# ==============================================================================
# step_2_grid_simulation.R
#
# Based on the working step_2_bc_multiset.R pattern.
# Added dimensions vs original:
#   - R0 sweep         : 1.25 to 3.00 by 0.25 (9 levels)
#   - seed_symptomatic : TRUE / FALSE (2 levels)
#
# Grid dimensions:
#   pathogen            : H1N1, SC2                          (2)
#   R0                  : 1.25, 1.50, ..., 3.00              (9)
#   antiviral_start     : 0, 5, 10, 20, 40                   (5)
#   efficacy_multiplier : 0.0, 0.2, 0.4, 0.6, 0.8, 1.0      (6)
#   intervention        : mild, moderate, intensive, maximum  (4)
#   seed_symptomatic    : TRUE, FALSE                         (2)
#   Total scenarios     : 2 x 9 x 5 x 6 x 4 x 2 = 4,320
#   Replicates/scenario : 500
#
# Outputs:
#   output/grid_results.parquet
# ==============================================================================

library(dplyr)
library(parallel)
library(arrow)

source("functions/network_bp_sim.R")

# ==============================================================================
# [Configuration Block]
# ==============================================================================

net_data_path <- "output/net_data_sf.rds"
output_dir    <- "output"
progress_dir  <- "output/.progress"
if (!dir.exists(output_dir)) dir.create(output_dir)

hh    <- 1
nonhh <- 0.1

# Pathogen definitions — functions inline as in original working code
pathogen_params <- list(
  
  H1N1 = list(
    label             = "H1N1",
    prop_asymptomatic = 0.50,
    generation_time_fn         = function(n) rgamma(n, shape = 4, rate = 1.538),
    infection_to_onset_fn      = function(n) rlnorm(n, meanlog = log(1.4), sdlog = log(1.51)),
    infectious_before_onset_fn = function(n) runif(n, 0.5, 0.7)
  ),
  
  SC2 = list(
    label             = "SARS-CoV-2",
    prop_asymptomatic = 0.15,
    generation_time_fn         = function(n) rweibull(n, shape = 2.826, scale = 5.665),
    infection_to_onset_fn      = function(n) rlnorm(n, meanlog = 1.621, sdlog = 0.418),
    infectious_before_onset_fn = function(n) runif(n, 0.4, 0.6)
  )
)

# Antiviral efficacy base values
base_eff_inf_e_max         <- 0.50
antiviral_eff_inf_kappa    <- 6.05
antiviral_eff_inf_tau      <- 0.32

base_eff_inf_hh_e_max      <- 0.71
antiviral_eff_inf_hh_kappa <- 6.91
antiviral_eff_inf_hh_tau   <- 1.10 

base_eff_trans_e_max       <- 0.70
antiviral_eff_trans_kappa  <- 0.91
antiviral_eff_trans_tau    <- 1.30

quarantine_efficacy     <- 0.9
n_sim                   <- 200L
epidemic_prob_threshold <- 0.999
seeding_cases           <- 1L
# n_cores                 <- max(1L, detectCores() - 1L)
n_cores                 <- 100

# ==============================================================================
# [Section 1] Load network
# ==============================================================================

cat(strrep("=", 62), "\n")
cat(" SF Branching Process — Scenario Grid\n")
cat(strrep("=", 62), "\n\n")

if (!file.exists(net_data_path))
  stop(sprintf("Network file not found: %s", net_data_path))

net_data <- readRDS(net_data_path)
n_people <- nrow(net_data$nodes)

cat(sprintf("Network      : %s\n",   net_data$network_type))
cat(sprintf("Population   : %d\n",   n_people))
cat(sprintf("Households   : %d\n",   net_data$n_households))
cat(sprintf("Cores        : %d\n\n", n_cores))

cat("Precomputing neighbor lists...\n")
neighbors_by_type <- precompute_neighbors(net_data)
cat("Done.\n\n")

# ==============================================================================
# [Section 2] Scenario grid
# ==============================================================================

intervention_levels <- list(
  mild = list(
    label                     = "Minimum",
    prob_treat_self           = 0.9,
    prob_treat_household      = 0.9,
    prob_treat_community      = 0.0,
    prob_quarantine_self      = 0.0,
    prob_quarantine_household = 0.0,
    prob_quarantine_community = 0.0
  ),
  moderate = list(
    label                     = "Moderate",
    prob_treat_self           = 0.9,
    prob_treat_household      = 0.9,
    prob_treat_community      = 0.9,
    prob_quarantine_self      = 0.0,
    prob_quarantine_household = 0.0,
    prob_quarantine_community = 0.0
  ),
  intensive = list(
    label                     = "Intensive",
    prob_treat_self           = 0.9,
    prob_treat_household      = 0.9,
    prob_treat_community      = 0.9,
    prob_quarantine_self      = 0.9,
    prob_quarantine_household = 0.9,
    prob_quarantine_community = 0.0
  ),
  maximum = list(
    label                     = "Maximum",
    prob_treat_self           = 0.9,
    prob_treat_household      = 0.9,
    prob_treat_community      = 0.9,
    prob_quarantine_self      = 0.9,
    prob_quarantine_household = 0.9,
    prob_quarantine_community = 0.9
  )
)

scenario_grid <- expand.grid(
  pathogen            = names(pathogen_params),
  R0                  = seq(1.00, 3.00, by = 0.5),
  # antiviral_start     = c(0, 5, 10, 20, 40),
  antiviral_start     = c(0, 14, 28),
  e_max_mult          = seq(0.0, 1.0, by = 0.25),
  intervention        = names(intervention_levels),
  seed_symptomatic    = c(TRUE, FALSE),
  stringsAsFactors    = FALSE
)
scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))
n_scenarios <- nrow(scenario_grid)

cat(sprintf("Total scenarios : %d\n",   n_scenarios))
cat(sprintf("Replicates      : %d\n",   n_sim))
cat(sprintf("Total runs      : %d\n\n", n_scenarios * n_sim))

# ==============================================================================
# [Section 3] run_scenario — same signature pattern as original working code
# ==============================================================================

run_scenario <- function(sc_row,
                         pathogen_params,
                         intervention_levels,
                         neighbors_by_type,
                         node_df,
                         hh, nonhh,
                         epidemic_prob_threshold,
                         base_eff_inf_e_max,
                         antiviral_eff_inf_kappa, antiviral_eff_inf_tau,
                         base_eff_inf_hh_e_max,
                         antiviral_eff_inf_hh_kappa, antiviral_eff_inf_hh_tau,
                         base_eff_trans_e_max,
                         antiviral_eff_trans_kappa, antiviral_eff_trans_tau,
                         quarantine_efficacy,
                         seeding_cases, n_sim,
                         progress_dir) {
  
  pg       <- pathogen_params[[sc_row$pathogen]]
  iv       <- intervention_levels[[sc_row$intervention]]
  emx      <- sc_row$e_max_mult
  av_start <- sc_row$antiviral_start
  R0       <- sc_row$R0
  
  prop_asymptomatic          <- pg$prop_asymptomatic
  generation_time_fn         <- pg$generation_time_fn
  infection_to_onset_fn      <- pg$infection_to_onset_fn
  infectious_before_onset_fn <- pg$infectious_before_onset_fn
  
  p_inf_household <- R0 / (2.14 * hh + 9.16 * nonhh)
  p_inf_community <- p_inf_household * nonhh
  
  eff_inf_e_max    <- base_eff_inf_e_max    * emx
  eff_inf_hh_e_max <- base_eff_inf_hh_e_max * emx
  eff_trans_e_max  <- base_eff_trans_e_max  * emx
  
  drug_eff_inf_fn <- function(t)
    eff_inf_e_max / (1 + exp(antiviral_eff_inf_kappa * (t - antiviral_eff_inf_tau)))
  
  drug_eff_inf_hh_fn <- function(t)
    eff_inf_hh_e_max / (1 + exp(antiviral_eff_inf_hh_kappa * (t - antiviral_eff_inf_hh_tau)))
  
  drug_eff_trans_fn <- function(t)
    eff_trans_e_max / (1 + exp(antiviral_eff_trans_kappa * (t - antiviral_eff_trans_tau)))
  
  delay_fn <- function(n) runif(n, 0, 3)
  
  reps <- lapply(seq_len(n_sim), function(rep_i) {
    
    res <- tryCatch(
      network_bp_sim(
        neighbors_by_type          = neighbors_by_type,
        node_df                    = node_df,
        R0                         = R0,
        epidemic_prob_threshold    = epidemic_prob_threshold,
        generation_time_fn         = generation_time_fn,
        infection_to_onset_fn      = infection_to_onset_fn,
        infectious_before_onset_fn = infectious_before_onset_fn,
        prop_asymptomatic          = prop_asymptomatic,
        p_inf_household            = p_inf_household,
        p_inf_community            = p_inf_community,
        antiviral_start            = av_start,
        logistical_delay_fn        = delay_fn,
        drug_eff_inf_fn            = drug_eff_inf_fn,
        drug_eff_inf_hh_fn         = drug_eff_inf_hh_fn,
        drug_eff_trans_fn          = drug_eff_trans_fn,
        eff_inf_e_max              = eff_inf_e_max,
        eff_inf_hh_e_max           = eff_inf_hh_e_max,
        eff_trans_e_max            = eff_trans_e_max,
        prob_treat_self            = iv$prob_treat_self,
        prob_treat_household       = iv$prob_treat_household,
        prob_treat_community       = iv$prob_treat_community,
        time_to_quarantine_fn      = delay_fn,
        prob_quarantine_self       = iv$prob_quarantine_self,
        prob_quarantine_household  = iv$prob_quarantine_household,
        prob_quarantine_community  = iv$prob_quarantine_community,
        quarantine_efficacy        = quarantine_efficacy,
        seeding_cases              = seeding_cases,
        seed_symptomatic           = sc_row$seed_symptomatic,
        t0                         = 0,
        seed                       = rep_i
      ),
      error = function(e) NULL
    )
    
    if (is.null(res)) {
      return(data.frame(
        scenario_id       = sc_row$scenario_id,
        rep               = rep_i,
        is_epidemic       = NA,
        n_infected        = NA_integer_,
        max_generation    = NA_integer_,
        epidemic_duration = NA_real_,
        active_at_end     = NA_integer_,
        error             = TRUE
      ))
    }
    
    inf_df <- res$infected
    n_inf  <- nrow(inf_df)
    
    data.frame(
      scenario_id       = sc_row$scenario_id,
      rep               = rep_i,
      is_epidemic       = res$is_epidemic,
      n_infected        = n_inf,
      max_generation    = if (n_inf > 0) max(inf_df$generation, na.rm = TRUE) else 0L,
      epidemic_duration = if (n_inf > 0)
        max(inf_df$time_infection, na.rm = TRUE) -
        min(inf_df$time_infection, na.rm = TRUE) else 0.0,
      active_at_end     = res$active_at_end,
      error             = FALSE
    )
  })
  
  writeLines("done", file.path(progress_dir,
                               sprintf("sc_%04d.done", sc_row$scenario_id)))
  
  do.call(rbind, reps)
}

# ==============================================================================
# [Section 4] Cluster — exact same export pattern as original working code
# ==============================================================================

if (dir.exists(progress_dir)) unlink(progress_dir, recursive = TRUE)
dir.create(progress_dir)

cat(sprintf("Launching cluster with %d workers...\n", n_cores))
cl <- makeCluster(n_cores)

clusterExport(cl, varlist = c(
  "run_scenario",
  "network_bp_sim", "precompute_neighbors",
  "pathogen_params", "intervention_levels",
  "neighbors_by_type",
  "scenario_grid",
  "hh", "nonhh",
  "epidemic_prob_threshold",
  "base_eff_inf_e_max",
  "antiviral_eff_inf_kappa", "antiviral_eff_inf_tau",
  "base_eff_inf_hh_e_max",
  "antiviral_eff_inf_hh_kappa", "antiviral_eff_inf_hh_tau",
  "base_eff_trans_e_max",
  "antiviral_eff_trans_kappa", "antiviral_eff_trans_tau",
  "quarantine_efficacy",
  "seeding_cases", "n_sim", "progress_dir"
), envir = environment())

clusterExport(cl, varlist = "node_df",
              envir = list2env(list(node_df = net_data$nodes)))

cat(sprintf("Running %d scenarios x %d reps...\n", n_scenarios, n_sim))
t_start <- proc.time()

chunk_size   <- n_cores * 4
n_chunks     <- ceiling(n_scenarios / chunk_size)
results_list <- vector("list", n_scenarios)

for (chunk_i in seq_len(n_chunks)) {
  
  idx_start <- (chunk_i - 1) * chunk_size + 1
  idx_end   <- min(chunk_i * chunk_size, n_scenarios)
  chunk_ids <- idx_start:idx_end
  
  chunk_results <- parLapply(cl, chunk_ids, function(i) {
    run_scenario(
      sc_row                     = scenario_grid[i, ],
      pathogen_params            = pathogen_params,
      intervention_levels        = intervention_levels,
      neighbors_by_type          = neighbors_by_type,
      node_df                    = node_df,
      hh                         = hh,
      nonhh                      = nonhh,
      epidemic_prob_threshold    = epidemic_prob_threshold,
      base_eff_inf_e_max         = base_eff_inf_e_max,
      antiviral_eff_inf_kappa    = antiviral_eff_inf_kappa,
      antiviral_eff_inf_tau      = antiviral_eff_inf_tau,
      base_eff_inf_hh_e_max      = base_eff_inf_hh_e_max,
      antiviral_eff_inf_hh_kappa = antiviral_eff_inf_hh_kappa,
      antiviral_eff_inf_hh_tau   = antiviral_eff_inf_hh_tau,
      base_eff_trans_e_max       = base_eff_trans_e_max,
      antiviral_eff_trans_kappa  = antiviral_eff_trans_kappa,
      antiviral_eff_trans_tau    = antiviral_eff_trans_tau,
      quarantine_efficacy        = quarantine_efficacy,
      seeding_cases              = seeding_cases,
      n_sim                      = n_sim,
      progress_dir               = progress_dir
    )
  })
  
  results_list[chunk_ids] <- chunk_results
  
  # Progress bar
  n_done  <- length(list.files(progress_dir, pattern = "\\.done$"))
  pct     <- n_done / n_scenarios * 100
  elapsed <- (proc.time() - t_start)[["elapsed"]]
  eta     <- if (n_done > 0) elapsed / n_done * (n_scenarios - n_done) else NA
  eta_str <- if (!is.na(eta))
    sprintf("ETA: %.0f sec", eta) else "ETA: --"
  bar     <- round(pct / 2)
  cat(sprintf("\r  |%s%s| %5.1f%%  (%d/%d)  %s     ",
              strrep("=", bar), strrep(" ", 50 - bar),
              pct, n_done, n_scenarios, eta_str))
  flush.console()
}

stopCluster(cl)
unlink(progress_dir, recursive = TRUE)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nDone. %.1f sec (%.1f min)\n\n", elapsed, elapsed / 60))

# ==============================================================================
# [Section 5] Assemble and save
# ==============================================================================

results_df <- do.call(rbind, results_list)

results_full <- results_df %>%
  left_join(
    scenario_grid %>% select(
      scenario_id, pathogen, R0, antiviral_start,
      e_max_mult, intervention, seed_symptomatic
    ),
    by = "scenario_id"
  ) %>%
  mutate(
    intervention_label = sapply(intervention, function(x)
      intervention_levels[[x]]$label)
  ) %>%
  select(
    scenario_id, rep,
    pathogen, R0, antiviral_start,
    e_max_mult, intervention, intervention_label,
    seed_symptomatic,
    is_epidemic, n_infected, max_generation,
    epidemic_duration, active_at_end, error
  )

out_path <- file.path(output_dir, "grid_results.parquet")
write_parquet(results_full, out_path)
cat(sprintf("Saved : %s\n",   out_path))
cat(sprintf("Rows  : %d\n",   nrow(results_full)))
cat(sprintf("Errors: %d\n\n", sum(results_full$error, na.rm = TRUE)))

# ==============================================================================
# [Section 6] Quick summary
# ==============================================================================

cat(strrep("-", 62), "\n")
cat(" Summary by pathogen x intervention\n")
cat(strrep("-", 62), "\n")

results_full %>%
  filter(!error) %>%
  group_by(pathogen, intervention_label) %>%
  summarise(
    p_epidemic      = mean(is_epidemic,    na.rm = TRUE),
    mean_n_infected = mean(n_infected,     na.rm = TRUE),
    mean_max_gen    = mean(max_generation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  { for (i in seq_len(nrow(.)))
    cat(sprintf("  %-14s | %-10s : P(ep)=%.3f  mean_inf=%.1f  mean_gen=%.1f\n",
                .[[i,"pathogen"]], .[[i,"intervention_label"]],
                .[[i,"p_epidemic"]], .[[i,"mean_n_infected"]],
                .[[i,"mean_max_gen"]])) }

cat(strrep("=", 62), "\n\n")
cat("Done.\n")
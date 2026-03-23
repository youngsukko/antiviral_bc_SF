# ==============================================================================
# step_2_bc_multiset
#
# Purpose:
#   Run a full scenario grid on the SF-Prem network with parallel Monte Carlo.
#   Uses makeCluster + parLapply for Windows-compatible multicore execution.
#
# Grid dimensions:
#   - Pathogen            : 2 (H1N1 / SARS-CoV-2)
#   - Intervention level  : 4 (mild / moderate / intensive / maximum)
#   - e_max multiplier    : 6 (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)  — applied equally to
#                           antiviral_eff_inf_e_max, antiviral_eff_inf_hh_e_max,
#                           and antiviral_eff_trans_e_max
#   - Delay upper bound   : 4 (Uniform(0, 1/2/3/4 days)) — applied equally to
#                           logistical_delay_fn and quarantine_fn
#   - Antiviral start day : 5 (0, 7, 14, 21, 28)
#   Total scenarios       : 2 x 4 x 6 x 4 x 5 = 960
#   Each scenario         : n_sim replications (parallel across scenarios)
#
# Pathogen definitions:
#   H1N1        : 2009 pandemic influenza A/H1N1
#                   R0 = 1.5
#                   Generation time : Gamma(shape=4, rate=1.538), mean ~2.6 d
#                     anchored to symptom onset (White et al. 2009)
#                   Incubation      : LogNormal(meanlog=log(1.4), sdlog=log(1.51))
#                     (Lessler et al. 2009)
#                   Presymptomatic infectious proportion : Uniform(0.5, 0.7)
#                   prop_asymptomatic = 0.50
#
#   SARS-CoV-2  : ancestral Wuhan lineage
#                   R0 = 2.79
#                   Generation time : Weibull(shape=2.826, scale=5.665), mean ~5.0 d
#                     anchored to symptom onset
#                     (Ferretti et al., Science 2020;
#                      https://www.science.org/doi/10.1126/science.abb6936)
#                   Incubation      : LogNormal(meanlog=1.621, sdlog=0.418),
#                                     mean ~5.5 d, median ~5.1 d
#                     (Lauer et al., Ann Intern Med 2020;
#                      https://www.acpjournals.org/doi/10.7326/M20-0504)
#                   Presymptomatic infectious proportion : Uniform(0.4, 0.6)
#                   prop_asymptomatic = 0.15
#                   Note: antiviral efficacy parameters (kappa, tau, e_max) are
#                   held identical to H1N1 across the grid. The pathogen dimension
#                   isolates the effect of Tg, presymptomatic fraction, and R0
#                   on intervention effectiveness under identical assumed drug
#                   properties.
#
# Intervention level definitions (prob values are 0 or 0.9):
#   Mild      : AV self only
#   Moderate  : AV self+HH,              quarantine self
#   Intensive : AV self+HH+community,    quarantine self+HH
#   Maximum   : AV self+HH+community,    quarantine self+HH+community
#
# Epidemic termination rule (Whittle 1955 — asymptomatic-only threshold):
#   The Whittle threshold is evaluated against ASYMPTOMATIC active infections
#   only. Symptomatic cases trigger ring PEP and self-isolation; their
#   transmission potential is substantially curtailed by the intervention model
#   and they are excluded from the R0-based epidemic probability calculation.
#
# Progress tracking:
#   Each worker writes a small file to progress_dir when a scenario completes.
#   The main process polls the file count between chunk batches and prints
#   a progress bar to the console.
#
# Integer encoding (defined in network_bp_sim.R):
#   status            : 1=S, 2=I, 3=R
#   contact_type      : 0=seed, 1=household, 2=community
#   treat_reason      : 1=self, 2=household, 3=community
#   quarantine_reason : 1=self, 2=household, 3=community
#
# Outputs (output/ folder):
#   scenario_results.rds : full results list (one row per replication)
#   scenario_summary.csv : aggregated summary per scenario
# ==============================================================================

library(dplyr)
library(parallel)
library(arrow)

source("functions/network_bp_sim.R")


# ==============================================================================
# [Configuration Block] — Modify only this section
# ==============================================================================

# --- Paths --------------------------------------------------------------------
net_data_path <- "output/net_data_sf.rds"
output_dir    <- "output"
progress_dir  <- "output/.progress"   # temporary per-scenario flag files
if (!dir.exists(output_dir))  dir.create(output_dir)

# --- Network degree parameters (shared across pathogens) ----------------------
hh    <- 1      # relative household transmission weight
nonhh <- 0.1   # relative community transmission weight
# p_inf_household and p_inf_community are computed inside run_scenario()
# using each pathogen's R0, so they are not set globally here.

# --- Pathogen parameter sets --------------------------------------------------
#
# Each pathogen list contains all parameters that vary between pathogens.
# Antiviral efficacy parameters (kappa, tau, e_max) are shared across
# pathogens — the grid sweep isolates the effect of pathogen biology
# (Tg, incubation, R0, prop_asymptomatic, presymptomatic fraction)
# on intervention outcome.
#
# Generation time is anchored to symptom onset (or time_infectious as proxy
# for asymptomatic sources), consistent with network_bp_sim.R:
#   t_idx = onset for symptomatic, t_infectious for asymptomatic
#
# infectious_before_onset_fn(n): draws the proportion of the incubation period
#   elapsed before the person becomes infectious (i.e., latent fraction).
#   time_infectious = t_infection + proportion * incubation_period

pathogen_params <- list(
  
  H1N1 = list(
    label             = "H1N1",
    R0                = 1.5,
    prop_asymptomatic = 0.50,
    # Gamma generation time anchored to onset: mean ~2.6 d (White et al. 2009)
    generation_time_fn         = function(n) rgamma(n, shape = 4, rate = 1.538),
    # Log-normal incubation (Lessler et al. 2009)
    infection_to_onset_fn      = function(n)
      rlnorm(n, meanlog = log(1.4), sdlog = log(1.51)),
    # ~10% of H1N1 transmission is presymptomatic; latent fraction ~0.5–0.7
    infectious_before_onset_fn = function(n) runif(n, 0.5, 0.7)
  ),
  
  SC2 = list(
    label             = "SARS-CoV-2",
    R0                = 2.79,
    prop_asymptomatic = 0.15,
    # Weibull generation time anchored to onset: mean ~5.0 d, SD ~1.9 d
    # Ferretti et al., Science 2020
    # (https://www.science.org/doi/10.1126/science.abb6936)
    generation_time_fn         = function(n) rweibull(n, shape = 2.826, scale = 5.665),
    # Log-normal incubation: mean ~5.5 d, median ~5.1 d
    # Lauer et al., Ann Intern Med 2020
    # (https://www.acpjournals.org/doi/10.7326/M20-0504)
    infection_to_onset_fn      = function(n) rlnorm(n, meanlog = 1.621, sdlog = 0.418),
    # ~46% of SC2 transmission is presymptomatic; latent fraction ~0.4–0.6
    infectious_before_onset_fn = function(n) runif(n, 0.4, 0.6)
  )
  
)

# --- Antiviral efficacy — community infection prevention (base values) ---------
# Multiplied by e_max_mult in the scenario grid sweep.
base_eff_inf_e_max      <- 0.50
antiviral_eff_inf_kappa <- 7.51
antiviral_eff_inf_tau   <- 1.34

# --- Antiviral efficacy — household infection prevention (base values) ---------
# Separate logistic; household PEP trials report higher efficacy.
# Also multiplied by e_max_mult (proportional scaling assumption).
base_eff_inf_hh_e_max      <- 0.63
antiviral_eff_inf_hh_kappa <- 6.06
antiviral_eff_inf_hh_tau   <- 0.89

# --- Antiviral efficacy — transmission reduction (base values) ----------------
base_eff_trans_e_max      <- 0.70
antiviral_eff_trans_kappa <- 0.91
antiviral_eff_trans_tau   <- 1.30

quarantine_efficacy <- 0.9

# --- Simulation control -------------------------------------------------------
n_sim                   <- 500
epidemic_prob_threshold <- 0.999
seeding_cases           <- 1
# n_cores                 <- max(1L, detectCores() - 1L)
n_cores                 <- 4


# ==============================================================================
# [Section 1] Load network and precompute neighbors (done once)
# ==============================================================================

cat(strrep("=", 62), "\n")
cat(" SF Branching Process — Scenario Grid\n")
cat(strrep("=", 62), "\n\n")

if (!file.exists(net_data_path))
  stop(sprintf("Network file not found: %s", net_data_path))

net_data <- readRDS(net_data_path)
n_people <- nrow(net_data$nodes)

cat(sprintf("Network         : %s\n", net_data$network_type))
cat(sprintf("Population      : %d\n", n_people))
cat(sprintf("Households      : %d\n", net_data$n_households))
cat(sprintf("Pathogens       : %s\n",
            paste(sapply(pathogen_params, `[[`, "label"), collapse = ", ")))
cat(sprintf("Replications    : %d per scenario\n", n_sim))
cat(sprintf("Cores           : %d\n\n", n_cores))

cat("Precomputing neighbor lists...\n")
neighbors_by_type <- precompute_neighbors(net_data)
cat("Done.\n\n")


# ==============================================================================
# [Section 2] Build scenario grid
# ==============================================================================

intervention_levels <- list(
  mild = list(
    label                     = "Mild",
    prob_treat_self           = 0.9,
    prob_treat_household      = 0.0,
    prob_treat_community      = 0.0,
    prob_quarantine_self      = 0.0,
    prob_quarantine_household = 0.0,
    prob_quarantine_community = 0.0
  ),
  moderate = list(
    label                     = "Moderate",
    prob_treat_self           = 0.9,
    prob_treat_household      = 0.9,
    prob_treat_community      = 0.0,
    prob_quarantine_self      = 0.9,
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

e_max_multipliers  <- seq(0, 1.0, by = 0.25)
# delay_upper_bounds <- c(1, 2, 3, 4)
delay_upper_bounds <- c(1, 4)
# av_start_days      <- c(0, 7, 14, 21, 28)
av_start_days      <- c(0, 14, 28)

scenario_grid <- expand.grid(
  pathogen        = names(pathogen_params),
  intervention    = names(intervention_levels),
  e_max_mult      = e_max_multipliers,
  delay_max       = delay_upper_bounds,
  antiviral_start = av_start_days,
  stringsAsFactors = FALSE
)
scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))

n_scenarios <- nrow(scenario_grid)
cat(sprintf(
  "Scenario grid   : %d scenarios (%d pathogens x %d interventions x %d e_max x %d delays x %d AV starts)\n",
  n_scenarios,
  length(pathogen_params),
  length(intervention_levels),
  length(e_max_multipliers),
  length(delay_upper_bounds),
  length(av_start_days)
))
cat(sprintf("Total runs      : %d\n\n", n_scenarios * n_sim))


# ==============================================================================
# [Section 3] Run simulations (Windows-compatible cluster)
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
  
  # Resolve pathogen-specific parameters
  pg       <- pathogen_params[[sc_row$pathogen]]
  iv       <- intervention_levels[[sc_row$intervention]]
  emx      <- sc_row$e_max_mult
  dmax     <- sc_row$delay_max
  av_start <- sc_row$antiviral_start
  
  R0                         <- pg$R0
  prop_asymptomatic          <- pg$prop_asymptomatic
  generation_time_fn         <- pg$generation_time_fn
  infection_to_onset_fn      <- pg$infection_to_onset_fn
  infectious_before_onset_fn <- pg$infectious_before_onset_fn
  
  # Per-contact SAR derived from R0 and network degree structure.
  # Household degree ~2.14, community degree ~9.16 (SF-Prem network).
  p_inf_household <- R0 / (2.14 * hh + 9.16 * nonhh)
  p_inf_community <- p_inf_household * nonhh
  
  # Scale all e_max values by the same multiplier.
  # This keeps the household/community efficacy ratio fixed across the grid
  # while allowing a uniform sweep of overall antiviral strength.
  eff_inf_e_max    <- base_eff_inf_e_max    * emx
  eff_inf_hh_e_max <- base_eff_inf_hh_e_max * emx
  eff_trans_e_max  <- base_eff_trans_e_max  * emx
  
  drug_eff_inf_fn <- function(t)
    eff_inf_e_max / (1 + exp(antiviral_eff_inf_kappa * (t - antiviral_eff_inf_tau)))
  
  drug_eff_inf_hh_fn <- function(t)
    eff_inf_hh_e_max / (1 + exp(antiviral_eff_inf_hh_kappa * (t - antiviral_eff_inf_hh_tau)))
  
  drug_eff_trans_fn <- function(t)
    eff_trans_e_max / (1 + exp(antiviral_eff_trans_kappa * (t - antiviral_eff_trans_tau)))
  
  delay_fn <- function(n) runif(n, 0, dmax)
  
  reps <- lapply(seq_len(n_sim), function(rep_i) {
    
    res <- network_bp_sim(
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
      t0                         = 0,
      seed                       = rep_i
    )
    
    data.frame(
      scenario_id     = sc_row$scenario_id,
      pathogen        = sc_row$pathogen,
      intervention    = sc_row$intervention,
      e_max_mult      = emx,
      delay_max       = dmax,
      antiviral_start = av_start,
      rep             = rep_i,
      is_epidemic     = res$is_epidemic,
      n_infected      = nrow(res$infected),
      active_at_end   = res$active_at_end
    )
  })
  
  # Write a completion flag file so the main process can track progress
  writeLines("done", file.path(progress_dir, sprintf("sc_%04d.done", sc_row$scenario_id)))
  
  do.call(rbind, reps)
}

# --------------------------------------------------------------------------
# Set up cluster, export everything workers need, then run
# --------------------------------------------------------------------------

# Fresh progress directory
if (dir.exists(progress_dir)) unlink(progress_dir, recursive = TRUE)
dir.create(progress_dir)

cat(sprintf("Launching cluster with %d workers...\n", n_cores))
cl <- makeCluster(n_cores)

# Export all objects the worker function references.
# infectious_before_onset_fn lives inside pathogen_params, so exporting
# pathogen_params is sufficient — no separate export needed.
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

# Chunk the scenario list into n_cores-sized batches.
# Progress is reported between chunks — parLapply is blocking so
# real-time polling within a chunk is not possible in base R.
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
  
  # Progress report after each batch
  n_done   <- length(list.files(progress_dir, pattern = "\\.done$"))
  pct      <- n_done / n_scenarios * 100
  elapsed  <- (proc.time() - t_start)[["elapsed"]]
  eta      <- if (n_done > 0) elapsed / n_done * (n_scenarios - n_done) else NA
  eta_str  <- if (!is.na(eta)) sprintf("ETA: %s", format(as.difftime(eta, units = "secs"), digits = 3)) else "ETA: --"
  bar_done <- round(pct / 2)
  bar_str  <- paste0(strrep("=", bar_done), strrep(" ", 50 - bar_done))
  cat(sprintf("\r  |%s| %5.1f%%  (%d/%d)  %s     ",
              bar_str, pct, n_done, n_scenarios, eta_str))
  flush.console()
}

stopCluster(cl)
unlink(progress_dir, recursive = TRUE)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nDone. Total elapsed: %.1f sec\n\n", elapsed))


# ==============================================================================
# [Section 4] Aggregate and save results
# ==============================================================================

results_df <- do.call(rbind, results_list)

summary_df <- results_df %>%
  group_by(scenario_id, pathogen, intervention, e_max_mult, delay_max, antiviral_start) %>%
  summarise(
    n_sim           = n(),
    p_epidemic      = mean(is_epidemic),
    mean_infected   = mean(n_infected),
    median_infected = median(n_infected),
    q25_infected    = quantile(n_infected, 0.25),
    q75_infected    = quantile(n_infected, 0.75),
    .groups         = "drop"
  ) %>%
  mutate(
    pathogen = factor(pathogen,
                      levels = c("H1N1", "SC2"),
                      labels = c("H1N1", "SARS-CoV-2")),
    intervention = factor(intervention,
                          levels = c("mild", "moderate", "intensive", "maximum"),
                          labels = c("Mild", "Moderate", "Intensive", "Maximum"))
  ) %>%
  arrange(scenario_id)

path_rds <- file.path(output_dir, "scenario_results.rds")
path_csv <- file.path(output_dir, "scenario_summary.csv")
saveRDS(results_df, path_rds)
write.csv(summary_df, path_csv, row.names = FALSE)

cat(sprintf("Saved: %s\n", path_rds))
cat(sprintf("Saved: %s\n\n", path_csv))


# ==============================================================================
# [Section 5] Console summary
# ==============================================================================

cat(strrep("-", 62), "\n")
cat(" Epidemic probability by pathogen\n")
cat(" (averaged across all intervention, e_max, delay, and AV start)\n")
cat(strrep("-", 62), "\n")

summary_df %>%
  group_by(pathogen) %>%
  summarise(
    mean_p_epidemic = mean(p_epidemic),
    mean_n_infected = mean(mean_infected),
    .groups = "drop"
  ) %>%
  { for (i in seq_len(nrow(.))) {
    cat(sprintf("  %-14s : P(epidemic) = %.3f | Mean infected = %.1f\n",
                as.character(.[[i, "pathogen"]]),
                .[[i, "mean_p_epidemic"]],
                .[[i, "mean_n_infected"]]))
  }}

cat("\n")
cat(strrep("-", 62), "\n")
cat(" Epidemic probability by pathogen x intervention level\n")
cat(" (averaged across all e_max, delay, and AV start combinations)\n")
cat(strrep("-", 62), "\n")

summary_df %>%
  group_by(pathogen, intervention) %>%
  summarise(
    mean_p_epidemic = mean(p_epidemic),
    mean_n_infected = mean(mean_infected),
    .groups = "drop"
  ) %>%
  { for (i in seq_len(nrow(.))) {
    cat(sprintf("  %-14s | %-12s : P(epidemic) = %.3f | Mean infected = %.1f\n",
                as.character(.[[i, "pathogen"]]),
                as.character(.[[i, "intervention"]]),
                .[[i, "mean_p_epidemic"]],
                .[[i, "mean_n_infected"]]))
  }}

cat("\n")
cat(strrep("-", 62), "\n")
cat(" Epidemic probability by AV start day\n")
cat(" (averaged across all pathogen, intervention, e_max, and delay)\n")
cat(strrep("-", 62), "\n")

summary_df %>%
  group_by(antiviral_start) %>%
  summarise(
    mean_p_epidemic = mean(p_epidemic),
    mean_n_infected = mean(mean_infected),
    .groups = "drop"
  ) %>%
  { for (i in seq_len(nrow(.))) {
    cat(sprintf("  AV start day %-3d : P(epidemic) = %.3f | Mean infected = %.1f\n",
                .[[i, "antiviral_start"]],
                .[[i, "mean_p_epidemic"]],
                .[[i, "mean_n_infected"]]))
  }}

cat(strrep("-", 62), "\n\n")
cat(strrep("=", 62), "\n")
cat(" Done\n")
cat(strrep("=", 62), "\n\n")
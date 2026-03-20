# ==============================================================================
# step_1_bc_singlerun.R
#
# Purpose:
#   Load the SF-Prem network (output/net_data_sf.rds) built by
#   scp_build_network_sf.R, run one realization of the stochastic
#   branching process simulation, and visualize the results.
#
# Run this file standalone.
# Prerequisites:
#   - scp_build_network_sf.R must have been run first
#   - functions/network_bp_sim.R must exist
#
# Recent changes reflected:
#   - infectious_before_onset_fn added (latent period as proportion of incubation)
#   - pathogen block added (H1N1 / SARS-CoV-2 switchable)
#   - generation time now anchored to onset (t_idx = onset for symptomatic,
#     t_infectious for asymptomatic), consistent with network_bp_sim.R
#
# Integer encoding (defined in network_bp_sim.R):
#   status            : 1=S, 2=I, 3=R
#   contact_type      : 0=seed, 1=household, 2=community
#   treat_reason      : 1=self, 2=household, 3=community
#   quarantine_reason : 1=self, 2=household, 3=community
# ==============================================================================

library(dplyr)
library(ggplot2)

source("functions/network_bp_sim.R")


# ==============================================================================
# [Configuration Block] — Modify only this section
# ==============================================================================

# --- Pathogen -----------------------------------------------------------------
# Switch between "H1N1" and "SC2" to change pathogen.
pathogen <- "H1N1"   # "H1N1" or "SC2"

pathogen_params <- list(
  H1N1 = list(
    label             = "H1N1 (2009 pandemic)",
    R0                = 1.5,
    prop_asymptomatic = 0.50,
    # Gamma generation time anchored to onset: mean ~2.6 d (White et al. 2009)
    generation_time_fn        = function(n) rgamma(n, shape = 4, rate = 1.538),
    # Log-normal incubation (Lessler et al. 2009)
    infection_to_onset_fn     = function(n) rlnorm(n, meanlog = log(1.4), sdlog = log(1.51)),
    # Proportion of incubation period elapsed before becoming infectious
    # ~10% of H1N1 transmission is presymptomatic (Cowling et al. 2010)
    infectious_before_onset_fn = function(n) runif(n, 0.5, 0.7)
  ),
  SC2 = list(
    label             = "SARS-CoV-2 (ancestral)",
    R0                = 2.79,
    prop_asymptomatic = 0.15,
    # Weibull generation time anchored to onset: mean ~5.0 d, SD ~1.9 d
    # Ferretti et al., Science 2020
    # (https://www.science.org/doi/10.1126/science.abb6936)
    generation_time_fn        = function(n) rweibull(n, shape = 2.826, scale = 5.665),
    # Log-normal incubation: mean ~5.5 d, median ~5.1 d
    # Lauer et al., Ann Intern Med 2020
    # (https://www.acpjournals.org/doi/10.7326/M20-0504)
    infection_to_onset_fn     = function(n) rlnorm(n, meanlog = 1.621, sdlog = 0.418),
    # ~46% of SC2 transmission is presymptomatic (Ferretti et al. 2020)
    infectious_before_onset_fn = function(n) runif(n, 0.4, 0.6)
  )
)

pg <- pathogen_params[[pathogen]]

# --- Network degree parameters ------------------------------------------------
hh    <- 1
nonhh <- 0.1

# --- Derived parameters -------------------------------------------------------
R0                         <- pg$R0
p_inf_household            <- R0 / (2.14 * hh + 9.16 * nonhh)
p_inf_community            <- p_inf_household * nonhh
prop_asymptomatic          <- pg$prop_asymptomatic
generation_time_fn         <- pg$generation_time_fn
infection_to_onset_fn      <- pg$infection_to_onset_fn
infectious_before_onset_fn <- pg$infectious_before_onset_fn

# --- Antiviral efficacy — community infection prevention ----------------------
# Logistic decay: E(t) = e_max / (1 + exp(kappa * (t - tau)))
antiviral_eff_inf_e_max   <- 0.50
antiviral_eff_inf_kappa   <- 7.51
antiviral_eff_inf_tau     <- 1.34

# --- Antiviral efficacy — household infection prevention ----------------------
# Separate logistic; household PEP trials report higher efficacy.
antiviral_eff_inf_hh_e_max  <- 0.63
antiviral_eff_inf_hh_kappa  <- 6.06
antiviral_eff_inf_hh_tau    <- 0.89

# --- Antiviral efficacy — transmission reduction (source side) ----------------
# Anchored to symptom onset (or time_infectious proxy for asymptomatic).
antiviral_eff_trans_e_max <- 0.70
antiviral_eff_trans_kappa <- 0.91
antiviral_eff_trans_tau   <- 1.30

# Efficacy functions
drug_eff_inf_fn <- function(t)
  antiviral_eff_inf_e_max / (1 + exp(antiviral_eff_inf_kappa * (t - antiviral_eff_inf_tau)))

drug_eff_inf_hh_fn <- function(t)
  antiviral_eff_inf_hh_e_max / (1 + exp(antiviral_eff_inf_hh_kappa * (t - antiviral_eff_inf_hh_tau)))

drug_eff_trans_fn <- function(t)
  antiviral_eff_trans_e_max / (1 + exp(antiviral_eff_trans_kappa * (t - antiviral_eff_trans_tau)))

# --- Intervention setting -----------------------------------------------------
antiviral_start           <- 14

prob_treat_self           <- 0.5
prob_treat_household      <- 0.5
prob_treat_community      <- 0.5

prob_quarantine_self      <- 0.0
prob_quarantine_household <- 0.0
prob_quarantine_community <- 0.0

quarantine_efficacy <- 0.0

logistical_delay_fn <- function(n) runif(n, 0, 3)
quarantine_fn       <- function(n) runif(n, 0, 3)

# --- Simulation control -------------------------------------------------------
epidemic_prob_threshold <- 0.999   # evaluated on asymptomatic active cases only
seeding_cases           <- 1
sim_seed                <- 42

# --- Paths --------------------------------------------------------------------
net_data_path <- "output/net_data_sf.rds"
figure_dir    <- "figures"
if (!dir.exists(figure_dir)) dir.create(figure_dir)


# ==============================================================================
# [Section 2] Load network and precompute neighbors
# ==============================================================================

cat(strrep("=", 62), "\n")
cat(" SF Branching Process Simulation — Single Run\n")
cat(strrep("=", 62), "\n\n")

if (!file.exists(net_data_path))
  stop(sprintf(
    "Network file not found: %s\nPlease run scp_build_network_sf.R first.",
    net_data_path
  ))

net_data <- readRDS(net_data_path)
n_people <- nrow(net_data$nodes)

epidemic_queue_threshold <- if (R0 > 1) {
  ceiling(log(1 - epidemic_prob_threshold) / log(1 / R0))
} else Inf

cat(sprintf("Network loaded  : %s\n",   net_data_path))
cat(sprintf("Network type    : %s\n",   net_data$network_type))
cat(sprintf("Population      : %d\n",   n_people))
cat(sprintf("Households      : %d\n",   net_data$n_households))
cat(sprintf("Pathogen        : %s\n",   pg$label))
cat(sprintf("R0              : %.2f\n", R0))
cat(sprintf(
  "Epidemic cutoff : %d asymptomatic active infections (Whittle threshold at %.1f%%)\n\n",
  epidemic_queue_threshold, epidemic_prob_threshold * 100
))

cat("Precomputing neighbor lists...\n")
neighbors_by_type <- precompute_neighbors(net_data)
cat("Done.\n\n")


# ==============================================================================
# [Section 3] Run simulation
# ==============================================================================

cat("Parameters\n")
cat(sprintf("  p_inf_household             : %.4f\n", p_inf_household))
cat(sprintf("  p_inf_community             : %.4f\n", p_inf_community))
cat(sprintf("  prop_asymptomatic           : %.2f\n",  prop_asymptomatic))
cat(sprintf("  antiviral_start             : day %d\n", antiviral_start))
cat(sprintf("  eff_inf (community) e_max   : %.2f\n",  antiviral_eff_inf_e_max))
cat(sprintf("  eff_inf (household) e_max   : %.2f\n",  antiviral_eff_inf_hh_e_max))
cat(sprintf("  eff_trans e_max             : %.2f\n",  antiviral_eff_trans_e_max))
cat(sprintf("  seeding_cases               : %d\n",    seeding_cases))
cat(sprintf("  sim_seed                    : %d\n\n",  sim_seed))

t_start <- proc.time()

result <- network_bp_sim(
  neighbors_by_type          = neighbors_by_type,
  node_df                    = net_data$nodes,
  R0                         = R0,
  epidemic_prob_threshold    = epidemic_prob_threshold,
  generation_time_fn         = generation_time_fn,
  infection_to_onset_fn      = infection_to_onset_fn,
  infectious_before_onset_fn = infectious_before_onset_fn,
  prop_asymptomatic          = prop_asymptomatic,
  p_inf_household            = p_inf_household,
  p_inf_community            = p_inf_community,
  antiviral_start            = antiviral_start,
  logistical_delay_fn        = logistical_delay_fn,
  drug_eff_inf_fn            = drug_eff_inf_fn,
  drug_eff_inf_hh_fn         = drug_eff_inf_hh_fn,
  drug_eff_trans_fn          = drug_eff_trans_fn,
  eff_inf_e_max              = antiviral_eff_inf_e_max,
  eff_inf_hh_e_max           = antiviral_eff_inf_hh_e_max,
  eff_trans_e_max            = antiviral_eff_trans_e_max,
  prob_treat_self            = prob_treat_self,
  prob_treat_household       = prob_treat_household,
  prob_treat_community       = prob_treat_community,
  time_to_quarantine_fn      = quarantine_fn,
  prob_quarantine_self       = prob_quarantine_self,
  prob_quarantine_household  = prob_quarantine_household,
  prob_quarantine_community  = prob_quarantine_community,
  quarantine_efficacy        = quarantine_efficacy,
  seeding_cases              = seeding_cases,
  t0                         = 0,
  seed                       = sim_seed
)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("Simulation time : %.1f sec\n\n", elapsed))


# ==============================================================================
# [Section 4] Results summary
# ==============================================================================

tdf_full     <- result$full
tdf_infected <- result$infected
n_infected   <- nrow(tdf_infected)
is_epidemic  <- result$is_epidemic

ct_labels <- c("0" = "seed", "1" = "household", "2" = "community")
tr_labels <- c("1" = "self", "2" = "household", "3" = "community")

cat(strrep("-", 62), "\n")
cat(" Results\n")
cat(strrep("-", 62), "\n")
cat(sprintf("  Total infected    : %d / %d (%.1f%%)\n",
            n_infected, n_people, n_infected / n_people * 100))
cat(sprintf("  Outcome           : %s\n",
            ifelse(is_epidemic,
                   sprintf("EPIDEMIC (asymptomatic active queue reached %d)", result$active_at_end),
                   "Contained (natural extinction)")))

if (n_infected > 0) {
  cat(sprintf("  Epidemic duration : %.1f - %.1f days\n",
              min(tdf_infected$time_infection),
              max(tdf_infected$time_infection)))
  cat(sprintf("  Max generation    : %d\n",
              max(tdf_infected$generation, na.rm = TRUE)))
  
  ct <- table(ct_labels[as.character(tdf_infected$contact_type)], useNA = "ifany")
  cat("\n  Infection route breakdown:\n")
  for (nm in names(ct))
    cat(sprintf("    %-15s : %4d  (%.1f%%)\n", nm, ct[nm], ct[nm] / n_infected * 100))
  
  n_asymp <- sum(tdf_infected$asymptomatic == 1L, na.rm = TRUE)
  cat(sprintf("\n  Symptomatic       : %d (%.1f%%)\n",
              n_infected - n_asymp, (n_infected - n_asymp) / n_infected * 100))
  cat(sprintf("  Asymptomatic      : %d (%.1f%%)\n",
              n_asymp, n_asymp / n_infected * 100))
  
  n_treated <- sum(tdf_infected$treated == 1L, na.rm = TRUE)
  cat(sprintf("\n  Treated (any)     : %d (%.1f%%)\n",
              n_treated, n_treated / n_infected * 100))
  tr <- table(tr_labels[as.character(tdf_infected$treat_reason)], useNA = "ifany")
  for (nm in names(tr))
    cat(sprintf("    %-15s : %4d\n", ifelse(is.na(nm), "untreated", nm), tr[nm]))
  
  hh_ar <- tdf_full %>%
    group_by(household) %>%
    summarise(
      hh_size    = n(),
      n_infected = sum(status %in% c(2L, 3L)),
      .groups    = "drop"
    ) %>%
    filter(n_infected > 0) %>%
    mutate(ar = n_infected / hh_size)
  
  cat(sprintf("\n  Households with >=1 infection : %d / %d\n",
              nrow(hh_ar), net_data$n_households))
  cat(sprintf("  Mean household AR             : %.1f%%\n",
              mean(hh_ar$ar) * 100))
  cat(sprintf("  Fully infected households     : %d\n",
              sum(hh_ar$ar == 1.0)))
}
cat(strrep("-", 62), "\n\n")


# ==============================================================================
# [Section 5] Visualization
# ==============================================================================

if (n_infected == 0) {
  cat("No infections — skipping visualization.\n")
} else {
  
  # --------------------------------------------------------------------------
  # 5-1. Epidemic Curve
  # Stacked bar of daily new infections by route.
  # Red dashed line marks the antiviral program start day.
  # --------------------------------------------------------------------------
  epi_df <- tdf_infected %>%
    mutate(
      day          = floor(time_infection),
      contact_type = factor(
        ct_labels[as.character(contact_type)],
        levels = c("seed", "household", "community"),
        labels = c("Seed", "Household", "Community")
      )
    ) %>%
    count(day, contact_type, name = "n_new")
  
  p_epi <- ggplot(epi_df, aes(x = day, y = n_new, fill = contact_type)) +
    geom_col(alpha = 0.85, width = 0.85) +
    geom_vline(xintercept = antiviral_start,
               linetype = "dashed", color = "red2", linewidth = 0.8) +
    annotate("text",
             x     = antiviral_start + 0.3,
             y     = max(epi_df$n_new) * 0.95,
             label = sprintf("AV Start\n(Day %d)", antiviral_start),
             hjust = 0, color = "red2", size = 3.5) +
    scale_fill_manual(values = c("Seed"      = "gray50",
                                 "Household" = "steelblue",
                                 "Community" = "tomato")) +
    labs(
      title    = sprintf("Daily New Infections (Epidemic Curve) — %s", pg$label),
      subtitle = sprintf("Total %d / %d | p_hh=%.4f, p_comm=%.4f | AV start day %d | %s",
                         n_infected, n_people, p_inf_household, p_inf_community,
                         antiviral_start,
                         ifelse(is_epidemic, "EPIDEMIC", "Contained")),
      x    = "Day of infection",
      y    = "New infections",
      fill = "Route"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  path_epi <- file.path(figure_dir, "plot_singlerun_epi_curve.png")
  ggsave(path_epi, p_epi, width = 10, height = 5, dpi = 150)
  cat(sprintf("Saved: %s\n", path_epi))
  
  # --------------------------------------------------------------------------
  # 5-2. Transmission Tree
  # x-axis: absolute time of infection
  # y-axis: chronological rank (earliest = bottom)
  # Lines connect each infector to its infectees.
  # Point color = infection route; shape = symptom status.
  # --------------------------------------------------------------------------
  tree_nodes <- tdf_infected %>%
    arrange(time_infection) %>%
    mutate(
      y_pos        = row_number(),
      contact_type = factor(
        ct_labels[as.character(contact_type)],
        levels = c("seed", "household", "community"),
        labels = c("Seed", "Household", "Community")
      ),
      shape_type = ifelse(asymptomatic == 1L, "Asymptomatic", "Symptomatic")
    )
  
  y_lookup <- setNames(tree_nodes$y_pos, tree_nodes$person_id)
  
  tree_edges <- tree_nodes %>%
    filter(!is.na(ancestor_id)) %>%
    mutate(
      x_start = tdf_infected$time_infection[match(ancestor_id, tdf_infected$person_id)],
      y_start = y_lookup[as.character(ancestor_id)],
      x_end   = time_infection,
      y_end   = y_pos
    ) %>%
    filter(!is.na(x_start))
  
  p_tree <- ggplot() +
    geom_segment(
      data = tree_edges,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "gray70", linewidth = 0.3, alpha = 0.6
    ) +
    geom_point(
      data = tree_nodes,
      aes(x = time_infection, y = y_pos,
          color = contact_type, shape = shape_type),
      size = 2.0, alpha = 0.85
    ) +
    geom_vline(xintercept = antiviral_start,
               linetype = "dashed", color = "red2", linewidth = 0.7) +
    scale_color_manual(values = c("Seed"      = "gray40",
                                  "Household" = "steelblue",
                                  "Community" = "tomato")) +
    scale_shape_manual(values = c("Symptomatic" = 16, "Asymptomatic" = 17)) +
    labs(
      title    = sprintf("Transmission Tree — %s", pg$label),
      subtitle = sprintf("Red dashed line = AV start (Day %d) | %s",
                         antiviral_start,
                         ifelse(is_epidemic, "EPIDEMIC", "Contained")),
      x        = "Time of infection (days)",
      y        = "Infected individuals (chronological)",
      color    = "Route",
      shape    = "Symptom status"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold")
    )
  
  path_tree <- file.path(figure_dir, "plot_singlerun_transmission_tree.png")
  ggsave(path_tree, p_tree,
         width  = 10,
         height = max(4, min(14, n_infected / 15)),
         dpi    = 150)
  cat(sprintf("Saved: %s\n", path_tree))
}
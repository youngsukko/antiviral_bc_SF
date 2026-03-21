# ==============================================================================
# step_1_bc_fixedrun.R
#
# Purpose:
#   Run n_sim Monte Carlo replications of the branching process under a single
#   fixed parameter configuration and summarize the distribution of outcomes.
#
# Use this file to: 
#   - Quickly inspect the outcome distribution for one specific scenario
#   - Validate model behavior before launching the full grid (multirun)
#   - Compare two pathogens side-by-side under identical intervention settings
#
# Prerequisites:
#   - scp_build_network_sf.R must have been run first
#   - functions/network_bp_sim.R must exist
#
# Outputs (figures/ folder):
#   plot_fixedrun_outbreak_size.png   : histogram of total infected per rep
#   plot_fixedrun_epidemic_prob.png   : bar showing P(epidemic) vs P(contained)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(parallel)
library(arrow)

source("functions/network_bp_sim.R")


# ==============================================================================
# [Configuration Block] — Modify only this section
# ==============================================================================

# --- Paths --------------------------------------------------------------------
net_data_path <- "output/net_data_sf.rds"
figure_dir    <- "figures"
if (!dir.exists(figure_dir)) dir.create(figure_dir)

# --- Pathogen -----------------------------------------------------------------
# Switch between "H1N1" and "SC2" (SARS-CoV-2) to compare pathogens.
# pathogen <- "H1N1"   # "H1N1" or "SC2"
pathogen <- "SC2"   # "H1N1" or "SC2"

pathogen_params <- list(
  H1N1 = list(
    label             = "H1N1 (2009 pandemic)",
    R0                = 1.5,
    prop_asymptomatic = 0.50,
    # Gamma generation time: mean ~2.6 d (White et al. 2009)
    generation_time_fn    = function(n) rgamma(n, shape = 4, rate = 1.538),
    infectious_before_onset_fn = function(n) runif(n, 0.5, 0.7),
    # Log-normal incubation (Lessler et al. 2009)
    infection_to_onset_fn = function(n)
      rlnorm(n, meanlog = log(1.4), sdlog = log(1.51))
  ),
  SC2 = list(
    label             = "SARS-CoV-2 (ancestral)",
    R0                = 2.79,
    prop_asymptomatic = 0.15,
    # Weibull generation time: mean ~5.0 d, SD ~1.9 d
    # Ferretti et al., Science 2020
    # (https://www.science.org/doi/10.1126/science.abb6936)
    generation_time_fn    = function(n) rweibull(n, shape = 2.826, scale = 5.665),
    infectious_before_onset_fn = function(n) runif(n, 0.4, 0.6),
    # Log-normal incubation: mean ~5.5 d, median ~5.1 d
    # Lauer et al., Ann Intern Med 2020
    # (https://www.acpjournals.org/doi/10.7326/M20-0504)
    infection_to_onset_fn = function(n) rlnorm(n, meanlog = 1.621, sdlog = 0.418)
  )
)

pg <- pathogen_params[[pathogen]]

# --- Network degree parameters ------------------------------------------------
hh    <- 1
nonhh <- 0.1

# --- Derived risk -------------------------------------------------------------
R0              <- pg$R0
p_inf_household <- R0 / (2.14 * hh + 9.16 * nonhh)
p_inf_community <- p_inf_household * nonhh

prop_asymptomatic     <- pg$prop_asymptomatic
generation_time_fn    <- pg$generation_time_fn
infection_to_onset_fn <- pg$infection_to_onset_fn
infectious_before_onset_fn <- pg$infectious_before_onset_fn  

# --- Antiviral efficacy — community infection prevention ----------------------
antiviral_eff_inf_e_max   <- 0.50
antiviral_eff_inf_kappa   <- 7.51
antiviral_eff_inf_tau     <- 1.34

# --- Antiviral efficacy — household infection prevention ----------------------
antiviral_eff_inf_hh_e_max  <- 0.63
antiviral_eff_inf_hh_kappa  <- 6.06
antiviral_eff_inf_hh_tau    <- 0.89

# --- Antiviral efficacy — transmission reduction (source side) ----------------
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
# Adjust these to the scenario you want to inspect.
antiviral_start           <- 28    # day AV program begins (0 = immediate)

prob_treat_self           <- 0.9
prob_treat_household      <- 0.9
prob_treat_community      <- 0.9

prob_quarantine_self      <- 0.0
prob_quarantine_household <- 0.0
prob_quarantine_community <- 0.0

quarantine_efficacy <- 0.9

logistical_delay_fn <- function(n) runif(n, 0, 2)   # Uniform(0, 2 days)
quarantine_fn       <- function(n) runif(n, 0, 2)

# --- Simulation control -------------------------------------------------------
n_sim                   <- 500
epidemic_prob_threshold <- 0.999
seeding_cases           <- 1
# n_cores                 <- max(1L, detectCores() - 1L)
n_cores                 <- 10


# ==============================================================================
# [Section 1] Load network
# ==============================================================================

cat(strrep("=", 62), "\n")
cat(" SF Branching Process — Fixed-Setting Monte Carlo\n")
cat(strrep("=", 62), "\n\n")

if (!file.exists(net_data_path))
  stop(sprintf("Network file not found: %s\nRun scp_build_network_sf.R first.", net_data_path))

net_data <- readRDS(net_data_path)
n_people <- nrow(net_data$nodes)

epidemic_queue_threshold <- if (R0 > 1)
  ceiling(log(1 - epidemic_prob_threshold) / log(1 / R0)) else Inf

cat(sprintf("Network         : %s\n",   net_data$network_type))
cat(sprintf("Population      : %d\n",   n_people))
cat(sprintf("Pathogen        : %s\n",   pg$label))
cat(sprintf("R0              : %.2f\n", R0))
cat(sprintf("p_inf_household : %.4f\n", p_inf_household))
cat(sprintf("p_inf_community : %.4f\n", p_inf_community))
cat(sprintf("prop_asymp      : %.2f\n", prop_asymptomatic))
cat(sprintf("AV start day    : %d\n",   antiviral_start))
cat(sprintf("Epidemic cutoff : %d asymptomatic active (Whittle @ %.1f%%)\n",
            epidemic_queue_threshold, epidemic_prob_threshold * 100))
cat(sprintf("Replications    : %d\n",   n_sim))
cat(sprintf("Cores           : %d\n\n", n_cores))

cat("Precomputing neighbor lists...\n")
neighbors_by_type <- precompute_neighbors(net_data)
cat("Done.\n\n")


# ==============================================================================
# [Section 2] Run n_sim replications in parallel
# ==============================================================================

cat(sprintf("Running %d replications...\n", n_sim))
t_start <- proc.time()

cl <- makeCluster(n_cores)

clusterExport(cl, varlist = c(
  "network_bp_sim", "precompute_neighbors",
  "neighbors_by_type", "net_data",
  "R0", "epidemic_prob_threshold",
  "p_inf_household", "p_inf_community",
  "prop_asymptomatic",
  "generation_time_fn", "infection_to_onset_fn", "infectious_before_onset_fn",
  "antiviral_start",
  "logistical_delay_fn", "quarantine_fn",
  "drug_eff_inf_fn", "drug_eff_inf_hh_fn", "drug_eff_trans_fn",
  "antiviral_eff_inf_e_max",   "antiviral_eff_inf_kappa",   "antiviral_eff_inf_tau",     # 추가
  "antiviral_eff_inf_hh_e_max","antiviral_eff_inf_hh_kappa","antiviral_eff_inf_hh_tau",  # 추가
  "antiviral_eff_trans_e_max", "antiviral_eff_trans_kappa", "antiviral_eff_trans_tau",   # 추가
  "prob_treat_self", "prob_treat_household", "prob_treat_community",
  "prob_quarantine_self", "prob_quarantine_household", "prob_quarantine_community",
  "quarantine_efficacy", "seeding_cases"
), envir = environment())

results_list <- parLapply(cl, seq_len(n_sim), function(rep_i) {
  
  res <- network_bp_sim(
    neighbors_by_type         = neighbors_by_type,
    node_df                   = net_data$nodes,
    R0                        = R0,
    epidemic_prob_threshold   = epidemic_prob_threshold,
    generation_time_fn        = generation_time_fn,
    infection_to_onset_fn     = infection_to_onset_fn,
    infectious_before_onset_fn= infectious_before_onset_fn,
    prop_asymptomatic         = prop_asymptomatic,
    p_inf_household           = p_inf_household,
    p_inf_community           = p_inf_community,
    antiviral_start           = antiviral_start,
    logistical_delay_fn       = logistical_delay_fn,
    drug_eff_inf_fn           = drug_eff_inf_fn,
    drug_eff_inf_hh_fn        = drug_eff_inf_hh_fn,
    drug_eff_trans_fn         = drug_eff_trans_fn,
    eff_inf_e_max             = antiviral_eff_inf_e_max,
    eff_inf_hh_e_max          = antiviral_eff_inf_hh_e_max,
    eff_trans_e_max           = antiviral_eff_trans_e_max,
    prob_treat_self           = prob_treat_self,
    prob_treat_household      = prob_treat_household,
    prob_treat_community      = prob_treat_community,
    time_to_quarantine_fn     = quarantine_fn,
    prob_quarantine_self      = prob_quarantine_self,
    prob_quarantine_household = prob_quarantine_household,
    prob_quarantine_community = prob_quarantine_community,
    quarantine_efficacy       = quarantine_efficacy,
    seeding_cases             = seeding_cases,
    t0                        = 0,
    seed                      = rep_i
  )
  
  data.frame(
    rep           = rep_i,
    is_epidemic   = res$is_epidemic,
    n_infected    = nrow(res$infected),
    active_at_end = res$active_at_end
  )
})

stopCluster(cl)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("Done. Elapsed: %.1f sec\n\n", elapsed))

results_df <- do.call(rbind, results_list)


# ==============================================================================
# [Section 3] Console summary
# ==============================================================================

n_epidemic  <- sum(results_df$is_epidemic)
n_contained <- n_sim - n_epidemic
p_epidemic  <- n_epidemic / n_sim

contained_df <- results_df[!results_df$is_epidemic, ]
epidemic_df  <- results_df[results_df$is_epidemic, ]

cat(strrep("-", 62), "\n")
cat(" Summary\n")
cat(strrep("-", 62), "\n")
cat(sprintf("  Replications      : %d\n", n_sim))
cat(sprintf("  Epidemic          : %d (%.1f%%)\n", n_epidemic,  p_epidemic * 100))
cat(sprintf("  Contained         : %d (%.1f%%)\n", n_contained, (1 - p_epidemic) * 100))

cat(sprintf("\n  --- All runs (n=%d) ---\n", n_sim))
cat(sprintf("  Mean infected     : %.1f\n", mean(results_df$n_infected)))
cat(sprintf("  Median infected   : %.1f\n", median(results_df$n_infected)))
cat(sprintf("  SD infected       : %.1f\n", sd(results_df$n_infected)))
cat(sprintf("  [Q25, Q75]        : [%.1f, %.1f]\n",
            quantile(results_df$n_infected, 0.25),
            quantile(results_df$n_infected, 0.75)))
cat(sprintf("  [P5, P95]         : [%.1f, %.1f]\n",
            quantile(results_df$n_infected, 0.05),
            quantile(results_df$n_infected, 0.95)))

if (nrow(contained_df) > 0) {
  cat(sprintf("\n  --- Contained runs (n=%d) ---\n", nrow(contained_df)))
  cat(sprintf("  Mean infected     : %.1f\n", mean(contained_df$n_infected)))
  cat(sprintf("  Median infected   : %.1f\n", median(contained_df$n_infected)))
  cat(sprintf("  [Q25, Q75]        : [%.1f, %.1f]\n",
              quantile(contained_df$n_infected, 0.25),
              quantile(contained_df$n_infected, 0.75)))
}

if (nrow(epidemic_df) > 0) {
  cat(sprintf("\n  --- Epidemic runs (n=%d) ---\n", nrow(epidemic_df)))
  cat(sprintf("  Mean infected     : %.1f\n", mean(epidemic_df$n_infected)))
  cat(sprintf("  Median infected   : %.1f\n", median(epidemic_df$n_infected)))
  cat(sprintf("  [Q25, Q75]        : [%.1f, %.1f]\n",
              quantile(epidemic_df$n_infected, 0.25),
              quantile(epidemic_df$n_infected, 0.75)))
}

cat(strrep("-", 62), "\n\n")


# ==============================================================================
# [Section 4] Visualization
# ==============================================================================

results_df$outcome <- ifelse(results_df$is_epidemic, "Epidemic", "Contained")

# --------------------------------------------------------------------------
# 4-1. Outbreak size distribution (histogram)
# Split by outcome; epidemic runs stacked/overlaid in a different color.
# Log scale on x-axis because contained runs cluster near 0.
# --------------------------------------------------------------------------
p_hist <- ggplot(results_df, aes(x = n_infected, fill = outcome)) +
  geom_histogram(bins = 60, alpha = 0.80, position = "stack", color = NA) +
  scale_fill_manual(values = c("Contained" = "steelblue", "Epidemic" = "tomato")) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title    = sprintf("Outbreak Size Distribution — %s", pg$label),
    subtitle = sprintf(
      "%d reps | AV start day %d | P(epidemic) = %.1f%% | R0 = %.2f",
      n_sim, antiviral_start, p_epidemic * 100, R0
    ),
    x    = "Total infected per replication",
    y    = "Count",
    fill = "Outcome"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

path_hist <- file.path(figure_dir, "plot_fixedrun_outbreak_size.png")
ggsave(path_hist, p_hist, width = 10, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", path_hist))

# --------------------------------------------------------------------------
# 4-2. Contained-only outbreak size (zoomed in)
# Epidemic runs are excluded so the small-outbreak tail is visible.
# --------------------------------------------------------------------------
if (nrow(contained_df) > 1) {
  p_contained <- ggplot(contained_df, aes(x = n_infected)) +
    geom_histogram(bins = 40, fill = "steelblue", alpha = 0.85, color = NA) +
    labs(
      title    = sprintf("Contained Outbreak Size — %s", pg$label),
      subtitle = sprintf(
        "Contained runs only (n=%d of %d) | AV start day %d | R0 = %.2f",
        nrow(contained_df), n_sim, antiviral_start, R0
      ),
      x = "Total infected per replication",
      y = "Count"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  path_contained <- file.path(figure_dir, "plot_fixedrun_contained_size.png")
  ggsave(path_contained, p_contained, width = 8, height = 4, dpi = 150)
  cat(sprintf("Saved: %s\n", path_contained))
}

# --------------------------------------------------------------------------
# 4-3. Outcome bar: P(epidemic) vs P(contained)
# --------------------------------------------------------------------------
outcome_df <- data.frame(
  outcome     = c("Epidemic", "Contained"),
  probability = c(p_epidemic, 1 - p_epidemic)
)

p_bar <- ggplot(outcome_df, aes(x = outcome, y = probability, fill = outcome)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", probability * 100)),
            vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Contained" = "steelblue", "Epidemic" = "tomato")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08)) +
  labs(
    title    = sprintf("Epidemic Probability — %s", pg$label),
    subtitle = sprintf(
      "%d reps | AV start day %d | R0 = %.2f | Whittle threshold = %d asymptomatic",
      n_sim, antiviral_start, R0, epidemic_queue_threshold
    ),
    x    = NULL,
    y    = "Probability"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

path_bar <- file.path(figure_dir, "plot_fixedrun_epidemic_prob.png")
ggsave(path_bar, p_bar, width = 6, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", path_bar))

cat("\n")
cat(strrep("=", 62), "\n")
cat(" Done\n")
cat(strrep("=", 62), "\n\n")
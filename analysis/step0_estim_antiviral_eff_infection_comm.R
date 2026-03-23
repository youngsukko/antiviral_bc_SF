################################################################################
# Logistic Antiviral Efficacy Model — Infection Prevention (Community)
# Bayesian estimation via Stan (HMC/NUTS)
#
# Model: E(t) = Emax / (1 + exp(kappa * (t - tau)))
# Emax is fixed; kappa and tau are freely estimated.
################################################################################

library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ============================================================
# *** USER CONFIGURATION ***
# ============================================================

# [1] Observed data
obs <- data.frame(
  a     = c(-10,  2),
  b     = c(0,    NA),
  E_obs = c(0.5,  0.00),
  type  = c("interval", "point"),
  stringsAsFactors = FALSE
)

# [2] Fixed parameter
EMAX_FIXED  <- 0.5   # Emax is fixed, not estimated

# [3] Parameter bounds for free parameters (uniform prior support)
KAPPA_LOWER <- 0.01
KAPPA_UPPER <- 10.0
TAU_LOWER   <- -5.0
TAU_UPPER   <- 10.0

# [4] Likelihood measurement error SD
SIGMA <- 0.03

# [5] MCMC settings
N_ITER   <- 4000
BURN_IN  <- 2000
CHAINS   <- 4
SEED     <- 42

# [6] Plot range
PLOT_T_RANGE <- c(-2, 4)
PLOT_Y_RANGE <- c(0, 0.6)

# ============================================================
# 1. Stan Model
# ============================================================

stan_code <- "
data {
  int<lower=0> N_int;
  vector[N_int] a_int;
  vector[N_int] b_int;
  vector[N_int] E_obs_int;

  int<lower=0> N_point;
  vector[N_point] t_point;
  vector[N_point] E_obs_point;

  real<lower=0> sigma;
  real Emax;              // fixed — passed from R, not sampled

  real kappa_lower; real kappa_upper;
  real tau_lower;   real tau_upper;
}
parameters {
  real<lower=kappa_lower, upper=kappa_upper> kappa;
  real<lower=tau_lower,   upper=tau_upper>   tau;
}
model {
  // Interval-averaged observations (analytic integral)
  for (i in 1:N_int) {
    real Fb = Emax * (b_int[i] - log1p_exp(kappa * (b_int[i] - tau)) / kappa);
    real Fa = Emax * (a_int[i] - log1p_exp(kappa * (a_int[i] - tau)) / kappa);
    real E_pred = (Fb - Fa) / (b_int[i] - a_int[i]);
    E_obs_int[i] ~ normal(E_pred, sigma);
  }
  // Instantaneous point observations
  for (i in 1:N_point) {
    real E_pred = Emax / (1 + exp(kappa * (t_point[i] - tau)));
    E_obs_point[i] ~ normal(E_pred, sigma);
  }
}
"

# ============================================================
# 2. Prepare Data
# ============================================================
set.seed(SEED)

is_int   <- obs$type == "interval"
is_point <- !is_int

stan_data <- list(
  N_int       = sum(is_int),
  N_point     = sum(is_point),
  a_int       = array(obs$a[is_int],      dim = sum(is_int)),
  b_int       = array(obs$b[is_int],      dim = sum(is_int)),
  E_obs_int   = array(obs$E_obs[is_int],  dim = sum(is_int)),
  t_point     = array(obs$a[is_point],    dim = sum(is_point)),
  E_obs_point = array(obs$E_obs[is_point],dim = sum(is_point)),
  sigma       = SIGMA,
  Emax        = EMAX_FIXED,
  kappa_lower = KAPPA_LOWER, kappa_upper = KAPPA_UPPER,
  tau_lower   = TAU_LOWER,   tau_upper   = TAU_UPPER
)

# ============================================================
# 3. Fit (with caching)
# ============================================================
OUTPUT_DIR <- "output"
CACHE_FILE <- file.path(OUTPUT_DIR, "antiviral_eff_estim_infection_comm.rds")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

if (file.exists(CACHE_FILE)) {
  cat("--- Cached samples found. Loading from:", CACHE_FILE, "---\n")
  fit <- readRDS(CACHE_FILE)
} else {
  cat("--- No cache found. Compiling and running Stan (HMC/NUTS) ---\n")
  fit <- stan(
    model_code = stan_code,
    data       = stan_data,
    iter       = N_ITER,
    warmup     = BURN_IN,
    chains     = CHAINS,
    seed       = SEED,
    refresh    = max(1, N_ITER / 10)
  )
  saveRDS(fit, CACHE_FILE)
  cat("--- Posterior samples saved to:", CACHE_FILE, "---\n")
}

# ============================================================
# 4. Post-processing
# ============================================================
post       <- rstan::extract(fit)
chain_full <- cbind(EMAX_FIXED, post$kappa, post$tau)
colnames(chain_full) <- c("Emax", "kappa", "tau")

post_mean <- colMeans(chain_full)
post_sd   <- apply(chain_full, 2, sd)
n_post    <- nrow(chain_full)

# ============================================================
# 5. Posterior Summary
# ============================================================
cat("\n=============================================================\n")
cat(sprintf(" Posterior Summary  (Emax fixed = %.2f)\n", EMAX_FIXED))
cat("=============================================================\n")
print(fit, pars = c("kappa", "tau"))

# Utility functions
E_instant <- function(t, Emax, kappa, tau) {
  Emax / (1 + exp(kappa * (t - tau)))
}
E_interval_avg <- function(a, b, Emax, kappa, tau) {
  Fb <- Emax * (b - log1p(exp(kappa * (b - tau))) / kappa)
  Fa <- Emax * (a - log1p(exp(kappa * (a - tau))) / kappa)
  (Fb - Fa) / (b - a)
}
E_predicted <- function(obs_row, Emax, kappa, tau) {
  if (obs_row$type == "interval") {
    E_interval_avg(obs_row$a, obs_row$b, Emax, kappa, tau)
  } else {
    E_instant(obs_row$a, Emax, kappa, tau)
  }
}

cat("\n--- Goodness of fit (posterior mean) ---\n")
for (i in seq_len(nrow(obs))) {
  E_pred <- E_predicted(obs[i, ], post_mean["Emax"], post_mean["kappa"], post_mean["tau"])
  if (obs$type[i] == "interval") {
    cat(sprintf("  interval [%g,%g]: obs=%.4f  pred=%.4f  resid=%+.6f\n",
                obs$a[i], obs$b[i], obs$E_obs[i], E_pred, obs$E_obs[i] - E_pred))
  } else {
    cat(sprintf("  point    t=%g  : obs=%.4f  pred=%.4f  resid=%+.6f\n",
                obs$a[i], obs$E_obs[i], E_pred, obs$E_obs[i] - E_pred))
  }
}

# ============================================================
# 6. Plots
# ============================================================
cat("\n--- Generating plots ---\n")
if (!dir.exists("figures")) dir.create("figures")

ts <- seq(PLOT_T_RANGE[1], PLOT_T_RANGE[2], length.out = 300)

# Pre-compute posterior predictive band
E_matrix <- sapply(seq_len(n_post), function(k)
  E_instant(ts, chain_full[k, "Emax"], chain_full[k, "kappa"], chain_full[k, "tau"]))
cilo <- apply(E_matrix, 1, quantile, 0.025)
cihi <- apply(E_matrix, 1, quantile, 0.975)

# Spaghetti sample
idx_draw <- sample(n_post, min(300, n_post))
spaghetti_df <- do.call(rbind, lapply(idx_draw, function(k) {
  data.frame(
    t    = ts,
    E    = E_instant(ts, chain_full[k, "Emax"], chain_full[k, "kappa"], chain_full[k, "tau"]),
    draw = k
  )
}))

# Posterior mean + CI band
fit_df <- data.frame(
  t    = ts,
  E    = E_instant(ts, post_mean["Emax"], post_mean["kappa"], post_mean["tau"]),
  cilo = cilo,
  cihi = cihi
)

# Interval observations: clip a to plot range for display
obs_interval <- obs %>%
  filter(type == "interval") %>%
  mutate(a_clip = pmax(a, PLOT_T_RANGE[1]))   # clip -10 to -2

# Point observations
obs_point <- obs %>%
  filter(type == "point") %>%
  mutate(t_plot = a)

# --------------------------------------------------------------
# Figure 1 — Posterior fit
# --------------------------------------------------------------
p_fit <- ggplot() +
  geom_line(data = spaghetti_df,
            aes(x = t, y = E, group = draw),
            color = "gray60", alpha = 0.04, linewidth = 0.3) +
  geom_ribbon(data = fit_df,
              aes(x = t, ymin = cilo, ymax = cihi),
              fill = "darkgreen", alpha = 0.15) +
  geom_line(data = fit_df, aes(x = t, y = cilo),
            color = "darkgreen", linewidth = 0.5, linetype = "dashed") +
  geom_line(data = fit_df, aes(x = t, y = cihi),
            color = "darkgreen", linewidth = 0.5, linetype = "dashed") +
  geom_line(data = fit_df, aes(x = t, y = E),
            color = "darkgreen", linewidth = 1.0) +
  # Interval observation: dashed red horizontal line clipped to plot range
  geom_segment(data = obs_interval,
               aes(x = a_clip, xend = b, y = E_obs, yend = E_obs),
               color = "red", linewidth = 1.0, linetype = "dashed") +
  # Point observation: filled circle
  geom_point(data = obs_point,
             aes(x = t_plot, y = E_obs),
             shape = 16, color = "red", size = 3) +
  scale_y_continuous(limits = PLOT_Y_RANGE) +
  scale_x_continuous(limits = PLOT_T_RANGE) +
  labs(x = "Time since index's onset (days)",
       y = "Infection prevention efficacy E(t)") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(panel.grid.minor = element_blank())

ggsave("figures/antiviral_infection_comm_fit.png",
       p_fit, width = 5, height = 4, dpi = 500)
cat("Saved: figures/antiviral_infection_comm_fit.png\n")

# --------------------------------------------------------------
# Figure 2 — Diagnostics (trace + marginal posteriors, kappa & tau only)
# --------------------------------------------------------------
chain_free <- as.data.frame(chain_full[, c("kappa", "tau")]) %>%
  mutate(sample = seq_len(n()))

param_levels <- c("kappa", "tau")

long_df <- chain_free %>%
  pivot_longer(all_of(param_levels), names_to = "param", values_to = "value") %>%
  mutate(param = factor(param, levels = param_levels))

post_mean_df <- data.frame(
  param = factor(param_levels, levels = param_levels),
  mean  = post_mean[param_levels]
)

ci_df <- data.frame(
  param = factor(param_levels, levels = param_levels),
  lo    = sapply(param_levels, function(p) quantile(chain_full[, p], 0.025)),
  hi    = sapply(param_levels, function(p) quantile(chain_full[, p], 0.975))
)

# Trace plot
p_trace <- ggplot(long_df, aes(x = sample, y = value)) +
  geom_line(color = "gray40", alpha = 0.4, linewidth = 0.2) +
  geom_hline(data = post_mean_df, aes(yintercept = mean),
             color = "blue", linewidth = 0.7) +
  facet_wrap(~ param, ncol = 2, scales = "free_y") +
  labs(x = "Sample", y = "Value") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

# Marginal posterior histograms
p_hist <- ggplot(long_df, aes(x = value)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(data = post_mean_df, aes(xintercept = mean),
             color = "red", linewidth = 0.8) +
  geom_vline(data = ci_df, aes(xintercept = lo),
             color = "red", linewidth = 0.5, linetype = "dashed") +
  geom_vline(data = ci_df, aes(xintercept = hi),
             color = "red", linewidth = 0.5, linetype = "dashed") +
  facet_wrap(~ param, ncol = 2, scales = "free") +
  labs(x = "Value", y = "Count") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

p_diag <- p_trace / p_hist

ggsave("figures/antiviral_infection_comm_diagnostics.png",
       p_diag, width = 8, height = 8, dpi = 500)
cat("Saved: figures/antiviral_infection_comm_diagnostics.png\n")

cat("=============================================================\n")
cat(" Finished. Check the 'figures' folder.\n")
cat("=============================================================\n")
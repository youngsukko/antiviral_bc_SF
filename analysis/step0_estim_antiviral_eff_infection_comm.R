################################################################################
# Logistic Antiviral Efficacy Model — Infection Prevention
# Bayesian estimation via Stan (Hamiltonian Monte Carlo / NUTS)
#
# Model: E(t) = Emax / (1 + exp(kappa * (t - tau)))
# * Closed-form analytic integral used for interval data (no numerical integration)
################################################################################
# install.packages("rstan")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # Enable parallel processing

# ============================================================
# *** USER CONFIGURATION ***
# ============================================================

# [1] Observed data
obs <- data.frame(
  a     = c(0,    2),
  b     = c(1,    NA),
  E_obs = c(0.5, 0.00),
  type  = c("interval", "point"),
  stringsAsFactors = FALSE
)

# [2] Fix parameters?
FIX_KAPPA   <- FALSE
KAPPA_FIXED <- 1.0

FIX_EMAX    <- TRUE
EMAX_FIXED  <- 0.5

# [3] Parameter bounds (used as uniform prior support)
EMAX_LOWER  <- 0.01
EMAX_UPPER  <- 0.80
KAPPA_LOWER <- 0.01
KAPPA_UPPER <- 10.0
TAU_LOWER   <- -5.0
TAU_UPPER   <- 10.0

# [4] Likelihood: measurement error SD
SIGMA <- 0.03

# [5] MCMC settings (Stan is highly efficient, so fewer iterations are needed)
N_ITER    <- 4000
BURN_IN   <- 2000
CHAINS    <- 4
SEED      <- 42

# [6] Plot range
PLOT_T_RANGE <- c(0, 4)
PLOT_Y_RANGE <- c(0, 0.9)

# (Note: INIT_FREE and PROPOSAL_SD are no longer needed — Stan handles this automatically.)

# ============================================================
# 1. Stan Model Definition (Inline)
# ============================================================

stan_code <- "
data {
  // Interval data
  int<lower=0> N_int;
  vector[N_int] a_int;
  vector[N_int] b_int;
  vector[N_int] E_obs_int;
  
  // Point data
  int<lower=0> N_point;
  vector[N_point] t_point;
  vector[N_point] E_obs_point;

  real<lower=0> sigma;

  // Fixed flags & values
  int<lower=0, upper=1> fix_Emax;
  int<lower=0, upper=1> fix_kappa;
  real Emax_fixed;
  real kappa_fixed;

  // Bounds
  real Emax_lower; real Emax_upper;
  real kappa_lower; real kappa_upper;
  real tau_lower; real tau_upper;
}
parameters {
  // Conditionally active parameters (size 0 if fixed, size 1 if free)
  array[fix_Emax ? 0 : 1] real<lower=Emax_lower, upper=Emax_upper> Emax_free;
  array[fix_kappa ? 0 : 1] real<lower=kappa_lower, upper=kappa_upper> kappa_free;
  
  // tau is always free
  real<lower=tau_lower, upper=tau_upper> tau;
}
transformed parameters {
  // Map fixed/free values to the actual parameters used in the model
  real Emax = fix_Emax ? Emax_fixed : Emax_free[1];
  real kappa = fix_kappa ? kappa_fixed : kappa_free[1];
}
model {
  // Priors: Implicitly uniform over the defined lower/upper bounds

  // Likelihood for interval-averaged data (using analytic CDF trick)
  for (i in 1:N_int) {
    // F(t) = Emax * (t - log1p_exp(kappa * (t - tau)) / kappa)
    real Fb = Emax * (b_int[i] - log1p_exp(kappa * (b_int[i] - tau)) / kappa);
    real Fa = Emax * (a_int[i] - log1p_exp(kappa * (a_int[i] - tau)) / kappa);
    real E_pred = (Fb - Fa) / (b_int[i] - a_int[i]);
    
    E_obs_int[i] ~ normal(E_pred, sigma);
  }

  // Likelihood for instantaneous point data
  for (i in 1:N_point) {
    real E_pred = Emax / (1 + exp(kappa * (t_point[i] - tau)));
    E_obs_point[i] ~ normal(E_pred, sigma);
  }
}
"

# ============================================================
# 2. Prepare Data & Fit Stan Model
# ============================================================
set.seed(SEED)

is_int <- obs$type == "interval"
stan_data <- list(
  N_int       = sum(is_int),
  N_point     = sum(!is_int),
  
  # Force 1-D arrays to prevent R from passing single values as scalars
  a_int       = array(obs$a[is_int], dim = sum(is_int)),
  b_int       = array(obs$b[is_int], dim = sum(is_int)),
  E_obs_int   = array(obs$E_obs[is_int], dim = sum(is_int)),
  
  t_point     = array(obs$a[!is_int], dim = sum(!is_int)),
  E_obs_point = array(obs$E_obs[!is_int], dim = sum(!is_int)),
  
  sigma       = SIGMA,
  
  fix_Emax    = as.integer(FIX_EMAX),
  fix_kappa   = as.integer(FIX_KAPPA),
  Emax_fixed  = EMAX_FIXED,
  kappa_fixed = KAPPA_FIXED,
  
  Emax_lower  = EMAX_LOWER,  Emax_upper  = EMAX_UPPER,
  kappa_lower = KAPPA_LOWER, kappa_upper = KAPPA_UPPER,
  tau_lower   = TAU_LOWER,   tau_upper   = TAU_UPPER
)

# ============================================================
# Cache settings — load existing samples if available,
# otherwise run Stan and save the result for future runs.
# ============================================================
OUTPUT_DIR  <- "output"
CACHE_FILE  <- file.path(OUTPUT_DIR, "antiviral_eff_estim_infection.rds")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

if (file.exists(CACHE_FILE)) {
  cat("--- Cached samples found. Loading from:", CACHE_FILE, "---\n")
  fit <- readRDS(CACHE_FILE)
} else {
  cat("--- No cache found. Compiling and Running Stan (HMC/NUTS) ---\n")
  fit <- stan(
    model_code = stan_code,
    data       = stan_data,
    iter       = N_ITER,
    warmup     = BURN_IN,
    chains     = CHAINS,
    seed       = SEED,
    refresh    = max(1, N_ITER / 10)
  )
  saveRDS(fit, file = CACHE_FILE)
  cat("--- Posterior samples saved to:", CACHE_FILE, "---\n")
}

# ============================================================
# 3. Post-processing (Mimicking original structures for plots)
# ============================================================

# Extract posterior samples
post <- rstan::extract(fit)
n_post <- length(post$tau)

# Build chain_full matrix
chain_full <- cbind(post$Emax, post$kappa, post$tau)
colnames(chain_full) <- c("Emax", "kappa", "tau")

post_mean <- colMeans(chain_full)
post_sd   <- apply(chain_full, 2, sd)

# Reconstruct 'free' parameter mapping for trace plots
free_idx   <- c()
free_names <- c()
if (!FIX_EMAX)  { free_idx <- c(free_idx, 1); free_names <- c(free_names, "Emax") }
if (!FIX_KAPPA) { free_idx <- c(free_idx, 2); free_names <- c(free_names, "kappa") }
free_idx   <- c(free_idx, 3); free_names <- c(free_names, "tau")
n_free     <- length(free_names)

# Build raw chain matrix (free parameters only)
chain <- chain_full[, free_idx, drop=FALSE]
colnames(chain) <- free_names

# ============================================================
# 4. Posterior Summary & Goodness of Fit
# ============================================================

cat("\n=============================================================\n")
cat(" Posterior Summary — Stan Fit\n")
cat("=============================================================\n")
print(fit, pars=c("Emax", "kappa", "tau"))

# Utility functions for predictions (used in GOF check and plots)
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
for (i in 1:nrow(obs)) {
  E_pred <- E_predicted(obs[i, ], post_mean[1], post_mean[2], post_mean[3])
  if (obs$type[i] == "interval") {
    cat(sprintf("  interval [%g,%g]: obs=%.4f  pred=%.4f  resid=%+.6f\n",
                obs$a[i], obs$b[i], obs$E_obs[i], E_pred, obs$E_obs[i] - E_pred))
  } else {
    cat(sprintf("  point    t=%g  : obs=%.4f  pred=%.4f  resid=%+.6f\n",
                obs$a[i], obs$E_obs[i], E_pred, obs$E_obs[i] - E_pred))
  }
}

# ============================================================
# 5. Plots (100% compatible with your existing layout)
# ============================================================

cat("\n--- Generating plots ---\n")

# Build filename suffix based on which parameters are fixed
fix_tag <- ""
if (FIX_KAPPA && FIX_EMAX) {
  fix_tag <- "_both_fixed"
} else if (FIX_KAPPA) {
  fix_tag <- "_kappa_fixed"
} else if (FIX_EMAX) {
  fix_tag <- "_Emax_fixed"
}

for (dtype in c("pdf", "png")) {
  # Create output directory if it does not already exist
  if(!dir.exists("figures")) dir.create("figures")
  
  fname <- paste0("figures/antiviral_infection_prevention_stan", fix_tag, ".", dtype)
  if (dtype == "pdf") pdf(fname, width = 12, height = 10)
  else png(fname, width = 1200, height = 1000, res = 120)
  
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE))
  par(mar = c(4.5, 4.5, 3, 1))
  
  ts <- seq(PLOT_T_RANGE[1], PLOT_T_RANGE[2], length.out = 300)
  
  # Build title suffix showing any fixed parameter values
  fix_str <- c()
  if (FIX_EMAX)  fix_str <- c(fix_str, sprintf("Emax=%.2f", EMAX_FIXED))
  if (FIX_KAPPA) fix_str <- c(fix_str, sprintf("kappa=%.2f", KAPPA_FIXED))
  title_suffix <- if (length(fix_str) > 0)
    paste0(" (", paste(fix_str, collapse = ", "), " fixed)") else ""
  
  # (a) Posterior predictive curves + data
  Ec <- E_instant(ts, post_mean[1], post_mean[2], post_mean[3])
  
  plot(ts, Ec, type = "n",
       xlab = "Time since index's onset (t, days)",
       ylab = "Infection prevention efficacy E(t)",
       main = paste0("(a) Posterior fit (Stan)", title_suffix),
       ylim = PLOT_Y_RANGE)
  
  # Draw a random subset of posterior curves as a spaghetti plot
  idx_draw <- sample(n_post, min(300, n_post))
  for (k in idx_draw) {
    Ek <- E_instant(ts, chain_full[k, 1], chain_full[k, 2], chain_full[k, 3])
    lines(ts, Ek, col = rgb(0.5, 0.5, 0.5, 0.03))
  }
  
  # Compute 95% credible interval band across all posterior samples
  E_matrix <- matrix(NA, n_post, length(ts))
  for (k in 1:n_post)
    E_matrix[k, ] <- E_instant(ts, chain_full[k, 1], chain_full[k, 2], chain_full[k, 3])
  cilo <- apply(E_matrix, 2, quantile, 0.025)
  cihi <- apply(E_matrix, 2, quantile, 0.975)
  polygon(c(ts, rev(ts)), c(cilo, rev(cihi)),
          col = rgb(0, 0.5, 0, 0.15), border = NA)
  lines(ts, Ec, lwd = 2, col = "darkgreen")
  lines(ts, cilo, lwd = 1, col = "darkgreen", lty = 2)
  lines(ts, cihi, lwd = 1, col = "darkgreen", lty = 2)
  
  # Overlay observed data points (rectangles for interval, triangles for point)
  for (i in 1:nrow(obs)) {
    if (obs$type[i] == "interval") {
      rect(obs$a[i], 0, obs$b[i], obs$E_obs[i],
           col = rgb(1, 0, 0, 0.2), border = "red", lwd = 1.5)
      points((obs$a[i] + obs$b[i]) / 2, obs$E_obs[i],
             pch = 16, col = "red", cex = 1.5)
    } else {
      points(obs$a[i], obs$E_obs[i], pch = 17, col = "red", cex = 1.5)
    }
  }
  legend("topright", bty = "n",
         legend = c("Posterior mean", "95% CrI", "Observed"),
         col = c("darkgreen", rgb(0, 0.5, 0, 0.3), "red"),
         lty = c(1, 1, NA), lwd = c(2, 8, NA), pch = c(NA, NA, 16))
  
  # (b, c) Trace plots for free parameters
  if (n_free >= 1) {
    plot(chain[, 1], type = "l", col = rgb(0, 0, 0, 0.3),
         xlab = "Sample", ylab = free_names[1],
         main = sprintf("(b) Trace — %s", free_names[1]))
    abline(h = post_mean[free_names[1]], col = "blue", lwd = 1.5)
  }
  if (n_free >= 2) {
    plot(chain[, 2], type = "l", col = rgb(0, 0, 0, 0.3),
         xlab = "Sample", ylab = free_names[2],
         main = sprintf("(c) Trace — %s", free_names[2]))
    abline(h = post_mean[free_names[2]], col = "blue", lwd = 1.5)
  } else {
    plot.new()
  }
  
  # (d, e) Posterior marginal histograms
  hist_colors <- c(Emax = "lightblue", kappa = "lightgreen", tau = "lightyellow")
  all_names <- c("Emax", "kappa", "tau")
  
  if (n_free >= 1) {
    p1 <- free_idx[1]
    ci1 <- quantile(chain_full[, p1], c(0.025, 0.975))
    hist(chain_full[, p1], breaks = 50, col = hist_colors[p1], border = "white",
         main = sprintf("(d) Posterior — %s", all_names[p1]),
         xlab = all_names[p1], prob = TRUE)
    abline(v = post_mean[p1], col = "red", lwd = 2)
    abline(v = ci1, col = "red", lty = 2)
  }
  if (n_free >= 2) {
    p2 <- free_idx[2]
    ci2 <- quantile(chain_full[, p2], c(0.025, 0.975))
    hist(chain_full[, p2], breaks = 50, col = hist_colors[p2], border = "white",
         main = sprintf("(e) Posterior — %s", all_names[p2]),
         xlab = all_names[p2], prob = TRUE)
    abline(v = post_mean[p2], col = "red", lwd = 2)
    abline(v = ci2, col = "red", lty = 2)
  } else {
    plot.new()
  }
  
  # (f) Joint posterior scatter of the first two free parameters
  if (n_free >= 2) {
    p1 <- free_idx[1]; p2 <- free_idx[2]
    plot(chain_full[, p1], chain_full[, p2],
         pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.1),
         xlab = all_names[p1], ylab = all_names[p2],
         main = "(f) Joint posterior")
    points(post_mean[p1], post_mean[p2], pch = 4, col = "red", cex = 2, lwd = 2)
  } else {
    plot.new()
  }
  
  dev.off()
}

cat("=============================================================\n")
cat(" Finished! Check the 'figures' folder.\n")
cat("=============================================================\n")
###############################################################################
# STAGE 1 - SCENARIO B
###############################################################################

source("~/modular-bayesian-calibration/Stage1/stage1_core_full_fixed_m.R")


SHOW_PLOTS <- TRUE
SAVE_PLOTS <- TRUE
N_CHAINS <- 4
N_ITER <- 8000
BURN <- 2000
PROPOSAL_SD_K <- 0.8
PROPOSAL_SD_GAMMA <- 0.9

INIT_K <- c(0.03, 0.15, 0.45, 0.8)
INIT_GAMMA <- c(0.15, 0.3, 0.6, 0.9)
CHAIN_SEEDS <- 1201:1204

RESULT_DIR <- file.path("results", "stage1_scenario_B_confounding")
FIGURE_DIR <- file.path(RESULT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 1. DATA GENERATOR
###############################################################################

generate_confounding_data <- function(n = 45, t_max = 2.77, theta_true = c(g = 9.8, h0 = 46.45), lambda2_true = 0.01, m_true = 15, delta_g_conf = 1, seed = 321) {
  set.seed(seed)
  
  t_raw <- seq(0, t_max, length.out = n)
  prep <- prepare_time(t_raw)
  physical_true <- ball_drop_mean(t_raw, theta_true)
  
  delta_true <- rep(0, n)
  idx_D <- seq.int(m_true + 1, n)
  t_m <- t_raw[m_true]
  
  # Deterministic discrepancy aligned with the physical g direction.
  delta_true[idx_D] <- 0.5 * delta_g_conf * (t_raw[idx_D]^2 - t_m^2)
  epsilon_true <- rnorm(n, sd = sqrt(lambda2_true))

  y <- physical_true + delta_true + epsilon_true
  
  data <- data.frame(
    i = seq_len(n),
    t_raw = t_raw,
    t_norm = prep$u,
    y = y,
    physical_true = physical_true,
    delta_true = delta_true,
    epsilon_true = epsilon_true,
    segment = ifelse(seq_len(n) <= m_true, "physical", "discrepancy")
  )
  
  list(
    data = data,
    y = y,
    t_raw = t_raw,
    truth = list(
      theta = theta_true,
      lambda2 = lambda2_true,
      m = m_true,
      delta_g_conf = delta_g_conf,
      k = NULL,
      gamma_delta = NULL
    )
  )
}

confounding_projection <- function(dat) {
  X <- make_X(dat$t_raw)
  idx_D <- seq.int(dat$truth$m + 1, nrow(dat$data))
  X_D <- X[idx_D, , drop = FALSE]
  delta_D <- dat$data$delta_true[idx_D]
  
  beta <- as.vector(safe_solve(crossprod(X_D)) %*% crossprod(X_D, delta_D))
  names(beta) <- c("g", "h0")
  
  data.frame(
    quantity = c("projection_on_g", "projection_on_h0", "cor_delta_with_g_column"),
    value = c(beta["g"], beta["h0"], cor(delta_D, X_D[, "g"]))
  )
}

dat <- generate_confounding_data()
print(dat$truth)
projection <- confounding_projection(dat)
print(projection)

write.csv(dat$data, file.path(RESULT_DIR, "simulated_data.csv"), row.names = FALSE)
write.csv(projection, file.path(RESULT_DIR, "confounding_projection.csv"), row.names = FALSE)

if (SHOW_PLOTS) {
  plot_simulated_data(dat, "Stage 1 Scenario B: confounding data")
  plot_true_discrepancy(dat, "Stage 1 Scenario B: DGP true discrepancy")
}

###############################################################################
## 2. RUN FOUR INDEPENDENT CHAINS
###############################################################################

chains <- vector("list", N_CHAINS)
for (j in seq_len(N_CHAINS)) {
  chains[[j]] <- run_full_fixed_chain(
    y = dat$y,
    t_raw = dat$t_raw,
    m = dat$truth$m,
    n_iter = N_ITER,
    burn = BURN,
    init_k = INIT_K[j],
    init_gamma = INIT_GAMMA[j],
    proposal_sd_k = PROPOSAL_SD_K,
    proposal_sd_gamma = PROPOSAL_SD_GAMMA,
    seed = CHAIN_SEEDS[j]
  )
}

pooled_draws <- combine_chains(chains)
summary_table <- parameter_summary(pooled_draws)
diagnostics <- continuous_diagnostics(chains)
acceptance <- acceptance_summary(chains)

print(round(summary_table, 5))
print(diagnostics)
print(acceptance)

write.csv(pooled_draws, file.path(RESULT_DIR, "posterior_draws.csv"), row.names = FALSE)
write.csv(summary_table, file.path(RESULT_DIR, "posterior_summary.csv"))
write.csv(diagnostics, file.path(RESULT_DIR, "diagnostics.csv"), row.names = FALSE)
write.csv(acceptance, file.path(RESULT_DIR, "acceptance.csv"), row.names = FALSE)

###############################################################################
## 4. PHYSICAL-PARAMETER RECOVERY
###############################################################################

physical_parameters <- c("g", "h0", "lambda2")
truth <- c(
  g = unname(dat$truth$theta["g"]),
  h0 = unname(dat$truth$theta["h0"]),
  lambda2 = dat$truth$lambda2
)

recovery <- data.frame(
  parameter = physical_parameters,
  truth = as.numeric(truth[physical_parameters]),
  mean = summary_table[physical_parameters, "mean"],
  q025 = summary_table[physical_parameters, "q025"],
  q975 = summary_table[physical_parameters, "q975"],
  row.names = NULL
)
recovery$error <- recovery$mean - recovery$truth
recovery$contains_truth <- recovery$truth >= recovery$q025 & recovery$truth <= recovery$q975
print(recovery)
write.csv(recovery, file.path(RESULT_DIR, "physical_parameter_recovery.csv"), row.names = FALSE)

cat(
  "\nNOTE: Scenario B has no data-generating truth for k or gamma_delta.\n",
  "Those parameters describe the fitted GP discrepancy model, not the deterministic DGP discrepancy.\n"
)

###############################################################################
## 5. PLOTS
###############################################################################

plot_parameter <- function(x, truth = NULL, xlab, main) {
  plot(density(x), lwd = 2, xlab = xlab, main = main)
  if (!is.null(truth)) abline(v = truth, lty = 2, lwd = 2)
}

if (SHOW_PLOTS) {
  par(mfrow = c(2, 3))
  plot_parameter(pooled_draws$g, dat$truth$theta["g"], "g", "Posterior of g")
  plot_parameter(pooled_draws$h0, dat$truth$theta["h0"], "h0", "Posterior of h0")
  plot_parameter(pooled_draws$lambda2, dat$truth$lambda2, expression(lambda^2), "Posterior of lambda2")
  plot_parameter(pooled_draws$k, NULL, "k", "Posterior of k")
  plot_parameter(pooled_draws$gamma_delta, NULL, expression(gamma[delta]), "Posterior of gamma")
  par(mfrow = c(1, 1))
}

if (SAVE_PLOTS) {
  png(file.path(FIGURE_DIR, "posterior_parameters.png"), width = 1400, height = 900, res = 140)
  par(mfrow = c(2, 3))
  plot_parameter(pooled_draws$g, dat$truth$theta["g"], "g", "Posterior of g")
  plot_parameter(pooled_draws$h0, dat$truth$theta["h0"], "h0", "Posterior of h0")
  plot_parameter(pooled_draws$lambda2, dat$truth$lambda2, expression(lambda^2), "Posterior of lambda2")
  plot_parameter(pooled_draws$k, NULL, "k", "Posterior of k")
  plot_parameter(pooled_draws$gamma_delta, NULL, expression(gamma[delta]), "Posterior of gamma")
  par(mfrow = c(1, 1))
  dev.off()
}

cat("\nStage 1 Scenario B complete. Results written to:\n", RESULT_DIR, "\n")

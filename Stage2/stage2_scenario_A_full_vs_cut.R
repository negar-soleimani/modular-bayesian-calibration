###############################################################################
# STAGE 2 - SCENARIO A
###############################################################################

source("~/modular-bayesian-calibration/Stage2/stage2_core_cut_fixed_m.R")

SHOW_PLOTS <- TRUE
SAVE_PLOTS <- TRUE
N_CHAINS <- 4

# Outer Module-P
OUTER_ITER <- 1400
OUTER_BURN <- 400
OUTER_THIN <- 5

# Module-D chain for every retained Module-P.
INNER_ITER <- 600
INNER_BURN <- 200
PROPOSAL_SD_K <- 0.02
PROPOSAL_SD_GAMMA <- 0.15

INIT_LAMBDA2 <- c(0.004, 0.01, 0.03, 0.08)
INIT_K <- c(0.03, 0.15, 0.45, 0.8)
INIT_GAMMA <- c(0.15, 0.3, 0.6, 0.9)
CHAIN_SEEDS <- 2101:2104

RESULT_DIR <- file.path("results", "stage2_scenario_A_well_specified_gp_cut")
FIGURE_DIR <- file.path(RESULT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 1. DATA GENERATOR
###############################################################################

generate_well_specified_gp_data <- function(n = 45, t_max = 2.77, theta_true = c(g = 9.8, h0 = 46.45),
                                   lambda2_true = 0.01, m_true = 15, k_true = 0.1, gamma_true = 0.3, seed = 123) {
  set.seed(seed)
  
  t_raw <- seq(0, t_max, length.out = n)
  prep <- prepare_time(t_raw)
  physical_true <- ball_drop_mean(t_raw, theta_true)
  
  delta_true <- rep(0, n)
  idx_D <- seq.int(m_true + 1, n)
  Sigma_delta <- GP_covariance(prep$u[idx_D], sigma_sq_delta = lambda2_true / k_true, gamma_delta = gamma_true)
  delta_true[idx_D] <- as.vector(rmvnorm(1, mean = rep(0, length(idx_D)), sigma = Sigma_delta))
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
      k = k_true,
      gamma_delta = gamma_true
    )
  )
}

###############################################################################
## 2. SIMULATE DATA
###############################################################################

dat <- generate_well_specified_gp_data()
print(dat$truth)
write.csv(dat$data, file.path(RESULT_DIR, "simulated_data.csv"), row.names = FALSE)

if (SHOW_PLOTS) {
  plot_simulated_data(dat, "Stage 2 Scenario A: well-specified GP data")
  plot_true_discrepancy(dat, "Stage 2 Scenario A: DGP true discrepancy")
}

###############################################################################
## 3. CUT CHAINS
###############################################################################

chains <- vector("list", N_CHAINS)
for (j in seq_len(N_CHAINS)) {
  chains[[j]] <- run_cut_fixed_chain(
    y = dat$y,
    t_raw = dat$t_raw,
    m = dat$truth$m,
    outer_iter = OUTER_ITER,
    outer_burn = OUTER_BURN,
    outer_thin = OUTER_THIN,
    inner_iter = INNER_ITER,
    inner_burn = INNER_BURN,
    init_lambda2 = INIT_LAMBDA2[j],
    init_k = INIT_K[j],
    init_gamma = INIT_GAMMA[j],
    proposal_sd_k = PROPOSAL_SD_K,
    proposal_sd_gamma = PROPOSAL_SD_GAMMA,
    seed = CHAIN_SEEDS[j],
    downstream_pick = "random",
    verbose = TRUE
  )
}

pooled_draws <- combine_cut_chains(chains)
summary_table <- parameter_summary(pooled_draws)
upstream_diagnostics <- module_P_diagnostics(chains)
inner_acceptance <- inner_acceptance_summary(chains)

print(summary_table)
print(upstream_diagnostics)
print(inner_acceptance)

write.csv(pooled_draws, file.path(RESULT_DIR, "cut_posterior_draws.csv"), row.names = FALSE)
write.csv(summary_table, file.path(RESULT_DIR, "cut_posterior_summary.csv"))
write.csv(upstream_diagnostics, file.path(RESULT_DIR, "module_P_diagnostics.csv"), row.names = FALSE)
write.csv(inner_acceptance, file.path(RESULT_DIR, "module_D_acceptance.csv"), row.names = FALSE)

###############################################################################

truth <- c(
  g = unname(dat$truth$theta["g"]),
  h0 = unname(dat$truth$theta["h0"]),
  lambda2 = dat$truth$lambda2,
  k = dat$truth$k,
  gamma_delta = dat$truth$gamma_delta
)

recovery <- data.frame(
  parameter = rownames(summary_table),
  truth = as.numeric(truth[rownames(summary_table)]),
  mean = summary_table$mean,
  q025 = summary_table$q025,
  q975 = summary_table$q975,
  row.names = NULL
)
recovery$error <- recovery$mean - recovery$truth
recovery$contains_truth <- recovery$truth >= recovery$q025 & recovery$truth <= recovery$q975
print(recovery)
write.csv(recovery, file.path(RESULT_DIR, "parameter_recovery.csv"), row.names = FALSE)

###############################################################################
## 4. PLOTS
###############################################################################

plot_parameter <- function(x, truth, xlab, main) {
  plot(density(x), lwd = 2, xlab = xlab, main = main)
  abline(v = truth, lty = 2, lwd = 2)
}

if (SHOW_PLOTS) {
  par(mfrow = c(2, 3))
  plot_parameter(pooled_draws$g, dat$truth$theta["g"], "g", "Cut posterior of g")
  plot_parameter(pooled_draws$h0, dat$truth$theta["h0"], "h0", "Cut posterior of h0")
  plot_parameter(pooled_draws$lambda2, dat$truth$lambda2, expression(lambda^2), "Cut posterior of lambda2")
  plot_parameter(pooled_draws$k, dat$truth$k, "k", "Cut posterior of k")
  plot_parameter(pooled_draws$gamma_delta, dat$truth$gamma_delta, expression(gamma[delta]), "Cut posterior of gamma")
  par(mfrow = c(1, 1))
}

if (SAVE_PLOTS) {
  png(file.path(FIGURE_DIR, "posterior_parameters.png"), width = 1400, height = 900, res = 140)
  par(mfrow = c(2, 3))
  plot_parameter(pooled_draws$g, dat$truth$theta["g"], "g", "Cut posterior of g")
  plot_parameter(pooled_draws$h0, dat$truth$theta["h0"], "h0", "Cut posterior of h0")
  plot_parameter(pooled_draws$lambda2, dat$truth$lambda2, expression(lambda^2), "Cut posterior of lambda2")
  plot_parameter(pooled_draws$k, dat$truth$k, "k", "Cut posterior of k")
  plot_parameter(pooled_draws$gamma_delta, dat$truth$gamma_delta, expression(gamma[delta]), "Cut posterior of gamma")
  par(mfrow = c(1, 1))
  dev.off()
}

cat("\nStage 2 Scenario A complete. Results written to:\n", RESULT_DIR, "\n")


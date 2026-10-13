###############################################################################
# ARTICLE 2 - STAGE 1 - SCENARIO A
# WELL-SPECIFIED GP DISCREPANCY WITH FIXED, KNOWN m
#
# Purpose:
#   Validate the Stage-1 Full fixed-m implementation in a benign setting where
#   the post-boundary discrepancy is generated from the SAME GP family used
#   for inference.
# This is a control / code-validation scenario. It is not the confounding
# stress test.
###############################################################################

source("stage1_core_full_fixed_m.R")

SHOW_PLOTS <- TRUE

###############################################################################
## A1. DATA GENERATOR
###############################################################################

generate_well_specified_gp_data <- function(
    n = 45,
    t_max = 2.77,
    theta_true = c(g = 9.8, h0 = 46.45),
    lambda2_true = 0.01,
    m_true = 15,
    k_true = 0.1,
    gamma_true = 0.3,
    seed = 123
) {
  set.seed(seed)
  
  t_raw <- seq(0, t_max, length.out = n)
  prep <- prepare_time(t_raw)
  t <- prep$t
  
  physical_true <- balldropg(t, theta = theta_true, t_min = prep$t_min, t_range = prep$t_range)
  
  delta_true <- rep(0, n)
  idx_D <- (m_true + 1):n
  
  # The GP is defined on normalized time, consistently with make_C().
  Sigma_delta <- GP_covariance(t[idx_D],
                               sigma_sq_delta = lambda2_true / k_true,
                               gamma_delta = gamma_true
  )
  
  delta_true[idx_D] <- as.vector(rmvnorm(1, mean = rep(0, length(idx_D)), sigma = Sigma_delta))
  
  epsilon <- rnorm(n, mean = 0, sd = sqrt(lambda2_true))
  
  y <- physical_true + delta_true + epsilon
  
  data <- data.frame(
    i = 1:n,
    t_raw = t_raw,
    t_norm = t,
    y = y,
    physical_true = physical_true,
    delta_true = delta_true,
    epsilon_true = epsilon,
    segment = ifelse(
      1:n <= m_true,
      "physical",
      "discrepancy"
    )
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
## A2. GENERATE ONE DATASET
###############################################################################

dat <- generate_well_specified_gp_data(seed = 123, m_true = 15)

print(dat$truth)
print(dat$data[12:19, ])

if (any(dat$data$physical_true <= 0)) {
  warning("Physical trajectory reaches zero; inspect balldropg/make_X consistency.")
}


###############################################################################
## A2b. SIMULATED DATA AND TRUE DISCREPANCY IN R
###############################################################################

if (SHOW_PLOTS) {
  plot_stage1_data(dat, main = "Stage 1 - Scenario A: well-specified GP data")
  plot_stage1_true_delta(dat, main = "Stage 1 - Scenario A: true GP discrepancy")
}

###############################################################################
## A3. AUDIT THE FIXED-m COLLAPSED CALCULATION AT THE TRUE GP HYPERPARAMETERS
###############################################################################

prep <- prepare_time(dat$t_raw)

st_true <- full_fixed_stats(
  y = dat$y,
  t = prep$t,
  t_min = prep$t_min,
  t_range = prep$t_range,
  m = dat$truth$m,
  k = dat$truth$k,
  gamma_delta = dat$truth$gamma_delta
)

print(st_true$theta_hat)
print(st_true$SSE)
print(st_true$lambda_shape)
print(st_true$lambda_rate)
print(st_true$log_marginal)

# With k_true = 0.1 and R_ii = 1:
#   C_ii = 1 + 1/k = 11.
print(diag(st_true$C)[1:5])

###############################################################################
## A4. SHORT SMOKE TEST
###############################################################################

test_chain <- run_full_fixed_chain(
  y = dat$y,
  t_raw = dat$t_raw,
  m = dat$truth$m,
  n_iter = 300,
  burn = 100,
  init_k = 0.1,
  init_gamma = 0.3,
  proposal_sd_k = 0.2,
  proposal_sd_gamma = 0.2,
  seed = 101
)

print(head(test_chain$chain))
print(test_chain$acceptance)

###############################################################################
## A5. FOUR FULL FIXED-m CHAINS
###############################################################################

init_k <- c(0.03, 0.15, 0.45, 0.8)
init_gamma <- c(0.15, 0.3, 0.6, 0.9)

chains <- lapply(
  1:4,
  function(j) {
    cat("Running Stage-1 Scenario-A chain", j, "\n")
    
    run_full_fixed_chain(
      y = dat$y,
      t_raw = dat$t_raw,
      m = dat$truth$m,
      n_iter = 8000,
      burn = 2000,
      init_k = init_k[j],
      init_gamma = init_gamma[j],
      proposal_sd_k = 0.20,
      proposal_sd_gamma = 0.20,
      seed = 100 + j
    )
  }
)

###############################################################################
## A6. DIAGNOSTICS AND PARAMETER SUMMARY
###############################################################################

diagnostics <- summarize_full_fixed_chains(chains)
print(diagnostics$gelman)
print(diagnostics$ess)
print(diagnostics$acceptance)

draws <- pool_full_fixed_draws(chains)
parameter_summary <- make_parameter_summary(draws)
print(round(parameter_summary, 5))

###############################################################################
## A7. POSTERIOR RECONSTRUCTION OF delta
###############################################################################

delta_summary <- reconstruct_delta_fixed(
  draws = draws,
  y = dat$y,
  t_raw = dat$t_raw,
  m = dat$truth$m,
  n_reconstruct = 300,
  seed = 2026
)


###############################################################################
## A7b. SHOW POSTERIOR DIAGNOSTIC FIGURES IN R
###############################################################################

if (SHOW_PLOTS) {
  plot_stage1_trace(chains, main_prefix = "Scenario A")
  plot_stage1_posterior(draws, dat$truth, "Scenario A", show_gp_truth = TRUE)
  plot_stage1_delta_reconstruction(
    dat, delta_summary,
    main = "Scenario A: posterior discrepancy reconstruction"
  )
}

###############################################################################
## A8. SAVE OUTPUTS 
###############################################################################

result_dir <- file.path(
  "results",
  "stage1_scenario_A_well_specified_gp"
)

save_stage1_outputs(
  result_dir = result_dir,
  truth = dat$truth,
  chains = chains,
  draws = draws,
  parameter_summary = parameter_summary,
  diagnostics = diagnostics,
  delta_summary = delta_summary
)

write.csv(
  dat$data,
  file.path(result_dir, "simulated_data.csv"),
  row.names = FALSE
)


figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_png(
  function() plot_stage1_data(dat, main = "Stage 1 - Scenario A: well-specified GP data"),
  file.path(figure_dir, "01_simulated_data.png")
)
save_plot_png(
  function() plot_stage1_true_delta(dat, main = "Stage 1 - Scenario A: true GP discrepancy"),
  file.path(figure_dir, "02_true_discrepancy.png")
)
save_plot_png(
  function() plot_stage1_trace(chains, main_prefix = "Scenario A"),
  file.path(figure_dir, "03_trace_plots.png")
)
save_plot_png(
  function() plot_stage1_posterior(draws, dat$truth, "Scenario A", show_gp_truth = TRUE),
  file.path(figure_dir, "04_posterior_distributions.png")
)
save_plot_png(
  function() plot_stage1_delta_reconstruction(
    dat, delta_summary,
    main = "Scenario A: posterior discrepancy reconstruction"
  ),
  file.path(figure_dir, "05_delta_reconstruction.png")
)

save_stage1_text_report(result_dir, dat$truth, parameter_summary, diagnostics)

cat("\nScenario A outputs written to:\n", normalizePath(result_dir), "\n")
cat("Figures saved in:\n", normalizePath(figure_dir), "\n")
cat("Use the RStudio Plots-pane arrows to browse plots shown during the run.\n")


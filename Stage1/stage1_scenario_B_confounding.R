###############################################################################
# ARTICLE 2 - STAGE 1 - SCENARIO B
# FIXED-m FULL MODEL UNDER A CONFOUNDING 
###############################################################################

source("stage1_core_full_fixed_m.R")

SHOW_PLOTS <- TRUE

###############################################################################
## B1. CONFOUNDING DATA GENERATOR
###############################################################################

generate_confounding_data <- function(
    n = 45,
    t_max = 2.77,
    theta_true = c(g = 9.8, h0 = 46.45),
    lambda2_true = 0.01,
    m_true = 15,
    delta_g_conf = 1.0,
    seed = 321
) {
  set.seed(seed)
  
  t_raw <- seq(0, t_max, length.out = n)
  prep <- prepare_time(t_raw)
  t <- prep$t
  
  physical_true <- balldropg(t, theta = theta_true, t_min = prep$t_min, t_range = prep$t_range)
  
  delta_true <- rep(0, n)
  idx_D <- (m_true + 1):n
  
  # Last clean physical time.
  t_boundary_ref <- t_raw[m_true]
  
  # Deterministic discrepancy aligned with the g direction:
  #   delta_true(t_i)
  #     = 0.5 * delta_g_conf * (t_i^2 - t_m^2),  i > m_true.
  #
  # It is zero at the last clean point and increases afterwards.
  delta_true[idx_D] <- 0.5 * delta_g_conf * (t_raw[idx_D]^2 - t_boundary_ref^2)
  
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
      delta_g_conf = delta_g_conf
    )
  )
}

confounding_projection <- function(dat) {
  prep <- prepare_time(dat$t_raw)
  X <- make_X(prep$t, prep$t_min, prep$t_range)
  
  m <- dat$truth$m
  idx_D <- (m + 1):nrow(dat$data)
  
  X_D <- X[idx_D, , drop = FALSE]
  delta_D <- dat$data$delta_true[idx_D]
  
  beta_projection <- as.vector(safe_solve(t(X_D) %*% X_D) %*% t(X_D) %*% delta_D)
  names(beta_projection) <- c("g", "h0")
  
  data.frame(
    quantity = c(
      "projection_on_g",
      "projection_on_h0",
      "cor_delta_with_g_column"
    ),
    value = c(
      beta_projection["g"],
      beta_projection["h0"],
      cor(delta_D, X_D[, "g"])
    )
  )
}

###############################################################################
## B2. GENERATE ONE CONFOUNDING DATASET
###############################################################################

dat <- generate_confounding_data(seed = 321, m_true = 15, delta_g_conf = 1.0)

print(dat$truth)
print(dat$data[12:20, ])

if (any(dat$data$physical_true <= 0)) {
  warning("Physical trajectory reaches zero; inspect balldropg/make_X consistency.")
}

projection_result <- confounding_projection(dat)
print(projection_result)


###############################################################################
## B2b. SIMULATED DATA AND TRUE DISCREPANCY
###############################################################################

if (SHOW_PLOTS) {
  plot_stage1_data(dat, main = "Stage 1 - Scenario B: confounding data")
  plot_stage1_true_delta(dat, main = "Stage 1 - Scenario B: true confounding discrepancy")
}

###############################################################################
## B3. AUDIT THE FIXED-m FULL CALCULATION AT REFERENCE GP VALUES
###############################################################################

reference_k <- 0.10
reference_gamma <- 0.30

prep <- prepare_time(dat$t_raw)

st_reference <- full_fixed_stats(
  y = dat$y,
  t = prep$t,
  t_min = prep$t_min,
  t_range = prep$t_range,
  m = dat$truth$m,
  k = reference_k,
  gamma_delta = reference_gamma
)

print(st_reference$theta_hat)
print(st_reference$SSE)
print(st_reference$lambda_shape)
print(st_reference$lambda_rate)

###############################################################################
## B4. SHORT SMOKE TEST
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
## B5. FOUR FULL FIXED-m CHAINS
###############################################################################

init_k <- c(0.03, 0.15, 0.45, 0.8)
init_gamma <- c(0.15, 0.3, 0.6, 0.9)

chains <- lapply(
  1:4,
  function(j) {
    cat("Running Stage-1 Scenario-B chain", j, "\n")
    
    run_full_fixed_chain(
      y = dat$y,
      t_raw = dat$t_raw,
      m = dat$truth$m,
      n_iter = 8000,
      burn = 2000,
      init_k = init_k[j],
      init_gamma = init_gamma[j],
      proposal_sd_k = 0.2,
      proposal_sd_gamma = 0.2,
      seed = 100 + j
    )
  }
)

###############################################################################
## B6. DIAGNOSTICS AND PARAMETER SUMMARY
###############################################################################

diagnostics <- summarize_full_fixed_chains(chains)
print(diagnostics$gelman)
print(diagnostics$ess)
print(diagnostics$acceptance)

draws <- pool_full_fixed_draws(chains)
parameter_summary <- make_parameter_summary(draws)
print(round(parameter_summary, 5))

###############################################################################
## B7. POSTERIOR RECONSTRUCTION OF delta
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
## B7b. POSTERIOR DIAGNOSTIC FIGURES 
###############################################################################

if (SHOW_PLOTS) {
  plot_stage1_trace(chains, main_prefix = "Scenario B")
  plot_stage1_posterior(draws, dat$truth, "Scenario B", show_gp_truth = FALSE)
  plot_stage1_delta_reconstruction(
    dat, delta_summary,
    main = "Scenario B: posterior discrepancy reconstruction"
  )
}

###############################################################################
## B8. SAVE OUTPUTS
###############################################################################

result_dir <- file.path(
  "results",
  "stage1_scenario_B_confounding"
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

write.csv(
  projection_result,
  file.path(result_dir, "confounding_projection.csv"),
  row.names = FALSE
)


figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_png(
  function() plot_stage1_data(dat, main = "Stage 1 - Scenario B: confounding data"),
  file.path(figure_dir, "01_simulated_data.png")
)
save_plot_png(
  function() plot_stage1_true_delta(dat, main = "Stage 1 - Scenario B: true confounding discrepancy"),
  file.path(figure_dir, "02_true_discrepancy.png")
)
save_plot_png(
  function() plot_stage1_trace(chains, main_prefix = "Scenario B"),
  file.path(figure_dir, "03_trace_plots.png")
)
save_plot_png(
  function() plot_stage1_posterior(draws, dat$truth, "Scenario B", show_gp_truth = FALSE),
  file.path(figure_dir, "04_posterior_distributions.png")
)
save_plot_png(
  function() plot_stage1_delta_reconstruction(
    dat, delta_summary,
    main = "Scenario B: posterior discrepancy reconstruction"
  ),
  file.path(figure_dir, "05_delta_reconstruction.png")
)

save_stage1_text_report(result_dir, dat$truth, parameter_summary, diagnostics)

cat("\nScenario B outputs written to:\n", normalizePath(result_dir), "\n")
cat("Figures saved in:\n", normalizePath(figure_dir), "\n")
cat("Use the RStudio Plots-pane arrows to browse plots shown during the run.\n")


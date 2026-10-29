###############################################################################
# STAGE 2 - SCENARIO B
# FIXED-m FULL vs CUT COMPARISON
###############################################################################

source("stage2_core_cut_fixed_m_corrected.R")

SHOW_PLOTS <- TRUE
SAVE_PLOTS <- TRUE

# Main Stage-2 computational settings.
N_CUT_PHYSICAL <- 10000

# The physical Cut block above is exact direct sampling.
# The following settings are only for the downstream GP block.
N_CUT_OUTER <- 300
INNER_ITER <- 1200
INNER_BURN <- 400
PROPOSAL_SD_K <- 0.02
PROPOSAL_SD_GAMMA <- 0.15

RESULT_DIR <- file.path("results", "stage2_scenario_B_full_vs_cut")

FIGURE_DIR <- file.path(RESULT_DIR, "figures")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\nStage-2 Scenario B output directory:\n")
cat(normalizePath(RESULT_DIR, mustWork = FALSE), "\n")

###############################################################################
## 1. DATASET AND FULL POSTERIOR DRAWS
###############################################################################

stage1_result_dir <- file.path("results", "stage1_scenario_B_confounding")

data_file <- file.path(stage1_result_dir, "simulated_data.csv")

full_draw_file <- file.path(stage1_result_dir, "pooled_draws.rds")

if (!file.exists(data_file) || !file.exists(full_draw_file)) {
  stop(
    paste0(
      "Stage-1 saved outputs were not found. ",
      "Run the corresponding Stage-1 scenario first.\n",
      "Expected files:\n",
      data_file, "\n",
      full_draw_file
    )
  )
}

data_stage1 <- read.csv(data_file, stringsAsFactors = FALSE)

full_fixed_draws <- readRDS(full_draw_file)

truth <- list(
  theta = c(g = 9.8, h0 = 46.45),
  lambda2 = 0.01,
  m = 15,
  delta_g_conf = 1
)

dat <- list(
  data = data_stage1,
  y = data_stage1$y,
  t_raw = data_stage1$t_raw,
  truth = truth
)

# Basic checks before the Cut analysis.
stopifnot(
  length(dat$y) == nrow(dat$data),
  truth$m >= 3,
  truth$m < length(dat$y)
)

cat("\nLoaded Stage-1 Scenario B data and Full draws.\n")
cat("Number of Full retained draws:", nrow(full_fixed_draws), "\n")
print(truth)

###############################################################################
## 2. AUDIT MODULE P AT THE FIXED m
###############################################################################

prep <- prepare_time(dat$t_raw)

module_P_check <- module_P_stats(y = dat$y, t = prep$t, t_min = prep$t_min,
               t_range = prep$t_range, m = truth$m, prior_type = PRIOR_TYPE)

cat("\nModule-P center theta_hat_P:\n")
print(module_P_check$theta_hat)

cat("\nModule-P S_P:\n")
print(module_P_check$SSE)

cat("\nModule-P lambda^2 posterior parameters:\n")
print(
  c(
    shape = module_P_check$lambda_shape,
    rate = module_P_check$lambda_rate
  )
)

###############################################################################
## 3. PRIMARY CUT RESULT: DIRECT DRAWS FOR theta AND lambda^2
###############################################################################

cut_physical_draws <- run_cut_physical_fixed(y = dat$y,
  t_raw = dat$t_raw, m = truth$m,
  n_draws = N_CUT_PHYSICAL,
  prior_type = PRIOR_TYPE,
  seed = 2999
)

cut_physical_summary <- rbind(
  g = summarize_parameter(cut_physical_draws$g),
  h0 = summarize_parameter(cut_physical_draws$h0),
  lambda2 = summarize_parameter(cut_physical_draws$lambda2)
)

cat("\nCut physical-module posterior summary:\n")
print(round(cut_physical_summary, 5))

###############################################################################
## 4. FULL vs CUT COMPARISON FOR THE PHYSICAL PARAMETERS
###############################################################################

physical_comparison <- make_full_cut_physical_comparison(
  full_draws = full_fixed_draws,
  cut_physical_draws = cut_physical_draws,
  truth = truth
)

physical_comparison_print <- physical_comparison
num_cols <- sapply(physical_comparison_print, is.numeric)
physical_comparison_print[, num_cols] <-
  round(
    physical_comparison_print[, num_cols, drop = FALSE],
    5
  )

cat("\nFull vs Cut physical comparison:\n")
print(physical_comparison_print)

if (SHOW_PLOTS) {
  plot_full_cut_physical(
    full_draws = full_fixed_draws,
    cut_draws = cut_physical_draws,
    truth = truth,
    scenario_label = "Scenario B"
  )
}

saveRDS(cut_physical_draws, file.path(RESULT_DIR, "cut_physical_draws.rds"))
write.csv(cut_physical_summary, file.path(RESULT_DIR, "cut_physical_summary.csv"), row.names = TRUE)
write.csv(physical_comparison, file.path(RESULT_DIR, "full_cut_physical_comparison.csv"), row.names = FALSE)

if (SAVE_PLOTS) {
  save_plot_png(
    function() plot_full_cut_physical(
      full_draws = full_fixed_draws,
      cut_draws = cut_physical_draws,
      truth = truth,
      scenario_label = "Scenario B"
    ),
    file.path(
      FIGURE_DIR,
      "01_full_vs_cut_physical.png"
    ),
    width = 1800,
    height = 700
  )
}

cat("\nPrimary physical Full-vs-Cut results saved in:\n")
cat(normalizePath(RESULT_DIR, mustWork = FALSE), "\n")

###############################################################################
## 5. SHORT DOWNSTREAM CUT SMOKE TEST
###############################################################################

cut_smoke <- run_cut_fixed_m(
  y = dat$y,
  t_raw = dat$t_raw,
  m = truth$m,
  n_outer = 20,
  inner_iter = 500,
  inner_burn = 150,
  proposal_sd_k = PROPOSAL_SD_K,
  proposal_sd_gamma = PROPOSAL_SD_GAMMA,
  prior_type = PRIOR_TYPE,
  seed = 3001,
  verbose = TRUE
)

cat("\nDownstream smoke-test acceptance means:\n")
print(
  colMeans(
    cut_smoke[, c("acc_k", "acc_gamma"), drop = FALSE]
  )
)

###############################################################################
## 6. k AND gamma_delta
###############################################################################

cut_joint_draws <- run_cut_fixed_m(
  y = dat$y,
  t_raw = dat$t_raw,
  m = truth$m,
  n_outer = N_CUT_OUTER,
  inner_iter = INNER_ITER,
  inner_burn = INNER_BURN,
  proposal_sd_k = PROPOSAL_SD_K,
  proposal_sd_gamma = PROPOSAL_SD_GAMMA,
  prior_type = PRIOR_TYPE,
  seed = 3002,
  verbose = TRUE
)

cut_joint_summary <- make_parameter_summary(
  cut_joint_draws
)

cat("\nJoint Cut summary:\n")
print(round(cut_joint_summary, 5))

downstream_acceptance <- data.frame(
  parameter = c("k", "gamma_delta"),
  mean_acceptance = c(
    mean(cut_joint_draws$acc_k),
    mean(cut_joint_draws$acc_gamma)
  ),
  min_acceptance = c(
    min(cut_joint_draws$acc_k),
    min(cut_joint_draws$acc_gamma)
  ),
  max_acceptance = c(
    max(cut_joint_draws$acc_k),
    max(cut_joint_draws$acc_gamma)
  )
)

cat("\nCut downstream inner-MH acceptance summary:\n")
print(round_numeric_columns(downstream_acceptance, digits = 4))

###############################################################################
## 7. SECONDARY FULL vs CUT COMPARISON FOR k AND gamma_delta
###############################################################################

# No GP truth exists in Scenario B because the true discrepancy is deterministic.
truth_gp <- NULL

if (SHOW_PLOTS) {
  plot_full_cut_gp(
    full_draws = full_fixed_draws,
    cut_joint_draws = cut_joint_draws,
    scenario_label = "Scenario B",
    truth_gp = truth_gp
  )
}


###############################################################################
## 8. RECONSTRUCT THE DISCREPANCY UNDER FULL AND CUT
###############################################################################

delta_full <- reconstruct_delta_fixed(
  draws = full_fixed_draws,
  y = dat$y,
  t_raw = dat$t_raw,
  m = truth$m,
  n_reconstruct = 300,
  seed = 3111
)

delta_cut <- reconstruct_delta_fixed(
  draws = cut_joint_draws,
  y = dat$y,
  t_raw = dat$t_raw,
  m = truth$m,
  n_reconstruct = min(300, nrow(cut_joint_draws)),
  seed = 3222
)

idx_D <- (truth$m + 1):nrow(dat$data)
full_delta_ISE <- mean((delta_full$mean[idx_D] - dat$data$delta_true[idx_D])^2)
cut_delta_ISE <- mean((delta_cut$mean[idx_D] - dat$data$delta_true[idx_D])^2)

delta_error <- data.frame(method = c("Full", "Cut"),
  downstream_mean_squared_error = c(
    full_delta_ISE,
    cut_delta_ISE
  )
)

cat("\nDiscrepancy posterior-mean error for this one dataset:\n")
print(delta_error)

if (SHOW_PLOTS) {
  plot_full_cut_delta_means(
    dat = dat,
    delta_full = delta_full,
    delta_cut = delta_cut,
    scenario_label = "Scenario B"
  )
  
  plot_full_cut_delta_intervals(
    dat = dat,
    delta_full = delta_full,
    delta_cut = delta_cut,
    scenario_label = "Scenario B"
  )
}

###############################################################################
## 9. SAVE
###############################################################################

result_dir <- RESULT_DIR

save_stage2_outputs(
  result_dir = result_dir,
  full_draws = full_fixed_draws,
  cut_physical_draws = cut_physical_draws,
  cut_joint_draws = cut_joint_draws,
  physical_comparison = physical_comparison,
  cut_physical_summary = cut_physical_summary,
  cut_joint_summary = cut_joint_summary,
  downstream_acceptance = downstream_acceptance,
  delta_full = delta_full,
  delta_cut = delta_cut,
  delta_error = delta_error
)

settings <- data.frame(
  quantity = c(
    "fixed_m",
    "full_retained_draws_loaded_from_stage1",
    "cut_physical_direct_draws",
    "cut_outer_draws",
    "inner_iterations_per_outer_draw",
    "inner_burn_per_outer_draw",
    "proposal_sd_k",
    "proposal_sd_gamma"
  ),
  value = c(
    truth$m,
    nrow(full_fixed_draws),
    N_CUT_PHYSICAL,
    N_CUT_OUTER,
    INNER_ITER,
    INNER_BURN,
    PROPOSAL_SD_K,
    PROPOSAL_SD_GAMMA
  )
)

write.csv(
  settings,
  file.path(result_dir, "algorithm_settings.csv"),
  row.names = FALSE
)

###############################################################################

figure_dir <- FIGURE_DIR

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (SAVE_PLOTS) {
  save_plot_png(
    function() plot_full_cut_physical(
      full_draws = full_fixed_draws,
      cut_draws = cut_physical_draws,
      truth = truth,
      scenario_label = "Scenario B"
    ),
    file.path(
      figure_dir,
      "01_full_vs_cut_physical.png"
    ),
    width = 1800,
    height = 700
  )
  
  save_plot_png(
    function() plot_full_cut_gp(
      full_draws = full_fixed_draws,
      cut_joint_draws = cut_joint_draws,
      scenario_label = "Scenario B",
      truth_gp = truth_gp
    ),
    file.path(
      figure_dir,
      "02_full_vs_cut_gp.png"
    ),
    width = 1400,
    height = 700
  )
  
  save_plot_png(
    function() plot_full_cut_delta_means(
      dat = dat,
      delta_full = delta_full,
      delta_cut = delta_cut,
      scenario_label = "Scenario B"
    ),
    file.path(
      figure_dir,
      "03_full_vs_cut_delta_means.png"
    ),
    width = 1400,
    height = 900
  )
  
  save_plot_png(
    function() plot_full_cut_delta_intervals(
      dat = dat,
      delta_full = delta_full,
      delta_cut = delta_cut,
      scenario_label = "Scenario B"
    ),
    file.path(
      figure_dir,
      "04_full_vs_cut_delta_intervals.png"
    ),
    width = 1800,
    height = 800
  )
}


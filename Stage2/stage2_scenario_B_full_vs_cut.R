###############################################################################
# ARTICLE 2 - STAGE 2 - SCENARIO B
# FIXED-m FULL vs CUT COMPARISON
#
# IMPORTANT:
# - This script uses the SAME simulated dataset and SAME Full posterior draws
#   already produced in Stage 1.
# - It does NOT regenerate a different dataset.
# - It does NOT rerun the Full chain.
# - It adds the fixed-m Cut analysis and Full-vs-Cut comparison.
###############################################################################

source("stage2_core_cut_fixed_m_corrected.R")

SHOW_PLOTS <- TRUE
SAVE_PLOTS <- TRUE

# Main Stage-2 computational settings.
N_CUT_PHYSICAL <- 10000L

# The physical Cut block above is exact direct sampling.
# The following settings are only for the downstream GP block.
N_CUT_OUTER <- 300L
INNER_ITER <- 1200L
INNER_BURN <- 400L
PROPOSAL_SD_K <- 0.02
PROPOSAL_SD_GAMMA <- 0.15

# Explicit Stage-2 output directory for this scenario.
# It is created before the downstream calculation so that the primary
# physical comparison is saved even if a later secondary calculation stops.
RESULT_DIR <- file.path(
  "results",
  "stage2_scenario_B_full_vs_cut"
)

FIGURE_DIR <- file.path(
  RESULT_DIR,
  "figures"
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\nStage-2 Scenario B output directory:\n")
cat(normalizePath(RESULT_DIR, mustWork = FALSE), "\n")

# Scenario B is the confounding stress test.
# The inference model and priors are unchanged from Scenario A.
# Only the DGP discrepancy is changed.
# This is the central fixed-m experiment for showing the effect of cutting
# downstream feedback on the physical calibration parameters.

###############################################################################
## 1. LOAD THE EXACT STAGE-1 DATASET AND FULL POSTERIOR DRAWS
###############################################################################

stage1_result_dir <- file.path(
  "results",
  "stage1_scenario_B_confounding"
)

data_file <- file.path(
  stage1_result_dir,
  "simulated_data.csv"
)

full_draw_file <- file.path(
  stage1_result_dir,
  "pooled_draws.rds"
)

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

data_stage1 <- read.csv(
  data_file,
  stringsAsFactors = FALSE
)

full_fixed_draws <- readRDS(
  full_draw_file
)
seed = 2999L
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

# Save the PRIMARY Stage-2 physical results immediately.
# These files remain available even if a later downstream diagnostic stops.
saveRDS(
  cut_physical_draws,
  file.path(RESULT_DIR, "cut_physical_draws.rds")
)

write.csv(
  cut_physical_summary,
  file.path(RESULT_DIR, "cut_physical_summary.csv"),
  row.names = TRUE
)

write.csv(
  physical_comparison,
  file.path(RESULT_DIR, "full_cut_physical_comparison.csv"),
  row.names = FALSE
)

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
  n_reconstruct = 300L,
  seed = 3111L
)

delta_cut <- reconstruct_delta_fixed(
  draws = cut_joint_draws,
  y = dat$y,
  t_raw = dat$t_raw,
  m = truth$m,
  n_reconstruct = min(300L, nrow(cut_joint_draws)),
  seed = 3222L
)

idx_D <- (truth$m + 1):nrow(dat$data)

full_delta_ISE <- mean(
  (
    delta_full$mean[idx_D] -
      dat$data$delta_true[idx_D]
  )^2
)

cut_delta_ISE <- mean(
  (
    delta_cut$mean[idx_D] -
      dat$data$delta_true[idx_D]
  )^2
)

delta_error <- data.frame(
  method = c("Full", "Cut"),
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
## 10. SAVE FIGURES AND SHOW THEM IN R
###############################################################################

figure_dir <- FIGURE_DIR

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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

cat("\n============================================================\n")
cat("STAGE 2 - SCENARIO B COMPLETE\n")
cat("============================================================\n")
cat("Results directory:\n", normalizePath(result_dir), "\n")
cat("Figures directory:\n", normalizePath(figure_dir), "\n")
cat("\nThe MAIN Stage-2 physical comparison uses", N_CUT_PHYSICAL,
    "direct Cut draws.\n")
cat("The downstream Cut uses", N_CUT_OUTER,
    "outer draws; each outer draw uses an inner MH run of",
    INNER_ITER, "iterations with", INNER_BURN, "burn-in.\n")
cat("In RStudio, use the Plots pane arrows to browse the displayed figures.\n")

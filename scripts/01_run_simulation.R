run_one_scenario <- function(scenario, scenario_index) {
  
  scenario_number <- as.integer(
    SCENARIO_NUMBERS[scenario]
  )
  
  scenario_label <- as.character(
    SCENARIO_LABELS[scenario]
  )
  
  cat(
    "\nScenario ",
    scenario_number,
    ": ",
    scenario_label,
    "\n",
    sep = ""
  )
  
  set.seed(SEED_DATA[scenario_index])
  
  sim <- simulate_one_dataset(scenario)
  
  y <- sim$y
  delta_true <- sim$delta_true
  
  y_P <- y[I_P]
  t_P <- t[I_P]
  
  y_D <- y[I_D]
  t_D <- t[I_D]
  
  delta_true_D <- delta_true[I_D]
  
  # -------------------------------------------------------------------------
  # Module P
  # -------------------------------------------------------------------------
  
  set.seed(SEED_MODULE_P[scenario_index])
  
  phys_chain <- mcmc_module_P(
    y_P = y_P,
    t_P = t_P,
    n_keep = N_KEEP_PHYS,
    burn = BURN_PHYS,
    init = init_phys,
    verbose = FALSE
  )
  
  #moduleP_ess <- safe_effective_size(phys_chain)
  
  #names(moduleP_ess) <- colnames(phys_chain)
  
  # Same physical posterior draws are used for the two Cut versions.
  set.seed(SEED_PHYS_IDS[scenario_index])
  
  phys_ids <- sample.int(
    nrow(phys_chain),
    size = N_OUTER_PHYS,
    replace = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Full Mixture
  # -------------------------------------------------------------------------
  
  set.seed(SEED_FULL[scenario_index])
  
  full_res <- mcmc_full_mixture(
    y = y,
    t_all = t,
    n_keep = N_KEEP_FULL,
    burn = BURN_FULL,
    init_phys_local = init_phys,
    init_mix_local = init_mix,
    verbose = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Conditional Cut Mixture
  # -------------------------------------------------------------------------
  
  set.seed(SEED_CUT_MIX[scenario_index])
  
  cut_mix_res <- mcmc_cut_mixture_conditional(
    y_D = y_D,
    t_D = t_D,
    phys_chain = phys_chain,
    phys_ids = phys_ids,
    n_keep_D = N_KEEP_D_PER_PHYS,
    burn_D = BURN_D_PER_PHYS,
    init = init_mix,
    verbose = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Conditional Cut No-Mixture
  # -------------------------------------------------------------------------
  
  set.seed(SEED_CUT_NOMIX[scenario_index])
  
  cut_nomix_res <- mcmc_cut_nomix_conditional(
    y_D = y_D,
    t_D = t_D,
    phys_chain = phys_chain,
    phys_ids = phys_ids,
    n_keep_D = N_KEEP_D_PER_PHYS,
    burn_D = BURN_D_PER_PHYS,
    init = init_D,
    verbose = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Plug-in No-Mixture
  # -------------------------------------------------------------------------
  
  set.seed(SEED_PLUGIN[scenario_index])
  
  plugin_res <- mcmc_plugin_nomix(
    y_D = y_D,
    t_D = t_D,
    phys_chain = phys_chain,
    n_keep = N_KEEP_PLUGIN,
    burn = BURN_PLUGIN,
    init = init_D,
    verbose = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Summaries
  # -------------------------------------------------------------------------
  
  physical_summary <- rbind(
    physical_summary_from_chain(
      full_res$chain,
      "Full Mixture"
    ),
    physical_summary_from_chain(
      phys_chain,
      "Conditional Cut Mixture"
    ),
    physical_summary_from_chain(
      phys_chain,
      "Conditional Cut No-Mixture"
    ),
    physical_summary_from_point(
      plugin_res$physical_hat,
      "Plug-in No-Mixture"
    )
  )
  
  physical_summary$scenario <- scenario
  physical_summary$replicate <- REPLICATE_ID
  
  physical_summary$method <- factor(
    physical_summary$method,
    levels = METHOD_LEVELS
  )
  
  hyper_summary <- rbind(
    hyper_summary_from_chain(
      full_res$chain,
      "Full Mixture"
    ),
    hyper_summary_from_chain(
      cut_mix_res$hyper,
      "Conditional Cut Mixture"
    ),
    hyper_summary_from_chain(
      cut_nomix_res$hyper,
      "Conditional Cut No-Mixture"
    ),
    hyper_summary_from_chain(
      plugin_res$hyper,
      "Plug-in No-Mixture"
    )
  )
  
  alpha_summary <- rbind(
    alpha_summary_from_chain(
      full_res$chain,
      "Full Mixture"
    ),
    alpha_summary_from_chain(
      cut_mix_res$hyper,
      "Conditional Cut Mixture"
    )
  )
  
  full_delta_D <-
    full_res$delta_eff[
      ,
      I_D,
      drop = FALSE
    ]
  
  delta_full <- summarize_delta_draws(
    delta_draws = full_delta_D,
    method = "Full Mixture",
    scenario = scenario,
    delta_true_D = delta_true_D
  )
  
  delta_cut_mix <- summarize_delta_draws(
    delta_draws = cut_mix_res$delta_eff,
    method = "Conditional Cut Mixture",
    scenario = scenario,
    delta_true_D = delta_true_D
  )
  
  delta_cut_nomix <- summarize_delta_draws(
    delta_draws = cut_nomix_res$delta,
    method = "Conditional Cut No-Mixture",
    scenario = scenario,
    delta_true_D = delta_true_D
  )
  
  delta_plugin <- summarize_delta_draws(
    delta_draws = plugin_res$delta,
    method = "Plug-in No-Mixture",
    scenario = scenario,
    delta_true_D = delta_true_D
  )
  
  delta_summary <- rbind(
    delta_full,
    delta_cut_mix,
    delta_cut_nomix,
    delta_plugin
  )
  
  delta_summary$method <- factor(
    delta_summary$method,
    levels = METHOD_LEVELS
  )
  
  delta_metrics <- rbind(
    delta_metrics_from_summary(
      delta_full
    ),
    delta_metrics_from_summary(
      delta_cut_mix
    ),
    delta_metrics_from_summary(
      delta_cut_nomix
    ),
    delta_metrics_from_summary(
      delta_plugin
    )
  )
  
  zeta_df <- rbind(
    data.frame(
      method = "Full Mixture",
      index = I_D,
      t = t_D,
      prob_zeta2 = full_res$zeta2_prob[I_D]
    ),
    data.frame(
      method = "Conditional Cut Mixture",
      index = I_D,
      t = t_D,
      prob_zeta2 = cut_mix_res$zeta2_prob
    )
  )
  
  plots <- list(
    g = plot_physical_parameter(
      physical_summary,
      parameter = "g",
      title = paste0(
        "Scénario ",
        scenario_number,
        " - g"
      ),
      y_label = "g"
    ),
    h0 = plot_physical_parameter(
      physical_summary,
      parameter = "h0",
      title = paste0(
        "Scénario ",
        scenario_number,
        " - h0"
      ),
      y_label = "h0"
    ),
    delta = plot_delta(
      delta_summary,
      title = paste0(
        "Scénario ",
        scenario_number,
        " - discrépance"
      )
    )
  )
  
  return(
    list(
      scenario = scenario,
      data = data.frame(
        t = t,
        y = y,
        f_true = sim$f_true,
        delta_true = delta_true,
        module = ifelse(
          seq_len(n) <= m_phys,
          "P",
          "D"
        )
      ),
      phys_chain = phys_chain,
      full = full_res,
      cut_mixture = cut_mix_res,
      cut_nomix = cut_nomix_res,
      plugin = plugin_res,
      physical_summary = physical_summary,
      hyper_summary = hyper_summary,
      alpha_summary = alpha_summary,
      delta_summary = delta_summary,
      delta_metrics = delta_metrics,
      zeta = zeta_df,
      #moduleP_ESS = moduleP_ess,
      plots = plots
    )
  )
}

###############################################################################
# Run
###############################################################################

# For all six scenarios:
SCENARIOS_TO_RUN <- SCENARIOS

# For a quick check of only one scenario, use for example:
# SCENARIOS_TO_RUN <- "moderate_smooth_positive"

results <- list()

for (scenario in SCENARIOS_TO_RUN) {
  
  scenario_index <- match(
    scenario,
    SCENARIOS
  )
  
  results[[scenario]] <- run_one_scenario(
    scenario = scenario,
    scenario_index = scenario_index
  )
  
  print(
    results[[scenario]]$physical_summary
  )
  
  print(
    results[[scenario]]$delta_metrics
  )
  
  print(
    results[[scenario]]$plots$g
  )
  
  print(
    results[[scenario]]$plots$delta
  )
}

# Example:
# result3 <- results[["moderate_smooth_positive"]]
# result3$physical_summary
# result3$delta_metrics
# print(result3$plots$g)
# print(result3$plots$delta)
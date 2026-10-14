###############################################################################
# ARTICLE 2 - STAGE 2
# CUT BAYESIAN MODEL WITH A FIXED, KNOWN CHANGE POINT
#
# Stage 1:
#   Full posterior. Both Y_P and Y_D can update theta and lambda^2.
#
# Stage 2:
#   Cut posterior.
#   Module P:
#       (theta, lambda^2) are learned from Y_P ONLY.
#   Module D:
#       (k, gamma_delta) are learned from Y_D conditional on the upstream draw.
#
# Fixed-m Cut target:
#   pi_cut^m(phi, eta | Y) = pi_P(phi | Y_P, m) * pi_D(eta | Y_D, phi, m),
#
#   phi = (theta, lambda^2)
#   eta = (k, gamma_delta)
###############################################################################

source("stage1_core_full_fixed_m.R")

###############################################################################
## 1. MODULE P: PHYSICAL PARAMETERS AND COMMON VARIANCE
###############################################################################
#
#   A_P = X_P' X_P
#   B_P = X_P' Y_P
#   theta_hat_P = A_P^{-1} B_P
#   S_P = ||Y_P - X_P theta_hat_P||^2
#
#   lambda^2 | Y_P,m ~ IG(m/2, S_P/2)
#
#   theta | lambda^2,Y_P,m
#     ~ N(theta_hat_P, lambda^2 A_P^{-1})
#
module_P_stats <- function(y, t, t_min, t_range, m, prior_type = PRIOR_TYPE) {
  X <- make_X(t, t_min, t_range)
  y_P <- y[1:m]
  X_P <- X[1:m, , drop = FALSE]
  
  if (qr(X_P)$rank < ncol(X_P)) {
    stop("Module P design matrix is not full rank for this m.")
  }
  
  A_data <- t(X_P) %*% X_P
  B_data <- t(X_P) %*% y_P
  
  # ---------------------------------------------------------------------------
  if (prior_type == "jeffreys") {
    A_inv <- safe_solve(A_data)
    
    theta_hat <- as.vector(A_inv %*% B_data)
    names(theta_hat) <- c("g", "h0")
    
    residual <- y_P - as.vector(X_P %*% theta_hat)
    SSE <- sum(residual^2)
    
    if (!is.finite(SSE) || SSE <= 0) {
      stop("Non-positive Module-P residual sum of squares.")
    }
    
    return(list(
      theta_hat = theta_hat,
      A_inv = A_inv,
      A = A_data,
      B = B_data,
      SSE = SSE,
      lambda_shape = m / 2,
      lambda_rate = SSE / 2
    ))
  }
  
  # ---------------------------------------------------------------------------
  # OPTIONAL PROPER-PRIOR SENSITIVITY ONLY
  # ---------------------------------------------------------------------------
  if (prior_type == "proper") {
    if (!exists("mu0") || !exists("V0") ||
        !exists("a0") || !exists("b0")) {
      stop("Proper prior selected, but mu0, V0, a0, b0 are not defined.")
    }
    
    V0_inv <- safe_solve(V0)
    
    A <- A_data + V0_inv
    A_inv <- safe_solve(A)
    B <- B_data + V0_inv %*% mu0
    
    theta_hat <- as.vector(A_inv %*% B)
    names(theta_hat) <- c("g", "h0")
    
    S <- sum(y_P^2) + as.numeric(t(mu0) %*% V0_inv %*% mu0) -
      as.numeric(t(B) %*% A_inv %*% B)
    
    return(list(
      theta_hat = theta_hat,
      A_inv = A_inv,
      A = A,
      B = B,
      SSE = S,
      lambda_shape = a0 + m / 2,
      lambda_rate = b0 + S / 2
    ))
  }
  
  stop("Unknown prior_type. Use 'jeffreys' or 'proper'.")
}

# One exact composition draw from the normalized Module-P posterior:
# first lambda^2, then theta conditional on lambda^2.
sample_module_P <- function(y, t, t_min, t_range, m, prior_type = PRIOR_TYPE) {
  st <- module_P_stats(
    y = y,
    t = t,
    t_min = t_min,
    t_range = t_range,
    m = m,
    prior_type = prior_type
  )
  lambda2 <- rinvgamma(1, shape = st$lambda_shape, rate = st$lambda_rate)
  theta <- as.vector(rmvnorm(1, mean = st$theta_hat, sigma = lambda2 * st$A_inv))
  names(theta) <- c("g", "h0")
  
  list(
    theta = theta,
    lambda2 = lambda2
  )
}

# These are DIRECT posterior draws. There is no MCMC chain for
# (theta, lambda^2) in Module P because the conditional distributions are known.
run_cut_physical_fixed <- function(y, t_raw, m, n_draws = 10000, prior_type = PRIOR_TYPE, seed = 999) {
  set.seed(seed)
  
  prep <- prepare_time(t_raw)
  
  out <- data.frame(
    draw = 1:n_draws,
    lambda2 = numeric(n_draws),
    g = numeric(n_draws),
    h0 = numeric(n_draws)
  )
  
  for (r in 1:n_draws) {
    up <- sample_module_P(
      y = y,
      t = prep$t,
      t_min = prep$t_min,
      t_range = prep$t_range,
      m = m,
      prior_type = prior_type
    )
    
    out$lambda2[r] <- up$lambda2
    out$g[r] <- up$theta["g"]
    out$h0[r] <- up$theta["h0"]
  }
  
  out
}

###############################################################################
## 2. MODULE D: GP HYPERPARAMETERS CONDITIONAL ON THE UPSTREAM DRAW
###############################################################################

# Conditional downstream log target, up to an additive constant:
#
#   -0.5 log|C_D|
#   - [r_D(theta)' C_D^{-1} r_D(theta)] / (2 lambda^2)
#
# theta and lambda^2 are fixed conditioning values here.

log_downstream_target <- function(y, t, t_min, t_range, m, theta, lambda2, k, gamma_delta) {
  if (k <= k_lower || k >= k_upper) return(-Inf)
  if (gamma_delta <= gamma_lower || gamma_delta >= gamma_upper) return(-Inf)
  
  n <- length(y)
  X <- make_X(t, t_min, t_range)
  idx_D <- (m + 1):n
  
  y_D <- y[idx_D]
  X_D <- X[idx_D, , drop = FALSE]
  t_D <- t[idx_D]
  
  C <- make_C(t_D, k, gamma_delta)
  C_inv <- safe_solve(C)
  
  residual <- y_D - as.vector(X_D %*% theta)
  
  q <- as.numeric(t(residual) %*% C_inv %*% residual)
  
  # Uniform priors are constant inside their supports.
  -0.5 * safe_logdet(C) - q / (2 * lambda2)
}

# Truncated-normal MH update for k.
update_k_downstream <- function(k_current, y, t, t_min, t_range, m, theta, lambda2, gamma_delta, proposal_sd = 0.02) {
  k_prop <- rtruncnorm(1, a = k_lower, b = k_upper, mean = k_current, sd = proposal_sd)
  
  log_current <- log_downstream_target(y, t, t_min, t_range, m, theta, lambda2, k_current, gamma_delta)
  
  log_prop <- log_downstream_target(y, t, t_min, t_range, m, theta, lambda2, k_prop, gamma_delta)
  
  log_q_reverse <- log(dtruncnorm(k_current, a = k_lower, b = k_upper, mean = k_prop, sd = proposal_sd))
  
  log_q_forward <- log(dtruncnorm(k_prop, a = k_lower, b = k_upper, mean = k_current, sd = proposal_sd))
  
  log_ratio <- (log_prop - log_current) + (log_q_reverse - log_q_forward)
  
  if (
    is.finite(log_ratio) &&
    log(runif(1)) < min(0, log_ratio)
  ) {
    list(value = k_prop, accepted = 1)
  } else {
    list(value = k_current, accepted = 0)
  }
}

# Truncated-normal MH update for gamma_delta.
update_gamma_downstream <- function(gamma_current, y, t, t_min, t_range, m, theta, lambda2, k, proposal_sd = 0.15) {
  gamma_prop <- rtruncnorm(1, a = gamma_lower, b = gamma_upper, mean = gamma_current, sd = proposal_sd)
  
  log_current <- log_downstream_target(y, t, t_min, t_range, m, theta, lambda2, k, gamma_current)
  
  log_prop <- log_downstream_target(y, t, t_min, t_range, m, theta, lambda2, k, gamma_prop)
  
  log_q_reverse <- log(dtruncnorm(gamma_current, a = gamma_lower, b = gamma_upper, mean = gamma_prop, sd = proposal_sd))
  
  log_q_forward <- log(dtruncnorm(gamma_prop, a = gamma_lower, b = gamma_upper, mean = gamma_current, sd = proposal_sd))
  
  log_ratio <- (log_prop - log_current) + (log_q_reverse - log_q_forward)
  
  if (
    is.finite(log_ratio) &&
    log(runif(1)) < min(0, log_ratio)
  ) {
    list(value = gamma_prop, accepted = 1)
  } else {
    list(value = gamma_current, accepted = 0)
  }
}

# A short inner MCMC is run for the CURRENT upstream draw
# (theta, lambda^2), and one retained state is selected at random.
# This is the computationally expensive part of the full Cut composition.
sample_downstream <- function(y, t, t_min, t_range, m, theta, lambda2, n_iter = 1200,
                              burn = 400, init_k = 0.1, init_gamma = 0.3, proposal_sd_k = 0.02, proposal_sd_gamma = 0.15) {
  
  if (burn >= n_iter) stop("burn must be smaller than n_iter.")
  
  k <- init_k
  gamma_delta <- init_gamma
  
  keep_k <- numeric(n_iter - burn)
  keep_gamma <- numeric(n_iter - burn)
  
  accept_k <- 0
  accept_gamma <- 0
  j <- 0
  
  for (iter in 1:n_iter) {
    out_k <- update_k_downstream(
      k_current = k,
      y = y,
      t = t,
      t_min = t_min,
      t_range = t_range,
      m = m,
      theta = theta,
      lambda2 = lambda2,
      gamma_delta = gamma_delta,
      proposal_sd = proposal_sd_k
    )
    
    k <- out_k$value
    accept_k <- accept_k + out_k$accepted
    
    out_gamma <- update_gamma_downstream(
      gamma_current = gamma_delta,
      y = y,
      t = t,
      t_min = t_min,
      t_range = t_range,
      m = m,
      theta = theta,
      lambda2 = lambda2,
      k = k,
      proposal_sd = proposal_sd_gamma
    )
    
    gamma_delta <- out_gamma$value
    accept_gamma <- accept_gamma + out_gamma$accepted
    
    if (iter > burn) {
      j <- j + 1
      keep_k[j] <- k
      keep_gamma[j] <- gamma_delta
    }
  }
  
  pick <- sample.int(length(keep_k), 1)
  
  list(
    k = keep_k[pick],
    gamma_delta = keep_gamma[pick],
    acceptance = c(
      k = accept_k / n_iter,
      gamma_delta = accept_gamma / n_iter
    )
  )
}

###############################################################################
## 3. COMPLETE FIXED-m CUT COMPOSITION SAMPLER
###############################################################################

# For each outer draw:
#   1) draw (theta, lambda^2) exactly from Module P using Y_P only;
#   2) conditional on that draw, approximate one draw of (k, gamma_delta)
#      from Module D using an inner MH chain.
#
# The physical block of the Cut is therefore exact composition sampling.
# The GP-hyperparameter block is approximate because of the finite inner MCMC.
run_cut_fixed_m <- function(
    y,
    t_raw,
    m,
    n_outer = 300,
    inner_iter = 1200,
    inner_burn = 400,
    proposal_sd_k = 0.02,
    proposal_sd_gamma = 0.15,
    prior_type = PRIOR_TYPE,
    seed = 1002,
    verbose = TRUE
) {
  set.seed(seed)
  
  prep <- prepare_time(t_raw)
  
  out <- data.frame(
    draw = 1:n_outer,
    m = rep(m, n_outer),
    k = numeric(n_outer),
    gamma_delta = numeric(n_outer),
    lambda2 = numeric(n_outer),
    g = numeric(n_outer),
    h0 = numeric(n_outer),
    acc_k = numeric(n_outer),
    acc_gamma = numeric(n_outer)
  )
  
  for (r in 1:n_outer) {
    # Module P: upstream draw from Y_P only.
    up <- sample_module_P(
      y = y,
      t = prep$t,
      t_min = prep$t_min,
      t_range = prep$t_range,
      m = m,
      prior_type = prior_type
    )
    
    # Module D: downstream draw conditional on the upstream draw.
    down <- sample_downstream(
      y = y,
      t = prep$t,
      t_min = prep$t_min,
      t_range = prep$t_range,
      m = m,
      theta = up$theta,
      lambda2 = up$lambda2,
      n_iter = inner_iter,
      burn = inner_burn,
      proposal_sd_k = proposal_sd_k,
      proposal_sd_gamma = proposal_sd_gamma
    )
    
    out$k[r] <- down$k
    out$gamma_delta[r] <- down$gamma_delta
    out$lambda2[r] <- up$lambda2
    out$g[r] <- up$theta["g"]
    out$h0[r] <- up$theta["h0"]
    out$acc_k[r] <- down$acceptance["k"]
    out$acc_gamma[r] <- down$acceptance["gamma_delta"]
    
    if (
      verbose &&
      r %% max(1, floor(n_outer / 10)) == 0
    ) {
      cat("Fixed-m Cut outer draw", r, "of", n_outer, "\n")
    }
  }
  
  out
}

###############################################################################
## 4. FULL vs CUT 
###############################################################################

make_full_cut_physical_comparison <- function(
    full_draws,
    cut_physical_draws,
    truth
) {
  full_summary <- make_parameter_summary(full_draws)
  
  cut_summary <- rbind(
    g = summarize_parameter(cut_physical_draws$g),
    h0 = summarize_parameter(cut_physical_draws$h0),
    lambda2 = summarize_parameter(cut_physical_draws$lambda2)
  )
  
  truth_vec <- c(
    g = unname(truth$theta["g"]),
    h0 = unname(truth$theta["h0"]),
    lambda2 = truth$lambda2
  )
  
  out <- data.frame(
    parameter = c("g", "h0", "lambda2"),
    truth = truth_vec[c("g", "h0", "lambda2")],
    full_mean = c(
      full_summary["g", "mean"],
      full_summary["h0", "mean"],
      full_summary["lambda2", "mean"]
    ),
    full_sd = c(
      full_summary["g", "sd"],
      full_summary["h0", "sd"],
      full_summary["lambda2", "sd"]
    ),
    full_q025 = c(
      full_summary["g", "q025"],
      full_summary["h0", "q025"],
      full_summary["lambda2", "q025"]
    ),
    full_q975 = c(
      full_summary["g", "q975"],
      full_summary["h0", "q975"],
      full_summary["lambda2", "q975"]
    ),
    cut_mean = c(
      cut_summary["g", "mean"],
      cut_summary["h0", "mean"],
      cut_summary["lambda2", "mean"]
    ),
    cut_sd = c(
      cut_summary["g", "sd"],
      cut_summary["h0", "sd"],
      cut_summary["lambda2", "sd"]
    ),
    cut_q025 = c(
      cut_summary["g", "q025"],
      cut_summary["h0", "q025"],
      cut_summary["lambda2", "q025"]
    ),
    cut_q975 = c(
      cut_summary["g", "q975"],
      cut_summary["h0", "q975"],
      cut_summary["lambda2", "q975"]
    ),
    row.names = NULL
  )
  
  out$full_error <- out$full_mean - out$truth
  out$cut_error <- out$cut_mean - out$truth
  out$full_abs_error <- abs(out$full_error)
  out$cut_abs_error <- abs(out$cut_error)
  out$full_width <- out$full_q975 - out$full_q025
  out$cut_width <- out$cut_q975 - out$cut_q025
  
  # Full minus Cut.
  # For g this is a direct descriptive measure of the feedback-related shift
  # in the fixed-m experiment, not a causal estimand.
  out$feedback_shift <- out$full_mean - out$cut_mean
  
  out
}

###############################################################################
## 5. PLOT
###############################################################################

plot_full_cut_physical <- function(
    full_draws,
    cut_draws,
    truth,
    scenario_label = "Stage 2"
) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(mfrow = c(1, 3))
  
  params <- c("g", "h0", "lambda2")
  truth_values <- c(
    g = unname(truth$theta["g"]),
    h0 = unname(truth$theta["h0"]),
    lambda2 = truth$lambda2
  )
  
  xlabels <- list(
    "g",
    expression(h[0]),
    expression(lambda^2)
  )
  
  for (j in seq_along(params)) {
    v <- params[j]
    
    xlim <- range(
      full_draws[[v]],
      cut_draws[[v]],
      finite = TRUE
    )
    
    plot(
      density(full_draws[[v]]),
      lwd = 2,
      xlim = xlim,
      xlab = xlabels[[j]],
      main = paste(scenario_label, "-", v)
    )
    
    lines(
      density(cut_draws[[v]]),
      lty = 2,
      lwd = 2
    )
    
    abline(
      v = truth_values[v],
      lty = 3,
      lwd = 2
    )
    
    legend(
      "topright",
      legend = c("Full", "Cut", "Truth"),
      lty = c(1, 2, 3),
      lwd = c(2, 2, 2),
      bty = "n"
    )
  }
}

# Secondary comparison for the GP hyperparameters.
# For Scenario B, truth_gp must be NULL because k and gamma_delta do not have
# DGP truth values under the deterministic confounding discrepancy.
plot_full_cut_gp <- function(
    full_draws,
    cut_joint_draws,
    scenario_label = "Stage 2",
    truth_gp = NULL
) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(mfrow = c(1, 2))
  
  params <- c("k", "gamma_delta")
  xlabels <- list(
    "k",
    expression(gamma[delta])
  )
  
  for (j in seq_along(params)) {
    v <- params[j]
    
    xlim <- range(
      full_draws[[v]],
      cut_joint_draws[[v]],
      finite = TRUE
    )
    
    plot(
      density(full_draws[[v]]),
      lwd = 2,
      xlim = xlim,
      xlab = xlabels[[j]],
      main = paste(scenario_label, "-", v)
    )
    
    lines(
      density(cut_joint_draws[[v]]),
      lty = 2,
      lwd = 2
    )
    
    if (!is.null(truth_gp) && is.finite(truth_gp[[v]])) {
      abline(
        v = truth_gp[[v]],
        lty = 3,
        lwd = 2
      )
      
      legend(
        "topright",
        legend = c("Full", "Cut", "Truth"),
        lty = c(1, 2, 3),
        lwd = c(2, 2, 2),
        bty = "n"
      )
    } else {
      legend(
        "topright",
        legend = c("Full", "Cut"),
        lty = c(1, 2),
        lwd = c(2, 2),
        bty = "n"
      )
    }
  }
}

# Direct comparison of posterior discrepancy means.
plot_full_cut_delta_means <- function(
    dat,
    delta_full,
    delta_cut,
    scenario_label = "Stage 2"
) {
  true_boundary_time <- mean(
    dat$data$t_raw[
      c(dat$truth$m, dat$truth$m + 1)
    ]
  )
  
  ylim <- range(
    dat$data$delta_true,
    delta_full$mean,
    delta_cut$mean,
    finite = TRUE
  )
  
  plot(
    dat$data$t_raw,
    dat$data$delta_true,
    type = "l",
    lwd = 2,
    ylim = ylim,
    xlab = "Time (seconds)",
    ylab = expression(delta(t)),
    main = paste(scenario_label, "- discrepancy posterior means")
  )
  
  lines(
    dat$data$t_raw,
    delta_full$mean,
    lty = 2,
    lwd = 2
  )
  
  lines(
    dat$data$t_raw,
    delta_cut$mean,
    lty = 3,
    lwd = 2
  )
  
  abline(
    v = true_boundary_time,
    lty = 4,
    lwd = 1
  )
  
  legend(
    "topleft",
    legend = c(
      "True discrepancy",
      "Full posterior mean",
      "Cut posterior mean",
      "True boundary"
    ),
    lty = c(1, 2, 3, 4),
    lwd = c(2, 2, 2, 1),
    bty = "n"
  )
}

# Two-panel reconstruction with pointwise credible limits.
plot_full_cut_delta_intervals <- function(
    dat,
    delta_full,
    delta_cut,
    scenario_label = "Stage 2"
) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  true_boundary_time <- mean(
    dat$data$t_raw[
      c(dat$truth$m, dat$truth$m + 1)
    ]
  )
  
  common_ylim <- range(
    dat$data$delta_true,
    delta_full$q025,
    delta_full$q975,
    delta_cut$q025,
    delta_cut$q975,
    finite = TRUE
  )
  
  par(mfrow = c(1, 2))
  
  plot(
    dat$data$t_raw,
    dat$data$delta_true,
    type = "l",
    lwd = 2,
    ylim = common_ylim,
    xlab = "Time (seconds)",
    ylab = expression(delta(t)),
    main = paste(scenario_label, "- Full")
  )
  lines(dat$data$t_raw, delta_full$mean, lty = 2, lwd = 2)
  lines(dat$data$t_raw, delta_full$q025, lty = 3)
  lines(dat$data$t_raw, delta_full$q975, lty = 3)
  abline(v = true_boundary_time, lty = 4)
  
  plot(
    dat$data$t_raw,
    dat$data$delta_true,
    type = "l",
    lwd = 2,
    ylim = common_ylim,
    xlab = "Time (seconds)",
    ylab = expression(delta(t)),
    main = paste(scenario_label, "- Cut")
  )
  lines(dat$data$t_raw, delta_cut$mean, lty = 2, lwd = 2)
  lines(dat$data$t_raw, delta_cut$q025, lty = 3)
  lines(dat$data$t_raw, delta_cut$q975, lty = 3)
  abline(v = true_boundary_time, lty = 4)
}

###############################################################################
## 6. OUTPUT
###############################################################################

save_stage2_outputs <- function(
    result_dir,
    full_draws,
    cut_physical_draws,
    cut_joint_draws,
    physical_comparison,
    cut_physical_summary,
    cut_joint_summary,
    downstream_acceptance,
    delta_full,
    delta_cut,
    delta_error
) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(
    cut_physical_draws,
    file.path(result_dir, "cut_physical_draws.rds")
  )
  
  saveRDS(
    cut_joint_draws,
    file.path(result_dir, "cut_joint_draws.rds")
  )
  
  saveRDS(
    delta_full,
    file.path(result_dir, "delta_full_reconstruction.rds")
  )
  
  saveRDS(
    delta_cut,
    file.path(result_dir, "delta_cut_reconstruction.rds")
  )
  
  write.csv(
    make_parameter_summary(full_draws),
    file.path(result_dir, "full_parameter_summary.csv"),
    row.names = TRUE
  )
  
  write.csv(
    cut_physical_summary,
    file.path(result_dir, "cut_physical_summary.csv"),
    row.names = TRUE
  )
  
  write.csv(
    cut_joint_summary,
    file.path(result_dir, "cut_joint_summary.csv"),
    row.names = TRUE
  )
  
  write.csv(
    physical_comparison,
    file.path(result_dir, "full_cut_physical_comparison.csv"),
    row.names = FALSE
  )
  
  write.csv(
    downstream_acceptance,
    file.path(result_dir, "cut_downstream_acceptance.csv"),
    row.names = FALSE
  )
  
  write.csv(
    delta_error,
    file.path(result_dir, "delta_error_metrics.csv"),
    row.names = FALSE
  )
}


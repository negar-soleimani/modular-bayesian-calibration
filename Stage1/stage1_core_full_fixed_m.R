###############################################################################
# STAGE 1
# FULL BAYESIAN MODEL WITH A FIXED, KNOWN CHANGE POINT
#
# Model:
#   i <= m : y_i = f(t_i, theta) + epsilon_i
#   i >  m : y_i = f(t_i, theta) + delta_i + epsilon_i
#
#   epsilon_i ~ N(0, lambda^2)
#   delta_D  ~ GP(0, (lambda^2 / k) R_gamma)
#
# Main priors:
#   pi(theta, lambda^2) proportional to (lambda^2)^(-(d/2 + 1))
#   k           ~ Uniform(0.015, 1)
#   gamma_delta ~ Uniform(0.1, 1)
#
# IMPORTANT:
# - m is FIXED and KNOWN in Stage 1.
# - delta_D is analytically integrated out in the main sampler.
# - delta_D can be reconstructed after inference from its Gaussian conditional.
###############################################################################

library(mvtnorm)
library(invgamma)
library(truncnorm)
library(coda)

###############################################################################
## 1. PRIOR SETTINGS
###############################################################################

PRIOR_TYPE <- "jeffreys"
d <- 2

k_lower <- 0.015
k_upper <- 1.0

gamma_lower <- 0.10
gamma_upper <- 1.0

# -----------------------------------------------------------------------------
# OPTIONAL PROPER-PRIOR SENSITIVITY ONLY
# -----------------------------------------------------------------------------
# PRIOR_TYPE <- "proper"
# mu0 <- c(g = 0, h0 = 0)
# V0 <- diag(c(1e4, 1e4))
# a0 <- 0.01
# b0 <- 0.01
#
# theta | lambda^2 ~ N(mu0, lambda^2 V0)
# lambda^2         ~ IG(a0, b0)
# -----------------------------------------------------------------------------

###############################################################################
## 2. MODEL FUNCTIONS
###############################################################################

safe_solve <- function(M, jitter = 1e-10) {
  M <- 0.5 * (M + t(M))
  solve(M + jitter * diag(nrow(M)))
}

safe_logdet <- function(M, jitter = 1e-10) {
  M <- 0.5 * (M + t(M))
  M <- M + jitter * diag(nrow(M))
  out <- determinant(M, logarithm = TRUE)
  if (out$sign <= 0) return(-Inf)
  as.numeric(out$modulus)
}

# GP covariance.
# IMPORTANT: in the current implementation the input t supplied here is the
# normalized time u in [0,1], not the physical time in seconds.
GP_covariance <- function(t, sigma_sq_delta, gamma_delta) {
  outer(t, t, function(ti, tj)
    sigma_sq_delta * exp(-abs(ti - tj) / gamma_delta)
  )
}

balldropg <- function(t, theta, t_min, t_range) {
  g <- theta[1]
  h0 <- theta[2]
  
  t_physical <- t * t_range + t_min
  h <- h0 - 0.5 * g * t_physical^2
  
  h[h < 0] <- 0
  as.vector(h)
}

# Linear design matrix used by the Stage-1 inference.
make_X <- function(t, t_min, t_range) {
  t_physical <- t * t_range + t_min
  cbind(
    g = -0.5 * t_physical^2,
    h0 = 1
  )
}

# Normalize time to [0,1] for the GP correlation calculation.
prepare_time <- function(t_raw) {
  t_min <- min(t_raw)
  t_range <- max(t_raw) - min(t_raw)
  
  if (t_range <= 0) {
    stop("All time points are identical.")
  }
  
  t <- (t_raw - t_min) / t_range
  
  list(
    t_raw = t_raw,
    t = t,
    t_min = t_min,
    t_range = t_range
  )
}

# After integrating delta_D:
#   C_D = I + R_gamma / k
# and
#   Y_D | theta,lambda2,k,gamma,m ~ N(X_D theta, lambda2 C_D).
make_C <- function(t_D, k, gamma_delta) {
  R_gamma <- GP_covariance(
    t_D,
    sigma_sq_delta = 1,
    gamma_delta = gamma_delta
  )
  
  diag(length(t_D)) + R_gamma / k
}

###############################################################################
## 3. FIXED-m FULL COLLAPSED CALCULATION
###############################################################################
#
#   A_m = X_P' X_P + X_D' C_D^{-1} X_D
#   B_m = X_P' Y_P + X_D' C_D^{-1} Y_D
#   theta_hat_m = A_m^{-1} B_m
#
#   S_m = r_P' r_P + r_D' C_D^{-1} r_D
#
# Under the baseline prior:
#   lambda2 | Y,k,gamma,m ~ IG(n/2, S_m/2)
#   theta | lambda2,Y,k,gamma,m
#       ~ N(theta_hat_m, lambda2 A_m^{-1})
#
# The collapsed log target for (k,gamma) is, up to common constants:
#   -0.5 log|C_D| -0.5 log|A_m|
#   + lgamma(n/2) -(n/2) log(S_m/2).
full_fixed_stats <- function(y, t, t_min, t_range, m, k, gamma_delta, prior_type = PRIOR_TYPE) {
  n <- length(y)
  X <- make_X(t, t_min, t_range)
  
  if (m <= 0 || m >= n) {
    return(list(valid = FALSE, log_marginal = -Inf))
  }
  
  idx_P <- 1:m
  idx_D <- (m + 1):n
  
  y_P <- y[idx_P]
  y_D <- y[idx_D]
  
  X_P <- X[idx_P, , drop = FALSE]
  X_D <- X[idx_D, , drop = FALSE]
  
  t_D <- t[idx_D]
  
  C <- make_C(t_D, k, gamma_delta)
  C_inv <- safe_solve(C)
  logdet_C <- safe_logdet(C)
  
  A_data <- t(X_P) %*% X_P + t(X_D) %*% C_inv %*% X_D
  
  B_data <- t(X_P) %*% y_P + t(X_D) %*% C_inv %*% y_D
  
  if (prior_type == "jeffreys") {
    A <- A_data
    A_inv <- safe_solve(A)
    
    theta_hat <- as.vector(A_inv %*% B_data)
    names(theta_hat) <- c("g", "h0")
    
    res_P <- y_P - as.vector(X_P %*% theta_hat)
    res_D <- y_D - as.vector(X_D %*% theta_hat)
    
    SSE <- sum(res_P^2) + as.numeric(t(res_D) %*% C_inv %*% res_D)
    
    if (!is.finite(SSE) || SSE <= 0) {
      return(list(valid = FALSE, log_marginal = -Inf))
    }
    
    logdet_A <- safe_logdet(A)
    
    log_marginal <- (-0.5) * logdet_C - 0.5 * logdet_A + lgamma(n / 2) - (n / 2) * log(SSE / 2)
    
    return(list(
      valid = TRUE,
      log_marginal = log_marginal,
      theta_hat = theta_hat,
      A_inv = A_inv,
      SSE = SSE,
      lambda_shape = n / 2,
      lambda_rate = SSE / 2,
      C = C,
      C_inv = C_inv,
      idx_P = idx_P,
      idx_D = idx_D
    ))
  }
  
  # ---------------------------------------------------------------------------
  # OPTIONAL PROPER-PRIOR SENSITIVITY BRANCH
  # ---------------------------------------------------------------------------
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
  
  yWy <- sum(y_P^2) + as.numeric(t(y_D) %*% C_inv %*% y_D)
  
  mu_quad <- as.numeric(t(mu0) %*% V0_inv %*% mu0)
  
  S <- yWy + mu_quad - as.numeric(t(B) %*% A_inv %*% B)
  
  a_n <- a0 + n / 2
  b_n <- b0 + S / 2
  
  log_marginal <- -0.5 * logdet_C - 0.5 * safe_logdet(A) -
    0.5 * safe_logdet(V0) + lgamma(a_n) - lgamma(a0) +
    a0 * log(b0) - a_n * log(b_n)
  
  list(
    valid = TRUE,
    log_marginal = log_marginal,
    theta_hat = theta_hat,
    A_inv = A_inv,
    SSE = S,
    lambda_shape = a_n,
    lambda_rate = b_n,
    C = C,
    C_inv = C_inv,
    idx_P = idx_P,
    idx_D = idx_D
  )
}

###############################################################################
## 4. COLLAPSED TARGET FOR k AND gamma_delta
###############################################################################

log_full_fixed_target <- function(y, t, t_min, t_range, m, k, gamma_delta, prior_type = PRIOR_TYPE) {
  
  if (k <= k_lower || k >= k_upper) return(-Inf)
  if (gamma_delta <= gamma_lower || gamma_delta >= gamma_upper) return(-Inf)
  
  st <- full_fixed_stats(y, t, t_min, t_range, m, k, gamma_delta, prior_type)
  
  if (!isTRUE(st$valid)) return(-Inf)
  
  # Uniform priors on k and gamma_delta are constants inside their support.
  st$log_marginal
}

# Truncated-normal random-walk MH update for k.
update_k_full_fixed <- function(k_current, y, t, t_min, t_range, m, gamma_delta,
                                proposal_sd = 0.20, prior_type = PRIOR_TYPE) {
  k_prop <- rtruncnorm(1, a = k_lower, b = k_upper, mean = k_current, sd = proposal_sd)
  
  log_current <- log_full_fixed_target(y, t, t_min, t_range, m, k_current, gamma_delta, prior_type)
  
  log_prop <- log_full_fixed_target(y, t, t_min, t_range, m, k_prop, gamma_delta, prior_type)
  
  log_q_reverse <- log(dtruncnorm(k_current, a = k_lower, b = k_upper, mean = k_prop, sd = proposal_sd))
  
  log_q_forward <- log(dtruncnorm(k_prop, a = k_lower, b = k_upper, mean = k_current, sd = proposal_sd))
  
  log_ratio <-(log_prop - log_current) + (log_q_reverse - log_q_forward)
  
  if (
    is.finite(log_ratio) &&
    log(runif(1)) < min(0, log_ratio)
  ) {
    list(value = k_prop, accepted = 1)
  } else {
    list(value = k_current, accepted = 0)
  }
}

# Truncated-normal random-walk MH update for gamma_delta.
update_gamma_full_fixed <- function(gamma_current, y, t, t_min, t_range, m, k, proposal_sd = 0.20, prior_type = PRIOR_TYPE) {
  gamma_prop <- rtruncnorm(1, a = gamma_lower, b = gamma_upper, mean = gamma_current, sd = proposal_sd)
  
  log_current <- log_full_fixed_target(y, t, t_min, t_range, m, k, gamma_current, prior_type)
  
  log_prop <- log_full_fixed_target(y, t, t_min, t_range, m, k, gamma_prop, prior_type)
  
  log_q_reverse <- log(dtruncnorm(gamma_current, a = gamma_lower, b = gamma_upper, mean = gamma_prop, sd = proposal_sd))
  
  log_q_forward <- log(dtruncnorm(gamma_prop, a = gamma_lower, b = gamma_upper, mean = gamma_current, sd = proposal_sd))
  
  log_ratio <- (log_prop - log_current) + (log_q_reverse - log_q_forward)
  
  if (
    is.finite(log_ratio) &&
    log(runif(1)) < min(0, log_ratio)
  ) {
    list(value = gamma_prop, accepted = 1L)
  } else {
    list(value = gamma_current, accepted = 0L)
  }
}

###############################################################################
## 5. STAGE-1 FIXED-m FULL SAMPLER
###############################################################################

# One MCMC iteration:
#   1) MH update k using the collapsed target.
#   2) MH update gamma_delta using the collapsed target.
#   3) Exact draw lambda2 | Y,k,gamma,m.
#   4) Exact draw theta | lambda2,Y,k,gamma,m.
#
# delta_D is NOT an MCMC state variable here.
run_full_fixed_chain <- function(
    y,
    t_raw,
    m,
    n_iter = 8000,
    burn = 2000,
    init_k = 0.1,
    init_gamma = 0.3,
    proposal_sd_k = 0.2,
    proposal_sd_gamma = 0.2,
    prior_type = PRIOR_TYPE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (burn >= n_iter) stop("burn must be smaller than n_iter.")
  
  prep <- prepare_time(t_raw)
  t <- prep$t
  
  k <- init_k
  gamma_delta <- init_gamma
  
  chain <- data.frame(
    iter = 1:n_iter,
    k = numeric(n_iter),
    gamma_delta = numeric(n_iter),
    lambda2 = numeric(n_iter),
    g = numeric(n_iter),
    h0 = numeric(n_iter)
  )
  
  accept_k <- 0
  accept_gamma <- 0
  
  for (iter in 1:n_iter) {
    out_k <- update_k_full_fixed(
      k_current = k,
      y = y,
      t = t,
      t_min = prep$t_min,
      t_range = prep$t_range,
      m = m,
      gamma_delta = gamma_delta,
      proposal_sd = proposal_sd_k,
      prior_type = prior_type
    )
    
    k <- out_k$value
    accept_k <- accept_k + out_k$accepted
    
    out_gamma <- update_gamma_full_fixed(
      gamma_current = gamma_delta,
      y = y,
      t = t,
      t_min = prep$t_min,
      t_range = prep$t_range,
      m = m,
      k = k,
      proposal_sd = proposal_sd_gamma,
      prior_type = prior_type
    )
    
    gamma_delta <- out_gamma$value
    accept_gamma <- accept_gamma + out_gamma$accepted
    
    st <- full_fixed_stats(y, t, prep$t_min, prep$t_range, m, k, gamma_delta, prior_type)
    
    if (!isTRUE(st$valid)) {
      stop("Invalid Stage-1 Full state encountered.")
    }
    
    lambda2 <- rinvgamma(1, shape = st$lambda_shape, rate = st$lambda_rate)
    
    theta <- as.vector(rmvnorm(1, mean = st$theta_hat, sigma = lambda2 * st$A_inv))
    names(theta) <- c("g", "h0")
    
    chain$k[iter] <- k
    chain$gamma_delta[iter] <- gamma_delta
    chain$lambda2[iter] <- lambda2
    chain$g[iter] <- theta["g"]
    chain$h0[iter] <- theta["h0"]
  }
  
  keep <- (burn + 1):n_iter
  
  list(
    chain = chain[keep, , drop = FALSE],
    acceptance = c(
      k = accept_k / n_iter,
      gamma_delta = accept_gamma / n_iter
    )
  )
}

###############################################################################
## 6. POSTERIOR RECONSTRUCTION OF delta_D
###############################################################################

sample_delta_conditional <- function(y, t_raw, m, theta, lambda2, k, gamma_delta) {
  prep <- prepare_time(t_raw)
  t <- prep$t
  X <- make_X(t, prep$t_min, prep$t_range)
  n <- length(y)
  
  idx_D <- (m + 1):n
  y_D <- y[idx_D]
  X_D <- X[idx_D, , drop = FALSE]
  t_D <- t[idx_D]
  
  R_gamma <- GP_covariance(t_D, sigma_sq_delta = 1, gamma_delta = gamma_delta)
  
  K_D <- R_gamma / k
  C <- diag(length(idx_D)) + K_D
  C_inv <- safe_solve(C)
  
  residual <- y_D - as.vector(X_D %*% theta)
  
  mean_delta <- as.vector(K_D %*% C_inv %*% residual)
  
  Cov_delta <- lambda2 * (K_D - K_D %*% C_inv %*% K_D)
  
  Cov_delta <- 0.5 * (Cov_delta + t(Cov_delta))
  
  delta_draw <- as.vector(rmvnorm(1, mean = mean_delta, sigma = Cov_delta + 1e-10 * diag(length(idx_D))))
  
  list(
    index = idx_D,
    mean = mean_delta,
    covariance = Cov_delta,
    draw = delta_draw
  )
}

reconstruct_delta_fixed <- function(draws, y, t_raw, m, n_reconstruct = 300, seed = 2026) {
  set.seed(seed)
  
  n <- length(y)
  B <- min(n_reconstruct, nrow(draws))
  
  id <- sample(seq_len(nrow(draws)),B,replace = FALSE)
  
  # The first block remains exactly zero by model construction.
  delta_mat <- matrix(0, nrow = B, ncol = n)
  
  for (b in seq_len(B)) {
    dr <- draws[id[b], ]
    
    rec <- sample_delta_conditional(
      y = y,
      t_raw = t_raw,
      m = m,
      theta = c(
        g = dr$g,
        h0 = dr$h0
      ),
      lambda2 = dr$lambda2,
      k = dr$k,
      gamma_delta = dr$gamma_delta
    )
    
    delta_mat[b, rec$index] <- rec$draw
  }
  
  list(
    draws = delta_mat,
    mean = colMeans(delta_mat),
    q025 = apply(delta_mat, 2, quantile, probs = 0.025),
    q975 = apply(delta_mat, 2, quantile, probs = 0.975)
  )
}

###############################################################################
## 7. SUMMARIES AND CONVERGENCE DIAGNOSTICS
###############################################################################

summarize_parameter <- function(x) {
  c(
    mean = mean(x),
    sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    median = median(x),
    q975 = unname(quantile(x, 0.975))
  )
}

make_parameter_summary <- function(draws) {
  rbind(
    g = summarize_parameter(draws$g),
    h0 = summarize_parameter(draws$h0),
    lambda2 = summarize_parameter(draws$lambda2),
    k = summarize_parameter(draws$k),
    gamma_delta = summarize_parameter(draws$gamma_delta)
  )
}

summarize_full_fixed_chains <- function(chains) {
  vars <- c("k", "gamma_delta", "lambda2", "g", "h0")
  
  mcmc_list <- mcmc.list(
    lapply(
      chains,
      function(x)
        mcmc(
          as.matrix(
            x$chain[, vars, drop = FALSE]
          )
        )
    )
  )
  
  list(
    gelman = gelman.diag(
      mcmc_list,
      multivariate = FALSE
    ),
    ess = effectiveSize(mcmc_list),
    acceptance = do.call(
      rbind,
      lapply(chains, function(x) x$acceptance)
    )
  )
}

pool_full_fixed_draws <- function(chains) {
  do.call(
    rbind,
    lapply(chains, function(x) x$chain)
  )
}
############################################################################
# save
############################################################################
save_stage1_outputs <- function(
    result_dir,
    truth,
    chains,
    draws,
    parameter_summary,
    diagnostics,
    delta_summary
) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(chains, file.path(result_dir, "chains.rds"))
  saveRDS(draws, file.path(result_dir, "pooled_draws.rds"))
  saveRDS(delta_summary, file.path(result_dir, "delta_reconstruction.rds"))
  
  write.csv(
    parameter_summary,
    file.path(result_dir, "parameter_summary.csv"),
    row.names = TRUE
  )
  
  write.csv(
    diagnostics$acceptance,
    file.path(result_dir, "acceptance_rates.csv"),
    row.names = TRUE
  )
  
  write.csv(
    data.frame(parameter = names(diagnostics$ess),
               ESS = as.numeric(diagnostics$ess)),
    file.path(result_dir, "effective_sample_size.csv"),
    row.names = FALSE
  )
  
  gelman_mat <- diagnostics$gelman$psrf
  write.csv(
    data.frame(
      parameter = rownames(gelman_mat),
      Rhat = gelman_mat[, 1],
      upper_CI = gelman_mat[, 2]
    ),
    file.path(result_dir, "rhat.csv"),
    row.names = FALSE
  )
  
  delta_df <- data.frame(
    index = seq_along(delta_summary$mean),
    mean = delta_summary$mean,
    q025 = delta_summary$q025,
    q975 = delta_summary$q975
  )
  write.csv(
    delta_df,
    file.path(result_dir, "delta_summary.csv"),
    row.names = FALSE
  )
  
  capture.output(
    print(truth),
    file = file.path(result_dir, "truth.txt")
  )
}

save_plot_png <- function(plot_fun, filename,
                          width = 1400, height = 900, res = 150) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

plot_stage1_data <- function(dat,
                             main = "Stage 1 simulated data") {
  true_boundary_time <- mean(
    dat$data$t_raw[c(dat$truth$m, dat$truth$m + 1)]
  )
  
  plot(
    dat$data$t_raw,
    dat$data$y,
    pch = 19,
    xlab = "Time (seconds)",
    ylab = "Height",
    main = main
  )
  lines(dat$data$t_raw, dat$data$physical_true, lwd = 2)
  abline(v = true_boundary_time, lty = 2, lwd = 2)
  legend(
    "bottomleft",
    legend = c("Observed data", "True physical model", "True boundary"),
    pch = c(19, NA, NA),
    lty = c(NA, 1, 2),
    bty = "n"
  )
}

plot_stage1_true_delta <- function(dat,
                                   main = "True discrepancy") {
  true_boundary_time <- mean(
    dat$data$t_raw[c(dat$truth$m, dat$truth$m + 1)]
  )
  
  plot(
    dat$data$t_raw,
    dat$data$delta_true,
    type = "l",
    lwd = 2,
    xlab = "Time (seconds)",
    ylab = expression(delta(t)),
    main = main
  )
  abline(h = 0, lty = 3)
  abline(v = true_boundary_time, lty = 2, lwd = 2)
}

plot_stage1_trace <- function(chains,
                              main_prefix = "Stage 1 Full fixed-m") {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(3, 2))
  
  variables <- c("g", "h0", "lambda2", "k", "gamma_delta")
  ylabels <- list("g", expression(h[0]), expression(lambda^2),
                  "k", expression(gamma[delta]))
  
  for (j in seq_along(variables)) {
    v <- variables[j]
    matplot(
      sapply(chains, function(x) x$chain[[v]]),
      type = "l",
      lty = 1,
      xlab = "Post-burn-in iteration",
      ylab = ylabels[[j]],
      main = paste(main_prefix, "-", v)
    )
  }
  
  plot.new()
  legend("center", legend = paste("Chain", seq_along(chains)), lty = 1, bty = "n")
}

plot_stage1_posterior <- function(draws,
                                  truth,
                                  scenario_label = "Stage 1",
                                  show_gp_truth = TRUE) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(3, 2))
  
  params <- c("g", "h0", "lambda2", "k", "gamma_delta")
  truth_values <- list(
    g = unname(truth$theta["g"]),
    h0 = unname(truth$theta["h0"]),
    lambda2 = truth$lambda2,
    k = if (show_gp_truth && !is.null(truth$k)) truth$k else NA_real_,
    gamma_delta = if (show_gp_truth && !is.null(truth$gamma_delta))
      truth$gamma_delta else NA_real_
  )
  
  for (v in params) {
    x <- draws[[v]]
    hist(x, breaks = 35, probability = TRUE,
         main = paste(scenario_label, "-", v), xlab = v)
    if (length(unique(x)) > 1L) lines(density(x), lwd = 2)
    tv <- truth_values[[v]]
    if (is.finite(tv)) abline(v = tv, lty = 2, lwd = 2)
  }
  
  plot.new()
  legend(
    "center",
    legend = c("Posterior density", "Data-generating truth (when defined)"),
    lty = c(1, 2),
    lwd = c(2, 2),
    bty = "n"
  )
}

plot_stage1_delta_reconstruction <- function(dat,
                                             delta_summary,
                                             main = "Posterior reconstruction of discrepancy") {
  true_boundary_time <- mean(
    dat$data$t_raw[c(dat$truth$m, dat$truth$m + 1)]
  )
  
  ylim <- range(c(dat$data$delta_true,
                  delta_summary$q025,
                  delta_summary$q975), finite = TRUE)
  
  plot(
    dat$data$t_raw,
    delta_summary$mean,
    type = "l",
    lwd = 2,
    xlab = "Time (seconds)",
    ylab = expression(delta(t)),
    ylim = ylim,
    main = main
  )
  lines(dat$data$t_raw, delta_summary$q025, lty = 3)
  lines(dat$data$t_raw, delta_summary$q975, lty = 3)
  lines(dat$data$t_raw, dat$data$delta_true, lty = 2, lwd = 2)
  abline(h = 0, lty = 3)
  abline(v = true_boundary_time, lty = 2)
  legend(
    "topleft",
    legend = c("Posterior mean", "95% pointwise credible limits",
               "True discrepancy", "True boundary"),
    lty = c(1, 3, 2, 2),
    lwd = c(2, 1, 2, 1),
    bty = "n"
  )
}

save_stage1_text_report <- function(result_dir,
                                    truth,
                                    parameter_summary,
                                    diagnostics) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  capture.output(
    {
      cat("============================================================\n")
      cat("STAGE 1 NUMERICAL SUMMARY\n")
      cat("============================================================\n\n")
      cat("TRUTH / DGP SETTINGS\n")
      print(truth)
      cat("\nPOSTERIOR PARAMETER SUMMARY\n")
      print(parameter_summary)
      cat("\nGELMAN-RUBIN / R-HAT\n")
      print(diagnostics$gelman)
      cat("\nEFFECTIVE SAMPLE SIZE\n")
      print(diagnostics$ess)
      cat("\nMH ACCEPTANCE RATES\n")
      print(diagnostics$acceptance)
    },
    file = file.path(result_dir, "numerical_report.txt")
  )
}

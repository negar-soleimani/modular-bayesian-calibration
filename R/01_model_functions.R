GP_covariance <- function(t, sigma_sq_delta, psi_delta) {
  Sigma <- outer(t, t, function(ti, tj) {
    sigma_sq_delta * exp(-abs(ti - tj) / psi_delta)
  }
  )
  return(Sigma)
}

balldropg <- function(t, theta) {
  g <- theta[1]
  h0 <- theta[2]
  theta_vec <- rbind(h0, g)
  x_vec <- cbind(1, -0.5 * (t * t_range + t_min)^2)
  h <- x_vec %*% theta_vec
  h[h < 0] <- 0
  return(as.vector(h))
}

make_X <- function(t_local) { 
  cbind(1, -0.5 * (t_local * t_range + t_min)^2)
}

safe_solve <- function(M, jitter = 1e-8) {
  M <- 0.5 * (M + t(M))
  solve(M + jitter * diag(nrow(M)))
}
repair_covariance_if_needed <- function(M, tolerance = 1e-8) {
  M <- 0.5 * (M + t(M))
  
  min_eig <- min(
    eigen(M, symmetric = TRUE, only.values = TRUE)$values
  )
  
  repaired <- FALSE
  
  if (min_eig < -tolerance) {
    M <- M + (-min_eig + 1e-10) * diag(nrow(M))
    repaired <- TRUE
  }
  
  if (min_eig >= -tolerance && min_eig < 1e-12) {
    M <- M + 1e-12 * diag(nrow(M))
  }
  
  return(
    list(
      matrix = M,
      repaired = repaired
    )
  )
}

safe_effective_size <- function(chain_matrix) {
  if (is.null(dim(chain_matrix))) {
    chain_matrix <- matrix(chain_matrix, ncol = 1)
  }
  
  if (nrow(chain_matrix) < 3) {
    values <- rep(NA_real_, ncol(chain_matrix))
    names(values) <- colnames(chain_matrix)
    return(values)
  }
  
  result <- tryCatch(
    coda::effectiveSize(coda::mcmc(chain_matrix)),
    error = function(e) {
      values <- rep(NA_real_, ncol(chain_matrix))
      names(values) <- colnames(chain_matrix)
      values
    }
  )
  
  return(as.numeric(result))
}

simulate_one_dataset <- function(scenario) {
  f_true <- balldropg(t, c(g_true, h0_true))
  delta_true <- rep(0, n)
  xD <- seq(0, 1, length.out = n_D)
  
  if (scenario == "no_discrepancy") {
    delta_true[] <- 0
    
  } else if (scenario == "weak_smooth_positive") {
    delta_D <- 1.2 * xD^1.4
    delta_true[I_D] <- delta_D
    
  } else if (scenario == "moderate_smooth_positive") {
    delta_D <- delta_amplitude * xD^1.4
    delta_true[I_D] <- delta_D
    
  } else if (scenario == "delayed_smooth_positive") {
    start <- 0.25
    x_shift <- pmax((xD - start) / (1 - start), 0)
    delta_D <- delta_amplitude * x_shift^1.5
    delta_true[I_D] <- delta_D
    
  } else if (scenario == "moderate_local_bump") {
    center <- 0.55
    width <- 0.16
    bump <- exp(-0.5 * ((xD - center) / width)^2)
    delta_D <- delta_amplitude * bump
    delta_true[I_D] <- delta_D
    
    # Scenario 6: quadratic discrepancy 2.5*x_D^2, intentionally close in shape to the quadratic physical trajectory and therefore prone to confounding.
  } else if (scenario == "moderate_quadratic_confounding") {
    delta_D <- delta_amplitude * xD^2
    delta_true[I_D] <- delta_D
    
  } else {
    stop(paste0("Unknown scenario: ", scenario))
  }
  y <- f_true + delta_true + rnorm(n, mean = 0, sd = sqrt(sigma_sq_err_true))
  
  return(list(y = y, f_true = f_true, delta_true = delta_true, scenario = scenario))
}

###############################################################################
# MODULE P
###############################################################################

mcmc_module_P <- function(y_P, t_P, n_keep = N_KEEP_PHYS, burn = BURN_PHYS, init = init_phys, verbose = FALSE) {
  n_P_local <- length(y_P)
  total_iter <- burn + n_keep
  X_P <- make_X(t_P)
  d_free <- 2
  
  theta <- init
  names(theta) <- c("g", "h0", "sigma_sq_err")
  
  chain <- matrix(NA_real_, nrow = total_iter, ncol = 3)
  colnames(chain) <- c("g", "h0", "sigma_sq_err")
  
  for (iter in seq_len(total_iter)) {
    sigma_sq_err <- as.numeric(theta["sigma_sq_err"])
    
    A <- t(X_P) %*% X_P
    B <- t(X_P) %*% matrix(y_P, ncol = 1)
    A_inv <- safe_solve(A)
    
    mu_theta <- A_inv %*% B
    Sigma_theta <- sigma_sq_err * A_inv
    
    theta_sample <- mvtnorm::rmvnorm(1, mean = as.vector(mu_theta), sigma = Sigma_theta)
    
    h0 <- as.numeric(theta_sample[1])
    g <- as.numeric(theta_sample[2])
    
    f_P <- balldropg(t_P, c(g, h0))
    rss_P <- sum((y_P - f_P)^2)
    
    shape_lambda <- n_P_local / 2 + d_free / 2
    rate_lambda <- max(0.5 * rss_P, 1e-14)
    
    sigma_sq_err <- invgamma::rinvgamma(1, shape = shape_lambda, rate = rate_lambda)
    
    theta <- c(g = g, h0 = h0, sigma_sq_err = sigma_sq_err)
    
    if (any(!is.finite(theta))) {
      stop(paste0("Non-finite Module-P draw at iteration ", iter))
    }
    
    chain[iter, ] <- theta
    
    if (verbose && iter %% 5000 == 0) {
      cat(
        "Module P:", iter, "/", total_iter,
        "| g =", round(g, 4),
        "| h0 =", round(h0, 4),
        "| lambda2 =", round(sigma_sq_err, 6),
        "\n"
      )
    }
  }
  
  keep <- seq.int(from = burn + 1, to = total_iter)
  return(chain[keep, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

update_psi_delta <- function(delta, t_local, sigma_sq_err, k, psi_delta) {
  n_local <- length(delta)
  
  psi_prop <- truncnorm::rtruncnorm(1, a = PSI_LOWER, b = PSI_UPPER, mean = psi_delta, sd = PSI_PROPOSAL_SD)
  Sigma_current <- GP_covariance(t = t_local, sigma_sq_delta = sigma_sq_err / k, psi_delta = psi_delta)
  Sigma_prop <- GP_covariance(t = t_local, sigma_sq_delta = sigma_sq_err / k, psi_delta = psi_prop)
  log_like_current <- tryCatch(mvtnorm::dmvnorm(delta, mean = rep(0, n_local), sigma = Sigma_current, log = TRUE), error = function(e) -Inf)
  log_like_prop <- tryCatch(mvtnorm::dmvnorm(delta, mean = rep(0, n_local), sigma = Sigma_prop, log = TRUE), error = function(e) -Inf)
  log_prior_current <- dbeta(psi_delta, shape1 = 7, shape2 = 13, log = TRUE)
  log_prior_prop <- dbeta(psi_prop, shape1 = 7, shape2 = 13, log = TRUE)
  log_q_current_given_prop <- log(truncnorm::dtruncnorm(psi_delta, a = PSI_LOWER, b = PSI_UPPER, mean = psi_prop, sd = PSI_PROPOSAL_SD))
  log_q_prop_given_current <- log(truncnorm::dtruncnorm(psi_prop, a = PSI_LOWER, b = PSI_UPPER, mean = psi_delta, sd = PSI_PROPOSAL_SD))
  log_ratio <- (log_like_prop + log_prior_prop) - (log_like_current + log_prior_current) +
    (log_q_current_given_prop - log_q_prop_given_current)
  
  accepted <- FALSE
  
  if (is.finite(log_ratio) && log(runif(1)) < log_ratio) {
    psi_delta <- psi_prop
    accepted <- TRUE
  }
  
  return(
    list(
      psi_delta = psi_delta,
      accepted = accepted
    )
  )
}

# -----------------------------------------------------------------------------

update_k <- function(delta, t_local, sigma_sq_err, psi_delta) {
  
  n_local <- length(delta)
  
  R <- GP_covariance(t = t_local, sigma_sq_delta = 1, psi_delta = psi_delta)
  
  R_inv <- safe_solve(R)
  
  quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
  
  alpha_k <- (n_local / 2) + 1
  beta_k <- max(quad_form_delta / (2 * sigma_sq_err), 1e-14)
  F0 <- pgamma(K_LOWER, shape = alpha_k, rate = beta_k)
  F1 <- pgamma(K_UPPER, shape = alpha_k, rate = beta_k)
  u <- runif(1, F0, F1)
  k <- qgamma(u, shape = alpha_k, rate = beta_k)
  
  return(
    list(
      k = k,
      quad_form_delta = quad_form_delta
    )
  )
}

# -----------------------------------------------------------------------------

sample_delta_mixture <- function(residual, zeta, t_local, sigma_sq_err, psi_delta, k) {
  
  n_local <- length(residual)
  sigma_sq_delta <- sigma_sq_err / k
  Sigma_delta <- GP_covariance(t = t_local, sigma_sq_delta = sigma_sq_delta, psi_delta = psi_delta)
  
  zeta_2_indices <- which(zeta == 2)
  
  if (length(zeta_2_indices) > 0) {
    
    Sigma_delta_ymym <- sigma_sq_err * diag(length(zeta_2_indices)) + Sigma_delta[zeta_2_indices, zeta_2_indices, drop = FALSE]
    Sigma_delta_ym <- Sigma_delta[ , zeta_2_indices, drop = FALSE]
    Sigma_inv <- safe_solve(Sigma_delta_ymym)
    mu_delta_hat <- Sigma_delta_ym %*% Sigma_inv %*% residual[zeta_2_indices]
    Sigma_delta_hat <- Sigma_delta - Sigma_delta_ym %*% Sigma_inv %*% t(Sigma_delta_ym)
    Sigma_delta_hat <- repair_covariance_if_needed(Sigma_delta_hat)$matrix
    delta <- as.vector(mvtnorm::rmvnorm(1, mean = as.vector(mu_delta_hat), sigma = Sigma_delta_hat))
    
  } else {
    
    delta <- as.vector(mvtnorm::rmvnorm(1, mean = rep(0, n_local), sigma = Sigma_delta)
    )
  }
  
  return(delta)
}
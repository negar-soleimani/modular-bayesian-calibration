###############################################################################
# STAGE 3 Main
# FULL BAYESIAN CALIBRATION WITH AN UNKNOWN DISCRETE CHANGE POINT m
###############################################################################

library(mvtnorm)
library(invgamma)
library(truncnorm)
library(coda)

###############################################################################
## 1. HELPERS FunCTIONS
###############################################################################

d <- 2

k_lower <- 0.015
k_upper <- 1

gamma_lower <- 0.1
gamma_upper <- 1

DEFAULT_M_MIN <- d + 1
DEFAULT_MIN_D <- 3


make_m_prior <- function(n, m_min = DEFAULT_M_MIN, min_D = DEFAULT_MIN_D) {
  m_min <- as.integer(m_min)
  m_max <- as.integer(n - min_D)
  
  if (m_min < 1 || m_min > m_max || m_max >= n) {
    stop("Invalid admissible support for m.")
  }
  values <- seq.int(m_min, m_max)
  prob <- rep(1 / length(values), length(values))
  names(prob) <- as.character(values)
  
  list(
    values = values,
    prob = prob,
    m_min = m_min,
    m_max = m_max,
    min_D = min_D
  )
}


safe_solve <- function(M, jitter = 1e-10) {
  M <- 0.5 * (M + t(M))
  solve(M + jitter * diag(nrow(M)))
}


safe_logdet <- function(M, jitter = 1e-10) {
  M <- 0.5 * (M + t(M))
  M <- M + jitter * diag(nrow(M))
  
  out <- determinant(M, logarithm = TRUE)
  
  if (out$sign <= 0) {
    return(-Inf)
  }
  as.numeric(out$modulus)
}


prepare_time <- function(t_raw) {
  t_min <- min(t_raw)
  t_range <- max(t_raw) - t_min
  
  if (!is.finite(t_range) || t_range <= 0) {
    stop("The observation times must not all be identical.")
  }
  
  list(
    t_raw = as.numeric(t_raw),
    u = (t_raw - t_min) / t_range,
    t_min = t_min,
    t_range = t_range
  )
}


make_X <- function(t_raw) {
  cbind(g = -0.5 * t_raw^2, h0 = 1)
}


ball_drop_mean <- function(t_raw, theta) {
  X <- make_X(t_raw)
  h <- as.vector(X %*% as.numeric(theta))
  h[h < 0] <- 0
  h
}


GP_correlation <- function(u, gamma_delta) {
  outer(u, u, function(ui, uj) {
      exp(-abs(ui - uj) / gamma_delta)
    }
  )
}

GP_covariance <- function(u, sigma_sq_delta, gamma_delta) {
  sigma_sq_delta * GP_correlation(u = u, gamma_delta = gamma_delta)
}


make_C_D <- function(u_D, k, gamma_delta) {
  if (k <= k_lower || k >= k_upper) {
    stop("k is outside its prior support.")
  }
  
  if (gamma_delta <= gamma_lower || gamma_delta >= gamma_upper) {
    stop("gamma_delta is outside its prior support.")
  }
  
  R_gamma <- GP_correlation(u = u_D, gamma_delta = gamma_delta)
  diag(length(u_D)) + R_gamma / k
}

log_q_tnorm_rw <- function(y, x, sd, lower, upper) {
  if (y <= lower || y >= upper) {
    return(-Inf)
  }
  
  z <- pnorm((upper - x) / sd) - pnorm((lower - x) / sd)
  
  if (!is.finite(z) || z <= 0) {
    return(-Inf)
  }
  
  dnorm(y, mean = x, sd = sd, log = TRUE) - log(z)
}


log_downstream_target <- function(Y_D, X_D, u_D, theta, lambda2, k, gamma_delta) {
  if (!is.finite(lambda2) || lambda2 <= 0) {
    return(-Inf)
  }
  
  if (k <= k_lower || k >= k_upper) {
    return(-Inf)
  }
  
  if (gamma_delta <= gamma_lower || gamma_delta >= gamma_upper) {
    return(-Inf)
  }
  
  C_D <- make_C_D(u_D = u_D, k = k, gamma_delta = gamma_delta)
  C_D_inv <- safe_solve(C_D)
  logdet_C_D <- safe_logdet(C_D)
  
  if (!is.finite(logdet_C_D)) {
    return(-Inf)
  }
  
  r_D <- Y_D - as.vector(X_D %*% as.numeric(theta))
  quadratic_D <- as.numeric(crossprod(r_D, C_D_inv %*% r_D)) 
  if (!is.finite(quadratic_D) || quadratic_D < 0) {
    return(-Inf)
  }
  
  log_target <- -0.5 * logdet_C_D - quadratic_D / (2 * lambda2)
  
  as.numeric(log_target)
}

###############################################################################
## 2. METROPOLIS-WITHIN-GIBBS CHAIN
###############################################################################

run_full_unknown_m_chain <- function(
    y,
    t_raw,
    n_iter = 8000,
    burn = 2000,
    m_min = DEFAULT_M_MIN,
    min_D = DEFAULT_MIN_D,
    init_m = NULL,
    init_theta = NULL,
    init_lambda2 = NULL,
    init_k = 0.1,
    init_gamma = 0.3,
    proposal_sd_k = 0.2,
    proposal_sd_gamma = 0.2,
    seed = NULL
) {
  
  ###########################################################################
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(y)
  
  if (length(t_raw) != n) {
    stop("y and t_raw must have the same length.")
  }
  
  if (burn >= n_iter) {
    stop("burn must be smaller than n_iter.")
  }
  
  if (init_k <= k_lower || init_k >= k_upper) {
    stop("Invalid init_k.")
  }
  
  if (init_gamma <= gamma_lower || init_gamma >= gamma_upper) {
    stop("Invalid init_gamma.")
  }
  
  prep <- prepare_time(t_raw)
  X <- make_X(t_raw)
  
  m_prior <- make_m_prior(n = n, m_min = m_min, min_D = min_D)
  
  ###########################################################################
  ## 2.2 INITIAL VALUES
  ###########################################################################
  
  if (is.null(init_m)) {
    init_m <-
      m_prior$values[
        ceiling(
          length(m_prior$values) / 2
        )
      ]
  }
  
  m <- as.integer(init_m)
  
  if (!(m %in% m_prior$values)) {
    stop("init_m is outside the admissible support.")
  }
  
  k <- init_k
  gamma_delta <- init_gamma
  
  # Build the initial partition only to initialize theta and lambda2.
  idx_P <- seq_len(m)
  idx_D <- seq.int(m + 1, n)
  
  Y_P <- y[idx_P]
  Y_D <- y[idx_D]
  
  X_P <- X[idx_P, , drop = FALSE]
  X_D <- X[idx_D, , drop = FALSE]
  
  u_D <- prep$u[idx_D]
  C_D <- make_C_D(u_D = u_D, k = k, gamma_delta = gamma_delta)
  C_D_inv <- safe_solve(C_D)
  A_m <- crossprod(X_P) + crossprod(X_D, C_D_inv %*% X_D)
  B_m <- crossprod(X_P, Y_P) + crossprod(X_D, C_D_inv %*% Y_D)
  A_m_inv <- safe_solve(A_m)
  theta_hat_m <- as.vector(A_m_inv %*% B_m)
  names(theta_hat_m) <- c("g", "h0")
  
  if (is.null(init_theta)) {
    theta <- theta_hat_m
  } else {
    theta <- as.numeric(init_theta)
    
    if (
      length(theta) != d ||
      any(!is.finite(theta))
    ) {
      stop("Invalid init_theta.")
    }
    
    names(theta) <- c("g", "h0")
  }
  
  if (is.null(init_lambda2)) {
    
    r_P_init <- Y_P - as.vector(X_P %*% theta)
    r_D_init <- Y_D - as.vector(X_D %*% theta)
    Q_m_init <- as.numeric(crossprod(r_P_init) + crossprod(r_D_init, C_D_inv %*% r_D_init))
    
    if (!is.finite(Q_m_init) || Q_m_init <= 0) {
      stop("Initial Q_m(theta) must be positive and finite.")
    }
    
    lambda2 <- rinvgamma(n = 1, shape = (n + d) / 2, rate = Q_m_init / 2)
    
  } else {
    
    lambda2 <- as.numeric(init_lambda2)
    
    if (
      length(lambda2) != 1 || !is.finite(lambda2) || lambda2 <= 0
    ) {
      stop("Invalid init_lambda2.")
    }
  }
  
  
  ###########################################################################
  raw <- data.frame(
    iter = seq_len(n_iter),
    m = integer(n_iter),
    g = numeric(n_iter),
    h0 = numeric(n_iter),
    lambda2 = numeric(n_iter),
    k = numeric(n_iter),
    gamma_delta = numeric(n_iter),
    Q_m_for_lambda2 = numeric(n_iter)
  )
  
  accept_k <- 0
  accept_gamma <- 0
  
  m_prob_sum <- rep(0, length(m_prior$values))
  names(m_prob_sum) <- as.character(m_prior$values)
  m_prob_count <- 0
  
  ###########################################################################
  ## 2.3 MAIN MCMC LOOP
  ###########################################################################
  
  for (s in seq_len(n_iter)) {
    
    #########################################################################
    ## STEP 1: EXACT DISCRETE GIBBS UPDATE OF m
    ##
    ## For every candidate j:
    ## log w_j = log pi(j) - 0.5 log|C_D(j;k,gamma_delta)| - Q_j(theta)/(2 lambda2)
    #########################################################################
    
    log_weights <- rep(-Inf, length(m_prior$values))
  
    for (j_index in seq_along(m_prior$values)) {
      
      candidate_m <- m_prior$values[j_index]
      candidate_idx_P <- seq_len(candidate_m)
      candidate_idx_D <- seq.int(candidate_m + 1, n)
      candidate_Y_P <- y[candidate_idx_P]
      candidate_Y_D <- y[candidate_idx_D]
      candidate_X_P <- X[candidate_idx_P, , drop = FALSE]
      candidate_X_D <- X[candidate_idx_D, , drop = FALSE]
      candidate_u_D <- prep$u[candidate_idx_D]
      candidate_C_D <- make_C_D(u_D = candidate_u_D, k = k, gamma_delta = gamma_delta)
      candidate_C_D_inv <- safe_solve(candidate_C_D)
      candidate_logdet_C_D <- safe_logdet(candidate_C_D)
      
      if (!is.finite(candidate_logdet_C_D)) {
        next
      }

      candidate_r_P <- candidate_Y_P - as.vector(candidate_X_P %*% theta)
      candidate_r_D <- candidate_Y_D - as.vector(candidate_X_D %*% theta)
      candidate_Q <- as.numeric(crossprod(candidate_r_P) + crossprod(candidate_r_D, candidate_C_D_inv %*% candidate_r_D))
      
      log_weights[j_index] <- log(m_prior$prob[j_index]) - 0.5 * candidate_logdet_C_D - candidate_Q / (2 * lambda2)
    }
    
    if (!any(is.finite(log_weights))) {
      stop("All candidate m log-weights are invalid.")
    }
    
    max_log_weight <- max(log_weights)
    weights <- exp(log_weights - max_log_weight)
    m_probabilities <- weights / sum(weights)
    names(m_probabilities) <- as.character(m_prior$values)
    m <- as.integer(sample(m_prior$values, size = 1, prob = m_probabilities))
    
    # Accumulate the complete-conditional probabilities only for retained
    # iterations. (Rao-Blackwellized PMF)
    if (s > burn) {m_prob_sum <- m_prob_sum + m_probabilities
      m_prob_count <- m_prob_count + 1
    }
    
    #########################################################################
    ## STEP 2: GIBBS UPDATE OF theta AT THE NEW m
    #########################################################################
    
    idx_P <- seq_len(m)
    idx_D <- seq.int(m + 1, n)
    Y_P <- y[idx_P]
    Y_D <- y[idx_D]
    X_P <- X[idx_P, , drop = FALSE]
    X_D <- X[idx_D, , drop = FALSE]
    u_D <- prep$u[idx_D]
    C_D <- make_C_D(u_D = u_D, k = k, gamma_delta = gamma_delta)
    C_D_inv <- safe_solve(C_D)
    logdet_C_D <- safe_logdet(C_D)
    
    if (!is.finite(logdet_C_D)) {
      stop("C_D is not numerically positive definite.")
    }
    
    A_m <- crossprod(X_P) + crossprod(X_D, C_D_inv %*% X_D)
    B_m <- crossprod(X_P, Y_P) + crossprod(X_D, C_D_inv %*% Y_D)
    A_m_inv <- safe_solve(A_m)
    theta_hat_m <- as.vector(A_m_inv %*% B_m)
    names(theta_hat_m) <- c("g", "h0")
    theta <- as.vector(rmvnorm(n = 1, mean = theta_hat_m, sigma = lambda2 * A_m_inv))
    names(theta) <- c("g", "h0")
    
    #########################################################################
    ## STEP 3: GIBBS UPDATE OF lambda2 AT THE NEW m AND NEW theta
    ##
    ## Q_m(theta) = r_P' r_P + r_D' C_D^{-1} r_D
    ##
    ## lambda2 | m,theta,k,gamma_delta,Y ~ IG((n+d)/2, Q_m(theta)/2)
    #########################################################################
    
    r_P <- Y_P - as.vector(X_P %*% theta)
    r_D <- Y_D - as.vector(X_D %*% theta)
    Q_m_theta <- as.numeric(crossprod(r_P) + crossprod(r_D, C_D_inv %*% r_D))
    
    if (!is.finite(Q_m_theta) || Q_m_theta <= 0) {
      stop("Q_m(theta) must be positive and finite.")
    }
    
    lambda2 <- rinvgamma(n = 1, shape = (n + d) / 2, rate = Q_m_theta / 2)
    
    #########################################################################
    ## STEP 4: METROPOLIS-HASTINGS UPDATE OF k
    #########################################################################
    
    k_star <- rtruncnorm(n = 1, a = k_lower, b = k_upper, mean = k, sd = proposal_sd_k)
    log_target_k_current <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta, lambda2 = lambda2, k = k, gamma_delta = gamma_delta)
    log_target_k_star <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta, lambda2 = lambda2, k = k_star, gamma_delta = gamma_delta)
    log_hastings_k <- log_q_tnorm_rw(y = k, x = k_star, sd = proposal_sd_k, lower = k_lower, upper = k_upper) 
                    - log_q_tnorm_rw(y = k_star, x = k, sd = proposal_sd_k, lower = k_lower, upper = k_upper)
    log_alpha_k <- log_target_k_star - log_target_k_current + log_hastings_k
    accepted_k <- log(runif(1)) < min(0, log_alpha_k)
    
    if (accepted_k) {
      k <- k_star
      accept_k <- accept_k + 1
    }
    
    #########################################################################
    ## STEP 5: METROPOLIS-HASTINGS UPDATE OF gamma_delta
    #########################################################################
    
    gamma_star <- rtruncnorm(n = 1, a = gamma_lower, b = gamma_upper, mean = gamma_delta, sd = proposal_sd_gamma)
    # k is the newly updated value from STEP 4.
    log_target_gamma_current <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta, lambda2 = lambda2, k = k, gamma_delta = gamma_delta)
    log_target_gamma_star <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta, lambda2 = lambda2, k = k, gamma_delta = gamma_star)
    log_hastings_gamma <- log_q_tnorm_rw(y = gamma_delta, x = gamma_star, sd = proposal_sd_gamma, lower = gamma_lower, upper = gamma_upper) 
                        - log_q_tnorm_rw(y = gamma_star, x = gamma_delta, sd = proposal_sd_gamma, lower = gamma_lower, upper = gamma_upper)
    log_alpha_gamma <- log_target_gamma_star - log_target_gamma_current + log_hastings_gamma
    accepted_gamma <- log(runif(1)) < min(0, log_alpha_gamma)
    
    if (accepted_gamma) {
      gamma_delta <- gamma_star
      accept_gamma <- accept_gamma + 1
    }
    
    #########################################################################
    
    raw$m[s] <- m
    raw[s, c("g", "h0")] <- theta
    raw$lambda2[s] <- lambda2
    raw$k[s] <- k
    raw$gamma_delta[s] <- gamma_delta
    raw$Q_m_for_lambda2[s] <- Q_m_theta
  }
  
  
  ###########################################################################
  ## 2.5 REMOVE BURN-IN AND BUILD POSTERIOR m SUMMARY
  ###########################################################################
  
  keep <- seq.int(burn + 1, n_iter)
  chain <- raw[keep, , drop = FALSE]

  if (m_prob_count > 0) {
    m_probability_mean <- m_prob_sum / m_prob_count
  } else {
    m_probability_mean <- rep(NA_real_, length(m_prob_sum))
  }
  
  names(m_probability_mean) <- names(m_prob_sum)
  
  ###########################################################################

  list(
    chain = chain,
    raw_chain = raw,
    m_prior = m_prior,
    m_probability_mean = m_probability_mean,
    acceptance = c( k = accept_k / n_iter, gamma_delta = accept_gamma / n_iter),
    scan_order = c("m_Gibbs", "theta_Gibbs", "lambda2_Gibbs", "k_MH", "gamma_MH"),
    settings = list(
      n_iter = n_iter,
      burn = burn,
      m_min = m_min,
      min_D = min_D,
      proposal_sd_k = proposal_sd_k,
      proposal_sd_gamma = proposal_sd_gamma,
      seed = seed
    )
  )
}


###############################################################################
## 3. SUMMARIES 
###############################################################################

parameter_summary <- function(draws, parameters = c("g", "h0", "lambda2", "k", "gamma_delta")) {
  out <- t(
    vapply(
      parameters,
      function(p) {
        x <- draws[[p]]
        c(
          mean = mean(x),
          sd = sd(x),
          q025 = unname(quantile(x, 0.025)),
          median = unname(quantile(x, 0.5)),
          q975 = unname(quantile(x, 0.975))
        )
      },
      numeric(5)
    )
  )
  
  as.data.frame(out)
}

combine_chains <- function(chains) {
  do.call(rbind, lapply(seq_along(chains),
      function(j) {
        x <- chains[[j]]$chain
        x$chain_id <- j
        x
      }
    )
  )
}


continuous_diagnostics <- function(chains, parameters = c("g","h0","lambda2","k","gamma_delta")) {
  mlist <- mcmc.list(lapply(chains,
      function(z) {
        mcmc(z$chain[, parameters])
      }
    )
  )
  
  rhat <- gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  ess <- effectiveSize(mlist)
  data.frame(
    parameter = parameters,
    Rhat = as.numeric(rhat[parameters]),
    ESS = as.numeric(ess[parameters]),
    row.names = NULL
  )
}


acceptance_summary <- function(chains) {
  data.frame(chain = seq_along(chains),
    k = vapply(chains, function(z) {
        z$acceptance["k"]
      },
      numeric(1)
    ),
    
    gamma_delta = vapply(chains,
      function(z) {
        z$acceptance["gamma_delta"]
      },
      numeric(1)
    )
  )
}


m_empirical_pmf <- function(draws, m_values) {
  tab <- table(factor(draws$m, levels = m_values))
  probability <- as.numeric(tab) / sum(tab)
  data.frame(m = m_values, probability = probability)
}

m_rb_pmf <- function(chains) {
  probs <- do.call(
    rbind,
    lapply(
      chains,
      function(z) {
        z$m_probability_mean
      }
    )
  )
  
  data.frame(
    m = as.integer(
      colnames(probs)
    ),
    probability = colMeans(
      probs,
      na.rm = TRUE
    )
  )
}


m_chain_summary <- function(chains) {
  do.call(
    rbind,
    lapply(
      seq_along(chains),
      function(j) {
        x <- chains[[j]]$chain$m
        
        data.frame(
          chain = j,
          mean_m = mean(x),
          median_m = median(x),
          mode_m = as.integer(
            names(
              which.max(
                table(x)
              )
            )
          ),
          change_rate =
            mean(
              diff(x) != 0
            )
        )
      }
    )
  )
}


###############################################################################
## 4. PLOT HELPERS
###############################################################################

plot_simulated_data <- function(dat, main = "Simulated data") {
  plot(dat$data$t_raw, dat$data$y, pch = 19, xlab = "Time (seconds)", ylab = "Height", main = main)
  lines(dat$data$t_raw, dat$data$physical_true, lwd = 2)
  abline(v = dat$data$t_raw[dat$truth$m], lty = 2, lwd = 2)
  legend("bottomleft", legend = c("Observed data", "True physical model", "True change point"), pch = c(19, NA, NA), lty = c(NA, 1, 2), bty = "n")
}

plot_true_discrepancy <- function(dat, main = "DGP true discrepancy") {
  plot(dat$data$t_raw, dat$data$delta_true, type = "l", lwd = 2, xlab = "Time (seconds)", ylab = expression(delta(t)), main = main)
  abline(h = 0, lty = 3)
  abline(v = dat$data$t_raw[dat$truth$m], lty = 2, lwd = 2)
}

plot_m_pmf <- function(pmf, true_m = NULL, main = "Posterior PMF of m") {
  plot(pmf$m, pmf$probability, type = "h", lwd = 3, xlab = "Change point m", ylab = "Posterior probability", main = main)
  points(pmf$m, pmf$probability, pch = 19)
  if (!is.null(true_m)) {
    abline(v = true_m, lty = 2, lwd = 2)
  }
}

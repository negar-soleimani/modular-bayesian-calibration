###############################################################################
# STAGE 1 CORE
# FULL BAYESIAN CALIBRATION WITH A FIXED, KNOWN CHANGE POINT (m)
#
# Model:
#   i <= m : y_i = x_i' theta + epsilon_i
#   i >  m : y_i = x_i' theta + delta_i + epsilon_i
#
#   epsilon_i ~ N(0, lambda2)
#   delta_D  ~ N(0, (lambda2/k) R_gamma)
#
# After integrating delta_D:
#   C_D = I + R_gamma/k
#   Y | theta,lambda2,k,gamma,m ~ N(X theta, lambda2 V_m)
#
# Priors:
#   pi(theta,lambda2) proportional to lambda2^(-(d/2+1))
#   k     ~ Uniform(0.015,1)
#   gamma ~ Uniform(0.1,1)
###############################################################################

library(mvtnorm)
library(invgamma)
library(truncnorm)
library(coda)

###############################################################################
## 1. HELPERS FUNCTIONS
###############################################################################

d <- 2

k_lower <- 0.015
k_upper <- 1

gamma_lower <- 0.1
gamma_upper <- 1


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


# Exponential GP correlation on normalized time u in [0,1].
GP_correlation <- function(u, gamma_delta) {
  outer(u, u, function(ui, uj) {exp(-abs(ui - uj) / gamma_delta)})
}

# Used by the scenario scripts for data generation.
GP_covariance <- function(u, sigma_sq_delta, gamma_delta) {
  sigma_sq_delta * GP_correlation(u = u, gamma_delta = gamma_delta)
}


# C_D = I + R_gamma/k.
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


# Log density of a truncated-Normal random-walk proposal q(y|x).
# It is used for the Hastings correction in both MH updates.
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


# Shared downstream log target:
#
# ell_D(k,gamma_delta)
#   = -0.5 log|C_D|
#     -0.5/lambda2 * r_D' C_D^{-1} r_D
#
# This is exactly the target appearing in the LaTeX document.
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
  
  r_D <- Y_D - as.vector( X_D %*% as.numeric(theta))
  
  q_D <- as.numeric(crossprod(r_D, C_D_inv %*% r_D))
  
  -0.5 * logdet_C_D - q_D / (2 * lambda2)
}


###############################################################################
## 2. METROPOLIS-WITHIN-GIBBS CHAIN
###############################################################################

run_full_fixed_chain <- function(
    y,
    t_raw,
    m,
    n_iter = 8000,
    burn = 2000,
    init_theta = NULL,
    init_lambda2 = NULL,
    init_k = 0.1,
    init_gamma = 0.3,
    proposal_sd_k = 0.2,
    proposal_sd_gamma = 0.2,
    seed = NULL
) {
  
  ###########################################################################
  ## 2.1 BASIC CHECKS AND FIXED QUANTITIES
  ###########################################################################
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(y)
  
  if (length(t_raw) != n) {
    stop("y and t_raw must have the same length.")
  }
  
  if (m < 1 || m >= n) {
    stop("m must satisfy 1 <= m < n.")
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
  
  # -------------------------------------------------------------------------
  # These quantities are computed ONCE because m is fixed and known.
  #
  # Since m does not move in Stage 1:
  #   P = {1,...,m}
  #   D = {m+1,...,n}
  # remain the same during all iterations.
  # -------------------------------------------------------------------------
  
  prep <- prepare_time(t_raw)
  X <- make_X(t_raw)
  
  idx_P <- seq_len(m)
  idx_D <- seq.int(m + 1, n)
  
  Y_P <- y[idx_P]
  Y_D <- y[idx_D]
  
  X_P <- X[idx_P, , drop = FALSE]
  X_D <- X[idx_D, , drop = FALSE]
  
  u_D <- prep$u[idx_D]
  
  ###########################################################################
  ## 2.2 INITIAL VALUES
  ###########################################################################
  
  k <- init_k
  gamma_delta <- init_gamma
  
  C_D <- make_C_D(u_D = u_D, k = k, gamma_delta = gamma_delta)
  C_D_inv <- safe_solve(C_D)
  A_m <- crossprod(X_P) + crossprod(X_D,C_D_inv %*% X_D)
  B_m <- crossprod(X_P, Y_P) + crossprod(X_D, C_D_inv %*% Y_D)
  A_m_inv <- safe_solve(A_m)
  theta_hat_m <- as.vector(A_m_inv %*% B_m)
  
  names(theta_hat_m) <- c("g", "h0")
  
  if (is.null(init_theta)) {
    theta <- theta_hat_m
  } else {
    theta <- as.numeric(init_theta)
    
    if (length(theta) != d || any(!is.finite(theta))) {
      stop("Invalid init_theta.")
    }
    
    names(theta) <- c("g", "h0")
  }
  
  # If lambda2 is not supplied, create only an INITIAL value.
  # This is not part of the recorded MCMC scan.
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
      length(lambda2) != 1 ||
      !is.finite(lambda2) ||
      lambda2 <= 0
    ) {
      stop("Invalid init_lambda2.")
    }
  }
  
  
  ###########################################################################
  ## 2.3 STORAGE
  ###########################################################################
  
  raw <- data.frame(
    iter = seq_len(n_iter),
    g = numeric(n_iter),
    h0 = numeric(n_iter),
    lambda2 = numeric(n_iter),
    k = numeric(n_iter),
    gamma_delta = numeric(n_iter),
    Q_m_for_lambda2 = numeric(n_iter)
  )
  
  accept_k <- 0
  accept_gamma <- 0
  
  
  ###########################################################################
  ## 2.4 MAIN MCMC LOOP
  ###########################################################################
  
  for (s in seq_len(n_iter)) {
    
    #########################################################################
    ## GIBBS UPDATE OF theta
    ##
    ## theta | lambda2,k,gamma_delta,Y,m
    ##   ~ N(theta_hat_m, lambda2 A_m^{-1})
    #########################################################################
    
    # C_D depends on the CURRENT values of k and gamma_delta.
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
    ## GIBBS UPDATE OF lambda2
    ##
    ## Q_m(theta) = r_P' r_P + r_D' C_D^{-1} r_D
    ##
    ## lambda2 | theta,k,gamma_delta,Y,m
    ##   ~ IG((n+d)/2, Q_m(theta)/2)
    #########################################################################
    
    r_P <- Y_P - as.vector(X_P %*% theta)
    r_D <- Y_D - as.vector(X_D %*% theta)
    
    Q_m_theta <- as.numeric(crossprod(r_P) + crossprod(r_D, C_D_inv %*% r_D))
    
    if (!is.finite(Q_m_theta) || Q_m_theta <= 0) {
      stop("Q_m(theta) must be positive and finite.")
    }
    
    shape_lambda2 <- (n + d) / 2
    rate_lambda2 <- Q_m_theta / 2
    
    lambda2 <- rinvgamma(n = 1, shape = shape_lambda2, rate = rate_lambda2)
    
    #########################################################################
    ## METROPOLIS-HASTINGS UPDATE OF k
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
    ## METROPOLIS-HASTINGS UPDATE OF gamma_delta
    #########################################################################
    
    gamma_star <- rtruncnorm(n = 1, a = gamma_lower, b = gamma_upper, mean = gamma_delta, sd = proposal_sd_gamma)
  
    # k is the newly updated value
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
    ## STORE THE CURRENT ITERATION
    #########################################################################
    
    raw[s, c("g", "h0")] <- theta
    raw$lambda2[s] <- lambda2
    raw$k[s] <- k
    raw$gamma_delta[s] <- gamma_delta
    raw$Q_m_for_lambda2[s] <- Q_m_theta
  }
  
  
  ###########################################################################
  ## 2.5 REMOVE BURN-IN AND RETURN RESULTS
  ###########################################################################
  
  keep <- seq.int(burn + 1, n_iter)
  chain <- raw[keep, , drop = FALSE]
  list(
    chain = chain,
    raw_chain = raw,
    acceptance = c(k = accept_k / n_iter, gamma_delta = accept_gamma / n_iter),
    scan_order = c("theta_Gibbs", "lambda2_Gibbs", "k_MH", "gamma_MH"),
    settings = list(
      m = m,
      n_iter = n_iter,
      burn = burn,
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
  out <- t(vapply(parameters, function(p) {
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
  do.call(
    rbind,
    lapply(
      seq_along(chains),
      function(j) {
        x <- chains[[j]]$chain
        x$chain_id <- j
        x
      }
    )
  )
}


continuous_diagnostics <- function(chains, parameters = c("g", "h0", "lambda2", "k", "gamma_delta")) {
  mcmc_list <- mcmc.list(
    lapply(
      chains,
      function(z) {
        mcmc(
          z$chain[, parameters]
        )
      }
    )
  )
  
  rhat <- gelman.diag( mcmc_list, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  ess <- effectiveSize(mcmc_list)
  data.frame(
    parameter = parameters,
    Rhat = as.numeric(rhat[parameters]),
    ESS = as.numeric(ess[parameters]),
    row.names = NULL
  )
}


acceptance_summary <- function(chains) {
  data.frame(
    chain = seq_along(chains),
    
    k = vapply(
      chains,
      function(z) {
        z$acceptance["k"]
      },
      numeric(1)
    ),
    
    gamma_delta = vapply(
      chains,
      function(z) {
        z$acceptance["gamma_delta"]
      },
      numeric(1)
    )
  )
}


###############################################################################
## 4. SIMPLE PLOTTING HELPERS FOR THE SCENARIO SCRIPTS
###############################################################################

plot_simulated_data <- function(dat, main = "Simulated data") {
  
  plot( dat$data$t_raw, dat$data$y, pch = 19, xlab = "Time (seconds)", ylab = "Height", main = main)
  lines(dat$data$t_raw, dat$data$physical_true, lwd = 2)
  abline(v = dat$data$t_raw[dat$truth$m], lty = 2, lwd = 2)
  legend("bottomleft", legend = c("Observed data", "True physical model", "True change point"), pch = c(19, NA, NA), lty = c(NA, 1, 2), bty = "n")

}


plot_true_discrepancy <- function(dat,main = "DGP true discrepancy") {
  
  plot(dat$data$t_raw, dat$data$delta_true, type = "l", lwd = 2, xlab = "Time (seconds)", ylab = expression(delta(t)), main = main)
  abline(h = 0, lty = 3)
  abline(v = dat$data$t_raw[dat$truth$m], lty = 2, lwd = 2)
  
}

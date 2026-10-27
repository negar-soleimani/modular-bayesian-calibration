###############################################################################
# STAGE 2 CORE
# CUT BAYESIAN CALIBRATION WITH A FIXED, KNOWN CHANGE POINT
#
#   pi_cut^m(phi,eta | Y)
#     = pi_P(phi | Y_P,m) * pi_D(eta | Y_D,phi,m),
#
#   phi = (theta,lambda2)
#   eta = (k,gamma_delta).
###############################################################################

library(mvtnorm)
library(invgamma)
library(truncnorm)
library(coda)

###############################################################################
## 1. HELPERS Functions
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
  cbind(
    g = -0.5 * t_raw^2,
    h0 = 1
  )
}


ball_drop_mean <- function(t_raw, theta) {
  X <- make_X(t_raw)
  h <- as.vector(X %*% as.numeric(theta))
  h[h < 0] <- 0
  h
}

GP_correlation <- function(u, gamma_delta) {
  outer(u, u, function(ui, uj) {exp(-abs(ui - uj) / gamma_delta)})
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


# Downstream conditional log target for fixed theta and lambda2:
#
# ell_D(k,gamma_delta)
#   = -0.5 log|C_D|
#     - r_D' C_D^{-1} r_D / (2 lambda2)
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
  q_D <- as.numeric(crossprod(r_D, C_D_inv %*% r_D)) 
  if (!is.finite(quadratic_D) || quadratic_D < 0) {
    return(-Inf)
  }
  
  log_target <-
    -0.5 * logdet_C_D -
    q_D / (2 * lambda2)
  
  as.numeric(log_target)
}

###############################################################################
## CUT SAMPLER
###############################################################################

run_cut_fixed_chain <- function(y, t_raw, m, outer_iter = 1400, outer_burn = 400, outer_thin = 5,
                       inner_iter = 600, inner_burn = 200, init_theta = NULL, init_lambda2 = NULL,
                       init_k = 0.1, init_gamma = 0.3, proposal_sd_k = 0.02, proposal_sd_gamma = 0.15,
                       seed = 2001, downstream_pick = c("random", "last"), verbose = TRUE) {
  
  downstream_pick <- match.arg(downstream_pick)
  set.seed(seed)
  
  ###########################################################################
  
  n <- length(y)
  if (length(t_raw) != n) {
    stop("y and t_raw must have the same length.")
  }
  if (m < d + 1L || m >= n) {
    stop("m must leave a full-rank physical prefix and downstream observations.")
  }
  if (outer_burn >= outer_iter) {
    stop("outer_burn must be smaller than outer_iter.")
  }
  if (outer_thin < 1L) {
    stop("outer_thin must be at least 1.")
  }
  if (inner_burn >= inner_iter) {
    stop("inner_burn must be smaller than inner_iter.")
  }
  if (init_k <= k_lower || init_k >= k_upper) {
    stop("Invalid init_k.")
  }
  if (init_gamma <= gamma_lower || init_gamma >= gamma_upper) {
    stop("Invalid init_gamma.")
  }
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
  ## MODULE P: FIXED QUANTITIES
  ##
  ## A_P = X_P' X_P
  ## B_P = X_P' Y_P
  ## theta_hat_P = A_P^{-1} B_P
  ###########################################################################
  
  A_P <- crossprod(X_P)
  B_P <- crossprod(X_P, Y_P)
  A_P_inv <- safe_solve(A_P)
  
  theta_hat_P <- as.vector(A_P_inv %*% B_P)
  
  names(theta_hat_P) <- c("g", "h0")
  
  ###########################################################################
  ## MODULE P: INITIAL VALUES
  ###########################################################################
  
  if (is.null(init_theta)) {
    theta <- theta_hat_P
  } else {
    theta <- as.numeric(init_theta)
    
    if (length(theta) != d || any(!is.finite(theta))) {
      stop("Invalid init_theta.")
    }
    
    names(theta) <- c("g", "h0")
  }
  
  if (is.null(init_lambda2)) {
    
    r_P_init <- Y_P - as.vector(X_P %*% theta)
    Q_P_init <- as.numeric(crossprod(r_P_init))
    
    if (!is.finite(Q_P_init) || Q_P_init <= 0) {
      stop("Initial Q_P(theta) must be positive and finite.")
    }
    lambda2 <- rinvgamma(n = 1, shape = (m + d) / 2, rate = Q_P_init / 2)
    
  } else {
    
    lambda2 <- as.numeric(init_lambda2)
    
    if (length(lambda2) != 1 || !is.finite(lambda2) || lambda2 <= 0) {
      stop("Invalid init_lambda2.")
    }
  }
  
  
  ###########################################################################
  ## MODULE-P GIBBS CHAIN
  ###########################################################################
  
  outer_raw <- data.frame(
    iter = seq_len(outer_iter),
    g = numeric(outer_iter),
    h0 = numeric(outer_iter),
    lambda2 = numeric(outer_iter),
    Q_P_for_lambda2 = numeric(outer_iter)
  )
  
  for (s in seq_len(outer_iter)) {
    
    #########################################################################
    ## MODULE P - theta
    ##
    ## theta | lambda2,Y_P,m
    ##   ~ N(theta_hat_P, lambda2 A_P^{-1})
    #########################################################################
    
    theta <- as.vector(rmvnorm(n = 1, mean = theta_hat_P, sigma = lambda2 * A_P_inv))
    names(theta) <- c("g", "h0")
    
    #########################################################################
    ## MODULE P - lambda2
    ##
    ## Q_P(theta) = (Y_P-X_P theta)'(Y_P-X_P theta)
    ##
    ## lambda2 | theta,Y_P,m
    ##   ~ IG((m+d)/2, Q_P(theta)/2)
    #########################################################################
    
    r_P <- Y_P - as.vector(X_P %*% theta)
    Q_P_theta <- as.numeric(crossprod(r_P))
    
    if (!is.finite(Q_P_theta) || Q_P_theta <= 0) {
      stop("Q_P(theta) must be positive and finite.")
    }
    
    lambda2 <- rinvgamma(n = 1, shape = (m + d) / 2, rate = Q_P_theta / 2)
    
    #########################################################################
    ## MODULE-P ITERATION
    #########################################################################
    
    outer_raw[s, c("g", "h0")] <- theta
    outer_raw$lambda2[s] <- lambda2
    outer_raw$Q_P_for_lambda2[s] <- Q_P_theta
  }
  
  outer_keep <- seq.int(outer_burn + 1, outer_iter, by = outer_thin)
  
  upstream <- outer_raw[outer_keep, , drop = FALSE]
  
  if (nrow(upstream) < 1) {
    stop("No retained Module-P draws.")
  }
  
  st_P <- list(
    m = as.integer(m),
    Y_P = Y_P,
    X_P = X_P,
    A_P = A_P,
    B_P = B_P,
    A_P_inv = A_P_inv,
    theta_hat_P = theta_hat_P
  )
  
  module_P <- list(
    chain = upstream,
    raw_chain = outer_raw,
    st_P = st_P,
    scan_order = c(
      "theta_Gibbs",
      "lambda2_Gibbs"
    ),
    settings = list(
      n_iter = outer_iter,
      burn = outer_burn,
      thin = outer_thin,
      seed = seed
    )
  )
  
  
  ###########################################################################
  ## CONDITIONAL MODULE-D SAMPLING
  ###########################################################################
  
  T_keep <- nrow(upstream)
  
  cut_draws <- data.frame(
    draw = seq_len(T_keep),
    outer_iter = upstream$iter,
    m = rep(m, T_keep),
    g = upstream$g,
    h0 = upstream$h0,
    lambda2 = upstream$lambda2,
    Q_P_for_lambda2 = upstream$Q_P_for_lambda2,
    k = numeric(T_keep),
    gamma_delta = numeric(T_keep),
    inner_acc_k = numeric(T_keep),
    inner_acc_gamma = numeric(T_keep)
  )
  
  for (r in seq_len(T_keep)) {
    
    theta_r <- c(g = upstream$g[r], h0 = upstream$h0[r])
    lambda2_r <- upstream$lambda2[r]
    
    # Every retained upstream draw receives its own conditional downstream chain.
    k <- init_k
    gamma_delta <- init_gamma
    
    inner_raw <- data.frame(
      iter = seq_len(inner_iter),
      k = numeric(inner_iter),
      gamma_delta = numeric(inner_iter),
      log_target = numeric(inner_iter)
    )
    
    accept_k <- 0
    accept_gamma <- 0
    
    for (q in seq_len(inner_iter)) {
      
      #######################################################################
      ## MODULE D - k
      #######################################################################
      
      k_star <- rtruncnorm(n = 1, a = k_lower, b = k_upper, mean = k, sd = proposal_sd_k)
      log_target_k_current <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta_r, lambda2 = lambda2_r, k = k, gamma_delta = gamma_delta)
      log_target_k_star <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta_r, lambda2 = lambda2_r, k = k_star, gamma_delta = gamma_delta)
      log_hastings_k <- log_q_tnorm_rw(y = k, x = k_star, sd = proposal_sd_k, lower = k_lower, upper = k_upper)
                      - log_q_tnorm_rw(y = k_star, x = k, sd = proposal_sd_k, lower = k_lower, upper = k_upper)
      log_alpha_k <- log_target_k_star - log_target_k_current + log_hastings_k
      
      accepted_k <- log(runif(1)) < min(0, log_alpha_k)
      
      if (accepted_k) {
        k <- k_star
        accept_k <- accept_k + 1
      }
      
      #######################################################################
      ## MODULE D - gamma_delta
      #######################################################################
      
      gamma_star <- rtruncnorm(n = 1, a = gamma_lower, b = gamma_upper, mean = gamma_delta, sd = proposal_sd_gamma)
      # k is the newly updated value from the previous step.
      log_target_gamma_current <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta_r, lambda2 = lambda2_r, k = k, gamma_delta = gamma_delta)
      log_target_gamma_star <- log_downstream_target(Y_D = Y_D, X_D = X_D, u_D = u_D, theta = theta_r, lambda2 = lambda2_r, k = k, gamma_delta = gamma_star)
      log_hastings_gamma <- log_q_tnorm_rw(y = gamma_delta, x = gamma_star, sd = proposal_sd_gamma, lower = gamma_lower, upper = gamma_upper) 
                          - log_q_tnorm_rw(y = gamma_star, x = gamma_delta, sd = proposal_sd_gamma, lower = gamma_lower, upper = gamma_upper)
      
      log_alpha_gamma <- log_target_gamma_star - log_target_gamma_current + log_hastings_gamma
      accepted_gamma <- log(runif(1)) < min(0, log_alpha_gamma)
      
      if (accepted_gamma) {
        gamma_delta <- gamma_star
        accept_gamma <- accept_gamma + 1
      }
      
      #######################################################################
      ## CONDITIONAL MODULE-D ITERATION
      #######################################################################
      
      inner_raw$k[q] <- k
      inner_raw$gamma_delta[q] <- gamma_delta
      
      inner_raw$log_target[q] <- log_downstream_target(Y_D = Y_D, X_D = X_D,
        u_D = u_D, theta = theta_r,
        lambda2 = lambda2_r,
        k = k,
        gamma_delta = gamma_delta
      )
    }
    
    inner_keep <- seq.int(inner_burn + 1, inner_iter)
    inner_chain <- inner_raw[inner_keep, , drop = FALSE]
    
    pick <- if (downstream_pick == "random") {
      sample.int(
        nrow(inner_chain),
        size = 1
      )
    } else {
      nrow(inner_chain)
    }
    
    cut_draws$k[r] <- inner_chain$k[pick]
    cut_draws$gamma_delta[r] <- inner_chain$gamma_delta[pick]
    
    cut_draws$inner_acc_k[r] <- accept_k / inner_iter
    
    cut_draws$inner_acc_gamma[r] <- accept_gamma / inner_iter
    
    if (
      verbose &&
      r %% max(1, floor(T_keep / 10)) == 0
    ) {
      cat(
        "Stage 2: downstream conditional",
        r,
        "of",
        T_keep,
        "complete\n"
      )
    }
  }
  
  
  ###########################################################################
  ## 2.6 RETURN COMPLETE CUT OUTPUT
  ###########################################################################
  
  list(
    module_P = module_P,
    cut_draws = cut_draws,
    mean_inner_acceptance = c(
      k = mean(cut_draws$inner_acc_k),
      gamma_delta = mean(cut_draws$inner_acc_gamma)
    ),
    
    settings = list(
      m = m,
      outer_iter = outer_iter,
      outer_burn = outer_burn,
      outer_thin = outer_thin,
      inner_iter = inner_iter,
      inner_burn = inner_burn,
      proposal_sd_k = proposal_sd_k,
      proposal_sd_gamma = proposal_sd_gamma,
      downstream_pick = downstream_pick,
      seed = seed
    )
  )
}


###############################################################################
## SUMMARIES 
###############################################################################

parameter_summary <- function(
    draws,
    parameters = c(
      "g",
      "h0",
      "lambda2",
      "k",
      "gamma_delta"
    )
) {
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


combine_cut_chains <- function(chains) {
  do.call(
    rbind,
    lapply(
      seq_along(chains),
      function(j) {
        x <- chains[[j]]$cut_draws
        x$chain_id <- j
        x
      }
    )
  )
}


# Rhat and ESS are evaluated on Module P because Module P is the actual
# outer MCMC chain whose marginal must remain pi_P(phi | Y_P,m).
module_P_diagnostics <- function(
    chains,
    parameters = c(
      "g",
      "h0",
      "lambda2"
    )
) {
  mlist <- mcmc.list(
    lapply(
      chains,
      function(z) {
        mcmc(
          z$module_P$chain[, parameters]
        )
      }
    )
  )
  
  rhat <-
    gelman.diag(
      mlist,
      autoburnin = FALSE,
      multivariate = FALSE
    )$psrf[, 1]
  
  ess <- effectiveSize(mlist)
  
  data.frame(
    parameter = parameters,
    Rhat = as.numeric(rhat[parameters]),
    ESS = as.numeric(ess[parameters]),
    row.names = NULL
  )
}


inner_acceptance_summary <- function(chains) {
  data.frame(
    chain = seq_along(chains),
    
    k = vapply(
      chains,
      function(z) {
        z$mean_inner_acceptance["k"]
      },
      numeric(1)
    ),
    
    gamma_delta = vapply(
      chains,
      function(z) {
        z$mean_inner_acceptance["gamma_delta"]
      },
      numeric(1)
    )
  )
}


###############################################################################
## 4. SIMPLE PLOTTING HELPERS
###############################################################################

plot_simulated_data <- function(dat,main = "Simulated data") {
  plot(dat$data$t_raw, dat$data$y, pch = 19, xlab = "Time (seconds)", ylab = "Height", main = main)
  lines(dat$data$t_raw, dat$data$physical_true,lwd = 2)
  abline(v = dat$data$t_raw[dat$truth$m], lty = 2, lwd = 2)
  legend("bottomleft",legend = c("Observed data", "True physical model", "True change point"),pch = c(19, NA, NA), lty = c(NA, 1, 2), bty = "n")
}

plot_true_discrepancy <- function(dat,main = "DGP true discrepancy") {
  plot(dat$data$t_raw, dat$data$delta_true, type = "l", lwd = 2, xlab = "Time (seconds)", ylab = expression(delta(t)), main = main)
  abline(h = 0, lty = 3)
  abline(v = dat$data$t_raw[dat$truth$m], lty = 2, lwd = 2)
}

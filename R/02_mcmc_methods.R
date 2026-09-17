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

###############################################################################
# FULL MIXTURE
###############################################################################

mcmc_full_mixture <- function(y, t_all, n_keep = N_KEEP_FULL, burn = BURN_FULL, init_phys_local = init_phys, init_mix_local = init_mix, verbose = TRUE) {
  n_local <- length(y)
  total_iter <- burn + n_keep
  d_free <- 2
  X <- make_X(t_all)
  
  g <- as.numeric(init_phys_local["g"])
  h0 <- as.numeric(init_phys_local["h0"])
  sigma_sq_err <- as.numeric(init_phys_local["sigma_sq_err"])
  psi_delta <- as.numeric(init_mix_local["psi_delta"])
  k <- as.numeric(init_mix_local["k"])
  alpha <- as.numeric(init_mix_local["alpha"])
  
  delta <- rep(0, n_local)
  zeta <- ifelse(seq_len(n_local) %% 2 == 0, 2L, 1L)
  
  chain_main <- matrix(NA_real_, nrow = n_keep, ncol = 6)
  
  colnames(chain_main) <- c(
    "g",
    "h0",
    "sigma_sq_err",
    "psi_delta",
    "k",
    "alpha"
  )
  
  delta_eff_keep <- matrix(NA_real_, nrow = n_keep, ncol = n_local)
  colnames(delta_eff_keep) <- paste0("delta_eff_", seq_len(n_local))
  zeta2_sum <- numeric(n_local)
  accept_psi <- 0
  psd_repair_count <- 0
  
  # <- t(X) %*% X
  #A_inv <- safe_solve(A)
  
  for (iter in seq_len(total_iter)) {
    
    # --------------- THETA UPDATE
    
    # THETA UPDATE
    
    idx_1 <- which(zeta == 1)
    idx_2 <- which(zeta == 2)
    
    X_1 <- X[idx_1, , drop = FALSE]
    X_2 <- X[idx_2, , drop = FALSE]
    
    y_1 <- y[idx_1]
    y_2 <- y[idx_2]
    
    delta_2 <- delta[idx_2]
    
    A <- t(X_1) %*% X_1 + t(X_2) %*% X_2
    
    B <- t(X_1) %*% matrix(y_1, ncol = 1) + t(X_2) %*% matrix(y_2 - delta_2, ncol = 1)
    
    A_inv <- safe_solve(A)
    
    mu_theta <- A_inv %*% B
    Sigma_theta <- sigma_sq_err * A_inv
    
    theta_sample <- mvtnorm::rmvnorm(1, mean = as.vector(mu_theta), sigma = Sigma_theta)
    
    h0 <- as.numeric(theta_sample[1])
    g  <- as.numeric(theta_sample[2])
    
    # ------------- DELTA UPDATE
    
    f <- balldropg(t_all, c(g, h0))
    residual <- y - f
    delta <- sample_delta_mixture(residual = residual, zeta = zeta, t_local = t_all, sigma_sq_err = sigma_sq_err, psi_delta = psi_delta, k = k)
    
    # ------------- ZETA UPDATE
    
    log_w1 <- log(alpha) + dnorm(y, mean = f, sd = sqrt(sigma_sq_err), log = TRUE)
    log_w2 <- log(1 - alpha) + dnorm(y, mean = f + delta, sd = sqrt(sigma_sq_err), log = TRUE)
    log_max <- pmax(log_w1, log_w2)
    log_den <- log_max + log(exp(log_w1 - log_max) + exp(log_w2 - log_max))
    prob_zeta2 <- exp(log_w2 - log_den)
    
    zeta <- ifelse(runif(n_local) < prob_zeta2, 2, 1)
    
    # ------------- ALPHA UPDATE
    alpha <- rbeta(1, sum(zeta == 1) + 0.5, sum(zeta == 2) + 0.5)
    
    # ------------ PSI_DELTA UPDATE
    
    psi_res <- update_psi_delta(delta = delta, t_local = t_all, sigma_sq_err = sigma_sq_err, k = k, psi_delta = psi_delta)
    psi_delta <- psi_res$psi_delta
    if (isTRUE(psi_res$accepted)) accept_psi <- accept_psi + 1
    
    # ----------- K UPDATE
    
    k_res <- update_k(delta = delta, t_local = t_all, sigma_sq_err = sigma_sq_err, psi_delta = psi_delta)
    k <- k_res$k
    
    # ----------- SIGMA_SQ_ERR UPDATE
    
    f <- balldropg(t_all, c(g, h0))
    
    rss1 <- sum((y[zeta == 1] - f[zeta == 1])^2)
    rss2 <- sum((y[zeta == 2] - f[zeta == 2] - delta[zeta == 2])^2)
    
    R <- GP_covariance(t = t_all, sigma_sq_delta = 1, psi_delta = psi_delta)
    
    R_inv <- safe_solve(R)
    quad_form_delta <- as.numeric(t(delta) %*% R_inv %*% delta)
    quad_form_delta <- max(quad_form_delta, 0)
    rate_err <- 0.5 * (rss1 + rss2 + k * quad_form_delta)
    rate_err <- max(rate_err, 1e-14)
    shape_err <- n_local + d_free / 2
    sigma_sq_err <- invgamma::rinvgamma(1, shape = shape_err, rate = rate_err)
    
    # -------------------------------------------------------------------------
    
    if (iter > burn) {
      keep_id <- iter - burn
      
      chain_main[keep_id, ] <- c(g = g, h0 = h0, sigma_sq_err = sigma_sq_err, psi_delta = psi_delta, k = k, alpha = alpha)
      
      active_keep <- as.numeric(zeta == 2)
      delta_eff_keep[keep_id, ] <- delta * active_keep
      zeta2_sum <- zeta2_sum + active_keep
    }
    
    if (verbose && iter %% 5000 == 0) {
      cat(
        "Full Mixture:", iter, "/", total_iter,
        "| g =", round(g, 4),
        "| h0 =", round(h0, 4),
        "| lambda2 =", round(sigma_sq_err, 6),
        "| alpha =", round(alpha, 3),
        "\n"
      )
    }
  }
  
  ess_values <- safe_effective_size(chain_main)
  names(ess_values) <- colnames(chain_main)
  
  return(
    list(
      chain = chain_main,
      delta_eff = delta_eff_keep,
      zeta2_prob = zeta2_sum / n_keep,
      accept_rate_psi = accept_psi / total_iter,
      psd_repair_count = psd_repair_count,
      ESS = ess_values
    )
  )
}

###############################################################################
# MODULE D - NO MIXTURE, CONDITIONALLY ON PHYSICAL PARAMETERS
###############################################################################

mcmc_module_D_nomix_given_phys <- function(y_D, t_D, g, h0, sigma_sq_err, n_keep = N_KEEP_D_PER_PHYS, burn = BURN_D_PER_PHYS, init = init_D, verbose = FALSE) {
  n_D_local <- length(y_D)
  total_iter <- burn + n_keep
  
  psi_delta <- as.numeric(init["psi_delta"])
  k <- as.numeric(init["k"])
  delta <- rep(0, n_D_local)
  
  chain_hyper <- matrix(NA_real_, nrow = n_keep, ncol = 2)
  colnames(chain_hyper) <- c("psi_delta", "k")
  delta_keep <- matrix(NA_real_, nrow = n_keep, ncol = n_D_local)
  colnames(delta_keep) <- paste0("delta_", I_D)
  
  accept_psi <- 0
  psd_repair_count <- 0
  
  f_D <- balldropg(t_D, c(g, h0))
  residual_D <- y_D - f_D
  
  for (iter in seq_len(total_iter)) { 
    sigma_sq_delta <- sigma_sq_err / k
    
    Sigma_delta <- GP_covariance( t = t_D, sigma_sq_delta = sigma_sq_delta, psi_delta = psi_delta)
    Sigma_y <- Sigma_delta + sigma_sq_err * diag(n_D_local)
    Sigma_y_inv <- safe_solve(Sigma_y)
    mu_delta <- Sigma_delta %*% Sigma_y_inv %*% residual_D
    V_delta <- Sigma_delta - Sigma_delta %*% Sigma_y_inv %*% Sigma_delta
    psd_result <- repair_covariance_if_needed(V_delta)
    V_delta <- psd_result$matrix
    if (isTRUE(psd_result$repaired)) psd_repair_count <- psd_repair_count + 1
    delta <- as.vector(mvtnorm::rmvnorm(1, mean = as.vector(mu_delta), sigma = V_delta))
    psi_res <- update_psi_delta(delta = delta, t_local = t_D, sigma_sq_err = sigma_sq_err, k = k, psi_delta = psi_delta)
    psi_delta <- psi_res$psi_delta
    if (isTRUE(psi_res$accepted)) accept_psi <- accept_psi + 1
    
    k_res <- update_k(delta = delta, t_local = t_D, sigma_sq_err = sigma_sq_err, psi_delta = psi_delta)
    
    k <- k_res$k
    
    if (iter > burn) {
      keep_id <- iter - burn
      chain_hyper[keep_id, ] <- c(
        psi_delta = psi_delta,
        k = k
      )
      delta_keep[keep_id, ] <- delta
    }
    
    if (verbose && iter %% 500 == 0) {
      cat(
        "  D No-Mixture:", iter, "/", total_iter,
        "| fixed g =", round(g, 4),
        "| psi =", round(psi_delta, 4),
        "| k =", round(k, 4),
        "\n"
      )
    }
  }
  
  ess_values <- safe_effective_size(chain_hyper)
  names(ess_values) <- colnames(chain_hyper)
  
  return(
    list(
      hyper = chain_hyper,
      delta = delta_keep,
      accept_rate_psi = accept_psi / total_iter,
      psd_repair_count = psd_repair_count,
      ESS = ess_values,
      max_physical_change = 0
    )
  )
}

###############################################################################
# CONDITIONAL CUT - NO MIXTURE
###############################################################################

mcmc_cut_nomix_conditional <- function(y_D, t_D, phys_chain, phys_ids, n_keep_D = N_KEEP_D_PER_PHYS, burn_D = BURN_D_PER_PHYS, init = init_D, verbose = TRUE) {
  n_outer_phys <- length(phys_ids)
  
  outer_phys_draws <- phys_chain[phys_ids, c("g", "h0", "sigma_sq_err"), drop = FALSE]
  
  hyper_list <- vector("list", n_outer_phys)
  delta_list <- vector("list", n_outer_phys)
  diagnostics_list <- vector("list", n_outer_phys)
  
  for (outer_id in seq_len(n_outer_phys)) {
    g_fixed <- as.numeric(outer_phys_draws[outer_id, "g"])
    h0_fixed <- as.numeric(outer_phys_draws[outer_id, "h0"])
    sigma_sq_err_fixed <- as.numeric(outer_phys_draws[outer_id, "sigma_sq_err"])
    
    if (verbose) {
      cat(
        "Cut No-Mixture outer", outer_id, "/", n_outer_phys,
        "| g =", round(g_fixed, 4),
        "| h0 =", round(h0_fixed, 4),
        "| lambda2 =", round(sigma_sq_err_fixed, 6),
        "\n"
      )
    }
    
    inner_res <- mcmc_module_D_nomix_given_phys(
      y_D = y_D,
      t_D = t_D,
      g = g_fixed,
      h0 = h0_fixed,
      sigma_sq_err = sigma_sq_err_fixed,
      n_keep = n_keep_D,
      burn = burn_D,
      init = init,
      verbose = FALSE
    )
    
    hyper_list[outer_id] <- list(inner_res$hyper)
    delta_list[outer_id] <- list(inner_res$delta)
    
    diagnostics_list[outer_id] <- list(
      data.frame(
        outer_id = outer_id,
        phys_chain_row = phys_ids[outer_id],
        g = g_fixed,
        h0 = h0_fixed,
        sigma_sq_err = sigma_sq_err_fixed,
        accept_rate_psi = inner_res$accept_rate_psi,
        ESS_psi_delta = as.numeric(inner_res$ESS["psi_delta"]),
        ESS_k = as.numeric(inner_res$ESS["k"]),
        psd_repair_count = inner_res$psd_repair_count,
        max_physical_change = inner_res$max_physical_change,
        row.names = NULL
      )
    )
  }
  
  hyper_all <- do.call(rbind, hyper_list)
  delta_all <- do.call(rbind, delta_list)
  inner_diagnostics <- do.call(rbind, diagnostics_list)
  
  rownames(hyper_all) <- NULL
  rownames(delta_all) <- NULL
  rownames(inner_diagnostics) <- NULL
  
  return(
    list(
      hyper = hyper_all,
      delta = delta_all,
      outer_phys_draws = outer_phys_draws,
      phys_ids = phys_ids,
      inner_diagnostics = inner_diagnostics,
      accept_rate_psi = mean(inner_diagnostics$accept_rate_psi),
      psd_repair_count = sum(inner_diagnostics$psd_repair_count)
    )
  )
}

###############################################################################
# MODULE D - MIXTURE VERSION, CONDITIONALLY ON PHYSICAL PARAMETERS
###############################################################################

mcmc_module_D_mixture_given_phys <- function(
    y_D,
    t_D,
    g,
    h0,
    sigma_sq_err,
    n_keep = N_KEEP_D_PER_PHYS,
    burn = BURN_D_PER_PHYS,
    init = init_mix,
    verbose = FALSE
) {
  n_D_local <- length(y_D)
  total_iter <- burn + n_keep
  
  psi_delta <- as.numeric(init["psi_delta"])
  k <- as.numeric(init["k"])
  alpha <- as.numeric(init["alpha"])
  
  delta <- rep(0, n_D_local)
  zeta <- ifelse(seq_len(n_D_local) %% 2 == 0, 2L, 1L)
  
  chain_hyper <- matrix(
    NA_real_,
    nrow = n_keep,
    ncol = 3
  )
  
  colnames(chain_hyper) <- c("psi_delta", "k", "alpha")
  
  delta_eff_keep <- matrix(
    NA_real_,
    nrow = n_keep,
    ncol = n_D_local
  )
  
  colnames(delta_eff_keep) <- paste0("delta_eff_", I_D)
  
  zeta2_sum <- numeric(n_D_local)
  accept_psi <- 0L
  psd_repair_count <- 0L
  
  f_D <- balldropg(t_D, c(g, h0))
  residual_D <- y_D - f_D
  
  for (iter in seq_len(total_iter)) {
    
    delta <- sample_delta_mixture(
      residual = residual_D,
      zeta = zeta,
      t_local = t_D,
      sigma_sq_err = sigma_sq_err,
      psi_delta = psi_delta,
      k = k
    )
    
    log_w1 <- log(alpha) +
      dnorm(y_D, mean = f_D, sd = sqrt(sigma_sq_err), log = TRUE)
    
    log_w2 <- log(1 - alpha) +
      dnorm(y_D, mean = f_D + delta, sd = sqrt(sigma_sq_err), log = TRUE)
    
    log_max <- pmax(log_w1, log_w2)
    
    log_den <- log_max +
      log(
        exp(log_w1 - log_max) +
          exp(log_w2 - log_max)
      )
    
    prob_zeta2 <- exp(log_w2 - log_den)
    
    zeta <- ifelse(
      runif(n_D_local) < prob_zeta2,
      2L,
      1L
    )
    
    alpha <- rbeta(
      1,
      sum(zeta == 1L) + 0.5,
      sum(zeta == 2L) + 0.5
    )
    
    psi_res <- update_psi_delta(
      delta = delta,
      t_local = t_D,
      sigma_sq_err = sigma_sq_err,
      k = k,
      psi_delta = psi_delta
    )
    
    psi_delta <- psi_res$psi_delta
    if (isTRUE(psi_res$accepted)) accept_psi <- accept_psi + 1L
    
    k_res <- update_k(
      delta = delta,
      t_local = t_D,
      sigma_sq_err = sigma_sq_err,
      psi_delta = psi_delta
    )
    
    k <- k_res$k
    
    if (iter > burn) {
      keep_id <- iter - burn
      active_keep <- as.numeric(zeta == 2L)
      
      chain_hyper[keep_id, ] <- c(
        psi_delta = psi_delta,
        k = k,
        alpha = alpha
      )
      
      delta_eff_keep[keep_id, ] <- delta * active_keep
      zeta2_sum <- zeta2_sum + active_keep
    }
    
    if (verbose && iter %% 500 == 0) {
      cat(
        "  D Mixture:", iter, "/", total_iter,
        "| fixed g =", round(g, 4),
        "| psi =", round(psi_delta, 4),
        "| k =", round(k, 4),
        "| alpha =", round(alpha, 3),
        "\n"
      )
    }
  }
  
  ess_values <- safe_effective_size(chain_hyper)
  names(ess_values) <- colnames(chain_hyper)
  
  return(
    list(
      hyper = chain_hyper,
      delta_eff = delta_eff_keep,
      zeta2_prob = zeta2_sum / n_keep,
      accept_rate_psi = accept_psi / total_iter,
      psd_repair_count = psd_repair_count,
      ESS = ess_values,
      max_physical_change = 0
    )
  )
}

###############################################################################
# CONDITIONAL CUT - MIXTURE VERSION
###############################################################################

mcmc_cut_mixture_conditional <- function(
    y_D,
    t_D,
    phys_chain,
    phys_ids,
    n_keep_D = N_KEEP_D_PER_PHYS,
    burn_D = BURN_D_PER_PHYS,
    init = init_mix,
    verbose = TRUE
) {
  n_outer_phys <- length(phys_ids)
  
  outer_phys_draws <- phys_chain[
    phys_ids,
    c("g", "h0", "sigma_sq_err"),
    drop = FALSE
  ]
  
  hyper_list <- vector("list", n_outer_phys)
  delta_eff_list <- vector("list", n_outer_phys)
  zeta_prob_list <- vector("list", n_outer_phys)
  diagnostics_list <- vector("list", n_outer_phys)
  
  for (outer_id in seq_len(n_outer_phys)) {
    g_fixed <- as.numeric(outer_phys_draws[outer_id, "g"])
    h0_fixed <- as.numeric(outer_phys_draws[outer_id, "h0"])
    sigma_sq_err_fixed <- as.numeric(outer_phys_draws[outer_id, "sigma_sq_err"])
    
    if (verbose) {
      cat(
        "Cut Mixture outer", outer_id, "/", n_outer_phys,
        "| g =", round(g_fixed, 4),
        "| h0 =", round(h0_fixed, 4),
        "| lambda2 =", round(sigma_sq_err_fixed, 6),
        "\n"
      )
    }
    
    inner_res <- mcmc_module_D_mixture_given_phys(
      y_D = y_D,
      t_D = t_D,
      g = g_fixed,
      h0 = h0_fixed,
      sigma_sq_err = sigma_sq_err_fixed,
      n_keep = n_keep_D,
      burn = burn_D,
      init = init,
      verbose = FALSE
    )
    
    hyper_list[outer_id] <- list(inner_res$hyper)
    delta_eff_list[outer_id] <- list(inner_res$delta_eff)
    zeta_prob_list[outer_id] <- list(inner_res$zeta2_prob)
    
    diagnostics_list[outer_id] <- list(
      data.frame(
        outer_id = outer_id,
        phys_chain_row = phys_ids[outer_id],
        g = g_fixed,
        h0 = h0_fixed,
        sigma_sq_err = sigma_sq_err_fixed,
        accept_rate_psi = inner_res$accept_rate_psi,
        ESS_psi_delta = as.numeric(inner_res$ESS["psi_delta"]),
        ESS_k = as.numeric(inner_res$ESS["k"]),
        ESS_alpha = as.numeric(inner_res$ESS["alpha"]),
        psd_repair_count = inner_res$psd_repair_count,
        max_physical_change = inner_res$max_physical_change,
        row.names = NULL
      )
    )
  }
  
  hyper_all <- do.call(rbind, hyper_list)
  delta_eff_all <- do.call(rbind, delta_eff_list)
  inner_diagnostics <- do.call(rbind, diagnostics_list)
  
  zeta_prob_matrix <- do.call(rbind, zeta_prob_list)
  zeta2_prob <- colMeans(zeta_prob_matrix)
  
  rownames(hyper_all) <- NULL
  rownames(delta_eff_all) <- NULL
  rownames(inner_diagnostics) <- NULL
  
  return(
    list(
      hyper = hyper_all,
      delta_eff = delta_eff_all,
      zeta2_prob = zeta2_prob,
      outer_phys_draws = outer_phys_draws,
      phys_ids = phys_ids,
      inner_diagnostics = inner_diagnostics,
      accept_rate_psi = mean(inner_diagnostics$accept_rate_psi),
      psd_repair_count = sum(inner_diagnostics$psd_repair_count)
    )
  )
}

###############################################################################
# PLUG-IN - NO MIXTURE
###############################################################################

mcmc_plugin_nomix <- function(
    y_D,
    t_D,
    phys_chain,
    n_keep = N_KEEP_PLUGIN,
    burn = BURN_PLUGIN,
    init = init_D,
    verbose = FALSE
) {
  g_hat <- mean(phys_chain[, "g"])
  h0_hat <- mean(phys_chain[, "h0"])
  sigma_sq_err_hat <- mean(phys_chain[, "sigma_sq_err"])
  
  plugin_res <- mcmc_module_D_nomix_given_phys(
    y_D = y_D,
    t_D = t_D,
    g = g_hat,
    h0 = h0_hat,
    sigma_sq_err = sigma_sq_err_hat,
    n_keep = n_keep,
    burn = burn,
    init = init,
    verbose = verbose
  )
  
  plugin_res$physical_hat <- c(
    g = g_hat,
    h0 = h0_hat,
    sigma_sq_err = sigma_sq_err_hat
  )
  
  return(plugin_res)
}
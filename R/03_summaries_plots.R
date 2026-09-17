summarize_draw_vector <- function(values) {
  c(
    mean = mean(values),
    sd = sd(values),
    q025 = as.numeric(quantile(values, 0.025, names = FALSE)),
    q50 = as.numeric(quantile(values, 0.50, names = FALSE)),
    q975 = as.numeric(quantile(values, 0.975, names = FALSE))
  )
}

true_value_for_parameter <- function(parameter) {
  if (parameter == "g") return(g_true)
  if (parameter == "h0") return(h0_true)
  if (parameter == "sigma_sq_err") return(sigma_sq_err_true)
  return(NA_real_)
}

physical_summary_from_chain <- function(chain, method) {
  parameters <- c("g", "h0", "sigma_sq_err")
  rows <- vector("list", length(parameters))
  
  for (parameter_id in seq_along(parameters)) {
    parameter <- parameters[parameter_id]
    s <- summarize_draw_vector(chain[, parameter])
    truth <- true_value_for_parameter(parameter)
    
    rows[parameter_id] <- list(
      data.frame(
        method = method,
        parameter = parameter,
        mean = as.numeric(s["mean"]),
        sd = as.numeric(s["sd"]),
        q025 = as.numeric(s["q025"]),
        q50 = as.numeric(s["q50"]),
        q975 = as.numeric(s["q975"]),
        truth = truth,
        bias = as.numeric(s["mean"]) - truth,
        covered = as.integer(
          as.numeric(s["q025"]) <= truth &&
            as.numeric(s["q975"]) >= truth
        ),
        row.names = NULL
      )
    )
  }
  
  return(do.call(rbind, rows))
}

physical_summary_from_point <- function(point_values, method) {
  parameters <- c("g", "h0", "sigma_sq_err")
  rows <- vector("list", length(parameters))
  
  for (parameter_id in seq_along(parameters)) {
    parameter <- parameters[parameter_id]
    value <- as.numeric(point_values[parameter])
    truth <- true_value_for_parameter(parameter)
    
    rows[parameter_id] <- list(
      data.frame(
        method = method,
        parameter = parameter,
        mean = value,
        sd = 0,
        q025 = value,
        q50 = value,
        q975 = value,
        truth = truth,
        bias = value - truth,
        covered = NA_integer_,
        row.names = NULL
      )
    )
  }
  
  return(do.call(rbind, rows))
}

hyper_summary_from_chain <- function(chain, method) {
  parameters <- intersect(c("psi_delta", "k"), colnames(chain))
  
  if (length(parameters) == 0) {
    return(data.frame())
  }
  
  rows <- vector("list", length(parameters))
  
  for (parameter_id in seq_along(parameters)) {
    parameter <- parameters[parameter_id]
    s <- summarize_draw_vector(chain[, parameter])
    
    rows[parameter_id] <- list(
      data.frame(
        method = method,
        parameter = parameter,
        mean = as.numeric(s["mean"]),
        sd = as.numeric(s["sd"]),
        q025 = as.numeric(s["q025"]),
        q50 = as.numeric(s["q50"]),
        q975 = as.numeric(s["q975"]),
        row.names = NULL
      )
    )
  }
  
  return(do.call(rbind, rows))
}

alpha_summary_from_chain <- function(chain, method) {
  if (!("alpha" %in% colnames(chain))) return(data.frame())
  
  s <- summarize_draw_vector(chain[, "alpha"])
  
  data.frame(
    method = method,
    parameter = "alpha",
    mean = as.numeric(s["mean"]),
    sd = as.numeric(s["sd"]),
    q025 = as.numeric(s["q025"]),
    q50 = as.numeric(s["q50"]),
    q975 = as.numeric(s["q975"]),
    row.names = NULL
  )
}

summarize_delta_draws <- function(
    delta_draws,
    method,
    scenario,
    delta_true_D
) {
  means <- colMeans(delta_draws)
  q025 <- apply(delta_draws, 2, quantile, probs = 0.025, names = FALSE)
  q50 <- apply(delta_draws, 2, quantile, probs = 0.50, names = FALSE)
  q975 <- apply(delta_draws, 2, quantile, probs = 0.975, names = FALSE)
  
  data.frame(
    scenario = scenario,
    method = method,
    index = I_D,
    t = t[I_D],
    xD = seq(0, 1, length.out = n_D),
    delta_true = delta_true_D,
    mean = means,
    q025 = q025,
    q50 = q50,
    q975 = q975,
    covered = as.integer(q025 <= delta_true_D & q975 >= delta_true_D),
    row.names = NULL
  )
}

delta_metrics_from_summary <- function(delta_summary) {
  error <- delta_summary$mean - delta_summary$delta_true
  
  data.frame(
    scenario = unique(delta_summary$scenario),
    method = unique(delta_summary$method),
    RMSE_delta = sqrt(mean(error^2)),
    MAE_delta = mean(abs(error)),
    coverage_delta = mean(delta_summary$covered),
    mean_abs_estimated_delta = mean(abs(delta_summary$mean)),
    row.names = NULL
  )
}

plot_physical_parameter <- function(
    physical_summary,
    parameter,
    title,
    y_label
) {
  
  d <- physical_summary[
    physical_summary$parameter == parameter,
    ,
    drop = FALSE
  ]
  
  ggplot(
    d,
    aes(x = method, y = mean, color = method)
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = q025, ymax = q975),
      width = 0.12
    ) +
    geom_hline(
      yintercept = unique(d$truth),
      linetype = "dashed"
    ) +
    labs(
      title = title,
      x = NULL,
      y = y_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      )
    )
}

plot_delta <- function(delta_summary, title) {
  
  ggplot(
    delta_summary,
    aes(x = t)
  ) +
    geom_ribbon(
      aes(
        ymin = q025,
        ymax = q975,
        fill = method
      ),
      alpha = 0.20,
      color = NA
    ) +
    geom_line(
      aes(
        y = mean,
        color = method
      ),
      linewidth = 0.9
    ) +
    geom_line(
      aes(y = delta_true),
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    facet_wrap(
      ~ method,
      ncol = 1
    ) +
    labs(
      title = title,
      x = "t",
      y = expression(delta(t))
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}
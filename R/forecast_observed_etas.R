source("R/etas_likelihood.R")
source("R/loglinear_incompleteness.R")
source("R/fit_loglinear_mle.R")

etas_lambda_at_vectorised <- function(t_eval, event_ts, event_magnitudes,
                                      theta, M_cut, chunk_size = 500L) {
  validate_etas_theta(theta)
  productivity <- etas_productivity(event_magnitudes, theta, M_cut)
  out <- numeric(length(t_eval))
  chunks <- split(seq_along(t_eval), ceiling(seq_along(t_eval) / chunk_size))
  for (idx in chunks) {
    dt <- outer(t_eval[idx], event_ts, "-")
    positive <- dt > 0
    kernel <- ((theta[["p"]] - 1) / theta[["c"]]) *
      (1 + pmax(dt, 0) / theta[["c"]])^(-theta[["p"]])
    kernel[!positive] <- 0
    out[idx] <- theta[["mu"]] + drop(kernel %*% productivity)
  }
  out
}

multi_observed_loglik <- function(theta, catalogue, start, end, M0, b = 1,
                                  model = c("naive", "plugin"),
                                  G = NA, H = NA, M_trigger = 5,
                                  grid_n = 1000L,
                                  fixed_trigger_time = NULL,
                                  fixed_trigger_magnitude = NULL) {
  model <- match.arg(model)
  x <- catalogue[order(catalogue$ts), , drop = FALSE]
  target <- x$ts > start & x$ts <= end
  event_rate <- etas_lambda_at_vectorised(
    x$ts[target], x$ts, x$magnitudes, theta, M0
  )
  use_fixed_trigger <- model == "plugin" &&
    !is.null(fixed_trigger_time) && !is.null(fixed_trigger_magnitude)
  trigger_ts <- if (use_fixed_trigger) {
    fixed_trigger_time
  } else {
    x$ts[x$magnitudes >= M_trigger]
  }
  trigger_magnitudes <- if (use_fixed_trigger) {
    fixed_trigger_magnitude
  } else {
    x$magnitudes[x$magnitudes >= M_trigger]
  }
  psi_event <- if (model == "naive") rep(1, sum(target)) else
    loglinear_detectable_fraction(
      x$ts[target], trigger_ts, trigger_magnitudes, M0, G, H, b
    )
  detection_active <- if (model == "naive" || !length(trigger_ts)) {
    FALSE
  } else {
    recovery_end <- loglinear_recovery_end(
      trigger_ts, trigger_magnitudes, M0, G, H
    )
    any(trigger_ts < end & recovery_end > start)
  }
  integral <- if (model == "naive" || !detection_active) {
    etas_integrated_intensity(start, end, x$ts, x$magnitudes, theta, M0)
  } else {
    q <- make_event_split_quadrature(start, end, x$ts, grid_n)
    grid_rate <- etas_lambda_at_vectorised(
      q$t, x$ts, x$magnitudes, theta, M0
    )
    psi_grid <- loglinear_detectable_fraction(
      q$t, trigger_ts, trigger_magnitudes, M0, G, H, b
    )
    sum(grid_rate * psi_grid * q$weight)
  }
  sum(log(pmax(event_rate * psi_event, 1e-300))) - integral
}

simulate_observed_history_etas <- function(theta, history, start, end, M0,
                                           b = 1, model = c("naive", "plugin"),
                                           G = NA, H = NA, M_trigger = 5,
                                           max_events = 2000L,
                                           fixed_trigger_time = NULL,
                                           fixed_trigger_magnitude = NULL) {
  model <- match.arg(model)
  events <- history[history$ts <= start, c("ts", "magnitudes"), drop = FALSE]
  events <- events[order(events$ts), , drop = FALSE]
  event_ts <- events$ts
  event_magnitudes <- events$magnitudes
  current <- start
  future_ts <- numeric(max_events)
  future_magnitudes <- numeric(max_events)
  n_future <- 0L
  overflow <- FALSE
  beta <- b * log(10)
  eps <- max(.Machine$double.eps^0.5, 1e-10)
  repeat {
    eval_t <- min(current + eps, end)
    bound <- etas_lambda_at(eval_t, event_ts, event_magnitudes, theta, M0)
    if (!is.finite(bound) || bound <= 0) break
    proposal <- current + rexp(1, bound)
    if (proposal > end) break
    latent <- etas_lambda_at(proposal, event_ts, event_magnitudes, theta, M0)
    use_fixed_trigger <- model == "plugin" &&
      !is.null(fixed_trigger_time) && !is.null(fixed_trigger_magnitude)
    trigger_ts <- if (use_fixed_trigger) {
      fixed_trigger_time
    } else {
      event_ts[event_magnitudes >= M_trigger]
    }
    trigger_magnitudes <- if (use_fixed_trigger) {
      fixed_trigger_magnitude
    } else {
      event_magnitudes[event_magnitudes >= M_trigger]
    }
    Mc <- if (model == "naive") M0 else compute_loglinear_mct_at(
      proposal, trigger_ts, trigger_magnitudes, M0, G, H
    )
    psi <- if (model == "naive") 1 else 10^(-b * (Mc - M0))
    if (runif(1) <= min(1, psi * latent / bound)) {
      magnitude <- Mc + rexp(1, beta)
      n_future <- n_future + 1L
      future_ts[n_future] <- proposal
      future_magnitudes[n_future] <- magnitude
      event_ts <- c(event_ts, proposal)
      event_magnitudes <- c(event_magnitudes, magnitude)
      if (n_future >= max_events) {
        overflow <- TRUE
        break
      }
    }
    current <- proposal
  }
  future <- data.frame(
    ts = future_ts[seq_len(n_future)],
    magnitudes = future_magnitudes[seq_len(n_future)]
  )
  attr(future, "overflow") <- overflow
  future
}

log_mean_exp <- function(x) {
  m <- max(x); m + log(mean(exp(x - m)))
}

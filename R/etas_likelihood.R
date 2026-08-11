# Shared ETAS intensity and likelihood helpers.

# Canonical model used everywhere in this project:
#   lambda(t | H_t) = mu + sum K 10^(alpha (M_i - M_cut)) g(t - t_i)
#   g(u) = ((p - 1) / c) (1 + u / c)^(-p), u > 0
# Here K is the expected number of direct offspring from an M_cut event over
# an infinite time horizon. Do not mix these parameters with the unnormalised,
# natural-exponential parameterisation used by ETAS.inlabru::Temporal.ETAS().

validate_etas_theta <- function(theta) {
  required <- c("mu", "K", "alpha", "c", "p")
  if (!all(required %in% names(theta))) {
    stop("theta must contain: ", paste(required, collapse = ", "))
  }
  if (theta[["mu"]] < 0 || theta[["K"]] < 0 || theta[["alpha"]] < 0 ||
      theta[["c"]] <= 0 || theta[["p"]] <= 1) {
    stop("Require mu, K, alpha >= 0, c > 0, and p > 1.")
  }
  invisible(TRUE)
}

etas_productivity <- function(magnitudes, theta, M_cut) {
  theta[["K"]] * 10^(theta[["alpha"]] * (magnitudes - M_cut))
}

etas_omori_density <- function(dt, theta) {
  out <- numeric(length(dt))
  positive <- dt > 0
  out[positive] <- ((theta[["p"]] - 1) / theta[["c"]]) *
    (1 + dt[positive] / theta[["c"]])^(-theta[["p"]])
  out
}

etas_omori_cdf <- function(dt, theta) {
  out <- numeric(length(dt))
  positive <- dt > 0
  out[positive] <- 1 -
    (1 + dt[positive] / theta[["c"]])^(1 - theta[["p"]])
  out
}

etas_lambda_at <- function(t_eval, event_ts, event_magnitudes, theta, M_cut) {
  validate_etas_theta(theta)
  intensity <- rep(theta[["mu"]], length(t_eval))

  for (i in seq_along(t_eval)) {
    past <- which(event_ts < t_eval[i])
    if (length(past) > 0) {
      intensity[i] <- intensity[i] + sum(
        etas_productivity(event_magnitudes[past], theta, M_cut) *
          etas_omori_density(t_eval[i] - event_ts[past], theta)
      )
    }
  }
  intensity
}

etas_integrated_intensity <- function(T1, T2, event_ts, event_magnitudes,
                                      theta, M_cut) {
  validate_etas_theta(theta)
  if (T2 <= T1) stop("T2 must be greater than T1.")

  total <- theta[["mu"]] * (T2 - T1)
  possible_parents <- which(event_ts < T2)
  if (length(possible_parents) == 0) return(total)

  lower <- pmax(T1, event_ts[possible_parents])
  upper_mass <- etas_omori_cdf(T2 - event_ts[possible_parents], theta)
  lower_mass <- etas_omori_cdf(lower - event_ts[possible_parents], theta)
  total + sum(
    etas_productivity(event_magnitudes[possible_parents], theta, M_cut) *
      (upper_mass - lower_mass)
  )
}

etas_to_natural <- function(eta) {
  c(
    mu = exp(eta[1]),
    K = exp(eta[2]),
    alpha = exp(eta[3]),
    c = exp(eta[4]),
    p = 1 + exp(eta[5])
  )
}

etas_to_internal <- function(theta) {
  c(
    log(theta[["mu"]]),
    log(theta[["K"]]),
    log(theta[["alpha"]]),
    log(theta[["c"]]),
    log(theta[["p"]] - 1)
  )
}

baseline_lambda <- function(ts, magnitudes, eta, M_cut) {
  theta <- etas_to_natural(eta)
  etas_lambda_at(ts, ts, magnitudes, theta, M_cut)
}

compute_mct_at <- function(t_eval,
                           event_ts,
                           event_magnitudes,
                           M_cut,
                           M_trigger,
                           gamma,
                           tau,
                           mode = "sum",
                           include_current = TRUE) {
  Mc_eval <- rep(M_cut, length(t_eval))
  trigger_idx <- which(event_magnitudes >= M_trigger)

  if (length(trigger_idx) == 0) {
    return(Mc_eval)
  }

  trigger_ts <- event_ts[trigger_idx]
  trigger_magnitudes <- event_magnitudes[trigger_idx]

  for (i in seq_along(t_eval)) {
    past <- if (include_current) {
      which(trigger_ts <= t_eval[i])
    } else {
      which(trigger_ts < t_eval[i])
    }

    if (length(past) > 0) {
      dt <- t_eval[i] - trigger_ts[past]
      jumps <- gamma * pmax(trigger_magnitudes[past] - M_cut, 0)
      effects <- jumps * exp(-dt / tau)

      if (mode == "max") {
        Mc_eval[i] <- M_cut + max(effects)
      } else {
        Mc_eval[i] <- M_cut + sum(effects)
      }
    }
  }

  Mc_eval
}

tvc_lambda <- function(ts, magnitudes, Mc_t, eta, M_cut) {
  latent_intensity <- baseline_lambda(ts, magnitudes, eta, M_cut)
  detection_probability <- pmin(1, 10^(-(Mc_t - M_cut)))
  latent_intensity * detection_probability
}

tvc_lambda_at <- function(t_grid, ts, magnitudes, Mc_grid, eta, M_cut) {
  theta <- etas_to_natural(eta)
  intensity <- etas_lambda_at(t_grid, ts, magnitudes, theta, M_cut)

  detection_probability <- pmin(1, 10^(-(Mc_grid - M_cut)))
  intensity * detection_probability
}

baseline_log_lik <- function(ts, magnitudes, eta, T1, T2, M_cut) {
  lambdas <- baseline_lambda(ts, magnitudes, eta, M_cut)
  lambdas[lambdas <= 0] <- 1e-6

  target_idx <- ts > T1 & ts <= T2
  log_sum <- sum(log(lambdas[target_idx]))

  theta <- etas_to_natural(eta)
  log_sum - etas_integrated_intensity(
    T1, T2, ts, magnitudes, theta, M_cut
  )
}

tvc_log_lik <- function(ts,
                        magnitudes,
                        Mc_t,
                        eta,
                        T1,
                        T2,
                        M_cut,
                        M_trigger,
                        gamma,
                        tau,
                        mode = "sum",
                        grid_n = 400) {
  lambdas <- tvc_lambda(ts, magnitudes, Mc_t, eta, M_cut)
  lambdas[lambdas <= 0] <- 1e-6

  target_idx <- ts > T1 & ts <= T2
  log_sum <- sum(log(lambdas[target_idx]))

  t_grid <- seq(T1, T2, length.out = grid_n)
  Mc_grid <- compute_mct_at(
    t_eval = t_grid,
    event_ts = ts,
    event_magnitudes = magnitudes,
    M_cut = M_cut,
    M_trigger = M_trigger,
    gamma = gamma,
    tau = tau,
    mode = mode
  )

  lambda_grid <- tvc_lambda_at(t_grid, ts, magnitudes, Mc_grid, eta, M_cut)
  dt_grid <- diff(t_grid)
  integrated_intensity <- sum(
    0.5 * (lambda_grid[-1] + lambda_grid[-length(lambda_grid)]) * dt_grid
  )

  log_sum - integrated_intensity
}

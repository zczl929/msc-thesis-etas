# Completeness-curve estimation and quadrature helpers.

source("R/etas_likelihood.R")
source("R/loglinear_incompleteness.R")

estimate_loglinear_plugin <- function(catalogue,
                                      mainshock_time,
                                      mainshock_magnitude,
                                      M0,
                                      b = 1,
                                      H_bounds = c(0.2, 1.5),
                                      G_bounds = c(3, 7),
                                      grid_n = 2001L) {
  if (b <= 0) stop("b must be positive.")
  x <- catalogue[catalogue$ts > mainshock_time, , drop = FALSE]
  if (nrow(x) == 0L) stop("No post-mainshock events available.")
  dt <- x$ts - mainshock_time
  mag <- x$magnitudes
  beta <- b * log(10)

  evaluate_H <- function(H) {
    # Under deterministic threshold censoring, every observed mark must lie
    # above Mc(t).  For fixed H this gives a lower bound on G.  The conditional
    # truncated-GR mark likelihood is maximised at the smallest feasible G.
    required_G <- mainshock_magnitude - H * log10(dt) - mag
    G <- max(G_bounds[1], max(required_G) + 1e-10)
    if (G > G_bounds[2]) return(c(G = G, log_lik = -Inf))
    Mc <- pmax(M0, mainshock_magnitude - G - H * log10(dt))
    if (any(mag < Mc - 1e-9)) return(c(G = G, log_lik = -Inf))
    log_lik <- sum(log(beta) - beta * (mag - Mc))
    c(G = G, log_lik = log_lik)
  }

  H_grid <- seq(H_bounds[1], H_bounds[2], length.out = grid_n)
  grid_eval <- t(vapply(H_grid, evaluate_H, numeric(2)))
  best <- which.max(grid_eval[, "log_lik"])
  if (!is.finite(grid_eval[best, "log_lik"])) {
    stop("No feasible plug-in G,H estimate within the supplied bounds.")
  }

  # Refine locally while preserving the global grid search as protection
  # against the non-smooth support boundary of a hard-threshold model.
  left <- H_grid[max(1L, best - 2L)]
  right <- H_grid[min(length(H_grid), best + 2L)]
  refined <- optimize(
    function(H) -evaluate_H(H)[["log_lik"]],
    interval = c(left, right),
    tol = 1e-8
  )
  final <- evaluate_H(refined$minimum)

  list(
    estimate = c(G = unname(final[["G"]]), H = refined$minimum),
    log_lik = unname(final[["log_lik"]]),
    bounds = list(G = G_bounds, H = H_bounds),
    method = "conditional truncated-GR hard-threshold plug-in MLE",
    n_post_mainshock = nrow(x)
  )
}

estimate_loglinear_multitrigger_plugin <- function(catalogue,
                                                   M0,
                                                   M_trigger = 5,
                                                   b = 1,
                                                   H_bounds = c(0.2, 1.5),
                                                   G_bounds = c(3, 7),
                                                   grid_n = 2001L) {
  x <- catalogue[order(catalogue$ts), , drop = FALSE]
  triggers <- x[x$magnitudes >= M_trigger, c("ts", "magnitudes"), drop = FALSE]
  if (!nrow(triggers)) stop("No observed multi-trigger events available.")
  beta <- b * log(10)

  evaluate_H <- function(H) {
    required <- numeric(0)
    for (i in seq_len(nrow(x))) {
      past <- triggers$ts < x$ts[i]
      if (!any(past)) next
      dt <- x$ts[i] - triggers$ts[past]
      required <- c(
        required,
        max(triggers$magnitudes[past] - H * log10(dt) - x$magnitudes[i])
      )
    }
    G <- max(G_bounds[1], max(required, -Inf) + 1e-10)
    if (!is.finite(G) || G > G_bounds[2]) return(c(G = G, log_lik = -Inf))
    Mc <- compute_loglinear_mct_at(
      x$ts, triggers$ts, triggers$magnitudes, M0, G, H,
      include_current = FALSE
    )
    if (any(x$magnitudes < Mc - 1e-9)) return(c(G = G, log_lik = -Inf))
    c(G = G, log_lik = sum(log(beta) - beta * (x$magnitudes - Mc)))
  }

  H_grid <- seq(H_bounds[1], H_bounds[2], length.out = grid_n)
  grid_eval <- t(vapply(H_grid, evaluate_H, numeric(2)))
  best <- which.max(grid_eval[, "log_lik"])
  if (!is.finite(grid_eval[best, "log_lik"])) stop("No feasible multi-trigger G,H estimate.")
  left <- H_grid[max(1L, best - 2L)]
  right <- H_grid[min(length(H_grid), best + 2L)]
  refined <- optimize(
    function(H) -evaluate_H(H)[["log_lik"]], c(left, right), tol = 1e-8
  )
  final <- evaluate_H(refined$minimum)
  list(
    estimate = c(G = unname(final[["G"]]), H = refined$minimum),
    log_lik = unname(final[["log_lik"]]),
    method = "multi-trigger truncated-GR hard-threshold plug-in MLE",
    M_trigger = M_trigger,
    n_triggers = nrow(triggers)
  )
}

make_event_split_quadrature <- function(T1, T2, event_ts, grid_n = 3000L) {
  breaks <- sort(unique(c(T1, event_ts[event_ts > T1 & event_ts < T2], T2)))
  max_step <- (T2 - T1) / grid_n
  nodes <- c(-sqrt(3 / 5), 0, sqrt(3 / 5))
  node_weights <- c(5 / 9, 8 / 9, 5 / 9)
  grid_parts <- vector("list", length(breaks) - 1L)
  weight_parts <- vector("list", length(breaks) - 1L)
  for (i in seq_len(length(breaks) - 1L)) {
    sub_breaks <- seq(
      breaks[i], breaks[i + 1L],
      length.out = max(2L, ceiling((breaks[i + 1L] - breaks[i]) / max_step) + 1L)
    )
    left <- sub_breaks[-length(sub_breaks)]
    right <- sub_breaks[-1L]
    mid <- (left + right) / 2
    half <- (right - left) / 2
    grid_parts[[i]] <- as.vector(outer(half, nodes, `*`) + mid)
    weight_parts[[i]] <- as.vector(outer(half, node_weights, `*`))
  }
  list(
    t = unlist(grid_parts, use.names = FALSE),
    weight = unlist(weight_parts, use.names = FALSE)
  )
}

loglinear_time_log_lik <- function(ts,
                                   magnitudes,
                                   eta,
                                   T1,
                                   T2,
                                   M0,
                                   mainshock_time,
                                   mainshock_magnitude,
                                   G,
                                   H,
                                   b = 1,
                                   correction = c("full", "kamranzad"),
                                   history_ts = ts,
                                   history_magnitudes = magnitudes,
                                   quadrature = NULL,
                                   grid_n = 3000L) {
  correction <- match.arg(correction)
  ensure_fast_rate_kernel()
  theta <- etas_to_natural(eta)
  target <- ts > T1 & ts <= T2

  latent_event <- weighted_etas_rate_cpp(
    ts, history_ts, history_magnitudes, rep(1, length(history_ts)),
    theta[["mu"]], theta[["K"]], theta[["alpha"]],
    theta[["c"]], theta[["p"]], M0
  )
  psi_event <- loglinear_detectable_fraction(
    ts, mainshock_time, mainshock_magnitude, M0, G, H, b
  )
  observed_event <- if (correction == "full") {
    latent_event * psi_event
  } else {
    theta[["mu"]] + (latent_event - theta[["mu"]]) * psi_event
  }

  if (is.null(quadrature)) {
    quadrature <- make_event_split_quadrature(T1, T2, history_ts, grid_n)
  }
  latent_grid <- weighted_etas_rate_cpp(
    quadrature$t, history_ts, history_magnitudes, rep(1, length(history_ts)),
    theta[["mu"]], theta[["K"]], theta[["alpha"]],
    theta[["c"]], theta[["p"]], M0
  )
  psi_grid <- loglinear_detectable_fraction(
    quadrature$t, mainshock_time, mainshock_magnitude, M0, G, H, b
  )
  observed_grid <- if (correction == "full") {
    latent_grid * psi_grid
  } else {
    theta[["mu"]] + (latent_grid - theta[["mu"]]) * psi_grid
  }

  sum(log(pmax(observed_event[target], 1e-300))) -
    sum(observed_grid * quadrature$weight)
}

fit_loglinear_fixed_mle <- function(catalogue,
                                    T1,
                                    T2,
                                    M0,
                                    mainshock_time,
                                    mainshock_magnitude,
                                    G,
                                    H,
                                    b = 1,
                                    correction = c("full", "kamranzad"),
                                    start,
                                    history_catalogue = NULL,
                                    grid_n = 3000L) {
  correction <- match.arg(correction)
  x <- catalogue[order(catalogue$ts, catalogue$ID), , drop = FALSE]
  history <- if (is.null(history_catalogue)) x else {
    history_catalogue[order(history_catalogue$ts, history_catalogue$ID), , drop = FALSE]
  }
  quadrature <- make_event_split_quadrature(T1, T2, history$ts, grid_n)
  objective <- function(eta) {
    value <- loglinear_time_log_lik(
      x$ts, x$magnitudes, eta, T1, T2, M0,
      mainshock_time, mainshock_magnitude, G, H, b,
      correction, history$ts, history$magnitudes, quadrature, grid_n
    )
    if (!is.finite(value)) return(1e100)
    -value
  }
  fit <- optim(
    etas_to_internal(start), objective, method = "L-BFGS-B",
    lower = log(c(1e-4, 1e-4, 0.05, 1e-5, 0.01)),
    upper = log(c(20, 2, 0.98 * b, 2, 3)),
    control = list(maxit = 300, factr = 1e8)
  )
  list(
    model = paste0("loglinear_", correction),
    estimate = etas_to_natural(fit$par),
    G = G,
    H = H,
    log_lik = -fit$value,
    convergence = fit$convergence,
    message = fit$message,
    n = nrow(x)
  )
}

fit_loglinear_joint_laplace <- function(catalogue,
                                        T1,
                                        T2,
                                        M0,
                                        mainshock_time,
                                        mainshock_magnitude,
                                        b = 1,
                                        start,
                                        G_start = 4.5,
                                        H_start = 0.75,
                                        G_bounds = c(3, 7),
                                        H_bounds = c(0.2, 1.5),
                                        n_draws = 2000L,
                                        grid_n = 2000L,
                                        seed = 2026L,
                                        fail_on_non_pd = FALSE) {
  x <- catalogue[order(catalogue$ts, catalogue$ID), , drop = FALSE]
  quadrature <- make_event_split_quadrature(T1, T2, x$ts, grid_n)
  lower_natural <- c(
    mu = 1e-4, K = 1e-4, alpha = 0.05, c = 0.001,
    p_minus_1 = 0.01, G = G_bounds[1], H = H_bounds[1]
  )
  upper_natural <- c(
    mu = 2, K = 2, alpha = 0.98 * b, c = 0.2,
    p_minus_1 = 1.5, G = G_bounds[2], H = H_bounds[2]
  )
  z0 <- c(etas_to_internal(start), log(G_start), log(H_start))

  log_prior_and_jacobian <- function(z, theta, G, H) {
    lp <- dgamma(theta[["mu"]], shape = 2, rate = 20, log = TRUE) +
      dlnorm(theta[["K"]], meanlog = log(0.21), sdlog = 0.8, log = TRUE) +
      dnorm(G, mean = 4.5, sd = 0.5, log = TRUE) +
      dlnorm(H, meanlog = log(0.75), sdlog = 0.35, log = TRUE)
    # Uniform priors for alpha, c and p contribute constants inside bounds.
    # Add the Jacobian for z -> (mu,K,alpha,c,p,G,H).
    lp + sum(z[c(1, 2, 3, 4, 6, 7)]) + z[5]
  }

  objective <- function(z) {
    theta <- etas_to_natural(z[1:5])
    G <- exp(z[6])
    H <- exp(z[7])
    ll <- loglinear_time_log_lik(
      x$ts, x$magnitudes, z[1:5], T1, T2, M0,
      mainshock_time, mainshock_magnitude, G, H, b,
      correction = "full", quadrature = quadrature, grid_n = grid_n
    )
    value <- ll + log_prior_and_jacobian(z, theta, G, H)
    if (!is.finite(value)) return(1e100)
    -value
  }

  fit <- optim(
    z0, objective, method = "L-BFGS-B",
    lower = log(lower_natural),
    upper = log(upper_natural),
    control = list(maxit = 800, factr = 1e5, pgtol = 1e-7)
  )
  hessian <- optimHess(fit$par, objective)
  hessian <- 0.5 * (hessian + t(hessian))
  eig <- eigen(hessian, symmetric = TRUE)
  positive_definite <- all(eig$values > 1e-8)
  if (!positive_definite) {
    negative_direction <- eig$vectors[, which.min(eig$values)]
    profile_steps <- c(-0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05)
    profile_objective <- vapply(
      profile_steps,
      function(step) {
        proposal <- pmin(
          pmax(fit$par + step * negative_direction, log(lower_natural)),
          log(upper_natural)
        )
        objective(proposal)
      },
      numeric(1)
    )
    diagnostic <- list(
      method = "joint Bayesian log-linear ETAS Laplace identifiability pilot",
      approximation = "Laplace approximation unavailable: non-positive Hessian",
      map = c(
        etas_to_natural(fit$par[1:5]),
        G = exp(fit$par[6]), H = exp(fit$par[7])
      ),
      convergence = fit$convergence,
      message = fit$message,
      log_posterior_at_map = -fit$value,
      hessian_eigenvalues = eig$values,
      negative_eigenvector = setNames(
        negative_direction,
        c("log_mu", "log_K", "log_alpha", "log_c", "log_p_minus_1", "log_G", "log_H")
      ),
      negative_direction_profile = data.frame(
        step = profile_steps,
        negative_log_posterior = profile_objective,
        delta_from_map = profile_objective - fit$value
      ),
      positive_definite = FALSE,
      internal_map = fit$par,
      warning = paste(
        "The local Gaussian joint posterior is invalid.",
        "This may indicate weak identification, a non-smooth mode, or",
        "insufficient optimisation and must be resolved before replication."
      )
    )
    if (fail_on_non_pd) {
      stop(
        "Joint Laplace Hessian is not positive definite; minimum eigenvalue = ",
        min(eig$values)
      )
    }
    return(diagnostic)
  }
  covariance <- eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
  covariance <- 0.5 * (covariance + t(covariance))
  set.seed(seed)
  standard <- matrix(rnorm(n_draws * 7L), ncol = 7L)
  chol_cov <- chol(covariance)
  z_draws <- sweep(standard %*% chol_cov, 2, fit$par, `+`)

  natural_draws <- cbind(
    mu = exp(z_draws[, 1]),
    K = exp(z_draws[, 2]),
    alpha = exp(z_draws[, 3]),
    c = exp(z_draws[, 4]),
    p = 1 + exp(z_draws[, 5]),
    G = exp(z_draws[, 6]),
    H = exp(z_draws[, 7])
  )
  within_bounds <- apply(
    natural_draws,
    1,
    function(v) all(v[c("mu", "K", "alpha", "c", "p", "G", "H")] >= c(
      lower_natural[c("mu", "K", "alpha", "c")],
      1 + lower_natural[["p_minus_1"]], G_bounds[1], H_bounds[1]
    )) && all(v[c("mu", "K", "alpha", "c", "p", "G", "H")] <= c(
      upper_natural[c("mu", "K", "alpha", "c")],
      1 + upper_natural[["p_minus_1"]], G_bounds[2], H_bounds[2]
    ))
  )
  natural_draws <- natural_draws[within_bounds, , drop = FALSE]
  summaries <- t(apply(natural_draws, 2, function(v) c(
    mean = mean(v),
    median = median(v),
    q025 = quantile(v, 0.025, names = FALSE),
    q975 = quantile(v, 0.975, names = FALSE)
  )))

  map <- c(etas_to_natural(fit$par[1:5]), G = exp(fit$par[6]), H = exp(fit$par[7]))
  list(
    method = "joint Bayesian log-linear ETAS Laplace identifiability pilot",
    approximation = "Laplace approximation to the marginal-time posterior",
    map = map,
    posterior_summary = summaries,
    posterior_draws = natural_draws,
    convergence = fit$convergence,
    message = fit$message,
    log_posterior_at_map = -fit$value,
    hessian_eigenvalues = eig$values,
    positive_definite = TRUE,
    covariance_internal = covariance,
    retained_draw_fraction = nrow(natural_draws) / n_draws,
    warning = paste(
      "This is an identifiability gate, not the final joint model.",
      "It uses the unmarked marginal-time censorship likelihood and ignores",
      "triggering by latent missing events."
    )
  )
}

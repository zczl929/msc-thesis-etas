# Exact-posterior MCMC targets for the frozen Ridgecrest analysis.
#
# The fitted process remains the observed-history Plug-in TVC-ETAS model.  It
# does not augment unobserved earthquakes.  The compensator is evaluated as an
# exact ETAS compensator plus a high-order quadrature correction restricted to
# the short intervals in which the detection probability is below one.

source("R/mcmc_etas.R")
source("R/loglinear_incompleteness.R")

make_multitrigger_adjustment_quadrature <- function(
    T1, T2, event_ts, trigger_ts, trigger_magnitudes, M0, G, H, b = 1,
    nodes_per_interval = 16L) {
  recovery_end <- loglinear_recovery_end(
    trigger_ts, trigger_magnitudes, M0, G, H
  )
  active <- recovery_end > T1 & trigger_ts < T2
  if (!any(active)) {
    return(list(t = numeric(), weight = numeric(), psi = numeric()))
  }
  trigger_ts <- trigger_ts[active]
  trigger_magnitudes <- trigger_magnitudes[active]
  recovery_end <- recovery_end[active]

  in_active_window <- vapply(
    event_ts,
    function(t) any(t > trigger_ts & t < recovery_end),
    logical(1)
  )
  breaks <- sort(unique(c(
    T1, T2,
    pmax(T1, trigger_ts),
    pmin(T2, recovery_end),
    event_ts[in_active_window & event_ts > T1 & event_ts < T2]
  )))
  rule <- gauss_legendre_rule(nodes_per_interval)
  t_parts <- list()
  weight_parts <- list()
  part <- 0L
  for (i in seq_len(length(breaks) - 1L)) {
    left <- breaks[i]
    right <- breaks[i + 1L]
    if (right <= left) next
    midpoint <- (left + right) / 2
    midpoint_psi <- loglinear_detectable_fraction(
      midpoint, trigger_ts, trigger_magnitudes, M0, G, H, b
    )
    if (midpoint_psi >= 1 - 1e-14) next
    half_width <- (right - left) / 2
    part <- part + 1L
    t_parts[[part]] <- midpoint + half_width * rule$nodes
    weight_parts[[part]] <- half_width * rule$weights
  }
  if (!part) {
    return(list(t = numeric(), weight = numeric(), psi = numeric()))
  }
  t <- unlist(t_parts, use.names = FALSE)
  list(
    t = t,
    weight = unlist(weight_parts, use.names = FALSE),
    psi = loglinear_detectable_fraction(
      t, trigger_ts, trigger_magnitudes, M0, G, H, b
    )
  )
}

ensure_multitrigger_plugin_loglik_cpp <- local({
  loaded <- FALSE
  function() {
    if (loaded) return(invisible(TRUE))
    ensure_exact_etas_loglik_cpp()
    if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
    Rcpp::cppFunction(code = '
      double multitrigger_plugin_loglik_cpp(
          Rcpp::NumericVector event_t, Rcpp::NumericVector mag,
          Rcpp::LogicalVector include_likelihood,
          Rcpp::NumericVector event_psi,
          Rcpp::NumericVector quad_t, Rcpp::NumericVector quad_weight,
          Rcpp::NumericVector quad_psi,
          double T1, double T2, double M0,
          double mu, double K, double alpha, double c, double p) {
        int n = event_t.size();
        if (!(mu > 0.0 && K > 0.0 && alpha >= 0.0 &&
              c > 0.0 && p > 1.0)) return R_NegInf;
        if (event_psi.size() != n ||
            quad_t.size() != quad_weight.size() ||
            quad_t.size() != quad_psi.size()) return R_NegInf;

        double norm = (p - 1.0) / c;
        Rcpp::NumericVector productivity(n);
        for (int j = 0; j < n; ++j) {
          productivity[j] =
            K * std::pow(10.0, alpha * (mag[j] - M0));
        }

        double event_term = 0.0;
        for (int i = 0; i < n; ++i) {
          if (!(include_likelihood[i] && event_t[i] > T1 &&
                event_t[i] <= T2)) continue;
          double rate = mu;
          for (int j = 0; j < i; ++j) {
            if (event_t[j] >= event_t[i]) break;
            double dt = event_t[i] - event_t[j];
            rate += productivity[j] * norm *
              std::pow(1.0 + dt / c, -p);
          }
          if (!R_finite(rate) || rate <= 0.0 ||
              !R_finite(event_psi[i]) || event_psi[i] <= 0.0) {
            return R_NegInf;
          }
          event_term += std::log(rate) + std::log(event_psi[i]);
        }

        double integral = mu * (T2 - T1);
        for (int j = 0; j < n; ++j) {
          if (event_t[j] >= T2) continue;
          double lower = std::max(T1, event_t[j]);
          double upper_dt = T2 - event_t[j];
          double lower_dt = lower - event_t[j];
          double upper_mass =
            1.0 - std::pow(1.0 + upper_dt / c, 1.0 - p);
          double lower_mass = lower_dt > 0.0 ?
            1.0 - std::pow(1.0 + lower_dt / c, 1.0 - p) : 0.0;
          integral += productivity[j] * (upper_mass - lower_mass);
        }

        for (int q = 0; q < quad_t.size(); ++q) {
          double rate = mu;
          for (int j = 0; j < n; ++j) {
            if (event_t[j] >= quad_t[q]) break;
            double dt = quad_t[q] - event_t[j];
            rate += productivity[j] * norm *
              std::pow(1.0 + dt / c, -p);
          }
          integral +=
            quad_weight[q] * (quad_psi[q] - 1.0) * rate;
        }
        if (!R_finite(integral)) return R_NegInf;
        return event_term - integral;
      }', env = globalenv())
    loaded <<- TRUE
    invisible(TRUE)
  }
})

make_multitrigger_plugin_log_posterior <- function(
    catalogue, T1, T2, M0, prior, G, H, b = 1, M_trigger = 5,
    nodes_per_interval = 16L) {
  ensure_multitrigger_plugin_loglik_cpp()
  x <- catalogue[order(catalogue$ts, catalogue$ID), , drop = FALSE]
  include <- if ("include_likelihood" %in% names(x)) {
    as.logical(x$include_likelihood)
  } else {
    rep(TRUE, nrow(x))
  }
  triggers <- x$magnitudes >= M_trigger
  event_psi <- loglinear_detectable_fraction(
    x$ts, x$ts[triggers], x$magnitudes[triggers], M0, G, H, b
  )
  quadrature <- make_multitrigger_adjustment_quadrature(
    T1, T2, x$ts, x$ts[triggers], x$magnitudes[triggers],
    M0, G, H, b, nodes_per_interval
  )
  log_posterior <- function(z) {
    theta <- etas_z_to_theta(z, prior)
    lp <- etas_log_prior_natural(theta, prior) +
      etas_log_jacobian(z, prior)
    if (!is.finite(lp)) return(-Inf)
    ll <- multitrigger_plugin_loglik_cpp(
      x$ts, x$magnitudes, include, event_psi,
      quadrature$t, quadrature$weight, quadrature$psi,
      T1, T2, M0,
      theta[["mu"]], theta[["K"]], theta[["alpha"]],
      theta[["c"]], theta[["p"]]
    )
    ll + lp
  }
  attr(log_posterior, "quadrature") <- quadrature
  attr(log_posterior, "event_psi") <- event_psi
  log_posterior
}

make_fixedtrigger_plugin_log_posterior <- function(
    catalogue, T1, T2, M0, prior, G, H,
    trigger_time, trigger_magnitude, b = 1,
    nodes_per_interval = 16L) {
  ensure_multitrigger_plugin_loglik_cpp()
  if (!all(is.finite(c(
    G, H, trigger_time, trigger_magnitude, b
  )))) {
    stop("G, H, trigger_time, trigger_magnitude and b must be finite.")
  }
  x <- catalogue[order(catalogue$ts, catalogue$ID), , drop = FALSE]
  include <- if ("include_likelihood" %in% names(x)) {
    as.logical(x$include_likelihood)
  } else {
    rep(TRUE, nrow(x))
  }
  event_psi <- loglinear_detectable_fraction(
    x$ts, trigger_time, trigger_magnitude, M0, G, H, b
  )
  quadrature <- make_multitrigger_adjustment_quadrature(
    T1, T2, x$ts, trigger_time, trigger_magnitude,
    M0, G, H, b, nodes_per_interval
  )
  log_posterior <- function(z) {
    theta <- etas_z_to_theta(z, prior)
    lp <- etas_log_prior_natural(theta, prior) +
      etas_log_jacobian(z, prior)
    if (!is.finite(lp)) return(-Inf)
    ll <- multitrigger_plugin_loglik_cpp(
      x$ts, x$magnitudes, include, event_psi,
      quadrature$t, quadrature$weight, quadrature$psi,
      T1, T2, M0,
      theta[["mu"]], theta[["K"]], theta[["alpha"]],
      theta[["c"]], theta[["p"]]
    )
    ll + lp
  }
  attr(log_posterior, "quadrature") <- quadrature
  attr(log_posterior, "event_psi") <- event_psi
  attr(log_posterior, "fixed_trigger") <- c(
    time = trigger_time, magnitude = trigger_magnitude
  )
  log_posterior
}

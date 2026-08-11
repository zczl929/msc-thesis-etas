# Exact-likelihood random-walk MCMC for the final synthetic and Ridgecrest
# analyses.

source("R/etas_likelihood.R")

ensure_exact_etas_loglik_cpp <- local({
  loaded <- FALSE
  function() {
    if (loaded) return(invisible(TRUE))
    if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
    Rcpp::cppFunction(code = '
      double exact_etas_loglik_cpp(
          Rcpp::NumericVector event_t, Rcpp::NumericVector mag,
          Rcpp::LogicalVector include_likelihood,
          double T1, double T2, double M0,
          double mu, double K, double alpha, double c, double p) {
        int n = event_t.size();
        if (!(mu > 0.0 && K > 0.0 && alpha >= 0.0 &&
              c > 0.0 && p > 1.0)) return R_NegInf;
        double norm = (p - 1.0) / c;
        Rcpp::NumericVector productivity(n);
        for (int j = 0; j < n; ++j) {
          productivity[j] =
            K * std::pow(10.0, alpha * (mag[j] - M0));
        }
        double event_term = 0.0;
        for (int i = 0; i < n; ++i) {
          if (!(include_likelihood[i] && event_t[i] > T1 &&
                event_t[i] < T2)) continue;
          double rate = mu;
          for (int j = 0; j < i; ++j) {
            if (event_t[j] >= event_t[i]) break;
            double dt = event_t[i] - event_t[j];
            rate += productivity[j] * norm *
              std::pow(1.0 + dt / c, -p);
          }
          if (!R_finite(rate) || rate <= 0.0) return R_NegInf;
          event_term += std::log(rate);
        }
        double integral = mu * (T2 - T1);
        for (int j = 0; j < n; ++j) {
          if (event_t[j] >= T2) continue;
          double lower = std::max(T1, event_t[j]);
          double upper_dt = T2 - event_t[j];
          double lower_dt = lower - event_t[j];
          double upper_mass = 1.0 - std::pow(1.0 + upper_dt / c, 1.0 - p);
          double lower_mass = lower_dt > 0.0 ?
            1.0 - std::pow(1.0 + lower_dt / c, 1.0 - p) : 0.0;
          integral += productivity[j] * (upper_mass - lower_mass);
        }
        if (!R_finite(integral)) return R_NegInf;
        return event_term - integral;
      }', env = globalenv())
    loaded <<- TRUE
    invisible(TRUE)
  }
})

bounded_forward <- function(z, lower, upper) {
  lower + (upper - lower) * plogis(z)
}

bounded_inverse <- function(x, lower, upper) {
  qlogis((x - lower) / (upper - lower))
}

etas_z_to_theta <- function(z, prior) {
  setNames(c(
    exp(z[1]),
    exp(z[2]),
    bounded_forward(z[3], prior$a_alpha, prior$b_alpha),
    bounded_forward(z[4], prior$a_c, prior$b_c),
    bounded_forward(z[5], prior$a_p, prior$b_p)
  ), c("mu", "K", "alpha", "c", "p"))
}

etas_theta_to_z <- function(theta, prior) {
  setNames(c(
    log(theta[["mu"]]),
    log(theta[["K"]]),
    bounded_inverse(theta[["alpha"]], prior$a_alpha, prior$b_alpha),
    bounded_inverse(theta[["c"]], prior$a_c, prior$b_c),
    bounded_inverse(theta[["p"]], prior$a_p, prior$b_p)
  ), c("mu", "K", "alpha", "c", "p"))
}

etas_log_jacobian <- function(z, prior) {
  bounded_jacobian <- function(value, lower, upper) {
    s <- plogis(value)
    log(upper - lower) + log(s) + log1p(-s)
  }
  z[1] + z[2] +
    bounded_jacobian(z[3], prior$a_alpha, prior$b_alpha) +
    bounded_jacobian(z[4], prior$a_c, prior$b_c) +
    bounded_jacobian(z[5], prior$a_p, prior$b_p)
}

etas_log_prior_natural <- function(theta, prior) {
  dgamma(theta[["mu"]], shape = prior$a_mu, rate = prior$b_mu, log = TRUE) +
    dlnorm(theta[["K"]], meanlog = prior$a_K, sdlog = prior$b_K, log = TRUE) +
    dunif(
      theta[["alpha"]], prior$a_alpha, prior$b_alpha, log = TRUE
    ) +
    dunif(theta[["c"]], prior$a_c, prior$b_c, log = TRUE) +
    dunif(theta[["p"]], prior$a_p, prior$b_p, log = TRUE)
}

make_exact_etas_log_posterior <- function(catalogue, T1, T2, M0, prior) {
  ensure_exact_etas_loglik_cpp()
  x <- catalogue[order(catalogue$ts, catalogue$ID), , drop = FALSE]
  include <- if ("include_likelihood" %in% names(x)) {
    as.logical(x$include_likelihood)
  } else if ("gen" %in% names(x)) {
    x$gen != -1L
  } else {
    rep(TRUE, nrow(x))
  }
  function(z) {
    theta <- etas_z_to_theta(z, prior)
    lp <- etas_log_prior_natural(theta, prior) +
      etas_log_jacobian(z, prior)
    if (!is.finite(lp)) return(-Inf)
    ll <- exact_etas_loglik_cpp(
      x$ts, x$magnitudes, include, T1, T2, M0,
      theta[["mu"]], theta[["K"]], theta[["alpha"]],
      theta[["c"]], theta[["p"]]
    )
    ll + lp
  }
}

gauss_legendre_rule <- function(n = 16L) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least 2.")
  index <- seq_len(n - 1L)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  eig <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(eig$values)
  list(
    nodes = eig$values[ordering],
    weights = 2 * eig$vectors[1L, ordering]^2
  )
}

make_detection_adjustment_quadrature <- function(
    T1, T2, event_ts, mainshock_time, mainshock_magnitude,
    M0, G, H, nodes_per_interval = 16L) {
  recovery_end <- mainshock_time +
    10^((mainshock_magnitude - G - M0) / H)
  lower <- max(T1, mainshock_time)
  upper <- min(T2, recovery_end)
  if (!is.finite(upper) || upper <= lower) {
    return(list(t = numeric(), weight = numeric()))
  }
  breaks <- sort(unique(c(
    lower, event_ts[event_ts > lower & event_ts < upper], upper
  )))
  rule <- gauss_legendre_rule(nodes_per_interval)
  t_parts <- vector("list", length(breaks) - 1L)
  weight_parts <- vector("list", length(breaks) - 1L)
  for (i in seq_len(length(breaks) - 1L)) {
    midpoint <- (breaks[i] + breaks[i + 1L]) / 2
    half_width <- (breaks[i + 1L] - breaks[i]) / 2
    t_parts[[i]] <- midpoint + half_width * rule$nodes
    weight_parts[[i]] <- half_width * rule$weights
  }
  list(
    t = unlist(t_parts, use.names = FALSE),
    weight = unlist(weight_parts, use.names = FALSE)
  )
}

ensure_plugin_etas_loglik_cpp <- local({
  loaded <- FALSE
  function() {
    if (loaded) return(invisible(TRUE))
    ensure_exact_etas_loglik_cpp()
    Rcpp::cppFunction(code = '
      double plugin_etas_loglik_cpp(
          Rcpp::NumericVector event_t, Rcpp::NumericVector mag,
          Rcpp::LogicalVector include_likelihood,
          Rcpp::NumericVector quad_t, Rcpp::NumericVector quad_weight,
          double T1, double T2, double M0,
          double mainshock_time, double mainshock_magnitude,
          double G, double H, double b,
          double mu, double K, double alpha, double c, double p) {
        int n = event_t.size();
        if (!(mu > 0.0 && K > 0.0 && alpha >= 0.0 &&
              c > 0.0 && p > 1.0 && G > 0.0 && H > 0.0 && b > 0.0)) {
          return R_NegInf;
        }
        double norm = (p - 1.0) / c;
        Rcpp::NumericVector productivity(n);
        for (int j = 0; j < n; ++j) {
          productivity[j] =
            K * std::pow(10.0, alpha * (mag[j] - M0));
        }
        double event_term = 0.0;
        for (int i = 0; i < n; ++i) {
          if (!(include_likelihood[i] && event_t[i] > T1 &&
                event_t[i] < T2)) continue;
          double rate = mu;
          for (int j = 0; j < i; ++j) {
            if (event_t[j] >= event_t[i]) break;
            double dt = event_t[i] - event_t[j];
            rate += productivity[j] * norm *
              std::pow(1.0 + dt / c, -p);
          }
          double psi = 1.0;
          if (event_t[i] > mainshock_time) {
            double dt_main = event_t[i] - mainshock_time;
            double mc = std::max(
              M0, mainshock_magnitude - G - H * std::log10(dt_main)
            );
            psi = std::pow(10.0, -b * std::max(mc - M0, 0.0));
          }
          if (!R_finite(rate) || rate <= 0.0 || psi <= 0.0) {
            return R_NegInf;
          }
          event_term += std::log(rate) + std::log(psi);
        }
        double integral = mu * (T2 - T1);
        for (int j = 0; j < n; ++j) {
          if (event_t[j] >= T2) continue;
          double lower = std::max(T1, event_t[j]);
          double upper_dt = T2 - event_t[j];
          double lower_dt = lower - event_t[j];
          double upper_mass = 1.0 - std::pow(1.0 + upper_dt / c, 1.0 - p);
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
          double dt_main = quad_t[q] - mainshock_time;
          double mc = std::max(
            M0, mainshock_magnitude - G - H * std::log10(dt_main)
          );
          double psi = std::pow(10.0, -b * std::max(mc - M0, 0.0));
          integral += quad_weight[q] * (psi - 1.0) * rate;
        }
        if (!R_finite(integral)) return R_NegInf;
        return event_term - integral;
      }', env = globalenv())
    loaded <<- TRUE
    invisible(TRUE)
  }
})

make_plugin_etas_log_posterior <- function(
    catalogue, T1, T2, M0, prior, mainshock_time, mainshock_magnitude,
    G, H, b = 1, nodes_per_interval = 16L) {
  ensure_plugin_etas_loglik_cpp()
  x <- catalogue[order(catalogue$ts, catalogue$ID), , drop = FALSE]
  include <- if ("include_likelihood" %in% names(x)) {
    as.logical(x$include_likelihood)
  } else if ("gen" %in% names(x)) {
    x$gen != -1L
  } else {
    rep(TRUE, nrow(x))
  }
  quadrature <- make_detection_adjustment_quadrature(
    T1, T2, x$ts, mainshock_time, mainshock_magnitude,
    M0, G, H, nodes_per_interval
  )
  function(z) {
    theta <- etas_z_to_theta(z, prior)
    lp <- etas_log_prior_natural(theta, prior) +
      etas_log_jacobian(z, prior)
    if (!is.finite(lp)) return(-Inf)
    ll <- plugin_etas_loglik_cpp(
      x$ts, x$magnitudes, include,
      quadrature$t, quadrature$weight,
      T1, T2, M0, mainshock_time, mainshock_magnitude,
      G, H, b,
      theta[["mu"]], theta[["K"]], theta[["alpha"]],
      theta[["c"]], theta[["p"]]
    )
    ll + lp
  }
}

run_exact_etas_chain <- function(log_posterior, initial, proposal_cov,
                                 iterations = 12000L, warmup = 4000L,
                                 thin = 4L, seed = 1L) {
  set.seed(seed)
  d <- length(initial)
  proposal_cov <- proposal_cov + diag(1e-8, d)
  chol_cov <- chol(proposal_cov)
  current <- initial
  current_lp <- log_posterior(current)
  if (!is.finite(current_lp)) stop("Initial MCMC state has zero density.")
  kept <- matrix(
    NA_real_,
    nrow = floor((iterations - warmup) / thin),
    ncol = d,
    dimnames = list(NULL, names(initial))
  )
  accepted <- 0L
  accepted_window <- 0L
  scale <- 2.38 / sqrt(d)
  keep_index <- 0L
  for (iteration in seq_len(iterations)) {
    proposal <- current +
      as.numeric(rnorm(d) %*% chol_cov) * scale
    proposal_lp <- log_posterior(proposal)
    accept <- is.finite(proposal_lp) &&
      log(runif(1)) < proposal_lp - current_lp
    if (accept) {
      current <- proposal
      current_lp <- proposal_lp
      accepted <- accepted + 1L
      accepted_window <- accepted_window + 1L
    }
    if (iteration <= warmup && iteration %% 100L == 0L) {
      rate <- accepted_window / 100
      scale <- scale * exp(min(0.25, max(-0.25, rate - 0.234)))
      accepted_window <- 0L
    }
    if (iteration > warmup && (iteration - warmup) %% thin == 0L) {
      keep_index <- keep_index + 1L
      kept[keep_index, ] <- current
    }
  }
  list(
    z = kept,
    acceptance = accepted / iterations,
    final_scale = scale
  )
}

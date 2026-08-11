# ============================================================

if (!exists("etas_productivity", mode = "function")) {
  source("R/etas_likelihood.R")
}
# simulate_etas.R
#
# Custom implementation for temporal ETAS simulation.
# Generates a complete synthetic catalogue under known parameter
# values theta = (mu, K, c, p, alpha), as specified in the proposal
# Stage 1: Simulation-Based Bias Quantification.
#
# Simulation approach:
# Hawkes / ETAS branching-process representation.
#
# Conditional intensity:
#   lambda(t | H_t) = mu + sum_{ti < t} kappa(M_i) g(t - t_i)
#
# Productivity:
#   kappa(M_i) = K * 10^{alpha * (M_i - M_cut)}
#
# Normalised Omori kernel:
#   g(u) = ((p - 1) / c) * (1 + u / c)^(-p), u > 0
# ============================================================

#' Simulate a complete temporal ETAS catalogue
#'
#' @param theta named list or one-row data.frame with elements mu, K, alpha, c, p
#' @param T1,T2 numeric, simulation time window in days
#' @param M_cut numeric, magnitude cutoff, lower bound of catalogue
#' @param beta numeric, Gutenberg-Richter parameter, beta = b * log(10)
#' @param Ht optional data.frame with columns ts and magnitudes; seeded events
#'           such as known mainshocks. Pass NULL for none.
#' @param seed integer, optional, for reproducibility
#' @param max_events integer, safety cap to avoid explosive simulations
#'
#' @return data.frame with columns ts, magnitudes, gen, parent_id, ID
simulate_etas <- function(theta,
                          T1, T2,
                          M_cut,
                          beta = log(10),
                          Ht = NULL,
                          seed = NULL,
                          max_events = 1e6) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Basic checks ----------------------------------------------------
  if (T2 <= T1) {
    stop("T2 must be greater than T1.")
  }

  if (beta <= 0) {
    stop("beta must be positive.")
  }

  # Allow theta to be supplied as a one-row data.frame
  if (is.data.frame(theta)) {
    theta <- as.list(theta[1, , drop = FALSE])
  }

  required_pars <- c("mu", "K", "alpha", "c", "p")
  missing_pars <- setdiff(required_pars, names(theta))

  if (length(missing_pars) > 0) {
    stop("theta is missing: ", paste(missing_pars, collapse = ", "))
  }


  # ---- Unpack parameters ----------------------------------------------
  # mu    : background event rate
  # K     : productivity scale
  # alpha : magnitude effect on productivity
  # c_par : Omori time offset parameter
  # p_par : Omori temporal decay exponent

  mu    <- as.numeric(theta$mu)
  K     <- as.numeric(theta$K)
  alpha <- as.numeric(theta$alpha)
  c_par <- as.numeric(theta$c)
  p_par <- as.numeric(theta$p)

  if (mu < 0 || K < 0 || c_par <= 0 || p_par <= 1) {
    stop("Require mu >= 0, K >= 0, c > 0, and p > 1.")
  }


  # ---- Helper functions ------------------------------------------------

  # Productivity function:
  # This is where M_i and M_cut enter the model.
  theta_natural <- c(mu = mu, K = K, alpha = alpha, c = c_par, p = p_par)
  validate_etas_theta(theta_natural)

  productivity <- function(mag) {
    etas_productivity(mag, theta_natural, M_cut)
  }

  # CDF of the normalised Omori kernel:
  # g(u) = ((p - 1) / c) * (1 + u / c)^(-p)
  #
  # Therefore:
  # G(u) = integral_0^u g(s) ds
  #      = 1 - (1 + u / c)^(1 - p)
  omori_cdf <- function(u) {
    etas_omori_cdf(u, theta_natural)
  }

  # Draw waiting times from the Omori kernel conditional on being
  # within [0, max_u]. This is needed because the simulation window
  # ends at T2.
  sample_truncated_omori <- function(n, max_u) {

    if (n <= 0 || max_u <= 0) {
      return(numeric(0))
    }

    G_max <- omori_cdf(max_u)

    if (G_max <= 0) {
      return(numeric(0))
    }

    u <- runif(n, min = 0, max = G_max)

    # Inverse CDF:
    # u = 1 - (1 + x / c)^(1 - p)
    # x = c * ((1 - u)^(1 / (1 - p)) - 1)
    waits <- c_par * ((1 - u)^(1 / (1 - p_par)) - 1)

    return(waits)
  }

  # Gutenberg-Richter magnitudes above M_cut.
  # M = M_cut + Exponential(beta)
  sample_magnitudes <- function(n) {
    M_cut + rexp(n, rate = beta)
  }


  # ---- Step 1: Background events --------------------------------------
  # Homogeneous Poisson process with rate mu on [T1, T2].
  # Magnitudes are drawn from Gutenberg-Richter distribution above M_cut.

  duration <- T2 - T1
  n_bg <- rpois(1, mu * duration)

  if (n_bg > 0) {
    bg_events <- data.frame(
      ts         = runif(n_bg, min = T1, max = T2),
      magnitudes = sample_magnitudes(n_bg),
      gen        = 0L,
      parent_id  = NA_integer_
    )
  } else {
    bg_events <- data.frame(
      ts         = numeric(0),
      magnitudes = numeric(0),
      gen        = integer(0),
      parent_id  = integer(0)
    )
  }


  # ---- Step 2: Add seeded events Ht -----------------------------------
  # User-provided events, such as known mainshocks, are added to the
  # catalogue with gen = -1 to mark them as externally seeded.

  if (!is.null(Ht) && nrow(Ht) > 0) {

    if (!("ts" %in% names(Ht))) {
      stop("Ht must contain a column named 'ts'.")
    }

    if (!("magnitudes" %in% names(Ht))) {
      stop("Ht must contain a column named 'magnitudes'.")
    }

    seeded_events <- data.frame(
      ts         = as.numeric(Ht$ts),
      magnitudes = as.numeric(Ht$magnitudes),
      gen        = -1L,
      parent_id  = NA_integer_
    )

    # Keep only seeded events inside the simulation window.
    seeded_events <- seeded_events[
      seeded_events$ts >= T1 & seeded_events$ts <= T2,
      ,
      drop = FALSE
    ]

  } else {
    seeded_events <- data.frame(
      ts         = numeric(0),
      magnitudes = numeric(0),
      gen        = integer(0),
      parent_id  = integer(0)
    )
  }


  # ---- Step 3: Combine initial catalogue -------------------------------

  catalogue <- rbind(bg_events, seeded_events)

  if (nrow(catalogue) == 0) {
    catalogue$ID <- integer(0)
    catalogue <- catalogue[, c("ts", "magnitudes", "gen", "parent_id", "ID")]
    return(catalogue)
  }

  # Assign unique event IDs.
  catalogue$ID <- seq_len(nrow(catalogue))


  # ---- Step 4: Branching loop -----------------------------------------
  # Each event can trigger offspring.
  #
  # For parent event i:
  #   expected number of offspring before T2
  #   = kappa(M_i) * G(T2 - t_i)
  #
  # where:
  #   kappa(M_i) = K * 10^{alpha * (M_i - M_cut)}
  #   G(.) is the CDF of the Omori kernel.

  next_parent_index <- 1L

  while (next_parent_index <= nrow(catalogue)) {

    if (nrow(catalogue) > max_events) {
      stop(
        "Simulation exceeded max_events = ", max_events,
        ". The ETAS process may be explosive. Try smaller K, alpha, or T2."
      )
    }

    parent <- catalogue[next_parent_index, ]

    parent_time <- parent$ts
    parent_mag  <- parent$magnitudes
    parent_id   <- parent$ID
    parent_gen  <- parent$gen

    max_wait <- T2 - parent_time

    if (max_wait > 0) {

      kappa_i <- productivity(parent_mag)
      expected_offspring <- kappa_i * omori_cdf(max_wait)

      n_child <- rpois(1, expected_offspring)

      if (n_child > 0) {

        child_waits <- sample_truncated_omori(n_child, max_wait)
        child_times <- parent_time + child_waits

        # Numerical safety.
        keep <- child_times >= T1 & child_times <= T2

        if (any(keep)) {

          n_keep <- sum(keep)

          # If parent is an externally seeded event with gen = -1,
          # its direct offspring are treated as generation 1.
          # If parent is a background event with gen = 0,
          # its direct offspring are also generation 1.
          child_gen <- ifelse(parent_gen < 0, 1L, parent_gen + 1L)

          child_events <- data.frame(
            ts         = child_times[keep],
            magnitudes = sample_magnitudes(n_keep),
            gen        = as.integer(child_gen),
            parent_id  = as.integer(parent_id)
          )

          child_events$ID <- seq(
            from = nrow(catalogue) + 1L,
            length.out = n_keep
          )

          catalogue <- rbind(catalogue, child_events)
        }
      }
    }

    next_parent_index <- next_parent_index + 1L
  }


  # ---- Step 5: Sort and finalise ---------------------------------------
  # Sort catalogue by event time.
  # IDs are kept as birth IDs, so parent_id still points to the original ID.

  catalogue <- catalogue[order(catalogue$ts, catalogue$ID), ]
  rownames(catalogue) <- NULL

  catalogue <- catalogue[, c("ts", "magnitudes", "gen", "parent_id", "ID")]

  return(catalogue)
}

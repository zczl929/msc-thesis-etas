# Helmstetter/Kamranzad short-term catalogue incompleteness.
#
# For a trigger j and t > t_j,
#   Mc_j(t) = m_j - G - H log10(t - t_j).
# The catalogue-level threshold is max(M0, Mc_j(t)) over the selected
# triggers.  Kamranzad et al. (2025) use the seeded mainshock only in their
# primary likelihood; the general multi-trigger form is retained here for a
# later misspecification experiment.

select_loglinear_triggers <- function(catalogue,
                                      trigger_mode = c("seeded", "threshold", "all"),
                                      M_trigger = 5.5) {
  trigger_mode <- match.arg(trigger_mode)
  if (!all(c("ts", "magnitudes") %in% names(catalogue))) {
    stop("catalogue must contain ts and magnitudes.")
  }

  keep <- switch(
    trigger_mode,
    seeded = {
      if (!"gen" %in% names(catalogue)) {
        stop("trigger_mode='seeded' requires a gen column.")
      }
      catalogue$gen == -1L
    },
    threshold = catalogue$magnitudes >= M_trigger,
    all = rep(TRUE, nrow(catalogue))
  )
  catalogue[keep, c("ts", "magnitudes"), drop = FALSE]
}

compute_loglinear_mct_at <- function(t_eval,
                                     trigger_ts,
                                     trigger_magnitudes,
                                     M0,
                                     G,
                                     H,
                                     include_current = FALSE,
                                     cap_at_trigger = FALSE) {
  if (G <= 0 || H <= 0) stop("G and H must be positive.")
  if (length(trigger_ts) != length(trigger_magnitudes)) {
    stop("trigger_ts and trigger_magnitudes must have equal length.")
  }

  out <- rep(M0, length(t_eval))
  if (length(trigger_ts) == 0L || length(t_eval) == 0L) return(out)

  for (i in seq_along(t_eval)) {
    eligible <- if (include_current) {
      trigger_ts <= t_eval[i]
    } else {
      trigger_ts < t_eval[i]
    }
    if (!any(eligible)) next

    dt <- t_eval[i] - trigger_ts[eligible]
    positive <- dt > 0
    if (!any(positive)) next
    mags <- trigger_magnitudes[eligible][positive]
    candidates <- mags - G - H * log10(dt[positive])
    if (cap_at_trigger) candidates <- pmin(candidates, mags)
    out[i] <- max(M0, candidates)
  }
  out
}

loglinear_detectable_fraction <- function(t_eval,
                                          trigger_ts,
                                          trigger_magnitudes,
                                          M0,
                                          G,
                                          H,
                                          b = 1,
                                          include_current = FALSE,
                                          cap_at_trigger = FALSE) {
  if (b <= 0) stop("b must be positive.")
  Mc <- compute_loglinear_mct_at(
    t_eval, trigger_ts, trigger_magnitudes, M0, G, H,
    include_current = include_current,
    cap_at_trigger = cap_at_trigger
  )
  pmin(1, 10^(-b * pmax(Mc - M0, 0)))
}

apply_loglinear_incompleteness <- function(catalogue,
                                           M0,
                                           G = 4.5,
                                           H = 0.75,
                                           b = 1,
                                           trigger_mode = c("seeded", "threshold", "all"),
                                           M_trigger = 5.5,
                                           cap_at_trigger = FALSE,
                                           retain_seeded = TRUE,
                                           retain_triggers = FALSE) {
  trigger_mode <- match.arg(trigger_mode)
  required <- c("ts", "magnitudes")
  if (!all(required %in% names(catalogue))) {
    stop("catalogue must contain ts and magnitudes.")
  }

  order_cols <- if ("ID" %in% names(catalogue)) {
    order(catalogue$ts, catalogue$ID)
  } else {
    order(catalogue$ts)
  }
  catalogue <- catalogue[order_cols, , drop = FALSE]
  rownames(catalogue) <- NULL
  triggers <- select_loglinear_triggers(catalogue, trigger_mode, M_trigger)

  Mc <- compute_loglinear_mct_at(
    t_eval = catalogue$ts,
    trigger_ts = triggers$ts,
    trigger_magnitudes = triggers$magnitudes,
    M0 = M0,
    G = G,
    H = H,
    include_current = FALSE,
    cap_at_trigger = cap_at_trigger
  )
  detected <- catalogue$magnitudes >= Mc
  if (retain_seeded && "gen" %in% names(catalogue)) {
    detected[catalogue$gen == -1L] <- TRUE
  }
  if (retain_triggers) {
    trigger_key <- paste(triggers$ts, triggers$magnitudes)
    event_key <- paste(catalogue$ts, catalogue$magnitudes)
    detected[event_key %in% trigger_key] <- TRUE
  }

  audit <- catalogue
  audit$Mc_t <- Mc
  audit$Psi_t <- pmin(1, 10^(-b * pmax(Mc - M0, 0)))
  audit$detected <- detected
  observed <- audit[detected, , drop = FALSE]
  rownames(observed) <- NULL

  attr(observed, "detection_audit") <- audit
  attr(observed, "observation_model") <- "helmstetter_loglinear_hard_threshold"
  attr(observed, "incompleteness_parameters") <- list(
    M0 = M0, G = G, H = H, b = b,
    trigger_mode = trigger_mode,
    M_trigger = M_trigger,
    cap_at_trigger = cap_at_trigger,
    retain_triggers = retain_triggers
  )
  observed
}

loglinear_recovery_end <- function(trigger_time, trigger_magnitude, M0, G, H) {
  trigger_time + 10^((trigger_magnitude - G - M0) / H)
}

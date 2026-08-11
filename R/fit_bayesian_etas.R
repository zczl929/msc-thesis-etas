# Shared Bayesian ETAS fitting helpers using ETAS.inlabru/INLA.

library(ETAS.inlabru)
library(inlabru)
library(INLA)
library(foreach)

source("R/etas_likelihood.R")
source("R/loglinear_incompleteness.R")

make_inlabru_catalog <- function(catalog, T1, T2) {
  catalog <- catalog[order(catalog$ts), ]
  catalog <- catalog[catalog$ts >= T1 & catalog$ts <= T2, ]
  include_likelihood <- if ("include_likelihood" %in% names(catalog)) {
    as.logical(catalog$include_likelihood)
  } else if ("gen" %in% names(catalog)) {
    # Synthetic mainshocks with gen == -1 are externally seeded.  They remain
    # in the conditioning history and can trigger offspring, but their own
    # occurrence must not be treated as a stochastic ETAS event.
    catalog$gen != -1L
  } else {
    rep(TRUE, nrow(catalog))
  }

  list(
    catalog = catalog,
    catalog_bru = data.frame(
      ts = catalog$ts,
      magnitudes = catalog$magnitudes,
      idx.p = seq_len(nrow(catalog)),
      include_likelihood = include_likelihood
    )
  )
}

default_etas_prior_spec <- function() {
  list(
    a_mu = 2,
    b_mu = 4,
    a_K = log(0.5),
    b_K = 1.0,
    a_alpha = 0.05,
    b_alpha = 2.0,
    a_c = 0.001,
    b_c = 0.2,
    a_p = 1.01,
    b_p = 2.5,
    mu_init = 0.5,
    K_init = 0.5,
    alpha_init = 0.8,
    c_init = 0.02,
    p_init = 1.3
  )
}

make_etas_link_functions <- function(prior_spec) {
  list(
    mu = function(x) gamma_t(x, prior_spec$a_mu, prior_spec$b_mu),
    K = function(x) loggaus_t(x, prior_spec$a_K, prior_spec$b_K),
    alpha = function(x) unif_t(x, prior_spec$a_alpha, prior_spec$b_alpha),
    c_ = function(x) unif_t(x, prior_spec$a_c, prior_spec$b_c),
    p = function(x) unif_t(x, prior_spec$a_p, prior_spec$b_p)
  )
}

make_bru_options <- function(prior_spec, verbose = 3, max_iter = 20,
                             rel_tol = 0.1) {
  if (!is.finite(rel_tol) || rel_tol <= 0 || rel_tol >= 1) {
    stop("rel_tol must lie strictly between zero and one.")
  }
  initial_values <- list(
    th.mu = inv_gamma_t(prior_spec$mu_init, prior_spec$a_mu, prior_spec$b_mu),
    th.K = inv_loggaus_t(prior_spec$K_init, prior_spec$a_K, prior_spec$b_K),
    th.alpha = inv_unif_t(prior_spec$alpha_init, prior_spec$a_alpha, prior_spec$b_alpha),
    th.c = inv_unif_t(prior_spec$c_init, prior_spec$a_c, prior_spec$b_c),
    th.p = inv_unif_t(prior_spec$p_init, prior_spec$a_p, prior_spec$b_p)
  )

  list(
    bru_verbose = verbose,
    bru_max_iter = max_iter,
    bru_initial = initial_values,
    bru_method = list(rel_tol = rel_tol),
    num.threads = "1:1"
  )
}

apply_resume_fit <- function(bru_options, resume_path = NULL) {
  if (is.null(resume_path) || !nzchar(resume_path) || !file.exists(resume_path)) {
    return(bru_options)
  }
  previous <- readRDS(resume_path)
  if (is.null(previous$fit)) stop("Resume file has no fit object: ", resume_path)
  states <- previous$fit$bru_iinla$states
  if (length(states) == 0) stop("Resume file has no nonlinear states: ", resume_path)
  bru_options$bru_initial <- tail(states, 1)[[1]]
  cat("Resuming nonlinear iteration from:", resume_path, "\n")
  bru_options
}

summarise_posterior_draws <- function(posterior_draws, true_values = NULL) {
  posterior_summary <- t(vapply(
    names(posterior_draws),
    function(param_name) {
      x <- posterior_draws[[param_name]]
      out <- c(
        mean = mean(x),
        median = median(x),
        q025 = unname(quantile(x, 0.025)),
        q975 = unname(quantile(x, 0.975))
      )
      if (!is.null(true_values)) {
        out <- c(out, true = unname(true_values[param_name]))
      }
      out
    },
    numeric(if (is.null(true_values)) 4 else 5)
  ))

  rownames(posterior_summary) <- names(posterior_draws)
  posterior_summary
}

compact_bayesian_result <- function(result) {
  keep <- c(
    "method", "model_version", "experiment_id", "run_config", "convergence", "note",
    "input_path", "oracle_history_path", "oracle_detection_path",
    "input_md5", "true_values", "posterior_summary", "posterior_draws",
    "posterior_param"
  )
  result[intersect(keep, names(result))]
}

save_bayesian_fit <- function(result, output_dir, output_stem,
                              compact = Sys.getenv("COMPACT_OUTPUT", unset = "0") == "1") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(output_dir, paste0(output_stem, ".rds"))
  csv_path <- file.path(output_dir, paste0(output_stem, "_summary.csv"))
  rds_tmp <- tempfile(paste0(output_stem, "_"), tmpdir = output_dir)
  csv_tmp <- tempfile(paste0(output_stem, "_summary_"), tmpdir = output_dir)
  on.exit(unlink(c(rds_tmp, csv_tmp)), add = TRUE)

  saved_result <- if (compact) compact_bayesian_result(result) else result
  saveRDS(saved_result, rds_tmp)
  write.csv(result$posterior_summary, csv_tmp, row.names = TRUE)
  if (!file.rename(rds_tmp, rds_path)) {
    stop("Could not atomically install fit result: ", rds_path)
  }
  if (!file.rename(csv_tmp, csv_path)) {
    stop("Could not atomically install posterior summary: ", csv_path)
  }

  list(rds_path = rds_path, csv_path = csv_path)
}

summarise_bru_convergence <- function(fit, rel_tol = 0.1) {
  messages <- fit$bru_iinla$log$log$message
  deviation_messages <- messages[grepl("Max deviation from previous", messages)]
  last_deviation <- NA_real_
  if (length(deviation_messages) > 0) {
    last_deviation <- suppressWarnings(as.numeric(sub(
      ".*previous: ([0-9.]+)%.*",
      "\\1",
      tail(deviation_messages, 1)
    )))
  }
  max_reached <- any(grepl("Maximum iterations reached", messages))
  tolerance_percent_sd <- 100 * rel_tol
  list(
    converged = !max_reached && is.finite(last_deviation) &&
      last_deviation < tolerance_percent_sd,
    max_iterations_reached = max_reached,
    iterations = length(fit$bru_iinla$states),
    last_deviation_percent_sd = last_deviation,
    tolerance_percent_sd = tolerance_percent_sd
  )
}

fit_inlabru_baseline <- function(config) {
  prior_spec <- config$prior_spec
  link_functions <- make_etas_link_functions(prior_spec)
  bru_options <- make_bru_options(
    prior_spec,
    verbose = config$bru_verbose,
    max_iter = config$bru_max_iter,
    rel_tol = config$bru_rel_tol
  )
  bru_options <- apply_resume_fit(bru_options, config$resume_path)

  catalog_raw <- readRDS(config$input_path)
  catalog_data <- make_inlabru_catalog(catalog_raw, config$T1, config$T2)

  cat("Number of conditioning history events:", sum(catalog_data$catalog_bru$ts <= config$T1), "\n")
  cat("Number of target events:", sum(catalog_data$catalog_bru$ts > config$T1 & catalog_data$catalog_bru$ts < config$T2), "\n")

  cat("Starting canonical baseline ETAS/INLA fit...\n")
  fit <- Temporal.ETAS.Canonical(
    total.data = catalog_data$catalog_bru,
    M0 = config$M0,
    T1 = config$T1,
    T2 = config$T2,
    link.functions = link_functions,
    coef.t. = config$coef_t,
    delta.t. = config$delta_t,
    N.max. = config$N_max,
    bru.opt = bru_options,
    use_detection = FALSE,
    M_trigger = Inf,
    gamma = 0,
    tau = 1,
    mct_mode = "max",
    b_GR = 1
  )
  cat("Finished canonical baseline ETAS/INLA fit.\n")

  input_list <- list(
    catalog = catalog_data$catalog,
    catalog.bru = catalog_data$catalog_bru,
    model.fit = fit,
    T12 = c(config$T1, config$T2),
    M0 = config$M0,
    link.functions = link_functions,
    prior_spec = prior_spec,
    bru.opt.list = bru_options,
    coef.t = config$coef_t,
    delta.t = config$delta_t,
    Nmax = config$N_max
  )

  posterior_param <- get_posterior_param(input_list)
  posterior_draws <- post_sampling(input_list, n.samp = config$n_draws, max.batch = config$max_batch)
  posterior_summary <- summarise_posterior_draws(posterior_draws, config$true_values)

  result <- list(
    method = "Canonical normalised temporal ETAS/INLA baseline",
    model_version = "canonical-base10-normalised-v1",
    experiment_id = config$experiment_id,
    run_config = list(
      n_draws = config$n_draws,
      max_batch = config$max_batch,
      bru_max_iter = config$bru_max_iter,
      bru_rel_tol = config$bru_rel_tol,
      coef_t = config$coef_t,
      delta_t = config$delta_t,
      N_max = config$N_max,
      prior_spec = config$prior_spec
    ),
    convergence = summarise_bru_convergence(fit, config$bru_rel_tol),
    note = paste(
      "Same base-10 productivity and normalised Omori kernel as the simulator",
      "and TVC fit, with detection probability fixed to one."
    ),
    input_path = config$input_path,
    input_md5 = unname(tools::md5sum(config$input_path)),
    true_values = config$true_values,
    posterior_summary = posterior_summary,
    posterior_draws = posterior_draws,
    posterior_param = posterior_param$post.df,
    fit = fit,
    input_list = input_list
  )

  paths <- save_bayesian_fit(result, config$output_dir, config$output_stem)
  list(result = result, paths = paths)
}

Temporal.ETAS.Canonical <- function(total.data,
                              history.data = NULL,
                              detection.data = NULL,
                              M0,
                              T1,
                              T2,
                              link.functions,
                              coef.t.,
                              delta.t.,
                              N.max.,
                              bru.opt,
                              use_detection = TRUE,
                              detection_model = "exponential_recovery",
                              M_trigger,
                              gamma,
                              tau,
                              mct_mode,
                              b_GR = 1,
                              loglinear_G = NA_real_,
                              loglinear_H = NA_real_,
                              mainshock_time = NA_real_,
                              mainshock_magnitude = NA_real_) {
  compute_psi_at <- function(t_eval, event_ts, event_magnitudes,
                             include_current = TRUE) {
    if (!use_detection) return(rep(1, length(t_eval)))

    Mc_eval <- if (identical(detection_model, "loglinear_mainshock")) {
      if (!all(is.finite(c(
        loglinear_G, loglinear_H, mainshock_time, mainshock_magnitude
      )))) {
        stop("Finite G, H, mainshock_time and mainshock_magnitude are required.")
      }
      compute_loglinear_mct_at(
        t_eval = t_eval,
        trigger_ts = mainshock_time,
        trigger_magnitudes = mainshock_magnitude,
        M0 = M0,
        G = loglinear_G,
        H = loglinear_H,
        include_current = include_current
      )
    } else if (identical(detection_model, "loglinear_multitrigger")) {
      if (!all(is.finite(c(loglinear_G, loglinear_H, M_trigger)))) {
        stop("Finite G, H and M_trigger are required.")
      }
      trigger <- event_magnitudes >= M_trigger
      compute_loglinear_mct_at(
        t_eval = t_eval,
        trigger_ts = event_ts[trigger],
        trigger_magnitudes = event_magnitudes[trigger],
        M0 = M0, G = loglinear_G, H = loglinear_H,
        include_current = include_current
      )
    } else if (identical(detection_model, "exponential_recovery")) {
      compute_mct_at(
        t_eval = t_eval,
        event_ts = event_ts,
        event_magnitudes = event_magnitudes,
        M_cut = M0,
        M_trigger = M_trigger,
        gamma = gamma,
        tau = tau,
        mode = mct_mode,
        include_current = include_current
      )
    } else {
      stop("Unknown detection_model: ", detection_model)
    }

    pmin(1, 10^(-b_GR * (Mc_eval - M0)))
  }

  clean_catalog <- function(x, label) {
    if (sum(is.na(x)) > 0) {
      cat("Some", label, "data are NA; removing rows\n")
      x <- na.omit(x)
    }
    idx.after <- x$ts > T2
    if (sum(idx.after) > 0) {
      x <- x[!idx.after, ]
      warning("Removing ", label, " events after T2")
    }
    x[order(x$ts), , drop = FALSE]
  }

  target.data <- clean_catalog(total.data, "target")
  if (is.null(history.data)) history.data <- target.data
  if (is.null(detection.data)) detection.data <- target.data
  history.data <- clean_catalog(history.data, "history")
  detection.data <- clean_catalog(detection.data, "detection")

  if (!"include_likelihood" %in% names(target.data)) {
    target.data$include_likelihood <- TRUE
  }
  idx.sample <- target.data$ts > T1 & target.data$ts < T2 &
    target.data$include_likelihood
  sample.s <- target.data[idx.sample, ]

  sample.s$psi_t <- compute_psi_at(
    t_eval = sample.s$ts,
    event_ts = detection.data$ts,
    event_magnitudes = detection.data$magnitudes,
    include_current = FALSE
  )

  psi_grid <- seq(T1, T2, length.out = 1000)
  psi_vals <- compute_psi_at(
    t_eval = psi_grid,
    event_ts = detection.data$ts,
    event_magnitudes = detection.data$magnitudes
  )
  integral_psi <- sum(0.5 * (psi_vals[-1] + psi_vals[-length(psi_vals)]) * diff(psi_grid))

  df.0 <- data.frame(counts = 0, exposures = 1, part = "background")

  cat("Start creating TVC grid...\n")
  time.g.st <- Sys.time()
  df.j <- foreach(idx = seq_len(nrow(history.data)), .combine = rbind) %do% {
    time_grid(
      data.point = history.data[idx, ],
      coef.t = coef.t.,
      delta.t = delta.t.,
      T1. = T1,
      T2. = T2,
      N.exp. = N.max.
    )
  }

  df.j$counts <- 0
  df.j$exposures <- 1
  df.j$part <- "triggered"

  t.lower <- pmax(df.j$t.start, df.j$ts)
  t.mid <- 0.5 * (t.lower + df.j$t.end)
  df.j$psi_mid <- compute_psi_at(
    t_eval = t.mid,
    event_ts = detection.data$ts,
    event_magnitudes = detection.data$magnitudes
  )

  cat("Finished creating TVC grid, time", Sys.time() - time.g.st, "\n")

  It_df_tvc <- function(param_, time.df) {
    tth <- as.numeric(time.df$ts)
    T1b <- as.numeric(time.df$t.start)
    T2b <- as.numeric(time.df$t.end)
    param_c <- param_[4]
    param_p <- param_[5]
    T.l <- pmax(tth, T1b)
    fun.l <- (1 + (T.l - tth) / param_c)^(1 - param_p)
    fun.u <- (1 + (T2b - tth) / param_c)^(1 - param_p)
    base_integral <- fun.l - fun.u
    base_integral * time.df$psi_mid
  }

  cond_lambda_mle_param <- function(theta, t, th, mh, M0) {
    if (is.null(th) || length(th) == 0 || all(th >= t)) {
      return(theta$mu)
    }

    past <- th < t
    dt <- t - th[past]
    productivity <- theta$K * 10^(theta$alpha * (mh[past] - M0))
    decay <- ((theta$p - 1) / theta$c) * (1 + dt / theta$c)^(-theta$p)
    theta$mu + sum(productivity * decay)
  }

  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, list.input_) {
    theta_ <- c(
      0,
      link.functions$K(th.K[1]),
      link.functions$alpha(th.alpha[1]),
      link.functions$c_(th.c[1]),
      link.functions$p(th.p[1])
    )

    comp. <- It_df_tvc(param_ = theta_, time.df = list.input_$df_grid)
    log(10) * theta_[3] * (list.input_$df_grid$magnitudes - list.input_$M0) +
      log(theta_[2] + 1e-100) +
      log(comp. + 1e-100)
  }

  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = "SL")

  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, psi_event, M0) {
    th.p <- list(
      mu = link.functions$mu(th.mu[1]),
      K = link.functions$K(th.K[1]),
      alpha = link.functions$alpha(th.alpha[1]),
      c = link.functions$c_(th.c[1]),
      p = link.functions$p(th.p[1])
    )

    mean(vapply(
      seq_along(tt),
      function(i) {
        x <- tt[i]
        th_x <- th < x
        log(cond_lambda_mle_param(theta = th.p, t = x, th = th[th_x], mh = mh[th_x], M0 = M0)) +
          log(pmax(psi_event[i], 1e-100))
      },
      numeric(1)
    ))
  }

  list.input <- list(
    df_grid = df.j,
    M0 = M0,
    sample.s = sample.s,
    target.data = target.data,
    history.data = history.data,
    detection.data = detection.data,
    integral_psi = integral_psi
  )

  data.input <- dplyr::bind_rows(df.0, df.s, df.j)
  list.input <- append(
    list.input,
    list(
      idx.bkg = data.input$part == "background",
      idx.trig = data.input$part == "triggered",
      idx.sl = data.input$part == "SL"
    )
  )

  predictor.fun <- function(th.mu, th.K, th.alpha, th.c, th.p, list.input, T1, T2, M0) {
    out <- rep(0, nrow(data.input))
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) +
      log(list.input$integral_psi + 1e-100)
    out[list.input$idx.trig] <- logLambda.h.inla(
      th.K = th.K,
      th.alpha = th.alpha,
      th.c = th.c,
      th.p = th.p,
      list.input_ = list.input
    )
    out[list.input$idx.sl] <- loglambda.inla(
      th.mu = th.mu,
      th.K = th.K,
      th.alpha = th.alpha,
      th.c = th.c,
      th.p = th.p,
      tt = list.input$sample.s$ts,
      th = list.input$history.data$ts,
      mh = list.input$history.data$magnitudes,
      psi_event = list.input$sample.s$psi_t,
      M0 = M0
    )
    out
  }

  merged.form <- counts ~ predictor.fun(
    th.mu = th.mu,
    th.K = th.K,
    th.alpha = th.alpha,
    th.c = th.c,
    th.p = th.p,
    list.input = list.input,
    T1 = T1,
    T2 = T2,
    M0 = M0
  )

  cmp.part <- counts ~ -1 +
    th.mu(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.K(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.c(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.p(1, model = "linear", mean.linear = 0, prec.linear = 1)

  inlabru::bru(
    components = cmp.part,
    inlabru::bru_obs(
      formula = merged.form,
      data = data.input,
      family = "poisson",
      E = data.input$exposures
    ),
    options = bru.opt
  )
}

fit_inlabru_tvc <- function(config) {
  prior_spec <- config$prior_spec
  link_functions <- make_etas_link_functions(prior_spec)
  bru_options <- make_bru_options(
    prior_spec,
    verbose = config$bru_verbose,
    max_iter = config$bru_max_iter,
    rel_tol = config$bru_rel_tol
  )
  bru_options <- apply_resume_fit(bru_options, config$resume_path)

  catalog_raw <- readRDS(config$input_path)
  catalog_data <- make_inlabru_catalog(catalog_raw, config$T1, config$T2)
  history_data <- NULL
  detection_data <- NULL
  if (!is.null(config$oracle_history_path) && nzchar(config$oracle_history_path)) {
    history_raw <- readRDS(config$oracle_history_path)
    history_data <- make_inlabru_catalog(history_raw, config$T1, config$T2)$catalog_bru
  }
  if (!is.null(config$oracle_detection_path) && nzchar(config$oracle_detection_path)) {
    detection_raw <- readRDS(config$oracle_detection_path)
    detection_data <- make_inlabru_catalog(detection_raw, config$T1, config$T2)$catalog_bru
  }

  cat("Number of conditioning history events:", sum(catalog_data$catalog_bru$ts <= config$T1), "\n")
  cat("Number of target events:", sum(catalog_data$catalog_bru$ts > config$T1 & catalog_data$catalog_bru$ts < config$T2), "\n")
  cat("Oracle history:", !is.null(history_data), "| oracle detection:", !is.null(detection_data), "\n")
  cat("TVC mode:", config$tvc_settings$mct_mode, "\n")

  cat("Starting TVC ETAS.inlabru/INLA fit...\n")
  fit <- Temporal.ETAS.Canonical(
    total.data = catalog_data$catalog_bru,
    history.data = history_data,
    detection.data = detection_data,
    M0 = config$M0,
    T1 = config$T1,
    T2 = config$T2,
    link.functions = link_functions,
    coef.t. = config$coef_t,
    delta.t. = config$delta_t,
    N.max. = config$N_max,
    bru.opt = bru_options,
    use_detection = TRUE,
    detection_model = if (is.null(config$tvc_settings$detection_model)) {
      "exponential_recovery"
    } else {
      config$tvc_settings$detection_model
    },
    M_trigger = config$tvc_settings$M_trigger,
    gamma = config$tvc_settings$gamma,
    tau = config$tvc_settings$tau,
    mct_mode = config$tvc_settings$mct_mode,
    b_GR = config$tvc_settings$b_GR,
    loglinear_G = if (is.null(config$tvc_settings$G)) NA_real_ else config$tvc_settings$G,
    loglinear_H = if (is.null(config$tvc_settings$H)) NA_real_ else config$tvc_settings$H,
    mainshock_time = if (is.null(config$tvc_settings$mainshock_time)) {
      NA_real_
    } else {
      config$tvc_settings$mainshock_time
    },
    mainshock_magnitude = if (is.null(config$tvc_settings$mainshock_magnitude)) {
      NA_real_
    } else {
      config$tvc_settings$mainshock_magnitude
    }
  )
  cat("Finished TVC ETAS.inlabru/INLA fit.\n")

  input_list <- list(
    catalog = catalog_data$catalog,
    catalog.bru = catalog_data$catalog_bru,
    model.fit = fit,
    T12 = c(config$T1, config$T2),
    M0 = config$M0,
    link.functions = link_functions,
    prior_spec = prior_spec,
    tvc_settings = config$tvc_settings,
    bru.opt.list = bru_options,
    coef.t = config$coef_t,
    delta.t = config$delta_t,
    Nmax = config$N_max
  )

  posterior_param <- get_posterior_param(input_list)
  posterior_draws <- post_sampling(input_list, n.samp = config$n_draws, max.batch = config$max_batch)
  posterior_summary <- summarise_posterior_draws(posterior_draws, config$true_values)

  result <- list(
    method = if (is.null(history_data) && is.null(detection_data)) {
      "Canonical normalised temporal TVC-ETAS/INLA"
    } else {
      "Simulation-only oracle canonical temporal TVC-ETAS/INLA"
    },
    model_version = "canonical-base10-normalised-v1",
    experiment_id = config$experiment_id,
    run_config = list(
      n_draws = config$n_draws,
      max_batch = config$max_batch,
      bru_max_iter = config$bru_max_iter,
      bru_rel_tol = config$bru_rel_tol,
      coef_t = config$coef_t,
      delta_t = config$delta_t,
      N_max = config$N_max,
      prior_spec = config$prior_spec,
      tvc_settings = config$tvc_settings
    ),
    convergence = summarise_bru_convergence(fit, config$bru_rel_tol),
    note = "Custom Temporal.ETAS.TVC with lambda_obs(t) = lambda_ETAS(t) * psi(Mc(t)).",
    input_path = config$input_path,
    input_md5 = unname(tools::md5sum(config$input_path)),
    oracle_history_path = config$oracle_history_path,
    oracle_detection_path = config$oracle_detection_path,
    true_values = config$true_values,
    posterior_summary = posterior_summary,
    posterior_draws = posterior_draws,
    posterior_param = posterior_param$post.df,
    fit = fit,
    input_list = input_list
  )

  paths <- save_bayesian_fit(result, config$output_dir, config$output_stem)
  list(result = result, paths = paths)
}

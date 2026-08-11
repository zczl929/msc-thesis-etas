# Shared orchestration and integrity checks for scripts/submission/.

source("config/submission_experiment.R")

submission_rscript <- function() {
  file.path(R.home("bin"), "Rscript")
}

submission_fit_expected <- function(cfg, input_path, model,
                                    tvc_settings = NULL) {
  list(
    experiment_id = cfg$experiment_id,
    input_md5 = unname(tools::md5sum(input_path)),
    model = model,
    n_draws = cfg$fit$n_draws,
    rel_tol = cfg$fit$rel_tol,
    prior_spec = cfg$primary_prior,
    fit_spec = cfg$fit,
    tvc_settings = tvc_settings
  )
}

submission_fit_is_current <- function(path, expected) {
  if (!file.exists(path)) return(FALSE)
  fit <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(fit)) return(FALSE)
  run <- fit$run_config
  prior <- run$prior_spec
  method_matches <- switch(
    expected$model,
    naive = identical(
      fit$method,
      "Canonical normalised temporal ETAS/INLA baseline"
    ),
    plugin = identical(
      fit$method,
      "Canonical normalised temporal TVC-ETAS/INLA"
    ),
    FALSE
  )
  prior_names <- names(expected$prior_spec)
  prior_matches <- all(vapply(
    prior_names,
    function(name) {
      !is.null(prior[[name]]) &&
        isTRUE(all.equal(prior[[name]], expected$prior_spec[[name]]))
    },
    logical(1)
  ))
  fit_spec_matches <- all(c(
    isTRUE(all.equal(run$max_batch, expected$fit_spec$max_batch)),
    isTRUE(all.equal(run$bru_max_iter, expected$fit_spec$max_iter)),
    isTRUE(all.equal(run$coef_t, expected$fit_spec$coef_t)),
    isTRUE(all.equal(run$delta_t, expected$fit_spec$delta_t)),
    isTRUE(all.equal(run$N_max, expected$fit_spec$N_max))
  ))
  tvc_matches <- if (is.null(expected$tvc_settings)) {
    TRUE
  } else {
    saved <- run$tvc_settings
    all(vapply(
      names(expected$tvc_settings),
      function(name) {
        !is.null(saved[[name]]) &&
          isTRUE(all.equal(saved[[name]], expected$tvc_settings[[name]]))
      },
      logical(1)
    ))
  }
  required_parameters <- c("mu", "K", "alpha", "c", "p")
  draws <- fit$posterior_draws
  draws_valid <- is.data.frame(draws) &&
    all(required_parameters %in% names(draws)) &&
    nrow(draws) >= expected$n_draws &&
    all(vapply(
      draws[required_parameters],
      function(x) is.numeric(x) && all(is.finite(x)),
      logical(1)
    )) &&
    all(draws$mu > 0) &&
    all(draws$K > 0) &&
    all(draws$alpha >= expected$prior_spec$a_alpha) &&
    all(draws$alpha <= expected$prior_spec$b_alpha) &&
    all(draws$c >= expected$prior_spec$a_c) &&
    all(draws$c <= expected$prior_spec$b_c) &&
    all(draws$p >= expected$prior_spec$a_p) &&
    all(draws$p <= expected$prior_spec$b_p)
  checks <- c(
    method_matches,
    identical(fit$experiment_id, expected$experiment_id),
    identical(fit$input_md5, expected$input_md5),
    isTRUE(run$n_draws >= expected$n_draws),
    isTRUE(all.equal(run$bru_rel_tol, expected$rel_tol)),
    prior_matches,
    fit_spec_matches,
    tvc_matches,
    draws_valid,
    isTRUE(fit$convergence$converged)
  )
  all(checks)
}

submission_tvc_expected_from_env <- function(extra_env) {
  get_value <- function(name, numeric = TRUE) {
    prefix <- paste0(name, "=")
    value <- sub(prefix, "", extra_env[startsWith(extra_env, prefix)][1L])
    if (!length(value) || is.na(value)) return(NULL)
    if (numeric) as.numeric(value) else value
  }
  out <- list(
    detection_model = get_value("DETECTION_MODEL", numeric = FALSE),
    M_trigger = get_value("M_TRIGGER"),
    b_GR = get_value("B_GR"),
    G = get_value("LOGLINEAR_G"),
    H = get_value("LOGLINEAR_H"),
    mainshock_time = get_value("MAINSHOCK_TIME"),
    mainshock_magnitude = get_value("MAINSHOCK_MAGNITUDE")
  )
  out[!vapply(out, is.null, logical(1))]
}

submission_common_env <- function(cfg, output_dir, prior = cfg$primary_prior,
                                  n_draws = cfg$fit$n_draws) {
  truth <- cfg$theta_true
  env <- c(
    paste0("EXPERIMENT_ID=", cfg$experiment_id),
    paste0("OUTPUT_DIR=", output_dir),
    paste0("M0=", cfg$M0),
    paste0("T1=", cfg$T1),
    paste0("T2=", cfg$T2),
    paste0("N_DRAWS=", n_draws),
    paste0("MAX_BATCH=", cfg$fit$max_batch),
    paste0("BRU_MAX_ITER=", cfg$fit$max_iter),
    paste0("BRU_REL_TOL=", cfg$fit$rel_tol),
    "FUTURE_GLOBALS_MAX_GB=2",
    "BRU_VERBOSE=1",
    "COMPACT_OUTPUT=1",
    paste0("N_MAX=", cfg$fit$N_max),
    paste0("COEF_T=", cfg$fit$coef_t),
    paste0("DELTA_T=", cfg$fit$delta_t),
    paste0("MU_PRIOR_SHAPE=", prior$a_mu),
    paste0("MU_PRIOR_RATE=", prior$b_mu),
    paste0("MU_INIT=", prior$mu_init),
    paste0("K_PRIOR_LOGMEAN=", prior$a_K),
    paste0("K_PRIOR_LOGSD=", prior$b_K),
    paste0("K_INIT=", prior$K_init),
    paste0("ALPHA_PRIOR_LOWER=", prior$a_alpha),
    paste0("ALPHA_PRIOR_UPPER=", prior$b_alpha),
    paste0("ALPHA_INIT=", prior$alpha_init),
    paste0("C_PRIOR_LOWER=", prior$a_c),
    paste0("C_PRIOR_UPPER=", prior$b_c),
    paste0("C_INIT=", prior$c_init),
    paste0("P_PRIOR_LOWER=", prior$a_p),
    paste0("P_PRIOR_UPPER=", prior$b_p),
    paste0("P_INIT=", prior$p_init)
  )
  if (is.null(truth)) {
    c(env, "REAL_DATA=1")
  } else {
    c(
      env,
      paste0("TRUE_MU=", truth[["mu"]]),
      paste0("TRUE_K=", truth[["K"]]),
      paste0("TRUE_ALPHA=", truth[["alpha"]]),
      paste0("TRUE_C=", truth[["c"]]),
      paste0("TRUE_P=", truth[["p"]])
    )
  }
}

run_submission_fit <- function(cfg, model, input_path, output_dir,
                               output_stem, extra_env = character(),
                               force = FALSE) {
  output_path <- file.path(output_dir, paste0(output_stem, ".rds"))
  expected <- submission_fit_expected(
    cfg,
    input_path,
    model,
    tvc_settings = if (model == "plugin") {
      submission_tvc_expected_from_env(extra_env)
    } else {
      NULL
    }
  )
  if (!force && submission_fit_is_current(output_path, expected)) {
    cat("Verified and reused:", output_path, "\n")
    return(output_path)
  }
  script <- if (model == "naive") {
    "scripts/submission/internal/fit_baseline_initialisation.R"
  } else {
    "scripts/submission/internal/fit_plugin_initialisation.R"
  }
  status <- system2(
    submission_rscript(),
    c("--vanilla", script),
    env = c(
      paste0("INPUT_PATH=", input_path),
      paste0("OUTPUT_STEM=", output_stem),
      submission_common_env(cfg, output_dir),
      extra_env
    )
  )
  if (!identical(status, 0L)) {
    # Some INLA/future backends can return a non-zero shutdown status after a
    # complete compact result has already been written.  Accept only if the
    # saved object passes every frozen metadata, checksum, draw-count and
    # convergence check; otherwise fail hard.
    if (submission_fit_is_current(output_path, expected)) {
      warning(
        "Backend returned status ", status,
        " after writing a fully validated fit: ", output_stem
      )
      return(output_path)
    }
    stop("Fit failed and no valid result was produced: ", output_stem)
  }
  if (!submission_fit_is_current(output_path, expected)) {
    stop("Fit completed but failed integrity checks: ", output_path)
  }
  output_path
}

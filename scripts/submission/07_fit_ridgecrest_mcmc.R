# Final four-chain MCMC fits for the frozen Ridgecrest Naive and Plug-in models.

source("config/submission_ridgecrest.R")
source("R/ridgecrest_mcmc.R")
source("R/fit_loglinear_mle.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_ridgecrest_config()
train <- readRDS(cfg$train_path)
stopifnot(identical(attr(train, "experiment_id"), cfg$experiment_id))
condition_mainshock <- identical(
  tolower(Sys.getenv("CONDITION_MAINSHOCK", unset = "false")), "true"
)
train$include_likelihood <- TRUE
completeness_mode <- tolower(Sys.getenv(
  "COMPLETENESS_MODE", unset = "single"
))
if (!completeness_mode %in% c("single", "multi")) {
  stop("COMPLETENESS_MODE must be 'single' or 'multi'.")
}
is_mainshock <- train$ID == cfg$mainshock_id
if (sum(is_mainshock) != 1L) {
  stop("Configured Ridgecrest mainshock is not unique in the training data.")
}
mainshock_time <- train$ts[is_mainshock]
mainshock_magnitude <- train$magnitudes[is_mainshock]
if (condition_mainshock) {
  train$include_likelihood[is_mainshock] <- FALSE
}
cfg$T1 <- min(train$ts) - 1e-6
cfg$T2 <- cfg$train_end
input_md5 <- unname(tools::md5sum(cfg$train_path))

gh_path <- file.path(
  cfg$results_dir,
  if (completeness_mode == "single") "plugin_GH_single.rds" else "plugin_GH.rds"
)
if (completeness_mode == "single") {
  gh_object <- estimate_loglinear_plugin(
    train, mainshock_time, mainshock_magnitude, cfg$M0, cfg$b
  )
  attr(gh_object, "experiment_id") <- cfg$experiment_id
  attr(gh_object, "input_md5") <- input_md5
  attr(gh_object, "completeness_mode") <- completeness_mode
  saveRDS(gh_object, gh_path)
} else {
  gh_object <- readRDS(gh_path)
  stopifnot(
    identical(attr(gh_object, "experiment_id"), cfg$experiment_id),
    identical(attr(gh_object, "input_md5"), input_md5)
  )
}
gh <- gh_object$estimate

models <- strsplit(
  Sys.getenv("MCMC_MODELS", unset = "naive,plugin"), ",", fixed = TRUE
)[[1]]
if (!all(models %in% c("naive", "plugin"))) {
  stop("MCMC_MODELS must contain only naive and/or plugin.")
}
settings <- list(
  iterations = as.integer(Sys.getenv("MCMC_ITER", unset = "24000")),
  warmup = as.integer(Sys.getenv("MCMC_WARMUP", unset = "8000")),
  thin = as.integer(Sys.getenv("MCMC_THIN", unset = "4")),
  n_chains = as.integer(Sys.getenv("MCMC_CHAINS", unset = "4")),
  nodes_per_interval = as.integer(Sys.getenv("QUAD_NODES", unset = "16"))
)
force <- identical(tolower(Sys.getenv("MCMC_FORCE")), "true")
recalibrate <- identical(
  tolower(Sys.getenv("MCMC_RECALIBRATE")), "true"
)
analysis_name <- paste0(
  "mcmc",
  if (condition_mainshock) "_conditioned" else "",
  if (completeness_mode == "single") "_single" else ""
)
out <- file.path(cfg$results_dir, analysis_name)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

for (model in models) {
  output_path <- file.path(out, paste0(model, ".rds"))
  existing <- if (file.exists(output_path)) readRDS(output_path) else NULL
  if (file.exists(output_path) && !force) {
    reusable <- identical(existing$experiment_id, cfg$experiment_id) &&
      identical(existing$input_md5, input_md5) &&
      identical(existing$model, model) &&
      identical(existing$prior, cfg$primary_prior) &&
      identical(existing$condition_mainshock, condition_mainshock) &&
      identical(existing$completeness_mode, completeness_mode) &&
      existing$settings$iterations >= settings$iterations &&
      existing$settings$warmup >= settings$warmup &&
      identical(existing$settings$thin, settings$thin) &&
      identical(existing$settings$n_chains, settings$n_chains) &&
      identical(
        existing$settings$nodes_per_interval,
        settings$nodes_per_interval
      )
    if (reusable) {
      cat("Reusing Ridgecrest MCMC fit:", model, "\n")
      next
    }
  }

  use_exact_pilot <- recalibrate && !is.null(existing) &&
    identical(existing$experiment_id, cfg$experiment_id) &&
    identical(existing$input_md5, input_md5) &&
    identical(existing$model, model)
  draws <- if (use_exact_pilot) {
    as.data.frame(do.call(rbind, existing$chains))
  } else {
    approximation <- readRDS(
      file.path(
        cfg$results_dir, "initialisation", paste0(model, ".rds")
      )
    )
    approximation$posterior_draws
  }
  draws <- draws[c("mu", "K", "alpha", "c", "p")]
  transformed <- t(apply(
    draws, 1, etas_theta_to_z, prior = cfg$primary_prior
  ))
  proposal_cov <- cov(transformed)
  initial_index <- unique(round(seq(
    1, nrow(draws), length.out = settings$n_chains + 2L
  )))[seq_len(settings$n_chains) + 1L]
  initial <- lapply(initial_index, function(index) transformed[index, ])

  log_posterior <- if (model == "naive") {
    make_exact_etas_log_posterior(
      train, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior
    )
  } else if (completeness_mode == "single") {
    make_fixedtrigger_plugin_log_posterior(
      train, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior,
      gh[["G"]], gh[["H"]],
      mainshock_time, mainshock_magnitude, cfg$b,
      settings$nodes_per_interval
    )
  } else {
    make_multitrigger_plugin_log_posterior(
      train, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior,
      gh[["G"]], gh[["H"]], cfg$b, cfg$M_trigger,
      settings$nodes_per_interval
    )
  }

  run_one <- function(chain) {
    result <- run_exact_etas_chain(
      log_posterior, initial[[chain]], proposal_cov,
      settings$iterations, settings$warmup, settings$thin,
      seed = 960700L + (if (model == "plugin") 100L else 0L) + chain
    )
    values <- t(apply(
      result$z, 1, etas_z_to_theta, prior = cfg$primary_prior
    ))
    list(
      draws = coda::mcmc(values),
      acceptance = result$acceptance,
      final_scale = result$final_scale
    )
  }
  raw_chains <- if (.Platform$OS.type == "unix") {
    parallel::mclapply(
      seq_len(settings$n_chains), run_one,
      mc.cores = settings$n_chains
    )
  } else {
    lapply(seq_len(settings$n_chains), run_one)
  }

  saveRDS(
    list(
      experiment_id = cfg$experiment_id,
      model = model,
      input_md5 = input_md5,
      condition_mainshock = condition_mainshock,
      completeness_mode = completeness_mode,
      conditioned_event_id = if (condition_mainshock) {
        cfg$mainshock_id
      } else {
        NULL
      },
      n_history_events = nrow(train),
      n_likelihood_targets = sum(train$include_likelihood),
      prior = cfg$primary_prior,
      gh_estimate = if (model == "plugin") gh else NULL,
      completeness_trigger = if (model == "plugin" &&
        completeness_mode == "single") {
        c(time = mainshock_time, magnitude = mainshock_magnitude)
      } else {
        NULL
      },
      M_trigger = if (model == "plugin" &&
        completeness_mode == "multi") cfg$M_trigger else NULL,
      chains = lapply(raw_chains, `[[`, "draws"),
      acceptance = vapply(raw_chains, `[[`, numeric(1), "acceptance"),
      final_scale = vapply(raw_chains, `[[`, numeric(1), "final_scale"),
      settings = settings,
      seed_rule = "960700 + 100 * I(plugin) + chain",
      proposal_source = if (use_exact_pilot) {
        "previous exact-MCMC pilot covariance"
      } else {
        "INLA/Laplace approximation covariance"
      }
    ),
    output_path
  )
  cat("Finished Ridgecrest MCMC fit:", model, "\n")
}

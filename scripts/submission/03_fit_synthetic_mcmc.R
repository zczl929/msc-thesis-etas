# Final exact-posterior MCMC fits for the 100 paired synthetic catalogues.
#
# This script does not regenerate catalogues or re-estimate the primary INLA
# fits.  The latter are used only to initialise dispersed chains and construct
# a proposal covariance.  Every saved draw targets the exact Naive or Plug-in
# posterior implemented in R/mcmc_etas.R.
#
# Optional environment variables:
#   MCMC_IDS=1,2,3        catalogue IDs (takes precedence)
#   MCMC_START=1 MCMC_END=50
#   MCMC_MODELS=naive,plugin
#   MCMC_ITER=12000       total iterations per chain
#   MCMC_WARMUP=4000
#   MCMC_THIN=4
#   MCMC_CHAINS=4
#   QUAD_NODES=16
#   MCMC_FORCE=true       overwrite existing compatible fits

source("config/submission_experiment.R")
source("R/mcmc_etas.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
ids <- if (nzchar(Sys.getenv("MCMC_IDS"))) {
  as.integer(strsplit(Sys.getenv("MCMC_IDS"), ",", fixed = TRUE)[[1]])
} else {
  seq.int(
    as.integer(Sys.getenv("MCMC_START", unset = "1")),
    as.integer(Sys.getenv("MCMC_END", unset = as.character(cfg$n_sim)))
  )
}
models <- strsplit(
  Sys.getenv("MCMC_MODELS", unset = "naive,plugin"), ",", fixed = TRUE
)[[1]]
if (anyNA(ids) || !all(ids %in% seq_len(cfg$n_sim))) {
  stop("MCMC_IDS contains an invalid catalogue ID.")
}
if (!all(models %in% c("naive", "plugin"))) {
  stop("MCMC_MODELS must contain only naive and/or plugin.")
}

settings <- list(
  iterations = as.integer(Sys.getenv("MCMC_ITER", unset = "12000")),
  warmup = as.integer(Sys.getenv("MCMC_WARMUP", unset = "4000")),
  thin = as.integer(Sys.getenv("MCMC_THIN", unset = "4")),
  n_chains = as.integer(Sys.getenv("MCMC_CHAINS", unset = "4")),
  nodes_per_interval = as.integer(Sys.getenv("QUAD_NODES", unset = "16"))
)
if (settings$warmup >= settings$iterations || settings$thin < 1L ||
    settings$n_chains < 2L) {
  stop("Invalid MCMC settings.")
}
force <- identical(tolower(Sys.getenv("MCMC_FORCE")), "true")
out <- file.path(cfg$results_dir, "mcmc_primary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fit_is_reusable <- function(path, input_md5, model, settings) {
  if (!file.exists(path) || force) return(FALSE)
  fit <- tryCatch(readRDS(path), error = function(e) NULL)
  settings_compatible <- !is.null(fit) &&
    fit$settings$iterations >= settings$iterations &&
    fit$settings$warmup >= settings$warmup &&
    identical(fit$settings$thin, settings$thin) &&
    identical(fit$settings$n_chains, settings$n_chains) &&
    identical(
      fit$settings$nodes_per_interval, settings$nodes_per_interval
    )
  !is.null(fit) &&
    identical(fit$experiment_id, cfg$experiment_id) &&
    identical(fit$input_md5, input_md5) &&
    identical(fit$model, model) &&
    settings_compatible &&
    length(fit$chains) == settings$n_chains
}

import_legacy_fit <- function(path, input_md5, id, model, gh) {
  legacy_path <- file.path(
    cfg$results_dir, "diagnostics", "incomplete_inla_vs_mcmc",
    sprintf("mcmc_%s_%03d.rds", model, id)
  )
  if (!file.exists(legacy_path) || force) return(FALSE)
  legacy <- readRDS(legacy_path)
  if (!identical(legacy$input_md5, input_md5) ||
      !identical(legacy$model, model) ||
      !identical(legacy$settings, settings)) {
    return(FALSE)
  }
  saveRDS(
    list(
      experiment_id = cfg$experiment_id,
      sim_id = id,
      model = model,
      input_md5 = input_md5,
      prior = cfg$primary_prior,
      gh_estimate = gh,
      chains = legacy$chains,
      settings = settings,
      seed_rule = "730000 + 1000 * sim_id + 100 * I(plugin) + chain",
      imported_from = legacy_path
    ),
    path
  )
  TRUE
}

for (id in ids) {
  input_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", id))
  if (!file.exists(input_path)) stop("Missing catalogue: ", input_path)
  input_md5 <- unname(tools::md5sum(input_path))
  catalogue <- readRDS(input_path)

  for (model in models) {
    output_path <- file.path(out, sprintf("%s_%03d.rds", model, id))
    inla_path <- file.path(
      cfg$results_dir, "initialisation", sprintf("%s_%03d.rds", model, id)
    )
    if (!file.exists(inla_path)) stop("Missing primary fit: ", inla_path)
    inla <- readRDS(inla_path)
    stopifnot(
      identical(inla$input_md5, input_md5),
      isTRUE(inla$convergence$converged)
    )

    gh <- NULL
    if (model == "plugin") {
      gh_path <- file.path(
        cfg$results_dir, "initialisation", sprintf("plugin_GH_%03d.rds", id)
      )
      gh <- readRDS(gh_path)$estimate
    }
    if (fit_is_reusable(output_path, input_md5, model, settings)) {
      cat("Reusing MCMC fit:", model, id, "\n")
      next
    }
    if (import_legacy_fit(output_path, input_md5, id, model, gh)) {
      cat("Imported validated 20-catalogue MCMC fit:", model, id, "\n")
      next
    }

    inla_draws <- inla$posterior_draws[c("mu", "K", "alpha", "c", "p")]
    transformed <- t(apply(
      inla_draws, 1, etas_theta_to_z, prior = cfg$primary_prior
    ))
    proposal_cov <- cov(transformed)
    initial_index <- unique(round(seq(
      1, nrow(inla_draws), length.out = settings$n_chains + 2L
    )))[seq_len(settings$n_chains) + 1L]
    initial <- lapply(initial_index, function(index) transformed[index, ])

    log_posterior <- if (model == "naive") {
      make_exact_etas_log_posterior(
        catalogue, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior
      )
    } else {
      make_plugin_etas_log_posterior(
        catalogue, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior,
        cfg$mainshock$ts, cfg$mainshock$magnitudes,
        gh[["G"]], gh[["H"]], cfg$b, settings$nodes_per_interval
      )
    }

    run_one <- function(chain) {
      result <- run_exact_etas_chain(
        log_posterior, initial[[chain]], proposal_cov,
        settings$iterations, settings$warmup, settings$thin,
        seed = 730000L + 1000L * id +
          (if (model == "plugin") 100L else 0L) + chain
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
    chains <- lapply(raw_chains, `[[`, "draws")

    saveRDS(
      list(
        experiment_id = cfg$experiment_id,
        sim_id = id,
        model = model,
        input_md5 = input_md5,
        prior = cfg$primary_prior,
        gh_estimate = gh,
        chains = chains,
        acceptance = vapply(raw_chains, `[[`, numeric(1), "acceptance"),
        final_scale = vapply(raw_chains, `[[`, numeric(1), "final_scale"),
        settings = settings,
        seed_rule = "730000 + 1000 * sim_id + 100 * I(plugin) + chain",
        imported_from = NULL
      ),
      output_path
    )
    cat("Finished MCMC fit:", model, id, "\n")
  }
}

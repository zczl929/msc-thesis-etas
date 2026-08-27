# Alternative-prior sensitivity for the 100 paired synthetic catalogues.
# Both Naive and Plug-in models are refitted so that the principal comparison
# is evaluated under the same alternative prior.  Primary outputs are never
# overwritten.

source("config/submission_experiment.R")
source("R/fit_loglinear_mle.R")
source("R/mcmc_etas.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
prior <- cfg$sensitivity_prior
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
  stop("Invalid catalogue ID.")
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
    settings$n_chains < 2L) stop("Invalid MCMC settings.")
force <- identical(tolower(Sys.getenv("MCMC_FORCE")), "true")
out <- file.path(cfg$results_dir, "prior_sensitivity", "synthetic")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pilot <- lapply(c("naive", "plugin"), function(model) {
  path <- file.path(cfg$results_dir, "initialisation", paste0(model, "_001.rds"))
  fit <- readRDS(path)
  stopifnot(isTRUE(fit$convergence$converged))
  draws <- fit$posterior_draws[c("mu", "K", "alpha", "c", "p")]
  z <- t(apply(draws, 1, etas_theta_to_z, prior = prior))
  indices <- unique(round(seq(
    1, nrow(z), length.out = settings$n_chains + 2L
  )))[seq_len(settings$n_chains) + 1L]
  list(
    path = path,
    proposal_cov = cov(z),
    initial = lapply(indices, function(index) z[index, ])
  )
})
names(pilot) <- c("naive", "plugin")

reusable <- function(path, input_md5, model) {
  if (!file.exists(path) || force) return(FALSE)
  fit <- tryCatch(readRDS(path), error = function(e) NULL)
  isTRUE(!is.null(fit) &&
    identical(fit$experiment_id, cfg$experiment_id) &&
    identical(fit$analysis, "alternative-prior-sensitivity") &&
    identical(fit$model, model) &&
    identical(fit$input_md5, input_md5) &&
    identical(fit$prior, prior) &&
    fit$settings$iterations >= settings$iterations &&
    fit$settings$warmup >= settings$warmup &&
    identical(fit$settings$thin, settings$thin) &&
    identical(fit$settings$n_chains, settings$n_chains) &&
    identical(fit$settings$nodes_per_interval, settings$nodes_per_interval) &&
    length(fit$chains) == settings$n_chains)
}

for (id in ids) {
  input_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", id))
  input_md5 <- unname(tools::md5sum(input_path))
  catalogue <- readRDS(input_path)
  gh_fit <- estimate_loglinear_plugin(
    catalogue, cfg$mainshock$ts, cfg$mainshock$magnitudes,
    cfg$M0, b = cfg$b
  )
  for (model in models) {
    output_path <- file.path(out, sprintf("%s_%03d.rds", model, id))
    if (reusable(output_path, input_md5, model)) {
      cat("Reusing alternative-prior fit:", model, id, "\n")
      next
    }
    log_posterior <- if (model == "naive") {
      make_exact_etas_log_posterior(
        catalogue, cfg$T1, cfg$T2, cfg$M0, prior
      )
    } else {
      gh <- gh_fit$estimate
      make_plugin_etas_log_posterior(
        catalogue, cfg$T1, cfg$T2, cfg$M0, prior,
        cfg$mainshock$ts, cfg$mainshock$magnitudes,
        gh[["G"]], gh[["H"]], cfg$b, settings$nodes_per_interval
      )
    }
    run_one <- function(chain) {
      result <- run_exact_etas_chain(
        log_posterior, pilot[[model]]$initial[[chain]],
        pilot[[model]]$proposal_cov,
        settings$iterations, settings$warmup, settings$thin,
        seed = 980000L + 1000L * id +
          (if (model == "plugin") 100L else 0L) + chain
      )
      values <- t(apply(result$z, 1, etas_z_to_theta, prior = prior))
      list(
        draws = coda::mcmc(values),
        acceptance = result$acceptance,
        final_scale = result$final_scale
      )
    }
    raw <- if (.Platform$OS.type == "unix") {
      parallel::mclapply(
        seq_len(settings$n_chains), run_one, mc.cores = settings$n_chains
      )
    } else {
      lapply(seq_len(settings$n_chains), run_one)
    }
    saveRDS(
      list(
        experiment_id = cfg$experiment_id,
        analysis = "alternative-prior-sensitivity",
        sim_id = id,
        model = model,
        input_md5 = input_md5,
        assumed_b = cfg$b,
        prior = prior,
        primary_prior = cfg$primary_prior,
        gh_estimate = if (model == "plugin") gh_fit$estimate else NULL,
        gh_fit = if (model == "plugin") gh_fit else NULL,
        chains = lapply(raw, `[[`, "draws"),
        acceptance = vapply(raw, `[[`, numeric(1), "acceptance"),
        final_scale = vapply(raw, `[[`, numeric(1), "final_scale"),
        settings = settings,
        initialisation_path = pilot[[model]]$path,
        seed_rule = paste0(
          "980000 + 1000 * sim_id + 100 * I(plugin) + chain"
        )
      ),
      output_path
    )
    cat("Finished alternative-prior fit:", model, id, "\n")
  }
}

# Fixed-b sensitivity analysis for the synthetic Plug-in ETAS fits.
#
# The generated catalogues are unchanged.  For each assumed b value this
# script re-estimates (G,H) from the observed catalogue and then targets the
# corresponding exact Plug-in ETAS posterior.  Outputs are isolated from the
# primary b=1 analysis.
#
# Required environment variable:
#   ASSUMED_B=0.8
# Optional environment variables:
#   MCMC_IDS=1,2,3 (takes precedence over MCMC_START/MCMC_END)
#   MCMC_START=1 MCMC_END=100
#   MCMC_ITER=12000 MCMC_WARMUP=4000 MCMC_THIN=4 MCMC_CHAINS=4
#   QUAD_NODES=16 MCMC_FORCE=true

source("config/submission_experiment.R")
source("R/fit_loglinear_mle.R")
source("R/mcmc_etas.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
assumed_b <- suppressWarnings(as.numeric(Sys.getenv("ASSUMED_B")))
if (!is.finite(assumed_b) || assumed_b <= 0) {
  stop("ASSUMED_B must be a positive number.")
}
ids <- if (nzchar(Sys.getenv("MCMC_IDS"))) {
  as.integer(strsplit(Sys.getenv("MCMC_IDS"), ",", fixed = TRUE)[[1]])
} else {
  seq.int(
    as.integer(Sys.getenv("MCMC_START", unset = "1")),
    as.integer(Sys.getenv("MCMC_END", unset = as.character(cfg$n_sim)))
  )
}
if (anyNA(ids) || !all(ids %in% seq_len(cfg$n_sim))) {
  stop("MCMC_IDS contains an invalid catalogue ID.")
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

b_tag <- paste0("b_", gsub("\\.", "p", format(assumed_b, trim = TRUE)))
out <- file.path(cfg$results_dir, "b_sensitivity", "synthetic", b_tag)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# The retained primary INLA pilot supplies only dispersed starts and proposal
# geometry.  It does not determine the posterior targeted by these chains.
initialisation_id <- as.integer(Sys.getenv("SENSITIVITY_PILOT_ID", unset = "1"))
initialisation_path <- file.path(
  cfg$results_dir, "initialisation", sprintf("plugin_%03d.rds", initialisation_id)
)
if (!file.exists(initialisation_path)) {
  stop("Missing pilot initialisation: ", initialisation_path)
}
initialisation <- readRDS(initialisation_path)
initialisation_input <- file.path(
  cfg$data_dir, sprintf("observed_%03d.rds", initialisation_id)
)
stopifnot(
  identical(
    initialisation$input_md5,
    unname(tools::md5sum(initialisation_input))
  ),
  isTRUE(initialisation$convergence$converged)
)
pilot_draws <- initialisation$posterior_draws[c("mu", "K", "alpha", "c", "p")]
pilot_z <- t(apply(
  pilot_draws, 1, etas_theta_to_z, prior = cfg$primary_prior
))
proposal_cov <- cov(pilot_z)
initial_index <- unique(round(seq(
  1, nrow(pilot_z), length.out = settings$n_chains + 2L
)))[seq_len(settings$n_chains) + 1L]
initial <- lapply(initial_index, function(index) pilot_z[index, ])

fit_is_reusable <- function(path, input_md5) {
  if (!file.exists(path) || force) return(FALSE)
  fit <- tryCatch(readRDS(path), error = function(e) NULL)
  compatible <- !is.null(fit) &&
    identical(fit$experiment_id, cfg$experiment_id) &&
    identical(fit$model, "plugin") &&
    identical(fit$input_md5, input_md5) &&
    isTRUE(all.equal(fit$assumed_b, assumed_b, tolerance = 1e-12)) &&
    fit$settings$iterations >= settings$iterations &&
    fit$settings$warmup >= settings$warmup &&
    identical(fit$settings$thin, settings$thin) &&
    identical(fit$settings$n_chains, settings$n_chains) &&
    identical(fit$settings$nodes_per_interval, settings$nodes_per_interval) &&
    length(fit$chains) == settings$n_chains
  isTRUE(compatible)
}

for (id in ids) {
  input_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", id))
  if (!file.exists(input_path)) stop("Missing catalogue: ", input_path)
  input_md5 <- unname(tools::md5sum(input_path))
  output_path <- file.path(out, sprintf("plugin_%03d.rds", id))
  if (fit_is_reusable(output_path, input_md5)) {
    cat("Reusing fixed-b Plug-in fit:", assumed_b, id, "\n")
    next
  }
  catalogue <- readRDS(input_path)
  gh_fit <- estimate_loglinear_plugin(
    catalogue,
    cfg$mainshock$ts,
    cfg$mainshock$magnitudes,
    cfg$M0,
    b = assumed_b
  )
  gh <- gh_fit$estimate
  log_posterior <- make_plugin_etas_log_posterior(
    catalogue, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior,
    cfg$mainshock$ts, cfg$mainshock$magnitudes,
    gh[["G"]], gh[["H"]], assumed_b, settings$nodes_per_interval
  )

  run_one <- function(chain) {
    result <- run_exact_etas_chain(
      log_posterior, initial[[chain]], proposal_cov,
      settings$iterations, settings$warmup, settings$thin,
      seed = 930000L + round(10000 * assumed_b) + 1000L * id + chain
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
      seq_len(settings$n_chains), run_one, mc.cores = settings$n_chains
    )
  } else {
    lapply(seq_len(settings$n_chains), run_one)
  }

  saveRDS(
    list(
      experiment_id = cfg$experiment_id,
      analysis = "fixed-b-sensitivity",
      sim_id = id,
      model = "plugin",
      input_md5 = input_md5,
      assumed_b = assumed_b,
      data_generating_b = cfg$b,
      prior = cfg$primary_prior,
      gh_estimate = gh,
      gh_fit = gh_fit,
      chains = lapply(raw_chains, `[[`, "draws"),
      acceptance = vapply(raw_chains, `[[`, numeric(1), "acceptance"),
      final_scale = vapply(raw_chains, `[[`, numeric(1), "final_scale"),
      settings = settings,
      initialisation_path = initialisation_path,
      seed_rule = paste0(
        "930000 + round(10000 * assumed_b) + 1000 * sim_id + chain"
      )
    ),
    output_path
  )
  cat("Finished fixed-b Plug-in fit:", assumed_b, id, "\n")
}

# Fixed-b sensitivity fits for the conditioned, single-trigger Ridgecrest
# Plug-in analysis.  The Naive temporal posterior is invariant to b and is
# therefore reused at the forecast-evaluation stage.

source("config/submission_ridgecrest.R")
source("R/ridgecrest_mcmc.R")
source("R/fit_loglinear_mle.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_ridgecrest_config()
assumed_b <- suppressWarnings(as.numeric(Sys.getenv("ASSUMED_B")))
if (!is.finite(assumed_b) || assumed_b <= 0) {
  stop("ASSUMED_B must be a positive number.")
}
settings <- list(
  iterations = as.integer(Sys.getenv("MCMC_ITER", unset = "24000")),
  warmup = as.integer(Sys.getenv("MCMC_WARMUP", unset = "8000")),
  thin = as.integer(Sys.getenv("MCMC_THIN", unset = "4")),
  n_chains = as.integer(Sys.getenv("MCMC_CHAINS", unset = "4")),
  nodes_per_interval = as.integer(Sys.getenv("QUAD_NODES", unset = "16"))
)
force <- identical(tolower(Sys.getenv("MCMC_FORCE")), "true")

train <- readRDS(cfg$train_path)
stopifnot(identical(attr(train, "experiment_id"), cfg$experiment_id))
train$include_likelihood <- TRUE
is_mainshock <- train$ID == cfg$mainshock_id
stopifnot(sum(is_mainshock) == 1L)
train$include_likelihood[is_mainshock] <- FALSE
mainshock_time <- train$ts[is_mainshock]
mainshock_magnitude <- train$magnitudes[is_mainshock]
cfg$T1 <- min(train$ts) - 1e-6
cfg$T2 <- cfg$train_end
input_md5 <- unname(tools::md5sum(cfg$train_path))

b_tag <- paste0("b_", gsub("\\.", "p", format(assumed_b, trim = TRUE)))
out <- file.path(cfg$results_dir, "b_sensitivity", b_tag)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(out, "plugin.rds")
if (file.exists(output_path) && !force) {
  existing <- readRDS(output_path)
  reusable <- identical(existing$experiment_id, cfg$experiment_id) &&
    identical(existing$analysis, "fixed-b-sensitivity") &&
    identical(existing$input_md5, input_md5) &&
    isTRUE(all.equal(existing$assumed_b, assumed_b, tolerance = 1e-12)) &&
    existing$settings$iterations >= settings$iterations &&
    existing$settings$warmup >= settings$warmup &&
    identical(existing$settings$thin, settings$thin) &&
    identical(existing$settings$n_chains, settings$n_chains) &&
    identical(existing$settings$nodes_per_interval, settings$nodes_per_interval)
  if (isTRUE(reusable)) {
    cat("Reusing Ridgecrest fixed-b fit:", assumed_b, "\n")
    quit(save = "no", status = 0L)
  }
}

gh_fit <- estimate_loglinear_plugin(
  train, mainshock_time, mainshock_magnitude, cfg$M0, b = assumed_b
)
gh <- gh_fit$estimate

primary_path <- file.path(
  cfg$results_dir, "mcmc_conditioned_single", "plugin.rds"
)
if (!file.exists(primary_path)) stop("Missing primary exact fit: ", primary_path)
primary <- readRDS(primary_path)
stopifnot(
  identical(primary$experiment_id, cfg$experiment_id),
  identical(primary$input_md5, input_md5),
  identical(primary$model, "plugin"),
  isTRUE(primary$condition_mainshock),
  identical(primary$completeness_mode, "single")
)
pilot_draws <- as.data.frame(do.call(rbind, primary$chains))
pilot_draws <- pilot_draws[c("mu", "K", "alpha", "c", "p")]
pilot_z <- t(apply(
  pilot_draws, 1, etas_theta_to_z, prior = cfg$primary_prior
))
proposal_cov <- cov(pilot_z)
initial_index <- unique(round(seq(
  1, nrow(pilot_z), length.out = settings$n_chains + 2L
)))[seq_len(settings$n_chains) + 1L]
initial <- lapply(initial_index, function(index) pilot_z[index, ])

log_posterior <- make_fixedtrigger_plugin_log_posterior(
  train, cfg$T1, cfg$T2, cfg$M0, cfg$primary_prior,
  gh[["G"]], gh[["H"]], mainshock_time, mainshock_magnitude,
  assumed_b, settings$nodes_per_interval
)
run_one <- function(chain) {
  result <- run_exact_etas_chain(
    log_posterior, initial[[chain]], proposal_cov,
    settings$iterations, settings$warmup, settings$thin,
    seed = 970700L + round(10000 * assumed_b) + chain
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
    model = "plugin",
    input_md5 = input_md5,
    assumed_b = assumed_b,
    primary_b = cfg$b,
    condition_mainshock = TRUE,
    completeness_mode = "single",
    conditioned_event_id = cfg$mainshock_id,
    n_history_events = nrow(train),
    n_likelihood_targets = sum(train$include_likelihood),
    prior = cfg$primary_prior,
    gh_estimate = gh,
    gh_fit = gh_fit,
    completeness_trigger = c(
      time = mainshock_time, magnitude = mainshock_magnitude
    ),
    chains = lapply(raw_chains, `[[`, "draws"),
    acceptance = vapply(raw_chains, `[[`, numeric(1), "acceptance"),
    final_scale = vapply(raw_chains, `[[`, numeric(1), "final_scale"),
    settings = settings,
    proposal_source = primary_path,
    seed_rule = "970700 + round(10000 * assumed_b) + chain"
  ),
  output_path
)
cat("Finished Ridgecrest fixed-b Plug-in fit:", assumed_b, "\n")

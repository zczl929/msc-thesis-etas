# MCMC diagnostics and fixed-period posterior-predictive Ridgecrest evaluation.

source("config/submission_ridgecrest.R")
source("R/forecast_observed_etas.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_ridgecrest_config()
assumed_b <- suppressWarnings(as.numeric(Sys.getenv(
  "ASSUMED_B", unset = as.character(cfg$b)
)))
if (!is.finite(assumed_b) || assumed_b <= 0) stop("Invalid ASSUMED_B.")
cfg$b <- assumed_b
full <- readRDS(cfg$full_path)
train <- readRDS(cfg$train_path)
condition_mainshock <- identical(
  tolower(Sys.getenv("CONDITION_MAINSHOCK", unset = "false")), "true"
)
completeness_mode <- tolower(Sys.getenv(
  "COMPLETENESS_MODE", unset = "single"
))
if (!completeness_mode %in% c("single", "multi")) {
  stop("COMPLETENESS_MODE must be 'single' or 'multi'.")
}
train$include_likelihood <- TRUE
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
history <- full[full$ts <= cfg$train_end, , drop = FALSE]
test <- full[
  full$ts > cfg$train_end & full$ts <= cfg$forecast_end,
  ,
  drop = FALSE
]
gh_path <- file.path(
  cfg$results_dir,
  if (completeness_mode == "single") "plugin_GH_single.rds" else "plugin_GH.rds"
)
gh <- readRDS(gh_path)$estimate
analysis_name <- paste0(
  "mcmc",
  if (condition_mainshock) "_conditioned" else "",
  if (completeness_mode == "single") "_single" else ""
)
fit_dir <- file.path(cfg$results_dir, analysis_name)
naive_fit_path <- Sys.getenv(
  "NAIVE_FIT_PATH", unset = file.path(fit_dir, "naive.rds")
)
plugin_fit_path <- Sys.getenv(
  "PLUGIN_FIT_PATH", unset = file.path(fit_dir, "plugin.rds")
)
fit_paths <- c(naive = naive_fit_path, plugin = plugin_fit_path)
plugin_fit_for_gh <- readRDS(plugin_fit_path)
if (!is.null(plugin_fit_for_gh$gh_estimate)) {
  gh <- plugin_fit_for_gh$gh_estimate
}
out <- Sys.getenv("FORECAST_OUT", unset = file.path(fit_dir, "forecast"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

diagnostic_rows <- list()
summary_rows <- list()
posterior_draws <- list()
for (model in c("naive", "plugin")) {
  fit <- readRDS(fit_paths[[model]])
  stopifnot(
    identical(fit$experiment_id, cfg$experiment_id),
    identical(fit$input_md5, input_md5),
    identical(fit$model, model),
    identical(fit$condition_mainshock, condition_mainshock),
    identical(fit$completeness_mode, completeness_mode),
    identical(fit$n_history_events, nrow(train)),
    identical(fit$n_likelihood_targets, sum(train$include_likelihood))
  )
  chains <- coda::mcmc.list(fit$chains)
  gelman <- coda::gelman.diag(
    chains, autoburnin = FALSE, multivariate = FALSE
  )$psrf[, "Point est."]
  ess <- coda::effectiveSize(chains)
  combined <- as.data.frame(do.call(rbind, fit$chains))
  posterior_draws[[model]] <- combined
  for (parameter in names(gelman)) {
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      model = model,
      parameter = parameter,
      rhat = gelman[[parameter]],
      effective_size = ess[[parameter]],
      acceptance_min = min(fit$acceptance),
      acceptance_max = max(fit$acceptance)
    )
    interval <- quantile(combined[[parameter]], c(0.025, 0.5, 0.975))
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      model = model,
      parameter = parameter,
      mean = mean(combined[[parameter]]),
      median = interval[[2L]],
      q025 = interval[[1L]],
      q975 = interval[[3L]]
    )
  }
}
diagnostics <- do.call(rbind, diagnostic_rows)
posterior_summary <- do.call(rbind, summary_rows)
write.csv(diagnostics, file.path(out, "mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(
  posterior_summary, file.path(out, "posterior_summary.csv"), row.names = FALSE
)
failures <- diagnostics[
  diagnostics$rhat > 1.05 | diagnostics$effective_size < 400,
]
if (nrow(failures)) {
  print(failures)
  stop("Ridgecrest MCMC convergence gate failed.")
}

n_predictive <- min(
  cfg$n_predictive,
  min(vapply(posterior_draws, nrow, integer(1)))
)
rows <- list()
predictive <- list()
for (model in c("naive", "plugin")) {
  draws <- posterior_draws[[model]]
  set.seed(if (model == "naive") 6072301L else 6072302L)
  index <- sample(seq_len(nrow(draws)), n_predictive, replace = FALSE)
  log_lik <- numeric(n_predictive)
  counts <- numeric(n_predictive)
  max_magnitude <- rep(NA_real_, n_predictive)
  overflow <- logical(n_predictive)

  for (k in seq_len(n_predictive)) {
    theta <- as.list(draws[index[k], c("mu", "K", "alpha", "c", "p")])
    log_lik[k] <- multi_observed_loglik(
      theta, full, cfg$train_end, cfg$forecast_end, cfg$M0, cfg$b,
      model, gh[["G"]], gh[["H"]], cfg$M_trigger,
      fixed_trigger_time = if (completeness_mode == "single") {
        mainshock_time
      } else NULL,
      fixed_trigger_magnitude = if (completeness_mode == "single") {
        mainshock_magnitude
      } else NULL
    )
    simulated <- simulate_observed_history_etas(
      theta, history, cfg$train_end, cfg$forecast_end, cfg$M0, cfg$b,
      model, gh[["G"]], gh[["H"]], cfg$M_trigger,
      fixed_trigger_time = if (completeness_mode == "single") {
        mainshock_time
      } else NULL,
      fixed_trigger_magnitude = if (completeness_mode == "single") {
        mainshock_magnitude
      } else NULL
    )
    counts[k] <- nrow(simulated)
    if (nrow(simulated)) max_magnitude[k] <- max(simulated$magnitudes)
    overflow[k] <- isTRUE(attr(simulated, "overflow"))
    if (k %% 100L == 0L) cat(model, k, "of", n_predictive, "\n")
  }

  count_q <- quantile(counts, c(0.025, 0.5, 0.975))
  magnitude_q <- quantile(
    max_magnitude, c(0.025, 0.5, 0.975), na.rm = TRUE
  )
  branching <- ifelse(
    draws$alpha < cfg$b,
    draws$K * cfg$b / (cfg$b - draws$alpha),
    Inf
  )
  rows[[model]] <- data.frame(
    model = model,
    assumed_b = cfg$b,
    condition_mainshock = condition_mainshock,
    completeness_mode = completeness_mode,
    n_history_events = nrow(train),
    n_likelihood_targets = sum(train$include_likelihood),
    n_test = nrow(test),
    log_score = log_mean_exp(log_lik),
    log_score_per_event = log_mean_exp(log_lik) / nrow(test),
    count_q025 = count_q[[1L]],
    count_median = count_q[[2L]],
    count_q975 = count_q[[3L]],
    count_covered =
      count_q[[1L]] <= nrow(test) && nrow(test) <= count_q[[3L]],
    overflow_rate = mean(overflow),
    observed_max = max(test$magnitudes),
    maxmag_q025 = magnitude_q[[1L]],
    maxmag_median = magnitude_q[[2L]],
    maxmag_q975 = magnitude_q[[3L]],
    branching_median = median(branching),
    probability_supercritical = mean(branching >= 1),
    probability_alpha_ge_b = mean(draws$alpha >= cfg$b)
  )
  predictive[[model]] <- list(
    loglik = log_lik,
    counts = counts,
    max_magnitude = max_magnitude,
    overflow = overflow,
    posterior_index = index
  )
}

result <- do.call(rbind, rows)
result$information_gain_vs_naive <- (
  result$log_score - result$log_score[result$model == "naive"]
) / result$n_test
write.csv(result, file.path(out, "forecast_summary.csv"), row.names = FALSE)
forecast_count_draws <- do.call(rbind, lapply(names(predictive), function(model) {
  data.frame(
    model = model,
    draw = seq_along(predictive[[model]]$counts),
    count = predictive[[model]]$counts,
    overflow = predictive[[model]]$overflow
  )
}))
write.csv(
  forecast_count_draws,
  file.path(out, "forecast_count_draws.csv"),
  row.names = FALSE
)
saveRDS(predictive, file.path(out, "forecast_draws.rds"))

set.seed(6072303L)
B <- 10000L
information_gain <- replicate(B, {
  naive <- sample(predictive$naive$loglik, replace = TRUE)
  plugin <- sample(predictive$plugin$loglik, replace = TRUE)
  (log_mean_exp(plugin) - log_mean_exp(naive)) / nrow(test)
})
mc <- data.frame(
  n_posterior_draws = length(predictive$naive$loglik),
  information_gain_per_event =
    result$information_gain_vs_naive[result$model == "plugin"],
  mc_q025 = unname(quantile(information_gain, 0.025)),
  mc_q975 = unname(quantile(information_gain, 0.975)),
  probability_positive = mean(information_gain > 0),
  note = paste(
    "Interval quantifies posterior Monte Carlo integration error,",
    "not generalisation across earthquake sequences."
  )
)
write.csv(mc, file.path(out, "forecast_mc_error.csv"), row.names = FALSE)

print(diagnostics, row.names = FALSE)
print(posterior_summary, row.names = FALSE)
print(result, row.names = FALSE)
print(mc, row.names = FALSE)

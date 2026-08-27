# Evaluate the Oracle Plug-in and complete-data synthetic benchmarks without
# overwriting the frozen Naive and estimated Plug-in evaluation outputs.

source("config/submission_experiment.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
fit_dir <- file.path(cfg$results_dir, "mcmc_primary")
out <- file.path(fit_dir, "summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
models <- c("naive", "plugin", "oracle", "complete")

wilson_interval <- function(successes, n, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  estimate <- successes / n
  denominator <- 1 + z^2 / n
  centre <- (estimate + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(
    estimate * (1 - estimate) / n + z^2 / (4 * n^2)
  ) / denominator
  c(lower = centre - half_width, upper = centre + half_width)
}

posterior_rows <- list()
diagnostic_rows <- list()
k <- 0L
d <- 0L
for (id in seq_len(cfg$n_sim)) {
  for (model in c("oracle", "complete")) {
    catalogue <- if (model == "complete") "complete" else "observed"
    input_path <- file.path(cfg$data_dir, sprintf("%s_%03d.rds", catalogue, id))
    fit_path <- file.path(fit_dir, sprintf("%s_%03d.rds", model, id))
    if (!file.exists(fit_path)) stop("Missing benchmark fit: ", fit_path)
    fit <- readRDS(fit_path)
    stopifnot(
      identical(fit$experiment_id, cfg$experiment_id),
      identical(fit$input_md5, unname(tools::md5sum(input_path))),
      identical(fit$model, model),
      length(fit$chains) == fit$settings$n_chains
    )
    chains <- coda::mcmc.list(fit$chains)
    combined <- as.data.frame(do.call(rbind, fit$chains))
    rhat <- coda::gelman.diag(
      chains, autoburnin = FALSE, multivariate = FALSE
    )$psrf[, "Point est."]
    ess <- coda::effectiveSize(chains)
    acceptance <- fit$acceptance
    for (parameter in names(cfg$theta_true)) {
      interval <- quantile(combined[[parameter]], c(0.025, 0.5, 0.975))
      truth <- cfg$theta_true[[parameter]]
      k <- k + 1L
      posterior_rows[[k]] <- data.frame(
        sim_id = id, model = model, parameter = parameter, true = truth,
        median = interval[[2L]], q025 = interval[[1L]], q975 = interval[[3L]],
        width = interval[[3L]] - interval[[1L]],
        covered = interval[[1L]] <= truth && truth <= interval[[3L]]
      )
      d <- d + 1L
      diagnostic_rows[[d]] <- data.frame(
        sim_id = id, model = model, parameter = parameter,
        rhat = rhat[[parameter]], effective_size = ess[[parameter]],
        acceptance_min = min(acceptance), acceptance_max = max(acceptance)
      )
    }
  }
}

benchmark_posterior <- do.call(rbind, posterior_rows)
benchmark_diagnostics <- do.call(rbind, diagnostic_rows)
primary_posterior <- read.csv(file.path(out, "posterior_long.csv"))
primary_diagnostics <- read.csv(file.path(out, "mcmc_diagnostics.csv"))
posterior <- rbind(primary_posterior, benchmark_posterior)
diagnostics <- rbind(primary_diagnostics, benchmark_diagnostics)

model_summary <- do.call(rbind, lapply(
  split(posterior, list(posterior$model, posterior$parameter), drop = TRUE),
  function(x) {
    coverage <- mean(x$covered)
    coverage_ci <- wilson_interval(sum(x$covered), nrow(x))
    model <- x$model[1L]
    parameter <- x$parameter[1L]
    dx <- diagnostics[
      diagnostics$model == model & diagnostics$parameter == parameter,
    ]
    data.frame(
      model = model, parameter = parameter, n_catalogues = nrow(x),
      bias = mean(x$median - x$true),
      rmse = sqrt(mean((x$median - x$true)^2)),
      coverage = coverage,
      coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
      coverage_ci_low = coverage_ci[["lower"]],
      coverage_ci_high = coverage_ci[["upper"]],
      mean_interval_width = mean(x$width),
      median_interval_width = median(x$width),
      max_rhat = max(dx$rhat), min_effective_size = min(dx$effective_size)
    )
  }
))

wide <- reshape(
  posterior[c("sim_id", "parameter", "true", "model", "median", "covered", "width")],
  idvar = c("sim_id", "parameter", "true"), timevar = "model",
  direction = "wide"
)
decomposition <- do.call(rbind, lapply(split(wide, wide$parameter), function(x) {
  rmse <- vapply(models, function(m) {
    sqrt(mean((x[[paste0("median.", m)]] - x$true)^2))
  }, numeric(1))
  coverage <- vapply(models, function(m) {
    mean(x[[paste0("covered.", m)]])
  }, numeric(1))
  data.frame(
    parameter = x$parameter[1L], n_catalogues = nrow(x),
    rmse_naive = rmse[["naive"]], rmse_plugin = rmse[["plugin"]],
    rmse_oracle = rmse[["oracle"]], rmse_complete = rmse[["complete"]],
    coverage_naive = coverage[["naive"]],
    coverage_plugin = coverage[["plugin"]],
    coverage_oracle = coverage[["oracle"]],
    coverage_complete = coverage[["complete"]],
    plugin_minus_oracle_rmse = rmse[["plugin"]] - rmse[["oracle"]],
    oracle_minus_complete_rmse = rmse[["oracle"]] - rmse[["complete"]]
  )
}))

failures <- benchmark_diagnostics[
  benchmark_diagnostics$rhat > 1.05 |
    benchmark_diagnostics$effective_size < 400,
]
write.csv(benchmark_posterior, file.path(out, "benchmark_posterior_long.csv"), row.names = FALSE)
write.csv(benchmark_diagnostics, file.path(out, "benchmark_mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(model_summary, file.path(out, "benchmark_model_summary.csv"), row.names = FALSE)
write.csv(wide, file.path(out, "benchmark_paired_results.csv"), row.names = FALSE)
write.csv(decomposition, file.path(out, "benchmark_decomposition.csv"), row.names = FALSE)
write.csv(failures, file.path(out, "benchmark_convergence_failures.csv"), row.names = FALSE)

print(decomposition, row.names = FALSE)
cat("\nBenchmark convergence failures (Rhat > 1.05 or ESS < 400):",
    nrow(failures), "\n")

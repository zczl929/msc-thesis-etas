# Evaluate the final 100-catalogue MCMC experiment.

source("config/submission_experiment.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
fit_dir <- file.path(cfg$results_dir, "mcmc_primary")
out <- file.path(fit_dir, "summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

posterior_rows <- list()
diagnostic_rows <- list()
index <- 0L
diagnostic_index <- 0L

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

for (id in seq_len(cfg$n_sim)) {
  input_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", id))
  input_md5 <- unname(tools::md5sum(input_path))
  for (model in c("naive", "plugin")) {
    path <- file.path(fit_dir, sprintf("%s_%03d.rds", model, id))
    if (!file.exists(path)) stop("Missing MCMC fit: ", path)
    fit <- readRDS(path)
    stopifnot(
      identical(fit$experiment_id, cfg$experiment_id),
      identical(fit$input_md5, input_md5),
      identical(fit$model, model),
      length(fit$chains) == fit$settings$n_chains
    )
    chains <- coda::mcmc.list(fit$chains)
    combined <- as.data.frame(do.call(rbind, fit$chains))
    gelman <- coda::gelman.diag(
      chains, autoburnin = FALSE, multivariate = FALSE
    )$psrf[, "Point est."]
    ess <- coda::effectiveSize(chains)
    acceptance <- if (!is.null(fit$acceptance)) fit$acceptance else NA_real_
    acceptance_min <- if (all(is.na(acceptance))) {
      NA_real_
    } else {
      min(acceptance, na.rm = TRUE)
    }
    acceptance_max <- if (all(is.na(acceptance))) {
      NA_real_
    } else {
      max(acceptance, na.rm = TRUE)
    }

    for (parameter in names(cfg$theta_true)) {
      values <- combined[[parameter]]
      interval <- quantile(values, c(0.025, 0.5, 0.975))
      truth <- cfg$theta_true[[parameter]]
      index <- index + 1L
      posterior_rows[[index]] <- data.frame(
        sim_id = id,
        model = model,
        parameter = parameter,
        true = truth,
        median = interval[[2L]],
        q025 = interval[[1L]],
        q975 = interval[[3L]],
        width = interval[[3L]] - interval[[1L]],
        covered = interval[[1L]] <= truth && truth <= interval[[3L]]
      )
      diagnostic_index <- diagnostic_index + 1L
      diagnostic_rows[[diagnostic_index]] <- data.frame(
        sim_id = id,
        model = model,
        parameter = parameter,
        rhat = gelman[[parameter]],
        effective_size = ess[[parameter]],
        acceptance_min = acceptance_min,
        acceptance_max = acceptance_max
      )
    }
  }
}

posterior <- do.call(rbind, posterior_rows)
diagnostics <- do.call(rbind, diagnostic_rows)

aggregate_results <- do.call(rbind, lapply(
  split(
    posterior,
    list(posterior$model, posterior$parameter),
    drop = TRUE
  ),
  function(x) {
    coverage <- mean(x$covered)
    coverage_interval <- wilson_interval(sum(x$covered), nrow(x))
    data.frame(
      model = x$model[1L],
      parameter = x$parameter[1L],
      n_catalogues = nrow(x),
      bias = mean(x$median - x$true),
      rmse = sqrt(mean((x$median - x$true)^2)),
      coverage = coverage,
      coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
      coverage_ci_low = coverage_interval[["lower"]],
      coverage_ci_high = coverage_interval[["upper"]],
      mean_interval_width = mean(x$width),
      median_interval_width = median(x$width),
      max_rhat = max(diagnostics$rhat[
        diagnostics$model == x$model[1L] &
          diagnostics$parameter == x$parameter[1L]
      ]),
      min_effective_size = min(diagnostics$effective_size[
        diagnostics$model == x$model[1L] &
          diagnostics$parameter == x$parameter[1L]
      ])
    )
  }
))

paired <- merge(
  posterior[posterior$model == "naive", ],
  posterior[posterior$model == "plugin", ],
  by = c("sim_id", "parameter", "true"),
  suffixes = c("_naive", "_plugin")
)
paired$absolute_error_change <-
  abs(paired$median_plugin - paired$true) -
  abs(paired$median_naive - paired$true)
paired$plugin_improved_absolute_error <- paired$absolute_error_change < 0

paired_results <- do.call(rbind, lapply(
  split(paired, paired$parameter),
  function(x) {
    improved_interval <- wilson_interval(
      sum(x$plugin_improved_absolute_error), nrow(x)
    )
    rmse_naive <- sqrt(mean((x$median_naive - x$true)^2))
    rmse_plugin <- sqrt(mean((x$median_plugin - x$true)^2))
    data.frame(
      parameter = x$parameter[1L],
      n_catalogues = nrow(x),
      rmse_reduction_fraction = 1 - rmse_plugin / rmse_naive,
      mean_absolute_error_change = mean(x$absolute_error_change),
      mean_absolute_error_change_mcse =
        sd(x$absolute_error_change) / sqrt(nrow(x)),
      median_absolute_error_change = median(x$absolute_error_change),
      proportion_plugin_improved = mean(x$plugin_improved_absolute_error),
      proportion_improved_ci_low = improved_interval[["lower"]],
      proportion_improved_ci_high = improved_interval[["upper"]],
      plugin_minus_naive_coverage =
        mean(x$covered_plugin) - mean(x$covered_naive),
      mean_width_ratio_plugin_over_naive =
        mean(x$width_plugin / x$width_naive)
    )
  }
))

convergence_failures <- diagnostics[
  diagnostics$rhat > 1.05 | diagnostics$effective_size < 400,
]

write.csv(posterior, file.path(out, "posterior_long.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(
  aggregate_results, file.path(out, "model_summary.csv"), row.names = FALSE
)
write.csv(paired, file.path(out, "paired_results.csv"), row.names = FALSE)
write.csv(
  paired_results, file.path(out, "paired_summary.csv"), row.names = FALSE
)
write.csv(
  convergence_failures,
  file.path(out, "convergence_failures.csv"),
  row.names = FALSE
)

print(aggregate_results, row.names = FALSE)
cat("\nPaired comparison\n")
print(paired_results, row.names = FALSE)
cat("\nConvergence failures (Rhat > 1.05 or ESS < 400):",
    nrow(convergence_failures), "\n")
if (nrow(convergence_failures)) print(convergence_failures, row.names = FALSE)

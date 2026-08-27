# Evaluate the alternative-prior Naive and Plug-in synthetic fits and compare
# them with the frozen primary-prior summaries.

source("config/submission_experiment.R")
if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
fit_dir <- file.path(cfg$results_dir, "prior_sensitivity", "synthetic")
out <- file.path(fit_dir, "summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

posterior_rows <- list()
diagnostic_rows <- list()
for (id in seq_len(cfg$n_sim)) {
  input_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", id))
  input_md5 <- unname(tools::md5sum(input_path))
  for (model in c("naive", "plugin")) {
    path <- file.path(fit_dir, sprintf("%s_%03d.rds", model, id))
    if (!file.exists(path)) stop("Missing alternative-prior fit: ", path)
    fit <- readRDS(path)
    stopifnot(
      identical(fit$experiment_id, cfg$experiment_id),
      identical(fit$analysis, "alternative-prior-sensitivity"),
      identical(fit$model, model),
      identical(fit$input_md5, input_md5),
      identical(fit$prior, cfg$sensitivity_prior)
    )
    chains <- coda::mcmc.list(fit$chains)
    combined <- as.data.frame(do.call(rbind, fit$chains))
    rhat <- coda::gelman.diag(
      chains, autoburnin = FALSE, multivariate = FALSE
    )$psrf[, "Point est."]
    ess <- coda::effectiveSize(chains)
    for (parameter in names(cfg$theta_true)) {
      q <- quantile(combined[[parameter]], c(0.025, 0.5, 0.975))
      truth <- cfg$theta_true[[parameter]]
      posterior_rows[[length(posterior_rows) + 1L]] <- data.frame(
        prior = "alternative", sim_id = id, model = model,
        parameter = parameter, true = truth, median = q[[2L]],
        q025 = q[[1L]], q975 = q[[3L]], width = q[[3L]] - q[[1L]],
        covered = q[[1L]] <= truth && truth <= q[[3L]]
      )
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        prior = "alternative", sim_id = id, model = model,
        parameter = parameter, rhat = rhat[[parameter]],
        effective_size = ess[[parameter]],
        acceptance_min = min(fit$acceptance),
        acceptance_max = max(fit$acceptance)
      )
    }
  }
}
alternative <- do.call(rbind, posterior_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
primary <- read.csv(file.path(
  cfg$results_dir, "mcmc_primary", "summary", "posterior_long.csv"
))
primary$prior <- "primary"
posterior <- rbind(
  primary[c("prior", "sim_id", "model", "parameter", "true", "median",
            "q025", "q975", "width", "covered")],
  alternative
)

summarise_group <- function(x) data.frame(
  prior = x$prior[1L], model = x$model[1L], parameter = x$parameter[1L],
  n_catalogues = nrow(x),
  bias = mean(x$median - x$true),
  rmse = sqrt(mean((x$median - x$true)^2)),
  coverage = mean(x$covered),
  coverage_mcse = sqrt(mean(x$covered) * (1 - mean(x$covered)) / nrow(x)),
  mean_interval_width = mean(x$width),
  median_interval_width = median(x$width)
)
model_summary <- do.call(rbind, lapply(
  split(
    posterior,
    list(posterior$prior, posterior$model, posterior$parameter),
    drop = TRUE
  ),
  summarise_group
))

paired_rows <- list()
for (prior_name in c("primary", "alternative")) {
  x <- posterior[posterior$prior == prior_name, ]
  paired <- merge(
    x[x$model == "naive", ], x[x$model == "plugin", ],
    by = c("prior", "sim_id", "parameter", "true"),
    suffixes = c("_naive", "_plugin")
  )
  paired_rows[[prior_name]] <- do.call(rbind, lapply(
    split(paired, paired$parameter),
    function(y) {
      rmse_naive <- sqrt(mean((y$median_naive - y$true)^2))
      rmse_plugin <- sqrt(mean((y$median_plugin - y$true)^2))
      data.frame(
        prior = prior_name, parameter = y$parameter[1L],
        n_catalogues = nrow(y),
        naive_rmse = rmse_naive, plugin_rmse = rmse_plugin,
        rmse_reduction_fraction = 1 - rmse_plugin / rmse_naive,
        naive_coverage = mean(y$covered_naive),
        plugin_coverage = mean(y$covered_plugin),
        plugin_minus_naive_coverage =
          mean(y$covered_plugin) - mean(y$covered_naive)
      )
    }
  ))
}
paired_summary <- do.call(rbind, paired_rows)
convergence_failures <- diagnostics[
  diagnostics$rhat > 1.05 | diagnostics$effective_size < 400,
]

write.csv(posterior, file.path(out, "posterior_long.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(model_summary, file.path(out, "model_summary.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(out, "paired_summary.csv"), row.names = FALSE)
write.csv(
  convergence_failures, file.path(out, "convergence_failures.csv"),
  row.names = FALSE
)
print(model_summary, row.names = FALSE)
cat("\nPrimary comparison under each prior\n")
print(paired_summary, row.names = FALSE)
cat("\nConvergence failures:", nrow(convergence_failures), "\n")

# Summarise fixed-b sensitivity of the synthetic Plug-in ETAS analysis.

source("config/submission_experiment.R")

if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")

cfg <- submission_experiment_config()
b_values <- as.numeric(strsplit(
  Sys.getenv("B_VALUES", unset = "0.8,1,1.2"), ",", fixed = TRUE
)[[1]])
if (any(!is.finite(b_values)) || any(b_values <= 0)) stop("Invalid B_VALUES.")

root <- file.path(cfg$results_dir, "b_sensitivity", "synthetic")
out <- file.path(root, "summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

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
gh_rows <- list()
pi <- 0L
di <- 0L
gi <- 0L

for (assumed_b in b_values) {
  if (isTRUE(all.equal(assumed_b, cfg$b, tolerance = 1e-12))) {
    primary_path <- file.path(
      cfg$results_dir, "mcmc_primary", "summary", "posterior_long.csv"
    )
    primary_diagnostics_path <- file.path(
      cfg$results_dir, "mcmc_primary", "summary", "mcmc_diagnostics.csv"
    )
    primary <- subset(read.csv(primary_path), model == "plugin")
    primary_diag <- subset(read.csv(primary_diagnostics_path), model == "plugin")
    primary$assumed_b <- assumed_b
    primary_diag$assumed_b <- assumed_b
    posterior_rows[[length(posterior_rows) + 1L]] <- primary
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- primary_diag
    next
  }

  b_tag <- paste0("b_", gsub("\\.", "p", format(assumed_b, trim = TRUE)))
  fit_dir <- file.path(root, b_tag)
  for (id in seq_len(cfg$n_sim)) {
    input_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", id))
    input_md5 <- unname(tools::md5sum(input_path))
    path <- file.path(fit_dir, sprintf("plugin_%03d.rds", id))
    if (!file.exists(path)) stop("Missing sensitivity fit: ", path)
    fit <- readRDS(path)
    stopifnot(
      identical(fit$experiment_id, cfg$experiment_id),
      identical(fit$analysis, "fixed-b-sensitivity"),
      identical(fit$model, "plugin"),
      identical(fit$input_md5, input_md5),
      isTRUE(all.equal(fit$assumed_b, assumed_b, tolerance = 1e-12)),
      length(fit$chains) == fit$settings$n_chains
    )
    chains <- coda::mcmc.list(fit$chains)
    combined <- as.data.frame(do.call(rbind, fit$chains))
    rhat <- coda::gelman.diag(
      chains, autoburnin = FALSE, multivariate = FALSE
    )$psrf[, "Point est."]
    ess <- coda::effectiveSize(chains)
    for (parameter in names(cfg$theta_true)) {
      interval <- quantile(combined[[parameter]], c(0.025, 0.5, 0.975))
      truth <- cfg$theta_true[[parameter]]
      pi <- pi + 1L
      posterior_rows[[length(posterior_rows) + 1L]] <- data.frame(
        sim_id = id,
        model = "plugin",
        parameter = parameter,
        true = truth,
        median = interval[[2L]],
        q025 = interval[[1L]],
        q975 = interval[[3L]],
        width = interval[[3L]] - interval[[1L]],
        covered = interval[[1L]] <= truth && truth <= interval[[3L]],
        assumed_b = assumed_b
      )
      di <- di + 1L
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        sim_id = id,
        model = "plugin",
        parameter = parameter,
        rhat = rhat[[parameter]],
        effective_size = ess[[parameter]],
        acceptance_min = min(fit$acceptance),
        acceptance_max = max(fit$acceptance),
        assumed_b = assumed_b
      )
    }
    gi <- gi + 1L
    gh_rows[[gi]] <- data.frame(
      sim_id = id,
      assumed_b = assumed_b,
      G = fit$gh_estimate[["G"]],
      H = fit$gh_estimate[["H"]]
    )
  }
}

posterior <- do.call(rbind, posterior_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
gh <- if (length(gh_rows)) do.call(rbind, gh_rows) else data.frame()

summary_rows <- lapply(
  split(posterior, list(posterior$assumed_b, posterior$parameter), drop = TRUE),
  function(x) {
    coverage <- mean(x$covered)
    interval <- wilson_interval(sum(x$covered), nrow(x))
    relevant_diag <- diagnostics[
      diagnostics$assumed_b == x$assumed_b[1L] &
        diagnostics$parameter == x$parameter[1L],
    ]
    data.frame(
      assumed_b = x$assumed_b[1L],
      parameter = x$parameter[1L],
      n_catalogues = nrow(x),
      bias = mean(x$median - x$true),
      rmse = sqrt(mean((x$median - x$true)^2)),
      coverage = coverage,
      coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(x)),
      coverage_ci_low = interval[["lower"]],
      coverage_ci_high = interval[["upper"]],
      mean_interval_width = mean(x$width),
      median_interval_width = median(x$width),
      max_rhat = max(relevant_diag$rhat),
      min_effective_size = min(relevant_diag$effective_size)
    )
  }
)
model_summary <- do.call(rbind, summary_rows)
model_summary <- model_summary[order(model_summary$parameter, model_summary$assumed_b), ]
convergence_failures <- diagnostics[
  diagnostics$rhat > 1.05 | diagnostics$effective_size < 400,
]

write.csv(posterior, file.path(out, "posterior_long.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(model_summary, file.path(out, "model_summary.csv"), row.names = FALSE)
write.csv(gh, file.path(out, "gh_estimates.csv"), row.names = FALSE)
write.csv(
  convergence_failures, file.path(out, "convergence_failures.csv"),
  row.names = FALSE
)

print(model_summary, row.names = FALSE)
cat("\nConvergence failures (Rhat > 1.05 or ESS < 400):",
    nrow(convergence_failures), "\n")

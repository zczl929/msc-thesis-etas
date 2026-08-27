#!/usr/bin/env Rscript

# Fast, read-only verification of the dissertation claims and submission
# structure. Run from the repository root after generating the final summaries.

fail <- function(...) stop(..., call. = FALSE)
assert_close <- function(actual, expected, tolerance = 5e-4, label = "value") {
  if (length(actual) != 1L || is.na(actual) ||
      abs(actual - expected) > tolerance) {
    fail(label, ": expected ", expected, ", found ", actual)
  }
}

required_files <- c(
  "results/submission_v1/mcmc_primary/summary/model_summary.csv",
  "results/submission_v1/mcmc_primary/summary/paired_summary.csv",
  "results/submission_v1/mcmc_primary/summary/mcmc_diagnostics.csv",
  "results/submission_v1/mcmc_primary/summary/benchmark_decomposition.csv",
  "results/submission_v1/mcmc_primary/summary/benchmark_model_summary.csv",
  paste0(
    "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/",
    c("forecast_summary.csv", "forecast_mc_error.csv", "mcmc_diagnostics.csv")
  ),
  "results/submission_v1/ridgecrest/completeness_summary.csv",
  paste0(
    "results/submission_v1/figures/thesis/",
    c(
      "02_synthetic_observation_curve.csv",
      "02_synthetic_observation_mechanism_data.csv",
      "03_synthetic_parameter_recovery_data.csv",
      "04_synthetic_empirical_coverage_data.csv",
      "05_synthetic_benchmark_rmse_data.csv",
      "05_ridgecrest_completeness_curve.csv",
      "06_ridgecrest_prediction_summary_data.csv",
      "S1_synthetic_mcmc_diagnostics_data.csv",
      "S2_fixed_b_synthetic_sensitivity_data.csv",
      "S3_ridgecrest_fixed_b_forecast_data.csv",
      "S4_ridgecrest_fixed_b_posterior_data.csv",
      "ridgecrest_fixed_b_predictive_scores_data.csv"
    )
  ),
  paste0(
    "results/submission_v1/figures/thesis/",
    c(
      "02_synthetic_observation_mechanism.pdf",
      "03_synthetic_parameter_recovery.pdf",
      "04_synthetic_empirical_coverage.pdf",
      "05_ridgecrest_catalogue_completeness.pdf",
      "06_ridgecrest_prediction_summary.pdf",
      "S1_synthetic_mcmc_diagnostics.pdf",
      "S2_fixed_b_synthetic_sensitivity.pdf",
      "S3_ridgecrest_fixed_b_forecast.pdf"
    )
  )
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  fail("Missing required submission files:\n", paste(missing_files, collapse = "\n"))
}

model <- read.csv(
  "results/submission_v1/mcmc_primary/summary/model_summary.csv",
  check.names = FALSE
)
paired <- read.csv(
  "results/submission_v1/mcmc_primary/summary/paired_summary.csv",
  check.names = FALSE
)
synthetic_diagnostics <- read.csv(
  "results/submission_v1/mcmc_primary/summary/mcmc_diagnostics.csv",
  check.names = FALSE
)
forecast <- read.csv(
  paste0(
    "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/",
    "forecast_summary.csv"
  ),
  check.names = FALSE
)
forecast_mc <- read.csv(
  paste0(
    "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/",
    "forecast_mc_error.csv"
  ),
  check.names = FALSE
)
ridgecrest_diagnostics <- read.csv(
  paste0(
    "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/",
    "mcmc_diagnostics.csv"
  ),
  check.names = FALSE
)
completeness <- read.csv(
  "results/submission_v1/ridgecrest/completeness_summary.csv",
  check.names = FALSE
)

get_model <- function(model_name, parameter_name, column) {
  row <- model$model == model_name & model$parameter == parameter_name
  if (sum(row) != 1L) fail("Missing model summary row: ", model_name, "/", parameter_name)
  model[[column]][row]
}
get_paired <- function(parameter_name, column) {
  row <- paired$parameter == parameter_name
  if (sum(row) != 1L) fail("Missing paired summary row: ", parameter_name)
  paired[[column]][row]
}
get_forecast <- function(model_name, column) {
  row <- forecast$model == model_name
  if (sum(row) != 1L) fail("Missing forecast row: ", model_name)
  forecast[[column]][row]
}

# Table 1 and the main paired-comparison claims.
expected <- list(
  mu = c(naive_bias = -0.0074, naive_rmse = 0.0148, naive_coverage = 0.96,
         plugin_bias = 0.0054, plugin_rmse = 0.0143, plugin_coverage = 0.94),
  K = c(naive_bias = 0.0948, naive_rmse = 0.1069, naive_coverage = 0.29,
        plugin_bias = -0.0154, plugin_rmse = 0.0323, plugin_coverage = 0.90),
  alpha = c(naive_bias = -0.1044, naive_rmse = 0.1093, naive_coverage = 0.04,
            plugin_bias = 0.0202, plugin_rmse = 0.0345, plugin_coverage = 0.90),
  c = c(naive_bias = 0.0076, naive_rmse = 0.0101, naive_coverage = 0.69,
        plugin_bias = 0.0024, plugin_rmse = 0.0047, plugin_coverage = 0.92),
  p = c(naive_bias = 0.0052, naive_rmse = 0.0510, naive_coverage = 0.97,
        plugin_bias = 0.0266, plugin_rmse = 0.0555, plugin_coverage = 0.93)
)
for (parameter_name in names(expected)) {
  values <- expected[[parameter_name]]
  for (model_name in c("naive", "plugin")) {
    for (column in c("bias", "rmse", "coverage")) {
      key <- paste(model_name, column, sep = "_")
      assert_close(
        get_model(model_name, parameter_name, column), values[[key]],
        label = paste(model_name, parameter_name, column)
      )
    }
  }
}
for (parameter_name in c("K", "alpha", "c")) {
  expected_reduction <- c(K = 0.698, alpha = 0.684, c = 0.530)[[parameter_name]]
  expected_improved <- c(K = 0.85, alpha = 0.94, c = 0.86)[[parameter_name]]
  assert_close(
    get_paired(parameter_name, "rmse_reduction_fraction"),
    expected_reduction, tolerance = 5e-4,
    label = paste(parameter_name, "RMSE reduction")
  )
  assert_close(
    get_paired(parameter_name, "proportion_plugin_improved"),
    expected_improved, label = paste(parameter_name, "paired improvement")
  )
}
assert_close(max(synthetic_diagnostics$rhat), 1.038, tolerance = 5e-4,
             label = "maximum synthetic R-hat")
assert_close(min(synthetic_diagnostics$effective_size), 432, tolerance = 0.5,
             label = "minimum synthetic ESS")

# Ridgecrest counts, scores, instability, and diagnostics.
assert_close(completeness$G, 5.5751, label = "Ridgecrest G")
assert_close(completeness$H, 0.7923, label = "Ridgecrest H")
assert_close(completeness$return_to_baseline_days, 0.0159,
             label = "Ridgecrest completeness recovery")
stopifnot(
  all(forecast$n_history_events == 811L),
  all(forecast$n_likelihood_targets == 810L),
  all(forecast$n_test == 224L)
)
assert_close(get_forecast("naive", "log_score_per_event"), 2.4021,
             label = "Naive temporal LPD")
assert_close(get_forecast("plugin", "log_score_per_event"), 2.4294,
             label = "Plug-in temporal LPD")
assert_close(forecast_mc$information_gain_per_event, 0.0273,
             label = "information gain")
assert_close(forecast_mc$mc_q025, 0.0258, label = "information-gain lower bound")
assert_close(forecast_mc$mc_q975, 0.0292, label = "information-gain upper bound")
for (model_name in c("naive", "plugin")) {
  targets <- if (model_name == "naive") {
    c(count_q025 = 68.975, count_median = 176, count_q975 = 752.4,
      overflow_rate = 0.004, probability_supercritical = 0.0174791666666667)
  } else {
    c(count_q025 = 149, count_median = 263, count_q975 = 1501.7,
      overflow_rate = 0.020, probability_supercritical = 0.83971875)
  }
  for (column in names(targets)) {
    assert_close(get_forecast(model_name, column), targets[[column]],
                 label = paste(model_name, column))
  }
}
assert_close(max(ridgecrest_diagnostics$rhat), 1.004, tolerance = 5e-4,
             label = "maximum Ridgecrest R-hat")
assert_close(min(ridgecrest_diagnostics$effective_size), 621, tolerance = 0.6,
             label = "minimum Ridgecrest ESS")

manuscript_files <- c("writing/first_draft_feedback.tex", "writing/refs.bib")
manuscript_available <- all(file.exists(manuscript_files))
if (manuscript_available) {
  # Every citation key and cross-reference used by the manuscript must resolve.
  tex <- paste(readLines(manuscript_files[[1L]], warn = FALSE), collapse = "\n")
  bib <- paste(readLines(manuscript_files[[2L]], warn = FALSE), collapse = "\n")
  citation_matches <- regmatches(
    tex, gregexpr("\\\\cite(?:t|p)?(?:\\[[^]]*\\])?\\{[^}]+\\}", tex, perl = TRUE)
  )[[1L]]
  citation_keys <- unique(trimws(unlist(strsplit(
    sub(".*\\{([^}]+)\\}", "\\1", citation_matches), ",", fixed = TRUE
  ))))
  bib_keys <- regmatches(
    bib, gregexpr("(?<=@\\w{1,20}\\{)[^,]+", bib, perl = TRUE)
  )[[1L]]
  missing_citations <- setdiff(citation_keys, bib_keys)
  if (length(missing_citations)) {
    fail("Citation keys missing from refs.bib: ", paste(missing_citations, collapse = ", "))
  }
  labels <- regmatches(tex, gregexpr("(?<=\\\\label\\{)[^}]+", tex, perl = TRUE))[[1L]]
  references <- regmatches(
    tex, gregexpr("(?<=\\\\(?:eqref|ref)\\{)[^}]+", tex, perl = TRUE)
  )[[1L]]
  missing_labels <- setdiff(unique(references), unique(labels))
  if (length(missing_labels)) {
    fail("Cross-references without labels: ", paste(missing_labels, collapse = ", "))
  }
  if (anyDuplicated(labels)) fail("Duplicate LaTeX labels detected")
}

# Parse every active R source file without executing expensive analyses.
r_files <- c(
  Sys.glob("R/*.R"), Sys.glob("config/*.R"),
  read.csv("results/submission_v1/code_manifest.csv", stringsAsFactors = FALSE)$path
)
for (path in sort(unique(r_files))) {
  tryCatch(parse(path), error = function(error) {
    fail("R parse failure in ", path, ": ", conditionMessage(error))
  })
}

manual_items <- character()
if (manuscript_available) {
  if (grepl("Acknowledgements will be completed", tex, fixed = TRUE)) {
    manual_items <- c(manual_items, "replace the acknowledgements placeholder")
  }
  class_available <- file.exists("writing/statsmsc.cls")
  if (!class_available && nzchar(Sys.which("kpsewhich"))) {
    class_available <- length(suppressWarnings(
      system2("kpsewhich", "statsmsc.cls", stdout = TRUE, stderr = FALSE)
    )) > 0L
  }
  if (!class_available) {
    manual_items <- c(manual_items, "provide statsmsc.cls in the LaTeX environment")
  }
}

cat("Submission verification passed.\n")
cat("Checked analysis claims, figures, and R syntax.\n")
if (manuscript_available) {
  cat("Checked manuscript citations and cross-references.\n")
} else {
  cat("Manuscript checks skipped because writing files are not distributed.\n")
}
if (length(manual_items)) {
  cat("Manual items remaining:\n- ", paste(manual_items, collapse = "\n- "), "\n", sep = "")
}

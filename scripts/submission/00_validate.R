source("config/submission_experiment.R")
source("R/etas_likelihood.R")
source("R/simulate_etas.R")
source("R/loglinear_incompleteness.R")
source("R/fit_loglinear_mle.R")
source("R/fit_bayesian_etas.R")

cfg <- submission_experiment_config()
dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)

theta <- as.list(cfg$theta_true)
branching <- theta$K * cfg$b / (cfg$b - theta$alpha)
stopifnot(abs(branching - cfg$branching_ratio) < 1e-12)

test <- simulate_etas(
  theta, cfg$T1, 10, cfg$M0, cfg$beta,
  Ht = data.frame(ts = 5, magnitudes = 6.5),
  seed = 1
)
test_repeat <- simulate_etas(
  theta, cfg$T1, 10, cfg$M0, cfg$beta,
  Ht = data.frame(ts = 5, magnitudes = 6.5),
  seed = 1
)
stopifnot(identical(test, test_repeat))
children <- which(!is.na(test$parent_id))
if (length(children)) {
  parent_row <- match(test$parent_id[children], test$ID)
  stopifnot(!anyNA(parent_row))
  stopifnot(all(test$ts[parent_row] < test$ts[children]))
  stopifnot(all(test$gen[children] >= 1L))
}
formatted <- make_inlabru_catalog(test, cfg$T1, 10)$catalog_bru
seed_row <- which(test$gen == -1L)
stopifnot(length(seed_row) == 1L)
stopifnot(!formatted$include_likelihood[seed_row])
stopifnot(all(formatted$include_likelihood[-seed_row]))

# Independent analytic-versus-numerical check of the normalised ETAS
# compensator used by both simulation and fitting.
check_ts <- c(0, 0.2, 1)
check_mag <- c(6.5, 3, 4)
grid <- seq(0, 30, length.out = 200000L)
lambda <- etas_lambda_at(grid, check_ts, check_mag, theta, cfg$M0)
numeric_integral <- sum(
  0.5 * (lambda[-1L] + lambda[-length(lambda)]) * diff(grid)
)
analytic_integral <- etas_integrated_intensity(
  0, 30, check_ts, check_mag, theta, cfg$M0
)
integral_relative_error <- abs(numeric_integral - analytic_integral) /
  analytic_integral
stopifnot(integral_relative_error < 0.01)

recovery <- loglinear_recovery_end(
  cfg$mainshock$ts, cfg$mainshock$magnitudes, cfg$M0,
  cfg$incompleteness$G, cfg$incompleteness$H
)
stopifnot(recovery > cfg$mainshock$ts)
recovery_threshold <- compute_loglinear_mct_at(
  recovery,
  cfg$mainshock$ts,
  cfg$mainshock$magnitudes,
  cfg$M0,
  cfg$incompleteness$G,
  cfg$incompleteness$H
)
stopifnot(abs(recovery_threshold - cfg$M0) < 1e-10)

# The hard-threshold catalogue operation and the likelihood's detectable
# fraction must encode the same Mc(t). Seeded events are retained explicitly,
# while every other retained event must satisfy its contemporaneous threshold.
observed_test <- apply_loglinear_incompleteness(
  test,
  cfg$M0,
  cfg$incompleteness$G,
  cfg$incompleteness$H,
  cfg$b,
  trigger_mode = "seeded",
  retain_seeded = TRUE
)
audit <- attr(observed_test, "detection_audit")
non_seed <- audit$gen != -1L
stopifnot(
  all(audit$detected[non_seed] == (
    audit$magnitudes[non_seed] >= audit$Mc_t[non_seed]
  )),
  all(audit$detected[audit$gen == -1L])
)
expected_psi <- 10^(-cfg$b * pmax(audit$Mc_t - cfg$M0, 0))
stopifnot(max(abs(audit$Psi_t - expected_psi)) < 1e-12)

# Event-split Gauss-Legendre quadrature must integrate constants exactly.
quadrature <- make_event_split_quadrature(0, 10, c(1, 1, 3, 8), 100L)
stopifnot(
  all(quadrature$t > 0 & quadrature$t < 10),
  abs(sum(quadrature$weight) - 10) < 1e-12
)

versions <- data.frame(
  package = c(
    "R", "ETAS.inlabru", "inlabru", "INLA", "Rcpp", "coda", "digest"
  ),
  version = c(
    R.version.string,
    as.character(packageVersion("ETAS.inlabru")),
    as.character(packageVersion("inlabru")),
    as.character(packageVersion("INLA")),
    as.character(packageVersion("Rcpp")),
    as.character(packageVersion("coda")),
    as.character(packageVersion("digest"))
  )
)
write.csv(
  versions,
  file.path(cfg$results_dir, "software_versions.csv"),
  row.names = FALSE
)
code_files <- c(
  "config/submission_experiment.R",
  "config/submission_ridgecrest.R",
  "R/etas_likelihood.R",
  "R/simulate_etas.R",
  "R/loglinear_incompleteness.R",
  "R/fit_loglinear_mle.R",
  "R/fit_bayesian_etas.R",
  "R/mcmc_etas.R",
  "R/ridgecrest_mcmc.R",
  "R/submission_workflow.R",
  "R/forecast_observed_etas.R",
  file.path(
    "scripts/submission",
    c(
      sprintf("%02d_%s.R", 0:9, c(
        "validate", "simulate", "prepare_synthetic_initialisation",
        "fit_synthetic_mcmc", "evaluate_synthetic",
        "prepare_ridgecrest", "prepare_ridgecrest_initialisation",
        "fit_ridgecrest_mcmc", "evaluate_ridgecrest",
        "make_thesis_figures"
      )),
      "04_evaluate_synthetic_benchmarks.R",
      "12_fit_synthetic_b_sensitivity.R",
      "13_evaluate_synthetic_b_sensitivity.R",
      "14_fit_ridgecrest_b_sensitivity.R",
      "15_fit_synthetic_prior_sensitivity.R",
      "16_evaluate_synthetic_prior_sensitivity.R",
      "19_make_sensitivity_figures.R"
    )
  ),
  "scripts/submission/internal/fit_baseline_initialisation.R",
  "scripts/submission/internal/fit_plugin_initialisation.R"
)
code_files <- sort(unique(code_files[file.exists(code_files)]))
write.csv(
  data.frame(
    path = code_files,
    md5 = unname(tools::md5sum(code_files))
  ),
  file.path(cfg$results_dir, "code_manifest.csv"),
  row.names = FALSE
)
saveRDS(cfg, file.path(cfg$results_dir, "frozen_config.rds"))

cat("Submission validation passed.\n")
cat("Experiment:", cfg$experiment_id, "\n")
cat("Branching ratio:", branching, "\n")
cat("Compensator relative error:", integral_relative_error, "\n")
cat("Recovery duration (hours):", 24 * (recovery - cfg$mainshock$ts), "\n")
print(versions)

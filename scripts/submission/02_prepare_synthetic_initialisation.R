# Preliminary fits used only to initialise the final synthetic MCMC chains.
source("config/submission_experiment.R")
source("R/submission_workflow.R")
source("R/fit_loglinear_mle.R")

cfg <- submission_experiment_config()
start <- as.integer(Sys.getenv("SIM_START", unset = "1"))
end <- as.integer(Sys.getenv("SIM_END", unset = as.character(cfg$n_sim)))
force <- Sys.getenv("FORCE_REFIT", unset = "0") == "1"
if (start < 1L || end > cfg$n_sim || start > end) {
  stop("Require 1 <= SIM_START <= SIM_END <= N_SIM.")
}

out <- file.path(cfg$results_dir, "initialisation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
detection_env <- c(
  "DETECTION_MODEL=loglinear_mainshock",
  paste0("MAINSHOCK_TIME=", cfg$mainshock$ts),
  paste0("MAINSHOCK_MAGNITUDE=", cfg$mainshock$magnitudes),
  paste0("B_GR=", cfg$b)
)

for (i in seq.int(start, end)) {
  input <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", i))
  if (!file.exists(input)) stop("Missing input: ", input)
  gh_path <- file.path(out, sprintf("plugin_GH_%03d.rds", i))
  gh <- estimate_loglinear_plugin(
    readRDS(input),
    cfg$mainshock$ts,
    cfg$mainshock$magnitudes,
    cfg$M0,
    cfg$b
  )
  saveRDS(gh, gh_path)

  run_submission_fit(
    cfg, "naive", input, out, sprintf("naive_%03d", i), force = force
  )
  run_submission_fit(
    cfg, "plugin", input, out, sprintf("plugin_%03d", i),
    extra_env = c(
      detection_env,
      paste0("LOGLINEAR_G=", gh$estimate[["G"]]),
      paste0("LOGLINEAR_H=", gh$estimate[["H"]])
    ),
    force = force
  )
  cat("Primary pair complete:", i, "\n")
}

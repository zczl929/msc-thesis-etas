# Preliminary fits used only to initialise the final Ridgecrest MCMC chains.
source("config/submission_ridgecrest.R")
source("R/submission_workflow.R")
source("R/fit_loglinear_mle.R")

cfg <- submission_ridgecrest_config()
train <- readRDS(cfg$train_path)
if (!identical(attr(train, "experiment_id"), cfg$experiment_id)) {
  stop("Run 05_prepare_ridgecrest.R first.")
}
cfg$T1 <- min(train$ts) - 1e-6
cfg$T2 <- cfg$train_end
initialisation_dir <- file.path(cfg$results_dir, "initialisation")
dir.create(initialisation_dir, recursive = TRUE, showWarnings = FALSE)

gh <- estimate_loglinear_multitrigger_plugin(
  train, cfg$M0, cfg$M_trigger, cfg$b
)
attr(gh, "experiment_id") <- cfg$experiment_id
attr(gh, "input_md5") <- unname(tools::md5sum(cfg$train_path))
saveRDS(gh, file.path(initialisation_dir, "plugin_GH.rds"))

force <- Sys.getenv("FORCE_REFIT", unset = "0") == "1"
run_submission_fit(
  cfg, "naive", cfg$train_path, initialisation_dir, "naive", force = force
)
run_submission_fit(
  cfg, "plugin", cfg$train_path, initialisation_dir, "plugin",
  extra_env = c(
    "DETECTION_MODEL=loglinear_multitrigger",
    paste0("M_TRIGGER=", cfg$M_trigger),
    paste0("B_GR=", cfg$b),
    paste0("LOGLINEAR_G=", gh$estimate[["G"]]),
    paste0("LOGLINEAR_H=", gh$estimate[["H"]])
  ),
  force = force
)
print(gh)

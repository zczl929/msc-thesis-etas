# Internal Naive ETAS approximation used only for MCMC initialisation.

set.seed(2026)
future_limit_gb <- as.numeric(
  Sys.getenv("FUTURE_GLOBALS_MAX_GB", unset = "2")
)
options(future.globals.maxSize = future_limit_gb * 1024^3)
source("R/fit_bayesian_etas.R")

real_data_mode <- Sys.getenv("REAL_DATA", unset = "0") == "1"

required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop(name, " must be set by an active workflow script.")
  value
}

prior_spec <- default_etas_prior_spec()
if (nzchar(Sys.getenv("MU_PRIOR_SHAPE", unset = ""))) {
  prior_spec$a_mu <- as.numeric(Sys.getenv("MU_PRIOR_SHAPE"))
}
if (nzchar(Sys.getenv("MU_PRIOR_RATE", unset = ""))) {
  prior_spec$b_mu <- as.numeric(Sys.getenv("MU_PRIOR_RATE"))
}
if (nzchar(Sys.getenv("MU_INIT", unset = ""))) {
  prior_spec$mu_init <- as.numeric(Sys.getenv("MU_INIT"))
}
if (nzchar(Sys.getenv("K_PRIOR_LOGMEAN", unset = ""))) {
  prior_spec$a_K <- as.numeric(Sys.getenv("K_PRIOR_LOGMEAN"))
}
if (nzchar(Sys.getenv("K_PRIOR_LOGSD", unset = ""))) {
  prior_spec$b_K <- as.numeric(Sys.getenv("K_PRIOR_LOGSD"))
}
if (nzchar(Sys.getenv("K_INIT", unset = ""))) {
  prior_spec$K_init <- as.numeric(Sys.getenv("K_INIT"))
}
for (entry in list(
  c("a_alpha", "ALPHA_PRIOR_LOWER"),
  c("b_alpha", "ALPHA_PRIOR_UPPER"),
  c("alpha_init", "ALPHA_INIT"),
  c("a_c", "C_PRIOR_LOWER"),
  c("b_c", "C_PRIOR_UPPER"),
  c("c_init", "C_INIT"),
  c("a_p", "P_PRIOR_LOWER"),
  c("b_p", "P_PRIOR_UPPER"),
  c("p_init", "P_INIT")
)) {
  value <- Sys.getenv(entry[[2]], unset = "")
  if (nzchar(value)) prior_spec[[entry[[1]]]] <- as.numeric(value)
}

config <- list(
  experiment_id = Sys.getenv("EXPERIMENT_ID", unset = "unspecified"),
  input_path = required_env("INPUT_PATH"),
  output_dir = required_env("OUTPUT_DIR"),
  output_stem = required_env("OUTPUT_STEM"),
  resume_path = Sys.getenv("RESUME_PATH", unset = ""),
  M0 = as.numeric(Sys.getenv("M0", unset = "2.5")),
  T1 = as.numeric(Sys.getenv("T1", unset = "0")),
  T2 = as.numeric(Sys.getenv("T2", unset = "30")),
  true_values = if (real_data_mode) {
    NULL
  } else {
    c(
      mu = as.numeric(Sys.getenv("TRUE_MU", unset = "0.2")),
      K = as.numeric(Sys.getenv("TRUE_K", unset = "0.4")),
      alpha = as.numeric(Sys.getenv("TRUE_ALPHA", unset = "0.7")),
      c = as.numeric(Sys.getenv("TRUE_C", unset = "0.01")),
      p = as.numeric(Sys.getenv("TRUE_P", unset = "1.2"))
    )
  },
  prior_spec = prior_spec,
  bru_verbose = as.integer(Sys.getenv("BRU_VERBOSE", unset = "3")),
  bru_max_iter = as.integer(Sys.getenv("BRU_MAX_ITER", unset = "20")),
  bru_rel_tol = as.numeric(Sys.getenv("BRU_REL_TOL", unset = "0.1")),
  coef_t = as.numeric(Sys.getenv("COEF_T", unset = "0.15")),
  delta_t = as.numeric(Sys.getenv("DELTA_T", unset = "0.001")),
  N_max = as.integer(Sys.getenv("N_MAX", unset = "250")),
  n_draws = as.integer(Sys.getenv("N_DRAWS", unset = "2000")),
  max_batch = as.integer(Sys.getenv("MAX_BATCH", unset = "1000"))
)

fit_out <- fit_inlabru_baseline(config)

print(round(fit_out$result$posterior_summary, 4))
cat("Saved results to:", fit_out$paths$rds_path, "\n")

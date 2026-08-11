# ============================================================
# setup_environment.R
# Install and verify the R environment for the project.
# Run this once per machine. Subsequent scripts assume these
# packages are installed and load-able.
# ============================================================

# ---- 1. CRAN packages -------------------------------------------------------
# Packages required by the final MCMC, data-integrity and initialisation code.
cran_pkgs <- c(
  "remotes",      # install from GitHub
  "coda",         # MCMC diagnostics
  "Rcpp",         # fast ETAS likelihood
  "digest"        # frozen Ridgecrest checksum
)

to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install)
}

# ---- 2. R-INLA (from INLA's own repository, not CRAN) -----------------------
# Using the testing channel as recommended by ETAS.inlabru.
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(
    "INLA",
    repos = c(getOption("repos"),
              INLA = "https://inla.r-inla-download.org/R/testing"),
    dep = TRUE
  )
}

# ---- 3. inlabru -------------------------------------------------------------
if (!requireNamespace("inlabru", quietly = TRUE)) {
  install.packages("inlabru")
}

# ---- 4. ETAS.inlabru (GitHub) -----------------------------------------------
if (!requireNamespace("ETAS.inlabru", quietly = TRUE)) {
  remotes::install_github("edinburgh-seismicity-hub/ETAS.inlabru")
}

# ---- 5. Sanity check: load everything and print versions --------------------
suppressPackageStartupMessages({
  library(INLA)
  library(inlabru)
  library(ETAS.inlabru)
  library(coda)
  library(Rcpp)
  library(digest)
})

cat("\n========== Environment Check ==========\n")
cat(sprintf("R version       : %s\n", R.version.string))
cat(sprintf("INLA version    : %s\n", as.character(inla.version("version"))))
cat(sprintf("inlabru version : %s\n", as.character(packageVersion("inlabru"))))
cat(sprintf("ETAS.inlabru    : %s\n", as.character(packageVersion("ETAS.inlabru"))))
cat(sprintf("Working dir     : %s\n", getwd()))
cat("=======================================\n")

# ---- 6. Minimal INLA smoke test --------------------------------------------
# Fits a trivial model to confirm INLA's backend C++ binary actually runs.
# (Installation can succeed while the binary fails on some macOS setups.)
cat("\nRunning minimal INLA fit (smoke test)...\n")
set.seed(1)
n <- 50
df <- data.frame(x = rnorm(n))
df$y <- 2 + 1.5 * df$x + rnorm(n, sd = 0.3)

fit <- tryCatch(
  inla(y ~ x, data = df, family = "gaussian"),
  error = function(e) {
    cat("INLA smoke test FAILED:\n")
    cat(conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(fit)) {
  cat("INLA smoke test PASSED.\n")
  cat("Posterior mean of intercept:", round(fit$summary.fixed["(Intercept)", "mean"], 3),
      "(true value: 2)\n")
  cat("Posterior mean of x slope  :", round(fit$summary.fixed["x", "mean"], 3),
      "(true value: 1.5)\n")
}

cat("\nSetup check complete.\n")

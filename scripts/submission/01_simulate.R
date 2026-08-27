source("config/submission_experiment.R")
source("R/etas_likelihood.R")
source("R/simulate_etas.R")
source("R/loglinear_incompleteness.R")

cfg <- submission_experiment_config()
dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)
force <- Sys.getenv("FORCE_RESIMULATE", unset = "0") == "1"
rows <- vector("list", cfg$n_sim)

is_current <- function(path, seed, type) {
  if (!file.exists(path)) return(FALSE)
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  identical(attr(x, "experiment_id"), cfg$experiment_id) &&
    identical(attr(x, "simulation_seed"), seed) &&
    identical(attr(x, "catalogue_type"), type)
}

for (i in seq_len(cfg$n_sim)) {
  seed <- cfg$seed_start + i
  complete_path <- file.path(cfg$data_dir, sprintf("complete_%03d.rds", i))
  observed_path <- file.path(cfg$data_dir, sprintf("observed_%03d.rds", i))

  if (!force &&
      is_current(complete_path, seed, "complete") &&
      is_current(observed_path, seed, "observed")) {
    complete <- readRDS(complete_path)
    observed <- readRDS(observed_path)
  } else {
    complete <- simulate_etas(
      theta = as.list(cfg$theta_true),
      T1 = cfg$T1,
      T2 = cfg$T2,
      M_cut = cfg$M0,
      beta = cfg$beta,
      Ht = cfg$mainshock,
      seed = seed,
      max_events = cfg$max_events
    )
    observed <- apply_loglinear_incompleteness(
      complete,
      M0 = cfg$M0,
      G = cfg$incompleteness$G,
      H = cfg$incompleteness$H,
      b = cfg$b,
      trigger_mode = cfg$incompleteness$trigger_mode,
      M_trigger = cfg$incompleteness$M_trigger,
      cap_at_trigger = cfg$incompleteness$cap_at_trigger,
      retain_seeded = TRUE
    )
    attr(complete, "experiment_id") <- cfg$experiment_id
    attr(complete, "simulation_seed") <- seed
    attr(complete, "catalogue_type") <- "complete"
    attr(observed, "experiment_id") <- cfg$experiment_id
    attr(observed, "simulation_seed") <- seed
    attr(observed, "catalogue_type") <- "observed"
    saveRDS(complete, complete_path)
    saveRDS(observed, observed_path)
  }

  stopifnot(sum(complete$gen == -1L) == 1L)
  stopifnot(sum(observed$gen == -1L) == 1L)
  stopifnot(all(observed$magnitudes >= observed$Mc_t))
  rows[[i]] <- data.frame(
    sim_id = i,
    seed = seed,
    complete_md5 = unname(tools::md5sum(complete_path)),
    observed_md5 = unname(tools::md5sum(observed_path)),
    n_complete = nrow(complete),
    n_observed = nrow(observed),
    observed_fraction = nrow(observed) / nrow(complete),
    first_day_complete = sum(
      complete$ts > cfg$mainshock$ts &
        complete$ts <= cfg$mainshock$ts + 1
    ),
    first_day_observed = sum(
      observed$ts > cfg$mainshock$ts &
        observed$ts <= cfg$mainshock$ts + 1
    )
  )
}

manifest <- do.call(rbind, rows)
write.csv(manifest, file.path(cfg$data_dir, "catalogue_manifest.csv"), row.names = FALSE)
cat("Generated/verified", nrow(manifest), "paired catalogues.\n")
print(summary(manifest[c("n_complete", "n_observed", "observed_fraction")]))

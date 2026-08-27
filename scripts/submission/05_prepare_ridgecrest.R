# Freeze and verify the Ridgecrest training and test catalogues.
source("config/submission_ridgecrest.R")

cfg <- submission_ridgecrest_config()
if (!file.exists(cfg$raw_path)) stop("Missing Ridgecrest source catalogue.")
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required to validate the Ridgecrest source.")
}
raw_sha256 <- digest::digest(
  file = cfg$raw_path, algo = "sha256", serialize = FALSE
)
if (!identical(raw_sha256, cfg$raw_sha256)) {
  stop("Ridgecrest source SHA-256 does not match the frozen query download.")
}
x <- read.csv(cfg$raw_path, stringsAsFactors = FALSE)
required <- c("id", "time", "latitude", "longitude", "mag")
if (!all(required %in% names(x))) {
  stop("Ridgecrest source is missing required columns.")
}
x$time_posix <- as.POSIXct(
  x$time, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"
)
if (anyNA(x$time_posix)) stop("Failed to parse one or more event times.")
if (anyDuplicated(x$id)) stop("Ridgecrest source contains duplicate event IDs.")
main <- x[x$id == cfg$mainshock_id, , drop = FALSE]
if (nrow(main) != 1L) {
  stop("Configured Ridgecrest mainshock ID is not unique in the source.")
}
expected_main_time <- as.POSIXct(
  cfg$mainshock_time_utc, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"
)
query_start <- as.POSIXct(
  cfg$query_start_utc, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"
)
query_end <- as.POSIXct(
  cfg$query_end_utc, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"
)
stopifnot(
  isTRUE(all.equal(main$time_posix, expected_main_time)),
  isTRUE(all.equal(main$latitude, cfg$centre_latitude)),
  isTRUE(all.equal(main$longitude, cfg$centre_longitude)),
  main$mag == 7.1,
  min(x$time_posix) >= query_start,
  max(x$time_posix) <= query_end,
  min(x$mag, na.rm = TRUE) >= cfg$M0
)
rad <- pi / 180
dlat <- (x$latitude - cfg$centre_latitude) * rad
dlon <- (x$longitude - cfg$centre_longitude) * rad
a <- sin(dlat / 2)^2 +
  cos(cfg$centre_latitude * rad) * cos(x$latitude * rad) * sin(dlon / 2)^2
x$distance_km <- 6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
x$ts <- as.numeric(difftime(x$time_posix, main$time_posix, units = "days"))

keep <- x$time_posix >= query_start &
  x$time_posix <= query_end &
  x$distance_km <= cfg$radius_km &
  x$mag >= cfg$M0 &
  x$ts <= cfg$forecast_end
full <- data.frame(
  ID = x$id[keep],
  ts = x$ts[keep],
  magnitudes = x$mag[keep]
)
full <- full[order(full$ts, full$ID), , drop = FALSE]
rownames(full) <- NULL
train <- full[full$ts <= cfg$train_end, , drop = FALSE]
T1 <- min(train$ts) - 1e-6

attr(full, "experiment_id") <- cfg$experiment_id
attr(full, "design") <- list(
  M0 = cfg$M0,
  radius_km = cfg$radius_km,
  train_start = T1,
  train_end = cfg$train_end,
  forecast_end = cfg$forecast_end
)
attr(train, "experiment_id") <- cfg$experiment_id
attr(train, "design") <- attr(full, "design")

dir.create(dirname(cfg$full_path), recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(full, cfg$full_path)
saveRDS(train, cfg$train_path)
manifest <- data.frame(
  source_name = cfg$source_name,
  source_url = cfg$source_url,
  retrieved_at_utc = cfg$retrieved_at_utc,
  query_start_utc = cfg$query_start_utc,
  query_end_utc = cfg$query_end_utc,
  query_centre_latitude = cfg$query_centre_latitude,
  query_centre_longitude = cfg$query_centre_longitude,
  centre_latitude = cfg$centre_latitude,
  centre_longitude = cfg$centre_longitude,
  radius_km = cfg$radius_km,
  minimum_magnitude = cfg$M0,
  raw_sha256 = raw_sha256,
  raw_md5 = unname(tools::md5sum(cfg$raw_path)),
  full_md5 = unname(tools::md5sum(cfg$full_path)),
  train_md5 = unname(tools::md5sum(cfg$train_path)),
  n_source_events = nrow(x),
  n_filtered_events = nrow(full),
  mainshock_id = main$id,
  mainshock_magnitude = main$mag,
  train_start = T1,
  n_pre_mainshock = sum(train$ts < 0),
  n_post_mainshock_train = sum(train$ts >= 0),
  n_holdout = sum(full$ts > cfg$train_end)
)
write.csv(
  manifest,
  file.path(cfg$results_dir, "catalogue_manifest.csv"),
  row.names = FALSE
)
print(manifest)

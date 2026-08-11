source("config/submission_experiment.R")

submission_ridgecrest_config <- function() {
  primary <- submission_experiment_config()
  list(
    experiment_id = "submission-v1-ridgecrest-strict-inla",
    theta_true = NULL,
    source_name = "USGS ComCat FDSN Event Web Service",
    source_url = paste0(
      "https://earthquake.usgs.gov/fdsnws/event/1/query.csv?",
      "starttime=2019-07-04T00%3A00%3A00Z&",
      "endtime=2019-07-16T03%3A19%3A53.040Z&",
      "latitude=35.7695&longitude=-117.5993&",
      "maxradiuskm=75&minmagnitude=2.95&orderby=time-asc"
    ),
    retrieved_at_utc = "2026-07-25T04:06:47Z",
    raw_path = "data/real/ridgecrest_usgs_exact.csv",
    raw_sha256 = "e4bec8199e7dfdf83856c5c8fabf7c79d9e7a45126f35a96d311b179af1e1923",
    full_path = "data/submission_v1/ridgecrest_full.rds",
    train_path = "data/submission_v1/ridgecrest_train.rds",
    results_dir = "results/submission_v1/ridgecrest",
    query_start_utc = "2019-07-04T00:00:00Z",
    query_end_utc = "2019-07-16T03:19:53.040Z",
    mainshock_id = "ci38457511",
    mainshock_time_utc = "2019-07-06T03:19:53.040Z",
    query_centre_latitude = 35.7695,
    query_centre_longitude = -117.5993,
    centre_latitude = 35.7695,
    centre_longitude = -117.5993333,
    M0 = 2.95,
    b = 1,
    radius_km = 75,
    T1 = NA_real_,
    T2 = 1,
    train_end = 1,
    forecast_end = 10,
    M_trigger = 5,
    n_predictive = as.integer(Sys.getenv("N_PREDICTIVE", unset = "1000")),
    primary_prior = primary$primary_prior,
    fit = list(
      n_draws = primary$fit$n_draws,
      max_batch = primary$fit$max_batch,
      # The short, dense Ridgecrest sequence converges more slowly than the
      # synthetic catalogues.  Keep the same 1% criterion, but allow enough
      # iterations for that criterion to be met instead of accepting a
      # max-iteration result.
      max_iter = 200L,
      rel_tol = primary$fit$rel_tol,
      N_max = 150L,
      coef_t = 0.15,
      delta_t = 0.0002
    )
  )
}

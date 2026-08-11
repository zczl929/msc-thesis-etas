# Frozen, submission-facing experiment specification.
#
# The primary prior is fixed independently of the simulation truth and follows
# the broad ETAS.inlabru defaults used in Naylor et al. (2023).  Truth values
# are used only to generate data and score recovery.

submission_experiment_config <- function() {
  b <- 1
  alpha <- 0.7
  branching_ratio <- 0.7
  K <- branching_ratio * (b - alpha) / b

  list(
    experiment_id = "submission-v1-conditioned-mainshock-strict-inla",
    n_sim = as.integer(Sys.getenv("N_SIM", unset = "100")),
    seed_start = 51000L,
    data_dir = "data/submission_v1",
    results_dir = "results/submission_v1",
    theta_true = c(mu = 0.1, K = K, alpha = alpha, c = 0.01, p = 1.2),
    branching_ratio = branching_ratio,
    b = b,
    beta = b * log(10),
    T1 = 0,
    T2 = 1500,
    M0 = 2.5,
    max_events = 200000L,
    mainshock = data.frame(ts = 500, magnitudes = 6.5),
    incompleteness = list(
      G = 4.5,
      H = 0.75,
      trigger_mode = "seeded",
      M_trigger = 5.5,
      cap_at_trigger = FALSE
    ),
    fit = list(
      n_draws = as.integer(Sys.getenv("N_DRAWS", unset = "2000")),
      max_batch = 500L,
      max_iter = 80L,
      rel_tol = 0.01,
      N_max = 200L,
      coef_t = 0.15,
      delta_t = 0.0005
    ),
    primary_prior = list(
      a_mu = 0.5,
      b_mu = 0.5,
      a_K = -1,
      b_K = 1,
      a_alpha = 0.05,
      b_alpha = 2,
      a_c = 0.001,
      b_c = 0.2,
      a_p = 1.01,
      b_p = 2.5,
      mu_init = 0.5,
      K_init = exp(-1),
      alpha_init = 0.8,
      c_init = 0.02,
      p_init = 1.3
    ),
    sensitivity_prior = list(
      a_mu = 1,
      b_mu = 2,
      a_K = log(0.1),
      b_K = 1.5,
      a_alpha = 0.05,
      b_alpha = 2,
      a_c = 0.001,
      b_c = 0.2,
      a_p = 1.01,
      b_p = 2.5,
      mu_init = 0.25,
      K_init = 0.1,
      alpha_init = 0.8,
      c_init = 0.02,
      p_init = 1.3
    )
  )
}

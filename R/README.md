# Final R code map

## Scientific model

- `etas_likelihood.R`: normalised temporal ETAS intensity, productivity,
  Omori kernel, and compensator.
- `simulate_etas.R`: branching-process simulation of complete catalogues.
- `loglinear_incompleteness.R`: \(M_c(t)\), detectable fraction, and hard
  threshold used to create incomplete catalogues.
- `fit_loglinear_mle.R`: conditional magnitude-likelihood estimate of
  \(G,H\), plus numerical-integration helpers.

## Bayesian fitting

- `mcmc_etas.R`: parameter transformations, priors, Naive/Plug-in posterior
  targets, and adaptive random-walk Metropolis sampler.
- `ridgecrest_mcmc.R`: Ridgecrest Plug-in posterior target.

## Prediction

- `forecast_observed_etas.R`: test-period temporal log predictive density and
  future observed-history catalogue simulation.

## Computational setup only

- `fit_bayesian_etas.R`, `submission_workflow.R`: preliminary approximation
  and integrity helpers used to initialise MCMC.

The scientific posterior summaries come from `mcmc_etas.R` and
`ridgecrest_mcmc.R`, not from the preliminary approximation.

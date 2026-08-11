# What to understand in the final analysis

## 1. The research logic

The thesis does not compare many competing scientific models. It compares two
ways of fitting the same temporal ETAS structure to an incompletely recorded
catalogue:

- **Naive ETAS** treats the recorded catalogue as complete above \(M_0\).
- **Plug-in ETAS** first estimates the completeness curve \(G,H\), fixes it,
  and multiplies the observed-history ETAS intensity by the detectable
  fraction \(\psi(t)\).

There are then two evaluations:

1. **Synthetic evaluation:** truth is known, so bias, RMSE, and empirical
   credible-interval coverage can be measured.
2. **Ridgecrest evaluation:** truth is unknown, so the fitted Bayesian
   posterior is evaluated through prediction on data not used for fitting.

Ridgecrest is a complementary posterior-predictive analysis, not an external
validation of simulation coverage.

## 2. The model in four equations

The temporal ETAS intensity is

\[
\lambda(t)=\mu+\sum_{t_i<t}
K10^{\alpha(M_i-M_0)}
\frac{p-1}{c}\left(1+\frac{t-t_i}{c}\right)^{-p}.
\]

After a fixed mainshock \((T_m,M_m)\), the working completeness curve is

\[
M_c(t)=\max\{M_0,M_m-G-H\log_{10}(t-T_m)\}.
\]

Under the Gutenberg--Richter magnitude model, the detectable fraction is

\[
\psi(t;G,H)=10^{-b\{M_c(t)-M_0\}}.
\]

The Plug-in observed-event intensity is

\[
\lambda_{\mathrm{plug}}(t)
=\psi(t;\widehat G,\widehat H)
\lambda_{\mathrm{ETAS}}(t\mid\mathcal H_t^{\mathrm{obs}}).
\]

The last expression is an approximation because the ETAS sum contains only
recorded events. Missing parents and their offspring are not reconstructed.

## 3. How one synthetic replicate is created

`01_simulate.R` performs the following:

1. Fix the true values
   \((\mu,K,\alpha,c,p)=(0.10,0.21,0.70,0.01,1.20)\).
2. Seed an \(M6.5\) mainshock at day 500.
3. Simulate a complete ETAS catalogue on days 0--1500.
4. Apply the hard threshold with \(G=4.5,H=0.75\).
5. Save both the latent complete catalogue and its observed incomplete
   counterpart.

The complete catalogue is used to generate and diagnose the experiment. Both
Naive and Plug-in are fitted to the same incomplete counterpart, making the
comparison paired.

## 4. What is fitted

For each incomplete catalogue:

1. `02_prepare_synthetic_initialisation.R` estimates \(G,H\) from observed
   post-mainshock magnitudes by constrained maximum likelihood.
2. `03_fit_synthetic_mcmc.R` builds either the Naive or Plug-in posterior.
3. Four random-walk Metropolis chains explore that posterior.
4. `04_evaluate_synthetic.R` forms posterior medians and equal-tailed 95%
   credible intervals.

The preliminary INLA/Laplace fit in step 1 supplies only dispersed starting
states and a proposal covariance. The accepted MCMC draws are evaluated using
the explicit likelihood and priors in `R/mcmc_etas.R`.

## 5. What the MCMC code is doing

ETAS parameters have constraints, so MCMC works with unconstrained transformed
parameters:

- \(\mu,K\): log scale;
- \(\alpha,c,p\): logit transformation of their prior bounds.

At each iteration the algorithm proposes a multivariate-normal step centred at
the current transformed value. It accepts or rejects the proposal using the
posterior-density ratio. Warm-up adapts only the overall proposal scale; saved
draws are produced after adaptation stops.

The essential diagnostics are:

- \(\widehat R\) close to 1: chains give compatible distributions;
- effective sample size: number of approximately independent draws;
- acceptance rate: diagnostic of proposal movement, not a measure of model
  quality.

You do not need to explain the C++ implementation line by line. It is a speed
implementation of the likelihood written in the Methods.

## 6. How the synthetic results answer RQ1

`04_evaluate_synthetic.R` calculates:

- **Bias:** average posterior-median error;
- **RMSE:** size of posterior-median error;
- **Coverage:** fraction of the 100 credible intervals containing the known
  truth;
- **paired improvement:** how often Plug-in has smaller absolute error on the
  same catalogue.

The authoritative values are in:

- `results/submission_v1/mcmc_primary/summary/model_summary.csv`
- `results/submission_v1/mcmc_primary/summary/paired_summary.csv`

The main finding is large Plug-in improvement for \(K,\alpha,c\), but not a
uniform improvement for every parameter.

## 7. How Ridgecrest answers RQ2

`05_prepare_ridgecrest.R` verifies the frozen USGS query and creates:

- 811 training-history events through day 1;
- 810 likelihood targets because the given M7.1 mainshock is conditioned on;
- 224 test events during days 1--10.

`07_fit_ridgecrest_mcmc.R` fits Naive and single-mainshock Plug-in models to
the training period. `08_evaluate_ridgecrest.R` then uses posterior draws in
two ways:

1. **Temporal LPD:** scores the observed test event path while propagating
   parameter uncertainty.
2. **Free future simulation:** generates a future catalogue without inserting
   the actual test events, producing a posterior predictive distribution for
   total event count.

The authoritative output is:

`results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/forecast_summary.csv`.

The Plug-in temporal score is better and its count median is closer to 224,
but its count distribution is wider with a heavier upper tail. This supports
a mixed, not uniformly superior, forecast conclusion.

## 8. What to be able to explain orally

You should be able to explain:

1. why simulation is required to measure bias and coverage;
2. why Ridgecrest can assess prediction but not parameter truth;
3. how \(M_c(t)\), \(\psi(t)\), Naive, and Plug-in differ;
4. why fixing \(\widehat G,\widehat H\) is called Plug-in;
5. why the Plug-in model is not a latent-event ETAS model;
6. what MCMC posterior draws, \(\widehat R\), and ESS mean;
7. what RMSE, coverage, LPD, and predictive counts measure;
8. why the conclusion is strong for selected synthetic parameters but mixed
   for real-data prediction.

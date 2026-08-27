# Final analysis workflow

Run all commands from the repository root. Scripts are ordered by the data
flow rather than by the history in which they were developed.

## Synthetic experiment

```sh
Rscript scripts/submission/00_validate.R
Rscript scripts/submission/01_simulate.R
Rscript scripts/submission/02_prepare_synthetic_initialisation.R
Rscript scripts/submission/03_fit_synthetic_mcmc.R
Rscript scripts/submission/04_evaluate_synthetic.R
```

Steps 1--2 define and generate the 100 paired complete/incomplete catalogues.
Step 3 estimates catalogue-specific \(G,H\) and creates preliminary
approximate fits used only for MCMC starting values and proposal covariance.
Step 4 targets the Naive and Plug-in posteriors by four-chain MCMC. Step 5
calculates bias, RMSE, empirical 95% coverage, and convergence diagnostics.

The simulation-only benchmark decomposition adds two further exact MCMC
targets. Oracle Plug-in fits the incomplete catalogue using the known
data-generating completeness parameters, while Complete-data ETAS fits the
original catalogue before incompleteness is imposed. They can reuse the pilot
initialisation from replicate 1; this affects only starting values and proposal
covariance, not either posterior target.

```sh
MCMC_MODELS=oracle,complete \
Rscript scripts/submission/03_fit_synthetic_mcmc.R
Rscript scripts/submission/04_evaluate_synthetic_benchmarks.R
```

The benchmark evaluator combines these fits with the frozen Naive and Plug-in
summaries, and writes separate `benchmark_*.csv` files without replacing the
primary two-model evaluation.

The five Naive fits that required longer chains were generated with:

```sh
MCMC_IDS=15,39,42,43,63 \
MCMC_MODELS=naive \
MCMC_ITER=24000 \
MCMC_WARMUP=8000 \
MCMC_FORCE=true \
Rscript scripts/submission/03_fit_synthetic_mcmc.R
```

The final acceptance gate is parameter-wise \(\widehat R\leq1.05\) and
effective sample size at least 400.

## Ridgecrest analysis

```sh
Rscript scripts/submission/05_prepare_ridgecrest.R
Rscript scripts/submission/06_prepare_ridgecrest_initialisation.R

# Final Plug-in fit: INLA/Laplace covariance.
CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
MCMC_MODELS=plugin \
MCMC_ITER=48000 \
MCMC_WARMUP=16000 \
MCMC_THIN=4 \
MCMC_FORCE=true \
Rscript scripts/submission/07_fit_ridgecrest_mcmc.R

# Exact-MCMC pilot used to recalibrate the final Naive proposal covariance.
CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
MCMC_MODELS=naive \
MCMC_ITER=24000 \
MCMC_WARMUP=8000 \
MCMC_THIN=4 \
MCMC_FORCE=true \
Rscript scripts/submission/07_fit_ridgecrest_mcmc.R

# Final Naive fit: recalibrate from the pilot above, then run the longer
# chain required by its K and alpha effective sample sizes.
CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
MCMC_MODELS=naive \
MCMC_ITER=64000 \
MCMC_WARMUP=16000 \
MCMC_THIN=4 \
MCMC_FORCE=true \
MCMC_RECALIBRATE=true \
Rscript scripts/submission/07_fit_ridgecrest_mcmc.R

CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
N_PREDICTIVE=1000 \
Rscript scripts/submission/08_evaluate_ridgecrest.R
```

The M7.1 mainshock is retained in the triggering history but excluded from
likelihood targets. The completeness curve is anchored at that mainshock
only. Days 1--10 are not used to fit ETAS parameters or \(G,H\).

## Outputs

- `results/submission_v1/mcmc_primary/summary/model_summary.csv`
- `results/submission_v1/mcmc_primary/summary/paired_summary.csv`
- `results/submission_v1/mcmc_primary/summary/mcmc_diagnostics.csv`
- `results/submission_v1/mcmc_primary/summary/benchmark_decomposition.csv`
- `results/submission_v1/mcmc_primary/summary/benchmark_model_summary.csv`
- `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/forecast_summary.csv`
- `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/forecast_count_draws.csv`
- `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/mcmc_diagnostics.csv`

Existing compatible fit objects are reused unless a force flag is supplied.
Do not force a full refit unless the model, priors, or data-generating
mechanism has intentionally changed.

## Supplementary sensitivity analyses

The fixed-\(b\) synthetic sensitivity re-estimates \(G,H\) and refits the
Plug-in model for the two non-primary assumed values.  The primary \(b=1\)
summaries are reused by the evaluator.

```sh
ASSUMED_B=0.8 Rscript scripts/submission/12_fit_synthetic_b_sensitivity.R
ASSUMED_B=1.2 Rscript scripts/submission/12_fit_synthetic_b_sensitivity.R
Rscript scripts/submission/13_evaluate_synthetic_b_sensitivity.R
```

For Ridgecrest, the Naive temporal posterior is invariant to the assumed
\(b\)-value.  Each alternative Plug-in posterior is paired with the primary
Naive posterior when regenerating its forecast summaries.

```sh
ASSUMED_B=0.8 Rscript scripts/submission/14_fit_ridgecrest_b_sensitivity.R
ASSUMED_B=0.8 \
CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
NAIVE_FIT_PATH=results/submission_v1/ridgecrest/mcmc_conditioned_single/naive.rds \
PLUGIN_FIT_PATH=results/submission_v1/ridgecrest/b_sensitivity/b_0p8/plugin.rds \
FORECAST_OUT=results/submission_v1/ridgecrest/b_sensitivity/b_0p8/forecast \
Rscript scripts/submission/08_evaluate_ridgecrest.R

ASSUMED_B=1.2 Rscript scripts/submission/14_fit_ridgecrest_b_sensitivity.R
ASSUMED_B=1.2 \
CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
NAIVE_FIT_PATH=results/submission_v1/ridgecrest/mcmc_conditioned_single/naive.rds \
PLUGIN_FIT_PATH=results/submission_v1/ridgecrest/b_sensitivity/b_1p2/plugin.rds \
FORECAST_OUT=results/submission_v1/ridgecrest/b_sensitivity/b_1p2/forecast \
Rscript scripts/submission/08_evaluate_ridgecrest.R
```

The targeted prior sensitivity refits both primary synthetic models under the
alternative priors and then compares the resulting repeated-sampling metrics.

```sh
Rscript scripts/submission/15_fit_synthetic_prior_sensitivity.R
Rscript scripts/submission/16_evaluate_synthetic_prior_sensitivity.R
```

Files under `scripts/submission/internal/` are called automatically by steps 2
and 6 to prepare MCMC initialisation. They are not separate analyses and should
not be run directly.

## Figures and final verification

Generate the final figures and then run the read-only submission check:

```sh
Rscript scripts/submission/09_make_thesis_figures.R
Rscript scripts/submission/19_make_sensitivity_figures.R
Rscript scripts/submission/10_verify_submission.R
```

The verification script checks numerical values against the authoritative CSV
files, verifies required figures, and parses all active R code. If the
separately maintained manuscript and bibliography are present locally, it also
checks citations and cross-references. It does not rerun MCMC.

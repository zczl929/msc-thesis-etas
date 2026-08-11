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

CONDITION_MAINSHOCK=true \
COMPLETENESS_MODE=single \
MCMC_ITER=48000 \
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
- `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/forecast_summary.csv`
- `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/mcmc_diagnostics.csv`

Existing compatible fit objects are reused unless a force flag is supplied.
Do not force a full refit unless the model, priors, or data-generating
mechanism has intentionally changed.

Files under `scripts/submission/internal/` are called automatically by steps 2
and 6 to prepare MCMC initialisation. They are not separate analyses and should
not be run directly.

## Figures and final verification

Generate the final figures and then run the read-only submission check:

```sh
Rscript scripts/submission/09_make_thesis_figures.R
Rscript scripts/submission/10_verify_submission.R
```

The verification script checks manuscript values against the authoritative
CSV files, verifies citations, cross-references and required figures, and
parses all active R code. It does not rerun MCMC.

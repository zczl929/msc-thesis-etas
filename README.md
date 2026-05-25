# MSc Thesis: Parameter Bias in ETAS Models under Short-Term Incompleteness

A Bayesian study of how short-term aftershock incompleteness affects parameter
recovery and credible-interval coverage in temporal ETAS models.

Seojeong Hong · MSc Statistics, Imperial College London · 2026

See [`docs/proposal.pdf`](docs/proposal.pdf) for the full proposal.

## Research Questions

1. To what extent does short-term aftershock incompleteness bias ETAS parameter estimates?
2. Can modelling time-varying completeness within a Bayesian ETAS framework improve
   parameter recovery and uncertainty calibration?

## Repository Structure

```
.
├── R/              # Reusable functions (sourced by scripts)
├── scripts/        # Numbered scripts that run experiments end-to-end
├── notebooks/      # Exploratory / visualisation .Rmd files
├── data/
│   ├── raw/        # Original SCEDC catalogues (not tracked)
│   └── processed/  # Cleaned catalogues for analysis
├── results/
│   ├── simulations/  # Synthetic ETAS catalogues (.rds)
│   ├── fits/         # Fitted inlabru model objects (.rds)
│   ├── figures/      # Plots for report
│   └── tables/       # Summary tables (bias, RMSE, coverage)
├── logs/           # Simulation run logs
└── docs/           # Proposal, notes, references
```

## Reproducibility

Built in R using `inlabru` / `R-INLA` and the
[`ETAS.inlabru`](https://github.com/edinburgh-seismicity-hub/ETAS.inlabru) package.

1. Clone this repository.
2. Run `scripts/00_setup_check.R` to install dependencies.
3. Run scripts in numerical order.

## Data

- **Synthetic catalogues:** Generated via Ogata thinning (`R/simulate_etas.R`).
- **Empirical case study:** 2019 Ridgecrest sequence from the
  [SCEDC catalogue](https://scedc.caltech.edu/). Raw files not committed.

## Status

Work in progress (started Jan 2026).
# MSc Thesis: Parameter Bias in ETAS Models under Short-Term Incompleteness

A Bayesian study of how short-term aftershock incompleteness affects parameter
recovery and credible-interval coverage in temporal ETAS models.

**Author:** Seojeong Hong (CID 06058915)
**Programme:** MSc Statistics (Data Science and Machine Learning), Imperial College London
**Module:** MATH70088 — Research Project
**Supervisor:** [TBD]

## Research Questions

1. To what extent does short-term aftershock incompleteness bias ETAS parameter estimates?
2. Can modelling time-varying completeness within a Bayesian ETAS framework
   improve parameter recovery and uncertainty calibration (frequentist coverage
   of Bayesian credible intervals)?

See `docs/proposal.pdf` for the full proposal.

## Repository Structure

.
├── R/              # Reusable functions (sourced by scripts)
├── scripts/        # Numbered scripts that run the experiments end-to-end
├── notebooks/      # Exploratory / visualisation .Rmd files
├── data/
│   ├── raw/        # Original SCEDC catalogues (not tracked; see below)
│   └── processed/  # Cleaned catalogues for analysis
├── results/
│   ├── simulations/  # Synthetic ETAS catalogues (.rds)
│   ├── fits/         # Fitted inlabru model objects (.rds)
│   ├── figures/      # Plots for report
│   └── tables/       # Summary tables (bias, RMSE, coverage)
├── logs/           # Simulation run logs
└── docs/           # Proposal, notes, references

## Reproducibility

Built in R using `inlabru` / `R-INLA` and the
[`ETAS.inlabru`](https://github.com/edinburgh-seismicity-hub/ETAS.inlabru) package.

To reproduce:
1. Clone this repository.
2. Run `scripts/00_setup_check.R` to install dependencies and verify the environment.
3. Run scripts in numerical order.

Package versions are pinned via `renv` (see `renv.lock`).

## Data Sources

- **Synthetic catalogues:** Generated via Ogata thinning. Code in `R/simulate_etas.R`.
- **Empirical case study:** 2019 Ridgecrest sequence from the
  [SCEDC catalogue](https://scedc.caltech.edu/). Download instructions in
  `scripts/06_ridgecrest_case.R`. Raw files are not committed.

## Status

Work in progress (started May 2026).

## License

Thesis text: All rights reserved.
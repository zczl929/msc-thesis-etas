# Bayesian ETAS under short-term catalogue incompleteness

Code, data, and analysis outputs for an MSc dissertation comparing
Naive ETAS with a Plug-in correction for short-term earthquake catalogue
incompleteness. The repository contains the synthetic simulation study and
the 2019 Ridgecrest case study. The dissertation source and bibliography are
maintained separately.

## Repository structure

```text
R/                      Model, simulation, MCMC, and evaluation functions
config/                 Analysis settings
data/real/              Ridgecrest catalogue
data/submission_v1/     Synthetic catalogue manifest
scripts/submission/     Analysis scripts in execution order
results/submission_v1/  Final CSV summaries and dissertation figures
```

## Running the analysis

Run commands from the repository root. Install/check the R environment once:

```sh
Rscript scripts/setup_environment.R
```

The main workflow follows the numbered scripts in `scripts/submission/`:

```text
00                         validate inputs and configuration
01--04                     synthetic simulation, MCMC, and evaluation
05--08                     Ridgecrest preparation, MCMC, and prediction
09                         main dissertation figures
12--16 and 19              supplementary sensitivity analyses and figures
```

Analysis settings, including MCMC chain lengths and sensitivity specifications,
are defined in the configuration files and analysis scripts.

## Main outputs

- Synthetic summaries: `results/submission_v1/mcmc_primary/summary/`
- Ridgecrest summaries: `results/submission_v1/ridgecrest/`
- Final figures: `results/submission_v1/figures/thesis/`

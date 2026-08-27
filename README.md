# MSc thesis: Bayesian ETAS under short-term incompleteness

This repository contains the final analysis for two related questions:

1. In controlled simulations, does a Plug-in short-term completeness
   correction improve ETAS parameter recovery and 95% credible-interval
   coverage relative to Naive ETAS?
2. In the Ridgecrest sequence, does the same correction improve Bayesian
   posterior prediction on a fixed test period?

The two primary models are fitted by four-chain MCMC. Preliminary
INLA/Laplace fits are retained only to initialise and tune the MCMC proposal;
they are not a scientific comparison and do not contribute posterior draws to
the reported results.

## Start here

- [`notes/analysis_guide.md`](notes/analysis_guide.md): conceptual and code
  guide for the thesis analysis.
- [`scripts/submission/README.md`](scripts/submission/README.md): exact
  execution order.

## Active project structure

```text
config/                 Frozen synthetic and Ridgecrest settings
R/                      Functions used by the final workflow
scripts/submission/     Ordered analysis, figure, and verification scripts
data/submission_v1/     Manifest for generated synthetic catalogues
data/real/              Frozen USGS source catalogue and provenance
results/submission_v1/  Final summaries and dissertation figures
```

The dissertation source and bibliography are maintained separately and are
not distributed in this code-and-results repository.

Only `scripts/submission/` defines the final reproduction path. PETAI,
Modular, multi-trigger, earlier INLA coverage, and rejected model experiments
are not part of the submitted analysis.

`scripts/setup_environment.R` is a one-time package installation and
environment check; it is not an experiment.

## Authoritative results

- Synthetic summaries:
  `results/submission_v1/mcmc_primary/summary/`
- Ridgecrest summaries:
  `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/`

The CSV summary files are the values used in the dissertation. Regeneratable
catalogues, preliminary fits, and posterior-draw objects are intentionally
excluded from the compact submission; the ordered workflow recreates them.

## Claim boundary

The Plug-in model corrects the expected rate of recorded events but uses only
the observed triggering history. It does not impute undetected earthquakes or
restore their triggering contributions. The analysis therefore does not claim
to be a latent-event ETAS model or to outperform PETAI/ETASI.

## Submission verification

The fast verification command checks reported numerical claims against the
authoritative CSV summaries, confirms the required figures exist, and parses
every active R source file:

```sh
Rscript scripts/submission/10_verify_submission.R
```

See [`SUBMISSION_CHECKLIST.md`](SUBMISSION_CHECKLIST.md) for the remaining
manual checks. Generated LaTeX auxiliary files are not part of the code
submission.

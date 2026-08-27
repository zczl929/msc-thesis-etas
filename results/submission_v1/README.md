# Final result map

- `mcmc_primary/summary/`: authoritative synthetic CSV results used in the
  dissertation.
- `ridgecrest/completeness_summary.csv`: fitted single-mainshock completeness
  parameters and recovery time.
- `ridgecrest/mcmc_conditioned_single/forecast/`: authoritative Ridgecrest
  diagnostics, posterior-predictive summaries, and event-count draws used in
  Figure 5.
- `figures/thesis/`: the six figures included in the dissertation body, four
  supplementary figures, and their compact plotting data where applicable.

Regeneratable catalogues, preliminary fits, posterior draws, duplicate PNG
renders, and obsolete result trees are excluded from the compact submission.
The committed PDF figures and compact CSV files are sufficient to inspect the
reported results; omitted computational objects can be recreated by the
ordered scripts in `scripts/submission/`.

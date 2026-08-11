# Submission checklist

Run all commands from the repository root.

## Automated checks

```sh
Rscript scripts/submission/00_validate.R
Rscript scripts/submission/10_verify_submission.R
```

The first command validates the frozen experiment, numerical integration,
catalogue observation rule, software versions, and code checksums. The second
checks every numerical claim in the dissertation against the final CSV
summaries, verifies citations and cross-references, confirms the required
figures, and parses all active R code.

## Manual checks before submission

- Replace the acknowledgements placeholder with an accurate statement.
- Add any required generative-AI disclosure using the programme's prescribed
  wording; do not invent or paraphrase institutional policy.
- Compile `writing/first_draft_feedback.tex` in the official template
  environment containing `statsmsc.cls`.
- Read the compiled PDF once for float placement, page breaks, overfull boxes,
  and unresolved `??` references.
- Confirm the title-page date, candidate ID, and final dissertation filename.
- Confirm that the submitted files exclude `.git`, `.DS_Store`, `Rplots.pdf`,
  and LaTeX auxiliary files.

## Final scientific outputs

- Synthetic results:
  `results/submission_v1/mcmc_primary/summary/`
- Ridgecrest results:
  `results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/`
- Thesis figures:
  `results/submission_v1/figures/thesis/`

The compact submission excludes regeneratable intermediate fits and posterior
draws. The CSV files above are the authoritative values used in the
dissertation; the ordered workflow recreates omitted computational objects.

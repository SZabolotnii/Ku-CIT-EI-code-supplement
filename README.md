# Ku-CIT-EI Code Supplement

Public verification supplement for the manuscript:

**A Transferability Criterion for Null-Optimized Variance Reduction in Cumulant-Based Error-Independence Testing**

Repository: <https://github.com/SZabolotnii/Ku-CIT-EI-code-supplement>

## Scope

The current manuscript evaluates whether a null-optimized PMM2
variance-reduction statistic transfers safely to a cumulant-based
error-independence test. The supplement contains the R code, generated
CSV artifacts, figures, and verification checks needed to reproduce the
headline numerical claims in the manuscript.

The active manuscript scope is:

- W&S-style third-order naive cumulant estimator;
- PMM2 variance-reduced estimator and its loss of alternative-side
  consistency;
- PMM3-style fourth-order diagnostic probe;
- targeted revision checks: Wilson intervals, bootstrap sensitivity,
  Tukey g-and-h heavy-tail sensitivity, and raw distance-covariance
  sanity baseline;
- a real-data empirical illustration on the PHQ-8 depression scale
  (BRFSS 2010, 70 four-plus-four half-splits): the null-optimized
  correction attenuates the statistic in every split and erases a
  one-directional subset of the naive detections;
- PATP only as a conceptual basis-adaptivity direction, not as a
  positive empirical method in this manuscript.

Some older exploratory GSA/DSGE/PATP-DSGE files remain in the repository
for provenance of the broader research track. They are not part of the
current manuscript's headline verification.

## Contents

- `scripts/02_naive_estimators.R`, `scripts/05_dgp.R`,
  `scripts/08_pmm2_estimator.R` - core DGP and estimator helpers.
- `scripts/09_sim_comparative.R` - PMM2 versus naive Monte Carlo driver.
- `scripts/18_pmm3_symmetric_probe.R` - PMM3-style diagnostic probe.
- `scripts/21_revision_experiments.R` - targeted revision experiments.
- `scripts/01_download_brfss2010.R` - downloads and processes the public
  BRFSS 2010 PHQ-8 high-risk subset (regenerates the local `.rds`; raw
  files are not redistributed).
- `scripts/03_phq8_70_splits.R` - naive Delta-c3 baseline over the 70
  PHQ-8 half-splits.
- `scripts/22_phq8_criterion.R` - real-data transferability-criterion
  analysis; writes `output/tables/phq8_criterion.csv` and
  `output/figures/fig_phq8_criterion.{pdf,png}`.
- `scripts/verify_reported_values.R` - fast verification of the
  manuscript-level numerical claims (including the PHQ-8 Section 6
  results) from generated CSV artifacts.
- `output/tables/` - generated CSV tables used for manuscript tables and
  diagnostics.
- `output/figures/` - generated figures from the broader supplement.
- `output/session_info/` - R session snapshots for main follow-up runs.
- `docs/DATA_POLICY.md` - data and redistribution boundary.
- `RUNBOOK.md` - execution-oriented reproduction notes.

## Quick Verification

Run from the repository root:

```sh
Rscript scripts/verify_reported_values.R
```

Expected result:

```text
All current-manuscript verification checks passed.
```

This check uses the generated CSV artifacts and does not rerun the full
Monte Carlo workflow.

## Reproduction Notes

The generated CSV files are included so reviewers can inspect the exact
artifacts used in the manuscript. Full simulation reruns can take
nontrivial time. A focused reproduction path is:

```sh
Rscript scripts/09_sim_comparative.R
Rscript scripts/18_pmm3_symmetric_probe.R
Rscript scripts/21_revision_experiments.R
Rscript scripts/verify_reported_values.R
```

The real-data Section 6 illustration (regenerates the committed
`phq8_criterion.csv` and figure) requires the processed BRFSS subset,
which is produced locally from public CDC data and is not redistributed:

```sh
Rscript scripts/01_download_brfss2010.R   # downloads + processes BRFSS 2010 PHQ-8
Rscript scripts/03_phq8_70_splits.R       # naive baseline (optional cross-check)
Rscript scripts/22_phq8_criterion.R       # criterion table + figure
```

## R Dependencies

The active manuscript scripts use base R plus:

- `moments`
- `here`
- `future`
- `future.apply`
- `ggplot2`
- `dplyr`
- `tidyr`
- `testthat` for local estimator checks

The revision distance-covariance sanity check is implemented locally and
does not require `energy` or `Hmisc`.

## Data Boundary

The manuscript's core results are simulation-based; Section 6 adds a
real-data illustration on the public BRFSS 2010 PHQ-8 data. The processed
subset is regenerated locally by `scripts/01_download_brfss2010.R` from
the public CDC source. This repository ships the generated per-split
criterion table (`output/tables/phq8_criterion.csv`) and figure, but does
not redistribute raw BRFSS files, the processed `.rds`, or third-party
articles.

## License

No explicit open-source license has been declared in this repository
yet. Until a license is added, the code and generated artifacts are
available for review and verification, but reuse rights should be
clarified by the author.

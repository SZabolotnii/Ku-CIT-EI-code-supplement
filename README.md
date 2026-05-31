# Ku-CIT-EI Code Supplement

This repository contains the open verification supplement for the manuscript:

**Kunchenko Stochastic Polynomials for Cumulant-Based Error-Independence Testing**

The supplement is intended to let reviewers and readers inspect the R code, generated Monte Carlo tables, figures, and headline numerical checks used in the manuscript. It intentionally excludes raw BRFSS files, processed `.rds` files, third-party PDFs, and private scratch material.

## Contents

- `scripts/` - R scripts for the W&S-style naive estimators, PMM2 comparison, PMM3 diagnostic, GSA-LLR diagnostics, DSGE/PATP-DSGE hybrid classifiers, figure generation, and headline-value verification.
- `output/tables/` - generated CSV tables used to report Monte Carlo diagnostics, ablations, bootstrap intervals, and verdict checks.
- `output/figures/` - generated figures used in the manuscript and supplement.
- `output/session_info/` - R session snapshots for the main follow-up runs.
- `docs/DATA_POLICY.md` - data and redistribution boundary for this public supplement.
- `RUNBOOK.md` - execution-oriented reproduction notes.

## Quick Verification

To verify the main numerical claims from the included CSV artifacts:

```sh
Rscript scripts/verify_reported_values.R
```

Expected result: all checks pass.

## Reproduction Notes

The included CSV files are the exact generated artifacts used for manuscript tables and diagnostics. The full Monte Carlo scripts are also included, but rerunning every full simulation can take nontrivial time.

Minimal rerun path:

```sh
Rscript scripts/09_sim_comparative.R
Rscript scripts/15_gsa_llr_repair_sweep.R
Rscript scripts/18_pmm3_symmetric_probe.R
Rscript scripts/17_dsge_hybrid_cit_classifier.R
Rscript scripts/19_patp_dsge_hybrid_sweep.R
Rscript scripts/20_make_patp_dsge_hybrid_figure.R
Rscript scripts/verify_reported_values.R
```

Some scripts support quick modes through command-line flags or environment variables; see `RUNBOOK.md`.

## R Dependencies

The core scripts use base R plus:

- `moments`
- `here`
- `future`
- `future.apply`
- `ggplot2`
- `dplyr`
- `tidyr`
- `testthat` for local estimator checks

## Data Boundary

The manuscript is simulation-based. BRFSS 2010 PHQ-8 is used only as public-data context for the Wiedermann-Shi application setting. This repository does not redistribute BRFSS data or third-party articles.

## License

No explicit open-source license has been declared in this repository yet. Until a license is added, the code and generated artifacts are available for review and verification, but reuse rights should be clarified by the author.


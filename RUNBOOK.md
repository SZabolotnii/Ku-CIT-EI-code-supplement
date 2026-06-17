# RUNBOOK

Last updated: 2026-06-17

Scope: public verification supplement for the current Ku-CIT-EI
manuscript, **A Transferability Criterion for Null-Optimized Variance
Reduction in Cumulant-Based Error-Independence Testing**.

## Headline Verification

Run:

```sh
Rscript scripts/verify_reported_values.R
```

This checks the manuscript-level numerical claims against generated CSV
artifacts in `output/tables/`:

- ARE table for PMM2 versus naive;
- bias/attenuation and power-loss table;
- targeted revision diagnostics;
- PMM3-style diagnostic values;
- the real-data PHQ-8 Section 6 results (`phq8_criterion.csv`).

Expected final line:

```text
All current-manuscript verification checks passed.
```

## PMM2 Comparative Monte Carlo

Full run:

```sh
Rscript scripts/09_sim_comparative.R
```

Generated artifacts:

- `output/tables/are_pmm2_vs_naive.csv`
- `output/tables/sim_comparative_summary.csv`
- `data/mc_comparative.rds`

Figure:

```sh
Rscript scripts/10_fig_comparative.R
```

## PMM3 Symmetric Diagnostic Probe

```sh
Rscript scripts/18_pmm3_symmetric_probe.R
```

Expected current result: PMM3 reduces variance but fails the
nuisance-aware testing gate; it is a diagnostic probe, not a new test.

Generated artifacts include:

- `output/tables/pmm3_symmetric_probe_results.csv`
- `output/tables/pmm3_symmetric_probe_verdict.csv`
- `output/tables/revision_pmm3_ci.csv`

## Targeted Revision Experiments

```sh
Rscript scripts/21_revision_experiments.R
```

Generated artifacts:

- `output/tables/revision_ci_existing.csv`
- `output/tables/revision_bootstrap_sensitivity.csv`
- `output/tables/revision_heavytail_sensitivity.csv`
- `output/tables/revision_heavytail_ratios.csv`
- `output/tables/revision_dcov_sanity.csv`
- `output/tables/revision_pmm3_ci.csv`

The script adds:

- Wilson confidence intervals for Type-I and power rates;
- bootstrap sensitivity for `B=200` versus `B=1000` and percentile
  versus BCa intervals;
- a Tukey g-and-h heavy-tail alternative for the skew-heavy confounder;
- a local permutation distance-covariance sanity check.

## PHQ-8 Real-Data Illustration (Section 6)

This regenerates the committed `output/tables/phq8_criterion.csv` and
`output/figures/fig_phq8_criterion.{pdf,png}`. It needs the processed
BRFSS 2010 PHQ-8 subset, which is produced locally from public CDC data
and is not redistributed here.

```sh
Rscript scripts/01_download_brfss2010.R   # download + process BRFSS 2010 PHQ-8 -> data/*.rds
Rscript scripts/03_phq8_70_splits.R       # naive Delta-c3 baseline (optional cross-check)
Rscript scripts/22_phq8_criterion.R       # per-split criterion table + figure
```

Expected current result: the H0-optimal weight estimated on real data is
non-zero (median K* ~ 0.21, matching the simulation value ~ 0.24), the
correction attenuates the statistic toward zero in all 70 splits, and it
erases 6 of the 58 naive detections one-directionally (creating none).
The committed CSV lets `verify_reported_values.R` check these without
the BRFSS `.rds`.

## Legacy Exploratory Branches

The repository still contains older GSA/DSGE/PATP-DSGE exploratory
scripts and generated artifacts from the broader research track. They
are retained for provenance but are not part of the current manuscript's
headline verification.

Do not cite those exploratory positive-candidate results as claims of
the current manuscript.

## Exclusions

This public supplement excludes:

- raw BRFSS files;
- processed `.rds` files from BRFSS;
- third-party source articles or PDFs;
- private scratch files;
- the manuscript source itself.

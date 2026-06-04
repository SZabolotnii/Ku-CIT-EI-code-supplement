# RUNBOOK

Last updated: 2026-06-04

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
- PMM3-style diagnostic values.

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

## Legacy Exploratory Branches

The repository still contains older GSA/DSGE/PATP-DSGE exploratory
scripts and generated artifacts from the broader research track. They
are retained for provenance but are not part of the current manuscript's
headline verification.

Do not cite those exploratory positive-candidate results as claims of
the current JMASM submission.

## Exclusions

This public supplement excludes:

- raw BRFSS files;
- processed `.rds` files from BRFSS;
- third-party source articles or PDFs;
- private scratch files;
- the manuscript source itself.

# Ku-CIT-EI Code Supplement

Public verification supplement for the manuscript:

**A Transferability Criterion for Null-Optimized Variance Reduction in Cumulant-Based Error-Independence Testing**

Repository: <https://github.com/SZabolotnii/Ku-CIT-EI-code-supplement>

Archived versions (Zenodo):

- v1.0.0, the version cited by the SORT revision 1:
  <https://doi.org/10.5281/zenodo.22933610>
- all versions (concept DOI, resolves to the latest):
  <https://doi.org/10.5281/zenodo.22933609>

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
- the null side of the same criterion: when the two idiosyncratic error
  components differ in shape, the auxiliary means are no longer zero and
  the PMM2 statistic loses size under a true null, while the naive
  statistic does not;
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
- `scripts/23_revision_w_shape_size.R` - null-side size experiment: the
  same design with separate skewness for the two error components;
  writes `output/tables/revision_w_shape_size{,_reps}.csv`.
- `scripts/24_revision_kstar_plim.R` - the weight used in Corollary 3:
  probability limits of the sample PMM2 weights in the R = 2 cell and the
  H0 cell of the Section-5 design; writes
  `output/tables/revision_kstar_plim.csv`.
- `scripts/01_download_brfss2010.R` - downloads and processes the public
  BRFSS 2010 PHQ-8 high-risk subset (regenerates the local `.rds`; raw
  files are not redistributed).
- `scripts/03_phq8_70_splits.R` - naive Delta-c3 baseline over the 70
  PHQ-8 half-splits.
- `scripts/22_phq8_criterion.R` - real-data transferability-criterion
  analysis; writes `output/tables/phq8_criterion.csv` and
  `output/figures/fig_phq8_criterion.{pdf,png}`.
- `scripts/22b_phq8_criterion_figure_rev1.R` - redraws that figure for
  the current manuscript from the committed per-split table; recomputes
  nothing.
- `scripts/verify_reported_values.R` - fast verification of the
  manuscript-level numerical claims (including the PHQ-8 Section 6
  results) from generated CSV artifacts.
- `output/tables/` - generated CSV tables used for manuscript tables and
  diagnostics.
- `output/figures/` - generated figures from the broader supplement.
- `output/session_info/` - R session snapshots for main follow-up runs.
- `docs/DATA_POLICY.md` - data and redistribution boundary.
- `RUNBOOK.md` - execution-oriented reproduction notes.

## Manuscript map (SORT revision 1)

Every table, figure and computed number of the revised manuscript, with the
script that produces it and the artifact it is read from. Section and table
numbers are those of the revised manuscript.

| Manuscript item | Script | Artifact in `output/` |
|---|---|---|
| Table 1 (basis functions) | — (definitions) | — |
| §4, remark after Proposition 2; Corollary 3; Appendix A (K* ≈ 0.24, attenuation 84 %) | `24_revision_kstar_plim.R` | `tables/revision_kstar_plim.csv` |
| §5.1, Tukey g-and-h skewness 2.06 / excess kurtosis 14.5 | closed-form g-and-h moments | — |
| §5.1, reference-sample values 2.05 / 14.8 | `21_revision_experiments.R` (object `tgh_ref`) | — |
| §5.2, Table 2 (ARE under H0) | `09_sim_comparative.R`, aggregated by `10_fig_comparative.R` | `tables/are_pmm2_vs_naive.csv` → `tables/sim_comparative_summary.csv` |
| §5.3, Table 3 (bias and attenuation under H1) | `09_sim_comparative.R` | `tables/are_pmm2_vs_naive.csv` |
| §5.4, Table 4 (power, aggregated by (R, γ_T)) | `09_sim_comparative.R`, `10_fig_comparative.R` | `tables/sim_comparative_summary.csv` |
| §5.4, Type-I rates under H0 (twelve cells, 200 replications each) | `09_sim_comparative.R` | `tables/are_pmm2_vs_naive.csv` (rows with R = 1, columns `power_naive`, `power_pmm2`) |
| §5.4, Table 5 and its paragraph (bootstrap sensitivity, g-and-h alternative, dCov) | `21_revision_experiments.R` | `tables/revision_bootstrap_sensitivity.csv`, `revision_heavytail_sensitivity.csv`, `revision_heavytail_ratios.csv`, `revision_dcov_sanity.csv` |
| §5.5, Table 6 (size under H0 with unequal error shapes) | `23_revision_w_shape_size.R` | `tables/revision_w_shape_size.csv` (summary), `revision_w_shape_size_reps.csv` (per replication) |
| §6, data (N = 2136 subsample) | `01_download_brfss2010.R` | local `.rds` (not redistributed) |
| §6, Table 7 and the per-split counts | `03_phq8_70_splits.R`, `22_phq8_criterion.R` | `tables/phq8_70splits_results.csv`, `tables/phq8_criterion.csv` |
| §6, Figure 1 | `22b_phq8_criterion_figure_rev1.R` | `figures/fig_phq8_criterion_rev1.{pdf,png}` |
| §7.1, PMM3-style fourth-order probe | `18_pmm3_symmetric_probe.R`; intervals by `21_revision_experiments.R` | `tables/pmm3_symmetric_probe_*.csv`, `tables/revision_pmm3_ci.csv` |
| Headline numbers of §§5–7 (cross-check) | `verify_reported_values.R` | reads the tables above |

## Notation map (code names to manuscript symbols)

The revised manuscript renamed its symbols; the code keeps its original
variable and column names so that every committed table stays
reproducible. The correspondence is:

| Code | Manuscript | Meaning |
|---|---|---|
| `phi1` | φ₁ = X₁X₂² − X₁²X₂ | target of the naive Δc₃ statistic |
| `phi4` | φ₂ = X₁³ − X₂³ | third-order marginal auxiliary (the single auxiliary of Proposition 2) |
| `phi7` | φ₃ = X₁³X₂ − X₁X₂³ | fourth-order cross auxiliary |
| `phi8` | φ₄ = X₁⁴ − X₂⁴ | fourth-order marginal auxiliary |
| `a1`, `a2` in `05_dgp.R` | α₁ = 1, α₂ = R | confounder effects on X₁, X₂ |
| `lambda` in `05_dgp.R` | λ | true-score loading |
| `lambda` in `08_pmm2_estimator.R` | — | ridge added to the auxiliary covariance matrix (0.01) |
| `K_aug` | q | number of auxiliaries (the basis has q + 1 functions) |
| `K_star` | K* | control-variate weight |
| `naive_phi1`, `pmm2_phi4` in `phq8_criterion.csv` | φ̄₁, φ̄₁ + K*φ̄₂ | naive and single-auxiliary statistics per split |

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
Rscript scripts/23_revision_w_shape_size.R   # ~17 min on 8 cores
Rscript scripts/24_revision_kstar_plim.R     # ~10 s
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

MIT License, Copyright (c) 2026 Serhii Zabolotnii - see `LICENSE`. This
covers the code in `scripts/` and the generated artifacts in `output/`.
It does not extend to the underlying BRFSS 2010 data, which is public
CDC material and is not redistributed here (see `docs/DATA_POLICY.md`),
or to any third-party article cited by the manuscript.

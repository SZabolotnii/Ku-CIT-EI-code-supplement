# RUNBOOK

Last updated: 2026-05-31

Scope: public verification supplement for the Ku-CIT-EI manuscript.

## Headline Verification

Run:

```sh
Rscript scripts/verify_reported_values.R
```

This checks the headline values reported in the manuscript against the generated CSV artifacts in `output/tables/`.

## PMM2 Comparative Monte Carlo

Full run:

```sh
Rscript scripts/09_sim_comparative.R
```

Generated artifacts:

- `output/tables/are_pmm2_vs_naive.csv`
- `output/tables/sim_comparative_summary.csv`
- `data/mc_comparative.rds`

Figures:

```sh
Rscript scripts/10_fig_comparative.R
```

## GSA-LLR Diagnostics

Compact Monte Carlo:

```sh
Rscript scripts/14_gsa_llr_mc_compact.R
```

Smoke run:

```sh
Rscript scripts/14_gsa_llr_mc_compact.R --quick
```

Repair sweep:

```sh
Rscript scripts/15_gsa_llr_repair_sweep.R
```

Expected current result: GSA-LLR controls null behavior and detects skewness-sensitive alternatives, but the strict PHQ-like pure-kurtosis cell remains below the manuscript power gate.

## PMM3 Symmetric Probe

```sh
Rscript scripts/18_pmm3_symmetric_probe.R
```

Expected current result: variance reduction is present, but the practical power gain and Gaussian nuisance guard do not pass the manuscript gate.

## DSGE / PATP-DSGE Hybrid Branch

Pure DSGE probe:

```sh
Rscript scripts/16_dsge_cit_probe.R
```

Hybrid DSGE classifier:

```sh
Rscript scripts/17_dsge_hybrid_cit_classifier.R
```

PATP-DSGE sweep:

```sh
Rscript scripts/19_patp_dsge_hybrid_sweep.R
```

Figure:

```sh
Rscript scripts/20_make_patp_dsge_hybrid_figure.R
```

Expected current result: the positive result is the existing-feature plus PATP-DSGE hybrid at `alpha = 0.75`, not pure DSGE and not a distribution-free test.

## Exclusions

This public supplement excludes:

- raw BRFSS files;
- processed `.rds` files from BRFSS;
- third-party source articles or PDFs;
- private scratch files;
- the manuscript source itself.


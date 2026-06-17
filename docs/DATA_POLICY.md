# Data Policy

This repository is a public code and generated-artifact supplement for a manuscript whose core results are simulation-based and whose Section 6 adds a real-data illustration on public BRFSS 2010 PHQ-8 data.

## Included

- R scripts used for simulation diagnostics and classifier probes.
- Generated CSV tables used to verify reported Monte Carlo results and
  the real-data Section 6 criterion results (`phq8_criterion.csv`).
- Generated manuscript-support figures (including `fig_phq8_criterion`).
- R session snapshots for the main follow-up runs.

## Excluded

- Raw BRFSS 2010 files.
- Processed BRFSS `.rds` files.
- Third-party source articles or PDFs.
- Private scratch files or local checkout metadata.

## BRFSS Boundary

Section 6 of the manuscript uses the public BRFSS 2010 PHQ-8 data as a real-data illustration of the transferability criterion. The supplement ships the generated per-split criterion table (`output/tables/phq8_criterion.csv`) and figure, plus the scripts that reproduce them. The processed `.rds` subset is regenerated locally from the public CDC source by `scripts/01_download_brfss2010.R`; raw BRFSS files and the processed `.rds` are not redistributed here. No individual-level respondent records are included — only aggregate per-split statistics.


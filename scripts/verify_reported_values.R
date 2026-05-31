#!/usr/bin/env Rscript

tol <- 5e-4

read_table <- function(name) {
  path <- file.path("output", "tables", name)
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE)
}

near <- function(actual, expected, label, tolerance = tol) {
  if (!isTRUE(abs(actual - expected) <= tolerance)) {
    stop(sprintf("%s: expected %.6f, got %.6f", label, expected, actual), call. = FALSE)
  }
  cat(sprintf("PASS %-48s %.6f\n", label, actual))
}

is_true <- function(value, label) {
  if (!isTRUE(value)) stop(sprintf("%s: expected TRUE", label), call. = FALSE)
  cat(sprintf("PASS %-48s TRUE\n", label))
}

is_false <- function(value, label) {
  if (!identical(value, FALSE)) stop(sprintf("%s: expected FALSE", label), call. = FALSE)
  cat(sprintf("PASS %-48s FALSE\n", label))
}

sim <- read_table("sim_comparative_summary.csv")
row <- sim[sim$R == 2 & sim$gamma_T == 2.25, ]
near(row$power_gain_pp, -51.9, "PMM2 power loss, R=2 gamma_T=2.25", 0.05)

are <- read_table("are_pmm2_vs_naive.csv")
row <- are[are$R == 2 & are$gamma_U == 2.25 & are$gamma_T == 0 & are$N == 2000, ]
true_dc3 <- row$R * (row$R - 1) * row$gamma_U
attenuation <- 1 - row$mean_pmm2 / true_dc3
near(true_dc3, 4.5, "true Delta c3, R=2 gamma_U=2.25")
near(attenuation, 0.812, "PMM2 attenuation, R=2 gamma_U=2.25 N=2000", 0.001)

pmm3 <- read_table("pmm3_symmetric_probe_results.csv")
var_ratio <- pmm3$var_H0_ng[pmm3$method == "naive_delta_c4"] /
  pmm3$var_H0_ng[pmm3$method == "pmm3_delta_c4"]
power_gain <- pmm3$power_H1_ng[pmm3$method == "pmm3_delta_c4"] -
  pmm3$power_H1_ng[pmm3$method == "naive_delta_c4"]
nuisance <- pmm3$nuisance_H1_gauss[pmm3$method == "pmm3_delta_c4"]
near(var_ratio, 1.127, "PMM3 variance ratio", 0.001)
near(power_gain, 0.02125, "PMM3 power gain", 0.0001)
near(nuisance, 0.295, "PMM3 Gaussian nuisance rejection")

gsa <- read_table("gsa_llr_phase11q_verdict.csv")
is_true(gsa$pass[gsa$criterion == "typeI_acceptable_in_90pct_null_cells"],
        "GSA Type-I verdict")
is_false(gsa$pass[gsa$criterion == "order4_power_gt_0.40_in_phq_like_pure_kurtosis_cell"],
         "GSA PHQ-like pure-kurtosis gate")

hybrid <- read_table("dsge_hybrid_cit_results.csv")
row <- hybrid[hybrid$method == "hybrid", ]
near(row$typeI_H0_ng, 0.030, "Existing + DSGE Type-I")
near(row$nuisance_H1_gauss, 0.0783333333, "Existing + DSGE nuisance", 1e-6)
near(row$power_H1_ng, 0.9016666667, "Existing + DSGE power", 1e-6)
near(row$accuracy, 0.9311111111, "Existing + DSGE accuracy", 1e-6)

hybrid_boot <- read_table("dsge_hybrid_cit_bootstrap.csv")
row <- hybrid_boot[hybrid_boot$candidate == "hybrid", ]
near(row$delta, 0.095, "Existing + DSGE bootstrap delta")
near(row$ci_low, 0.0788888889, "Existing + DSGE bootstrap CI low", 1e-6)
near(row$ci_high, 0.111125, "Existing + DSGE bootstrap CI high", 1e-6)

patp <- read_table("patp_dsge_hybrid_sweep_best.csv")
near(patp$alpha, 0.75, "PATP hybrid alpha")
near(patp$typeI_H0_ng, 0.0233333333, "PATP hybrid Type-I", 1e-6)
near(patp$nuisance_H1_gauss, 0.055, "PATP hybrid nuisance")
near(patp$power_H1_ng, 0.8866666667, "PATP hybrid power", 1e-6)
near(patp$accuracy, 0.9361111111, "PATP hybrid accuracy", 1e-6)

patp_boot <- read_table("patp_dsge_hybrid_sweep_bootstrap.csv")
row <- patp_boot[patp_boot$candidate == "patp_hybrid", ]
near(row$delta, 0.0994444444, "PATP hybrid bootstrap delta", 1e-6)
near(row$ci_low, 0.0838888889, "PATP hybrid bootstrap CI low", 1e-6)
near(row$ci_high, 0.1155555556, "PATP hybrid bootstrap CI high", 1e-6)

patp_verdict <- read_table("patp_dsge_hybrid_sweep_verdict.csv")
is_true(patp_verdict$pass[patp_verdict$criterion == "patp_hybrid_passes_positive_candidate_gates"],
        "PATP hybrid positive-candidate verdict")

cat("\nAll headline verification checks passed.\n")


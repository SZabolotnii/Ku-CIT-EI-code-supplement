#!/usr/bin/env Rscript

tol <- 5e-4

read_table <- function(name) {
  path <- file.path("output", "tables", name)
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE)
}

near <- function(actual, expected, label, tolerance = tol) {
  if (length(actual) != 1L || !is.finite(actual) ||
      !isTRUE(abs(actual - expected) <= tolerance)) {
    stop(sprintf("%s: expected %.6f, got %.6f", label, expected, actual),
         call. = FALSE)
  }
  cat(sprintf("PASS %-58s %.6f\n", label, actual))
}

row1 <- function(dat, expr, label) {
  rows <- dat[expr, , drop = FALSE]
  if (nrow(rows) != 1L) {
    stop(sprintf("%s: expected one matching row, got %d", label, nrow(rows)),
         call. = FALSE)
  }
  rows[1, , drop = FALSE]
}

sim <- read_table("sim_comparative_summary.csv")
expected_are <- data.frame(
  R = c(1, 1, 1.25, 1.25, 2, 2),
  gamma_T = c(0, 2.25, 0, 2.25, 0, 2.25),
  ARE_mean = c(1.89, 4.63, 1.23, 1.34, 5.16, 4.80),
  ARE_median = c(1.89, 4.67, 1.21, 1.32, 4.64, 4.63),
  power_naive_mean = c(0.061, 0.065, 0.157, 0.159, 0.864, 0.853),
  power_pmm2_mean = c(0.062, 0.068, 0.085, 0.066, 0.366, 0.334),
  power_gain_pp = c(0.1, 0.2, -7.2, -9.3, -49.8, -51.9)
)

for (i in seq_len(nrow(expected_are))) {
  er <- expected_are[i, ]
  sr <- row1(sim, sim$R == er$R & sim$gamma_T == er$gamma_T,
             sprintf("summary row R=%s gamma_T=%s", er$R, er$gamma_T))
  near(sr$ARE_mean, er$ARE_mean, sprintf("Table ARE mean R=%s gamma_T=%s", er$R, er$gamma_T), 0.005)
  near(sr$ARE_median, er$ARE_median, sprintf("Table ARE median R=%s gamma_T=%s", er$R, er$gamma_T), 0.005)
  near(sr$power_naive_mean, er$power_naive_mean, sprintf("Table power naive R=%s gamma_T=%s", er$R, er$gamma_T), 0.0005)
  near(sr$power_pmm2_mean, er$power_pmm2_mean, sprintf("Table power PMM2 R=%s gamma_T=%s", er$R, er$gamma_T), 0.0005)
  near(sr$power_gain_pp, er$power_gain_pp, sprintf("Table power loss R=%s gamma_T=%s", er$R, er$gamma_T), 0.05)
}

are <- read_table("are_pmm2_vs_naive.csv")
bias_cells <- data.frame(
  R = c(1.25, 1.25, 1.25, 1.25, 2, 2, 2, 2),
  gamma_U = c(0.75, 0.75, 2.25, 2.25, 0.75, 0.75, 2.25, 2.25),
  N = c(500, 2000, 500, 2000, 500, 2000, 500, 2000),
  true_dc3 = c(0.234, 0.234, 0.703, 0.703, 1.500, 1.500, 4.500, 4.500),
  mean_naive = c(0.297, 0.239, 0.649, 0.723, 1.475, 1.514, 4.431, 4.452),
  mean_pmm2 = c(0.172, 0.134, 0.219, 0.309, 0.187, 0.193, 0.694, 0.846),
  attenuation = c(0.265, 0.427, 0.688, 0.560, 0.875, 0.871, 0.846, 0.812),
  power_gain_pp = c(-2.0, -1.0, -5.5, -9.0, -29.5, -74.0, -74.0, -34.5)
)

for (i in seq_len(nrow(bias_cells))) {
  bc <- bias_cells[i, ]
  ar <- row1(are, are$R == bc$R & are$gamma_U == bc$gamma_U &
               are$gamma_T == 0 & are$N == bc$N,
             sprintf("bias row R=%s gamma_U=%s N=%s", bc$R, bc$gamma_U, bc$N))
  true_dc3 <- ar$R * (ar$R - 1) * ar$gamma_U
  attenuation <- 1 - ar$mean_pmm2 / true_dc3
  near(true_dc3, bc$true_dc3, sprintf("Table bias true Dc3 R=%s gamma_U=%s N=%s", bc$R, bc$gamma_U, bc$N), 0.0006)
  near(ar$mean_naive, bc$mean_naive, sprintf("Table bias naive mean R=%s gamma_U=%s N=%s", bc$R, bc$gamma_U, bc$N), 0.0006)
  near(ar$mean_pmm2, bc$mean_pmm2, sprintf("Table bias PMM2 mean R=%s gamma_U=%s N=%s", bc$R, bc$gamma_U, bc$N), 0.0006)
  near(attenuation, bc$attenuation, sprintf("Table bias attenuation R=%s gamma_U=%s N=%s", bc$R, bc$gamma_U, bc$N), 0.004)
  near(ar$power_gain_pp, bc$power_gain_pp, sprintf("Table bias power loss R=%s gamma_U=%s N=%s", bc$R, bc$gamma_U, bc$N), 0.05)
}

boot <- read_table("revision_bootstrap_sensitivity.csv")
strong_pct_naive <- row1(boot, boot$scenario == "H1_strong" &
                           boot$estimator == "naive_delta_c3" &
                           boot$B == 1000 & boot$ci_method == "percentile",
                         "revision bootstrap strong naive percentile")
strong_pct_pmm2 <- row1(boot, boot$scenario == "H1_strong" &
                          boot$estimator == "pmm2_delta_c3" &
                          boot$B == 1000 & boot$ci_method == "percentile",
                        "revision bootstrap strong PMM2 percentile")
strong_bca_naive <- row1(boot, boot$scenario == "H1_strong" &
                           boot$estimator == "naive_delta_c3" &
                           boot$B == 1000 & boot$ci_method == "bca",
                         "revision bootstrap strong naive BCa")
strong_bca_pmm2 <- row1(boot, boot$scenario == "H1_strong" &
                          boot$estimator == "pmm2_delta_c3" &
                          boot$B == 1000 & boot$ci_method == "bca",
                        "revision bootstrap strong PMM2 BCa")
near(strong_pct_naive$rejection_rate, 0.8833333333, "Revision bootstrap percentile naive power", 1e-6)
near(strong_pct_pmm2$rejection_rate, 0.2500000000, "Revision bootstrap percentile PMM2 power", 1e-6)
near(strong_bca_naive$rejection_rate, 0.9000000000, "Revision bootstrap BCa naive power", 1e-6)
near(strong_bca_pmm2$rejection_rate, 0.3000000000, "Revision bootstrap BCa PMM2 power", 1e-6)

heavy <- read_table("revision_heavytail_sensitivity.csv")
h500_naive <- row1(heavy, heavy$R == 2 & heavy$N == 500 & heavy$estimator == "naive_delta_c3",
                   "heavy-tail R=2 N=500 naive")
h500_pmm2 <- row1(heavy, heavy$R == 2 & heavy$N == 500 & heavy$estimator == "pmm2_delta_c3",
                  "heavy-tail R=2 N=500 PMM2")
h2000_naive <- row1(heavy, heavy$R == 2 & heavy$N == 2000 & heavy$estimator == "naive_delta_c3",
                    "heavy-tail R=2 N=2000 naive")
h2000_pmm2 <- row1(heavy, heavy$R == 2 & heavy$N == 2000 & heavy$estimator == "pmm2_delta_c3",
                   "heavy-tail R=2 N=2000 PMM2")
near(h500_naive$power, 0.6625, "Heavy-tail R=2 N=500 naive power", 1e-6)
near(h500_pmm2$power, 0.0750, "Heavy-tail R=2 N=500 PMM2 power", 1e-6)
near(h2000_naive$power, 0.9875, "Heavy-tail R=2 N=2000 naive power", 1e-6)
near(h2000_pmm2$power, 0.2500, "Heavy-tail R=2 N=2000 PMM2 power", 1e-6)

heavy_ratios <- read_table("revision_heavytail_ratios.csv")
hr <- row1(heavy_ratios, heavy_ratios$scenario == "Tukey_g0.35_h0.10_R2_N2000",
           "heavy-tail R=2 N=2000 signal ratio")
near(hr$signal_ratio, 0.1143563359, "Heavy-tail R=2 N=2000 PMM2/naive signal ratio", 1e-9)

dcov <- read_table("revision_dcov_sanity.csv")
d0 <- row1(dcov, dcov$scenario == "H0_common_true_score", "raw dCov H0 row")
near(d0$raw_dcov_rejection_rate, 1.0, "Raw dCov H0 observed-score rejection", 1e-12)
near(d0$ci_lo, 0.9010990077, "Raw dCov H0 Wilson CI low", 1e-9)

pmm3 <- read_table("pmm3_symmetric_probe_results.csv")
var_ratio <- pmm3$var_H0_ng[pmm3$method == "naive_delta_c4"] /
  pmm3$var_H0_ng[pmm3$method == "pmm3_delta_c4"]
power_gain <- pmm3$power_H1_ng[pmm3$method == "pmm3_delta_c4"] -
  pmm3$power_H1_ng[pmm3$method == "naive_delta_c4"]
nuisance <- pmm3$nuisance_H1_gauss[pmm3$method == "pmm3_delta_c4"]
near(var_ratio, 1.127, "PMM3 diagnostic variance ratio", 0.001)
near(power_gain, 0.02125, "PMM3 diagnostic power gain", 0.0001)
near(nuisance, 0.295, "PMM3 diagnostic nuisance rejection", 0.0001)

pmm3_ci <- read_table("revision_pmm3_ci.csv")
pmm3_type1 <- row1(pmm3_ci, pmm3_ci$method == "pmm3_delta_c4" &
                     pmm3_ci$metric == "typeI_H0_ng", "PMM3 type-I CI")
near(pmm3_type1$rate, 0.07125, "PMM3 diagnostic Type-I rate", 1e-6)

# --- Section 6: real-data PHQ-8 / BRFSS 2010 transferability-criterion checks ---
# Self-contained: reads the committed per-split table (no BRFSS .rds needed).
crit <- read_table("phq8_criterion.csv")
near(nrow(crit), 70, "PHQ-8 criterion: number of half-splits", 1e-9)
near(sum(crit$rel_attenuation < 0), 70, "PHQ-8: splits attenuated toward zero", 1e-9)
near(sum(crit$sig_naive), 58, "PHQ-8: naive detections (alpha=0.05)", 1e-9)
near(sum(crit$sig_pmm2), 52, "PHQ-8: PMM2 detections (alpha=0.05)", 1e-9)
near(sum(crit$sig_naive == 1 & crit$sig_pmm2 == 0), 6,
     "PHQ-8: detections erased by PMM2", 1e-9)
near(sum(crit$sig_naive == 0 & crit$sig_pmm2 == 1), 0,
     "PHQ-8: detections created by PMM2", 1e-9)
near(sum(crit$criterion_safe), 10, "PHQ-8: no-bias precondition met", 1e-9)
med_K <- median(crit$K_cv)
if (!(med_K > 0.13 && med_K < 0.25))
  stop(sprintf("PHQ-8 median K* out of expected range: %.4f", med_K), call. = FALSE)
cat(sprintf("PASS %-58s %.6f\n", "PHQ-8: real-data H0-optimal K* (median ~ sim 0.24)", med_K))
med_att <- median(abs(crit$rel_attenuation))
if (!(med_att > 0.50 && med_att < 0.60))
  stop(sprintf("PHQ-8 median attenuation out of expected range: %.4f", med_att), call. = FALSE)
cat(sprintf("PASS %-58s %.6f\n", "PHQ-8: median |relative attenuation|", med_att))

cat("\nAll current-manuscript verification checks passed.\n")

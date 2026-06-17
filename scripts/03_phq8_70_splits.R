# 03_phq8_70_splits.R
# ---------------------------------------------------------------------------
# Replicates W&S (2026) Table 5 / Figures 8-10 computation:
# Computes Δ̂c₃, Δ̂c₄ + 95% percentile bootstrap CIs for all C(8,4) = 70
# half-splits of the PHQ-8 items on the BRFSS 2010 high-risk subset.
#
# Requires:  source-documents/data/brfss2010_phq8_highrisk.rds
#            scripts/02_naive_estimators.R
# Produces:  output/tables/phq8_70splits_results.csv
#
# Run: Rscript scripts/03_phq8_70_splits.R
# ---------------------------------------------------------------------------

library(here)
source(here("scripts", "02_naive_estimators.R"))

# ── 0. Setup ────────────────────────────────────────────────────────────────
phq_rds  <- here("source-documents", "data", "brfss2010_phq8_highrisk.rds")
out_csv  <- here("output", "tables", "phq8_70splits_results.csv")
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

SEED_BASE <- 20260525L
B_BOOT    <- 2000L
ALPHA     <- 0.05

# ── 1. Load data ─────────────────────────────────────────────────────────────
cat("Loading BRFSS PHQ-8 high-risk subset...\n")
phq_data <- readRDS(phq_rds)

phq_vars <- c("ADPLEASR", "ADDOWN", "ADSLEEP", "ADENERGY",
               "ADEAT1",   "ADFAIL",  "ADTHINK",  "ADMOVE")
phq_vars <- phq_vars[phq_vars %in% names(phq_data)]
cat(sprintf("PHQ-8 items available: %d/8: %s\n", length(phq_vars),
            paste(phq_vars, collapse = ", ")))
if (length(phq_vars) < 8L) stop("Not all 8 PHQ-8 items found. Re-run 01_download_brfss2010.R.")

data_mat <- as.matrix(phq_data[, phq_vars])
N <- nrow(data_mat)
cat(sprintf("N = %d (W&S report N = 2136)\n", N))

# ── 2. All C(8,4) = 70 splits ───────────────────────────────────────────────
splits <- combn(8L, 4L)          # 4 × 70 matrix: each column = first-half item indices
n_splits <- ncol(splits)
cat(sprintf("Processing %d splits (C(8,4) = %d)...\n", n_splits, n_splits))

# ── 3. Compute results ───────────────────────────────────────────────────────
results <- vector("list", n_splits)

for (s in seq_len(n_splits)) {
  idx1 <- splits[, s]
  idx2 <- setdiff(1L:8L, idx1)

  x1 <- rowSums(data_mat[, idx1])
  x2 <- rowSums(data_mat[, idx2])

  # Pearson correlation and Spearman-Brown
  cor_x1x2 <- cor(x1, x2)
  sb_rel    <- spearman_brown(cor_x1x2)

  # Δc₃
  res_c3 <- bootstrap_ci(delta_c3_naive, x1, x2, B = B_BOOT, alpha = ALPHA,
                          seed = SEED_BASE + s)

  # Δc₄
  res_c4 <- bootstrap_ci(delta_c4_naive, x1, x2, B = B_BOOT, alpha = ALPHA,
                          seed = SEED_BASE + 100L + s)

  results[[s]] <- data.frame(
    split_id      = s,
    items_1       = paste(phq_vars[idx1], collapse = "+"),
    items_2       = paste(phq_vars[idx2], collapse = "+"),
    N             = N,
    cor_X1X2      = round(cor_x1x2, 4),
    sb_reliability= round(sb_rel,    4),
    delta_c3      = round(res_c3["estimate"], 5),
    ci_lo_c3      = round(res_c3["ci_lo"],    5),
    ci_hi_c3      = round(res_c3["ci_hi"],    5),
    sig_c3        = res_c3["sig"],
    delta_c4      = round(res_c4["estimate"], 5),
    ci_lo_c4      = round(res_c4["ci_lo"],    5),
    ci_hi_c4      = round(res_c4["ci_hi"],    5),
    sig_c4        = res_c4["sig"],
    stringsAsFactors = FALSE
  )

  if (s %% 10L == 0L || s == n_splits) {
    cat(sprintf("  Split %2d/%d done.\n", s, n_splits))
  }
}

df <- do.call(rbind, results)
rownames(df) <- NULL

# ── 4. Save ───────────────────────────────────────────────────────────────────
write.csv(df, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s  (%d rows × %d cols)\n", out_csv, nrow(df), ncol(df)))

# ── 5. Acceptance-criteria checks ────────────────────────────────────────────
cat("\n=== Acceptance-criteria checks ===\n")

# A1: 70 rows
ac1 <- nrow(df) == 70L
cat(sprintf("A1: 70 rows:                 %s (actual %d)\n",
            if (ac1) "PASS" else "FAIL", nrow(df)))

# A2: ≥10 columns
ac2 <- ncol(df) >= 10L
cat(sprintf("A2: ≥10 columns:             %s (actual %d)\n",
            if (ac2) "PASS" else "FAIL", ncol(df)))

# A3: Pearson cor range [0.77, 0.86]
cor_range <- range(df$cor_X1X2)
ac3 <- cor_range[1] >= 0.77 & cor_range[2] <= 0.86
cat(sprintf("A3: cor ∈ [0.77, 0.86]:      %s (actual [%.3f, %.3f])\n",
            if (ac3) "PASS" else "WARN", cor_range[1], cor_range[2]))

# A4: Spearman-Brown reliability [0.87, 0.93]
sb_range <- range(df$sb_reliability)
ac4 <- sb_range[1] >= 0.87 & sb_range[2] <= 0.93
cat(sprintf("A4: SB ∈ [0.87, 0.93]:       %s (actual [%.3f, %.3f])\n",
            if (ac4) "PASS" else "WARN", sb_range[1], sb_range[2]))

# A5: sig_c3 ≈ 58/70 (82.9%) ± 3
n_sig_c3 <- sum(df$sig_c3)
ac5 <- abs(n_sig_c3 - 58L) <= 3L
cat(sprintf("A5: sig Δc₃ ≈ 58 ± 3:        %s (actual %d, %.1f%%)\n",
            if (ac5) "PASS" else "WARN", n_sig_c3, 100 * n_sig_c3 / 70))

# A6: sig_c4 ≈ 48/70 (68.6%) ± 3
n_sig_c4 <- sum(df$sig_c4)
ac6 <- abs(n_sig_c4 - 48L) <= 3L
cat(sprintf("A6: sig Δc₄ ≈ 48 ± 3:        %s (actual %d, %.1f%%)\n",
            if (ac6) "PASS" else "WARN", n_sig_c4, 100 * n_sig_c4 / 70))

cat("\n=== Summary ===\n")
cat(sprintf("Correlations: mean=%.3f, range=[%.3f, %.3f]\n",
            mean(df$cor_X1X2), cor_range[1], cor_range[2]))
cat(sprintf("SB reliability: mean=%.3f, range=[%.3f, %.3f]\n",
            mean(df$sb_reliability), sb_range[1], sb_range[2]))
cat(sprintf("Significant Δc₃: %d/70 (%.1f%%)  [W&S: 58/70 = 82.9%%]\n",
            n_sig_c3, 100 * n_sig_c3 / 70))
cat(sprintf("Significant Δc₄: %d/70 (%.1f%%)  [W&S: 48/70 = 68.6%%]\n",
            n_sig_c4, 100 * n_sig_c4 / 70))

cat("\n=== DONE ===\n")

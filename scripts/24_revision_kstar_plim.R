# 24_revision_kstar_plim.R
# ---------------------------------------------------------------------------
# Provenance of K* ~ 0.24 in Corollary 3 and Appendix A (SORT revision 1,
# 2026-09-24). The submitted text called it "estimated from the Monte Carlo
# runs of Section 5", but no table in output/ carried it. This script measures
# it on the unchanged Section-5 DGP (scripts/05_dgp.R).
#
# What is measured, per cell, over n_draws independent samples of size N:
#   * the weights the PMM2 estimator computes from the sample (as in
#     scripts/08_pmm2_estimator.R, without the ridge, which is negligible at
#     this N), single auxiliary (phi4 = X1^3 - X2^3, phi_2 in the manuscript)
#     and the three-auxiliary basis of Table 1;
#   * the sample means of phi1 and of the auxiliaries;
#   * the resulting statistic phi1_bar + K' aux_bar.
# At N = 1e6 the draw-to-draw spread is small, so the means are the
# probability limits to the reported precision.
#
# Cells: R = 2, gamma_U = 2.25, cor(E1, E2) = 0.38 (the Corollary-3 cell) and
# the Section-5 H0 cell (R = 1, no confounder, sigma_W = 1).
# Output: output/tables/revision_kstar_plim.csv
#         output/session_info/revision_kstar_plim_session_info.txt
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(here))
options(pmm2.skip_tests = TRUE)
source(here("scripts", "05_dgp.R"))

N       <- 1e6
N_DRAWS <- 20L
SEED    <- 20260924L

sample_weights <- function(x1, x2, basis) {
  x1 <- x1 - mean(x1); x2 <- x2 - mean(x2)
  phi1 <- x1 * x2^2 - x1^2 * x2
  aug  <- list(phi4 = x1^3 - x2^3,
               phi7 = x1^3 * x2 - x1 * x2^3,
               phi8 = x1^4 - x2^4)[basis]
  A  <- do.call(cbind, aug)
  Ac <- sweep(A, 2, colMeans(A))
  n1 <- length(x1) - 1
  K  <- -solve(crossprod(Ac) / n1,
               as.vector(crossprod(Ac, phi1 - mean(phi1)) / n1))
  names(K) <- paste0("K_", basis)
  m <- colMeans(A); names(m) <- paste0("mean_", basis)
  c(K, mean_phi1 = mean(phi1), m, statistic = mean(phi1) + sum(K * colMeans(A)))
}

cells <- list(
  list(cell = "H1_R2_gU2.25", R = 2, cor = 0.38, gamma_U = 2.25),
  list(cell = "H0_section5",  R = 1, cor = 0.00, gamma_U = 2.25)
)
bases <- list(single = "phi4", three = c("phi4", "phi7", "phi8"))

set.seed(SEED)
rows <- list()
for (cc in cells) {
  draws <- list(single = NULL, three = NULL)
  for (b in seq_len(N_DRAWS)) {
    d <- generate_ws_dgp(N = N, gamma_U = cc$gamma_U, gamma_T = 0, gamma_W = 0,
                         R = cc$R, cor_E1E2 = cc$cor)
    for (bn in names(bases))
      draws[[bn]] <- rbind(draws[[bn]], sample_weights(d$x1, d$x2, bases[[bn]]))
  }
  for (bn in names(bases)) {
    X <- draws[[bn]]
    rows[[length(rows) + 1]] <- data.frame(
      cell = cc$cell, R = cc$R, cor_E1E2 = cc$cor, gamma_U = cc$gamma_U,
      basis = bn, quantity = colnames(X),
      mean = colMeans(X), sd_over_draws = apply(X, 2, sd),
      n_draws = N_DRAWS, N = N, row.names = NULL)
  }
}
out <- do.call(rbind, rows)

# Closed-form check of Corollary 3 with the measured single-auxiliary weight.
k1  <- out$mean[out$cell == "H1_R2_gU2.25" & out$basis == "single" &
                out$quantity == "K_phi4"]
st1 <- out$mean[out$cell == "H1_R2_gU2.25" & out$basis == "single" &
                out$quantity == "statistic"]
pred <- 2.25 * (2 - 1) * (2 - k1 * 7)
cat(sprintf("R = 2 cell: single-auxiliary weight %.4f; predicted centre %.3f, measured %.3f; attenuation %.1f %% (target 4.50)\n",
            k1, pred, st1, 100 * (1 - st1 / 4.5)))

dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, here("output", "tables", "revision_kstar_plim.csv"), row.names = FALSE)
dir.create(here("output", "session_info"), showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()),
           here("output", "session_info", "revision_kstar_plim_session_info.txt"))
cat("written: output/tables/revision_kstar_plim.csv\n")

# 23_revision_w_shape_size.R
# ---------------------------------------------------------------------------
# SORT revision 1 (2026-09-17), concern "do the zero H0 means of the
# auxiliaries need W1 and W2 to have the same distributional shape?"
#
# Null-side SIZE experiment. The errors are INDEPENDENT in every cell
# (no confounder U at all, R = 1, cor(E1,E2) = 0), so Delta c3 = 0 and a
# level-alpha test should reject at rate alpha. What varies is the SHAPE of
# the idiosyncratic errors: W1 and W2 get separate skewness arguments.
#
#   E[phi4] = E[X1^3 - X2^3] = kappa3(W1) - kappa3(W2)            (R = 1)
#   E[phi7] = 0                        (equal error variances)
#   E[phi8] = E[X1^4 - X2^4] = kappa4(W1) - kappa4(W2)  (equal error variances)
#
# so the PMM2 statistic mean(phi1) + K' mean(Phi_aug) carries the fixed
# population offset K' mu_aux although Delta c3 = 0; the naive statistic
# mean(phi1) does not (E[phi1] does not involve W at all).
#
# NOTHING is re-implemented here: the DGP component generator, the two
# estimators and the percentile bootstrap are the repository's own
# functions, sourced from
#   scripts/05_dgp.R            generate_component(), generate_ws_dgp()
#   scripts/02_naive_estimators.R  delta_c3_naive(), bootstrap_ci()
#   scripts/08_pmm2_estimator.R    delta_c3_pmm2()  (full basis and basis="phi4")
# and wilson_ci() is evaluated verbatim out of scripts/21_revision_experiments.R
# (that file cannot be source()d without re-running all of its experiments).
#
# The only new code is (i) a wrapper DGP that is line-for-line the
# independent-errors branch of generate_ws_dgp() with gamma_W split into
# gamma_W1 / gamma_W2, (ii) a DIAGNOSTIC that recomputes the weights K* the
# way 08 does, only to report them -- every rep asserts that these weights
# reproduce delta_c3_pmm2()'s own output to 1e-9 -- and (iii) bookkeeping.
#
# How this matches Section 5 (scripts/09_sim_comparative.R, power pass):
#   alpha = 0.05; percentile bootstrap via bootstrap_ci(); B = 200;
#   200 replications; seed_r = SEED_BASE + 500000 + cond_idx*10000 + r;
#   data seed = seed_r, naive bootstrap seed = seed_r + 1, PMM2 bootstrap
#   seed = seed_r + 2 (single-auxiliary PMM2: seed_r + 3, new);
#   H0 cells: no U, lambda = 1, sigma_T = 1, sigma_W = sigma_T = 1.
#   NOTE: the manuscript text says sigma_T = sigma_W = 1/sqrt(2); the code
#   that produced Section 5 uses 1 (09 never passes sigma_T, 05 defaults to
#   1 and sets sigma_W = sigma_T in the independent-errors branch). The
#   primary runs follow the CODE. A third run repeats the 200-rep run with
#   the SAME seeds at 1/sqrt(2); the data are then exactly the primary data
#   times 1/sqrt(2), and the only thing that is not scale-equivariant is the
#   Tikhonov term lambda = 0.01 inside delta_c3_pmm2().
#   SEED_BASE follows the repo convention 20260525 + <script number>.
#
# Runs (never mixed; column `run`):
#   sec5_match_200        200 reps, B = 200, sigma = 1         offset  500000
#   extended_1000        1000 reps, B = 200, sigma = 1         offset 1000000
#                        (independent of the 200-rep run: disjoint seeds)
#   scale_invsqrt2_200    200 reps, B = 200, sigma = 1/sqrt(2) offset  500000
#                        (same seeds as sec5_match_200, paired)
#
# Run from the repository root:
#   Rscript scripts/23_revision_w_shape_size.R
# Quick smoke test (20 reps, writes *_QUICK files only):
#   Rscript scripts/23_revision_w_shape_size.R --quick
# ---------------------------------------------------------------------------

options(pmm2.skip_tests = TRUE)
SKIP_PMM2_TESTS <- TRUE

suppressPackageStartupMessages({
  library(here)
  library(parallel)
})

source(here("scripts", "05_dgp.R"))
source(here("scripts", "02_naive_estimators.R"))
source(here("scripts", "08_pmm2_estimator.R"))

# Evaluate one top-level function definition out of a script that cannot be
# sourced as a whole. The definition is used verbatim.
extract_fn <- function(path, name) {
  exprs <- parse(path, keep.source = FALSE)
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name(name))) {
      eval(e, envir = globalenv())
      return(invisible(TRUE))
    }
  }
  stop("definition of ", name, " not found in ", path)
}
extract_fn(here("scripts", "21_revision_experiments.R"), "wilson_ci")
# Single-auxiliary point estimator of scripts/22 (no Tikhonov term); used
# only for a parity check against delta_c3_pmm2(basis = "phi4").
extract_fn(here("scripts", "22_phq8_criterion.R"), "pmm2_phi4_point")

args  <- commandArgs(trailingOnly = TRUE)
QUICK <- "--quick" %in% args

# ── Configuration ───────────────────────────────────────────────────────────
SEED_BASE <- 20260525L + 23L
ALPHA     <- 0.05
B_BOOT    <- 200L
N_CORES   <- max(1L, min(8L, parallel::detectCores() - 2L))

N_vals <- c(500L, 2000L, 5000L)
shape_conditions <- data.frame(
  condition = c("control_gaussian", "control_identical_skew",
                "treatment_skew_vs_gaussian", "treatment_mild_skew_vs_gaussian"),
  gamma_W1  = c(0, 2.25, 2.25, 0.75),
  gamma_W2  = c(0, 2.25, 0,    0),
  stringsAsFactors = FALSE
)
GAMMA_T <- 0

# cond_idx: N outer, shape inner (1..12)
conditions <- do.call(rbind, lapply(seq_along(N_vals), function(i) {
  cbind(N = N_vals[i], shape_conditions)
}))
conditions$cond_idx <- seq_len(nrow(conditions))

runs <- list(
  list(run = "sec5_match_200",     n_reps = 200L,  sigma = 1,           offset = 500000L),
  list(run = "extended_1000",      n_reps = 1000L, sigma = 1,           offset = 1000000L),
  list(run = "scale_invsqrt2_200", n_reps = 200L,  sigma = 1 / sqrt(2), offset = 500000L)
)
if (QUICK) runs <- lapply(runs, function(r) { r$n_reps <- 20L; r })

# ── Wrapper DGP ─────────────────────────────────────────────────────────────
# Line-for-line the `else` (independent-errors) branch of generate_ws_dgp():
# same RNG call order (seed, T, E1, E2), lambda * T + E_i, sigma_W = sigma_T.
# The only change: W1 and W2 take separate skewness arguments.
generate_ws_h0_wshape <- function(N, gamma_T = 0, gamma_W1 = 0, gamma_W2 = 0,
                                  sigma_T = 1.0, lambda = 1.0, seed = NULL) {
  stopifnot(N >= 2L)
  if (!is.null(seed)) set.seed(seed)
  T_raw   <- generate_component(N, gamma_T, 0, seed = NULL) * sigma_T
  sigma_W <- sigma_T
  E1 <- generate_component(N, gamma_W1, 0, seed = NULL) * sigma_W
  E2 <- generate_component(N, gamma_W2, 0, seed = NULL) * sigma_W
  data.frame(x1 = lambda * T_raw + E1, x2 = lambda * T_raw + E2)
}

# Single-auxiliary PMM2: the estimator's own single-basis path.
delta_c3_pmm2_phi4 <- function(x1, x2) delta_c3_pmm2(x1, x2, basis = "phi4")

# ── Diagnostic: the weights K* as 08_pmm2_estimator.R computes them ─────────
# Reported only. `recon` is what these weights give; the caller asserts that
# it equals delta_c3_pmm2()'s own return value.
pmm2_weights_diag <- function(x1, x2, basis = c("phi4", "phi7", "phi8"),
                              lambda = 0.01) {
  N  <- length(x1)
  x1 <- x1 - mean(x1); x2 <- x2 - mean(x2)
  phi1 <- x1 * x2^2 - x1^2 * x2
  aug  <- list(phi4 = x1^3 - x2^3,
               phi7 = x1^3 * x2 - x1 * x2^3,
               phi8 = x1^4 - x2^4)[basis]
  Phi   <- do.call(cbind, aug)
  Phi_c <- sweep(Phi, 2, colMeans(Phi), "-")
  F_mat <- crossprod(Phi_c) / (N - 1)
  B_vec <- as.vector(crossprod(Phi_c, phi1 - mean(phi1)) / (N - 1))
  K     <- -solve(F_mat + lambda * diag(ncol(Phi)), B_vec)
  names(K) <- basis
  list(K = K, aux_means = colMeans(Phi),
       recon = mean(phi1) + sum(K * colMeans(Phi)))
}

# ── Closed forms for a standardized shifted-gamma component scaled by sigma ─
# skewness gamma, excess kurtosis 6/shape = 1.5 * gamma^2 (shape = 4/gamma^2)
kappa3_W <- function(gamma, sigma) gamma * sigma^3
kappa4_W <- function(gamma, sigma) 1.5 * gamma^2 * sigma^4

# ── One replication ─────────────────────────────────────────────────────────
one_rep <- function(r, cond, sigma, offset) {
  seed_r <- SEED_BASE + offset + cond$cond_idx * 10000L + r
  d <- generate_ws_h0_wshape(N = cond$N, gamma_T = GAMMA_T,
                             gamma_W1 = cond$gamma_W1, gamma_W2 = cond$gamma_W2,
                             sigma_T = sigma, lambda = 1.0, seed = seed_r)
  ci_n <- bootstrap_ci(delta_c3_naive,     d$x1, d$x2, B = B_BOOT, alpha = ALPHA, seed = seed_r + 1L)
  ci_p <- bootstrap_ci(delta_c3_pmm2,      d$x1, d$x2, B = B_BOOT, alpha = ALPHA, seed = seed_r + 2L)
  ci_s <- bootstrap_ci(delta_c3_pmm2_phi4, d$x1, d$x2, B = B_BOOT, alpha = ALPHA, seed = seed_r + 3L)
  wf <- pmm2_weights_diag(d$x1, d$x2)
  ws <- pmm2_weights_diag(d$x1, d$x2, basis = "phi4")
  c(rep = r, seed = seed_r,
    est_naive = ci_n[["estimate"]], est_pmm2_full = ci_p[["estimate"]],
    est_pmm2_phi4 = ci_s[["estimate"]],
    sig_naive = ci_n[["sig"]], sig_pmm2_full = ci_p[["sig"]],
    sig_pmm2_phi4 = ci_s[["sig"]],
    Kf_phi4 = wf$K[["phi4"]], Kf_phi7 = wf$K[["phi7"]], Kf_phi8 = wf$K[["phi8"]],
    Ks_phi4 = ws$K[["phi4"]],
    aux_phi4 = wf$aux_means[["phi4"]], aux_phi7 = wf$aux_means[["phi7"]],
    aux_phi8 = wf$aux_means[["phi8"]],
    diag_err_full = abs(wf$recon - ci_p[["estimate"]]),
    diag_err_phi4 = abs(ws$recon - ci_s[["estimate"]]),
    parity22_phi4 = abs(pmm2_phi4_point(d$x1, d$x2) - ci_s[["estimate"]]))
}

run_cell <- function(cond, run) {
  t0  <- Sys.time()
  out <- parallel::mclapply(seq_len(run$n_reps), one_rep, cond = cond,
                            sigma = run$sigma, offset = run$offset,
                            mc.cores = N_CORES, mc.preschedule = TRUE)
  bad <- vapply(out, function(o) inherits(o, "try-error") || is.null(o), logical(1))
  if (any(bad)) stop("worker failure in cell ", cond$cond_idx, " run ", run$run)
  reps <- as.data.frame(do.call(rbind, out))
  reps <- reps[order(reps$rep), ]
  stopifnot(nrow(reps) == run$n_reps, all(is.finite(as.matrix(reps))))
  attr(reps, "elapsed") <- as.numeric(Sys.time() - t0, units = "secs")
  reps
}

# Population weights from one very large sample (reference only).
population_weights <- function(cond, sigma) {
  d <- generate_ws_h0_wshape(N = 2000000L, gamma_T = GAMMA_T,
                             gamma_W1 = cond$gamma_W1, gamma_W2 = cond$gamma_W2,
                             sigma_T = sigma,
                             seed = SEED_BASE + 9000000L +
                               match(cond$condition, shape_conditions$condition))
  list(full = pmm2_weights_diag(d$x1, d$x2)$K,
       phi4 = pmm2_weights_diag(d$x1, d$x2, basis = "phi4")$K)
}

summarise_cell <- function(reps, cond, run, kpop) {
  s  <- run$sigma
  dk3 <- kappa3_W(cond$gamma_W1, s) - kappa3_W(cond$gamma_W2, s)
  dk4 <- kappa4_W(cond$gamma_W1, s) - kappa4_W(cond$gamma_W2, s)
  n  <- nrow(reps)
  mk <- function(est_col, sig_col, estimator, K4, K7, K8, Kp) {
    w  <- wilson_ci(sum(reps[[sig_col]]), n)
    m4 <- if (is.null(K4)) NA_real_ else mean(K4)
    m7 <- if (is.null(K7)) NA_real_ else mean(K7)
    m8 <- if (is.null(K8)) NA_real_ else mean(K8)
    off_k3 <- if (is.null(K4)) 0 else m4 * dk3
    off_all <- if (is.null(K4)) 0 else off_k3 + (if (is.null(K8)) 0 else m8 * dk4)
    off_pop <- if (is.null(Kp)) 0 else
      Kp[["phi4"]] * dk3 + (if ("phi8" %in% names(Kp)) Kp[["phi8"]] * dk4 else 0)
    data.frame(
      run = run$run, n_reps = n, B = B_BOOT, alpha = ALPHA,
      sigma_T = s, sigma_W = s, lambda = 1, gamma_T = GAMMA_T,
      cond_idx = cond$cond_idx, condition = cond$condition,
      gamma_W1 = cond$gamma_W1, gamma_W2 = cond$gamma_W2, N = cond$N,
      estimator = estimator,
      n_reject = sum(reps[[sig_col]]),
      rejection_rate = w[["rate"]], ci_lo = w[["ci_lo"]], ci_hi = w[["ci_hi"]],
      mean_estimate = mean(reps[[est_col]]), sd_estimate = stats::sd(reps[[est_col]]),
      mcse_mean_estimate = stats::sd(reps[[est_col]]) / sqrt(n),
      mean_K_phi4 = m4, mean_K_phi7 = m7, mean_K_phi8 = m8,
      sd_K_phi4 = if (is.null(K4)) NA_real_ else stats::sd(K4),
      kappa3_W1_minus_W2 = dk3, kappa4_W1_minus_W2 = dk4,
      offset_theory_k3_only = off_k3,
      offset_theory_all_aux = off_all,
      Kpop_phi4 = if (is.null(Kp)) NA_real_ else Kp[["phi4"]],
      Kpop_phi7 = if (is.null(Kp) || !"phi7" %in% names(Kp)) NA_real_ else Kp[["phi7"]],
      Kpop_phi8 = if (is.null(Kp) || !"phi8" %in% names(Kp)) NA_real_ else Kp[["phi8"]],
      offset_theory_Kpop = off_pop,
      seed_base = SEED_BASE, seed_first = min(reps$seed), seed_last = max(reps$seed),
      cell_elapsed_sec = attr(reps, "elapsed"),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    mk("est_naive",     "sig_naive",     "naive_delta_c3",      NULL, NULL, NULL, NULL),
    mk("est_pmm2_full", "sig_pmm2_full", "pmm2_full_3aux",
       reps$Kf_phi4, reps$Kf_phi7, reps$Kf_phi8, kpop$full),
    mk("est_pmm2_phi4", "sig_pmm2_phi4", "pmm2_single_aux_phi4",
       reps$Ks_phi4, NULL, NULL, kpop$phi4)
  )
}

# ── Pre-flight parity checks ────────────────────────────────────────────────
cat("=== 23_revision_w_shape_size.R ===\n")
cat(sprintf("Started: %s | cores: %d | quick: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), N_CORES, QUICK))
for (N in N_vals) {
  a <- generate_ws_dgp(N = N, gamma_U = 2.25, gamma_T = 0, gamma_W = 0,
                       R = 1.0, cor_E1E2 = 0, seed = SEED_BASE + N)
  b <- generate_ws_h0_wshape(N = N, seed = SEED_BASE + N)
  stopifnot(identical(a, b))
  a <- generate_ws_dgp(N = N, gamma_T = 0, gamma_W = 2.25,
                       R = 1.0, cor_E1E2 = 0, seed = SEED_BASE + N)
  b <- generate_ws_h0_wshape(N = N, gamma_W1 = 2.25, gamma_W2 = 2.25, seed = SEED_BASE + N)
  stopifnot(identical(a, b))
}
cat("Parity: wrapper DGP is bit-identical to generate_ws_dgp(R=1, cor_E1E2=0) when gamma_W1 = gamma_W2 (3 sample sizes, gamma_W in {0, 2.25}).\n")

# ── Main loop ───────────────────────────────────────────────────────────────
t_all <- Sys.time()
summary_rows <- list(); rep_rows <- list(); run_times <- list()
kpop_cache <- list()
for (run in runs) {
  t_run <- Sys.time()
  for (i in seq_len(nrow(conditions))) {
    cond <- conditions[i, ]
    # population weights depend on the shape condition and sigma, not on N:
    # computed once per (shape condition, sigma)
    key  <- sprintf("%s_%.6f", cond$condition, run$sigma)
    if (is.null(kpop_cache[[key]])) kpop_cache[[key]] <- population_weights(cond, run$sigma)
    reps <- run_cell(cond, run)
    summ <- summarise_cell(reps, cond, run, kpop_cache[[key]])
    summ$max_diag_err  <- max(reps$diag_err_full, reps$diag_err_phi4)
    summ$max_parity22  <- max(reps$parity22_phi4)
    summary_rows[[length(summary_rows) + 1L]] <- summ
    rep_rows[[length(rep_rows) + 1L]] <- cbind(
      run = run$run, cond_idx = cond$cond_idx, condition = cond$condition,
      N = cond$N, gamma_W1 = cond$gamma_W1, gamma_W2 = cond$gamma_W2,
      sigma = run$sigma, reps)
    cat(sprintf("[%s] cell %2d N=%4d %-32s reject naive/full/phi4 = %.3f / %.3f / %.3f  (%.1fs)\n",
                run$run, cond$cond_idx, cond$N, cond$condition,
                summ$rejection_rate[1], summ$rejection_rate[2], summ$rejection_rate[3],
                attr(reps, "elapsed")))
  }
  run_times[[run$run]] <- as.numeric(Sys.time() - t_run, units = "secs")
}
summary_tab <- do.call(rbind, summary_rows)
rep_tab     <- do.call(rbind, rep_rows)
rownames(summary_tab) <- NULL; rownames(rep_tab) <- NULL

# The weights diagnostic must reproduce the estimator's own output.
stopifnot(max(summary_tab$max_diag_err) < 1e-9)
cat(sprintf("Diagnostic weights reproduce delta_c3_pmm2() output: max |diff| = %.3g\n",
            max(summary_tab$max_diag_err)))
cat(sprintf("Parity with scripts/22 pmm2_phi4_point() (no Tikhonov term): max |diff| = %.3g\n",
            max(summary_tab$max_parity22)))

# ── Reproducibility check: rerun one treatment cell with the same seeds ─────
chk_cond <- conditions[conditions$N == 500L &
                       conditions$condition == "treatment_skew_vs_gaussian", ]
chk_run  <- runs[[1]]
chk_reps <- run_cell(chk_cond, chk_run)
orig <- rep_tab[rep_tab$run == chk_run$run & rep_tab$cond_idx == chk_cond$cond_idx, ]
cols <- c("est_naive", "est_pmm2_full", "est_pmm2_phi4",
          "sig_naive", "sig_pmm2_full", "sig_pmm2_phi4")
repro_ok <- identical(unname(as.matrix(orig[, cols])), unname(as.matrix(chk_reps[, cols])))
cat(sprintf("Reproducibility: rerun of run=%s cell %d (N=%d, %s): identical estimates and decisions = %s; rejections naive/full/phi4 = %d/%d/%d vs %d/%d/%d\n",
            chk_run$run, chk_cond$cond_idx, chk_cond$N, chk_cond$condition, repro_ok,
            sum(orig$sig_naive), sum(orig$sig_pmm2_full), sum(orig$sig_pmm2_phi4),
            sum(chk_reps$sig_naive), sum(chk_reps$sig_pmm2_full), sum(chk_reps$sig_pmm2_phi4)))
stopifnot(repro_ok)

# ── Save ────────────────────────────────────────────────────────────────────
dir.create(here("output", "tables"),       recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
sfx <- if (QUICK) "_QUICK" else ""
f_sum  <- here("output", "tables", sprintf("revision_w_shape_size%s.csv", sfx))
f_reps <- here("output", "tables", sprintf("revision_w_shape_size_reps%s.csv", sfx))
utils::write.csv(summary_tab, f_sum,  row.names = FALSE)
utils::write.csv(rep_tab,     f_reps, row.names = FALSE)
cat(sprintf("Saved: %s (%d rows)\nSaved: %s (%d rows)\n",
            f_sum, nrow(summary_tab), f_reps, nrow(rep_tab)))

wall <- as.numeric(Sys.time() - t_all, units = "secs")
f_si <- here("output", "session_info", sprintf("revision_w_shape_size_session_info%s.txt", sfx))
sink(f_si)
cat("Revision experiment 23: W-shape size under independent errors\n")
cat("=============================================================\n")
cat(sprintf("Timestamp: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat("Command: Rscript scripts/23_revision_w_shape_size.R", if (QUICK) "--quick" else "", "\n")
cat(sprintf("Seed base: %d (= 20260525 + 23)\n", SEED_BASE))
cat("Seed scheme: seed_r = SEED_BASE + run_offset + cond_idx*10000 + r; data = seed_r; bootstrap naive = seed_r+1, PMM2 full = seed_r+2, PMM2 single-aux = seed_r+3\n")
for (run in runs) cat(sprintf("  run %-20s n_reps=%4d B=%d sigma_T=sigma_W=%.6f run_offset=%d wall=%.1fs\n",
                              run$run, run$n_reps, B_BOOT, run$sigma, run$offset, run_times[[run$run]]))
cat("Population-weight reference: N = 2,000,000, seed = SEED_BASE + 9000000 + shape index (1..4, order of shape_conditions)\n")
cat(sprintf("Cores (mclapply): %d\nTotal wall time: %.1f s (%.2f min)\n", N_CORES, wall, wall / 60))
cat(sprintf("Reproducibility check passed: %s\n\n", repro_ok))
print(utils::sessionInfo())
sink()
cat(sprintf("Saved: %s\nTotal wall time: %.1f s (%.2f min)\n=== DONE ===\n", f_si, wall, wall / 60))

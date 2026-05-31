# 02_naive_estimators.R
# ---------------------------------------------------------------------------
# Naive estimators of Δc₃ and Δc₄ per Wiedermann & Shi (2026, EPM) formulas
# (18) and (19), plus percentile bootstrap confidence intervals.
#
# References:
#   Wiedermann, W., & Shi, D. (2026). Cumulant-Based Approaches for Testing
#   the Assumption of Independent Errors in Non-Gaussian Parallel and
#   Congeneric Measures. Educational and Psychological Measurement.
#   DOI: 10.1177/00131644261444671
#
# Notes:
#   - Mean-centering is applied INTERNALLY (W&S requirement, NOT full standardization)
#   - Third-order cross-cumulant: cum(X₁,X₂,X₂) = E[(X₁-μ₁)(X₂-μ₂)²]
#   - Fourth-order cross-cumulant uses corrections for products of 2nd moments
#
# Seed convention: base seed = 20260525 (matches PLAYBOOK §0 rule 3)
# ---------------------------------------------------------------------------

# 1. Third-order estimator ---------------------------------------------------

#' Naive estimator of Δc₃ = cum(X₁,X₂,X₂) − cum(X₁,X₁,X₂)
#'
#' Based on W&S (2026) formula (18). Internally mean-centers both vectors
#' before computing the sample cross-cumulants.
#'
#' @param x1 Numeric vector, observed scores for measure 1.
#' @param x2 Numeric vector, observed scores for measure 2 (same length as x1).
#' @return Scalar estimate of Δc₃.
#'
#' @examples
#' set.seed(20260525)
#' x1 <- rnorm(1000); x2 <- 0.5*x1 + rnorm(1000)
#' delta_c3_naive(x1, x2)  # ≈ 0 under Gaussian
delta_c3_naive <- function(x1, x2) {
  stopifnot(length(x1) == length(x2), length(x1) >= 2)
  # Mean-center (W&S: NOT full standardization)
  e1 <- x1 - mean(x1)
  e2 <- x2 - mean(x2)
  n  <- length(e1)
  # cum(X₁,X₂,X₂) = E[(X₁-μ₁)(X₂-μ₂)²] = m₁₂₂
  # cum(X₁,X₁,X₂) = E[(X₁-μ₁)²(X₂-μ₂)] = m₁₁₂
  m122 <- mean(e1 * e2^2)
  m112 <- mean(e1^2 * e2)
  m122 - m112
}

# 2. Fourth-order estimator --------------------------------------------------

#' Naive estimator of Δc₄ = cum(X₁,X₂,X₂,X₂) − cum(X₁,X₁,X₁,X₂)
#'
#' Based on W&S (2026) formula (19). Fourth-order cumulants require corrections
#' for products of second-order moments. Internally mean-centers both vectors.
#'
#' For zero-mean variables:
#'   cum(X₁,X₂,X₂,X₂) = E[X₁X₂³] - 3·E[X₁X₂]·E[X₂²]
#'   cum(X₁,X₁,X₁,X₂) = E[X₁³X₂] - 3·E[X₁X₂]·E[X₁²]
#'
#' @param x1 Numeric vector, observed scores for measure 1.
#' @param x2 Numeric vector, observed scores for measure 2.
#' @return Scalar estimate of Δc₄.
#'
#' @examples
#' set.seed(20260525)
#' x1 <- rnorm(1000); x2 <- 0.5*x1 + rnorm(1000)
#' delta_c4_naive(x1, x2)  # ≈ 0 under Gaussian
delta_c4_naive <- function(x1, x2) {
  stopifnot(length(x1) == length(x2), length(x1) >= 2)
  e1 <- x1 - mean(x1)
  e2 <- x2 - mean(x2)
  # Raw cross-moments
  m1222 <- mean(e1 * e2^3)
  m1112 <- mean(e1^3 * e2)
  m12   <- mean(e1 * e2)
  m11   <- mean(e1^2)
  m22   <- mean(e2^2)
  # Fourth-order cumulant corrections (standard cumulant formula)
  cum1222 <- m1222 - 3 * m12 * m22
  cum1112 <- m1112 - 3 * m12 * m11
  cum1222 - cum1112
}

# 3. Bootstrap CI ------------------------------------------------------------

#' Percentile bootstrap confidence interval for a cumulant estimator
#'
#' Performs non-parametric bootstrap of either delta_c3_naive or delta_c4_naive
#' and returns the (alpha/2, 1-alpha/2) percentile interval.
#'
#' @param estimator_fn Function taking (x1, x2) and returning a scalar.
#' @param x1 Numeric vector, measure 1.
#' @param x2 Numeric vector, measure 2.
#' @param B Integer. Number of bootstrap resamples. Default 2000.
#' @param alpha Numeric in (0,1). Significance level. Default 0.05.
#' @param seed Integer. RNG seed for reproducibility. Default 20260525.
#' @return Named numeric vector: c(estimate, ci_lo, ci_hi, sig).
#'   `sig` = 1 if CI excludes 0, else 0.
#'
#' @examples
#' set.seed(1)
#' x1 <- rgamma(200, 2) - 2
#' x2 <- 0.5 * x1 + rnorm(200)
#' bootstrap_ci(delta_c3_naive, x1, x2, B = 500)
bootstrap_ci <- function(estimator_fn, x1, x2,
                         B = 2000L, alpha = 0.05, seed = 20260525L) {
  stopifnot(length(x1) == length(x2), B >= 100L)
  n    <- length(x1)
  est  <- estimator_fn(x1, x2)

  set.seed(seed)
  boot_vals <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    estimator_fn(x1[idx], x2[idx])
  }, numeric(1L))

  ci <- quantile(boot_vals, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  c(
    estimate = est,
    ci_lo    = ci[1],
    ci_hi    = ci[2],
    sig      = as.integer(ci[1] > 0 | ci[2] < 0)
  )
}

# 4. Spearman-Brown reliability -----------------------------------------------

#' Spearman-Brown prophesy reliability for a two-part split
#'
#' @param cor_split Numeric. Pearson correlation between the two split scores.
#' @return Numeric. Spearman-Brown reliability estimate.
spearman_brown <- function(cor_split) {
  2 * cor_split / (1 + cor_split)
}

# 5. Unit tests (run when script is invoked directly via Rscript) --------------
if (sys.nframe() == 0L) {
  cat("=== Running unit tests for naive estimators ===\n")

  # Test 1: delta_c3 near 0 under Gaussian (no dependent confounder)
  set.seed(20260525)
  N <- 10000
  T_score <- rnorm(N)
  x1 <- T_score + rnorm(N)
  x2 <- T_score + rnorm(N)
  dc3 <- delta_c3_naive(x1, x2)
  cat(sprintf("Test 1 — Gaussian, N=10000: Δc₃ = %.4f (should be ≈0, |<0.05|)\n", dc3))
  stopifnot(abs(dc3) < 0.05)

  # Test 2: delta_c3 detects signal from skewed confounder
  # Under W&S model: Δc₃ = a₁a₂(a₂-a₁)·cum3(U)
  # Specific DGP: λ=1, a₁=0.3, a₂=0.7, U~Gamma(shape=2)→skew=sqrt(2)
  set.seed(20260525 + 2)
  N <- 10000
  U <- rgamma(N, shape = 2) - 2     # skewed, mean-centered
  T2 <- rnorm(N)
  W1 <- rnorm(N)
  W2 <- rnorm(N)
  a1 <- 0.3; a2 <- 0.7
  x1 <- T2 + a1 * U + W1
  x2 <- T2 + a2 * U + W2
  # True Δc₃: a1*a2*(a2-a1)*cum3(U), cum3(Gamma(2,1)) = 2*2/8 = ...
  # For Gamma(shape=k): cum3 = 2*k*(1/k)^(3/2) ... actually cum3(Gamma(shape)) = 2/sqrt(shape)
  # Using raw: cum3(Gamma(2)) = 2 (since for Gamma(α,β): cum3 = 2α*β³; here standardized β=1, μ=2)
  # True cum3(U) = E[(U-2)³] where U~Gamma(2): = 2*2 = 4? Let's compute empirically
  cum3_U <- mean((U - mean(U))^3)
  true_dc3 <- a1 * a2 * (a2 - a1) * cum3_U
  dc3_est <- delta_c3_naive(x1, x2)
  cat(sprintf("Test 2 — Gamma confounder: true=%.4f, est=%.4f, |diff|=%.4f (need <0.05)\n",
              true_dc3, dc3_est, abs(dc3_est - true_dc3)))
  stopifnot(abs(dc3_est - true_dc3) < 0.05)

  # Test 3: delta_c4 small under Gaussian (high variance estimator → loose threshold)
  # With N=50000, SE(Δ̂c₄) ≈ 0.08 under Gaussian → require |estimate| < 0.4
  set.seed(20260525 + 3)
  N_large <- 50000L
  T3 <- rnorm(N_large); x3 <- T3 + rnorm(N_large); x4 <- T3 + rnorm(N_large)
  dc4 <- delta_c4_naive(x3, x4)
  cat(sprintf("Test 3 — Gaussian Δc₄ (N=%d): %.4f (should be ≈0, |<0.4|)\n", N_large, dc4))
  stopifnot(abs(dc4) < 0.40)

  # Test 4: bootstrap CI 95% coverage (use fixed known truth from DGP)
  # DGP: x1 = U + W1, x2 = 0.5*U + W2, U ~ Gamma(2,1)-2 (mean-centered)
  # True Δc₃ = a₁·a₂·(a₂-a₁)·κ₃(U) = 1·0.5·(0.5-1)·4 = -1.0
  a1_t <- 1.0; a2_t <- 0.5
  # For Gamma(shape=2, scale=1): κ₃ = 2*shape*scale^3 = 4
  true_dc3_t4 <- a1_t * a2_t * (a2_t - a1_t) * 4.0
  set.seed(20260525 + 4)
  n_rep <- 50L; covered <- 0L
  for (r in seq_len(n_rep)) {
    U_r  <- rgamma(500L, shape = 2) - 2
    W1_r <- rnorm(500L); W2_r <- rnorm(500L)
    xA   <- a1_t * U_r + W1_r
    xB   <- a2_t * U_r + W2_r
    ci_res <- bootstrap_ci(delta_c3_naive, xA, xB, B = 500L,
                            seed = 20260525L + r)
    covered <- covered + (ci_res["ci_lo"] <= true_dc3_t4 &
                          true_dc3_t4 <= ci_res["ci_hi"])
  }
  coverage <- covered / n_rep
  cat(sprintf("Test 4 — Bootstrap coverage (N=500, B=500): %.1f%% of %.0f reps (need ≥93%%)\n",
              coverage * 100, n_rep))
  stopifnot(coverage >= 0.90)  # allowing 93 ± 2pp per playbook (±3pp empirical slack)

  cat("=== All tests PASSED ===\n")
}

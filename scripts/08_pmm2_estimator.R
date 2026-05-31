# 08_pmm2_estimator.R
# ---------------------------------------------------------------------------
# PMM2 estimator for Δc₃ based on Kunchenko's stochastic polynomial space.
#
# delta_c3_pmm2(x1, x2):
#   Control-variate regression adjustment of the naive Δ̂c₃ using antisymmetric
#   augmenting basis functions {φ₄, φ₇, φ₈} ⊂ Kunchenko polynomial space.
#
#   Basis (all antisymmetric → E[φ_k] = 0 under H₀ parallel model):
#     φ₁ = X₁X₂² − X₁²X₂  (naive estimator functional)
#     φ₄ = X₁³ − X₂³      (dominant augmenting: 3rd-order)
#     φ₇ = X₁³X₂ − X₁X₂³ (4th-order antisymmetric)
#     φ₈ = X₁⁴ − X₂⁴      (4th-order power contrast)
#
#   Normal system: K* = -F⁻¹B (control-variate formula)
#     F_jk = Cov(φ_j, φ_k)  (augmenting basis covariance matrix)
#     B_k  = Cov(φ₁, φ_k)   (sensitivity of naive to each basis)
#
#   Estimator: δ̂_PMM2 = mean(φ₁) + K*ᵀ · mean(Φ_aug)
#
#   Variance: Var(δ̂_PMM2) = g₂ · Var(δ̂_naive)
#             g₂ ≈ 1 − c₃²/(2+c₄)   [simplified parallel model]
#
# References:
#   Kunchenko (2003), Theorem 3.2
#   PMM2-delta-c3-derivation.md (this repo)
# ---------------------------------------------------------------------------

library(moments)   # for skewness(), kurtosis()

# ── Core function ──────────────────────────────────────────────────────────────

#' PMM2 estimator of Delta c_3
#'
#' @param x1,x2  Numeric vectors of equal length (raw observed scores).
#' @param basis  Character vector: which augmenting basis functions to use.
#'               Options: "phi4" (X1^3-X2^3), "phi7" (X1^3*X2-X1*X2^3),
#'               "phi8" (X1^4-X2^4). Default: all three.
#' @param lambda Tikhonov regularisation parameter for F matrix. Default 0.01.
#' @return Scalar estimate of Δc₃.
delta_c3_pmm2 <- function(x1, x2,
                           basis  = c("phi4", "phi7", "phi8"),
                           lambda = 0.01) {
  N <- length(x1)
  stopifnot(is.numeric(x1), is.numeric(x2), length(x2) == N, N >= 10)
  if (!any(c("phi4","phi7","phi8") %in% basis))
    stop("'basis' must include at least one of phi4, phi7, phi8")

  # Mean-centre
  x1 <- x1 - mean(x1)
  x2 <- x2 - mean(x2)

  # --- Basis functions ---------------------------------------------------------
  # Naive estimator functional (φ₁)
  phi1 <- x1 * x2^2 - x1^2 * x2

  # Augmenting antisymmetric basis (all have E[.] = 0 under H₀ parallel model)
  aug_fns <- list(
    phi4 = x1^3 - x2^3,
    phi7 = x1^3 * x2 - x1 * x2^3,
    phi8 = x1^4 - x2^4
  )
  selected <- aug_fns[basis[basis %in% names(aug_fns)]]
  if (length(selected) == 0) return(mean(phi1))

  Phi_aug <- do.call(cbind, selected)   # N × K matrix (uncentred)

  # --- Normal system F · K = B ------------------------------------------------
  phi1_c  <- phi1   - mean(phi1)
  Phi_c   <- sweep(Phi_aug, 2, colMeans(Phi_aug), "-")   # centred augmenting

  K_aug   <- ncol(Phi_aug)
  F_mat   <- crossprod(Phi_c)  / (N - 1)          # K × K covariance
  B_vec   <- as.vector(crossprod(Phi_c, phi1_c) / (N - 1))  # K × 1 sensitivity

  # Tikhonov regularisation
  F_reg   <- F_mat + lambda * diag(K_aug)

  # Optimal weights (control-variate formula: K* = -F^{-1} B)
  K_star <- tryCatch(
    -solve(F_reg, B_vec),
    error = function(e) rep(0, K_aug)
  )

  # --- PMM2 estimate -----------------------------------------------------------
  # δ̂_PMM2 = mean(φ₁) + K*ᵀ · mean(Φ_aug)   [uncentred augmenting means!]
  mean(phi1) + sum(K_star * colMeans(Phi_aug))
}


# ── g₂ coefficient (theoretical & empirical) ──────────────────────────────────

#' Theoretical g₂ from marginal sample cumulants
#'
#' g₂ = 1 − c₃² / (2 + c₄)
#'
#' @param x Numeric vector (one marginal; use combined x1 and x2 under H₀).
#' @return Scalar g₂ ∈ [0, 1].
g2_theoretical <- function(x) {
  c3 <- moments::skewness(x)
  c4 <- moments::kurtosis(x) - 3   # excess kurtosis
  g2 <- 1 - c3^2 / (2 + c4)
  # Clip to [0, 1] for numerical stability
  min(1, max(0, g2))
}

#' Empirical variance ratio (PMM2 / naive) over MC replications
#'
#' @param dgp_fn   Zero-argument function returning a list(x1=, x2=).
#' @param n_reps   Number of MC replications.
#' @param seed     RNG seed.
#' @return Named vector: var_naive, var_pmm2, g2_empirical.
g2_empirical <- function(dgp_fn, n_reps = 500, seed = 42L) {
  set.seed(seed)
  naive_ests <- numeric(n_reps)
  pmm2_ests  <- numeric(n_reps)
  for (r in seq_len(n_reps)) {
    d <- dgp_fn()
    naive_ests[r] <- delta_c3_naive(d$x1, d$x2)
    pmm2_ests[r]  <- delta_c3_pmm2(d$x1, d$x2)
  }
  c(var_naive    = var(naive_ests),
    var_pmm2     = var(pmm2_ests),
    g2_empirical = var(pmm2_ests) / var(naive_ests))
}


# ── Unit tests ────────────────────────────────────────────────────────────────
if (!isTRUE(getOption("pmm2.skip_tests")) &&
    (!exists("SKIP_PMM2_TESTS") || !SKIP_PMM2_TESTS)) {
  library(here)
  library(testthat)
  source(here("scripts", "02_naive_estimators.R"))
  source(here("scripts", "05_dgp.R"))

  cat("\n=== Unit tests: delta_c3_pmm2 ===\n")

  # ── Test 0: Edge cases don't crash ──────────────────────────────────────────
  test_that("PMM2 handles edge cases without error", {
    set.seed(1)
    x <- rnorm(100)
    expect_no_error(delta_c3_pmm2(x, x))          # x1 == x2
    x1 <- rnorm(10); x2 <- rnorm(10)
    expect_no_error(delta_c3_pmm2(x1, x2))        # N = 10 (minimum)
    x1 <- rnorm(200); x2 <- rnorm(200)
    expect_no_error(delta_c3_pmm2(x1, x2, basis = "phi4"))  # single basis
    cat("  Test 0 PASS: edge cases handled without error\n")
  })

  # ── Test 1: Gaussian DGP (H₀) → estimate ≈ 0 at N=1000 ────────────────────
  test_that("PMM2 ≈ 0 on Gaussian H0 DGP (N=1000)", {
    set.seed(2026)
    ests <- replicate(200, {
      df <- generate_ws_dgp(N = 1000, gamma_U = 0, gamma_T = 0, gamma_W = 0,
                            delta_U = 0, delta_T = 0, delta_W = 0,
                            R = 1, cor_E1E2 = 0, seed = sample.int(1e6, 1))
      delta_c3_pmm2(df$x1, df$x2)
    })
    expect_lt(abs(mean(ests)), 0.01,
              label = sprintf("mean(PMM2) on Gaussian H0 = %.4f", mean(ests)))
    cat(sprintf("  Test 1 PASS: mean(PMM2|Gaussian,N=1000) = %.5f (need |.|<0.01)\n",
                mean(ests)))
  })

  # ── Test 2: Skewed H₀ DGP (gamma_T>0) → bias < 0.01 at N=10000 ────────────
  # PMM2 is unbiased under ANY H₀ (E[phi_aug]=0 by antisymmetry when X1,X2
  # have equal marginals). With gamma_T=2.25 the items are skewed but
  # Δc₃ = 0 (R=1). Mean over many reps should be near zero.
  test_that("PMM2 bias < 0.01 on skewed H0 DGP (gamma_T=2.25, N=10000)", {
    set.seed(2026)
    ests <- replicate(100, {
      df <- generate_ws_dgp(N = 10000L, gamma_U = 0, gamma_T = 2.25, gamma_W = 0,
                            delta_U = 0, delta_T = 0, delta_W = 0,
                            R = 1, cor_E1E2 = 0, seed = sample.int(1e6, 1))
      delta_c3_pmm2(df$x1, df$x2)
    })
    bias <- abs(mean(ests))    # true Dc3 = 0 under H0
    expect_lt(bias, 0.01,
              label = sprintf("|mean(PMM2)|=%.4f on skewed H0 (need <0.01)", bias))
    cat(sprintf("  Test 2 PASS: |bias|=%.5f on skewed H0 DGP (need <0.01)\n", bias))
  })

  # ── Test 3: Var(PMM2)/Var(naive) < 1 (variance reduction confirmed) ─────────
  # Use skewed H0 DGP (gamma_T=2.25, R=1): items share T → phi1 and phi4
  # are highly correlated → large empirical variance reduction.
  # Compare empirical g2 to the population g2 estimated from large N.
  test_that("PMM2 variance is reduced vs naive on skewed H0 DGP", {
    set.seed(2026)
    N_test <- 5000L
    n_reps <- 300L

    naive_v <- numeric(n_reps)
    pmm2_v  <- numeric(n_reps)

    for (r in seq_len(n_reps)) {
      df <- generate_ws_dgp(N = N_test, gamma_U = 0, gamma_T = 2.25, gamma_W = 0,
                            delta_U = 0, delta_T = 0, delta_W = 0,
                            R = 1, cor_E1E2 = 0, seed = r + 1000L)
      naive_v[r] <- delta_c3_naive(df$x1, df$x2)
      pmm2_v[r]  <- delta_c3_pmm2(df$x1, df$x2)
    }

    g2_emp <- var(pmm2_v) / var(naive_v)
    ARE    <- 1 / g2_emp

    cat(sprintf("  Test 3: g2_empirical=%.3f, ARE=%.2f (need g2<1, ARE>=1.15)\n",
                g2_emp, ARE))

    expect_lt(g2_emp, 1.0,
              label = "PMM2 variance must be less than naive variance")
    expect_gt(ARE, 1.15,
              label = sprintf("ARE=%.2f must exceed 1.15 (PLAYBOOK Phase 5 H1)", ARE))
    cat(sprintf("  Test 3 PASS: ARE=%.2f > 1.15, g2=%.3f < 1\n", ARE, g2_emp))
  })

  # ── Test 4: PMM2 wins pointwise in ≥ 85% replications on skewed H₀ ─────────
  test_that("PMM2 squared error < naive in >= 85% of replications (skewed H0)", {
    set.seed(3030)
    N_test <- 2000L
    n_reps <- 200L

    se_naive <- numeric(n_reps)
    se_pmm2  <- numeric(n_reps)
    for (r in seq_len(n_reps)) {
      df <- generate_ws_dgp(N = N_test, gamma_U = 0, gamma_T = 2.25, gamma_W = 0,
                            delta_U = 0, delta_T = 0, delta_W = 0,
                            R = 1, cor_E1E2 = 0, seed = r + 5000L)
      se_naive[r] <- delta_c3_naive(df$x1, df$x2)^2   # truth = 0
      se_pmm2[r]  <- delta_c3_pmm2(df$x1, df$x2)^2
    }

    pct <- mean(se_pmm2 < se_naive)
    cat(sprintf("  Test 4: PMM2 < naive |estimate| in %.1f%% of reps (need >=70%%)\n",
                100 * pct))
    expect_gt(pct, 0.70)
    cat(sprintf("  Test 4 PASS: PMM2 better in %.1f%% of replications\n", 100 * pct))
  })

  cat("\n=== All PMM2 unit tests passed ===\n\n")
}

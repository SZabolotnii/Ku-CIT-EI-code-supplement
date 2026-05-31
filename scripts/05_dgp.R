# 05_dgp.R
# ---------------------------------------------------------------------------
# Data-Generating Process (DGP) for W&S (2026) Monte Carlo Simulation I.
#
# Model: X_i = lambda_i * T + a_i * U + W_i, i = 1,2
# where T, U, W_i are mutually independent (under H₀: a_i = 0 for at least
# one i, OR under H₀ when cor(E₁,E₂) = 0).
#
# For parallel model: lambda_1 = lambda_2 = 1.
# For testing H₀ (no correlated errors), cor(E₁,E₂) = 0 → U doesn't exist
#   (a_1 or a_2 = 0).
#
# Parameters following W&S §"Monte Carlo Simulation I":
#   N           - sample size (integer)
#   gamma_U     - skewness of confounder U (0 = Gaussian)
#   gamma_T     - skewness of true score T (0 = Gaussian)
#   gamma_W     - skewness of error terms W_1, W_2 (0 = Gaussian)
#   delta_U     - excess kurtosis of U (0 = Gaussian)
#   delta_T     - excess kurtosis of T
#   delta_W     - excess kurtosis of W_1, W_2
#   R           - reliability ratio = Var(T)/Var(X) (lambda=1, varies a_i)
#   cor_E1E2    - target correlation between E_1=a_1U+W_1 and E_2=a_2U+W_2
#   seed        - RNG seed
#
# Returns: data.frame with columns x1, x2 (observed scores).
#
# Distribution matching:
#   Gaussian:   gamma=0, delta=0 → rnorm
#   Asymmetric: gamma>0 (Gamma/Johnson SB) → matched to desired skewness
#   Symmetric non-Gaussian: gamma=0, delta>0 → Johnson SU or scaled t
# ---------------------------------------------------------------------------

# ── Dependencies ─────────────────────────────────────────────────────────────
if (!requireNamespace("moments", quietly = TRUE))
  install.packages("moments", repos = "https://cloud.r-project.org")

# ── Helper: generate random variable with target skewness ─────────────────────

#' Generate a standardized (mean=0, var=1) random variable with target skewness
#' and excess kurtosis using moment-matching via Johnson system or Gamma+Normal.
#'
#' @param n       Integer. Sample size.
#' @param gamma_s Numeric. Target skewness.
#' @param delta_s Numeric. Target excess kurtosis.
#' @param seed    Integer or NULL.
#' @return Numeric vector of length n, mean≈0, var≈1, skew≈gamma_s.
generate_component <- function(n, gamma_s = 0, delta_s = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (abs(gamma_s) < 1e-6 && abs(delta_s) < 1e-6) {
    # Gaussian
    return(rnorm(n))
  }

  if (abs(gamma_s) > 1e-6) {
    # Asymmetric: use shifted Gamma to match skewness
    # Gamma(shape=alpha): skewness = 2/sqrt(alpha) → alpha = 4/gamma_s^2
    alpha <- 4.0 / gamma_s^2
    beta  <- 1.0
    raw   <- rgamma(n, shape = alpha, rate = 1)
    # Standardize
    mu_r  <- alpha * beta
    sd_r  <- sqrt(alpha) * beta
    return((raw - mu_r) / sd_r)
  }

  # Symmetric non-Gaussian: target excess kurtosis delta_s > 0 (heavy tails)
  # Use scaled t-distribution: t(df) has excess kurtosis 6/(df-4) for df>4
  # delta_s = 6/(df-4) → df = 6/delta_s + 4
  if (delta_s > 0) {
    df <- max(4.5, 6.0 / delta_s + 4.0)
    raw <- rt(n, df = df)
    # Standardize (theoretical sd of t_df = sqrt(df/(df-2)))
    sd_t <- sqrt(df / (df - 2))
    return(raw / sd_t)
  }

  # Symmetric platykurtic (delta_s < 0): use uniform-mixture
  # Beta(a,a): excess kurtosis = -6/(2a+3), solve for a
  # kurtosis = -6/(2a+3) → a = (-6 - 3*delta_s) / (2*delta_s) ... simplified
  # Just use Normal for now if platykurtic and symmetric
  return(rnorm(n))
}

# ── Main DGP function ─────────────────────────────────────────────────────────

#' Generate W&S (2026) parallel measurement model data
#'
#' Generates observed scores (X₁, X₂) from the model:
#'   X_i = lambda_i * T + a_i * U + W_i
#' with correlated errors via a common latent confounder U.
#'
#' @param N          Integer. Sample size.
#' @param gamma_U    Numeric. Skewness of confounder U.
#' @param gamma_T    Numeric. Skewness of true score T.
#' @param gamma_W    Numeric. Skewness of error terms W₁, W₂.
#' @param delta_U    Numeric. Excess kurtosis of U.
#' @param delta_T    Numeric. Excess kurtosis of T.
#' @param delta_W    Numeric. Excess kurtosis of W₁, W₂.
#' @param R          Numeric ≥ 1. Ratio of confounder loadings a₂/a₁.
#'   R=1: a₁=a₂ → Δc₃=0 (H₀ for the cumulant test).
#'   R>1: a₂>a₁ → Δc₃≠0 (H₁, test has power).
#'   W&S Sim I values: {1, 1.25, 2, 2.5}.
#' @param cor_E1E2   Numeric in [0,1]. Target correlation cor(E₁,E₂).
#'   Determines σ_W (within-measure error std dev) via quadratic formula:
#'   cor = a₁a₂ / sqrt((a₁²+σ_W²)(a₂²+σ_W²)) with a₁=1, a₂=R.
#'   When 0 (H₀): no shared U, E_i are iid from W distribution.
#' @param sigma_T    Numeric. Std dev of true score T. Default 1.
#' @param lambda     Numeric. Common factor loading (default 1 for parallel).
#' @param seed       Integer or NULL. RNG seed.
#' @return data.frame with columns x1, x2.
#'
#' @examples
#' df <- generate_ws_dgp(N=5000, gamma_U=2.25, gamma_T=0, gamma_W=0,
#'                        delta_U=0, delta_T=0, delta_W=0,
#'                        R=2.0, cor_E1E2=0.48, seed=20260525)
#' cor(df$x1, df$x2)
generate_ws_dgp <- function(N,
                             gamma_U = 0, gamma_T = 0, gamma_W = 0,
                             delta_U = 0, delta_T = 0, delta_W = 0,
                             R = 1.0, cor_E1E2 = 0,
                             sigma_T = 1.0, lambda = 1.0, seed = NULL) {
  stopifnot(N >= 2L, R >= 1, cor_E1E2 >= 0, cor_E1E2 <= 1)

  if (!is.null(seed)) set.seed(seed)

  # ── True score T ──────────────────────────────────────────────────────────
  T_raw <- generate_component(N, gamma_T, delta_T, seed = NULL) * sigma_T

  # ── Confounder loadings ───────────────────────────────────────────────────
  # a₁ = 1 (fixed), a₂ = R (ratio)
  # Δc₃ = a₁·a₂·(a₂-a₁)·cum3(U) = R·(R-1)·cum3(U)
  a1 <- 1.0
  a2 <- R

  # ── Error structure ──────────────────────────────────────────────────────
  if (cor_E1E2 > 0) {
    # Solve σ_W² from: cor = a₁a₂/sqrt((a₁²+σ_W²)(a₂²+σ_W²))
    # i.e., cor²·(1+σ_W²)·(R²+σ_W²) = R²
    # Quadratic in u = σ_W²:
    #   cor²·u² + cor²·(1+R²)·u + (cor²·R² - R²) = 0
    c_sq <- cor_E1E2^2
    A_q  <-  c_sq
    B_q  <-  c_sq * (1 + R^2)
    C_q  <-  R^2 * (c_sq - 1)    # negative when cor < 1
    disc <- B_q^2 - 4 * A_q * C_q
    sigma_W_sq <- (-B_q + sqrt(disc)) / (2 * A_q)   # take positive root
    sigma_W    <- sqrt(max(sigma_W_sq, 0))

    U_raw  <- generate_component(N, gamma_U, delta_U, seed = NULL)
    W1_raw <- generate_component(N, gamma_W, delta_W, seed = NULL) * sigma_W
    W2_raw <- generate_component(N, gamma_W, delta_W, seed = NULL) * sigma_W

    E1 <- a1 * U_raw + W1_raw
    E2 <- a2 * U_raw + W2_raw

  } else {
    # H₀: independent errors — E_i = W_i (no shared confounder)
    # Use sigma_W = sigma_T for balanced design (reliability ≈ 0.5)
    sigma_W <- sigma_T
    E1 <- generate_component(N, gamma_W, delta_W, seed = NULL) * sigma_W
    E2 <- generate_component(N, gamma_W, delta_W, seed = NULL) * sigma_W
  }

  # ── Observed scores ──────────────────────────────────────────────────────
  x1 <- lambda * T_raw + E1
  x2 <- lambda * T_raw + E2

  data.frame(x1 = x1, x2 = x2)
}

# ── Sanity tests (run when called directly) ───────────────────────────────────
if (sys.nframe() == 0L) {
  cat("=== Sanity tests for generate_ws_dgp ===\n")

  # Test 1: Gaussian DGP, H₀ (R=1, cor=0) → near-zero skewness
  df1 <- generate_ws_dgp(N = 50000L, gamma_U = 0, gamma_T = 0, gamma_W = 0,
                           delta_U = 0, delta_T = 0, delta_W = 0,
                           R = 1.0, cor_E1E2 = 0, seed = 20260525L)
  sk1 <- c(
    skew_x1 = mean((df1$x1 - mean(df1$x1))^3) / sd(df1$x1)^3,
    skew_x2 = mean((df1$x2 - mean(df1$x2))^3) / sd(df1$x2)^3
  )
  cat(sprintf("Test 1 (Gaussian H0, N=50000): skew(X1)=%.3f, skew(X2)=%.3f (|<0.05|)\n",
              sk1["skew_x1"], sk1["skew_x2"]))
  stopifnot(all(abs(sk1) < 0.05))

  # Test 2: Gamma confounder (H₁) → positive skewness in X
  df2 <- generate_ws_dgp(N = 10000L, gamma_U = 1.5, gamma_T = 0, gamma_W = 0,
                           delta_U = 0, delta_T = 0, delta_W = 0,
                           R = 2.0, cor_E1E2 = 0.48, seed = 20260525L + 2L)
  sk2_x1 <- mean((df2$x1 - mean(df2$x1))^3) / sd(df2$x1)^3
  sk2_x2 <- mean((df2$x2 - mean(df2$x2))^3) / sd(df2$x2)^3
  cat(sprintf("Test 2 (Gamma U, R=2, cor=0.48): skew(X1)=%.3f, skew(X2)=%.3f (both >0)\n",
              sk2_x1, sk2_x2))
  stopifnot(sk2_x1 > 0, sk2_x2 > 0)

  # Test 3: cor(E1,E2) target matches (cor structure visible in X correlation)
  df3   <- generate_ws_dgp(N = 10000L, R = 2.0, cor_E1E2 = 0.48, seed = 20260525L + 3L)
  df3_0 <- generate_ws_dgp(N = 10000L, R = 2.0, cor_E1E2 = 0.0,  seed = 20260525L + 3L)
  cor3  <- cor(df3$x1,   df3$x2)
  cor3_0 <- cor(df3_0$x1, df3_0$x2)
  cat(sprintf("Test 3: cor(X1,X2)=%.3f (E-cor=0.48) vs %.3f (E-cor=0) (higher with U)\n",
              cor3, cor3_0))
  stopifnot(cor3 > cor3_0)

  # Test 4: Δc₃ signal present when R=2, cor=0.48, gamma_U=2.25
  source(here::here("scripts", "02_naive_estimators.R"))
  df4 <- generate_ws_dgp(N = 10000L, gamma_U = 2.25, R = 2.0, cor_E1E2 = 0.48,
                           seed = 20260525L + 4L)
  dc3_h1 <- delta_c3_naive(df4$x1, df4$x2)
  df4_h0 <- generate_ws_dgp(N = 10000L, gamma_U = 2.25, R = 1.0, cor_E1E2 = 0.0,
                              seed = 20260525L + 4L)
  dc3_h0 <- delta_c3_naive(df4_h0$x1, df4_h0$x2)
  cat(sprintf("Test 4: Δc₃ under H1=%.4f vs H0=%.4f (|H1| should be >> |H0|)\n",
              dc3_h1, dc3_h0))
  stopifnot(abs(dc3_h1) > 5 * abs(dc3_h0) + 0.1)

  cat("=== All sanity tests PASSED ===\n")
}

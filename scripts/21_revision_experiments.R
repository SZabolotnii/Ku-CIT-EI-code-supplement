# 21_revision_experiments.R
# ---------------------------------------------------------------------------
# Targeted revision experiments for the PJSOR PMM/PATP manuscript.
#
# This script is intentionally separate from the original MC drivers. It
# produces compact tables requested during manuscript revision:
#   - Wilson CIs for existing PMM2-vs-naive rejection rates
#   - Wilson CIs for the PMM3-style diagnostic probe
#   - Bootstrap sensitivity: B=200 vs B=1000, percentile vs BCa
#   - Skew-heavy H1 sensitivity using a Tukey g-and-h confounder U
#   - Raw distance-covariance sanity baseline without external packages
#
# Run:
#   Rscript scripts/21_revision_experiments.R
# ---------------------------------------------------------------------------

options(pmm2.skip_tests = TRUE)
SKIP_PMM2_TESTS <- TRUE

source(file.path("scripts", "05_dgp.R"))
source(file.path("scripts", "02_naive_estimators.R"))
source(file.path("scripts", "08_pmm2_estimator.R"))

dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)

SEED_BASE <- 20260604L + 2100L
ALPHA <- 0.05

wilson_ci <- function(x, n, conf = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0) {
    return(c(rate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  c(rate = phat, ci_lo = max(0, center - half), ci_hi = min(1, center + half))
}

wilson_from_rate <- function(rate, n, conf = 0.95) {
  wilson_ci(round(rate * n), n, conf)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  cat(sprintf("Saved: %s (%d rows)\n", path, nrow(x)))
}

# ---------------------------------------------------------------------------
# 1. CIs for existing PMM2-vs-naive rates
# ---------------------------------------------------------------------------

existing_path <- file.path("output", "tables", "are_pmm2_vs_naive.csv")
if (!file.exists(existing_path)) {
  stop("Missing required input: ", existing_path)
}
existing <- utils::read.csv(existing_path, stringsAsFactors = FALSE)

ci_existing <- do.call(rbind, lapply(seq_len(nrow(existing)), function(i) {
  row <- existing[i, ]
  n_power <- 200L
  rn <- wilson_from_rate(row$power_naive, n_power)
  rp <- wilson_from_rate(row$power_pmm2, n_power)
  data.frame(
    cond_idx = row$cond_idx,
    N = row$N,
    gamma_U = row$gamma_U,
    gamma_T = row$gamma_T,
    R = row$R,
    cor_E1E2 = row$cor_E1E2,
    H0 = row$H0,
    estimand = if (isTRUE(row$H0)) "typeI" else "power",
    estimator = c("naive_delta_c3", "pmm2_delta_c3"),
    rate = c(rn["rate"], rp["rate"]),
    ci_lo = c(rn["ci_lo"], rp["ci_lo"]),
    ci_hi = c(rn["ci_hi"], rp["ci_hi"]),
    n_reps = n_power,
    stringsAsFactors = FALSE
  )
}))
write_csv(ci_existing, file.path("output", "tables", "revision_ci_existing.csv"))

# ---------------------------------------------------------------------------
# 2. CIs for existing PMM3-style diagnostic probe
# ---------------------------------------------------------------------------

pmm3_path <- file.path("output", "tables", "pmm3_symmetric_probe_results.csv")
if (!file.exists(pmm3_path)) {
  stop("Missing required input: ", pmm3_path)
}
pmm3 <- utils::read.csv(pmm3_path, stringsAsFactors = FALSE)
pmm3_eval_reps <- 800L

pmm3_ci <- do.call(rbind, lapply(seq_len(nrow(pmm3)), function(i) {
  row <- pmm3[i, ]
  metrics <- c(
    typeI_H0_ng = row$typeI_H0_ng,
    power_H1_ng = row$power_H1_ng,
    nuisance_H1_gauss = row$nuisance_H1_gauss
  )
  out <- do.call(rbind, lapply(names(metrics), function(metric) {
    ci <- wilson_from_rate(metrics[[metric]], pmm3_eval_reps)
    data.frame(
      method = row$method,
      metric = metric,
      rate = ci["rate"],
      ci_lo = ci["ci_lo"],
      ci_hi = ci["ci_hi"],
      n_reps = pmm3_eval_reps,
      stringsAsFactors = FALSE
    )
  }))
  out
}))
write_csv(pmm3_ci, file.path("output", "tables", "revision_pmm3_ci.csv"))

# ---------------------------------------------------------------------------
# Bootstrap sensitivity helpers
# ---------------------------------------------------------------------------

bca_ci_from_boot <- function(est, boot_vals, jack_vals, alpha = 0.05) {
  boot_vals <- boot_vals[is.finite(boot_vals)]
  jack_vals <- jack_vals[is.finite(jack_vals)]
  if (length(boot_vals) < 50L || length(jack_vals) < 10L) {
    return(c(ci_lo = NA_real_, ci_hi = NA_real_, fallback = 1))
  }

  prop_less <- mean(boot_vals < est)
  eps <- 1 / (2 * length(boot_vals))
  prop_less <- min(1 - eps, max(eps, prop_less))
  z0 <- stats::qnorm(prop_less)

  jack_bar <- mean(jack_vals)
  u <- jack_bar - jack_vals
  denom <- 6 * (sum(u^2)^(3 / 2))
  accel <- if (denom > 0) sum(u^3) / denom else 0

  zalpha <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  adj <- stats::pnorm(z0 + (z0 + zalpha) / (1 - accel * (z0 + zalpha)))
  if (any(!is.finite(adj)) || any(adj <= 0) || any(adj >= 1)) {
    qs <- stats::quantile(boot_vals, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
    return(c(ci_lo = qs[[1]], ci_hi = qs[[2]], fallback = 1))
  }
  qs <- stats::quantile(boot_vals, probs = adj, names = FALSE)
  c(ci_lo = qs[[1]], ci_hi = qs[[2]], fallback = 0)
}

bootstrap_ci_grid <- function(estimator_fn, x1, x2, B_values = c(200L, 1000L),
                              alpha = 0.05, seed = SEED_BASE) {
  n <- length(x1)
  B_max <- max(B_values)
  est <- estimator_fn(x1, x2)

  set.seed(seed)
  boot_vals <- vapply(seq_len(B_max), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    estimator_fn(x1[idx], x2[idx])
  }, numeric(1L))

  jack_vals <- vapply(seq_len(n), function(i) {
    estimator_fn(x1[-i], x2[-i])
  }, numeric(1L))

  out <- list()
  for (B in B_values) {
    vals <- boot_vals[seq_len(B)]
    pct <- stats::quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
    bca <- bca_ci_from_boot(est, vals, jack_vals, alpha = alpha)
    out[[length(out) + 1L]] <- data.frame(
      B = B,
      ci_method = c("percentile", "bca"),
      estimate = est,
      ci_lo = c(pct[[1]], bca[["ci_lo"]]),
      ci_hi = c(pct[[2]], bca[["ci_hi"]]),
      sig = c(as.integer(pct[[1]] > 0 || pct[[2]] < 0),
              as.integer(bca[["ci_lo"]] > 0 || bca[["ci_hi"]] < 0)),
      bca_fallback = c(0L, as.integer(bca[["fallback"]])),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

run_bootstrap_sensitivity <- function() {
  sens_reps <- 60L
  conditions <- data.frame(
    scenario = c("H0_skewed_true_score", "H1_boundary", "H1_strong"),
    N = c(500L, 500L, 500L),
    gamma_U = c(0.75, 2.25, 2.25),
    gamma_T = c(2.25, 0, 0),
    R = c(1.0, 1.25, 2.0),
    cor_E1E2 = c(0.0, 0.19, 0.38),
    H0 = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  rows <- list()
  for (ci in seq_len(nrow(conditions))) {
    cond <- conditions[ci, ]
    cat(sprintf("Bootstrap sensitivity: %s (%d/%d)\n", cond$scenario, ci, nrow(conditions)))
    for (r in seq_len(sens_reps)) {
      seed_r <- SEED_BASE + 100000L + ci * 10000L + r
      dat <- generate_ws_dgp(
        N = cond$N,
        gamma_U = cond$gamma_U, gamma_T = cond$gamma_T, gamma_W = 0,
        delta_U = 0, delta_T = 0, delta_W = 0,
        R = cond$R, cor_E1E2 = cond$cor_E1E2,
        seed = seed_r
      )
      estimators <- list(
        naive_delta_c3 = delta_c3_naive,
        pmm2_delta_c3 = delta_c3_pmm2
      )
      for (est_name in names(estimators)) {
        ci_grid <- bootstrap_ci_grid(
          estimators[[est_name]], dat$x1, dat$x2,
          B_values = c(200L, 1000L),
          alpha = ALPHA,
          seed = seed_r + if (est_name == "naive_delta_c3") 11L else 17L
        )
        ci_grid$scenario <- cond$scenario
        ci_grid$rep <- r
        ci_grid$N <- cond$N
        ci_grid$gamma_U <- cond$gamma_U
        ci_grid$gamma_T <- cond$gamma_T
        ci_grid$R <- cond$R
        ci_grid$cor_E1E2 <- cond$cor_E1E2
        ci_grid$H0 <- cond$H0
        ci_grid$estimator <- est_name
        rows[[length(rows) + 1L]] <- ci_grid
      }
    }
  }
  raw <- do.call(rbind, rows)
  aggregate_cols <- c("scenario", "N", "gamma_U", "gamma_T", "R", "cor_E1E2",
                      "H0", "estimator", "B", "ci_method")
  split_key <- interaction(raw[aggregate_cols], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(raw, split_key), function(d) {
    rate <- mean(d$sig)
    ci <- wilson_ci(sum(d$sig), nrow(d))
    data.frame(
      scenario = d$scenario[[1]],
      N = d$N[[1]],
      gamma_U = d$gamma_U[[1]],
      gamma_T = d$gamma_T[[1]],
      R = d$R[[1]],
      cor_E1E2 = d$cor_E1E2[[1]],
      H0 = d$H0[[1]],
      estimand = if (isTRUE(d$H0[[1]])) "typeI" else "power",
      estimator = d$estimator[[1]],
      B = d$B[[1]],
      ci_method = d$ci_method[[1]],
      rejection_rate = rate,
      ci_lo = ci[["ci_lo"]],
      ci_hi = ci[["ci_hi"]],
      n_reps = nrow(d),
      mean_estimate = mean(d$estimate),
      sd_estimate = stats::sd(d$estimate),
      bca_fallbacks = sum(d$bca_fallback),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$scenario, out$estimator, out$B, out$ci_method), ]
}

bootstrap_sensitivity <- run_bootstrap_sensitivity()
write_csv(bootstrap_sensitivity, file.path("output", "tables", "revision_bootstrap_sensitivity.csv"))

# ---------------------------------------------------------------------------
# 3. Tukey g-and-h skew-heavy H1 sensitivity
# ---------------------------------------------------------------------------

tgh_raw <- function(n, g = 0.35, h = 0.10) {
  z <- stats::rnorm(n)
  if (abs(g) < 1e-12) {
    out <- z * exp(h * z^2 / 2)
  } else {
    out <- (exp(g * z) - 1) / g * exp(h * z^2 / 2)
  }
  out
}

tgh_ref <- local({
  set.seed(SEED_BASE + 300000L)
  raw <- tgh_raw(300000L)
  mu <- mean(raw)
  sig <- stats::sd(raw)
  z <- (raw - mu) / sig
  list(
    mu = mu,
    sig = sig,
    gamma3 = mean(z^3),
    gamma4 = mean(z^4) - 3
  )
})

tgh_component <- function(n) {
  (tgh_raw(n) - tgh_ref$mu) / tgh_ref$sig
}

generate_ws_tgh_dgp <- function(N, R, cor_E1E2, seed = NULL,
                                sigma_T = 1.0, lambda = 1.0) {
  if (!is.null(seed)) set.seed(seed)
  T_raw <- stats::rnorm(N) * sigma_T
  a1 <- 1.0
  a2 <- R
  c_sq <- cor_E1E2^2
  A_q <- c_sq
  B_q <- c_sq * (1 + R^2)
  C_q <- R^2 * (c_sq - 1)
  disc <- B_q^2 - 4 * A_q * C_q
  sigma_W_sq <- (-B_q + sqrt(disc)) / (2 * A_q)
  sigma_W <- sqrt(max(sigma_W_sq, 0))

  U_raw <- tgh_component(N)
  W1 <- stats::rnorm(N) * sigma_W
  W2 <- stats::rnorm(N) * sigma_W
  data.frame(
    x1 = lambda * T_raw + a1 * U_raw + W1,
    x2 = lambda * T_raw + a2 * U_raw + W2
  )
}

bootstrap_ci_percentile_only <- function(estimator_fn, x1, x2, B = 200L,
                                         alpha = 0.05, seed = SEED_BASE) {
  n <- length(x1)
  est <- estimator_fn(x1, x2)
  set.seed(seed)
  boot_vals <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    estimator_fn(x1[idx], x2[idx])
  }, numeric(1L))
  qs <- stats::quantile(boot_vals, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  c(estimate = est, ci_lo = qs[[1]], ci_hi = qs[[2]],
    sig = as.integer(qs[[1]] > 0 || qs[[2]] < 0))
}

run_heavytail_sensitivity <- function() {
  heavy_reps <- 80L
  conditions <- expand.grid(
    N = c(500L, 2000L),
    R = c(1.25, 2.0),
    stringsAsFactors = FALSE
  )
  conditions$cor_E1E2 <- ifelse(conditions$R == 1.25, 0.19, 0.38)
  conditions$scenario <- paste0("Tukey_g0.35_h0.10_R", conditions$R, "_N", conditions$N)

  rows <- list()
  for (ci in seq_len(nrow(conditions))) {
    cond <- conditions[ci, ]
    cat(sprintf("Heavy-tail sensitivity: %s (%d/%d)\n", cond$scenario, ci, nrow(conditions)))
    for (r in seq_len(heavy_reps)) {
      seed_r <- SEED_BASE + 400000L + ci * 10000L + r
      dat <- generate_ws_tgh_dgp(
        N = cond$N, R = cond$R, cor_E1E2 = cond$cor_E1E2, seed = seed_r
      )
      vals <- list(
        naive_delta_c3 = bootstrap_ci_percentile_only(delta_c3_naive, dat$x1, dat$x2,
                                                       B = 200L, seed = seed_r + 31L),
        pmm2_delta_c3 = bootstrap_ci_percentile_only(delta_c3_pmm2, dat$x1, dat$x2,
                                                      B = 200L, seed = seed_r + 37L)
      )
      for (est_name in names(vals)) {
        v <- vals[[est_name]]
        rows[[length(rows) + 1L]] <- data.frame(
          scenario = cond$scenario,
          rep = r,
          N = cond$N,
          R = cond$R,
          cor_E1E2 = cond$cor_E1E2,
          estimator = est_name,
          estimate = v[["estimate"]],
          sig = v[["sig"]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  raw <- do.call(rbind, rows)
  split_key <- interaction(raw$scenario, raw$estimator, drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(raw, split_key), function(d) {
    ci <- wilson_ci(sum(d$sig), nrow(d))
    data.frame(
      scenario = d$scenario[[1]],
      N = d$N[[1]],
      R = d$R[[1]],
      cor_E1E2 = d$cor_E1E2[[1]],
      U_family = "Tukey g-and-h",
      U_g = 0.35,
      U_h = 0.10,
      U_ref_skewness = tgh_ref$gamma3,
      U_ref_excess_kurtosis = tgh_ref$gamma4,
      estimator = d$estimator[[1]],
      power = mean(d$sig),
      ci_lo = ci[["ci_lo"]],
      ci_hi = ci[["ci_hi"]],
      n_reps = nrow(d),
      mean_estimate = mean(d$estimate),
      sd_estimate = stats::sd(d$estimate),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL

  # Add PMM2-to-naive signal and power ratios per condition.
  naive <- out[out$estimator == "naive_delta_c3", ]
  pmm2 <- out[out$estimator == "pmm2_delta_c3", ]
  key <- paste(naive$scenario)
  idx <- match(paste(pmm2$scenario), key)
  ratio <- data.frame(
    scenario = pmm2$scenario,
    estimator = "pmm2_vs_naive_ratio",
    signal_ratio = pmm2$mean_estimate / naive$mean_estimate[idx],
    power_diff_pp = 100 * (pmm2$power - naive$power[idx]),
    stringsAsFactors = FALSE
  )
  list(
    table = out[order(out$scenario, out$estimator), ],
    ratio = ratio[order(ratio$scenario), ]
  )
}

heavy_result <- run_heavytail_sensitivity()
heavy <- heavy_result$table
heavy_ratio <- heavy_result$ratio
write_csv(heavy, file.path("output", "tables", "revision_heavytail_sensitivity.csv"))
write_csv(heavy_ratio, file.path("output", "tables", "revision_heavytail_ratios.csv"))

# ---------------------------------------------------------------------------
# 4. Raw dCov sanity baseline
# ---------------------------------------------------------------------------

dcov_stat <- function(x, y) {
  A <- abs(outer(x, x, "-"))
  B <- abs(outer(y, y, "-"))
  A <- sweep(sweep(A, 1L, rowMeans(A), "-"), 2L, colMeans(A), "-") + mean(A)
  B <- sweep(sweep(B, 1L, rowMeans(B), "-"), 2L, colMeans(B), "-") + mean(B)
  mean(A * B)
}

dcov_perm_test <- function(x, y, B_perm = 99L, seed = SEED_BASE) {
  set.seed(seed)
  A <- abs(outer(x, x, "-"))
  Y <- abs(outer(y, y, "-"))
  A <- sweep(sweep(A, 1L, rowMeans(A), "-"), 2L, colMeans(A), "-") + mean(A)
  Y <- sweep(sweep(Y, 1L, rowMeans(Y), "-"), 2L, colMeans(Y), "-") + mean(Y)
  obs <- mean(A * Y)
  perm <- vapply(seq_len(B_perm), function(i) {
    idx <- sample.int(length(y))
    mean(A * Y[idx, idx])
  }, numeric(1L))
  pval <- (1 + sum(perm >= obs)) / (B_perm + 1)
  c(stat = obs, p_value = pval, reject = as.integer(pval < ALPHA))
}

run_dcov_sanity <- function() {
  dcov_reps <- 35L
  conditions <- data.frame(
    scenario = c("H0_common_true_score", "H1_boundary", "H1_strong"),
    N = c(300L, 300L, 300L),
    gamma_U = c(0.75, 2.25, 2.25),
    gamma_T = c(0, 0, 0),
    R = c(1.0, 1.25, 2.0),
    cor_E1E2 = c(0.0, 0.19, 0.38),
    H0_error_independence = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  rows <- list()
  for (ci in seq_len(nrow(conditions))) {
    cond <- conditions[ci, ]
    cat(sprintf("dCov sanity: %s (%d/%d)\n", cond$scenario, ci, nrow(conditions)))
    for (r in seq_len(dcov_reps)) {
      seed_r <- SEED_BASE + 500000L + ci * 10000L + r
      dat <- generate_ws_dgp(
        N = cond$N,
        gamma_U = cond$gamma_U, gamma_T = cond$gamma_T, gamma_W = 0,
        delta_U = 0, delta_T = 0, delta_W = 0,
        R = cond$R, cor_E1E2 = cond$cor_E1E2,
        seed = seed_r
      )
      test <- dcov_perm_test(dat$x1, dat$x2, B_perm = 99L, seed = seed_r + 41L)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = cond$scenario,
        rep = r,
        N = cond$N,
        R = cond$R,
        cor_E1E2 = cond$cor_E1E2,
        H0_error_independence = cond$H0_error_independence,
        stat = test[["stat"]],
        p_value = test[["p_value"]],
        reject = test[["reject"]],
        stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  split_key <- interaction(raw$scenario, drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(raw, split_key), function(d) {
    ci <- wilson_ci(sum(d$reject), nrow(d))
    data.frame(
      scenario = d$scenario[[1]],
      N = d$N[[1]],
      R = d$R[[1]],
      cor_E1E2 = d$cor_E1E2[[1]],
      H0_error_independence = d$H0_error_independence[[1]],
      raw_dcov_rejection_rate = mean(d$reject),
      ci_lo = ci[["ci_lo"]],
      ci_hi = ci[["ci_hi"]],
      n_reps = nrow(d),
      permutations = 99L,
      mean_p_value = mean(d$p_value),
      median_p_value = stats::median(d$p_value),
      interpretation = "Raw dCov tests observed-score independence, not W&S error-independence; common T creates dependence under H0.",
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$scenario), ]
}

dcov <- run_dcov_sanity()
write_csv(dcov, file.path("output", "tables", "revision_dcov_sanity.csv"))

# ---------------------------------------------------------------------------
# Session info
# ---------------------------------------------------------------------------

sink(file.path("output", "session_info", "revision_experiments_session_info.txt"))
cat("Revision experiment session info\n")
cat("================================\n")
cat(sprintf("Timestamp: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("Seed base: %d\n", SEED_BASE))
cat(sprintf("Tukey g-and-h U reference skewness: %.4f\n", tgh_ref$gamma3))
cat(sprintf("Tukey g-and-h U reference excess kurtosis: %.4f\n", tgh_ref$gamma4))
print(utils::sessionInfo())
sink()

cat("=== Revision experiments completed ===\n")

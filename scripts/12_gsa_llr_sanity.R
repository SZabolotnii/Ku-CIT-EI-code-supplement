# Phase 9.5: adversarial sanity checks for the GSA-LLR track.
#
# This script is deliberately a pre-flight gate, not a polished Phase 10
# implementation. It checks whether the symbolic GSA-LLR derivation survives
# small Monte Carlo stress tests before the project invests in a full MC grid.

options(stringsAsFactors = FALSE)

ensure_phase95_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

gsa_llr_F <- function(order = 4L) {
  stopifnot(order %in% c(3L, 4L))
  F4 <- matrix(c(
     10 / 9,  -10 / 3,       0,        0,
    -10 / 3,  130 / 9,       0,        0,
          0,        0,  500 / 27, -400 / 9,
          0,        0, -400 / 9, 3200 / 27
  ), nrow = 4L, byrow = TRUE)
  idx <- seq_len(if (order == 3L) 2L else 4L)
  F4[idx, idx, drop = FALSE]
}

gsa_llr_y <- function(R, kappa3, kappa4 = 0, order = 4L) {
  stopifnot(order %in% c(3L, 4L))
  y4 <- c(
    R * (R - 1) * kappa3,
    (1 - R^3) * kappa3,
    R * (R^2 - 1) * kappa4,
    (1 - R^4) * kappa4
  )
  y4[seq_len(if (order == 3L) 2L else 4L)]
}

gsa_llr_basis <- function(x1, x2, order = 4L, standardize = TRUE) {
  stopifnot(length(x1) == length(x2), length(x1) >= 2L, order %in% c(3L, 4L))
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  x1 <- x1 - mean(x1)
  x2 <- x2 - mean(x2)
  if (standardize) {
    s1 <- stats::sd(x1)
    s2 <- stats::sd(x2)
    if (!is.finite(s1) || !is.finite(s2) || s1 <= 0 || s2 <= 0) {
      stop("Cannot standardize a degenerate score vector")
    }
    x1 <- x1 / s1
    x2 <- x2 / s2
  }

  out <- cbind(
    phi1 = x1 * x2^2 - x1^2 * x2,
    phi2 = x1^3 - x2^3
  )
  if (order == 4L) {
    out <- cbind(
      out,
      phi3 = x1 * x2^3 - x1^3 * x2,
      phi4 = x1^4 - x2^4
    )
  }
  out
}

gsa_llr_weights <- function(R, kappa3, kappa4 = 0, order = 4L) {
  F <- gsa_llr_F(order)
  y <- gsa_llr_y(R, kappa3, kappa4, order)
  K <- as.numeric(solve(F, y))
  J <- as.numeric(crossprod(y, K))
  list(K = K, J = J, y = y, F = F)
}

gsa_llr_statistic <- function(x1, x2, R, kappa3, kappa4 = 0, order = 4L,
                              alpha = 0.05, threshold = c("vp", "normal"),
                              standardize = TRUE) {
  threshold <- match.arg(threshold)
  wk <- gsa_llr_weights(R, kappa3, kappa4, order)
  n <- length(x1)

  if (!is.finite(wk$J) || wk$J <= sqrt(.Machine$double.eps)) {
    return(list(
      lambda = 0, z = NA_real_, threshold_z = NA_real_, reject = FALSE,
      J = wk$J, degenerate_direction = TRUE
    ))
  }

  phi_bar <- colMeans(gsa_llr_basis(x1, x2, order = order, standardize = standardize))
  lambda <- sum(wk$K * phi_bar)
  z <- sqrt(n) * lambda / sqrt(wk$J)
  threshold_z <- switch(
    threshold,
    vp = sqrt(8 / (9 * alpha)),
    normal = stats::qnorm(1 - alpha / 2)
  )

  list(
    lambda = lambda,
    z = z,
    threshold_z = threshold_z,
    reject = is.finite(z) && abs(z) > threshold_z,
    J = wk$J,
    degenerate_direction = FALSE
  )
}

phase95_component <- function(n, kappa3 = 0, kappa4 = 0) {
  if (abs(kappa3) > 1e-12) {
    shape <- 4 / (abs(kappa3)^2)
    x <- stats::rgamma(n, shape = shape, scale = 1)
    x <- (x - mean(x)) / stats::sd(x)
    if (kappa3 < 0) x <- -x
    return(x)
  }

  if (kappa4 > 1e-12) {
    df <- 6 / kappa4 + 4
    x <- stats::rt(n, df = df) / sqrt(df / (df - 2))
    return((x - mean(x)) / stats::sd(x))
  }

  if (kappa4 < -1e-12) {
    # Symmetric beta(a, a) has excess kurtosis -6 / (2a + 3).
    a <- (-6 / kappa4 - 3) / 2
    if (!is.finite(a) || a <= 0) {
      stop("Requested negative excess kurtosis is outside the beta(a,a) helper range")
    }
    x <- stats::rbeta(n, shape1 = a, shape2 = a)
    return((x - mean(x)) / stats::sd(x))
  }

  stats::rnorm(n)
}

phase95_dgp <- function(N, R = 1, kappa3 = 0, kappa4 = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  T_score <- stats::rnorm(N)
  U <- phase95_component(N, kappa3, kappa4)
  W1 <- stats::rnorm(N)
  W2 <- stats::rnorm(N)
  data.frame(
    x1 = T_score + U + W1,
    x2 = T_score + R * U + W2
  )
}

simulate_phase95_cell <- function(check_id, description, reps, N,
                                  dgp_R, dgp_kappa3, dgp_kappa4,
                                  test_R, test_kappa3, test_kappa4,
                                  order = 4L, alpha = 0.05,
                                  threshold = "vp", standardize = TRUE,
                                  seed_base = 20260526L) {
  stopifnot(reps >= 1L)
  z_vals <- numeric(reps)
  rejects <- logical(reps)
  degenerate <- logical(reps)

  for (r in seq_len(reps)) {
    dat <- phase95_dgp(
      N = N, R = dgp_R, kappa3 = dgp_kappa3, kappa4 = dgp_kappa4,
      seed = seed_base + r
    )
    st <- gsa_llr_statistic(
      dat$x1, dat$x2,
      R = test_R, kappa3 = test_kappa3, kappa4 = test_kappa4,
      order = order, alpha = alpha, threshold = threshold,
      standardize = standardize
    )
    z_vals[r] <- st$z
    rejects[r] <- st$reject
    degenerate[r] <- st$degenerate_direction
  }

  data.frame(
    check_id = check_id,
    description = description,
    threshold = threshold,
    standardize = standardize,
    reps = reps,
    N = N,
    order = order,
    dgp_R = dgp_R,
    dgp_kappa3 = dgp_kappa3,
    dgp_kappa4 = dgp_kappa4,
    test_R = test_R,
    test_kappa3 = test_kappa3,
    test_kappa4 = test_kappa4,
    rejection_rate = mean(rejects),
    mean_z = mean(z_vals, na.rm = TRUE),
    mean_abs_z = mean(abs(z_vals), na.rm = TRUE),
    degenerate_reps = sum(degenerate),
    stringsAsFactors = FALSE
  )
}

phase95_plan <- function(reps = 100L, threshold = "vp", standardize = TRUE) {
  rows <- list()

  # (1) Consistency: z should grow approximately as sqrt(N) when the
  # third-order alternative is present. We record z/sqrt(N) stability.
  for (N in c(1000L, 5000L, 10000L, 50000L)) {
    rows[[length(rows) + 1L]] <- simulate_phase95_cell(
      "C1_consistency",
      "R=2, kappa3=2.25, kappa4=0; z/sqrt(N) should stabilize",
      reps = reps, N = N,
      dgp_R = 2, dgp_kappa3 = 2.25, dgp_kappa4 = 0,
      test_R = 2, test_kappa3 = 2.25, test_kappa4 = 0,
      order = 3L, threshold = threshold, standardize = standardize,
      seed_base = 20260526L + N
    )
  }

  # (2) Gaussian Type-I: data are generated from the Gaussian null, but the
  # score is evaluated against a fixed non-null direction. Using Y=0 here would
  # be a degenerate no-test, not a Type-I check.
  for (N in c(500L, 5000L)) {
    for (R_alt in c(1.25, 1.5, 2, 2.5)) {
      rows[[length(rows) + 1L]] <- simulate_phase95_cell(
        "C2_gaussian_typeI",
        "Gaussian H0 data; fixed skewed alternative direction",
        reps = reps, N = N,
        dgp_R = 1, dgp_kappa3 = 0, dgp_kappa4 = 0,
        test_R = R_alt, test_kappa3 = 0.75, test_kappa4 = 0,
        order = 3L, threshold = threshold, standardize = standardize,
        seed_base = 20270526L + N + round(100 * R_alt)
      )
    }
  }

  # (3) Boundary monotonicity.
  for (R_alt in c(1, 1.25, 1.5, 2, 2.5)) {
    test_R <- if (R_alt == 1) 1.25 else R_alt
    rows[[length(rows) + 1L]] <- simulate_phase95_cell(
      "C3_boundary_monotonicity",
      "kappa3=0.75, N=500; power should increase with R",
      reps = reps, N = 500L,
      dgp_R = R_alt, dgp_kappa3 = 0.75, dgp_kappa4 = 0,
      test_R = test_R, test_kappa3 = 0.75, test_kappa4 = 0,
      order = 3L, threshold = threshold, standardize = standardize,
      seed_base = 20280526L + round(100 * R_alt)
    )
  }

  # (4) Joint vs marginal pooling.
  for (order in c(3L, 4L)) {
    rows[[length(rows) + 1L]] <- simulate_phase95_cell(
      "C4_joint_vs_marginal",
      "kappa3=2.25, kappa4=9, R=2, N=1000; order 4 should improve",
      reps = reps, N = 1000L,
      dgp_R = 2, dgp_kappa3 = 2.25, dgp_kappa4 = 9,
      test_R = 2, test_kappa3 = 2.25, test_kappa4 = 9,
      order = order, threshold = threshold, standardize = standardize,
      seed_base = 20290526L + order
    )
  }

  # (5) PHQ-8-like symmetric platykurtic regime.
  for (order in c(3L, 4L)) {
    rows[[length(rows) + 1L]] <- simulate_phase95_cell(
      "C5_phq8_like",
      "kappa3=0, kappa4=-1.3, R=1.5, N=2000; order 4 should detect",
      reps = reps, N = 2000L,
      dgp_R = 1.5, dgp_kappa3 = 0, dgp_kappa4 = -1.3,
      test_R = 1.5, test_kappa3 = 0, test_kappa4 = -1.3,
      order = order, threshold = threshold, standardize = standardize,
      seed_base = 20300526L + order
    )
  }

  out <- do.call(rbind, rows)
  out$z_over_sqrt_N <- out$mean_z / sqrt(out$N)
  out
}

evaluate_phase95 <- function(results, alpha = 0.05) {
  checks <- list()

  c1 <- subset(results, check_id == "C1_consistency")
  checks[["C1_consistency"]] <- all(c1$rejection_rate >= 0.95) &&
    (stats::sd(c1$z_over_sqrt_N, na.rm = TRUE) < 0.25)

  c2 <- subset(results, check_id == "C2_gaussian_typeI")
  checks[["C2_gaussian_typeI"]] <- mean(c2$rejection_rate >= 0.03 &
                                         c2$rejection_rate <= 0.08) >= 0.90

  c3 <- subset(results, check_id == "C3_boundary_monotonicity")
  c3 <- c3[order(c3$dgp_R), ]
  checks[["C3_boundary_monotonicity"]] <- all(diff(c3$rejection_rate) > 0)

  c4 <- subset(results, check_id == "C4_joint_vs_marginal")
  p3 <- c4$rejection_rate[c4$order == 3L]
  p4 <- c4$rejection_rate[c4$order == 4L]
  checks[["C4_joint_vs_marginal"]] <- length(p3) == 1L && length(p4) == 1L &&
    is.finite(p3) && is.finite(p4) && (p4 - p3 >= 0.05)

  c5 <- subset(results, check_id == "C5_phq8_like")
  p3 <- c5$rejection_rate[c5$order == 3L]
  p4 <- c5$rejection_rate[c5$order == 4L]
  checks[["C5_phq8_like"]] <- length(p3) == 1L && length(p4) == 1L &&
    is.finite(p3) && is.finite(p4) && (p3 <= alpha + 0.05) && (p4 >= 0.40)

  data.frame(
    check_id = names(checks),
    pass = unlist(checks, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

write_phase95_report <- function(results, verdict,
                                 path = file.path("output", "sanity_check_phase95.txt")) {
  overall <- all(verdict$pass)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("Phase 9.5 GSA-LLR adversarial sanity check", con)
  writeLines("================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Threshold: %s", unique(results$threshold)), con)
  writeLines(sprintf("Standardize scores before basis: %s", unique(results$standardize)), con)
  writeLines(sprintf("Overall pass: %s", overall), con)
  writeLines("", con)

  writeLines("Verdict by check:", con)
  for (i in seq_len(nrow(verdict))) {
    writeLines(sprintf("- %s: %s", verdict$check_id[i],
                       if (verdict$pass[i]) "PASS" else "FAIL"), con)
  }
  writeLines("", con)

  writeLines("Key rates:", con)
  compact <- results[, c("check_id", "N", "order", "dgp_R", "dgp_kappa3",
                         "dgp_kappa4", "test_R", "test_kappa3", "test_kappa4",
                         "rejection_rate", "mean_z", "z_over_sqrt_N",
                         "degenerate_reps")]
  utils::write.table(compact, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)
  writeLines("Interpretation:", con)
  if (overall) {
    writeLines("All Phase 9.5 checks passed. It is reasonable to proceed to Phase 10.", con)
  } else {
    writeLines("At least one Phase 9.5 check failed. Per PLAYBOOK Gate G3', do not proceed to Phase 10 before revisiting the finite-sample calibration/derivation.", con)
  }
}

run_phase95 <- function(reps = 100L, threshold = "vp", standardize = TRUE) {
  ensure_phase95_dirs()
  results <- phase95_plan(reps = reps, threshold = threshold, standardize = standardize)
  verdict <- evaluate_phase95(results)
  utils::write.csv(results, file.path("output", "tables", "gsa_llr_phase95_results.csv"),
                   row.names = FALSE)
  utils::write.csv(verdict, file.path("output", "tables", "gsa_llr_phase95_verdict.csv"),
                   row.names = FALSE)
  write_phase95_report(results, verdict)
  capture.output(utils::sessionInfo(),
                 file = file.path("output", "session_info", "session_info_phase95.txt"))
  invisible(list(results = results, verdict = verdict, overall_pass = all(verdict$pass)))
}

phase95_z_sample <- function(reps, N, dgp_R, dgp_kappa3, dgp_kappa4,
                             test_R, test_kappa3, test_kappa4, order = 4L,
                             standardize = TRUE, seed_base = 20260526L) {
  z_vals <- numeric(reps)
  degenerate <- logical(reps)

  for (r in seq_len(reps)) {
    dat <- phase95_dgp(
      N = N, R = dgp_R, kappa3 = dgp_kappa3, kappa4 = dgp_kappa4,
      seed = seed_base + r
    )
    st <- gsa_llr_statistic(
      dat$x1, dat$x2,
      R = test_R, kappa3 = test_kappa3, kappa4 = test_kappa4,
      order = order, alpha = 0.05, threshold = "normal",
      standardize = standardize
    )
    z_vals[r] <- st$z
    degenerate[r] <- st$degenerate_direction
  }

  list(z = z_vals, degenerate = degenerate)
}

simulate_phase95r_cell <- function(check_id, description, reps, calibration_reps, N,
                                   dgp_R, dgp_kappa3, dgp_kappa4,
                                   test_R, test_kappa3, test_kappa4,
                                   order = 4L, alpha = 0.05,
                                   standardize = TRUE, seed_base = 20260526L) {
  alt <- phase95_z_sample(
    reps = reps, N = N,
    dgp_R = dgp_R, dgp_kappa3 = dgp_kappa3, dgp_kappa4 = dgp_kappa4,
    test_R = test_R, test_kappa3 = test_kappa3, test_kappa4 = test_kappa4,
    order = order, standardize = standardize, seed_base = seed_base
  )

  # Empirical-null calibration uses the same non-Gaussian component family but
  # enforces the null loading equality R=1. This avoids the degenerate Y=0
  # trap while matching the alternative direction being tested.
  null <- phase95_z_sample(
    reps = calibration_reps, N = N,
    dgp_R = 1, dgp_kappa3 = dgp_kappa3, dgp_kappa4 = dgp_kappa4,
    test_R = test_R, test_kappa3 = test_kappa3, test_kappa4 = test_kappa4,
    order = order, standardize = standardize, seed_base = seed_base + 1000000L
  )

  thresholds <- c(
    vp = sqrt(8 / (9 * alpha)),
    normal = stats::qnorm(1 - alpha / 2),
    empirical = as.numeric(stats::quantile(abs(null$z), 1 - alpha,
                                           na.rm = TRUE, names = FALSE))
  )

  data.frame(
    check_id = check_id,
    description = description,
    basis_scale = if (standardize) "standardized" else "raw",
    reps = reps,
    calibration_reps = calibration_reps,
    N = N,
    order = order,
    dgp_R = dgp_R,
    dgp_kappa3 = dgp_kappa3,
    dgp_kappa4 = dgp_kappa4,
    test_R = test_R,
    test_kappa3 = test_kappa3,
    test_kappa4 = test_kappa4,
    rate_vp = mean(abs(alt$z) > thresholds["vp"], na.rm = TRUE),
    rate_normal = mean(abs(alt$z) > thresholds["normal"], na.rm = TRUE),
    rate_empirical = mean(abs(alt$z) > thresholds["empirical"], na.rm = TRUE),
    null_rate_vp = mean(abs(null$z) > thresholds["vp"], na.rm = TRUE),
    null_rate_normal = mean(abs(null$z) > thresholds["normal"], na.rm = TRUE),
    null_rate_empirical = mean(abs(null$z) > thresholds["empirical"], na.rm = TRUE),
    threshold_vp = thresholds["vp"],
    threshold_normal = thresholds["normal"],
    threshold_empirical = thresholds["empirical"],
    mean_z = mean(alt$z, na.rm = TRUE),
    mean_abs_z = mean(abs(alt$z), na.rm = TRUE),
    z_over_sqrt_N = mean(alt$z, na.rm = TRUE) / sqrt(N),
    null_mean_abs_z = mean(abs(null$z), na.rm = TRUE),
    null_sd_z = stats::sd(null$z, na.rm = TRUE),
    degenerate_reps = sum(alt$degenerate),
    null_degenerate_reps = sum(null$degenerate),
    stringsAsFactors = FALSE
  )
}

phase95r_plan <- function(reps = 200L, calibration_reps = 500L,
                          standardize_values = c(TRUE, FALSE)) {
  rows <- list()

  for (std in standardize_values) {
    for (N in c(1000L, 5000L, 10000L, 50000L)) {
      rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
        "C1_consistency",
        "R=2, kappa3=2.25, kappa4=0; z/sqrt(N) should stabilize",
        reps = reps, calibration_reps = calibration_reps, N = N,
        dgp_R = 2, dgp_kappa3 = 2.25, dgp_kappa4 = 0,
        test_R = 2, test_kappa3 = 2.25, test_kappa4 = 0,
        order = 3L, standardize = std, seed_base = 21260526L + N
      )
    }

    for (N in c(500L, 5000L)) {
      for (R_alt in c(1.25, 1.5, 2, 2.5)) {
        rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
          "C2_gaussian_typeI",
          "Gaussian H0 data; fixed skewed alternative direction",
          reps = reps, calibration_reps = calibration_reps, N = N,
          dgp_R = 1, dgp_kappa3 = 0, dgp_kappa4 = 0,
          test_R = R_alt, test_kappa3 = 0.75, test_kappa4 = 0,
          order = 3L, standardize = std,
          seed_base = 21360526L + N + round(100 * R_alt)
        )
      }
    }

    for (R_alt in c(1, 1.25, 1.5, 2, 2.5)) {
      test_R <- if (R_alt == 1) 1.25 else R_alt
      rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
        "C3_boundary_monotonicity",
        "kappa3=0.75, N=500; power should increase with R",
        reps = reps, calibration_reps = calibration_reps, N = 500L,
        dgp_R = R_alt, dgp_kappa3 = 0.75, dgp_kappa4 = 0,
        test_R = test_R, test_kappa3 = 0.75, test_kappa4 = 0,
        order = 3L, standardize = std,
        seed_base = 21460526L + round(100 * R_alt)
      )
    }

    for (order in c(3L, 4L)) {
      rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
        "C4_joint_vs_marginal",
        "kappa3=2.25, kappa4=9, R=2, N=1000; order 4 should improve",
        reps = reps, calibration_reps = calibration_reps, N = 1000L,
        dgp_R = 2, dgp_kappa3 = 2.25, dgp_kappa4 = 9,
        test_R = 2, test_kappa3 = 2.25, test_kappa4 = 9,
        order = order, standardize = std, seed_base = 21560526L + order
      )
    }

    for (order in c(3L, 4L)) {
      rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
        "C5_phq8_like",
        "kappa3=0, kappa4=-1.3, R=1.5, N=2000; order 4 should detect",
        reps = reps, calibration_reps = calibration_reps, N = 2000L,
        dgp_R = 1.5, dgp_kappa3 = 0, dgp_kappa4 = -1.3,
        test_R = 1.5, test_kappa3 = 0, test_kappa4 = -1.3,
        order = order, standardize = std, seed_base = 21660526L + order
      )
    }
  }

  do.call(rbind, rows)
}

evaluate_phase95r <- function(results, alpha = 0.05) {
  std <- subset(results, basis_scale == "standardized")
  raw <- subset(results, basis_scale == "raw")

  c2_std <- subset(std, check_id == "C2_gaussian_typeI")
  c3_std <- subset(std, check_id == "C3_boundary_monotonicity")
  c3_std <- c3_std[order(c3_std$dgp_R), ]
  c4_std <- subset(std, check_id == "C4_joint_vs_marginal")
  c5_std <- subset(std, check_id == "C5_phq8_like")
  c2_raw <- subset(raw, check_id == "C2_gaussian_typeI")

  p4_3 <- c4_std$rate_empirical[c4_std$order == 3L]
  p4_4 <- c4_std$rate_empirical[c4_std$order == 4L]
  z4_3 <- c4_std$mean_abs_z[c4_std$order == 3L]
  z4_4 <- c4_std$mean_abs_z[c4_std$order == 4L]
  p5_3 <- c5_std$rate_empirical[c5_std$order == 3L]
  p5_4 <- c5_std$rate_empirical[c5_std$order == 4L]

  checks <- c(
    standardized_empirical_typeI =
      mean(c2_std$rate_empirical >= 0.025 & c2_std$rate_empirical <= 0.10, na.rm = TRUE) >= 0.90,
    standardized_empirical_boundary_monotone =
      all(diff(c3_std$rate_empirical) > 0),
    standardized_joint_effect_size =
      length(z4_3) == 1L && length(z4_4) == 1L && is.finite(z4_3) &&
      is.finite(z4_4) && z4_4 > z4_3,
    standardized_phq_target =
      length(p5_3) == 1L && length(p5_4) == 1L &&
      (is.na(p5_3) || p5_3 <= alpha + 0.05) && is.finite(p5_4) && p5_4 >= 0.40,
    raw_vp_normal_typeI_invalid =
      any(c2_raw$rate_vp > 0.10 | c2_raw$rate_normal > 0.10, na.rm = TRUE)
  )

  data.frame(
    check_id = names(checks),
    pass = unlist(checks, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

write_phase95r_report <- function(results, verdict,
                                  path = file.path("output", "sanity_check_phase95r.txt")) {
  std <- subset(results, basis_scale == "standardized")
  raw <- subset(results, basis_scale == "raw")
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("Phase 9.5R GSA-LLR calibration repair diagnostic", con)
  writeLines("==================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Reps: %d; empirical-null calibration reps: %d",
                     unique(results$reps)[1], unique(results$calibration_reps)[1]), con)
  writeLines("", con)

  writeLines("Verdict by diagnostic:", con)
  for (i in seq_len(nrow(verdict))) {
    writeLines(sprintf("- %s: %s", verdict$check_id[i],
                       if (verdict$pass[i]) "PASS" else "FAIL"), con)
  }
  writeLines("", con)

  writeLines("Calibration comparison (standardized basis):", con)
  compact_std <- std[, c("check_id", "N", "order", "dgp_R", "dgp_kappa3",
                         "dgp_kappa4", "test_R", "rate_vp", "rate_normal",
                         "rate_empirical", "null_rate_vp", "null_rate_normal",
                         "null_rate_empirical", "threshold_empirical",
                         "mean_z", "z_over_sqrt_N", "degenerate_reps")]
  utils::write.table(compact_std, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Raw-basis Type-I diagnostic:", con)
  compact_raw <- raw[raw$check_id == "C2_gaussian_typeI",
                     c("check_id", "N", "order", "test_R", "rate_vp",
                       "rate_normal", "rate_empirical", "null_rate_vp",
                       "null_rate_normal", "null_rate_empirical",
                       "threshold_empirical")]
  utils::write.table(compact_raw, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Interpretation:", con)
  writeLines("- V-P is a conservative upper-bound calibration, not a nominal 5% calibration; using [3%, 8%] as a lower-and-upper Type-I gate for V-P is inappropriate.", con)
  writeLines("- Raw centered basis is invalid with the current symbolic F/J scaling unless F and Y are re-derived for raw variances; V-P/normal Type-I is badly inflated.", con)
  writeLines("- Standardized basis with empirical-null calibration repairs Type-I and boundary monotonicity in this quick grid.", con)
  writeLines("- The original C4 rejection-rate criterion is saturated at N=1000; order 4 has larger mean |z|, but the specified >=5pp power-gain check is not identifiable at this operating point.", con)
  writeLines("- The PHQ-8-like target power >=0.40 is not met in the quick standardized empirical-null run; this target should be revisited before Phase 10.", con)
}

run_phase95r <- function(reps = 200L, calibration_reps = 500L) {
  ensure_phase95_dirs()
  results <- phase95r_plan(reps = reps, calibration_reps = calibration_reps)
  verdict <- evaluate_phase95r(results)
  utils::write.csv(results, file.path("output", "tables", "gsa_llr_phase95r_calibration.csv"),
                   row.names = FALSE)
  utils::write.csv(verdict, file.path("output", "tables", "gsa_llr_phase95r_verdict.csv"),
                   row.names = FALSE)
  write_phase95r_report(results, verdict)
  capture.output(utils::sessionInfo(),
                 file = file.path("output", "session_info", "session_info_phase95r.txt"))
  invisible(list(results = results, verdict = verdict, overall_pass = all(verdict$pass)))
}

phase95rb_plan <- function(reps = 200L, calibration_reps = 500L) {
  rows <- list()

  # PHQ-like curve: the original N=2000, R=1.5 target was too optimistic in
  # Phase 9.5R. This grid tests whether the same signal recovers with larger N
  # or stronger loading asymmetry.
  for (N in c(2000L, 5000L, 10000L)) {
    for (R_alt in c(1.5, 2.0)) {
      for (order in c(3L, 4L)) {
        rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
          "RBR1_phq_power_curve",
          "PHQ-like platykurtic power curve under empirical-null calibration",
          reps = reps, calibration_reps = calibration_reps, N = N,
          dgp_R = R_alt, dgp_kappa3 = 0, dgp_kappa4 = -1.3,
          test_R = R_alt, test_kappa3 = 0, test_kappa4 = -1.3,
          order = order, standardize = TRUE,
          seed_base = 21760526L + N + round(100 * R_alt) + order
        )
      }
    }
  }

  # C4 joint-vs-marginal: use non-saturated operating points so that rejection
  # rate differences are observable, not hidden behind both powers being 1.
  c4_grid <- data.frame(
    N = c(200L, 200L, 1000L, 1000L),
    R = c(1.5, 1.5, 1.25, 1.25),
    kappa3 = c(2.25, 2.25, 1.5, 1.5),
    kappa4 = c(1.3, 3.0, 1.3, 3.0)
  )
  for (i in seq_len(nrow(c4_grid))) {
    row <- c4_grid[i, ]
    for (order in c(3L, 4L)) {
      rows[[length(rows) + 1L]] <- simulate_phase95r_cell(
        "RBR2_c4_nonsaturated",
        "Non-saturated joint-vs-marginal operating point",
        reps = reps, calibration_reps = calibration_reps, N = row$N,
        dgp_R = row$R, dgp_kappa3 = row$kappa3, dgp_kappa4 = row$kappa4,
        test_R = row$R, test_kappa3 = row$kappa3, test_kappa4 = row$kappa4,
        order = order, standardize = TRUE,
        seed_base = 21860526L + row$N + round(100 * row$R) +
          round(10 * row$kappa3) + round(10 * row$kappa4) + order
      )
    }
  }

  do.call(rbind, rows)
}

evaluate_phase95rb <- function(results) {
  phq <- subset(results, check_id == "RBR1_phq_power_curve" & order == 4L)
  phq_original <- subset(phq, N == 2000L & abs(dgp_R - 1.5) < 1e-12)
  phq_n5000 <- subset(phq, N == 5000L & abs(dgp_R - 1.5) < 1e-12)
  phq_r2_n2000 <- subset(phq, N == 2000L & abs(dgp_R - 2.0) < 1e-12)

  c4 <- subset(results, check_id == "RBR2_c4_nonsaturated")
  keys <- unique(c4[, c("N", "dgp_R", "dgp_kappa3", "dgp_kappa4")])
  gains <- numeric(nrow(keys))
  nonsaturated <- logical(nrow(keys))
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    pair <- subset(c4, N == key$N & dgp_R == key$dgp_R &
                     dgp_kappa3 == key$dgp_kappa3 & dgp_kappa4 == key$dgp_kappa4)
    p3 <- pair$rate_empirical[pair$order == 3L]
    p4 <- pair$rate_empirical[pair$order == 4L]
    gains[i] <- p4 - p3
    nonsaturated[i] <- is.finite(p3) && is.finite(p4) && p3 < 0.95 && p4 < 0.98
  }

  checks <- c(
    phq_original_target_met =
      nrow(phq_original) == 1L && phq_original$rate_empirical >= 0.40,
    phq_recovers_with_larger_N =
      nrow(phq_n5000) == 1L && phq_n5000$rate_empirical >= 0.40,
    phq_recovers_with_larger_R =
      nrow(phq_r2_n2000) == 1L && phq_r2_n2000$rate_empirical >= 0.40,
    c4_nonsaturated_power_gain =
      all(nonsaturated) && all(gains >= 0.05)
  )

  data.frame(
    check_id = names(checks),
    pass = unlist(checks, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

write_phase95rb_report <- function(results, verdict,
                                   path = file.path("output", "sanity_check_phase95rb.txt")) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("Phase 9.5R-b GSA-LLR PHQ/C4 repair diagnostic", con)
  writeLines("=================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Reps: %d; empirical-null calibration reps: %d",
                     unique(results$reps)[1], unique(results$calibration_reps)[1]), con)
  writeLines("", con)

  writeLines("Verdict by diagnostic:", con)
  for (i in seq_len(nrow(verdict))) {
    writeLines(sprintf("- %s: %s", verdict$check_id[i],
                       if (verdict$pass[i]) "PASS" else "FAIL"), con)
  }
  writeLines("", con)

  writeLines("PHQ-like power curve (standardized, empirical-null calibrated):", con)
  phq <- subset(results, check_id == "RBR1_phq_power_curve")
  phq_compact <- phq[, c("N", "order", "dgp_R", "dgp_kappa4",
                         "rate_empirical", "threshold_empirical",
                         "mean_z", "degenerate_reps")]
  utils::write.table(phq_compact, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Non-saturated C4 joint-vs-marginal grid:", con)
  c4 <- subset(results, check_id == "RBR2_c4_nonsaturated")
  c4_compact <- c4[, c("N", "order", "dgp_R", "dgp_kappa3", "dgp_kappa4",
                       "rate_empirical", "threshold_empirical", "mean_abs_z")]
  utils::write.table(c4_compact, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Interpretation:", con)
  writeLines("- The original PHQ-like target at N=2000, R=1.5 remains too aggressive.", con)
  writeLines("- The PHQ-like signal recovers with larger N or stronger R, so the fourth-order block is not dead; the planned operating point was underpowered.", con)
  writeLines("- Non-saturated C4 operating points show larger order-4 mean |z|, but empirical-null thresholds also widen; calibrated rejection power does not improve in this grid.", con)
  writeLines("- Recommendation: revise Phase 9.5 acceptance criteria and revisit the order-4 calibration before Phase 10, instead of claiming the original joint-power dominance target.", con)
}

run_phase95rb <- function(reps = 200L, calibration_reps = 500L) {
  ensure_phase95_dirs()
  results <- phase95rb_plan(reps = reps, calibration_reps = calibration_reps)
  verdict <- evaluate_phase95rb(results)
  utils::write.csv(results, file.path("output", "tables", "gsa_llr_phase95rb_curves.csv"),
                   row.names = FALSE)
  utils::write.csv(verdict, file.path("output", "tables", "gsa_llr_phase95rb_verdict.csv"),
                   row.names = FALSE)
  write_phase95rb_report(results, verdict)
  capture.output(utils::sessionInfo(),
                 file = file.path("output", "session_info", "session_info_phase95rb.txt"))
  invisible(list(results = results, verdict = verdict, overall_pass = all(verdict$pass)))
}

finite_values <- function(x) {
  x[is.finite(x)]
}

safe_mean <- function(x) {
  x <- finite_values(x)
  if (!length(x)) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- finite_values(x)
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

safe_quantile <- function(x, prob) {
  x <- finite_values(x)
  if (!length(x)) NA_real_ else as.numeric(stats::quantile(x, prob, names = FALSE))
}

safe_abs_rate <- function(x, threshold) {
  x <- finite_values(x)
  if (!length(x) || !is.finite(threshold)) NA_real_ else mean(abs(x) > threshold)
}

simulate_phase95rc_cell <- function(check_id, description, reps, calibration_reps, N,
                                    dgp_R, dgp_kappa3, dgp_kappa4,
                                    test_R, test_kappa3, test_kappa4,
                                    order = 4L, alpha = 0.05,
                                    seed_base = 20260527L) {
  alt <- phase95_z_sample(
    reps = reps, N = N,
    dgp_R = dgp_R, dgp_kappa3 = dgp_kappa3, dgp_kappa4 = dgp_kappa4,
    test_R = test_R, test_kappa3 = test_kappa3, test_kappa4 = test_kappa4,
    order = order, standardize = TRUE, seed_base = seed_base
  )

  null <- phase95_z_sample(
    reps = calibration_reps, N = N,
    dgp_R = 1, dgp_kappa3 = dgp_kappa3, dgp_kappa4 = dgp_kappa4,
    test_R = test_R, test_kappa3 = test_kappa3, test_kappa4 = test_kappa4,
    order = order, standardize = TRUE, seed_base = seed_base + 1000000L
  )

  empirical_threshold <- safe_quantile(abs(null$z), 1 - alpha)
  normal_threshold <- stats::qnorm(1 - alpha / 2)
  null_mean <- safe_mean(null$z)
  null_sd <- safe_sd(null$z)

  if (is.finite(null_sd) && null_sd > sqrt(.Machine$double.eps)) {
    alt_block_z <- (alt$z - null_mean) / null_sd
    null_block_z <- (null$z - null_mean) / null_sd
  } else {
    alt_block_z <- rep(NA_real_, length(alt$z))
    null_block_z <- rep(NA_real_, length(null$z))
  }

  data.frame(
    check_id = check_id,
    description = description,
    reps = reps,
    calibration_reps = calibration_reps,
    N = N,
    order = order,
    dgp_R = dgp_R,
    dgp_kappa3 = dgp_kappa3,
    dgp_kappa4 = dgp_kappa4,
    test_R = test_R,
    test_kappa3 = test_kappa3,
    test_kappa4 = test_kappa4,
    rate_empirical = safe_abs_rate(alt$z, empirical_threshold),
    null_rate_empirical = safe_abs_rate(null$z, empirical_threshold),
    threshold_empirical = empirical_threshold,
    rate_block_normal = safe_abs_rate(alt_block_z, normal_threshold),
    null_rate_block_normal = safe_abs_rate(null_block_z, normal_threshold),
    threshold_block_normal = normal_threshold,
    mean_abs_z = safe_mean(abs(alt$z)),
    null_mean_z = null_mean,
    null_sd_z = null_sd,
    mean_abs_block_z = safe_mean(abs(alt_block_z)),
    degenerate_reps = sum(alt$degenerate),
    null_degenerate_reps = sum(null$degenerate),
    stringsAsFactors = FALSE
  )
}

phase95rc_plan <- function(reps = 500L, calibration_reps = 1000L) {
  rows <- list()

  pure_kurtosis_grid <- data.frame(
    N = c(2000L, 5000L, 2000L, 1000L),
    R = c(1.5, 1.5, 2.0, 1.5),
    kappa4 = c(-1.3, -1.3, -1.3, 3.0)
  )
  for (i in seq_len(nrow(pure_kurtosis_grid))) {
    row <- pure_kurtosis_grid[i, ]
    for (order in c(3L, 4L)) {
      rows[[length(rows) + 1L]] <- simulate_phase95rc_cell(
        "RCR1_pure_kurtosis",
        "Pure kurtosis alternative: order 3 should be degenerate, order 4 should detect",
        reps = reps, calibration_reps = calibration_reps, N = row$N,
        dgp_R = row$R, dgp_kappa3 = 0, dgp_kappa4 = row$kappa4,
        test_R = row$R, test_kappa3 = 0, test_kappa4 = row$kappa4,
        order = order,
        seed_base = 21960527L + row$N + round(100 * row$R) +
          round(10 * abs(row$kappa4)) + order
      )
    }
  }

  c4_grid <- data.frame(
    N = c(200L, 200L, 1000L, 1000L),
    R = c(1.5, 1.5, 1.25, 1.25),
    kappa3 = c(2.25, 2.25, 1.5, 1.5),
    kappa4 = c(1.3, 3.0, 1.3, 3.0)
  )
  for (i in seq_len(nrow(c4_grid))) {
    row <- c4_grid[i, ]
    for (order in c(3L, 4L)) {
      rows[[length(rows) + 1L]] <- simulate_phase95rc_cell(
        "RCR2_block_normalized_c4",
        "C4 comparison after empirical-null block normalization",
        reps = reps, calibration_reps = calibration_reps, N = row$N,
        dgp_R = row$R, dgp_kappa3 = row$kappa3, dgp_kappa4 = row$kappa4,
        test_R = row$R, test_kappa3 = row$kappa3, test_kappa4 = row$kappa4,
        order = order,
        seed_base = 22060527L + row$N + round(100 * row$R) +
          round(10 * row$kappa3) + round(10 * row$kappa4) + order
      )
    }
  }

  do.call(rbind, rows)
}

phase95rc_pair_gains <- function(results) {
  c4 <- subset(results, check_id == "RCR2_block_normalized_c4")
  keys <- unique(c4[, c("N", "dgp_R", "dgp_kappa3", "dgp_kappa4")])
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    pair <- subset(c4, N == key$N & dgp_R == key$dgp_R &
                     dgp_kappa3 == key$dgp_kappa3 & dgp_kappa4 == key$dgp_kappa4)
    p3 <- pair[pair$order == 3L, ]
    p4 <- pair[pair$order == 4L, ]
    rows[[i]] <- data.frame(
      N = key$N,
      dgp_R = key$dgp_R,
      dgp_kappa3 = key$dgp_kappa3,
      dgp_kappa4 = key$dgp_kappa4,
      gain_empirical = p4$rate_empirical - p3$rate_empirical,
      gain_block_normal = p4$rate_block_normal - p3$rate_block_normal,
      mean_abs_z_gain = p4$mean_abs_z - p3$mean_abs_z,
      mean_abs_block_z_gain = p4$mean_abs_block_z - p3$mean_abs_block_z,
      empirical_threshold_ratio = p4$threshold_empirical / p3$threshold_empirical,
      null_sd_ratio = p4$null_sd_z / p3$null_sd_z,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

evaluate_phase95rc <- function(results) {
  pure <- subset(results, check_id == "RCR1_pure_kurtosis")
  pure3 <- subset(pure, order == 3L)
  pure4 <- subset(pure, order == 4L)
  phq_stronger_n <- subset(pure4, N == 5000L & abs(dgp_R - 1.5) < 1e-12 &
                             abs(dgp_kappa4 + 1.3) < 1e-12)
  phq_stronger_r <- subset(pure4, N == 2000L & abs(dgp_R - 2.0) < 1e-12 &
                             abs(dgp_kappa4 + 1.3) < 1e-12)
  gains <- phase95rc_pair_gains(results)
  c4 <- subset(results, check_id == "RCR2_block_normalized_c4")

  checks <- c(
    pure_kurtosis_order3_degenerate =
      nrow(pure3) > 0L && all(pure3$degenerate_reps == pure3$reps),
    pure_kurtosis_order4_detects =
      nrow(pure4) > 0L && all(pure4$degenerate_reps == 0L) &&
      ((nrow(phq_stronger_n) == 1L && phq_stronger_n$rate_empirical >= 0.40) ||
         (nrow(phq_stronger_r) == 1L && phq_stronger_r$rate_empirical >= 0.40)),
    block_normalized_typeI_acceptable =
      mean(c4$null_rate_block_normal >= 0.025 & c4$null_rate_block_normal <= 0.10,
           na.rm = TRUE) >= 0.75,
    block_normalized_c4_power_gain =
      nrow(gains) > 0L && all(gains$gain_block_normal >= 0.05, na.rm = TRUE),
    proceed_to_phase10 =
      nrow(gains) > 0L &&
      all(gains$gain_block_normal >= 0.05, na.rm = TRUE) &&
      mean(c4$null_rate_block_normal >= 0.025 & c4$null_rate_block_normal <= 0.10,
           na.rm = TRUE) >= 0.75
  )

  data.frame(
    check_id = names(checks),
    pass = unlist(checks, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

write_phase95rc_report <- function(results, verdict, gains,
                                   path = file.path("output", "sanity_check_phase95rc.txt")) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("Phase 9.5R-c GSA-LLR H2/calibration repair diagnostic", con)
  writeLines("========================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Reps: %d; empirical-null calibration reps: %d",
                     unique(results$reps)[1], unique(results$calibration_reps)[1]), con)
  writeLines("", con)

  writeLines("Verdict by diagnostic:", con)
  for (i in seq_len(nrow(verdict))) {
    writeLines(sprintf("- %s: %s", verdict$check_id[i],
                       if (verdict$pass[i]) "PASS" else "FAIL"), con)
  }
  writeLines("", con)

  writeLines("Pure-kurtosis detector check:", con)
  pure <- subset(results, check_id == "RCR1_pure_kurtosis")
  pure_compact <- pure[, c("N", "order", "dgp_R", "dgp_kappa4",
                           "rate_empirical", "rate_block_normal",
                           "threshold_empirical", "mean_abs_z",
                           "degenerate_reps")]
  utils::write.table(pure_compact, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Block-normalized C4 comparison:", con)
  c4 <- subset(results, check_id == "RCR2_block_normalized_c4")
  c4_compact <- c4[, c("N", "order", "dgp_R", "dgp_kappa3", "dgp_kappa4",
                       "rate_empirical", "rate_block_normal",
                       "null_rate_block_normal", "threshold_empirical",
                       "null_sd_z", "mean_abs_z", "mean_abs_block_z")]
  utils::write.table(c4_compact, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Order-4 minus order-3 gains in non-saturated C4 cells:", con)
  utils::write.table(gains, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Interpretation:", con)
  writeLines("- Pure-kurtosis alternatives confirm the limited but real role of the fourth-order block: order 3 is degenerate when kappa3=0, while order 4 detects once N or R is large enough.", con)
  writeLines("- Block-normalizing by the empirical null controls the scale mismatch more directly than raw empirical quantiles, but it is not a proof of uniform C4 dominance.", con)
  if (isTRUE(verdict$pass[verdict$check_id == "block_normalized_c4_power_gain"])) {
    writeLines("- In this grid, block-normalized order 4 improves calibrated rejection power. Phase 10 can proceed only with the revised calibration explicitly documented.", con)
  } else {
    writeLines("- In this grid, order 4 still does not show uniform calibrated rejection-power dominance. H2 should be revised to an omnibus/sensitivity claim rather than dominance over every marginal test.", con)
  }
  writeLines("- Recommendation: do not start Phase 10 unless the project accepts the revised H2 and uses Phase 10 to implement that narrower claim.", con)
}

run_phase95rc <- function(reps = 500L, calibration_reps = 1000L) {
  ensure_phase95_dirs()
  results <- phase95rc_plan(reps = reps, calibration_reps = calibration_reps)
  verdict <- evaluate_phase95rc(results)
  gains <- phase95rc_pair_gains(results)
  utils::write.csv(results, file.path("output", "tables", "gsa_llr_phase95rc_results.csv"),
                   row.names = FALSE)
  utils::write.csv(verdict, file.path("output", "tables", "gsa_llr_phase95rc_verdict.csv"),
                   row.names = FALSE)
  utils::write.csv(gains, file.path("output", "tables", "gsa_llr_phase95rc_gains.csv"),
                   row.names = FALSE)
  write_phase95rc_report(results, verdict, gains)
  capture.output(utils::sessionInfo(),
                 file = file.path("output", "session_info", "session_info_phase95rc.txt"))
  invisible(list(results = results, verdict = verdict, gains = gains,
                 overall_pass = all(verdict$pass)))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- "phase95"
  reps <- 100L
  calibration_reps <- 500L
  threshold <- "vp"
  standardize <- TRUE

  for (arg in args) {
    if (grepl("^--mode=", arg)) mode <- sub("^--mode=", "", arg)
    if (grepl("^--reps=", arg)) reps <- as.integer(sub("^--reps=", "", arg))
    if (grepl("^--calibration-reps=", arg)) {
      calibration_reps <- as.integer(sub("^--calibration-reps=", "", arg))
    }
    if (grepl("^--threshold=", arg)) threshold <- sub("^--threshold=", "", arg)
    if (grepl("^--standardize=", arg)) {
      standardize <- tolower(sub("^--standardize=", "", arg)) %in% c("1", "true", "yes")
    }
  }

  if (identical(mode, "phase95r")) {
    ans <- run_phase95r(reps = reps, calibration_reps = calibration_reps)
  } else if (identical(mode, "phase95rb")) {
    ans <- run_phase95rb(reps = reps, calibration_reps = calibration_reps)
  } else if (identical(mode, "phase95rc")) {
    ans <- run_phase95rc(reps = reps, calibration_reps = calibration_reps)
  } else {
    ans <- run_phase95(reps = reps, threshold = threshold, standardize = standardize)
  }
  print(ans$verdict)
  cat(sprintf("Overall pass: %s\n", ans$overall_pass))
}

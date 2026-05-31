# GSA-LLR repair sprint: empirical high-order basis, Gaussian fixed-direction
# Type-I, and moment audit.
#
# This is intentionally diagnostic. Closed-form F/Y are available only for the
# s=3/s=4 derivation; for s=5/s=6 we estimate F and Y from calibration samples.
#
# Full run:
#   Rscript scripts/15_gsa_llr_repair_sweep.R
#
# Smoke run:
#   Rscript scripts/15_gsa_llr_repair_sweep.R --quick

source(file.path("scripts", "13_gsa_llr_detector.R"))

repair_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

repair_arg_value <- function(args, name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", name, "="), "", hit[[1L]])
}

repair_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  quick <- "--quick" %in% args || identical(Sys.getenv("GSA_REPAIR_QUICK"), "1")
  list(
    quick = quick,
    reps = as.integer(Sys.getenv(
      "GSA_REPAIR_REPS",
      repair_arg_value(args, "--reps", if (quick) 80L else 300L)
    )),
    calibration_reps = as.integer(Sys.getenv(
      "GSA_REPAIR_CALIBRATION_REPS",
      repair_arg_value(args, "--calibration-reps", if (quick) 160L else 600L)
    )),
    train_n = as.integer(Sys.getenv(
      "GSA_REPAIR_TRAIN_N",
      repair_arg_value(args, "--train-n", if (quick) 5000L else 20000L)
    )),
    audit_n = as.integer(Sys.getenv(
      "GSA_REPAIR_AUDIT_N",
      repair_arg_value(args, "--audit-n", if (quick) 10000L else 50000L)
    )),
    alpha = 0.05,
    ridge_factor = 1e-6,
    seed_base = 20260527L
  )
}

repair_exponents <- function(max_degree) {
  max_degree <- as.integer(max_degree)
  if (max_degree < 3L) stop("max_degree must be >= 3")
  rows <- list()
  for (degree in seq.int(3L, max_degree)) {
    for (a in seq.int(0L, floor((degree - 1L) / 2L))) {
      b <- degree - a
      if (a < b) {
        rows[[length(rows) + 1L]] <- data.frame(
          degree = degree,
          a = a,
          b = b,
          name = sprintf("d%d_a%d_b%d", degree, a, b),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

repair_standardize_pair <- function(x1, x2) {
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  x1 <- x1 - mean(x1)
  x2 <- x2 - mean(x2)
  s1 <- stats::sd(x1)
  s2 <- stats::sd(x2)
  if (!is.finite(s1) || !is.finite(s2) || s1 <= 0 || s2 <= 0) {
    stop("Degenerate score vector in repair basis")
  }
  list(x1 = x1 / s1, x2 = x2 / s2)
}

repair_basis <- function(x1, x2, max_degree = 6L) {
  xs <- repair_standardize_pair(x1, x2)
  ex <- repair_exponents(max_degree)
  out <- matrix(NA_real_, nrow = length(xs$x1), ncol = nrow(ex))
  for (j in seq_len(nrow(ex))) {
    a <- ex$a[j]
    b <- ex$b[j]
    out[, j] <- (xs$x1^a) * (xs$x2^b) - (xs$x1^b) * (xs$x2^a)
  }
  colnames(out) <- ex$name
  out
}

repair_dgp_with_u <- function(N, R = 1, kappa3 = 0, kappa4 = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  T_score <- stats::rnorm(N)
  U <- gsa_llr_component(N, kappa3 = kappa3, kappa4 = kappa4)
  W1 <- stats::rnorm(N)
  W2 <- stats::rnorm(N)
  data.frame(
    U = U,
    T_score = T_score,
    W1 = W1,
    W2 = W2,
    x1 = T_score + U + W1,
    x2 = T_score + R * U + W2
  )
}

repair_skew <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  z <- (x - mean(x)) / stats::sd(x)
  mean(z^3)
}

repair_excess <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  z <- (x - mean(x)) / stats::sd(x)
  mean(z^4) - 3
}

repair_train <- function(max_degree, R, alt_kappa3, alt_kappa4,
                         null_kappa3, null_kappa4, train_n, seed,
                         ridge_factor = 1e-6) {
  h0 <- repair_dgp_with_u(
    N = train_n, R = 1, kappa3 = null_kappa3, kappa4 = null_kappa4,
    seed = seed + 1L
  )
  h1 <- repair_dgp_with_u(
    N = train_n, R = R, kappa3 = alt_kappa3, kappa4 = alt_kappa4,
    seed = seed + 2L
  )
  b0 <- repair_basis(h0$x1, h0$x2, max_degree = max_degree)
  b1 <- repair_basis(h1$x1, h1$x2, max_degree = max_degree)

  Y <- colMeans(b1) - colMeans(b0)
  F <- stats::cov(b0) + stats::cov(b1)
  diag_mean <- mean(diag(F))
  ridge <- ridge_factor * ifelse(is.finite(diag_mean) && diag_mean > 0, diag_mean, 1)
  F_reg <- F + diag(ridge, nrow(F))
  K <- as.numeric(solve(F_reg, Y))
  J <- as.numeric(crossprod(Y, K))
  cond_F <- kappa(F_reg, exact = FALSE)

  list(
    max_degree = max_degree,
    basis_dim = ncol(b0),
    exponents = repair_exponents(max_degree),
    F = F,
    F_reg = F_reg,
    Y = Y,
    K = K,
    J = J,
    ridge = ridge,
    cond_F = cond_F,
    degenerate_direction = !is.finite(J) || J <= sqrt(.Machine$double.eps)
  )
}

repair_score <- function(x1, x2, fit) {
  if (isTRUE(fit$degenerate_direction)) {
    return(NA_real_)
  }
  b <- repair_basis(x1, x2, max_degree = fit$max_degree)
  lambda <- sum(fit$K * colMeans(b))
  sqrt(length(x1)) * lambda / sqrt(fit$J)
}

repair_z_sample <- function(reps, N, R, kappa3, kappa4, fit, seed) {
  z <- rep(NA_real_, reps)
  for (r in seq_len(reps)) {
    dat <- repair_dgp_with_u(
      N = N, R = R, kappa3 = kappa3, kappa4 = kappa4,
      seed = seed + r
    )
    z[r] <- repair_score(dat$x1, dat$x2, fit)
  }
  z
}

repair_rate <- function(z, threshold) {
  if (!is.finite(threshold) || all(!is.finite(z))) return(NA_real_)
  mean(abs(z) > threshold, na.rm = TRUE)
}

repair_cell <- function(check_id, max_degree, N, R, alt_kappa3, alt_kappa4,
                        null_kappa3, null_kappa4, reps, calibration_reps,
                        train_n, alpha, seed, ridge_factor) {
  fit <- repair_train(
    max_degree = max_degree, R = R,
    alt_kappa3 = alt_kappa3, alt_kappa4 = alt_kappa4,
    null_kappa3 = null_kappa3, null_kappa4 = null_kappa4,
    train_n = train_n, seed = seed, ridge_factor = ridge_factor
  )

  null_z_cal <- repair_z_sample(
    reps = calibration_reps, N = N, R = 1,
    kappa3 = null_kappa3, kappa4 = null_kappa4,
    fit = fit, seed = seed + 100000L
  )
  threshold <- if (all(!is.finite(null_z_cal))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(abs(null_z_cal), 1 - alpha,
                               na.rm = TRUE, names = FALSE))
  }

  null_z_eval <- repair_z_sample(
    reps = reps, N = N, R = 1,
    kappa3 = null_kappa3, kappa4 = null_kappa4,
    fit = fit, seed = seed + 200000L
  )
  alt_z <- repair_z_sample(
    reps = reps, N = N, R = R,
    kappa3 = alt_kappa3, kappa4 = alt_kappa4,
    fit = fit, seed = seed + 300000L
  )

  data.frame(
    check_id = check_id,
    max_degree = max_degree,
    basis_dim = fit$basis_dim,
    N = N,
    R = R,
    alt_kappa3 = alt_kappa3,
    alt_kappa4 = alt_kappa4,
    null_kappa3 = null_kappa3,
    null_kappa4 = null_kappa4,
    reps = reps,
    calibration_reps = calibration_reps,
    train_n = train_n,
    threshold_empirical = threshold,
    null_rate_empirical = repair_rate(null_z_eval, threshold),
    alt_rate_empirical = repair_rate(alt_z, threshold),
    mean_abs_alt_z = mean(abs(alt_z), na.rm = TRUE),
    mean_abs_null_z = mean(abs(null_z_eval), na.rm = TRUE),
    null_sd_z = stats::sd(null_z_eval, na.rm = TRUE),
    J = fit$J,
    cond_F = fit$cond_F,
    ridge = fit$ridge,
    degenerate_direction = fit$degenerate_direction,
    stringsAsFactors = FALSE
  )
}

repair_moment_audit <- function(audit_n, seed_base) {
  grid <- data.frame(
    family = c("gaussian", "skew_only", "pure_kurtosis", "pure_kurtosis",
               "mixed_ad_hoc", "mixed_ad_hoc"),
    R = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.25),
    requested_kappa3 = c(0, 0.75, 0, 0, 0.75, 2.25),
    requested_kappa4 = c(0, 0, -1.3, 3.0, 1.3, 3.0)
  )
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dat <- repair_dgp_with_u(
      N = audit_n,
      R = grid$R[i],
      kappa3 = grid$requested_kappa3[i],
      kappa4 = grid$requested_kappa4[i],
      seed = seed_base + i
    )
    rows[[i]] <- data.frame(
      family = grid$family[i],
      R = grid$R[i],
      requested_kappa3 = grid$requested_kappa3[i],
      requested_kappa4 = grid$requested_kappa4[i],
      U_skew = repair_skew(dat$U),
      U_excess = repair_excess(dat$U),
      x1_skew = repair_skew(dat$x1),
      x1_excess = repair_excess(dat$x1),
      x2_skew = repair_skew(dat$x2),
      x2_excess = repair_excess(dat$x2),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

repair_plan <- function(cfg) {
  rows <- list()

  # Gap 1 from review: Gaussian data, fixed nonzero alternative directions.
  fixed_dirs <- data.frame(
    label = c("gaussian_fixed_skew", "gaussian_fixed_platykurt",
              "gaussian_fixed_leptokurt"),
    R = c(1.5, 1.5, 1.5),
    alt_kappa3 = c(0.75, 0, 0),
    alt_kappa4 = c(0, -1.3, 1.3)
  )
  for (i in seq_len(nrow(fixed_dirs))) {
    for (degree in c(4L, 5L, 6L)) {
      rows[[length(rows) + 1L]] <- repair_cell(
        check_id = fixed_dirs$label[i],
        max_degree = degree,
        N = 2000L,
        R = fixed_dirs$R[i],
        alt_kappa3 = fixed_dirs$alt_kappa3[i],
        alt_kappa4 = fixed_dirs$alt_kappa4[i],
        null_kappa3 = 0,
        null_kappa4 = 0,
        reps = cfg$reps,
        calibration_reps = cfg$calibration_reps,
        train_n = cfg$train_n,
        alpha = cfg$alpha,
        seed = cfg$seed_base + 1000L * i + degree,
        ridge_factor = cfg$ridge_factor
      )
    }
  }

  # Repair target: strict PHQ-like pure-kurtosis cell that failed Phase 11Q.
  for (degree in c(4L, 5L, 6L)) {
    rows[[length(rows) + 1L]] <- repair_cell(
      check_id = "phq_like_pure_kurtosis_repair",
      max_degree = degree,
      N = 2000L,
      R = 1.5,
      alt_kappa3 = 0,
      alt_kappa4 = -1.3,
      null_kappa3 = 0,
      null_kappa4 = -1.3,
      reps = cfg$reps,
      calibration_reps = cfg$calibration_reps,
      train_n = cfg$train_n,
      alpha = cfg$alpha,
      seed = cfg$seed_base + 9000L + degree,
      ridge_factor = cfg$ridge_factor
    )
  }

  do.call(rbind, rows)
}

repair_verdict <- function(results) {
  type1 <- subset(results, grepl("^gaussian_fixed_", check_id))
  phq <- subset(results, check_id == "phq_like_pure_kurtosis_repair")
  phq56 <- subset(phq, max_degree %in% c(5L, 6L))
  checks <- c(
    gaussian_fixed_typeI_control =
      nrow(type1) > 0L && all(type1$null_rate_empirical >= 0.025 &
                                type1$null_rate_empirical <= 0.10,
                              na.rm = TRUE),
    high_order_phq_power_gt_0.40 =
      nrow(phq56) > 0L && any(phq56$alt_rate_empirical > 0.40, na.rm = TRUE),
    high_order_conditioning_acceptable =
      nrow(phq56) > 0L && all(phq56$cond_F < 1e8, na.rm = TRUE)
  )
  data.frame(
    criterion = names(checks),
    pass = unlist(checks, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

repair_write_report <- function(results, verdict, audit, cfg,
                                path = file.path("output", "sanity_check_gsa_repair_sweep.txt")) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("GSA-LLR repair sweep: s=5/s=6, Type-I gaps, moment audit", con)
  writeLines("================================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Quick mode: %s", cfg$quick), con)
  writeLines(sprintf("Reps: %d; calibration reps: %d; train_n: %d; audit_n: %d",
                     cfg$reps, cfg$calibration_reps, cfg$train_n, cfg$audit_n), con)
  writeLines("", con)

  writeLines("Verdict:", con)
  for (i in seq_len(nrow(verdict))) {
    writeLines(sprintf("- %s: %s", verdict$criterion[i],
                       if (verdict$pass[i]) "PASS" else "FAIL"), con)
  }
  writeLines(sprintf("- overall: %s", all(verdict$pass)), con)
  writeLines("", con)

  writeLines("Gaussian fixed-direction Type-I checks:", con)
  type1 <- subset(results, grepl("^gaussian_fixed_", check_id))
  utils::write.table(
    type1[, c("check_id", "max_degree", "N", "R", "alt_kappa3", "alt_kappa4",
              "null_rate_empirical", "alt_rate_empirical", "cond_F", "basis_dim")],
    con, sep = "\t", quote = FALSE, row.names = FALSE
  )
  writeLines("", con)

  writeLines("Strict PHQ-like pure-kurtosis repair target:", con)
  phq <- subset(results, check_id == "phq_like_pure_kurtosis_repair")
  utils::write.table(
    phq[, c("max_degree", "basis_dim", "N", "R", "alt_kappa4",
            "null_rate_empirical", "alt_rate_empirical", "threshold_empirical",
            "mean_abs_alt_z", "null_sd_z", "J", "cond_F")],
    con, sep = "\t", quote = FALSE, row.names = FALSE
  )
  writeLines("", con)

  writeLines("Moment audit:", con)
  utils::write.table(audit, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Interpretation:", con)
  if (any(subset(results, check_id == "phq_like_pure_kurtosis_repair" &
                 max_degree %in% c(5L, 6L))$alt_rate_empirical > 0.40, na.rm = TRUE)) {
    writeLines("- Empirical s=5/s=6 repair reaches the PHQ-like power gate in this diagnostic run. This is a candidate for a follow-up calibration sprint, not yet a symbolic theorem.", con)
  } else {
    writeLines("- Empirical s=5/s=6 repair does not reach the PHQ-like power gate in this diagnostic run. GSA-LLR should remain diagnostic/operating-domain unless a different basis or model repair is introduced.", con)
  }
  writeLines("- Mixed C3+C4 rows in prior Phase 11Q should remain exploratory unless the moment audit confirms requested and realized cumulants are close.", con)
}

run_repair_sweep <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- repair_config(args)
  repair_dirs()
  results <- repair_plan(cfg)
  audit <- repair_moment_audit(audit_n = cfg$audit_n, seed_base = cfg$seed_base + 500000L)
  verdict <- repair_verdict(results)

  utils::write.csv(
    results,
    file.path("output", "tables", "gsa_llr_repair_sweep_results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    subset(results, grepl("^gaussian_fixed_", check_id)),
    file.path("output", "tables", "gsa_llr_repair_gaussian_typeI.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    audit,
    file.path("output", "tables", "gsa_llr_repair_moment_audit.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    verdict,
    file.path("output", "tables", "gsa_llr_repair_verdict.csv"),
    row.names = FALSE
  )
  repair_write_report(results, verdict, audit, cfg)
  capture.output(
    utils::sessionInfo(),
    file = file.path("output", "session_info", "session_info_gsa_repair_sweep.txt")
  )

  invisible(list(results = results, audit = audit, verdict = verdict,
                 overall_pass = all(verdict$pass), config = cfg))
}

if (sys.nframe() == 0L) {
  ans <- run_repair_sweep()
  if (!isTRUE(ans$overall_pass)) {
    quit(status = 1L)
  }
}

# Phase 11Q compact Monte Carlo for the revised GSA-LLR H2.
#
# Default run matches NEXT_ACTION_SPEC_2026-05-27.md:
#   N in {500, 2000, 5000}
#   R in {1, 1.25, 1.5, 2}
#   kappa3 in {0, 0.75, 2.25}
#   kappa4 in {-1.3, 0, 1.3, 3}
#   reps = 500
#   calibration_reps = 1000
#
# Smoke run:
#   Rscript scripts/14_gsa_llr_mc_compact.R --quick

source(file.path("scripts", "13_gsa_llr_detector.R"))

phase11q_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

phase11q_arg_value <- function(args, name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", name, "="), "", hit[[1L]])
}

phase11q_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  quick <- "--quick" %in% args || identical(Sys.getenv("PHASE11Q_QUICK"), "1")
  default_reps <- if (quick) 80L else 500L
  default_calibration_reps <- if (quick) 160L else 1000L

  reps <- as.integer(Sys.getenv(
    "PHASE11Q_REPS",
    phase11q_arg_value(args, "--reps", default_reps)
  ))
  calibration_reps <- as.integer(Sys.getenv(
    "PHASE11Q_CALIBRATION_REPS",
    phase11q_arg_value(args, "--calibration-reps", default_calibration_reps)
  ))

  list(
    quick = quick,
    reps = reps,
    calibration_reps = calibration_reps,
    alpha = 0.05,
    seed_base = 20260527L,
    N_values = c(500L, 2000L, 5000L),
    R_values = c(1, 1.25, 1.5, 2),
    kappa3_values = c(0, 0.75, 2.25),
    kappa4_values = c(-1.3, 0, 1.3, 3)
  )
}

phase11q_scale_cols <- function(x) {
  x <- sweep(x, 2L, colMeans(x), "-")
  s <- sqrt(colSums(x^2) / (nrow(x) - 1L))
  bad <- !is.finite(s) | s <= 0
  if (any(bad)) {
    stop("Degenerate simulated column encountered")
  }
  sweep(x, 2L, s, "/")
}

phase11q_component_matrix <- function(N, reps, kappa3 = 0, kappa4 = 0) {
  n_total <- N * reps

  if (abs(kappa3) > 1e-12 && abs(kappa4) <= 1e-12) {
    shape <- 4 / (abs(kappa3)^2)
    x <- matrix(stats::rgamma(n_total, shape = shape, scale = 1), nrow = N)
    x <- phase11q_scale_cols(x)
    if (kappa3 < 0) x <- -x
    return(x)
  }

  if (abs(kappa3) <= 1e-12 && kappa4 > 1e-12) {
    df <- 6 / kappa4 + 4
    x <- matrix(stats::rt(n_total, df = df) / sqrt(df / (df - 2)), nrow = N)
    return(phase11q_scale_cols(x))
  }

  if (abs(kappa3) <= 1e-12 && kappa4 < -1e-12) {
    a <- (-6 / kappa4 - 3) / 2
    if (!is.finite(a) || a <= 0) {
      stop("Requested negative excess kurtosis is outside beta(a,a) helper range")
    }
    x <- matrix(stats::rbeta(n_total, shape1 = a, shape2 = a), nrow = N)
    return(phase11q_scale_cols(x))
  }

  if (abs(kappa3) > 1e-12 && abs(kappa4) > 1e-12) {
    skew <- phase11q_component_matrix(N, reps, kappa3 = kappa3, kappa4 = 0)
    kurt <- phase11q_component_matrix(N, reps, kappa3 = 0, kappa4 = kappa4)
    w3 <- abs(kappa3)
    w4 <- sqrt(abs(kappa4))
    return(phase11q_scale_cols(w3 * skew + w4 * kurt))
  }

  matrix(stats::rnorm(n_total), nrow = N)
}

phase11q_pair_matrix <- function(N, reps, R, kappa3, kappa4, seed) {
  set.seed(seed)
  T_score <- matrix(stats::rnorm(N * reps), nrow = N)
  U <- phase11q_component_matrix(N, reps, kappa3 = kappa3, kappa4 = kappa4)
  W1 <- matrix(stats::rnorm(N * reps), nrow = N)
  W2 <- matrix(stats::rnorm(N * reps), nrow = N)
  list(
    x1 = T_score + U + W1,
    x2 = T_score + R * U + W2
  )
}

phase11q_z_matrix <- function(pair, test_R, kappa3, kappa4, order = 4L,
                              standardize = TRUE) {
  order <- gsa_llr_match_order(order)
  wk <- gsa_llr_weights(R = test_R, kappa3 = kappa3, kappa4 = kappa4, order = order)
  reps <- ncol(pair$x1)
  N <- nrow(pair$x1)

  if (isTRUE(wk$degenerate_direction)) {
    return(rep(NA_real_, reps))
  }

  x1 <- pair$x1
  x2 <- pair$x2
  if (isTRUE(standardize)) {
    x1 <- phase11q_scale_cols(x1)
    x2 <- phase11q_scale_cols(x2)
  } else {
    x1 <- sweep(x1, 2L, colMeans(x1), "-")
    x2 <- sweep(x2, 2L, colMeans(x2), "-")
  }

  phi_means <- rbind(
    phi1 = colMeans(x1 * x2^2 - x1^2 * x2),
    phi2 = colMeans(x1^3 - x2^3)
  )
  if (order == 4L) {
    phi_means <- rbind(
      phi_means,
      phi3 = colMeans(x1 * x2^3 - x1^3 * x2),
      phi4 = colMeans(x1^4 - x2^4)
    )
  }

  lambda <- as.numeric(crossprod(wk$K, phi_means))
  sqrt(N) * lambda / sqrt(wk$J)
}

phase11q_max_abs <- function(z3, z4) {
  abs_mat <- cbind(abs(z3), abs(z4))
  apply(abs_mat, 1L, function(row) {
    if (all(!is.finite(row))) NA_real_ else max(row, na.rm = TRUE)
  })
}

phase11q_rate <- function(z, threshold) {
  if (!is.finite(threshold) || all(!is.finite(z))) {
    return(NA_real_)
  }
  mean(abs(z) > threshold, na.rm = TRUE)
}

phase11q_binary_rate <- function(rejects) {
  if (length(rejects) == 0L || all(is.na(rejects))) NA_real_ else mean(rejects, na.rm = TRUE)
}

phase11q_cell_family <- function(R, kappa3, kappa4) {
  if (R == 1 && (abs(kappa3) > 1e-12 || abs(kappa4) > 1e-12)) {
    return("null_symmetric_loading")
  }
  if (abs(kappa3) <= 1e-12 && abs(kappa4) <= 1e-12) {
    return("gaussian_no_direction")
  }
  if (R > 1 && abs(kappa3) > 1e-12 && abs(kappa4) <= 1e-12) {
    return("skew_only")
  }
  if (R > 1 && abs(kappa3) <= 1e-12 && abs(kappa4) > 1e-12) {
    return("pure_kurtosis")
  }
  if (R > 1 && abs(kappa3) > 1e-12 && abs(kappa4) > 1e-12) {
    return("mixed")
  }
  "other"
}

phase11q_cell <- function(cell_id, N, R, kappa3, kappa4, cfg) {
  start <- proc.time()[["elapsed"]]
  test_R <- if (R <= 1) 1.25 else R
  seed <- cfg$seed_base + 10000L * cell_id

  alt_pair <- phase11q_pair_matrix(
    N = N, reps = cfg$reps, R = R,
    kappa3 = kappa3, kappa4 = kappa4, seed = seed + 1L
  )
  null_pair <- phase11q_pair_matrix(
    N = N, reps = cfg$calibration_reps, R = 1,
    kappa3 = kappa3, kappa4 = kappa4, seed = seed + 2L
  )

  z3_alt <- phase11q_z_matrix(alt_pair, test_R, kappa3, kappa4, order = 3L)
  z4_alt <- phase11q_z_matrix(alt_pair, test_R, kappa3, kappa4, order = 4L)
  z3_null <- phase11q_z_matrix(null_pair, test_R, kappa3, kappa4, order = 3L)
  z4_null <- phase11q_z_matrix(null_pair, test_R, kappa3, kappa4, order = 4L)

  threshold_normal <- gsa_llr_threshold(cfg$alpha, "normal_approx")
  threshold_vp <- gsa_llr_threshold(cfg$alpha, "vp_bound")
  threshold_normal_bonf <- gsa_llr_threshold(cfg$alpha / 2, "normal_approx")
  threshold_vp_bonf <- gsa_llr_threshold(cfg$alpha / 2, "vp_bound")

  threshold_emp3 <- if (all(!is.finite(z3_null))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(abs(z3_null), 1 - cfg$alpha, na.rm = TRUE, names = FALSE))
  }
  threshold_emp4 <- if (all(!is.finite(z4_null))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(abs(z4_null), 1 - cfg$alpha, na.rm = TRUE, names = FALSE))
  }

  max_alt <- phase11q_max_abs(z3_alt, z4_alt)
  max_null <- phase11q_max_abs(z3_null, z4_null)
  threshold_joint <- if (all(!is.finite(max_null))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(max_null, 1 - cfg$alpha, na.rm = TRUE, names = FALSE))
  }

  omni_normal_alt <- (is.finite(z3_alt) & abs(z3_alt) > threshold_normal_bonf) |
    (is.finite(z4_alt) & abs(z4_alt) > threshold_normal_bonf)
  omni_vp_alt <- (is.finite(z3_alt) & abs(z3_alt) > threshold_vp_bonf) |
    (is.finite(z4_alt) & abs(z4_alt) > threshold_vp_bonf)
  omni_normal_null <- (is.finite(z3_null) & abs(z3_null) > threshold_normal_bonf) |
    (is.finite(z4_null) & abs(z4_null) > threshold_normal_bonf)
  omni_vp_null <- (is.finite(z3_null) & abs(z3_null) > threshold_vp_bonf) |
    (is.finite(z4_null) & abs(z4_null) > threshold_vp_bonf)

  elapsed <- proc.time()[["elapsed"]] - start

  data.frame(
    cell_id = cell_id,
    family = phase11q_cell_family(R, kappa3, kappa4),
    N = N,
    R = R,
    test_R = test_R,
    kappa3 = kappa3,
    kappa4 = kappa4,
    reps = cfg$reps,
    calibration_reps = cfg$calibration_reps,
    alpha = cfg$alpha,
    order3_degenerate = all(!is.finite(z3_alt)),
    order4_degenerate = all(!is.finite(z4_alt)),
    threshold_normal = threshold_normal,
    threshold_vp = threshold_vp,
    threshold_empirical_order3 = threshold_emp3,
    threshold_empirical_order4 = threshold_emp4,
    threshold_empirical_joint = threshold_joint,
    order3_rate_normal = phase11q_rate(z3_alt, threshold_normal),
    order3_rate_vp = phase11q_rate(z3_alt, threshold_vp),
    order3_rate_empirical = phase11q_rate(z3_alt, threshold_emp3),
    order4_rate_normal = phase11q_rate(z4_alt, threshold_normal),
    order4_rate_vp = phase11q_rate(z4_alt, threshold_vp),
    order4_rate_empirical = phase11q_rate(z4_alt, threshold_emp4),
    omnibus_rate_normal_bonferroni = phase11q_binary_rate(omni_normal_alt),
    omnibus_rate_vp_bonferroni = phase11q_binary_rate(omni_vp_alt),
    omnibus_rate_empirical_joint = phase11q_rate(max_alt, threshold_joint),
    order3_null_rate_normal = phase11q_rate(z3_null, threshold_normal),
    order3_null_rate_vp = phase11q_rate(z3_null, threshold_vp),
    order3_null_rate_empirical = phase11q_rate(z3_null, threshold_emp3),
    order4_null_rate_normal = phase11q_rate(z4_null, threshold_normal),
    order4_null_rate_vp = phase11q_rate(z4_null, threshold_vp),
    order4_null_rate_empirical = phase11q_rate(z4_null, threshold_emp4),
    omnibus_null_rate_normal_bonferroni = phase11q_binary_rate(omni_normal_null),
    omnibus_null_rate_vp_bonferroni = phase11q_binary_rate(omni_vp_null),
    omnibus_null_rate_empirical_joint = phase11q_rate(max_null, threshold_joint),
    mean_abs_z3 = mean(abs(z3_alt), na.rm = TRUE),
    mean_abs_z4 = mean(abs(z4_alt), na.rm = TRUE),
    null_sd_z3 = stats::sd(z3_null, na.rm = TRUE),
    null_sd_z4 = stats::sd(z4_null, na.rm = TRUE),
    runtime_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

phase11q_summarize <- function(results) {
  safe_median <- function(x) {
    if (all(!is.finite(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  }
  safe_max <- function(x) {
    if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE)
  }
  families <- sort(unique(results$family))
  rows <- lapply(families, function(fam) {
    x <- results[results$family == fam, ]
    data.frame(
      family = fam,
      cells = nrow(x),
      median_order3_empirical = safe_median(x$order3_rate_empirical),
      median_order4_empirical = safe_median(x$order4_rate_empirical),
      median_omnibus_empirical_joint = safe_median(x$omnibus_rate_empirical_joint),
      median_omnibus_null_empirical_joint = safe_median(x$omnibus_null_rate_empirical_joint),
      max_omnibus_empirical_joint = safe_max(x$omnibus_rate_empirical_joint),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

phase11q_verdict <- function(results) {
  null_cells <- results[results$family == "null_symmetric_loading", ]
  null_ok <- if (nrow(null_cells) == 0L) {
    FALSE
  } else {
    mean(null_cells$omnibus_rate_empirical_joint >= 0.025 &
           null_cells$omnibus_rate_empirical_joint <= 0.10,
         na.rm = TRUE) >= 0.90
  }

  skew <- results[results$family == "skew_only", ]
  skew_groups <- split(skew, paste(skew$N, skew$kappa3, skew$kappa4, sep = "_"))
  skew_monotone <- if (length(skew_groups) == 0L) {
    FALSE
  } else {
    all(vapply(skew_groups, function(g) {
      g <- g[order(g$R), ]
      all(diff(g$order3_rate_empirical) >= -0.03)
    }, logical(1L)))
  }

  phq_cell <- results[results$family == "pure_kurtosis" &
                        results$N == 2000L &
                        abs(results$R - 1.5) < 1e-12 &
                        abs(results$kappa3) < 1e-12 &
                        abs(results$kappa4 + 1.3) < 1e-12, ]
  pure_kurtosis_ok <- nrow(phq_cell) == 1L &&
    is.finite(phq_cell$order4_rate_empirical) &&
    phq_cell$order4_rate_empirical > 0.40

  skew_alt <- results[results$family == "skew_only" & results$R > 1, ]
  kurt_alt <- results[results$family == "pure_kurtosis" & results$R > 1, ]
  omnibus_covers_families <- nrow(skew_alt) > 0L && nrow(kurt_alt) > 0L &&
    stats::median(skew_alt$omnibus_rate_empirical_joint, na.rm = TRUE) >
      stats::median(null_cells$omnibus_rate_empirical_joint, na.rm = TRUE) &&
    stats::median(kurt_alt$omnibus_rate_empirical_joint, na.rm = TRUE) >
      stats::median(null_cells$omnibus_rate_empirical_joint, na.rm = TRUE)

  verdict <- data.frame(
    criterion = c(
      "typeI_acceptable_in_90pct_null_cells",
      "order3_power_monotone_in_R_for_skewness",
      "order4_power_gt_0.40_in_phq_like_pure_kurtosis_cell",
      "omnibus_covers_skewness_and_kurtosis_families_without_dominance_claim"
    ),
    pass = c(
      null_ok,
      skew_monotone,
      pure_kurtosis_ok,
      omnibus_covers_families
    ),
    stringsAsFactors = FALSE
  )
  verdict
}

phase11q_write_report <- function(results, summary, verdict, cfg,
                                  path = file.path("output", "sanity_check_phase11q.txt")) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("Phase 11Q compact Monte Carlo for revised GSA-LLR H2", con)
  writeLines("=====================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Quick mode: %s", cfg$quick), con)
  writeLines(sprintf("Reps: %d; empirical-null calibration reps: %d", cfg$reps, cfg$calibration_reps), con)
  writeLines(sprintf("Grid cells: %d", nrow(results)), con)
  writeLines(sprintf("Total runtime seconds: %.2f", sum(results$runtime_sec, na.rm = TRUE)), con)
  writeLines("", con)

  writeLines("Verdict:", con)
  for (i in seq_len(nrow(verdict))) {
    writeLines(sprintf("- %s: %s", verdict$criterion[i],
                       if (verdict$pass[i]) "PASS" else "FAIL"), con)
  }
  writeLines(sprintf("- overall: %s", all(verdict$pass)), con)
  writeLines("", con)

  writeLines("Summary by family:", con)
  utils::write.table(summary, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  null_cells <- results[results$family == "null_symmetric_loading", ]
  if (nrow(null_cells) > 0L) {
    ok_share <- mean(null_cells$omnibus_rate_empirical_joint >= 0.025 &
                       null_cells$omnibus_rate_empirical_joint <= 0.10,
                     na.rm = TRUE)
    writeLines(sprintf("Null-cell acceptable share: %.3f", ok_share), con)
  }

  phq_cell <- results[results$family == "pure_kurtosis" &
                        results$N == 2000L &
                        abs(results$R - 1.5) < 1e-12 &
                        abs(results$kappa4 + 1.3) < 1e-12, ]
  if (nrow(phq_cell) == 1L) {
    writeLines(sprintf(
      "PHQ-like pure-kurtosis cell order-4 empirical power: %.3f",
      phq_cell$order4_rate_empirical
    ), con)
  }

  writeLines("", con)
  writeLines("Interpretation:", con)
  if (all(verdict$pass)) {
    writeLines("- Revised H2 has a go signal for manuscript outlining: empirical-null omnibus calibration is usable in this compact grid.", con)
    writeLines("- The result supports an omnibus/sensitivity story, not a uniform dominance story.", con)
  } else {
    writeLines("- Revised H2 is not ready as a manuscript claim under this compact grid. Demote or narrow the GSA-LLR track before PHQ/manuscript expansion.", con)
  }
}

run_phase11q <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- phase11q_config(args)
  phase11q_dirs()

  grid <- expand.grid(
    N = cfg$N_values,
    R = cfg$R_values,
    kappa3 = cfg$kappa3_values,
    kappa4 = cfg$kappa4_values,
    KEEP.OUT.ATTRS = FALSE
  )
  grid <- grid[order(grid$N, grid$R, grid$kappa3, grid$kappa4), ]
  rownames(grid) <- NULL

  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    cat(sprintf(
      "Phase11Q cell %03d/%03d: N=%d R=%.2f k3=%.2f k4=%.2f\n",
      i, nrow(grid), grid$N[i], grid$R[i], grid$kappa3[i], grid$kappa4[i]
    ))
    results[[i]] <- phase11q_cell(
      cell_id = i,
      N = grid$N[i],
      R = grid$R[i],
      kappa3 = grid$kappa3[i],
      kappa4 = grid$kappa4[i],
      cfg = cfg
    )
    if (i %% 12L == 0L) {
      gc(verbose = FALSE)
    }
  }

  results <- do.call(rbind, results)
  summary <- phase11q_summarize(results)
  verdict <- phase11q_verdict(results)

  utils::write.csv(
    results,
    file.path("output", "tables", "gsa_llr_phase11q_results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path("output", "tables", "gsa_llr_phase11q_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    verdict,
    file.path("output", "tables", "gsa_llr_phase11q_verdict.csv"),
    row.names = FALSE
  )
  phase11q_write_report(results, summary, verdict, cfg)
  capture.output(
    utils::sessionInfo(),
    file = file.path("output", "session_info", "session_info_phase11q.txt")
  )

  invisible(list(results = results, summary = summary, verdict = verdict,
                 overall_pass = all(verdict$pass), config = cfg))
}

if (sys.nframe() == 0L) {
  ans <- run_phase11q()
  if (!isTRUE(ans$overall_pass)) {
    quit(status = 1L)
  }
}

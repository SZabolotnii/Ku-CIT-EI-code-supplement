# PMM3-style symmetric probe for the W&S pure-kurtosis cell.
#
# This is a deliberately gated sanity check, not a manuscript-ready PMM3
# derivation. It asks whether a PMM3-style antisymmetric control-variate
# correction can improve the naive Delta c4 statistic in the symmetric
# platykurtic setting:
#   N = 2000, R = 1.5, kappa3 = 0, kappa4 = -1.3.
#
# The lesson from PMM2 is encoded as a hard gate: H0 variance reduction alone is
# insufficient; the estimator must also retain H1 signal and power.

source(file.path("scripts", "02_naive_estimators.R"))
source(file.path("scripts", "13_gsa_llr_detector.R"))

pmm3_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

pmm3_arg_value <- function(args, name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", name, "="), "", hit[[1L]])
}

pmm3_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  quick <- "--quick" %in% args || identical(Sys.getenv("PMM3_SYMMETRIC_QUICK"), "1")
  list(
    quick = quick,
    N = as.integer(pmm3_arg_value(args, "N", "2000")),
    R_alt = as.numeric(pmm3_arg_value(args, "R_alt", "1.5")),
    kappa3 = 0,
    kappa4 = as.numeric(pmm3_arg_value(args, "kappa4", "-1.3")),
    alpha = 0.05,
    calibration_reps = if (quick) 160L else 800L,
    eval_reps = if (quick) 160L else 800L,
    reference_n = if (quick) 50000L else 300000L,
    seed_base = 20260527L + 1800L,
    ridge = 1e-4
  )
}

pmm3_delta4_terms <- function(x1, x2) {
  stopifnot(length(x1) == length(x2), length(x1) >= 10L)
  e1 <- as.numeric(x1) - mean(x1)
  e2 <- as.numeric(x2) - mean(x2)
  m12 <- mean(e1 * e2)
  target <- e1 * e2^3 - e1^3 * e2 - 3 * m12 * (e2^2 - e1^2)

  aug <- cbind(
    psi3 = e1^3 - e2^3,
    psi5 = e1^5 - e2^5,
    psi41 = e1^4 * e2 - e1 * e2^4,
    psi32 = e1^3 * e2^2 - e1^2 * e2^3
  )
  storage.mode(aug) <- "double"
  list(target = target, aug = aug)
}

pmm3_delta_c4_cv <- function(x1, x2, ridge = 1e-4) {
  terms <- pmm3_delta4_terms(x1, x2)
  target <- terms$target
  aug <- terms$aug
  keep <- vapply(seq_len(ncol(aug)), function(j) {
    stats::sd(aug[, j], na.rm = TRUE) > 1e-12
  }, logical(1L))
  aug <- aug[, keep, drop = FALSE]
  if (ncol(aug) == 0L) {
    return(mean(target))
  }

  target_c <- target - mean(target)
  aug_c <- sweep(aug, 2L, colMeans(aug), "-")
  F <- crossprod(aug_c) / max(1L, length(target) - 1L)
  B <- as.numeric(crossprod(aug_c, target_c) / max(1L, length(target) - 1L))
  ridge_abs <- ridge * max(1, mean(diag(F)))
  K <- tryCatch(
    -solve(F + ridge_abs * diag(ncol(F)), B),
    error = function(e) rep(0, ncol(F))
  )
  mean(target) + sum(K * colMeans(aug))
}

pmm3_sample_cumulants <- function(x) {
  z <- as.numeric(x) - mean(x)
  m2 <- mean(z^2)
  m3 <- mean(z^3)
  m4 <- mean(z^4)
  m6 <- mean(z^6)
  gamma3 <- m3 / (m2^(3 / 2))
  gamma4 <- m4 / (m2^2) - 3
  gamma6 <- m6 / (m2^3) - 15 * (m4 / m2^2) + 30
  denom <- 6 + 9 * gamma4 + gamma6
  g3 <- if (is.finite(denom) && denom > 0) 1 - gamma4^2 / denom else NA_real_
  data.frame(
    gamma3 = gamma3,
    gamma4 = gamma4,
    gamma6 = gamma6,
    g3_theory_like = g3,
    stringsAsFactors = FALSE
  )
}

pmm3_generate <- function(cfg, R, kappa4, seed) {
  gsa_llr_dgp(
    N = cfg$N, R = R, kappa3 = cfg$kappa3, kappa4 = kappa4, seed = seed
  )
}

pmm3_estimator_values <- function(cfg, reps, R, kappa4, seed_offset) {
  out <- data.frame(
    rep = seq_len(reps),
    scenario_R = R,
    scenario_kappa4 = kappa4,
    naive_delta_c4 = NA_real_,
    pmm3_delta_c4 = NA_real_,
    stringsAsFactors = FALSE
  )
  for (r in seq_len(reps)) {
    dat <- pmm3_generate(
      cfg = cfg, R = R, kappa4 = kappa4,
      seed = cfg$seed_base + seed_offset + r
    )
    out$naive_delta_c4[[r]] <- delta_c4_naive(dat$x1, dat$x2)
    out$pmm3_delta_c4[[r]] <- pmm3_delta_c4_cv(dat$x1, dat$x2, ridge = cfg$ridge)
  }
  out
}

pmm3_thresholds <- function(calibration, cfg) {
  data.frame(
    method = c("naive_delta_c4", "pmm3_delta_c4"),
    threshold = c(
      as.numeric(stats::quantile(abs(calibration$naive_delta_c4), 1 - cfg$alpha, names = FALSE)),
      as.numeric(stats::quantile(abs(calibration$pmm3_delta_c4), 1 - cfg$alpha, names = FALSE))
    ),
    stringsAsFactors = FALSE
  )
}

pmm3_eval_method <- function(method, threshold, null_eval, alt_eval, nuisance_eval) {
  data.frame(
    method = method,
    threshold = threshold,
    typeI_H0_ng = mean(abs(null_eval[[method]]) > threshold),
    power_H1_ng = mean(abs(alt_eval[[method]]) > threshold),
    nuisance_H1_gauss = mean(abs(nuisance_eval[[method]]) > threshold),
    mean_H0_ng = mean(null_eval[[method]]),
    mean_H1_ng = mean(alt_eval[[method]]),
    mean_H1_gauss = mean(nuisance_eval[[method]]),
    var_H0_ng = stats::var(null_eval[[method]]),
    var_H1_ng = stats::var(alt_eval[[method]]),
    stringsAsFactors = FALSE
  )
}

pmm3_reference <- function(cfg) {
  ref_h0 <- gsa_llr_dgp(
    N = cfg$reference_n, R = 1, kappa3 = cfg$kappa3,
    kappa4 = cfg$kappa4, seed = cfg$seed_base + 8000000L
  )
  ref_h1 <- gsa_llr_dgp(
    N = cfg$reference_n, R = cfg$R_alt, kappa3 = cfg$kappa3,
    kappa4 = cfg$kappa4, seed = cfg$seed_base + 8100000L
  )
  ref_gauss <- gsa_llr_dgp(
    N = cfg$reference_n, R = cfg$R_alt, kappa3 = 0,
    kappa4 = 0, seed = cfg$seed_base + 8200000L
  )
  marg <- pmm3_sample_cumulants(c(ref_h0$x1, ref_h0$x2))
  data.frame(
    scenario = c("H0_ng", "H1_ng", "H1_gauss_nuisance"),
    R = c(1, cfg$R_alt, cfg$R_alt),
    kappa4 = c(cfg$kappa4, cfg$kappa4, 0),
    naive_delta_c4_reference = c(
      delta_c4_naive(ref_h0$x1, ref_h0$x2),
      delta_c4_naive(ref_h1$x1, ref_h1$x2),
      delta_c4_naive(ref_gauss$x1, ref_gauss$x2)
    ),
    pmm3_delta_c4_reference = c(
      pmm3_delta_c4_cv(ref_h0$x1, ref_h0$x2, ridge = cfg$ridge),
      pmm3_delta_c4_cv(ref_h1$x1, ref_h1$x2, ridge = cfg$ridge),
      pmm3_delta_c4_cv(ref_gauss$x1, ref_gauss$x2, ridge = cfg$ridge)
    ),
    marginal_gamma3_H0_ng = marg$gamma3,
    marginal_gamma4_H0_ng = marg$gamma4,
    marginal_gamma6_H0_ng = marg$gamma6,
    marginal_g3_theory_like_H0_ng = marg$g3_theory_like,
    reference_n = cfg$reference_n,
    stringsAsFactors = FALSE
  )
}

pmm3_verdict <- function(results, reference, cfg) {
  naive <- results[results$method == "naive_delta_c4", ]
  pmm3 <- results[results$method == "pmm3_delta_c4", ]
  ref_alt <- reference[reference$scenario == "H1_ng", ]
  alt_ratio <- pmm3$mean_H1_ng / naive$mean_H1_ng
  reference_ratio <- ref_alt$pmm3_delta_c4_reference / ref_alt$naive_delta_c4_reference
  variance_reduction <- naive$var_H0_ng / pmm3$var_H0_ng
  checks <- data.frame(
    criterion = c(
      "H0 variance reduction vs naive",
      "Type-I controlled for PMM3-style estimator",
      "H1 signal retained vs naive mean",
      "Large-N H1 reference retained",
      "Power beats naive by at least 5pp",
      "Gaussian nuisance not worse than 10pct"
    ),
    value = c(
      variance_reduction,
      pmm3$typeI_H0_ng,
      alt_ratio,
      reference_ratio,
      pmm3$power_H1_ng - naive$power_H1_ng,
      pmm3$nuisance_H1_gauss
    ),
    pass = c(
      is.finite(variance_reduction) && variance_reduction > 1.0,
      is.finite(pmm3$typeI_H0_ng) && pmm3$typeI_H0_ng <= 0.10,
      is.finite(alt_ratio) && alt_ratio >= 0.80 && alt_ratio <= 1.20,
      is.finite(reference_ratio) && reference_ratio >= 0.80 && reference_ratio <= 1.20,
      is.finite(pmm3$power_H1_ng - naive$power_H1_ng) &&
        pmm3$power_H1_ng - naive$power_H1_ng >= 0.05,
      is.finite(pmm3$nuisance_H1_gauss) && pmm3$nuisance_H1_gauss <= 0.10
    ),
    stringsAsFactors = FALSE
  )
  overall <- data.frame(
    criterion = "pmm3_symmetric_probe_passes_all_gates",
    pass = all(checks$pass),
    stringsAsFactors = FALSE
  )
  list(checks = checks, overall = overall)
}

pmm3_write_report <- function(results, thresholds, verdict, reference, cfg) {
  con <- file(file.path("output", "sanity_check_pmm3_symmetric_probe.txt"), open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("PMM3-style symmetric Delta c4 probe", con)
  writeLines("================================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Quick mode: %s", cfg$quick), con)
  writeLines(sprintf(
    "Target: N=%d; R_alt=%.2f; kappa3=%.2f; kappa4=%.2f; alpha=%.2f",
    cfg$N, cfg$R_alt, cfg$kappa3, cfg$kappa4, cfg$alpha
  ), con)
  writeLines(sprintf(
    "Calibration/eval reps=%d/%d; reference_n=%d; ridge=%g",
    cfg$calibration_reps, cfg$eval_reps, cfg$reference_n, cfg$ridge
  ), con)
  writeLines("", con)

  writeLines("Thresholds:", con)
  utils::write.table(thresholds, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Results:", con)
  utils::write.table(results, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Large-N reference / cumulant audit:", con)
  utils::write.table(reference, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Gate checks:", con)
  utils::write.table(verdict$checks, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)
  writeLines("Overall verdict:", con)
  utils::write.table(verdict$overall, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  if (isTRUE(verdict$overall$pass[[1L]])) {
    writeLines("Interpretation:", con)
    writeLines("- The PMM3-style symmetric control-variate probe passes the practical gates in this cell. It still needs a formal PMM3 cross-cumulant derivation before manuscript use.", con)
  } else {
    writeLines("Interpretation:", con)
    writeLines("- The PMM3-style symmetric control-variate probe does not pass all gates. Treat PMM3 here as diagnostic unless redesigned.", con)
  }
}

pmm3_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  pmm3_dirs()
  cfg <- pmm3_config(args)

  calibration <- pmm3_estimator_values(
    cfg, cfg$calibration_reps, R = 1, kappa4 = cfg$kappa4,
    seed_offset = 1000000L
  )
  thresholds <- pmm3_thresholds(calibration, cfg)
  null_eval <- pmm3_estimator_values(
    cfg, cfg$eval_reps, R = 1, kappa4 = cfg$kappa4,
    seed_offset = 2000000L
  )
  alt_eval <- pmm3_estimator_values(
    cfg, cfg$eval_reps, R = cfg$R_alt, kappa4 = cfg$kappa4,
    seed_offset = 3000000L
  )
  nuisance_eval <- pmm3_estimator_values(
    cfg, cfg$eval_reps, R = cfg$R_alt, kappa4 = 0,
    seed_offset = 4000000L
  )

  results <- do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) {
    pmm3_eval_method(
      method = thresholds$method[[i]],
      threshold = thresholds$threshold[[i]],
      null_eval = null_eval,
      alt_eval = alt_eval,
      nuisance_eval = nuisance_eval
    )
  }))
  reference <- pmm3_reference(cfg)
  verdict <- pmm3_verdict(results, reference, cfg)

  utils::write.csv(results, file.path("output", "tables", "pmm3_symmetric_probe_results.csv"), row.names = FALSE)
  utils::write.csv(thresholds, file.path("output", "tables", "pmm3_symmetric_probe_thresholds.csv"), row.names = FALSE)
  utils::write.csv(reference, file.path("output", "tables", "pmm3_symmetric_probe_reference.csv"), row.names = FALSE)
  utils::write.csv(verdict$checks, file.path("output", "tables", "pmm3_symmetric_probe_checks.csv"), row.names = FALSE)
  utils::write.csv(verdict$overall, file.path("output", "tables", "pmm3_symmetric_probe_verdict.csv"), row.names = FALSE)
  pmm3_write_report(results, thresholds, verdict, reference, cfg)
  utils::capture.output(
    sessionInfo(),
    file = file.path("output", "session_info", "session_info_pmm3_symmetric_probe.txt")
  )

  invisible(list(results = results, thresholds = thresholds, reference = reference,
                 verdict = verdict, config = cfg))
}

if (sys.nframe() == 0L) {
  out <- pmm3_main()
  if (!isTRUE(out$verdict$overall$pass[[1L]])) {
    quit(status = 1L)
  }
}

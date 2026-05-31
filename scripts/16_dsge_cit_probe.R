# DSGE-CIT feasibility probe.
#
# This is not an analytic Wiedermann-Shi replacement test. It is a gated
# reconstruction-error diagnostic for the same strict PHQ-like cell that blocked
# the GSA-LLR track:
#   N = 2000, R = 1.5, kappa3 = 0, kappa4 = -1.3.
#
# Detector: fit class-specific bidirectional reconstructions under H0 and H1,
# compute S = log(MSED_H0) - log(MSED_H1), calibrate threshold under held-out H0.

source(file.path("scripts", "13_gsa_llr_detector.R"))

dsge_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

dsge_arg_value <- function(args, name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", name, "="), "", hit[[1L]])
}

dsge_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  quick <- "--quick" %in% args || identical(Sys.getenv("DSGE_CIT_QUICK"), "1")
  list(
    quick = quick,
    N = 2000L,
    R_alt = 1.5,
    kappa3 = 0,
    kappa4 = -1.3,
    alpha = 0.05,
    train_n = if (quick) 8000L else 30000L,
    calibration_reps = if (quick) 160L else 600L,
    eval_reps = if (quick) 80L else 300L,
    ridge = 1e-6,
    seed_base = 20260527L,
    gsa_reference_power = 0.296,
    gaussian_nuisance_max = 0.10,
    bases = c("poly3", "poly4", "signed_frac", "robust")
  )
}

dsge_standardize <- function(x, label) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(sprintf("%s contains non-finite values", label))
  }
  mu <- mean(x)
  sig <- stats::sd(x)
  if (!is.finite(sig) || sig <= 0) {
    stop(sprintf("Cannot standardize degenerate vector: %s", label))
  }
  list(z = (x - mu) / sig, mean = mu, sd = sig)
}

dsge_basis <- function(z, basis = "poly3") {
  z <- as.numeric(z)
  if (any(!is.finite(z))) {
    stop("basis input contains non-finite values")
  }
  signed_pow <- function(p) sign(z) * abs(z)^p
  mag_pow <- function(p) abs(z)^p

  out <- switch(
    basis,
    poly3 = cbind(z = z, z2 = z^2, z3 = z^3),
    poly4 = cbind(z = z, z2 = z^2, z3 = z^3, z4 = z^4),
    signed_frac = cbind(
      sp05 = signed_pow(0.5),
      z = z,
      sp15 = signed_pow(1.5),
      z2 = z^2,
      sp25 = signed_pow(2.5),
      abs05 = mag_pow(0.5),
      abs15 = mag_pow(1.5)
    ),
    robust = cbind(
      z = z,
      tanh = tanh(z),
      sp05 = signed_pow(0.5),
      sp15 = signed_pow(1.5),
      logabs = log1p(abs(z)),
      z_clip = pmax(pmin(z, 2.5), -2.5),
      z2_clip = pmax(pmin(z, 2.5), -2.5)^2
    ),
    stop(sprintf("Unknown DSGE basis: %s", basis))
  )

  out <- as.matrix(out)
  storage.mode(out) <- "double"
  out
}

dsge_prepare_dat <- function(dat) {
  x1 <- dsge_standardize(dat$x1, "window_x1")$z
  x2 <- dsge_standardize(dat$x2, "window_x2")$z
  data.frame(x1 = x1, x2 = x2)
}

dsge_scale_basis <- function(phi) {
  center <- colMeans(phi)
  phi_c <- sweep(phi, 2L, center, "-")
  scale <- sqrt(colSums(phi_c^2) / max(1L, nrow(phi_c) - 1L))
  bad <- !is.finite(scale) | scale <= 0
  if (any(bad)) {
    keep <- !bad
    phi_c <- phi_c[, keep, drop = FALSE]
    center <- center[keep]
    scale <- scale[keep]
  }
  if (ncol(phi_c) == 0L) {
    stop("All DSGE basis columns are degenerate")
  }
  list(phi = sweep(phi_c, 2L, scale, "/"), center = center, scale = scale)
}

dsge_apply_basis_scale <- function(phi, center, scale) {
  common <- intersect(colnames(phi), names(center))
  if (length(common) != length(center)) {
    stop("Basis columns do not match fitted model")
  }
  phi <- phi[, names(center), drop = FALSE]
  sweep(sweep(phi, 2L, center, "-"), 2L, scale, "/")
}

dsge_fit_direction <- function(x, y, basis = "poly3", ridge = 1e-6) {
  xs <- dsge_standardize(x, "x")
  ys <- dsge_standardize(y, "y")
  raw_phi <- dsge_basis(xs$z, basis = basis)
  sc <- dsge_scale_basis(raw_phi)
  phi <- sc$phi

  F <- crossprod(phi) / max(1L, nrow(phi) - 1L)
  B <- as.numeric(crossprod(phi, ys$z) / max(1L, nrow(phi) - 1L))
  ridge_abs <- ridge * max(1, mean(diag(F)))
  K <- as.numeric(solve(F + ridge_abs * diag(ncol(F)), B))
  cond_F <- kappa(F + ridge_abs * diag(ncol(F)), exact = TRUE)

  list(
    basis = basis,
    x_mean = xs$mean,
    x_sd = xs$sd,
    y_mean = ys$mean,
    y_sd = ys$sd,
    phi_center = sc$center,
    phi_scale = sc$scale,
    K = K,
    F = F,
    B = B,
    ridge_abs = ridge_abs,
    cond_F = cond_F
  )
}

dsge_direction_mse <- function(model, x, y) {
  xs <- (as.numeric(x) - model$x_mean) / model$x_sd
  ys <- (as.numeric(y) - model$y_mean) / model$y_sd
  phi <- dsge_apply_basis_scale(
    dsge_basis(xs, basis = model$basis),
    center = model$phi_center,
    scale = model$phi_scale
  )
  pred <- as.numeric(phi %*% model$K)
  mean((ys - pred)^2)
}

dsge_fit_class <- function(dat, basis = "poly3", ridge = 1e-6) {
  forward <- dsge_fit_direction(dat$x1, dat$x2, basis = basis, ridge = ridge)
  backward <- dsge_fit_direction(dat$x2, dat$x1, basis = basis, ridge = ridge)
  list(
    basis = basis,
    forward = forward,
    backward = backward,
    cond_max = max(forward$cond_F, backward$cond_F)
  )
}

dsge_msed <- function(model, dat) {
  e12 <- dsge_direction_mse(model$forward, dat$x1, dat$x2)
  e21 <- dsge_direction_mse(model$backward, dat$x2, dat$x1)
  mean(c(e12, e21))
}

dsge_score <- function(fit, dat) {
  e0 <- dsge_msed(fit$h0, dat)
  e1 <- dsge_msed(fit$h1, dat)
  log(e0) - log(e1)
}

dsge_train_models <- function(cfg, basis, seed) {
  h0 <- gsa_llr_dgp(
    N = cfg$train_n, R = 1, kappa3 = cfg$kappa3, kappa4 = cfg$kappa4,
    seed = seed + 1L
  )
  h1 <- gsa_llr_dgp(
    N = cfg$train_n, R = cfg$R_alt, kappa3 = cfg$kappa3, kappa4 = cfg$kappa4,
    seed = seed + 2L
  )
  h0 <- dsge_prepare_dat(h0)
  h1 <- dsge_prepare_dat(h1)
  list(
    h0 = dsge_fit_class(h0, basis = basis, ridge = cfg$ridge),
    h1 = dsge_fit_class(h1, basis = basis, ridge = cfg$ridge)
  )
}

dsge_score_sample <- function(fit, reps, N, R, kappa3, kappa4, seed) {
  scores <- numeric(reps)
  for (r in seq_len(reps)) {
    dat <- gsa_llr_dgp(
      N = N, R = R, kappa3 = kappa3, kappa4 = kappa4,
      seed = seed + r
    )
    dat <- dsge_prepare_dat(dat)
    scores[[r]] <- dsge_score(fit, dat)
  }
  scores
}

dsge_cell <- function(cfg, basis, cell_id) {
  seed <- cfg$seed_base + 100000L * cell_id
  start <- proc.time()[["elapsed"]]
  fit <- dsge_train_models(cfg, basis = basis, seed = seed)

  cal_null <- dsge_score_sample(
    fit = fit, reps = cfg$calibration_reps, N = cfg$N,
    R = 1, kappa3 = cfg$kappa3, kappa4 = cfg$kappa4,
    seed = seed + 10000L
  )
  threshold <- as.numeric(stats::quantile(
    cal_null, probs = 1 - cfg$alpha, names = FALSE, na.rm = TRUE
  ))

  eval_null <- dsge_score_sample(
    fit = fit, reps = cfg$eval_reps, N = cfg$N,
    R = 1, kappa3 = cfg$kappa3, kappa4 = cfg$kappa4,
    seed = seed + 20000L
  )
  eval_alt <- dsge_score_sample(
    fit = fit, reps = cfg$eval_reps, N = cfg$N,
    R = cfg$R_alt, kappa3 = cfg$kappa3, kappa4 = cfg$kappa4,
    seed = seed + 30000L
  )
  nuisance <- dsge_score_sample(
    fit = fit, reps = cfg$eval_reps, N = cfg$N,
    R = cfg$R_alt, kappa3 = 0, kappa4 = 0,
    seed = seed + 40000L
  )

  elapsed <- proc.time()[["elapsed"]] - start
  data.frame(
    basis = basis,
    N = cfg$N,
    R_alt = cfg$R_alt,
    kappa3 = cfg$kappa3,
    kappa4 = cfg$kappa4,
    train_n = cfg$train_n,
    calibration_reps = cfg$calibration_reps,
    eval_reps = cfg$eval_reps,
    threshold_empirical = threshold,
    null_rate = mean(eval_null > threshold),
    alt_rate = mean(eval_alt > threshold),
    gaussian_nuisance_rate = mean(nuisance > threshold),
    mean_score_null = mean(eval_null),
    mean_score_alt = mean(eval_alt),
    mean_score_gaussian_nuisance = mean(nuisance),
    sd_score_null = stats::sd(eval_null),
    sd_score_alt = stats::sd(eval_alt),
    sd_score_gaussian_nuisance = stats::sd(nuisance),
    cond_h0 = fit$h0$cond_max,
    cond_h1 = fit$h1$cond_max,
    cond_max = max(fit$h0$cond_max, fit$h1$cond_max),
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

dsge_verdict <- function(results, cfg) {
  type1_ok <- is.finite(results$null_rate) &
    results$null_rate >= 0.035 & results$null_rate <= 0.075
  power_ok <- is.finite(results$alt_rate) & results$alt_rate > 0.40
  beats_gsa <- is.finite(results$alt_rate) & results$alt_rate > cfg$gsa_reference_power
  nuisance_ok <- is.finite(results$gaussian_nuisance_rate) &
    results$gaussian_nuisance_rate <= cfg$gaussian_nuisance_max
  cond_ok <- is.finite(results$cond_max) & results$cond_max < 1e6
  per_basis <- data.frame(
    basis = results$basis,
    type1_ok = type1_ok,
    phq_power_gt_0.40 = power_ok,
    beats_gsa_order4_reference = beats_gsa,
    gaussian_nuisance_ok = nuisance_ok,
    conditioning_ok = cond_ok,
    basis_pass = type1_ok & power_ok & beats_gsa & nuisance_ok & cond_ok,
    stringsAsFactors = FALSE
  )
  overall <- data.frame(
    criterion = c(
      "any_basis_passes_all_gates",
      "any_basis_phq_power_gt_0.40",
      "any_basis_beats_gsa_order4_reference",
      "any_basis_passes_gaussian_nuisance_guard",
      "all_conditioning_ok"
    ),
    pass = c(
      any(per_basis$basis_pass),
      any(per_basis$phq_power_gt_0.40),
      any(per_basis$beats_gsa_order4_reference),
      any(per_basis$gaussian_nuisance_ok),
      all(per_basis$conditioning_ok)
    ),
    stringsAsFactors = FALSE
  )
  list(per_basis = per_basis, overall = overall)
}

dsge_write_report <- function(results, verdict, cfg) {
  con <- file(file.path("output", "sanity_check_dsge_cit_probe.txt"), open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("DSGE-CIT feasibility probe", con)
  writeLines("================================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Quick mode: %s", cfg$quick), con)
  writeLines(sprintf(
    "Target: N=%d; R_alt=%.2f; kappa3=%.2f; kappa4=%.2f",
    cfg$N, cfg$R_alt, cfg$kappa3, cfg$kappa4
  ), con)
  writeLines(sprintf(
    "Train_n=%d; calibration_reps=%d; eval_reps=%d",
    cfg$train_n, cfg$calibration_reps, cfg$eval_reps
  ), con)
  writeLines("All windows are column-standardized before reconstruction.", con)
  writeLines("", con)

  writeLines("Verdict:", con)
  utils::write.table(verdict$overall, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Per-basis gates:", con)
  utils::write.table(verdict$per_basis, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Results:", con)
  show_cols <- c(
    "basis", "null_rate", "alt_rate", "gaussian_nuisance_rate",
    "threshold_empirical", "mean_score_null", "mean_score_alt",
    "mean_score_gaussian_nuisance", "cond_max", "elapsed_sec"
  )
  utils::write.table(results[, show_cols], con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  if (any(verdict$per_basis$basis_pass)) {
    writeLines("Interpretation:", con)
    writeLines("- At least one DSGE reconstruction basis passes the feasibility gates. Treat this as a candidate method requiring a larger calibration and ablation study, not as a finished manuscript claim.", con)
  } else {
    writeLines("Interpretation:", con)
    writeLines("- No DSGE reconstruction basis passes all feasibility gates in this probe. DSGE should remain exploratory unless a stronger basis-selection or CF-style modeling idea is introduced.", con)
  }
}

dsge_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  dsge_dirs()
  cfg <- dsge_config(args)
  results <- do.call(rbind, lapply(seq_along(cfg$bases), function(i) {
    dsge_cell(cfg, basis = cfg$bases[[i]], cell_id = i)
  }))
  verdict <- dsge_verdict(results, cfg)

  utils::write.csv(
    results,
    file.path("output", "tables", "dsge_cit_probe_results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    verdict$per_basis,
    file.path("output", "tables", "dsge_cit_probe_basis_verdict.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    verdict$overall,
    file.path("output", "tables", "dsge_cit_probe_verdict.csv"),
    row.names = FALSE
  )
  dsge_write_report(results, verdict, cfg)
  utils::capture.output(
    sessionInfo(),
    file = file.path("output", "session_info", "session_info_dsge_cit_probe.txt")
  )

  invisible(list(results = results, verdict = verdict, config = cfg))
}

if (sys.nframe() == 0L) {
  out <- dsge_main()
  if (!any(out$verdict$per_basis$basis_pass)) {
    quit(status = 1L)
  }
}

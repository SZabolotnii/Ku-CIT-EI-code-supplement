# DSGE hybrid CIT classifier.
#
# Positive-result candidate:
#   y=1: H1_ng, asymmetric loading with platykurtic confounder
#   y=0: H0_ng and H1_gauss_nuisance
#
# Feature sets:
#   existing  = cumulant / moment / GSA features
#   dsge      = reconstruction log-MSED features
#   hybrid    = existing + DSGE
#
# Thresholds are calibrated on validation negative windows only.

source(file.path("scripts", "02_naive_estimators.R"))
source(file.path("scripts", "13_gsa_llr_detector.R"))
source(file.path("scripts", "16_dsge_cit_probe.R"))

hybrid_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

hybrid_arg_value <- function(args, name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", name, ""), "", hit[[1L]])
}

hybrid_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  quick <- "--quick" %in% args || identical(Sys.getenv("DSGE_HYBRID_QUICK"), "1")
  list(
    quick = quick,
    N = 2000L,
    alpha = 0.05,
    train_per_scenario = if (quick) 80L else 400L,
    validation_per_scenario = if (quick) 80L else 400L,
    test_per_scenario = if (quick) 120L else 600L,
    dsge_train_n = if (quick) 8000L else 30000L,
    bootstrap_R = if (quick) 500L else 2000L,
    seed_base = 20260527L,
    ridge = 1e-6,
    bases = c("poly3", "poly4", "signed_frac", "robust")
  )
}

hybrid_scenarios <- function() {
  data.frame(
    scenario = c("H0_ng", "H1_gauss_nuisance", "H1_ng"),
    y = c(0L, 0L, 1L),
    R = c(1.0, 1.5, 1.5),
    kappa3 = c(0, 0, 0),
    kappa4 = c(-1.3, 0, -1.3),
    stringsAsFactors = FALSE
  )
}

hybrid_generate_window <- function(N, R, kappa3, kappa4, seed) {
  dsge_prepare_dat(gsa_llr_dgp(
    N = N, R = R, kappa3 = kappa3, kappa4 = kappa4, seed = seed
  ))
}

hybrid_skew <- function(x) {
  z <- (x - mean(x)) / stats::sd(x)
  mean(z^3)
}

hybrid_excess <- function(x) {
  z <- (x - mean(x)) / stats::sd(x)
  mean(z^4) - 3
}

hybrid_gsa_z4 <- function(dat) {
  out <- tryCatch(
    gsa_llr_score(
      dat$x1, dat$x2,
      R = 1.5, kappa3 = 0, kappa4 = -1.3,
      order = 4L, standardize = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(out) || !is.finite(out$z)) 0 else out$z
}

hybrid_existing_features <- function(dat) {
  dc3 <- delta_c3_naive(dat$x1, dat$x2)
  dc4 <- delta_c4_naive(dat$x1, dat$x2)
  sx1 <- stats::sd(dat$x1)
  sx2 <- stats::sd(dat$x2)
  c(
    cor_x1x2 = stats::cor(dat$x1, dat$x2),
    sd_ratio = sx2 / sx1,
    var_ratio = stats::var(dat$x2) / stats::var(dat$x1),
    delta_c3 = dc3,
    delta_c4 = dc4,
    abs_delta_c3 = abs(dc3),
    abs_delta_c4 = abs(dc4),
    skew_diff = hybrid_skew(dat$x2) - hybrid_skew(dat$x1),
    excess_diff = hybrid_excess(dat$x2) - hybrid_excess(dat$x1),
    abs_skew_diff = abs(hybrid_skew(dat$x2) - hybrid_skew(dat$x1)),
    abs_excess_diff = abs(hybrid_excess(dat$x2) - hybrid_excess(dat$x1)),
    gsa_z4_phq = hybrid_gsa_z4(dat),
    abs_gsa_z4_phq = abs(hybrid_gsa_z4(dat))
  )
}

hybrid_dsge_features <- function(dsge_fits, dat) {
  vals <- unlist(lapply(names(dsge_fits), function(basis) {
    fit <- dsge_fits[[basis]]
    e0 <- dsge_msed(fit$h0, dat)
    e1 <- dsge_msed(fit$h1, dat)
    c(
      log_msed_h0 = log(e0),
      log_msed_h1 = log(e1),
      score = log(e0) - log(e1)
    )
  }), use.names = FALSE)
  names(vals) <- as.vector(vapply(names(dsge_fits), function(basis) {
    paste0("dsge_", basis, "_", c("log_msed_h0", "log_msed_h1", "score"))
  }, character(3L)))
  vals
}

hybrid_fit_dsge_models <- function(cfg) {
  base_cfg <- dsge_config(if (cfg$quick) "--quick" else character())
  base_cfg$train_n <- cfg$dsge_train_n
  base_cfg$ridge <- cfg$ridge
  base_cfg$N <- cfg$N
  base_cfg$R_alt <- 1.5
  base_cfg$kappa3 <- 0
  base_cfg$kappa4 <- -1.3
  names(cfg$bases) <- cfg$bases
  lapply(seq_along(cfg$bases), function(i) {
    dsge_train_models(base_cfg, basis = cfg$bases[[i]],
                      seed = cfg$seed_base + 700000L + i * 1000L)
  }) |>
    stats::setNames(cfg$bases)
}

hybrid_build_dataset <- function(cfg, dsge_fits, split_name, per_scenario, seed_offset) {
  scenarios <- hybrid_scenarios()
  rows <- vector("list", nrow(scenarios) * per_scenario)
  meta <- vector("list", length(rows))
  idx <- 1L

  for (s in seq_len(nrow(scenarios))) {
    sc <- scenarios[s, ]
    for (r in seq_len(per_scenario)) {
      seed <- cfg$seed_base + seed_offset + s * 100000L + r
      dat <- hybrid_generate_window(
        N = cfg$N, R = sc$R, kappa3 = sc$kappa3,
        kappa4 = sc$kappa4, seed = seed
      )
      rows[[idx]] <- c(
        hybrid_existing_features(dat),
        hybrid_dsge_features(dsge_fits, dat)
      )
      meta[[idx]] <- data.frame(
        split = split_name,
        scenario = sc$scenario,
        y = sc$y,
        seed = seed,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  X <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
  X[] <- lapply(X, as.numeric)
  cbind(do.call(rbind, meta), X)
}

hybrid_feature_columns <- function(dat, feature_set) {
  dsge_cols <- grep("^dsge_", names(dat), value = TRUE)
  existing_cols <- setdiff(names(dat), c("split", "scenario", "y", "seed", dsge_cols))
  switch(
    feature_set,
    existing = existing_cols,
    dsge = dsge_cols,
    hybrid = c(existing_cols, dsge_cols),
    stop(sprintf("Unknown feature set: %s", feature_set))
  )
}

hybrid_fit_classifier <- function(train, feature_cols) {
  x <- train[, feature_cols, drop = FALSE]
  nzv <- vapply(x, function(col) stats::sd(col, na.rm = TRUE) > 1e-12, logical(1L))
  feature_cols <- feature_cols[nzv]
  form <- stats::as.formula(paste("y ~", paste(sprintf("`%s`", feature_cols), collapse = " + ")))
  fit <- suppressWarnings(stats::glm(form, data = train, family = stats::binomial()))
  list(fit = fit, feature_cols = feature_cols)
}

hybrid_predict_score <- function(model, dat) {
  p <- suppressWarnings(stats::predict(model$fit, newdata = dat, type = "response"))
  as.numeric(p)
}

hybrid_eval_feature_set <- function(feature_set, train, validation, test, cfg) {
  feature_cols <- hybrid_feature_columns(train, feature_set)
  model <- hybrid_fit_classifier(train, feature_cols)
  val_score <- hybrid_predict_score(model, validation)
  val_negative <- validation$y == 0L
  threshold <- as.numeric(stats::quantile(
    val_score[val_negative], 1 - cfg$alpha, names = FALSE, na.rm = TRUE
  ))
  test_score <- hybrid_predict_score(model, test)
  pred <- as.integer(test_score > threshold)

  scenario_rates <- aggregate(
    pred,
    by = list(scenario = test$scenario),
    FUN = mean
  )
  names(scenario_rates)[2] <- "rejection_rate"
  rate_of <- function(scenario) {
    hit <- scenario_rates$rejection_rate[scenario_rates$scenario == scenario]
    if (length(hit) == 0L) NA_real_ else hit[[1L]]
  }
  accuracy <- mean(pred == test$y)
  balanced_accuracy <- mean(c(
    mean(pred[test$y == 1L] == 1L),
    mean(pred[test$y == 0L] == 0L)
  ))

  list(
    method = feature_set,
    model = model,
    threshold = threshold,
    predictions = data.frame(
      method = feature_set,
      scenario = test$scenario,
      y = test$y,
      score = test_score,
      pred = pred,
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      method = feature_set,
      n_features = length(model$feature_cols),
      threshold = threshold,
      typeI_H0_ng = rate_of("H0_ng"),
      nuisance_H1_gauss = rate_of("H1_gauss_nuisance"),
      power_H1_ng = rate_of("H1_ng"),
      accuracy = accuracy,
      balanced_accuracy = balanced_accuracy,
      stringsAsFactors = FALSE
    )
  )
}

hybrid_bootstrap_compare <- function(pred_base, pred_candidate, candidate,
                                     R = 2000L, seed = 20260527L) {
  stopifnot(nrow(pred_base) == nrow(pred_candidate))
  y <- pred_base$y
  base_correct <- as.integer(pred_base$pred == y)
  candidate_correct <- as.integer(pred_candidate$pred == y)
  delta_obs <- mean(candidate_correct) - mean(base_correct)

  set.seed(seed)
  deltas <- replicate(R, {
    idx <- sample.int(length(y), length(y), replace = TRUE)
    mean(candidate_correct[idx]) - mean(base_correct[idx])
  })
  data.frame(
    candidate = candidate,
    baseline = "existing",
    metric = "accuracy",
    delta = delta_obs,
    ci_low = as.numeric(stats::quantile(deltas, 0.025, names = FALSE)),
    ci_high = as.numeric(stats::quantile(deltas, 0.975, names = FALSE)),
    p_value_one_sided = mean(deltas <= 0),
    R = R,
    stringsAsFactors = FALSE
  )
}

hybrid_verdict <- function(results, bootstrap) {
  existing <- results[results$method == "existing", ]
  candidates <- c("dsge", "hybrid")
  candidate_rows <- do.call(rbind, lapply(candidates, function(candidate) {
    row <- results[results$method == candidate, ]
    bs <- bootstrap[bootstrap$candidate == candidate, ]
    data.frame(
      candidate = candidate,
      typeI_H0_ng_ok = row$typeI_H0_ng <= 0.10,
      nuisance_H1_gauss_ok = row$nuisance_H1_gauss <= 0.10,
      power_H1_ng_ok = row$power_H1_ng > 0.40,
      power_beats_existing_ok = row$power_H1_ng - existing$power_H1_ng >= 0.05,
      accuracy_bootstrap_ok = bs$p_value_one_sided < 0.05 && bs$ci_low > 0,
      stringsAsFactors = FALSE
    )
  }))
  candidate_rows$candidate_pass <- apply(candidate_rows[, -1, drop = FALSE], 1L, all)

  overall <- data.frame(
    criterion = c("any_dsge_candidate_passes_all_gates"),
    pass = c(any(candidate_rows$candidate_pass)),
    stringsAsFactors = FALSE
  )
  list(candidates = candidate_rows, overall = overall)
}

hybrid_write_report <- function(results, bootstrap, verdict, cfg) {
  con <- file(file.path("output", "sanity_check_dsge_hybrid_cit.txt"), open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("DSGE hybrid CIT classifier", con)
  writeLines("================================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Quick mode: %s", cfg$quick), con)
  writeLines(sprintf(
    "N=%d; train/validation/test per scenario=%d/%d/%d; bootstrap_R=%d",
    cfg$N, cfg$train_per_scenario, cfg$validation_per_scenario,
    cfg$test_per_scenario, cfg$bootstrap_R
  ), con)
  writeLines("", con)

  writeLines("Method results:", con)
  utils::write.table(results, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Candidate vs existing bootstrap:", con)
  utils::write.table(bootstrap, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Candidate verdict:", con)
  utils::write.table(verdict$candidates, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)
  writeLines("Overall verdict:", con)
  utils::write.table(verdict$overall, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  if (any(verdict$candidates$candidate_pass)) {
    writeLines("Interpretation:", con)
    writeLines("- At least one DSGE candidate passes the nuisance-aware positive-result gates. This is a candidate positive result, still requiring wording as a calibrated simulation-family detector rather than a distribution-free analytic test.", con)
  } else {
    writeLines("Interpretation:", con)
    writeLines("- No DSGE candidate passes all positive-result gates. Do not use as manuscript-positive result without redesign.", con)
  }
}

hybrid_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  hybrid_dirs()
  cfg <- hybrid_config(args)
  dsge_fits <- hybrid_fit_dsge_models(cfg)

  train <- hybrid_build_dataset(cfg, dsge_fits, "train", cfg$train_per_scenario, 1000000L)
  validation <- hybrid_build_dataset(cfg, dsge_fits, "validation", cfg$validation_per_scenario, 2000000L)
  test <- hybrid_build_dataset(cfg, dsge_fits, "test", cfg$test_per_scenario, 3000000L)

  methods <- c("existing", "dsge", "hybrid")
  evals <- lapply(methods, hybrid_eval_feature_set, train = train, validation = validation, test = test, cfg = cfg)
  results <- do.call(rbind, lapply(evals, `[[`, "summary"))
  predictions <- do.call(rbind, lapply(evals, `[[`, "predictions"))

  pred_existing <- evals[[which(methods == "existing")]]$predictions
  pred_dsge <- evals[[which(methods == "dsge")]]$predictions
  pred_hybrid <- evals[[which(methods == "hybrid")]]$predictions
  bootstrap <- rbind(
    hybrid_bootstrap_compare(
      pred_base = pred_existing,
      pred_candidate = pred_dsge,
      candidate = "dsge",
      R = cfg$bootstrap_R,
      seed = cfg$seed_base + 9000000L
    ),
    hybrid_bootstrap_compare(
      pred_base = pred_existing,
      pred_candidate = pred_hybrid,
      candidate = "hybrid",
      R = cfg$bootstrap_R,
      seed = cfg$seed_base + 9100000L
    )
  )
  verdict <- hybrid_verdict(results, bootstrap)

  utils::write.csv(results, file.path("output", "tables", "dsge_hybrid_cit_results.csv"), row.names = FALSE)
  utils::write.csv(predictions, file.path("output", "tables", "dsge_hybrid_cit_predictions.csv"), row.names = FALSE)
  utils::write.csv(bootstrap, file.path("output", "tables", "dsge_hybrid_cit_bootstrap.csv"), row.names = FALSE)
  utils::write.csv(verdict$candidates, file.path("output", "tables", "dsge_hybrid_cit_candidate_verdict.csv"), row.names = FALSE)
  utils::write.csv(verdict$overall, file.path("output", "tables", "dsge_hybrid_cit_verdict.csv"), row.names = FALSE)
  hybrid_write_report(results, bootstrap, verdict, cfg)
  utils::capture.output(
    sessionInfo(),
    file = file.path("output", "session_info", "session_info_dsge_hybrid_cit.txt")
  )

  invisible(list(results = results, predictions = predictions, bootstrap = bootstrap, verdict = verdict, config = cfg))
}

if (sys.nframe() == 0L) {
  out <- hybrid_main()
  if (!any(out$verdict$candidates$candidate_pass)) {
    quit(status = 1L)
  }
}

# PATP alpha sweep for the DSGE hybrid positive-result track.
#
# This reuses the nuisance-aware hybrid design from script 17, but replaces the
# four discrete DSGE bases with a single PATP basis selected by validation
# performance. The key question is whether a continuous Kunchenko transition
# basis can preserve or improve the positive hybrid result with fewer DSGE
# features.

source(file.path("scripts", "17_dsge_hybrid_cit_classifier.R"))

patp_dirs <- function() {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "session_info"), recursive = TRUE, showWarnings = FALSE)
}

patp_arg_value <- function(args, name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", name, "="), "", hit[[1L]])
}

patp_parse_grid <- function(x) {
  if (length(x) != 1L || !nzchar(x)) return(NULL)
  vals <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  vals <- vals[is.finite(vals) & vals >= 0 & vals <= 1]
  sort(unique(vals))
}

patp_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- hybrid_config(args)
  quick <- cfg$quick
  grid_arg <- patp_arg_value(args, "alpha_grid", "")
  alpha_grid <- patp_parse_grid(grid_arg)
  if (is.null(alpha_grid) || length(alpha_grid) == 0L) {
    step <- as.numeric(patp_arg_value(args, "alpha_step", if (quick) "0.25" else "0.05"))
    alpha_grid <- seq(0, 1, by = step)
  }
  cfg$alpha_grid <- round(alpha_grid, 6)
  cfg$patp_n <- as.integer(patp_arg_value(args, "patp_n", "3"))
  cfg$ridge <- as.numeric(patp_arg_value(args, "ridge", "1e-4"))
  cfg$seed_base <- 20260527L + 1900L
  cfg$reference_discrete_file <- file.path("output", "tables", "dsge_hybrid_cit_results.csv")
  cfg
}

patp_power <- function(i, alpha) {
  1 / i + (4 - i - 3 / i) * alpha + (2 * i - 4 + 2 / i) * alpha^2
}

patp_basis_matrix <- function(z, alpha, n = 3L) {
  z <- as.numeric(z)
  signed_pow <- function(p) sign(z) * abs(z)^p
  raw_powers <- vapply(seq.int(2L, n + 1L), patp_power, numeric(1L), alpha = alpha)
  keep <- !duplicated(round(raw_powers, 6L)) & abs(raw_powers - 1) > 1e-6
  powers <- raw_powers[keep]
  if (length(powers) == 0L) {
    out <- cbind(z = z)
  } else {
    out <- cbind(z = z, do.call(cbind, lapply(seq_along(powers), function(j) {
      signed_pow(powers[[j]])
    })))
    colnames(out) <- c("z", sprintf("patp_p%.4f", powers))
  }
  storage.mode(out) <- "double"
  out
}

patp_basis_name <- function(alpha, n) {
  sprintf("patp_n%d_a%0.2f", n, alpha)
}

patp_parse_basis <- function(basis) {
  m <- regexec("^patp_n([0-9]+)_a([0-9]+\\.[0-9]+)$", basis)
  hit <- regmatches(basis, m)[[1L]]
  if (length(hit) != 3L) return(NULL)
  list(n = as.integer(hit[[2L]]), alpha = as.numeric(hit[[3L]]))
}

patp_install_dsge_basis <- function() {
  base_dsge_basis <- dsge_basis
  dsge_basis <<- function(z, basis = "poly3") {
    parsed <- patp_parse_basis(basis)
    if (!is.null(parsed)) {
      return(patp_basis_matrix(z, alpha = parsed$alpha, n = parsed$n))
    }
    base_dsge_basis(z, basis = basis)
  }
  invisible(TRUE)
}

patp_bootstrap_compare <- function(pred_base, pred_candidate, candidate,
                                   R = 2000L, seed = 20260527L) {
  hybrid_bootstrap_compare(
    pred_base = pred_base,
    pred_candidate = pred_candidate,
    candidate = candidate,
    R = R,
    seed = seed
  )
}

patp_alpha_cell <- function(alpha, cfg, train = NULL, validation = NULL, test = NULL) {
  basis <- patp_basis_name(alpha, cfg$patp_n)
  cfg_alpha <- cfg
  cfg_alpha$bases <- basis
  dsge_fits <- hybrid_fit_dsge_models(cfg_alpha)

  if (is.null(train)) {
    train <- hybrid_build_dataset(cfg_alpha, dsge_fits, "train", cfg_alpha$train_per_scenario, 1000000L)
  }
  if (is.null(validation)) {
    validation <- hybrid_build_dataset(cfg_alpha, dsge_fits, "validation", cfg_alpha$validation_per_scenario, 2000000L)
  }
  if (is.null(test)) {
    test <- hybrid_build_dataset(cfg_alpha, dsge_fits, "test", cfg_alpha$test_per_scenario, 3000000L)
  }

  methods <- c("existing", "dsge", "hybrid")
  evals <- lapply(methods, hybrid_eval_feature_set, train = train, validation = validation, test = test, cfg = cfg_alpha)
  summaries <- do.call(rbind, lapply(evals, `[[`, "summary"))
  summaries$alpha <- alpha
  summaries$patp_basis <- basis
  summaries$cond_h0 <- dsge_fits[[basis]]$h0$cond_max
  summaries$cond_h1 <- dsge_fits[[basis]]$h1$cond_max
  summaries$cond_max <- max(summaries$cond_h0, summaries$cond_h1)

  predictions <- do.call(rbind, lapply(evals, `[[`, "predictions"))
  predictions$alpha <- alpha
  predictions$patp_basis <- basis
  list(alpha = alpha, basis = basis, summaries = summaries, predictions = predictions,
       evals = evals, dsge_fits = dsge_fits)
}

patp_select_alpha <- function(sweep_results) {
  hybrid_rows <- sweep_results[sweep_results$method == "hybrid", ]
  guard <- hybrid_rows$typeI_H0_ng <= 0.10 &
    hybrid_rows$nuisance_H1_gauss <= 0.10 &
    hybrid_rows$cond_max < 1e6
  candidates <- hybrid_rows[guard, , drop = FALSE]
  if (nrow(candidates) == 0L) {
    candidates <- hybrid_rows
  }
  ord <- order(
    -candidates$balanced_accuracy,
    -candidates$power_H1_ng,
    candidates$nuisance_H1_gauss,
    candidates$typeI_H0_ng
  )
  candidates[ord[[1L]], , drop = FALSE]
}

patp_stability <- function(sweep_results, alpha_opt) {
  hybrid_rows <- sweep_results[sweep_results$method == "hybrid", ]
  neighbors <- hybrid_rows[abs(hybrid_rows$alpha - alpha_opt) <= 0.0500001, ]
  data.frame(
    alpha_opt = alpha_opt,
    n_neighbors_pm005 = nrow(neighbors),
    min_neighbor_balanced_accuracy = if (nrow(neighbors)) min(neighbors$balanced_accuracy) else NA_real_,
    max_neighbor_balanced_accuracy = if (nrow(neighbors)) max(neighbors$balanced_accuracy) else NA_real_,
    min_neighbor_power = if (nrow(neighbors)) min(neighbors$power_H1_ng) else NA_real_,
    max_neighbor_nuisance = if (nrow(neighbors)) max(neighbors$nuisance_H1_gauss) else NA_real_,
    stringsAsFactors = FALSE
  )
}

patp_discrete_reference <- function(cfg) {
  if (!file.exists(cfg$reference_discrete_file)) {
    return(data.frame(
      method = "hybrid",
      accuracy = NA_real_,
      balanced_accuracy = NA_real_,
      power_H1_ng = NA_real_,
      source = cfg$reference_discrete_file,
      stringsAsFactors = FALSE
    ))
  }
  ref <- utils::read.csv(cfg$reference_discrete_file)
  row <- ref[ref$method == "hybrid", , drop = FALSE]
  if (nrow(row) == 0L) {
    return(data.frame(
      method = "hybrid",
      accuracy = NA_real_,
      balanced_accuracy = NA_real_,
      power_H1_ng = NA_real_,
      source = cfg$reference_discrete_file,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    method = row$method,
    accuracy = row$accuracy,
    balanced_accuracy = row$balanced_accuracy,
    power_H1_ng = row$power_H1_ng,
    source = cfg$reference_discrete_file,
    stringsAsFactors = FALSE
  )
}

patp_verdict <- function(best_row, bootstrap, stability, discrete_reference, cfg) {
  bs <- bootstrap[bootstrap$candidate == "patp_hybrid", ]
  acc_ref <- discrete_reference$accuracy[[1L]]
  checks <- data.frame(
    criterion = c(
      "alpha_opt_not_boundary",
      "typeI_controlled",
      "gaussian_nuisance_guard",
      "power_gt_0.40",
      "bootstrap_accuracy_beats_existing",
      "conditioning_ok",
      "stable_neighborhood_pm005",
      "within_2pp_of_saved_discrete_hybrid_accuracy"
    ),
    value = c(
      best_row$alpha,
      best_row$typeI_H0_ng,
      best_row$nuisance_H1_gauss,
      best_row$power_H1_ng,
      bs$delta,
      best_row$cond_max,
      stability$min_neighbor_balanced_accuracy,
      if (is.finite(acc_ref)) best_row$accuracy - acc_ref else NA_real_
    ),
    pass = c(
      best_row$alpha > 0 && best_row$alpha < 1,
      best_row$typeI_H0_ng <= 0.10,
      best_row$nuisance_H1_gauss <= 0.10,
      best_row$power_H1_ng > 0.40,
      nrow(bs) == 1L && bs$p_value_one_sided < 0.05 && bs$ci_low > 0,
      best_row$cond_max < 1e6,
      is.finite(stability$min_neighbor_balanced_accuracy) &&
        stability$min_neighbor_balanced_accuracy >= best_row$balanced_accuracy - 0.03,
      !is.finite(acc_ref) || best_row$accuracy >= acc_ref - 0.02
    ),
    stringsAsFactors = FALSE
  )
  overall <- data.frame(
    criterion = "patp_hybrid_passes_positive_candidate_gates",
    pass = all(checks$pass),
    stringsAsFactors = FALSE
  )
  list(checks = checks, overall = overall)
}

patp_write_report <- function(sweep_results, best_row, bootstrap, stability,
                              discrete_reference, verdict, cfg) {
  con <- file(file.path("output", "sanity_check_patp_dsge_hybrid_sweep.txt"), open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("PATP DSGE hybrid alpha sweep", con)
  writeLines("================================================================", con)
  writeLines(sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("Quick mode: %s", cfg$quick), con)
  writeLines(sprintf(
    "N=%d; train/validation/test per scenario=%d/%d/%d; bootstrap_R=%d",
    cfg$N, cfg$train_per_scenario, cfg$validation_per_scenario,
    cfg$test_per_scenario, cfg$bootstrap_R
  ), con)
  writeLines(sprintf("PATP n=%d; alpha_grid=%s", cfg$patp_n, paste(cfg$alpha_grid, collapse = ",")), con)
  writeLines("", con)

  writeLines("Best PATP hybrid row:", con)
  utils::write.table(best_row, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Bootstrap vs existing at selected alpha:", con)
  utils::write.table(bootstrap, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Stability:", con)
  utils::write.table(stability, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Saved discrete-basis hybrid reference:", con)
  utils::write.table(discrete_reference, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Gate checks:", con)
  utils::write.table(verdict$checks, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("Overall verdict:", con)
  utils::write.table(verdict$overall, con, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("", con)

  writeLines("All sweep rows:", con)
  utils::write.table(sweep_results, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

patp_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  patp_dirs()
  cfg <- patp_config(args)
  patp_install_dsge_basis()

  cells <- lapply(seq_along(cfg$alpha_grid), function(i) {
    alpha <- cfg$alpha_grid[[i]]
    message(sprintf("PATP alpha %.2f (%d/%d)", alpha, i, length(cfg$alpha_grid)))
    patp_alpha_cell(alpha, cfg)
  })
  sweep_results <- do.call(rbind, lapply(cells, `[[`, "summaries"))
  predictions <- do.call(rbind, lapply(cells, `[[`, "predictions"))
  best_row <- patp_select_alpha(sweep_results)
  alpha_opt <- best_row$alpha[[1L]]
  selected <- cells[[which(abs(cfg$alpha_grid - alpha_opt) < 1e-9)[[1L]]]]

  methods <- c("existing", "dsge", "hybrid")
  pred_existing <- selected$evals[[which(methods == "existing")]]$predictions
  pred_dsge <- selected$evals[[which(methods == "dsge")]]$predictions
  pred_hybrid <- selected$evals[[which(methods == "hybrid")]]$predictions
  bootstrap <- rbind(
    patp_bootstrap_compare(
      pred_existing, pred_dsge, "patp_dsge",
      R = cfg$bootstrap_R, seed = cfg$seed_base + 9000000L
    ),
    patp_bootstrap_compare(
      pred_existing, pred_hybrid, "patp_hybrid",
      R = cfg$bootstrap_R, seed = cfg$seed_base + 9100000L
    )
  )
  stability <- patp_stability(sweep_results, alpha_opt)
  discrete_reference <- patp_discrete_reference(cfg)
  verdict <- patp_verdict(best_row, bootstrap, stability, discrete_reference, cfg)

  utils::write.csv(sweep_results, file.path("output", "tables", "patp_dsge_hybrid_sweep_results.csv"), row.names = FALSE)
  utils::write.csv(predictions, file.path("output", "tables", "patp_dsge_hybrid_sweep_predictions.csv"), row.names = FALSE)
  utils::write.csv(best_row, file.path("output", "tables", "patp_dsge_hybrid_sweep_best.csv"), row.names = FALSE)
  utils::write.csv(bootstrap, file.path("output", "tables", "patp_dsge_hybrid_sweep_bootstrap.csv"), row.names = FALSE)
  utils::write.csv(stability, file.path("output", "tables", "patp_dsge_hybrid_sweep_stability.csv"), row.names = FALSE)
  utils::write.csv(discrete_reference, file.path("output", "tables", "patp_dsge_hybrid_discrete_reference.csv"), row.names = FALSE)
  utils::write.csv(verdict$checks, file.path("output", "tables", "patp_dsge_hybrid_sweep_checks.csv"), row.names = FALSE)
  utils::write.csv(verdict$overall, file.path("output", "tables", "patp_dsge_hybrid_sweep_verdict.csv"), row.names = FALSE)
  patp_write_report(sweep_results, best_row, bootstrap, stability, discrete_reference, verdict, cfg)
  utils::capture.output(
    sessionInfo(),
    file = file.path("output", "session_info", "session_info_patp_dsge_hybrid_sweep.txt")
  )

  invisible(list(sweep_results = sweep_results, predictions = predictions,
                 best = best_row, bootstrap = bootstrap, stability = stability,
                 discrete_reference = discrete_reference, verdict = verdict,
                 config = cfg))
}

if (sys.nframe() == 0L) {
  out <- patp_main()
  if (!isTRUE(out$verdict$overall$pass[[1L]])) {
    quit(status = 1L)
  }
}

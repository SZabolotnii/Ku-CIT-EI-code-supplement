# Narrowed GSA-LLR detector API for Mini-paper 2.
#
# This file implements the revised H2 from:
#   03-mini-paper-2-GSA-LLR/derivation/H2_revision_2026-05-27.md
#
# Scope:
# - order-3 skewness-sensitive detector;
# - order-4 detector, including the pure-kurtosis direction;
# - empirical-joint omnibus rule over max(abs(z3), abs(z4)).
#
# It deliberately does not encode or test any uniform dominance claim.

options(stringsAsFactors = FALSE)

gsa_llr_match_order <- function(order, default = 4L) {
  if (length(order) != 1L) {
    order <- default
  }
  order <- as.integer(order)
  if (!order %in% c(3L, 4L)) {
    stop("order must be 3 or 4")
  }
  order
}

gsa_llr_match_calibration <- function(calibration_mode) {
  match.arg(calibration_mode, c("empirical_null", "normal_approx", "vp_bound"))
}

gsa_llr_F <- function(order = c(3L, 4L)) {
  order <- gsa_llr_match_order(order)
  F4 <- matrix(c(
     10 / 9,  -10 / 3,       0,        0,
    -10 / 3,  130 / 9,       0,        0,
          0,        0,  500 / 27, -400 / 9,
          0,        0, -400 / 9, 3200 / 27
  ), nrow = 4L, byrow = TRUE)
  idx <- seq_len(if (order == 3L) 2L else 4L)
  F4[idx, idx, drop = FALSE]
}

gsa_llr_Y <- function(R, kappa3, kappa4 = 0, order = c(3L, 4L)) {
  order <- gsa_llr_match_order(order)
  y4 <- c(
    R * (R - 1) * kappa3,
    (1 - R^3) * kappa3,
    R * (R^2 - 1) * kappa4,
    (1 - R^4) * kappa4
  )
  y4[seq_len(if (order == 3L) 2L else 4L)]
}

gsa_llr_standardize <- function(x, label) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(sprintf("%s contains non-finite values", label))
  }
  x <- x - mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) {
    stop(sprintf("Cannot standardize a degenerate score vector: %s", label))
  }
  x / s
}

gsa_llr_basis <- function(x1, x2, order = c(3L, 4L), standardize = TRUE) {
  order <- gsa_llr_match_order(order)
  if (length(x1) != length(x2) || length(x1) < 2L) {
    stop("x1 and x2 must have the same length >= 2")
  }
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  if (any(!is.finite(x1)) || any(!is.finite(x2))) {
    stop("x1 and x2 must contain only finite values")
  }

  if (standardize) {
    x1 <- gsa_llr_standardize(x1, "x1")
    x2 <- gsa_llr_standardize(x2, "x2")
  } else {
    x1 <- x1 - mean(x1)
    x2 <- x2 - mean(x2)
    if (stats::sd(x1) <= 0 || stats::sd(x2) <= 0) {
      stop("Degenerate centered score vector")
    }
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

gsa_llr_weights <- function(R, kappa3, kappa4 = 0, order = c(3L, 4L),
                            eps = sqrt(.Machine$double.eps)) {
  order <- gsa_llr_match_order(order)
  F <- gsa_llr_F(order)
  Y <- gsa_llr_Y(R = R, kappa3 = kappa3, kappa4 = kappa4, order = order)
  K <- as.numeric(solve(F, Y))
  J <- as.numeric(crossprod(Y, K))
  degenerate_direction <- !is.finite(J) || J <= eps
  list(
    K = K,
    J = J,
    F = F,
    Y = Y,
    order = order,
    degenerate_direction = degenerate_direction
  )
}

gsa_llr_score <- function(x1, x2, R, kappa3, kappa4 = 0,
                          order = c(3L, 4L), standardize = TRUE,
                          eps = sqrt(.Machine$double.eps)) {
  order <- gsa_llr_match_order(order)
  wk <- gsa_llr_weights(
    R = R, kappa3 = kappa3, kappa4 = kappa4,
    order = order, eps = eps
  )
  n <- length(x1)

  if (isTRUE(wk$degenerate_direction)) {
    return(list(
      lambda = 0,
      z = NA_real_,
      J = wk$J,
      n = n,
      K = wk$K,
      F = wk$F,
      Y = wk$Y,
      order = order,
      degenerate_direction = TRUE
    ))
  }

  phi_bar <- colMeans(gsa_llr_basis(
    x1 = x1, x2 = x2, order = order, standardize = standardize
  ))
  lambda <- sum(wk$K * phi_bar)
  z <- sqrt(n) * lambda / sqrt(wk$J)

  list(
    lambda = lambda,
    z = z,
    J = wk$J,
    n = n,
    K = wk$K,
    F = wk$F,
    Y = wk$Y,
    order = order,
    degenerate_direction = FALSE
  )
}

gsa_llr_threshold <- function(alpha = 0.05,
                              calibration_mode = c("normal_approx", "vp_bound")) {
  calibration_mode <- match.arg(calibration_mode)
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }
  switch(
    calibration_mode,
    normal_approx = stats::qnorm(1 - alpha / 2),
    vp_bound = sqrt(8 / (9 * alpha))
  )
}

gsa_llr_component <- function(n, kappa3 = 0, kappa4 = 0) {
  if (abs(kappa3) > 1e-12 && abs(kappa4) <= 1e-12) {
    shape <- 4 / (abs(kappa3)^2)
    x <- stats::rgamma(n, shape = shape, scale = 1)
    x <- (x - mean(x)) / stats::sd(x)
    if (kappa3 < 0) x <- -x
    return(x)
  }

  if (abs(kappa3) <= 1e-12 && kappa4 > 1e-12) {
    df <- 6 / kappa4 + 4
    x <- stats::rt(n, df = df) / sqrt(df / (df - 2))
    return((x - mean(x)) / stats::sd(x))
  }

  if (abs(kappa3) <= 1e-12 && kappa4 < -1e-12) {
    a <- (-6 / kappa4 - 3) / 2
    if (!is.finite(a) || a <= 0) {
      stop("Requested negative excess kurtosis is outside the beta(a,a) helper range")
    }
    x <- stats::rbeta(n, shape1 = a, shape2 = a)
    return((x - mean(x)) / stats::sd(x))
  }

  if (abs(kappa3) > 1e-12 && abs(kappa4) > 1e-12) {
    skew_part <- gsa_llr_component(n, kappa3 = kappa3, kappa4 = 0)
    kurt_part <- gsa_llr_component(n, kappa3 = 0, kappa4 = kappa4)
    w3 <- abs(kappa3)
    w4 <- sqrt(abs(kappa4))
    x <- w3 * skew_part + w4 * kurt_part
    return((x - mean(x)) / stats::sd(x))
  }

  stats::rnorm(n)
}

gsa_llr_dgp <- function(N, R = 1, kappa3 = 0, kappa4 = 0,
                        seed = NULL, dgp_family = "phase95_component") {
  if (dgp_family != "phase95_component") {
    stop("Only dgp_family='phase95_component' is currently implemented")
  }
  if (!is.null(seed)) set.seed(seed)
  T_score <- stats::rnorm(N)
  U <- gsa_llr_component(N, kappa3 = kappa3, kappa4 = kappa4)
  W1 <- stats::rnorm(N)
  W2 <- stats::rnorm(N)
  data.frame(
    x1 = T_score + U + W1,
    x2 = T_score + R * U + W2
  )
}

gsa_llr_z_sample <- function(reps, N, dgp_R, dgp_kappa3, dgp_kappa4,
                             test_R, test_kappa3, test_kappa4,
                             order = c(3L, 4L), standardize = TRUE,
                             seed = 20260527L,
                             dgp_family = "phase95_component") {
  order <- gsa_llr_match_order(order)
  z_vals <- rep(NA_real_, reps)
  degenerate <- rep(FALSE, reps)

  for (r in seq_len(reps)) {
    dat <- gsa_llr_dgp(
      N = N, R = dgp_R, kappa3 = dgp_kappa3, kappa4 = dgp_kappa4,
      seed = seed + r, dgp_family = dgp_family
    )
    st <- gsa_llr_score(
      dat$x1, dat$x2,
      R = test_R, kappa3 = test_kappa3, kappa4 = test_kappa4,
      order = order, standardize = standardize
    )
    z_vals[r] <- st$z
    degenerate[r] <- isTRUE(st$degenerate_direction)
  }

  list(z = z_vals, degenerate = degenerate)
}

gsa_llr_calibrate_empirical_null <- function(N, R, kappa3, kappa4 = 0,
                                             order = c(3L, 4L),
                                             alpha = 0.05,
                                             calibration_reps = 1000L,
                                             seed = 20260527L,
                                             standardize = TRUE,
                                             dgp_family = "phase95_component",
                                             return_z = FALSE) {
  order <- gsa_llr_match_order(order)
  if (calibration_reps < 10L) {
    stop("calibration_reps must be >= 10")
  }
  null <- gsa_llr_z_sample(
    reps = calibration_reps, N = N,
    dgp_R = 1, dgp_kappa3 = kappa3, dgp_kappa4 = kappa4,
    test_R = R, test_kappa3 = kappa3, test_kappa4 = kappa4,
    order = order, standardize = standardize, seed = seed,
    dgp_family = dgp_family
  )
  abs_z <- abs(null$z)
  if (all(!is.finite(abs_z))) {
    threshold <- NA_real_
    null_rejection_estimate <- NA_real_
  } else {
    threshold <- as.numeric(stats::quantile(
      abs_z, probs = 1 - alpha, na.rm = TRUE, names = FALSE
    ))
    null_rejection_estimate <- mean(abs_z > threshold, na.rm = TRUE)
  }

  out <- list(
    calibration_mode = "empirical_null",
    alpha = alpha,
    calibration_reps = calibration_reps,
    dgp_family = dgp_family,
    seed = seed,
    N = N,
    R = R,
    kappa3 = kappa3,
    kappa4 = kappa4,
    order = order,
    standardize = standardize,
    threshold = threshold,
    threshold_z = threshold,
    null_mean = mean(null$z, na.rm = TRUE),
    null_sd = stats::sd(null$z, na.rm = TRUE),
    null_mean_abs = mean(abs_z, na.rm = TRUE),
    null_rejection_estimate = null_rejection_estimate,
    null_degenerate_reps = sum(null$degenerate)
  )
  if (isTRUE(return_z)) {
    out$null_z <- null$z
  }
  out
}

gsa_llr_component_test <- function(x1, x2, R, kappa3, kappa4 = 0,
                                   order = c(3L, 4L),
                                   alpha = 0.05,
                                   calibration_mode = c("empirical_null", "normal_approx", "vp_bound"),
                                   calibration_reps = 1000L,
                                   seed = 20260527L,
                                   standardize = TRUE,
                                   dgp_family = "phase95_component") {
  order <- gsa_llr_match_order(order)
  calibration_mode <- gsa_llr_match_calibration(calibration_mode)
  score <- gsa_llr_score(
    x1, x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
    order = order, standardize = standardize
  )

  if (isTRUE(score$degenerate_direction)) {
    return(list(
      component = sprintf("order%d", order),
      score = score,
      calibration = NULL,
      calibration_mode = calibration_mode,
      threshold_z = NA_real_,
      reject = FALSE
    ))
  }

  if (calibration_mode == "empirical_null") {
    calibration <- gsa_llr_calibrate_empirical_null(
      N = length(x1), R = R, kappa3 = kappa3, kappa4 = kappa4,
      order = order, alpha = alpha, calibration_reps = calibration_reps,
      seed = seed, standardize = standardize, dgp_family = dgp_family
    )
    threshold_z <- calibration$threshold_z
  } else {
    calibration <- list(
      calibration_mode = calibration_mode,
      alpha = alpha,
      threshold_z = gsa_llr_threshold(alpha, calibration_mode),
      seed = NA_integer_
    )
    threshold_z <- calibration$threshold_z
  }

  list(
    component = sprintf("order%d", order),
    score = score,
    calibration = calibration,
    calibration_mode = calibration_mode,
    threshold_z = threshold_z,
    reject = is.finite(score$z) && is.finite(threshold_z) && abs(score$z) > threshold_z
  )
}

gsa_llr_test_order3 <- function(x1, x2, R, kappa3, kappa4 = 0,
                                alpha = 0.05,
                                calibration_mode = c("empirical_null", "normal_approx", "vp_bound"),
                                calibration_reps = 1000L,
                                seed = 20260527L,
                                standardize = TRUE,
                                dgp_family = "phase95_component") {
  gsa_llr_component_test(
    x1 = x1, x2 = x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
    order = 3L, alpha = alpha, calibration_mode = calibration_mode,
    calibration_reps = calibration_reps, seed = seed,
    standardize = standardize, dgp_family = dgp_family
  )
}

gsa_llr_test_order4 <- function(x1, x2, R, kappa3, kappa4 = 0,
                                alpha = 0.05,
                                calibration_mode = c("empirical_null", "normal_approx", "vp_bound"),
                                calibration_reps = 1000L,
                                seed = 20260527L,
                                standardize = TRUE,
                                dgp_family = "phase95_component") {
  gsa_llr_component_test(
    x1 = x1, x2 = x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
    order = 4L, alpha = alpha, calibration_mode = calibration_mode,
    calibration_reps = calibration_reps, seed = seed,
    standardize = standardize, dgp_family = dgp_family
  )
}

gsa_llr_row_max_abs <- function(z3, z4) {
  abs_mat <- cbind(abs(z3), abs(z4))
  apply(abs_mat, 1L, function(row) {
    if (all(!is.finite(row))) NA_real_ else max(row, na.rm = TRUE)
  })
}

gsa_llr_calibrate_omnibus_empirical_null <- function(N, R, kappa3, kappa4 = 0,
                                                     alpha = 0.05,
                                                     calibration_reps = 1000L,
                                                     seed = 20260527L,
                                                     standardize = TRUE,
                                                     dgp_family = "phase95_component",
                                                     return_z = FALSE) {
  z3 <- rep(NA_real_, calibration_reps)
  z4 <- rep(NA_real_, calibration_reps)
  deg3 <- rep(FALSE, calibration_reps)
  deg4 <- rep(FALSE, calibration_reps)

  for (r in seq_len(calibration_reps)) {
    dat <- gsa_llr_dgp(
      N = N, R = 1, kappa3 = kappa3, kappa4 = kappa4,
      seed = seed + r, dgp_family = dgp_family
    )
    st3 <- gsa_llr_score(
      dat$x1, dat$x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
      order = 3L, standardize = standardize
    )
    st4 <- gsa_llr_score(
      dat$x1, dat$x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
      order = 4L, standardize = standardize
    )
    z3[r] <- st3$z
    z4[r] <- st4$z
    deg3[r] <- isTRUE(st3$degenerate_direction)
    deg4[r] <- isTRUE(st4$degenerate_direction)
  }

  max_abs <- gsa_llr_row_max_abs(z3, z4)
  threshold <- if (all(!is.finite(max_abs))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(max_abs, probs = 1 - alpha, na.rm = TRUE, names = FALSE))
  }

  component_thresholds <- c(
    order3 = if (all(!is.finite(z3))) NA_real_ else as.numeric(stats::quantile(abs(z3), 1 - alpha, na.rm = TRUE, names = FALSE)),
    order4 = if (all(!is.finite(z4))) NA_real_ else as.numeric(stats::quantile(abs(z4), 1 - alpha, na.rm = TRUE, names = FALSE))
  )

  out <- list(
    calibration_mode = "empirical_joint_null",
    alpha = alpha,
    calibration_reps = calibration_reps,
    dgp_family = dgp_family,
    seed = seed,
    N = N,
    R = R,
    kappa3 = kappa3,
    kappa4 = kappa4,
    standardize = standardize,
    threshold = threshold,
    threshold_z = threshold,
    component_thresholds = component_thresholds,
    null_rejection_estimate = mean(max_abs > threshold, na.rm = TRUE),
    null_mean_abs_order3 = mean(abs(z3), na.rm = TRUE),
    null_mean_abs_order4 = mean(abs(z4), na.rm = TRUE),
    null_sd_order3 = stats::sd(z3, na.rm = TRUE),
    null_sd_order4 = stats::sd(z4, na.rm = TRUE),
    null_degenerate_order3 = sum(deg3),
    null_degenerate_order4 = sum(deg4)
  )
  if (isTRUE(return_z)) {
    out$null_z_order3 <- z3
    out$null_z_order4 <- z4
    out$null_max_abs <- max_abs
  }
  out
}

gsa_llr_omnibus <- function(x1, x2, R, kappa3, kappa4 = 0,
                            alpha = 0.05,
                            calibration_mode = c("empirical_null", "normal_approx", "vp_bound"),
                            calibration_reps = 1000L,
                            seed = 20260527L,
                            standardize = TRUE,
                            dgp_family = "phase95_component") {
  calibration_mode <- gsa_llr_match_calibration(calibration_mode)
  score3 <- gsa_llr_score(
    x1, x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
    order = 3L, standardize = standardize
  )
  score4 <- gsa_llr_score(
    x1, x2, R = R, kappa3 = kappa3, kappa4 = kappa4,
    order = 4L, standardize = standardize
  )

  if (calibration_mode == "empirical_null") {
    calibration <- gsa_llr_calibrate_omnibus_empirical_null(
      N = length(x1), R = R, kappa3 = kappa3, kappa4 = kappa4,
      alpha = alpha, calibration_reps = calibration_reps, seed = seed,
      standardize = standardize, dgp_family = dgp_family
    )
    threshold_z <- calibration$threshold_z
    component_thresholds <- calibration$component_thresholds
    combined_rule <- "max_abs_joint_empirical_null"
  } else {
    component_alpha <- alpha / 2
    threshold_z <- gsa_llr_threshold(component_alpha, calibration_mode)
    component_thresholds <- c(order3 = threshold_z, order4 = threshold_z)
    calibration <- list(
      calibration_mode = calibration_mode,
      alpha = alpha,
      component_alpha = component_alpha,
      threshold_z = threshold_z,
      note = "Bonferroni component alpha split for non-empirical calibration"
    )
    combined_rule <- "any_component_reject_bonferroni"
  }

  z3 <- score3$z
  z4 <- score4$z
  component_reject <- c(
    order3 = is.finite(z3) && is.finite(component_thresholds["order3"]) &&
      abs(z3) > component_thresholds["order3"],
    order4 = is.finite(z4) && is.finite(component_thresholds["order4"]) &&
      abs(z4) > component_thresholds["order4"]
  )
  max_abs_z <- gsa_llr_row_max_abs(z3, z4)
  omnibus_reject <- if (calibration_mode == "empirical_null") {
    is.finite(max_abs_z) && is.finite(threshold_z) && max_abs_z > threshold_z
  } else {
    any(component_reject)
  }

  list(
    order3 = list(score = score3, reject = unname(component_reject["order3"])),
    order4 = list(score = score4, reject = unname(component_reject["order4"])),
    component_reject = component_reject,
    omnibus = list(
      reject = omnibus_reject,
      max_abs_z = as.numeric(max_abs_z),
      threshold_z = threshold_z,
      rule = combined_rule
    ),
    calibration = calibration,
    calibration_mode = calibration_mode,
    alpha = alpha,
    dominance_claim = FALSE
  )
}

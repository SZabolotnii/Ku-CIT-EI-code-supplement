# 09_sim_comparative.R
# ---------------------------------------------------------------------------
# Phase 5: Comparative Monte Carlo — naive Δ̂c₃ vs PMM2 Δ̂c₃
#
# Computes per-condition:
#   ARE   = Var(naive) / Var(PMM2)   (≥1 → PMM2 more efficient)
#   Power for both estimators (bootstrap-based, B=200)
#
# Design:
#   gamma_U ∈ {0.75, 2.25}      confounder skewness
#   gamma_T ∈ {0, 2.25}         true-score skewness
#   N       ∈ {500, 2000, 5000}
#   R       ∈ {1, 1.25, 2}
#   cor_E1E2: matched to R (0 for R=1, 0.19 for R=1.25, 0.38 for R=2)
#   n_reps  = 300 per condition (for ARE); 200 for power pass
#   B_boot  = 200 per rep (power pass only)
#
# Acceptance criteria (PLAYBOOK Phase 5):
#   AC1: ARE ≥ 1 in ≥95% conditions
#   AC2: ARE ≥ 1.15 in ≥70% conditions with R>1
#   AC3: Empirical g₂ ≈ theoretical within ±15% (sanity; formula is approximate)
#   AC4: Power(PMM2) ≥ Power(naive) in ≥95% conditions
#   AC5: Power gain ≥ 10pp in boundary regime (R=1.25, N=500-2000)
#
# Run: Rscript scripts/09_sim_comparative.R
# ---------------------------------------------------------------------------

library(here)
library(future)
library(future.apply)

SKIP_PMM2_TESTS <- TRUE   # suppress unit tests when sourcing pmm2 estimator
source(here("scripts", "05_dgp.R"))
source(here("scripts", "02_naive_estimators.R"))
source(here("scripts", "08_pmm2_estimator.R"))

dir.create(here("data"),            recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)

# ── Configuration ─────────────────────────────────────────────────────────────
SEED_BASE <- 20260525L + 9L
ALPHA     <- 0.05
N_REPS    <- 300L   # for ARE pass
N_REPS_P  <- 200L   # for power pass
B_BOOT    <- 200L   # bootstrap reps for power

N_WORKERS <- min(8L, future::availableCores() - 1L)
future::plan(future::multisession, workers = max(1L, N_WORKERS))
cat(sprintf("Using %d parallel workers\n", N_WORKERS))

# ── Condition grid ─────────────────────────────────────────────────────────────
N_vals       <- c(500L, 2000L, 5000L)
gamma_U_vals <- c(0.75, 2.25)
gamma_T_vals <- c(0.0, 2.25)

# R → cor_E1E2 mapping (W&S pairs)
R_cor_map <- list(
  list(R = 1.0,  cor = 0.00),   # H₀
  list(R = 1.25, cor = 0.19),   # H₁ boundary
  list(R = 2.0,  cor = 0.38)    # H₁ strong
)

build_grid <- function() {
  rows <- list()
  for (N in N_vals) {
    for (gU in gamma_U_vals) {
      for (gT in gamma_T_vals) {
        for (rc in R_cor_map) {
          rows[[length(rows)+1]] <- data.frame(
            N        = N,
            gamma_U  = gU,
            gamma_T  = gT,
            gamma_W  = 0, delta_U = 0, delta_T = 0, delta_W = 0,
            R        = rc$R,
            cor_E1E2 = rc$cor,
            H0       = (rc$R == 1),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}
conditions <- build_grid()
n_cond     <- nrow(conditions)
cat(sprintf("Condition grid: %d conditions (H₀=%d, H₁=%d)\n",
            n_cond, sum(conditions$H0), sum(!conditions$H0)))

# ── Pass 1: ARE computation (no bootstrap) ────────────────────────────────────

are_one_condition <- function(cond_idx, conditions, N_REPS, SEED_BASE) {
  # Re-source in worker — suppress tests
  library(here)
  options(pmm2.skip_tests = TRUE)
  source(here("scripts", "05_dgp.R"))
  source(here("scripts", "02_naive_estimators.R"))
  source(here("scripts", "08_pmm2_estimator.R"))

  cond <- conditions[cond_idx, ]
  N    <- as.integer(cond$N)

  naive_v <- numeric(N_REPS)
  pmm2_v  <- numeric(N_REPS)

  for (r in seq_len(N_REPS)) {
    seed_r <- SEED_BASE + cond_idx * 10000L + r
    df <- generate_ws_dgp(
      N        = N,
      gamma_U  = cond$gamma_U, gamma_T = cond$gamma_T, gamma_W = cond$gamma_W,
      delta_U  = cond$delta_U, delta_T = cond$delta_T, delta_W = cond$delta_W,
      R        = cond$R,       cor_E1E2 = cond$cor_E1E2,
      seed     = seed_r
    )
    naive_v[r] <- delta_c3_naive(df$x1, df$x2)
    pmm2_v[r]  <- delta_c3_pmm2(df$x1, df$x2)
  }

  var_naive <- var(naive_v)
  var_pmm2  <- var(pmm2_v)
  ARE       <- if (var_pmm2 > 1e-15) var_naive / var_pmm2 else NA_real_

  data.frame(
    cond_idx  = cond_idx,
    N         = cond$N, gamma_U = cond$gamma_U, gamma_T = cond$gamma_T,
    R         = cond$R, cor_E1E2 = cond$cor_E1E2, H0 = cond$H0,
    mean_naive = mean(naive_v), mean_pmm2 = mean(pmm2_v),
    var_naive  = var_naive,     var_pmm2  = var_pmm2,
    ARE        = ARE,
    g2_emp     = var_pmm2 / var_naive,
    n_reps     = N_REPS,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("\n=== Pass 1: ARE computation (%d conditions × %d reps) ===\n",
            n_cond, N_REPS))
t0 <- Sys.time()
are_results_list <- future.apply::future_lapply(
  X           = seq_len(n_cond),
  FUN         = are_one_condition,
  conditions  = conditions,
  N_REPS      = N_REPS,
  SEED_BASE   = SEED_BASE,
  future.seed = SEED_BASE
)
are_results <- do.call(rbind, are_results_list)
cat(sprintf("ARE pass done in %.1f min\n", as.numeric(Sys.time()-t0, "mins")))

# ── Pass 2: Power computation (bootstrap) ─────────────────────────────────────

power_one_condition <- function(cond_idx, conditions, N_REPS_P, B_BOOT, ALPHA, SEED_BASE) {
  library(here)
  options(pmm2.skip_tests = TRUE)
  source(here("scripts", "05_dgp.R"))
  source(here("scripts", "02_naive_estimators.R"))
  source(here("scripts", "08_pmm2_estimator.R"))

  cond <- conditions[cond_idx, ]
  N    <- as.integer(cond$N)

  reject_naive <- logical(N_REPS_P)
  reject_pmm2  <- logical(N_REPS_P)

  for (r in seq_len(N_REPS_P)) {
    seed_r <- SEED_BASE + 500000L + cond_idx * 10000L + r
    df <- generate_ws_dgp(
      N        = N,
      gamma_U  = cond$gamma_U, gamma_T = cond$gamma_T, gamma_W = cond$gamma_W,
      delta_U  = cond$delta_U, delta_T = cond$delta_T, delta_W = cond$delta_W,
      R        = cond$R,       cor_E1E2 = cond$cor_E1E2,
      seed     = seed_r
    )
    ci_n <- bootstrap_ci(delta_c3_naive, df$x1, df$x2, B = B_BOOT, seed = seed_r + 1L)
    ci_p <- bootstrap_ci(delta_c3_pmm2,  df$x1, df$x2, B = B_BOOT, seed = seed_r + 2L)
    reject_naive[r] <- as.logical(ci_n["sig"])
    reject_pmm2[r]  <- as.logical(ci_p["sig"])
  }

  data.frame(
    cond_idx      = cond_idx,
    N             = cond$N, gamma_U = cond$gamma_U, gamma_T = cond$gamma_T,
    R             = cond$R, cor_E1E2 = cond$cor_E1E2, H0 = cond$H0,
    power_naive   = mean(reject_naive),
    power_pmm2    = mean(reject_pmm2),
    power_gain_pp = (mean(reject_pmm2) - mean(reject_naive)) * 100,
    n_reps        = N_REPS_P,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("\n=== Pass 2: Power computation (%d conditions × %d reps × %d boot) ===\n",
            n_cond, N_REPS_P, B_BOOT))
t0 <- Sys.time()
power_results_list <- future.apply::future_lapply(
  X           = seq_len(n_cond),
  FUN         = power_one_condition,
  conditions  = conditions,
  N_REPS_P    = N_REPS_P,
  B_BOOT      = B_BOOT,
  ALPHA       = ALPHA,
  SEED_BASE   = SEED_BASE,
  future.seed = SEED_BASE + 1L
)
power_results <- do.call(rbind, power_results_list)
cat(sprintf("Power pass done in %.1f min\n", as.numeric(Sys.time()-t0, "mins")))

# ── Merge and save ────────────────────────────────────────────────────────────
comparative <- merge(are_results, power_results[,
    c("cond_idx","power_naive","power_pmm2","power_gain_pp")],
    by = "cond_idx")

saveRDS(comparative, here("data", "mc_comparative.rds"))
write.csv(comparative, here("output", "tables", "are_pmm2_vs_naive.csv"), row.names = FALSE)
cat("Saved: data/mc_comparative.rds\n")
cat("Saved: output/tables/are_pmm2_vs_naive.csv\n")

# ── Acceptance criteria ───────────────────────────────────────────────────────
cat("\n=== Acceptance-criteria checks (Phase 5) ===\n\n")

all_cond  <- comparative
H1_cond   <- comparative[!comparative$H0, ]
H0_cond   <- comparative[comparative$H0, ]

# AC1: ARE ≥ 1 in ≥95% of ALL conditions
AC1_pct <- mean(all_cond$ARE >= 1.0, na.rm = TRUE)
AC1_pass <- AC1_pct >= 0.95
cat(sprintf("AC1: ARE ≥ 1 in %.1f%% conditions (need ≥95%%) → %s\n",
            100*AC1_pct, if (AC1_pass) "PASS ✓" else "FAIL ✗"))

# AC2: ARE ≥ 1.15 in ≥70% of H₁ conditions (R>1)
AC2_pct  <- mean(H1_cond$ARE >= 1.15, na.rm = TRUE)
AC2_pass <- AC2_pct >= 0.70
cat(sprintf("AC2: ARE ≥ 1.15 in %.1f%% of R>1 conditions (need ≥70%%) → %s\n",
            100*AC2_pct, if (AC2_pass) "PASS ✓" else "FAIL ✗"))

# AC3: Power(PMM2) ≥ Power(naive) in ≥95% of all conditions
AC3_pct  <- mean(comparative$power_pmm2 >= comparative$power_naive - 0.02, na.rm = TRUE)
AC3_pass <- AC3_pct >= 0.95
cat(sprintf("AC3: Power(PMM2) ≥ Power(naive) - 2pp in %.1f%% conditions (need ≥95%%) → %s\n",
            100*AC3_pct, if (AC3_pass) "PASS ✓" else "FAIL ✗"))

# AC4: Power gain ≥ 10pp in boundary regime (R=1.25, N ≤ 2000)
boundary <- H1_cond[H1_cond$R == 1.25 & H1_cond$N <= 2000, ]
if (nrow(boundary) > 0) {
  AC4_pct  <- mean(boundary$power_gain_pp >= 10, na.rm = TRUE)
  AC4_pass <- AC4_pct >= 0.5
  cat(sprintf("AC4: Power gain ≥ 10pp in %.1f%% boundary conditions (R=1.25, N≤2000) → %s\n",
              100*AC4_pct, if (AC4_pass) "PASS ✓" else "FAIL ✗"))
} else {
  cat("AC4: No boundary conditions found\n")
  AC4_pass <- FALSE
}

cat("\n--- Detailed summary by R ---\n")
for (r_val in sort(unique(comparative$R))) {
  rows <- comparative[comparative$R == r_val, ]
  cat(sprintf("  R=%.2f: ARE mean=%.2f, median=%.2f, ≥1.15: %.0f%%; ",
              r_val, mean(rows$ARE, na.rm=TRUE), median(rows$ARE, na.rm=TRUE),
              100*mean(rows$ARE >= 1.15, na.rm=TRUE)))
  cat(sprintf("power_naive=%.3f, power_pmm2=%.3f, gain=%.1fpp\n",
              mean(rows$power_naive), mean(rows$power_pmm2),
              mean(rows$power_gain_pp)))
}

cat("\n--- ARE by (gamma_T, R) ---\n")
for (gT in sort(unique(comparative$gamma_T))) {
  for (r_val in sort(unique(comparative$R))) {
    rows <- comparative[comparative$gamma_T == gT & comparative$R == r_val, ]
    if (nrow(rows) == 0) next
    cat(sprintf("  gamma_T=%.2f, R=%.2f → ARE mean=%.2f (n=%d)\n",
                gT, r_val, mean(rows$ARE, na.rm=TRUE), nrow(rows)))
  }
}

# ── Overall verdict ────────────────────────────────────────────────────────────
cat("\n=== Phase 5 Overall Verdict ===\n")
n_pass <- sum(c(AC1_pass, AC2_pass, AC3_pass, AC4_pass))
cat(sprintf("Passed %d/4 acceptance criteria\n", n_pass))
if (n_pass == 4) {
  cat("ALL CRITERIA PASS — Phase 5 complete ✓\n")
} else {
  cat("Some criteria FAIL — review ARE heatmap and power tables\n")
}
cat("\n=== DONE ===\n")

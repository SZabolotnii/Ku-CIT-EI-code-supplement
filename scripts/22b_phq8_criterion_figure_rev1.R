# 22b_phq8_criterion_figure_rev1.R
# ---------------------------------------------------------------------------
# Redraws fig_phq8_criterion for the SORT revision 1 (2026-09-17) from the
# per-split artifact output/tables/phq8_criterion.csv produced by
# scripts/22_phq8_criterion.R. The numbers are NOT recomputed here.
#
# Changes vs the figure in scripts/22 (which stays as it was for the submitted
# version):
#   * subtitle removed — its "fall toward y=0" wording was false for the two
#     splits (ids 2, 69) whose corrected statistic overshoots through zero, and
#     its "55%" (median |shift/naive|) differed from the manuscript's 54%
#     (1 - median |PMM2/naive|);
#   * y-axis label in the revision's notation: phi_2 = X1^3 - X2^3 (the code
#     column is still named pmm2_phi4);
#   * legend wording: "detection" is a rejection of H0, so the erased class
#     reads "Rejected by naive only";
#   * title shortened so it is not clipped at 7.5 in;
#   * PDF written with cairo_pdf (2026-09-24): the default pdf() device left
#     Helvetica and Symbol unembedded, which fails the embedded-font check.
# Output: output/figures/fig_phq8_criterion_rev1.{pdf,png}
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here); library(ggplot2) })

df <- read.csv(here("output", "tables", "phq8_criterion.csv"))
stopifnot(nrow(df) == 70L)

df$decision <- factor(
  ifelse(df$sig_naive == 1L & df$sig_pmm2 == 0L, "Rejected by naive only",
  ifelse(df$sig_naive == 1L & df$sig_pmm2 == 1L, "Rejected by both",
                                                  "Rejected by neither")),
  levels = c("Rejected by both", "Rejected by naive only", "Rejected by neither"))

n_overshoot <- sum(sign(df$pmm2_phi4) != sign(df$naive_phi1) |
                   abs(df$pmm2_phi4) > abs(df$naive_phi1))
cat(sprintf("splits: %d | rejected by naive %d | by PMM2 %d | naive only %d | overshoot %d\n",
            nrow(df), sum(df$sig_naive), sum(df$sig_pmm2),
            sum(df$decision == "Rejected by naive only"), n_overshoot))

lim <- max(abs(c(df$naive_phi1, df$pmm2_phi4)))
p <- ggplot(df, aes(x = naive_phi1, y = pmm2_phi4, color = decision, shape = decision)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = c("Rejected by both" = "#4D4D4D",
                                "Rejected by naive only" = "#D73027",
                                "Rejected by neither" = "#9ECAE1"),
                     name = NULL) +
  scale_shape_manual(values = c("Rejected by both" = 16,
                                "Rejected by naive only" = 17,
                                "Rejected by neither" = 1),
                     name = NULL) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    title = "PHQ-8 (BRFSS 2010), 70 half-splits",
    x = expression("Naive " * Delta * hat(c)[3] * "  (" * bar(varphi)[1] * ")"),
    y = expression("Single-auxiliary PMM2 " * Delta * hat(c)[3] * "  (" * bar(varphi)[1] + K^"*" * bar(varphi)[2] * ")")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(here("output", "figures", "fig_phq8_criterion_rev1.pdf"), p, width = 7.5, height = 6,
       device = cairo_pdf)
ggsave(here("output", "figures", "fig_phq8_criterion_rev1.png"), p, width = 7.5, height = 6, dpi = 150)
cat("written: output/figures/fig_phq8_criterion_rev1.{pdf,png}\n")

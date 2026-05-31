# Rebuild the PATP-DSGE hybrid alpha-sweep figure used in the manuscript.

fig_dirs <- function() {
  dir.create(file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)
}

read_patp_hybrid_curve <- function(path = file.path("output", "tables", "patp_dsge_hybrid_sweep_results.csv")) {
  dat <- utils::read.csv(path)
  dat <- dat[dat$method == "hybrid", ]
  dat[order(dat$alpha), ]
}

draw_patp_hybrid_curve <- function(dat) {
  layout(matrix(c(1, 2), nrow = 2), heights = c(4.4, 1.05))
  on.exit({
    layout(1)
  }, add = TRUE)

  cols <- c(
    accuracy = "#0072B2",
    balanced = "#009E73",
    power = "#D55E00",
    nuisance = "#CC79A7",
    typeI = "#4D4D4D"
  )

  par(mar = c(4.3, 4.8, 1.0, 1.2), las = 1)
  plot(
    dat$alpha, dat$accuracy,
    type = "o", pch = 16, lwd = 2.4, col = cols[["accuracy"]],
    ylim = c(0, 1.02), xlim = c(0, 1),
    xlab = expression(PATP~alpha), ylab = "Rate / accuracy",
    cex.lab = 1.12, cex.axis = 0.95
  )
  lines(dat$alpha, dat$balanced_accuracy, type = "o", pch = 17, lwd = 2.2, col = cols[["balanced"]])
  lines(dat$alpha, dat$power_H1_ng, type = "o", pch = 15, lwd = 2.2, col = cols[["power"]])
  lines(dat$alpha, dat$nuisance_H1_gauss, type = "o", pch = 1, lwd = 2.0, col = cols[["nuisance"]])
  lines(dat$alpha, dat$typeI_H0_ng, type = "o", pch = 2, lwd = 2.0, col = cols[["typeI"]])

  abline(h = 0.10, lty = 3, lwd = 1.3, col = "grey45")
  abline(v = 0.75, lty = 2, lwd = 1.2, col = "grey35")
  text(0.755, 0.99, expression(alpha == 0.75), pos = 4, cex = 0.95, col = "grey20")
  text(0.98, 0.118, "0.10 guard", pos = 2, cex = 0.78, col = "grey35")
  grid(nx = NA, ny = NULL, lty = 3, col = "grey88")
  box()

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    horiz = FALSE,
    ncol = 3,
    bty = "n",
    cex = 0.86,
    x.intersp = 0.75,
    y.intersp = 1.05,
    legend = c(
      "Accuracy",
      "Balanced accuracy",
      "Power: H1 non-Gaussian",
      "Gaussian nuisance",
      "Type-I: H0 non-Gaussian"
    ),
    col = unname(cols),
    lty = 1,
    lwd = c(2.4, 2.2, 2.2, 2.0, 2.0),
    pch = c(16, 17, 15, 1, 2)
  )
}

write_patp_hybrid_figure <- function() {
  fig_dirs()
  dat <- read_patp_hybrid_curve()

  grDevices::pdf(file.path("output", "figures", "patp_dsge_hybrid_alpha_sweep.pdf"),
                 width = 7.4, height = 5.4, onefile = FALSE)
  draw_patp_hybrid_curve(dat)
  grDevices::dev.off()

  grDevices::png(file.path("output", "figures", "patp_dsge_hybrid_alpha_sweep.png"),
                 width = 1600, height = 1160, res = 220)
  draw_patp_hybrid_curve(dat)
  grDevices::dev.off()
}

if (identical(environment(), globalenv())) {
  write_patp_hybrid_figure()
}

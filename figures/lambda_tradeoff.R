# Accepted six-configuration A3 plot: replay validated numerical values only.
source("R/load.R"); load_mvgc()
source("R/utils/reproduction.R")
require_packages(c("ggplot2", "patchwork", "ragg"))
d <- read.csv(reproduction_input("lambda_tradeoff"))
stopifnot(identical(d$lambda, c(0, .025, .125, .25, .5, 1)),
          all(is.finite(d$mean_hellinger)), all(is.finite(d$fragmentation_pct)),
          d$lambda[which.min(d$mean_hellinger)] == .025)
main <- which(d$lambda == .125)
out <- mvgc_output_dir("figures/spatial_regularisation")
library(ggplot2); library(patchwork)
# Scale the reference's 29 cm export at 0.75*17.7 cm manuscript inclusion.
# Draw this review at prospective full text width; Figure A uses the same
# physical type/line/marker sizes at its existing 0.85*textwidth inclusion.
reference_scale <- 1 # Physical publication-size typography, shared with A1
blue <- "#1f77b4"; red <- "#d62728"
common <- theme_minimal(base_size = 8.5 * reference_scale, base_family = "sans") +
  theme(panel.grid.major = element_line(linewidth = .35 * reference_scale, colour = "#EBEBEB"),
    panel.grid.minor = element_line(linewidth = .2 * reference_scale, colour = "#F2F2F2"),
    axis.title = element_text(size = 8.5 * reference_scale, colour = "#202020"),
    axis.text = element_text(size = 7.5 * reference_scale, colour = "#202020"),
    plot.title = element_text(size = 9 * reference_scale, face = "plain", colour = "#202020", margin = margin(b = 5)),
    plot.margin = margin(5, 6, 4, 5), axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)))
x_axis <- scale_x_continuous(breaks = c(0, .125, .25, .5, .75, 1), labels = c("0", "0.125", "0.25", "0.5", "0.75", "1"),
                             minor_breaks = c(.0625, .1875, .375, .625, .875),
                             limits = c(-.025, 1.025), expand = expansion(mult = 0))
layers <- list(geom_vline(xintercept = .125, linetype = "dotted", colour = "grey55", linewidth = .3 * reference_scale),
  geom_line(colour = blue, linewidth = .45 * reference_scale),
  geom_point(colour = blue, fill = "white", shape = 21, size = 2.1 * reference_scale, stroke = .55 * reference_scale),
  geom_point(data = d[main, ], shape = 21, colour = red, fill = red, size = 2.6 * reference_scale, stroke = .5 * reference_scale),
  x_axis, common)
# D-MMM (0.2479196869) stays documented in the unchanged provenance note.
# Show absolute GC score variation, without rescaling the scores themselves.
pa <- ggplot(d, aes(lambda, mean_hellinger)) +
  layers + labs(title = "(a) Reconstruction accuracy", x = expression(paste("Smoothness parameter (", lambda, ")")),
                y = "Area-weighted Hellinger distance") +
  scale_y_continuous(breaks = seq(.185, .210, .005), minor_breaks = seq(.1875, .2075, .005),
                     limits = c(.1835, .2115), labels = function(x) sprintf("%.3f", x), expand = expansion(mult = 0))
pb <- ggplot(d, aes(lambda, fragmentation_pct)) + layers +
  labs(title = "(b) Label-map fragmentation", x = expression(paste("Smoothness parameter (", lambda, ")")), y = "Fragmentation (%)") +
  scale_y_continuous(breaks = seq(0, 30, 10), minor_breaks = seq(5, 35, 10), limits = c(0, 34), expand = expansion(mult = 0))
p <- pa + pb + plot_layout(widths = c(1, 1))
for (panel in list(pa, pb)) {
  built <- ggplot_build(panel)$data
  y <- if (identical(panel$labels$y, "Area-weighted Hellinger distance")) d$mean_hellinger else d$fragmentation_pct
  stopifnot(identical(built[[3]]$x, d$lambda), identical(built[[3]]$y, y),
            identical(built[[4]]$x, d$lambda[main]), identical(built[[4]]$y, y[main]),
            identical(built[[4]]$fill, red))
}
ggsave(file.path(out, "lambda_tradeoff_sensitivity.pdf"), p, device = cairo_pdf, width = 17.7, height = 7.5, units = "cm", bg = "white")
ggsave(file.path(out, "lambda_tradeoff_sensitivity.png"), p, device = ragg::agg_png, width = 17.7, height = 7.5, units = "cm", dpi = 450, bg = "white")

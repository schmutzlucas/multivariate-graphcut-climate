# Fig. 5 voting variants: connected observations and simple boxplots.
# Run from repository root. Both variants use the identical reviewed dataset.
# PME_REFERENCE_DATA_FILE optionally overrides the reviewed minimal plotting CSV.
# PME_FINAL_BOXPLOT_OUTPUT_DIR optionally overrides the output directory.
source("R/load.R")
load_mvgc()
source("R/utils/reproduction.R")
if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
library(ggplot2)
input <- reproduction_input("pme_boxplot")
out <- Sys.getenv("PME_FINAL_BOXPLOT_OUTPUT_DIR", mvgc_output_dir("figures"))
df <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("reference", "period_id", "period", "metric", "metric_label", "ell_LDR",
              "H_D_MMM", "H_GC", "relative_improvement_pct")
periods <- c("1976_2000", "2001_2025", "2026_2050", "2051_2075", "2076_2100")
period_labels <- c("1976\u20132000", "2001\u20132025", "2026\u20132050",
                   "2051\u20132075", "2076\u20132100")
stopifnot(identical(names(df), required), nrow(df) == 150L,
          setequal(df$metric, c("full", "0.01")), setequal(df$period_id, periods),
          length(unique(df$reference)) == 15L,
          !anyDuplicated(df[c("reference", "period_id", "metric")]),
          all(is.finite(df$relative_improvement_pct)), all(is.finite(df$H_D_MMM)),
          all(is.finite(df$H_GC)), all(df$H_D_MMM > 0),
          all(is.na(df$ell_LDR[df$metric == "full"])),
          all(df$ell_LDR[df$metric == "0.01"] == 0.01),
          all(table(df$period_id, df$metric) == 15L),
          max(abs(df$relative_improvement_pct -
            100 * (df$H_D_MMM - df$H_GC) / df$H_D_MMM)) < 1e-10)

refs <- unique(df$reference)
df$metric <- factor(df$metric, levels = c("full", "0.01"))
df$x <- match(df$period_id, periods) + ifelse(df$metric == "full", -.18, .18)
df$x_point <- df$x + (((match(df$reference, refs) * 7L) %% 15L) - 7) / 7 * .065
mean_df <- aggregate(relative_improvement_pct ~ metric + x, df, mean)

# Match the final PME overlays' 13-pt sans type, pale grids, white background,
# default 5.5-pt margins and 29 x 18 cm format. Green and violet distinguish
# diagnostics; points and trajectories exactly match the box outline colours.
fills <- c(full = "#A6CEA8", "0.01" = "#C2A5CF")
outlines <- c(full = "#2F6B43", "0.01" = "#75418C")
labels <- expression("Hellinger distance", "LDR Hellinger score, " * "\u2113"[plain(LDR)] == 0.01)
p <- ggplot(df, aes(x, relative_improvement_pct,
                    group = interaction(period_id, metric))) +
  geom_hline(yintercept = 0, colour = "#6A6A6A", linewidth = .65, linetype = "dashed") +
  geom_line(aes(x = x_point, colour = metric, group = interaction(reference, metric)),
             linewidth = .30, alpha = .24, show.legend = FALSE) +
  geom_boxplot(aes(fill = metric, colour = metric), width = .265,
               linewidth = .65, outlier.shape = NA, coef = 1.5) +
  geom_point(aes(x = x_point, colour = metric), size = 1.6, alpha = .68,
             show.legend = FALSE) +
  geom_point(data = mean_df, aes(x, relative_improvement_pct), inherit.aes = FALSE,
             shape = 23, fill = "white", colour = "#303030", size = 2.6, stroke = .65) +
  scale_fill_manual(values = fills, labels = labels, name = NULL, drop = FALSE) +
  scale_colour_manual(values = outlines, labels = labels, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = 1:5, labels = period_labels, limits = c(.5, 5.5),
                     minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(0, 60, 10), minor_breaks = seq(5, 55, 10),
                     limits = c(-4, 70), expand = expansion(mult = c(0, .02))) +
  labs(x = "Target period", y = "Relative improvement (%)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(linewidth = .55, colour = "#EBEBEB"),
    panel.grid.minor = element_line(linewidth = .35, colour = "#EBEBEB"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title = element_text(size = 13),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "inside", legend.position.inside = c(.97, .96),
    legend.justification = c(1, 1), legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", colour = "#6A6A6A", linewidth = .4),
    legend.margin = margin(5, 7, 5, 7),
    legend.text = element_text(size = 12.2, colour = "#303030"),
    legend.key.width = grid::unit(.6, "cm"),
    legend.key.height = grid::unit(.55, "cm"),
    legend.spacing.x = grid::unit(.15, "cm")) +
  guides(fill = guide_legend(nrow = 1), colour = guide_legend(nrow = 1))

# Verify the plotted boxes themselves, not just the source population counts.
built <- ggplot_build(p)
boxes <- built$data[[3]]
lines <- built$data[[2]]
stopifnot(nrow(boxes) == 10L, nrow(built$data[[4]]) == 150L,
          nrow(built$data[[5]]) == 10L, nrow(lines) == 150L,
          length(unique(lines$group)) == 30L, all(table(lines$group) == 5L))
for (g in unique(lines$group)) {
  z <- lines[lines$group == g, ]
  stopifnot(all(diff(z$x) > 0), length(unique(z$colour)) == 1L)
}
stopifnot(setequal(lines$colour, outlines),
          setequal(built$data[[4]]$colour, outlines), setequal(boxes$colour, outlines))
for (i in seq_len(nrow(boxes))) {
  values <- df$relative_improvement_pct[abs(df$x - boxes$x[i]) < 1e-12]
  quartiles <- quantile(values, c(.25, .5, .75), type = 7, names = FALSE)
  iqr <- quartiles[3] - quartiles[1]
  whisker_values <- values[values >= quartiles[1] - 1.5 * iqr &
                           values <= quartiles[3] + 1.5 * iqr]
  stopifnot(length(values) == 15L,
    max(abs(c(boxes$lower[i], boxes$middle[i], boxes$upper[i]) - quartiles)) < 1e-12,
    max(abs(c(boxes$ymin[i], boxes$ymax[i]) - range(whisker_values))) < 1e-12)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base <- file.path(out, Sys.getenv("PME_FINAL_BOXPLOT_BASENAME", "pme_reference_improvement"))
ggsave(paste0(base, ".png"), p, width = 29, height = 18, units = "cm", dpi = 300, bg = "white")
ggsave(paste0(base, ".pdf"), p, width = 29, height = 18, units = "cm", device = cairo_pdf, bg = "white")
# Byte-for-byte preservation of the reviewed plotting dataset.
stopifnot(file.copy(input, paste0(base, "_data.csv"), overwrite = TRUE))
caption <- paste0(
 "Reference-wise relative improvement of GC-0.125 over D-MMM in the perfect-model experiment ",
 "for Hellinger distance (green) and LDR Hellinger score (violet; $\\ell_{\\mathrm{LDR}}=0.01$), computed as ",
 "$100(H_{\\mathrm{D\\mbox{-}MMM}}-H_{\\mathrm{GC}})/H_{\\mathrm{D\\mbox{-}MMM}}$ ",
 "from reference-wise area-weighted distances. Positive values favour GC-0.125. ",
 "Each box summarizes 15 equally weighted held-out reference models: the central line is the median, ",
 "the box spans the interquartile range, and whiskers extend to the most distant values within ",
 "1.5 interquartile ranges. Points show individual references; white diamonds mark arithmetic means. ",
 "Faint lines connect each held-out reference across periods separately for each diagnostic, ",
 "using the same colours as the corresponding boxplot outlines and points.")
writeLines(c("% Final Fig. 5 caption proposal", paste0("\\caption{", caption, "}")),
           paste0(base, "_caption.tex"))

cat('Saved accepted connected reference-improvement figure to', base, '\n')

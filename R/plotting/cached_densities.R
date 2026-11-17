# Final overlay layout, separated from density estimation for Level 1 reproduction.
#' Replay a five-period PME density overlay from archived curve coordinates
#' @param recipe Numeric plot-recipe list with coordinates (x,y,group,colour),
#'   five-row annotations (period_label,mmm,gc,mmm_color,gc_color,y) and a finite
#'   positive y_limit. No ggplot object or full experiment workspace is needed.
#' @param metric_type Either full for ordinary Hellinger or partial for LDR.
#' @param tau_key Numeric-character ell_LDR threshold in (0,1) for partial; NULL
#'   for full. Formatting is preserved in the plotmath table heading.
#' @param inset_layout Standard main-text or enlarged appendix inset typography.
#' @return ggplot object with the original curves, table, colours and layout.
#'   Does not re-estimate densities, average references, save files or run graph cut.
plot_cached_overlay <- function(recipe, metric_type = "full", tau_key = NULL, inset_layout = c("standard", "appendix_enlarged")) {
  if (!is.character(metric_type) || length(metric_type) != 1L || is.na(metric_type) ||
      !metric_type %in% c("full", "partial")) stop("metric_type must be full or partial.")
  if (metric_type == "partial") {
    ell_ldr <- suppressWarnings(as.numeric(tau_key))
    if (length(ell_ldr) != 1L || !is.finite(ell_ldr) || ell_ldr <= 0 || ell_ldr >= 1)
      stop("Partial overlays require one finite tau_key in (0,1).")
  }
  curves <- recipe$coordinates
  table_df <- recipe$annotations
  common_y_limit <- recipe$y_limit
  if (!is.data.frame(curves) || !nrow(curves) ||
      !all(c("x", "y", "group", "colour") %in% names(curves)) ||
      !is.numeric(curves$x) || !is.numeric(curves$y) || any(!is.finite(curves$x)) ||
      any(!is.finite(curves$y)) || anyNA(curves$group) || anyNA(curves$colour) ||
      !is.data.frame(table_df) || nrow(table_df) != 5L ||
      !all(c("period_label", "mmm", "gc", "mmm_color", "gc_color", "y") %in% names(table_df)) ||
      anyNA(table_df$period_label) || anyDuplicated(table_df$period_label) ||
      !all(vapply(table_df[c("mmm", "gc", "y")], function(x) is.numeric(x) && all(is.finite(x)), logical(1))) ||
      !is.numeric(common_y_limit) || length(common_y_limit) != 1L ||
      !is.finite(common_y_limit) || common_y_limit <= 0) {
    stop("Expected finite cached density coordinates, five distinct period summaries and a positive y_limit.")
  }
  gc_label <- "GC-0.125"
  inset_layout <- match.arg(inset_layout)
  enlarged <- inset_layout == "appendix_enlarged"
  left <- if(enlarged) .45 else .57
  right <- .985; lower <- if(enlarged) .25 else .405; upper <- .975
  col_x <- if(enlarged) left + (right-left)*c(.05,.60,.88) else c(.595,.81,.935)
  title_y <- if(enlarged) upper-(upper-lower)*.105 else .915
  header_y <- if(enlarged) upper-(upper-lower)*.23 else .85
  separator_y <- if(enlarged) upper-(upper-lower)*.29 else .825
  if(enlarged) table_df$y <- common_y_limit*(upper-(upper-lower)*seq(.38,.88,length.out=nrow(table_df)))
  table_box <- data.frame(
    xmin = left, xmax = right,
    ymin = common_y_limit * lower, ymax = common_y_limit * upper
  )

  density_panel <- ggplot2::ggplot(
    curves,
    ggplot2::aes(x = x, y = y, group = group, color = colour)
  ) +
    ggplot2::geom_line(linewidth = 0.85, alpha = 1) +
    ggplot2::geom_rect(
      data = table_box,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white",
      color = "#6A6A6A",
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "text", x = col_x[1], y = common_y_limit * title_y,
      label = if (metric_type == "full") "bolditalic('Area-weighted H')" else
        paste0("bolditalic('Area-weighted')~bolditalic(H)[plain(LDR)]~'('*'\u2113'[plain(LDR)]==", tau_key, "*')'"),
      parse = TRUE,
      hjust = 0, fontface = "bold.italic",
      size = if(enlarged) 6.0 else 4.7, color = "#303030"
    ) +
    ggplot2::annotate(
      "text", x = col_x, y = common_y_limit * header_y,
      label = c("Period", "D-MMM", gc_label), hjust = c(0, 0.5, 0.5),
      fontface = "bold", size = if(enlarged) 6.0 else 4.4,
      color = c("#303030", "#1f77b4", "#d62728")
    ) +
    ggplot2::annotate(
      "segment", x = if(enlarged) left+.02 else .59, xend = right-.015,
      y = common_y_limit * separator_y, yend = common_y_limit * separator_y,
      linewidth = 0.35, color = "#A0A0A0"
    ) +
    ggplot2::geom_text(
      data = table_df,
      ggplot2::aes(x = col_x[1], y = y, label = period_label),
      inherit.aes = FALSE, hjust = 0, size = if(enlarged) 6.0 else 4.3, color = "#303030"
    ) +
    ggplot2::geom_text(
      data = table_df,
      ggplot2::aes(x = col_x[2], y = y, label = sprintf("%.3f", mmm), color = mmm_color),
      inherit.aes = FALSE, hjust = 0.5, size = if(enlarged) 6.0 else 4.3, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = table_df,
      ggplot2::aes(x = col_x[3], y = y, label = sprintf("%.3f", gc), color = gc_color),
      inherit.aes = FALSE, hjust = 0.5, size = if(enlarged) 6.0 else 4.3, fontface = "bold"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, common_y_limit), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(
      x = if (metric_type == "full") "Hellinger distance" else "LDR Hellinger score",
      y = "Area-weighted density",
      title = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.55),
      panel.grid.minor = ggplot2::element_line(linewidth = 0.35),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.title = ggplot2::element_text(size = 13),
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  density_panel
}

# Same supplementary facet layout as the production plotting script, using the
# precomputed density coordinates and reference-specific summary markers.
#' Replay a reference-faceted PME density figure from its numeric recipe
#' @param recipe Plot-recipe list with coordinates (x,y,method,group,reference),
#'   annotations (x_value,method,statistic,reference) and labels (title,subtitle).
#'   MMM is the internal D-MMM method identifier; GC-0.125 is the other method.
#' @return ggplot object using five facet columns, original curve coordinates and
#'   weighted-mean/median markers. The visible legend says D-MMM. No recomputation
#'   of densities or summaries and no output-file or graphics-device side effects.
plot_cached_reference_density <- function(recipe) {
  ggplot2::ggplot(recipe$coordinates,
    ggplot2::aes(x=x,y=y,colour=method,group=group)) +
    ggplot2::geom_line(linewidth=.8) +
    ggplot2::geom_vline(data=recipe$annotations,
      ggplot2::aes(xintercept=x_value,colour=method,linetype=statistic),
      linewidth=.55,alpha=.9) +
    ggplot2::facet_wrap(~reference,ncol=5) +
    ggplot2::scale_color_manual(values=c("MMM"="#1f77b4","GC-0.125"="#d62728"),
      labels=c("MMM"="D-MMM","GC-0.125"="GC-0.125")) +
    ggplot2::scale_x_continuous(limits=c(0,1)) +
    ggplot2::scale_linetype_manual(values=c("Weighted mean"="dashed","Weighted median"="solid")) +
    ggplot2::labs(title=recipe$labels$title,subtitle=recipe$labels$subtitle,
      x="Hellinger distance",y="Weighted density",color="Method",linetype="Statistic") +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(legend.position="top",plot.title=ggplot2::element_text(face="bold"),
                   strip.text=ggplot2::element_text(face="bold"))
}

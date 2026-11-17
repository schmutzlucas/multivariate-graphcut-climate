#' Compute a scalar weighted empirical quantile of finite samples
#' @param x,w Numeric values and corresponding nonnegative weights of identical
#'   length and shape. Nonfinite values/weights and zero weights are omitted.
#' @param prob Scalar finite probability in [0,1], default median 0.5.
#' @return First sorted positive-weight value reaching prob, or NA_real_ if
#'   no positive finite weight remains. No interpolation or sample correction.
weighted_quantile <- function(x, w, prob = 0.5) {
  if (!is.numeric(x) || !is.numeric(w) || length(x) != length(w) ||
      !identical(dim(x), dim(w)) || any(w[is.finite(w)] < 0) ||
      !is.numeric(prob) || length(prob) != 1L || !is.finite(prob) || prob < 0 || prob > 1) {
    stop("Expected paired numeric values/nonnegative weights and one probability in [0,1].")
  }
  ok <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(ok)) {
    return(NA_real_)
  }

  x_ok <- x[ok]
  w_ok <- w[ok]
  ord <- order(x_ok)
  x_ok <- x_ok[ord]
  w_ok <- w_ok[ord]

  w_sum <- sum(w_ok)
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(NA_real_)
  }

  cdf <- cumsum(w_ok) / w_sum
  cdf[length(cdf)] <- 1
  idx <- which(cdf >= prob)[1]
  if (is.na(idx)) {
    return(tail(x_ok, 1))
  }
  x_ok[idx]
}

#' Estimate and plot Hellinger densities for each held-out reference
#' @param df_metric Data frame with reference, method, value and weight columns;
#'   method uses the internal MMM/GC-0.125 identifiers. Values are scores in [0,1]
#'   and weights are cosine-latitude areas for each regular-grid cell.
#' @param summary_table Data frame with reference, method, weighted_mean and
#'   weighted_median columns, computed from the same cells and weights.
#' @param metric_label Scalar figure title identifying the diagnostic.
#' @param method_colors Named colour vector for the internal method identifiers.
#' @param period_label Target-period text included in the calibration/evaluation subtitle.
#' @return ggplot object with boundary-corrected weighted KDEs (adjust 1, bounds
#'   [0,1]), five facet columns and supplied mean/median markers. No files written.
#'   Visible legend labels translate MMM to D-MMM; raw numerical inputs are retained.
plot_weighted_density_by_reference <- function(
  df_metric,
  summary_table,
  metric_label,
  method_colors,
  period_label
) {
  summary_mean <- summary_table |>
    dplyr::mutate(statistic = "Weighted mean", x_value = weighted_mean)
  summary_median <- summary_table |>
    dplyr::mutate(statistic = "Weighted median", x_value = weighted_median)
  summary_markers <- dplyr::bind_rows(summary_mean, summary_median)
  summary_markers$statistic <- factor(
    summary_markers$statistic,
    levels = c("Weighted mean", "Weighted median")
  )

  subtitle <- paste(
    "Split by PME reference; curves are cosine-latitude weighted kernel densities.",
    "GraphCut data term: 1951-1975 calibration; smoothness/evaluation:",
    period_label,
    "."
  )

  ggplot2::ggplot(
    df_metric,
    ggplot2::aes(x = value, weight = weight, color = method)
  ) +
    ggplot2::geom_density(
      linewidth = 0.8,
      adjust = 1,
      bounds = c(0, 1)
    ) +
    ggplot2::geom_vline(
      data = summary_markers,
      ggplot2::aes(xintercept = x_value, color = method, linetype = statistic),
      linewidth = 0.55,
      alpha = 0.9
    ) +
    ggplot2::facet_wrap(~reference, ncol = 5) +
    ggplot2::scale_color_manual(values = method_colors,
      labels = function(x) ifelse(x == "MMM", "D-MMM", x)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_linetype_manual(values = c("Weighted mean" = "dashed", "Weighted median" = "solid")) +
    ggplot2::labs(
      title = metric_label,
      subtitle = subtitle,
      x = "Hellinger distance",
      y = "Weighted density",
      color = "Method",
      linetype = "Statistic"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

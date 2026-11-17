#' Resolve one raw PME density map without reinterpreting its spatial axes
#' @param map_obj Numeric matrix of the expected shape, or a dimensionless
#'   legacy vector with longitude-fastest column-major ordering.
#' @param expected_dim Positive integer c(nlon,nlat) dimensions.
#' @param label Scientific map identity included in an invalid-shape error.
#' @return Lon x lat matrix. Wrong-shaped arrays are rejected even if total
#'   length agrees; no transposition/regridding is guessed from array size.
coerce_to_map_matrix <- function(map_obj, expected_dim, label) {
  if (is.matrix(map_obj) || is.array(map_obj)) {
    d <- dim(map_obj)
    if (length(d) == 2L && identical(as.integer(d), as.integer(expected_dim))) {
      return(as.matrix(map_obj))
    }
  }
  if (is.numeric(map_obj) && is.null(dim(map_obj)) && length(map_obj) == prod(expected_dim)) {
    return(matrix(as.numeric(map_obj), nrow = expected_dim[1], ncol = expected_dim[2]))
  }
  stop(
    label,
    " could not be coerced to [",
    paste(expected_dim, collapse = ", "),
    "]. Class=",
    paste(class(map_obj), collapse = "/"),
    ", dim=",
    if (is.null(dim(map_obj))) "NULL" else paste(dim(map_obj), collapse = ", "),
    "."
  )
}

#' Compute the raw-density helper's finite-pair weighted mean
#' @param x,w Numeric observation/weight vectors with identical lengths.
#' @return Weighted arithmetic mean over finite pairs with nonnegative weights;
#'   NA_real_ for no positive valid total. References are not normalised separately.
weighted_mean <- function(x, w) {
  if (!is.numeric(x) || !is.numeric(w) || length(x) != length(w))
    stop("Values and weights must be numeric vectors of the same length.")
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok)) {
    return(NA_real_)
  }
  w_sum <- sum(w[ok])
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(NA_real_)
  }
  sum(x[ok] * w[ok]) / w_sum
}

#' Build regular-grid cosine-latitude weights for pooled PME densities
#' @param lon Nonempty unique finite longitude vector, determining row count.
#' @param lat Nonempty unique finite latitude vector in degrees north, within [-90,90].
#' @return Unnormalised lon x lat cos(latitude) matrix. Assumes equal-width
#'   longitudes; this is not a general irregular-grid area calculator.
build_weight_matrix <- function(lon, lat) {
  if (!is.numeric(lon) || !length(lon) || any(!is.finite(lon)) || anyDuplicated(lon) ||
      !is.numeric(lat) || !length(lat) || any(!is.finite(lat)) ||
      any(abs(lat)>90) || anyDuplicated(lat))
    stop("Expected unique finite longitude/latitude vectors with latitudes in [-90,90].")
  matrix(
    rep(cos(lat * pi / 180), each = length(lon)),
    nrow = length(lon),
    ncol = length(lat),
    byrow = FALSE
  )
}

#' Require complete PME reference identities before pooling raw density maps
#' @param results Compact PME result list with metadata and named method maps.
#' @param period_id Evaluation-period key.
#' @param lambda_key Character GC smoothness key.
#' @param tau_key Character LDR threshold key, ignored for full Hellinger.
#' @param metric_type Exactly full or partial, the internal LDR identifier.
#' @return Metadata model order after verifying both methods contain exactly
#'   that set once. Does not silently intersect or drop missing references.
get_reference_order <- function(results, period_id, lambda_key, tau_key = NULL, metric_type) {
  if (length(metric_type)!=1L || !metric_type %in% c("full","partial"))
    stop("metric_type must be full or partial.")
  if (metric_type == "full") {
    ref_sets <- list(
      names(results$MMM_hdist_future_list[[period_id]]),
      names(results$GC_hdist_future_list[[period_id]][[lambda_key]])
    )
  } else {
    ref_sets <- list(
      names(results$MMM_partial_hdist_future_list[[period_id]][[tau_key]]),
      names(results$GC_partial_hdist_future_list[[period_id]][[tau_key]][[lambda_key]])
    )
  }
  ref_order <- results$meta$model_names
  if (!is.character(ref_order) || !length(ref_order) || anyNA(ref_order) ||
      any(!nzchar(ref_order)) || anyDuplicated(ref_order) ||
      any(vapply(ref_sets,function(x)anyDuplicated(x)>0L || !setequal(x,ref_order),logical(1))))
    stop("Incomplete or unexpected held-out reference set for ",period_id,
      "; both D-MMM and GC must contain every metadata model exactly once.")
  ref_order
}

#' Stack raw spatial PME distance fields for one pooled density panel
#' @param results Compact PME result list containing both methods for every reference.
#' @param period_id Evaluation-period key.
#' @param lambda_key GC smoothness key used to construct the internal GC method label.
#' @param tau_key LDR threshold key, ignored for full Hellinger.
#' @param metric_type Full Hellinger or partial internal LDR identifier.
#' @param lon,lat Coordinate vectors determining expected map shape.
#' @param weight_matrix Matching lon x lat cosine-latitude weight matrix.
#' @return Data frame with method factor, value and weight. Internal MMM means
#'   D-MMM; nonfinite values/weights and negative weights are removed. This is
#'   a pooled grid-cell density, not the equal-reference distribution of relative
#'   improvements; no per-reference weight renormalisation is applied.
build_period_density_data <- function(results, period_id, lambda_key, tau_key, metric_type, lon, lat, weight_matrix) {
  expected_dim <- c(length(lon), length(lat))
  ref_order <- get_reference_order(results, period_id, lambda_key, tau_key, metric_type)
  weight_vector <- as.vector(weight_matrix)
  rows <- vector("list", length(ref_order) * 2L)
  cursor <- 1L

  if (metric_type == "full") {
    mmm_maps <- results$MMM_hdist_future_list[[period_id]]
    gc_maps <- results$GC_hdist_future_list[[period_id]][[lambda_key]]
  } else {
    mmm_maps <- results$MMM_partial_hdist_future_list[[period_id]][[tau_key]]
    gc_maps <- results$GC_partial_hdist_future_list[[period_id]][[tau_key]][[lambda_key]]
  }

  for (ref_name in ref_order) {
    mmm_map <- coerce_to_map_matrix(mmm_maps[[ref_name]], expected_dim, paste0("MMM map ", period_id, "/", ref_name))
    gc_map <- coerce_to_map_matrix(gc_maps[[ref_name]], expected_dim, paste0("GC map ", period_id, "/", ref_name))

    rows[[cursor]] <- data.frame(
      method = "MMM",
      value = as.vector(mmm_map),
      weight = weight_vector,
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L

    rows[[cursor]] <- data.frame(
      method = paste0("GC-", lambda_key),
      value = as.vector(gc_map),
      weight = weight_vector,
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }

  out <- do.call(rbind, rows)
  out <- out[is.finite(out$value) & is.finite(out$weight) & out$weight >= 0, , drop = FALSE]
  out$method <- factor(out$method, levels = c("MMM", paste0("GC-", lambda_key)))
  out
}

#' Summarise method-wise area-weighted means of stacked PME density samples
#' @param df_period Data frame with method, value and weight columns.
#' @return Tibble with one row per method and weighted_mean, computed over finite
#'   value/weight pairs. Requires dplyr; does not compute reference-wise improvements.
summarise_means <- function(df_period) {
  df_period |>
    dplyr::group_by(method) |>
    dplyr::summarise(weighted_mean = weighted_mean(value, weight), .groups = "drop")
}

#' Find the maximum kernel-density height for a common overlay axis
#' @param df_period Stacked method/value/weight samples for one period.
#' @param method_colors Named colour vector keyed by internal method labels.
#' @return Scalar maximum of ggplot2's evaluated area-weighted density curves,
#'   using adjust=1 and bounds/limits [0,1]. Does not save or display the plot.
density_ymax <- function(df_period, method_colors) {
  p <- ggplot2::ggplot(
    df_period,
    ggplot2::aes(x = value, weight = weight, color = method)
  ) +
    ggplot2::geom_density(linewidth = 0.95, adjust = 1, bounds = c(0, 1)) +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::scale_x_continuous(limits = c(0, 1))

  y_values <- ggplot2::ggplot_build(p)$data[[1]]$y
  max(y_values, na.rm = TRUE)
}

#' Construct the five-period Article 1 PME density overlay from raw compact maps
#'
#' Blue shades denote D-MMM and red shades GC, with period-specific intensity.
#' Includes an area-weighted mean table, no vertical mean markers. This raw-map
#' renderer supports re-curating numerical curve coordinates; Level 1 publication
#' figures normally use the compact cached curves instead of all spatial maps.
#' @param results Compact PME result list with metadata, periods and complete maps.
#' @param lambda_key GC smoothness key, normally "0.125".
#' @param output_base Output stem for PDF/PNG, without extension. Parent directory
#'   must exist; saving can overwrite these filenames, so select a new stem for review.
#' @param save_dpi Raster resolution passed to ggsave.
#' @param metric_type Full ordinary Hellinger or partial internal LDR identifier.
#' @param tau_key Character LDR threshold, used for maps and the table label;
#'   ignored in full-Hellinger mode.
#' @param save_outputs Whether to save the 29 x 18 cm PDF and PNG; FALSE constructs
#'   the plot for numerical curve caching without writing figure files.
#' @return ggplot object invisibly, with pooled density curves and a summary table.
#'   Does not run PMF estimation, graph cut or select/repair an input dataset.
make_overlay_density_figure <- function(results, lambda_key, output_base, save_dpi,
                                        metric_type = "full", tau_key = NULL, save_outputs = TRUE) {
  lon <- as.numeric(results$meta$lon)
  lat <- as.numeric(results$meta$lat)
  periods <- results$periods
  period_order <- as.character(periods$period_id)
  period_labels <- as.character(periods$label)
  names(period_labels) <- period_order

  weight_matrix <- build_weight_matrix(lon, lat)
  gc_label <- paste0("GC-", lambda_key)
  mmm_shades <- grDevices::colorRampPalette(c("#9ECAE1", "#08519C"))(length(period_order))
  gc_shades <- grDevices::colorRampPalette(c("#FCAE91", "#99000D"))(length(period_order))

  period_data <- vector("list", length(period_order))
  period_means <- vector("list", length(period_order))
  for (idx in seq_along(period_order)) {
    period_id <- period_order[[idx]]
    period_label <- period_labels[[period_id]]
    period_data[[idx]] <- build_period_density_data(
      results = results,
      period_id = period_id,
      lambda_key = lambda_key,
      tau_key = tau_key,
      metric_type = metric_type,
      lon = lon,
      lat = lat,
      weight_matrix = weight_matrix
    ) |>
      dplyr::mutate(
        period_id = period_id,
        period_label = period_label,
        curve_id = paste(method, period_id, sep = "__"),
        curve_color = dplyr::if_else(method == "MMM", mmm_shades[[idx]], gc_shades[[idx]])
      )
    period_means[[idx]] <- period_data[[idx]] |>
      dplyr::group_by(method, period_id, period_label, curve_color) |>
      dplyr::summarise(weighted_mean = weighted_mean(value, weight), .groups = "drop")
  }

  overlay_data <- do.call(rbind, period_data)
  overlay_means <- do.call(rbind, period_means)
  common_y_limit <- max(
    vapply(
      period_data,
      density_ymax,
      numeric(1),
      method_colors = c("MMM" = "#1f77b4", setNames("#d62728", gc_label))
    ),
    na.rm = TRUE
  ) * 1.08

  table_df <- data.frame(
    period_label = period_labels[period_order],
    mmm = vapply(period_order, function(id) {
      overlay_means$weighted_mean[overlay_means$period_id == id & overlay_means$method == "MMM"]
    }, numeric(1)),
    gc = vapply(period_order, function(id) {
      overlay_means$weighted_mean[overlay_means$period_id == id & overlay_means$method == gc_label]
    }, numeric(1)),
    mmm_color = mmm_shades,
    gc_color = gc_shades,
    y = seq(common_y_limit * 0.79, common_y_limit * 0.485, length.out = length(period_order)),
    stringsAsFactors = FALSE
  )

  table_box <- data.frame(
    xmin = 0.57, xmax = 0.985,
    ymin = common_y_limit * 0.405, ymax = common_y_limit * 0.975
  )

  density_panel <- ggplot2::ggplot(
    overlay_data,
    ggplot2::aes(x = value, weight = weight, group = curve_id, color = curve_color)
  ) +
    ggplot2::geom_density(linewidth = 0.85, adjust = 1, bounds = c(0, 1), alpha = 1) +
    ggplot2::geom_rect(
      data = table_box,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white",
      color = "#6A6A6A",
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "text", x = 0.595, y = common_y_limit * 0.915,
      label = if (metric_type == "full") "bolditalic('Area-weighted H')" else
        paste0("bolditalic('Area-weighted')~bolditalic(H)[plain(LDR)]~'('*'\u2113'[plain(LDR)]==", tau_key, "*')'"),
      parse = TRUE,
      hjust = 0, fontface = "bold.italic",
      size = 4.7, color = "#303030"
    ) +
    ggplot2::annotate(
      "text", x = c(0.595, 0.81, 0.935), y = common_y_limit * 0.85,
      label = c("Period", "D-MMM", gc_label), hjust = c(0, 0.5, 0.5),
      fontface = "bold", size = 4.4,
      color = c("#303030", "#1f77b4", "#d62728")
    ) +
    ggplot2::annotate(
      "segment", x = 0.59, xend = 0.97,
      y = common_y_limit * 0.825, yend = common_y_limit * 0.825,
      linewidth = 0.35, color = "#A0A0A0"
    ) +
    ggplot2::geom_text(
      data = table_df,
      ggplot2::aes(x = 0.595, y = y, label = period_label),
      inherit.aes = FALSE, hjust = 0, size = 4.3, color = "#303030"
    ) +
    ggplot2::geom_text(
      data = table_df,
      ggplot2::aes(x = 0.81, y = y, label = sprintf("%.3f", mmm), color = mmm_color),
      inherit.aes = FALSE, hjust = 0.5, size = 4.3, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = table_df,
      ggplot2::aes(x = 0.935, y = y, label = sprintf("%.3f", gc), color = gc_color),
      inherit.aes = FALSE, hjust = 0.5, size = 4.3, fontface = "bold"
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

  if (save_outputs) ggplot2::ggsave(paste0(output_base, ".pdf"), plot = density_panel, width = 29, height = 18, units = "cm", dpi = save_dpi, device = grDevices::cairo_pdf)
  if (save_outputs) ggplot2::ggsave(paste0(output_base, ".png"), plot = density_panel, width = 29, height = 18, units = "cm", dpi = save_dpi)
  invisible(density_panel)
}

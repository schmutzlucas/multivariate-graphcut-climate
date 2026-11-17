# Pure definitions extracted from the final Article 1 implementation.
a4_plot_width_cm <- 17.5

#' Build a map-row title with the manuscript precipitation-unit typography
#' @param row_title Character title or plotmath language/expression object. The
#'   exact precipitation-bias string receives superscript inverse-day formatting.
#' @return grid text grob at the original 10.5-point size; no device is opened.
make_row_title_grob <- function(row_title) {
  if (identical(row_title, "Precipitation bias [mm day^-1]")) {
    return(grid::textGrob(
      expression(bold("Precipitation bias [mm day"^-1*"]")),
      gp = grid::gpar(fontsize = 10.5)
    ))
  }
  grid::textGrob(row_title, gp = grid::gpar(fontsize = 10.5, fontface = "bold"))
}

#' Build the shared two-method header for a four-map figure
#' @param left_title,right_title Scalar character method labels. MMM denotes the
#'   mean-field baseline; distributional callers explicitly supply D-MMM.
#' @return grid tree containing two 13-point titles at the fixed manuscript positions.
make_column_title_grob <- function(left_title = "MMM", right_title = "GC-0.125") {
  grid::gTree(children = grid::gList(
    grid::textGrob(left_title, x = grid::unit(0.31, "npc"), y = grid::unit(0.5, "npc"),
                   gp = grid::gpar(fontsize = 13, fontface = "plain")),
    grid::textGrob(right_title, x = grid::unit(0.68, "npc"), y = grid::unit(0.5, "npc"),
                   gp = grid::gpar(fontsize = 13, fontface = "plain"))
  ))
}

#' Validate paired regular-grid maps before constructing plotting coordinates
#' @param baseline,candidate Numeric longitude x latitude matrices with identical
#'   dimensions. Missing values are permitted; infinities are rejected.
#' @param lon,lat Finite distinct grid centres matching the matrix dimensions.
#'   Longitude lies in [-180,360) without duplicated periodic cells; latitude is
#'   in [-90,90]. Callers supply the regular cylindrical Article 1 grid.
#' @param limits Two finite strictly increasing colour-scale endpoints.
#' @return TRUE invisibly, or an error for incompatible/ambiguous map identities.
validate_era5_map_pair <- function(baseline, candidate, lon, lat, limits) {
  if (!is.matrix(baseline) || !is.numeric(baseline) || !is.matrix(candidate) ||
      !is.numeric(candidate) || !identical(dim(baseline), dim(candidate)) ||
      !identical(dim(baseline), c(length(lon), length(lat))) ||
      any(is.infinite(baseline)) || any(is.infinite(candidate))) {
    stop("Expected numeric lon x lat maps with identical coordinate-matched dimensions and no infinities.")
  }
  if (!is.numeric(lon) || !length(lon) || any(!is.finite(lon)) ||
      any(lon < -180 | lon >= 360) || anyDuplicated(ifelse(lon >= 180, lon - 360, lon)) ||
      !is.numeric(lat) || !length(lat) || any(!is.finite(lat)) ||
      any(abs(lat) > 90) || anyDuplicated(lat)) {
    stop("Expected finite distinct cylindrical longitude/latitude coordinates in range.")
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) || limits[1] >= limits[2]) {
    stop("limits must have two finite strictly increasing endpoints.")
  }
  invisible(TRUE)
}

#' Render a pair of ordinary or restricted Hellinger maps
#' @param h_mmm,h_gc Numeric longitude x latitude D-MMM/GC score matrices in [0,1];
#'   missing cells remain missing and values above limits are visually squished.
#' @param lon,lat Regular cylindrical grid centres matching the matrices.
#' @param limits Two finite increasing shared colour-scale endpoints.
#' @param palette_high Colour at the upper endpoint of a white-to-colour gradient.
#' @param left_title,right_title Scalar character method panel titles; empty hides.
#' @param row_title Character or plotmath title for the shared row; empty hides.
#' @param show_x_title Logical; include the Longitude axis title when TRUE.
#' @return patchwork plot with equal-size maps, shared scale/legend and the
#'   manuscript coastline, crop and typography. No scores or files are modified.
build_hdist3_pair_plot <- function(
  h_mmm,
  h_gc,
  lon,
  lat,
  limits = c(0, 0.5),
  palette_high = "#015a8c",
  left_title = "D-MMM",
  right_title = "GC-0.125",
  row_title = "Hellinger distance",
  show_x_title = TRUE
) {
  required_pkgs <- c("ggplot2", "reshape2", "patchwork", "maps", "scales", "grid")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }

  validate_era5_map_pair(h_mmm, h_gc, lon, lat, limits)
  if (any(h_mmm[is.finite(h_mmm)] < 0 | h_mmm[is.finite(h_mmm)] > 1 + 1e-8) ||
      any(h_gc[is.finite(h_gc)] < 0 | h_gc[is.finite(h_gc)] > 1 + 1e-8)) {
    stop("Hellinger score maps must lie in [0,1] up to numerical tolerance.")
  }

  #' Convert the enclosing figure field into geographic plot rows
  #' Internal closure of build_hdist3_pair_plot; not a standalone entry point.
  #' @param field_matrix Longitude-by-latitude diagnostic or bias matrix.
  #' @param method_label Method label attached to all Hellinger plot rows.
  #' @return Data frame of grid indices, diagnostic values and wrapped longitudes/latitudes from the enclosing coordinate vectors.
  to_map_df <- function(field_matrix, method_label) {
    map_df <- reshape2::melt(
      field_matrix,
      varnames = c("lon_idx", "lat_idx"),
      value.name = "H_dist"
    )
    map_df$lon <- lon[map_df$lon_idx]
    map_df$lat <- lat[map_df$lat_idx]
    map_df$lon_wrapped <- ifelse(map_df$lon > 180, map_df$lon - 360, map_df$lon)
    map_df$method <- method_label
    map_df
  }

  mmm_df <- to_map_df(h_mmm, "MMM")
  gc_df  <- to_map_df(h_gc, "GC-0.125")

  #' Draw one panel of the enclosing paired geographic figure
  #' Internal closure of build_hdist3_pair_plot; not a standalone entry point.
  #' @param map_df Long-form geographic diagnostic data frame.
  #' @param panel_title Character title for this panel.
  #' @param right_panel Whether to use the paired right-panel axis treatment.
  #' @return ggplot object using the enclosing limits, labels, units and colour scale. Call through the enclosing paired-plot function.
  build_map_panel <- function(map_df, panel_title, right_panel = FALSE) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = map_df,
        ggplot2::aes(x = lon_wrapped, y = lat, fill = H_dist)
      ) +
      ggplot2::borders("world", colour = "grey45", linewidth = 0.08) +
      ggplot2::coord_fixed(
        ratio = 1.3,
        xlim = c(-180, 180),
        ylim = c(-84, 90),
        expand = FALSE
      ) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0),
        breaks = c(-100, 0, 100)
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        breaks = c(-50, 0, 50)
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = palette_high,
        limits = limits,
        breaks = seq(limits[1], limits[2], length.out = 3),
        oob = scales::squish,
        name = "H",
        guide = ggplot2::guide_colorbar(
          barheight = grid::unit(3.8, "cm"),
          barwidth  = grid::unit(0.35, "cm"),
          frame.colour = NA,
          ticks.colour = "grey70"
        )
      ) +
      ggplot2::labs(
        title = if (nzchar(panel_title)) panel_title else NULL,
        x = if (show_x_title) "Longitude" else " ",
        y = "Latitude"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 14,
          face = "plain",
          hjust = 0.5,
          margin = ggplot2::margin(t = 2, r = 0, b = 2, l = 0)
        ),

        axis.title.x = ggplot2::element_text(size = 7, margin = ggplot2::margin(t = 2)),
        axis.title.y = ggplot2::element_text(size = 7, margin = ggplot2::margin(r = 2)),
        axis.text.x  = ggplot2::element_text(size = 5),
        axis.text.y  = ggplot2::element_text(size = 5),

        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "grey35", linewidth = 0.35, fill = NA),

        legend.position = "right",
        legend.title = ggplot2::element_text(size = 8.5),
        legend.text  = ggplot2::element_text(size = 7.5),
        legend.key.height = grid::unit(0.9, "cm"),
        legend.key.width  = grid::unit(0.35, "cm"),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        legend.margin = ggplot2::margin(0, 0, 0, -7),

        # Small top margin, but not negative
        plot.margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0)
      )

    if (right_panel) {
      p <- p +
        ggplot2::labs(y = NULL) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()
        )
    }

    p
  }

  p_mmm <- build_map_panel(mmm_df, left_title, right_panel = FALSE)
  p_gc  <- build_map_panel(gc_df, right_title, right_panel = TRUE)

  p_pair <- (p_mmm + p_gc) +
    patchwork::plot_layout(
      guides = "collect",
      widths = c(1, 1)
    ) &
    ggplot2::theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      # Important: keep this at zero, not negative
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  if (is.language(row_title) || is.expression(row_title) || nzchar(row_title)) {
    row_title_grob <- patchwork::wrap_elements(
      full = make_row_title_grob(row_title)
    )
    p_pair <- row_title_grob / p_pair +
      patchwork::plot_layout(heights = c(0.045, 1))
  }

  p_pair
}

#' Render a pair of signed mean-field bias maps
#' @param bias_mmm,bias_gc Numeric longitude x latitude mean-field bias matrices.
#'   Missing cells remain missing; values beyond limits are visually squished.
#' @param lon,lat Regular cylindrical grid centres matching the matrices.
#' @param limits Two finite increasing shared colour-scale endpoints.
#' @param fill_name Colour-scale title, or NULL to omit.
#' @param left_title,right_title Scalar character method titles; empty hides.
#' @param row_title Scalar character shared row title; empty hides.
#' @param show_x_title Logical; include the Longitude axis title when TRUE.
#' @return patchwork plot retaining the original diverging palette, shared legend,
#'   crop and typography. No unit conversion, numerical evaluation or file writing.
build_bias_pair_plot <- function(
  bias_mmm,
  bias_gc,
  lon,
  lat,
  limits = c(-8, 8),
  fill_name = NULL,
  left_title = "MMM",
  right_title = "GC-0.125",
  row_title = "",
  show_x_title = TRUE
) {
  required_pkgs <- c("ggplot2", "reshape2", "patchwork", "maps", "scales", "grid", "RColorBrewer")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }
  validate_era5_map_pair(bias_mmm, bias_gc, lon, lat, limits)

  #' Convert the enclosing figure field into geographic plot rows
  #' Internal closure of build_bias_pair_plot; not a standalone entry point.
  #' @param field_matrix Longitude-by-latitude diagnostic or bias matrix.
  #' @return Data frame of grid indices, diagnostic values and wrapped longitudes/latitudes from the enclosing coordinate vectors.
  to_map_df <- function(field_matrix) {
    map_df <- reshape2::melt(
      field_matrix,
      varnames = c("lon_idx", "lat_idx"),
      value.name = "Bias"
    )
    map_df$lon <- lon[map_df$lon_idx]
    map_df$lat <- lat[map_df$lat_idx]
    map_df$lon_wrapped <- ifelse(map_df$lon > 180, map_df$lon - 360, map_df$lon)
    map_df
  }

  #' Draw one panel of the enclosing paired geographic figure
  #' Internal closure of build_bias_pair_plot; not a standalone entry point.
  #' @param map_df Long-form geographic diagnostic data frame.
  #' @param panel_title Character title for this panel.
  #' @param right_panel Whether to use the paired right-panel axis treatment.
  #' @return ggplot object using the enclosing limits, labels, units and colour scale. Call through the enclosing paired-plot function.
  build_map_panel <- function(map_df, panel_title, right_panel = FALSE) {
    bias_colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
    bias_breaks <- as.numeric(format(seq(limits[1], limits[2], length.out = 6), digits = 3))

    p <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = map_df,
        ggplot2::aes(x = lon_wrapped, y = lat, fill = Bias)
      ) +
      ggplot2::borders("world", colour = "grey45", linewidth = 0.08) +
      ggplot2::coord_fixed(
        ratio = 1.3,
        xlim = c(-180, 180),
        ylim = c(-84, 90),
        expand = FALSE
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0), breaks = c(-100, 0, 100)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), breaks = c(-50, 0, 50)) +
      ggplot2::scale_fill_gradientn(
        colours = bias_colors,
        limits = limits,
        breaks = bias_breaks,
        oob = scales::squish,
        name = fill_name,
        guide = ggplot2::guide_colorbar(
          barheight = grid::unit(3.8, "cm"),
          barwidth = grid::unit(0.35, "cm"),
          frame.colour = NA,
          ticks.colour = "grey70"
        )
      ) +
      ggplot2::labs(
        title = if (nzchar(panel_title)) panel_title else NULL,
        x = if (show_x_title) "Longitude" else " ",
        y = "Latitude"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 14,
          face = "plain",
          hjust = 0.5,
          margin = ggplot2::margin(t = 2, r = 0, b = 2, l = 0)
        ),
        axis.title.x = ggplot2::element_text(size = 7, margin = ggplot2::margin(t = 2)),
        axis.title.y = ggplot2::element_text(size = 7, margin = ggplot2::margin(r = 2)),
        axis.text.x = ggplot2::element_text(size = 5),
        axis.text.y = ggplot2::element_text(size = 5),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "grey35", linewidth = 0.35, fill = NA),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 8.5),
        legend.text = ggplot2::element_text(size = 7.5),
        legend.key.height = grid::unit(0.9, "cm"),
        legend.key.width = grid::unit(0.35, "cm"),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        legend.margin = ggplot2::margin(0, 0, 0, -7),
        plot.margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0)
      )

    if (right_panel) {
      p <- p + ggplot2::labs(y = NULL) + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    }
    p
  }

  p_mmm <- build_map_panel(to_map_df(bias_mmm), left_title, right_panel = FALSE)
  p_gc <- build_map_panel(to_map_df(bias_gc), right_title, right_panel = TRUE)

  p_pair <- (p_mmm + p_gc) +
    patchwork::plot_layout(guides = "collect", widths = c(1, 1)) &
    ggplot2::theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  if (nzchar(row_title)) {
    row_title_grob <- patchwork::wrap_elements(
      full = make_row_title_grob(row_title)
    )
    p_pair <- row_title_grob / p_pair +
      patchwork::plot_layout(heights = c(0.045, 1))
  }

  p_pair
}

#' Compute mean-field ensemble and selected-model biases against a reference
#' @param mean_future Named list of numeric longitude x latitude x model mean
#'   arrays. Each array has reference at index 1 and candidate models afterwards.
#' @param gc_labels Finite integer longitude x latitude matrix indexing candidates
#'   from 1, so label k selects full-array model k+1.
#' @param var_name Scalar character variable key in mean_future.
#' @return List bias_mmm and bias_gc of longitude x latitude matrices, each method
#'   minus reference. The MMM is an equal-model mean; missing model values propagate.
#'   No precipitation/unit conversion is performed and invalid labels raise an error.
compute_mmm_gc_bias_fields <- function(mean_future, gc_labels, var_name) {
  if (!is.character(var_name) || length(var_name) != 1L || is.na(var_name) ||
      !var_name %in% names(mean_future)) {
    stop("Variable `", var_name, "` is missing in `mean_future`.", call. = FALSE)
  }
  var_array <- mean_future[[var_name]]
  if (!is.numeric(var_array) || length(dim(var_array)) != 3 || dim(var_array)[3] < 2) {
    stop("`mean_future[[var_name]]` must be [lon, lat, model] with at least 2 models.", call. = FALSE)
  }
  grid <- dim(var_array)[1:2]
  if (!is.numeric(gc_labels) || !identical(dim(gc_labels), grid) ||
      any(!is.finite(gc_labels)) || any(gc_labels != floor(gc_labels)) ||
      any(gc_labels < 1L | gc_labels >= dim(var_array)[3])) {
    stop("GraphCut labels must be finite integer candidate indices on the mean-field grid.", call. = FALSE)
  }

  ref_field <- matrix(var_array[, , 1], nrow = grid[1], ncol = grid[2])
  mmm_field <- apply(var_array[, , 2:dim(var_array)[3], drop = FALSE], c(1, 2), mean)

  gc_field <- matrix(NA_real_, nrow = dim(var_array)[1], ncol = dim(var_array)[2])
  label_values <- sort(unique(as.integer(gc_labels)))
  for (label in label_values) {
    model_idx <- label + 1
    is_label <- gc_labels == label
    candidate <- matrix(var_array[, , model_idx], nrow = grid[1], ncol = grid[2])
    gc_field[is_label] <- candidate[is_label]
  }

  list(
    bias_mmm = mmm_field - ref_field,
    bias_gc = gc_field - ref_field
  )
}


# Pure rendering of the approved cached Article 1 fields. No experiment code.
# Geometry is in mm at final publication size, shared by both figures.
a4_plot_width_cm <- 17.5
#' Describe fixed physical geometry of the accepted ERA5 four-panel figure
#' @return Named list of millimetre dimensions, degree extents, aspect ratio and axis breaks.
era5_map_geometry <- function() {
  widths <- c(outer_margin=.3, latitude_axis=4.7, map_left=78, column_gap=1.5,
              map_right=78, colourbar_gap=1.5, colourbar=11)
  map_height <- 78 * 1.3 * diff(c(-84,90)) / diff(c(-180,180))
  # Reuse the existing header bands in method-then-diagnostic order. The numeric
  # heights and all map/colourbar/axis positions remain unchanged.
  heights <- c(top_margin=.6, column_heading=4.2, title_top=4, heading_gap=1.3,
    map_top=map_height, row_gap=1.7, title_bottom=4.5, title_gap=1.3,
    map_bottom=map_height, longitude_axis=4.8, bottom_margin=.5)
  list(widths=widths, heights=heights, width_mm=sum(widths),height_mm=sum(heights),
    map_width_mm=78,map_height_mm=map_height,map_rows=c(5L,9L),title_rows=c(3L,7L),method_row=2L,
    xlim=c(-180,180),ylim=c(-84,90),aspect=1.3,
    longitude_breaks=c(-100,0,100),latitude_breaks=c(-50,0,50))
}
fourmap_total_height_cm <- era5_map_geometry()$height_mm / 10

#' Select diagnostic row headings
#' @param kind Diagnostic kind: distributional or bias for layouts, hellinger or bias for scales.
#' @return List of two text or plotmath labels; unknown kinds fail.
era5_row_labels <- function(kind) {
  switch(kind,
    distributional=list('Hellinger distance',
      expression('LDR Hellinger score ('*'\u2113'[plain(LDR)]==0.01*')')),
    bias=list('Temperature bias','Precipitation bias'),stop('Unknown ERA5 figure kind.'))
}

#' Construct a nine-point diagnostic heading
#' @param name Grid grob identifier.
#' @param label Diagnostic text or plotmath label.
#' @return Grid text grob; no device is opened.
era5_row_title_grob <- function(label, name='row-heading') {
  grid::textGrob(label,x=.5,y=.5,
    gp=grid::gpar(fontfamily='sans',fontsize=9,fontface='plain'),name=name)
}

#' Select mean-bias colourbar units
#' @param kind Diagnostic kind: distributional or bias for layouts, hellinger or bias for scales.
#' @param row One-based diagnostic row (1 or 2).
#' @return Kelvin or mm/day label, or NULL for distributional maps.
era5_colourbar_unit <- function(kind,row) {
  if(kind!='bias') return(NULL)
  if(row==1) '[K]' else expression('['*plain(mm)~plain(day)^{-1}*']')
}

#' Construct a ten-point method heading
#' @param name Grid grob identifier.
#' @param title Method heading text.
#' @return Grid text grob positioned in its map column.
era5_column_title_grob <- function(title, name='column-heading') {
  # Placed in the SAME gtable column as its geographic panel, never a page fraction.
  grid::textGrob(title,x=.5,y=.5,
    gp=grid::gpar(fontfamily='sans',fontsize=10,fontface='plain'),name=name)
}

#' Lay out cached fields with wrapped longitude coordinates
#' @param field Numeric longitude x latitude matrix; supplied values are not normalised or imputed.
#' @param lon,lat Longitude and latitude cell centres in degrees, matching matrix axes.
#' @return Data frame with longitude, latitude, wrapped longitude and unchanged values; shape mismatch fails.
era5_map_data <- function(field,lon,lat) {
  stopifnot(is.matrix(field),identical(dim(field),c(length(lon),length(lat))))
  data.frame(lon=rep(lon,length(lat)),lat=rep(lat,each=length(lon)),
    lon_wrapped=rep(ifelse(lon>180,lon-360,lon),length(lat)),value=as.vector(field))
}

#' Construct the accepted diagnostic colour scale
#' @param kind Diagnostic kind: distributional or bias for layouts, hellinger or bias for scales.
#' @param limits Numeric two-element colour limits in the diagnostic units.
#' @return ggplot2 fill scale with supplied limits and squishing outside them; no values modified.
era5_fill_scale <- function(kind,limits) {
  if(kind=='hellinger') return(ggplot2::scale_fill_gradient(
    low='white',high='#015a8c',limits=limits,
    breaks=seq(limits[1],limits[2],length.out=3),oob=scales::squish,name=NULL))
  ggplot2::scale_fill_gradientn(colours=rev(RColorBrewer::brewer.pal(11,'RdBu')),
    limits=limits,breaks=as.numeric(format(seq(limits[1],limits[2],length.out=6),digits=3)),
    oob=scales::squish,name=NULL)
}

#' Render a cached regular-grid diagnostic panel
#' @param kind Diagnostic kind: distributional or bias for layouts, hellinger or bias for scales.
#' @param field Numeric longitude x latitude matrix; supplied values are not normalised or imputed.
#' @param lon,lat Longitude and latitude cell centres in degrees, matching matrix axes.
#' @param limits Numeric two-element colour limits in the diagnostic units.
#' @param left,bottom Logical flags controlling visible latitude and longitude axes.
#' @param render_mode Raster (default) or tile primitive; changes rendering only, not arrays, colours or geometry.
#' @return ggplot object with unchanged data and scales; raster fields prevent PDF seams, while overlays stay vector.
era5_map_panel <- function(field,lon,lat,kind,limits,left,bottom,render_mode = c("raster", "tile")) {
  render_mode <- match.arg(render_mode)
  field_geom <- if (render_mode == "raster") ggplot2::geom_raster else ggplot2::geom_tile
  geometry <- era5_map_geometry()
  ggplot2::ggplot(era5_map_data(field,lon,lat),ggplot2::aes(lon_wrapped,lat,fill=value))+
    field_geom()+
    ggplot2::borders('world',colour='grey45',linewidth=.08)+
    ggplot2::coord_fixed(ratio=geometry$aspect,xlim=geometry$xlim,ylim=geometry$ylim,expand=FALSE)+
    ggplot2::scale_x_continuous(expand=c(0,0),breaks=geometry$longitude_breaks)+
    ggplot2::scale_y_continuous(expand=c(0,0),breaks=geometry$latitude_breaks)+
    era5_fill_scale(kind,limits)+ggplot2::labs(x=NULL,y=NULL)+
    ggplot2::theme_bw(base_family='sans')+
    ggplot2::theme(axis.title=ggplot2::element_blank(),
      axis.text=ggplot2::element_text(size=7,colour='grey25'),
      axis.text.x=if(bottom) ggplot2::element_text(margin=ggplot2::margin(t=1)) else ggplot2::element_blank(),
      axis.text.y=if(left) ggplot2::element_text(margin=ggplot2::margin(r=1)) else ggplot2::element_blank(),
      axis.ticks=ggplot2::element_line(linewidth=.25,colour='grey35'),
      axis.ticks.x=if(bottom) ggplot2::element_line() else ggplot2::element_blank(),
      axis.ticks.y=if(left) ggplot2::element_line() else ggplot2::element_blank(),
      axis.ticks.length=grid::unit(.65,'mm'),
      panel.grid=ggplot2::element_blank(),panel.background=ggplot2::element_blank(),
      panel.border=ggplot2::element_rect(colour='grey35',linewidth=.35,fill=NA),
      legend.position='none',plot.margin=ggplot2::margin(0,0,0,0))
}

#' Build a physical colourbar from the map scale
#' @param row One-based diagnostic row (1 or 2).
#' @param name Grid grob identifier.
#' @param scale ggplot2 fill scale already trained for the panel.
#' @param unit_label Optional text/plotmath physical-unit label.
#' @return Grid tree with raster colour strip, vector ticks and unit label; no device opened.
era5_colourbar <- function(scale,row,name,unit_label=NULL) {
  # Exact plot scale mapping, limits, breaks and direction. No data rescaling.
  g <- era5_map_geometry(); h <- g$map_height_mm-12
  limits <- scale$get_limits(); breaks <- scale$get_breaks()
  colours <- scale$map(seq(limits[1],limits[2],length.out=300))
  children <- list(grid::rasterGrob(matrix(rev(colours),ncol=1),
    x=grid::unit(0,'mm'),y=.5,width=grid::unit(2.6,'mm'),height=grid::unit(h,'mm'),
    just=c('left','centre'),interpolate=FALSE,name=paste0(name,'-field')))
  if(!is.null(unit_label)) children <- c(children,list(grid::textGrob(unit_label,
    # Centre the short Kelvin label over the strip itself. The longer compound
    # unit uses the existing legend+gap width, without reducing either map.
    x=grid::unit(if(row==1) 1.3 else 4.7,'mm'),
    y=grid::unit(.5,'npc')+grid::unit(h/2+2,'mm'),just=c('centre','bottom'),
    gp=grid::gpar(fontsize=7,fontfamily='sans'),name=paste0(name,'-unit'))))
  for(i in seq_along(breaks)) {
    y <- grid::unit(.5,'npc')+grid::unit((breaks[i]-mean(limits))/diff(limits)*h,'mm')
    children <- c(children,list(
      grid::segmentsGrob(x0=grid::unit(2.1,'mm'),x1=grid::unit(2.6,'mm'),y0=y,y1=y,
        gp=grid::gpar(col='grey70',lwd=.5)),
      grid::textGrob(scale$get_labels(breaks)[i],x=grid::unit(3.1,'mm'),y=y,just='left',
        gp=grid::gpar(fontsize=7,fontfamily='sans',col='grey25'),name=paste0(name,'-tick-',i))))
  }
  grid::gTree(children=do.call(grid::gList,children),name=name)
}

#' Describe geographic panel positions
#' @return Four-row data frame of row, column and millimetre panel positions/sizes.
era5_panel_geometry <- function() {
  g <- era5_map_geometry(); x <- c(0,cumsum(g$widths)); y <- c(0,cumsum(g$heights))
  data.frame(row=rep(1:2,each=2),column=rep(1:2,2),
    x_mm=rep(x[c(3,5)],2),y_mm=rep(y[g$map_rows],each=2),
    width_mm=g$map_width_mm,height_mm=g$map_height_mm)
}

#' Assemble the accepted four-map geometry
#' @param kind Diagnostic kind: distributional or bias for layouts, hellinger or bias for scales.
#' @param fields List of four numeric longitude x latitude matrices in row-major panel order.
#' @param lon,lat Longitude and latitude cell centres in degrees, matching matrix axes.
#' @return Gtable with era5_layout metadata containing plots, fields, scales and coordinates; no numerical experiment.
assemble_era5_fourmap <- function(fields,lon,lat,kind) {
  stopifnot(length(fields)==4L)
  g <- era5_map_geometry(); stopifnot(g$width_mm==175)
  labels <- era5_row_labels(kind)
  headings <- if(kind=='distributional') c('D-MMM','GC-0.125') else c('MMM','GC-0.125')
  limits <- if(kind=='distributional') list(c(0,.5),c(0,.5)) else list(c(-8,8),c(-4,4))
  scale_kind <- if(kind=='distributional') 'hellinger' else 'bias'
  table <- gtable::gtable(widths=grid::unit(g$widths,'mm'),heights=grid::unit(g$heights,'mm'))
  table <- gtable::gtable_add_grob(table,grid::rectGrob(gp=grid::gpar(fill='white',col=NA)),
    t=1,l=1,b=length(g$heights),r=7,z=-Inf,name='white-background')
  for(column in 1:2) table <- gtable::gtable_add_grob(table,
    era5_column_title_grob(headings[column],paste0('method-',column)),t=g$method_row,l=c(3,5)[column],clip='off')
  plots <- list(); built <- list()
  for(row in 1:2) {
    table <- gtable::gtable_add_grob(table,era5_row_title_grob(labels[[row]],paste0('diagnostic-',row)),
      t=g$title_rows[row],l=3,r=5,clip='off',name=paste0('diagnostic-',row))
    for(column in 1:2) {
      i <- (row-1)*2+column
      plot <- era5_map_panel(fields[[i]],lon,lat,scale_kind,limits[[row]],column==1,row==2)
      plots[[i]] <- plot; built[[i]] <- ggplot2::ggplot_build(plot)
      gt <- ggplot2::ggplotGrob(plot)
      get <- function(name) gt$grobs[[match(name,gt$layout$name)]]
      table <- gtable::gtable_add_grob(table,get('panel'),t=g$map_rows[row],l=c(3,5)[column],
        clip='on',name=paste0('geo-panel-',row,'-',column))
      if(column==1) table <- gtable::gtable_add_grob(table,get('axis-l'),t=g$map_rows[row],l=2,
        clip='off',name=paste0('latitude-axis-',row))
      if(row==2) table <- gtable::gtable_add_grob(table,get('axis-b'),t=10,l=c(3,5)[column],
        clip='off',name=paste0('longitude-axis-',column))
    }
    scale <- built[[(row-1)*2+1]]$plot$scales$get_scales('fill')
    table <- gtable::gtable_add_grob(table,era5_colourbar(scale,row,paste0('colourbar-',row),era5_colourbar_unit(kind,row)),
      t=g$map_rows[row],l=7,clip='off',name=paste0('colourbar-',row))
  }
  attr(table,'era5_layout') <- list(kind=kind,geometry=g,panels=era5_panel_geometry(),
    headings=headings,row_labels=labels,plots=plots,built=built,fields=fields,lon=lon,lat=lat,limits=limits)
  table
}

#' Assemble full-distribution and LDR ERA5 maps
#' @param lon,lat Longitude and latitude cell centres in degrees, matching matrix axes.
#' @param h_mmm,h_gc Longitude x latitude full-distribution Hellinger distance matrices.
#' @param partial_h_mmm,partial_h_gc Longitude x latitude non-renormalised LDR Hellinger score matrices.
#' @param tau_label Numeric-character LDR probability-mass threshold; accepted layout requires 0.01.
#' @return Gtable for full Hellinger distance and non-renormalised LDR Hellinger score; threshold must be 0.01.
build_hdist3_partial_4map_plot <- function(h_mmm,h_gc,partial_h_mmm,partial_h_gc,lon,lat,tau_label) {
  if(as.numeric(tau_label)!=.01) stop('This publication figure requires the unchanged LDR threshold 0.01.')
  assemble_era5_fourmap(list(h_mmm,h_gc,partial_h_mmm,partial_h_gc),lon,lat,'distributional')
}

#' Assemble signed ERA5 mean-bias maps
#' @param lon,lat Longitude and latitude cell centres in degrees, matching matrix axes.
#' @param bias_tas_mmm,bias_tas_gc Signed temperature bias matrices (K), longitude x latitude.
#' @param bias_pr_mmm,bias_pr_gc Signed precipitation bias matrices (mm/day), longitude x latitude.
#' @return Gtable using unchanged temperature [-8,8] K and precipitation [-4,4] mm/day scales.
build_bias_4map_plot <- function(bias_tas_mmm,bias_tas_gc,bias_pr_mmm,bias_pr_gc,lon,lat) {
  assemble_era5_fourmap(list(bias_tas_mmm,bias_tas_gc,bias_pr_mmm,bias_pr_gc),lon,lat,'bias')
}

#' Fingerprint cached scientific arrays and coordinates
#' @param fields List of four numeric longitude x latitude matrices in row-major panel order.
#' @param lon,lat Longitude and latitude cell centres in degrees, matching matrix axes.
#' @return List of dimensions, ranges, missingness and SHA-256 identities; no fields are modified.
era5_field_signature <- function(fields,lon,lat) {
  list(coordinates=list(lon_sha256=digest::digest(lon,algo='sha256'),
                       lat_sha256=digest::digest(lat,algo='sha256')),
    fields=lapply(fields,function(x) list(dim=dim(x),range=range(x,na.rm=TRUE),
      missing=sum(is.na(x)),missing_sha256=digest::digest(is.na(x),algo='sha256'),
      numerical_sha256=digest::digest(x,algo='sha256'))))
}

#' Validate the four-panel layout and numerical colour mapping
#' @param plot Gtable returned by the accepted four-map builder, with era5_layout metadata.
#' @return Audit list of dimensions, labels, colour breaks and scientific signatures; mismatches raise errors. Requires an open device for text metrics.
verify_era5_fourmap <- function(plot) {
  meta <- attr(plot,'era5_layout'); g <- meta$geometry; panel <- meta$panels
  stopifnot(nrow(panel)==4,all(panel$width_mm==78),all(panel$height_mm==g$map_height_mm),
    all(diff(panel$x_mm[c(1,3)])==0),all(diff(panel$x_mm[c(2,4)])==0),
    sum(g$widths)==175,length(grep('^geo-panel',plot$layout$name))==4,
    length(grep('^colourbar-',plot$layout$name))==2,
    length(grep('^latitude-axis-',plot$layout$name))==2,
    length(grep('^longitude-axis-',plot$layout$name))==2)
  text_boxes <- list(); add_text <- function(grob,name) {
    text_boxes[[name]] <<- c(width_mm=grid::convertWidth(grid::grobWidth(grob),'mm',valueOnly=TRUE),
                            height_mm=grid::convertHeight(grid::grobHeight(grob),'mm',valueOnly=TRUE))
  }
  for(i in seq_along(meta$plots)) {
    p <- meta$plots[[i]]; b <- meta$built[[i]]; row <- (i-1)%/%2+1; col <- (i-1)%%2+1
    stopifnot(is.null(p$labels$x),is.null(p$labels$y),
      identical(p$data$value,as.vector(meta$fields[[i]])),
      identical(p$data$lon_wrapped,rep(ifelse(meta$lon>180,meta$lon-360,meta$lon),length(meta$lat))),
      identical(b$layout$panel_params[[1]]$x.range,g$xlim),
      identical(b$layout$panel_params[[1]]$y.range,g$ylim),
      identical(b$plot$scales$get_scales('fill')$get_limits(),meta$limits[[row]]),
      identical(p$coordinates$ratio,1.3))
    reference <- era5_fill_scale(if(meta$kind=='bias') 'bias' else 'hellinger',meta$limits[[row]])
    stopifnot(identical(b$data[[1]]$fill,reference$map(as.vector(meta$fields[[i]]))))
    if(col==2) stopifnot(inherits(p$theme$axis.text.y,'element_blank'),inherits(p$theme$axis.ticks.y,'element_blank'))
    if(row==1) stopifnot(inherits(p$theme$axis.text.x,'element_blank'),inherits(p$theme$axis.ticks.x,'element_blank'))
  }
  for(row in 1:2) {
    label <- era5_row_title_grob(meta$row_labels[[row]])
    add_text(label,paste0('row',row))
    stopifnot(tail(text_boxes,1)[[1]][['width_mm']] < 157.5,
      tail(text_boxes,1)[[1]][['height_mm']] <= g$heights[g$title_rows[row]])
    bar <- era5_colourbar(meta$built[[(row-1)*2+1]]$plot$scales$get_scales('fill'),row,'check',
      era5_colourbar_unit(meta$kind,row))
    ticks <- Filter(function(x) inherits(x,'text'),as.list(bar$children))
    for(i in seq_along(ticks)) {
      add_text(ticks[[i]],paste0('bar',row,'_text',i))
      width <- tail(text_boxes,1)[[1]][['width_mm']]
      if(grepl('-unit$',ticks[[i]]$name)) {
        centre <- grid::convertX(ticks[[i]]$x,'mm',valueOnly=TRUE)
        stopifnot(centre-width/2 >= -g$widths[['colourbar_gap']],
          centre+width/2 <= g$widths[['colourbar']])
      } else stopifnot(width<=g$widths[['colourbar']])
    }
  }
  for(col in 1:2) {
    heading <- era5_column_title_grob(meta$headings[col]);add_text(heading,paste0('column',col))
    stopifnot(tail(text_boxes,1)[[1]][['height_mm']]<=g$heights[['column_heading']])
  }
  list(status='passed',kind=meta$kind,dimensions_mm=c(width=g$width_mm,height=g$height_mm),
    geographic_panels=panel,column_heading_centres_mm=panel$x_mm[1:2]+39,
    column_heading_y_mm=cumsum(g$heights)[g$method_row]-g$heights[g$method_row]/2,
    row_heading_centres_y_mm=cumsum(g$heights)[g$title_rows]-g$heights[g$title_rows]/2,
    row_heading_x_mm=mean(panel$x_mm[1:2]+39),
    colourbar_centres_y_mm=panel$y_mm[c(1,3)]+g$map_height_mm/2,
    row_label_gutter_mm=0,axis_text_pt=7,row_text_pt=9,column_text_pt=10,
    text_bounds_mm=text_boxes,colour_limits=meta$limits,
    colour_breaks=lapply(c(1,3),function(i) meta$built[[i]]$plot$scales$get_scales('fill')$get_breaks()),
    method_order=meta$headings,scientific_signature=era5_field_signature(meta$fields,meta$lon,meta$lat))
}

#' Save and verify accepted ERA5 map layouts
#' @param plot Gtable returned by the accepted four-map builder, with era5_layout metadata.
#' @param base Output filename stem without extension; parent directory must exist.
#' @param cache_path Existing compact RDS input whose bytes must remain unchanged.
#' @return Audit list invisibly; writes PDF, PNG and verification JSON at base, overwriting prior files. Requires Cairo and ragg; verifies cache hash remains unchanged.
save_era5_fourmap <- function(plot,base,cache_path) {
  original_file_hash <- digest::digest(file=cache_path,algo='sha256')
  signature <- with(attr(plot,'era5_layout'),era5_field_signature(fields,lon,lat))
  audits <- list()
  for(format in c('pdf','png')) {
    path <- paste0(base,'.',format)
    if(format=='pdf') grDevices::cairo_pdf(path,width=a4_plot_width_cm/2.54,
      height=fourmap_total_height_cm/2.54,family='sans',bg='white',onefile=FALSE)
    else ragg::agg_png(path,width=175,height=era5_map_geometry()$height_mm,units='mm',res=300,background='white')
    audits[[format]] <- tryCatch({
      grid::grid.newpage();grid::grid.draw(plot);verify_era5_fourmap(plot)
    },finally=grDevices::dev.off())
  }
  stopifnot(identical(original_file_hash,digest::digest(file=cache_path,algo='sha256')),
    identical(signature,with(attr(plot,'era5_layout'),era5_field_signature(fields,lon,lat))))
  audits$cache <- list(path=cache_path,sha256_before=original_file_hash,sha256_after=original_file_hash)
  jsonlite::write_json(audits,paste0(base,'_verification.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
  invisible(audits)
}

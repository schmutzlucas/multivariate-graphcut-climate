# Use an external temporary device for layout metrics before explicit export.
options(device = function(...) grDevices::pdf(file = tempfile("mvgc-layout-", fileext=".pdf")))
source("R/load.R"); load_mvgc()
source("R/plotting/era5_maps.R")
x <- readRDS("results/era5/distributional_fields.rds")
plots <- lapply(c("raster", "tile"), function(mode)
  era5_map_panel(x$h_gc,x$lon,x$lat,"hellinger",c(0,.5),TRUE,TRUE,mode))
stopifnot(inherits(plots[[1]]$layers[[1]]$geom,"GeomRaster"),
          inherits(plots[[2]]$layers[[1]]$geom,"GeomTile"),
          identical(plots[[1]]$data,plots[[2]]$data),
          identical(plots[[1]]$labels,plots[[2]]$labels))
b <- lapply(plots,ggplot2::ggplot_build)
stopifnot(identical(b[[1]]$data[[1]]$fill,b[[2]]$data[[1]]$fill),
          identical(b[[1]]$layout$panel_params[[1]]$x.range,b[[2]]$layout$panel_params[[1]]$x.range),
          identical(b[[1]]$layout$panel_params[[1]]$y.range,b[[2]]$layout$panel_params[[1]]$y.range))
p <- build_hdist3_partial_4map_plot(x$h_d_mmm,x$h_gc,x$ldr_d_mmm,x$ldr_gc,x$lon,x$lat,"0.01")
meta <- attr(p,"era5_layout")
stopifnot(all(meta$panels$width_mm==78),meta$geometry$width_mm==175,
          identical(meta$headings,c("D-MMM","GC-0.125")),
          identical(meta$row_labels,era5_row_labels("distributional")))
cat("PASS explicit primitives, unchanged arrays/colours/extents, accepted panel sizes/titles.\n")

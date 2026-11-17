# Fig. 6 candidates. First run experiments/sensitivity/spatial_regularisation.R.
# Every panel displays an actual graph-cut solution, with identical model colours.
source("R/load.R")
load_mvgc()
source("R/utils/reproduction.R")
require_packages(c("ggplot2", "maps"))
suppressPackageStartupMessages(library(ggplot2))
out <- mvgc_output_dir("figures/lambda")
bundle_file <- reproduction_input("label_maps")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(bundle_file)) {
  b <- readRDS(bundle_file)
} else {
  # Allow review of a complete candidate while other independent solutions run.
  checkpoint_root <- Sys.getenv("MVGC_LABEL_CHECKPOINT_DIR", unset = "")
  if (!nzchar(checkpoint_root)) {
    checkpoint_root <- file.path(mvgc_output_dir("spatial_regularisation"), "checkpoints")
    identities <- list.dirs(checkpoint_root, recursive = FALSE, full.names = TRUE)
    identities <- identities[grepl("^identity_", basename(identities))]
    if (length(identities) > 1L) stop("Several spatial input identities exist; select MVGC_LABEL_CHECKPOINT_DIR explicitly.")
    if (length(identities) == 1L) checkpoint_root <- identities
  }
  files <- list.files(checkpoint_root, pattern = "^lambda_.*\\.rds$", full.names = TRUE)
  if (!length(files)) stop("Compute the lambda label maps first: ", bundle_file)
  solutions <- lapply(files, readRDS)
  solutions <- solutions[order(vapply(solutions, function(s) s$lambda, numeric(1)))]
  names(solutions) <- vapply(solutions, function(s) format(s$lambda, scientific = FALSE, trim = TRUE), "")
  signature <- solutions[[1]]$input_signature
  stopifnot(all(vapply(solutions, function(s) identical(s$input_signature, signature), logical(1))))
  b <- list(solutions = solutions, lon = -180:179, lat = -90:90,
            model_names = signature$models, input_signature = signature)
}
stopifnot(identical(b$lon, -180:179), identical(b$lat, -90:90),
          is.character(b$model_names), length(b$model_names) == 22L,
          !anyNA(b$model_names), all(nzchar(b$model_names)), !anyDuplicated(b$model_names),
          is.list(b$solutions), length(b$solutions) > 0L,
          !is.null(names(b$solutions)), !anyDuplicated(names(b$solutions)))
for (i in seq_along(b$solutions)) {
  s <- b$solutions[[i]]
  stopifnot(is.numeric(s$lambda), length(s$lambda) == 1L, is.finite(s$lambda), s$lambda >= 0,
            identical(names(b$solutions)[i], format(s$lambda, scientific = FALSE, trim = TRUE)),
            identical(s$input_signature, b$input_signature))
  if (!is.null(s$state) && !identical(s$state, "complete")) stop("Cannot render an incomplete spatial solution.")
}

# Same model order and categorical palette as ERA5_experiment_figures.R.
palette <- setNames(c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#5254a3"), b$model_names)
#' Format a numeric lambda for lookup in a compact label-map bundle
#' @param x Scalar lambda already validated by the bundle/candidate-set checks.
#' @return Decimal character list key with whitespace trimmed and no scientific notation.
key <- function(x) format(x, scientific = FALSE, trim = TRUE)
#' Build a concise lambda facet title with the original manuscript formatting
#' @param x Numeric lambda; values 1 and 2 retain a trailing decimal zero.
#' @return Scalar character title with a lambda glyph and decimal value.
lambda_title <- function(x) {
  value <- if (x == 1 || x == 2) sprintf("%.1f", x) else key(x)
  paste0("\u03bb = ", value)
}

# Count each unordered four-connected edge once. Longitude wraps; latitude
# does not. This is an unweighted edge fraction, NOT the Hellinger/latitude-
# weighted smoothness energy. Polar grid rows are retained as in gcoWrapR.
#' Count differently labelled edges on the graph-cut cylindrical grid
#' @param labels Numeric matrix of positive finite integer labels indexed longitude
#'   x latitude, with at least three longitude and two latitude cells.
#' @return One-row data frame with longitude/latitude/different/total edge counts
#'   and fragmentation_pct. Each four-connected edge is counted once; longitude
#'   wraps, latitude does not, polar rows remain, and no area/energy weighting occurs.
fragmentation <- function(labels) {
  stopifnot(is.matrix(labels), is.numeric(labels), all(is.finite(labels)),
            all(labels >= 1 & labels == floor(labels)), nrow(labels) >= 3L, ncol(labels) >= 2L)
  nx <- nrow(labels); ny <- ncol(labels)
  longitude_different <- sum(labels != labels[c(2:nx, 1L), , drop = FALSE])
  latitude_different <- sum(labels[, 1:(ny - 1L), drop = FALSE] !=
                              labels[, 2:ny, drop = FALSE])
  edges <- nx * ny + nx * (ny - 1L)
  data.frame(longitude_edges = nx * ny, latitude_edges = nx * (ny - 1L),
             longitude_different = longitude_different, latitude_different = latitude_different,
             different_edges = longitude_different + latitude_different,
             total_edges = edges,
             fragmentation_pct = 100 * (longitude_different + latitude_different) / edges)
}
# Small independent checks catch boundary/double-counting mistakes.
stopifnot(fragmentation(matrix(1L, 3, 2))$fragmentation_pct == 0)
toy <- cbind(c(1L, 1L, 2L), c(1L, 2L, 2L))
stopifnot(fragmentation(toy)$different_edges == 5,
          fragmentation(toy)$total_edges == 9)
stats <- do.call(rbind, lapply(b$solutions, function(s) {
  labs <- s$result$label_attribution
  stopifnot(identical(dim(labs), c(360L, 181L)),
            all(labs %in% seq_along(b$model_names)))
  cbind(data.frame(lambda = s$lambda, origin = s$origin,
                   selected_models = length(unique(as.vector(labs)))), fragmentation(labs))
}))
rownames(stats) <- NULL
stopifnot(all(stats$total_edges == 129960))
sets <- list(lambda_label_maps_balanced = c(0, .025, .125, .25, .5, 1),
             lambda_label_maps_balanced_3x2 = c(0, .025, .125, .25, .5, 1),
             lambda_label_maps_balanced_fragmentation = c(0, .025, .125, .25, .5, 1),
             lambda_label_maps_extreme = c(0, .025, .125, .25, 1, 2),
             lambda_label_maps_transition_zoom = c(0, .025, .0625, .125, .25, .5))
balanced_values <- sets[[1]]
selected <- Sys.getenv("MVGC_LAMBDA_FIGURES", unset = "lambda_label_maps_balanced_3x2")
if (nzchar(selected)) {
  selected <- strsplit(selected, ",", fixed = TRUE)[[1]]
  stopifnot(all(selected %in% names(sets)))
  sets <- sets[selected]
}
if (!all(unlist(sets) %in% stats$lambda)) stop("Requested candidate still has missing lambda solutions.")
if (all(balanced_values %in% stats$lambda)) {
  write.csv(stats[match(balanced_values, stats$lambda), ],
            file.path(out, "lambda_label_maps_balanced_fragmentation.csv"), row.names = FALSE)
}
write.csv(stats, file.path(out, "lambda_label_maps_all_statistics.csv"), row.names = FALSE)
write.csv(data.frame(model_index = seq_along(palette), model = names(palette), colour = palette),
          file.path(out, "lambda_label_maps_model_palette.csv"), row.names = FALSE)

# Draw coastlines, excluding internal country boundaries directly in the maps
# database. NA separators prevent connecting distinct coast segments/islands.
world <- maps::map("world", plot = FALSE, interior = FALSE)
coast_id <- cumsum(is.na(world$x) | is.na(world$y) |
                    world$x < -180.5 | world$x > 179.5)
coast <- data.frame(lon = world$x, lat = world$y, group = coast_id)
coast <- coast[is.finite(coast$lon) & is.finite(coast$lat) &
                 coast$lon >= -180.5 & coast$lon <= 179.5, ]
rownames(coast) <- NULL

caption <- paste0(
  "Effect of the smoothness parameter $\\lambda$ on graph-cut model selection in the ERA5 experiment ",
  "(calibration: 1951--1975; projection: 1999--2023). Increasing $\\lambda$ favours spatially coherent ",
  "model patches over independent local selection. The main configuration ($\\lambda=0.125$) and ",
  "stronger-regularisation sensitivity case ($\\lambda=0.25$) are included. Colours identify the same ",
  "22 CMIP6 models in every panel; only $\\lambda$ changes between solutions. The grid uses ",
  "four-connected neighbours with periodic longitude and no additional polar connection.")
fragment_caption <- paste0(" Panel annotations give the percentage of neighbouring grid-cell pairs ",
  "assigned different models, counting each edge once (129,960 edges) without area weighting; ",
  "this fragmentation diagnostic is distinct from the weighted smoothness energy.")

for (filename in names(sets)) {
  values <- sets[[filename]]
  show_stat <- grepl("fragmentation$", filename)
  balanced_maps_only <- filename %in% c("lambda_label_maps_balanced", "lambda_label_maps_balanced_3x2")
  panel_columns <- if (filename == "lambda_label_maps_balanced_3x2") 2L else 3L
  # Compact only the portrait layout to leave room for the manuscript caption.
  # Keep the landscape balanced layout at its earlier manuscript aspect.
  map_ratio <- if (panel_columns == 2L) 1.2 else if (balanced_maps_only) 1.3 else 1
  frames <- lapply(values, function(x) {
    labs <- b$solutions[[key(x)]]$result$label_attribution
    data.frame(lon = rep(b$lon, length(b$lat)), lat = rep(b$lat, each = length(b$lon)),
               model = factor(b$model_names[as.vector(labs)], levels = b$model_names),
               panel = lambda_title(x))
  })
  df <- do.call(rbind, frames)
  df$panel <- factor(df$panel, levels = vapply(values, lambda_title, ""))
  annotations <- data.frame(panel = levels(df$panel), lon = 0, lat = -103,
    label = sprintf("Different neighbours: %.1f%%", stats$fragmentation_pct[match(values, stats$lambda)]))
  annotations$panel <- factor(annotations$panel, levels = levels(df$panel))
  p <- ggplot(df, aes(lon, lat)) +
    geom_raster(aes(fill = model), interpolate = FALSE) +
    geom_path(data = coast, aes(lon, lat, group = group), inherit.aes = FALSE,
              colour = "#303030", linewidth = .12) +
    facet_wrap(vars(panel), ncol = panel_columns) +
    scale_fill_manual(values = palette, limits = b$model_names, drop = FALSE, name = NULL) +
    coord_fixed(ratio = map_ratio, xlim = c(-180.5, 179.5), ylim = c(-90.5, 90.5),
                 expand = FALSE, clip = if (show_stat) "off" else "on") +
    guides(fill = "none") +
    theme_void(base_size = 9, base_family = "sans") +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "#999999", linewidth = .25),
          panel.spacing.x = grid::unit(if (balanced_maps_only) 4 else 2.5, "mm"),
          panel.spacing.y = grid::unit(if (panel_columns == 2L) 2 else if (balanced_maps_only) 5 else 3, "mm"),
          strip.text = element_text(size = 10, colour = "#303030", margin = margin(0, 0, 3, 0)),
          legend.position = "none",
          plot.margin = margin(5, 5, if (show_stat) 14 else 5, 5))
  if (show_stat) {
    # Keep every map cell visible: place the statistic below, not over Antarctica.
    p <- p + geom_text(data = annotations, aes(lon, lat, label = label),
                        inherit.aes = FALSE, hjust = .5, size = 2.45, colour = "#555555") +
      theme(panel.spacing.y = grid::unit(6, "mm"))
  }
  # Verify exact panel populations, identity mapping and row-major order.
  built <- ggplot_build(p)
  stopifnot(nrow(built$layout$layout) == 6L,
            max(built$layout$layout$COL) == panel_columns,
            max(built$layout$layout$ROW) == 6L / panel_columns,
            all(table(built$data[[1]]$PANEL) == 65160L),
            identical(as.character(built$layout$layout$panel), levels(df$panel)),
            identical(unname(built$plot$scales$get_scales("fill")$map(b$model_names)), unname(palette)))
  base <- file.path(out, if (filename == "lambda_label_maps_balanced_3x2") "lambda_label_maps" else filename)
  figure_height <- if (panel_columns == 2L) 18.6 else if (balanced_maps_only) 9.6 else if (show_stat) 8.3 else 7.6
  ggsave(paste0(base, ".png"), p, width = 19, height = figure_height, units = "cm", dpi = 600, bg = "white")
  ggsave(paste0(base, ".pdf"), p, width = 19, height = figure_height, units = "cm",
         device = cairo_pdf, fallback_resolution = 600, bg = "white")
  writeLines(paste0("\\caption{", caption, if (show_stat) fragment_caption else "", "}"),
             paste0(base, "_caption.tex"))
  cat("Verified and saved", filename, "\n")
}
note <- c("Fig. 6 candidate provenance", "",
  paste("Source:", b$input_signature$workspace),
  "ERA5 calibration 1951-1975; projection 1999-2023; full year; pr, tas, psl; 8^3 bins.",
  "Data cost uses calibration Hellinger distances; smoothness uses projection-period model PMFs.",
  "Unchanged settings: weight_data=1, seed=1, original saved optimizer defaults, swap(-1).",
  paste("gcoWrapR commit:", b$input_signature$gco_commit),
  "Longitude-periodic, four-connected grid; latitude/polar wrap absent.",
  "Loaded existing final-period map: lambda=0.125.",
  paste("Newly computed actual maps: lambda=", paste(stats$lambda[stats$origin == "new graph-cut computation"], collapse = ", ")),
  "Legacy Fig. 6 PNGs are not numerical label-map inputs and were not interpolated.",
  "Colours are the unchanged ERA5 categorical palette, indexed by the saved model_names.",
  "Fragmentation counts each of 129,960 unique grid edges once; no latitude/area weighting.",
  "Candidates only: main.tex and the current manuscript figure have not been replaced.", "",
  "Reproduce from the code-repository root:",
  "Set ERA5_RESULTS_FILE to the saved final ERA5 workspace, then run:",
  "Rscript experiments/sensitivity/spatial_regularisation.R", "Rscript figures/spatial_regularisation.R")
writeLines(note, file.path(out, "lambda_label_maps_provenance.txt"))
print(stats[, c("lambda", "fragmentation_pct", "selected_models", "origin")], row.names = FALSE)

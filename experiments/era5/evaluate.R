# Evaluation reads saved production results; it does not refit labels.
# Missing/incompatible source objects fail explicitly. Outputs are derived scores
# under the configured output root and may overwrite existing evaluation files.
source("R/load.R")
load_mvgc()
# Global summary metrics for the ERA5 GraphCut experiment.
# This script is table-only: it does not generate or regenerate figures.

workspace_file <- production_result_file("ERA5_RESULTS_FILE", "the final ERA5 RData workspace")

if (!nzchar(workspace_file)) {
  candidates <- list.files(
    "workspaces",
    pattern = "_workspace_multiRes_3v_gc_results_ERA5_fullyear_cal1951-1975_proj1999-2023\\.RData$",
    full.names = TRUE
  )
  if (!length(candidates)) {
    stop("No final ERA5 workspace found. Run experiments/era5/run.R first.")
  }
  workspace_file <- candidates[which.max(file.info(candidates)$mtime)]
}

if (!file.exists(workspace_file)) {
  stop("ERA5 workspace not found: ", workspace_file)
}

source("R/metrics/ldr.R")
source("R/metrics/partial_hellinger.R")

required_objects <- c(
  "lon", "lat", "year_present", "year_future", "model_names", "GC_result", "mean_future",
  "pdf1_future", "pdf2_future", "pdf3_ref_fut", "pdf3_models_fut",
  "pdf3_MMM_fut", "pdf3_GC_fut", "MMM_hdist_fut", "GC_hdist_fut"
)
# A workspace must supply its own scientific inputs: an interactive session must
# not fill missing fields, and archived helper definitions must not replace code.
workspace <- new.env(parent = emptyenv())
load(workspace_file, envir = workspace)
missing_objects <- setdiff(required_objects, ls(workspace, all.names = TRUE))
if (length(missing_objects) > 0L) {
  stop("ERA5 workspace is missing required object(s): ", paste(missing_objects, collapse = ", "))
}
available_objects <- c(required_objects, intersect("smooth_cost", ls(workspace)))
list2env(mget(available_objects, envir = workspace, inherits = FALSE), envir = environment())
rm(workspace)
cfg <- article1_config()
if (!is.numeric(year_present) || !is.numeric(year_future) ||
    !identical(as.numeric(year_present), as.numeric(cfg$calibration_years)) ||
    !identical(as.numeric(year_future), as.numeric(cfg$era5_target_years))) {
  stop("ERA5 summary requires calibration 1951-1975 and projection 1999-2023.")
}
if (!is.numeric(lon) || !length(lon) || any(!is.finite(lon)) || anyDuplicated(lon) ||
    !is.numeric(lat) || !length(lat) || any(!is.finite(lat)) || any(abs(lat) > 90) || anyDuplicated(lat)) {
  stop("ERA5 summary requires finite distinct longitude/latitude coordinates and latitudes in [-90,90].")
}
if (!is.character(model_names) || !length(model_names) || anyNA(model_names) ||
    any(!nzchar(model_names)) || anyDuplicated(model_names)) {
  stop("ERA5 summary requires distinct nonempty candidate model names.")
}
# New whitelisted workspaces intentionally do not save incidental summary weights.
# Reconstruct the published equal-longitude cosine-latitude convention explicitly.
area_weight_matrix <- matrix(rep(cos(lat * pi / 180), each = length(lon)), nrow = length(lon))

tau_values <- cfg$era5_ldr_thresholds
method_names <- c("MMM", "GC-0.125")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- mvgc_output_dir("evaluation/era5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#' Compute ordinary Hellinger distance for two valid PMF vectors
#' @param pdf_ref,pdf_model Numeric normalised PMFs on identical ordered bins.
#' @return Scalar Hellinger distance in [0,1]. Invalid masses raise an error;
#'   they are never clipped or replaced with zeros.
hellinger_vector <- function(pdf_ref, pdf_model) {
  hellinger_distance(pdf_ref, pdf_model)
}

#' Remove only a singleton model axis from a Hellinger map
#' @param x Numeric longitude x latitude matrix or longitude x latitude x 1 array.
#' @return Numeric matrix retaining both spatial dimensions, including singleton
#'   spatial axes. Missing scores remain missing. Other shapes raise an error.
coerce_h_map <- function(x) {
  dims <- dim(x)
  if (!is.numeric(x) || !(length(dims) == 2L || (length(dims) == 3L && dims[3] == 1L))) {
    stop("Expected a lon x lat map with at most a singleton model axis.")
  }
  matrix(x, nrow = dims[1], ncol = dims[2])
}

#' Select finite positive-weight samples without recycling mismatched fields
#' @param field Numeric field, usually a longitude x latitude matrix.
#' @param weights Numeric weights with exactly the same length and dimensions.
#'   Nonfinite or zero weights are excluded; finite negative weights are invalid.
#' @param transform Function returning one numeric transformed value per input cell.
#' @return List with valid transformed values, their weights and a logical valid
#'   mask in column-major order. Missing field cells are excluded with their weights.
weighted_values <- function(field, weights, transform = identity) {
  if (!is.numeric(field) || !is.numeric(weights) || length(field) != length(weights) ||
      !identical(dim(field), dim(weights)) || any(weights[is.finite(weights)] < 0)) {
    stop("Field and nonnegative weights must have identical dimensions and length.")
  }
  values <- transform(as.numeric(field))
  if (!is.numeric(values) || length(values) != length(field)) {
    stop("transform must return one numeric value per field cell.")
  }
  w <- as.numeric(weights)
  ok <- is.finite(values) & is.finite(w) & w > 0
  list(values = values[ok], weights = w[ok], valid = ok)
}

#' Compute a finite-cell area-weighted mean
#' @param field,weights Numeric field and nonnegative weights with identical shape.
#' @param transform Elementwise function applied before selecting finite values.
#' @return Scalar weighted mean or NA_real_ if no positive valid weight remains.
weighted_mean_value <- function(field, weights, transform = identity) {
  xw <- weighted_values(field, weights, transform)
  if (length(xw$values) == 0L) {
    return(NA_real_)
  }
  sum(xw$values * xw$weights) / sum(xw$weights)
}

#' Compute weighted empirical quantiles without interpolation
#' @param field,weights Numeric field and nonnegative weights with identical shape.
#' @param probs Finite numeric probabilities in [0,1], in the desired output order.
#' @param transform Elementwise function applied before selecting finite values.
#' @return Numeric vector, the first sorted value whose cumulative normalised
#'   weight reaches each probability. Returns NA entries when no valid weight remains.
weighted_quantile_value <- function(field, weights, probs, transform = identity) {
  if (!is.numeric(probs) || any(!is.finite(probs)) || any(probs < 0 | probs > 1)) {
    stop("probs must be finite probabilities in [0,1].")
  }
  xw <- weighted_values(field, weights, transform)
  if (length(xw$values) == 0L) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(xw$values)
  values <- xw$values[ord]
  w <- xw$weights[ord]
  cw <- cumsum(w) / sum(w)
  cw[length(cw)] <- 1 # Protect the exact 100% quantile against summation roundoff.
  vapply(probs, function(p) values[which(cw >= p)[1]], numeric(1))
}

#' Compute the weighted population standard deviation of a field
#' @param field,weights Numeric field and nonnegative weights with identical shape.
#' @param transform Elementwise function applied before selecting finite values.
#' @return Scalar population SD (no sample correction), or NA_real_ with no valid area.
weighted_sd_value <- function(field, weights, transform = identity) {
  xw <- weighted_values(field, weights, transform)
  if (length(xw$values) == 0L) {
    return(NA_real_)
  }
  mu <- sum(xw$values * xw$weights) / sum(xw$weights)
  sqrt(sum(xw$weights * (xw$values - mu) ^ 2) / sum(xw$weights))
}

#' Measure the fraction of positive-weight area with finite field values
#' @param field,weights Numeric field and nonnegative weights with identical shape.
#' @return Scalar in [0,1], or NA_real_ if total positive finite weight is zero.
valid_area_fraction <- function(field, weights) {
  weighted_values(field, weights)
  ok_weight <- is.finite(weights) & weights > 0
  ok <- is.finite(field) & ok_weight
  if (!any(ok_weight)) {
    return(NA_real_)
  }
  sum(weights[ok]) / sum(weights[ok_weight])
}

#' Compute the area-weighted root mean squared error of a bias field
#' @param field Numeric signed bias field.
#' @param weights Nonnegative weights with identical shape to field.
#' @return Scalar sqrt(weighted mean(field^2)) on finite cells, or NA_real_.
weighted_rmse_against_zero <- function(field, weights) {
  sqrt(weighted_mean_value(field ^ 2, weights))
}

#' Compute the weighted Pearson spatial correlation on common finite cells
#' @param field,ref Numeric spatial fields with identical shape.
#' @param weights Nonnegative weights of the same shape.
#' @return Scalar correlation; NA_real_ with no common area or either field constant.
weighted_spatial_correlation <- function(field, ref, weights) {
  weighted_values(field, weights)
  weighted_values(ref, weights)
  ok <- is.finite(field) & is.finite(ref) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  x <- field[ok]
  y <- ref[ok]
  w <- weights[ok] / sum(weights[ok])
  mx <- sum(w * x)
  my <- sum(w * y)
  vx <- sum(w * (x - mx) ^ 2)
  vy <- sum(w * (y - my) ^ 2)
  if (vx <= 0 || vy <= 0) {
    return(NA_real_)
  }
  sum(w * (x - mx) * (y - my)) / sqrt(vx * vy)
}

#' Measure the area with a strict reduction in absolute error
#' @param err_gc,err_mmm Numeric GC and D-MMM errors with identical shape.
#' @param weights Nonnegative weights of the same shape.
#' @return Weighted fraction abs(err_gc) < abs(err_mmm) on common finite cells,
#'   or NA_real_ when none remain. Exact ties are not improvements.
weighted_fraction_improved <- function(err_gc, err_mmm, weights) {
  weighted_values(err_gc, weights)
  weighted_values(err_mmm, weights)
  ok <- is.finite(err_gc) & is.finite(err_mmm) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(weights[ok] * (abs(err_gc[ok]) < abs(err_mmm[ok]))) / sum(weights[ok])
}

#' Summarise a spatial diagnostic with the published area-weighting convention
#' @param field Numeric field, typically a Hellinger score map.
#' @param weights Nonnegative weights with identical shape.
#' @return One-row data frame containing mean, median, population SD, IQR,
#'   90th/95th percentiles and finite-area fraction.
summarize_map <- function(field, weights) {
  qs <- weighted_quantile_value(field, weights, probs = c(0.25, 0.50, 0.75, 0.90, 0.95))
  data.frame(
    weighted_mean = weighted_mean_value(field, weights),
    weighted_median = qs[2],
    weighted_sd = weighted_sd_value(field, weights),
    weighted_iqr = qs[3] - qs[1],
    weighted_p90 = qs[4],
    weighted_p95 = qs[5],
    valid_area_fraction = valid_area_fraction(field, weights),
    stringsAsFactors = FALSE
  )
}

#' Append non-NULL table fragments in order
#' @param target Existing data frame, empty data frame or NULL.
#' @param rows List of conformable data frames; NULL elements are omitted.
#' @return Row-bound data frame, or target unchanged when no fragments remain.
append_rows <- function(target, rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(target)
  }
  if (is.null(target) || nrow(target) == 0L) {
    return(do.call(rbind, rows))
  }
  rbind(target, do.call(rbind, rows))
}

#' Select every support bin at every reference grid cell
#' @param pdf_ref Reference PMF array indexed longitude x latitude x bins.
#' @return Nested longitude/latitude lists of integer bin sequences for the
#'   ordinary Hellinger calculation. Caller supplies a valid PMF array.
build_all_indices <- function(pdf_ref) {
  lon_size <- dim(pdf_ref)[1]
  lat_size <- dim(pdf_ref)[2]
  n_bins <- dim(pdf_ref)[3]
  lapply(seq_len(lon_size), function(i) {
    lapply(seq_len(lat_size), function(j) seq_len(n_bins))
  })
}

#' Select reference-defined LDR bins for every grid cell
#' @param pdf_ref Reference PMF array indexed longitude x latitude x bins.
#' @param tau LDR cumulative-mass threshold in (0,1), called ell_LDR in the paper.
#' @return Nested longitude/latitude lists following select_ldr_indices(),
#'   including zero-reference bins and its complete cutoff-tie rule.
build_ldr_indices <- function(pdf_ref, tau) {
  lon_size <- dim(pdf_ref)[1]
  lat_size <- dim(pdf_ref)[2]
  lapply(seq_len(lon_size), function(i) {
    lapply(seq_len(lat_size), function(j) select_ldr_indices(pdf_ref[i, j, ], tau = tau))
  })
}

#' Evaluate D-MMM and GC PMFs on identical selected bins
#' @param pdf_ref,pdf_mmm,pdf_gc Valid PMF arrays indexed longitude x latitude x bins.
#' @param selected_indices Nested longitude/latitude bin-index lists.
#' @return Named list MMM and GC-0.125 of two-dimensional score maps. Restricted
#'   masses are not renormalised; empty selections retain missing scores.
compute_h_maps_from_pdfs <- function(pdf_ref, pdf_mmm, pdf_gc, selected_indices) {
  list(
    MMM = coerce_h_map(compute_partial_hdist(pdf_ref, pdf_mmm, selected_indices)),
    `GC-0.125` = coerce_h_map(compute_partial_hdist(pdf_ref, pdf_gc, selected_indices))
  )
}

#' Evaluate every candidate model against reference-selected bins
#' @param pdf_ref Valid reference PMF array longitude x latitude x bins.
#' @param pdf_models Valid candidate array longitude x latitude x bins x models.
#' @param selected_indices Nested longitude/latitude bin-index lists.
#' @return Longitude x latitude x models score array, retaining missing scores
#'   and a singleton model dimension. Does not normalise restricted PMFs.
compute_model_h_maps_from_pdfs <- function(pdf_ref, pdf_models, selected_indices) {
  compute_partial_hdist(pdf_ref, pdf_models, selected_indices)
}

#' Extract the spatially selected candidate PMF from a reference-first array
#' @param pdf_full Numeric array longitude x latitude x bins x models, with
#'   reference at index 1 and at least one candidate afterwards.
#' @param gc_labels Finite integer longitude x latitude matrix indexing candidates
#'   from 1, so candidate label k selects full model index k+1.
#' @return Longitude x latitude x bins array, retaining singleton spatial axes.
#'   Invalid shapes, missing labels or labels outside candidate range error.
extract_gc_pdf_from_full_array <- function(pdf_full, gc_labels) {
  dims <- dim(pdf_full)
  if (!is.numeric(pdf_full) || length(dims) != 4L || dims[4] < 2L) {
    stop("Expected full PDF array with dimensions lon x lat x bins x model.")
  }
  n_lon <- dims[1]
  n_lat <- dims[2]
  n_bins <- dims[3]
  n_models_excluding_ref <- dims[4] - 1L
  if (!is.numeric(gc_labels) || !identical(dim(gc_labels), c(n_lon, n_lat))) {
    stop("GraphCut label dimensions do not match PDF grid.")
  }
  if (any(!is.finite(gc_labels)) || any(gc_labels != floor(gc_labels)) ||
      any(gc_labels < 1L | gc_labels > n_models_excluding_ref)) {
    stop("GraphCut labels are outside the expected model index range.")
  }

  pdf_gc <- array(NA_real_, dim = c(n_lon, n_lat, n_bins))
  for (label in seq_len(n_models_excluding_ref)) {
    mask <- gc_labels == label
    if (!any(mask, na.rm = TRUE)) {
      next
    }
    full_index <- label + 1L
    for (b in seq_len(n_bins)) {
      target <- matrix(pdf_gc[, , b], nrow = n_lon, ncol = n_lat)
      source <- matrix(pdf_full[, , b, full_index], nrow = n_lon, ncol = n_lat)
      target[mask] <- source[mask]
      pdf_gc[, , b] <- target
    }
  }
  pdf_gc
}

#' Build an equal-model distributional mean without filling missing PMF masses
#' @param pdf_full Numeric longitude x latitude x bins x models array, reference
#'   first and at least one candidate model afterwards.
#' @return Longitude x latitude x bins candidate-mean PMF array. Missing values
#'   propagate: bin-wise available-model mixing could produce a nonnormalised PMF.
compute_mmm_pdf_from_full_array <- function(pdf_full) {
  if (!is.numeric(pdf_full) || length(dim(pdf_full)) != 4L || dim(pdf_full)[4] < 2L) {
    stop("Expected a reference-first PMF array with at least one candidate model.")
  }
  apply(pdf_full[, , , -1L, drop = FALSE], c(1, 2, 3), mean)
}

#' Create one spatial distance-summary row with archived column identities
#' @param method Method identifier; MMM internally denotes D-MMM.
#' @param dimension,variable_set Dimension and variable-set labels.
#' @param metric_type Archived metric identifier, e.g. full_H or partial_LDR_H.
#' @param tau LDR threshold; NULL for ordinary Hellinger.
#' @param field,weights Numeric score map and nonnegative weights of identical shape.
#' @return One-row data frame of identifiers, H units and summarize_map() statistics.
make_distance_row <- function(method, dimension, variable_set, metric_type, tau, field, weights) {
  cbind(
    data.frame(
      method = method,
      dimension = dimension,
      variable_set = variable_set,
      metric_type = metric_type,
      tau = if (is.null(tau)) NA_real_ else tau,
      units = "H",
      stringsAsFactors = FALSE
    ),
    summarize_map(field, weights)
  )
}

#' Compare D-MMM and GC on the same finite spatial area
#' @param metric_family,dimension,variable_set,metric_type Archived identity labels.
#' @param tau LDR threshold or NULL for ordinary Hellinger.
#' @param mmm_field,gc_field Numeric error/score fields of identical shape.
#' @param weights Nonnegative weights of that shape.
#' @param local_error_already_absolute Retained compatibility argument. Fraction
#'   improved always compares absolute values, which also applies to nonnegative H.
#' @return One-row data frame with paired means, GC-minus-D-MMM delta, relative
#'   improvement (NA for a nonpositive/missing baseline) and improved-area fraction.
#'   Both means use the intersection of finite fields with positive finite weights.
make_comparison_row <- function(metric_family,
                                dimension,
                                variable_set,
                                metric_type,
                                tau,
                                mmm_field,
                                gc_field,
                                weights,
                                local_error_already_absolute = FALSE) {
  weighted_values(mmm_field, weights)
  weighted_values(gc_field, weights)
  paired_weights <- weights
  paired_weights[!is.finite(mmm_field) | !is.finite(gc_field)] <- 0
  mmm_mean <- weighted_mean_value(mmm_field, paired_weights)
  gc_mean <- weighted_mean_value(gc_field, paired_weights)
  fraction_improved <- weighted_fraction_improved(gc_field, mmm_field, paired_weights)
  data.frame(
    metric_family = metric_family,
    dimension = dimension,
    variable_set = variable_set,
    metric_type = metric_type,
    tau = if (is.null(tau)) NA_real_ else tau,
    mmm_weighted_mean = mmm_mean,
    gc_weighted_mean = gc_mean,
    delta_vs_MMM = gc_mean - mmm_mean,
    improvement_pct_vs_MMM = if (is.finite(mmm_mean) && mmm_mean > 0 && is.finite(gc_mean))
      100 * (mmm_mean - gc_mean) / mmm_mean else NA_real_,
    fraction_area_improved_vs_MMM = fraction_improved,
    stringsAsFactors = FALSE
  )
}

#' Rank finite candidate and ensemble scores with smaller values preferred
#' @param scores Numeric vector with distinct nonempty candidate model names.
#' @param mmm_value,gc_value Scalar D-MMM and GC scores; nonfinite values are unranked.
#' @return One-row data frame of best/worst model identities and scores, unweighted
#'   type-7 ensemble quartiles, and minimum ranks among models plus both methods.
#'   No finite candidate scores yields missing summary fields.
summarize_ensemble_scores <- function(scores, mmm_value, gc_value) {
  if (!is.numeric(scores) || (length(scores) && (is.null(names(scores)) ||
      anyNA(names(scores)) || any(!nzchar(names(scores))) || anyDuplicated(names(scores)))) ||
      any(names(scores) %in% c("MMM", "GC-0.125")) || !is.numeric(mmm_value) ||
      length(mmm_value) != 1L || !is.numeric(gc_value) || length(gc_value) != 1L) {
    stop("Expected named distinct candidate scores and scalar numeric method scores.")
  }
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0L) {
    return(data.frame(
      best_model = NA_character_,
      best_model_value = NA_real_,
      ensemble_q25 = NA_real_,
      ensemble_median = NA_real_,
      ensemble_q75 = NA_real_,
      worst_model = NA_character_,
      worst_model_value = NA_real_,
      MMM_rank_among_models_and_methods = NA_integer_,
      GC_rank_among_models_and_methods = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  score_names <- names(scores)
  ord <- order(scores)
  combined <- c(scores, MMM = mmm_value, `GC-0.125` = gc_value)
  combined <- combined[is.finite(combined)]
  combined_rank <- rank(combined, ties.method = "min")
  data.frame(
    best_model = score_names[ord[1]],
    best_model_value = scores[ord[1]],
    ensemble_q25 = as.numeric(stats::quantile(scores, 0.25, na.rm = TRUE, names = FALSE)),
    ensemble_median = as.numeric(stats::quantile(scores, 0.50, na.rm = TRUE, names = FALSE)),
    ensemble_q75 = as.numeric(stats::quantile(scores, 0.75, na.rm = TRUE, names = FALSE)),
    worst_model = score_names[ord[length(ord)]],
    worst_model_value = scores[ord[length(ord)]],
    MMM_rank_among_models_and_methods = unname(combined_rank["MMM"]),
    GC_rank_among_models_and_methods = unname(combined_rank["GC-0.125"]),
    stringsAsFactors = FALSE
  )
}

distance_rows <- data.frame()
comparison_rows <- data.frame()
ensemble_rows <- data.frame()

#' Build ordinary and LDR distance summaries for one variable set
#' @param dimension,variable_set Archived dimension/variable labels.
#' @param pdf_ref,pdf_mmm,pdf_gc Normalised longitude x latitude x bins arrays.
#' @param pdf_models Normalised longitude x latitude x bins x candidates array.
#' @param full_mmm_map,full_gc_map Optional precomputed ordinary Hellinger maps;
#'   if either is NULL both are recomputed from PMFs.
#' @return List of distance, comparison and ensemble summary data frames for
#'   ordinary Hellinger plus each tau_values threshold. Uses caller-environment
#'   area_weight_matrix, model_names and tau_values; no output files are written.
#'   Shapes/model count and complete finite normalised PMFs are validated first.
add_distance_metric <- function(dimension,
                                variable_set,
                                pdf_ref,
                                pdf_models,
                                pdf_mmm,
                                pdf_gc,
                                full_mmm_map = NULL,
                                full_gc_map = NULL) {
  expected_grid <- dim(area_weight_matrix)
  ref_dims <- dim(pdf_ref)
  if (!is.numeric(pdf_ref) || length(ref_dims) != 3L ||
      !identical(ref_dims[1:2], expected_grid) || !identical(dim(pdf_mmm), ref_dims) ||
      !identical(dim(pdf_gc), ref_dims) ||
      !identical(dim(pdf_models), c(ref_dims, length(model_names)))) {
    stop("PMF grids, support bins and candidate count must match ERA5 evaluation inputs.")
  }
  # Bin-wise omission of missing model masses changes the distributional baseline.
  # Validate complete PMFs before either full or restricted distance evaluation.
  validate_target_pmfs(pdf_models)
  for (pmf in list(pdf_ref, pdf_mmm, pdf_gc)) {
    dim(pmf) <- c(ref_dims, 1L)
    validate_target_pmfs(pmf)
  }
  all_indices <- build_all_indices(pdf_ref)
  if (is.null(full_mmm_map) || is.null(full_gc_map)) {
    full_maps <- compute_h_maps_from_pdfs(pdf_ref, pdf_mmm, pdf_gc, all_indices)
    full_mmm_map <- full_maps$MMM
    full_gc_map <- full_maps$`GC-0.125`
  }
  full_mmm_map <- coerce_h_map(full_mmm_map)
  full_gc_map <- coerce_h_map(full_gc_map)
  if (!identical(dim(full_mmm_map), expected_grid) || !identical(dim(full_gc_map), expected_grid) ||
      any(!is.finite(full_mmm_map)) || any(full_mmm_map < 0) ||
      any(!is.finite(full_gc_map)) || any(full_gc_map < 0)) {
    stop("Precomputed ordinary Hellinger maps must be complete nonnegative maps on the evaluation grid.")
  }

  distance_rows_local <- list(
    make_distance_row("MMM", dimension, variable_set, "full_H", NULL, full_mmm_map, area_weight_matrix),
    make_distance_row("GC-0.125", dimension, variable_set, "full_H", NULL, full_gc_map, area_weight_matrix)
  )
  comparison_rows_local <- list(
    make_comparison_row(
      metric_family = "hellinger",
      dimension = dimension,
      variable_set = variable_set,
      metric_type = "full_H",
      tau = NULL,
      mmm_field = full_mmm_map,
      gc_field = full_gc_map,
      weights = area_weight_matrix
    )
  )

  model_maps <- compute_model_h_maps_from_pdfs(pdf_ref, pdf_models, all_indices)
  model_scores <- vapply(seq_len(dim(model_maps)[3]), function(m) {
    weighted_mean_value(matrix(model_maps[, , m], nrow = expected_grid[1], ncol = expected_grid[2]), area_weight_matrix)
  }, numeric(1))
  names(model_scores) <- model_names[seq_along(model_scores)]
  ensemble_rows_local <- list(cbind(
    data.frame(
      metric_family = "hellinger",
      dimension = dimension,
      variable_set = variable_set,
      metric_type = "full_H",
      tau = NA_real_,
      stringsAsFactors = FALSE
    ),
    summarize_ensemble_scores(
      scores = model_scores,
      mmm_value = weighted_mean_value(full_mmm_map, area_weight_matrix),
      gc_value = weighted_mean_value(full_gc_map, area_weight_matrix)
    )
  ))

  for (tau in tau_values) {
    ldr_indices <- build_ldr_indices(pdf_ref, tau)
    partial_maps <- compute_h_maps_from_pdfs(pdf_ref, pdf_mmm, pdf_gc, ldr_indices)
    distance_rows_local <- c(
      distance_rows_local,
      list(
        make_distance_row("MMM", dimension, variable_set, "partial_LDR_H", tau, partial_maps$MMM, area_weight_matrix),
        make_distance_row("GC-0.125", dimension, variable_set, "partial_LDR_H", tau, partial_maps$`GC-0.125`, area_weight_matrix)
      )
    )
    comparison_rows_local <- c(
      comparison_rows_local,
      list(
        make_comparison_row(
          metric_family = "hellinger",
          dimension = dimension,
          variable_set = variable_set,
          metric_type = "partial_LDR_H",
          tau = tau,
          mmm_field = partial_maps$MMM,
          gc_field = partial_maps$`GC-0.125`,
          weights = area_weight_matrix
        )
      )
    )

    model_partial_maps <- compute_model_h_maps_from_pdfs(pdf_ref, pdf_models, ldr_indices)
    model_partial_scores <- vapply(seq_len(dim(model_partial_maps)[3]), function(m) {
      weighted_mean_value(matrix(model_partial_maps[, , m], nrow = expected_grid[1], ncol = expected_grid[2]), area_weight_matrix)
    }, numeric(1))
    names(model_partial_scores) <- model_names[seq_along(model_partial_scores)]
    ensemble_rows_local <- c(
      ensemble_rows_local,
      list(cbind(
        data.frame(
          metric_family = "hellinger",
          dimension = dimension,
          variable_set = variable_set,
          metric_type = "partial_LDR_H",
          tau = tau,
          stringsAsFactors = FALSE
        ),
        summarize_ensemble_scores(
          scores = model_partial_scores,
          mmm_value = weighted_mean_value(partial_maps$MMM, area_weight_matrix),
          gc_value = weighted_mean_value(partial_maps$`GC-0.125`, area_weight_matrix)
        )
      ))
    )
  }

  list(
    distance = do.call(rbind, distance_rows_local),
    comparison = do.call(rbind, comparison_rows_local),
    ensemble = do.call(rbind, ensemble_rows_local)
  )
}

cat("Computing 3D Hellinger summaries...\n")
metric_result <- add_distance_metric(
  dimension = "3D",
  variable_set = "pr_tas_psl",
  pdf_ref = pdf3_ref_fut,
  pdf_models = pdf3_models_fut,
  pdf_mmm = pdf3_MMM_fut,
  pdf_gc = pdf3_GC_fut,
  full_mmm_map = MMM_hdist_fut,
  full_gc_map = GC_hdist_fut
)
distance_rows <- append_rows(distance_rows, list(metric_result$distance))
comparison_rows <- append_rows(comparison_rows, list(metric_result$comparison))
ensemble_rows <- append_rows(ensemble_rows, list(metric_result$ensemble))

gc_labels <- GC_result$label_attribution

cat("Computing 2D Hellinger summaries...\n")
for (pair_name in names(pdf2_future)) {
  pdf_full <- pdf2_future[[pair_name]]
  pdf_ref <- array(pdf_full[, , , 1L], dim = dim(pdf_full)[1:3])
  pdf_models <- pdf_full[, , , -1L, drop = FALSE]
  pdf_mmm <- compute_mmm_pdf_from_full_array(pdf_full)
  pdf_gc <- extract_gc_pdf_from_full_array(pdf_full, gc_labels)

  metric_result <- add_distance_metric(
    dimension = "2D",
    variable_set = pair_name,
    pdf_ref = pdf_ref,
    pdf_models = pdf_models,
    pdf_mmm = pdf_mmm,
    pdf_gc = pdf_gc
  )
  distance_rows <- append_rows(distance_rows, list(metric_result$distance))
  comparison_rows <- append_rows(comparison_rows, list(metric_result$comparison))
  ensemble_rows <- append_rows(ensemble_rows, list(metric_result$ensemble))
}

cat("Computing 1D Hellinger summaries...\n")
for (var_name in names(pdf1_future)) {
  pdf_full <- pdf1_future[[var_name]]
  pdf_ref <- array(pdf_full[, , , 1L], dim = dim(pdf_full)[1:3])
  pdf_models <- pdf_full[, , , -1L, drop = FALSE]
  pdf_mmm <- compute_mmm_pdf_from_full_array(pdf_full)
  pdf_gc <- extract_gc_pdf_from_full_array(pdf_full, gc_labels)

  metric_result <- add_distance_metric(
    dimension = "1D",
    variable_set = var_name,
    pdf_ref = pdf_ref,
    pdf_models = pdf_models,
    pdf_mmm = pdf_mmm,
    pdf_gc = pdf_gc
  )
  distance_rows <- append_rows(distance_rows, list(metric_result$distance))
  comparison_rows <- append_rows(comparison_rows, list(metric_result$comparison))
  ensemble_rows <- append_rows(ensemble_rows, list(metric_result$ensemble))
}

cat("Computing mean-field bias summaries...\n")
bias_rows <- data.frame()
bias_comparison_rows <- data.frame()
bias_ensemble_rows <- data.frame()

#' Extract selected candidate means from a reference-first mean array
#' @param mean_array Numeric longitude x latitude x models array, reference first.
#' @param gc_labels Finite integer longitude x latitude candidate labels in
#'   1..(number of full models minus one).
#' @return Longitude x latitude matrix selecting full model index label+1.
#'   Missing field values propagate; invalid labels or mismatched grids error.
extract_gc_mean_field <- function(mean_array, gc_labels) {
  dims <- dim(mean_array)
  if (!is.numeric(mean_array) || length(dims) != 3L || dims[3] < 2L) {
    stop("Expected mean array with dimensions lon x lat x model.")
  }
  n_models_excluding_ref <- dims[3] - 1L
  if (!is.numeric(gc_labels) || !identical(dim(gc_labels), dims[1:2]) ||
      any(!is.finite(gc_labels)) || any(gc_labels != floor(gc_labels)) ||
      any(gc_labels < 1L | gc_labels > n_models_excluding_ref)) {
    stop("GraphCut labels must be finite candidate indices on the mean-field grid.")
  }
  out <- matrix(NA_real_, nrow = dims[1], ncol = dims[2])
  for (label in seq_len(n_models_excluding_ref)) {
    mask <- gc_labels == label
    if (!any(mask, na.rm = TRUE)) {
      next
    }
    full_index <- label + 1L
    source <- matrix(mean_array[, , full_index], nrow = dims[1], ncol = dims[2])
    out[mask] <- source[mask]
  }
  out
}

#' Summarise mean-field error against the reference
#' @param method Internal method identifier; MMM is the distributional baseline.
#' @param variable Variable identifier pr, tas or psl.
#' @param method_field,ref_field Numeric spatial means on an identical grid.
#' @param weights Nonnegative weights of that grid shape.
#' @return One-row data frame of signed/absolute bias, RMSE, absolute-error
#'   quantiles, spatial correlation and valid-area fraction. Historical units are
#'   K for tas and explicitly unqualified native units for pr/psl; no conversion.
make_bias_row <- function(method, variable, method_field, ref_field, weights) {
  weighted_values(method_field, weights)
  weighted_values(ref_field, weights)
  bias <- method_field - ref_field
  abs_bias <- abs(bias)
  qs_abs <- weighted_quantile_value(abs_bias, weights, probs = c(0.50, 0.90, 0.95))
  data.frame(
    method = method,
    variable = variable,
    units = ifelse(variable == "tas", "K", ifelse(variable == "pr", "native_pr_units", "native_psl_units")),
    weighted_mean_bias = weighted_mean_value(bias, weights),
    weighted_mean_absolute_bias = weighted_mean_value(abs_bias, weights),
    weighted_rmse = weighted_rmse_against_zero(bias, weights),
    weighted_median_absolute_bias = qs_abs[1],
    weighted_p90_absolute_bias = qs_abs[2],
    weighted_p95_absolute_bias = qs_abs[3],
    weighted_spatial_correlation = weighted_spatial_correlation(method_field, ref_field, weights),
    valid_area_fraction = valid_area_fraction(bias, weights),
    stringsAsFactors = FALSE
  )
}

for (var_name in names(mean_future)) {
  mean_array <- mean_future[[var_name]]
  ref_field <- matrix(mean_array[, , 1L], nrow = length(lon), ncol = length(lat))
  model_fields <- mean_array[, , -1L, drop = FALSE]
  mmm_field <- apply(model_fields, c(1, 2), mean)
  gc_field <- extract_gc_mean_field(mean_array, gc_labels)

  bias_rows <- append_rows(bias_rows, list(
    make_bias_row("MMM", var_name, mmm_field, ref_field, area_weight_matrix),
    make_bias_row("GC-0.125", var_name, gc_field, ref_field, area_weight_matrix)
  ))

  mmm_abs_bias <- abs(mmm_field - ref_field)
  gc_abs_bias <- abs(gc_field - ref_field)
  bias_comparison <- make_comparison_row(
    metric_family = "mean_field_bias",
    dimension = "map",
    variable_set = var_name,
    metric_type = "mean_absolute_bias",
    tau = NA_real_,
    mmm_field = mmm_abs_bias,
    gc_field = gc_abs_bias,
    weights = area_weight_matrix,
    local_error_already_absolute = TRUE
  )
  mmm_mae <- bias_comparison$mmm_weighted_mean
  gc_mae <- bias_comparison$gc_weighted_mean
  bias_comparison_rows <- append_rows(bias_comparison_rows, list(bias_comparison))

  model_mae <- vapply(seq_len(dim(model_fields)[3]), function(m) {
    model_field <- matrix(model_fields[, , m], nrow = length(lon), ncol = length(lat))
    weighted_mean_value(abs(model_field - ref_field), area_weight_matrix)
  }, numeric(1))
  names(model_mae) <- model_names[seq_along(model_mae)]
  bias_ensemble_rows <- append_rows(bias_ensemble_rows, list(cbind(
    data.frame(
      metric_family = "mean_field_bias",
      dimension = "map",
      variable_set = var_name,
      metric_type = "mean_absolute_bias",
      tau = NA_real_,
      stringsAsFactors = FALSE
    ),
    summarize_ensemble_scores(
      scores = model_mae,
      mmm_value = mmm_mae,
      gc_value = gc_mae
    )
  )))
}

comparison_rows <- append_rows(comparison_rows, list(bias_comparison_rows))
ensemble_rows <- append_rows(ensemble_rows, list(bias_ensemble_rows))

config <- data.frame(
  key = c(
    "source_workspace",
    "created_at",
    "methods",
    "tau_values",
    "weighting",
    "smooth_cost",
    "n_lon",
    "n_lat",
    "n_models_excluding_reference"
  ),
  value = c(
    normalizePath(workspace_file, winslash = "/", mustWork = FALSE),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    paste(method_names, collapse = ", "),
    paste(tau_values, collapse = ", "),
    "cos(latitude), finite-cell means; method comparisons use common finite cells",
    if (exists("smooth_cost", inherits = TRUE)) as.character(smooth_cost) else NA_character_,
    as.character(length(lon)),
    as.character(length(lat)),
    as.character(length(model_names))
  ),
  stringsAsFactors = FALSE
)

distance_rows$metric_family <- "hellinger"
distance_rows <- distance_rows[
  c(
    "metric_family", "method", "dimension", "variable_set", "metric_type", "tau", "units",
    "weighted_mean", "weighted_median", "weighted_sd", "weighted_iqr",
    "weighted_p90", "weighted_p95", "valid_area_fraction"
  )
]

long_distance <- data.frame(
  metric_family = distance_rows$metric_family,
  method = distance_rows$method,
  dimension = distance_rows$dimension,
  variable_set = distance_rows$variable_set,
  metric_type = distance_rows$metric_type,
  tau = distance_rows$tau,
  metric_name = "weighted_mean",
  metric_value = distance_rows$weighted_mean,
  units = distance_rows$units,
  stringsAsFactors = FALSE
)

long_bias <- do.call(rbind, lapply(
  c(
    "weighted_mean_bias", "weighted_mean_absolute_bias", "weighted_rmse",
    "weighted_median_absolute_bias", "weighted_p90_absolute_bias",
    "weighted_p95_absolute_bias", "weighted_spatial_correlation"
  ),
  function(metric_name) {
    data.frame(
      metric_family = "mean_field_bias",
      method = bias_rows$method,
      dimension = "map",
      variable_set = bias_rows$variable,
      metric_type = metric_name,
      tau = NA_real_,
      metric_name = metric_name,
      metric_value = bias_rows[[metric_name]],
      units = bias_rows$units,
      stringsAsFactors = FALSE
    )
  }
))

long_comparison <- data.frame(
  metric_family = comparison_rows$metric_family,
  method = "GC-0.125_vs_MMM",
  dimension = comparison_rows$dimension,
  variable_set = comparison_rows$variable_set,
  metric_type = comparison_rows$metric_type,
  tau = comparison_rows$tau,
  metric_name = "improvement_pct_vs_MMM",
  metric_value = comparison_rows$improvement_pct_vs_MMM,
  units = "%",
  stringsAsFactors = FALSE
)

long_all_metrics <- rbind(long_distance, long_bias, long_comparison)

if (nrow(distance_rows) == 0L || nrow(bias_rows) == 0L ||
    nrow(comparison_rows) == 0L || nrow(ensemble_rows) == 0L ||
    nrow(long_all_metrics) == 0L) {
  stop("One or more ERA5 summary tables are unexpectedly empty.")
}

csv_file <- file.path(output_dir, paste0(stamp, "_ERA5_global_metrics_summary_long.csv"))

write.csv(long_all_metrics, csv_file, row.names = FALSE)

cat("ERA5 evaluation CSV:", csv_file, "\n")

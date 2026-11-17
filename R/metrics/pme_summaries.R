#' Compute a finite-pair weighted mean for a spatial PME score
#' @param x,w Numeric value/weight vectors of identical length.
#' @return Weighted arithmetic mean over finite pairs with nonnegative weights;
#'   NA_real_ if no positive valid total weight remains. Negative/nonfinite weights
#'   are omitted, not converted to positive weights. Does not weight references.
#' @examples
#' stopifnot(weighted_mean(c(1,3),c(1,1)) == 2)
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

#' Construct Article 1 cosine-latitude weights on a regular global grid
#' @param lon Unique finite longitude centres in degrees; their number determines
#'   row count. Equal-width longitude cells are assumed, not checked here.
#' @param lat Unique finite latitude centres in degrees north, in [-90,90].
#' @return Lon x lat matrix with cos(latitude) repeated along longitude. Weights
#'   are unnormalised; this is not an irregular-grid cell-area calculation.
#' @examples
#' build_weight_matrix(0:2,c(-30,30))
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

#' Resolve one PME map without silently reinterpreting a different grid shape
#' @param map_obj Numeric matrix of the expected shape, or a dimensionless legacy
#'   vector already in longitude-fastest column-major order.
#' @param expected_dim Positive integer c(nlon,nlat) shape.
#' @param label Map identity used in error messages.
#' @return Lon x lat matrix. A matrix/array with a different shape fails even
#'   when its element count matches; only a raw vector can be reshaped.
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
    "]."
  )
}

#' Require the complete held-out PME reference set for both methods
#' @param results Compact PME result list, including meta$model_names and method maps.
#' @param period_id One named evaluation-period identifier.
#' @param lambda_key Character GC smoothness key, usually "0.125".
#' @param tau_key Character LDR threshold key, e.g. "0.01"; ignored for full Hellinger.
#' @param metric_type Exactly full or partial (the internal identifier for LDR).
#' @return Character vector in metadata model order. Missing, extra or duplicate
#'   references on either method raise an error instead of reducing the sample.
#'   Publication entry points separately require the configured set of 15 models.
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

#' Summarise paired spatial PME scores separately for each reference and period
#'
#' First average each method's distance field with the supplied spatial weights,
#' then compute 100 * (H_D-MMM - H_GC) / H_D-MMM for that reference. This is not
#' an average of local percentage ratios and does not pool references. Both
#' methods must have the same valid-cell mask; matching missing cells retain
#' the historical common-area omission rule. Positive improvement favours GC.
#' @param results Compact PME result list with periods, metadata and method maps.
#' @param metric_type Full ordinary Hellinger or partial internal LDR identifier.
#' @param metric_label Human-readable label copied into each output row.
#' @param lambda_key GC smoothness map key, normally "0.125".
#' @param tau_key LDR threshold map key, ignored for full Hellinger.
#' @param weight_matrix Finite nonnegative lon x lat weights with positive total.
#' @return Data frame with metric_label, period_id, period, reference, H_MMM,
#'   H_GC, delta_H and rel_improvement in percent. H_MMM denotes D-MMM internally.
#'   Requires complete reference sets, distinct periods, correct map shapes and
#'   finite paired means with a strictly positive D-MMM denominator.
build_metric_means_from_rds <- function(results, metric_type, metric_label, lambda_key, tau_key, weight_matrix) {
  periods <- results$periods
  if (!is.data.frame(periods) || !nrow(periods) ||
      !all(c("period_id","label") %in% names(periods)) ||
      anyNA(periods$period_id) || any(!nzchar(as.character(periods$period_id))) ||
      anyDuplicated(periods$period_id)) stop("Expected distinct nonempty PME period identities.")
  if (!is.matrix(weight_matrix) || !is.numeric(weight_matrix) ||
      any(!is.finite(weight_matrix)) || any(weight_matrix<0) || sum(weight_matrix)<=0)
    stop("Expected a finite nonnegative spatial weight matrix with positive total.")
  period_order <- as.character(periods$period_id)
  period_labels <- as.character(periods$label)
  names(period_labels) <- period_order
  expected_dim <- dim(weight_matrix)

  rows <- list()
  cursor <- 1L
  for (period_id in period_order) {
    refs <- get_reference_order(results, period_id, lambda_key, tau_key, metric_type)
    if (metric_type == "full") {
      mmm_maps <- results$MMM_hdist_future_list[[period_id]]
      gc_maps <- results$GC_hdist_future_list[[period_id]][[lambda_key]]
    } else {
      mmm_maps <- results$MMM_partial_hdist_future_list[[period_id]][[tau_key]]
      gc_maps <- results$GC_partial_hdist_future_list[[period_id]][[tau_key]][[lambda_key]]
    }

    for (ref_name in refs) {
      mmm_map <- coerce_to_map_matrix(mmm_maps[[ref_name]], expected_dim, paste0("MMM ", metric_type, " ", period_id, "/", ref_name))
      gc_map <- coerce_to_map_matrix(gc_maps[[ref_name]], expected_dim, paste0("GC ", metric_type, " ", period_id, "/", ref_name))
      # Never compare methods averaged over different geographic subsets.
      if (!identical(is.finite(mmm_map),is.finite(gc_map)))
        stop("D-MMM/GC valid-cell masks differ for ",period_id,"/",ref_name)
      H_MMM <- weighted_mean(as.vector(mmm_map), as.vector(weight_matrix))
      H_GC <- weighted_mean(as.vector(gc_map), as.vector(weight_matrix))
      rows[[cursor]] <- data.frame(
        metric_label = metric_label,
        period_id = period_id,
        period = period_labels[[period_id]],
        reference = ref_name,
        H_MMM = H_MMM,
        H_GC = H_GC,
        delta_H = H_MMM - H_GC,
        rel_improvement = relative_improvement(H_MMM, H_GC),
        stringsAsFactors = FALSE
      )
      cursor <- cursor + 1L
    }
  }

  do.call(rbind, rows)
}

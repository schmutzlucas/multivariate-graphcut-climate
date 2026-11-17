# Public vector API. Production array calculations retain compute_partial_hdist().
#' Validate a probability mass function without modifying it
#' @param p Numeric nonempty vector of finite nonnegative masses.
#' @param tolerance Allowed absolute deviation of total mass from one.
#' @return p invisibly; raises an error for an invalid PMF. No renormalisation occurs.
#' @examples
#' validate_pmf(c(0, 0.25, 0.75))
validate_pmf <- function(p, tolerance = 1e-10) {
  if(!is.numeric(tolerance)||length(tolerance)!=1L||!is.finite(tolerance)||tolerance<0)
    stop("tolerance must be one finite non-negative number.")
  if (!is.numeric(p) || !length(p) || any(!is.finite(p)) || any(p < 0) ||
      abs(sum(p) - 1) > tolerance) stop("Expected a finite, non-negative PMF summing to one.")
  invisible(p)
}

#' Compute the Hellinger distance between two PMFs
#' @param p,q Numeric normalised PMFs on the same ordered bins.
#' @return Scalar sqrt(sum((sqrt(p) - sqrt(q))^2) / 2), between zero and one.
#' @examples
#' hellinger_distance(c(1, 0), c(0, 1))
hellinger_distance <- function(p, q) {
  validate_pmf(p); validate_pmf(q)
  if (length(p) != length(q)) stop("PMFs must share the same bins.")
  sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)
}

#' Compute non-renormalised LDR Hellinger score
#' @param reference,candidate Numeric normalised PMFs on identical ordered bins.
#' @param ell_ldr Scalar threshold strictly between zero and one, usually 0.01.
#' @return Scalar restricted Hellinger score on reference-selected LDR bins.
#' @details The selected subvectors are not renormalised. This measures an
#'   absolute contribution to squared Hellinger distance, not distance between
#'   conditional distributions. Ties and zero-reference bins follow select_ldr_indices().
#' @examples
#' ldr_hellinger_score(c(0, 0.1, 0.9), c(0.2, 0.1, 0.7), ell_ldr = 0.01)
ldr_hellinger_score <- function(reference, candidate, ell_ldr = .01) {
  validate_pmf(reference); validate_pmf(candidate)
  if (length(reference) != length(candidate)) stop("PMFs must share the same bins.")
  if (!is.numeric(ell_ldr) || length(ell_ldr) != 1L || !is.finite(ell_ldr) ||
      ell_ldr <= 0 || ell_ldr >= 1) stop("ell_ldr must lie strictly between zero and one.")
  bins <- select_ldr_indices(reference, ell_ldr)
  compute_partial_hdist(reference, candidate, bins)
}

#' Average a global field with cosine-latitude weights
#' @param field Numeric matrix indexed longitude x latitude, on equal-width
#'   longitude cells. Nonfinite field values are excluded with their weights.
#' @param lat Latitude centres in degrees north, one per field column, in [-90, 90].
#' @return Finite weighted mean, or NA_real_ if no positive valid weight remains.
#' @details Weights are proportional to cos(latitude), as used in Article 1.
#'   This is not a general irregular-grid cell-area calculator.
#' @examples
#' area_weighted_mean(matrix(c(1, 1, 3, 3), 2), c(0, 60))
area_weighted_mean <- function(field, lat) {
  if (!is.matrix(field) || ncol(field) != length(lat) || any(!is.finite(lat)) ||
      any(abs(lat) > 90)) stop("Expected a lon x lat matrix and matching latitudes in [-90,90].")
  weights <- matrix(rep(cos(lat * pi / 180), each = nrow(field)), nrow(field))
  valid <- is.finite(field) & is.finite(weights) & weights >= 0
  if (!any(valid) || sum(weights[valid]) <= 0) return(NA_real_)
  sum(field[valid] * weights[valid]) / sum(weights[valid])
}

#' Compute paired relative improvement over D-MMM
#' @param baseline Finite positive D-MMM scores, one per comparison.
#' @param candidate Finite GC scores in the same order and of the same length.
#' @return Numeric percentages 100 * (baseline - candidate) / baseline.
#' @details Positive values favour GC. Reference-wise paper summaries apply
#'   this after obtaining each held-out reference's area-weighted score; they
#'   do not average grid-cell percentage ratios or weight references by model size.
#' @examples
#' relative_improvement(c(4, 2), c(3, 1.5))
relative_improvement <- function(baseline, candidate) {
  if (length(baseline) != length(candidate) || any(!is.finite(baseline)) ||
      any(!is.finite(candidate)) || any(baseline <= 0)) {
    stop("Relative improvement requires finite paired scores and a positive baseline.")
  }
  100 * (baseline - candidate) / baseline
}

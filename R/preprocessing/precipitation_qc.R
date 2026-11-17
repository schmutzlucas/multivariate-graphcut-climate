#' Enforce nonnegative finite precipitation without modifying missing values
#' @param values Numeric precipitation vector/array in already prepared units.
#' @return Copy with every finite negative observation replaced by zero. Existing
#'   zeros/positives, dimensions, NA and NaN are preserved. Infinities are retained
#'   as invalid observations, not silently repaired; production validation rejects
#'   them. No unit conversion, file editing or logarithm is performed here.
#' @examples
#' stopifnot(identical(precipitation_qc(c(-1, 0, 2, NA_real_)), c(0, 0, 2, NA_real_)))
precipitation_qc <- function(values) {
  if (!is.numeric(values) || is.complex(values)) stop("Precipitation must be real numeric data.")
  negative <- is.finite(values) & values < 0
  if (any(negative)) values[negative] <- 0
  values
}

#' Apply the Article 1 nonnegative precipitation logarithm
#' @param values Numeric precipitation vector/array in prepared mm/day units.
#' @return log1p(precipitation_qc(values)), retaining dimensions and missing values.
#'   Infinite inputs stop explicitly. The shared transform is idempotent with
#'   respect to QC, not with respect to the logarithm: call once on physical data.
#' @examples
#' stopifnot(identical(precipitation_transform(c(-1, 0, 1)), c(0, 0, log1p(1))))
precipitation_transform <- function(values) {
  if (any(is.infinite(values))) stop("Infinite precipitation is invalid; inspect the input/QC report.")
  log1p(precipitation_qc(values))
}

#' Compute original-scale precipitation means before the logarithm
#' @param cube Numeric longitude x latitude x time array in prepared mm/day units.
#' @return Longitude x latitude means of nonnegative corrected observations,
#'   omitting missing values. All-missing cells return NA; infinities stop. No
#'   log/exp round trip is used, so means use corrected physical observations.
precipitation_mean_grid <- function(cube) {
  dimensions <- dim(cube)
  if (length(dimensions) != 3L || any(dimensions < 1L)) stop("Expected a nonempty lon x lat x time precipitation cube.")
  if (any(is.infinite(cube))) stop("Infinite precipitation in mean calculation.")
  cube <- precipitation_qc(cube)
  dim(cube) <- c(prod(dimensions[1:2]), dimensions[3])
  means <- rowMeans(cube, na.rm = TRUE)
  means[is.nan(means)] <- NA_real_
  matrix(means, dimensions[1], dimensions[2])
}

#' Collect additive raw-precipitation quality-control statistics
#' @param values Numeric physical precipitation vector/array before QC.
#' @return Named list with counts, extrema, raw/corrected sums and total absolute
#'   correction. Finite negatives alone are changed. Missing and infinite values
#'   are separately counted; an empty finite sample has missing extrema. These
#'   sufficient statistics can be merged across chunks without pooling means.
#'   Counters use doubles so global 150-year grids do not overflow R's 32-bit
#'   integers; integer counts remain exact below 2^53 observations.
precipitation_qc_statistics <- function(values) {
  corrected <- precipitation_qc(values)
  finite <- is.finite(values)
  negative <- finite & values < 0
  list(observations = as.double(length(values)), finite_observations = as.double(sum(finite)),
    missing_observations = as.double(sum(is.na(values))), unexpected_invalid_values = as.double(sum(is.infinite(values))),
    negative_observations = as.double(sum(negative)), values_changed = as.double(sum(negative)),
    raw_minimum = if (any(finite)) min(values[finite]) else NA_real_,
    corrected_minimum = if (any(finite)) min(corrected[finite]) else NA_real_,
    raw_sum = sum(values[finite]), corrected_sum = sum(corrected[finite]),
    correction_sum = -sum(values[negative]),
    maximum_absolute_correction = if (any(negative)) -min(values[negative]) else 0)
}

#' Merge raw-precipitation statistics over disjoint observation chunks
#' @param accumulated NULL for the initial chunk, or a prior statistics list.
#' @param chunk Statistics list from precipitation_qc_statistics(). Chunks must
#'   represent disjoint samples within one model/variable/scope.
#' @return Merged additive statistics with finite extrema and maximum correction.
#'   Derived fractions/means must be computed from merged sums/counts, not by
#'   averaging chunk means. Input observations are not retained.
merge_precipitation_qc_statistics <- function(accumulated, chunk) {
  if (is.null(accumulated)) return(chunk)
  extrema <- c("raw_minimum", "corrected_minimum", "maximum_absolute_correction")
  for (name in setdiff(names(chunk), extrema)) accumulated[[name]] <- as.double(accumulated[[name]]) + as.double(chunk[[name]])
  for (name in c("raw_minimum", "corrected_minimum")) {
    valid <- c(accumulated[[name]], chunk[[name]])
    accumulated[[name]] <- if (all(is.na(valid))) NA_real_ else min(valid, na.rm = TRUE)
  }
  accumulated$maximum_absolute_correction <- max(accumulated$maximum_absolute_correction, chunk$maximum_absolute_correction)
  accumulated
}

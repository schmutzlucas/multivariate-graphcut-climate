#' Select a tie-inclusive high-density region
#'
#' Sort reference-bin probabilities from largest to smallest and retain every
#' bin at least as probable as the first cutoff reaching mass 1 - tau.
#' @param pdf_vector Numeric vector of reference-bin masses. Positive total
#'   mass is normalised internally; nonfinite and negative entries are treated
#'   as zero for legacy compatibility. Use validate_pmf() for strict validation.
#' @param tau Scalar in [0, 1); nominal mass outside the HDR.
#' @return Sorted one-based integer bin indices, or integer(0) for zero total mass.
#' @details Boundary ties are all included; the attained mass may exceed 1 - tau.
#' @examples
#' select_hdr_indices(c(0, 0.1, 0.9), tau = 0.1)
select_hdr_indices <- function(pdf_vector, tau = 0.10) {
  if (!is.numeric(pdf_vector) || length(pdf_vector) == 0L) {
    return(integer(0))
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0 || tau >= 1) {
    stop("tau must be a single finite numeric value in [0, 1).")
  }

  p <- as.numeric(pdf_vector)
  p[!is.finite(p)] <- 0
  p[p < 0] <- 0

  total <- sum(p)
  if (!is.finite(total) || total <= 0) {
    return(integer(0))
  }

  # Work on normalized probabilities to make target mass explicit.
  p <- p / total
  target_mass <- 1 - tau

  sorted_indices_desc <- order(p, decreasing = TRUE)
  cumulative_sum_desc <- cumsum(p[sorted_indices_desc])
  threshold_index <- which(cumulative_sum_desc >= target_mass)[1]

  if (is.na(threshold_index)) {
    return(sort(sorted_indices_desc))
  }

  kappa <- p[sorted_indices_desc[threshold_index]]
  sort(which(p >= kappa))
}

#' Select the manuscript's tie-inclusive low-density region
#'
#' Select the smallest-probability reference bins up to the first cumulative
#' cutoff reaching tau, then include every bin tied at that cutoff.
#' @param pdf_vector Numeric reference masses. Positive total mass is normalised
#'   internally. Legacy handling maps negative/nonfinite entries to zero; the
#'   strict public score API validates PMFs before calling this selector.
#' @param tau Scalar in [0, 1), called ell_LDR in the manuscript. Tau is nominal,
#'   not a guarantee that selected mass equals the requested value.
#' @return Sorted one-based bin indices; integer(0) if total reference mass is zero.
#' @details Zero-reference bins are retained. Ties may make the selected mass
#'   larger than tau. Selection uses the reference only, never candidate masses.
#'   At tau = 0 the legacy selector still retains zero bins, or all minimum-mass
#'   ties if no zero exists. ldr_hellinger_score() requires a strictly positive threshold.
#' @examples
#' select_ldr_indices(c(0, 0.005, 0.005, 0.005, 0.985), tau = 0.01)
select_ldr_indices <- function(pdf_vector, tau = 0.10) {
  if (!is.numeric(pdf_vector) || length(pdf_vector) == 0L) {
    return(integer(0))
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0 || tau >= 1) {
    stop("tau must be a single finite numeric value in [0, 1).")
  }

  p <- as.numeric(pdf_vector)
  p[!is.finite(p)] <- 0
  p[p < 0] <- 0

  total <- sum(p)
  if (!is.finite(total) || total <= 0) {
    return(integer(0))
  }

  # LDR is the lowest-probability set whose cumulative reference mass reaches tau.
  p <- p / total
  sorted_indices_asc <- order(p, decreasing = FALSE)
  cumulative_sum_asc <- cumsum(p[sorted_indices_asc])
  threshold_index <- which(cumulative_sum_asc >= tau)[1]

  if (is.na(threshold_index)) {
    return(sort(sorted_indices_asc))
  }

  kappa <- p[sorted_indices_asc[threshold_index]]
  sort(which(p <= kappa))
}

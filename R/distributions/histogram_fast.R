# Fast, numerically equivalent equal-width histogram helpers.
#
# Scientific invariants:
# - supports are the existing local physical/transformed min/max values;
# - bins are axis aligned and equal width;
# - precipitation is transformed upstream as log(1 + pr), as in production;
# - flattened states use the same column-major convention as compute_histND().

#' Assign equal-width bins with explicit floating-edge handling
#' @param x Numeric observations; NA/NaN produce missing indices. Infinite and
#'   out-of-range observations saturate into the outermost bins.
#' @param xmin,xmax Finite scalar support limits with xmax strictly larger.
#' @param nbins Positive integer number of bins along this one variable.
#' @param zero_based Logical; return 0 through nbins - 1 rather than R's 1-based bins.
#' @return Integer bin indices of length length(x), preserving missing entries.
#' @details Values within eight scaled machine epsilons of an internal edge are
#'   snapped to that edge before flooring. This stabilises nested-bin derivation
#'   but is not bit-identical to findInterval for every value infinitesimally below
#'   an edge. Ordinary exact edges belong to the bin on their right; xmax belongs
#'   to the final bin. No physical support bounds are re-estimated here.
#' @examples
#' equal_width_bin_index(c(0, 0.25, 0.5, 1), 0, 1, 4)
equal_width_bin_index <- function(x, xmin, xmax, nbins, zero_based = TRUE) {
  if (!is.numeric(x) || !is.numeric(c(xmin, xmax)) ||
      length(xmin) != 1L || length(xmax) != 1L ||
      !is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) {
    stop("x must be numeric and xmin/xmax must be finite scalars with xmax > xmin.")
  }
  if (!is.numeric(nbins) || length(nbins) != 1L || is.na(nbins) ||
      !is.finite(nbins) ||
      nbins < 1L || nbins != as.integer(nbins)) {
    stop("nbins must be a single positive integer.")
  }
  nbins <- as.integer(nbins)

  scaled <- nbins * (x - xmin) / (xmax - xmin)

  # Values constructed from an internal edge should enter the bin on its right,
  # exactly as findInterval(..., rightmost.closed = TRUE). Floating arithmetic
  # can place such a value infinitesimally below an integer, so snap only values
  # within a few machine epsilons of an integer before flooring.
  nearest <- round(scaled)
  tol <- 8 * .Machine$double.eps * pmax(1, abs(scaled))
  snap <- is.finite(scaled) & abs(scaled - nearest) <= tol
  scaled[snap] <- nearest[snap]

  bin0 <- floor(scaled)
  bin0[is.infinite(x) & x < 0] <- 0
  bin0[is.infinite(x) & x > 0] <- nbins - 1L
  bin0[!is.na(bin0) & bin0 < 0] <- 0
  bin0[!is.na(bin0) & bin0 >= nbins] <- nbins - 1L
  bin0 <- as.integer(bin0)

  if (isTRUE(zero_based)) bin0 else bin0 + 1L
}

#' Coarsen zero-based indices on exactly nested equal-width grids
#' @param finest_bin0 Zero-based indices on the finest grid, possibly missing.
#' @param finest_nbins Positive integer finest-grid resolution.
#' @param target_nbins Positive integer divisor of finest_nbins.
#' @return Zero-based integer indices obtained by integer division. No samples
#'   are rebinned independently; all grids must have the same physical bounds.
#' @examples
#' derive_nested_bin_index(0:7, 8, 4)
derive_nested_bin_index <- function(finest_bin0, finest_nbins, target_nbins) {
  finest_nbins <- as.integer(finest_nbins)
  target_nbins <- as.integer(target_nbins)
  if (length(finest_nbins) != 1L || length(target_nbins) != 1L ||
      finest_nbins < 1L || target_nbins < 1L ||
      finest_nbins %% target_nbins != 0L) {
    stop("finest_nbins must be an integer multiple of target_nbins.")
  }
  as.integer(finest_bin0 %/% (finest_nbins %/% target_nbins))
}

#' Flatten and count multivariate bin states in R array order
#' @param bin_indices0 Observations x variables matrix of valid zero-based bin
#'   indices; incomplete rows are excluded. Callers supply indices below nbins.
#' @param nbins Number of bins per variable, identical along every axis.
#' @return Integer counts of length nbins raised to the number of variables.
#' @details The first variable varies fastest: for pr/tas/psl the state is
#'   b_pr + nbins*b_tas + nbins^2*b_psl. Returned counts are not normalised PMFs.
#' @examples
#' tabulate_linear_states(rbind(c(0, 1), c(1, 0)), 2)
tabulate_linear_states <- function(bin_indices0, nbins) {
  if(!is.numeric(nbins)||length(nbins)!=1L||!is.finite(nbins)||nbins<1L||nbins!=as.integer(nbins))
    stop("nbins must be a positive integer.")
  if (!is.matrix(bin_indices0)) {
    bin_indices0 <- as.matrix(bin_indices0)
  }
  n_vars <- ncol(bin_indices0)
  if (n_vars < 1L) {
    stop("bin_indices0 must contain at least one variable column.")
  }
  valid <- stats::complete.cases(bin_indices0)
  if (!any(valid)) {
    return(tabulate(integer(0), nbins = as.integer(nbins)^n_vars))
  }

  bins <- bin_indices0[valid, , drop = FALSE]
  if(!is.numeric(bins)||any(!is.finite(bins))||any(bins<0|bins>=nbins)||any(bins!=trunc(bins)))
    stop("Bin indices must be finite integers between zero and nbins - 1.")
  multipliers <- as.double(nbins) ^ (seq_len(n_vars) - 1L)
  # Variable order is pr, tas, psl. Thus pr varies fastest:
  # state0 = b_pr + k*b_tas + k^2*b_psl.
  state1 <- as.vector(bins %*% multipliers) + 1L
  tabulate(state1, nbins = as.integer(nbins)^n_vars)
}

#' Count an equal-width histogram using vectorised indices
#' @param data Numeric observations x variables matrix or matrix-coercible data.
#' @param range_var Variables x 2 matrix of finite minimum/maximum bounds.
#' @param nbins Positive integer bins per axis, not total multivariate states.
#' @return Integer count vector of length nbins^ncol(data), in column-major bin order.
#' @details Applies equal_width_bin_index() separately to each variable and
#'   excludes rows with missing bin assignments. One-observation and empty inputs
#'   retain their correct multivariate dimensionality. See floating-edge semantics
#'   in equal_width_bin_index(); this function does not normalise counts.
#' @examples
#' compute_histND_equal_width_fast(matrix(c(0.2, 0.8), 1), cbind(c(0, 0), c(1, 1)), 8)
compute_histND_equal_width_fast <- function(data, range_var, nbins) {
  if (!is.matrix(data)) data <- as.matrix(data)
  if (!is.numeric(data)) stop("data must be numeric.")
  if (!is.matrix(range_var) || ncol(range_var) != 2L ||
      nrow(range_var) != ncol(data) || !is.numeric(range_var) ||
      any(!is.finite(range_var)) || any(range_var[, 2] <= range_var[, 1])) {
    stop("range_var must be finite [ncol(data), 2] with max > min.")
  }

  bins0 <- vapply(
    seq_len(ncol(data)),
    function(v) equal_width_bin_index(
      data[, v], range_var[v, 1], range_var[v, 2], nbins,
      zero_based = TRUE
    ),
    integer(nrow(data))
  )
  # vapply drops the observation axis when exactly one observation is present.
  bins0 <- matrix(bins0, nrow = nrow(data), ncol = ncol(data))
  tabulate_linear_states(bins0, nbins)
}

#' Build the Article 1 univariate, bivariate and trivariate counts together
#' @param pr,tas,psl Equally sized numeric observation vectors. Precipitation
#'   must already have the production log(1 + precipitation) transform applied.
#' @param range_3d A 3 x 2 matrix of transformed/physical bounds ordered pr,tas,psl.
#' @param nbins3d,nbins2d,nbins1d Per-axis resolutions for three-, two-, and
#'   one-dimensional histograms. Article 1 fixes these to 8, 16 and 32.
#' @param profile Logical; append indexing/tabulation timings in seconds.
#' @return List pdf3 (joint counts), pdf2 (pr_tas, pr_psl, tas_psl count vectors),
#'   and pdf1 (pr, tas, psl counts), with optional timing. The pdf names are legacy:
#'   values are counts until the calling PMF builder normalises each histogram.
#' @details Nested resolutions reuse finest-grid indices. Non-nested settings
#'   compute their grids separately. Missing values are excluded per diagnostic,
#'   so marginal sample counts can differ when raw variables contain missing data.
#' @examples
#' x <- c(0.1, 0.5, 0.9)
#' counts <- compute_multires_hist_counts_fast(x, x, x, cbind(rep(0, 3), rep(1, 3)))
#' sum(counts$pdf3)
compute_multires_hist_counts_fast <- function(
  pr,
  tas,
  psl,
  range_3d,
  nbins3d = 8L,
  nbins2d = 16L,
  nbins1d = 32L,
  profile = FALSE
) {
  if (!(length(pr) == length(tas) && length(tas) == length(psl))) {
    stop("pr, tas, and psl must have equal lengths.")
  }
  if (!is.matrix(range_3d) || !identical(dim(range_3d), c(3L, 2L))) {
    stop("range_3d must be a 3 x 2 matrix ordered pr, tas, psl.")
  }

  choices<-list(nbins3d,nbins2d,nbins1d)
  if(any(lengths(choices)!=1L)||!all(vapply(choices,is.numeric,TRUE)))
    stop("All bin counts must be single positive integers.")
  requested<-unlist(choices)
  if(any(!is.finite(requested))||any(requested<1)||any(requested!=trunc(requested)))
    stop("All bin counts must be positive integers.")
  resolutions <- as.integer(requested)
  if (any(resolutions < 1L)) stop("All bin counts must be positive integers.")
  finest <- max(resolutions)
  nested <- all(finest %% resolutions == 0L)
  values <- list(pr = pr, tas = tas, psl = psl)

  t_index <- proc.time()[["elapsed"]]
  bins <- lapply(seq_along(values), function(v) {
    fine <- equal_width_bin_index(
      values[[v]], range_3d[v, 1], range_3d[v, 2], finest,
      zero_based = TRUE
    )
    if (nested) {
      list(
        b3 = derive_nested_bin_index(fine, finest, nbins3d),
        b2 = derive_nested_bin_index(fine, finest, nbins2d),
        b1 = derive_nested_bin_index(fine, finest, nbins1d)
      )
    } else {
      list(
        b3 = equal_width_bin_index(values[[v]], range_3d[v, 1], range_3d[v, 2], nbins3d),
        b2 = equal_width_bin_index(values[[v]], range_3d[v, 1], range_3d[v, 2], nbins2d),
        b1 = equal_width_bin_index(values[[v]], range_3d[v, 1], range_3d[v, 2], nbins1d)
      )
    }
  })
  indexing_seconds <- proc.time()[["elapsed"]] - t_index

  t_tabulate <- proc.time()[["elapsed"]]
  counts <- list(
    pdf3 = tabulate_linear_states(cbind(bins[[1]]$b3, bins[[2]]$b3, bins[[3]]$b3), nbins3d),
    pdf2 = list(
      pr_tas = tabulate_linear_states(cbind(bins[[1]]$b2, bins[[2]]$b2), nbins2d),
      pr_psl = tabulate_linear_states(cbind(bins[[1]]$b2, bins[[3]]$b2), nbins2d),
      tas_psl = tabulate_linear_states(cbind(bins[[2]]$b2, bins[[3]]$b2), nbins2d)
    ),
    pdf1 = list(
      pr = tabulate_linear_states(matrix(bins[[1]]$b1, ncol = 1L), nbins1d),
      tas = tabulate_linear_states(matrix(bins[[2]]$b1, ncol = 1L), nbins1d),
      psl = tabulate_linear_states(matrix(bins[[3]]$b1, ncol = 1L), nbins1d)
    )
  )
  tabulation_seconds <- proc.time()[["elapsed"]] - t_tabulate

  if (isTRUE(profile)) {
    counts$timing <- c(
      bin_index = indexing_seconds,
      tabulation = tabulation_seconds
    )
  }
  counts
}

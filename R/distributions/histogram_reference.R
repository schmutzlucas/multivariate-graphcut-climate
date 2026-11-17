#' Compute N-dimensional Histogram
#'
#' This function computes an n-dimensional histogram from a given dataset, where each variable
#' is assigned to a specific bin based on the provided range and number of bins. The function
#' returns a flattened vector representing the histogram counts for each combination of bins
#' across all variables.
#'
#' @param data A matrix where each column represents a variable and each row is an observation.
#' @param range_var A matrix of dimension `[n_vars, 2]` that defines the minimum and maximum values
#' for each variable. The first column corresponds to the minimum values, and the second column
#' corresponds to the maximum values.
#' @param nbins A numeric value specifying the number of bins for each variable.
#'
#' @return A numeric vector representing the flattened n-dimensional histogram. The vector has length
#' equal to `nbins^n_vars`, where `n_vars` is the number of variables.
#' @details These are integer counts, not normalised PMFs. The first variable's
#'   bin varies fastest in the flattened array. Out-of-support values saturate
#'   into edge bins; rows containing missing observations do not contribute.
#'
#' @examples
#' # Example data with two variables
#' data <- matrix(seq(-1, 1, length.out = 100), ncol = 2) # 50 paired observations
#' range_var <- cbind(c(-3, -3), c(3, 3)) # Min/max columns
#' nbins <- 10  # Number of bins per variable
#'
#' # Compute the n-dimensional histogram
#' hist_vector <- compute_histND(data, range_var, nbins)
#'
#' @export
compute_histND <- function(data, range_var, nbins) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("data must be numeric.")
  }

  n_vars <- ncol(data)
  n_obs <- nrow(data)

  if (!is.numeric(nbins) || length(nbins) != 1L || is.na(nbins) || nbins < 1L || nbins != as.integer(nbins)) {
    stop("nbins must be a single positive integer.")
  }
  nbins <- as.integer(nbins)

  if (!is.matrix(range_var) || ncol(range_var) != 2L || nrow(range_var) != n_vars) {
    stop("range_var must be a numeric matrix with ncol=2 and nrow=ncol(data).")
  }
  if (!is.numeric(range_var) || any(!is.finite(range_var))) {
    stop("range_var must contain finite numeric values.")
  }
  if (any(range_var[, 2] <= range_var[, 1])) {
    stop("Each variable range must satisfy max > min in range_var.")
  }

  if (n_obs == 0L) {
    return(tabulate(integer(0), nbins = nbins^n_vars))
  }

  # Compute bin edges for each variable
  bin_edges <- lapply(seq_len(n_vars), function(v) {
    seq(range_var[v, 1], range_var[v, 2], length.out = nbins + 1L)
  })

  # For each variable, determine the bin index for each observation
  bin_indices <- matrix(sapply(seq_len(n_vars), function(v) {
    inds <- findInterval(data[, v], vec = bin_edges[[v]], rightmost.closed = TRUE)
    inds[inds < 1L] <- 1L
    inds[inds > nbins] <- nbins
    inds
  }), nrow = n_obs, ncol = n_vars)

  # Compute linear indices for the n-dimensional histogram.
  dims <- rep(nbins, n_vars)
  multiplier <- cumprod(c(1L, dims[-length(dims)]))
  idx_linear <- as.vector((as.matrix(bin_indices) - 1L) %*% multiplier) + 1L

  # Build histogram vector using tabulate; its length is exactly nbins^n_vars.
  hist_vector <- tabulate(idx_linear, nbins = prod(dims))

  hist_vector
}

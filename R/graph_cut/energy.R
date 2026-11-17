# Transparent double-precision reference for testing the production float callbacks.
#' Enumerate four-neighbour edges on a longitude-periodic grid
#' @param nlon Integer longitude count, at least three.
#' @param nlat Integer latitude count, at least one.
#' @return Two-column integer matrix p/q of one-based flattened pixel IDs, one
#'   row per undirected edge. Longitude varies fastest. Every row wraps east-west;
#'   latitude never wraps and no diagonal edge is added.
#' @examples
#' cylindrical_edges(3, 2)
cylindrical_edges <- function(nlon, nlat) {
  if (length(nlon) != 1L || length(nlat) != 1L || nlon < 3 || nlat < 1 ||
      nlon != as.integer(nlon) || nlat != as.integer(nlat)) {
    stop("Use at least three longitude cells and one latitude row.")
  }
  ids <- matrix(seq_len(nlon * nlat), nlon, nlat)
  edges <- cbind(as.vector(ids), as.vector(ids[c(2:nlon, 1L), , drop = FALSE]))
  if (nlat > 1L) edges <- rbind(edges,
    cbind(as.vector(ids[, -nlat, drop = FALSE]), as.vector(ids[, -1L, drop = FALSE])))
  colnames(edges) <- c("p", "q")
  edges
}

#' Evaluate the Article 1 graph-cut energy independently in double precision
#' @param labels One-based integer-valued longitude x latitude label matrix.
#' @param calibration_h Nonnegative longitude x latitude x models unary distances.
#' @param target_pmfs Normalised longitude x latitude x bins x models PMFs.
#' @param lat Latitude centres in degrees north, matching the second array axis.
#' @param lambda Nonnegative smoothness multiplier, normally 0.125.
#' @param weight_data Nonnegative multiplier of the unary term, normally one.
#' @return Named numeric vector data, smooth and total, including their multipliers.
#' @details The seam is the mean endpoint cosine weight times the sum of endpoint
#'   Hellinger distances between the edge labels' candidate PMFs. This R reference
#'   does not optimise labels and is not bit-identical to native float arithmetic.
#' @examples
#' graph_cut_energy(matrix(1, 3, 2), array(0.1, c(3, 2, 2)),
#'                  array(0.5, c(3, 2, 2, 2)), c(-30, 30))
graph_cut_energy <- function(labels, calibration_h, target_pmfs, lat, lambda = .125,
                             weight_data = 1) {
  dims <- dim(target_pmfs)
  if (length(dims) != 4L || !is.matrix(labels) ||
      !identical(as.integer(dim(labels)), as.integer(dims[1:2])) ||
      !identical(as.integer(dim(calibration_h)), as.integer(dims[c(1, 2, 4)])) ||
      length(lat) != dims[2] || any(!is.finite(lat)) || any(abs(lat) > 90) ||
      any(!is.finite(target_pmfs)) || any(target_pmfs < 0) ||
      any(!is.finite(calibration_h)) || any(calibration_h < 0) ||
      anyNA(labels) || any(labels != as.integer(labels)) || any(labels < 1 | labels > dims[4]) ||
      length(lambda) != 1L || !is.finite(lambda) || lambda < 0 ||
      length(weight_data) != 1L || !is.finite(weight_data) || weight_data < 0) {
    stop("Incompatible labels, PMFs, unary distances, latitude or energy weights.")
  }
  if (any(abs(apply(target_pmfs, c(1, 2, 4), sum) - 1) > 1e-10)) stop("Target PMFs must sum to one.")
  ncell <- prod(dims[1:2]); labs <- as.vector(labels)
  weights <- rep(abs(cos(lat * pi / 180)), each = dims[1])
  unary <- matrix(calibration_h, ncell, dims[4])
  pmfs <- array(target_pmfs, c(ncell, dims[3], dims[4]))
  data <- weight_data * sum(weights * unary[cbind(seq_len(ncell), labs)])
  edges <- cylindrical_edges(dims[1], dims[2])
  seam <- vapply(seq_len(nrow(edges)), function(e) {
    p <- edges[e, 1]; q <- edges[e, 2]; m <- labs[p]; n <- labs[q]
    hp <- sqrt(sum((sqrt(pmfs[p, , m]) - sqrt(pmfs[p, , n]))^2)) / sqrt(2)
    hq <- sqrt(sum((sqrt(pmfs[q, , m]) - sqrt(pmfs[q, , n]))^2)) / sqrt(2)
    (weights[p] + weights[q]) / 2 * (hp + hq)
  }, numeric(1))
  smooth <- lambda * sum(seam)
  c(data = data, smooth = smooth, total = data + smooth)
}

#' Measure unweighted label fragmentation with cylindrical neighbours
#' @param labels Complete longitude x latitude matrix of categorical model labels.
#' @return Percentage of undirected neighbouring cell pairs with different labels.
#' @details Includes periodic longitude edges and ordinary north-south edges,
#'   with no diagonal or polar wrap. This is an edge fraction, not area-weighted
#'   patch size or the graph-cut seam energy.
#' @examples
#' label_fragmentation(matrix(1, 3, 2))
label_fragmentation <- function(labels) {
  if (!is.matrix(labels) || anyNA(labels)) stop("Expected a complete label matrix.")
  edges <- cylindrical_edges(nrow(labels), ncol(labels))
  100 * mean(as.vector(labels)[edges[, 1]] != as.vector(labels)[edges[, 2]])
}

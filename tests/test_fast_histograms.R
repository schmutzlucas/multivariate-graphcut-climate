source("R/load.R")
load_mvgc()
cfg <- article1_config()

assert_identical <- function(actual, expected, label) {
  if (!identical(actual, expected)) {
    stop(label, " failed. max_abs_diff=", max(abs(actual - expected)))
  }
}

compare_one <- function(data, ranges, nbins, label) {
  valid <- stats::complete.cases(data)
  legacy <- compute_histND(data[valid, , drop = FALSE], ranges, nbins)
  fast <- compute_histND_equal_width_fast(data, ranges, nbins)
  assert_identical(fast, legacy, label)
}

set.seed(20260924)
edges <- seq(-2, 3, length.out = 9L)
boundary_values <- c(
  edges,
  edges[-length(edges)] - 1e-12,
  edges[-1L] + 1e-12,
  -Inf, Inf, NA_real_
)
compare_one(
  matrix(boundary_values, ncol = 1L),
  matrix(c(-2, 3), nrow = 1L),
  8L,
  "1D boundaries"
)

ranges3 <- rbind(
  pr = c(0, log1p(100)),
  tas = c(230, 330),
  psl = c(85000, 108000)
)
n <- 5000L
data3 <- cbind(
  pr = log1p(runif(n, -0.2, 120)),
  tas = runif(n, 220, 340),
  psl = runif(n, 84000, 109000)
)
data3[sample.int(n, 20), sample.int(3, 1)] <- NA_real_

for (k in c(4L, 7L, 8L, 12L, 16L, 32L)) {
  compare_one(data3, ranges3, k, paste0("random 3D k=", k))
  compare_one(data3[, 1:2, drop = FALSE], ranges3[1:2, , drop = FALSE], k, paste0("random 2D k=", k))
}

compare_multires <- function(nbins3d, nbins2d, nbins1d, label) {
  got <- compute_multires_hist_counts_fast(
    data3[, 1], data3[, 2], data3[, 3], ranges3,
    nbins3d, nbins2d, nbins1d
  )
  expected3 <- compute_histND(data3[complete.cases(data3), , drop = FALSE], ranges3, nbins3d)
  assert_identical(got$pdf3, expected3, paste(label, "3D"))

  pair_cols <- list(pr_tas = 1:2, pr_psl = c(1, 3), tas_psl = 2:3)
  for (name in names(pair_cols)) {
    cols <- pair_cols[[name]]
    valid <- complete.cases(data3[, cols, drop = FALSE])
    expected <- compute_histND(data3[valid, cols, drop = FALSE], ranges3[cols, , drop = FALSE], nbins2d)
    assert_identical(got$pdf2[[name]], expected, paste(label, name))
  }
  for (v in seq_len(3L)) {
    valid <- !is.na(data3[, v])
    expected <- compute_histND(matrix(data3[valid, v], ncol = 1L), ranges3[v, , drop = FALSE], nbins1d)
    assert_identical(got$pdf1[[v]], expected, paste(label, rownames(ranges3)[v]))
  }
}

compare_multires(8L, 16L, 32L, "nested")
compare_multires(7L, 11L, 13L, "generic fallback")

cat("All fast histogram equivalence tests passed.\n")

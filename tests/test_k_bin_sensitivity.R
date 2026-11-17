source("R/load.R")
load_mvgc()
cfg <- article1_config()
#!/usr/bin/env Rscript

load_pmf_sensitivity_cpp()

set.seed(5102)
n <- 240L
pr <- log(runif(n, 0, 20) + 1)
tas <- rnorm(n, 280, 7)
psl <- rnorm(n, 101000, 900)
ranges <- rbind(range(pr), range(tas), range(psl))
splits <- list(rep(c(TRUE, FALSE), length.out = n), rep(c(FALSE, TRUE), length.out = n))
k_values <- 4:12

actual <- cell_pmfs_for_k(pr, tas, psl, ranges, splits, k_values)
for (k in k_values) {
  for (s in 1:2) {
    expected_counts <- compute_histND(
      cbind(pr, tas, psl)[splits[[s]], , drop = FALSE], ranges, k
    )
    expected <- expected_counts / sum(expected_counts)
    stopifnot(identical(actual[[as.character(k)]][[s]], expected))
    stopifnot(identical(hellinger_distance_pmf(expected, expected), 0))
  }
}

raw_pr <- array(pr, dim = c(1L, 1L, n))
raw_tas <- array(tas, dim = c(1L, 1L, n))
raw_psl <- array(psl, dim = c(1L, 1L, n))
range_block <- array(as.numeric(ranges), dim = c(1L, 1L, 3L, 2L))
compiled <- multik_hist_counts_cpp(raw_pr, raw_tas, raw_psl, range_block, splits[[1]], k_values)
for (ki in seq_along(k_values)) {
  k <- k_values[ki]
  nstates <- k^3
  offset <- compiled$offsets[ki]
  for (s in 1:2) {
    from <- offset + (s - 1L) * nstates + 1L
    observed_counts <- compiled$counts[from:(from + nstates - 1L)]
    expected_counts <- compute_histND(
      cbind(pr, tas, psl)[splits[[s]], , drop = FALSE], ranges, k
    )
    stopifnot(identical(observed_counts, expected_counts))
  }
}
compiled_h <- multik_hellinger_cpp(compiled, compiled)
stopifnot(all(compiled_h == 0))

# Rank stability must be invariant to strictly monotone score transforms.
a <- c(0.1, 0.4, 0.2, 0.5, 0.3)
b <- a^3
stopifnot(isTRUE(all.equal(unname(cor(a, b, method = "spearman")), 1)))

cat("k-sensitivity unit tests passed.\n")

source("R/load.R"); load_mvgc()
set.seed(42)
for (i in 1:200) {
  p <- rexp(20); p <- p / sum(p); q <- rexp(20); q <- q / sum(q)
  h <- hellinger_distance(p, q)
  stopifnot(hellinger_distance(p, p) == 0, h >= 0, h <= 1,
            identical(h, hellinger_distance(q, p)),
            identical(h, compute_partial_hdist(p, q, seq_along(p))))
}
stopifnot(abs(hellinger_distance(c(1,0), c(0,1)) - 1) < 1e-15)
p <- c(0, .005, .005, .005, .005, .98); q <- c(.01,.002,.004,.008,.006,.97)
stopifnot(identical(select_ldr_indices(p, .01), 1:5),
          sum(p[select_ldr_indices(p, .01)]) == .02,
          abs(ldr_hellinger_score(p,q) - sqrt(sum((sqrt(p[1:5])-sqrt(q[1:5]))^2)/2)) < 1e-15)
stopifnot(abs(ldr_hellinger_score(c(0,.1,.9), c(.2,.1,.7))^2 - .1) < 1e-15)
for (i in 1:300) {
  p <- sample(0:20, 512, TRUE); p <- p / sum(p); q <- rexp(512); q <- q / sum(q)
  for (ell in c(.005,.01,.025,.05,.1)) {
    bins <- select_ldr_indices(p, ell)
    cutoff <- max(p[bins]); smaller <- which(p < cutoff)
    stopifnot(identical(bins, which(p <= cutoff)), all(which(p == 0) %in% bins),
              sum(p[bins]) >= ell - 1e-12, sum(p[smaller]) < ell + 1e-12,
              ldr_hellinger_score(p,q,ell) <= hellinger_distance(p,q) + 1e-12)
  }
}
stopifnot(identical(select_ldr_indices(rep(0, 4), .01), integer(0)),
          abs(area_weighted_mean(matrix(c(1,1,3,3),2),c(0,60))-5/3) < 1e-14,
          area_weighted_mean(matrix(c(NA,NA,3,3),2),c(0,60)) == 3,
          is.na(area_weighted_mean(matrix(NA_real_,2,2),c(0,60))),
          relative_improvement(4,3) == 25)
cat("PASS: Hellinger properties, 1,500 tied/zero LDR cases, non-renormalised scores and area weights.\n")

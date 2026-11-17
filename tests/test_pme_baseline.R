source("R/load.R");load_mvgc();source("R/experiments/perfect_model.R")
set.seed(74)
for (shape in list(c(1L,1L),c(2L,1L),c(1L,3L),c(3L,2L))) {
  p <- array(runif(prod(shape)*8),c(shape,8L))
  p <- sweep(p,c(1,2),rowSums(p,dims=2),"/")
  q <- array(runif(prod(shape)*8*3),c(shape,8L,3L))
  q <- sweep(q,c(1,2,4),apply(q,c(1,2,4),sum),"/")
  expected_mixture <- apply(q,1:3,mean)
  mixture <- compute_mmm_pdf(q)
  stopifnot(identical(mixture,expected_mixture))
  expected <- matrix(NA_real_,shape[1],shape[2])
  bins <- rep(list(rep(list(1:8),shape[2])),shape[1])
  for(i in seq_len(shape[1]))for(j in seq_len(shape[2]))
    expected[i,j] <- sqrt(sum((sqrt(expected_mixture[i,j,])-sqrt(p[i,j,]))^2))/sqrt(2)
  full <- compute_mmm_hdist(p,q)
  partial <- compute_mmm_partial_hdist(p,q,bins)
  stopifnot(identical(full,expected),identical(partial,expected),
    identical(compute_mmm_hdist(p,q,mixture),full),
    identical(compute_mmm_partial_hdist(p,q,bins,mixture),partial))
  ldr <- build_ldr_indices(p,.01,shape[1],shape[2])
  expected_ldr <- matrix(compute_partial_hdist(p,mixture,ldr),shape[1],shape[2])
  stopifnot(identical(compute_mmm_partial_hdist(p,q,ldr,mixture),expected_ldr))
}
for (bad_value in c(NA_real_,NaN,Inf,-.1,0)) {
  p <- array(c(.2,.8),c(1,1,2));q <- array(.5,c(1,1,2,2))
  bad <- p;bad[1] <- bad_value
  stopifnot(inherits(try(compute_mmm_hdist(bad,q),silent=TRUE),"try-error"),
    inherits(try(compute_mmm_partial_hdist(bad,q,list(list(1:2))),silent=TRUE),"try-error"))
  bad_q <- q;bad_q[1] <- bad_value
  stopifnot(inherits(try(compute_mmm_pdf(bad_q),silent=TRUE),"try-error"))
}
p <- array(c(.2,.8),c(1,1,2));q <- array(.5,c(1,1,2,2))
stopifnot(inherits(try(compute_mmm_partial_hdist(p,q,list(list(integer()))),silent=TRUE),"try-error"),
  inherits(try(compute_mmm_hdist(p,q,array(.5,c(1,2,1))),silent=TRUE),"try-error"))
cat("PASS equal-model D-MMM mixture, cached/uncached full/LDR scores, singleton maps and invalid-PMF rejection\n")

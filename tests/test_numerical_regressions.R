source("R/load.R");load_mvgc()
# A single observation must produce one count, not d^2 counts or d univariate bins.
for(n in 0:2)for(d in 1:3) {
 x<-matrix(.5,n,d);bounds<-cbind(rep(0,d),rep(1,d))
 expected<-numeric(8^d)
 if(n>0)expected[1+sum(4*8^(0:(d-1))) ]<-n
 for(fun in list(compute_histND,compute_histND_equal_width_fast)) {
  y<-fun(x,bounds,8)
  stopifnot(length(y)==8^d,sum(y)==n,identical(as.numeric(y),expected))
 }
}
p<-array(c(.2,.8),c(1,1,2));q<-array(c(.3,.7),c(1,1,2,1))
indices<-list(list(1:2))
one<-compute_partial_hdist(p,q,indices)
three<-compute_partial_hdist(p,array(q,c(1,1,2)),indices)
stopifnot(identical(dim(one),c(1L,1L,1L)),identical(one,three),
 abs(one[1]-hellinger_distance(as.vector(p),as.vector(q)))<1e-15)
for(tolerance in list(-1,NA_real_,Inf,c(0,1)))
 stopifnot(inherits(try(validate_pmf(c(.5,.5),tolerance),silent=TRUE),"try-error"))
stopifnot(isTRUE(validate_target_pmfs(array(.5,c(3,2,2,2)))))
for(value in c(-.1,NA_real_,Inf,.9)) {
 bad<-array(.5,c(3,2,2,2));bad[1]<-value
 stopifnot(inherits(try(validate_target_pmfs(bad),silent=TRUE),"try-error"))
}
cat("PASS audited singleton histogram dimensions and singleton-model Hellinger arrays\n")

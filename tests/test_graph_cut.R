source("R/load.R"); load_mvgc()
edges <- cylindrical_edges(360L, 181L)
stopifnot(nrow(edges) == 129960L,
          !anyDuplicated(cbind(pmin(edges[,1],edges[,2]),pmax(edges[,1],edges[,2]))),
          sum(abs(edges[,1]-edges[,2]) == 359L) == 181L)
if (!requireNamespace("gcoWrapR", quietly = TRUE) || !requireNamespace("RcppXPtrUtils", quietly = TRUE)) {
  cat("SKIP: native solver unavailable; cylindrical edge construction passed.\n"); quit(status=0)
}
require_production_solver()
set.seed(11)
pmfs <- array(runif(4*2*8*2),c(4,2,8,2))
pmfs <- sweep(pmfs,c(1,2,4),apply(pmfs,c(1,2,4),sum),"/")
unary <- array(runif(16),c(4,2,2)); lat <- c(-30,30)
comparisons <- list()
for (lambda in c(0,.125,.25)) {
  fit <- GraphCutHellinger_nD_lat(pmfs,unary,1,lambda,8L,lat,1L,FALSE,FALSE)
  e <- graph_cut_energy(fit$label_attribution,unary,pmfs,lat,lambda)
  native <- unlist(fit[["Data and smooth cost"]])
  comparisons[[length(comparisons)+1L]] <- data.frame(lambda=lambda,
    R_data=e[["data"]],native_data=unname(native[1]),R_smooth=e[["smooth"]],
    native_smooth=unname(native[2]),max_abs_difference=max(abs(unname(native)-e[c("data","smooth")])))
  stopifnot(max(abs(unname(native)-e[c("data","smooth")])) < 2e-6)
  candidates <- as.matrix(expand.grid(rep(list(1:2),8)))
  energies <- apply(candidates,1,function(z) graph_cut_energy(matrix(z,4,2),unary,pmfs,lat,lambda)[["total"]])
  stopifnot(abs(e[["total"]]-min(energies)) < 2e-6)
  again <- GraphCutHellinger_nD_lat(pmfs,unary,1,lambda,8L,lat,1L,FALSE,FALSE)
  stopifnot(identical(fit$label_attribution,again$label_attribution))
}
write.csv(do.call(rbind,comparisons),file.path(mvgc_output_dir("verification"),"native_energy.csv"),row.names=FALSE)
cat("PASS: periodic edges, float callbacks versus independent R energy, exhaustive binary optimum and seed repeatability.\n")

source("R/load.R"); load_mvgc()
x <- readRDS("tests/fixtures/era5_production_patch.rds")
# This cropped historical patch validates distribution components, not the
# optimum of the full accepted experiment. No native solver is needed here.
stopifnot(identical(dim(x$calibration_models),c(4L,2L,512L,22L)))
for (key in c("calibration_reference","target_reference","target_d_mmm","target_gc"))
 stopifnot(max(abs(apply(x[[key]],c(1,2),sum)-1)) < 1e-10)
for (model in seq_len(22L)) for (i in seq_len(4L)) for (j in seq_len(2L)) {
 p <- x$calibration_reference[i,j,]
 q <- x$calibration_models[i,j,,model]
 independent <- sqrt(sum((sqrt(p)-sqrt(q))^2)/2)
 stopifnot(abs(hellinger_distance(p,q)-independent)<1e-14,
           abs(independent-x$calibration_h[i,j,model])<1e-12)
}
for (i in seq_len(4L)) for (j in seq_len(2L)) {
 independent <- function(q) sqrt(sum((sqrt(x$target_reference[i,j,])-sqrt(q))^2)/2)
 stopifnot(abs(independent(x$target_d_mmm[i,j,])-x$full_d_mmm[i,j])<1e-12,
           abs(independent(x$target_gc[i,j,])-x$full_gc[i,j])<1e-12)
}
cat("PASS: 176 real calibration distances and 16 real target distances from independently evaluated PMFs.\n")

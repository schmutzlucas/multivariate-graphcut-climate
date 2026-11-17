source("R/load.R");load_mvgc();source("R/preprocessing/support_builder.R")
x<-array(c(-2,-1e-12,0,2,NA_real_,NaN),c(2,1,3))
original<-x;expected<-x;expected[1:2]<-0
stopifnot(identical(precipitation_qc(x),expected),identical(x,original),
 identical(precipitation_transform(x),log1p(expected)),
 identical(precipitation_qc(precipitation_qc(x)),expected),
 identical(precipitation_mean_grid(x),matrix(c(0,1),2,1)))
q<-precipitation_qc_statistics(x)
stopifnot(q$negative_observations==2,q$values_changed==2,q$missing_observations==2,
 q$raw_minimum==-2,q$corrected_minimum==0,q$maximum_absolute_correction==2,
 q$unexpected_invalid_values==0,abs(q$correction_sum-(2+1e-12))<1e-15)
split<-merge_precipitation_qc_statistics(precipitation_qc_statistics(x[1:3]),precipitation_qc_statistics(x[4:6]))
stopifnot(all(vapply(names(q),function(name)abs(q[[name]]-split[[name]])<1e-14,logical(1))))
# Production PME pools exceed 2^31-1 observations. Test the boundary without
# allocating billions of values, including legacy integer-valued chunk counts.
large<-precipitation_qc_statistics(c(0,1))
large$observations<-large$finite_observations<-2147483600L
chunk<-precipitation_qc_statistics(c(0,1))
chunk$observations<-chunk$finite_observations<-100L
combined<-merge_precipitation_qc_statistics(large,chunk)
stopifnot(is.double(combined$observations),combined$observations==2147483700,
 combined$finite_observations==2147483700,
 all(vapply(q[c("observations","finite_observations","missing_observations",
  "unexpected_invalid_values","negative_observations","values_changed")],is.double,logical(1))))
global<-precipitation_qc_statistics(c(0,1));global$observations<-global$finite_observations<-3571005600
stopifnot(merge_precipitation_qc_statistics(global,global)$observations==7142011200)
stopifnot(precipitation_qc_statistics(c(Inf,-Inf,NA_real_))$unexpected_invalid_values==2,
 inherits(try(precipitation_transform(c(Inf,-Inf)),silent=TRUE),"try-error"))
positive<-array(c(2,3,4,5),c(1,1,4));p<-compute_cube_min_max(positive,"pr")
stopifnot(p$min==log1p(2),p$max==log1p(5),p$min>0)
for(first in c(0,-1e-6)) {
 zero<-positive;zero[1]<-first;r<-compute_cube_min_max(zero,"pr")
 stopifnot(r$min==0,r$max==log1p(5))
}
for(value in c(0,log1p(3))) {
 r<-regularize_empirical_support(matrix(value,1,1),matrix(value,1,1),"pr")
 stopifnot(r$min==value,r$max>value,all(r$constant),r$expansion>0,
 sum(compute_histND(matrix(rep(value,5),ncol=1),matrix(c(r$min,r$max),nrow=1),32))==5)
}
missing<-compute_cube_min_max(array(NA_real_,c(1,1,3)),"pr")
stopifnot(is.na(missing$min),is.na(missing$max),
 is.na(regularize_empirical_support(missing$min,missing$max,"pr")$min))
# The real disk readers must agree despite permuted axes/repeated time requests.
path<-tempfile(fileext=".nc")
dims<-list(ncdf4::ncdim_def("time","days since 2000-01-01",0:2),
 ncdf4::ncdim_def("lat","degrees_north",0:1),ncdf4::ncdim_def("lon","degrees_east",0:1))
n<-ncdf4::nc_create(path,ncdf4::ncvar_def("pr","mm day-1",dims,-999,prec="double"))
raw<-array(c(-1e-6,NA_real_,0,2,3,4,-1e-12,5,6,7,8,9),c(3,2,2))
# ncdf4 may replace missing values in its supplied buffer; isolate that write
# buffer so the independent expected observations retain their NA entries.
ncdf4::ncvar_put(n,"pr",raw+0);ncdf4::nc_close(n)
n<-ncdf4::nc_open(path);info<-get_variable_dimension_info(n,"pr","time")
block<-read_variable_time_blocks(n,"pr",info,1:2,1:2,1:3)
slab<-read_variable_time_chunk(n,"pr",info,1:2,1:2,1,3,report_qc=TRUE)
stopifnot(identical(block,precipitation_qc(aperm(raw,c(3,2,1)))),identical(block,slab$data),
 slab$qc$values_changed==2,identical(ncdf4::ncvar_get(n,"pr",collapse_degen=FALSE),raw))
ncdf4::nc_close(n)
cat("PASS precipitation QC, missing/invalid values, raw means, empirical positive/zero/constant/missing supports and both disk readers\n")

# Exercise the real multi-model dispatcher in separate worker processes. The
# existing fixture includes finite negatives, zeros, positives and missing pr,
# unequal variable file start years and exact legacy/full/blocked comparisons.
source("tests/test_variable_date_alignment.R")
require_packages(c("future","future.apply"))
first_model<-model;model<-"FIXTURE_2"
write_variable("pr",1950);write_variable("tas",1949);write_variable("psl",1950)
models<-c(first_model,model)
for(backend in c("legacy","optimized")) {
 got<-compute_nd_pdf_multi_seasonal(variables=variables,model_names=models,data_dir=root,
  year_present=1951L,year_future=1952L,lon=0:1,lat=0:1,range_var=bounds,
  season_start_md=c(1,1),season_end_md=c(12,31),nbins3d=8L,nbins2d=16L,nbins1d=32L,
  workers=2L,pmf_backend=backend,spatial_block_size=1L,diagnostics=FALSE,verbose=FALSE)
 for(m in seq_along(models))for(period in c("present","future")) {
  stopifnot(identical(got$pdf3[[period]][,,,m],a[[paste0("pdf3_",period)]]))
  for(pair in c("pr_tas","pr_psl","tas_psl"))
   stopifnot(identical(got$pdf2[[period]][[pair]][,,,m],a[[paste0("pdf2_",period)]][[pair]]))
  for(variable in variables) {
   stopifnot(identical(got$pdf1[[period]][[variable]][,,,m],a[[paste0("pdf1_",period)]][[variable]]),
    identical(got$mean[[period]][[variable]][,,m],a[[paste0("mean_",period)]][[variable]]))
  }
 }
}
cat("PASS two-process legacy/optimized multi-model dispatch: QC PMFs and physical means identical to serial reference\n")

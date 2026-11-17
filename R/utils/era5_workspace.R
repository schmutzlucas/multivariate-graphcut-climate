# Save explicit scientific inputs/results, not the caller's entire interactive state.
#' Save the explicit ERA5 scientific-workspace whitelist
#' @param file Destination RData path; normal R save() overwrite behaviour applies.
#' @param envir Environment containing experiment objects, normally the global
#'   Rscript environment. Private enclosing environments of saved functions can
#'   capture extra objects; the validation driver deliberately avoids that case.
#' @param stage Either all (the historical scientific whitelist) or pmfs (only
#'   pre-optimisation inputs, raw PMF arrays and means). The pmfs stage excludes
#'   stale score/label objects remaining in an interactive caller environment.
#' @return Character names actually saved, invisibly. Missing whitelist objects
#'   are omitted. Writes an uncompressed workspace, potentially many GB, including
#'   PMFs/support/configuration/maps and the solver function but not arbitrary
#'   interactive objects. Does not create the parent directory or validate a run.
save_era5_workspace <- function(file,envir=parent.frame(),stage=c("all","pmfs")) {
  stage <- match.arg(stage)
  objects <- c("lon","lat","lon_size","lat_size","variables","year_present","year_future",
    "reference_name","model_names","ref_index","nbins3d","nbins2d","nbins1d","nbins_total3d",
    "smooth_cost","run_metadata","range_var_gc","pdf3_present","pdf3_future",
    "pdf3_ref_pres","pdf3_ref_fut","pdf3_models_pres","pdf3_models_fut",
    "pdf3_GC_pres","pdf3_GC_fut","pdf3_MMM_pres","pdf3_MMM_fut",
    "pdf1_present","pdf1_future","pdf2_present","pdf2_future","mean_present","mean_future",
    "h_dist_pres","h_dist_fut","GC_result","GC_hdist_pres","GC_hdist_fut",
    "MMM_hdist_pres","MMM_hdist_fut","GraphCutHellinger_nD_lat")
  if(stage=="pmfs")objects<-intersect(objects,c("lon","lat","lon_size","lat_size","variables",
    "year_present","year_future","reference_name","model_names","ref_index","nbins3d","nbins2d","nbins1d",
    "run_metadata","range_var_gc","pdf3_present","pdf3_future","pdf1_present","pdf1_future",
    "pdf2_present","pdf2_future","mean_present","mean_future"))
  available <- intersect(objects,ls(envir=envir,all.names=TRUE))
  save(list=available,envir=envir,file=file,compress=FALSE)
  invisible(available)
}

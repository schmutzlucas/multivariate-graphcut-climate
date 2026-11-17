# Source from the repository root. No formal R package or package attachment is needed.
#' Load the reusable Article 1 functions without running an experiment
#' @param envir Environment in which to define configuration and reusable functions;
#'   defaults to the caller. Existing definitions of the same names are replaced.
#' @return Character vector of sourced module paths, invisibly. Paths are relative
#'   to the repository root. Does not download data or source experiment entry points.
#' @examples
#' e <- new.env()
#' load_mvgc(e)
#' stopifnot(exists("hellinger_distance", envir=e, inherits=FALSE))
load_mvgc <- function(envir = parent.frame()) {
  modules <- c(
    "config/article1.R", "R/utils/runtime.R", "R/utils/input_registry.R", "R/utils/help.R", "R/utils/era5_workspace.R",
    "R/distributions/histogram_reference.R", "R/distributions/histogram_fast.R",
    "R/preprocessing/calendar.R", "R/preprocessing/cf_calendar_validation.R", "R/preprocessing/precipitation_qc.R", "R/preprocessing/netcdf_io.R", "R/preprocessing/era5_utc.R",
    "R/distributions/pmf_pipeline.R", "R/distributions/bin_sensitivity.R",
    "R/metrics/ldr.R", "R/metrics/partial_hellinger.R", "R/metrics/scores.R",
    "R/graph_cut/energy.R", "R/graph_cut/solver.R")
  for (path in modules) sys.source(path, envir = envir)
  invisible(modules)
}

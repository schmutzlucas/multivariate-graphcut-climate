# Current CRAN setup is opt-in, not a restoration of the tested environment.
args <- commandArgs(TRUE)
if (!"--accept-current-cran" %in% args)
  stop("Use DEPENDENCIES.md for restoration, or explicitly pass --accept-current-cran.")
repos <- "https://cloud.r-project.org"
packages <- c("ncdf4", "abind", "PCICt", "ncdf4.helpers", "future", "future.apply",
  "matrixStats", "Rcpp", "RcppXPtrUtils", "remotes", "digest", "ggplot2", "dplyr",
  "rlang", "patchwork", "maps", "png", "jsonlite", "openssl", "reshape2",
  "RColorBrewer", "scales", "ragg", "gtable")
if ("--with-docs" %in% args) packages <- c(packages, "roxygen2")
# Packages removed from current CRAN are installed from the CRAN Archive at the
# tested version recorded in config/software_versions.csv.
archived <- c(ncdf4.helpers = "0.3-7", PCICt = "0.5-4.4")
#' Check availability of a selected dependency without attaching it
#' @param p Package name.
#' @return Logical namespace availability; package startup is suppressed.
available <- function(p) requireNamespace(p, quietly = TRUE)
missing <- packages[!vapply(packages, available, logical(1))]
if (length(missing)) install.packages(missing, repos = repos)
for (p in intersect(names(archived), packages)) {
  if (available(p)) next
  url <- sprintf("%s/src/contrib/Archive/%s/%s_%s.tar.gz", repos, p, p, archived[[p]])
  message("Installing ", p, " ", archived[[p]], " from the CRAN Archive: ", url)
  install.packages(url, repos = NULL, type = "source")
}
missing <- packages[!vapply(packages, available, logical(1))]
if (length(missing)) stop("Dependency installation incomplete. Missing required packages: ",
                          paste(missing, collapse = ", "))
cat("PASS: selected dependencies available in this environment.\n")
cat("Current CRAN setup does not restore historical versions. Native solver remains separate.\n")

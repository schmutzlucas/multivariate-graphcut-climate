# Exercise failure verification without network or package installation.
env <- new.env(parent=globalenv())
env$commandArgs <- function(...) c("--accept-current-cran")
env$requireNamespace <- function(...) FALSE
env$install.packages <- function(...) invisible(NULL)
failure <- tryCatch({sys.source("tools/install_dependencies.R",env);NULL},error=identity)
stopifnot(inherits(failure,"error"),
          grepl("Missing required packages",conditionMessage(failure),fixed=TRUE),
          grepl("ragg",conditionMessage(failure),fixed=TRUE),
          !grepl("roxygen2",conditionMessage(failure),fixed=TRUE))
cat("PASS: incomplete mandatory installation fails; docs optional by default.\n")

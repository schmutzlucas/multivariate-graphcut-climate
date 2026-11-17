# Runtime checks never install packages or silently select a different dataset.
#' Require installed R dependencies without installing or attaching them
#' @param packages Character vector of package names.
#' @return TRUE invisibly if all namespaces load; otherwise an error listing
#'   missing packages and the documented installation entry point.
#' @examples
#' require_packages(c("stats", "utils"))
require_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing R packages: ", paste(absent, collapse = ", "),
                           ". See DEPENDENCIES.md or tools/install_dependencies.R.")
  invisible(TRUE)
}

#' Create and resolve a configured local output directory
#' @param component Relative output subdirectory, or an empty string for the root.
#' @return Path under article1_config()$output_root. Creates missing directories
#'   recursively; callers should pass repository-controlled relative components.
mvgc_output_dir <- function(component = "") {
  path <- file.path(article1_config()$output_root, component)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' Check whether the configured production data collection is discoverable
#' @param config Configuration list with input_registry and data_root fields.
#' @return Registry directory or data root invisibly. A registry is structurally
#'   checked; a directory-only root is checked for existence. This is not a payload
#'   checksum, coverage, unit, or complete-input validation. No computation starts.
require_climate_data <- function(config = article1_config()) {
  if(nzchar(config$input_registry)) {
    read_input_registry(config$input_registry)
    return(invisible(dirname(config$input_registry)))
  }
  if (!nzchar(config$data_root) || !dir.exists(config$data_root)) {
    stop("Climate data unavailable. Set MVGC_DATA_ROOT to the prepared daily ",
         "ERA5/CMIP6 NetCDF hierarchy. Level 1 and the synthetic demonstration ",
         "do not need these data. No production computation was started.")
  }
  invisible(config$data_root)
}

#' Require the pinned custom longitude-periodic graph-cut solver
#' @param config Configuration containing solver_repository and solver_commit.
#' @return TRUE invisibly if gcoWrapR/RcppXPtrUtils load and gcoWrapR RemoteSha
#'   matches the pin. Missing/unknown/wrong versions raise an error; nothing is installed.
require_production_solver <- function(config = article1_config()) {
  require_packages(c("gcoWrapR", "RcppXPtrUtils"))
  installed <- packageDescription("gcoWrapR")$RemoteSha
  if (!identical(installed, config$solver_commit)) {
    stop("Article 1 requires ", config$solver_repository, "@", config$solver_commit,
         "; installed RemoteSha is ", if (is.null(installed)) "unknown" else installed,
         ". The custom longitude-periodic solver is required.")
  }
  invisible(TRUE)
}

#' Translate archived support-path aliases without loading their contents
#' @param path Character vector of recorded support paths.
#' @return Character vector replacing the two historical ranges/ aliases with
#'   bundled data/supports paths; other entries are unchanged. This does not
#'   select corrected supports or establish equivalence of their numerical values.
#' @examples
#' canonical_support_path("data/supports/pme_1950-2100.rds")
canonical_support_path <- function(path) {
  aliases <- c(
    "ranges/range_var_final_allModelsPar_1950-2100_90deg_3v_PME_fullyear.rds" = "data/supports/pme_1950-2100.rds",
    "ranges/range_var_final_ERA5_1950-2023_3v.rds" = "data/supports/era5_1950-2023.rds")
  ifelse(path %in% names(aliases), unname(aliases[path]), path)
}

#' Require an explicitly selected external production result file
#' @param variable Environment-variable name holding the path.
#' @param description Human-readable description included in a missing-file error.
#' @return Existing path as a character string. Does not scan directories, choose
#'   a newest file, load a workspace, or verify its scientific/input identity.
production_result_file <- function(variable, description) {
  path <- Sys.getenv(variable, "")
  if (!nzchar(path) || !file.exists(path)) {
    stop("Set ", variable, " to ", description, ". Compact bundled outputs can ",
         "be rendered without loading a production workspace.")
  }
  path
}

#' Fingerprint scientific source without invalidating R caches for comments alone
#' @param paths Nonempty distinct existing source paths, normally repository-relative.
#' @return Named SHA-256 character vector. R files are parsed without source
#'   references and canonically deparsed, so ordinary comments/formatting are
#'   ignored but executable/default changes are retained. Other files (including
#'   native C++) use complete byte hashes. No source statements are evaluated.
#' @examples
#' scientific_code_fingerprint("R/metrics/ldr.R")
scientific_code_fingerprint <- function(paths) {
  if(!is.character(paths)||!length(paths)||anyNA(paths)||any(!nzchar(paths))||
     anyDuplicated(paths)||any(!file.exists(paths)))stop("Expected distinct existing scientific source paths.")
  require_packages("openssl")
  hashes<-vapply(paths,function(path) {
    if(tolower(tools::file_ext(path))!="r")return(mvgc_sha256(path))
    code<-paste(deparse(parse(path,keep.source=FALSE),width.cutoff=500L),collapse="\n")
    paste(as.character(openssl::sha256(charToRaw(code))),collapse="")
  },"")
  hashes
}

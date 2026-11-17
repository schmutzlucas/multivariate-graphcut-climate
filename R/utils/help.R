#' Open a generated function reference without installing an R package
#' @param topic Function name, as a single character string.
#' @param source Optional repository-relative source file to disambiguate helpers
#'   with the same name in multiple modules.
#' @param open Logical; open the HTML page in a browser. Defaults to interactive().
#' @return Absolute HTML path invisibly. Errors if no unique generated page exists.
#' @details Run Rscript tools/build_reference_docs.R after editing roxygen blocks.
#'   Plain sourcing does not register standard R help databases, so ?topic is not
#'   guaranteed to find a repository-local function. This helper is explicit and
#'   leaves the scientific project as a script-based repository.
#' @examples
#' if (file.exists("docs/reference/topics.csv")) {
#'   mvgc_help("hellinger_distance", open = FALSE)
#' }
mvgc_help <- function(topic, source=NULL, open=interactive()) {
 if(length(topic)!=1L||!is.character(topic)||is.na(topic))stop("Provide one function name")
 index<-"docs/reference/topics.csv"
 if(!file.exists(index))stop("Build the function reference first: Rscript tools/build_reference_docs.R")
 entries<-read.csv(index,stringsAsFactors=FALSE)
 rows<-entries[entries$function_name==topic&entries$documented,,drop=FALSE]
 if(!is.null(source))rows<-rows[rows$source==source,,drop=FALSE]
 if(nrow(rows)!=1L)stop("No unique documented topic; supply source= for duplicate names")
 path<-normalizePath(rows$html,winslash="/",mustWork=TRUE)
 if(isTRUE(open))utils::browseURL(path)
 invisible(path)
}

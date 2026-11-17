# Zero-copy, exact-file access. Never fall back to the old directory on failure.
#' Compute a payload SHA-256 through a binary file connection
#' @param path Path to a readable existing file, including institutional SMB paths.
#' @return Lowercase hexadecimal SHA-256 string. Requires openssl, reads the
#'   complete file, closes its connection on exit, and never changes the payload.
mvgc_sha256 <- function(path) {
  require_packages("openssl")
  con<-file(path,"rb");on.exit(close(con))
  paste(as.character(openssl::sha256(con)),collapse="")
}

#' Read and validate the frozen input registry without reading climate payloads
#' @param path Registry JSON path, defaulting to MVGC_INPUT_REGISTRY.
#' @param required_models Optional distinct model identifiers whose complete
#'   variable coverage must be present before starting an experiment.
#' @param required_variables Distinct expected variables for each required model,
#'   default pr/tas/psl; only used when required_models is supplied.
#' @return Parsed list with roots/inputs and any producer metadata unchanged.
#'   Checks nonempty schema, unique identities, safe relative paths, byte sizes,
#'   optional checksum syntax and requested model-variable completeness. Does not
#'   open NetCDFs, establish complete dates/units or freshly hash payloads.
read_input_registry <- function(path=Sys.getenv("MVGC_INPUT_REGISTRY", ""),
                                required_models=NULL,
                                required_variables=c("pr","tas","psl")) {
  require_packages("jsonlite")
  if(!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)||!file.exists(path))
    stop("Frozen input registry unavailable: ",paste(path,collapse=", "))
  r<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  #' Recognise one nonempty registry string
  #' Internal closure of read_input_registry; not a standalone entry point.
  #' @param x Potential scalar character field from an input registry.
  #' @return TRUE only for one nonmissing nonempty character value.
  scalar_string<-function(x)is.character(x)&&length(x)==1L&&!is.na(x)&&nzchar(x)
  if(!is.list(r)||!is.list(r$roots)||!length(r$roots)||is.null(names(r$roots))||
     anyNA(names(r$roots))||any(!nzchar(names(r$roots)))||anyDuplicated(names(r$roots))||
     !all(vapply(r$roots,scalar_string,logical(1)))||!is.list(r$inputs)||!length(r$inputs))
    stop("Input registry must contain nonempty named roots and input records.")
  for(i in r$inputs) {
    if(!is.list(i)||!all(vapply(i[c("model","variable","root_id","relative_path")],scalar_string,logical(1)))||
       !all(c("model","variable","root_id","relative_path","bytes")%in%names(i))||
       !i$root_id%in%names(r$roots))stop("Malformed input registry record or unknown root.")
    if(!is.numeric(i$bytes)||length(i$bytes)!=1L||!is.finite(i$bytes)||i$bytes<0||i$bytes!=floor(i$bytes))
      stop("Registry byte sizes must be finite nonnegative integers.")
    if(grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)",i$relative_path)||
       grepl("^[/\\\\]|^[A-Za-z]:",i$relative_path))stop("Unsafe registry path")
    if(!is.null(i$sha256)&&(!scalar_string(i$sha256)||!grepl("^[0-9a-fA-F]{64}$",i$sha256)))
      stop("Malformed registry SHA-256 for ",i$model,"/",i$variable)
  }
  keys<-vapply(r$inputs,function(i)paste(i$model,i$variable,sep="/"),"")
  if(anyDuplicated(keys))stop("Duplicate model/variable in input registry")
  if(!is.null(required_models)) {
    for(x in list(required_models,required_variables))
      if(!is.character(x)||!length(x)||anyNA(x)||any(!nzchar(x))||anyDuplicated(x))
        stop("Required model/variable identities must be distinct nonempty strings.")
    expected<-as.vector(outer(required_models,required_variables,paste,sep="/"))
    selected<-keys[vapply(r$inputs,function(i)i$model%in%required_models,logical(1))]
    if(!setequal(expected,selected))stop("Frozen registry does not contain exactly the required model-variable inputs.")
  }
  r
}

#' Resolve an exact model-variable payload from the active frozen registry
#' @param model,variable Single model/variable identifiers, matched exactly.
#' @param required_years Nonempty finite integer years; optional recorded endpoint
#'   years must bracket the request. This endpoint test does not prove no dates are missing.
#' @return Exact existing payload path after schema/path/size/endpoint checks.
#'   Never falls back to directory discovery. Full SHA-256 verification belongs
#'   to preparation and the explicit prepared-run launcher, not every point read.
registry_input_file <- function(model,variable,required_years) {
  if(!is.numeric(required_years)||!length(required_years)||any(!is.finite(required_years))||
     any(required_years!=floor(required_years)))stop("Required input years must be finite integers.")
  r<-read_input_registry()
  rows<-Filter(function(i)identical(i$model,model)&&identical(i$variable,variable),r$inputs)
  if(length(rows)!=1L)stop("No unique frozen input for ",model,"/",variable)
  row<-rows[[1]];root<-r$roots[[row$root_id]]
  if(is.null(root)||grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)",row$relative_path)||
     grepl("^[/\\\\]|^[A-Za-z]:",row$relative_path))stop("Unsafe registry path")
  path<-file.path(root,row$relative_path)
  if(!file.exists(path)||file.info(path)$size!=row$bytes)stop("Frozen input missing or size changed: ",path)
  if(!is.null(row$first_date)&&!is.null(row$last_date)) {
    if(min(required_years)<as.integer(substr(row$first_date,1,4))||
       max(required_years)>as.integer(substr(row$last_date,1,4)))stop("Requested years outside frozen input: ",path)
  }
  path
}

#' Define exact model and year unions for experiment-specific common supports
#' @param mode ERA5 or PME profile.
#' @param cfg Shared Article 1 configuration.
#' @return List model_names/years. ERA5 uses its reference plus 22 candidates and
#'   only calibration/evaluation windows; PME uses 15 models, calibration and all
#'   five target periods in one support, never a support per held-out reference.
#' @examples
#' stopifnot(length(support_scope("pme")$model_names)==15L)
support_scope <- function(mode=c("era5","pme"),cfg=article1_config()) {
  mode<-match.arg(mode)
  list(model_names=scan(if(mode=="era5")cfg$era5_models else cfg$pme_models,what="",quiet=TRUE),
    years=if(mode=="era5")sort(unique(c(cfg$calibration_years,cfg$era5_target_years)))else
      sort(unique(c(cfg$calibration_years,unlist(Map(seq.int,cfg$periods$start_year,cfg$periods$end_year))))))
}

#' Require the configured support and verify corrected-data readiness/provenance
#' @param mode ERA5 or PME support profile.
#' @return Support path. With an active registry, requires the readiness receipt,
#'   unchanged registry/support SHA-256 and exact support model/year scope. Without
#'   a registry, only requires the archived support file to exist. Does not itself
#'   freshly hash climate inputs or approve scientific/source-version qualifications.
resolve_experiment_support <- function(mode=c("era5","pme")) {
  mode<-match.arg(mode);cfg<-article1_config();path<-cfg[[paste0(mode,"_support")]]
  if(!file.exists(path))stop("Histogram support unavailable: ",path)
  support_metadata<-readRDS(path)$meta
  if(!identical(support_metadata$precipitation_preprocessing,cfg$precipitation_preprocessing))
    stop("Incompatible precipitation preprocessing in support; regenerate QC supports instead of reusing old bounds.")
  if(nzchar(cfg$input_registry)) {
    readiness<-file.path(dirname(cfg$input_registry),"readiness.json")
    if(!file.exists(readiness))stop("Data/support preparation has not completed; readiness.json missing")
    ready<-jsonlite::fromJSON(readiness)
    if(ready$state!="ready_for_rerun_with_documented_qualifications")stop("Data/support preparation is not ready")
    hash<-mvgc_sha256(cfg$input_registry)
    if(!identical(ready$input_registry_sha256,hash))stop("Input registry changed after preparation")
    obj<-readRDS(path);scope<-support_scope(mode,cfg)
    if(!identical(obj$meta$input_registry_sha256,hash)||
       !identical(as.character(obj$meta$model_names),as.character(scope$model_names))||
       !is.numeric(obj$meta$years)||
       !identical(as.numeric(obj$meta$years),as.numeric(scope$years)))
      stop("Support is not tied to the current registry and exact experiment scope")
    expected<-ready$supports[[mode]]$sha256
    if(!identical(mvgc_sha256(path),expected))stop("Support checksum differs from preparation receipt")
  }
  path
}

#' Record the registry and support identity for a result or checkpoint
#' @param support Path to one histogram-support RDS.
#' @return List containing paths, freshly calculated registry/support SHA-256 and
#'   an input-identity description. Payload hashes are not recalculated here;
#'   directory-only legacy selection is explicitly weaker than a frozen registry.
experiment_input_provenance <- function(support) {
  path<-Sys.getenv("MVGC_INPUT_REGISTRY","")
  list(input_registry=path,input_registry_sha256=if(nzchar(path))mvgc_sha256(path)else NA_character_,
    support_file=support,support_sha256=mvgc_sha256(support),
    input_identity=if(nzchar(path))"Frozen registry; see full payload checks in preparation receipt"else"Legacy directory selection")
}

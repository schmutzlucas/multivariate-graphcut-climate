# Level 1 only: selecting compact results never changes production configuration.

#' Resolve the compact numerical inputs for a named article reproduction profile
#' @param profile Current (completed October reruns)
#'   (the accepted publication cohort). MVGC_REPRODUCTION_PROFILE selects it.
#' @return Named list with profile, paths, environment override names and cohort
#'   manifests. No file is read, directory created or experiment configured.
reproduction_profile <- function(profile = Sys.getenv("MVGC_REPRODUCTION_PROFILE", "current")) {
  if (length(profile) != 1L || is.na(profile) || !profile %in% "current")
    stop("MVGC_REPRODUCTION_PROFILE must be current.")
  paths <- c(
    era5_fields = "results/era5/distributional_fields.rds",
    era5_bias = "results/era5/mean_bias_fields.rds",
    era5_scores = "results/era5/summary_long.csv",
    pme_densities = "results/perfect_model/curated/densities",
    pme_boxplot = "results/perfect_model/curated/boxplot_data.csv",
    pme_scores = "results/perfect_model/evaluation/reference_scores.csv",
    label_maps = "results/spatial_regularisation/label_maps.rds",
    lambda_tradeoff = "results/spatial_regularisation/lambda_tradeoff_metrics.csv",
    histogram = "results/histogram_resolution/cellwise.rds"
  )
  variables <- c(era5_fields = "MVGC_ERA5_FIELDS_FILE", era5_bias = "MVGC_ERA5_BIAS_FIELDS_FILE",
    era5_scores = "MVGC_ERA5_SCORES_FILE", pme_densities = "MVGC_PME_DENSITY_ROOT",
    pme_boxplot = "PME_REFERENCE_DATA_FILE", pme_scores = "MVGC_PME_SCORES_FILE",
    label_maps = "MVGC_LABEL_MAPS_FILE", lambda_tradeoff = "MVGC_LAMBDA_TRADEOFF_FILE", histogram = "MVGC_HISTOGRAM_RESULT_FILE")
  manifests <- paste0("results/", c("era5",
    "perfect_model", "spatial_regularisation",
    "histogram_resolution"), "/manifest.csv")
  list(profile = profile, paths = paths, variables = variables, manifests = manifests)
}

#' Select one Level 1 input while preserving explicit environment overrides
#' @param key One named input from reproduction_profile()$paths.
#' @param profile Resolved compact-result profile, defaulting to the selected profile.
#' @return Path from the named override, or the selected profile. An explicitly
#'   empty override fails; missing files never trigger selection of another cohort.
reproduction_input <- function(key, profile = reproduction_profile()) {
  if (length(key) != 1L || is.na(key) || !key %in% names(profile$paths))
    stop("Unknown reproduction input key.")
  path <- Sys.getenv(profile$variables[[key]], unset = NA_character_)
  if (is.na(path)) return(profile$paths[[key]])
  if (!nzchar(trimws(path))) stop("Empty reproduction override: ", profile$variables[[key]])
  path
}

#' Verify and inventory the compact inputs selected for Level 1 rendering
#' @param scripts Nonempty selected figure/table entry-point paths.
#' @param profile Resolved current publication compact-result profile.
#' @return Data frame of exact selected input paths, sizes, SHA-256, overrides and
#'   manifest matches. Default inputs must match their archived manifests. Custom
#'   overrides are hashed but not certified as the publication cohort. Climate
#'   inputs, supports, checkpoints and large production workspaces are never read.
reproduction_inventory <- function(scripts, profile = reproduction_profile()) {
  require_packages("digest")
  mapping <- list(era5_distributional = "era5_fields", era5_mean_bias = "era5_bias",
    pme_overlays = "pme_densities", pme_reference_improvements = "pme_boxplot",
    spatial_regularisation = "label_maps", lambda_tradeoff = "lambda_tradeoff", histogram_resolution = "histogram",
    render = c("era5_scores", "pme_scores"))
  stems <- tools::file_path_sans_ext(basename(scripts))
  if (!length(stems) || any(!stems %in% names(mapping))) stop("Unknown reproduction script.")
  keys <- unique(unlist(mapping[stems], use.names = FALSE))
  manifests <- do.call(rbind, lapply(profile$manifests, function(path) {
    m <- read.csv(path, stringsAsFactors = FALSE)
    data.frame(path = if ("file" %in% names(m)) file.path(dirname(path), m$file) else m$filename,
      bytes = as.numeric(m$bytes), sha256 = m$sha256, stringsAsFactors = FALSE)
  }))
  rows <- list()
  for (key in keys) {
    path <- reproduction_input(key, profile)
    overridden <- !is.na(Sys.getenv(profile$variables[[key]], unset = NA_character_))
    if (key == "pme_densities") {
      required <- c("full_overlay", paste0("ldr_", c("0p005", "0p01", "0p025", "0p05", "0p10"), "_overlay"),
        paste0("reference_", c("1976_2000", "2001_2025", "2026_2050", "2051_2075", "2076_2100")))
      selected <- Sys.getenv("MVGC_PME_FIGURES", "all")
      if (selected != "all") {
        requested <- strsplit(selected, ",", fixed = TRUE)[[1]]
        if (!length(requested) || any(!requested %in% required)) stop("Unknown PME figure selection.")
        required <- required[required %in% requested]
      }
      paths <- file.path(path, paste0(required, ".rds"))
    } else paths <- path
    for (input in paths) {
      if (!file.exists(input) || dir.exists(input)) stop("Missing compact reproduction input: ", input)
      bytes <- file.info(input)$size
      hash <- digest::digest(file = input, algo = "sha256")
      target <- normalizePath(input, winslash = "/", mustWork = TRUE)
      known <- normalizePath(manifests$path, winslash = "/", mustWork = FALSE)
      match <- which(known == target & manifests$bytes == bytes & manifests$sha256 == hash)
      if (!overridden && !length(match)) stop("Compact input does not match its cohort manifest: ", input)
      rows[[length(rows) + 1L]] <- data.frame(profile = profile$profile, key = key,
        environment_variable = profile$variables[[key]], explicit_override = overridden,
        path = target, bytes = bytes, sha256 = hash,
        verification = if (length(match)) "matches_profile_manifest" else "custom_override_not_publication_verified",
        stringsAsFactors = FALSE)
    }
  }
  if ("render" %in% stems) {
    metadata <- c("data/metadata/models.csv", "data/metadata/model_references.bib",
                  "config/era5_models.txt", "config/pme_models.txt")
    approved <- read.csv("publication/release_manifest.csv",stringsAsFactors=FALSE)
    for (input in metadata) {
      hit <- which(approved$path == input)
      if (length(hit) != 1L || !file.exists(input)) stop("Missing table metadata identity: ",input)
      bytes <- file.info(input)$size
      hash <- digest::digest(file=input,algo="sha256")
      if (bytes != approved$bytes[hit] || hash != approved$sha256[hit])
        stop("Table metadata does not match the publication manifest: ",input)
      rows[[length(rows)+1L]] <- data.frame(profile=profile$profile,key=basename(input),
        environment_variable="",explicit_override=FALSE,
        path=normalizePath(input,winslash="/",mustWork=TRUE),bytes=bytes,sha256=hash,
        verification="matches_profile_manifest",stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

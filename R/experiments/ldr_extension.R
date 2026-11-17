# Legacy PME diagnostic extension. These definitions do not start an experiment.
#' Validate saved PME identities and labels before extending LDR diagnostics
#' @param bundle Compact PME list with meta, periods, GC_result_list and existing
#'   ordinary/LDR result collections. A complete one-percent baseline is required.
#' @param cfg Shared Article 1 configuration; compact grids/model lists may be
#'   supplied explicitly by tests. Publication defaults remain 15 references/8^3 bins.
#' @return TRUE invisibly after exact configured model/grid/variable/year/period
#'   checks, integer candidate-label checks and finite existing one-percent score checks.
#'   Does not approve climate payloads or recompute scores.
validate_ldr_extension_bundle <- function(bundle, cfg = article1_config()) {
  if (!is.list(bundle) || !is.list(bundle$meta) || !is.data.frame(bundle$periods))
    stop("Expected a compact PME bundle with scientific metadata and periods.")
  m <- bundle$meta; models <- scan(cfg$pme_models, what = "", quiet = TRUE)
  if (length(models) < 3L || anyDuplicated(models) || !identical(m$model_names, models) ||
      !identical(m$variables, cfg$variables) || !identical(m$lon, cfg$lon) || !identical(m$lat, cfg$lat) ||
      !is.numeric(m$nbins3d) || !identical(as.numeric(m$nbins3d), as.numeric(cfg$nbins3d)) ||
      !is.numeric(m$calibration_years) || !identical(as.numeric(m$calibration_years), as.numeric(cfg$calibration_years)) ||
      !identical(bundle$periods$period_id, cfg$periods$period_id) ||
      !is.numeric(bundle$periods$start_year) || !is.numeric(bundle$periods$end_year) ||
      !identical(as.numeric(bundle$periods$start_year), as.numeric(cfg$periods$start_year)) ||
      !identical(as.numeric(bundle$periods$end_year), as.numeric(cfg$periods$end_year)))
    stop("PME extension model/grid/variable/bin/calibration/period identities differ from configuration.")
  if (!is.null(m$run_state) && !identical(m$run_state, "complete")) stop("Cannot extend an incomplete PME run.")
  lambdas <- m$lambda_keys
  if (!is.character(lambdas) || !length(lambdas) || anyNA(lambdas) || anyDuplicated(lambdas) ||
      !"0.125" %in% lambdas || !is.numeric(m$lambda) ||
      any(!is.finite(m$lambda)) || any(m$lambda < 0) || anyDuplicated(m$lambda) ||
      !identical(suppressWarnings(as.numeric(lambdas)), as.numeric(m$lambda)))
    stop("PME extension requires distinct numeric lambda identities including 0.125.")
  grid <- c(length(m$lon), length(m$lat))
  for (pid in bundle$periods$period_id) {
    # Existing 1% scores are a mandatory numerical anchor, even for cache reuse.
    anchor <- list("0.01" = list(MMM = bundle$MMM_partial_hdist_future_list[[pid]][["0.01"]],
      GC = bundle$GC_partial_hdist_future_list[[pid]][["0.01"]]))
    validate_ldr_extension_block(anchor, "0.01", lambdas, models, grid)
    for (lambda in lambdas) {
      fits <- bundle$GC_result_list[[pid]][[lambda]]
      if (!is.list(fits) || anyDuplicated(names(fits)) || !setequal(names(fits), models))
        stop("Missing/ambiguous saved labels: ", pid, "/", lambda)
      for (model in models) {
        labels <- fits[[model]]$label_attribution
        if (!is.numeric(labels) || !identical(dim(labels), grid) || any(!is.finite(labels)) ||
            any(labels != floor(labels)) || any(labels < 1 | labels > length(models) - 1L))
          stop("Invalid saved candidate labels: ", pid, "/", lambda, "/", model)
      }
    }
  }
  invisible(TRUE)
}

#' Validate a complete numerical LDR latitude-block payload
#' @param maps Named threshold list, each with MMM reference maps and GC lambda/reference maps.
#' @param keys,lambda_keys Exact threshold and lambda key vectors expected in this block.
#' @param models Ordered held-out reference names.
#' @param grid Integer longitude/block-latitude dimensions.
#' @return TRUE invisibly; missing/duplicate identities, wrong shapes and nonfinite,
#'   negative or above-one scores raise an error. Missing values are never approved.
validate_ldr_extension_block <- function(maps, keys, lambda_keys, models, grid) {
  #' Validate exact named-list identities inside an LDR block
  #' Internal closure of validate_ldr_extension_block; not a standalone entry point.
  #' @param x Object expected to be a uniquely named list.
  #' @param expected Character set of required names.
  #' @return TRUE for a uniquely named list with exactly the requested identities.
  named <- function(x, expected) is.list(x) && !is.null(names(x)) &&
    !anyNA(names(x)) && !anyDuplicated(names(x)) && setequal(names(x), expected)
  #' Validate an LDR field against the enclosing grid
  #' Internal closure of validate_ldr_extension_block; not a standalone entry point.
  #' @param x Numeric matrix to validate.
  #' @return TRUE for finite nonnegative scores no larger than one within 1e-8 on the enclosing grid.
  valid <- function(x) is.numeric(x) && is.matrix(x) && identical(dim(x), as.integer(grid)) &&
    all(is.finite(x)) && all(x >= 0 & x <= 1 + 1e-8)
  if (!named(maps, keys)) stop("LDR checkpoint has missing/ambiguous threshold keys.")
  for (key in keys) {
    if (!named(maps[[key]]$MMM, models) || !named(maps[[key]]$GC, lambda_keys) ||
        !all(vapply(maps[[key]]$MMM, valid, logical(1)))) stop("Invalid LDR baseline block: ", key)
    for (lambda in lambda_keys) if (!named(maps[[key]]$GC[[lambda]], models) ||
        !all(vapply(maps[[key]]$GC[[lambda]], valid, logical(1))))
      stop("Invalid LDR GC block: ", key, "/", lambda)
  }
  invisible(TRUE)
}

#' Evaluate reference-defined LDR scores for all saved lambda label maps
#' @param pmfs Numeric array bins x cells x models of finite nonnegative PMFs
#'   normalised per cell/model on the single common PME support.
#' @param labels Named lambda list of named reference matrices on this latitude block;
#'   labels index the ordered candidates obtained after excluding that reference.
#' @param thresholds Named numeric ell_LDR values strictly in (0,1).
#' @param models Ordered reference/model identifiers matching the PMF model axis.
#' @param lon_size Positive integer number of longitude cells; cells must divide exactly.
#' @return Named threshold list with MMM reference matrices and GC lambda/reference
#'   matrices. Uses equal-model means, all zero-reference/cutoff-tied bins and the
#'   non-renormalised Hellinger formula. No labels are optimised or modified.
compute_ldr_extension_block <- function(pmfs, labels, thresholds, models, lon_size) {
  dims <- dim(pmfs)
  if (!is.numeric(pmfs) || length(dims) != 3L || any(dims < 1L) || dims[3] != length(models) ||
      length(models) < 3L || anyDuplicated(models) || any(!is.finite(pmfs)) || any(pmfs < 0) ||
      !is.numeric(lon_size) || length(lon_size) != 1L || !is.finite(lon_size) ||
      lon_size < 1 || lon_size != floor(lon_size) || dims[2] %% lon_size != 0 ||
      !is.numeric(thresholds) || !length(thresholds) || any(!is.finite(thresholds)) ||
      any(thresholds <= 0 | thresholds >= 1) || is.null(names(thresholds)) || anyDuplicated(names(thresholds)))
    stop("Invalid PMF/threshold/model/block dimensions for LDR evaluation.")
  ncell <- dims[2]; grid <- c(as.integer(lon_size), as.integer(ncell / lon_size))
  for (m in seq_along(models)) if (any(abs(colSums(matrix(pmfs[, , m], dims[1], ncell)) - 1) > 1e-8))
    stop("LDR input PMFs must sum to one in every cell/model.")
  if (!is.list(labels) || !length(labels) || is.null(names(labels)) || anyDuplicated(names(labels)))
    stop("Expected distinct named lambda label collections.")
  for (lambda in names(labels)) {
    if (!is.list(labels[[lambda]]) || anyDuplicated(names(labels[[lambda]])) ||
        !setequal(names(labels[[lambda]]), models)) stop("Missing reference labels for lambda ", lambda)
    for (model in models) {
      x <- labels[[lambda]][[model]]
      if (!is.numeric(x) || !identical(dim(x), grid) || any(!is.finite(x)) ||
          any(x != floor(x)) || any(x < 1 | x > length(models) - 1L)) stop("Invalid LDR candidate labels.")
    }
  }
  maps <- setNames(lapply(thresholds, function(tau) list(MMM = setNames(vector("list", length(models)), models),
    GC = setNames(lapply(labels, function(x) setNames(vector("list", length(models)), models)), names(labels)))), names(thresholds))
  for (m in seq_along(models)) {
    ref <- matrix(pmfs[, , m], dims[1], ncell); candidates <- setdiff(seq_along(models), m)
    mmm <- apply(pmfs[, , candidates, drop = FALSE], c(1, 2), mean)
    selected <- lapply(labels, function(x) candidates[as.vector(x[[models[m]]])])
    for (key in names(thresholds)) {
      h_mmm <- numeric(ncell); h_gc <- matrix(NA_real_, ncell, length(labels))
      for (cell in seq_len(ncell)) {
        bins <- select_ldr_indices(ref[, cell], thresholds[[key]])
        if (!length(bins)) stop("LDR reference bin selection is empty.")
        h_mmm[cell] <- sqrt(sum((sqrt(mmm[bins, cell]) - sqrt(ref[bins, cell]))^2)) / sqrt(2)
        for (l in seq_along(labels)) h_gc[cell, l] <- sqrt(sum((sqrt(pmfs[bins, cell, selected[[l]][cell]]) -
          sqrt(ref[bins, cell]))^2)) / sqrt(2)
      }
      maps[[key]]$MMM[[models[m]]] <- matrix(h_mmm, grid[1], grid[2])
      for (l in seq_along(labels)) maps[[key]]$GC[[l]][[models[m]]] <- matrix(h_gc[, l], grid[1], grid[2])
    }
  }
  validate_ldr_extension_block(maps, names(thresholds), names(labels), models, grid)
  maps
}

#' Bind LDR extension caches to saved results, support, settings and input discovery
#' @param input,range_path Existing source-bundle and common-support file paths.
#' @param bundle Validated source bundle, including saved labels and scientific settings.
#' @param data_dir Actual prepared data root used for NetCDF discovery.
#' @param cfg Shared configuration, including the requested LDR threshold set.
#' @return Signature list with fresh bundle/support/registry hashes, exact scientific
#'   settings, canonical code/native hashes and selected input file size/mtime records.
#'   Raw climate payloads are not rehashed here; legacy directory identity is
#'   explicitly weaker than a frozen registry plus the prepared-run payload checks.
ldr_extension_signature <- function(input, range_path, bundle, data_dir, cfg = article1_config()) {
  years <- sort(unique(unlist(lapply(seq_len(nrow(bundle$periods)), function(i)
    seq.int(bundle$periods$start_year[i], bundle$periods$end_year[i])))))
  if (nzchar(cfg$input_registry))
    read_input_registry(cfg$input_registry, required_models = bundle$meta$model_names, required_variables = bundle$meta$variables)
  records <- lapply(bundle$meta$model_names, function(model) lapply(bundle$meta$variables, function(variable) {
    path <- find_matching_netcdf_file(data_dir, model, variable, years); info <- file.info(path)
    list(model = model, variable = variable, path = normalizePath(path, winslash = "/", mustWork = TRUE),
      bytes = unname(info$size), mtime = as.numeric(info$mtime))
  }))
  list(schema = "article1-ldr-extension-v2", bundle_sha256 = mvgc_sha256(input), support_sha256 = mvgc_sha256(range_path),
    registry_sha256 = if (nzchar(cfg$input_registry)) mvgc_sha256(cfg$input_registry) else NA_character_,
    input_identity = if (nzchar(cfg$input_registry)) "Frozen registry plus selected-file metadata; raw payloads not rehashed here" else
      "Legacy directory selection: file path/size/mtime only; no raw-payload checksum approval",
    input_records = records, models = bundle$meta$model_names, variables = bundle$meta$variables,
    lon = bundle$meta$lon, lat = bundle$meta$lat, nbins3d = bundle$meta$nbins3d,
    periods = bundle$periods, thresholds = cfg$pme_ldr_thresholds, lambda_keys = bundle$meta$lambda_keys,
    selection = bundle$meta[c("season_start_md", "season_end_md", "strict_time_checks", "season_selection_mode",
      "time_decode_mode", "requested_year_policy", "partial_h_definition")],
    source_sha256 = scientific_code_fingerprint(c("R/experiments/ldr_extension.R", "experiments/sensitivity/ldr_threshold.R",
      "R/distributions/bin_sensitivity.R", "src/bin_sensitivity.cpp", "R/preprocessing/calendar.R",
      "R/preprocessing/cf_calendar_validation.R", "R/preprocessing/netcdf_io.R", "R/metrics/ldr.R")),
    R_version = as.character(getRversion()), platform = R.version$platform,
    software_versions = vapply(c("Rcpp", "ncdf4"), function(x) as.character(utils::packageVersion(x)), ""))
}

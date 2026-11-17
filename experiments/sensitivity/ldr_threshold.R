# Sensitivity entry point: changing the tested resolution or regularisation
# changes the scientific configuration. Keep external inputs/supports fixed when
# comparing settings. Resource controls affect execution only; retain distinct
# output/checkpoint namespaces to avoid mixing incompatible experiments.
# Extend diagnostics from existing PME labels; never optimise graph cuts here.
source("R/load.R"); load_mvgc()
source("R/experiments/ldr_extension.R")
cfg <- article1_config()
input <- production_result_file("PME_FINAL_RESULTS_FILE", "the final PME compact bundle")
output <- Sys.getenv("PME_LDR_OUTPUT_FILE", unset = file.path(mvgc_output_dir("workspaces"), "pme_final_periods_extended_ldr_compact.rds"))
if (file.exists(output)) stop("Output exists; choose a new PME_LDR_OUTPUT_FILE instead of overwriting a result.")
bundle <- readRDS(input)
validate_ldr_extension_bundle(bundle, cfg)
if (nzchar(cfg$input_registry)) {
  current_support <- resolve_experiment_support("pme")
  if (!identical(bundle$meta$input_provenance, experiment_input_provenance(current_support)))
    stop("Saved PME result does not belong to the current corrected input/support preparation.")
}
data_dir <- Sys.getenv("MVGC_DATA_ROOT", unset = "")
if (!nzchar(data_dir)) data_dir <- bundle$meta$data_dir
if (!is.character(data_dir) || length(data_dir) != 1L || !dir.exists(data_dir))
  stop("Source NetCDF storage is unavailable; set MVGC_DATA_ROOT.")
lon <- bundle$meta$lon; lat <- bundle$meta$lat; models <- bundle$meta$model_names
variables <- bundle$meta$variables; lambda_keys <- bundle$meta$lambda_keys
range_path <- canonical_support_path(bundle$meta$range_sources)
if (length(range_path) != 1L || is.na(range_path)) stop("Expected one common PME support file.")
if (!file.exists(range_path)) range_path <- file.path(Sys.getenv("MVGC_SOURCE_WORKSPACE", unset = dirname(dirname(input))), range_path)
if (!file.exists(range_path)) stop("Saved PMF range file is unavailable: ", range_path)
range_source <- readRDS(range_path)
lon_order <- match(ifelse(lon < 0, lon + 360, lon), 0:359)
lat_order <- match(lat, -90:90)
if (anyNA(lon_order) || anyNA(lat_order)) stop("LDR extension requires the common one-degree support grid.")
ranges <- array(NA_real_, c(length(lon), length(lat), length(variables), 2L))
for (v in seq_along(variables)) {
  bounds <- range_source$ranges[[variables[v]]]
  if (!is.numeric(bounds) || !identical(dim(bounds), c(360L, 181L, 2L))) stop("Malformed common support bounds.")
  ranges[, , v, ] <- bounds[lon_order, lat_order, , drop = FALSE]
}
tau_values <- cfg$pme_ldr_thresholds
tau_keys <- vapply(tau_values, function(x) format(x, trim = TRUE, scientific = FALSE, decimal.mark = ".", nsmall = 2), "")
signature <- ldr_extension_signature(input, range_path, bundle, data_dir, cfg)
require_packages("openssl")
identity_hash <- as.character(openssl::sha256(serialize(signature, NULL, version = 2)))
checkpoint_dir <- file.path(mvgc_output_dir("checkpoints/pme_ldr"), paste0("identity_", substr(identity_hash, 1, 24)))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
checked_periods <- character()

for (period_id in bundle$periods$period_id) {
  row <- bundle$periods[bundle$periods$period_id == period_id, , drop = FALSE]
  mmm_existing <- bundle$MMM_partial_hdist_future_list[[period_id]]
  gc_existing <- bundle$GC_partial_hdist_future_list[[period_id]]
  if (length(setdiff(union(names(mmm_existing), names(gc_existing)), tau_keys)))
    stop("Source has LDR thresholds outside the configured set; preserve/review them explicitly.")
  # Existing collections may be absent for a whole threshold/lambda, but an
  # existing collection with missing references/invalid values is corruption.
  for (key in union(names(mmm_existing), names(gc_existing))) {
    collections <- c(if (!is.null(mmm_existing[[key]])) list(mmm_existing[[key]]) else list(), gc_existing[[key]])
    for (collection in collections) {
      if (!is.list(collection) || anyDuplicated(names(collection)) || !setequal(names(collection), models))
        stop("Incomplete existing LDR reference collection: ", period_id, "/", key)
      for (model in models) {
        x <- collection[[model]]
        if (!is.numeric(x) || !identical(dim(x), c(length(lon), length(lat))) ||
            any(!is.finite(x)) || any(x < 0 | x > 1 + 1e-8))
          stop("Invalid existing LDR score map: ", period_id, "/", key, "/", model)
      }
    }
    if (anyDuplicated(names(gc_existing[[key]])) || length(setdiff(names(gc_existing[[key]]), lambda_keys)))
      stop("Unknown/duplicate lambda in existing LDR results: ", period_id, "/", key)
  }
  missing_mmm <- tau_keys[vapply(tau_keys, function(key) is.null(mmm_existing[[key]]), logical(1))]
  missing_gc <- setNames(lapply(tau_keys, function(key) setdiff(lambda_keys, names(gc_existing[[key]]))), tau_keys)
  missing_keys <- union(missing_mmm, tau_keys[lengths(missing_gc) > 0L])
  if (!length(missing_keys)) next
  calculate_keys <- unique(c("0.01", missing_keys))
  plans <- prepare_k_sensitivity_plans(models, variables, data_dir, row$start_year:row$end_year, lon, lat,
    season_start_md = bundle$meta$season_start_md, season_end_md = bundle$meta$season_end_md,
    strict_time_checks = bundle$meta$strict_time_checks, season_selection_mode = bundle$meta$season_selection_mode,
    time_decode_mode = bundle$meta$time_decode_mode, requested_year_policy = bundle$meta$requested_year_policy)
  for (key in missing_mmm) bundle$MMM_partial_hdist_future_list[[period_id]][[key]] <-
    setNames(lapply(models, function(x) matrix(NA_real_, length(lon), length(lat))), models)
  for (key in missing_keys) for (lambda in missing_gc[[key]])
    bundle$GC_partial_hdist_future_list[[period_id]][[key]][[lambda]] <-
      setNames(lapply(models, function(x) matrix(NA_real_, length(lon), length(lat))), models)
  blocks <- split(seq_along(lat), ceiling(seq_along(lat) / 4L))
  load_pmf_sensitivity_cpp()
  for (idx in blocks) {
    checkpoint <- file.path(checkpoint_dir, sprintf("%s_lat%03d-%03d.rds", period_id, min(idx), max(idx)))
    block_identity <- list(signature = signature, period = row, latitude_indices = idx, calculate_keys = calculate_keys)
    reused <- file.exists(checkpoint)
    if (reused) {
      record <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
      if (!is.list(record) || !identical(record$identity, block_identity))
        stop("Incompatible or corrupt LDR checkpoint: ", checkpoint)
      maps <- record$maps
    } else {
      cat(period_id, "latitude rows", min(idx), "to", max(idx), "\n")
      ncell <- length(lon) * length(idx); bins <- as.integer(bundle$meta$nbins3d^3)
      pmfs <- array(NA_real_, c(bins, ncell, length(models)))
      for (m in seq_along(models)) {
        raw <- read_k_sensitivity_block(plans[[models[m]]], variables, idx)
        counts <- multik_hist_counts_cpp(raw$pr, raw$tas, raw$psl, ranges[, idx, , , drop = FALSE],
          rep(TRUE, dim(raw$pr)[3]), as.integer(bundle$meta$nbins3d))$counts
        counts <- matrix(counts[seq_len(bins * ncell)], nrow = bins)
        total <- colSums(counts)
        if (any(!is.finite(total)) || any(total <= 0) || any(!is.finite(counts)) || any(counts < 0))
          stop("LDR histogram has invalid counts or zero complete observations.")
        pmfs[, , m] <- sweep(counts, 2, total, "/")
        rm(raw, counts); gc(FALSE)
      }
      labels <- setNames(lapply(lambda_keys, function(lambda) setNames(lapply(models, function(model)
        bundle$GC_result_list[[period_id]][[lambda]][[model]]$label_attribution[, idx, drop = FALSE]), models)), lambda_keys)
      maps <- compute_ldr_extension_block(pmfs, labels, setNames(tau_values[match(calculate_keys, tau_keys)], calculate_keys),
        models, length(lon))
      rm(pmfs); gc(FALSE)
    }
    validate_ldr_extension_block(maps, calculate_keys, lambda_keys, models, c(length(lon), length(idx)))
    # This check is deliberately outside the cache branch: reusing a checkpoint
    # never bypasses its numerical anchor or comparisons with other saved scores.
    for (key in calculate_keys) for (model in models) {
      if (!is.null(mmm_existing[[key]]) && max(abs(maps[[key]]$MMM[[model]] - mmm_existing[[key]][[model]][, idx, drop = FALSE])) > 1e-12)
        stop("Existing D-MMM LDR diagnostic differs: ", period_id, "/", key, "/", model)
      for (lambda in intersect(lambda_keys, names(gc_existing[[key]])))
        if (max(abs(maps[[key]]$GC[[lambda]][[model]] - gc_existing[[key]][[lambda]][[model]][, idx, drop = FALSE])) > 1e-12)
          stop("Existing GC LDR diagnostic differs: ", period_id, "/", key, "/", lambda, "/", model)
    }
    if (!reused) {
      temporary <- tempfile("ldr_block_", tmpdir = checkpoint_dir, fileext = ".rds")
      saveRDS(list(identity = block_identity, maps = maps), temporary, compress = "xz")
      if (!file.rename(temporary, checkpoint)) stop("Could not publish LDR checkpoint.")
    }
    for (key in missing_mmm) for (model in models)
      bundle$MMM_partial_hdist_future_list[[period_id]][[key]][[model]][, idx] <- maps[[key]]$MMM[[model]]
    for (key in missing_keys) for (lambda in missing_gc[[key]]) for (model in models)
      bundle$GC_partial_hdist_future_list[[period_id]][[key]][[lambda]][[model]][, idx] <- maps[[key]]$GC[[lambda]][[model]]
  }
  checked_periods <- c(checked_periods, period_id)
}
for (pid in bundle$periods$period_id) {
  maps <- setNames(lapply(tau_keys, function(key) list(MMM = bundle$MMM_partial_hdist_future_list[[pid]][[key]],
    GC = bundle$GC_partial_hdist_future_list[[pid]][[key]])), tau_keys)
  validate_ldr_extension_block(maps, tau_keys, lambda_keys, models, c(length(lon), length(lat)))
}
if (!identical(signature, ldr_extension_signature(input, range_path, readRDS(input), data_dir, cfg)))
  stop("Source bundle/support/input metadata or scientific code changed during LDR extension.")
bundle$meta$partial_tau <- tau_values; bundle$meta$partial_tau_keys <- tau_keys
bundle$meta$ldr_extension <- list(source_bundle = basename(input), source_sha256 = signature$bundle_sha256,
  scientific_settings = signature, lambda_keys = lambda_keys, reused_graphcut_labels = TRUE,
  checked_existing_tau = if (length(checked_periods)) .01 else NULL, checked_periods = checked_periods,
  regression_tolerance = 1e-12, complete_threshold_lambda_coverage = TRUE, completed_at = as.character(Sys.time()))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) stop("Output appeared during processing; choose a new PME_LDR_OUTPUT_FILE.")
temporary <- tempfile("pme_ldr_", tmpdir = dirname(output), fileext = ".rds")
saveRDS(bundle, temporary, compress = "xz")
if (!file.rename(temporary, output)) stop("Could not publish extended PME bundle.")
cat("Saved", output, "\n")

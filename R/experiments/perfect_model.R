#' Read a positive-integer environment setting for PME orchestration
#' @param var_name Environment-variable name.
#' @param default_value Positive integer fallback when unset or invalid.
#' @return Integer, warning and using the default for fractional/nonfinite/
#'   nonpositive/out-of-range values. This helper does not choose scientific parameters.
read_positive_int_env <- function(var_name, default_value) {
  raw_value <- Sys.getenv(var_name, unset = "")
  if (raw_value == "") {
    return(as.integer(default_value))
  }

  parsed <- suppressWarnings(as.numeric(raw_value))
  if (!is.finite(parsed) || parsed < 1L || parsed != floor(parsed) ||
      parsed > .Machine$integer.max) {
    warning(
      "Ignoring invalid value for ",
      var_name,
      ": '",
      raw_value,
      "'. Using default ",
      default_value,
      ".",
      call. = FALSE
    )
    return(as.integer(default_value))
  }

  as.integer(parsed)
}

#' Format a graph-cut lambda for the archived PME nested-list keys
#' @param value Numeric smoothness value, normally a scalar from the configured grid.
#' @return Fixed-notation character key without padding, e.g. "0.125". Not a
#'   numerical fingerprint; complete scientific settings separately identify caches.
lambda_to_key <- function(value) {
  format(as.numeric(value), trim = TRUE, scientific = FALSE)
}

#' Format an LDR threshold for the archived PME nested-list keys
#' @param value Numeric ell_LDR threshold (legacy argument name tau).
#' @return Fixed-notation character key with at least two decimal places, e.g.
#'   "0.01", "0.10" or "0.005". This does not change the threshold used in scoring.
tau_to_key <- function(value) {
  format(as.numeric(value), trim = TRUE, scientific = FALSE, nsmall = 2)
}

#' Make an ASCII model identifier suitable for a checkpoint filename
#' @param value Character name or vector of names.
#' @return Names replacing non-letter/digit/underscore/hyphen characters with
#'   underscores. Not an identity-preserving encoding; the checkpoint metadata
#'   retains exact model names and ordered model indices.
sanitize_name <- function(value) {
  gsub("[^A-Za-z0-9_-]", "_", value)
}

#' Read the common PME support and reorder longitude to the graph-cut grid
#' @param range_file_path RDS with a ranges list in 0:359 longitude / -90:90 latitude order.
#' @param variables Requested ordered variable identifiers, normally pr/tas/psl.
#' @return Array [360,181,variables,2] in -180:179 longitude order with dimnames.
#'   Reads existing bounds, never estimates per-reference supports. Requires each
#'   variable entry to be exactly [360,181,2]; later PMF validation checks finite,
#'   strictly increasing bounds and the prepared-support resolver checks provenance.
build_range_cube <- function(range_file_path, variables) {
  range_obj <- readRDS(range_file_path)

  if (!is.list(range_obj) || is.null(range_obj$ranges)) {
    stop("Range file has unexpected structure: ", range_file_path)
  }

  if (!all(variables %in% names(range_obj$ranges))) {
    missing_vars <- variables[!(variables %in% names(range_obj$ranges))]
    stop(
      "Missing variables in range file ",
      range_file_path,
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  lon_file <- 0:359
  lon_adj <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
  lon_order_adj <- order(lon_adj)

  range_cube <- array(
    NA_real_,
    dim = c(360, 181, length(variables), 2),
    dimnames = list(
      lon = -180:179,
      lat = -90:90,
      var = variables,
      bound = c("min", "max")
    )
  )

  for (k in seq_along(variables)) {
    var_name <- variables[k]
    mat <- range_obj$ranges[[var_name]]

    if (!identical(dim(mat), c(360L, 181L, 2L))) {
      stop(
        "Unexpected dimension for variable ",
        var_name,
        " in ",
        range_file_path,
        ". Expected [360, 181, 2], got [",
        paste(dim(mat), collapse = ", "),
        "]."
      )
    }

    range_cube[, , k, ] <- mat[lon_order_adj, , , drop = FALSE]
  }

  range_cube
}

#' Combine two legacy support cubes using local minima and maxima
#' @param range_a,range_b Identically shaped [lon,lat,variable,2] support arrays.
#' @return Cube of the same shape/dimnames, taking the available endpoint when
#'   only one is NA. Errors if combined min exceeds max; does not approve remaining
#'   missing/zero-width bounds. Retained for archived workflows, not used to build
#'   separate held-out supports in the current common-support PME.
combine_two_ranges <- function(range_a, range_b) {
  if (!identical(dim(range_a), dim(range_b))) {
    stop("Range cubes must have identical dimensions before combining.")
  }

  min_a <- range_a[, , , 1]
  min_b <- range_b[, , , 1]
  max_a <- range_a[, , , 2]
  max_b <- range_b[, , , 2]

  combined_min <- ifelse(
    is.na(min_a), min_b,
    ifelse(is.na(min_b), min_a, pmin(min_a, min_b))
  )
  combined_max <- ifelse(
    is.na(max_a), max_b,
    ifelse(is.na(max_b), max_a, pmax(max_a, max_b))
  )

  out <- array(
    NA_real_,
    dim = dim(range_a),
    dimnames = dimnames(range_a)
  )
  out[, , , 1] <- combined_min
  out[, , , 2] <- combined_max

  invalid <- which(out[, , , 1] > out[, , , 2], arr.ind = TRUE)
  if (nrow(invalid) > 0L) {
    stop("Found ", nrow(invalid), " grid cells where combined min > combined max.")
  }

  out
}

#' Average candidate PMFs with equal model weight to form D-MMM
#' @param pdf_models Normalised nonnegative finite PMFs, lon x lat x bins x models.
#' @return Lon x lat x bins arithmetic mean. Models, not their sample counts,
#'   have equal weight. Invalid candidate PMFs raise an error before averaging.
#' @examples
#' compute_mmm_pdf(array(c(.2,.8,.4,.6),c(1,1,2,2)))
compute_mmm_pdf <- function(pdf_models) {
  validate_target_pmfs(pdf_models)
  apply(pdf_models, c(1, 2, 3), mean)
}

#' Validate matched reference and D-MMM PMF arrays
#' @param pdf_ref,pdf_d_mmm Numeric lon x lat x bins arrays on identical supports.
#' @return TRUE invisibly; invalid shape, masses or normalisation raise an error.
#'   No missing value replacement, renormalisation or score adjustment is performed.
validate_mmm_comparison <- function(pdf_ref, pdf_d_mmm) {
  d <- dim(pdf_ref)
  if (!is.numeric(pdf_ref) || length(d) != 3L || any(d < 1L) ||
      !is.numeric(pdf_d_mmm) || !identical(d, dim(pdf_d_mmm))) {
    stop("Reference and D-MMM must have matching longitude x latitude x bins dimensions.")
  }
  for (p in list(pdf_ref, pdf_d_mmm)) {
    if (any(!is.finite(p)) || any(p < 0) || any(abs(rowSums(p,dims=2L)-1) > 1e-8)) {
      stop("Reference and D-MMM must contain finite nonnegative PMFs summing to one.")
    }
  }
  invisible(TRUE)
}

#' Evaluate ordinary Hellinger distance from a reference to D-MMM
#' @param pdf_ref Reference PMFs, lon x lat x bins, normalised on the common support.
#' @param pdf_models Candidate PMFs, lon x lat x bins x models, used if mmm_pdf is NULL.
#' @param mmm_pdf Optional precomputed equal-model PMF average. Supplying it avoids
#'   recomputing the same mixture for each diagnostic/threshold; caller owns model identity.
#' @return Lon x lat distance matrix, retaining singleton spatial axes. Invalid
#'   masses fail rather than converting a NaN distance to a perfect score of zero.
compute_mmm_hdist <- function(pdf_ref, pdf_models, mmm_pdf = NULL) {
  lon_size <- dim(pdf_ref)[1]
  lat_size <- dim(pdf_ref)[2]

  if (is.null(mmm_pdf)) mmm_pdf <- compute_mmm_pdf(pdf_models)
  validate_mmm_comparison(pdf_ref, mmm_pdf)
  hdist <- matrix(NA_real_, nrow = lon_size, ncol = lat_size)

  for (i in seq_len(lon_size)) {
    for (j in seq_len(lat_size)) {
      hdist[i, j] <- sqrt(sum((sqrt(mmm_pdf[i, j, ]) - sqrt(pdf_ref[i, j, ]))^2)) / sqrt(2)
    }
  }

  hdist
}

#' Evaluate non-renormalised LDR Hellinger score from a reference to D-MMM
#' @param pdf_ref Reference PMFs, lon x lat x bins, normalised on the common support.
#' @param pdf_models Candidate PMFs, lon x lat x bins x models, used if mmm_pdf is NULL.
#' @param selected_indices Nested list indexed longitude then latitude, containing
#'   valid one-based reference-selected bin indices. Empty selections are not valid here.
#' @param mmm_pdf Optional precomputed equal-model average, as in compute_mmm_hdist().
#' @return Lon x lat score matrix with no renormalisation on selected bins. Invalid
#'   inputs/nonfinite scores raise an error; singleton spatial dimensions are retained.
compute_mmm_partial_hdist <- function(pdf_ref, pdf_models, selected_indices, mmm_pdf = NULL) {
  if (is.null(mmm_pdf)) mmm_pdf <- compute_mmm_pdf(pdf_models)
  validate_mmm_comparison(pdf_ref, mmm_pdf)
  partial <- compute_partial_hdist(
    pdf_ref = pdf_ref,
    pdf_model = mmm_pdf,
    selected_indices = selected_indices
  )
  partial <- matrix(partial, nrow=dim(pdf_ref)[1], ncol=dim(pdf_ref)[2])
  if (any(!is.finite(partial))) stop("D-MMM LDR scores are nonfinite; check selected bins and PMFs.")
  partial
}

#' Extract calibration/evaluation score maps for graph-cut-selected models
#' @param gc_result Solver list with integer-valued one-based label_attribution matrix.
#' @param h_dist_pres,h_dist_fut,partial_fut Matching [lon,lat,model] calibration,
#'   ordinary evaluation and LDR evaluation score arrays, respectively.
#' @return List of lon x lat matrices present, future and partial_future. Does not
#'   reoptimise labels; each cell uses its selected candidate's already computed
#'   score. Invalid label identity/shape raises an error instead of leaving NA cells.
compute_gc_selected_maps <- function(gc_result, h_dist_pres, h_dist_fut, partial_fut) {
  d<-dim(h_dist_pres);labels<-gc_result$label_attribution
  if(length(d)!=3L||!identical(d,dim(h_dist_fut))||!identical(d,dim(partial_fut))||
     !is.numeric(labels)||!identical(dim(labels),d[1:2])||any(!is.finite(labels))||
     any(labels!=floor(labels))||any(labels<1|labels>d[3]))stop("Invalid graph-cut labels or score array dimensions.")
  lon_size <- dim(h_dist_pres)[1]
  lat_size <- dim(h_dist_pres)[2]
  n_models <- dim(h_dist_pres)[3]

  out_pres <- matrix(NA_real_, nrow = lon_size, ncol = lat_size)
  out_fut <- matrix(NA_real_, nrow = lon_size, ncol = lat_size)
  out_partial_fut <- matrix(NA_real_, nrow = lon_size, ncol = lat_size)

  for (l in seq_len(n_models)) {
    mask <- gc_result$label_attribution == l
    out_pres[mask] <- h_dist_pres[, , l][mask]
    out_fut[mask] <- h_dist_fut[, , l][mask]
    out_partial_fut[mask] <- partial_fut[, , l][mask]
  }

  list(
    present = out_pres,
    future = out_fut,
    partial_future = out_partial_fut
  )
}

#' Select tie-inclusive reference LDR bins independently at every grid cell
#' @param pdf_ref Valid normalised reference PMFs [lon,lat,bins] on the common support.
#' @param tau Nominal ell_LDR mass threshold, passed to select_ldr_indices().
#' @param lon_size,lat_size Spatial sizes matching the reference array.
#' @return Nested list indexed longitude then latitude of one-based selected bin
#'   indices. Includes zero-probability bins and all cutoff ties. Uses reference
#'   probabilities only; no selected-mass renormalisation is performed.
build_ldr_indices <- function(pdf_ref, tau, lon_size, lat_size) {
  lapply(seq_len(lon_size), function(i) {
    lapply(seq_len(lat_size), function(j) {
      select_ldr_indices(pdf_ref[i, j, ], tau = tau)
    })
  })
}

#' Return the five configured PME target periods without a duplicate date table
#' @return The period_id/label/start_year/end_year data frame from article1_config().
#'   Does not alter calibration or targets; current windows remain 1976--2100.
make_periods <- function() article1_config()$periods

#' Select an explicit nonempty subset of PME periods for a limited run
#' @param periods Configured period table with unique period_id values.
#' @return All rows when PME_PERIODS is unset/empty, otherwise rows in the
#'   comma-separated requested order. Unknown, duplicate or blank-only selections
#'   raise an error; a limited run must not silently execute zero/repeated periods.
filter_periods_from_env <- function(periods) {
  raw_value <- Sys.getenv("PME_PERIODS", unset = "")
  if (raw_value == "") {
    return(periods)
  }

  requested <- trimws(strsplit(raw_value, ",", fixed = TRUE)[[1]])
  requested <- requested[requested != ""]
  if(!length(requested)||anyDuplicated(requested))stop("PME_PERIODS must select distinct nonempty period ids.")
  missing <- requested[!(requested %in% periods$period_id)]
  if (length(missing) > 0L) {
    stop(
      "PME_PERIODS contains unknown period ids: ",
      paste(missing, collapse = ", "),
      ". Valid ids are: ",
      paste(periods$period_id, collapse = ", ")
    )
  }

  periods[match(requested, periods$period_id), , drop = FALSE]
}

#' Build the period/reference PME checkpoint filename
#' @param checkpoint_dir Checkpoint root.
#' @param period_id Target period key.
#' @param ref_index One-based reference position in the full ordered model list.
#' @param ref_name Exact model identifier; sanitised only for the filename.
#' @return Path under checkpoint_dir/period_id with index and model name. Does
#'   not create directories or prove compatibility of any existing payload.
checkpoint_path_for <- function(checkpoint_dir, period_id, ref_index, ref_name) {
  file.path(
    checkpoint_dir,
    period_id,
    paste0(period_id, "_ref_", sprintf("%02d", ref_index), "_", sanitize_name(ref_name), ".rds")
  )
}

#' Require a PME map's exact spatial array dimensions
#' @param map_obj Matrix/array to inspect.
#' @param expected_dim Expected dimension vector; no transposition or reshaping is guessed.
#' @param label Scientific map identity included in errors.
#' @return NULL invisibly on success; errors on wrong rank/shape. Does not itself
#'   certify finite values, successful optimisation, units or statistical identity.
validate_map_dim <- function(map_obj, expected_dim, label) {
  if (!is.matrix(map_obj) && !is.array(map_obj)) {
    stop(label, " is not matrix-like.")
  }
  if (!identical(as.integer(dim(map_obj)), as.integer(expected_dim))) {
    stop(
      label,
      " has dim [",
      paste(dim(map_obj), collapse = ", "),
      "] but expected [",
      paste(expected_dim, collapse = ", "),
      "]."
    )
  }
}

#' Construct a reproducible signature of the PME scientific implementation/settings
#' @param cfg Shared experiment configuration, optionally a compact test grid.
#' @return List of exact scientific parameters, canonical source fingerprints,
#'   R/platform and relevant installed software versions. Ignores ordinary R
#'   comments/formatting but not executable changes.
#'   Does not hash climate payloads; input/support provenance is checked separately.
pme_scientific_settings <- function(cfg=article1_config()) {
  fields<-c("variables","lon","lat","nbins3d","nbins2d","nbins1d",
    "pme_lambdas","pme_ldr_thresholds","season_start_md","season_end_md","seed",
    "strict_time_checks","season_selection_mode","time_decode_mode","requested_year_policy",
    "solver_repository","solver_commit","precipitation_preprocessing")
  files<-c("R/distributions/histogram_reference.R","R/distributions/histogram_fast.R",
    "R/distributions/pmf_pipeline.R","R/preprocessing/calendar.R","R/preprocessing/cf_calendar_validation.R",
    "R/preprocessing/netcdf_io.R","R/preprocessing/precipitation_qc.R","R/metrics/ldr.R","R/metrics/partial_hellinger.R",
    "R/graph_cut/solver.R","R/graph_cut/energy.R","R/experiments/perfect_model.R",
    "R/utils/runtime.R","R/utils/input_registry.R","experiments/perfect_model/run.R")
  list(schema="article1-pme-checkpoint-v2",parameters=cfg[fields],
    source_sha256=scientific_code_fingerprint(files),R_version=as.character(getRversion()),
    platform=R.version$platform,
    software_versions=vapply(c("Rcpp","gcoWrapR","RcppXPtrUtils","ncdf4","future","future.apply"),
      function(p)tryCatch(as.character(utils::packageVersion(p)),error=function(e)NA_character_),""))
}

#' Check scientific identity and successful numerical payload of a PME checkpoint
#' @param checkpoint_file Existing checkpoint path, or a missing candidate path.
#' @param period_row One-row table describing the requested target period.
#' @param calibration_years Exact requested integer calibration years.
#' @param period_role Scientific calibration/projection role identifier.
#' @param lambda_keys,tau_keys Ordered smoothness and ell_LDR nested-list keys.
#' @param expected_model_names Full ordered held-out/candidate model list.
#' @param range_sources Exact support path used by this computation.
#' @param partial_h_definition Recorded LDR evaluation-definition string.
#' @param expected_settings Signature from pme_scientific_settings(), updated with
#'   the actual numerical lambda/threshold grids when run-specific overrides apply.
#' @return TRUE only for matching, successful finite score/label payloads. Missing,
#'   old-signature, failed, malformed or different-input checkpoints return FALSE.
#'   Ordinary comments alone do not invalidate otherwise identical checkpoints.
checkpoint_is_compatible <- function(
  checkpoint_file,
  period_row,
  calibration_years,
  period_role,
  lambda_keys,
  tau_keys,
  expected_model_names,
  range_sources,
  partial_h_definition,
  expected_settings=pme_scientific_settings()
) {
  if (!file.exists(checkpoint_file)) {
    return(FALSE)
  }

  checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
  if (!is.list(checkpoint) || !is.list(checkpoint$meta)) {
    return(FALSE)
  }

  period_id <- period_row$period_id
  projection_years <- period_row$start_year:period_row$end_year
  meta <- checkpoint$meta
  if(!identical(meta$scientific_settings,expected_settings)||!identical(meta$checkpoint_state,"complete"))return(FALSE)
  expected<-experiment_input_provenance(range_sources)
  if(!identical(meta$input_provenance,expected))return(FALSE)

  if (!identical(meta$period_role, period_role)) {
    return(FALSE)
  }
  for(field in c("calibration_years","year_present","year_future"))
    if(!is.numeric(meta[[field]])||!length(meta[[field]])||any(!is.finite(meta[[field]]))||
       any(meta[[field]]!=floor(meta[[field]])))return(FALSE)
  if (!identical(as.numeric(meta$calibration_years), as.numeric(calibration_years))) {
    return(FALSE)
  }
  if (!identical(as.numeric(meta$year_present), as.numeric(calibration_years))) {
    return(FALSE)
  }
  if (!identical(as.numeric(meta$year_future), as.numeric(projection_years))) {
    return(FALSE)
  }
  if (!identical(as.character(meta$lambda_keys), as.character(lambda_keys))) {
    return(FALSE)
  }
  if (!identical(as.character(meta$partial_tau_keys), as.character(tau_keys))) {
    return(FALSE)
  }
  if (!identical(as.character(meta$model_names), as.character(expected_model_names))) {
    return(FALSE)
  }
  if (!identical(as.character(meta$range_sources), as.character(range_sources))) {
    return(FALSE)
  }
  if (!identical(meta$partial_h_definition, partial_h_definition)) {
    return(FALSE)
  }
  if (is.null(checkpoint$period_id) || !identical(checkpoint$period_id, period_id)) {
    return(FALSE)
  }

  # Metadata is not proof that a solver succeeded: reject failed/NA payloads.
  tryCatch({
  index<-meta$ref_index
  if(!is.numeric(index)||length(index)!=1L||!is.finite(index)||index!=floor(index)||
     index<1||index>length(expected_model_names))return(FALSE)
  ref<-expected_model_names[index]
  if(!identical(meta$ref_name,ref)||!identical(checkpoint$ref_name,ref)||
     !identical(basename(checkpoint_file),basename(checkpoint_path_for("",period_id,index,ref))))return(FALSE)
  expected_dim<-c(length(expected_settings$parameters$lon),length(expected_settings$parameters$lat))
  #' Validate one checkpoint score field against its configured grid
  #' Internal closure of checkpoint_is_compatible; not a standalone entry point.
  #' @param x Potential numeric longitude-by-latitude score array.
  #' @return TRUE when dimensions match the enclosing expected_dim and every score is finite and nonnegative.
  valid_score<-function(x)is.numeric(x)&&identical(dim(x),expected_dim)&&all(is.finite(x))&&all(x>=0)
  if(!valid_score(checkpoint$MMM_hdist_future_list[[period_id]][[ref]]))return(FALSE)
  for(tau in tau_keys)
    if(!valid_score(checkpoint$MMM_partial_hdist_future_list[[period_id]][[tau]][[ref]]))return(FALSE)
  for(lambda in lambda_keys) {
    fit<-checkpoint$GC_result_list[[period_id]][[lambda]][[ref]]
    if(!is.list(fit))return(FALSE)
    labels<-fit$label_attribution
    if(!is.numeric(labels)||!identical(dim(labels),expected_dim)||any(!is.finite(labels))||
       any(labels!=floor(labels))||any(labels<1|labels>=length(expected_model_names)))return(FALSE)
    if(!valid_score(checkpoint$GC_hdist_future_list[[period_id]][[lambda]][[ref]]))return(FALSE)
    for(tau in tau_keys)
      if(!valid_score(checkpoint$GC_partial_hdist_future_list[[period_id]][[tau]][[lambda]][[ref]]))return(FALSE)
  }
  TRUE
  },error=function(e)FALSE)
}

#' Initialise the nested containers for one PME period
#' @param lambda_keys Ordered GC smoothness keys.
#' @param tau_keys Ordered ell_LDR keys.
#' @return Empty score/solver/reference/checkpoint containers retaining the
#'   archived internal MMM/partial naming; no files or scientific results are created.
empty_period_result <- function(lambda_keys, tau_keys) {
  list(
    GC_result_list = setNames(lapply(lambda_keys, function(x) list()), lambda_keys),
    GC_hdist_future_list = setNames(lapply(lambda_keys, function(x) list()), lambda_keys),
    GC_partial_hdist_future_list = setNames(
      lapply(tau_keys, function(tau_key) setNames(lapply(lambda_keys, function(x) list()), lambda_keys)),
      tau_keys
    ),
    MMM_hdist_future_list = list(),
    MMM_partial_hdist_future_list = setNames(lapply(tau_keys, function(x) list()), tau_keys),
    checkpoint_files = character(0),
    completed_references = character(0)
  )
}

#' Assemble one period from scientifically compatible successful reference checkpoints
#' @param period_row One-row requested period table.
#' @param model_names Ordered references to look for.
#' @param checkpoint_dir Checkpoint root using the standard period/reference naming.
#' @param lambda_keys,tau_keys Ordered GC smoothness and ell_LDR keys.
#' @param expected_dim Exact spatial map dimensions.
#' @param calibration_years Integer requested calibration years.
#' @param period_role Recorded scientific period-role identifier.
#' @param expected_model_names Full model-list identity required in checkpoint metadata.
#' @param range_sources Requested support path, verified by hash as well as name.
#' @param partial_h_definition Recorded evaluation-definition string.
#' @param expected_settings Expected complete scientific signature.
#' @return Period result containers with only compatible completed references;
#'   absent/incompatible/failed references are skipped for later recomputation.
#'   Does not run graph cut or turn missing cases into zero scores.
read_period_from_checkpoints <- function(
  period_row,
  model_names,
  checkpoint_dir,
  lambda_keys,
  tau_keys,
  expected_dim,
  calibration_years,
  period_role,
  expected_model_names,
  range_sources,
  partial_h_definition,
  expected_settings=pme_scientific_settings()
) {
  period_id <- period_row$period_id
  out <- empty_period_result(lambda_keys, tau_keys)

  for (m in seq_along(model_names)) {
    ref_name <- model_names[m]
    checkpoint_file <- checkpoint_path_for(checkpoint_dir, period_id, m, ref_name)
    if (!checkpoint_is_compatible(
      checkpoint_file = checkpoint_file,
      period_row = period_row,
      calibration_years = calibration_years,
      period_role = period_role,
      lambda_keys = lambda_keys,
      tau_keys = tau_keys,
      expected_model_names = expected_model_names,
      range_sources = range_sources,
      partial_h_definition = partial_h_definition,
      expected_settings = expected_settings
    )) {
      next
    }

    checkpoint <- readRDS(checkpoint_file)

    validate_map_dim(
      checkpoint$MMM_hdist_future_list[[period_id]][[ref_name]],
      expected_dim,
      paste0("MMM future checkpoint map for ", period_id, "/", ref_name)
    )
    for (tau_key in tau_keys) {
      validate_map_dim(
        checkpoint$MMM_partial_hdist_future_list[[period_id]][[tau_key]][[ref_name]],
        expected_dim,
        paste0("MMM partial checkpoint map for ", period_id, "/", ref_name, " @ tau ", tau_key)
      )
    }
    for (lambda_key in lambda_keys) {
      validate_map_dim(
        checkpoint$GC_hdist_future_list[[period_id]][[lambda_key]][[ref_name]],
        expected_dim,
        paste0("GC future checkpoint map for ", period_id, "/", ref_name, " @ ", lambda_key)
      )
      for (tau_key in tau_keys) {
        validate_map_dim(
          checkpoint$GC_partial_hdist_future_list[[period_id]][[tau_key]][[lambda_key]][[ref_name]],
          expected_dim,
          paste0("GC partial checkpoint map for ", period_id, "/", ref_name, " @ ", lambda_key, ", tau ", tau_key)
        )
      }

      out$GC_result_list[[lambda_key]][[ref_name]] <-
        checkpoint$GC_result_list[[period_id]][[lambda_key]][[ref_name]]
      out$GC_hdist_future_list[[lambda_key]][[ref_name]] <-
        checkpoint$GC_hdist_future_list[[period_id]][[lambda_key]][[ref_name]]
      for (tau_key in tau_keys) {
        out$GC_partial_hdist_future_list[[tau_key]][[lambda_key]][[ref_name]] <-
          checkpoint$GC_partial_hdist_future_list[[period_id]][[tau_key]][[lambda_key]][[ref_name]]
      }
    }
    out$MMM_hdist_future_list[[ref_name]] <-
      checkpoint$MMM_hdist_future_list[[period_id]][[ref_name]]
    for (tau_key in tau_keys) {
      out$MMM_partial_hdist_future_list[[tau_key]][[ref_name]] <-
        checkpoint$MMM_partial_hdist_future_list[[period_id]][[tau_key]][[ref_name]]
    }
    out$checkpoint_files <- c(out$checkpoint_files, checkpoint_file)
    out$completed_references <- c(out$completed_references, ref_name)
  }

  out
}

#' Save a compact PME period result with explicit metadata
#' @param period_result Nested score/solver/checkpoint containers for one period.
#' @param period_row One-row target period descriptor.
#' @param period_output_dir Destination directory, created recursively if absent.
#' @param run_stamp Run identifier used in the filename; select a new run directory
#'   to avoid normal saveRDS overwrite behaviour.
#' @param meta_base Shared scientific/input metadata.
#' @return Written compact RDS path. Adds save time and completed-reference list;
#'   does not independently certify that all required references have completed.
write_period_result <- function(period_result, period_row, period_output_dir, run_stamp, meta_base) {
  dir.create(period_output_dir, recursive = TRUE, showWarnings = FALSE)
  period_id <- period_row$period_id
  period_file <- file.path(period_output_dir, paste0(run_stamp, "_", period_id, "_pme_period_compact.rds"))

  saveRDS(
    list(
      meta = c(
        meta_base,
        list(
          period = period_row,
          period_saved_at = as.character(Sys.time()),
          completed_references = period_result$completed_references,
          checkpoint_files = period_result$checkpoint_files
        )
      ),
      period_id = period_id,
      GC_result_list = period_result$GC_result_list,
      GC_hdist_future_list = period_result$GC_hdist_future_list,
      GC_partial_hdist_future_list = period_result$GC_partial_hdist_future_list,
      MMM_hdist_future_list = period_result$MMM_hdist_future_list,
      MMM_partial_hdist_future_list = period_result$MMM_partial_hdist_future_list
    ),
    file = period_file,
    compress = "xz"
  )

  period_file
}

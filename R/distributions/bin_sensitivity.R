# Bin-count stability diagnostics for the trivariate calibration PMF.
#
# This module intentionally stops before graph cut, future-period processing,
# and LDR. It streams one latitude block and one model at a time and retains
# only Hellinger scores long enough to calculate split-sample rank stability.

#' Load the histogram-sensitivity native routines for the current process/source
#' @param cpp_file Path to the Rcpp histogram/count-distance implementation.
#' @param cache_dir Rcpp build-cache directory; created when needed.
#' @param force_load Whether to reload even when source identity and functions match.
#' @return TRUE invisibly after native functions are available. Reuse requires
#'   matching canonical path/SHA-256 and process ID, not merely a global loaded flag.
#'   Multisession workers load their own DLL. Compilation requires Rcpp/toolchain.
load_pmf_sensitivity_cpp <- function(
  cpp_file = "src/bin_sensitivity.cpp",
  cache_dir = file.path(tempdir(), "MVGraphCutRcppCache"),
  force_load = FALSE
) {
  identity<-list(file=normalizePath(cpp_file,winslash="/",mustWork=TRUE),
    sha256=mvgc_sha256(cpp_file),process_id=Sys.getpid())
  if (!isTRUE(force_load) && isTRUE(getOption("mvgc.pmf_cpp_loaded", FALSE)) &&
      identical(getOption("mvgc.pmf_cpp_identity"),identity) &&
      exists("multik_hist_counts_cpp", mode = "function", inherits = TRUE) &&
      exists("multik_hellinger_cpp", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp is required for the global k-sensitivity experiment.")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(cpp_file, rebuild = FALSE, cacheDir = cache_dir, showOutput = FALSE)
  options(mvgc.pmf_cpp_loaded = TRUE)
  options(mvgc.pmf_cpp_identity = identity)
  invisible(TRUE)
}

#' Prepare matched odd/even calibration-date read plans across model variables
#' @param model_names Distinct prepared model identifiers, including the reference.
#' @param variables Ordered pr/tas/psl identifiers.
#' @param data_dir Prepared NetCDF hierarchy root unless a frozen registry is active.
#' @param years Requested calibration years; record year parity defines the split.
#' @param lon,lat Requested rectilinear coordinate vectors, without interpolation.
#' @param season_start_md,season_end_md Inclusive non-leap reference month/day boundaries.
#' @param strict_time_checks Whether to enforce legacy daily-mapping assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for unavailable/incomplete requested years.
#' @param verbose Whether to print per-model metadata progress.
#' @return Named model list with per-variable extraction plans, selected_year,
#'   odd/even logical masks. Reads time/coordinate metadata, not climate cubes;
#'   compares date identities across variables rather than their file-local offsets.
prepare_k_sensitivity_plans <- function(
  model_names,
  variables,
  data_dir,
  years,
  lon,
  lat,
  season_start_md = c(1, 1),
  season_end_md = c(12, 31),
  strict_time_checks = TRUE,
  season_selection_mode = "auto",
  time_decode_mode = "auto",
  requested_year_policy = "warn",
  verbose = TRUE
) {
  required_years <- sort(unique(as.integer(years)))
  plans <- setNames(vector("list", length(model_names)), model_names)

  for (model in model_names) {
    if (isTRUE(verbose)) cat(format(Sys.time(), "%H:%M:%S"), "- metadata", model, "\n")
    model_plans <- setNames(vector("list", length(variables)), variables)
    selected_indices <- NULL
    selected_years <- NULL

    for (variable in variables) {
      file_path <- find_matching_netcdf_file(data_dir, model, variable, required_years)
      extraction <- extract_fractional_variable_data(
        file_path = file_path,
        variable = variable,
        lon = lon,
        lat = lat,
        year_present = required_years,
        year_future = required_years,
        season_start_md = season_start_md,
        season_end_md = season_end_md,
        strict_time_checks = strict_time_checks,
        season_selection_mode = season_selection_mode,
        time_decode_mode = time_decode_mode,
        requested_year_policy = requested_year_policy,
        read_data = FALSE
      )

      if (is.null(selected_indices)) {
        selected_indices <- extraction$selected_present_dates
        selected_years <- extraction$selected_present_year
      } else if (!identical(selected_indices, extraction$selected_present_dates) ||
                 !identical(selected_years, extraction$selected_present_year)) {
        stop("Calibration selection differs across variables for model ", model, ".")
      }
      model_plans[[variable]] <- extraction$plan
    }

    available_years <- sort(unique(selected_years))
    if (!identical(available_years, required_years)) {
      warning(
        "Model ", model, " uses years ", paste(available_years, collapse = ","),
        " instead of the complete requested set.", call. = FALSE
      )
    }
    plans[[model]] <- list(
      variables = model_plans,
      selected_year = selected_years,
      odd = selected_years %% 2L == 1L,
      even = selected_years %% 2L == 0L
    )
  }
  plans
}

#' Read one model's calibration latitude block for resolution sensitivity
#' @param model_plan Model entry from prepare_k_sensitivity_plans().
#' @param variables Ordered pr/tas/psl identifiers.
#' @param lat_positions Positions in the requested latitude vector, not file-axis indices.
#' @return Named list of lon x block-lat x time cubes; pr uses the production
#'   shared nonnegative QC then log1p; tas/psl remain in prepared units. Opens/closes each
#'   variable file independently, including on read errors. Does not regrid data.
read_k_sensitivity_block <- function(model_plan, variables, lat_positions) {
  result <- setNames(vector("list", length(variables)), variables)
  for (variable in variables) {
    plan <- model_plan$variables[[variable]]
    nc_obj <- ncdf4::nc_open(plan$file_path)
    result[[variable]] <- tryCatch(
      read_variable_time_blocks(
        nc_obj = nc_obj,
        variable = variable,
        dimension_info = plan$dimension_info,
        lon_indices = plan$lon_indices,
        lat_indices = plan$lat_indices[lat_positions],
        time_indices = plan$selected_present
      ),
      finally = ncdf4::nc_close(nc_obj)
    )
  }
  result$pr <- precipitation_transform(result$pr)
  result
}

#' Construct one cell's odd/even trivariate PMFs at explicit resolutions
#' @param pr,tas,psl Aligned numeric calibration samples; pr is already log transformed.
#' @param range_3d Three x two common-support matrix in pr/tas/psl order.
#' @param split_masks List of logical sample masks, normally odd and even years.
#' @param k_values Positive integer per-variable resolutions to evaluate.
#' @return Named resolution list, each containing split PMF vectors of length k^3.
#'   Empty complete-case splits retain NULL slots rather than deleting/shifting
#'   list entries. Uses the same equal-width floating-edge convention as production.
cell_pmfs_for_k <- function(pr, tas, psl, range_3d, split_masks, k_values) {
  output <- setNames(vector("list", length(k_values)), as.character(k_values))
  for (k in k_values) {
    indices <- list(
      equal_width_bin_index(pr, range_3d[1, 1], range_3d[1, 2], k),
      equal_width_bin_index(tas, range_3d[2, 1], range_3d[2, 2], k),
      equal_width_bin_index(psl, range_3d[3, 1], range_3d[3, 2], k)
    )
    valid <- Reduce(function(x, y) x & !is.na(y), indices, init = rep(TRUE, length(pr)))
    split_pmfs <- vector("list", length(split_masks))
    for (s in seq_along(split_masks)) {
      keep <- valid & split_masks[[s]]
      if (!any(keep)) {
        split_pmfs[s] <- list(NULL)
      } else {
        counts <- tabulate_linear_states(
          cbind(indices[[1]][keep], indices[[2]][keep], indices[[3]][keep]),
          k
        )
        split_pmfs[[s]] <- counts / sum(counts)
      }
    }
    output[[as.character(k)]] <- split_pmfs
  }
  output
}

#' Compute the sensitivity reference helper's ordinary Hellinger distance
#' @param p,q Valid normalised PMF vectors on identical ordered bins, or NULL for
#'   an empty split. Unlike the strict public API this internal helper assumes validity.
#' @return Scalar ordinary Hellinger distance, or NA_real_ when either split is NULL.
#'   No subset renormalisation, baseline percentage or latitude weighting is applied.
#' @examples
#' stopifnot(hellinger_distance_pmf(c(1,0),c(0,1)) == 1)
hellinger_distance_pmf <- function(p, q) {
  if (is.null(p) || is.null(q)) return(NA_real_)
  sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)
}

#' Build the R sensitivity reference PMFs for a latitude block
#' @param raw Named pr/tas/psl cubes, lon x local-lat x time, with transformed pr.
#' @param model_plan Reference selected-year odd/even masks.
#' @param range_var Full requested-grid lon x lat x variable x min/max support array.
#' @param lat_positions Full-grid latitude positions corresponding to raw's local rows.
#' @param k_values Integer per-variable resolutions.
#' @return Named resolution list of lon x local-lat x two split list arrays.
#'   Every entry is a k^3 PMF or NULL. R verification/reference path; global
#'   production rank sensitivity uses compact native count arrays instead.
compute_reference_pmfs_block <- function(raw, model_plan, range_var, lat_positions, k_values) {
  nlon <- dim(raw$pr)[1]
  nlat_block <- length(lat_positions)
  output <- setNames(vector("list", length(k_values)), as.character(k_values))
  for (k in k_values) {
    output[[as.character(k)]] <- array(vector("list", nlon * nlat_block * 2L),
                                       dim = c(nlon, nlat_block, 2L))
  }
  splits <- list(model_plan$odd, model_plan$even)
  for (i in seq_len(nlon)) {
    for (j_local in seq_len(nlat_block)) {
      j <- lat_positions[j_local]
      range_3d <- rbind(range_var[i, j, 1, ], range_var[i, j, 2, ], range_var[i, j, 3, ])
      pmfs <- cell_pmfs_for_k(
        raw$pr[i, j_local, ], raw$tas[i, j_local, ], raw$psl[i, j_local, ],
        range_3d, splits, k_values
      )
      for (k in k_values) {
        output[[as.character(k)]][i, j_local, 1] <- list(pmfs[[as.character(k)]][[1]])
        output[[as.character(k)]][i, j_local, 2] <- list(pmfs[[as.character(k)]][[2]])
      }
    }
  }
  output
}

#' Evaluate one candidate's split Hellinger scores through the R reference path
#' @param raw Named pr/tas/psl calibration cubes with already transformed pr.
#' @param model_plan Candidate selected-year odd/even masks.
#' @param reference_pmfs Reference list arrays from compute_reference_pmfs_block().
#' @param range_var Common full-grid support array [lon,lat,3,2].
#' @param lat_positions Full-grid positions for this latitude block.
#' @param k_values Integer per-variable resolutions in desired output order.
#' @return Numeric array lon x local-lat x resolution x two splits. Empty sample
#'   comparisons are NA. This is a verification helper, not graph-cut optimisation.
compute_model_h_scores_block <- function(
  raw, model_plan, reference_pmfs, range_var, lat_positions, k_values
) {
  nlon <- dim(raw$pr)[1]
  nlat_block <- length(lat_positions)
  scores <- array(NA_real_, dim = c(nlon, nlat_block, length(k_values), 2L))
  splits <- list(model_plan$odd, model_plan$even)
  for (i in seq_len(nlon)) {
    for (j_local in seq_len(nlat_block)) {
      j <- lat_positions[j_local]
      range_3d <- rbind(range_var[i, j, 1, ], range_var[i, j, 2, ], range_var[i, j, 3, ])
      pmfs <- cell_pmfs_for_k(
        raw$pr[i, j_local, ], raw$tas[i, j_local, ], raw$psl[i, j_local, ],
        range_3d, splits, k_values
      )
      for (k_index in seq_along(k_values)) {
        key <- as.character(k_values[k_index])
        for (s in 1:2) {
          ref <- reference_pmfs[[key]][i, j_local, s][[1]]
          scores[i, j_local, k_index, s] <- hellinger_distance_pmf(pmfs[[key]][[s]], ref)
        }
      }
    }
  }
  scores
}

#' Compute odd/even candidate-ranking stability for one latitude block
#' @param lat_positions Positions of the block in the requested latitude vector.
#' @param plans Matched model read plans from prepare_k_sensitivity_plans().
#' @param model_names Distinct model identifiers including the reference.
#' @param reference_name Model treated as reference, excluded from ranked candidates.
#' @param variables Ordered pr/tas/psl identifiers.
#' @param range_var Common full-grid support array [lon,lat,3,2].
#' @param k_values Integer per-variable histogram resolutions.
#' @param verbose Whether to print block progress.
#' @return List lat_positions, rho, n_models_valid and elapsed_seconds. rho and
#'   counts are lon x local-lat x resolution arrays. Spearman correlation requires
#'   at least three finite candidate pairs and nonconstant ranks in both splits;
#'   otherwise rho is NA. Streams models, retaining only Hellinger scores, not
#'   all PMFs. Does not run future periods, graph cut, LDR or tune production k=8.
compute_k_sensitivity_block <- function(
  lat_positions, plans, model_names, reference_name, variables,
  range_var, k_values, verbose = FALSE
) {
  # Multisession workers receive serialized R wrappers but must load the cached
  # native DLL in their own process.
  load_pmf_sensitivity_cpp()
  started <- proc.time()[["elapsed"]]
  if (isTRUE(verbose)) {
    cat(format(Sys.time(), "%H:%M:%S"), "- latitude block",
        min(lat_positions), "to", max(lat_positions), "\n")
  }
  range_block <- range_var[, lat_positions, , , drop = FALSE]
  reference_raw <- read_k_sensitivity_block(plans[[reference_name]], variables, lat_positions)
  reference_counts <- multik_hist_counts_cpp(
    reference_raw$pr, reference_raw$tas, reference_raw$psl,
    range_block, plans[[reference_name]]$odd, k_values
  )
  rm(reference_raw)

  candidates <- setdiff(model_names, reference_name)
  nlon <- dim(range_var)[1]
  h_a <- array(NA_real_, dim = c(nlon, length(lat_positions), length(candidates), length(k_values)))
  h_b <- h_a
  for (m in seq_along(candidates)) {
    model <- candidates[m]
    raw <- read_k_sensitivity_block(plans[[model]], variables, lat_positions)
    model_counts <- multik_hist_counts_cpp(
      raw$pr, raw$tas, raw$psl, range_block, plans[[model]]$odd, k_values
    )
    scores <- multik_hellinger_cpp(model_counts, reference_counts)
    dim(scores) <- c(nlon, length(lat_positions), length(k_values), 2L)
    h_a[, , m, ] <- scores[, , , 1]
    h_b[, , m, ] <- scores[, , , 2]
    rm(raw, model_counts, scores)
    gc(verbose = FALSE)
  }

  rho <- array(NA_real_, dim = c(nlon, length(lat_positions), length(k_values)))
  n_models_valid <- array(0L, dim = dim(rho))
  for (i in seq_len(nlon)) {
    for (j_local in seq_along(lat_positions)) {
      for (k_index in seq_along(k_values)) {
        a <- h_a[i, j_local, , k_index]
        b <- h_b[i, j_local, , k_index]
        valid <- is.finite(a) & is.finite(b)
        n_models_valid[i, j_local, k_index] <- sum(valid)
        if (sum(valid) >= 3L && length(unique(a[valid])) > 1L && length(unique(b[valid])) > 1L) {
          rho[i, j_local, k_index] <- suppressWarnings(
            stats::cor(a[valid], b[valid], method = "spearman")
          )
        }
      }
    }
  }
  list(
    lat_positions = lat_positions,
    rho = rho,
    n_models_valid = n_models_valid,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
}

#' Validate a freshly computed or cached histogram-resolution block
#' @param value Block list with rho, n_models_valid, lat_positions and elapsed_seconds.
#' @param dimensions Expected longitude x local-latitude x resolution dimensions.
#' @param lat_positions Exact ordered full-grid latitude positions.
#' @param candidate_count Total distinct ranked candidates, excluding the reference.
#' @return TRUE invisibly or an error. Undefined correlations may be NA, but not
#'   NaN/Inf; finite correlations must lie in [-1,1] and have at least three valid
#'   candidates. Counts must be integer-valued and within the candidate population.
#'   This is structural/physical validation, not a cryptographic payload checksum.
validate_k_sensitivity_block <- function(value, dimensions, lat_positions, candidate_count) {
  if (!is.list(value) || !identical(value$lat_positions, lat_positions) ||
      !is.numeric(value$rho) || !identical(dim(value$rho), dimensions) ||
      any(is.nan(value$rho) | is.infinite(value$rho)) || any(abs(value$rho) > 1, na.rm = TRUE) ||
      !is.numeric(value$n_models_valid) || !identical(dim(value$n_models_valid), dimensions) ||
      any(!is.finite(value$n_models_valid)) || any(value$n_models_valid < 0 | value$n_models_valid > candidate_count) ||
      any(value$n_models_valid != floor(value$n_models_valid)) ||
      any(is.finite(value$rho) & value$n_models_valid < 3) ||
      !is.numeric(value$elapsed_seconds) || length(value$elapsed_seconds) != 1L ||
      !is.finite(value$elapsed_seconds) || value$elapsed_seconds < 0)
    stop("Invalid sensitivity block scores, candidate counts, grid positions or timing.")
  invisible(TRUE)
}

#' Run calibration-only split-sample model-rank stability across histogram sizes
#'
#' Uses odd/even calibration years and the same fixed local physical support at
#' each k. Does not choose a new production discretisation, run graph cut or use
#' future-period/reference-specific supports. Exact scientific input/code identity
#' namespaces checkpoints; ordinary R comments/formatting do not invalidate them.
#' Legacy directory identity uses file stats/plans, not full payload hashes; use
#' the prepared-run launcher for fresh checksums before scientific reruns.
#' @param model_names Distinct model identifiers: one reference and at least three candidates.
#' @param reference_name Single model identifier excluded from candidate ranks.
#' @param variables Exactly c("pr","tas","psl").
#' @param data_dir Prepared daily NetCDF root unless registry resolution is active.
#' @param years Integer calibration years containing odd and even years.
#' @param lon,lat Nonempty unique finite requested coordinates; latitudes in [-90,90].
#' @param range_var Finite strictly increasing bounds [lon,lat,3,2], with pr transformed.
#' @param k_values Integer bins per variable, sorted/deduplicated; default 4 through 12.
#' @param workers Positive integer worker count; restores the prior future plan on exit.
#' @param latitude_block_size Positive integer requested latitude positions per task.
#' @param strict_time_checks Whether to enforce legacy daily-time assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional full-year selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for unavailable/incomplete years.
#' @param verbose Whether to print metadata/block progress.
#' @param checkpoint_dir Optional cache root. Namespaces exact identities, atomically
#'   installs new blocks, and rejects malformed/incompatible existing records.
#' @return List summary (area/unweighted mean rho and valid coverage), rho and
#'   n_models_valid arrays [lon,lat,k], and metadata including input/code cache
#'   identity and per-block reuse/timings. Correlations average over valid cells;
#'   this measures rank stability, not GC relative improvement over D-MMM.
run_k_sensitivity <- function(
  model_names,
  reference_name,
  variables,
  data_dir,
  years,
  lon,
  lat,
  range_var,
  k_values = 4:12,
  workers = 1L,
  latitude_block_size = 2L,
  strict_time_checks = TRUE,
  season_selection_mode = "auto",
  time_decode_mode = "auto",
  requested_year_policy = "warn",
  verbose = TRUE,
  checkpoint_dir = NULL
) {
  if(!is.character(model_names)||length(model_names)<4L||anyNA(model_names)||
     any(!nzchar(model_names))||anyDuplicated(model_names)||
     length(reference_name)!=1L||!reference_name%in%model_names||
     !identical(variables,c("pr","tas","psl")))
    stop("Rank stability requires one named reference and at least three distinct candidates with pr/tas/psl.")
  for(grid in list(lon,lat))
    if(!is.numeric(grid)||!length(grid)||any(!is.finite(grid))||anyDuplicated(grid))
      stop("Expected nonempty unique finite coordinates.")
  if(any(abs(lat)>90))stop("Latitudes must lie in [-90,90].")
  if(!is.numeric(years)||!length(years)||any(!is.finite(years))||any(years!=floor(years)))
    stop("Calibration years must be finite integers.")
  years<-sort(unique(as.integer(years)))
  if(!any(years%%2L==0L)||!any(years%%2L==1L))stop("Odd and even calibration years are both required.")
  if(!is.numeric(range_var)||!identical(dim(range_var),c(length(lon),length(lat),3L,2L))||
     any(!is.finite(range_var))||any(range_var[,,,2]<=range_var[,,,1]))
    stop("Expected finite increasing support bounds [longitude,latitude,3,2].")
  if(!is.numeric(k_values)||!length(k_values)||any(!is.finite(k_values))||any(k_values!=trunc(k_values)))
    stop("k_values must contain finite integer resolutions.")
  k_values <- sort(unique(as.integer(k_values)))
  if (any(k_values < 2L)) stop("All k values must be at least 2.")
  for(value in list(workers,latitude_block_size))
    if(!is.numeric(value)||length(value)!=1L||!is.finite(value)||value<1||
       value!=floor(value)||value>.Machine$integer.max)
      stop("workers and latitude_block_size must be finite positive integers.")
  workers <- as.integer(workers)
  latitude_block_size <- as.integer(latitude_block_size)

  plans <- prepare_k_sensitivity_plans(
    model_names, variables, data_dir, years, lon, lat,
    strict_time_checks = strict_time_checks,
    season_selection_mode = season_selection_mode,
    time_decode_mode = time_decode_mode,
    requested_year_policy = requested_year_policy,
    verbose = verbose
  )
  lat_blocks <- split(seq_along(lat), ceiling(seq_along(lat) / latitude_block_size))

  # A count/min/max filename is not a scientific identity: different models,
  # interior years, coordinates or k sets can have exactly those same summaries.
  require_packages("digest")
  registry<-Sys.getenv("MVGC_INPUT_REGISTRY","")
  source_files<-c("R/distributions/bin_sensitivity.R","src/bin_sensitivity.cpp",
    "R/distributions/histogram_fast.R","R/preprocessing/calendar.R",
    "R/preprocessing/cf_calendar_validation.R","R/preprocessing/netcdf_io.R")
  input_files<-unique(unlist(lapply(plans,function(p)
    vapply(p$variables,`[[`,"","file_path")),use.names=FALSE))
  input_info<-file.info(input_files)[,c("size","mtime"),drop=FALSE]
  cache_configuration<-list(model_names=model_names,reference_name=reference_name,
    variables=variables,years=years,lon=lon,lat=lat,k_values=k_values,
    range_sha256=digest::digest(range_var,algo="sha256"),plans=plans,
    registry_sha256=if(nzchar(registry))mvgc_sha256(registry)else NA_character_,
    legacy_file_info=if(nzchar(registry))NULL else input_info,
    strict_time_checks=strict_time_checks,season_selection_mode=season_selection_mode,
    time_decode_mode=time_decode_mode,requested_year_policy=requested_year_policy,
    source_sha256=scientific_code_fingerprint(source_files),R_version=as.character(getRversion()))
  cache_identity<-digest::digest(cache_configuration,algo="sha256")

  checkpoint_key <- paste0(
    "lon", length(lon), "_lat", length(lat), "_models", length(model_names),
    "_years", min(years), "-", max(years),
    "_k", min(k_values), "-", max(k_values)
  )
  if (!is.null(checkpoint_dir)) {
    checkpoint_dir<-file.path(checkpoint_dir,paste0("identity_",substr(cache_identity,1,24)))
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  #' Compute or reuse one scientifically compatible latitude block
  #' Internal closure of run_k_sensitivity; not a standalone entry point.
  #' @param block Integer positions in the enclosing run latitude vector.
  #' @return Validated block list containing rho, valid-model counts and checkpoint identity. Uses enclosing prepared plans and output directory.
  compute_or_load_block <- function(block) {
    checkpoint_path <- if (is.null(checkpoint_dir)) NULL else file.path(
      checkpoint_dir,
      sprintf("%s_lat%04d-%04d.rds", checkpoint_key, min(block), max(block))
    )
    if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
      cached <- readRDS(checkpoint_path)
      expected_dim<-c(length(lon),length(block),length(k_values))
      if (!identical(cached$cache_identity,cache_identity)||!identical(cached$lat_positions,block)||
          !identical(dim(cached$rho),expected_dim)||!identical(dim(cached$n_models_valid),expected_dim))
        stop("Invalid or incompatible sensitivity checkpoint; choose a new checkpoint directory: ",checkpoint_path)
      validate_k_sensitivity_block(cached, expected_dim, block, length(model_names) - 1L)
      cached$checkpoint_reused<-TRUE
      return(cached)
    }
    value <- compute_k_sensitivity_block(
      block, plans, model_names, reference_name, variables,
      range_var, k_values, verbose
    )
    validate_k_sensitivity_block(value, c(length(lon), length(block), length(k_values)), block, length(model_names) - 1L)
    value$cache_identity<-cache_identity
    value$checkpoint_reused<-FALSE
    if (!is.null(checkpoint_path)) {
      temporary_path <- paste0(checkpoint_path, ".", Sys.getpid(), ".tmp")
      saveRDS(value, temporary_path, compress = FALSE)
      if (!file.rename(temporary_path, checkpoint_path)) {
        stop("Could not atomically install checkpoint: ", checkpoint_path)
      }
    }
    value
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (workers == 1L) {
    future::plan(future::sequential)
  } else if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  block_results <- future.apply::future_lapply(
    lat_blocks,
    compute_or_load_block,
    future.packages = c("ncdf4", "stats", "Rcpp"),
    future.seed = TRUE,
    future.scheduling = Inf
  )

  rho <- array(NA_real_, dim = c(length(lon), length(lat), length(k_values)),
               dimnames = list(lon = lon, lat = lat, k = k_values))
  n_models_valid <- array(0L, dim = dim(rho), dimnames = dimnames(rho))
  for (block_result in block_results) {
    idx <- block_result$lat_positions
    rho[, idx, ] <- block_result$rho
    n_models_valid[, idx, ] <- block_result$n_models_valid
  }

  weights_lat <- pmax(0, cos(lat * pi / 180))
  summary <- data.frame(
    k = k_values,
    rho_area_weighted = NA_real_,
    rho_unweighted = NA_real_,
    valid_cells = integer(length(k_values)),
    total_cells = length(lon) * length(lat),
    valid_cell_fraction = NA_real_,
    valid_area_fraction = NA_real_,
    stringsAsFactors = FALSE
  )
  for (k_index in seq_along(k_values)) {
    values <- rho[, , k_index]
    valid <- is.finite(values)
    weights <- matrix(rep(weights_lat, each = length(lon)), nrow = length(lon))
    summary$valid_cells[k_index] <- sum(valid)
    summary$valid_cell_fraction[k_index] <- mean(valid)
    summary$valid_area_fraction[k_index] <- sum(weights[valid]) / sum(weights)
    if (!any(valid)) next # Undefined rank stability is NA, not a normalised NaN.
    summary$rho_area_weighted[k_index] <- stats::weighted.mean(values[valid], weights[valid])
    summary$rho_unweighted[k_index] <- mean(values[valid])
  }

  list(
    summary = summary,
    rho = rho,
    n_models_valid = n_models_valid,
    metadata = list(
      reference_name = reference_name,
      candidate_models = setdiff(model_names, reference_name),
      variables = variables,
      calibration_years = years,
      split_a = "odd calibration years",
      split_b = "even calibration years",
      k_values = k_values,
      workers = workers,
      latitude_block_size = latitude_block_size,
      cache_identity = cache_identity,
      source_sha256 = cache_configuration$source_sha256,
      checkpoint_reused = vapply(block_results,`[[`,logical(1),"checkpoint_reused"),
      block_elapsed_seconds = vapply(block_results, `[[`, numeric(1), "elapsed_seconds")
    )
  )
}

#' Build the calibration rank-stability plot with truthful scope annotations
#' @param result Result list from run_k_sensitivity() with summary, rho array and
#'   metadata. Requires a finite k=8 summary point; other undefined scores may be NA.
#' @return ggplot without writing files. Uses actual candidate/cell counts and
#'   calibration years; archived global input retains its published annotations.
#'   Rejects invalid/duplicate resolutions, impossible correlations or inconsistent
#'   scope metadata before rendering. Does not change any numerical result.
build_k_sensitivity_plot <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to render the article-style sensitivity figure.")
  }
  plot_data <- result$summary
  if (!is.data.frame(plot_data) || !all(c("k", "rho_area_weighted", "total_cells") %in% names(plot_data)) ||
      !nrow(plot_data) || !is.numeric(plot_data$k) || any(!is.finite(plot_data$k)) ||
      any(plot_data$k < 2 | plot_data$k != floor(plot_data$k)) || anyDuplicated(plot_data$k) ||
      is.unsorted(plot_data$k) || !is.numeric(plot_data$rho_area_weighted) ||
      any(is.nan(plot_data$rho_area_weighted) | is.infinite(plot_data$rho_area_weighted)) ||
      any(abs(plot_data$rho_area_weighted) > 1 + 1e-12, na.rm = TRUE))
    stop("Expected ordered unique resolutions and valid rank-stability summary scores.")
  cells <- dim(result$rho)
  if (length(cells) != 3L || cells[3] != nrow(plot_data) ||
      !is.numeric(plot_data$total_cells) || anyNA(plot_data$total_cells) ||
      any(plot_data$total_cells != prod(cells[1:2])))
    stop("Summary cell counts do not match the recorded sensitivity grid.")
  candidates <- result$metadata$candidate_models; years <- result$metadata$calibration_years
  if (!is.character(candidates) || length(candidates) < 3L || anyNA(candidates) ||
      any(!nzchar(candidates)) || anyDuplicated(candidates) || !is.numeric(years) ||
      !length(years) || any(!is.finite(years)) || any(years != floor(years)) || anyDuplicated(years))
    stop("Missing or invalid candidate-model/calibration-year metadata.")
  production <- plot_data[plot_data$k == 8L, , drop = FALSE]
  if (nrow(production) != 1L || !is.finite(production$rho_area_weighted))
    stop("A finite production k=8 score is required to render this diagnostic; inspect valid-cell coverage.")
  y_lower <- floor((min(plot_data$rho_area_weighted, na.rm = TRUE) - 0.001) * 200) / 200
  y_upper <- ceiling((max(plot_data$rho_area_weighted, na.rm = TRUE) + 0.001) * 200) / 200
  year_label <- if (length(years) > 1L && identical(sort(as.numeric(years)), as.numeric(seq(min(years), max(years)))))
    paste0(min(years), "-", max(years)) else paste(sort(years), collapse = ", ")

  sensitivity_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = k, y = rho_area_weighted)
  ) +
    ggplot2::geom_vline(xintercept = 8, color = "grey55", linetype = "dotted", linewidth = .3) +
    ggplot2::geom_line(color = "#1f77b4", linewidth = .45) +
    ggplot2::geom_point(color = "#1f77b4", fill = "white", shape = 21, size = 2.1, stroke = .55) +
    ggplot2::geom_point(
      data = production,
      color = "#d62728", fill = "#d62728", shape = 21, size = 2.6, stroke = .5
    ) +
    ggplot2::scale_x_continuous(breaks = plot_data$k) +
    ggplot2::scale_y_continuous(
      breaks = seq(y_lower, y_upper, by = 0.005),
      minor_breaks = seq(y_lower, y_upper, by = 0.0025),
      labels = function(x) sprintf("%.3f", x)
    ) +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
    ggplot2::labs(
      x = "Bins per variable (k)",
      y = "Area-weighted Spearman correlation"
    ) +
    ggplot2::theme_minimal(base_size = 8.5) +
    ggplot2::theme(panel.grid.major=ggplot2::element_line(colour="#EBEBEB",linewidth=.35),
      panel.grid.minor=ggplot2::element_line(colour="#F2F2F2",linewidth=.2),
      axis.title=ggplot2::element_text(size=8.5),axis.text=ggplot2::element_text(size=7.5,colour="#202020"),
      plot.margin=ggplot2::margin(5,6,4,5))

  sensitivity_plot
}

#' Save the Article 1 histogram-resolution diagnostic and its numerical outputs
#' @param result Result list from run_k_sensitivity(), including a finite k=8 score.
#' @param output_dir Output directory to create if absent. Fixed diagnostic
#'   filenames may overwrite earlier renders; production controllers require a new
#'   run directory. Validation precedes all writes, including CSV/RDS output.
#' @return Paths to global CSV, cellwise RDS, PNG and PDF, invisibly. Writes a
#'   29 x 18 cm figure at 300 dpi, with actual scope annotations. Does not change
#'   production k, launch an experiment, or imply scientific approval of a run.
write_k_sensitivity_outputs <- function(result, output_dir) {
  sensitivity_plot <- build_k_sensitivity_plot(result)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(output_dir, "k_sensitivity_global.csv")
  rds_path <- file.path(output_dir, "k_sensitivity_cellwise.rds")
  png_path <- file.path(output_dir, "pmf_resolution_sensitivity.png")
  pdf_path <- file.path(output_dir, "pmf_resolution_sensitivity.pdf")
  utils::write.csv(result$summary, csv_path, row.names = FALSE)
  saveRDS(result, rds_path, compress = "xz")
  ggplot2::ggsave(png_path, plot = sensitivity_plot, width = 13.275, height = 8.3, units = "cm", dpi = 300, bg = "white")
  ggplot2::ggsave(pdf_path, plot = sensitivity_plot, width = 13.275, height = 8.3, units = "cm", dpi = 300, device = grDevices::cairo_pdf)

  invisible(list(csv = csv_path, rds = rds_path, png = png_path, pdf = pdf_path))
}

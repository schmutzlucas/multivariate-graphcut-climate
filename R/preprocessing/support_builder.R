# Full-year support construction shared by ERA5 and PME.
#' Read a positive-integer environment setting for support construction
#' @param var_name Environment-variable name.
#' @param default_value Default positive integer when unset or invalid.
#' @return Integer, warning and using the default on fractional/nonfinite/
#'   nonpositive/out-of-range settings. Does not alter the process environment.
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

#' Read an enumerated environment setting for support construction
#' @param var_name Environment-variable name.
#' @param default_value Value to return when unset or empty; caller supplies a valid choice.
#' @param allowed Character vector of exact accepted values; matching is case sensitive.
#' @return Default or selected string; unrecognised nonempty settings raise an error.
read_choice_env <- function(var_name, default_value, allowed) {
  raw_value <- Sys.getenv(var_name, unset = "")
  if (raw_value == "") {
    return(default_value)
  }
  if (!(raw_value %in% allowed)) {
    stop(var_name, " must be one of: ", paste(allowed, collapse = ", "))
  }
  raw_value
}

#' Read a finite positive numeric resource setting
#' @param var_name Environment-variable name, usually a support memory-budget setting.
#' @param default_value Finite positive numeric default supplied by the caller.
#' @return Numeric value; unset/empty settings use the default, invalid nonempty
#'   settings warn and use it. Does not modify the environment or allocate memory.
read_positive_numeric_env <- function(var_name, default_value) {
  raw_value <- Sys.getenv(var_name, unset = "")
  if (raw_value == "") {
    return(as.numeric(default_value))
  }

  parsed <- suppressWarnings(as.numeric(raw_value))
  if (!is.finite(parsed) || parsed <= 0) {
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
    return(as.numeric(default_value))
  }

  parsed
}

#' Estimate the payload size of one dense spatial-time cube
#' @param nlon,nlat,ntime Nonnegative dimensions of the requested cube.
#' @param bytes_per_value Storage bytes per numeric value, default eight for R doubles.
#' @return Numeric size in GiB (2^30 bytes), despite the historical gb suffix.
#'   Does not include copies, object overhead, decompression or NetCDF caches.
#' @examples
#' stopifnot(estimate_selected_cube_gb(1024,1024,128) == 1)
estimate_selected_cube_gb <- function(nlon, nlat, ntime, bytes_per_value = 8) {
  (as.numeric(nlon) * as.numeric(nlat) * as.numeric(ntime) * bytes_per_value) / (1024^3)
}

#' Apply the support builder's heuristic per-worker memory multiplier
#' @param variable Prepared variable name; pr requires transformation temporaries.
#' @param cube_gb Dense selected-cube size in GiB.
#' @return Estimated GiB, using multiplier 2.2 for pr and 1.2 for other variables.
#'   This is a scheduling heuristic, not a measured or guaranteed peak-memory bound.
estimate_task_memory_gb <- function(variable, cube_gb) {
  memory_multiplier <- if (identical(variable, "pr")) 2.2 else 1.2
  cube_gb * memory_multiplier
}

#' Limit support workers using the heuristic per-task memory estimate
#' @param variable Variable name passed to the task-memory estimator.
#' @param slab_gb Estimated raw cube/slab size per task in GiB; must be positive.
#' @param requested_workers Positive integer upper limit from the resource configuration.
#' @param memory_budget_gb Finite positive scheduling budget in GiB.
#' @return Integer worker count between one and requested_workers. The minimum
#'   of one does not guarantee a task fits when its estimate exceeds the budget.
derive_variable_worker_count <- function(variable, slab_gb, requested_workers, memory_budget_gb) {
  task_memory_gb <- estimate_task_memory_gb(variable, slab_gb)
  memory_limited_workers <- max(1L, floor(memory_budget_gb / task_memory_gb))
  as.integer(max(1L, min(requested_workers, memory_limited_workers)))
}

#' Split contiguous NetCDF blocks into bounded-size chunks
#' @param blocks Data frame with finite integer start >= 1 and count >= 0 columns.
#' @param max_count Finite positive integer maximum records per output chunk.
#' @return Data frame of start/count chunks, preserving input block order and
#'   coverage. Zero-count blocks are omitted; empty input returns a typed empty table.
#' @examples
#' split_blocks_to_max_count(data.frame(start=4L,count=8L),3L)
split_blocks_to_max_count <- function(blocks, max_count) {
  if (!is.numeric(max_count) || length(max_count) != 1L ||
      !is.finite(max_count) || max_count < 1 || max_count != floor(max_count) ||
      max_count > .Machine$integer.max) {
    stop("max_count must be one finite positive integer.")
  }
  if (!is.data.frame(blocks) || !all(c("start", "count") %in% names(blocks)) ||
      !is.numeric(blocks$start) || !is.numeric(blocks$count) ||
      any(!is.finite(blocks$start)) || any(!is.finite(blocks$count)) ||
      any(blocks$start < 1 | blocks$start != floor(blocks$start)) ||
      any(blocks$count < 0 | blocks$count != floor(blocks$count))) {
    stop("blocks must have finite integer start >= 1 and count >= 0 columns.")
  }
  out <- vector("list", 0)
  cursor <- 1L

  for (b in seq_len(nrow(blocks))) {
    remaining <- blocks$count[b]
    start <- blocks$start[b]

    while (remaining > 0L) {
      count <- min(remaining, max_count)
      out[[cursor]] <- data.frame(start = start, count = count)
      cursor <- cursor + 1L
      start <- start + count
      remaining <- remaining - count
    }
  }

  if (!length(out)) return(data.frame(start = integer(), count = integer()))
  do.call(rbind, out)
}

#' Read one contiguous time slab and requested rectilinear spatial positions
#' @param nc_obj Open NetCDF handle, owned and closed by the caller.
#' @param variable Three-dimensional variable name.
#' @param dimension_info Original axis names, positions and sizes from get_variable_dimension_info().
#' @param lon_indices,lat_indices Nonempty integer positions in the original file axes.
#' @param time_start One-based first record of a contiguous time interval.
#' @param time_count Positive integer number of time records.
#' @param report_qc Whether to return list(data, qc) with raw precipitation
#'   statistics alongside corrected data; FALSE retains the ordinary array API.
#' @return Lon x lat x time array, retaining singleton axes and requested spatial
#'   order, or list(data, qc) when requested. Finite negative pr is corrected
#'   immediately after reading; statistics describe only the selected cells.
#'   Reads a bounding spatial rectangle; no regridding or unit conversion.
read_variable_time_chunk <- function(
  nc_obj,
  variable,
  dimension_info,
  lon_indices,
  lat_indices,
  time_start,
  time_count,
  report_qc = FALSE
) {
  lon_start <- min(lon_indices)
  lon_count <- max(lon_indices) - lon_start + 1L
  lat_start <- min(lat_indices)
  lat_count <- max(lat_indices) - lat_start + 1L

  lon_local <- lon_indices - lon_start + 1L
  lat_local <- lat_indices - lat_start + 1L

  start <- rep.int(1L, length(dimension_info$variable_size))
  count <- dimension_info$variable_size

  start[dimension_info$lon_pos] <- lon_start
  count[dimension_info$lon_pos] <- lon_count

  start[dimension_info$lat_pos] <- lat_start
  count[dimension_info$lat_pos] <- lat_count

  start[dimension_info$time_pos] <- time_start
  count[dimension_info$time_pos] <- time_count

  raw_cube <- ncdf4::ncvar_get(
    nc = nc_obj,
    varid = variable,
    start = start,
    count = count,
    collapse_degen = FALSE
  )

  raw_selected <- aperm(
    raw_cube,
    perm = c(dimension_info$lon_pos, dimension_info$lat_pos, dimension_info$time_pos)
  )[lon_local, lat_local, , drop = FALSE]
  statistics <- if (isTRUE(report_qc) && identical(variable, "pr")) precipitation_qc_statistics(raw_selected) else NULL
  if (identical(variable, "pr")) raw_selected <- precipitation_qc(raw_selected)
  if (isTRUE(report_qc)) return(list(data = raw_selected, qc = statistics))
  raw_selected
}

#' Verify a historical retained-negative exception, not current ingestion QC
#' @param values Nonempty numeric precipitation array/vector in mm per day.
#' @param model Prepared model identifier.
#' @param sha256 Verified whole-input SHA-256 supplied by the calling scanner.
#' @param exception NULL for strict nonnegativity, or a reviewed list containing
#'   model, sha256, minimum and action equal to retain_unchanged. The minimum
#'   must lie strictly above -1 and no higher than zero.
#' @return Number of negative observations, invisibly. Never modifies values,
#'   clips, converts units or changes the production precipitation transform.
#'   Nonfinite values, mismatched exceptions and out-of-bound negatives stop.
#'   Acceptance is a numerical preparation qualification, not physical validation.
#'   Retained only to interpret historical receipts. Current production reads use
#'   precipitation_qc() instead; this function does not implement the QC rerun.
validate_support_precipitation <- function(values, model, sha256, exception = NULL) {
  if (!is.numeric(values) || !length(values) || any(!is.finite(values)))
    stop("Prepared precipitation must contain finite numeric observations.")
  minimum <- 0
  if (!is.null(exception)) {
    if (!is.list(exception) || !identical(exception$model, model) ||
        !identical(exception$sha256, sha256) ||
        !identical(exception$action, "retain_unchanged") ||
        !is.numeric(exception$minimum) || length(exception$minimum) != 1L ||
        !is.finite(exception$minimum) || exception$minimum <= -1 || exception$minimum > 0)
      stop("Invalid or source-mismatched precipitation exception.")
    minimum <- exception$minimum
  }
  if (any(values < minimum))
    stop("Precipitation below approved source-specific minimum; no clipping permitted.")
  invisible(sum(values < 0))
}

#' Reduce a climate cube to local histogram-support extrema
#' @param cube Numeric lon x lat x time array in prepared physical units.
#' @param variable Variable identifier. pr is transformed with log1p before extrema;
#'   tas and psl remain unchanged. Finite negative pr is corrected before log1p.
#' @return List with lon x lat min/max matrices, omitting missing samples;
#'   nonfinite extrema become NA. This is not the strict corrected-input approval
#'   gate: infinities are invalid; missing observations are omitted. Cells with
#'   no valid observations retain NA bounds rather than an invented support.
compute_cube_min_max <- function(cube, variable) {
  if (identical(variable, "pr")) {
    cube <- precipitation_transform(cube)
  }

  cube_dims <- dim(cube)
  dim(cube) <- c(cube_dims[1] * cube_dims[2], cube_dims[3])
  ranges <- matrixStats::rowRanges(cube, na.rm = TRUE, useNames = FALSE)

  min_values <- matrix(ranges[, 1], nrow = cube_dims[1], ncol = cube_dims[2])
  max_values <- matrix(ranges[, 2], nrow = cube_dims[1], ncol = cube_dims[2])
  min_values[!is.finite(min_values)] <- NA_real_
  max_values[!is.finite(max_values)] <- NA_real_

  list(min = min_values, max = max_values)
}

#' Make constant empirical supports valid without forcing minima to zero
#' @param lower,upper Numeric arrays of empirical lower/upper bounds with the
#'   same shape. Jointly missing bounds are retained, not fabricated.
#' @param variable Variable name; precipitation minima must be nonnegative.
#' @return List min, max, constant and expansion. The empirical minimum is
#'   unchanged. For finite constant cells only, the upper bound is increased by
#'   max(1e-12, 64 * machine epsilon * max(1, abs(min))). On the log-precipitation
#'   scale this is an explicit numerical envelope around a point mass; identical
#'   observations remain in the first bin. Missing/constant cells are reported
#'   by production aggregation. Inconsistent, infinite or reversed bounds fail.
regularize_empirical_support <- function(lower, upper, variable) {
  if (!is.numeric(lower) || !is.numeric(upper) ||
      !identical(dim(lower), dim(upper)) || length(lower) != length(upper) ||
      any(is.infinite(lower)) || any(is.infinite(upper)) ||
      !identical(is.na(lower), is.na(upper)) || any(lower > upper, na.rm = TRUE))
    stop("Invalid empirical histogram bounds.")
  if (identical(variable, "pr") && any(lower < 0, na.rm = TRUE))
    stop("QC precipitation support cannot have a negative empirical minimum.")
  constant <- !is.na(lower) & lower == upper
  expansion <- lower * 0
  expansion[is.na(lower)] <- NA_real_
  expansion[constant] <- pmax(1e-12, 64 * .Machine$double.eps * pmax(1, abs(lower[constant])))
  upper[constant] <- upper[constant] + expansion[constant]
  list(min = lower, max = upper, constant = constant, expansion = expansion)
}

#' Plan a full-year support selection on the experiment's requested grid
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param file_path Input path used for coverage hints and diagnostic messages.
#' @param variable Three-dimensional prepared variable name.
#' @param years Requested year vector; may contain disjoint experiment windows.
#' @param lon,lat Requested coordinate vectors, matched without interpolation.
#' @param strict_time_checks Whether to enforce legacy daily-mapping assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional. The selected
#'   seasonal boundaries are always January 1 through December 31.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for unavailable/incomplete years under
#'   the legacy mapping; this is not independent proof of complete calendar dates.
#' @return List with selected indices, per-year summary, calendar/selection mode,
#'   original dimension metadata and matched file-coordinate indices.
select_full_year_indices <- function(
  nc_obj,
  file_path,
  variable,
  years,
  lon,
  lat,
  strict_time_checks,
  season_selection_mode,
  time_decode_mode,
  requested_year_policy
) {
  time_dim_name <- find_time_dimension_name(nc_obj)
  calendar_raw <- ncdf4::ncatt_get(nc_obj, time_dim_name, "calendar")$value
  calendar <- normalize_calendar_name(calendar_raw)

  time_structure <- infer_time_structure_daily(
    nc_obj = nc_obj,
    file_path = file_path,
    time_dim_name = time_dim_name,
    calendar = calendar,
    required_years = years,
    strict_time_checks = strict_time_checks,
    time_decode_mode = time_decode_mode,
    requested_year_policy = requested_year_policy
  )

  dimension_info <- get_variable_dimension_info(
    nc_obj = nc_obj,
    variable = variable,
    time_dim_name = time_dim_name
  )

  coord_indices <- match_requested_coordinates(
    nc_obj = nc_obj,
    dimension_info = dimension_info,
    lon = lon,
    lat = lat
  )

  season_start_md <- c(1, 1)
  season_end_md <- c(12, 31)
  start_ref_doy <- month_day_to_reference_doy(season_start_md)
  end_ref_doy <- month_day_to_reference_doy(season_end_md)
  wraps_year_end <- FALSE

  selection_mode_used <- season_selection_mode
  exact_time_components <- NULL
  if (!identical(season_selection_mode, "fractional")) {
    time_values <- read_time_values(nc_obj, time_dim_name)
    exact_time_components <- decode_time_components_for_exact_selection(
      nc_obj = nc_obj,
      time_dim_name = time_dim_name,
      time_values = time_values,
      time_units = time_structure$time_units,
      calendar = calendar,
      time_decode_mode = time_decode_mode
    )

    if (isTRUE(exact_time_components$ok)) {
      selection_mode_used <- "exact_month_day"
    } else if (identical(season_selection_mode, "exact_month_day")) {
      stop(
        "Exact month/day selection failed for file ",
        basename(file_path),
        ". ",
        exact_time_components$error_message
      )
    } else {
      selection_mode_used <- "fractional"
      warning(
        "Exact month/day selection unavailable for ",
        basename(file_path),
        ". Falling back to fractional selection. Reason: ",
        exact_time_components$error_message,
        call. = FALSE
      )
    }
  }

  if (identical(selection_mode_used, "exact_month_day")) {
    selection <- build_exact_month_day_indices_for_years(
      year_table = time_structure$year_table,
      timestep_year = time_structure$timestep_year,
      timestep_slot = time_structure$timestep_slot,
      timestep_month = exact_time_components$month,
      timestep_day = exact_time_components$day,
      selected_years = years,
      season_start_md = season_start_md,
      season_end_md = season_end_md,
      requested_year_policy = requested_year_policy
    )
  } else {
    selection <- build_fractional_indices_for_years(
      year_table = time_structure$year_table,
      timestep_year = time_structure$timestep_year,
      timestep_slot = time_structure$timestep_slot,
      selected_years = years,
      start_ref_doy = start_ref_doy,
      end_ref_doy = end_ref_doy,
      wraps_year_end = wraps_year_end,
      requested_year_policy = requested_year_policy
    )
  }

  list(
    indices = selection$indices,
    year_summary = selection$year_summary,
    calendar = calendar,
    selection_mode_used = selection_mode_used,
    dimension_info = dimension_info,
    coord_indices = coord_indices
  )
}

#' Stream one model-variable's experiment-scoped support extrema
#' @param model Single prepared model identifier, resolved through registry or data_dir.
#' @param variable Prepared variable name; pr extrema use log1p scale.
#' @param data_dir Prepared input hierarchy root when no input registry is active.
#' @param years Exact requested year union, not automatically the file's entire coverage.
#' @param lon,lat Requested longitude/latitude vectors in desired output order.
#' @param strict_time_checks Whether to enforce legacy daily-mapping assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional full-year selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for unavailable/incomplete years.
#' @param time_chunk_size Maximum records per chunk in chunked mode.
#' @param read_mode Full, slab, chunked, or auto. Auto compares selected-cube size
#'   with max_full_read_gb, choosing full or slab.
#' @param max_full_read_gb Positive GiB threshold for automatic full-cube reading.
#' @param slab_time_count Maximum records per slab in slab mode.
#' @return List with lon x lat min/max matrices and input/calendar/read diagnostics.
#'   Opens/closes its own file and omits missing samples; this legacy helper alone
#'   does not certify corrected-input provenance, complete dates or physical units.
compute_range_for_model_variable <- function(
  model,
  variable,
  data_dir,
  years,
  lon,
  lat,
  strict_time_checks,
  season_selection_mode,
  time_decode_mode,
  requested_year_policy,
  time_chunk_size,
  read_mode,
  max_full_read_gb,
  slab_time_count
) {
  file_path <- find_matching_netcdf_file(
    data_dir = data_dir,
    model = model,
    variable = variable,
    required_years = years
  )

  nc_obj <- ncdf4::nc_open(file_path)
  on.exit(ncdf4::nc_close(nc_obj), add = TRUE)

  selection <- select_full_year_indices(
    nc_obj = nc_obj,
    file_path = file_path,
    variable = variable,
    years = years,
    lon = lon,
    lat = lat,
    strict_time_checks = strict_time_checks,
    season_selection_mode = season_selection_mode,
    time_decode_mode = time_decode_mode,
    requested_year_policy = requested_year_policy
  )

  selected_cube_gb <- estimate_selected_cube_gb(
    nlon = length(lon),
    nlat = length(lat),
    ntime = length(selection$indices)
  )

  effective_read_mode <- read_mode
  if (identical(read_mode, "auto")) {
    effective_read_mode <- if (selected_cube_gb <= max_full_read_gb) "full" else "slab"
  }

  if (identical(effective_read_mode, "full")) {
    cube <- read_variable_time_blocks(
      nc_obj = nc_obj,
      variable = variable,
      dimension_info = selection$dimension_info,
      lon_indices = selection$coord_indices$lon_indices,
      lat_indices = selection$coord_indices$lat_indices,
      time_indices = selection$indices
    )
    stats <- compute_cube_min_max(cube, variable)
    min_values <- stats$min
    max_values <- stats$max
    rm(cube, stats)
    gc(verbose = FALSE)
  } else {
    max_count <- if (identical(effective_read_mode, "slab")) slab_time_count else time_chunk_size
    blocks <- split_blocks_to_max_count(
      split_contiguous_blocks(selection$indices),
      max_count = max_count
    )

    min_values <- matrix(Inf, nrow = length(lon), ncol = length(lat))
    max_values <- matrix(-Inf, nrow = length(lon), ncol = length(lat))

    for (b in seq_len(nrow(blocks))) {
      cube <- read_variable_time_chunk(
        nc_obj = nc_obj,
        variable = variable,
        dimension_info = selection$dimension_info,
        lon_indices = selection$coord_indices$lon_indices,
        lat_indices = selection$coord_indices$lat_indices,
        time_start = blocks$start[b],
        time_count = blocks$count[b]
      )

      stats <- compute_cube_min_max(cube, variable)
      min_values <- ifelse(is.na(stats$min), min_values, pmin(min_values, stats$min, na.rm = TRUE))
      max_values <- ifelse(is.na(stats$max), max_values, pmax(max_values, stats$max, na.rm = TRUE))

      rm(cube, stats)
      gc(verbose = FALSE)
    }
  }

  min_values[!is.finite(min_values)] <- NA_real_
  max_values[!is.finite(max_values)] <- NA_real_

  list(
    min = min_values,
    max = max_values,
    diagnostics = list(
      model = model,
      variable = variable,
      file_path = file_path,
      calendar = selection$calendar,
      selection_mode_used = selection$selection_mode_used,
      read_mode = effective_read_mode,
      selected_cube_gb = selected_cube_gb,
      selected_timesteps = length(selection$indices),
      year_summary = selection$year_summary
    )
  )
}

#' Build a separately reviewed candidate support using exact experiment scope
#'
#' Requires explicit MVGC_REBUILD_SUPPORTS=1 and external climate data. Configured
#' models/year unions come from support_scope(): ERA5 includes its reference and
#' 22 candidates in its two windows; PME uses its 15 models and one shared union
#' of calibration plus all five targets, never leave-one-out-specific supports.
#' Resource/read/decode settings are PME_RANGE_* environment overrides; output is
#' MVGC_SUPPORT_OUTPUT_FILE or an outputs/supports candidate. Existing files are
#' never overwritten. Restores the prior future plan when exiting.
#' @param mode ERA5 or PME experiment profile.
#' @return No scientific return value; writes an RDS containing named variable
#'   support arrays [0..359 longitude,-90..90 latitude,min/max] and rich metadata,
#'   and prints progress. Does not establish identity with the accepted support.
#'   See DATA_PREPROCESSING.md for the corrected-input workflow.
run_support_builder <- function(mode = c("era5", "pme")) {
  mode <- match.arg(mode)
  cfg <- article1_config()
  require_climate_data(cfg)
  require_packages(c("ncdf4", "future", "future.apply", "matrixStats"))
  if (Sys.getenv("MVGC_REBUILD_SUPPORTS", "0") != "1") stop("Set MVGC_REBUILD_SUPPORTS=1 explicitly to build a new support. Exact accepted supports are external inputs; see docs/DATA_PREPROCESSING.md.")
  variables <- c("pr", "tas", "psl")
  lon <- 0:359
  lat <- -90:90
  scope <- support_scope(mode,cfg)
  years <- scope$years
  data_dir <- Sys.getenv(
    "PME_RANGE_DATA_DIR",
    unset = Sys.getenv("MVGC_DATA_ROOT", unset = "")
  )
  if (!nzchar(data_dir)) stop("Set PME_RANGE_DATA_DIR or MVGC_DATA_ROOT to the CMIP6 data directory.")
  model_list_file <- if (mode == "era5") cfg$era5_models else cfg$pme_models
  model_names <- scan(model_list_file, what = "", quiet = TRUE)
  output_file <- Sys.getenv(
    "MVGC_SUPPORT_OUTPUT_FILE",
    unset = file.path(mvgc_output_dir("supports"), paste0(mode, "_candidate.rds"))
  )

  strict_time_checks <- TRUE
  season_selection_mode <- Sys.getenv("PME_RANGE_SEASON_SELECTION_MODE", unset = "auto")
  time_decode_mode <- Sys.getenv("PME_RANGE_TIME_DECODE_MODE", unset = "auto")
  requested_year_policy <- read_choice_env(
    "PME_RANGE_REQUESTED_YEAR_POLICY",
    default_value = "warn",
    allowed = c("warn", "stop")
  )
  workers <- read_positive_int_env("PME_RANGE_WORKERS", default_value = 6L)
  read_mode <- read_choice_env(
    "PME_RANGE_READ_MODE",
    default_value = "slab",
    allowed = c("auto", "full", "slab", "chunked")
  )
  max_full_read_gb <- read_positive_numeric_env("PME_RANGE_MAX_FULL_READ_GB", default_value = 12)
  memory_budget_gb <- read_positive_numeric_env("PME_RANGE_MEMORY_BUDGET_GB", default_value = 128)
  time_chunk_size <- read_positive_int_env("PME_RANGE_TIME_CHUNK", default_value = 366L)
  slab_time_count <- read_positive_int_env("PME_RANGE_SLAB_TIME_COUNT", default_value = 8192L)
  output_compress <- read_choice_env(
    "PME_RANGE_OUTPUT_COMPRESS",
    default_value = "gzip",
    allowed = c("gzip", "xz", "none")
  )

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  if(file.exists(output_file))stop("Refusing to overwrite an existing support: ",output_file)
  cat("Computing ", mode, " ranges for exact requested years: ", paste(years,collapse=","), ".\n", sep = "")
  cat("Models: ", length(model_names), "; variables: ", paste(variables, collapse = ", "), "\n", sep = "")
  cat("Output: ", output_file, "\n", sep = "")
  cat(
    "Read mode: ",
    read_mode,
    "; workers: ",
    workers,
    "; max full-read cube per worker: ",
    signif(max_full_read_gb, 4),
    " GB; slab time count: ",
    slab_time_count,
    "; memory budget: ",
    signif(memory_budget_gb, 4),
    " GB\n",
    sep = ""
  )

  final_ranges <- setNames(vector("list", length(variables)), variables)
  diagnostics <- list()
  default_selected_cube_gb <- estimate_selected_cube_gb(
    nlon = length(lon),
    nlat = length(lat),
    ntime = length(years) * 366L
  )
  default_slab_gb <- estimate_selected_cube_gb(
    nlon = length(lon),
    nlat = length(lat),
    ntime = slab_time_count
  )
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  for (variable in variables) {
    cat("\nVariable: ", variable, "\n", sep = "")

    global_min <- matrix(Inf, nrow = length(lon), ncol = length(lat))
    global_max <- matrix(-Inf, nrow = length(lon), ncol = length(lat))
    workers_for_variable <- derive_variable_worker_count(
      variable = variable,
      slab_gb = if (identical(read_mode, "full")) default_selected_cube_gb else default_slab_gb,
      requested_workers = workers,
      memory_budget_gb = memory_budget_gb
    )
    if (.Platform$OS.type == "unix") {
      future::plan(future::multicore, workers = workers_for_variable)
    } else {
      future::plan(future::multisession, workers = workers_for_variable)
    }
    cat(
      "  -> processing ",
      length(model_names),
      " models with ",
      workers_for_variable,
      " workers (estimated task memory ",
      signif(
        estimate_task_memory_gb(
          variable,
          if (identical(read_mode, "full")) default_selected_cube_gb else default_slab_gb
        ),
        4
      ),
      " GB)\n",
      sep = ""
    )

    model_ranges <- future.apply::future_lapply(
      X = model_names,
      FUN = function(model) {
        compute_range_for_model_variable(
        model = model,
        variable = variable,
        data_dir = data_dir,
        years = years,
        lon = lon,
        lat = lat,
        strict_time_checks = strict_time_checks,
        season_selection_mode = season_selection_mode,
        time_decode_mode = time_decode_mode,
        requested_year_policy = requested_year_policy,
        time_chunk_size = time_chunk_size,
        read_mode = read_mode,
        max_full_read_gb = max_full_read_gb,
        slab_time_count = slab_time_count
        )
      },
      future.packages = c("ncdf4", "stats", "matrixStats")
    )
    names(model_ranges) <- model_names

    diagnostics[[variable]] <- list()
    for (model in model_names) {
      model_range <- model_ranges[[model]]
      global_min <- ifelse(is.na(model_range$min), global_min, pmin(global_min, model_range$min, na.rm = TRUE))
      global_max <- ifelse(is.na(model_range$max), global_max, pmax(global_max, model_range$max, na.rm = TRUE))
      diagnostics[[variable]][[model]] <- model_range$diagnostics
    }
    rm(model_ranges)
    gc(verbose = FALSE)

    global_min[!is.finite(global_min)] <- NA_real_
    global_max[!is.finite(global_max)] <- NA_real_
    if(anyNA(global_min)||anyNA(global_max))stop("No valid pooled support observations for ",variable,
      "; missing cells must be investigated, not assigned fabricated bounds.")
    adjusted <- regularize_empirical_support(global_min,global_max,variable)

    range_array <- array(NA_real_, dim = c(length(lon), length(lat), 2))
    range_array[, , 1] <- adjusted$min
    range_array[, , 2] <- adjusted$max
    diagnostics[[variable]]$constant_cells <- which(adjusted$constant,arr.ind=TRUE)
    final_ranges[[variable]] <- range_array
  }

  saveRDS(
    list(
      ranges = final_ranges,
      meta = list(
        script = paste0("experiments/",if(mode=="pme")"perfect_model"else"era5","/build_supports.R"),
        input_registry_sha256 = if(nzchar(cfg$input_registry))mvgc_sha256(cfg$input_registry)else NA_character_,
        created_at = as.character(Sys.time()),
        precipitation_preprocessing = cfg$precipitation_preprocessing,
        years = years,
        season_start_md = c(1, 1),
        season_end_md = c(12, 31),
        lon = lon,
        lat = lat,
        variables = variables,
        model_list_file = model_list_file,
        model_names = model_names,
        data_dir = data_dir,
        requested_year_policy = requested_year_policy,
        season_selection_mode = season_selection_mode,
        time_decode_mode = time_decode_mode,
        workers = workers,
        memory_budget_gb = memory_budget_gb,
        read_mode = read_mode,
        max_full_read_gb = max_full_read_gb,
        slab_time_count = slab_time_count,
        time_chunk_size = time_chunk_size,
        diagnostics = diagnostics
      )
    ),
    file = output_file,
    compress = if (identical(output_compress, "none")) FALSE else output_compress
  )

  cat("\nSaved PME full-year range file: ", output_file, "\n", sep = "")
}

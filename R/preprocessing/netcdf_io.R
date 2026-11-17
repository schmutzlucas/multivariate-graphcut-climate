#' Group unique file indices into contiguous disk-read blocks
#' @param indices Numeric vector of positive finite integer indices. Unsorted
#'   input and repetitions are allowed; invalid positions raise an error.
#' @return Data frame with integer start/count columns in increasing index order.
#'   Empty input returns a zero-row table. Caller order is not represented here.
#' @examples
#' split_contiguous_blocks(c(5L, 1L, 3L, 3L, 4L))
split_contiguous_blocks <- function(indices) {
  if (!is.numeric(indices) || any(!is.finite(indices)) ||
      any(indices < 1 | indices != floor(indices)) ||
      any(indices > .Machine$integer.max)) {
    stop("indices must be finite positive integers within R's integer range.")
  }
  if (length(indices) == 0L) {
    return(data.frame(start = integer(0), count = integer(0), stringsAsFactors = FALSE))
  }

  idx <- sort(unique(as.integer(indices)))
  cut_points <- c(1L, which(diff(idx) != 1L) + 1L)
  ends <- c(cut_points[-1L] - 1L, length(idx))

  data.frame(
    start = idx[cut_points],
    count = idx[ends] - idx[cut_points] + 1L,
    stringsAsFactors = FALSE
  )
}

#' Read selected NetCDF records in longitude-latitude-time order
#'
#' Reads a bounding spatial rectangle per contiguous unique-time block, permutes
#' original file axes, and subsets the requested spatial coordinates. It does not
#' regrid, convert units, log-transform precipitation, or infer missing values.
#' Finite negative precipitation is set to zero immediately after each disk read.
#' @param nc_obj Open NetCDF handle, owned and closed by the caller.
#' @param variable Name of the three-dimensional variable to read.
#' @param dimension_info Axis metadata from \code{get_variable_dimension_info()}.
#' @param lon_indices,lat_indices Nonempty positive integer positions in original
#'   file axes. Requested order and repetitions are retained.
#' @param time_indices Positive integer time positions, possibly unordered or
#'   repeated. Empty input is allowed; positions outside the file are rejected.
#' @return Numeric array longitude x latitude x requested time. Singleton axes
#'   remain present; time ordering/repetitions exactly follow the request.
read_variable_time_blocks <- function(
  nc_obj,
  variable,
  dimension_info,
  lon_indices,
  lat_indices,
  time_indices
) {
  blocks <- split_contiguous_blocks(time_indices)
  requested_time <- as.integer(time_indices)
  selected_time <- sort(unique(requested_time))
  for (axis in c("lon", "lat")) {
    index <- if (axis == "lon") lon_indices else lat_indices
    limit <- dimension_info$variable_size[dimension_info[[paste0(axis, "_pos")]]]
    if (!is.numeric(index) || !length(index) || any(!is.finite(index)) ||
        any(index < 1 | index != floor(index) | index > limit)) {
      stop(axis, "_indices must be nonempty integer positions within the file axis.")
    }
  }
  if (length(selected_time) && max(selected_time) >
      dimension_info$variable_size[dimension_info$time_pos]) {
    stop("time_indices extend beyond the file time axis.")
  }

  lon_start <- min(lon_indices)
  lon_count <- max(lon_indices) - lon_start + 1L
  lat_start <- min(lat_indices)
  lat_count <- max(lat_indices) - lat_start + 1L

  lon_local <- lon_indices - lon_start + 1L
  lat_local <- lat_indices - lat_start + 1L

  output <- array(
    NA_real_,
    dim = c(length(lon_indices), length(lat_indices), length(selected_time))
  )

  write_from <- 1L
  for (b in seq_len(nrow(blocks))) {
    start <- rep.int(1L, length(dimension_info$variable_size))
    count <- dimension_info$variable_size

    start[dimension_info$lon_pos] <- lon_start
    count[dimension_info$lon_pos] <- lon_count

    start[dimension_info$lat_pos] <- lat_start
    count[dimension_info$lat_pos] <- lat_count

    start[dimension_info$time_pos] <- blocks$start[b]
    count[dimension_info$time_pos] <- blocks$count[b]

    raw_cube <- ncdf4::ncvar_get(
      nc = nc_obj,
      varid = variable,
      start = start,
      count = count,
      collapse_degen = FALSE
    )
    if (identical(variable, "pr")) raw_cube <- precipitation_qc(raw_cube)

    cube_lon_lat_time <- aperm(
      raw_cube,
      perm = c(dimension_info$lon_pos, dimension_info$lat_pos, dimension_info$time_pos)
    )

    block_slice <- cube_lon_lat_time[lon_local, lat_local, , drop = FALSE]

    write_to <- write_from + blocks$count[b] - 1L
    output[, , write_from:write_to] <- block_slice
    write_from <- write_to + 1L
  }

  # Read each unique time once, then restore the caller's order and repetitions.
  # Normal experiment plans are already sorted/unique: avoid a second cube copy.
  if (!identical(requested_time, selected_time)) {
    return(output[, , match(requested_time, selected_time), drop = FALSE])
  }
  output
}

#' Tabulate the fractional seasonal mapping for distinct year lengths
#' @param year_lengths Nonempty vector of positive integer model/empirical year lengths.
#' @param start_ref_doy Inclusive start day on the 365-day reference calendar.
#' @param end_ref_doy Inclusive end day on that reference calendar.
#' @param wraps_year_end Whether selection includes both ends of each year.
#' @return Data frame with one row per sorted unique year length and columns
#'   year_length, start_mod, end_mod and wraps_window; NULL for empty input.
#' @examples
#' build_mapping_by_year_length(c(360,365,366),1,365,FALSE)
build_mapping_by_year_length <- function(year_lengths, start_ref_doy, end_ref_doy, wraps_year_end) {
  unique_lengths <- sort(unique(as.integer(year_lengths)))

  mapping <- lapply(unique_lengths, function(length_value) {
    mapped <- map_fractional_window_to_year_length(
      year_length = length_value,
      start_ref_doy = start_ref_doy,
      end_ref_doy = end_ref_doy,
      wraps_year_end = wraps_year_end
    )

    data.frame(
      year_length = length_value,
      start_mod = mapped$start_mod,
      end_mod = mapped$end_mod,
      wraps_window = wraps_year_end,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, mapping)
}

#' Plan or read calibration and evaluation records for one climate variable
#'
#' Opens and closes its own file. Auto season selection prefers decoded inclusive
#' month/day boundaries and warns before falling back to fractional year slots.
#' Precipitation is nonnegative-QC'd by the shared disk reader; its logarithmic
#' transformation belongs to the PMF layer. Other physical variables are unchanged.
#' Selection diagnostics distinguish filenames, observed record slots and decoded
#' date identities. Callers aligning variables must compare selected date keys,
#' not file-local indices, because file start years may differ.
#' @param file_path Path to one prepared rectilinear daily NetCDF file.
#' @param variable Name of its three-dimensional climate variable.
#' @param lon,lat Requested coordinate vectors in degrees, without interpolation.
#' @param year_present,year_future Nonempty vectors of calibration/evaluation years.
#' @param season_start_md,season_end_md Inclusive month/day pairs on a 365-day
#'   reference calendar; wrapped intervals include both ends of each selected year.
#' @param strict_time_checks Whether to enforce the legacy daily mapping assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional.
#' @param time_decode_mode Auto, helpers_cf, or empirical date-decoding backend.
#' @param requested_year_policy Warn or stop for unavailable/incomplete requested
#'   years under the legacy mapping. Warn is retained from production.
#' @param read_data Logical; FALSE inspects metadata and creates plans without
#'   reading climate-variable cubes (coordinate/time axes are still read).
#' @return List with optional data_present/data_future cubes, selected file
#'   indices and date keys, per-record year/slot metadata, diagnostics, year_table
#'   and a reusable block-reading plan. Data arrays are lon x lat x selected time.
#'   Known supported calendars use independent CF date keys; Julian/unknown
#'   fallback identities are explicitly labelled unverified year/slot strings.
extract_fractional_variable_data <- function(
  file_path,
  variable,
  lon,
  lat,
  year_present,
  year_future,
  season_start_md,
  season_end_md,
  strict_time_checks = TRUE,
  season_selection_mode = c("auto", "exact_month_day", "fractional"),
  time_decode_mode = c("auto", "helpers_cf", "empirical"),
  requested_year_policy = c("warn", "stop"),
  read_data = TRUE
) {
  extraction_started <- proc.time()[["elapsed"]]
  season_selection_mode <- match.arg(season_selection_mode)
  time_decode_mode <- match.arg(time_decode_mode)
  requested_year_policy <- match.arg(requested_year_policy)

  nc_obj <- ncdf4::nc_open(file_path)
  on.exit(ncdf4::nc_close(nc_obj), add = TRUE)

  time_dim_name <- find_time_dimension_name(nc_obj)
  calendar_raw <- ncdf4::ncatt_get(nc_obj, time_dim_name, "calendar")$value
  calendar_raw_string <- if (
    !is.character(calendar_raw) || length(calendar_raw) == 0L || is.na(calendar_raw[1])
  ) {
    "missing"
  } else {
    as.character(calendar_raw[1])
  }
  calendar <- normalize_calendar_name(calendar_raw)

  required_years <- sort(unique(c(year_present, year_future)))

  time_structure <- infer_time_structure_daily(
    nc_obj = nc_obj,
    file_path = file_path,
    time_dim_name = time_dim_name,
    calendar = calendar,
    required_years = required_years,
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

  start_ref_doy <- month_day_to_reference_doy(season_start_md)
  end_ref_doy <- month_day_to_reference_doy(season_end_md)
  wraps_year_end <- start_ref_doy > end_ref_doy

  selection_mode_used <- season_selection_mode
  selection_fallback_reason <- NA_character_
  exact_decode_source <- NA_character_
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
      exact_decode_source <- exact_time_components$source
    } else if (identical(season_selection_mode, "exact_month_day")) {
      stop(
        "Exact month/day seasonal selection failed for file ",
        basename(file_path),
        ". ",
        exact_time_components$error_message
      )
    } else {
      selection_mode_used <- "fractional"
      selection_fallback_reason <- exact_time_components$error_message
      warning(
        "Exact month/day selection unavailable for ",
        basename(file_path),
        ". Falling back to fractional selection. Reason: ",
        selection_fallback_reason,
        call. = FALSE
      )
    }
  }

  if (identical(selection_mode_used, "exact_month_day")) {
    present_selection <- build_exact_month_day_indices_for_years(
      year_table = time_structure$year_table,
      timestep_year = time_structure$timestep_year,
      timestep_slot = time_structure$timestep_slot,
      timestep_month = exact_time_components$month,
      timestep_day = exact_time_components$day,
      selected_years = year_present,
      season_start_md = season_start_md,
      season_end_md = season_end_md,
      requested_year_policy = requested_year_policy
    )

    future_selection <- build_exact_month_day_indices_for_years(
      year_table = time_structure$year_table,
      timestep_year = time_structure$timestep_year,
      timestep_slot = time_structure$timestep_slot,
      timestep_month = exact_time_components$month,
      timestep_day = exact_time_components$day,
      selected_years = year_future,
      season_start_md = season_start_md,
      season_end_md = season_end_md,
      requested_year_policy = requested_year_policy
    )
  } else {
    present_selection <- build_fractional_indices_for_years(
      year_table = time_structure$year_table,
      timestep_year = time_structure$timestep_year,
      timestep_slot = time_structure$timestep_slot,
      selected_years = year_present,
      start_ref_doy = start_ref_doy,
      end_ref_doy = end_ref_doy,
      wraps_year_end = wraps_year_end,
      requested_year_policy = requested_year_policy
    )

    future_selection <- build_fractional_indices_for_years(
      year_table = time_structure$year_table,
      timestep_year = time_structure$timestep_year,
      timestep_slot = time_structure$timestep_slot,
      selected_years = year_future,
      start_ref_doy = start_ref_doy,
      end_ref_doy = end_ref_doy,
      wraps_year_end = wraps_year_end,
      requested_year_policy = requested_year_policy
    )
  }

  metadata_seconds <- proc.time()[["elapsed"]] - extraction_started
  read_started <- proc.time()[["elapsed"]]
  data_present <- NULL
  data_future <- NULL
  if (isTRUE(read_data)) {
    data_present <- read_variable_time_blocks(
      nc_obj = nc_obj,
      variable = variable,
      dimension_info = dimension_info,
      lon_indices = coord_indices$lon_indices,
      lat_indices = coord_indices$lat_indices,
      time_indices = present_selection$indices
    )

    data_future <- read_variable_time_blocks(
      nc_obj = nc_obj,
      variable = variable,
      dimension_info = dimension_info,
      lon_indices = coord_indices$lon_indices,
      lat_indices = coord_indices$lat_indices,
      time_indices = future_selection$indices
    )
  }
  read_seconds <- proc.time()[["elapsed"]] - read_started

  mapping_by_length <- if (identical(selection_mode_used, "fractional")) {
    build_mapping_by_year_length(
      year_lengths = time_structure$year_table$year_length,
      start_ref_doy = start_ref_doy,
      end_ref_doy = end_ref_doy,
      wraps_year_end = wraps_year_end
    )
  } else {
    data.frame(
      year_length = sort(unique(as.integer(time_structure$year_table$year_length))),
      start_mod = NA_integer_,
      end_mod = NA_integer_,
      wraps_window = wraps_year_end,
      stringsAsFactors = FALSE
    )
  }

  year_length_lookup <- setNames(
    as.integer(time_structure$year_table$year_length),
    as.character(time_structure$year_table$year)
  )
  selected_present_year <- time_structure$timestep_year[present_selection$indices]
  selected_future_year <- time_structure$timestep_year[future_selection$indices]
  selected_present_slot <- time_structure$timestep_slot[present_selection$indices]
  selected_future_slot <- time_structure$timestep_slot[future_selection$indices]
  selected_present_year_length <- as.integer(year_length_lookup[as.character(selected_present_year)])
  selected_future_year_length <- as.integer(year_length_lookup[as.character(selected_future_year)])

  selected_present_month <- NULL
  selected_present_day <- NULL
  selected_future_month <- NULL
  selected_future_day <- NULL
  if (isTRUE(!is.null(exact_time_components) && exact_time_components$ok)) {
    selected_present_month <- exact_time_components$month[present_selection$indices]
    selected_present_day <- exact_time_components$day[present_selection$indices]
    selected_future_month <- exact_time_components$month[future_selection$indices]
    selected_future_day <- exact_time_components$day[future_selection$indices]
  }

  diagnostics <- list(
    variable = variable,
    file = file_path,
    time_dimension = time_dim_name,
    calendar_raw = calendar_raw_string,
    calendar = calendar,
    calendar_mode = time_structure$calendar_mode,
    year_length_mode = time_structure$year_length_mode,
    decode_backend_requested = time_structure$decode_backend_requested,
    decode_backend = time_structure$decode_backend,
    fallback_used = time_structure$fallback_used,
    fallback_reason = time_structure$fallback_reason,
    helper_result_class = time_structure$helper_result_class,
    requested_year_policy = time_structure$requested_year_policy,
    season_selection_mode_requested = season_selection_mode,
    season_selection_mode_used = selection_mode_used,
    season_selection_fallback_reason = selection_fallback_reason,
    exact_decode_source = exact_decode_source,
    season_start_md = season_start_md,
    season_end_md = season_end_md,
    wraps_year_end = wraps_year_end,
    declared_calendar_mismatch = time_structure$declared_calendar_mismatch,
    theoretical_year_lengths_unique = time_structure$theoretical_year_lengths_unique,
    observed_year_lengths_unique = time_structure$observed_year_lengths_unique,
    inferred_year_length = time_structure$inferred_year_length,
    time_units = time_structure$time_units,
    n_timesteps = time_structure$n_time,
    step_days_min = time_structure$step_days_min,
    step_days_median = time_structure$step_days_median,
    step_days_max = time_structure$step_days_max,
    max_offset_round_error = time_structure$max_offset_round_error,
    time_phase = time_structure$time_phase,
    filename_year_start = time_structure$filename_year_start,
    filename_year_end = time_structure$filename_year_end,
    filename_dates_source = time_structure$filename_dates_source,
    filename_start_md = time_structure$filename_start_md,
    filename_end_md = time_structure$filename_end_md,
    starts_on_jan1 = time_structure$starts_on_jan1,
    full_year_coverage = time_structure$full_year_coverage,
    daily_time_step = time_structure$daily_time_step,
    has_duplicates_or_reverse = time_structure$has_duplicates_or_reverse,
    expected_length_matches = time_structure$expected_length_matches,
    time_structure_matches_fast_method = time_structure$assumptions_ok,
    anomaly_years = time_structure$anomaly_years,
    anomaly_years_requested = time_structure$anomaly_years_requested,
    anomaly_years_outside = time_structure$anomaly_years_outside,
    inferred_year_lengths = sort(unique(time_structure$year_table$year_length)),
    year_table = time_structure$year_table,
    mapping_by_year_length = mapping_by_length,
    present_selection = present_selection$year_summary,
    future_selection = future_selection$year_summary,
    selected_timesteps_present = length(present_selection$indices),
    selected_timesteps_future = length(future_selection$indices),
    timing = c(
      metadata_and_time_axis = metadata_seconds,
      netcdf_variable_reads = read_seconds
    )
  )

  # Compare daily identities, not file-local offsets: corrected precipitation
  # can start in 1950 while temperature/pressure files still start in 1940.
  if(calendar %in% c("standard","proleptic_gregorian","360_day","365_day","366_day")) {
    decoded<-audit_cf_dates(read_time_values(nc_obj,time_dim_name),time_structure$time_units,calendar)
    date_keys<-sprintf("%04d-%02d-%02d",decoded$year,decoded$month,decoded$day)
    date_identity_source<-"independent_CF_calendar"
  } else {
    date_keys<-paste(time_structure$timestep_year,time_structure$timestep_slot,sep="/")
    date_identity_source<-"inferred_year_slot_not_verified_calendar_dates"
  }

  list(
    selected_present_dates = date_keys[present_selection$indices],
    selected_future_dates = date_keys[future_selection$indices],
    date_identity_source = date_identity_source,
    data_present = data_present,
    data_future = data_future,
    selected_present = present_selection$indices,
    selected_future = future_selection$indices,
    selected_present_year = selected_present_year,
    selected_future_year = selected_future_year,
    selected_present_slot = selected_present_slot,
    selected_future_slot = selected_future_slot,
    selected_present_year_length = selected_present_year_length,
    selected_future_year_length = selected_future_year_length,
    selected_present_month = selected_present_month,
    selected_present_day = selected_present_day,
    selected_future_month = selected_future_month,
    selected_future_day = selected_future_day,
    year_table = time_structure$year_table,
    diagnostics = diagnostics,
    plan = list(
      file_path = file_path,
      variable = variable,
      dimension_info = dimension_info,
      lon_indices = coord_indices$lon_indices,
      lat_indices = coord_indices$lat_indices,
      selected_present = present_selection$indices,
      selected_future = future_selection$indices,
      selected_present_year = selected_present_year,
      selected_future_year = selected_future_year
    )
  )
}

#' Execute a prepared extraction plan for a subset of latitude positions
#' @param plan List returned as plan by \code{extract_fractional_variable_data()}.
#' @param lat_positions Nonempty integer positions within the plan's requested
#'   latitude selection, not positions in the original file axis.
#' @return List with present/future lon x selected-latitude x selected-time cubes.
#'   Opens/closes the file for this block; retains original data units.
read_fractional_variable_plan_block <- function(plan, lat_positions) {
  if (!is.list(plan) || is.null(plan$file_path) || is.null(plan$dimension_info)) {
    stop("Invalid prepared variable extraction plan.")
  }
  if (!is.numeric(lat_positions) || length(lat_positions) == 0L ||
      any(!is.finite(lat_positions)) || any(lat_positions != floor(lat_positions)) ||
      any(lat_positions < 1L) ||
      any(lat_positions > length(plan$lat_indices))) {
    stop("lat_positions are outside the prepared coordinate selection.")
  }
  lat_positions <- as.integer(lat_positions)

  nc_obj <- ncdf4::nc_open(plan$file_path)
  on.exit(ncdf4::nc_close(nc_obj), add = TRUE)
  list(
    present = read_variable_time_blocks(
      nc_obj = nc_obj,
      variable = plan$variable,
      dimension_info = plan$dimension_info,
      lon_indices = plan$lon_indices,
      lat_indices = plan$lat_indices[lat_positions],
      time_indices = plan$selected_present
    ),
    future = read_variable_time_blocks(
      nc_obj = nc_obj,
      variable = plan$variable,
      dimension_info = plan$dimension_info,
      lon_indices = plan$lon_indices,
      lat_indices = plan$lat_indices[lat_positions],
      time_indices = plan$selected_future
    )
  )
}

#' Count complete multivariate samples with an explicit empty-sample result
#' @param data_matrix Observation x variable matrix; vectors become one-column matrices.
#' @param range_matrix Variable x two matrix of lower/upper common histogram bounds.
#' @param nbins Positive integer number of bins per variable.
#' @param compute_hist_nd_fun Histogram callback accepting data, range_var and nbins.
#' @return Flattened count vector of length nbins to the number of variables.
#'   Rows with NA/NaN are removed; an entirely empty sample returns all zero counts.
#'   Infinite values are not removed by complete.cases and follow the callback's
#'   boundary policy. Corrected production inputs are independently checked finite.
#' @examples
#' safe_hist_counts(matrix(c(0.2,NA,0.8),ncol=1),matrix(c(0,1),nrow=1),8,compute_histND)
safe_hist_counts <- function(data_matrix, range_matrix, nbins, compute_hist_nd_fun) {
  if (is.null(dim(data_matrix))) {
    data_matrix <- matrix(data_matrix, ncol = 1)
  }

  valid_rows <- stats::complete.cases(data_matrix)
  if (!any(valid_rows)) {
    return(rep.int(0, nbins^ncol(data_matrix)))
  }

  compute_hist_nd_fun(
    data = data_matrix[valid_rows, , drop = FALSE],
    range_var = range_matrix,
    nbins = nbins
  )
}

#' Average observed values while retaining an explicit all-missing result
#' @param x Numeric observations; NA and NaN are omitted, not infinite values.
#' @return Scalar arithmetic mean, or NA_real_ for empty/all-missing input.
#'   Precipitation callers back-transform log values before invoking this helper.
#' @examples
#' stopifnot(is.na(safe_mean_value(c(NA_real_,NA_real_))))
safe_mean_value <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

#' Normalise finite nonnegative histogram counts separately at every grid cell
#' @param arr Numeric longitude x latitude x bin count array; axes must be nonempty.
#' @return Array of identical shape with each nonempty cell summing to one.
#'   All-zero cells remain all zero under the legacy convention: they are not
#'   valid PMFs, and native graph-cut validation rejects them. Invalid input
#'   counts raise errors rather than being silently replaced with zeros.
#' @examples
#' normalize_histogram_array(array(c(1,3),c(1,1,2)))
normalize_histogram_array <- function(arr) {
  if(!is.numeric(arr)||length(dim(arr))!=3L||any(dim(arr)<1L)||
     any(!is.finite(arr))||any(arr<0))
    stop("Histogram counts must be a finite nonnegative longitude x latitude x bins array.")
  totals <- rowSums(arr, dims = 2)
  arr_norm <- sweep(arr, c(1, 2), totals, "/")
  arr_norm[!is.finite(arr_norm)] <- 0
  arr_norm
}
#' Print the time/season mapping diagnostics from one extraction
#' @param diag_obj Diagnostics list returned by \code{extract_fractional_variable_data()}.
#' @return No scientifically meaningful return value; prints backend/fallback,
#'   empirical completeness, per-year selection tables and a completeness preview.
print_fractional_selection_diagnostics <- function(diag_obj) {
  cat("\n=== Seasonal Selection Diagnostics ===\n")
  cat("File:       ", diag_obj$file, "\n", sep = "")
  cat("Variable:   ", diag_obj$variable, "\n", sep = "")
  cat("Cal raw:    ", diag_obj$calendar_raw, "\n", sep = "")
  cat("Calendar:   ", diag_obj$calendar, "\n", sep = "")
  cat("Cal mode:   ", diag_obj$calendar_mode, "\n", sep = "")
  cat("Decoder req:", diag_obj$decode_backend_requested, "\n", sep = "")
  cat("Decoder use:", diag_obj$decode_backend, "\n", sep = "")
  cat("Fallback:   ", diag_obj$fallback_used, "\n", sep = "")
  if (isTRUE(diag_obj$fallback_used)) {
    cat("Fallback rsn:", diag_obj$fallback_reason, "\n", sep = "")
  }
  if (!is.null(diag_obj$helper_result_class) && !is.na(diag_obj$helper_result_class)) {
    cat("Helper cls: ", diag_obj$helper_result_class, "\n", sep = "")
  }
  cat("Year policy:", diag_obj$requested_year_policy, "\n", sep = "")
  if (!is.null(diag_obj$season_selection_mode_requested)) {
    cat("Season req: ", diag_obj$season_selection_mode_requested, "\n", sep = "")
  }
  if (!is.null(diag_obj$season_selection_mode_used)) {
    cat("Season use: ", diag_obj$season_selection_mode_used, "\n", sep = "")
  }
  if (!is.null(diag_obj$season_selection_fallback_reason) && !is.na(diag_obj$season_selection_fallback_reason)) {
    cat("Season fbk: ", diag_obj$season_selection_fallback_reason, "\n", sep = "")
  }
  if (!is.null(diag_obj$exact_decode_source) && !is.na(diag_obj$exact_decode_source)) {
    cat("Exact src:  ", diag_obj$exact_decode_source, "\n", sep = "")
  }
  if (!is.null(diag_obj$filename_dates_source) && !is.na(diag_obj$filename_dates_source)) {
    cat("Date src:   ", diag_obj$filename_dates_source, "\n", sep = "")
  }
  if (!is.null(diag_obj$year_length_mode) && !is.na(diag_obj$year_length_mode)) {
    cat("Year mode:  ", diag_obj$year_length_mode, "\n", sep = "")
  }
  if (isTRUE(diag_obj$declared_calendar_mismatch)) {
    cat("Cal mismatch: TRUE\n")
    cat(
      "Observed year lengths:   ",
      paste(diag_obj$observed_year_lengths_unique, collapse = ", "),
      "\n",
      sep = ""
    )
    cat(
      "Theoretical year lengths:",
      paste(diag_obj$theoretical_year_lengths_unique, collapse = ", "),
      "\n",
      sep = ""
    )
  } else {
    cat("Cal mismatch: FALSE\n")
  }
  if (!is.null(diag_obj$inferred_year_length) && !is.na(diag_obj$inferred_year_length)) {
    cat("Cal infer:  ", diag_obj$inferred_year_length, "\n", sep = "")
  }
  cat("Time dim:   ", diag_obj$time_dimension, "\n", sep = "")
  cat("Timesteps:  ", diag_obj$n_timesteps, "\n", sep = "")
  cat("Step min:   ", diag_obj$step_days_min, "\n", sep = "")
  cat("Step med:   ", diag_obj$step_days_median, "\n", sep = "")
  cat("Step max:   ", diag_obj$step_days_max, "\n", sep = "")
  cat("Step round: ", diag_obj$max_offset_round_error, "\n", sep = "")
  if (!is.null(diag_obj$time_phase) && !is.na(diag_obj$time_phase)) {
    cat("Step phase: ", diag_obj$time_phase, "\n", sep = "")
  }
  cat("File years: ", diag_obj$filename_year_start, "-", diag_obj$filename_year_end, "\n", sep = "")
  cat(
    "Start MD:   ",
    sprintf("%02d-%02d", diag_obj$filename_start_md[1], diag_obj$filename_start_md[2]),
    "\n",
    sep = ""
  )
  cat(
    "End MD:     ",
    sprintf("%02d-%02d", diag_obj$filename_end_md[1], diag_obj$filename_end_md[2]),
    "\n",
    sep = ""
  )
  cat("Starts Jan1:", diag_obj$starts_on_jan1, "\n")
  cat("Daily step: ", diag_obj$daily_time_step, "\n", sep = "")
  cat("Fast OK:    ", diag_obj$time_structure_matches_fast_method, "\n", sep = "")
  cat("Anomaly yrs:", paste(diag_obj$anomaly_years, collapse = ", "), "\n")
  cat("Anom req:   ", paste(diag_obj$anomaly_years_requested, collapse = ", "), "\n")
  cat("Anom out:   ", paste(diag_obj$anomaly_years_outside, collapse = ", "), "\n")

  cat("\nMapped window by year length:\n")
  print(diag_obj$mapping_by_year_length, row.names = FALSE)

  cat("\nSelected timesteps per PRESENT year:\n")
  print(diag_obj$present_selection, row.names = FALSE)

  cat("\nSelected timesteps per FUTURE year:\n")
  print(diag_obj$future_selection, row.names = FALSE)

  cat("\nYear completeness preview:\n")
  print(
    utils::head(
      diag_obj$year_table[
        ,
        c(
          "year",
          "year_length",
          "theoretical_year_length",
          "expected_timesteps_for_full_year",
          "observed_timesteps",
          "observed_unique_doy",
          "is_full_year"
        )
      ],
      10L
    ),
    row.names = FALSE
  )
}

#' Inspect one point's seasonal extraction and its time diagnostics
#' @param nc_file Path to a prepared daily NetCDF file.
#' @param variable Variable name to inspect.
#' @param year_span Requested years, sorted/deduplicated before extraction.
#' @param season_start_md,season_end_md Inclusive non-leap reference month/day boundaries.
#' @param lon_target,lat_target Scalar coordinates; missing coordinates default
#'   independently to the first corresponding file coordinate.
#' @param strict_time_checks Whether to enforce daily-mapping assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for unavailable/incomplete years.
#' @param verbose Whether to print the diagnostic tables.
#' @return Extraction diagnostics list. This convenience inspector reads the
#'   selected climate values at one point; it is not metadata-only.
inspect_fractional_seasonal_selection <- function(
  nc_file,
  variable,
  year_span,
  season_start_md,
  season_end_md,
  lon_target = NULL,
  lat_target = NULL,
  strict_time_checks = TRUE,
  season_selection_mode = c("auto", "exact_month_day", "fractional"),
  time_decode_mode = c("auto", "helpers_cf", "empirical"),
  requested_year_policy = c("warn", "stop"),
  verbose = TRUE
) {
  season_selection_mode <- match.arg(season_selection_mode)
  time_decode_mode <- match.arg(time_decode_mode)
  requested_year_policy <- match.arg(requested_year_policy)
  year_span <- sort(unique(as.integer(year_span)))

  if (is.null(lon_target) || is.null(lat_target)) {
    nc_obj <- ncdf4::nc_open(nc_file)
    on.exit(ncdf4::nc_close(nc_obj), add = TRUE)

    time_dim_name <- find_time_dimension_name(nc_obj)
    dimension_info <- get_variable_dimension_info(
      nc_obj = nc_obj,
      variable = variable,
      time_dim_name = time_dim_name
    )

    lon_values <- read_dimension_values(nc_obj, dimension_info$lon_dim_name)
    lat_values <- read_dimension_values(nc_obj, dimension_info$lat_dim_name)

    if (is.null(lon_target)) {
      lon_target <- ifelse(lon_values[1] >= 180, lon_values[1] - 360, lon_values[1])
    }
    if (is.null(lat_target)) {
      lat_target <- lat_values[1]
    }
  }

  result <- extract_fractional_variable_data(
    file_path = nc_file,
    variable = variable,
    lon = lon_target,
    lat = lat_target,
    year_present = year_span,
    year_future = year_span,
    season_start_md = season_start_md,
    season_end_md = season_end_md,
    strict_time_checks = strict_time_checks,
    season_selection_mode = season_selection_mode,
    time_decode_mode = time_decode_mode,
    requested_year_policy = requested_year_policy
  )

  if (isTRUE(verbose)) {
    print_fractional_selection_diagnostics(result$diagnostics)
  }

  result$diagnostics
}

#' Extract one location's raw seasonal series with calendar metadata
#' @param nc_file Path to a prepared daily NetCDF file.
#' @param variable Variable name, in its prepared units without further conversion.
#' @param lon_target,lat_target Scalar geographic coordinates matched without interpolation.
#' @param year_span Requested years, sorted and deduplicated.
#' @param season_start_md,season_end_md Inclusive non-leap reference month/day boundaries.
#' @param strict_time_checks Whether to enforce legacy daily-mapping assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for unavailable/incomplete years.
#' @return Data frame with year, year_length, doy and value. Here doy is the
#'   active mapping's slot (possibly empirical rank), not necessarily true DOY.
#'   Decodable exact-selection dates additionally provide month and day columns.
extract_fractional_point_series <- function(
  nc_file,
  variable,
  lon_target,
  lat_target,
  year_span,
  season_start_md,
  season_end_md,
  strict_time_checks = TRUE,
  season_selection_mode = c("auto", "exact_month_day", "fractional"),
  time_decode_mode = c("auto", "helpers_cf", "empirical"),
  requested_year_policy = c("warn", "stop")
) {
  season_selection_mode <- match.arg(season_selection_mode)
  time_decode_mode <- match.arg(time_decode_mode)
  requested_year_policy <- match.arg(requested_year_policy)
  year_span <- sort(unique(as.integer(year_span)))

  extraction <- extract_fractional_variable_data(
    file_path = nc_file,
    variable = variable,
    lon = lon_target,
    lat = lat_target,
    year_present = year_span,
    year_future = year_span,
    season_start_md = season_start_md,
    season_end_md = season_end_md,
    strict_time_checks = strict_time_checks,
    season_selection_mode = season_selection_mode,
    time_decode_mode = time_decode_mode,
    requested_year_policy = requested_year_policy
  )

  values <- as.numeric(extraction$data_present[1, 1, ])
  if (length(values) != length(extraction$selected_present)) {
    stop(
      "Point-series extraction length mismatch: values=",
      length(values),
      ", selected_present=",
      length(extraction$selected_present),
      "."
    )
  }

  out <- data.frame(
    year = extraction$selected_present_year,
    year_length = extraction$selected_present_year_length,
    doy = extraction$selected_present_slot,
    value = values,
    stringsAsFactors = FALSE
  )

  if (!is.null(extraction$selected_present_month) && !is.null(extraction$selected_present_day)) {
    out$month <- extraction$selected_present_month
    out$day <- extraction$selected_present_day
  }

  out
}

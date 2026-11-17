#' Validate a seasonal boundary on the non-leap reference calendar
#'
#' Seasonal boundaries use a 365-day reference year, not a particular model
#' calendar. February 29 is therefore not an admissible boundary. Whole-year
#' extraction can nevertheless include leap days in the data.
#' @param md Numeric length-two vector containing integer month and day.
#' @param arg_name Argument label used in error messages.
#' @return Integer vector \code{c(month, day)}; invalid boundaries raise an error.
#' @examples
#' validate_month_day(c(12, 31), "season_end")
validate_month_day <- function(md, arg_name) {
  if (!is.numeric(md) || length(md) != 2L || any(!is.finite(md)) ||
      any(md != floor(md))) {
    stop(arg_name, " must contain two finite integers: c(month, day).")
  }

  month <- as.integer(md[1])
  day <- as.integer(md[2])

  if (month < 1L || month > 12L) {
    stop(arg_name, " month must be in [1, 12].")
  }
  if (day < 1L || day > 31L) {
    stop(arg_name, " day must be in [1, 31].")
  }

  ref_date <- as.Date(sprintf("2001-%02d-%02d", month, day), format = "%Y-%m-%d")
  if (is.na(ref_date)) {
    stop(
      arg_name,
      " = c(",
      month,
      ", ",
      day,
      ") is not a valid month/day in a 365-day reference year."
    )
  }

  c(month, day)
}

#' Convert a seasonal boundary to its one-based reference day
#' @param md Integer month/day pair valid in a 365-day year.
#' @return Integer day of year in 1 through 365, using year 2001.
#' @examples
#' stopifnot(month_day_to_reference_doy(c(3, 1)) == 60L)
month_day_to_reference_doy <- function(md) {
  md <- validate_month_day(md, "month_day")
  ref_date <- as.Date(sprintf("2001-%02d-%02d", md[1], md[2]), format = "%Y-%m-%d")
  as.integer(ref_date - as.Date("2001-01-01")) + 1L
}

#' Canonicalise the supported CF calendar names
#' @param calendar_raw Raw character calendar attribute; only its first value is used.
#' @return One of \code{standard}, \code{proleptic_gregorian}, \code{julian},
#'   \code{360_day}, \code{365_day}, \code{366_day}, or \code{unknown}.
#'   Missing, blank or non-character attributes default to \code{standard}.
#' @examples
#' stopifnot(normalize_calendar_name("noleap") == "365_day")
normalize_calendar_name <- function(calendar_raw) {
  if (!is.character(calendar_raw) || length(calendar_raw) == 0L || is.na(calendar_raw[1]) || trimws(calendar_raw[1]) == "") {
    return("standard")
  }

  cal <- tolower(trimws(calendar_raw[1]))

  if (cal %in% c("standard", "gregorian")) {
    return("standard")
  }
  if (cal %in% c("proleptic_gregorian", "proleptic gregorian")) {
    return("proleptic_gregorian")
  }
  if (cal %in% c("julian")) {
    return("julian")
  }
  if (cal %in% c("365_day", "365", "noleap", "no_leap")) {
    return("365_day")
  }
  if (cal %in% c("366_day", "366", "all_leap", "allleap")) {
    return("366_day")
  }
  if (cal %in% c("360_day", "360")) {
    return("360_day")
  }

  "unknown"
}

#' Estimate a fixed year length for an unsupported calendar
#'
#' This legacy fallback divides the number of records by the inclusive filename
#' year span and rounds to an integer. It cannot establish that individual dates
#' or complete years are actually present; use the independent CF audit for that.
#' @param parsed_dates List with integer \code{start_year} and \code{end_year}.
#' @param n_time Number of time records in the file.
#' @return Estimated year length, restricted to 200 through 500 days.
#' @examples
#' infer_unknown_year_length(list(start_year=2000L, end_year=2001L), 720L)
infer_unknown_year_length <- function(parsed_dates, n_time) {
  n_file_years <- parsed_dates$end_year - parsed_dates$start_year + 1L
  if (!is.finite(n_file_years) || n_file_years < 1L) {
    stop("Cannot infer unknown calendar: invalid filename year span.")
  }

  inferred <- as.integer(round(n_time / n_file_years))
  if (!is.finite(inferred) || inferred < 200L || inferred > 500L) {
    stop(
      "Cannot robustly infer unknown calendar year length from file length. ",
      "n_time=",
      n_time,
      ", n_file_years=",
      n_file_years,
      ", inferred_year_length=",
      inferred,
      "."
    )
  }

  inferred
}

#' Test the proleptic Gregorian leap-year rule
#' @param year Numeric vector of years.
#' @return Logical vector of equal length; missing years propagate.
#' @examples
#' stopifnot(identical(is_leap_gregorian(c(1900, 2000)), c(FALSE, TRUE)))
is_leap_gregorian <- function(year) {
  ((year %% 4L == 0L) & (year %% 100L != 0L)) | (year %% 400L == 0L)
}

#' Test the Julian leap-year rule
#' @param year Numeric vector of years.
#' @return Logical vector indicating divisibility by four; missing years propagate.
#' @examples
#' stopifnot(is_leap_julian(1900))
is_leap_julian <- function(year) {
  year %% 4L == 0L
}

#' Return theoretical year lengths under a calendar rule
#' @param year Numeric vector of years.
#' @param calendar One canonical name from \code{normalize_calendar_name()}.
#' @return Vector of year lengths. Fixed calendars return 360, 365 or 366;
#'   Julian uses divisibility by four; other names use Gregorian rules.
#' @examples
#' stopifnot(all(days_in_year_by_calendar(1999:2001, "360_day") == 360L))
days_in_year_by_calendar <- function(year, calendar) {
  if (calendar == "360_day") {
    return(rep.int(360L, length(year)))
  }
  if (calendar == "365_day") {
    return(rep.int(365L, length(year)))
  }
  if (calendar == "366_day") {
    return(rep.int(366L, length(year)))
  }
  if (calendar == "julian") {
    return(ifelse(is_leap_julian(year), 366L, 365L))
  }

  ifelse(is_leap_gregorian(year), 366L, 365L)
}

#' Return theoretical month lengths with explicit scalar broadcasting
#' @param year Finite integer years within R's integer range.
#' @param month Finite integer months in 1 through 12. Length must equal
#'   \code{year}, or one of the two arguments must have length one.
#' @param calendar Canonical CF calendar name; unrecognised names use Gregorian rules.
#' @return Numeric or integer month lengths in input order; an empty argument
#'   yields an empty integer vector. Partial recycling is rejected.
#' @examples
#' stopifnot(all(days_in_month_by_calendar(1999:2001, 2, "standard") == c(28,29,28)))
days_in_month_by_calendar <- function(year, month, calendar) {
  if (!is.numeric(year) || !is.numeric(month) || any(!is.finite(year)) ||
      any(!is.finite(month)) || any(year != floor(year)) ||
      any(month != floor(month)) || any(abs(year) > .Machine$integer.max) ||
      any(month < 1 | month > 12)) {
    stop("year and month must be finite integers, with month in [1, 12].")
  }
  if (!length(year) || !length(month)) return(integer())
  n <- max(length(year), length(month))
  if (!(length(year) %in% c(1L, n)) || !(length(month) %in% c(1L, n))) {
    stop("year and month must have equal lengths, or one must be a scalar.")
  }
  # Broadcast only scalars: R's partial recycling can silently pair wrong dates.
  year <- rep_len(as.integer(year), n)
  month <- as.integer(month)
  month <- rep_len(month, n)

  if (calendar == "360_day") {
    return(rep.int(30L, length(month)))
  }

  if (calendar == "365_day") {
    month_days <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
    return(month_days[month])
  }

  if (calendar == "366_day") {
    month_days <- c(31L, 29L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
    return(month_days[month])
  }

  if (calendar == "julian") {
    feb_days <- ifelse(is_leap_julian(year), 29L, 28L)
  } else {
    feb_days <- ifelse(is_leap_gregorian(year), 29L, 28L)
  }

  month_days <- cbind(
    31L,
    feb_days,
    31L,
    30L,
    31L,
    30L,
    31L,
    31L,
    30L,
    31L,
    30L,
    31L
  )
  month_days[cbind(seq_along(month), month)]
}

#' Read the last date-range token from a climate filename
#'
#' Performs lexical range validation, not full calendar validation: February 30
#' can be legitimate for a 360-day model. Filename coverage is a discovery hint,
#' never independent proof of the actual time-coordinate coverage.
#' @param file_path Single path whose basename may contain YYYYMMDD-YYYYMMDD.
#' @return List containing \code{ok} and start/end year, month and day integers.
#'   Missing or invalid tokens return \code{ok=FALSE} rather than throwing.
#' @examples
#' parse_filename_date_range("pr_MODEL_19500101-21001230.nc")
parse_filename_date_range <- function(file_path) {
  file_name <- basename(file_path)
  matches <- regmatches(file_name, gregexpr("[0-9]{8}-[0-9]{8}", file_name, perl = TRUE))[[1]]

  if (length(matches) == 0L || !nzchar(matches[1])) {
    return(list(
      ok = FALSE,
      start_year = NA_integer_,
      start_month = NA_integer_,
      start_day = NA_integer_,
      end_year = NA_integer_,
      end_month = NA_integer_,
      end_day = NA_integer_
    ))
  }

  token <- matches[length(matches)]
  parts <- strsplit(token, "-", fixed = TRUE)[[1]]
  start_token <- parts[1]
  end_token <- parts[2]

  parsed <- list(
    ok = TRUE,
    start_year = as.integer(substr(start_token, 1, 4)),
    start_month = as.integer(substr(start_token, 5, 6)),
    start_day = as.integer(substr(start_token, 7, 8)),
    end_year = as.integer(substr(end_token, 1, 4)),
    end_month = as.integer(substr(end_token, 5, 6)),
    end_day = as.integer(substr(end_token, 7, 8))
  )

  if (any(!is.finite(unlist(parsed[2:7])))) {
    parsed$ok <- FALSE
  }
  if (
    parsed$start_month < 1L || parsed$start_month > 12L ||
    parsed$end_month < 1L || parsed$end_month > 12L ||
    parsed$start_day < 1L || parsed$start_day > 31L ||
    parsed$end_day < 1L || parsed$end_day > 31L
  ) {
    parsed$ok <- FALSE
  }
  # Numeric YYYYMMDD ordering is valid for each supported calendar, without
  # imposing Gregorian month lengths on e.g. a legitimate 360-day February 30.
  if (as.numeric(start_token) > as.numeric(end_token)) {
    parsed$ok <- FALSE
  }

  parsed
}

#' Resolve one prepared NetCDF input for a model and variable
#'
#' An active \code{MVGC_INPUT_REGISTRY} takes precedence and never falls back to
#' directory discovery. Otherwise search \code{data_dir/model/variable}; with
#' multiple files choose the smallest filename year span covering the request.
#' A sole file is returned without proving its dates; extraction audits it later.
#' @param data_dir Root directory used only without an active input registry.
#' @param model Single model identifier, e.g. \code{ERA5}.
#' @param variable Prepared variable identifier, e.g. \code{pr}.
#' @param required_years Nonempty integer vector of requested years.
#' @return One path. Missing/ambiguous coverage raises an error; choosing among
#'   multiple covering files emits a warning.
find_matching_netcdf_file <- function(data_dir, model, variable, required_years) {
  if (nzchar(Sys.getenv("MVGC_INPUT_REGISTRY", "")))
    return(registry_input_file(model, variable, required_years))
  var_dir <- file.path(data_dir, model, variable)
  if (!dir.exists(var_dir)) {
    stop("Directory not found for model/variable: ", var_dir)
  }

  pattern <- glob2rx(paste0(variable, "_", model, "*.nc"))
  files <- sort(list.files(var_dir, pattern = pattern, full.names = TRUE))

  if (length(files) == 0L) {
    stop(
      "No NetCDF file found for model='", model,
      "', variable='", variable,
      "' in ", var_dir,
      ". Pattern used: ", pattern
    )
  }

  if (length(files) == 1L) {
    return(files[1])
  }

  year_min <- min(required_years)
  year_max <- max(required_years)

  file_table <- data.frame(
    file = files,
    parse_ok = FALSE,
    start_year = NA_integer_,
    end_year = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(files)) {
    parsed <- parse_filename_date_range(files[i])
    if (isTRUE(parsed$ok)) {
      file_table$parse_ok[i] <- TRUE
      file_table$start_year[i] <- parsed$start_year
      file_table$end_year[i] <- parsed$end_year
    }
  }

  can_cover <- with(
    file_table,
    parse_ok & (start_year <= year_min) & (end_year >= year_max)
  )

  if (any(can_cover)) {
    candidates <- file_table[can_cover, , drop = FALSE]
    span <- candidates$end_year - candidates$start_year
    ord <- order(span, candidates$start_year, candidates$file)
    chosen <- candidates$file[ord[1]]

    if (nrow(candidates) > 1L) {
      warning(
        "Multiple files match model='", model,
        "', variable='", variable,
        "' and cover requested years ", year_min, "-", year_max,
        ". Using file with smallest date span: ", basename(chosen),
        call. = FALSE
      )
    }

    return(chosen)
  }

  detail <- paste0("  - ", basename(files), collapse = "\n")
  stop(
    "Multiple files found for model='", model,
    "', variable='", variable,
    "', but none unambiguously covers requested years ",
    year_min, "-", year_max, ".\n",
    "Candidate files:\n", detail
  )
}

#' Find a supported time dimension in an open NetCDF file
#' @param nc_obj Open \code{ncdf4} handle, owned and closed by the caller.
#' @return \code{time}, or \code{valid_time} when \code{time} is absent;
#'   raises an error if neither dimension exists.
find_time_dimension_name <- function(nc_obj) {
  dim_names <- names(nc_obj$dim)

  if ("time" %in% dim_names) {
    return("time")
  }
  if ("valid_time" %in% dim_names) {
    return("valid_time")
  }

  stop("No time dimension named 'time' or 'valid_time' found in NetCDF file.")
}

#' Read raw numeric NetCDF time coordinates
#' @param nc_obj Open \code{ncdf4} handle, not closed by this function.
#' @param time_dim_name Name of the time coordinate/dimension.
#' @return Numeric vector, preferring a variable to cached dimension values.
#'   No decoding or time-zone conversion is performed.
read_time_values <- function(nc_obj, time_dim_name) {
  if (time_dim_name %in% names(nc_obj$var)) {
    return(as.numeric(ncdf4::ncvar_get(nc_obj, time_dim_name)))
  }

  if (time_dim_name %in% names(nc_obj$dim) && !is.null(nc_obj$dim[[time_dim_name]]$vals)) {
    return(as.numeric(nc_obj$dim[[time_dim_name]]$vals))
  }

  stop("Unable to read values for time dimension '", time_dim_name, "'.")
}

#' Read a rectilinear NetCDF coordinate vector
#' @param nc_obj Open \code{ncdf4} handle, not closed by this function.
#' @param dim_name Coordinate dimension name.
#' @return Numeric coordinates, preferring cached dimension values to a variable.
#'   Raises an error if neither representation is available.
read_dimension_values <- function(nc_obj, dim_name) {
  if (dim_name %in% names(nc_obj$dim) && !is.null(nc_obj$dim[[dim_name]]$vals)) {
    return(as.numeric(nc_obj$dim[[dim_name]]$vals))
  }

  if (dim_name %in% names(nc_obj$var)) {
    return(as.numeric(ncdf4::ncvar_get(nc_obj, dim_name)))
  }

  stop("Unable to read coordinate values for dimension '", dim_name, "'.")
}

#' Convert supported CF time units to a scale in days
#' @param time_units Character units attribute of the form days/hours/seconds since origin.
#' @return Scalar multiplier (1, 1/24, or 1/86400); unsupported units raise an error.
#' @examples
#' stopifnot(resolve_time_units_scale_days("hours since 1950-01-01") == 1/24)
resolve_time_units_scale_days <- function(time_units) {
  if (!is.character(time_units) || length(time_units) == 0L || is.na(time_units[1])) {
    stop("Missing NetCDF time units attribute.")
  }

  units_low <- tolower(trimws(time_units[1]))

  if (grepl("^days\\s+since\\s+", units_low)) {
    return(1)
  }
  if (grepl("^hours\\s+since\\s+", units_low)) {
    return(1 / 24)
  }
  if (grepl("^seconds\\s+since\\s+", units_low)) {
    return(1 / 86400)
  }

  stop("Unsupported time units '", time_units[1], "'. Expected days/hours/seconds since <origin>.")
}

#' Extract a Gregorian origin date from a CF units string
#'
#' Used for the explicitly estimated fallback range, not for exact non-Gregorian
#' date decoding. Time of day is discarded. Intended inputs have UTC origins;
#' arbitrary numeric time-zone offsets are not supported by this parser.
#' @param time_units One units string containing \code{since} and a date origin.
#' @return Scalar R \code{Date}; an unparseable origin raises an error.
#' @examples
#' parse_time_origin_date("days since 1950-01-01 00:00:00 UTC")
parse_time_origin_date <- function(time_units) {
  origin_str <- sub(".*since\\s+", "", time_units, ignore.case = TRUE)
  origin_str <- trimws(origin_str)
  origin_str <- sub(" UTC$", "", origin_str, ignore.case = TRUE)
  origin_str <- sub("Z$", "", origin_str)

  # Keep YYYY-MM-DD or YYYY-MM-DD HH:MM:SS forms.
  origin_date <- suppressWarnings(as.Date(origin_str))
  if (is.na(origin_date)) {
    origin_date <- suppressWarnings(as.Date(substr(origin_str, 1L, 10L)))
  }
  if (is.na(origin_date)) {
    stop("Could not parse origin date from time units: ", time_units)
  }

  origin_date
}

#' Extract a UTC date-time origin from a CF units string
#' @param time_units One units string containing \code{since}, an ISO date and
#'   optional time. UTC or Z suffixes are accepted; numeric zone offsets are not supported.
#' @return Scalar UTC \code{POSIXct}; missing or unparseable origins raise an error.
#' @examples
#' parse_time_origin_datetime("hours since 1950-01-01T12:00:00Z")
parse_time_origin_datetime <- function(time_units) {
  origin_str <- sub(".*since\\s+", "", time_units, ignore.case = TRUE)
  origin_str <- trimws(origin_str)
  origin_str <- sub(" UTC$", "", origin_str, ignore.case = TRUE)
  origin_str <- sub("Z$", "", origin_str)

  formats <- c(
    "%Y-%m-%d %H:%M:%OS",
    "%Y-%m-%dT%H:%M:%OS",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%dT%H:%M",
    "%Y-%m-%d"
  )

  for (fmt in formats) {
    parsed <- suppressWarnings(as.POSIXct(origin_str, format = fmt, tz = "UTC"))
    if (!is.na(parsed)) {
      return(parsed)
    }
  }

  stop("Could not parse origin datetime from time units: ", time_units)
}

#' Decode a Gregorian time coordinate to UTC dates
#'
#' This routine uses R's proleptic Gregorian/UTC arithmetic. It is appropriate
#' for the modern Article 1 Gregorian dates, not 360-day or no-leap calendars
#' or historical mixed-calendar reform transitions.
#' @param time_values Numeric offsets in the declared time units.
#' @param time_units One days/hours/seconds-since units string with UTC origin.
#' @return R \code{Date} vector; records are assigned to their UTC calendar day.
#' @examples
#' stopifnot(all(as.character(decode_time_to_dates_gregorian(c(0,24),
#'   "hours since 2000-01-01")) == c("2000-01-01","2000-01-02")))
decode_time_to_dates_gregorian <- function(time_values, time_units) {
  units_low <- tolower(trimws(time_units))
  origin_dt <- parse_time_origin_datetime(time_units)

  seconds <- if (grepl("^days\\s+since\\s+", units_low)) {
    time_values * 86400
  } else if (grepl("^hours\\s+since\\s+", units_low)) {
    time_values * 3600
  } else if (grepl("^seconds\\s+since\\s+", units_low)) {
    time_values
  } else {
    stop("Unsupported time units for Gregorian decoding: ", time_units)
  }

  dt <- origin_dt + seconds
  date_values <- as.Date(dt, tz = "UTC")
  if (anyNA(date_values)) {
    stop("Gregorian date decoding produced NA values from time axis.")
  }

  date_values
}

#' Identify the calendars eligible for the modern Gregorian decoder
#' @param calendar Character vector of canonical calendar names.
#' @return Logical vector, TRUE only for standard or proleptic_gregorian.
#' @examples
#' is_gregorian_like_calendar(c("standard", "360_day"))
is_gregorian_like_calendar <- function(calendar) {
  calendar %in% c("standard", "proleptic_gregorian")
}

#' Return the most frequent finite integer in a vector
#' @param x Numeric values, converted to integers before counting.
#' @return Scalar integer mode, smallest value on a frequency tie, or NA for
#'   no finite observations. Used for legacy empirical year-length diagnostics.
#' @examples
#' stopifnot(mode_integer(c(365L,366L)) == 365L)
mode_integer <- function(x) {
  x <- as.integer(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_integer_)
  }

  tab <- sort(table(x), decreasing = TRUE)
  as.integer(names(tab)[1])
}

#' Adapt optional ncdf4.helpers time decoding to R Date
#'
#' Tries three helper calling conventions and several Date coercions. Conversion
#' can fail for dates unrepresentable by R's Gregorian Date (e.g. February 30).
#' The caller chooses whether failure is fatal or permits a documented fallback.
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param time_dim_name Time dimension name supplied to the optional helper.
#' @return List with Date vector \code{date_values}, \code{helper_method} and
#'   \code{helper_result_class}; raises an error if all attempts fail.
decode_time_with_ncdf4_helpers <- function(nc_obj, time_dim_name) {
  if (!requireNamespace("ncdf4.helpers", quietly = TRUE)) {
    stop("Package 'ncdf4.helpers' is required for helpers_cf decode mode.")
  }

  attempts <- list(
    list(fun = ncdf4.helpers::nc.get.time.series, args = list(nc_obj, time.dim.name = time_dim_name)),
    list(fun = ncdf4.helpers::nc.get.time.series, args = list(nc_obj)),
    list(fun = ncdf4.helpers::nc.get.time.series, args = list(nc_obj, v = time_dim_name))
  )

  last_error <- NULL
  for (i in seq_along(attempts)) {
    candidate <- attempts[[i]]
    result <- tryCatch(
      do.call(candidate$fun, candidate$args),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      last_error <- result$message
      next
    }

    if (length(result) == 0L) {
      next
    }

    date_values <- tryCatch(
      suppressWarnings(as.Date(result)),
      error = function(e) rep(as.Date(NA), length(result))
    )
    if (anyNA(date_values)) {
      result_chr <- tryCatch(
        as.character(result),
        error = function(e) rep(NA_character_, length(result))
      )
      date_values <- suppressWarnings(as.Date(result_chr))
    }
    if (anyNA(date_values)) {
      result_posix <- tryCatch(
        as.POSIXct(result, tz = "UTC"),
        error = function(e) rep(as.POSIXct(NA, tz = "UTC"), length(result))
      )
      date_values <- suppressWarnings(as.Date(result_posix))
    }
    if (anyNA(date_values)) {
      result_fmt <- tryCatch(
        format(result, "%Y-%m-%d"),
        error = function(e) rep(NA_character_, length(result))
      )
      date_values <- suppressWarnings(as.Date(result_fmt))
    }
    if (anyNA(date_values)) {
      last_error <- paste0(
        "ncdf4.helpers returned an object that could not be coerced to Date. class=",
        paste(class(result), collapse = "/")
      )
      next
    }

    return(list(
      date_values = date_values,
      helper_method = paste0("nc.get.time.series attempt ", i),
      helper_result_class = paste(class(result), collapse = "/")
    ))
  }

  stop(
    "ncdf4.helpers failed to decode the time axis for dimension '",
    time_dim_name,
    "'. Last error: ",
    if (is.null(last_error)) "unknown" else last_error
  )
}

#' Derive fallback endpoints when no filename date range can be read
#'
#' Prefers the optional CF helper, then the modern Gregorian decoder. Last resort:
#' add rounded min/max offsets to a Gregorian origin and explicitly mark the range
#' estimated. The latter is not exact CF decoding for non-Gregorian calendars.
#' @param file_path Path used in diagnostic errors.
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param time_dim_name Time dimension name.
#' @param time_values Nonempty numeric time offsets.
#' @param time_units Character CF units attribute.
#' @param scale_days Multiplier converting offsets to days.
#' @param calendar Canonical calendar name.
#' @param time_decode_mode Requested auto, helpers_cf, or empirical backend.
#' @return List with endpoint fields in \code{parsed}, a \code{source} label,
#'   and optional caution \code{note}. Errors if no origin can be recovered.
derive_date_range_from_time_axis <- function(
  file_path,
  nc_obj,
  time_dim_name,
  time_values,
  time_units,
  scale_days,
  calendar,
  time_decode_mode
) {
  parsed <- list(
    ok = FALSE,
    start_year = NA_integer_,
    start_month = NA_integer_,
    start_day = NA_integer_,
    end_year = NA_integer_,
    end_month = NA_integer_,
    end_day = NA_integer_
  )
  source <- NA_character_
  note <- NULL

  helper_try <- NULL
  if (!identical(time_decode_mode, "empirical")) {
    helper_try <- tryCatch(
      decode_time_with_ncdf4_helpers(nc_obj = nc_obj, time_dim_name = time_dim_name),
      error = function(e) e
    )
  }

  if (!is.null(helper_try) && !inherits(helper_try, "error")) {
    dates <- helper_try$date_values
    parsed <- list(
      ok = TRUE,
      start_year = as.integer(format(dates[1L], "%Y")),
      start_month = as.integer(format(dates[1L], "%m")),
      start_day = as.integer(format(dates[1L], "%d")),
      end_year = as.integer(format(dates[length(dates)], "%Y")),
      end_month = as.integer(format(dates[length(dates)], "%m")),
      end_day = as.integer(format(dates[length(dates)], "%d"))
    )
    source <- "time_axis_helpers_cf"
    return(list(parsed = parsed, source = source, note = note))
  }

  if (is_gregorian_like_calendar(calendar)) {
    dates <- decode_time_to_dates_gregorian(time_values, time_units[1L])
    parsed <- list(
      ok = TRUE,
      start_year = as.integer(format(dates[1L], "%Y")),
      start_month = as.integer(format(dates[1L], "%m")),
      start_day = as.integer(format(dates[1L], "%d")),
      end_year = as.integer(format(dates[length(dates)], "%Y")),
      end_month = as.integer(format(dates[length(dates)], "%m")),
      end_day = as.integer(format(dates[length(dates)], "%d"))
    )
    source <- "time_axis_gregorian"
    return(list(parsed = parsed, source = source, note = note))
  }

  origin_date <- tryCatch(parse_time_origin_date(time_units[1L]), error = function(e) e)
  if (!inherits(origin_date, "error")) {
    first_offset_days <- as.integer(round(min(time_values) * scale_days))
    last_offset_days <- as.integer(round(max(time_values) * scale_days))
    start_date <- origin_date + first_offset_days
    end_date <- origin_date + last_offset_days

    parsed <- list(
      ok = TRUE,
      start_year = as.integer(format(start_date, "%Y")),
      start_month = as.integer(format(start_date, "%m")),
      start_day = as.integer(format(start_date, "%d")),
      end_year = as.integer(format(end_date, "%Y")),
      end_month = as.integer(format(end_date, "%m")),
      end_day = as.integer(format(end_date, "%d"))
    )
    source <- "time_units_origin_estimated"
    note <- "Date range estimated from time units origin because calendar decode was unavailable."
    return(list(parsed = parsed, source = source, note = note))
  }

  helper_error <- NULL
  if (!is.null(helper_try) && inherits(helper_try, "error")) {
    helper_error <- helper_try$message
  }
  stop(
    "Could not parse start/end dates from filename and failed to derive them from time axis for file ",
    basename(file_path),
    ". helper_error=",
    if (is.null(helper_error)) "none" else helper_error,
    ", origin_parse_error=",
    origin_date$message,
    "."
  )
}

#' Build legacy daily time-selection metadata and completeness diagnostics
#'
#' Supports helper-based or empirical mapping, with Gregorian dates decoded from
#' UTC offsets and non-Gregorian empirical years inferred from filename start
#' plus day offsets. Full-year non-Gregorian files must start on January 1.
#'
#' Important: Gregorian/helper \code{is_full_year} uses modal observed counts
#' within leap/non-leap groups, not only theoretical calendar lengths. A repeated
#' missing-day pattern can therefore be classified full under the legacy mapping.
#' \code{declared_calendar_mismatch} exposes this; it is not proof of complete
#' dates. Corrected-data preparation separately verifies exact CF date coverage.
#' Warn policy preserves the historical experiment behaviour; it does not repair
#' or approve missing data. \code{time_phase} is a legacy NA placeholder.
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param file_path Path containing a date-range hint, used in diagnostics.
#' @param time_dim_name Time dimension name.
#' @param calendar Canonical CF calendar name.
#' @param required_years Integer vector of requested years.
#' @param strict_time_checks Whether invalid ordering/fractional non-Gregorian
#'   offsets raise errors; gaps may still warn and receive per-year checks.
#' @param time_decode_mode One of auto, helpers_cf (no decode fallback), empirical.
#' @param requested_year_policy Warn or stop when requested years fail the active
#'   mapping's completeness test. Defaults to warn, unchanged from production.
#' @return List containing backend/fallback diagnostics, \code{year_table},
#'   file coverage, and per-record \code{timestep_year}, \code{timestep_slot}
#'   and \code{timestep_doy}. Slots can be observed-record ranks, not true DOYs.
infer_time_structure_daily <- function(
  nc_obj,
  file_path,
  time_dim_name,
  calendar,
  required_years,
  strict_time_checks = TRUE,
  time_decode_mode = c("auto", "helpers_cf", "empirical"),
  requested_year_policy = c("warn", "stop")
) {
  time_decode_mode <- match.arg(time_decode_mode)
  requested_year_policy <- match.arg(requested_year_policy)

  time_units <- ncdf4::ncatt_get(nc_obj, time_dim_name, "units")$value
  scale_days <- resolve_time_units_scale_days(time_units)

  time_values <- read_time_values(nc_obj, time_dim_name)
  n_time <- length(time_values)
  if (n_time == 0L) {
    stop("Time axis is empty in file: ", file_path)
  }
  if (any(!is.finite(time_values))) {
    stop("Time axis contains nonfinite values in file: ", file_path)
  }

  # Tolerances for floating-point time axes.
  offset_round_tol <- 1e-3
  step_days_min <- NA_real_
  step_days_max <- NA_real_
  step_days_median <- NA_real_
  has_duplicates_or_reverse <- FALSE
  daily_time_step <- TRUE
  step_daily_tol <- 1e-3
  max_offset_round_error <- NA_real_
  time_phase <- NA_real_
  day_number <- rep.int(1L, n_time)

  if (n_time > 1L) {
    step_days <- diff(time_values) * scale_days
    step_days_min <- min(step_days)
    step_days_max <- max(step_days)
    step_days_median <- stats::median(step_days)
    has_duplicates_or_reverse <- any(step_days <= 0)
    daily_time_step <- all(abs(step_days - 1) <= step_daily_tol)
  }

  # Map raw time to integer day numbers from the first timestep.
  # This remains robust if rare gaps exist later in the file.
  time_day_offset <- (time_values - time_values[1]) * scale_days
  day_offset_rounded <- round(time_day_offset)
  max_offset_round_error <- max(abs(time_day_offset - day_offset_rounded))
  day_number <- as.integer(day_offset_rounded) + 1L
  if (any(day_number < 1L)) {
    stop("Invalid time mapping produced day numbers < 1 in file: ", basename(file_path))
  }

  parsed_dates <- parse_filename_date_range(file_path)
  filename_dates_source <- "filename"
  if (!isTRUE(parsed_dates$ok)) {
    derived_range <- derive_date_range_from_time_axis(
      file_path = file_path,
      nc_obj = nc_obj,
      time_dim_name = time_dim_name,
      time_values = time_values,
      time_units = time_units,
      scale_days = scale_days,
      calendar = calendar,
      time_decode_mode = time_decode_mode
    )
    parsed_dates <- derived_range$parsed
    filename_dates_source <- derived_range$source
    warning(
      "Could not parse start/end dates from filename ",
      basename(file_path),
      ". Using ",
      filename_dates_source,
      " start/end dates instead.",
      call. = FALSE
    )
    if (!is.null(derived_range$note) && nzchar(derived_range$note)) {
      warning(derived_range$note, call. = FALSE)
    }
  }

  starts_on_jan1 <- parsed_dates$start_month == 1L && parsed_dates$start_day == 1L

  calendar_mode <- if (identical(calendar, "unknown")) "inferred_fixed" else "known_cf"
  inferred_year_length <- NA_integer_
  decode_backend <- "empirical"
  decode_backend_requested <- time_decode_mode
  fallback_used <- FALSE
  fallback_reason <- NA_character_
  helper_result_class <- NA_character_

  #' Derive observed native-calendar year slots from decoded daily dates
  #' Internal closure of infer_time_structure_daily; not a standalone entry point.
  #' @param date_values Date-like decoded time vector supporting year formatting.
  #' @param backend_tag Decoder identity used to label the returned year-length mode.
  #' @return List of timestep years/slots, observed years/counts and theoretical native-calendar year lengths using the enclosing calendar settings.
  build_vectors_from_dates <- function(date_values, backend_tag) {
    timestep_year_local <- as.integer(format(date_values, "%Y"))
    year_vector_local <- sort(unique(timestep_year_local))
    if (length(year_vector_local) == 0L) {
      stop("Decoded date vector did not produce any year values.")
    }

    year_index_of_timestep <- match(timestep_year_local, year_vector_local)
    year_lengths_local <- as.integer(tabulate(year_index_of_timestep, nbins = length(year_vector_local)))

    timestep_slot_local <- integer(length(timestep_year_local))
    for (k in seq_along(year_vector_local)) {
      idx <- which(year_index_of_timestep == k)
      if (length(idx) > 0L) {
        timestep_slot_local[idx] <- seq_along(idx)
      }
    }

    theoretical_local <- if (identical(calendar_mode, "known_cf")) {
      as.integer(days_in_year_by_calendar(year_vector_local, calendar))
    } else {
      year_lengths_local
    }

    list(
      timestep_year = timestep_year_local,
      timestep_slot = timestep_slot_local,
      year_vector = year_vector_local,
      year_lengths = year_lengths_local,
      theoretical_year_lengths = theoretical_local,
      year_length_mode = if (backend_tag == "helpers_cf") "helpers_cf_observed_by_year" else "empirical_by_year"
    )
  }

  #' Decode daily year slots using the enclosing empirical calendar controls
  #' Internal closure of infer_time_structure_daily; not a standalone entry point.
  #' @return List of timestep years/slots and year lengths. Uses enclosing time_values, units and calendar; updates inferred_year_length for unknown fixed calendars.
  empirical_decode <- function() {
    if (is_gregorian_like_calendar(calendar)) {
      date_values <- decode_time_to_dates_gregorian(time_values, time_units[1L])
      return(build_vectors_from_dates(date_values, backend_tag = "empirical"))
    }

    year_vector_local <- integer(0)
    year_lengths_local <- integer(0)
    theoretical_local <- integer(0)
    cumulative_year_end <- integer(0)
    max_day_number <- max(day_number)
    current_year <- parsed_dates$start_year

    if (identical(calendar_mode, "inferred_fixed")) {
      inferred_year_length <<- infer_unknown_year_length(parsed_dates, n_time)
    }

    while (length(cumulative_year_end) == 0L || tail(cumulative_year_end, 1) < max_day_number) {
      target_len <- if (identical(calendar_mode, "known_cf")) {
        as.integer(days_in_year_by_calendar(current_year, calendar))
      } else {
        inferred_year_length
      }
      year_vector_local <- c(year_vector_local, current_year)
      year_lengths_local <- c(year_lengths_local, target_len)
      theoretical_local <- c(theoretical_local, target_len)
      cumulative_year_end <- c(cumulative_year_end, sum(year_lengths_local))
      current_year <- current_year + 1L
    }

    cumulative_year_start <- cumulative_year_end - year_lengths_local + 1L
    year_index_of_timestep <- findInterval(day_number - 1L, cumulative_year_end) + 1L
    timestep_year_local <- year_vector_local[year_index_of_timestep]
    timestep_slot_local <- day_number - cumulative_year_start[year_index_of_timestep] + 1L

    list(
      timestep_year = timestep_year_local,
      timestep_slot = timestep_slot_local,
      year_vector = year_vector_local,
      year_lengths = year_lengths_local,
      theoretical_year_lengths = theoretical_local,
      year_length_mode = if (identical(calendar_mode, "known_cf")) "calendar_rule" else "inferred_fixed"
    )
  }

  #' Select the requested daily calendar decoder and record fallbacks
  #' Internal closure of infer_time_structure_daily; not a standalone entry point.
  #' @return Decoded year-slot list from the helper or empirical backend. Updates enclosing backend/fallback metadata and errors when a forced helper fails.
  choose_decode <- function() {
    if (identical(calendar_mode, "inferred_fixed") && !identical(time_decode_mode, "helpers_cf")) {
      decode_backend <<- "empirical"
      return(empirical_decode())
    }

    if (identical(time_decode_mode, "empirical")) {
      decode_backend <<- "empirical"
      return(empirical_decode())
    }

    helper_result <- tryCatch(
      decode_time_with_ncdf4_helpers(nc_obj = nc_obj, time_dim_name = time_dim_name),
      error = function(e) e
    )

    if (inherits(helper_result, "error")) {
      if (identical(time_decode_mode, "helpers_cf")) {
        stop(
          "helpers_cf decoding failed for file ",
          basename(file_path),
          ": ",
          helper_result$message
        )
      }

      fallback_used <<- TRUE
      fallback_reason <<- paste0("helpers_cf failed: ", helper_result$message)
      decode_backend <<- "empirical"
      return(empirical_decode())
    }

    helper_vectors <- build_vectors_from_dates(
      date_values = helper_result$date_values,
      backend_tag = "helpers_cf"
    )
    helper_result_class <<- helper_result$helper_result_class

    helper_declared_mismatch <- identical(calendar_mode, "known_cf") &&
      any(helper_vectors$year_lengths != helper_vectors$theoretical_year_lengths)

    if (identical(time_decode_mode, "auto") && helper_declared_mismatch) {
      fallback_used <<- TRUE
      fallback_reason <<- paste0(
        "helpers_cf produced year lengths inconsistent with declared calendar: observed=",
        paste(sort(unique(helper_vectors$year_lengths)), collapse = "/"),
        ", theoretical=",
        paste(sort(unique(helper_vectors$theoretical_year_lengths)), collapse = "/")
      )
      decode_backend <<- "empirical"
      return(empirical_decode())
    }

    decode_backend <<- "helpers_cf"
    helper_vectors
  }

  decode_vectors <- choose_decode()
  timestep_year <- decode_vectors$timestep_year
  timestep_slot <- decode_vectors$timestep_slot
  year_vector <- decode_vectors$year_vector
  year_lengths <- decode_vectors$year_lengths
  theoretical_year_lengths <- decode_vectors$theoretical_year_lengths
  year_length_mode <- decode_vectors$year_length_mode

  years_from_filename <- seq.int(parsed_dates$start_year, parsed_dates$end_year)
  if (identical(calendar_mode, "known_cf")) {
    expected_year_lengths <- as.integer(days_in_year_by_calendar(years_from_filename, calendar))
    expected_n_time <- sum(expected_year_lengths)
    end_day_expected <- days_in_month_by_calendar(parsed_dates$end_year, 12L, calendar)
    ends_on_year_end <- parsed_dates$end_month == 12L && parsed_dates$end_day == end_day_expected
  } else {
    inferred_year_length <- infer_unknown_year_length(parsed_dates, n_time)
    expected_year_lengths <- rep.int(inferred_year_length, length(years_from_filename))
    expected_n_time <- sum(expected_year_lengths)
    ends_on_year_end <- (n_time %% inferred_year_length) == 0L
  }
  full_year_coverage <- starts_on_jan1 && ends_on_year_end

  expected_length_matches <- identical(as.integer(expected_n_time), as.integer(n_time))

  year_table <- data.frame(
    year = year_vector,
    year_length = year_lengths,
    theoretical_year_length = theoretical_year_lengths,
    expected_timesteps_for_full_year = NA_integer_,
    start_index = NA_integer_,
    end_index = NA_integer_,
    observed_timesteps = 0L,
    observed_unique_doy = 0L,
    observed_min_doy = NA_integer_,
    observed_max_doy = NA_integer_,
    is_full_year = FALSE,
    stringsAsFactors = FALSE
  )

  for (k in seq_along(year_vector)) {
    idx <- which(timestep_year == year_vector[k])
    if (length(idx) == 0L) {
      next
    }

    slot_vals <- timestep_slot[idx]
    year_table$start_index[k] <- min(idx)
    year_table$end_index[k] <- max(idx)
    year_table$observed_timesteps[k] <- length(idx)
    year_table$observed_unique_doy[k] <- length(unique(slot_vals))
    year_table$observed_min_doy[k] <- min(slot_vals)
    year_table$observed_max_doy[k] <- max(slot_vals)
  }

  if (year_length_mode %in% c("empirical_by_year", "helpers_cf_observed_by_year")) {
    leap_mask <- year_table$theoretical_year_length == 366L
    expected_obs <- integer(nrow(year_table))

    for (flag in c(FALSE, TRUE)) {
      group_idx <- which(leap_mask == flag)
      if (length(group_idx) == 0L) {
        next
      }
      expected_obs[group_idx] <- mode_integer(year_table$observed_timesteps[group_idx])
    }

    year_table$expected_timesteps_for_full_year <- expected_obs
    contiguous_idx <- (
      !is.na(year_table$start_index) &
      !is.na(year_table$end_index) &
      (year_table$end_index - year_table$start_index + 1L == year_table$observed_timesteps)
    )
    year_table$is_full_year <- (
      contiguous_idx &
      year_table$observed_timesteps == year_table$expected_timesteps_for_full_year &
      year_table$observed_unique_doy == year_table$year_length &
      year_table$observed_min_doy == 1L &
      year_table$observed_max_doy == year_table$year_length
    )
  } else {
    year_table$expected_timesteps_for_full_year <- year_table$year_length
    year_table$is_full_year <- (
      year_table$observed_timesteps == year_table$year_length &
      year_table$observed_unique_doy == year_table$year_length &
      year_table$observed_min_doy == 1L &
      year_table$observed_max_doy == year_table$year_length
    )
  }

  declared_calendar_mismatch <- identical(calendar_mode, "known_cf") &&
    any(year_table$observed_timesteps != year_table$theoretical_year_length)
  theoretical_year_lengths_unique <- sort(unique(year_table$theoretical_year_length))
  observed_year_lengths_unique <- sort(unique(year_table$observed_timesteps))

  if (
    strict_time_checks &&
    !is_gregorian_like_calendar(calendar) &&
    max_offset_round_error > offset_round_tol
  ) {
    stop(
      "Fast fractional seasonal extraction assumptions failed for file ", basename(file_path), ". ",
      "time axis is not close to integer-day offsets from first timestep. ",
      "max_offset_round_error=", signif(max_offset_round_error, 7),
      ", offset_round_tol=", offset_round_tol,
      "."
    )
  }

  if (strict_time_checks && has_duplicates_or_reverse) {
    stop(
      "Fast fractional seasonal extraction assumptions failed for file ", basename(file_path), ". ",
      "duplicates_or_reverse=", has_duplicates_or_reverse,
      ", step_days_min=", signif(step_days_min, 7),
      ", step_days_median=", signif(step_days_median, 7),
      ", step_days_max=", signif(step_days_max, 7),
      ", step_daily_tol=", step_daily_tol,
      "."
    )
  }

  if (strict_time_checks && !starts_on_jan1 && identical(requested_year_policy, "stop")) {
    stop(
      "Fast fractional seasonal extraction assumptions failed for file ",
      basename(file_path),
      ". starts_on_jan1=",
      starts_on_jan1,
      "."
    )
  }

  if (
    !strict_time_checks &&
    (
      has_duplicates_or_reverse ||
      !starts_on_jan1 ||
      (!is_gregorian_like_calendar(calendar) && max_offset_round_error > offset_round_tol)
    )
  ) {
    warning(
      "Fast fractional extraction assumptions not fully satisfied for ", basename(file_path),
      ". Proceeding because strict_time_checks=FALSE.",
      call. = FALSE
    )
  }

  if (strict_time_checks && !daily_time_step) {
    warning(
      "Detected non-unit day steps in ", basename(file_path),
      " (min/median/max day step: ",
      signif(step_days_min, 7), "/",
      signif(step_days_median, 7), "/",
      signif(step_days_max, 7),
      "). Continuing because requested years are validated via year/day mapping.",
      call. = FALSE
    )
  }

  if (!starts_on_jan1 && identical(requested_year_policy, "warn")) {
    warning(
      "File ",
      basename(file_path),
      " does not start on January 1. Proceeding because requested_year_policy='warn'.",
      call. = FALSE
    )
  }

  if (declared_calendar_mismatch) {
    mismatch_note <- if (!starts_on_jan1 || !ends_on_year_end) {
      " This may be expected for pre-subset seasonal files that do not contain full years."
    } else {
      ""
    }
    warning(
      "Declared calendar '",
      calendar,
      "' in ",
      basename(file_path),
      " is inconsistent with observed per-year timestep counts. ",
      "Observed year lengths: ",
      paste(observed_year_lengths_unique, collapse = "/"),
      "; theoretical under declared calendar: ",
      paste(theoretical_year_lengths_unique, collapse = "/"),
      ". Using empirical year lengths for fractional seasonal selection.",
      mismatch_note,
      call. = FALSE
    )
  }

  year_table$requested_year <- year_table$year %in% required_years
  year_table$anomaly_class <- ifelse(
    !year_table$is_full_year & year_table$requested_year,
    "requested_year_anomaly",
    ifelse(!year_table$is_full_year, "out_of_range_anomaly", "none")
  )

  assumptions_ok <-
    daily_time_step &&
    !has_duplicates_or_reverse &&
    starts_on_jan1 &&
    expected_length_matches &&
    full_year_coverage

  if (!strict_time_checks && !assumptions_ok) {
    warning(
      "Filename-derived coverage checks differ from actual timestep count for ",
      basename(file_path),
      ". full_year_coverage=",
      full_year_coverage,
      ", expected_n_time=",
      expected_n_time,
      ", actual_n_time=",
      n_time,
      ".",
      ", max_offset_round_error=",
      signif(max_offset_round_error, 7),
      call. = FALSE
    )
  }

  full_years_available <- year_table$year[year_table$is_full_year]
  anomaly_years <- year_table$year[!year_table$is_full_year]
  anomaly_years_requested <- sort(intersect(anomaly_years, required_years))
  anomaly_years_outside <- sort(setdiff(anomaly_years, required_years))

  if (length(full_years_available) == 0L) {
    year_diag_preview <- head(
      year_table[
        ,
        c(
          "year",
          "year_length",
          "theoretical_year_length",
          "expected_timesteps_for_full_year",
          "observed_timesteps",
          "observed_unique_doy",
          "observed_min_doy",
          "observed_max_doy"
        )
      ],
      10L
    )
    msg <- paste0(
      "File ",
      basename(file_path),
      " does not contain any fully covered model year under calendar '",
      calendar,
      "' with year_length_mode='",
      year_length_mode,
      "'. Example yearly diagnostics:\n",
      paste(capture.output(print(year_diag_preview, row.names = FALSE)), collapse = "\n")
    )
    if (identical(requested_year_policy, "stop")) {
      stop(msg)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  missing_required <- setdiff(required_years, full_years_available)
  if (length(missing_required) > 0L) {
    msg <- if (identical(calendar_mode, "inferred_fixed")) {
      paste0(
        "Requested years are not fully available after inferring unknown calendar for file ",
        basename(file_path),
        ". inferred_year_length=",
        inferred_year_length,
        ", missing requested years=",
        paste(missing_required, collapse = ", "),
        "."
      )
    } else {
      paste0(
        "Requested years are not fully available in file ",
        basename(file_path),
        ": ",
        paste(missing_required, collapse = ", "),
        ". Fully available years are ",
        if (length(full_years_available) == 0L) "none" else paste0(min(full_years_available), "-", max(full_years_available)),
        ". year_length_mode=",
        year_length_mode,
        "."
      )
    }

    if (identical(requested_year_policy, "stop")) {
      stop(msg)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  list(
    time_units = as.character(time_units[1]),
    n_time = n_time,
    daily_time_step = daily_time_step,
    has_duplicates_or_reverse = has_duplicates_or_reverse,
    step_days_min = step_days_min,
    step_days_median = step_days_median,
    step_days_max = step_days_max,
    step_daily_tol = step_daily_tol,
    max_offset_round_error = max_offset_round_error,
    time_phase = time_phase,
    starts_on_jan1 = starts_on_jan1,
    full_year_coverage = full_year_coverage,
    ends_on_year_end = ends_on_year_end,
    expected_length_matches = expected_length_matches,
    assumptions_ok = assumptions_ok,
    decode_backend = decode_backend,
    decode_backend_requested = decode_backend_requested,
    fallback_used = fallback_used,
    fallback_reason = fallback_reason,
    helper_result_class = helper_result_class,
    calendar_mode = calendar_mode,
    year_length_mode = year_length_mode,
    declared_calendar_mismatch = declared_calendar_mismatch,
    theoretical_year_lengths_unique = theoretical_year_lengths_unique,
    observed_year_lengths_unique = observed_year_lengths_unique,
    inferred_year_length = inferred_year_length,
    year_table = year_table,
    requested_year_policy = requested_year_policy,
    anomaly_years = anomaly_years,
    anomaly_years_requested = anomaly_years_requested,
    anomaly_years_outside = anomaly_years_outside,
    filename_year_start = parsed_dates$start_year,
    filename_year_end = parsed_dates$end_year,
    filename_dates_source = filename_dates_source,
    filename_start_md = c(parsed_dates$start_month, parsed_dates$start_day),
    filename_end_md = c(parsed_dates$end_month, parsed_dates$end_day),
    timestep_year = timestep_year,
    timestep_slot = timestep_slot,
    timestep_doy = timestep_slot
  )
}
#' Locate the first recognised spatial axis
#' @param dimension_names Character vector of variable dimension names.
#' @param candidates Accepted names for one spatial axis.
#' @param label Axis label used if no candidate is found.
#' @return One-based index of the first match; missing matches raise an error.
resolve_spatial_dimension_index <- function(dimension_names, candidates, label) {
  idx <- which(dimension_names %in% candidates)
  if (length(idx) == 0L) {
    stop(
      "Could not find ", label,
      " dimension. Available dimensions: ",
      paste(dimension_names, collapse = ", "),
      "."
    )
  }
  idx[1]
}

#' Describe the axes of a three-dimensional rectilinear climate variable
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param variable NetCDF variable name, required to have exactly three axes.
#' @param time_dim_name Name of that variable's time axis.
#' @return List with original dimension names/sizes, longitude/latitude/time
#'   positions and names. Accepted spatial aliases are lon/longitude/x and lat/latitude/y.
get_variable_dimension_info <- function(nc_obj, variable, time_dim_name) {
  var_obj <- nc_obj$var[[variable]]
  if (is.null(var_obj)) {
    stop("Variable '", variable, "' not found in NetCDF file.")
  }

  dimension_names <- vapply(var_obj$dim, function(dim_obj) dim_obj$name, character(1))
  if (length(dimension_names) != 3L) {
    stop(
      "Variable '", variable,
      "' must have exactly 3 dimensions (lon, lat, time). Found: ",
      paste(dimension_names, collapse = ", "),
      "."
    )
  }

  lon_pos <- resolve_spatial_dimension_index(
    dimension_names,
    candidates = c("lon", "longitude", "x"),
    label = "longitude"
  )
  lat_pos <- resolve_spatial_dimension_index(
    dimension_names,
    candidates = c("lat", "latitude", "y"),
    label = "latitude"
  )
  time_pos <- which(dimension_names == time_dim_name)

  if (length(time_pos) == 0L) {
    stop(
      "Time dimension '", time_dim_name,
      "' not found in variable '", variable,
      "'. Variable dimensions: ",
      paste(dimension_names, collapse = ", "),
      "."
    )
  }

  list(
    dimension_names = dimension_names,
    variable_size = as.integer(var_obj$varsize),
    lon_pos = lon_pos,
    lat_pos = lat_pos,
    time_pos = time_pos[1],
    lon_dim_name = dimension_names[lon_pos],
    lat_dim_name = dimension_names[lat_pos],
    time_dim_name = time_dim_name
  )
}

#' Match coordinates exactly or by one unambiguous nearby value
#' @param requested Numeric coordinate vector in desired order.
#' @param available Numeric coordinates on the file axis.
#' @param tolerance Absolute tolerance, in coordinate units, for fallback matching.
#' @return Integer indices into available. Unmatched or ambiguous nonexact
#'   requests remain NA; exact matching uses R's first-match convention.
#' @examples
#' stopifnot(match_values_with_tolerance(1+1e-8, 1:3) == 1L)
match_values_with_tolerance <- function(requested, available, tolerance = 1e-6) {
  idx <- match(requested, available)
  missing <- which(is.na(idx))

  if (length(missing) > 0L) {
    for (k in missing) {
      candidates <- which(abs(available - requested[k]) <= tolerance)
      if (length(candidates) == 1L) {
        idx[k] <- candidates[1]
      }
    }
  }

  idx
}

#' Match requested rectilinear coordinates without regridding
#'
#' Longitudes at least 180 are shifted once by -360 in both request and file.
#' This handles the Article 1 0..359/-180..179 conventions, not arbitrary turns
#' of longitude. Output indices preserve requested order, including latitude order.
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param dimension_info Axis description from \code{get_variable_dimension_info()}.
#' @param lon Requested longitude vector in degrees east.
#' @param lat Requested latitude vector in degrees north.
#' @return List with integer \code{lon_indices} and \code{lat_indices} into
#'   original file axes. Missing/ambiguous coordinates raise errors.
match_requested_coordinates <- function(nc_obj, dimension_info, lon, lat) {
  lon_values <- read_dimension_values(nc_obj, dimension_info$lon_dim_name)
  lat_values <- read_dimension_values(nc_obj, dimension_info$lat_dim_name)

  lon_file_adj <- ifelse(lon_values >= 180, lon_values - 360, lon_values)
  lon_req_adj <- ifelse(lon >= 180, lon - 360, lon)

  lon_order <- order(lon_file_adj)
  lon_sorted <- lon_file_adj[lon_order]

  lon_idx_sorted <- match_values_with_tolerance(lon_req_adj, lon_sorted)
  if (anyNA(lon_idx_sorted)) {
    missing_lons <- unique(lon[is.na(lon_idx_sorted)])
    stop(
      "Requested longitude(s) not found in file grid: ",
      paste(head(missing_lons, 8), collapse = ", "),
      if (length(missing_lons) > 8L) " ..." else ""
    )
  }

  lon_idx_file <- lon_order[lon_idx_sorted]

  lat_idx_file <- match_values_with_tolerance(lat, lat_values)
  if (anyNA(lat_idx_file)) {
    missing_lats <- unique(lat[is.na(lat_idx_file)])
    stop(
      "Requested latitude(s) not found in file grid: ",
      paste(head(missing_lats, 8), collapse = ", "),
      if (length(missing_lats) > 8L) " ..." else ""
    )
  }

  list(
    lon_indices = as.integer(lon_idx_file),
    lat_indices = as.integer(lat_idx_file)
  )
}

#' Map a non-leap seasonal interval proportionally onto a model year
#'
#' Uses 1 + floor((start-1)*L/365) and floor(end*L/365), clamped to 1..L.
#' A whole-year request always maps to every slot. A wrapped interval is the
#' union of the end and start segments within each selected calendar year,
#' not a season reassigned to the following year's label.
#' @param year_length Positive integer model/empirical year length L.
#' @param start_ref_doy One-based start day on the 365-day reference calendar.
#' @param end_ref_doy One-based end day on the same reference calendar.
#' @param wraps_year_end Logical scalar indicating an interval across December/January.
#' @return List of integer-like scalars \code{start_mod} and \code{end_mod}.
#' @examples
#' stopifnot(map_fractional_window_to_year_length(360,1,365,FALSE)$end_mod == 360)
map_fractional_window_to_year_length <- function(year_length, start_ref_doy, end_ref_doy, wraps_year_end) {
  year_length <- as.integer(year_length)
  start_mod <- 1L + floor((start_ref_doy - 1L) * year_length / 365)
  end_mod <- floor(end_ref_doy * year_length / 365)

  start_mod <- max(1L, min(year_length, as.integer(start_mod)))
  end_mod <- max(1L, min(year_length, as.integer(end_mod)))

  if (!wraps_year_end && end_mod < start_mod) {
    end_mod <- start_mod
  }

  list(start_mod = start_mod, end_mod = end_mod)
}

#' Select seasonal file indices under the fractional year-slot convention
#' @param year_table Year diagnostics from \code{infer_time_structure_daily()}.
#' @param timestep_year Integer year label for every file record.
#' @param timestep_slot Corresponding one-based slot under the active time mapping.
#' @param selected_years Requested years, sorted and deduplicated before selection.
#' @param start_ref_doy Seasonal start day on a 365-day reference year.
#' @param end_ref_doy Seasonal end day on that reference year.
#' @param wraps_year_end Whether to take both ends of each selected year.
#' @param requested_year_policy Warn or stop for absent/incomplete/empty years.
#'   Warn mode uses observed maximum slot as effective length for incomplete years.
#' @return List with integer \code{indices} in chronological selected-year
#'   order and a \code{year_summary} table. Errors if nothing can be selected.
build_fractional_indices_for_years <- function(
  year_table,
  timestep_year,
  timestep_slot,
  selected_years,
  start_ref_doy,
  end_ref_doy,
  wraps_year_end,
  requested_year_policy = c("warn", "stop")
) {
  requested_year_policy <- match.arg(requested_year_policy)
  selected_years <- sort(unique(as.integer(selected_years)))

  summaries <- vector("list", length(selected_years))
  all_indices <- vector("list", length(selected_years))

  for (k in seq_along(selected_years)) {
    year_value <- selected_years[k]
    row_idx <- match(year_value, year_table$year)

    if (is.na(row_idx)) {
      if (identical(requested_year_policy, "stop")) {
        stop("Requested year ", year_value, " is not available in year table.")
      } else {
        warning(
          "Requested year ",
          year_value,
          " is not present in the year table. Skipping because requested_year_policy='warn'.",
          call. = FALSE
        )
        next
      }
    }

    year_length <- as.integer(year_table$year_length[row_idx])
    if (!isTRUE(year_table$is_full_year[row_idx])) {
      if (identical(requested_year_policy, "stop")) {
        stop(
          "Requested year ",
          year_value,
          " is not fully covered in the time axis under the active year-length mapping; ",
          "cannot apply strict fractional seasonal selection."
        )
      } else {
        warning(
          "Requested year ",
          year_value,
          " is not fully covered in the time axis. Proceeding because requested_year_policy='warn'.",
          call. = FALSE
        )
      }
    }

    effective_year_length <- if (isTRUE(year_table$is_full_year[row_idx])) {
      year_length
    } else {
      max(1L, as.integer(year_table$observed_max_doy[row_idx]))
    }

    mapped <- map_fractional_window_to_year_length(
      year_length = effective_year_length,
      start_ref_doy = start_ref_doy,
      end_ref_doy = end_ref_doy,
      wraps_year_end = wraps_year_end
    )

    sel_mask <- timestep_year == year_value
    if (wraps_year_end) {
      sel_mask <- sel_mask & (timestep_slot >= mapped$start_mod | timestep_slot <= mapped$end_mod)
    } else {
      sel_mask <- sel_mask & (timestep_slot >= mapped$start_mod & timestep_slot <= mapped$end_mod)
    }

    global_indices <- which(sel_mask)
    if (length(global_indices) == 0L) {
      if (identical(requested_year_policy, "stop")) {
        stop(
          "No timesteps selected for requested year ",
          year_value,
          " after fractional mapping. This indicates an inconsistent year/day mapping."
        )
      } else {
        warning(
          "No timesteps selected for requested year ",
          year_value,
          " after fractional mapping. Skipping this year because requested_year_policy='warn'.",
          call. = FALSE
        )
        next
      }
    }
    all_indices[[k]] <- global_indices

    summaries[[k]] <- data.frame(
      year = year_value,
      year_length = year_length,
      effective_year_length = effective_year_length,
      start_mod = mapped$start_mod,
      end_mod = mapped$end_mod,
      selected_timesteps = length(global_indices),
      wraps_window = wraps_year_end,
      year_full = isTRUE(year_table$is_full_year[row_idx]),
      stringsAsFactors = FALSE
    )
  }

  summaries <- Filter(Negate(is.null), summaries)
  all_indices <- Filter(Negate(is.null), all_indices)
  combined_indices <- unlist(all_indices, use.names = FALSE)
  if (length(combined_indices) == 0L) {
    stop("No timesteps were selected for requested years.")
  }

  list(
    indices = as.integer(combined_indices),
    year_summary = do.call(rbind, summaries)
  )
}

#' Test membership in an inclusive seasonal month/day interval
#' @param month Numeric month vector of decoded dates.
#' @param day Numeric day vector parallel to month.
#' @param start_md Validated month/day start boundary.
#' @param end_md Validated month/day end boundary.
#' @param wraps_year_end Whether the interval joins December to January.
#' @return Logical membership vector. Uses date ordering, not calendar-day arithmetic.
#' @examples
#' stopifnot(all(month_day_in_window(c(12,1),c(31,1),c(12,1),c(2,28),TRUE)))
month_day_in_window <- function(month, day, start_md, end_md, wraps_year_end) {
  md_ge_start <- (month > start_md[1L]) | (month == start_md[1L] & day >= start_md[2L])
  md_le_end <- (month < end_md[1L]) | (month == end_md[1L] & day <= end_md[2L])

  if (wraps_year_end) {
    md_ge_start | md_le_end
  } else {
    md_ge_start & md_le_end
  }
}

#' Try to obtain year/month/day vectors for exact seasonal selection
#' @param nc_obj Open NetCDF handle, not closed here.
#' @param time_dim_name Time dimension supplied to optional CF helpers.
#' @param time_values Numeric offsets for Gregorian fallback decoding.
#' @param time_units CF time units with UTC origin.
#' @param calendar Canonical calendar name.
#' @param time_decode_mode Auto/helpers_cf/empirical; empirical skips optional helpers.
#' @return List with \code{ok}, decoded \code{year}, \code{month}, \code{day},
#'   \code{source}, and \code{error_message}. Non-Gregorian helper failure
#'   returns ok=FALSE; this routine does not invent Gregorian model dates.
decode_time_components_for_exact_selection <- function(
  nc_obj,
  time_dim_name,
  time_values,
  time_units,
  calendar,
  time_decode_mode
) {
  helper_error <- NULL
  if (!identical(time_decode_mode, "empirical")) {
    helper_result <- tryCatch(
      decode_time_with_ncdf4_helpers(nc_obj = nc_obj, time_dim_name = time_dim_name),
      error = function(e) e
    )

    if (!inherits(helper_result, "error")) {
      date_values <- helper_result$date_values
      return(list(
        ok = TRUE,
        year = as.integer(format(date_values, "%Y")),
        month = as.integer(format(date_values, "%m")),
        day = as.integer(format(date_values, "%d")),
        source = "helpers_cf",
        error_message = NA_character_
      ))
    }
    helper_error <- helper_result$message
  }

  if (is_gregorian_like_calendar(calendar)) {
    date_values <- decode_time_to_dates_gregorian(time_values, time_units[1L])
    return(list(
      ok = TRUE,
      year = as.integer(format(date_values, "%Y")),
      month = as.integer(format(date_values, "%m")),
      day = as.integer(format(date_values, "%d")),
      source = "gregorian_time_units",
      error_message = NA_character_
    ))
  }

  list(
    ok = FALSE,
    year = NULL,
    month = NULL,
    day = NULL,
    source = NA_character_,
    error_message = paste0(
      "Exact month/day selection requires decodable dates. calendar=",
      calendar,
      ", time_decode_mode=",
      time_decode_mode,
      if (!is.null(helper_error)) paste0(", helpers_cf_error=", helper_error) else ""
    )
  )
}

#' Select inclusive seasonal dates from decoded month/day components
#'
#' Wrapped seasons include both ends of each selected year; no cross-year
#' reassignment is made. Completeness diagnostics come from year_table; this
#' routine checks missing years and empty selections, not every expected date.
#' @param year_table Year diagnostics from \code{infer_time_structure_daily()}.
#' @param timestep_year Decoded year for each file record.
#' @param timestep_slot Parallel slot/rank vector used in diagnostics only.
#' @param timestep_month Parallel decoded month vector.
#' @param timestep_day Parallel decoded day vector.
#' @param selected_years Requested years, sorted and deduplicated.
#' @param season_start_md Inclusive start boundary on the 365-day reference calendar.
#' @param season_end_md Inclusive end boundary on that calendar.
#' @param requested_year_policy Warn or stop for missing years or empty selections.
#' @return List with integer file \code{indices} and \code{year_summary} table;
#'   raises an error if the combined selection is empty.
build_exact_month_day_indices_for_years <- function(
  year_table,
  timestep_year,
  timestep_slot,
  timestep_month,
  timestep_day,
  selected_years,
  season_start_md,
  season_end_md,
  requested_year_policy = c("warn", "stop")
) {
  requested_year_policy <- match.arg(requested_year_policy)
  selected_years <- sort(unique(as.integer(selected_years)))
  wraps_year_end <- month_day_to_reference_doy(season_end_md) < month_day_to_reference_doy(season_start_md)

  summaries <- vector("list", length(selected_years))
  all_indices <- vector("list", length(selected_years))

  for (k in seq_along(selected_years)) {
    year_value <- selected_years[k]
    row_idx <- match(year_value, year_table$year)

    if (is.na(row_idx)) {
      if (identical(requested_year_policy, "stop")) {
        stop("Requested year ", year_value, " is not available in year table.")
      } else {
        warning(
          "Requested year ",
          year_value,
          " is not present in the year table. Skipping because requested_year_policy='warn'.",
          call. = FALSE
        )
        next
      }
    }

    in_year <- timestep_year == year_value
    in_window <- month_day_in_window(
      month = timestep_month,
      day = timestep_day,
      start_md = season_start_md,
      end_md = season_end_md,
      wraps_year_end = wraps_year_end
    )
    global_indices <- which(in_year & in_window)

    if (length(global_indices) == 0L) {
      if (identical(requested_year_policy, "stop")) {
        stop(
          "No timesteps selected for requested year ",
          year_value,
          " in exact month/day mode. Check decoded calendar dates."
        )
      } else {
        warning(
          "No timesteps selected for requested year ",
          year_value,
          " in exact month/day mode. Skipping because requested_year_policy='warn'.",
          call. = FALSE
        )
        next
      }
    }

    year_length <- as.integer(year_table$year_length[row_idx])
    effective_year_length <- if (isTRUE(year_table$is_full_year[row_idx])) {
      year_length
    } else {
      max(1L, as.integer(year_table$observed_max_doy[row_idx]))
    }

    selected_slots <- timestep_slot[global_indices]
    all_indices[[k]] <- global_indices
    summaries[[k]] <- data.frame(
      year = year_value,
      year_length = year_length,
      effective_year_length = effective_year_length,
      start_mod = min(selected_slots),
      end_mod = max(selected_slots),
      selected_timesteps = length(global_indices),
      wraps_window = wraps_year_end,
      year_full = isTRUE(year_table$is_full_year[row_idx]),
      selection_mode = "exact_month_day",
      stringsAsFactors = FALSE
    )
  }

  summaries <- Filter(Negate(is.null), summaries)
  all_indices <- Filter(Negate(is.null), all_indices)
  combined_indices <- unlist(all_indices, use.names = FALSE)
  if (length(combined_indices) == 0L) {
    stop("No timesteps were selected for requested years.")
  }

  list(
    indices = as.integer(combined_indices),
    year_summary = do.call(rbind, summaries)
  )
}

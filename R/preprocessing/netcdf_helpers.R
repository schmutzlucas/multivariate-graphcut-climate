# Reusable NetCDF metadata, output creation and numerical equality helpers.
#' Hash an input through a binary connection, including Windows SMB files
#' @param path Single existing readable file path.
#' @return Lowercase hexadecimal SHA-256 string. Requires openssl and reads the
#'   entire file without modifying it; always closes its connection, including on error.
repair_sha256 <- function(path) {
  # On the Windows/SMB validation host file.access(path,4) returns -1 for files
  # that can actually be opened. digest(file=...) refuses those paths. Hash an
  # explicitly opened binary connection instead; never skip checksum validation.
  if (!requireNamespace("openssl",quietly=TRUE)) stop("Recovery needs the openssl R package")
  connection <- file(path,open="rb")
  on.exit(close(connection))
  paste(as.character(openssl::sha256(connection)),collapse="")
}

#' Require and read a NetCDF attribute during preparation
#' @param nc Open ncdf4 handle, not closed here.
#' @param variable Variable name or zero for a global attribute.
#' @param name Attribute name.
#' @return Attribute value coerced to character; missing attributes raise an error.
repair_attr <- function(nc, variable, name) {
  value <- ncdf4::ncatt_get(nc, variable, name)
  if (!isTRUE(value$hasatt)) stop("Missing NetCDF attribute: ", name)
  as.character(value$value)
}

#' Decode and validate an ordered, unique daily-date CF axis
#' @param nc Open NetCDF handle with supported CF time units and calendar attributes.
#' @return List of time name/values/units/calendar, independently decoded date
#'   components and YYYY-MM-DD keys. Duplicate or non-increasing dates fail;
#'   completeness and a daily timestamp phase are checked separately by the preparation workflow.
repair_axis <- function(nc) {
  name <- find_time_dimension_name(nc)
  values <- read_time_values(nc, name)
  units <- repair_attr(nc, name, "units")
  calendar <- normalize_calendar_name(repair_attr(nc, name, "calendar"))
  dates <- audit_cf_dates(values, units, calendar)
  keys <- sprintf("%04d-%02d-%02d", dates$year, dates$month, dates$day)
  if (anyDuplicated(keys) || is.unsorted(keys, strictly = TRUE))
    stop("Duplicate or unordered source dates")
  list(name=name, values=values, units=units, calendar=calendar, dates=dates, keys=keys)
}

#' Create a new rectilinear NetCDF container for a recovered input
#' @param nc Open template handle with two rectilinear spatial axes followed by time.
#' @param variable Float/double source variable name; packed integer variables are unsupported.
#' @param path New output path; existing files are refused.
#' @param times Full numeric output time coordinates supplied by the recovery plan.
#' @param units CF time-coordinate units, not climate-variable units.
#' @param calendar Canonical CF calendar for the output time axis.
#' @param output_units Optional climate-variable units override; defaults to its template attribute.
#' @return Open writable ncdf4 handle, which the caller must close. Writes a
#'   NetCDF4 definition with compression 1, full-spatial/one-time chunks, original
#'   precision/missing value, and selected source/member/experiment/grid attributes.
#'   Does not fill the climate data or copy every template metadata attribute.
repair_new_file <- function(nc, variable, path, times, units, calendar, output_units=NULL) {
  if (file.exists(path)) stop("Refusing to overwrite NetCDF: ", path)
  definition <- nc$var[[variable]]
  dimensions <- definition$dim
  if (length(dimensions) != 3L || dimensions[[3]]$name != find_time_dimension_name(nc))
    stop("Recovery currently requires two spatial dimensions followed by time")
  spatial <- lapply(dimensions[1:2], function(d) {
    # A regular or Gaussian 1-D longitude/latitude grid is required. Never invent
    # a curvilinear grid from its dimension indices.
    if (!d$name %in% c("lon","lat","longitude","latitude"))
      stop("Unsupported spatial coordinate: ",d$name)
    ncdf4::ncdim_def(d$name, d$units, d$vals)
  })
  time <- ncdf4::ncdim_def(dimensions[[3]]$name, units, times, calendar=calendar)
  precision <- definition$prec
  if (!precision %in% c("float","double")) stop("Unsupported packed data precision")
  if (is.null(output_units)) output_units <- repair_attr(nc,variable,"units")
  var <- ncdf4::ncvar_def(variable, output_units, c(spatial,list(time)),
    missval=definition$missval, prec=precision,compression=1,
    chunksizes=c(dimensions[[1]]$len,dimensions[[2]]$len,1L))
  result <- ncdf4::nc_create(path,var,force_v4=TRUE)
  for (name in c("source_id","variant_label","experiment_id","grid_label","Conventions")) {
    value <- ncdf4::ncatt_get(nc,0,name)
    if (isTRUE(value$hasatt)) ncdf4::ncatt_put(result,0,name,value$value)
  }
  result
}

#' Require exact preservation of existing decoded climate values
#' @param actual,expected Numeric arrays to compare, with matching dimensions and NA pattern.
#' @param context Diagnostic description included in failure messages.
#' @return Number of expected values invisibly. Errors on shape/missing-pattern
#'   mismatch, non-missing nonfinite actual values, or any nonzero absolute difference.
#'   This is exact value comparison, not a configurable numerical tolerance.
#' @examples
#' repair_equal(array(c(1,NA),c(1,1,2)),array(c(1,NA),c(1,1,2)),"example")
repair_equal <- function(actual, expected, context) {
  if (!identical(dim(actual),dim(expected)) || !identical(is.na(actual),is.na(expected)))
    stop("Shape or missing-value pattern differs: ",context)
  difference <- abs(actual-expected)
  if (any(!is.finite(actual) & !is.na(actual)) || any(difference > 0,na.rm=TRUE))
    stop("Existing values differ: ",context)
  invisible(length(expected))
}

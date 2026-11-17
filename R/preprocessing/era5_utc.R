# Small, testable contracts shared by the monthly and merged UTC validators.

#' Fingerprint the complete decoded following-midnight ERA5 precipitation map
#' @param values Numeric longitude x latitude x one-hour array, in native metres.
#'   Coordinates/units/calendar must already have passed the native-grid checks.
#' @param utc_seconds Single finite UTC timestamp, in seconds since the Unix epoch.
#' @return Lowercase SHA-256 over dimensions, decoded double values and timestamp,
#'   using R's version-2 XDR serialisation. This hashes the actual required midnight,
#'   not unused later hours or container metadata in the following monthly file.
#'   Finite negative hourly increments are fingerprinted, not silently modified;
#'   the resulting daily product separately requires nonnegative precipitation.
era5_midnight_fingerprint <- function(values, utc_seconds) {
  if (!is.numeric(values) || length(dim(values)) != 3L || dim(values)[3] != 1L ||
      !length(values) || any(!is.finite(values)) || !is.numeric(utc_seconds) ||
      length(utc_seconds) != 1L || !is.finite(utc_seconds))
    stop("Expected one finite hourly map and UTC timestamp for ERA5 identity.")
  payload <- list(dimensions = dim(values), values = as.double(values), utc_seconds = as.double(utc_seconds))
  unclass(as.character(openssl::sha256(serialize(payload, NULL, version = 2, xdr = TRUE))))
}

#' Require physically valid daily ERA5 precipitation without clipping
#' @param values Nonempty numeric daily precipitation vector/array in mm per day.
#' @return List minimum, maximum and values_scanned after complete finiteness and
#'   nonnegativity checks. Any negative value, even a small one, fails explicitly;
#'   source investigation must precede any separately justified tolerance policy.
#'   Does not write or modify the supplied values or the original data archive.
#' @examples
#' stopifnot(validate_era5_precipitation_values(c(0, 1, 2))$maximum == 2)
validate_era5_precipitation_values <- function(values) {
  if (!is.numeric(values) || !length(values) || any(!is.finite(values)))
    stop("Non-finite or empty daily ERA5 precipitation; no filling permitted.")
  if (any(values < 0)) stop("Negative daily ERA5 precipitation; no clipping permitted.")
  list(minimum = min(values), maximum = max(values), values_scanned = length(values))
}

#' Check every merged ERA5 daily value against its signed monthly input
#' @param final_nc Open caller-owned ncdf4 handle for the merged daily pr product.
#' @param expected_keys Complete ordered YYYY-MM-DD keys already verified for final_nc.
#' @param sources Ordered monthly records with YYYYMM month, file, bytes and sha256.
#' @return Total compared values invisibly. Requires every expected month exactly
#'   once, fresh monthly payload hashes and matching dates/coordinates/units.
#'   Uses exact decoded-value equality, not a numerical tolerance or summary mean.
#'   Opens/closes monthly handles and reads at most one monthly cube at a time;
#'   leaves final_nc open. Requires the NetCDF helpers repair_axis/repair_equal
#'   definitions, loaded explicitly by the UTC validation entry point.
verify_era5_monthly_merge <- function(final_nc, expected_keys, sources) {
  expected_months <- gsub("-", "", unique(substr(expected_keys, 1, 7)), fixed = TRUE)
  if (!is.list(sources) || !length(sources) ||
      !identical(vapply(sources, function(x) x$month, ""), expected_months))
    stop("Merged ERA5 verification requires every monthly source exactly once in order.")
  compared <- 0
  for (item in sources) {
    if (!file.exists(item$file) || file.info(item$file)$size != item$bytes || mvgc_sha256(item$file) != item$sha256)
      stop("Monthly ERA5 source checksum/size mismatch: ", item$month)
    n <- ncdf4::nc_open(item$file)
    tryCatch({
      axis <- repair_axis(n)
      keys <- expected_keys[gsub("-", "", substr(expected_keys, 1, 7), fixed = TRUE) == item$month]
      if (!identical(axis$keys, keys) || repair_attr(n, "pr", "units") != "mm day-1")
        stop("Monthly ERA5 date/unit mismatch: ", item$month)
      for (name in c("lon", "lat")) if (!identical(as.numeric(ncdf4::ncvar_get(n, name)), as.numeric(ncdf4::ncvar_get(final_nc, name))))
        stop("Monthly ERA5 grid mismatch: ", item$month)
      indices <- match(keys, expected_keys)
      source <- ncdf4::ncvar_get(n, "pr", collapse_degen = FALSE)
      actual <- ncdf4::ncvar_get(final_nc, "pr", start = c(1, 1, indices[1]), count = c(-1, -1, length(indices)), collapse_degen = FALSE)
      validate_era5_precipitation_values(source)
      compared <- compared + repair_equal(actual, source, paste("ERA5 monthly merge", item$month))
    }, finally = ncdf4::nc_close(n))
  }
  invisible(compared)
}

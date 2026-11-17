# Independent date labels for metadata audits only; never coerces 360-day dates
# (e.g. February 30) through Gregorian Date/POSIXct. Production selection unchanged.
#' Independently decode CF dates without Gregorian coercion of fixed calendars
#'
#' Gregorian inputs use UTC arithmetic. Fixed calendars use year-length and
#' month-offset arithmetic, including origin time-of-day and negative offsets.
#' A 1e-8-day boundary tolerance is applied before flooring fixed-calendar days.
#' This decoder labels records; it does not itself certify ordering, uniqueness,
#' complete years or the correctness of the climate-variable values.
#' @param values Finite numeric CF time offsets.
#' @param units Days/hours/seconds-since units string with a supported UTC origin.
#' @param calendar Calendar name or alias: standard/proleptic Gregorian, 360_day,
#'   365_day or 366_day. Julian and unknown calendars are not supported here.
#' @return List with ok=TRUE, integer year/month/day vectors, a source label and
#'   NA error_message. Invalid coordinates/origins or unsupported calendars error.
#' @examples
#' audit_cf_dates(59,"days since 2000-01-01","360_day") # February 30
audit_cf_dates <- function(values,units,calendar) {
  if(any(!is.finite(values)))stop("Non-finite time coordinate.")
  scale <- resolve_time_units_scale_days(units)
  calendar <- normalize_calendar_name(calendar)
  if(calendar%in%c("standard","gregorian","proleptic_gregorian")) {
    dates <- as.Date(parse_time_origin_datetime(units)+values*scale*86400,tz="UTC")
    return(list(ok=TRUE,year=as.integer(format(dates,"%Y")),
      month=as.integer(format(dates,"%m")),day=as.integer(format(dates,"%d")),
      source="independent_gregorian_UTC",error_message=NA_character_))
  }
  if(!calendar%in%c("360_day","365_day","366_day"))
    stop("Unsupported audit calendar: ",calendar)
  origin <- sub(".*since[[:space:]]+","",units,ignore.case=TRUE)
  pattern <- "^([0-9]+)-([0-9]{1,2})-([0-9]{1,2})(?:[ T]([0-9]{1,2}):([0-9]{1,2})(?::([0-9.]+))?)?(?:[ ]?UTC|Z)?$"
  fields <- regmatches(trimws(origin),regexec(pattern,trimws(origin),perl=TRUE))[[1]]
  if(!length(fields))stop("Unsupported fixed-calendar time origin: ",origin)
  parts <- suppressWarnings(as.numeric(fields[-1]));parts[is.na(parts)]<-0
  y<-parts[1];m<-parts[2];d<-parts[3]
  lengths <- if(calendar=="360_day")rep(30,12)else c(31,if(calendar=="366_day")29 else 28,31,30,31,30,31,31,30,31,30,31)
  if(m<1||m>12||d<1||d>lengths[m])stop("Invalid fixed-calendar origin.")
  cumulative <- c(0,cumsum(lengths));year_length <- sum(lengths)
  ordinal <- y*year_length+cumulative[m]+d-1+
    floor(values*scale+parts[4]/24+parts[5]/1440+parts[6]/86400+1e-8)
  years <- floor(ordinal/year_length);doy <- ordinal-years*year_length
  months <- findInterval(doy,cumulative)
  list(ok=TRUE,year=as.integer(years),month=as.integer(months),
    day=as.integer(doy-cumulative[months]+1),source="independent_fixed_calendar_arithmetic",
    error_message=NA_character_)
}

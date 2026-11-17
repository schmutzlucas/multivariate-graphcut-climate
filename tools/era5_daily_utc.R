# Read-only source checks / independent point validation for the CDO UTC builder.
source("R/load.R");load_mvgc();source("R/preprocessing/cf_calendar_validation.R")
source(Sys.getenv("MVGC_REPAIR_MODULE","R/preprocessing/netcdf_helpers.R"))
require_packages(c("ncdf4","jsonlite","openssl"))
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
job<-jsonlite::fromJSON(args[2]);mode<-args[1]
if(mode=="verify"&&file.exists(job$report))stop("Monthly verification report already exists")
#' Validate the native ERA5 hourly precipitation grid and decode its UTC stamps
#' @param n Open caller-owned ncdf4 handle containing tp in metres, accumulation
#'   step type, on the exact 0.25-degree longitude/descending-latitude/valid_time grid.
#' @return Numeric UTC seconds since the Unix epoch. Checks declared dimensions,
#'   units and Gregorian calendar; caller separately verifies complete hourly
#'   sequences/month boundaries. Reads coordinates only and does not close n.
hourly_axis<-function(n) {
 v<-n$var$tp;stopifnot(!is.null(v),length(v$dim)==3L,
  repair_attr(n,"tp","units")=="m",repair_attr(n,"tp","GRIB_stepType")=="accum")
 stopifnot(identical(vapply(v$dim,`[[`,"","name"),c("longitude","latitude","valid_time")),
  identical(as.numeric(v$dim[[1]]$vals),seq(0,359.75,by=.25)),
  identical(as.numeric(v$dim[[2]]$vals),seq(90,-90,by=-.25)))
 t<-v$dim[[3]];stopifnot(normalize_calendar_name(repair_attr(n,t$name,"calendar"))%in%
  c("standard","proleptic_gregorian"))
 as.numeric(parse_time_origin_datetime(t$units))+t$vals*resolve_time_units_scale_days(t$units)*86400
}
start<-as.Date(paste0(substr(job$month,1,4),"-",substr(job$month,5,6),"-01"))
next_start<-seq(start,by="month",length.out=2)[2];days<-as.integer(next_start-start)
expected<-as.numeric(as.POSIXct(start,tz="UTC"))+seq(0,days*24-1)*3600
if(mode=="source_identity") {
 # A resume reads all current source bytes in Python, but only the actually used
 # first map of the following month here. No unused following-month payload scan.
 n<-ncdf4::nc_open(job$raw);times<-hourly_axis(n)
 stopifnot(length(times)==length(expected),max(abs(times-expected))<1e-6)
 ncdf4::nc_close(n)
 n<-ncdf4::nc_open(job$following);times<-hourly_axis(n)
 stopifnot(times[1]==as.numeric(as.POSIXct(next_start,tz="UTC")))
 values<-ncdf4::ncvar_get(n,"tp",start=c(1,1,1),count=c(-1,-1,1),collapse_degen=FALSE)
 fingerprint<-era5_midnight_fingerprint(values,times[1]);ncdf4::nc_close(n)
 if(file.exists(job$fingerprint_report))stop("Source-identity report already exists")
 writeLines(jsonlite::toJSON(list(month=job$month,following_midnight_value_sha256=fingerprint),auto_unbox=TRUE),job$fingerprint_report)
} else if(mode=="boundary") {
 n<-ncdf4::nc_open(job$following);times<-hourly_axis(n)
 stopifnot(times[1]==as.numeric(as.POSIXct(next_start,tz="UTC")))
 out<-repair_new_file(n,"tp",job$boundary,n$var$tp$dim[[3]]$vals[1],
  n$var$tp$dim[[3]]$units,repair_attr(n,"valid_time","calendar"))
 values<-ncdf4::ncvar_get(n,"tp",start=c(1,1,1),count=c(-1,-1,1),collapse_degen=FALSE)
 stopifnot(all(is.finite(values)))
 ncdf4::ncvar_put(out,"tp",values)
 ncdf4::nc_close(out);ncdf4::nc_close(n)
 n<-ncdf4::nc_open(job$raw);times<-hourly_axis(n)
 stopifnot(length(times)==length(expected),max(abs(times-expected))<1e-6)
 ncdf4::nc_close(n);cat("PASS complete hourly source and next midnight",job$month,"\n")
} else if(mode=="verify") {
 n<-ncdf4::nc_open(job$raw);times<-hourly_axis(n)
 stopifnot(length(times)==length(expected),max(abs(times-expected))<1e-6)
 b<-ncdf4::nc_open(job$boundary)
 stopifnot(repair_attr(b,"tp","units")=="m",
  identical(as.numeric(b$var$tp$dim[[1]]$vals),seq(0,359.75,by=.25)),
  identical(as.numeric(b$var$tp$dim[[2]]$vals),seq(90,-90,by=-.25)))
 t<-b$var$tp$dim[[3]]
 boundary_time<-as.numeric(parse_time_origin_datetime(t$units))+t$vals*resolve_time_units_scale_days(t$units)*86400
 stopifnot(length(boundary_time)==1L,boundary_time==as.numeric(as.POSIXct(next_start,tz="UTC")))
 boundary_values<-ncdf4::ncvar_get(b,"tp",start=c(1,1,1),count=c(-1,-1,1),collapse_degen=FALSE)
 boundary_value_sha256<-era5_midnight_fingerprint(boundary_values,boundary_time)
 out<-ncdf4::nc_open(job$product);a<-repair_axis(out)
 stopifnot(identical(a$keys,format(seq(start,next_start-1,by="day"),"%Y-%m-%d")),
  repair_attr(out,"pr","units")=="mm day-1",
  identical(as.numeric(ncdf4::ncvar_get(out,"lon")),as.numeric(0:359)),
  identical(as.numeric(ncdf4::ncvar_get(out,"lat")),as.numeric(-90:90)))
 # Whole output scan; do not silently clip negatives or replace missing values.
 minimum<-Inf;maximum<--Inf
 for(d in seq_len(days)) {
  x<-ncdf4::ncvar_get(out,"pr",start=c(1,1,d),count=c(-1,-1,1))
  quality<-validate_era5_precipitation_values(x)
  minimum<-min(minimum,quality$minimum);maximum<-max(maximum,quality$maximum)
 }
 points<-list(c(6,46),c(120,0),c(250,-30));checks<-list()
 for(p in points) {
  ix<-as.integer(p[1]*4+1);iy<-as.integer((90-p[2])*4+1)
  h<-as.numeric(ncdf4::ncvar_get(n,"tp",start=c(ix,iy,1),count=c(1,1,-1)))
  boundary<-as.numeric(ncdf4::ncvar_get(b,"tp",start=c(ix,iy,1),count=c(1,1,1)))
  correct<-colSums(matrix(c(h[-1],boundary),nrow=24))*1000
  old<-colSums(matrix(h,nrow=24))*1000
  actual<-as.numeric(ncdf4::ncvar_get(out,"pr",start=c(p[1]+1,p[2]+91,1),count=c(1,1,-1)))
  tolerance<-1e-5+1e-6*abs(correct)
  stopifnot(all(is.finite(correct)),all(abs(actual-correct)<=tolerance))
  checks[[length(checks)+1L]]<-list(lon=p[1],lat=p[2],days=days,
   max_abs_vs_independent_hourly_sum=max(abs(actual-correct)),
   max_abs_change_from_old_window=max(abs(actual-old)))
 }
 ncdf4::nc_close(out);ncdf4::nc_close(n);ncdf4::nc_close(b)
 result<-list(state="utc_month_grid_and_independent_points_verified",month=job$month,days=days,
  source_sha256=job$source_sha256,following_source=job$following,
  boundary_sha256=repair_sha256(job$boundary),product_sha256=repair_sha256(job$product),
  product_bytes=file.info(job$product)$size,
  following_midnight_value_sha256=boundary_value_sha256,
  minimum_mm=minimum,maximum_mm=maximum,negative_values=0L,nonfinite_values=0L,
  independent_points=checks,validation_scope="All output cells scanned; three exact native/output grid points independently summed for every day. Not exhaustive source/output equality.")
 writeLines(jsonlite::toJSON(result,pretty=TRUE,auto_unbox=TRUE),job$report)
 cat("PASS verified UTC month",job$month,"\n")
} else stop("Unknown mode")

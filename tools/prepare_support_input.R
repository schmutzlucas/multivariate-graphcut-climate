# One exact prepared input: checksum, physical/time checks and scoped extrema.
source("R/load.R");load_mvgc();source("R/preprocessing/cf_calendar_validation.R")
source(Sys.getenv("MVGC_REPAIR_MODULE","R/preprocessing/netcdf_helpers.R"))
source("R/preprocessing/support_builder.R")
require_packages(c("ncdf4","jsonlite","openssl","matrixStats"))
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
job_sha256<-mvgc_sha256(args[1])
job<-jsonlite::fromJSON(args[1],simplifyVector=FALSE)
if(!identical(job$precipitation_preprocessing,article1_config()$precipitation_preprocessing))
 stop("A current nonnegative-precipitation QC job is required; old caches are incompatible")
if(!is.numeric(job$negative_review_threshold_mm_day)||length(job$negative_review_threshold_mm_day)!=1L||
 !is.finite(job$negative_review_threshold_mm_day)||job$negative_review_threshold_mm_day<=0)
 stop("Expected a documented positive precipitation review threshold")
if(!is.character(job$variable)||length(job$variable)!=1L||!job$variable%in%c("pr","tas","psl"))
 stop("Expected one prepared pr/tas/psl variable")
if(!is.list(job$scopes)||!length(job$scopes)||is.null(names(job$scopes))||
 anyNA(names(job$scopes))||anyDuplicated(names(job$scopes))||
 !all(names(job$scopes)%in%c("era5","pme")))stop("Expected explicit ERA5/PME year scopes")
for(scope in job$scopes) {
 if(any(lengths(scope)!=1L))stop("Missing or nonscalar support-year entry")
 years<-unlist(scope,use.names=FALSE)
 if(!is.numeric(years)||!length(years)||any(!is.finite(years))||any(years!=floor(years))||
  any(abs(years)>.Machine$integer.max)||anyDuplicated(years)||is.unsorted(years))
  stop("Support years must be sorted, distinct, finite integers; no implicit truncation")
}
for(coordinate in list(job$lon,job$lat)) {
 if(any(lengths(coordinate)!=1L))stop("Missing or nonscalar support coordinate")
 values<-unlist(coordinate,use.names=FALSE)
 if(!is.numeric(values)||!length(values)||any(!is.finite(values))||anyDuplicated(values))
  stop("Expected explicit unique finite support coordinates")
}
if(identical(normalizePath(job$output,winslash="/",mustWork=FALSE),
 normalizePath(job$report,winslash="/",mustWork=FALSE)))stop("Support cache/report paths must differ")
if(file.exists(job$output)||file.exists(job$report))stop("Use new support cache/report paths")
before<-file.info(job$file)
if(before$size!=job$bytes||repair_sha256(job$file)!=job$sha256)stop("Input SHA-256 mismatch")
n<-ncdf4::nc_open(job$file)
v<-job$variable;a<-repair_axis(n)
lon<-as.numeric(unlist(job$lon));lat<-as.numeric(unlist(job$lat))
stopifnot(identical(as.numeric(ncdf4::ncvar_get(n,"lon")),lon),
 identical(as.numeric(ncdf4::ncvar_get(n,"lat")),lat))
declared<-repair_attr(n,v,"units")
expected_units<-c(pr="mm day-1",tas="K",psl="Pa")[[v]]
if(v=="pr") {
 if(!declared%in%c("mm","mm day-1","mm/day")&&!isTRUE(job$unit_controls_passed))
  stop("Precipitation units need independently verified conversion controls")
} else stopifnot(declared==expected_units)
scopes<-job$scopes;wanted<-sort(unique(as.integer(unlist(scopes))))
selected<-which(a$dates$year%in%wanted)
year_counts<-vapply(wanted,function(y)sum(a$dates$year==y),0L)
stopifnot(identical(as.integer(year_counts),as.integer(days_in_year_by_calendar(wanted,a$calendar))))
step<-resolve_time_units_scale_days(a$units)
for(year in wanted) {
 ix<-which(a$dates$year==year)
 stopifnot(all(diff(ix)==1L),all(abs(diff(a$values[ix])*step-1)<1e-7),
  a$dates$month[ix[1]]==1L,a$dates$day[ix[1]]==1L,
  tail(a$dates$month[ix],1)==12L,
  tail(a$dates$day[ix],1)==days_in_month_by_calendar(year,12,a$calendar))
}
info<-get_variable_dimension_info(n,v,find_time_dimension_name(n))
coordinates<-match_requested_coordinates(n,info,lon,lat)
# Compare real production selection against independent CF date labels.
for(scope in names(scopes)) {
 years<-as.integer(unlist(scopes[[scope]]))
 selection<-select_full_year_indices(n,job$file,v,years,lon,lat,TRUE,"auto","auto","warn")
 stopifnot(identical(as.integer(selection$indices),as.integer(which(a$dates$year%in%years))))
}
ranges<-lapply(scopes,function(s)list(min=matrix(Inf,length(lon),length(lat)),
 max=matrix(-Inf,length(lon),length(lat))))
stats<-list(raw_min=Inf,raw_max=-Inf,negative_values=0,nonfinite_values=0,values_scanned=0)
qc<-NULL;qc_scopes<-setNames(vector("list",length(scopes)),names(scopes))
# Keep each chunk inside one CF year so every scope receives exact raw QC sums.
blocks<-do.call(rbind,lapply(wanted,function(year)
 split_blocks_to_max_count(split_contiguous_blocks(which(a$dates$year==year)),366L)))
for(k in seq_len(nrow(blocks))) {
 b<-blocks[k,];ix<-seq.int(b$start,length.out=b$count)
 read<-read_variable_time_chunk(n,v,info,coordinates$lon_indices,coordinates$lat_indices,b$start,b$count,report_qc=TRUE)
 cube<-read$data
 stats$values_scanned<-stats$values_scanned+length(cube)
 stats$nonfinite_values<-stats$nonfinite_values+sum(!is.finite(cube))
 if(any(is.infinite(cube)))stop("Infinite prepared climate values; see input ",job$file)
 if(any(is.finite(cube))) {
  stats$raw_min<-min(stats$raw_min,min(cube,na.rm=TRUE));stats$raw_max<-max(stats$raw_max,max(cube,na.rm=TRUE))
 }
 if(v=="pr") {
  qc<-merge_precipitation_qc_statistics(qc,read$qc)
  if(read$qc$maximum_absolute_correction>job$negative_review_threshold_mm_day) {
   writeLines(jsonlite::toJSON(list(model=job$model,file=job$file,scope_years=wanted,
    state="large_negative_precipitation_requires_investigation",statistics=qc),auto_unbox=TRUE,pretty=TRUE,null="null",na="null"),
    paste0(job$report,".qc_failure.json"))
   stop("Large negative precipitation requires source/units investigation; see QC failure report")
  }
 }
 for(scope in names(scopes)) {
  keep<-which(a$dates$year[ix]%in%as.integer(unlist(scopes[[scope]])))
  if(!length(keep))next
  if(v=="pr")qc_scopes[[scope]]<-merge_precipitation_qc_statistics(qc_scopes[[scope]],read$qc)
  s<-compute_cube_min_max(cube[,,keep,drop=FALSE],v)
  ranges[[scope]]$min<-pmin(ranges[[scope]]$min,s$min,na.rm=TRUE)
  ranges[[scope]]$max<-pmax(ranges[[scope]]$max,s$max,na.rm=TRUE)
 }
 rm(cube);if(k%%10L==0L){gc(FALSE);cat("Scanned blocks",k,"/",nrow(blocks),"\n");flush.console()}
}
ncdf4::nc_close(n)
after<-file.info(job$file)
stopifnot(before$size==after$size,before$mtime==after$mtime)
if(mvgc_sha256(args[1])!=job_sha256)stop("Support job changed during execution")
for(scope in names(ranges)) {
 for(bound in c("min","max"))ranges[[scope]][[bound]][!is.finite(ranges[[scope]][[bound]])]<-NA_real_
 s<-ranges[[scope]]
 stopifnot(identical(is.na(s$min),is.na(s$max)),all(s$min<=s$max,na.rm=TRUE))
}
if(v=="pr") {stats$raw_min<-qc$raw_minimum;stats$negative_values<-qc$negative_observations}
meta<-list(model=job$model,variable=v,file=job$file,sha256=job$sha256,bytes=job$bytes,
 job_sha256=job_sha256,preparation_identity=job$preparation_identity,
 software=list(R=as.character(getRversion()),ncdf4=as.character(packageVersion("ncdf4")),
  matrixStats=as.character(packageVersion("matrixStats"))),
 calendar=a$calendar,declared_units=declared,effective_units=expected_units,
 transform=if(v=="pr")"log1p(pmax(pr in mm day-1, 0)); preserve missing"else"identity",statistics=stats,
 precipitation_exception=NULL,precipitation_preprocessing=job$precipitation_preprocessing,
 precipitation_values_modified=if(v=="pr")qc$values_changed>0 else FALSE,
 qc_statistics=qc,qc_by_scope=if(v=="pr")qc_scopes else NULL,
 scopes=scopes,required_years=wanted,days_by_year=as.list(setNames(year_counts,wanted)),
 input_registry_sha256=job$input_registry_sha256,source_identity=job$source_identity,
 source_qualifications=job$qualifications)
saveRDS(list(ranges=ranges,meta=meta),job$output)
meta$state<-"payload_calendar_grid_units_and_scoped_ranges_verified"
meta$cache_sha256<-repair_sha256(job$output)
writeLines(jsonlite::toJSON(meta,auto_unbox=TRUE,pretty=TRUE,null="null",na="null",digits=NA),job$report)
cat("PASS prepared support input",job$model,v,"\n")

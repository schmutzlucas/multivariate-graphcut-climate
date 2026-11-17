source("R/load.R"); load_mvgc()
source("R/preprocessing/support_builder.R")

# Scalar broadcasting must agree with independent scalar calls for every calendar.
calendars <- c("standard", "proleptic_gregorian", "julian", "360_day", "365_day", "366_day")
for (calendar in calendars) {
  stopifnot(identical(as.integer(days_in_month_by_calendar(2000, 1:12, calendar)),
    vapply(1:12, function(m) as.integer(days_in_month_by_calendar(2000, m, calendar)), 0L)))
  years <- c(1900L, 1999L, 2000L, 2004L, 2100L)
  stopifnot(identical(as.integer(days_in_month_by_calendar(years, 2, calendar)),
    vapply(years, function(y) as.integer(days_in_month_by_calendar(y, 2, calendar)), 0L)))
  stopifnot(identical(as.integer(days_in_month_by_calendar(years, 1:5, calendar)),
    vapply(seq_along(years), function(i)
      as.integer(days_in_month_by_calendar(years[i], i, calendar)), 0L)))
}
stopifnot(identical(as.integer(days_in_month_by_calendar(1999:2001, 2, "standard")), c(28L,29L,28L)),
  identical(days_in_month_by_calendar(numeric(), 2, "standard"), integer()))
for (md in list(c(1.9,2.9), c(NA,1), c(1,Inf), c(2,29), c(4,31)))
  stopifnot(inherits(try(validate_month_day(md,"season"),silent=TRUE),"try-error"))
stopifnot(identical(validate_month_day(c(12,31),"season"), c(12L,31L)),
  month_day_to_reference_doy(c(12,31))==365L)
for (args in list(list(2000,1.5),list(NA_real_,2),list(2000,13),list(1:2,1:3)))
  stopifnot(inherits(try(do.call(days_in_month_by_calendar,c(args,list(calendar="standard"))),
    silent=TRUE),"try-error"))
stopifnot(!parse_filename_date_range("x_20001231-20000101.nc")$ok,
  parse_filename_date_range("x_20000230-20001230.nc")$ok)

# Bad chunk size previously made the while loop fail to advance forever.
blocks <- data.frame(start=c(2L,20L),count=c(7L,3L))
expected <- data.frame(start=c(2L,5L,8L,20L),count=c(3L,3L,1L,3L))
stopifnot(isTRUE(all.equal(split_blocks_to_max_count(blocks,3),expected)),
  identical(split_blocks_to_max_count(blocks[FALSE,],3),data.frame(start=integer(),count=integer())),
  nrow(split_blocks_to_max_count(data.frame(start=1L,count=0L),3))==0L)
for (n in list(0,-1,NA_real_,Inf,1.5,c(1,2)))
  stopifnot(inherits(try(split_blocks_to_max_count(blocks,n),silent=TRUE),"try-error"))
for (bad in list(data.frame(start=0L,count=1L),data.frame(start=1L,count=-1L),
                 data.frame(start=1.5,count=1L),data.frame(start=1L,count=NA_real_)))
  stopifnot(inherits(try(split_blocks_to_max_count(bad,3),silent=TRUE),"try-error"))

local({
  path <- tempfile(fileext=".nc")
  on.exit(unlink(path))
  td <- ncdf4::ncdim_def("time","days since 2000-01-01",as.numeric(0:2))
  nc <- ncdf4::nc_create(path,ncdf4::ncvar_def("dummy","1",list(td),-999))
  ncdf4::ncvar_put(nc,"time",c(0,NA_real_,2))
  ncdf4::nc_close(nc)
  nc <- ncdf4::nc_open(path); on.exit(ncdf4::nc_close(nc),add=TRUE)
  failure <- try(infer_time_structure_daily(nc,path,"time","standard",2000,
    time_decode_mode="empirical"),silent=TRUE)
  stopifnot(inherits(failure,"try-error"),grepl("nonfinite",failure,fixed=TRUE))
})
cat("PASS calendar scalar broadcasting, date contracts and bounded chunk splitting\n")

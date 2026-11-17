source("R/load.R");load_mvgc()
local({
  path <- tempfile(fileext=".nc");on.exit(unlink(path))
  # Deliberately store axes in a different order to the returned lon/lat/time cube.
  dims <- list(ncdf4::ncdim_def("time","days since 2000-01-01",0:4),
    ncdf4::ncdim_def("lat","degrees_north",c(-10,10)),
    ncdf4::ncdim_def("lon","degrees_east",0:2))
  n <- ncdf4::nc_create(path,ncdf4::ncvar_def("tas","K",dims,-999,prec="double"))
  values <- array(NA_real_,c(5,2,3))
  for(t in 1:5)for(j in 1:2)for(i in 1:3)values[t,j,i] <- 1000*i+100*j+t
  ncdf4::ncvar_put(n,"tas",values);ncdf4::nc_close(n)
  n <- ncdf4::nc_open(path);on.exit(ncdf4::nc_close(n),add=TRUE)
  info <- get_variable_dimension_info(n,"tas","time")
  selected <- c(5L,1L,3L,3L)
  got <- read_variable_time_blocks(n,"tas",info,c(3L,1L),c(2L,1L),selected)
  expected <- aperm(values[selected,c(2L,1L),c(3L,1L),drop=FALSE],c(3,2,1))
  stopifnot(identical(got,expected),all(is.finite(got)),
    identical(dim(read_variable_time_blocks(n,"tas",info,1L,1L,integer())),c(1L,1L,0L)))
  # Sorted unique production requests preserve the previous output exactly.
  stopifnot(identical(read_variable_time_blocks(n,"tas",info,1:3,1:2,c(1L,3L,4L)),
    aperm(values[c(1L,3L,4L),,,drop=FALSE],c(3,2,1))))
  for(bad in list(NA_real_,Inf,0,-1,1.5,6))
    stopifnot(inherits(try(read_variable_time_blocks(n,"tas",info,1L,1L,bad),silent=TRUE),"try-error"))
  for(bad in list(NA_real_,Inf,0,1.5,4,numeric()))
    stopifnot(inherits(try(read_variable_time_blocks(n,"tas",info,bad,1L,1L),silent=TRUE),"try-error"))
  plan <- list(file_path=path,dimension_info=info,lat_indices=1:2)
  for(bad in list(NA_real_,1.5,Inf))
    stopifnot(inherits(try(read_fractional_variable_plan_block(plan,bad),silent=TRUE),"try-error"))
})
stopifnot(identical(split_contiguous_blocks(c(5L,1L,3L,3L,4L)),
  data.frame(start=c(1L,3L),count=c(1L,3L))))
cat("PASS NetCDF read order, repeated times, permuted axes and invalid-index rejection\n")

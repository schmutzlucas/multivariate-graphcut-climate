source("R/load.R");load_mvgc()
local({
  variable <- "MVGC_AUDIT_INTEGER_SETTING"
  before <- Sys.getenv(variable,unset=NA_character_)
  on.exit(if(is.na(before))Sys.unsetenv(variable)else do.call(Sys.setenv,setNames(list(before),variable)))
  for(path in c("R/preprocessing/support_builder.R","R/experiments/perfect_model.R")) {
    e <- new.env(parent=globalenv());sys.source(path,e)
    Sys.unsetenv(variable)
    stopifnot(identical(e$read_positive_int_env(variable,4L),4L))
    for(value in c("1.5","Inf","NaN","NA","0","-2","2147483648","text")) {
      do.call(Sys.setenv,setNames(list(value),variable));warned <- FALSE
      result <- withCallingHandlers(e$read_positive_int_env(variable,4L),warning=function(w) {
        warned <<- TRUE;invokeRestart("muffleWarning")
      })
      stopifnot(warned,identical(result,4L))
    }
    for(value in c("1","6","1e2")) {
      do.call(Sys.setenv,setNames(list(value),variable))
      stopifnot(identical(e$read_positive_int_env(variable,4L),as.integer(as.numeric(value))))
    }
  }
  e <- new.env(parent=globalenv());sys.source("R/preprocessing/support_builder.R",e)
  do.call(Sys.setenv,setNames(list("Inf"),variable));warned <- FALSE
  result <- withCallingHandlers(e$read_positive_numeric_env(variable,12),warning=function(w) {
    warned <<- TRUE;invokeRestart("muffleWarning")
  })
  stopifnot(warned,identical(result,12))
})
cat("PASS finite nonfractional integer settings and finite positive memory settings\n")

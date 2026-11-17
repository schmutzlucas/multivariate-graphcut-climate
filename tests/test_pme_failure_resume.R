# Reduced controller integration: PMF/solver stubs, real checkpoint/result logic.
# This is not a native solver or full-climate numerical-equivalence claim.
source("R/load.R");load_mvgc()
local({
  root<-tempfile("pme_failure_resume_");dir.create(root)
  variables<-c("MVGC_INPUT_REGISTRY","PME_MAX_MODELS","PME_RUN_STAMP","PME_CHECKPOINT_DIR",
    "PME_OUTPUT_DIR","PME_PERIODS","PME_LAMBDAS","PME_FORCE_RECOMPUTE")
  previous<-Sys.getenv(variables,unset=NA_character_)
  on.exit(for(i in seq_along(variables))if(is.na(previous[i]))Sys.unsetenv(variables[i])else
    do.call(Sys.setenv,setNames(list(previous[i]),variables[i])))
  Sys.unsetenv(variables)
  Sys.setenv(PME_MAX_MODELS="3",PME_RUN_STAMP="test_resume")
  cfg<-article1_config();cfg$lon<-0:2;cfg$lat<-c(-30,30);cfg$nbins3d<-2L
  cfg$periods<-cfg$periods[1,,drop=FALSE];cfg$pme_lambdas<-c(.125,.25)
  cfg$data_root<-root;cfg$output_root<-root;cfg$pme_support<-file.path(root,"support.rds")
  saveRDS(list(synthetic=TRUE),cfg$pme_support)
  e<-new.env(parent=globalenv());load_mvgc(e);sys.source("R/experiments/perfect_model.R",e)
  e$article1_config<-function()cfg
  e$source<-function(file,...) {
    stopifnot(file%in%c("R/load.R","R/experiments/perfect_model.R"));invisible(NULL)
  }
  e$load_mvgc<-function(...)invisible(NULL)
  e$require_climate_data<-function(...)invisible(TRUE)
  e$require_production_solver<-function(...)invisible(TRUE)
  e$resolve_experiment_support<-function(...)cfg$pme_support
  e$build_range_cube<-function(...)array(rep(c(0,1),each=18),c(3,2,3,2))
  e$compute_nd_pdf_multi_seasonal<-function(...)list(
    pdf3=list(present=array(1/8,c(3,2,8,3)),future=array(1/8,c(3,2,8,3))),
    calendars=rep("365_day",3))
  counter<-new.env(parent=emptyenv());counter$calls<-0L;counter$fail_first<-TRUE
  e$GraphCutHellinger_nD_lat<-function(...) {
    counter$calls<-counter$calls+1L
    if(counter$fail_first&&counter$calls==1L)stop("deliberate synthetic solver failure")
    list(label_attribution=matrix(1L,3,2))
  }
  first<-try(suppressWarnings(sys.source("experiments/perfect_model/run.R",e)),silent=TRUE)
  stopifnot(inherits(first,"try-error"),counter$calls==6L,
    file.exists(file.path(e$output_dir,"test_resume_pme_incomplete_periods_compact.rds")),
    !file.exists(file.path(e$output_dir,"test_resume_pme_final_periods_compact.rds")))
  cp<-readRDS(e$checkpoint_paths[1])
  stopifnot(cp$meta$checkpoint_state=="incomplete_solver_failure",!anyDuplicated(names(cp$meta)))
  counter$fail_first<-FALSE
  suppressWarnings(sys.source("experiments/perfect_model/run.R",e))
  stopifnot(counter$calls==8L,file.exists(e$final_bundle_file),
    e$final_meta$run_state=="complete",e$final_meta$n_completed_period_references==3L,
    identical(e$partial_tau,cfg$pme_ldr_thresholds))
})
cat("PASS failed PME solver is not marked final; resume reuses only successful references and reruns failed cases\n")

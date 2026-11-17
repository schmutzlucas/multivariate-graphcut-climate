source("R/load.R");load_mvgc()
summary_env <- new.env(parent=globalenv());sys.source("R/metrics/pme_summaries.R",summary_env)
plot_env <- new.env(parent=globalenv());sys.source("R/plotting/pme_density.R",plot_env)
refs <- scan("config/pme_models.txt",what="",quiet=TRUE)
periods <- article1_config()$periods
result <- list(meta=list(model_names=refs),periods=periods,
  MMM_hdist_future_list=list(),GC_hdist_future_list=list(),
  MMM_partial_hdist_future_list=list(),GC_partial_hdist_future_list=list())
for (pid in periods$period_id) {
  maps <- setNames(lapply(seq_along(refs),function(i)matrix(seq_len(6)/100+i/1000,3,2)),refs)
  gc_maps <- lapply(maps,function(x).75*x)
  result$MMM_hdist_future_list[[pid]] <- maps
  result$GC_hdist_future_list[[pid]] <- list("0.125"=gc_maps[rev(refs)])
  result$MMM_partial_hdist_future_list[[pid]] <- list("0.01"=maps)
  result$GC_partial_hdist_future_list[[pid]] <- list("0.01"=list("0.125"=gc_maps))
}
weights <- summary_env$build_weight_matrix(0:2,c(-30,60))
for (type in c("full","partial")) {
  rows <- summary_env$build_metric_means_from_rds(result,type,type,"0.125","0.01",weights)
  stopifnot(nrow(rows)==75L,all(table(rows$period_id)==15L),
    identical(rows$reference,rep(refs,5)),max(abs(rows$rel_improvement-25))<1e-12)
}
for (e in list(summary_env,plot_env)) {
  stopifnot(identical(e$coerce_to_map_matrix(1:6,c(3L,2L),"vector"),matrix(as.numeric(1:6),3,2)),
    inherits(try(e$coerce_to_map_matrix(matrix(1:6,2,3),c(3L,2L),"transposed"),silent=TRUE),"try-error"),
    inherits(try(e$weighted_mean(1:3,1:2),silent=TRUE),"try-error"),
    inherits(try(e$build_weight_matrix(0:2,c(0,100)),silent=TRUE),"try-error"))
  bad <- result;bad$GC_hdist_future_list[[1]][[1]][[refs[1]]] <- NULL
  stopifnot(inherits(try(e$get_reference_order(bad,periods$period_id[1],"0.125",metric_type="full"),silent=TRUE),"try-error"))
}
bad <- result;bad$GC_hdist_future_list[[1]][[1]][[refs[1]]][1] <- NA_real_
stopifnot(inherits(try(summary_env$build_metric_means_from_rds(bad,"full","H","0.125",NULL,weights),silent=TRUE),"try-error"))
# Matching missing cells retain the old common-area averaging convention.
bad$MMM_hdist_future_list[[1]][[refs[1]]][1] <- NA_real_
matched <- summary_env$build_metric_means_from_rds(bad,"full","H","0.125",NULL,weights)
stopifnot(max(abs(matched$rel_improvement-25))<1e-12)
cat("PASS 15-reference summaries, exact identities, matrix orientation and paired valid-cell masks\n")

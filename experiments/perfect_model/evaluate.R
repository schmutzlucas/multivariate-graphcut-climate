# Evaluation reads saved production results; it does not refit labels.
# Missing/incompatible source objects fail explicitly. Outputs are derived scores
# under the configured output root and may overwrite existing evaluation files.
# Evaluate existing compact PME maps; never runs graph cut or builds climate PMFs.
source("R/load.R"); load_mvgc()
source("R/metrics/pme_summaries.R")
bundle <- readRDS(production_result_file("PME_FINAL_RESULTS_FILE","the compact final PME bundle"))
stopifnot(identical(as.integer(bundle$meta$calibration_years),1951:1975),
          identical(bundle$periods$period_id,article1_config()$periods$period_id),
          identical(bundle$meta$model_names,scan("config/pme_models.txt",what="",quiet=TRUE)))
weights <- build_weight_matrix(bundle$meta$lon,bundle$meta$lat)
rows <- lapply(c("full","0.005","0.01","0.025","0.05","0.10"),function(key) {
  z <- build_metric_means_from_rds(bundle,if(key=="full")"full"else"partial",
    if(key=="full")"Hellinger distance"else"LDR Hellinger score","0.125",
    if(key=="full")NULL else key,weights)
  data.frame(reference=z$reference,period_id=z$period_id,period=z$period,
    metric=key,ell_LDR=if(key=="full")NA_real_ else as.numeric(key),
    H_D_MMM=z$H_MMM,H_GC=z$H_GC,relative_improvement_pct=z$rel_improvement)
})
write.csv(do.call(rbind,rows),file.path(mvgc_output_dir("evaluation/perfect_model"),"reference_scores.csv"),row.names=FALSE)

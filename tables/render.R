# Level 1 tables, generated only from unrounded archived CSVs and model lists.
source("R/load.R"); load_mvgc()
source("R/utils/reproduction.R")
out <- mvgc_output_dir("tables")
era <- read.csv(reproduction_input("era5_scores"))
pme <- read.csv(reproduction_input("pme_scores"))
periods <- article1_config()$periods
#' Validate the archived numerical inputs required by Article 1 tables
#' @param era ERA5 long-statistics data frame, with unique scientific identities
#'   and finite numeric metric_value. Missing tau is allowed for ordinary metrics.
#' @param pme Reference-score data frame with reference, period_id, metric,
#'   H_D_MMM, H_GC and relative_improvement_pct columns.
#' @param periods Configured period data frame containing unique period_id values.
#' @param expected_references Exact configured held-out reference set, normally 15 models.
#' @return TRUE invisibly, or an error for invalid/duplicate data, inconsistent
#'   relative improvements or incomplete full-H/LDR-one-percent reference-period populations.
#'   No table is written and no numerical input is modified.
validate_table_inputs <- function(era, pme, periods, expected_references) {
  era_columns <- c("metric_family", "method", "dimension", "variable_set", "metric_type", "tau",
                   "metric_name", "metric_value", "units")
  pme_columns <- c("reference", "period_id", "metric", "H_D_MMM", "H_GC", "relative_improvement_pct")
  if (!is.data.frame(era) || !nrow(era) || !all(era_columns %in% names(era)) ||
      !is.numeric(era$metric_value) || any(!is.finite(era$metric_value)) ||
      anyNA(era[setdiff(era_columns, c("tau", "metric_value"))]) ||
      anyDuplicated(era[setdiff(era_columns, "metric_value")]))
    stop("ERA5 table inputs require finite values and complete unique scientific identities.")
  if (!is.data.frame(pme) || !nrow(pme) || !all(pme_columns %in% names(pme)) ||
      anyNA(pme[c("reference", "period_id", "metric")]) ||
      !all(vapply(pme[c("H_D_MMM", "H_GC", "relative_improvement_pct")], function(x)
        is.numeric(x) && all(is.finite(x)), logical(1))) ||
      any(pme$H_D_MMM <= 0) || any(pme$H_GC < 0) ||
      anyDuplicated(pme[c("reference", "period_id", "metric")]) ||
      !setequal(unique(pme$reference), expected_references) ||
      !setequal(unique(pme$period_id), periods$period_id))
    stop("PME table inputs require finite unique scores on the configured reference/period set.")
  expected_gain <- relative_improvement(pme$H_D_MMM, pme$H_GC)
  if (any(abs(expected_gain - pme$relative_improvement_pct) > 1e-10))
    stop("PME stored relative improvements disagree with their unrounded D-MMM/GC scores.")
  for (metric in c("full", "0.01")) for (pid in periods$period_id) {
    refs <- pme$reference[pme$metric == metric & pme$period_id == pid]
    if (length(refs) != length(expected_references) || !setequal(refs, expected_references))
      stop("Incomplete PME table population: ", metric, "/", pid)
  }
  invisible(TRUE)
}
expected_references <- scan(article1_config()$pme_models, what = "", quiet = TRUE)
validate_table_inputs(era, pme, periods, expected_references)

#' Look up one unrounded ERA5 statistic by its scientific identity
#' @param method,vars,metric Scalar internal method, variable-set and metric-type keys.
#' @param tau Scalar LDR threshold in (0,1), or NA_real_ for an ordinary metric.
#' @param name Scalar metric-name key; default weighted_mean.
#' @return One finite numeric value from caller-environment era. Missing or
#'   ambiguous rows raise an error; no first-row fallback or rounding occurs.
lookup <- function(method,vars,metric,tau=NA_real_,name="weighted_mean") {
  if (!all(vapply(list(method, vars, metric, name), function(x)
      is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x), logical(1))) ||
      !is.numeric(tau) || length(tau) != 1L ||
      (!is.na(tau) && (!is.finite(tau) || tau <= 0 || tau >= 1)))
    stop("Lookup requires scalar scientific keys and an LDR threshold or NA.")
  z <- era[era$method==method & era$variable_set==vars & era$metric_type==metric &
             era$metric_name==name & if(is.na(tau))is.na(era$tau)else !is.na(era$tau)&abs(era$tau-tau)<1e-12,]
  if(nrow(z)!=1L || !is.numeric(z$metric_value) || !is.finite(z$metric_value))
    stop("Ambiguous, absent or nonfinite ERA5 table input: ",method," / ",vars," / ",metric)
  z$metric_value
}
#' Write a manuscript-ready booktabs table from already formatted LaTeX fragments
#' @param filename Output filename below caller-environment out; existing generated
#'   files follow writeLines() overwrite behaviour.
#' @param label,caption LaTeX label and caption text, supplied without escaping.
#' @param header Character header row without the terminating line break.
#' @param rows Character vector of body rows without terminating line breaks.
#' @param columns LaTeX tabular column specification, e.g. lrrr.
#' @return NULL invisibly after writing a complete table environment. Caller must
#'   provide valid LaTeX and booktabs in the manuscript; no document compilation.
write_table <- function(filename,label,caption,header,rows,columns) {
  writeLines(c("% Generated from unrounded archived numerical inputs.","\\begin{table}[t]","\\centering","\\small",
    paste0("\\caption{",caption,"}"),paste0("\\label{",label,"}"),
    paste0("\\begin{tabular}{",columns,"}"),"\\toprule",paste0(header," \\\\"),"\\midrule",
    paste0(rows," \\\\"),"\\bottomrule","\\end{tabular}","\\end{table}"),file.path(out,filename))
}
#' Format paired ERA5 scores and their relative improvement as a LaTeX row
#' @param label Display label, already valid LaTeX.
#' @param vars,metric Variable-set and metric-type keys for lookup().
#' @param tau LDR threshold or NA_real_ for ordinary Hellinger/mean-field error.
#' @param name Metric-name key for lookup(), default weighted_mean.
#' @return Character row with baseline/GC scores to three decimals and percentage
#'   gain to one decimal. Improvement is computed from unrounded values first.
score_row <- function(label,vars,metric,tau=NA_real_,name="weighted_mean") {
  b <- lookup("MMM",vars,metric,tau,name);g <- lookup("GC-0.125",vars,metric,tau,name)
  sprintf("%s & %.3f & %.3f & %.1f\\,\\%%",label,b,g,relative_improvement(b,g))
}
rows <- c(score_row("Hellinger distance","pr_tas_psl","full_H"),
  vapply(c(.01,.05,.1),function(t)score_row(paste0("LDR score, $\\ell_{\\mathrm{LDR}}=",format(t),"$"),
    "pr_tas_psl","partial_LDR_H",t),""),
  vapply(c("tas","pr","psl"),function(v)score_row(paste0("Mean absolute error, ",v),v,
    "weighted_mean_absolute_bias",name="weighted_mean_absolute_bias"),""))
write_table("era5_global_scores.tex","tab:era5_global_scores",
  "ERA5 target-period scores. Distributional baselines are D-MMM; mean-field baselines are MMM. Positive relative improvement favours GC-0.125. Mean errors are in K (tas), mm/day (pr) and Pa (psl). Scores are rounded to three decimals and gains to one decimal, calculated from unrounded scores.",
  "Metric & Baseline & GC-0.125 & Improvement",rows,"lrrr")
sets <- c("pr_tas_psl","pr_tas","pr_psl","tas_psl","pr","tas","psl")
rows <- vapply(sets,function(vars) {
  full <- c(lookup("MMM",vars,"full_H"),lookup("GC-0.125",vars,"full_H"))
  ldr <- c(lookup("MMM",vars,"partial_LDR_H",.01),lookup("GC-0.125",vars,"partial_LDR_H",.01))
  label <- gsub("_",", ",vars,fixed=TRUE)
  sprintf("%s & %.3f & %.3f & %.1f & %.3f & %.3f & %.1f",label,full[1],full[2],relative_improvement(full[1],full[2]),
          ldr[1],ldr[2],relative_improvement(ldr[1],ldr[2]))
},"")
write_table("era5_variable_subsets.tex","tab:era5_distributional_sensitivity",
  "Independent 32-bin univariate, 16-bin bivariate and 8-bin trivariate diagnostics, retaining the trivariate graph-cut labels. LDR uses $\\ell_{\\mathrm{LDR}}=0.01$. Improvements are percentages relative to D-MMM, calculated before display rounding (scores: three decimals; gains: one decimal).",
  "Variables & H D-MMM & H GC & Gain & LDR D-MMM & LDR GC & Gain",rows,"lrrrrrr")
refs <- sort(expected_references)
for(key in c("full","0.01")) {
  z <- pme[pme$metric==key,]
  # Match by both identifiers explicitly; never depend on source row ordering.
  values <- sapply(periods$period_id,function(id){a<-z[z$period_id==id,];a$relative_improvement_pct[match(refs,a$reference)]})
  rows <- vapply(seq_along(refs),function(i)paste(c(refs[i],sprintf("%.1f",values[i,])),collapse=" & "),"")
  rows <- c(rows,paste(c("Mean",sprintf("%.1f",colMeans(values))),collapse=" & "),
    paste(c("Median",sprintf("%.1f",apply(values,2,median))),collapse=" & "),
    paste(c("IQR (Q1--Q3)",apply(values,2,function(v)paste(sprintf("%.1f",quantile(v,c(.25,.75),type=7)),collapse="--"))),collapse=" & "))
  write_table(paste0("pme_reference_",if(key=="full")"hellinger"else"ldr",".tex"),
    if(key=="full")"tab:pme_reference_improvement"else"tab:pme_reference_improvement_ldr",
    paste0("Reference-wise relative improvement of GC-0.125 over D-MMM for ",if(key=="full")"Hellinger distance"else"the LDR Hellinger score at $\\ell_{\\mathrm{LDR}}=0.01$",
      ". All 15 held-out references receive equal weight. Values are percentages."),
    paste(c("Held-out reference",gsub("-","--",periods$label,fixed=TRUE)),collapse=" & "),rows,"lrrrrr")
}
models <- scan("config/era5_models.txt",what="",quiet=TRUE)
pme_models <- scan("config/pme_models.txt",what="",quiet=TRUE)
metadata <- read.csv("data/metadata/models.csv",stringsAsFactors=FALSE)
stopifnot(identical(metadata$model,models),
          setequal(metadata$model[metadata$pme=="yes"],pme_models))
rows <- paste(metadata$model,metadata$institution,
              paste0("\\citet{",metadata$bibtex_key,"}"),
              ifelse(metadata$pme=="yes","X",ifelse(metadata$model=="ERA5","--","")),sep=" & ")
write_table("cmip6_models.tex","tab:cmip6_models","Reference dataset and CMIP6 models used in the experiments. All listed CMIP6 models are candidates in the ERA5 experiment; PME marks the duplicate-filtered subset. References are in data/metadata/model_references.bib.",
  "Dataset or model & Institution & Reference & PME",rows,"llll")
write.csv(data.frame(reference=refs),file.path(out,"reference_order.csv"),row.names=FALSE)
cat("PASS: generated five manuscript table components from archived unrounded inputs.\n")

# Level 1 orchestration. Default: completed, corrected article cohort.
scripts <- c("figures/era5_distributional.R","figures/era5_mean_bias.R",
             "figures/pme_overlays.R","figures/pme_reference_improvements.R",
             "figures/spatial_regularisation.R","figures/histogram_resolution.R","figures/lambda_tradeoff.R","tables/render.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) {
  requested <- strsplit(args[[1]],",",fixed=TRUE)[[1]]
  known <- tools::file_path_sans_ext(basename(scripts))
  if (!length(requested) || any(!requested %in% known)) stop("Unknown reproduction entry point.")
  scripts <- scripts[tools::file_path_sans_ext(basename(scripts)) %in% requested]
  if(!length(scripts)) stop("No matching reproduction entry point.")
}
source("R/load.R");load_mvgc()
source("R/utils/reproduction.R")
profile <- reproduction_profile()
# An independent output root preserves the frozen publication assets.
# Explicit input/output overrides are recorded in the external execution log.
if (!nzchar(Sys.getenv("MVGC_OUTPUT_ROOT", "")))
  Sys.setenv(MVGC_OUTPUT_ROOT = file.path("outputs/reproduction", profile$profile))
inventory <- reproduction_inventory(scripts, profile)
private <- Sys.getenv("MVGC_PRIVATE_VERIFICATION_DIR", tempfile("mvgc-render-"))
dir.create(private, recursive = TRUE, showWarnings = FALSE)
cat("Reproduction profile:", profile$profile, "\n")
if (any(inventory$explicit_override))
  message("Explicit input overrides are active; see reproduction_inputs.csv for exact identities.")
write.csv(inventory, file.path(private, "reproduction_inputs.csv"), row.names = FALSE)
write.csv(as.data.frame(installed.packages()[,c("Package","Version"),drop=FALSE]),
          file.path(private,"software_versions.csv"),row.names=FALSE)
rs <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
log_dir <- file.path(private,"render_logs")
dir.create(log_dir,recursive=TRUE,showWarnings=FALSE)
rows <- list()
for(script in scripts) {
  cat("Reproducing",script,"\n")
  log <- file.path(log_dir,paste0(tools::file_path_sans_ext(basename(script)),".txt"))
  start <- proc.time()[["elapsed"]]
  code <- system2(rs,shQuote(script),stdout=log,stderr=log)
  rows[[length(rows)+1L]] <- data.frame(script=script,exit_status=code,
    elapsed_seconds=proc.time()[["elapsed"]]-start,log=log)
  if(code!=0L) {cat(readLines(log,warn=FALSE),sep="\n");stop("Reproduction failed: ",script)}
}
write.csv(do.call(rbind,rows),file.path(private,"render_execution.csv"),row.names=FALSE)
outputs <- list.files(article1_config()$output_root,recursive=TRUE,full.names=TRUE)
outputs <- outputs[file.exists(outputs)&!dir.exists(outputs)]
write.csv(data.frame(file=outputs,bytes=file.info(outputs)$size,
  sha256=vapply(outputs,function(p)digest::digest(file=p,algo="sha256"),"")),
  file.path(private,"render_outputs.csv"),row.names=FALSE)
cat("PASS: all selected figure/table entry points completed.\n")

if (length(scripts) == 8L) {
  figures <- read.csv("publication/figure_manifest.csv")
  tables <- read.csv("publication/table_manifest.csv")
  expected <- c(paste0(rep(figures$renderer_stem, each=2), c(".pdf", ".png")),
                basename(tables$accepted_latex))
  actual <- basename(outputs)
  if (length(expected) != 39L || any(vapply(expected,function(x) sum(actual==x),integer(1)) != 1L))
    stop("Expected exactly one of each of 17 figure pairs and five table components.")
}

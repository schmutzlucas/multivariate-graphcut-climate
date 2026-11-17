#!/usr/bin/env Rscript

# Regenerate the article-style bin-sensitivity figure from the completed result
# without rereading the climate archive.

source("R/load.R")
load_mvgc()
source("R/utils/reproduction.R")
require_packages("ggplot2")

result_path <- reproduction_input("histogram")
if (!file.exists(result_path)) stop("Sensitivity result not found: ", result_path)
result <- readRDS(result_path)
write_k_sensitivity_outputs(result, mvgc_output_dir("figures/histogram_resolution"))

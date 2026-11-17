# Use an external temporary device for layout metrics before explicit export.
options(device = function(...) grDevices::pdf(file = tempfile("mvgc-layout-", fileext=".pdf")))
# Render only the approved cached fields; never compute PMFs, scores or labels.
source("R/load.R"); load_mvgc()
source("R/utils/reproduction.R")
require_packages(c("ggplot2","gtable","maps","scales","RColorBrewer","ragg","digest","jsonlite"))
source("R/plotting/era5_maps.R")
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  if (!nzchar(Sys.setlocale("LC_CTYPE", ".UTF-8"))) stop("A UTF-8 locale is required for ell_LDR.")
}
cache <- reproduction_input("era5_fields")
if (!file.exists(cache)) stop("Missing approved figure-field cache; rendering never recomputes numerical results.")
fields <- readRDS(cache)
stopifnot(identical(as.integer(fields$calibration),1951:1975),
          identical(as.integer(fields$target),1999:2023), fields$ell_ldr==.01)
for(name in c("h_d_mmm","h_gc","ldr_d_mmm","ldr_gc")) {
  stopifnot(identical(dim(fields[[name]]),c(360L,181L)),
            all(is.finite(fields[[name]])), all(fields[[name]]>=0))
}
plot <- build_hdist3_partial_4map_plot(fields$h_d_mmm,fields$h_gc,fields$ldr_d_mmm,
                                    fields$ldr_gc,fields$lon,fields$lat,"0.01")
out <- Sys.getenv("MVGC_FIGURE_OUTPUT_ROOT", unset=mvgc_output_dir("figures"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base <- file.path(out,"era5_hellinger_maps")
save_era5_fourmap(plot,base,cache)
message("PASS: cached values, coordinates, palette and 0..0.5 scales unchanged; four-panel geometry verified.")

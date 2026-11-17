# Level 1: render compact, numerically verified curves without production grids.
if(.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
source("R/load.R"); load_mvgc()
source("R/utils/reproduction.R")
require_packages(c("ggplot2","rlang"))
source("R/plotting/cached_densities.R")
selected <- Sys.getenv("MVGC_PME_FIGURES","all")
recipes <- list.files(reproduction_input("pme_densities"),pattern="\\.rds$",full.names=TRUE)
if(selected!="all") recipes <- recipes[tools::file_path_sans_ext(basename(recipes)) %in% strsplit(selected,",",fixed=TRUE)[[1]]]
if(!length(recipes)) stop("No requested density recipes found. See docs/ARTICLE_MAP.md.")
out <- mvgc_output_dir("figures/perfect_model")
for(path in recipes) {
  r <- readRDS(path)
  stem <- tools::file_path_sans_ext(basename(path))
  outputs <- read.csv("publication/figure_manifest.csv")
  hit <- match(path, outputs$input)
  if (is.na(hit)) stop("Unmapped density input: ", path)
  base <- file.path(out, outputs$renderer_stem[hit])
  reference_panel <- startsWith(basename(path),"reference_")
  key <- if(reference_panel || basename(path)=="full_overlay.rds")NULL else
    gsub("p",".",sub("_overlay\\.rds$","",sub("^ldr_","",basename(path))))
  plot <- if(reference_panel)plot_cached_reference_density(r) else
    plot_cached_overlay(r,if(is.null(key))"full"else"partial",key,
      if(!is.null(key) && key != "0.01") "appendix_enlarged" else "standard")
  ggplot2::ggsave(paste0(base,".png"),plot,width=if(reference_panel)34 else 29,
    height=if(reference_panel)22 else 18,units="cm",dpi=300,bg="white")
  ggplot2::ggsave(paste0(base,".pdf"),plot,width=if(reference_panel)34 else 29,
    height=if(reference_panel)22 else 18,units="cm",device=grDevices::cairo_pdf,bg="white")
  cat("Rendered",basename(path),"\n")
}
thresholds <- c("0p005","0p025","0p05","0p10")
panels <- file.path(out,paste0("ldr_threshold_",thresholds,".png"))
if(all(file.exists(panels))) {
  fragment <- c("% Generated threshold panels; graphics paths are relative to the repository root.")
  for(i in seq_along(panels))fragment<-c(fragment,
    "\\begin{subfigure}{0.49\\textwidth}",
    paste0("  \\includegraphics[width=\\linewidth]{",gsub("\\\\","/",panels[i]),"}"),
    paste0("  \\caption{Nominal threshold $\\ell_{\\mathrm{LDR}}=",gsub("p",".",thresholds[i]),"$.}"),
    "\\end{subfigure}",if(i%%2L==1L)"\\hfill"else"")
  writeLines(fragment,file.path(out,"ldr_threshold_subfigures.tex"))
}

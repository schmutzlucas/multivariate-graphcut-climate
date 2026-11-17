source("R/load.R"); load_mvgc()
# The bundled supports are path-redacted copies of the production supports.
# Check identity, array layout and the exact model/year scope used by the article.
manifest <- read.csv("data/supports/manifest.csv", stringsAsFactors = FALSE)
stopifnot(identical(sort(manifest$file), c("era5_support.rds", "pme_support.rds")))
cfg <- article1_config()
for (i in seq_len(nrow(manifest))) {
  path <- file.path("data/supports", manifest$file[i])
  stopifnot(file.size(path) == manifest$bytes[i], identical(mvgc_sha256(path), manifest$sha256[i]))
  mode <- sub("_support.rds", "", manifest$file[i], fixed = TRUE)
  x <- readRDS(path)
  scope <- support_scope(mode, cfg)
  stopifnot(identical(names(x$ranges), c("pr", "tas", "psl")),
            all(vapply(x$ranges, function(a) identical(dim(a), c(360L, 181L, 2L)), logical(1))),
            all(vapply(x$ranges, function(a) all(is.finite(a)) && all(a[, , 2] >= a[, , 1]), logical(1))),
            identical(as.character(x$meta$model_names), as.character(scope$model_names)),
            identical(as.numeric(x$meta$years), as.numeric(scope$years)),
            identical(x$meta$precipitation_preprocessing, "nonnegative_then_log1p_v1"),
            identical(x$meta$input_registry_sha256,
                      "c34724582dfbe7ed2f1ef73886bf83f88d758482e8f81922f70d6ca7a2a653b6"),
            length(x$meta$input_records) == 3L * length(scope$model_names),
            !any(grepl("[/\\\\]", vapply(x$meta$input_records, `[[`, "", "file"))))
  cat(sprintf("PASS: %s support, %d datasets, %d years, arrays 360x181x2 for pr/tas/psl.\n",
              mode, length(scope$model_names), length(scope$years)))
}

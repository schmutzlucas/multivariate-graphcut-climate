source("R/load.R"); load_mvgc(); require_packages(c("jsonlite", "openssl"))
source("R/preprocessing/support_builder.R")
local({
  e <- new.env(parent = globalenv())
  for (expr in parse("tools/aggregate_prepared_supports.R"))
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("aggregate_prepared_supports"))) eval(expr, e)
  root <- tempfile("support_aggregation_"); dir.create(root)
  models <- c("REF", "A", "B", "C"); variables <- c("pr", "tas", "psl")
  scopes <- list(era5 = list(models = models, years = c(1951L, 1999L)),
    pme = list(models = models[-1], years = 1951:1952))
  identity <- paste(rep("a", 64), collapse = ""); registry <- paste(rep("b", 64), collapse = "")
  records <- list()
  for (m in seq_along(models)) for (v in seq_along(variables)) {
    bounds <- list(min = matrix(10 * v + m, 360, 181), max = matrix(10 * v + m + 1, 360, 181))
    cache <- file.path(root, paste0(models[m], "_", variables[v], ".rds"))
    meta <- list(model = models[m], variable = variables[v], job_sha256 = identity,
      precipitation_preprocessing=article1_config()$precipitation_preprocessing,
      preparation_identity = identity, input_registry_sha256 = registry,
      scopes = lapply(scopes, `[[`, "years"))
    saveRDS(list(ranges = list(era5 = bounds, pme = bounds), meta = meta), cache)
    records[[length(records) + 1L]] <- list(model = models[m], variable = variables[v], cache = cache,
      cache_sha256 = mvgc_sha256(cache), job_sha256 = identity, state = "payload_calendar_grid_units_and_scoped_ranges_verified")
  }
  archive <- file.path(root, "archive.rds")
  archived <- setNames(lapply(seq_along(variables), function(v) array(c(rep(10 * v, 65160), rep(10 * v + 6, 65160)), c(360, 181, 2))), variables)
  saveRDS(list(ranges = archived), archive)
  make_job <- function(case) list(scopes = scopes, records = records, preparation_identity = identity,
    input_registry_sha256 = registry, outputs = list(era5 = file.path(root, case, "era5.rds"), pme = file.path(root, case, "pme.rds")),
    archived = list(era5 = archive, pme = archive), report = file.path(root, case, "validation.json"))
  job <- make_job("valid"); manifest <- e$aggregate_prepared_supports(job)
  for (mode in names(scopes)) {
    x <- readRDS(job$outputs[[mode]])
    for (v in seq_along(variables)) {
      stopifnot(all(x$ranges[[v]][, , 1] == 10 * v + if (mode == "era5") 1 else 2),
        all(x$ranges[[v]][, , 2] == 10 * v + 5))
    }
    stopifnot(identical(x$meta$years, scopes[[mode]]$years), identical(x$meta$model_names, scopes[[mode]]$models),
      x$meta$discretisation$trivariate == 8L, manifest[[mode]]$sha256 == mvgc_sha256(job$outputs[[mode]]),
      length(x$meta$input_records) == length(scopes[[mode]]$models) * 3L)
  }
  original <- vapply(unlist(job$outputs), mvgc_sha256, "")
  stopifnot(inherits(try(e$aggregate_prepared_supports(job), silent = TRUE), "try-error"),
    identical(original, vapply(unlist(job$outputs), mvgc_sha256, "")))
  for (value in c(0,log1p(3))) {
    constant <- make_job(paste0("point_mass_",value))
    for(index in seq(1L,length(records),by=3L)) {
      x<-readRDS(records[[index]]$cache)
      for(mode in names(scopes))x$ranges[[mode]]$min[,]<-x$ranges[[mode]]$max[,]<-value
      cache<-file.path(root,paste0("point_mass_",value,"_",index,".rds"));saveRDS(x,cache)
      constant$records[[index]]$cache<-cache;constant$records[[index]]$cache_sha256<-mvgc_sha256(cache)
    }
    e$aggregate_prepared_supports(constant)
    x<-readRDS(constant$outputs$pme)
    stopifnot(all(x$ranges$pr[,,1]==value),all(x$ranges$pr[,,2]>value),
      nrow(x$meta$constant_cells$pr)==65160)
  }
  # Failed PME validation must not leave an apparently successful ERA5 support.
  for (kind in c("fractional_year", "scope", "duplicate", "identity", "shape", "missing", "archive")) {
    bad <- make_job(kind)
    if (kind == "fractional_year") bad$scopes$pme$years[1] <- 1951.5
    if (kind == "scope") bad$scopes$pme$years <- 1952:1953
    if (kind == "duplicate") bad$records <- c(bad$records, bad$records[1])
    if (kind %in% c("identity", "shape")) {
      # Change only PME content in one cache, then refresh its byte checksum.
      index <- 4L; source <- readRDS(records[[index]]$cache)
      if (kind == "identity") source$meta$model <- "wrong"
      if (kind == "shape") source$ranges$pme$min <- matrix(1, 1, 181)
      cache <- file.path(root, paste0(kind, ".rds")); saveRDS(source, cache)
      bad$records[[index]]$cache <- cache; bad$records[[index]]$cache_sha256 <- mvgc_sha256(cache)
    }
    if (kind == "missing") for (index in seq(1L, length(records), by = 3L)) {
      source <- readRDS(records[[index]]$cache)
      source$ranges$pme$min[,] <- source$ranges$pme$max[,] <- NA_real_
      cache <- file.path(root, paste0("constant_", index, ".rds")); saveRDS(source, cache)
      bad$records[[index]]$cache <- cache; bad$records[[index]]$cache_sha256 <- mvgc_sha256(cache)
    }
    if (kind == "archive") { p <- file.path(root, "bad_archive.rds"); saveRDS(list(ranges = list()), p); bad$archived$pme <- p }
    stopifnot(inherits(try(e$aggregate_prepared_supports(bad), silent = TRUE), "try-error"),
      !any(file.exists(c(unlist(bad$outputs), bad$report))))
  }
  # A publish failure rolls back only files this invocation created, and keeps
  # all source caches. It must not leave the first support looking approved.
  counter <- 0L
  e$file.rename <- function(from, to) { counter <<- counter + 1L; if (counter == 2L) FALSE else base::file.rename(from, to) }
  failed <- make_job("publish_failure")
  stopifnot(inherits(try(e$aggregate_prepared_supports(failed), silent = TRUE), "try-error"),
    !any(file.exists(c(unlist(failed$outputs), failed$report))),
    !length(list.files(dirname(failed$report))), all(file.exists(vapply(records, `[[`, "", "cache"))))
})
cat("PASS support aggregation: exact population/extrema, complete provenance, all-scopes validation, identity/dimension guards and no overwrite\n")

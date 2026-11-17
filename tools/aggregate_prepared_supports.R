source("R/load.R"); load_mvgc(); require_packages(c("jsonlite", "openssl"))
source("R/preprocessing/support_builder.R")

#' Pool checked per-input extrema into separate ERA5 and common PME supports
#' @param job Parsed aggregation job with era5/pme scopes, signed cache records,
#'   input_registry_sha256, preparation_identity, new output paths, archived
#'   support paths and a new JSON report path. Bounds use the fixed 360 x 181 grid.
#' @return Named era5/pme manifest invisibly. Each support contains exactly its
#'   ordered model/year scope, transforms, numerical settings and input metadata.
#'   Validates both complete supports before publishing any output; refuses
#'   overwrite. Checks cache model/variable/job/source identities and dimensions,
#'   not just their filenames. No raw inputs are changed or experiments started.
#'   Missing pooled cells fail with explicit locations. Constant empirical bounds
#'   preserve their minimum and use the documented upper-only numerical envelope.
aggregate_prepared_supports <- function(job) {
  modes <- c("era5", "pme"); variables <- c("pr", "tas", "psl")
  if (!is.list(job) || !setequal(names(job$scopes), modes) ||
      !is.character(job$preparation_identity) || length(job$preparation_identity) != 1L ||
      is.na(job$preparation_identity) || !grepl("^[[:xdigit:]]{64}$", job$preparation_identity) ||
      !is.character(job$input_registry_sha256) || length(job$input_registry_sha256) != 1L ||
      is.na(job$input_registry_sha256) || !grepl("^[[:xdigit:]]{64}$", job$input_registry_sha256))
    stop("Expected two explicit support scopes and signed preparation/registry identities.")
  destinations <- c(unlist(job$outputs[modes], use.names = FALSE), job$report)
  if (length(destinations) != 3L || !is.character(destinations) || anyNA(destinations) ||
      any(!nzchar(destinations)) || anyDuplicated(normalizePath(destinations, winslash = "/", mustWork = FALSE)) ||
      any(file.exists(destinations))) stop("Choose distinct new support and report output paths.")
  # Detect accidental output-to-input aliases even if a caller bypassed Python.
  inputs <- c(vapply(job$records, function(x) x$cache, ""), unlist(job$archived[modes], use.names = FALSE))
  if (any(normalizePath(destinations, winslash = "/", mustWork = FALSE) %in%
          normalizePath(inputs, winslash = "/", mustWork = TRUE))) stop("Support output aliases an input.")
  objects <- manifest <- list()
  for (mode in modes) {
    scope <- job$scopes[[mode]]; models <- unlist(scope$models, use.names = FALSE); years <- unlist(scope$years, use.names = FALSE)
    if (any(lengths(scope$models) != 1L) || any(lengths(scope$years) != 1L) ||
        !is.character(models) || !length(models) || anyNA(models) || any(!nzchar(models)) || anyDuplicated(models) ||
        !is.numeric(years) || !length(years) || any(!is.finite(years)) || any(years != floor(years)) ||
        any(abs(years) > .Machine$integer.max) || anyDuplicated(years) || is.unsorted(years))
      stop("Expected distinct models and sorted, finite integer years for ", mode, ".")
    years <- as.integer(years)
    ranges <- setNames(vector("list", 3L), variables); records <- list(); constant_cells <- list()
    for (v in variables) {
      low <- matrix(Inf, 360, 181); high <- matrix(-Inf, 360, 181)
      for (model in models) {
        item <- Filter(function(r) identical(r$model, model) && identical(r$variable, v), job$records)
        if (length(item) != 1L) stop("Expected exactly one cache for ", model, "/", v, ".")
        item <- item[[1]]
        if (!identical(item$state, "payload_calendar_grid_units_and_scoped_ranges_verified") ||
            !identical(mvgc_sha256(item$cache), item$cache_sha256)) stop("Unverified or changed support cache: ", item$cache)
        x <- readRDS(item$cache)
        if (!identical(x$meta$model, model) || !identical(x$meta$variable, v) ||
            !identical(x$meta$job_sha256, item$job_sha256) || is.null(item$job_sha256) ||
            !identical(x$meta$preparation_identity, job$preparation_identity) ||
            !identical(x$meta$precipitation_preprocessing, article1_config()$precipitation_preprocessing) ||
            !identical(x$meta$input_registry_sha256, job$input_registry_sha256) ||
            !identical(as.numeric(unlist(x$meta$scopes[[mode]])), as.numeric(years)))
          stop("Support cache source/job/scope identity mismatch: ", item$cache)
        bounds <- x$ranges[[mode]]
        if (!is.numeric(bounds$min) || !is.numeric(bounds$max) ||
            !identical(dim(bounds$min), c(360L, 181L)) || !identical(dim(bounds$max), c(360L, 181L)) ||
            any(is.infinite(bounds$min)) || any(is.infinite(bounds$max)) ||
            !identical(is.na(bounds$min),is.na(bounds$max)) || any(bounds$min > bounds$max,na.rm=TRUE))
          stop("Invalid per-input support extrema: ", item$cache)
        low <- pmin(low, bounds$min,na.rm=TRUE); high <- pmax(high, bounds$max,na.rm=TRUE)
        records[[paste(model, v, sep = "/")]] <- x$meta
      }
      absent <- which(!is.finite(low)|!is.finite(high),arr.ind=TRUE)
      if(nrow(absent)) {
        dir.create(dirname(job$report),recursive=TRUE,showWarnings=FALSE)
        utils::write.csv(data.frame(scope=mode,variable=v,lon=absent[,1]-1L,lat=absent[,2]-91L),
          paste0(job$report,".missing_cells.csv"),row.names=FALSE)
        stop("No valid pooled observations; locations reported. No support invented for ",mode,"/",v)
      }
      adjusted <- regularize_empirical_support(low,high,v)
      cells <- which(adjusted$constant,arr.ind=TRUE)
      constant_cells[[v]] <- data.frame(lon=cells[,1]-1L,lat=cells[,2]-91L,
        empirical_value=low[cells],upper_expansion=adjusted$expansion[cells])
      ranges[[v]] <- array(c(adjusted$min, adjusted$max), c(360, 181, 2))
    }
    objects[[mode]] <- list(ranges = ranges, meta = list(script = "tools/aggregate_prepared_supports.R", mode = mode,
      created_at = as.character(Sys.time()), years = years, model_names = models, variables = variables,
      lon = 0:359, lat = -90:90, season_start_md = c(1, 1), season_end_md = c(12, 31),
      support_definition = "Per-cell min/max pooled over exactly the listed models and years; one common PME support for all held-out references and periods",
      precipitation_preprocessing = article1_config()$precipitation_preprocessing,
      empirical_minimum_policy = "Observed min(log1p(pmax(pr,0))); strictly positive minima are retained",
      constant_cells = constant_cells,
      constant_policy = "Preserve empirical minimum; upper-only envelope max(1e-12,64*epsilon*max(1,abs(min))) for point masses",
      discretisation = list(trivariate = 8L, bivariate = 16L, univariate = 32L),
      preparation_identity = job$preparation_identity, input_registry_sha256 = job$input_registry_sha256, input_records = records))
    manifest[[mode]] <- list(path=job$outputs[[mode]],models=models,years=years)
    if (!is.null(job$archived[[mode]])) {
    old <- readRDS(job$archived[[mode]])
    comparison <- lapply(variables, function(v) {
      if (!is.numeric(old$ranges[[v]]) || !identical(dim(old$ranges[[v]]), c(360L, 181L, 2L)) ||
          any(!is.finite(old$ranges[[v]]))) stop("Invalid archived comparison support: ", mode, "/", v)
      d <- ranges[[v]] - old$ranges[[v]]
      list(variable = v, changed_bounds = sum(d != 0), total_bounds = length(d), max_abs_bound_change = max(abs(d)),
        interpretation = "Input and scope changed intentionally; not a numerical-equivalence test")
    })
    manifest[[mode]] <- list(path = job$outputs[[mode]], models = models, years = years, comparison_with_archive = comparison)
    }
  }
  # No output exists until all populations, cache identities, bound dimensions
  # and archived comparisons have passed. Stage each complete payload locally.
  staged <- installed <- character(); success <- FALSE
  on.exit({unlink(staged); if (!success) unlink(installed)}, add = TRUE)
  for (mode in modes) {
    path <- job$outputs[[mode]]; dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("support_", tmpdir = dirname(path), fileext = ".rds"); staged <- c(staged, temporary)
    saveRDS(objects[[mode]], temporary, compress = "gzip")
    manifest[[mode]]$sha256 <- mvgc_sha256(temporary); manifest[[mode]]$bytes <- file.info(temporary)$size
  }
  dir.create(dirname(job$report), recursive = TRUE, showWarnings = FALSE)
  report <- tempfile("support_report_", tmpdir = dirname(job$report), fileext = ".json"); staged <- c(staged, report)
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), report)
  for (i in seq_along(destinations)) {
    if (file.exists(destinations[i]) || !file.rename(staged[i], destinations[i])) stop("Could not publish support output: ", destinations[i])
    installed <- c(installed, destinations[i])
  }
  success <- TRUE
  invisible(manifest)
}

args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Supply exactly one aggregation job JSON.")
aggregate_prepared_supports(jsonlite::fromJSON(args[1], simplifyVector = FALSE))
cat("PASS two separately scoped common supports\n")

# Sensitivity entry point: changing the tested resolution or regularisation
# changes the scientific configuration. Keep external inputs/supports fixed when
# comparing settings. Resource controls affect execution only; retain distinct
# output/checkpoint namespaces to avoid mixing incompatible experiments.
#!/usr/bin/env Rscript
source("R/load.R")
load_mvgc()
source("R/experiments/perfect_model.R") # Shared, name-aware support loader; definitions only.
cfg <- article1_config()

# Calibration-only trivariate bin-count stability diagnostic.
# This script does not run graph cut, future periods, or LDR and does not alter
# the production setting k = 8.

if (!requireNamespace("ncdf4", quietly = TRUE) ||
    !requireNamespace("future", quietly = TRUE) ||
    !requireNamespace("future.apply", quietly = TRUE) ||
    !requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Install ncdf4, future, future.apply, and Rcpp before running this script.")
}


require_climate_data(cfg)
data_dir <- cfg$data_root

# Validate before casting: as.integer() would silently turn 1.5 into 1 and
# malformed limits into a full, potentially expensive production run.
settings <- c(MVGC_SENSITIVITY_WORKERS = 4, MVGC_SENSITIVITY_LAT_BLOCK_SIZE = 2,
  MVGC_SENSITIVITY_LON_LIMIT = 0, MVGC_SENSITIVITY_LAT_LIMIT = 0, MVGC_SENSITIVITY_MODEL_LIMIT = 0)
for (name in names(settings)) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(settings[[name]]))))
  minimum <- if (grepl("_LIMIT$", name)) 0 else 1
  if (length(value) != 1L || !is.finite(value) || value < minimum || value != floor(value) ||
      value > .Machine$integer.max) stop(name, " must be an integer >= ", minimum, ".")
  settings[[name]] <- value
}
workers <- as.integer(settings[["MVGC_SENSITIVITY_WORKERS"]])
block_size <- as.integer(settings[["MVGC_SENSITIVITY_LAT_BLOCK_SIZE"]])

model_names <- scan(cfg$era5_models, what = "", quiet = TRUE)
reference_name <- "ERA5"
if (length(model_names) < 4L || anyDuplicated(model_names) || !reference_name %in% model_names)
  stop("Expected ERA5 and at least three distinct candidate models.")
variables <- cfg$variables
calibration_years <- cfg$calibration_years
k_values <- 4:12
lon <- cfg$lon
lat <- cfg$lat

range_file <- resolve_experiment_support("era5")
range_cube <- build_range_cube(range_file, variables)
lon_index <- match(lon, -180:179); lat_index <- match(lat, -90:90)
if (anyNA(lon_index) || anyNA(lat_index)) stop("Requested grid is outside the graph-cut support grid.")
range_var <- range_cube[lon_index, lat_index, , , drop = FALSE]
rm(range_cube)

# Optional limits support quick integration checks; unset values run the full experiment.
lon_limit <- as.integer(settings[["MVGC_SENSITIVITY_LON_LIMIT"]])
lat_limit <- as.integer(settings[["MVGC_SENSITIVITY_LAT_LIMIT"]])
model_limit <- as.integer(settings[["MVGC_SENSITIVITY_MODEL_LIMIT"]])
if (is.finite(lon_limit) && lon_limit > 0L) {
  keep <- seq_len(min(lon_limit, length(lon)))
  lon <- lon[keep]
  range_var <- range_var[keep, , , , drop = FALSE]
}
if (is.finite(lat_limit) && lat_limit > 0L) {
  keep <- seq_len(min(lat_limit, length(lat)))
  lat <- lat[keep]
  range_var <- range_var[, keep, , , drop = FALSE]
}
if (is.finite(model_limit) && model_limit > 0L) {
  if (model_limit < 3L) stop("Rank stability needs at least three candidate models.")
  candidates <- setdiff(model_names, reference_name)
  model_names <- c(reference_name, head(candidates, model_limit))
}

output_dir <- mvgc_output_dir("histogram_resolution")
if (any(file.exists(file.path(output_dir, c("k_sensitivity_global.csv", "k_sensitivity_cellwise.rds",
    "k_sensitivity_rank_stability.png", "k_sensitivity_rank_stability.pdf")))))
  stop("Histogram-resolution outputs already exist; choose a new MVGC_OUTPUT_ROOT.")
set.seed(cfg$seed)

cat("Trivariate PMF bin sensitivity\n")
cat("  calibration:", min(calibration_years), "-", max(calibration_years), "\n")
cat("  split A: odd years; split B: even years\n")
cat("  k:", paste(k_values, collapse = ", "), "\n")
cat("  grid:", length(lon), "x", length(lat), "\n")
cat("  models:", length(model_names) - 1L, "candidates +", reference_name, "reference\n")

started <- Sys.time()
result <- run_k_sensitivity(
  model_names = model_names,
  reference_name = reference_name,
  variables = variables,
  data_dir = data_dir,
  years = calibration_years,
  lon = lon,
  lat = lat,
  range_var = range_var,
  k_values = k_values,
  workers = workers,
  latitude_block_size = block_size,
  strict_time_checks = cfg$strict_time_checks,
  season_selection_mode = cfg$season_selection_mode,
  time_decode_mode = cfg$time_decode_mode,
  requested_year_policy = cfg$requested_year_policy,
  verbose = TRUE,
  checkpoint_dir = mvgc_output_dir("checkpoints/histogram_resolution")
)
result$metadata$started_at <- format(started, "%Y-%m-%d %H:%M:%S %Z")
result$metadata$finished_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
result$metadata$data_dir <- data_dir
result$metadata$range_file <- range_file
result$metadata$input_provenance <- experiment_input_provenance(range_file)
result$metadata$production_k_unchanged <- cfg$nbins3d
result$metadata$seed <- cfg$seed
result$metadata$reduced_run <- any(c(lon_limit, lat_limit, model_limit) > 0L)

paths <- write_k_sensitivity_outputs(result, output_dir)
print(result$summary)
cat("Outputs:\n")
print(paths)
cat("Elapsed minutes:", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), "\n")

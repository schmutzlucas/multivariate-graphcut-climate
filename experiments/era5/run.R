# Production ERA5 workflow: load exact external prepared fields and shared
# supports, construct calibration/target PMFs, then invoke the pinned cylindrical
# solver. Scientific periods, models, season, bins, lambda and seed come from
# config; changing them defines a different experiment. Worker/block settings
# bound resources only. Input/solver identity checks fail before expensive work.
# The saved workspace retains PMFs, reference diagnostics, selected labels and
# provenance for evaluation; do not point output at an archived result directory.
source("R/load.R")
load_mvgc()
source("R/experiments/perfect_model.R") # Shared range-cube validation/reordering only.
cfg <- article1_config()
require_climate_data(cfg)
require_packages(c("ncdf4", "future", "future.apply"))
require_production_solver(cfg)
future::plan(future::sequential)
workspace_dir <- mvgc_output_dir("workspaces")
support_file <- resolve_experiment_support("era5")
range_var_gc <- build_range_cube(support_file, cfg$variables)

# Setting variables
lon <- cfg$lon
lat <- cfg$lat
lon_size <- length(lon)
lat_size <- length(lat)
# Final paper temporal ranges. These match the 25-year PME calibration and
# the latest complete 25-year projection window available in ERA5.
year_present <- cfg$calibration_years
year_future <- cfg$era5_target_years

if (length(year_present) != 25L || length(year_future) != 25L) {
  stop("The final ERA5 calibration and projection periods must each contain 25 years.")
}
if (max(year_present) >= min(year_future)) {
  stop("The final ERA5 calibration and projection periods must not overlap.")
}
period_tag <- paste0(
  "cal", min(year_present), "-", max(year_present),
  "_proj", min(year_future), "-", max(year_future)
)

# Seasonal window (user-settable): month/day
# Example summer: April 15 to October 14
season_start_md <- cfg$season_start_md
season_end_md <- cfg$season_end_md

# data directory
data_dir <- cfg$data_root
if (!nzchar(data_dir)) stop("Set MVGC_DATA_ROOT to the ERA5/CMIP6 data directory.")


# List of the variable used
variables <- cfg$variables

## 1.  bin-resolution choices
nbins3d <- cfg$nbins3d
nbins2d <- cfg$nbins2d
nbins1d <- cfg$nbins1d

## 2.  model list
model_names <- scan(cfg$era5_models, what = "", quiet = TRUE)
if (length(model_names) < 3L || model_names[1] != "ERA5" || anyDuplicated(model_names))
  stop("ERA5 requires its reference first and at least two distinctly named candidate models.")

## 3.  number of parallel workers
resource_defaults <- c(MVGC_WORKERS = 4L, MVGC_LAT_BLOCK_SIZE = 8L)
resource_settings <- lapply(names(resource_defaults), function(setting_name) {
  value <- suppressWarnings(as.numeric(Sys.getenv(setting_name, unset = as.character(resource_defaults[[setting_name]]))))
  if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value) || value > .Machine$integer.max)
    stop(setting_name, " must be a positive integer.")
  as.integer(value)
})
workers <- resource_settings[[1]]
spatial_block_size <- resource_settings[[2]]
strict_time_checks <- cfg$strict_time_checks
season_selection_mode <- cfg$season_selection_mode
time_decode_mode <- cfg$time_decode_mode
requested_year_policy <- cfg$requested_year_policy
seasonal_diagnostics <- TRUE
seasonal_verbose <- FALSE

run_metadata <- list(
  script = "experiments/era5/run.R",
  run_state = "building_pmfs",
  started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  calibration_years = year_present,
  projection_years = year_future,
  period_tag = period_tag,
  period_length_years = length(year_present),
  season_start_md = season_start_md,
  season_end_md = season_end_md,
  spatial_block_size = spatial_block_size,
  workers = workers,
  seed = cfg$seed,
  nbins3d = nbins3d,
  nbins2d = nbins2d,
  nbins1d = nbins1d,
  requested_year_policy = requested_year_policy,
  model_list_file = cfg$era5_models,
  range_file = support_file,
  input_provenance = experiment_input_provenance(support_file),
  source_sha256 = scientific_code_fingerprint(c("experiments/era5/run.R", "config/article1.R",
    "R/distributions/pmf_pipeline.R", "R/distributions/histogram_fast.R",
    "R/distributions/histogram_reference.R", "R/preprocessing/calendar.R",
    "R/preprocessing/cf_calendar_validation.R", "R/preprocessing/netcdf_io.R",
    "R/metrics/partial_hellinger.R", "R/graph_cut/solver.R", "R/graph_cut/energy.R",
    "R/experiments/perfect_model.R", "R/utils/era5_workspace.R"))
)

cat("Final ERA5 paper experiment\n")
cat("  calibration:", min(year_present), "-", max(year_present), "\n")
cat("  projection: ", min(year_future), "-", max(year_future), "\n")
cat("  output tag: ", period_tag, "\n")

## 4.  call the multi-resolution histogram builder
cat("-> building PDFs and means ...\n")
t0 <- Sys.time()
tmp <- compute_nd_pdf_multi_seasonal(
  variables, model_names,
  data_dir,
  year_present, year_future,
  lon, lat,
  range_var_gc,
  season_start_md = season_start_md,
  season_end_md   = season_end_md,
  nbins3d = nbins3d,
  nbins2d = nbins2d,
  nbins1d = nbins1d,
  workers = workers,
  strict_time_checks = strict_time_checks,
  season_selection_mode = season_selection_mode,
  time_decode_mode = time_decode_mode,
  requested_year_policy = requested_year_policy,
  spatial_block_size = spatial_block_size,
  diagnostics = seasonal_diagnostics,
  verbose = seasonal_verbose
)
cat("   done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")

if (isTRUE(seasonal_diagnostics)) {
  cat("Detected calendars by model:\n")
  print(tmp$calendars)

  sample_model <- names(tmp$diagnostics)[1]
  sample_var <- variables[1]
  cat("\nSample seasonal selection counts for ", sample_model, " / ", sample_var, ":\n", sep = "")
  print(tmp$diagnostics[[sample_model]]$variables[[sample_var]]$present_selection)
  cat("\nSample decode backend for ", sample_model, " / ", sample_var, ":\n", sep = "")
  print(
    tmp$diagnostics[[sample_model]]$variables[[sample_var]][
      c("decode_backend_requested", "decode_backend", "fallback_used", "fallback_reason")
    ]
  )
}

gc()

# --------------------------------------------------------------------
#  The configured ERA5 reference is first; later indices identify candidates.
# --------------------------------------------------------------------
ref_index      <- 1L
reference_name <- model_names[ref_index]
model_names    <- model_names[-ref_index]

pdf3_present <- tmp$pdf3$present          # 3-D joint PDFs
pdf3_future  <- tmp$pdf3$future
expected_pdf_dimensions <- c(lon_size, lat_size, as.integer(nbins3d^3), length(model_names) + 1L)
if (!identical(dim(pdf3_present), expected_pdf_dimensions) || !identical(dim(pdf3_future), expected_pdf_dimensions))
  stop("ERA5 joint PMF dimensions do not match grid, bins and configured models.")
# Validate both periods, including the reference, before scores or native callbacks.
# A zero-observation PMF must not become a perfect calibration match via NaN -> 0.
validate_target_pmfs(pdf3_present)
validate_target_pmfs(pdf3_future)

pdf2_present <- tmp$pdf2$present          # list of 3 two-D vectors
pdf2_future  <- tmp$pdf2$future

pdf1_present <- tmp$pdf1$present          # list of 3 one-D vectors
pdf1_future  <- tmp$pdf1$future

mean_present <- tmp$mean$present          # cell-wise means
mean_future  <- tmp$mean$future

# drop the tmp list if memory is tight
rm(tmp)
gc()

# --------------------------------------------------------------------
#  Save the workspace
# --------------------------------------------------------------------
stamp    <- format(Sys.time(), "%Y%m%d%H%M")
filename <- file.path(
  workspace_dir,
  paste0(stamp, "_workspace_multiRes_3v_", period_tag, ".RData")
)
run_metadata$run_state <- "pmfs_complete"
if (file.exists(filename)) stop("PMF workspace exists; choose a new MVGC_OUTPUT_ROOT: ", filename)
save_era5_workspace(filename, stage = "pmfs")
cat("workspace saved to", filename, "\n")


# --------------------------------------------------------------------
#  1) Split out reference vs. models
# --------------------------------------------------------------------
# assume: pdf3_present/future and ref_index are in scope
pdf3_ref_pres    <- pdf3_present[,,, ref_index]
pdf3_models_pres <- pdf3_present[,,,-ref_index]

pdf3_ref_fut    <- pdf3_future[,,, ref_index]
pdf3_models_fut <- pdf3_future[,,,-ref_index]

rm(pdf3_future, pdf3_present)

# --------------------------------------------------------------------
#  2) build all-bins index list for full Hellinger
# --------------------------------------------------------------------
nbins_total3d <- nbins3d^3
selected_indices_all <- lapply(seq_len(lon_size), function(i) {
  lapply(seq_len(lat_size), function(j) {
    seq_len(nbins_total3d)
  })
})

# --------------------------------------------------------------------
#  3) compute full Hellinger distance maps
# --------------------------------------------------------------------
h_dist_pres <- compute_partial_hdist(
  pdf_ref      = pdf3_ref_pres,
  pdf_model    = pdf3_models_pres,
  selected_indices = selected_indices_all
)

h_dist_fut  <- compute_partial_hdist(
  pdf_ref      = pdf3_ref_fut,
  pdf_model    = pdf3_models_fut,
  selected_indices = selected_indices_all
)

# --------------------------------------------------------------------
#  4) invalid scores fail rather than being treated as perfect matches
# --------------------------------------------------------------------
if (any(!is.finite(h_dist_pres)) || any(h_dist_pres < 0) ||
    any(!is.finite(h_dist_fut)) || any(h_dist_fut < 0)) stop("ERA5 Hellinger score arrays are invalid.")

# quick check
hist(h_dist_pres, main="H-dist Present (all bins)")
hist(h_dist_fut,  main="H-dist Future  (all bins)")

# --------------------------------------------------------------------
#  5) run GraphCut on the 3-D cost map
# --------------------------------------------------------------------
smooth_cost <- cfg$lambda
run_metadata$smooth_cost <- smooth_cost

GC_result <- tryCatch({
  fit <- GraphCutHellinger_nD_lat(
    pdf_models_future = pdf3_models_fut,   # target PMFs define pairwise seam costs
    h_dist            = h_dist_pres,        # datacost = Hellinger(pres)
    weight_data       = 1,
    weight_smooth     = smooth_cost,
    nBins             = nbins_total3d,
    lat               = lat,
    seed              = cfg$seed,
    verbose           = TRUE,
    rebuild           = TRUE
  )
  labels <- fit$label_attribution
  if (!is.numeric(labels) || !identical(dim(labels), c(lon_size, lat_size)) ||
      any(!is.finite(labels)) || any(labels != floor(labels)) || any(labels < 1 | labels > length(model_names)))
    stop("GraphCut returned invalid ERA5 candidate labels.")
  fit
}, error = function(e) {
  run_metadata$run_state <<- "failed_solver"
  run_metadata$solver_error <<- conditionMessage(e)
  failure_file <- file.path(workspace_dir, paste0(stamp, "_ERA5_solver_failure.rds"))
  saved <- try(saveRDS(list(metadata = run_metadata, pmf_workspace = filename), failure_file), silent = TRUE)
  if (inherits(saved, "try-error")) warning("Could not save solver failure metadata: ", failure_file, call. = FALSE)
  stop("ERA5 GraphCut failed: ", conditionMessage(e), ". PMFs retained in ", filename, call. = FALSE)
})
gc()



# Compute the pdf3 for MMM and GC | Compute the hdist_3
{
  # Initialize arrays for GraphCut PDFs
  pdf3_GC_pres <- array(NA, dim = c(length(lon), length(lat), dim(pdf3_models_pres)[3]))
  pdf3_GC_fut <- array(NA, dim = c(length(lon), length(lat), dim(pdf3_models_fut)[3]))

  # Assign PDFs based on GraphCut label attribution
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Get the GraphCut-selected label (indexing in R is 1-based)
      label <- GC_result$label_attribution[i, j]  # [longitude, latitude], one-based candidate index

      # Assign the PDF of the selected model
      pdf3_GC_pres[i, j, ] <- pdf3_models_pres[i, j, , label]
      pdf3_GC_fut[i, j, ] <- pdf3_models_fut[i, j, , label]
    }
  }


  # Initialize a lon x lat matrix for each smooth cost
  GC_hdist_pres <- matrix(NA, nrow = length(lon), ncol = length(lat))
  GC_hdist_fut <- matrix(NA, nrow = length(lon), ncol = length(lat))

  for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
    islabel <- which(GC_result$label_attribution == l)
    GC_hdist_pres[islabel] <- h_dist_pres[,,l][islabel]
    GC_hdist_fut[islabel] <- h_dist_fut[,, l][islabel]
  }


  hist(GC_hdist_fut)


  mean(GC_hdist_fut)

  # Compute Multi-Model Mean for Present
  pdf3_MMM_pres <- apply(pdf3_models_pres, c(1, 2, 3), mean)

  # Compute Multi-Model Mean for Future
  pdf3_MMM_fut <- apply(pdf3_models_fut, c(1, 2, 3), mean)


  # Initialize arrays to store the Hellinger distance for MMM
  MMM_hdist_pres <- array(NA, dim = c(length(lon), length(lat)))
  MMM_hdist_fut <- array(NA, dim = c(length(lon), length(lat)))

  # Compute Hellinger distance for the Multi-Model Mean
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance for present
      MMM_hdist_pres[i, j] <- sqrt(sum((sqrt(pdf3_MMM_pres[i, j, ]) - sqrt(pdf3_ref_pres[i, j, ]))^2)) / sqrt(2)

      # Compute Hellinger distance for future
      MMM_hdist_fut[i, j] <- sqrt(sum((sqrt(pdf3_MMM_fut[i, j, ]) - sqrt(pdf3_ref_fut[i, j, ]))^2)) / sqrt(2)
    }
  }

  if (any(!is.finite(MMM_hdist_pres)) || any(!is.finite(MMM_hdist_fut)) ||
      any(!is.finite(GC_hdist_pres)) || any(!is.finite(GC_hdist_fut)))
    stop("ERA5 final Hellinger maps contain nonfinite values.")

  # Visualize or analyze the results
  hist(MMM_hdist_pres)
  hist(GC_hdist_pres)
  hist(MMM_hdist_fut)
  hist(GC_hdist_fut)
}

# --------------------------------------------------------------------
#  Save the workspace
# --------------------------------------------------------------------
stamp    <- format(Sys.time(), "%Y%m%d%H%M")
run_metadata$finished_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
run_metadata$run_state <- "complete"
filename <- file.path(
  workspace_dir,
  paste0(
    stamp,
    "_workspace_multiRes_3v_gc_results_ERA5_fullyear_",
    period_tag,
    ".RData"
  )
)
if (file.exists(filename)) stop("Final workspace exists; choose a new MVGC_OUTPUT_ROOT: ", filename)
save_era5_workspace(filename)
cat("workspace saved to", filename, "\n")

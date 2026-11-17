# Production PME workflow: external prepared climate files/supports and pinned
# cylindrical solver are required; Level 1 reproduction never enters this script.
# Scientific controls: model/reference membership, periods, season, PMF bins,
# supports, lambdas and seed define the experiment. PME_MAX_MODELS truncates the
# scientific population and must not be confused with a worker/resource limit.
# Worker, block-size and memory settings only bound execution resources.
# Each reference is excluded from its candidate set; shared PMFs are prepared per
# target period before calibration unary and target-period seam costs reach GCO.
# Resume accepts only complete, identity-compatible period/reference checkpoints.
# PME_FORCE_RECOMPUTE=1 intentionally bypasses reuse and overwrites checkpoint
# paths. Use a distinct PME_RUN_STAMP/output directory to retain an earlier run.
# Saved checkpoints contain metadata, input/source identities and lambda results;
# the final compact bundle combines periods, references, diagnostics and failures.
# Solver failures retain successful checkpoints but make the final run fail.
source("R/load.R")
load_mvgc()
cfg <- article1_config()
require_climate_data(cfg)
require_packages(c("ncdf4", "future", "future.apply"))
require_production_solver(cfg)
future::plan(future::sequential)
workspace_dir <- mvgc_output_dir("workspaces")
source("R/experiments/perfect_model.R")
# Global configuration
variables <- cfg$variables
lon <- cfg$lon
lat <- cfg$lat
lon_size <- length(lon)
lat_size <- length(lat)

calibration_years <- cfg$calibration_years
period_role <- "fixed_calibration_hybrid_projection_smoothness"
year_present <- calibration_years
periods <- filter_periods_from_env(cfg$periods)
if (any((periods$end_year - periods$start_year + 1L) != 25L)) {
  stop("All projection periods must have exactly 25 years.")
}

season_start_md <- cfg$season_start_md
season_end_md <- cfg$season_end_md

data_dir <- cfg$data_root
if (!nzchar(data_dir)) stop("Set MVGC_DATA_ROOT to the ERA5/CMIP6 data directory.")

nbins3d <- cfg$nbins3d
nbins2d <- cfg$nbins2d
nbins1d <- cfg$nbins1d
nbins_total3d <- nbins3d^3

smooth_costs <- cfg$pme_lambdas
custom_lambdas <- Sys.getenv("PME_LAMBDAS", "")
if (nzchar(custom_lambdas)) {
  smooth_costs <- as.numeric(strsplit(custom_lambdas, ",", fixed = TRUE)[[1]])
  stopifnot(length(smooth_costs) > 0L, all(is.finite(smooth_costs)), all(smooth_costs >= 0), !anyDuplicated(smooth_costs))
}
smooth_keys <- vapply(smooth_costs, lambda_to_key, character(1))
if(anyDuplicated(smooth_keys))stop("Distinct lambda values collide after checkpoint-key formatting.")
partial_tau <- cfg$pme_ldr_thresholds
partial_tau_keys <- vapply(partial_tau, tau_to_key, character(1))
partial_h_definition <- "Corrected 2D MMM-vs-reference partial H over reference projection-period LDR bins for multiple tau values."

# Execution controls: bound workers without changing scientific membership.
workers <- read_positive_int_env("PME_WORKERS", default_value = 3L)
strict_time_checks <- cfg$strict_time_checks
season_selection_mode <- cfg$season_selection_mode
time_decode_mode <- cfg$time_decode_mode
requested_year_policy <- cfg$requested_year_policy
seasonal_diagnostics <- FALSE
seasonal_verbose <- FALSE
force_recompute <- Sys.getenv("PME_FORCE_RECOMPUTE", unset = "0") == "1"

model_list_file <- "config/pme_models.txt"
model_names <- scan(model_list_file, what = "", quiet = TRUE)

max_models <- suppressWarnings(as.integer(Sys.getenv("PME_MAX_MODELS", unset = "0")))
if (!is.na(max_models) && max_models > 0 && max_models < length(model_names)) {
  model_names <- model_names[seq_len(max_models)]
  cat("Limiting run to", max_models, "models due to PME_MAX_MODELS.\n")
}

if (length(model_names) < 3L) {
  stop("Native PME graph cut requires at least three models: one held out and at least two candidate labels.")
}

run_stamp <- Sys.getenv("PME_RUN_STAMP", unset = "")
if (run_stamp == "") {
  run_stamp <- format(Sys.time(), "%Y%m%d%H%M")
}
output_dir <- Sys.getenv("PME_OUTPUT_DIR", unset = workspace_dir)
checkpoint_dir <- Sys.getenv(
  "PME_CHECKPOINT_DIR",
  unset = file.path(output_dir, "pme_final_period_checkpoints", run_stamp)
)
period_output_dir <- file.path(output_dir, "pme_final_period_results", run_stamp)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(period_output_dir, recursive = TRUE, showWarnings = FALSE)

# External supports are shared across references and periods, not refitted here.
range_file <- resolve_experiment_support("pme")
if (!file.exists(range_file)) {
  stop(
    "PME full-year range file is missing: ",
    range_file,
    ". Run experiments/perfect_model/build_supports.R before this experiment."
  )
}

cat("Loading PME full-year range cube...\n")
range_var_gc <- build_range_cube(range_file, variables = variables)
gc()

all_bins <- seq_len(nbins_total3d)
shared_lat_indices <- rep(list(all_bins), lat_size)
selected_indices_all <- rep(list(shared_lat_indices), lon_size)

scientific_settings<-pme_scientific_settings(cfg)
scientific_settings$parameters$pme_lambdas<-smooth_costs
scientific_settings$parameters$pme_ldr_thresholds<-partial_tau
meta_base <- list(
  scientific_settings = scientific_settings,
  script = "experiments/perfect_model/run.R",
  run_stamp = run_stamp,
  created_at = as.character(Sys.time()),
  model_list_file = model_list_file,
  model_names = model_names,
  variables = variables,
  lon = lon,
  lat = lat,
  year_present = year_present,
  calibration_years = calibration_years,
  periods = periods,
  projection_periods = periods,
  period_role = period_role,
  period_length_years = 25L,
  graphcut_data_term = "calibration-period model-vs-reference H",
  graphcut_smoothness_term = "projection-period model PDF similarity",
  evaluation_period = "projection period shown by period_id",
  MMM_definition = "local ensemble mean PDF across non-reference models",
  MMM_full_hdist = "H(MMM_projection_pdf, ref_projection_pdf)",
  MMM_partial_hdist = "H(MMM_projection_pdf, ref_projection_pdf) over reference projection LDR bins",
  season_start_md = season_start_md,
  season_end_md = season_end_md,
  nbins3d = nbins3d,
  nbins2d = nbins2d,
  nbins1d = nbins1d,
  lambda = smooth_costs,
  lambda_key = smooth_keys,
  lambda_keys = smooth_keys,
  partial_tau = partial_tau,
  partial_tau_keys = partial_tau_keys,
  partial_h_region = "LDR",
  workers = workers,
  strict_time_checks = strict_time_checks,
  season_selection_mode = season_selection_mode,
  time_decode_mode = time_decode_mode,
  requested_year_policy = requested_year_policy,
  seasonal_diagnostics = seasonal_diagnostics,
  data_dir = data_dir,
  range_sources = range_file,
  range_combination = "single common support pooled across exactly the configured PME models and required years",
  input_provenance = experiment_input_provenance(range_file),
  checkpoint_dir = checkpoint_dir,
  period_output_dir = period_output_dir,
  force_recompute = force_recompute,
  partial_h_definition = partial_h_definition
)

GC_result_list <- list()
GC_hdist_future_list <- list()
GC_partial_hdist_future_list <- list()
MMM_hdist_future_list <- list()
MMM_partial_hdist_future_list <- list()
calendar_by_period <- list()
period_files <- character(0)
checkpoint_files <- character(0)
completed_period_references <- data.frame(
  period_id = character(0),
  ref_name = character(0),
  stringsAsFactors = FALSE
)

for (p in seq_len(nrow(periods))) {
  period_row <- periods[p, , drop = FALSE]
  period_id <- period_row$period_id
  year_future <- period_row$start_year:period_row$end_year

  cat("\n=== Projection period ", period_row$label, " (", period_id, ") ===\n", sep = "")

  checkpoint_paths <- vapply(
    seq_along(model_names),
    function(m) checkpoint_path_for(checkpoint_dir, period_id, m, model_names[m]),
    character(1)
  )
  compatible_checkpoints <- vapply(
    checkpoint_paths,
    checkpoint_is_compatible,
    logical(1),
    period_row = period_row,
    calibration_years = calibration_years,
    period_role = period_role,
    lambda_keys = smooth_keys,
    tau_keys = partial_tau_keys,
    expected_model_names = model_names,
    range_sources = range_file,
    partial_h_definition = partial_h_definition,
    expected_settings = scientific_settings
  )
  period_already_complete <- all(compatible_checkpoints) && !force_recompute

  if (period_already_complete) {
    cat("All checkpoints already exist for ", period_id, "; loading compact period result.\n", sep = "")
    period_result <- read_period_from_checkpoints(
      period_row = period_row,
      model_names = model_names,
      checkpoint_dir = checkpoint_dir,
      lambda_keys = smooth_keys,
      tau_keys = partial_tau_keys,
      expected_dim = c(lon_size, lat_size),
      calibration_years = calibration_years,
      period_role = period_role,
      expected_model_names = model_names,
      range_sources = range_file,
      partial_h_definition = partial_h_definition,
      expected_settings = scientific_settings
    )
  } else {
    cat("-> building PDFs for PME model set and future period ", period_row$label, " ...\n", sep = "")
    t0 <- Sys.time()
    tmp <- compute_nd_pdf_multi_seasonal(
      variables = variables,
      model_names = model_names,
      data_dir = data_dir,
      year_present = year_present,
      year_future = year_future,
      lon = lon,
      lat = lat,
      range_var = range_var_gc,
      season_start_md = season_start_md,
      season_end_md = season_end_md,
      nbins3d = nbins3d,
      nbins2d = nbins2d,
      nbins1d = nbins1d,
      workers = workers,
      strict_time_checks = strict_time_checks,
      season_selection_mode = season_selection_mode,
      time_decode_mode = time_decode_mode,
      requested_year_policy = requested_year_policy,
      diagnostics = seasonal_diagnostics,
      verbose = seasonal_verbose
    )
    cat("PDF build done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes.\n")

    pdf3_present <- tmp$pdf3$present
    pdf3_future <- tmp$pdf3$future
    calendar_by_period[[period_id]] <- tmp$calendars
    rm(tmp)
    gc()
    # Include the held-out reference in the validity gate, not only the 14
    # candidates passed to the solver. Empty reference histograms are not PMFs.
    validate_target_pmfs(pdf3_present)
    validate_target_pmfs(pdf3_future)

    period_result <- read_period_from_checkpoints(
      period_row = period_row,
      model_names = model_names,
      checkpoint_dir = checkpoint_dir,
      lambda_keys = smooth_keys,
      tau_keys = partial_tau_keys,
      expected_dim = c(lon_size, lat_size),
      calibration_years = calibration_years,
      period_role = period_role,
      expected_model_names = model_names,
      range_sources = range_file,
      partial_h_definition = partial_h_definition,
      expected_settings = scientific_settings
    )

    for (m in seq_along(model_names)) {
      ref_name <- model_names[m]
      checkpoint_file <- checkpoint_path_for(checkpoint_dir, period_id, m, ref_name)
      dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)

      existing_checkpoint_ok <- checkpoint_is_compatible(
        checkpoint_file = checkpoint_file,
        period_row = period_row,
        calibration_years = calibration_years,
        period_role = period_role,
        lambda_keys = smooth_keys,
        tau_keys = partial_tau_keys,
        expected_model_names = model_names,
        range_sources = range_file,
        partial_h_definition = partial_h_definition,
        expected_settings = scientific_settings
      )
      if (existing_checkpoint_ok && !force_recompute) {
        cat("[", period_id, " ", m, "/", length(model_names), "] reference = ", ref_name, " already checkpointed; skipping.\n", sep = "")
        next
      }

      cat("[", period_id, " ", m, "/", length(model_names), "] reference = ", ref_name, "\n", sep = "")
      ref_start <- Sys.time()

      pdf3_ref_pres <- pdf3_present[, , , m]
      pdf3_models_pres <- pdf3_present[, , , -m, drop = FALSE]
      pdf3_ref_fut <- pdf3_future[, , , m]
      pdf3_models_fut <- pdf3_future[, , , -m, drop = FALSE]

      h_dist_pres <- compute_partial_hdist(
        pdf_ref = pdf3_ref_pres,
        pdf_model = pdf3_models_pres,
        selected_indices = selected_indices_all
      )
      h_dist_fut <- compute_partial_hdist(
        pdf_ref = pdf3_ref_fut,
        pdf_model = pdf3_models_fut,
        selected_indices = selected_indices_all
      )

      if (any(!is.finite(h_dist_pres)) || any(!is.finite(h_dist_fut)))
        stop("Nonfinite reference/candidate Hellinger distances; inspect PMFs rather than replacing scores with zero.")

      # One identical equal-model mixture serves ordinary Hellinger and every
      # LDR threshold. Averaging 14 global PMFs again for each threshold is wasteful.
      mmm_pdf_fut <- compute_mmm_pdf(pdf3_models_fut)
      mmm_hdist_fut <- compute_mmm_hdist(pdf_ref = pdf3_ref_fut,
        pdf_models = pdf3_models_fut, mmm_pdf = mmm_pdf_fut)

      validate_map_dim(mmm_hdist_fut, c(lon_size, lat_size), paste0("MMM future map for ", period_id, "/", ref_name))

      ldr_indices_by_tau <- setNames(vector("list", length(partial_tau_keys)), partial_tau_keys)
      mmm_partial_by_tau <- setNames(vector("list", length(partial_tau_keys)), partial_tau_keys)
      partial_fut_by_tau <- setNames(vector("list", length(partial_tau_keys)), partial_tau_keys)

      for (tau_idx in seq_along(partial_tau)) {
        tau_value <- partial_tau[tau_idx]
        tau_key <- partial_tau_keys[tau_idx]

        ldr_indices <- build_ldr_indices(
          pdf_ref = pdf3_ref_fut,
          tau = tau_value,
          lon_size = lon_size,
          lat_size = lat_size
        )

        mmm_partial_hdist_fut <- compute_mmm_partial_hdist(
          pdf_ref = pdf3_ref_fut,
          pdf_models = pdf3_models_fut,
          selected_indices = ldr_indices,
          mmm_pdf = mmm_pdf_fut
        )

        partial_fut <- compute_partial_hdist(
          pdf_ref = pdf3_ref_fut,
          pdf_model = pdf3_models_fut,
          selected_indices = ldr_indices
        )
        if (any(!is.finite(partial_fut)))
          stop("Nonfinite reference/candidate LDR scores; inspect selected bins and PMFs.")

        validate_map_dim(
          mmm_partial_hdist_fut,
          c(lon_size, lat_size),
          paste0("MMM partial map for ", period_id, "/", ref_name, " @ tau ", tau_key)
        )

        ldr_indices_by_tau[[tau_key]] <- ldr_indices
        mmm_partial_by_tau[[tau_key]] <- mmm_partial_hdist_fut
        partial_fut_by_tau[[tau_key]] <- partial_fut
      }

      gc_result_by_lambda <- setNames(vector("list", length(smooth_keys)), smooth_keys)
      gc_future_by_lambda <- setNames(vector("list", length(smooth_keys)), smooth_keys)
      gc_partial_by_tau_lambda <- setNames(
        lapply(partial_tau_keys, function(tau_key) setNames(vector("list", length(smooth_keys)), smooth_keys)),
        partial_tau_keys
      )

      for (lambda_idx in seq_along(smooth_costs)) {
        smooth_cost <- smooth_costs[lambda_idx]
        smooth_key <- smooth_keys[lambda_idx]

        # Calibration unary fit and target-period candidate compatibility enter
        # the pinned float callback solver; target-reference scores evaluate only.
        gc_result <- tryCatch({
          GraphCutHellinger_nD_lat(
            pdf_models_future = pdf3_models_fut,
            h_dist = h_dist_pres,
            weight_data = 1,
            weight_smooth = smooth_cost,
            nBins = nbins_total3d,
            lat = lat,
            seed = cfg$seed,
            verbose = FALSE,
            rebuild = FALSE
          )
        }, error = function(e) {
          warning(
            "GraphCut failed for period '",
            period_id,
            "', reference '",
            ref_name,
            "' at lambda=",
            smooth_cost,
            ": ",
            e$message,
            call. = FALSE
          )
          NULL
        })

        if (is.null(gc_result)) {
          gc_future_map <- matrix(NA_real_, nrow = lon_size, ncol = lat_size)
        } else {
          gc_maps <- compute_gc_selected_maps(
            gc_result = gc_result,
            h_dist_pres = h_dist_pres,
            h_dist_fut = h_dist_fut,
            partial_fut = partial_fut_by_tau[[partial_tau_keys[1]]]
          )
          gc_future_map <- gc_maps$future
        }

        validate_map_dim(gc_future_map, c(lon_size, lat_size), paste0("GC future map for ", period_id, "/", ref_name, " @ ", smooth_key))

        gc_result_by_lambda[smooth_key] <- list(gc_result)
        gc_future_by_lambda[[smooth_key]] <- gc_future_map

        for (tau_key in partial_tau_keys) {
          if (is.null(gc_result)) {
            gc_partial_map <- matrix(NA_real_, nrow = lon_size, ncol = lat_size)
          } else {
            gc_maps_partial <- compute_gc_selected_maps(
              gc_result = gc_result,
              h_dist_pres = h_dist_pres,
              h_dist_fut = h_dist_fut,
              partial_fut = partial_fut_by_tau[[tau_key]]
            )
            gc_partial_map <- gc_maps_partial$partial_future
          }

          validate_map_dim(
            gc_partial_map,
            c(lon_size, lat_size),
            paste0("GC partial map for ", period_id, "/", ref_name, " @ ", smooth_key, ", tau ", tau_key)
          )
          gc_partial_by_tau_lambda[[tau_key]][[smooth_key]] <- gc_partial_map
        }
      }

      reference_successful<-all(vapply(gc_result_by_lambda,function(x)!is.null(x),logical(1)))&&
        all(vapply(gc_future_by_lambda,function(x)all(is.finite(x)),logical(1)))
      checkpoint_obj <- list(
        meta = utils::modifyList(
          meta_base,
          list(
            checkpoint_state = if(reference_successful)"complete"else"incomplete_solver_failure",
            checkpoint_created_at = as.character(Sys.time()),
            period = period_row,
            period_role = period_role,
            calibration_years = calibration_years,
            year_present = calibration_years,
            year_future = year_future,
            lambda = smooth_costs,
            lambda_key = smooth_keys,
            lambda_keys = smooth_keys,
            partial_tau = partial_tau,
            partial_tau_keys = partial_tau_keys,
            ref_index = m,
            ref_name = ref_name
          )
        ),
        period_id = period_id,
        ref_name = ref_name,
        GC_result_list = setNames(
          list(lapply(gc_result_by_lambda, function(x) setNames(list(x), ref_name))),
          period_id
        ),
        GC_hdist_future_list = setNames(
          list(lapply(gc_future_by_lambda, function(x) setNames(list(x), ref_name))),
          period_id
        ),
        GC_partial_hdist_future_list = setNames(
          list(lapply(
            gc_partial_by_tau_lambda,
            function(by_lambda) lapply(by_lambda, function(x) setNames(list(x), ref_name))
          )),
          period_id
        ),
        MMM_hdist_future_list = setNames(
          list(setNames(list(mmm_hdist_fut), ref_name)),
          period_id
        ),
        MMM_partial_hdist_future_list = setNames(
          list(lapply(mmm_partial_by_tau, function(x) setNames(list(x), ref_name))),
          period_id
        )
      )

      # Save each reference promptly, including an explicit failure state. A file
      # existing on disk alone is insufficient evidence for safe resume.
      saveRDS(checkpoint_obj, file = checkpoint_file, compress = "xz")

      for (smooth_key in smooth_keys) {
        period_result$GC_result_list[[smooth_key]][[ref_name]] <- gc_result_by_lambda[[smooth_key]]
        period_result$GC_hdist_future_list[[smooth_key]][[ref_name]] <- gc_future_by_lambda[[smooth_key]]
        for (tau_key in partial_tau_keys) {
          period_result$GC_partial_hdist_future_list[[tau_key]][[smooth_key]][[ref_name]] <- gc_partial_by_tau_lambda[[tau_key]][[smooth_key]]
        }
      }
      period_result$MMM_hdist_future_list[[ref_name]] <- mmm_hdist_fut
      for (tau_key in partial_tau_keys) {
        period_result$MMM_partial_hdist_future_list[[tau_key]][[ref_name]] <- mmm_partial_by_tau[[tau_key]]
      }
      period_result$checkpoint_files <- unique(c(period_result$checkpoint_files, checkpoint_file))
      if(reference_successful)
        period_result$completed_references <- unique(c(period_result$completed_references, ref_name))

      rm(
        pdf3_ref_pres, pdf3_models_pres, pdf3_ref_fut, pdf3_models_fut,
        h_dist_pres, h_dist_fut, mmm_pdf_fut, mmm_hdist_fut, ldr_indices_by_tau,
        mmm_partial_by_tau, partial_fut_by_tau, gc_result_by_lambda,
        gc_future_by_lambda, gc_partial_by_tau_lambda, checkpoint_obj
      )
      gc()

      cat(
        "Reference ",
        ref_name,
        " for ",
        period_id,
        " done in ",
        round(difftime(Sys.time(), ref_start, units = "mins"), 2),
        " minutes.\n",
        sep = ""
      )
    }

    rm(pdf3_present, pdf3_future)
    gc()
  }

  period_file <- write_period_result(
    period_result = period_result,
    period_row = period_row,
    period_output_dir = period_output_dir,
    run_stamp = run_stamp,
    meta_base = meta_base
  )

  period_files <- c(period_files, period_file)
  checkpoint_files <- c(checkpoint_files, period_result$checkpoint_files)
  completed_period_references <- rbind(
    completed_period_references,
    data.frame(
      period_id = period_id,
      ref_name = period_result$completed_references,
      stringsAsFactors = FALSE
    )
  )

  GC_result_list[[period_id]] <- period_result$GC_result_list
  GC_hdist_future_list[[period_id]] <- period_result$GC_hdist_future_list
  GC_partial_hdist_future_list[[period_id]] <- period_result$GC_partial_hdist_future_list
  MMM_hdist_future_list[[period_id]] <- period_result$MMM_hdist_future_list
  MMM_partial_hdist_future_list[[period_id]] <- period_result$MMM_partial_hdist_future_list

  cat("Period compact result saved to: ", period_file, "\n", sep = "")
  rm(period_result)
  gc()
}

for (p in seq_len(nrow(periods))) {
  period_id <- periods$period_id[p]
  for (ref_name in model_names) {
    validate_map_dim(
      MMM_hdist_future_list[[period_id]][[ref_name]],
      c(lon_size, lat_size),
      paste0("Final MMM future map for ", period_id, "/", ref_name)
    )
    validate_map_dim(
      MMM_partial_hdist_future_list[[period_id]][[partial_tau_keys[1]]][[ref_name]],
      c(lon_size, lat_size),
      paste0("Final MMM partial map for ", period_id, "/", ref_name, " @ tau ", partial_tau_keys[1])
    )
    validate_map_dim(
      GC_hdist_future_list[[period_id]][[smooth_keys[1]]][[ref_name]],
      c(lon_size, lat_size),
      paste0("Final GC future map for ", period_id, "/", ref_name, " @ ", smooth_keys[1])
    )
    for (tau_key in partial_tau_keys) {
      validate_map_dim(
        MMM_partial_hdist_future_list[[period_id]][[tau_key]][[ref_name]],
        c(lon_size, lat_size),
        paste0("Final MMM partial map for ", period_id, "/", ref_name, " @ tau ", tau_key)
      )
    }
    for (smooth_key in smooth_keys) {
      validate_map_dim(
        GC_hdist_future_list[[period_id]][[smooth_key]][[ref_name]],
        c(lon_size, lat_size),
        paste0("Final GC future map for ", period_id, "/", ref_name, " @ ", smooth_key)
      )
      for (tau_key in partial_tau_keys) {
        validate_map_dim(
          GC_partial_hdist_future_list[[period_id]][[tau_key]][[smooth_key]][[ref_name]],
          c(lon_size, lat_size),
          paste0("Final GC partial map for ", period_id, "/", ref_name, " @ ", smooth_key, ", tau ", tau_key)
        )
      }
    }
  }
}

final_meta <- c(
  meta_base,
  list(
    finished_at = as.character(Sys.time()),
    calendar_by_period = calendar_by_period,
    completed_period_references = completed_period_references,
    n_completed_period_references = nrow(completed_period_references),
    run_state = if(nrow(completed_period_references)==nrow(periods)*length(model_names))"complete"else"incomplete_solver_failure",
    checkpoint_files = unique(checkpoint_files),
    period_files = period_files
  )
)

final_bundle <- list(
  meta = final_meta,
  periods = periods,
  GC_result_list = GC_result_list,
  GC_hdist_future_list = GC_hdist_future_list,
  GC_partial_hdist_future_list = GC_partial_hdist_future_list,
  MMM_hdist_future_list = MMM_hdist_future_list,
  MMM_partial_hdist_future_list = MMM_partial_hdist_future_list
)

complete<-identical(final_meta$run_state,"complete")
final_bundle_file <- file.path(output_dir,paste0(run_stamp,
  if(complete)"_pme_final_periods_compact.rds"else"_pme_incomplete_periods_compact.rds"))
# Aggregate compact period/reference diagnostics, lambda results and identities.
# An incomplete bundle remains diagnostic evidence and triggers failure below.
saveRDS(final_bundle, file = final_bundle_file, compress = "xz")
if(!complete)stop("PME has failed solver cases. Successful checkpoints were retained; incomplete diagnostic bundle: ",final_bundle_file)

cat("\nFinal PME multi-period compact bundle saved to: ", final_bundle_file, "\n", sep = "")
cat("Per-period compact files saved to: ", period_output_dir, "\n", sep = "")
cat("Period/reference checkpoints saved to: ", checkpoint_dir, "\n", sep = "")

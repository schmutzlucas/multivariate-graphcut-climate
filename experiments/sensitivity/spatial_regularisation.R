# Sensitivity entry point: changing the tested resolution or regularisation
# changes the scientific configuration. Keep external inputs/supports fixed when
# comparing settings. Resource controls affect execution only; retain distinct
# output/checkpoint namespaces to avoid mixing incompatible experiments.
source("R/load.R")
load_mvgc()
cfg <- article1_config()
require_production_solver(cfg)
# Compute the Fig. 6 lambda candidates from the final matched-period ERA5 run.
# Run from the repository root. No climate inputs are rebuilt or modified.
# ERA5_RESULTS_FILE selects the saved experiment; MVGC_OUTPUT_ROOT selects output.
options(warn = 1)
workspace <- Sys.getenv("ERA5_RESULTS_FILE", unset = "")
if (!nzchar(workspace) || !file.exists(workspace)) {
  stop("Set ERA5_RESULTS_FILE to the final ERA5 result workspace.")
}
out <- mvgc_output_dir("spatial_regularisation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
all_lambdas <- cfg$era5_lambda_sensitivity
# Optional disjoint workers can precompute slow independent solutions. Each
# solution has a dedicated checkpoint; a full bundle is written only when ready.
requested <- Sys.getenv("MVGC_LAMBDA_VALUES", unset = "")
lambdas <- if (nzchar(requested)) suppressWarnings(as.numeric(strsplit(requested, ",", fixed = TRUE)[[1]])) else all_lambdas
stopifnot(length(lambdas) > 0, !anyNA(lambdas), !anyDuplicated(lambdas), all(lambdas %in% all_lambdas))
#' Format a spatial-sensitivity lambda as a decimal list key
#' @param x One finite nonnegative numeric lambda.
#' @return Character scalar with scientific notation disabled and whitespace trimmed.
key <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) stop("Expected one finite nonnegative lambda.")
  format(x, scientific = FALSE, trim = TRUE)
}
#' Resolve a lambda checkpoint in the current scientific-identity namespace
#' @param x One finite nonnegative numeric lambda.
#' @return Path below caller-environment checkpoint_dir. Does not create a file;
#'   each full scientific input/settings signature has its own parent directory.
checkpoint_file <- function(x) file.path(checkpoint_dir, paste0("lambda_", key(x), ".rds"))

#' Check a spatial-sensitivity checkpoint for scientific identity and success
#' @param record Deserialised checkpoint list, or malformed value to reject.
#' @param lambda Expected scalar nonnegative numeric lambda.
#' @param signature Exact expected scientific/input signature list.
#' @param grid Expected longitude/latitude dimensions, a two-element integer vector.
#' @param n_models Positive integer candidate-model count.
#' @return TRUE only for a complete matching checkpoint with finite integer labels
#'   in range and finite nonnegative native data/smooth energies. Malformed or
#'   unsuccessful payloads return FALSE, never partial numerical approval.
spatial_checkpoint_is_complete <- function(record, lambda, signature, grid, n_models) {
  tryCatch({
    if (!is.list(record) || !identical(record$state, "complete") ||
        !is.numeric(record$lambda) || length(record$lambda) != 1L ||
        !is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0 ||
        !identical(as.numeric(record$lambda), as.numeric(lambda)) ||
        !is.list(signature) || !identical(record$input_signature, signature) ||
        !is.numeric(grid) || length(grid) != 2L || any(!is.finite(grid)) || any(grid < 1 | grid != floor(grid)) ||
        !is.numeric(n_models) || length(n_models) != 1L || !is.finite(n_models) ||
        n_models < 1 || n_models != floor(n_models) || !is.list(record$result)) return(FALSE)
    labels <- record$result$label_attribution
    energies <- record$result[["Data and smooth cost"]]
    if (!is.numeric(labels) || !is.matrix(labels) || !identical(dim(labels), as.integer(grid)) ||
        any(!is.finite(labels)) || any(labels != floor(labels)) || any(labels < 1 | labels > n_models) ||
        !is.list(energies) || !identical(names(energies), c("Data cost", "Smooth cost")) ||
        !all(vapply(energies, function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0, logical(1))))
      return(FALSE)
    TRUE
  }, error = function(e) FALSE)
}

cat("Loading final ERA5 inputs: ", workspace, "\n", sep = "")
e <- new.env(parent = emptyenv())
load(workspace, envir = e)
if(nzchar(cfg$input_registry)) {
  current_support<-resolve_experiment_support("era5")
  if(!identical(e$run_metadata$input_provenance,experiment_input_provenance(current_support)))
    stop("Saved ERA5 workspace does not belong to the current corrected input/support preparation")
}
models <- scan(cfg$era5_models, what = "", quiet = TRUE)
expected_models <- models[-1L]
grid <- c(length(cfg$lon), length(cfg$lat))
if (!is.list(e$run_metadata) || !is.numeric(e$run_metadata$seed) ||
    !identical(as.numeric(e$run_metadata$seed), as.numeric(cfg$seed)))
  stop("Saved ERA5 metadata must record the configured seed before its main solution can be reused.")
stopifnot(length(models) >= 3L, models[1] == "ERA5", !anyDuplicated(models),
          is.numeric(e$year_present), is.numeric(e$year_future),
          identical(as.numeric(e$year_present), as.numeric(cfg$calibration_years)),
          identical(as.numeric(e$year_future), as.numeric(cfg$era5_target_years)),
          identical(e$reference_name, "ERA5"), identical(e$model_names, expected_models),
          identical(e$variables, cfg$variables), is.numeric(e$nbins3d),
          identical(as.numeric(e$nbins3d), as.numeric(cfg$nbins3d)),
          identical(e$lon, cfg$lon), identical(e$lat, cfg$lat),
          is.numeric(e$smooth_cost), identical(as.numeric(e$smooth_cost), as.numeric(cfg$lambda)),
          identical(dim(e$h_dist_pres), c(grid, length(expected_models))),
          identical(dim(e$pdf3_models_fut), c(grid, as.integer(cfg$nbins3d^3L), length(expected_models))),
          is.numeric(e$h_dist_pres), all(is.finite(e$h_dist_pres)), all(e$h_dist_pres >= 0))
validate_target_pmfs(e$pdf3_models_fut)
inputs <- list(pdf = e$pdf3_models_fut, h_dist = e$h_dist_pres,
               lon = e$lon, lat = e$lat, model_names = e$model_names,
               existing = e$GC_result, metadata = e$run_metadata)
# Verify the saved optimizer has the same numerical implementation as the
# paper-owned function, including seed-1 initialization and swap(-1).
saved_optimizer <- e$GraphCutHellinger_nD_lat
saved_text <- paste(deparse(body(saved_optimizer)), collapse = "\n")
current_text <- paste(deparse(body(GraphCutHellinger_nD_lat)), collapse = "\n")
if (!identical(saved_text, current_text)) {
  # Saved workspaces may include the paper-2 extension. The owned deterministic
  # implementation has patch-level numerical comparisons; no global equivalence
  # is inferred from the textual comparison alone.
  message("Saved optimizer differs; using the paper-owned deterministic implementation. See tests/test_production_fixture.R for the tested numerical comparison.")
}
optimizer <- GraphCutHellinger_nD_lat
rm(e, saved_optimizer, saved_text, current_text)
gc()
# Path/size/mtime alone cannot establish content identity. Hash once, then bind
# the numerical implementation/settings as well as the saved experiment bytes.
input_signature <- list(schema_version = 2L,
  workspace = normalizePath(workspace, winslash = "/", mustWork = TRUE),
  workspace_sha256 = mvgc_sha256(workspace), bytes = unname(file.info(workspace)$size),
  calibration = cfg$calibration_years, projection = cfg$era5_target_years,
  seed = cfg$seed, models = inputs$model_names, lon = inputs$lon, lat = inputs$lat,
  nbins3d = cfg$nbins3d, weight_data = 1, main_lambda = cfg$lambda,
  lambda_grid = all_lambdas, gco_commit = cfg$solver_commit,
  input_provenance = inputs$metadata$input_provenance,
  source_sha256 = scientific_code_fingerprint(c("R/graph_cut/solver.R", "R/graph_cut/energy.R",
    "R/utils/runtime.R", "experiments/sensitivity/spatial_regularisation.R")),
  r_version = as.character(getRversion()), platform = R.version$platform,
  package_versions = vapply(c("Rcpp", "gcoWrapR", "RcppXPtrUtils"), function(x)
    as.character(utils::packageVersion(x)), ""))
require_packages("openssl")
identity_hash <- as.character(openssl::sha256(serialize(input_signature, NULL, version = 2)))
checkpoint_dir <- file.path(out, "checkpoints", paste0("identity_", substr(identity_hash, 1, 24)))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
solutions <- setNames(vector("list", length(lambdas)), vapply(lambdas, key, ""))
for (lambda in lambdas) {
  path <- checkpoint_file(lambda)
  if (file.exists(path)) {
    result <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!spatial_checkpoint_is_complete(result, lambda, input_signature, grid, length(inputs$model_names)))
      stop("Invalid or incomplete spatial checkpoint; inspect it before resuming: ", path)
    cat("Reusing checkpoint lambda =", lambda, "\n")
  } else {
    started <- Sys.time()
    cat("Starting lambda =", lambda, "at", format(started), "\n")
    if (lambda == cfg$lambda) {
      gc_result <- inputs$existing
      origin <- "loaded final ERA5 workspace"
    } else {
      gc_result <- optimizer(pdf_models_future = inputs$pdf,
                             h_dist = inputs$h_dist, weight_data = 1,
                             weight_smooth = lambda, nBins = cfg$nbins3d^3L,
                             lat = inputs$lat, seed = cfg$seed,
                             verbose = TRUE, rebuild = FALSE)
      origin <- "new graph-cut computation"
    }
    result <- list(lambda = lambda, result = gc_result, origin = origin, state = "complete",
                   input_signature = input_signature,
                   elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    if (!spatial_checkpoint_is_complete(result, lambda, input_signature, grid, length(inputs$model_names)))
      stop("Invalid graph-cut solution; no successful checkpoint written for lambda = ", lambda)
    # Disjoint workers never expose a partially written checkpoint to a reader.
    temporary <- tempfile("solution_", tmpdir = checkpoint_dir, fileext = ".rds")
    saveRDS(result, temporary)
    if (!file.rename(temporary, path)) stop("Could not publish spatial checkpoint: ", path)
    cat("Saved lambda =", lambda, "after", round(result$elapsed_seconds, 1), "seconds\n")
  }
  solutions[[key(lambda)]] <- result
}
if (!all(file.exists(vapply(all_lambdas, checkpoint_file, "")))) {
  cat("Requested solutions checkpointed; other lambdas still pending.\n")
  quit(save = "no", status = 0)
}
solutions <- setNames(lapply(all_lambdas, function(x) readRDS(checkpoint_file(x))),
                      vapply(all_lambdas, key, ""))
if (!all(vapply(seq_along(all_lambdas), function(i)
    spatial_checkpoint_is_complete(solutions[[i]], all_lambdas[i], input_signature, grid, length(inputs$model_names)), logical(1))))
  stop("Cannot assemble final bundle: a spatial checkpoint is incomplete or incompatible.")
bundle <- list(solutions = solutions, lon = inputs$lon, lat = inputs$lat,
               model_names = inputs$model_names, input_signature = input_signature,
               metadata = inputs$metadata, state = "complete")
temporary_bundle <- tempfile("bundle_", tmpdir = out, fileext = ".rds")
saveRDS(bundle, temporary_bundle)
if (!file.rename(temporary_bundle, file.path(out, "lambda_label_maps_compact.rds")))
  stop("Could not atomically publish the complete spatial bundle.")
cat("Complete: compact label maps saved in", out, "\n")

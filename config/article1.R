# Scientific defaults are explicit and shared across the publication entry points.
#' Read the shared Article 1 scientific and runtime configuration
#'
#' Scientific dates, model-list paths, 8/16/32 resolutions, lambda/threshold grids
#' and seed 1 are explicit in this function. Only input registry, data/output/
#' workspace roots and support paths use MVGC_* environment overrides. Loading
#' configuration does not prove external inputs or supports are ready for a rerun.
#' @return List of scientific parameters, model-list/support paths, data roots
#'   and the required custom solver repository/commit. No files are opened or
#'   directories created. Support paths must be supplied explicitly;
#'   no historical support is silently substituted.
#' @examples
#' stopifnot(article1_config()$nbins3d == 8L, article1_config()$seed == 1L)
article1_config <- function() {
  list(
    variables = c("pr", "tas", "psl"), lon = -180:179, lat = -90:90,
    calibration_years = 1951:1975, era5_target_years = 1999:2023,
    periods = data.frame(
      period_id = c("1976_2000", "2001_2025", "2026_2050", "2051_2075", "2076_2100"),
      label = c("1976-2000", "2001-2025", "2026-2050", "2051-2075", "2076-2100"),
      start_year = c(1976L, 2001L, 2026L, 2051L, 2076L),
      end_year = c(2000L, 2025L, 2050L, 2075L, 2100L)),
    nbins3d = 8L, nbins2d = 16L, nbins1d = 32L,
    lambda = .125, era5_lambda_sensitivity = c(0, .025, .0625, .125, .25, .5, 1, 2),
    pme_lambdas = c(.125, .25, .5, 1),
    era5_ldr_thresholds = c(.01, .05, .10),
    pme_ldr_thresholds = c(.005, .01, .025, .05, .10),
    season_start_md = c(1, 1), season_end_md = c(12, 31),
    seed = 1L, strict_time_checks = TRUE, season_selection_mode = "auto",
    precipitation_preprocessing = "nonnegative_then_log1p_v1",
    time_decode_mode = "auto", requested_year_policy = "warn",
    input_registry = Sys.getenv("MVGC_INPUT_REGISTRY", ""),
    data_root = Sys.getenv("MVGC_DATA_ROOT", if(nzchar(Sys.getenv("MVGC_INPUT_REGISTRY", "")))
      dirname(Sys.getenv("MVGC_INPUT_REGISTRY")) else ""),
    output_root = Sys.getenv("MVGC_OUTPUT_ROOT", "outputs"),
    source_workspace = Sys.getenv("MVGC_SOURCE_WORKSPACE", ""),
    era5_support = Sys.getenv("MVGC_ERA5_SUPPORT", ""),
    pme_support = Sys.getenv("MVGC_PME_SUPPORT", ""),
    era5_models = "config/era5_models.txt", pme_models = "config/pme_models.txt",
    solver_repository = "schmutzlucas/gcoWrapR",
    solver_commit = "6b655d5feff7e7912d2ac0a7de95d0f8edd6d434")
}

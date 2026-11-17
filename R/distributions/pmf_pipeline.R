#' Build one model's marginal and joint PMFs using the reference histograms
#'
#' Internal worker for the public multi-model entry point. Reads all selected
#' variables for both periods before looping over cells. Precipitation uses the
#' shared nonnegative QC then log1p; precipitation means use corrected physical
#' observations before transformation. Missing rows are omitted separately for each marginal/joint
#' diagnostic. Empty local samples retain zero arrays, not valid PMFs.
#' @param model Single prepared model identifier.
#' @param variables Exactly c("pr", "tas", "psl"), in this order.
#' @param data_dir Prepared data hierarchy root, unless registry resolution is active.
#' @param year_present,year_future Requested calibration/evaluation year vectors.
#' @param lon,lat Ordered longitude/latitude coordinate vectors.
#' @param range_var Common local support array [lon,lat,variable,min/max], with
#'   precipitation bounds already on the log(1+pr) scale.
#' @param season_start_md,season_end_md Inclusive non-leap reference month/day boundaries.
#' @param pdf2_names Pair identifiers pr_tas, pr_psl and tas_psl.
#' @param nbins3d,nbins2d,nbins1d Positive integer bins per axis for joint,
#'   bivariate and univariate histograms respectively, normally 8/16/32.
#' @param strict_time_checks Whether to enforce legacy daily-time assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical time decoding.
#' @param requested_year_policy Warn or stop for incomplete/unavailable years.
#' @param diagnostics Whether detailed extraction diagnostics may be printed.
#' @param verbose Whether to print progress and, if diagnostics is TRUE, tables.
#' @param compute_hist_nd_fun Reference-compatible count callback with data,
#'   range_var and nbins arguments.
#' @param profile Compatibility profiling flag; timing is collected either way.
#' @param spatial_block_size Accepted for backend compatibility; ignored by this worker.
#' @return List of present/future three-dimensional PMF arrays, named lists of
#'   bivariate/univariate arrays and mean matrices, plus diagnostics/timing/performance.
#'   The first variable varies fastest in flattened histogram states.
compute_model_pdfs_fractional_legacy <- function(
  model,
  variables,
  data_dir,
  year_present,
  year_future,
  lon,
  lat,
  range_var,
  season_start_md,
  season_end_md,
  pdf2_names,
  nbins3d,
  nbins2d,
  nbins1d,
  strict_time_checks,
  season_selection_mode,
  time_decode_mode,
  requested_year_policy,
  diagnostics,
  verbose,
  compute_hist_nd_fun,
  profile = FALSE,
  spatial_block_size = NULL
) {
  model_started <- proc.time()[["elapsed"]]
  if (isTRUE(verbose)) {
    cat(format(Sys.time(), "%H:%M:%S"), "- reading", model, "\n")
  }

  var_present <- vector("list", length(variables))
  var_future <- vector("list", length(variables))
  var_diag <- vector("list", length(variables))
  names(var_present) <- variables
  names(var_future) <- variables
  names(var_diag) <- variables

  selected_present_ref <- NULL
  selected_future_ref <- NULL

  required_years <- sort(unique(c(year_present, year_future)))

  transform_seconds <- 0
  for (v in seq_along(variables)) {
    variable <- variables[v]

    file_path <- find_matching_netcdf_file(
      data_dir = data_dir,
      model = model,
      variable = variable,
      required_years = required_years
    )

    extraction <- extract_fractional_variable_data(
      file_path = file_path,
      variable = variable,
      lon = lon,
      lat = lat,
      year_present = year_present,
      year_future = year_future,
      season_start_md = season_start_md,
      season_end_md = season_end_md,
      strict_time_checks = strict_time_checks,
      season_selection_mode = season_selection_mode,
      time_decode_mode = time_decode_mode,
      requested_year_policy = requested_year_policy
    )

    if (is.null(selected_present_ref)) {
      selected_present_ref <- extraction$selected_present_dates
      selected_future_ref <- extraction$selected_future_dates
    } else {
      if (!identical(selected_present_ref, extraction$selected_present_dates)) {
        stop(
          "Selected PRESENT calendar dates differ across variables for model ", model,
          ". This indicates inconsistent time axes between variable files."
        )
      }
      if (!identical(selected_future_ref, extraction$selected_future_dates)) {
        stop(
          "Selected FUTURE calendar dates differ across variables for model ", model,
          ". This indicates inconsistent time axes between variable files."
        )
      }
    }

    data_present <- extraction$data_present
    data_future <- extraction$data_future

    if (variable == "pr") {
      transform_started <- proc.time()[["elapsed"]]
      pr_means_present <- precipitation_mean_grid(data_present)
      pr_means_future <- precipitation_mean_grid(data_future)
      data_present <- precipitation_transform(data_present)
      data_future <- precipitation_transform(data_future)
      transform_seconds <- transform_seconds + proc.time()[["elapsed"]] - transform_started
    }

    var_present[[variable]] <- data_present
    var_future[[variable]] <- data_future
    var_diag[[variable]] <- extraction$diagnostics

    gc(verbose = FALSE)
  }

  calendar_set <- unique(vapply(var_diag, function(x) x$calendar, character(1)))
  if (length(calendar_set) > 1L) {
    stop(
      "Calendar mismatch across variables for model ", model,
      ": ", paste(calendar_set, collapse = ", "),
      "."
    )
  }

  nlon <- length(lon)
  nlat <- length(lat)

  pdf3_present <- array(0, dim = c(nlon, nlat, nbins3d^3))
  pdf3_future <- array(0, dim = c(nlon, nlat, nbins3d^3))

  pdf2_present <- lapply(pdf2_names, function(x) array(0, dim = c(nlon, nlat, nbins2d^2)))
  pdf2_future <- lapply(pdf2_names, function(x) array(0, dim = c(nlon, nlat, nbins2d^2)))
  names(pdf2_present) <- pdf2_names
  names(pdf2_future) <- pdf2_names

  pdf1_present <- lapply(variables, function(x) array(0, dim = c(nlon, nlat, nbins1d)))
  pdf1_future <- lapply(variables, function(x) array(0, dim = c(nlon, nlat, nbins1d)))
  names(pdf1_present) <- variables
  names(pdf1_future) <- variables

  mean_present <- lapply(variables, function(x) matrix(NA_real_, nrow = nlon, ncol = nlat))
  mean_future <- lapply(variables, function(x) matrix(NA_real_, nrow = nlon, ncol = nlat))
  names(mean_present) <- variables
  names(mean_future) <- variables

  histogram_started <- proc.time()[["elapsed"]]
  for (i in seq_len(nlon)) {
    for (j in seq_len(nlat)) {
      pr_present <- as.numeric(var_present$pr[i, j, ])
      tas_present <- as.numeric(var_present$tas[i, j, ])
      psl_present <- as.numeric(var_present$psl[i, j, ])

      pr_future <- as.numeric(var_future$pr[i, j, ])
      tas_future <- as.numeric(var_future$tas[i, j, ])
      psl_future <- as.numeric(var_future$psl[i, j, ])

      range_pr <- matrix(range_var[i, j, 1, ], ncol = 2, byrow = TRUE)
      range_tas <- matrix(range_var[i, j, 2, ], ncol = 2, byrow = TRUE)
      range_psl <- matrix(range_var[i, j, 3, ], ncol = 2, byrow = TRUE)

      range_pr_tas <- rbind(range_var[i, j, 1, ], range_var[i, j, 2, ])
      range_pr_psl <- rbind(range_var[i, j, 1, ], range_var[i, j, 3, ])
      range_tas_psl <- rbind(range_var[i, j, 2, ], range_var[i, j, 3, ])
      range_3d <- rbind(range_var[i, j, 1, ], range_var[i, j, 2, ], range_var[i, j, 3, ])

      data_present <- cbind(pr_present, tas_present, psl_present)
      data_future <- cbind(pr_future, tas_future, psl_future)

      pdf3_present[i, j, ] <- safe_hist_counts(data_present, range_3d, nbins3d, compute_hist_nd_fun)
      pdf3_future[i, j, ] <- safe_hist_counts(data_future, range_3d, nbins3d, compute_hist_nd_fun)

      pdf2_present$pr_tas[i, j, ] <- safe_hist_counts(data_present[, 1:2, drop = FALSE], range_pr_tas, nbins2d, compute_hist_nd_fun)
      pdf2_future$pr_tas[i, j, ] <- safe_hist_counts(data_future[, 1:2, drop = FALSE], range_pr_tas, nbins2d, compute_hist_nd_fun)

      pdf2_present$pr_psl[i, j, ] <- safe_hist_counts(data_present[, c(1, 3), drop = FALSE], range_pr_psl, nbins2d, compute_hist_nd_fun)
      pdf2_future$pr_psl[i, j, ] <- safe_hist_counts(data_future[, c(1, 3), drop = FALSE], range_pr_psl, nbins2d, compute_hist_nd_fun)

      pdf2_present$tas_psl[i, j, ] <- safe_hist_counts(data_present[, 2:3, drop = FALSE], range_tas_psl, nbins2d, compute_hist_nd_fun)
      pdf2_future$tas_psl[i, j, ] <- safe_hist_counts(data_future[, 2:3, drop = FALSE], range_tas_psl, nbins2d, compute_hist_nd_fun)

      pdf1_present$pr[i, j, ] <- safe_hist_counts(matrix(pr_present, ncol = 1), range_pr, nbins1d, compute_hist_nd_fun)
      pdf1_future$pr[i, j, ] <- safe_hist_counts(matrix(pr_future, ncol = 1), range_pr, nbins1d, compute_hist_nd_fun)

      pdf1_present$tas[i, j, ] <- safe_hist_counts(matrix(tas_present, ncol = 1), range_tas, nbins1d, compute_hist_nd_fun)
      pdf1_future$tas[i, j, ] <- safe_hist_counts(matrix(tas_future, ncol = 1), range_tas, nbins1d, compute_hist_nd_fun)

      pdf1_present$psl[i, j, ] <- safe_hist_counts(matrix(psl_present, ncol = 1), range_psl, nbins1d, compute_hist_nd_fun)
      pdf1_future$psl[i, j, ] <- safe_hist_counts(matrix(psl_future, ncol = 1), range_psl, nbins1d, compute_hist_nd_fun)

      mean_present$pr[i, j] <- pr_means_present[i, j]
      mean_future$pr[i, j] <- pr_means_future[i, j]
      mean_present$tas[i, j] <- safe_mean_value(tas_present)
      mean_future$tas[i, j] <- safe_mean_value(tas_future)
      mean_present$psl[i, j] <- safe_mean_value(psl_present)
      mean_future$psl[i, j] <- safe_mean_value(psl_future)
    }
  }
  histogram_seconds <- proc.time()[["elapsed"]] - histogram_started

  normalization_started <- proc.time()[["elapsed"]]
  pdf3_present <- normalize_histogram_array(pdf3_present)
  pdf3_future <- normalize_histogram_array(pdf3_future)

  for (name in pdf2_names) {
    pdf2_present[[name]] <- normalize_histogram_array(pdf2_present[[name]])
    pdf2_future[[name]] <- normalize_histogram_array(pdf2_future[[name]])
  }

  for (name in variables) {
    pdf1_present[[name]] <- normalize_histogram_array(pdf1_present[[name]])
    pdf1_future[[name]] <- normalize_histogram_array(pdf1_future[[name]])
  }
  normalization_seconds <- proc.time()[["elapsed"]] - normalization_started

  model_diag <- list(
    model = model,
    calendar = calendar_set[1],
    variables = var_diag
  )

  if (isTRUE(diagnostics) && isTRUE(verbose)) {
    cat("Diagnostics summary for model", model, "\n")
    for (v in variables) {
      print_fractional_selection_diagnostics(var_diag[[v]])
    }
  }

  timing <- c(
    metadata_and_time_axis = sum(vapply(var_diag, function(x) x$timing[["metadata_and_time_axis"]], numeric(1))),
    netcdf_variable_reads = sum(vapply(var_diag, function(x) x$timing[["netcdf_variable_reads"]], numeric(1))),
    precipitation_transform = transform_seconds,
    legacy_histogram_and_assembly = histogram_seconds,
    normalization = normalization_seconds,
    total = proc.time()[["elapsed"]] - model_started
  )

  list(
    pdf3_present = pdf3_present,
    pdf3_future = pdf3_future,
    pdf2_present = pdf2_present,
    pdf2_future = pdf2_future,
    pdf1_present = pdf1_present,
    pdf1_future = pdf1_future,
    mean_present = mean_present,
    mean_future = mean_future,
    diagnostics = model_diag,
    timing = timing,
    performance = list(
      blocked_reads = FALSE,
      latitude_block_size = nlat,
      latitude_blocks = 1L,
      approximate_peak_raw_input_bytes = as.numeric(
        sum(vapply(c(var_present, var_future), object.size, numeric(1)))
      )
    )
  )
}

#' Build one model's PMFs using shared nested bin indices and latitude blocks
#'
#' Uses exactly the supplied support, never estimates it per reference. Fast
#' histogram grids are nested where their resolutions permit; floating-edge
#' snapping follows equal_width_bin_index(). Date keys across variables must agree.
#' A latitude block bounds raw input memory, not the final retained PMF arrays.
#' @param model Single prepared model identifier.
#' @param variables Exactly c("pr", "tas", "psl"), in this order.
#' @param data_dir Prepared data root unless an input registry is active.
#' @param year_present,year_future Requested calibration/evaluation year vectors.
#' @param lon,lat Ordered coordinate vectors, matched without interpolation.
#' @param range_var Common [lon,lat,3,2] support; pr is on log(1+pr) scale.
#' @param season_start_md,season_end_md Inclusive non-leap reference month/day boundaries.
#' @param pdf2_names Pair identifiers pr_tas, pr_psl and tas_psl.
#' @param nbins3d,nbins2d,nbins1d Bins per variable at each diagnostic dimensionality,
#'   normally 8/16/32. No statistical discretisation choice is made here.
#' @param strict_time_checks Whether to enforce legacy daily-time assumptions.
#' @param season_selection_mode Auto, exact_month_day, or fractional selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical decoding.
#' @param requested_year_policy Warn or stop for incomplete/unavailable years.
#' @param diagnostics Whether detailed extraction diagnostics may be printed.
#' @param verbose Whether to print progress and requested diagnostic tables.
#' @param compute_hist_nd_fun Unused compatibility argument; fast count helpers are used.
#' @param profile Compatibility profiling flag; timing is collected either way.
#' @param spatial_block_size Maximum requested latitude positions per raw-input
#'   block. NULL or a value at least the latitude count reads complete selected cubes.
#' @return List matching the legacy worker's PMFs/means/diagnostics, with additional
#'   block/bin-index timing and approximate raw-input peak-memory statistics.
#'   PMFs omit incomplete rows per diagnostic; empty cells retain legacy zero arrays.
#'   pr uses shared nonnegative QC then log1p. Its means use physical corrected
#'   observations, not a logarithm/exponential round trip.
compute_model_pdfs_fractional_optimized <- function(
  model,
  variables,
  data_dir,
  year_present,
  year_future,
  lon,
  lat,
  range_var,
  season_start_md,
  season_end_md,
  pdf2_names,
  nbins3d,
  nbins2d,
  nbins1d,
  strict_time_checks,
  season_selection_mode,
  time_decode_mode,
  requested_year_policy,
  diagnostics,
  verbose,
  compute_hist_nd_fun = NULL,
  profile = FALSE,
  spatial_block_size = NULL
) {
  if (!exists("compute_multires_hist_counts_fast", mode = "function")) {
    stop("compute_multires_hist_counts_fast() is unavailable. Source R/distributions/histogram_fast.R.")
  }
  model_started <- proc.time()[["elapsed"]]
  if (isTRUE(verbose)) cat(format(Sys.time(), "%H:%M:%S"), "- reading", model, "\n")

  nlon <- length(lon)
  nlat <- length(lat)
  if (!is.null(spatial_block_size)) {
    if (length(spatial_block_size) != 1L || is.na(spatial_block_size) || spatial_block_size < 1L) {
      stop("spatial_block_size must be NULL or a positive integer.")
    }
    spatial_block_size <- as.integer(spatial_block_size)
  }
  use_blocked_reads <- !is.null(spatial_block_size) && spatial_block_size < nlat

  var_present <- setNames(vector("list", length(variables)), variables)
  var_future <- setNames(vector("list", length(variables)), variables)
  var_plan <- setNames(vector("list", length(variables)), variables)
  var_diag <- setNames(vector("list", length(variables)), variables)
  selected_present_ref <- NULL
  selected_future_ref <- NULL
  required_years <- sort(unique(c(year_present, year_future)))
  transform_seconds <- 0

  for (variable in variables) {
    file_path <- find_matching_netcdf_file(data_dir, model, variable, required_years)
    extraction <- extract_fractional_variable_data(
      file_path = file_path,
      variable = variable,
      lon = lon,
      lat = lat,
      year_present = year_present,
      year_future = year_future,
      season_start_md = season_start_md,
      season_end_md = season_end_md,
      strict_time_checks = strict_time_checks,
      season_selection_mode = season_selection_mode,
      time_decode_mode = time_decode_mode,
      requested_year_policy = requested_year_policy,
      read_data = !use_blocked_reads
    )

    if (is.null(selected_present_ref)) {
      selected_present_ref <- extraction$selected_present_dates
      selected_future_ref <- extraction$selected_future_dates
    } else {
      if (!identical(selected_present_ref, extraction$selected_present_dates)) {
        stop("Selected PRESENT calendar dates differ across variables for model ", model, ".")
      }
      if (!identical(selected_future_ref, extraction$selected_future_dates)) {
        stop("Selected FUTURE calendar dates differ across variables for model ", model, ".")
      }
    }

    var_plan[[variable]] <- extraction$plan
    if (!use_blocked_reads) {
      data_present <- extraction$data_present
      data_future <- extraction$data_future
      if (identical(variable, "pr")) {
        transform_started <- proc.time()[["elapsed"]]
        pr_means_present <- precipitation_mean_grid(data_present)
        pr_means_future <- precipitation_mean_grid(data_future)
        data_present <- precipitation_transform(data_present)
        data_future <- precipitation_transform(data_future)
        transform_seconds <- transform_seconds + proc.time()[["elapsed"]] - transform_started
      }
      var_present[[variable]] <- data_present
      var_future[[variable]] <- data_future
    }
    var_diag[[variable]] <- extraction$diagnostics
    rm(extraction)
    gc(verbose = FALSE)
  }

  calendar_set <- unique(vapply(var_diag, function(x) x$calendar, character(1)))
  if (length(calendar_set) > 1L) {
    stop("Calendar mismatch across variables for model ", model, ": ", paste(calendar_set, collapse = ", "), ".")
  }

  pdf3_present <- array(0, dim = c(nlon, nlat, nbins3d^3))
  pdf3_future <- array(0, dim = c(nlon, nlat, nbins3d^3))
  pdf2_present <- setNames(lapply(pdf2_names, function(x) array(0, c(nlon, nlat, nbins2d^2))), pdf2_names)
  pdf2_future <- setNames(lapply(pdf2_names, function(x) array(0, c(nlon, nlat, nbins2d^2))), pdf2_names)
  pdf1_present <- setNames(lapply(variables, function(x) array(0, c(nlon, nlat, nbins1d))), variables)
  pdf1_future <- setNames(lapply(variables, function(x) array(0, c(nlon, nlat, nbins1d))), variables)
  mean_present <- setNames(lapply(variables, function(x) matrix(NA_real_, nlon, nlat)), variables)
  mean_future <- setNames(lapply(variables, function(x) matrix(NA_real_, nlon, nlat)), variables)

  bin_index_seconds <- 0
  tabulation_seconds <- 0
  block_read_seconds <- 0
  peak_raw_input_bytes <- if (use_blocked_reads) 0 else {
    sum(vapply(c(var_present, var_future), object.size, numeric(1)))
  }
  lat_blocks <- if (use_blocked_reads) {
    split(seq_len(nlat), ceiling(seq_len(nlat) / spatial_block_size))
  } else {
    list(seq_len(nlat))
  }

  assembly_started <- proc.time()[["elapsed"]]
  for (lat_block in lat_blocks) {
    if (use_blocked_reads) {
      block_present <- setNames(vector("list", length(variables)), variables)
      block_future <- setNames(vector("list", length(variables)), variables)
      read_started <- proc.time()[["elapsed"]]
      for (variable in variables) {
        block_data <- read_fractional_variable_plan_block(var_plan[[variable]], lat_block)
        block_present[[variable]] <- block_data$present
        block_future[[variable]] <- block_data$future
      }
      block_read_seconds <- block_read_seconds + proc.time()[["elapsed"]] - read_started
      peak_raw_input_bytes <- max(
        peak_raw_input_bytes,
        sum(vapply(c(block_present, block_future), object.size, numeric(1)))
      )
      transform_started <- proc.time()[["elapsed"]]
      pr_mean_present_block <- precipitation_mean_grid(block_present$pr)
      pr_mean_future_block <- precipitation_mean_grid(block_future$pr)
      block_present$pr <- precipitation_transform(block_present$pr)
      block_future$pr <- precipitation_transform(block_future$pr)
      transform_seconds <- transform_seconds + proc.time()[["elapsed"]] - transform_started
    } else {
      block_present <- var_present
      block_future <- var_future
      pr_mean_present_block <- pr_means_present
      pr_mean_future_block <- pr_means_future
    }

    for (i in seq_len(nlon)) {
      for (j_local in seq_along(lat_block)) {
        j <- lat_block[j_local]
      range_3d <- rbind(
        range_var[i, j, 1, ],
        range_var[i, j, 2, ],
        range_var[i, j, 3, ]
      )
      present <- compute_multires_hist_counts_fast(
        as.numeric(block_present$pr[i, j_local, ]),
        as.numeric(block_present$tas[i, j_local, ]),
        as.numeric(block_present$psl[i, j_local, ]),
        range_3d, nbins3d, nbins2d, nbins1d,
        profile = profile
      )
      future <- compute_multires_hist_counts_fast(
        as.numeric(block_future$pr[i, j_local, ]),
        as.numeric(block_future$tas[i, j_local, ]),
        as.numeric(block_future$psl[i, j_local, ]),
        range_3d, nbins3d, nbins2d, nbins1d,
        profile = profile
      )

      pdf3_present[i, j, ] <- present$pdf3
      pdf3_future[i, j, ] <- future$pdf3
      for (name in pdf2_names) {
        pdf2_present[[name]][i, j, ] <- present$pdf2[[name]]
        pdf2_future[[name]][i, j, ] <- future$pdf2[[name]]
      }
      for (name in variables) {
        pdf1_present[[name]][i, j, ] <- present$pdf1[[name]]
        pdf1_future[[name]][i, j, ] <- future$pdf1[[name]]
      }
      if (isTRUE(profile)) {
        bin_index_seconds <- bin_index_seconds + present$timing[["bin_index"]] + future$timing[["bin_index"]]
        tabulation_seconds <- tabulation_seconds + present$timing[["tabulation"]] + future$timing[["tabulation"]]
      }

      pr_present <- as.numeric(block_present$pr[i, j_local, ])
      pr_future <- as.numeric(block_future$pr[i, j_local, ])
      mean_present$pr[i, j] <- pr_mean_present_block[i, j_local]
      mean_future$pr[i, j] <- pr_mean_future_block[i, j_local]
      mean_present$tas[i, j] <- safe_mean_value(block_present$tas[i, j_local, ])
      mean_future$tas[i, j] <- safe_mean_value(block_future$tas[i, j_local, ])
      mean_present$psl[i, j] <- safe_mean_value(block_present$psl[i, j_local, ])
      mean_future$psl[i, j] <- safe_mean_value(block_future$psl[i, j_local, ])
      }
    }

    if (use_blocked_reads) {
      rm(block_present, block_future, block_data)
      gc(verbose = FALSE)
    }
  }
  assembly_elapsed <- proc.time()[["elapsed"]] - assembly_started
  assembly_seconds <- max(
    0,
    assembly_elapsed - block_read_seconds - if (use_blocked_reads) transform_seconds else 0
  )

  normalization_started <- proc.time()[["elapsed"]]
  pdf3_present <- normalize_histogram_array(pdf3_present)
  pdf3_future <- normalize_histogram_array(pdf3_future)
  for (name in pdf2_names) {
    pdf2_present[[name]] <- normalize_histogram_array(pdf2_present[[name]])
    pdf2_future[[name]] <- normalize_histogram_array(pdf2_future[[name]])
  }
  for (name in variables) {
    pdf1_present[[name]] <- normalize_histogram_array(pdf1_present[[name]])
    pdf1_future[[name]] <- normalize_histogram_array(pdf1_future[[name]])
  }
  normalization_seconds <- proc.time()[["elapsed"]] - normalization_started

  model_diag <- list(model = model, calendar = calendar_set[1], variables = var_diag)
  if (isTRUE(diagnostics) && isTRUE(verbose)) {
    cat("Diagnostics summary for model", model, "\n")
    for (v in variables) print_fractional_selection_diagnostics(var_diag[[v]])
  }

  timing <- c(
    metadata_and_time_axis = sum(vapply(var_diag, function(x) x$timing[["metadata_and_time_axis"]], numeric(1))),
    netcdf_variable_reads = sum(vapply(var_diag, function(x) x$timing[["netcdf_variable_reads"]], numeric(1))) + block_read_seconds,
    precipitation_transform = transform_seconds,
    bin_index = bin_index_seconds,
    histogram_tabulation = tabulation_seconds,
    histogram_and_assembly = assembly_seconds,
    normalization = normalization_seconds,
    total = proc.time()[["elapsed"]] - model_started
  )

  list(
    pdf3_present = pdf3_present,
    pdf3_future = pdf3_future,
    pdf2_present = pdf2_present,
    pdf2_future = pdf2_future,
    pdf1_present = pdf1_present,
    pdf1_future = pdf1_future,
    mean_present = mean_present,
    mean_future = mean_future,
    diagnostics = model_diag,
    timing = timing,
    performance = list(
      blocked_reads = use_blocked_reads,
      latitude_block_size = if (use_blocked_reads) spatial_block_size else nlat,
      latitude_blocks = length(lat_blocks),
      approximate_peak_raw_input_bytes = as.numeric(peak_raw_input_bytes)
    )
  )
}

#' Dispatch a single-model PMF computation to the selected backend
#' @param ... Named arguments forwarded to the corresponding single-model worker.
#' @param backend Optimized shared-index/block implementation or legacy reference.
#' @param profile Forwarded compatibility profiling flag.
#' @return The selected worker's list of PMFs, means and diagnostics. This dispatcher
#'   does not choose periods, supports, models, or a different histogram resolution.
compute_model_pdfs_fractional <- function(..., backend = c("optimized", "legacy"), profile = FALSE) {
  backend <- match.arg(backend)
  fun <- if (identical(backend, "optimized")) {
    compute_model_pdfs_fractional_optimized
  } else {
    compute_model_pdfs_fractional_legacy
  }
  fun(..., profile = profile)
}

#' Build calibration and evaluation PMFs on one fixed support for multiple models
#'
#' Validates common grid/support/model configuration, dispatches model workers,
#' then assembles arrays with a final model axis. This is a memory-intensive
#' computation: latitude blocking limits raw input cubes, not retained outputs.
#' Missing observations are excluded per marginal/joint sample; all-zero local
#' histograms remain zero under the legacy convention and must not be fed to a
#' solver as normalised PMFs. No new support or reference-specific support is made.
#' @param variables Exactly c("pr", "tas", "psl"), in that order.
#' @param model_names Distinct nonempty model identifiers in desired output order.
#' @param data_dir Prepared daily NetCDF root; registry resolution takes precedence.
#' @param year_present,year_future Nonempty finite integer calibration/evaluation years.
#' @param lon,lat Nonempty unique finite longitude/latitude coordinate vectors.
#' @param range_var Finite strictly increasing bounds [lon,lat,3,2], aligned with
#'   the requested grid; pr bounds must already represent log(1+pr).
#' @param season_start_md,season_end_md Inclusive month/day boundaries on a 365-day
#'   reference year. Helper defaults are seasonal; Article 1 entry points explicitly
#'   pass January 1 through December 31, so do not infer production dates from defaults.
#' @param nbins3d,nbins2d,nbins1d Positive integer bins per variable, default 8/16/32.
#' @param workers Positive integer model-worker count. Uses multicore on Unix and
#'   multisession elsewhere, restoring the previous future plan on exit.
#' @param strict_time_checks Whether to enforce the legacy daily-mapping assumptions.
#' @param diagnostics Whether detailed model extraction diagnostics may be printed.
#' @param verbose Whether workers print progress and requested diagnostics.
#' @param season_selection_mode Auto, exact_month_day, or fractional seasonal selection.
#' @param time_decode_mode Auto, helpers_cf, or empirical time decoding.
#' @param requested_year_policy Warn or stop for missing/incomplete requested years;
#'   warn preserves the original experiment behaviour, not an input repair policy.
#' @param pmf_backend Optimized shared-index/blocked worker or legacy reference worker.
#' @param profile Compatibility profiling flag forwarded to each model worker.
#' @param spatial_block_size Latitude positions per optimized raw-input block; NULL
#'   disables blocked reads. Defaults to eight, independently of the histogram bins.
#' @return List pdf3/pdf2/pdf1/mean, each with present/future outputs; arrays have
#'   longitude x latitude x bins x models, mean arrays omit the bin axis.
#'   Bivariate/univariate/mean entries are named lists. Also returns diagnostics,
#'   calendars, timing/performance lists by model, and the selected pmf_backend.
compute_nd_pdf_multi_seasonal <- function(
  variables,
  model_names,
  data_dir,
  year_present,
  year_future,
  lon,
  lat,
  range_var,
  season_start_md = c(4, 15),
  season_end_md = c(10, 14),
  nbins3d = 8,
  nbins2d = 16,
  nbins1d = 32,
  workers = 4,
  strict_time_checks = TRUE,
  diagnostics = FALSE,
  verbose = FALSE,
  season_selection_mode = c("auto", "exact_month_day", "fractional"),
  time_decode_mode = c("auto", "helpers_cf", "empirical"),
  requested_year_policy = c("warn", "stop"),
  pmf_backend = c("optimized", "legacy"),
  profile = FALSE,
  spatial_block_size = 8L
) {
  if (!exists("compute_histND", mode = "function")) {
    stop("compute_histND() is not available. Source R/distributions/histogram_reference.R before calling this function.")
  }

  if (!is.character(variables) || length(variables) != 3L || !all(variables == c("pr", "tas", "psl"))) {
    stop("variables must be exactly c('pr', 'tas', 'psl').")
  }
  if(!is.character(model_names)||!length(model_names)||anyNA(model_names)||
     any(!nzchar(model_names))||anyDuplicated(model_names))stop("model_names must be distinct nonempty names.")
  for(years in list(year_present,year_future))
    if(!is.numeric(years)||!length(years)||any(!is.finite(years))||any(years!=trunc(years)))
      stop("Experiment years must be finite integers.")
  for(grid in list(lon,lat))
    if(!is.numeric(grid)||!length(grid)||any(!is.finite(grid))||anyDuplicated(grid))
      stop("Longitude and latitude must be finite nonempty unique coordinate vectors.")
  if(!is.numeric(range_var)||!identical(as.integer(dim(range_var)),c(length(lon),length(lat),3L,2L))||
     any(!is.finite(range_var))||any(range_var[,,,2]<=range_var[,,,1]))
    stop("range_var must be finite [longitude, latitude, 3, 2] with strictly increasing bounds.")
  bins<-list(nbins3d,nbins2d,nbins1d)
  if(any(lengths(bins)!=1L)||!all(vapply(bins,is.numeric,TRUE))||
     any(!is.finite(unlist(bins)))||any(unlist(bins)<1)||any(unlist(bins)!=trunc(unlist(bins))))
    stop("Histogram resolutions must be single positive integers.")

  season_start_md <- validate_month_day(season_start_md, "season_start_md")
  season_end_md <- validate_month_day(season_end_md, "season_end_md")
  season_selection_mode <- match.arg(season_selection_mode)
  time_decode_mode <- match.arg(time_decode_mode)
  requested_year_policy <- match.arg(requested_year_policy)
  pmf_backend <- match.arg(pmf_backend)

  year_present <- sort(unique(as.integer(year_present)))
  year_future <- sort(unique(as.integer(year_future)))

  if (length(year_present) == 0L || length(year_future) == 0L) {
    stop("year_present and year_future must each contain at least one year.")
  }

  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) || workers < 1L || workers!=trunc(workers)) {
    stop("workers must be a positive integer.")
  }
  if (!is.null(spatial_block_size) &&
      (!is.numeric(spatial_block_size) || length(spatial_block_size) != 1L ||
       is.na(spatial_block_size) || spatial_block_size < 1L)) {
    stop("spatial_block_size must be NULL or a positive integer.")
  }
  if (!is.null(spatial_block_size)) spatial_block_size <- as.integer(spatial_block_size)

  n_models <- length(model_names)
  nlon <- length(lon)
  nlat <- length(lat)
  pdf2_names <- c("pr_tas", "pr_psl", "tas_psl")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }

  model_list <- future.apply::future_lapply(
    X = model_names,
    FUN = function(model_name) {
      compute_model_pdfs_fractional(
        model = model_name,
        variables = variables,
        data_dir = data_dir,
        year_present = year_present,
        year_future = year_future,
        lon = lon,
        lat = lat,
        range_var = range_var,
        season_start_md = season_start_md,
        season_end_md = season_end_md,
        pdf2_names = pdf2_names,
        nbins3d = nbins3d,
        nbins2d = nbins2d,
        nbins1d = nbins1d,
        strict_time_checks = strict_time_checks,
        season_selection_mode = season_selection_mode,
        time_decode_mode = time_decode_mode,
        requested_year_policy = requested_year_policy,
        diagnostics = diagnostics,
        verbose = verbose,
        compute_hist_nd_fun = compute_histND,
        backend = pmf_backend,
        profile = profile,
        spatial_block_size = spatial_block_size
      )
    },
    future.packages = c("ncdf4", "stats")
  )

  pdf3_present <- array(NA_real_, dim = c(nlon, nlat, nbins3d^3, n_models))
  pdf3_future <- array(NA_real_, dim = c(nlon, nlat, nbins3d^3, n_models))

  pdf2_present <- lapply(pdf2_names, function(x) array(NA_real_, dim = c(nlon, nlat, nbins2d^2, n_models)))
  pdf2_future <- lapply(pdf2_names, function(x) array(NA_real_, dim = c(nlon, nlat, nbins2d^2, n_models)))
  names(pdf2_present) <- pdf2_names
  names(pdf2_future) <- pdf2_names

  pdf1_present <- lapply(variables, function(x) array(NA_real_, dim = c(nlon, nlat, nbins1d, n_models)))
  pdf1_future <- lapply(variables, function(x) array(NA_real_, dim = c(nlon, nlat, nbins1d, n_models)))
  names(pdf1_present) <- variables
  names(pdf1_future) <- variables

  mean_present <- lapply(variables, function(x) array(NA_real_, dim = c(nlon, nlat, n_models)))
  mean_future <- lapply(variables, function(x) array(NA_real_, dim = c(nlon, nlat, n_models)))
  names(mean_present) <- variables
  names(mean_future) <- variables

  diagnostics_by_model <- vector("list", n_models)

  for (m in seq_len(n_models)) {
    pdf3_present[, , , m] <- model_list[[m]]$pdf3_present
    pdf3_future[, , , m] <- model_list[[m]]$pdf3_future

    for (name in pdf2_names) {
      pdf2_present[[name]][, , , m] <- model_list[[m]]$pdf2_present[[name]]
      pdf2_future[[name]][, , , m] <- model_list[[m]]$pdf2_future[[name]]
    }

    for (name in variables) {
      pdf1_present[[name]][, , , m] <- model_list[[m]]$pdf1_present[[name]]
      pdf1_future[[name]][, , , m] <- model_list[[m]]$pdf1_future[[name]]
      mean_present[[name]][, , m] <- model_list[[m]]$mean_present[[name]]
      mean_future[[name]][, , m] <- model_list[[m]]$mean_future[[name]]
    }

    diagnostics_by_model[[m]] <- model_list[[m]]$diagnostics
  }

  names(diagnostics_by_model) <- model_names

  calendars <- vapply(
    diagnostics_by_model,
    function(x) x$calendar,
    character(1)
  )
  names(calendars) <- model_names

  timing_by_model <- lapply(model_list, function(x) x$timing)
  names(timing_by_model) <- model_names
  performance_by_model <- lapply(model_list, function(x) x$performance)
  names(performance_by_model) <- model_names

  list(
    pdf3 = list(present = pdf3_present, future = pdf3_future),
    pdf2 = list(present = pdf2_present, future = pdf2_future),
    pdf1 = list(present = pdf1_present, future = pdf1_future),
    mean = list(present = mean_present, future = mean_future),
    diagnostics = diagnostics_by_model,
    calendars = calendars,
    timing = timing_by_model,
    performance = performance_by_model,
    pmf_backend = pmf_backend
  )
}

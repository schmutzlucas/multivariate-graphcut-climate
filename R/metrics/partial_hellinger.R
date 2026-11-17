#' Compute restricted Hellinger scores for vectors or global PMF arrays
#' @param pdf_ref Reference PMF vector, or array longitude x latitude x bins.
#' @param pdf_model Candidate vector; bins x models matrix for a vector reference;
#'   or longitude x latitude x bins array with an optional fourth model axis.
#' @param selected_indices One-based bin indices for a single location, or nested
#'   lists selected_indices[[longitude]][[latitude]] for a global array.
#' @return Scalar or model vector for vector inputs; longitude x latitude x models
#'   array for grid inputs, retaining a singleton model dimension. Empty selections
#'   return NA rather than silently returning a zero score.
#' @details Divides the square-root sum by sqrt(2) without renormalising selected
#'   masses. Selecting all bins gives ordinary Hellinger distance. Callers must
#'   supply common bins, valid masses and consistent selected indices.
#' @examples
#' compute_partial_hdist(c(0.2, 0.8), c(0.3, 0.7), 1:2)
compute_partial_hdist <- function(pdf_ref, pdf_model, selected_indices) {
  ref_dims <- dim(pdf_ref)
  model_dims <- dim(pdf_model)

  # Case 1: Single grid point, single model (vectors).
  if (is.null(ref_dims) && is.null(model_dims)) {
    if (length(selected_indices) > 0L) {
      return(sqrt(sum((sqrt(pdf_model[selected_indices]) - sqrt(pdf_ref[selected_indices]))^2)) / sqrt(2))
    }
    return(NA_real_)
  }

  # Case 2: Single grid point, multiple models (matrix: bins x models).
  if (is.null(ref_dims) && !is.null(model_dims) && length(model_dims) == 2L) {
    n_models <- model_dims[2]
    h_dist_partial <- rep(NA_real_, n_models)

    for (m in seq_len(n_models)) {
      if (length(selected_indices) > 0L) {
        h_dist_partial[m] <- sqrt(
          sum((sqrt(pdf_model[selected_indices, m]) - sqrt(pdf_ref[selected_indices]))^2)
        ) / sqrt(2)
      }
    }
    return(h_dist_partial)
  }

  # Case 3: Full grid (ref: lon x lat x bins, model: lon x lat x bins [x models]).
  if (is.null(ref_dims) || length(ref_dims) < 3L) {
    stop("pdf_ref must be either a vector or an array with dimensions lon x lat x bins.")
  }

  lon_size <- ref_dims[1]
  lat_size <- ref_dims[2]
  n_models <- if (!is.null(model_dims) && length(model_dims) == 4L) model_dims[4] else 1L

  h_dist_partial <- array(NA_real_, dim = c(lon_size, lat_size, n_models))

  for (j in seq_len(lat_size)) {
    for (i in seq_len(lon_size)) {
      selected_indices_grid <- selected_indices[[i]][[j]]

      if (length(selected_indices_grid) > 0L) {
        for (m in seq_len(n_models)) {
          # A singleton model axis is still a fourth array dimension in R.
          model_pdf <- if (length(model_dims) == 4L) {
            pdf_model[i, j, selected_indices_grid, m]
          } else {
            pdf_model[i, j, selected_indices_grid]
          }

          h_dist_partial[i, j, m] <- sqrt(
            sum((sqrt(model_pdf) - sqrt(pdf_ref[i, j, selected_indices_grid]))^2)
          ) / sqrt(2)
        }
      }
    }
  }

  h_dist_partial
}

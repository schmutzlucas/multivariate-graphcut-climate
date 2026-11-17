#' Validate target PMFs before handing their memory to native callbacks
#' @param pdf_models_future Numeric longitude x latitude x bins x models array.
#' @param tolerance Absolute tolerance for each local PMF's total mass.
#' @return TRUE invisibly; invalid masses or dimensions raise an R error.
#' @details Scans one latitude slab at a time to bound extra memory. Missing,
#'   negative and nonfinite masses are rejected before C++ square roots or graph
#'   optimisation can consume them. Valid PMFs are not modified or renormalised.
#' @examples
#' validate_target_pmfs(array(0.5, c(3, 2, 2, 2)))
validate_target_pmfs <- function(pdf_models_future,tolerance=1e-8) {
 dims<-dim(pdf_models_future)
 if(!is.numeric(pdf_models_future)||length(dims)!=4L||any(dims<1L)||
  !is.numeric(tolerance)||length(tolerance)!=1L||!is.finite(tolerance)||tolerance<0)
  stop("Expected a four-dimensional PMF array and a non-negative finite tolerance")
 for(j in seq_len(dims[2])) {
  slab<-pdf_models_future[,j,,,drop=FALSE]
  if(any(!is.finite(slab))||any(slab<0))stop("Target PMFs contain negative or nonfinite masses")
  probabilities<-matrix(aperm(slab,c(3,1,2,4)),nrow=dims[3])
  if(any(abs(colSums(probabilities)-1)>tolerance))stop("Target PMFs must sum to one in every cell/model")
 }
 invisible(TRUE)
}

#' Optimise latitude-weighted cylindrical graph-cut model labels
#'
#' Fit labels using calibration Hellinger unary costs and target-period PMF seam
#' costs. Article 1 uses the pinned custom gcoWrapR solver and swap(-1).
#' @param pdf_models_future Normalised target PMFs, dimensions longitude x
#'   latitude x joint bins x candidate models. Longitude must be periodic and
#'   ordered consistently; the native grid requires at least 3 x 2 cells and 2 labels.
#' @param h_dist Nonnegative calibration Hellinger distances, longitude x latitude x models.
#' @param weight_data Nonnegative scalar multiplying the unary term, normally 1.
#' @param weight_smooth Nonnegative scalar lambda multiplying the seam term.
#' @param nBins Total number of joint bins, normally 512 for 8 x 8 x 8,
#'   not the number of bins along one variable.
#' @param lat Latitude centres in degrees north, one per array latitude row.
#' @param seed Seed for the pseudo-random initial labels; Article 1 fixes seed 1.
#' @param verbose Logical; print native solver and progress information.
#' @param rebuild Logical; force recompilation of the native callback functions.
#' @return List with label_attribution, a one-based longitude x latitude matrix,
#'   and 'Data and smooth cost', the native unary/seam energies. Runtime is printed
#'   when verbose but is not included in the returned list.
#' @details R labels become zero-based for C++ and one-based on return. Unary
#'   costs use cos(latitude); an edge uses mean endpoint latitude weight times
#'   the sum of endpoint Hellinger differences between its two selected models.
#'   Longitude wraps east-west; latitude has no polar wrap. The pinned custom
#'   solver revision is enforced by production entry points, not this low-level API.
#'   Native callbacks deliberately accumulate in float precision. This function
#'   resets the process RNG kind/state and uses a locally converged swap solution;
#'   it does not promise a global optimum for arbitrary multi-label problems.
#' @examplesIf requireNamespace("gcoWrapR", quietly = TRUE) && requireNamespace("RcppXPtrUtils", quietly = TRUE)
#' p <- array(0.5, c(3, 2, 2, 2))
#' h <- array(0.1, c(3, 2, 2))
#' fit <- GraphCutHellinger_nD_lat(p, h, 1, 0.125, 2, c(-30, 30),
#'                                seed = 1, verbose = FALSE, rebuild = FALSE)
#' dim(fit$label_attribution)
GraphCutHellinger_nD_lat <- function(
  pdf_models_future,
  h_dist,
  weight_data,
  weight_smooth,
  nBins,
  lat,
  seed,
  verbose,
  rebuild
) {
  dims <- dim(pdf_models_future)
  if (length(dims) != 4L || dims[1] < 3L || dims[2] < 2L || dims[4] < 2L ||
      !identical(as.integer(dim(h_dist)),as.integer(dims[c(1,2,4)])) ||
      !is.numeric(lat) || length(lat) != dims[2] || any(!is.finite(lat)) || any(abs(lat)>90) ||
      length(nBins)!=1L || is.na(nBins) || nBins != dims[3] ||
      any(!is.finite(h_dist)) || any(h_dist < 0) ||
      length(weight_data)!=1L || !is.finite(weight_data) || weight_data<0 ||
      length(weight_smooth)!=1L || !is.finite(weight_smooth) || weight_smooth<0 ||
      length(seed)!=1L || !is.finite(seed)) {
    stop("Graph cut requires [lon>=3,lat>=2,bins,models>=2] PMFs, matching unary distances, valid latitudes, seed and non-negative weights.")
  }
  validate_target_pmfs(pdf_models_future)
  require_packages(c("gcoWrapR", "RcppXPtrUtils"))
  GCoptimizationGridGraph <- gcoWrapR::GCoptimizationGridGraph
  cppXPtr <- RcppXPtrUtils::cppXPtr
  n_labs <- dim(h_dist)[3]
  height <- ncol(pdf_models_future[,,,1])
  width <- nrow(pdf_models_future[,,,1])

  if (isTRUE(verbose)) {
    cat("GraphCut grid:", width, "x", height, ", labels:", n_labs, "\n")
  }

  if (length(lat) != height) {
    stop("Error: Latitude vector length does not match grid height.")
  }

  lat_weights <- matrix(
    rep(abs(cos(lat * pi / 180)), each = width),
    nrow = width,
    ncol = height,
    byrow = FALSE
  )

  if (isTRUE(verbose)) {
    cat("Latitude weights dim:", paste(dim(lat_weights), collapse = " x "), "\n")
  }

  # Joint-bin index varies fastest in the seam callback; unary arrays instead
  # use pixel-major cells followed by label. Changing either layout changes costs.
  h_dist_cpp <- c(aperm(h_dist, c(1, 2, 3)))
  pdf_models_cpp <- c(aperm(pdf_models_future, c(3, 1, 2, 4)))
  lat_weights_cpp <- c(aperm(lat_weights, c(1, 2)))

  # Instantiate the GraphCut environment
  gco <- new(GCoptimizationGridGraph, width, height, n_labs)

  if (isTRUE(verbose)) {
    print(gco)
  }

  # DataCost function using 2D latitude weights
  ptrDataCost <- cppXPtr(
    code = 'float dataFn(int p, int l, Rcpp::List extraData)
{
    int width         = extraData["width"];
    int height        = extraData["height"];
    int numPix        = width * height;
    float weight_global = extraData["weight"];
    Rcpp::NumericVector data        = extraData["data"];
    Rcpp::NumericVector lat_weights = extraData["lat_weights"];

    float lat_weight = lat_weights[p];
    float dval       = data[p + numPix * l];

    return(weight_global * dval * lat_weight);
}',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = rebuild, showOutput = FALSE, verbose = FALSE
  )

  if (isTRUE(verbose)) cat("Creating SmoothCost function...  ")
  ptrSmoothCost <- cppXPtr(
    code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
{
    int width = extraData["width"];
    int height = extraData["height"];
    int numPix = width * height;
    float weight_global = extraData["weight"];
    NumericVector data = extraData["data"];
    NumericVector lat_weights = extraData["lat_weights"];
    int nBins = extraData["nBins"];

    // Ensure p1 and p2 are within the valid range
    if (p1 < 0 || p1 >= numPix || p2 < 0 || p2 >= numPix) {
        Rcpp::stop("Pixel index out of bounds in smoothFn.");
    }

    // Retrieve latitude weights safely
    float lat_weight = (lat_weights[p1] + lat_weights[p2]) / 2.0;

    float cost = 0.0f;
    float tmp1 = 0.0f;
    float tmp2 = 0.0f;
    float diff1, diff2;
    int index1, index2;

    int offset_p1_l1 = (p1 + numPix * l1) * nBins;
    int offset_p1_l2 = (p1 + numPix * l2) * nBins;
    int offset_p2_l1 = (p2 + numPix * l1) * nBins;
    int offset_p2_l2 = (p2 + numPix * l2) * nBins;

    for (int i = 0; i < nBins; i++) {
        index1 = offset_p1_l1 + i;
        index2 = offset_p1_l2 + i;
        diff1 = sqrt(data[index1]) - sqrt(data[index2]);
        tmp1 += diff1 * diff1;

        index1 = offset_p2_l1 + i;
        index2 = offset_p2_l2 + i;
        diff2 = sqrt(data[index1]) - sqrt(data[index2]);
        tmp2 += diff2 * diff2;
    }

    cost = (sqrt(tmp1) + sqrt(tmp2)) / sqrt(2.0f);
    return weight_global * lat_weight * cost;
}
',
    includes = c("#include <math.h>", "#include <Rcpp.h>"),
    rebuild = rebuild, showOutput = isTRUE(verbose), verbose = FALSE
  )

  # Set DataCost and SmoothCost
  gco$setDataCost(ptrDataCost, list(
    numPix = width * height,
    width = width,
    height = height,
    data = h_dist_cpp,
    weight = weight_data,
    lat_weights = lat_weights_cpp
  ))

  gco$setSmoothCost(ptrSmoothCost, list(
    width = width,
    height = height,
    data = pdf_models_cpp,
    weight = weight_smooth,
    nBins = nBins,
    lat_weights = lat_weights_cpp
  ))

  # Initialize labels pseudo-randomly
  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion") # For reproducibility across OS
  set.seed(seed)

  init_labels <- sample.int(n_labs, width * height, replace = TRUE) - 1
  for (z in 0:((width * height) - 1)) {
    gco$setLabel(z, init_labels[z + 1])
  }

  # Perform graph cut optimization
  if (isTRUE(verbose)) {
    cat("Starting GraphCut optimization...  ")
    print(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  }
  begin <- Sys.time()
  gco$swap(-1)  # Run until convergence
  time_spent <- Sys.time() - begin
  if (isTRUE(verbose)) {
    print(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    cat("GraphCut optimization done :  ")
    print(time_spent)
  }

  # Retrieve results
  data_cost <- gco$giveDataEnergy()
  smooth_cost <- gco$giveSmoothEnergy()
  data_smooth_list <- list("Data cost" = data_cost, "Smooth cost" = smooth_cost)

  # Extract label attribution
  label_attribution <- matrix(0, nrow = width, ncol = height)
  for (j in 1:height) {
    for (i in 1:width) {
      label_attribution[i, j] <- gco$whatLabel((i - 1) + width * (j - 1))
    }
  }
  label_attribution <- label_attribution + 1

  # Clean up and force garbage collection
  rm(gco)
  gc()

  # Return results
  list(
    "label_attribution" = label_attribution,
    "Data and smooth cost" = data_smooth_list
  )
}

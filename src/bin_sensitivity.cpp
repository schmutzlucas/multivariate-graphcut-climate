#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <climits>
using namespace Rcpp;

inline int equal_width_bin_cpp(const double x, const double xmin, const double xmax, const int k) {
  if (NumericVector::is_na(x) || std::isnan(x)) return -1;
  if (x == R_NegInf) return 0;
  if (x == R_PosInf) return k - 1;
  double scaled = static_cast<double>(k) * (x - xmin) / (xmax - xmin);
  const double nearest = std::round(scaled);
  const double tolerance = 8.0 * DBL_EPSILON * std::fmax(1.0, std::fabs(scaled));
  if (std::isfinite(scaled) && std::fabs(scaled - nearest) <= tolerance) scaled = nearest;
  int bin = static_cast<int>(std::floor(scaled));
  if (bin < 0) bin = 0;
  if (bin >= k) bin = k - 1;
  return bin;
}

// [[Rcpp::export]]
List multik_hist_counts_cpp(
    NumericVector pr,
    NumericVector tas,
    NumericVector psl,
    NumericVector range_block,
    LogicalVector split_a,
    IntegerVector k_values) {
  IntegerVector raw_dim = pr.attr("dim");
  IntegerVector tas_dim = tas.attr("dim");
  IntegerVector psl_dim = psl.attr("dim");
  IntegerVector range_dim = range_block.attr("dim");
  if (raw_dim.size() != 3 || range_dim.size() != 4 || tas_dim.size() != 3 || psl_dim.size() != 3)
    stop("Invalid array dimensions.");
  for (int d = 0; d < 3; ++d) {
    if (raw_dim[d] <= 0 || tas_dim[d] != raw_dim[d] || psl_dim[d] != raw_dim[d])
      stop("Variable array shapes differ or contain empty dimensions.");
  }
  const int nlon = raw_dim[0];
  const int nlat = raw_dim[1];
  const int ntime = raw_dim[2];
  const double cells_d = static_cast<double>(nlon) * nlat;
  if (cells_d > INT_MAX / 6) stop("Spatial grid exceeds integer support indexing capacity.");
  const int ncell = static_cast<int>(cells_d);
  if (tas.size() != pr.size() || psl.size() != pr.size() || split_a.size() != ntime)
    stop("Variable arrays and split mask have inconsistent dimensions.");
  if (range_dim[0] != nlon || range_dim[1] != nlat || range_dim[2] != 3 || range_dim[3] != 2)
    stop("range_block must have dimensions [lon, lat, 3, 2].");

  const int nk = k_values.size();
  if (nk == 0) stop("At least one histogram resolution is required.");
  for (int t = 0; t < ntime; ++t)
    if (split_a[t] == NA_LOGICAL) stop("Split mask must not contain NA.");
  for (R_xlen_t r = 0; r < range_block.size(); ++r)
    if (!std::isfinite(range_block[r])) stop("Histogram support must be finite.");
  for (int cell = 0; cell < ncell; ++cell)
    for (int v = 0; v < 3; ++v)
      if (range_block[cell + ncell * v + ncell * 3] <= range_block[cell + ncell * v])
        stop("Every histogram maximum must be greater than its minimum.");
  IntegerVector offsets(nk);
  double total_length_d = 0.0;
  for (int ki = 0; ki < nk; ++ki) {
    if (k_values[ki] == NA_INTEGER || k_values[ki] < 2)
      stop("Histogram resolutions must be integers at least two.");
    const double states = std::pow(static_cast<double>(k_values[ki]), 3.0);
    const double next_length = total_length_d + states * static_cast<double>(ncell) * 2.0;
    // Offsets/state arithmetic below is intentionally integer-based. Reject an
    // unsupported long layout BEFORE narrowing, rather than silently overflowing.
    if (!std::isfinite(next_length) || next_length > INT_MAX)
      stop("Count layout exceeds integer offset capacity; reduce the spatial block or resolutions.");
    offsets[ki] = static_cast<int>(total_length_d);
    total_length_d = next_length;
  }
  if (total_length_d > static_cast<double>(R_XLEN_T_MAX)) stop("Count array is too large.");
  IntegerVector counts(static_cast<R_xlen_t>(total_length_d));

  for (int cell = 0; cell < ncell; ++cell) {
    const double pr_min = range_block[cell + ncell * 0 + ncell * 3 * 0];
    const double tas_min = range_block[cell + ncell * 1 + ncell * 3 * 0];
    const double psl_min = range_block[cell + ncell * 2 + ncell * 3 * 0];
    const double pr_max = range_block[cell + ncell * 0 + ncell * 3 * 1];
    const double tas_max = range_block[cell + ncell * 1 + ncell * 3 * 1];
    const double psl_max = range_block[cell + ncell * 2 + ncell * 3 * 1];

    for (int t = 0; t < ntime; ++t) {
      const R_xlen_t raw_index = cell + static_cast<R_xlen_t>(ncell) * t;
      const double x0 = pr[raw_index];
      const double x1 = tas[raw_index];
      const double x2 = psl[raw_index];
      if (NumericVector::is_na(x0) || std::isnan(x0) ||
          NumericVector::is_na(x1) || std::isnan(x1) ||
          NumericVector::is_na(x2) || std::isnan(x2)) continue;
      const int split = split_a[t] == TRUE ? 0 : 1;

      for (int ki = 0; ki < nk; ++ki) {
        const int k = k_values[ki];
        const int b0 = equal_width_bin_cpp(x0, pr_min, pr_max, k);
        const int b1 = equal_width_bin_cpp(x1, tas_min, tas_max, k);
        const int b2 = equal_width_bin_cpp(x2, psl_min, psl_max, k);
        if (b0 < 0 || b1 < 0 || b2 < 0) continue;
        const int state = b0 + k * b1 + k * k * b2;
        const int nstates = k * k * k;
        const R_xlen_t out_index = offsets[ki] + state +
          static_cast<R_xlen_t>(nstates) * (cell + ncell * split);
        counts[out_index] += 1;
      }
    }
  }
  return List::create(
    _["counts"] = counts,
    _["offsets"] = offsets,
    _["k_values"] = k_values,
    _["ncell"] = ncell
  );
}

// Count lists are also callable from R. Validate their complete memory layout
// before using unchecked native indexing in the Hellinger loop.
inline void validate_count_layout(const List& input) {
  const char* required[] = {"counts", "offsets", "k_values", "ncell"};
  for (const char* name : required)
    if (!input.containsElementNamed(name)) stop("Count layout is missing a required field.");
  SEXP raw_counts = input["counts"], raw_offsets = input["offsets"], raw_ks = input["k_values"];
  if (TYPEOF(raw_counts) != INTSXP || TYPEOF(raw_offsets) != INTSXP || TYPEOF(raw_ks) != INTSXP)
    stop("Count vectors, offsets and resolutions must be stored as integers.");
  IntegerVector counts = input["counts"];
  IntegerVector offsets = input["offsets"];
  IntegerVector ks = input["k_values"];
  const double cells_d = as<double>(input["ncell"]);
  if (!std::isfinite(cells_d) || cells_d < 1 || cells_d != std::floor(cells_d) || cells_d > INT_MAX)
    stop("Count layout needs a positive integer cell count.");
  if (ks.size() == 0 || offsets.size() != ks.size()) stop("Invalid count-layout offsets/resolutions.");
  double length = 0;
  for (int i = 0; i < ks.size(); ++i) {
    if (ks[i] == NA_INTEGER || ks[i] < 2 || offsets[i] == NA_INTEGER || offsets[i] != length)
      stop("Count-layout resolution or offset is invalid.");
    length += std::pow(static_cast<double>(ks[i]), 3.0) * cells_d * 2;
    if (!std::isfinite(length) || length > INT_MAX) stop("Count layout exceeds integer offset capacity.");
  }
  if (length != counts.size()) stop("Count vector length does not match its layout.");
  for (R_xlen_t i = 0; i < counts.size(); ++i)
    if (counts[i] == NA_INTEGER || counts[i] < 0) stop("Counts must be nonnegative integers without NA.");
}

// [[Rcpp::export]]
NumericVector multik_hellinger_cpp(List model, List reference) {
  validate_count_layout(model);
  validate_count_layout(reference);
  IntegerVector model_counts = model["counts"];
  IntegerVector ref_counts = reference["counts"];
  IntegerVector offsets = model["offsets"];
  IntegerVector ref_offsets = reference["offsets"];
  IntegerVector k_values = model["k_values"];
  IntegerVector ref_k_values = reference["k_values"];
  const int ncell = as<int>(model["ncell"]);
  const int ref_ncell = as<int>(reference["ncell"]);
  if (ncell != ref_ncell || k_values.size() != ref_k_values.size())
    stop("Model and reference count layouts differ.");
  for (int i = 0; i < k_values.size(); ++i)
    if (k_values[i] != ref_k_values[i]) stop("Model and reference count layouts differ.");
  const int nk = k_values.size();
  NumericVector output(static_cast<R_xlen_t>(ncell) * nk * 2, NA_REAL);
  output.attr("dim") = IntegerVector::create(ncell, nk, 2);

  for (int split = 0; split < 2; ++split) {
    for (int ki = 0; ki < nk; ++ki) {
      const int k = k_values[ki];
      const int nstates = k * k * k;
      for (int cell = 0; cell < ncell; ++cell) {
        const R_xlen_t model_start = offsets[ki] +
          static_cast<R_xlen_t>(nstates) * (cell + ncell * split);
        const R_xlen_t ref_start = ref_offsets[ki] +
          static_cast<R_xlen_t>(nstates) * (cell + ncell * split);
        double model_total = 0.0;
        double ref_total = 0.0;
        for (int state = 0; state < nstates; ++state) {
          model_total += model_counts[model_start + state];
          ref_total += ref_counts[ref_start + state];
        }
        if (model_total <= 0.0 || ref_total <= 0.0) continue;
        double squared = 0.0;
        for (int state = 0; state < nstates; ++state) {
          const double p = model_counts[model_start + state] / model_total;
          const double q = ref_counts[ref_start + state] / ref_total;
          const double delta = std::sqrt(p) - std::sqrt(q);
          squared += delta * delta;
        }
        output[cell + ncell * ki + static_cast<R_xlen_t>(ncell) * nk * split] =
          std::sqrt(squared) / std::sqrt(2.0);
      }
    }
  }
  return output;
}

# Multivariate graph-cut combination of climate models

Code, compact numerical results and documentation for the article

> **Beyond mean fields: a multivariate graph-cut approach to climate-model combination**
> Lucas Schmutz and co-authors. Manuscript prepared for *Earth System Dynamics* (ESD).
> The article DOI will be added here when it is available.

The repository regenerates every figure and table of the article from included
numerical results. It demonstrates the graph-cut solver on a small example and
provides entry points for re-running the full ERA5 and perfect-model experiments
with external climate data.

## Scientific background

Multi-model means average climate-model fields and can smooth away
distributional features. The article instead combines models **locally, using
their full multivariate distributions**. At every grid point, one candidate
CMIP6 model's empirical joint distribution of precipitation (`pr`), near-surface
temperature (`tas`) and sea-level pressure (`psl`) is selected. A graph-cut
energy is minimised over the global grid:

- **data term:** Hellinger distance between each candidate's calibration-period
  (1951–1975) joint PMF and the reference, weighted by cos(latitude);
- **smoothness term:** λ times the sum of inter-model Hellinger distances
  between the two neighbouring labels, computed from target-period candidate PMFs.
  The grid is four-connected and periodic in longitude.

The minimisation uses alpha-beta swap with a longitude-periodic build of the GCO
library (main configuration λ = 0.125, seed 1). The result is a spatial
patchwork of model distributions. It is evaluated against D-MMM, the
equal-weight mixture of candidate PMFs, using the Hellinger distance and a
low-density-region (LDR) Hellinger score. Mean-field biases are compared with
the conventional multi-model mean (MMM). See [docs/METHODOLOGY.md](docs/METHODOLOGY.md).

| Experiment | Reference | Candidates | Calibration | Evaluation |
| --- | --- | --- | --- | --- |
| ERA5 | ERA5 reanalysis | 22 CMIP6 models (`config/era5_models.txt`) | 1951–1975 | 1999–2023 |
| Perfect model (PME) | each of 15 models in turn | the other 14 (`config/pme_models.txt`) | 1951–1975 | 1976–2000, 2001–2025, 2026–2050, 2051–2075, 2076–2100 |

The PME uses a 15-model subset of the 22 candidates to limit near-duplicate
models. CMCC-ESM2, CNRM-ESM2-1, HadGEM3-GC31-LL, HadGEM3-GC31-MM, INM-CM4-8,
KACE-1-0-G and UKESM1-0-LL are excluded. The article describes this screening
as pragmatic and partly empirical, not a formal independence test. Trivariate PMFs use 8 bins per
variable on fixed per-cell supports (`data/supports/`). Precipitation is clipped
at zero and transformed with `log1p`. Scientific defaults live in `config/article1.R`.

## Repository contents

```text
R/            distributions, metrics, calendars/NetCDF I/O, graph-cut energy, plotting
src/          C++ kernel for the PMF-resolution sensitivity
config/       scientific settings, model lists, grid, tested software versions
experiments/  full ERA5/PME runs, evaluation, sensitivities and the compact demo
figures/      one renderer per article figure family
tables/       generator for all article tables
results/      compact numerical results used by figures/tables (corrected cohort)
data/         exact histogram supports, input provenance, public data manifest
publication/  figure/table assets exactly as published, manifests, checksums
tests/        mathematical, regression and solver tests; small real-data fixture
tools/        installation, reproduction, data preparation and manifest tools
docs/         methodology, data, provenance, known issues, offline R reference
```

`publication/article_inventory.csv` maps every manuscript figure, table and
in-text result to its LaTeX label, generation script, numerical inputs,
experiment, configuration, verification status and limitations.

## Installation

R ≥ 4.4 with a C++ toolchain (Rtools on Windows); tested with R 4.4.3. Python ≥ 3.10
is used only for checksum and manifest tests. Linux needs the NetCDF, curl,
OpenSSL, libpng and font development libraries listed in `.github/workflows/`.

```text
Rscript --vanilla tools/install_dependencies.R --accept-current-cran
```

This installs current CRAN packages. `ncdf4.helpers` and `PCICt`, now archived
on CRAN, are installed from the CRAN Archive at the tested versions. To restore
the exact tested versions instead, use `renv.lock`; see [DEPENDENCIES.md](DEPENDENCIES.md).

The graph-cut solver is only needed for Levels 2 and 3. It is a separate package,
not included here:

```r
remotes::install_github("schmutzlucas/gcoWrapR@6b655d5feff7e7912d2ac0a7de95d0f8edd6d434")
```

Its bundled GCO library has its own research-use and redistribution terms; read
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) before installing it. The code
checks for this exact commit. A standard rectangular-grid GCO build changes the
graph and is not an equivalent substitute.

## Level 1 — reproduce all figures and tables (no climate data)

```text
Rscript --vanilla tools/reproduce.R
```

Takes about one minute. Every figure family and table is rendered from
`results/` into `outputs/` (set `MVGC_OUTPUT_ROOT` to change this). Input
checksums are verified first, and the archived article assets in `publication/`
are never overwritten. Single families can be run directly, e.g.
`Rscript --vanilla figures/pme_overlays.R` or `Rscript --vanilla tables/render.R`.

Verified against the manuscript (Overleaf revision `02e1166`; Windows, R 4.4.3):

- **All five tables and the quoted in-text values** (PME pooled gains, A3 values)
  are identical to the manuscript.
- **All 17 figure pairs** regenerate pixel-identical to the figures in the
  manuscript; the manuscript's figure files are the repository's renders.

Other systems and graphics library versions may change PDF bytes or anti-aliasing
without changing the plotted numbers.

## Level 2 — graph-cut demonstration and tests

```text
Rscript --vanilla experiments/demo.R
Rscript --vanilla tests/run_tests.R
python -m unittest discover -s tests -p "test_*.py" -v
```

The demo runs the pinned solver on synthetic PMFs on a 4×2 periodic grid. It
reports labels and the data/smoothness energy. The tests check the Hellinger/LDR
definitions, PMF construction, calendars and precipitation QC; the native energy
callbacks against an independent R energy, the exhaustive binary optimum, periodic
edges and seed repeatability; distances on a real production PMF patch; and
the checksums of all compact inputs and supports. Tests that need the solver
report `SKIP` when it is absent.

## Level 3 — full experiments (external data required)

The full experiments need the prepared daily 1° climate collection: 23 datasets
× 3 variables, about 803 GB. It is **not** included. The run takes many hours
and tens of GiB of memory/disk.

```text
MVGC_DATA_ROOT=/path/to/prepared
MVGC_ERA5_SUPPORT=data/supports/era5_support.rds
MVGC_PME_SUPPORT=data/supports/pme_support.rds
MVGC_OUTPUT_ROOT=/path/to/new/output
Rscript --vanilla experiments/era5/run.R
Rscript --vanilla experiments/perfect_model/run.R      # PME_LAMBDAS=0.125 for the main setting only
Rscript --vanilla experiments/sensitivity/histogram_resolution.R
Rscript --vanilla experiments/sensitivity/spatial_regularisation.R
Rscript --vanilla experiments/sensitivity/ldr_threshold.R
```

Read [docs/FULL_EXPERIMENTS.md](docs/FULL_EXPERIMENTS.md),
[docs/DATA_REQUIREMENTS.md](docs/DATA_REQUIREMENTS.md) (directory layout,
grid, calendars, units) and [docs/DATA_PREPROCESSING.md](docs/DATA_PREPROCESSING.md)
(download, UTC aggregation, regridding, QC) first. **The full experiments have
not been independently re-run from raw data for this release.** Results built
from newly downloaded data may differ in detail from the published values.

## Data availability and provenance

[data/DATA_MANIFEST.csv](data/DATA_MANIFEST.csv) lists every data asset in four
categories, explained in [docs/DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md):

1. **included** here, with SHA-256: compact results, the two exact histogram
   supports, the article figure/table assets and a test fixture;
2. **external public** source collections: CMIP6 projections on the Copernicus
   Climate Data Store (https://doi.org/10.24381/cds.c866074c) and ERA5 hourly
   single levels (https://doi.org/10.24381/cds.adbb2d47). These are collection
   DOIs, not DOIs of the prepared files;
3. **required but not archived:** the 69 prepared daily files, with recorded
   bytes and SHA-256, and the large intermediate workspaces;
4. **pending verification:** redistribution terms for prepared files and
   unresolved upstream version qualifications.

Exact model, member and version identities are in `data/metadata/source_manifest.csv`.

## Limitations and known issues

See [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md). In brief: the histogram
supports pool evaluation-period data. Some upstream source versions remain qualified. Full
experiments require non-archived prepared data. LDR scores are distributional
diagnostics, not validation of physical extremes.

## Citation

Please cite the article and this software. Citation metadata are in
[CITATION.cff](CITATION.cff); the archived release DOI will be added after
deposit. The graph-cut method builds on Thao et al. (2022). The GCO notice
requests citation of Boykov, Veksler and Zabih (2001), Kolmogorov and Zabih
(2004) and Boykov and Kolmogorov (2004).

## Licence

- **Code** (`R/`, `src/`, `experiments/`, `figures/*.R`, `tables/`, `tests/`,
  `tools/`, `config/`): MIT, see [LICENSE](LICENSE).
- **Results, data records, figures and documentation** (`results/`, `data/`,
  `publication/`, `docs/`, Markdown files): CC BY 4.0, see [LICENSE-DATA](LICENSE-DATA).

Neither licence covers the external graph-cut solver, R/CRAN dependencies, or the
ERA5/CMIP6 source data, which keep their own terms. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and `config/dependency_licences.csv`.

Contact: Lucas Schmutz, University of Lausanne (lucas.schmutz@unil.ch).

# Reproducing the article

Three levels are distinguished. Only Levels 1 and 2 have been executed for this
release; Level 3 is documented but has not been independently re-run.

| Level | Command | Needs | Verified |
| --- | --- | --- | --- |
| 1 Figures and tables | `Rscript --vanilla tools/reproduce.R` | this repository + CRAN packages | yes (local Windows/R 4.4.3; Linux CI) |
| 2 Algorithm | `Rscript --vanilla experiments/demo.R`, `tests/run_tests.R` | + pinned solver | yes (local; Linux CI job installs solver) |
| 3 Full experiments | `experiments/*/run.R` | + ~803 GB prepared data | not re-run for this release |

## Level 1

`tools/reproduce.R` runs eight independent subprocesses (`figures/*.R`,
`tables/render.R`). Before rendering, it checks the SHA-256 of every compact
input listed in `results/*/manifest.csv`, and fails if any is missing or changed.
Output goes to `MVGC_OUTPUT_ROOT` (default `outputs/`): 17 PDF/PNG figure pairs
and five LaTeX tables. The archived article assets in `publication/` are never
written. Select one family with, e.g., `Rscript --vanilla tools/reproduce.R pme_overlays`.

Explicit input overrides through `MVGC_*`/`PME_*` variables are recorded, but
their outputs are not article data. Clear these variables for ordinary reproduction.

Inputs are compact: map arrays (tas K, pr mm/day, psl Pa, with coordinates),
unrounded scores and the exact density-curve coordinates of the published
figures. Level 1 lays out these numbers; it does not re-estimate densities or
rerun the optimiser. PME tables derive from all 450 reference/period/diagnostic
score pairs. The histogram figure derives from saved cellwise correlations.

### Agreement with the final manuscript (Overleaf `f22ca16`)

| Check | Result |
| --- | --- |
| archived figure assets vs manuscript files (34 files) | byte-identical |
| regenerated figures vs published (pixel comparison of PNG and 150-dpi PDF rasters) | 14/17 identical; Fig. 5 differs in 54 anti-aliased pixels (≤ 2/765), same data; Figs. 1–2 differ (older cohort, see below) |
| generated table values vs manuscript tables (Tables 1–3, B1–B2) | identical row by row |
| in-text values (PME pooled gains 25.9/25.7/24.8/20.5/13.5 %; A3 Hellinger/fragmentation values; D-MMM 0.24792) | identical |

PDF bytes always differ between renders because of embedded metadata. Pixel and
numerical comparisons are the meaningful checks. Other operating systems or
graphics library versions can change anti-aliasing and fonts.

### ERA5 maps (Figs. 1–2)

The published ERA5 maps were rendered from the research caches that preceded the
October 2026 data correction. The tables and text use the corrected cohort. At the
author's request, the archive keeps the published figure files unchanged. Running
the renderers on the included corrected arrays yields maps with the same layout
and slightly different values. Details: [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md),
[docs/ERA5_RECONCILIATION.md](docs/ERA5_RECONCILIATION.md).

## Level 2

`experiments/demo.R` solves a synthetic 4×2 periodic problem with the pinned
solver (`schmutzlucas/gcoWrapR@6b655d5`) and prints labels and energies. The
installed solver commit is checked at runtime. `tests/run_tests.R` runs 22 R test
scripts in isolated processes. Native tests compare callbacks with an independent
R energy and check the exhaustive binary optimum, periodic edges and seed
repeatability. Without the solver they report `SKIP`, never `PASS`.

## Level 3

See [docs/FULL_EXPERIMENTS.md](docs/FULL_EXPERIMENTS.md). Exact-input reproduction
requires the prepared files listed in `data/DATA_MANIFEST.csv` (not publicly
archived), the included supports and the pinned solver. Replication from newly
downloaded ERA5/CMIP6 data is possible with
[docs/DATA_PREPROCESSING.md](docs/DATA_PREPROCESSING.md). Different source
versions, preprocessing, libraries or hardware may change individual labels and
percentages. Quantitative equivalence is not guaranteed and has not been tested.
The article's sensitivity analyses cover only PMF resolution, LDR threshold and λ.

## Preserved originals

Generated outputs never go into `publication/`, `results/` or `data/`. Test and
runtime logs go to a temporary directory, or to `MVGC_PRIVATE_VERIFICATION_DIR`.

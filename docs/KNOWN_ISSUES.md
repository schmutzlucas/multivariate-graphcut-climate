# Known issues and limitations

This file lists every known limitation that affects how the article's results
can be reproduced or interpreted. It separates deficiencies of the **original**
(pre-October 2026) data collection from the status of the **corrected cohort**
whose numbers are reported in the final manuscript and included here.
Nothing below changes a published number.

## 1. ERA5 maps (Figs. 1–2): resolved

All article figures, tables and quoted values now come from the same corrected
cohort. Earlier manuscript versions (up to Overleaf `f22ca16`) showed ERA5 maps
rendered from pre-correction research caches. Those maps were replaced in the
manuscript by renders of the included corrected arrays (Overleaf `a19496f`);
the archived assets are byte-identical to them. The superseded figures are not
in this repository. Their identities and the cell-wise differences (area means
change only in the fourth decimal; GC labels changed in a few per cent of the
globe) are kept in `data/metadata/era5_asset_provenance.json`,
`era5_field_comparison.csv` and `era5_bias_comparison.csv`. All text and table
values were checked and did not need to change.

## 2. Histogram supports pool evaluation-period information

Supports (per-cell bin ranges) are data-informed and shared within each experiment:

- **ERA5:** one support pools ERA5 and the 22 candidates over 1951–1975 **and
  1999–2023**. ERA5 target-period extrema therefore influence bin edges. They
  do not enter the unary or seam energies.
- **PME:** one support pools all 15 models over **1951–2100**. Each held-out
  reference's future data, and all future periods, therefore influence bin edges
  for every reference/period case. The held-out model is excluded from candidate
  labels and from D-MMM.

The results are conditional on this discretisation. Performance with a
calibration-only support has not been demonstrated. Both exact supports are
included (`data/supports/`). The manuscript Methods describe this pooling.

## 3. Original data deficiencies (historical, corrected for the reported cohort)

These affected the **original** collection and the older research caches. The
corrected cohort used re-acquired or re-prepared inputs. Each prepared input
is recorded in `data/metadata/source_manifest.csv` and `data/input_provenance.json`.

- **Missing day-31 samples.** The original CDS CMIP6 daily requests selected days
  1–30 only. As a result, 18 models lacked 175 daily samples per audited 25-year
  window. The corrected collection is calendar-complete: each support records
  per-year day counts, and the 360-day models correctly have no day 31.
- **Precipitation units/metadata.** Prepared files are interpreted as effective
  mm/day. Four native 360-day products (HadGEM3-GC31-LL/MM, KACE-1-0-G,
  UKESM1-0-LL) contain mm/day values despite stale flux attributes; they must not
  be re-converted. ERA5 precipitation was rebuilt as UTC-day accumulations in
  mm/day. Finite negative values are clipped to zero before `log1p`. The supports
  record QC statistics per file.
- **Support provenance.** The original archived supports had incomplete creation
  metadata. The old ERA5 support builder selected the PME model list. The corrected
  supports record exact model/year scopes, contributing file hashes and the
  input-registry hash.

## 4. Source identity qualifications (unresolved)

- **ACCESS-CM2 historical:** published v20191108 was used. A publisher registration
  advertises v20191118, but no retrievable release was found.
- **BCC-CSM2-MR historical:** r1i1p1f1 v20181126 has an unresolved
  replacement/unpublication history. The final pairing is member-coherent r1i1p1f1.
- **ERA5:** original CDS request receipts are not retained for every hourly month.
  The tas/psl daily aggregation receipts are incomplete.
- Exact member/version identities are recorded where known. Unknown values are
  left empty rather than inferred.

## 5. Reproduction scope

- **Level 1** (included data): all 17 figure pairs and five tables render. 16 pairs
  are pixel-identical to the article. Fig. 5 differs only by anti-aliasing in 54
  pixels, with identical plotted data. Table numbers equal the manuscript row by row.
- **Level 2** (synthetic demo and tests with the pinned solver): verified locally
  on Windows/R 4.4.3. A Linux CI job installs the solver from its public repository.
- **Level 3** (full experiments): requires the ~803 GB prepared collection, which
  is not publicly archived. No full ERA5/PME rerun was performed for this release,
  and no independent end-to-end reproduction has been demonstrated.
- The native solver uses float arithmetic in callbacks and finds a local minimum
  (alpha-beta swap, seed 1). Bit-identical labels across platforms are not promised.

## 6. Diagnostic interpretation

- LDR bins are chosen from the reference PMF in ascending probability until the
  nominal mass is reached. Zero bins and all boundary ties are included, and
  restricted PMFs are not renormalised. Attained mass can exceed the threshold.
  LDR scores are distributional diagnostics, not validation of physical extremes.
- Pooled area-mean gains differ from the mean of the 15 reference-wise gains;
  both are reported where relevant.
- The graph-cut product is a spatial patchwork of distributions, not a coherent
  daily space-time simulation.

## 7. Licensing and third-party software

- Code: MIT (`LICENSE`). Results, data records, figures and documentation:
  CC BY 4.0 (`LICENSE-DATA`). Scope in `THIRD_PARTY_NOTICES.md`.
- The longitude-periodic GCO solver (`schmutzlucas/gcoWrapR@6b655d5`) is **not**
  included. Its bundled GCO v3 notice imposes research-only use and redistribution
  conditions, while its wrapper metadata declares GPL (>=2). This repository only
  references the solver; users obtain it from its own repository under its terms.
- ERA5/CMIP6 source data remain under their providers' terms.

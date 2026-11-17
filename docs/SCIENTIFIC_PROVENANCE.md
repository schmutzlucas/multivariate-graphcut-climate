# Scientific source and result provenance

Scientific modules and entry points originate in
`schmutzlucas/MV_graphcut_paper`, selected revision
`69a6faf8b20451a013425a7d72290e9b30a8221f`, with final plotting code from
`paper1-submission` `0671145c239d6b8ee9e53760bef46ac4d569614b`. This public
repository has its own history and omits internal development records. The
manuscript verified against is `schmutzlucas/MV_graphcut_overleaf`
`a19496fb5ce37f967c974a30b2eb9053ce97be62`, including the separate supplement.
Both research repositories are private and are cited for provenance only.

The accepted original computations used different frozen source records:

- PME source record: `3abff1f5f8d0a965f18e61ef3834c953cdca9698`, 15 models,
  75 reference-period cases, nonnegative precipitation and one 1951–2100 support.
- Histogram sensitivity source record:
  `e5d08073971484622622cacc6a98f2a2ee83ec88`, 22 candidates, calibration
  1951–1975, k=4:12.
- ERA5 corrected workspace identity:
  `ed0fb12586b7fd06b1d60e4668d85736394a313ea5b25f2c62d0e78546fcec1c`.
  Its compact source records do not establish one unambiguous producing Git
  commit. Exact original source lineage remains qualified rather than inferred
  from a nearby legacy run.
- Spatial solutions are tied to the corrected ERA5 workspace, seed 1 and six
  shown lambdas; the .125 solution matches the main selected labels.

The registry/support identities in DATA_REQUIREMENTS.md select corrected inputs.
Input records preserve calendars, effective units, exact model/member/version
identities where known and original prepared checksums. BCC/ACCESS historical
version qualifications and incomplete ERA5 request receipts remain unresolved.
Producer checksums are historical identities, not a fresh rehash of all climate
payloads. The four native 360-day precipitation attributes require effective-unit
interpretation; their source-map controls supported values already in mm/day.

Publication RDS metadata removes private storage locations while retaining
scientific arrays, coordinates, diagnostic labels and meaningful source/input
hashes. These sanitized files have **new file hashes** in the compact manifests;
original source/result hashes must not be confused with their publication-copy
hashes. Source identifiers are not download URLs. Reusable date/NetCDF helper
definitions were retained under descriptive module names; historical repair
orchestration and investigation scripts are excluded. Plotters retain accepted
layout/numerical data, with only the connected Fig. 5 published.

The compact workflow renders 17 figure pairs and five tables, all from the
corrected cohort (ERA5 map history in ERA5_RECONCILIATION.md). It does not establish
new raw-climate/PMF/global-optimizer execution. Original production source
identity and source-to-input acquisition remain qualified where evidence is
incomplete. No claim of universal bit-identical solver labels across platforms
or exact original climate-file reconstruction is made.

## Final manuscript plotting selection

Current manuscript: `schmutzlucas/MV_graphcut_overleaf@d0e224f60e3bbb0de442f5a83b20ee7cb496c709`.
Updated ERA5/A1/A2 plotting source: `schmutzlucas/MV_graphcut_paper@0671145c239d6b8ee9e53760bef46ac4d569614b`.
A3 plotting and validated six-lambda CSV: research revision `3f8cce7a3da597bfe34b938ca51669404c9330f3`,
`analysis/lambda_tradeoff_review/analyse.R` and `lambda_tradeoff_metrics.csv`.
CSV SHA-256: `49554fdb3dbe773d223f45058415a257ba920309008ecbaff63ee624c064bc66`.
The CSV preserves original checkpoint, workspace, registry and support identities.
It adds derived diagnostics to the existing corrected spatial cohort; it replaces
no earlier numerical input. Its minimum tested Hellinger distance is at 0.025.
These plotting revisions do not replace the original experiment-producing
source identities above. The final software identity is the Git commit cited
with the archive; it is distinct from every source revision.

## Manuscript/cohort consistency

All figures, tables and quoted values of manuscript revision `a19496f` come from
the corrected cohort. Earlier ERA5 maps from pre-correction caches were replaced;
see ERA5_RECONCILIATION.md.

## Histogram supports

`data/supports/` holds the production supports built on 2026-10-24 from the
corrected inputs (original SHA-256 `5530eab3...` ERA5, `ca9394ce...` PME; full
values in `data/supports/manifest.csv`). Only the private storage path of each
contributing file was reduced to its basename. `ranges` and all other metadata
are identical (checked with `identical()` in R), and every per-year daily count
matches its declared calendar.

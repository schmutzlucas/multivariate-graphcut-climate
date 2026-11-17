# Prepared-data contract

Expected external hierarchy:

```text
prepared/
  ERA5/pr/pr_ERA5_<dates>.nc
  ERA5/tas/tas_ERA5_<dates>.nc
  ERA5/psl/psl_ERA5_<dates>.nc
  ACCESS-CM2/pr/pr_ACCESS-CM2_<dates>.nc
  ... all configured models and variables ...
```

Use the source manifest for exact original filenames; abbreviated names above
illustrate the directory contract, not exact-file selection. Each model-variable
file must contain the requested complete daily periods and common 360-by-181
grid, compatible native calendar and effective units. Discovery rejects ambiguous
files. Required-year warn policy alone does not certify complete dates.

Set `MVGC_DATA_ROOT` to the external prepared hierarchy and `MVGC_ERA5_SUPPORT` /
`MVGC_PME_SUPPORT` to their appropriate support RDS. The exact production
supports used by the article are included:

```text
MVGC_ERA5_SUPPORT=data/supports/era5_support.rds   # ERA5 + 22 candidates, 1951-1975 and 1999-2023
MVGC_PME_SUPPORT=data/supports/pme_support.rds     # 15 PME models, 1951-2100
```

They are not selected silently; set the variables explicitly. Bounds are arrays
[360,181,2] for each variable, indexed [0:359 longitude,-90:90 latitude,min/max]
before production reordering. Each file also records its contributing prepared
files (basename, bytes, SHA-256), calendars, effective units, precipitation QC
statistics and source identities. Only the private storage path of each
contributing file was reduced to its basename; `ranges` and all other metadata
are identical to the originals (`identical()` checked in R). The original and
public file hashes are listed in `data/supports/manifest.csv`. Both supports pool
target-period years; see [KNOWN_ISSUES.md](KNOWN_ISSUES.md).

An optional exact-file runtime registry contains `roots` (named external roots)
and `inputs` (model, variable, root_id, safe relative_path, bytes, sha256,
optional first_date/last_date). The public catalogue is source metadata, not a
drop-in private-path runtime registry. For an active registry, the original
readiness contract additionally requires its SHA, appropriate support model/year
scope, matching support SHA and completed preparation state. Reconstruct these
from actual preparation; do not fabricate success receipts to bypass checks.

The accepted identities are:

- registry `c34724582dfbe7ed2f1ef73886bf83f88d758482e8f81922f70d6ca7a2a653b6`
- ERA5 support `5530eab3e4b67b5a26121614e08c52654ce55e7ab58790e097370ac0dd9e37ca`
- PME support `ca9394ce74fe78084e224116e37a196821782f31c1b6681592330ae3802058f4`

These identify original prepared artifacts including their metadata. They are
not access URLs. The support hashes identify the original objects, not the
path-redacted public copies in `data/supports/`, which carry the same arrays.
A readiness receipt built for the original registry therefore rejects the public
copies by checksum; a new registry/receipt must be produced from actual
preparation. The prepared climate files themselves are not included.

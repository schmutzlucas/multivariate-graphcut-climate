# Preparation of daily inputs

Use `data/metadata/source_manifest.csv` for expected source identities, effective
units, native calendars and original prepared-file checksums. All prepared files
must expose daily `pr`, `tas` or `psl` on 0:359 longitude and -90:90 latitude.
Production readers reorder longitude to -180:179 and align dates across variables.
No model-to-model daily synchronization is imposed. Preserve Gregorian, noleap,
all-leap and 360-day calendars; never manufacture a 31st in a 360-day year.

1. Select native daily CMIP6 historical and SSP585 data, matching variables/member
   and source versions. Concatenate historical/scenario segments with native
   calendar/time units, rejecting overlaps, duplicate or missing selected dates.
2. Form ERA5 daily means for tas/psl. For precipitation, sum accumulations stamped
   01–23 UTC of a day plus 00 UTC of the next day, expressed in mm/day.
3. Bilinearly remap to `config/grid_1degree.txt` using the original CDO `remapbil`
   convention. Convert tas to K, psl to Pa and pr to effective mm/day. For CMIP6
   flux in kg m-2 s-1 the daily factor is 86400 when not already applied.
4. Independently check variables, dimensions, dates, calendars, finite/missing
   values and effective units. Four native 360-day model precipitation products
   in the original collection already contain mm/day values despite stale flux
   attributes: HadGEM3-GC31-LL/MM, KACE-1-0-G and UKESM1-0-LL. Never reconvert
   merely because their original attribute says flux.
5. At PMF ingestion apply finite-negative precipitation QC to zero, then log1p;
   mean diagnostics use the nonnegative physical values. Preserve missing values
   and investigate large negative values/source units rather than hiding them.
6. Build separate common support arrays and PMFs with the exact scopes in
   `config/article1.R`; no per-reference PME support is substituted.

The implemented UTC month recipe is the following CDO operator chain (right to
left), with the native month followed by one midnight map from the next month:

```text
cdo -L -f nc4 -z zip_1 -setattribute,pr@units="mm day-1" -chname,tp,pr -remapbil,config/grid_1degree.txt -mulc,1000 -seldate,FIRST_DAY,LAST_DAY -daysum -shifttime,-1hour -mergetime MONTH.nc NEXT_MIDNIGHT.nc DAILY.nc
```

Replace the placeholders with the selected month/date range. The preserved
`tools/era5_daily_utc.R` provides `boundary`, `source_identity` and `verify` modes
on an explicit JSON job (`month` YYYYMM, `raw`, `following`, `boundary`, `product`,
`report`, `fingerprint_report`, `source_sha256`). `boundary` extracts only the
next midnight; `verify` compares daily sums independently at native/output grid
points and checks the whole product's dates/units/nonnegative finite values.
It never downloads source data. Full historical tas/psl preparation cannot be
recovered byte-for-byte from incomplete original receipts; this remains a
limitation of exact reconstruction, not evidence that current products are equal.

## Support construction

With prepared daily files in the hierarchy in DATA_REQUIREMENTS.md:

```text
Rscript --vanilla experiments/era5/build_supports.R
Rscript --vanilla experiments/perfect_model/build_supports.R
```

Set `MVGC_DATA_ROOT`, `MVGC_REBUILD_SUPPORTS=1`, `MVGC_OUTPUT_ROOT`, and a new
`MVGC_SUPPORT_OUTPUT_FILE` for each command. The shared support scope is
1951–1975 plus 1999–2023 for ERA5, and all 15 PME models over 1951–2100.
Common empirical minima/maxima are pooled per cell/variable; constant cells
retain their empirical minimum and use the documented upper-only numerical
envelope. Nonnegative pr is transformed with log1p. A new support produced from
different data or metadata is not the original support merely by assertion.

For signed per-input preparation, `tools/prepare_support_input.R JOB.json`
checks an explicit `file/model/variable/bytes/sha256`, grid/scopes, units and time
coverage and generates a new extrema cache/report. `tools/aggregate_prepared_supports.R`
pools those caches using explicit ERA5/PME scopes and signed identities. Its
comparison against older support is optional when no old support exists; the
pooling and constant-bound rules are unchanged. These expert job interfaces
retain their Roxygen2 contracts and refuse overwrite. Runtime paths/reports
belong in external data/output storage, not the publication repository.

## What the compact software archive distributes

Included: derived compact CSV/RDS fields and curves, accepted graphics/tables,
histogram supports, source-identity metadata and checksums. Not included: the original prepared
climate collection and large workspaces. The two exact production supports are
included in `data/supports/` (path-redacted copies with unchanged arrays).
Public upstream access routes do not establish exact reconstruction of historical
prepared products: missing acquisition/processing receipts, unavailable original
versions and effective-unit qualifications remain documented. Storage burden,
absence of an exact original copy, missing receipts and source-right restrictions
are separate issues. This archive does not assert a blanket legal prohibition
on sharing prepared data; specific rights must be checked against each provider.

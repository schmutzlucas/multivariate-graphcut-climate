# Full experiments from prepared inputs

Exact-input reproduction depends on the prepared files, the included supports
and the pinned solver; replication may use newly prepared source data. Neither is
executed by Level 1 reproduction. Configure the external paths in DATA_REQUIREMENTS.md,
choose new output directories, then run from the repository root:

```text
Rscript --vanilla experiments/era5/run.R
Rscript --vanilla experiments/era5/evaluate.R
Rscript --vanilla experiments/perfect_model/run.R
Rscript --vanilla experiments/perfect_model/evaluate.R
```

For the ERA5 evaluator set `ERA5_RESULTS_FILE` explicitly to the completed final
workspace. For PME evaluation set `PME_FINAL_RESULTS_FILE` to the completed
compact bundle. Do not select the newest arbitrary workspace. The ERA5 default
calibration/target are 1951–1975/1999–2023 with all 22 candidates. PME uses exactly
15 references, 14 other candidates each and five targets. The main article uses
.125; production defaults retain four lambdas .125/.25/.5/1. Set
`PME_LAMBDAS=0.125` for a main-only calculation, recording this deliberate
configuration reduction; all-lambda sensitivity is a different run scope.

`MVGC_WORKERS` / `MVGC_LAT_BLOCK_SIZE` bound ERA5 parallel extraction.
`PME_WORKERS`, `PME_OUTPUT_DIR`, `PME_CHECKPOINT_DIR` and `PME_RUN_STAMP` control
PME output/checkpoint management. Preserve original run stamps/settings when
resuming; invalid/mismatched scientific checkpoint identities fail. Never weaken
the gate or force overwrite to reuse incompatible PMFs or labels. Failed attempts
and logs stay in external storage. Reductions such as `PME_MAX_MODELS` and
`PME_PERIODS` are diagnostic subsets, not full article reproductions.

Sensitivities:

```text
Rscript --vanilla experiments/sensitivity/spatial_regularisation.R
Rscript --vanilla experiments/sensitivity/histogram_resolution.R
Rscript --vanilla experiments/sensitivity/ldr_threshold.R
```

Spatial sensitivity needs the explicitly selected completed ERA5 workspace and
computes actual solutions. Histogram sensitivity compares odd/even calibration
years for k=4:12 on the experiment support. LDR extension evaluates reference
bins on saved PME labels and prepared data, not a new optimizer. Use new output
paths; reduced grid/model/year settings do not reproduce the whole article.

The historical collection is about 803 GB prepared; a saved ERA5 intermediate
workspace was about 32 GiB. PMFs and multiple workers require substantial
additional memory/disk, and full experiments can take many hours or longer.
Those observations are not minimum RAM/time guarantees. Start with one worker
and a documented small prepared-input experiment on suitable hardware. The
included 4-by-2 native demo avoids those costs but starts from synthetic PMFs.
Fresh exact-input global reruns and complete source-to-prepared reconstruction
are not established by this publication package's compact rendering checks.

## Exact configuration

| Setting | Value (`config/article1.R`) |
| --- | --- |
| variables | `pr`, `tas`, `psl`; 1° grid 360×181 (`config/grid_1degree.txt`) |
| ERA5 candidates | 22 CMIP6 models, `config/era5_models.txt` |
| PME models | 15, `config/pme_models.txt`; each held out in turn, 14 candidates |
| calibration | 1951–1975 |
| ERA5 evaluation | 1999–2023 |
| PME targets | 1976–2000, 2001–2025, 2026–2050, 2051–2075, 2076–2100 |
| PMF bins per variable | 8 (trivariate); 16 (bivariate) and 32 (univariate) diagnostics |
| supports | `data/supports/era5_support.rds`, `data/supports/pme_support.rds` |
| λ | 0.125 (main); ERA5 sensitivity 0, 0.025, 0.125, 0.25, 0.5, 1 |
| LDR thresholds | ERA5 0.01, 0.05, 0.10; PME 0.005, 0.01, 0.025, 0.05, 0.10 |
| solver | `schmutzlucas/gcoWrapR@6b655d5feff7e7912d2ac0a7de95d0f8edd6d434`, alpha-beta swap, seed 1 |
| tested software | R 4.4.3, `config/software_versions.csv`, `renv.lock` |

Set `MVGC_ERA5_SUPPORT` and `MVGC_PME_SUPPORT` to the included supports (they are
not chosen implicitly). For claims and limits of each reproduction level see
[REPRODUCIBILITY.md](../REPRODUCIBILITY.md) and [KNOWN_ISSUES.md](KNOWN_ISSUES.md).

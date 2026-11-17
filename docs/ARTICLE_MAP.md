# Article figures and tables

The selected article is **Beyond mean fields: a multivariate graph-cut approach
to climate-model combination**, manuscript revision
`02e11666e50bfde8c6976c1a0abd654caa3ff2b8`, including its separate Supplement.

| Article component | Generation entry point | Compact numerical input |
| --- | --- | --- |
| Main Fig. 1: ERA5 distributional maps | `figures/era5_distributional.R` | ERA5 `distributional_fields.rds` |
| Main Fig. 2: ERA5 mean biases | `figures/era5_mean_bias.R` | ERA5 `mean_bias_fields.rds` |
| Main Fig. 3: PME full Hellinger densities | `figures/pme_overlays.R` | PME `full_overlay.rds` |
| Main Fig. 4: PME LDR densities, threshold .01 | `figures/pme_overlays.R` | PME `ldr_0p01_overlay.rds` |
| Main Fig. 5: connected reference improvements | `figures/pme_reference_improvements.R` | PME `boxplot_data.csv` |
| Main Fig. 6: spatial regularization | `figures/spatial_regularisation.R` | Spatial `label_maps.rds` |
| Appendix: histogram stability | `figures/histogram_resolution.R` | Histogram `cellwise.rds` |
| Appendix: four additional LDR thresholds | `figures/pme_overlays.R` | PME `ldr_*_overlay.rds` |
| Appendix A3: λ trade-off | `figures/lambda_tradeoff.R` | Spatial `lambda_tradeoff_metrics.csv` |
| Supplement Figs. S1–S5: reference-resolved densities | `figures/pme_overlays.R` | PME `reference_*.rds` |

Full project-relative input and accepted output paths are in
`publication/figure_manifest.csv`. The compact-input cohort directories are
selected by `reproduction_profile()`; only `current` is supported. The complete
workflow renders 17 PDF/PNG pairs.

`tables/render.R` generates the model-membership table, ERA5 global-score and
variable-subset tables, and two PME reference-wise tables. Inputs are the model
lists and `data/metadata/models.csv`, the ERA5 unrounded `summary_long.csv`, and
PME `reference_scores.csv`. Bibliographic entries for the model table are in
`data/metadata/model_references.bib`. The five accepted table fragments and their
manuscript labels are identified in `publication/table_manifest.csv`.

For the distinction between rendering compact results and executing the climate
experiments, see [REPRODUCIBILITY.md](../REPRODUCIBILITY.md).

Appendix A3: `figures/lambda_tradeoff.R` renders the six-row
`lambda_tradeoff_metrics.csv` without optimisation. Lowest tested Hellinger
distance is at lambda=0.025; the selected configuration remains 0.125.

All mapped figure bytes are identical to this manuscript revision. The complete
machine-readable mapping, including LaTeX labels, inputs, configuration,
verification status and limitations, is `publication/article_inventory.csv`.

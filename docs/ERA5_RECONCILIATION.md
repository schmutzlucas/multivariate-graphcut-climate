# ERA5 maps: history and resolution

Status: **resolved** in manuscript revision `a19496f` (MV_graphcut_overleaf).

Until Overleaf revision `f22ca16`, the manuscript's ERA5 maps (Figs. 1–2) were
rendered from research caches computed before the October 2026 data correction,
while Tables 2–3, Fig. 6, Fig. A3 and all quoted values came from the corrected
cohort. The maps were then re-rendered with the same renderers
(`figures/era5_distributional.R`, `figures/era5_mean_bias.R`) from the corrected
arrays in `results/era5/`. Page size, layout and captions are unchanged.

Differences between the superseded and current maps (area-weighted):

| Field | Superseded | Current |
| --- | --- | --- |
| D-MMM Hellinger distance (mean) | 0.24714 | 0.24792 |
| GC Hellinger distance (mean) | 0.18867 | 0.18936 |
| MMM temperature MAE | 0.868 K | 0.866 K |
| GC temperature MAE | 0.611 K | 0.604 K |
| GC precipitation MAE | 0.481 mm/day | 0.481 mm/day |

Baseline maps are essentially unchanged (spatial correlation ≥ 0.999). The GC
maps differ in patches where a different model was selected: 3,460 changed labels
at λ = 0.125 and 21 instead of 20 selected models. Cells changed by up to 0.45 in
Hellinger distance, 8.3 K and 9.0 mm/day. Regional means in every region named
in the text (tropical South America, Andes, East Africa, Maritime Continent, High
Asia, ITCZ) differ by at most 0.002 in Hellinger distance, so the qualitative
description in the article holds for both versions.

Evidence: `data/metadata/era5_asset_provenance.json` (current and superseded asset
hashes and cache identities), `era5_field_comparison.csv`, `era5_bias_comparison.csv`.

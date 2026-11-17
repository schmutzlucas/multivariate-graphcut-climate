# Changelog

## 1.0.0 (unreleased)

First public version accompanying the ESD manuscript *Beyond mean fields: a
multivariate graph-cut approach to climate-model combination*.

- Scientific code from research revision `69a6faf`, with final plotting code from
  `paper1-submission` `0671145` and the A3 renderer from `3f8cce7`.
- Compact numerical results of the corrected October 2026 cohort; all tables and
  in-text values equal the manuscript (`MV_graphcut_overleaf` `f22ca16`).
- 17 published figure PDF/PNG pairs and five table fragments, byte-identical to
  the manuscript assets.
- Exact production histogram supports (path-redacted copies).
- Machine-readable article inventory and public data manifest.
- Dependency installer falls back to the CRAN Archive for `ncdf4.helpers`/`PCICt`.
- CI: compact reproduction, tests and manifest checks; separate job installing the
  pinned solver for the Level 2 demonstration.

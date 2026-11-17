# Third-party software and data

The MIT licence applies to original author-owned research code and accompanying
documentation. It does not relicense external software, ERA5/CMIP6 datasets,
institutional/provider material or the native solver. Original published result
graphics/numerical products are supplied by the author for scientific reproduction;
underlying climate-data rights remain with their providers.

## Native graph-cut solver

The experiments require `schmutzlucas/gcoWrapR` at
`6b655d5feff7e7912d2ac0a7de95d0f8edd6d434`. The wrapper DESCRIPTION declares
GPL (>=2). Its bundled GCO version 3 notice credits Olga Veksler and Andrew
Delong and imposes research-use and redistribution restrictions. That notice is
not superseded here by the wrapper metadata. This repository does **not**
include or redistribute the solver's sources or binaries and cannot grant
permissions for them. Users obtain it from its own public repository
(`remotes::install_github("schmutzlucas/gcoWrapR@6b655d5feff7e7912d2ac0a7de95d0f8edd6d434")`)
and are responsible for complying with its terms. The continuous-integration job
installs it the same way for testing only.

The pinned fork implements periodic longitude. Substituting an ordinary rectangular
grid build changes the scientific graph and is unsupported.

The GCO notice requests these core citations:

- Boykov, Veksler and Zabih (2001), *Fast Approximate Energy Minimization via Graph Cuts*.
- Kolmogorov and Zabih (2004), *What Energy Functions Can Be Minimized via Graph Cuts?*
- Boykov and Kolmogorov (2004), *An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision*.

See the [pinned upstream notice](https://github.com/schmutzlucas/gcoWrapR/blob/6b655d5feff7e7912d2ac0a7de95d0f8edd6d434/inst/GCO_README.TXT).
Label-cost-feature-specific citations apply if that separate feature is used.
No claim about current patent validity is made from this historical notice.

## Other dependencies and data

Declared package licences/versions are in `config/dependency_licences.csv`.
R, Rcpp, CRAN packages and CDO retain their own terms; only calls/interfaces and
author-owned scientific kernels are included. Maps/coastline data are supplied
by the separately installed `maps` package under its distribution terms.

ERA5 and CMIP6 source collections and exact source identities are listed in
`docs/DATA_DOWNLOAD.md` and `data/metadata/source_manifest.csv`. Follow each
provider/model's actual access and licence terms before acquiring or redistributing
data. BCC/ACCESS source-version qualifications remain unresolved and explicit.

## Scope of the MIT licence

MIT covers original author-owned code and documentation in this repository. It
is not a licence for GCO, external climate files or provider material. The
compact numerical results and figures are derived from ERA5/CMIP6 data; users
should also acknowledge and comply with the source collections' terms.

Redistribution permission for the cylindrical GCO modification has not been
sought from the GCO rights holders. This repository avoids the question by not
redistributing the solver. Any future bundling would need that permission.

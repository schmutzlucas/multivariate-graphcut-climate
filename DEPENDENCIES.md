# Dependencies and installation

The scientific implementation was tested with R 4.4.3 on Windows; Python 3.12
supports the small Python tests. Rtools44 or another compatible C++ toolchain is
needed for Rcpp callbacks and the sensitivity kernel. Linux/macOS require the
usual R development tools and NetCDF libraries. Cross-platform native-label
bit identity is not promised for floating ties.

`config/software_versions.csv` records the original tested direct packages.
`renv.lock` records the available tested non-base direct/transitive R package
versions, without the restricted native solver. `config/dependency_licences.csv`
records their declared licence identifiers and versions. A snapshot is not
proof that all historical source distributions remain accessible today.

For version restoration, install `renv` explicitly, then in R:

```r
renv::restore(lockfile="renv.lock", library=Sys.getenv("MVGC_R_LIBRARY"), prompt=FALSE)
```

Set `MVGC_R_LIBRARY` to an empty external library directory and put it on
`R_LIBS_USER` for subsequent runs. Do not pass an empty path. If a recorded
package is unavailable, retain the failure and state the environment deviation;
do not claim that newer packages restore the original environment.
Alternatively, the opt-in `tools/install_dependencies.R --accept-current-cran`
installs missing CRAN packages into the caller's configured library and declares
that it creates a new environment. `ncdf4.helpers` (0.3-7) and `PCICt` (0.5-4.4)
were removed from CRAN in 2025; the installer takes these tested versions from
the CRAN Archive.

Rendering requires ggplot2, patchwork, maps and related packages. Climate workflows
add ncdf4, ncdf4.helpers, PCICt, future/future.apply, matrixStats, jsonlite and
openssl. Rcpp/RcppXPtrUtils support native scientific computation. roxygen2 is
for documentation generation; png supports pixel comparisons. The declared
licences of external packages govern their use; no external package source is
bundled under this project's MIT licence.

The exact longitude-periodic solver is **not installed automatically**. After
reading the terms in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md), install
`schmutzlucas/gcoWrapR` revision `6b655d5feff7e7912d2ac0a7de95d0f8edd6d434`
from its public repository with
`remotes::install_github("schmutzlucas/gcoWrapR@6b655d5feff7e7912d2ac0a7de95d0f8edd6d434")`. Runtime checks require that exact `RemoteSha`;
a rectangular-grid build with the same version number is not an equivalent
substitute. No restricted solver sources or binaries are in this archive.

Data acquisition requires the upstream providers' credentials/permissions,
separately managed by the researcher. Regridding recipes use CDO `remapbil`.
Exact original preprocessing installation/command receipts are incomplete;
fresh downloads/preparation are methodological replication until compared to
the recorded scientific input identities. No credentials belong in Git.

Figure rendering explicitly requires ragg and gtable. Documentation generation is
optional for scientific use; select it with `--with-docs`. CI validates both.
The installer checks every selected package after installation and fails on missing
packages. Linux source builds require libuv, libpng, FreeType, Fontconfig, HarfBuzz,
FriBidi, JPEG/TIFF and the NetCDF/OpenSSL/curl/XML development libraries declared
in the workflow. Tested package versions and renv.lock remain historical records;
current CRAN installation neither pins nor restores them.

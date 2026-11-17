# Scientific method and numerical conventions

Daily observations are summarized as local empirical joint PMFs. `pr` is
interpreted as effective mm/day, finite negative values are clipped to zero,
and `log1p` is applied only for distribution construction. Missing samples and
invalid/empty PMFs follow the checked original reader/normalization contracts.
Variables are ordered `pr,tas,psl`; longitude varies fastest in R grid arrays,
and bins vary fastest in native pairwise-cost arrays. Each joint PMF sums to one.

For normalized PMFs, Hellinger is
`sqrt(sum((sqrt(p)-sqrt(q))^2)/2)`, bounded by zero and one. Separate
univariate/bivariate diagnostics use 32/16 bins per variable, rather than
marginalizing the trivariate 8-bin PMF. Labels remain selected by the trivariate
objective. The distributional baseline D-MMM averages candidate PMFs with equal
model weights. The mean-state MMM averages candidate climatological means.

Unary cost is calibration-reference Hellinger times cosine latitude. A seam
between labels is lambda times the average endpoint cosine weight times the
sum of inter-model Hellinger distances at the two endpoints, using target-period
candidate PMFs. Four-neighbor longitude is periodic; latitude is not. The custom
native solver uses alpha-beta swap, initialized with seed 1 and L'Ecuyer-CMRG.
It finds an approximate local minimum for multiple labels, not a promised global
optimum. The graph retains polar rows at effectively zero area; it is cylindrical,
not a tessellation merging all polar vertices.

LDR bins are selected from the reference alone by ascending probability until
nominal mass is reached, including zeros and all cutoff ties. Restricted PMFs
are not renormalized. Main threshold is .01; sensitivities are .005/.025/.05/.10
in PME and .05/.10 in ERA5. Attained mass can exceed the threshold. These are
distributional diagnostics, not physical-event definitions of extremes.

Global scores use cosine-latitude weighting. Percentage gain is
`100*(baseline-GC)/baseline`, with a positive baseline. A ratio of pooled
area-mean scores differs from the arithmetic mean of 15 reference-wise gains;
both are supplied where appropriate. Reference quartiles use R type 7.
Scores are rounded only for display, not before gains are calculated.

ERA5 common support pools 23 datasets over 1951–1975 and 1999–2023 only.
PME support pools all 15 models over 1951–2100. Thus the held-out model is
excluded from candidate labels/baseline, but its future extrema can affect bin
edges. This shared data-informed support conditions the transfer results.

No daily time synchronization across models is assumed by distribution comparison.
A spatial patchwork of PMFs does not automatically supply coherent daily
space-time fields. Related CMIP6 models are not independent random replicates.

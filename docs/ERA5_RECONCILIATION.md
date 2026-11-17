# ERA5 article assets and numerical provenance

Archive-asset selection: **CONFIRMED BY THE AUTHOR**. The archive must contain
exactly the two ERA5 PDF/PNG pairs used in the article. Those accepted files are
retained in `publication/figures/main/era5_hellinger_maps.*` and `era5_bias_maps.*`.
The source revision and four SHA-256 identities are recorded in
`data/metadata/era5_asset_provenance.json` and checked by the publication tests.

Regenerated figures match 14/17 published pairs pixel for pixel; Fig. 5 differs
only by anti-aliasing; the two ERA5 pairs differ for the reason below. This asset
selection does not authorize replacing the corrected arrays or claim that they
produce the archived ERA5 figures. Their historical cohort difference remains
explicit below.

Both accepted ERA5 PDF rasters, every embedded image and their PNG previews match
renders from the older research caches exactly. With the corrected October
publication inputs, dimensions and text positions agree, but numerical fields
and pixels differ. The full Hellinger GC mean is 0.188669922715227 for the older
cache versus 0.189363224696502 for the corrected cohort. D-MMM is
0.247142215813498 versus 0.247919686851482. Maximum cell differences are
0.449670338510979 in GC Hellinger distance, 8.2978668820341 K in GC temperature
bias and 8.980547686517115 mm/day in GC precipitation bias.

Evidence: `data/metadata/era5_field_comparison.csv`, `era5_bias_comparison.csv`
and `era5_asset_provenance.json` identify exact arrays, sources and differences.
The accepted graphics trace to the older research cache bytes, not a newly
inferred producing workspace commit. The corrected workspace identity remains
in SCIENTIFIC_PROVENANCE.md. Underlying producing lineage is qualified.

Diagnosis: numerical arrays and historical workspace/cohort differ (B/C), with
model labelling differences evidenced by the historical review (E). The review
records 3,460 changed selected labels at lambda=0.125 and 20 versus 21 selected
models. Preprocessing/support effects cannot be isolated from these derived
fields alone (D remains unproven); this is not merely cosmetic (A). No new climate
experiment or label optimisation was run. A3 matches the corrected values.

Author decision: preserve the accepted article graphics in the prospective
Zenodo archive. Corrected-cohort renders remain separate generated outputs;
reproduction never overwrites the accepted assets. No Zenodo deposit is created
by this decision. Any future numerical-cohort reconciliation or manuscript change
is a separate scientific review, not an implicit consequence of archive selection.

Neither accepted graphics nor corrected scientific arrays have been replaced.
Rendering all 17 figure families does not by itself establish 17/17 agreement.

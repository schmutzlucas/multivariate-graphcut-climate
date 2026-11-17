The small `era5_production_patch.rds` contains cropped empirical PMFs from a
historical 4-by-2-cell, 22-model ERA5 experiment (source workspace dated
2026-09-13), together with stored Hellinger components. It is a reusable
component fixture, not a crop that reproduces the full global graph optimum,
and is not the accepted October result cohort.

`test_real_pmf_fixture.R` independently recomputes 192 stored distances. Source
location strings were removed; scientific arrays are unchanged. All accepted
figure inputs are under `results/`, selected through the `current` profile.

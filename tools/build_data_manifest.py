"""Write data/DATA_MANIFEST.csv: included, external-public, required-not-archived
and verification-pending data assets. Deterministic; run from any directory."""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/DATA_MANIFEST.csv"
FIELDS = ["category", "asset", "description", "bytes", "sha256", "provenance", "access", "status"]
COHORTS = {c["canonical_path"]: c for c in json.loads(
    (ROOT / "data/cohort_provenance.json").read_text(encoding="utf-8"))["cohorts"]}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def provenance(rel):
    for prefix, c in COHORTS.items():
        if rel.startswith(prefix + "/"):
            rev = c.get("research_revision", "")
            where = f"research {rev[:12]}:{c['research_path']}" if rev else c["research_path"]
            return f"{c['historical_cohort']}; {where}"
    if rel.startswith("publication/figures/"):
        return "manuscript asset, schmutzlucas/MV_graphcut_overleaf (see publication/figure_manifest.csv)"
    if rel.startswith("publication/tables/"):
        return "manuscript table fragment (see publication/table_manifest.csv)"
    if rel.startswith("tests/fixtures/"):
        return "cropped production PMF patch for regression tests (see tests/fixtures/README.md)"
    return "curated scientific metadata (see docs/SCIENTIFIC_PROVENANCE.md)"

DESCRIPTIONS = {
    "results/era5/": "ERA5 experiment scores/fields (corrected cohort)",
    "results/perfect_model/": "Perfect-model scores, boxplot data and density coordinates",
    "results/histogram_resolution/": "PMF-resolution (k=4..12) sensitivity results",
    "results/spatial_regularisation/": "Six-lambda label maps, fragmentation and A3 metrics",
    "data/supports/": "Exact production histogram support (path-redacted copy)",
    "data/metadata/": "Source identities, model metadata and ERA5 cohort comparison",
    "data/input_provenance.json": "Prepared-input identities of the corrected cohort",
    "data/cohort_provenance.json": "Research-source identities of compact cohorts",
    "publication/figures/": "Figure asset exactly as published in the manuscript",
    "publication/tables/": "Table fragment exactly as published in the manuscript",
    "tests/fixtures/": "Small production-derived regression fixture",
}

def describe(rel):
    return next((d for p, d in DESCRIPTIONS.items() if rel.startswith(p)), "")

rows = []
included = []
for base in ("results", "data/supports", "data/metadata", "publication/figures", "publication/tables", "tests/fixtures"):
    included += [p for p in (ROOT / base).rglob("*") if p.is_file() and p.name != "README.md"]
included += [ROOT / "data/input_provenance.json", ROOT / "data/cohort_provenance.json"]
for path in sorted(set(included)):
    rel = path.relative_to(ROOT).as_posix()
    rows.append(dict(category="1_included", asset=rel, description=describe(rel), bytes=path.stat().st_size,
                     sha256=sha(path), provenance=provenance(rel), access="this repository", status="verified checksum"))

rows += [
    dict(category="2_external_public", asset="CMIP6 climate projections (CDS)",
         description="Source collection for CMIP6 historical and SSP5-8.5 daily data", bytes="", sha256="",
         provenance="https://doi.org/10.24381/cds.c866074c", access="Copernicus Climate Data Store; registration and licence acceptance",
         status="collection DOI; not a DOI of the prepared files"),
    dict(category="2_external_public", asset="ERA5 hourly data on single levels (CDS)",
         description="Source collection for ERA5 pr/tas/psl", bytes="", sha256="",
         provenance="https://doi.org/10.24381/cds.adbb2d47", access="Copernicus Climate Data Store; registration and licence acceptance",
         status="collection DOI; not a DOI of the prepared files"),
]
with (ROOT / "data/metadata/source_manifest.csv").open(encoding="utf-8", newline="") as stream:
    for r in csv.DictReader(stream):
        rows.append(dict(category="3_required_not_archived", asset=r["prepared_filename"],
                         description=f"Prepared daily {r['variable']} for {r['model']} on the 1-degree grid ({r['calendar']}, {r['effective_units']})",
                         bytes=r["prepared_bytes"], sha256=r["prepared_sha256"],
                         provenance="derived from the source identities in data/metadata/source_manifest.csv",
                         access="not publicly archived; access route not yet decided by the author",
                         status="required for Level 3 exact-input reproduction"))
rows += [
    dict(category="3_required_not_archived", asset="ERA5 experiment workspace and PME result bundles/checkpoints",
         description="Global PMFs, labels and intermediate optimisation outputs (tens of GiB)", bytes="", sha256="",
         provenance="identities in docs/SCIENTIFIC_PROVENANCE.md", access="not publicly archived",
         status="needed only to re-evaluate without rerunning; compact results are included"),
    dict(category="4_verification_pending", asset="Redistribution of the 69 prepared climate files",
         description="Provider terms (CDS/CMIP6 model licences) for redistributing derived daily files", bytes="", sha256="",
         provenance="", access="", status="permission and persistent hosting not established"),
    dict(category="4_verification_pending", asset="ACCESS-CM2 and BCC-CSM2-MR historical source versions",
         description="Version/replacement history qualifications", bytes="", sha256="",
         provenance="data/metadata/source_manifest.csv", access="", status="unresolved upstream; documented"),
    dict(category="4_verification_pending", asset="ERA5 original CDS request receipts",
         description="Incomplete hourly-month request/aggregation receipts", bytes="", sha256="",
         provenance="data/input_provenance.json", access="", status="cannot be recovered; documented"),
]
with OUT.open("w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
counts = {}
for r in rows:
    counts[r["category"]] = counts.get(r["category"], 0) + 1
print("PASS: data manifest written", counts)

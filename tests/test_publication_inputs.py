"""Reusable scientific-input and frozen-publication regression contracts."""
import csv
import hashlib
import json
from pathlib import Path
import statistics
import unittest

ROOT = Path(__file__).resolve().parents[1]


def records(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


class PublicationInputs(unittest.TestCase):
    def test_compact_manifests(self):
        manifests = [
            "results/era5/manifest.csv",
            "results/perfect_model/manifest.csv",
            "results/spatial_regularisation/manifest.csv",
            "results/histogram_resolution/manifest.csv",
            "data/supports/manifest.csv",
        ]
        for manifest in manifests:
            for row in records(manifest):
                path = ROOT / Path(manifest).parent / row["file"]
                raw = path.read_bytes()
                self.assertEqual(len(raw), int(row["bytes"]), str(path))
                self.assertEqual(hashlib.sha256(raw).hexdigest(), row["sha256"])

    def test_pme_population_and_gain_formula(self):
        rows = records("results/perfect_model/evaluation/reference_scores.csv")
        models = (ROOT / "config/pme_models.txt").read_text().split()
        periods = ("1976_2000", "2001_2025", "2026_2050", "2051_2075", "2076_2100")
        metrics = ("full", "0.005", "0.01", "0.025", "0.05", "0.10")
        self.assertEqual(len(models), 15)
        self.assertEqual(len(rows), 450)
        self.assertEqual({(r["reference"], r["period_id"], r["metric"]) for r in rows},
                         {(m, p, d) for m in models for p in periods for d in metrics})
        for row in rows:
            baseline, candidate = float(row["H_D_MMM"]), float(row["H_GC"])
            self.assertGreater(baseline, 0)
            self.assertAlmostEqual(100 * (baseline-candidate)/baseline,
                                   float(row["relative_improvement_pct"]), places=10)

    def test_pooled_and_reference_aggregation(self):
        rows = records("results/perfect_model/evaluation/reference_scores.csv")
        sample = [r for r in rows if r["period_id"] == "1976_2000" and r["metric"] == "full"]
        baseline = statistics.mean(float(r["H_D_MMM"]) for r in sample)
        candidate = statistics.mean(float(r["H_GC"]) for r in sample)
        pooled = 100 * (baseline-candidate)/baseline
        average_gain = statistics.mean(float(r["relative_improvement_pct"]) for r in sample)
        self.assertAlmostEqual(pooled, 25.9482662133142, places=10)
        self.assertAlmostEqual(average_gain, 24.6714787661649, places=10)
        self.assertNotAlmostEqual(pooled, average_gain, places=2)

    def test_frozen_figure_selection(self):
        rows = records("publication/figure_manifest.csv")
        self.assertEqual(len(rows), 17)
        self.assertEqual(len({r["publication_pdf"] for r in rows}), 17)
        for row in rows:
            pdf = ROOT / row["publication_pdf"]
            self.assertEqual(hashlib.sha256(pdf.read_bytes()).hexdigest(), row["sha256_pdf"])
            png = ROOT / row["publication_png"]
            self.assertEqual(hashlib.sha256(png.read_bytes()).hexdigest(), row["sha256_png"])
        self.assertTrue(any("pme_reference_improvement.pdf" in r["publication_pdf"] for r in rows))
        self.assertFalse(any("boxplots_only" in r["publication_pdf"] for r in rows))

    def test_era5_article_assets(self):
        ledger = json.loads((ROOT / "data/metadata/era5_asset_provenance.json").read_text(encoding="utf-8"))
        assets = ledger["current_assets"]["assets"]
        self.assertEqual(len(assets), 4)
        expected = {f"publication/figures/main/{stem}.{ext}" for stem in
                    ("era5_hellinger_maps", "era5_bias_maps") for ext in ("pdf", "png")}
        self.assertEqual({a["archive_path"] for a in assets}, expected)
        for asset in assets:
            self.assertEqual(hashlib.sha256((ROOT / asset["archive_path"]).read_bytes()).hexdigest(),
                             asset["sha256"], asset["archive_path"])

    def test_five_tables_and_raw_scores(self):
        metadata = records("data/metadata/models.csv")
        self.assertEqual([r["model"] for r in metadata],
                         (ROOT / "config/era5_models.txt").read_text().split())
        self.assertEqual({r["model"] for r in metadata if r["pme"] == "yes"},
                         set((ROOT / "config/pme_models.txt").read_text().split()))
        bibliography = (ROOT / "data/metadata/model_references.bib").read_text(encoding="utf-8")
        for model in metadata:
            self.assertIn("{" + model["bibtex_key"] + ",", bibliography)
        rows = records("publication/table_manifest.csv")
        self.assertEqual(len(rows), 5)
        for row in rows:
            self.assertTrue((ROOT / row["accepted_latex"]).is_file())
            self.assertIn(row["table_label"], (ROOT / row["accepted_latex"]).read_text(encoding="utf-8"))
        era = records("results/era5/summary_long.csv")
        self.assertEqual(len(era), 129)
        gain = next(float(r["metric_value"]) for r in era if
                    r["metric_name"] == "improvement_pct_vs_MMM" and r["variable_set"] == "pr" and r["metric_type"] == "full_H")
        self.assertAlmostEqual(gain, -21.4499687270578, places=10)

    def test_input_scope_and_qualifications(self):
        ledger = json.loads((ROOT / "data/input_provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(len(ledger["inputs"]), 69)
        self.assertEqual(ledger["registry_sha256"], "c34724582dfbe7ed2f1ef73886bf83f88d758482e8f81922f70d6ca7a2a653b6")
        self.assertEqual(ledger["supports"]["era5"]["sha256"], "5530eab3e4b67b5a26121614e08c52654ce55e7ab58790e097370ac0dd9e37ca")
        self.assertEqual(ledger["supports"]["pme"]["sha256"], "ca9394ce74fe78084e224116e37a196821782f31c1b6681592330ae3802058f4")
        self.assertEqual(ledger["support_scopes"]["pme"]["years"], list(range(1951, 2101)))
        for item in ledger["inputs"]:
            if item["model"] in ("BCC-CSM2-MR", "ACCESS-CM2"):
                self.assertTrue(item["provenance_qualifications"])


if __name__ == "__main__":
    unittest.main()

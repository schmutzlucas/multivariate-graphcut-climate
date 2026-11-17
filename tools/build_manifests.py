"""Deterministic inventories; exclude both inventory files from their own hashes.
Run after documentation generation, before tests/reproduction. Scientific payloads
are never written. Add intended archive files to Git before running this tool.
"""
from pathlib import Path
import csv
import hashlib
import subprocess
ROOT = Path(__file__).resolve().parents[1]
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
figures = list(csv.DictReader((ROOT / "publication/figure_manifest.csv").open()))
for row in figures:
    for extension in ("pdf", "png"):
        assert sha(ROOT / row["publication_" + extension]) == row["sha256_" + extension], row
# Existing scientific manifests are verified, never silently rebaselined.
for manifest in sorted((ROOT / "results").glob("*/manifest.csv")) + [ROOT / "data/supports/manifest.csv"]:
    for row in csv.DictReader(manifest.open()):
        path = manifest.parent / row["file"]
        assert path.stat().st_size == int(row["bytes"]) and sha(path) == row["sha256"], path
old = {row["path"]: row["category"] for row in csv.DictReader((ROOT / "publication/release_manifest.csv").open())}
files = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True).splitlines()
rows = []
for name in sorted(files):
    if name in ("publication/release_manifest.csv", "publication/SHA256SUMS"):
        continue
    path = ROOT / name
    category = old.get(name, "scientific_metadata" if name.startswith("data/") else "maintained_workflow")
    rows.append(dict(path=name, category=category, bytes=path.stat().st_size, sha256=sha(path)))
write_csv(ROOT / "publication/release_manifest.csv", rows)
(ROOT / "publication/SHA256SUMS").write_text("".join(f"{r['sha256']}  {r['path']}\n" for r in rows), encoding="utf-8", newline="\n")
print(f"PASS: {len(rows)} archive identities; self-inventories excluded")

"""Static release audit and deterministic SHA-256 inventory."""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import re
from pathlib import Path

from _bootstrap import METADATA, ROOT
from run_pipeline import STAGES

EXCLUDED_PARTS = {".git", "__pycache__", "library", "staging", "run_logs"}
GENERATED_TOP_LEVEL = {"results", "figures", "metadata"}
PUBLISHED_DATA_FILES = {
    "data/README.md",
    "data/input_manifest.tsv",
    "data/GPL13497_probe2symbol.tsv",
}
FORBIDDEN = ("D:/R/projects/piezo1_oa", "E:/codex/2026", "C:/Users/Administrator")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def inventory() -> int:
    output = METADATA / "SHA256SUMS.txt"
    rows: list[str] = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == output:
            continue
        relative = path.relative_to(ROOT)
        if any(part in EXCLUDED_PARTS for part in relative.parts):
            continue
        if relative.parts[0] in GENERATED_TOP_LEVEL:
            continue
        if relative.parts[0] == "data" and relative.as_posix() not in PUBLISHED_DATA_FILES:
            continue
        rows.append(f"{digest(path)}  {relative.as_posix()}")
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return len(rows)


def main() -> int:
    checks: list[tuple[str, bool, str]] = []
    authoritative = sorted({name for values in STAGES.values() for name in values})
    for name in authoritative:
        path = ROOT / "scripts" / name
        checks.append((f"script:{name}", path.exists(), "authoritative script exists"))
        if path.exists():
            text = path.read_text(encoding="utf-8-sig")
            checks.append((f"portable:{name}", not any(token in text for token in FORBIDDEN), "no workstation absolute path"))
            if path.suffix == ".py":
                try:
                    ast.parse(text)
                    valid_python = True
                except SyntaxError:
                    valid_python = False
                checks.append((f"python-syntax:{name}", valid_python, "AST parses"))
    for name in (".zenodo.json", "codemeta.json"):
        try:
            json.loads((ROOT / name).read_text(encoding="utf-8"))
            valid_json = True
        except (OSError, json.JSONDecodeError):
            valid_json = False
        checks.append((f"json:{name}", valid_json, "valid JSON"))
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    checks.append(("citation-version", "version: 2.0.0" in cff, "CFF version matches release"))
    checks.append(("citation-orcid", "0009-0003-5907-7037" in cff, "ORCID present"))
    checks.append(("doi-not-fabricated", not re.search(r"10\.5281/zenodo\.\d+", cff), "no unassigned DOI"))
    published = ROOT / "reference" / "figures" / "published"
    checks.append(("main-figures", len(list(published.glob("Figure[1-5].png"))) == 5, "five main figures"))
    checks.append(("supp-figures", len(list(published.glob("Figure_S*.png"))) == 15, "fifteen supplementary figures"))
    checks.append(("reference-results", len(list((ROOT / "reference" / "results").glob("*.csv"))) >= 70, "frozen result tables present"))
    with (ROOT / "data" / "input_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    checks.append(("input-manifest", len(manifest) >= 20, "raw/derived input inventory present"))
    bundled = ROOT / "data" / "GPL13497_probe2symbol.tsv"
    expected = next(row["sha256"] for row in manifest if row["local_path"] == "data/GPL13497_probe2symbol.tsv")
    checks.append(("bundled-probe-map", bundled.exists() and digest(bundled).lower() == expected.lower(), "bundled derivative checksum"))
    count = inventory()
    checks.append(("sha256-inventory", count >= 100, f"{count} files inventoried"))
    failures = [row for row in checks if not row[1]]
    lines = ["PIEZO1 OA reproducible repository v2.0.0 audit", ""]
    lines.extend(f"{'PASS' if ok else 'FAIL'}\t{name}\t{detail}" for name, ok, detail in checks)
    lines.extend(["", f"FINAL STATUS: {'PASS' if not failures else 'FAIL'}"])
    METADATA.mkdir(parents=True, exist_ok=True)
    (METADATA / "repository_audit.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

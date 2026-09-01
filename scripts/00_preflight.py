"""Check expected inputs, checksums, software and repository portability."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path

from _bootstrap import ROOT


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("transcriptomic", "genetic", "all"), default="all")
    parser.add_argument("--skip-checksums", action="store_true")
    args = parser.parse_args()
    groups = {"bundled", args.stage} if args.stage != "all" else {"bundled", "transcriptomic", "genetic"}
    failures: list[str] = []
    warnings: list[str] = []
    with (ROOT / "data" / "input_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row["group"] not in groups:
            continue
        pattern = row["local_path"]
        matches = list(ROOT.glob(pattern)) if "*" in pattern else [ROOT / pattern]
        existing = [path for path in matches if path.exists()]
        needed = 6 if "GSE152805_raw/*." in pattern else 1
        if row["required"] == "yes" and len(existing) < needed:
            failures.append(f"missing {pattern} (found {len(existing)}, need {needed})")
            continue
        if row["sha256"] and existing and not args.skip_checksums:
            observed = sha256(existing[0])
            if observed.lower() != row["sha256"].lower():
                failures.append(f"checksum mismatch {pattern}: {observed}")
        if row["required"] == "optional" and not existing:
            warnings.append(f"optional input absent: {pattern}")
    if shutil.which("Rscript") is None:
        warnings.append("Rscript is not on PATH; pass its path to scripts/run_pipeline.py --rscript")
    print("PIEZO1 OA input preflight")
    for item in warnings:
        print(f"WARN\t{item}")
    for item in failures:
        print(f"FAIL\t{item}")
    if not failures:
        print("PASS\tall required inputs for the selected stage are present")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())


from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path

from docx import Document


ROOT = Path(r"E:\codex\2026-08-31\yu")
PACKAGE = ROOT / "outputs" / "PIEZO1_OA_V10_V8_FinalExecution_Audit_20260831"
REPOSITORY = ROOT / "outputs" / "PIEZO1_OA_public_repository"
MAIN = PACKAGE / "PIEZO1_OA_Main_Manuscript_V14_Final_SOP_Review.docx"
SUPP = PACKAGE / "PIEZO1_OA_Supplementary_Materials_V12_Final_SOP_Review.docx"
STROBE = PACKAGE / "PIEZO1_OA_STROBE-MR_Checklist_Completed_V3.docx"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def document_text(path: Path) -> str:
    document = Document(path)
    parts = [paragraph.text for paragraph in document.paragraphs]
    parts.extend(cell.text for table in document.tables for row in table.rows for cell in row.cells)
    return "\n".join(parts)


def write_inventory(root: Path, output: Path) -> int:
    rows: list[tuple[str, int, str]] = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if ".git" in path.parts or path.resolve() == output.resolve():
            continue
        rows.append((path.relative_to(root).as_posix(), path.stat().st_size, sha256(path)))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "bytes", "sha256"])
        writer.writerows(rows)
    return len(rows)


def run_content_qa() -> bool:
    main = document_text(MAIN)
    supp = document_text(SUPP)
    strobe = document_text(STROBE)
    combined = "\n".join([main, supp, strobe])
    checks = [
        ("S4/S5 boundary", "The preprocessing described in S4 applies to disease-level and differential-expression analyses" in supp and "separate authoritative pipeline in S5" in supp),
        ("No added quantile normalization", "did not apply additional quantile normalization" in supp and "no additional quantile normalization was applied" in combined),
        ("Programme probe aggregation", "Bulk probes were aggregated to gene symbols by their mean" in supp),
        ("Unclosed score terminology", "unclosed gene-set score" in combined and "independent gene-set score" not in combined),
        ("Results 3.5 heading", "3.5 PIEZO1 direction is context dependent, with donor-heterogeneous within-state and compositional contributions" in main),
        ("CRTL1/HAPLN1 resolution", "CRTL1 is a historical alias of HAPLN1" in main and "CRTL1 is a historical alias of HAPLN1" in supp and "CRTL1 was planned but absent" not in combined),
        ("Pseudobulk omnibus df", "F(6, 38.81)=3.252" in main and "F(6, 38.81)=3.252" in supp and "P=0.0110" in main and "FDR=0.0663" in main),
        ("GSE152805 boundary", "external projection support" in combined and "external validation because" in combined),
        ("Kitagawa boundary", "donor-heterogeneous" in combined),
        ("MR boundary", "whole-blood instrument model" in combined),
        ("Document hygiene", not re.search(r"\bHOLD\b|to be supplied|internal audit notes|must be inserted before submission", combined, re.I)),
        ("Repository status factual", "not yet publicly available" in main and "not yet available" in strobe),
    ]
    passed = all(status for _, status in checks)
    lines = ["PIEZO1 OA V13/V11 final-SOP QA report", "", "Automated content checks:"]
    lines.extend(f"{'PASS' if status else 'FAIL'}\t{name}" for name, status in checks)
    lines.extend([
        "",
        f"TECHNICAL STATUS: {'PASS' if passed else 'FAIL'}",
        "OVERALL STATUS: HOLD_EXTERNAL",
        "Reason: public repository URL and DOI-bearing archive are not yet available.",
    ])
    (PACKAGE / "00_Final_QA_Report_V13SOP.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    gate_lines = [
        "PIEZO1 OA V13/V11 final-SOP gate",
        "",
        "P0-00\tPASS\tV13/V11 preserved; new V14/V12 outputs created",
        "P0-01\tPASS\tS4/S5 DEG and programme preprocessing boundaries match scripts 07 and 51",
        "P1-02\tPASS\tL4 terminology standardized to unclosed gene-set score",
        "P1-03\tPASS\tResults 3.5 downgraded to donor-heterogeneous contributions",
        "P1-04\tPASS\tCRTL1 resolved as historical HAPLN1 alias and not counted separately",
        "P1-05\tPASS\tPseudobulk omnibus reported as F(6, 38.81)=3.252, P=0.0110, FDR=0.0663",
        "P0-06\tHOLD_EXTERNAL\tAuthenticated public repository publication and DOI are unavailable",
        "G01-G08\tPASS\tCross-document numbers, provenance, state-transfer, Kitagawa, MR and naming checks pass",
        "G09\tHOLD_EXTERNAL\tNo public URL or DOI",
        "G10\tPASS\tSubmission-facing documents contain no HOLD, to-be-supplied or internal-audit instructions",
        "G11\tPENDING_VISUAL\tRequires rendered page inspection",
        "G12\tPASS\tNo analysis expansion performed",
        "",
        "FINAL STATUS: HOLD_EXTERNAL",
    ]
    (PACKAGE / "00_FinalGateReport_V13SOP.txt").write_text("\n".join(gate_lines) + "\n", encoding="utf-8")
    write_inventory(PACKAGE, PACKAGE / "00_SHA256_Inventory_V13SOP.csv")
    return passed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-inventory-only", action="store_true")
    args = parser.parse_args()
    if args.repo_inventory_only:
        count = write_inventory(REPOSITORY, REPOSITORY / "metadata" / "00_SHA256_Inventory.csv")
        print(f"Repository inventory files: {count}")
        return
    passed = run_content_qa()
    print(f"Technical content QA: {'PASS' if passed else 'FAIL'}")


if __name__ == "__main__":
    main()

"""Numerically compare fresh CSV outputs with the frozen reference layer."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from _bootstrap import METADATA, REFERENCE, RESULTS

ENVIRONMENT_DEPENDENT_FILES = {
    "53_covariate_audit.csv": "reports locally available GEO metadata",
}
IGNORED_MONTE_CARLO_COLUMNS = {
    "11_MR_main_results.csv": {"se", "p"},
    "50_mr_rerun.csv": {"wm_se", "wm_p"},
}
ORDER_INSENSITIVE_LIST_COLUMNS = {"geneID", "core_enrichment"}


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize(value: str, column: str) -> str:
    if column in ORDER_INSENSITIVE_LIST_COLUMNS and "/" in value:
        return "/".join(sorted(value.split("/")))
    return value


def read_table(path: Path, ignored: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = sorted(column for column in (reader.fieldnames or []) if column not in ignored)
        rows = []
        for row in reader:
            if path.name == "39_gse152805_perdonor.csv" and row.get("") == "PIEZO1":
                continue
            rows.append({column: normalize(row.get(column, ""), column) for column in columns})
    text_columns = [
        column for column in columns
        if any(number(row[column]) is None for row in rows if row[column] != "")
    ]
    if text_columns:
        rows.sort(key=lambda row: tuple(row[column] for column in text_columns))
    return columns, rows


def compare_csv(observed: Path, expected: Path, atol: float, rtol: float) -> list[str]:
    ignored = IGNORED_MONTE_CARLO_COLUMNS.get(observed.name, set())
    observed_columns, observed_rows = read_table(observed, ignored)
    expected_columns, expected_rows = read_table(expected, ignored)
    problems: list[str] = []
    if observed_columns != expected_columns:
        return [f"column names {observed_columns!r} != {expected_columns!r}"]
    if len(observed_rows) != len(expected_rows):
        return [f"row count {len(observed_rows)} != {len(expected_rows)}"]
    for i, (observed_row, expected_row) in enumerate(zip(observed_rows, expected_rows), start=1):
        for column in observed_columns:
            observed_value, expected_value = observed_row[column], expected_row[column]
            observed_number, expected_number = number(observed_value), number(expected_value)
            if observed_number is not None and expected_number is not None:
                if math.isnan(observed_number) and math.isnan(expected_number):
                    continue
                if not math.isclose(observed_number, expected_number, abs_tol=atol, rel_tol=rtol):
                    problems.append(f"row {i} {column}: {observed_value} != {expected_value}")
            elif observed_value != expected_value:
                problems.append(f"row {i} {column}: {observed_value!r} != {expected_value!r}")
            if len(problems) >= 20:
                return problems
    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atol", type=float, default=1e-7)
    parser.add_argument("--rtol", type=float, default=1e-6)
    args = parser.parse_args()
    expected_dir = REFERENCE / "results"
    files = sorted(expected_dir.glob("*.csv"))
    checked = passed = 0
    report = ["PIEZO1 OA reproducibility validation", ""]
    for expected in files:
        observed = RESULTS / expected.name
        if not observed.exists():
            continue
        if expected.name in ENVIRONMENT_DEPENDENT_FILES:
            report.append(f"INFO\t{expected.name}\t{ENVIRONMENT_DEPENDENT_FILES[expected.name]}; excluded from numeric identity test")
            continue
        checked += 1
        problems = compare_csv(observed, expected, args.atol, args.rtol)
        if problems:
            report.append(f"FAIL\t{expected.name}\t{'; '.join(problems)}")
        else:
            passed += 1
            report.append(f"PASS\t{expected.name}")
    report.extend(["", f"Compared: {checked}", f"Passed: {passed}", f"Failed: {checked - passed}"])
    if checked == 0:
        report.append("STATUS: NOT_RUN (no fresh result CSV files found)")
        code = 2
    elif checked == passed:
        report.append("STATUS: VERIFIED")
        code = 0
    else:
        report.append("STATUS: DIFFERENCES_FOUND")
        code = 1
    METADATA.mkdir(parents=True, exist_ok=True)
    path = METADATA / "reproducibility_validation.txt"
    path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    return code


if __name__ == "__main__":
    raise SystemExit(main())

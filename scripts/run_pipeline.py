"""Deterministic, fail-fast orchestrator for the authoritative analysis."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from datetime import datetime, timezone

from _bootstrap import METADATA, ROOT

STAGES = {
    "transcriptomic": [
        "10_scrna_piezo1.R",
        "07_deg_paired.R",
        "33_deconvolution_mechano.R",
        "35_replication_gse57218.R",
        "39_gse152805.R",
        "40_P0_donorlevel_inference.R",
        "51_P1_compositional_sensitivity.R",
        "52_P1_gse152805_within_state.R",
        "53_P1_covariate_audit.py",
        "54_P1_gse51588_covariate_sensitivity.R",
    ],
    "genetic": [
        "11_mr_piezo1.R",
        "11b_coloc_piezo1.R",
        "41_P0_coloc_PPH3.R",
        "50_P1_mr_exposure_audit.R",
        "50b_P1_gtex_scale_validation.py",
        "32_power_primary_mr.py",
    ],
    "validation": [
        "71_v10_controlled_validations.R",
        "80_final_execution_nested_validation.R",
        "83_final_audit_tables.R",
    ],
    "figures": [
        "70_rebuild_figure2_v9.R",
        "75_rebuild_figure3_v10.R",
        "72_rebuild_figure4_v10.R",
        "73_rebuild_figure5_v10.R",
        "74_rebuild_supplementary_figures_v8.R",
        "81_rebuild_figure_s15_from_audit.R",
        "82_rebuild_figure1_final.R",
        "86_rebuild_supp_figure_s1_final.R",
    ],
}


def command_for(name: str, rscript: str) -> list[str]:
    path = ROOT / "scripts" / name
    return [rscript, str(path)] if path.suffix.lower() == ".r" else [sys.executable, str(path)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=(*STAGES, "all"), default="all")
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    selected = list(STAGES) if args.stage == "all" else [args.stage]
    scripts = [(stage, name) for stage in selected for name in STAGES[stage]]
    for stage, name in scripts:
        print(f"{stage:14s} {' '.join(command_for(name, args.rscript))}")
    if args.dry_run:
        return 0
    env = os.environ.copy()
    env["PIEZO1_PROJECT_ROOT"] = str(ROOT)
    log_dir = METADATA / "run_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    for stage, name in scripts:
        log_path = log_dir / f"{stamp}_{name}.log"
        print(f"[run] {stage}/{name}")
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                command_for(name, args.rscript), cwd=ROOT, env=env,
                stdout=log, stderr=subprocess.STDOUT, text=True, check=False,
            )
        if completed.returncode:
            print(f"[fail] {name}; see {log_path}", file=sys.stderr)
            return completed.returncode
        print(f"[pass] {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


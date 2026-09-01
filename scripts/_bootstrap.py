"""Shared path resolution for portable PIEZO1 OA scripts."""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("PIEZO1_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA = Path(os.environ.get("PIEZO1_DATA_DIR", ROOT / "data")).resolve()
RESULTS = Path(os.environ.get("PIEZO1_RESULTS_DIR", ROOT / "results")).resolve()
FIGURES = Path(os.environ.get("PIEZO1_FIGURES_DIR", ROOT / "figures")).resolve()
METADATA = Path(os.environ.get("PIEZO1_METADATA_DIR", ROOT / "metadata")).resolve()
REFERENCE = ROOT / "reference"

for path in (RESULTS, FIGURES, METADATA):
    path.mkdir(parents=True, exist_ok=True)

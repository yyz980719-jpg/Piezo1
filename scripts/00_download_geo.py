"""Download the public GEO inputs used by the authoritative workflow."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import tarfile
import urllib.request
from pathlib import Path

from _bootstrap import DATA

BASE_SERIES = "https://ftp.ncbi.nlm.nih.gov/geo/series"
BASE_PLATFORM = "https://ftp.ncbi.nlm.nih.gov/geo/platforms"

DOWNLOADS = {
    DATA / "GSE51588_series_matrix.txt.gz": f"{BASE_SERIES}/GSE51nnn/GSE51588/matrix/GSE51588_series_matrix.txt.gz",
    DATA / "GSE55235_series_matrix.txt.gz": f"{BASE_SERIES}/GSE55nnn/GSE55235/matrix/GSE55235_series_matrix.txt.gz",
    DATA / "GSE57218_series_matrix.txt.gz": f"{BASE_SERIES}/GSE57nnn/GSE57218/matrix/GSE57218_series_matrix.txt.gz",
    DATA / "GSE82107_series_matrix.txt.gz": f"{BASE_SERIES}/GSE82nnn/GSE82107/matrix/GSE82107_series_matrix.txt.gz",
    DATA / "GSE46750_series_matrix.txt.gz": f"{BASE_SERIES}/GSE46nnn/GSE46750/matrix/GSE46750_series_matrix.txt.gz",
    DATA / "GPL6947.annot.gz": f"{BASE_PLATFORM}/GPL6nnn/GPL6947/annot/GPL6947.annot.gz",
    DATA / "GSE104782_series_matrix.txt.gz": f"{BASE_SERIES}/GSE104nnn/GSE104782/matrix/GSE104782_series_matrix.txt.gz",
    DATA / "scrna" / "GSE104782_allcells_UMI_count.txt.gz": f"{BASE_SERIES}/GSE104nnn/GSE104782/suppl/GSE104782_allcells_UMI_count.txt.gz",
    DATA / "scrna" / "GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx": f"{BASE_SERIES}/GSE104nnn/GSE104782/suppl/GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx",
}
LARGE_TARGET = DATA / "GSE152805_RAW.tar"
LARGE_URL = f"{BASE_SERIES}/GSE152nnn/GSE152805/suppl/GSE152805_RAW.tar"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size:
        print(f"[keep] {target.relative_to(DATA.parent)}")
        return
    temp = target.with_suffix(target.suffix + ".part")
    print(f"[download] {url}")
    with urllib.request.urlopen(url, timeout=120) as response, temp.open("wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    temp.replace(target)


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            target = (root / member.name).resolve()
            if root != target and root not in target.parents:
                raise RuntimeError(f"unsafe archive member: {member.name}")
        bundle.extractall(root, filter="data")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-large", action="store_true", help="also download and extract the 314-MB GSE152805 archive")
    args = parser.parse_args()
    for target, url in DOWNLOADS.items():
        download(url, target)
        print(f"[sha256] {target.name} {sha256(target)}")
    if args.include_large:
        download(LARGE_URL, LARGE_TARGET)
        print(f"[sha256] {LARGE_TARGET.name} {sha256(LARGE_TARGET)}")
        safe_extract(LARGE_TARGET, DATA / "GSE152805_raw")
        print("[extract] data/GSE152805_raw")
    else:
        print("[skip] GSE152805_RAW.tar; rerun with --include-large for the full workflow")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


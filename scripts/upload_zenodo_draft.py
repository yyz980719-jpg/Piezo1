#!/usr/bin/env python3
"""Verify and upload repository artifacts to an existing Zenodo draft.

This program intentionally has no publish operation. The Zenodo token only needs
the ``deposit:write`` scope; final publication remains a manual author action.
"""

from __future__ import annotations

import hashlib
import os
import sys
import time
from pathlib import Path
from urllib.parse import quote

import requests


ZENODO_ROOT = "https://zenodo.org"
UPLOAD_DIRECTORY = Path("zenodo_upload_parts")
MANIFEST_PATH = UPLOAD_DIRECTORY / "SHA256SUMS.txt"
SUPPORT_FILES = (
    UPLOAD_DIRECTORY / "ZENODO_UPLOAD_PARTS_README.md",
    MANIFEST_PATH,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_expected_files() -> list[Path]:
    expected: list[Path] = []
    for line in MANIFEST_PATH.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected_hash, filename = line.split(maxsplit=1)
        path = UPLOAD_DIRECTORY / filename.strip()
        if not path.is_file():
            raise RuntimeError(f"Missing upload file: {path}")
        actual_hash = sha256(path)
        if actual_hash != expected_hash.lower():
            raise RuntimeError(f"SHA-256 mismatch for {path.name}")
        expected.append(path)

    for path in SUPPORT_FILES:
        if not path.is_file():
            raise RuntimeError(f"Missing support file: {path}")
        if path not in expected:
            expected.append(path)
    return expected


def request_with_retry(
    session: requests.Session,
    method: str,
    url: str,
    *,
    attempts: int = 4,
    **kwargs: object,
) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.request(method, url, timeout=(30, 600), **kwargs)
            if response.status_code < 500 and response.status_code != 429:
                return response
            last_error = RuntimeError(
                f"Zenodo returned HTTP {response.status_code}: {response.text[:300]}"
            )
        except requests.RequestException as exc:
            last_error = exc
        if attempt < attempts:
            time.sleep(2**attempt)
    raise RuntimeError(f"Request failed after {attempts} attempts: {last_error}")


def upload_legacy(
    session: requests.Session, deposition: dict[str, object], paths: list[Path]
) -> None:
    links = deposition.get("links")
    if not isinstance(links, dict) or not isinstance(links.get("bucket"), str):
        raise RuntimeError("Zenodo draft response did not contain a bucket link")
    bucket_url = str(links["bucket"]).rstrip("/")

    existing = deposition.get("files")
    existing_names = {
        item.get("filename") or item.get("key")
        for item in existing
        if isinstance(existing, list) and isinstance(item, dict)
    } if isinstance(existing, list) else set()

    for path in paths:
        if path.name in existing_names:
            print(f"Already present; skipping: {path.name}")
            continue
        target = f"{bucket_url}/{quote(path.name)}"
        print(f"Uploading {path.name} ({path.stat().st_size} bytes)")
        with path.open("rb") as stream:
            response = request_with_retry(session, "PUT", target, data=stream)
        if response.status_code not in (200, 201):
            raise RuntimeError(
                f"Upload failed for {path.name}: HTTP {response.status_code} "
                f"{response.text[:500]}"
            )


def verify_legacy(
    session: requests.Session, deposition_id: str, expected_paths: list[Path]
) -> None:
    url = f"{ZENODO_ROOT}/api/deposit/depositions/{deposition_id}/files"
    response = request_with_retry(session, "GET", url)
    if response.status_code != 200:
        raise RuntimeError(
            f"Could not verify Zenodo files: HTTP {response.status_code} "
            f"{response.text[:500]}"
        )
    records = response.json()
    remote = {
        item.get("filename"): int(item.get("filesize", -1))
        for item in records
        if isinstance(item, dict)
    }
    errors = [
        path.name
        for path in expected_paths
        if remote.get(path.name) != path.stat().st_size
    ]
    if errors:
        raise RuntimeError("Remote size verification failed: " + ", ".join(errors))
    print(f"Verified {len(expected_paths)} Zenodo files. Draft was not published.")


def main() -> int:
    token = os.environ.get("ZENODO_TOKEN", "").strip()
    deposition_id = os.environ.get("ZENODO_DEPOSITION_ID", "").strip()
    if not token:
        raise RuntimeError("ZENODO_TOKEN secret is missing")
    if not deposition_id.isdigit():
        raise RuntimeError("ZENODO_DEPOSITION_ID must contain digits only")

    paths = load_expected_files()
    print(f"Local integrity checks passed for {len(paths)} files")

    session = requests.Session()
    session.headers.update({"Authorization": f"Bearer {token}"})
    deposition_url = f"{ZENODO_ROOT}/api/deposit/depositions/{deposition_id}"
    response = request_with_retry(session, "GET", deposition_url)
    if response.status_code != 200:
        raise RuntimeError(
            f"Cannot access Zenodo draft {deposition_id}: HTTP "
            f"{response.status_code} {response.text[:500]}"
        )

    deposition = response.json()
    if deposition.get("submitted") is True:
        raise RuntimeError("The Zenodo record is already published; refusing to modify it")

    upload_legacy(session, deposition, paths)
    verify_legacy(session, deposition_id, paths)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # noqa: BLE001 - concise CI failure reporting
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

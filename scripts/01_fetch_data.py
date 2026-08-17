#!/usr/bin/env python3
"""Mirror the shared Xenium Drive folder into ``data/raw/<section>/``.

The folder is world-readable, so no credentials are needed; file listings are
scraped from the Drive web view because the Drive API would require a key.

    python scripts/01_fetch_data.py                # the core set (~150 MB)
    python scripts/01_fetch_data.py --all          # incl. transcripts (~2.2 GB)
    python scripts/01_fetch_data.py --sections M493_2 F536_1
    python scripts/01_fetch_data.py --list-only    # audit what is on Drive

By default only the files the pipeline actually needs are fetched;
``transcripts.parquet`` (~0.5 GB per section) is pulled only for re-segmentation
and subcellular QC, so it is opt-in via ``--all``.

Files already present with a non-zero size are skipped, so the script is safe to
re-run after a partial download.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
# "Xenium Lea" -- full output bundles for the four ROI sections.
ROOT_FOLDER = "11stRgPPOnkQSuFG6LQtrgdD9YSriA5Tz"
FOLDER_MIME = "application/vnd.google-apps.folder"

# Everything the pipeline needs.  transcripts.parquet is ~0.5 GB per section and
# is only used for re-segmentation, so it is excluded unless --all is given.
CORE_FILES = {
    "cell_feature_matrix.h5",
    # Some sections ship the matrix only as the mtx trio, with no .h5.
    "barcodes.tsv.gz",
    "matrix.mtx.gz",
    "features.tsv.gz",
    "cells.parquet",
    "nucleus_boundaries.parquet",
    "cell_boundaries.parquet",
    "experiment.xenium",
    "gene_panel.json",
    "metrics_summary.csv",
}
# Redundant with the files above, or too bulky to mirror by default.
SKIP_ALWAYS = {
    "cell_feature_matrix.zarr.zip",
    "analysis.zarr.zip",
    "analysis_summary.html",
    "cell_boundaries.csv.gz",
}

_ITEM_RE = re.compile(r'\[\[null,"([0-9A-Za-z_-]{20,})"\],null,null,null,"([^"]+)"')
_NAME_RE = re.compile(r'\[\[\["([^"]{1,120})",null,1\]\]\]')
_SIZE_RE = re.compile(r'\[\[\["([\d.]+ (?:KB|MB|GB|bytes))"\]\]\]')


def _fetch_folder_html(folder_id: str, attempts: int = 4) -> str:
    url = f"https://drive.google.com/drive/folders/{folder_id}"
    for attempt in range(attempts):
        proc = subprocess.run(
            ["curl", "-sL", "--max-time", "120", url], capture_output=True, text=True
        )
        if proc.returncode == 0 and len(proc.stdout) > 1000:
            return proc.stdout
        time.sleep(2**attempt)
    raise RuntimeError(f"could not list Drive folder {folder_id}")


def list_folder(folder_id: str) -> list[dict]:
    """Return ``[{id, name, mime, size}]`` for one Drive folder."""
    html = _fetch_folder_html(folder_id)
    items, seen = [], set()
    for match in _ITEM_RE.finditer(html):
        file_id, mime = match.group(1), match.group(2)
        if file_id in seen:
            continue
        window = html[match.end() : match.end() + 1500]
        name = _NAME_RE.search(window)
        if not name:
            continue
        size = _SIZE_RE.search(window)
        seen.add(file_id)
        items.append(
            {
                "id": file_id,
                "name": name.group(1),
                "mime": mime,
                "size": size.group(1) if size else None,
            }
        )
    return items


def walk(folder_id: str, prefix: str = "", depth: int = 0, max_depth: int = 3) -> list[dict]:
    out = []
    for item in list_folder(folder_id):
        path = f"{prefix}/{item['name']}".lstrip("/")
        out.append({**item, "path": path})
        if item["mime"] == FOLDER_MIME and depth < max_depth:
            out.extend(walk(item["id"], path, depth + 1, max_depth))
    return out


# Sections live either at the root or one level down under "additional slides",
# and each may nest files inside cell_feature_matrix/ or morphology_focus/.
# Both are flattened to data/raw/<section>/ so downstream code sees one layout.
SUBDIRS = {"cell_feature_matrix", "morphology_focus"}
CONTAINER_DIRS = {"additional slides"}


def _section_of(path: str) -> str:
    parts = [p for p in path.split("/") if p and p not in CONTAINER_DIRS]
    parts = [p for p in parts if p not in SUBDIRS]
    return parts[-2] if len(parts) >= 2 else "_root"


def _local_path(path: str) -> Path:
    section = _section_of(path)
    name = path.split("/")[-1]
    # Keep the mtx trio together so scanpy can read the directory directly.
    if name in {"barcodes.tsv.gz", "matrix.mtx.gz"} or (
        name == "features.tsv.gz" and "cell_feature_matrix" in path
    ):
        return Path(section) / "cell_feature_matrix" / name
    return Path(section) / name if section != "_root" else Path(name)


def download(file_id: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  skip   {dest.name} (present)")
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  get    {dest.name}")
    proc = subprocess.run(
        [sys.executable, "-m", "gdown", "-q", file_id, "-O", str(dest)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        print(f"  FAILED {dest.name}: {proc.stderr.strip()[:300]}")
        dest.unlink(missing_ok=True)
        return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sections", nargs="*", help="section ids, default: all")
    parser.add_argument("--list-only", action="store_true", help="print the tree, download nothing")
    parser.add_argument("--all", action="store_true", help="also fetch transcripts.parquet (~0.5 GB/section)")
    parser.add_argument("--images", action="store_true", help="also fetch morphology_focus DAPI (~0.7 GB/section)")
    parser.add_argument("--max-depth", type=int, default=4)
    args = parser.parse_args()

    tree = walk(ROOT_FOLDER, max_depth=args.max_depth)
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / "drive_manifest.json").write_text(json.dumps(tree, indent=1))

    if args.list_only:
        for entry in sorted(tree, key=lambda e: e["path"]):
            kind = "dir " if entry["mime"] == FOLDER_MIME else "file"
            print(f"{kind} {entry['path']:60s} {entry['size'] or '-':>10s}")
        return 0

    wanted = set(CORE_FILES)
    if args.all:
        wanted.add("transcripts.parquet")
    if args.images:
        wanted.add("morphology_focus_0000.ome.tif")

    files = [
        e
        for e in tree
        if e["mime"] != FOLDER_MIME
        and e["name"] not in SKIP_ALWAYS
        and (e["name"] in wanted or "/" not in e["path"])
    ]
    failures = []
    for entry in sorted(files, key=lambda e: e["path"]):
        section = _section_of(entry["path"])
        if args.sections and section not in args.sections:
            continue
        print(f"{section}/{entry['name']}")
        if not download(entry["id"], RAW / _local_path(entry["path"])):
            failures.append(entry["path"])

    if failures:
        print(f"\n{len(failures)} download(s) failed:", *failures, sep="\n  ")
        return 1
    print(f"\nMirrored {len(files)} file(s) into {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

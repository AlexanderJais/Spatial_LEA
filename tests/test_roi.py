"""ROI selection must reproduce the cell counts recorded in the ROI export.

These tests need the per-sample ``cells.parquet`` files under ``data/raw`` (see
``scripts/01_fetch_data.py``).  They are skipped when the data is absent.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from spatial_lea.roi import cells_in_roi, load_rois

REPO = Path(__file__).resolve().parents[1]
ROI_FILE = REPO / "config" / "roi_coordinates_male_2vs2.json"
RAW = REPO / "data" / "raw"

# ROI keys use G073_2 while the Drive folder / experiment.xenium use G_073_2.
SAMPLE_DIRS = {"G073_2": "G_073_2"}


def _cells_path(sample: str) -> Path:
    return RAW / SAMPLE_DIRS.get(sample, sample) / "cells.parquet"


ROIS = load_rois(ROI_FILE) if ROI_FILE.exists() else {}


@pytest.mark.parametrize("sample", sorted(ROIS))
def test_roi_reproduces_exported_cell_count(sample: str) -> None:
    path = _cells_path(sample)
    if not path.exists():
        pytest.skip(f"missing {path}; run scripts/01_fetch_data.py")
    roi = ROIS[sample]
    cells = pd.read_parquet(path, columns=["x_centroid", "y_centroid"])
    assert int(cells_in_roi(cells, roi).sum()) == roi.n_cells_expected


def test_rotation_is_required_where_specified() -> None:
    """Ignoring the transform must change the selection, else it is untested."""
    rotated = [s for s, r in ROIS.items() if r.is_rotated]
    assert rotated, "expected at least one rotated ROI in the male 2v2 export"

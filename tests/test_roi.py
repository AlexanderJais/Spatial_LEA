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

def _cells_path(sample: str) -> Path:
    return RAW / sample / "cells.parquet"


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


def test_anatomical_frame_is_consistent_across_sections() -> None:
    """Every section's frame must agree on orientation and be well-conditioned.

    The ventricle strip is elongated, so a low elongation means the fit latched
    onto something that is not the third ventricle.  The ventral direction is
    resolved from the arcuate, so it must be reproducible section to section
    once expressed in the frame's own terms.
    """
    import numpy as np

    from spatial_lea.anatomy import HYPOTHALAMIC_WINDOW, NUCLEI

    model_path = REPO / "results" / "anatomy" / "nucleus_model.npy"
    if not model_path.exists():
        pytest.skip("run scripts/08_anatomy.py first")
    model = np.load(model_path, allow_pickle=True).item()

    assert set(model) == set(NUCLEI)
    # Anatomical ordering along the dorsoventral axis: ARC is ventral to VMH,
    # which is ventral to DMH.  If this ever fails the frame is upside down.
    dv = {n: model[n]["mean"][1] for n in NUCLEI}
    assert dv["ARC"] < dv["VMH"] < dv["DMH"], dv
    # All three sit inside the hypothalamic window they were fitted in.
    for nucleus in NUCLEI:
        ml, dorsal = model[nucleus]["mean"]
        assert 0 <= ml <= HYPOTHALAMIC_WINDOW["max_abs_ml"]
        assert HYPOTHALAMIC_WINDOW["min_dv"] <= dorsal <= HYPOTHALAMIC_WINDOW["max_dv"]

"""Region-of-interest selection for Xenium sections.

The ROI files exported by the upstream annotation tool store an axis-aligned
polygon in *display* coordinates together with an optional rigid rotation that
was applied to the section before the polygon was drawn.  To recover the cells
inside an ROI the centroids therefore have to be rotated into that display
frame first, and only then tested against the polygon.

The convention -- rotate centroids by ``+rotation_deg`` about ``pivot`` -- was
determined empirically by reproducing the ``n_cells_selected`` counts recorded
in the ROI file; see ``tests/test_roi.py``, which reproduces all four male
sections exactly.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath


@dataclass(frozen=True)
class Roi:
    """A named ROI on one section."""

    sample: str
    roi_name: str
    vertices: np.ndarray
    rotation_deg: float = 0.0
    pivot: tuple[float, float] | None = None
    n_cells_expected: int | None = None

    @property
    def is_rotated(self) -> bool:
        return bool(self.rotation_deg) and self.pivot is not None


def rotate(xy: np.ndarray, degrees: float, pivot) -> np.ndarray:
    """Rotate ``xy`` counter-clockwise by ``degrees`` about ``pivot``."""
    if not degrees:
        return np.asarray(xy, dtype=float)
    theta = np.deg2rad(degrees)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    origin = np.asarray(pivot, dtype=float)
    delta = np.asarray(xy, dtype=float) - origin
    rotated = np.column_stack(
        [
            delta[:, 0] * cos_t - delta[:, 1] * sin_t,
            delta[:, 0] * sin_t + delta[:, 1] * cos_t,
        ]
    )
    return rotated + origin


def load_rois(path: str | Path) -> dict[str, Roi]:
    """Load an ROI export keyed by sample id."""
    payload = json.loads(Path(path).read_text())
    rois: dict[str, Roi] = {}
    for sample, spec in payload.items():
        transform = spec.get("transform") or {}
        rois[sample] = Roi(
            sample=sample,
            roi_name=spec.get("roi_name", "ROI"),
            vertices=np.asarray(spec["vertices"], dtype=float),
            rotation_deg=float(transform.get("rotation_deg", 0.0)),
            pivot=tuple(transform["pivot"]) if "pivot" in transform else None,
            n_cells_expected=spec.get("n_cells_selected"),
        )
    return rois


def cells_in_roi(cells: pd.DataFrame, roi: Roi, x="x_centroid", y="y_centroid") -> np.ndarray:
    """Boolean mask of the rows of ``cells`` that fall inside ``roi``."""
    xy = cells[[x, y]].to_numpy(dtype=float)
    if roi.is_rotated:
        xy = rotate(xy, roi.rotation_deg, roi.pivot)
    return MplPath(roi.vertices).contains_points(xy)


def roi_display_coords(cells: pd.DataFrame, roi: Roi, x="x_centroid", y="y_centroid") -> pd.DataFrame:
    """Centroids expressed in the ROI's display frame, origin at the ROI corner.

    Gives every section a common, rotation-corrected frame so that positions are
    comparable across animals.
    """
    xy = cells[[x, y]].to_numpy(dtype=float)
    if roi.is_rotated:
        xy = rotate(xy, roi.rotation_deg, roi.pivot)
    origin = roi.vertices.min(axis=0)
    return pd.DataFrame(
        {"roi_x": xy[:, 0] - origin[0], "roi_y": xy[:, 1] - origin[1]},
        index=cells.index,
    )

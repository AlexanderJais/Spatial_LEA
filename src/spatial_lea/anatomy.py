"""An intrinsic anatomical frame for the mediobasal hypothalamus.

Replaces the hand-drawn ROI rectangle.  A rectangle is a drawing decision: the
four supplied boxes differ in area by 40%, so how much DMH versus LHA each one
contains varies between animals, and any abundance expressed as a fraction of
"cells in the box" inherits that.

Instead every section is given a coordinate system built from its own anatomy:

``dv``
    distance dorsal from the ventral brain surface, along the third-ventricle
    axis, in micrometres.
``ml``
    signed perpendicular distance from the third-ventricle midline.  Nuclei are
    bilateral, so ``abs_ml`` is what the nucleus model is fitted on.

Both landmarks are found from the data, with no manual input:

* the **third ventricle** is traced from its lining -- tanycytes (``Gpr50``) and
  ependymal cells (``Spag16``).  Scattered false positives are removed by taking
  the largest spatial cluster, then the axis is the principal direction of that
  strip.
* the **ventral direction** is fixed by the arcuate nucleus, which lies ventral
  to the ventricle strip: the vector from the strip centroid to the ARC anchor
  centroid points ventrally.  This makes the frame independent of how the
  section happened to be mounted or rotated.

Nuclei are then defined once, in that shared frame, as robust 2-D Gaussians
fitted to marker-defined anchor cells pooled over all sections -- so the same
definition applies to every section, including ones that were never annotated
by hand.

Marker rules below were chosen on precision/recall against the annotated ROI
cells; see ``scripts/07_calibrate_anchors.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import linalg
from sklearn.cluster import DBSCAN

# Anchors favour precision over recall: a landmark is only wrong if the cells it
# is built from are in the wrong place, and recall costs nothing here.
VENTRICLE_RULE = [("Gpr50", 2), ("Spag16", 2)]   # union
ARC_RULE = [("Agrp", 2), ("Pomc", 5)]            # union; 89% / 77% precision
VMH_RULE = [("Adcyap1", 8)]                      # 74% precision, 69% recall
DMH_RULE = [("Grp", 3)]                          # 52% precision -- weakest anchor

# Lining cells further than this from the ARC centroid are not the third
# ventricle.  Generous relative to MBH dimensions, tight enough to exclude the
# lateral ventricles.
SEED_RADIUS_UM = 1800.0

NUCLEI = ("ARC", "VMH", "DMH")

# The MBH is, by anatomical definition, the tissue within roughly 1.5 mm of the
# third ventricle and 1.8 mm of the ventral surface.  These are fixed anatomical
# constants applied identically to every section -- unlike a drawn rectangle,
# they do not vary between animals.  The window matters because Adcyap1 and Grp
# are expressed brain-wide: without it the VMH and DMH anchors are dominated by
# cortex and thalamus rather than hypothalamus.
HYPOTHALAMIC_WINDOW = {"max_abs_ml": 1500.0, "min_dv": -100.0, "max_dv": 1800.0}


def in_window(coords: pd.DataFrame, window: dict | None = None) -> np.ndarray:
    """Cells inside the coarse hypothalamic window of the anatomical frame."""
    w = window or HYPOTHALAMIC_WINDOW
    return (
        (coords["ml"].abs() <= w["max_abs_ml"])
        & (coords["dv"] >= w["min_dv"])
        & (coords["dv"] <= w["max_dv"])
    ).to_numpy()


@dataclass(frozen=True)
class Frame:
    """The anatomical frame of one section."""

    section: str
    origin: np.ndarray      # ventral surface point on the 3V axis
    ventral_dir: np.ndarray # unit vector pointing ventrally
    lateral_dir: np.ndarray # unit vector, perpendicular, right-handed
    n_lining: int
    quality: float          # elongation of the ventricle strip; higher is better


def _marker_mask(counts: pd.DataFrame, rules) -> np.ndarray:
    mask = np.zeros(len(counts), dtype=bool)
    for gene, cut in rules:
        if gene in counts.columns:
            mask |= counts[gene].to_numpy() >= cut
    return mask


def find_frame(xy: np.ndarray, counts: pd.DataFrame, section: str = "") -> Frame:
    """Locate the third-ventricle axis and the ventral direction."""
    lining = _marker_mask(counts, VENTRICLE_RULE)
    if lining.sum() < 50:
        raise ValueError(f"{section}: only {lining.sum()} ventricle-lining cells found")

    arc = _marker_mask(counts, ARC_RULE)
    if arc.sum() < 20:
        raise ValueError(f"{section}: only {arc.sum()} ARC anchor cells found")
    # ARC anchors are hypothalamus-specific, so their centroid seeds the search.
    # Ependymal cells line every ventricle; without this seed the strip fit can
    # latch onto the lateral ventricle instead of the third.
    arc_centre = np.median(xy[arc], axis=0)

    near_arc = np.linalg.norm(xy[lining] - arc_centre, axis=1) <= SEED_RADIUS_UM
    pts = xy[lining][near_arc]
    if len(pts) < 50:
        raise ValueError(f"{section}: only {len(pts)} lining cells near the ARC seed")

    # The lining forms one elongated strip; scattered positives elsewhere are
    # noise.  eps in micrometres, so the value carries physical meaning.
    labels = DBSCAN(eps=60.0, min_samples=10).fit_predict(pts)
    if (labels >= 0).any():
        biggest = pd.Series(labels[labels >= 0]).value_counts().idxmax()
        pts = pts[labels == biggest]

    centre = pts.mean(axis=0)
    centred = pts - centre
    # Total-least-squares line fit: first principal direction of the strip.
    _, sing, vt = linalg.svd(centred, full_matrices=False)
    axis = vt[0]
    quality = float(sing[0] / max(sing[1], 1e-9))

    # ARC sits ventral to the ventricle strip, which fixes the sign of the axis.
    to_arc = arc_centre - centre
    ventral = axis if np.dot(axis, to_arc) > 0 else -axis

    # Origin at the ventral end of the strip (median eminence / tuberal surface),
    # taken as a high percentile rather than the extreme so one stray cell cannot
    # move the whole coordinate system.
    proj = (pts - centre) @ ventral
    origin = centre + ventral * np.percentile(proj, 99)

    lateral = np.array([-ventral[1], ventral[0]])
    return Frame(section, origin, ventral, lateral, int(lining.sum()), quality)


def to_frame(xy: np.ndarray, frame: Frame) -> pd.DataFrame:
    """Express centroids in the anatomical frame (micrometres)."""
    delta = xy - frame.origin
    return pd.DataFrame(
        {
            "dv": -(delta @ frame.ventral_dir),  # positive = dorsal to the surface
            "ml": delta @ frame.lateral_dir,
        }
    )


def _robust_gaussian(points: np.ndarray, trim: float = 0.1):
    """Mean and covariance after trimming the most peripheral ``trim`` fraction."""
    mu = np.median(points, axis=0)
    cov = np.cov(points.T)
    for _ in range(5):
        inv = linalg.pinv(cov)
        d = np.einsum("ij,jk,ik->i", points - mu, inv, points - mu)
        keep = d <= np.quantile(d, 1 - trim)
        mu = points[keep].mean(axis=0)
        cov = np.cov(points[keep].T)
    return mu, cov


def fit_nucleus_model(coords: pd.DataFrame, counts: pd.DataFrame) -> dict:
    """Fit one 2-D Gaussian per nucleus in (abs_ml, dv), pooled over sections.

    Fitted on ``abs_ml`` because ARC, VMH and DMH are bilateral -- folding the
    two sides together doubles the anchor cells and forces the definition to be
    symmetric, which it should be.
    """
    points = np.column_stack([coords["ml"].abs().to_numpy(), coords["dv"].to_numpy()])
    window = in_window(coords)
    model = {}
    for nucleus, rule in (("ARC", ARC_RULE), ("VMH", VMH_RULE), ("DMH", DMH_RULE)):
        mask = _marker_mask(counts, rule) & window
        if mask.sum() < 30:
            raise ValueError(f"{nucleus}: only {mask.sum()} anchor cells to fit")
        mu, cov = _robust_gaussian(points[mask])
        model[nucleus] = {"mean": mu, "cov": cov, "n_anchor": int(mask.sum())}
    return model


def assign_nuclei(coords: pd.DataFrame, model: dict, max_sd: float = 2.0) -> pd.DataFrame:
    """Assign each cell to the nearest nucleus in Mahalanobis distance.

    Cells beyond ``max_sd`` of every nucleus are left unassigned rather than
    forced into the closest one -- the MBH has a boundary, and cells outside it
    should not silently enter a per-nucleus denominator.
    """
    points = np.column_stack([coords["ml"].abs().to_numpy(), coords["dv"].to_numpy()])
    dists = {}
    for nucleus, spec in model.items():
        inv = linalg.pinv(spec["cov"])
        delta = points - spec["mean"]
        dists[nucleus] = np.sqrt(np.einsum("ij,jk,ik->i", delta, inv, delta))
    dist_df = pd.DataFrame(dists, index=coords.index)

    nearest = dist_df.idxmin(axis=1)
    best = dist_df.min(axis=1)
    out = pd.DataFrame(
        {
            "nucleus": np.where(best <= max_sd, nearest, "outside"),
            "nucleus_sd": best,
            # Margin to the runner-up: small values flag boundary-ambiguous cells.
            "nucleus_margin": dist_df.apply(lambda r: np.diff(np.sort(r.to_numpy())[:2])[0], axis=1),
        },
        index=coords.index,
    )
    for nucleus in dist_df.columns:
        out[f"sd_{nucleus}"] = dist_df[nucleus]
    return out


# The annotated cell types that define each nucleus.  Fitting on these rather
# than on raw marker thresholds matters: Agrp/Pomc anchors sit in the ventral
# core of the arcuate, so a Gaussian fitted to them alone is too tight and
# excludes the TIDA and KNDy populations, which are genuinely arcuate.
NUCLEUS_TYPES = {
    "ARC": ["ARC Agrp/Npy", "ARC Pomc", "ARC Th/Slc6a3 (TIDA)", "ARC Tac2/Esr1 (KNDy-like)"],
    "VMH": ["VMH-like Glut Rasgrf2", "VMH-like Glut Calb1", "VMH-like Glut Tac1"],
    "DMH": ["DMH Grp/Ppp1r17"],
}


def fit_nucleus_model_from_labels(
    coords: pd.DataFrame, cell_types: pd.Series, trim: float = 0.05
) -> dict:
    """Fit the nucleus Gaussians from annotated cell types.

    ``coords`` and ``cell_types`` must be aligned.  Only the sections that were
    clustered and annotated are needed -- the fitted model is global and is then
    applied to every section, annotated or not.
    """
    points = np.column_stack([coords["ml"].abs().to_numpy(), coords["dv"].to_numpy()])
    window = in_window(coords)
    labels = cell_types.to_numpy().astype(str)

    model = {}
    for nucleus, types in NUCLEUS_TYPES.items():
        mask = np.isin(labels, types) & window
        if mask.sum() < 30:
            raise ValueError(f"{nucleus}: only {mask.sum()} labelled anchor cells")
        mu, cov = _robust_gaussian(points[mask], trim=trim)
        model[nucleus] = {"mean": mu, "cov": cov, "n_anchor": int(mask.sum())}
    return model

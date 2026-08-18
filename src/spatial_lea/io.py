"""Loading Xenium sections into AnnData, with ROI membership attached.

A Xenium ``cell_feature_matrix.h5`` holds the panel genes *and* the control
channels (negative-control probes, negative-control codewords, genomic controls,
unassigned/deprecated codewords) in one matrix.  The controls are what make it
possible to say whether a low count for a gene like ``Galr1`` is signal or
background, so they are split out into ``obs``/``obsm`` rather than discarded.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from spatial_lea.roi import Roi, cells_in_roi, load_rois, roi_display_coords

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
ROI_FILE = REPO / "config" / "roi_coordinates_male_2vs2.json"

GENE_TYPE = "Gene Expression"
CONTROL_TYPES = (
    "Negative Control Probe",
    "Negative Control Codeword",
    "Genomic Control",
    "Unassigned Codeword",
    "Deprecated Codeword",
)

# section -> animal, and the study design.  Kept here so a loaded object always
# carries its group labels and can never be analysed with them detached.
SECTION_ANIMAL = {
    # Cohort A -- male, panel 7ZBFXR, nucleus-expansion segmentation, 3 sections each.
    "F536_1": "F536", "F536_2": "F536", "F536_3": "F536",
    "G073_1": "G_073", "G073_2": "G_073", "G073_3": "G_073",
    "M399_1": "M399", "M399_2": "M399", "M399_3": "M399",
    "M493_1": "M493", "M493_2": "M493", "M493_3": "M493",
    # Cohort B -- female, panel NCY734, stain-kit segmentation, 1 section each.
    "K238_2": "K238", "P953_1": "P953", "F739_2": "F739", "Q378_2": "Q378",
}
# The four sections with a hand-drawn ROI; kept only to reproduce the original
# analysis and to validate the anatomical frame against it.
ROI_SECTIONS = ("F536_1", "G073_2", "M399_3", "M493_2")
# Flagged by the lab as physically distorted; analysed but reported separately.
SUSPECT_SECTIONS = ("M399_2",)
# The section the lab selected as best for each animal.
SELECTED_SECTIONS = ("F536_1", "G073_2", "M399_3", "M493_2",
                     "K238_2", "P953_1", "F739_2", "Q378_2")

# Eight animals, four aged and four adult, in four blocks.  Sex is completely
# confounded with panel and with segmentation chemistry -- all four males ran on
# 7ZBFXR with nucleus expansion, all four females on NCY734 with the stain kit --
# so no male/female difference can be interpreted.  **Age is not confounded**:
# every block is one aged and one adult animal on the same slide run, panel and
# chemistry, so an age contrast formed within blocks never crosses those
# boundaries.  Sex is a blocking factor here, not a variable.
ANIMAL_META = {
    "F536": {"age_group": "aged", "age_weeks": 70, "batch": "B1", "sex": "male",
             "panel": "7ZBFXR", "detection": "nucleus_expansion"},
    "M493": {"age_group": "adult", "age_weeks": 29, "batch": "B1", "sex": "male",
             "panel": "7ZBFXR", "detection": "nucleus_expansion"},
    "G_073": {"age_group": "aged", "age_weeks": 65, "batch": "B2", "sex": "male",
              "panel": "7ZBFXR", "detection": "nucleus_expansion"},
    "M399": {"age_group": "adult", "age_weeks": 30, "batch": "B2", "sex": "male",
             "panel": "7ZBFXR", "detection": "nucleus_expansion"},
    "K238": {"age_group": "aged", "age_weeks": 73, "batch": "B3", "sex": "female",
             "panel": "NCY734", "detection": "stain_kit"},
    "P953": {"age_group": "adult", "age_weeks": 25, "batch": "B3", "sex": "female",
             "panel": "NCY734", "detection": "stain_kit"},
    "F739": {"age_group": "aged", "age_weeks": 68, "batch": "B4", "sex": "female",
             "panel": "NCY734", "detection": "stain_kit"},
    "Q378": {"age_group": "adult", "age_weeks": 24, "batch": "B4", "sex": "female",
             "panel": "NCY734", "detection": "stain_kit"},
}
# Aged/adult pair per block; the unit every age contrast is formed on.
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399"),
          "B3": ("K238", "P953"), "B4": ("F739", "Q378")}
AGED = tuple(a for a, m in ANIMAL_META.items() if m["age_group"] == "aged")
ADULT = tuple(a for a, m in ANIMAL_META.items() if m["age_group"] == "adult")


def shared_genes(sections=tuple(SECTION_ANIMAL), raw: Path = RAW) -> list[str]:
    """Genes present on every panel in `sections`, in a stable order.

    The two panels overlap in 315 of 325 targets.  Mixing cohorts without
    intersecting them would leave a gene missing from one panel looking like a
    gene absent from those animals, which is a technical zero, not a biological
    one.
    """
    common, order = None, []
    for section in sections:
        genes = _panel_genes(raw / section / "gene_panel.json")
        if common is None:
            common, order = set(genes), list(genes)
        else:
            common &= set(genes)
    return [g for g in order if g in common]


def _panel_genes(path: Path) -> list[str]:
    payload = json.loads(path.read_text()).get("payload", {})
    out = []
    for target in payload.get("targets", []):
        name = target.get("type", {}).get("data", {}).get("name")
        if name:
            out.append(name)
    return out


def _read_matrix(folder: Path) -> ad.AnnData:
    """Read the count matrix, whichever form the bundle shipped it in.

    Most sections carry ``cell_feature_matrix.h5``; a few carry only the mtx
    trio in ``cell_feature_matrix/``.
    """
    h5 = folder / "cell_feature_matrix.h5"
    if h5.exists():
        return sc.read_10x_h5(h5, gex_only=False)
    mtx_dir = folder / "cell_feature_matrix"
    if (mtx_dir / "matrix.mtx.gz").exists():
        return sc.read_10x_mtx(mtx_dir, gex_only=False, make_unique=False)
    raise FileNotFoundError(f"{folder}: no cell_feature_matrix.h5 and no mtx directory")


def load_section(section: str, raw: Path = RAW, roi_file: Path = ROI_FILE) -> ad.AnnData:
    """Load one section: counts, cell metadata, design labels and ROI mask."""
    folder = raw / section
    adata = _read_matrix(folder)
    adata.var_names_make_unique()

    cells = pd.read_parquet(folder / "cells.parquet").set_index("cell_id")
    missing = adata.obs_names.difference(cells.index)
    if len(missing):
        raise ValueError(f"{section}: {len(missing)} cells in matrix absent from cells.parquet")
    adata.obs = cells.reindex(adata.obs_names).copy()

    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy(dtype=float)

    animal = SECTION_ANIMAL[section]
    adata.obs["section"] = pd.Categorical([section] * adata.n_obs)
    adata.obs["animal"] = pd.Categorical([animal] * adata.n_obs)
    for key, value in ANIMAL_META[animal].items():
        adata.obs[key] = pd.Categorical([value] * adata.n_obs)

    # Only four sections were hand-drawn; the rest rely on the anatomical frame
    # (spatial_lea.anatomy), so a missing ROI is normal, not an error.
    rois = load_rois(roi_file)
    roi = rois.get(section)
    if roi is not None:
        adata.obs["in_roi"] = cells_in_roi(adata.obs, roi)
        adata.obs[["roi_x", "roi_y"]] = roi_display_coords(adata.obs, roi)
        adata.uns["roi"] = {
            "name": roi.roi_name,
            "rotation_deg": roi.rotation_deg,
            "n_cells_expected": roi.n_cells_expected,
        }
    else:
        adata.obs["in_roi"] = False

    _split_controls(adata)
    return adata


def _split_controls(adata: ad.AnnData) -> None:
    """Move control channels out of ``X`` into ``obs``/``uns``.

    Per-cell control totals become obs columns (used for background calibration);
    the per-probe control matrix is kept in ``uns`` so the null distribution of a
    single probe can be compared against a single gene.
    """
    feature_type = adata.var["feature_types"].astype(str)
    is_gene = feature_type == GENE_TYPE
    is_control = feature_type.isin(CONTROL_TYPES)

    control = adata[:, is_control]
    counts = np.asarray(control.X.sum(axis=1)).ravel()
    adata.obs["control_counts"] = counts
    for ctype in CONTROL_TYPES:
        sel = feature_type == ctype
        if sel.any():
            key = ctype.lower().replace(" ", "_") + "_counts"
            adata.obs[key] = np.asarray(adata[:, sel].X.sum(axis=1)).ravel()

    # Kept in obsm, not uns, so it is subset along with the cells.  Storing it in
    # uns silently survives a cell filter and then misaligns with X.
    adata.obsm["control_counts_matrix"] = control.X.copy()
    adata.uns["control_names"] = control.var_names.to_numpy()
    adata.uns["control_types"] = feature_type[is_control].to_numpy()

    adata._inplace_subset_var(is_gene.to_numpy())
    adata.obs["gene_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()


def load_all(sections=tuple(SECTION_ANIMAL), raw: Path = RAW, roi_only: bool = False) -> ad.AnnData:
    """Concatenate sections. ``roi_only`` restricts to cells inside the MBH ROI."""
    parts = []
    for section in sections:
        adata = load_section(section, raw=raw)
        if roi_only:
            adata = adata[adata.obs["in_roi"]].copy()
        parts.append(adata)
    merged = ad.concat(parts, label="section_key", keys=list(sections), index_unique="-", merge="same")
    merged.uns["sections"] = list(sections)
    return merged


def counts_vector(adata, gene: str) -> np.ndarray:
    """Raw counts for one gene as a dense 1-D array.

    ``adata.X`` is z-scored after ``sc.pp.scale``, so ``X > 0`` is meaningless as
    a detection test.  Counts always come from the ``counts`` layer.
    """
    if "counts" not in adata.layers:
        raise KeyError("no 'counts' layer; load with load_section or 03_cluster.py")
    col = adata[:, gene].layers["counts"]
    return np.asarray(col.todense() if hasattr(col, "todense") else col).ravel()


def counts_matrix(adata) -> np.ndarray:
    mat = adata.layers["counts"]
    return np.asarray(mat.todense() if hasattr(mat, "todense") else mat)

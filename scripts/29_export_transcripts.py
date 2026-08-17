#!/usr/bin/env python3
"""Export transcripts for re-segmentation, cropped to the hypothalamic window.

Why re-segment at all.  Every result in this repository rests on the vendor
segmentation: a DAPI nucleus dilated by a fixed 5 um.  Measuring the objects it
produces showed the dilation is not neutral --

  * neuronal nuclei are 8.4 um across, glial nuclei 6.5 um
  * the realised expansion is 3.8 um, not 5, because dilation stops at neighbours
  * so a "neuron" ends up 15.8 um wide against a real soma of 15-25 um
    (under-captured), and a "glial cell" 13.9 um against a real 8-12 um
    (over-expanded, ~78% of its area is dilation territory)

Under-capturing neurons loses their cytoplasmic transcripts; over-expanding glia
imports the neighbours'.  Both push measurements the same way -- toward every
cell looking like the local average -- which is the most likely explanation for
the ~40% of cells carrying both GABAergic and glutamatergic markers.

Baysor does not assume a shape.  It assigns each molecule to a cell using the
local molecular composition as well as position, so a neuron can claim the
transcripts in its processes and a glial cell is not forced to swallow a fixed
disc of neuropil.

What is exported
----------------
One CSV per section, in the section's own micrometre coordinates, containing
only molecules inside the anatomical hypothalamic window that every downstream
analysis already uses (|ml| <= 1500 um, dv in [-100, 1800] um).  The window is
derived from the section's own third ventricle and ventral surface, so the same
tissue is cut out of each animal -- this is a crop, not a subsample, and it is
the same crop the existing results are computed on.

The prior column holds the Xenium cell id for molecules **inside a nucleus**,
and 0 elsewhere.  That hands Baysor the DAPI segmentation -- which is real
measurement -- while leaving the cytoplasmic extent, which is the part we do not
trust, for it to infer.

    python scripts/29_export_transcripts.py G073_1 M399_3
    python scripts/29_export_transcripts.py G073_1 --patch 600   # timing pilot
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.anatomy import HYPOTHALAMIC_WINDOW, find_frame, to_frame  # noqa: E402
from spatial_lea.io import RAW, counts_matrix, load_section  # noqa: E402

OUT = REPO / "data" / "segmentation"
MIN_COUNTS, MIN_GENES = 10, 5
ANCHOR_GENES = ["Gpr50", "Spag16", "Agrp", "Pomc", "Adcyap1", "Grp", "Ppp1r17"]
MIN_QV = 20.0                 # 10x's own call threshold; drops 4.4% of molecules
BATCH = 2_000_000
COLS = ["cell_id", "overlaps_nucleus", "feature_name",
        "x_location", "y_location", "z_location", "qv", "is_gene"]


def section_frame(section: str):
    """The section's anatomical frame, from its own ventricle and ARC."""
    adata = load_section(section)
    keep = ((adata.obs["gene_counts"] >= MIN_COUNTS)
            & (adata.obs["n_genes"] >= MIN_GENES)
            & (adata.obs["nucleus_count"] == 1))
    adata = adata[keep.to_numpy()].copy()
    adata.layers["counts"] = adata.X.copy()   # load_section returns raw counts in X
    present = [g for g in ANCHOR_GENES if g in adata.var_names]
    counts = pd.DataFrame(counts_matrix(adata[:, present]), columns=present,
                          index=adata.obs_names)
    xy = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    return find_frame(xy, counts, section), adata


def export(section: str, patch: float | None) -> dict:
    frame, adata = section_frame(section)
    src = RAW / section / "transcripts.parquet"
    pf = pq.ParquetFile(src)

    # Centre of the window in world coordinates, used only for --patch.
    w = HYPOTHALAMIC_WINDOW
    centre_dv = (w["min_dv"] + w["max_dv"]) / 2
    centre = frame.origin - centre_dv * frame.ventral_dir

    kept, n_total, n_gene = [], 0, 0
    for batch in pf.iter_batches(batch_size=BATCH, columns=COLS):
        d = batch.to_pandas()
        n_total += len(d)
        d = d[d["is_gene"].to_numpy() & (d["qv"].to_numpy() >= MIN_QV)]
        n_gene += len(d)
        if not len(d):
            continue
        xy = d[["x_location", "y_location"]].to_numpy()
        f = to_frame(xy, frame)
        inside = ((f["ml"].abs() <= w["max_abs_ml"])
                  & (f["dv"] >= w["min_dv"]) & (f["dv"] <= w["max_dv"])).to_numpy()
        if patch:
            inside &= (np.abs(xy - centre) <= patch / 2).all(axis=1)
        if not inside.any():
            continue
        d = d[inside].copy()
        d["ml"] = f["ml"].to_numpy()[inside]
        d["dv"] = f["dv"].to_numpy()[inside]
        kept.append(d)

    mol = pd.concat(kept, ignore_index=True)
    # Xenium writes feature_name and cell_id as bytes in some output versions.
    for col in ("feature_name", "cell_id"):
        if isinstance(mol[col].iloc[0], bytes):
            mol[col] = mol[col].str.decode("utf-8")

    # Prior: the DAPI nucleus each molecule falls inside, 0 for the rest.
    # Baysor requires integers, so the string ids are factorised; the mapping is
    # written alongside so a Baysor cell can be traced back to a Xenium cell.
    in_nuc = (mol["overlaps_nucleus"].to_numpy() == 1) & \
             (mol["cell_id"].to_numpy() != "UNASSIGNED")
    codes, uniques = pd.factorize(mol.loc[in_nuc, "cell_id"])
    mol["prior"] = 0
    mol.loc[in_nuc, "prior"] = codes + 1

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / f"{section}_molecules.csv"
    mol[["x_location", "y_location", "z_location", "feature_name", "prior",
         "ml", "dv", "overlaps_nucleus"]].rename(
        columns={"x_location": "x", "y_location": "y", "z_location": "z",
                 "feature_name": "gene"}).to_csv(dest, index=False)
    pd.Series(uniques, name="xenium_cell_id").to_csv(
        OUT / f"{section}_prior_ids.csv", index_label="prior")

    # Cells whose centroid is in the same window, for the like-for-like count.
    cf = to_frame(adata.obs[["x_centroid", "y_centroid"]].to_numpy(), frame)
    n_cells = int(((cf["ml"].abs() <= w["max_abs_ml"])
                   & (cf["dv"] >= w["min_dv"])
                   & (cf["dv"] <= w["max_dv"])).sum())

    return {
        "section": section, "molecules_total": n_total, "molecules_qc": n_gene,
        "molecules_window": len(mol), "pct_of_section": round(len(mol) / n_gene * 100, 1),
        "in_nucleus": int(in_nuc.sum()),
        "pct_in_nucleus": round(in_nuc.mean() * 100, 1),
        "prior_cells": len(uniques), "xenium_cells_in_window": n_cells,
        "genes": mol["feature_name"].nunique(),
        "median_per_prior_cell": int(np.median(np.bincount(mol.loc[in_nuc, "prior"])[1:])),
        "area_mm2": round((2 * w["max_abs_ml"]) * (w["max_dv"] - w["min_dv"]) / 1e6, 2),
        "file_mb": round(dest.stat().st_size / 1e6, 1),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sections", nargs="+")
    p.add_argument("--patch", type=float, default=None,
                   help="crop further to a square of this side (um) for timing")
    args = p.parse_args()

    rows = []
    for section in args.sections:
        print(f"--- {section} ---", flush=True)
        rows.append(export(section, args.patch))
        print(pd.Series(rows[-1]).to_string(), flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "export_summary.csv", index=False)
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Cluster the anatomically-defined MBH across all 12 sections and label it.

Cell types are not re-annotated from scratch.  The four hand-drawn sections were
already annotated against marker evidence (04_annotate.py), so those labels are
carried onto the full set: cluster all 12 sections together, then give each
cluster the majority label of the annotated cells inside it.  That keeps the
naming reviewable and stops the label set drifting between runs.

Clusters containing too few annotated cells to call are reported rather than
silently named.
"""

from __future__ import annotations

import sys
from pathlib import Path

import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "mbh_atlas"
RESOLUTION = 1.5          # finer than the ROI run: more cells, more sections
MIN_ANNOTATED = 15        # annotated cells needed to name a cluster
SEED = 0


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 4

    full = sc.read_h5ad(PROC / "mbh_anatomical.h5ad")
    adata = full[full.obs["in_mbh"].to_numpy()].copy()
    del full
    print(f"MBH cells across {adata.obs['section'].nunique()} sections: {adata.n_obs:,}")
    print(pd.crosstab(adata.obs["animal"], adata.obs["nucleus"]).to_string(), "\n")

    adata.layers["counts"] = adata.X.copy()
    target = float(np.median(np.asarray(adata.layers["counts"].sum(axis=1)).ravel()))
    sc.pp.normalize_total(adata, target_sum=target)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)

    harmony = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["animal"], random_state=SEED)
    corrected = np.asarray(harmony.Z_corr)
    adata.obsm["X_pca_harmony"] = corrected if corrected.shape[0] == adata.n_obs else corrected.T

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
    sc.tl.leiden(adata, resolution=RESOLUTION, key_added="leiden",
                 flavor="igraph", n_iterations=2, random_state=SEED)
    print(f"{adata.obs['leiden'].nunique()} clusters at resolution {RESOLUTION}")

    labels = _annotated_labels(adata)
    adata.obs["ref_label"] = labels.to_numpy()

    # Majority vote per cluster, with the vote's own purity carried alongside so
    # a weakly-supported name is visible rather than implied to be solid.
    rows = []
    mapping = {}
    for cluster, sub in adata.obs.groupby("leiden", observed=True):
        known = sub["ref_label"].dropna()
        if len(known) < MIN_ANNOTATED:
            mapping[cluster] = "Unresolved"
            rows.append({"cluster": cluster, "n_cells": len(sub), "n_annotated": len(known),
                         "label": "Unresolved", "purity_pct": np.nan})
            continue
        counts = known.value_counts()
        mapping[cluster] = counts.index[0]
        rows.append({"cluster": cluster, "n_cells": len(sub), "n_annotated": len(known),
                     "label": counts.index[0],
                     "purity_pct": round(counts.iloc[0] / len(known) * 100, 1)})

    adata.obs["cell_type"] = pd.Categorical(adata.obs["leiden"].map(mapping))
    transfer = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    transfer.to_csv(OUT / "label_transfer.csv", index=False)

    print("\n=== Label transfer (majority annotated label per cluster) ===")
    with pd.option_context("display.width", 200):
        print(transfer.to_string(index=False))
    weak = transfer[(transfer.purity_pct < 60) | transfer.purity_pct.isna()]
    print(f"\n{len(weak)} of {len(transfer)} clusters have purity < 60% or too few "
          f"annotated cells; they carry {int(weak.n_cells.sum()):,} cells.")

    print("\n=== Cell types by nucleus (all 12 sections) ===")
    tab = pd.crosstab(adata.obs["cell_type"], adata.obs["nucleus"])
    tab["total"] = tab.sum(axis=1)
    print(tab.sort_values("total", ascending=False).to_string())

    adata.write_h5ad(PROC / "mbh_atlas_12.h5ad")
    tab.to_csv(OUT / "celltype_by_nucleus.csv")
    print(f"\nWrote {PROC / 'mbh_atlas_12.h5ad'}")
    return 0


def _annotated_labels(adata) -> pd.Series:
    ref = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    key = ref.obs["section"].astype(str) + "|" + ref.obs_names.str.split("-").str[0]
    ref_map = pd.Series(ref.obs["cell_type"].astype(str).to_numpy(), index=key)
    ref_map = ref_map[~ref_map.index.duplicated()]
    own = adata.obs["section"].astype(str) + "|" + adata.obs_names.str.split("-").str[0]
    out = pd.Series(ref_map.reindex(own).to_numpy(), index=adata.obs_names)
    print(f"{out.notna().sum():,} of {len(out):,} MBH cells carry a prior annotation")
    return out


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Cluster the MBH ROI cells and dump per-cluster markers for annotation.

Deliberately exploratory: it produces the clustering and the evidence needed to
name the clusters, but does not name them.  Annotation happens in 04_annotate.py
against a written marker table, so the naming step is reviewable rather than
buried in a heuristic.

With a 297-gene panel every gene is informative, so no HVG selection is done.
Integration is over animal, and an un-integrated embedding is kept alongside so
the effect of integration on the age signal can be inspected rather than assumed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import harmonypy
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import SECTION_ANIMAL, load_section  # noqa: E402

OUT = REPO / "results" / "cluster"
PROC = REPO / "data" / "processed"
MIN_COUNTS, MIN_GENES = 10, 5
RESOLUTIONS = (0.5, 1.0, 1.5, 2.0)
SEED = 0


def build() -> sc.AnnData:
    parts = []
    for section in SECTION_ANIMAL:
        adata = load_section(section)
        roi = adata[adata.obs["in_roi"]].copy()
        keep = (
            (roi.obs["gene_counts"] >= MIN_COUNTS)
            & (roi.obs["n_genes"] >= MIN_GENES)
            & (roi.obs["nucleus_count"] == 1)
        )
        parts.append(roi[keep.to_numpy()].copy())
    merged = sc.concat(parts, label="section", keys=list(SECTION_ANIMAL), index_unique="-", merge="same")
    merged.obs["section"] = merged.obs["section"].astype("category")
    return merged


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 4

    adata = build()
    print(f"MBH ROI cells after QC: {adata.n_obs:,} across {adata.obs['animal'].nunique()} animals")
    print(adata.obs.groupby(["animal", "age_group"], observed=True).size().to_string())

    adata.layers["counts"] = adata.X.copy()
    # Median-count target keeps the scale near the data instead of an arbitrary 1e4.
    target = float(np.median(np.asarray(adata.layers["counts"].sum(axis=1)).ravel()))
    sc.pp.normalize_total(adata, target_sum=target)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)
    adata.obsm["X_pca_unintegrated"] = adata.obsm["X_pca"].copy()

    # Called directly rather than via sc.external.pp.harmony_integrate: harmonypy
    # 2.0 returns Z_corr already as (cells, PCs), and the scanpy wrapper still
    # transposes it, which raises a shape error.
    harmony = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["animal"], random_state=SEED)
    corrected = np.asarray(harmony.Z_corr)
    if corrected.shape[0] != adata.n_obs:
        corrected = corrected.T
    adata.obsm["X_pca_harmony"] = corrected

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    for res in RESOLUTIONS:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2, random_state=SEED)
        print(f"  leiden res={res}: {adata.obs[key].nunique()} clusters")

    # Integration diagnostic: how mixed are animals within clusters, and is the
    # age split still visible?  Reported, not assumed away.
    key = "leiden_1.0"
    mix = (
        pd.crosstab(adata.obs[key], adata.obs["animal"], normalize="index")
        .pipe(lambda d: -(d * np.log(d + 1e-12)).sum(axis=1) / np.log(d.shape[1]))
        .rename("animal_entropy")
    )
    mix.to_csv(OUT / "cluster_animal_entropy.csv")
    print(f"\nAnimal-mixing entropy per cluster (1.0 = perfectly mixed): "
          f"median {mix.median():.3f}, min {mix.min():.3f}")

    sc.tl.rank_genes_groups(adata, key, method="wilcoxon", n_genes=25)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(OUT / "cluster_markers.csv", index=False)

    sizes = adata.obs[key].value_counts().sort_index()
    comp = pd.crosstab(adata.obs[key], adata.obs["animal"])
    summary = pd.DataFrame({
        "n_cells": sizes,
        "pct": (sizes / sizes.sum() * 100).round(2),
        "animal_entropy": mix.round(3),
        "top_markers": markers.groupby("group", observed=True)
                              .head(8).groupby("group", observed=True)["names"]
                              .apply(lambda s: ", ".join(s)),
    }).join(comp)
    summary.to_csv(OUT / "cluster_summary.csv")

    print(f"\n=== {adata.obs[key].nunique()} clusters at resolution 1.0 ===")
    with pd.option_context("display.width", 250, "display.max_colwidth", 95):
        print(summary[["n_cells", "pct", "animal_entropy", "top_markers"]].to_string())

    adata.write_h5ad(PROC / "mbh_roi_clustered.h5ad")
    print(f"\nWrote {PROC / 'mbh_roi_clustered.h5ad'} and tables to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

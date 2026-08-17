#!/usr/bin/env python3
"""Sweep every Gal- or Galr1-expressing population for an age difference.

No modelling.  The criterion is the one that survived audit: **complete
separation of the six aged sections from the six adult sections**, with the
direction consistent inside both blocks.  It needs no AP correction, no
covariate and no distributional assumption -- either the ranges overlap or they
do not.

Populations searched:
  * nucleus x annotated cell type
  * subclusters of the Gal-positive cells, on their own
  * subclusters of the Galr1-positive cells, on their own

Calibration: the identical criterion is applied to all 297 panel genes in each
population, so every hit comes with the number of panel genes that separate just
as cleanly there.  A population where 40 genes separate is telling you about the
population, not about galanin.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "gal_sweep"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
GENES = ["Gal", "Galr1"]
MIN_PER_SECTION = 15
MIN_SECTIONS = 12          # all twelve, so separation is over the full set
SEED = 0


def subcluster(adata, gene: str, resolution=1.0) -> pd.Series:
    """Leiden subclusters of the cells positive for one gene."""
    gi = int(np.where(adata.var_names.to_numpy() == gene)[0][0])
    pos = counts_matrix(adata[:, [gi]]).ravel() > 0
    # Build from the counts layer.  adata.X in this object is already
    # normalised and scaled, so re-normalising it produces a degenerate matrix.
    sub = sc.AnnData(counts_matrix(adata[pos]).astype(np.float32))
    sub.obs_names = adata.obs_names[pos]
    sub.var_names = adata.var_names
    # Genes with no counts inside the subset have zero variance, and scaling
    # them yields NaN, which PCA rejects.  Drop them for the embedding only.
    sub = sub[:, np.asarray(sub.X.sum(axis=0)).ravel() > 0].copy()
    sc.pp.normalize_total(sub)
    sc.pp.log1p(sub)
    sc.pp.scale(sub, max_value=10)
    sub.X = np.nan_to_num(np.asarray(sub.X))
    sc.tl.pca(sub, n_comps=min(30, sub.n_vars - 1), svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(sub, n_neighbors=15, random_state=SEED)
    sc.tl.leiden(sub, resolution=resolution, key_added="sub",
                 flavor="igraph", n_iterations=2, random_state=SEED)
    labels = pd.Series(f"{gene}+ sub", index=adata.obs_names, dtype=object)
    labels[:] = np.nan
    labels[sub.obs_names] = [f"{gene}+ c{c}" for c in sub.obs["sub"]]
    print(f"  {gene}+ cells: {int(pos.sum()):,} -> {sub.obs['sub'].nunique()} subclusters")
    return labels


def separation_screen(counts: np.ndarray, sections: np.ndarray, genes: np.ndarray,
                      mask: np.ndarray) -> pd.DataFrame:
    """For every gene in one population: does it separate aged from adult sections?"""
    sec = sections[mask]
    keep = [s for s in pd.unique(sec) if (sec == s).sum() >= MIN_PER_SECTION]
    if len(keep) < MIN_SECTIONS:
        return pd.DataFrame()
    sub = counts[mask]
    mat = np.vstack([sub[sec == s].sum(axis=0) for s in keep])
    tot = mat.sum(axis=1, keepdims=True)
    cpm = np.log2(mat / tot * 1e6 + 1)

    grp = np.array([ANIMAL_META[SECTION_ANIMAL[s]]["age_group"] for s in keep])
    animals = np.array([SECTION_ANIMAL[s] for s in keep])
    aged, adult = cpm[grp == "aged"], cpm[grp == "adult"]
    sep_up = aged.min(axis=0) > adult.max(axis=0)
    sep_dn = aged.max(axis=0) < adult.min(axis=0)

    # Direction must also hold inside each block, so one deviant animal cannot
    # carry the result on its own.
    block_ok = np.ones(cpm.shape[1], dtype=bool)
    for a, d in BLOCKS.values():
        va, vd = cpm[animals == a], cpm[animals == d]
        if not len(va) or not len(vd):
            return pd.DataFrame()
        up = va.min(axis=0) > vd.max(axis=0)
        dn = va.max(axis=0) < vd.min(axis=0)
        block_ok &= (up & sep_up) | (dn & sep_dn)

    lfc = np.mean([cpm[animals == a].mean(axis=0) - cpm[animals == d].mean(axis=0)
                   for a, d in BLOCKS.values()], axis=0)
    return pd.DataFrame({"gene": genes, "separates": (sep_up | sep_dn) & block_ok,
                         "lfc": lfc, "expressed": mat.sum(axis=0) >= 200})


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 4
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(
        ["ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH"])].copy()

    print("Subclustering the ligand- and receptor-positive cells")
    schemes = {"nucleus x cell type":
               (adata.obs["nucleus_ext"].astype(str) + " | "
                + adata.obs["cell_type"].astype(str)).astype(object)}
    for gene in GENES:
        schemes[f"{gene}+ subclusters"] = subcluster(adata, gene)

    counts = counts_matrix(adata)
    sections = adata.obs["section"].astype(str).to_numpy()
    genes = adata.var_names.to_numpy()

    rows = []
    for scheme, labels in schemes.items():
        labels = pd.Series(labels).reindex(adata.obs_names)
        for pop in labels.dropna().unique():
            mask = (labels == pop).to_numpy()
            if mask.sum() < MIN_PER_SECTION * MIN_SECTIONS:
                continue
            res = separation_screen(counts, sections, genes, mask)
            if res.empty:
                continue
            n_panel = int((res["separates"] & res["expressed"]).sum())
            n_tested = int(res["expressed"].sum())
            for gene in GENES:
                r = res[res.gene == gene]
                if r.empty or not bool(r["expressed"].iloc[0]):
                    continue
                rows.append({
                    "scheme": scheme, "population": pop, "gene": gene,
                    "n_cells": int(mask.sum()),
                    "separates": bool(r["separates"].iloc[0]),
                    "lfc": round(float(r["lfc"].iloc[0]), 2),
                    "panel_genes_separating": n_panel,
                    "panel_genes_tested": n_tested,
                    "panel_sep_rate_pct": round(n_panel / max(n_tested, 1) * 100, 1),
                })

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "separation_sweep.csv", index=False)
    print(f"\nScreened {res['population'].nunique()} populations "
          f"across {res['scheme'].nunique()} schemes\n")

    hits = res[res["separates"]].sort_values("panel_sep_rate_pct")
    print("=== Populations where Gal or Galr1 fully separates aged from adult ===")
    print("   (all 6 aged sections one side of all 6 adult, consistent in both blocks)")
    if len(hits):
        with pd.option_context("display.width", 240, "display.max_colwidth", 34):
            print(hits[["scheme", "population", "gene", "n_cells", "lfc",
                        "panel_genes_separating", "panel_genes_tested",
                        "panel_sep_rate_pct"]].to_string(index=False))
        print("\n  panel_sep_rate_pct is the false-positive rate for that population:")
        print("  a hit where 30% of the panel also separates is not about galanin.")
    else:
        print("  none")

    clean = hits[hits["panel_sep_rate_pct"] < 10]
    print(f"\n=== Hits in populations where the panel rate is under 10%: {len(clean)} ===")
    if len(clean):
        print(clean[["population", "gene", "n_cells", "lfc",
                     "panel_sep_rate_pct"]].to_string(index=False))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

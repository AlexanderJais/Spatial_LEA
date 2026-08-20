#!/usr/bin/env python3
"""Delineate hypothalamic nuclei from tissue structure, not from one marker.

The Gaussian model in 08_anatomy.py fits each nucleus to marker-defined anchor
cells.  That works where several cell types anchor a nucleus (ARC: 1856 cells
across 4 types; VMH: 5525 across 3) but not for the DMH, which had a single
anchor type -- 343 Grp/Ppp1r17 cells.  Grp marks a restricted DMH subpopulation,
so the fitted ellipse is the extent of the Grp+ core, not of the nucleus, and it
comes out ~477 um dorsoventrally against a published 500-700.

This takes the standard alternative: describe every cell by the cell-type
composition of its spatial neighbourhood, then cluster those descriptions.
Domains emerge with the shape the tissue actually has instead of an ellipse, and
a nucleus anchored by many minor types is delineated as well as one anchored by
a single abundant marker.

Runs over the whole hypothalamic window, not just the previously-assigned MBH,
so a domain is free to extend past where the ellipses stopped.
"""

from __future__ import annotations

import sys
from pathlib import Path

import harmonypy
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.anatomy import in_window  # noqa: E402
from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "domains"
N_NEIGHBOURS = 30       # spatial neighbours defining a niche, ~90 um radius here
CELL_RES = 1.5          # Leiden resolution for cell types
# Niches are clustered with KMeans at a fixed k rather than Leiden.  Leiden at
# any resolution that resolved the DMH also split the meningeal rim into dozens
# of 40-cell fragments; a fixed k keeps domains at nucleus scale, which is the
# scale the question is asked at.  K_DOMAINS is swept and reported below.
K_DOMAINS = 14
K_SWEEP = (8, 10, 12, 14, 16, 20)
SEED = 0


def build_window(adata):
    coords = pd.DataFrame({"ml": adata.obs["ml"].to_numpy(), "dv": adata.obs["dv"].to_numpy()})
    keep = in_window(coords)
    sub = adata[keep].copy()
    print(f"Hypothalamic window: {sub.n_obs:,} cells across {sub.obs['section'].nunique()} sections")
    return sub


def cluster_cells(adata):
    adata.layers["counts"] = adata.X.copy()
    target = float(np.median(np.asarray(adata.layers["counts"].sum(axis=1)).ravel()))
    sc.pp.normalize_total(adata, target_sum=target)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)
    harmony = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["animal"], random_state=SEED)
    z = np.asarray(harmony.Z_corr)
    adata.obsm["X_pca_harmony"] = z if z.shape[0] == adata.n_obs else z.T
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
    sc.tl.leiden(adata, resolution=CELL_RES, key_added="leiden",
                 flavor="igraph", n_iterations=2, random_state=SEED)
    print(f"  {adata.obs['leiden'].nunique()} cell clusters")
    return adata


def transfer_labels(adata) -> pd.Series:
    """Name each cluster from the annotated reference cells inside it.

    The reference is the hand-drawn ROI annotation, which exists only for the
    male sections.  Clusters that contain no reference cell -- which the female
    cohort can produce, since it ran a different panel and segmentation
    chemistry -- are named instead by correlating their mean profile against the
    reference cell-type profiles.  The two routes are counted separately so it
    stays visible how many labels rest on which.
    """
    ref = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    key = ref.obs["section"].astype(str) + "|" + ref.obs_names.str.split("-").str[0]
    ref_map = pd.Series(ref.obs["cell_type"].astype(str).to_numpy(), index=key)
    ref_map = ref_map[~ref_map.index.duplicated()]
    own = adata.obs["section"].astype(str) + "|" + adata.obs_names.str.split("-").str[0]
    prior = pd.Series(ref_map.reindex(own).to_numpy(), index=adata.obs_names)

    clusters = adata.obs["leiden"].astype(str)
    mapping, by_prior, by_profile = {}, 0, 0
    for cluster, sub in prior.groupby(clusters, observed=True):
        known = sub.dropna()
        if len(known) >= 15:
            mapping[cluster] = known.value_counts().index[0]
            by_prior += 1
        else:
            mapping[cluster] = None

    unnamed = [c for c, v in mapping.items() if v is None]
    if unnamed:
        profiles = _reference_profiles(ref, adata.var_names)
        expr = _cluster_profiles(adata, clusters, unnamed, list(profiles.columns))
        for cluster, vec in expr.items():
            corr = profiles.apply(lambda p: np.corrcoef(p, vec)[0, 1], axis=1)
            mapping[cluster] = (corr.idxmax() if corr.max() > 0.3
                                else f"unlabelled_{cluster}")
            by_profile += 1

    labelled = clusters.map(mapping)
    print(f"  {prior.notna().sum():,} cells carried a reference label; "
          f"{by_prior} clusters named from reference cells, "
          f"{by_profile} from profile correlation")
    return pd.Series(labelled.to_numpy(), index=adata.obs_names)


def _reference_profiles(ref, var_names) -> pd.DataFrame:
    """Mean centred log-CPM profile per annotated cell type, on the shared genes."""
    genes = [g for g in var_names if g in ref.var_names]
    counts = ref[:, genes].layers["counts"]
    counts = np.asarray(counts.todense() if hasattr(counts, "todense") else counts)
    cpm = np.log1p(counts / np.maximum(counts.sum(axis=1, keepdims=True), 1) * 1e4)
    df = pd.DataFrame(cpm, columns=genes)
    df["cell_type"] = ref.obs["cell_type"].astype(str).to_numpy()
    df = df[~df["cell_type"].str.startswith("unlabelled")]
    prof = df.groupby("cell_type").mean()
    return prof.sub(prof.mean(axis=1), axis=0)


def _cluster_profiles(adata, clusters, wanted, genes) -> dict:
    """Mean centred log-CPM per cluster, in the same gene order as the profiles."""
    counts = adata[:, genes].layers["counts"]
    counts = np.asarray(counts.todense() if hasattr(counts, "todense") else counts)
    out = {}
    for cluster in wanted:
        sel = (clusters == cluster).to_numpy()
        block = counts[sel]
        v = np.log1p(block.sum(axis=0) / max(block.sum(), 1) * 1e4)
        out[cluster] = v - v.mean()
    return out


def niche_composition(adata, labels: pd.Series) -> pd.DataFrame:
    """Cell-type fractions among each cell's spatial neighbours, within section."""
    types = pd.Categorical(labels)
    onehot = pd.get_dummies(types).to_numpy(dtype=float)
    out = np.zeros_like(onehot)
    for section in adata.obs["section"].unique():
        sel = (adata.obs["section"] == section).to_numpy()
        xy = adata.obs.loc[sel, ["ml", "dv"]].to_numpy()
        nn = NearestNeighbors(n_neighbors=min(N_NEIGHBOURS, sel.sum())).fit(xy)
        _, idx = nn.kneighbors(xy)
        block = onehot[sel]
        out[sel] = block[idx].mean(axis=1)
    return pd.DataFrame(out, index=adata.obs_names, columns=types.categories)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 4

    cache = PROC / "hypothalamus_window_typed.h5ad"
    niche_path = OUT / "niche_composition.parquet"
    if cache.exists() and niche_path.exists():
        print(f"Reusing {cache.name}")
        adata = sc.read_h5ad(cache)
        niche = pd.read_parquet(niche_path)
    else:
        full = sc.read_h5ad(PROC / "mbh_anatomical.h5ad")
        adata = build_window(full)
        del full
        adata = cluster_cells(adata)
        adata.obs["cell_type"] = transfer_labels(adata).to_numpy()
        print("\nBuilding spatial niches...")
        niche = niche_composition(adata, adata.obs["cell_type"])
        niche.to_parquet(niche_path)
        adata.write_h5ad(cache)

    from sklearn.decomposition import PCA
    # A niche is described by the cell-type composition of its neighbourhood, so
    # the feature count is however many cell types this run resolved -- not a
    # constant.  Asking for 30 components assumed at least 30 types exist and
    # raised outright when 27 were resolved, which is a legitimate outcome of
    # 04_annotate dropping clusters whose markers did not support any label.
    n_comp = min(30, niche.shape[1])
    if n_comp < 30:
        print(f"  {niche.shape[1]} cell types in the niche descriptors: "
              f"using {n_comp} components rather than 30")
    pcs = PCA(n_components=n_comp, random_state=SEED).fit_transform(niche.to_numpy())

    print("\n=== Choosing the number of domains ===")
    print("   smallest domain should stay at nucleus scale, not fragment the meninges")
    for k in K_SWEEP:
        lab = KMeans(n_clusters=k, n_init=10, random_state=SEED).fit_predict(pcs)
        sizes = pd.Series(lab).value_counts()
        print(f"   k={k:3d}  smallest {sizes.min():6,d} cells  "
              f"largest {sizes.max():6,d}  median {int(sizes.median()):6,d}")

    adata.obs["domain"] = pd.Categorical(
        KMeans(n_clusters=K_DOMAINS, n_init=20, random_state=SEED).fit_predict(pcs).astype(str)
    )
    print(f"\nUsing k = {K_DOMAINS}: {adata.obs['domain'].nunique()} spatial domains")

    # Name each domain by where it sits and what defines it.
    rows = []
    for dom, sub in adata.obs.groupby("domain", observed=True):
        comp = sub["cell_type"].value_counts(normalize=True)
        rows.append({
            "domain": dom, "n_cells": len(sub),
            "ml_median": round(float(sub["ml"].abs().median())),
            "dv_median": round(float(sub["dv"].median())),
            "dv_p10": round(float(sub["dv"].quantile(.10))),
            "dv_p90": round(float(sub["dv"].quantile(.90))),
            "ml_p90": round(float(sub["ml"].abs().quantile(.90))),
            "prev_nucleus": sub["nucleus"].value_counts(normalize=True).round(2).to_dict(),
            "top_types": ", ".join(f"{t} {p*100:.0f}%" for t, p in comp.head(4).items()),
        })
    summary = pd.DataFrame(rows).sort_values("dv_median")
    summary.to_csv(OUT / "domain_summary.csv", index=False)
    print("\n=== Spatial domains, ordered ventral to dorsal ===")
    with pd.option_context("display.width", 250, "display.max_colwidth", 62):
        print(summary[["domain", "n_cells", "ml_median", "dv_median", "dv_p10",
                       "dv_p90", "ml_p90", "top_types"]].to_string(index=False))

    print("\n=== Where the previous ellipses placed each domain ===")
    cross = pd.crosstab(adata.obs["domain"], adata.obs["nucleus"], normalize="index").round(2)
    print(cross.to_string())

    adata.write_h5ad(PROC / "hypothalamus_domains.h5ad")
    _plot(adata, summary)
    print(f"\nWrote {PROC / 'hypothalamus_domains.h5ad'} and {OUT}")
    return 0


def _plot(adata, summary) -> None:
    obs = adata.obs
    doms = sorted(obs["domain"].unique(), key=int)
    cmap = plt.get_cmap("tab20")
    colours = {d: cmap(i % 20) for i, d in enumerate(doms)}

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    ax = axes[0]
    sample = obs.sample(min(120000, len(obs)), random_state=0)
    ax.scatter(sample["ml"], sample["dv"], s=.8, alpha=.55,
               c=[colours[d] for d in sample["domain"]], rasterized=True)
    for _, r in summary.iterrows():
        ax.annotate(r["domain"], (r["ml_median"], r["dv_median"]), fontsize=11,
                    fontweight="bold", ha="center",
                    bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=.75))
    ax.set_xlim(-1500, 1500); ax.set_ylim(-150, 1800); ax.set_aspect("equal")
    ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv above ventral surface (µm)")
    ax.set_title("Spatial domains from neighbourhood composition\n(all 12 sections pooled)")

    ax = axes[1]
    prev = {"ARC": "#B4531A", "VMH": "#0E5A61", "DMH": "#6A3D9A", "outside": "#DCE2E2"}
    ax.scatter(sample["ml"], sample["dv"], s=.8, alpha=.55,
               c=[prev[n] for n in sample["nucleus"]], rasterized=True)
    ax.set_xlim(-1500, 1500); ax.set_ylim(-150, 1800); ax.set_aspect("equal")
    ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv above ventral surface (µm)")
    ax.set_title("Previous Gaussian ellipses, for comparison")
    fig.tight_layout()
    fig.savefig(OUT / "domains_vs_ellipses.png", dpi=145)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

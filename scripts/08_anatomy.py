#!/usr/bin/env python3
"""Build the anatomical MBH frame for all 12 sections and assign nuclei.

Replaces the hand-drawn rectangles with a definition every section shares.  Each
section gets a coordinate system from its own third ventricle and ventral
surface; ARC, VMH and DMH are then fitted once, in that shared frame, from
marker-defined anchor cells pooled over all sections.

Validation is against the four hand-drawn sections: the anatomical ARC has to
contain the annotated ARC neurons, and so on.  That check is printed, not
assumed.

Writes results/anatomy/ and data/processed/mbh_anatomical.h5ad.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.anatomy import (  # noqa: E402
    ARC_RULE, DMH_RULE, NUCLEI, VMH_RULE, assign_nuclei, find_frame,
    fit_nucleus_model, fit_nucleus_model_from_labels, to_frame,
)
from spatial_lea.io import (  # noqa: E402
    ANIMAL_META, SECTION_ANIMAL, SUSPECT_SECTIONS, counts_matrix, load_section,
)

OUT = REPO / "results" / "anatomy"
PROC = REPO / "data" / "processed"
MIN_COUNTS, MIN_GENES = 10, 5
ANCHOR_GENES = ["Gpr50", "Spag16", "Agrp", "Pomc", "Adcyap1", "Grp", "Ppp1r17"]


def load_qc(section: str):
    adata = load_section(section)
    keep = (
        (adata.obs["gene_counts"] >= MIN_COUNTS)
        & (adata.obs["n_genes"] >= MIN_GENES)
        & (adata.obs["nucleus_count"] == 1)
    )
    adata = adata[keep.to_numpy()].copy()
    adata.layers["counts"] = adata.X.copy()
    present = [g for g in ANCHOR_GENES if g in adata.var_names]
    counts = pd.DataFrame(
        counts_matrix(adata[:, present]), columns=present, index=adata.obs_names
    )
    return adata, counts


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    sections = list(SECTION_ANIMAL)

    parts, frames, coord_parts, count_parts = [], [], [], []
    print("=== Anatomical frame per section ===")
    print(f"{'section':9s} {'animal':7s} {'group':6s} {'cells':>7s} {'lining':>7s} "
          f"{'elong':>6s} {'ventral direction':>20s}")
    for section in sections:
        adata, counts = load_qc(section)
        frame = find_frame(adata.obsm["spatial"], counts, section)
        coords = to_frame(adata.obsm["spatial"], frame)
        coords.index = adata.obs_names

        adata.obs["dv"] = coords["dv"].to_numpy()
        adata.obs["ml"] = coords["ml"].to_numpy()
        parts.append(adata)
        coord_parts.append(coords)
        count_parts.append(counts)
        frames.append(frame)

        meta = ANIMAL_META[SECTION_ANIMAL[section]]
        flag = "  (flagged distorted)" if section in SUSPECT_SECTIONS else ""
        print(f"{section:9s} {SECTION_ANIMAL[section]:7s} {meta['age_group']:6s} "
              f"{adata.n_obs:7d} {frame.n_lining:7d} {frame.quality:6.2f} "
              f"({frame.ventral_dir[0]:+.2f},{frame.ventral_dir[1]:+.2f}){flag}")

    merged_pre = ad.concat(parts, label="section", keys=sections, index_unique="-", merge="same")
    coords_pre = pd.DataFrame(
        {"ml": merged_pre.obs["ml"].to_numpy(), "dv": merged_pre.obs["dv"].to_numpy()},
        index=merged_pre.obs_names,
    )

    # Fit on the annotated sections if they exist; fall back to raw marker
    # anchors otherwise, so the pipeline still runs on a fresh checkout.
    labels = _annotated_labels(merged_pre)
    if labels is not None and labels.notna().sum() > 1000:
        have = labels.notna().to_numpy()
        model = fit_nucleus_model_from_labels(coords_pre[have], labels[have])
        source = f"annotated cell types ({int(have.sum()):,} labelled cells)"
    else:
        counts_all = pd.concat(count_parts, keys=sections, names=["section", "cell"])
        model = fit_nucleus_model(coords_pre.reset_index(drop=True), counts_all.reset_index(drop=True))
        source = "marker anchors"

    print(f"\n=== Nucleus model from {source}; shared frame (|ml|, dv) in micrometres ===")
    for nucleus in NUCLEI:
        mu, cov = model[nucleus]["mean"], model[nucleus]["cov"]
        sd = np.sqrt(np.diag(cov))
        print(f"  {nucleus:4s} centre |ml|={mu[0]:6.0f}  dv={mu[1]:6.0f}   "
              f"sd |ml|={sd[0]:5.0f}  dv={sd[1]:5.0f}   (n anchor = {model[nucleus]['n_anchor']})")
    np.save(OUT / "nucleus_model.npy", model, allow_pickle=True)

    merged = merged_pre
    nuc = assign_nuclei(coords_pre, model)
    for col in nuc.columns:
        merged.obs[col] = nuc[col].to_numpy()
    merged.obs["in_mbh"] = merged.obs["nucleus"] != "outside"

    print("\n=== Cells per nucleus, by section ===")
    tab = pd.crosstab(merged.obs["section"], merged.obs["nucleus"])
    tab["MBH_total"] = tab[[n for n in NUCLEI if n in tab.columns]].sum(axis=1)
    tab["section_cells"] = merged.obs["section"].value_counts().reindex(tab.index)
    print(tab.to_string())

    _validate(merged)
    _plot(merged, model, frames, sections)

    merged.write_h5ad(PROC / "mbh_anatomical.h5ad")
    tab.to_csv(OUT / "cells_per_nucleus.csv")
    print(f"\nWrote {PROC / 'mbh_anatomical.h5ad'} and tables to {OUT}")
    return 0


def _annotated_labels(merged):
    """Cell-type labels from 04_annotate.py, aligned to ``merged`` (NaN if absent)."""
    annotated = PROC / "mbh_roi_annotated.h5ad"
    if not annotated.exists():
        return None
    import scanpy as sc

    ref = sc.read_h5ad(annotated)
    ref_key = pd.Series(
        ref.obs["cell_type"].astype(str).to_numpy(),
        index=ref.obs["section"].astype(str) + "|" + ref.obs_names.str.split("-").str[0],
    )
    ref_key = ref_key[~ref_key.index.duplicated()]
    merged_key = merged.obs["section"].astype(str) + "|" + merged.obs_names.str.split("-").str[0]
    return pd.Series(ref_key.reindex(merged_key).to_numpy(), index=merged.obs_names)


def _validate(merged) -> None:
    """The anatomical nuclei must contain the cell types that define them."""
    annotated = PROC / "mbh_roi_annotated.h5ad"
    if not annotated.exists():
        print("\n(skipping validation: run 04_annotate.py first)")
        return
    import scanpy as sc

    ref = sc.read_h5ad(annotated)
    # Match on the original per-section cell id, which both objects carry.
    ref_key = pd.Series(
        ref.obs["cell_type"].astype(str).to_numpy(),
        index=ref.obs["section"].astype(str) + "|" + ref.obs_names.str.split("-").str[0],
    )
    merged_key = merged.obs["section"].astype(str) + "|" + merged.obs_names.str.split("-").str[0]
    matched = ref_key.reindex(merged_key).to_numpy()

    df = pd.DataFrame({"cell_type": matched, "nucleus": merged.obs["nucleus"].to_numpy()}).dropna()
    expect = {
        "ARC Agrp/Npy": "ARC", "ARC Pomc": "ARC", "ARC Th/Slc6a3 (TIDA)": "ARC",
        "ARC Tac2/Esr1 (KNDy-like)": "ARC",
        "VMH-like Glut Rasgrf2": "VMH", "VMH-like Glut Calb1": "VMH",
        "VMH-like Glut Tac1": "VMH", "DMH Grp/Ppp1r17": "DMH",
    }
    print("\n=== Validation: where do the annotated marker populations land? ===")
    rows = []
    for cell_type, want in expect.items():
        sub = df[df.cell_type == cell_type]
        if not len(sub):
            continue
        share = (sub.nucleus == want).mean() * 100
        rows.append({"cell_type": cell_type, "expected": want, "n": len(sub),
                     "pct_in_expected": round(share, 1),
                     "pct_outside_MBH": round((sub.nucleus == "outside").mean() * 100, 1)})
    print(pd.DataFrame(rows).to_string(index=False))


def _plot(merged, model, frames, sections) -> None:
    obs = merged.obs
    colours = {"ARC": "#B4531A", "VMH": "#0E5A61", "DMH": "#6A3D9A", "outside": "#D3DADA"}

    ncol = 4
    nrow = int(np.ceil(len(sections) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 3.8 * nrow), squeeze=False)
    for ax, section in zip(axes.ravel(), sections):
        sub = obs[obs["section"] == section]
        ax.scatter(sub["ml"], sub["dv"], s=.7, alpha=.5,
                   c=[colours[n] for n in sub["nucleus"]], rasterized=True)
        ax.axvline(0, color="#444", lw=.8, ls="--")
        ax.set_title(f"{section} ({ANIMAL_META[SECTION_ANIMAL[section]]['age_group']})", fontsize=10)
        ax.set_xlim(-1400, 1400); ax.set_ylim(-150, 1700)
        ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv above ventral surface (µm)")
        ax.set_aspect("equal")
    for ax in axes.ravel()[len(sections):]:
        ax.axis("off")
    fig.suptitle("MBH nuclei in the anatomical frame — 3V at ml = 0, ventral surface at dv = 0", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "nuclei_per_section.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6.5))
    sample = obs.sample(min(40000, len(obs)), random_state=0)
    ax.scatter(sample["ml"].abs(), sample["dv"], s=.6, alpha=.25, color="#C3CCCC", rasterized=True)
    for nucleus in NUCLEI:
        mu, cov = model[nucleus]["mean"], model[nucleus]["cov"]
        vals, vecs = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
        for n_sd in (1, 2):
            ax.add_patch(Ellipse(mu, 2 * n_sd * np.sqrt(vals[-1]), 2 * n_sd * np.sqrt(vals[0]),
                                 angle=angle, fill=False, lw=1.8 if n_sd == 1 else 1.1,
                                 edgecolor=colours[nucleus],
                                 label=f"{nucleus}" if n_sd == 1 else None))
        ax.annotate(nucleus, mu, color=colours[nucleus], fontweight="bold", fontsize=12)
    ax.set_xlabel("|ml| from 3rd ventricle (µm)"); ax.set_ylabel("dv above ventral surface (µm)")
    ax.set_title("Nucleus model, pooled over 12 sections\n(1 and 2 SD contours)")
    ax.set_xlim(0, 1400); ax.set_ylim(-100, 1700); ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "nucleus_model.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

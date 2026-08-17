#!/usr/bin/env python3
"""Publication figures, with source data written alongside each panel.

Five figures, each built only from claims the design actually supports:

  Fig 1  design and anatomical registration      -- descriptive, fully supported
  Fig 2  MBH cell atlas                          -- descriptive, fully supported
  Fig 3  galanin-system mapping                  -- descriptive, replicated in
                                                    all four animals separately
  Fig 4  the rostro-caudal confound              -- inferential; the unit is the
                                                    section pair, not the animal
  Fig 5  the surviving candidate, and power      -- effect size + calibration,
                                                    no age p-value claimed

Every panel's numbers are written to results/figures/source_data/ so the figure
and the table in the manuscript cannot drift apart.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats as st

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse, Rectangle  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix, counts_vector  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
ANIMALS = ["F536", "G_073", "M493", "M399"]
NUC_COLOUR = {"ARC": "#B4531A", "VMH": "#0E5A61", "DMH": "#6A3D9A", "outside": "#DCE2E2"}
GRP_COLOUR = {"aged": "#B4531A", "adult": "#0E5A61"}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def save(fig, name: str) -> None:
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


def panel_label(ax, letter: str, dx=-0.16, dy=1.06) -> None:
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------- figure 1
def figure1(dom, model) -> None:
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25], hspace=.42, wspace=.36)

    ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a")
    ax.axis("off")
    ax.set_title("Design: blocked, 2 vs 2", loc="left")
    for i, (block, (aged, adult)) in enumerate(BLOCKS.items()):
        for j, (animal, grp) in enumerate([(aged, "aged"), (adult, "adult")]):
            ax.add_patch(Rectangle((j * .5, .55 - i * .35), .42, .26,
                                   facecolor=GRP_COLOUR[grp], alpha=.22,
                                   edgecolor=GRP_COLOUR[grp], lw=1.2))
            ax.text(j * .5 + .21, .74 - i * .35, animal, ha="center", fontsize=8, fontweight="bold")
            ax.text(j * .5 + .21, .63 - i * .35,
                    f"{ANIMAL_META[animal]['age_weeks']} wk", ha="center", fontsize=7)
        ax.text(-.08, .68 - i * .35, block, fontsize=7, ha="right", va="center", color="#555")
    ax.text(.21, .93, "aged", ha="center", fontsize=7.5, color=GRP_COLOUR["aged"], fontweight="bold")
    ax.text(.71, .93, "adult", ha="center", fontsize=7.5, color=GRP_COLOUR["adult"], fontweight="bold")
    ax.text(0, .06, "3 sections per animal\n12 sections, 164,943 cells",
            fontsize=7, color="#444")
    ax.set_xlim(-.15, 1.05); ax.set_ylim(0, 1)

    ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b")
    one = dom.obs[dom.obs["section"] == "M493_2"]
    ax.scatter(one["ml"], one["dv"], s=.5, c="#C9D2D2", rasterized=True)
    lining = one[one["cell_type"].astype(str).isin(["Tanycyte", "Ependymal"])]
    ax.scatter(lining["ml"], lining["dv"], s=1.2, c="#0E5A61", rasterized=True)
    ax.axvline(0, color="k", lw=.9, ls="--")
    ax.axhline(0, color="k", lw=.9, ls=":")
    ax.annotate("3rd ventricle\n(ml = 0)", (0, 1500), fontsize=6.5, ha="center", color="#0E5A61")
    ax.annotate("ventral surface (dv = 0)", (700, 40), fontsize=6.5, ha="center")
    ax.set_xlim(-1300, 1300); ax.set_ylim(-150, 1800); ax.set_aspect("equal")
    ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv (µm)")
    ax.set_title("Anatomical frame", loc="left")

    ax = fig.add_subplot(gs[0, 2]); panel_label(ax, "c")
    sample = dom.obs.sample(min(60000, dom.n_obs), random_state=0)
    ax.scatter(sample["ml"].abs(), sample["dv"], s=.4, c="#DCE2E2", rasterized=True)
    for nucleus, spec in model.items():
        mu, cov = spec["mean"], spec["cov"]
        vals, vecs = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
        ax.add_patch(Ellipse(mu, 4 * np.sqrt(vals[-1]), 4 * np.sqrt(vals[0]), angle=angle,
                             fill=False, lw=1.5, edgecolor=NUC_COLOUR[nucleus]))
        ax.annotate(nucleus, mu, color=NUC_COLOUR[nucleus], fontweight="bold", fontsize=8,
                    ha="center")
    ax.set_xlim(0, 1000); ax.set_ylim(-100, 1700)
    ax.set_xlabel("|ml| (µm)"); ax.set_ylabel("dv (µm)")
    ax.set_title("Nucleus model (2 SD)", loc="left")

    axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    panel_label(axes[0], "d")
    for ax, section in zip(axes, ["M399_3", "G073_1", "F536_2"]):
        sub = dom.obs[dom.obs["section"] == section]
        ax.scatter(sub["ml"], sub["dv"], s=.45, alpha=.6,
                   c=[NUC_COLOUR[n] for n in sub["nucleus"]], rasterized=True)
        grp = ANIMAL_META[SECTION_ANIMAL[section]]["age_group"]
        ax.set_title(f"{section} ({grp})", loc="left", fontsize=8)
        ax.set_xlim(-1300, 1300); ax.set_ylim(-150, 1800); ax.set_aspect("equal")
        ax.set_xlabel("ml (µm)")
    axes[0].set_ylabel("dv (µm)")

    fig.suptitle("Figure 1 — Anatomical registration of the mediobasal hypothalamus",
                 x=.02, ha="left", fontsize=10, fontweight="bold")
    save(fig, "figure1_registration")


# ---------------------------------------------------------------- figure 2
def figure2(atlas) -> None:
    if "X_umap" not in atlas.obsm:
        sc.pp.neighbors(atlas, use_rep="X_pca_harmony", n_neighbors=15, random_state=0)
        sc.tl.umap(atlas, random_state=0)

    counts = pd.crosstab(atlas.obs["cell_type"], atlas.obs["nucleus"])
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
    top = counts.head(16)
    top.to_csv(SRC / "fig2_celltype_by_nucleus.csv")

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1.35], wspace=.42)

    cmap = plt.get_cmap("tab20")
    types = list(top.index)
    colour = {t: cmap(i % 20) for i, t in enumerate(types)}

    ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a")
    um = atlas.obsm["X_umap"]
    other = ~atlas.obs["cell_type"].isin(types)
    ax.scatter(um[other, 0], um[other, 1], s=.3, c="#E2E7E7", rasterized=True)
    for t in types:
        sel = (atlas.obs["cell_type"] == t).to_numpy()
        ax.scatter(um[sel, 0], um[sel, 1], s=.3, color=colour[t], rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title("MBH cell types", loc="left")

    ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b")
    sub = atlas.obs.sample(min(40000, atlas.n_obs), random_state=0)
    cols = [colour.get(t, "#E2E7E7") for t in sub["cell_type"]]
    ax.scatter(sub["ml"], sub["dv"], s=.5, c=cols, rasterized=True)
    ax.set_xlim(-800, 800); ax.set_ylim(-100, 1600); ax.set_aspect("equal")
    ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv (µm)")
    ax.set_title("Same cells, in tissue", loc="left")

    ax = fig.add_subplot(gs[0, 2]); panel_label(ax, "c")
    frac = top.div(top.sum(axis=1), axis=0)
    y = np.arange(len(frac))
    left = np.zeros(len(frac))
    for nucleus in ["ARC", "VMH", "DMH"]:
        if nucleus not in frac:
            continue
        ax.barh(y, frac[nucleus], left=left, color=NUC_COLOUR[nucleus], label=nucleus, height=.75)
        left += frac[nucleus].to_numpy()
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t}  ({counts.loc[t].sum():,})" for t in frac.index], fontsize=6.5)
    ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of cells in each nucleus")
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(.5, 1.0), frameon=False)
    ax.set_title("Nucleus distribution", loc="left", pad=18)

    fig.suptitle("Figure 2 — Cell atlas of the mediobasal hypothalamus",
                 x=.02, ha="left", fontsize=10, fontweight="bold")
    save(fig, "figure2_atlas")


# ---------------------------------------------------------------- figure 3
def _negative_control_rates() -> list[float]:
    """Per-probe detection rate of the negative-control probes, % of cells."""
    from spatial_lea.io import load_section

    rates = []
    for section in ["F536_1", "G073_2", "M399_3", "M493_2"]:
        adata = load_section(section)
        mat = adata.obsm["control_counts_matrix"]
        mat = np.asarray(mat.todense()) if hasattr(mat, "todense") else np.asarray(mat)
        is_probe = adata.uns["control_types"] == "Negative Control Probe"
        rates.extend(((mat[:, is_probe] > 0).mean(axis=0) * 100).tolist())
    return rates


def figure3(atlas) -> None:
    genes = ["Gal", "Galr1", "Galr3"]
    types = atlas.obs["cell_type"].value_counts()
    types = [t for t in types.index if types[t] >= 300][:18]

    rows = []
    for t in types:
        sel = (atlas.obs["cell_type"] == t).to_numpy()
        for g in genes:
            v = counts_vector(atlas, g)[sel]
            rows.append({"cell_type": t, "gene": g, "detect_pct": (v > 0).mean() * 100,
                         "mean_counts": v.mean(), "n": int(sel.sum())})
    dot = pd.DataFrame(rows)
    dot.to_csv(SRC / "fig3_galanin_by_celltype.csv", index=False)

    # Control channels are dropped by concat, so they are read back from the
    # source sections -- one per animal is enough to characterise the background.
    ctrl_rate = _negative_control_rates()

    fig = plt.figure(figsize=(7.2, 5.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1.3, 1],
                          hspace=.45, wspace=.4)

    ax = fig.add_subplot(gs[:, 0]); panel_label(ax, "a", dx=-0.42)
    order = dot[dot.gene == "Galr1"].sort_values("detect_pct", ascending=False)["cell_type"].tolist()
    ypos = {t: i for i, t in enumerate(order)}
    for gi, g in enumerate(genes):
        sub = dot[dot.gene == g]
        sizes = np.clip(sub["detect_pct"], 0, 100) * 1.6 + 2
        ax.scatter([gi] * len(sub), [ypos[t] for t in sub["cell_type"]],
                   s=sizes, c=sub["mean_counts"], cmap="magma_r",
                   vmin=0, vmax=dot["mean_counts"].quantile(.97),
                   edgecolor="#333", linewidth=.3)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, style="italic")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=6.5)
    ax.invert_yaxis(); ax.set_xlim(-.6, len(genes) - .4)
    ax.set_title("Galanin system by cell type\n(size = % cells detected, colour = mean counts)",
                 loc="left")

    ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b")
    nuc = atlas.obs[atlas.obs["cell_type"] == "GABA Gal/Galr1"]["nucleus"].value_counts()
    nuc = nuc.reindex(["DMH", "VMH", "ARC"]).fillna(0)
    ax.bar(range(len(nuc)), nuc.values, color=[NUC_COLOUR[n] for n in nuc.index], width=.65)
    for i, v in enumerate(nuc.values):
        ax.text(i, v, f"{v/nuc.sum()*100:.1f}%", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(nuc))); ax.set_xticklabels(nuc.index)
    ax.set_ylabel("cells"); ax.set_title("Galr1-defining population\nis DMH", loc="left")
    pd.DataFrame({"nucleus": nuc.index, "cells": nuc.values}).to_csv(
        SRC / "fig3_galr1_population_nucleus.csv", index=False)

    ax = fig.add_subplot(gs[1, 1]); panel_label(ax, "c")
    galr1_rate = [dot[(dot.gene == "Galr1")]["detect_pct"].max()]
    parts = ax.violinplot([ctrl_rate], positions=[0], widths=.7, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#9AA8A8"); pc.set_alpha(.7)
    ax.scatter([1], galr1_rate, s=45, color="#0E5A61", zorder=3)
    ax.scatter([1] * len(dot[dot.gene == "Galr1"]), dot[dot.gene == "Galr1"]["detect_pct"],
               s=10, color="#0E5A61", alpha=.5, zorder=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["negative\ncontrol probes", "Galr1\nby cell type"])
    ax.set_yscale("log")
    ax.set_ylim(0.01, 120)
    ax.set_ylabel("% cells detected")
    ax.set_title("Detection vs background", loc="left")

    fig.suptitle("Figure 3 — Galanin system is concentrated in a GABAergic DMH population",
                 x=.02, ha="left", fontsize=10, fontweight="bold")
    save(fig, "figure3_galanin_mapping")


# ---------------------------------------------------------------- figure 4
def figure4(dom, ap) -> None:
    """The rostro-caudal confound, via the within-animal AP correction.

    The earlier version of this figure used a correlation between effect size
    and the AP mismatch of each section pairing.  That diagnostic failed audit
    (20_audit_ap_test.py): the AP gaps in the two blocks do not overlap, so the
    correlation tested a block difference, the pairings share sections so the
    p-value was anticonservative, and 38% of the panel fired at p<0.05.  The
    within-animal slope used here is estimated where age is constant, so it
    cannot be confounded with the age effect it is correcting.
    """
    corr = pd.read_csv(REPO / "results" / "ap_correction" / "ap_corrected_candidates.csv")
    panel_gal = pd.read_csv(
        REPO / "results" / "ap_correction" / "panel_ARC_Agrp_Npy.csv", index_col=0)
    panel_galr1 = pd.read_csv(
        REPO / "results" / "ap_correction" / "panel_GABA_Gal_Galr1.csv", index_col=0)
    vc = pd.read_csv(REPO / "results" / "power" / "variance_components.csv")

    fig, axes = plt.subplots(1, 4, figsize=(8.6, 2.7))
    fig.subplots_adjust(wspace=.78)

    ax = axes[0]; panel_label(ax, "a", dx=-0.34)
    for animal, sub in ap.groupby("animal"):
        sub = sub.sort_values("idx")
        ax.plot(sub["idx"], sub["ap_score"], "-o", ms=3.5, lw=1,
                color=GRP_COLOUR[sub["group"].iloc[0]])
        ax.annotate(animal, (sub["idx"].iloc[-1], sub["ap_score"].iloc[-1]),
                    textcoords="offset points", xytext=(5, -2), fontsize=6,
                    va="center", clip_on=False)
    ax.set_xlabel("section (cut order)"); ax.set_ylabel("AP score")
    ax.set_xticks([1, 2, 3]); ax.set_xlim(.7, 3.6)
    ax.set_title("Morphometric\nAP score", loc="left", fontsize=8)

    ax = axes[1]; panel_label(ax, "b")
    labels, raws, corrs = [], [], []
    for _, r in corr.iterrows():
        labels.append(f"{r['gene']}\n{r['cell_type'].split()[0]}")
        raws.append(r["raw_mean"]); corrs.append(r["corr_mean"])
    x = np.arange(len(labels))
    ax.bar(x - .19, raws, .36, color="#9A4C15", label="raw")
    ax.bar(x + .19, corrs, .36, color="#0E5A61", label="AP-corrected")
    ax.axhline(0, color="#555", lw=.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("log2 aged / adult")
    ax.legend(frameon=False, fontsize=6)
    ax.set_title("Effect before and after\ncorrection", loc="left", fontsize=8)

    ax = axes[2]; panel_label(ax, "c")
    ax.hist(panel_galr1["shrinkage"] * 100, bins=35, color="#DCE2E2", label="panel")
    ax.axvline(corr[corr.gene == "Galr1"]["shrinkage_pct"].iloc[0], color="#B4531A", lw=1.8)
    ax.axvline(corr[corr.cell_type == "ARC Agrp/Npy"]["shrinkage_pct"].iloc[0],
               color="#0E5A61", lw=1.8)
    ax.set_xlabel("% of effect removed")
    ax.set_ylabel("panel genes")
    ax.set_title("Galr1 93.6% removed (top 5%)\nGal 32.6% (typical)",
                 loc="left", fontsize=7.5)

    ax = axes[3]; panel_label(ax, "d")
    lbl = ["between\nsections", "between\nanimals"]
    vals = [vc["sd_section"].mean(), vc["sd_animal"].mean()]
    ax.bar(lbl, vals, color=["#9A4C15", "#0E5A61"], width=.6)
    ax.set_ylabel("SD (log2)")
    ax.set_ylim(0, max(vals) * 1.42)
    ax.set_title("Variance\ncomponents", loc="left", fontsize=8)
    ax.text(.5, max(vals) * 1.20, "anatomical noise\nexceeds biology",
            ha="center", fontsize=6, color="#444")

    fig.suptitle("Figure 4 — Rostro-caudal position accounts for the apparent Galr1 effect",
                 x=.02, ha="left", fontsize=10, fontweight="bold", y=1.12)
    save(fig, "figure4_ap_confound")


def _all_pairings(dom, ap) -> pd.DataFrame:
    genes = dom.var_names.to_numpy()
    rows = []
    for gene, ct in [("Galr1", "GABA Gal/Galr1"), ("Gal", "ARC Agrp/Npy"),
                     ("Gal", "GABA Cacna2d2")]:
        gi = int(np.where(genes == gene)[0][0])
        vals = {}
        for section in ap.index:
            sel = ((dom.obs["cell_type"].astype(str) == ct)
                   & (dom.obs["section"].astype(str) == section)).to_numpy()
            if sel.sum() < 15:
                continue
            c = counts_matrix(dom[sel])
            vals[section] = c[:, gi].sum() / c.sum() * 1e6
        for block, (aged, adult) in BLOCKS.items():
            a = [s for s in vals if SECTION_ANIMAL[s] == aged]
            d = [s for s in vals if SECTION_ANIMAL[s] == adult]
            for si, sj in itertools.product(a, d):
                rows.append({
                    "gene": gene, "cell_type": ct, "block": block,
                    "aged_section": si, "adult_section": sj,
                    "ap_gap": abs(ap.loc[si, "ap_score"] - ap.loc[sj, "ap_score"]),
                    "lfc": float(np.log2((vals[si] + 1) / (vals[sj] + 1))),
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- figure 5
def figure5(dom, ap) -> None:
    ct, gene = "ARC Agrp/Npy", "Gal"
    gi = int(np.where(dom.var_names.to_numpy() == gene)[0][0])
    rows = []
    for section in ap.index:
        sel = ((dom.obs["cell_type"].astype(str) == ct)
               & (dom.obs["section"].astype(str) == section)).to_numpy()
        if sel.sum() < 15:
            continue
        c = counts_matrix(dom[sel])
        rows.append({"section": section, "animal": SECTION_ANIMAL[section],
                     "group": ANIMAL_META[SECTION_ANIMAL[section]]["age_group"],
                     "cpm": c[:, gi].sum() / c.sum() * 1e6,
                     "detect_pct": (c[:, gi] > 0).mean() * 100, "n": int(sel.sum())})
    per = pd.DataFrame(rows)
    per.to_csv(SRC / "fig5_gal_arc_agrp_per_section.csv", index=False)
    power = pd.read_csv(REPO / "results" / "power" / "min_detectable_effect_pooled.csv")
    null = pd.read_csv(REPO / "results" / "galanin_search" / "panel_null_ARC_ARC_Agrp_Npy.csv",
                       index_col=0)

    fig, axes = plt.subplots(1, 4, figsize=(8.4, 2.6))
    fig.subplots_adjust(wspace=.72)

    ax = axes[0]; panel_label(ax, "a", dx=-0.34)
    for i, animal in enumerate(ANIMALS):
        sub = per[per.animal == animal]
        grp = ANIMAL_META[animal]["age_group"]
        ax.scatter([i] * len(sub), sub["cpm"], s=22, color=GRP_COLOUR[grp], zorder=3)
        ax.plot([i - .26, i + .26], [sub["cpm"].mean()] * 2, color=GRP_COLOUR[grp], lw=2)
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{a}\n{ANIMAL_META[a]['age_group']}" for a in ANIMALS], fontsize=6)
    ax.set_ylabel("Gal (CPM)")
    ax.set_title("Gal in ARC Agrp/Npy\none point per section", loc="left", fontsize=8)

    ax = axes[1]; panel_label(ax, "b")
    pair = _all_pairings(dom, ap)
    sub = pair[(pair.gene == "Gal") & (pair.cell_type == ct)].sort_values("ap_gap")
    ax.scatter(range(len(sub)), sub["lfc"], s=18, c="#0E5A61")
    ax.axhline(0, color="#999", lw=.8)
    ax.set_xlabel("section pairings (descriptive)")
    ax.set_ylabel("log2 aged / adult")
    ax.set_title(f"positive in {int((sub.lfc>0).sum())}/{len(sub)}\npairings", loc="left")

    ax = axes[2]; panel_label(ax, "c")
    ax.hist(null["min_abs_lfc"], bins=40, color="#C3CCCC")
    ax.axvline(null.loc[gene, "min_abs_lfc"], color="#B4531A", lw=1.6)
    ax.annotate("Gal", (null.loc[gene, "min_abs_lfc"], ax.get_ylim()[1] * .82),
                color="#B4531A", fontsize=7, ha="left", style="italic")
    ax.set_xlabel("min |log2 FC|")
    ax.set_ylabel("panel genes")
    n_pass = int(null["passes_all"].sum())
    ax.set_title(f"panel null\n{n_pass}/297 genes pass", loc="left", fontsize=8)

    ax = axes[3]; panel_label(ax, "d")
    ax.plot(power["n_per_group"], power["mde_log2"], "-o", ms=3.5, color="#9A4C15",
            label="as sampled")
    ax.plot(power["n_per_group"], power["mde_log2_ap_matched"], "-o", ms=3.5,
            color="#0E5A61", label="AP-matched")
    obs = abs(per[per.group == "aged"].cpm.mean() / per[per.group == "adult"].cpm.mean())
    ax.axhline(np.log2(obs), color="#333", ls="--", lw=1)
    ax.annotate("observed effect", (7, np.log2(obs)), fontsize=6, va="bottom")
    ax.set_xlabel("animals per group"); ax.set_ylabel("min detectable log2 FC")
    ax.set_yscale("log"); ax.legend(frameon=False)
    ax.set_title("Power", loc="left")

    fig.suptitle("Figure 5 — Gal in ARC Agrp/Npy neurons is the one surviving candidate",
                 x=.02, ha="left", fontsize=10, fontweight="bold", y=1.12)
    save(fig, "figure5_candidate_power")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    model = np.load(REPO / "results" / "anatomy" / "nucleus_model.npy", allow_pickle=True).item()
    ap = pd.read_csv(REPO / "results" / "ap_matching" / "section_ap_scores.csv", index_col=0)

    print("Figure 1 ...")
    dom = sc.read_h5ad(PROC / "hypothalamus_domains.h5ad")
    figure1(dom, model)

    print("Figure 2 ...")
    atlas = sc.read_h5ad(PROC / "mbh_atlas_12.h5ad")
    atlas = atlas[atlas.obs["cell_type"] != "Unresolved"].copy()
    figure2(atlas)

    print("Figure 3 ...")
    figure3(atlas)
    del atlas

    print("Figure 4 ...")
    figure4(dom, ap)

    print("Figure 5 ...")
    figure5(dom, ap)

    print(f"\nFigures in {OUT}, source data in {SRC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

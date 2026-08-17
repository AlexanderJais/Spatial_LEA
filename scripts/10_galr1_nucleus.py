#!/usr/bin/env python3
"""Galr1 by nucleus and cell type, over 12 sections with per-nucleus denominators.

Two things change relative to 05_galr1.py:

* Abundance is expressed as a share of its own **nucleus**, not of a drawn box.
  The ROI rectangles differed in area by 40%, so a "% of ROI cells" abundance
  carried the drawing decision with it; a "% of DMH cells" abundance does not.
* Three sections per animal instead of one.  Sections are technical replicates,
  so they are pooled per animal for the group contrast -- they buy precision,
  not degrees of freedom.  Their spread is reported separately, as the check on
  whether an effect reproduces within an animal.

The blocked design is unchanged: B1 = F536 vs M493, B2 = G_073 vs M399.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, SUSPECT_SECTIONS, counts_vector  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "galr1_nucleus"
GENES = ["Galr1", "Galr3", "Gal"]
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
MIN_CELLS, MIN_POS = 30, 15


def per_unit(adata, gene: str, by) -> pd.DataFrame:
    counts = counts_vector(adata, gene)
    df = adata.obs[by].copy()
    df["count"] = counts
    df["total"] = adata.obs["gene_counts"].to_numpy()
    g = df.groupby(by, observed=True)
    return g.apply(
        lambda x: pd.Series({
            "n_cells": len(x),
            "n_pos": int((x["count"] > 0).sum()),
            "detect_pct": (x["count"] > 0).mean() * 100,
            "per10k": x["count"].sum() / x["total"].sum() * 1e4,
        }),
        include_groups=False,
    ).reset_index()


def blocked(wide: pd.DataFrame, n: pd.DataFrame, n_pos: pd.DataFrame, metric: str) -> pd.DataFrame:
    keep = (n >= MIN_CELLS).all(axis=1) & (n_pos >= MIN_POS).all(axis=1)
    wide, n, n_pos = wide[keep], n[keep], n_pos[keep]
    out = pd.DataFrame(index=wide.index)
    for block, (aged, adult) in BLOCKS.items():
        out[f"lfc_{block}"] = np.log2((wide[aged] + 0.01) / (wide[adult] + 0.01))
    out["mean_lfc"] = out[[f"lfc_{b}" for b in BLOCKS]].mean(axis=1)
    out["consistent"] = np.sign(out["lfc_B1"]) == np.sign(out["lfc_B2"])
    out["min_abs_lfc"] = out[[f"lfc_{b}" for b in BLOCKS]].abs().min(axis=1)
    return out.join(wide.round(2).add_prefix(f"{metric}_")).sort_values("mean_lfc")


def section_spread(adata, gene: str, key_cols) -> pd.DataFrame:
    """Per-section values, and whether the group gap exceeds within-animal noise.

    With three sections per animal the section-to-section spread is measurable,
    which turns "the two groups differ" into "the two groups differ by more than
    the same animal differs from itself".
    """
    stats = per_unit(adata, gene, key_cols + ["section"])
    stats["animal"] = stats["section"].map(SECTION_ANIMAL)
    stats["age_group"] = stats["animal"].map(lambda a: ANIMAL_META[a]["age_group"])
    rows = []
    for keys, sub in stats.groupby(key_cols, observed=True):
        sub = sub[sub["n_cells"] >= MIN_CELLS]
        if sub["animal"].nunique() < 4:
            continue
        aged = sub[sub.age_group == "aged"]["detect_pct"]
        adult = sub[sub.age_group == "adult"]["detect_pct"]
        within = sub.groupby("animal", observed=True)["detect_pct"].std().mean()
        rows.append({
            **dict(zip(key_cols, keys if isinstance(keys, tuple) else (keys,))),
            "n_sections": len(sub),
            "aged_mean": round(aged.mean(), 2),
            "adult_mean": round(adult.mean(), 2),
            "gap": round(aged.mean() - adult.mean(), 2),
            "within_animal_sd": round(within, 2),
            # How big the group gap is relative to how much one animal's own
            # sections disagree.  Below ~1 the gap is inside the noise.
            "gap_over_noise": round(abs(aged.mean() - adult.mean()) / max(within, 1e-9), 2),
        })
    return pd.DataFrame(rows).sort_values("gap_over_noise", ascending=False)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "mbh_atlas_12.h5ad")
    adata = adata[adata.obs["cell_type"] != "Unresolved"].copy()
    print(f"{adata.n_obs:,} MBH cells, 12 sections, {adata.obs['cell_type'].nunique()} cell types\n")

    # ---- composition, now per nucleus -------------------------------------
    print("=== Composition: cell type as % of its OWN nucleus (was % of drawn ROI) ===")
    comp = (
        adata.obs.groupby(["nucleus", "cell_type", "animal"], observed=True).size()
        .rename("n").reset_index()
    )
    tot = adata.obs.groupby(["nucleus", "animal"], observed=True).size().rename("nucleus_n")
    comp = comp.join(tot, on=["nucleus", "animal"])
    comp["pct"] = comp["n"] / comp["nucleus_n"] * 100
    comp.to_csv(OUT / "composition_per_nucleus.csv", index=False)

    wide = comp.pivot_table(index=["nucleus", "cell_type"], columns="animal", values="pct")
    n_w = comp.pivot_table(index=["nucleus", "cell_type"], columns="animal", values="n")
    comp_res = blocked(wide, n_w, n_w, "pct")
    comp_res.to_csv(OUT / "composition_blocked.csv")
    focus = [i for i in comp_res.index if i[1] in
             ("DMH Grp/Ppp1r17", "GABA Gal/Galr1", "ARC Agrp/Npy", "ARC Pomc", "Tanycyte")]
    with pd.option_context("display.width", 220):
        print(comp_res.loc[focus][["pct_F536", "pct_G_073", "pct_M493", "pct_M399",
                                   "lfc_B1", "lfc_B2", "consistent"]].round(2).to_string())

    # ---- Galr1 by nucleus x cell type -------------------------------------
    all_res = {}
    for gene in GENES:
        stats = per_unit(adata, gene, ["nucleus", "cell_type", "animal"])
        stats.to_csv(OUT / f"{gene}_per_animal.csv", index=False)
        idx = ["nucleus", "cell_type"]
        for metric in ("detect_pct", "per10k"):
            w = stats.pivot_table(index=idx, columns="animal", values=metric)
            n = stats.pivot_table(index=idx, columns="animal", values="n_cells")
            npos = stats.pivot_table(index=idx, columns="animal", values="n_pos")
            res = blocked(w, n, npos, metric)
            res.to_csv(OUT / f"{gene}_{metric}_blocked.csv")
            all_res[(gene, metric)] = res

    print("\n=== Galr1 detection rate by nucleus x cell type (aged vs adult) ===")
    res = all_res[("Galr1", "detect_pct")]
    with pd.option_context("display.width", 230):
        print(res[["detect_pct_F536", "detect_pct_G_073", "detect_pct_M493", "detect_pct_M399",
                   "lfc_B1", "lfc_B2", "mean_lfc", "consistent"]].round(2).to_string())
    k, n = int(res["consistent"].sum()), len(res)
    print(f"  cross-block agreement: {k}/{n}")

    print("\n=== Shortlist: consistent in both blocks and >=1.5x in the weaker one ===")
    any_hit = False
    for (gene, metric), r in all_res.items():
        hits = r[r["consistent"] & (r["min_abs_lfc"] >= 0.58)]
        if len(hits):
            any_hit = True
            print(f"\n{gene} / {metric}:")
            print(hits[["lfc_B1", "lfc_B2", "mean_lfc"]].round(2).to_string())
    if not any_hit:
        print("  (none)")

    # ---- does anything survive the within-animal noise check? --------------
    print("\n=== Galr1: group gap vs within-animal section spread (top 12) ===")
    spread = section_spread(adata, "Galr1", ["nucleus", "cell_type"])
    spread.to_csv(OUT / "Galr1_section_spread.csv", index=False)
    with pd.option_context("display.width", 220):
        print(spread.head(12).to_string(index=False))

    _plot_focus(adata)
    print(f"\nWrote tables and figures to {OUT}")
    return 0


def _plot_focus(adata) -> None:
    """The GABA Gal/Galr1 population, per section, per animal."""
    sub = adata[adata.obs["cell_type"] == "GABA Gal/Galr1"].copy()
    counts = counts_vector(sub, "Galr1")
    df = sub.obs[["section", "animal", "age_group", "nucleus", "gene_counts"]].copy()
    df["pos"] = counts > 0
    df["count"] = counts

    per_section = df.groupby(["animal", "age_group", "section"], observed=True).agg(
        n=("pos", "size"), detect=("pos", "mean"),
        per10k=("count", lambda c: c.sum()), total=("gene_counts", "sum"),
    ).reset_index()
    per_section["detect"] *= 100
    per_section["per10k"] = per_section["per10k"] / per_section["total"] * 1e4
    per_section.to_csv(OUT / "gaba_gal_galr1_per_section.csv", index=False)

    order = ["F536", "G_073", "M493", "M399"]
    colours = {"aged": "#B4531A", "adult": "#0E5A61"}
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5))

    for ax, (col, label) in zip(
        axes[:2], [("detect", "% of cells Galr1+"), ("per10k", "Galr1 counts per 10k")]
    ):
        for i, animal in enumerate(order):
            rows = per_section[per_section.animal == animal]
            grp = rows["age_group"].iloc[0]
            ax.scatter([i] * len(rows), rows[col], s=64, color=colours[grp], zorder=3)
            ax.plot([i - .22, i + .22], [rows[col].mean()] * 2, color=colours[grp], lw=2.4)
        ax.set_xticks(range(4))
        ax.set_xticklabels([f"{a}\n{ANIMAL_META[a]['age_group']}" for a in order], fontsize=9)
        ax.set_ylabel(label)
        ax.set_title(f"GABA Gal/Galr1 (DMH) — {label}\none point per section")
        ax.grid(axis="y", alpha=.25)

    ax = axes[2]
    nuc = adata.obs[adata.obs["cell_type"] == "GABA Gal/Galr1"]["nucleus"].value_counts()
    ax.bar(nuc.index, nuc.values, color=["#6A3D9A", "#0E5A61", "#B4531A"][: len(nuc)])
    ax.set_ylabel("cells"); ax.set_title("Where the GABA Gal/Galr1 population sits")
    for i, v in enumerate(nuc.values):
        ax.text(i, v, f" {v/nuc.sum()*100:.1f}%", ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUT / "gaba_gal_galr1.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

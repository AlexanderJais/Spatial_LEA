#!/usr/bin/env python3
"""Galr1 (and Galr3, Gal) by cell type: the A/B/C decomposition.

"Galr1 changed with age" is three separable claims, and they are tested apart:

  A  a different *fraction of cells* express it        -> detection rate
  B  the *same cells* express more or less             -> counts per 10k, expressing cells
  C  the *abundance of expressing cell types* changed  -> composition

The unit of replication is the animal, so every quantity is computed per animal
first.  The design is blocked (B1 = F536 aged vs M493 adult, B2 = G_073 aged vs
M399 adult), so the estimator is the within-block log2 fold change, one per
block, and the evidence is whether the two blocks agree.

With n=2 per group the exact blocked permutation null has four assignments, so
no per-gene p-value below 0.5 (two-sided) is attainable.  Nothing here reports
one.  What is reported is effect size, per-animal transparency, and cross-block
consistency -- plus a global test of whether more cell types agree in direction
than chance would give.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, counts_vector  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "galr1"
GENES = ["Galr1", "Galr3", "Gal"]
MIN_CELLS_PER_ANIMAL = 30
MIN_POS_CELLS = 15  # expressing cells per animal; below this a fold change is noise
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}  # (aged, adult)


def per_animal_stats(adata, gene: str) -> pd.DataFrame:
    """Detection rate and expression per 10k counts, per cell type per animal."""
    counts = counts_vector(adata, gene)
    total = adata.obs["gene_counts"].to_numpy()
    df = pd.DataFrame(
        {
            "cell_type": adata.obs["cell_type"].to_numpy(),
            "animal": adata.obs["animal"].to_numpy(),
            "count": counts,
            "total": total,
        }
    )
    grouped = df.groupby(["cell_type", "animal"], observed=True)
    out = grouped.apply(
        lambda g: pd.Series(
            {
                "n_cells": len(g),
                "n_pos": int((g["count"] > 0).sum()),
                # A: what fraction of cells carry any transcript
                "detect_pct": (g["count"] > 0).mean() * 100,
                # B: expression normalised to sequencing depth of those cells
                "per10k": g["count"].sum() / g["total"].sum() * 1e4,
                # B, restricted to expressing cells: intensity, not prevalence
                "mean_in_pos": g.loc[g["count"] > 0, "count"].mean() if (g["count"] > 0).any() else np.nan,
            }
        ),
        include_groups=False,
    )
    return out.reset_index()


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson interval -- behaves at k=0, unlike the normal approximation."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (100 * max(0.0, centre - half), 100 * min(1.0, centre + half))


def blocked(stats_df: pd.DataFrame, metric: str) -> pd.DataFrame:
    wide = stats_df.pivot(index="cell_type", columns="animal", values=metric)
    n = stats_df.pivot(index="cell_type", columns="animal", values="n_cells")
    n_pos = stats_df.pivot(index="cell_type", columns="animal", values="n_pos")

    # Two gates.  Enough cells to estimate anything, and enough *expressing*
    # cells that a fold change is not driven by a near-zero denominator: a
    # 0-vs-4-cell contrast produces a spectacular log2 FC that means nothing.
    keep = (n >= MIN_CELLS_PER_ANIMAL).all(axis=1) & (n_pos >= MIN_POS_CELLS).all(axis=1)
    wide, n, n_pos = wide[keep], n[keep], n_pos[keep]

    out = pd.DataFrame(index=wide.index)
    for block, (aged, adult) in BLOCKS.items():
        out[f"lfc_{block}"] = np.log2((wide[aged] + 0.01) / (wide[adult] + 0.01))
    out["mean_lfc"] = out[[f"lfc_{b}" for b in BLOCKS]].mean(axis=1)
    out["consistent"] = np.sign(out["lfc_B1"]) == np.sign(out["lfc_B2"])
    out["min_abs_lfc"] = out[[f"lfc_{b}" for b in BLOCKS]].abs().min(axis=1)
    out = (out.join(wide.round(3).add_prefix(f"{metric}_"))
              .join(n.astype(int).add_prefix("n_"))
              .join(n_pos.astype(int).add_prefix("npos_")))
    return out.sort_values("mean_lfc")


def sign_test(result: pd.DataFrame) -> tuple[int, int, float]:
    """Do more cell types agree across blocks than chance would give?

    Under a null of no age effect the two blocks' signs are independent coin
    flips, so agreement is Binomial(n, 0.5).  This is a statement about the set
    of cell types, not about any one of them.
    """
    k, n = int(result["consistent"].sum()), len(result)
    return k, n, float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    print(f"{adata.n_obs:,} MBH cells, {adata.obs['cell_type'].nunique()} cell types\n")

    all_results = {}
    for gene in GENES:
        stats_df = per_animal_stats(adata, gene)
        stats_df.to_csv(OUT / f"{gene}_per_animal.csv", index=False)

        for metric, title in [("detect_pct", "A: detection rate (% cells)"),
                              ("per10k", "B: counts per 10k")]:
            res = blocked(stats_df, metric)
            res.to_csv(OUT / f"{gene}_{metric}_blocked.csv")
            all_results[(gene, metric)] = res

            k, n, p = sign_test(res)
            print(f"=== {gene} -- {title} ===")
            cols = [f"{metric}_{a}" for a in ["F536", "G_073", "M493", "M399"]] + \
                   ["lfc_B1", "lfc_B2", "mean_lfc", "consistent"]
            with pd.option_context("display.width", 230):
                print(res[cols].round(3).to_string())
            print(f"  cross-block agreement: {k}/{n} cell types "
                  f"(binomial p = {p:.3f}; aged/adult columns are F536,G_073 vs M493,M399)\n")

    _plot(all_results, adata)

    # The headline shortlist: consistent direction in both blocks, and an effect
    # of at least 1.5x in the weaker block.
    print("=== Shortlist: consistent across both blocks, >=0.58 log2 (1.5x) in BOTH ===")
    for (gene, metric), res in all_results.items():
        hits = res[res["consistent"] & (res["min_abs_lfc"] >= 0.58)]
        if len(hits):
            print(f"\n{gene} / {metric}:")
            print(hits[["lfc_B1", "lfc_B2", "mean_lfc"]].round(2).to_string())
    print(f"\nWrote tables and figures to {OUT}")
    return 0


def _plot(all_results, adata) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(19, 9))
    for ax, gene in zip(axes, GENES):
        res = all_results[(gene, "detect_pct")].sort_values("mean_lfc")
        y = np.arange(len(res))
        ax.axvline(0, color="#666", lw=1)
        for i, block in enumerate(BLOCKS):
            ax.scatter(res[f"lfc_{block}"], y + (i - .5) * .22, s=34,
                       color=["#0E5A61", "#B4531A"][i], label=f"block {block}", zorder=3)
        for yi, (_, row) in enumerate(res.iterrows()):
            ax.plot([row["lfc_B1"], row["lfc_B2"]], [yi - .11, yi + .11],
                    color="#B0BDBD" if not row["consistent"] else "#4A5A5A", lw=1.4, zorder=2)
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{'* ' if c else '  '}{t}" for t, c in zip(res.index, res["consistent"])], fontsize=8)
        ax.set_xlabel("log2 FC, aged / adult (detection rate)")
        ax.set_title(f"{gene}\n* = both blocks agree in direction")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(axis="x", alpha=.25)
    fig.tight_layout()
    fig.savefig(OUT / "blocked_lfc_by_celltype.png", dpi=155)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

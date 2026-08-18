#!/usr/bin/env python3
"""The segmentation figure: what a 5 um dilation does, and what changes without it.

Panel (a) is the argument.  Everything else quantifies it.

  a  one 160 um field of DMH tissue, every molecule drawn twice -- coloured by
     the vendor cell it was assigned to, then by the Baysor cell -- with the
     measured DAPI outlines on both.  The nuclei are identical; only what counts
     as belonging to them differs.
  b  cell size, as equivalent diameter, against the sizes real somata have
  c  the mixed GABAergic/glutamatergic rate, which is the artefact
  d  where the Galr1 signal sits, before and after
  e  cell-type composition, before and after

    python scripts/33_baysor_figure.py G073_1 M399_3
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.anatomy import find_frame  # noqa: E402
from spatial_lea.io import RAW, counts_matrix, load_section  # noqa: E402

SEG = REPO / "data" / "segmentation"
PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
RES = REPO / "results" / "baysor"
CROP_UM = 160.0
GLIA = ("Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Tanycyte", "Ependymal")
VENDOR_C, BAYSOR_C = "#9AA8A8", "#0E5A61"
ANCHOR_GENES = ["Gpr50", "Spag16", "Agrp", "Pomc", "Adcyap1", "Grp", "Ppp1r17"]

plt.rcParams.update({
    "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def label(ax, letter, dx=-0.16, dy=1.12):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def frame_of(section: str):
    adata = load_section(section)
    keep = ((adata.obs["gene_counts"] >= 10) & (adata.obs["n_genes"] >= 5)
            & (adata.obs["nucleus_count"] == 1))
    adata = adata[keep.to_numpy()].copy()
    adata.layers["counts"] = adata.X.copy()
    present = [g for g in ANCHOR_GENES if g in adata.var_names]
    counts = pd.DataFrame(counts_matrix(adata[:, present]), columns=present,
                          index=adata.obs_names)
    return find_frame(adata.obs[["x_centroid", "y_centroid"]].to_numpy(), counts, section)


def crop_data(section: str, centre_ml: float, centre_dv: float):
    """Molecules in one field, with both segmentations, plus the DAPI outlines."""
    frame = frame_of(section)
    centre = frame.origin - centre_dv * frame.ventral_dir + centre_ml * frame.lateral_dir
    lo, hi = centre - CROP_UM / 2, centre + CROP_UM / 2

    cols = ["cell_id", "feature_name", "x_location", "y_location", "qv", "is_gene"]
    parts = []
    for b in pq.ParquetFile(RAW / section / "transcripts.parquet").iter_batches(
            batch_size=2_000_000, columns=cols):
        d = b.to_pandas()
        xy = d[["x_location", "y_location"]].to_numpy()
        m = ((xy >= lo).all(axis=1) & (xy <= hi).all(axis=1)
             & d["is_gene"].to_numpy() & (d["qv"].to_numpy() >= 20))
        if m.any():
            parts.append(d[m])
    mol = pd.concat(parts, ignore_index=True)
    for c in ("cell_id", "feature_name"):
        if len(mol) and isinstance(mol[c].iloc[0], bytes):
            mol[c] = mol[c].str.decode("utf-8")

    bay = pd.read_csv(SEG / f"{section}_baysor" / "segmentation.csv",
                      usecols=["x", "y", "cell"])
    bxy = bay[["x", "y"]].to_numpy()
    bay = bay[(bxy >= lo).all(axis=1) & (bxy <= hi).all(axis=1)]

    nb = pd.read_parquet(RAW / section / "nucleus_boundaries.parquet")
    xcol = "vertex_x" if "vertex_x" in nb else nb.columns[1]
    ycol = "vertex_y" if "vertex_y" in nb else nb.columns[2]
    nb = nb[(nb[xcol] >= lo[0]) & (nb[xcol] <= hi[0])
            & (nb[ycol] >= lo[1]) & (nb[ycol] <= hi[1])]
    return mol, bay, nb, (xcol, ycol), lo, hi


def scatter_seg(ax, x, y, ids, lo, hi, title, unassigned):
    """Molecules coloured by cell id; unassigned ones grey."""
    ids = pd.Series(ids).astype(str).to_numpy()
    free = ids == str(unassigned)
    ax.scatter(x[free], y[free], s=1.1, c="#D8DEDE", rasterized=True, linewidths=0)
    codes = pd.factorize(ids[~free])[0]
    rng = np.random.default_rng(0)
    cmap = plt.get_cmap("tab20")(rng.permutation(20))
    ax.scatter(x[~free], y[~free], s=1.6, c=cmap[codes % 20], rasterized=True,
               linewidths=0)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(.6); s.set_color("#B8C2C2")
    ax.set_title(f"{title}\n{(~free).sum() / len(ids) * 100:.0f}% of molecules assigned",
                 loc="left", fontsize=7)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sections", nargs="+")
    p.add_argument("--crop-section", default=None)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    comp_stats = pd.read_csv(RES / "segmentation_comparison.csv")
    cells = {s: pd.read_csv(RES / f"{s}_baysor_cells.csv", index_col=0)
             for s in args.sections}
    ref = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")

    crop_section = args.crop_section or args.sections[0]
    dmh = ref.obs[(ref.obs["section"].astype(str) == crop_section)
                  & (ref.obs["nucleus_ext"] == "DMH")]
    centre_ml, centre_dv = float(dmh["ml"].median()), float(dmh["dv"].median())
    mol, bay, nb, (xcol, ycol), lo, hi = crop_data(crop_section, centre_ml, centre_dv)

    matched = pd.concat([pd.read_csv(RES / f"{s}_mixed_matched.csv")
                         for s in args.sections], ignore_index=True)

    fig = plt.figure(figsize=(7.2, 7.8))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.05, .80, .82],
                          width_ratios=[1, 1, 1.02], hspace=.80, wspace=.58)
    W3 = .26
    C_VEN, C_BAY, C_ALL = VENDOR_C, BAYSOR_C, "#BFD4D4"
    xs = np.arange(len(comp_stats))
    sec_labels = [f"{r.section}\n{r.age_group}" for _, r in comp_stats.iterrows()]

    # (a) the same tissue, twice
    for i, (x, y, ids, title, un) in enumerate([
            (mol["x_location"].to_numpy(), mol["y_location"].to_numpy(),
             mol["cell_id"], "Vendor: nucleus + 5 µm", "UNASSIGNED"),
            (bay["x"].to_numpy(), bay["y"].to_numpy(),
             bay["cell"].fillna(""), "Baysor", "")]):
        ax = fig.add_subplot(gs[0, i])
        if i == 0:
            label(ax, "a", dx=-0.08, dy=1.26)
        scatter_seg(ax, x, y, ids, lo, hi, title, un)
        for _, g in nb.groupby(nb.columns[0]):
            ax.plot(np.r_[g[xcol], g[xcol].iloc[:1]],
                    np.r_[g[ycol], g[ycol].iloc[:1]],
                    lw=.45, color="#1A1A1A", alpha=.85)
        if i == 0:
            ax.plot([lo[0] + 12, lo[0] + 32], [lo[1] + 12] * 2, lw=1.6, color="k")
            ax.text(lo[0] + 22, lo[1] + 17, "20 µm", ha="center", fontsize=5.5)
            ax.annotate(f"{crop_section}, dorsomedial hypothalamus. Black outlines are\n"
                        "the measured DAPI nuclei — identical in both panels;\n"
                        "only what counts as belonging to them differs.",
                        xy=(0, -.09), xycoords="axes fraction", fontsize=6,
                        color="#5A6A6A", va="top", linespacing=1.5)

    # (b) mixed identity, three ways
    ax = fig.add_subplot(gs[0, 2]); label(ax, "b", dx=-0.34, dy=1.26)
    for off, vals, col in ((-W3, comp_stats["vendor_mixed_pct"], C_VEN),
                           (0, matched["baysor_matched"], C_BAY),
                           (W3, matched["baysor_all"], C_ALL)):
        ax.bar(xs + off, vals, W3, color=col)
        for j, v in enumerate(vals):
            ax.text(j + off, v + .5, f"{v:.0f}", ha="center", fontsize=5.5)
    ax.set_xticks(xs); ax.set_xticklabels(sec_labels, fontsize=6)
    ax.set_ylim(0, 31)
    ax.set_ylabel("% of marker⁺ cells")
    ax.set_title("Cells carrying GABAergic and\nglutamatergic markers at once", loc="left")
    ax.annotate("Middle bar is the comparison that means anything:\n"
                "matched objects, matched depth.",
                xy=(0, -.24), xycoords="axes fraction", fontsize=5.8,
                color="#5A6A6A", va="top")

    # (c) cell size
    ax = fig.add_subplot(gs[1, 0]); label(ax, "c", dx=-0.24)
    s0 = args.sections[0]
    ven = ref[(ref.obs["section"].astype(str) == s0)]
    ven_d = 2 * np.sqrt(ven.obs["cell_area"].to_numpy() / np.pi)
    bay_d = 2 * np.sqrt(cells[s0].loc[cells[s0]["prior_cell"] > 0, "area"].to_numpy() / np.pi)
    bins = np.linspace(0, 36, 55)
    ax.hist(ven_d, bins=bins, color=C_VEN, alpha=.85, density=True)
    ax.hist(bay_d, bins=bins, histtype="step", lw=1.5, color=C_BAY, density=True)
    ytop = ax.get_ylim()[1]
    ax.axvspan(15, 25, color="#B4531A", alpha=.10)
    ax.axvspan(8, 12, color="#6A3D9A", alpha=.10)
    ax.text(20, ytop * .97, "neuron somata", ha="center", va="top", fontsize=5.5,
            color="#B4531A", rotation=90)
    ax.text(10, ytop * .97, "glia", ha="center", va="top", fontsize=5.5,
            color="#6A3D9A", rotation=90)
    ax.set_xlabel("equivalent cell diameter (µm)"); ax.set_ylabel("density")
    ax.set_title(f"Cell size, {s0}", loc="left")

    # (d) depth -- why the pale bar in (b) is lower
    ax = fig.add_subplot(gs[1, 1]); label(ax, "d", dx=-0.24)
    for off, vals, col in ((-W3, matched["median_counts_vendor"], C_VEN),
                           (0, matched["median_counts_matched"], C_BAY),
                           (W3, comp_stats["baysor_median_counts"], C_ALL)):
        ax.bar(xs + off, vals, W3, color=col)
    ax.set_xticks(xs); ax.set_xticklabels(sec_labels, fontsize=6)
    ax.set_ylabel("median transcripts / cell")
    ax.set_title("Depth of those same cells", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (C_VEN, C_BAY, C_ALL)]
    ax.legend(handles, ["vendor cells",
                        "Baysor cells on a measured nucleus",
                        "all Baysor cells, incl. nucleus-free"],
              frameon=False, fontsize=5.8, loc="upper center",
              bbox_to_anchor=(.5, -.26), ncol=1, handlelength=1.1,
              handletextpad=.4, labelspacing=.35)

    # (e) Galr1 by nucleus
    ax = fig.add_subplot(gs[1, 2]); label(ax, "e", dx=-0.28)
    nucs = [c[len("galr1_pct_"):-len("_vendor")] for c in comp_stats.columns
            if c.startswith("galr1_pct_") and c.endswith("_vendor")]
    v = comp_stats[[f"galr1_pct_{n}_vendor" for n in nucs]].mean().to_numpy()
    b = comp_stats[[f"galr1_pct_{n}_baysor" for n in nucs]].mean().to_numpy()
    order = np.argsort(-v)[:7]
    nucs = [nucs[i] for i in order]; v, b = v[order], b[order]
    nx = np.arange(len(nucs))
    ax.bar(nx - .19, v, .38, color=C_VEN)
    ax.bar(nx + .19, b, .38, color=C_BAY)
    ax.set_xticks(nx); ax.set_xticklabels(nucs, fontsize=6, rotation=40, ha="right")
    ax.set_ylabel("% of all Galr1 signal")
    ax.set_title("Galr1 localisation is unchanged", loc="left")
    pd.DataFrame({"nucleus": nucs, "vendor": v.round(2),
                  "baysor": b.round(2)}).to_csv(SRC / "baysor_fig_e_galr1.csv", index=False)

    # (f) composition, matched set
    ax = fig.add_subplot(gs[2, :2]); label(ax, "f", dx=-0.09)
    comp = pd.read_csv(RES / f"{s0}_composition_matched.csv", index_col=0)
    comp = comp.rename(columns={"delta_matched": "delta"})
    tops = comp.reindex(comp["delta"].abs().sort_values(ascending=False).index).head(10)
    tops = tops.sort_values("delta")
    ys = np.arange(len(tops))
    ax.barh(ys, tops["delta"], color=[C_BAY if d > 0 else "#B4531A" for d in tops["delta"]])
    ax.axvline(0, color="k", lw=.6)
    ax.set_yticks(ys); ax.set_yticklabels(tops.index, fontsize=6)
    ax.set_xlim(-11, 11)
    ax.set_xlabel("change in share of cells (percentage points), Baysor − vendor")
    ax.set_title(f"Composition, matched set ({s0})", loc="left")
    ax.text(10.4, .2, "every shift under 1.3 pp.\nUsing all Baysor cells\ninstead moves"
            "\nastrocytes by +9.8 pp —\nthat gain is fragments.",
            fontsize=5.5, color="#7A8A8A", ha="right", va="bottom")

    # (g) mixed rate against depth
    ax = fig.add_subplot(gs[2, 2]); label(ax, "g", dx=-0.34, dy=1.22)
    d = pd.read_csv(RES / f"{s0}_mixed_by_depth.csv", index_col=0)
    dx = np.arange(len(d))
    ax.plot(dx, d["vendor_pct_mixed"], "o-", color=C_VEN, ms=3, lw=1.2, label="vendor")
    ax.plot(dx, d["baysor_pct_mixed"], "o-", color=C_BAY, ms=3, lw=1.2, label="Baysor")
    ax.set_xticks(dx)
    ax.set_xticklabels([i.split(",")[0].strip("[") for i in d.index], fontsize=5.5)
    ax.set_xlabel("transcripts per cell (bin start)")
    ax.set_ylabel("% of marker⁺ cells mixed")
    ax.set_title("Mixing rises steeply with depth;\nat matched depth Baysor is lower", loc="left")
    ax.legend(frameon=False, loc="upper left")

    fig.suptitle("Re-segmentation with Baysor", x=.02, ha="left", y=1.00,
                 fontsize=10, fontweight="bold")
    fig.savefig(OUT / "baysor_segmentation.png")
    fig.savefig(OUT / "baysor_segmentation.pdf")
    plt.close(fig)
    comp_stats.to_csv(SRC / "baysor_fig_summary.csv", index=False)
    print(f"wrote {OUT / 'baysor_segmentation.png'} (+ .pdf)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

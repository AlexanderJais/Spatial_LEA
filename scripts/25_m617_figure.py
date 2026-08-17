#!/usr/bin/env python3
"""The Xenium figure for the M617 paper.

One main figure.  In a treatment paper the spatial data is supporting material:
it establishes what M617 can act on and predicts the sign of the effect, so the
functional result can be interpreted.  Every panel is descriptive, replicated
across four animals, and free of the n=2 limitation that blocks any age claim.

  a  where the receptor is, in tissue
  b  receptor share by nucleus -- DMH and LHA dominate, ARC is nearly empty
  c  receptor by cell type, against the negative-control floor
  d  transmitter identity of the receptor-bearing cells -> predicts disinhibition
  e  Galr3 co-expression -> where Galr1 selectivity matters least
  f  target availability is comparable between ages
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

from spatial_lea.io import ANIMAL_META, counts_matrix, load_section  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
ANIMALS = ["F536", "G_073", "M493", "M399"]
NUC_COLOUR = {"DMH": "#6A3D9A", "LHA": "#2E7D32", "VMH": "#0E5A61",
              "DHA_PH": "#5D6D7E", "ZI": "#C2185B", "ARC": "#B4531A", "ME_3V": "#E0A468"}
GABA = ["Gad1", "Gad2", "Slc32a1"]
GLUT = ["Slc17a6", "Slc17a7"]

plt.rcParams.update({
    "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def label(ax, letter, dx=-0.22, dy=1.10):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def control_floor() -> float:
    """Median per-probe detection rate of the negative-control probes."""
    rates = []
    for section in ["F536_1", "G073_2", "M399_3", "M493_2"]:
        ad = load_section(section)
        mat = ad.obsm["control_counts_matrix"]
        mat = np.asarray(mat.todense() if hasattr(mat, "todense") else mat)
        is_probe = ad.uns["control_types"] == "Negative Control Probe"
        rates.extend(((mat[:, is_probe] > 0).mean(axis=0) * 100).tolist())
    return float(np.median(rates))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(NUCLEI)].copy()

    need = ["Galr1", "Galr3", "Gal"] + GABA + GLUT
    present = [g for g in need if g in adata.var_names]
    m = pd.DataFrame(counts_matrix(adata[:, present]), columns=present)
    obs = adata.obs.reset_index(drop=True).copy()
    obs["galr1"] = m["Galr1"]; obs["galr3"] = m["Galr3"]; obs["gal"] = m["Gal"]

    gaba = m[[g for g in GABA if g in m]].sum(axis=1)
    glut = m[[g for g in GLUT if g in m]].sum(axis=1)
    nt = pd.Series("ambiguous", index=obs.index)
    nt[(gaba >= 2) & (glut == 0)] = "GABAergic"
    nt[(glut >= 2) & (gaba == 0)] = "glutamatergic"
    nt[(gaba < 2) & (glut < 2)] = "other"
    obs["nt"] = nt
    total_galr1 = obs["galr1"].sum()
    floor = control_floor()

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1],
                          width_ratios=[1.25, .95, 1.0], hspace=.62, wspace=.60)

    # (a) where the receptor is
    ax = fig.add_subplot(gs[0, 0]); label(ax, "a")
    bg = obs.sample(min(45000, len(obs)), random_state=0)
    ax.scatter(bg["ml"], bg["dv"], s=.5, c="#E4E9E9", rasterized=True)
    pos = obs[obs.galr1 > 0].sample(min(14000, int((obs.galr1 > 0).sum())), random_state=0)
    ax.scatter(pos["ml"], pos["dv"], s=.7, alpha=.75,
               c=[NUC_COLOUR[n] for n in pos["nucleus_ext"]], rasterized=True)
    ax.set_xlim(-1450, 1450); ax.set_ylim(-150, 1800); ax.set_aspect("equal")
    ax.set_xlabel("mediolateral (µm)"); ax.set_ylabel("dorsoventral (µm)")
    ax.set_title("Galr1⁺ cells, 12 sections", loc="left")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=3.5, color=NUC_COLOUR[n], label=n)
               for n in ("DMH", "LHA", "VMH", "ZI", "DHA_PH", "ARC", "ME_3V")]
    ax.legend(handles=handles, ncol=2, frameon=False, fontsize=5.5,
              loc="upper center", bbox_to_anchor=(.5, -.24), handletextpad=.2,
              columnspacing=.8)

    # (b) share by nucleus
    ax = fig.add_subplot(gs[0, 1]); label(ax, "b")
    share = (obs.groupby("nucleus_ext", observed=True)["galr1"].sum()
             / total_galr1 * 100).sort_values()
    ax.barh(range(len(share)), share.values,
            color=[NUC_COLOUR[n] for n in share.index])
    ax.set_yticks(range(len(share))); ax.set_yticklabels(share.index)
    ax.set_xlabel("% of all Galr1 signal")
    ax.set_title("DMH and LHA carry 60%;\nARC carries 0.7%", loc="left")
    share.round(2).to_csv(SRC / "m617_fig_b_nucleus_share.csv")

    # (c) by cell type, against the control floor
    ax = fig.add_subplot(gs[0, 2]); label(ax, "c")
    obs["pop"] = obs["nucleus_ext"].astype(str) + " | " + obs["cell_type"].astype(str)
    tab = obs.groupby("pop", observed=True).agg(
        n=("galr1", "size"), pct=("galr1", lambda x: (x > 0).mean() * 100),
        share=("galr1", lambda x: x.sum() / total_galr1 * 100))
    tab = tab[tab.n >= 150].nlargest(12, "share").sort_values("share")
    ax.barh(range(len(tab)), tab["pct"], color="#0E5A61")
    ax.axvline(floor, color="#B4531A", lw=1.2, ls="--")
    ax.annotate(f"neg-control floor {floor:.2f}%", (floor, -1.4), fontsize=5.5,
                color="#B4531A", ha="left", annotation_clip=False)
    ax.set_ylim(-2.2, len(tab) - .3)
    ax.set_yticks(range(len(tab)))
    ax.set_yticklabels([p[:30] for p in tab.index], fontsize=5.5)
    ax.set_xlabel("% of cells Galr1⁺"); ax.set_xscale("log")
    ax.set_title("Top targets vs background", loc="left")
    tab.round(2).to_csv(SRC / "m617_fig_c_top_targets.csv")

    # (d) transmitter identity
    ax = fig.add_subplot(gs[1, 0]); label(ax, "d")
    nt_share = (obs.groupby("nt", observed=True)["galr1"].sum() / total_galr1 * 100)
    order = [k for k in ["GABAergic", "glutamatergic", "ambiguous", "other"] if k in nt_share]
    cols = {"GABAergic": "#B4531A", "glutamatergic": "#0E5A61",
            "ambiguous": "#9AA8A8", "other": "#DCE2E2"}
    ax.bar(range(len(order)), [nt_share[k] for k in order],
           color=[cols[k] for k in order], width=.65)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["GABA", "Glut", "ambig.", "other"][:len(order)], fontsize=6.5)
    ax.set_ylabel("% of all Galr1 signal")
    ratio = nt_share.get("GABAergic", 0) / max(nt_share.get("glutamatergic", 1e-9), 1e-9)
    ax.set_title(f"Gi-coupled receptor sits {ratio:.1f}×\nmore on GABAergic cells", loc="left")

    # (e) Galr3 overlap
    ax = fig.add_subplot(gs[1, 1]); label(ax, "e")
    ov = obs.groupby("pop", observed=True).agg(
        n=("galr1", "size"),
        g1=("galr1", lambda x: (x > 0).mean() * 100),
        both=("galr1", "size"))
    both = obs.groupby("pop", observed=True).apply(
        lambda x: ((x.galr1 > 0) & (x.galr3 > 0)).mean() * 100, include_groups=False)
    ov["both"] = both
    ov = ov[(ov.n >= 150) & (ov.g1 >= 10)].nlargest(6, "both").sort_values("both")
    ax.barh(range(len(ov)), ov["g1"], color="#DCE2E2", label="Galr1⁺")
    ax.barh(range(len(ov)), ov["both"], color="#C2185B", label="Galr1⁺Galr3⁺")
    ax.set_yticks(range(len(ov)))
    ax.set_yticklabels([p[:26] for p in ov.index], fontsize=5.5)
    ax.set_xlabel("% of cells"); ax.legend(frameon=False, loc="lower right")
    ax.set_title("Where Galr1 selectivity\nmatters least", loc="left")
    ov.round(2).to_csv(SRC / "m617_fig_e_galr3_overlap.csv")

    # (f) target availability by age
    ax = fig.add_subplot(gs[1, 2]); label(ax, "f")
    per = (obs.groupby(["nucleus_ext", "animal"], observed=True)["galr1"]
           .apply(lambda x: (x > 0).mean() * 100).unstack())
    per = per.loc[share.index[::-1]]
    y = np.arange(len(per))
    for i, animal in enumerate(ANIMALS):
        if animal not in per:
            continue
        grp = ANIMAL_META[animal]["age_group"]
        ax.scatter(per[animal], y + (i - 1.5) * .17, s=18,
                   color="#B4531A" if grp == "aged" else "#0E5A61",
                   label=grp if i in (0, 2) else None)
    ax.set_yticks(y); ax.set_yticklabels(per.index, fontsize=6)
    ax.set_xlabel("% of cells Galr1⁺")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.22), ncol=2)
    ax.set_xlim(0, 48)
    ax.set_title("Target availability,\none point per animal", loc="left")
    ax.grid(axis="x", alpha=.25)
    per.round(2).to_csv(SRC / "m617_fig_f_availability_by_age.csv")

    fig.suptitle("Xenium — the Galr1 target landscape M617 acts on",
                 x=.02, ha="left", fontsize=10, fontweight="bold")
    fig.savefig(OUT / "m617_figure_xenium.png")
    fig.savefig(OUT / "m617_figure_xenium.pdf")
    plt.close(fig)
    print(f"wrote {OUT / 'm617_figure_xenium.png'} (+ .pdf), source data in {SRC}")
    print(f"\nnegative-control floor: {floor:.3f}% | GABA:Glut Galr1 share ratio {ratio:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Figure 1: which MBH cells upregulate Galr1 with age.

The question the study asks, in the order a reader needs it answered:

    Galr1 is present in the mediobasal hypothalamus.  Which cells carry it, and
    which of them change between adult and aged mice?

    Of every cell type the design can test, one moves: a glutamatergic
    Otp+/Cbln1+/Prdm8+ population of the DMH, which gains Galr1 and Ghsr with
    age.  The population that carries the most Galr1 -- GABA Gal/Galr1 -- does
    not move, so this is not "the aged hypothalamus has more Galr1".

Panels:
  a  one representative section: every cell, the Galr1+ cells among them, and
     this population outlined on top, inside the registered subregions
  b  which cell types carry Galr1 at all, as the fraction of their cells that
     are positive -- the landscape the screen is run over
  c  the screen: Galr1 aged versus adult in every cell type with enough cells,
     the same statistic and thresholds everywhere.  Whatever passes is named
  d  where the population that moved sits: enrichment by subregion, DMH highest
  e  what it is: its markers against the rest of the window, and its HypoMap
     correspondence
  f  Galr1 in that population, one point per animal, paired within block
  g  Ghsr, the second receptor these neurons carry, behaves the same way

Choices made by rule rather than by eye:
  * representative section = the section with the highest whole-section bilateral
    balance, i.e. the most intact tissue.  Balance is the ratio of cells left
    and right of the midline; a torn or folded section scores low.  Restricting
    that test to the analysed window is not enough -- a section can be perfectly
    symmetric inside the window and badly damaged outside it.
  * a cell type enters the screen when every animal contributes at least
    MIN_CELLS_PER_ANIMAL of it.  The comparison is blocked within animal pairs,
    so a type one animal lacks cannot be tested at all, and a type thin in one
    animal is tested mostly on that animal's noise.  The rule is a property of
    the design, not of any result: it is applied before the statistics are
    looked at.
  * panel a tints the subregions largest first.  Whichever region is drawn last
    takes the shared pixels, and in the order the colour table happens to list
    them the DMH -- the region panel d makes a claim about -- was painted over
    by LHA, ZI and DHA_PH, its three larger neighbours, at exactly the borders
    the parcellation is least sure of.  Size order is a rule that applies to
    every region rather than a thumb on the scale for this one.
  * panels c, d and e are computed here, not read from a stored table, so the
    whole figure regenerates from the data it claims to show.

Adult is plotted first throughout, as the reference condition.  Fos falls
steeply in these neurons in aged mice, but these animals are untreated, so
immediate-early gene expression is hard to interpret without a stimulus; it is
reported in Extended Data rather than as a main-figure claim.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.figstyle import (  # noqa: E402
    ADULT, AGED as C_AGED, FULL, GALR1, GALR1_LOW, INK, NUCLEUS_COLOUR,
    NUCLEUS_LABEL, POP, TISSUE, bare, panel, scalebar, use_style,
)
from spatial_lea.io import (  # noqa: E402
    ADULT as A_ADULT, AGED as A_AGED, BLOCKS, ONE_PER_MOUSE, counts_matrix,
)

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
HYPOMAP = REPO / "results" / "hypomap" / "hypomap_C185_named_all_celltypes.csv"
POPNAME = "Glut Prdm8/Cbln1"
# The type that carries the most Galr1 in this tissue.  Named in the figure as
# the comparison the result has to survive: if ageing simply raised Galr1, this
# is the population it would show up in first.
GALR1_RICH = "GABA Gal/Galr1"
# Adult first: the reference condition.
ANIMALS = list(A_ADULT) + list(A_AGED)
# Priority order for this study; Galr1 leads, then the peptide that names the
# population, then the other two receptors it carries.
KEY_GENES = ["Galr1", "Ghrh", "Ghsr", "Galr3"]
MAIN_GENES = ["Galr1", "Ghsr"]
RECEPTORS = ("Galr1", "Ghsr", "Galr3")
MARKERS_SHOWN = 8
LANDSCAPE_SHOWN = 12
# A cell type enters the screen only when every animal contributes at least this
# many of it.  The test is blocked within animal pairs, so a type missing from
# one animal cannot be tested at all, and a type that is thin in one animal is
# tested mostly on that animal's noise.
MIN_CELLS_PER_ANIMAL = 30


def blocked_stats(cpm: pd.DataFrame) -> pd.DataFrame:
    lfc = pd.DataFrame({b: cpm.loc[x] - cpm.loc[y] for b, (x, y) in BLOCKS.items()}).T
    mean = lfc.mean()
    agree = (np.sign(lfc) == np.sign(mean)).sum()
    vals = cpm.loc[ANIMALS].to_numpy()
    n_ad = len(A_ADULT)
    idx = list(range(len(ANIMALS)))
    obs = vals[n_ad:].mean(axis=0) - vals[:n_ad].mean(axis=0)
    null = np.array([vals[list(c)].mean(axis=0)
                     - vals[[i for i in idx if i not in c]].mean(axis=0)
                     for c in combinations(idx, n_ad)])
    p = (np.abs(null) >= np.abs(obs)).mean(axis=0)
    aged_v, adult_v = vals[n_ad:], vals[:n_ad]
    margin = np.maximum(aged_v.min(axis=0) - adult_v.max(axis=0),
                        adult_v.min(axis=0) - aged_v.max(axis=0))
    out = pd.DataFrame({"lfc": mean, "blocks": agree, "p": p, "margin": margin})
    for b in BLOCKS:
        out[f"lfc_{b}"] = lfc.loc[b]
    return out


def regional_enrichment(is_pop: np.ndarray, nucleus: np.ndarray) -> pd.DataFrame:
    """Each subregion's share of this population, over the population's share of
    the whole window.

    The denominator is the entire analysed window, unnamed tuberal territory
    included: the question is where the population concentrates within the
    tissue this study looked at, not within the named nuclei only.
    """
    overall = float(is_pop.mean())
    rows = []
    for name in NUCLEUS_LABEL:
        m = nucleus == name
        if not m.any():
            continue
        n_pop, n_region = int((is_pop & m).sum()), int(m.sum())
        rows.append({"nucleus": name, "n_pop": n_pop, "n_region": n_region,
                     "pct_of_region": round(n_pop / n_region * 100, 3),
                     "enrichment": round((n_pop / n_region) / overall, 3)})
    return pd.DataFrame(rows).sort_values("enrichment", ascending=False)


def galr1_by_celltype(counts: np.ndarray, var: np.ndarray, cell_types: np.ndarray,
                      animals: np.ndarray) -> pd.DataFrame:
    """Galr1, measured and tested by the same rule in every cell type the design
    can support.

    ``pct_pos`` and ``cpm`` are the landscape -- who carries the receptor.
    ``lfc``/``blocks``/``margin``/``p`` are the age test.  Both come from one
    pass so the panel that asks "which cells carry Galr1" and the panel that
    asks "which of them change" cannot be built on different cell sets.
    """
    gene = int(np.flatnonzero(var == "Galr1")[0])
    rows = []
    for name in pd.unique(cell_types):
        if str(name).startswith(("unlabelled", "unresolved")):
            continue
        sel = cell_types == name
        per_animal = [int((sel & (animals == a)).sum()) for a in ANIMALS]
        if min(per_animal) < MIN_CELLS_PER_ANIMAL:
            continue
        mat = pd.DataFrame({a: counts[sel & (animals == a)].sum(axis=0)
                            for a in ANIMALS}, index=var).T
        cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
        r = blocked_stats(cpm).loc["Galr1"]
        pooled = counts[sel]
        rows.append({
            "cell_type": name, "n_cells": int(sel.sum()),
            "min_cells": min(per_animal),
            "pct_pos": round(float((pooled[:, gene] > 0).mean() * 100), 1),
            "cpm": round(float(np.log2(pooled[:, gene].sum() / pooled.sum() * 1e6 + 1)), 2),
            "lfc": round(float(r.lfc), 3), "fold": round(float(2 ** r.lfc), 2),
            "blocks": int(r.blocks), "p": round(float(r.p), 4),
            "margin": round(float(r.margin), 3),
        })
    return pd.DataFrame(rows).sort_values("lfc", ascending=False)


def population_markers(counts: np.ndarray, var: np.ndarray,
                       is_pop: np.ndarray) -> pd.Series:
    """Log2 enrichment of every gene in the population against the rest of the
    window, which is what the identity claim rests on."""
    inside = counts[is_pop].sum(axis=0) / counts[is_pop].sum() * 1e6
    outside = counts[~is_pop].sum(axis=0) / counts[~is_pop].sum() * 1e6
    return pd.Series(np.log2((inside + 1) / (outside + 1)), index=var)


def hypomap_match() -> tuple[str, float] | None:
    """Best HypoMap C185 correspondence for this population, if 44_query_hypomap
    has been run."""
    if not HYPOMAP.exists():
        return None
    hm = pd.read_csv(HYPOMAP, index_col=0)
    if POPNAME not in hm.columns:
        return None
    best = hm[POPNAME].sort_values(ascending=False)
    return best.index[0].split(": ", 1)[-1], float(best.iloc[0])


def paired_panel(ax, values: dict, ylabel: str, title: str) -> None:
    """One point per animal; adult on the left, joined to its block partner."""
    for a_aged, a_adult in BLOCKS.values():
        ax.plot([0, 1], [values[a_adult], values[a_aged]], color="#C9C9C9",
                lw=.6, zorder=1)
    for j, (members, colour) in enumerate(((A_ADULT, ADULT), (A_AGED, C_AGED))):
        ys = [values[a] for a in members]
        ax.scatter([j] * len(ys), ys, s=13, color=colour, zorder=3,
                   linewidths=0, clip_on=False)
        ax.plot([j - .22, j + .22], [np.mean(ys)] * 2, color=colour, lw=1.3,
                solid_capstyle="butt", zorder=2)
    ax.set_xlim(-.45, 1.45); ax.set_xticks([0, 1])
    ax.set_xticklabels(["adult", "aged"])
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=3)


def main() -> int:
    use_style()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    win = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    win = win[win.obs["section"].astype(str).isin(ONE_PER_MOUSE)].copy()
    cell_types = win.obs["cell_type"].astype(str).to_numpy()
    if POPNAME not in set(cell_types):
        raise SystemExit(
            f"'{POPNAME}' is not among this run's cell types. 04_annotate.py "
            "matches clusters to labels on marker evidence; check its report at "
            "results/annotation/cluster_label_matching.csv before going further.")
    is_pop = cell_types == POPNAME
    animals = win.obs["animal"].astype(str).to_numpy()
    sections = win.obs["section"].astype(str).to_numpy()
    ml, dv = win.obs["ml"].to_numpy(), win.obs["dv"].to_numpy()
    nuc = win.obs["nucleus_ext"].astype(str).to_numpy()
    counts = counts_matrix(win)
    var = win.var_names.to_numpy()
    galr1 = counts[:, int(np.flatnonzero(var == "Galr1")[0])]

    mat = pd.DataFrame({a: counts[is_pop & (animals == a)].sum(axis=0)
                        for a in ANIMALS}, index=var).T
    expressed = mat.sum(axis=0) >= 200
    cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
    stats = blocked_stats(cpm)[expressed]
    stats.to_csv(SRC / "fig_main_stats.csv")
    # Fos is reported in Extended Data, so its per-animal values travel too.
    cpm.loc[ANIMALS, [g for g in ("Fos",) if g in cpm.columns]].to_csv(
        SRC / "fig_main_fos_per_animal.csv")

    # Abundance is reported in Extended Data; its per-animal values travel too.
    pct = {a: (is_pop & (animals == a)).sum() / (animals == a).sum() * 100
           for a in ANIMALS}
    pd.DataFrame({"abundance": pct}).loc[ANIMALS].to_csv(SRC / "fig_main_abundance.csv")

    # Panels c, d and e are derived here rather than read from a stored table,
    # so the whole figure regenerates from the data it claims to show.
    enr = regional_enrichment(is_pop, nuc)
    enr.to_csv(SRC / "population_regional_enrichment.csv", index=False)
    by_ct = galr1_by_celltype(counts, var, cell_types, animals)
    by_ct.to_csv(SRC / "galr1_by_celltype.csv", index=False)
    markers = population_markers(counts, var, is_pop)
    markers.sort_values(ascending=False).to_csv(SRC / "population_markers.csv")

    # --- representative section: the most intact tissue, whole section --------
    whole = sc.read_h5ad(PROC / "mbh_anatomical.h5ad", backed="r")
    wsec = whole.obs["section"].astype(str).to_numpy()
    wml = whole.obs["ml"].to_numpy()
    balance = {}
    for s in ONE_PER_MOUSE:
        x = wml[wsec == s]
        left, right = (x < -200).sum(), (x > 200).sum()
        balance[s] = min(left, right) / max(left, right)
    rep = max(balance, key=balance.get)
    whole.file.close()

    fig = plt.figure(figsize=(FULL, 5.1))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, .95], hspace=.62)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.55, 1.05, 1.15], wspace=.55)
    bot = outer[1].subgridspec(1, 4, width_ratios=[1.0, 1.25, .85, .85], wspace=.75)

    # (a) one section: the tissue, the receptor, and the population that moved
    ax = fig.add_subplot(top[0]); panel(ax, "a", dx=-0.04, dy=1.13)
    on = sections == rep
    # Subregions as a faint ground, largest first so no region is buried by a
    # bigger one drawn after it.
    present = [k for k in NUCLEUS_LABEL if (on & (nuc == k)).any()]
    ax.scatter(ml[on], dv[on], s=.9, c=TISSUE, linewidths=0, rasterized=True)
    for key in sorted(present, key=lambda k: int((on & (nuc == k)).sum()), reverse=True):
        m = on & (nuc == key)
        ax.scatter(ml[m], dv[m], s=.9, c=NUCLEUS_COLOUR[key], alpha=.30,
                   linewidths=0, rasterized=True)
    # The receptor itself, in two tints: one transcript, or more than one.
    for lo, hi, colour, size in ((1, 1, GALR1_LOW, 1.9), (2, 10 ** 6, GALR1, 3.1)):
        m = on & (galr1 >= lo) & (galr1 <= hi)
        ax.scatter(ml[m], dv[m], s=size, c=colour, linewidths=0, rasterized=True)
    m = on & is_pop
    ax.scatter(ml[m], dv[m], s=15, facecolors="none", edgecolors=POP,
               linewidths=.55, rasterized=True)
    for key in present:
        m = on & (nuc == key)
        if m.sum() < 40:
            continue
        side = ml[m] > 0 if (ml[m] > 0).sum() > 30 else ml[m] < 0
        ax.annotate(NUCLEUS_LABEL[key],
                    (np.median(ml[m][side]), np.median(dv[m][side])),
                    fontsize=5.4, color=INK, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none",
                              alpha=.70))
    ax.set_xlim(-1500, 1500); ax.set_ylim(-100, 1800)
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 500, "500 µm")
    for j, (text, colour) in enumerate((
            ("$\\it{Galr1}$ 1 transcript", GALR1_LOW),
            ("$\\it{Galr1}$ 2+", GALR1),
            (f"{POPNAME}", POP))):
        ax.annotate(text, xy=(-1470, 1760 - j * 105), fontsize=5.4, color=colour,
                    va="top", ha="left", fontweight="bold")
    ax.set_title(f"one section ({rep})", loc="left", pad=2)

    # (b) which cell types carry Galr1 at all
    ax = fig.add_subplot(top[1]); panel(ax, "b", dx=-0.46, dy=1.13)
    land = by_ct.sort_values("pct_pos", ascending=False).head(LANDSCAPE_SHOWN)
    land = land.sort_values("pct_pos")
    ys = np.arange(len(land))
    colours = [POP if n == POPNAME else INK if n == GALR1_RICH else "#C4C4C4"
               for n in land.cell_type]
    ax.barh(ys, land.pct_pos, height=.68, color=colours)
    ax.set_yticks(ys); ax.set_yticklabels(land.cell_type, fontsize=5.2)
    for tick, n in zip(ax.get_yticklabels(), land.cell_type):
        tick.set_fontweight("bold" if n in (POPNAME, GALR1_RICH) else "normal")
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("% of cells $\\it{Galr1}^+$")
    ax.set_title("which cells carry it", loc="left", pad=2)

    # (c) the screen: which of them change with age
    ax = fig.add_subplot(top[2]); panel(ax, "c", dx=-0.30, dy=1.13)
    scr = by_ct.sort_values("lfc")
    ys = np.arange(len(scr))
    passes = ((scr.blocks == 4) & (scr.margin > 0)).to_numpy()
    ax.axvline(0, color=INK, lw=.5)
    ax.scatter(scr.lfc[~passes], ys[~passes], s=10, color="#B0B0B0",
               linewidths=0, zorder=3)
    # Whatever passes the criterion is drawn and named.  Marking every passing
    # type "this population" would be true only while exactly one passes, and
    # the panel's claim is that the criterion was applied blind.
    is_pop_row = (scr.cell_type == POPNAME).to_numpy()
    for mask, colour in ((passes & is_pop_row, POP), (passes & ~is_pop_row, INK)):
        if not mask.any():
            continue
        ax.scatter(scr.lfc[mask], ys[mask], s=24, color=colour, linewidths=0,
                   zorder=4)
        for y, row in zip(ys[mask], scr[mask].itertuples()):
            label = (POPNAME if row.cell_type == POPNAME
                     else f"{row.cell_type} ({row.pct_pos:.0f}% pos.)")
            ax.annotate(label, (row.lfc, y), xytext=(5, 0),
                        textcoords="offset points", fontsize=5.2, va="center",
                        ha="left", color=colour, fontweight="bold")
    # The population that carries the most Galr1 is named whether it moves or
    # not: if ageing simply raised Galr1, this is where it would show.
    rich = scr[scr.cell_type == GALR1_RICH]
    if len(rich):
        y = int(ys[(scr.cell_type == GALR1_RICH).to_numpy()][0])
        ax.scatter(rich.lfc, [y], s=16, color=INK, linewidths=0, zorder=4)
        ax.annotate(f"{GALR1_RICH}\n({float(rich.pct_pos.iloc[0]):.0f}% pos., unchanged)",
                    (float(rich.lfc.iloc[0]), y), xytext=(-5, 0),
                    textcoords="offset points", fontsize=5.2, va="center",
                    ha="right", color=INK)
    ax.set_yticks([]); ax.set_ylim(-1, len(scr))
    # Limits follow the data: a fixed window silently drops any cell type that
    # moves further than the window was drawn for.
    ax.set_xlim(min(-0.9, float(scr.lfc.min()) - .15),
                max(1.9, float(scr.lfc.max()) + .15))
    ax.set_xlabel("$\\it{Galr1}$, aged / adult (log$_2$)")
    ax.set_ylabel(f"{len(scr)} cell types")
    ax.set_title(f"{int(passes.sum())} of {len(scr)} change", loc="left", pad=2)

    # (d) where the population that moved sits
    ax = fig.add_subplot(bot[0]); panel(ax, "d", dx=-0.44)
    reg = enr.sort_values("enrichment")
    ys = np.arange(len(reg))
    ax.barh(ys, reg.enrichment, height=.68,
            color=[POP if n == "DMH" else "#C4C4C4" for n in reg.nucleus])
    ax.axvline(1, color=INK, lw=.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([NUCLEUS_LABEL.get(n, n) for n in reg.nucleus])
    for tick, n in zip(ax.get_yticklabels(), reg.nucleus):
        tick.set_fontweight("bold" if n == "DMH" else "normal")
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("enrichment over the whole window")
    ax.set_title("most enriched in DMH", loc="left", pad=2)

    # (e) what it is
    ax = fig.add_subplot(bot[1]); panel(ax, "e", dx=-0.34)
    top_marks = markers.nlargest(MARKERS_SHOWN).sort_values()
    ys = np.arange(len(top_marks))
    ax.barh(ys, top_marks.values, height=.68,
            color=[GALR1 if g in RECEPTORS else "#C4C4C4" for g in top_marks.index])
    ax.set_yticks(ys)
    ax.set_yticklabels(top_marks.index, fontsize=5.6, style="italic")
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("log$_2$ vs rest of window")
    hm = hypomap_match()
    if hm is not None:
        ax.annotate(f"HypoMap C185: {hm[0]}\nSpearman ρ = {hm[1]:.2f}",
                    xy=(.97, .06), xycoords="axes fraction", fontsize=5.4,
                    color=INK, ha="right", va="bottom")
    ax.set_title("what it is", loc="left", pad=2)

    # (f, g) the receptors in that population
    for j, gene in enumerate(MAIN_GENES):
        ax = fig.add_subplot(bot[j + 2]); panel(ax, "fg"[j], dx=-0.42)
        r = stats.loc[gene]
        paired_panel(ax, cpm[gene].to_dict(), f"$\\it{{{gene}}}$ (log$_2$ CPM)",
                     f"{2 ** r.lfc:.2f}×  $P$ = {r.p:.3f}")

    fig.savefig(OUT / "figure_main.pdf")
    fig.savefig(OUT / "figure_main.png")
    plt.close(fig)

    print("whole-section bilateral balance (1.0 = symmetric):")
    for s, v in sorted(balance.items(), key=lambda kv: -kv[1]):
        print(f"   {s:8s} {v:.3f}{'   <- shown' if s == rep else ''}")
    print("\nsubregions in the window (pooled over the eight sections):")
    in_window = [k for k in NUCLEUS_LABEL if (nuc == k).any()]
    for key in sorted(in_window, key=lambda k: int((nuc == k).sum()), reverse=True):
        m = nuc == key
        a_ml, d = np.abs(ml[m]), dv[m]
        print(f"   {key:7s} {int(m.sum()):6d} cells   "
              f"|ml| {np.percentile(a_ml, 5):4.0f}-{np.percentile(a_ml, 95):4.0f}"
              f"   dv {np.percentile(d, 5):4.0f}-{np.percentile(d, 95):4.0f} um")

    print(f"\n{int(is_pop.sum())} {POPNAME} neurons across 8 animals, "
          f"{float(by_ct.loc[by_ct.cell_type == POPNAME, 'pct_pos'].iloc[0]):.0f}% "
          "of them Galr1+")
    for g in KEY_GENES:
        if g in stats.index:
            r = stats.loc[g]
            print(f"  {g:6s} {2 ** r.lfc:.2f}x  blocks {int(r.blocks)}/4  P {r.p:.4f}")

    print(f"\nGalr1 screen over {len(by_ct)} cell types "
          f"(>= {MIN_CELLS_PER_ANIMAL} cells in every animal):")
    shown = by_ct[["cell_type", "n_cells", "pct_pos", "fold", "blocks", "p", "margin"]]
    with pd.option_context("display.width", 200):
        print(shown.to_string(index=False))
    print(f"\nwrote {OUT / 'figure_main.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

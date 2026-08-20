#!/usr/bin/env python3
"""Extended Data: the technical validation behind the main figure.

Nothing here is a biological result.  Each panel answers one objection a
reviewer would reasonably raise about the main figure.

  a  design -- eight animals in four blocks, one aged and one adult each
  b  the two cohorts ran different gene panels; how much they share
  c  anatomical registration succeeded in every section
  d  the age effect does not depend on which section was chosen per mouse
  e  how often the significance criterion fires across the whole panel
  f  the conclusion is unchanged when segmentation is re-derived from raw
     transcripts rather than taken from the vendor
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.figstyle import (  # noqa: E402
    ADULT, AGED as C_AGED, FULL, INK, MUTED, POP, panel, use_style,
)
from spatial_lea.io import (  # noqa: E402
    ANIMAL_META, BLOCKS, ONE_PER_MOUSE, RAW, SECTION_ANIMAL,
)

OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
GREY = "#9A9A9A"


def panel_genes(section: str) -> set:
    payload = json.loads((RAW / section / "gene_panel.json").read_text())["payload"]
    return {t["type"]["data"]["name"] for t in payload["targets"]
            if t.get("type", {}).get("data", {}).get("name")}


def main() -> int:
    use_style()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(FULL, 4.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=.75, wspace=.48)

    # (a) design
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a")
    for i, (b, (a_aged, a_adult)) in enumerate(BLOCKS.items()):
        y = len(BLOCKS) - 1 - i
        ax.plot([0, 1], [y, y], color="#D0D0D0", lw=.6, zorder=1)
        ax.scatter([0], [y], s=22, color=C_AGED, zorder=3, linewidths=0)
        ax.scatter([1], [y], s=22, color=ADULT, zorder=3, linewidths=0)
        # Names sit beside their own marker, not above the row, so a label can
        # never be read as belonging to the block below it.
        ax.text(-0.10, y, f"{b}   {a_aged}, {ANIMAL_META[a_aged]['age_weeks']} wk",
                ha="right", va="center", fontsize=5, color=INK)
        ax.text(1.10, y, f"{a_adult}, {ANIMAL_META[a_adult]['age_weeks']} wk",
                ha="left", va="center", fontsize=5, color=INK)
    ax.set_xlim(-1.25, 2.05); ax.set_ylim(-.5, len(BLOCKS) - .5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["aged", "adult"])
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_title("one aged and one adult mouse per block,\n"
                 "matched for slide run, panel and chemistry",
                 loc="left", pad=4, color=MUTED, fontsize=5.8)

    # (b) panel overlap
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b")
    g1, g2 = panel_genes("M493_1"), panel_genes("K238_2")
    shared, only1, only2 = len(g1 & g2), len(g1 - g2), len(g2 - g1)
    ax.barh([0], [shared], color=GREY, height=.55)
    ax.barh([0], [only1], left=[shared], color="#D8D8D8", height=.55)
    ax.barh([1], [shared], color=GREY, height=.55)
    ax.barh([1], [only2], left=[shared], color="#D8D8D8", height=.55)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["cohort 1", "cohort 2"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("panel targets")
    ax.annotate(f"{shared} shared targets", xy=(shared / 2, 1.42), ha="center",
                fontsize=5.5, color=INK)
    ax.annotate(f"{only1}", xy=(shared + only1 / 2, 0), xytext=(0, -13),
                textcoords="offset points", ha="center", fontsize=5, color=MUTED)
    ax.annotate(f"{only2}", xy=(shared + only2 / 2, 1), xytext=(0, 11),
                textcoords="offset points", ha="center", fontsize=5, color=MUTED)
    ax.annotate("cohort-specific", xy=(shared + 14, .5), ha="left", va="center",
                fontsize=5, color=MUTED)
    ax.set_xlim(0, shared + 100); ax.set_ylim(-.55, 1.75)
    ax.set_title("analysis uses the shared set only", loc="left", pad=4,
                 color=MUTED, fontsize=5.8)

    # (c) registration QC
    ax = fig.add_subplot(gs[0, 2]); panel(ax, "c")
    qc = pd.read_csv(REPO / "results" / "anatomy" / "frame_qc.csv", index_col=0)
    qc = qc.loc[[s for s in SECTION_ANIMAL if s in qc.index]]
    used = [s in ONE_PER_MOUSE for s in qc.index]
    ax.scatter(np.arange(len(qc))[~np.array(used)],
               qc.pct_arc_in_window[~np.array(used)], s=11, facecolors="none",
               edgecolors=MUTED, linewidths=.5)
    ax.scatter(np.arange(len(qc))[np.array(used)],
               qc.pct_arc_in_window[np.array(used)], s=13, color=INK, linewidths=0)
    ax.axhline(55, color=MUTED, lw=.5, ls=(0, (3, 2)))
    ax.annotate("acceptance threshold", xy=(0.3, 57), fontsize=5.2, color=MUTED)
    ax.set_ylim(0, 105); ax.set_xticks([])
    ax.set_xlabel(f"{len(qc)} sections")
    ax.set_ylabel("ARC anchors inside\nfitted window (%)")
    ax.set_title("filled, sections analysed", loc="left", pad=4, color=MUTED,
                 fontsize=5.8)

    # (d) section-choice sensitivity
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "d")
    sens = pd.read_csv(REPO / "results" / "gal_n8" / "section_choice_sensitivity.csv")
    ax.hist(sens.mean_lfc, bins=18, color=GREY)
    chosen = sens[sens.combo == "F536_3+G073_2+M493_1+M399_1"]
    if len(chosen):
        ax.axvline(float(chosen.mean_lfc.iloc[0]), color=INK, lw=1.0)
        ax.annotate("set used", xy=(float(chosen.mean_lfc.iloc[0]), ax.get_ylim()[1]),
                    xytext=(2, -2), textcoords="offset points", fontsize=5.5,
                    color=INK, va="top")
    ax.axvline(0, color=MUTED, lw=.5)
    ax.set_xlabel("$\\it{Galr1}$ aged / adult (log$_2$)")
    ax.set_ylabel(f"section sets ({len(sens)})")
    ax.set_title(f"positive in {int((sens.mean_lfc > 0).sum())} of {len(sens)} "
                 "possible\nsection choices", loc="left", pad=4, color=MUTED,
                 fontsize=5.8)

    # (e) panel-wide calibration
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "e")
    cal = pd.read_csv(REPO / "results" / "gal_n8" / "calibration_by_cell_type.csv")
    bars = [("all four blocks\nagree", cal.pct_blocks_agree_4.mean(), 12.5),
            ("groups separate\ncompletely", cal.pct_separating.mean(), 2.9),
            ("both, with\neffect size", cal.pct_strict.mean(), 0.36)]
    xs = np.arange(len(bars))
    ax.bar(xs, [b[1] for b in bars], width=.5, color=GREY)
    for x, (_, obs, exp) in zip(xs, bars):
        ax.plot([x - .3, x + .3], [exp, exp], color=INK, lw=1.0)
    ax.set_xticks(xs); ax.set_xticklabels([b[0] for b in bars], fontsize=5.2)
    ax.set_ylabel("panel genes meeting\ncriterion (%)")
    ax.set_title("bars, observed; rules, expected by chance", loc="left", pad=4,
                 color=MUTED, fontsize=5.8)

    # (f) segmentation validation
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "f")
    mm = pd.concat([pd.read_csv(p) for p in
                    sorted((REPO / "results" / "baysor").glob("*_mixed_matched.csv"))],
                   ignore_index=True)
    xs = np.arange(len(mm)); w = .34
    ax.bar(xs - w / 2, mm.vendor, w, color="#D8D8D8")
    ax.bar(xs + w / 2, mm.baysor_matched, w, color=GREY)
    ax.set_xticks(xs); ax.set_xticklabels(mm.section, fontsize=5.2)
    ax.set_ylabel("cells with conflicting\ntransmitter markers (%)")
    ax.annotate("vendor", xy=(-w / 2, mm.vendor.iloc[0] + .8), ha="center",
                fontsize=5.2, color=MUTED)
    ax.annotate("re-segmented", xy=(w / 2, mm.baysor_matched.iloc[0] + .8),
                ha="center", fontsize=5.2, color=INK)
    ax.set_ylim(0, mm.vendor.max() * 1.35)
    ax.set_title("matched objects at matched depth", loc="left", pad=4,
                 color=MUTED, fontsize=5.8)

    fig.savefig(OUT / "figure_extended_data.pdf")
    fig.savefig(OUT / "figure_extended_data.png")
    plt.close(fig)
    print(f"panels shared {shared}, cohort-specific {only1}/{only2}")
    print(f"registration: {qc.pct_arc_in_window.min():.0f}-{qc.pct_arc_in_window.max():.0f}% "
          f"across {len(qc)} sections")
    print(f"wrote {OUT / 'figure_extended_data.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

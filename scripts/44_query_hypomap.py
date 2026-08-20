#!/usr/bin/env python3
"""Match Xenium cell types to HypoMap clusters.  Runs on your machine, not here.

    python scripts/44_query_hypomap.py \
        --hypomap ~/Desktop/hypoMap.h5ad \
        --profiles results/hypomap/xenium_celltype_profiles.csv \
        --population "Glut Prdm8/Cbln1"

What it does, and what it deliberately does not do
--------------------------------------------------
The Xenium panel measures 287 genes; HypoMap is whole transcriptome.  Anything
that embeds the two together -- ingest, scArches, symphony -- has to reconstruct
20,000 missing genes per Xenium cell, and the reconstruction is what then drives
the assignment.  On a 287-gene panel that is mostly imputation, so this script
does not do it.

Instead it compares pseudobulk profiles directly on the shared genes: each
Xenium cell type against each HypoMap cluster, by Spearman correlation of
log-CPM.  Rank correlation rather than Pearson, because the two platforms have
different detection efficiencies per gene, which shifts magnitudes but largely
preserves order.

That gives a ranked shortlist and a margin between the best and second-best
cluster.  It is a defensible correspondence, not a label transfer, and the
output says so.  A correspondence with a small margin means the panel cannot
separate those clusters -- which is a fact about the panel, not a failure of
the match.

Three things are reported:
  1. best-matching HypoMap clusters, at every annotation level HypoMap carries
  2. how the population's own defining markers behave in those clusters, so the
     match can be checked by eye rather than trusted
  3. Galr1 and Gal expression across HypoMap, which says whether the receptor
     result sits in a cluster HypoMap also considers galanoceptive

Requires: scanpy, pandas, numpy, scipy.  HypoMap is large; the cluster means are
accumulated in chunks so it never needs to be held in memory uncompressed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import spearmanr

CHUNK = 20000
MARKERS = ["Cbln1", "Slc17a6", "Otp", "Prdm8", "Bdnf", "Galr1", "Ghsr",
           "Galr3", "Gal", "Ghrh"]


def resolve_symbols(hm, panel: list[str]) -> pd.Index:
    """Find where the gene symbols live and return them aligned to hm.var.

    HypoMap distributions vary: the index may be symbols, Ensembl ids, or an
    internal id, with symbols in a var column under any of several names.  Rather
    than guess, every candidate is scored by how many panel genes it recovers and
    the best is taken; a case-insensitive pass catches objects that upper-case
    symbols.
    """
    want = {g.lower() for g in panel}
    best, best_hits, best_where = None, 0, ""
    candidates = [("var_names", pd.Index(hm.var_names.astype(str)))]
    for col in hm.var.columns:
        values = hm.var[col]
        if values.dtype.name in ("category", "object"):
            candidates.append((f"var['{col}']", pd.Index(values.astype(str))))
    for where, values in candidates:
        hits = sum(v.lower() in want for v in values)
        if hits > best_hits:
            best, best_hits, best_where = values, hits, where
    if best is None or best_hits == 0:
        preview = ", ".join(map(str, hm.var_names[:4]))
        cols = ", ".join(hm.var.columns) or "(none)"
        raise SystemExit(
            "Could not find gene symbols in this HypoMap object.\n"
            f"  var_names look like: {preview}\n"
            f"  var columns: {cols}\n"
            "Open it and check which column holds symbols such as Agrp or Pomc.")
    print(f"  gene symbols taken from {best_where} ({best_hits} panel genes matched)")
    return pd.Index(best)


def annotation_levels(obs: pd.DataFrame) -> list[str]:
    """HypoMap's nested cluster columns, coarse to fine.

    Named C<n>_named in the published object.  Anything else that looks like a
    categorical annotation with a sensible number of levels is offered too, so
    the script still works if the object has been renamed.
    """
    named = [c for c in obs.columns
             if c.startswith("C") and c.endswith("_named")]
    if named:
        return sorted(named, key=lambda c: int(c.split("_")[0][1:]))
    return [c for c in obs.columns
            if obs[c].dtype.name in ("category", "object")
            and 3 <= obs[c].nunique() <= 600]


def cluster_means(adata, key: str, genes: list[str],
                  symbols: pd.Index) -> pd.DataFrame:
    """Mean log-CPM per cluster, accumulated in chunks."""
    lookup = {}
    for i, s in enumerate(symbols):
        lookup.setdefault(str(s).lower(), i)
    idx = [lookup[g.lower()] for g in genes]
    labels = adata.obs[key].astype(str).to_numpy()
    groups = pd.unique(labels)
    total = pd.DataFrame(0.0, index=groups, columns=genes)
    n = pd.Series(0, index=groups, dtype=int)

    for start in range(0, adata.n_obs, CHUNK):
        stop = min(start + CHUNK, adata.n_obs)
        block = adata[start:stop]
        x = block.X
        x = x.toarray() if issparse(x) else np.asarray(x)
        # HypoMap ships log-normalised counts; if this object holds raw counts,
        # normalise it the same way rather than mixing scales.
        if x.max() > 50:
            libsize = np.maximum(x.sum(axis=1, keepdims=True), 1)
            x = np.log1p(x / libsize * 1e4)
        x = x[:, idx]
        lab = labels[start:stop]
        frame = pd.DataFrame(x, columns=genes)
        frame["_g"] = lab
        agg = frame.groupby("_g").sum()
        cnt = frame.groupby("_g").size()
        total.loc[agg.index] += agg
        n.loc[cnt.index] += cnt
    return total.div(n.replace(0, np.nan), axis=0).dropna(how="all")


def match(profile: pd.Series, means: pd.DataFrame) -> pd.DataFrame:
    rho = {}
    for cluster, row in means.iterrows():
        ok = row.notna() & profile.notna()
        if ok.sum() >= 20:
            rho[cluster] = spearmanr(profile[ok], row[ok]).statistic
    out = pd.Series(rho).sort_values(ascending=False).to_frame("spearman")
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hypomap", type=Path,
                   default=Path.home() / "Desktop" / "hypoMap.h5ad")
    p.add_argument("--profiles", type=Path, default=None,
                   help="defaults to xenium_celltype_profiles.csv beside this "
                        "script, or in results/hypomap/")
    p.add_argument("--population", default="Glut Prdm8/Cbln1")
    p.add_argument("--top", type=int, default=8)
    p.add_argument("--out", type=Path, default=Path("results/hypomap"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # Look beside the script and in the repo layout before giving up, so the
    # script works whether it was copied to a desktop or run from the repo.
    name = "xenium_celltype_profiles.csv"
    here = Path(__file__).resolve().parent
    candidates = ([args.profiles] if args.profiles else
                  [here / name, here.parent / "results" / "hypomap" / name,
                   Path.cwd() / name, Path("results/hypomap") / name])
    found = next((c for c in candidates if c and c.exists()), None)
    if found is None:
        raise SystemExit(
            f"Could not find {name}. Looked in:\n  "
            + "\n  ".join(str(c) for c in candidates)
            + "\nPass it explicitly with --profiles /path/to/" + name)
    print(f"Profiles: {found}")
    prof = pd.read_csv(found, index_col=0)
    if args.population not in prof.index:
        raise SystemExit(f"{args.population!r} not in {args.profiles}. "
                         f"Available: {', '.join(prof.index[:10])} ...")

    print(f"Reading {args.hypomap} (backed, chunked)")
    hm = sc.read_h5ad(args.hypomap, backed="r")
    print(f"  {hm.n_obs:,} cells x {hm.n_vars:,} genes")

    symbols = resolve_symbols(hm, list(prof.columns))
    have = {str(s).lower() for s in symbols}
    shared = [g for g in prof.columns if g.lower() in have]
    print(f"  {len(shared)} of {prof.shape[1]} panel genes present in HypoMap")
    missing = [g for g in prof.columns if g.lower() not in have]
    if missing:
        print(f"  absent: {', '.join(missing[:12])}"
              + (" ..." if len(missing) > 12 else ""))
    if len(shared) < 50:
        raise SystemExit("Too few shared genes; check that gene symbols match "
                         "(HypoMap uses mouse symbols such as Agrp, Pomc).")

    levels = annotation_levels(hm.obs)
    print(f"  annotation levels found: {', '.join(levels)}\n")

    target = prof.loc[args.population, shared]
    all_hits = []
    for key in levels:
        means = cluster_means(hm, key, shared, symbols)
        hits = match(target, means)
        margin = (hits["spearman"].iloc[0] - hits["spearman"].iloc[1]
                  if len(hits) > 1 else np.nan)
        print(f"=== {key}: {len(means)} clusters ===")
        print(hits.head(args.top).round(3).to_string())
        print(f"  margin over runner-up: {margin:.3f}"
              + ("   (small -- the panel cannot separate these)"
                 if margin < .02 else ""))
        best = hits.index[0]
        show = [g for g in MARKERS if g in means.columns]
        print(f"  markers in {best}:")
        print("   " + "  ".join(f"{g} {means.loc[best, g]:.2f}" for g in show))
        print()
        hits["level"] = key
        all_hits.append(hits.assign(cluster=hits.index))

        # Every Xenium cell type against this level, for the record.
        table = pd.DataFrame({t: match(prof.loc[t, shared], means)["spearman"]
                              for t in prof.index})
        table.to_csv(args.out / f"hypomap_{key}_all_celltypes.csv")

    pd.concat(all_hits).to_csv(args.out / "hypomap_best_matches.csv", index=False)

    # Where does HypoMap put Galr1 and Gal?  Reported at C185, the level that
    # is interpretable as neuron types; the finer levels split the same
    # populations further without adding an identity a reader can use.
    finest = next((l for l in levels if l.startswith("C185")), levels[-1])
    means = cluster_means(hm, finest,
                          [g for g in ("Galr1", "Gal", "Galr3")
                           if g.lower() in have], symbols)
    print(f"=== Galr1 and Gal across HypoMap {finest} ===")
    ranked = means.sort_values("Galr1", ascending=False)
    print("  highest Galr1:")
    print(ranked.head(10).round(3).to_string())
    if "Gal" in means.columns:
        print("\n  highest Gal, for contrast -- the ligand and the receptor are "
              "carried by different clusters:")
        print(means.sort_values("Gal", ascending=False).head(10).round(3).to_string())
    means.to_csv(args.out / f"hypomap_{finest}_galanin.csv")

    print(f"\nWrote tables to {args.out}")
    print("These are profile correspondences on a 287-gene panel, not label "
          "transfer.\nQuote the correlation and the margin, not the cluster name "
          "alone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

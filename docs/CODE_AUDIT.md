# Code audit

Six independent lenses over the pipeline, each finding adversarially verified,
plus first-hand verification of the highest-stakes claims by re-running the
analysis. Everything below was checked against `data/processed/*.h5ad` on this
branch. Findings are ranked by what threatens the result, not by how loudly
they fail.

Verification status is marked per finding: **[verified here]** means the number
was reproduced by running the pipeline during this audit; **[reported]** means
an audit agent found it and it survived adversarial review but was not
re-derived independently.

---

## 1. The headline depends on one undocumented section substitution — HIGH

**[verified here]** `src/spatial_lea/io.py:51-59`, consumed at
`scripts/42_figure_main.py`.

Two "one section per mouse" sets exist, both commented as lab-chosen, and they
disagree on three of four cohort-A animals:

    SELECTED_SECTIONS = ("F536_1", "G073_2", "M399_3", "M493_2", ...)
    ONE_PER_MOUSE     = ("F536_3", "G073_2", "M399_1", "M493_1", ...)

Only `SELECTED_SECTIONS` has an artifact behind it: `config/samples.yaml` marks
exactly those four with `in_roi_export: true`. `ONE_PER_MOUSE` — the set the
figure uses — matches no rule recorded anywhere in the repository.

Re-running the screen on each set:

| sections | Galr1 in Glut Prdm8/Cbln1 | P | margin | types significant |
|---|---|---|---|---|
| `ONE_PER_MOUSE` (used) | 1.46x | 0.0286 | +0.112 | 1 of 12 |
| `SELECTED_SECTIONS` | 1.34x | 0.1143 | -0.074 | **0 of 13** |
| `ONE_PER_MOUSE`, M493_2 for M493_1 | 1.41x | 0.0857 | -0.074 | **0 of 13** |

One substitution, M493_1 in place of M493_2, is necessary and sufficient for the
result. Nothing in the repository states that, and the Extended Data panel that
would disclose it (`scripts/43_figure_extended.py:126`) reads
`results/gal_n8/section_choice_sensitivity.csv`, a file no committed script
writes — script 43 dies there with `FileNotFoundError`.

**Fix.** Derive the set from a single named criterion already recorded (e.g.
argmax of `quality` in `results/anatomy/frame_qc.csv`) instead of asserting a
literal, and write the section-choice sweep that panel already expects.

## 2. Half the design carries group labels with no source — HIGH

**[verified here]** `src/spatial_lea/io.py:77-84` vs `config/samples.yaml:42-45`.

The manifest records cohort B as `age_weeks: null, age_group: null` with the
comment "ages and groups not yet supplied. Confirm before use". `io.py` asserts
K238 aged/73w, P953 adult/25w, F739 aged/68w, Q378 adult/24w. No source for
those values exists in the repository, the raw bundles or the Drive manifest.

The labels are load-bearing: swapping cohort B's assignment turns Galr1 from
+0.542 (1.46x) to -0.087. Neither cohort alone can produce a significant result
— a 2-vs-2 permutation has a floor of P = 0.333, and both cohorts sit exactly on
it. The result exists only in the 4-vs-4 combination, and that combination
exists only because of these labels.

`config/samples.yaml` also does not parse (`yaml.safe_load` raises at line 42 on
the bare `?` in `sex: female?`), and no code reads it, so the divergence between
the manifest and `io.py` cannot be detected by anything.

**Fix.** Get the ages from the lab. Then quote the `female?` values so the file
parses, and replace the `ANIMAL_META` literal with a loader that reads the
manifest and raises on a null `age_group`.

## 3. The significance rule is one condition stated as three — HIGH

**[verified here]** `scripts/42_figure_main.py`, `blocked_stats` and the
`significant` column.

The rule reads as three hurdles: all four run-pairs agreeing in sign, the two
groups not overlapping, and a permutation P <= 0.05. Over 7,077 gene x cell-type
rows: all 260 rows with non-overlapping groups have P exactly 2/70 = 0.0286 and
all four signs agreeing, and the lowest P among the other 6,817 rows is 0.0571.
The full rule and `margin > 0` alone select identical sets.

So "P = 0.029" carries no information beyond "the four aged values are all above
the four adult values". With eight animals and an exact permutation test, 2/70
is the smallest attainable P, and complete separation forces it.

**Fix.** State the criterion as what it is — complete separation of 4 vs 4 —
and report P as the exact permutation value at its floor, or drop two of the
three conjuncts from the description.

## 4. Multiplicity is not accounted for, and it changes the headline — HIGH

**[verified here]** Null distributions built from the 36 balanced regroupings
(2 aged + 2 adult per side, carrying no age contrast):

- **Cell-type screen.** The real split gives 1 of 12 types passing. Six of 36
  null splits also give at least one. **P ~ 0.19.** Finding one cell type out of
  twelve is close to what the criterion does by chance.
- **Gene list in that population.** The real split gives 22 of 169 genes; the
  null splits give a median of 3 and never reach 22. **P = 0.027**, with roughly
  14% of the 22 expected to be false.

The population-level transcriptional change is supported; "Galr1 changes in
exactly one cell type" is marginal once the twelve tests are counted.

**Fix.** Report the screen with its family-wise error, and lead on the
population-level result, which is the one that survives.

## 5. `docs/RESULTS.md` states the negation of the headline — HIGH

**[verified here]** `:28` "Galr1 shows no age difference that survives the
additional sections", `:164` "No `Galr1` population survives. Not one, in any of
the four carvings." Both predate the eight-animal design. The word *Prdm8*
appears nowhere in `docs/` outside `M617_IMPLICATIONS.md`.

**[reported]** Its anatomy tables no longer reproduce: `:298`/`:314` give the DMH
as 1,254 cells/section and 1000 x 477 um; the current run gives 13,042 cells and
903 x 431 um, and all seven nuclei disagree. Six of the eight KMeans domain ids
cited in prose now denote different tissue.

**Fix.** Banner the file as cohort A, n = 2/group, add a current section, and
regenerate the anatomy tables from `results/extended_nuclei/nucleus_extents.csv`
rather than by hand. Stop citing bare KMeans ids in prose — they are not stable
across runs.

## 6. The domain-8 split is documented, plotted and reported, but never runs — HIGH

**[verified here]** `scripts/18_extended_nuclei.py:39` defines `SPLIT_ML_UM`;
it is referenced only inside `_plot` (`:216`, `:218`). `assign()` maps one label
per domain with no `ml` test, so `nucleus_ext` never contains the fringe. The
docstring, `docs/RESULTS.md:289-301` and a shipped figure panel all describe the
split as applied.

**[reported]** The hardcoded id is also now the wrong domain — domain 8 in this
run is 62% VLMC/fibroblast — which is the exact failure the file's own header
says it removed.

**[verified here]** `:137` computes `per_section` as `len(sub) / 12` while the
object holds 16 sections, so every published per-section count is 33% high. The
same literal appears at `scripts/22_m617_target_atlas.py:107` and `:151`.

**Fix.** Delete the split or implement it inside `assign()` selecting the domain
by evidence rather than by the literal `"8"`. Pass `obs["section"].nunique()`
into `extent()` instead of the 12.

## 7. `docs/M617_IMPLICATIONS.md` names the wrong nucleus — HIGH

**[reported]** The doc gives DMH 37.2% of all Galr1 against LHA 23.0%.
Recomputed on the same cohort-A sections: LHA 39.7%, DMH 21.8% — inverted. The
top Galr1-carrying population becomes `LHA | GABA Cacna2d2`, which the doc's
table does not contain. This document exists to direct a wet-lab experiment.

**Fix.** Re-run `22_m617_target_atlas.py` with the cohort filter it describes and
the section-count fix, rewrite the section-1 table, and commit its outputs so a
reader can check them.

## 8. Figure fidelity — MEDIUM

**[verified here]** Panel e's UMAP is the reference clustering: 37,483 cells
from the four ROI sections, not the 115,972 cells from eight sections the rest
of the figure uses. The panel no longer says so.

**[verified here]** Panel b draws one section's tissue with all eight animals'
population cells on it. Both facts belong in the legend.

**[reported]** Panel f's axis starts at 0.6, which exaggerates the gap between
the top HypoMap match and its rivals.

## 9. Smaller findings — MEDIUM to LOW

**[reported]**

- `transfer_labels` in `13_spatial_domains.py` names a cluster by bare plurality
  of reference labels, with no margin requirement.
- The domain cache is reused on file existence alone, so an upstream change does
  not invalidate it.
- `NUCLEUS_TYPES` in `anatomy.py` names cell types that no longer exist; the fit
  silently ignores them.
- `docs/DATA_AUDIT.md` describes 14 sections; 16 exist.
- `_robust_gaussian` returns a trimmed covariance with no consistency
  correction, so fitted nuclei are slightly too small.
- `docs/STATISTICS.md` is stale and its power table hardcodes 3 sections.

## What the audit checked and found sound

- The anatomical frame: ventricle tracing, ARC-seeded ventral direction, the
  rank-based anchor rule that fixed the second cohort's inverted frames.
- `04_annotate.py`'s marker-based cluster matching, including the case of more
  clusters than table entries.
- Cell-level plumbing: no step silently drops or duplicates cells; obs and
  count-matrix ordering hold throughout.
- Panel c's row labels, bold flags and colours track the same ordering.
- The reconstruction of M493_2's missing `cells.parquet`, now superseded by the
  real file, which it matched to 0.04 um on centroids and exactly on counts.
- The figure regenerates byte-identically from the data in one command.

## Fix these three first

1. **The section set (1).** The headline is currently a property of a choice
   nobody recorded. Either justify `ONE_PER_MOUSE` against a stated criterion or
   report the sweep — until then the P value is not interpretable.
2. **The cohort-B ages (2).** Four of eight animals carry invented labels, and
   the result does not exist without them.
3. **The reporting of significance (3, 4).** State the criterion as complete
   separation, and report the screen's family-wise error. The population-level
   result (22 genes, P = 0.027) is the defensible headline.

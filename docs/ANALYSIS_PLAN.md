# Analysis plan — Galanin receptor signalling in the ageing mediobasal hypothalamus

**Question.** In which mediobasal hypothalamic (MBH) cell types — neuronal
subtypes in particular — does **Galr1** expression differ between mature adult
(~7 mo) and middle-aged (~15–16 mo) male mice, and is that change restricted to
specific nuclei such as the DMH?

**Data.** 10x Xenium, 297-gene `mBrain_50g` panel, 4 male animals (2 aged, 2
adult) with 2–3 sections each; MBH rectangles supplied for one section per
animal (40 863 cells). A second cohort of 4 animals exists on a different
platform version.

**Status.** The panel contains `Galr1`, `Galr3` and `Gal`, the ROI convention is
solved exactly, and the design is properly blocked. One thing is missing before
any of this can run: **the cell × gene expression matrix is not in the shared
Drive folder** (see §1).

Full factual basis: [`DATA_AUDIT.md`](DATA_AUDIT.md).

---

## 1. Blocker and open decisions

### 1.1 Blocking — expression matrices needed

Each Drive sample folder has only `cells.parquet`, `features.tsv.gz` and
`experiment.xenium`. `cells.parquet` holds *per-cell transcript totals*, not
per-gene counts, so Galr1 cannot be measured and cell types cannot be called.

Please add, for every section, from the Xenium output bundle:

| priority | file | ≈ size / section | needed for |
|---|---|---|---|
| **required** | `cell_feature_matrix.h5` | 50–150 MB | everything |
| high | `nucleus_boundaries.parquet`, `cell_boundaries.parquet` | 10–40 MB | segmentation QC, figures, nucleus-only sensitivity analysis |
| high | `morphology_focus/` (DAPI OME-TIFF) | 0.5–3 GB | Allen CCF registration → nucleus labels; anatomy figure panels |
| medium | `analysis/` (10x clustering, UMAP) | small | independent cross-check of clustering |
| optional | `transcripts.parquet` | 1–3 GB | Baysor/Proseg re-segmentation, subcellular QC |

The required tier is ~1.5 GB for all 14 sections — comfortably within this
environment. `scripts/01_fetch_data.py` will mirror whatever is present, so
dropping the files into the existing Drive folders is enough.

### 1.2 Decisions worth making now

1. **Cohort B metadata.** K238, P953, F739, Q378 (Nov 2025) — sex, ages, groups.
   If that cohort is also adult vs aged it becomes an independent replication
   set, which matters a great deal for §9.
2. **ROIs for the remaining sections.** The export covers 1 section per animal;
   6 more exist (F536 ×2, M493 ×2, G_073 ×1, M399 ×1). Adding them stabilises
   each animal's estimate. Either you draw them, or we propagate the existing
   ROI to sibling sections by landmark registration and you approve the result.
3. **Nucleus-label ground truth.** With DAPI images we can register to the Allen
   CCF. Without them we fall back on marker + geometry based domains (§5), which
   is defensible but weaker for a "the effect is in the DMH" claim.

---

## 2. Design and the statistics that follow from it

| batch | aged | adult |
|---|---|---|
| B1 (26 Jun, Run1) | F536 — 70 w | M493 — 29 w |
| B2 (01 Jul, Run2) | G_073 — 65 w | M399 — 30 w |

A randomised complete block design: **age is orthogonal to batch, run day, slide
and cassette.** Nothing needs to be deconfounded, and each batch yields one
independent aged−adult contrast.

**Unit of replication is the animal, always.** 40 863 cells is a precision
statement about 4 mice, not 40 863 independent observations. Section- and
cell-level data will be aggregated to per-animal values before any test
([pseudoreplication inflates FDR](https://www.sc-best-practices.org/conditions/differential_gene_expression.html);
pseudobulk is the accepted remedy).

### The honest power position

With 2 vs 2 in a blocked design, the exact permutation null has only four
label assignments (swap within B1 × swap within B2). The observed assignment and
its mirror are indistinguishable, so **the smallest attainable exact two-sided
p-value is 0.5, and one-sided 0.25 — for any gene, however large the effect.**

No amount of analysis changes that. So the plan does not pretend otherwise:

- **Cohort A delivers** a rigorously annotated spatial atlas, per-animal effect
  sizes with confidence intervals, replication across two independent batches
  and across sections, and a ranked, pre-registered shortlist of cell type ×
  nucleus hypotheses.
- **A high-impact claim additionally needs** one or more of: cohort B as an
  independent replication; more animals (n ≥ 5/group makes standard pseudobulk
  testing meaningful); orthogonal validation of the top hits by RNAscope/smFISH
  — which such journals expect for an in-situ receptor claim regardless of n.

Formal statistics that *are* legitimate here: moderated tests that borrow
strength across genes (limma-voom / edgeR with `~ age + batch`), reported with
the residual-df limitation stated; sign-consistency across the two blocks;
per-animal dot plots over every bar. Effect sizes lead, p-values follow.

---

## 3. Stage 0 — environment and provenance

Pinned `requirements.txt`, fixed random seeds, one config file
(`config/samples.yaml`) as the single source of truth for group assignment.
Every stage writes a versioned `.h5ad` plus a QC report, so any figure can be
traced to an input file and a commit.

## 4. Stage 1 — QC and cell filtering

Per cell: `total_counts`, genes detected, `cell_area`, `nucleus_area`,
`nucleus_count`, and the control channels already present in `cells.parquet`
(`control_probe_counts`, `genomic_control_counts`, `control_codeword_counts`,
`unassigned_codeword_counts`).

- **Background calibration first.** The 27 negative-control probes and 41
  negative-control codewords give a per-section false-detection rate. This sets
  the threshold for what counts as a real `Galr1` transcript later (§7) — the
  question a reviewer always asks about a 1–2 count GPCR.
- Filters: `total_counts ≥ 10`, `nucleus_count == 1` (drops merged
  segmentations), area within section-wise MAD limits, control fraction below a
  section-specific cut. Every threshold re-run at ±1 step as a sensitivity check.
- Because segmentation is 5 µm nucleus expansion, **transcript bleed-through
  from neighbours is the main artefact risk.** Reported as: correlation of
  neighbouring-cell marker contamination, and a nucleus-only recount if
  boundaries arrive.
- Output: per-section QC panel, retained-cell table, and a go/no-go note per
  section.

## 5. Stage 2 — MBH extraction and an anatomical coordinate frame

ROI selection is implemented and verified (`src/spatial_lea/roi.py`,
`tests/test_roi.py` reproduce all four recorded counts exactly).

The supplied rectangles are generous — they contain ARC, VMH, DMH and ME **plus**
surrounding LHA/PH/TU. Sub-nucleus assignment is what the DMH question rests on,
so it gets three independent routes and a consensus:

1. **Allen CCFv3 registration** (DeepSlice → ABBA/VisuAlign, or `brainreg`)
   using DAPI, or a rendered nucleus-density image if DAPI is unavailable.
   Yields per-cell CCF structure labels (ARH, VMH, DMH, ME, TU, LHA, PH).
2. **Reference mapping to the Allen ABC atlas MERFISH** component (Yao/Zhang
   2023; S3 access verified), which carries both cluster identity and CCF
   coordinates — transfers anatomy and cell type at once.
3. **Unsupervised spatial domains** (Banksy) annotated by marker enrichment plus
   geometry: the 3rd ventricle as a cell-sparse midline, `Gpr50`/`Rfx4`
   tanycytes lining it, ME ventrally, `Adcyap1`-dense ovoid VMH, `Grp`/`Pdyn`
   DMH dorsal to VMH.

Cells are kept where **≥2 routes agree**; the agreement matrix is reported and
boundary-ambiguous cells are excluded, with a dilation/erosion sensitivity run.

**In parallel, a continuous frame** — distance to the 3rd ventricle and
normalised dorsoventral depth, anchored per section. Continuous position is
immune to registration error and lets Galr1 be modelled as a spatial gradient;
discrete nuclei then serve as the interpretable summary. This directly
pre-empts the standard criticism of nucleus-specific claims.

## 6. Stage 3–4 — normalisation, integration, cell typing

Counts normalised per cell and log1p-transformed for embedding; **raw counts
retained for all testing**. With only 297 genes, all genes are used rather than
an HVG subset. Integration over animal with Harmony or scVI, explicitly checked
so that the age signal is not regressed away (integration diagnostics reported
alongside a no-integration run).

Cell types are assigned by three independent routes and reconciled:

- **HypoMap** (Steuernagel 2022, 384 925 cells, hierarchical C7→C465) — the
  reference built for exactly this tissue; CELLxGENE access verified. Label
  transfer restricted to the 297 shared genes via scArches/scANVI.
- **Allen ABC whole-brain taxonomy** for broad classes and cross-checking.
- **Unsupervised Leiden** at several resolutions with manual marker annotation.

Consensus labels plus a per-cell confidence; concordance reported as a matrix.
Expected MBH populations: ARC `Agrp`/`Npy`, ARC `Pomc`/`Cartpt`, ARC `Ghrh`,
ARC `Th`/`Slc6a3` (TIDA), VMH `Adcyap1`/`Esr1`, DMH `Grp`, DMH `Gal`, DMH `Npy`,
`Sst` and `Pdyn` interneuron sets, α/β tanycytes, astrocytes, oligodendrocytes,
OPC, microglia, endothelium, mural, ependyma.

**Sanity check before anything else — does the tissue look older?** `Gfap`,
`Cd68`, `Trem2`, `Spp1`, `Cd44`, `Ly6a`, `Igfbp5` should rise in the aged
animals. If the canonical ageing signature is absent, that is reported first and
it reframes every downstream result. This runs early, not at the end.

## 7. Stage 5 — the Galr1 analysis

"Galr1 changed" is three different biological claims, and they are separated:

| # | claim | test |
|---|---|---|
| **A** | a different *fraction of cells* express Galr1 | per-animal detection rate per cell type × nucleus; background-calibrated threshold from §4; blocked comparison |
| **B** | the *same cells express more/less* Galr1 | pseudobulk per animal × cell type × nucleus, negative-binomial with `log(total counts)` offset, `~ age + batch` |
| **C** | the *abundance of Galr1-high cell types* changed | compositional analysis (scCODA / propeller), which is not the same finding at all |

Reporting a bulk shift without this decomposition is the most common way an
in-situ receptor paper gets sent back. All three are run for **`Galr1`, and in
parallel for `Galr3` and the ligand `Gal`** — a change in the ligand–receptor
*balance* is a stronger and more mechanistic result than a receptor count.

Multiple testing: Benjamini–Hochberg across cell type × nucleus combinations,
with the combination list fixed before testing. Combinations with fewer than a
pre-set minimum cells per animal are excluded up front rather than filtered
after seeing results.

## 8. Stage 6 — spatial and circuit-level analyses

These are what separate a spatial paper from a dissociated-tissue paper:

- **Gal → Galr1 spatial coupling.** Neighbourhood enrichment and co-occurrence
  (squidpy) between `Gal`-expressing and `Galr1`-expressing cells; distance from
  each Galr1⁺ cell to the nearest Gal⁺ cell, adult vs aged. Tests whether the
  ageing change is in receptor abundance or in the *geometry* of the signalling
  relationship.
- **Galr1 as a spatial gradient** along distance-to-3V and dorsoventral depth,
  per cell type, fitted with GAMs and compared between groups.
- **Niche composition** around Galr1⁺ DMH neurons — which cell types are their
  neighbours, and does that change with age.
- **Local ageing microenvironment** — does Galr1 change track local microglial
  or astrocytic activation within the same section.

## 9. Stage 7 — robustness

Run as a fixed battery, reported whether or not it is flattering:

- Leave-one-animal-out for every headline result.
- Section-level replication — each animal's 2–3 sections analysed independently
  (needs the extra ROIs, §1.2).
- Exact blocked permutation, with its 0.25 floor stated rather than hidden.
- Threshold sensitivity: detection cut, QC cuts, nucleus-boundary dilation.
- Segmentation sensitivity: nucleus-only counts, and Baysor/Proseg
  re-segmentation if `transcripts.parquet` is supplied.
- Negative-control genes: the same pipeline on panel genes with no expected age
  effect, plus on negative-control probes, to show specificity.
- Power curve: the effect size detectable at n = 2, 3, 5, 8 per group — this
  doubles as the justification for the next cohort.

## 10. Deliverables

Figure-level plan, built to be dropped into a manuscript:

1. **Design + QC** — cohort table with ages, ROI placement per section, QC
   metrics, negative-control calibration.
2. **MBH cell atlas** — UMAP and spatial maps of consensus cell types, per
   animal; marker dot plot; HypoMap concordance.
3. **Anatomy** — nucleus segmentation with three-route agreement; the continuous
   3V-distance / DV frame; ageing signature confirming the age contrast.
4. **Galr1 map** — spatial expression per section, per cell type, per nucleus;
   detection rate vs level; the A/B/C decomposition with per-animal points.
5. **Nucleus × cell type specificity** — the DMH question answered directly, with
   both blocks shown separately.
6. **Gal–Galr1 circuit** — ligand/receptor spatial coupling, gradients, niches.
7. **Supplement** — robustness battery, sensitivity analyses, power curve, full
   statistics tables.

Plus: versioned `.h5ad` objects, per-cell annotation tables (cell type, nucleus,
coordinates, Galr1 status), a reproducible pipeline, and a written results
report.

---

## 11. Risk register

| risk | severity | handling |
|---|---|---|
| Expression matrix absent from Drive | **blocking** | §1.1 — nothing proceeds until supplied |
| n = 2 vs 2 caps attainable p-values | **high** | §2 — discovery framing, replication in cohort B, orthogonal validation, power curve for the next cohort |
| "Aged" is 15–16 mo = middle-aged by NIA convention | medium | state ages explicitly; frame as mid-life onset; add an 18–24 mo group if possible |
| No `Nr5a1`/`Tbx3`/`Trh` in panel for nucleus ID | medium | three-route consensus + continuous coordinates (§5) |
| Galr1 is a low-abundance GPCR | medium | negative-control-calibrated detection (§4, §7A); report count distributions |
| 5 µm nucleus expansion → transcript bleed-through | medium | nucleus-only recount, contamination diagnostics, optional re-segmentation |
| Cohorts A and B differ in panel design, software and segmentation | medium | never pooled; analysed separately as independent replication |
| ROIs cover only 1 of 2–3 sections per animal | low | request or propagate the remaining ROIs |

---

## 12. Sequence of work

| step | depends on | rough effort |
|---|---|---|
| Fetch matrices, integrity check | files uploaded | < 1 h |
| QC + background calibration | ↑ | 0.5 d |
| ROI extraction *(already implemented and verified)* | ↑ | done |
| Ageing-signature sanity check | QC | 0.5 d |
| Cell typing (3 routes + consensus) | QC | 1–2 d |
| Nucleus assignment (3 routes + consensus) | cell typing, ideally DAPI | 1–2 d |
| Galr1 A/B/C analysis | above | 1 d |
| Spatial / circuit analyses | above | 1 d |
| Robustness battery | above | 1 d |
| Figures + report | above | 1–2 d |

The ageing-signature check is deliberately early: it is the cheapest way to find
out whether the age contrast is detectable at all before investing in the rest.

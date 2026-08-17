# Data audit — Xenium MBH ageing study

Everything below was read directly from the shared Drive folder
(`18VEea5Ath_wbyXeIz_PYZcgwfx6lYJZt`) on 2026-08-17. Nothing here is assumed.

---

## 1. Blocking issue: the expression matrix is not on Drive

Every one of the 14 sample folders contains exactly three files:

| file | size | what it holds |
|---|---|---|
| `cells.parquet` | 1.8–2.6 MB | cell id, centroid, total/control counts, cell & nucleus area, segmentation method |
| `features.tsv.gz` | 4 KB | the panel definition (gene symbols + control codewords) |
| `experiment.xenium` | 2–3 KB | run metadata |

There is no `cell_feature_matrix.h5`, no `cell_feature_matrix/` directory, no
`transcripts.parquet`, no boundary parquets and no morphology images. The
`Xenium_BL6` folder is empty.

`cells.parquet` carries **per-cell transcript totals only — not per-gene
counts**. Galr1 cannot be quantified from it, and neither can cell types. Every
downstream step in the analysis plan depends on the cell × gene matrix.

What is needed, per section (~50–150 MB each, ~1.5 GB for all 14):

```
cell_feature_matrix.h5            # required
# or the equivalent directory:
cell_feature_matrix/barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
```

Useful but not blocking: `cell_boundaries.parquet`, `nucleus_boundaries.parquet`,
`analysis/` (10x clustering), `morphology_focus/` (needed for Allen CCF
registration), `transcripts.parquet` (needed only for re-segmentation).

---

## 2. The panel contains Galr1 — the core question is answerable

Panel `mBrain_50g` = 247 predesigned (`mBrain_v1.1`) + 50 custom = **297 genes**,
plus 244 control codewords (41 negative-control codeword, 27 negative-control
probe, 175 unassigned, 1 deprecated).

Galanin system coverage:

| target | in panel |
|---|---|
| **Galr1** | **yes** |
| **Galr3** | **yes** |
| `Gal` (the ligand) | yes |
| Galr2 | **no** |

So both a receptor pair and the ligand can be measured in the same cells — that
supports a ligand–receptor story, not just a receptor count.

**Relevant gap:** the panel has no `Nr5a1`/`Sf1`, `Tbx3`, `Trh`, `Pmch`, `Kiss1`,
`Foxb1` or `Rax`. Those are the canonical transcription-factor markers used to
call VMH, ARC and mammillary territory. Nucleus assignment therefore cannot rely
on marker genes alone and needs spatial registration (plan §5).

Usable anchors that *are* present — ARC: `Agrp`, `Npy`, `Pomc`, `Cartpt`,
`Ghrh`, `Th`, `Slc6a3`, `Ghsr`, `Lepr`; VMH-leaning: `Adcyap1`, `Nts`, `Cckar`,
`Esr1`, `Fezf2`; DMH-leaning: `Grp`, `Pdyn`, `Gal`, `Cartpt`, `Npy`, `Prlh`;
tanycyte/ME: `Gpr50`, `Rfx4`, `Slc13a4`, `Col1a1`; LHA: `Hcrt`; PVN: `Oxt`,
`Avp`, `Crh`, `Otp`. Glia and vasculature are well covered (`Aqp4`, `Gfap`,
`Slc1a2`, `Ntsr2`, `Sox10`, `Opalin`, `Pdgfra`, `Cspg4`, `Trem2`, `Cd68`,
`Siglech`, `Spi1`, `Cldn5`, `Pecam1`, `Kdr`, `Acta2`, `Dcn`).

---

## 3. Two acquisition batches, and they are not interchangeable

| | cohort A (the ROI cohort) | cohort B |
|---|---|---|
| animals | F536, M493, G_073, M399 | K238, P953, F739, Q378 |
| sections | 10 (2–3 per animal) | 4 (1 per animal) |
| dates | 26 Jun / 01 Jul 2025 | 03 Nov / 06 Nov 2025 |
| panel design id | `7ZBFXR` | `NCY734` |
| analysis software | `xenium-3.3.0.1` | `xenium-4.0.1.0` |
| segmentation | 5 µm nucleus expansion (100%) | multimodal (nuc-expansion fraction 0.0) |
| transcripts / 100 µm² | 149–190 | 330–430 |
| fraction transcripts assigned | 0.79–0.82 | 0.69–0.72 |

Cohort B has roughly **twice** the transcript density and a different
segmentation algorithm. Pooling A and B into one differential test would put a
~2× technical gain on the same axis as the biological contrast. They can be used
as **independent cohorts that replicate each other**, which is scientifically
stronger than pooling anyway — but not as one merged dataset.

### Run structure within cohort A — the design is properly blocked

| batch | run | slide | animal | age | group | sections |
|---|---|---|---|---|---|---|
| B1 | 20250626 Run1 | 0014602 | F536 | 70 w (16.1 mo) | **aged** | F536_1, F536_2, F536_3 |
| B1 | 20250626 Run1 | 0014610 | M493 | 29 w (6.7 mo) | **adult** | M493_1, M493_2, M493_3 |
| B2 | 20250701 Run2 | 0014290 | G_073 | 65 w (15.0 mo) | **aged** | G_073_1, G_073_2 |
| B2 | 20250701 Run2 | 0014292 | M399 | 30 w (6.9 mo) | **adult** | M399_1, M399_3 |

All four are male. **Each batch contains one aged and one adult animal**, so age
is orthogonal to run day, slide and cassette. This is a randomised complete
block design: batch enters the model as a blocking factor, and every technical
drift shared by a run day cancels within the block. It also means the study
yields **two independent within-batch aged−adult differences**, and agreement
between them is the strongest evidence the design can produce at this n.

This is the best-case outcome for a 2 vs 2 layout, and it should be stated
explicitly in the manuscript — reviewers look for exactly this confound first.

### One framing caveat on "aged"

By the NIA/JAX convention for C57BL/6 (adult 3–6 mo, middle-aged 10–14 mo, old
18–24 mo), the two "aged" animals at **15.0 and 16.1 months are middle-aged**,
and the "adult" animals at 6.7–6.9 months sit at the upper edge of mature adult.
The contrast is real and biologically interesting — many hypothalamic ageing
phenotypes begin in middle age — but calling it "aged vs adult" without stating
the ages invites a reviewer objection at a high-impact journal. Recommended
framing: *mature adult (~7 months) vs middle-aged (~15–16 months)*, with the
ages given in the abstract or first figure, and any claim about ageing scoped to
"early/mid-life onset". If a genuinely old group (18–24 mo) can be added later,
that converts a caveat into a trajectory.

---

## 4. The ROI convention is solved and verified exactly

The ROI export stores an axis-aligned rectangle plus, for two sections, a manual
rotation. Rotating centroids by **+`rotation_deg` about `pivot`** and then
testing against the rectangle reproduces the recorded `n_cells_selected` with
**zero error on all four sections**:

| section | rotation | expected | no transform | **+deg** | −deg |
|---|---|---|---|---|---|
| M493_2 | +1.0° | 10 291 | 10 287 | **10 291** ✓ | 10 299 |
| F536_1 | — | 11 614 | **11 614** ✓ | — | — |
| G073_2 | — | 8 592 | **8 592** ✓ | — | — |
| M399_3 | −2.5° | 10 366 | 10 410 | **10 366** ✓ | 10 610 |

Implemented in `src/spatial_lea/roi.py`, locked in by `tests/test_roi.py`.

Note the ROI file uses key `G073_2` while Drive and `experiment.xenium` use
`G_073_2`; the mapping is handled in code.

### ROI geometry

| section | width | height | area | cells | density |
|---|---|---|---|---|---|
| M493_2 | 2386 µm | 1818 µm | 4.34 mm² | 10 291 | 2372 /mm² |
| F536_1 | 2679 µm | 1923 µm | 5.15 mm² | 11 614 | 2254 /mm² |
| G073_2 | 2373 µm | 1550 µm | 3.68 mm² | 8 592 | 2336 /mm² |
| M399_3 | 2350 µm | 1756 µm | 4.13 mm² | 10 366 | 2512 /mm² |
| **total** | | | 17.3 mm² | **40 863** | |

Cell density is consistent across sections (2254–2512 /mm², ~10% spread), which
argues against gross tissue-quality differences between animals.

The ROIs are generous rectangles over ventral tissue, so they will contain ARC,
VMH, DMH and ME **plus** neighbouring territory (LHA, PH, TU, and possibly
ventral PVN at the dorsal edge). Sub-nucleus assignment inside the ROI is
therefore mandatory, not optional — it is what the "is it the DMH?" question
turns on.

---

## 5. Sample-level QC from `experiment.xenium`

| section | cells | tx/cell | frac assigned | area mm² |
|---|---|---|---|---|
| F536_1 | 123 654 | 172 | 0.787 | 59.0 |
| F536_2 | 125 523 | 171 | 0.793 | 55.5 |
| F536_3 | 127 434 | 188 | 0.790 | 57.7 |
| M493_1 | 129 490 | 184 | 0.821 | 53.8 |
| M493_2 | 109 539 | 210 | 0.807 | 50.0 |
| M493_3 | 121 441 | 207 | 0.798 | 53.3 |
| G_073_1 | 116 145 | 217 | 0.785 | 56.0 |
| G_073_2 | 114 629 | 218 | 0.789 | 50.7 |
| M399_1 | 103 258 | 173 | 0.803 | 46.4 |
| M399_3 | 104 259 | 190 | 0.808 | 45.0 |
| K238_2 | 89 279 | 160 | 0.724 | 38.9 |
| P953_1 | 90 799 | 191 | 0.712 | 41.2 |
| F739_2 | 123 104 | 240 | 0.723 | 55.0 |
| Q378_2 | 127 025 | 262 | 0.692 | 59.0 |

Cohort A is tight and healthy (172–218 tx/cell, 0.79–0.82 assigned). Within
cohort A, transcripts/cell tracks run day somewhat (Run1 median ~188, Run2
~204), one more reason to keep run in the model.

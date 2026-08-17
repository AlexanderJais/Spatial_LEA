# Spatial_LEA — Galanin receptor signalling in the ageing mediobasal hypothalamus

Xenium spatial transcriptomics analysis of the mediobasal hypothalamus (MBH) in
mature adult (~7 months) vs middle-aged (~15–16 months) male mice, focused on
**Galr1** — which cell types express it, and whether it changes with age in
specific hypothalamic nuclei such as the DMH.

## Start here

| document | what it is |
|---|---|
| [`docs/ANALYSIS_PLAN.md`](docs/ANALYSIS_PLAN.md) | the analysis plan, statistics, deliverables and risks |
| [`docs/DATA_AUDIT.md`](docs/DATA_AUDIT.md) | what is actually in the data, verified directly |
| [`config/samples.yaml`](config/samples.yaml) | single source of truth for sample → animal → group |

## Current status

**Blocked on data.** The shared Drive folder contains `cells.parquet`,
`features.tsv.gz` and `experiment.xenium` for each section, but **not** the
cell × gene expression matrix (`cell_feature_matrix.h5`). Galr1 cannot be
quantified without it. See `docs/DATA_AUDIT.md` §1.

Already established and verified:

- The panel contains **`Galr1`, `Galr3` and `Gal`** (`Galr2` is absent) — the
  core question is answerable once the matrices arrive.
- The design is **properly blocked**: each batch holds one aged and one adult
  animal, so age is orthogonal to run day, slide and cassette.
- The **ROI transform convention is solved exactly** — rotate centroids by
  `+rotation_deg` about `pivot`, then apply the polygon. This reproduces the
  recorded `n_cells_selected` for all four sections with zero error.

## Layout

```
config/     sample manifest, ROI export, panel definition
docs/       analysis plan and data audit
src/        spatial_lea package (ROI selection implemented)
scripts/    01_fetch_data.py — mirrors the Drive folder into data/raw/
tests/      ROI counts locked against the recorded values
```

## Usage

```bash
pip install -r requirements.txt
pip install -e .

python scripts/01_fetch_data.py --list-only   # audit what is on Drive
python scripts/01_fetch_data.py               # mirror into data/raw/
pytest                                        # verifies ROI selection
```

`data/` is git-ignored; the Drive folder is the source of record.

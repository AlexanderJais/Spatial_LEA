# What statistics this design supports

Cohort A as it stands: **2 aged and 2 adult male mice, blocked** (each batch
holds one of each), **3 sections per animal**, 164 943 cells in the
hypothalamic window.

The unit of replication is the **animal**. Cells and sections are
pseudoreplicates for any age contrast, and aggregating them before testing is
what keeps the false-discovery rate honest. Everything below follows from that.

---

## The one hard limit

With 2 animals per group in a blocked design, the exact permutation null has
four label assignments (swap within B1 × swap within B2). The observed
assignment and its mirror are indistinguishable, so:

> **The smallest attainable exact two-sided p-value for any age contrast is
> 0.5 — for any gene, at any effect size.**

No test, model or transformation changes this. So the paper should **not report
a p-value for any age effect.** It should report effect sizes, calibrated
rankings, and the reproducibility criteria below. Stating this limit explicitly
is a strength in review, not a weakness.

---

## Tier 1 — fully supported, no age contrast involved

These are descriptive claims with enormous n and no pseudoreplication problem,
because the comparison is within-animal or across cells rather than across
groups. **Report each one replicated independently in all four animals.**

| claim | test | where |
|---|---|---|
| Galr1 is detected above background | detection rate vs the 27 negative-control probes; binomial CI, n = cells | Fig 3c |
| The Galr1-defining population is DMH | proportion with Wilson CI, per animal | Fig 3b |
| Galanin genes differ between cell types | within-animal comparison across cell types; χ² or Mann–Whitney on cells, repeated per animal | Fig 3a |
| Nuclei are correctly delineated | % of marker-defined cells landing in the expected nucleus | Fig 1 |
| Cell atlas composition | proportions with CIs, per animal | Fig 2 |

**Reporting form:** "in each of four animals independently, X" is stronger here
than any pooled p-value, and immune to the n = 2 problem.

## Tier 2 — legitimate inference, about nuisance rather than age

These have a valid unit of replication that is *not* the animal, so they carry
real p-values.

**Rostro-caudal confound test.** Effect size versus AP mismatch across all 18
possible section pairings. The unit is the section pair; the question is whether
the apparent effect tracks anatomical mismatch. Pearson correlation.

- `Galr1` in DMH Gal/Galr1 neurons: **r = −0.69, p = 0.002** → the effect grows
  as matching worsens. Artifact.
- `Gal` in ARC Agrp/Npy neurons: **r = +0.12, p = 0.63** → indifferent to
  matching. Survives.

This is the single most defensible inferential result in the study, and it is
also a methodological contribution: it quantifies a confound that spatial
hypothalamic studies generally do not control.

**Panel-null calibration.** Applying identical filters (consistent direction in
both blocks, ≥1.5× in the weaker block, complete section-level separation within
both blocks) to all 297 panel genes gives an empirical false-positive rate with
the same cells, animals and design.

- ~9% of the panel clears the first two filters in a typical population.
- Only **9 of 297 genes (3.0%)** clear all three in ARC Agrp/Npy, and `Gal` is
  one of them.

This is a legitimate multiplicity statement — an empirical FDR — and it does not
require an age p-value.

**Variance components.** Three sections per animal make the nesting estimable.

| component | SD (log2) | fixed by |
|---|---|---|
| between sections within an animal | **0.265** | AP matching |
| between animals within a group | **0.152** (2 df) | more animals |

**Anatomical sampling noise exceeds biological variation between animals.** That
single number justifies both the AP-matching protocol and the next cohort's
design. The animal term has only 2 df, so quote it as an estimate with that
caveat.

## Tier 3 — effect sizes, reported without p-values

Legitimate to report, provided no significance is claimed:

- **Blocked log2 fold change**, one per block, plus their mean. Always show both
  blocks separately — agreement between two independent batches is the evidence.
- **Per-animal values with per-animal points shown** (Fig 5a). Error bars should
  be labelled as within-animal measurement precision, not biological variance.
- **Cross-block sign consistency**, with the panel-wide pass rate as its
  denominator (Tier 2).
- **Complete section-level separation** — e.g. `Gal` in ARC Agrp/Npy is positive
  in **18 of 18** pairings. Descriptive reproducibility, not a p-value, because
  sections within an animal are not independent.

## Tier 4 — not available, do not attempt

- Any per-gene p-value, q-value or FDR for the **age** effect.
- Mixed-effects models with animal as a random effect: 4 animals cannot support
  a variance estimate for inference.
- Cell-level or section-level tests of the age effect (DESeq2/edgeR/Wilcoxon on
  cells or sections). These will produce tiny p-values and they will be wrong —
  this is textbook pseudoreplication.
- Claims that any effect is "significant", "not significant", or absent. The
  correct negative statement is bounded by the power analysis below.

---

## Power, and what the negative results mean

Minimum detectable |log2 FC| at 80% power, α = 0.05, blocked design, using the
pooled between-animal SD of 0.152:

| animals/group | as sampled | AP-matched |
|---|---|---|
| **2 (current)** | **3.47 (11.1×)** | 2.47 (5.5×) |
| 3 | 0.99 (2.0×) | 0.70 (1.6×) |
| 5 | 0.51 (1.4×) | 0.36 (1.3×) |
| 8 | 0.35 (1.3×) | 0.25 (1.2×) |
| 10 | 0.30 (1.2×) | 0.21 (1.2×) |

So the honest negative statement is:

> "This design had 80% power to detect only effects larger than ~11-fold. The
> absence of a detectable Galr1 age difference therefore excludes very large
> effects, and nothing smaller."

And the honest positive statement:

> "`Gal` in ARC Agrp/Npy neurons differs by 1.59× (log2 +0.67), consistent in
> both blocks, positive in all 18 section pairings, independent of rostro-caudal
> matching (r = +0.12, p = 0.63), and one of 9 of 297 panel genes passing all
> reproducibility filters. Confirming it at 80% power requires n ≥ 5 per group."

---

## What the paper is

Given the above, the defensible framing is a **spatially registered MBH atlas
with galanin-system mapping, plus a methodological result about rostro-caudal
confounding, plus one pre-registered candidate.** That is a complete, honest
paper. What it is not is a demonstration that galanin signalling changes with
age — and attempting that claim at n = 2 is what a reviewer will catch.

Three things would change the tier of the age claim, in order of value:

1. **n ≥ 5 per group with AP-matched sections** — moves the `Gal` candidate from
   Tier 3 to a real test.
2. **RNAscope on `Gal` in `Agrp`⁺ ARC neurons** in AP-matched sections — an
   orthogonal method, and the standard expectation for an in-situ claim.
3. **An 18–24-month group** — the current "aged" animals are middle-aged by the
   NIA convention (15.0 and 16.1 months), which should be stated in the abstract
   regardless.

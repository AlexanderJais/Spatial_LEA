#!/usr/bin/env bash
# Re-segment exported molecules with Baysor.
#
#   scripts/30_run_baysor.sh <section> [extra baysor args...]
#
# Baysor v0.7.1 (the Julia line, Petukhov et al. 2022) is run from source rather
# than the released binary: this container cannot reach GitHub releases, but it
# can clone the repository and resolve the Julia dependencies, so the package is
# installed from a v0.7.1 worktree.  scripts/setup_baysor.sh builds that worktree
# and applies patches/baysor-0.7.1-skip-ncv-colors.patch, which is required at
# this data size -- without it the run finishes segmenting and is then OOM-killed
# in a plotting-only step, before writing anything.  BAYSOR_HOME must point at
# the patched worktree.
#
# Parameters and why:
#
#   -m 30   minimum molecules for a cell to be called real.  Nuclei alone carry a
#           median of 119 molecules here, so 30 only removes fragments.
#   -s 6.0 --scale-std 50%
#           expected cell *radius* in um, and how wide a spread of sizes one run
#           may cover.  Chosen by scripts/32_baysor_scale_sweep.py, not by eye --
#           see results/baysor/scale_sweep.csv.  Scored against the measured DAPI
#           nuclei, 4 um has the fewest errors only because it never leaves the
#           nucleus (8.1 um cells against a 7.1 um nucleus, half the counts per
#           cell), while 8-10 um merges 21-29% of nuclei into shared cells.  6 um
#           is the only setting whose cell count its own cell size can account
#           for, and it is not estimated from the prior: the prior is nuclei only,
#           so estimating from it would bake the under-capture back in.
#   --prior-segmentation-confidence 0.5
#           the DAPI nuclei are real measurement, so they are trusted, but not
#           absolutely -- 0.5 lets Baysor split or merge where the molecular
#           composition clearly disagrees.
#   --n-clusters 8
#           coarse molecule clusters used as the composition prior.  The Baysor
#           docs recommend ~10 for large high-gene-count Xenium panels.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SECTION="${1:?usage: 30_run_baysor.sh <section> [extra args]}"
shift || true

: "${BAYSOR_HOME:?set BAYSOR_HOME to the Baysor v0.7.1 worktree}"
: "${JULIA_BIN:?set JULIA_BIN to the julia executable}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$BAYSOR_HOME/../.julia}"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-$(nproc)}"

IN="$REPO/data/segmentation/${SECTION}_molecules.csv"
OUT="$REPO/data/segmentation/${SECTION}_baysor"
[ -f "$IN" ] || { echo "missing $IN -- run scripts/29_export_transcripts.py first"; exit 1; }
mkdir -p "$OUT"

echo "=== Baysor: $SECTION ($(wc -l < "$IN") molecules, $JULIA_NUM_THREADS threads) ==="
time "$JULIA_BIN" --project="$BAYSOR_HOME" -e 'using Baysor; Baysor.command_main()' -- \
    run -x x -y y -z z -g gene -m 30 -s 6.0 --scale-std 50% \
    --prior-segmentation-confidence 0.5 \
    --n-clusters 8 \
    --count-matrix-format tsv \
    --polygon-format none \
    -o "$OUT" "$@" \
    "$IN" :prior

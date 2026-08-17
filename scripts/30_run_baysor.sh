#!/usr/bin/env bash
# Re-segment exported molecules with Baysor.
#
#   scripts/30_run_baysor.sh <section> [extra baysor args...]
#
# Baysor v0.7.1 (the Julia line, Petukhov et al. 2022) is run from source rather
# than the released binary: this container cannot reach GitHub releases, but it
# can clone the repository and resolve the Julia dependencies, so the package is
# installed from a v0.7.1 worktree.  BAYSOR_HOME must point at that worktree.
#
# Parameters and why:
#
#   -m 30   minimum molecules for a cell to be called real.  Nuclei alone carry a
#           median of 119 molecules here, so 30 only removes fragments.
#   -s 6.0  expected cell *radius* in um.  Not estimated from the prior, because
#           the prior is nuclei only (radius ~4 um) and estimating from it would
#           bake the under-capture we are trying to remove back in.  6 um is the
#           centre of a mixed neuron/glia population: real somata run 15-25 um
#           across for neurons and 8-12 um for glia.  Baysor treats it as a prior
#           mean with a 25% spread and adapts per cell.
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
    run -x x -y y -z z -g gene -m 30 -s 6.0 \
    --prior-segmentation-confidence 0.5 \
    --n-clusters 8 \
    --count-matrix-format tsv \
    --polygon-format none \
    -o "$OUT" "$@" \
    "$IN" :prior

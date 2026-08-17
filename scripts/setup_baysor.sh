#!/usr/bin/env bash
# Install Baysor v0.7.1 from source, into a directory of your choosing.
#
#   BAYSOR_ROOT=~/baysor scripts/setup_baysor.sh
#   export BAYSOR_HOME=~/baysor/Baysor JULIA_BIN=~/baysor/julia-1.10.5/bin/julia
#
# Why from source rather than the released binary: this container can reach
# github.com over git and can reach the Julia package registry, but not GitHub's
# releases endpoint, so the precompiled Baysor binary cannot be downloaded.
# Cloning the repository and resolving the dependencies through the registry
# works, and gives the same v0.7.1 code.
#
# Note that Baysor's master branch is now a C++ rewrite (0.8.x).  v0.7.1 is the
# last Julia release and the line the published method describes
# (Petukhov et al., Nat Biotechnol 2022), so that is the tag used here.
#
# One patch is applied, patches/baysor-0.7.1-skip-ncv-colors.patch.  It guards a
# plotting-only step behind the --plot flag; without it a 5M-molecule run
# allocates roughly 12 GB after the segmentation has already finished and is
# killed before writing its output.  The segmentation itself is unmodified.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${BAYSOR_ROOT:-$REPO/.baysor}"
JULIA_VERSION="${JULIA_VERSION:-1.10.5}"
TAG="${BAYSOR_TAG:-v0.7.1}"

mkdir -p "$ROOT"
cd "$ROOT"

if [ ! -x "julia-$JULIA_VERSION/bin/julia" ]; then
    echo "== fetching julia $JULIA_VERSION"
    curl -sSL -o julia.tar.gz \
        "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
    tar xzf julia.tar.gz && rm julia.tar.gz
fi
JULIA="$ROOT/julia-$JULIA_VERSION/bin/julia"

if [ ! -d Baysor ]; then
    echo "== cloning Baysor"
    git clone https://github.com/kharchenkolab/Baysor.git Baysor
fi
git -C Baysor fetch --tags origin
git -C Baysor checkout --quiet "$TAG"

if ! grep -q "PATCHED (Spatial_LEA)" Baysor/src/processing/utils/cli_wrappers.jl; then
    echo "== applying $(basename "$REPO/patches/baysor-0.7.1-skip-ncv-colors.patch")"
    git -C Baysor apply "$REPO/patches/baysor-0.7.1-skip-ncv-colors.patch"
fi

echo "== resolving dependencies (this takes ~10 min, CairoMakie dominates)"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$ROOT/.julia}"
"$JULIA" --project="$ROOT/Baysor" -e '
using Pkg
Pkg.instantiate()
Pkg.precompile()
using Baysor
println("Baysor ", pkgversion(Baysor), " ready")'

cat <<EOF

Done.  Export these before running scripts/30_run_baysor.sh:

  export BAYSOR_HOME=$ROOT/Baysor
  export JULIA_BIN=$JULIA
  export JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH
EOF

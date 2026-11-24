#!/usr/bin/env bash
# Compile lean/YMStatement/Reduction.lean against the pinned Mathlib olean cache.
# Usage: bash lean/check-reduction.sh
#
# Same v4.33.0-rc1 pin and LEAN_PATH scheme as lean/check-conditional.sh. The
# file imports Mathlib only (it is the self-contained reduction ledger), so no
# YMStatement/Handoff olean is needed first.
#
# Expected result: exit 0, with exactly FIVE `declaration uses 'sorry'` warnings
# (H1_holds, CP_holds, W_geo_nonneg_fact, P3_strong_coupling,
#  no_dissolution_kinematic). Every (A)/(C) theorem must show the standard axiom
# footprint only — no `sorryAx`.
set -u

ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="$MP$DEPS"

cd "$(dirname "$0")/.."

echo "=== lean/YMStatement/Reduction.lean"
~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean lean/YMStatement/Reduction.lean

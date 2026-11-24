#!/usr/bin/env bash
# Compile lean/YMStatement/H1Sharpening.lean against the pinned Mathlib olean cache.
# Usage: bash lean/check-h1sharpening.sh
#
# Same v4.33.0-rc1 pin and LEAN_PATH scheme as lean/check-reduction.sh. The
# file imports Mathlib only.
#
# Expected result: exit 0, ZERO `sorry` warnings, and every theorem's axiom
# footprint is only [propext, Classical.choice, Quot.sound] (no sorryAx).
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

echo "=== lean/YMStatement/H1Sharpening.lean"
~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean lean/YMStatement/H1Sharpening.lean

#!/usr/bin/env bash
# Reproducible check for the glue-d (GL7) deliverable.
# Usage: bash research/glue-d-dual-handoff/check-gl7.sh
#
# (1) Compiles GL7Dual.lean against the pinned Mathlib olean cache (same pin as
#     lean/check-handoff.sh, v4.33.0-rc1) and enforces the standard-axiom
#     footprint via #guard_msgs / #print axioms.
# (2) Runs the exact-rational collapse check, the frontier equivalence, and the
#     falsification cases for the two different mu*.
set -u

ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="$MP$DEPS"

cd "$(dirname "$0")/../.." || exit 1

echo "=== (1) Lean: research/glue-d-dual-handoff/GL7Dual.lean"
OUT=$(~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean research/glue-d-dual-handoff/GL7Dual.lean 2>&1)
echo "$OUT"
if echo "$OUT" | grep -q "sorryAx"; then
  echo "FAIL: sorryAx in axiom footprint"; exit 1
fi
if echo "$OUT" | grep -q "error"; then
  echo "FAIL: Lean error"; exit 1
fi

echo
echo "=== (2) numerical checks: research/glue-d-dual-handoff/check_gl7.py"
python -I research/glue-d-dual-handoff/check_gl7.py || exit 1

echo
echo "ALL GL7 CHECKS PASS"

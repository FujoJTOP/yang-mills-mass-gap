#!/usr/bin/env bash
# Reproducible check for the glue-a (GL3/GL4/GL9) deliverable.
# Usage: bash research/glue-a-limits/check-glue-a.sh
#
# (1) Compiles the four Lean files against the pinned Mathlib olean cache
#     (same pin as lean/check-handoff.sh, v4.33.0-rc1) and enforces the
#     standard-axiom footprint.  GL4L2.lean and GL3Limits.lean use
#     `#guard_msgs`; the other two print the axioms, which this script checks.
# (2) Runs the Gaussian sharpness computation for GL4-C (a numerical
#     counterexample to "weak limits preserve clustering without uniformity").
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

for f in \
  research/glue-a-limits/GL9Decay.lean \
  research/glue-a-limits/GL4Transfer.lean \
  research/glue-a-limits/GL4L2.lean \
  research/glue-a-limits/GL3Limits.lean
do
  echo "=== (1) Lean: $f"
  OUT=$(~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean "$f" 2>&1)
  echo "$OUT"
  if echo "$OUT" | grep -q "sorryAx"; then
    echo "FAIL: sorryAx in axiom footprint of $f"; exit 1
  fi
  if echo "$OUT" | grep -q "error"; then
    echo "FAIL: compile errors in $f"; exit 1
  fi
  if ! grep -q "#guard_msgs" "$f"; then
    if ! echo "$OUT" | grep -q "propext, Classical.choice, Quot.sound"; then
      echo "FAIL: unexpected axiom footprint in $f"; exit 1
    fi
  fi
  echo
done

echo "=== (2) GL4-C sharpness: massive 4D Gaussian, m -> 0"
python -I research/glue-a-limits/gl4_gaussian_counterexample.py || exit 1

echo
echo "ALL GLUE-A CHECKS PASS"

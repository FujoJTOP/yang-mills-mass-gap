#!/usr/bin/env bash
# Compile SYMDim4.lean against the pinned Mathlib olean cache and verify that
# every headline theorem has the standard axiom footprint.
# Usage: bash research/os2-rotation/check-symdim4.sh
#
# Pinned to v4.33.0-rc1 for the same reason as lean/check-handoff.sh: the
# vendored OSforGFF/Mathlib oleans are built with that toolchain. Uses the
# prebuilt cache only; it does not run lake and does not rebuild Mathlib.
set -u

ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="$MP$DEPS"

F="$(dirname "$0")/lean/SYMDim4.lean"
OUT=$(~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean "$F" 2>&1) || {
  echo "$OUT"; echo "COMPILE FAILED"; exit 1
}
echo "$OUT"

# The file ends with three `#guard_msgs in #print axioms` blocks that freeze the
# axiom footprint: if any theorem gained an axiom (including sorryAx), the guard
# fails and the compile above already failed. Reject any stray report as well.
BAD=$(printf '%s\n' "$OUT" | grep "depends on axioms" \
  | grep -v "\[propext, Classical.choice, Quot.sound\]" || true)
if [ -n "$BAD" ]; then
  echo "NON-STANDARD AXIOMS:"; echo "$BAD"; exit 1
fi
echo "SYMDIM4 OK: compiles; 3 #guard_msgs axiom guards passed (standard axioms only)"

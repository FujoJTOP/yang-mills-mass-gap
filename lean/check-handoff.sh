#!/usr/bin/env bash
# Compile all three handoff Lean files against the pinned Mathlib olean cache.
# Usage: bash lean/check-handoff.sh
#
# Pinned to v4.33.0-rc1 because the vendored OSforGFF/Mathlib oleans are built
# with that toolchain. The machine also has Lean 4.36; do not switch the pin
# without rebuilding Mathlib. `HandoffMargin.lean` carries #guard_msgs blocks,
# so this script also enforces the standard-axiom footprint.
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

for f in \
  lean/YMStatement/Handoff.lean \
  lean/YMStatement/RGHandoff.lean \
  research/handoff-proof/HandoffMargin.lean
do
  echo "=== $f"
  ~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean "$f" || exit 1
done

echo "ALL THREE COMPILE"

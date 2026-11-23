#!/usr/bin/env bash
# Compile the minimal-handoff lemma against the pinned Mathlib olean cache.
# Usage: bash research/minimal-handoff/check-minimal.sh
# Same pin and LEAN_PATH pattern as lean/check-conditional.sh (v4.33.0-rc1).
set -u

ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="D:/Dev/math/lean;D:/Dev/math/research/minimal-handoff;$MP$DEPS"

cd "$(dirname "$0")/../.."

LEAN=~/.elan/bin/elan
PIN=leanprover/lean4:v4.33.0-rc1

echo "=== research/minimal-handoff/MinimalHandoff.lean"
"$LEAN" run "$PIN" lean research/minimal-handoff/MinimalHandoff.lean || exit 1

echo "MINIMAL HANDOFF COMPILES"

#!/usr/bin/env bash
# Compile the conditional-theorem skeleton against the pinned Mathlib olean cache.
# Usage: bash lean/check-conditional.sh
#
# Same pin and LEAN_PATH as lean/check-handoff.sh (v4.33.0-rc1), plus the two
# directories whose modules Conditional.lean imports.  It builds the imported
# modules to .olean first so `import YMStatement.Handoff` / `import HandoffMargin`
# resolve; `Conditional.lean` carries #guard_msgs blocks, so this script also
# enforces the standard-axiom footprint.
set -u

ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="D:/Dev/math/lean;D:/Dev/math/research/handoff-proof;D:/Dev/math/research/glue-a-limits;$MP$DEPS"

cd "$(dirname "$0")/.."

LEAN=~/.elan/bin/elan
PIN=leanprover/lean4:v4.33.0-rc1

echo "=== lean/YMStatement/Handoff.lean (build olean)"
"$LEAN" run "$PIN" lean -o lean/YMStatement/Handoff.olean lean/YMStatement/Handoff.lean || exit 1

echo "=== research/handoff-proof/HandoffMargin.lean (build olean)"
"$LEAN" run "$PIN" lean -o research/handoff-proof/HandoffMargin.olean \
  research/handoff-proof/HandoffMargin.lean || exit 1

echo "=== research/glue-a-limits/GL9Decay.lean (build olean)"
"$LEAN" run "$PIN" lean -o research/glue-a-limits/GL9Decay.olean \
  research/glue-a-limits/GL9Decay.lean || exit 1

echo "=== lean/YMStatement/Conditional.lean"
"$LEAN" run "$PIN" lean lean/YMStatement/Conditional.lean || exit 1

echo "CONDITIONAL COMPILES (axioms frozen in #guard_msgs)"

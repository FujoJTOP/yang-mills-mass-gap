#!/usr/bin/env bash
# Compile the GL2a / GL2b' Lean artefacts against the pinned Mathlib olean cache.
# Usage: bash research/glue-b-blocking/check.sh
#
# Pinned to v4.33.0-rc1 to match the vendored Mathlib oleans (same pin as
# lean/check-handoff.sh; do not switch without rebuilding Mathlib).
# Enforces: (1) no sorry/admit/axiom in the sources (comments stripped),
#           (2) compilation, (3) the standard axiom footprint via #guard_msgs.
set -u

ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="$MP$DEPS"

cd "$(dirname "$0")"

fail=0

# 1. source scan with comments stripped (block `/- ... -/` and line `--` tails)
if ! python -I - "$@" <<'PY'
import re, sys
bad = re.compile(r'\b(sorry|admit|axiom)\b')
found = False
for path in ("lean/GL2a.lean", "lean/GL2bPrime.lean"):
    src = open(path, encoding='utf-8').read()
    # strip nested block comments (iterate to handle nesting)
    prev = None
    while prev != src:
        prev = src
        src = re.sub(r'/-.*?-/', '', src, flags=re.S)
    src = re.sub(r'--[^\n]*', '', src)
    for i, line in enumerate(src.splitlines(), 1):
        if bad.search(line):
            print(f"!!! {path}:{i}: {line.strip()}")
            found = True
    if found:
        break
sys.exit(1 if found else 0)
PY
then
  echo "!!! source scan failed"; fail=1
fi

# 2 + 3. compile; #guard_msgs freezes [propext, Classical.choice, Quot.sound]
for f in lean/GL2a.lean lean/GL2bPrime.lean; do
  echo "=== $f"
  if ! ~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean "$f"; then
    fail=1
  fi
done

if [ "$fail" -ne 0 ]; then
  echo "FAILED"
  exit 1
fi
echo "GL2a + GL2bPrime COMPILE (no sorry, standard axioms frozen)"
exit 0

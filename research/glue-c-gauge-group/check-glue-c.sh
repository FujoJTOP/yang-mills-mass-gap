#!/usr/bin/env bash
# GL1/GL6/GL8 evidence path for glue-c.  Reproduces:
#   - Lean: the group-blind assembly lemma (GL8)
#   - Python: the GL1 normalization dictionary and the GL6 orbit identity
# No files outside this directory are written.
#
# Usage: bash check-glue-c.sh
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$ROOT/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="$MP$DEPS"

cd "$HERE"

echo "=== Lean: GroupUniform.lean (GL8)"
~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean GroupUniform.lean || exit 1

echo
echo "=== Lean: GL6Fubini.lean (GL6 orbit Fubini + S3 enumeration)"
~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean GL6Fubini.lean || exit 1

echo
echo "=== GL1: normalization dictionary"
python gl1_check.py --group S3 --lattice 2x1 --beta 2.0 || exit 1

echo
echo "=== GL6: gauge-fixing independence (S3, Gribov copies)"
python gl6_check.py --group S3 --lattice 2x1 || exit 1

echo
echo "=== GL8: group-uniform constants"
python gl8_check.py || exit 1

echo
echo "ALL GLUE-C CHECKS PASS"

#!/usr/bin/env bash
# Compile-check files against the vendored mathlib olean cache (toolchain v4.33.0-rc1).
# Usage: ./check.sh YangMillsInfra/Some/File.lean
# Builds YangMillsInfra/Basic.olean first (the package's only root dependency).
set -u
ROOTPKG="D:/Dev/math/research/lean-infra"
DEPSRC="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$DEPSRC/mathlib/.lake/build.badcache/lib/lean"
DEPS=""
for p in batteries aesop Qq proofwidgets LeanSearchClient importGraph plausible; do
  DP="$DEPSRC/$p/.lake/build/lib/lean"
  [ -d "$DP" ] && DEPS="$DEPS;$DP"
done
export LEAN_PATH="$ROOTPKG;$MP$DEPS"
LEAN=~/.elan/bin/elan
cd "$ROOTPKG"
build() { # $1 = module path without .lean
  local f="$1.lean" o="$1.olean"
  if [ ! -f "$o" ] || [ "$f" -nt "$o" ]; then
    "$LEAN" run leanprover/lean4:v4.33.0-rc1 lean -o "$o" "$f" || return 1
  fi
}
build YangMillsInfra/Basic || exit 1
for target in "$@"; do
  "$LEAN" run leanprover/lean4:v4.33.0-rc1 lean -o "${target%.lean}.olean" "$target"
done

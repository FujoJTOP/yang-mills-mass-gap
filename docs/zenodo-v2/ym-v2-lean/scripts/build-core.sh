#!/usr/bin/env bash
# Build the two modules advertised in the paper (Reduction, H1Sharpening)
# plus the conditional-theorem module (Conditional), against the pinned Mathlib.
set -euo pipefail
cd "$(dirname "$0")/.."

lake update
# The official Mathlib cache may or may not exist for this revision. If it does
# not, this line fails and `lake build` below compiles Mathlib from source
# (tens of minutes to hours).
lake exe cache get || echo "WARNING: Mathlib cache unavailable for this revision; building from source."

lake build YMStatement.Reduction YMStatement.H1Sharpening YMStatement.Conditional
echo "CORE OK: Reduction, H1Sharpening, Conditional built."

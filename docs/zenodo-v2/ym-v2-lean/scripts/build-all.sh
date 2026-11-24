#!/usr/bin/env bash
# Build every library in the supplement. This is not needed for the paper's
# claims; it is for readers who want to check the surrounding development.
set -euo pipefail
cd "$(dirname "$0")/.."

lake update
lake exe cache get || echo "WARNING: Mathlib cache unavailable for this revision; building from source."

lake build YMStatement YangMillsInfra HandoffMargin GL9Decay GL2a GL6Fubini GL7Dual MinimalHandoff NonTrivialCore SYMDim4
echo "ALL OK."

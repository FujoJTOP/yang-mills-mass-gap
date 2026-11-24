import Lake
open Lake DSL

package «ym-v2-lean» where

require mathlib from git "https://github.com/leanprover-community/mathlib4.git" @ "79d0395a1825a6264ad5d269e35e60537518955e"

/-- The campaign's statement/reduction modules.
`YMStatement.Basic` is intentionally not a root: it imports the external
`OSforGFF` package, which is not redistributed here. -/
lean_lib «YMStatement» where
  roots := #[`YMStatement.Handoff, `YMStatement.Reduction,
             `YMStatement.H1Sharpening, `YMStatement.RGHandoff,
             `YMStatement.Conditional]

lean_lib «YangMillsInfra» where
  srcDir := "research/lean-infra"

lean_lib «HandoffMargin» where
  srcDir := "research/handoff-proof"

lean_lib «GL9Decay» where
  srcDir := "research/glue-a-limits"

lean_lib «GL2a» where
  srcDir := "research/glue-b-blocking/lean"

lean_lib «GL6Fubini» where
  srcDir := "research/glue-c-gauge-group"

lean_lib «GL7Dual» where
  srcDir := "research/glue-d-dual-handoff"

lean_lib «MinimalHandoff» where
  srcDir := "research/minimal-handoff"

lean_lib «NonTrivialCore» where
  srcDir := "research/nontriviality"

lean_lib «SYMDim4» where
  srcDir := "research/os2-rotation/lean"

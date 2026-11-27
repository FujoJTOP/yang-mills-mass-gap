import Lake
open Lake DSL

/-- OS1 minimal subset: four standalone Lean files compiled against the
mathlib already vendored at `vendor/OSforGFF` (v4.33.0-rc1).
`.lake/packages` is a junction to the vendored package store, so no mathlib
is fetched twice. -/
package «os1»

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.33.0-rc1"

@[default_target]
lean_lib «OS1» where
  roots := #[`OS1Tight, `OS1JointShift, `OS1Statement, `OS1Schwartz]

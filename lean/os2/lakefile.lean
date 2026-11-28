import Lake
open Lake DSL

/-- OS2 subset: the four OS2 files compiled against the mathlib already
vendored at `vendor/OSforGFF` (v4.33.0-rc1). `.lake/packages` is a junction
to the shared package store, so mathlib is not fetched twice. -/
package «os2» where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.33.0-rc1"

@[default_target]
lean_lib «OS2» where
  roots := #[`OS2Clip, `OS2Topology, `OS2Theta, `OS2Wiring]

/-- GL4 依赖文件（OS2Wiring 的 `rp_transfer` / `rp_l2_closure` 出处），随包一起编译。 -/
lean_lib «GL4» where
  roots := #[`GL4Transfer, `GL4L2]

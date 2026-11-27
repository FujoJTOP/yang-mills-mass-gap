/-
OS2Wiring.lean — ①②③ 并成一条：`rp_full_cone`

`import` 说明（**本文件未编译**；编译需要把两个交付目录放进 `LEAN_PATH`）：
  export LEAN_PATH="D:/Dev/math/vendor/OSforGFF/.lake/packages/mathlib/.lake/build.badcache/lib/lean\
;<…batteries/aesop/Qq/…>\;research/glue-a-limits;research/v3/species/ooo2"
  ~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean research/v3/species/ooo2/OS2Wiring.lean
（`research/glue-a-limits/check-glue-a.sh` 的 LEAN_PATH 再加两个冒号项即可。）

这条定理把 NNN2-OS2.md §3「三定理并一条」写出来：
* ③ 把 `hθ : MeasurePreserving (⇑θ) μ μ` 从**假设**降为**推论**（`OS2Theta.measurePreserving_of_weak_limit`）；
* ② 保证 `Ω` 有拓扑、`θ` 连续、`μ` 侧强收敛/唯一性都在位；
* ① 把 `GL4L2.trunc` 换成 `clipₙ`，于是 `rp_transfer`（要 `Ω →ᵇ ℝ`）能吃每一个 `clipₙ ∘ A`。
剩下两条输入，写在签名里、不再藏：
* `hRP : ∀ n i, 0 ≤ ∫ clipₙ(A(θω))·clipₙ(Aω) ∂μ_a`（μ_a 侧的 clip 层 RP；见 `OS2Clip` 末尾说明）；
* `hA2 : MemLp A 2 μ`（无界层要极限侧 `L²`，即 REPORT.md:238 的 F4-3 边界）。

依赖核验（mathlib4 v4.33.0-rc1；本文件**未编译**）
--------------------------------------------------
`GlueA.rp_transfer`                research/glue-a-limits/GL4Transfer.lean:48
`GlueA.rp_l2_closure`              research/glue-a-limits/GL4L2.lean:86
`GlueA.OS2.clipBCF_apply`          research/v3/species/ooo2/OS2Clip.lean
`GlueA.OS2.clip_memLp` / `clip_tendsto` / `clip_norm_le`   同上
`GlueA.OS2.measurePreserving_of_weak_limit`                research/v3/species/ooo2/OS2Theta.lean
-/
import Mathlib.Topology.Order.OrderClosed
import Mathlib.MeasureTheory.Measure.Portmanteau
import GL4Transfer
import GL4L2
import OS2Clip
import OS2Theta

open MeasureTheory Filter Topology
open scoped BoundedContinuousFunction Topology

namespace GlueA.OS2

variable {Ω : Type*} [TopologicalSpace Ω] [MeasurableSpace Ω] [BorelSpace Ω]
  [PseudoMetrizableSpace Ω]

/-- **GL4-RP 全锥，单条机检目标。**
`A` 只假设连续（不必有界）+ `L²(μ)`；测度保持**不再假设**，由 ③ 从 `μ_a` 侧推上来。 -/
theorem rp_full_cone
    (θe : Ω ≃ᵐ Ω) (hθc : Continuous (⇑θe))
    {ι : Type*} {F : Filter ι} [NeBot F]
    (μs : ι → ProbabilityMeasure Ω) (μ : ProbabilityMeasure Ω)
    (hweak : Tendsto μs F (𝓝 μ))
    (hinv : ∀ i, Measure.map (⇑θe) (μs i : Measure Ω) = (μs i : Measure Ω))
    (A : Ω → ℝ) (hAc : Continuous A) (hA2 : MemLp A 2 (μ : Measure Ω))
    (hRP : ∀ n i, 0 ≤ ∫ ω, clip n (A (⇑θe ω)) * clip n (A ω) ∂(μs i : Measure Ω)) :
    0 ≤ ∫ ω, A (⇑θe ω) * A ω ∂(μ : Measure Ω) := by
  -- ③：把 `rp_l2_closure` 的 `hθ` 由假设变成定理
  have hθinv : MeasurePreserving (⇑θe) (μ : Measure Ω) (μ : Measure Ω) :=
    measurePreserving_of_weak_limit θe hθc μs μ hweak hinv
  refine rp_l2_closure θe hθinv A hA2 (fun n ω => clip n (A ω)) ?_ ?_ ?_ ?_
  · -- `hAR`：`clipₙ ∘ A ∈ L²`（由 `|clipₙ| ≤ |·|`）
    intro n
    exact clip_memLp hAc hA2 n
  · -- `hpt`：逐点收敛
    intro ω
    exact clip_tendsto (A ω)
  · -- `hbound`：`‖clipₙ (A ω)‖ ≤ ‖A ω‖`
    intro n ω
    exact clip_norm_le n (A ω)
  · -- `hRP` 在极限测度上：①+`rp_transfer` 逐 `n` 供给
    intro n
    have h := rp_transfer (⇑θe) hθc μs μ hweak (clipBCF n A hAc) (fun i => ?_)
    · simpa only [clipBCF_apply] using h
    · simpa only [clipBCF_apply] using hRP n i

end GlueA.OS2

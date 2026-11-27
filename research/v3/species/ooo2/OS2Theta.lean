/-
OS2Theta.lean — ③ `θ_* μ = μ`：声明 + 它到底需要什么

仓内现状
--------
`grep "Measure.map|map_eq|map_apply|MeasurePreserving|eq_map" @ research/glue-a-limits/*.lean`
→ 3 命中，全是 `MeasurePreserving` **假设**（`GL4L2.lean:74,87,211`）；定理 **0** 条。
即：`rp_l2_closure` 把「θ 测度保持」当输入（`hθ : MeasurePreserving (⇑θ) μ μ`），
而 `REPORT.md:162` 第 1 步（`θ_*μ = μ`）只是三行散文。

本文件
------
* `map_eq_of_measurePreserving`：`MeasurePreserving ⇒ Measure.map = μ`。它是**字段选择**（`h.map_eq`），
  所以 GL4L2 的 `hθ` **字面上就是结论**——真正要证的是「极限测度也保持」，见下一条。
* `map_eq_of_weak_limit`：**真内容**——逐 a 不变 + 弱收敛 ⇒ 极限不变。
* `measurePreserving_of_weak_limit`：同一事实的 `MeasurePreserving` 形式，
  直接填 `rp_l2_closure` 的 `hθ` 槽位。

需要什么 / 不需要什么（这是本条判词的核心）
------------------------------------------
需要：(i) `Continuous θ`（⇒ 推前映射 `ν ↦ ν.map θ` 弱连续；
          `ProbabilityMeasure.tendsto_map_of_tendsto_of_continuous`）；
      (ii) `Tendsto μs F (𝓝 μ)`（弱收敛，**不是**抽象的 `hconv`）；
      (iii) 逐 a 的 `θ_*μ_a = μ_a`（物理输入；F4-4 是它的失效口：反射面/分块不对齐）；
      (iv) `ProbabilityMeasure Ω` 的 Hausdorff 性（极限唯一）← `[TopologicalSpace.PseudoMetrizableSpace Ω] [BorelSpace Ω]`。
不需要：`Ω` 紧、tightness、Prokhorov、任何 `a`-一致矩界、`A` 的 `L²` 界。
         也不要求 `θ` 是**对合**或双射——`measure.map θ μ = μ` 本身不需要可逆性；
         但 `rp_l2_closure` 的 `Ω ≃ᵐ Ω` 槽位要，那由 `OS2Topology.coordSwapMeasEquiv` 供给。

依赖核验（mathlib4 v4.33.0-rc1，逐个核到源；本文件**未编译**）
--------------------------------------------------------------
`MeasurePreserving.map_eq`（字段）  Dynamics/Ergodic/MeasurePreserving.lean:47
`ProbabilityMeasure.tendsto_map_of_tendsto_of_continuous`
                                    MeasureTheory/Measure/ProbabilityMeasure.lean:639
`ProbabilityMeasure.toMeasure_map`（`@[simp]`）   同上 :614
`ProbabilityMeasure.toMeasure_injective`          同上 :128
`ProbabilityMeasure.t2Space`                      同上 :422（要 `[HasOuterApproxClosed Ω] [BorelSpace Ω]`）
`HasOuterApproxClosed` 全局实例（由 `PseudoMetrizableSpace`）  MeasureTheory/Measure/HasOuterApproxClosed.lean:236 附近
`tendsto_nhds_unique`               Topology/Separation/Hausdorff.lean:179
`Continuous.measurable`（`[OpensMeasurableSpace]` ⇒）  MeasureTheory/Constructions/BorelSpace/Basic.lean:183 给实例
-/
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure
import Mathlib.MeasureTheory.Measure.HasOuterApproxClosed
import Mathlib.MeasureTheory.Constructions.BorelSpace.Basic
import Mathlib.Topology.Separation.Hausdorff
import Mathlib.Topology.Metrizable.Basic

open MeasureTheory Filter Topology
open scoped ENNReal

namespace GlueA.OS2

variable {Ω : Type*} [TopologicalSpace Ω] [MeasurableSpace Ω] [BorelSpace Ω]
  [TopologicalSpace.PseudoMetrizableSpace Ω]

/-! ## 1. 平凡版：`MeasurePreserving` 的定义就是 `map = μ` -/

/-- `MeasurePreserving θ μ μ` 的字段 `map_eq` 逐字就是 `Measure.map θ μ = μ`。
这条**不构成内容**：它说明 `GL4L2.lean:87` 的 `hθ` 是「把结论当假设」。 -/
theorem map_eq_of_measurePreserving (θ : Ω → Ω) (μ : Measure Ω)
    (h : MeasurePreserving θ μ μ) : Measure.map θ μ = μ :=
  h.map_eq

/-! ## 2. 真内容：逐 `a` 不变 + 弱收敛 ⇒ 极限不变 -/

/-- **③ 的主定理**：`μs ⇀ μ`、每个 `μs i` 都 `θ`-不变 ⇒ `μ` 也 `θ`-不变。 -/
theorem map_eq_of_weak_limit
    (θ : Ω → Ω) (hθ : Continuous θ)
    {ι : Type*} {F : Filter ι} [NeBot F]
    (μs : ι → ProbabilityMeasure Ω) (μ : ProbabilityMeasure Ω)
    (hweak : Tendsto μs F (𝓝 μ))
    (hinv : ∀ i, Measure.map θ (μs i : Measure Ω) = (μs i : Measure Ω)) :
    Measure.map θ (μ : Measure Ω) = (μ : Measure Ω) := by
  -- (1) 推前在弱拓扑下连续：`θ` 连续是唯一结构条件
  have hcont : Tendsto (fun i => ProbabilityMeasure.map (μs i) hθ.measurable.aemeasurable)
      F (𝓝 (ProbabilityMeasure.map μ hθ.measurable.aemeasurable)) :=
    ProbabilityMeasure.tendsto_map_of_tendsto_of_continuous μs μ hweak hθ
  -- (2) 由逐 a 不变性，左边那族就是 `μs` 本身
  have hconst : (fun i => ProbabilityMeasure.map (μs i) hθ.measurable.aemeasurable) = μs := by
    funext i
    refine ProbabilityMeasure.toMeasure_injective ?_
    rw [ProbabilityMeasure.toMeasure_map]
    exact hinv i
  rw [hconst] at hcont
  -- (3) Hausdorff ⇒ 极限唯一
  have hmap : ProbabilityMeasure.map μ hθ.measurable.aemeasurable = μ :=
    tendsto_nhds_unique hcont hweak
  have hcoe : (ProbabilityMeasure.map μ hθ.measurable.aemeasurable : Measure Ω)
      = (μ : Measure Ω) := congrArg (fun ν : ProbabilityMeasure Ω => (ν : Measure Ω)) hmap
  rwa [ProbabilityMeasure.toMeasure_map] at hcoe

/-- 同一事实的 `MeasurePreserving` 形式：直接填 `GL4L2.lean:87` 的 `hθ` 槽位。 -/
theorem measurePreserving_of_weak_limit
    (θe : Ω ≃ᵐ Ω) (hθc : Continuous (⇑θe))
    {ι : Type*} {F : Filter ι} [NeBot F]
    (μs : ι → ProbabilityMeasure Ω) (μ : ProbabilityMeasure Ω)
    (hweak : Tendsto μs F (𝓝 μ))
    (hinv : ∀ i, Measure.map (⇑θe) (μs i : Measure Ω) = (μs i : Measure Ω)) :
    MeasurePreserving (⇑θe) (μ : Measure Ω) (μ : Measure Ω) :=
  ⟨θe.measurable, map_eq_of_weak_limit (⇑θe) hθc μs μ hweak hinv⟩

end GlueA.OS2

/-
判词要点（写在这里，供 OOO2-LEAN-OS2.md §3 引用）：
1. `θ_*μ = μ` **不需要** `Ω` 紧、不需要 tightness、不需要任何 `a`-一致矩界；
   它严格弱于 GL4-RP 的其它两段，只消费「θ 连续 + 逐 a 不变 + 弱收敛 + 极限唯一」。
2. 它对 `θ` 的要求只有**连续 + 可测**；不要求对合、不要求双射。
   `rp_l2_closure` 另外要的 `Ω ≃ᵐ Ω`（供 `nrm_comp` 用）由坐标置换的自同胚性单独供给。
3. 真正的物理输入是 `hinv : ∀ i, θ_*μ_a = μ_a`（F4-4）：反射面与分块不对齐时它先坏，
   `θ` 仍连续但不变性失效；这与 os3 的 `(B1)` 判据（`A ∈ W ∩ RWR⁻¹`）同源。
4. 因此 GL4L2 的 `hθ : MeasurePreserving (⇑θ) μ μ` 应当**降级**为定理，
   而不是继续当假设；降级后全锥链条里只剩 `hRP`（clip 层）与 `A ∈ L²` 两条输入。
-/

/-
OS2Topology.lean — ② 给 `Ω` 加拓扑（`GL4L2.lean:25` 的 `variable` 段没有 `TopologicalSpace`）

结论（三层）
------------
(L0) **抽象层**：`Ω` 带 `[TopologicalSpace] [TopologicalSpace.PseudoMetrizableSpace] [MeasurableSpace] [BorelSpace]`
     即足以让 GL4 全套需要的实例自动出现（下面用 `example … := inferInstance` 逐个钉死）：
       `OpensMeasurableSpace Ω`              ← `BorelSpace.opensMeasurable`
       `HasOuterApproxClosed Ω`              ← `PseudoMetrizableSpace` 的全局实例
       `T2Space (ProbabilityMeasure Ω)`      ← `ProbabilityMeasure.t2Space`
       `CompactSpace (ProbabilityMeasure Ω)` ← `Prokhorov`（另需 `[CompactSpace Ω] [T2Space Ω]`）
(L1) **具体层**（REPORT.md §3.1「公共空间」）：`Ω := κ → [−1,1]`（归一化 Wilson 环的坐标），
     拓扑 = 积拓扑，紧 = Tychonoff，Borel = `Pi.borelSpace`（要 `[Countable κ]`）；
     `θ` = 反射诱导的**坐标置换** `x ↦ x ∘ σ`：连续（`continuous_pi`），且是 `≃ₜ`、`≃ᵐ`。
(L2) **对 GL4L2 的补丁**：`variable` 段换成带拓扑的版本（见文件末注释；不改既有文件）。

依赖核验（mathlib4 v4.33.0-rc1，逐个核到源；本文件**未编译**）
--------------------------------------------------------------
`BorelSpace.opensMeasurable`       MeasureTheory/Constructions/BorelSpace/Basic.lean:183
`Pi.borelSpace`（要 `[Countable ι]` + 纤维第二可数） 同上 :645
`Subtype.borelSpace`               同上 :187
`Subtype.secondCountableTopology`  Topology/Bases.lean:875
`HasOuterApproxClosed` 全局实例    MeasureTheory/Measure/HasOuterApproxClosed.lean:236 附近
`ProbabilityMeasure.t2Space`       MeasureTheory/Measure/ProbabilityMeasure.lean:422
`instCompactSpaceProbabilityMeasure`  MeasureTheory/Measure/Prokhorov.lean:168
（该节变量：`[MeasurableSpace E] [TopologicalSpace E] [T2Space E] [BorelSpace E]`，Prokhorov.lean:65）
`compactSpace_Icc`                 Topology/Order/Compact.lean:112
`Function.compactSpace`            Topology/Compactness/Compact.lean:1170
`pseudoMetrizableSpace_pi`         Topology/Metrizable/Basic.lean:100
`PseudoMetrizableSpace.subtype`    同上 :96
`continuous_pi` / `continuous_apply`  Topology/Constructions.lean:757,761
`Homeomorph` 字段                   Topology/Homeomorph/Defs.lean:43–49
`Homeomorph.toMeasurableEquiv`     MeasureTheory/Constructions/BorelSpace/Basic.lean:563
-/
import Mathlib.Topology.Order.Compact
import Mathlib.Topology.Metrizable.Basic
import Mathlib.Topology.Metrizable.Real
import Mathlib.MeasureTheory.Constructions.BorelSpace.Basic
import Mathlib.MeasureTheory.Measure.HasOuterApproxClosed
import Mathlib.MeasureTheory.Measure.Prokhorov

open MeasureTheory Topology

universe u

namespace GlueA.OS2

/-! ## 1. 抽象层：四个类换 GL4 全套实例 -/

/-- `BorelSpace` ⇒ `OpensMeasurableSpace`（`GL4Transfer.lean:39` 的变量段要它）。 -/
example {Ω : Type*} [TopologicalSpace Ω] [MeasurableSpace Ω] [BorelSpace Ω] :
    OpensMeasurableSpace Ω := inferInstance

/-- 拟可度量化 ⇒ 闭集示性函数可由有界连续函数外逼近（`OS2Theta` 的 T2 用它）。 -/
example {Ω : Type*} [TopologicalSpace Ω] [TopologicalSpace.PseudoMetrizableSpace Ω] :
    HasOuterApproxClosed Ω := inferInstance

/-- 弱拓扑在 `ProbabilityMeasure Ω` 上是 Hausdorff ⇒ **极限唯一**（`θ_*μ = μ` 的唯一性步）。 -/
example {Ω : Type*} [TopologicalSpace Ω] [TopologicalSpace.PseudoMetrizableSpace Ω] [MeasurableSpace Ω]
    [BorelSpace Ω] : T2Space (ProbabilityMeasure Ω) := inferInstance

/-- 紧空间上的概率测度空间紧（Prokhorov）⇒ 弱极限**存在**（OS0/tightness 那一侧）。 -/
example {Ω : Type*} [TopologicalSpace Ω] [T2Space Ω] [MeasurableSpace Ω] [BorelSpace Ω]
    [CompactSpace Ω] : CompactSpace (ProbabilityMeasure Ω) := inferInstance

/-! ## 2. 具体层：`Ω := κ → [−1,1]`（Wilson 环坐标），`θ` = 坐标置换 -/

/-- 单个归一化 Wilson 环坐标，取值在 `[−1,1]`。 -/
abbrev LoopCoord : Type := Set.Icc (-1 : ℝ) 1

/-- 公共紧空间：坐标族的积。 -/
abbrev LoopSpace (κ : Type u) : Type u := κ → LoopCoord

example (κ : Type*) : TopologicalSpace (LoopSpace κ) := inferInstance
example (κ : Type*) : CompactSpace (LoopSpace κ) := inferInstance            -- Tychonoff
-- 拟可度量化需要可数坐标族（不可数积非度量化，故此处加 `[Countable κ]`）
example (κ : Type*) [Countable κ] : TopologicalSpace.PseudoMetrizableSpace (LoopSpace κ) :=
  inferInstance
example (κ : Type*) : T2Space (LoopSpace κ) := inferInstance
example (κ : Type*) [Countable κ] : SecondCountableTopology (LoopSpace κ) := inferInstance
example (κ : Type*) [Countable κ] : BorelSpace (LoopSpace κ) := inferInstance
example (κ : Type*) [Countable κ] : OpensMeasurableSpace (LoopSpace κ) := inferInstance
/-- OS0/Prokhorov 侧：概率测度空间紧。 -/
example (κ : Type*) [Countable κ] :
    CompactSpace (ProbabilityMeasure (LoopSpace κ)) := inferInstance

variable {κ : Type*}

/-- 反射诱导的坐标置换：`θ x i = x (σ i)`（`σ : κ ≃ κ` 由「反射把一条 Wilson 环映到另一条」给出）。 -/
def coordSwap (σ : κ ≃ κ) (x : LoopSpace κ) : LoopSpace κ := fun i => x (σ i)

/-- `θ` 连续 —— `GL4Transfer.lean:49` 的 `hθ : Continuous θ`，这就是它的出处。 -/
lemma continuous_coordSwap (σ : κ ≃ κ) : Continuous (coordSwap σ) :=
  continuous_pi fun i => continuous_apply (σ i)

/-- `θ` 是自同胚（不需要 `σ` 是对合；反射是对合只是物理侧的额外事实）。 -/
noncomputable def coordSwapHomeo (σ : κ ≃ κ) : LoopSpace κ ≃ₜ LoopSpace κ where
  toFun := coordSwap σ
  invFun := coordSwap σ.symm
  left_inv x := by funext i; simp [coordSwap]
  right_inv x := by funext i; simp [coordSwap]
  continuous_toFun := continuous_coordSwap σ
  continuous_invFun := continuous_coordSwap σ.symm

/-- 进 `GL4L2.lean:87` 的 `Ω ≃ᵐ Ω` 槽位（`rp_l2_closure` / `nrm_comp` 要它）。 -/
noncomputable def coordSwapMeasEquiv [Countable κ] (σ : κ ≃ κ) :
    LoopSpace κ ≃ᵐ LoopSpace κ :=
  (coordSwapHomeo σ).toMeasurableEquiv

end GlueA.OS2

/-
(L2) 对 `GL4L2.lean:25` 的补丁（**不改既有文件**，只把该行应当长成的样子记在此处）：

  -- 原文: variable {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]
  variable {Ω : Type*} [TopologicalSpace Ω] [TopologicalSpace.PseudoMetrizableSpace Ω]
    [MeasurableSpace Ω] [BorelSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]

效果：`OpensMeasurableSpace Ω`（`BorelSpace` 给）、`HasOuterApproxClosed Ω`（`PseudoMetrizableSpace` 给）、
`clipBCF : Ω →ᵇ ℝ` 可用 ⇒ `trunc` 换成 `clipₙ` 后 `rp_transfer` 能直接吃。
注意 `BorelSpace` 而非仅 `OpensMeasurableSpace`：`ProbabilityMeasure.t2Space` 与
`ProbabilityMeasure.tendsto_map_of_tendsto_of_continuous` 都要求目标侧 `BorelSpace`。
-/

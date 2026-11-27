/-
  OOO3 · CompactTransfer.lean   —   接口层（只写类型，不写模型）
  ================================================================

  口径：固定 `(L, a)`、Wilson 作用、**紧**群 `G`、`0 ≤ β < ∞`、`a > 0`；
  `T_a = TransferData.transfer a`、`Ω = vacuum`、`λ₀ = 1`、`λ₁ = ‖T_a|_{Ω^⊥}‖`、
  `Δ_a = −(1/a)·log λ₁`。

  本文件只提供**类型**：

  1. `CompactTransfer`：`TransferData` 的四条加料（`compact / quadPos / kernel_pos /
     row_normalized`）。`TransferData` 本体（`Reconstruction/OS.lean:98–116`）**零改动**。
  2. `RepresentsKernel`：B 级桥（`L²(μ)` 表示）——**单列一条**，不塞进四字段（理由见 md §1）。
  3. `gapOp / lambdaOne / Delta / StepGap`：谱侧**定义**。

  未编译（本轮禁区：无 `lake build`、无 `#check`、无 `#eval`）。
  每个 `sorry` 带 `SORRY n` 标签，台账见 `OOO3-LEAN-FIXED.md` §4。

  逐行依赖（模块 → 用到的名字）：
    YangMillsInfra.Reconstruction.OS      `TransferData` 及其 7 个字段
    Mathlib.Analysis.InnerProductSpace.Basic      `inner`、`‖·‖`
    Mathlib.Analysis.InnerProductSpace.Orthogonal `Submodule.ᗮ`、`Submodule.mem_orthogonal`
    Mathlib.Analysis.InnerProductSpace.Adjoint    `ContinuousLinearMap.IsSymmetric`
    Mathlib.Topology.Algebra.Module.ContinuousLinearMap.Restrict  `ContinuousLinearMap.restrict`
    Mathlib.Analysis.Normed.Operator.Compact.Basic `IsCompactOperator`
    Mathlib.Analysis.SpecialFunctions.Log.Basic    `Real.log`
    Mathlib.MeasureTheory.Integral.Bochner.Basic   `∫ · ∂μ`、`Integrable`
    Mathlib.MeasureTheory.Measure.ProbabilityMeasure `IsProbabilityMeasure`
-/
import YangMillsInfra.Reconstruction.OS
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.InnerProductSpace.Orthogonal
import Mathlib.Analysis.InnerProductSpace.Adjoint
import Mathlib.Topology.Algebra.Module.ContinuousLinearMap.Restrict
import Mathlib.Analysis.Normed.Operator.Compact.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure

noncomputable section

open scoped InnerProductSpace
open MeasureTheory

namespace YangMillsInfra
namespace Reconstruction

/-! ## 1 · 真空线与压缩（纯 Hilbert 空间，不含 `TransferData`） -/

section Hilbert

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]

/-- 真空线 `ℂ·Ω`。 -/
def vacuumLine (Ω : E) : Submodule ℂ E := Submodule.span ℂ {Ω}

/-- `A` 保持真空线（只要 `A Ω = Ω` 与线性）。 -/
theorem mapsTo_vacuumLine (A : E →L[ℂ] E) (Ω : E) (hΩ : A Ω = Ω) :
    ∀ v ∈ vacuumLine Ω, A v ∈ vacuumLine Ω := by
  intro v hv
  induction hv using Submodule.span_induction with
  | mem x hx =>
      rw [Set.mem_singleton_iff] at hx
      subst hx
      rw [hΩ]
      exact Submodule.subset_span (by simp)
  | zero => rw [map_zero]; exact Submodule.zero_mem _
  | add x y _ _ hx hy => rw [map_add]; exact Submodule.add_mem _ hx hy
  | smul c x _ hx => rw [map_smul]; exact Submodule.smul_mem _ c hx

/-- **压缩**：`A` 限制到 `Ω^⊥`。自伴性给出"`Ω^⊥` 不变"，喂给
`ContinuousLinearMap.restrict`（`Topology/Algebra/Module/ContinuousLinearMap/Restrict.lean:31`）。 -/
def gapOp (A : E →L[ℂ] E) (hAsym : A.IsSymmetric) (Ω : E) (hΩ : A Ω = Ω) :
    (vacuumLine Ω)ᗮ →L[ℂ] (vacuumLine Ω)ᗮ :=
  A.restrict (by
    intro x hx
    -- 用无撇版 `mem_orthogonal'`（`⟪x, v⟫` 序），与 `hAsym x v` 的朝向一致
    rw [Submodule.mem_orthogonal'] at hx ⊢
    intro v hv
    rw [hAsym x v]
    have hv' : A v = v := by
      -- SORRY 0 · `A` 在 `ℂ·Ω` 上是恒等（`hΩ` + `Submodule.span_induction`，约 5 行）。
      sorry
    rw [hv']
    exact hx v hv)

end Hilbert

/-! ## 2 · `CompactTransfer`：四字段 -/

/-- **四条不可去的输入。** `FFF4-LAM1.md:24` 的 5×5 见证：去 `kernel_pos` / 去 `compact` /
去 `row_normalized` 各给 `λ₁ = 1`；去 `quadPos` 则 `λ₁` 可为负（`Δ_a` 无定义）。

`T : TransferData` 作**参数**而非 `extends`：`TransferData` 的
`instNormed / instInner / instComplete` 是普通字段、不是 instance binder，`extends` 之后
新字段类型里的 `T.transfer` 找不到 `NormedAddCommGroup T.H`；写成参数、每个字段类型自带
`letI` 即通过。四字段全 `Prop`，故整个结构取 `Prop`。 -/
structure CompactTransfer (T : TransferData) (X : Type*) [MeasurableSpace X]
    (μ : Measure X) [IsProbabilityMeasure μ] (κ : ℝ → X → X → ℝ) : Prop where
  /-- ① `T_a` 紧。 -/
  compact : letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ a : ℝ, 0 < a → IsCompactOperator (T.transfer a)
  /-- ② `T_a` 正（Loewner 二次型非负）。 -/
  quadPos : letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ a : ℝ, 0 ≤ a → ∀ f : T.H, 0 ≤ (inner ℂ f (T.transfer a f)).re
  /-- ③ 核严格正。 -/
  kernel_pos : ∀ a : ℝ, 0 < a → ∀ x y : X, 0 < κ a x y
  /-- ④ 每行归一（`T_a 1 = 1` 的函数形）。 -/
  row_normalized : ∀ a : ℝ, 0 < a → ∀ x : X, ∫ y, κ a x y ∂μ = 1

/-! ## 3 · B 级桥：核表示（单列，不塞进四字段） -/

/-- **`L²(μ)` 表示，只保留论证真正用到的那两条推论。**

（a）`Ω^⊥` 内的不动向量可换成一个**实值** `L²` 函数 `u`：均值 `0`、满足核方程；
（b）该表示**忠实**：`u = 0 a.e.` ⇒ `ψ = 0`。

由真正的 `L²(μ)` 模型推出这两条 = B 级（本仓 `T_a` 尚不存在：`Lattice.lean:38–42`
只有 `V / L / Plaq`、无时间层）。单列成假设，是为了让谱侧（A 级）保持可证。 -/
def RepresentsKernel (T : TransferData) (X : Type*) [MeasurableSpace X] (μ : Measure X)
    (κ : ℝ → X → X → ℝ) : Prop :=
  letI := T.instNormed; letI := T.instInner; letI := T.instComplete
  ∀ a : ℝ, 0 < a → ∀ ψ : T.H, T.transfer a ψ = ψ →
    (∀ v ∈ vacuumLine T.vacuum, inner ℂ ψ v = 0) →
    ∃ u : X → ℝ, Integrable u μ ∧ (∫ x, u x ∂μ = 0) ∧
      (∀ᵐ x ∂μ, u x = ∫ y, κ a x y * u y ∂μ) ∧
      ((∀ᵐ x ∂μ, u x = 0) → ψ = 0)

/-! ## 4 · 谱侧定义：`λ₁`、`Δ_a`、单步 gap -/

namespace TransferData

/-- `T_a` 自伴。（A 级胶水：`IsSelfAdjoint → IsSymmetric` 的 CLM 版；线性映射版见
`Adjoint.lean:710` 的 `isSymmetric_iff_isSelfAdjoint`。名字本机无法 `#check`，故 `sorry`。） -/
theorem isSymmetric_transfer (T : TransferData) {a : ℝ} (ha : 0 ≤ a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    (T.transfer a).IsSymmetric := by
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  -- SORRY 5 · `T.selfAdjoint a ha` 给出 `IsSelfAdjoint (T.transfer a)`；
  -- 需要 CLM 版 `IsSelfAdjoint → IsSymmetric`（`star = adjoint` + `adjoint_inner_left`）。
  sorry

/-- `T_a` 压缩到 `Ω^⊥`；`a ≥ 0` 只为拿到 `T_a Ω = Ω`（`vacuum_invariant`）。

注意类型里也要 `letI`：`ᗮ` 与 `→L[ℂ]` 都要 `T.H` 上的实例，而 `TransferData` 的实例
是字段不是 instance binder。 -/
def gapOp (T : TransferData) (a : ℝ) (ha : 0 ≤ a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    (vacuumLine T.vacuum)ᗮ →L[ℂ] (vacuumLine T.vacuum)ᗮ :=
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  Reconstruction.gapOp (T.transfer a) (T.isSymmetric_transfer ha) T.vacuum
    (T.vacuum_invariant a ha)

/-- **`λ₁`**：`Ω^⊥` 上压缩算子的算子范数。对正紧自伴算子它等于 `Ω^⊥` 上谱的上确界，
即"次一点"；这就是 `λ₁` 的无坐标定义（不用 Hilbert–Schmidt、不用本征向量基）。 -/
def lambdaOne (T : TransferData) (a : ℝ) (ha : 0 ≤ a) : ℝ :=
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  ‖T.gapOp a ha‖

/-- **`Δ_a = −(1/a)·log λ₁`**，`λ₁ = 0` 分支钉成 `1`。

钉死的原因：Lean 的 `Real.log 0 = 0`，于是 `−(1/a)·log 0 = 0`，会给出 `Δ_a = 0` 而
不是 `> 0`；而 `λ₁ = 0`（`β = 0` 的自由情形）时任意正 `Δ` 都对，取 `1`。 -/
def Delta (T : TransferData) (a : ℝ) (ha : 0 ≤ a) : ℝ :=
  if T.lambdaOne a ha = 0 then 1 else -(1 / a) * Real.log (T.lambdaOne a ha)

/-- **单步 gap**（时间 `a`、单一 `Δ`）。这是四字段真正给的东西。

与 `Spectral.SemigroupGap`（`Spectral/Transfer.lean:37`）的区别：后者要求**全体**
`t ≥ 0` 共用同一个 `Δ`。四字段对**每个** `a > 0` 各给一个 `Δ_a`，但 `t ↦ Δ_t` 的
一致性（半群插值 + cutoff 一致）是另一条输入，不在四字段内 —— 见 md §0 判词第 3 行。 -/
structure StepGap (T : TransferData) (a Δ : ℝ) : Prop where
  Δ_pos : 0 < Δ
  decay : letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ ψ : T.H, (∀ v ∈ vacuumLine T.vacuum, inner ℂ ψ v = 0) →
      ‖T.transfer a ψ‖ ≤ Real.exp (-Δ * a) * ‖ψ‖

end TransferData

end Reconstruction
end YangMillsInfra

/-
  OOO3 · LambdaOne.lean   —   主定理骨架：四字段 ⇒ 顶端本征值单重 ⇒ `λ₁ < 1` ⇒ `Δ_a > 0`
  =====================================================================================

  三层，界线画死：

  * **A 级（谱，本文件 §A/§C/§D）**：`gapOp_norm_lt_one` —— 紧 + 自伴 + `0 ≤ T ≤ 1` +
    "`1` 的本征空间只有 `Ω`" ⇒ `‖T|_{Ω^⊥}‖ < 1`。只用 `IsCompactOperator`，
    **不碰 Hilbert–Schmidt**（Mathlib 全库 0 命中，见 md §3）。
  * **B 级（分析，本文件 §B）**：`fixed_point_eq_const`（极值原理：正核 + 行归一 ⇒
    不动向量 a.e. 常值）＋ `RepresentsKernel` 的忠实性。**未证，`sorry`。**
  * **收尾（§D）**：`Δ_a = −(1/a)log λ₁ > 0`，以及单步 `StepGap`。

  未编译（本轮禁区）。`SORRY n` 台账见 `OOO3-LEAN-FIXED.md` §4。

  依赖（除 `CompactTransfer.lean` 外，逐行）：
    Mathlib.Analysis.InnerProductSpace.Spectrum
        `ContinuousLinearMap.orthogonalComplement_iSup_eigenspaces_eq_bot`(:443)、
        `ContinuousLinearMap.finite_dimensional_eigenspace`(:463)、
        `IsSymmetric.conj_eigenvalue_eq_self`、`IsSymmetric.restrict_invariant`
    Mathlib.Analysis.InnerProductSpace.Rayleigh
        `ContinuousLinearMap.spectralRadius_eq_nnnorm`(:182)、
        `IsSelfAdjoint.hasEigenvector_of_isMaxOn`(:276)
    Mathlib.Analysis.Normed.Operator.Compact.FredholmAlternative
        `IsCompactOperator.hasEigenvalue_iff_mem_spectrum`(:221)
    Mathlib.Analysis.Normed.Operator.Compact.Basic
        `IsCompactOperator.restrict'`(:349)、`IsCompactOperator.restrict`(:336)
    Mathlib.Analysis.Normed.Algebra.Spectrum
        `spectrum.isCompact`(:150)、`spectrum.exists_nnnorm_eq_spectralRadius_of_nonempty`(:229)
    Mathlib.Analysis.InnerProductSpace.Orthogonal
        `Submodule.isClosed_orthogonal`(:114)、`IsClosed.completeSpace_coe`
-/
import CompactTransfer
import Mathlib.Analysis.InnerProductSpace.Spectrum
import Mathlib.Analysis.InnerProductSpace.Rayleigh
import Mathlib.Analysis.Normed.Operator.Compact.FredholmAlternative
import Mathlib.Analysis.Normed.Algebra.Spectrum

noncomputable section

open scoped InnerProductSpace
open MeasureTheory

namespace YangMillsInfra
namespace Reconstruction

/-! ## §A · 谱侧核心（A 级，Mathlib 全覆盖，除两处胶水） -/

section Spectral

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E] [CompleteSpace E]

/-- **核心引理。** `A` 紧、自伴、`‖A‖ ≤ 1`、二次型非负、`A Ω = Ω`，而 `1` 的本征空间
在 `Ω^⊥` 上是零 ⇒ 压缩范数 `< 1`。

证明路线（每步都点名 Mathlib）：
 假设 `1 ≤ ‖S‖`（`S := gapOp A …`）⇒ `spectralRadius ℂ S = ‖S‖₊ ≥ 1`
 [`spectralRadius_eq_nnnorm`, Rayleigh:182] ⇒ 谱紧 [`spectrum.isCompact`] 且非空
 ⇒ 取到 `μ ∈ spectrum` 使 `‖μ‖ = spectralRadius` [`exists_nnnorm_eq_…
 _of_nonempty`, Spectrum:229] ⇒ `μ ≠ 0` ⇒ `μ` 是本征值
 [`hasEigenvalue_iff_mem_spectrum`, FredholmAlternative:221] ⇒ 取本征向量 `ψ ≠ 0`；
 压缩不增范 ⇒ `‖μ‖ ≤ 1`；`0 ≤ S` ⇒ `μ ≥ 0`；本征值实 ⇒ `μ = 1`；
 抬回 `E` 得 `A ψ = ψ`、`ψ ⊥ Ω`，与 `hsimple` 矛盾。 -/
theorem gapOp_norm_lt_one
    (A : E →L[ℂ] E) (hAsym : A.IsSymmetric) (hAk : IsCompactOperator A)
    (hA1 : ‖A‖ ≤ 1) (hApos : ∀ x : E, 0 ≤ (inner ℂ x (A x)).re)
    (Ω : E) (hΩ : A Ω = Ω) (hΩne : Ω ≠ 0)
    (hsimple : ∀ ψ : E, A ψ = ψ → ψ ∈ (vacuumLine Ω)ᗮ → ψ = 0) :
    ‖gapOp A hAsym Ω hΩ‖ < 1 := by
  set p := (vacuumLine Ω)ᗮ with hp
  by_cases hpbot : p = ⊥
  · -- 真空补退化：压缩算子落在零空间上，范数 `0`
    -- SORRY 9 · `‖gapOp A hAsym Ω hΩ‖ = 0`（`Subsingleton.elim` + `norm_zero`，约 3 行）
    sorry
  · haveI : Nontrivial p := Submodule.nontrivial_iff_ne_bot.mpr hpbot
    haveI : CompleteSpace p := p.isClosed_orthogonal.completeSpace_coe
    set S : p →L[ℂ] p := gapOp A hAsym Ω hΩ with hSdef
    -- 压缩仍紧：`IsCompactOperator.restrict'`（Compact/Basic.lean:349）以 `M →ₗ M` 陈述，
    -- 需在 CLM 与线性映射之间来回。
    have hSk : IsCompactOperator S := by
      -- SORRY 1 · 约 10 行
      sorry
    -- 压缩仍自伴：`IsSymmetric.restrict_invariant`（Symmetric.lean:137）的 CLM 包装。
    have hSsym : S.IsSymmetric := by
      -- SORRY 2 · 约 10 行（模式见 InnerProductSpace/Spectrum.lean:451）
      sorry
    have hSsa : IsSelfAdjoint S := hSsym.isSelfAdjoint
    by_contra hcon
    have hge : (1 : ℝ) ≤ ‖S‖ := not_lt.mp hcon
    -- 谱半径 = 范数
    have hr : spectralRadius ℂ S = ‖S‖₊ := S.spectralRadius_eq_nnnorm hSsa
    -- 谱紧 + 非空 ⇒ 取到
    obtain ⟨μ, hμmem, hμn⟩ :=
      spectrum.exists_nnnorm_eq_spectralRadius_of_nonempty (spectrum.nonempty S)
    -- `‖μ‖ = ‖S‖ ≥ 1`，故 `μ ≠ 0`
    have hμ0 : μ ≠ 0 := by
      rintro rfl
      have h0 : ‖S‖₊ = 0 := by rw [← hr, ← hμn]; simp
      have h0' : ‖S‖ = 0 := by simpa using h0
      linarith
    -- Fredholm：非零谱点 = 本征值
    have heig : HasEigenvalue S μ :=
      (IsCompactOperator.hasEigenvalue_iff_mem_spectrum hSk hμ0).mpr hμmem
    obtain ⟨ψ, hψ⟩ := heig.exists_hasEigenvector
    -- `HasEigenvector` 是 `x ∈ eigenspace μ ∧ x ≠ 0` 的 abbrev（Eigenspace/Basic.lean:417）
    have hψne : ψ ≠ 0 := hψ.2
    have hψe : S ψ = μ • ψ := hψ.apply_eq_smul
    -- 压缩不增范 ⇒ `‖μ‖ ≤ 1`
    have hμ1 : ‖μ‖ ≤ 1 := by
      -- SORRY 3 · `‖Sψ‖ = ‖A ψ‖ ≤ ‖A‖‖ψ‖ ≤ ‖ψ‖` 且 `‖Sψ‖ = ‖μ‖‖ψ‖`，约 10 行
      sorry
    -- `μ` 实、`μ ≥ 0`、`‖μ‖ = 1` ⇒ `μ = 1`
    have hμone : μ = 1 := by
      -- SORRY 4 · 三件事：`conj_eigenvalue_eq_self`（实）、二次型非负（`≥ 0`）、
      -- `‖μ‖ = ‖S‖` + `hge` + `hμ1`（`‖μ‖ = 1`）。约 25 行。
      sorry
    -- 抬回 `E`：`A ψ = ψ`
    have hfix : A (ψ : E) = (ψ : E) := by
      -- SORRY 6 · `↑(S ψ) = A ↑ψ`（`ContinuousLinearMap.restrict_apply`）+ `hμone`
      sorry
    exact hψne (Subtype.ext (hsimple (ψ : E) hfix ψ.property))

end Spectral

/-! ## §B · 极值原理 ＋ 核表示（B 级，两条 `sorry`） -/

/-- **极值原理**（`FFF4-LAM1.md:39` 的实值版）：正核 + 行归一 ⇒ 不动向量 a.e. 常值。

证明计划（约 200–400 行，本仓 `MeasureTheory` 零件未拼装）：`M := essSup u < ∞`
（`hint`）；`hfix` ⇒ `u ≤ M` a.e.；若 `{u < M − ε}` 有正测度，则用 `hpos` 的严格正性
与 `hrow` 得 `u ≤ M − δ` a.e.，与 `essSup = M` 矛盾；故 `u ≡ M` a.e.，`M = ∫u`。 -/
theorem fixed_point_eq_const {X : Type*} [MeasurableSpace X] {μ : Measure X}
    [IsProbabilityMeasure μ] {κ : X → X → ℝ}
    (hpos : ∀ x y, 0 < κ x y) (hrow : ∀ x, ∫ y, κ x y ∂μ = 1)
    {u : X → ℝ} (hint : Integrable u μ)
    (hfix : ∀ᵐ x ∂μ, u x = ∫ y, κ x y * u y ∂μ) :
    ∀ᵐ x ∂μ, u x = ∫ y, u y ∂μ := by
  -- SORRY 8 · 见上（B 级；`essSup` + 正性 + 行归一）
  sorry

/-- **顶端本征值单重**：`Ω^⊥` 内 `T_a` 的不动向量只有 `0`。四字段里 `kernel_pos` 与
`row_normalized` 在这一步被用掉；`quadPos` 用在 §A 的 `μ ≥ 0`。 -/
theorem TransferData.top_eigenspace_trivial (T : TransferData) {X : Type*} [MeasurableSpace X]
    {μ : Measure X} [IsProbabilityMeasure μ] {κ : ℝ → X → X → ℝ}
    (CT : CompactTransfer T X μ κ) (hrep : RepresentsKernel T X μ κ) {a : ℝ} (ha : 0 < a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ ψ : T.H, T.transfer a ψ = ψ →
      (∀ v ∈ vacuumLine T.vacuum, inner ℂ ψ v = 0) → ψ = 0 := by
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  intro ψ hfix horth
  obtain ⟨u, hint, hmean, hkernel, hfaith⟩ := hrep a ha ψ hfix horth
  have hconst : ∀ᵐ x ∂μ, u x = ∫ y, u y ∂μ :=
    fixed_point_eq_const (CT.kernel_pos a ha) (CT.row_normalized a ha) hint hkernel
  have hzero : ∀ᵐ x ∂μ, u x = 0 := by
    filter_upwards [hconst] with x hx
    rw [hx, hmean]
  exact hfaith hzero

/-! ## §C · `λ₁ < 1`（把四字段喂进 §A） -/

/-- **`λ₁ < 1`。** 这就是"固定 `(L,a)`"的全部内容：`Δ_a > 0` 由它一行推出。 -/
theorem TransferData.lambdaOne_lt_one (T : TransferData) {X : Type*} [MeasurableSpace X]
    {μ : Measure X} [IsProbabilityMeasure μ] {κ : ℝ → X → X → ℝ}
    (CT : CompactTransfer T X μ κ) (hrep : RepresentsKernel T X μ κ) {a : ℝ} (ha : 0 < a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    T.lambdaOne a ha.le < 1 := by
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  -- `top_eigenspace_trivial` 给的是"内积对每个 `v ∈ Ω` 为零"，这里换成成员形式喂给 §A
  have hsimple : ∀ ψ : T.H, T.transfer a ψ = ψ → ψ ∈ (vacuumLine T.vacuum)ᗮ → ψ = 0 := by
    intro ψ hfix hmem
    rw [Submodule.mem_orthogonal'] at hmem
    exact T.top_eigenspace_trivial CT hrep ha ψ hfix hmem
  exact gapOp_norm_lt_one (T.transfer a) (T.isSymmetric_transfer ha.le) (CT.compact a ha)
    (T.contraction a ha.le) (CT.quadPos a ha.le) T.vacuum (T.vacuum_invariant a ha.le)
    T.vacuum_ne hsimple

/-! ## §D · `Δ_a > 0` 与单步 decay -/

/-- **`Δ_a > 0`。** 分支 `λ₁ = 0` 走 `Delta` 的钉死值 `1`：Lean 里 `Real.log 0 = 0`，
不加这个分支会得到 `Δ_a = 0`。 -/
theorem TransferData.Delta_pos (T : TransferData) {X : Type*} [MeasurableSpace X]
    {μ : Measure X} [IsProbabilityMeasure μ] {κ : ℝ → X → X → ℝ}
    (CT : CompactTransfer T X μ κ) (hrep : RepresentsKernel T X μ κ) {a : ℝ} (ha : 0 < a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    0 < T.Delta a ha.le := by
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  have hq : T.lambdaOne a ha.le < 1 := T.lambdaOne_lt_one CT hrep ha
  have hq0 : 0 ≤ T.lambdaOne a ha.le := norm_nonneg _
  unfold TransferData.Delta
  split_ifs with h0
  · norm_num
  · have hpos : 0 < T.lambdaOne a ha.le := lt_of_le_of_ne hq0 (Ne.symm h0)
    have hlog : Real.log (T.lambdaOne a ha.le) < 0 := Real.log_neg hpos hq
    have hinv : -(1 / a : ℝ) < 0 := by
      have : (0 : ℝ) < 1 / a := by positivity
      linarith
    exact mul_pos_of_neg_of_neg hinv hlog

/-- **单步 decay**（`t = a`，单一 `Δ_a`）：四字段真正给的东西。

**不是** `SemigroupGap`：后者要全体 `t ≥ 0` 共用一个 `Δ`。从四字段只能得
`n·a` 这一条算术级数上的 `λ₁^n`（见 `stepGap_nat`），全体 `t` 需要半群插值 +
cutoff 一致的 `Δ₀` —— 那是靶，不是四字段。 -/
theorem TransferData.stepGap (T : TransferData) {X : Type*} [MeasurableSpace X]
    {μ : Measure X} [IsProbabilityMeasure μ] {κ : ℝ → X → X → ℝ}
    (CT : CompactTransfer T X μ κ) (hrep : RepresentsKernel T X μ κ) {a : ℝ} (ha : 0 < a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    T.StepGap a (T.Delta a ha.le) := by
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  refine ⟨T.Delta_pos CT hrep ha, ?_⟩
  intro ψ horth
  have hmem : ψ ∈ (vacuumLine T.vacuum)ᗮ := by
    rw [Submodule.mem_orthogonal']
    exact horth
  have hcomp : T.transfer a ψ = T.gapOp a ha.le ⟨ψ, hmem⟩ := by
    -- SORRY 10 · `gapOp` 与 `T_a` 在 `Ω^⊥` 上一致（`ContinuousLinearMap.restrict_apply`）
    sorry
  calc ‖T.transfer a ψ‖
      = ‖T.gapOp a ha.le ⟨ψ, hmem⟩‖ := by rw [hcomp]
    _ ≤ ‖T.gapOp a ha.le‖ * ‖ψ‖ := ContinuousLinearMap.le_opNorm _ _
    _ = T.lambdaOne a ha.le * ‖ψ‖ := rfl
    _ = Real.exp (-(T.Delta a ha.le) * a) * ‖ψ‖ := by
        -- SORRY 11 · `exp(−Δ_a·a) = λ₁`：`λ₁ > 0` 用 `Real.exp_log`（注意 `Δ_a` 的定义
        -- 让 `−Δ_a·a = log λ₁`）；`λ₁ = 0` 分支 `exp·‖ψ‖ ≥ 0` 也够。
        sorry

/-- **算术级数上的 decay**：`t = n·a`。半群律（`T.semigroup`）+ `Ω^⊥` 不变性
（自伴 + `vacuum_invariant`）⇒ `λ₁^n`。这是四字段能推到的"全体 `t`"的边界。 -/
theorem TransferData.stepGap_nat (T : TransferData) {X : Type*} [MeasurableSpace X]
    {μ : Measure X} [IsProbabilityMeasure μ] {κ : ℝ → X → X → ℝ}
    (CT : CompactTransfer T X μ κ) (hrep : RepresentsKernel T X μ κ) {a : ℝ} (ha : 0 < a) :
    letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ n : ℕ, ∀ ψ : T.H, (∀ v ∈ vacuumLine T.vacuum, inner ℂ ψ v = 0) →
      ‖T.transfer ((n : ℝ) * a) ψ‖ ≤ T.lambdaOne a ha.le ^ n * ‖ψ‖ := by
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  -- SORRY 12 · 归纳：`n → n+1` 用 `T.semigroup`；每步要 `Ω^⊥` 不变
  -- （自伴 + `vacuum_invariant`，约 15 行）。
  sorry

end Reconstruction
end YangMillsInfra

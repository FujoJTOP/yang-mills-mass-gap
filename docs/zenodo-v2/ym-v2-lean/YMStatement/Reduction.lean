/-
# 条件归约链的机器检验 · Machine-checking the conditional reduction chain

本文件把 4D 纯 SU(N) Yang–Mills 质量间隙项目里 **"H1 ∧ (CP) ∧ P3 ∧ … ⇒ m_gap > 0"**
这条**归约逻辑**形式化，并把其中**真正可证**的结构性引理（Feshbach–Bogoliubov Schur 引理、
susceptibility 挤压、`G″ = E[H] − Var + W` 恒等式）以 **0 sorry** 证完。

## 诚实边界（对应 CLAUDE.md 的铁律）

* 本文件 **不证明 Clay 问题**，也不证明质量间隙。
* `(A)` 段全部是 **PROVEN（0 sorry）**：纯线性代数 + 有序域代数。
* `(B)` 段是 **STATED（sorry）**：它们是**未证的解析输入**，每条都在注释里写明精确缺失内容。
* `(C)` 段的 `clay_gap_conditional` 是**机器检验的核心**：它把 `(B)` 的未证输入作为
  **假设**，只证明 **`假设 ⇒ m_gap>0`** 这一步逻辑。它 **不调用任何 (B) 的 sorry 定理**，
  因而恒为 **0 sorry**。`H1` 始终是假设，绝不被声称已证。

## 记号

* `f₁²` = 边原子残数（edge residue），`χ` = susceptibility，`λ` = 粗模 Poincaré 常数。
* `H1`：`Var_ν(φ) ≤ (1/λ)∫|∇φ|²`，等价于 `χ ≤ 1/λ`（本文件用后者）。
* `(CP)`：禁闭相 + Kato-保持 + 一致隔离，这里取其实质内容 = **乘性残数搬运**
  `f₁² = Z_O²·f₁²_strong`，`Z_O>0`（见 `docs/new-tools-registry.md`、INVENT-TRANSPORT）。
* `P3`：强耦合残数下界 `f₁²_strong ≥ r₀ > 0`（`β < 1/48`，`PROOF-P3-UNCONDITIONAL.md`）。
* `W_geo = E[Φ″] − Var(Φ′)`：Haar 几何项（Prékopa 凸性，`PROOF-JACOBIAN-SAVES.md`）。
* `G″ = E[H] − Var + W`：粗模曲率分解（G2.0 恒等式层）。
-/

import Mathlib

noncomputable section

namespace YMStatement.Reduction

open Matrix
open scoped Matrix

/-! ## (A) 完全证明 · FULLY PROVEN (0 sorry) -/

/-! ### A1. Feshbach–Bogoliubov 双块 Schur 引理

块算子 `[[A, B], [Bᵀ, D]]`，`A ≻ 0`。Schur 补 `Σ = D − Bᵀ A⁻¹ B`。
若 `Σ ⪰ 0`，则全算子 `⪰ 0`。这是纯线性代数（Mathlib 的 `Matrix.PosDef.fromBlocks₁₁`）。
`A ≻ 0` 蕴含 `A ⪰ aI`（对某个 `a>0`），即把 `G″ = S − D` 的"用 `A` 配平方"那一步严格化：
`Q(x,y) = (x + A⁻¹By)ᴴ A (x + A⁻¹By) + yᴴ Σ y ≥ 0`。 -/

/-- **Feshbach–Bogoliubov / Schur 补引理**（实矩阵，0 sorry）。

`A ≻ 0`、`Σ := D − Bᵀ A⁻¹ B ⪰ 0` ⇒ `[[A, B],[Bᵀ, D]] ⪰ 0`。 -/
theorem feshbach_bogoliubov_schur {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m]
    {A : Matrix m m ℝ} {B : Matrix m n ℝ} {D : Matrix n n ℝ} [Invertible A]
    (hA : A.PosDef) (hSigma : (D - Bᵀ * A⁻¹ * B).PosSemidef) :
    (Matrix.fromBlocks A B Bᵀ D).PosSemidef := by
  -- 把 `Bᵀ` 写成 Mathlib 引理中的 `Bᴴ`（实数上二者相等）。
  have hSigma' : (D - Bᴴ * A⁻¹ * B).PosSemidef := by
    rwa [← Matrix.conjTranspose_eq_transpose_of_trivial] at hSigma
  have h := (Matrix.PosDef.fromBlocks₁₁ (A := A) B D hA).mpr hSigma'
  rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-! ### A2. `G″ = E[H] − Var + W` 恒等式与正性胶合

文档 G2.0 的精确恒等式 `G″(k→0) = E[H] − Var(∇S) + W`，以及 `G″ = S − D`（`S=E[H]`，
`D=Var−W`）。这是**定义性/代数**恒等式；只要 `Var ≤ E[H]` 且 `W ≥ 0` 就有 `G″ ≥ 0`。
自由场（`E[H]=Var=W=0`）上 `G″=0`——过"自由场酸测"。 -/

/-- 粗模曲率分解 `G″ = E[H] − Var + W`。 -/
def gDouble (EH Var W : ℝ) : ℝ := EH - Var + W

/-- 定义性恒等式（`rfl`）。 -/
theorem gDouble_eq (EH Var W : ℝ) : gDouble EH Var W = EH - Var + W := rfl

/-- 等价形式 `G″ = S − D`，其中 `S = E[H]`、`D = Var − W`。 -/
theorem gDouble_eq_S_sub_D (EH Var W : ℝ) :
    gDouble EH Var W = EH - (Var - W) := by
  unfold gDouble; ring

/-- 正性胶合：`Var ≤ E[H] ∧ 0 ≤ W ⇒ 0 ≤ G″`。 -/
theorem gDouble_nonneg {EH Var W : ℝ} (hEV : Var ≤ EH) (hW : 0 ≤ W) :
    0 ≤ gDouble EH Var W := by
  unfold gDouble; linarith

/-- 自由场酸测：`E[H]=Var=W=0 ⇒ G″=0`（工具在自由场上失效，非 vacuous）。 -/
theorem gDouble_free : gDouble 0 0 0 = 0 := by unfold gDouble; norm_num

/-! ### A3. 方差 / susceptibility 挤压（有序域代数）

REDUCTION 的精确形式：`f₁²/χ ≤ inf supp ν ≤ M₀/χ`。
充分条件：`f₁² ≥ r₀ ∧ χ ≤ 1/λ ⇒ m₁² = f₁²/χ ≥ r₀·λ > 0`。
全部是线性有序域上的代数，0 sorry。 -/

/-- 下侧：`f ≤ χ·lower ⇒ f/χ ≤ lower`。 -/
theorem div_le_of_le_mul {f χ lower : ℝ} (hχ : 0 < χ) (h : f ≤ χ * lower) :
    f / χ ≤ lower := by
  rw [div_le_iff₀ hχ]; linarith

/-- 上侧：`χ·upper ≤ M₀ ⇒ upper ≤ M₀/χ`。 -/
theorem le_div_of_mul_le {upper M₀ χ : ℝ} (hχ : 0 < χ) (h : χ * upper ≤ M₀) :
    upper ≤ M₀ / χ := by
  rw [le_div_iff₀ hχ]; linarith

/-- **两步挤压** `f₁²/χ ≤ inf supp ν ≤ M₀/χ`（0 sorry）。 -/
theorem susceptibility_squeeze {f2 M0 infv χ : ℝ} (hχ : 0 < χ)
    (hlow : f2 ≤ χ * infv) (hup : χ * infv ≤ M0) :
    f2 / χ ≤ infv ∧ infv ≤ M0 / χ :=
  ⟨div_le_of_le_mul hχ hlow, le_div_of_mul_le hχ hup⟩

/-- **充分条件**：`0<f`、`χ ≤ 1/λ` ⇒ `f·λ ≤ f/χ`（即 `m₁² = f/χ ≥ f·λ`）。 -/
theorem m_sq_ge_f_lam {f χ lam : ℝ} (hf : 0 < f) (hχ : 0 < χ) (hlam : 0 < lam)
    (hχle : χ ≤ 1 / lam) : f * lam ≤ f / χ := by
  have hχmul : χ * lam ≤ 1 := by
    have := mul_le_mul_of_nonneg_right hχle (le_of_lt hlam)
    rwa [one_div, inv_mul_cancel₀ (ne_of_gt hlam)] at this
  have h2 : lam ≤ 1 / χ := by
    rw [le_div_iff₀ hχ]; nlinarith [hχmul]
  calc f * lam ≤ f * (1 / χ) := mul_le_mul_of_nonneg_left h2 (le_of_lt hf)
    _ = f / χ := by ring

/-- **精确挤压陈述**（与 REDUCTION 字面对应，0 sorry）：
`0<f`、`f ≤ χ·lower`、`χ ≤ 1/λ` ⇒ `f·λ ≤ m = f/χ ≤ lower`。 -/
theorem variance_susceptibility_squeeze {f lower χ lam : ℝ} (hf : 0 < f) (hχ : 0 < χ)
    (hlam : 0 < lam) (hle : f ≤ χ * lower) (hχle : χ ≤ 1 / lam) :
    f * lam ≤ f / χ ∧ f / χ ≤ lower :=
  ⟨m_sq_ge_f_lam hf hχ hlam hχle, div_le_of_le_mul hχ hle⟩

/-- **REDUCTION 的充分条件版**：`r₀>0`、`χ ≤ 1/λ`、`r₀ ≤ f₁²` ⇒ `r₀·λ ≤ f₁²/χ`。 -/
theorem m_sq_ge_of_squeeze {r0 lam χ f2 : ℝ} (hr0 : 0 < r0) (hlam : 0 < lam)
    (hχ : 0 < χ) (hχle : χ ≤ 1 / lam) (hf : r0 ≤ f2) :
    r0 * lam ≤ f2 / χ := by
  calc r0 * lam ≤ r0 / χ := m_sq_ge_f_lam hr0 hχ hlam hχle
    _ ≤ f2 / χ := div_le_div_of_nonneg_right hf (le_of_lt hχ)

/-! ## 理论数据包 · the bundled theory data

把归约链**消费**的实数量打包。`IsPhysical` 是诚实标记：这个 bundle 是物理的
4D 纯 SU(N) YM 连续理论。它**不是** `True` 占位——它是一个抽象的 `Prop` 字段，**只**
在 `(B)` 的未证事实里作为守卫出现，**不****在** `(C)` 的归约定理里出现（归约定理直接以
`H1`、`W_geo≥0` 等**假设**为输入，与"能否构造物理 bundle / 断言其物理性"解耦）。 -/

structure YMTheory where
  /- 强耦合残数（P3 一侧）。 -/
  f1sqStrong : ℝ
  r0 : ℝ
  betaStrong : ℝ
  /- 乘性残数搬运因子（CP 一侧）。 -/
  ZO : ℝ
  /- 物理残数。 -/
  f1sq : ℝ
  /- susceptibility 与粗模 Poincaré 常数（H1）。 -/
  χ : ℝ
  lam : ℝ
  /- Haar 几何项（Prékopa）。 -/
  Wgeo : ℝ
  /- 曲率分解 `G″ = E[H] − Var + W` 的两项。 -/
  EH : ℝ
  Var : ℝ
  /- 物理质量间隙（速率）。 -/
  mgap : ℝ
  /- 诚实标记：本 bundle 是物理理论（抽象 Prop，不在此断言为真）。 -/
  IsPhysical : Prop

/-- `P3`：强耦合残数下界 `0 < r₀ ∧ β<1/48 ∧ r₀ ≤ f₁²_strong`。 -/
def P3 (T : YMTheory) : Prop :=
  0 < T.r0 ∧ T.betaStrong < 1 / 48 ∧ T.r0 ≤ T.f1sqStrong

/-- `(CP)`：单相/Kato-保持的实质内容 = 乘性残数搬运 `f₁² = Z_O²·f₁²_strong`，`Z_O>0`。 -/
def CP (T : YMTheory) : Prop :=
  0 < T.ZO ∧ T.f1sq = T.ZO ^ 2 * T.f1sqStrong

/-- `H1`：粗模 Poincaré 的 susceptibility 形式 `0<λ ∧ χ≤1/λ`。 -/
def H1 (T : YMTheory) : Prop :=
  0 < T.lam ∧ 0 < T.χ ∧ T.χ ≤ 1 / T.lam

/-- `W_geo ≥ 0`（Prékopa 几何项非负）。 -/
def WGeoNonneg (T : YMTheory) : Prop := 0 ≤ T.Wgeo

/-- 搬运后的物理残数下界 `ρ = Z_O²·r₀`（P3 的种子经 CP 搬到物理标度）。 -/
def transportedRho (T : YMTheory) : ℝ := T.ZO ^ 2 * T.r0

/-- 谱质量平方 `m₁² = f₁²/χ`。 -/
def m1sq (T : YMTheory) : ℝ := T.f1sq / T.χ

/-- 无溶解的运动学裕度：`m₁²>0 ⇒ 物理质量间隙>0`（原子不溶入多粒子连续谱）。 -/
def NoDissolution (T : YMTheory) : Prop := 0 < m1sq T → 0 < T.mgap

/-- 质量间隙：`m_gap > 0`。 -/
def MassGap (T : YMTheory) : Prop := 0 < T.mgap

/-- `G″` 的 bundle 版本。 -/
def gDoubleOf (T : YMTheory) : ℝ := gDouble T.EH T.Var T.Wgeo

/-- **`W_geo≥0` 给出 `G″` 的下界一侧**：`E[H] − Var ≤ G″`（0 sorry）。

这是 `G″ = E[H] − Var + W` 与 `W_geo≥0` 的直接代数后果，对应项目 G2.0 里
"`G″` 形式与 H1 形式是同一对象"的下界半边。 -/
theorem gdouble_lower_bound (T : YMTheory) (hW : WGeoNonneg T) :
    T.EH - T.Var ≤ gDoubleOf T := by
  unfold gDoubleOf gDouble WGeoNonneg at *
  linarith

/-! ## (B) 以 sorry 陈述的未证解析输入 · STATED (sorry)

每条都在注释里写明 **精确缺失的分析内容**（CLAUDE.md：绝不把未证之物称定理）。
它们**不**被 `(C)` 调用；`(C)` 只把它们当作假设。 -/

/-- **H1（未证）** — 粗模 Poincaré / 方差上界 `Var_ν(φ) ≤ (1/λ)∫|∇φ|²`，
等价 `χ ≤ 1/λ`。

缺失内容：在 `(L,a)` 上一致的 `∇_∥`（`k→0`）Poincaré 常数下界——这就是质量间隙本身
（`docs/YANG-MILLS-MASTER.md` §6.1、`PROOF-SKELETON.md` G2.5）。 -/
theorem H1_holds (T : YMTheory) (_hT : T.IsPhysical) : H1 T := by
  sorry

/-- **(CP)（未证）** — 单 T=0 禁闭相 + 0⁺⁺ 边 Kato-保持 + 一致隔离，这里是其乘性搬运内核
`f₁² = Z_O²·f₁²_strong`，`Z_O>0`。缺失内容：跨 Savvidy 交叉的 Kato/隔离与连续极限
（`CONTINUUM-CP`）。 -/
theorem CP_holds (T : YMTheory) (_hT : T.IsPhysical) : CP T := by
  sorry

/-- **W_geo ≥ 0（未证）** — Haar 几何项 `W_geo = E[Φ″] − Var(Φ′) ≥ 0`。

缺失内容：Haar 纤维上有效作用的 **Prékopa 对数凹/凸性**（`PROOF-JACOBIAN-SAVES.md`
有 β=0 闭合与数值，但全 β 的无条件证明未落定）。 -/
theorem W_geo_nonneg_fact (T : YMTheory) (_hT : T.IsPhysical) : WGeoNonneg T := by
  sorry

/-- **P3 强耦合残数定理（未证）** — `β < 1/48 ⇒ 0 < r₀ ∧ r₀ ≤ f₁²_strong`。

缺失内容：强耦合 surface/character 展开的 graded-strip KP 引理与 `a₀>0` 的**一致余项**
（`PROOF-P3-UNCONDITIONAL.md`；目前只到 `β<1/48`，未搬运到 `β=2.5`）。 -/
theorem P3_strong_coupling (T : YMTheory) (_hT : T.IsPhysical) (hβ : T.betaStrong < 1 / 48) :
    0 < T.r0 ∧ T.r0 ≤ T.f1sqStrong := by
  sorry

/-- **无溶解运动学裕度（未证）** — 边原子的运动学裕度使其不自溶入多粒子连续谱：
`m₁²>0 ⇒ m_gap>0`。缺失内容：在连续极限下把边原子（非对角、usc）的零点与物理谱边
**一致地**等同（`PROOF-NO-DISSOLUTION-KINEMATIC.md` 给了固定 `a`/固定盒的运动学，连续
`a→0` 的一致版本未证）。 -/
theorem no_dissolution_kinematic (T : YMTheory) (_hT : T.IsPhysical) : NoDissolution T := by
  sorry

/-! ## (C) 归约定理 · THE REDUCTION THEOREM (0 sorry)

逻辑核心：`P3` 经 `(CP)` 搬运成物理残数下界 `ρ = Z_O²r₀ > 0`；`H1` 给 `χ ≤ 1/λ`;
A3 的挤压给 `m₁² = f₁²/χ ≥ ρ·λ > 0`;无溶解裕度把它转成 `m_gap>0`。
**所有未证输入都是假设；证明不触碰任何 sorry 定理。** -/

/-- **搬运后的残数下界**：`P3 ∧ CP ⇒ 0 < ρ ∧ ρ ≤ f₁²`。
（`ρ = Z_O² r₀`；`Z_O>0` 保号，`r₀ ≤ f₁²_strong` 经 `f₁² = Z_O² f₁²_strong` 搬运。） -/
theorem transported_residue (T : YMTheory) (hCP : CP T) (hP3 : P3 T) :
    0 < transportedRho T ∧ transportedRho T ≤ T.f1sq := by
  obtain ⟨hZO, heq⟩ := hCP
  obtain ⟨hr0, _hβ, hr0le⟩ := hP3
  have hZO2 : 0 < T.ZO ^ 2 := pow_pos hZO 2
  refine ⟨?_, ?_⟩
  · unfold transportedRho; positivity
  · unfold transportedRho
    rw [heq]
    exact mul_le_mul_of_nonneg_left hr0le (le_of_lt hZO2)

/-- **谱质量正性**：`H1 ∧ CP ∧ P3 ∧ 无溶解 ⇒ m_gap>0`。
只用到 (A) 的 `m_sq_ge_of_squeeze` 与 (C) 的 `transported_residue`，**0 sorry**。 -/
theorem mass_gap_of_spectral (T : YMTheory) (hH1 : H1 T) (hCP : CP T) (hP3 : P3 T)
    (hno : NoDissolution T) : MassGap T := by
  obtain ⟨hlam, hχ, hχle⟩ := hH1
  obtain ⟨hρpos, hρle⟩ := transported_residue T hCP hP3
  have hm : 0 < m1sq T := by
    unfold m1sq
    exact lt_of_lt_of_le (mul_pos hρpos hlam)
      (m_sq_ge_of_squeeze hρpos hlam hχ hχle hρle)
  exact hno hm

/-- **THE REDUCTION THEOREM（机器检验的核心，0 sorry）**

`H1 ∧ (CP) ∧ (W_geo≥0) ∧ P3 ∧ 无溶解裕度 ⇒ m_gap > 0`。

`W_geo≥0` 经 `gdouble_lower_bound` 给出 `G″ = E[H] − Var + W ≥ E[H] − Var`（`G″` 分解的
下界一侧，对应项目中 `G″` 形式与 H1 形式是同一对象 G2.0）；结论 `m_gap>0` 的逻辑由 A3
挤压承担。**`W_geo≥0` 对"`m_gap>0`"这一方向不是必需**——这里保留它作为归约链输入清单的
忠实记录（它服务于等价的 `G″` 形式），并显式驱动 `gdouble_lower_bound`。
**注意**：本定理只检验"假设 ⇒ 结论"的**逻辑有效性**，`H1` 等仍是**假设**，不是已证事实。 -/
theorem clay_gap_conditional (T : YMTheory)
    (hH1 : H1 T) (hCP : CP T) (hWgeo : WGeoNonneg T) (hP3 : P3 T)
    (hno : NoDissolution T) : MassGap T := by
  -- `W_geo≥0` ⇒ `G″ ≥ E[H] − Var`（G2.0 等价形式的下界一侧；此处显式记录）。
  have _hG := gdouble_lower_bound T hWgeo
  -- 主逻辑：H1 + CP + P3 + 无溶解 ⇒ m_gap>0。
  exact mass_gap_of_spectral T hH1 hCP hP3 hno

/-! ## 公理审计 · axiom audit

下面打印主定理的公理足迹。`(A)` 与 `(C)` 的定理**不得**含 `sorryAx`。
（`(B)` 的定理会含 `sorryAx`，因为它们就是 sorry 的输入。） -/

#print axioms feshbach_bogoliubov_schur
#print axioms gDouble_nonneg
#print axioms gdouble_lower_bound
#print axioms susceptibility_squeeze
#print axioms variance_susceptibility_squeeze
#print axioms m_sq_ge_of_squeeze
#print axioms transported_residue
#print axioms mass_gap_of_spectral
#print axioms clay_gap_conditional

end YMStatement.Reduction

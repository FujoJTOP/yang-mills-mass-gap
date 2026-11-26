/-
  L04PropC.lean — 命题 C（I2 ⟺ I1 ∧ ∃次小谱点）与引理 A1（s-加权平均 ≥ m₁）的 Lean 形式化

  工具状态：lean 4.36.0（elan）；mathlib **部分预编译**
  （D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build\lib\lean 下 2136/2357 模块）。
  可用：Mathlib.Data.Real.Basic、Mathlib.Algebra.Order.BigOperators.Group.Finset。
  **不可用**：Mathlib.Tactic（linarith/ring/norm_num 缺失）、MeasureTheory.*（无 olean）。
  ⇒ 不 import MeasureTheory；A1 的积分版只在 .md 里陈述。本文件证明不用数值自动化。

  编译（在 D:\Dev\math\vendor\OSforGFF 下）：
    lake env lean D:\Dev\math\research\v3\lean\L04PropC.lean
-/
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset

open scoped BigOperators

namespace L04

/-! ## 第 1 部分：引理 A1（有限和形式）—— s-加权平均 ≥ 最小点 -/

/-- 标量从和里提出来：`c * ∑ w = ∑ c * w`（mathlib 的 `Finset.mul_sum` 未预编译，自证）。-/
theorem mul_sum_eq {ι : Type*} (s : Finset ι) (w : ι → ℝ) (c : ℝ) :
    c * s.sum w = s.sum (fun i => c * w i) := by
  classical
  induction s using Finset.induction_on with
  | empty => rw [Finset.sum_empty, Finset.sum_empty, mul_zero]
  | insert a s ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, mul_add, ih]

/-- **引理 A1（家族形式）**：非负权、总权 1、每点 `m i ≥ mmin` ⇒ 加权平均 `≥ mmin`。
    即凸性引理（最小点 ≤ 加权平均）。 -/
theorem a1 {ι : Type*} [Fintype ι] (w m : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (h1 : ∑ i, w i = 1) (mmin : ℝ) (h : ∀ i, mmin ≤ m i) :
    mmin ≤ ∑ i, w i * m i := by
  have hgen : mmin * ∑ i, w i ≤ ∑ i, w i * m i := by
    rw [mul_sum_eq Finset.univ w mmin]
    refine Finset.sum_le_sum (fun i _ => ?_)
    have h2 : w i * mmin ≤ w i * m i := mul_le_mul_of_nonneg_left (h i) (hw i)
    rw [mul_comm mmin (w i)]
    exact h2
  rw [h1, mul_one] at hgen
  exact hgen

/-! ## 第 2 部分：命题 C 的精确陈述

支撑取 `List ℝ`（有限支撑、可重列）；权重为 `ℝ → ℝ`。
`Λ★` 不引入可计算函数，而用「存在且最小的次小点」谓词 `LtStar` 表达。

注意：`Spec` **不**内置「`pts` = 正权重支撑」这一条。是否内置决定 `I1` 是否自动成立
（见 `I1_automatic` 与反例 `counterexample_I2_not_I1`），这是命题 C 的唯一实质分歧点。
-/

/-- 有限支撑的谱：点表 `pts`（非空）、权重函数 `wt` 非负。 -/
structure Spec where
  pts : List ℝ
  hne : pts ≠ []
  wt : ℝ → ℝ
  hnn : ∀ x, 0 ≤ wt x

/-- `m1` 是谱的最小点。 -/
def IsMin (S : Spec) (m1 : ℝ) : Prop := m1 ∈ S.pts ∧ ∀ x ∈ S.pts, m1 ≤ x

/-- **I1**：`f₁² > 0`，即最小点上的权重（基态残数）为正。 -/
def I1 (S : Spec) (m1 : ℝ) : Prop := 0 < S.wt m1

/-- **次小谱点**：去掉最小点后支撑非空 ⇔ 存在严格大于 `m1` 的谱点。 -/
def SecExists (S : Spec) (m1 : ℝ) : Prop := ∃ x, x ∈ S.pts ∧ m1 < x

/-- `Λ★` = `supp ∖ {m1}` 的最小点（存在且最小，谓词形式）。 -/
def LtStar (S : Spec) (m1 : ℝ) (Λ : ℝ) : Prop :=
  m1 < Λ ∧ Λ ∈ S.pts ∧ ∀ x ∈ S.pts, m1 < x → Λ ≤ x

/-- **I2**：`Λ★ > m1²`。 -/
def I2 (S : Spec) (m1 : ℝ) : Prop := ∃ Λ, LtStar S m1 Λ ∧ m1 ^ 2 < Λ

/-! ### 辅助：有限点集最小元的存在性（纯归纳，不用 `Finset.min'`） -/

/-- 非空有限实数集有最小元。 -/
theorem exists_min : ∀ l : List ℝ, l ≠ [] → ∃ m, m ∈ l ∧ ∀ x ∈ l, m ≤ x
  | [], h => absurd rfl h
  | [a], _ => ⟨a, by simp, by
      intro x hx
      simp only [List.mem_singleton] at hx
      exact le_of_eq hx.symm⟩
  | a :: b :: t, _ => by
      have hne : b :: t ≠ [] := List.cons_ne_nil b t
      obtain ⟨m, hm, hmin⟩ := exists_min (b :: t) hne
      rcases le_total a m with hle | hle
      · refine ⟨a, by simp, ?_⟩
        intro x hx
        simp only [List.mem_cons] at hx
        rcases hx with rfl | hx
        · exact le_rfl
        · exact le_trans hle (hmin x (by simp only [List.mem_cons]; exact hx))
      · refine ⟨m, by simp [hm], ?_⟩
        intro x hx
        simp only [List.mem_cons] at hx
        rcases hx with rfl | hx
        · exact hle
        · exact hmin x (by simp only [List.mem_cons]; exact hx)

/-- 去掉最小点后的子集若有元素，则有最小元（次小谱点）。 -/
theorem exists_min_above (l : List ℝ) (m1 : ℝ) (h : ∃ x, x ∈ l ∧ m1 < x) :
    ∃ Λ, m1 < Λ ∧ Λ ∈ l ∧ ∀ x ∈ l, m1 < x → Λ ≤ x := by
  induction l with
  | nil => obtain ⟨x, hx, _⟩ := h; simp at hx
  | cons a t ih =>
      by_cases ha : m1 < a
      · by_cases ht : ∃ x, x ∈ t ∧ m1 < x
        · obtain ⟨Λ, hΛ, hmem, hmin⟩ := ih ht
          rcases le_total a Λ with hle | hle
          · exact ⟨a, ha, by simp, by
              intro x hx hxlt
              simp only [List.mem_cons] at hx
              rcases hx with rfl | hx
              · exact le_rfl
              · exact le_trans hle (hmin x hx hxlt)⟩
          · exact ⟨Λ, hΛ, by simp [hmem], by
              intro x hx hxlt
              simp only [List.mem_cons] at hx
              rcases hx with rfl | hx
              · exact hle
              · exact hmin x hx hxlt⟩
        · push_neg at ht
          exact ⟨a, ha, by simp, by
            intro x hx hxlt
            simp only [List.mem_cons] at hx
            rcases hx with rfl | hx
            · exact le_rfl
            · exact absurd hxlt (not_lt.mpr (ht x hx))⟩
      · push_neg at ha
        obtain ⟨Λ, hΛ, hmem, hmin⟩ := ih (by
          obtain ⟨x, hx, hlt⟩ := h
          simp only [List.mem_cons] at hx
          rcases hx with rfl | hx
          · exact absurd hlt (not_lt.mpr ha)
          · exact ⟨x, hx, hlt⟩)
        exact ⟨Λ, hΛ, by simp [hmem], by
          intro x hx hxlt
          simp only [List.mem_cons] at hx
          rcases hx with rfl | hx
          · exact absurd hxlt (not_lt.mpr ha)
          · exact hmin x hx hxlt⟩

/-- 窗口算术：`0 < m1`、`m1 < 1` ⇒ `m1² < m1`（无 `norm_num`/`nlinarith`，手证）。 -/
theorem sq_lt_self {m1 : ℝ} (h0 : 0 < m1) (h1 : m1 < 1) : m1 ^ 2 < m1 := by
  rw [pow_two]
  calc m1 * m1 < m1 * 1 := mul_lt_mul_of_pos_left h1 h0
    _ = m1 := mul_one m1

/-! ### 命题 C 方向 1（无条件）：`I2 ⇒ 次小谱点` -/

theorem propC_fwd (S : Spec) (m1 : ℝ) (h : I2 S m1) : SecExists S m1 := by
  obtain ⟨Λ, ⟨hlt, hmem, _⟩, _⟩ := h
  exact ⟨Λ, hmem, hlt⟩

/-! ### 命题 C 方向 2（需窗口 `0 < m1 < 1`）：`次小谱点 ⇒ I2` -/

theorem propC_bwd (S : Spec) (m1 : ℝ) (h0 : 0 < m1) (h1 : m1 < 1)
    (h : SecExists S m1) : I2 S m1 := by
  obtain ⟨Λ, hΛlt, hmem, hmin⟩ := exists_min_above S.pts m1 h
  exact ⟨Λ, ⟨hΛlt, hmem, hmin⟩, lt_trans (sq_lt_self h0 h1) hΛlt⟩

/-- **命题 C（核）**：`I2 ⟺ 次小谱点`，窗口 `0 < m1 < 1`。
    右端不含 `I1`：若 `pts` 就是正权重支撑，则 `I1` 自动成立（`I1_automatic`）。 -/
theorem propC (S : Spec) (m1 : ℝ) (h0 : 0 < m1) (h1 : m1 < 1) :
    I2 S m1 ↔ SecExists S m1 :=
  ⟨propC_fwd S m1, propC_bwd S m1 h0 h1⟩

/-- `I1` 在「`pts` = 正权重支撑」下自动成立。 -/
theorem I1_automatic (S : Spec) (hpos : ∀ x ∈ S.pts, 0 < S.wt x)
    (m1 : ℝ) (hmem : m1 ∈ S.pts) : I1 S m1 :=
  hpos m1 hmem

/-- **命题 C 的字面形式**（`I1` 留在右端）：在正支撑假设下与 `propC` 等价。 -/
theorem propC_literal (S : Spec) (hpos : ∀ x ∈ S.pts, 0 < S.wt x)
    (m1 : ℝ) (hmem : m1 ∈ S.pts) (h0 : 0 < m1) (h1 : m1 < 1) :
    I2 S m1 ↔ I1 S m1 ∧ SecExists S m1 := by
  rw [propC S m1 h0 h1]
  exact ⟨fun h => ⟨I1_automatic S hpos m1 hmem, h⟩, fun h => h.2⟩

/-! ### 反例：若 `pts` 是「网格点表」而非正支撑，`I2` **不蕴含** `I1`

⇒ 字面 `I2 ⟺ I1 ∧ ∃次小谱点` 的 `I1` 一支不能从 `I2` 推出。 -/

/-- 点表 `[0, 1]`、权重全压在 `1` 上、`0` 点权重为 0。 -/
noncomputable def demoWt : ℝ → ℝ := fun x => if x = 0 then 0 else 1

noncomputable def demoSpec : Spec where
  pts := [0, 1]
  hne := List.cons_ne_nil 0 [1]
  wt := demoWt
  hnn := by
    intro x
    unfold demoWt
    split_ifs
    · exact le_rfl
    · exact zero_le_one

/-- 反例：`I2 demoSpec 0` 成立（`Λ★ = 1 > 0² = 0`），但 `I1 demoSpec 0` 不成立。 -/
theorem counterexample_I2_not_I1 : I2 demoSpec 0 ∧ ¬ I1 demoSpec 0 := by
  have hstar : LtStar demoSpec 0 1 := by
    refine ⟨zero_lt_one, ?_, ?_⟩
    · show (1 : ℝ) ∈ [0, 1]
      simp
    · intro x hx hxlt
      change x ∈ [0, 1] at hx
      simp only [List.mem_cons] at hx
      rcases hx with h0 | h1
      · subst h0
        exact absurd hxlt (not_lt.mpr le_rfl)
      · rcases h1 with h1 | h2
        · subst h1
          exact le_rfl
        · simp at h2
  refine ⟨⟨1, hstar, ?_⟩, ?_⟩
  · rw [pow_two, zero_mul]; exact zero_lt_one
  · have h0 : demoSpec.wt (0 : ℝ) = 0 := by
      change demoWt 0 = 0
      unfold demoWt
      rw [if_pos rfl]
    show ¬ (0 < demoSpec.wt 0)
    rw [h0]
    exact not_lt.mpr le_rfl

end L04

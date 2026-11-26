/-
  L3B — Assembly:  I2 ∧ (★) ⇒ m > 0   （有限谱模型，sorry-free）

  I2      : Λ★ − m₁² ≥ δ > 0    （0⁺⁺ 道隔离：谱上除基态 m₀ 外全部 ≥ m₀ + δ）
  (★)     : 其他道的最轻质量 ≥ c > 0
  结论    : ∃ g > 0, 谱上除 m₀ 外的所有点 ≥ g （全谱隙）

  约定：`spec : Finset ℝ` 是有限谱（物理上是平方质量谱 m²），`m₀` 是 0⁺⁺ 基态质量，
  `m₀ ≥ 0` 只在一处用到：把 `x ≥ m₀ + δ` 降成 `x ≥ δ`。见 `assembly_of_I2`。
-/

import Mathlib.Data.Real.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Order.MinMax
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic.Linarith

set_option autoImplicit false

namespace MassGap.L3B

/-! ## 1. 主装配定理（用户给定的裸签名） -/

/-- **L3B 主定理**：I2 与 (★) 一起给出全谱隙，显式常数 `g := min δ c`。 -/
theorem assembly (spec : Finset ℝ) (m₀ δ c : ℝ)
    (hm : 0 ≤ m₀) (hδ : 0 < δ) (hc : 0 < c)
    (h1 : ∀ x ∈ spec, x = m₀ ∨ x ≥ m₀ + δ)
    (hstar : ∀ y ∈ spec, y ≠ m₀ → y ≥ c) :
    ∃ g > 0, ∀ x ∈ spec, x ≠ m₀ → x ≥ g := by
  refine ⟨min δ c, lt_min hδ hc, ?_⟩
  intro x hx hne
  rcases h1 x hx with h | h
  · exact absurd h hne
  · exact le_trans (min_le_left δ c) (le_trans (by linarith) h)

/-! ## 2. 诚实备注：本形式下 (★) 是冗余的 -/

/-- **只用 I2**：`g := δ` 已经够了。也就是说第 1 节的 `min δ c` 没有用到 (★)。
    (★) 真正不可替代的场合是 `m₀` 未知/未定为非负的情况（见第 4 节）。 -/
theorem assembly_of_I2 (spec : Finset ℝ) (m₀ δ : ℝ)
    (hm : 0 ≤ m₀) (hδ : 0 < δ)
    (h1 : ∀ x ∈ spec, x = m₀ ∨ x ≥ m₀ + δ) :
    ∃ g > 0, ∀ x ∈ spec, x ≠ m₀ → x ≥ g :=
  ⟨δ, hδ, fun x hx hne => by
    rcases h1 x hx with h | h
    · exact absurd h hne
    · linarith⟩

/-! ## 3. 两扇区版本（(★) 在这里才真的出力） -/

/-- I2 只作用于 0⁺⁺ 扇区 `A`，(★) 作用于其它道扇区 `B`，两者并起来才是全谱。 -/
theorem assembly_two_sector (A B : Finset ℝ) (m₀ δ c : ℝ)
    (hm : 0 ≤ m₀) (hδ : 0 < δ) (hc : 0 < c)
    (h1 : ∀ x ∈ A, x = m₀ ∨ x ≥ m₀ + δ)
    (hstar : ∀ y ∈ B, y ≥ c) :
    ∃ g > 0, ∀ x ∈ A ∪ B, x ≠ m₀ → x ≥ g := by
  refine ⟨min δ c, lt_min hδ hc, ?_⟩
  intro x hx hne
  rcases Finset.mem_union.mp hx with hxA | hxB
  · rcases h1 x hxA with h | h
    · exact absurd h hne
    · exact le_trans (min_le_left δ c) (le_trans (by linarith) h)
  · exact le_trans (min_le_right δ c) (hstar x hxB)

/-! ## 4. 若不假定 `m₀ ≥ 0`：(★) 独自就给出绝对地板 -/

/-- 去掉 `0 ≤ m₀` 后 `min δ c` 不再自动为正，此时 (★) 的 `c` 接管地板。 -/
theorem assembly_no_ground_sign (spec : Finset ℝ) (m₀ c : ℝ) (hc : 0 < c)
    (hstar : ∀ y ∈ spec, y ≠ m₀ → y ≥ c) :
    ∃ g > 0, ∀ x ∈ spec, x ≠ m₀ → x ≥ g :=
  ⟨c, hc, fun x hx hne => hstar x hx hne⟩

/-! ## 5. 从谱隙到质量隙：`m := √g > 0` -/

/-- 平方质量谱上的隙 `g > 0` 直接给出质量 `m > 0`，`m² ≤ x`。 -/
theorem assembly_mass (spec : Finset ℝ) (m₀ δ c : ℝ)
    (hm : 0 ≤ m₀) (hδ : 0 < δ) (hc : 0 < c)
    (h1 : ∀ x ∈ spec, x = m₀ ∨ x ≥ m₀ + δ)
    (hstar : ∀ y ∈ spec, y ≠ m₀ → y ≥ c) :
    ∃ m > 0, ∀ x ∈ spec, x ≠ m₀ → m ^ 2 ≤ x := by
  obtain ⟨g, hg, hb⟩ := assembly spec m₀ δ c hm hδ hc h1 hstar
  refine ⟨Real.sqrt g, Real.sqrt_pos.mpr hg, ?_⟩
  intro x hx hne
  rw [Real.sq_sqrt (le_of_lt hg)]
  exact hb x hx hne

/-! ## 6. `Set ℝ` 版本：不需要有限性/离散性 -/

/-- 同一陈述在任意 `Set ℝ` 上成立（证明逐字相同）：
    `∃ g > 0, ∀ x ∈ spec, x ≠ m₀ → x ≥ g` 只是"一致下界"，不需要谱有限。 -/
theorem assembly_set (spec : Set ℝ) (m₀ δ c : ℝ)
    (hm : 0 ≤ m₀) (hδ : 0 < δ) (hc : 0 < c)
    (h1 : ∀ x ∈ spec, x = m₀ ∨ x ≥ m₀ + δ)
    (hstar : ∀ y ∈ spec, y ≠ m₀ → y ≥ c) :
    ∃ g > 0, ∀ x ∈ spec, x ≠ m₀ → x ≥ g := by
  refine ⟨min δ c, lt_min hδ hc, ?_⟩
  intro x hx hne
  rcases h1 x hx with h | h
  · exact absurd h hne
  · exact le_trans (min_le_left δ c) (le_trans (by linarith) h)

/-! ## 7. 非空性：假设体系可满足（定理不 vacuous） -/

/-- 显式见证：`spec = {0, 1}`, `m₀ = 0`, `δ = 1`, `c = 1`。 -/
theorem witness_satisfiable :
    ∃ spec : Finset ℝ, ∃ m₀ δ c : ℝ,
      0 ≤ m₀ ∧ 0 < δ ∧ 0 < c ∧
        (∀ x ∈ spec, x = m₀ ∨ x ≥ m₀ + δ) ∧
        (∀ y ∈ spec, y ≠ m₀ → y ≥ c) := by
  classical
  refine ⟨{0, 1}, 0, 1, 1, le_refl 0, one_pos, one_pos, ?_, ?_⟩
  · intro x hx
    simp only [Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with h | h
    · subst h; exact Or.inl rfl
    · subst h; exact Or.inr (by linarith)
  · intro y hy hne
    simp only [Finset.mem_insert, Finset.mem_singleton] at hy
    rcases hy with h | h
    · subst h; exact absurd rfl hne
    · subst h; linarith

/-! ## 7. 公理审计 -/

#print axioms MassGap.L3B.witness_satisfiable
#print axioms MassGap.L3B.assembly
#print axioms MassGap.L3B.assembly_of_I2
#print axioms MassGap.L3B.assembly_two_sector
#print axioms MassGap.L3B.assembly_set
#print axioms MassGap.L3B.assembly_no_ground_sign
#print axioms MassGap.L3B.assembly_mass

end MassGap.L3B

/-
M03 — 基准反例（有限维 / 2×2）

目的（环境裁决 + 数学内容一并落盘）：
  证明 "残数型量 f₁² 消失 ⇒ 无隙 (m = 0)" 这一步是**非法推理**。
  在最省的有限维模型里造一个显式反例：反交换子矩阵元为零，而谱隙 m > 0。

物理对位（2D YM 型链条）：
  trF² = H 守恒  ⇒  [H, O] = 0
              ⇒  ⟨𝔾|[H,O]|n⟩ = 0
              ⇒  (H 的本征值不退化时)  f₁² ∝ O_{01} = 0
  但 **不等于** m = 0 —— 本文件给出同时满足两者（O 非对角元 = 0 且 谱隙 > 0）的
  2×2 见证。故 "f₁² = 0 ⇒ m = 0" 在这一类模型上被证伪（[T]，Lean 内核检查）。

编译方式（已实测，本仓唯一可用路径，见 M03-LEANPATH.md）：
  cd vendor\OSforGFF ; lake env lean ..\..\research\v3\lean\M03-COUNTER.lean
-/
import Mathlib.Data.Matrix.Basic
import Mathlib.Tactic.NormNum

namespace M03

open Matrix

/-! ### A 层：算术核（无量词技巧，纯 `Rat` 域运算）

`⟨𝔾|[H,O]|n⟩` 在 H 对角、谱隙 `a ≠ b` 时等于 `(a - b) * O_{01}`；
该项为零且 `a ≠ b` 强制 `O_{01} = 0`。这是链条里唯一有内容的算术步。 -/

/-- 交换子矩阵元为零 + 谱隙非退化 ⇒ 非对角元为零。 -/
theorem f1sq_eq_zero (a b o01 : ℚ) (hne : a ≠ b) (hcomm : a * o01 = o01 * b) :
    o01 = 0 := by
  have h : (a - b) * o01 = 0 := by
    rw [sub_mul, hcomm, mul_comm b o01, sub_self]
  exact (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr hne)

/-- 谱隙为正 ⇒ 本征值不退化（把 `a ≠ b` 换成 `a < b` 的接口）。 -/
theorem ne_of_lt_q {a b : ℚ} (h : a < b) : a ≠ b := fun he => by
  rw [he] at h
  exact lt_irrefl b h

/-! ### B 层：2×2 显式见证

H = diag(0, 1)（谱隙 = 1 > 0），O = diag(3, -1)（非标量，非对角元 = 0）。
两者都是对角矩阵 ⇒ 交换；隙为正；残数项 `O₀₁` 为零。 -/

/-- 有隙哈密顿量。 -/
def H : Matrix (Fin 2) (Fin 2) ℚ := diagonal ![0, 1]

/-- 与 H 交换的可观测量（trF² 型守恒量的替身）。 -/
def O : Matrix (Fin 2) (Fin 2) ℚ := diagonal ![3, -1]

/-- 残数型量 `f₁²` 的替身：非对角矩阵元。 -/
def f1sq : ℚ := O 0 1

/-- `trF² = H` 守恒的对位：`[H, O] = 0`。 -/
theorem comm : H * O = O * H := by simp [H, O, Matrix.diagonal_mul_diagonal]

/-- 矩阵元为零：`⟨𝔾|[H,O]|n⟩ = 0 ⇒ f₁² = 0`。 -/
theorem f1sq_zero : f1sq = 0 := by decide

/-- 但仍**有**谱隙：`m = H₁₁ - H₀₀ = 1 > 0`。 -/
theorem gap_pos : (0 : ℚ) < H 1 1 - H 0 0 := by
  rw [show H 1 1 = 1 from by decide, show H 0 0 = 0 from by decide]
  norm_num

/-- H 的两条本征值不退化。 -/
theorem spectrum_nondegenerate : H 0 0 ≠ H 1 1 := by
  rw [show H 1 1 = 1 from by decide, show H 0 0 = 0 from by decide]
  decide

/-- **基准反例**：存在交换的一对 (H, O)，其残数项为零，而谱隙严格为正。 -/
theorem counterexample :
    ∃ H O : Matrix (Fin 2) (Fin 2) ℚ,
      H * O = O * H ∧ O 0 1 = 0 ∧ (0 : ℚ) < H 1 1 - H 0 0 :=
  ⟨H, O, comm, f1sq_zero, gap_pos⟩

/-- 反例的否定式表述：`f₁² = 0` 不能推出 `m = 0`。 -/
theorem not_f1sq_zero_implies_gapless :
    ¬ (∀ H O : Matrix (Fin 2) (Fin 2) ℚ,
        H * O = O * H → O 0 1 = 0 → H 1 1 - H 0 0 = 0) := by
  intro h
  have hz : H 1 1 - H 0 0 = 0 := h H O comm f1sq_zero
  have hg := gap_pos
  rw [hz] at hg
  exact lt_irrefl (0 : ℚ) hg

/-! ### 自检：无 `sorry` -/

#print axioms M03.counterexample
#print axioms M03.not_f1sq_zero_implies_gapless

end M03


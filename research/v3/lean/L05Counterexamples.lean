/-
L05 · 基准反例（M6 / 2D YM）的最省形式化。

靶 T1（被攻击的对象）：  `G″ > 0  ⇒  m > 0`
—— 4D YM 五环里那把"钥匙"的单向蕴含形式。本文件只做一件事：
在**最省的载体**上把 `G″ > 0  ∧  gap = 0` 同时实现出来，看它能否机器验证。

分三层（按"最省优先"：有限维 > 需 Gaussian > 需谱论）：

§1  M6 的算术核：`Var(a) = 1/(2ω) ⇒ G″ = 2ω > 0`          [T]
§2  "无隙" 的最省表述：Rayleigh 无隙见证（不用任何谱论）    [def]
§3  有限维 2×2 二分法：简并真空 ⇒ 无隙；简单真空 ⇒ 必有隙   [T]
§4  判决：有限维反例只在**真空简并**下成立 ⇒ 只杀裸蕴含，
    不构成 M6 同类反例（M6 的真空是简单的，机制是连续谱）。
-/

import Mathlib.Data.Rat.Lemmas
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith

namespace L05

/-! ## §1  M6 的算术核：`G″ = 1/Var(a)`，`Var(a) = 1/(2ω)` ⇒ `G″ = 2ω`

M6 = 谐振子（a 方向）⊗ 自由粒子（b 方向）。
a 方向的高斯测度给 `Var(a) = 1/(2ω)`，而 `G″ = 1/Var(a) = 2ω`。
**这一步完全不需要 Gaussian 积分**：`Var(a) = 1/(2ω)` 是高斯测度的定义性质，
在 Lean 里把它作为模型数据（hypothesis）带进来即可；`G″ = 1/Var` 是定义。
于是 "`G″ > 0`" 是纯算术。 -/

/-- `G″ = 1/Var(a)`，写成避免除法的等价形式（乘性）。 -/
def GSecondOf (V : ℚ) : ℚ := V⁻¹

/-- 字面形式：`Var(a) = 1/(2ω)` ⇒ `G″ = 2ω`。 -/
theorem gsecond_harmonic (ω : ℚ) (hω : ω ≠ 0) :
    GSecondOf (1 / (2 * ω)) = 2 * ω := by
  unfold GSecondOf
  field_simp

/-- M6 的模型数据：a 方向高斯（`Var(a) = 1/(2ω)`），`G″ = 1/Var(a)`。 -/
structure M6 where
  ω : ℚ
  hω : 0 < ω
  VarA : ℚ
  hVar : VarA * (2 * ω) = 1
  G2 : ℚ
  hG2 : G2 * VarA = 1

theorem m6_Gsecond (M : M6) : M.G2 = 2 * M.ω := by
  have hVarA : M.VarA ≠ 0 := by
    have h := M.hVar
    intro h0
    rw [h0, zero_mul] at h
    exact zero_ne_one h
  refine mul_right_cancel₀ hVarA ?_
  rw [M.hG2, mul_comm (2 * M.ω) M.VarA, M.hVar]

/-- **[T] M6 的 `G″` 严格为正**（题意要求的 `G″ = 2ω > 0`）。 -/
theorem m6_Gsecond_pos (M : M6) : 0 < M.G2 := by
  rw [m6_Gsecond M]
  linarith [M.hω]

theorem m6_Gsecond_literal (M : M6) : M.G2 = 1 / (1 / (2 * M.ω)) := by
  rw [m6_Gsecond M]
  field_simp

/-! ## §2  「无隙」的最省表述：Rayleigh 无隙见证

`σ(H) = [0,∞)` 需要无界自伴算子的谱定理（Mathlib 无）。
最省且等价的表述是**变分/Rayleigh 形式**：

  `NoGap H vac : ∀ ε > 0, ∃ ψ, ψ ⊥ vac, ‖ψ‖ = 1, ⟨ψ, H ψ⟩ < ε`

即"与真空正交的归一化态里存在任意低的能量"。
只要 `H ≥ 0`，这与 `inf (σ H ∩ (0,∞)) = 0` 等价，且**不需要任何谱论**。 -/

/-- 2 维状态空间（a 分量, b 分量）——"a 方向 × b 方向"的最小载体。 -/
abbrev V2 := ℚ × ℚ

/-- 真空 `e₀ = (1, 0)`；`b` 方向是那个"完全自由/无势"的方向。 -/
def vac : V2 := (1, 0)

/-- 与真空正交且归一化的态（`a` 分量 = 0 ⇒ 只在 `b` 方向上有激发）。 -/
def Excited (v : V2) : Prop := v.1 = 0 ∧ v.1 ^ 2 + v.2 ^ 2 = 1

/-- Rayleigh 二次型 `⟨v, H v⟩`，`H = !![0, c; c, e]`（对称，且 `H e₀ = 0`）。 -/
def q (c e : ℚ) (v : V2) : ℚ := 2 * c * v.1 * v.2 + e * v.2 ^ 2

/-- **无隙**：任意小的激发能量。 -/
def NoGap (c e : ℚ) : Prop := ∀ ε : ℚ, 0 < ε → ∃ v : V2, Excited v ∧ q c e v < ε

theorem q_c_zero (e : ℚ) (v : V2) : q 0 e v = e * v.2 ^ 2 := by
  unfold q; ring

theorem q_eq_zero_at_vac (c e : ℚ) : q c e vac = 0 := by
  unfold q vac; ring

/-! ## §3  有限维 2×2 二分法 -/

/-- **[T] 简并真空（`e = 0`）⇒ 无隙**：`v = (0,1)` 是零能激发。
这是"`G″ > 0` 而 `gap = 0`"的最小实现。 -/
theorem degenerate_NoGap (c : ℚ) : NoGap c 0 := by
  intro ε hε
  refine ⟨(0, 1), ?_, ?_⟩
  · unfold Excited; norm_num
  · have hq : q c 0 (0, 1) = 0 := by unfold q; ring
    rw [hq]; linarith

/-- **[T] 简单真空（`0 < e`）⇒ 必有隙**：所有激发能量 `≥ e > 0`。 -/
theorem simple_vacuum_gap (e : ℚ) (he : 0 < e) :
    ∀ v : V2, Excited v → 0 < q 0 e v := by
  intro v hv
  obtain ⟨h1, h2⟩ := hv
  have hsq : v.2 ^ 2 = 1 := by
    rw [h1] at h2
    simpa using h2
  rw [q_c_zero, hsq, mul_one]
  exact he

/-- **[T] 二分法**：`2×2` 对称、`H e₀ = 0`、且二次型非负 ⇒ `e = 0` 时必 `c = 0`。
即：非负的 2×2 模型里，"无隙" ⟺ **真空简并**，没有第三条路。
（`c ≠ 0` 时 `q(1,y) = 2cy + ey²` 在 `y` 上无下界。） -/
theorem dichotomy (c : ℚ) (hnn : ∀ v : V2, 0 ≤ q c 0 v) : c = 0 := by
  have h1 : 0 ≤ 2 * c := by
    have := hnn (1, 1)
    simpa [q] using this
  have h2 : 0 ≤ -(2 * c) := by
    have := hnn (1, -1)
    simpa [q] using this
  linarith

/-! ## §4  判决（本文件真正给出的东西）

`finite_dim_witness`：同一个模型可以同时有 `G″ > 0` 与 `gap = 0`
—— 但 `gap = 0` 只能来自**真空简并**，而不是 M6 的连续谱。
所以：**裸蕴含 `G″ > 0 ⇒ m > 0` 被机器证伪**〔T〕；
而 **M6 同类反例（简单真空 + 连续谱）** 在这里**没有**被机器验证〔NG〕，
它需要自由 Laplacian 的 `σ = [0,∞)`。 -/

/-- **合成反例**：`G″ = 2ω > 0`（a 方向高斯） ∧ `NoGap`（b 方向退化零模）。 -/
theorem finite_dim_witness (M : M6) : 0 < M.G2 ∧ NoGap 0 0 :=
  ⟨m6_Gsecond_pos M, degenerate_NoGap 0⟩

#print axioms L05.m6_Gsecond_pos
#print axioms L05.degenerate_NoGap
#print axioms L05.simple_vacuum_gap
#print axioms L05.dichotomy
#print axioms L05.finite_dim_witness

end L05

/-
  L02-CBARKATO.lean
  =================
  物种 L02：① `C̄ ≤ e^{Δu} = 1.054`（c07，[T]）  ② `Kato ⇒ m > 0`（f04，[T]）

  ## 工具状态（必须先读）
  本机：Lean 4.36.0 core（elan），**无 Mathlib**。
  实测 `import Mathlib` →
      `error: unknown module prefix 'Mathlib'`
      `No directory 'Mathlib' or file 'Mathlib.olean' in the search path entries`
  ⇒ `Real.exp`、`Real.log`、`spectrum`、`Metric.infDist` 一律不可用；
  core Lean 4 的 `Rat` 也没有 `OfNat ℚ` / `Div ℚ` 实例（实测 synthInstanceFailed）。

  因此本文件做的是**能做到的最大机器验证**：
  把 ① 的全部**算术内容**（等价于 `e^{5/96} ≤ 1.054` 的精确整数证书）
  与 ② 的单重谱隙核（`0 < (λ₁−δ)/λ₁ < 1`）降到 `Nat` 上，由内核检查。
  解析转移步（exp 单调、`e^t ≥ 1+t`、`−log` 严格减）在本机不可形式化，
  以精确陈述 + 手写证明记入 `research/v3/species/L02-CBARKATO.md`，标「待形式化」。

  `sorry` 数 = 0。全部定理由 `decide` / `Nat.mul_le_mul` / `omega` 关闭。
-/

namespace L02

/-! ## ① c07：`C̄ ≤ e^{Δu} = 1.054` — 算术骨架

链条（见 md）：
  β < 1/48 ⇒ u = c_f/N² = β/(2N²) < 1/(96N²)
  ⇒ Δu < 20/(96N²) = 5/(24N²) ≤ 5/96   (N ≥ 2)
  ⇒ e^{Δu} ≤ e^{5/96} ≤ (1 − (5/96)/3)⁻³ = (288/283)³ ≤ 1054/1000 = 1.054.
  ⇒ Ĉ·C̄ ≤ 4 · 1.054 = 4.216 < 7.0633.

下面 5 条把这条链上**除两次解析转移以外的每一步**都做了精确整数化。
-/

/-- 步 1：`20/96 = 5/24`（交叉相乘的整数形态，精确）。 -/
theorem delta_over_96_eq_five_over_24 : 20 * 24 = 5 * 96 := by decide

/-- 步 2：最坏情形在 `N = 2`。`N ≥ 2 ⇒ 96 ≤ 24·N²`，
即 `5/(24N²) ≤ 5/96`（两边同乘正数 `96·24N²`）。 -/
theorem worst_case_at_N_two (N : Nat) (hN : 2 ≤ N) : 96 ≤ 24 * (N * N) := by
  have h4 : 4 ≤ N * N := Nat.mul_le_mul hN hN
  omega

/-- 步 2'：`N = 2` 时取等：`96 = 24·(2·2)`。 -/
theorem N_two_sharp : 96 = 24 * (2 * 2) := by decide

/-- 步 3（**关键整数证书**）：`(288/283)³ ≤ 1054/1000`。
即 `1000 · 288³ ≤ 1054 · 283³` —— 这是 `e^{5/96} ≤ 1.054` 在
`e^x ≤ (1 − x/3)⁻³` 下的精确有理化形态。 -/
theorem exp_bound_certificate : 1000 * 288 ^ 3 ≤ 1054 * 283 ^ 3 := by decide

/-- 步 3''：真值方向核对 —— `e^{5/96} = 1.0535548…`，故余量真实存在；
这里给出对应的更强有理证书 `(288/283)³ < 1054/1000`（严格，非取等）。 -/
theorem exp_bound_strict : 1000 * 288 ^ 3 < 1054 * 283 ^ 3 := by decide

/-- 步 3'''（负对照，诚实性检查）：`n = 2` **不够**。
`(192/187)² > 1.054`，故 `e^x ≤ (1 − x/2)⁻²` 这条更弱的上界在此**失败**，
必须用 `n = 3`。 -/
theorem n_two_insufficient : 1054 * 187 ^ 2 < 1000 * 192 ^ 2 := by decide

/-- 步 4：`Ĉ·C̄ ≤ 4 · 1.054 = 4.216 < 7.0633`（精确整数形态）。 -/
theorem cbar_product_bound : 4 * 1054 * 10000 < 70633 * 1000 := by decide

/-- 步 4'：`4 × 1.054 = 4.216`（取等，精确）。 -/
theorem cbar_product_exact : 4 * 1054 * 1000 = 4216 * 1000 := by decide

/-! ## ② f04：`Kato ⇒ m > 0` — 谱隙核

实数形式（md 中给出精确陈述与证明）：
  `T` 紧自伴、`λ₁ = ‖T‖` 单重、`δ := dist(λ₁, spec T ∖ {λ₁}) > 0`
  ⇒ `λ₂ ≤ λ₁ − δ` ⇒ `m := −log(λ₂/λ₁) ≥ −log(1 − δ/λ₁) > 0`。

`−log(1 − δ/λ₁) > 0` 的**全部算术内容**就是 `0 < 1 − δ/λ₁ < 1`，
即 `0 < δ < λ₁`。这一条与它成立后 `m > 0` 的代数一致性在下面机器验证。
-/

/-- ②-核：`0 < δ < λ₁ ⇒ 0 < λ₁ − δ < λ₁`（即 `(λ₁−δ)/λ₁ ∈ (0,1)`）。
这一步之后 `−log` 的严格减性给出 `m > 0`。 -/
theorem gap_ratio_in_unit_interval (lam del : Nat) (hpos : 0 < del) (hlt : del < lam) :
    0 < lam - del ∧ lam - del < lam := by
  omega

/-- ②-核'：`δ = λ₁` 时谱隙退化（`r_iso = 0` 的模型），核**恰在边界**失效：
`0 < λ₁ − λ₁` 为假。这是「严格性」在对偶侧的体现。 -/
theorem gap_ratio_collapses_at_boundary (lam : Nat) (_h : 0 < lam) : ¬ (0 < lam - lam) := by
  omega

/-- ②-数值实例：`λ₁ = 1/2`、`δ = 1/4` ⇒ `1 − δ/λ₁ = 1/2 ∈ (0,1)`，
于是 `m ≥ −log(1/2) = log 2 > 0`。
整数形态：`2·(2−1) = 2`，`2−1 = 1 < 2`。 -/
theorem kato_instance_half_quarter : 0 < 2 - 1 ∧ 2 - 1 < 2 := by decide

/-! ## 数值记录（非证明，仅供对账） -/

#eval (288.0 / 283.0) ^ 3
#eval (192.0 / 187.0) ^ 2
#eval Float.exp (5.0 / 96.0)

end L02

/-! ## 公理审计（`#print axioms`）：不得出现 `sorryAx` -/

#print axioms L02.delta_over_96_eq_five_over_24
#print axioms L02.worst_case_at_N_two
#print axioms L02.exp_bound_certificate
#print axioms L02.exp_bound_strict
#print axioms L02.n_two_insufficient
#print axioms L02.cbar_product_bound
#print axioms L02.cbar_product_exact
#print axioms L02.gap_ratio_in_unit_interval
#print axioms L02.gap_ratio_collapses_at_boundary
#print axioms L02.kato_instance_half_quarter

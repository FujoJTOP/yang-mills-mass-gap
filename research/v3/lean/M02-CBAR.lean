/-
  M02-CBAR.lean
  =============
  物种 M02：`C̄ ≡ N` 的**算术核**（mathlib 形式化）。口径逐字取 `f07` 的裁决
  （`research/v3/species/F07-ADJCBAR.md` §0 行 3、行 5、行 7）。

  ## 环境（实测，非据称）
  * `lean 4.33.0-rc1`（`C:\Users\hooya\.elan\toolchains\leanprover--lean4---v4.33.0-rc1\bin\lean.exe`）
    —— mathlib 本体就是用这个 toolchain 编的；`v4.36.0` 读它的 `.olean` 会因版本不符失败。
  * `mathlib v4.33.0-rc1`（rev `79d0395a`）的**完整**构建在
    `vendor\OSforGFF\.lake\packages\mathlib\.lake\build.badcache\lib\lean`（8279 个 `.olean`，
    含根 `Mathlib.olean`）。同一包的 `...\.lake\build\lib\lean` 只有 2136 个 `.olean`
    且**无根 `Mathlib.olean`** —— `import Mathlib` 在那里必然失败。
  * `LEAN_PATH` 必须**显式包含 toolchain 的 `lib\lean`**（`LEAN_PATH` 一旦设置就顶替默认
    搜索路径），否则 `import Mathlib` 在 `Std.Http.Data.Body.Basic` 上报
    `failed to read file …Basic.olean`。

  跑法与实测输出见 `research/v3/species/M02-LEANC BAR.md`。

  ## 命题（`f07` 口径）
  ```
  |w|^{1/n}/u = c_f/u = N            （cube ∏mult = 1）
  => Ĉ·C̄ = 4N                        （Ĉ = 4 包络）
     N = 2  => 8                     而阈值 2N²·48/eΔ = 19.2/e = 7.0632853…
     窗口  β < 2N²/(eΔ·C̄)            N=2 时 = 1/(5e) = 3.5316426×(1/48)
  ```

  `sorry` 数 = 0。全部由 `norm_num` / `field_simp` / `ring` / `nlinarith` /
  `Real.exp_one_gt_d9` / `Real.exp_one_lt_d9` 关闭；`#print axioms` 见文件末（无 `sorryAx`）。
-/

import Mathlib

namespace M02

/-! ## 1 · 归一：`c_f := β/(2N)`、`u := c_f/N`（`f07 §1.4`；`Z06-HOLEA` 的 `u = c_f/N`） -/

/-- `c_f = β/(2N)`。`SU(2)`、`β = 1/12` 时 `= 1/48 = 0.0208333…`；而 `f07` 实测
`⟨(1/N)tr U_p⟩ = 0.0208273`（差 0.03%）。 -/
noncomputable def cf (β : ℝ) (N : ℕ) : ℝ := β / (2 * N)

/-- `u = c_f/N = β/(2N²)`，即 `f07` / `c03` 用的那个归一。 -/
noncomputable def u (β : ℝ) (N : ℕ) : ℝ := cf β N / N

/-- `SU(2)` 的 `c_f`：`β = 1/12 ⇒ c_f = 1/48`（精确有理数）。 -/
theorem cf_SUN2 : cf (1 / 12) 2 = 1 / 48 := by
  unfold cf; norm_num

/-- `SU(2)` 的 `u`：`β = 1/12 ⇒ u = 1/96`。 -/
theorem u_SUN2 : u (1 / 12) 2 = 1 / 96 := by
  unfold u cf; norm_num

/-! ## 2 · `C̄ ≡ N`：定义恒等（不是估计） -/

/-- **`c_f/u = N`** —— 一般 `N`、一般 `β ≠ 0`。这是 `C̄` 的**值**：`C̄` 不是 `1`，
也不是 `Δ = 20`，而是 `N`（`f07` §0 行 5，`[T]`）。 -/
theorem cbar_eq_N (β : ℝ) (N : ℕ) (hβ : β ≠ 0) (hN : (N : ℝ) ≠ 0) :
    cf β N / u β N = (N : ℝ) := by
  simp only [u, cf]
  field_simp

/-- `N = 2`：`C̄ = c_f/u = 2`。 -/
theorem cbar_eq_N_two (β : ℝ) (hβ : β ≠ 0) : cf β 2 / u β 2 = 2 := by
  simpa using cbar_eq_N β 2 hβ (by norm_num)

/-- `Φ(0)/u = 2N²`（`Φ(0) = Σ_{r≠0} c_r dim r = β = 2N²u`，`f`/`f̄` 各 `β/2`）。
故 `C̄ ∈ [N, 2N²]`，`N = 2` 时 `[2, 8]`。 -/
theorem phi0_over_u (β : ℝ) (N : ℕ) (hβ : β ≠ 0) (hN : (N : ℝ) ≠ 0) :
    (2 * (N : ℝ) ^ 2 * u β N) / u β N = 2 * (N : ℝ) ^ 2 := by
  have hu : u β N ≠ 0 := by
    simp only [u, cf]
    exact div_ne_zero (div_ne_zero hβ (mul_ne_zero (by norm_num) hN)) hN
  field_simp

/-- `N = 2`：`Φ(0)/u = 8` —— 与 `Ĉ·C̄ = 8` 数值相同但**是另一个量**
（`f07` §3 的「同名不同物」：`Δ = 20` 属树计数，`2N²` 属权重比）。 -/
theorem phi0_over_u_two (β : ℝ) (hβ : β ≠ 0) :
    (2 * (2 : ℝ) ^ 2 * u β 2) / u β 2 = 8 := by
  have h := phi0_over_u β 2 hβ (by norm_num)
  norm_num at h ⊢
  exact h

/-! ## 3 · `Ĉ·C̄ = 4N`（`Ĉ = 4` 包络） -/

/-- `Ĉ·C̄ = 4·2 = 8`。 -/
theorem chat_times_cbar : (4 : ℝ) * 2 = 8 := by norm_num

/-- 通式 `Ĉ·C̄ = 4N`；`N = 2` ⇒ `8`。 -/
theorem chat_times_cbar_N (N : ℕ) : (4 : ℝ) * (N : ℝ) = 4 * N := by ring

/-! ## 4 · 阈值 `2N²·48/eΔ`：与 `Ĉ·C̄ = 8` 比较 -/

/-- 阈值的**精确形状**：`2N²·48/(eΔ) = 192/(10e) = 19.2/e`（`N = 2, Δ = 20`）。 -/
theorem thresh_eq :
    2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20) = 192 / (10 * Real.exp 1) := by
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  field_simp
  ring

/-- **`8 > 7.0632853…`**：核心不等式，`sorry`-free。 -/
theorem thresh : (4 : ℝ) * 2 > 2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20) := by
  have he : (2.4 : ℝ) < Real.exp 1 := by linarith [Real.exp_one_gt_d9]
  have hp : (0 : ℝ) < Real.exp 1 * 20 := by positivity
  rw [gt_iff_lt, div_lt_iff₀ hp]
  nlinarith

/-- 该比较**恰好等价于** `e > 2.4`：`e = 12/5` 时两边相等 ⇒ 余量 = `e/(12/5) = 5e/12`。
把「一个数值巧合」升级成「一条等价」：`f07` 的 `7.0633` 不是孤立的数。 -/
theorem thresh_iff :
    ((4 : ℝ) * 2 > 2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20)) ↔ (2.4 : ℝ) < Real.exp 1 := by
  have hp : (0 : ℝ) < Real.exp 1 * 20 := by positivity
  rw [gt_iff_lt, div_lt_iff₀ hp]
  constructor <;> intro h <;> nlinarith [Real.exp_pos 1]

/-- 阈值 = `Ĉ·C̄` 除以 `β = 1/12` 处的 `x = eΔ·N·u`。 -/
theorem thresh_eq_eight_div_x :
    2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20) = 8 / (5 * Real.exp 1 / 12) := by
  norm_num
  field_simp
  ring

/-- `x = eΔ·N·u` 在 `β = 1/12` 处的精确值：`= 5e/12 = 1.1326174…`
（`f07 §1.4` 的 `eΔ·N·u = 1.132617`，逐位一致）。 -/
theorem x_at_beta_twelfth :
    Real.exp 1 * 20 * 2 * (1 / 96 : ℝ) = 5 * Real.exp 1 / 12 := by
  field_simp
  ring

/-- 比值 `Ĉ·C̄ / 阈值 = 5e/12 = 1.1326174…`。 -/
theorem ratio :
    ((4 : ℝ) * 2) / (2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20)) = 5 * Real.exp 1 / 12 := by
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  field_simp
  ring

/-- **`x > 1`**（`f07` §0 行 3 的 ★ 缺口）：`β = 1/12` 处 `eΔ·N·u = 1.1326 > 1`，
故 `24(eΔ·N·u)^α/(1 − eΔ·N·u)` 的几何级数**发散** —— 界在此失效，窗口必须缩到
`β < 1/(5e) = 0.0735759`。 -/
theorem x_gt_one : (1 : ℝ) < ((4 : ℝ) * 2) / (2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20)) := by
  rw [ratio]
  have he : (2.4 : ℝ) < Real.exp 1 := by linarith [Real.exp_one_gt_d9]
  linarith

/-! ## 5 · 窗口 `β < 2N²/(eΔ·C̄)`：`N=2, C̄=2` ⇒ `β < 1/(5e) = 3.5316426×(1/48)`

`f07 §0 行 7` 的窗口是 `0.073576 = 3.532×(1/48)`（无约定口径 `0.0735925`）。
**这里修正任务书给的式子**：`48/(8/(2·2²·48/(e·20))) = 42.37`，不是 `3.53` ——
那个嵌套把「比值」又除了 `48`。正确的窗口比是 `48·2N²/(eΔ·C̄)`。 -/

/-- 窗口通式：`48·2N²/(eΔ·C̄)`，`C̄ = N` ⇒ `48·2N/(20e)`。 -/
theorem window_general (N : ℝ) (hN : N ≠ 0) :
    (48 : ℝ) * (2 * N ^ 2) / (Real.exp 1 * 20 * N) = 48 * (2 * N) / (20 * Real.exp 1) := by
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  field_simp

/-- **窗口比**（`N = 2, C̄ = 2, Δ = 20`）：`48·8/(20e·2) = 48/(5e) = 3.5316426…`。 -/
theorem window :
    (48 : ℝ) * (2 * (2 : ℝ) ^ 2) / (Real.exp 1 * 20 * 2) = 48 / (5 * Real.exp 1) := by
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  field_simp
  ring

/-- 窗口的 `β` 值：`2N²/(eΔ·C̄) = 1/(5e)`。 -/
theorem window_beta :
    2 * (2 : ℝ) ^ 2 / (Real.exp 1 * 20 * 2) = 1 / (5 * Real.exp 1) := by
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  field_simp
  ring

/-- 窗口比落在 `(3.53, 3.54)`：`3.5316426…`。这就是任务书写下的 `3.53`。 -/
theorem window_bounds :
    (3.53 : ℝ) < 48 / (5 * Real.exp 1) ∧ 48 / (5 * Real.exp 1) < 3.54 := by
  have h1 : (2.7182818283 : ℝ) < Real.exp 1 := Real.exp_one_gt_d9
  have h2 : Real.exp 1 < 2.7182818286 := Real.exp_one_lt_d9
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  constructor
  · rw [lt_div_iff₀ (by positivity : (0 : ℝ) < 5 * Real.exp 1)]
    nlinarith
  · rw [div_lt_iff₀ (by positivity : (0 : ℝ) < 5 * Real.exp 1)]
    nlinarith

/-- 窗口 `β < 0.073576`（`f07` 的 `0.073576 = 3.532×(1/48)`），且 `> 0.073575`。 -/
theorem window_beta_bounds :
    (0.073575 : ℝ) < 1 / (5 * Real.exp 1) ∧ 1 / (5 * Real.exp 1) < 0.073576 := by
  have h1 : (2.7182818283 : ℝ) < Real.exp 1 := Real.exp_one_gt_d9
  have h2 : Real.exp 1 < 2.7182818286 := Real.exp_one_lt_d9
  have hp : (0 : ℝ) < Real.exp 1 := Real.exp_pos 1
  constructor
  · rw [lt_div_iff₀ (by positivity : (0 : ℝ) < 5 * Real.exp 1)]
    nlinarith
  · rw [div_lt_iff₀ (by positivity : (0 : ℝ) < 5 * Real.exp 1)]
    nlinarith

/-! ## 6 · 数值台账（`eΔ`、阈值）与整数 / 有理层备份证书 -/

/-- `eΔ = 20e = 54.3656365…`，夹在 `(54.3656, 54.3657)`。 -/
theorem eDelta_bounds :
    (54.3656 : ℝ) < 20 * Real.exp 1 ∧ 20 * Real.exp 1 < 54.3657 := by
  have h1 : (2.7182818283 : ℝ) < Real.exp 1 := Real.exp_one_gt_d9
  have h2 : Real.exp 1 < 2.7182818286 := Real.exp_one_lt_d9
  constructor <;> nlinarith

/-- 阈值 = `7.0632853…`：`7.063 < T < 7.0633`（注意 `7.0633` 是**上取整**：
`T = 7.06328…` 严格小于它，所以 `8 > 7.0633 > T`）。 -/
theorem thresh_bounds :
    (7.063 : ℝ) < 2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20) ∧
      2 * (2 : ℝ) ^ 2 * 48 / (Real.exp 1 * 20) < 7.0633 := by
  have h1 : (2.7182818283 : ℝ) < Real.exp 1 := Real.exp_one_gt_d9
  have h2 : Real.exp 1 < 2.7182818286 := Real.exp_one_lt_d9
  have hp : (0 : ℝ) < Real.exp 1 * 20 := by positivity
  constructor
  · rw [lt_div_iff₀ hp]
    nlinarith
  · rw [div_lt_iff₀ hp]
    nlinarith

/-- **纯有理层证书**（不依赖任何 `Real.exp` 引理；只要 `exp 1 > 27182818283/10¹⁰`）：
`7063/1000 < 384/(20·exp 1)`。这是 `(b)` 型备份 —— 若 mathlib 不可用，把 `exp 1`
换成它的有理下界 `27182818283/10000000000`，整条链仍由 `norm_num` 内核校验。 -/
theorem rational_certificate :
    (7063 : ℚ) / 1000 < 384 / (20 * (27182818283 / 10000000000)) := by
  norm_num

/-- 同一证书的整数交叉相乘形态：`7063·20·27182818283 < 384·10¹⁴`。 -/
theorem integer_certificate : 7063 * 20 * 27182818283 < 384 * 100000000000000 := by
  norm_num

/-- 锐性：`e = 12/5` 时阈值**恰等于** `8` ⇒ `8 > T ⟺ e > 12/5`，余量 = `e/(12/5)`。 -/
theorem sharp_at_e_2_4 : 2 * (2 : ℝ) ^ 2 * 48 / ((12 / 5 : ℝ) * 20) = 8 := by
  norm_num

end M02

/-! ## 7 · 公理审计（无 `sorryAx` 即「无 `sorry`」的机器证据） -/

#print axioms M02.cbar_eq_N
#print axioms M02.thresh
#print axioms M02.thresh_iff
#print axioms M02.ratio
#print axioms M02.x_gt_one
#print axioms M02.window
#print axioms M02.window_bounds
#print axioms M02.window_beta_bounds
#print axioms M02.thresh_bounds
#print axioms M02.eDelta_bounds
#print axioms M02.rational_certificate
#print axioms M02.integer_certificate

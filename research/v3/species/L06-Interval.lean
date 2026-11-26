/-
  L06-INTERVAL.lean
  =================
  问题：Lean/Mathlib 能否机器验证 q11 的 `Q11-CERT.md` 证书（定向舍入手算的
  S ∈ [..], D ∈ [..], S_lo − D_up > 0）？

  本文件只做**最小可编译例子**，回答两件事：
    (A) 实数超越不等式（`Real.exp`）能否在本机工具链上闭合？          → 能（下方 t1/t2）
    (B) `intervalIntegral` 能否用？                                   → 不能（import 直接缺 olean）
    (C) 证书的**代数骨架**能否形式化？                                → 能（cert_skeleton / q11_gap_shape）

  工具链（实测，非推断）：
    Lean 4.33.0-rc1 · Lake 5.0.0-src+62eed1d · Mathlib（vendor/OSforGFF 部分构建，2136 个 olean）
    运行方式：在 vendor/OSforGFF 里 `lake env lean <本文件>`。
  说明：默认 PATH 上的 `lean` 是 4.36.0，读不了 4.33.0-rc1 的 olean，
        必须走 `lake env` 才拿到匹配的 toolchain + LEAN_PATH。

  (B) 的实测报错（把下面这行取消注释即可复现）：
    -- import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
    error: object file '…\.lake\packages\mathlib\.lake\build\lib\lean\
            Mathlib\MeasureTheory\Integral\IntervalIntegral\Basic.olean'
            of module Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic does not exist
  本机 Mathlib 的 MeasureTheory（含 Integral/IntervalIntegral）**一个 olean 都没构建**，
  故 `∫ x in (0:ℝ)..1, x` 这类语句在本机无法编译。
-/

import Mathlib.Analysis.Complex.Exponential
-- ↑ 该 olean 本机已构建；`Real.exp`、`Real.exp_bound'`、`Real.exp_bound` 都在这里。
-- （`Mathlib.Analysis.SpecialFunctions.Exp` 也可，但 `exp_bound'` 定义在 Complex.Exponential。）

namespace L06

/-! ## 例子 A：纯有理/超越不等式 · 可编译、可闭合 -/

/-- 泰勒 + 显式余项：`exp x ≤ Σ_{m<n} x^m/m! + x^n (n+1)/(n!·n)`，取 `x = 5/96, n = 4`。
    数值：Σ = 1 + 5/96 + (5/96)²/2 + (5/96)³/6 ≈ 1.0534632，
          余项 ≈ (5/96)⁴·5/96 ≈ 3.83e-7，合计 ≈ 1.0534636 < 1.054。-/
theorem t1 : Real.exp (5 / 96 : ℝ) ≤ 1.054 := by
  have h := Real.exp_bound' (x := (5 / 96 : ℝ)) (by norm_num) (by norm_num) (n := 4) (by norm_num)
  have hnum :
      (∑ m ∈ Finset.range 4, (5 / 96 : ℝ) ^ m / (m.factorial : ℝ)) +
        (5 / 96 : ℝ) ^ 4 * ((4 : ℝ) + 1) / ((Nat.factorial 4 : ℝ) * 4) ≤ 1.054 := by
    norm_num [Finset.sum_range_succ, Nat.factorial]
  linarith

/-- 由 t1 直接放大：`4 * 1.054 = 4.216 < 7.0633`。 -/
theorem t2 : (4 : ℝ) * Real.exp (5 / 96) < 7.0633 := by
  have h := t1
  nlinarith [h]

/-! ## 例子 C：q11 证书的**代数骨架**（不含任何超越函数） -/

/-- 证书逻辑本身：区间端点比较 `D_up < S_lo` ⇒ `S − D > 0`。 -/
theorem cert_skeleton (S D S_lo D_up : ℝ) (hS : S_lo ≤ S) (hD : D ≤ D_up)
    (hgap : D_up < S_lo) : 0 < S - D := by
  linarith

/-- q11 的两个具体端点相减（`Q11-CERT.md` 的数）：`S_lo − D_up > 0`。
    `norm_num` 在有理数上直接闭合，不需要定向舍入库。 -/
theorem q11_gap_shape :
    (0.20511529712307314798 : ℝ) - (0.06912058608108493695 : ℝ) > 0 := by
  norm_num

end L06

-- 自检：#print axioms t1 应只依赖 [propext, Classical.choice, Quot.sound]，无 sorry。
#print axioms L06.t1
#print axioms L06.t2
#print axioms L06.q11_gap_shape

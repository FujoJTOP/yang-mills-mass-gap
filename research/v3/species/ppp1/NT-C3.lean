/-
# PPP1 / NT-C3.lean — 直路（`κ₃ ≠ 0 ⇒ ¬Gaussian` 用在 `q` 上）：为什么它是 `[NG]`

三件，按 §2 的判词组织：

1. **复合观测量不是 free-vanishing 见证。** `q = ¼tr F²`（=`1 − ½Re tr U_□` 的连续形）是
   `F` 的**二次**泛函（Wick 平方）；Gaussian `F` 的 `:F²:` 其**连通三点**本身非零
   （`research/nontriviality/REPORT.md:45` 逐字：`⟨:F²:(x):F²:(y):F²:(z)⟩_c ≠ 0`）。
   ⇒ 「`κ₃(q) ≠ 0 ⇒ 源非 Gaussian`」**是假推理**。本文件的有限精确形：
   `composite_third_cumulant_not_a_witness`（`ℚ` 上，`sorry` 为零）。
2. **直路真正需要的具名输入**：`C3Uniform`（`(a,L)`-一致非零下界）。
3. **`C3Uniform ⇒ NT`：已证**（`c3_uniform_implies_NT`）；**裸读数 ⇒ `C3Uniform`：`sorry`**
   （唯一一处；它标记的是**不存在的一步**，不是待证定理）。

## 依赖（逐行标注；本轮**未编译**）

* `import Mathlib` — 用 `abs`、`linarith`、`norm_num`（第 2、3 节的证明）。
* `import YangMillsInfra.Statement.NonTriviality` — `EuclideanTheory`/`TestFunction`/
  `IsGaussianTheory`/`nthCumulant`/`isGaussianTheory_cumulants_vanish`；
  lake 命令与 measure 侧对照见 `NT.lean` 头。
* 下面的 `NT` 与 `kappa3_implies_NT` 与 `NT.lean` 重复 5 行：`ppp1/` 不在任何 lake 包内，
  相对 `import` 不可解析（`lake` 只认包内模块名），故两个文件各自作**独立模块**。
  若日后并入 `YangMillsInfra`，删这两行、改 `import YangMillsInfra.Species.PPP1.NT`。
-/
import Mathlib
import YangMillsInfra.Statement.NonTriviality

noncomputable section

namespace PPP1

open YangMillsInfra

variable {d : ℕ}

/-- 见 `NT.lean`（本地副本，理由见文件头）。 -/
def NT (T : EuclideanTheory d) : Prop := ¬ IsGaussianTheory T

/-- 见 `NT.lean` 的 `nthCumulant_three_ne_zero_implies_NT`（本地副本）。 -/
theorem kappa3_implies_NT {T : EuclideanTheory d} {J : TestFunction d}
    (h : nthCumulant T J 3 ≠ 0) : NT T :=
  fun hG => h (isGaussianTheory_cumulants_vanish hG 3 (by decide) J)

/-! ## 1. 复合观测量：三阶累积量不是 free-vanishing 见证

`ℚ` 上的精确反例（无浮点、无采样）。源取偶律 `X ∈ {−2,−1,1,2}`，权重 `(1/8,3/8,3/8,1/8)`
（`m₁ = m₃ = 0` ⇒ `κ₃(源) = 0`，正是「Gaussian 侧应当通过」的那一半）；观测量取二次泛函
`Φ = X²`。`κ₃(Φ) = 81/32 ≠ 0`。这与 `REPORT.md:42–45` 的 `:F²:` 逐项同型。 -/

/-- 单变量三阶累积量：`κ₃ = m₃ − 3 m₂ m₁ + 2 m₁³`。 -/
def κ3 (m1 m2 m3 : ℚ) : ℚ := m3 - 3 * m2 * m1 + 2 * m1 ^ 3

/-- 偶源（`m₁ = m₃ = 0`）的三阶累积量恒为零。 -/
theorem even_source_kappa3 (m2 : ℚ) : κ3 0 m2 0 = 0 := by
  simp [κ3]

/-- 源的一阶矩：`m₁ = 0`。 -/
theorem source_m1 :
    (1 / 8 * (-2 : ℚ) + 3 / 8 * (-1 : ℚ) + 3 / 8 * (1 : ℚ) + 1 / 8 * (2 : ℚ)) = 0 := by
  norm_num

/-- 源的二阶矩：`m₂ = 7/4`。 -/
theorem source_m2 :
    (1 / 8 * (4 : ℚ) + 3 / 8 * (1 : ℚ) + 3 / 8 * (1 : ℚ) + 1 / 8 * (4 : ℚ)) = 7 / 4 := by
  norm_num

/-- 源的三阶矩：`m₃ = 0`。 -/
theorem source_m3 :
    (1 / 8 * (-8 : ℚ) + 3 / 8 * (-1 : ℚ) + 3 / 8 * (1 : ℚ) + 1 / 8 * (8 : ℚ)) = 0 := by
  norm_num

/-- 源的三阶累积量为零（`m₁ = m₃ = 0`）。 -/
theorem source_kappa3_zero : κ3 0 (7 / 4) 0 = 0 :=
  even_source_kappa3 (7 / 4)

/-- 二次观测量 `Φ = X²` 的矩：`m₁' = 7/4`、`m₂' = 19/4`、`m₃' = 67/4`。 -/
theorem composite_moments :
    (1 / 8 * (4 : ℚ) + 3 / 8 * (1 : ℚ) + 3 / 8 * (1 : ℚ) + 1 / 8 * (4 : ℚ)) = 7 / 4
  ∧ (1 / 8 * (16 : ℚ) + 3 / 8 * (1 : ℚ) + 3 / 8 * (1 : ℚ) + 1 / 8 * (16 : ℚ)) = 19 / 4
  ∧ (1 / 8 * (64 : ℚ) + 3 / 8 * (1 : ℚ) + 3 / 8 * (1 : ℚ) + 1 / 8 * (64 : ℚ)) = 67 / 4 := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

/-- 复合观测量 `Φ = X²` 的三阶累积量：`κ₃ = 81/32`。 -/
theorem composite_kappa3 : κ3 (7 / 4) (19 / 4) (67 / 4) = 81 / 32 := by
  norm_num [κ3]

/-- **空洞性（有限精确形，本文件主结论一）**：源偶（`κ₃ = 0`），复合观测量 `Φ = X²` 的
`κ₃ = 81/32 ≠ 0`。⇒ 「`κ₃(Φ) ≠ 0 ⇒ 源非 Gaussian`」**假**。

对应 `q = ¼tr F²`：`κ₃(q̄) ≠ 0` 是**自由/Wick 平方**也会给出的读数，故不能作 `NT` 的见证
（`REPORT.md:42–45,55–57`；`Basic.lean:39–43` 对连通 4 点写了同一句警告）。 -/
theorem composite_third_cumulant_not_a_witness :
    κ3 0 (7 / 4) 0 = 0 ∧ κ3 (7 / 4) (19 / 4) (67 / 4) ≠ 0 := by
  refine ⟨source_kappa3_zero, ?_⟩
  rw [composite_kappa3]
  norm_num

/-! ## 2. 直路真正需要的具名输入

`AA6-NT:45` 的幂次计数：`q = c·a⁴trF² + O(a⁶)` ⇒ `κ₃(q̄)` 的首项带 `a¹²`；
`AA6-NT:46` 说要用它必须先有「重整化三阶累积的 `(a,L)`-一致下界」。
把这句话写成类型： -/

/-- **具名输入（未证）**：重整化三阶累积量存在 `(a,L)`-一致的非零下界。

语义：`c ≤ |κ₃(⟨φ,J⟩)|` 对所有 `J` 与所有 `(a,L)` 成立（`c > 0` 与 `(a,L)` 无关）。
按 `AA6-NT:46`，其强度已接近 `NT` 本身。 -/
structure C3Uniform (T : EuclideanTheory d) : Prop where
  /-- 一致下界。 -/
  c : ℝ
  /-- 下界严格为正。 -/
  c_pos : 0 < c
  /-- 对每个线性测试泛函都成立。 -/
  witness : ∃ J : TestFunction d, c ≤ |nthCumulant T J 3|

/-- **`C3Uniform ⇒ NT`（已证，`sorry` 为零）**：一致下界 ⇒ 每个 `κ₃ ≠ 0` ⇒ 引理 C 逆否。 -/
theorem c3_uniform_implies_NT {T : EuclideanTheory d} (h : C3Uniform T) : NT T := by
  obtain ⟨c, hc, J, hJ⟩ := h
  refine kappa3_implies_NT (J := J) ?_
  intro hzero
  rw [hzero, abs_zero] at hJ
  linarith

/-! ## 3. 唯一 `sorry`：裸读数 ⇒ 一致下界

裸读数（`SU(2)`、`a=1`、`β=2.5`；`QQ2-ALOCAL:51–53`，**数据不是公理**）：

| `L` | 几何 | `κ₃(q̄)` | `±` |
|---|---|---|---|
| 6 | iso | `+8.620e−9` | `3.1e−9` |
| 4 | iso | `+1.037e−7` | `2.4e−8` |
| 6 | `6³×12` | `+1.772e−9` | `8.4e−10` |
| 6 | `6³×20` | `+5.878e−10` | `1.1e−10` |
| 8 | `8³×12` | `+1.076e−9` | `2.3e−10` |
| 0 | 8 副本合并 | `+6.272e−10` | `6.8e−11` |

六格**全正**但级差 `≈176×`（`1.037e−7 / 5.878e−10`），无标度律、无 `a→0` 族 ⇒ 给不出下界。 -/

/-- **直路的缺口（本文件唯一 `sorry`）**：把固定 `(L,β)` 的裸读数升成 `(a,L)`-一致下界。

**本定理不被主张为真**：它记录的是直路必须补上的**形状**，而不是一个待证命题；
§2 已判其前件不成立（`q` 是 Wick 平方，见 `composite_third_cumulant_not_a_witness`）。
保留 `sorry` 的目的是让缺口在 `#print axioms` 下显形（会打印 `sorryAx`），不作公理使用。

核验命令（本轮**未执行**）：`#print axioms PPP1.hole_bare_to_uniform`。 -/
theorem hole_bare_to_uniform (T : EuclideanTheory d) : C3Uniform T := by
  sorry

end PPP1

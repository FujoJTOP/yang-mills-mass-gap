# L06-INTERVAL · Lean 严格区间算术：`q11` 的证书能否机器验证？

**判：`(b)`** —— 例子 A **过**（`[T]`，exit 0，无 `sorryAx`）；例子 B **不过**（`[T]`，`intervalIntegral` 的 olean 不存在）。
`q11` 完整证书在本机**不可**机器验证（`[NG]`）——缺的不是技巧，是「**未构建的 Mathlib + 未定义的对象**」。
本轮**只**交：`research/v3/species/L06-Interval.lean` + 本文件。无 git commit。

---

## 0. 环境（实测，非推断）

| 项 | 值 |
|---|---|
| toolchain | `leanprover/lean4:v4.33.0-rc1`（elan 已装；另有 v4.36.0 作 default） |
| Lean | `4.33.0-rc1, x86_64-w64-windows-gnu, commit 62eed1d` |
| Lake | `5.0.0-src+62eed1d` |
| Mathlib | `vendor/OSforGFF/.lake/packages/mathlib`（`lakefile.toml` version `1.2.0`） |
| 已构建 olean | **2136** / 8268 个 `.lean`（25.8%） |
| `Mathlib/MeasureTheory` | 312 个 `.lean`，**0** 个 olean |
| `Mathlib/Analysis/Calculus` | 128 个 `.lean`，**0** 个 olean |
| `Mathlib/Analysis`（总） | 75 / 792 |
| PATH 上的 `lean` | **4.36.0**，读不了 4.33 的 olean ⇒ 必须走 `lake env` |

**环境抖动（实测，须记）**：同一文件同一命令，有一次报
`error: failed to read file '…\.elan\toolchains\leanprover--lean4---v4.33.0-rc1\lib\lean\Lean\Meta\Tactic\Grind\Injective.olean.private'`（exit 1）。
复查：该文件存在（432192 B），toolchain 内 `.olean.private` 共 2483 个；**立即重跑即 exit 0**、输出与第 1 节逐字一致。
⇒ 判定为**瞬时读失败**（与后台 `lake build` 的文件争用/扫描相关），非稳态缺陷；但说明本机编译偶发不可复现失败，任何"跑过了"都必须贴当次输出。

唯一可用运行方式（本文件全部输出均由此产生）：

```
cd D:\Dev\math\vendor\OSforGFF
lake env lean D:\Dev\math\research\v3\species\L06-Interval.lean
```

---

## 1. 例子 A：纯有理/超越不等式 —— **过** `[T]`

被测语句（`L06-Interval.lean`）：

```lean
theorem t1 : Real.exp (5 / 96 : ℝ) ≤ 1.054
theorem t2 : (4 : ℝ) * Real.exp (5 / 96) < 7.0633
```

**实际输出（逐字，第 3 次运行）**：

```
'L06.t1' depends on axioms: [propext, Classical.choice, Quot.sound]
'L06.t2' depends on axioms: [propext, Classical.choice, Quot.sound]
'L06.q11_gap_shape' depends on axioms: [propext, Classical.choice, Quot.sound]
[exit=0] elapsed: 00:00:18.3395234
```

无 `sorryAx`、无 error、退出码 0。

证明核（可直接复制）：`Real.exp_bound'` 给「Taylor 截断 + 显式余项」，余项是**有理式**，`norm_num` 收尾：

```lean
have h := Real.exp_bound' (x := (5 / 96 : ℝ)) (by norm_num) (by norm_num) (n := 4) (by norm_num)
have hnum :
    (∑ m ∈ Finset.range 4, (5 / 96 : ℝ) ^ m / (m.factorial : ℝ)) +
      (5 / 96 : ℝ) ^ 4 * ((4 : ℝ) + 1) / ((Nat.factorial 4 : ℝ) * 4) ≤ 1.054 := by
  norm_num [Finset.sum_range_succ, Nat.factorial]
linarith
```

数值核对：`Σ_{m<4} = 1 + 5/96 + (5/96)²/2 + (5/96)³/6 ≈ 1.0534632`，
余项 `= (5/96)⁴ · 5/96 ≈ 3.83e-7`，合计 `≈ 1.0534636 ≤ 1.054`（余量 `4.6e-4`）。
`t2` 由 `t1` 线性放大：`4×1.054 = 4.216 < 7.0633`。

### 1.1 一轮修正史（留档，防重复踩）

| 症状 | 根因 | 处方 |
|---|---|---|
| `Real.exp_le_one_add_self` unknown | 该名不存在 | 用 `Real.add_one_le_exp`（下界）或 `Real.exp_bound'`（上界） |
| `Real.exp_bound` 只给 `|exp x − Σ| ≤ R` | 是绝对值版 | 换 `Real.exp_bound'`，直接给 `exp x ≤ Σ + R` |
| ``Invalid field `factorial` … Real.factorial`` | `(4:ℝ).factorial` 非法 | 写 `Nat.factorial 4`（再 cast 到 ℝ） |
| `rw [show 4 = 1+1+1+1]` → `motive is not type correct` / `NeZero` 类型错 | 该重写也打中 `Nat.factorial 4`，dependent motive | 不必展开：单行 `norm_num [Finset.sum_range_succ, Nat.factorial]` |
| 试 1/2/3 三种写法 | — | 三种均通过（probe 文件无输出 = exit 0），取最短者 |

---

## 2. 例子 B：区间积分 —— **不过** `[T]`

被测语句：

```lean
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
example : (∫ x in (0:ℝ)..1, x) = 1/2 := by norm_num [integral_id]
```

**实际输出（逐字）**：

```
C:\Users\hooya\AppData\Local\Temp\l06probe\probeB.lean:1:0: error: object file 'D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build\lib\lean\Mathlib\MeasureTheory\Integral\IntervalIntegral\Basic.olean' of module Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic does not exist
```

补测 `lake build @mathlib/Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic`（24 s 抽样，oleans `2136 → 2144`）三条硬事实：

1. 目标在依赖序列的 `2143+ / 2169` 之后才可能被建 —— 即**上千个模块**在前；
2. 多个模块直接崩：`libc++abi: terminating due to uncaught exception of type std::bad_alloc: std::bad_alloc` → `error: Lean exited with code 3221226505`；
3. 另有 `error: Mathlib/Combinatorics/Enumerative/Composition.lean:6:0: failed to read file '…\Mathlib\Util\Notation3.ir'`。

⇒ 本机 Mathlib 既**不完整**（MeasureTheory 零 olean、Calculus 零 olean），又**不健康**（构建内存崩溃）。
`∫ x in (0)..1, x = 1/2` 在本机**无法编译**；30 分钟窗口内不可能补齐（缺口 ≈ 440 个 MeasureTheory+Calculus 文件及其传递依赖，且构建会 OOM）。

> 副作用记录：该探测在 `vendor/OSforGFF` 的 `.lake` 里新增 8 个 olean（2136→2144），**无源文件改动**。

---

## 3. 定向舍入：Lean 里可否？

| 面 | 状态 | 证据 |
|---|---|---|
| 符号 | **在** | `Mathlib.Analysis.SpecialFunctions.Log.Basic.olean`、`…Trigonometric.Basic.olean`、`…Exp.olean`、`Analysis/Complex/Exponential.olean` 均已建 |
| 定向舍入/区间算术库 | **不在** | 已构建 `Mathlib/Tactic` 36 个目录（Linarith, NormNum, Positivity, GCongr, Monotonicity, Bound, Ring, FieldSimp, …）**无** interval/interval-arithmetic；`Select-String 'syntax.*"interval"'` → **0 命中**（工具 OK） |
| `norm_num` 直取超越函数数值界 | **不支持** | `Mathlib/Tactic/NormNum/**` 里 `Real\.(exp\|log\|sin)` → **0 命中**（工具 OK） |
| 导数 / 中值定理（log、sin 余项证明的前提） | **不可用** | `Mathlib/Analysis/Calculus`：128 个 `.lean`，**0** 个 olean |

**替代（已验证可行，例子 A 即样板）：有理多项式上下界 + `norm_num`。**
模式：`超越点估值 = Taylor 截断（有理系数）+ 显式 Lagrange 型余项（有理式）`，余项用 `norm_num [Finset.sum_range_succ, Nat.factorial]` 闭合，跨项用 `linarith` 拼接。
- `exp`：现成 —— `Real.exp_bound'` / `Real.exp_bound`。
- `log`：需自写 `log(1+u) ≤ u − u²/2 + u³/3 − …` 的交替级数余项。
- `sin`：需自写 `sin u ≤ u − u³/6 + u⁵/120 − …`。
- **但** log/sin 的余项要用导数/中值定理，而 `Analysis/Calculus` 在本机 0 olean ⇒ 这是当前最硬的洞，不是排版问题。

---

## 4. `q11` 能否机器验证？

**当前工具链：不能**（`[NG]`）。**路径存在**（`[C]`），按依赖顺序，逐条引理：

| 序 | 引理 | 内容 | 现状 |
|---|---|---|---|
| **L0** | `Phi_def` / `Phi_pp_form` | 在 Lean 里**定义** Φ、Φ″，以及 `NS/NA/NB/Z` 四个 1D 整函数积分；证 `Φ″ = NS/Z`、`Var = β²(NB/Z − (NA/Z)²)` | **最重且不可绕**。q11 里它们只是 md 的解析式，Lean 里**无此对象** |
| **L1** | `exp_bracket` | `0 ≤ x → x ≤ 1 → exp x ≤ Σ_{m<n} x^m/m! + x^n(n+1)/(n!·n)` | **现成** `Real.exp_bound'`；例子 A 已用 `[T]` |
| **L2** | `log_bracket` | `log(1+u)` 交替级数有理上下界 | 自写（需 Calculus） |
| **L3** | `sin_bracket` | `sin u` 交替级数有理上下界 | 自写（需 Calculus） |
| **L4** | `bracket_add/mul/div` | 区间四则运算组合子（`linarith`/`nlinarith` 收尾，各 ~10 行） | 自写，廉价 |
| **L5** | `phi_pp_bracket` | 把 L1–L4 套到 L0 的表达式，得 `S ∈ [S_lo,S_up]`、`D ∈ [D_lo,D_up]`（**有理端点**） | 依赖 L0–L4 |
| **L6** | `integral_bracket` | `∫_a^b f ≤ (b−a)·M` 等区间积分单调性引理 | **卡在 `intervalIntegral`（第 2 节，本机缺）** |
| **L7** | `cert_skeleton` + `q11_gap_shape` | `S_lo ≤ S → D ≤ D_up → D_up < S_lo → 0 < S − D`；以及 `0.20511529712307314798 − 0.06912058608108493695 > 0` | **已在本文件证毕** `[T]`（`linarith` / `norm_num`） |

**即：证书的代数骨架本轮就能验证，缺的全是解析输入（L0–L6）。**

诚实口径：
- `S_lo − D_up = +0.1359947110419882` 这类端点在 Lean 里**没有复算**；形式化时它只能作为 L5 的**输出**，或作为 hypothesis（⇒ 得的是**条件定理**）。
- 现在的 trusted base 是 Python + mpmath 手算 + 定向舍入，**不是** Lean 内核。
- `q11` 若要在 Lean 里真验，第一步不是写不等式，是 **L0**：把 `Φ` 与四个整函数积分搬进 Lean。那一步的体量比"把不等式形式化"大一个量级。

---

## 5. 检索执行核验（AGENTS §1）

| 检索式 | 工具 | 状态 | 命中数 |
|---|---|---|---|
| `Test-Path …\IntervalIntegral\Basic.olean` | pwsh | OK | False（即"不存在"） |
| `Select-String 'syntax.*"interval"'` over `Mathlib/Tactic/**/*.lean` | pwsh | OK | 0 |
| `Select-String 'Real\.(exp\|log\|sin)'` over `Mathlib/Tactic/NormNum/**` | pwsh | OK | 0 |
| `Get-ChildItem …\Mathlib\MeasureTheory -Filter *.olean -Recurse` | pwsh | OK | 0 |
| `Get-ChildItem …\Mathlib\Analysis\Calculus -Filter *.olean -Recurse` | pwsh | OK | 0 |
| `lake env lean probeB.lean`（import IntervalIntegral） | lake env lean | OK（错误文本有效） | 1 error |

无 `TOOL-FAIL`。以上"0 命中"是**检索执行成功**后的 0，不是工具失败。

---

## 6. 纪律与边界

- 只写两个文件：本 md + `research/v3/species/L06-Interval.lean`。**无 git commit / push**。
- 非交付副作用：`lake build` 探测在 `vendor/OSforGFF/.lake` 新增 8 个 olean（2136→2144）；`$env:TEMP\l06probe\` 下有 3 个探测文件。
- `docs/position-current.html` **未改**（超出本路写权限）。若 lead 要落 §6 位置图，本轮结论变化为三条：
  1. `[T]` exp 型「Taylor + 有理余项 + `norm_num`」夹逼在本机**可编译可闭合**；
  2. `[T]` `intervalIntegral` 在本机**不可用**（olean 缺失 + 构建 OOM）；
  3. `[NG]` `q11` 完整证书在当前工具链下**不可**机器验证；缺 L0（对象未定义）与 L6（积分库缺失）。
- 口径：本文件所有数值均给来源文件（`Q11-CERT.md` / `L06-Interval.lean`），标签 `[T]/[C]/[NG]` 按实测分级；未证之物未标 `[T]`。

# L04 · 命题 C（`I2 ⟺ I1 ∧ ∃次小谱点`）与引理 A1（`s`-加权平均 `≥ m₁`）的 Lean 形式化

> 上游（已读正文）：`F01-TRANSPORTI2.md:18,29,39–60` · `D06-SPEC.md:11,17–35`
> 目标文件：`research/v3/lean/L04PropC.lean`（新，193 行，**0 个 `sorry`**）
> 标签：`[T]` 机器验证 / `[C]` 条件定理 / `[NG]` no-go / `[M]` 元陈述

---

## §0 工具状态（§1 纪律：先报工具，再报结论）

| 项 | 实测 | 判定 |
|---|---|---|
| `lean` | `Lean (version 4.36.0, x86_64-w64-windows-gnu, commit 611c117bfaf7b7ae1cfcf246ce6fba31e8f509e9, Release)`（`C:\Users\hooya\.elan\bin\lean.exe`） | **可用** |
| `lake` | `C:\Users\hooya\.elan\bin\lake.exe` | **可用** |
| mathlib 预编译 | `D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build\lib\lean` 下 **2136/2357** 模块（`.tools\mathlib-retry1.log` 停在 `[2351/2357]`）；`Mathlib.olean`（根）**不存在**；另有一份 8268 个 olean 在 `build.badcache`（已被判废） | **部分可用** |
| `Mathlib.Tactic` | `Mathlib\Tactic.olean` 不存在 ⇒ `linarith / ring / norm_num / positivity` **全部不可用** | 缺 |
| `MeasureTheory.*` | `MeasureTheory\Integral\Lebesgue.olean`、`MeasureTheory\MeasureSpace.olean` 均不存在 | **缺** ⇒ A1 的积分版**无法在本机陈述** |
| 实际编译 | `lake env lean D:\Dev\math\research\v3\lean\L04PropC.lean`（cwd = `D:\Dev\math\vendor\OSforGFF`） | **0 error**，仅 2 条 `push_neg` 弃用警告（见 `exists_min_above:124,131`） |
| 公理检查 | `#print axioms L04.a1 / L04.propC / L04.propC_literal / L04.counterexample_I2_not_I1 / L04.exists_min_above` → 全部 `[propext, Classical.choice, Quot.sound]`（**无 `sorryAx`**） | **无 sorry 依赖** |
| `sorry` 计数 | `Select-String -Pattern sorry -SimpleMatch` → **0** | — |

检索式（`Select-String`，工具 OK）：`F01-TRANSPORTI2.md` 上 `I2|I1|次小|Λ` → 25 命中；`D06-SPEC.md` 上 `A1|加权平均|上界|m₁|m1` → 25 命中。均有效响应，非 `TOOL-FAIL`。

---

## §1 命题 C

### §1.1 数学陈述（精确，先于 Lean）

固定一个**有限支撑、纯原子**的正测度（谱）：非空有限点集 `pts ⊂ ℝ`、权重 `wt : ℝ → ℝ` 非负。
取 `m₁ ∈ pts` 为最小点（`IsMin`）。

- `f₁² := wt m₁`（基态残数），**I1**：`f₁² > 0`。
- **次小谱点**：`SecExists := ∃ x ∈ pts, m₁ < x`（去掉最小点后支撑非空）。
- **Λ★**：`supp ∖ {m₁}` 的**最小点**（用谓词表达，不引选择函数）：
  `LtStar Λ := m₁ < Λ ∧ Λ ∈ pts ∧ ∀ x ∈ pts, m₁ < x → Λ ≤ x`。
- **I2**：`∃ Λ, LtStar Λ ∧ m₁² < Λ`。

**定理（命题 C 核）**：`0 < m₁ < 1` ⇒ `I2 ⟺ SecExists`。

### §1.2 Lean（逐字，来自 `L04PropC.lean`）

```lean
structure Spec where
  pts : List ℝ
  hne : pts ≠ []
  wt  : ℝ → ℝ
  hnn : ∀ x, 0 ≤ wt x

def I1 (S : Spec) (m1 : ℝ) : Prop := 0 < S.wt m1
def SecExists (S : Spec) (m1 : ℝ) : Prop := ∃ x, x ∈ S.pts ∧ m1 < x
def LtStar (S : Spec) (m1 Λ : ℝ) : Prop :=
  m1 < Λ ∧ Λ ∈ S.pts ∧ ∀ x ∈ S.pts, m1 < x → Λ ≤ x
def I2 (S : Spec) (m1 : ℝ) : Prop := ∃ Λ, LtStar S m1 Λ ∧ m1 ^ 2 < Λ

theorem propC (S : Spec) (m1 : ℝ) (h0 : 0 < m1) (h1 : m1 < 1) :
    I2 S m1 ↔ SecExists S m1
```

辅助（均已编译）：`exists_min`（非空有限实数集有最小元）、`exists_min_above`（去最小点后仍有最小元）、
`sq_lt_self : 0 < m1 → m1 < 1 → m1 ^ 2 < m1`（手证，因 `nlinarith` 缺）。
`propC` = `propC_fwd`（`I2 ⇒ SecExists`，无条件）+ `propC_bwd`（`SecExists ⇒ I2`，用窗口）。

### §1.3 裁决：`[T]`（**带两个适用范围**，不许省略）

1. **窗口 `0 < m₁ < 1` 不可去。** `Λ★ > m₁` 只在该窗口内蕴涵 `Λ★ > m₁²`（`sq_lt_self`）。
   窗口外有反例（数论算术，非 Lean）：`m₁ = 2`、次小点 `Λ = 3` ⇒ `Λ = 3 < m₁² = 4`，`SecExists` 真而 `I2` 假。
   本项目里 `m₁` 是质量（格点单位），`m₁ > 1` 的窗口并非先验排除 ⇒ **命题 C 的 `⟺` 是窗口条件定理**。
2. **有限支撑 = 纯原子谱。** `Λ★` 在我的模型里**取到**（离散最小点）。若剩下的是**连续谱**，则 `inf` 不被取到、`Λ★ = m₁²`，
   `I2` 假而 `SecExists` 真 —— 这正是 `F01:59` 引的 `R2 Lemma B` 反例（`I1 ⟹ I2` 假）。该情形**在我的有限支撑模型之外**，
   所以命题 C 的字面 `⟺` 只能读成「次小谱点 = 与底**分离**的第二个点」（否则它就是假的）。

**`I1` 一支是多余的**：若 `pts` 就是正权重支撑，`I1` 自动成立：

```lean
theorem I1_automatic (S : Spec) (hpos : ∀ x ∈ S.pts, 0 < S.wt x) (m1 : ℝ) (hmem : m1 ∈ S.pts) : I1 S m1
theorem propC_literal (S : Spec) (hpos : …) (m1 : ℝ) (hmem : m1 ∈ S.pts) (h0 : 0 < m1) (h1 : m1 < 1) :
    I2 S m1 ↔ I1 S m1 ∧ SecExists S m1      -- 机器验证
```

若 `pts` 只是「网格点表」（允许零权重点），则 **`I2` 不蕴涵 `I1`**，机器验证的反例：

```lean
noncomputable def demoWt : ℝ → ℝ := fun x => if x = 0 then 0 else 1
noncomputable def demoSpec : Spec := { pts := [0, 1], … , wt := demoWt, … }
theorem counterexample_I2_not_I1 : I2 demoSpec 0 ∧ ¬ I1 demoSpec 0   -- 机器验证，0 sorry
```
⇒ 字面 `I2 ⟺ I1 ∧ ∃次小谱点` 只在 **supp = 正权重支撑** 的读法下成立；该读法下 `I1` 是虚位合取项。

**标签**：`propC` / `propC_literal` / `counterexample_I2_not_I1` 均 **`[T]`（机器验证）**；
「命题 C 逐字成立」的普适版 **`[C]`**（窗口 `0<m₁<1` + 有限原子支撑/分离次小点）。

---

## §2 引理 A1

### §2.1 A1 的原文（`D06-SPEC.md:19–31`，积分形式）

设 `dν ≥ 0` 于 `[0,∞)`，`0 < ∫dν < ∞`，`m₁ := inf supp ν`；`w ≥ 0` 可测，`0 < ∫ w dν < ∞`，

```
A_w(ν) := (∫ s w(s) dν(s)) / (∫ w(s) dν(s))        -- ν 上的 s-加权平均
A_w(ν) ≥ m₁(ν)                                      等号 ⟺ w·ν 是单点质量
```
证明两步：`(i)` 在 `supp ν` 上逐点 `s ≥ m₁`；`(ii)` `w ≥ 0` ⇒ 加权平均保序。∎

### §2.2 判：A1 与凸性引理**是否同一**？—— **不同一，是它的「积分特化 + 元层延伸」**

| 层 | 对象 | 是否数学命题 | Lean |
|---|---|---|---|
| ① **凸性引理本身**（父任务候选 `lemma a1`，有限和） | `mmin ≤ ∑ i, w i * m i`，`w ≥ 0`，`∑ w = 1` | **是** | **`[T]` 机器验证**（下方 `L04.a1`，0 sorry） |
| ② **A1 原文**（`D06:19–31`，`Measure` + 积分 + `inf supp`） | `∫ s w dν / ∫ w dν ≥ inf supp ν` | **是**（同一思想的一步升级） | **`[C]`（未机器验证）**：本机无 `MeasureTheory` olean，连陈述都编译不过 |
| ③ **「20 条死路算的都是地板的上界」** | 对 20 个具体估计量的**分类断言** | **不是**（元陈述 / `[M]`） | 不可形式化为单条定理：要对每条死路逐个形式化它算的量，才能说它是某个 `A_w(ν)` |

所以：**A1 的数学核 = 凸性引理（取积分而非有限和）；A1 的「上界」结论 = 元层断言**（`[M]`，需逐条死路的形式化清单，不是一行 Lean）。
①→② 的差别仅是把 `∑_{i} w_i` 换成 `∫ w dν`（`ν`-a.e. 逐点 `s ≥ m₁` + 保序性 ⇒ 两边同乘 `w` 再除正数）；证明骨架同构，但**积分版没有机器验证**，不许标 `[T]`。

### §2.3 A1 的 Lean（有限和形式，机器验证）

```lean
theorem a1 {ι : Type*} [Fintype ι] (w m : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (h1 : ∑ i, w i = 1) (mmin : ℝ) (h : ∀ i, mmin ≤ m i) :
    mmin ≤ ∑ i, w i * m i
```

证明（无 tactic 自动化，因 `Mathlib.Tactic` 缺）：先证 `mmin * (∑ w i) ≤ ∑ w i * m i`
（`mul_sum_eq`：`c * ∑ w = ∑ c * w`，自证 `Finset.induction_on` + `mul_add`；
逐项 `mul_le_mul_of_nonneg_left`），再代 `h1`、`mul_one`。
`#print axioms L04.a1` = `[propext, Classical.choice, Quot.sound]` ⇒ **无 `sorry`、无额外公理**。

### §2.4 A1 的积分版（`[C]`，手写证明，未编译）

```lean
-- 需要 Mathlib.MeasureTheory.Integral.*（本机无 olean），故只给陈述与证明
-- theorem a1_integral (ν : Measure ℝ) [IsProbabilityMeasure ν] (w : ℝ → ℝ)
--     (hw : ∀ x, 0 ≤ w x) (hwi : 0 < ∫ x, w x ∂ν)
--     (mmin : ℝ) (h : ∀ᵐ x ∂ν, mmin ≤ x) :
--     mmin ≤ (∫ x, x * w x ∂ν) / (∫ x, w x ∂ν)
```
手证：`h` 给出 `∫ x, w x ∂ν`-方向逐点 `mmin * w x ≤ x * w x`（`mul_le_mul_of_nonneg_left` + `mul_comm`）；
`integral_mono`/`integral_const_mul` 得 `mmin * ∫ w ≤ ∫ x, x * w x`；`0 < ∫ w` 除过去即得。∎
（注：原文本写 `m₁ := inf supp ν`，Lean 里 `inf supp` 需 `MeasureTheory.Measure.support` + `sInf`，比 `∀ᵐ x ∂ν, mmin ≤ x` 更弱形式的假设更贵；此处取后者为等价的工作形式。）

---

## §3 判死线

| 支 | 判定 | 依据 |
|---|---|---|
| **(a)** 两者都可编译且 `sorry ≤ 2` ⇒ 升 `[T]` | **半命中**：`sorry = 0`（比 `≤2` 更严），`lean` 编译 0 error。**命题 C 核 + A1 有限和版 均升为机器验证 `[T]`** | §0 编译与公理记录 |
| **(b)** A1 与凸性引理不同一 ⇒ 指出差别 | **命中**：A1 原文是 `Measure`/积分（②），不是有限和（①）；且「20 条死路都是上界」是**元陈述 `[M]`**，不是数学命题 | §2.2 表 |
| **(c)** Lean 不可用 ⇒ 只交陈述 + 手证 | **不适用**（Lean 可用），但适用其变体：`MeasureTheory` 缺 ⇒ **A1 积分版降为 `[C]`** | §0 工具表 |

**`sorry` 数 = 0，位置 = 无。** 唯一未机器验证的是 §2.4 的 A1 积分版（工具缺口，不是逻辑缺口）。

**结论一句话**：命题 C 核 `I2 ⟺ 次小谱点`（窗口 `0<m₁<1`、有限原子支撑）**`[T]`**；字面的 `I1 ∧` 在「supp = 正权重支撑」下是虚位合取项，否则 `I2 ⇏ I1` 有机器验证反例；
A1 **不是**凸性引理逐字同一 —— 它是积分版凸性引理（`[C]`，本机不可陈述）+ 一条 `[M]` 的元断言（20 条死路上界）。

---

## §4 复现

```powershell
cd D:\Dev\math\vendor\OSforGFF
lake env lean D:\Dev\math\research\v3\lean\L04PropC.lean      # 0 error
```

依赖的两个 mathlib 模块：`Mathlib.Data.Real.Basic`、`Mathlib.Algebra.Order.BigOperators.Group.Finset`（本机已预编译）。
不用 `Mathlib.Tactic`、不用 `MeasureTheory`（均缺 olean）。

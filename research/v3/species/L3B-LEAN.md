# L3B — Lean 装配：`I2 ∧ (★) ⇒ m > 0`（有限谱模型）

**状态**：`sorry = 0`，`exit = 0`，7 个声明全部 `#print axioms` 干净（仅 `propext / Classical.choice / Quot.sound`，**无 `sorryAx`**）。
**文件**：`research/v3/lean/L3B-Assembly.lean`（133 行）。
**纪律**：`cd vendor\OSforGFF; lake env lean <file>`，toolchain `leanprover/lean4:v4.33.0-rc1`（`lake env lean --version` → `Lean (version 4.33.0-rc1, ...)`）。olean 走 `mathlib\.lake\build\lib\lean`（`Mathlib.Data.Real.Basic / Finset.Basic / Analysis.SpecialFunctions.Pow.Real / Order.MinMax / Tactic.Linarith` 均命中；`build.badcache` 未用到，`build` 树已有本次所需 olean）。

---

## 1. 精确陈述（逐字，`L3B-Assembly.lean:25`）

```lean
theorem assembly (spec : Finset ℝ) (m₀ δ c : ℝ)
    (hm : 0 ≤ m₀) (hδ : 0 < δ) (hc : 0 < c)
    (h1 : ∀ x ∈ spec, x = m₀ ∨ x ≥ m₀ + δ)          -- I2：0⁺⁺ 道隔离
    (hstar : ∀ y ∈ spec, y ≠ m₀ → y ≥ c) :          -- (★)：其他道地板
    ∃ g > 0, ∀ x ∈ spec, x ≠ m₀ → x ≥ g := by        -- 结论：全谱隙
  refine ⟨min δ c, lt_min hδ hc, ?_⟩
  intro x hx hne
  rcases h1 x hx with h | h
  · exact absurd h hne
  · exact le_trans (min_le_left δ c) (le_trans (by linarith) h)
```

与任务给的签名相比只多两条显式常数正性（`hδ, hc`，本来就在 `I2`/`(★)` 的定义里）和一条 `0 ≤ m₀`（见 §4，这是**唯一**用到 `m₀` 符号的地方）。`g := min δ c`，证明三行。

## 2. 编译结果

```
D:\Dev\math\research\v3\lean\L3B-Assembly.lean:28:5: warning: Variable name `hstar` is not explicitly referenced.
'MassGap.L3B.witness_satisfiable' depends on axioms: [propext, Classical.choice, Quot.sound]
'MassGap.L3B.assembly' depends on axioms: [propext, Classical.choice, Quot.sound]
'MassGap.L3B.assembly_of_I2' depends on axioms: [propext, Classical.choice, Quot.sound]
'MassGap.L3B.assembly_two_sector' depends on axioms: [propext, Classical.choice, Quot.sound]
'MassGap.L3B.assembly_set' depends on axioms: [propext, Classical.choice, Quot.sound]
'MassGap.L3B.assembly_no_ground_sign' depends on axioms: [propext, Classical.choice, Quot.sound]
'MassGap.L3B.assembly_mass' depends on axioms: [propext, Classical.choice, Quot.sound]
[exit=0]
```

`sorry` 出现次数：**0**（全文唯一匹配是第 2 行注释里的字 "sorry-free"）。冷启 `olean.private` 坑本轮未触发。

## 3. 判定

**(a) `sorry = 0` ⇒ L3 的装配已机器验证** —— 但要说清"验证了什么"：机器验证的是**装配这一步**（隔离 + 地板 ⇒ 全谱隙，常数 `min δ c`），**不是** `I2` 或 `(★)` 本身。`I2`、`(★)` 在本文件里仍是假设（`h1`、`hstar`），从测度/转移矩阵推出它们不在本文件范围。故：

- `[T]` `assembly`：`I2 ∧ (★) ⇒ ∃ g > 0, 其余谱 ≥ g`，无条件成立，`sorry`-free。
- `[T]` `assembly_mass`：由 `g` 得质量 `m := √g > 0`，`m² ≤ x`（`Real.sq_sqrt`）。
- `[C]` `spec` 取为"平方质量谱"的物理映射（本文件只用实数序结构，未建物理模型）。

**(b) 是否需要 `(★)` 的更强形式**：**不需要。** 现有 `∀ y ∈ spec, y ≠ m₀ → y ≥ c`（`c > 0`）已足够，且**不需要有限性**——`assembly_set`（`L3B-Assembly.lean:91`）在任意 `Set ℝ` 上逐字复现同一证明。两点精确的边界：

1. **本形式下 `(★)` 甚至冗余**：`assembly_of_I2`（:40）证明只用 `I2` 就给出 `g := δ`。机器证据就是编译首行的 linter 警告 `Variable name 'hstar' is not explicitly referenced` —— `(★)` 在主定理里没被引用。
2. **`(★)` 真正出力的是 `m₀` 符号未知的场合**：`min δ c` 要为正需要 `0 ≤ m₀`（把 `x ≥ m₀ + δ` 降成 `x ≥ δ`）。去掉 `hm` 后 `min δ c` 不保证正，此时靠 `(★)` 的 `c` 作绝对地板——`assembly_no_ground_sign`（:68）。所以"更强"的 `(★)` 无意义；需要补的是 `0 ≤ m₀`（平方质量非负，物理上免费）或直接以 `c` 为地板。
3. 两扇区拆开时 `(★)` 才按物理原意出手：`assembly_two_sector`（:52），`I2` 管 `A`、`(★)` 管 `B`，结论 `min δ c`。

**(c) Lean 环境**：通。`lake env lean` 返回 `exit=0`；唯一注意点是 import 面（`Finset`/`Real.sqrt` 不在 `Mathlib.Data.Real.Basic` 里，需 `Mathlib.Data.Finset.Basic` + `Mathlib.Analysis.SpecialFunctions.Pow.Real`）。首轮不加这两条 import 会得到 `Unknown identifier 'Finset'`——那是 import 面问题，不是装配问题。

## 4. 诚实备注

1. 装配是**簿记**，不是物理难点。真正的缺口仍是 `I2`、`(★)` 的无条件推导；本文件把"给定隔离 ⇒ 全谱隙"这一步钉死，防止后续在常数上打架。
2. `witness_satisfiable`（:105）给出 `spec = {0,1}, m₀ = 0, δ = 1, c = 1` 的具体见证，证明假设体系可满足、定理不 vacuous。
3. 局限：`Finset` 版本不涉及谱的重数、不涉及算子自伴性，也不涉及 `Λ★` 与 `m₁²` 的关系；`Λ★ − m₁² ≥ δ` 只以 `x ≥ m₀ + δ` 的序形式进入。

## 5. 复现

```
cd D:\Dev\math\vendor\OSforGFF
lake env lean D:\Dev\math\research\v3\lean\L3B-Assembly.lean   # exit=0
```

标签汇总：`[T]` assembly / assembly_of_I2 / assembly_two_sector / assembly_set / assembly_no_ground_sign / assembly_mass / witness_satisfiable；`[C]` spec = 平方质量谱的物理读法。

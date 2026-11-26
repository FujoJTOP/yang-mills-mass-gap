# L05 · 两个基准反例（M6 / 2D YM）的形式化 —— "靶没定对"的判据

> 轮次：第 N 轮（时间戳 2026-11-25 16:32–17:00 America/Toronto）
> 文件：`research/v3/lean/L05Counterexamples.lean`（唯一新增 `.lean`）
> 工具：Lean 4.36.0 + Mathlib（部分 olean）—— `lake env lean` 实际输出见 §3
> 标签：`[T]` 机器验证 / `[C]` 条件（需说明前提）/ `[NG]` 未证（缺什么已写清）

---

## 0. 靶（被攻击的对象）

```
T1:  G″ > 0  ⇒  m > 0        （4D YM 五环"钥匙"的单向蕴含形式）
```

两个基准反例各自给一对"一个正、一个零"的数值：

| 反例 | 正 | 零 | 机制 |
|---|---|---|---|
| **M6**（谐振子 ⊗ 自由方向） | `G″ = 2ω > 0` | `m = 0` | `G″ = 1/Var(a)`，`Var(a)=1/(2ω)`；b 方向**无势** ⇒ 连续谱 `[0,∞)` |
| **2D YM** | `m = 3g²/8 > 0` | `f₁² ≡ 0` | `tr F² = H` 守恒 |

**"最省"排序（本文件遵守）**：有限维 > 需 Gaussian > 需谱论。

---

## 1. 结论先说（判死线 (a)/(b)/(c)）

**(c) 成立、且是最省的载体，但它只杀裸蕴含。**
`2×2` 对称模型 `H = !![0,c; c,e]`（`H e₀ = 0`，`e₀ = (1,0)`）上，本文件机器验证了**二分法** `[T]`：

- `e = 0`（**真空简并**）⇒ `NoGap`：`∀ε>0` 存在归一化 `v ⊥ e₀` 使 `⟨v,Hv⟩ < ε`；
- `0 < e`（**简单真空**）⇒ **必有隙**：所有激发能量 `> 0`；
- 且非负性强制 `e = 0 ⇒ c = 0`，**没有第三条路**。

> 它算不算合法反例？**对裸蕴含 `G″>0 ⇒ m>0`：算**（`[T]`）。
> **对 M6 的物理机制：不算**（`[C]`→实际是 `[NG]`）。理由：M6 的 `m = 0` 来自 b 方向
> **连续谱**配**简单真空**；有限维里 `gap = 0` 只能来自**真空简并**，是被 fine-tune 出来的
> 退化，不是"自由方向"的 IR 机制。两者"`gap=0`"同名不同源。

**(b) Gaussian：本 checkout 不可用（不是缺定理，是缺 olean）。**
- 源码**有** `ProbabilityTheory.gaussianReal`（`Mathlib/Probability/Distributions/Gaussian/Real.lean:222`）`[T-源码]`；
- 源码**没有**名为 `variance_gaussian` 的引理 —— 检索式 `def gaussianReal|variance_gaussian|def gaussian \(|theorem gaussianReal`，工具 grep OK，命中 1（即 `gaussianReal` 本身），`variance_gaussian` 命中 **0**；
- 实际可用的是 `variance_id_gaussianReal : Var[id; gaussianReal μ v] = v`（`Real.lean:575`）、`variance_fun_id_gaussianReal`（`:550`）；
- **但 Probability 一棵树的 olean = 0 个**，`lake env lean` 报
  `object file .../Mathlib/Probability/Moments/Variance.olean does not exist`。
  ⇒ 要用 `gaussianReal` 必须整构建 Mathlib（超出本轮）。
- **因此 M6 的 `Var(a) = 1/(2ω)` 不能从 Gaussian 积分推出**；本轮把它当**模型数据**（假设）带入 ——
  这一步是诚实的，因为 `Var(a)=1/(2ω)` 是谐振子高斯测度的定义性质，算术核不需要它。

**(a) 升为机器验证的部分**：`G″ = 2ω > 0` 的算术核 + 有限维 `gap = 0` 见证 +
二分法，全部 `[T]`，`sorry` 数 = **0**，公理只有 `propext / Classical.choice / Quot.sound`。

---

## 2. "无隙"的最省表述（判：Rayleigh 形式，不用谱）

`σ(H) = [0,∞)` 需要无界自伴算子谱定理（Mathlib 无）⇒ **不取**。
取**变分/Rayleigh 无隙见证**：

```lean
def NoGap (c e : ℚ) : Prop := ∀ ε : ℚ, 0 < ε → ∃ v : V2, Excited v ∧ q c e v < ε
-- Excited v  := v.1 = 0 ∧ v.1^2 + v.2^2 = 1     -- 与真空 e₀=(1,0) 正交且归一
-- q c e v    := 2*c*v.1*v.2 + e*v.2^2           -- ⟨v, H v⟩, H = !![0,c; c,e]
```

理由：只要 `H ≥ 0`，"与真空正交的态里存在任意低能量" ⟺ `inf(σ H ∩ (0,∞)) = 0`，
而这条**完全不需要谱论**，只需二次型算术 —— 这正是"最省"的定义选择。

---

## 3. Lean 实际输出（工具状态：`lake env lean` OK，exit 0）

工作目录 `vendor/OSforGFF/.lake/packages/mathlib`（唯一带 Mathlib olean 的包），命令
`lake env lean D:\Dev\math\research\v3\lean\L05Counterexamples.lean`：

```
'L05.m6_Gsecond_pos' depends on axioms: [propext, Classical.choice, Quot.sound]
'L05.degenerate_NoGap' depends on axioms: [propext, Classical.choice, Quot.sound]
'L05.simple_vacuum_gap' depends on axioms: [propext, Classical.choice, Quot.sound]
'L05.dichotomy' depends on axioms: [propext, Classical.choice, Quot.sound]
'L05.finite_dim_witness' depends on axioms: [propext, Classical.choice, Quot.sound]
```

（无 `sorryAx` ⇒ `sorry` 数 0；无报错，仅一条 unused-variable linter 警告，已修。）

### 已验证条目

| 条目 | 内容 | 状态 |
|---|---|---|
| `gsecond_harmonic` | `G″(1/(2ω)) = 2ω`（字面形式，`ω≠0`） | `[T]` |
| `m6_Gsecond` | `Var(a)·2ω = 1 ∧ G″·Var(a) = 1 ⇒ G″ = 2ω` | `[T]` |
| `m6_Gsecond_pos` | `G″ > 0`（`ω>0`） | `[T]` |
| `m6_Gsecond_literal` | `G″ = 1/(1/(2ω))` | `[T]` |
| `degenerate_NoGap` | 简并真空 ⇒ `NoGap`（`v=(0,1)`，零能激发） | `[T]` |
| `simple_vacuum_gap` | 简单真空 ⇒ 所有激发能量 `> 0` | `[T]` |
| `dichotomy` | `2×2` 对称非负 + `H e₀=0` ⇒ `e=0` 时 `c=0`（无第三条路） | `[T]` |
| `finite_dim_witness` | 同一模型：`G″ > 0 ∧ NoGap` | `[T]` |

---

## 4. 明说缺什么（`[NG]`，逐条）

1. **Gaussian 积分**：`ProbabilityTheory.gaussianReal` 的 olean 未构建（Probability 全树 0 olean）；
   `variance_gaussian` 这个**名字在 Mathlib 源码里不存在**（0 命中），要的是 `variance_id_gaussianReal`。
   ⇒ `Var(a)=1/(2ω)` 无法从积分算出，本轮作假设。
2. **M6 真正的 `m = 0`**：需要自由 Laplacian 在 `L²(ℝ)` 上 `σ = [0,∞)`（或"连续谱"），
   Mathlib 无无界自伴算子谱定理 ⇒ **无法机器验证**。这是 M6 反例**唯一**的实质缺口。
3. **2D YM**：`[H,O]=0 ⇒ f₁² = 0` 的最省形式（算符对易 + 矩阵元 `⟨𝔾|O|𝔾⟩`）
   需先在 Lean 里把 2D 格点 `tr F²` 建成有限维自伴算子，本轮 25 分钟内**未做**。
   判：它比 M6 更"便宜"（有限维矩阵，无测度），但**比 (c) 贵**（要建格点模型），
   且它的"`f₁² ≡ 0`"是**恒等式**而非"`G″` 与隙的错配"，逻辑形态与 M6 不同。
4. **`G″ = 1/Var(a)` 本身**：本轮当定义用。若它也是待证命题，则 L05 的两条反例都只到条件层面。

---

## 5. 对"靶没定对"的裁决

- **`G″ > 0 ⇒ m > 0` 作为无条件蕴含：假**，`[T]`（有限维反例，`finite_dim_witness`）。
- **但**有限维反例的 `gap=0` 靠真空简并，M6 的 `gap=0` 靠连续谱。
  ⇒ 若要"物理同类"的反例，必须加**非简并真空**假设；加了之后
  `2×2` 非负对称类里 `G″>0 ⇒ gap>0`：`e=0` 分支 `[T]`（`dichotomy` ⇒ `c=0` ⇒ `q≡0` 才是无隙的那种），
  `e>0` 分支 `[T]`（`simple_vacuum_gap`）；`c≠0 ∧ e>0` 与 q ≥ 0 矛盾为 [C]（手工一行：q(1,-c/e) = -c²/e < 0，未机验）。
- 因此靶需要改成：`G″ > 0 ∧ 真空非简并 ∧ <XX> ⇒ m > 0`，缺的 `<XX>` 就是
  **"b 方向那种连续谱在 4D YM 里不出现"**（= 无质量自由度不存在）。这正是本轮要指出的结构缺口：
  `G″>0` 管不了 IR 有没有零模，它只管 a 方向涨落；`m>0` 管的是**有没有任意低激发**。
  两个量测的不是同一件事 —— **这就是"靶没定对"的机器可验证证据（有限维版）**。

---

## 6. 复现

```powershell
cd D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib
lake env lean D:\Dev\math\research\v3\lean\L05Counterexamples.lean
```

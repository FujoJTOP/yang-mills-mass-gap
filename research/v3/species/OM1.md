# OM1 — `Ω` 到底是什么？（`gap2` 引 `R1:4.1` 的"缺的是 `Ω` 不是数"）

> 本轮唯一产出文件。只读核 + 定性，没有任何新数值。
> 纪律：`R1:4.1` 逐字是第一问；不编；标签 `[T]/[C]/[NG]/[O]`。

---

## §0 检索执行记录（`AGENTS §1`）

| 项 | 值 |
|---|---|
| 工具状态 | `pwsh` OK；`Select-String`（等价 grep）OK；`read` OK。无 `TOOL-FAIL`，无 `ENOENT`，无 `ripgrep provider failure` |
| 范围 | `D:\Dev\math\research` 递归 `*.md`，文件数 **1422** |
| 检索式 1 | `Select-String -Pattern 'Ω'`（全 `research/**/*.md`） ⇒ **命中 1917 行 / 345 文件** |
| 检索式 2 | `Get-ChildItem research -Recurse -Include *.md \| ? Name -match '^R1'` ⇒ 命中 5：`research/unconditional/R1-small-S1-bions.md`、`research/v3/species/R1-BRIDGE.md`、`research/v3/R1-RESIDUE-TRANSPORT.md`、`R10/R11/R12-*`（非 R1） |
| 检索式 3 | `Select-String -Pattern 'R1[:：]' research\v3\species\LOOP-RUNBOOK.md` ⇒ **命中 5**（6666、6671、17839、17877、17897） |
| 检索式 4 | `Select-String -Path research\v3\R1-RESIDUE-TRANSPORT.md -Pattern '^#+.*4\.1\|Ω'` ⇒ **命中 §4.1 于 L181**；`(4.1)` 在 L188、`(4.2)` 在 L189 |
| 消歧 | `R1:4.1` 的两个候选：`v3/species/R1-BRIDGE.md §4.1` 是"判定"节（转移矩阵桥，讲 `\|Ω\|²` 真空测度），`v3/R1-RESIDUE-TRANSPORT.md §4.1` 是 "The honest identity"（讲搬运算子 `Z_O` + 态重叠）。`gap2` 的上下文是"**不搬 β** ⇒ 不搬 `Z_O` ⇒ 缺的因子"⇒ 只有 **`R1-RESIDUE-TRANSPORT.md`** 对得上。 |

---

## §1 第一问：`R1:4.1` 的逐字

`research/v3/R1-RESIDUE-TRANSPORT.md:181` 起（`§4 · The exact decomposition: the transport carries the *value*, the missing factor is the *overlap*`）：

> **### 4.1 The honest identity**
>
> Let `f₁²_{lat}(β)` be the lightest `0⁺⁺` atom of the *lattice* theory at coupling `β` (a residue of the
> finite-`L` transfer matrix, which is always `>0` by compactness). To compare `β_eff<1/48` and `β_tgt=2.5`
> we must identify the lattice theories as continuum regularizations and the states as the *same* physical
> glueball. Inserting the two identifications:
> ```
>    f₁²_{phys}(β_tgt) = Z_O²(μ_tgt,μ_eff) · Ω(β_tgt,β_eff) · f₁²_{lat}(β_eff) ,        (4.1)
>    Ω(β_tgt,β_eff) := |⟨G(β_tgt)|G(β_eff)⟩_O|²   (O-weighted overlap of the two lightest atoms).  (4.2)
> ```
> Here `Z_O²` is the exact §1 covariance (it acts on the operator), and `Ω` is the **state** factor — the
> squared `O`-overlap between the lightest `0⁺⁺` eigenstate at the two couplings. The domain of (4.1)
> requires both ends to be in the scaling window and the state to be the same. The literature identity
> `f₁²=Z_O² f₁²` used by the corpus is the special case `Ω=1`.

同文件 `:196–210`（`§4.2 What Ω is, and why it is exactly the missing step`）：

> `Ω` is `1` **iff** the strong-coupling "glueball" is the physical glueball, carried without mixing across
> the crossover. … `Ω<1` and can be `0`. …
> - `Ω=1` ⇔ the lightest atom is the *same* atom, with the same quantum numbers, isolated on `[β_eff,β_tgt]`
>   ⇒ **uniform isolation `I₂`** (`Λ_β−m₁²(β)≥δ>0`);
> - `Ω=0` ⇔ the `O`-matrix element of the lightest state vanishes (an **overlap zero**, `a(β_tgt)=0`), or
>   the atom has dissolved — `f₁²` touches `0`.
>
> So **`Ω` is isolation `I₂`, is the single-`T=0`-phase statement `(CP)`, and is the gap.**

同文件 `:293–296`（`§8`）：

> **Missing step (exact).** A non-perturbative lower bound
> ```
>    Ω(β_tgt,β_eff) = |⟨G(β_tgt)|O|G(β_eff)⟩|² ≥ ω₀ > 0 ,        (P-Ω)
> ```

`gap2` 的逐字在 `research/v3/species/LOOP-RUNBOOK.md:17877`：`**缺的是 `Ω`，不是数**（`R1:4.1`）`。**引文对得上，不是错引。** `[T]`

---

## §2 第二问：哪个 `Ω`？

**答案是 (4.2)：`Ω` 是一个标量因子 —— 两个耦合下最轻 `0⁺⁺` 本征态之间的 `O`-加权重叠模方，取值 `∈[0,1]`。** `[T]`（定义即此）

**不是**真空态向量，**不是**正交补，**不是**代数，**不是**构型空间。见 §3 撞名表。

对齐关系（`R1:191–194` `[T]`）：

- `Z_O` = **算子**侧的协变（作用在 `O` 上）；对 `f₁²` 是乘性的 `Z_O²`。
- `Ω` = **态**侧的因子（作用在态上）；语料一直隐含取 `Ω=1`。
- 所以 `f₁² = Z_O²·f₁²`（语料式） = `(4.1)` 在 `Ω=1` 的**特例**。

---

## §3 ★ 撞名表（`Ω` 在语料里至少 4 义，必须带限定词才能写）

| # | 义 | 逐字 / 出处 | 类型 |
|---|---|---|---|
| **Ω-1（本义）** | **态重叠因子** `Ω(β_tgt,β_eff)=\|⟨G(β_tgt)\|G(β_eff)⟩_O\|²∈[0,1]` | `v3/R1-RESIDUE-TRANSPORT.md:189` `(4.2)`；`LOOP-RUNBOOK.md:6340` `Ω = \|⟨G_{2.5} \| G_{1/48}⟩\|² ∈ [0,1]` | 标量；`[T]` 定义 / `[O]` 下界 |
| **Ω-2** | **真空态向量** `\|Ω⟩`（及 `a_n=⟨n\|A_L\|Ω⟩`、`f₁²=\|⟨Ω\|O\|𝔾⟩\|²`） | `confinement/PROOF.md:175,192`；`LOOP-RUNBOOK.md:4547,5055`；`confinement/REPORT.md:84` | Hilbert 空间单位向量；`[T]` |
| **Ω-3** | **真空波函数 / 真空测度** `Ω>0` 为 PF 基态函数，`μ = Ω² dU = \|Ω\|²` | `v3/species/R1-BRIDGE.md:19,99,102,106,151`；`gauge-invariance/REPORT.md:143`（`TΩ=Ω`、`真空 Ω=1`） | 构型空间上的正函数；`[T]` |
| **Ω-4** | **正交补** `Ω^⊥` / `(ℂΩ)^⊥`（谱隙的定义域） | `confinement/PROOF.md:192`（`Δ:=inf spec(H\|_{(ℂΩ)^⊥})`）；`gauge-invariance/REPORT.md:150,419`；`LOOP-RUNBOOK.md:190,5135,6388` | 子空间；`[T]` |
| **Ω-5** | **构型空间** `Ω = G^{links}` | `bianchi/NA-3-T-spec.md:14` | 流形；`[T]` |
| （旁义，异字形） | `ω` = 2-cocycle 相位 `e^{2πikq_iq_j/N}`；`ω = s/6+s³/360` | `defect-bootstrap/*`; `LOOP-RUNBOOK.md:4524` | 不冲突，但排版易混 |

⇒ **`gap2` 用的绝对是 Ω-1**。`Ω-2/3/4` 是**同一个字**在做三件不同的事，写任何一句含 `Ω` 的结论时必须显式限定（建议：`Ω_st`（state）/ `Ω_vac`（vacuum）/ `Ω^⊥`（perp）/ `Ω_cfg`（config））。这是 `LOOP-RUNBOOK.md:17899` 预判的"第 32 处撞名"，**成立**。`[T]`

---

## §4 第三问：★ 为什么"缺的是 `Ω`"

### 4.1 `Ω` 在 `G″ = S − D` 里扮演什么角色？——**它不在那条恒等式里，这正是不撞车的原因** `[T]`

`G″ = S − D`（精确方差分解，`β=2.5`，`D/S≈2%`）+ Fredholm 是**单点、单一 β** 的**无展开恒等式**：它只用到该 `β` 处测度的方差结构，因此**既不需要 `Z_O`、也不需要 `Ω`** —— 不需要搬任何东西。

这与 `gap2` 的三段式完全自洽：

```
不搬 β  ⇒  不需要 Z_O（算子协变）
        ⇒  也不需要 (4.1) 里那个跨耦合的 Ω
但  ⇒  "在 2.5 就地" 得到的是 f₁²_{lat}(2.5)，不是 f₁²_{phys}
```

**⇒ `Ω` 的角色是"换算器/同一性证书"，不是"数值常数"**：

- `Z_O` 回答"**算子**怎么变"（保形协变，§1 已精确）。
- `Ω` 回答"**这是不是同一个原子**"（态的不混合、不溶解、不交叉）。
- (4.1) 说：`f₁²` 的搬运 = **算子因子 × 态因子**。语料 60+ 轮只在做**前者**（`Z_O>0` 平凡且免费），把**后者偷偷设成 1**。

### 4.2 为什么"缺的是 `Ω` 不是数" —— 三条精确理由

1. **`0` 是吸收点（N1，`R1:218`）**：`f₁² ↦ Z_O² f₁²` 对任何 `Z_O>0` 都在 `0` 处消失 ⇒ **任何"数"（任何正标量）都搬不动一个可能为 `0` 的量**。`[NG]`
2. **`Ω` 只能从上界（N2，`R1:226`）**：Cauchy–Schwarz / 正性给 `|⟨G'|O|G⟩|² ≤ …`，Schur-凹 ⇒ **一切 `ν≥0` / RP / Schur 工具只出上界，永不出下界**。所以缺的不是"更好的数"，是**一个方向上根本不同的界**。`[NG]`
3. **`Ω` 就是隙本身（`R1:206`）**：`Ω=1 ⇔ I₂ ⇔ (CP) ⇔ m₁≥μ`。所以"缺 `Ω`"= "缺隙"，只是换了**残数服装**；这也解释了为什么它不能是"数"——**数是结论，`Ω` 是结论的另一件衣服**。`[T]`（对象等同）/ `[O]`（下界未证）

### 4.3 一句话

> `Z_O` 免费，`Ω` 是全部工作；`(P-Ω): Ω≥ω₀>0` **原封未动**（`R1:310` 逐字：*"The transport (`Z_O`) is free; `(P-Ω)` is the work, and it is unmoved."*）`[O]`

---

## §5 判：`Ω` 的定义清楚吗？

**定义清楚，且是显式可写的闭式；但"可构造到下界"不成立。**

| 层级 | 断言 | 标签 | 出处 |
|---|---|---|---|
| 定义 | `Ω(β₂,β₁)=\|⟨G(β₂)\|O\|G(β₁)⟩\|²`，显式、无自由参数 | `[T]` | `R1:189` |
| 分解式 | `f₁²_phys = Z_O²·Ω·f₁²_lat`，且 `Ω=1` 时退化为语料式 | `[T]` | `R1:188,194` |
| 对象等同 | `Ω = I₂ = (CP) = gap` | `[T]`（等同）/ `[O]`（未证） | `R1:206,328` |
| 上界 | `Ω ≤ 1`（Cauchy–Schwarz），且**只**有上界 | `[NG]` | `R1:226` |
| 下界 | `Ω ≥ ω₀ > 0` | **`[O]` 未证** | `R1:295` |
| 可测性 | "a direct lattice measurement of `Ω` with certified isolation `I₂` … **not available**" | `[O]` | `R1:306–307` |

### 裁决

- **(a) 成立（限定版）**：`Ω` 定义**清楚**、**显式可写**、**原则可算**（`Ω` 是两态重叠模方，给定两端本征态就能算）。所以 `gap2` 的引用**不是空名**，逐字见 §1。`[T]`
  - **但** (a) 的"可构造"**只到"能写出/能数值测"**，**不到"能证下界"**：`(P-Ω)` 是 `[O]`，且 §4.2 的 N2 说明**正性类工具结构上不可能**给出它。⇒ **`Ω` 是"名字给对了"，不是"路已经通了"。**
- **(b) 成立**：`Ω` 有 **≥4 个不相干的意思**（§3），是本项目第 32 处撞名。引用 `R1:4.1` 时必须写 `Ω-1/Ω_st`。
- **(c) 不成立**：`R1:4.1` **找到了**，且**就是**这个意思（`R1-RESIDUE-TRANSPORT.md §4.1`，非 `R1-BRIDGE.md §4.1`）。引文无错。

**⇒ 对 `gap2` 那句话的判决**：`"缺的是 Ω 不是数"` 是**正确且有力**的一句 —— 它把墙从"数不够好"改述成"对象没被识别"。但它**不构成路**：`Ω` 的可构造性止于"可写可测"，其下界 `(P-Ω)` 与隙**恒等**，因此在"在 `2.5` 就地"的图景里，**`G″=S−D` 绕开了 `Z_O`，却没有绕开 `Ω`** —— 只是让 `Ω` 从"跨耦合的搬运因子"退化成"**单点 2.5 上 `0⁺⁺` 道是否被 `O` 看见**"这一个更小的问句。`[C]`

---

## §6 留给下一手的精确问句（不改靶，只缩靶）

1. **缩靶后的 `Ω`**：`Ω_st(2.5) := |⟨G^{0⁺⁺}_{lightest}(2.5)|O|G^{0⁺⁺}_{lightest}(2.5)⟩|²`——**单点版**是否仍是 `[O]`？若是，则 `gap2` 的"就地"路线**等价于**"在 2.5 上测 `0⁺⁺` 道的重叠"，与 `D/S≈2%` 那条恒等式**不共享任何一步**。`[C]`
2. **是否已有 `Ω` 的数值**：`Ω` 只有上界 `≤1`（`RUNBOOK:6689` 逐字："只有上界 `Ω ≤ 1`，无下界"）⇒ 直接 lattice 测 `⟨G|O|G⟩`（`R1:306` 的 (b) 项）**是否已在 `w7` 的 `gevp_4d`/`t06_su3_gap` 里跑过**？这是可查的，不是可证的。
3. **撞名清理**：`Ω_vac / Ω_st / Ω^⊥ / Ω_cfg` 四名分立，写进 `docs/position-current.html`（本轮受"只写 OM1.md"约束，未动）。

---

## §7 本轮状态

- **新增（本轮唯一）**：`Ω` 的 5 义撞名表（§3）+ `Ω` 在 `G″=S−D` 中**不在场**这一结构事实（§4.1）+ 裁决（a-限定 + b，非 c）。
- **未动**：任何数值、任何 `[T]` 的既有断言、`gap2` 的结论。
- **未合 `AGENTS §6`**：`docs/position-current.html` 未更新（用户限"只写 OM1.md"）；此为**本轮唯一**未合项。

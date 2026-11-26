# Z09-HOLED — 洞 D 裁决：`O` 的「密度 / 整片求和」二义、`f₁²` 的广延性、`a₀` 的 `L`-无关性

> 2026-11-25 · 只处理洞 D（`Z05-AUDITI1.md:78,91`）。**判据 = `docs/paper-v2/ch05-theorems.tex:60`。**
> **源（实读，非转述）**：`docs/paper-v2/ch05-theorems.tex`（:20–99、:415–470、全文 grep）· `research/unconditional/P3-DRESSED-MASS.md`（:14–79、:110–139）· `P3-CYLINDER-RESIDUE.md`（:20–74）· `P3-FOURWALL-DERIVATION.md`（:74–87）· `PROOF-P3-UNCONDITIONAL.md`（全 36 行）· `research/v3/species/U04-MEASF1.md`（:57–72、:203–211）· `U08-CORR.md`（:37–53）· `research/v3/R4-RESIDUE-PROTOCOL.md`（:74–83）· `LOOP-RUNBOOK.md`（:6205–6235）。
> 检索核验见 §8。标签：`[T]` 真/逐字事实 · `[C]` 条件 · `[R]` 缺失步 · `[M]` 数值 · `[J]` 判断。

---

## 0 · 裁决（一句话）

**洞 D = 真现象 + 假诊断。**

- **现象真**：`O` 在**同一语料里被写成三档**（密度 / 整片求和 / 归一化平均），三档差**恰好一个 `V_s=L³` 的幂**；按主链文件（`P3-DRESSED-MASS`、`P3-FOURWALL`、`ch05`）的记账，`f₁²_L` **是广延量 `∝V_s`**，而 `a₀`、`r₀` 是 `L`-无关的强度量 ⇒ **两个不同归一化被直接比大小**。
- **诊断假**：`a₀` 的 `L`-无关性**已被证明**（`P3-DRESSED-MASS.md:38` 几何和、`:79` "finite backbone positions"，`:77` wrap `e^{−cL/a}`）；Z05 的"`a₀` 的 `L`-无关性未证"把**归一化缺失**误报成**解析缺失**。

**判死线：走 (c)**（两定义并存，分配见 §1/§2）+ **按 (b) 的后果处理**（`r₀` 必须带 `V_s`，或 `f₁²` 必须归一化）。**不是 (a)**——洞不"消失"，裸写的定理陈述仍是错的。`[J]`

---

## 1 · ★ 逐字核准：`O` 的五处定义（`[T]`，原文抄）

| # | 文件:行 | 逐字 | 归一化档（`L` 幂） |
|---|---|---|---|
| 1 | `PROOF-P3-UNCONDITIONAL.md:11` | "the lightest `0⁺⁺` atom of **`O=¼tr F²`**" | **标签=密度**（未给求和；幂 0，强度） |
| 2 | `P3-DRESSED-MASS.md:25` | "`O(x) = Σ_{spatial plaquettes p ⊂ x-slice} (2N)^{-1} Re tr U_p = Σ_p (2N)^{-1}(χ_f(U_p)+χ_{\bar f}(U_p))`" | **整片求和**（`(x)` 是**时间片**号，不是格点；幂 +1，广延） |
| 3 | `P3-CYLINDER-RESIDUE.md:61` | "Let `O := (2N)^{−1} Σ_{p∈A1g}(χ_f+χ_{\bar f})(U_p)` be the `A1g`/zero-momentum gauge-invariant **single-plaquette density** on one time slice" | **同句二义**：字面 density、"zero-momentum" 要的却是片内求和 |
| 4 | `P3-FOURWALL-DERIVATION.md:79` | "`C(t)=Σ_x Σ_{y∈spatial} ⟨O(x,0)O(y,t)⟩_c` (zero-momentum sum over sink; source fixed by translation invariance)" | **双空间求和**（幂 +1，广延） |
| 5 | `U04-MEASF1.md:57` / `:67`；`U08-CORR.md:37` / `:49` | 逐点 "`O(x) ≡ ¼ tr F_{μν}F_{μν}`"；片算符 "`O(t) = (1/V_s) Σ_{x_s} O(x_s,t)`, `V_s = L³`" | 密度（幂 0）/ **归一化平均**（幂 −1） |

**结论 `[T]`**：`O` 的定义**不是一个，是三档**；`Z05` 写的 `O = Σ_x ¼trF²(x)` 只是 #2 的口语化，而 `U04/U08` 用的是 #5（两种都写全了）。

---

## 2 · ★ 判据文件 `ch05` 内部：把两档缝在同一节

| 行 | 逐字 | 属哪档 |
|---|---|---|
| `:32` | `f_{1,L}^2:=\nu_L(\{E_{1,L}\})=\lim_{t\to\infty}e^{E_{1,L}t}C_L(t)` | 残数（未指定 `O`） |
| `:52` | "`\operatorname{tr}F^2` is a finite A1g sum of local plaquette observables, so `C(t)` is a sum over connected polymers" | **求和档**（#2/#4） |
| `:60` | "`f_1^2=|\langle\Omega|O|\mathbb{G}\rangle|^2` is a modulus-square, so no sign can undo it." | **强度档**（模方 + 与 `:58` 的 `L`-无关 `a₀` 并列） |
| `:62` | "the remainder is `O(u)` **uniformly in `t,V,L`**" | 求和档 |
| `:64` | "A genuine **volume- and `L`-uniform** lower bound on the residue `f_1^2`" | 求和档 |
| `:532` | "`f_1^2\ge r_0(\beta)>0`, **uniform in `V,L`**"；条件栏仍写 `T^4_L` | 求和档（且标签未同步） |

**两式不能同时成立**：若 `O` 是 `:52` 的片内 A1g 求和，则 `|⟨Ω|O|𝔾⟩|² ∝ V_s`，`:64/:532` 的 "**volume-uniform**" 对增长量而言**平凡**（无信息）；若 `O` 是密度，则 `:60` 与 `a₀` 的 `L`-无关性成立，但 `:52` 的求和记账、`P3-DRESSED-MASS.md:38` 的横移求和**都无对象**。

**判据侧硬事实 `[T]`**：`ch05-theorems.tex` 全文检索 `V_s` / `\sum_{x}` / `per volume` —— **命中 0**（§8）。⇒ 判据文件**根本没做**这次归一化，`f₁²` 被当成既广延又 `V`-一致的量使用。

---

## 3 · `f₁²` 是广延量：`f₁²_L = V_s·a₀(1+O(u))`（`[T]` mod 环 1/2/4 既有输入）

管权 `w(B_T)=a₀u^{4T}`（`P3-DRESSED-MASS.md:36`），而 `B_T` 的截面是**固定**的源空间面 `p₀`（`:44` "cross-section the source spatial plaquette `p₀`"）。片算符要对 `p₀` 求和，平移不变 ⇒ 每个 `p₀` 同权：

```
C_L(t) = Σ_{p₀} Σ_{p_t} w(tube) = V_s · a₀ u^{4T} (1+O(u))     ⇒   f₁²_L = V_s · a₀ (1+O(u)),   V_s=L³
```

- `P3-DRESSED-MASS.md:38` 的 `Σ_x u^{|x|_link}=O(1)` 是**相对位移**（kink 链）之和——正确；**缺的是绝对位置和 `V_s`**。
- `A1g` 面数（3 或 6）是 `O(1)`，并入 `a₀`；幂次结论 `∝V_s=L³` 不受影响。`[T]`

### 语料三处**独立**同判（逐字，非推演）

1. `R4-RESIDUE-PROTOCOL.md:81–82`：**"The physical (per-volume) residue is `f₁²/V`; on the zero-momentum operator the amplitude is extensive, so a genuine atom has `f₁² ∝ V` while the mass `m_G` is volume-stable."** `[T]`
2. `U04-MEASF1.md:208–209`：`O(t)` 用 `1/V_s` 平均，"`C(0) ∝ 1/V_s` ⇒ `f₁²_avg` **不是**物理残数（物理残数对应 `Σ_{x_s}` 求和、`f²_sum = V_s² f²_avg`）"。`[T]`
3. `U08-CORR.md:51–52`：照抄同一警告，"物理残数对应 `Σ_{x_s}` 求和"。`[T]`

**且 `P3-DRESSED-MASS.md:79` 自己写了解法**：
> "The finitely many spatial positions of the backbone (folded into `a₀`, §0) are `L`-independent **up to the overall `L³` measure that cancels in the normalized residue**." `[T]`

⇒ 那个 `V_s` **已被点名**（`:79`），但**没有进任何定理陈述**：`PROOF:11–15`、`P3-CYLINDER:73–74`、`ch05:32/60/64/532` 一律裸写 `f₁²_L`。**洞 D 的全部内容就是这一句没上升为定义。**`[J]`

---

## 4 · `a₀` 的 `L`-无关性：**已证，不是洞**

| 步骤 | 内容 | 标签 |
|---|---|---|
| `P3-DRESSED-MASS.md:38` | 横移（kink）求和 `Σ_x u^{|x|_link}` 几何收敛 ⇒ `O(1)`；并**显式**排除熵解读："it is **not** an entropy `μ^T` competing with the energy" | `[T]`（每移一格 ≥1 个 `u`；指数是约定，任何 `u^{c|x|}` 均几何收敛） |
| `P3-DRESSED-MASS.md:79` | backbone 位置 "finitely many … `L`-independent up to the overall `L³`" | `[T]` |
| `P3-DRESSED-MASS.md:77` | wrap 修饰 `u^{L/a}=e^{−cL}` ⇒ `f_dec(L)=f_dec(∞)+O(e^{−cL/a})` | `[T]` |
| 干净上界（独立于 `L`） | `a₀(L)=Σ_{Δ∈T³_L}u^{c|Δ|}(recoupling) ≤ Σ_{Δ∈ℤ³}u^{c|Δ|}<∞`，且 `≥` `Δ=0` 项 `>0` | `[T]`（正项级数可延拓到无限格） |

**⇒ Z05 §3（`:78`、表 D 的"`a₀` 的 `L`-无关性未证"）为假**：真缺陷只有"`f₁²_L` 与 `a₀` 不属同一归一化"。`[J]`

---

## 5 · 与 `r₀ = ½c₁u^p` 的关系（任务 2）

`r₀` 的两个版本**都 `L`-无关**：`½c₁u^p`（`PROOF:31`）与 `¼(2N)^{-2}`（`P3-DRESSED-MASS:124`）。而 `a₀=(2N)^{-p}s`（`:119`）同样 `L`-无关 ⇒ **`r₀` 属强度/归一化档**（裸端点正系数乘积的量级；`p` 在两处指不同东西，那是**洞 E**，与 D 正交）。`[T]`

于是 `f₁²_L ≥ r₀` 在广延归一化下是

```
V_s · a₀ ≥ r₀   ⟺   V_s ≥ r₀/a₀ = O(1)
```

——**只要 `L` 够大就平凡成立**（任务 2 的观察成立：`r₀` 足够小即平凡）。后果：`inf_L` 落在 `L_min`，**"`L`-一致"退化为单点**，而 `Σ_x` 的 `L³` 从未显示。

**要与 `r₀=½c₁u^p` 相容，只有两条路，必须二选一明写**：

- **(广延版)** `r₀^{ext}(β,L) := V_s·½c₁u^p = L³·½c₁u^p` —— `r₀` 带 `V_s`，"`r₀(β)` 与 `L` 无关"这句**作废**；
- **(强度版)** `f̄₁² := f₁²_L/V_s ≥ ½c₁u^p` —— 此时 `a₀` 的 `L`-无关性**正是**要证的强度量性质，且已证（§4）。

两版数学等价（同一条不等式乘/除 `V_s`）；**不可写的只有现文这一种**：裸 `f₁²_L` 对 `L`-无关 `r₀`。`[T]`

---

## 6 · 对 `I1` 的影响（三条）

1. **真值不变**：`a₀>0`、KP（环 2）、`[KL]`（环 4）、`t→∞` at fixed `L` 全不动。归一化是**记账**，不进 `[C]→[T]` 的障碍集（A/B 才是硬伤）。`[J]`
2. **陈述必须改写**三处全裸写的位置：`PROOF:13–15`、`ch05:32/60/64/532`、`P3-CYLINDER:73–74`。（这是"修洞最小集"第 (3) 条的可执行形式。）
3. **"`L`-一致的来源"被认错**：归一化后 `f̄₁²(L)=a₀(L)(1+O(u))`，其 `L`-一致性来自 **`a₀` 的正项几何和**（§4，`[T]`，**不依赖** wrap 项）；wrap `e^{−cL/a}` 管的是 `f_dec`/质量，不是振幅地板。Z05 把 `a₀` 一支判为"方向可能反"，**方向反了**。`[J]`

---

## 7 · 判（判死线逐条）

- **(a) 不成立**：不能说"只需更正 `U04`/`U08`"。`U04/U08` **已经**写对了（密度 + 归一化 + 警告），错的是 `DR/FW/ch05` 的记账与陈述。
- **(c) 成立（两定义并存）**，分配已钉死：
  **标签档 = 密度**（`PROOF:11`；`CYLINDER:61` 的字面；`U04:57`；`U08:37`）·
  **记账档 = 整片求和**（`DR:25`；`FW:79`；`CYLINDER:61` 的 `Σ_{p∈A1g}`；`ch05:52`）·
  **数值档 = 归一化平均**（`U04:67`；`U08:49`）。
- **(b) 的后果对主链成立**：`DR/FW/ch05` 按整片求和记 ⇒ `f₁² ∝ V_s` ⇒ **`r₀` 必须依赖 `V_s`**（§5 两条路择一），否则 `I1` 的"`L`-一致"平凡化。
- **是否真洞**：**现象真、诊断假**。真=差一个 `V_s`；假=`a₀` 的 `L`-无关性未证。**修法一行**：给 `f₁²_L` 加 `1/V_s`（或给 `r₀` 乘 `L³`），并把 `P3-DRESSED-MASS.md:79` 的"cancels in the normalized residue"提升为**定义**。`[J]`
- **不是 `[NG]`**：无一步给出 `f₁²/V_s→0`。`[J]`
- **同类漂移旁证**：`ch05:532` 摘要表仍写 `T^4_L`（定理在柱面 `T³_L×ℝ`）—— 标签/记账未同步，与 D 同源。`[T]`

---

## 8 · 检索核验（AGENTS §1：工具状态 + 检索式 + 命中数）

| 检索式 | 目标 | 工具状态 | 命中 |
|---|---|---|---|
| `V_s\|\\sum_\{x\|per volume\|volume factor\|residue` | `docs/paper-v2/ch05-theorems.tex` | grep **OK** | 26 行，其中 `V_s`/`\sum_{x}`/`per volume` **0** ⇒ ch05 无归一化 `[采信]` |
| `\\frac\{1\}\{4\}\|frac14\|⅓¼ ?tr` | 同上 | grep **OK** | **0** ⇒ ch05 不写算符显式式 `[采信]` |
| `O\(x\)\|Σ_|密度\|density` | `U04` / `U08` / `DR` / `FW` | grep **OK** | 10 / 14 / 4 / 多（逐字已抄）`[采信]` |
| `广延\|强度量\|extensive\|intensive` | `research/**/*.md` | grep **OK** | 208（R4`:81–82`、U04`:208–209` 逐字）`[采信]` |

**无 `TOOL-FAIL`。** 本轮未读 `Y06`/`Z01`/`Y05` 的结论（只读原始文件与其逐字内容）。

*源路径（绝对）*：`D:\Dev\math\docs\paper-v2\ch05-theorems.tex`、`D:\Dev\math\research\unconditional\{PROOF-P3-UNCONDITIONAL,P3-DRESSED-MASS,P3-CYLINDER-RESIDUE,P3-FOURWALL-DERIVATION}.md`、`D:\Dev\math\research\v3\species\{U04-MEASF1,U08-CORR,Z05-AUDITI1,LOOP-RUNBOOK}.md`、`D:\Dev\math\research\v3\R4-RESIDUE-PROTOCOL.md`。
*无 git commit。*

# GPP-STRUCT — `G''(0)` 是什么对象：结构类型与四条候选不等式的判决

**日期**: 2026-11-25 · **状态**: 已判定（结构类型 [T]；四条候选 [NG]）· **禁 git commit**（遵守）· 只写本文件

---

## 0. 判死线（先读这段）

**`G''(0)` 是"条件化（钉扎）自由能的二阶导"**：`G(a) = −log ρ_v(a)`，`ρ_v` 是粗模 `ξ=⟨v,X⟩` 的条件边缘密度在平直图 `dU=J dX`、`V=S_W+Φ`、`Φ=−log J` 下的密度；`G''(0)` 是它在 `a=0` 的曲率。展开式是**差**：`E[∂_v²V|0] − Var(∂_vV|0)` = **条件平均曲率 − 条件力方差**。它**不是方差**，也**不是 −方差**。

- **(a) 前提不成立，条款不触发。** `G''` 含正部 `E[H_vv+Φ''|0]`，而 `Φ''⪰0` 点态（Lemma B：`E[Φ'']≥N/12`），所以 `δ≤0` 恒成立是假的：`δ>0` 有四处精确值（`+0.135995` 单 plaquette `β=2.5,θ=π`；`+0.2909…+2.587` 单 plaquette `θ=0`；`+2.0101` 2D YM；`+0.22196` P4 周期）。⇒ **H1 不能"只以 `δ≤0` 成立"**；它的内容就是 `δ>0` 的定量下界，见 §4(a)。
- **(b) 它是差，`≥0` 不平凡；但"某条已知不等式定它的号"是否定的。** 四条候选里三条（Chebyshev/FKG/Griffiths）方向**反**（给 `Cov≥0` 下界，而 `G''≥0` 要**方差上界**）；唯一方向对的是 **Ginibre/Lebowitz 上界族**，其前提是**正四次项 / `Hess V⪰0`**——正是 Savvidy 墙（`λ_min(Hess S_W)=−g|B|<0`）。⇒ **H1 的真身 = 1D 边缘的 log-凹/BL 型"条件方差 ≤ 条件曲率"不等式**，φ⁴ 上有（Lebowitz/Ginibre），SU(N) Wilson 上缺席。不是四条候选中的任何一条。
- **(c) 恒等式本身是换名（`[NG]` 作为路线），但"`G''≥0` ⟺ H1"作为陈述不成立**：精确等价的是**整族** `H1(δ) ⟺ G''≥δ`（[T]），而 `G''≥0` 是 `δ=0` 成员，**严格弱于**物理点 `H1(δ=a²Δ₀²+W_full>0)`（`W_full=W_geo+X`）。判：**凸区已证 [T]（`β=0`：`G''≥N/12>0`）、某处已证假 [T]（`β>β*=4.30661`）、物理点 `β=2.5` 未证（可证未证，缺的工具不在这四条里）**。

---

## 1. `G` 的精确定义（逐字核准，不凭记忆）

### 1.1 论文原文（两处，互相独立）

`docs/paper-v2/ch05-theorems.tex:230–240`（式 `eq:split`，标签 `\label{eq:split}`）：

```
G''(0) = E[H_vv + Φ'' | ξ=0] − Var(G_v + Φ' | ξ=0) = Δ_dyn + W_geo + X,
Δ_dyn = E[H] − Var(G),      X = −2 Cov(G, Φ').
```

`docs/paper-v2/ch15-uvface.tex:20–28`（式 `eq:v3object`，**给出对象的名字**）：

```
G''(v;a₀) = E_a[∂_v²V] − Var_a(∂_vV) = Δ_dyn + W_geo + X,   X = −2Cov(G,Φ'),
V = S_W + Φ,  Φ = −log J 的 Haar Jacobian,  条件化在粗场 a = a₀.
```
同处 `ch15:20–21` 明写该对象是 **“the coarse-field-conditioned marginal curvature”**。
`ch08-immune.tex:138` 明写它是 **“the Poincaré constant of the 1D marginal `=1/χ=G''`=H1”**。
`app.tex:173–175` 明写 **`G''·Var(ξ)=1` 当且仅当边缘是 Gaussian**（实测偏离 `2–5%`，不得当恒等式用）。

### 1.2 因此 `G` 是什么（可核的写法）

设 `dμ = e^{−V(X)}dX`（`V=S_W+Φ`，`Φ=−log J`，`dU=J dX`），粗坐标 `ξ(X)=⟨v,X⟩`，条件边缘密度

```
ρ_v(a) := ∫ e^{−V(X)} δ(⟨v,X⟩ − a) dX ,        G(v;a) := −log ρ_v(a) .
```

则（把 `Z(a)=ρ_v(a)` 对 `a` 求导两次，标准条件化公式）

```
G''(v;a) = ‖v‖⁻⁴·{ E[ ∂_v²V | ξ=a ] − Var( ∂_vV | ξ=a ) } = ‖v‖⁻⁴·{ E[H_vv+Φ''|a] − Var(G_v+Φ'|a) } ,   （`∂_v:=⟨v,∇⟩` 未归一；`‖v‖=1` 时即 `eq:split` 逐项一致，`‖v‖²=1/L` 的 `L`-link 平均下差 `L²`）
```

与 `eq:split` / `eq:v3object` 逐项一致。**`G` 本身只以二阶导出场**：`G` 可差一个仿射项，论文从未写出 `G(a)`，只写 `G''(0)`；可核的定义是 `G''`，不是 `G`。

**记号冲突（必须记，防止复用错误）**：`eq:split` 里的 **`G_v` 是力** `=∂_vS_W=⟨v,∇S_W⟩`，自由能是 **`G(a)`**；`H1A-EXACT.md §1.2` 把力直接写成 `G`（"`G=⟨v,∇S_W⟩`"，`X=−2Cov(G,Φ')`）。同一字母 `G` 在两个文件里指**两个不同对象**。写任何式子前先分清：**`G(a)` = 自由能，`G_v` = 力**。

### 1.3 `G''(0)` 的精确表达式

```
G''(0) = ∂²_a[−log ρ_v(a)]_{a=0}
       = E_μ[H_vv + Φ'' | ξ=0] − Var_μ(G_v + Φ' | ξ=0)          （eq:split 第一行，[T]）
       = Δ_dyn + W_geo + X ,                                      （eq:split 第二行，[T]）
  Δ_dyn = E[H|0] − Var(G|0)，  W_geo = E[Φ''|0] − Var(Φ'|0)，  X = −2Cov(G_v,Φ'|0)。
```

`H := Hess V`（**不是** Hamiltonian），`' := ∂_v`，`v` = 粗横向 `k→0` 模，`ξ=0` = 条件化点。三个块中只有 `W_geo` 有已证符号：`W_geo≥0`（Prékopa，纯 Haar 边缘；`eq:wgeo`，`[T]`；点态 `Φ''⪰0` 给 Lemma B `E[Φ'']≥N/12`）。`Δ_dyn` 与 `X` 符号不定（`X∈[−0.09,+0.16]` 逐模变号，`β=5.7` 翻号）。

---

## 2. 它是方差 / 协方差 / 自由能二阶导 / 差？

**答：自由能二阶导，且写成差（条件平均曲率 − 条件方差）。不是方差，不是协方差，也不是 −方差。**

| 读法 | 是否成立 | 判据（精确反例，均来自 `H1A-EXACT.md` §2，数值已与存档对 6 位） |
|---|---|---|
| 方差型（天然 `≥0`） | **否** | `β=6, θ=π`：`G''=−0.193019<0`（`Var=0.398135`、`E[H]=0`、`W_full=0.205115`） |
| `−方差`型（天然 `≤0`） | **否** | `β=2.5, θ=π`：`+0.135995`；`θ=0, β=2.5`：`+0.729924`；2D YM：`+2.0101`；P4 `L→∞`：`+0.22196` |
| 协方差型 | **否** | `G''` 是曲率，不是 `Cov(f,g)`；`X` 才是协方差（且只是加项之一） |
| 自由能二阶导 | **是** | `G=−log ρ_v`；`G''=∂²_a(−log ρ_v)`；Gaussian 情形 `G''Var(ξ)=1`（`app.tex:174`） |
| 差 `E[Hess]−Var` | **是** | `eq:split` 第一行；`−Var` 只吃**一部分**，正部 `E[H_vv+Φ''\|0]` 可正可负（`θ=π` 时精确 `0`，`θ=0` 时 `+0.553`） |

**结构等式（本文件的判定核心）**

```
G''(0) ≥ 0   ⟺   Var_μ(∂_vV | ξ=0)  ≤  E_μ(∂_v²V | ξ=0)
             ⟺   条件边缘 ρ_v 在 a=0 处局部 log-凹（BL / Poincaré，λ=0 的弱形式）
             ⟺   G 在 0 处凸（局部）；"G 不凸"就是 G''<0，不是"G'' 是 −Var"
```

Gaussian 边缘是等式情形（`G''=1/Var(ξ)`），所以 `G''` 的正确名字是**局部精度 / log-凹模量 / 粗模有效质量平方**，而不是任何方差。`H1A-EXACT.md` 记的"`G` 不凸"（单 plaquette `β>4.30661`）正是这条局部 log-凹性的失效，**不能**升级为"`G''` 型为 `−Var`"——因为那样 `θ=0` 与 2D YM 的四个正值无法解释。

---

## 3. 四条候选不等式，逐个判

判据统一为：**要 `G''≥0` 就必须给出"条件方差的上界"（或"条件平均曲率的下界"）；任何只给 `Cov≥0` 的不等式方向反。**

| # | 候选 | 它到底给什么 | 在什么测度上成立 | 对 SU(N) Wilson 已知？ | 能否给 `G''≥0` |
|---|---|---|---|---|---|
| 1 | **Chebyshev 关联不等式** | 同向单调函数 `Cov(f,g)≥0`（单实变量无假设；多维需**共同标量/全序**，重排 = Fréchet–Hoeffding = **Cauchy–Schwarz**） | 任意 1D 测度（平凡）；多维要格序 | **不适用**：纤维 `2(N²−1)V` 维，无共同标量、无全序；`Chebyshev/TP₂ 要共同标量（无）`（`C6-DEPENDENCE-CLASS.md:19`、`C5-NEG-ASSOC.md:16–17`） | **不能**。方向：`Cov` 下界。多维版退化为 CS，给的是**双边上界** `\|X\|≤2√(VarG·VarΦ')≈0.40` 对 slack `0.32` ⇒ 只到 `G''≥−0.08`（**窗口：`SU(3)`、`β=5.7`、`L=8`（CS 退化，非 `β=2.5`）；`W_geo=+0.318±0.006`**；`ABC2-CONSTRUCTIVE.md §4`、`H5-HEAD-NUMERICS.md:156`、`B4-VARIANCE-MAG.md:15`；`docs/PROOF-SKELETON.md:44` 逐字"Chebyshev/Fréchet=CS"） |
| 2 | **FKG / 正关联** | 有限**分配格**上增函数正相关：`Cov(f,g)≥0` | 乘积链 + log-超模（铁磁）密度 | **未证/前提空**：`SU(N)` 无平移不变格序（torsion）、`e^{−S_W}` 非逐坐标超模、Weingarten 面权带号（`C6:17,31`；`C8-SPECIES-NOGO.md:6`；`PROOF-COMPARISON-FKG §1`） | **不能**。且方向反：`FKG/Griffiths: Cov≥0 ⇒ X≤0 ⇒ G''≤W_geo+Δ_dyn`——**降** `G''`（`C6:17,31` 原文"WRONG direction"） |
| 2′ | **RP（格点规范确实有，但 `RP≠FKG`）** | `⟨θF,F⟩≥0`、转移矩阵正性、OS 重构 | **任意紧群、任意 `β≥0`**（Osterwalder–Seiler 1978） | **是**（定理级），但它是**反射内积正定**，不是协方差不等式 | **不能**。RP 的输出是谱正性/上界，语料已钉："positivity is blind to the edge"（`docs/YANG-MILLS-CONDITIONAL-REDUCTION.md:480`）。`RP` 与 `FKG` 互不蕴含：前者是双线性型正定，后者是协方差号；把 RP 当 FKG 用是**类型错误** |
| 3 | **Griffiths 不等式（GKS）** | `⟨σ_Aσ_B⟩ ≥ ⟨σ_A⟩⟨σ_B⟩`（GKS 1970）、`⟨σ_A⟩≥0`（Griffiths 1967） | 铁磁 Ising / 向量自旋，`J≥0` | **未证/前提空**：Wilson 面权对 `N≥3` 特征标系数变号，无铁磁耦合结构（`C6:17`；`Q4-LOGCONVEX.md:46` 引 `GKS-GAUGE §4`——**该文件在仓库中不存在**，见 §5 悬空引用） | **不能**。同样是 `Cov` 下界，方向反 |
| 4 | **Ginibre / Lebowitz** | **上界**：`m=2` 时 `U₄=⟨φ⁴⟩−3⟨φ²⟩²≤0`（Lebowitz 不等式）；一般 `n`：截断函数 `\|S_n\|` 上界；Ginibre 型 `⟨:φ²:(x):φ²:(y)⟩^T ≤ 2⟨φ(x)φ(y)⟩²` | **`λφ⁴`（`λ≥0`）**、铁磁向量自旋（Ginibre 的 Griffiths 普遍形式） | **部分**：**abelian** 格点规范（Villain，经 Wilson/Manton 极限）的 Ginibre 关联不等式**已证**（Chevyrev–Garban）；**非阿贝尔 SU(N) 未查到（命中 0）** | **不能，但它是唯一方向对的一条**。它的前提就是"正四次项/凸势"，即 `Hess V⪰0` = Savvidy 墙（`λ_min(Hess S_W)=−g\|B\|<0`，`H1E-SAVVIDY`/`PROOF-SAVVIDY-WGEO-DOMINATION`） |

**四条合并的结论**：1、2、2′、3 是**正关联族**（`Cov≥0`），方向与 `G''≥0` 相反，**原理上不可能**给出所需方差上界；4 是**上界族**，方向对，但它的可证性依赖凸性，而凸性正是墙。⇒ **`G''≥0` 不对应任何一条已知相关不等式；它对应的是"条件方差 ≤ 条件曲率"的 BL/log-凹型不等式**——这与 `docs/YANG-MILLS-CONDITIONAL-REDUCTION.md:585` 记的族 14"lower-bound machine, upper-bound target"同一条，也与 `docs/PROOF-SKELETON.md:44` 的"BL 需 `Hess⪰0`=墙；FKG 反号"逐字一致。

### 出处（可核）

- Chebyshev 积分不等式：**Hardy–Littlewood–Pólya, _Inequalities_, 2nd ed., Thm 236**（`[T-standard]`；**书目级，未逐字核**）。概率版（Harris 型）：`ABC-NEW-INEQ`/`ABC-FINAL-LEDGER` 参考表列 "Harris 1960"（**未逐字核**）。
- FKG：**Fortuin–Kasteleyn–Ginibre, Comm. Math. Phys. 22 (1971) 89–103**，"Correlation inequalities on some partially ordered sets"。**摘要级已核**（原摘要："increasing functions on a finite distributive lattice are positively correlated by positive measures satisfying a suitable convexity property"）。
- Griffiths：**J. Math. Phys. 8 (1967) 478**（II 篇 484）；GKS：**Griffiths–Hurst–Sherman, J. Math. Phys. 11 (1970) 790**（`[T-standard]`；**书目级，DOI 未逐条核**）。
- Lebowitz：**Comm. Math. Phys. 35 (1974) 87–92**，"GHS and other inequalities"，DOI `10.1007/BF01646608`（NSTL 书目页）／另一书目页记 `10.1007/BF01646068`——**两个 DOI 不一致，标 [E]，未逐字核**。`U₄≤0` 即 Lebowitz 不等式、`m>2` 无严格结果：**摘要级已核**（J-STAGE 检索命中原文）。
- Ginibre：**Comm. Math. Phys. 16 (1970) 310**（"General formulation of Griffiths' inequalities"）；不等式陈述见 **Simon, Caltech preprint（"Theorem 2.4 (Ginibre Theorem 4)"）** 与 **arXiv:1905.09605, Prop. 11.3（Lebowitz inequality）**，`[C]` 摘要/正文片段级。
- 非阿贝尔/abelian 分野：**Chevyrev–Garban, "Villain action in lattice gauge theory", arXiv:2404.09928**（J. Stat. Phys. 2025）——"In the special case of **Abelian** symmetry … the validity of **Ginibre correlation inequality** for lattice gauge theory with Villain interaction"。**摘要级已核**。
- RP：**Osterwalder–Seiler, Ann. Phys. 110 (1978) 440–471, DOI 10.1016/0003-4916(78)90039-8**（**书目级已核**；内部评估 `research/lit-v3/LIT-ERA-04.md:16,198`）。

---

## 4. 判：已证 / 可证未证 / 循环

### (a) `G''=−Var` 型？——**不触发**，且要写清"若不是那它是什么"

`G''` 不是 `−Var`：它有正部 `E[H_vv+Φ''|0]`。`Φ''⪰0` 点态（`[T]`，Lemma B），`H_vv` 可正可负但**可为正**（`θ=0` 切片 `E[H]=0.553`…`2.419`；2D YM 的 `E[H]=W_geo−δ`）。所以：

> **`δ≤0` 恒成立是假的。H1 的内容 = "`G''(0) ≥ a²Δ₀²+W_full > 0`"（`δ:=a²Δ₀²+W_full`，`W_full=W_geo+X`），即两件事同时成立：(i) 条件化点落在**稳定真空**上（否则 `E[H|0]→0` 甚至为负，如 `θ=π` 的精确 `E[H]≡0`），(ii) **力方差被条件曲率压住**（`Var(G_v+Φ'|0) ≤ E[H_vv+Φ''|0] − a²Δ₀² − W_full`）。两半都是动力学的；没有一半是"自由的方差估计"。**

这也解释了 `H1A-EXACT.md` 的证伪机制：`δ<0 ⟺ 真空不在条件化点`（Savvidy 凹性的精确 3 维版本），**不是**"方差天然为负"。

### (b) 它是差 ⇒ `≥0` 不平凡 ⇒ 那条"真身"不等式是什么？

`G'' = E[Hess] − Var`，不是 `Var(A)−Var(B)`（`E[Hess]` 不是方差；无条件 Stein 恒等式 `E[Hess V]=E[∇V⊗∇V]⪰0` 只对**无条件**期望成立，`app.tex:167–172`，对条件化后的 `E[H|0]` 无效——这正是"条件化修正"是缺口的原因）。

**若一条已知不等式能定它的号，那条不等式就是 H1 的真身。判定：这样的不等式在已知名录里不存在，但有唯一"正确的类型"。**

- 正确的类型 = **BL / log-凹 / 上界族**：`Var(∂_vV|0) ≤ E[∂_v²V|0]`。它的凸性版本是 Brascamp–Lieb（`Hess V ⪰ λ>0`），它的 φ⁴ 可计算版本是 **Lebowitz/Ginibre 的上界**。两者都要求"势凸/正四次项"。
- Wilson 作用量在 `β=2.5` 上**不凸**：`λ_min(Hess S_W)=−g|B|<0`（Savvidy，实测 7.1% 负模；`CLAUDE.md §5`、`PROOF-SAVVIDY-WGEO-DOMINATION`）。⇒ BL 前提假，Ginibre/Lebowitz 无对应物。
- ⇒ **(b) 的答案：H1 的真身 = "1D 条件边缘的 BL/log-凹不等式"，即第六物种；它不是四条候选中的任何一条，且它在 SU(N) Wilson 上的**非阿贝尔**版本未查到（命中 0，见 §5）。**

### (c) `G''≥0` 等价于 H1？

分两层，必须分开判：

1. **作为恒等式/换名**：`H1(δ) ⟺ G''(0)≥δ` 是**代数恒等**（`[T]`，`eq:scalar-rule`，`ch18:42–48`），逐点、逐模、不依赖模型。作为**路线**，它**不产生新估计** ⇒ `[NG]`（`ch18:59–63` 的"no independent variance inequality"就是这条）。这一点上说"换名"是对的。
2. **作为陈述**：`G''≥0` **不等于** H1。`H1(δ)` 的物理点是 `δ=a²Δ₀²+W_full>0`（`app.tex:190–191` 逐字 `δ:=a²Δ₀²+W_full`；`W_full=W_geo+X`，`app.tex:796`；另见 `app.tex:182–202`，`app:sharp`）；`G''≥0` 是族里 `δ=0` 的那一个成员，**严格更弱**。故"`G''≥0` ⟺ H1"**不成立**；成立的是"整族 `⟺` 对 `G''` 的整族下界"。

**最终判词**：

| 陈述 | 判 |
|---|---|
| `G''(0)` 是自由能二阶导 = 差 `E[Hess]−Var` | `[T]`（`eq:split`/`eq:v3object`，逐字核） |
| `G''≥0` 的**类型** = 条件方差上界 = BL/log-凹 | `[T]`（恒等式改写） |
| `G''≥0` 在凸区（`β=0`，`V=Φ`） | `[T]`：`G''≥N/12>0`（Prékopa + Lemma B；`ch05:209–212`） |
| `G''≥0` 对一切 `(β,真空)` | `[T]` **假**：`β>β*=4.30661` 时 `−0.193019` |
| `G''≥0` 在物理点 `β=2.5`（`SU(N)`，`(L,a)`-一致） | **可证未证**；缺的工具 = 非阿贝尔 BL/Ginibre 型上界（未查到，命中 0） |
| "`G''≥0` 就是 H1" | `[NG]`（δ=0 成员，弱于物理点；作为路线是换名） |

---

## 5. AGENTS.md §1 检索执行核验

**工具状态**：`pwsh`（`Get-ChildItem`/`Select-String`）= **OK**（返回目录与带行号命中）；DSH `grep` = **OK**（全仓 435 命中，250 条落盘到 spill 文件，可复读）；`read` = **OK**；`web_search` = **OK**（8 条查询，均返回有效外部源）；`web_fetch` = **未调用**（外部断言全部标"摘要级/书目级"，未声称读到正文）；`glob` = 未调用（用 `pwsh Get-ChildItem` 等价替代，OK）。**无 `TOOL-FAIL`。**

| # | 检索式 | 范围 | 命中 | 结论 |
|---|---|---|---|---|
| 1 | `G''\(0\)\|G\(a\)\|G_v\|pinning\|constrained free energy\|free energy` | `docs/paper-v2/*.tex` | 27 | 定义定位到 `ch05:230–240`、`ch15:20–28`、`ch18:27–48`、`ch03:129–134`、`app:167–202` |
| 2 | `FKG\|Griffiths\|Ginibre\|Lebowitz\|Chebyshev\|reflection positiv\|positive association` | 全仓 `*.md` | 435（截断 250） | 内部判词：`C6:17,19,31`；`C5:16–17`；`C8:6`；`PROOF-SKELETON:44`；`YANG-MILLS-CONDITIONAL-REDUCTION:551,585`；`BASE-INFINITE-VOL:47–48` |
| 3 | `Test-Path research/v3/species/GKS-GAUGE.md` | 仓库 | **0（文件不存在）** | `Q4-LOGCONVEX.md:46` 引用的 `GKS-GAUGE §4` 是**悬空内部引用**，不得当出处用 |
| 4 | `Lebowitz "GHS and other inequalities" Communications in Mathematical Physics 35 1974 87 DOI` | web | 2 源 | CMP 35(2) 87–92；DOI 两版本冲突（`BF01646608` vs `BF01646068`）⇒ 标 `[E]` |
| 5 | `Ginibre inequality phi^4 fourth cumulant negative Lebowitz inequality proof statement` | web | 3 源 | `U₄≤0` = Lebowitz 不等式（m=2；m>2 无严格结果）；Ginibre Thm 4（Simon 讲义） |
| 6 | `Fortuin Kasteleyn Ginibre "Correlation inequalities on some partially ordered sets"` | web | 1 源 | CMP **22**(2) 89–103 (1971)，摘要原文已核 |
| 7 | `"Ginibre inequality" "lattice gauge" correlation inequality SU(N) Wilson action` | web | 10 源，**非阿贝尔 0** | 仅 **abelian**（Villain/Wilson/Manton）有：Chevyrev–Garban arXiv:2404.09928 |
| 8 | `non-abelian lattice gauge theory FKG inequality Griffiths inequality SU(N) Wilson action correlation inequalities` | web | 10 源，**SU(N)/非阿贝尔定理 0** | 全为 Ising/Z₂/FK 渗流/abelian；**未查到（检索式：见本行；工具：OK；命中：0）** |
| 9 | `Lebowitz inequality lattice gauge theory non-abelian` | web | 10 源，**非阿贝尔规范 0** | 全为 φ⁴/Ising/铁磁自旋；**未查到（工具：OK；命中：0）** |
| 10 | `Osterwalder Seiler 1978 gauge field theories on a lattice reflection positivity Annals of Physics 110 440 DOI` | web | 1 源 | `10.1016/0003-4916(78)90039-8`（书目级） |
| 11 | `Chebyshev integral inequality Hardy Littlewood Polya Theorem 236` | web | 10 源，**精确 Thm 号 0** | HLP 存在且是标准出处；**Thm 236 编号未核到（工具：OK；命中：0）**，故标 `[T-standard]/书目级` |

---

## 6. 与 `H1A-EXACT.md` 的对账（不改其结论）

- `H1A` 用的 `δ_max=G''(0)` 与 `Var−E[H]−W_full=−G''(0)` 本文件**逐字复算一致**（`H1A §1.3` ↔ `ch18:32–37`）。符号表四项（`+0.135995 / −0.193019 / +2.0101 / +0.22196`）本文件直接引用，未重算。
- `H1A` 把力写成 `G`，本文件按 `eq:split` 区分 `G(a)`（自由能）与 `G_v`（力）。**这是唯一的记号冲突，不影响任何数值。**
- `H1A §2.3` 的"`δ<0 ⟺ 真空不在条件化点"在本文件被重述为"`G` 在 0 处非局部 log-凹"，**同一件事的两种说法**（`[T]` 恒等式 + `[C]` 数值）。

## 7. 未闭合（诚实标注）

1. "非阿贝尔 `Hess V⪰0` 失败 ⇒ 无 BL/Ginibre 型上界"是**否定性检索结论**（命中 0），不是不可能性定理；要把它升成 `[NG]` 还需一条"任何上界族都必须过凸性"的分类定理——本文件没有。
2. `W_geo` 在 `eq:split` 里用**条件化 Wilson 期望**书写，而 `W_geo≥0` 的证明（Prékopa）在**纯 Haar** 边缘上；两者在 `Φ'≡0` 的切片（单 plaquette `θ=π`、P4 周期、Dirichlet 弦）上重合，在一般切片上未必。这是**论文级**的记号缝隙，标记待审（不影响本文件的符号判决：`G''` 的符号由 `H1A` 的精确值直接钉死）。
3. 本文件全部 `G''` 均为格点单位曲面下界，未做 `a→0` 标度（`ch03 §sec:dim`）。

**标签**：§1 定义与 `eq:split`/`eq:v3object` `[T]`（逐字核，含行号）；§2 结构等价 `[T]`（恒等式）+ 反例数值 `[C]`；§3 四条判词 `[NG]`（方向论证 `[T]`，文献 `[C]/[E]` 分级已标）；§4(a)(c) `[T]`，§4(b) `[NG]/[C]`；§5 工具状态与命中数为实测。

---

## 勘误（O14 审计）

来源：`O14-AUDIT-GPP.md` `§1.4`、`§4.1–4.2`、`§2`/`§0 判决表 #2`。三处均为**本文件**行文错误；`o14` 的三条结论未错。

### E1 `‖v‖⁻⁴` 缺失（§1.2 恒等式）

- **原（逐字，§1.2 原行）**：`G''(v;a) = E[ ∂_v²V | ξ=a ] − Var( ∂_vV | ξ=a ) = E[H_vv+Φ''|a] − Var(G_v+Φ'|a) ,`
- **正**：`G''(v;a) = ‖v‖⁻⁴·{ E[ ∂_v²V | ξ=a ] − Var( ∂_vV | ξ=a ) }`
- **幂次核准（按论文，非按记忆）**：`ch05-theorems.tex:230–240`（`eq:split`）与 `ch15-uvface.tex:20–28`（`eq:v3object`）**都无显式 `‖v‖⁻⁴`** —— 两式直接写 `G''(v;a₀)=E_a[∂_v²V]−Var_a(∂_vV)`，即**隐式 `‖v‖=1` 约定**；同一章 `ch05:208`（Lemma B `E[Φ'']≥(N/12)‖v‖²`）与 `ch05:214–228`（`eq:wgeo-quant`，"an `L`-link average gives `‖v‖²=1/L`"，`c=(N/12)/L`）证明 `∂_v:=⟨v,∇⟩` **未归一**。故：幂次确为 **`−4`**（不是 `−2`，也不是无因子）；`‖v‖=1` 时退化为论文原式，`‖v‖²=1/L` 时本文件与 `eq:split` 差 `L²`。
- **两条独立核准**：(i) 纤维分解 `dX=‖v‖⁻²da dY` ⇒ `G''=‖v‖⁻⁴(E[∂_v²V]−Var(∂_vV))`（`O14 §1.2 (i)–(iv)`）；(ii) 标度：`ρ_{λv}(a)=λ⁻¹ρ_v(a/λ) ⇒ G''_{λv}(a)=λ⁻²G_v''(a/λ)`，而 `{E[∂_v²V]−Var(∂_vV)}_{λv}=λ²{…}_v` ⇒ 因子 `λ⁻⁴`。
- **未动**：§1.3 的 `eq:split` 逐字引用行、§4 的 `G''≥0 ⟺` 结构等式（后者是 `δ=0` 的符号陈述，独立于归一化）。**`o14` 侧无错**。

### E2 漏 `W_full`（物理点）

- **原**：物理点写作 `H1(δ=a²Δ₀²>0)` / `G''(0) ≥ a²Δ₀² > 0`（§0(c)、§4(a) 引文、§4(c)）
- **正**：`δ := a²Δ₀² + W_full`，`W_full = W_geo + X` ⇒ `G''(0) ≥ a²Δ₀² + W_full`；等价的方差式 `Var(G_v+Φ'|0) ≤ E[H_vv+Φ''|0] − a²Δ₀² − W_full`
- **出处（逐字）**：`docs/paper-v2/app.tex:190–191` —— "with $W_{\mathrm{full}}=W_{\mathrm{geo}}+X$ and $\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality"；`app.tex:796` 符号表 "$W_{\mathrm{full}}$: $W_{\mathrm{geo}}+X$"；`app.tex:182–188`（`app:sharp` 的 (M) 式：`Var(⟨v,∇V⟩|0) ≤ E[H_vv|0] − a²Δ₀²`）。
- **`a²Δ₀²` 是什么（语料核）**：它就是"`G''` 在格点单位下的目标值"——`Δ₀>0` 是 `(L,a)`-一致的**无量纲**下界，`G''_a/a² ≥ Δ₀²` ⟺ `m_lat ≥ aΔ₀`（`ch03-existence.tex:158–163`：`λ_a ∼ a²Δ₀²`，"固定 `a` 下 `G''_a≥λ`（常数）"是量纲错；`ch13-machine.tex:65,70`；`ch18-scalar-h1.tex:90`）。**语料从未给 `Δ₀` 数值**（全仓 `Δ₀|Δ₀²` 命中 41 条，全为符号），故此处无可"数不对"。
- **后果限定**：`W_geo≥0` 已证 `[T]`（`eq:wgeo`），`X` 符号不定（实测 `−0.09…+0.16`）⇒ **`W_full≥0` 未证**；本文件原先把一个待证项当已知用了。判词方向不变（`δ>0` 仍严格强于 `δ=0`）。
- **同错范围（不在本文件，仅记录）**：`ch18-scalar-h1.tex:90`、`H1A-EXACT.md §4.1`、`LOOP-RUNBOOK.md:2158` 同样写 `δ=a²Δ₀²`；本文件**未改**这些文件。

### E3 窗口标签（`−0.08` 的适用窗口）

- **原**：`只到 G''≥−0.08（ABC2-CONSTRUCTIVE.md §4；…）` —— 无窗口标签（沿用 `B3/B4/B5` 的 `β=2.5` 读法）
- **正**：窗口是 **`SU(3)`、`β=5.7`、`L=8`**，且该链条是 **Cauchy–Schwarz（Fréchet–Hoeffding 退化）**，不是"Chebyshev 给出的界"
- **出处**：`H5-HEAD-NUMERICS.md:156`（`k0 SU3 | 5.7 | 8 (12) | … | W_geo=+0.318±0.006`）、`:27`、`:169`（`W_geo=+0.254`(β=2.4) / `+0.318`(β=5.7)，SE≈0.005）；`B4-VARIANCE-MAG.md:15`（"single β=5.7"）、`:6`/`:30`（`VarG≈2`、`VarΦ'≈0.02`、`W_geo≈0.32` ⇒ `|X|≤0.40` ⇒ `G″≳−0.08`）；退化判据 `C6-DEPENDENCE-CLASS.md:29`、`docs/PROOF-SKELETON.md:44`。
- **判词不变**：CS 要闭合需 `VarG<W_geo²/(4VarΦ')≈0.8–1.3<2`，在**任何**窗口都不成立。**`o14` 侧无错**（它就是这么记的）。

### 自检（`pwsh` `Select-String -SimpleMatch`，工具：OK；范围=本文件，共 233 行；"正文"=1–194 行，`## 勘误` 在第 195 行）

| 检索式 | 命中 正文 / 全文 | 说明 |
|---|---|---|
| `‖v‖⁻⁴` | 1 / 6 | 新式已在 §1.2（余 5 处在勘误节内引用） |
| `= E[ ∂_v²V \| ξ=a ]`（旧式首，无因子） | 0 / 1 | 正文旧式已无（1 处为 E1 的"原（逐字）"引文） |
| `a²Δ₀²>0` | 0 / 2 | 正文旧物理点已无（2 处为勘误节引原式） |
| `a²Δ₀²+W_full` | 3 / 4 | §0(c)、§4(a)、§4(c) |
| `− a²Δ₀² − W_full` | 1 / 3 | §4(a) 的方差式 |
| `（CS 退化，非 `β=2.5`）` | 1 / 2 | 新窗口标签已在 §3 行 1（第 2 处是本表自身的引文） |
| `只到 `G''≥−0.08`（`ABC2`（旧无窗口形） | 0 / 0 | 错标签已无 |

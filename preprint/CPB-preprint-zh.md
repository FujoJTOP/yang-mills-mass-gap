# Crossover Preservation and the Weak–Strong Handoff for Four-Dimensional Pure Yang–Mills: A Reduction and a Map of Obstructions

# 四维纯 Yang–Mills 的交叉保持与弱—强交接：一个归约及其障碍图谱

> 预印本草稿 · 2026-11-22 · 面向 Zenodo

---

## Abstract

We reduce the existence and mass-gap problem for four-dimensional pure Yang–Mills to a **weak–strong handoff inequality** `q₀ ≥ q₁`, formulate the two-sided *crossover preservation* statements `CPB-W` and `CPB-S` that would discharge it, and prove the handoff algebra in a **machine-checked** form. We then report, with proofs and falsification tests, that every mechanism currently available for `CPB-W`/`CPB-S` fails at the crossover `q ≈ 1`: Bianchi/sign cancellation is absent (`κ_c = 1` provable in the strong-coupling representation; no additive defect charge in 4D non-abelian gauge theory), Gaussian-weighted class preservation is blocked by a cluster constant (`C_cum ≥ 2` against a required `≤ 1.6`), and the numerical proxies that appeared to suggest otherwise are Gaussian tails of a blocked field that is itself Gaussian.

The outcome is not a proof of the mass gap. It is a precise statement of the obstruction: at `q ≈ 1` the single-block activity must be small, either because the defect action is large or because signs cancel, and no known structure supplies either. We also record the surrounding **conditional theorem**: ten named inputs, nine glue lemmas, nine machine-checked blocks, and four explicit external declarations.

**Keywords:** Yang–Mills, mass gap, lattice gauge theory, renormalization group, crossover, Bianchi identity, Lean 4, numerical falsification.

---

## 摘要

本文将四维纯 Yang–Mills 的存在性与质量隙问题归约为一个**弱—强交接不等式** `q₀ ≥ q₁`，给出两侧的*交叉保持*命题 `CPB-W` 与 `CPB-S`，并以**机器验证**的形式证明交接代数。随后我们给出证明与证伪测试，说明现有可用于 `CPB-W`/`CPB-S` 的每一条机制都在交叉点 `q ≈ 1` 失效：Bianchi/带号抵消不存在（强耦合表象中 `κ_c = 1` 可证；4D 非阿贝尔规范理论没有可加缺陷荷），Gaussian 计权的类保持被团簇常数卡住（`C_cum ≥ 2`，而需要 `≤ 1.6`），而看似支持相反结论的数值代理，其实是**本身为高斯**的阻塞场的高斯尾。

这不是质量隙的证明。它是障碍的精确陈述：`q ≈ 1` 处的单块活动度必须小，要么因为缺陷作用量大，要么因为符号相消，而现有结构两者都不提供。我们也记录了周边的**条件定理**：十个具名输入、九条粘合引理、九个机器验证块、四个显式外部声明。

---

## 1. 引言

四维纯 Yang–Mills 的存在性与质量隙是 Clay 千禧问题之一。近年一个被认真评估的构造（2026 年提交、经独立委员会审查）给出了**子列无限体积极限**与 OS 重构，但没有给出质量隙、唯一性或禁闭。本文的工作起点是：把"质量隙"与"存在性"之间的**连接部分**（即所谓 handoff）写成可检查的不等式，并逐一判定所需输入的状态。

本文**不**证明质量隙。它做三件事：

1. **归约**：把质量隙归约为一个显式的、量纲自洽的交接不等式（§2），并证明其代数核心（机器验证）。
2. **两侧的精确陈述**：把交接所需的两个类保持命题 `CPB-W`（弱侧）与 `CPB-S`（强侧）写成可证伪的形式（§3），并给出证明骨架与承重引理（§4）。
3. **障碍图谱**：证明/判定每一条已知机制为何在交叉点失效（§5–§6），并列出剩余路径（§7）。

最后（§8）我们记录围绕此归约的**条件定理**及其机器验证状态，以便读者判断"若输入成立，则命题成立"这一步是否可信。

---

## 2. 交接不等式

### 2.1 设定与常数

固定紧、连通、单连通单李群 `G`；`h^∨` 为对偶 Coxeter 数，`dim f` 为定义表示维数，`λ_G := h^∨ g²` 为 RG 不变量。物理归一化 `T(f) = 1/2` 下

```
b₀ = 11 h^∨ / (48π²) ,      B := 2b₀ = 11 h^∨ / (24π²) ,
1/g²(μ) = B ln(μ/Λ) + O(g²) ,      C_RG := 24π²/11 .
```

**记号纪律**（本文全程遵守）：测度权重系数是 `γ := β/N_c = 2/g²`，`N_c := dim f`；跑动耦合用 `λ_G = h^∨ g²`。二者仅在 `SU(N)` 时重合（`h^∨ = dim f = N`）。混用会在 `g² = 1` 处使 `q₀` 差 `e^{3.59} ≈ 36` 倍。

### 2.2 定理 1（交接不等式的等价形式）

令 `q(μ) := m/μ`（`m` 为物理质量隙，RG 不变），`μ_up` 由 `g²(μ_up) = g*²` 定义，`q₀ := q(μ_up)`，`q₁` 为强侧阈值。则以下等价：

```
(i)   q₁ ≤ q₀ ;
(ii)  ln(m/Λ) ≥ t_up(λ*) + ln q₁ ;
(iii) M := ln(m/Λ) − t_up(λ*) − ln q₁ ≥ 0 ;
(iv)  t_λ ≤ t_S   （(t, ln q) 平面中两域相接）。
```

其中 `t_up(λ*) = 24π²/(11λ*)`（一圈）。**几何读法**：弱侧覆盖 `t ≥ t_λ`，强侧覆盖 `t ≤ t_S`；覆盖全尺度当且仅当 `t_λ ≤ t_S`。这是纯实分析，已在 Lean 4 中形式化（`handoff_iff_margin`、`required_lambda`，公理为 `propext, Classical.choice, Quot.sound`）。

### 2.3 定理 2（射程 ⇒ 显式质量下界）

设弱侧类保持对 `λ ≤ λ_c` 成立（记射程 `W`），强侧 signed-gas 对 `q ≥ q_c` 成立（记射程 `S`）。则交接成立当且仅当

```
m/Λ ≥ q_c · exp(t_up(λ_c)) .
```

**这是本文的核心归约**：质量隙的证明被拆成"两个方法射程的估计"与"物理质量 `m/Λ`"两个可独立攻击的量。§3 的两条命题正是 `W` 与 `S`；§4–§6 说明它们为何在 `q ≈ 1` 处难以获得。

### 2.4 两圈闭式与数表

两圈截断流的精确解：

```
t_up^{2L}(λ) = C_RG/λ − (51/121) ln(1 + 121 C_RG/(51 λ)) .
```

`SU(3)`、`λ = 6`（`g² = 2`）：`t_up^{2L} = 2.639`；取 `m/Λ = 8` 得 `q₀ = 0.571`；`q₀ = 1` 在 `λ ≈ 7.28`（`g² ≈ 2.43`）。一圈给 `3.589`（`q₀ = 0.221`），三圈/数值给约 `2.36/2.0`（方案依赖，**不作严格输入**）。

**结论：`CPB-W` 必须给出至少两圈的严格界**；只有一圈控制时 `q₁` 至多 `0.22`，交接无法闭合。

---

## 3. 两侧类保持

### 3.1 `CPB-W`（弱侧）

**猜想.** 存在 `λ_c ≈ 7.3`（两圈；一圈保守版 `≈ 12`）与弱侧类 `C_W(ε₀)`（以质量 `m` 的 Gaussian 为参考 + 有界相互作用，相对误差 `ε₀ ≤ 1/2`），使有效作用量在全部 `λ ≤ λ_c` 上保持在 `C_W` 中，且一步分块后 `C_W` 不变；同时给出 `t_up(λ_c)` 的两圈严格界 `RB2`。

`RB2` 的精确内容是

```
|β(λ) − β_{≤2}(λ)| ≤ C₃ λ⁴ ,      0 < λ ≤ 7.3 ,
```

`C₃` 显式、对 UV 截断与体积一致；或单边符号 `β ≤ −Aλ² − Bλ³`。

**状态：开放。** 定量的核心 `NPB` 是非微扰 β 界；Balaban 型严格 RG 按构造止于 `λ ~ 1`（见 §4.3）。

### 3.2 `CPB-S`（强侧）

把有效作用量写成 `S_k = S_k^G + V_k`（massive Gaussian + 相互作用），在 signed-gas（polymer）基下定义**交叉活动度**

```
z_k(q) := sup_x Σ_{C∋x} |w_q(C)| · e^{+κ q |C|} ,     κ = 1/2 ,
```

其中 `w_q(C)` 是连通团簇 `C` 的带号权重，`q = m a_k` 是块单位的无量纲质量。**猜想.** 存在绝对常数 `q_c ∈ [1, 3/2]`、`z_c = 1/2`，使 `q_k ≥ q_c` 时 `z_k(q_k) ≤ z_c`，且类在一步分块后保持。

**关键结构**：`q_{k+1} = 2 q_k` 是精确运动学，且 `z` 随 `q` 单调不增（在 `c₀ > κ` 时），故**最坏点在窗口最小 `q`，即 `q_c`**；因此 `CPB-S` 只需在 `q_c` 处证明一次。

### 3.3 为什么标准路线在交叉点失效

- **Balaban**：小场展开围绕平凡联络，展开参数是耦合，按构造在 `λ ~ 1` 停止；`λ ~ 6` 处大场不被 `β` 压制。
- **强耦合展开**：收敛半径由 `β/(2 dim f)` 控制；`q₁ ~ 1` 需 `β ~ 2.2`，而严格 cluster expansion 的阈值远低于此。
- **Bakry–Émery**：Wilson 作用的 Hessian 含 `−[F,F]` 型负项，只有质量项能压住它，而质量未知。

---

## 4. 证明骨架与承重墙

### 4.1 九引理骨架

`CPB-W` 的类保持部分可分解为十条引理（`[B]` = Balaban/标准，`[新]`，`[未解]`）：

```
1  分块与协方差分解         [B]        G = B^T G' B + G_⊥
2  Gaussian 条件期望/Wick   标准
3  Cumulant 展开            [B]+[新]    余项 O(ε³L^{3/2})，需 εL^{1/2} ≤ 1/2
4  涨落和（一圈 log）        [新]       L(q) ≤ C₃ ln(1/q) + C₄
5  Lipschitz 收缩           [B]+[新]    sup 不收缩，C¹/Lip 才看得到
6  大场 Lipschitz/C² 界      [B] 仅 g² 小；q~1 [未解]   ← 承重墙
7  规范无关/等变分块         [B]+[新]
8  体积/UV 一致             [B] 部分
9  一步类保持               [新]       U(1) 条件性成立
10 迭代 + RB2/NPB           [未解]
```

### 4.2 引理 6：大场界在 `q ≈ 1` 失败

以 `C²` 范数（sup + Lip + Hessian）在球 `‖A‖ ≤ R` 上控制 `V`，非退化要求

```
S'' ≥ (c₀ m² − Cε)·Id > 0 ,     ε < c₀ m² / C .
```

非阿贝尔交换子给出 `‖V''‖ ≲ R³/g² + R²/g`；在 `q ≈ 1`、`g² ≈ 2`、典型场幅 `R ~ g/q ≈ 1.4` 处这是 `≈ 2.8`，**超过 `c₀ m²/C` 的容许值**。Balaban 的大场估计属于同类，但只对小的 `g`（UV）证明；没有 `q ~ 1` 的版本。

**这是 `CPB-W` 的类保持核心，不是 `RB2/NPB`。**

### 4.3 Gaussian 计权类与团簇常数

一个自然的补救是改用局部 `L²` 类

```
C_W^G(ε) = { S = S^G + V : sup_x ∫ |V_x|² dμ_G ≤ ε² } ,
```

`μ_G` 为协方差 `g²(−Δ+q²)^{-1}` 的 massive Gaussian。它把展开参数从 `ε ≈ 2`（必败）移到

```
εL^{1/2} = 0.24 + C_cum·R₃ ,      R₃ = 0.16–0.46 ,
```

过线需 `εL^{1/2} ≤ 1/2`，即 `C_cum ≲ 0.6–1.6`。而标准的 Kotecký–Preiss 界在 `g² ~ 2` 处给 **`C_cum ≥ 2`**：

```
C_cum = 1 :  0.40–0.70  （仅低端过）
C_cum = 2 :  0.56–1.16  （不过）
C_cum = 3 :  0.72–1.62  （不过）
```

**所以 Gaussian 计权不绕过引理 6，它把墙移到 `C_cum`，而标准常数已经出局。**

---

## 5. Bianchi 机制及其失败

一个反复出现的希望是：Bianchi 恒等式 `dF = 0` 会通过符号抵消把活动度压下去。我们逐一判定。

### 5.1 BL-1：`dF=0` 是拓扑陈述

`dF = 0` 精确等价于「通量 2-余循环」（每个立方体上 `∏_{p∈∂c} U_p = 1`）。在阿贝尔 + 量子化通量的情形它才推出「缺陷 2-链闭合」。它约束通量守恒，**不约束单块权重的符号或大小**。

### 5.2 BL-3：闭 2-链可以是 codim-1

`x₄ = 0` 切片是闭 2-链，但它是 **codim-1 平板**，其逐块代价与 `q` 无关。因此「缺陷是闭 2-曲面」**不蕴含**「面积型作用量 `σ(q) ~ q²`」。

### 5.3 BL-4：改进只在熵，不在作用量密度

`dF = 0` 作用在连通 2-链的**计数率** `γ_all → γ_cyc` 上，一个**与 `q` 无关的纯数**；它不作用在 `|w_q(C)|` 上。在加法形式 `|w| ≤ exp(−c₀q|C| − c₁q²|C|)` 中，**`c₁ = 0` 可证**。

### 5.4 BL-5：`q ≈ 1` 处二次标度是空的

在 `q = 1`，任何 `q^α`（`α > 1`）恒等于 `q`。**任何"二次标度"的改进，在需要它的那个点上增益恰为零。** 真正的门槛是 `q = 1` 处的绝对效率 `σ(1)`，不是标度。

### 5.5 BL-6：强耦合下 `κ_c = 1`

唯一能同时服务两侧的机制是**固定边界上的带号抵消**

```
|Σ_{∂C=∂C₀} w_q(C)| ≤ e^{−c₀q|C₀|} · κ_c^{|C₀|} ,      κ_c < 1 ,
```

其中 `κ_c` 是每块抵消因子（`κ_c = e^{−δ}`）。`CPB-W` 需 `κ_c ≤ 0.80`；`CPB-S` 需更强。

**判定（BL-6，非定理）：默认 `κ_c ≈ 1`，无已知机制给出 `κ_c < 1`。** 三条方向一致的证据：

1. **强耦合 surface/character 表象中 `κ_c = 1` 可精确证**：`SU(2)` 系数 `2I_{2r+1}(β)/β > 0` 全为正，曲面和各项非负，带号和等于绝对和；
2. **权重局部因子化** ⟹ 符号不随面积反相；`κ_c < 1` 需要逐块的 `e^{iπ·可加荷}`；
3. **`dF = 0` 不提供可加缺陷荷**：2D `O(N)` 靠精确约束 + 拓扑指标，4D `U(1)` 靠对偶/Poisson，**4D 非阿贝尔纯 YM 两者皆无**（瞬子荷密度不是缺陷荷，`θ` 只做整体移相）。

`U(1)` 给出一个干净的反证：`F = dA` 使 `dF = 0` 自动成立，却没有非平凡的带号相消。**故 `dF = 0` 与 `κ_c < 1` 逻辑独立。**

> **这是关于路线的负结果，不是对 `CPB-W`/`CPB-S` 为假的证明。** 被否定的是"用 Bianchi / 符号抵消证明它们"。

---

## 6. 数值检验

### 6.1 装置与验证

`SU(3)` Wilson 作用量，Cabibbo–Marinari 热浴（精确 `SU(2)` 子群采样器），与 Metropolis、强耦合级数、已知 `β = 6` plaquette 交叉验证。格点 `L = 4, 6, 8, 10`；主网格 `β = 2.6, 2.8, 3.0, 3.3`（`g² = 1.8–2.3`），各 8000 测量。

### 6.2 `R₄` 与 `κ_eff`

`R₄`（阻塞 plaquette 的连通四点比）在所有 `β, L, b` 上 `|R₄| ≲ 0.13`（`b ≥ 2` 时 `≲ 0.01`），无液滴结构。

`κ_eff`（投影 Hessian 最低特征值）在 `β = 1.7…6` 全部为负，**且在 `β = 40` 的弱耦合对照处仍为 `−1.74`**。其标度 `|κ_eff| ≈ (0.85–1.6)·β(1−⟨P⟩)` 把它识别为**微观** Wilson 作用的 Savvidy/Nielsen–Olesen 负模——任何非凸作用量都有，包括连续统 YM。**因此它不是 `CPB-W` 的判别量**，本文报告但不用。

### 6.3 特征标系数全正

plaquette 分布的字符系数 `c_R = ⟨χ_R(U_p)⟩`：

| `β` | `c₃` | `c₆` | `c₈` | `c₁₀` | `r₁` |
|---|---|---|---|---|---|
| 2.6 | +0.52107(19) | +0.13495(18) | +0.23809(26) | +0.02323(16) | 1.0000 |
| 3.0 | +0.61476(19) | +0.18739(20) | +0.32620(27) | +0.03800(17) | 1.0000 |
| 3.3 | +0.68861(19) | +0.23549(20) | +0.40444(29) | +0.05317(18) | 1.0000 |

`L = 6` 与 `L = 8` 一致。**没有 `c_R < 0`**，单块带号/绝对比 `r₁ = 1.0000`：`c₁` 层面没有符号抵消。这与 BL-6 的 `κ_c ≈ 1` 一致（是**确认**，不是失败），但它只是 `c₁` 代理，不是 `κ_c` 的测量。

### 6.4 缺陷比例是高斯尾

一个看似强大的代理是用阈值缺陷比例估计缺陷作用量 `c₀ = −ln(fraction)`。该代理**被证伪**：阻塞场本身是**高斯**的，任何阈值比例只是高斯尾。

```
β=2.6, b=2, L=8:  μ=0.17355, σ=0.02805
高斯尾预测:  f(⟨P⟩/2)=0.00099 ,  f(0.75⟨P⟩)=0.0610
实测:        f(⟨P⟩/2)=0.00089 ,  f(0.75⟨P⟩)=0.0602
```

`L=6` 与 `L=8` 的方差一致到 `<1%`（`0.000781` vs `0.000787`），故非有限体积效应。换阈值与耦合，`c₀` 在 `2.4–10.1` 之间变动——**它由阈值与分布宽度决定，不是独立作用量**。

### 6.5 `ξ/a ≈ 7` 的可行性

σ 证书要求 `q₁^σ ≤ 0.143`，对应 `ξ/a ≈ 7`。实测 `ξ/a(β=6) = 0.93(8)`；`β ≤ 3.3` 处 `ξ/a ≲ 0.3`（关联子在一个格距内消失）。`ξ/a = 7` 需要 `L ≈ 28–32`。

**数值不构成证书。** `q₁^σ ≤ 0.143` 是一条严格的数学不等式（定理或输入），点值或置信区间都不能证明它；数值只能提供可行性/一致性证据。若作证据，正确的统计对象是 **2σ 下界**：5% 精度下需测到 `ξ/a ≈ 8.5`（`8.5 − 2×0.425 = 7.65 > 7`），且需 `L/ξ ≥ 4`，即 `L ≳ 34`。`L = 12` 在 `ξ ≈ 7` 时 `L/ξ = 1.7`，只能定出 `ξ/a ≈ 3–4`，**不能提供 `ξ/a ≥ 7` 的证据**。

纯 NumPy 单核成本 `∝ L⁴`：`L=32` 约 `38 s/测量`，`10⁴` 测量约 `105 h`，`10⁵` 测量约 `1050 h`。换 kernel（C/OpenMP，16 核）与各向异性体积（`16³×64` 或 `24³×64`）后，按 `10⁶` 量级独立组态估计仍需**数天**。**这是可行性限制，不是数值噪声。**

---

## 7. 剩余路径

以本文的判定，`CPB-W`/`CPB-S` 在现有机制下无已知路线。剩余候选：

1. **逆蒙特卡洛（IMC）**：拟合阻塞有效作用量的耦合 → 得到 `w_q(C)` → 直接给 `c₀` 与 `κ_c`。这是**唯一剩下的数值途径**，代价是一个独立项目。
2. **新的可加符号/指标结构**：`κ_c < 1` 需要的对象。2D `O(N)` 与 4D `U(1)` 各有其版本，4D 非阿贝尔是否可能存在，目前未知。
3. **改变交接架构**：不同的弱/强匹配、不同的 RG 方案，或放弃交接改走全耦合构造。
4. **换 certificate**：用弦张量 `σ` 而非质量 `m`，或换信号更干净的可观测量；`ξ/a ≈ 7` 这道门槛是否可绕，未知。

---

## 8. 条件定理与机器验证

本文的归约嵌套在一个更大的**条件定理**中：十个具名输入 `IN-0 … IN-9`（含 2026 构造的认证输入、UV 极限、`O(4)`、`CPB-W`、`CPB-S`、重叠、IR 唯一性、非平凡性、群一致、谱预解收敛）⇒ Clay 陈述。

其"九变一"部分是九条粘合引理 `GL1–GL9`，其中：

```
无条件 ✅（6 块）：GL1（审计）· GL2a · GL4 · GL6 · GL5-KIN · GL8-A
✅-蕴含（6 块）   ：GL3 · GL2b′ · GL5-DYN · GL8-S · GL9 · GL7
无可形式化内容    ：GL5-DYN（IN-2 ⟹ OS2 是记账）
审计（非定理）    ：GL1（跨文档合规）
```

**机器验证状态**：11 块中 9 块已 Lean 4 化，无 `sorry`，公理冻结为 `[propext, Classical.choice, Quot.sound]`。四件"显式声明而非机器证明"的外部件被推到**类型层面**：

1. Minlos（GL3 的 `hMinlos` 假设，Mathlib 缺失）；
2. 截断的 Borel 可测性（GL4 的显式假设，Mathlib 缺 `comp_stronglyMeasurable`）；
3. 变分隙 ↔ 谱隙的等同（GL2b′ 的声明性 scope）；
4. Reed–Simon VIII.24b（GL2b′，Mathlib 缺失）。

**这四件是外部假设，不是证明。** 它们被写进签名或文件头，正是为了让读者能看见——而不是被散文藏起来。

---

## 9. 结论

本文把四维纯 Yang–Mills 的质量隙问题精确地归约为一个交接不等式，并系统判定：每一条已知的、能在 `q ≈ 1` 处压低单块活动度的机制，都在需要它的那个点上失效。

- Bianchi 是拓扑的，作用于熵而非权重；`q²` 缩放在 `q = 1` 处为空；强耦合下 `κ_c = 1` 可证。
- Gaussian 计权把墙移到团簇常数 `C_cum`，标准估计已经出局。
- 数值上，`κ_c` 不可直接测，σ 的缺陷比例代理是高斯尾，`ξ/a ≈ 7` 是可行性边界。

**这不是失败的报告，而是一张地图**：它告诉后来者哪些路已经走过、每条为什么封死、以及在什么条件下会重新打开（`q ≈ 1` 处实测 `κ_c < 1`，或出现新的可加符号/指标结构）。

---

## 附录 A：常数与符号

| 符号 | 含义 | 不可混用 |
|---|---|---|
| `λ_G` | `h^∨ g²`（RT 不变量） | 不是 `dim f·g²` |
| `b₀` | `11h^∨/(48π²)` | 不是 `11 dim f/(48π²)`（仅 SU(N) 重合） |
| `B` | `2b₀ = 11h^∨/(24π²)` | 不是 `b₀` |
| `γ` | `β/N_c = 2/g²`（测度权重） | 与 `λ_G` 是两套变量 |
| `κ` | `1/2`（Kotecký–Preiss 因子） | 不是 `κ_c` |
| `κ_c` | 带号抵消因子 `<1` | 不是 `κ` |
| `σ(q)` | 缺陷作用量密度 | 不是弦张量 `σ` |

## 附录 B：复现

```
bash lean/check-handoff.sh            # 交接代数，公理冻结
bash lean/check-conditional.sh        # 条件定理骨架
bash research/glue-a-limits/check-glue-a.sh
bash research/glue-b-blocking/check.sh
bash research/glue-c-gauge-group/check-glue-c.sh
bash research/glue-d-dual-handoff/check-gl7.sh
bash research/os2-rotation/check-symdim4.sh
python research/os2-rotation/check-d14.py
```

数值：`research/lattice-cpb/`（热浴、分析、`t_up.csv`、`block_summary.txt`）。

## 参考文献（部分待核）

1. K. Wilson, *Confinement of quarks*, Phys. Rev. D 10 (1974).
2. K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. 110 (1978).
3. T. Balaban, *Renormalization group approach to lattice gauge field theories*, Comm. Math. Phys. 102 (1985), 109 (1986), 122 (1989).
4. R. Kotecký, D. Preiss, *Cluster expansion for abstract polymer models*, Comm. Math. Phys. 103 (1986).
5. S. Cao, M. Park, S. Sheffield, *Random surfaces and lattice Yang–Mills*（2024，待核）.
6. A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Stochastic quantisation of Yang–Mills–Higgs in 3D*, Invent. Math. 237 (2024).
7. 2026 年提交的 4D 纯 YM 构造（技术报告，经独立委员会评估；未发表）。
8. hexagonmath 重构 `hexagon:2611.00492v1`（站点 ID，**不是 arXiv 编号**；AI 辅助预印本）。
9. 2D `O(N)` 构造性结果（具体出处待核）。

---

> **免责声明.** 本文不声称证明 Yang–Mills 存在性或质量隙。所有"开放"的输入在正文中逐条标注；所有负结果均附证伪测试。若读者发现某条判定有误，请以正文中的可复现命令与数值数据为准。

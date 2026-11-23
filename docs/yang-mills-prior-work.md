# 杨–米尔斯存在性与质量间隙：前人方向与当前坐标

> 调研日期：2026-11-22
> 目的：在用 Lean 4 介入之前，先把问题的精确陈述、已有路线、2026 年的新进展和形式化现状对齐。

---

## 0. 先给结论

1. **问题仍未解决，但坐标已经变了。** 2026 年 9–11 月发生了三件事：
   - OpenAI 用 AI 生成了 2D O(N) sigma model 的五篇预印本，其中"指数衰减"一篇同步做了 Lean 形式化；
   - OpenAI 提交了 4D 纯 Yang–Mills **存在性**的 Lean 构造（9,470 个模块，Lean 4.32），2026-11-10 由 David Gross 主持的五人委员会出具验证报告；
   - 一份独立重构（hexagonmath, `hexagon:2611.00492v1`，站内 ID，非 arXiv）把构造拆成 16 个已认证输入，并证明**这 16 个输入推不出质量间隙**；质量间隙被归约为一条"体积一致、标号一致的指数聚类速率"（记作 MG）。

2. **因此正确的介入方式不是"从零证明 YM"，而是站在现有 Lean 工件上做质量间隙。** 存在性部分已经有可编译、可审计的代码；真正的空白是：
   - 质量间隙的**非空洞陈述层**；
   - 强耦合 → 弱耦合的 **RG 交接**；
   - 无限体积极限的**唯一性**；
   - 4D 中无穷多个 counterterm 的控制。

3. **最现实的第一阶段成果**：把"质量间隙"形式化为一个可机器检验的不等式，证明它与谱隙、指数聚类、转移矩阵间隙等价，并形式化已知的强耦合质量间隙。这不是千禧年奖，但是后续任何证明都必须踩的台阶。

4. **关于"黎曼猜想已被 AI 证明"**：截至 2026-11-22，这是 Levent Alpöge 于 11-20 发布的**声称**，附带约 1700 万行 Lean 代码，尚未通过同行评审，也尚未完成独立的全量检查。Lean 编译通过不等于陈述写对了，这一点对杨–米尔斯同样适用。

---

## 1. 问题的精确陈述

### 1.1 Clay 官方表述

Jaffe–Witten 的原文（Clay 官网）：

> Prove that for any compact simple gauge group G, a non-trivial quantum Yang–Mills theory exists on ℝ⁴ and has a mass gap Δ > 0. Existence includes establishing axiomatic properties at least as strong as those cited in Streater–Wightman (1964) and Osterwalder–Schrader (1973, 1975).

拆开是两件事：

| 子问题 | 内容 |
|---|---|
| **存在性** | 对任意紧单群 G，构造一个非平凡的、满足 Wightman 或 OS 公理的 4D 量子场论 |
| **质量间隙** | 该理论的哈密顿量 H 在真空之外有正的下界：`spec(H) ∩ (0, Δ) = ∅`，Δ > 0 |

"非平凡"排除了自由场。对 SU(3)，质量间隙意味着胶球质量有正下界。

### 1.2 Wightman 公理（Minkowski 侧）

- **W0**：可分 Hilbert 空间 ℋ 上有幺正的 Poincaré 表示；谱条件 `spec(P) ⊂ V₊`（`P₀ ≥ 0`, `P₀² − P² ≥ 0`）；存在唯一的 Poincaré 不变真空 Ω。
- **W1**：场是算子值分布，`φ(f)` 与测试函数做 smear，定义在公共稠密域 D 上并连续。
- **W2**：协变性 `U(a,Λ) φ(f) U(a,Λ)⁻¹ = φ(f_{a,Λ})`。
- **W3**：微观因果性 —— 类空间隔的场对易（或反对易）。
- **W4**（常并入）：真空循环性（Reeh–Schlieder）。

### 1.3 Osterwalder–Schrader 公理（Euclidean 侧）

设 μ 是 `𝒮′(ℝ⁴)` 上的概率测度，Schwinger 函数 `Sₙ(x₁,…,xₙ) = ∫ φ(x₁)⋯φ(xₙ) dμ`。Douglas 等人在 Lean 中采用的是 Glimm–Jaffe 形式：

| 公理 | 内容 |
|---|---|
| **OS0（解析性）** | `(z₁,…,zₙ) ↦ S(Σ zⱼ fⱼ)` 在 ℂᴺ 上整 |
| **OS1（正则性/温性）** | `|S(f)| ≤ exp[c(‖f‖₁ + ‖f‖ₚᵖ)]`，`1 ≤ p ≤ 2` |
| **OS2（Euclidean 不变性）** | 对每个 Euclidean 等距 E，`S(f) = S(Ef)`，即 `μ = E∗μ` |
| **OS3（反射正性）** | `⟨θA, A⟩_{L²(μ)} ≥ 0`，A 属于正半空间上场的指数所生成的类 |
| **OS4（遍历性/聚类）** | 时间平均收敛到均值；等价地 `S(f + T(s)g) → S(f)S(g)` |

**OS 重构定理**：OS0–OS4 ⇒ 存在 Wightman 理论。这是把 Euclidean 概率测度翻译成 Minkowski 量子场论的标准桥梁。

### 1.4 质量间隙的等价表述

这是形式化时最该先钉死的部分。以下三条等价（在 OS + 重构的前提下）：

| 表述 | 内容 | 形式化难度 |
|---|---|---|
| **谱隙** | `H ≥ Δ > 0` 在 `(ℂΩ)⊥` 上 | 中：需谱理论、GNS |
| **指数聚类** | 对局部算子 O，`|⟨O(x)O(y)⟩_c| ≤ C e^{-m|x−y|}` | 中：需连通关联函数 |
| **转移矩阵间隙** | `λ₂/λ₁ ≤ e^{-aΔ}`，a 为格距 | 低–中：有限维矩阵，格点友好 |

**关键警告（hexagonmath `hexagon:2611.00492v1`）**：仅有 OS4 的"聚类"（不带速率）**推不出**谱隙。他们给出一个满足全部 16 个输入的显式 OS 模型，其 H 没有谱隙。因此质量间隙必须表述为**带一致速率的指数聚类**，而不能只写"聚类"。

这也是为什么 Google DeepMind `formal-conjectures` 的 Yang–Mills PR 被审稿人指出"字段是 `True`、定义空洞"：陈述层一旦写松，形式化就失去意义。

---

## 2. 七条主要路线

### 2.1 强耦合格点展开

**核心思想**：在格点上取 β 小（耦合大），把 `e^{β tr U(plaquette)}` 展开成 Taylor 级数，用群积分的正交关系逐项计算 Wilson 环期望。展开项组织成二维"弦世界面"，收敛半径正。

**已证**：
- Wilson 面积律 `|E[W(L_A)]| ≤ C e^{-σA}`；
- 质量间隙（以 `‖x−y‖₁` 表述）；
- Osterwalder–Seiler (1978)、Seiler 专著 (1982) 的公理化；
- Cao–Park–Sheffield, *Random surfaces and lattice Yang–Mills* (arXiv:2307.06790, Comm. AMS 2025) 给出任意亏格嵌入曲面的求和表示；
- Cao–Nissim–Sheffield 把面积律推到更大参数区（arXiv:2505.16585, 2509.04688）。

**卡点**：强耦合是 β 小。连续极限要求 β → ∞（渐近自由：`β = 2N/g²`, `g → 0`）。**强耦合展开到不了连续极限。**

### 2.2 Balaban 重整化群

**核心思想**：在格点上做分块自旋（block spin），逐尺度积分掉高频自由度，保持有效作用量在一个可控的类里。小场区域做鞍点/Gaussian 展开，大场区域证明其稀少。

**已证**（1980 年代 CMP 系列）：
- 4D 格点 YM 的 **UV 稳定性**：配分函数在 `a → 0` 时不爆炸，沿渐近自由轨迹重整化参数；
- 小盒子内有效作用量的控制。

**未证**：
- 规范不变的局部关联函数族；
- 无限体积极限；
- 完整 OS 公理；
- 质量间隙。

**评价**：Douglas 2026 说 Balaban 是"UV 构造的原理性证明"，Dimock 认为"有缺口但可补"。2026 年 OpenAI 的 sigma model 工作正是把 Balaban 的分块方案**改到能跑通**：不围绕约束能量极小点展开，而围绕分片常数的 block spin 展开；粗键用带号气体（signed gas）而非 polymer expansion；有效密度保持在一个固定类里。

### 2.3 Magnen–Rivasseau–Sénéor 与构造性 QFT

**核心思想**：用多尺度分析 / 相空间展开直接构造连续统理论。

**已证**：
- Magnen–Rivasseau–Sénéor (1993)：带红外截断的 YM₄ 构造（规范固定，SU(2)）；
- Glimm–Jaffe–Spencer (1974)：P(φ)₂ 的 Wightman 公理与粒子结构；
- Feldman–Osterwalder (1976)：弱耦合 (φ⁴)₃ 的 Wightman 公理与质量间隙；
- Aizenman 等：4D φ⁴ 的平凡性。

**卡点**：
- MRG 是规范固定的，不是规范不变的物理理论；
- 4D 非阿贝尔规范理论的无穷多 counterterm；
- 没有质量间隙。

### 2.4 随机量子化

**核心思想**（Parisi–Wu）：把 QFT 测度定义为 Langevin 动力学的唯一不变测度：

```
∂φ/∂t = −δS/δφ + ξ,    ξ = 白噪声
```

不动点就是 `e^{−S}` 路径积分。发散通过逐项减除 counterterm 处理。

**已证**：
- 2D YM：Chandra–Chevyrev–Hairer–Shen, *Langevin dynamic for the 2D Yang–Mills measure*, Publ. IHÉS (2022)；
- 3D YM–Higgs：Chandra–Chevyrev–Hairer–Shen, Invent. Math. 237 (2024) 541–696，用 regularity structures 构造重整化流，证明 gauge-covariant in law；
- 任意维格点强耦合：Shen–Zhu–Zhu, CMP 400 (2023) 805–851，用 Langevin 动力学的指数遍历性给出唯一不变测度、无限体积极限、**正质量间隙**；
- 4D：**未解决**。4D 不是超可重整的，需要无穷多个 counterterm，regularity structures 不能直接用。

**评价**：这是概率学派目前最有活力的路线。Douglas 认为它相对传统 CQFT 的最大优势是"不需要复杂的 cluster expansion 就能处理大场问题"。但 4D 仍是硬墙。

### 2.5 二维精确解

**核心思想**：2D YM 可精确求解。

**已证**：
- Migdal (1975)、Gross–Witten (1980) 的精确公式；
- Osterwalder–Seiler 的公理化；
- Driver (1989)：连续统期望、格点收敛、lasso；
- Lévy (2003)、Sengupta (1997)：紧曲面上的 YM 测度。

**评价**：2D 是**形式化的理想第一个目标**：测度显式、OS 公理可证、面积律精确。但 2D 没有传播自由度，质量间隙是退化的，不能直接推广到 4D。

### 2.6 转移矩阵 / 哈密顿量路线

**核心思想**：用 Kogut–Susskind 哈密顿量或 Euclidean 转移矩阵 `T = e^{−aH}`，把质量间隙变成 `λ₂/λ₁` 的界。

**已证**：强耦合下转移矩阵有隙。

**卡点**：需要**对格距 a 一致**的界，且要活到弱耦合。这等价于 RG 路线。

### 2.7 其他路线

| 路线 | 状态 | 备注 |
|---|---|---|
| Seiberg–Witten / monopole condensation | 软破缺 N=2 已解，纯 YM 未严格化 | Douglas 称"迄今最好的解析处理" |
| 张量网络 RG | 数值强，无严格定理 | 凝聚态方法 |
| 大 N / 全息 | 无严格构造 | 启发式 |
| 可积性 | 仅低维特殊模型 | 不适用于 4D |
| 算子代数 / AQFT | 公理化框架成熟，构造缺 | 可作为陈述层参考 |

---

## 3. 2026 年的转折点

### 3.1 OpenAI 的 2D O(N) sigma model 五篇

2026-09-23 至 10-05，OpenAI 从一个未发布模型生成 722 篇手稿（后撤回 3 篇），其中 family 215 是五篇 sigma model 论文：

1. **Exponential decay in 2D classical O(n) models**（9-23）：对 n ≥ 3、任意有限温度，证明两点关联指数衰减；**唯一做了 Lean 形式化的一篇**（约 1,500 行，117 个 Mathlib 文件，标准公理）。
2. **Sharp mass bounds for the 2D O(4) model**（9-23）。
3. **The canonical massive continuum limit of the 2D O(3) model**（10-04，122 页）：构造交互连续极限，满足 OS 公理，唯一真空，正哈密顿量间隙。
4. **An isolated particle pole for the 2D O(3) spin field**（10-04）。
5. **Exact mass asymptotics for the 2D O(4) lattice model**（10-05）。

**注意**：Lean 只覆盖第 1 篇的有限方格衰减定理；无限体积极限、O(4) 谱隙、连续极限、OS 公理都**没有**形式化。社区重构在 `github.com/bauerscr/nlsm`。

### 3.2 模板：可复用的四步

Douglas 把 sigma model 的证明提炼成一个对 YM 也适用的模板：

1. **衰减定理**（格点、所有耦合）：证明关联指数衰减。这是"close the gap"的粗估计。
2. **分块 RG**：Balaban 式，但围绕分片常数 block spin 展开；粗键用 signed gas；有效密度保持固定类；步进定理复现一圈、两圈 running coupling，并带误差界。
3. **间隙检验**：定义无量纲量
   ```
   Δ(n) = 4 log Z(n×n) − log Z(2n×2n)
   ```
   它在精确分块下不变；反射正性说明，若它在某个盒子上小，则转移矩阵有间隙 `≥ log 2 / n`。
4. **一致性**：RG 说明 Δ(n) 在固定物理尺寸的盒子上与截断无关；一次裸尺度估计（Ward 恒等式迭代）说明它通过一次。于是每个截断都有间隙。

### 3.3 OpenAI 的 4D 存在性 + Gross 委员会

- 2026-10-16：OpenAI 完成 Lean 形式化；
- 2026-10-28：提交，代码在 `yang-mills-group/OpenAI-YM-20261028`；
- 2026-11-09/10：Gross（主席）、Douglas、Jaffe、Klebanov、Nekrasov 五人委员会发布 35 页验证报告。

**委员会认可的**：
- 4D 纯 YM 可由 Wilson 格点规范理论的标度极限构造；
- 满足 Wightman 公理；
- 对任意紧单群成立；
- 9,470 个模块重建无错；17 个主定理只用 `propext`, `Classical.choice`, `Quot.sound`，无 `sorryAx`；
- 非平凡性、守恒应力张量、OPE 界、一圈 running coupling。

**明确未认可的**：
- **质量间隙**；
- 无限体积极限的**唯一性**（只有子列极限）；
- 不同子列理论的一致性；
- 收敛而非仅子列收敛；
- 禁闭/面积律；
- 与微扰论逐阶一致。

**委员会自己的话**："We do not claim to fully understand the proof yet."（我们尚未声称完全理解该证明。）报告还说明他们用 OpenAI 的 Astra 模型辅助消化 Lean 代码。

### 3.4 16 个输入与 MG 归约

hexagonmath 的独立重构（`hexagon:2611.00492v1`, 2026-11-15；站内 ID，非 arXiv）把构造拆成：

```
RG1–RG5, UV1, PH1–PH3, AN1, Y1–Y4, B1, C1   （共 16 个已认证输入）
```

在这 16 个输入下，它复现了：固定体积 UV 极限、子列无限体积极限、相分解与聚类、OS 重构、应力张量与 Ward 恒等式、迹反常、一圈 running coupling、OPE。

它的两个逻辑边界最有价值：

1. **反例**：存在一个满足全部 16 个输入、E4 聚类成立的 OS 模型，但其哈密顿量**没有谱隙**。
2. **等价**：质量间隙 ⇔ E4 聚类具有**体积一致、标号一致的指数速率**（MG）。

因此：**16 个输入不够；必须额外证明 MG。**

**审计更新（2026-11-22）**：16 个外部输入（RG1–5、UV1、PH1–3、AN1、Y1–4、B1、C1；PH4 由 AN1 推出）中**没有**交接项 `q₀ ≥ q₁`。委员会的检查范围表明确不包含质量隙、无限体积唯一性、禁闭或全阶微扰匹配；A5 的 `MassLimit` 注释也写明"no intrinsic scale or spectral gap is asserted"。因此交接不等式是一个**新增输入**，不是隐藏的已有输入。详见 [`research/handoff-audit/16-input-audit.md`](../research/handoff-audit/16-input-audit.md)。

### 3.5 仍然缺什么

| 缺口 | 性质 | 可能的工具 |
|---|---|---|
| MG（一致指数聚类） | 核心未解 | RG + 强耦合交接；或对偶/单极凝聚 |
| 无限体积唯一性 | 未解 | 聚类 + 相唯一性；或 Dobrushin 型准则 |
| 弱→强耦合 RG 交接 | 核心未解 | Balaban 方案的现代化 |
| 4D counterterm 控制 | 技术硬点 | 渐近自由 + 多尺度 |
| 禁闭/面积律 | 未解 | 强耦合 + RG |
| 16 个输入的独立审计 | 工程/数学 | 人工 + Lean |

---

## 4. Lean 形式化现状

### 4.1 mathlib

Douglas 2026 的估计：约 180 万行，覆盖大部分本科数学（测度论、基础概率、Lie 群），但**缺**：
- 多元 Gaussian（截至 2026-10 仍缺）；
- Hilbert 空间上的算子理论；
- Sobolev / 加权 Sobolev 空间；
- 复分析与 Riemann 面的一部分。

Mathlib 有 Schwartz 函数与缓增分布（`Mathlib.Analysis.Distribution.SchwartzSpace`），这是 OS 陈述层的基础。

### 4.2 PhysLean

Joseph Tooby-Smith 主持的物理形式化项目（前身 HepLean）。已覆盖：
- 经典力学、相对论、凝聚态；
- 场论指标记法（arXiv:2411.07667）；
- Maxwell 方程从拉氏量导出；
- Wick 定理（arXiv:2505.07939, 2025）；
- 规划中的规范理论：主丛上的 Ehresmann 联络。

**与 YM 的关系**：提供几何侧（联络、曲率、规范变换）的词汇表，但不涉及构造性 QFT 的测度论。

### 4.3 OSforGFF（最直接可复用的工件）

- 仓库：`github.com/mrdouglasny/OSforGFF`
- 论文：Douglas, Hoback, Mei, Nissim, *Formalization of QFT*, arXiv:2603.15770
- 规模：约 25,500 行 Lean，52 个文件，**0 sorry, 0 axioms**
- 内容：在 `𝒮′(ℝᵈ)` 上构造有质量 Gaussian 自由场测度，证明 OS0–OS4（含聚类与遍历性），d = 2,3,4,5 及所有 d ≥ 2
- 依赖：
  - `bochner`：Minlos 定理、Hilbert–核结构、正定性
  - `gaussian-field`：Schwartz 空间可分性与核性
  - `kolmogorov_extension4`：Kolmogorov 扩张定理
- 关键设计：用 `𝒮` 而非 `𝒟`，因为 `𝒟 ⊂ 𝒮` 且 `𝒮′ ⊂ 𝒟′`，所以 OS 公理蕴含 Glimm–Jaffe 版本
- 已解决的技术点：双重 mollifier 收敛、Fubini 引理、Bessel 协方差、两套 Fourier 约定
- 历史：最初依赖 Minlos、Schwartz 核性、Goursat 三个假设；后来全部证明或绕开

**这是我们的起点。** 它已经把"OS 公理在 Lean 中怎么写"这件事做完了，而且是非空洞的。

### 4.4 4D YM 仓库

- `github.com/yang-mills-group/OpenAI-YM-20261028`（原始提交 + 验证包）
- `github.com/yang-mills-group/YM`（评估、复现脚本、Comparator 配置）

规模与工具链：
- 9,470 个模块；
- Lean 4.32.0；
- Comparator 检查（改动陈述、改定义、塞禁用公理会拒绝）；
- 17 个主定理，标准公理。

**务实的警告**：完整构建大概率需要数小时到数十小时、几十 GB 磁盘。第一阶段不必全量重建，先读陈述签名与 `docs/assessment-1.pdf`。

### 4.5 反面教材：formal-conjectures 的空洞陈述

Google DeepMind `formal-conjectures` PR #7636 试图形式化 Clay YM 陈述，但审稿人发现：
- OS 字段直接写成 `True`；
- 真空方程右边缺失；
- 曲率恒等式断言过强；
- 动力学缺失；
- 群的范围允许任意拓扑。

**教训**：陈述层的每一个字段都必须有非空洞的数学内容，并且要有"反例测试"——自由场应通过，平凡理论应失败，质量间隙必须真的约束谱。

### 4.6 还缺什么（形式化侧）

| 缺口 | 状态 | 优先级 |
|---|---|---|
| OS 重构定理（OS → Wightman） | Douglas 称 2026-10 已有验证包；公开搜索未找到独立确认 | 高 |
| 质量间隙的非空洞陈述 | 未完成 | 最高 |
| 谱隙 ↔ 指数聚类 ↔ 转移矩阵 的等价 | 未完成 | 最高 |
| 格点 YM 强耦合质量间隙 | 未形式化 | 高 |
| 2D YM 精确测度 | 未形式化 | 中 |
| 3D YM / YM–Higgs | 未形式化 | 中 |
| RG 步进定理 | 未形式化（Douglas 说"有地图"） | 长期 |
| Sobolev / 重整化 | mathlib 缺失 | 长期 |

---

## 5. 我们的可行路线

### 5.1 三条候选路线

| 路线 | 内容 | 可行性 | 风险 |
|---|---|---|---|
| **A. 审计** | 读 OpenAI 的 4D Lean 代码，验证 16 个输入是否真的被证明、陈述是否匹配 Clay | 高 | 低；但主要是理解，不产生新定理 |
| **B. 归约 + 强耦合** | 形式化质量间隙陈述，证明三种表述等价，形式化强耦合质量间隙与 2D YM | 中高 | 中；有明确里程碑 |
| **C. 4D 正面进攻** | 直接攻 MG / RG 交接 | 低 | 极高；可能多年无果 |

### 5.2 推荐

**以 B 为主线，A 并行，C 作为长期观察。**

理由：
- B 的每一步都有明确的数学内容和可验收的 Lean 工件；
- B 产出的"陈述层 + 归约定理"是 C 的前提；
- A 能让我们判断现有 4D 构造的真实强度，避免在错误的假设上盖楼；
- C 不应作为近期承诺。

### 5.3 分阶段里程碑

#### M0：环境与工件审计（1–2 周）

- 安装 `elan` + Lean 4.32.0；
- 克隆并构建 `OSforGFF`（先做，规模可控）；
- 读 `yang-mills-group/YM` 的评估报告与 Lean 陈述签名，不做全量构建；
- 产出：`docs/m0-audit.md`，回答：
  - 4D 主定理的精确类型是什么？
  - 16 个输入在 Lean 里是定理还是结构字段？
  - `MassGap` 是否出现在陈述中？如果没有，缺什么？
  - OS 重构是否真的形式化？

验收：能在本地 `lake build` 通过 `OSforGFF`；能引用 4D 仓库的主定理签名。

#### M1：质量间隙陈述层（2–4 周）

建立独立 Lean 包 `YangMillsStatement`：

```lean
-- 草图
structure WightmanTheory where
  H : Type
  hH : HilbertSpace H
  Ω : H
  U : PoincareGroup → H ≃ₗᵢ[ℂ] H
  -- ...
  spectralCondition : ...
  uniqueVacuum : ...
  massGap : ∃ Δ > 0, ∀ ψ ⟂ Ω, ⟪ψ, H ψ⟫ ≥ Δ * ‖ψ‖^2
```

要证明：
1. `massGap ↔ exponentialClustering`（在 OS + 重构下）；
2. `massGap ↔ transferMatrixGap`（有限格点）；
3. **非空洞性**：
   - 有质量自由场满足 `massGap`（复用 OSforGFF）；
   - 无质量自由场不满足（反例测试）；
   - hexagonmath 的"有聚类但无谱隙"模型不满足（如果可形式化）。

验收：0 sorry；三个定理；两个反例测试。

#### M2：强耦合格点质量间隙（1–2 个月）

- 形式化有限格点 YM 的 Wilson 测度；
- 形式化强耦合展开的收敛性（或 Shen–Zhu–Zhu 的 Langevin 路线）；
- 证明：存在 `β₀ > 0`，对 `β < β₀`，转移矩阵有一致间隙。

验收：`strong_coupling_mass_gap : ∃ β₀, ∀ β < β₀, ...`，0 sorry。

#### M3：2D YM（2–4 个月）

- 形式化 2D YM 测度（Migdal 公式 / Driver 构造）；
- 证明 OS 公理；
- 证明面积律与转移矩阵间隙。

验收：2D YM 满足 M1 的 `massGap` 定义。

#### M4：3D YM / YM–Higgs（3–6 个月）

- 走随机量子化或 Balaban 3D；
- 目标不是完整 3D 质量间隙，而是把已知构造形式化到可复用。

#### M5：4D RG 交接（长期，高风险）

- 形式化 RG 步进定理；
- 形式化间隙检验 Δ(n)；
- 证明"Δ(n) 与截断无关"；
- 攻"一次裸尺度估计"；
- 核心难点：从弱耦合流到 O(1) 耦合。

#### M6：整合与审计

- 把 M1–M3 的定理接到 4D 构造的 Wightman 结构上；
- 明确列出"若 MG 成立，则 Clay 问题解决"的形式化条件定理。

### 5.4 诚实预期

- **短期内（3–6 个月）**：可以做出一个有价值的形式化工件——非空洞的质量间隙陈述层 + 强耦合质量间隙 + 2D YM。这足以成为论文或开源贡献。
- **中期（1–2 年）**：可以把 3D 或 2D + Higgs 的构造形式化，并审计 4D 存在性证明。
- **千禧年奖**：需要解决 MG 与弱→强 RG 交接。这是公开难题，不能承诺。我们最可能的贡献是把问题**归约成一个可机器检验的不等式**，并排除陈述层的空洞化。

---

## 附录 A：关键文献

### 官方与综述
- A. Jaffe, E. Witten, *Quantum Yang–Mills Theory*, Clay Mathematics Institute.
- M. R. Douglas, *The Yang–Mills Millennium Prize Problem*, Yang Symposium, Stony Brook, 2026-10-13. [slides](https://media.scgp.stonybrook.edu/presentations/2026/20261013_Douglas.pdf)
- M. R. Douglas, *The Yang–Mills Millennium problem*, Nature Reviews Physics 8 (2026) 86–97.

### 构造性 QFT
- J. Glimm, A. Jaffe, T. Spencer, *The Wightman axioms and particle structure in the P(φ)₂ quantum field model*, Ann. Math. 100 (1974) 585–632.
- J. Feldman, K. Osterwalder, *The Wightman axioms and the mass gap for weakly coupled (φ⁴)₃ quantum field theories*, Ann. Physics 97 (1976) 80–135.
- K. Osterwalder, R. Schrader, *Axioms for Euclidean Green's functions* I, Comm. Math. Phys. 31 (1973) 83–112; II, 42 (1975) 281–305.
- E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, Lecture Notes in Physics 159, Springer, 1982.
- T. Balaban, *Renormalization group approach to lattice gauge field theories* I–IV, Comm. Math. Phys., 1980s.
- J. Magnen, V. Rivasseau, R. Sénéor, *Construction of YM₄ with an infrared cutoff*, Comm. Math. Phys. 155 (1993) 417–430.
- J. Dimock, *The renormalization group according to Balaban* I–III.

### 随机量子化与概率
- G. Parisi, Y.-S. Wu, *Perturbation theory without gauge fixing*, Sci. Sinica 24 (1981) 483.
- A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Langevin dynamic for the 2D Yang–Mills measure*, Publ. IHÉS (2022).
- A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Stochastic quantisation of Yang–Mills–Higgs in 3D*, Invent. Math. 237 (2024) 541–696.
- H. Shen, R. Zhu, X. Zhu, *A stochastic analysis approach to lattice Yang–Mills at strong coupling*, Comm. Math. Phys. 400 (2023) 805–851.
- S. Chatterjee, *Gauge Theory* (book) / *Yang–Mills for probabilists*.

### 面积律与格点
- S. Cao, M. Park, S. Sheffield, *Random surfaces and lattice Yang–Mills*, arXiv:2307.06790, Comm. AMS (2025).
- S. Cao, R. Nissim, S. Sheffield, *Expanded regimes of area law for lattice Yang–Mills theories*, arXiv:2505.16585 (2025).
- S. Cao, R. Nissim, S. Sheffield, *Dynamical approach to area law for lattice Yang–Mills*, arXiv:2509.04688 (2025).

### 二维精确解
- A. A. Migdal, *Recursion equations in gauge theories*, Sov. Phys. JETP 42 (1975) 413.
- D. Gross, E. Witten, *Possible third-order phase transition in the large-N lattice gauge theory*, Phys. Rev. D 21 (1980) 446.
- B. Driver, *YM₂: continuum expectations, lattice convergence, and lassos*, Comm. Math. Phys. 123 (1989) 575–616.
- T. Lévy, *Yang–Mills measure on compact surfaces*, Mem. AMS (2003).
- A. Sengupta, *Gauge theory on compact surfaces*, Mem. AMS (1997).

### 2026 年 AI 相关
- OpenAI, *openai/math* (722 manuscripts), GitHub, 2026-10.
- M. R. Douglas, S. Hoback, A. Mei, R. Nissim, *Formalization of QFT*, arXiv:2603.15770 (2026).
- `github.com/mrdouglasny/OSforGFF`
- `github.com/yang-mills-group/OpenAI-YM-20261028`
- `github.com/yang-mills-group/YM`
- hexagonmath (Chunping Zeng), *From Four-Dimensional Wilson Lattice Gauge Theory to Quantum Yang–Mills Theory*, `hexagon:2611.00492v1` (2026-11-15), hexagonmath 站内预印本；**不是 arXiv**（`arXiv:2611.00492` 是另一篇语言学论文）。
- Yang–Mills Working Group, *Validation of the OpenAI Yang–Mills construction*, 2026-11-10.

### Lean 生态
- J. Tooby-Smith, PhysLean / HepLean.
- *Digitalizing Wick's theorem*, arXiv:2505.07939 (2025).
- *Formalization of physics index notation in Lean 4*, arXiv:2411.07667.
- Google DeepMind, `formal-conjectures` PR #7636（反面教材）.

---

## 附录 B：一页纸速查

```
问题 = 存在性 + 质量间隙
存在性 = 2026 年 Lean 化，条件于 16 个输入，子列极限，未独立审计
质量间隙 = 未解；等价于 MG（体积一致、标号一致的指数聚类速率）
历史主线 = 强耦合（到不了连续极限）→ Balaban RG（UV 可控，IR 缺）→ 随机量子化（2D/3D 成，4D 缺）
2026 模板 = 衰减定理 + 分块 RG + 间隙检验 Δ(n) + 强耦合交接
Lean 起点 = OSforGFF（GFF + OS0–OS4，25.5k 行，0 sorry）
Lean 缺口 = 非空洞质量间隙陈述、等价性定理、强耦合质量间隙、RG 步进
近期目标 = M1 陈述层 + M2 强耦合 + M3 2D YM
```

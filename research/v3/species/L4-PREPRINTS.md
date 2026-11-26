# L4 — 预印本横扫：4D Yang–Mills 质量隙 / 连续极限 / 格距一致下界（2023–2026）

> 检索时间：2026-11-25（本地）。方法：arXiv API `export.arxiv.org/api/query`（逐条 abs 页核 ID）+ `arxiv.org/abs/<ID>` 元数据核验 + PDF 全文抽取（pdftotext）查关键步。
> 纪律：所有 ID 均实际打开过（HTTP 200）；引文逐字取自 arXiv 摘要或 PDF 正文，未凭记忆。
> 结论先行：**(a) 命中——有 3 条预印本自称证明 4D YM 质量隙，全部有致命缺口**；其中 1 条已被 arXiv 管理员撤回。**(b) 同时命中**：多条"提高半径/改进口袋"的严肃结果（SZZ 体系），最远点为强耦合 `'t Hooft` 区 + `U(N)` 随机环境。**(c) 不成立**，故不写"未查到"。

---

## 0. 实际检索式（14 条 arXiv API 查询 + 4 条网页检索，均可复现）

| # | 检索式（原样） | 命中(2023+) |
|---|---|---|
| Q01 | `abs:"mass gap" AND abs:"Yang-Mills"` | 40+（取 80 条按日期倒序） |
| Q02 | `abs:"mass gap" AND abs:"lattice gauge"` | 6 |
| Q03 | `abs:"continuum limit" AND abs:"Yang-Mills" AND abs:"mass gap"` | 3 |
| Q04 | `abs:"uniform in the lattice spacing"` | 1（**假阳性**：讲"uniform range of lattice spacings"） |
| Q05 | `abs:"exponential clustering" AND abs:"Yang-Mills"` | 1 |
| Q06 | `abs:"mass gap" AND abs:"thermodynamic limit" AND abs:"continuum limit"` | 0 |
| Q07 | `au:Chatterjee_S AND abs:"gauge"` | 0 |
| Q08 | `au:Cao_S AND abs:"gauge"` | 0（署名格式不匹配，改用 Q-补） |
| Q09 | `au:Chandra_A AND abs:"gauge"` | 1（无关：4d N=2 Schur 指标） |
| Q10 | `au:Chevyrev` | 18（3D/2D YM，无 4D 隙） |
| Q11 | `au:Hairer_M AND abs:"gauge"` | 0（2023+ 无新条目） |
| Q12 | `au:Bauerschmidt` | 10（全是 φ⁴/sinh-Gordon，**非规范场**） |
| Q13 | `au:Shen_H AND abs:"gauge"` | 4（含 YM–Higgs） |
| Q14 | `abs:"Yang-Mills" AND abs:"spectral gap" AND abs:"lattice"` | 1 |
| W1 | `arXiv 2023 2024 2025 Yang-Mills mass gap proof` | 二手/Zenodo 混杂 |
| W2 | `"mass gap" "4D" Yang-Mills lattice "continuum limit" arXiv 2024` | 二手 |
| W3 | `"mass gap" Yang-Mills arxiv 2025 "uniform in the lattice spacing"` | 1（Zenodo，见 §4） |
| W4 | `"mass gap" site:arxiv.org 2025 gauge`（等价 API 复检 Q01/Q14） | 同上 |

> 注：`export.arxiv.org` 有 429 限流；Q15–Q39 中多条被限流丢弃（已在 §4 用 abs 页逐一补核），未把限流误记为"0 命中"。

---

## 1. (a) 自称"4D YM 质量隙已证"的预印本 —— 逐字引用 + 缺口

### A1. arXiv:2606.19362v1 — Mir Faizal, Arshid Shabir (2026)，hep-lat
`Reflection-Positive Construction of a Four-Dimensional SU(N) Yang-Mills Theory with Mass Gap and Confinement`

**摘要逐字**：
> "the Osterwalder-Schrader reconstruction turns these Euclidean facts into a Minkowski theory with a self-adjoint Hamiltonian, the spectral gap lying above the vacuum and the linear potential for static charges appearing, which gives a concrete picture of confinement"

**声称**：从反射正定格点 Wilson 作用出发，(i) 转移算子有 uniform gap；(ii) 大 Wilson 圈由收敛的 character/polymer 展开给出面积律；(iii) 有限程协变多尺度分析把特征逐尺度搬运，缺陷可求和 ⇒ 连续极限中指数聚类 + 严格正弦张力保持；(iv) OS 重构得自伴 Hamiltonian 与谱隙；(v) 对"容许 slicing"做 telescoping 得唯一极限测度（universality）；(vi) 弱耦合入口 ⇒ 强弱耦合到达同一理论。220+ 页，附 A–G 全部证明；支持性计算另发 Int. J. Geom. Meth. Mod. Phys. (2026) 2650111–2650114。

**缺口（致命，三处）**：
1. **核心一步是"假设"不是"定理"**。PDF §10.1 原文：
   > `Hypothesis 1 (Stated for convenience; proved as Theorem ((10.4)) below))`
   紧接着的一步 RG 估计是
   > `∥Kk+1 ∥ ≤ θ ∥Kk ∥ + c1 gk2 , gk+1 = gk − β0 gk3 + Rk`
   即：**"单步 RG 在 (g,K)=(0,0) 邻域内解析收缩"被当作公设**，随后又以"Theorem (10.4)"自我引用闭合。所有 uniform-in-a 估计最终都回落到这个"admissible scheme"类的公设上（PDF 附录 F 原文：`The estimate Eq.(F.85) is uniform in a and Λ because the constants Cn and Lk are determined solely by the admissible locality/positivity data`）——**uniformity 是类定义的一部分，不是被证明的性质**。
2. **它证明的是"到 Gaussian 不动点的收缩"，不是"非平凡隙"**。`∥Kk∥→0`、`gk→0` 恰是渐近自由 UV 极限；在 `a→0` 下要得到**物理单位**的 `Δ>0`，需要的正是把强耦合格点单位的隙（`~1/a` 发散）搬到弱耦合窗口的**非微扰搬运**。该文从未给出任何 `Δ_phys > 0` 的显式数值下界（无 `Δ ≥ c·Λ_QCD` 形式的定理）。这正是本项目 CLAUDE.md §5 记录的"尚不存在的那台机器"。
3. **面积律/弦张力同样只在强耦合 `β≪1` 有内容**；"弱耦合入口"用的是 `(g,K)=(0,0)` 邻域的小场展开，与 Savvidy 大场凹性交叉点无任何定量控制。

**判断**：不是定理，是"公设化的 RG 程序 + 未证的可加性"。在任何同行评审意义上都**不构成** Clay 解答；其 uniform-in-a 声明正是本题死亡线所描述的那 99%。

---

### A2. arXiv:2307.00788v7 — Adrian P. C. Lim (2023-07 首发，v7 = 2025-01)，math-ph，95 页
`Positive mass gap of quantum Yang-Mills Fields`

**摘要逐字**：
> "The eigenvalue of the corresponding momentum operator will be shown to be strictly less than the eigenvalue of the Hamiltonian, hence showing the existence of a positive mass gap in each component Hilbert space."

同摘要另一句（缺口的直接来源）：
> "the corresponding Hamiltonian eigenvalue is given by the quadratic Casimir operator"

**声称**：构造一个满足 Wightman 公理的 4 维理论：`H_YM(g) = {1} ⊕ ⊕_{n≥1} H(ρ_n)`，按单 Lie 代数 `g` 的不等价不可约表示直和；态由"类空矩形面 S + Lie 代数丛的可测截面 + Minkowski 标架"给出，内积与 S 的**面积**挂钩；用 YM 路径积分 + 重正化 + 渐近自由量子化动量与能量；每块 `H(ρ_n)` 中 `E > |p|`，且各块隙的下确界 `>0`；末节给出隙 ⇒ Clustering。

**缺口**：
1. **隙 = 二次 Casimir 本征值**，是**表示标签（无量纲数）**，不是质量。要变成 `GeV`，必须引入 `Λ_QCD` 或格距 `a`；文中从未把 `a→0` 与 `β(a)`（渐近自由跑动）联系到 Casimir 上。因此"隙 > 0"是无量纲陈述，**没有物理标度**。
2. **理论是按公理"拼"出来的，不是 Wilson 格点 YM 的连续极限**：`H(ρ)` 由面的面积内积直接定义（非可分 Hilbert 空间直和），文中**没有**任何"格点测度在 `a→0` 弱收敛到该理论"的定理，也没有 OS 重构所需的测度存在性论证。即：证明了"存在某个满足 Wightman 且带隙的理论"，**没有**证明"它就是 4D 纯 YM"。
3. 路径积分量子化那一步（§7）用的是未加控制的 `YM` 泛函积分 + "renormalization techniques and asymptotic freedom will be used"的口号式引用，无 `ε–δ` 估计。
4. 版本史可疑：v7 已到 2025-01，`Journal ref` 仍为空（abs 页 `J:` 空），无同行评审记录。

**判断**：**不触及** Clay 问题。它把"隙"偷换成一个离散标签，把"存在理论"偷换成"构造 Hilbert 空间"。可作为"隙的定义必须带标度"的反面教材。

---

### A3. arXiv:2506.00284v2 — D. C. Jacobsen (2025-05)，physics.gen-ph —— **已被 arXiv 撤回**
`A Constructive Proof of Existence and Mass Gap for Pure SU(3) Yang-Mills in Four-Dimensional Space-Time`

**摘要逐字（最狂的一句）**：
> "All steps rest on explicit epsilon-delta estimates and combinatorial bounds for SU(3), leaving no remaining gap between heuristic physics and mathematical proof."

另有两句暴露结构：
> "Our approach embeds the four-dimensional gauge theory as the zero-mode sector of a five-dimensional orbifold regulator that preserves gauge invariance and reflection positivity."
> "a Sturm-Liouville analysis of five-dimensional fluctuations, combined with transfer-matrix spectral projections, isolates a strictly positive glueball mass m0."

**硬事实（比缺口更硬）**：abs 页 Comments 字段现行原文：
> "arXiv admin note: This submission has been withdrawn because it does not meet arXiv's research content quality standards"

**缺口**：
1. **已被撤回**（arXiv 管理员，非作者自撤），且学科分类是 `physics.gen-ph` 而非 `math-ph`/`hep-lat`。任何引用都必须带撤回标记。
2. 方法学：把 4D 理论当作 **5D orbifold 正则化子的零模扇区**——这是**换了一个模型**，不是 Wilson 格点 YM 的 `a→0` 极限；"convergent joint polymer expansion … uniform control of the continuum (a→0) and infinite-volume (L→∞) limits" 是断言而非构造（无 Kotecký–Preiss 型显式判据、无 `β(a)` 窗口）。这与 A1 是同一类洞。
3. "Sturm-Liouville + 转移矩阵谱投影"给出 `m0>0`，但**没有任何 `m0` 与格距无关的量纲化下界**；5D 涨落谱与 4D 胶球谱之间的对应从未被证明。

**判断**：这是 (a) 类声明里最干净的标本——**自称"no remaining gap"，被平台判定质量不达标而撤回**。

---

## 2. (b) 提高半径 / 改进口袋：严肃但只到"强耦合域"

### B1. arXiv:2510.22788v1 — Ron Nissim (2025)，math.PR
`U(N) lattice Yang-Mills in the 't Hooft regime`
**逐字**：
> "We establish a mass gap, prove the existence of a unique infinite volume limit, and give a new proof of the large $N$ limit for $\mathrm{U}(N)$ lattice Yang-Mills theory in the 't Hooft regime."
**证明**：`U(N)` 的 Bakry–Émery 判据失效（Ricci 曲率非一致正），改写成随机环境 `SU(N)` 模型（随机性来自 `U(1)` 场），用 cluster expansion + Langevin 动力学处理 `U(1)×SU(N)`。
**没证明**：仍是 **固定格距**的质量隙；`β` 在强耦合/‘t Hooft 窗口 `|β| < 1/(16(d−1))` 型区间。**最远点 = `U(N)` 群 + ‘t Hooft 标度**，与 `a→0` 无关。

### B2. arXiv:2509.04688v2 — Sky Cao, Ron Nissim, Scott Sheffield (2025)，math.PR（8 页）
`Dynamical approach to area law for lattice Yang-Mills`
**逐字**：
> "The main point is to verify the mass gap condition from [DF80], from which area law directly follows."
**证明**：把 [SZZ23] 的 Langevin 方法用于验证 `[DF80]` 的面积律前提，覆盖 `G ∈ {U(N), SU(N), SO(2N)}`（中心非平凡是 [DF80] 的关键假设）。
**没证明**：**以质量隙为前提**（条件定理），且面积律 ⇒ 弦张力 ⇒ 仍是固定 `a`、强耦合域。**最远点 = 把面积律从"假定隙"降为"验证隙的前提"**。

### B3. arXiv:2401.13299v3 — Hao Shen, Rongchan Zhu, Xiangchan Zhu (2024)，math.PR（60 页）
`Langevin dynamics of lattice Yang-Mills-Higgs and applications`
**逐字**：
> "We establish exponential decay of correlations for a broad class of observables, namely, the infinite volume measure exhibits a strictly positive mass gap."
**证明**：把 [SZZ22/SZZ23] 纯 YM 结果扩到**耦合 Higgs**（紧/非紧靶空间），用分解论证 + 相关性精细分析控制无界 Higgs。
**没证明**：`(β, κ)` 小的强耦合区；无 `a→0`。**最远点 = 加入物质场后仍保持强耦合隙**。

### B4. arXiv:2509.03513v2 — Paulo A. F. da Veiga, Michael O'Carroll (2025)，hep-th/math-ph（20 页）
`On Charge Conjugation, Correlations, Elitzur's Theorem and the Mass Gap Problem in Lattice SU(N) Yang-Mills Models in d=4`
**逐字**：
> "the action associated with the exponent of the exponentiated Haar measure density has a positive local mass term which is proportional to the square of the gauge coupling $β=g^{-2}>0$ and the $SU(N)$ quadratic Casimir operator eigenvalue"
**证明**：`d=4` Wilson 格点上定义电荷共轭算子并给出 `N≠2` 时物理 Hilbert 空间的正交分解；Elitzur 型定理；两个不同场强张量期望为零。
**没证明**：**无任何隙**。仅指出"Haar 测度指数里出现正定局域质量项"（且 `∝1/g²` 反而是**紫外**行为），并由此**预期** `N≠2`、小 `β` 时出现二重简并的单粒子胶球态。属"结构观察"，不是结果。

### B5. arXiv:2301.06996v7 — Puskar Mondal (2023)，hep-th（42 页）
`A Geometric Approach to the Yang-Mills Mass Gap`
**逐字**：
> "assuming the existence of a quantized Yang-Mills theory on $(\mathbb{R}^{1+2},η)$ and $(\mathbb{R}^{1+3},η)$, respectively"
**证明**：轨道空间上的 Riemann 度量 + 正截面曲率 ⇒ 正则化 Bakry–Émery Ricci 曲率 ⇒ 隙（`2+1` 与 `3+1`），**条件于**理论已存在；作者自认 `3+1` 维需额外长度标度、且"remains to be understood in a rigorous way"。
**没证明**：理论存在性；`3+1` 的标度。**最远点 = 几何化条件定理**（与项目"第 6 种物种"思路同族，但恰是条件式）。

### B6. arXiv:2301.04561v5 — V. Gogokhia, G. G. Barnaföldi (2023，v5=2026)
`The Mass Gap Approach to QCD.I.`
**逐字**：
> "Assuming that a non-trivial quantum Yang-Mills theory exists, we have proved that it should have a mass gap $Δ^2 > 0$, indeed."
**证明**：Slavnov–Taylor + 胶子传播子运动方程给出新约束，两解在低能分叉，差异源于 tadpole 常数项。
**没证明**：存在性（自认假设）；且是 QCD 传播子层面的解析论证，非构造性场论。

---

## 3. 明确"不声称"或邻近工作的条目（用于排他）

| ID | 作者/年 | 逐字关键句 | 定位 |
|---|---|---|---|
| arXiv:2611.32690v1 | Yu. N. Berdinsky, 2026, hep-th | "No continuum limit, mass gap, confinement or Makeenko-Migdal equation is claimed; the Giles reconstruction theorem is cited only." | Wilson 圈表述 + Lean 4 机检 9 条陈述；**自承不碰隙**。诚实标本。 |
| arXiv:2611.31309v1 | Victor Atanasov, 2026, physics.gen-ph | "recasting the existence and mass-gap question as a uniform estimate on this cocycle." | 曲率轨道量子化 + Radon–Nikodym 上链；**把问题改写**为一致估计，未证。 |
| arXiv:2608.10133v2 | V. P. Nair, 2026, hep-th | "We give the outlines of a proof for the existence of a mass gap for the 3d Yang-Mills theory in terms of an inequality on the eigenvalues of the Laplacian on the gauge-invariant configuration space." | **3D**，且是 outline。 |
| arXiv:2503.15539v3 | Shoshauna Gauvin, 2025, gen-ph | "The interacting four-dimensional construction, nontrivial renormalized local fields and a uniform physical mass gap remain unproved." | 给出固定 cutoff 下 `≥8a/3`（`b/a≤1/648`）的 Hamiltonian 隙 + "conditional continuum-transfer criteria"；**自承缝未闭合**，是本轮最诚实的一条。 |
| arXiv:2605.02373v2 | Alexander Migdal, 2026 | "providing an analytic mechanism for a purely gluonic mass gap" | twistor-string / Makeenko–Migdal 大绕数 WKB；非严格，属"机制"。 |
| arXiv:2608.02452v1 | Ning Bao, 2026 | "This distinguishes the suppression of global flux sectors from the local spectral information needed to discuss a mass gap." | 量子纠错视角，明确区分"通量抑制"与"谱"。 |

**窗口外的地基文献（2023 前，仅作参照，不计数）**：arXiv:2204.12737（Shen–Zhu–Zhu，强耦合 `|β|<1/(16(d−1))` 下无穷体积极限有严格正质量隙）——所有 B1–B4 的源头。

---

## 4. 被排除的候选（各一句理由）

1. **Zenodo 24117822**（"…a spectral gap ϖ for resampling cubes of physical size l, uniform in the lattice spacing, gives every nontrivial continuum limit the mass gap `16ϖ/(625 e² l)`"）—— 正是最贴近死亡线的"格距一致"声明，但**仅在 Zenodo，无 arXiv ID，无法核可**，按纪律不列为可核条目；若确为 Chatterjee 体系，其"cube gap 在 `β→0`（`a→0`）区是否成立"仍是未证的口袋。
2. **Zenodo 23484057**（"Spacing-Uniform Lattice Inequalities Cannot Produce the Yang-Mills Mass Gap"）—— 与项目"新 no-go"同方向，但同样**无 arXiv ID**，且标题即结论、无可核证明链；仅作线索。
3. **Zenodo 17173707 / 19595333 / engage.cambridge 684bd4ec（"Completion of Agawa's Proof"）** —— 自称无条件格点证明 / "Agawa 补全"，均**无 arXiv ID**，且属典型"自出版 Clay 解答"，按 §1 同标准即判死。
4. **arXiv:2609.31845 / 2609.10841 / 2607.13140 / 2603.17943 / 2405.13696 / 2512.03217 / 2306.00935 / 2412.01313** —— 摘要含 "mass gap" 但分别是黑洞 sphaleron 云、胶球 Van der Waals、`O(6)` σ 模型 transseries、超对称振幅、半经典中心涡旋、半经典双迹形变、矩阵量子力学、QED₃ 有效弦；**均非 4D 连续极限的隙**。
5. **arXiv:2509.09625 / 2510.27611 / 2504.04714 / 2305.04437 / 2303.11016 / 2411.18402 / 2411.14022 / 2412.20794 / 2510.08006** —— 分别是 `U(1)` 单极环、截断 Hilbert 空间压力测试、SU(2)-Higgs 相图、Schwinger 模型、格点 QCD 重子质量、`θ`-依赖谱、Casimir 效应、`SU(N)` 拓扑磁化率；**主题为数值/模型，无 4D 隙定理**。
6. **arXiv:2505.19632（Takook）** —— Krein 空间量子化 + 量子算法，YM 只在背景提及；无定理。
7. **arXiv:2606.06221（Alvarez et al.）** —— `su(2,2|3)` 统一规范理论的手征破缺，是**别的模型**。
8. **arXiv:2512.02107（Chandra–Mukhi–Singh，4d N=2 Schur 指标/MLDE）与 Chevyrev 一批（2608.27828 3D U(1) Langevin、2607.22236 2D YM 测度、2503.03060 3D YM–Higgs 重整唯一性、2302.12160 2D YM 不变测度）** —— 3D/2D YM 或超对称指标，**维度或对象不对**。
9. **arXiv:2506.00284** —— 已列入 §1 A3（撤回），不重复计数。

---

## 5. 裁决

- **(a) 命中。** 三条自称"4D YM 质量隙已证"的预印本：`2606.19362`（Faizal–Shabir，uniform-in-a 靠 `Hypothesis 1` 公设）、`2307.00788`（Lim，隙 = Casimir 标签，无标度）、`2506.00284`（Jacobsen，**arXiv 管理员撤回**）。三条的缺口是同一个：**强耦合格点隙 → 弱耦合连续窗口的非微扰搬运未做**，而它们都把一个"容许类/邻域收缩"的假设当成定理。
- **(b) 同时命中。** 严肃口袋结果的最远点是：`U(N)` + ‘t Hooft 区（`2510.22788`）、YM–Higgs 强耦合（`2401.13299`）、面积律条件化（`2509.04688`）。**全部固定格距**，无一触及 `a→0` 一致下界。
- **(c) 不成立**，故不写"未查到"。
- **对我们项目的直接含义**：2023–2026 的 arXiv 上**不存在**任何"格距一致下界 ⇒ 连续极限质量隙"的真定理；存在的只是三种把假设包装成定理的写法。这从反面确认 CLAUDE.md §5 的判断：缺的正是跨 Savvidy 交叉点的那台机器，而不是又一条固定 `β` 的格点不等式。

### 未能核验项（诚实标注）
- `export.arxiv.org` 429 限流导致 Q15–Q39 部分查询未执行完（Q37–Q39 全丢）；已用 `arxiv.org/abs/<ID>` 逐条补核关键条目，但**不排除**有漏网的 2023–2026 预印本未被上述 14 条式覆盖。
- Faizal–Shabir 的 `journal_ref`（API 报 `Fortsch. Phys. 74 (2026) 4, e70097`, DOI `10.1002/prop.70097`）在 abs 页未显示，未独立核 DOI 真伪。
- Lim 的 `2307.00788` 未见任何期刊版本，亦未见公开的同行评议意见。

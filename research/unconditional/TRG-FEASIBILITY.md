# TRG-FEASIBILITY — 严格认证张量网络能否在 β≈2.5 提取 4D SU(N) 有限体积质量隙？

> 2026-11-24 · 任务：诚实评估 Douglas(2026) 点名的"计算/张量 RG"路线，是否能给出 **逐字认证（interval / error-bounded）**的 4D 纯 SU(2)/SU(3) 格点 Yang–Mills 有限体积质量隙。
> 标签（宪章 §2）：**[T]** 真定理 / **[C]** 条件定理 / **[O]** 精确缺失步（缺但可设想）/ **[NG]** 具名类 no-go / **[V-here]** 本轮新算 / **[R]** 外部据报（本轮已核一手）/ **[unverified]** 未能核实的自出版。**不在未证上敲 QED。**

---

## VERDICT（一句话）

**不可行——不是被某条 no-go 定理禁止，而是卡在"认证"而非"张量网络近似"：4D 网络收缩在 β=2.5 处的 *近似* 计算本身就正好处在当今前沿（D≈10–30，4D 纯规范 SU(2)/SU(3) 尚无人做完，只做完 Z₂/U(1)），而 *逐字认证* 版本额外要求 (i) 一个 d≥2 尚无严格面积律/χ-误差定理的截断界，(ii) 区间误差 < 指数小的 `e^{−ΔL}` 信号，(iii) β=2.5 的规范表示截断尾界——而 (iii) 的严格化本身就是强耦型下界 = 已有的那堵墙 H1。 判定 **[O]**（缺但可设想的计算，非 [NG]），但完整认证流水线在 L=6,8（看标度所必需）处超出当前技术。**

**最可达情形**：SU(2)、`L=4`（或时间向 `L_t=4`）、`j_max=2`、`D≈20`：ATRG 成本 `~5×10¹¹` 运算/步、内存 `~26 MB`，认证误差乐观 `~10⁻²`，信号 `e^{−ΔL_t=4}≈0.13` —— 原理可行。 `L=6,8` 不可行。

**它究竟证什么**：只证 **固定 (a,L) 的有限体积转移矩阵谱隙 Δ(L,a)>0**（等价于该体积的关联长度上界）。**这不是 Clay 隙**：还需 (a) `a→0, L→∞` 且 Δ 一致（=H1）、(b) 连续统重构（=L3 / G1）。有限体积证书是必要非充分。

**酸测**：**过**（自由无质量场无隙 ⇒ 正确证书必须且确实报"无隙"；本轮精确算得自由场隙 `2−2cos(2π/N)→4π²/N²→0`）。

---

## 1. 研究现状（逐条核实；区分"数值精确"与"严格认证"）

### 1.1 严格 MPS/PEPS 近似（有真定理，但只在 1D 或受限类）
- **[R] Hastings 2007**，*An area law for one-dimensional quantum systems*, J. Stat. Mech. (2007) P08024：**1D** 有隙局域哈密顿量基态满足面积律（纠缠有界）——**严格**。
- **[R] Verstraete–Cirac 2006**，*Matrix product states represent ground states faithfully*, Phys. Rev. B **73**, 094423 (2006)：有此面积律则 MPS 忠实表示，`χ ~ poly(1/ε)`——**严格（1D）**。
- **缺口**：**d≥2 无一般严格面积律**（仅自由费米子/特殊类有；一般有隙 2D 系统面积律是公开问题）。故 **2D/3D PEPS 无严格 "χ vs 截断误差" 定理**。这是任何 d≥2 张量网络认证的**根本空白**，不是 4D 规范特有。

### 1.2 收缩复杂度（真定理：精确收缩不可行，必须截断+认证）
- **[R] Schuch, Wolf, Verstraete, Cirac**，*Computational Complexity of Projected Entangled Pair States*, Phys. Rev. Lett. **98**, 140506 (2007), arXiv:quant-ph/0611050：**一般张量网络收缩与经典模拟 PEPS 均为 #P-完全**。⇒ 精确收缩不可行；唯一出路是利用平移不变+有隙结构做**近似+误差证书**，而该结构（面积律/谱隙）恰恰在 d≥2 未证。

### 1.3 TRG 实用算法（全部"数值"、非"认证"）
- **[R] Levin–Nave**，*Tensor renormalization group approach to 2D classical lattice models*, PRL **99**, 120601 (2007)：TRG 原始算法（数值，无认证界）。
- **[R] HOTRG**：Xie, Chen, Qin, Zhu, Yang, Xiang, *Coarse-graining renormalization by higher-order SVD*, Phys. Rev. B **86**, 045139 (2012)。d 维成本 `O(χ^{4d−1})`、内存 `O(χ^{2d})`（本轮据 Lattice2021/Akiyama 2025 报告核）。
- **[R] ATRG**：Adachi, Okubo, Todo, *Anisotropic Tensor Renormalization Group*, arXiv:1906.02007 (Phys. Rev. B **102**, 054432 (2020))。d 维成本 `O(χ^{2d+1})`、内存 `O(χ^5)`（d=4）。
- **4D 数值边界**：HOTRG 4D 实际只到 `χ≈13`；4D 主要工作集中在 **Z₂、U(1)、Z₂ gauge-Higgs**（Akiyama–Kuramashi 等）；**4D 纯 SU(2)/SU(3) 规范 TRG 尚未完成**；3D SU(2)/SU(3) 已做（Atis 等）。[R，据检索]

### 1.4 严格/认证张量数值（正在起步，但不在 4D 规范）
- **[R] arXiv:2611.03241**（2026-11），*Establishing Néel order in the square-lattice AFM Heisenberg model: a bootstrap proof*：用 **interval arithmetic + 验证浮点** 给出**认证变分上界**（ε₀ ≤ −0.668047）。这是"认证张量数值"的**真实范例**，但：(a) 2D 自旋模型（非 4D 规范）；(b) 认证的是**能量上界**（变分），**不是质量隙**；(c) 预印本、未见同行评审。
- **[R] Eisert 等**，*Lower bounds to variational problems with guarantees*, Phys. Rev. A（2026 接受），DOI 10.1103/cz6k-y46r：Anderson 界+SDP 给平移不变周期格点哈密顿量的**下界**——与上界配对可成"认证夹逼"，但非张量网络隙。
- **诚实结论**：**"严格认证的 4D 规范场有限体积质量隙"在同行评审文献中不存在**。

### 1.5 必须排除的自出版"认证"声明（宪章 §2：不传播）
以下均为 **Zenodo 自出版、未见评审**，且多自承"Clay 仍开放"或明言不构造 YM 隙，**不得作为结果引用**：
- "Kozić, *Rigorous Rearrangement Tensor Renormalization in Four Dimensions*"（Zenodo 23016981/23284903）：自承**不构造连续统 QFT、不证 YM 质量隙**；仅高温度平凡不动点附近的局部 RG 定理。
- "Package B — Full QCD Conjecture"（Zenodo 17470707）、"Pichardo 2026 Constructive RG YM S³"（Zenodo 19204109）、"Complete normalization… lattice SU(2)"（Zenodo 23935825）、"Spectral Certificates… O(3)/SU(2)"（Zenodo 23784620）。
**[unverified]** ——列出仅为示警，不作证据。

**1 节小结**：
| 层次 | 对象 | 状态 |
|---|---|---|
| 严格近似 | 1D MPS | **[T]**（Hastings / VC） |
| 严格近似 | 2D/3D PEPS 误差界 | **空白**（无面积律定理） |
| 精确收缩 | 一般 TN | **[T]** #P-完全 ⇒ 不可行 |
| 数值 TRG | 4D 纯 SU(2)/SU(3) | **未做**（只 Z₂/U(1)/Z₂-Higgs） |
| 认证 TN 数值 | 2D 自旋能量上界 | **[R]** 起步（arXiv:2611.03241） |
| 认证 4D 规范隙 | — | **不存在** |

---

## 2. 规范约束障碍（**订正：任务所引"d≥2 no-go"作为具名定理未能核实，且对有限群为假**）

**任务前提**："精确规范不变（Gauss 定律）与有限键维 PEPS 在 d≥2 不相容（Zohar–Cirac–Burrello / Haegeman–Verstraete）"。

**本轮核实结果（宪章 §2：引用逐个核对，抓到错误立即更正）**：
- **[R] Kull, Molnár, Zohar, Cirac**，*Classification of MPS with a Local (Gauge) Symmetry*, **Annals of Physics 386, 199–241 (2017)**（Zbl 1375.81179）：严格**分类**局部规范不变态——但对象是 **MPS（1D）**，**不是 PEPS、不是 d≥2**。
- **[R] Zohar–Burrello**，*Building projected entangled pair states with a local gauge symmetry*, **New J. Phys. 18, 043008 (2016)**：**构造**出有限键维、**精确**规范不变的 PEPS（有限群）。即"有限 D 规范不变 PEPS 在 d≥2 **存在**"——与"d≥2 no-go"字面**相反**。
- **[R] Haegeman, Van Acoleyen, Schuch, Cirac, Verstraete**，*Gauging Quantum States: From Global to Local Symmetries*, **Phys. Rev. X 5, 011024 (2015)**, arXiv:1407.1025：给出**"gauging"构造**（把全局对称态变局域规范态），**不是** no-go 定理。
- **[R] Z2/Z3 纯规范 PEPS 数值（2+1d）**：Emonts, Bañuls, Cirac, Zohar, PRD **102**, 074501 (2020)；Emonts 等 PRD **107**, 014505 (2023)；Kelman 等 PRD **110**, 054511 (2024)：有限 D 规范不变 PEPS 实际可用。

**因此——精确的正确命题（本轮给出，非"具名 no-go"）**：
- **(I) [事实，d-无关]** **连续**规范群的链变量 Hilbert 空间是**无限维**（`L²(G)`，或无穷多个不可约表示）。有限 χ 张量网络的物理腿**有限维**、只能装有限多个表示 ⇒ 对**连续群**的**精确**规范不变在**任何**有限键维下都不可能（**任意 d**，这是维数计数，平凡）。
- **(II) [R]** 对**有限群**（Z_N，或 SU(2) 的有限子群如二十面体群），有限 D 精确规范不变 PEPS **存在**（§2 构造与数值）⇒ "d≥2 不相容"作为**一般**命题**为假**。
- **(III) [R]** 连续群 + 精确规范不变**需要无限键维**（连续 PEPS/CPEPS 的虚拟空间趋于无限维；如 JHEP 07(2026)044 讨论连续 PEPS 虚拟空间"应为无限"；Itou 2025 讲稿："gauge d.o.f. needs an infinite-dimensional Hilbert space"）。

**它禁止什么、允许什么**：
- **禁止**：连续群 SU(2)/SU(3) 的**精确**规范不变在有限 D 下（平凡维数事实）。
- **允许**（须付代价，且代价需认证）：
  (a) **表示截断**：`j ≤ j_max`（quantum-link 截断）→ 对截断群**精确**、对完整 SU(2) **近似**，误差 = 尾部 `j>j_max` 权重，**必须认证**；
  (b) **规范固定**：放弃 Gauss 定律做无规范/固定规范表述 → 需控制 Gribov 副本/FP 行列式，**必须认证**。
- **影响（诚实）**：**这不是硬 no-go，而是一个额外的、必须认证的截断误差源**。但在 β=2.5，`(a)` 的尾部严格界本身就是一个**强耦型下界**（没有小参数），即 **= H1**；`(b)` 的 Gribov 控制在弱耦亦未解决。故规范障碍在本问题里**不是"禁止"，而是"把认证负担再加一层，且这层在 β=2.5 严格化= 那堵墙"**。

---

## 3. "严格的有限体积隙证书"到底是什么（须认证什么）

**目标对象**：固定 `(β, a, L)` 的 4D 环面转移矩阵 `T` 的**谱隙** `Δ(L,a) = −lim_{t→∞}(1/t)log C(t)`，认证 `Δ ≥ Δ₀ > 0`。三者关系：
- **(a) 转移矩阵谱隙** 与 **(b) 关联长度上界**：在有限体积等价（`ξ = 1/Δ`），给出其一即得另一（[T]，有限维线性代数）。
- **(c) 热力学外推 `L→∞`**：与 (a),(b) **不同**——需 Δ(L) 关于 L 的一致下界，**这正是 H1 那一步**（不是有限体积证书能给的）。

**必须认证的四个误差源**（逐个，缺一即非证书）：

| # | 误差源 | 严格界现状（β=2.5） |
|---|---|---|
| 1 | **键维截断** `‖T−T_χ‖ ≤ ε_χ`（算子范数，用于 Weyl 型谱扰动） | **无 d≥2 定理**（§1.1 面积律空白）。**[O]** |
| 2 | **规范表示截断** `‖(j≤j_max)−exact‖ ≤ ε_{j_max}` | 需尾部界 = 强耦型 ⇒ **=H1**。**[O]** |
| 3 | **区间算术收缩** 给出 `T_χ` 的严格围道 | 需误差 `< e^{−ΔL}`（见 §4.3）。技术可达但随 L 恶化 |
| 4 | **`L→∞` 外推** Δ(L) 的一致下界 | = H1。本项目 T1 单调曲线 `A_m̃(T)` 是相关可证伪对象 |

**证书能证的**：`Δ(L,a) ≥ Δ₀ > 0`（有限体积、固定格距）。**证不了的**：`a→0` 与 `L→∞` 下 Δ 一致（H1）、连续统重构（L3/G1）。**故即便全流水线跑通，也只是把 H1 在"第三种技术方言"里再确认一次**，与 `TRANSPORT-RIGOROUS-NUMERICS.md` §4 的预判一致。

---

## 4. 在 β=2.5 的可行性 + [V-here] 计算

脚本：`D:\Dev\math\research\unconditional\w7\trg_feasibility.py`（`py -3.12`，单作业，纯计数/精确谱，无大数组、无 OOM 风险）。β 取标准 SU(2) Wilson 归一 `S=βΣ_p(1−½Tr U_p)`，bulk 转变 `β_c≈2.298`，故 `β=2.5` 正落在标度/禁闭区；`g²=4/β≈1.6=O(1)`（精确值依归一，但"无小参数"在任何归一成立）。

### 4.1 [V-here] 维数账本（周期 L⁴；Kogut-Susskind 链空间，pre-Gauss-law）
局部链维 `dim H_link = Σ_{j≤j_max}(2j+1)²`：`j_max=½→5`，`1→14`，`3/2→30`，`2→55`。

| L | #links | #plaq | slice links | log₁₀ dim（j≤½） | log₁₀ dim（j≤1） | log₁₀ dim（j≤2） |
|---|---|---|---|---|---|---|
| 2 | 64 | 96 | 32 | 22.4 | 36.7 | 55.7 |
| 3 | 324 | 486 | 108 | 75.5 | 123.8 | 188.0 |
| 4 | 1024 | 1536 | 256 | **178.9** | 293.4 | 445.5 |
| 6 | 5184 | 7776 | 864 | 603.9 | 990.3 | 1503.7 |
| 8 | 16384 | 24576 | 2048 | 1431.5 | 2347.3 | 3564.3 |

**读法**：即便 `L=4, j≤½`，pre-constraint 切片 Hilbert 维 `~10^179`；**稠密对角化绝无可能**，TRG 是唯一路，且**必须同时截断键维与规范表示**。

### 4.2 [V-here] 4D 收缩成本（已核标度）
HOTRG d=4：时间 `O(χ^15)`、内存 `O(χ^8)`；ATRG d=4：时间 `O(χ^9)`、内存 `O(χ^5)`。

| χ | HOTRG 时间 | HOTRG 内存 | ATRG 时间 | ATRG 内存 |
|---|---|---|---|---|
| 10 | 1.0×10¹⁵ | 0.80 GB | 1.0×10⁹ | ~0 |
| 13 | 5.1×10¹⁶ | 6.5 GB | 1.1×10¹⁰ | ~0 |
| 20 | 3.3×10¹⁹ | **205 GB** | 5.1×10¹¹ | 0.03 GB |
| 30 | 1.4×10²² | **5.2 TB** | 2.0×10¹³ | 0.19 GB |
| 50 | 3.1×10²⁵ | 312 TB | 2.0×10¹⁵ | 2.5 GB |

**成本墙 #1（HOTRG 内存）**：χ=20 即 205 GB、χ=30 即 5.2 TB ⇒ 工作站不可行；实际 4D HOTRG 只到 χ≈13。**ATRG 内存可行**（χ=20 ~26 MB），但 `χ^9≈5×10¹¹` 运算/步 × `log₂L` 步，**且这只是近似（未认证）成本**；再加：
- `× n_reps`（SU(2) 特征展开每条键带表示标号，`n_reps~2j_max+1`，抬高有效键维）；
- `×2–3`（区间端点存两倍内存）；
- `×~10`（区间舍入/相关性的运算记账）。

### 4.3 [V-here] 认证误差 vs 指数小信号（**最关键的一道墙**）
有限环面时间向 `L_t` 上的隙信号 `C(L_t)~e^{−ΔL_t}`。取 `Δa≈0.5`（β≈2.5 交叉区 O(1) 质量）：

| L_t | `e^{−ΔL_t}`（=所需认证绝对误差上限） |
|---|---|
| 4 | 1.35×10⁻¹ |
| 6 | 4.98×10⁻² |
| 8 | **1.83×10⁻²** |
| 10 | 6.74×10⁻³ |
| 12 | 2.48×10⁻³ |

朴素区间舍入增长 `err ~ (#ops)·u`，`u=2⁻⁵³≈1.1×10⁻¹⁶`，按 χ⁹ 运算/步、8 次扫 × 4 步：

| χ | 每步误差界 | 总（×8×4） |
|---|---|---|
| 10 | 1.1×10⁻⁷ | 3.6×10⁻⁶ |
| 13 | 1.2×10⁻⁶ | 3.8×10⁻⁵ |
| 20 | 5.7×10⁻⁵ | 1.8×10⁻³ |
| 30 | 2.2×10⁻³ | **7.0×10⁻²** |

**读法**：即便**最好情形线性**区间界，χ=20–30 总误差 `~10⁻³–10⁻²`，与 `L_t=8` 的信号 `1.8×10⁻²` **同量级**；且区间界在 4D 收缩里通常**超线性**增长（相关性/包裹效应），而靶子是**两个指数小配分函数之差的指数小量**。⇒ **`L_t=4` 尚可（信号 0.13 ≫ 误差 10⁻²）；`L_t≥8` 处在/超出区间认证边缘**，而看标度正需要 L=6,8。**这是"认证"版相对"近似"版多出来的、且更硬的一道墙。**

### 4.4 ξ/a 与所需 χ
β=2.5 SU(2)：`m_{0++}a` O(0.5–0.7) ⇒ `ξ/a = 1/(m a)` **O(1–3)**（交叉区，与任务估的 2–4 同量级）。有隙边界态所需 `χ ~ 10–30`（PEPS 经验；**最坏精确 `χ~e^{cξ}`**，但以 1% 误差论 PEPS 惯例 10–100）。**这恰好等于现今 4D-TRG 前沿 D=10–30**：**近似问题已在前沿，认证问题在其之外。**

---

## 5. 酸测（自由无质量场必须无隙）

**本轮精确计算**（周期环长 N，转移/拉普拉斯谱）：
- 无质量：`λ_n = 2−2cos(2πn/N)`，隙 `= 2−2cos(2π/N) → 4π²/N² → 0`（N→∞），`ξ=1/√λ₁ ∝ N →∞`。实测：`N=8→0.586, 32→0.0384, 128→2.41×10⁻³, 512→1.51×10⁻⁴`（**无隙**）。
- 有质量 m：`λ_n = m²+2−2cos(…)`，隙 `= m²+λ₁ → m²>0`（`m=0.5` 时 `N=512→0.2502`，**有隙**）。
- **自由 Maxwell（4D U(1)）**：横光子 `ω=|k|`，`k→0` 无隙（精确；无禁闭则无质量生成）。

**酸测结论 [V-here]**：**过**。正确证书对自由无质量场必须且确实报"无隙"（TRG 在自由场处是**临界不动点**，RG 映射不收缩 → 基于收缩的"隙证书"自然返回无隙）。**但须诚实**：过酸测只证明工具**非 vacuous**（对自由场失效），**不**证明它能在 β=2.5 **认证**出隙——后者是 §4 的四道墙。

---

## 6. 最终裁决 + 精确障碍

**问题：β≈2.5 处严格认证的有限 (a,L) 隙，实际可行吗？**
**答：以现有技术，不可行。** 不是被具名 no-go 禁止（**[O]**，非 [NG]），而是卡在"认证"而非"近似"：

**精确障碍（按先后）**：
1. **收缩成本/体积墙**：4D HOTRG 内存 `χ^8`（χ=20→205 GB 不可行）；ATRG `χ^9` 运算/步（χ=20→5×10¹¹，小时–天级 × 区间因子）. 且 **4D 纯 SU(2)/SU(3) 近似 TRG 本身尚未有人做完**（只 Z₂/U(1)），所以不存在可"加认证"的成熟数值基座。
2. **d≥2 无严格面积律/χ-误差定理**（§1.1）：键维截断误差 `ε_χ` **无严格界**——这是任何 d≥2 张量网络**认证**的**原则性空白**。**[O]**
3. **认证误差 vs 指数小信号**（§4.3）：认证区间误差 `~10⁻²` 在最乐观情形与 `e^{−Δ·8}≈1.8×10⁻²` 同量级；`L=6,8`（看标度必需）处在/超出区间边缘。**这是"认证"版独有的、比"近似"版更硬的墙。**
4. **规范表示截断尾界 @β=2.5 = H1**（§2）：`j_max` 尾部严格界本身是强耦型下界；β=2.5 无小参数 ⇒ **=H1**。规范固定路线则需 Gribov/FP 控制，弱耦亦未解。
5. **`L→∞` 一致外推**（§3-4）= H1；**`a→0` 连续统** = L3/G1。

**何谓 [O] 而非 [NG]**：无任何定理禁止"认证 4D 规范隙"；它是**缺但可设想的计算**。**但**其中第 4 项（β=2.5 规范尾界）与第 5 项（L→∞ 一致）**各自都是已命名的开放硬问题**，其严格化即 H1 与 L3。故**完整认证流水线 = [O]，但其两个子步骤是"过不了的已知墙"，而非本路线独有的新障碍。**

**最小可设想情形（若要试）**：SU(2)、`L=4`、`j_max=2`、ATRG `D≈20`：成本 `~5×10¹¹` 运算/步（×区间 ~10 ⇒ 天级）、内存 `~26 MB`；认证误差乐观 `~10⁻²` < 信号 `e^{−Δ·4}≈0.13` ⇒ 原理可行；**但它只证 Δ(L=4,a)>0，非 Clay 隙**，且仍需第 4、5 项才升格。

**一句话定位**：**严格的 4D 规范有限体积 *近似* TRG = 当今计算前沿；严格的 *认证* 版本 = 前沿之外，且其两个子步骤各自撞上 H1/L3 这两堵已命名的墙。**这与 `TRANSPORT-RIGOROUS-NUMERICS.md` §4 的预判（"用第三种技术再测一次 H1，不削弱它"）**完全一致**。**Clay(G2) 仍 0%。**

---

## 7. 外部引用（本轮逐一核实；DOI/arXiv）

- **[R]** Douglas, M. R., *The Yang–Mills Millennium problem*, **Nature Rev. Phys. 8, 86–97 (2026)**, DOI 10.1038/s42254-025-00909-2（2026-01-12）。点名 computational/tensor RG 为物理学家可贡献方向；"no clear path known"。
- **[R]** Hastings, M. B., *An area law for one-dimensional quantum systems*, J. Stat. Mech. (2007) P08024。
- **[R]** Verstraete, F.; Cirac, J. I., *Matrix product states represent ground states faithfully*, Phys. Rev. B **73**, 094423 (2006)。
- **[R]** Schuch, N.; Wolf, M. M.; Verstraete, F.; Cirac, J. I., *Computational Complexity of Projected Entangled Pair States*, Phys. Rev. Lett. **98**, 140506 (2007), arXiv:quant-ph/0611050（#P-完全）。
- **[R]** Levin, M.; Nave, C. P., *Tensor renormalization group approach to 2D classical lattice models*, Phys. Rev. Lett. **99**, 120601 (2007)。
- **[R]** Xie, Z. Y. et al., *Coarse-graining renormalization by higher-order SVD*（HOTRG）, Phys. Rev. B **86**, 045139 (2012)。
- **[R]** Adachi, Y.; Okubo, T.; Todo, S., *Anisotropic Tensor Renormalization Group*（ATRG）, arXiv:1906.02007 (Phys. Rev. B **102**, 054432 (2020))。d 维成本 `O(χ^{2d+1})`。
- **[R]** Haegeman, J.; Van Acoleyen, K.; Schuch, N.; Cirac, J. I.; Verstraete, F., *Gauging Quantum States: From Global to Local Symmetries in Many-Body Systems*, **Phys. Rev. X 5, 011024 (2015)**, arXiv:1407.1025。
- **[R]** Kull, I.; Molnár, A.; Zohar, E.; Cirac, J. I., *Classification of MPS with a Local (Gauge) Symmetry*, **Ann. Phys. 386, 199–241 (2017)**, Zbl 1375.81179（**MPS/1D**，非 PEPS）。
- **[R]** Zohar, E.; Burrello, M., *Building projected entangled pair states with a local gauge symmetry*, **New J. Phys. 18, 043008 (2016)**（有限群有限 D 规范不变 PEPS 构造）。
- **[R]** Zohar, E.; Burrello, M.; Wahl, T. B.; Cirac, J. I., Ann. Phys. **363**, 385 (2015); Wahl, Burrello, Cirac, Ann. Phys. **374**, 84 (2016)。
- **[R]** Emonts, P.; Bañuls, M. C.; Cirac, J. I.; Zohar, E., PRD **102**, 074501 (2020); Emonts 等 PRD **107**, 014505 (2023); Kelman 等 PRD **110**, 054511 (2024)（2+1d Z₂/Z₃ 纯规范 PEPS 数值）。
- **[R]** arXiv:2611.03241（2026-11），*Establishing Néel order … bootstrap proof*（interval arithmetic 认证变分**上界**；2D 自旋；预印本）。
- **[R]** Eisert 等，*Lower bounds to variational problems with guarantees*, Phys. Rev. A（2026 接受），DOI 10.1103/cz6k-y46r。
- **[R]** JHEP 07(2026)044（连续 PEPS，虚拟空间"应为无限维"，与 §2(III) 一致）；Itou 2025 Stony Brook 讲稿（"gauge d.o.f. needs an infinite-dimensional Hilbert space"）。
- **[unverified]** 见 §1.5 的 Zenodo 自出版条目（列出仅为示警，不作证据）。
- **内部**：`TRANSPORT-RIGOROUS-NUMERICS.md`（§4 预判）、`TRANSPORT-FIVE-WAY-VERDICT.md`、`BRIEF-FOR-MODEL.md`、`A5-RIGOROUS-NUMERICS.md`、`PROOF-SKELETON.md`、`CLAUDE.md §5`。
- **[V-here]** `research/unconditional/w7/trg_feasibility.py`（本轮账本/成本/酸测，`py -3.12` 单作业）。

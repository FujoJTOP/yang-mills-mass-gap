# SURVEY — 外部文献横扫：2025–2026 有没有一条真正通向 4D 隙的路？

> 2026-11-24 · 广度扫 arXiv + Zenodo + 期刊 + 预印本平台；对每一条核核实引用、区分 **真定理 / 条件 / 猜测 / no-go**，并逐条过"自由场酸测"。
> 标签：**[T]** 真定理（已证）/ **[C]** 条件定理（含未证前提）/ **[O]** 精确缺失步 / **[NG]** 具名类的 no-go / **[R]** 外部据报（我未逐字核全文）/ **[J]** 判决（本册意见，非定理）。
> **不在未证上敲 QED。** 本节所有 arXiv 条目经 `export.arxiv.org/api` 核实（id/作者/日期/摘要）；所有 Zenodo/DOI 经页面核实。

---

## 0. TOP-LINE（先给结论）

**没有一条给出 4D 隙的新路。** 本轮扫到的 2025–2026 候选里，**零个**是"既真、又新、又落到 4D 物理点"的证明。它们分成五个死法，且**每一个死法都落在本项目的已知墙上**：

1. **固定 cutoff / 只到强耦合一端 + "连续搬运待证"**（Mihajlovic、Gauvin、Cao–Nissim–Sheffield、Nissim、SZZ 线）——把"可证的强耦合端"再往前推一点，但跨交叉点的**体积一致符号/谱下界**仍显式列为未证 = **H1 / G2.5**。
2. **轨道空间曲率（Bakry–Émery 类）**（Mondal 2023、Sakib 2026）——需要轨道空间上**正 Ricci 曲率**（3+1D 还要用胶球质量固定长度尺度 = **循环**），作者自认"heuristic"；这是物种 #2（谱隙-曲率 Γ₂）的轨道空间实例，不是新物种。
3. **动力学学错道**（Lim：**从不用 YM 作用**，内积是面积而非 OS 测度，H=手设 Casimir，且"modified Wightman"；Comitini：胶子质量 = 规范依赖，非规范不变隙）——**未证 Clay，且过不了自由场酸测**（见 §3）。
4. **整条流水线式"全证明"**（Platek Zenodo 114pp、One Cube Zenodo、Schultz–Errante–Hosain G2/SUSY）——要么**文件受限、且自述"从不进入强耦合"**（Platek，= 物种 #13 梯度流 + H1 方向），要么**把立方体估计当公理**（One Cube，"Clay 仍未解"），要么**SUSY 保护在错道 + F₃₇ 数术**（Schultz，无一处承认从 N=1 SYM 到纯 YM 的缺口）。
5. **明确撤稿 / 边缘 / 重构**（Jacobsen 2506.00284 **被 arXiv 管理员撤稿**；"kaleidoscopic/zeta"、"history dependence" 等）。

**唯一在概念上真正"新"的东西**：Bao, *Confinement as Decoding*（arXiv:2608.02452）——把禁闭/格点 YM 读成**高阶形式量子纠错码的解码问题**（Nishimori 系综 + 相对权重）。这是**一个新的类比/透镜**，但**当前只是一个词典，没有任何 4D 隙的 handle**；其"4D ℤ_N memory ↔ confining PSU(N)" 是"possible relation"（猜测）。记为 [O]，**不构成路线**。

**排序（按"真·前景"，全部 LOW）**：Mihajlovic（最诚实、最有量化内容，但自认未建制）> Gauvin（固定 cutoff 证书）> Cao–Nissim–Sheffield / Nissim（'t Hooft 正则律）> Mondal / Sakib（轨道曲率，条件/启发）> 其余（结构/错道/撤稿/边缘）。

**一句话**：2026 年的外部产出**没有越过墙**；它把"墙两侧各自能证到什么"又画细了一圈，但墙本身（在 α_s~O(1) 交叉点给 `0⁺⁺` 符号控制、非循环、(L,a) 一致的下界）**仍是同一堵墙**。与本项目 TRANSPORT-FIVE-WAY 的"五路各证明同一不可约性"一致；Douglas, *Nature Rev. Phys.* **8**, 86–97 (2026) 独立确认"no clear path known"。

---

## 1. 总表（citation | 证了什么 | 新 handle? | 死法 | 判决）

| # | 引用（已核） | 实际证了什么（[T]/[C]/[O] 分级） | 新 handle? | 精确死法 | 判决 |
|---|---|---|---|---|---|
| 1 | **Mihajlovic**, *SU(2) Yang-Mills: Spectral Bounds, Vacuum Concentration, and Criteria for the Continuum Limit*, arXiv:**2611.11534** (2026-11-08), 414pp；数据 Zenodo **10.5281/zenodo.23986334** | **固定空间 cutoff** 的谱界、非平凡局部扇区、真空集中熵界、硬电 cutoff 的**保真度障碍**、固定丢弃链扇区的完整界；把"体积一致隙所需的额外输入"**隔离出来**。**显式自认**：4D cutoff 构造与沿之的**正谱界未建立；不断言精确连续质量** [C]/[O] | **否**——"真空集中/保真度障碍"是**诊断**（告诉你缺什么），不是**下界**；缺的输入 = 体积一致谱界 = **H1** | 体积一致（L,a）符号/谱下界未证 | **最诚实、最有量化内容**；[C] 无假 QED。**不开路**（把墙说得更细） |
| 2 | **Gauvin**, *Pinched Multi Affine Geometry and Confinement*, arXiv:**2503.15539** (2025-03, v3 2026-09-14) | 固定空间 cutoff 下物理 SU(3) Hamiltonian 有**体积无关隙 ≥ 8a/3**（磁/电比 b/a≤1/648）；热力学基态 + 存活有限能 Wilson 激发；正谱表示。**"conditional continuum-transfer criteria"**，且明说"相互作用 4D 构造、重整化局域场、一致物理质量隙**未证**" [C] | **否**——固定 cutoff 隙本质是**自由盒 UV 隙**（对自由场也成立，见 §3）；搬运条件是墙 | 连续搬运与 IR 隙需"uniform comparison"，未给 | 诚实 [C]；**regime-limited**，不开路 |
| 3 | **Cao, Nissim, Sheffield**, *Dynamical approach to area law for lattice Yang-Mills*, arXiv:**2509.04688** (2025-09) | 在 **'t Hooft 参数区**（大 N）验证 DF80 的 mass-gap 条件 ⇒ Wilson 面积律；群 U(N)/SU(N)/SO(2N)（需非平凡中心）[T-regime] | **否**——SZZ 动力学的应用；'t Hooft 正则律 | **大 N / 't Hooft**，非固定 N 物理点 | 真定理但**regime-limited**（大 N）；不碰 N=3 物理点 |
| 4 | **Nissim**, *U(N) lattice Yang-Mills in the 't Hooft regime*, arXiv:**2510.22788** (2025-10) | U(N) 't Hooft 区**质量隙 + 唯一无限体积极限 + 新大 N 极限证明**，用 SU(N) 随机环境 + cluster + Langevin(Bakry–Émery) [T-regime] | **否**——SZZ/Bakry-Émery 物种的 U(N) 修补 | **'t Hooft/大 N**；Bakry–Émery 仍是物种 #2 | 真定理**regime-limited**（已知谱系）；不开路 |
| 5 | **Mondal**, *A Geometric Approach to the Yang-Mills Mass Gap*, arXiv:**2301.06996** (2023, v7 2023-11), 42pp | **轨道空间**（连接模规范）上由约化动能定 Riemann 度量，**若**正则化 Bakry–Émery Ricci 曲率**正**，则给 2+1D、3+1D 隙（**假定已存在量子化理论**）。作者自认"**at least as a heuristic one**"；3+1D 需引入长度尺度"**由最低胶球能量固定**"、"remains to be understood rigorously" [C]/[O] | **否**——轨道空间 Bakry–Émery = 物种 #2 的实例（Nair 2608.10133 的 3D 轨道 Laplacian 之 4D 前驱） | 需**正 Ricci 曲率**（= 未证符号）；3+1D 长度尺度由胶球质量固定 = **循环** | 诚实 [C]；把 3+1D 说成"heuristic"；**不开路**（同一物种） |
| 6 | **Sakib**, *Towards the Yang-Mills Mass Gap via Gauge-Orbit Curvature: A Strong-Coupling Bound and Data-Driven Lattice Evidence*, Cambridge Engage **DOI 10.33774/coe-2026-4q17x** (2026-11-21, 未同行评审) | 对规范不变观测量，Bakry–Émery 曲率比 SZZ 用的 Ricci 曲率**大 (N/2)Q_d**；把可证强耦合区 **β<1/48 → β<195/7168**（4D，因子 1.31）；胶球质量比数据再分析（作者自警"percent-level agreement alone is not evidence"）[C] | **否**——"改进 Bakry–Émery 曲率常数"=物种 #2 的**数值加宽** | **仍只在强耦合端**；无 crossover 论证；factor 1.31 | 真·regime 加宽，但**不开路**（同一物种，微调常数） |
| 7 | **Faria da Veiga, O'Carroll**, *On Charge Conjugation, Correlations, Elitzur's Theorem and the Mass Gap Problem in Lattice SU(N) YM in d=4*, arXiv:**2509.03513** (2025-09, v2 2026-01), 20pp | 物理 Hilbert 空间按 C=±1 **正交分解**（N≠2）；Elitzur 版本；两"不同"场强张量期望 = 0；并指出**指数化 Haar 测度密度**的展开有**正局部质量项 ∝ β·Casimir** [T-structural]/[O] | **否/弱**——"Haar 测度给质量项"是 Lagrangian 二次项，**非谱隙**；结论是"small β 下**expected** 双重态"（猜测） | **层次错**（作用量二次项≠谱原子）+ 只在 small β（强耦合）"expected" | 结构 [T] 但**不证隙**；"Haar 质量项"值得知晓，非新钥匙 |
| 8 | **Bao**, *Confinement as Decoding: Higher Form Codes and Lattice YM*, arXiv:**2608.02452** (2026-08) | 把（Abelian 同调码）解码表述为高阶形式规范场；正局部噪声下逻辑类=Nishimori 系综拓扑扇区；Fourier 关系、"contour/fractional moment"阈值判据；定义 4D ℤ_N memory 与 confining PSU(N) 的"**possible relation**"；有限曲率中心片模型 + 强耦合展开 [O 类比] | **概念新，但无 gap handle**——词典而已；4D 关联是猜测 | 类比/词典；无 4D 隙下界 | **本轮唯一真正新透镜**，记 [O]；**不作路线** |
| 9 | **Lim**, *Positive mass gap of quantum Yang-Mills Fields*, arXiv:**2307.00788** (2023, v7 2025-01), 95pp | 声称构造 4D QFT 满足（**modified**）Wightman 公理 + 每扇区正隙 + inf>0 + Clustering。**见 §2 深度审计** [C→fails] | **否**——**从不用 YM 作用**；内积=面积（非 OS 测度）；H=手设 dim·C₂(ρ) | **动力学学错道 + 自由场酸测不过**（见 §2/§3） | **不是 Clay 证明**（作者自承未用 YM 作用）；精确缺失步见 §2 |
| 10 | **Platek**, *The Yang Mills Mass Gap from Finite-Resolution Coercivity*, Zenodo **10.5281/zenodo.20114011** (2026-05-10, v5, **文件受限**) | 声称**完整证明** 4D SU(N) 有正隙：固定梯度流时刻 s₀（弱耦）coercivity ⇒ Poincaré 隙 + Stroock–Zegarlinski 聚类 + 桥定理 m≥min{λ_SZ,E₀/2} + Mosco 收敛 + subsequential OS。**见 §2 第二审计** [C/fails] | **否**——梯度流 = 物种 #13；Bakry–Émery = #2；Bogomolny/coarea 只给**上界** | **自述"strong-coupling regime is never entered"** ⇒ 要么 UV 自由盒隙（酸测败），要么 Bogomolny 偷给符号（=H1） | **不可核（文件受限）+ 结构性跳过交叉点**；**不开路** |
| 11 | **Bhattacharjee, Bhattacharya, Lockwood, Romanov**, *One Cube Estimate Implies the Yang-Mills Mass Gap*, Zenodo **10.5281/zenodo.23992883** (2026-11-06, v4.37.0, 类型=Software) | 条件链：单立方 Dobrushin–Shlosman 估计 ⇒ 连续极限中的物理隙；Lean4 核验组合/有限算术。**自认未证**：标签的平均协方差、跨尺度有效作用、非平凡连续极限、"**cube estimate 本身是前提（if）**"；"**Clay Millennium problem is open**" [C] | **否**——cube estimate = 一个**单侧/平均协方差界** = H1 方向 | 核心前提未证（=H1）；版本号 4.37 的"软件" | 诚实自标 [C]，**自己说未解 Clay**；辅助工具或可借鉴，**不开路** |
| 12 | **Schultz, Errante, Hosain**, *The Yang-Mills Mass Gap: A Constructive Proof via G2 Holonomy and Supersymmetric Soft Breaking*, Zenodo **10.5281/zenodo.18840503** (2026-03-02) | 声称 **Δ_YM = 6Λ_YM**：G2 holonomy KK 隙 + Veneziano–Yankielowicz + N=1 SYM 有隙（反常匹配）+ "soft SUSY breaking bridge" + OS 重建；C=6 来自 emergent so(6) 在 **𝔽₃₇** 上的 Dynkin 指数。无一处承认缺口 [R/NG-错道] | **否**——SUSY 保护在**错道**（N=1 SYM 隙 ≠ 纯 YM 隙；soft breaking 不保证隙存活） | 从 N=1 SYM 到纯 YM 的**未证跳步**；有限域数术；SUSY 保护在错道 | **不是证明**；引用只为存档，[R]/低可信 |
| 13 | **Jacobsen**, *A Constructive Proof of Existence and Mass Gap for Pure SU(3) YM in 4D*, arXiv:**2506.00284** (2025-05) | 声称 5D orbifold 正则化零模嵌入 + 联合聚合物展开 + OS 重建 + Sturm–Liouville 隔离正胶球质量。**arXiv 管理员注：已撤稿（不符合研究内容质量标准）** [R/withdrawn] | 否 | **撤稿**；5D 零模"隔离正质量"是可疑步（KK 塔 vs 4D 物理隙） | **已撤稿**；不采信 |
| 14 | **Comitini**, *Dynamical gluon mass generation, the Schwinger mechanism... color confinement in linear covariant gauges*, arXiv:**2611.09631** (2026-11-06) | 重审 Kugo–Ojima，指出其推导错误；论证 Schwinger 机制 = Higgs 型（自发破缺全局规范荷），ABFP 束缚态=非物理复合 Goldstone；PT-BFM 零动量胶子传播子给 RG 不变质量标度 [T-side] | **否**——**胶子质量 = 规范依赖**，非规范不变 `f₁²`/谱隙 | **错道**（胶子质量/传播子 ≠ 0⁺⁺ 谱原子） | 非隙证明；**错道** |
| 15 | **Farinelli**, *4D Quantum YM for Weak Coupling: Mass Gap Implies Quark Confinement*, arXiv:**2406.16881** (2024) | 推广 Elitzur，证**若**有 mass gap 则小裸耦下夸克禁闭 [C] | 否 | **假定隙**（前提，非结论）；弱耦 | 非隙证明（条件）；[C] |

> 另有 **R**（外部据报、未逐字核）：Preprints.org 202610.1678「Defect-Platform … conditional SU(3) YM」（摘要自标"conditional、非微扰连续存在仍未解"，**全文 403 未取**）；"Kaleidoscopic Partition Filter / zeta 零点 Δ≈14.13"（ResearchGate，谱同构于 ζ 零点，**边缘**，不采信）；"A Geometric Resolution … History Dependence"（aixiv，**重构非解决**）；academia.edu「stochastic quantization + monopole」（**不可信/AI 生成**类）。

---

## 2. 最深审计（"单条最严肃候选"）

### 2A. 【主审】Lim, arXiv:2307.00788 — 声称 Wightman + 正隙，**实为kinematic错道**

**为什么选它做最深审计**：它是本批里**唯一未被撤稿、非受限、以常规学术口吻声称"满足 Wightman 公理 + 正质量隙 + Clustering"的 arXiv 正刊格式论文**（95pp），即"最像 Clay 级声明"的一条。下列引文均取自 arXiv HTML 正文（v7，CC BY 4.0，逐字核实）。

**它真正构造了什么（[T] 部分）**：一个 Hilbert 空间 `ℍ_YM(𝔤) = ⨁_{n} ℋ(ρ_n)`（每个不可约表示 ρ_n 一个扇区，n=0 为真空线）；每个态是"时空矩形面 S 的可测截面 ⊗ 矩阵 ⊗ Minkowski 标架"的三元组；内积（Eq 2.6）是**对面积 d|ρ'| 的面积分**，仅在两面标架一致处非零；每个 ℋ(ρ) **不可分**（Prop 2.11）。

**精确缺失步（这是审计的价值所在）——三重、且作者自曝**：

1. **从不用 Yang-Mills 作用**。正文原话："**nowhere did we use the Yang-Mills action**"（Rmk 2.10 一带）。一个"Yang-Mills 谱隙证明"**从不使用 YM 作用** ⇒ 它证的不是 YM 的隙。
2. **内积不是物理 OS 测度**。作者自承"**area is not invariant under the action of the Poincare group**"；态只因"测度**motivated** by path-integral flux"才叫 YM 场（Rmk 2.10）。而 Clay 隙 = **OS 重构的物理 Hamiltonian** 在**由格点/连续 YM 测度诱导的** Hilbert 空间上的谱隙。此处 Hilbert 空间是**面积内积的 kinematic 直和**，与 OS 重构**没有建立任何同一性**（正文可见部分**无反射正性/无 OS 重构**）。
3. **Hamiltonian 是手设的，且"隙 = Casimir"**。`Ĥ(ρ), P̂(ρ)` 只作为平移群相位出现，作者自承"**it is not clear how we should choose these numbers**"；随后断言"the corresponding Hamiltonian eigenvalue is **given by the quadratic Casimir operator**"（实为 dim·C₂(ρ)）。即：**"质量"是群表示论标签（Casimir），不是 4D 动力学的结果**。`m₀=inf_n m(ρ_n)>0` 是**表示谱**的下确界，与 β、与交叉点、与 Savvidy 无一关系。
4. **条件性**：摘要明写"**assumes that asymptotic freedom holds … a positive mass gap is implied**"——即便接受其框架，也是**条件于 AF**，而非无条件证明。
5. **"modified Wightman"**：作者承认其构造"**satisfies a modified version of Wightman's axioms**"——**被修改的公理**，不是 Clay 要求的 OS/Wightman。

**判决 [J]**：这**不是** 4D 纯 YM 质量隙的证明。精确缺失步 = **把该面积内积、手设 `H=dim·C₂(ρ)`、modified-Wightman 的直和空间，与"4D YM 测度的 OS 重构"同一化**——这一步**完全缺席**，且作者明说"从不用 YM 作用"使其**原理上不可能**由该构造补上。⇒ **错道（kinematic）/ 未证 Clay**。**并经 §3 自由场酸测判为 vacuous**（H 与 g 无关）。

---

### 2B. 【次审】Platek, Zenodo 10.5281/zenodo.20114011 — 声称"完整证明流水线"，**结构上跳过交叉点**

**为什么次审**：它是本批**唯一提出端到端机器**（coercivity → SZ 聚类 → 桥 → Mosco → OS）的"全证明"声明，114pp。**但文件受限（restricted）——这本身是可复现性硬伤**（[R] 边界：我只能核摘要页，不能逐字核证明）。摘要页给出的结构已被我定位到至少三处与本项目墙正面相撞：

1. **"固定梯度流时刻 s₀、弱耦合"且"strong-coupling regime is **never entered**"**（摘要页原话）。**这就是死穴**：物理质量在 β≈2.5（α_s~O(1)）生成；一个**按自身声明从不进入强耦合**、只在弱耦固定 s₀ 建立 Poincaré/聚类率的方案，**要么**得到的是**自由/UV 盒隙**（格点伪影，不作为物理质量存活于连续极限），**要么**它在"disordered region"用 Bogomolny/coarea 偷渡了符号——**Bogomolny 界是作用量的下界 ⇒ `e^{-S}` 的上界 ⇒ 只给"病态构型权重之上界" = H1 方向，无法下界原子**。
2. **梯度流 = 物种 #13**（本项目已入库：完全单调 = 反射正性，只给上界；酸测败——自由无质量场 `t²E≡1/64π²>0` 也"有隙"）。用它当"coercivity 引擎"**原理上**给不出非自由场的符号下界。
3. **Bakry–Émery 只在"kinetic window 子空间"上有 CD(κ_w,∞)，并自承"absence of a global log-Sobolev inequality on disjoint topological sectors at weak coupling"**。在"kinetic window"外（= 典型场 × k→0 = 唯一符号格）它没有控制；那正是 G2.5/H1 的开格。
4. **常数靠数值脚本验 SU(2)–SU(5)**（计算证据，非数学证明）；**OS 是 subsequential，不断言极限唯一**（OS4 下游）。

**判决 [J]**：即便接受其自述，这仍**不开路**：它把"弱耦 Poincaré 隙"当作物理隙的种子，而**种子在交叉点之外**；跨交叉点的体积一致符号控制（= 禁闭）被"从不进入强耦合"**明确排除**。**要么 vacuous（UV 盒隙/酸测败），要么循环（Bogomolny/coarea 偷给符号 = H1）**。精确缺失步 = **弱耦 s₀ 的 SZ 率到物理 `0⁺⁺` 残数下界的、跨交叉点的非微扰桥**——与 TRANSPORT-FIVE-WAY 的"输运机"是**同一个缺失对象**。注：**文件受限**使其连"逐行找错"都做不了；按宪章 §2 记为 [R]/不可核。

> **审稿方法论对照**：两条"全证明"声明（Lim readable、Platek restricted）各自独立地**在同一处**断裂——Lim 断裂在"空间与测度不是物理 OS 的"（连作用都没用），Platek 断裂在"从不进入质量生成区"。这不是巧合：**它们是"跳过交叉点"的两种语法**。

---

## 3. 自由场酸测（宪章硬条件 4）

**判据**：合法工具必须对**自由无质量场**失效（自由：`f₁²=0`、`G″=0`）。若某"隙工具"在自由场上也报正隙 ⇒ **vacuous**。

| 工具 / 候选 | 自由场上是否也报"隙"？ | 酸测判定 | 说明 |
|---|---|---|---|
| **轨道空间 Bakry–Émery**（Mondal / Sakib；= Nair 物种） | 曲率随 g 缩放，g→0 时 →0 | **过**（对自由失效） | 酸测合格，但**物理点 g~O(1) 的正性未证**；3+1D 长度尺度由胶球质量固定 = 循环。**合格却够不到**。 |
| **梯度流 / finite-resolution coercivity**（Platek；物种 #13） | 自由场 flowed Dirichlet form 也有正 Poincaré 隙（自由盒隙≠0） | **败**（vacuous） | 自由无质量场"有隙"⇒ 该引擎不区分（不）相互作用。 |
| **Lim 的面积内积 + H=dim·C₂(ρ)** | H 与 g 无关，任何单 Lie 代数都赋正隙 | **败**（vacuous；连自由/有质量/相互作用都不区分） | 未用 YM 作用 ⇒ 原理上无法区分自由场。 |
| **固定 cutoff 证书**（Gauvin、Mihajlovic） | 固定 a 下自由盒隙恒正 | **已知·自承**（他们明说连续物理隙未证） | 非"工具"，是 fixed-a 量；诚实；酸测对它们的物理主张**不适用**（它们不主张物理隙）。 |
| **Haar 测度局部质量项**（Faria da Veiga–O'Carroll） | 项 ∝ β·Casimir，β→0 时 →0 | **过**（对自由消失） | 但这是**Lagrangian 二次项**，非谱隙；且 small-β"expected"= 猜测。 |
| **SUSY 保护**（Schultz et al.） | 纯自由 YM gapless；其隙来自 SUSY | **条件/vacuous**（对自由纯 YM 沉默） | N=1 SYM 隙 ≠ 纯 YM 隙；soft breaking 未证保隙。 |

**结论**：**两条"全证明"声明（Lim、Platek）都过不了酸测**——Lim 是 `H` 与 g 无关（vacuous），Platek 是梯度流在自由场也给隙（vacuous）。唯一**过酸测**的真工具是轨道曲率（g→0 失效），但它**在物理点正性未证**且 3+1D 循环。这与本项目"正性类工具天生只能证上界/无序"的元结论一致。

---

## 4. 是否有一条"新物种"？（对照 11 物种死亡图）

逐一映射到已知物种：

| 候选 | 落在哪个已知物种 / 死法 | 是否新物种？ |
|---|---|---|
| Mondal / Sakib（轨道曲率） | 物种 #2 谱隙-曲率（Bakry–Émery Γ₂）的**轨道空间实例** | **否**（同一物种，Nair 3D 的 4D 前驱/加宽） |
| Platek（梯度流 coercivity） | 物种 **#13 梯度流**（本项目 TRANSPORT 已入库）+ #2 Bakry–Émery | **否** |
| Gauvin / Mihajlovic（固定 cutoff 证书 + 搬运判据） | "只到强耦合一端 / 连续搬运待证" | **否**（= 墙体本身） |
| Cao–Nissim–Sheffield / Nissim | "大 N / 't Hooft / 强耦合" | **否**（已知 regime） |
| One Cube（cube estimate → 隙） | "单侧/平均协方差界" = H1 方向 | **否** |
| Faria da Veiga–O'Carroll（Haar 质量项） | "作用量二次项 ≠ 谱原子"（层次错） | **否**（但 Haar 质量项值得记一笔） |
| **Bao（confinement as decoding）** | **量子纠错/高阶形式码**——**不在 11 物种内** | **概念上算"第 12 类透镜"**，但**只是词典、无 handle、4D 关联是猜测** ⇒ 记 [O]，**不构成可证路** |
| Lim（面积内积 + Casimir） | 非标准 kinematic 构造 | **否**（错道，非新数学物种） |
| Schultz（G2/SUSY） | BPS/SUSY 保护在错道 | **否** |

**结论**：本批**没有**产出项目宪章 §1 所要求的"第 6/12 种能下界原子、对自由场失效"的新物种。**唯一概念新意（Bao 解码）尚未发育成工具**；其余全是旧物种的实例/加宽/错道/撤稿。

---

## 5. 对本项目的净输入（诚实清单）

- **[T] 真定理（可借用）**：Cao–Nissim–Sheffield / Nissim 的 't Hooft 正则律（大 N 面积律 + 质量隙），以及 Faria da Veiga–O'Carroll 的 **C=±1 正交分解**与 **Elitzur 版本**（结构定理，可入 G1 侧；注意其"双重态 expected"仍是猜测）。
- **[O] 精确缺失步（与我们一致，无新信息）**：Mihajlovic 的"**体积一致谱界**"、Gauvin 的"**uniform comparison**"、One Cube 的"**cube estimate**"、Platek 的"**跨 s₀/弱耦到物理隙的桥**"——**全是 H1 / G2.5 / 输运机**的别名。外部独立工作把它们各自命名了一遍，**没有提供新钥匙**。
- **[J] 负面确认**：Douglas (2026) 与多篇独立作者（含 Mihajlovic、Gauvin 的自述）一致：**没有已建立的 4D 连续纯 YM 质量隙**。这加强（不替代）本项目的"11 物种×同一不可约性"证据。
- **[R] 可复现性警示**：两条"全证明"（Platek 文件受限；One Cube 类型=Software、版本 4.37）与一条撤稿（Jacobsen）——**声称全证明 ≠ 可核证明**。按宪章 §2，一律不入 [T]。

---

## 6. 一句话判决

**2025–2026 的外部文献没有越过墙：** 唯一概念新意（Bao 的"禁闭=解码"）尚未成为工具；唯一两条"全证明"声明（Lim、Platek）各自以不同语法**跳过交叉点**且**过不了自由场酸测**；其余多条是轨道曲率/梯度流/强耦合/'t Hooft 等**已知物种**的实例或加宽。**墙仍是同一堵墙：在 α_s~O(1) 交叉点给 `0⁺⁺` 残数的、非循环、(L,a) 一致的下界 = 禁闭 = 那台尚不存在的输运机。** **Clay 仍 0%。**

---

### 附：核实与边界

- **已核（arXiv API：id/作者/日期/摘要逐字）**：2611.11534, 2503.15539, 2509.04688, 2510.22788, 2301.06996, 2307.00788, 2509.03513, 2608.02452, 2611.09631, 2406.16881, 2506.00284（撤稿注）。
- **已核（Zenodo/DOI 页面）**：23992883（One Cube，4.37.0，Software，自认 Clay open）、20114011（Platek，v5，restricted）、23986334（Mihajlovic 数据）、18840503（Schultz G2/SUSY）；Sakib DOI 10.33774/coe-2026-4q17x（Cambridge Engage，未同行评审）。
- **已知·不重复（据任务清单）**：Nair 2608.10133（3D 轨道 Laplacian）、Bringmann–Cao–Hairer–Zhao 2611.03640（2D G1）、Kirk Zenodo 18773860、CNS 2505.16585、Federbush 1995、RCS、Wasserstein/Ollivier、Battle–Federbush、rigorous TRG、2026 "2D O(3)" AI 争议。本册仅在**谱系对照**处引用它们。
- **未取全文**：Platek（restricted）、Preprints.org 202610.1678（403）、One Cube（我只核了描述页，未逐字核 Lean 源码）。相关判断标 [R]/[J]，未在其上敲 QED。

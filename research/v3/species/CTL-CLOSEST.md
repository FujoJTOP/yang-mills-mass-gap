# CTL-CLOSEST — 从已发表定理到 Clay 的最短链：现状清单、唯一断环、判定

> 批次 **B3 路**（"绕过新物种"），工作目录 `D:\Dev\math`。
> 日期 2026-11-25。**纯文献 + 推理，未跑任何重计算。**
> 标签：**[R] 已核外部文献** / **[M] 元数据级（Crossref / arXiv API 全记录）** / **[J] 判断（非定理）** / **[C] 条件** / **[?] 声称未核** / **[NG] no-go**。
> 纪律：**不凭记忆写引用**；每条给可点来源；查不到的写"未查到（检索式：…）"；区分"读摘要"与"读全文"。
> 上游：`LIT-CHECK.md`（§Q3、Q3′ 及其引用）、`docs/paper-v2/ch13-machine.tex` §"Why the existing rigorous tools stop"、`docs/paper-v2/ch11b-related.tex`、`research/unconditional/HISTORICAL-TOOL-SWEEP.md`、`research/v3/species/LOOP-RUNBOOK.md` §1。

---

## §1 方法与证据分级（本节决定后面每一句话的强度）

### §1.1 检索做了什么

全部检索于 **2026-11-25** 由 `web_search` / `web_fetch` 实做，端点三类：

1. **arXiv API**（`export.arxiv.org/api/query`）——按 `abs:` / `au:` / `id_list` 抓取，拿到的是**摘要原文 + 元数据**（`journal_ref`、`doi`、`comment`）。用了 6 组查询（详见 §7.2）。
2. **Crossref API**（`api.crossref.org`）——拿到**卷/期/页/日期/DOI/参考文献表**。本轮最有价值的一次调用是 Douglas 2026 综述的 Crossref 全记录：它返回了该文的 **74 条参考文献表**（`reference-count: 74`），这让我能对一个"读不到正文"的综述做**书目级的硬事实陈述**（§2.11）。
3. **`web_search`** 定位于书目补齐与"有没有人反驳"。

### §1.2 阅读深度（三级，逐条标注）

| 级别 | 含义 | 本轮篇数 |
|---|---|---|
| **[F]** | 读了全文 | **0** |
| **[A]** | 读了 arXiv 摘要原文（整段） | 18 |
| **[M]** | 只读了书目记录（卷期页/DOI/参考文献表） | 6 |

**这一条必须写在最前面**：本轮**没有任何一篇读了全文**。因此：

- 凡"这篇的结论精确到什么程度"——我只复述**摘要原文**的措辞，不补细节；
- 凡"这篇的证明内部有没有洞"——一律 **[J]**，且必须给出间接证据链；
- 凡"某定理是否 (L,a)-一致"——若摘要未声明，我写"**摘要未声明**"，不写"不一致"。

### §1.3 "离 Clay 的距离"怎么量（排序泛函）

Clay（Jaffe–Witten 表述）要的是：**4 维、纯 Yang–Mills（无物质）、连续极限存在且满足 OS/Wightman 公理、Hamiltonian 谱在真空之上有严格正下界**。据此定义五个坐标轴，逐条打分（✔ / 部分 / ✘）：

1. **维度** = 4；
2. **群** = 连续紧单群（SU(N)），且**结论与 N 无关**；
3. **耦合区间** = 覆盖物理点（`β≈2.5`，SU(3) 下 `g²=2.4`、`λ_tH≈7.2`）；
4. **体积一致** = 常数不随 `L` 退化（含 `N`-一致）；
5. **含连续极限** = 有 `a→0` 的陈述（不是格点单位的隙）。

排序规则：**(5) 优先于 (3)**。理由：一个只到强耦合但**含 `a→0`** 的结果，比一个覆盖物理耦合但**只给格点单位隙**的结果更接近 Clay——因为前者已经解决了"取极限"这一半，后者还要解决"隙在极限里活不活得下来"，而后者恰恰是格点单位隙在 `a→0` 时**发散**（`m_phys=m_lat/a`）的那个陷阱。这条规则会让 §2 的排序和"常识印象"不同，先声明。

---

## §2 T1 现状清单：已发表的 10 条最强结果（按离 Clay 的距离升序）+ 综述状态（§2.11）

> 排序 = 距离**近**到**远**。每条给：作者/年份/载体 · 结论精确程度（维度/群/耦合/`L`-一致/`N`-一致/连续极限）· 差什么（一句）。

### 2.1 【最近】Shen–Zhu–Zhu 2023 — 4D 纯 YM 的**唯一**有隙定理

- **载体**：*Comm. Math. Phys.* **400** (2023) 805–851；arXiv:[2204.12737](https://arxiv.org/abs/2204.12737)；DOI [10.1007/s00220-022-04609-1](https://doi.org/10.1007/s00220-022-04609-1)。[A]+[M]
- **精确程度**：`d>1` 任意维；`G ∈ {SU(N), SO(N)}`；'t Hooft 标度 `βN`；条件 **`|β| < 1/(16(d−1))`（SU(N)）** 或 `|β| < (N−2)/(32(d−1)N)`（SO(N)）——**与 `N` 无关**。结论：Langevin 动力学在 `Z^d` 上唯一不变测度 + 指数遍历；**LSI + Poincaré**；**无穷体积极限唯一**；**关联指数衰减（严格正质量隙）**；Wilson 环大 `N` 极限可在无穷体积直接取。
- **维度** 4 ✔ · **群** SU(N) ✔ · **耦合** `d=4 ⇒ β<1/48=0.02083`，到 `2.5` 差 **120×** ✘ · **体积一致** ✔（且 `N`-一致）· **连续极限** ✘（隙在格点单位；固定 `β` 下 `a→0` 不构成 AF 连续极限）。
- **差什么**：整条定理住在强耦合相，**从来没跨过 `β≈2.5`，也从未说 `a→0`**。
- **注意**：`LIT-CHECK` §Q3′ 已核定——语料里流通的 `λ_tH<1/24` **在源文中不存在**，正确读法是 `λ_tH = βN < N/48`（N-依赖）。本轮沿用该勘误，未再触发。

### 2.2 Nissim 2025 — `U(N)` 版本（修掉 Bakry–Émery 的失效）

- **载体**：arXiv:[2510.22788](https://arxiv.org/abs/2510.22788)（2025-10-26，25 页，**无 `journal_ref`**）。[A]
- **精确程度**：`U(N)`、't Hooft regime；结论同 2.1（质量隙 + 唯一无穷体积极限 + 大 `N` 极限新证明）。方法：把 `U(N)` 重写成**随机环境 `SU(N)`**（随机性来自 `U(1)` 场），cluster expansion + Langevin 混合。
- **差什么**：**摘要原文自己说明了为什么需要这条弯路**——*"this approach fails because its Ricci curvature is not uniformly positive"*。这句话本身是 §5 的关键证据（Bakry–Émery 在该 regime **不是松，是假**）。
- **同窗口**：仍是 't Hooft regime（强耦合侧），未过 `2.5`，无连续极限。

### 2.3 Cao–Nissim–Sheffield 2025b — 动力学面积律（同窗口，**结论更弱**）

- **载体**：arXiv:[2509.04688](https://arxiv.org/abs/2509.04688)（v2 2025-09-28，8 页，"Submitted"，**无 `journal_ref`**）。[A]
- **精确程度**：把 SZZ 的动力学方法用到 't Hooft regime 证 **Wilson 面积律**；`G ∈ {U(N), SU(N), SO(2N)}`。摘要原话：*"The main point is to verify the mass gap condition from [DF80], from which area law directly follows."*
- **差什么**：**这是面积律，不是隙**。它做的是"把 DF80 的 mass-gap **条件**验证成前提"，不是证出隙。**这条是 §5 R4 的关键证据**：2025 年的进展是横向的。

### 2.4 Cao–Nissim–Sheffield 2025a — 面积律 regime 扩张（改进 OS78）

- **载体**：arXiv:[2505.16585](https://arxiv.org/abs/2505.16585)（v3 2025-09-26，35 页，**无 `journal_ref`**）。[A]
- **精确程度**：纯 `U(N)` 格点 YM 的**面积律**在扩张的参数区成立，**"in particular when N is large"**；改进 Osterwalder–Seiler 1978。方法：把 master loop equation 当**线性非齐次方程**看，证解的先验界；核心novelty 是处理 merger 项（引入 truncated model + 逼近）。
- **差什么**：结论是面积律；`d=4` 大 `λ` / 小 `N` 仍开；与 2.1 同一"强耦合/大 N"侧。

### 2.5 Adhikari–Cao 2025 — 本轮检索到的**唯一**在**弱耦合**侧拿到相关衰减的结果（但群有限）

- **载体**：*The Annals of Probability* **53** (2025) 140–174；DOI [10.1214/24-aop1702](https://doi.org/10.1214/24-aop1702)；arXiv:[2202.10375](https://arxiv.org/abs/2202.10375)v3（"Accepted version. To appear in AOP"）。[A]+[M]
- **精确程度**：**有限**（可非交换）规范群、**弱耦合（大 `β`）**，对一大类规范不变函数（**包含 Wilson 环的任意函数**）证**指数衰减**。方法：把边配置看成格点基本群的同态 + "swapping map"（不用 cluster expansion）。
- **差什么**：**群有限** ⇒ 无渐近自由、无连续极限；且"对每个足够大的 `β` 有隙"**没有跨 `β` 的一致速率**（摘要不声明、结构上也不可能声明）。这条的价值是把"弱耦合半边也有人占着"钉死——**但占的方式不能通向 Clay**。

### 2.6 Chatterjee–Yakir 2025 — Chatterjee 线的最新落点（仍是 abelian + Higgs）

- **载体**：arXiv:[2509.19176](https://arxiv.org/abs/2509.19176)（2025-09-23，38 页，**无 `journal_ref`**）。[A]
- **精确程度**：格点 **YM–Higgs**、**`U(1)`**、`d≥2`、"complete breakdown of symmetry" regime；对**任意** `m>0`，若 `α=mβ` 且 `β` **足够大**，证指数衰减。改进 OS78（后者额外要求 `m` 足够大）；验证 Fradkin–Shenker 1979 的相图。方法：Glimm–Jaffe–Spencer cluster expansion around a massive Gaussian field，**following Balaban et al. (1984)**。
- **差什么**：abelian + 有 Higgs + 大 `β`。**但它是本轮唯一一条"把某个门槛真的往前挪了"的结果**（挪的是 `m` 的门槛，不是 `β` 的门槛）——这个对比在 §5 R4 里很重要。

### 2.7 Chandra–Chevyrev–Hairer–Shen 2024 — 3D YM–Higgs 随机量化（同行评议，Invent. Math.）

- **载体**：*Invent. Math.* **237** (2024) 541–696；DOI [10.1007/s00222-024-01264-2](https://doi.org/10.1007/s00222-024-01264-2)。[M]（Crossref 全记录）+ 二手摘要 [A]
- **精确程度**：3 维 YM–**Higgs** 随机量化方程：定义状态空间与关联 Markov 过程。**不声称隙**。
- **差什么**：维度 3、带 Higgs、无隙声明。是"3D 侧最硬的构造结果"，但不是 3D 的隙。

### 2.8 Cao–Chatterjee 2024 — 3D 欧氏 YM 的**状态空间**

- **载体**：*Comm. Math. Phys.* **405** (2024) 3；DOI [10.1007/s00220-023-04870-y](https://doi.org/10.1007/s00220-023-04870-y)。[M]
- **精确程度**：给出 3D 欧氏 YM 的一个状态空间（测度定义域），**不含隙**。
- **差什么**：维度 3；只解决"住哪儿"，不解决"有没有隙"。
- **3D 总结（任务点名要）**：**未查到 3D 纯 `SU(N)` YM 的严格质量隙定理**。arXiv API 查 `abs:"three-dimensional" AND abs:"Yang-Mills" AND abs:"mass gap"` 共 6 条（全部全息/唯象，无一为定理）；`abs:"lattice gauge" AND abs:"mass gap" AND cat:math.PR` 共 **0 条**。3D 的严格资产是 **状态空间**（2.8）+ **YM–Higgs 随机量化**（2.7）+ **SZZ 的 `d≥3` 分支**（`β<1/32`，同样强耦合）。**"3D 更容易"是物理直觉，不是文献事实。**

### 2.9 Balaban 1987 → Dimock 2013 → Dybalski–Stottmeister–Tanimoto 2024/26（**Balaban 程序的最新状态**）

- **载体**：Bałaban, *Comm. Math. Phys.* **109** (1987) 249–301，DOI [10.1007/BF01215223](https://doi.org/10.1007/BF01215223)；Dimock, "The renormalization group according to Balaban — I. Small fields", *Rev. Math. Phys.* **25** (2013) 1330010，DOI [10.1142/S0129055X13300100](https://doi.org/10.1142/S0129055X13300100)；**DST**, "The Bałaban variational problem in the non-linear sigma model", *Rev. Math. Phys.* **38**(8) (2026) 2461003，DOI [10.1142/S0129055X24610038](https://doi.org/10.1142/S0129055X24610038)，arXiv:[2403.09800](https://arxiv.org/abs/2403.09800)（online 2024-10-30）。[A]（DST）+[M]
- **精确程度**：Balaban = 小场近似下的有效作用生成 + 耦合常数重整化 + 大场界；**UV 侧**，IR 未闭合。DST 2024 = 把 Balaban 程序里最关键的一步（**block averaging 约束下的作用量极小化**，即"变分问题"）第一次在 Balaban 之外严格化：用 Lie 群上的变分法导出临界点方程，再用 **Banach 压缩映射**解它，代价是对格点 Green 函数及其积分核做随机游走展开的精细控制。
- **差什么 / 关键负面发现**：**DST 做的是 `O(4)` 二维非线性 σ 模型，不是 4D YM**。所以任务问的"有没有把 `λ~1` 往前推的 2020s 续作"——答案是：**续作存在，但它往"更简单的模型上把这一步做严格"走，没有往"YM4 的 `λ` 推进"走**。这是一个方向性的、可核的事实（摘要原文：*"We analyze this step for the `O(4)` non-linear sigma model in two dimensions"*）。
- 另一条 2020s 的 Balaban 支线：**Chatterjee–Yakir 2025**（2.6）也是"following the approach of Balaban et al. (1984)"——同样是**往 abelian + Higgs 走**，不是往 `λ~1` 走。

### 2.10 Magnen–Rivasseau–Sénéor 1993 — constructive YM₄（固定 IR cutoff）

- **载体**：*Comm. Math. Phys.* **155** (1993) 325–384，DOI [10.1007/BF02097397](https://doi.org/10.1007/BF02097397)。[M]（Crossref 参考文献表核对；摘要级二手转述见 [NII-REO 记录](https://reo.nii.ac.jp/oja/ART1000345241)）
- **精确程度**：纯 `SU(2)` `YM₄` 在**固定红外截断、无紫外截断**、正则化轴规范、平凡拓扑扇区下的 Schwinger 函数构造基础。
- **差什么**：**质量隙对 cutoff 不一致**——这正是 `ch13-machine.tex` 第 113–114 行那句的文献依据。constructive 线的后继（含 MRS 之后的 Balaban 型工作）**没有**给出 `cutoff`-一致的隙。

### 2.11 综述状态：Douglas 2026（**我读的是书目，不是正文**）

- **载体**：Michael R. Douglas, "The Yang–Mills Millennium problem", *Nature Reviews Physics* **8**(2) (2026) 86–97；DOI [10.1038/s42254-025-00909-2](https://doi.org/10.1038/s42254-025-00909-2)；accepted **2025-11-14**，first online **2026-01-12**；单一作者；**`reference-count: 74`**。[M]
- **正文读不到**：`nature.com` 与 PDF 均被 `idp.nature.com` 重定向拦下（本轮独立复现了 `LIT-CHECK` §2.3 的同一失败）。**因此本文件不引用它任何一句结论。**
- **但书目是可核的硬事实**（这是本轮的方法学收获）：在那 74 条参考文献里，**没有任何一条是 4D 纯 YM 质量隙的证明**。我能逐条点出的"最前沿"是：SZZ、CNS 面积律、Park–Pfeffer–Sheffield–Yu（*Probab. Math. Phys.* **7** (2025) 37–121）、Chandra et al.（Invent. Math. 2024）、Cao–Chatterjee（CMP 2024）、Bałaban（CMP 1987）、Dimock（RMP 2013）、**DST（RMP 2026）**、**Rai et al.（Quantum 2026）**、**Ebel–Kennedy–Rychkov（arXiv:2506.03247）**、Meurice 的张量 RG 综述（*Rev. Mod. Phys.* **94** (2022) 025005）、以及 Chatterjee 2025 的 CMSA 讲稿。
- **[J] 由此可下的唯一推断**：一份 2026-01 的权威综述，其参考文献表里**没有"已证"这一项**；它列的是**候选路线**。任务书里"Douglas 明说 no clear path known"这句我**仍然无法核实**（且 `LIT-CHECK` 已记录 Meurice 的二手转述呈**乐观**语气）——**不要引那句话**。
- **一个额外收获**：Douglas 把 **Rai et al. 的谱隙证书层级**、**EKR 的张量 RG**、**DST 的 Balaban 变分问题**同时列进参考文献，说明**这三条正是圈内认为"还没被用尽"的东西**。这直接支撑 §6 的 T4 选择（不是我拍脑袋挑的）。

---

## §3 T1′ 2024–2026 的"声称"清单 + 同行评议状态（任务硬要求，逐条核）

> 口径：**"声称"≠"结果"**。本节只回答两件事：(a) 有没有经过同行评议；(b) 我核到什么深度。

| # | 声称 | 载体 | 同行评议？ | 我的核实深度 | 判定 |
|---|---|---|---|---|---|
| 1 | 4D `SU(N)` YM 的反射正构造 + 质量隙 + 禁闭 | **Faizal & Shabir**, *Fortsch. Phys.* **74** (2026) 4, e70097；DOI [10.1002/prop.70097](https://doi.org/10.1002/prop.70097)；arXiv:[2606.19362](https://arxiv.org/abs/2606.19362) | **是**（Fortschritte der Physik 是同行评议刊；Scopus 有记录） | **[A] 只读摘要**；PDF 无 HTML 版（`arxiv.org/html/2606.19362v1` 返回 404） | **[?] 声称未核**。三条间接反证：(i) 摘要**通篇没有**给出任何显式 `β` 区间、常数、或 `(L,a)` 一致性陈述，只有"interlaced inequalities whose small defects can be summed"这类定性语；(ii) Douglas 2026 的 74 条参考文献**不含它**；(iii) Zenodo 上有第三方审计（灰文献）明说其未给出 Wilson-specific 估计。**不采信，但也不宣布为假**（我没读全文）。 |
| 2 | 纯 `SU(3)` YM₄ 的构造性存在 + 质量隙 | **Jacobsen**, arXiv:[2506.00284](https://arxiv.org/abs/2506.00284)v2 | **否** | [A] | **已撤稿**。arXiv 管理员注记原文：*"This submission has been withdrawn because it does not meet arXiv's research content quality standards"*。**决定性否定**。 |
| 3 | 4D 量子 YM 的正质量隙（Wightman 公理） | **Lim**, arXiv:[2307.00788](https://arxiv.org/abs/2307.00788)v7（95 页，v7 2025-01-17） | **未查到**（arXiv 元数据**无 `journal_ref`**） | [A] | [?]。构造依赖"可数多个不可分 Hilbert 空间的直和"与"面积给内积"，摘要未给 `(L,a)` 一致陈述；**无同行评议记录**。 |
| 4 | QCD 基态必有 mass gap（ST 恒等式路线） | **Gogokhia & Barnaföldi**, arXiv:[2301.04561](https://arxiv.org/abs/2301.04561)v5 | **未查到**（无 `journal_ref`） | [A] | [?]。前提是"假设非平凡量子 YM 存在"，属 **[C] 条件**陈述。 |
| 5 | 固定 cutoff 下 `SU(3)` Hamiltonian 的**体积无关**隙 | **Gauvin**, arXiv:[2503.15539](https://arxiv.org/abs/2503.15539)v3（physics.gen-ph） | **否**（gen-ph，无 `journal_ref`） | [A] | **诚实的非声称**。摘要原文：`b/a ≤ 1/648` 时隙 `≥ 8a/3`，且**自陈**："The interacting four-dimensional construction, nontrivial renormalized local fields and a uniform physical mass gap **remain unproved**." ⇒ 这是**固定 `a`、小 `b/a`** 的结果，正好落在 §1.3 的第 5 轴 ✘。 |
| 6 | 3D YM 质量隙的"证明纲要" | **Nair**, arXiv:[2608.10133](https://arxiv.org/abs/2608.10133)v2 | **未查到** | [A] | **不是证明**。摘要自称 *"outlines of a proof"*，工具是"规范不变构型空间上 Laplacian 特征值的不等式"，**未给该不等式**。 |
| 7 | C-宇称扇区分解 + Elitzur + "预期"二重简并 glueball | **Faria da Veiga & O'Carroll**, arXiv:[2509.03513](https://arxiv.org/abs/2509.03513)v2 | **未查到** | [A] | **结构性结果，非隙**。真结论：`N≠2` 时 Hilbert 空间按 `C=±1` 正交分解；两个不同场强张量的期望为零；指数化 Haar 密度带**正定局域质量项** `∝ β·C₂`。隙只是 *"is expected to be present"*。 |
| 8 | 灰文献（Zenodo / PhilArchive / viXra）的 4D 证明 | 多条 | **否** | 未逐条读 | **不计入证据**（沿用 `LIT-CHECK` §2.3 的判定）。本轮新增两条同性质记录（Zenodo 上"Yang–Mills Mass Gap Follows from One Cube"自陈"the Clay Millennium problem is open"；另一条自陈"proves no mass gap … conditional theorem"）。**灰文献自己都在自陈未闭合**。 |
| 9 | "某机构 2026-10 已证 4D gauge+fermion 存在性" | 仅博客级页面 | **否** | 只见到二手叙述 | **无 arXiv 记录、无期刊记录、无可核预印本**（检索式：`OpenAI Yang-Mills proof October 2026 arXiv`）。按 IRON RULE **不计入**。 |

**§3 小结（一句话）**：2024–2026 关于 4D 纯 YM 质量隙的**声称**很多，**经同行评议且被独立采信的一条也没有**；唯一进入同行评议期刊的那一条（#1）在其摘要层面**不含任何可检验的定量断言**，且被领域综述的书目表排除在外。

---

## §4 T2 最短链：`A(已证) → B(已证) → … → Clay`，逐环标注

### §4.1 选哪条

T1 里只有**一条**链的每一环都是**关于 4D 纯 YM 测度本身的定理**：**SZZ 线**（2.1 → 2.2 → 2.3/2.4），配 Balaban 线的 UV 半边（2.9）。其余各条（2.5 有限群、2.6 abelian+Higgs、2.7/2.8 三维、2.10 固定 IR cutoff）在**第一步**就换掉了问题（换群 / 换维度 / 换 cutoff 结构），所以它们不是"更远的路"，而是"**别的路**"。

**结论：最短链 = SZZ 线 + Balaban UV 半边。**

### §4.2 链（六环）

```
L1  [已证·定理]  OS78 + SZZ23
     d=4, G=SU(N)/SO(N), |β|<1/48：LSI + Poincaré（L-一致、N-一致）
     ⇒ 唯一无穷体积测度 + 指数簇聚 ⇒ 格点单位下的严格正质量隙
        │
        ▼
L2  [已证·定理]  Nissim25（U(N)，随机环境 SU(N)）+ CNS25b（同窗口面积律）
     把"群"这一维补齐；把 DF80 的 mass-gap 条件在 't Hooft regime 验证
        │
        ▼
L3  [已证·定理]  Balaban87 + Dimock13 + DST24/26
     UV/小场多重标度控制；block-averaging 变分步已在 2D O(4) NLSM 上严格化
     ⇒ 覆盖 AF-UV 锥（g_eff→0）；**不覆盖典型场**
        │
        ▼
L4  [★唯一断环★]  典型场 × k→0 单元的单向方差界（= H1）
     存在 λ>0，使对 (L,a) 一致、每个 k→0 的粗横向 v：
        Var_μ(⟨v,∇V⟩ | ξ=0) ≤ E_μ[H_vv | ξ=0] + W_full − δ
     等价形式：D ≤ S − λ 且 S ≥ λ；等价形式：G''_a / a² ≥ Δ₀²
        │
        ▼
L5  [条件·由 L4 推出]  L4 ⇒ 指数簇聚 ⇒ OS 重构 ⇒ Wightman 理论 + 谱隙
     （OS 反射正性对 Wilson 作用是定理；a→0 由 L4 的 (L,a) 一致性与 AF 锥给出）
        │
        ▼
L6  Clay
```

### §4.3 逐环标注（**该环是否已证 / 若未证，要的是"更强估计"还是"新定理"**）

| 环 | 状态 | 若未证，需要什么 | 依据 |
|---|---|---|---|
| **L1** | **已证** | — | SZZ 主定理（2.1） |
| **L2** | **已证** | — | Nissim 摘要、CNS25b 摘要（2.2、2.3） |
| **L3** | **已证（在其声明范围内）** | 其"未闭合"部分是 IR，不是 UV；UV 侧的缺件**不是瓶颈**，因为沿 AF 流 `g_eff→0`，那里**有小参数** | Balaban/Dimock/DST（2.9）；`ch13-machine.tex` 第 54–56 行"the AF-UV cone by `g→0`" |
| **L4** | **未证 = 唯一断环** | 见 §5：**新定理（新机制）**，不是更强估计 | `ch13-machine.tex` 第 60–73 行 eq.(sharpest)、第 144–158 行 (H1)/(NA)/(TCI) 三式（eq.(sharpest) 与 (H1) 是同一条不等式的两种写法） |
| **L5** | **未证，但是 L4 的条件推论** | 给定 L4，只需 OS 重构（定理）+ `a→0` 的一致性（L4 提供） | `ch13-machine.tex` §"What building it would close"：*"A single one-sided variance upper bound … would close both G1's L3 (UV) and G2's H1 (IR) simultaneously."* |
| **L6** | — | — | — |

### §4.4 为什么"唯一"是 L4，而不是"L3 也没闭合"（这一条必须自己拆掉）

很容易反驳说：L3 也没闭合（Balaban 从未完成 4D YM 的构造），所以至少有**两个**断环。我的回应，以及它的诚实边界：

1. **L3 的缺件与 L4 的缺件是同一个缺件的两个坐标。** L3 缺的是"大场/典型场的控制"；L4 缺的是"典型场 × `k→0` 的方差上界"。Balaban 的小场展开在其适用域内是**完备**的（给到任意阶的界），它在典型场处失效**不是因为算得不够，而是因为展开参数是裸耦合 `g²=O(1)`**。所以 L3 的洞**就是** L4 的洞在 UV 坐标下的像。
2. **`ch13-machine.tex` 的单元分解把这一点说死了**：整个收敛被归约到 `(场 × 动量)` 平面上的**一个单元**——"typical field × high momentum"；其余单元"low field by `Φ''⪰N/12+` BL、high-field cone by WSTC、AF-UV cone by `g→0`"。若接受这个分解，断环唯一。
3. **诚实边界 [J]**：上面第 1、2 条**都是我接受的内部论证，不是我的定理**。"UV 端的缺件可以由 L4 补上"这件事本身**在文献里没有被在 `β=2.5` 处做过**——Balaban 的所有严格陈述都在小裸耦合。所以我给的判定是：**逻辑上唯一断环 = L4；工程上第二风险 = L5 的 `a→0` 子步**。两者不是并列的"两个断环"，因为后者没有独立的数学内容（它是 L4 的一致性在极限里的读法）。

### §4.5 一个必须记录的逻辑旁路（以及它为什么也死在 L4）

有人会说：不必跨过 `β=2.5`，可以证"有隙的 `β` 集合在 `(0,∞)` 上**既开又闭**"——开性是 cluster expansion 的定理，闭性用"隙不能连续地关掉"（Kato/解析性），于是 `β=2.5` 自动被覆盖。**这条路是死的，且死法已知**：`HISTORICAL-TOOL-SWEEP.md` §2 第 7 行（"complex analysis"）已经把它钉死——*"Analyticity (Kato) transports the **value** of `a(β)`, not its **non-vanishing**: a real-analytic `a≥0` positive at both `1/48` and `2.5` may still touch `0` between."* 也就是说：**闭性步骤需要"`a` 在中间不为零"，而这正是要证的东西。** 旁路不通 ⇒ 只能正面攻 L4。

---

## §5 T3 判定：L4 需要的是"已有技术的更强估计"，还是"目前不存在的东西"？

### §5.1 判定

> **新定理（新机制）。** 但必须加限定：**不存在性是技术性的（没有载体），不是原理性的（没有形状）。**

下面给理由，不给结论式断言。

### §5.2 R1 — 参数结构：三条主力的控制参数在 `[1/48, 2.5]` 上都是 `O(1)` 或更差，**不是"松"**

- **簇/聚合物展开**：SZZ 的阈值 `1/(16(d−1))` 是**组合常数**给的（`β ×` 单 plaquette 上聚合物计数的上界 `< 1`），与 `N` 无关。`d=4` 给 `1/48 = 0.02083`，到 `2.5` 差 **120×**。这不是"界松了 120 倍"，而是"**级数在 `β=2.5` 处根本不收敛**"——展开参数本身是 `O(1)`。
- **Bakry–Émery Γ₂**：需要 `Hess S ⪰ ρ > 0`。Savvidy Hessian 是 `−D² + 2ad F_B`，**不定**（负模横向磁）。在该方向 `ρ` 不是"变小"，是**变号**。`CLAUDE.md §5` 记录 β=2.5 时负模占比 **7.1%**（且随 β 恶化）。**外部佐证**：Nissim 摘要原文——`U(N)` 的 Ricci 曲率 **not uniformly positive**，所以要绕道。**"按定义失效"和"界不够紧"是两件事。**
- **Balaban 型 RG**：展开参数是**裸耦合**。`β=2.5` 在 SU(3) 下 `g²=6/β=2.4`、`g≈1.55`、`λ_tH≈7.2`——**没有小场区域**，典型场就是 `O(1)`。
- **大 `N`（Chatterjee 线）**：`1/N` 控制的是**涨落**，不是固定点问题；且结论仍在强耦合（2.6 的最新落点还退到 abelian + Higgs）。

⇒ 三条主力在目标 regime 的控制参数**没有一条**处在"再优化一个数量级就能过"的位置。

### §5.3 R2 — 方向：L4 要的是**单向下界**，而现有工具的产出全是别的方向

- **二阶/方差型**（Bakry–Émery、Brascamp–Lieb、LSI、Poincaré、Christoffel、Szegő）：给的是**上界**（对方差）或**下界**（对倒数），换算到 `inf supp ν` 的方向都是错的。
- **计数/整数型**（谱流、index、Birman–Krein、Levinson）：给的是**整数**，被连通 `β`-线上的**刚性**冻结成自由场的值 `0`。
- **等价式**（Tauberian/Bernstein）：`−(1/t)log G(t) → inf supp ν` 是**恒等式**，不是可改进的不等式。
- **内部已完成的穷尽性证据**：`HISTORICAL-TOOL-SWEEP.md` 扫了 13 类、~1800→今的工具，结论是 **Family I 全部因 Schur 凹（线性数据对 `inf supp` 无下界信息）死；Family II 方向对、过酸测不空洞，但全部因输入循环死（`E₂` 回归 / Haag / 刚性 / 无穷体积连续性）**。
⇒ 要补 L4，必须**翻转不等式的方向**且**不引入循环输入**。这是**换不等式**，不是**换常数**。

### §5.4 R3 — 形状已知、载体未知

L4 的五个条件（非微扰 / 单向 / 非循环 / 反自由场 / `(L,a)` 一致）已经把对象钉成一个**具体的数学对象**（`ch13-machine.tex` 第 12–30 行、第 144–158 行），且它的六个等价形式（H1、(CP)、edge atom、`G''`、行列式零点、Savvidy 符号）**互相等价**——这说明形状是**稳的**，不是我在拼凑。但 `LIT-CHECK` §2.2 的五机制核查给出 **0/5 被 gap 语境占用**（其中 2 条在别处被占用且占用手法恰好抹掉信号，1 条原料本身未证，1 条方向反）。

⇒ 缺的是**载体**（哪条不等式能同时满足五条件），不是**形状**。这正是"技术性不存在"的含义。

### §5.5 R4 — 经验反驳（**最关键的一条，因为它是可核的事实而不是论证**）

**若 L4 只是"更强估计"，那么 2023→2026 应该看到阈值的纵向推进。实际一次也没有。**

- SZZ 的 `|β| < 1/(16(d−1))` 从 arXiv v1（2022-04-27）到 CMP 正刊（2023）到本轮检索日（2026-11-25），**一步未动**。
- 2025 年的四篇同线工作**全部是横向的**：
 - **群**方向：`U(N)`（Nissim）、`SO(2N)`（CNS25b）；
 - **结论**方向：面积律（CNS25a/b）——**面积律 ≠ 隙**，CNS25b 自己说它的作用是 *"verify the mass gap condition from [DF80]"*，即把一个**假设**验成条件；
 - **耦合坐标**方向：'t Hooft regime（Nissim、CNS25b）；
 - **维度**方向：`d≥2` 统一（SZZ 本来就 `d>1`）；
 - **加物质**方向：`U(1)`+Higgs（Chatterjee–Yakir，2.6）。
- **对照实验（同一批人、同一个方法族）**：Chatterjee–Yakir 2025 确实把一个门槛挪了——挪的是 **Higgs 质量 `m`** 的门槛（从"`m` 足够大"到"任意 `m>0`"），**不是 `β` 的门槛**。这说明这条技术线**有能力**推进门槛，但它能推进的那个门槛**不是**我们需要的那个。
⇒ 三年、零纵向进展、120 倍的 `β` 差距：这是"机制而非常数"的**经验证据**。

### §5.6 R5 — 什么会推翻我这个判定（可证伪性）

如果出现下列**任一**，R1 被削弱、L4 退回"更强估计"，我的判定就要改：

1. 在 `[1/48, 2.5]` 的**任何**子区间上给出 `LSI` 或 `Poincaré` 的**一致常数**（哪怕只 `d=4, SU(2)`、哪怕常数极小）；
2. 把 Balaban 的小场展开在 `g²=O(1)` 处**重新组织**成收敛的（例如用"场强而非场"作展开参数）；
3. 在有限群情形（2.5）证明速率对 `β` **有正的一致下界**——那会说明"大 `β` 侧的衰减"可以脱离"有限群"这个限制。

**目前为止：三条全部"未查到"**（检索式见 §7.2）。

### §5.7 判定的一句话形式

> L4 需要的是**一条在 `O(1)` 耦合处仍成立、且方向朝下、且输入不循环的不等式**。现有技术在 `[1/48, 2.5]` 上的控制参数**全部**是 `O(1)` 或更差（簇展开的 120× 收敛失败、Bakry–Émery 的 `ρ` 变号、Balaban 的裸耦合 `g²=2.4`），因此这**不是**把常数调紧能到达的；但它的**形状**已被五条件钉死、**候选载体**在 2020s 文献里至少有 3 个（§6），因此它也**不是**"原理上不知道要找什么"。

---

## §6 T4 三条"最可能被忽略的已有工具"

> 选择标准（三条同时满足）：(a) 2020s 的、**同行评议或高可信**的结果；(b) 与本链的**某一具体环**形状对得上；(c) 有一条**可指名的**"为什么还没人用"的结构性理由。
> 一个外部约束：本节 T4-2（Rai et al.）与 T4-3（EKR）**都出现在 Douglas 2026 的 74 条参考文献表里**（DST 也在，见 §2.9），说明"圈内已认为它们可能有用、但没人做"——这不是我的偏好。

### T4-1 **Polchinski 方程的 LSI 判据**（Bauerschmidt–Dagallier）

- **来源**：Roland Bauerschmidt & Benoit Dagallier, "Log-Sobolev inequality for near critical Ising models", *Comm. Pure Appl. Math.* **77** (2024) 2568–2576；DOI [10.1002/cpa.22172](https://doi.org/10.1002/cpa.22172)；arXiv:[2202.02301](https://arxiv.org/abs/2202.02301)。[A]
- **它做什么**（摘要原文要点）：对耦合矩阵谱半径有界的铁磁 Ising 模型，证 log-Sobolev 常数有一个**只依赖磁化率**的界；由此得 LSI 常数**在体积上一致、直到临界点**，"without using any mixing conditions"；若磁化率满足 mean-field 界，则 LSI 常数对到临界点的距离**多项式依赖**（`d>4` 时适用）。证明用：**Polchinski（重整化群）方程表述的 LSI 判据** + 一条 Ising 关联不等式 + Perron–Frobenius + 乘积 Bernoulli 的 LSI。
- **它可能补哪一环**：**L1 → L4 的 LSI 搬运**。SZZ 的 LSI 只到 `β<1/48`，Bakry–Émery 在 `O(1)` 处**按定义失效**；这是目前唯一一条"**不用曲率、用 RG 方程**拿 LSI"的现成机器，形状正好对上"把 LSI 搬过 crossover"。
- **为什么还没人用**：它的核心输入是一条**正关联（FKG 型）关联不等式**；`SU(N)` Wilson 测度的 Weingarten 权重**变号**（`LIT-CHECK` §2.1(c)，Collins–Matsumoto–Novak 综述为证），FKG 的前提**在第一步就不成立**。**但**：判据只用到"某条相关不等式"，若能在 YM 上找到一条**符号无关**的替代——例如只对**中心 `Z_N` 部分**用 FKG（`Z_N` 是 abelian，FKG 可用），把 `SU(N)` 的非交换性留在 cluster 项里——这一步可能被绕开。**这正是"值得试、但没人试"的确切位置。**

### T4-2 **谱隙证书的 SDP 层级**（Rai–Kull–Emonts–Tura–Schuch–Baccari）

- **来源**：*Quantum* **10**, 2065 (2026)；DOI [10.22331/q-2026-04-13-2065](https://doi.org/10.22331/q-2026-04-13-2065)；arXiv:[2411.03680](https://arxiv.org/abs/2411.03680)v3（v3 2026-04-01）。[A]
- **它做什么**：把"**热力学极限下**的谱隙下界"写成**一列优化问题（半定规划）**，层级越高证书越强；把 Knabe 界及其后续改进作为**特例**包含进来（摘要原话：*"guaranteed to either match or surpass them"*）；在 1D 自旋链上，下界精度与"能检出隙的参数范围"都提高**数个量级**。
- **它可能补哪一环**：**L4 的"证书形式"**。本仓库的 `F3∧F4=∅` 说"可计算证书类为空 / 满类等于隙本身"（见 `LIT-CHECK` §Q4，以及正确的收紧形态"**族上不可一致可计算**"）。这条层级给了一个**中间态**：不是一条不等式，而是**一列单调改进、可计算、可验证的下界**。如果能把 YM 的转移算子塞进去（哪怕只是固定 `a` 的 `SU(2)`、或 1+1D 截断），"可计算但无用"与"有用但不可计算"之间就有了**定量刻度**——这直接是 B3 路要的"绕过"素材。
- **为什么还没人用**：层级要求 (i) **frustration-free**、(ii) 有限局域维数、(iii) 热力学极限的局域拓扑。Wilson 转移算子**不是** frustration-free（Kogut–Susskind 的动能项与磁能项不对易，规范约束把局域 Hilbert 空间变成无穷维）；论文的实证只在 1D。⇒ 移植前必须先做一个"**规范固定 + 有限截断 + frustration-free 化**"的模型——**这个活本身没人做过**。

### T4-3 **严格张量 RG + 计算机辅助**（Ebel–Kennedy–Rychkov）

- **来源**：Nikolay Ebel, Tom Kennedy, Slava Rychkov, "Tensor Renormalization Group Meets Computer Assistance", arXiv:[2506.03247](https://arxiv.org/abs/2506.03247)v2（v1 2025-06-03，v2 2026-06-22，45+5 页，含开源代码）。[A]
- **它做什么**：引入 **2×1 重整化映射**（一个方向粗粒化 2 倍 + 90° 旋转），发展一套图形语言把映射作用翻译成**张量分量上的不等式系统**（Hilbert–Schmidt 范数下的严格估计），定义有限维"hat-tensor"包围盒与主函数；**数值迭代 + 严格误差界**证明：落在可量化邻域（63 个正交扇区内偏差 `≤0.02`）内的张量流到高温不动点；对 2D Ising / XY 给出显式高温相区界。
- **它可能补哪一环**：**"搬运"这件事的严格形态**。L4 的本质是"把强耦合端的结论搬到另一端"；EKR 示范了"**RG 映射 + 数值迭代 + 严格包围**"可以在**没有小参数**的地方产出定理——这正是 L4 需要的**形态**（而不是又一条不等式）。
- **为什么还没人用**：4D 的键维数/截断误差**完全不受控**（`0.02` 的邻域是 2D 的），且已证的收敛只到**高温（=强耦合）不动点**，也就是我们**已经站着的那个端点**。⇒ 它是"**形态可用、4D 不可用**"的工具。但它把"严格搬运"从口号变成了一份**可执行模板**，而且 Douglas 2026 把它列进了参考文献——说明圈内注意到了、没人做 4D。

### §6.4 备选（记录但不推荐）

**Bao, "Confinement as Decoding: Higher Form Codes and Lattice Yang–Mills Theory"**, arXiv:[2608.02452](https://arxiv.org/abs/2608.02452)（2026-08-03，hep-th，无 `journal_ref`）。[A]
它给的是"contour 与 fractional moment 的**阈值判据**"，且**自己点出了与 L4 的差别**：*"This distinguishes the suppression of global flux sectors from the local spectral information needed to discuss a mass gap."* ⇒ 它能补的是"全局通量扇区抑制"，**不是**"局域谱信息"。**不推荐**它作为 L4 载体，但推荐它作为"阈值 vs 隙"这条界线的一份清晰表述。

---

## §7 结论、未查到清单、核对记录

### §7.1 结论（三句）

1. **最近的路** = **SZZ 线 + Balaban UV 半边**：`OS78/SZZ23（β<1/48 的 LSI+Poincaré+唯一无穷体积+指数簇聚）→ Nissim25（U(N)）→ CNS25a/b（面积律）→ Balaban87/Dimock13/DST24（UV 小场）` ⇒ **`L4`：典型场 × `k→0` 的单向方差界（H1）** ⇒ OS 重构 + `a→0` ⇒ Clay。
2. **断在 `L4`**，且**只断在 `L4`**（`L5` 是它的条件推论；`L3` 的洞是 `L4` 的洞在 UV 坐标下的像）。
3. **`L4` 需要的是"目前不存在的东西"（新定理 / 新机制），不是"已有技术的更强估计"**——但它的不存在性是**技术性**的（没有载体），不是**原理性**的（形状已被五条件钉死）。最强证据是**经验性的**：SZZ 阈值 2022→2026 **一步未动**，同期四篇同线工作**全部横向扩展**（群 / 面积律 / 't Hooft / 加 Higgs），而面积律**不等于**隙。

### §7.2 未查到清单（明写检索式，**不补**）

| 未查到的东西 | 检索式 |
|---|---|
| Douglas 2026 综述**正文**（含"no clear path"原句） | `nature.com/articles/s42254-025-00909-2`（→ `idp.nature.com` 重定向）；`Douglas "The Yang-Mills Millennium problem" Nature Reviews Physics full text` |
| **3D 纯 `SU(N)` YM 的严格质量隙定理** | arXiv API `abs:"three-dimensional" AND abs:"Yang-Mills" AND abs:"mass gap"`（6 条，全为全息/唯象）；`abs:"lattice gauge" AND abs:"mass gap" AND cat:math.PR`（**0 条**） |
| **把 SZZ 阈值从 `1/48` 往前推的任何论文** | `Shen Zhu Zhu strong coupling threshold improvement`；`lattice Yang-Mills mass gap beta larger than 1/48`；`improved strong coupling radius lattice gauge mass gap 2025` |
| **Faizal–Shabir 的同行评议裁决**（referee report / erratum / retraction） | `Faizal Shabir Yang-Mills Fortschritte der Physik erratum retraction`；只查到 Zenodo 灰文献的第三方审计（明说其未证 Wilson-specific 估计），**非同行评议，不作证据** |
| **"2026-10 某机构已证 4D"的可核预印本** | `OpenAI Yang-Mills proof October 2026 arXiv`；只有博客级页面，**无 arXiv / 无期刊记录** |
| **"`F3∧F4=∅` 的正确形态（族上不可一致可计算）"的等价文献** | 沿用 `LIT-CHECK` §Q4 的检索式，本轮未新增 |

### §7.3 核对记录与风险（可审计）

- **本轮读全文的篇数：0。** 全部为摘要级（18 篇）+ 书目级（6 篇），篇目见 §2/§3 逐条的 `[A]` / `[M]` 标注。凡"证明内部有洞"的判断均为 **[J]**，含对 Faizal–Shabir 的判定——依据是三条**间接**证据（摘要无显式区间/常数、Douglas 书目不含它、灰文献审计），**不是**我读了全文发现的洞。
- **沿用而非重核的两条**：`λ_tH<1/24` 在源文中不存在（`LIT-CHECK` §Q3′）；`(NA)` 已被双向否证、`(TCI)` 双读法判死、只剩 H1（`LOOP-RUNBOOK.md` §1）。本轮**没有**重新触发这两条。
- **本轮新增的可核事实（4 条）**：
 1. Douglas 2026 的 Crossref 记录给出 `reference-count: 74`，且**该 74 条里没有一篇 4D 纯 YM 质量隙的证明**（§2.11）；
 2. **Balaban 程序的 2020s 续作（DST 2024/26）做的是 `O(4)` 2D NLSM，不是 YM4**（§2.9）——这直接回答了任务的一个问题，且答案是负面的；
 3. **Adhikari–Cao 2025（AOP 53, 140–174）** 是本轮检索到的唯一在**弱耦合**侧拿到相关衰减的结果，但**群有限**（§2.5）；
 4. **arXiv:2506.00284 已被 arXiv 管理员撤稿**（§3 #2）——这是 2024–2026 声称潮里唯一一条**决定性的**状态记录。
- **一条自我限制**：§1.3 的排序泛函把"含连续极限"排在"覆盖物理耦合"之前。这个偏好会让 SZZ 排第一；若改成"覆盖物理耦合"优先，排第一的会是 **Adhikari–Cao**（弱耦合但有限群）——但那条路**结构性**地到不了 Clay（有限群无 AF），所以我认为当前排序是对的。这条偏好是 **[J]**，写出来供反驳。

### §7.4 一句话

**今天离 Clay 最近的路是 SZZ 那条（`β<1/48` 的 4D 纯 YM 有隙定理）+ Balaban 的 UV 半边；链在第四环 `L4` 断掉，断点是 `(L,a)` 一致的"典型场 × `k→0`"单向方差界（H1）；那一环要的不是把某个常数调紧 120 倍，而是一条在 `O(1)` 耦合处仍然朝下、且输入不循环的新不等式——而 2023–2026 三年里，这条链的阈值一步没动，动的全是横向。**

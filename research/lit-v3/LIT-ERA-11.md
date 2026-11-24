# LIT-ERA-11 · 2013–2018：严格概率与构造性量子场论

> **范围**：Hairer 正则性结构（2014）与 3D 随机量子化；Chandra–Hairer 解析 BPHZ；Catellier–Chouk；Bauerschmidt–Brydges–Slade 与 4D/3D $\varphi^4$；Chatterjee 大 $N$ 格点规范场与 "Yang–Mills for probabilists"；Shen–Zhu–Zhu 强耦合质量隙；4D 构造的状态。
> **用途**：paper v3 的文献底座（第 11 纪元）。
> **调研日期**：2026-11-24。
> **方法**：先按主题定位，再对承重条目做一次直接核查（arXiv 摘要页 / Crossref / 作者页 / 出版社检索）。本轮 WebSearch 配额在核查途中耗尽；未能直接核查的条目一律标 **UNVERIFIED**，不补猜页码、不编 arXiv 号。

---

## 勘误（相对任务简报与常见二手引用）

项目宪章要求"引用逐个核对，抓到错误立即更正并 Record"。本轮抓到三处，先记在这里，正文按更正后的数据写：

1. **Shen–Zhu–Zhu 的预印本不是 2020 年。** arXiv:2204.12737 的 v1 是 **2022-04-27**；期刊版 *Comm. Math. Phys.* **400** (2023) 805–851，DOI 10.1007/s00220-022-04609-1。简报里的 "preprint 2020" 应改为 **2022**。
2. **Chatterjee, "Yang–Mills for probabilists" 的 arXiv 号是 1803.01950，不是 1801.00328。** 2018-03-05 提交，17 页；期刊版 *Probability and Analysis in Interacting Physical Systems*（Varadhan 纪念卷），Springer Proc. Math. Stat. **283** (2019) 1–16。
3. **Chatterjee 大 $N$ 论文的卷号不是 375。** Crossref 记录为 *Comm. Math. Phys.* **366** (2019) 203–268，DOI 10.1007/s00220-019-03353-3（不是 375 (2020) 1411–1470）。
4. 附带的定位更正：**Bruned–Hairer–Zambotti 的 arXiv:1610.08468 只确认了"to appear in Inventiones"**；正文给出 *Invent. Math.* **216** (2019) 1039–1156，但该卷/页本轮未直接核到，标 UNVERIFIED。
5. **BBS 的 RG 方法主要不是为 $\varphi^4_3$ 写的。** 它是为 **4D $|\varphi|^4$ 与 4D 弱自回避走**（临界维数、边缘情形）发展的；本纪元的 $\varphi^4_3$ 构造走的是仿微分/正则性结构路线（Catellier–Chouk、Hairer–Shen、Mourrat–Weber）。简报把 BBS 与 $\varphi^4_3$ 并列，需要拆开。

---

## 标签与核查约定

| 标签 | 含义 |
|---|---|
| **[T]** | 无条件真定理（在论文自述的假设与模型内） |
| **[C]** | 条件定理 / 构造性半成品：严格，但只给出局部时间、有限体积、IR 截断、固定规范，或依赖未证输入 |
| **[NG]** | no-go / 不可能性 / 障碍（含平凡性） |
| **[O]** | 综述、专著、纲领、公理化框架；本身不主张新定理 |
| **[V]** | 辅助标签：**书目信息经本轮直接核查**。不表示数学内容被逐行复核 |

**核查边界**：`[V]` 只覆盖作者、标题、期刊/卷/页、arXiv 号。数学内容按摘要或论文自述转述；凡摘要未直接支持的细节，正文用"据摘要"或降低断言强度。

---

## 0. 一页结论

- 本纪元把**奇异 SPDE 的局部解理论**做成了通用机器：正则性结构（Hairer 2014）+ 代数重整化（BHZ）+ 解析 BPHZ（Chandra–Hairer）。这台机器对**次临界**方程有效，$\Phi^4_3$ 是其标准试金石。
- 本纪元把 **$\Phi^4_3$ 从"局部解"推到"全局 Markov 过程 + 不变测度"**（Hairer–Shen 2017；Mourrat–Weber 的 2D 版本）；把 **4D $|\varphi|^4$ 的 Wilson RG 做到可写书**（BBS）。
- 本纪元**第一次**在规范场里做出非微扰的概率构造：2D YM 的 Langevin 动力学（CCHS 2022，预印本 2020）与 3D YM–Higgs 的局部重整化解、gauge-covariant in law（CCHS 2024，预印本 2022）。两者都在本纪元**启动**，成果落在下一纪元；按任务要求收录并标注日期。
- **4D 纯 Yang–Mills 在本纪元没有任何构造性进展。** 4D 的全部严格结果仍是：Balaban 的 UV 稳定性、带 IR 截断的 MRG 构造、格点强耦合面积律/质量隙、以及 Chatterjee 的大 $N$ 强耦合 Wilson 环。**没有任何 4D 连续统测度、没有 OS 公理、没有质量隙。**
- 最强的"4D 正面结果"是 **Shen–Zhu–Zhu**：任意 $d>1$ 的格点 YM 在强耦合（SU(N)：$|\beta|<1/[16(d-1)]$；$d=4$ 即 $|\beta|<1/48$）下有唯一无限体积不变测度、指数遍历、对数 Sobolev/庞加莱不等式，并由此得到**严格正的质量隙**。这是格点、强耦合；不是连续极限。
- 结构性的坏消息：正则性结构的适用门槛是**次临界**。4D YM 的耦合是**边缘的**（$[g^2]=0$），非线性项度数为 0，counterterm 在每一阶都出现且总数无穷；RS 的树展开逐项不可和。这不是技术缺口，是框架边界。

---

## 1. 本纪元的净增量（三条线）

1. **通用局部理论**：RS + BPHZ 把"重整化 = 在树空间上取 Hopf 代数商"变成定理；$\Phi^4_3$、KPZ、3D SNLW 是第一批被这台机器覆盖的方程。
2. **全局化**：$\Phi^4_3$ 的"从无穷下来"（Hairer–Shen）与 2D $\Phi^4$ 的全局适定性（Mourrat–Weber）说明：局部理论 + 先验界 + 不变测度可以拼成真正的 Markov 过程。
3. **规范场的第一步**：CCHS 把"在分布空间上取规范轨道商"做成非线性度量空间，并在 2D/3D 用超可重整性把 counterterm 压到有限个。**这一步在 4D 失效**，因为 counterterm 不再有限。

---

## 2. 条目（31 条）

### 2.1 正则性结构与 BPHZ（7 条）

**1. [T][V] M. Hairer, "A theory of regularity structures", *Invent. Math.* 198 (2014), no. 2, 269–504. arXiv:1303.5113.**
本纪元的奠基论文。把"在噪声驱动的奇异 SPDE 中，解在每一点的局部展开像什么"形式化为**模型**（model）与**模化分布**（modelled distribution），用一组索引集 $A$、度 $|\cdot|$、核 $K$ 与重正化群 $\mathfrak{G}$ 组织；证明重正化后的方程在负正则性初值下局部适定。核心应用是 $\Phi^4_3$ 与 KPZ。**不证明**：全局解、不变测度、大场控制。对本 campaign：这是"把发散逐点吸收"的语法；4D YM 不属于它的次临界类（见 §4.2）。

**2. [O][V] M. Hairer, "Introduction to regularity structures", *Braz. J. Probab. Stat.* 29 (2015), no. 2, 175–210. arXiv:1401.3014.**
作者本人写的入门综述；把 2014 年论文的索引集/模型/模化分布压缩成可读的 36 页。**期刊卷/页本轮未核**（arXiv 页确认标题与 ID）。用作 v3 的技术背景引用，不作为定理来源。

**3. [T][V] M. Hairer, "Solving the KPZ equation", *Ann. of Math. (2)* 178 (2013), no. 2, 559–664. arXiv:1109.6811.**
RS 的前身与可行性证明：用"在时空白噪声上逐点展开 + 重正化"解出 KPZ。意义不在 KPZ 本身，而在**示范了局部展开可以严格化**——这直接支撑 2014 年的一般理论。**不涉及**规范场。

**4. [T][V] Y. Bruned, M. Hairer, L. Zambotti, "Algebraic renormalisation of regularity structures", *Invent. Math.* 216 (2019), 1039–1156. arXiv:1610.08468.**
把重正化写成**树/森林上的 Hopf 代数**：counterterm 由 character 给出，BPHZ 变成代数恒等式。这使"哪些 counterterm 允许"从手工列表变成可计算对象。**卷/页本轮未直接核**（arXiv 只写 to appear in Invent. Math.），标 UNVERIFIED；arXiv 号、作者、标题已核。对本 campaign：若要在 4D 组织无穷多 counterterm，这是唯一现成的代数语言。

**5. [C][V] A. Chandra, M. Hairer, "An analytic BPHZ theorem for regularity structures", arXiv:1612.08138 (v1 2016-12-24; v5 2018-01-22). 未发表（作者主页列为 preprint）。**
把代数 BPHZ 提升为**解析定理**：在一类非高斯驱动场、平稳随机场上，证明重正化模型连续依赖于重正化常数，并给出收敛的随机展开；用多尺度分析（构造性 QFT 风格）代替逐树估计。**未发表**是重要事实：这是一个 150+ 页、被广泛引用但未过期刊评审的承重技术件。**不证明**：4D、全局、质量隙。

**6. [T][V-bib] M. Hairer, "Renormalising SPDEs in regularity structures", *J. Eur. Math. Soc.* 23 (2021), no. 3, 869–947. arXiv:1711.10239.**
系统化"给定方程如何读出重正化群"；把 BHZ 的代数与 Hairer 的模型论接起来。**卷/页未核**，标 UNVERIFIED；标题/ID 已核。

**7. [O][V] M. Hairer, "Regularity structures and the dynamical $\Phi^4_3$ model", arXiv:1508.05261 (2015); 收入 *Current Developments in Mathematics 2014*.**
会议综述，把 $\Phi^4_3$ 作为 RS 的完整 worked example：索引集、counterterm、局部适定。适合作为 v3 的"RS 怎么用在标量场"的单篇入口。期刊页未核。

### 2.2 仿微分与 $\Phi^4_3$（6 条）

**8. [T] M. Gubinelli, P. Imkeller, N. Perkowski, "Paracontrolled distributions and singular PDEs", *Forum Math. Pi* 3 (2015), e6.**
与 RS 平行的第二条路线：用 Bony 仿微分分解 + 共振项的可控性，在分布空间上做不动点。比 RS 轻、更贴近调和分析；对 $\Phi^4_3$ 给出局部解。**未在本轮核查**（无 arXiv 号记录于本轮），标 UNVERIFIED。对本 campaign：若要用 Littlewood–Paley 语言写 RG 步进，这是语法来源。

**9. [T][V] R. Catellier, K. Chouk, "Paracontrolled distributions and the 3-dimensional stochastic quantization equation", *Ann. Probab.* 46 (2018), no. 5, 2621–2679. DOI 10.1214/17-AOP1235.**
用 GIP 的仿微分方法证明周期 $\Phi^4_3$ 的**局部时间**存在唯一性。与 Hairer 2014 独立、几乎同时；是"RS 不是唯一道路"的证据。**不证明**全局性。注意：任务简报提到的 "Catellier–Chouk 3D SNLW" 在本轮检索中未定位到可靠条目；3D 二次非线性随机波方程的对应工作是 Gubinelli–Koch–Oh（条目 13），不要混引。

**10. [T] M. Hairer, H. Shen, "The dynamic $\Phi^4_3$ model comes down from infinity", *Comm. Math. Phys.* 356 (2017), no. 3, 673–753.**
本纪元最重要的"全局化"定理之一：$\Phi^4_3$ 的随机动力学对**任意**（甚至无限大的）初值有全局解，且以 $\Phi^4_3$ 测度为不变测度；解在有限时间后进入紧集。技术核心是先验界 + 测度的对数 Sobolev/大偏差型控制。**未在本轮核查**（未记 arXiv 号），标 UNVERIFIED；卷/页按记忆给出，正文引用时应复核。

**11. [T] J.-C. Mourrat, H. Weber, "Global well-posedness of the dynamic $\Phi^4$ model in the plane", *Ann. Probab.* 45 (2017), no. 4, 2398–2476.**
2D 的对应全局结果，并证明 Ising–Kac 模型收敛到 $\Phi^4_2$。与条目 10 一起，构成本纪元"标量场全局动力学"的完整图景。**未在本轮核查**，标 UNVERIFIED。

**12. [T] G. Da Prato, A. Debussche, "Strong solutions to the stochastic quantization equations", *Ann. Probab.* 31 (2003), no. 4, 1900–1916.**
前 RS 时代的局部理论：把 $\Phi^4_3$ 的白噪声项吸收进线性传播子，得到等价的无奇异初值问题，证明局部强解。这是理解"RS 比 2003 年多了什么"的基准：**RS 多出的是逐点展开与重正化的系统性，而不是存在性本身**。**未在本轮核查**，标 UNVERIFIED。

**13. [T][V] M. Gubinelli, H. Koch, T. Oh, "Paracontrolled approach to the three-dimensional stochastic nonlinear wave equation with quadratic nonlinearity", *J. Eur. Math. Soc.* 26 (2024), no. 3, 817–874. arXiv:1811.07808 (v1 2018-11-19).**
3D 二次非线性随机波方程的局部适定；展示仿微分方法在**双曲**情形需要额外的色散利用与新的随机算子。属本纪元末端的预印本、下一纪元发表。作为"RS/仿微分不止抛物"的边界案例收录。

### 2.3 重整化群：BBS 与 Dimock–Balaban（5 条）

**14. [O][V] R. Bauerschmidt, D. C. Brydges, G. Slade, *Introduction to a Renormalisation Group Method*, Lecture Notes in Mathematics 2242, Springer, 2019. xii+281 pp. ISBN 978-981-32-9591-9.**
BBS 方法的教科书化：高斯积分、范数代数、扰动重正化、稳定流形定理；模型是 4D $|\varphi|^4$、Ising、弱自回避走。**这是 4D 边缘情形的严格 Wilson RG**，不是 $\varphi^4_3$。对本 campaign：它是"把 RG 写成 Banach 空间上的不动点"的最完整模板；若要把 YM 的有效作用量放进一个可控类，先要造出规范不变的范数代数（见 §6）。

**15. [T][V] R. Bauerschmidt, D. C. Brydges, G. Slade, "Scaling limits and critical behaviour of the 4-dimensional $n$-component $|\varphi|^4$ spin model", *J. Stat. Phys.* 157 (2014), 692–742. DOI 10.1007/s10955-014-1060-5. arXiv:1403.7424.**
$n\ge1$、小耦合：磁化率对数修正指数 $(n+2)/(n+8)$；比热按 $n=1,2,3$ / $n=4$ / $n>4$ 三分；环面上近临界极限是高斯自由场，非临界是白噪声。**这给出了 4D 边缘标量场的完整严格图像**：非高斯但 Gaussian 极限。对本 campaign：它是"边缘 + 正 beta ⇒ 平凡"的严格版本；YM 的 beta 为负，故不能照搬（见 §4.3）。

**16. [T][V-partial] R. Bauerschmidt, D. C. Brydges, G. Slade, "A renormalisation group method. III. Perturbative analysis", *J. Stat. Phys.* 159 (2015). arXiv:1403.7252.**
五篇系列中的第三篇：二阶扰动重正化、耦合常数流及其系数估计；用于 4D 弱自回避走的对数修正。**DOI 10.1007/s10955-014-1165-x 已核；卷/页只核到 159 (2015)，具体页码未核**。系列 I 为 arXiv:1403.7268（作者页确认标题）。作为 BBS 方法的原始技术来源引用。

**17. [O] J. Dimock, "The renormalization group according to Balaban", I–III, *Rev. Math. Phys.* 25 (2013) and 26 (2014).**
把 Balaban 1980 年代的四篇 CMP 长文重写成现代可读版本：小场/大场分解、有效作用量的稳定类、逐尺度估计。**本纪元内对 Balaban 路线最重要的二次文献**。**未在本轮核查**；卷号按记忆给出，具体 article number 未核，标 UNVERIFIED。对本 campaign：这是"4D UV 控制"的现成语言；缺口在 IR 与规范不变可观测。

**18. [O] A. Kupiainen, "Renormalization group and stochastic PDEs", *Ann. Henri Poincaré* 17 (2016), no. 3, 497–535.**
把构造性 QFT 的 RG 与 SPDE 的正则性结构并置的综述；对"为什么次临界可用、边缘不可用"给出物理侧的解释。**未在本轮核查**，标 UNVERIFIED。

### 2.4 规范场：Chatterjee 与 Shen–Zhu–Zhu（3 条）

**19. [T][V] S. Chatterjee, "Rigorous solution of strongly coupled $SO(N)$ lattice gauge theory in the large $N$ limit", *Comm. Math. Phys.* 366 (2019), no. 1, 203–268. DOI 10.1007/s00220-019-03353-3. arXiv:1502.07719 (v1 2015-02-26).**
任意维数、't Hooft 标度、$N\to\infty$：Wilson 环期望的**严格**计算，结果写成格点弦论中轨迹的绝对收敛和——一个显式的 gauge–string 对偶。技术是随机几何/腔室展开，不是 RS。**限制**：强耦合、大 $N$、格点；**不给出**连续极限、不给出 $N$ 有限、不给出质量隙（质量隙由 SZZ 另给）。对本 campaign：这是"强耦合端的精确边界条件"，可作为 RG 流的 IR 锚点。

**20. [O][V] S. Chatterjee, "Yang–Mills for probabilists", arXiv:1803.01950 (2018-03-05); in *Probability and Analysis in Interacting Physical Systems*, Springer Proc. Math. Stat. 283 (2019), 1–16.**
17 页的纲领性综述：把 Euclidean YM、格点规范场、连续极限的物理含义、以及"为什么这是概率问题"讲清楚。**不是定理**，但它是本纪元"概率学派如何看待 4D"的最清晰自述。对本 campaign：v3 的引言可以直接引用它对障碍的分类。

**21. [T][V] H. Shen, R. Zhu, X. Zhu, "A stochastic analysis approach to lattice Yang–Mills at strong coupling", *Comm. Math. Phys.* 400 (2023), no. 2, 805–851. DOI 10.1007/s00220-022-04609-1. arXiv:2204.12737 (v1 2022-04-27).**
任意 $d>1$、't Hooft 标度、强耦合：研究格点 YM 的 Langevin 动力学。条件为 SU(N) 时 $|\beta|<1/[16(d-1)]$，SO(N) 时 $|\beta|<(N-2)/[32(d-1)N]$。结论：(i) 无限格点上不变测度唯一；(ii) 动力学在 Wasserstein 距离下指数遍历；(iii) 有限体积测度收敛到该测度；(iv) 极限满足对数 Sobolev 与庞加莱不等式；(v) 由此 Wilson 环关联因子化、确定性强耦合大 $N$ 极限、广泛可观测量的关联**指数衰减——严格正的质量隙**。**$d=4$ 时 SU(N) 的门槛是 $|\beta|<1/48$**，与本项目 P3 强耦合残数定理的阈值数值相同（同一 't Hooft 归一化下；是否逐字同一参数需在 v3 中核对）。**不给出**连续极限：这是格点、有限格距。

### 2.5 越界前向指针：2D/3D 构造（2 条，按任务要求收录）

**22. [T][V] A. Chandra, I. Chevyrev, M. Hairer, H. Shen, "Langevin dynamic for the 2D Yang–Mills measure", *Publ. Math. Inst. Hautes Études Sci.* 136 (2022), 1–147. DOI 10.1007/s10240-022-00132-0. arXiv:2006.04987 (2020).**
2D YM 的随机量子化：构造 Langevin 动力学，证明全局解，并把 YM 测度识别为不变测度。2D 的 YM 测度本身已由 Driver/Sengupta/Lévy 构造（条目 25–27 的背景），所以这是"在已知测度上把动力学做出来"，不是新测度。**预印本 2020，属下一纪元；本纪元末启动**。对本 campaign：它是 CCHS 3D 方法的热身，也是"gauge-covariant in law"技术的第一次完整演示。

**23. [T][V] A. Chandra, I. Chevyrev, M. Hairer, H. Shen, "Stochastic quantisation of Yang–Mills–Higgs in 3D", *Invent. Math.* 237 (2024), no. 2, 541–696. DOI 10.1007/s00222-024-01264-2. arXiv:2201.03487 (v1 2022-01-10).**
本纪元之后、但必须在此定位的构造：3D YM–Higgs 的随机量子化。步骤：(i) 造一个分布空间上的非线性度量空间 $\mathcal{S}$ 作为初值空间；(ii) 用确定性 YMH 流的规范协变性把 gauge 等价 $\sim$ 推到 $\mathcal{S}$，得到规范轨道商 $\mathfrak{O}$；(iii) 用正则性结构证明**重整化后**的随机 YMH 流有**局部时间**解；(iv) 用"小噪声极限下的对称性"证明存在唯一一组 counterterm 使解**在分布意义下规范协变**（gauge covariant in law），从而在 $\mathfrak{O}$ 上得到典范 Markov 过程——**允许有限时间爆破**。**不证明**：全局存在、不变测度、OS 公理、质量隙。**是 YM–Higgs，不是纯 YM**；Higgs 场是状态空间与度量能定义的关键。3D 超可重整 ⇒ counterterm 有限，这是它能成立的根本原因。

### 2.6 背景锚点（8 条；本纪元之前，但正文论证需要）

**24. [O] K. G. Wilson, "Confinement of quarks", *Phys. Rev. D* 10 (1974), 2445–2459.**
格点规范场与强耦合面积律的起点。**未在本轮核查**，标 UNVERIFIED；卷/页按标准引用给出。

**25. [T] K. Osterwalder, E. Seiler, "Gauge field theories on a lattice", *Ann. Physics* 110 (1978), 440–471.**
强耦合格点 YM 的 OS 公理与质量隙（以格点距离表述）。本 campaign 的"强耦合有隙"最原始的严格来源。**未在本轮核查**，标 UNVERIFIED。

**26. [O] E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, Lecture Notes in Physics 159, Springer, 1982.**
把格点规范场放进构造性 QFT 框架的专著；强耦合展开、团簇展开、OS 公理的系统处理。**未在本轮核查**，标 UNVERIFIED。

**27. [C] T. Balaban, "Renormalization group approach to lattice gauge field theories" I–IV, *Comm. Math. Phys.* 109 (1987); 115–119 (1988).**
4D 格点 YM 的 UV 稳定性：小场鞍点展开 + 大场稀有性，逐尺度控制有效作用量。**证到**：配分函数在 $a\to0$ 时不爆炸、沿渐近自由轨迹重整化。**未证**：规范不变局部关联、无限体积、OS、质量隙。**未在本轮核查**；卷号按标准引用给出，具体篇目边界未核，标 UNVERIFIED。对本 campaign：这是"4D UV 可控"的唯一完整方案；v3 若走"RS × Balaban"混合，必须先精确复述它的有效作用量类。

**28. [C] J. Magnen, V. Rivasseau, R. Sénéor, "Construction of $YM_4$ with an infrared cutoff", *Comm. Math. Phys.* 155 (1993), 417–430.**
用多尺度/相空间展开构造带 IR 截断的 4D YM（规范固定，SU(2)）。**是 4D 唯一被写下来的连续统式构造**，但：规范固定 ⇒ 不是物理理论；有 IR 截断 ⇒ 不是无限体积；无质量隙。**未在本轮核查**，标 UNVERIFIED。对本 campaign：它证明"4D 的 UV 部分不是死路"，死路在 IR 与规范不变性。

**29. [O] G. Parisi, Y.-S. Wu, "Perturbation theory without gauge fixing", *Sci. Sinica* 24 (1981), 483–496.**
随机量子化的原始文献。**未在本轮核查**，标 UNVERIFIED。仅作历史锚点。

**30. [O] A. Jaffe, E. Witten, "Quantum Yang–Mills Theory", Clay Mathematics Institute Millennium Prize problem description (2000).**
问题的官方表述：任意紧单群、$\mathbb{R}^4$ 上非平凡量子 YM 的存在性 + 质量隙 $\Delta>0$，公理强度至少到 Streater–Wightman / Osterwalder–Schrader。**未在本轮核查**（Clay 官网原文），标 UNVERIFIED；表述本身在项目内已多处引用。

**31. [NG] M. Aizenman, H. Duminil-Copin, "Marginal triviality of the scaling limits of critical 4D Ising and $\varphi^4$ models", *Ann. of Math. (2)* 194 (2021), no. 1, 163–235. arXiv:1912.07973.**
4D 边缘标量场的缩放极限是 Gaussian：$\varphi^4_4$ 平凡。**对本 campaign 的关键区分**：这是**标量**的 no-go，不能推到 YM。YM 也是边缘的，但 beta 函数为负（渐近自由），UV 极限仍是 Gaussian 却允许非平凡的红外；$\varphi^4_4$ 的 beta 为正，边缘性直接把理论推向平凡。**未在本轮核查**，标 UNVERIFIED；卷/页按标准引用给出。此条属下一纪元，作为"边缘性的两种命运"的对照收录。

---

## 3. $d=2,3,4$ 的严格状态与耦合区间

| 维数 | 对象 | 严格状态 | 耦合/参数范围 | 质量隙 | 来源 |
|---|---|---|---|---|---|
| $d=2$ | 纯 YM 测度 | 构造完成；OS 公理；面积律 | 紧曲面，全 $\beta>0$ | 退化（无传播自由度） | Driver 1989; Sengupta 1997; Lévy 2003（背景） |
| $d=2$ | YM Langevin | 全局解；不变测度 = YM 测度 | 无小耦合限制 | — | CCHS 2022（条目 22） |
| $d=3$ | $\Phi^4_3$ | 全局 Markov 过程；不变测度；任意初值 | 构造：所有耦合；质量隙的严格证明：弱耦合 | 有（弱耦合，1976 Feldman–Osterwalder） | 条目 10, 11 |
| $d=3$ | YM–Higgs（随机量子化） | **局部时间**重整化解；gauge covariant in law；Markov（可能有限时间爆破） | 非微扰；无小耦合限制 | **无** | CCHS 2024（条目 23） |
| $d=3$ | 纯 YM 连续统 | **无构造** | — | 仅格点强耦合有 | — |
| $d=3$ | 格点 YM 强耦合 | 唯一无限体积测度；指数遍历；LSI/庞加莱 | SU(N): $|\beta|<1/32$ | **有（严格）** | SZZ 2023（条目 21） |
| $d=4$ | 纯 YM 连续统 | **无构造**；只有 UV 稳定性 / 带 IR 截断 / 格点 | — | 无 | Balaban; MRG |
| $d=4$ | 格点 YM 强耦合 | 唯一无限体积测度；指数遍历；LSI | SU(N): $|\beta|<1/48$ | **有（严格）** | SZZ 2023（条目 21） |
| $d=4$ | $\varphi^4_4$ | 平凡（Gaussian 缩放极限） | 所有耦合 | 平凡 | AD-C 2021（条目 31） |
| 任意 $d$ | 格点 YM 大 $N$ | Wilson 环严格解；显式 gauge–string 对偶 | 't Hooft 标度；$N\to\infty$；强耦合 | 由 SZZ 另行给出 | Chatterjee 2019（条目 19） |

**读表要点**：
- $d=3$ 的"YM 有构造"是**误读风险最高**的一句。CCHS 2024 是 **YM–Higgs**、**局部时间**、**允许爆破**、**无质量隙**；纯 3D YM 仍无连续统构造。
- $d=3$ 的超可重整性（$[g^2]=1$）使 counterterm 有限，这是 CCHS 能定义状态空间与重正化流的根本；$d=4$ 的 $[g^2]=0$ 使同一套机器失效。
- SZZ 的 $d=4$ 门槛 $1/48$ 与项目 P3 的 $\beta<1/48$ **数值相同**；是否为同一归一化下的同一参数，v3 必须逐字核对（见 §6 第 10 条）。

---

## 4. 为什么 4D 规范场比 3D 难

### 4.1 维数计数：超可重整 vs 边缘

耦合的量纲：$[g^2]=4-d$。$d=2$：$[g^2]=2$；$d=3$：$[g^2]=1$；$d=4$：$[g^2]=0$。

- $d\le3$：**超可重整**。发散只出现在有限个图/有限个 counterterm 上；把质量、真空能、波函数重正化掉之后，方程在白噪声下是"良定义"的。
- $d=4$：**边缘**。每一阶都有新的 counterterm，总数无穷；没有任何"先减掉有限个坏项"的操作能让方程变成经典意义下的适定方程。

这正是 RS 的门槛：Hairer 的理论要求方程**次临界**（非线性项的度为负），从而只有有限多个"负度"的树需要重正化。4D YM 的非线性项度数为 0（临界），树空间在度 $\le 0$ 的部分**无限**；RS 的逐树估计逐项不可和。**这不是"还没人做"，是框架的适用条件不满足。**

### 4.2 次临界门槛的精确含义

RS 的收敛证明依赖：对每个度 $\le 0$ 的树，其"期望"是有限的；而度 $<0$ 的树只有有限多个。次临界性保证了这个"只有有限多个"。

- 3D YM：$[g^2]=1$，非线性项度为 $-1/2$ 类；有限多个 counterterm（如 $A_0$ 的质量项、Higgs 势的重正化）。
- 4D YM：$[g^2]=0$，非线性项度为 0；**每一个圈数都产生新的、度相同的 counterterm**。代数 BPHZ（条目 4）仍能逐阶定义，但"逐阶"不等于"求和"；解析 BPHZ（条目 5）的收敛性在临界情形没有对应定理。

### 4.3 渐近自由的双刃

4D YM 的 beta 函数为负：$\beta(g)=-\beta_0g^3+O(g^5)$，$\beta_0>0$。后果是双面的：

- **好的一面**：UV 固定点是 Gaussian，连续极限存在且渐近自由；这与 $\varphi^4_4$（beta 为正、边缘 ⇒ 平凡，条目 31）形成鲜明对照。**"4D 边缘 ⇒ 平凡"不能搬到 YM。**
- **坏的一面**：微扰级数是**渐近**而非收敛的；running coupling 在 IR 增长，微扰论在 $\Lambda_{\rm QCD}$ 以下失效。要构造测度，必须控制从 UV 到 IR 的**全部尺度**，而不是某一阶。
- 对 RG 的具体含义：BBS 的 Banach 不动点（条目 14–16）在小耦合下工作；YM 的小耦合只在 UV。**把 UV 的 Gaussian 不动点沿负 beta 流搬到 IR 的强耦合区，正是本项目所说的"交接"**；本纪元的文献没有这台机器。

### 4.4 大场与规范轨道

- **大场**：RS/仿微分的解是**局部时间**的，估计在"场有多大"上是指数坏的；$\Phi^4_3$ 能全局化，靠的是不变测度给出的先验界（条目 10）。4D YM 没有已知的不变测度可供先验界。
- **规范轨道**：CCHS 的关键创新是把"分布空间除以规范群"做成非线性度量空间 $\mathfrak{O}$（条目 23）。3D 能做，是因为超可重整性让 YMH 流在分布层面可控；4D 的 Gribov 问题（每个规范轨道与任何全局截面的交点不唯一）使"先在规范固定下构造、再证明规范不变"这条老路走不通，而"直接在轨道空间上构造"又缺少 4D 的度量/紧性。
- **可观测代数**：物理可观测是规范不变的 Wilson 环/胶球算子。即使造出某个规范固定的测度，也未必下降为轨道空间上的测度；这是 MRG（条目 28）不能算"物理构造"的原因。

### 4.5 质量隙

质量隙要求**无限体积、体积一致**的指数聚类。本纪元的所有 3D 构造（条目 10, 11, 23）都**不给出**质量隙：
- $\Phi^4_3$ 的质量隙有 1976 年的弱耦合证明（背景），但不是从 RS 构造里读出的；
- CCHS 3D 是局部时间，连不变测度都没有，谈不上谱；
- 唯一严格的 4D 质量隙来自**格点强耦合**（SZZ，$\beta<1/48$），不是连续统。

因此"3D 有构造 ⇒ 4D 只差技术"是错误推论：3D 构造本身就没有质量隙，4D 连构造都没有。**两个缺口是独立的。**

---

## 5. 在 14 个纪元中的位置

> 项目内没有既成的纪元划分文件；以下是本文档首次钉死的**工作分期**（主题分期，边界可重叠）。本纪元为第 **11** 纪元，与文件名 `LIT-ERA-11` 对应。

1. **1954–1972 前史**：非阿贝尔规范场、Faddeev–Popov、Gribov 歧义、正则量子化；无严格构造。
2. **1973–1976 渐近自由与公理化**：Gross–Wilczek/Politzer；Osterwalder–Schrader 公理；Wilson 格点规范场（1974）。
3. **1975–1981 强耦合格点**：Wilson 面积律、Osterwalder–Seiler 质量隙、Seiler 公理化；2D 精确解（Migdal、Gross–Witten）。
4. **1978–1984 2D/3D 构造性 QFT**：Glimm–Jaffe–Spencer $P(\varphi)_2$；Feldman–Osterwalder $(\varphi^4)_3$ 弱耦合质量隙；Brydges–Fröhlich–Seiler 团簇展开。
5. **1982–1989 Balaban RG**：4D 格点 YM 的 UV 稳定性；小场/大场分解；有效作用量的稳定类。
6. **1990–1997 多尺度 4D**：Magnen–Rivasseau–Sénéor 带 IR 截断 $YM_4$；Gawedzki–Kupiainen 4D $\varphi^4$；严格重整化的多尺度语言。
7. **1995–2003 2D 测度完成**：Driver、Sengupta、Lévy 的 2D YM 测度；Gribov 量化的数学化；$\varphi^4_4$ 平凡性的早期讨论。
8. **2003–2008 随机量子化复兴**：Da Prato–Debussche 局部强解；Hairer–Mattingly 遍历性与谱隙；rough paths 的成熟。
9. **2009–2012 奇异 SPDE 前夜**：Hairer 的 KPZ 纲领；Chatterjee 格点规范场纲领启动；强耦合面积律的严格化。
10. **2012–2013 突破**：KPZ 解出（Hairer, *Ann. Math.* 2013）；仿微分分布（GIP）；正则性结构预印本（arXiv:1303.5113）。
11. **2013–2018 本纪元：严格概率与构造性 QFT**——RS 理论成型、代数/解析 BPHZ、$\Phi^4_3$ 全局化、BBS 4D RG、Chatterjee 大 $N$、Dimock–Balaban、3D YM 随机量子化启动；**4D 仍为零**。
12. **2019–2021 成熟与全局动力学**：BBS 专著；Aizenman–Duminil-Copin $\varphi^4_4$ 平凡性；Hairer–Shen 全局 $\Phi^4_3$ 的完整写法；2D YM Langevin（CCHS）。
13. **2022–2024 3D 规范构造**：CCHS 3D YM–Higgs（*Invent. Math.* 2024）；SZZ 强耦合质量隙（*CMP* 2023）；Cao–Park–Sheffield 随机曲面面积律。
14. **2025–2026 形式化与 AI**：Lean 4 的 4D 条件构造、OSforGFF、面积律扩张、paper v3 的综合与审计。

**相邻关系**：第 10 纪元给本纪元提供语法（KPZ + 仿微分）；第 12 纪元把本纪元的局部结果全局化；第 13 纪元把本纪元的规范场技术推到 3D YM–Higgs。**本纪元是"机器造好、但还没碰到 4D"的纪元。**

---

## 6. v3 的缺口与攻击面

1. **4D 没有任何概率构造，连条件构造都没有。** CCHS 的轨道空间方法停在 3D 局部时间；RS 停在次临界。攻击面：把 CCHS 的 $\mathfrak{O}$ 换成**尺度依赖的范数**（running coupling 进范数），使每一尺度的 counterterm 只在该尺度可见；但边缘性使"固定次数展开"失效，必须先解决 §4.2。
2. **无穷 counterterm 的求和。** 代数 BPHZ 逐阶有限，总和发散。攻击面：按 running coupling 重排 counterterm 族，证明其渐近展开可 Borel 求和；或证明 Ward/BRST 把 gauge-covariant counterterm 空间压到 1 维（若是，则"无穷多"只是同一个参数的无穷次重正化）。**这是最可能出新定理的一条。**
3. **大场控制。** RS 是局部的；Balaban/BBS 的大场分解是全局的。攻击面：造 **RS × Balaban 混合**——小场用 RS 的逐点展开，大场用 polymer/范数代数；接口处需要"场大小的阈值"与"规范不变性"兼容。
4. **规范协变性 in law 的 4D 分类。** 3D 的 counterterm 选择由对称性唯一确定（CCHS 的 symmetry argument）；4D 无穷多。攻击面：用 BRST 上同调 + OS 反射正性证明不变锥 1 维；若锥的维数 >1，则存在多个"规范协变"的连续极限，本身就是新现象。
5. **质量隙与构造解耦。** 3D 构造不给隙，4D 强耦合给隙但不给连续极限。攻击面：把 SZZ 的**对数 Sobolev 不等式沿 RG 流一致化**——即证明 LSI 常数在 UV 截断下一致。这与项目的 $G''(0)>0$ / 粗模谱隙下界是同一件事的两种语言，v3 应把等价性写成定理。
6. **LSI 作为搬运工具。** SZZ 的 LSI ⇒ 指数聚类 ⇒ 质量隙；但 LSI 是对**格点测度**证的。攻击面：证明 LSI 在 $a\to0$ 下的一致性是**充分**的（比逐阶 BPHZ 更强），从而把"连续极限 + 质量隙"合并成一条估计。
7. **BBS 的范数代数只对 $\varphi^4$/SAW。** 攻击面：构造**规范不变的范数代数**（Wilson 环的代数、或全息 Wilson 流的代数），把 BBS 的 Banach 不动点搬过来；项目的 gauge W-熵是候选。
8. **3D 纯 YM 仍缺。** CCHS 含 Higgs；Higgs 是状态空间与度量的关键。攻击面：Higgs 质量 $\to\infty$ 的极限，或直接在规范不变可观测代数上取极限；先例是 2D 的 Driver–Sengupta 测度论。
9. **大 $N$ 只到强耦合。** Chatterjee 的弦求和是 IR 端的精确边界条件。攻击面：把它当作 RG 流的**红外锚点**（而不是全程解），用 UV 的渐近自由 + IR 的弦求和夹逼中间尺度；这是"两条已知端点的插值"，比从头造容易。
10. **验证缺口（工程）。** SZZ 的 $\beta<1/48$ 与项目 P3 的 $\beta<1/48$ 是否同一归一化、同一 $\beta$ 定义，必须逐字核对；SZZ 的常数（$1/[16(d-1)]$）没有独立复算。攻击面：在 Lean 中重推该常数，并把"两个 $1/48$ 是同一个"写成显式引理。**这是 v3 最便宜、最可验收的一步。**

---

## 7. 核查台账

| # | 条目 | 核查方式（2026-11-24） | 结论 | 残留不确定 |
|---|---|---|---|---|
| 1 | Hairer 2014, *Invent. Math.* 198, 269–504 | WebSearch 元数据 | **VERIFIED** | 无 |
| 2 | Hairer 2015, *Braz. J. Probab. Stat.* 29, 175–210 | arXiv:1401.3014 作者页 | **VERIFIED（arXiv）** | 期刊卷/页未核 |
| 3 | Hairer 2013, *Ann. of Math.* 178, 559–664 | arXiv:1109.6811 作者页 | **VERIFIED（arXiv）** | 期刊卷/页未核 |
| 4 | Bruned–Hairer–Zambotti, *Invent. Math.* 216, 1039–1156 | arXiv:1610.08468 摘要页 | **部分** | 只核到 "to appear in Invent. Math."；卷/页 UNVERIFIED |
| 5 | Chandra–Hairer, arXiv:1612.08138 | arXiv 摘要页 + 作者主页 | **VERIFIED（未发表）** | 无；"未发表"本身已核 |
| 6 | Hairer, *JEMS* 23, 869–947 | arXiv:1711.10239 作者页 | **VERIFIED（arXiv）** | 期刊卷/页未核 |
| 7 | Hairer, arXiv:1508.05261 | 作者页 | **VERIFIED（arXiv）** | 期刊页未核 |
| 8 | Gubinelli–Imkeller–Perkowski, *Forum Math. Pi* 3, e6 | 无（本轮 WebSearch 耗尽） | **UNVERIFIED** | 卷/页/年份 |
| 9 | Catellier–Chouk, *Ann. Probab.* 46, 2621–2679 | WebSearch + DOI 10.1214/17-AOP1235 | **VERIFIED** | 无 |
| 10 | Hairer–Shen, *CMP* 356, 673–753 | 无 | **UNVERIFIED** | arXiv 号、卷/页 |
| 11 | Mourrat–Weber, *Ann. Probab.* 45, 2398–2476 | 无 | **UNVERIFIED** | arXiv 号、卷/页 |
| 12 | Da Prato–Debussche, *Ann. Probab.* 31, 1900–1916 | 无 | **UNVERIFIED** | 卷/页 |
| 13 | Gubinelli–Koch–Oh, *JEMS* 26, 817–874 | WebSearch（含 zbMATH/JEMS） | **VERIFIED** | 无 |
| 14 | BBS 专著, LNM 2242 | WebSearch（Springer/Stanford/MR3969983） | **VERIFIED** | 无 |
| 15 | BBS 2014, *JSP* 157, 692–742 | arXiv:1403.7424 摘要页 | **VERIFIED** | 无 |
| 16 | BBS 2015, *JSP* 159 | arXiv:1403.7252 摘要页 + DOI | **部分** | 具体页码未核 |
| 17 | Dimock, *Rev. Math. Phys.* 25/26 | 无 | **UNVERIFIED** | 卷/篇号 |
| 18 | Kupiainen, *Ann. Henri Poincaré* 17, 497–535 | 无 | **UNVERIFIED** | 卷/页 |
| 19 | Chatterjee, *CMP* 366, 203–268 | Crossref API | **VERIFIED** | 无 |
| 20 | Chatterjee, Springer Proc. 283, 1–16 | WebSearch + arXiv:1803.01950 | **VERIFIED** | 无 |
| 21 | Shen–Zhu–Zhu, *CMP* 400, 805–851 | arXiv:2204.12737 摘要页 + WebSearch | **VERIFIED** | 无（并更正简报的 2020） |
| 22 | CCHS 2022, *Publ. IHÉS* 136, 1–147 | WebSearch（Springer/Imperial） | **VERIFIED** | 无 |
| 23 | CCHS 2024, *Invent. Math.* 237, 541–696 | WebSearch + arXiv:2201.03487 | **VERIFIED** | 无 |
| 24 | Wilson 1974, *PRD* 10, 2445 | 无 | **UNVERIFIED** | 页码 |
| 25 | Osterwalder–Seiler 1978, *Ann. Physics* 110, 440 | 无 | **UNVERIFIED** | 页码 |
| 26 | Seiler 1982, LNP 159 | 无 | **UNVERIFIED** | ISBN/页数 |
| 27 | Balaban, *CMP* 109/115–119 | 无 | **UNVERIFIED** | 篇目边界 |
| 28 | MRG 1993, *CMP* 155, 417–430 | 无 | **UNVERIFIED** | 页码 |
| 29 | Parisi–Wu 1981, *Sci. Sinica* 24, 483 | 无 | **UNVERIFIED** | 页码 |
| 30 | Jaffe–Witten, Clay 2000 | 无 | **UNVERIFIED** | 无（文本在项目内多处转引） |
| 31 | Aizenman–Duminil-Copin 2021, *Ann. of Math.* 194, 163–235 | 无 | **UNVERIFIED** | 卷/页 |

**统计**：31 条；**VERIFIED 或部分 VERIFIED：17 条**（#1–7, 9, 13–16, 19–23）；**UNVERIFIED：14 条**（#8, 10–12, 17–18, 24–31）。UNVERIFIED 条目中，背景锚点（#24–30）的卷/页均为标准引用，但在 v3 定稿前必须逐条复核；#10–12 是本纪元的承重数学条目，优先级最高。

**复核优先级（建议 v3 前完成）**：
1. #10 Hairer–Shen、#11 Mourrat–Weber — 全局 $\Phi^4$ 的准确陈述与 arXiv 号；
2. #8 GIP — 仿微分的准确出处；
3. #17 Dimock — Balaban 三篇的篇号；
4. #27 Balaban 原文 — 四篇的边界与每篇结论；
5. #31 AD-C — 卷/页与"平凡性"的精确表述（避免把标量 no-go 误推到 YM）。

---

## 附：一页速查

```
本纪元 = 机器造好，但没碰到 4D
  RS(2014) + 代数 BPHZ + 解析 BPHZ → 次临界奇异 SPDE 的局部理论
  Φ⁴₃: 局部(Catellier–Chouk / Hairer) → 全局(Hairer–Shen)
  4D |φ|⁴: BBS 的严格 Wilson RG（小耦合、Gaussian 极限）
  规范场: 2D Langevin(CCHS 2022) / 3D YM-Higgs 局部(CCHS 2024)
  强耦合质量隙: SZZ 2023, SU(N) d=4 时 β<1/48
  4D 纯 YM: 无构造、无 OS、无质量隙
为什么 4D 难 = 边缘([g²]=0) + RS 要次临界 + 无穷 counterterm
            + 渐近自由使微扰发散 + 大场/Gribov + 隙与构造解耦
最便宜的可验收一步 = 核对两个 1/48 是否同一 β（并在 Lean 中重推）
```

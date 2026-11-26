# L5-CLASSICAL · 严格构造性 4D 规范理论这一线：有没有带 `a`-依赖的质量隙陈述？

任务：逐篇核摘要/正文，回答唯一问题——**这条线上是否有任何一篇给出 `m` 的界（哪怕只在某范围内），且该界的 `a`-依赖被明确写出？**

纪律：标题/摘要逐字；拿不到正文/摘要的，写"未取数（检索式：…）"，不凭记忆补。
取数方式：期刊/档案馆摘要页（Springer 镜像 NII-REO、CERN CDS、arXiv abs）为主；PDF 正文在本会话的抓取器里返回 `unsupported content type "application/pdf"`，因此凡标【snippet】者均为检索索引对 PDF 的逐字切片（附 URL 与锚点），非我本人读取正文所得。所有逐字引文均标来源。

---

## 0. 一句话结论

**（b）成立**：本线（Balaban / Magnen–Rivasseau–Sénéor / Dimock–Hurd / 格点红外）中，凡我取到摘要的每一篇，最强结论都是 **UV 稳定性 / 有效作用量构造 / 变分问题极小元存在 / 可观测代数构造 / 非零温红外界**——**没有一篇给出 `m`（质量隙）的界，遑论带 `a`-依赖的界**。
本线唯一真正的 `m` 界来自 **Osterwalder–Seiler 1978**（强耦合），而它的形式是 **`m·a ≥ c(β)`（格点单位下 `m ≥ c > 0`）**，即物理单位 `m ≥ c/a → ∞`——这是"强耦合相关长度 = O(a)"的平凡 UV 陈述，**不构成任何 `a→0` 一致的正质量隙**。这正是这条线没有入口的**结构性原因**：它的机器全部架在"UV 稳定 / 红外冻结"上，而质量隙是纯红外对象。

---

## 1. Dimock–Hurd（CMP，1980s）："observable algebra" / "scaling limits"

- **已核实存在的一篇（同作者、同年代簇）**：J. Dimock, T. Hurd, *A renormalization group analysis of the Kosterlitz-Thouless phase*, **Commun. Math. Phys. 137 (1991) 263–287**。出处：Dimock 本人论文目录（[publications.pdf](https://www.nsm.buffalo.edu/~dimock/publications.pdf)，本会话 PDF 不可读，仅由检索索引给出该行）。
- **"Construction of the observable algebra for Yang-Mills theory" 一文**：**未取数**。检索式：`au:Dimock au:Hurd "observable algebra" Yang-Mills`（zbMATH / MathSciNet / Springer Link 站内）。三次不同措辞检索（`"observable algebra for Yang-Mills" Dimock Hurd`、`Dimock Hurd "Yang-Mills" "Letters in Mathematical Physics" 1999`）均**未返回该标题的书目记录**，只返回 Dimock 1980 的 *Algebras of local observables on a manifold*（CMP 77, 219–228）与若干二手引用。⇒ 我不写标题，不写年份，不写期刊。
- 可确证的 Dimock 侧最近代的相关工作：*The renormalization group according to Balaban. I: Small fields*, **Rev. Math. Phys. 25 (2013) 1330010**（二手引用一致，未见摘要原文）。
- **是否有 `m` 的界**：否（KT 相 / dipole gas / Balaban 综述均不以 `m` 为结论对象）。
- **`a`-依赖**：不适用。

## 2. Magnen–Rivasseau–Sénéor, *CMP* **155** (1993) 325–383

- 作者/年/刊：J. Magnen, V. Rivasseau, R. Sénéor, **Commun. Math. Phys. 155** (1993), Issue 2, 325–383；**DOI 10.1007/BF02097397**。来源：[NII-REO 摘要页](https://reo.nii.ac.jp/oja/ART1000345241/en)。
- **标题逐字**：`Construction of YM4 with an infrared cutoff`（页面排版作 "Construction ofYM4"，为排版吞空格）。
- **最强结论（摘要逐字）**：
  > "We provide the basis for a rigorous construction of the Schwinger functions of the pure SU(2) Yang-Mills field theory in four dimensions (in the trivial topological sector) **with a fixed infrared cutoff but no ultraviolet cutoff**, in a regularized axial gauge. … We prove that this combined effect is **stabilizing** if we use cutoffs of a certain type in momentum space. We check the validity of the construction by showing that **Slavnov identities** … do hold non-perturbatively."
- **判读**：证的正是**UV 稳定性 + Slavnov 恒等式（非微扰规范不变性）**，**红外被固定为 cutoff**。质量隙是 `|x−y|→∞`（纯红外）现象，而此文把红外**切掉**了 ⇒ 结构上不可能含 `m` 界。
- **是否有 `m` 的界**：**无**。全文摘要无 `mass`、`gap`、`exponential decay` 任一语词。
- **`a`-依赖**：不适用（连续谱、动量空间 cutoff；作者明言无 UV cutoff ⇒ 根本不存在格点间距 `a`）。

## 3. Balaban, *CMP* 102 / 109 / 116 / 118 / 122 系列——逐篇核标题+摘要

| # | 出处 | 标题逐字（已核） | 摘要里的最强结论 | 有 `m` 界？ |
|---|---|---|---|---|
| 1 | **CMP 102** (1985) 277–309, DOI 10.1007/BF01229381 | `The variational problem and background fields in the renormalization group method for lattice gauge theories` | "We consider the action of a lattice gauge theory on a space of regular gauge field configurations with fixed averages, and we prove that **there exists a minimum of this action**."（NII-REO 摘要型页；CERN CDS record 161548 给卷页） | 否 |
| 2 | **CMP 102** (1985) 255–（同卷 Issue 2, p.175–347） | `Ultraviolet stability of three-dimensional lattice pure gauge field theories` | "We prove the **ultraviolet stability** for three-dimensional lattice gauge field theories. We consider only the Wilson lattice approximation for pure Yang-Mills field theories."（[NII-REO 卷目](https://reo.nii.ac.jp/oja/do/vol_issues?vol_id=ISS0010021541&lang=en)）注意：**3 维**，非 4 维 | 否 |
| 3 | **CMP 109** (1987) 249–301 | `Renormalization group approach to lattice gauge field theories: I. Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions` | "We study four-dimensional pure gauge field theories by the renormalization group approach. **The analysis is restricted to small field approximation.** In this region we construct a sequence of localized effective actions by cluster expansions in one step renormalization transformations."（[NII-REO 卷目 p.249](https://reo.nii.ac.jp/oja/vol_issues?vol_id=ISS0010021580&lang=en)） | 否 |
| 4 | **CMP 116** (1988) 1–22 | `Renormalization group approach to lattice gauge field theories: II. Cluster expansions` | "The fluctuation field integral, constructed in Part I, is represented by the exponentiated cluster expansion. It is proved that **the terms of the expansion satisfy the inductive assumptions**."（[ADS 1988CMaPh.116....1B](https://ui.adsabs.harvard.edu/abs/1988CMaPh.116....1B/exportcitation)；NII-REO vol_id ISS0010021607 同） | 否 |
| 5 | **CMP 118** (1988)、**CMP 122** (1989) 的 III / IV 篇 | **未取数** | — | — |

- 对 #5：检索式 `Balaban "Renormalization group approach to lattice gauge field theories" III 118`、`Balaban IV "Correlation functions" 122` 在本会话窗口内未返回可核的标题/摘要页。**不填记忆里的标题。**
- **模式判读（对已核 4 篇）**：Balaban 系列的结论形状恒为 **"small field 下的有效作用量 / cluster expansion 满足归纳假设 / UV 稳定性"**。**摘要中 `mass`、`gap`、`exponential decay` 均不出现。**
- **是否有 `m` 的界**：**无（已核 4 篇）**。**`a`-依赖**：不适用（结论是 `a`-一致的 UV 控制，不是 `m`）。

## 4. Tomita / King / Borgs——格点规范的红外行为

- **Borgs（+Seiler）**：C. Borgs, E. Seiler, `Lattice Yang-Mills theory at nonzero temperature and the confinement problem`, **CMP 91 (1983) 329–380**（[Google Scholar 记录](https://scholar.google.com.hk/citations?view_op=view_citation&hl=ko&user=WX0pI3wAAAAJ&citation_for_view=WX0pI3wAAAAJ:IjCSPb-OGe4C)）。摘要逐字（索引切片）：
  > "This same analogy also allows us to establish **infrared bounds for the gauge models that yield a lower bound for this critical coupling** and thereby show the existence of a **weak coupling regime without confinement at nonzero temperature in three or more space dimensions**."
  ⇒ 红外界是**关于临界耦合 `β_c` 的**（`β_c` 的下界），**不是关于 `m` 的**；且对象是**非零温**（有限温度 = 有限时间周长，质量隙在那里被温度尺度替代）。
- **Tomita**：**未取数**。检索式：`au:Tomita lattice gauge "infrared"`（CMP 1980–1990 卷目）、`"strong coupling limit" lattice gauge Tomita`。本会话检索被大量非同行评议的 zenodo 预印本淹没，未见可核书目记录 ⇒ 不写标题。
- **King**：**未取数**。检索式：`au:King C "strong coupling" lattice gauge CMP`、`King "mass gap" lattice gauge 1980s`。同上，未取到可核记录。
- **是否有 `m` 的界**：Borgs–Seiler = 无（是 `β_c` 界）；Tomita/King = 未取数。

## 5. Seiler, *Gauge Theories as a Problem of Constructive QFT*（书）

- 书目（已核）：E. Seiler, **Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics**, *Lecture Notes in Physics* **159**, Springer-Verlag, 1982, iv + 192 p.（[Chicago 目录](https://catalog.lib.uchicago.edu/vufind/Record/11067557/Similar#tabnav)、[CiNii BA03458450](https://ci.nii.ac.jp/ncid/BA03458450/refworks)）。**ISBN 未取数**（Springer 书页未取到；检索式：`site:link.springer.com "Gauge Theories as a Problem of Constructive Quantum Field Theory"`）。
- **书里对 4D YM 质量隙状态的逐字表述**：**未取数**。检索式：`Seiler LNP 159 "mass gap"`+Google Books 全文；`"Gauge theories as a problem" Seiler mass gap open`。我只取到二手旁证（见 §6 的 nLab/讲义），**不能当作书内逐字**。
- 可直接确证的相关硬结果在**同作者 1978 年论文**：K. Osterwalder, E. Seiler, `Gauge field theories on a lattice`, **Ann. Phys. 110 (1978) 440–471**, DOI 10.1016/0003-4916(78)90039-8（该 DOI 在 Balaban 1985 参考文献表中逐字出现）。**强耦合质量隙/指数衰减即出自此文**——但它的 `m` 是**格点单位下的常数**：`m a ≥ c(β) > 0`，`β` 小。**`a→0` 时 `m ≥ c/a → ∞`**：这是"相关长度 = O(a)"的平凡陈述，**不是** `a`-一致正隙。
- **是否有 `m` 的界**：书 = 未取数；书内所据的 Osterwalder–Seiler = 有，但形式为 `m a ≥ c(β)`（无物理 `a`-一致内容）。

## 6. Chatterjee, *Gauge Theory and Probability*（书/讲义）

- **书本身正文**：**未取数**（本会话未取到书页/PDF；检索式：`"Gauge Theory and Probability" Chatterjee book pdf mass gap`、`souravchatterjee.su.domains gauge`）。下列为**同一作者、同一论题的公开讲义/课程笔记**的逐字切片（【snippet】，PDF 正文不可读，故按索引切片引用并给锚点）：
  1. **讲义（Columbia CFT18）**：`Euclidean Yang-Mills theories are scaling limits of lattice gauge theories (probability-theoretic **open problem**).`【snippet，[Chatterjee.pdf](https://www.math.columbia.edu/department/probability/CFT18/Chatterjee.pdf)】
     同文另列：`Open problem (Area law) Take any compact non-Abelian Lie group G ⊆ U(N) … consider any infinite volume limit of 4D lattice gauge theory …`【snippet，同上】
  2. **报告（Utah FPD18）**：`Even physicists do not think they have a proof of mass gap.`【snippet，[Chatterjee.pdf](https://math.utah.edu/~firas/FPD18/Talks/Chatterjee.pdf)】同页：`There is a proof at small β by Osterwalder & Seiler (1978).`（对"面积律/质量隙"的现状定位：只有小 `β` 证明。）
  3. **二手（nLab 条目，转述其立场）**：`If one can show that this theory has a mass gap at all values of the coupling strength, that would explain why … glueballs … have mass. **All such questions remain open.**`（[nLab: Sourav Chatterjee](https://ncatlab.org/nlab/show/Sourav+Chatterjee/)）——二手，不作书内引文。
  4. **课程笔记（Stanford stats369）**：`This is also an open problem – both this and mass gap do not have convincing physics arguments beyond numerical simulations in four dimensions.`【snippet，[stats369.pdf](https://web.stanford.edu/~lindrew/stats369.pdf)】——二手教学材料，注明为转述。
- **是否有 `m` 的界**：**无**（其可核实结论全是"开问题"定位 + 引用 1978 强耦合结果）。

---

## 判死线裁决

**(a) 是否有任一篇给出 `m` 的界？**
本线**已核**各篇（MRS 1993；Balaban 102/109/116 及 102-3D；Borgs–Seiler 1983）：**没有任何 `m` 的界**，因此**没有可引的式号**。
唯一真实存在的 `m`-型断言在**本线之外一篇、却被全线反复引用为"唯一已知"**：Osterwalder–Seiler 1978，形式 `m a ≥ c(β)`，`β` 小。**它不是入口**——因为 `a→0` 时它给 `m ≥ c/a → ∞`，即"强耦合没有连续极限"，与"连续极限下 `m > 0`"是**相反方向**的陈述。任何以 `1/a` 为下界的"隙"在 `a→0` 极限里都自动成立而毫无内容；要的是**`a`-一致（且在 `β→β_c` 上一致）的正下界**，本线**零篇**给出。

**(b) 是否成立？** **成立。**
> **这条线上没有任何质量隙界。**

补强（为什么是结构性的，不是"还没做"）：MRS 把**红外切掉**（`with a fixed infrared cutoff`），质量隙住在被切掉的那一侧；Balaban 全系列的结论形状是 **`a`-一致的 UV 控制**（有效作用量 + cluster expansion 归纳假设 + UV 稳定性），其"小场近似"限制恰好排除了长程涨落；Borgs–Seiler 的红外界是**关于 `β_c` 的界**、且在**非零温**。三件事拼起来：这条线的机器在原理上不测量 `m`。

**(c) 领域自己的判词（逐字）**
- `Even physicists do not think they have a proof of mass gap.` —— Chatterjee, FPD18 报告，[math.utah.edu/~firas/FPD18/Talks/Chatterjee.pdf](https://math.utah.edu/~firas/FPD18/Talks/Chatterjee.pdf)【snippet】
- `Euclidean Yang-Mills theories are scaling limits of lattice gauge theories (probability-theoretic open problem).` —— Chatterjee, Columbia CFT18 讲义，[Chatterjee.pdf](https://www.math.columbia.edu/department/probability/CFT18/Chatterjee.pdf)【snippet】
- 书内逐字：**未取数**（Seiler LNP 159 / Chatterjee 书均未取到正文；检索式见 §5、§6）。**不以二手代替。**

---

## 未取数清单（供下一轮补齐，勿以记忆填充）

1. Dimock–Hurd `observable algebra` 一文：`au:Dimock au:Hurd "observable algebra"`。
2. Balaban III（CMP 118）、IV（CMP 122）：`Balaban "lattice gauge field theories" III / IV`。
3. Tomita：`au:Tomita lattice gauge infrared`（限 CMP/LNP 1980–1990）。
4. King：`au:King C lattice gauge strong coupling`。
5. Seiler 书 ISBN + 4D 质量隙段落逐字：Google Books 全文 / Springer 书页。
6. Chatterjee 书正文质量隙段落逐字：书 PDF。
7. PDF 正文抓取器在本会话返回 `unsupported content type "application/pdf"`；所有【snippet】引文需下一轮以可读正文复核后方可升格为正式引文。

## 对本项目的意义（一句话）

把"严格构造 4D 规范理论"当入口是**方向性错误**：这条线交付的是 UV 稳定与红外冻结，`m` 从未出现在其结论里；我们要的 `a`-一致正隙必须自建（`G″(0)>0` / 粗模谱隙下界那一支），**不能**指望从 Balaban–MRS 线上"漏掉的式子"里取。

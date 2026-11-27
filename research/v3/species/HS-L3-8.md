# HS-L3-8 — 文献扫：4D 纯 YM 的 `a`-一致矩界／紧性，有没有已发表结果？

> **轮**：HS-8 · 子代理 `agent-8` · 工作目录 `D:\Dev\math` · 只新建本文件（不改既有文件、不 commit、零 MC、不扩范围）。
> **任务**：只回答一件事 —— Balaban / Magnen–Rivasseau–Sénéor (`CMP` **155**, 1993) / Chatterjee / SZZ / Cao–Chatterjee 中，有没有**已发表**的 4D 纯 YM `a`-一致矩界或 `a→0` 紧性。
> **读级**：`[摘要]` = 读到官方/Crossref/INSPIRE 摘要逐字；`[元数据]` = 只核到卷-页-DOI-标题；`[评]` = 读到第三方评论（zbMATH）片段。**凡未标 `[正文]` 者正文均未逐字核**（Springer 付费墙；本轮 0 篇读全文）。
> **标签**：只交 `[T]`／`[NG]`。`[T]` = 书目事实或摘要逐字已核；`[NG]` = 具名缺步／对具名类 0 供给。

---

## §0 判词（先行；每行带「下一手」）

| # | 断言 | 判 | 依据（读级） | 下一手 |
|---|---|---|---|---|
| 0.1 | **五条具名候选，无一条**给出 4D 纯 YM 的 gauge-invariant `a`-一致矩界（`L3` 的 `(U1)`–`(U3)`），也**无一条**给出 `a→0` 的测度紧性。五条的"一致"全挂在别的东西上：`Z`/有效作用量、gauge-fixed Schwinger 函数、固定 `a` 的 `β`-窗口、3D 状态空间 | **`[NG]`** | §1 逐条 | 把 §3 的"谁供给哪一格"表并进 `HQ1-L3 §2.2` 的供给表，替换"Balaban 到不了大场"那一行 |
| 0.2 | **仓内错记**：`G1-L3.md:16`／`G1-ASSEMBLE-100.md:104` 写 Balaban "explicitly … short of the large-field problem"。**Balaban 已完成 4D 纯规范场 UV 稳定** —— `CMP` **122** (1989) 355 摘要逐字："This allows us to **complete the proof of the ultraviolet stability of four-dimensional pure gauge field theories**, as formulated in Theorem 1." | **`[T]`** `[摘要]` | DOI `10.1007/BF01238433`；同摘要见 NII-REO 卷期页、zbMATH Zbl 0704.58060 | 全仓"缺大场"表述改为"大场 `R` 操作已做（`CMP` 122.175 + 122.355）；缺的是 **`Z`→矩** 与 gauge 不变性" |
| 0.3 | **但 Balaban 的一致是 `Z` 级**：zbMATH 评 `CMP` 122,355 逐字："it follows the **ultraviolet stability bounds of the same type as for superrenormalizable models**"。即：配分函数/有效作用量密度的 `e^{O(1)}` 界，**不是** gauge-invariant smeared 矩的下界，**不是**测度 | **`[T]`** `[摘要]`+`[评]` | 同上 + zbMATH 0704.58060 评语 | `Z`→矩 那一步是**新引理**，不是文献搬运；先判它是否就是 `L3` 的 `(U1)` 本体（见 §3.1） |
| 0.4 | **唯一"`a`-一致矩界"的已发表实例是 2D 且 gauge-fixed**：Chandra–Chevyrev 摘要逐字 "we show that **moments … are bounded uniformly in the lattice spacing**"，对象是 2D 格点 Abelian Higgs 的 **gauge field marginal in a fixed gauge** | **`[T]`** `[摘要]` | arXiv:2207.05443 → `CMP` **405** (2024) art. 143，DOI `10.1007/s00220-024-05008-4` | 这是 `(U1)` 的**方法模板**（Besov 范数矩 + 抗磁不等式）；先查它的常数是否 `L`-一致、是否只对 gauge-fixed 边际成立 |
| 0.5 | **MRS 1993 不是矩界**：摘要逐字 "We **provide the basis for** a rigorous construction … with a **fixed infrared cutoff** but no ultraviolet cutoff, in a **regularized axial gauge**"，且"Trivial topological sector"／纯 `SU(2)` | **`[T]`** `[摘要]` | DOI `10.1007/BF02097397`，`CMP` **155**(2) (1993) 325–383（Crossref 核：页 325–383，非 417–430） | 继续按 `HQ1-L3 §2.2` 记"无 gauge-invariant OS"；**不要**把它的 UV 控制记成 `L3` 供给 |
| 0.6 | **SZZ / Cao–Chatterjee / Chatterjee 三条都不在轴上**：SZZ 固定 `a`、强耦合 `|β|<1/16(d−1)`（`SU(N)`，任意 `d>1`）；Cao–Chatterjee 是 **3D** 且条件于 GFF 假设；Chatterjee 4D 唯一 `a→0` 结果要 `g=O(ε^{50d})→0`，极限是**有质量高斯场**，自述非高斯开 | **`[T]`** `[摘要]` | arXiv:2204.12737 / `CMP` **400** (2023) 805–851，DOI `10.1007/s00220-022-04609-1`；arXiv:2111.12813；arXiv:2401.10507 / `Prob. Math. Phys.` **7** (2026) 339–381，DOI `10.2140/pmp.2026.7.339` | SZZ 只当"固定 `a` 反例"用；Chatterjee 4D 记为"`a→0` 有先例、但 `g→0` + Higgs ⇒ 不可搬" |
| 0.7 | ⇒ **`L3` 的 `(U1)`（gauge-invariant `a`-一致矩界）在已发表文献中 0 供给**；`a→0` 的 4D 紧性同样 0 供给。与 `G1-L3.md:21`（"no published result closes L3"）**一致**，本轮以具名 ID 逐条补上证 | **`[NG]`** | §1–§2 + §4 检索式（命中 0/1，失败标记 0） | 见 §3.2：下一步不是再扫，而是判 `Z` 级界能否升到 `(U1)`（一个可写的引理，或一条 `[NG]`） |

---

## §1 五个具名候选，逐条核过（全部只读摘要/元数据）

| # | 文献 | ID / DOI | 读级 | 它的"一致"是什么 | 给 `L3` 哪一格 |
|---|---|---|---|---|---|
| B1 | Balaban, *Renormalization group approach to lattice gauge field theories. 1. Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions*, `CMP` **109** (1987) 249–301 | DOI `10.1007/BF01215223` | `[摘要]` | 摘要逐字："**restricted to small field approximation**"；一阶 RG 变换的团簇展开 + `β`-函数 + 耦合常数重整 | 无（小场前半） |
| B2 | Balaban, *Renormalization group approach to lattice gauge field theories*（II）, `CMP` **116** (1988) 1–22 | DOI `10.1007/BF01239022` | `[元数据]`（Crossref 卷/页 OK；摘要未取到） | 系列续篇 | 无 |
| B3 | Balaban, *Large field renormalization. I. The basic step of the `R` operation*, `CMP` **122** (1989) 175–202 | DOI `10.1007/BF01257412` | `[摘要]` | 摘要逐字："removes the **main obstacle** to prove the ultraviolet stability of four-dimensional gauge field theories" | 上游（大场侧） |
| B4 | Balaban, *Large field renormalization. II*, `CMP` **122** (1989) 355–392 | DOI `10.1007/BF01238433` | `[摘要]`+`[评]` | 摘要逐字："**complete the proof of the ultraviolet stability of four-dimensional pure gauge field theories**"；评语：`Z` 级界（超可重整型） | **`(U2)` 的上游，非 `(U1)`**（见 §3.1） |
| B5 | Magnen–Rivasseau–Sénéor, *Construction of `YM₄` with an infrared cutoff*, `CMP` **155** (1993) 325–383 | DOI `10.1007/BF02097397` | `[摘要]`（INSPIRE 存 Springer 摘要全文） | "basis for" 构造：固定 IR、无 UV cutoff、正则化轴规范、平凡拓扑、`SU(2)`；大场正性 + 小场背景场 + Slavnov 检验 | 无（gauge-fixed、无 OS、无矩族） |
| C1 | Chatterjee, *Rigorous solution of strongly coupled `SO(N)` lattice gauge theory in the large `N` limit*, `CMP` **366** (2019) 203–268 | DOI `10.1007/s00220-019-03353-3`；arXiv:1502.07719 | `[元数据]`（Crossref 全记录）＋摘要经 SZZ 引文页 | 固定小 `β`、`N→∞` 的 Wilson 环弦和；**无 `a`** | 无 |
| S1 | Shen–Zhu–Zhu, arXiv:2204.12737 → `CMP` **400** (2023) 805–851 | DOI `10.1007/s00220-022-04609-1` | `[摘要]` | 固定 `a`、`|β|<1/16(d−1)`（`SU(N)`）；LSI/PI + 指数遍历 + 无穷体积测度 + 正质量隙 | 无（强耦合、`a` 固定） |
| CC1 | Cao–Chatterjee, *A state space for 3D Euclidean Yang–Mills theories*, `CMP` **405** (2024) art. 3 | arXiv:2111.12813；DOI `10.1007/s00220-023-04870-y`（承 `LIT-CHECK.md:103`） | `[摘要]` | **3D** 状态空间；"**provided that they indeed behave like the GFF**" ⇒ 条件 | 无 |
| C2 | Chatterjee, *A scaling limit of `SU(2)` lattice Yang–Mills–Higgs theory*, `Prob. Math. Phys.` **7** (2026) 339–381 | arXiv:2401.10507；DOI `10.2140/pmp.2026.7.339` | `[摘要]` | **唯一 4D `a→0` 标度极限**：`g=O(ε^{50d})`、`αg=cε`、unitary gauge ⇒ **有质量高斯场**；"non-Gaussian … remains open" | 4D `a→0` 有先例；`g→0` ⇒ 不可搬到纯 YM |

**一行判**：B1–B4 在 `a` 轴上（对），但给的是**有效作用量/`Z` 的一致**；B5 在 UV 轴上但 gauge-fixed 且"basis"；C1/S1/CC1/C2 全在**别的轴**（`N`、固定 `a`、3D、`g→0`）。

---

## §2 最近邻三条（不属于题目点名的五个，但必须记账）

| # | 文献 | ID / DOI | 读级 | 为什么必须记 |
|---|---|---|---|---|
| N1 | Chandra–Chevyrev, *Gauge field marginal of an Abelian Higgs model*, `CMP` **405** (2024) art. 143 | arXiv:2207.05443；DOI `10.1007/s00220-024-05008-4` | `[摘要]` | **已发表文献里唯一**"moments … bounded uniformly in the lattice spacing"的实例 ⇒ `(U1)` 的**方法模板**；但 2D、Abelian、**固定规范**、对象是 gauge field marginal 而非纯 YM gauge-invariant 矩 |
| N2 | Brydges–Fröhlich–Seiler, *On the construction of quantized gauge fields* I–III, `Ann. Phys.` **121** (1979) 227；`CMP` **71** (1980) 159；`CMP` **79** (1981) 353 | DOI `10.1016/0003-4916(79)90098-8`；`10.1007/BF01197918`；`10.1007/BF01208500` | `[元数据]`（经 Chatterjee `CMP` 366 参考文献表核到卷/页/DOI） | 2D／3D 构造 + 格点逼近收敛；**不是 4D**；`HQ1-L3 §2.2` 未列，但它是"`a`-一致矩 + 测度"的唯二先例族之一 |
| N3 | `Zenodo` 17247200（C. Maynard, 2025-10，单作者，题："Uniform Bounds for Flowed Gauge-Invariant Composites in 4D Yang-Mills (Lemma A)"） | DOI `10.5281/zenodo.17247200` | `[摘要]`（Zenodo 记录页） | **最贴近 `(U1)` 的公开声明**：4D `SU(N)`、flowed gauge-invariant composite 的矩一致有界 ⇒ tightness。但：**无 arXiv ID、无同行评议、单作者"三引理纲领"**，且只在 `t(a)=a²|log a|^p` 窗口成立 ⇒ 按 `L4-PREPRINTS §4` 纪律**不列为可核条目**，仅作线索 |

---

## §3 谁供给 `L3` 的哪一格（这是本轮的实际产出）

### 3.1 `Z` 级 ≠ 矩级（本文件最承重的一行）

Balaban 的"4D UV 稳定"给的是 `Z`（与有效作用量密度）在 RG 步数上的一致界，型如超可重整模型。`L3` 要的是 `|S^{(n)}_{L,a}(f_1..f_n)| ≤ C_n(L)∏‖f_i‖_{C^k}`，`∀a≤a₀`，**gauge-invariant**、**smeared**、带 `n!` 增长。二者的差距是三条，**每条都必须单独补**：

1. **从 `Z` 到矩**：`Z` 的一致界不含任何 smeared 观测量的高阶矩；`e^{O(1)}` 不会自动变成 `C_n(L)`。
2. **gauge 不变性**：Balaban 的机制建立在轴向/正则化规范 + 背景场分解上；把 gauge-fixed 的控制搬到 gauge-invariant 的 `S^{(n)}` 需要一次显式的"不变化"（或直接在不变量上做 RG）。
3. **`(U2)` 的方案常数**：`a^{-4}` 减除必须与 `L2` 的方案 C 匹配；文献里没有任何一条在**同一方案**下写出这个常数的一致递推（B4 的 counterterm 递推是 Balaban 方案，不是方案 C）。

⇒ 三条都不是"查漏"，是**新引理**。故 0.1 的 `[NG]` 不是"没搜到"，是"该格在已发表文献中确为空，且缺口已定位到三条具名步骤"。

### 3.2 酸测：这些文献能不能被"借"成 `L3`

| 借法 | 结论 | 一行理由 |
|---|---|---|
| 用 B4 的 UV 稳定 ⇒ `(U1)` | **不可**（差 §3.1 三条） | `Z` 级界 vs 矩级界；且 gauge-fixed |
| 用 N1 的 2D 模板搬到 4D | **不可直接搬**，可搬方法 | 2D 的 `a`-一致靠 2D 传播子的对数结构；4D 是 marginal，`β_j→∞` 塔在 2D 无对应物 |
| 用 C2 的 4D 极限当先例 | **只买"存在 4D `a→0` 极限"** | `g=O(ε^{50d})→0`，非阿贝尔相互作用被一起送走；极限是有质量**高斯** |
| 用 N3 当供给 | **不列为供给** | 无 arXiv/无评审；且 flow 窗口 `t(a)` 把"sharp moment bound"换成了"flowed operator"的界（仓内 `G1-L3-VIA-FLOW` 已判 flow 不改 sharp bound） |

---

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际检索式（逐条） | 命中 |
|---|---|---|---|---|
| 1 | `web_fetch` → Crossref | **OK**（HTTP 200，JSON） | `api.crossref.org/works/10.1007/BF02097397`；`…/10.1007/s00220-019-03353-3`；`…/10.1007/s00220-022-04609-1`；`works?query.author=Balaban&query.bibliographic=lattice+gauge+renormalization+group&rows=20` | 1 / 1 / 1 / total 73（返 20） |
| 2 | `web_fetch` → arXiv abs | **OK**（HTTP 200） | `arxiv.org/abs/2204.12737`；`/abs/2111.12813`；`/abs/2207.05443`；`/abs/2401.10507`；`/abs/2611.11534` | 5 页摘要逐字 |
| 3 | `web_fetch` → arXiv API | **OK**（XML） | ①`abs:"uniform in the lattice spacing"` ②`abs:"moment bounds" AND abs:"Yang-Mills"` ③`abs:"tightness" AND abs:"Yang-Mills"` ④`au:Dimock_J AND abs:Balaban` | ① **1**（唯一命中 arXiv:2510.08006，**假阳性**：拓扑磁化率数值文）② **0** ③ **7**（无一是 4D 纯 YM 的 `a`-一致紧性；含 Gauvin 2503.15539 只给**判据**且自述未证）④ **8**（Dimock 的 Balaban 复述系列是 3D `φ⁴`，非规范场） |
| 4 | `web_fetch` → INSPIRE API | **OK**（JSON） | ①`a Magnen and a Rivasseau and a Seneor and t infrared cutoff` ②`a Balaban, T. and t renormalization group and t lattice` ③`a Balaban, T. and t large field` ④`t "renormalization group approach to lattice gauge"` ⑤`doi:10.1007/BF01239022` | ①**1**（摘要全文到手）②**3** ③**3**（**B3/B4 摘要逐字**，本轮关键）④**1** ⑤**0**（INSPIRE 不收录 CMP 116，故 B2 只能 `[元数据]`） |
| 5 | `web_search` | **OK**（10 源/查询） | ①`"uniform in the lattice spacing" moment bounds Yang-Mills gauge invariant 4D` ②`Balaban ultraviolet stability … "uniform in the lattice spacing"` ③`a-uniform bounds Schwinger functions lattice Yang-Mills continuum limit compactness` ④`Dimock "The renormalization group according to Balaban" …` ⑤`Balaban four dimensional ultraviolet stability not proved only three dimensional` ⑥`"uniform in the lattice spacing" 4D pure Yang-Mills tightness … 2025` ⑦`Balaban 1989 "ultraviolet stability" … complete … assessment` | 用到的定位：N1 摘要句；B4 摘要（NII-REO／zbMATH Zbl 0704.58060 评语）；Zenodo 17247200；Gauvin 2503.15539 自述 |
| 6 | `grep`（仓内） | **OK**（无 provider failure） | `Magnen\|Rivasseau\|Sénéor` @ `research/**/*.md`；`Chatterjee` @ 同 | **64** / **154** |
| 7 | `read` | **OK** | `research/v3/species/HQ1-L3.md`(全 286)、《`A5-LIT-UNIFORM.md`》(全 102)、《`L4-PREPRINTS.md`》(全 187)、《`research/unconditional/G1-L3.md`》(全 29) | 全部返回正文 |
| 8 | **失败标记** | ①`web_fetch link.springer.com/article/10.1007/BF01215223` → `cross-origin redirect to idp.springer.com`，**SOFT-FAIL**（改走 INSPIRE 摘要，已补上）②`web_fetch http://export.arxiv.org/api/...` → cross-origin redirect **SOFT-FAIL**，改 `https://export.arxiv.org` 后 **OK** ③INSPIRE `doi:10.1007/BF01239022` → **0**（覆盖缺口，非工具故障） | — | — | `TOOL-FAIL` 硬失败 **0**；空响应 **0** |
| 9 | 未调用 | `pwsh` / `py` / Lean / MC / `git` / `skill`（本任务为文献裁决，catalog 的 `academic-*`/`deep-research` 不适用） | — | — |

**检索式的边界（诚实）**：以上 6 条 arXiv/Crossref/INSPIRE 式 + 7 条 web 式覆盖的是**摘要层**。所以本文件的"0 供给"读作：**在这 13 条式下、摘要层无一条 4D 纯 YM 的 gauge-invariant `a`-一致矩界或 `a→0` 紧性**。它不是"不存在"的证明，而是"以具名 ID 逐条点过、每条都说清它的一致挂在哪个轴上"。

---

## §5 诚实边界（6 条）

1. **0 篇读正文**：全部 `[摘要]`／`[元数据]`／`[评]`。B4 的"`Z` 级"是 zbMATH 评语 + 摘要措辞的合成判读，**未开正文**核 Theorem 1 的完整假设（尤其是 counterterm 是否沿 `β(a)` 递推、体积是否固定）。
2. **B2（`CMP` 116）摘要未取到**：故 B1–B4 的"系列分工"有一步是标题＋卷页推断，记 `[元数据]`，不作 `[T]`。
3. **0.2 是对仓内的更正，不是对 Balaban 的评价**：`HQ1-L3 §2.2` 的 M1 行、`G1-L3.md:16`、`G1-ASSEMBLE-100.md:104` 三处"缺大场"表述与 B3/B4 摘要直接冲突；本文件只指出冲突，不改那三份文件（纪律）。
4. **"`Z` 级 ⇒ 矩级"这一步本文件不证**：§3.1 把它拆成三条具名步骤，是**缺口定位**，不是三条的证明或否证。
5. **N3（Zenodo）不作引用**：按 `L4-PREPRINTS §4` 的既有纪律（无 arXiv ID、无评审）；记它只为让下一个代理不重复发现它。
6. **无 QED**：本文件不产一条估计、不升格任何 `[T]` 到更大范围；`0.1`/`0.7` 的 `[NG]` 是"具名格在已发表文献中为空"，不是"该格不可证"（仓内 `G1-UNIFORM-ZA.md:9` 的无内禀 no-go 判词不受影响）。

*HS-L3-8 · 只新建本文件；既有文件保持原样；未 `git commit`；零 MC；未动 `docs/position-current.html`；`≤130` 行。*

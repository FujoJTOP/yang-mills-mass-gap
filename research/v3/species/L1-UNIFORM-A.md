# L1-UNIFORM-A — `a`-一致（格距一致）的 4D 质量隙：文献裁决

> 角度只有一个：**有没有已发表/预印本结果形如「`a→0` 时 `m ≥ c·Λ`（`c` 纯数、与 `a` 无关）」或等价 `m a ≥ c'·aΛ`（`p=1`）**。
> 工具：`web_search`（多轮换词）、`web_fetch`（arXiv API / arXiv abs 页 / Crossref REST）。
> 证据等级：**[摘要]** = 逐字读到官方摘要；**[元数据]** = 只核到 Crossref 的题名/卷/页/DOI；**[二手]** = 只见于他人转述，未核一手。
> 检索窗口：2026-11-25（本会话）。**结论先给。**

---

## 0. 裁决

**主判：`(c)` 未查到 —— 同时含 `(b)` 维度不够。**

- 没有任何已发表/预印本定理给出 4D 纯 `SU(N)` 的 `m ≥ c·Λ`（`c` 与 `a` 无关）。**未查到。**
- 唯一两条**声称** `a`-一致（跨尺度搬运谱隙）的条目，一条是**已撤稿**（`2506.00284`），一条**摘要里没有任何显式下界/式号**（`2606.19362`）。二者都不能当定理用。
- 能逐字引到的「一致」全部挂在**别的轴上**：体积 `L`（`Osterwalder–Seiler`）、**有效作用量**而非谱隙（`Balaban`）、**IR cutoff 固定下的 UV 移除**（`Magnen–Rivasseau–Sénéor`）、**数值连续外推**（`Bonanno et al.`、`Loan et al.`）。
- 2D/3D 条目存在（数值 + 程序），但**维度不够**，且 2D 的 `g` 自带质量量纲、无 dimensional transmutation，`a m → 0`，与 4D 的 `p=1` 问题不同型（此句为我方推理，非文献）。

**与内部前序一致**：`research/v3/species/A5-LIT-UNIFORM.md` 的主判同为 `(b)/(c)`。本轮独立检索**未推翻**它，只补了 2025–2026 两条「声称型」新条目。

---

## 1. 候选逐条（含逐字引用 + 它统一的轴 + 有无 `m` 下界）

### C1 — Faizal & Shabir 2026：**声称**跨尺度搬运「uniform gap」，无显式界
- M. Faizal, A. Shabir, "Reflection-Positive Construction of a Four-Dimensional SU(N) Yang-Mills Theory with Mass Gap and Confinement", **arXiv:2606.19362**；`Fortschr. Phys.` **74** (2026) 4, e70097；DOI `10.1002/prop.70097`。**[摘要]**
- 摘录（逐字，一句连续）：
  > "starting from a reflection-positive lattice formulation of pure SU(N) Yang-Mills theory we obtain a transfer operator with a uniform gap, while large Wilson loops already show an area law by means of convergent character (polymer) expansions; a finite-range, gauge-covariant multiscale analysis then carries these features from one scale to the next with interlaced inequalities whose small defects can be summed, so that exponential clustering and a strictly positive string tension endure in the continuum"
- **统一的轴**：声称是**尺度/RG 步**（即 `a` 方向），且点名 4D 纯 `SU(N)`。
- **有 `m` 的下界吗**：**摘要里没有**。无 `m ≥ c·Λ`、无 `c`、无式号；只有 "uniform gap" 与 "strictly positive string tension" 的定性断言。正文（4 篇 `IJGMMP` 附录，2026 年 2650111–2650114）**未精读**，故**不能**声称它给出了可核的 `p=1` 不等式。
- **用法**：只能作为「有人声称过」的记录，**不可作定理引用**。

### C2 — Jacobsen 2025：逐字命中我们要的形式，但**已被 arXiv 撤稿**
- D. C. Jacobsen, "A Constructive Proof of Existence and Mass Gap for Pure SU(3) Yang-Mills in Four-Dimensional Space-Time", **arXiv:2506.00284**（v2）。**[摘要]**
- 摘录（逐字）：
  > "A convergent joint polymer expansion on the Wilson lattice provides uniform control of the continuum (a->0) and infinite-volume (L->inf) limits"
- **统一的轴**：`a` **与** `L` 双轴（正是我们想要的形状）。
- **有下界吗**：摘要称 "isolates a strictly positive glueball mass m0"，但**无式号、无常数**。
- **判死**：arXiv 注释逐字为
  > "This submission has been withdrawn because it does not meet arXiv's research content quality standards"
  故**排除**（撤稿、无有效内容）。

### C3 — Mihajlovic 2026：**明确宣布** `a`/体积一致谱隙未建立（最新 no-go 式自述）
- A. Mihajlovic, "SU(2) Yang-Mills: Spectral Bounds, Vacuum Concentration, and Criteria for the Continuum Limit", **arXiv:2611.11534** [math-ph], 414 pp, 2026-11-08。**[摘要]**
- 摘录（逐字）：
  > "A nontrivial four-dimensional spatial-cutoff construction and a positive spectral bound along it are not established; no exact continuum mass is asserted."
- 另（逐字，同一摘要）：
  > "collective projection and Schur estimates isolate the additional input needed for a volume-uniform gap."
- **统一的轴**：空间 cutoff `a` + 体积 `L`（双轴），对象是 transfer/物理时间动力学。
- **有下界吗**：**没有**，且自述「未建立」。**[二手] 线索**：检索索引里出现该文正文句 "Its lattice gap is m_j = a_j/j = 1/j^2, so m_j/a_j = 1/j → 0"——即它举了「格点隙随 `a` 归一后趋零」的例子；此句**我只从检索片段见到，未在正文核对**，故不作为引用，仅记为线索。
- **性质**：单作者预印本，非同行评议；价值在于它**自己承认** `a`-一致谱隙是缺口。

### C4 — Balaban 1984–88：轴是 `a`，但对象是**有效作用量**，不是谱隙
- T. Balaban, "Renormalization group approach to lattice gauge field theories", `Commun. Math. Phys.` **116** (1988) 1–22, DOI `10.1007/BF01239022`；同系列 `CMP` **95** (1984) 17–40 `10.1007/BF01215753`、**99** (1985) 389–434 `10.1007/BF01240355`、**102** (1985) 255–275 `10.1007/BF01229380`、**102** (1985) 277–309 `10.1007/BF01229381`。**[元数据，Crossref 核过]**
- **统一的轴**：**RG 步数 = 逐层缩小的 `a`**（这就是文献里 "uniform in the lattice spacing" 的实际落点）。
- **有 `m` 的下界吗**：**没有**。控制的是配分函数/有效作用的解析性与 UV 稳定性；谱隙从未被推出，也未接 OS 重构。

### C5 — Osterwalder–Seiler 1978：轴是**体积 `L`**，`a` 固定
- K. Osterwalder, E. Seiler, "Gauge field theories on a lattice", `Annals of Physics` **110** (1978) 440–471, DOI `10.1016/0003-4916(78)90039-8`。**[元数据，Crossref 核过：卷 110、页 440–471、1978-02]**
- **摘要逐字文本本轮未取到**（Elsevier 页未开；Crossref 该条无 `abstract` 字段）——按纪律**不补写、不凭记忆引**。
- **统一的轴**：体积 `L`（强耦合区），**不是 `a`**。
- **有下界吗**：在其成立范围内（强耦合、固定 `a`）有正隙；**不含任何 `a` 方向命题**。

### C6 — Magnen–Rivasseau–Sénéor 1993：轴是 UV cutoff 移除，但 IR cutoff 永不解除
- J. Magnen, V. Rivasseau, R. Sénéor, "Construction of YM₄ with an infrared cutoff", `Commun. Math. Phys.` **155** (1993) 325–383, DOI `10.1007/BF02097397`。**[元数据]**
- **统一的轴**：`a→0`（去 UV cutoff），方向对。
- **有下界吗**：**没有**。IR cutoff 固定、规范固定，不给规范不变 `0^{++}` 谱隙。

### C7 — Dimock（Balaban 型 RG，QED₃）：轴是 `a`，自述「只做第一步 / 不控余项」
- J. Dimock, arXiv:`math-ph/0210020`；arXiv:`math-ph/0407063`；arXiv:`1712.10029`（均 3D QED）。**[摘要]**
- 摘录（逐字，`math-ph/0210020`）：
  > "We start with the lattice gauge theory on a toroidal lattice and seek to control the singularities as the lattice spacing is taken to zero … In this paper we explore the first step in detail."
- **统一的轴**：`a`。
- **有下界吗**：**没有**；且 3D + QED（非 4D 纯 `SU(N)`）。

---

## 2. `(b)` 维度不够：2D/3D 的条目

| 条目 | ID / DOI | 逐字（摘要） | 统一的轴 | 有 `m` 下界？ |
|---|---|---|---|---|
| Loan–Brunner–Hamer 2002, 3D compact U(1) 胶球 | arXiv:`hep-lat/0209161` | "Convincing evidence of the scaling behaviour in the antisymmetric mass gap is observed over the range $1.4<β<2.25$." | `β`（↔`a`），**数值** | 否（拟合 scaling，非下界） |
| Balaban 1985, 3D 纯规范 UV 稳定性 | `CMP` 102, 255–275 | 见 C4 系列 | `a`（有效作用量） | 否 |
| Dimock 3D QED RG | 见 C7 | 见 C7 | `a` | 否 |

- **3D 严格结果**（Polyakov 1977 / Göpfert–Mack 1982 一类）**本轮未取到逐字摘要**，按纪律**不列为候选**。
- **2D**：本轮**未取到**可逐字引用的「`a`-一致」条目。结构上 2D 的 `g` 自带质量量纲、无 dimensional transmutation，`a m → 0`，与 4D 的 `p=1` 问题不同型（我方推理）。

---

## 3. 数值（非定理，但说明「实践上 `a`-一致」长什么样）

- C. Bonanno, C. Bonati, M. Papace, D. Vadacchino, "The imaginary-θ dependence of the SU(N) spectrum", **arXiv:2411.14022**（Lattice 2024）。**[摘要]**
  > "The coefficient of the $\mathcal{O}(θ^2)$ term in the Taylor expansion of the spectrum around θ=0 could be obtained in the continuum limit for N=3, and on two fairly fine lattices for N=6."
  → 轴 = `a`（连续外推）；**无下界**，纯数值。
- 同类（glueball 谱的连续外推，给出 `m/√σ` 平台）属数值证据，**不能替代 `m ≥ c·Λ`**。

---

## 4. 实际用过的检索式（可复现，≥8）

**web_search**（逐字）：
1. `"mass gap" "uniform in the lattice spacing"`
2. `"uniform mass gap" lattice gauge continuum limit`
3. `"continuum limit" "mass gap" Yang-Mills rigorous exponential decay uniform lattice spacing`
4. `"mass gap" "uniformly in a" lattice`
5. `site:arxiv.org mass gap uniform lattice spacing continuum`
6. `Osterwalder Seiler 1978 gauge field theories on a lattice mass gap exponential decay uniform in volume`
7. `Balaban renormalization group lattice gauge mass gap uniform cutoff`
8. `Brydges Frohlich Seiler convergence of the lattice approximation mass gap`
9. `"lower bound on the mass gap" lattice Yang-Mills arXiv`
10. `"mass gap survives the continuum limit" Yang-Mills`
11. `luminal mass gap lattice Yang-Mills bound`
12. `Bauerschmidt Dagallier upper bound mass gap lattice gauge theory arXiv`
13. `"uniform in the lattice spacing" "mass gap" exact phrase theorem`

**arXiv API**（`export.arxiv.org/api/query`）：
14. `abs:"mass gap" AND abs:"lattice spacing"` → **19** 命中，逐条看：无一条给 4D `a`-一致下界。
15. `abs:"mass gap" AND abs:"continuum limit" AND abs:"lattice gauge"` → **4** 命中（hep-lat/9704017、hep-th/9509032、hep-lat/0307029、hep-th/0006022），均为数值/方法，无定理。
16. `abs:"mass gap" AND abs:"uniform" AND abs:"lattice"` → **12** 命中，其中 2 条为 C1/C2，其余无关（石墨烯、非厄米、活性物质等）。

**Crossref**：17. `query.bibliographic=Osterwalder Seiler gauge field theories on a lattice` → 核到 C5 元数据。

**未完成的检索（诚实记录）**：`ti:"mass gap" AND abs:lattice`（按日期倒排）、`abs:"exponential decay" AND abs:"lattice gauge" AND abs:"uniform"`、`au:Dimock`/`au:Balaban` 全量列举 —— 这三条在发起时被 arXiv 返回 **HTTP 429（限流）**，本 10 分钟窗口内未重跑；C7/C4 的条目依赖内部前序 `A5-LIT-UNIFORM.md` 已核过的 Crossref 记录 + 本轮抽查。

---

## 5. 被排除的候选（各一句理由）

1. **arXiv:2506.00284**（Jacobsen）— **已撤稿**（arXiv admin note：不符合内容质量标准），且摘要无式号。
2. **arXiv:2606.19362**（Faizal–Shabir）— 声称 `a`-一致「uniform gap」，但摘要**无任何显式常数/不等式**，支撑文献为其同团队 4 篇 `IJGMMP`，独立复现未见；**降级为「声称」**（仍列 C1 备查，不作定理）。
3. **Zenodo / OSF / Academia 自印稿**（如 "Uniform Spectral Lower Bounds in Discrete Gauge Models"、"Volume-Uniform Poincaré Inequalities as a Mechanism for the Yang–Mills Mass Gap"、"A conditional continuum theorem…"、"Explicit strong-coupling mass gap and a robust landing interface"）— **非同行评议自印**，无 DOI 期刊归属，部分为 AI 生成文风，**不引用**。
4. **aiwiki / valency hub / computoergosum 等聚合站**— 二手且无出处链，**不引用**。
5. **MDPI `IJT` 卷 2 文章 00002**（`ijt-02-00002`，PDF 版本戳 2025；声称 "Δ ≥ C_G g²/log(1/g²) Λ_QCD"）— 期刊与该式均无严格推导可核，且这是**上界型**（`g²/log` 随 `g→0` 趋零），**结构上不是 `a`-一致正下界**，排除。
6. **`hep-lat/0509043`（Wolff et al.）**— O(N) σ 模型（非规范理论）+ 连续外推方法学，无下界。
7. **`hep-lat/9511004`（Hasenbusch–Horgan, SO(4) PCM）**— 2D 主手征模型，维度不够且明确报告外推可能错误。
8. **`cond-mat/9901184`（XXZ 链 scaling limit）**— 自旋链，非规范理论。
9. **`2610.34277`（Tomiya, Gross–Neveu Fisher 度量）**— 2D Gross–Neveu，维度不够且无质量隙下界。
10. **Bauerschmidt–Dagallier（LSI for U(1)₄）**— 仅见 [二手] 提及「proved LSI for U(1) in d=4」；**未核到一手 arXiv ID**，且 LSI 是泛函不等式、不是 `m ≥ cΛ`，**列而不引**。

---

## 6. 给 `L1` 的一句话

> 已被引用的「`a`-一致」只有两种落点：**有效作用量**（Balaban）与**数值连续外推**（Bonanno、Loan）；**谱隙**这一侧没有任何一条已发表定理，两条声称型预印本一条撤稿、一条无显式界。故 `m ≥ c·Λ`（`p=1`）在文献中**未查到**；`L1` 的缺口是领域边界，不是检索遗漏。

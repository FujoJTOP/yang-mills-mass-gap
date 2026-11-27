# HT6-LIT — 4D 纯 YM 的 `a`-一致矩界 / 紫外紧性：判词 `[NG]` + 一条待核声明

> **轮**：HT6 · 子代理 `ht6-lit` · 工作目录 `D:\Dev\math`。只新建本文件；不改既有文件、不 commit、零 MC。
> **靶（唯一）**：`L3^{UV}`（一行式见 `HQ1-L3.md:37–41`：`∃a₀∃{C_n}∀a≤a₀∀f: |S^{(n)}_{L,a}(f)|≤C_n(L)∏‖f_i‖` + `n!` 增长）在已发表文献里有没有人做过。
> **读级**：`【读摘要】` = arXiv/INSPIRE/Crossref 摘要逐字；**本轮 0 篇读正文**（Springer/Wiley 付费墙）。凡断言"某文没给矩界"者，均由该文摘要的量词形状判，非正文逐条核 —— 见 §2.3。
> **标签**：只交 `[T]`／`[NG]`。`[T]` = 书目事实或摘要逐字已核；`[NG]` = 具名格在已发表文献中 0 供给 / 具名缺步。
> **时间单位是轮，不是天。**

---

## §0 判词（先行；每行带「下一手」）

| # | 断言 | 判 | 依据（读级） | 下一手 |
|---|---|---|---|---|
| **0.1** | **4D 纯 YM 的 gauge-invariant `a`-一致矩界、以及 `a→0` 的测度紧性，在已发表文献中 0 条供给** ⇒ `L3^{UV}` 是**文献空位**，且这个空位**已被点过名**（Balaban 的 UV 稳定只到 `Z` 级，MRS 只到固定规范） | **`[NG]`** | §1（11 条具名）+ §3（18 条检索式，命中 0 条 4D 纯 YM 矩界） | 不再扫。下一手是判 `Z`→矩这一步能否写成引理（`HT6` 之后的活儿），见 0.2 |
| **0.2** | **Balaban 系列确实做到了 `a`-一致**，但一致的**对象是 `Z`/有效作用量**：`CMP` **122** (1989) 355 摘要逐字 "complete the proof of the ultraviolet stability of four-dimensional pure gauge field theories, as formulated in Theorem 1"；`CMP` **122** 175 摘要逐字 "removes the main obstacle to prove the ultraviolet stability"，R 操作 | **`[T]`** 【读摘要】 | DOI `10.1007/BF01238433`、`10.1007/BF01257412`（INSPIRE 逐字）；系列另见 `10.1007/BF01215223`(109,249)、`10.1007/BF01239022`(116,1)、`10.1007/BF01229380`(102,255, 3D) | 用 `CMP` **122** 355 替换仓内"Balaban 到不了大场"的表述（`G1-L3.md:16` 等）；`Z`→矩 三条缺口仍按 `HS-L3-8.md:56–62` 记 |
| **0.3** | **MRS `CMP` **155** (1993) 325–383 在 UV 轴上但不在矩轴上**：摘要逐字 "basis for a rigorous construction of the Schwinger functions … with a fixed infrared cutoff but no ultraviolet cutoff, in a regularized axial gauge"、"trivial topological sector"、纯 `SU(2)` | **`[T]`** 【读摘要】 | DOI `10.1007/BF02097397`（Crossref：页 325–383、issue 2、13 refs；INSPIRE abstract 逐字） | 只当"UV 极限在固定规范下走到过"的先例记；**不当** `L3^{UV}` 供给（无 gauge-invariant 矩族、无 OS） |
| **0.4** | **SZZ／Cao–Chatterjee／Chatterjee／Nissim 四条全不在 `a` 轴上**：SZZ 固定 `a`、`\|β\|<1/16(d−1)`（`SU(N)`，任意 `d>1`）；Cao–Chatterjee **3D** 且条件于 "provided that they indeed behave like the GFF"；Chatterjee 4D 要 `g=O(ε^{50d})→0`、极限是**有质量高斯**、"non-Gaussian … remains open"；Nissim 是 't Hooft 区固定 `a` 的隙/无穷体积 | **`[T]`** 【读摘要】 | arXiv:2204.12737 → DOI `10.1007/s00220-022-04609-1`（`CMP` **400**(2), 805–851）；arXiv:2111.12813 → DOI `10.1007/s00220-023-04870-y`（`CMP` **405**, art. 3）；arXiv:2401.10507 → DOI `10.2140/pmp.2026.7.339`；arXiv:2510.22788 | SZZ 仅作"固定 `a` 的强耦合参照"；Chatterjee 4D 只买"`a→0` 极限存在先例" |
| **0.5** | **"矩在格距上一致有界"这句已发表的话，只出现在 2D**：Chandra–Chevyrev 摘要逐字 "moments … are bounded uniformly in the lattice spacing"，对象是 **2D** Abelian Higgs 的 gauge field marginal、**固定规范**；He (2026) 摘要逐字 "tightness … as the lattice spacing tends to zero"、**2D** Abelian–Higgs Coulomb 规范 | **`[T]`** 【读摘要】 | arXiv:2207.05443 → DOI `10.1007/s00220-024-05008-4`（`CMP` **405**, art. 143）；arXiv:2611.16875 | 这两条是 `(U1)` 的**方法模板**（Besov 矩 + 抗磁不等式 / FBSDE + 高低频分解），维度是硬边界 |
| **0.6** | **★ 一条必须记账的强声明**：Faizal–Shabir，`Fortsch. Phys.` **74** (2026) 4, e70097，题 "Reflection-Positive Construction of a Four-Dimensional SU(N) Yang–Mills Theory With Mass Gap and Confinement"，摘要声称 "gauge-covariant multiscale analysis … interlaced inequalities whose small defects can be summed" ⇒ 连续极限 + "unique limiting measure"。**但摘要层无一个显式矩界**（无 `C_n`、无量词序、无 `a₀`），且被引 0 次 | **`[NG]`（列为 R1 待核声明，不列为供给）** 【读摘要】 | DOI `10.1002/prop.70097`（Crossref：2025-10-09 收、2026-03-17 受、2026-04-20 发；36 refs，含 Balaban 109/116、Mitoma 紧性、Billingsley） | **R1**：开正文核它是否真写出 `\|S^{(n)}_{L,a}(f)\|≤C_n(L)∏‖f_i‖`。若写出 ⇒ 本文件 0.1 翻成 `[T]`；若只写测度级紧性 ⇒ 仍 `[NG]` |
| **0.7** | **Bałaban–Jaffe 无联合定理**：两人的联合产出是两篇 NATO ASI 会议文集（`Constructive Gauge Theory`, 1986, 207–263；`Exact Renormalization Group for Gauge Theories`, 1984, 79–103），非新估计 | **`[T]`** 【读摘要】(Crossref 元数据) | DOI `10.1007/978-1-4757-0363-4_5`、`10.1007/978-1-4757-0280-4_4`；INSPIRE `a Balaban and a Jaffe` → 5 条，无期刊论文 | brief 里"Bałaban–Jaffe"这条线**结案**：它不供给任何一格 |

**一句话。** `L3^{UV}` 在已发表文献里**是空位**（`[NG]`），但空位两侧都已被踩过：左侧 Balaban 到 `Z`-一致，右侧 Chandra–Chevyrev 到 2D 矩一致；中间"4D + 纯 YM + 规范不变 + smeared 矩 + `a`-一致"这一格没人填。**唯一可能推翻这句话的是 R1（`10.1002/prop.70097`），下一手就是开它的正文。**

---

## §1 逐条核查（全部【读摘要】；0 篇读正文）

| # | 文献 | ID / DOI | 它的一致挂在什么上 | 到 `L3^{UV}` 的哪一格 | 判 |
|---|---|---|---|---|---|
| B1 | Balaban, `CMP` **109** (1987) 249–301 | `10.1007/BF01215223` | 小场近似下的有效作用量团簇展开、`β`-函数 | 上游（小场） | 固定截断外的准备 |
| B2 | Balaban, `CMP` **116** (1988) 1–22 | `10.1007/BF01239022` | 系列续篇（Crossref 只给卷页，摘要未取） | 上游 | 元数据 |
| B3 | Balaban, *Large field renorm. I*, `CMP` **122** (1989) 175–202 | `10.1007/BF01257412` | `R` 操作：大场区的**绝对值压制** | `(U2)` 上游 | `a`-一致，但非矩 |
| B4 | Balaban, *Large field renorm. II*, `CMP` **122** (1989) 355–392 | `10.1007/BF01238433` | **4D 纯规范场 UV 稳定的完成**；对象 = `Z`/有效作用量 | **`Z`-一致 ≠ 矩一致** | `[T]`（书目+摘要），**不是** `(U1)` |
| B5 | Bałaban–Jaffe（两篇 NATO ASI） | 见 0.7 | 会议综述 | 无 | `[T]` 无联合定理 |
| M1 | Magnen–Rivasseau–Sénéor, `CMP` **155** (1993) 325–383 | `10.1007/BF02097397` | 无 UV cutoff、固定 IR；轴规范；`SU(2)`；平凡拓扑；Slavnov 检验 | UV 轴走到过，**规范固定** | `[T]`，非 `(U1)` |
| S1 | Shen–Zhu–Zhu, `CMP` **400** (2023) 805–851 | `10.1007/s00220-022-04609-1` | 固定 `a` 的 Langevin/LSI/PI + 强耦合 `\|β\|<1/16(d−1)` | 无 `a` 轴 | `[T]`，固定截断 |
| C1 | Cao–Chatterjee, `CMP` **405** (2024) art. 3 | `10.1007/s00220-023-04870-y` | **3D** 状态空间，条件于 GFF 行为 | 无 | `[T]`，维度不符 |
| C2 | Chatterjee, `Prob. Math. Phys.` **7** (2026) 339–381 | `10.2140/pmp.2026.7.339` | 4D `a→0`，但 `g=O(ε^{50d})→0`、Higgs、极限有质量高斯 | 只买"极限存在" | `[T]`，`g→0` |
| C3 | Chatterjee, *Leading term of the YM free energy*, `JFA` **271** (2016) 2944 | arXiv:1602.01222 | 3D `a→0` 自由能；4D 仅弱耦合 | `Z`/自由能级 | `[T]`，非矩 |
| N1 | Chandra–Chevyrev, `CMP` **405** (2024) art. 143 | `10.1007/s00220-024-05008-4` | **2D**、固定规范 marginal、Besov 矩 `a`-一致 | `(U1)` 的**唯一已发表实例（2D）** | `[T]`，维度不符 |
| N2 | He, arXiv:2611.16875 (2026-11) | arXiv:2611.16875 | **2D** Abelian–Higgs Coulomb 规范、FBSDE、`a→0` 紧性 | 同上（2D） | `[T]`【读摘要】 |
| N3 | Nissim, arXiv:2510.22788 (2025-10) | arXiv:2510.22788 | `U(N)` 't Hooft 区、固定 `a` 的隙/无穷体积 | 无 `a` 轴 | `[T]` |
| N4 | Brydges–Fröhlich–Seiler I–III | `10.1016/0003-4916(79)90098-8`；`10.1007/BF01197918`；`10.1007/BF01208500` | 2D/3D 构造 + 格点逼近收敛 | 先例族，非 4D | 元数据（经 Cao–Chatterjee refs 核卷页） |
| N5 | Dimock, *RG according to Balaban* I–III | arXiv:1108.1335 / 1212.5562 / 1304.0705 | 3D `φ⁴` 的 Balaban 复述，末篇 "complete … a stability bound" | 说明"UV 稳定 = `Z` 界"的语义 | `[T]`【读摘要】 |
| R1 | Faizal–Shabir, `Fortsch. Phys.` **74** (2026) e70097 | `10.1002/prop.70097` | 声称 4D SU(N) 连续构造（多尺度 + 交错不等式） | 可能覆盖 `(U1)`，摘要未给算式 | **待核**，见 0.6 |
| X1 | arXiv:2506.00284（"4D SU(3) 构造 + 质量隙"） | arXiv:2506.00284 | **已被 arXiv 撤下**（"withdrawn … does not meet research content quality standards"） | 无 | 记号：勿引 |

**一行判**：`a`-一致**有**（B3/B4，`Z` 级）；矩一致**有**（N1/N2，2D）；**两者相交的那一格 = 0 篇**。R1 是唯一可能填它的声明。

---

## §2 空位判定

### 2.1 判定：**空位（`[NG]`）**，但不是"没人想过"

三条独立证据指向同一格为空：

1. **上游有、下游无**：Balaban 的 `Z`-一致（B4）与 Minlos 型紧性机器（N5 说明其语义）都在，但**没有任何一条已发表结果**把 `Z`-界升成 `C_n(L)∏‖f_i‖` 型的 smeared 矩界。这正是 `HQ1-L3.md:56–62` 拆出的三条缺口，本轮未找到第四条路。
2. **方法在、维度不在**：`a`-一致矩界的**技术**已在 2D 走通两遍（N1 Besov + 抗磁；N2 FBSDE + 高低频），且两篇摘要都用了同一句话级表述。4D 无对应物 —— 与 `G1-UNIFORM-ZA.md:9`（"No no-go blocks this"）不冲突：**空位 ≠ 有 no-go**。
3. **检索面到了摘要层全覆盖**：§3 的 18 条式覆盖 arXiv 全库摘要字段 + Crossref + INSPIRE；对 `abs:"uniform in the lattice spacing"`（命中 1，且是拓扑磁化率数值文）、`abs:"moment bounds" AND abs:"lattice spacing"`（命中 1，2D）两条**最贴题的式**，4D 纯 YM 命中均为 **0**。

### 2.2 与起点三句的关系

- `G1-L3.md:21` 逐字 "no published result closes L3. It is the frontier, not a cited lemma." ⇒ **本轮独立复核后维持**，并用 11 条具名 ID 补上书目。
- `G1-UNIFORM-ZA.md:9` 逐字 "No no-go blocks this" ⇒ **不冲突**：本轮判的是"文献里没人做"，不是"做不了"。
- `HS-L3-8.md:14` 的 `0.1 [NG]` 与本文件 0.1 **同判**；本轮新增的是 B5（Bałaban–Jaffe 结案）、N2（He 2026）、N5（Dimock 语义）、**R1（Faizal–Shabir）** 四条它没有的条目。

### 2.3 诚实边界（4 条）

1. **R1 是本判词唯一的活口**：0.6 的 `[NG]` 建立在"摘要未给显式矩界"上，不是"正文没有"。R1 正文若写出 `(U1)` 型不等式，本文件 0.1 立即作废并翻 `[T]`。
2. **0 篇读正文**：B1–B4、M1 的"没给矩界"是**摘要量词形状**的判读；其中 B4 的 Theorem 1 假设（counterterm 是否沿 `β(a)` 递推、体积是否固定）**未逐字核**。
3. **zbMATH 取不到**：`zbmath.org/?q=an:0704.58060` → HTTP 403（Cloudflare 挑战）⇒ 记 `TOOL-FAIL`，**不作**任何"评语未提及"的推断。
4. **无 QED**：本文件不产估计、不升任何 `[T]` 的范围；`[NG]` 读作"该格在已发表文献中 0 供给"，读作"不可证"是误读。

---

## §3 检索四要素（`AGENTS §1`）

| # | 工具 | 状态 | 实际检索式（逐条） | 命中 |
|---|---|---|---|---|
| 1 | `web_fetch`→arXiv API | **OK**（HTTP 200，Atom） | ①`abs:"uniform in the lattice spacing"` ②`abs:"Yang-Mills" AND abs:"tightness"` ③`abs:"Yang-Mills" AND abs:"uniform bounds"` ④`abs:"pure Yang-Mills" AND abs:"continuum limit"` ⑤`abs:"ultraviolet stability"` ⑥`abs:"lattice gauge" AND abs:"uniform" AND abs:"spacing"` ⑦`abs:"moment bounds" AND abs:"lattice spacing"` ⑧`abs:"gauge" AND abs:"tightness" AND abs:"lattice"` ⑨`abs:"Yang-Mills" AND abs:"mass gap"` | ①**1**（2510.08006，假阳性：拓扑磁化率数值）②**7**（0 条 4D 纯 YM 紧性）③**0** ④**7**（0 条构造性）⑤**17**（0 条 4D 纯 YM）⑥**9**（1 条相关：2611.16875，2D）⑦**1**（2611.16875，2D）⑧**108**（"tight-binding"噪声；1 条相关，2D）⑨**186**（扫前 25：得 R1、X1、N3） |
| 2 | `web_fetch`→arXiv API（作者/ID） | **OK** | ⑩`au:Balaban_T` ⑪`au:Dimock_J AND abs:Balaban` ⑫`au:Chatterjee_S AND abs:"Yang-Mills"` ⑬`au:Chatterjee AND abs:"lattice gauge"` ⑭`au:Faizal AND abs:"Yang-Mills"` ⑮`id_list=2204.12737,2111.12813,2401.10507,2207.05443` | ⑩**3**（全为化学家 A.T. Balaban，物理 Balaban 无 arXiv）⑪**8** ⑫**0**（字段式失败→改 ⑬）⑬**7** ⑭**11**（含 R1）⑮**4**（四条摘要逐字到手） |
| 3 | `web_fetch`→Crossref | **OK**（JSON） | ⑯`works/10.1007/BF02097397`、`…/s00220-022-04609-1`、`…/s00220-023-04870-y`、`…/prop.70097` ⑰`works?query.author=Balaban+Jaffe&rows=8` ⑱`works?query.author=Balaban&query.bibliographic=renormalization group lattice gauge field theories&rows=20` | ⑯**1/1/1/1**（卷-页-日期全核）⑰**8** 返（4 条相关：B5 两篇 + Imbrie 两篇）⑱**20** 返（Balaban 系列卷页全表） |
| 4 | `web_fetch`→INSPIRE | **OK**（JSON） | ⑲`q=a Balaban and a Jaffe` ⑳`q=a Balaban and t large field` ㉑`q=a Magnen and a Rivasseau and a Seneor and t infrared cutoff` | ⑲**5**（0 条期刊联合论文）⑳**3**（**B3/B4 摘要逐字**）㉑**1**（M1 摘要逐字） |
| 5 | `read`/`glob`（仓内） | **OK** | `HS-L3-8.md`(全 102)、`HQ1-L3.md`(1–80)、`G1-UNIFORM-ZA.md`(1–20)、`G1-L3.md`(1–28) | 4/4 返回正文 |
| 6 | **失败标记** | ①`zbmath.org/?q=an%3A0704.58060` → **HTTP 403**（Cloudflare "Just a moment…"）= **`TOOL-FAIL`**，故 B4 的第三方评语本轮未取（`HS-L3-8.md:16` 的评语仍有效，本文件不重复引用为独立证据）②`au:Chatterjee_S AND abs:"Yang-Mills"` → **0**，属**查询字段式产物**（改 ⑬ 后 7 条），不是覆盖缺口 | — | `TOOL-FAIL` **1**；查询式失败 **1**；空响应 **0** |

**检索面边界（诚实）**：以上 18 条式覆盖 **arXiv 摘要字段全文库 + Crossref + INSPIRE**，均为**摘要层**。所以 0.1 的 `[NG]` 读作：**在这 18 条式下、摘要层内，4D 纯 YM 的 gauge-invariant `a`-一致矩界与 `a→0` 紧性 0 条供给**。它不证明"不存在"，只把每一条候选钉到它实际所在的那根轴上；唯一未被钉死的是 **R1**（正文未开）。

*HT6-LIT · 只新建本文件；既有文件原样；未 `git commit`；零 MC；未动 `docs/position-current.html`；≤130 行。*

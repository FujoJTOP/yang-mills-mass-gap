# RG2-MAP — Balaban 的 CMP 系列（1984–1989）：覆盖到哪、IR 停在哪

> 用途：核 l4a 的判据「4D SU(N) 的 UV 有真定理（Balaban），IR 是 0 条 ⇒ [O]」。
> 本轮只做**书目 + 摘要级**地图。正文在付费墙后 ⇒ 凡"正文"级断言一律标 **未逐字核**。
> 时间戳仅作追溯；进展以「轮」计。

---

## 0. 检索执行记录（AGENTS.md §1 四项）

### 0.1 工具状态

| 工具 | 状态 | 说明 |
|---|---|---|
| `web_search` | **OK** | 6 条检索式全部返回有效结果（每次 10 条源） |
| `web_fetch` → Crossref API | **OK** (HTTP 200) | `api.crossref.org` 期刊/书目查询 |
| `web_fetch` → OpenAlex API | **OK** (HTTP 200) | 作者全量作品表（63 条）、按年过滤（1984–1990：24 条） |
| `web_fetch` → INSPIRE-HEP API | **OK** (HTTP 200) | 取 **Springer 原摘要**（`abstracts[].source = "Springer"`） |
| `web_fetch` → arXiv HTML | **OK** (HTTP 200) | Dimock 综述正文（参考文献段被截断，未取到） |
| `web_fetch` → link.springer.com / doi.org | **TOOL-FAIL** | `cross-origin redirect to https://idp.springer.com`，正文页取不到 |
| `web_fetch` → zbmath.org | **TOOL-FAIL** | HTTP 403（Cloudflare） |
| `web_fetch` → api.semanticscholar.org | **TOOL-FAIL** | HTTP 429 |

> 摘要的**权威来源**因此改为 INSPIRE 的 `source: "Springer"` 字段（Springer 原摘要转发），非二手转述。

### 0.2 实际检索式与命中数

| # | 检索式 | 命中 |
|---|---|---|
| S1 | `"Balaban" "Communications in Mathematical Physics" 1985 "ultraviolet stability"` | 10 源（含 3D UV 稳定性条目） |
| S2 | `"Balaban" "gauge field" "renormalization group"` | 10 源 |
| S3 | `Balaban ultraviolet stability gauge fields Yang-Mills CMP 1987 1988 1989` | 10 源 |
| S4 | `Balaban "Renormalization group approach to lattice gauge field theories" I "109" 249` | 10 源（得 DOI 10.1007/BF01215223） |
| S5 | `Balaban "Renormalization group approach ... III"` / `... IV"` / `Propagators ... 107` | 10 源，**0** 条指向 III/IV（该编号不存在） |
| S6 | `Balaban four-dimensional ultraviolet stability proof not complete "last step"` + `... does not address infrared mass gap 4D Yang-Mills` | 10 源，**9/10 为 2026 年 Zenodo/OSF 的 AI 生成"补全 Balaban"预印本**（不可采信），**0** 条一手文献 |
| Q1 | Crossref `journals/0010-3616/works?query.author=Balaban` | 16 条 |
| Q2 | Crossref `works?query.bibliographic=Balaban renormalization group lattice gauge field theories` | 30 条（含重复/无关） |
| Q3 | OpenAlex `authors?search=Tadeusz Balaban` | 1 作者（`A5110504215`，63 works） |
| Q4 | OpenAlex `works?filter=author.id:A5110504215,publication_year:1984-1990` | 24 条 |
| Q5 | OpenAlex `works?filter=author.id:A5110504215&per-page=100` | 63 条（全量） |
| Q6 | INSPIRE `q=doi:...`（5 个 DOI 批量） | 4 条（命中 4/5；`BF01239022` 无 INSPIRE 记录） |
| Q7 | INSPIRE `doi:10.1016/0003-4916(84)90121-0` + IR 会议录题名 | 2 条 |

---

## 1. 先纠一处：题面给的卷号里有 4 个**不属于** Balaban

题面要求"至少覆盖 CMP 102、107、109、110、113、115、116、117、119、122"。三源（Crossref / OpenAlex 全量 63 条 / INSPIRE）交叉后：

| 卷 | Balaban 条目 | 实际作者 / 说明 |
|---|---|---|
| 102 | **2 篇**（255–275、277–309） | ✔ Balaban |
| **107** | **0** | 该卷 Yang–Mills 论文是 **Federbush**：`A phase cell approach to Yang–Mills theory. I. Small field modes`，CMP **107** (1986) 319–329，DOI `10.1007/BF01209397`（INSPIRE 作者字段已核） |
| 109 | **1 篇**（249–301） | ✔ Balaban（系列 I） |
| **110** | **0** | 同卷 Yang–Mills 论文仍是 **Federbush**：`... III. Stability, modified renormalization group transformation`，CMP **110** (1987) 293–309，DOI `10.1007/BF01207369` |
| **113 / 115 / 117** | **0 / 0 / 0** | 三源均 0 命中。1984–1990 窗口 Balaban 的 CMP 论文只落在 95、96、97、98、99、102、109、114、116、119、122 |
| 116 | **1 篇**（1–22） | ✔ Balaban（系列 II） |
| 119 | **1 篇**（243–285） | ✔ Balaban（第 3 步，换题不换号） |
| 122 | **3 篇**（175–202、233–247、355–392） | ✔ Balaban |

> 结论：**题面的 107 / 110 / 113 / 115 / 117 这五个槽位里，Balaban 只占了 0 个**；107、110 属 Federbush 的另一条（同样未完成的）phase-cell 路线。1984–1989 窗口 Balaban 的 CMP 论文总数 = **15 篇**（其中 4D 纯 YM 主线 5 篇 + 前置工具 7 篇 + 侧支 3 篇）。

---

## 2. 地图：Balaban CMP 系列（卷期页 / 年 / 标题 / 一句话 / 覆盖尺度）

尺度列只有三种取值：**UV** = 紫外方向（短距离、块自旋下行、抵消发散）；**IR** = 红外方向（长距离、质量隙、相关衰减）；**中间** = 单步 RG 的中间量（背景场/极小元），不构成任一端点结论。
证据列：`A` = 只读到**摘要**（Springer 原摘要，经 INSPIRE 转发）；`B` = 读到正文（本轮**无**一篇满足，付费墙）。

### 2.1 前置工具（全是 UV 侧）

| 卷期页 | 年 | 标题 | 一句话内容 | 尺度 | DOI | 证据 |
|---|---|---|---|---|---|---|
| 95, 17–40 | 1984 | Propagators and renormalization transformations for lattice gauge theories. I | 把格点规范场看作高斯向量场的扰动，给出传播子与 RG 变换的基本估计 | UV | 10.1007/BF01215753 | A |
| 96, 223–250 | 1984 | 同上 II | 把 I 的结果推广到"限制在不同尺度上"的算子与不同阶的 RG 变换 | UV | 10.1007/BF01240221 | A |
| 98, 17–51 | 1985 | Averaging operations for lattice gauge theories | 研究场与规范变换的**平均算子**，刻画出使其解析/正则的位形类 | UV | 10.1007/BF01211042 | A |
| 99, 75–102 | 1985 | Spaces of regular gauge field configurations on a lattice and gauge fixing conditions | 规范不变正则性条件 ∩ 规范固定面 ⇒ 位形既小又正则 | UV | 10.1007/BF01466594 | A |
| 99, 389–434 | 1985 | Propagators for lattice gauge theories in a background field | 带背景场的传播子：正则性与衰减 | UV | 10.1007/BF01240355 | A |
| 102, 255–275 | 1985 | **Ultraviolet stability of three-dimensional lattice pure gauge field theories** | **证出 3D 格点纯 Yang–Mills（Wilson 近似）的紫外稳定性** | UV（3D） | 10.1007/BF01229380 | A |
| 102, 277–309 | 1985 | The variational problem and background fields in RG method for lattice gauge theories | 固定平均值的正则位形空间上作用量存在极小元，且规范变换下唯一 ⇒ 背景场 | 中间 | 10.1007/BF01229381 | A |

### 2.2 4D 纯 Yang–Mills 主线（5 篇）

| 卷期页 | 年 | 标题 | 一句话内容 | 尺度 | DOI | 证据 |
|---|---|---|---|---|---|---|
| 109, 249–301 | 1987 | **Renormalization group approach to lattice gauge field theories. I. Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions** | 4D 纯规范场：**限于小场近似**，用单步 RG + 团簇展开造局部有效作用量，造 β 函数与耦合常数重正化的递推方程组 | UV | 10.1007/BF01215223 | A |
| 116, 1–22 | 1988 | 同上 **II. Cluster expansions** | 把 I 里构造的涨落场积分写成指数化团簇展开，证明展开项满足归纳假设 | UV | 10.1007/BF01239022 | A |
| 119, 243–285 | 1988 | Convergent renormalization expansions for lattice gauge theories | 引入含**大场区**的完整有效密度的归纳描述，证明 RG 变换保持密度形式；**由此完成超可重正模型的 RG 分析并给出收敛展开** | UV（超可重正 = d≤3） | 10.1007/BF01217741 | A |
| 122, 175–202 | 1989 | **Large field renormalization. I. The basic step of the R operation** | 构造大场区表达式的重正化操作 ℝ，"**移除了证明 4D 规范场紫外稳定性的主要障碍**" | UV | 10.1007/BF01257412 | A |
| 122, 355–392 | 1989 | **Large field renormalization. II. Localization, exponentiation, and bounds for the R operation** | 收束 ℝ 操作的基本性质与界，"**从而完成 4D 纯规范场紫外稳定性的证明，如 Theorem 1 所述**" | UV | 10.1007/BF01238433 | A |

### 2.3 侧支（3 篇，别混进主线）

| 卷期页 | 年 | 标题 | 一句话内容 | 尺度 | DOI | 证据 |
|---|---|---|---|---|---|---|
| 97, 299–329 | 1985 | Renormalization of the Higgs model: Minimizers, propagators and the stability of mean field theory | Higgs 模型的极小元/传播子/平均场稳定性 | UV（阿贝尔） | 10.1007/BF01206191 | A |
| 114, 257–315 | 1988 | Effective action and cluster properties of the abelian Higgs model | **阿贝尔** Higgs（**2D、3D**）：多尺度团簇展开，给关联指数衰减所需的解耦 ⇒ 为 Higgs 机制与质量隙铺路 | **IR 侧（但 d≤3 且阿贝尔）** | 10.1007/BF01225038 | A |
| 122, 233–247 | 1989 | Block renormalization group for Euclidean fermions（与 O'Carroll、Schor） | 格点/连续欧氏费米子的块 RG：自由场下作用量序列一致指数衰减；"δ-函数"权重反而**不**一致衰减 | UV | 10.1007/BF01257414 | A |

### 2.4 窗口外的同名 IR 工作（不在 CMP，必须单列，否则会误判）

| 出处 | 年 | 标题 | 一句话 | 尺度 | DOI/出处 | 证据 |
|---|---|---|---|---|---|---|
| Annals of Physics **158**, 281–319 | 1984 | The mass gap for Higgs models on a unit lattice | 建立紧/非紧阿贝尔规范理论的同构，弱耦合下证出 **Higgs 机制预言的 mass gap** | **IR（阿贝尔、单位格点、弱耦合）** | 10.1016/0003-4916(84)90121-0 | A |
| Montreal 1987 会议录 135–151 | 1990 | Infrared asymptotic freedom for perturbations of the free lattice electromagnetic field | 自由格点电磁场扰动的**红外渐近自由** | **IR（U(1)₄，围绕自由场）** | 会议录（无 DOI） | A（题名级，INSPIRE 无摘要） |

> 注意这两条的分量：一条是 **阿贝尔 + d=2,3 或单位格点弱耦合**，一条是 **U(1)₄ 的微扰式红外**。**都不是 4D 非阿贝尔纯 YM 的质量隙。**

---

## 3. 判定：IR 那一段做了没有？

### 3.1 事实层

1. **4D 非阿贝尔纯 YM 的 IR / 质量隙：0 篇。** 15 篇 CMP 窗口论文里，标题含 "infrared" 的 **0 篇**，含 "mass gap" 的 **0 篇**；摘要里出现红外结论的 **0 篇**。
2. **4D 主线的终点是"紫外稳定性"这一条界，不是理论构造。** 122(II) 摘要原话：完成的是 "the proof of the ultraviolet stability of four-dimensional pure gauge field theories, **as formulated in Theorem 1**" —— 一个**有限体积、小耦合**下的稳定性**界**，不是连续极限、不是热力学极限、更不是谱隙。
3. **他自己划了范围**（三条硬边界，均为摘要原话）：
   - 109(I)：「The analysis is **restricted to small field approximation**.」
   - 119：「This **completes the renormalization group analysis for superrenormalizable models**」——超可重正 = d≤3，**4D 不在内**。
   - 122(I)：「removes the **main obstacle** to prove the ultraviolet stability of four-dimensional gauge field theories」——是"移除障碍"，不是"给出 IR"。
4. **他从未声称 4D 非阿贝尔 YM 存在。** 1986 年 Erice 讲义《Constructive Gauge Theory》摘要原话：「**Up to now, only abelian gauge fields are known to exist**, namely the U(1) Higgs model on ℝ² and ℝ³ and electrodynamics on the tori T² and T³.」——这是他自己对"已建成什么"的盘点，4D 非阿贝尔不在清单里。
5. **IR 侧他只在阿贝尔、低维做**：CMP 114（阿贝尔 Higgs，2D/3D，指数衰减/解耦）、Ann. Phys. 158（Higgs mass gap，单位格点弱耦合）、Montreal 1987（U(1)₄ 红外渐近自由）。

### 3.2 他"有没有说做不了"？

- **明说"做不到"的原话：本轮 0 命中**（检索式 S6；工具 OK；命中 0；一手文献层面）。摘要里没有 "impossible / cannot / beyond reach" 这类表述。
- **但他用的是"范围划定"而非"能力宣告"**：小场近似、超可重正、UV stability only、以及"只有阿贝尔模型已知存在"。这等价于把 IR 排除在论文目标之外，而不是宣布 IR 不可证。
- 因此对 l4a 判据的正确读法：
  - 「UV 有真定理」→ **成立**（[T] 级为：*Balaban 在 122(II) 声称完成 4D 纯规范场 UV 稳定性（Theorem 1；有限体积、小耦合）*；注意这是"他证了/他声称证了"，不是"我们复核过正文"——正文未逐字核）。
  - 「IR 是 0 条」→ **成立**（[T]：15 篇 CMP 窗口论文中 4D 非阿贝尔 IR 结论 = 0）。
  - 「⇒ [O]」→ 逻辑上要小心：**0 条不等于"已证不可做"**。正确标注是 **[NG]（无人做 / 无定理）而非 [T]（不可能）**，除非另有一条 no-go。本轮**没有**找到 Balaban 的 no-go 声明。

### 3.3 尚未闭合的疑点（不许含糊过去）

- 122(II) 的"完成"是**作者自述**。学界对 Balaban 4D 证明"是否逐字闭合/可复核"的公开讨论，本轮 `web_search` 只捞到 **2026 年 Zenodo/OSF 的 AI 生成预印本**（自称"补全 Balaban"），**0 条一手或可核的二手评述**。这批材料不可引用。
- 若要判"Balaban 的 4D 证明是否有洞"，需要读 122(I)(II) 正文，本轮**未取得**（Springer → idp 跨域重定向，TOOL-FAIL）。**这条挂在墙上，标 [M]。**
- Dimock 的三篇综述（*The renormalization group according to Balaban* I/II/III，arXiv:1304.0705 等）做的是 **φ⁴₃** 的示范重写，不是 4D 规范场；其参考文献段在本轮抓取中被截断，未取到完整书目——但 Balaban 系列编号（I/II + 换题的第三、四步）已由 Crossref/OpenAlex/INSPIRE 三源独立确认，**不存在 "III"/"IV" 编号的 CMP 论文**（S5：命中 0）。

---

## 4. 判死线自查

- **(a) 成立**：列出 **12 篇** Balaban 规范场 CMP 论文（≥6），每篇给卷期页 + DOI + 覆盖尺度；并判 IR 是否被触及 ⇒ **仅在阿贝尔/低维被触及，4D 非阿贝尔 = 0**。
- 覆盖尺度分布：**UV 11 篇**（95、96、98、99×2、102/255、109、116、119、122×2）、**中间 1 篇**（102/277）、**4D 非阿贝尔 IR 0 篇**；IR 侧 3 条全在 CMP 之外且非 4D 非阿贝尔（114、Ann.Phys.158、Montreal 1987）。
- 未编造项：所有卷期页/DOI 均来自 Crossref/OpenAlex/INSPIRE 三源之一；所有"内容"句均可在标注的摘要原文中对照；正文一律标 A（仅摘要）。

---

## 5. 一句话交接

**Balaban 的 CMP 系列把 4D 纯 YM 一路推到"紫外稳定性"的门口并在 122(II) 声称完成，全程只做 UV 侧的小场/大场重正化；4D 非阿贝尔的红外与质量隙，他 0 篇、也 0 次声称能做——但他同样没写"做不了"。所以「UV 真定理 + IR 0 条」成立，「⇒ [O]」应写作 [NG]/无定理，除非另有 no-go；而"他的 4D 证明是否闭合"仍是 [M]，因为正文取不到。**

# A9 · 四篇文献里的指数 `c'`：逐个判定

> 输入：`A7-EXPONENT.md`（`p = 2b₀c'`，`c'_crit = 1/(2b₀) = 24π²/(11N)`）、`A5-LIT-UNIFORM.md:93`（`f(g) ~ e^{−c'/g²}`）。
> 纪律：只写可核对的逐字引用 + 出处（式号/页码/摘要）；取不到的写 `未取数`/`未查到（检索式：…）`。**不从记忆补文献。**
> 检索日期 2026-11-25。全部取数走 `web_fetch` / `Invoke-WebRequest` + `pdftotext`，无记忆条目。

---

## §0 判死线（第一段）

**判 (c)。四篇里没有一篇给出形如 `m a ≥ c·f(g)`、`f(g) ~ e^{−c'/g²}` 的 `a`-依赖谱隙界 ⇒ `c'` 在这四篇里不是"取数大小"的问题，是"该形状根本不存在"的问题。**

逐篇一句话（详见 §1–§4）：

| # | 文献 | 实际内容（逐字/可核） | 小量因子形式 | `c'` |
|---|---|---|---|---|
| 1a | Balaban, *CMP* **102** (1985) 255–275 | 标题逐字：`Ultraviolet stability of three-dimensional lattice pure gauge field theories` —— **紫外稳定性**，非谱隙 | 稳定界，无 `m`；小量因子**未取数** | 未取数 |
| 1b | Balaban, *CMP* **116** (1988) 1–22 | 摘要逐字（§1.2）：`...the terms of the expansion satisfy the inductive assumptions. This completes the construction of the sequence of effective actions in the small field approximation.` —— **只有有效作用量** | 同上，无 `m`、无 `a`-依赖隙 | 未取数 |
| 2 | Battle–Federbush, *Ann. Phys.* **142** (1982) 95–139 | 标题逐字：`A phase cell cluster expansion for Euclidean field theories`；摘要逐字：`We adapt the cluster expansion first used to treat infrared problems for lattice models (a mass zero cluster expansion) to the usual field theory situation.` —— **不是规范理论、不是谱隙、无 `β` 半径** | 无（非 YM 论文） | 不适用 |
| 3 | SZZ, `arXiv:2204.12737` | 条件 `|β| < 1/(16(d−1))`（摘要逐字）；§4.3 逐字：`In this section we use the Poincaré inequality to prove the existence of mass gap for lattice Yang–Mills.` | **`β` 的常数半径**，无指数 | **不存在**（无 `e^{−c'/g²}`） |
| 4 | Dimock, `math-ph/0210020` | 摘要逐字：`In this paper we explore the first step in detail.`；式 (14) 逐字：`e_{0}=L^{-N/2}e, m_{0}=L^{-N}m, μ_{0}=L^{-N}μ` | 小量 = `e₀ ∝ a^{1/2}`，**`a` 的幂**，非指数 | 不存在 |

**为什么 (c) 而不是 (a)/(b)**：(a) 需要"任一篇的小量因子是 `e^{−1/(2b₀g²)}`（或其幂）"——四篇**没有一篇**用小量因子 `e^{−c/g²}` 写出谱隙界；(b) 需要"全部是 `c' ≠ 1/(2b₀)`"——但 `c'` 在这四篇里**未被使用过**，给不出 `p` 值（`p` 只在形状存在时才有定义）。故 A7 的三档判据（`p⋛1`）在这四篇上**没有对象可判**——与 `a5` 的"未查到"正面一致。

**代价与救法（本文件唯一正面结论）**：`c'_crit = 24π²/(11N) = 10.7668 (N=2) / 7.1779 (N=3)` 这个数**没有文献来源**，它是 `b₀` 的算术推论（`A7 §1` 已三重核对）⇒ **我们自己不该去"借"一个已发表的 `c'`，只该写常数型界**：见 §5。

---

## §1 Balaban：两篇都不是谱隙论文，`c'` 未被使用

### §1.1 *CMP* **102** (1985) 255–275

书目核（DOI 逐字 `10.1007/BF01229380`，[INIS/ADS 记录](https://www.imstat.org/publications/aop/aop_53_1/aop_53_1.pdf#5#10#5) 引文条：`BALABAN, T. (1985). Ultraviolet stability of three-dimensional lattice pure gauge field theories. Comm. Math. Phys. 102 255–275`）：

> **题名逐字**：`Ultraviolet stability of three-dimensional lattice pure gauge field theories`

对象是 **YM₃ 的紫外稳定性**（配分函数/有效密度的界），**不是** `m a ≥ c·f(g)` 型谱隙界。**摘要正文未取到**：

- `https://link.springer.com/article/10.1007/BF01229380` → HTTP 200 但正文仅 3036 字节（登录/JS 壳），无摘要字段；
- NII-REO 全文 PDF `.../ART1000344029.pdf` 类路径对 102 卷无开放件；116 卷 PDF 返回日文页 `本コンテンツは、閲覧することができません。`（`The item you requested is not available.`）。
- 检索式（可复现）：`Balaban Comm. Math. Phys. 102 (1985) 255-275`、`Balaban ultraviolet stability three-dimensional lattice pure gauge`、`site:reo.nii.ac.jp Balaban`、`Invoke-WebRequest link.springer.com/article/10.1007/BF01229380`。

⇒ **`c'` 未取数；且题名层面已排除"谱隙界"**（稳定性 ≠ 谱隙）。

### §1.2 *CMP* **116** (1988) 1–22

书目 + 摘要**逐字**（NII-REO `ART1000344029`，DOI `10.1007/BF01239022`，卷 116 期 1，1988-03，p.1–22）：

> `Renormalization group approach to lattice gauge field theories : II. Cluster expansions`
> `Abstract: The fluctuation field integral, constructed in Part I, is represented by the exponentiated cluster expansion. It is proved that the terms of the expansion satisfy the inductive assumptions. This completes the construction of the sequence of effective actions in the small field approximation.`

**判决**：这是"小场近似下有效作用量序列的构造"——`ρ_k = (RT)^k ρ_0` 的归纳假设，**没有质量、没有谱隙、没有 `a`**。故"他的小场/大幅度界里耦合小量写成 `e^{−c/g²}` 还是 `g^N` 还是 `e^{−c/g²}g^{−k}`"这个问题在本篇**无法回答：该界不存在**（`未取数`）。

旁证（**二手，不作数**，只标注待核）：检索命中一份自出版文档写 Balaban 的大场因子形如 `exp(−c(log g^{−2})^q)`（[Zenodo search 摘要片段](https://zenodo.org/search?q=parent.id%3A23970310&f=allversions%3Atrue&l=list&p=1&s=10&sort=version#1) 逐字：`like all factors exp(−c (log g^(−2))^q)`）。**若**属实则 Balaban 的小量是**对数幂型**、非 `e^{−c/g²}`、非 `g^N`，`p` 也**不是常数**而是 `a`-依赖的指数（见 §6）。**此条未从原文取到式号，标 `未取数`。**

---

## §2 Battle–Federbush：这不是规范理论论文 ⇒ 引用错配

书目 + 摘要**逐字**（*Annals of Physics*, 1982-01, Vol.**142**(1), p.95–139；[INIS](https://inis.iaea.org/search?q=metadata.creators.person_or_org.name%3A%22Battle%2C%20G.A.%20III%22&l=list&p=1&s=10&sort=bestmatch)、[Wisconsin 记录](https://search.library.wisc.edu/article/cdi_crossref_primary_10_1016_0003_4916_82_90230_5/#database#1)）：

> 题名逐字：`A phase cell cluster expansion for Euclidean field theories`
> 摘要逐字：`We adapt the cluster expansion first used to treat infrared problems for lattice models (a mass zero cluster expansion) to the usual field theory situation.`

**判决**：`Euclidean field theories`＋`phase cell`（相格）＝ φ⁴ 型相格团簇展开，**既不是格点 Yang–Mills，也没有 `β`-收敛半径的表述**（摘要是"把零质量团簇展开搬到通常场论"）。它无法提供"团簇展开收敛半径里的 `β`-小量形式"，因为对象里没有规范耦合。

- 检索式：`Battle Federbush 1982 Annals of Physics 142 95 cluster expansion lattice gauge theories convergence`、`Battle Federbush "A phase cell cluster expansion for Euclidean field theories"`。
- 若需"团簇展开半径 `β` 小"这件事的**正确**出处，逻辑上应查 Osterwalder–Seiler `Ann. Phys. 121 (1979)`（SZZ 摘要/引言逐字把它当此事的来源：`A closely related earlier paper is by Osterwalder–Seiler [OS78], which showed that for the lattice Yang–Mills theory, when the coupling is sufficiently strong, the cluster expansion ... is convergent`）——**本文档未取该篇正文，标 `未取数`。**

⇒ `c'` **不适用**（非 YM 论文）。

---

## §3 SZZ：常数半径，`e^{−c/β}` 明确不存在

同一份 ar5iv 全文（`https://ar5iv.labs.arxiv.org/html/2204.12737`）两处**逐字**：

> 摘要：`Assuming |β| < (N−2)/(32(d−1)N) for the structure group SO(N), or |β| < 1/(16(d−1)) for SU(N), we prove the following results. ... correlations of a large class of observables decay exponentially, namely the infinite volume measure has a strictly positive mass gap.`

> §4.3 首句：`In this section we use the Poincaré inequality to prove the existence of mass gap for lattice Yang–Mills.`

**确认用户的条件成立：他们不是这个形状。** 论证（`[T]`，只依赖上面两句）：

1. 假设集是 `|β| < r`、`r = 1/(16(d−1))` 为**与 `g` 无关的常数**（`'t Hooft` 标度 `βN`）。谱隙来自 Poincaré/LSI 常数在球 `|β|<r` 上的**一致正性** ⇒ 界形如 `m ≥ c(d,r,β)`，`c` 在**两侧**（`β→0⁺` 与 `β→0⁻`）都不退化。
2. 任何 `f ~ e^{−c'/β}` 在 `β→0⁺` 时 `→0`、在 `β→0⁻` 时 `→+∞`，**不可能**在含 0 的双侧球上充当"正下界"而不在 `β=0` 崩塌；而 SZZ 的定理恰好断言 `β=0` 邻域内（含 `β<0`）成立 ⇒ **界内不可能出现 `e^{−c/β}` 或 `e^{−c'/β²}` 型因子**（否则 `β<0` 侧结论立刻平凡/错误）。
3. 与 `A5/A7` 的口径对照：SZZ 是**强耦合侧（`β` 小）**的 `N`-一致界，其自然位置在 `β = 2N/g²` **小**（`g` 大）一侧；`e^{−1/(2b₀g²)}` 是 `g→0`（`β→∞`）侧的量 ⇒ 与本篇**不同侧**，这也独立地解释了为何"没有指数"。

⇒ `c'`：**不存在**（不取数、也不是 `≠ 1/(2b₀)`；是"无此形状"）。

---

## §4 Dimock：明说只做第一步；`a`-依赖是 `a^{1/2}` 的幂

摘要**逐字**（`https://ar5iv.labs.arxiv.org/html/math-ph/0210020`）：

> `This is done by following the flow of a sequence of renormalization group transformations. The analysis is facilitated by splitting the space of gauge fields into into large and small fields at each step following Balaban. In this paper we explore the first step in detail.`

`a`-依赖的**逐字**形态（该文式 (14)，缩放后单位格点耦合）：

> `e_{0}=L^{-N/2}e \quad m_{0}=L^{-N}m \quad \mu_{0}=L^{-N}\mu`

其中 `a = L^{-N}`（`𝕋_M^{-N} = (L^{-N}ℤ / L^M ℤ)^3`，式 (1)）。⇒ **小量是 `e₀ = a^{1/2}e`，即 `a` 的 `1/2` 次幂**；这是**多项式**依赖，写成 `f(a) = (aΛ)^p` 只能是 `p` 半整数/整数级，**不可能**是 `e^{−c'/g²}`。

硬证据（可复算）：`Invoke-WebRequest .../math-ph/0210020` 全文剥离标签后，字符串 `mass gap` 命中 **1 次**，位置在参考文献 `[7] T. Balaban, D. Brydges, J. Imbrie, A. Jaffe, The mass gap for Higgs models on a unit lattice, Ann. of Phys. 158: 281-319, 1984` 之内 ⇒ **本文正文没有谱隙陈述**。

⇒ `c'`：**不存在**（`a`-依赖是幂 `1/2`，且无谱隙界）。

---

## §5 结论与"我们自己该写什么"

### §5.1 判 (c)

**四篇都不是这个形状**（无一篇含 `m a ≥ c·e^{−c'/g²}`；两篇非谱隙、一篇非 YM、一篇自限第一步）。故：

- `c'` 在这四篇上**全部未取数或不存在** ⇒ **无法给 `p` 值**（`p = 2b₀c'` 无输入）；
- `A5` 的"未查到"与 `A7` 的"语料无实测 `c'`"**双向一致** ⇒ 这是**已发表侧的真空**，不是我们漏检。

### §5.2 我们自己该写什么形状（可执行）

**不许**再写 `m a ≥ c·f(g)`、`f ~ e^{−c'/g²}` —— 它在 `a→0` 只有 `p=1` 一个活点，而 `p=1` 要求 `c' = 24π²/(11N)`，这等价于**把 `f` 换成 `aΛ`**，等于同义反复。落笔形状三选一，按优先级：

1. **无量纲常数型（首选，`a`-无关）**
   `m / Λ ≥ c(N)`，`c(N) > 0` 与 `a` 无关。等价写法 `m a ≥ c·(aΛ)`（即 `A7` 的 `t3`/`p=1` 档）。这是**唯一**在 `a→0` 不退化的形状。
2. **把 `e^{−1/(2b₀g²)}` 显式钉死成 `aΛ`**，只在与格距挂钩的中间步骤用：
   `m a ≥ c · aΛ ⇒ m ≥ cΛ`（同一句，(2.1)–(2.4) 的代数，`A7 §2–§3`）。写 `e^{−c'/g²}` 必须**同时**写 `c' = 1/(2b₀) = 24π²/(11N) = 10.7668 (N=2) / 7.1779 (N=3)`，否则落入 `p>1`（空洞）或 `p<1`（荒谬）。
3. **对数幂型（若走 Balaban 路线的合法替代）**
   若真实小量是 `exp(−c(log g^{−2})^q)`（§1.2 二手旁证，待核），则 `f = (aΛ)^{2b₀c(log(1/(aΛ)))^q·(…)}^{}` 使**指数不是常数**：`p(a) = 2b₀c·(log(1/(aΛ)))^q` 随 `a→0` **发散** ⇒ `p(a) > 1` 对充分小 `a` 必然成立 ⇒ **仍判空洞**。这条要单列成一个 species（"log-幂小量"），因为它**逃出** `A7` 的"单参数 `c'` 三档"分类——`A7 §0` 的穷尽性只对**纯** `e^{−c'/g²}` 成立。

**§10.2 记账**：本文件对 (b) 的贡献是**否证其前提**——不存在"全部 `c' ≠ 1/(2b₀)`"，因为不存在 `c'`。故"形状空洞"在这里不是"指数选错"，而是**"已发表侧没有这个形状"**；我们要造的不是更好的 `c'`，是 §5.2 第 1 条的**常数型界**（`c(N)>0` 的来源仍是 `f₁²>0`/(CP)/`G″(0)>0`，见 `CLAUDE.md §5`）。

---

## §6 附：`p` 值表（本文件能给的只有这一列）

`A7 §4` 的表是"假设有一个 `c'`"的表；本文件证明**四篇里没有 `c'` 可填**。为免误用，重列判据本体（`[T]`，与文献无关）：

| `c'` | `p = 2b₀c'` (N=2) | `p` (N=3) | 档 |
|---|---|---|---|
| `1/(2b₀) = 24π²/(11N)` | `p = 1` | `p = 1` | 活（`m ≥ cΛ`） |
| `> 1/(2b₀)` | `> 1` | `> 1` | 空洞 |
| `< 1/(2b₀)` | `< 1` | `< 1` | 荒谬 |

数值核（`b₀ = 11N/(48π²)`，`π²=9.8696044`）：`24π²/22 = 10.766841`（N=2）、`24π²/33 = 7.177894`（N=3）——**与给定值逐位一致**，`b₀` 归一核准通过。

---

## §7 检索记录（可复现 / 未取数清单）

| 目标 | 动作 | 结果 |
|---|---|---|
| Balaban 102 摘要 | `link.springer.com/article/10.1007/BF01229380` | 3036 B 壳页，**未取数** |
| Balaban 116 摘要 | `reo.nii.ac.jp/oja/ART1000344029/en` | **逐字取到**（§1.2） |
| Balaban 116 正文 | `.../ART1000344029.pdf` | 需订阅，`閲覧することができません`，**未取数** |
| Battle–Federbush 摘要 | INIS 记录 + Wisconsin 检索记录 | **逐字取到**（§2） |
| SZZ | `ar5iv …/2204.12737` 全文 | **逐字取到**（§3） |
| Dimock | `ar5iv …/math-ph/0210020` 全文 | **逐字取到**（§4，式 (14)、`mass gap` 仅 1 次） |
| 任一 `e^{−1/(2b₀g²)}` | 全文正则 `e\^\{-`、`exp\(-`、`\exp`、（Balaban 两篇仅题名/摘要层） | **0 命中**；Balaban 正文**未取数**（订阅墙），结论不依赖它（§1 已由题名/摘要排除谱隙） |

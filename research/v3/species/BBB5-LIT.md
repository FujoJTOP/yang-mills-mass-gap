# BBB5-LIT — 4D 纯 SU(N)：`O_h` 十道里有没有任何**一道**的严格质量下界？

> 只判一件事：逐道（或某一特定道）的严格下界，文献里有没有。
> 阅读深度：**【读摘要】** = 读到官方摘要；**【元数据】** = 只核到题名-作者-卷-页-DOI，正文未逐字核。
> 工具：`web_search` OK（7 次，均返回有效来源表）；`web_fetch` OK（HTTP 200 ×3）。

## §0 判词

**[NG]。命中 0 条。** 未查到任何给出"4D 纯 SU(N) 的某一道（`A1± A2± E± T1± T2±` 任一）严格正下界"的文献。三件事必须分开记：

1. **全谱有界 ⇒ 每道平凡继承，这不等于"逐道证明"。** 强耦合区（小 `β`）的格点质量隙定理给的是整个 `H ⊖ Ω` 上的谱底 `> 0`；对子空间取 `inf` 只会 ≥ 全谱值，故十道自动得界。但这是**限制**，非逐道估计；只在 `β < β_0`，无 `a`-一致，不进连续极限。
2. **连续极限侧连全谱都没有。** arXiv:2611.11534 摘要明文：四维空间 cutoff 构造与其上的正谱界 **not established**。全谱未立，逐道无从谈起。
3. **"按道分解"最接近的一篇只到分解，不到界。** arXiv:2509.03513 给 `C = ±1` 两扇区正交分解（`N ≠ 2`；`N = 2` 单扇区），把"小 `β` 下二重一粒子胶球态"写成 **expected** —— 是预期，不是定理。

**结构判**：逐道下界逻辑上是全谱下界的**推论**（`inf` 于子集 ≥ `inf` 于全集）。故"逐道"要有独立内容，必须出现**某道单独可证而全谱不可证**的分离。文献里唯一带算符选择性的论证是 plaquette–plaquette 关联衰减，而 `tr U_p` 只耦合 `A1+`（及 `C=+`）—— 强耦合区给 `A1+` 道的界，但那里全谱界已有，**无增量**。

**下一手**：靶改为"**通道选择性衰减**"——找一个只耦合目标道、且在弱耦合/大 `β` 侧仍有严格衰减的规范不变算符；命中即得该道单独下界。首选 `A1+` 之外的 `T1±`/`E±`（需 3-link/6-link 算符），先核其反射正性与正性度量。

## §1 检索式与命中

| # | 检索式（原样） | 工具 | 命中 |
|---|---|---|---|
| 1 | `rigorous mass gap lattice gauge lower bound` | web_search | 10 源 |
| 2 | `glueball mass lower bound SU(2) rigorous` | web_search | 10 源 |
| 3 | `spectral gap Yang-Mills proof` | web_search | 9 源 |
| 4 | `Osterwalder Seiler mass gap` | web_search | 10 源 |
| 5 | `Balaban exponential decay plaquette correlation large beta SU(2) lattice gauge mass gap rigorous` + `Osterwalder Seiler 1978 … Annals of Physics DOI` | web_search | 10 源 |
| 6 | `Faria da Veiga O'Carroll glueball one-particle states … rigorous small beta` + `Guo Zheng exact ground state … variational upper bound arXiv` | web_search | 10 源 |
| 7 | `rigorous lower bound glueball mass cubic group irrep O_h lattice gauge channel` | web_search | 10 源 |
| 8 | `arxiv.org/abs/2509.03513` / `1406.4177` / `2611.11534` | web_fetch | HTTP 200 ×3 |

命中判定：**0** 条满足"某一道 × 严格下界 × 4D 纯 SU(N)"。全部命中落三类：强耦合全谱隙、变分**上**界、按道分解但无界。

## §2 逐条核

**C1 Osterwalder–Seiler 1978 — 全谱（强耦合），无道标签**
K. Osterwalder, E. Seiler, "Gauge field theories on a lattice", *Ann. Phys.* **110** (1978) 440–471；DOI `10.1016/0003-4916(78)90039-8` **【元数据，Crossref 核过】**。给的是反射正性/OS 重构 + 强耦合区聚类（⇒ 谱底正）；十道只能靠限制继承。严格，但只在 `β < β_0`，衰减率随 `β → 0` 退化。正文未逐字核。

**C2 arXiv:2509.03513 — 扇区分解（`C = ±1`），界未给**
Faria da Veiga, O'Carroll **【读摘要】**（v2, 2026-01-24）。道的粒度是 `C = ±1` 两扇区（`N ≠ 2`；`N = 2` 单扇区），**不是** `O_h` 十道。分解与 Elitzur 型结论是定理；"小 `β` 二重一粒子胶球态"原文为 **expected** ⇒ 非定理，无量值下界。

**C3 arXiv:2611.11534 — 明文承认 4D 正谱界未立**
Mihajlovic **【读摘要】**（414 页）。有限盒内界对"整个 physical vacuum complement"（`inf` 于全集，十道皆被覆盖而无一单独）；有限矩阵判据含 all omitted representations；认证到 `β_W = 0.354`。摘要明文：4D 空间 cutoff 构造与正谱界 **not established**，不主张连续质量。价值 = 本项目问句的**上限声明**。

**C4 arXiv:1406.4177 — 声称全谱正，未采信**
Farinelli **【读摘要】**（v17）声称 `R^{1,3}` 上 Wightman 场、谱正且远离 0。**不作 `[T]` 引用**：`math.GM`、17 版、无评议记录；自述 `g → 0` 时 mass gap → 0，与渐近自由下 transmutation 图像冲突。只记为"存在声称"，不入台账。

**C5 Guo–Zheng 1988 — 方向相反（**上**界）**
S. Guo, W. Zheng, *Phys. Rev. D* **38** (1988) 2591；DOI `10.1103/PhysRevD.38.2591` **【元数据 + 摘要片段】**。`(2+1)` 维 `U(1)`/`SU(2)`/`SU(3)`：精确基态 + 变分**上**界。非 4D、非下界、非逐道。

**C6 Zenodo 自印稿（`20566425`、`23957890` 等）**：称 SU(2) 全谱隙或"One Cube"到 `β_W = 0.348`。**未评议，不引用**；`β` 与 C1/C3 同侧，无独立增量。

## §3 台账（AGENTS §1 四要素）

| 项 | 内容 |
|---|---|
| **① 工具状态** | `web_search` **OK**（7 次，均返结构化来源表，非空响应）；`web_fetch` **OK**（3 次 HTTP 200，取到官方摘要）；`glob`/`read` **OK** |
| **② 实际检索式** | §1 表内 7 条原样列出，未做"已检索"式概括 |
| **③ 命中数** | 满足判据者 **0**；各检索式返回 10/10/9/10/10/10/10 源 |
| **④ 失败标记** | 本次**无 TOOL-FAIL**。唯一边界：C1 正文未逐字核（元数据 + 领域通行表述），条目内已标 |

**标签**：
- `[NG]` 逐道严格下界（4D 纯 SU(N)，`O_h` 十道任一）：命中 0。
- `[NG]` 单一特定道的独立严格下界：命中 0。
- `[T]`（弱形式）强耦合全谱谱底正 ⇒ 十道经限制皆有正下界；严格性来自 C1 一类定理，**不含** `a`-一致与连续极限。
- `[T]`（结构性，本次自证）逐道下界 = 全谱下界之子集限制；其独立内容只能来自"通道选择性衰减"。

**下一手（按序）**：① 核 `T1±`/`E±` 的 3-link/6-link 算符在 Wilson 作用下的反射正性——成立则其关联衰减给该道单独界，不依赖全谱。② 反向核"某道谱底 < 全谱界"的严格反例（即限制式论证不紧）；有则逐道才有独立靶。

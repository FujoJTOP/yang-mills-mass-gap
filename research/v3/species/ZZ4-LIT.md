# ZZ4-LIT — 4D SU(2) Wilson `⟨plaq⟩` 在 β≈2.5 的严格下界：文献有/无

目标量：`plaq(β) = ½⟨tr U_□⟩`，β = 4/g²，Wilson 作用。靶：`plaq(2.5) ≥ 0.625`；实测 0.6526；余 4.4%。

## §0 判词

**[NG]**：在 β=2.5 给 `plaq ≥ 0.625` 的**无条件严格下界**，文献 **0 条**（检索式与命中见 §1，工具状态见 §3）。

已发表材料按 β 窗分三档，没有一档覆盖 2.5：

| 档 | 代表 | 给什么 | β 窗 | 能否到 0.625@2.5 |
|---|---|---|---|---|
| 强耦合严格 | Osterwalder–Seiler 1978 | Griffiths 不等式、Wilson 环单调、面积律 | β 小（本仓 P3 用 β<1/48） | 否（值极小） |
| 弱耦合/面积律 | Cao–Nissim–Sheffield 2025 | Wilson 环**上**界 | β 大 | 否（方向上界） |
| 中耦合 | positivity bootstrap（KZ24、GLYZ25） | plaquette **两侧**界，物理区 0.1% 精度 | 含 λ=3.2↔β=2.5 | **数值上**是；性质上是条件界 |

关键结构事实（本仓可用的唯一严格搬运工具）：Osterwalder–Seiler 的 Griffiths 型单调性 ⇒ `plaq(β)` 随 β 不减。故 β₀ 处的严格下界 `L` 自动成为 β>β₀ 的下界。这条**只搬运数值，不生产数值**；β=2.5 的 L 必须来自 β₀<2.5 的可证区，而该区可证的 L 远低于 0.625。

故 0.625 目前不落在任何已发表严格定理的辖域内；4.4% 余量由 MC + bootstrap 数值界共同支撑，不是定理。

## §1 检索式与命中

| # | 检索式（原文） | 工具 | 命中 | 其中相关 | 失败标记 |
|---|---|---|---|---|---|
| S1 | `rigorous lower bound plaquette SU(2) lattice gauge` | web_search | 10 源 | 1（JHEP12(2025)033，经 springer 代理）；其余 Zenodo 自出版稿，非同行评议，不采信 | 无 |
| S2 | `Osterwalder Seiler plaquette bound` | arXiv API `all:"Osterwalder-Seiler"` | 32 | 1（2505.16585，面积律） | 无 |
| S3 | `mean field lower bound Wilson action` | arXiv API `abs:"mean plaquette" AND abs:"lower bound"` | **0** | 0 | 无（返回 totalResults=0，非工具失败） |
| S4 | `Griffiths inequality lattice gauge plaquette` | arXiv API `abs:"Griffiths" AND abs:"lattice gauge"` | 2 | 2（均为单调性，无 plaquette 值） | 无 |
| S5 | 补充：`abs:"plaquette" AND abs:"bootstrap"` | arXiv API | 3 | 3（KZ24 / GLGY25 / Li-Zhou24） | 无 |
| S6 | 补充：`all:"plaquette" AND all:"rigorous" AND all:"lattice gauge"` | arXiv API | 4 | 1（KZ24） | 无 |
| S7 | 补充：`abs:"plaquette" AND abs:"semidefinite"` | arXiv API | 2 | 0（均非 4D SU(2) 连续群） | 无 |

## §2 逐条核对（给什么界 / β 窗 / 读到哪）

1. **Osterwalder & Seiler, Ann. Phys. 110 (1978) 440**，DOI `10.1016/0003-4916(78)90039-8`。给：格子规范场框架 + Griffiths 型不等式 ⇒ Wilson 环期望随 β 单调、小 β 面积律。β 窗：小 β。**未给任何 β=2.5 数值下界。**【读摘要/元数据（Crossref），正文在付费墙后，未逐字核】
2. **Forsström & Viklund, arXiv:2502.19942**。给：Ising 格子规范理论的随机流展开，另证 Griffiths 第二不等式 ⇒ Wilson 环期望对 β 单调（任意 β）。只给单调性，无值。**【读摘要】**
3. **Cao, Nissim & Sheffield, arXiv:2505.16585**。给：扩大 U(N) 面积律成立区（改进 OS1978），用主回路方程的先验界。方向是 Wilson 环上界。**【读摘要】**
4. **Kazakov & Zheng, arXiv:2404.16925（v4）**。给：4D SU(2) 单迹 bootstrap，`u_P = Tr U_□` 的**上缘与下缘**；耦合 `λ = 2N²/β = 8/β`，故 β=2.5 ↔ λ=3.2；声称 3D/4D 物理耦合区上下缘精度 ~0.1%；回路长度 ≤24。性质：以"无限体积测度存在 + SDP 数值认证"为前提的**条件界**，不是显式常数的定理。**【读正文 §1、§1.1；图 1 为图形，λ=3.2 处具体数值未核】**
5. **Guo, Li, Yang & Zhu, arXiv:2502.14421 = JHEP12(2025)033**，DOI `10.1007/JHEP12(2025)033`。给：SU(3) 为主，含 SU(2)/SU(3) 在 2D–4D 的 plaquette 两侧 SDP 界（文中自述 "rigorous bounds"）；同为条件+SDP 性质。**【读摘要 + §2 正文；4D 图数值未核】**
6. **Li & Zhou, arXiv:2404.17071**。给：Z₂/U(1) 的 bootstrap 两侧界，并指出 bootstrap 界与 Griffiths 不等式的关系。无 SU(2) 4D。**【读摘要】**
7. **Vallée et al., arXiv:2610.06272**。给：SDP 层级对基态能与 Wilson 环的认证区间，示例为 (1+1)D Z₂ 与 (2+1)D。无 4D SU(2)。**【读摘要】**

## §3 台账（四要素 + 断言等级）

四要素：

| 要素 | 记录 |
|---|---|
| 工具状态 | web_search → 有效响应（10 源）；arXiv export API → HTTP 200，XML 完整；Crossref API → HTTP 200；web_fetch(arxiv.org/pdf/*) → `unsupported content type application/pdf`（**TOOL-FAIL**，故 PDF 正文未读，改用 ar5iv/HTML 与摘要）；link.springer.com → 跨域重定向未跟随（**TOOL-FAIL**，JHEP 正文数值未读） |
| 实际检索式 | S1–S7，逐条见 §1，原文照录 |
| 命中数 | S1=10 源 / S2=32 / S3=**0** / S4=2 / S5=3 / S6=4 / S7=2 |
| 失败标记 | 2 处 TOOL-FAIL（PDF 解析、Springer 重定向）；其后改用可读通道，未把失败写成"未查到" |

断言：

- **A1 [NG]** β=2.5 处 `plaq ≥ 0.625` 的无条件严格下界：文献 0 条（S1–S7，工具 OK）。窗口：4D SU(2) 基本表示，β=4/g²。
- **A2 [C]** bootstrap 供给：4D SU(2) 在 λ=3.2（β=2.5）有两侧界，数值与 MC 一致到 ~0.1%，故数值上满足 0.625；成立条件 = 无限体积 Gibbs 测度存在 + SDP 数值认证；本仓未逐点复算。窗口同上。
- **A3 [T]（外部定理）** Griffiths 型单调性：`plaq(β)` 随 β 不减 ⇒ β₀<2.5 的严格下界可直接搬到 2.5。这是唯一可用的严格搬运工，但只搬不造。

下一手：

1. 把靶改写为"β₀∈(0,1/48) 处最大可证 `L(β₀)`" + A3 搬运，查这条链在 2.5 的天花板值（预计远低于 0.625，但需给出数）。
2. 抄录 2404.16925 图 1 在 λ=3.2 的下缘数值（读图或取作者数据），把 A2 从 [C] 抬成可复算的固定数。
3. 查 Balaban 类弱耦合构造是否有显式常数——若有，可补一档"β 大"的定量下界，仍不覆盖 2.5。

## §4 限定

覆盖：4D SU(2)、Wilson 基本表示、超立方格子、`plaq = ½⟨tr U_□⟩`。未覆盖：其他群/表示、有限体积修正、有限格点间距。未核：JHEP12(2025)033 的 4D 图数值、OS1978 正文定理编号、Balaban 类常数。

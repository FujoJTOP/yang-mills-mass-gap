# RG4 · Dimock arXiv:1108.1335 核验（"Balaban 梯子"到底爬到哪）

> 任务：核 `l4a` 的断言「Dimock 1108.1335 **只有 φ⁴₃ 梯子**」。
> 时间戳：2026-11-25（第 N 轮，见 AGENTS.md §5）。工作目录 `D:\Dev\math`。
> 纪律：arXiv ID 已核；只读到**摘要 + 元数据（abs 页）**，正文 PDF 未逐字核。

---

## 0. 工具状态（AGENTS.md §1）

| 项 | 值 |
|---|---|
| 工具 | `web_fetch` |
| URL | `https://arxiv.org/abs/1108.1335` |
| 状态 | **OK**（HTTP 200，返回完整 abs 页） |
| 读到 | **摘要 + 书目元数据**；正文 52 页 PDF **未读** |
| 二次检索 | `web_search` 两式 → **OK**，命中系列 II 及 Dimock 后续 QED₃ 论文（见 §5） |

无 TOOL-FAIL。

---

## 1. 标题与摘要（逐字，abs 页）

**Title:** `The Renormalization Group According to Balaban - I. Small fields`

**Authors:** J. Dimock

**arXiv:1108.1335** (math-ph). v1: 5 Aug 2011；v2: 21 Dec 2012。
Comments: `52 pages. Some corrections, additions, reorganization`
MSC: 81T17 (Primary), 81T08, 81T10, 81T16 (Secondary)
Journal DOI: `10.1142/S0129055X13300100`（Rev. Math. Phys. 25 (2013) 1330010）

> **Abstract:** This is an expository account of Balaban's approach to the renormalization group. The method is illustrated with a treatment of the ultraviolet problem for the scalar phi^4 model on a toroidal lattice in dimension d=3. This yields another proof of the stability bound. In this first paper we analyze the small field contribution to the partition function.

---

## 2. 它证的是什么（`φ⁴₃` 的什么性质？）

- **对象**：**环面格点**（toroidal lattice）上的**标量 `φ⁴`**，**d = 3**，`d=3` 即 `φ⁴₃`。
- **问题**：**紫外问题**（UV problem）—— 去掉格点截断后配分函数是否仍有意义。
- **结论类型**：**stability bound**（紫外稳定性界，即 `Z_Λ` 的一致下界/上界，形如 `e^{-C|Λ|}`）。
- **性质归属**：这是**配分函数的 UV 稳定性**，**不是**质量间隙、**不是**谱隙、**不是**相变、**不是**连续极限的存在性（后者需 `d=4` 或超可重整模型）。
- **本征贡献**：`another proof` —— 复证（Balaban 原本的证明路径的**重述**版），不是新定理。
- 本册（I）**只做小场**部分：`In this first paper we analyze the small field contribution to the partition function.` 大场在 **II**（`arXiv:1212.5562`，J. Math. Phys. 54, 092301）。

⇒ `l4a` 断言**成立**：**ID 核过，内容确实只有 `φ⁴₃` 的 UV 稳定性复证**。`[T]`（读到摘要）/ 正文未核。

---

## 3. 方法：是不是 Balaban 式多尺度？

**是。** 摘要自称 "Balaban's approach to the renormalization group"，即 **BFM（block field renormalization / 块场重整化）多尺度展开**，其骨架为：

1. **逐尺度块平均**（block averaging，`a → 2a`，`Λ_j = Λ_0 2^{-j}`）；
2. **小场 / 大场分解**（small-field / large-field decomposition）——
   小场：有效作用量在**解析、有界、polymer 活动**函数空间内闭合并收缩（本册 I）；
   大场：用 polymer 气体 / 熵计数压住（册 II）；
3. **有效作用量的逐阶展开**（`T_{j+1}` 迭代），误差项按 `2^{-j}`（或 `a^δ`）幂次吸收；
4. **稳定界**：`Z ∈ [e^{-C|Λ|}, e^{+C|Λ|}]`，`C` 与截断无关。

关键点：**它是一套"框架/语言"**（解析聚合物活动 + 显式常数跟踪），不是一条单一定理。

---

## 4. 4D gauge 能借什么：方法，还是只有常数？

**方法（骨架）可借，具体估计不可借。** 分三层：

| 层 | 4D gauge 可否借 | 说明 |
|---|---|---|
| **L1 多尺度骨架** | **可借（已借，Balaban 本人做过）** | Balaban 1980s 系列对 **4D 纯格点 `SU(2)` 规范场**做过 UV 稳定性：块平均 + 小/大场 + polymer。Dimock 1108.1335 是**把这套语言写成可读版**，用的样本是 `φ⁴₃`（更简单、无规范）。 |
| **L2 估计与常数** | **部分可借** | 聚合物活动的 `C^n/n!` 型界、`2^{-j}` 幂次计数、显式常数链是**通用机器**，换成规范场后形式不变。 |
| **L3 规范特有的那一步** | **Dimock 本书不含** | 见 §6 的"不可移植步"。 |

Dimock 自己后续把同一语言推进到 **QED₃**（`arXiv:1502.02946`，Ann. Henri Poincaré；"Ultraviolet stability for QED in d=3"）—— 这证明：**这套语言确实能从标量搬到"带规范"的模型**，但只搬到 **d=3 的 QED**（超可重整），**没搬到 d=4 的非阿贝尔**。

---

## 5. 系列与邻近文献（`web_search` 命中，元数据级）

| 文献 | 内容 | 状态 |
|---|---|---|
| 1108.1335 (I) | `φ⁴₃` 小场 | 摘要已读 |
| **1212.5562** (II) | `φ⁴₃` **大场** | 摘要已读（摘要同句式）|
| 1502.02946 | **标量 QED₃** 非微扰重整化（用 Balaban 表述） | 摘要级 |
| Dimock, *UV stability for QED in d=3* | QED₃ UV 稳定界（有限体积） | 摘要级 |
| Balaban 1980s 系列 | **4D 纯格点 SU(2)** UV 稳定性（规范固定） | 二手（本次未直接核） |

> 注：Balaban 系列"建立了规范固定下纯格点规范配分函数的 UV 稳定性"，**不含** OS 重构、**不含**质量间隙、**不含**规范不变连续极限。此判据来自检索到的第三方综述性 PDF（zenodo 文档），**未逐字核 Balaban 原文**，标 `[C]`。

---

## 6. 判：这条"梯子"能爬到 4D gauge 吗？

**答：能爬到"4D 规范场的 UV 稳定性"，爬不到"4D 规范场的质量间隙"。**

### (a) 可移植的那一步

- **块平均 + 小/大场分解 + 聚合物活动闭合**（`T_j` 迭代、`C^n/n!` 界、`2^{-j}` 幂次跟踪）。
- **证据**：Balaban 已对 4D 纯 `SU(2)` 做过；Dimock 语言在 QED₃ 上复用过。⇒ 骨架**不是** `φ⁴₃` 专用。

### (b) 不可移植 / `φ⁴₃` 专用的那几步

1. **维数**：`φ⁴₃` 是**超可重整**（`[φ]=1/2`，`φ⁴` 边缘但 `d<4` 使 `φ⁴` 相关/边缘混合，UV 发散有限阶）。4D gauge 是**严格可重整**，`1/g²` 边缘、`F²` 边缘 —— 逐阶抵消结构不同，**幂次计数不能照抄**。
2. **规范**：`φ⁴₃` 无规范冗余。4D gauge 必须**规范固定**（Balaban 用轴向/Landau 规范 + 换规范），否则块平均破坏规范不变性；换规范产生 **Faddeev–Popov / Gribov** 类问题，`φ⁴₃` 完全没有。
3. **规范不变性 vs 稳定性**：即便拿到 `Z` 的 UV 稳定界，那只是**规范固定**下的配分函数界；提升到**规范不变 OS 理论**是**额外一步**，Dimock 本书不含，Balaban 系列亦未闭合。
4. **目标不同**：`φ⁴₃` 要的是"`Z` 有界"；质量间隙要的是"**转移矩阵谱隙 > 0 且一致**"。UV 稳定性对谱隙**零信息**。

### (c) TOOL-FAIL

**无。** `web_fetch` HTTP 200，`web_search` 有效。

---

## 7. 结论（供总册引用）

- **`l4a` 断言成立**：1108.1335 = `φ⁴₃`（d=3 环面标量 `φ⁴`）**UV 稳定界**的**复证**，册 I 只做小场。ID 已核，读到摘要。 `[T]`
- **方法 = Balaban 多尺度 BFM 的可读重述**，非新方法。 `[T]`
- **4D gauge 可借**：L1 骨架 + L2 聚合物计数机器（Balaban 已在 4D `SU(2)` 用过；Dimock 已把语言搬到 QED₃）。 `[C]`
- **不可借**：维数（超可重整 ↔ 可重整）、规范固定/FP、UV 稳定 → OS 规范不变理论、UV 稳定 → 谱隙。 `[T]`（结构性判断）
- **爬不上去的那一步**：**从"规范固定的 UV 稳定界"到"规范不变的质量间隙"**。这正是本项目宪章 §5 说的"那台机器尚不存在"。

---

## 8. 与本项目位置的关系

- 本文件**不改变** `docs/position-current.html` 的五环/五洞状态：它只是**给 RG 路线（R 环）加一条文献边界注**——RG 骨架已有 4D 规范先例，但先例止于 UV 稳定性，未触及谱隙。
- 死路台账新增一条（编号待总册分配）：
  **「Balaban/Dimock 梯子 → 质量间隙」= 断在"稳定性 ≠ 谱隙"**（`[T]`，结构性的）。

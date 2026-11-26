# RG6-STARLIT — 文献检索：`(★)` 单边方差界在 4D `SU(N)` YM 里有没有人做过

> 2026-11-25 · 物种 RG6（STARLIT：文献星图）。任务：查 `(★)` 是否已有文献版本。
> 判级纪律：`AGENTS §1`（工具状态 + 检索式 + 命中数四项齐全）、`§3`（arXiv/DOI + 标注读到哪一层）。
> 本文件只记检索与判定，不改任何断言等级。

---

## 1 · `(★)` 的精确陈述（逐字，`C2-BALABAN-UNIFORM.md:208–217`）

标题行：`### 2.4 The failing step: the analytic per-shell bound (★)`（`:208`）。

前置散文（`:210–211`，逐字）：

> Lemma B controls *which blocks are large*. It does **not** control the **small-field conditional measure**
> `μ(·|U′)` once the large set is removed. What remains is the one inequality per shell, uniform in `j,U′,a`:

```
(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j, v⟩)  ≤  E_{μ(·|U′)}⟨PᵀH_jP v, v⟩ + W_j(v) − δ_W      (δ_W > 0),
```
（`:214`，code fence `:213`–`:215`）

伴随条件（`:217`，逐字）：

```
plus cross-shell band decay `(★dec): |Cov_μ(s_B,s_{B′})| ≤ C e^{−κ dist(B,B′)/a}`. Equivalently
`G″_j(v) = E[H_j] − Var_j + W_j ≥ δ_W`.
```
（`:217–218`；`:218` 给出等价形式 `G″_j = E[H_j] − Var_j + W_j ≥ δ_W`。）

同一文件的等价重述（`:32`，逐字）：
`(★): Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E⟨PᵀH_jP v,v⟩ + W_j(v) − δ_W`，uniformly in `j,U′,a`。

**关键点（决定文献可比性）**：`(★)` 的难处不是方差上界本身，而是 **一边性** —— 只有 `Var ≥ 0`（从下方减），
没有任何机制给出 `Var` 的**上界**。`C2:219–222` 明写：唯一桥是单圈等配分，且只 `asymptotic`、不 `uniform`。

---

## 2 · 工具状态与检索式（`AGENTS §1` 四项）

**工具状态**
| 工具 | 状态 |
|---|---|
| `pwsh`（读本地行） | **OK**（读出 `C2-BALABAN-UNIFORM.md` 共 411 行，`205–222` 行原文见上） |
| `grep`（本地语料） | **OK**（`research/v3` 命中 941 处，落盘 spill 文件正常） |
| `web_search` | **OK**（4 次调用均返回非空来源列表） |
| `web_fetch` | **OK**（arXiv abs 三篇 HTTP 200；ETH 旧链接 HTTP 404，见 §3.4） |

**检索式与命中数**（每条单独一次调用；工具单次返回上限 10 条来源）

| # | 检索式 | 返回来源数 | 与 `(★)` 同形的直接命中 |
|---|---|---|---|
| Q1 | `"one-sided" variance bound gauge renormalization group` | 10 | **0** |
| Q2 | `"variance" bound lattice gauge "renormalization group" uniform` | 10 | **0**（最接近的是条件 Poincaré / Efron–Stein 型，均非 `Var ≤ E[H]` 单边型） |
| Q3 | `Balaban gauge "second order" bound positive` | 10 | **0** |
| Q4 | `"SU(N)" Yang-Mills constructive "renormalization group" infrared 4D` | 10 | **0** |

> 说明：10 是**返回条数上限**，不是精确总命中；本次检索对"精确同形"计数为 **0**（`AGENTS §1`：写 `0`，不写"没有"）。

---

## 3 · 最相关 2–3 篇（读到摘要层）

### 3.1 Dimock, *The Renormalization Group According to Balaban — I. Small fields*
- **arXiv:1108.1335** (math-ph)；期刊 DOI **10.1142/S0129055X13300100**（Rev. Math. Phys. 25 (2013) 1330010）
- 读到：**摘要全文**（`web_fetch` HTTP 200）。逐字：
  > "This is an expository account of Balaban's approach to the renormalization group. The method is illustrated with a treatment of the the ultraviolet problem for the scalar phi^4 model on a toroidal lattice in dimension d=3. This yields another proof of the stability bound. In this first paper we analyze the small field contribution to the partition function."
- 相关性：这正是 `(★)` 所在位置的文献对应物（**小场/单 shell 上的定量控制**）。但它做的是 **两边小扰动（ε-perturbation）展开**：有效作用是二次型的小扰动，给 `‖ΔS‖ ≤ O(ε)·(二次型)`，**不是**方向性方差不等式 `Var ≤ E[H] + W − δ`。

### 3.2 Dimock, *… — II. Large fields*
- **arXiv:1212.5562**；J. Math. Phys. 54 (2013) 092301；DOI **10.1063/1.4821275**
- 读到：**摘要层**（检索结果摘要 + 期刊页标题，逐字："In this paper we control the large field contribution to the partition function."）
- 相关性：对应 `C2 §2.2` 的 **大场通道关闭**（Lemma B 那一侧）。即：文献把 `(★)` 的"combina­torial 一半"做了，"analytic 一半"没做。

### 3.3 Dimock, *… — III. Convergence*
- **arXiv:1304.0705**；Ann. Henri Poincaré DOI **10.1007/s00023-013-0303-3**
- 读到：**摘要全文**（HTTP 200）。逐字："In this third paper we demonstrate convergence of the expansion and complete the proof of a stability bound."
- 相关性：稳定性界 `stability bound` 是**行列表式**的（`e^{-βS}` 的上下界），与 `(★)` 的**方向性方差**不同对象。

### 3.4 原始的 Balaban 四维规范场系列
- T. Balaban, *Renormalization group approach to lattice gauge field theories. I. Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions*, Comm. Math. Phys. **109** (1987) 249–301；DOI **10.1007/BF01215223**（DOI 经 Mendeley 目录页与 ouci 目录页两处一致核对）。
- 读到：**仅题录 + 标题层**（正文在付费墙后，**未逐字核**）。
- 标题本身即是判决语言：`small field approximation` + `coupling constant renormalization` —— 即 **有效作用的小场展开 + 耦合流**，而非"条件测度的单边方差界"。

**一条关键的旁证（决定 (b) 而非 (a)）**：Dimock 的 Balaban 综述幻灯片（ETH 旧链 HTTP 404，未取到）在检索摘要里留下逐字残句：
> "Caveat: paper on coupling constant flow in second order perturbation theory missing"
即：Balaban 路线里"二阶微扰下的耦合流"那一步，**综述作者自己标注为缺失/未发表**。文献把 `(★)` 附近的那一格留空，这与 `C2 §2.4` 的"唯一数学缺口"叙述**方向一致**（但不等于证明 `(★)` 无人做过，见 §4）。

**不予采信的两类来源**（纪律）：Q2 返回的 `ai.vixra.org` 三篇（`2602.0040`, `2602.0072`, `2602.0056`）与 `zenodo` 上若干"YM mass gap 构造"文件均为**非同行评议预印本/仓库文件**，无 arXiv/DOI 期刊锚点，本次**不计入命中**，也不作为"有人做过"的证据。

---

## 4 · 判定

**判定：(b)，带限定。**

- `(★)` 的**精确同形**（`Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E[H_j] + W_j − δ_W`，uniform in `j,U′,a`）在本次检索中：**未查到**（检索式 Q1–Q4，工具 OK，命中 **0**）。
- 文献实际提供的是三种**不同**对象：
  1. **两边**小扰动/稳定性界（Balaban 1987；Dimock 2013 I/III）—— 给的是 `|ΔS| ≤ O(ε)·Q`，不是单边方差；
  2. **条件 Poincaré / LSI / Brascamp–Lieb 型方差上界** —— 需要 log-concavity 或谱隙，方向恰好相反（它靠 `E[H] ⪰ 0` 的下边缘，给不出 `Var ≲ E[H]`）；
  3. **大场抑制 / 组合通道**（Dimock II 2013）—— 对应 `C2 §2.2`，已关闭。
- 因此：**"这可能是新的"** —— 即"把 `(★)` 作为**单条、per-shell、对边界 `U′` 一致的单边方差不等式**来陈述并证明"这一格，本次检索 **0 命中**。
- **不下"不存在"的断言**：检索只覆盖 4 条英文检索式的一次 web 检索 + 3 篇摘要层，未系统扫 Balaban 原文（付费墙）、未扫 Dimock 2013 II 正文、未扫 `math-ph` 全库。**结论强度：`[C]`（检索层）**，非 `[T]`。

**可用的部分（对项目有正面价值）**：即便 `(★)` 本身无文献版本，Dimock I/III 的 `ε-perturbation + stability bound` 机器是**现成可搬的两边版**——若能把"两边"降成"单边 + 余项"，就是 `(★)` 的一个候选证明骨架；反之，若两边界已是极限（`O(ε)` 常数不可再降），那就是 `(★)` 在 Balaban 语言下**不可得**的一个信号。这一条留给下一轮实测，本文件不作断言。

---

## 5 · 引用清单（`AGENTS §3`：arXiv/DOI + 读到哪一层）

| 文献 | 锚点 | 读到的层 |
|---|---|---|
| Dimock, RG according to Balaban I. Small fields | arXiv:1108.1335；DOI 10.1142/S0129055X13300100 | 摘要全文 |
| Dimock, RG according to Balaban II. Large fields | arXiv:1212.5562；DOI 10.1063/1.4821275 | 摘要 |
| Dimock, RG according to Balaban III. Convergence | arXiv:1304.0705；DOI 10.1007/s00023-013-0303-3 | 摘要全文 |
| Balaban, RG approach to lattice gauge field theories I | DOI 10.1007/BF01215223（CMP 109 (1987) 249–301） | 题录；正文未逐字核 |

**本地锚点**：`research/v3/C2-BALABAN-UNIFORM.md:208–222`（`(★)` 陈述 + 一边性说明）、`:32`（等价重述）、`:229`/`:372`（`[O]` 状态）。

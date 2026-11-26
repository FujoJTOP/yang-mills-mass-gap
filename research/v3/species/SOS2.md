# SOS2 — 算子 SOS（Helton–McCullough / Putinar / Lasserre）对**无界** `H` 与**无穷维** `L²(SU(2)^24)` 覆盖到哪一步

**日期**: 2026-11-25（第 N 轮）· **裁决**: **(b) 代数 SOS 用不上** + **(c) 只剩"截断 + 尾界"，而尾界就是 `d02` 的那堵墙** · **禁 git commit**（遵守）
**只写**: 本文件
**读入**: `D02-DIRECT.md`（`T` 载体 ∞ 维、`T` 紧、PF 不出数）、`Q12-SOS.md`（非负系数 SOS 被 `a=0.75` 反例杀）、`SDP-CERT.md §2.3`（moment-SOS ≠ 矩矩阵）、`B05-DUAL.md:256`、`B06-VAR.md:60`（P3 需**显式 `H`**）
**任务**: `P3 SOS` = 要证 `H ≥ c`（算子不等式，形式核上）。本文只判"这条定理库能不能供它"。

---

## 0. 判死线（先读）

| 问 | 答 | 标签 |
|---|---|---|
| 有"非交换 Positivstellensatz"吗 | **有**，Helton 2002 / Helton–McCullough 2004 / Putinar 1993 / Lasserre 2001，四条逐字见 §1 | `[T]` |
| 它们覆盖**无界算子**吗 | **不覆盖**。假设全部落在"有界自伴算子元组 / 紧半代数集 / Archimedean"上；Archimedean \(\Rightarrow\) 每个生成元**有界**（§2.1 两行证明） | `[T]` |
| 有**一个**真的无界情形定理吗 | **有且仅有一条**：Schmüdgen 2005 Weyl 代数（§1.5）。但它要 CCR + "在 Schrödinger 表示里已是严格正算子" ⇒ 对我们**循环**，且结论是 \(b_i c b_i^*\) 的**加权**平方和，不是 \(c\) 本身 | `[T]` |
| `L²(SU(2)^24)` 无穷维 ⇒ 无穷多平方，有收敛定理吗 | **代数意义的收敛定理不存在**；存在的只有 (i) 形式级恒等式（§3.1，恒真但**预设**结论）、(ii) 截断 + 尾界（§3.2，`[T]`，但要求尾在**算子序**下 \(\to0\)） | `[T]`/`[NG]` |
| 对 YM 可行吗 | **不可行**（(b)）；唯一存活形式 = (c)，其缺口 \(\varepsilon_N\) 正是 `d02` 的墙 | `[NG]` |

**一句话**：SOS 的两半不对称 —— **"SOS ⇒ 正"** 恒真但要求先有恒等式（且要求公共核；这本身就是证明义务）；**"正 ⇒ SOS"**（唯一能**产出**证书的方向）要求生成元有界 + 紧/Archimedean + 有限维。YM 侧 `H` 的无界部分（Laplacian/动量）不在任何有界生成元的 polynom 代数里，且 `H` 本身非显式 ⇒ 连命题都摆不出来。

---

## 1. 逐字陈述（每条：假设 / 结论 / 是否覆盖无界）

> 标注口径：**[摘要级]** = 只核到 arXiv/出版社摘要，未逐字核正文（付费墙）。**§3 的 DOI 均核到出版社/Crossref 页面**。

### 1.1 Putinar 1993 —— Archimedean（交换、有限维）
Putinar, *Positive polynomials on compact semi-algebraic sets*, **Indiana Univ. Math. J. 42(3) (1993) 969–984**, DOI [`10.1512/iumj.1993.42.42045`](https://doi.org/10.1512/iumj.1993.42.42045)。
- **假设**：\(S=\{x: g_1\ge0,\dots,g_m\ge0\}\subset\mathbb R^n\) **紧**，且二次模 \(M(g)=\sigma_0+\sum\sigma_ig_i\) **Archimedean**：\(\exists N,\ N-\sum_j x_j^2\in M(g)\)。
- **结论**：\(f>0\) on \(S\ \Rightarrow\ f\in M(g)\)（\(\sigma_i\) 为 SOS）。
- 覆盖无界？**否**。假设里两处硬限：\(n<\infty\)、\(S\) 紧（\(\Leftrightarrow\) 坐标有界）。`[T]` [摘要级]

### 1.2 Lasserre 2001 —— moment-SOS 层级（有限维）
Lasserre, *Global optimization with polynomials and the problem of moments*, **SIAM J. Optim. 11(3) (2001) 796–817**, DOI [`10.1137/S1052623400366802`](https://doi.org/10.1137/S1052623400366802)。
- **假设**：\(n<\infty\)、\(S\) 紧（或 Archimedean）、**次数截断** \(d\)。
- **结论**：半定松弛值 \(\uparrow \min_S f\)（单调收敛，有限步精确当 \(S\) 有限/离散）。
- 覆盖无界/无限维？**否**。收敛性完全建立在"矩矩阵 \(M_d\) 有限阶 + \(S\) 紧"上。`[T]` [摘要级]

### 1.3 Helton 2002 —— 非交换：矩阵正 \(\Rightarrow\) SOS
Helton, *"Positive" noncommutative polynomials are sums of squares*, **Ann. of Math. 156(2) (2002) 675–694**, DOI [`10.2307/3597203`](https://doi.org/10.2307/3597203)。
- **假设**：\(p\) 非交换多项式（\(d\) 个自伴生成元 \(X\)），**矩阵正**：对**一切**尺寸的一切自伴矩阵元组 \(A\)，\(p(A)\succeq0\)。
- **结论**：\(p=\sum_j r_j^*r_j\)（非交换 SOS）。
- 覆盖无界？**否**：假设是"对一切**有界**（矩阵）元组"。**注意方向**：这是"正 \(\Rightarrow\) SOS"，不是"逐元组正 \(\Rightarrow\) SOS"。`[T]` [摘要级]

### 1.4 Helton–McCullough 2004 —— 非交换 Positivstellensatz（有界半代数集）
Helton & McCullough, *A Positivstellensatz for non-commutative polynomials*, **Trans. AMS 356(9) (2004) 3721–3737**, DOI [`10.1090/S0002-9947-04-03433-6`](https://doi.org/10.1090/S0002-9947-04-03433-6)。
- **摘要逐字**（出版社摘要）："A non-commutative polynomial which is positive on a **bounded** semi-algebraic set of operators has a **weighted** sum of squares representation."
- 覆盖无界？**否**，且"**bounded**"写在摘要里；结论是**加权** SOS（\(s_0+\sum s_i g_i s_i^*\) 型），不是裸 SOS。`[T]` [摘要级]

### 1.5 Schmüdgen 2005 —— **唯一**的无界生成元定理（Weyl 代数）
Schmüdgen, *A strict Positivstellensatz for the Weyl algebra*, **Math. Ann. 331 (2005) 779–794**, arXiv [`math/0403076`](https://arxiv.org/abs/math/0403076)。
- **摘要逐字**："Let \(c\) be an element of the Weyl algebra \(W(d)\) which is given by a **strictly positive operator in the Schrödinger representation**. It is shown that, **under some conditions**, there exist elements \(b_1,\dots,b_d\) in \(W(d)\) such that \(b_1cb_1^*+\cdots+b_dcb_d^*\) is a finite sum of squares."
- 覆盖无界？**是**（\(W(d)\) 的生成元即无界 \(p,q\)，CCR）。**但对 YM 三条全不适用**：
  1. 要 **Weyl 代数**（CCR，\([p,q]=i\)）。格点 YM 的 `T` 载体上，位置是 `SU(2)` 上的**有界**函数（紧群），动量是 Laplacian；生成的不是 \(W(d)\)。
  2. 假设是"**在 Schrödinger 表示里已是严格正算子**" —— 这正是我们要证的 `H ≥ c`，**循环**。
  3. 结论是 \(b_icb_i^*\) 的**加权**平方和，不是 \(c=\sum r^*r\)；`b_i` 依赖表示无关的代数元素，不给出可算证书。
- `[T]` [摘要级]

### 1.6（旁证，2025）算子值非交换多项式的单平方分解
Jindal–Klep–McCullough, *Positive operator-valued noncommutative polynomials are squares*, arXiv [`2511.06487`](https://arxiv.org/abs/2511.06487)（IEOT 98(1):7, 2026）。
- **摘要逐字**："every positive operator-valued noncommutative polynomial \(p\) admits a **single-square** factorization \(p=r^*r\)"。
- 覆盖无界？**否**（算子值 = 系数取 \(\mathcal B(\mathcal H)\)，生成元仍有界）。`[T]` [摘要级] —— 但它证明**"单平方"在非交换侧已是极限形式**，与本项目 §3.1 的形式恒等式同型。

### 1.7 域要求（无界侧的正确定理，非 Positivstellensatz）
Nelson, *Analytic vectors*, **Ann. of Math. 70(3) (1959) 572–615**, DOI [`10.2307/1970331`](https://doi.org/10.2307/1970331)。
- **内容**（标准表述）：若 \(D\) 是 \(A_1,\dots,A_d\)（自伴）的**公共解析向量核**，则每个非交换多项式 \(p(A_1,\dots,A_d)\) 在 \(D\) 上本质自伴。
- **反面**：没有公共核时 \(p(A)\) 的**定义域**可能塌缩，"\(p(A)\succeq0\)"根本不是命题。⇒ **无界侧 SOS 的第一道门槛不是代数，而是公共核**。`[T]` [摘要级]

---

## 2. 关键：为什么"无界"把 Positivstellensatz 整条假设打掉

### 2.1 Archimedean \(\Rightarrow\) 生成元有界（两行，自证）
\(N-\sum_j x_j^2\in M\) 且 \(M\) 里每个元素在任何 \(*\)-表示下都 \(\succeq0\)。取任一自伴元组 \(A\)：
\[
0\preceq N-\textstyle\sum_j A_j^2\ \Longrightarrow\ \sum_j A_j^2\preceq N\ \Longrightarrow\ \|A_j\|\le\sqrt N .
\]
⇒ **Archimedean 条件按定义就是"生成元一致有界"**。无界生成元上 Archimedean 恒假，Putinar/Lasserre 的整条链条（紧性 \(\to\) 矩矩阵有限阶截断 \(\to\) 收敛）无处落脚。`[T]`（自证）

### 2.2 两个方向不对称（这是本文的主结论）
| 方向 | 是否真 | 代价 |
|---|---|---|
| \(\sum_j A_j^*A_j\succeq0\ \Rightarrow\ H\ge c\) | **恒真** | 需 (i) 显式恒等式 \(H-c=\sum A_j^*A_j\)；(ii) \(A_j\) 闭、**公共核**（Nelson §1.7）；(iii) 该核是形式核 |
| \(H\ge c\ \Rightarrow\ \exists\) 有限 SOS（**产出**证书的方向） | **只在有界 + 紧/Archimedean + 有限维**下真（§1.1–1.4） | YM 侧两条全缺 |

⇒ **SOS 作为"发现证书的机器"在 YM 上不存在**；作为"验证模板"它退化成"你自己把 `H−c` 配成平方"，即回到 `H1`（`Q12-SOS.md §5` 已判：SOS ⟺ `H1`）。

### 2.3 格点 `T` 的结构诊断（本题专用）
`T` 的载体是 \(L^2(SU(2)^{24})\)，坐标 \(U_l\) 取值于**紧**群 ⇒ 位置方向**有界**，Archimedean 在 \(U\) 方向**免费成立**。但：
- 免费的那一半**没用**：`T` 的无界性全在动量/Laplacian 方向（`H = −Δ_{SU(2)^{24}} + V(U)` 型），而 \(-\Delta\) **不是**有界生成元的多项式；
- 而 `H` 在转移矩阵语言里**根本不显式**（`B05-DUAL.md:256`、`B06-VAR.md:60`：4D 转移算子非显式 ⇒ `[C]` 不可及）。⇒ 连"把 `H−c` 写成什么"都写不出。
- `[NG]`

---

## 3. 无穷维：有收敛定理吗

### 3.1 形式级恒等式（恒真、但**预设**结论）
设 \(q\) 是稠定、闭、下有界的形式，\(q\ge c\)，\(A\) 为其 Friedrichs 自伴算子（Schmüdgen, *Unbounded Self-adjoint Operators on Hilbert Space*, **GTM 265**, DOI [`10.1007/978-94-007-4753-1`](https://doi.org/10.1007/978-94-007-4753-1)，形式/Friedrichs/KLMN 各章）。则
\[
q(u)-c\|u\|^2=\big\|(A-c)^{1/2}u\big\|^2,\qquad u\in\operatorname{Dom}(q).
\]
- 这是**对无界算子成立的"单平方 SOS"**，无穷维也不需要收敛定理。`[T]`
- **但它是循环的**：右端的 \((A-c)^{1/2}\) 只在已知 \(A\ge c\) 后存在。⇒ 它把 `H ≥ c` 变成"`H ≥ c`"，零信息。（与 §1.6 的 single-square 同型：非交换侧"一个平方"也是极限形式。）
- 唯一**非**循环的无界 SOS 是形式级 \(H=\sum_j A_j^*A_j\) 且 \(A_j\) 显式（如 \(-\Delta=\sum_j\partial_j^*\partial_j\)，公共核 \(C_c^\infty\)）。`[T]` —— 但自由场过不了酸测（`f₁²=0`），对 YM 无内容。

### 3.2 无限的"收敛"只能是截断 + 尾界（`[T]` 引理，缺口 \(\varepsilon_N\)）
**引理（截断–尾界）**：设 \(H,H_N\) 自伴，且存在 \(\varepsilon_N\ge0\) 使**算子序**
\[
-\varepsilon_N\,\mathbb 1\ \le\ H-H_N\ \le\ \varepsilon_N\,\mathbb 1 .
\]
若 \(H_N\ge c_N\)（有限维截断上的 SOS/矩阵证书给出），则 \(H\ge c_N-\varepsilon_N\)。
- 证明：\(\langle u,(H-c_N)u\rangle=\langle u,(H_N-c_N)u\rangle+\langle u,(H-H_N)u\rangle\ge0-\varepsilon_N\|u\|^2\)。`[T]`（自证）
- **缺口精确为 \(\varepsilon_N\)**：需要 \(\varepsilon_N\to0\) 且 \(c_N\to c\) 同速率。对 `T`：截断是 Peter–Weyl/`j`-截断（`D02-DIRECT.md §3.1` 路 I），尾项**没有**算子序界，只有核的正性与紧性（`D02 §2`）；而 `T` 的载体**本质无穷维**（`D02-DIRECT.md §1`，`L²(SU(2)^24)` 可分无穷维，`3L³=24` 只是**矩阵元的积分维数**）。
- ⇒ **`d02` 的那堵墙 = 引理里的 \(\varepsilon_N\)**。`[NG]`
- 附注（为何不能用形式收敛糊过去）：强预解收敛只给**谱的收敛**，不给**下界**；下界要**算子序/范数**级别的尾界。`[T]`

---

## 4. 判：SOS 对 YM 可行吗

**(b) 成立 —— 代数 SOS 用不上**，四条理由各带证：
1. `H` 的无界部分不在有界生成元多项式代数里；Archimedean 恒假（§2.1）。
2. Positivstellensatz 的"正 \(\Rightarrow\) SOS"方向要求有限维紧集或有界算子元组（§1.1–1.4）。
3. 唯一无界定理（§1.5）假设循环、要 CCR、结论加权（§1.5 三条）。
4. `H` 非显式（`B05/B06`），连命题都摆不出（§2.3）。

**(c) 成立 —— 只剩"截断 + 尾界"**，即 §3.2 的引理；其 \(\varepsilon_N\) 就是 `d02` 的墙，也是 `t→0/t→∞` 桥（`LOOP-RUNBOOK §70.3`）缺的那一段。

**(a) 部分成立但不救**：无界情形**有**定理（§1.5），但**不覆盖**本问题的结构（CCR／非循环假设／加权结论）。

**净结论**：`P3 SOS` 的"代数证书"分支 `[NG]`；可存活的只有**形式级**版本：写显式 `H = Σ_j A_j^*A_j + R`（Bochner/Lichnerowicz 型，YM 里 `R` = 曲率/交换子项 `ad F`），于是 SOS 那一步化为"**逐场位形** `R ≥ c`"。而 `R` 正是 Savvidy 负模所在（本项目实测 `Hess−D²+2adF_B` 7.1% 负，`CLAUDE.md §5`）⇒ 该分支持续负，和 `Q12-SOS` 的 `Var≤−0.0367` 是同一堵墙的两个坐标。

**路径（逐字）**：`P3` 若还要走，必须换**证的方向**：不做"正 \(\Rightarrow\) SOS"，而做"**显式 `H` + 公共核 + 形式级 Bochner 恒等式**"。前提是 `H` 显式 —— 这正是 `B05-DUAL.md:256`/`B06-VAR.md:60` 判为 `[C]` 不可及的那一步。**SOS 不提供它。**

---

## 5. 检索核验（AGENTS §1）

| 工具 | 状态 | 检索式 | 命中 |
|---|---|---|---|
| `grep` | **OK** | `d02\|SOS\|Helton\|Lasserre\|Putinar` @ `research\v3\**\*.md` | **54** |
| `glob` | **OK** | `research/v3/species/SOS*.md` | **0**（本文件为新建） |
| `web_search` ×3 批 | **OK** | `Helton "Positive" noncommutative polynomials Annals DOI`；`Helton McCullough Positivstellensatz Trans AMS DOI`；`Putinar 1993 DOI 10.1512`；`Lasserre 2001 DOI 10.1137/S1052623400366802`；`Schmüdgen Unbounded Self-adjoint GTM 265 DOI`；`Marshall Positive Polynomials and Sums of Squares 10.1090/surv/146`；`"Positivstellensatz" "unbounded operators" arXiv`；`sum of squares unbounded self-adjoint common core theorem`；`Nelson analytic vectors DOI` | 各批 10 条摘要 |
| `web_fetch` | **OK** | `arxiv.org/abs/0709.3170`、`/abs/math/0403076`、`/abs/2511.06487` | 200，摘要逐字取到 |
| `web_fetch` | **403** | `zbmath.org/pdf/02156244.pdf`（Cloudflare） | 改用 arXiv 摘要成功 |

**未查到（检索式 `"Positivstellensatz" "unbounded operators"`；工具 OK；命中 0 条直接命中的"紧群/带 Laplacian 的算子代数"Positivstellensatz）** —— 最接近的只有 §1.5（Weyl 代数）。**这是 0 命中，不是 TOOL-FAIL。**
**未逐字核**：§1.1–1.6 全部**仅摘要级**（正文付费墙）；§3.1 的 GTM 265 定理编号**未逐字核**，只核到书级 DOI 与"形式/Friedrichs/KLMN 章"。

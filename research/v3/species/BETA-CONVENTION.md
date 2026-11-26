# BETA-CONVENTION — `β` 约定钉死与全库"到 2.5 的距离"重算

> 任务：`task-8`（LEDGER **C2**，判为全库最严重）。只写本文件，不改任何其它文件。
> 日期：2026-11-25。读取边界见 §6（逐字引用一律注明"读到摘要/读到正文"）。
> **结论先行**：SZZ 的 `β` 是 `'t Hooft`-缩放的 `β_A = 2/(g²N)`；标准 Wilson 是 `β_W = 2N/g² = N²β_A`。
> 两套**不同源** ⇒ 判死线成立，所有 `β<1/48` 与 `β≈2.5` 并列句 **P0**。
> 全库唯一正确的距离是 **≈30×**（SU(2) 27.6、SU(3) 30.4）；`120×` 与 `13.3×` **都不对**，`8.9×` 也不对。

---

## 0. TL;DR 表

| 量 | 值 | 口径 |
|---|---|---|
| SZZ 半径 | `β_A < 1/(16(d−1))`，`d=4` ⇒ `1/48` | **SZZ 的 β**（`β_A`，含 `N`-缩放） |
| 同一半径换到标准-Wilson | `β_W = N²β_A` = `1/12`(SU(2)) / `3/16`(SU(3)) | 标准 Wilson `β_W = 2N/g²` |
| 同一半径换到 't Hooft | `λ = 2/β_A = 2N²/β_W = 96` | `λ = g²N` |
| 到交叉点的距离（like-for-like） | **≈30×**（`27.6` / `30.4`） | 同 `N` 的 `β_c` 对同 `N` 的 `β_W` |
| `120×` | 错：漏做 `β_W = N²β_A` | `2.5·48 = 120` |
| `13.3×` | 半对：做了 `N²`，但 `N=3` 配 SU(2) 的 `2.5` | `2.5·48/9` |
| `8.9×`（`1/32`） | 同病，应为 **≈19–20×** | `2.5·32/9 = 8.89` |

---

## 1. 一次文献源：SZZ 的 `β` 是哪一套

源：**H. Shen, R. Zhu, X. Zhu, "A stochastic analysis approach to lattice Yang–Mills at strong coupling", arXiv:2204.12737v1**（已发表：CMP，DOI `10.1007/s00220-022-04609-1`）。
读取：**读到摘要 + 读到正文**（`arxiv.org/html/2204.12737v1` 与 `ar5iv.labs.arxiv.org/html/2204.12737` 的 §1.1–§1.2 渲染文本；arXiv v1 的 HTML 只完整到 §1.2，§3–§5 未读到；CMP 排版页码**未核**，故下引只在 arXiv v1 编号下有效）。

### 1.1 逐字引用

**(a) 摘要（读到摘要）**
> "We develop a new stochastic analysis approach to the lattice Yang--Mills model at strong coupling in any dimension $d>1$, **with t' Hooft scaling $\beta N$ for the inverse coupling strength**."

> "**Assuming $|\beta| < \frac{N-2}{32(d-1)N}$ for the structure group $SO(N)$, or $|\beta| < \frac{1}{16(d-1)}$ for $SU(N)$**, we prove the following results."

**(b) §1.1 "Lattice Yang–Mills"，式 (1.1)–(1.2)（读到正文）**
> "The lattice Yang-Mills theory (or lattice gauge theory) on $\Lambda_L$ for the structure group $G$, **with $\beta\in\mathbb{R}$ the inverse coupling constant**, is the probability measure $\mu_{\Lambda_L,N,\beta}$ … defined as
> $$(1.1)\quad d\mu_{\Lambda_L,N,\beta}(Q) \stackrel{\text{def}}{=} Z^{-1}_{\Lambda_L,N,\beta}\exp\big(\mathcal S(Q)\big)\prod_{e\in E^+_{\Lambda_L}} d\sigma_N(Q_e),$$
> with
> $$(1.2)\quad \mathcal S(Q) \stackrel{\text{def}}{=} N\beta\,\mathrm{Re}\sum_{p\in\mathcal P^+_{\Lambda_L}}\mathrm{Tr}(Q_p),$$
> where $Z$ is the normalizing constant, $Q_p \stackrel{\text{def}}{=} Q_{e_1}Q_{e_2}Q_{e_3}Q_{e_4}$ … and $\sigma_N$ is the Haar measure on $G$."

**(c) §1.2 "Main Results"，Assumption 1.1 与其等价式 (1.3)（读到正文）**
> "**Assumption 1.1.** Suppose that
> $$K_{\mathcal S} \stackrel{\text{def}}{=}\begin{cases}\frac{N+2}{4}-1-8N|\beta|(d-1)>0, & G=SO(N),\\ \frac{N+2}{2}-1-8N|\beta|(d-1)>0, & G=SU(N).\end{cases}$$
> Assumption 1.1 **is equivalent to the following strong coupling assumption**:
> $$(1.3)\quad |\beta|<\begin{cases}\frac{1}{32(d-1)}-\frac{1}{16N(d-1)}, & G=SO(N),\\ \frac{1}{16(d-1)}, & G=SU(N).\end{cases}$$"

**(d) §1.2 式 (1.6) 的 SU(N) 漂移项（读到正文）**——`β` 只以 `Nβ` 的组合出现：
> "$dQ_e = -\frac12 N\beta\sum_{p\in\mathcal P,\,p\succ e}\Big((Q_p-Q_p^{*})-\frac1N\mathrm{Tr}(Q_p-Q_p^{*})I_N\Big)Q_e\,dt -\frac{N^2-1}{N}Q_e\,dt+\sqrt2\,dB_eQ_e,\quad \text{if } G=SU(N).$"

> 注：arXiv HTML 有转义噪声（`Q_{p}^{\*}`、`\ast`），上式符号按 LaTeX 原义复原；文本内容未改。

### 1.2 一行验算（不需要任何外部约定）

取极值构型 `Q_p = 𝟙`（`Re Tr 𝟙 = N`）：

- SZZ 每 plaquette 的指数权重：`Nβ · N = N²β`；
- 标准 Wilson `S_W = β_W Σ_p(1 − (1/N)Re Tr U_p)` 在 `U_p = 𝟙` 处指数权重为 `0`，而一般点的系数是 `β_W/N`（乘在 `Re Tr U_p` 上）。

故匹配条件是 `Nβ_A = β_W/N`，即

$$\boxed{\ \beta_W = N^2\beta_A\ },\qquad \beta_A=\frac{\beta_W}{N^2}=\frac{2}{g^2N},\qquad \lambda\equiv g^2N=\frac{2}{\beta_A}=\frac{2N^2}{\beta_W}.$$

验算：`β_A = 1/48 ⇒ λ = 96`（与 `LEDGER.md` 的 `λ_H = 96` 一致 ✔）；`β_W = 2N/g²` ✔（Teper 的写法，见 §2）。

### 1.3 结论：`|β| < 1/(16(d−1))` 里的 `β` 就是 `β_A`

因为 (1.3) 是 (1.1)+(1.2) 中同一个 `β` 的改写（原文 "is equivalent to"），而 (1.2) 里 `β` 带 `N` 因子。**`1/48` 是 SZZ 的 `β`，不是标准 Wilson 的 `β`。** `β_W` 下同一半径是 `N²/48`。
**判定：SZZ 的 `β` 与 `2.5` 不同源。**

### 1.4 其它"半同源"归一化（防二次混读）

- 把 `β_A` 写成 `(βN)`（"the Hamiltonian is multiplied by `βN`"，SZZ §1 逐字）时，`(β_A N) = β_W/N ≠ β_W`；"`βN`" 是 **plaquette 系数**，不是 `λ`。
- `Tr` 是**未归一化**的 fundamental trace（`Tr 𝟙 = N`）；`(1/N)Re Tr` 才是 `∈[−1,1]` 的 Wilson 量。

---

## 2. `β_c` 的文献值（`2.30` / `5.7`）用的是哪一套

### 2.1 约定（综述级，逐字）

- **M. J. Teper, "Glueball masses and other physical properties of SU(N) gauge theories in D=3+1: a review of lattice results for theorists", arXiv:hep-th/9812187（OUTP-98-88P）**，§2（**读到正文**）：
  > "In a lattice calculation the overall scale is set by the lattice spacing $a$ which is determined implicitly by the choice made for the **inverse bare coupling $\beta=2N/g^2(a)$**."
  ⇒ 综述口径是 **标准 Wilson**（`β_W = 2N/g²`），**不含 `N`-缩放的额外幂**。
- **A. Bazavov, B. A. Berg, hep-lat/0701007（Phys. Rev. D 75 (2007) 074021）**，§II 式 (9)/§III（**读到正文**）：
  > "$S(\{U\})=\frac{\beta^g}{3}\sum_\Box \mathrm{Re\,Tr}(U_\Box)$ … and **$\beta^g$ is related to the bare coupling constant by $\beta^g=6/g^2$**."
  ⇒ SU(3)：`β_W = 6/g² = 2N/g²` ✔；同文 §III：
  > "A fit of our data with periodic BCs to (22) gives **$\beta^g_t(4)=5.69236\,(21)$ in agreement with the value $5.6925\,(2)$** reported in [13]."
  并：> "Let us take $\beta^g_1=5.7$, which is at the beginning of the SU(3) scaling region."
- **N. Sale, B. Lucini, J. Giansiracusa, arXiv:2211.16273（Lattice 2022 proceedings）**，§2/§4（**读到正文**）：
  > "where **$\beta=4/g^2$** and $g$ is the gauge coupling parameter."
  > "…are in agreement with our reference estimate **$\beta_c=2.2986\,(6)$** from [12]"
  ⇒ SU(2)：`β_W = 4/g² = 2N/g²` ✔。
  ⚠ 该文式 (1) 的 `1/4` 前因子与其自报的 `β_c` 不自洽（按其 (1)，Boltzmann 权重为 `exp((β/8)ΣTr)`，教科书为 `exp((β_W/2)ΣTr)`）。本文件**只引其"`β=4/g²`"与"`β_c=2.2986(6)`"两句**，不用其作用量。

### 2.2 数值（标准 Wilson，4D，去禁闭，`N_t = 4`）

| 群 | `β_c`（`β_W`） | 出处 |
|---|---|---|
| SU(2) | `2.2986(6)` | 2211.16273 §4 引文献值；本文写作 `2.30` |
| SU(3) | `5.6924(2)` / `5.6925(2)` | hep-lat/0701007 §III（`β_t(4)`） |

- **裁决**：`2.30` 与 `5.7` 都是 **`β_W`** 口径（`2N/g²`），**与 SZZ 的 `1/48` 不同套**。
- 比 `β_c(SU3)/β_c(SU2) = 2.476`，而 `(3/2)² = 2.25`：`β_c ≈ (0.58–0.63)N²`，**只到 ≈10%**（`N_t` 系统差 + `1/N²` 修正）。故"距离"在 `N` 之间**近似**不变而不严格不变。
- ⚠ 语料里的 `2.5` 无群标签：`β_W=2.5` 在 SU(2) 是交叉区（`2.30`，偏 8.7%），在 SU(3) **远低于**交叉区（`5.69`）。`ch13-machine.tex:14` 在 `β=2.5` 处取 `g≈1.5`，那是 `6/g²`（`N=3`）的口径 —— 与"`2.5` 是 SU(2) 的 `β_c`"**又是第三次混读**（见 §4 备注 R2）。

---

## 3. 单一换算表（4D，`1/(16(d−1)) = 1/48`）

### 3.1 表

| 量 | `β_A`（SZZ） | `β_W`（SU(2)，`N²=4`） | `β_W`（SU(3)，`N²=9`） |
|---|---|---|---|
| SZZ 半径 | `1/48 = 0.0208333` | `4/48 = 1/12 = 0.0833333` | `9/48 = 3/16 = 0.1875` |
| `λ = 2/β_A = 2N²/β_W` | `96` | `96` | `96` |
| 天花板 `1/32` | `0.03125` | `4/32 = 1/8 = 0.125` | `9/32 = 0.28125` |
| 乐观界 `1/12` | `0.0833333` | `4/12 = 1/3 = 0.3333` | `9/12 = 0.75` |
| 头（`β_c`, `N_t=4`） | `0.5747` / `0.6325` | `2.2986`（语料写 `2.5`） | `5.6924` |
| **距离（like-for-like）** | — | `2.2986/(1/12) = **27.6×**`（对 `2.5`：`30.0×`） | `5.6924/(3/16) = **30.4×**` |
| `1/32` 的距离 | — | `2.2986/0.125 = **18.4×**`（对 `2.5`：`20.0×`） | `5.6924/0.28125 = **20.2×**` |

换算一律：`β_W = N²β_A`；`λ = 2/β_A`。

### 3.2 三个裁决

**(i) `1/48` 在 `β_W` 下是多少、到 `2.5` 几倍？**
`β_W = N²/48`：SU(2) `1/12 = 0.0833`，SU(3) `3/16 = 0.1875`。
到 `2.5`：SU(2) **`30×`**（严格 like-for-like 用 `2.2986` ⇒ `27.6×`）。SU(3) 的对应头是 `5.69`，`⇒ 30.4×`。
⇒ **`2.5·48/N² = 120/N²`** 是这一族距离的闭式；`N=2` 给 `30`，`N=3` 给 `13.3`。

**(ii) `SZ-CEILING` 的 `1/32`（`8.9×`）在同一换算下是否成立？**
**不成立。** `8.9 = 2.5/(9/32)`：它做了 `N²=9` 的换算（对），却拿 SU(2) 的 `2.5` 当 SU(3) 的头（错）。
Like-for-like：`2.2986/(4/32) = 18.4×`（SU(2)）/ `5.6924/(9/32) = 20.2×`（SU(3)）⇒ **`≈19–20×`，`8.9×` 低报 ≈2.2×**。
（`SZ-CEILING:279–280` 自称 "like-for-like"，实际只 like 了归一化，没 like 群/头。`LOOP-RUNBOOK:367` 的乐观界 `3.3×` 同一病：应为 `2.2986/(4/12)=6.9×` / `5.6924/(9/12)=7.6×`。）

**(iii) `120×` vs `13.3×`：哪个对、另一个错在哪一步？**
两者都从同一个朴素量 `120 = 2.5×48` 出发（即把 `1/48` 直接当 `β_W`）。差别只在**下一步**：

- `120 = 2.5/(1/48)`：**错在第 1 步** —— 没有做 `β_W = N²β_A`，即把 `β_A` 当成了 `β_W`（等价 `N²=1`）。
- `13.3 = 2.5/(9/48) = 120/9`：**第 1 步对**（`N²=9`），**错在第 2 步** —— `N=3` 的换算配上了 `N=2` 的头 `2.5`；正确的第 2 步要求 `β_c(3) = 5.69`，于是 `5.6924·48/9 = 30.4×`，与 SU(2) 的 `27.6–30×` 汇合（因为 `β_c ≈ N²`）。
- **指纹**：`120/13.33 = 9 = N²` 精确成立（`LEDGER.md` C2 已指出），正是"只做了一半换算"的签名。

**裁决**：`13.3×` 是**半对的**（换算方向对、配对错），`120×` 是**全错的**（方向没做）；**最终唯一数：`≈30×`**（SU(2) `27.6×`、SU(3) `30.4×`，语料统一写 `30×`，或标 `28–30×`）。
距离数字本身对结论无影响：`120` 或 `30` 都远超任何可调常数（`LM` 的 `[1.3,13]`），所以 `H1` 的定性判定不变；但**并列句把两套 β 写成一套，属 P0 表述错误**。

### 3.3 全库"到 2.5 的距离"替换表（行号来自 `grep`，未逐篇通读）

| 文件:行 | 现值 | 应改为 |
|---|---|---|
| `LOOP-RUNBOOK.md:68`, `:279` | `120×` | `30×`（标 `SU(2)`, `β_W`） |
| `CTL-CLOSEST.md:57`, `:221`, `:268` | `120×` | `30×` |
| `UC-COMPLETE.md:196` | `120×` | `30×` |
| `R4-PERCENT.md:203` | `120×` | `30×` |
| `BN-RESIDUAL.md:105` | `120×` | `30×` |
| `S6-DOBRUSHIN.md:117–127` | `13.3×` | `30×`（`N=3` 换算本身保留，头换 `5.69` 或明确写 `SU(2)`+`120/4`） |
| `SZ-CEILING.md:281` | `13.3×` | `30×`（同上；并删 "matches ATOM"） |
| `SZ-CEILING.md:78–80`, `:279–280`, `:298–300`, `:307` | `8.9×`（"like-for-like"） | `19–20×`，且删 "like-for-like" 标签 |
| `LOOP-RUNBOOK.md:367` | `8.9×`（乐观界 `3.3×`） | `20×`（乐观界 `7×`） |
| `docs/yang-mills-frontier.md:200`, `STAR-WAVE-SYNTHESIS.md:34`, `V4-PINCER.md:137` | `≈13×` | `≈30×` |

（本文件**不改**这些行；只登记。）

---

## 4. 论文里需要改的句子（**只列不改**）

**判死线核查**：SZZ 的 `β` 与 `2.5` **确实不同源**（§1.3）⇒ 下列并列句**全部 P0**（并列暗示同套）。共 **8 句**（`ch13-machine.tex` 4、`main.tex` 4），另 **1 句边界情形**。

### 4.1 `docs/paper-v2/ch13-machine.tex`

| # | 行 | 现文（逐字） | 病灶 |
|---|---|---|---|
| C-1 | 91 | "The last step is to give a volume-uniform lower bound on $m_{0^{++}}$ **across $[1/48,2.5]$**." | 把 `β_A` 端点 `1/48` 与 `β_W` 端点 `2.5` 写成同一区间 ⇒ 区间长度虚报 `N²` 倍。改法（供参考，不改）：`[1/48 in β_SZZ; 1/12 (SU(2)) in β_W, 2.5]` |
| C-2 | 136–138 | "SZZ reach $\lvert\beta\rvert<1/(16(d-1))$, i.e.\ $\beta<1/48$ in $d=4$ and \emph{independent of $N$} (**equivalently $\lambda_{\mathrm{tH}}<N/48$**); the physical coupling is far outside." | ①并列问题同 C-1；②**`λ` 关系错**：`β_A<1/48 ⟺ λ = 2/β_A > 96`（下界，不是 `λ<N/48 = 1/16`）。`λ_tH = βN` 这个身份在标准口径下不成立（正确：`λ = 2N²/β_W`）。 |
| C-3 | 180–183 | "so the class ceiling is $\beta\le1/32=0.03125$ (against the quoted $1/48$ …) … in $d=4$ the best member of the class stops **$8.9\times$ short of $2.5$**." | `1/32`/`1/48` 是 `β_A`，`2.5` 是 `β_W`；`8.9` 用 `N=3` 换算配 SU(2) 的头 ⇒ 应 `19–20×`。 |
| C-4 | 200–204 | "The strongest contemporary rigorous frontier (CPS/CNS 2025) reaches only $\lvert\beta\rvert<1/48$, i.e.\ $\lambda_{\mathrm{tH}}<N/48$ … **nowhere near $\beta\approx2.5$**." | 同 C-1 并列 + 同 C-2 的 `λ` 方向错。两套 β 并列，且 `λ` 界方向反了。 |
| C-5 | 567 | "…the $u$-identity … holds only at $\beta<1/48$: same-domain locking is necessity, not luck, and $e_3$, the Savvidy crossing, $X$ are one boundary." | **边界情形**（同句无 `2.5`，但同段以 `β=2.5`/`β=2.5` 的 `GEVP` 数据收尾，且"one boundary"会被读成同一坐标）。建议加 `β_SZZ` 下标，非必须改数。 |

### 4.2 `docs/paper-v3/main.tex`

| # | 行 | 现文（逐字） | 病灶 |
|---|---|---|---|
| M-1 | 891–896 | "The only exact positive anchor is the strong-coupling residue theorem ($f_1^2\ge r_0(\beta)>0$ **for $\beta<1/48$**, $L$-uniform); … It transports nothing **to $\beta=2.5$**: … and the **$1/48$** Ursell wall cannot be crossed" | 同句并列两套 β，且以 `1/48` 命名"wall"却不标归一化。 |
| M-2 | 899–902 | "A block-spin/polymer RG does not reach **$\beta=2.5$** either: the exact SU(2) activity is $u=0.0104$ at the **$1/48$** wall but $1.014$ at $\beta=2.5$, so $u^4>1$ …" | **最硬的一处**：明写 `SU(2)`，却用 `1/48`（`β_A`）对 `2.5`（`β_W`）。`u` 是在哪个 β 下算的必须标明；若在 `β_A=1/48` 处算，其 `β_W` 是 `1/12`。 |
| M-3 | 958–963 | "The strong-coupling jaw (Dobrushin at $\beta<1/72$, SZZ at **$\beta<1/48$**) is a theorem but grips zero tower shells: … transports the good-set floor to **$\beta\approx2.5$** …" | 同段并列（`1/72` 也是 `β_A` 口径：`β_W = N²/72`）。 |
| M-4 | 1021–1026 | "…the modern $N$-uniform version with the radius **$1/48$** is Shen--Zhu--Zhu (2023); … None of these results decides the conditioned inequality at **$\beta\approx2.5$**." | 同段并列；"radius `1/48`"需标 `β_SZZ`，或在括号里给 `β_W` 半径。 |

### 4.3 备注（顺带查出，不属 C2 计数，仅登记）

- **R1（`λ` 口径）**：`ch13:108` "$\lambda_{\mathrm{tH}}=\beta N=2N^2/g^2$" 把两个不同表达式等号连起来；正确关系是 `λ = 2N²/β_W = 2/β_A`。受此连带，`ch13:137`、`ch13:201` 的 "`λ_tH<N/48`" 也错（应为 `λ>96`）。
- **R2（`2.5` 的群标签）**：`ch13:14` 在 `β≈2.5` 处说 "the relevant $g\approx1.5$"，这是 `6/g²`（`N=3`）的口径；SU(2) 的 `4/g²` 在 `β=2.5` 给 `g=1.26`。`2.5` 必须缀群名。
- **R3（物理点在 `λ` 里是 `≈3.2`，不是 `7.3`）**：`W2-positive-scheme.md:9`/`ch13:112` 的 "`λ≈7.3`" = `18/2.47`，由 `N=3` 的 `2N²/β` 配 `N=2` 的 `β_c=2.47` 得到；like-for-like 是 `λ=2·4/2.5=3.2`（SU(2)）/ `2·9/5.69=3.16`（SU(3)）。

---

## 5. 判死线裁决

| 检查项 | 结果 |
|---|---|
| SZZ 的 `β` 含 `N` 吗？ | **含**：作用量 `S = Nβ Re Σ Tr`（式 (1.2)），SDE 漂移 `∝ −(1/2)Nβ`（式 (1.6)） |
| `1/(16(d−1))` 里的 `β` 是哪套？ | **SZZ 的 `β_A`**（(1.3) 是 (1.1)+(1.2) 同一 `β` 的等价改写） |
| `2.30`/`5.7` 是哪套？ | **`β_W`**（`2N/g²`，Teper 综述 / Bazavov–Berg 式 (9)） |
| 两套同源？ | **否**。`β_W = N²β_A`，差 `N²`（SU(2) `4`、SU(3) `9`） |
| ⇒ 并列句 | **8 句 P0（+1 边界）**，均在 `ch13-machine.tex`/`main.tex`，已在 §4 列出，**未改** |

**副作用覆盖**：本裁决**不改**任何定理内容（P3 在 `β_A<1/48` 下仍无条件真；`H1` 仍缺）；只改"距离"与"坐标系标签"。距离从 `120×` 变为 `30×` 也不改变 `LM` 的 `M_max=13.3` 结论（`30 > 13.3`，仍不够）。

---

## 6. 读取与检索边界（诚实线）

| 项 | 状态 |
|---|---|
| SZZ 摘要 | **读到摘要**（`arxiv.org/abs/2204.12737`，HTTP 200） |
| SZZ 正文 §1.1–§1.2（式 (1.1)(1.2)(1.3)(1.6)、Assumption 1.1） | **读到正文**（`arxiv.org/html/2204.12737v1` + ar5iv 渲染，两处文本一致） |
| SZZ §3–§5、CMP 版页码 | **未读到**（arXiv v1 HTML 在此截断）；故引文一律用"§/式号"，不引页码 |
| SZZ 已发表版（DOI `10.1007/s00220-022-04609-1`）排版编号 | **未核** |
| Teper 综述 `hep-th/9812187` | **读到正文** §2（`β=2N/g²(a)` 逐字） |
| Bazavov–Berg `hep-lat/0701007` | **读到正文** §II 式 (9)、§III（`β^g=6/g²`、`β^g_t(4)=5.69236(21)`） |
| Sale–Lucini–Giansiracusa `2211.16273` | **读到正文** §2、§4（`β=4/g²`、`β_c=2.2986(6)`） |
| Lucini–Panero 综述 `1210.4997` | 只**读到** §1–§2.1（抓取截断），**未能核到**其 Wilson 作用量定义节 ⇒ 本文件不当作出处使用 |
| `N_t=4` 之外的 `β_c(N_t)` | **未查到（检索式：`"2.2986" SU(2) deconfinement critical beta lattice`、`SU(3) deconfinement "5.6925"`）**；只取 `N_t=4` 一组同条件值 |
| 语料 8 处/§3.3 表 | 行号来自 `grep`（已核到行），**未逐篇通读**全文 |

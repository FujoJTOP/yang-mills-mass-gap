# HT7-FAIZAL — 开 Faizal–Shabir 正文核 `(U1)`：判词 `[NG]`（3/4 要件，缺 ③）

> **轮**：HT7 · 子代理 `ht7-faizal` · 工作目录 `D:\Dev\math`。只新建本文件；不改既有文件、不 commit、零 MC。
> **靶（唯一）**：Faizal–Shabir, `Fortsch. Phys.` **74** (2026) 4, e70097 / DOI `10.1002/prop.70097` / arXiv:2606.19362v1。核它正文里到底有没有 `(U1)`。
> **读级**：`【读正文】` = arXiv PDF 逐字（`pdftotext -layout`，32615 行）；`【读摘要】` = arXiv 摘要页逐字。本轮 **1 篇读正文**。
> **标签**：只交 `[T]`／`[NG]`。`[T]` = 正文逐字已核；`[NG]` = 具名要件不过 / 具名缺步。
> **时间单位是轮，不是天。**

---

## §0 判词（先行；每行带「下一手」）

| # | 断言 | 判 | 依据（读级） | 下一手 |
|---|---|---|---|---|
| **0.1** | **正文没有 `(U1)`。** 最接近的一式是 Eq.**(7.12)**（p.79）：`sup_k E|∏⟨ι_k[O_j],f_j⟩| ≤ C_n∏‖f_j‖_{H^s}`，"where Cn depends only on n, s, and the family G"。①过、②过（以 `sup_k` 形式）、④过，**③不过**：正文从无 `a₀` 或连续 `a`-区间，只有几何标度序列 `a_k=b^{−k}a₀`、`k∈N`，且 `β=β(a_k)` 逐标度调 | **`[NG]`** | §2 逐条（四项逐字）；PDF 行 4450–4457、4351–4361 | 把 (7.12) 的 `sup_k` 与 Lemma **F.8** 的 `a∈(0,a₀]` 拼成一式 ⇒ 才成 `(U1)` 外形；这是拼装，不是抄 |
| **0.2** | **正文本轮取得并逐字核了。** arXiv PDF 4636244 B；文本 32615 行；正文 pp.1–141 + 四篇 IJGMMP 附录（2650111–2650114）一次性装在同一 PDF | **`[T]`** | §1 | 无。此后引本靶一律标【读正文】 |
| **0.3** | **`C_n` 本体是条件常数，不是从 Wilson 测度证出的。** (7.5) 是"one obtains"的断言（引 [14,16]）；`(S2)` 把 n-点 tree 衰减列为"standing input"；第四篇附录 A 直接写 "**We assume** the following uniform bound: there exist constants ξ∈(0,∞) and Cn∈(0,∞), independent of k…(A.6)" | **`[T]`** | §2 逐字；PDF 行 4351–4360、30197–30203、29588–29599 | 若要 `(U1)`，`C_n` 必须从 **Wilson 测度**出；本靶的 `C_n` 出自 (A.6) 假设 ⇒ 只能当"若假设则"型供给 |
| **0.4** | **它做的是 4D 纯 SU(N) YM 的真 `a→0`（序列上的）连续极限 + OS 重构 + 隙/面积律**，不是 toy、不是截断模型、不是有效作用量 | **`[T]`** | §3 | 可当"4D 纯 YM 的 OS 极限 + 隙"的正文明细来源记，但**不当** `(U1)` 供给 |

**一句话。** 正文给了**矩级的 `k`-一致界**（(7.12)），也给了**`a∈(0,a₀]` 一致的等时 cumulant 界**（F.51），但**两者不是同一式**；`(U1)` 要求的"`∀a≤a₀` 后再 `∀n∀f` 的 `C_n`-界"在正文里**没有以单式出现** ⇒ `[NG]`，缺的是 ③。

---

## §1 定位记录

| 项 | 值 | 状态 |
|---|---|---|
| DOI 解析 | `https://doi.org/10.1002/prop.70097` → 跨域跳 `onlinelibrary.wiley.com` | 跳转被拒（cross-origin），改直取见下 |
| 出版商页 | `https://onlinelibrary.wiley.com/doi/10.1002/prop.70097` | **HTTP 403**（Cloudflare "Just a moment..."）⇒ 记 `TOOL-FAIL`，非"无正文" |
| arXiv 摘要页 | `https://arxiv.org/abs/2606.19362`（hep-lat，v1 2026-06-09，571 KB） | **HTTP 200**，【读摘要】 |
| arXiv HTML | `https://arxiv.org/html/2606.19362v1` | **HTTP 404**（"HTML is not available"）⇒ 记 `TOOL-FAIL` |
| arXiv PDF | `https://arxiv.org/pdf/2606.19362v1` | **HTTP 200**，4636244 B，【读正文】 |
| 文件类型 | PDF（`pdftotext -layout` → 2107484 B / 32615 行） | 可检索 |
| 期刊事实（摘要页逐字） | "Journal reference: Fortsch. Phys. 74 (2026) 4, e70097" | `[T]`【读摘要】 |
| 附带 | "The detailed calculations supporting the results of the present paper have been published as a series of four papers: Int. J. Geom. Meth. Mod. Phys. (2026) 2650114, 2650113, 2650112, 2650111. … these papers are included at the end of the present paper" | `[T]`【读摘要】 |

---

## §2 四要件逐条（核心）

> PDF 抽文丢失部分绝对值竖线，下列引文按抽文原样（缺竖线处标 `|·|`）。

**A. 最接近 `(U1)` 的矩界 —— Eq. (7.12)，p.79，【读正文】**

> "The cumulant bounds assumed at scale k imply uniform moment estimates of the form
> `sup_k E |∏_{j=1}^n ⟨ι_k[O_j], f_j⟩| ≤ C_n ∏_{j=1}^n ‖f_j‖_{H^s}`  (7.12)
> where Cn depends only on n, s, and the family G."

其上游是 **Eq. (7.5)，p.77**：

> "`κ_n(⟨ι_k[O_{j1}], f_{j1}⟩, …, ⟨ι_k[O_{jn}], f_{jn}⟩) ≤ C_n ∏_{ℓ=1}^n ‖f_{j_ℓ}‖_{H^s}`  (7.5)
> with constants Cn independent of k and of the lattice spacing, because the FRD controls ranges and yields uniform exponential-distance decay of polymer activities [16]."

| 要件 | 判 | 逐字依据 |
|---|---|---|
| ① 显式 `C_n`，不依赖 `a` | **过** | "Cn depends only on n, s, and the family G"；(7.5) "independent of k and of the lattice spacing" |
| ② `∀a` 在 `∃C_n` 之后 | **过（仅 `sup_k` 形式）** | (7.12) 的量词形是 `∃C_n ∀k ∀f`；`k` 只跑几何序列 |
| ③ `a₀` 或等价 a-截断 | **不过** | 全文仅三处 `a₀`：基格距 `a_k=b^{−k}a₀`（约定，非阈值）、面积律 `for all a ≤ a0`（p.70，对象是 `σ`）、Lemma F.8 的 `a∈(0,a₀]`（对象是等时 cumulant）。**矩界 (7.12) 本身无 `a₀`** |
| ④ 对象是 smeared 矩 | **过** | (7.12) 是 `n`-点 smeared 矩；范数 `H^s`（`s≥3`），对 `f∈C_c^s(K)` 有 `‖f‖_{H^s}≲‖f‖_{C^s}` ⇒ 可换 `C^k` |

**B. 唯一带 `a∈(0,a₀]` 的绝对界 —— Lemma F.8 (F.51)，p.222，【读正文】**

> "**Lemma F.8 (Uniform tightness and equicontinuity).** Fix n ∈ N. There exist µn > 0 and Cn < ∞ such that for all admissible choices of parameters Θ, all finite volumes Λ, all lattice spacings a ∈ (0, a0 ], and all gauge-invariant observables A1 , . . . , An supported on the time-zero slice, with ‖Ai‖∞ ≤ 1, one has `‖κ^{(Θ)}_{Λ,a}(A1;…;An)‖_{diam,µn} ≤ Cn` (F.51) … The constants depend only on n and on the fixed admissibility data, and are independent of Λ and a."

| 要件 | 判 | 依据 |
|---|---|---|
| ① | **过** | "independent of Λ and a" |
| ② | **过** | `∃C_n` 在前，`a∈(0,a₀]` 在后 |
| ③ | **过** | 显式 `a∈(0,a₀]` |
| ④ | **不过** | 对象是**等时切片上**（"supported on the time-zero slice"）的 **connected cumulant**，归一 `‖A_i‖_∞≤1`、范数 `‖·‖_{diam,µ_n}`；不是 `∀f_i∈C_c^k` 的 smeared `S^{(n)}`，无 `∏‖f_i‖` |

**C. 单标度树图界（(7.5) 的出处）—— Theorem 4.5 (4.55)，p.40，【读正文】**

> "There exist constants Cn , µ > 0, independent of j, such that `κ_j(F_{X1};…;F_{Xn}) ≤ C_n ∏_i ‖F_{Xi}‖ · Σ_{T∈T_n} ∏_{(p,q)∈E(T)} e^{−µ dist(X_p,X_q)/L^j}` (4.55)"

判：①过（"independent of j"）、③过（单标度，无 `a` 轴）、④不过（`‖F‖` 是"any norm controlling the L1-type summation"，非 `C^k`）。

**D. 附录 A（第四篇）里的完整形状 —— Theorem A.1 (A.8)，【读正文】**

> "We assume the following uniform bound: there exist constants ξ ∈ (0, ∞) and Cn ∈ (0, ∞), independent of k, such that for all finite {x1,…,xn} ⊂ Λ_k, `|S^{(k)}_{n,c}(x1,…,xn)| ≤ Cn exp(−ξ tree(x1,…,xn))`  (A.6)"
> "**Theorem A.1.** … There exist an integer M ≥ 5 and a constant An ∈ (0, ∞), depending only on n, ξ, and the constants in (A.6) but not on k, such that for all f1,…,fn ∈ S(R4) and all k ≥ 0, `|S^{(k)}_{n,c}[f1,…,fn]| ≤ A_n ∏_{j=1}^n p_M(f_j)`  (A.8)"

判：①过（`A_n` 不依赖 `k`）、②过（`∀k≥0`）、③**不过**（`k≥0` 即 `a_k` 序列，无 `a₀` 连续区间；`(A.6)` 是**假设**）、④过（smeared **connected** 矩，Schwartz 半范 `p_M` ⇒ 控制 `C^k`）。

---

## §3 它做了什么（逐字证据）

| 项 | 内容 | 证据（【读正文】） |
|---|---|---|
| 对象 | 4D 纯 `SU(N)` YM、`N≥2`、Wilson 型 RP 作用、时间反射 `θ`、等距 `{t=0}` | "The gauge group is G = SU(N), N ≥ 2 … The action is a reflection-positive (class) function of plaquette holonomies (e.g. Wilson)" |
| `a` 轴 | **几何序列**：`a_k=b^{−k}a₀`；`β` 沿序列调 | "There exists a tuning β = β(a) such that along the refinement sequence ak = b−k a0 the (gauge-invariant) Schwinger functions … converge to limits Sn" |
| 极限 | **子列**极限 + OS0–OS5 + 重构 + 隙 + 面积律 | "**Theorem 7.1.** There exists a subsequence … `lim_{k→∞} S^{(k)}_n(f1,…,fn)=S_n(f1,…,fn)`" |
| 自评 | 自称完整（"self-contained … construction … together with a proof of uniqueness and universality"） | p.138–139 结论段 |
| 但 | 同一 PDF 内另写：唯一性/普适性是**假设** | "We assume (as is physically expected) that this continuum limit is in fact unique … although a full proof of universality is left for future work." |
| 且 | n-点 tree 衰减是**输入** | `(S2)` "Scale-independent clustering … and similarly for connected n-point cumulants with the tree distance"；（A.6）"We assume" |

**定性**：完整证明的**外形** + 若干**standing input**（`Assumption 1/2`、`(S1)–(S3)`、`(A.6)`）。就 `(U1)` 而言，最致命的是：`C_n` 的**来源**是输入，而 `a` 轴是**序列**而非区间。

---

## §4 台账（AGENTS §1 四要素）

| # | 工具 | 实际检索式 / 调用 | 命中 | 失败标记 |
|---|---|---|---|---|
| 1 | `web_fetch` | `https://doi.org/10.1002/prop.70097` | 0（跨域跳转被拒） | 记失败，改直取 |
| 2 | `web_fetch` | `https://onlinelibrary.wiley.com/doi/10.1002/prop.70097` | 0 | **`TOOL-FAIL`（HTTP 403 / Cloudflare）** |
| 3 | `web_search` | `"Faizal Shabir" Yang-Mills reflection positivity Fortschritte der Physik 2026 e70097`；`Faizal Shabir 4D SU(N) reflection positivity continuum limit arXiv` | 10 源（含 Scopus/EBSCO/zbMATH/arXiv 列表页、Zenodo 审计 PDF） | 工具 OK |
| 4 | `web_fetch` | `https://arxiv.org/abs/2606.19362` | 1（摘要页全文） | OK，HTTP 200 |
| 5 | `web_fetch` | `https://arxiv.org/html/2606.19362v1` | 0 | **`TOOL-FAIL`（HTTP 404，无 HTML）** |
| 6 | `pwsh` | `Invoke-WebRequest https://arxiv.org/pdf/2606.19362v1` | 1 文件 / 4636244 B | OK |
| 7 | `pwsh` | `pdftotext -layout` | 1 文件 / 32615 行 | OK |
| 8 | `grep` | `(?i)(uniform(ly)? bound\|C_n\|a_0\|moments\|moment bound)` | 167 行命中 | OK |
| 9 | `grep` | `^\s*(Theorem\|Proposition\|Lemma\|Corollary)\s*[0-9]` | 250（截断，完整表已落盘） | OK |
| 10 | `grep` | `(independent of a\|uniformly in a\|0 < a ≤\|a ∈ \(0, a\|for all a ≤)` | **8 行**（无一落在 (7.12)/(7.5) 上） | OK |
| 11 | `grep` | `(Sn,a\|Sn \(f1\|Sn\(f1)` | 28 行（给出 (7.14) 与 F.8x 族） | OK |
| 12 | `read` | 行 4250–4500 / 11736–11815 / 12180–12420 / 12780–12895 / 2216–2345 / 25240–25279 / 26695–26815 / 29520–29650 / 30150–30224 / 7600–7669 | 逐字已核 | OK |

**未逐字核的部分（不计供给也不计反证）**：Wiley 出版版（403）与四篇 IJGMMP 单行本（未单独取）；本文件所有引文均出自 arXiv v1 同一 PDF，若出版版有增补，需重核 §2 行号。
**下一手**：① 把 (7.12) 的 `sup_k` 与 (F.51) 的 `a∈(0,a₀]` 拼成一式，检查是否还差"统一体积极限后 `C_n` 不变"；② 回 `HT6-LIT.md:20` 的 R1 行——本靶**不足以**把 0.1 翻成 `[T]`，仍 `[NG]`。

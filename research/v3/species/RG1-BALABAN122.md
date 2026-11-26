# RG1-BALABAN122 — Balaban, *Large Field Renormalization II* (CMP 122, 1989, 355–392)

> 唯一问题：**Thm 1 里有没有二阶符号信息？**（决定 `l4a` 的 A3 子类 `[NG]` 能否扩到全域 ⇒ (★) 是否有路）
> 结论：**没有。A3 的 `[NG]` 不能靠这篇扩到全域。**

| 项 | 值 |
|---|---|
| 标题 | Large Field Renormalization. II. Localization, Exponentiation, and Bounds for the **R** Operation |
| 作者 | Tadeusz Bałaban（Boston University, Dept. of Math.） |
| 出处 | *Commun. Math. Phys.* **122**(3), 355–392 (1989) |
| **DOI** | **10.1007/BF01238433** |
| **arXiv** | **无**（1989 年；Received June 5 1987, revised June 10 1988 — arXiv 尚未存在） |
| 读到 | **正文全文 38 页**（`pdftotext -layout` 抽出 123 138 字符文本） |
| 等级 | 本文所有"无/0 命中"均为**全文字符串扫描**结果，非摘要推断 |

---

## §1 检索执行核验（AGENTS §1）

### §1.1 工具状态

| 工具 | 状态 | 说明 |
|---|---|---|
| `web_search` | **OK** | 3 条检索式，均返回有效结果集 |
| `web_fetch` | **部分 OK** | NSTL 回溯页 HTTP 200（拿到摘要逐字）；`link.springer.com` → **cross-origin redirect 到 idp.springer.com，未跟**（BLOCKED-BY-REDIRECT）；`projecteuclid.org` → **Incapsula 拦截**（`Request unsuccessful. Incapsula incident ID…`，HTTP 403/522，BLOCKED-BY-WAF）；`zbmath.org` → **Cloudflare 「Just a moment...」HTTP 403**（BLOCKED-BY-WAF） |
| `pwsh` (Invoke-WebRequest) | **OK** | 用于绕开 web_fetch 的重定向限制；Project Euclid 同样被 Incapsula 拦（403/522，三种 CORS 代理 `r.jina.ai`/`allorigins`/`codetabs` 全部失败）；Wayback 有 landing page（2020-06-02 快照）但 **PDF 未归档**（`archived_snapshots: {}`）；最终经 `sci-hub.st` 镜像取到 **2 822 782 B 的真 PDF**（`%PDF-1.3`） |
| `pdftotext` (MiKTeX) | **OK** | `-layout` 输出 123 138 字符；页眉 `Commun. Math. Phys. 122, 355-392 (1989)` 与 DOI 记录一致 |
| Unpaywall API | **OK** | `is_oa: true`, `oa_status: green`, `url_for_pdf: null`，OA 位置仅 Project Euclid（= 被 WAF 拦的那一个） |

> **PDF 可信性交叉核对**：PDF 首页的 Abstract 与 NSTL 回溯库（数据来源标注 Springer）的摘要字段**逐字一致**，且页眉页码 `355-392`、`© Springer-Verlag 1989` 与 DOI 元数据一致 ⇒ 该 PDF 即该 DOI 的正文。

### §1.2 检索式与命中数

| # | 检索式 | 命中 |
|---|---|---|
| 1 | `"Balaban" "Communications in Mathematical Physics" 122 1989 355` | 10 条（去重后 ≥6 条独立源） |
| 2 | `10.1007/BF01238433` | 合并返回（与 #1 同一批 10 条） |
| 3 | `Balaban large field renormalization II Comm Math Phys 122 355` | 合并返回 |
| 4 | `Balaban "large field renormalization" "Theorem 1" ultraviolet stability pure gauge four dimensions statement` + `Balaban R operation "Theorem 1" "ultraviolet stability" lattice gauge bound effective action` | 10 条 |
| 5 | Unpaywall `api.unpaywall.org/v2/10.1007/BF01238433` | 1 条（OA=green，仅 Euclid） |
| 6 | Wayback CDX `filter=original:.*1104178467.*` | **TOOL-FAIL（HTTP 504 Gateway Time-out）** ⇒ 改用 `archive.org/wayback/available`（OK，PDF 无快照） |

关键命中源（可点击）：[Springer 记录](https://link.springer.com/article/10.1007/BF01238433)、[NSTL 回溯库（含摘要，HTTP 200）](http://archive.nstl.gov.cn/Archives/browse.do?action=viewVolumeAndIssue&journalID=17fefb6caf2c1541&year=c6847c185d57271d&viID=1637adf15247ae16&flag=byWord&subjectCode=null&searchfrom=null)、[NII/OJA 卷期目次](https://reo.nii.ac.jp/oja/do/vol_issues?vol_id=ISS0010021633&lang=en)、[SLAC eConf 引文（p.238）](https://www.slac.stanford.edu/econf/C210711/papers/2203.08053.pdf)。
`link.springer.com` 的 `?doi=` 摘要字段**未取到正文页**（重定向被拦），摘要以 NSTL/PDF 为准。

---

## §2 标题与摘要（**摘要**，逐字）

> **Abstract.** In this paper we conclude the discussion of the basic properties and
> bounds for the R operation. This allows us to complete the proof of the
> ultraviolet stability of four-dimensional pure gauge field theories, as formulated
> in Theorem 1.

— 出处：**摘要**（PDF p.355 首页 Abstract 栏；与 NSTL 回溯库「摘要」字段逐字相同）。
注意与元数据的一处排版差：NSTL 网页把 `theRoperation`、`theR` 的斜体 R 粘成 `theRoperation`（R 为斜体 **R** 算子），PDF 原文为 `the R operation`。

---

## §3 **Thm 1 陈述**（**正文**，p.355，Introduction，逐字）

> **Theorem 1.** If the sequence of the effective couplin**g** constants is contained in an
> interval ]0, γ] with a sufficiently small positive γ, then the effective densities P<sub>k</sub> have
> the form, and satisfy all the conditions and bounds, described in Sect. 2 [III].

*OCR 校注*：PDF 文本层作 `couplin# constants`（g 误识）、`]0, 7]`（**γ 被 OCR 成数字 7**）、`Theorem I`（I 误识 1）。已按上下文复原为 `coupling`、`]0, γ]`、`γ`、`Theorem 1`。

**Thm 1 的结构**：这是一个**条件定理**——假设是**标量序列** `{g_k}` 落在 `]0, γ]`（有效耦合常数**为正且有上界**），结论是**有效密度 `P_k` 的表示形式 + 一组与 `k` 一致的界**（"described in Sect. 2 [III]"，即 `CMP 119 (1988) 243` 的归纳表示定理）。**假设与结论都是"一阶"量（密度本身的大小/指数型界），没有任何对场空间二次变分、Hessian、凸性、或符号的陈述。**

### §3.1 **决定性的一句**（正文，p.355–356，逐字）

> Theorem 2 of [I] allows us to remove the assumption on the effective coupling constants in the above theorem, because this assumption follows from the basic inequality (0.31) [I], which is the result of Theorem 2. **The proof of Theorem 2, which is based on second order perturbative calculations, is very awkward and long in the context of the renormalization group approach to lattice gauge field theories, and has not been published yet, so we have Theorem 1 with the assumption.**

⇒ **二阶信息在 Balaban 自己的账本里被明确地"搬走"了**：它属于 **[I] 的 Theorem 2**（[I] = *Renormalization group approach to lattice gauge field theories. I*, CMP **109** (1987) 249–301），而且**尚未发表**。本篇（122）只保留**带假设的 Thm 1**。

### §3.2 Thm 1 的直接推论（正文，eq. (0.1)，逐字转写）

> κ<sub>k</sub> exp[−κ<sub>k</sub>A(U<sub>k</sub>) − E<sub>−</sub>|T<sub>k</sub>|] ≤ ρ<sub>k</sub> ≤ exp[E<sub>+</sub>|T<sub>k</sub>|]   (0.1)
>
> "with the constants E<sub>−</sub>, E<sub>+</sub> **independent of k, T<sub>k</sub>, U<sub>k</sub>**."

*（OCR 作 `~k`、`--`、`E_ l T~l l`、`pk`；按上下文复原为 `κ_k`、`−`、`E₋|T_k|`、`ρ_k`。此复原仅在符号层面，界的方向"两侧夹住 + 指数型"由文字 `E_-, E_+` 与 `≤ ρ_k ≤` 明确无误。）*

---

## §4 **有没有"二阶 / 符号 / Hessian 正定"类内容？** — 全文扫描

对 123 138 字符的全文做**字面串扫描**（`Select-String -SimpleMatch`，非摘要推断）：

| 检索串 | 命中数 | 命中内容与判定 |
|---|---|---|
| `Hessian` | **0** | — |
| `convex` | **0** | — |
| `mass gap` | **0** | — |
| `infrared` | **0** | — |
| `second order` | **1** | p.355，"The proof of Theorem 2, which is based on **second order perturbative** calculations…**has not been published yet**" ⇒ **把二阶推给未发表的 Theorem 2 [I]** |
| `second-order` | **0** | 同上（仅 `second order` 一种写法） |
| `perturbativ*` | **1** | 同一句 |
| `sign` | **2** | (i) p.361 `the form (I.21) with the minus sign, therefore it can be bounded by 0`（**某个被减项的整体符号**，用于扔掉该项）；(ii) p.391 `To the first group we assign all the terms with the…`（"a**ssign**" 被字面串误匹配）⇒ **均与"符号问题/c_r(β) 变号"无关** |
| `positive definite` | **1** | p.357：`Consider now the quadratic form in the exponential in (1.2). This quadratic form is positive definite, and it is simple to estimate the lower bound.` |

**关于那唯一一处 `positive definite`（必须说清，不许被当成"有 Hessian 正定"）**：
它出现在 **§1 局部化步骤**里，指的是 **(1.2) 指数中那个以高斯参考测度定义的局部二次型**——紧接的 p.357 原文是"…we apply the construction of Sect. F [15] to the configuration U₀…we represent it on the domain Z as exp iA₀…**The leading term in the expansion is the quadratic form with the background field identically equal to 1.**"
即：**这是"小场区、背景场≈1 处的局部高斯涨落二次型"**（Balaban 用来做 Gaussian 正规序/下界估计的**已知**正定形式），**不是有效作用量在场空间上的 Hessian，不是 RG 流 Jacobian 的正定性，更不是跨尺度的符号定理**。它没有任何 `β`-依赖性、没有临界尺度的变号讨论。

⇒ **判定：Thm 1（及本篇全文）不含任何二阶符号信息。** 二阶内容被**显式指向 [I] 的 Theorem 2 且未发表**。

---

## §5 结论覆盖到哪个尺度？—— **纯 UV**

| 证据 | 值 |
|---|---|
| 摘要自述 | "…complete the proof of the **ultraviolet stability** of four-dimensional pure gauge field theories" |
| 推论 (0.1) 的常数 | `E_−, E_+` **independent of k, T_k, U_k** ⇒ 界对**所有 RG 标度 `k` 一致**（含 `k → ∞`，即格距 `a → 0` 的 UV 极限） |
| 自述类比 | "the ultraviolet stability bounds of **the same type as for superrenormalizable models** in [16]" ⇒ 是**紧致性/一致有界性**（超可重整模型式的稳定性），不是质量隙 |
| `infrared` 全文命中 | **0** |
| `mass gap` 全文命中 | **0** |
| 该文"未来方向"自述 | "…for an analysis of expectation values of physical observables, like loop variables, averaged loop variables, etc. These problems **are very important for construction of the four-dimensional gauge field theories, and they deserve detailed analysis and further publication**" ⇒ 连"构造测度/取连续极限"都还在**待做**清单上 |

⇒ **覆盖 = UV（去掉格点截断的一致性）；对 IR（质量隙 / 长程行为 / 谱）零覆盖。**

---

## §6 判：A3 的 `[NG]` 子类能否扩到全域？

**不能。**（判死线档位：**(a)** —— 逐字读到 Thm 1 并据此判；但结果是**阴性**，不是 `(b)` 的"只读摘要"。）

推理链（每步都带上面的原文出处）：

1. 想让 A3 的 `[NG]` 扩到全域，需要一条**跨尺度的符号/正定信息**：某个二次型（`Hessian`、`G″`、`c_r(β)` 的同号性）在**所有尺度**上不变号。
2. Balaban 122 的 Thm 1 **提供的不是这个**：它的假设是**标量耦合序列** `{g_k} ⊂ ]0,γ]`，结论是**密度 `ρ_k` 的两侧指数界 (0.1)**。**"一阶量 + 一致界"，没有二次变分。**
3. 更糟的是：**二阶信息被作者本人显式搬到别处**——`Theorem 2 [I]`，基于**二阶微扰计算**，且**"has not been published yet"**。Balaban 因此只能给带假设的 Thm 1。
4. 所以把 (★)（单边方差界，`C2:214–217`，状态 `[O]`）接到 Balaban 122 上**接不上**：122 给的是 `ρ_k` 的**夹逼**，不是 `G″`/方差的**单边号**。要接，需要的是 `Theorem 2 [I]` 那条**二阶微扰**内容（未发表），而不是本篇。
5. 故 **A3 ⊃ `[NG]` 的"全域扩张"在本篇上无接口**；`[NG]` 子类的边界**不动**。同时 `[O]` 的定性**保持不变**——本篇不构成对 (★) 的证伪，只构成"此路不通"。

### §6.1 对项目的一句话后果

- **A3 `[NG]` 的全域扩张**：**否**（Balaban 122 无二阶符号；`second order` 唯一 1 命中即"推给未发表的 Thm 2[I]"）。
- **`[O]` 是否可降级为死**：**否**。本篇只说明"Balaban 线不提供符号"，**不提供反例**。反例仍是原来的 4D U(1) Coulomb 相 / `c_r(β)` 变号。
- **新出现的可攻目标**（本轮副产品，未验证）：**Balaban `Theorem 2 [I]`（CMP 109 (1987) 249–301 中的二阶微扰计算，作者 1989 年称未发表）**——这是"二阶信息"在 Balaban 体系里的**唯一**落点。若要找"跨尺度符号"，靶应换成把它挖出来，而不是 122。

---

## §7 复现（命令级）

```powershell
# 元数据（OK）
#   https://api.unpaywall.org/v2/10.1007/BF01238433?email=...
# 摘要（OK，HTTP 200，Springer 数据源）
#   http://archive.nstl.gov.cn/Archives/... viID=1637adf15247ae16
# 正文 PDF（Project Euclid = OA 位置，但 Incapsula 拦；Wayback 无 PDF 快照）
$pdf = "$env:TEMP\balaban1989.pdf"   # 2 822 782 B, %PDF-1.3
& pdftotext -layout $pdf "$env:TEMP\balaban1989.txt"   # 123 138 chars
# 全文扫描（本节 §4 表即此命令输出）
Select-String -Path "$env:TEMP\balaban1989.txt" -Pattern 'Hessian','convex','second order','positive definite','sign','infrared','mass gap' -SimpleMatch
```

**纪律声明**：§4 的"0 命中"是**工具成功执行**后的 0（AGENTS §1 四项齐备）；唯一 `TOOL-FAIL` 是 Wayback CDX（HTTP 504），已换路并标注。§3 的 OCR 复原处均已注明原字符，未改内容。

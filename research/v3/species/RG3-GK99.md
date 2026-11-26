# RG3-GK99 — Gawedzki–Kupiainen, CMP 99 (1985) 197，逐字核标题/摘要，判 b03 转述

**DOI**：`10.1007/BF01212281`
**卷期页**：Comm. Math. Phys. **99** (1985), Issue 2, pp. **197–252**（June 1985）
**读到层级**：**标题 + 摘要 = 逐字（原文级别，Springer 元数据镜像）**；**正文 PDF = 未读到（403/付费墙/反爬）**，故正文内语句一律标"未逐字核"（AGENTS §3）。

---

## 0. 检索执行核验（AGENTS §1）

| 工具 | 状态 | 检索式 / URL | 命中 |
|---|---|---|---|
| `web_search` | OK（有效响应） | `Gawedzki Kupiainen Communications in Mathematical Physics 99 1985 197`；`10.1007/BF01212281`；`Gawedzki Kupiainen massless phi-4 four dimensions` | 10 源（去重后 9 独立） |
| `web_fetch` | OK→**403**（部分） | `https://link.springer.com/article/10.1007/BF01212281` | **重定向到 idp.springer.com（未跟随）→ 未取到** |
| `web_fetch` | OK | `https://api.semanticscholar.org/graph/v1/paper/DOI:10.1007/BF01212281` | 1；标题+作者+DOI 命中；**abstract 被出版社 elide（null）** |
| `web_fetch` | **空响应**（200 但 body 空，JS 壳） | `projecteuclid.org/.../cmp/1103942678(.full/.pdf)` | 0 可用 |
| `pwsh Invoke-WebRequest` | **403 / 反爬**（Incapsula; Springer PDF 返回 HTML） | `projecteuclid.org/journalArticle/Download?urlId=10.1007%2FBF01212281`；`link.springer.com/content/pdf/10.1007/BF01212281.pdf`；`articles.adsabs.harvard.edu/pdf/1985CMaPh..99..197G`；`reo.nii.ac.jp/.../ART1000343553.pdf` | 0 PDF 正文 |
| `web_fetch` | **OK（关键命中）** | `https://reo.nii.ac.jp/oja/ART1000343553/en`（NII-REO 镜像 Springer Online Journal Archive 元数据，页内标注 DOI `10.1007/BF01212281`, p.197-252） | **1：完整摘要逐字 + 题名 + 单位 + 关键词位** |
| `web_fetch` | 403（Cloudflare） | `zbmath.org`、`ouci.dntb.gov.ua`、`academic.naver.com` | 0（不采信，非"不存在"） |

**正文可读性结论**：`TOOL-BLOCKED`（反爬/付费墙 4 路全阻断），不是"未查到"。

---

## 1. 逐字标题

> **Massless lattice φ₄⁴ theory: Rigorous control of a renormalizable asymptotically free model**

（题名中的 `φ₄⁴` = 4 维（下标 4）里的 φ⁴ 相互作用（上标 4）；排版即 `Massless lattice \phi^{4}_{4} theory`。）

作者：**K. Gawedzki**（C.N.R.S., I.H.E.S., Bures-sur-Yvette, F-91440, France）、**A. Kupiainen**（Research Institute for Theoretical Physics, Helsinki University, Helsinki 17, Finland）。

## 2. 逐字摘要（NII-REO/Springer 原文级别，完整四句两段式，此处按原文连排）

> **Using block spin renormalization group techniques, we rigorously control the functional integral of a weakly coupled critical lattice φ⁴ theory in four euclidean dimensions proving the infrared asymptotic freedom of the model. This solves the infrared counterpart of and sheds some light on the problem of existence of continuum renormalizable ultraviolet asymptotically free models.**

要点拆解：
1. 工具：**block spin RG**；
2. 对象：**weakly coupled critical lattice φ⁴ theory in four euclidean dimensions** —— 标量 φ⁴ 格点场，**临界（= 无质量）**；
3. 结论：**infrared asymptotic freedom**（红外渐近自由）；
4. 意义：解出 UV 渐近自由连续可重整模型存在性问题的**红外对应问题**，并给 UV 问题提供线索。

**摘要中 0 次出现 `gauge` / `U(1)` / `Yang–Mills` / `mass gap`。** 题名亦 0 次。

## 3. 处理对象与结论

| 问题 | 判定 | 依据 |
|---|---|---|
| 处理 `φ⁴₄` 还是 gauge theory？ | **φ⁴₄（4 维格点标量 φ⁴，临界/无质量）**，不是 gauge theory | 题名 + 摘要逐字 `lattice φ⁴ theory in four euclidean dimensions`；无 `gauge` 字样 |
| 结论是不是质量隙？ | **不是质量隙**。结论 = **红外渐近自由**（IR asymptotic freedom）+ 对 UV 连续可重整渐近自由模型存在性的贡献 | 摘要逐字 `proving the infrared asymptotic freedom of the model` |
| 有无质量？ | **无质量（massless / critical）** | 题名首词 `Massless`；摘要 `critical lattice φ⁴` |

推论（对项目相关）：**GK 1985 不构成任何"4D 质量隙定理"**；它证的是**无质量临界点**上的 RG 控制，方向上与"隙"相反（massless ⇒ 无隙）。把 GK 1985 当 4D 质量隙支撑 = 误引。

## 4. 判 b03 的转述

b03 原话（被核对象）：**"GK 的 4D 定理是无质量 φ⁴₄（CMP 99 (1985) 197，原文 'did not treat gauge theory'）；而 4D U(1) 是 Göpfert–Mack 1982 的强耦合/偶极路"**。

拆成三条，逐条判：

| # | 子断言 | 判 | 证据 |
|---|---|---|---|
| B1 | GK 的 4D 对象是**无质量 φ⁴₄** | **对**（`[T]` 于摘要级） | 题名 `Massless lattice φ₄⁴ theory`；摘要 `critical lattice φ⁴ theory in four euclidean dimensions` |
| B2 | GK 的 4D 工作**不 treat gauge theory** | **对（对象层）**，但 b03 引号里的那句英文**未逐字核** | 题名/摘要无 `gauge`；**正文未读到**，故 `"did not treat gauge theory"` 这一字面句 = **NG（未逐字核）**，其命题内容由摘要支持 |
| B3 | 4D U(1) 的强耦合/偶极路属 Göpfert–Mack 1982（非 GK CMP 99） | **对（不在 GK CMP 99 内）**；`[C]`，Göpfert–Mack 原文本轮未核 | GK CMP 99 摘要对象为标量 φ⁴，无 U(1) 规范场 |

**总裁决**：**b03 的转述对** —— 走 **(a)**（逐字引原文佐证标题+摘要），**但附两条限定**：
- (a-1) b03 引号内的英文句子 `"did not treat gauge theory"` **不是**标题/摘要原文，正文本轮**未逐字核** ⇒ 该引号按 (c) 处理（不得当逐字引文用）。
- (a-2) b03 的"4D 定理"若被读作"质量隙定理"则**错**：GK 1985 结论是**红外渐近自由**（无质量模型），**不是质量隙**。b03 用词"无质量 φ⁴₄"本身正确，但不可反推 GK 提供任何隙。

**旁注（防止过度外推）**：GK 确实做过与 4D U(1) **对偶**的偶极气体/`(∇φ)⁴` 线（Annals of Physics **147** (1983) 198 一系，见检索命中 `Lattice dipole gas and (∇φ)⁴ models at long distances`）。所以"GK 从未碰过与规范场对偶的对象"是过强说法；但 **CMP 99 (1985) 197 这一篇**确为标量 φ⁴₄，不含规范场 —— b03 的判据范围（特指该文）成立。

## 5. 引用块（可直接抄用）

> K. Gawedzki, A. Kupiainen, "Massless lattice φ₄⁴ theory: Rigorous control of a renormalizable asymptotically free model", *Comm. Math. Phys.* **99** (1985) 197–252, DOI: 10.1007/BF01212281. ［题名+摘要已逐字核；正文未逐字核（付费墙/反爬）］

## 6. 残留洞

1. 正文（§1 引言）是否出现 `"we do not treat gauge theory"` 之类原句 —— **未核**；需机构订阅 PDF 或 ADS 扫描件。
2. Göpfert–Mack 1982（4D U(1) 强耦合/偶极路）原文 —— **本轮未核**（另开一物种）。

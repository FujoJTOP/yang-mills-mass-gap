# L3-TRANSMUTATION · `m = C·Λ` 形状先例的文献裁决

日期：2026-11-25 · 范围：只查"已**证明**形如 `m = C·Λ`（`C` 纯数、与格距无关）的质量隙"先例
方法：INSPIRE-HEP API / OpenAlex / 出版商摘要页逐条核 ID；**所有引文均取自检索到的摘要原文**，未凭记忆补写
禁止事项遵守：未执行 git commit；仅写本文件

---

## 0. 裁决摘要

| 判死线 | 结论 |
|---|---|
| **(a) 格点上严格证明 `m = CΛ`** | **不存在**（未查到）。最接近的是 Kupiainen 1980 的**格点** n-矢量 σ 模型定理：证明的是"**mass gap 存在**"+"1/n 展开"，**不是**逐字的 `m = CΛ` |
| **(b) 仅可积模型/连续场论的精确结果** | **成立且证据充足**。`m = C·Λ` 的形状在 2D 渐近自由可积模型中是被**算出**的精确等式（`C` 如 `8/e`、`√(32π/e)`、`sin(πΔ)/(πΔ)`），但都不是格点定理 |
| **(c) 连形状先例都没有** | **不适用**。形状先例存在（≈6 条独立模型），只是全部属于 (b) |

**一句话**：`m ∝ Λ` 作为**形状**在 2D 渐近自由理论里有硬先例（精确、`C` 为纯数）；但"**在格点上严格证明** `m = CΛ`（`C` 与 `a` 无关、`Λ` 为维度嬗变尺度）"这个更强的命题，**在文献里没有先例** —— 这正是 4D YM 那台"机器"缺的东西。

---

## 1. 逐条核过的条目（全部为 (b)：精确/可积，非格点定理）

### E1 — O(3)、O(4) 非线性 σ 模型：`m = (8/e)Λ_MS`、`m = √(32π/e)·Λ_MS`
- **作者/年**：P. Hasenfratz, M. Maggiore, F. Niedermayer, 1990
- **ID**：Phys.Lett.B **245** (1990) 522–528；**DOI 10.1016/0370-2693(90)90685-Y**；INSPIRE 296882（report BUTP-90/12-BERN）
- **模型**：d=2, O(3) / O(4) 非线性 σ 模型（渐近自由、维度嬗变）
- **定理/结论的逐字引用**（Elsevier 摘要，经 INSPIRE API 取回）：
  > "By comparing the Bethe Ansatz solution with perturbation theory we obtain the exact results m = 8 e Λ MS and m = ( 32 π e ) 1 2 Λ MS for the O(3) and O(4) non-linear σ-models, respectively. Here Λ MS is the renormalization group invariant scale in the MS scheme."
  （纯文本把分式压平了：即 `m = (8/e)·Λ_MS` 与 `m = (32π/e)^{1/2}·Λ_MS`）
- **`C`**：O(3) `C = 8/e = 2.9430…`；O(4) `C = √(32π/e) = 6.0814…`
- **格点严格？** **否**。可积模型：Bethe Ansatz 精确解 + 微扰论匹配定标度 ⇒ 精确，但不是格点定理
- **形状价值**：★ 最高。`C` 是**纯数**，`Λ_MS` 是 RG 不变量标度 ⇒ 逐字就是 `m = C·Λ` 的**精确等式**（不是拟合）

### E2 — O(N) σ 模型任意 N≥3 的精确质量隙
- **作者/年**：P. Hasenfratz, F. Niedermayer, 1990
- **ID**：Phys.Lett.B **245** (1990) 529–532；**DOI 10.1016/0370-2693(90)90686-Z**；INSPIRE 296884
- **模型**：d=2, O(N) σ 模型（N≥3）
- **逐字引用**：
  > "Starting from the exact S-matrix of the O(N) model we derive an integral equation, which determines the free energy as the function of the chemical potential coupled to a Noether charge. … this integral equation predicts for the mass gap m = ( 8 e ) 1 (N−2) × [ 1 Γ ( 1 + 1 (N − 2 ) ) ] Λ MS ."
- **`C`**：`C(N) = (8/e)^{1/(N−2)} / Γ(1 + 1/(N−2))`（N=3 退化为 `8/e`，与 E1 自洽）
- **格点严格？** 否（精确 S-矩阵 + TBA）

### E3 — Gross–Neveu 模型精确质量隙（TBA + 1/N 双重核验）
- **作者/年**：P. Forgács, F. Niedermayer, P. Weisz, 1991
- **ID**：Nucl.Phys.B **367** (1991) 123–143，**DOI 10.1016/0550-3213(91)90044-X**（part 1: TBA）；Nucl.Phys.B **367** (1991) 144–157，**DOI 10.1016/0550-3213(91)90045-Y**（part 2: 1/N）
- **模型**：d=2, O(N) 对称 Gross–Neveu（渐近自由、维度嬗变）
- **逐字引用（part 1）**：
  > "Using the postulated S -matrix of the O( N ) symmetric Gross-Neveu model the mass gap is calculated in terms of the Λ-parameter in the MS scheme."
- **逐字引用（part 2）**：
  > "The ratio of the mass gap to the Λ-parameter for the O( N ) symmetric Gross-Neveu model is calculated to order 1/ N . The result agrees with an exact formula derived previously by us. Furthermore the exact result is compared to estimates obtained from latticesimulations."
- **`C`**：**摘要未给出 `C` 的显式数值/表达式**（正文给出 `m/Λ_MS` 的 N 依赖闭式；未取到 PDF 正文，故**不在此处写死数值**，遵"不许凭记忆"）
- **格点严格？** 否。但 part 2 明确做了"精确 vs 格点数值"的对照 ⇒ 形状在格点上被**数值支持**，非定理
- **形状价值**：`m/Λ` 是"比"关系 ⇒ `m = C·Λ`；且 1/N 与 TBA 两条独立路径给同一 `C`

### E4 — 手征 SU(n)×SU(n)（主手征模型）精确质量隙
- **作者/年**：J. Balog, S. Naik, F. Niedermayer, P. Weisz, 1992
- **ID**：Phys.Rev.Lett. **69** (1992) 873–876；**DOI 10.1103/PhysRevLett.69.873**；INSPIRE 334874
- **模型**：d=2 手征 SU(n)×SU(n)（= 主手征模型）
- **逐字引用**：
  > "Using the Bethe ansatz equation the mass gap of the chiral SU(n)×SU(n) model is calculated in terms of the Λ parameter in the modified minimal subtraction scheme. The result is in reasonable agreement with recent Monte Carlo data."
- **`C`**：摘要未给数值；但"in terms of Λ"+"与 Monte Carlo 相符" ⇒ `m = C·Λ` 形状 + 格点数值同向
- **格点严格？** 否

### E5 — SUSY O(N) σ 模型：`m/Λ = 2^{2Δ} sin(πΔ)/(πΔ)`
- **作者/年**：J. M. Evans, T. J. Hollowood, 1994/95
- **ID**：Phys.Lett.B **343** (1995) 189–197；**DOI 10.1016/0370-2693(94)01477-T**；**arXiv hep-th/9409141**
- **逐字引用**：
  > "A formula for the mass-gap of the supersymmetric O($N$) sigma model ($N>4$) in two dimensions is derived: $m/\Lambda_{\overline{\rm MS}}=2^{2\Delta}\sin(\pi\Delta)/(\pi\Delta)$, where $\Delta=1/(N-2)$ and $m$ is the mass of the fundamental vector particle in the theory."
- **`C`**：`C = 2^{2Δ} sin(πΔ)/(πΔ)`，`Δ = 1/(N−2)`
- **格点严格？** 否（TBA + 微扰自由能双向匹配，兼解 CDD 歧义）

### E6 — SUSY `CP^{N−1}`：`m/Λ = sin(πΔ)/(πΔ)`  ← **任务点名的 `CP^{N-1}`**
- **作者/年**：J. M. Evans, T. J. Hollowood, 1994/95
- **ID**：Phys.Lett.B **343** (1995) 198–206；**DOI 10.1016/0370-2693(94)01478-U**；**arXiv hep-th/9409142**
- **模型**：d=2 **超对称 CP^{n−1}**（注意：是 SUSY 版本）
- **逐字引用**：
  > "A formula for the mass-gap of the supersymmetric $\CP^{n-1}$ sigma model ($n > 1$) in two dimensions is derived: $m/\Lambda_{\overline{\rm MS}}=\sin(\pi\Delta)/(\pi\Delta)$ where $\Delta=1/n$ and $m$ is the mass of the fundamental particle multiplet."
- **`C`**：`C = sin(π/n)/(π/n)`，`n = N`（CP 的 N）；`n→∞` 时 `C→1`（大 N 极限 `m = Λ`）
- **格点严格？** 否（Köberle–Kurak 精确 S-矩阵 + TBA）
- **重要限制**：**非 SUSY 的 2D CP^{N−1} 没有已知精确 S-矩阵**（Köberle–Kurak 是 SUSY 的）⇒ 本条的"精确性"不能直接搬到纯 CP^{N−1}；纯 CP^{N−1} 只有大 N 的 `m ∝ Λ`（此处未取到可逐字引用的大 N `C` 值，故不写数）

### E7 — SUSY O(N) σ 模型（大 N 重整化路线）
- **作者/年**：J. A. Gracey, 1992/93
- **ID**：Phys.Lett.B **298** (1993) 116–119；**DOI 10.1016/0370-2693(93)91716-Z**
- **逐字引用**：
  > "We write down the exact mass gap of the supersymmetric O( N ) σ model in terms of the Λ-parameter and a function of N − 2 by extracting the physical mass from the next to leading order renormalization of the model in the large N expansion in MS ."
- **`C`**：摘要未给显式（"function of N−2"）；与 E5 互为校验
- **格点严格？** 否

### E8 — 手征 Gross–Neveu 精确质量隙
- **作者/年**：P. Forgács, S. Naik, F. Niedermayer, 1992
- **ID**：Phys.Lett.B **283** (1992) 282–286；**DOI 10.1016/0370-2693(92)90020-5**
- **逐字引用**：
  > "Using the postulated S -matrix of the U( n ) symmetric chiral Gross-Neveu model the mass gap is calculated in terms of the Λ- parameter."
- **`C`**：摘要未给数值
- **格点严格？** 否

---

## 2. 唯一沾"格点"的严格结果：不是 `m = CΛ`，只是"`m > 0`"

### E9 — Kupiainen 1980：**格点上**的 `m > 0` 定理（不是形状定理）
- **作者/年**：A. J. Kupiainen, 1980
- **ID**：Commun.Math.Phys. **73** (1980) 273–294；**DOI 10.1007/BF01197703**；INSPIRE 无收录（math-phys 期刊），OpenAlex W2047723697；Project Euclid: `euclid.cmp/1103907876`
- **模型**：**格点**上的 n-分量非线性 σ 模型（经典 Heisenberg 模型），任意维
- **逐字引用**（Springer 摘要片段，经检索索引两次取回、内容一致；**出版页被 SSO/Anubis 拦截，未取到全文摘要**，故只引片段不补写）：
  > "The 1/n expansion is considered for the n-component non-linear σ-model (classical Heisenberg model) on a lattice of arbitrary dimensions. We show that the expansion for … Furthermore, the existence of a mass gap is established for these temperatures and n sufficiently large."
- **定理内容**：低温度 + `n` 充分大 ⇒ **1/n 展开可控 + mass gap 存在**
- **`C`**：**摘要未给**。1/n 展开的领头项给出 `m ∝ Λ`（形状），但**"`m = CΛ` 作为定理"不在摘要陈述里** ⇒ 不冒充
- **格点严格？** **是（格点 + 严格）**，但**只到"隙 > 0"这一步**
- **为什么不升级为 (a)**：定理的形式是"存在性 + 渐近展开"，不是"`m = C·Λ` 且 `C` 与 `a` 无关"的等式；`C` 的 `a`-无关性需要**跨标度的 RG 控制**，1980 年这条定理没做这一步

### E10 — 2D Yang–Mills：精确可解，但**没有 Λ** ⇒ 形状先例不成立
- **作者/年**：E. Witten, 1991
- **ID**：Commun.Math.Phys. **141** (1991) 153–209；**DOI 10.1007/BF02100009**；INSPIRE 325641
- **逐字引用**（摘要）：
  > "Two dimensional quantum Yang-Mills theory is studied from three points of view: (i) by standard physical methods; (ii) by relating it to the largek limit of three dimensional Chern-Simons theory and two dimensional conformal field theory; (iii) by relating its weak coupling limit to the theory of Reidemeister-Ray-Singer torsion."
- **关键判断（推理，非引文）**：2D YM 是**超可重整**理论，`g` 本身带质量量纲 ⇒ 理论里**唯一的标度就是 `g²`**，不存在维度嬗变产生的 `Λ`。因此 `m = C·g²`（`C` 纯数）即使成立，也**不是** `m = C·Λ` 的形状先例 —— 它是"隙 ∝ 唯一量纲参数"，与 4D 那种"无标度经典理论自己长出标度"是两种东西。
- **任务点名的 `SU(N) 2D YM gap 4π²/N`**：**未查到**可逐字引用的定理/公式支持。检索到的 2D YM 精确结果集中在配分函数、Wilson 圈、模空间体积、大 N 相变（Douglas–Kazakov / Gross–Witten 一侧）；**"gap = 4π²/N" 这个具体表达式我没有核到出处**，按纪律**不写进结论**。
- **格点严格？** 2D YM 的格点精确解（Migdal 递推）存在，但同样受"无 Λ"限制 ⇒ 不构成 (a)

---

## 3. 为什么 4D YM 里这条形状失效

1. **没有可积性**：E1–E8 的 `C` 全部来自"精确 S-矩阵 / Bethe Ansatz / TBA + 微扰匹配"。4D 纯 YM 无精确 S-矩阵、无 TBA ⇒ 拿不到任何 `C` 的闭式，连"算出 `C`"这条路都不存在。
2. **格点严格结果给的是"格点单位的隙"，不是 `Λ` 单位的隙**：强耦合团簇展开（Osterwalder–Seiler 型）给 `m ≈ a^{-1}·f(β)`，`f` 在强耦合区是 `O(1)` ⇒ 相关长度 `O(a)`，此时 `Λ` 尚未定义（`Λ` 是 `a→0` 极限里才有的 RG 不变量）。**"`m = CΛ` 且 `C` 与 `a` 无关"要求把 `a^{-1}f(β)` 一路搬运到 `Λ` 标度** —— 这一步在 2D 由渐近自由 + 1/n 或可积性兜住，在 4D 没有任何定理兜。
3. **本项目自己的 no-go（内部结论，非外部先例）**：`research/v3/species/A5-LIT-UNIFORM.md` 一线已裁："格点单位下与 `a` 一致的隙 ⇒ 连续两点函数被逼为 0"。即**"`a`-一致"在格点单位里说 ⇒ 真空；必须换成 `Λ` 单位说**。这与本文件 (a) 的空白是同一件事的两面：4D 缺的正是"从 `a^{-1}f(β)` 到 `C·Λ`"的搬运机器。
4. **`C` 的普适性问题**：E1–E8 里 `C` 依赖 scheme（`Λ_MS` vs `Λ_{\overline{MS}}`），但**与 `a` 无关**、且与 `N` 有闭式依赖。4D 若要 `m = CΛ`，必须先固定 scheme 再证 `C` 的 `a`-一致性 —— 目前连"`m/Λ` 有极限"都未证。

---

## 4. 检索记录（可复现）

查询串（实际发出）：
- `"mass gap" "dimensional transmutation" rigorous 2D lattice sigma model proof`
- `2D lattice Yang-Mills exact solution mass gap heat kernel SU(N) glueball mass`
- `exact mass gap O(3) O(4) nonlinear sigma model Hasenfratz Niedermayer`
- `Kupiainen 1/n expansion nonlinear sigma model mass gap rigorous`
- `"uniform in the lattice spacing" mass gap continuum limit two-dimensional`
- `two-dimensional lattice Yang-Mills "mass gap" rigorous spectral gap heat kernel SU(N) proof`
- `Fateev Zamolodchikov exact mass spectrum sigma model TBA mass gap`
- `Gross-Neveu exact mass gap "m =" Lambda "2 pi" dimensional transmutation large N`
- INSPIRE API：`t "exact mass gap"`（12 条命中，本文件采信 6 条）、`a Kupiainen and t "mass gap"`（0 条）、`t "quantum gauge theories in two dimensions" and a Witten`（1 条）
- OpenAlex：`doi:10.1007/BF01197703`（摘要字段为 null，故 Kupiainen 只引 Springer 索引片段）

**搜到但明确不采信**（记录以免下游误引）：
- Zenodo `records/17190165`（"Heat Kernel Sandwich Yields Sheet Tension and a 4D Yang Mills Mass Gap"）、`records/19401375`（"Yang-Mills Existence and Mass Gap — Hanners (2026)"）、`records/19204109`（"Constructive RG YangMills S3"，其 "Numerical Result 9.1 (Dimensional transmutation verified)"）、`records/23484057`（与本项目 A5 同题的自发布记录）、`records/19501991`（"Principia Initialis"，胶球质量预测）。
  理由：均为未审预印本/自发布，声称解决 Clay 问题但非同行评议定理；**不作为"形状先例"引用**。（`records/23484057` 与本项目 A5/A8 同题，若为本项目自有产出，应按内部结论引用而非外部文献。）
- `arxiv 2605.02373`（"Geometric QCD III"）、`inspire 296882` 之外的胶球谱文章：与"严格 `m = CΛ`"无关，未采信。

**未查到（诚实空缺）**：
1. 任何**格点严格**证明 `m = C·Λ`（`C` 与 `a` 无关）的模型 —— 无。
2. `SU(N) 2D YM` 的 `gap = 4π²/N` —— 无出处。
3. 非 SUSY `CP^{N−1}` 的精确 `C` 闭式（只有 SUSY 版 E6）—— 无逐字可引者。
4. GN 模型 `C(N)` 的显式闭式（摘要层面）—— 无（正文有，未取到 PDF）。

---

## 5. 对 `m = C·Λ` 这条线的净结论

- **形状有先例**：`m = C·Λ` 在 2D 渐近自由可积模型中是被**证明/算出**的精确等式，`C` 为纯数、与格距无关，且 `C` 有闭式（`8/e`、`√(32π/e)`、`(8/e)^{1/(N−2)}/Γ(1+1/(N−2))`、`sin(πΔ)/(πΔ)`、`2^{2Δ}sin(πΔ)/(πΔ)`）。
- **格点严格版没有先例**：唯一格点严格的是 Kupiainen 的"隙存在"（E9），不含形状等式。
- **因此本物种（L3-TRANSMUTATION）的合法定位**：把 `m = C·Λ` 当**已确立的形状目标**去证，而**不是**当"已有格点定理可直接搬运"。任何写作中若称"`m = CΛ` 在格点上已有先例"⇒ **越界，属诚实红线违规**。
- **4D 的确缺失步**：把 `a^{-1}f(β)`（强耦合隙）搬到 `C·Λ`（`C` 与 `a` 无关）的那台"跨标度搬运机"。2D 靠可积性/1-n 绕过，4D 无此捷径。

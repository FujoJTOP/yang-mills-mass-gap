# JJJ2-LUSCH — 面积律扇区的无隙方向：`ΔE=π/r`（Lüscher 项）是**无质量横向模的有限尺寸隙**，不是质量

> **靶（逐字）**：`FFF3-MPHYS` ①——"面积律扇区自带无隙方向（Lüscher）"。
> **三问**：(a) 哪个扇区、哪个方向；(b) 物理模式还是技术对象；(c) 是否约束隙的来源。
> **口径**：扇区 = 带一对静态夸克的扇区（`3⊗3*`）；`ΔE := E₁−E₀`；`d=4` ⇒ `d−2=2`；`γ := −(π/24)(d−2) = −π/12`。
> **限定**：Lüscher 的算术是**有效弦内**的结论；"四维 `SU(N)` 的真静态夸克扇区如此"是模型假设 + 数值确认，不是定理——0.3 单列此界。
> **纪律**：只交 `[T]`/`[NG]`；无百分数；每结论带「下一手」；零新 MC；既有文件零改动；未 `commit`；§6 位置图交 `lead`。

## §0 判词

| # | 断言 | 判 | 依据 / 逐字 | 下一手 |
|---|---|---|---|---|
| **0.1** | **(a) 扇区**=静态夸克（通量管）扇区；**方向**=`d−2=2` 个**横向位移** `h_i(z)`（世界面的法向涨落） | **`[T]`**（模型内） | `LW02`（`hep-lat/0207003`）§3.1：`h` 有 "`d−2` non-zero transverse components"；`S_eff=∫½∂_a h ∂_a h+…` | 把 `h_i` 与 `T2-TUBE-OVERLAP:261` 的 `log T` 修正对齐，判是否同一对象 |
| **0.2** | **逐字结论**：`E_n=σr+μ+(π/r){n−(d−2)/24}`；正文原句 **"the energy gap `ΔE` is equal to `π/r`"**；`w₁=d−2` | **`[T]`**（模型内，正文） | 同文 (3.6)(3.7)、§3.2；(3.10) `ΔE=(π/r)(1+b/r)`；(1.1) `γ=−(π/24)(d−2)` | 三条数（`π/r`、`w₁=2`、`γ=−π/12`）入台账，作可核对预言 |
| **0.3** | 把 0.2 当**四维 `SU(N)` 定理**用 | **`[NG]`** | `LW02` §1 逐字 "Assuming this to be the case"；`LW04` 摘要 "likely to be described … by an effective string theory" | 引用必写"模型内 + 数值"；不得写成"Lüscher 证明了 QCD 弦无隙" |
| **0.4** | **(b) 物理模式**：通量管横向涨落（集体坐标），质量项被**平移对称**禁戒 | **`[T]`**（模型内）+ 数值 | 作用量只含 `∂h`（含边界项 `(∂₁h)²`）⇒ `h→h+c` 对称 ⇒ 无 `m²h²`；`Meyer05` 摘要："The Luscher term only depends on the number of **massless** bosonic degrees of freedom" | 反向判据：紧致一个横向维 ⇒ 该模应**生 gap**（KT，`σL_p²≃1/8π`，中心荷减一）；不出现则 (b) 降级 |
| **0.5** | **(c) 硬约束（扇区版）**：`δ(r)=π/r`，`inf_r δ=0`；且 `E₀=V(r)=σr+μ+γ/r→∞` | **`[T]`**（模型内 + 逻辑） | (3.6)；(1.1) 的 `σr` 项 | 任何"弦谱有隙 ⇒ 全域有隙"的论证，用 `δ=π/r` 一句否证（墓碑 W05-④） |
| **0.6** | **(c) 越界版**："面积律与隙无关"、"隙不存在"、"必须换理论" | **`[NG]`** | `LW04` 摘要：`1/r` 展开；`1/r²` 修正必为零；正文："free-string energies are only the leading terms in an **asymptotic expansion**" ⇒ **不是界** | 只保留可证半边：`m>0 ⇒ σ>0`（Chatterjee'21，`B06-VAR:98`） |
| **0.7** | **对"往哪找"**：`δ` 只能在**真空扇区**（零绕数、局域算符）最轻态上；`σ` 降为**标度**，不作机制 | **`[T]`**（逻辑）/ **`[NG]`**（作桥） | 扇区不同 ⇒ 0.5 不触及真空扇区；`m_{0++}/√σ≈3.5–4` 为测量（`AREALAW-GAP:14,18,27`）；与 `MU5-SIGMA-M:74,84` 的 L4 同一缺口 | 落成可证对象：`λ₂(T_a)` 的 `(L,β)`-一致 `O(a)` 系数下界（= `FFF3 §3` 的靶） |

**一句话**：`ΔE=π/r` 不是"某个质量的近似"，而是**无质量二维场在长度 `r` 的盒中的第一 Dirichlet 模**（`k=π/r`）——该扇区的"隙"是**有限尺寸效应**，`r→∞` 必归零。把它当 `δ`，是把有限尺寸隙误认成质量隙。**判 `[T]`（模型内）/ `[NG]`（作四维定理、作 `σ⇒m` 的桥）。**

## §1 文献

| # | 文献 | ID | 读到 | 逐字 / 用途 |
|---|---|---|---|---|
| L1 | Lüscher–Symanzik–Weisz, *Anomalies of the free loop wave equation in the WKB approximation*, NPB **173** (1980) 365–396 | DOI `10.1016/0550-3213(80)90009-7` | 【读摘要/书目】DESY 80-074（KEK 扫描 `8005/8005074.pdf`）首页；**系数未逐字核正文** | `γ` 的出处；本文件的 `γ/r` 逐字取自 L2 (1.1)，归因 L1 |
| L2 | Lüscher–Weisz, *Quark confinement and the bosonic string*, JHEP **0207** (2002) 049 | arXiv `hep-lat/0207003`；DOI `10.1088/1126-6708/2002/07/049` | **【读正文】** §1、§3.1–3.3（ar5iv 全文） | (1.1) `γ`；(3.1)(3.2) `h`、`d−2`；(3.6) `E_n`；**"the energy gap ΔE is equal to π/r"**；(3.7) `w₁=d−2`；(3.10) `ΔE=(π/r)(1+b/r)` |
| L3 | Lüscher–Weisz, *String excitation energies in SU(N) gauge theories beyond the free-string approximation*, JHEP **0407** (2004) 014 | arXiv `hep-th/0406205`；DOI `10.1088/1126-6708/2004/07/014` | **【读摘要】**+【读正文片段】 | 摘要："expansion of the excitation energies in powers of `1/r`"；正文："free-string energies are only the leading terms in an **asymptotic expansion** in powers of `1/r`" |
| L4 | H. B. Meyer, *Vortices on the worldsheet of the QCD string*, NPB **724** (2005) 432–450 | arXiv `hep-th/0506034`；DOI `10.1016/j.nuclphysb.2005.07.001` | **【读摘要】** | "The Luscher term only depends on the number of **massless** bosonic degrees of freedom"；紧致横向维 ⇒ 位移矢量的**周期分量** "develops a mass gap"、中心荷减一 |

**逐字核对（"gapless"/"massless" 各指什么）**：`LW02` 全文**不用 "gapless"**——它的原句是 gap `= π/r`（`r→∞` 归零）；"massless" 的逐字出处是 `Meyer05` 摘要，指**世界面上的玻色自由度**（`d−2` 个横向位移），不是指四维粒子。故本文件把"无隙方向"读作：**该扇区的隙 = 无质量模的有限尺寸 `π/r`**。

**勘误（本轮）**：`hep-lat/0409003` 是 Chroma（Edwards–Joo），**不是** L3；L3 的 ID 是 `hep-th/0406205`（已核摘要页）。L1 的 DOI 初猜 `…90209-6` 错，实为 `…90009-7`（两处索引一致）。

## §2 物理模式，还是技术对象

**是物理模式（三条）。**
1. **它是通量管的形状自由度**：世界面偏离极小面的法向位移；`LW02` §3.1 逐字 "The important degree of freedom … is the deviation of the fluctuating surface from the minimal surface in the directions orthogonal to the latter"。不是规范冗余，也不是拟合参数。
2. **无质量由对称性强制**：`S_eff` 只含 `∂h`（边界项 (3.8) 也只有 `(∂₁h)²`）⇒ `h→h+c` 是精确对称 ⇒ 不许出现 `m²h²`。这是横向平移破缺的 Goldstone 型。`Meyer05` 把它数成 "massless bosonic degrees of freedom"，而 `γ` 只依赖这个**数目**——`γ` 的普适性就是"数无质量模"。
3. **有可测效应**：`LW02` 摘要逐字 "we confirm the presence of a `γ/r` correction … with a coefficient `γ` as predicted by the bosonic string theory"。

**两条边界（不许越）。**
- (i) 它是**有效理论**的集体模，不是四维基本场；"无隙"是**该扇区**（外源/长弦）的激发间距，不是真空扇区的谱下沿。
- (ii) `π/r` 是**无质量场在长度 `r` 的盒中的第一模**（Dirichlet）⇒ 是**有限尺寸隙**，`r` 有限时它 `>0`，但**不是质量**；`r→∞` 归零。

**判**：**物理模式，不是技术对象；但它是"扇区内的无质量模 + 有限尺寸隙"，不是"四维谱里的质量"。**

## §3 约束：隙往哪找

**(c) 裁决，分两读。**
1. **扇区版 `[T]`**：`inf_r(E₁−E₀)=inf_r π/r=0`，且 `E₀=σr+…→∞` ⇒ 面积律扇区**自身既无质量隙、也不提供最轻态**。所以"弦扇区有隙"**不是**可用的前提。
2. **越界版 `[NG]`**：不能由此推出"面积律与隙无关"、"隙不存在"、"必须换理论"。`1/r` 展开是**渐近级数**（L3），**不是界**；`σ` 仍是标度（`m_{0++}/√σ` 的测量）。**同一理论的 UV 端**才是 `δ` 的所在。

**对"往哪找"的三条硬含义。**
- (i) **对象**：`δ` 只能落在**真空扇区**（零绕数、局域规范不变算符）的最轻态；弦激发是 `π/r`，不是 `δ`。
- (ii) **方法**：任何把 `1/r` 展开截断当"下界"的论证无效——它连符号都不稳（`b` 项）。
- (iii) **地位**：`σ` 从"桥"降为"刻度"；已证箭头只有 `m>0 ⇒ σ>0`，反向仍无定理，且现在知道反向**不能**靠弦谱补。

**下一手**：① `λ₂(T_a)` 的 `(L,β)`-一致 `O(a)` 系数下界立项（`FFF3 §3`）；② 用 `Meyer05` 的 KT 开关做"横向模是物理的"独立检验（紧致横向维 ⇒ 应生 gap）；③ 墓碑一条：凡"弦谱有隙 ⇒ 全域有隙"，引 (3.6) 一句否证。

## §4 检索四要素（`AGENTS §1`）与限定

| 工具 | 状态 | 检索式 / 目标（逐条） | 命中 | 失败标记 |
|---|---|---|---|---|
| `grep` | OK | `FFF3-MPHYS\|JJJ2\|Lüscher\|Luscher` @ `research/v3` `*.md` | **148** | — |
| `grep` | OK | `Lüscher\|Luscher\|gapless\|无隙\|1/R` @ `MU5-SIGMA-M.md` | **10** | — |
| `glob` | OK | `research/v3/species/*.md` | **957**（前 100 落盘） | — |
| `web_fetch` | OK | `arxiv.org/abs/hep-lat/0207003`；`ar5iv.labs.arxiv.org/html/hep-lat/0207003` | 全文（§4.1 后截断） | — |
| `web_fetch` | OK | `arxiv.org/abs/hep-th/0406205`；`arxiv.org/abs/hep-th/0506034` | 摘要 + 书目 | — |
| `web_fetch` | **TOOL-FAIL** | `doi.org/10.1016/0550-3213(80)90209-6`（初猜） | **0** | cross-origin redirect 到 linkinghub；改由 `web_search` 取 `…90009-7`（两处索引一致）⇒ **初猜 DOI 错，已更正** |
| `web_fetch` | 无效目标 | `arxiv.org/abs/hep-lat/0409003` | 命中 Chroma | 非 TOOL-FAIL：ID 误猜；L3 实为 `hep-th/0406205` |
| `web_search` | OK | ① `Lüscher Weisz "String excitation energies…" arXiv` ② `"effective string" transverse fluctuations massless "Lüscher term" gapless` ③ `"Anomalies of the free loop wave equation…" DOI` | 10 / 10 / 9 源 | — |

**限定**：①`[T]` 的窗口 = `LW02` 的自由有效弦（`T,r ≫ 1/√σ`，`d≥3`）；②四维 `SU(N)` 上"扇区隙 `π/r`"只有数值确认，无定理（0.3）；③边界耦合 `b` 与 `1/r³` 项未用；④本文件不重算任何数，`π/12`、`π/r` 为逐字转录。

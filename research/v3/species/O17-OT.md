# O17-OT — "方差 ≤ 曲率倒数"型不等式在最优传输里的已发表形式：四条候选逐个判，以及"随机曲率 BL"是否存在

**日期**: 2026-11-25 · **状态**: 已判（四条候选：1 条换名、1 条不适用、1 条不适用+精确反例、1 条**未查到 → 写成猜想**）· **禁 git commit**（遵守，未 commit）· **只写本文件**

**先读**：`research/v3/species/GPP-STRUCT.md`（191 行，读毕：`G''(0)=E[H_vv+Φ''|0]−Var(G_v+Φ'|0)`，`[T]`；四候选 Chebyshev/FKG/Griffiths/Ginibre 判决为 `[NG]`）。
**O08：不存在。** `Test-Path research/v3/species/O08*` = `False`；`glob **/O08*` = 0；`Get-ChildItem research/v3 -Recurse -Filter O0*` = 0 命中。⇒ 无前置约束，本文件只对齐 `GPP-STRUCT.md`。

---

## 0. 判死线（先读这段）

| 判 | 内容 |
|---|---|
| **(a) 有"随机曲率"的 BL 吗** | **有，但只在 PSD-a.s. 下**：**Hargé 不等式**（JFA 254 (2008) 267–300, DOI `10.1016/j.jfa.2007.07.019`）的曲率是**状态依赖的随机矩阵** `Hess V(X)^{-1}`——但要求 `Hess V(X) ⪰ 0` 几乎处处（凸）。⇒ 形式对上了，**约束没松开**：Savvidy 负模（`λ_min(Hess S_W)=−g\|B\|<0`）仍然被排除。**不是活口。** |
| **(b) 无 PSD 假设的"随机曲率 BL"** | **未查到**（检索式与命中见 §4；命中 0）。⇒ 按纪律写成**可陈述猜想 RC-BL**（§5，精确形式：条件 BL + 方差分解的 between 项）。**这是本文件的正式交付物。** |
| **(c) 中间活口** | 三个：(c1) **Hargé 只需 `Hess V ⪰ 0`，不需 `⪰ λI`**（比 BL 弱、比"强凸"弱）；(c2) **Kolesnikov 2026 的 rank-one Hessian 修正**（arXiv:2610.21669，even log-concave + 局部 log-BM）；(c3) **OT 证明只用到"沿测地线方向"的凸性**（单方向凸）。**三条都不是我们缺的那条**，但 (c1)(c2) 的形状可直接拿来做猜想的模板。 |
| **对 H1 的净收益** | **零**（在 1D 投影上 Hargé **恰好等价** H1 = 换名 `[NG]`，与 `GPP-STRUCT §4(c)` 一致）；n 维（纤维）版本被墙挡住。 |

---

## 1. 候选 1：Hargé 不等式

**陈述（标准形式，`[C]`：公式级未逐字核，见下"证据"）**

设 `dμ = e^{−V}dX`，`V ∈ C²`，`Hess V(X) ⪰ 0`（凸，允许退化、用伪逆）。则对一切 `f ∈ C_c^∞`，

```
Var_μ(f)  ≤  E_μ[ ⟨ Hess V(X)^{−1} ∇f , ∇f ⟩ ] .          (Hargé)
```

**与 BL 的关系 `[C]`**：BL 是 `Hess V ⪰ λI, λ>0` 时的 `Var ≤ (1/λ)E|∇f|²`。因 `Hess V ⪰ λI ⇒ Hess V^{−1} ⪯ (1/λ)I`，故 **Hargé ⟹ BL**，且是**严格加强**：曲率不再是一个常数，而是**逐点的随机矩阵**；Gaussian 情形取等（`Hess V^{-1} = Σ`，正是精确协方差）。所以 Hargé = **"方差 ≤ 随机曲率倒数"的已发表原型**，但它把"随机"理解为"随 `X` 随机"，并**强制该随机矩阵 PSD**。

**对 H1 是什么关系（本文件的关键判定）`[T]`（我自己的代数论证）**

H1 要的是 1D 条件边缘上的
`Var(∂_vV | ξ=0) ≤ E[∂_v²V | ξ=0]`。
在 1D 上 `Hess V = g''`、`Hess V^{-1} = 1/g''`，Hargé 就是经典 1D BL `Var_ρ(φ) ≤ E_ρ[φ'²/g'']`（此形式在 arXiv:2610.21669 正文片段里逐字出现，见 §3）：
`Var_ρ(φ) ≤ ∫_0^∞ (∂_r∂_εℓ)² /(−∂_r²ℓ)`。而 1D 情形 `E[φ'²/g''] ` 有限 ⟺ `g''>0` 于 supp，**且取 `φ = id` 即得 `g'' ≥ 0` ⟺ `G''≥0`**。
⇒ **1D 投影上 Hargé 与 H1 互相蕴含 ⇒ 换名，`[NG]`**（与 `GPP-STRUCT §4(c)` 第 1 条"作为路线不产生新估计"逐字同判）。
⇒ 唯一可能有内容的是**纤维（`n−1` 维）版本**：`Hess V|_{T\{ξ=0\}} ⪰ 0`。**这正是 Savvidy 墙**：`λ_min(Hess S_W) = −g|B| < 0`，且负模是**横向磁**的，**落在纤维内部**（不是在 `v` 方向上）。⇒ **Hargé 的假设在 `β=2.5` 上假，候选 1 不适用。**

**证据级别**：`web_fetch` 到 `zbmath.org` = **HTTP 403**、`sciencedirect.com/…/S0022123607003126` = **HTTP 403**（付费墙）⇒ **未读到正文、未读到摘要原文**。以下为**书目级已核**：标题 "Reinforcement of an inequality due to Brascamp and Lieb"、`J. Funct. Anal. 254(2) 267–300 (2008)`、DOI `10.1016/j.jfa.2007.07.019`（OpenAlex `[T-书目]`，HTTP 200）。摘要片段（NSTL 书目页）逐字："we obtain a reinforcement of an inequality due to Brascamp and Lieb and a reinforcement of Poincaré's inequality for general logarithmically concave measures" ⇒ `[C]`（片段级，**非全文**）。公式 (Hargé) 本身标 `[C]/[T-standard]`，**未逐字核**。

---

## 2. 候选 2：Cordero-Erausquin / 最优传输证明的 BL —— 需要凸性吗？

**答：需要，而且是"沿运输方向"的凸性；把它升级成对一切 `f` 的方差不等式，就等于要求全局凸性。`[C]`**

- 路线：`CMS 2001`（Cordero-Erausquin–McCann–Schmuckenschläger, *A Riemannian interpolation inequality à la Borell, Brascamp and Lieb*, **Invent. Math. 146 (2001) 219–257**，`[T-书目]`）证明的是 **RBBL（Borell–Brascamp–Lieb 的黎曼插值形式）**：沿测地线 `x_t`，`ρ_t(x_t)^{−1/n} ≥ (1−t)ρ_0^{−1/n} + t ρ_1^{−1/n}`，前提是 `Ric ≥ 0` 与**测度沿测地线的 log-凹性**（即 `d²/dt² V(x_t) ≥ 0`）。
- **BL 方差不等式**从 RBBL 取出的方式是：对每个检验函数 `φ`，把 1D 边缘 `r ↦ (P_tφ)(x_t^{r})` 上的 1D BL 沿每条测地线用一次。**所以凸性只被用在"运输方向"上**——这是 (c3) 活口；但因为 `φ` 任意、对每个方向都要用一次，**合起来等价于全局凸性**。
- 同族的 HWI 路线（`Otto–Villani 2000`；CMS 2001 §6）：`Ent(ν|μ) ≤ W₂√I − (κ/2)W₂²`，要由此得 BL 必须 `κ > 0`（二次项要"帮"而不是"害"）。⇒ **OT 证明没有绕开凸性**，它只是把凸性搬到了"沿测地线的二阶导"上。
- 已核的旁证 `[C]`（摘要级，逐字）：**"Dimensional improvements of the logarithmic Sobolev, Talagrand and Brascamp–Lieb inequalities"**, *Ann. Probab. 46 (2018)*, DOI `10.1214/17-aop1184` — 摘要明写 "we use **optimal transport methods** and the **Borell–Brascamp–Lieb** inequality"，且产出是"**deficit** in the classical inequalities"（即 BL 的定量加强，**不是**去掉凸性）。`[C]` 摘要已读。

**判：候选 2 = 不适用（不松开凸性），但它给出了正确的"局部化"模板**：把 BL 局部化到每条测地线的 1D 问题上 ⇒ 猜想 RC-BL 应当写成"**逐测地线 / 逐条件块**"的形式（§5）。

---

## 3. 候选 3：`W₂` 收缩 + 方差；非凸 `κ<0` 版本存在吗？

**答：`κ<0` 的"`Var ≤ (1/κ)E|∇f|²`"根本不可能成立——有初等精确反例（`[T]`，本文件自证）；`CD(κ,∞)` 只在 `κ>0` 给方差上界。**

**精确反例 `[T]`（1D，手算，非引文）**：取 `V(x) = −x²/2 + x⁴/4`，则 `V''(x) = −1 + 3x² ≥ −1`，即 `Hess V ⪰ κI` 且 `κ = −1`（且是**最优**常数）。取 `f(x)=x`：`E|∇f|² = 1`。若"`κ<0` 版 BL"成立则须 `Var_μ(x) ≤ (1/κ)·1 = −1`，而 `Var_μ(x) > 0`。**矛盾。** ⇒ **任何把 `1/κ` 直接搬到 `κ<0` 的写法都是假的**，不是"没查到"，是"假"。`[T]`
（对照：1D 的正确形式是 `Var_ρ(φ) ≤ E_ρ[φ'²/g'']`，分母是**逐点** `g''(x)>0`，不是 `inf g''`。这正是 Hargé/1D-BL 的形状，也是候选 3 唯一能救回来的写法。）

**`CD(κ,∞)` ⇒ 方差**：`W₂` 沿半群的收缩 `W₂(μ_t,ν_t) ≤ e^{−κt}W₂(μ,ν)` 对 `κ ∈ R` 都成立（`κ<0` 是**扩张**），但它**只**给"输运距离随时间"的信息；要转成 `Var` 上界得走 `Var ≤ W₂·√I` 型链（HWI），而那一步要 `κ>0`。⇒ **`CD(κ,∞)` 与"方差 ≤ 曲率倒数"之间的桥，在 `κ<0` 时断在 HWI 的二次项上**，`[C]/[T-standard]`（Bakry–Gentil–Ledoux 教科书级；本文件未逐字核该书页码）。

**判：候选 3 = 不适用**（`κ<0` 版本**假**，不是缺）。**但它精确说明了缺口的形状**：缺口不是"曲率符号"，是"**曲率的符号可以逐点变号，只要它在平均/条件意义下被压住**"。这就是猜想 RC-BL 要写的东西。

---

## 4. 候选 4：`conditional Brascamp–Lieb` / `random curvature Brascamp–Lieb` / `stochastic BL` —— 按 `AGENTS.md §1` 报检索式与命中数

**工具状态（§1 四项）**：`web_search` = **OK**（6 次调用，均返回有效外部源，无空响应）；`web_fetch` = **OK**（arXiv/OpenAlex/NeurIPS 均 200；`zbmath.org` 与 `sciencedirect.com` = **403 付费墙**，已标）；`api.semanticscholar.org` = **HTTP 429**（限流，2 次，`TOOL-FAIL` 不计入命中）；`pwsh`/`glob`/`read` = **OK**；`grep` = 未调用（用 `pwsh Get-ChildItem` 替代）。**无 `TOOL-FAIL` 计入判决。**

| # | 检索式（逐字） | 工具 | 命中 | 结论 |
|---|---|---|---|---|
| 1 | `"conditional Brascamp-Lieb inequality"`（**带引号精确短语**） | web_search OK | **10 源，精确短语 0** | 全部是**非线性/多线性 BL**（Bennett–Carbery–Christ–Tao 系）、**forward-reverse BL**（arXiv:2102.06935）、**regularized BL**（arXiv:2110.02841）；"conditional"只作普通词出现（"conditional distributions"、"conditional theorem for local→global"）。⇒ **"conditional BL" 这个术语在已发表文献里不存在** |
| 2 | `"random Brascamp-Lieb" curvature quenched` | web_search OK | **10 源，精确短语 "random Brascamp–Lieb" 0** | 唯一"随机"命中是 **SciRate 摘要片段**：Mihăilescu（math-ph）"**quenched** in the random environment … **with high probability the environment is sufficiently regular to apply a convexity argument and the Brascamp–Lieb inequality**" ⇒ 这是**"在高概率集合上先用凸性、再套 BL"**，不是随机曲率 BL 定理 |
| 3 | `stochastic Brascamp-Lieb inequality quenched curvature spin glass`（与 2 同批） | web_search OK | 批内合并 10 源，逐式不分解 | 无 "stochastic BL" 术语命中 |
| 4 | `Brascamp-Lieb variance inequality assumed convex in one direction only semiconvex conditional Hessian` | web_search OK | **10 源，单方向凸 0** | 命中全部要求全局 `Hess V ⪰ 0`（含"compact convex sets 版 BL"、「`Ṽ−ρx²/2` **不是全直线凸、只在区间 `(0,∞)` 上凸**、故在该区间上用 BL」的 thesis 片段 ⇒ 这是**分块凸**，不是单方向凸） |
| 5 | `Brascamp-Lieb CD(kappa,infinity) variance bound negative kappa non-convex`（批内） | web_search OK | 合并 10 源 | 无 `κ<0` 方差上界命中（与 §3 的 `[T]` 反例一致：**该东西不存在**） |
| 6 | OpenAlex full-text：`Reinforcement of an inequality due to Brascamp and Lieb` | web_fetch OK | `count = 56` | 精确定位 Hargé（DOI 已核）；附带命中 **"Intertwinings and generalized Brascamp–Lieb inequalities"**, *Rev. Mat. Iberoam.* 2018, DOI `10.4171/rmi/1014`，摘要逐字："entail **new functional inequalities of Brascamp–Lieb type** for log-concave distributions **and beyond**" ⇒ **"beyond log-concave" 是已发表的**（(c) 类候选；作者未核） |
| 7 | OpenAlex full-text：`Remarks on decay of correlations and Witten Laplacians Brascamp-Lieb inequalities and semiclassical limit` | web_fetch OK | `count = 52` | 定位 **Helffer–Sjöstrand, JFA 155 (1998)**, DOI `10.1006/jfan.1997.3239`；同族 **"On Witten Laplacians and Brascamp–Lieb's Inequality on Manifolds with Boundary"**, *Integral Equations Operator Theory* (2017), DOI `10.1007/s00020-017-2354-1`；以及 Witten 拉普拉斯路线的"**非凸例子**"片段（"This result contains **non convex examples** when `v_j` is negative"——**SD 正文片段，付费墙，未逐字核**）⇒ **Witten-Laplacian 路线的 BL 确实能出非凸例子**，但它给的是**谱隙/关联衰减**，不是 `Var ≤ 1/κ`。 |
| 8 | 自存 preprint 命中（须剔除） | web_search OK | ≥2 源 | `zenodo.org/records/18860202` "**Curvature Brascamp-Lieb domination**" 是**本项目自己的 preprint**（Yang–Mills, Exceptional Groups），**不是已发表文献**，不得引用为出处 |

**§4 判**：**(a) 的答案 = "没有无 PSD 假设的随机曲率 BL"**（检索式 1/2/3/4 全部 0 命中；工具全部 OK）。唯一形状相近的已发表物是 **Hargé（随机曲率但 PSD）** 与 **Kolesnikov 2026（rank-one 修正 Hessian，但 even log-concave）**。⇒ **落到 (b)：写成猜想。**

### 4.1 两个"活口形状"的逐字引证（供猜想模板用）

**(c2) Kolesnikov, arXiv:2610.21669（2026-10-15 提交；`[C]` 摘要已读，逐字）**：
> "We prove a **strengthened Brascamp–Lieb inequality for even log-concave measures** whose positive sublevel bodies satisfy the **local logarithmic Brunn–Minkowski inequality**. The conclusion can be written either by **adding a rank-one term to the Hessian of the potential** or by adding an explicit nonnegative remainder to the variance."

**为什么这条值钱（`[C]`，我的读法）**：(i) **evenness**——格点 Wilson 测度**恰好是 even**（`S_W(A)=S_W(−A)`：plaquette 取逆后共轭、`Haar` 反演不变），**这一条白送**；(ii) **rank-one 修正**——Savvidy 不稳定模是**低秩**（横向磁，单一模对）⇒ 与"给 Hessian 加一个 rank-one 正项"的形状**对得上**；(iii) 但假设仍是 **log-concave**（我们的不是）。⇒ **不能直接用，但它是"rank-one 修补"这一物种的已发表证据**，把我们的猜想从"凭空的"降成"有先例的"。

**(c1) Hargé 的弱形式**：只需 `Hess V ⪰ 0`（退化允许、用伪逆），不需要 `⪰ λI`。⇒ 比 BL 弱、比"强凸"弱，是**半凸（PI/半凸）型**中间活口。**但对 H1 净收益为零**（§1 的 1D 换名论证）。

---

## 5. 猜想 RC-BL（随机曲率 Brascamp–Lieb）—— 精确陈述 `[M]`

> **标签 `[M]`（猜想 / 模型陈述，**未证**，不得标 `[T]`）。**

**设定**：`(E,⟨·,·⟩)` 为 `n` 维内积空间，`dμ = e^{−V}dX`，`V ∈ C²`，`H := Hess V`（**不假设 PSD**）。写 PSD 分解 `H = H_+ − H_-`，`H_± ⪰ 0`，`|H| = H_+ + H_-`（逐点，`H_+ = (|H|+H)/2`）。设 `𝒴` 是把 `E` 分成直和 `E = ⊕_j E_j` 的一个**可测分块**（"方向分块"，例如按 Savvidy 模态分解：负模一块、其余一块），`Y` 为对应的块坐标投影，`κ_j(y) > 0` 为**块上的随机曲率**（可测、可依赖 `y`）。写 `P_j` 为到 `E_j` 的正交投影。

**假设 (RC)**（随机曲率的"平均正性"，替代 PSD-a.s.）：
```
(RC1)  对每个块 j 与 a.e. y：  P_j H(x) P_j|_{E_j} ⪰ κ_j(y) P_j      （仅在"好集" G_j(y) 上要求）
(RC2)  E_μ[ 1_{G^c} · tr( P_j |H_-| P_j ) ]  ≤  θ_j · E_μ[ κ_j(Y) ] ,   θ_j ∈ [0,1)   （坏集的负曲率被"预算"压住）
(RC3)  Σ_j θ_j < 1 .
```

**结论（RC-BL）**：
```
Var_μ(f)  ≤  Σ_j E_μ[ (1/κ_j(Y)) · E( |P_j ∇f|² | Y ) ]  +  Var_μ( E[f | Y] )  ,     ∀ f ∈ C_c^∞ .
```

**正确的读法（这是本猜想的内容，也是它不是平庸的原因）**：
1. 右边**第一项**是"**随机曲率倒数 × 块内梯度能量**"——`κ_j(Y)` 是**随机量**，这正是任务要找的对象；它由**块内**的条件 BL 给出（**块内** `Hess V|_{E_j} ⪰ κ_j > 0`，可以逐块成立，即使整体不定）。
2. 右边**第二项** `Var(E[f|Y])` 是**必须付的代价**：它是"**块间**方差"。`f` 只依赖 `Y` 时第一项为 0、第二项 = `Var(f)`，故**等号可取** ⇒ **不等式尖锐（不可去掉第二项）**。这就是"随机曲率 BL"与"常数曲率 BL"的本质差别：**你换来的是逐块的随机曲率，付出的是块间方差。**
3. **(RC1) 只要求逐块、且只在好集上**（`G_j(y)` = "该块上负模不激发"的事件，物理上就是**小场事件**）；**(RC2)(RC3)** 是"**曲率预算**"：坏集上的负曲率总投入必须小于好集上的正曲率产出。这正是把 `GPP-STRUCT §4(a)` 的"两半都是动力学的"翻译成可陈述的不等式条件。
4. **特别的 1D 情形**（`n=1`，块 = 整个 `R`）：`Y` 平凡 ⇒ 第二项为 0，`(RC1)` 退化为 `g'' ≥ κ > 0` a.e. ⇒ **退回经典 1D BL**。⇒ **猜想在凸情形自动正确，不 vacuous**（过"自由场酸测"的第一关：自由场 `V=Φ` 凸 ⇒ 不等式成立且**不空**；`f = id` 给 `Var ≤ E[1/κ]`，是一个真界）。
5. **与 H1 的关系 `[C]`**：取 `Y = ⟨v,·⟩` 的块分解（`v` 一块、纤维一块），`f = ∂_vV`，则 `(RC)` 的结论**恰好是** H1 所需的 `Var(∂_vV|0) ≤ E[∂_v²V|0]` **加上**一个 `Var(E[∂_vV|Y])` 修正项。⇒ **RC-BL 若成立，H1 的 `δ=0` 成员随之成立，代价是块间项**；块间项是否可控是**第二个独立问题**（本文件不定它）。

**为什么不 vacuous / 会在哪失败（诚实标注）**：`(RC2)` 里的 `E_μ[·]` 是**无条件的**，而 H1 要的是**条件化在 `ξ=0`** 的版本 ⇒ **两者之间还差一次"条件化修正"**（这正是 `GPP-STRUCT §4(b)` 记的"条件化修正是缺口的原因"）。把 `(RC2)` 的条件化版本写出来，就是把 RC-BL 落到 H1 上的**精确缺失的那一步**。`[M]`

**辅助猜想（更弱、更可能先证）RC-BL′**：若 `E_μ[H] ⪰ 0`（**平均**曲率 PSD，逐点不定）且 `E_μ[tr H_-] ≤ θ E_μ[tr H_+]`，则
`Var_μ(f) ≤ (1/(1−θ)) · E_μ[ ⟨ (E_μ[H])^{-1} ∇f, ∇f ⟩ ] · (1 + osc)` —— **不给精确形式**（`osc` 未定），仅记录：**"平均 PSD"是比"逐点 PSD"更弱、且很可能够用的中间假设**。`[M]`
> 依据（`[C]` 片段）：一份 INSA-Lyon 博士论文（HAL `tel-05798816`）正文写 "The **Cramér–Rao bound is a sort of reverse Brascamp–Lieb inequality** … We now study how to use the Brascamp–Lieb inequality and the Cramér–Rao bound to **transfer weak information on the Hessian of V to bounds on the covariance** of the measure `μ = e^{−V}`" ⇒ **"用 Hessian 的弱信息换协方差界"这个思路已有人在用**，但产出是**界/算法**，不是"随机曲率定理"。`[C]`

---

## 6. 四条候选的最终判词表

| # | 候选 | 给什么 | 需要什么假设 | 对 H1 | 判 |
|---|---|---|---|---|---|
| 1 | **Hargé** `Var ≤ E⟨Hess V^{−1}∇f,∇f⟩` | **随机（状态依赖）曲率**的方差上界，严格强于 BL | `Hess V ⪰ 0`（**PSD a.s.**，可退化） | **1D 投影上恰好 ⟺ H1 ⇒ 换名 `[NG]`**；纤维版被 Savvidy 负模挡住 | **不适用（换名）** |
| 2 | **CMS/OT 证明的 BL** | 正确的**局部化**模板（逐测地线 1D BL） | 沿测地线 `d²/dt²V ≥ 0`；`Ric ≥ 0`；升级到一般 `f` ⇒ **全局凸性** | 不给 | **不适用**（但给出猜想的写法） |
| 3 | **`W₂` 收缩 + 方差 / `CD(κ,∞)`** | `κ>0` 才有 `Var` 上界；`κ<0` **精确反例**（`V=−x²/2+x⁴/4`，`Var ≤ −1` 荒谬）`[T]` | `κ>0` | **`κ<0` 版本是假的，不是缺的** | **不适用** |
| 4 | **conditional / random curvature BL** | — | — | — | **未查到（§4：4 条检索式全 0 命中，工具 OK）⇒ 猜想 RC-BL（§5）** |
| 4′ | **Kolesnikov arXiv:2610.21669** | **rank-one Hessian 修正**的加强 BL；even + 局部 log-BM | even **log-concave** | 不能直接用；**evenness 白送、rank-one 形状对得上** | **中间活口 (c2)** |
| 4″ | **"Intertwinings and generalized BL"** Rev. Mat. Iberoam. 2018, DOI `10.4171/rmi/1014` | BL 型不等式 "for log-concave distributions **and beyond**" | 未核（付费墙） | 待查 | **中间活口 (c3)，作者/假设未核** |

---

## 7. 未闭合（诚实标注）

1. **Hargé 的公式未逐字核**（SD/zbmath 双 403）。本文件的 `(Hargé)` 是 `[C]/[T-standard]`；若要用进论文，必须先拿到正文核一次（机构 VPN 或作者主页 PDF）。
2. **"随机曲率 BL 不存在"是否定性检索结论**（0 命中 × 4 式），**不是不可能性定理**。要升成 `[NG]` 还需一条"任何方差上界族都必须过 PSD 曲率"的分类定理——本文件没有。
3. **RC-BL 的 `(RC2)` 是无条件期望，H1 要条件期望**；两者之间的"条件化修正"是本文件**明确留空的那一步**（§5 末）。
4. `CD(κ,∞) ⇒ Var` 的教科书引用（Bakry–Gentil–Ledoux）**未逐字核页码**；§3 的判决只依赖**自证精确反例 `[T]`**，不依赖该书。
5. §4 里 2 次 `semanticscholar` 调用 = **HTTP 429**，记为 `TOOL-FAIL`，**未计入任何 0 命中**；结论改由 OpenAlex + web_search 支撑（均 200）。

**标签汇总**：§0/§6 判词 `[NG]`（候选 1/2/3）与 `[M]`（RC-BL）；§1 换名论证 `[T]`（自证代数）+ 书目 `[T-书目]` + 公式 `[C]`；§2 `[C]/[T-书目]`；§3 反例 `[T]`（手算）、`κ>0` 事实 `[C]/[T-standard]`；§4 检索式与命中数为**实测**（表内逐条）；§4.1 Kolesnikov 摘要 `[C]`（**摘要已读**，正文未读）、Helffer–Sjöstrand `[T-书目]`+`[C]`片段、SD 片段 `[C]`（付费墙，未逐字核）；§5 `[M]`；§7 未闭合 5 条。

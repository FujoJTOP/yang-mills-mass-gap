# GX2-TRANSPORT — 输运-熵（Talagrand `T₂`）/ 局部化：`λ_min<0` 时还能不能给模数？

> 2026-11-26 · 子代理 `gx2-transport` · **只写本文件**（未改任何既有文件、未 `git commit`、**零 MC**）· **No QED** · **Clay 0%**。
> 题面：`GH1-SPECTRAL.md §2.6`（`:232–239`）已证带模数的输入只有四类、无第五类；其中 ② = 一致 Poincaré/LSI 是墙（Savvidy `λ_min = −g|B| < 0`，`GI4-FORMS.md:40,95`；`G4-GENERAL-BOUND.md:201–204`）。本路判：**能不能找到一类"模数不要求 `Hess V ⪰ 0`"的不等式**（赌注：Talagrand `T₂` / 局部化一维归约）。
> 靶（逐字，`GH1-SPECTRAL.md:172`）：`(E) ∃Δ₀>0 ∀(L,a): log R_t ≥ (aΛ)·Δ₀`。
> 标签：`[T]` 真定理/恒等式 · `[T-标准]` 教科书级（未开卷） · `[C]` 条件/片段级证据 · `[R]` 外部读到摘要或片段 · `[NG]` 具名 no-go · `[O]` 精确缺失步 · `[M]` 口径 · `[E]` 本轮精确算术（式在文内，可复跑） · `[X]` 抓错。
> **读入（全部读毕）**：`GH1-SPECTRAL.md:160–366`；`GI4-FORMS.md:38–97`；`G4-GENERAL-BOUND.md:186–206`；`S1-CONDITIONAL-POINCARE.md:1–48`；`O01-MARGINAL-BL.md:100–258`；`U09-M6.md`（grep 全 30 命中，含 `:11,:15,:125–127,:150`）；`O17-OT.md:36–139`；`B01-RP.md:1–14`。
> **匹配 skill**：本任务是数学判定 + 逐行 `file:line` 核对 + 外部文献定位；`<available_skills>` 里的学术类 skill（`academic-paper`/`deep-research`/`academic-pipeline` 等）面向**写作/评审流水线**，与"判一条不等式能否过 `(E)`"不匹配 ⇒ **未加载 skill**，用 `read`/`grep`/`pwsh`/`web_search`/`web_fetch` 直做。

---

## §0 判词（七条，先写死）

### 0.1 一句话

**本路的赌注含一个类别错误，而纠错之后赌注自动作废。**题面写"Poincaré/LSI 都要求 `Hess V ⪰ 0` 型条件"——这句**只在证明方法层为真**（Bakry–Émery / Holley–Stroock / Prékopa–BL 是**充分条件**），在**陈述层为假**：Poincaré、LSI、`T₂` 的**陈述里根本没有 `Hess`**。所以"有没有一类不等式，其模数不要求这个"实际是在问"**有没有一条不用凸性的、给常数的证明方法**"——这是一个**换过名字的同一个问题**，不是第五类。纠错后：`T₂` 有凸性无关的正定理（§1.2，**有**，不是没查到），但它产出的是 **Poincaré**，而 `(E)` 只需要 Poincaré 常数（`GH1-SPECTRAL.md:209`）⇒ **`T₂` 比墙更严，不是更松**。

### 0.2 Q1 · `T₂` 在 `λ_min<0` 时有没有已知定理给常数 `>0`？——**有**，但正定理的产出就是墙

| 问 | 判 | 依据 |
|---|---|---|
| 存在"弱于凸"的 `T₂` 正定理吗 | **存在（三条形状）**：判据/Lyapunov 型（Cattiaux–Guillin，`arXiv:math.PR/0312081`，摘要逐字）；有界扰动型（Holley–Stroock，经 `O01-MARGINAL-BL.md:202` `[T]`）；张量型（乘积测度，2D 即此） | §1.2 |
| 这些正定理把 `T₂` 放在哪一格 | **`LSI ⟹ T₂ ⟹ Poincaré`，两箭头皆严格**（`T₂` 无 LSI 的例子：Cattiaux–Guillin, *JMPA* **86** (2006) 341–361，摘要逐字） | §1.3 |
| ⇒ 它绕过墙了吗 | **没有**。它落在墙的**内部**：`T₂ ∈ (Poincaré, LSI)`，而墙 ② 覆盖整段；`(E)` 的下界只需 Poincaré 常数 ⇒ `T₂` 常数是**墙内的一个更难点** | §1.5 |
| 判 | **`[NG]`（类型级，不是文献级）** | — |

### 0.3 Q2 · ★ 局部化一维归约是否只要求"边际凸"？——**否**，且对本靶是**空操作**

1. **前提层（外部）**：Borell 局部化 / needle decomposition 的前提是**联合** log-凹（等价 `CD(0,∞)`），不是边际凸：KLS 1995 的 needle 是"**exponential needle**"、"consider any **log-concave** function `F`"（片段逐字）；Klartag 2017 的推广前提是"**Ricci curvature bounded from below**"（片段逐字）；`arXiv:2108.07193` 片段逐字："the condition `CD(0,∞)` is **equivalent to log-concavity** of the considered measure"。⇒ **`[R]/[C]`**，§2.2。
2. **反例层（本文件新建，精确）**：`p=q=Uniform(−1,1)`（两个边缘都是**最强**的 log-凹：密度常数）配联合密度 `f=1+ρ(x²−⅓)(y²−⅓)`（`ρ=½` 或 `1`）：`f>0`（网格 `min=0.888889`），边缘**精确**为 Uniform（`∫(y²−⅓)/2 dy = −8.3e−12`，`x=0.7` 处边缘 `=1.000000000`），而 log-密度 Hessian 在 `(0.5,0.5)` 有一个正特征值 `+0.411777` ⇒ **联合不 log-凹**。⇒ "边际凸 ⇏ 联合凸" `[T]+[E]`，§2.4。
3. **反例层（仓内，更贴靶）**：`O01-MARGINAL-BL.md:159–166` 的精确反例——单 plaquette `SU(2)`,`θ=π`,`v=e₃`：切片测度 **log-凹**（`λ_⊥ ≥ inf Φ'' > 0`）而**同一切片的边缘** `G''(0)=W₀−V₀β²<0`（`β>4.30660785`）⇒ `(CBL)` 前提成立、`(BL-1D)` 结论假 `[T]`。机制逐字："**BL 是切向算子，`G''` 的缺口是法向方差 —— 这是类型不匹配，不是常数不优**"。
4. **归约层（本文件命题 L′）**：`G''` 是**沿固定方向的边缘**泛函；局部化的收益全在"从 needle 条件测度**拼回**联合"那一步，而拼回 = Prékopa = 凸性。⇒ 对本靶，1D 归约把问题化成它自己（`[T]`，§2.3）。
5. **对象层**：Savvidy 的不稳定方向**本身就是一维的**（`k=0` 粗模、单个横向磁模，`GX2-WINDOW.md:9–11`）⇒ "n 维降到 1 维"对一个"已经 1 维"的不稳定方向**零信息**（§2.5）；且即便成功，产出的是 `G''` 的模数 = **边缘盲**对象（`U09-M6.md:125–127`），不是 `(E)`。
6. **判**：**局部化不是第五类**。它把维度从 `n` 降到 `1`，而墙在**曲率**坐标上；**两个坐标正交** ⇒ 墙一格没动。`[NG]`

### 0.4 Q3 · 即使有模数，能否给 `(E)`（`(L,a)`-一致）？——**三条独立死因，全部 `[T]`**

| # | 死因 | 数 | 出处 |
|---|---|---|---|
| i | **常数 = 逆隙本身**（Gauss 上精确相等）：`C_T₂^♯(N(0,Σ)) = 2λ_max(Σ) = 2/Δ_min = 1/Λ` | 自由无质量场 `C_T₂ = L²/(2π²)`：`L=4/16/64` → `0.8106 / 12.9691 / 207.5058` ⇒ **不 `L`-一致** | §1.4、§3.2（`[T]+[E]`） |
| ii | **`a`-一致 ⇒ 反靶**：`C_T₂·a²` 有界 ⟺ `Δ_min ≥ 2c/a` ⟺ `m·a ≥ c` | `B04.md:92–97`、`B01-RP.md:34`；自由有质量场见证 `a m = 2asinh(am/2) → 0` | §3.3（`[T]`） |
| iii | **去凸性的代价 = `e^{osc}`**：Holley–Stroock 在 `β=2.5, L=4` 给 `4βN_□ = 15360`，常数 `≈10^{6671}` | 第二座（广延量）**未付**，与 `GS2-PATTERN.md:5` 同座 | §3.4（`[T]`，承 `O01:202`） |

### 0.5 Q4 · 判词

**不是第五类。** 精确位置：**`T₂` 与局部化都不在 `GH1 §2.6` 的四类之外**。
- `T₂` = **② 的内点**（`Poincaré ⊂ T₂ ⊂ LSI`）；`GH1:209` 已判"只要 Poincaré 常数即 `(E)`" ⇒ `T₂` 是一条**更难的充分路线**，不是绕道。
- 局部化 = **② 的证明技术**（其前提是联合 log-凹 = `Hess V ⪰ 0`，与 ② 的前件逐字相同，`§2.2`）。
- ⇒ `GH1 §2.6` 的"四类"**不多不少**；`GH1:239` 的"本轮没有发现第五类"在本路**被独立复核**（且本路给出了它**为什么**不该有第五类：第四类与第五类的差别不在不等式，在"常数住不住在模数里"，`GH1:216`）。

### 0.6 酸测（§4）

| 模型 | `T₂` / 局部化 的输出 | 判 |
|---|---|---|
| **自由场（无质量）** | `T₂` **成立**（`C=L²/(2π²)` 有限），而 `m=0` | **假阳性**（判别器失败）；同时 `C` 随 `L` 发散 ⇒ 不 `L`-一致 |
| **2D YM** | `T₂` 经**张量**成立，`C=O(1)`；`m=3g²/8>0` | **退化通过**：机制是"乘积测度"（2D 无共享 link），4D 无对应物；且与自由场**同输出** ⇒ 零判别 |
| **M6** | `C_T₂ = 1/m_L`：`L=40` 给 `81.3206` ⇒ 不 `L`-一致；且 M6 是 Gauss ⇒ BL/局部化**取等** ⇒ 零 slack | **杀 `L`-一致** + 零模数 |
| **单 plaquette** | 无方向可归约（局部化空转）；`K⪰0` 在物理点 `θ=π` **假**（`d²S/ds²=−(β/4)sin(s/2)<0`，`O01:136`） | **静默** |

**签名** = `(假阳性, 与自由场同输出, L-不一致, 单片静默)` ⇒ **不过酸测**。与 `F06-KLATOM.md:82–85` 同型方法论、与 `GH1:291–304` 的四模型签名一致。

### 0.7 交付

本轮**新增**：（a）`T₂` 在 `[Poincaré, LSI]` 内点的**精确定位**（含 Gauss 上 `C_T₂ = 1/Λ` 的两行证）；（b）**命题 L′**（对本靶 1D 归约是空操作）；（c）新建**边际-联合反例**（Uniform 边缘 + 非 log-凹联合，精确数）；（d）外部三条正定理的 DOI/arXiv 与"读到哪一层"；（e）抓错一条（`GH1c-LIT` 的 "LSI ⟺ T₂"）。**未证明任何关于 4D YM 隙的东西。Clay = 0%。**

---

## §1 `T₂` 判定

### 1.1 先纠题面的一个措辞（本路判死的根源，`[T]` 逻辑）

题面：**"而 Poincaré/LSI 都要求 `Hess V ⪰ 0` 型条件。"**

这句话必须拆成两层，否则整条路都会被误判：

| 层 | 内容 | 真值 |
|---|---|---|
| **陈述层** | `Poincaré(C)`：`Var_μ(f) ≤ C·E_μ|∇f|²`；`LSI(C)`：`H(ν|μ) ≤ (C/2)I(ν|μ)`；`T₂(C)`：`W₂²(μ,ν) ≤ C·H(ν|μ)`。三者的**陈述里没有 `Hess V`**，只出现 `μ`、`f`、`∇` | **"要求 `Hess V ⪰ 0`" 为假** `[T]` |
| **证明层** | 要**证出常数**，已知机器：Bakry–Émery `Γ₂ ≥ ρΓ`（要 `Ric + Hess V ⪰ ρ g`）、Brascamp–Lieb/Prékopa（要 `Hess V ⪰ 0`）、Holley–Stroock（不要凸性，要 `‖osc S_W‖ < ∞`）、Otto–Villani `HWI`（要 `κ>0`，`O17-OT.md:52,66`） | **"需要凸性" 为真** `[T]` |

⇒ **"有没有一类不等式，其模数不要求这个"** = **"有没有一条不用凸性、却能给出常数的证明方法"**。这不是"第五类不等式"，是**同一个问题的第二种问法**。本路据此给判：**若答案只能是"换证明方法"，那它必须同时回答"新方法的常数是不是 `(L,a)`-一致"**——而这一条与不等式叫什么名字无关（§3）。

> 这正是 `GH1:216` 的判据在另一坐标上的复现：**"分界落在常数住不住在模数里，不落在用哪个定理。"**

### 1.2 外部已知定理（"弱于凸"确实存在，逐条给读到哪一层）

| # | 定理 | 陈述（逐字/要点） | 前提 | 读到哪层 | 对本靶 |
|---|---|---|---|---|---|
| **R1** | **Cattiaux–Guillin 判据**，*A criterion for Talagrand's quadratic transportation cost inequality*，**arXiv:math.PR/0312081**（HAL `hal-00133642`） | 摘要逐字：*"We show that the quadratic transportation cost inequality `T₂` is **equivalent to both a Poincaré inequality and a strong form of the Gaussian concentration property**. The main ingredient … a new family of inequalities, called **modified** quadratic transportation cost inequalities in the spirit of the modified logarithmic-Sobolev inequalities by Bobkov and Ledoux, that are shown to hold **as soon as a Poincaré inequality is satisfied**."* | **无凸性前提**（判据型） | **读到摘要全文** | **决定性**：正定理存在，但它把 `T₂` **等价到 Poincaré** ⇒ 见 §1.5 |
| **R2** | **Cattiaux–Guillin**，*On quadratic transportation cost inequalities*，*J. Math. Pures Appl.* (9) **86** (2006), no. 4, **341–361**，DOI `10.1016/j.matpur.2006.06.003` | 摘要逐字：*"we introduce new families of inequalities (for quadratic transportation cost and for relative entropy) that are shown to be **equivalent to the Poincaré inequality**. This allows us to give some **examples of measures satisfying `T₂` but not the logarithmic-Sobolev inequality**."* | 无凸性前提 | **读到摘要**（ScienceDirect 摘要片段；正文付费墙，**未逐字核**） | **`LSI ⟹ T₂` 不可逆** ⇒ `T₂` **严格弱于** LSI（这是本路唯一真正的"松动"，但只松一格，§1.3） |
| **R3** | **Gentil–Guillin–Miclo**，*Modified logarithmic Sobolev inequalities and transportation inequalities*，*Probab. Theory Related Fields* **133**(3) (2005) **409–436**（HAL `hal-00001609`） | 题面所说"曲率下界为负时的修正"在文献中确实以 **modified (log-)Sobolev / modified transport** 形态存在 | 修正项 = **可加尾项**（不是把负曲率变成正曲率） | **只读到书目记录**（PDF 抓取返回 `unsupported content type "application/pdf"`，**正文未读**）⇒ 不据此作判 | 修正项在格点上就是"第二座"（广延/`a`-发散），§3.4 |
| **R4** | **Otto–Villani**，*Generalization of an inequality by Talagrand and links with the logarithmic Sobolev inequality*，*J. Funct. Anal.* **173**(2) (2000) **361–400** | `LSI ⟹ T₂`（本仓已用：`GH1c-LIT.md:50,158`） | 无 | 书目级（卷/页多源一致；**DOI 未核**、正文未读） | 上箭头；仓内 `[R]` |
| **R5** | **Holley–Stroock**，*J. Stat. Phys.* **46** (1987) **1159–1194**，DOI `10.1007/BF01011161` | 有界扰动：`Var ≤ e^{2osc}λ_0^{-1}E|∇f|²` | **无凸性**，但常数 `e^{2osc}` | 书目级（`O01:236` 已核 DOI） | 见 §3.4：`10^{6671}` |

**⇒ Q1 的答案（`AGENTS §1` 四项见 §5.1）**：**有**。`T₂` 的"弱于凸"正定理**不是缺的**，是**现成的**（R1 判据 / R5 扰动 / 张量）。**所以本路不是"`[NG]`@文献"，是"`[NG]`@类型"**——见 §1.5。

### 1.3 等价链与严格性（`[T-标准]` + `[R]`）

```
LSI  ⟹  T₂  ⟹  Poincaré
 (R4)      (R1)
两个箭头都严格：
  · T₂ ⇏ LSI：R2 摘要逐字给出"满足 T₂ 而无 LSI"的例子 `[R]`
  · LSI ⇏? T₂：不成立（R4 是单向蕴含）
  · Poincaré ⇏ 经典 T₂：R2 的"新族（= **modified** 版）⟺ Poincaré"只说明 Poincaré 够到 **modified** 版本；
    `T₂` 是否被 Poincaré 蕴含，摘要未说 ⇒ **本文件不用这一步**（判词只依赖前两个箭头）
```
- 仓内已有同判：`GH1c-LIT.md:158` 逐字"`von Renesse–Sturm` 型结果是**等价链**……**等价意味着：任何一端都不能比另一端更容易拿到**"。本路的**修正**是：这条链不是"全等价"，而是**单向 + 严格**，且 `T₂` 落在**中间**。
- 负曲率侧的桥：仓内已判 `O17-OT.md:66` 逐字"`CD(κ,∞)` 与'方差 ≤ 曲率倒数'之间的桥，在 `κ<0` 时**断在 HWI 的二次项上**"；`:139` 更给出 `κ<0` 的**精确反例**（`V=−x²/2+x⁴/4` 使 `Var ≤ −1` 荒谬）`[T]`。⇒ **"负曲率修正"在"`Var ≤ 1/κ`"这一格上是假命题，不是缺定理。**

### 1.4 Gauss 上 `T₂` 常数的**精确值**（本路唯一的 `[T]` 新算，两行）

> **命题 T（`[T]`）** 设 `μ=N(0,Σ)` 非退化。则 `T₂` 的**锐**常数为
> ```
> C_T₂^♯(N(0,Σ)) = 2·λ_max(Σ) = 2/λ_min(Σ^{-1}) = 2/Δ_min = 1/Λ(μ) ,
> ```
> 其中 `Δ_min` = 精度矩阵（= `Hess V`）的最小本征值，`Λ` = `U09-M6.md:150` 的归一化（`Λ = ½∫|∇f|²/Var = Δ_min/2`，`B01-RP.md:16` 的 `2Λ`）。

*证.*
1. **下界（锐性）**：取 `ν=N(m,Σ)`（纯平移）。则 `W₂²(μ,ν)=|m|²`、`H(ν|μ)=½ mᵀΣ^{-1}m`。`T₂` 要求 `|m|² ≤ (C/2)mᵀΣ^{-1}m` 对**一切** `m`，取 `m` 沿 `Σ` 的最大本征方向 ⇒ `C ≥ 2λ_max(Σ)`。
2. **上界（成立）**：`μ = Σ^{1/2}#γ`（`γ` 标准 Gauss）。标准 Gauss 满足 `T₂(2)`（Talagrand 1996，`[T-标准]`），而 `W₂²(A#α,A#β) ≤ ‖A‖²W₂²(α,β)`、熵在推前下不变 ⇒ `T₂` 对 `μ` 以 `2‖Σ^{1/2}‖² = 2λ_max(Σ)` 成立。∎

**数值复核**（`python -I`，纯 `math`，无第三方包）：
```
(B) sup_{grid} W₂²/H = 1.999994  (理论 ≤ 2；m=0.001,s=1 给 2.000000)
    (m,s)=(0.001,1.0000) → 2.000000 ; (2.0,2.0) → 1.781355
```

**读法（本路最重要的一行）**：**`T₂` 的常数在 Gauss 上逐字等于"逆隙"。**它不是隙的**代理**、不是隙的**上界数据**，而是**同一个数**。⇒ 想用 `T₂` 给 `(E)`，等于想用 `1/Δ_min` 给 `Δ_min` 的下界——**命题与靶同一**。

### 1.5 判：`T₂` 在 `[Poincaré, LSI]` 的内点 ⇒ 不是第五类

| 事实 | 依据 | 后果 |
|---|---|---|
| `Poincaré ⊂ T₂ ⊂ LSI`（严格） | §1.3 | `T₂` 是**墙内的一个点**，不是墙外的路 |
| 墙 ② = "一致 Poincaré / LSI / 曲率" | `GH1-SPECTRAL.md:235` | 墙覆盖**整段** |
| `(E)` 只需 Poincaré 常数 | `GH1-SPECTRAL.md:209`（"(b) 变分 + Poincaré … **只要常数即 `(E)`**"） | 任何**强于** Poincaré 的一致陈述都**足够**推 `(E)` |
| ⇒ 证出一致的 `T₂` ⇒ 得 `(E)` | 三段论 | **是**一条**合法但更难**的路线 |
| ⇒ `T₂` 是否能**绕过**墙 | — | **不能**。赌注想"求弱"，实际是"求强" |

**一句话判词**：**赌注只把不等式从 `LSI` 降了一格，停在 `Poincaré` 的**上面**一格；而墙是 `Poincaré` 那一格（更准确：覆盖整段）。⇒ 没越过。**

### 1.6 抓错 `[X]`：`GH1c-LIT.md:50,158` 的 "`LSI ⟺ T₂`"

- 逐字（`GH1c-LIT.md:50`）："2010s 输运-熵（vRS/Kuwada/Erbar–Maas）| **是等价，不是工具** | `LSI ⟺ T₂` 是**双向**陈述"；`:158` 同旨。
- **按 R2 摘要逐字**（Cattiaux–Guillin 2006）："measures **satisfying `T₂` but not** the logarithmic-Sobolev inequality" ⇒ **`T₂` 不蕴含 LSI** ⇒ "双向"**假**。
- **正确写法**：`LSI ⟹ T₂ ⟹ Poincaré`，**两箭头皆严格**。
- **影响**：`GH1c-LIT` 的**结论**（"输运-熵不产生新界"）**不变**——甚至更强（一个更弱的中间格都不产生新界，全等价的说法反而更弱）。只改**措辞与等级**：`:50` 的 `⟺` 应为 `⟹`，标 `[R]` 而非 `[T]`。
- **不做什么**：本文件**未改** `GH1c-LIT.md`（写域限制），交 `lead` 记账。

---

## §2 ★ 局部化：一维归约是否只要求"边际凸"？

### 2.1 三个对象必须分开（否则会把三种不同的"局部化"混成一条）

| 名字 | 做什么 | 前提 | 对本靶 |
|---|---|---|---|
| **Borell 局部化 / needle decomposition**（Borell 1975；Kannan–Lovász–Simonovits 1995；Klartag 2017） | 把 `n`-维积分不等式**归约到一根"针"**上的一维不等式 | **联合 log-凹**（`CD(0,∞)`） | §2.2–2.5：**前件在 Savvidy 方向假 + 对本靶空操作** |
| **Prékopa / Brascamp–Lieb** | 联合 log-凹 ⇒ 边缘 log-凹、`Var` 界 | **联合** `Hess V ⪰ 0` | 仓内已判：`W_geo ≥ 0` 在 `β=0` **取等**（`G4-GENERAL-BOUND.md:196–199`），`β>0` 不转移 `[NG]`（`:201–204`） |
| **Bakry–Émery / `CD(κ,∞)`** | `Γ₂ ≥ ρΓ ⇒ 隙 ≥ ρ` | `Ric + Hess V ⪰ ρg` | 仓内已判 `ρ ≈ −2.3 < 0` ⇒ vacuous（`O18-CHEEGER.md:32`；`O08-CD.md:121`） |

**关键**：题面把前两个当成同一条路。它们**前提相同**（联合 log-凹），**产出不同**（一维归约 vs 边缘界），**对本靶都死在同一个字上：联合**。

### 2.2 一维归约的前提，逐字（外部）

| 源 | 逐字/片段 | 读到的层 |
|---|---|---|
| **Kannan–Lovász–Simonovits (1995)**，*Isoperimetric problems for convex bodies and a localization lemma*，*Discrete Comput. Geom.* **13**(3–4) **541–559**，DOI `10.1007/BF02574061` | *"The Localization Lemma … reduces the problem of proving inequalities between integrals in `ℝⁿ` to inequalities between **one-dimensional** integrals"*；needle 的定义是 **"exponential needle"**（区间 + 指数 `γ`）；证明中 *"consider any **log-concave** function `F` … `Fℓ^{n−1}` is also **log-concave**, and hence Lemma 2.4 implies that there exists an exponential needle `E` violating (b)"* | 片段级（PDF 抓取按摘要片段读到） |
| **Klartag (2017)**，*Needle decompositions in Riemannian geometry*，*Mem. Amer. Math. Soc.* **249** (2017), no. **1180** | *"The localization technique from convex geometry is generalized to the setting of Riemannian manifolds whose **Ricci curvature is bounded from below**"* | 片段级（书目 + 摘要片段） |
| **`arXiv:2108.07193`**（*Leaves decompositions in Euclidean spaces*） | *"Note also that the condition **`CD(0,∞)` is equivalent to log-concavity** of the considered measure."* | 片段级 |
| 一份博士论文（HAL `tel-04577476`，Bizeul 2023） | *"Almost surely `μ_ω` is a 'needle', that is **a one dimensional log concave probability** or a Dirac"*；*"Then by the Localization Lemma (Lemma 4.1) there is a **one-dimensional log-concave needle**"* | 片段级 |
| **Borell (1975)**，*Convex set functions in d-space*，*Period. Math. Hungar.* **6**(2) **111–136**，DOI `10.1007/BF02018814` | 局部化原理的原始出处 | 书目级（正文未读） |

**⇒ 前提是 `CD(0,∞)` = 联合 log-凹 = `Hess V ⪰ 0`（联合），不是"边际凸"。** `[R]/[C]`

### 2.3 命题 L′：对本靶，一维归约是**空操作** `[T]`

> **命题 L′.** 设泛函 `F` 只依赖测度沿固定方向 `v` 的**边缘** `ρ_v`（本靶：`G'' = −∂²log ρ_v`）。则任何"把 `n`-维测度化成 1-维族"的定理，对 `F` 的**信息增量 = 0**。
>
> *证.* 一维归约的产物是 **各 needle 上的条件测度** `{μ_ℓ}`，其结论形如"`∀ℓ`：`I(μ_ℓ)` 成立 ⇒ `I(μ)` 成立"。而 `F` 依赖的是**边缘** `ρ_v`；要把 `{μ_ℓ}` 拼回 `ρ_v` 必须走 Prékopa（联合 log-凹 ⇒ 边缘 log-凹），而**一旦**有联合 log-凹，`ρ_v` 的 log-凹性由 Prékopa **直接**给出，**根本不需要 needle**。⇒ needle 这台机器对 `F` 的**全部**贡献都在"拼回"那一步，而拼回的成本 = 凸性（= 墙）。∎
>
> **推论（决定性）**：局部化若要给 `G''` 一个模数，它必须**先**付出 `Hess V ⪰ 0`（联合）；付了这一步，`G''≥λ` 已由 Prékopa 得到，与不定 Hessian 的假设**直接矛盾**。⇒ **局部化在逻辑上不可能成为"绕过墙"的第五类。**

**注（不夸大）**：命题 L′ **不是**说局部化无用。对**联合型**泛函（等周 profile、集中函数、`W_geo` 这类），归约是真的、且给定量界（这正是 KLS 的用途）。本靶偏偏是**边缘型**。

### 2.4 "边际凸不够"：两个精确反例

**(a) 本文件新建（概率意义的"边际"）`[T]+[E]`**

构造：`p=q=Uniform(−1,1)`（**两个边缘都 log-凹，且是最强的一种：密度常数**），联合
```
f(x,y) = 1 + ρ(x²−⅓)(y²−⅓),   (x,y) ∈ [−1,1]²,  ρ ∈ {½, 1}
```
- **正性**：`(x²−⅓)(y²−⅓) ∈ [−2/9, 4/9]` ⇒ `f ≥ 1−ρ·2/9 > 0`（网格实测 `min f = 0.888889`（`ρ=½`）/ `0.777778`（`ρ=1`））。
- **边缘精确为 Uniform**：`∫_{−1}^{1}(y²−⅓)·½ dy = ⅓−⅓ = 0`（数值 `−8.3e−12`）；`x=0.7` 处边缘 `=1.000000000`。
- **联合不 log-凹**：`log f` 的 Hessian（解析式见下）在 `(0.5,0.5)` 有本征值 `+0.411777`（`ρ=½`）/ `+0.813888`（`ρ=1`）；`(0.6,0.2)` 亦 `+0.149793`。

```
g_xx = 2ρB/D − 4ρ²x²B²/D² ,  g_yy = 2ρA/D − 4ρ²y²A²/D² ,  g_xy = 4ρxy/D² ,
A = x²−⅓ ,  B = y²−⅓ ,  D = 1+ρAB .
```
**数值**（`python -I`）：
```
ρ=½: (0.00,0.00) eig=(−0.315789,−0.315789) log-凹 ✓ ; (0.50,0.50) eig=(−0.581315,+0.411777) ✗
     (0.80,0.80) eig=(−0.929621,+1.405601) ✗ ; (0.60,0.20) eig=(−0.448755,+0.149793) ✗
ρ=1: (0.50,0.50) eig=(−1.158621,+0.813888) ✗ ; min f = 0.777778 ; 边缘 = 1.000000000
```
⇒ **"两个边缘都 log-凹"对联合 log-凹一点约束都不给。** `[T]`
⇒ 这直接否掉"localization 只要求边际凸"的**概率读法**。

**(b) 仓内已有（方向/切片意义的"边际"，更贴靶）`[T]`**

`O01-MARGINAL-BL.md:159–166`（单 plaquette `SU(2)`、`θ=π`、`v=e₃`、`ξ=0`）：
- 切片 `{ξ₃=0}`：`S_W ≡ β` 常数 ⇒ 切片测度 `∝ J d²ξ`，`J=(sin h/h)²`、`Φ=−log J` 凸 ⇒ **切片 log-凹，`λ_⊥ > 0`**。
- **同一切片的边缘**：`G''(0) = W₀ − V₀β²`，`β > β* = 4.30660785` 时 `<0`（精确：`β=6` 给 `−0.193019`，`O01:148`）。
- ⇒ `(CBL)` 前提成立、`(BL-1D)` 结论假 `[T]`。

**⇒ 两种读法的"边际凸"都被精确反例否掉，且第二个反例就长在靶上（`G''`）。**

### 2.5 第二与第三死因（本题特有，不用文献）

**(i) 不稳定方向已经是一维的 `[T]`**
Savvidy 的负模取在 `k=0` 粗模、物理（横向）子空间上（`GX2-WINDOW.md:9–11`：`G″ = −D²+2adF_B`，粗糙模 `k=0`），`λ_min = −g|B|`（`GI4-FORMS.md:95`）。一个"把 `n` 维降到 1 维"的定理，对一个**本来就只有一维**的不稳定方向，**可降的维度 = 0**。
- **更形式化的说法**：`G''` 是**单个模坐标**上的二阶导（`U09-M6.md:50`：`G''(q) := −∂_q²log ρ(q)`）；它的"一维化"就是它自己。

**(ii) 局部化的自然产物是"边缘盲"对象 `[T]`**
`U09-M6.md:125–127` 的**边缘盲引理**逐字：设 `F` 是 `ρ_a` 的任意泛函（含整个边缘函数、含一切矩），则 `F(M6) = F(HO_ω)`；若 `F` 在谐振子上 `>0`，则 `F(M6) > 0`，与 `m(M6)=0` 矛盾（实测边缘差 `2.259e−14`）⇒ **任何只依赖单坐标边缘的量都不可能当靶**。
⇒ 就算局部化**成功**给出 `G'' ≥ λ`，那也是**边缘泛函**的模数 ⇒ `M6` 一票否掉 ⇒ **不是 `(E)`，是 `(H1)` 的一个更强形式**（`O01:217`：`(BL-1D)ᵤ ⇒ H1 [T]`，且严格强于 H1）。

### 2.6 Brakke（题面第三条）

题面列 "Brakke / 局部化"。**两种读法都判 `[NG]`**：
- **Brakke（平均曲率流）**：`Brakke` 不等式 / Huisken 单调性公式是**几何流**的单调性，前提是**平均凸**（`H ≥ 0`）——与 Savvidy 是**同一个符号条件**，不是它的替代；且它输出的是"加权面积"的单调性，**没有测度侧的熵-输运模数**。本靶里也没有可挂的 MCF 结构。检索：见 §5.3 `#2`，**相关命中 0**。
- **`Brascamp`（若为笔误）**：已由 `GI4-FORMS.md:40,95` 判 `[NG]`（前提 `Hess V ⪰ 0` 假）。

**判**：`[NG]`（本靶）/ `[M]`（口径）。

### 2.7 §2 判词

| 问 | 答 | 级 |
|---|---|---|
| 局部化一维归约只要求边际凸吗 | **否**：前提是**联合** `CD(0,∞)`（外部四源逐字） | `[R]/[C]` |
| "边际凸"能否替代 | **否**：两个精确反例（本文件 `Uniform` 构造；仓内 `O01 §4.3` 切片构造） | `[T]+[E]` |
| 对本靶，1D 归约有多少信息 | **0**（命题 L′：收益全在"拼回"= Prékopa = 凸性） | `[T]` |
| 负方向需要"降维"吗 | **不需要**：它已经是一维（`k=0` 粗模） | `[T]` |
| 即使成功，产出是 `(E)` 吗 | **不是**：是边缘盲对象 `G''`（`M6` 否掉） | `[T]` |
| **⇒ 是否第五类** | **不是** | `[NG]` |

---

## §3 一致性：即使有模数，能否给 `(E)`（逐条判 `(L,a)`-一致）

`(E)` 的要求逐字（`GH1-SPECTRAL.md:172`）：`∃Δ₀>0 ∀(L,a): log R_t ≥ (aΛ)Δ₀`。下面四条独立判。

### 3.1 谱量的同一性（`[T]`，§1.4 的推论）

`C_T₂^♯ = 2/Δ_min = 1/Λ`。⇒ **`(L,a)`-一致的 `T₂` 常数 ⟺ `(L,a)`-一致的 `Δ_min` 下界 ⟺ `(E)`**（末一等价承 `GH1-SPECTRAL.md:209`：只要 Poincaré 常数即 `(E)`）。
- **判**：这不是"路线的强度问题"，是**命题与靶同一**（换名）。`[T]`

### 3.2 `L`-一致：自由无质量场直接给反例（`[T]+[E]`）

`T⁴_L` 环上无质量自由场（去零模）：`Δ_min = (2π/L)²` ⇒ `C_T₂ = L²/(2π²)`。
```
L= 4  Δ_min=2.467401  C_T₂=  0.8106
L= 8  Δ_min=0.616850  C_T₂=  3.2423
L=16  Δ_min=0.154213  C_T₂= 12.9691
L=32  Δ_min=0.038553  C_T₂= 51.8764
L=64  Δ_min=0.009638  C_T₂=207.5058
```
⇒ **任何 `L`-一致的 `T₂` 常数都被自由场否掉。** 而自由场 `m=0`（`GH1:275–277` 的三读表）⇒ 自由场同时是"`T₂` 成立"的**假阳性**见证（§4.1）。

### 3.3 `a`-一致：`a`-一致 ⇒ 反靶（`[T]`）

`C_T₂·a² = 2a²/Δ_min`。要求 `C_T₂` 在**格点单位**上有界 ⟺ `Δ_min ≥ 2c/a` ⟺ `m·a ≥ const`。
- 这正是 `GH1-SPECTRAL.md:211`（(d) 的出口 (ii)）逐字判过的"**反靶**"：`B04.md:92–97` 的"格点单位地板 = 相关长度 `O(a)` 的平凡 UV 陈述"（`B01-RP.md:34`）。
- **见证它不存在**：自由有质量场 `a m = 2·asinh(am/2) → 0`（`GC1-NEW.md:44`，`GH1:275`）⇒ 若 `C_T₂` 有 `a`-一致界，则自由场也必须有 `m·a ≥ c` ⇒ 与 `2asinh(am/2)→0` 矛盾。
- **判**：`[T]`。**这条与 §1.4 合起来是本路最硬的一刀**：`T₂` 的常数**不是**"另一个量"，它把 `m` 的坐标**逐字**搬了过来。

### 3.4 体积/广延：去凸性的代价是"第二座"（`[T]`，承 `O01`）

| 去凸性的手段 | 常数 | 数（`β=2.5, L=4`） |
|---|---|---|
| Holley–Stroock 有界扰动 | `e^{2·osc(S_W)}`，`osc ≤ 2βN_□`，`N_□ = 6L⁴` | 指数 `4βN_□ = 15360` ⇒ `≈ 10^{6671}`（`O01-MARGINAL-BL.md:202`，`[T]+[O]`） |
| 张量（乘积测度） | `max_□ C_□` | **只在测度真分解时可用**：2D 成立，4D **不成立**（plaquette 共享 link，`§4.2`） |
| 判据/Lyapunov（R1） | 由 Lyapunov 函数定 | 需一个显式的 `V`-型漂移；在 `L`-一致意义下等价于 `§3.2` |

- `HS` 的常数是**广延量**（`∝ N_□ = 6L⁴`）⇒ 正是 `GS2-PATTERN.md:5` 的"**第二座**（`a`-发散/广延的量）"，且第二座**没付**。
- **判**：`[T]`（HS 的算术承 `O01`）。"把凸性换成体积指数"= 用第二座换第一座，不是消掉墙。

### 3.5 §3 判词

| 层 | 判 | 级 |
|---|---|---|
| `T₂` 常数与 `(E)` 的关系 | **同一命题**（`C_T₂^♯ = 1/Λ`） | `[T]` |
| `L`-一致 | **自由场反例**：`C = L²/(2π²) → ∞` | `[T]+[E]` |
| `a`-一致 | **⇔ 反靶**（`m·a ≥ c`）；自由场见证其不存在 | `[T]` |
| 去凸性的代价 | `e^{2osc}` = 广延第二座（`10^{6671}`） | `[T]` |
| **能否给 `(E)`** | **只有把 `(E)` 当假设才能给 `(E)`** ⇒ 循环 | `[NG]` |

---

## §4 酸测（四模型；精确数，无 MC）

### 4.1 自由场（**必须失效**）—— 判：**假阳性**

- `T₂` 在自由场上**成立**（Gauss，`C = 2/H_min` 有限；无质量时 `=L²/(2π²)`），而 `m = 0`（`GH1:276`）。
- ⇒ **`T₂` 这个"判别器"对 `m=0` 与 `m>0` 都给"是"** ⇒ 与 `(S)`（`GH1:190`）**同类**：`(S)` 的真值甚至不是 `m` 的函数，`T₂` 的**可证性**也不是。
- 更要紧的一层：**Gauss 上 BL/局部化是等式**（`[T]`：`V=½⟨x,Hx⟩` 时 `Var(f)=E⟨Σ∇f,∇f⟩`；极值族 `f=⟨c,∇V⟩` 取等）。数值（`python -I`，`V=x²/2`）：
```
f = x   : Var = 1.00000000 ,  E[f'²/V''] = 1        ⇒ ratio = 1.000000000000  （取等）
f = x²  : Var = 2.00000000 ,  E[f'²/V''] = 4.000000 ⇒ ratio = 0.500000        （严格）
```
⇒ 对**有质量**的 Gauss（`m>0`）它也给**零 slack** ⇒ **假阴性**：在 `(E)` 为真的模型上输出 `0`。

### 4.2 2D YM —— 判：**退化通过（零判别）**

- 2D 的 Wilson 测度在**独立 plaquette** 上分解（2D 无共享 link）⇒ `T₂`/LSI/`Poincaré` **张量**成立，常数 `= max_□ C_□`，**与 `L` 无关**。`m = 3g²/8 > 0`（`U09-M6.md:150`）。
- ⇒ 2D **通过**。但**通过的机制是"乘积测度"，4D 没有这个机制**（共享 link ⇒ 无分解）。
- 且 `4.1` 已判自由场也通过 ⇒ **2D 与自由场同输出 ⇒ 零判别力**（与 `GH1:252` 的"`R_t>1` 不是判别器（2D 与自由 `U(1)` 同输出）"同型）。
- **判**：`[M]`（退化通过）+ `[NG]`（判别器）。

### 4.3 M6 —— 判：**杀 `L`-一致 + 零模数**

- `m_L = min(ω, ½(2π/L)²)` = `0.999800 / 0.196759 / 0.012297`（`L=2/10/40`，`U09-M6.md:11`，`m_L·L² = 19.675` 恒定 ⇒ `∝ L^{-2}`）。
- ⇒ `C_T₂ ≥ 1/m_L` = `1.0002 / 5.0824 / 81.3206` ⇒ **不 `L`-一致**（与 `§3.2` 的自由场同律）。
- M6 是 **Gauss**（`ρ_a` 精确 Gaussian，`U09:11`）⇒ BL/局部化**取等** ⇒ 零 slack ⇒ **零模数**。
- 边缘盲（`U09:125–127`）：M6 与 `HO_ω` 的 `a`-边缘差 `2.259e−14` 而 `m` 一个 `→0`、一个 `=ω` ⇒ **任何边缘泛函的模数都过不了 M6**。
- **判**：`[T]+[E]`（数承 `U09`）。

### 4.4 单 plaquette —— 判：**静默**

- 单片 plaquette `SU(2)`、`θ=π`、`v=e₃`：`d²S/ds² = −(β/4)sin(s/2) < 0`（`O01:136` `[T]`）⇒ `K ⪰ 0` **假**，BL/局部化的**前件在物理点直接假**。
- 且单片的"模空间"只有 1 维角度 ⇒ needle 分解**空转**（`命题 L′` 的极端情形）。
- **判**：`[NG]`（静默）。

### 4.5 酸测签名

| 模型 | `T₂` 成立？ | 局部化成立？ | 给模数？ | 判 |
|---|---|---|---|---|
| 自由无质量场 | ✓（`C=L²/(2π²)`） | ✓（Gauss ⇒ 取等） | ✗（`C` 不 `L`-一致） | **假阳性** |
| 自由有质量场 | ✓（`C=2/H_min`） | ✓（取等） | ✗（零 slack） | **假阴性**（`m>0` 却输出 0） |
| 2D YM | ✓（张量） | ✓ | ✓（`O(1)`） | **退化通过**（机制 4D 不存在） |
| M6 | ✓ | ✓（Gauss） | ✗（`1/m_L→∞`） | **杀 `L`-一致** |
| 单 plaquette | ✗（`K⋡0`） | ✗ | — | **静默** |

**⇒ 签名 =（假阳性 / 假阴性 / 退化通过 / `L`-不一致 / 静默）⇒ 不过酸测。**

---

## §5 台账

### 5.1 检索执行核验（`AGENTS §1` 四要素）

**(1) 工具状态**：`grep` = **OK**（多批，无 `provider failure`，结果含 `spill` 文件路径）；`glob` = **OK**；`read` = **OK**（`GH1-SPECTRAL.md:160–279`、`O01-MARGINAL-BL.md:100–258`、`G4-GENERAL-BOUND.md:186–206`、`S1-CONDITIONAL-POINCARE.md:1–48`）；`pwsh` = **OK**（`Get-ChildItem`/`Get-Content` 均返回；一次 `python -I` 因无 `numpy` 报 `ModuleNotFoundError` ⇒ **改纯 `math` 重跑成功**，见 §5.1b）；`web_search` = **OK**（4 批、22 条查询，每批返回 10 源）；`web_fetch` = **OK**（HAL TEI `HTTP 200`）/ **部分失败**：`math.univ-lyon1.fr/.../logsob-poincare.pdf` 返回 `unsupported content type "application/pdf"` ⇒ 记 **TOOL-FAIL（仅该条）**，R3 因此只到书目级。

**(2) 实际检索式（逐条）**

| # | 检索式 | 范围 | 命中 | 用途 |
|---|---|---|---|---|
| 1 | `Talagrand\|T_2\|T₂\|W_2\|localization\|局部化\|Borell\|Brascamp\|Otto\|HWI\|位移凸\|displacement` | `research/v3/species/*.md` | **250 / 322**（截断，spill `...be2fd13c99f7-grep-results.txt`） | 定位仓内已有判词（承重点：`GH1c-LIT:50,158`、`O01:198–209`、`O17:42–139`、`G4:196–204`、`GI4:40,95`、`S1:22–30`、`U09:125–127`） |
| 2 | `M6` | `P04-M6NUM.md` | 1 | 确认 M6 命名 |
| 3 | `M6\|定义\|六片\|模型` | `u09-m6.md` | 30 | M6 定义（`H=(p_a²+p_b²)/2+½ω²a²`，`b` 为周长 `L` 的环，`U09:23`）与 `m_L` 值（`:11`） |
| 4 | `Otto Villani "Generalization of an inequality by Talagrand" JFA 2000 …` | web | 10 | R4 卷/页 `173(2) 361–400` |
| 5 | `Cattiaux Guillin "quadratic transportation cost inequalities" JMPA 2006` | web | 10 | R2 DOI `10.1016/j.matpur.2006.06.003`、`86(4) 341–361`（作者主页参考文献表，**书目级**） |
| 6 | `Zegarlinski log-Sobolev unbounded spin systems bounded perturbation` | web | 10 | 旁证：有界扰动族存在（未作承重） |
| 7 | `Fradelizi Guedon Pajor localization lemma log-concave needle decomposition` | web | 10 | §2.2（KLS/Klartag 前提） |
| 8 | `Talagrand T2 negative curvature modified LSI Gentil Guillin Miclo` | web | 10 | R3 书目（`PTRF 133(3) 409–436`） |
| 9 | `transportation-entropy T2 implies Poincare Gozlan Leonard survey` | web | 10 | 链方向 |
| 10 | `Talagrand T2 does not imply LSI counterexample` | web | 10 | **R2 摘要逐字**（"satisfying T2 but not LSI"） |
| 11 | `does T2 imply Poincare` | web | 10 | **R1 摘要逐字**（"equivalent to both a Poincaré inequality …"） |
| 12 | `Klartag needle decomposition CD(0,N) log-concave requires convexity` | web | 10 | Klartag `Mem. AMS 249(1180)`；`CD(0,∞) ⟺ log-concavity` 片段 |
| 13 | `Kannan Lovasz Simonovits localization lemma DOI` | web | 10 | KLS `DCG 13(3–4) 541–559`，DOI `10.1007/BF02574061` |
| 14 | `Borell "Convex set functions in d-space"` | web | 10 | Borell `PMH 6(2) 111–136`，DOI `10.1007/BF02018814` |
| 15 | `localization lemma fails non-log-concave` | web | 10 | **未查到**（见 §5.3） |
| 16 | `Brakke inequality mean curvature flow Huisken mean convex` | web | 10 | 只有 MCF 侧命中，**与本靶相关 0** |
| 17 | `Fradelizi Guedon Pajor "Localization of measures" 2004` / `Borell localization exponential needles` | web | 10 | 前提措辞 |

**(3) 命中数**：见上表；两处 `0` 见 §5.3。

**(4) 失败标记**：`TOOL-FAIL` = **1 条**（`web_fetch` 对 `application/pdf` 不支持 ⇒ R3 正文未读；`python -I` 一次因 `numpy` 缺失失败，已用纯 `math` 成功重跑，**不影响任何数**）。**"工具挂了"与"确实没有"已分开记**（§5.3 两项均为"命中 0"而非"TOOL-FAIL"）。

**§5.1b 本轮全部数值均为 `python -I` 直跑、纯 `math`**：§1.4(B)、§2.4(a)、§3.2、§4.1(E)。无第三方包、无随机数、无 MC。

### 5.2 外部引用（逐条：读到哪一层）

| 引 | 出处 | 读到的层 |
|---|---|---|
| R1 | Cattiaux–Guillin, *A criterion for Talagrand's quadratic transportation cost inequality*，**arXiv:math.PR/0312081**（HAL `hal-00133642`） | **摘要全文逐字**（HAL TEI `HTTP 200`）；正文未读；期刊卷页**未核** |
| R2 | Cattiaux–Guillin, *On quadratic transportation cost inequalities*，*JMPA* (9) **86** (2006) no.4 **341–361**，DOI `10.1016/j.matpur.2006.06.003` | **摘要**（ScienceDirect 摘要；DOI/卷页出自作者主页参考文献表 = 二手书目级）；正文付费墙，**未逐字核** |
| R3 | Gentil–Guillin–Miclo, *Modified LSI and transportation inequalities*，*PTRF* **133**(3) (2005) **409–436** | **只到 HAL 书目记录**（正文 `TOOL-FAIL`） |
| R4 | Otto–Villani, *JFA* **173**(2) (2000) **361–400** | 书目级（多源一致；DOI **未核**） |
| R5 | Holley–Stroock, *JSP* **46** (1987) **1159–1194**，DOI `10.1007/BF01011161` | 书目级（承 `O01:236`） |
| R6 | Kannan–Lovász–Simonovits, *DCG* **13**(3–4) (1995) **541–559**，DOI `10.1007/BF02574061` | 书目级 + 摘要片段逐字（needle/exponential needle/log-concave） |
| R7 | Borell, *PMH* **6**(2) (1975) **111–136**，DOI `10.1007/BF02018814` | 书目级（正文未读） |
| R8 | Klartag, *Mem. AMS* **249** (2017) no. **1180** | 书目级 + 摘要片段逐字（Ricci bounded below） |
| R9 | *Leaves decompositions in Euclidean spaces*，**arXiv:2108.07193** | 片段逐字（`CD(0,∞) ⟺ log-concavity`）；作者/卷页未核 |
| R10 | Bizeul 2023 博士论文（HAL `tel-04577476`） | 片段逐字（needle = 1D log-concave） |
| R11 | Talagrand, *Adv. Math.* **120** (1996) 109–145（高斯 `T₂(2)`） | `[T-标准]`，**本轮未开卷** |

### 5.3 未查到（四项齐，采信）

| # | 检索式 | 工具 | 命中 | 结论 |
|---|---|---|---|---|
| 1 | `localization lemma fails non-log-concave measure counterexample needle decomposition requires log-concavity hypothesis` | `web_search` | **0**（10 源全部只在 log-凹前提下陈述局部化） | **未查到**"边际凸即可局部化"的定理（检索式：上；工具：OK；命中：0） |
| 2 | `Brakke inequality mean curvature flow Huisken monotonicity mean convex` | `web_search` | 10 源，**与本靶相关 0** | **未查到**任何"Brakke ⇒ 测度侧熵-输运模数"的桥（工具：OK；命中：0） |

**注意（`AGENTS §1`）**：这两条的"0"是**命中 0**，不是 `TOOL-FAIL`；工具状态见 §5.1(1)。

### 5.4 与既有文件的关系（不重叠，只接续）

| 文件 | 它已判什么 | 本文件加什么 |
|---|---|---|
| `GH1-SPECTRAL.md:232–239` | 四类、无第五类 | **复核并给出"为什么不该有"**：`T₂ ∈ (Poincaré, LSI)`（§1.3/1.5）、局部化与墙正交（§2.3） |
| `GH1c-LIT.md:50,158` | "输运-熵 = 等价链，不产生新界" | **抓错** `⟺` → `⟹`（§1.6）；结论不变、等级降 `[R]` |
| `O01-MARGINAL-BL.md:159–209` | 切片凸 ⇏ 边缘 log-凹（`[T]`）；HS 常数 `10^{6671}`；`(CBL) ⇏ (BL-1D)` | 把它的**两句话**接到"局部化"名下：前提是联合（§2.2）、去凸性代价是第二座（§3.4） |
| `O17-OT.md:52,66,139` | `HWI` 在 `κ<0` 断在二次项；`κ<0` 的 `Var≤1/κ` 有精确反例 | 把"负曲率修正"定位到 R3 的 modified 形态，并判其在本靶 = 第二座 |
| `G4-GENERAL-BOUND.md:196–204` | `W_geo≥0` 在 `β=0` **取等**；`β>0` 不转移 `[NG]` | 用作 §4.1 的"零 slack ⇒ 零模数"的承重点（Gauss 取等） |
| `S1-CONDITIONAL-POINCARE.md:24–48`（§0 判词 2） | 具名类 `[NG]`："`U′`-一致 log-凹/BE/BL/PI/**LSI**/Prékopa" | 本文件证明"输运-熵"**不能加进这个具名类的例外**（§1.5、§2.7） |
| `U09-M6.md:11,125–127,150` | `Λ ≡ m`（Poincaré 常数 = 隙）；边缘盲引理 | **`C_T₂^♯ = 1/Λ`**（§1.4）把 `T₂` 逐字钉到这条等式上；边缘盲作 §2.5(ii) |
| `GS2-PATTERN.md:5` | "双座"模式 | `T₂` 的去凸性代价 = 第二座（广延 `10^{6671}`），§3.4 |

### 5.5 诚实边界（逐条带窗口，`AGENTS §2.1`）

1. **`C_T₂^♯ = 2λ_max(Σ) = 1/Λ` 的窗口 = Gauss 测度**（命题 T 的证只在 Gauss 上；锐性用平移族）。**非 Gauss 测度上只知 `T₂ ⟹ Poincaré`（R1），常数关系未核** ⇒ §3.1 的"同一命题"在非 Gauss 上应读作"**同一类**，不比 Poincaré 弱"，不是"逐字相等"。
2. **R1/R2 只读到摘要**。若正文的"Poincaré inequality"是**带权**版本，则 `T₂ ⟹ 经典 Poincaré` 这一步需重核；但 R2 摘要的"new families ⟺ Poincaré"与"`T₂` 而**非** LSI 的例子"两点一致地把它压在 `(Poincaré, LSI)` 之间 ⇒ 本路判词**不依赖**那一步的强弱。
3. **§2.4(a) 的反例只否掉"概率边际"读法**；"方向/切片"读法由 `O01:159–166`（**单 plaquette `SU(2)`,`θ=π`,`v=e₃`** 窗口）否掉，**不得跨窗引用**。
4. **`β` 窗口**：本文件所有关于 Savvidy 的断言都在 `λ_min = −g|B| < 0` 的**非零色磁场**窗口；`B=0`（`β=0`）时 `Hess ⪰ 0` 成立、Prékopa 取等（`G4:196–199`），本文件**不**在那一格主张任何 `[NG]`。
5. **`L` 窗口**：`m_L>0` 逐 `L` 恒成立（`U09:29` 的口径）；"无隙"是 `L→∞` 的命题。§3.2/§4.3 的"不 `L`-一致"是同一条口径的定量版。
6. **本文件没有证明任何关于 4D YM 隙的正面东西。** 它只判"这条路不通"，并把"为什么不通"落到三条可核的数上（`L²/(2π²)`、`10^{6671}`、边缘盲 `2.259e−14`）。
7. **`T₂` 是否可能在**非 Gauss **测度上严格弱于 Poincaré**：R1 摘要排除（`T₂ ⟺ Poincaré + 强 Gauss 集中`）。若正文与此有出入 ⇒ 本路判词需重开；本文件**没有**读到那段正文。

### 5.6 位置图交接（`AGENTS §6`）

本轮的**结论变化**触发 `AGENTS §6` 条件 2（"某条路线判死"）与条件 4（"某个洞关闭"）：
- 新增死路一条：**`T₂` / 局部化（本文件）**——理由三条数：`C_T₂^♯ = 1/Λ`（同命题）、`C=L²/(2π²)`（不 `L`-一致）、`10^{6671}`（第二座）。
- 新增**抓错**一条：`GH1c-LIT.md:50,158` 的 `LSI ⟺ T₂`。
- **但本代理写域限于本文件** ⇒ `docs/position-current.html` **未改**，交 `lead`。建议加一行：
  > 「第五类的两个候选（`T₂`、局部化）判 `[NG]`：`T₂ ⊂ (Poincaré, LSI)` 且 `C_T₂^♯=1/Λ`（同命题）；局部化的 1D 归约作用在**维度**坐标、墙在**曲率**坐标，两坐标正交；`GH1c-LIT` 的 `LSI⟺T₂` 应为 `⟹`。」

---

*Status：本路**判死**，`[NG]`，并给出它为什么**不该**活的理由（不是"文献没查到"）。**Q1**：`T₂` 的"弱于凸"正定理**存在**（Cattiaux–Guillin 判据 `arXiv:math.PR/0312081` 摘要逐字："`T₂` is **equivalent to both a Poincaré inequality and a strong form of the Gaussian concentration property**"；Holley–Stroock 有界扰动；张量），所以本路不是缺定理——**是缺类型**：`LSI ⟹ T₂ ⟹ Poincaré` 且**两箭头严格**（`T₂` 而**非** LSI 的例子：Cattiaux–Guillin *JMPA* **86** (2006) 341–361 摘要逐字），故 `T₂` 是墙 ② 的**内点**，而 `(E)` 只需 Poincaré 常数（`GH1:209`）⇒ **证 `T₂` 比证墙更严，不是更松**。**Q2（★）**：局部化的 1D 归约**不**只要边际凸——前提是**联合** `CD(0,∞)`（KLS `DCG 13(1995)541`、Klartag `Mem.AMS 249(2018)1180`、`arXiv:2108.07193` 片段："`CD(0,∞)` is equivalent to log-concavity"）；"边际凸"由**两个精确反例**否掉（本文件新建：`Uniform(−1,1)` 边缘 + 联合 `1+ρ(x²−⅓)(y²−⅓)`，`min f=0.888889`、边缘 `=1.000000000`、log-Hess 本征值 `+0.411777`；仓内 `O01:159–166`：切片 log-凹而边缘 `G''(0)=W₀−V₀β²<0`）；更有**命题 L′**：`G''` 是**边缘**泛函，局部化的全部收益在"拼回联合"那一步，而拼回 = Prékopa = 凸性 ⇒ **对本靶 1D 归约是空操作**，且 Savvidy 负方向**本来就是一维**（`k=0` 粗模）⇒ 可降维度 = 0。**Q3**：三条独立死因——(i) Gauss 上 `C_T₂^♯ = 2λ_max(Σ) = 2/Δ_min = 1/Λ`（两行证 + 数值 `sup W₂²/H = 1.999994`）⇒ 常数**就是逆隙**，自由场给 `C=L²/(2π²)`（`L=64: 207.5058`）⇒ 不 `L`-一致；(ii) `a`-一致 ⟺ `m·a ≥ c` = 反靶（`B04:92–97`），自由场见证 `2asinh(am/2)→0`；(iii) 去凸性代价 `e^{2osc}=10^{6671}` = 广延第二座（`O01:202`、`GS2:5`）。**Q4**：**不是第五类**；`GH1 §2.6` 的四类不动。**酸测**：自由场假阳性/有质量 Gauss 假阴性（BL 在极值族取等，ratio `1.000000000000`）/2D 退化通过（乘积机制 4D 不存在）/M6 杀 `L`-一致（`1/m_L=81.3206@L=40`）+ 边缘盲（`2.259e−14`）/单片静默 ⇒ 不过。**抓错**：`GH1c-LIT:50,158` 的 `LSI ⟺ T₂` → `⟹`，标 `[R]`。**No QED。Clay = 0%。***

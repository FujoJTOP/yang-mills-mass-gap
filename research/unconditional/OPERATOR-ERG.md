# OPERATOR-ERG — Wegner 连续酉流对角化但不造隙：`f₁²` 是酉不变量，"spectral flow" 是误称（谱不流，只流矩阵元）；Savvidy 凹性只被绕开为分母，未被绕开为控制

> **判词 (verdict).** Wegner CUT 流 `dH/dl=[η,H]`、`η=[H_d,H]` 把 off-diagonal 按 `e^{−(ΔE)²l}` 衰减：形式精确、非微扰、无背景、无圈。**代价是它每一步都酉**，故 `spec(H(l))=spec(H)` 恒成立——流的只有矩阵元，**本征值不流**（"spectral flow" 名不副实）。`0⁺⁺` 道的谱原子 `f₁²=ν({0})`（⟺ `G″(k→0)` 有正下限/有隙）是**酉不变量**，所以**精确 operator ERG 与精确对角化、与谱问题本身逐字等价**：它不给独立抓手，只把"求谱"改写成"解一个非线性无穷流方程"。**可算的（截断）operator ERG 破坏酉性**；要把它变成 `G″≥λ` 的下界，就必须有截断误差界，而该界在 Savvidy 大场区不可控——逐字是 H1。自由无质量场 `H` 已对角、`η=0`、流平凡、连续谱、`f₁²=0`：酸测**无假阳性（过）**，且这"过"由酉性结构性地保证，与"证不出 YM 隙"是同一句话。Savvidy 凹性作为 Wetterich 式**分母**被形式绕开（CUT 不选背景、不写 Hessian），作为**截断控制问题**被完整继承；且流在 `ΔE=m₁`（最小能量差）处最慢——恰在最该准的地方最慢。它**不是第 6 种物种**，而是已耗尽两范式（谱 × RG）的交集，其 no-go 正是"谱"范式被耗尽的根因。五点：**1 过（形式）/ 2 败 / 3 败 / 4 过（最干净）/ 5 败**。

## §1 流方程，与"精确即问题"
`H(l)=U(l)HU†(l)`，`η=(∂_lU)U†=−η†`，`∂_lH=[η,H]`。取 `η=[H_d,H]` 得 `∂_l(H_od)_{mn}≈−(E_m−E_n)²(H_od)_{mn}`，即 `H_od∼e^{−(ΔE)²l}`（Wegner 1994；Kehrein 2006）。band-diagonal 变体只衰减 `|E_m−E_n|>Δ` 的对，**故意保留近简并块**——而能隙 `ΔE=m₁` 正是"最近那条缝"。流收敛时 `H_d(∞)` 的对角元即本征值。**它不是不够精确，是太精确**：酉 ⇒ 谱不变。

## §2 Species no-go：`f₁²` 酉不变
`f₁²` 是 `0⁺⁺` 谱测度在原点的原子；谱测度在 `H↦UHU†` 下不变，故任何酉流逐字保持 `f₁²`。于是"从解耦读出隙" = 解出 `H_d(∞)` 的低能对角元 = 已知谱 = 已知 `f₁²`。精确流给**等式**（谱守恒），目标要**不等式**（`G″≥λ`）；酉相似不能把等式变下界——"流有隙 ⟺ H 有隙"是恒真式，不是定理。任何真正的下界内容只可能来自**额外**的变分/微扰估计（Temple/Kato），而非来自流本身。**命题（species no-go）.** 任何酉相似流都不能在**不求解谱问题本身**的前提下证明 `f₁²>0`；operator ERG 是酉相似的最纯形式。

## §3 Savvidy：绕开分母，继承控制
Wetterich FRG（E1-FRG）把负模当分母 `(Γ_k^{(2)}+R_k)⁻¹`，在 2.5 处 `Γ₂` 不定、给不出 `≥λ`。CUT 不选背景、不写 Hessian，Savvidy 凹性**不作为分母出现**——形式绕开。但 (i) `H≥0`，Savvidy 负模从来不是 `H` 的本征值而是背景 Hessian；酉流不碰它，也**不能**把它翻成正下界（翻谱即违酉）。(ii) 可算必需截断；截断破坏酉性，丢掉项的最大值恰在大场/负 Hessian 区。(iii) 流率 `(ΔE)²`：`m₁` 是最小非零差，低能态**最后**解耦。这与 G1-L3-VIA-FLOW 的"scale profile 错半"同构：把动量换成能量，把拉回 `e^{ct/a²}` 换成慢流 `e^{−(ΔE)²l}`，墙不动；E4-RMT-EDGE 是谱族的另一支，它把隙定位在反隙边，operator ERG 则连边都不碰。

## §4 唯一活路与其代价
非空洞路线只有"受控截断 + 丢弃项范数界"，该界给 `G″` 下界，而它在大场区不可控 = H1（`(★)` 的算子版）。另：保 Gauss 律要求 `H_d` 规范不变，规范不变的"对角"基就是胶球基，**对角化它** = `0⁺⁺` 谱问题本身；PROOF-GAUSS-INVARIANT-PF 已示 PF 只到真空+有限体积隙、到不了激发 `0⁺⁺`。

## Five-point
1. **Nonperturbative — PASS（形式）.** 精确流、无背景无圈；但连续极限 `H` 无界、`l→∞` 收敛需谱条件，酉流存在性未证。
2. **Sign-controlling `G″(k→0)≥λ` — FAIL.** 酉 ⇒ 谱守恒 ⇒ 只等式，无 `λ`。
3. **Non-circular — FAIL.** 解耦率 `(ΔE)²`，`ΔE=gap` 是待求量；读出即预设。
4. **Non-SUSY + fails free field — PASS（最干净）.** 自由场 `η=0`、连续谱、`f₁²=0`；酉性保证永不假造正性。
5. **Uniform `≥λ` — FAIL.** 无 `λ`、无 `V,g` 一致控制，低 `ΔE` 最慢。

**Acid test.** 自由无质量 operator-ERG 留连续谱（无隙）——**过**（无假阳性），且由 §2 酉不变性结构性保证，与"证不出 YM 隙"同源。

**Sharpest line.** Operator ERG 是全家族**最诚实**的物种：酉流从不伪造正性，故在自由场上干净失败；也正因为它是酉相似，能隙/`f₁²` 是它的**不变量**而非**产出**——"读出隙"只能靠已含隙的 `H_d(∞)`，而算出它所需的受控截断，误差界逐字是 H1。它把 Savvidy 凹性从 Wetterich 的**分母**搬成 CUT 的**慢流 + 截断误差**，同一堵墙的第三种外衣（E1 分母 → G1 拉回 → operator 慢流），不是第 6 种物种；值得银行的恰是这条"酉不变性 no-go"：凡酉者，皆换基。

*Refs.* F. Wegner, *Flow-equations for Hamiltonians*, Ann. Phys. (Leipzig) **506** (1994) 77; S. D. Głazek, K. G. Wilson, *Renormalization of Hamiltonians*, Phys. Rev. D **48** (1993) 5863; S. Kehrein, *The Flow Equation Approach to Many-Particle Systems*, Springer Tracts in Modern Physics **217** (2006); J. Kogut, L. Susskind, Phys. Rev. D **11** (1975) 395; G. Savvidy, Phys. Lett. B **71** (1977) 133. Internal: E1-FRG（Wetterich 分母、`Γ₂` 不定）；G1-L3-VIA-FLOW（scale profile 错半、`e^{ct/a²}`）；PROOF-GAUSS-INVARIANT-PF（PF 到真空+有限体积隙，不到 `0⁺⁺`）；C8-SPECIES-NOGO（species 判词范式）；E4-RMT-EDGE（谱族反隙边）；B4-VARIANCE-MAG（等分）；PROOF-SAVVIDY-WGEO-DOMINATION（负模物理）；PROOF-COARSE-POINCARE（H1、`(★)`）；G1-EXISTENCE-THEOREM（L3 ⟺ 一致算子界）；new-tools-registry（谱、RG 皆列 free-field-blind）。

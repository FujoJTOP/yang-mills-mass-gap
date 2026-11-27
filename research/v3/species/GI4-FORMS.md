# GI4-FORMS — 「恒等式形态」总普查：14 种形态、逐条判「能否作用在 `m` 或 `f₁²` 上」

> 2026-11-26 · 第 N 轮 · 靶：`m₁ = inf supp ν` 与 `f₁² = ν({m₁²})`（`docs/BRIEF-FOR-MODEL.md:15–17`）。
> 纪律：AGENTS §1 四要素（§4）、AGENTS §3 引用分级、宪章 §2 诚实铁律。标签 `[T]` 真定理 / `[C]` 条件 / `[NG]` 具名 no-go / `[O]` 精确缺失步 / `[M]` 模型或记账 / `[E]` 实测 / `[J]` 判断。**未证不标 `[T]`；无 QED；Clay 0%。**
> 本轮**不新跑 MC**、不改既有文件、不 commit。§1 的 14 行 = 题面 12 条 + 仓内增补 2 条（F13 求和规则、F14 Ward 塔，二者是仓内投入最大、题面未列的两族）。

---

## §0 · 判词

**判词 1（口径必须先改，否则整张表都会被误读）。** 本轮的判据是「能否给出 `m` 或 `f₁²` 的**精确表达式**」，即 **READ**，不是「能否逼出 `f₁²>0`」，即 **FORCE**。这不是措辞游戏：`research/unconditional/META-GAPFORCING-IDENTITY.md:11–20` 已把 **FORCE** 一侧整类关闭 `[NG]`——每条精确恒等式的解集都含无隙元（自由无质量场，公理相容、`inf supp ν=0`），故「无隙 ⇒ 恒等式被违反」这条蕴含对任何具名恒等式族都不成立（主 no-go 骨架 `:120–128` `[T]`）。**⇒ 「恒等式形态」这个方向只能在 READ 语义下有意义；按 FORCE 语义，表里 14 行全是 `[NG]`，本轮无需再普查。** 这个区分是 `:156–158` 的「READ ≠ FORCE」，也是下面所有排序的前提。`[T]`

**判词 2（上一轮的成功属于哪一类）。** `φ^U(Γ) = (−1)^{n−1}N(G_Γ)`、`N = #{NBC 生成树} = T_{G_Γ}(1,0) ∈ ℤ_{≥1}` 是一条**精确整数恒等式**，证明是删-缩归纳并等价于 Whitney 破圈定理，**只消费「`G_Γ` 是有限连通图」——与维度、`Δ=20`、plaquette 结构全无关**（`docs/PAPER-CPB-CLOSURE.md:40,43,100,209`；`research/v3/species/GO1-SIGNTHEOREM.md:14–15,260`）。它的成功机制是「**把带号交错和压成整数计数**」，可复用性是**形态级**的、不是物种级的——这正是本轮该按形态而非按物种攻的理由。`[T]`

**判词 3（可用形态总数 = 4 条；未试且可算 = 2 条）。**
- **A 类（能给出 `m` 或 `f₁²` 的表达式）= 4 条**：F3 Fredholm/谱行列式、F4 矩的生成函数（Stieltjes/Padé/Prony）、F2 删-缩/NBC、F1 Möbius 反演（Mayer–Ursell/KP）。**其中 F1、F2 在 `β=2.5` 上已被判死**（只活于 `β<1/48` 窗内），所以真正落在靶窗 `β=2.5` 上的只有 **F3、F4 两条**。`[T]/[NG]`
- **B 类（只能给界 / 循环 / 换名）= 6 条**：F5 变分（BL/Borell）、F6 积分表示（Feynman–Kac/Berezin）、F8 谱 ζ/热核、F12 概率耦合、F13 求和规则、F14 Ward 塔。
- **C 类（具名死或对 SU(2) 不可算）= 4 条**：F7 对偶/配对、F9 指标定理、F10 加法组合、F11 复分析。
- **「未试且可算」= 2 条**，且它们是**同一组 20 个数的两种读法**：
  1. **F4 在真实格点数据上从未执行**——ESPRIT/Prony 估计器已实现且在合成信号上验证（原子 0.7% 误差），但文件逐字写明「`[V-num]` on synthetic data, **not a lattice measurement**」（`research/unconditional/WI-REALTIME.md:11–12`）。`[V-num]`
  2. **F3b 迹恒等式 `Tr T_O^k = Z_{S¹_k×T³}` 作为测量通道从未被用过**——A3 已入册（`research/v3/species/D08-LOGDET.md:190` 转引 `INVENT-NEW-IDENTITY.md:5` 与 `registry:29`），但仓内一律走 correlator 侧，没人把 `D_O` 的系数当**可测配分函数**读。`[T]`（恒等式本身）/ 未试（作为通道）。**代价诚实**：需在 `k`-扭曲几何上新跑，`6³×20` 单点无此对象。
- **F3a（Hankel 铅笔直读）不算「未试」**：已解析执行且有**方向定理** `m₁^{det} := (1/a)log z₁ = inf{E_n : f_n²>0} ≥ inf{E_n} = m₁`（`research/v3/species/H08-Z1.md:103–109` `[T]`+`[E]`），合成数据上给过 `m₁^{det}(M)↓` 序列 `0.1625→0.0780→0.0554→0.051293`（`:96–101`）。它给的是 **`m₁` 的上界序列**，不是证明。`[T]`

**判词 4（今晚的教训要精确化，否则会误导下一轮）。** 「恒等式成、不等式败」这条经验**只在 READ 语义下成立**。反例就在表里：F5 变分恒等式（Brascamp–Lieb/Prékopa/Borell）**不是不等式失败**——它的失败是**前提**在 `β=2.5` 为假（需 `Hess V ⪰ 0`，而 Savvidy 精确给 `λ_min=−g|B|<0`），且即便成立也只返回方差上界（`research/unconditional/C6-DEPENDENCE-CLASS.md:3,16` `[NG]`）。同理 F8 热核与 F9 指标**都是等式**，却因「渐近/刚性」而给不出 `m`。**⇒ 等式性既不充分也不必要；真正的判据是「右端是否落在 `β=2.5, 6³×20, SU(2)` 上可测量的那 20 个数里」。** 这是本轮普查唯一的正面筛选原则。`[J]`

**判词 5（诚实上限）。** 两条 A 类活路都**不产生定理**，只产生**可复核的数**：F3a 受 `H08-Z1.md:107` 的方向定理限制（零点是 `m₁` 的**上界**，差 `m₁^{det}−m₁` 由「漏了几个态」决定，不随精度减小，`:113–119`）；F4 的收敛性条件「孤立极点 + `e₂−e₁≥Δ>0`」其 `(L,a)`-一致性**恰是 H1**（`WI-PLATEAU-RIGOR.md:36–40` `[O]/[NG]`）。故本轮**不改变 Clay 百分比：0%**。`[T]/[O]`

---

## §1 · 形态表（14 条）

### 1A 形态本体：名字 / 最小公式模板 / 它需要什么结构 / 仓内先例

| # | 形态 | 最小公式模板 | 需要什么结构 | 仓内先例（`file:line` → 判词） |
|---|---|---|---|---|
| **F1** | **Möbius 反演 / 偏序求和**（Mayer–Ursell、cluster expansion、KP） | `φ^U(Γ)=Σ_{A⊆E,(V,A)连通}(−1)^{\|A\|}`；`log Z = Σ_Γ φ^U(Γ)a^Γ` | 一个偏序格 + **一致小活动度**（收敛半径）；`Σ_{γ∋x}\|w(γ)\|e^{α\|γ\|}≤α` | `C01-CBAR1.md:13`（`φ^U`=Ursell/Mayer `[T]`）；`T1-POLYMER-RG.md:95`（判据 `[T]`）、`:23`（**不到 `β=2.5`** `[NG]`）；`V06-RADIUS.md:113`（`β=2.5,SU(2)` 要 `K<1.97`、实有 `1182` ⇒ 差 **×600** `[C]/[M]`） |
| **F2** | **删-缩 / 容斥**（NBC 生成树、Tutte/色多项式） | `φ^U(Γ)=(−1)^{n−1}N(G_Γ)`，`N=#\{NBC 生成树\}=T_G(1,0)`；`P_G(u)=P_{G−e}(u)−uP_{G/e}(u)` | **只要「`G_Γ` 是有限连通图」**（维度/`Δ`/plaquette 全不用） | `PAPER-CPB-CLOSURE.md:40,43,100,209`（`[T]` ∀n ∀维）；`GO1-SIGNTHEOREM.md:14–15,159–176,260`（五步删-缩，`[T]`）；`GR4-SPECIES7.md:78`（把它当**载体**：`β=2.5` 处 `u=1.014390>1`、根 `β*=2.4469` ⇒ `[NG]`） |
| **F3** | **行列式 / 迹恒等式**（Fredholm–Hadamard、Selberg/Gutzwiller、B6） | `f₁²=−(1/z₁)N_O(z₁)/D_O′(z₁)`，`D_O=det(1−zT_O)`，`log D_O=−Σ_k z^k Tr T_O^k/k` | 迹类转移算子 + 循环子空间 + **隔离**（读第二零点需 `S₂=o(λ₁)`） | `INVENT-NEW-IDENTITY.md:3,70–76`（`[T]` fixed `a`）；`LOC-DETERMINANT.md:9`（固定 `L` 是 **perfect probe**）、`:15`；`H08-Z1.md:17`（显式 Hankel 铅笔 `det(1−zS⁻¹A)`）、`:103–109`（方向定理 `[T]`）；`G2F-ATOM-PLUS-ISOLATION.md:3,9` |
| **F4** | **生成函数 / 连分数**（Stieltjes、Padé、Prony、Christoffel、approach-law） | `C(t)=∫e^{−st}dν(s)`；`G(z)=N/D`；`\hatν(M)` 由 `C(0..2M−1)` 的 Gauss 节点 | 精确矩 + `ν≥0` +（要**等式**须 `rank ν≤M`） | `TOOL-COMPLETELY-MONOTONE.md:51`（`r(t)` 双向界 `[T]`）；`A4-SPECTRAL-BOOTSTRAP.md:13`（有限矩**不确定** `inf supp` `[T]`）；`PADE-ROUCHE-DEEP.md:16,19,23`；`R2-ISOLATION.md:156–162,329`（有限塔 Prony `[T]`）；`WI-REALTIME.md:11–12,25`（ESPRIT 已实现，仅合成数据）；`WI-PLATEAU-RIGOR.md:22–27`（**原子/割判别律** `[T]`） |
| **F5** | **变分恒等式**（Brascamp–Lieb、Prékopa、Borell–TIS、Bakry–Émery） | `Var(f) ≤ E⟨c,Hess V⁺c⟩`（`D₊≤0`）；`W_geo≥N/12` | **`Hess V ⪰ 0`（log-凹）** | `C6-DEPENDENCE-CLASS.md:3,16`（需 `Hess⪰0`=Savvidy 墙；成立也只返 `VarG≤cW_geo`=H1，`[NG]`）；`C2-BALABAN-UNIFORM.md:173`（Borell–TIS）；`G4-GENERAL-BOUND.md:203`（`K⪰0` 失败 ⇒ 无 BL） |
| **F6** | **积分表示**（Feynman–Kac、转移矩阵半群、Berezin–Lieb） | `C(t)=⟨Ω\|O†e^{−tH}O\|Ω⟩=∫e^{−st}dν`；Berezin–Lieb `Σλ_j≤Σε_j` | OS 正性 + 正自伴 `T=e^{−aH}` | `U09-M6.md:48`（Feynman–Kac 唯一用法：`β→∞` 边缘 `=\|ψ₀\|²`）；`U11.md:18,42,54,229–231`（Berezin–Lieb=**和的上界**；全仓 `Berezin–Lieb` 用于原子的先例 **0 命中**，工具 OK）；`mass-spectrum/REPORT.md:60`（KL=谱定理作用于 `AΩ`） |
| **F7** | **对偶 / 配对**（Poincaré/Serre、Abelian duality、Wegner、Bianchi） | `Z(β)=\tilde Z(\tildeβ)`；`ker d₂=im d₁` | **阿贝尔 / 自对偶 / 2D / 拓扑面** | `HISTORICAL-TOOL-SWEEP.md:45`（FK/Wegner/渗流**只在阿贝尔**有正表示，SU(N) 带号 ⇒ 无 Gibbs 测度，`[NG]`）；`B05-DUAL.md:25`（「单方向/某边缘上的量是全空间常数的上界」；**经 `H08-Z1.md:119` 转引，本轮未直读该行**）；`BIANCHI-COAREA.md:16,74`（Bianchi=image，**零额外刚性**） |
| **F8** | **谱 ζ 与热核系数**（Seeley–DeWitt、Dirac 谱作用） | `Tr e^{−tD²}∼Σa_k t^{(k−d)/2}`；`ζ(s)=Σλ^{−s}` | 椭圆算子 + **`a→0` 渐近展开** | `CONTINUUM-FREDHOLM.md:17`；`H1-T5-NCG.md:5,15`（透射尺度 `Λ_QCD=μe^{−1/2b₀g²}` 指数在 `1/g²`，**任何有限 SD 阶看不见** `[NG]`）；`D08-LOGDET.md:59`（`log D_O` 的 `a→0` 发散是 SD–BPHZ 型）；`FORGE-DIRAC-ETA.md:11`（谱作用 `g`-连续 ⇒ 自由值 ⇒ Schur-凸，`atom≤total` 错锥） |
| **F9** | **指标定理**（Atiyah–Singer、APS、Lefschetz 不动点、spectral flow） | `ind D=∫ch`；`η∈ℝ`、`flow∈ℤ` | Dirac 型算子 + 拓扑背景 + **族（path）** | `BASE-ORDER-PARAMETER.md:45–46`（index/`η` 均 **floor-blind** `[NG]`）；`BASE-INDEX.md:22–23`（整数+同伦不变 ⇒ 刚性）；`INVENT-NEW-INDEX.md:8`（错道 + 只上界）；`HISTORICAL-TOOL-SWEEP.md:37`（**刚性冻结于自由值 0**）；Lefschetz 侧：`Q7-TOOL-FORGE.md:94`（thimble 求和权重=测度控制，`[O]`+`[NG]`） |
| **F10** | **加法组合**（Euler 变换、Borel–Padé 重求和、簇代数） | `E(f)(q)=∏_n(1−q^n)^{−a_n}`；`Σ_Γτ(G_Γ)a^Γ=Σ_T a^T∏_{e∉T}(1+a_e)` | **非零收敛半径 / 可和性结构** | `C01-CBAR1.md:100`（树恒等式 `[T]`，仓内唯一真正的加法组合恒等式）；`GU2-FIVERINGS.md:190`（Euler/Borel **`[NG]`**：`β=2.5` 处 `ū=135.9`，级数「没有定义」而非「难算」）；`C02-URSELL.md:16`（真实 4D 簇 `τ^{1/96}=16.361`） |
| **F11** | **复分析**（Rouché、Jensen、Hadamard 分解、Phragmén–Lindelöf） | `log|a(β₀)| ≤ (1/2π)∫log|a|`；零点计数 | 解析性 + **无零比较子** | **整类已判死**：`HISTORICAL-TOOL-SWEEP.md:42`（max-modulus **上界** / min-modulus **循环** / Rouché **需靶**，`[NG]` 类）；`GA1-ANALYTIC.md:91,141`（`P-CLS` 类级封口）；`MS-ANALYTIC-NOZERO.md:19`（逐件失效）；仅存 `TOOL-ZEROCOUNT.md:23,57`（零自由带 `[T]`，**读工具**） |
| **F12** | **概率耦合**（Coupling、Wasserstein、Stein、Dobrushin） | `W_p(μ,ν)≤…`；`E[∂_n g]=E[gf]`（Stein） | 两个测度 + 传输计划 + 凸序 | `BASE-FIELDMEASURE.md:40`（Wasserstein drift `=0` 或 `≡gap` ⇒ 同义反复，`[NG]`）；`H02-COUPLING.md:133`（Jensen 双向界，只给界）；`Y1-LSC-ATOM.md:111`（`W_p` 下原子**连续移动** `[T]`——真结论，但是稳定性不是表达式）；`C1-SZZ-NSCALING.md:65`（`W_2` 指数遍历，`β<1/48`） |
| **F13** | **求和规则 / 谱表示**（Källén–Lehmann、迹反常 SR、SVZ、OPE 矩）〔增补〕 | `∫₀^∞ρ(s)w(s)ds=C`；`M_{−2}∝⟨θ⟩=⟨(β/2g)¼trF²⟩∼Λ⁴` | `ρ≥0` + 一个**已知常数** | `INEQ-R10-sumrule-positive.md:3,7`（每族必落**齐次**（`f₁²` 抵消）或**聚合**（矩）；第三类 point-isolating 从未出现 `[T]`）；`A1-SUMRULE-A.md:79`（任何 KL 求和规则产不出 `a`）；`META-GAPFORCING-IDENTITY.md:44–49`（I3） |
| **F14** | **Ward / IBP 恒等式塔**（Schwinger–Dyson、Stein、BRST）〔增补〕 | `E[af]=1`；`E[H_vv]=E[f²]`；`E[a^{2j}(H−f²)]=−2j(2j−1)E[a^{2j−2}]` | 紧群上分部积分（边界项由紧性消失） | `META-GAPFORCING-IDENTITY.md:28–35`（I1 整座塔=同一句 IbP `[T]`，逐点形式 `G̃(a):=E[H−f²\|a]=−p″(a)/p(a)` 亦在此）；`BRIEF-FOR-MODEL.md:17`（`G″=1/χ` 是**逆磁化率不是隙**） |

### 1B 三判：（a）在 `β=2.5 / 6³×20 / SU(2)` 上可算？（b）等式还是不等式？（c）是否已判死？

| # | (a) 可算？〔有定义 ∧ 落在该窗的 20 个时隙/可测量上〕 | (b) 等式 / 不等式 | (c) 判死？ | 净判 |
|---|---|---|---|---|
| F1 | **△ 有定义、不可求值**：`u=0.625`、`ū=33.98`、`uK≈10⁴`；`β=2.5` 在收敛盘外（`T1-POLYMER-RG.md:23`、`V06-RADIUS.md:113`） | **等式**（Mayer 恒等式） | **是**，作 `f₁²` 载体：`[NG]`（`T1:23`）；仅在 `β<1/48` 活，且那里给出 `f₁²≥r₀>0`（`T2-TUBE-OVERLAP.md:61`） | **B′：窗内死、窗外是已用的路** |
| F2 | **✓ 可算**（整数精确枚举；`n≤6` 全量已算：`Στ = 1,20,570,18776,669215,25060608`，`ADJ7.md:13`） | **等式**（精确整数恒等式，∀n∀维） | **否**（号定理 `[T]`）；但**当载体时是** `[NG]`（`GR4-SPECIES7.md:78`，`u=1.014390>1`） | **A：可算、是等式，但不直接落在 `m/f₁²` 上** |
| F3 | **✓ 可算**：固定 `(L,a)` 上「perfect probe per box」(`LOC-DETERMINANT.md:9`)；`H08-Z1.md:17` 给显式 `det(1−zS⁻¹A)`，`S_{ij}=C(i+j)` ⇒ `6³×20` 的 `t≤19` 给 `M≤9` | **等式**（fixed `a` `[T]`；含 `f₁²` 的残数式） | **否（读）/ 是（一致化）**：`S₂=o(λ₁)⟺I2` ⇒ `[O]`（`LOC-DETERMINANT.md:9,15`）；**方向**：零点位置 **只上界** `m₁`（`H08-Z1.md:107`） | **A：最强，但只上界 + 残数需 `rank≤M`** |
| F4 | **✓ 可算**（同一组 `C(t)`；ESPRIT 已实现） | **截断等式**；真定理侧只上界（Christoffel `f₁²≤1/K_n`，`A4:13`） | **否**（`PADE-ROUCHE-DEEP.md:16` 过酸测）；一致性 `[O]/[NG]`（`WI-PLATEAU-RIGOR.md:40`） | **A：唯一能判「原子 vs 割」的读法** |
| F5 | ✓ | **不等式**（`Var≤E[H]`、`D₊≤0`） | **是**：前提在开胞失败 `[NG]`（`C6:3,16`） | **B：降级** |
| F6 | ✓ | KL **等式** / Berezin–Lieb **不等式** | KL 部分 = 定义（循环）；BL `[NG]`（`U11.md:54,246`） | **B：循环** |
| F7 | **✗ 不可算**（SU(2) 无对偶） | 等式（阿贝尔上） | **是**（`SWEEP:45`） | **C** |
| F8 | ✓（可算，甚至精确） | **渐近等式**（非精确） | **是**（`H1-T5-NCG.md:15`；`FORGE-DIRAC-ETA.md:11`） | **C** |
| F9 | ✓（整数可算） | **等式**（整数不变量） | **是**（刚性，`SWEEP:37`；`BASE-ORDER-PARAMETER.md:45`） | **C** |
| F10 | ✓（形式幂级数可算） | **形式等式** | **是**（`ū=135.9` ⇒ 无定义，`GU2-FIVERINGS.md:190`） | **C** |
| F11 | ✓（数值可算） | 混合（Jensen 公式是等式；用法只上界/循环） | **是**（整类 `[NG]`，`SWEEP:42`） | **C（已标死）** |
| F12 | ✓ | Stein **等式** / 耦合 **不等式** | 作「搬 `m`」`[NG]`（`BASE-FIELDMEASURE.md:40`）；`W_p` 连续性是 `[T]` 副产品 | **B：降级，留一条 `[T]`** |
| F13 | ✓（反常与矩可算） | **等式**（求和规则） | **是**：齐次/聚合 ⇒ 永不解出 `f₁²`（`INEQ-R10:3`） | **B：等式但错变量** |
| F14 | ✓（逐样本可测） | **等式**（IbP 塔） | 作逼隙 `[NG]`；作 `G″` 的表达式 `[T]`——但 `G″=1/χ≠m`（`BRIEF:17`） | **B：换名（H1）** |

**统计**：等式形 11/14；不等式形 3/14（F5、F6-半、F12-半）。**A 类 4（窗内 2）· B 类 6 · C 类 4。**

> **口径注（防两次混读，AGENTS §2.1）**：表里出现**两个不同的 `u`**，不可互换——F1 的 `u = β/4 = 0.625`、`ū = eΔ·u = 33.98` 是仓内 KP/`λ_Ur` 记账口径（`UR1.md:46`、`P4-OTHER.md:67`）；F2 的 `u = 1.014390` 是 `GR4-SPECIES7.md:78` 自己的 Tutte/NBC 载体归一。两者都按各自文件的窗口读，**不得**用来互相判死或互相救活。`[M]`

---

## §2 · 排序与攻击步骤（按「最可能给出 `m` 或 `f₁²` 的精确表达式」）

### ① F4 · 矩的生成函数：在 `6³×20` 上第一次真正执行 ESPRIT–Prony + approach-law 判别
**为什么第一**：唯一一条**能在 `β=2.5` 上区分「孤立原子」与「底部连续谱」的 `[T]` 判别器——孤立原子给指数逼近 `e₁+Ce^{−Δτ}`，底部有割给幂律 `e₁+1/τ`（`WI-PLATEAU-RIGOR.md:22–27` `[T]`）；而工具链已写好、只在合成数据上验过（`WI-REALTIME.md:11–12`、`:25`）。
**攻击步骤（一句话）**：用现有 `m1_4d_run`/GEVP 引擎在 `β=2.5, SU(2), 6³×20`、零动量 `A₁⁺⁺`（`A02-ROUTES.md:189,235`）上取 `C(t), t=0..19`（`N≳2×10⁴`，`GH3-OPERATOR.md:75` 估 `320–600 s/链`），对 `y_t=C(t)/C(0)` 建 Hankel `H[i,j]=y_{i+j}` 做 SVD 与 shift-invariance 得 `z_j`、`E_j=−ln z_j/a`、Vandermonde 解 `c_j`，然后**只报三个数**：谱秩 `r` 随 (窗口,`N`) 的稳定性、`m₁^{ESPRIT}`、`f₁²=c_{argmin}`，并用**逼近律的形状**（指数 vs `1/τ`）而不是数值大小判原子。

### ② F3 · Fredholm 行列式：Hankel 铅笔的残数直读（`f₁²=−λ₁Res`）
**为什么第二**：它是唯一**同时供给** `m₁` 与 `f₁²` 的精确对象（`META-GAPFORCING-IDENTITY.md:197` 逐字：不消费、而供给 `m₁²`），且固定 `L` 上无近似（`LOC-DETERMINANT.md:9`）。
**攻击步骤（一句话）**：同一组 `C(t)` 构造 `S_{ij}=C(i+j)`、`A_{ij}=C(i+j+1)`，`M=1..9` 逐阶算 `eig(S⁻¹A)` 得 `z₁^{(M)}` 与残数 `−(1/z₁^{(M)})·N_O(z₁^{(M)})/D_O′(z₁^{(M)})`；**先把方向写死**——`H08-Z1.md:107` 的 `m₁^{det}≥m₁` 是 `[T]`，所以产出是**上界序列**加一个候选 `f₁²`，只有当 `M=3..9` 的 `z₁` 与 ①的 `z₁` 重合、且 `m₁^{det}` 收敛时，才允许把它读成 `f₁²` 的数值（否则 `rank>M`，读的是 Gauss 节点而非原子，`H08-Z1.md:101,119`）。

### ③ F2+F1 · 组合计数：把 NBC 整数计数当**振幅的闭式**，不再当符号
**为什么第三**：这是上一轮唯一成功形态的直接延伸，`[T]` 且与维度无关（`PAPER-CPB-CLOSURE.md:43`）；它给的是**精确整数**，不是估计。
**攻击步骤（一句话）**：在 `β=2.5, SU(2)` 取 4 壁管族 `Γ_t`（`|P|≥4t/a`），逐 `Γ` 用删-缩算 `N(G_Γ)=T_{G_Γ}(1,0)`（整数精确；`n≤6` 的 `Στ` 已在 `ADJ7.md:13` 全量复现），把 `φ^U=(−1)^{n−1}N` 代入 `Σ_Γφ^U u^{|Γ|}` 得振幅 `c` 的**整数闭式**，再把它在 `u=0.625` 处求值——预期诚实产出是「`c` 的精确组合表达式 + 其收敛半径」，**不是** `f₁²` 的号（判死线已给：`u_*=1.014390>1`、`β*=2.4469`，`GR4-SPECIES7.md:78`）。

---

## §3 · 不可能清单：结构上给不出 `m` 的形态

1. **元层面（最重要的负结论，覆盖全表）**：任何精确恒等式的解集都含无隙元——自由无质量场公理相容、`inf supp ν=0`、`f₁²=0`，故落在每条恒等式的解集内（`META-GAPFORCING-IDENTITY.md:11–20`；主 no-go `[T]` `:120–128`；无隙见证表 `:98–109`）。**⇒ 恒等式只能 READ，不能 FORCE。** 这条不是某个形态的死因，而是**所有形态的天花板**。`[T]`
2. **F5 变分（BL/Prékopa/Borell/Bakry–Émery）**：要 `Hess V ⪰ 0`，而 Savvidy 精确给 `λ_min = −g|B| < 0`；即便前提成立也只返回 `VarG≤cW_geo` = H1（`C6-DEPENDENCE-CLASS.md:3,16`；`G4-GENERAL-BOUND.md:203`）。**只对上界说话。** `[NG]`
3. **F13 求和规则（KL/迹反常/SVZ/OPE）**：每族要么**齐次**（两端同乘 `f₁²`，恒等抵消，`META:63–68`），要么**聚合**（矩/总量，只上界）；第三类「点隔离 coercivity `total≤C·atom`」在文献与仓内**从未出现**（`INEQ-R10-sumrule-positive.md:3,7`）。**等式，但作用于错的变量。** `[T]`
4. **F14 Ward/IBP 塔**：`E[af]=1` 与 `G″·Var(a)=1` 的右端是 `G″ = 1/χ`——**逆磁化率，不是隙**（`BRIEF-FOR-MODEL.md:17`）；塔对每个测度每个 `β` 成立（`META:28–35`），故对自由场同等成立。**换名。** `[T]`
5. **F9 指标 / F8 谱 ζ**：整数不变量在连通的 `β`-族上**局部常值 ⇒ 冻结于自由值 0**（`SWEEP:37`；`BASE-ORDER-PARAMETER.md:45`）；谱 ζ 是 `a→0` 渐近，透射尺度 `Λ_QCD∝e^{−1/2b₀g²}` 指数在 `1/g²`，任何有限 Seeley–DeWitt 阶看不见（`H1-T5-NCG.md:15`）。**一个是刚性，一个是渐近。** `[NG]`
6. **F7 对偶/阿贝尔化**：SU(N) 的曲面展开带号（SU(2) 按 Euler 示性数、SU(N) 约半数负）⇒ 无 Gibbs 测度、无渗流判据可陈述；FK/Wegner 只对阿贝尔/2D/自对偶成立（`SWEEP:45`）。**只在阿贝尔情形成立。** `[NG]`
7. **F11 复分析整类**：max-modulus 只上界、min-modulus 需「已无零」（循环）、Rouché 需无零比较子（=靶）（`SWEEP:42`），类级判词已在 `GA1-ANALYTIC.md:91,141` 收口。**注意：这是题面点名要标清的一条，本轮确认仍死，无翻案。** `[NG]`
8. **F6 积分表示（除 KL 定义外）**：Feynman–Kac 在仓内唯一用法是「`β→∞` 边缘 `=|ψ₀|²`」（`U09-M6.md:48`），它给的是**基态占据**不是谱边；Berezin–Lieb 给「和的上界」（`U11.md:42`），且其假设在 YM 上就是靶。**循环或反向。** `[NG]`
9. **F10 加法组合（Euler/Borel–Padé 重求和）**：`β=2.5` 处 `ū=135.9`，级数**没有定义**而不是「难算」；Euler/Borel 需要一个非零收敛半径作前提，而缺的正是它（`GU2-FIVERINGS.md:190`）。**只在收敛盘内是等式。** `[NG]`
10. **通用陷阱（不属任何单一形态，但会伪装成"新形态"）**：凡是输出为**线性/凸泛函**的（矩、迹、熵、容量、`∫φdν`），由线性数据引理对任意有限固定线性数据都存在「原子在任意处」的 `ν'`（`META:23` 引 `LEVER-nonconvex-edge §0B` `[T]`）。**⇒ 新形态必须先自证"非线性"才有资格进 A 类。**

---

## §4 · 台账（AGENTS §1 四要素）

### 4.1 工具状态 / 检索式 / 命中数 / 失败标记

| # | 工具 | 状态 | 检索式 | 命中 |
|---|---|---|---|---|
| 1 | `pwsh` (`Get-ChildItem`) | **OK**（有效响应，含截断落盘 `dsh-spill-…\34b9742b089e-pwsh.txt`） | 列 `D:\Dev\math` 与 `research/` 树 | 目录级全量 |
| 2 | `glob` | **OK** | `**/AGENTS.md`；`research/**/{LOC-DETERMINANT,ABC2-DUALITY,PROOF-P3-UNCONDITIONAL,TOOL-FACTORIZATION-ALGEBRA,SUPPORT-NOT-ATOM,B05-DUAL,A1-SUMRULE-A,KEY-EXISTS-OR-IMPOSSIBLE}.md`；`research/v3/species/GI*.md` | 2 / 8 / **0**（`GI*` 0 命中 ⇒ 本文件是该系列首篇） |
| 3 | `grep` | **OK** | `varphi\^U\|phi\^U` @`*.md` | 14 |
| 4 | `grep` | **OK** | `Rouch\|Jensen\|Hadamard\|儒歇\|阿达马` @`research/**/*.md` | 162（显示 250/总量截断，全表落 spill） |
| 5 | `grep` | **OK** | `Weitzenb\|Brascamp\|Borell\|Berezin\|Feynman\|Atiyah\|Lefschetz\|Wasserstein\|Euler 变换\|Möbius\|Koteck` @`research/**/*.md` | **250/834**（截断，全表落 `e62e7a3e826f-grep-results.txt`） |
| 6 | `grep` | **OK** | `Berezin\|Feynman.Kac\|Seeley\|Poincaré 对偶\|Serre 对偶\|Euler 变换\|簇代数\|Källén\|求和规则\|反演` @`research/**/*.md` | 250/413（截断） |
| 7 | `grep` | **OK** | `6³\|6\^3\|×20\|24\^3\|16\^3` @`research/**/*.md` | 35（`GH3-OPERATOR.md:75` 给出 `6³×20` 的成本 320–600 s/链） |
| 8 | `grep` | **OK** | `生成树\|spanning tree\|删-缩\|deletion-contraction\|破圈\|broken circuit\|Kirchhoff` @全仓 `*.md` | **121** |
| 9 | `grep` | **OK** | `Feynman.Kac\|Berezin` @全仓 `*.md` | **25**（其中 `Feynman–Kac` 本名仅 `U09-M6.md:48` 一处；`Berezin–Lieb` 见 `U11.md`） |
| 10 | `grep` | **OK** | `Prony\|ESPRIT\|matrix.pencil\|矩阵束\|Hankel.*本征\|广义本征` @`research/**/*.md` | **32**（关键：`WI-REALTIME.md:11–12` 明写 "not a lattice measurement"） |
| 11 | `read` | **OK** | `META-GAPFORCING-IDENTITY.md`（全 233 行）、`HISTORICAL-TOOL-SWEEP.md`（全 113 行）、`BRIEF-FOR-MODEL.md`（前 66 行）、`LOC-DETERMINANT.md:3–15`、`WI-PLATEAU-RIGOR.md:20–47`、`WI-REALTIME.md:1–40`、`H08-Z1.md:88–127` | — |

**失败标记**：**本轮无 `TOOL-FAIL`**。所有 `grep` 返回的均为有效响应（非 `ENOENT` / 非 `ripgrep provider failure`），截断者已注明总数并给出落盘路径。
**明确"0 命中"的一条**：`glob research/v3/species/GI*.md` → **0**（检索式如上；工具 OK；⇒ 本文件为 `GI4` 首篇，非重贴）。

### 4.2 本轮结论的等级分布

| 等级 | 条目 |
|---|---|
| `[T]` | F2 号定理（`PAPER-CPB-CLOSURE.md:40,209`）；F3 fixed-`a` 残数式（`INVENT-NEW-IDENTITY.md:74`）与**方向定理** `m₁^{det}≥m₁`（`H08-Z1.md:107`）；F4 `r(t)` 双向界（`TOOL-COMPLETELY-MONOTONE.md:51`）与原子/割逼近律（`WI-PLATEAU-RIGOR.md:22–27`）；F10 树恒等式（`C01-CBAR1.md:100`）；F14 IbP 塔（`META:28–35`）；F13 求和规则三分类（`INEQ-R10:3`）；恒等式-逼隙主 no-go（`META:120–128`） |
| `[C]` | F1/F2 的窗内强度：`β<1/48` 的 `f₁²≥r₀>0`（`T2-TUBE-OVERLAP.md:61`）；`H1A-EXACT.md` 的 `β*` 窗口读数 |
| `[NG]` | F1（载体）、F5、F7、F8、F9、F10、F11、F12（搬 `m`）、F13、F14（逼隙） |
| `[O]` | F3 的 `(L,a)`-一致化（`S₂=o(λ₁)⟺I2`）；F4 的 `2M=O(1/(am₁))` 一致矩（`PADE-ROUCHE-DEEP.md:19`） |
| `[M]` | `V06-RADIUS.md:113` 的 `K<1.97` vs `1182`（×600）记账；`GR4-SPECIES7.md:78` 的 `u=1.014390`、`β*=2.4469` |
| `[E]` | `H08-Z1.md:96–101` 的 `m₁^{det}(M)↓` 序列；`ADJ7.md:13` 的 `Στ` 全枚举；`WI-REALTIME.md` 的合成数据误差 |
| `[J]` | 本文件 §0 判词 3/4 的口径选择与排序 |

### 4.3 局限（必写）

1. **`β=2.5, SU(2)` 的 `T=0` 胶球质量在仓内"未取数"**（`A4-SCALE-SEP.md:73`；`GN4-REDTEAM.md:77`）——语料里 `β=2.5` 的数都是有限温屏蔽质量。故 §1B 的「可算」一律指**有定义、可算**，**不等于已算**。`[M]`
2. `β=2.5` 只在「SU(2) + `β_W=4/g²`」这一格是交叉区入口（`ξ/a≈1.0–1.3`），在 SU(3)/`β_A` 格里不是（`A4-SCALE-SEP.md:59`；`BETA-CONVENTION.md:107,194`）。**本表所有 `[M]/[E]` 行都按该格读。** `[T]`
3. §2 的三步都**不产生定理**，只产生可复核的数；产出若与 `H08-Z1.md:107` 的方向冲突，必须先怀疑「漏态」而非怀疑定理。`[J]`
4. 本轮**未开外部文献**：正文出现的 Whitney / Kotecký–Preiss / Brascamp–Lieb / Berezin–Lieb / Seeley–DeWitt / Atiyah–Singer 均为**给形态命名**用，其权威判词全部承自仓内 `file:line`（`HISTORICAL-TOOL-SWEEP.md`、`META-GAPFORCING-IDENTITY.md`、`U11.md` 等）。**未编造任何 arXiv/DOI。**

---

*状态：形态普查（14 条）+ 三判 + 三条排序 + 不可能清单 10 条。可用 4 条（窗内 2 条），未试且可算 **2 条**（F4 真实数据首跑、F3b 迹恒等式测量通道）。不证 Clay；不改既有文件；无 MC；**Clay 0%**。*

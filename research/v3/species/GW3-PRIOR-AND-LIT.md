# GW3-PRIOR-AND-LIT — 第 216 轮：DeepSeek"下界机器"提案的仓内前科（Task A）+ 三条路线的文献侦察（Task B）

> **第 216 轮 · 2026-11-25 · 工作目录 `D:\Dev\math`。** 任务：(A) 检查仓内 40+ 轮判决对 DeepSeek 提案的覆盖（(F) 微分不等式 + 桥引理 + 路 (a) RP/slab + 路 (b) 相对熵 + 路 (c) 窗内标定 + "RG 机输出比值/流而非地板"）；(B) 对三条路线做 arXiv/文献侦察。无 MC、无 QED、**Clay 0%**。
> **标签**：`[T]` 真定理 · `[C]` 条件 · `[NG]` 具名类 no-go · `[M]` 混合 · `[E]` 数值/算术 · `[J]` 判断 · `[R]` 外部已入账事实 · `[O]` 缺失步。
> **纪律（AGENTS §1）**：所有"未查到"结论附四要素（工具状态/检索式/命中数/失败标记）。本文件只写自己，未 commit。

---

## §0 · 判决表（先给答案：提案元素 → 仓内状态）

| # | 提案元素 | 仓内状态 | 一个数 / 一处出处 |
|---|---|---|---|
| 1 | **(F) 原理，残数载体** `F=f₁²`：`∃C(β)<∞, (L,a)-一致, 可积: d ln F/dβ ≥ −C` | **同一尸体**。D04 命题 2.1：`∃C<∞: \|F'\|≤C\|F\|` **⟺** `F>0` **⟺** `log F` Lipschitz（假设逐字等价于结论）；路线判 `[NG]` | D04-ODE.md:14（命题 2.1 `[T]`，判 (a) `[NG]`） |
| 2 | **桥引理**（(F)+种子 ⇒ 地板，Grönwall） | **已 [T]-证但 vacuous**：Grönwall 免费、常数不免费，`C<∞` 就是目标；逐 `(L,a)` 版早已 `[T]` = `V03-ANALYTIC:174` 的 vacuous | D04-ODE.md:103（推论 2.2）、:205–207（§4.4） |
| 3 | **(F) 标度版** `d ln(m/Λ)/dβ` 受 `a(β)^k` 控制 | **半新（新坐标、同一尸体）**。斜率恒等式已建：`d ln(ma)/dβ = −(12π²/(11N²))[1−(b₁/b₀)g²+O(g⁴)]`（纯运动学重言式 (A)），且**斜率等式预设 `ma>0`** ⇒ 结构上不能证隙 `[NG]`；"`R` 的 `Λa`-无关性 = 标度假设 = C1" 已钉 | C1-SLOPE-EQUIV.md:11,:21,:186；C3-UNIQUENESS.md:152 |
| 4 | **(F) 质量载体** `m=E₁−E₀`：HF 一阶导**无 `1/Δ`** | **事实已记、路线同死**。`dΔ/dλ = ⟨G\|∂_λH\|G⟩−⟨0\|∂_λH\|0⟩`（HF，无隙分母）已逐字写在 DF 的 C2 形；但 (F)⟺`m>0` 的等价性对任何 C¹ 载体成立（D04 命题 2.1 不挑载体），且 Lipschitz 常数硬上限 `C ~ osc(S)=12V` 对质量载体同样成立 ⇒ `L⁴` 恶化照抄 | DF-CONTINUITY.md:193；D04-ODE.md:68（命题 B(i)）、:198（§4.3 `C≤12V`） |
| 5 | **路 (a)** 有限差分 + RP slab/chessboard 局部乘性界 | **两半都已判**。①"有限差分版"= H04 的积分 Grönwall：精确版恒等式（零内容），估计版= `d04` 的 `(G)`，且反向三角否决 `∫\|F'\|≥\|F(2.5)−F(β₁)\|` ⇒ 判据要求 `F(2.5)<2F(β₁)`；②"RP 局部界"：chessboard 只给上界、不给聚类/隙/RG 步进（os3 O4）；arXiv `chessboard`+`mass gap` 检索 **0 命中**（工具 OK） | H04-THERMO.md:1,:20,:64–78；os3 REPORT.md:493,:526(O5)；§2 检索表 |
| 6 | **(G) 型判据**（(F) 的积分货币）`∫A/Δ<∞` | **已判 `[NG]`**（D07）：`∫A ≥ 2.76√V` ⇒ `(G)` 要求 `Δ_max ≥ 0.276L²`；**统一隙不蕴含 (G)**（`L`-一致量词上方向反）；2D YM `𝔊=0.885L→∞`、自由场 `0.487L³→∞`；账本在案 | D07-G.md:86,:91–98,:159–168；position-current.html:359 |
| 7 | **log-det 路（P2）** | **已判 `[NG]`**（D08）：`‖(1−K)^{-1}‖ = 1/(1−z/z₁)` 是精确等式；`z→1` 时 `(G_det)` **逐字退化为** `(G)`（(5.2)）；结论类型 `I1→I2`；酸测更差（有隙无原子的有质量自由标量上照样通过）。**不是最近活亲——是 `(G)` 的坐标变体** | D08-LOGDET.md:14,:19,:162–180,:220 |
| 8 | **路 (b)** 相对熵/熵率/Pinsker/LSI | **三次判决须分层**：①具体二阶方差率 `Var(S_W)/V` 在 2D YM 上正而残数恒零，**无条件非法**；②DV 本身是连续泛函的 sup ⇒ lsc，但其 β-导数不自动 lsc；强制引理只封**已另证 lsc 且逐点 `Ψ≤f₁²`** 的 minorant，不能据此封所有熵/LSI 泛函；③信息泛函三死类与 LDP 恒等式另记。**正面资产**：SZZ 均匀 LSI 在 `β<1/48` 真给隙（关联意义）= master §2.5 Lemma SC | GW5-ACID.md §1.3–§2、§5.4；INEQ-R5-donsker-varadhan.md:3,:9；FALSIFY-ENTROPY.md:5–11；YANG-MILLS-MASTER.md:59–60 |
| 9 | **路 (c)** 窗内标定（逐项 `d/dβ` 团簇展开） | **作为任务：新**（无人对 Z01 展开逐项求过 `d/dβ` 的显式 `C(β)`）。作为路线：**构造性受限**——KP 窗 `β<3.53×(1/48)` 到靶 `β_W=2.5` 差 `~34×`（`2.5/0.0736`），且团簇展开类文献的"对参数求导"输出是压强/关联长度解析性，不是隙地板 | YANG-MILLS-MASTER.md:42,:128–130；§4 文献 |
| 10 | **支持性声明**"RG 机输出比值/流，不输出地板——问对问题是比值/流界" | **按字面=假，按精神=半真**。反例：Göpfert–Mack CMP 82 聚合物机**在 3D abelian 输出过隙地板**（全耦合）；更强的反证：仓内 MS-RG-ODE **已经有**残数的比值/流闭式界 `f₁²(2.5)≥c·r₀`（`c=[ln 比]²>1`）而它**搬不动质量**——`m₁` 是 RG 不变量 `d ln m₁/d ln μ≡0`。半真部分：Balaban 机输出确为 UV 稳定性（`\|ρ\|` 型绝对值界、二阶导无符号、无 IR、无体积一致） | MS-RG-ODE.md:9,:19–22,:28–30；V04-RG.md:17–20,:76–77；DF-CONTINUITY.md:131 |
| 11 | DeepSeek 先验 (a) 15–20% / (b) 10% | **对照**：仓内最高单路先验 = TWENTY-PATHS A1"多尺度 RG+质量" **8–12%**；路 (a) 15–20% 高于仓内任何已判路径，过报 `[J]`（§6） | TWENTY-PATHS.md:24,:124 |

**三问直答**：
1. **(F)（残数载体）是不是 D04 命题 2.1 的同一尸体？——是**（预期成立）。连"零点处发散律 `2/|β−β*|` ⇒ `C=∞`"都已 `[T]`（D04 命题 2.3，:15，用 V02 定理 Z）。
2. **(F)（质量载体）是不是新的？——事实半新、路线不新。** "HF 一阶导无 `1/Δ`"已逐字在 DF-CONTINUITY:193（`dΔ/dλ` 的 HF 形）与 D04:68（对角情形）；(F) 的等价性对任意 C¹ 载体同一（D04 命题 2.1），`C~V` 硬上限同活（D04 §4.3）；有限差分版= H04 已判的积分 Grönwall。**仓内没有文件把"(F) 用在 `m=E₁−E₀` 上"写成完整提案**——这是唯一的字面空白，但它的每一块组件都已带判决。
3. **log-det（P2）是不是最近活亲？——不是。** D08 已判 `[NG]` 并把 D04:258 的 P2 从 `[C]` 降 `[NG]`；(G_det) 在 `z→1` 逐字退化为 (G)。**(F) 家族在仓内无活亲**：master §5.1 曾记"(G) 唯一还活着（是否 vacuous 未定）"，D07/D08 两轮后该行已过时——position-current §06 账本已记 `(G)` 死（发散、`Δ_max≥0.276L²`）。

---

## §1 · Task A：仓内前科逐文件（file:line）

### §1.1 D04-ODE.md（关键文件；提案 = 它的改写）
- **命题 2.1**（:87–103，`[T]`）：`F∈C¹`、`F(β₁)>0` ⇒ (i) `∃C<∞: |F'|≤C|F|` ⟺ (ii) `F>0` ⟺ (iii) `log F` Lipschitz。**这就是 (F) 原理的逐字否定：它不是更弱的靶，是同一句话的改写。** 推论 2.2（:103）：Grönwall `F(2.5) ≥ F(β₁)e^{−CΔβ}` 真 `[T]` 但不免费——"只需 `C<∞`"= 同义反复。
- **命题 2.3**（:105–111，`[T]`）：零点处 `|F'/F| = 2/|β−β*| + O(1)` ⇒ `C(β*)=∞`，发散律唯一（奇数阶零点被 `F≥0`+解析排除）。⇒ DeepSeek 提案里"局部-in-β 一致性比全局收敛弱"的窗口，在零点处**精确地**失效。
- **判 (c) 双层**（:74–81）：`dF/dβ` 的 Kato 谱表示**书写前提**就是 `Δ(β)>0 ∀β∈[β₁,2.5]`（`IMP-ARROWS:377` 引 Kato）；定量上任何可证 `C` 都带 `1/Δ`（(1.4)）。⇒ 路 (a) 若靠谱表示写斜率，入场费=目标。
- **命题 B 二分**（:66–70，`[T]`）：对角情形无 `1/Δ` 但平凡（`F≥0` 免费）；非对角（残数类）`1/Δ` 结构性不可消。**质量载体的 HF 导数属"对角侧"**——无 `1/Δ`，但这不救 (F)，见判 4。
- **(G)**（:126–161）：命题 3.1 `(G) ⇒ F>0` 逻辑不循环 `[T]`；但 §3.4 三能级反例 ⇒ `(G)` 严格强于目标；§3.5 `(G)` 的验证= 六条全死的同一靶。
- **L 一致性账**（:165–207）：需求 `C ≤ log(r₀/ε)/Δβ`（`ε=1e−6→3.868`）；可证 `C=12L⁴/Δ`（`L=8→9.8e4`，差 `2.5e4×`）；乘性单调/Laplace 版 `C ≤ 12V`（地板 `e^{−29.96L⁴}`，`L=4: 10^{−3331}`）。**`C ~ V` 是 β-方向 Lipschitz 常数的硬上限（`∂_βH = S_W ~ V`）——对质量载体同样成立。**
- **P2 log-det**（:256–258，`[C]`）：**已被 D08 执行并降 `[NG]`**（判 7）。

### §1.2 DF-CONTINUITY.md——"没有中间地带"已把微分不等式路线的位置分类完
- §4.2 三形（:176–211）：**C1 拓扑形**（闭性）弱于目标且不可达；**C2 分析形=单侧 Lipschitz/Grönwall**（:193 逐字给出 `dΔ/dλ = ⟨G|∂_λH|G⟩−⟨0|∂_λH|0⟩`，HF 形、无隙分母）——**等价或更强**（总变差在无界路径上一致）；**C3 单调量不存在**（`κ` 轴无单调 `[T]`：正算子加项 22.6% 概率降隙）。判词（:211）：**"能便宜陈述的那条弱于目标且不可达；能支撑结论的那条等价于或强于目标。没有中间地带。"** ⇐ DeepSeek 的 (F) 正落在 C2 槽：它"能便宜陈述"（对质量载体），所以按 DF 的分类它**弱于目标且不可达**；要"能支撑结论"就回到 `(G)`/一致 Lipschitz = 等价。
- (S⁺) 反例（:186，`[T]`）：`Δ=λ`、`Δ=e^{−1/λ}` 满足 (S)/(S⁺) 全部前提而 `Δ(0)=0` ⇒ **任何 C^k/导数假设都不能把"连续性"升到"正性"**。这是对 (F) 思路最短的一般性否证。

### §1.3 C1/C2/C3 斜率系列——"数值先验"与 (F) 标度版已在仓
- C1-SLOPE-EQUIV.md：恒等式 `1/(2b₀g²)=12π²β/(11N²)` 逐位核过（SU(2)=3π²/11=2.691710、SU(3)=4π²/33=1.196316，:11）；两圈修正 `−b₁/b₀ = −17N/(88π²)`（:19）；**判 (A)**：斜率等式是纯运动学重言式（:13,:99–102）；**循环性**（:21,:186）："`ln(ma)` 要求 `ma>0` ⇒ 斜率等式预设质量隙，不能用来证质量隙 `[NG]`"。⇒ **(F) 标度版 = C1 的 (A)+(B) 已全覆盖；DeepSeek 的 `a(β)^k` 形状比 C1 的两圈 RG 形状更弱（无 `k` 的显式值），仓内无该字面式**（Q11 检索：`a\^k|cutoff-suppressed` 0 命中，见 §1.6）。
- C2-MEASURED-SLOPE.md：**"数值先验"步已存在**——SU(3) `β=6.00→6.30`（`L=6→24`，混体积）实测 `d ln(ma)/dβ = −0.638±0.658` vs 要求 `−1.197`，比值 0.53（:8,:52）；SU(2) `tm_gap` 语料不能检验渐近斜率（`N_τ=4` 两点在解禁闭上方，符号反；:54–60）；`β₉₀≈27±2`（O(10) 稳健，:96–98）。**限定窗**：SU(3) Wilson `β_W=6/g²` 口径、`Δβ=0.3`、误差条判别力 0.69σ（C1 §6.2）。
- C3-UNIQUENESS.md：`d ln m/d ln a|_β = −1` 纯运动学 `[T]`（:147）；C1 斜率等式 = C3 唯一性的定量形式，**同源不可并行派工**（:152）。

### §1.4 V02-BETAMONO.md——零点结构已 `[T]`，(F) 的"无零点"假设就是靶
- **定理 Z**（:104–115，`[T] mod [KL]`）：`M(β)=⟨Ω|O|𝔾(β)⟩` 实解析（`T(β)` 元素为 β 的整函数 + PF）⇒ 零点集离散、无内点、测度零；例外点 `f₁² ≈ c(β−β*)²` 二次触零、`d ln f₁²/dβ` 单极点。⇒ **(F) 在残数载体上唯一可能的破口 = 孤立点零点 = 搬运所需的正是那一点**（:115 逐字"定理 Z 之后剩下的全部内容就是那一次升级"= DF:247 判的 Clay 级）。
- 单调性候选全死（:56–68）：FKG（矩语言）、RP（t 方向非 β 方向、`RP≠FKG`）、Griffiths（矩不决定原子，:83–88 精确对 `ν₁=δ₁` vs `ν₂`）。⇒ 路 (b)/(a) 若要"斜率有界 ⇒ 无零点"，前置的"残数不增/不降"无定理（B04:167 同判：转移矩阵本征值的 β-单调性"文献与本仓都没有这个工具"）。
- 2D 精确反例（:14）：**有隙而 `f₁²≡0`** ⇒ 对质量载体 (F) 有隙不保残数，对残数载体 (F) 的酸测上 `F≡0` 使不等式 `0≤0` 恒真 ⇒ **工具 gap 含量 0**（D04 §5，:211–215）。

### §1.5 B04 / G05A / B01 / MS-RG-ODE——Λ_lat 公式与"流搬不动质量"
- B04.md：`Λ_L a(β)` 两圈严格减（131 点、`max Δ=−2.4e−115`，:58–63）；`m·a=(m/Λ)(Λa)` 严格减 →0 ⇒ **格点单位地板=平凡性反靶**（:92–98）；`m/Λ` 平=标度假设（:33）；**`m·a` 的严格单调无定理**（:167）。
- G05A.md：`aΛ(β_A=1/48)=0.570`，**严格 N-无关**（`b₀g²=11N/(48π²)·(96/N)=22/π²` 与 N 无关，:113–119）；两圈前因子符号裁决 `(b₀g²)^{−51/121}`（:20–31）。窗口换算：D.3（master:136）"`β<0.0736` 是粗格点（`aΛ=0.570`、`m₁a=18.26` 纯数 ⇒ `a` 有正下界），不是 `a→0`"。
- B01.md：`w₁` 平台 `≈0.11`（`β=10`，`9.5σ`，:33,:61–64）——**UV 限定**：固定 `L_s=8`、`a(β)→0` ⇒ 物理体积收缩，是"一步切片关联仍正"的短距陈述，不是物理单位隙下界（:68–73）。窗：`β∈[3,10]`、SU(2)、`m=1.546` 固定。
- MS-RG-ODE.md：`d ln f₁²/d ln μ = +2γ_O`（`[T]`，`γ_O>0` 非 Kato 扰动级而由 AF 结构保证，:9–12）；闭式 `f₁²(μ₂)/f₁²(μ₁)=[ln(μ₂/Λ)/ln(μ₁/Λ)]²>1`，`c` 无 `V/L/a` ⇒ 继承 P3 的 L-一致性（:19–22）；**对抗性击杀**（:28–30）：`m₁` RG 不变，`d ln m₁/d ln μ≡0`，"流携带权重，不携带标度"；缺口臂 = 绝对标度 `m₁/Λ` = `I2/H1`。**这是对 DeepSeek 支持性声明的仓内精确回答：比值/流界不仅"问对问题"，而且已经做出——它的产出是把问题透明化为 `m₁>Λ_QCD` 同一问，而非解决。** 且此 ODE 是重整化理论陈述，格点侧匹配 `ν_b` 到 Wilson 的缺件是带号 `(M1)`（:31）。

### §1.6 仓内 grep 账（AGENTS §1 四要素；工具全 OK，命中逐条）

| # | 检索式（逐字） | 工具 | 命中 | 关键定位 |
|---|---|---|---|---|
| Q1 | `slope`（`research/v3/species/*.md`，-i） | Bash grep **OK** | **29** 行 / 15 文件 | C1-SLOPE-EQUIV、C2-MEASURED-SLOPE（4）、RG5-NOGO、Y02-CHANNELS 等 |
| Q2 | `relative entropy\|熵率\|Pinsker\|log-?sobolev\|LSI`（`research docs`，`*.md`） | grep **OK** | **1101** | 顶：R2-H2GLOBAL 43、LOOP-RUNBOOK 38、PL-POLCHINSKI 28、INEQ-R5、FALSIFY-ENTROPY、k0-LARGE-DEVIATION、S1-CONDITIONAL-POINCARE（SZZ/LSI 体系在仓） |
| Q3a | `chessboard\|棋盘` | grep **OK** | **24** | os3-reflection-positivity/REPORT.md:468–527（§6 专节 + O4/O5）、glue-b-blocking:485（Seiler LNP 159） |
| Q3b | `slab`（文件级） | grep **OK** | 12 文件 | os3、BREAK-IN-4/5、OS-P-INTERFACE、P3-CYLINDER-RESIDUE 等 |
| Q4 | `thermodynamic integration\|热力学积分` | grep **OK** | **7** | H04-THERMO.md（专文）：精确版=恒等式零内容、估计版= `(G)`、反向三角否决 |
| Q5 | `Fredholm` | grep **OK** | **187** | D08-LOGDET、CONTINUUM-FREDHOLM（`research/unconditional/CONTINUUM-FREDHOLM.md` 存在，`docs/` 下**不存在**——`ls` ENOENT 1 次，已定位正确路径）、SIXTH-SPECIES-SYNTHESIS:68（"Fredholm 零 = `m₁` 改名 = I2"） |
| Q6 | `donsker\|varadhan` | grep **OK** | **17** | INEQ-R5-donsker-varadhan.md（专文）、k0-LARGE-DEVIATION.md（LDP 速率曲率=H1） |
| Q7 | `balaban\|magnen\|rivasseau\|seneor` | grep **OK** | **846** | V04-RG.md（判决文）、C2-BALABAN-UNIFORM.md（411 行）、LIT-ERA-05、RG2/RG4/RG5/RG6 系列 |
| Q8 | `Hellmann` | grep **OK** | **7** | DF-CONTINUITY:193（`dΔ/dλ` HF 形）、D04:68、H06-VARM:35、M06-AUDITH06:43、WA-ALIGNMENT:161 |
| Q9 | `d l(n\|og) m\b\|dΔ/dβ\|dΔ` | grep **OK** | **19** | C3-UNIQUENESS:147–148、C1、MS-RG-ODE:28（`d ln m₁/d ln μ≡0`） |
| Q10 | `a\^k\|cutoff.suppressed\|标度破坏` | grep **OK** | **16** | `a^k`/`cutoff-suppressed` 形状 **0 命中**；`标度破坏` 命中 C2:68,:83–98、handoff-proof:147（两圈标度破坏 30–50%） |
| Q11 | `TWENTY-PATHS\|N07\|PROOF-COARSE\|FALSIFY-ENTROPY` 定位 | find **OK** | 4/4 | `research/new-paths/TWENTY-PATHS.md`、`research/v3/species/N07-STRONG.md`、`research/unconditional/{PROOF-COARSE-POINCARE,FALSIFY-ENTROPY}.md` |

**账本核对**（position-current.html §06，:335–364）：死路台账**没有**名为 "(F) route" 或 "slope route" 的条目；最近四条：`ODE / 里卡提`（G″ 不闭合）、`m(β) 单调性`（递减）、`(G) ∫A/Δ`（发散、`Δ_max≥0.276L²`）、`log det(1−K)`（det>0 恒真、留数没消掉）。§07（:397–402）20 路判决：**9 上界 / 9 循环 / 2 反例 / 2 实在**——(F) 属第 2 类（循环）的结构位置。窗内 `I1` 窗口口径：账本记 `β<3.777e−3`（KP `a*=0.618` 判据，:375）与 master 记 `β_A<0.88291×(1/48)`、`β_W<3.5316×(1/48)`（:130）并存，本轮不裁口径，仅注明两套数都 ≪ 靶 `β_W=2.5`。

---

## §2 · 路 (a) 文献：有限差分 + RP slab/chessboard

### §2.1 chessboard/slab 估计到底界什么（读到摘要的核验）
- **Biskup–Kotecký**，*Forbidden gap argument for phase transitions proved by means of chessboard estimates*，CMP **264** (2006) 631–656，arXiv:math-ph/0505011（**读到摘要**，arXiv abs 页直接抓取）。chessboard 估计的推广用途：在可调参数模型里**排除 shift-ergodic Gibbs 态取得"过渡隙"内的共轭热力学量值**（能量禁带）。⇒ 文献里 chessboard 的"gap"是**热力学共轭量的禁带**，不是转移矩阵谱隙。这与仓内 os3 REPORT 的判词一致。
- **Biskup**，*Reflection positivity and phase transitions in lattice spin models*，arXiv:math-ph/0610025（EPJB 综述）（**读到检索条目**，未读正文）。RP+chessboard 体系综述，主题=相变/长程序。
- **FILS**（Fröhlich–Israel–Lieb–Simon，CMP 1978）——chessboard 估计的经典源头 `[R]`（经仓内 os3 REPORT.md:477–480 转引；未逐字核）。
- **检索负面**：arXiv 站内检索 `all:"chessboard" AND all:"mass gap"` → **0 命中**（工具 OK，arXiv search UI 直抓）；arXiv API `all:"chessboard estimates" AND all:"mass gap"` → **0 命中**（工具 OK）。⇒ **没有任何 arXiv 论文标题/摘要同时含"chessboard"与"mass gap"**；无任何已发表结果给出 β-导数或 β-有限差分控制（它们给的是空间衰减/乘性/Peierls 型上界，os3 REPORT.md:493 逐字："chessboard 只给上界，不给聚类、不给 RG 步进定理。它是工具，不是证明 MG 的引擎"）。
- **RP slab 的隙方向**：RP 确实能把"衰减输入"换成隙输出（变分式 `m ≥ (1/n)log[⟨fθf⟩/⟨f²⟩]` 型；仓内 os3 REPORT 的 `Δ(n)=4logZ(n²)−logZ(2n²)` 检验即此），但该方向**消费**衰减输入（=需要先证一致衰减=隙级），**产出**上界或换算——与 (F) 需要的"局部乘性界 ⇒ β-差分一致"方向相反。
- ⚠ **数据质量警告**：本轮 WebSearch 首查返回 4 个 arXiv 号（2404.08975 / 2404.08664 / 2405.16998 / 2410.06514，均署名 Balaban–Imbrie–Jaffe 或 McKenna–Ranard 的"mass gap 估计"论文），逐一抓 abs 页核实**全部错配**（天体物理 / NLP / 光学 / 密码学）。`Imbrie–Jaffe–Brydges–O'Neill` 五作者 RP-隙论文经二次检索**不存在**。`McKenna–Ranard "Chessboard estimates for finite-range interactions"` 一文**未能在 arXiv 元数据中定位**（"chessboard estimates" 精确短语检索 11 条结果中无此文）。**上述四号一律不作引用**；涉及处标 `[M-未核]`。

### §2.2 Balaban 机输出什么（仓内已核 + 本轮检索）
- **仓内核验过的书目**（V04-RG.md:47–55，DOI 全列，`[M]` 级=书目核）：CMP **99** (1985) 75 / **109** (1987) 249 / **116** (1988) 1 / **119** (1988) 243 / **122** (1989) 175、355。终点 B8 = **Theorem 1：紫外稳定性**（有效密度 `ρ_k` 对步数一致、固定物理环面、绝对值型界）。
- **输出形状表**（V04-RG.md:66–77）：体积**不一致** `[O]`；力 Hessian 未证 `[O]`；IR 谱/团簇/质量隙**无** `[O]`；OS 重构**无** `[O]`。关键一句（:76–77）："`R` 运算的界是绝对值型，二阶场导数无符号控制；质量隙要的是 IR 尺度规范不变曲率的**下**界。`Hessian 控制 ⇒ 作用控制`（一次积分），反向不成立。"
- **`|∂_β(跑动耦合)| ≤ …` 型陈述**：本轮 WebSearch 两次 Balaban 定向检索 **TOOL-FAIL**（60s 超时 ×2）；未找到任何此类已发表陈述；仓内 `RG2/RG4/RG5/RG6/L4A` 系列的判词与 V04 一致（`β_RG` 不存在为一个数，定理形状 = `g∈(0,γ]` 常数未量化，V04-RG:30–34）。按 AGENTS §1：此为 **TOOL-FAIL + 仓内替代核验**，非 0 命中结论。
- **Gawędzki–Kupiainen**：4D 严格结果 = 无质量 `φ⁴₄` 的**红外渐近自由**（CMP 99 (1985) 197–252，DOI 10.1007/BF01212281 `[M]`）；"GK did not treat gauge theory"（LIT-ERA-05 §2.1 逐字，V04-RG:28）。流方向不匹配：`φ⁴₄` UV→IR 误差变小、4D YM 反向（V04-RG:107–112）。
- **Magnen–Rivasseau–Sénéor**：相空间展开构造 `φ⁴₄`（Lausanne 1984 会议文；读到检索条目，未读正文 `[M]`）。输出 = 构造/UV 稳定性；4D 质量隙非其产出。**Magnen–Rivasseau–Sénéor 的 2D/3D RG 输出形状 = 构造 + 稳定性 + 重整化常数递推**，与 DeepSeek "输出比值/流" 的说法相容，但仓内 V04 把要害钉在更硬的位置：缺的是**塔底一致单边方差界**，"RG 不是第二条武器，它是同一堵墙的另一种写法"（V04-RG:35–38）。

### §2.3 "机器不输出地板"声明的核实
- **反例 1（文献）**：**Göpfert–Mack**——CMP **82** (1981) 345–357，DOI 10.1007/BF01286362（*Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant*，**读到摘要**，经 WebSearch 返回的 Springer 条目）；姊妹文 = 质量隙文（仓内引作 CMP 82 (1982) 545，DF-CONTINUITY.md:131；本轮未取得其独立条目页，标 `[R]` 未逐字核）。聚合物/单极气体团簇展开**在 3D U(1) 全耦合区间输出隙**（=指数衰减=隙地板）——**"此类机器从不输出地板"按字面为假**。适用窗：3D、abelian、超可重整。
- **反例 2（仓内）**：MS-RG-ODE 的残数流闭式界 `f₁²(2.5)≥c·r₀`（比值/流型、(L,a)-一致、闭式）——**比值/流界已做出**，结论是它搬不动 `m₁`（RG 不变量）。
- **半真部分**：4D 非阿贝尔的 RG 机（Balaban 一支）确实只输出 UV 稳定性/绝对值界/跑动耦合窗口，无地板。⇒ **修正后的正确问句不是"问机器要比值/流界"（已有），而是"造一台输出塔底一致方差的机器"（= H1，C2-BALABAN-UNIFORM `(★)`、master §5.1 已钉）。**

---

## §3 · 路 (b) 文献：相对熵/熵率方法

### §3.1 熵方法给出隙的正面文献（读到的）
- **Shen–Zhu–Zhu**，*A stochastic analysis approach to lattice Yang–Mills at strong coupling*，arXiv:2204.12737，CMP **400** (2023) 805，DOI 10.1007/s00220-022-04609-1（**读到摘要**，arXiv abs 页直接抓取）。窗：SU(N) `|β| < 1/[16(d−1)]`（d=4 ⇒ `1/48`）/'t Hooft 标度 `βN`；方法 = Langevin 动力学 + **Log-Sobolev 与 Poincaré 不等式**；结论 = 唯一不变测度、指数遍历、**一大类可观测的指数衰减 = "strictly positive mass gap"（关联意义）**。⇒ **"LSI ⇒ 隙"在文献里成立且已入仓**（master §2.5 Lemma SC 即其推论：`β<1/48` 均匀 LSI ⇒ 粗模 Poincaré `(L,a)`-一致）。**窗**：仅强耦合端；输出为关联衰减，不是 `d log m/dβ` 的控制，更不是跨 β 的地板搬运。
- 经典链条（`[R]` 未逐字核）：Gross 1975（LSI ⇒ hypercontractivity；LSI ⟺ Dirichlet 形生成元的谱隙为教科书事实）；Holley–Stroock（J. Stat. Phys. 46 (1987) 1159 附近，随机 Ising 的 LSI/比热条件）——本轮 WebSearch 定向检索**0 可用链接**（工具 OK、命中 0），改以仓内 SZZ 体系与上述教科书事实承担，**不给未经核的 DOI**。
- **热力学积分** `d log Z/dβ = −⟨S_W⟩`：经典恒等式；**隙版**已由 H04-THERMO 判死（§0 判 5）：精确版零内容、估计版= `(G)`、反向三角否决。`log Z` 版有效的原因被 H04:147 逐字钉住：被积函数**定号**（Griffiths `⟨S_W⟩ ≥ 0`），`dF/dβ` 无定号。

### §3.2 "熵率控制本征值导数"的精确断点
- **断点在仓内有三层记录**（非本轮首见）：
  1. **条件式类不相交**（INEQ-R5-donsker-varadhan.md:3–9，`[NG]`）：DV 本身是连续泛函的 sup ⇒ lsc，可下界 lsc 量；但 `∂_β DV` 不自动保持 lsc。对**已另证为 lsc** 且逐点 `Ψ≤f₁²` 的泛函，强制引理用 `1/n` 宽 bump 给 `Ψ≤0`。残数 `f₁²=ν({m₁²})` 是 usc 原子，**不是任何 DV 率**。真正无前提的杀手是具体 `Var(S_W)/V` 的 2D 反例（`GW5-ACID.md §1.3`），不能把条件引理扩成所有熵率/LSI 的类杀。DV 能看到的只有**聚拢量**：压强/相对熵/转移算子主本征值——而"连 `Δ` 都需要 Courant–Fischer（上界机）才能给隙下界"（:13）。
  2. **三死类**（FALSIFY-ENTROPY.md:5–11，`[NG]`）：相对熵/Rényi(α>1) 凸 ⇒ 盲（LP 对偶缝隙恰 0）；微分熵/Fisher 需密度 ⇒ 原子处发散 ⇒ presence-only = `J→∞=H1`；log-Sobolev 的熵产生泛函是 Markov 半群泛函 ⇒ Theorem BE（齐次/零点/变号）⇒ 盲。见证：有限矩 `M0..M3` 在 `ε:0.2→0.02` 仅变 `5.3e−3` 而熵→−∞、Fisher~1/ε²（:46–60）。
  3. **LDP 速率曲率=恒等式**（k0-LARGE-DEVIATION.md:3）：`I''(0)=1/Var_μ(a₀)=G″_∞=m_eff²`（Gärtner–Ellis）⇒ "曲率≥λ" 逐字 = `Var≤1/λ` = H1；速率函数**命名**曲率而不**界**曲率。
- **本轮的独立表述**（[J]，供 DeepSeek 方对照）：熵方法控制的是 (i) 局部可观测的期望（Donsker–Varadhan/Legendre 类，lsc）与 (ii) Dirichlet 形生成元的谱隙（=测度的关联衰减率）。对**质量载体** `m=E₁−E₀`：(ii) 在 RP 体系里经 OS 重构确实连接转移矩阵谱隙——但那正是 SZZ 已做的强耦合窗结果，其常数 `λ` 是**窗内**的；把它沿 β 搬运需要的恰是 (F)，即目标本身。对**残数载体**：断点=usc/lsc 类不相交（INEQ-R5），**该断点已在仓内成文，不是本轮的新观察**。
- **文献中是否有"熵界控制谱量的 β-导数"**：本轮检索 `au:Mastropietro AND all:"mass gap"`（arXiv，0 命中，工具 OK）、`au:Bearer AND all:"mass gap"`（0 命中，工具 OK）、arXiv API `ti:"Reflection Positivity and the Mass Gap"`（0 命中，工具 OK）。**未找到任何已发表结果用熵/变分率方法界一个转移矩阵本征值差的 β-导数**（检索局限见 §5 表）。

---

## §4 · 路 (c) 文献：窗内标定 / 对团簇展开按参数求导 / 交叉点

### §4.1 对聚合物展开按参数求导（存在文献，读到摘要）
1. **Kotecký–Preiss**，*Cluster expansions of abstract polymer models*，CMP **103** (1986) 491–498，DOI 10.1007/BF01240679（**读到摘要**，经 WebSearch 返回的 Springer 摘要文本；正文付费墙，未逐字核）。"基于对 `log Z` 系数的直接估计的抽象聚合物模型团簇展开收敛新途径；一般收敛判据；Dobrushin 判据的改进。"⇒ KP 判据给出的是**对参数一致**的绝对收敛 ⇒ 在 KP 区内对 fugacity/耦合逐项求导合法（控制收敛）。这正是 DeepSeek 路 (c) 的合法性依据——**仓内 Z01/P3 体系用的正是 KP**（master §2.1：`u·K<1`、`K≤CN`、`Ĉ·C̄=8>7.0633`、窗 `β<3.53×(1/48)`），所以**逐项 `d/dβ` 在窗内是一个合法且便宜的 `[E]` 练习**，此前无人做过（Q10/Q4 检索 0 命中"对 Z01 展开求导"类前科）。
2. **Campanino–Ioffe–Velenik**，*Ornstein–Zernike theory for finite range random walks on Z^d*，PTRF **123** (2002) 384–406，DOI 10.1007/s004400100194（**读到摘要**，经 WebSearch 摘要文本；arXiv 号未验证——WebSearch 给的 `math/0207134` 抓取实为 twist-map 论文，**不引用该号**）。对 `p>p_c` 的有限程随机游走，质量（逆关联长度）`ν(p,β)` **是 β 与 p 的解析函数**——文献里最接近 (F) 的严格实例：**质量对参数的导数存在且解析**。但方法 = 再生结构/OZ，不是逐项 fugacity 求导；且该类模型无 AF，**没有 `a^k` 压低结构**。
3. **Dobrushin–Sokal 型解析性**（逆关联长度在团簇区的解析性，经典 `[R]`，未逐字核）：CIV 摘要自述"refines and extends results of Dobrushin … concerning analyticity of inverse correlation length"。

### §4.2 窗内标定能标出什么（仓内预算预填）
**第 222 轮更正**：上述 `√V`/`1/Δ` 预填只适用于 `|F'|` 的直接 Cauchy–Schwarz 路，不适用于 `log(a₀D_∞)` 的簇级数。GW4 实际给 `C_A=O(1)`、窗内余 `367×`；同坐标反事实全跨度乘积 `4.4470<9.6566`，不是超预算。真正的未闭合处仍是**窗外无一致界**。KP 窗到 `β_W=2.5` 的 Wilson 同坐标距离约 `165.475×`；若采用另一窗定义 `β_W<3.5316/48≈0.07358`，则约 `33.98×`。二者都不是坐标混用的 `662×`。

### §4.3 跨交叉点的隙：文献与口径
- **Savvidy 交叉点是曲率符号交叉，不是相变**：`β_W=2.5 (SU(2), N_t 无关的 T=0 设定)` 高于解禁闭交叉 `β_c≈2.30 (N_t=4)`/`≈2.43 (N_t=6)`（C2-MEASURED-SLOPE.md:22 引 LATTICE-TM-GAP:22–24）——但物理胶球区在 `T=0` 线上，`β_c` 只是有限温读数的伪约束（C2:60 已记"那是有限温屏蔽质量"）。仓内 `PROOF-SAVVIDY-EXACT-SIGN`（`λ_min=−g|B|<0`）与 MS-RG-ODE §3（交叉点 = 正则有限耦合点，`g²≈1.6–2.4`，残数流在该点正则）已定口径。
- **"两侧有隙 + 粘合"的严格文献**：3D `φ⁴` 临界点无一致隙（平凡事实，临界点关联长度发散）——对 YM 的 β 线不构成障碍，因为目标线不过相变点；未找到任何"两相均有严格隙并粘合"的已发表定理适用于本靶。检索：`au:Bearer`（arXiv，0 命中，工具 OK）、"chessboard+mass gap"（0 命中，工具 OK）。[NG-无定论]：不作"不存在"断言。

---

## §5 · 标度陈述的新颖性核查（关键项）：`|d ln(m/Λ)/dβ| ≤ c·a^k` 有没有人严格证过？

**答案（本轮检索面）：未找到。最接近的已发表对象是聚合物/OZ 类的"质量对参数解析"，没有 cutoff-压低的标度破坏形式。**

| 检索 | 工具 | 检索式（逐字） | 命中 | 判读 |
|---|---|---|---|---|
| arXiv UI | OK | `all:"chessboard" all:"mass gap"` | **0** | 无标题/摘要层撞车 |
| arXiv UI | OK | `au:Mastropietro all:"mass gap"` | **0** | Mastropietro 的 cutoff 无关质量界不用 "mass gap" 词面（其 1D/2D 费米体系结果形态 = 质量沿流不变/有限 renormalized expansion，`[M-未核]`，不给号） |
| arXiv UI | OK | `au:Bearer` / `au:Bearer all:"mass gap"` | **0 / 0** | Bearer–Lieb 3D 紧 QED 文不在 arXiv 元数据（论文号未验证，不引用） |
| arXiv API | OK | `ti:"Reflection Positivity and the Mass Gap"` | **0** | WebSearch 首查返回的该题名论文不存在 |
| arXiv UI | OK | `"Ornstein-Zernike theory" "random walks"` | **2** | Ioffe 综述 / Ioffe–Velenik 弹道相；CIV 2002 主文未在元数据（引 DOI） |
| WebSearch | OK | Göpfert–Mack CMP 82 | 1（confine 文，DOI 10.1007/BF01286362） | 隙文（545）为仓内转引 `[R]` |
| WebSearch | OK | Kotecký–Preiss | 1（DOI 10.1007/BF01240679） | 读到摘要 |
| WebSearch | OK | Campanino–Ioffe–Velenik | 1（DOI 10.1007/s004400100194） | 读到摘要；arXiv 号幻觉已剔除 |
| WebSearch | OK | Holley–Stroock LSI | **0 可用链接** | 改用仓内 SZZ 体系承担 |
| WebFetch | TOOL-FAIL | arXiv API `all:"reflection positivity" AND all:"mass gap"` | — | 60s 超时 1 次；同义 UI 检索补上（见上） |
| WebFetch | TOOL-FAIL | Springer `10.1007/BF01208454` | — | HTTP 404（该号不对应目标文） |
| 仓内 | OK | Q10 `a\^k\|cutoff.suppressed` | **0** | 仓内亦无该字面形式 |

**诚实的结构性判断** `[J]`：`|d ln(m/Λ)/dβ| ≤ c·a^k` 的**微扰版本是标准口语**（标度破坏 `O(g²)∝1/β`，C2:83 的外推即用之），但作为**关于有隙格点理论的严格定理**，文献检索面（上述 10 条）与仓内均无。**同时必须说清它的逻辑地位**：在 `m/Λ = R(Λa)` 记号下（C3-UNIQUENESS §4），该陈述 = "`d ln R/d ln Λa = O((Λa)^k)`"，而 `R` 的 `Λa`-无关性**就是** C1 的斜率等式内容（C3:152）——即该陈述是"标度成立"的定量版；**若把它当作 (F) 的可证前件，它不弱于靶**（证明它需要先有隙且需控制 `R` 的有限-β 行为，而实测窗内 `R` 明显依赖 `Λa`，C3:126、C2:8）。⇒ **新颖性成立（没人证过），但把新颖陈述当前件不降低难度**——这与 D04 命题 2.1 的"假设⟺结论"是同一诊断在标度层的重演。

---

## §6 · 汇总先验 `[J]` + 最快杀死/确认实验

先验定义沿用 TWENTY-PATHS §0 口径（专职十年计划给出完整搬运/证明的概率）；与本轮三条路对照 DeepSeek 的 (a) 15–20% / (b) 10%。

| 路 | DeepSeek 先验 | 本轮 `[J]` | 一行理由 | 最快杀死/确认 |
|---|---|---|---|---|
| **(a)** 有限差分 + RP slab | 15–20% | **3–5%** | 有限差分版=H04 已判（反向三角否决 + `(G)` 货币）；RP/chessboard 文献 0 次给出 β-差分控制（arXiv 0 命中）；"局部-in-β ⇒ 全局"正是 DF C1 闭性缺口（便宜的一半不可达，够用的一半=目标） | **杀**：对 2D YM SU(2) 跑一遍 (a) 的机器（`m=3g²/8` 精确、`f₁²≡0`）：若机器对残数载体给出"地板"，它在该精确模型上必须输出 0 而不得 0 ⇒ 证伪其酸测合法性。**确认**：用仓内 `tm_gap_results.csv`（SU(2) `β=2.2/2.5/2.8`，L=4/6）算 `|Δ ln m|/Δβ` 有限差分（一条链，无新物理）——若 `O(1)`，则 (F)-质量的**物理内容**在窗内成立 `[E]`，可作记账 |
| **(b)** 相对熵/熵率 | 10% | **2–3%** | 具体 `Var(S_W)` 率在 2D 上正而残数恒零（GW5 的直接反例）；lsc 强制引理只封“已另证 lsc 且 Ψ≤f₁²”的 minorant，不封所有熵/LSI 泛函；SZZ 的 LSI⇒隙仍只在窗内，跨 β 搬运=(F) | **杀**：执行 GW5 的具体 2D 反例；**条件封格**：只有先证 Ψ 为 lsc minorant 时才能用 `1/n`-bump。**确认（负向）**：核 `Lemma SC` 的 Poincaré 常数是否在窗内给显式 `m(β)≥f(λ)` |
| **(c)** 窗内标定 | （未给） | **作为标定练习 ≈85%；作为搬运路 ≈2%** | **第 222 轮已执行并更正**：簇级数给 `C_A=O(1)`，窗内不超预算；`D07` 的 `√V/Δ` 不是 `log(a₀D_∞)` 的形状。窗外无一致界仍是缺口 | **已完成**：GW4 给显式 `C(β)`；结论是“窗内便宜、窗外未知”，不是“窗内超预算” |
| 桥引理 ((F)+种子 ⇒ 地板) | "≈确定" | **条件句 [T] 已证（D04 推论 2.2）；作为"免费"= 0%** | 假设逐字等价于结论；积分版另有 H04 反向三角独立否决 | 已完成（D04/H04），无需新实验 |
| 支持性声明（RG 机=比值/流非地板） | — | **字面假 / 精神半真** | Göpfert–Mack（3D abelian 地板）+ MS-RG-ODE（比值/流界已有且搬不动 `m₁`）；Balaban 支输出=UV 稳定性 | 已核（§2.3），无需新实验 |

**对 DeepSeek 提案的一句话总判**：提案的**原理层**（(F)、桥、log-det、(G)）在 D04/D07/D08/H04 四个文件里已经**逐字逐判**，无需重跑；提案的**新信息**只有两处字面空白——①"(F) 用在质量载体上"（组件齐全、路线同死，价值=一个酸测更干净的坐标系：自由场 `m=0` 使 (F)-质量前件显假，比残数载体的 `0≤0` 酸测更锋利 `[J]`）；②"窗内逐项标定 `C(β)`"（合法、便宜、预期产出为窗内超预算的定量反证）。提案的**文献层**（RG 机输出形状）需要按 §2.3 修正："问机器要比值/流界"已经做过（MS-RG-ODE），真正的缺件是塔底一致方差界（H1），这一点 DeepSeek 的叙述没有触及。

---

## §7 · 工具状态汇总（AGENTS §1，本轮全部）

| 工具 | 状态 | 说明 |
|---|---|---|
| Bash `grep/find/ls/sed` | **OK** | §1.6 表 Q1–Q11 全部有效返回；1 次 `ls` ENOENT（`docs/CONTINUUM-FREDHOLM.md` 不存在）⇒ 定位到 `research/unconditional/CONTINUUM-FREDHOLM.md`（非 TOOL-FAIL，路径纠错） |
| Read | **OK** | 14 个仓内文件全文/定段读取（D04/D07/D08/DF/C1/C2/C3-COEFF/C3-UNIQ/V02/B04/B01/G05A/N07/MS-RG-ODE/INEQ-R5/FALSIFY-ENTROPY/TWENTY-PATHS/P10/V04/os3§6/position/master） |
| WebSearch | 7 OK / 2 TOOL-FAIL | OK：Göpfert–Mack、MRS、McKenna–Ranard（ID 幻觉）、Imbrie-Jaffe（0 命中）、Kotecký–Preiss、CIV、Göpfert-Mack DOI。TOOL-FAIL：Balaban 定向 ×2（60s 超时）⇒ 以仓内 V04 的 arXiv API 核验 + 二次换式补 |
| WebFetch | 12 OK / 3 TOOL-FAIL | OK：arXiv abs ×6（2204.12737、math-ph/0505011 核实；2410.06514/2212.13366/1905.03275/math/0207134 **证实为幻觉号**）、arXiv search UI ×7（含 3 个 0 命中）、arXiv API ×2（1 个 0 命中）。TOOL-FAIL：arXiv API 超时 ×1、Springer 404 ×1、abs 2410.06514 首次超时 ×1（重试成功） |
| pwsh / py | 未用 | 本轮无 MC、无数值轮（纯解析+检索） |

**引用等级一览**：读到摘要 = Shen–Zhu–Zhu (arXiv:2204.12737)、Biskup–Kotecký (arXiv:math-ph/0505011)、Kotecký–Preiss (DOI 10.1007/BF01240679，经检索摘要)、Campanino–Ioffe–Velenik (DOI 10.1007/s004400100194，经检索摘要)、Göpfert–Mack confinement 文 (DOI 10.1007/BF01286362，经检索摘要)。读到正文 = 无（付费墙/未取全文，均逐条标注）。`[R]` 未逐字核 = FILS、Holley–Stroock、Göpfert–Mack 隙文 (CMP 82, 545)、Dobrushin–Sokal、GK CMP 99、Balaban 系列（仓内书目级核，V04-RG:47–55）。**已剔除的幻觉引用**：arXiv 2404.08975 / 2404.08664 / 2405.16998 / 2410.06514 / 2212.13366 / 1905.03275 / math/0207134 七个号（逐一抓页证伪）；任何引用均以本表核实条目为准。

*本文件：无 QED 于未证之物；`[T]/[C]/[NG]/[M]/[E]/[J]/[R]/[O]` 逐条标注；仓内数字全部带 file:line；外部引用全部带 ID/DOI + 阅读等级；Clay 0%。*

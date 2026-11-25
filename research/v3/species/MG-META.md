# MG-META — 元层生成器：把"**死因**"当轴，而不是把"主题"当轴

> 2026-11-25 · 元层轮（`meta-generator`）。**No QED. Clay 0%.**
> 只新建本文件（`research/v3/species/MG-META.md`），未改任何其它文件，未 `git commit`。
> 标签：`[T]` 真定理/恒等式 · `[C]` 条件 · `[NG]` 具名类 no-go · `[M]` 精确缺失步 · `[E]` 数值证据 · `[J]` 判断。
>
> **记账规则（全篇遵守，可被重算）**
> 1. 一个族**只计一次主类**。判词句里被点名为"决定性"的那条机制当主类，其余列进**副类**并单独计数（表 1 两列）。
> 2. 族清单 = `ch08-immune.tex` 的 **40** 行 + 精细记账缺口 **1**（`SIXTH-SPECIES-SYNTHESIS.md` §3.1：TCI/RB/DM
>    全计时是 41/42/43，本文件取 TCI 为第 41 格，RB/DM 已在表内另有格位）+ 本轮新路 **24**。
>    ⇒ **65 条判词，逐条归类，无抽样。**
> 3. **本文件能写的最强话是"在具名类内不存在"**；"我没找到"不写成"不存在"。
> 4. α/β/γ 按 `SIXTH-SPECIES-CHARTER.md` §5：**α** 第五数据模（可计算 ∧ 地板敏感）· **β** 新等价坐标
>    （`𝔉_object` 正性 ⟺ gap，且带**严格强于 H1** 的新机器）· **γ** 决定性否证（把在册逃生口判假 / 把候选压成 H1 改名）。
> 5. 五条件 = `SIXTH-SPECIES-CHARTER.md` §2（非微扰 / 控号下界 / 非循环 / 非 SUSY 且败自由场 / `(L,a)` 均匀）；
>    酸测五点 = 同文件 §3。

---

## §0 判词（12 行）

1. **死因类 17 个**（D1–D17，加 1 个准类 D18），覆盖 **65** 条判词；主因分布见 §1 表 1。
2. 主因**高度集中**：`D2 循环 21` + `D1 方向反 11` + `D6 空前提 8` = **40/65（62%）**。三者恰是
   `HISTORICAL-TOOL-SWEEP.md`「Verdict up front」的两族 + 真空性：*Family I 上界机 / Family II 输入循环 / 前提空*。
3. **三个类从未当过主因**（只作副因）：`D12 无尺度分离`（副 1）、`D15 定标不匹配`（副 6）、`D17 恒等式≠不等式`（副 2）。
   ⇒ "决定性死因"的种类比"死因表"窄得多——这是本文件第一个可核新事实。
4. **D7（常数/标度）在 41 族里 0 次主因**，直到 `NS-NOSUBTRACTION`（`Δ_dyn=βE−β²Var`）与 `KT-RHO2`（`C=2` 是类常数）
   才第一次当主因，且**一攻就死**。⇒ `β`-解耦量这一层是最后才被攻的层。
5. **T2 造出 17 行逃生表**（16 个死因类 + 准类 D18）；其中 **7 个候选机制"在具名类内未出现"**，
   5 个是弱候选（按构造立刻落回邻类），**真正值得测的只有 2 个 + 1 个名义未出现**。
6. **T3 只测 3 个**（纪律：先写判死线、先做方向账）。三个全死，且**死在三个不同的类**：C1→`D1`、C2→`D18+D6`、C3→`D5`。
7. ⇒ 那句"在'死因类'这个轴上，可逃的机制也已被穷尽" **本轮不成立**：`7 − 3 = 4` 个未测，
   且 `D7` 与 `D2` 两行是"满格 4/4 已出现"。**不许写那句话**，理由与补测清单见 §3.4。
8. T2 的真正强结论不是穷尽，而是：**同时逃过四道"机制形状"滤波（D1/D5/D6/D10）的对象，整个语料只有一个**——
   `A2` 的 `G″=S−D`（场空间、有向、自由场静默、无 `β²` 项）；**但它的一致化仍落回 `D2`（`reg` A2：残留 = `H1`），
   缺口是定量的 `D/S≈2%`，不是结构的**。⇒ 形状不是瓶颈，**定量才是**。
9. 新增**准类 D18（输出值域在阈值处稠密 ⇒ 无法排除）**，本轮由 C2 发现；它与 `D7` 是否独立，需下一轮裁决。
10. 顺手抓到**记账缺口**：`HORMANDER.md`（子椭圆/Fefferman–Phong）、`MICROLOCAL-BUILD.md`（解析微局部）、
    `SOS-CERTIFICATE.md`（moment-SOS）三条物种**已判死但未入免疫图**（同型问题见 `BC-BETHE` 被"并入族 3"）；
    疑似第四处：`O-MINIMAL-UNIFORMITY.md`（o-minimality/proof mining，自记 2/5）——**未核，不写成结论**。
11. **下一轴建议：等价性轴**（`H1` 的**充分性**），不是量轴、不是元数学轴。理由见 §4。
12. Clay 0%。无 QED。

---

## §1 T1 · 死因类清单（有限、带计数、带机理）

> 每一类的"机理"一行是**从判词句里抽象出来的**，不是我的猜测：引文在表 2 表 3 里逐条给出。
> **主/副** 两列：主 = 该族被它**决定性**杀死；副 = 同一判词句里并列点名、但不是决定性那一半。

| 类 | 名称 | 主机理（一句话） | 主 | 副 | 代表条目（族号见 §1′） |
|---|---|---|---|---|---|
| **D1** | 方向反 | 工具的输出住在目标量的**另一侧**：目标要 `inf` 的下界，而工具只能给"对 `min` 取上界"（变分/Ritz/宽度/交错皆然） | 11 | 7 | T6 变分上界、T14 下界机器（"下界机器、上界靶"=反隙）、T20 `κ_NG=−1/(6π)`、T40 矩、N19 `λ₁R_eff≤2` |
| **D2** | 循环 | 证书的**求值本身**需要目标（或它的等价物：`E₂`/阈值/`a→0` 连续统/`L→∞` 隔离/`1/χ`） | 21 | 9 | T5 改名=gap、T23 `F(0)=0`、T33 诊断、T39 `ρ²≤C² ⟺ e₁≥t/(1+C²)`、N2 X1 热模隙`=β×`转移隙 |
| **D3** | 尺度/regime 错 | 工具的有效域是 `β<1/48` 或 UV 半边或 `β→∞`；答案住在 `β≈2.5` 的 IR 边 | 3 | 2 | T4 错半边 UV、T7 强耦合专享、T19 长线构造 |
| **D4** | 道/旋错 | 读的不是 `0⁺⁺`／N-ality-0／`trF²`-**分辨**的那一道（跨道对称、N-ality 盲、operator-blind） | 1 | 6 | T2 跨道（`local trF²` 盲、BPS `Z=0`）；副：T9、T17、T18、T19、T24、T26 |
| **D5** | 弱-* 盲 | 输入是 `ν` 的弱-* 连续泛函（矩/有限数据/有限 `t`）⇒ 看不见**支撑边界的不连续点** | 3 | 1 | T29 复分析工具、T30 FALSIFY、N17 DM（`d_N~e^{−0.4595N}`，双向不可判） |
| **D6** | 空前提/空量化 | 前提类在目标对象上**为空**（`ℋ₋={0}`、FKG 缺、`(logρ)″<0`、无正测度），或自由场已满足全部前提 | 8 | 4 | T10 四把手全空、T26 一致性被自由场满足、T37 SW、T38 LX、T41 TCI |
| **D7** | 常数/标度障碍 | 展开里有 `β²` 负项，或承重常数只在**类内**成立 ⇒ 不可能 `β`-一致 | 2 | 1 | N21 `Δ_dyn=βE[H₁]−β²Var(G₁)`、N23 `C²(𝒦)`（DD=4 / ND=9 / 变密度≥5.71 / 双井=∞） |
| **D8** | 只 decoupling 非 floor | 证明"消去/解耦/中性"，而地板要的是"**耦合下**有下界" | 1 | 0 | N1 SAV-CANCEL（`φ₋·∇O=0`，耦合 `−0.32` 有界，但"negative mode was never the atom's obstruction"） |
| **D9** | 符号由 FKG 定而 `SU(N)` 无 FKG | 比较/符号的**前提**是正关联；`SU(N)` 上无全序、Weingarten 面权带号 ⇒ 前提空 | 1 | 1 | T3 非阿贝尔符号；副：N20（BC-BETHE 的"符号由 FKG 定、几何定不了"） |
| **D10** | 自由场不静默／同值 | 机制的输入量在 `g=0` 上**同值或已取极值**：要么给无隙场发正地板（误报），要么"两边都 0"（空真） | 3 | 7 | T1 自由场盲正性、N4 张量网络（`L=2` 有限体积隙过不了酸测）、N7 NCG（谱作用=经典 YM） |
| **D11** | 被反例**否证**（假，不是弱） | 不是"不够强"，而是一个数字把它打掉：假设在第一步就失败 404.8×，或大类上有反例 | 3 | 3 | T35 PS（`ε₁/dist=404.8`）、N12 FOSP-A（`SU(3)` ~26% 2×2 子式负）、N18 `(NA)` 双向为假 |
| **D12** | 无尺度分离 | 工具要的层级（`ℓ/L→0`、`log₂(ξ/a)≫1`）在答案所在的尺度上**不存在** | 0 | 1 | 副：N24 BN-N3（`log₂(ξ/a)∈[−0.9,0.4]`，0–1 步）；同型：`reg` §F1 相位单元 `(ℓ/L)²→1` |
| **D13** | 载体结构性恒零/退化 | 载体在对称切片上**恒等于 0** 或退化 ⇒ 无信号（不是"弱"，是"没有"） | 3 | 1 | T34 HC（Plancherel 绝对连续、`atom_at_floor=0`）、T36 TR（`Δ_trop≡0`）、N14 CCM（`H□` 斜自伴、`tr≡0`） |
| **D14** | 投射错位 | 读的是同一算子的**另一个泛函**：只交错特征值不交错重叠、"原子 ⟂ gap"、控核不控对齐 | 3 | 6 | N10 高阶 BE、N13 FOSP-B、N16 SL-TEST；副：T21 WI、T29、T40 |
| **D15** | 定标不匹配 | 量纲/体积强度/单位（`χ` vs `χ/V`、`(L,a)` 单位）与目标不一致 ⇒ **不可能**一致 | 0 | 6 | 副：T4（finite box=上界）、T12（单元不张量化）、T19、T37、N4、N24（N1 类型错：值 vs 变差模数） |
| **D16** | 类型错位（整数/拓扑 vs 实数下确界） | 输出是整数、局部常数或 UV 尺度，目标是 IR 的**实数下确界** ⇒ 刚性 | 2 | 3 | T8 指标=0、N6 TQFT/Langlands（metric-independent 类=自由值） |
| **D17** | 恒等式 ≠ 不等式（读 ≠ 力） | 输出是**精确恒等式**（或"读"），没有方向；要力就必须是不等式，而语料里每个不等式都单向且方向不对 | 0 | 2 | 副：T23（"every exact identity is a read"）、T33（诊断是读）；同型：`A3/A6/A8/A9` 全部"看得见原子、不给地板" |
| **D18**（准类） | 输出值域在阈值处**稠密** | 不是看不到，而是**看到了也不判**：不变量的取值范围在判定阈值处没有间隙（II₁ 因子的维数像 = `[0,1]`；`C²(𝒦)` 的像 = `[0,∞)`） | 1\* | 0 | C2（本轮 §3.2，**不计入 65**）；疑似同族：`KT-RHO2` 的 `C` 是类函数（那里记为 D7） |

**表 1 的两条读法.**

- **(i) 主因集中是结构性的，不是记账偶然.** `D2+D1+D6 = 40/65` 与 `HISTORICAL-TOOL-SWEEP.md` 的两族分类
  （Family I = `ν` 的线性/凸泛函 ⇒ 上界；Family II = 正确的逻辑形状但死在 regress/Haag/刚性/连续性）**独立吻合**，
  而那份 sweep 是按**历史年代**分的、本表是按**判词**分的 ⇒ 两条独立路径给出同一张图。`[J]`
- **(ii) `D7` 晚出现.** 41 族里 0 次主因，是因为承重常数 `C` 与 `β<1/48` 直到 `SIXTH-SPECIES-SYNTHESIS` §3.2/§3.3
  才被当成"定理级对象"审；一旦审，`C=2` 立刻降级为**类常数**（`KT-RHO2`）。⇒ **"常数是不是结构常数"是一条独立的杀死线**，
  此前没有被系统化。`[J]`

---

## §1′ 逐族归类（**41 + 24 = 65 行，一条不漏**）

> 出处缩写：`ch08` = `docs/paper-v2/ch08-immune.tex` 表 `tab:immune`；`zoo` = `research/v3/PARADIGM-ZOO-SUMMARY.md` §0；
> `R1` = `species/ROUND1-VERDICTS.md`；`reg` = `docs/new-tools-registry.md`；`SYN` = `species/SIXTH-SPECIES-SYNTHESIS.md`。
> 判词只引**决定性**那半句（原文英文不译，免得改味）。

### 表 2 · 免疫图 40 族 + 第 41 格（TCI）

| # | 族（出处） | 决定性判词（引文） | 主 | 副 |
|---|---|---|---|---|
| T1 | 自由场盲的正性/解析性（`ch08` #1） | "`g=0` free field has `f₁²=0`; Schur-convex master obstruction" | D10 | D1 |
| T2 | 跨道对称（BPS/SUSY、1-形式反常、非可逆对称、对偶 Meissner，`ch08` #2） | "wrong channel: line/flux/charged sector, or spin-1; local `tr F²` is blind" | D4 | D16 |
| T3 | 非阿贝尔符号（random flow、FK、FKG、spectral independence、Picard–Lefschetz，`ch08` #3） | "no lattice order on the quotient; `SU(N)` Weingarten face weights signed; FKG wrong-signed" | D9 | D6 |
| T4 | 错半边 UV（L3、`H2`-global、梯度流、Dirichlet 比较、集中、热核、验证数值、step scaling、有效弦、`S³`，`ch08` #4） | "orthogonal to the IR spectral edge; finite box `=` upper bound" | D3 | D15, D14 |
| T5 | 改名 = gap（H1、NEO/ESI、I2、束缚能、单圈 L3，`ch08` #5） | "the object is the gap restated, in another vocabulary" | D2 | — |
| T6 | 上界方向（Bogoliubov、Hylleraas、FSS、Temple/Weinstein、GEVP，`ch08` #6） | "variational methods give upper bounds; threshold oracle circular" | D1 | D2 |
| T7 | 强耦合专享（大 `N`/Chatterjee、SZZ、随机面/CPS `\lvert β\rvert<1/48`、中心涡旋，`ch08` #7） | "cannot reach `β≈2.5`; method radius is the obstruction" | D3 | — |
| T8 | 指标 `=0`（NCG、localization、BPS、Witten 旋量，`ch08` #8） | "rigidity no-go: discrete locally-constant index `=` free value zero" | D16 | — |
| T9 | 精确归约/对偶（ABC A-leg 八路，`ch08` #9） | "all eight `→` H1 or misdirection"（含 `S_N≥S_∞` **反号**） | D2 | D1, D3, D4 |
| T10 | 构造性强制（ABC B-leg，`ch08` #10） | "region-complete but sign-null: open cell has all four handles null" | D6 | D1 |
| T11 | 新不等式/容量（ABC C-leg、`C-BBS-*`，`ch08` #11） | "capacity floor is the Poincaré constant of the 1D marginal `=1/χ=G″=` H1, necessary not sufficient" | D2 | D6, D9 |
| T12 | 单元局部化（`LOC-*`，`ch08` #12） | "free massless satisfies every premise; single-cell control does not tensorize to `L→∞`" | D6 | D2, D15 |
| T13 | 搬运/单调（`MS-*`，`ch08` #13） | "every transport is one-sided (upper bounds) or local (a `β`-disk)" | D1 | D2 |
| T14 | 下界机器（`LOWER-AMPLITUDE-MECHANISMS` 八件，`ch08` #14） | "lower-bound machine, upper-bound target … a lower bound on correlation `=` lighter mode `=` **anti-gap**" | D1 | — |
| T15 | Feshbach/Krein/squeeze/Dirac-η/谱流（`FORGE-*`，`ch08` #15） | "two-block lemma reduces to TCI `=` H1; Krein count `≠` IR min; spectral-flow/eta index rigid" | D2 | D16 |
| T16 | `k=0` 有效曲率（`K0-*`，`ch08` #16） | "finite-`V` strict positivity is trivial … uniform strictness `=` gap" | D6 | D2 |
| T17 | Bethe–Salpeter 双胶子（`BS-*`，`ch08` #17） | "Temple needs the threshold oracle (circular); RP sign-blind; area law is the wrong sector" | D2 | D1, D4 |
| T18 | 弦张力 `→` 胶球（`TRANSMUTE-*`，`ch08` #18） | "the available inequality points the wrong way"（加扇区不符） | D1 | D4 |
| T19 | 长线构造（`CONSTRUCT-*` 八路，`ch08` #19） | "every construction is either strong-coupling, dual-transformed into the wrong channel, or one-datum phenomenological" | D3 | D4, D15 |
| T20 | 通量管 `κ`（`KAPPA-*`，`ch08` #20） | "NG determinant has the wrong sign (`−1/(6π)`); `κ=` H1 rescaled" | D1 | D11, D2 |
| T21 | WI/IR 原子签名（`WI-*`，`ch08` #21） | "detects single-exponential dominance, not atom `+` continuum; uniformity `=` H1" | D2 | D14 |
| T22 | Savvidy 元胞（`SAVVIDY-*`，`ch08` #22） | "the surviving margin is near-critical and IR-marginal `=` H1" | D2 | — |
| T23 | 不动点方程（`SLOT-*` 八路，`ch08` #23） | "`F(0)=0` in every coordinate; every exact identity is a read; the missing object is a `Q`-independent source at the zero `=` H1" | D2 | D17 |
| T24 | Roof（`ROOF-*` 八路，`ch08` #24） | "kinematic orthogonality; 0-form 't Hooft impossible; the strut is H1" | D2 | D4 |
| T25 | Base（`BASE-*` 八路，`ch08` #25） | "every base is a known type … and each type is eliminated by a theorem" | D2 | — |
| T26 | 不可能基（`IMPOSSIBLE-*` 八路，`ch08` #26） | "consistency satisfied by the free field; indices rigid; entanglement wrong channel" | D6 | D16, D4 |
| T27 | Star/RG 吸引子（`STAR-*`，`ch08` #27） | "it reduces the problem to `N_f=0` outside the conformal window, a dynamical `β`-function fact" | D2 | D3 |
| T28 | 共形窗 / `(G1)`（`G1-*`，`ch08` #28） | "the nonperturbative statement is `(G1)` itself" | D2 | — |
| T29 | 复函数论工具（`TOOL-*` 八路，`ch08` #29） | "values/moments are integrals (OPE-RHO blind) … only Hardy detects presence, not weight" | D5 | D14 |
| T30 | 否证 ONLY-EXIT（`FALSIFY-*`，`ch08` #30） | "finite data is weak-`*` continuous; atomless measures are dense" | D5 | — |
| T31 | 直接轰 H1（`BLAST-*`，`ch08` #31） | "every route collapses H1 to one marginal inequality … the atom half is untouched" | D2 | — |
| T32 | 证书组装（`CERT-*`、`ROUTE-NEW-*`，`ch08` #32） | "the certificate is H1 plus A3, not a new object" | D2 | — |
| T33 | 诊断（`WI-P3-*`、`WI-RESLOPE-*`，`ch08` #33） | "each is a read gated by `L→∞` `=` H1" | D2 | D17 |
| T34 | HC 球函数（`ch08` #34 = `R1` §1） | "the split-type Plancherel measure is absolutely continuous with no discrete series, so there is no atom" | D13 | D10 |
| T35 | PS 伪谱/非正规（`ch08` #35 = `R1` §2） | "Rouché needs `V<dist/ε₁=2.47×10⁻³<1` … a single site already fails by `404.8×`" | D11 | D10, D14 |
| T36 | TR 热带/min-plus（`ch08` #36 = `R1` §3） | "`S_W≡0` on the constant-modulus slice … `Δ_trop≡0` against the exact `Δ=0.6444`" | D13 | D10 |
| T37 | SW 集合型/小球（`ch08` #37 = `R1` §4） | "`(logρ)″<0` on all of `(0,2)`, so the hypothesis is empty for every `β` and every Gaussian" | D6 | D15, D1 |
| T38 | LX Lorentzian/M-凸（`ch08` #38 = `R1` §5） | "there is no positive measure and log-concavity has no referent" | D6 | D13 |
| T39 | RB 比值专用（`ch08` #39 = `reg` §B19） | "`ρ²≤C² ⟺ e₁≥t/(1+C²)`, so the window inherits H1 and is not an independent predicate" | D2 | D1, D7 |
| T40 | WA 源对齐四机制（`ch08` #40 = `reg` §B20） | "direction (the double commutator gives moments, hence upper bounds)" | D1 | D14, D6, D5 |
| T41 | **TCI**（第 41 格；`reg` §C、`SYN` §3.1 的"全计"读法） | "free massless field `ℋ₋={0}` ⇒ `∀φ∈ℋ₋` 空量化 ⇒ (TCI) 空洞通过"；条件版 "退化成 H1 本身" | D6 | D2 |

### 表 3 · 本轮新路 24 条（`zoo` + `reg` §B14–B22 + Round-2/3 派路）

| # | 路（出处） | 决定性判词（引文） | 主 | 副 |
|---|---|---|---|---|
| N1 | SAV-CANCEL 负模消去（`zoo`） | "exact decoupling … but it is a **decoupling, not a floor**; the negative mode was never the atom's obstruction" | D8 | — |
| N2 | X1 modular/代数 QFT（`zoo`） | "local algebras type III₁, Connes spectrum `ℝ`; global thermomodular gap `=β×`(transfer gap) `=` H1" | D2 | — |
| N3 | X2 factorization algebra（`zoo`） | "gap `⟺` locally constant (**opposite sign**)" | D1 | — |
| N4 | X3 张量网络（`zoo`） | "4D `L=2` finite-volume gap fails free-massless acid test" | D10 | D15 |
| N5 | Z1 holography（`zoo`） | "the IR cap's scale **is** the confinement scale; gap `=` cap/thermal scale `→0`" | D2 | — |
| N6 | Z2 TQFT/Langlands/cohomology（`zoo`） | "metric-independent class `=` free value (Witten index constant while gap spans ×16)" | D16 | — |
| N7 | Z3 NCG（`zoo`） | "spectral action `=` classical YM (acid-vacuous)" | D10 | — |
| N8 | Z4 MERA（`zoo`） | "scale-invariant MERA **is** the critical fixed point" | D1 | — |
| N9 | Y1 下半连续（`zoo`） | "Portmanteau gives **upper** semicontinuity; LSC false in weak-`*`" | D1 | — |
| N10 | Y2 高阶 Bakry–Émery（`zoo`） | "fixed-jet double well: jets `w`-independent while `λ₁~e^{−w²/2}→0`" | D14 | D2 |
| N11 | Y4 规范不变 Feshbach/RG（`zoo`） | "valid, exactly gauge-independent, non-vacuous, **but Gaussian thermodynamics**" | D2 | — |
| N12 | FOSP-A 融合序 TP2（`zoo`） | "`SU(3)` fails TP2 (~26% of 2×2 minors negative)" | D11 | — |
| N13 | FOSP-B 外平方（`zoo`） | "gap `=μ₂−μ₁` (insufficient); double well `Δ` rises while gap dies" | D14 | — |
| N14 | CCM 约束环记忆（`zoo`） | "all `E_i` orthogonal projections ⇒ `H□` skew-adjoint, `tr≡0`, kinematic only" | D13 | D10 |
| N15 | ICA ratchet（`zoo`） | "sign wrong (mass shrinks `D`)" | D1 | — |
| N16 | SL-TEST 谱身份簇（`zoo`） | "exact 2×2 battery shows atom `⟂` gap" | D14 | — |
| N17 | DM 第五数据模（`PROFILE-FACET-DATAMODE`） | "`(spectral datum D) ⇒` [floor-blind] XOR [circular]; the only exit leaves the spectral category" | D5 | D2 |
| N18 | `(NA)` 命名不等式（`R1` §7 / `NA2` / `CT-CHAIN`） | "`(NA)` 为**假（双向）**，不是'仅一处反号'"；唯一活口 `rank G=1 ∧ 单 link` 上恒零 | D11 | D6 |
| N19 | QG 图约化/有效电阻（`reg` §Q4） | "`λ₁·R_eff(a,b) ≤ 2`，常数 2 精确 ⇒ **电阻上界 `=` 隙上界**" | D1 | D2, D10 |
| N20 | BC 树比较（`reg` §B22） | "比较的符号由 **FKG** 决定，几何决定不了；`SU(N)` 无 FKG ⇒ 前提空"；`β=2.5` 处界**恒为 0** | D6 | D9, D10 |
| N21 | NS 无减法靶（`NS-NOSUBTRACTION`） | "`Δ_dyn=βE[H₁]−β²Var(G₁)` ⇒ 正项 `β⁰`、负项 `β²` ⇒ 本靶**不可能 `β`-一致**" | D7 | D11 |
| N22 | DT Dirichlet 塔（`DT-DIRICHLET`） | "`L_s=π/√e₁` ⇒ 退化；一侧比较**不给**比值界"（floor 半 `[C]` 活） | D2 | D10 |
| N23 | KT `C=2` 火烤（`KT-RHO2`） | "`C=2` 作为**结构常数**是死的，作为**边界条件类常数**是活的"（`C²`：DD=4、ND=9、变密度≥5.71、双井 `∞`） | D7 | D11 |
| N24 | BN beyond-ν 残余（`BN-RESIDUAL`） | "'`beyond-ν ≡ H1`' 必须撤销；总数是 `非∅ ⟺ Clay`"（N1 类型错 + N2 方向反 + N3 无尺度分离 + N4 纤维内外错位） | D2 | D15, D1, D12, D14 |

**表 2 / 表 3 的合计（可重算）**

```
主类:  D1 6+5=11   D2 16+5=21   D3 3       D4 1       D5 2+1=3    D6 7+1=8
       D7 0+2=2    D8 0+1=1    D9 1      D10 1+2=3  D11 1+2=3   D12 0
       D13 2+1=3   D14 0+3=3   D15 0     D16 1+1=2  D17 0       D18 1*（C2，§3.2，不计入）
       合计 = 65 = 41 + 24  ✓
副类:  D1 7  D2 9  D3 2  D4 6  D5 1  D6 4  D7 1  D9 1  D10 7  D11 3
       D12 1  D13 1  D14 6  D15 6  D16 3  D17 2          （副类合计 60）
```

---

## §2 T2 · 死因 × 逃生条件 × 候选机制

### §2.0 先说一个**元类 D0**（它解释了 D6 与 D10，故不计入 17 类）

`OUT-transplant.md`「Why no transplant — the genericity obstruction (sharpest line)」给出一条比 D6/D10 更狠的话：

> 每个已证的谱隙定理都是 "`∀` operator `L` satisfying hypothesis `H`: `spec(L)∩(0,g)=∅`"。
> 若 `H` **接纳**自由无质量场 ⇒ 它不能输出一个在那里为假的界；若 `H` **排除**自由场 ⇒ `H` 本身
> 就是一条非微扰的 YM 事实（`G1`/`G2` 的伪装）。

⇒ **凡"假设 ⇒ 无谱"型机制，都有一条二难：D6（前提空/空量化）与 D10（自由场同值）是同一二难的两支。**
这就是为什么 `D6+D10 = 11/65` 且没有任何族同时被两者当主因。**唯一的逃法是让机制不是二难型**：
输出必须是一个**数**（不是"有/无谱"的二分），且这个数在自由场上取**边界值**而在 YM 上取**内部值**——
即"输入是测度、输出是连续量"。语料里满足这条的**只有一个对象**：`A2` 的 `G″=S−D`（`D=0` 高斯 / `D>0` 纤维 Savvidy）。
这条观察是本文件对"为什么 150 条路同构"的最短解释。`[J]`

### §2.1 逃生表（17 行；**已出现 / 未出现** 逐格判）

| 死因类 | 逃生条件（精确到数学性质） | 候选机制（判） | 已出现? | 新候选? |
|---|---|---|---|---|
| **D1** 方向反 | 必须给 `inf` 的**下**界，且**不是**任何 `inf`-over-子空间的**对偶量**（否则可计算方向必然反向：造测试子空间 = 给上界） | (a) 形式比较/障碍函数（Sturm–Prüfer、最大原理）**已出现**（`DT-DIRICHLET` 的 floor 半，`[C]` 活）；(b) 反向 CLR/Bargmann 束缚态存在性 **已出现**（`HISTORICAL-TOOL-SWEEP` #3，`[NG]`）；(c) **n-宽度/逼近数（Kolmogorov, Karlin–Pinkus）**；(d) 指标/谱流计数 **已出现**（T15） | 3/4 | **C1**（n-宽度） |
| **D2** 循环 | 必须存在"**输入可测 / 输出为界**"的可判定分离（见证 vs 验证），且验证不调用 `inf supp`；或**改变问题本身** | (a) 区间算术/验证数值 **已出现**（T4 内）；(b) 证明挖掘/proof mining **已出现**（`O-MINIMAL-UNIFORMITY.md`，落在元数学轴）；(c) 符号级正性 FP/子椭圆 **已出现**（`HORMANDER.md`，**未入免疫图**）；(d) moment-SOS 对偶 **已出现**（`SOS-CERTIFICATE.md`：度 2 恰为 Popoviciu） | 4/4 | 无（该类的逃生只能是"换轴"，见 §4） |
| **D3** 尺度/regime | 核心不等式**不含 `β` 的幂**、不含小参数，且在 `β∈[1/48,2.5]` 与 `(L,a)` 双极限上闭包 | (a) 组合/簇展开 **已出现**（T7）；(b) **半经典 EBK/Bohr–Sommerfeld 量子化**；(c) 算术界（Selberg `3/16`、Deligne–Ramanujan）**已出现**（`OUT-transplant.md` §3：无算术格、无 Hecke 结构）；(d) 纯拓扑/指标 **已出现**（T8） | 3/4 | EBK（弱：= TR 的连续版） |
| **D4** 道/旋 | 同道（`0⁺⁺`、N-ality 0、规范不变）**且算子分辨**：能区分同一道里两个不同中性算子 | (a) Wilson 环形状导数 **已出现**（`reg` §B12）；(b) 道分辨指标 **已出现**（T15）；(c) **算子值 K-理论 / 完全正映射（Choi 秩）** | 2/3 | 弱候选 |
| **D5** 弱-* 盲 | **不是** `ν` 的弱-* 连续泛函 ⇒ 场空间值或算子族值（`PROFILE-FACET-DATAMODE` 已整类关闭谱类别） | (a) 场空间非线性不等式 **已出现**（`=H1`）；(b) 算子族一致性质 **已出现**；(c) **次因子 / L²-Betti / Jones 指标（von Neumann 维数论）**；(d) **迹方法 / 闭游走计数（Alon–Boppana 型）** | 2/4 | **C2**、**C3** |
| **D6** 空前提 | 前提必须在 YM 上**非空**、在自由场上**假**，且输出既不恒 0 也不误报正（⇒ 见 D0） | (a) 纤维 BL 缺陷 `D` **已出现**（A2，唯一过 D0 者）；(b) 三阶累积量/Stein 交换配对 **已出现**（T23）；(c) 非交换几何/Moyal 平面上的谱作用 **已出现**（`species/Z3-NCG.md`，判 acid-vacuous） | 3/3 | 无 |
| **D7** 常数/标度 | 不含 `β²` 负项 ⇒ 只能作用在**与 `β` 解耦**的量：纯组合/纯几何/纯拓扑/`β`-无关比值 | (a) 纯组合 `log(q−1)`/`d_N` **已出现**（A10）；(b) 纯几何 `W_geo≥N/12` **已出现**（N21 T0(c)）；(c) 纯拓扑 N-ality **已出现**（T2）；(d) `β`-无关比值 `ρ²`/`m₂/m₁` **已出现**（T39/N23） | 4/4 | **无：这一层已穷尽** `[J]` |
| **D8** 只 decoupling | 必须给"**耦合下**的下界"（溶解/重叠），而不是"消去/中性" | (a) 精确对易消去 **已出现**（N1）；(b) 避免交叉的定量化 **已出现**（无溶解结构定理 `[T]`） | 2/2 | 无 |
| **D9** FKG/符号 | 不依赖正关联的符号控制：负关联，或**代数**符号（显式 Weingarten 符号组合学） | (a) 反 FKG/负相关 **已出现且被否证**（N18）；(b) 符号交错求和/配对 **已出现**（登记为 `[M]`）；(c) Pfaffian/determinantal 结构 **已出现**（`FOSP-B.md`、`C6-DEPENDENCE-CLASS.md`） | 3/3 | 无 |
| **D10** 自由场不静默 | 输入必须依赖**非高斯性**（`g=0` 上取边界值），输出在自由场上静默 | (a) BL 缺陷 `D` **已出现**（A2）；(b) 高阶累积量 **已出现**；(c) 自由概率 **已出现**（`H1-T8-free-prob.md`） | 3/3 | 无 |
| **D12** 无尺度分离 | 不要求 `ℓ/L→0`：自相似/尺度无关重整化，或**图的粗化** | (a) 张量网络/MERA **已出现**；(b) Bethe 自相似 **已出现**（N20）；(c) 图/度量粗化（`RCD-SYNTHETIC-RICCI.md`、`WASSERSTEIN-OLLIVIER.md`）**已出现** | 3/3 | 无 |
| **D13** 载体恒零 | 载体不得在对称切片上退化：必须在**破缺对称的非零模切片**上定义 | (a) 加权框架（共形势）**已出现**（T35 [D]：加权下自伴 ⇒ 非正规性是规范假象）；(b) **非自治/时变转移算子** | 1/2 | 弱（按构造落回 D15） |
| **D14** 投射错位 | 必须读**重叠**而非仅谱：需要源分辨（bivariate）且**不是 `ν` 的泛函**的数据 | (a) 主角度/Kato–Temple **已出现**（T40）；(b) 双算子关联矩阵 + GEVP **已出现**（`GEVP-4D` 已补数据件）；(c) **frame/Parseval 界（非正交 frame 的谱界）** | 2/3 | 弱 |
| **D15** 定标不匹配 | 自带 `(L,a)` 一致的无量纲化，且体积强度正确（`χ` 而非 `χ/V`） | (a) 乘性残数搬运 `Z_O²` **已出现**（A1）；(b) step scaling/FSS **已出现**（T4/T6） | 2/2 | 无 |
| **D16** 类型错位 | 输出必须**实值且非局部常数**（非刚性）——**且**在阈值处有间隙（否则落进 D18） | (a) 谱移函数 `ξ(λ)` **已出现**（`HISTORICAL-SWEEP` #2）；(b) η-不变量 **已出现**（T15）；(c) 算子值/矩阵值不变量 **已出现**（`MICROLOCAL-BLIND.md` 用 operator-valued 分布） | 3/3 | 无 |
| **D17** 恒等式≠不等式 | 输出必须是有方向的**不等式**，且不是把恒等式两侧截断 | (a) 桥 `L2`/QI-1 **已出现**；(b) moment-SOS/Lasserre **已出现**（`SOS-CERTIFICATE.md`） | 2/2 | 无 |
| **D18**（准类）值域稠密 | 输出值域必须在阈值处**离散或有间隙**（如只取 `{0}∪[c,∞)`） | (a) II₁ 维数（像 `[0,1]` 稠密）**未出现** → 落 C2；(b) 指标（像 `ℤ`，有间隙）**已出现但刚性**（D16） | 1/2 | **C2 的死法** |

**三个可核计数.**

1. **未出现 = 7**：`n-宽度`（D1）、`EBK`（D3）、`算子值 K-理论`（D4）、`L²-Betti/Jones`（D5/D18）、
   `迹方法`（D5）、`非自治转移`（D13）、`frame/Parseval`（D14）。
2. **其中"弱" = 5**（EBK = TR 的连续版；算子 K-理论 / 非自治 / frame 三个按构造立刻落回 D15 或 D1 的邻域）。
   ⇒ **真正值得测的只有 2 个：`L²-Betti` 与 `迹方法`**——再**加**一个"名义未出现、实质与 min-max 同源"的 `n-宽度`。
   这三个就是 §3 的 C1/C2/C3。
3. **两行是"满格 4/4 已出现"**：`D7`（`β`-解耦量）与 `D2`（循环）。⇒ `D7` 这一层的**穷尽性**是
   本文件能给的**最强单项结论**：*"与 `β` 解耦的量"这一整层（组合/几何/拓扑/比值）在具名类内已无候选。* `[J]`

---

## §3 T3 · 三个"未出现"候选：判死线 → 方向账 → 五条件 → 酸测 → α/β/γ

> 顺序按 `LOOP-RUNBOOK.md` §4 第 3.5 步：**第一段判死线、第二段方向账**。三个都给了定义式（章程 §5：没公式的候选视为未交卷）。

### 3.1 C1 · n-宽度 / 逼近数（Kolmogorov, Karlin–Pinkus）

**对象.** `d_n(L) := inf{ sup_{‖x‖≤1} dist(Lx,V) : dim V = n }`，`a_n(L) := inf{ ‖L−R‖ : rank R < n }`；
经典 `a_n = s_n`（奇异值），`L≥0` 自伴时 `s_n = λ_n` ⇒ `a_n` 是**特征值的精确表征**。

**判死线（第一段）.** 若 `a_n` 的**唯一可执行**方向是"造一个 `n` 维子空间 / 一个秩 `<n` 的逼近"，则一切可算估计给的是 `a_n` 的**上**界 ⇒ 给 `λ_n` 的**上**界 ⇒ 在第一段就死（`D1`），不必跑数。**确认：命中**，且这是定义式的后果、不是巧合。

**方向账（第二段）.** 目标 = `λ₁` 的**下**界；宽度的定义是**两个 `inf`**（对子空间、对秩）⇒ 可执行方向 = 给 `inf` 的上界（造对象）⇒ `λ_n = a_n ≤ 造出来的值`。**方向反**——与 `RB-RATIO-BOUND` §0-0 的"对 `min` 取上界 = 造一个测试空间"是同一条推理，只换了名字。

**五条件 / 酸测.** 1 非微扰 PASS（无 `β` 展开、无小参数）· 2 控号下界 **FAIL**（方向反）· 3 非循环 PASS · 4/5 不可判（未到）⇒ **1/5 + 2 未判**；按"方向反了立刻停手"⇒ **不跑数**，酸测记 `N/A（方向先死）`。

**α/β/γ.** **无归属**：不是新坐标、也不是 H1 改名，是 `D1` 的一个**对偶面**。`[NG]`

**附带（C1 的唯一价值，可复用）.** **`D1` 在对偶下封闭**：任何"对某族子空间取 `inf`"的量——Kolmogorov 宽度、Ky Fan rank-`r` 迹（A7）、Courant–Fischer、交错压缩、Ritz/GEVP——都把**可计算方向**送到上界 ⇒ `D1` 不可能靠"再找一个变分量"逃掉；只能靠**非变分**不等式（形式比较/Bochner/指标），而那三者分别落进 `D2`/`D16`（`DT`、`HORMANDER`、`T8`）。`[T]`-级观察

### 3.2 C2 · 次因子 / L²-Betti / Jones 指标（von Neumann 维数论）

**对象.** `𝒩 := W*(L_{0⁺⁺}, 真空态)` 作用在 `H_phys = L²(𝒜/𝒢, μ̄)`；目标量 = 低谱投影 `P_λ := 1_{[0,λ]}(L_{0⁺⁺})` 的 `dim_𝒩 P_λ`（或 `L²`-Betti 数、Jones 指标 `[𝒩:ℳ]`）。要证的是**排除**命题 `P_λ = 0`。

**判死线（第一段）.** 两条，命中任一条即死：**(i)** 若 `𝒩` 是 type III₁（`X1-MODULAR.md`：local algebras type III₁、Connes spectrum `ℝ`、"gapless for all masses"），则**无忠实正规迹** ⇒ 维数论**无指称**（`D6` 的"无指称"面孔）；**(ii)** 若 `𝒩` 是 II₁，维数像 = `[0,1]` **连续稠密** ⇒ `dim ≠ 0 ⟺ P ≠ 0` 里**没有**"`dim<ε ⇒ P=0`"这一步 ⇒ 排除不可能（`D18`）。**确认：(i) 由 `X1` 的在册判词直接命中；(ii) 是本轮新发现的死法。两类穷尽了"有维数论的因子"⇒ 本候选无论落在哪一支都死。**

**方向账（第二段）.** 要证的是**存在性的否定**。维数论给"若 `P≠0` 则 `dim P>0`"（**存在性方向**）；"排除"需要一条**阈值间隙**。II₁ 的像在 `0` 处稠密、III₁ 没有像 ⇒ **不是"上界还是下界"的问题，而是"输出值域在判定阈值处没有间隙"**——17 类里没有这个形状 ⇒ 记为准类 **`D18`**。

**五条件 / 酸测.** 1 PASS · 2 **FAIL**（排除不了）· 3 PASS（不预设 gap）· 4 **FAIL**（自由场的 `L` 有从 `0` 起的连续谱 ⇒ `P_λ≠0` 对一切 `λ` ⇒ 与 YM **同值**，`D10` 同时命中）· 5 **FAIL** ⇒ **1/5**；薄原子/双井无分辨力，量纲虽无量纲却**不带 `g`**。

**α/β/γ.** **γ**：不是"H1 改名"，而是把 H1 的**可判定性**翻译进维数论后落进 `D18`，按章程 §5 属**决定性否证**。`[NG]` **新贡献 = `D18` 这条准类本身**（值域稠密 ⇒ 无法排除），与 `D7`（类常数）是否独立待裁决。

### 3.3 C3 · 迹方法 / 闭游走计数（Alon–Boppana 型）

**对象.** 设 `A` 为 `L_{0⁺⁺}` 的某离散化（或 `ν_L` 的矩）。经典手法 `Tr(A^k) = Σ_i λ_i^k ≥ d^k + λ₂^k ⇒ λ₂ ≥ (Tr(A^k) − d^k)^{1/k}`（`d`-正则，`λ₁=d` 已知）。**这是真下界，且不是变分**——这是它值得测的唯一原因。

**判死线（第一段）.** 若 `Tr(A^k)` 是 `ν` 的第 `k` 阶矩，则按 `DM`/`PROFILE-FACET-DATAMODE` 它是**弱-* 连续泛函** ⇒ 对 `inf supp` 结构性盲；只要输入是迹/矩，就在 `D5` 死。**确认：命中，且用的是语料自己那条"无隙射线"。**

**方向账（第二段）.** **PASS**：`Tr ≥ d^k + λ₂^k` 用的是"其余项为正"⇒ 给 `λ₂` 的**下**界，且**不经过** `inf`-over-子空间（逃 `D1`）⇒ 它是 7 个未出现候选里**唯一方向正确**的。但它逃不了 `D5`：`ROUND1-VERDICTS.md` §6 的**无隙射线** `ν_h=(1−h)ν+hδ₀` 上每个弱-* 连续数据随 `h` **线性移动**（`D₈^{L1}=0.183h`、`M₁→1.7819`、`χ→1.8370`），而地板指示 `1_{G_0.2}` 在 `h>0` 恒 `0`、`h=0` 为 `1` ⇒ **地板是射线端点处唯一的不连续量**；任何由矩造的、对数据连续的下界当 `h→0` 连续退回 ⇒ 给不出一致的正下界。

**五条件 / 酸测.** 1 PASS · 2 **形状 PASS / 对象错**（`λ₂` 是第二本征值，目标是 `μ₁=inf supp`，只在无限体积对应）· 3 PASS · 4 **FAIL**（自由连续谱也给矩，`bulk` 边到 `0` 的距离由带宽定、不由隙定 ⇒ **误报**，`D10`）· 5 **FAIL**（Alon–Boppana 要 `d→∞`，是 regime）⇒ **1.5/5**。

**α/β/γ.** **γ**（决定性否证）。附带一条**可复用的一般化**（C3 的真价值）：
> **矩-下界引理（一行，`[T]`）**：设 `𝔉` 是 `ν` 前 `N` 阶矩的连续泛函、`λ ≥ Φ(𝔉)` 单调。则存在 `ν_h→ν`（`h↓0`）使 `𝔉(ν_h)→𝔉(ν)` 而 `inf supp ν_h = 0`、`inf supp ν = m₁>0` ⇒ `Φ(𝔉(ν)) = lim Φ(𝔉(ν_h)) ≤ 0`。**⇒ 矩造不出正下界，不是"还没找到"，是"连续性上不可能"。** 与 `DM` 的 `d_N~e^{−0.4595N}` 独立同向（那里"加阶数指数无用"，这里"任何阶数都不行，只要 `Φ` 连续"）。

### 3.4 判决：三个全死，但"穷尽"这句话**不成立**

- **C1 → `D1`**（对偶封闭）· **C2 → `D18`（+`D6`,`D10`）** · **C3 → `D5`（+`D10`）**。三个死因**互不相同**，且全部死在**判死线/方向账**阶段（无一需要跑数）——这正是 §4-3.5 的正确用法。
- **但"在'死因类'这个轴上，可逃的机制也已被穷尽"本轮不成立**：未出现候选 **7** 个，只测了 **3** 个，**剩 4 个未测**——`EBK/Bohr–Sommerfeld`（D3）、`算子值 K-理论 / Choi 秩`（D4）、`非自治/时变转移算子`（D13）、`frame/Parseval 界`（D14）。
- **我能写的最强话（限定形式）**：在具名类 `{41 族 ∪ 24 新路 ∪ 7 个未出现候选}` 内，**已测的 65+3 条无一逃出 17 类**；且 4 个未测候选里每一个在**动手前**就能看到它落回邻类（`D15` 或 `D1` 邻域）——这是 `[J]`，**不是定理**；要让那句话成立，必须把那 4 条逐条跑完判死线。
- **最接近穷尽的一句（可写）**：`D7` 那一行是**满格 4/4 已出现** ⇒ ***"与 `β` 解耦的量"（纯组合/纯几何/纯拓扑/`β`-无关比值）这一整层，在具名类内已无候选***。`[J]`

---

## §4 T4 · 下一轮建议（10 行）

1. **第三个轴 = 等价性轴**（改问 `H1` 的**充分性**：`H1 ⟹ G`），不是量轴、不是元数学轴。
2. **不是量轴**：量轴在谱类别内已被 `DM`/`PROFILE` 关死（`∀` 谱数据 = 盲 XOR 循环）；它唯一剩下的形态就是场空间泛函 = `H1` ⇒ 量轴只是等价性轴里的一个坐标。
3. **是等价性轴**：`BN-RESIDUAL` §0 已把"唯一活口"从 `H1` 这个**对象**改成"`H1` 充分性"这个**命题**（`L0+L2` 只能读成"H1 必要（平凡）"或"H1 充分（即 Clay）"）；`PROFILE` 同样说"唯一真自由度是 `ν→μ`"。
4. **死因轴给的定量理由**：唯一同时逃过 `D0/D1/D5/D6/D10` **形状**滤波的是 `A2` 的 `G″=S−D`（**其一致化本身 `=H1`，不算逃生**），缺口 `D/S≈2%` 是**定量**的 ⇒ 下一轴的货币不是"再造机器"，而是"哪条蕴含能吃掉 2%"。
5. **不是元数学轴**：`O-MINIMAL-UNIFORMITY.md` 已试；地板是弱-* 不连续点（描述集合论事实），元数学只能出 no-go ⇒ 留作**副轴**（把"不可判"写成定理），不闭合 Clay。
6. **第一步（可执行）**：把 `H1 ⟹ G` 拆成 (i) `H1 ⟹ I`（`s₂ ≥ (Ē−w₀s₁)/(1−w₀)`，`SYN` §4(i)）与 (ii) `H1 ⟹ (CP)`（单相 + 0⁺⁺ Kato 保持），在 `(L,a)` 网格上测**哪一句先坏**——纯解析，不需新 MC。
7. **判死线**：若 (i)、(ii) 都与 `H1` **逐字等价**（无严格弱化）⇒ 本轴也退化为循环，记 `γ` 并回退。
8. **顺手活**：给 §3.4 的 4 个未测候选各写一段判死线（解析的、一行级），把"死因轴穷尽"这句话要么坐实、要么永久放弃。
9. **记账补漏**：把 `HORMANDER.md`、`MICROLOCAL-BUILD.md`、`SOS-CERTIFICATE.md` 三条已死物种**补进免疫图**（本文件只引不改）。
10. **纪律不变**：无 QED；"未找到"只写成"在具名类内未出现"；Clay 0%。

---

## §5 自否证条件（3 条）

(a) **某族的主类归属被判错** ⇒ 重读 65 行判词句，看"决定性"那半句是否引全 ⇒ 只推翻表 1 计数与那一行，不动 §2/§3。
(b) **某"未出现"候选其实已在语料**（grep 漏别名）⇒ 用 `Pinkus`、`s-number`、`von Neumann dimension`、`Atiyah conjecture`、`closed walk`、`non-autonomous` 重扫 ⇒ 该候选降为"已出现"，C1/C2/C3 之一作废。
(c) **`D18` 与 `D7` 不独立**（值域稠密只是"类常数"的另一说法）⇒ 找一条 `D7` 成员但值域不稠密、或 `D18` 成员但无常数可谈 ⇒ 撤 D18 行与 §3.2 的 γ 归属。

*（本文件只新建、不改任何其它文件；未 commit。No QED. Clay 0%。）*

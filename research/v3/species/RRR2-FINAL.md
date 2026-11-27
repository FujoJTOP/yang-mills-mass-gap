# RRR2-FINAL — `OS1` 最终陈述与 `𝒮′` 载体结账

> 轮 RRR2 · 子代理 `rrr2-final` · `D:\Dev\math`。写域：只新建本文件；既有文件一字未动；未 `commit`；零新 MC；零 Lean 编译。
> 口径：只 `[T]`/`[NG]`；不用百分数；每结论带「下一手」；`AGENTS §1` 四要素见 §4；`AGENTS §6` 的位置图由 `lead` 回填。
> 窗口：`L` 固定、`a→0`（环面）；`L→∞` = IN-6，另账。编号按题面（OS1 = 平移）。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | **OS1 最终陈述 = 一般 Polish 载体 + 三条前提（稠密／弱收敛／格点不变）⇒ 平移不变**；`hTight` 已删（由 ① 六步供） | **`[T]`** | §1；`OOO1-LEAN-OS1.md:13,32–38`；`NNN1-OS1.md:18–24` | 位置图把 OS1 前提写成三行，删 `hTight` |
| 0.2 | **承重件由四降三**：`平移等距` 可删 —— S3 只消费「`(h,x)↦τ_hx` 在紧 `S×K` 上一致连续」+ 积一致结构，`\|h_a−t\|≤2a` 是**数的界**，与 `d` 是否平移不变无关 | **`[T]`**（更正 `NNN1:22`） | §1 降格条；`HeineCantor.lean:47` | 台账 `NNN1:0.2` 的「四件」改「三件 + 载体是拓扑 `ℝ⁴`-空间」 |
| 0.3 | **载体：只有 `H^{-s}` 与「`𝒮′` 带强对偶拓扑」与 Polish 三件相容**；`𝒟′`、`C(K)′`、「可分子空间」三条全 `[NG]` | **`[T]`/`[NG]`** | §2 | 首选「换拓扑不换空间」（`𝒮′` 强对偶）；`H^{-s}` 作等价备份 |
| 0.4 | **换载体后 OS1 结论够 Clay 用**：`H^{-s}` 上的平移不变性经 `ι` 搬到 `𝒮′`，逐字无损耗 | **`[T]`** | §3 搬运引理 | OS1 从 open 账上清掉，写成定理 |
| 0.5 | **换载体不打开任何 Clay 开件**：`AX` 开义务仍 `{L3} ∧ {GAP}`；OS1 过不了酸测（无质量 GFF 亦满足）⇒ 闭条件、零物理内容 | **`[NG]`**（"载体换完更接近 Clay"为假） | §3；`G1-ASSEMBLE-100.md:5–6`；`NNN1:14` | 只留 `{L3,GAP}` 两行；载体问题降为记账 |

## §1 最终陈述（完整量词句）

**OS1（载体一般形式）。** 设 `E` 是 Polish 拓扑空间（第一可数、可分、完备），`τ: ℝ⁴×E → E` 连续（`E` 为拓扑 `ℝ⁴`-空间），`G_a := aℤ⁴` 经 `τ` 作用，`a ∈ (0,a₀]`；设 `(μ_a) ⊂ 𝒫(E)`、`μ ∈ 𝒫(E)`。若

1. **(P1) 群像稠密**：`cl(∪_a G_a) = ℝ⁴`，且取 `h_a := a⌊t/a⌋` 给显式模 `|h_a − t| ≤ 2a`（逐坐标）；
2. **(P2) 弱收敛**：`μ_a ⇒ μ`，即 `∀f ∈ C_b(E) : ∫f dμ_a → ∫f dμ`；
3. **(P3) 格点不变**：`∀a ∀h ∈ G_a : (τ_h)_*μ_a = μ_a`（可放宽为 `|⟨f∘τ_h⟩_{μ_a} − ⟨f⟩_{μ_a}| ≤ ε_a‖f‖_∞`，`ε_a→0`）；

则　`∀t ∈ ℝ⁴ : (τ_t)_*μ = μ`。

**`hTight` 不在前提里**：`IsTightMeasureSet (range μ_a)` 由 (P2) 在 Polish 载体上经 `isTightMeasureSet_range_of_tendsto` 六步组装得出（`Compact.lean:634` → `IsCompact.isClosed` → `closure_minimal` → `IsCompact.of_isClosed_subset:103` → `Prokhorov.lean:609` → `ProbabilityMeasure.t2Space:422`），`[T]`，0 假设。

**证明骨架（五步，逐件点名）**：S1 稠密 → 取 `h_a→t`；(P3) 把 `⟨f⟩_{μ_a}` 换成 `∫f∘τ_{h_a}dμ_a`；S3 取 tight 的 `K`，`S:={h_a}∪{t}` 紧，`(h,x)↦τ_hx` 在紧 `S×K` 上一致连续 ⇒ `ω_a := sup_{x∈K}|f(τ_{h_a}x) − f(τ_tx)| → 0`，余项 `2‖f‖_∞ μ_a(K^c)`；S4 两端极限拼出 `|⟨f⟩_μ − ⟨f⟩_{(τ_t)_*μ}| ≤ 2‖f‖_∞ ε` 对一切 `ε`；S5 `C_b` 分离 `𝒫(E)` ⇒ 等号。

**降格条（本轮新增，更正 `NNN1:22`）**：S3 的 `ω_a→0` 只用到 `(h,x)↦τ_hx` 的**联合连续 + 紧性**（`HeineCantor.lean:47` 只要 `[UniformSpace]`）：`d((h_a,x),(t,x)) = |h_a−t| ≤ 2a` 对 `x` 一致，因为它是**数的界**而非度量的界。⇒「平移等距」不是承重件；Mathlib 亦无平移不变度量构造（`BirkhoffKakutani|rInvariant` 全库 1 命中且无关，§4 ④）。

| 承重件 | 现在的位置 | 判 |
|---|---|---|
| 群像稠密 | (P1)；`NNN1:20`（模 `2a` 显式），`cl(∪_a aℤ⁴)=ℝ⁴` | `[T]` |
| 紧上一致连续 | S3；`HeineCantor.lean:47`，作用于紧 `S×K` | `[T]` |
| tight | 已从前提删除；`isTightMeasureSet_range_of_tendsto`（`OOO1:32–38`） | `[T]` |
| 弱收敛 + 分离 | (P2)/S4–S5；`FiniteMeasure.ext_of_forall_integral_eq`（`FiniteMeasure.lean:379`） | `[T]` |
| `hJoint`（移动试验函数） | 数学侧由「紧上一致连续 + tight」合成（S3–S4）；Lean 侧仍 ②b（~80 行，`OOO1:58–66`） | `[T]` / `[NG]`（分账） |

**骨架前提 = 三条**：`hConv`(P2) + `hInv`(P3) + `hJoint`；`hTight` 删除后只剩这三条（`os1_map_shift_eq_self`，`OOO1:96–101`）。

## §2 载体判（Polish 三件：可度量化／可分／完备 + `ℝ⁴`-作用联合连续）

| 候选 | (i) 可度量化 | (ii) 可分 | (iii) 完备 | `ℝ⁴`-作用 | 判 | 代价／下一手 |
|---|---|---|---|---|---|---|
| **(a)** `𝒮′→𝒟′(ℝ⁴)` | **✗** 严格 LF 的（弱-* 或强）对偶，非第一可数 | ✓ | ✓ | ✓ | **`[NG]`** | `𝒟′` 内测度类更宽、Minlos 更松；删 |
| **(b)** `H^{-s}(𝕋⁴)`, `s>1` | ✓ Hilbert | ✓ | ✓ | ✓ 酉（Fourier 乘子模 1）⇒ 联合连续 | **`[T]`** | 唯一税：`sup_a E‖φ‖²_{H^{-s}}<∞`（`s>1`：自由场 `[T]`；YM = UV 二点界，`UVS/L3` 物种）。`ℝ⁴` 版用 `H^{-s}(ℝ⁴)`；环面版只给 `t∈𝕋⁴` |
| **(b′)** `𝒮′` 带**强对偶**拓扑 | ✓ Fréchet（Birkhoff–Kakutani 给相容度量） | ✓ 有理格点 `δ` 的有理组合稠密（+ 序列一致，同下注） | ✓ Fréchet 完备 | ✓ `τ_t` = `𝒮` 自同构的转置；等度连续给联合连续 | **`[T]`**（最小改动） | Mathlib 自注逐字：*"we use the pointwise topology **for now** and not the strong topology … not barrelled, hence the uniform boundedness principle does not hold"*（`TemperedDistribution.lean:54–58`）；𝔖-收敛构造已在（`UniformConvergence.lean:21–24,63`）⇒ 现行 `[NG]` 是**缺省选择**，不是对象不相容 |
| **(c)** `C(K)′`, `K` 紧 | **✗** 弱-* 只在范数球上可度量化（Banach–Alaoglu + 可分），全空间非第一可数 | ✓ | ✓ | **✗** `K` 固定，`τ_t` 移出 `K` | **`[NG]`** | 两轴同时失败；删 |
| **(d)** 不换载体、只在可分子空间上证 | **✗** 无限维 `τ`-不变子空间上点态拓扑仍有不可数坐标 ⇒ 非第一可数 | — | — | ✓ | **`[NG]`** | 改走「可数分离族」只给 `⟨f_n⟩` 上的不变；补回全 `𝒮` 需 `Z` 连续 = OS0（`os2 REPORT:337–357`，`NNN1:0.4`）⇒ 成本搬到 OS0 |

**载体裁决**：`𝒮′` 的 `[NG]` 的唯一来源是**拓扑选择**（点态 = 序列够用、Polish 不够用），不是空间不相容；换到强对偶后三件全 `[T]`，且测试函数侧引理（`SchwartzMap.compSubConstCLM`，`SchwartzSpace/Basic.lean:1068`）全部保留。`H^{-s}` 与 `𝒮′` 强对偶在 Polish 三件上同判 `[T]`，差别只在：`H^{-s}` 多一条矩界、少一层「Borel 类重述」；`𝒮′` 强对偶多一层 Borel/`C_b` 重述、少一条矩界。**序列收敛两版相同**（Fréchet–Montel 强对偶上序列收敛 ⇔ 点态）⇒ (P2) 的序列内容不随载体变。

**Lean 侧**：(b′) 数学 `[T]`、Lean `[NG]`（`SecondCountableTopology`/`CompleteSpace`/`PseudoMetricSpace` 三个实例仍要写，构造件在 Mathlib）；(b) 定义件已存在 —— `TemperedDistribution.besselPotential (s : ℝ)`、`MemSobolev (s : ℝ) p`（`Analysis/Distribution/Sobolev.lean:71,149`）`[T]`，Hilbert 完备包络 `[NG]`。

## §3 够不够 Clay 用

**搬运引理 `[T]`。** `ι: H^{-s}(ℝ⁴) ↪ 𝒮′(ℝ⁴)` 连续线性（`f∈𝒮⊂H^s` ⇒ `|⟨φ,f⟩| ≤ ‖φ‖_{H^{-s}}‖f‖_{H^s}`），且 `ι∘κ_t = τ_t∘ι`（两边都是「试验函数平移」的对偶）。若 `(κ_t)_*μ = μ`，则
`(τ_t)_*(ι_*μ) = (τ_t∘ι)_*μ = (ι∘κ_t)_*μ = ι_*((κ_t)_*μ) = ι_*μ`。
故 `μ̃ := ι_*μ` 在 `𝒮′` 上 `ℝ⁴`-平移不变；Schwinger 函数（对 `f∈𝒮` 的矩）由推前保号不变 ⇒ **OS1 的结论逐字搬过来，无损耗**。反向不成立（`𝒮′` 上的不变性不给 `H^{-s}` 上的不变性），但 OS1 不需要反向。

**裁决两格。**
1. **OS1 这一格：`[T]`。** 换载体后 OS1 是定理（三前提 + 载体联合连续），且不依赖存在性（`(P2)` 是前提）——与 `NNN1` 的「条件于极限在手」一致。
2. **Clay 这一格：`[NG]`。** 换载体不动任何开件：`AX` 的开义务仍是 `{L3} ∧ {GAP}`（`G1-ASSEMBLE-100.md:5–6,35`），OS2=`SYM-4D`/`IN-2`、OS4=`G2` 照旧；且**新增一条** `H^{-s}` 矩界（属 `L3`/`UVS` 物种）。OS1 过不了酸测（无质量 GFF 满足，`NNN1:0.5`）⇒ 它是闭条件、零物理内容，载体换得再干净也不产生质量隙。

**下一手**：① 位置图把 OS1 记「`[T]`；载体 = `𝒮′` 强对偶（首选）／`H^{-s}`（等价，附矩界）」；② 公理表只留 `{L3, GAP}` 两行，OS1 从 open 账删除；③ `H^{-s}` 矩界回查 `UVS` 是否已含（若含 ⇒ 载体问题零新增前提）；④ 载体换定后再谈 `②b` 的 80 行 Lean（`hJoint` 降为定理）。

## §4 检索核验（`AGENTS §1` 四要素）

| 工具 | 状态 | 检索式（逐字） | 命中／失败 |
|---|---|---|---|
| `grep` | **OK**（有效响应，无 `provider failure`/`ENOENT`/空响应） | ①`OS1` @ `research` ②`强对偶\|F[ré]chet\|Montel\|Sobolev\|H^\{-s\}` @ `research/**/*.md` ③`StrongDual\|UniformConvergenceBounded\|strong dual\|Bornology.*Bounded\|WeakDual` @ mathlib ④`BirkhoffKakutani\|invariantMetric\|rInvariant` @ mathlib | ①**132** ②**142** ③**835**（显示 250）④**1**（`Algebra/Star/SelfAdjoint.lean:31`，无关 ⇒ 无平移不变度量构造） |
| `read` | **OK**（带行号，无截断） | `OOO1-LEAN-OS1.md`(170)、`NNN1-OS1.md`(62)、`AA4-CARRIER.md`(69)、`TemperedDistribution.lean:36–75`、`UniformConvergence.lean:15–64` | 有效；**0** `TOOL-FAIL` |
| `glob`/`pwsh` | **OK**（`exit=0`） | `Test-Path research/v3/species/RRR2-FINAL.md`；`Get-ChildItem research/v3/species -Filter *.md`；`Get-ChildItem mathlib -Recurse -Filter *Sobolev*` | 目标**不存在**（新建不覆盖）；`species/*.md` = **1019**；Sobolev 文件 **2**（`Distribution/Sobolev.lean`、`FunctionalSpaces/SobolevInequality.lean`） |
| 未用 | — | `web_*` 本轮未用（判据全在仓内 + Mathlib 源码闭合；外部文献本轮无新引用） | **0** |

**失败标记**：`TOOL-FAIL` **0**；所有「缺」均为**已执行检索 0/1 命中**，非工具失败。零 Lean 编译 ⇒ 本文件无机器验证断言，全部为源码级判读（`文件:行` 已标）。

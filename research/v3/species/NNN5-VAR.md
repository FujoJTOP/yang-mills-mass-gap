# NNN5-VAR — 变分形与阈值等价：①②③ 全 `[T]`；① 须限制到 `Ω⊥`，③ 的 ⟺ 只对指数阈值

> **2026-11-27 · NNN5 轮 · 只新建本文件**；既有文件零改动；未 commit；零新 MC；不编译 Lean。
> **口径**（承 `JJJ1-VAR.md:4`、`FFF1-UNITS.md:8`、`FFF4-LAM1.md:38`）：`T_a` 自伴、`0≤T_a≤1`、`T_a1=1`、`Ω:=ker(1−T_a)=ℝ1`（PF）、`λ₂:=sup σ(T_a|_{Ω⊥})`（紧时取到）、`Δ_a:=−(1/a)ln λ₂`（质量单位）、`g:=1−λ₂`；固定 `L<∞`、`a>0`、`0≤β<∞`、Wilson 周期。**不含** `(L,a)`-一致。
> 标签只用 `[T]`/`[NG]`。

## §0 判词

| # | 命题 | 判 | 下一手 |
|---|---|---|---|
| ① | `inf σ(1−T_a) = inf_{‖f‖=1}⟨f,(1−T_a)f⟩` | **`[T]`**；逐字在整个空间上两端皆 `0`（`f=1`，零信息）；有内容版须 `f⊥Ω`，此时 `=g` | 引用一律带 `⊥Ω`，§1.1 |
| ② | `⟨f,(1−T_a)f⟩ = ½∬(f−f′)²K_a` | **`[T]`**（`K_a` 对称、行和 1；核表示下逐字） | 形式化先建核层，§2 |
| ③ | `Δ_a≥δ ⟺ g≥1−e^{−aδ}`；`g≥aδ ⟹ Δ_a≥δ`，反向假 | **`[T]`**（⟺ 纯代数；充分用 `1−x≤e^{−x}`；反例 `a=δ=1, λ₂=e^{−1}`） | 记账只写指数阈值；线性阈值标“仅充分” |
| ★ | 三条能否无条件 `[T]` | **能**：不依赖 L3/GAP；只改写、不产正下界 | 下界仍归 Dirichlet 形（H1），§3 |

**一句话**：`Δ_a≥δ` 的精确变分书写是 **`∀f⊥Ω: E_a(f)≥1−e^{−aδ}`**，不是“谱隙 `≥aδ`”。

## §1 三条证明

### 1.1 ① 谱定理
`A:=(1−T_a)|_{Ω⊥}` 自伴有界（`Ω⊥` 不变：`h∈Ω ⇒ ⟨Tf,h⟩=⟨f,Th⟩=⟨f,h⟩=0`），`m:=inf σ(A)`（`σ` 紧 ⇒ `m∈σ`）。谱定理给 `A−m≥0` ⇒ `⟨f,Af⟩≥m‖f‖²` ⇒ `inf_form≥m`。若 `c:=inf_form>m`，则 `A−c≥0`；正元 ⇒ `σ(A−c)⊂[0,∞)`（`ε>0` 时 `A−c+ε≥ε` ⇒ 下有界 ⇒ 单射闭值；自伴 ⇒ 稠值 ⇒ 可逆），故 `σ(A)⊂[c,∞)`，与 `m<c` 矛盾。⇒ `inf σ(A)=inf_form`。又 `(1−T)−(1−λ)=−(T−λ)` ⇒ `σ((1−T)|_{Ω⊥})=1−σ(T|_{Ω⊥})` ⇒ `inf(1−S)=1−sup σ(S)=1−λ₂=g`（`S:=T|_{Ω⊥}`）。**`H` 上版**：`1∈ker(1−T_a)` ⇒ `0∈σ`；`1−T_a≥0` ⇒ `inf σ=0=⟨1,(1−T_a)1⟩`。∎

### 1.2 ② 代数
`(T_af)(x)=∫K_a(x,y)f(y)dμ(y)`（`f` 实值），`K_a(x,y)=K_a(y,x)`、`∫K_a(x,y)dμ(y)=1`（`T_a1=1`）。则
`⟨f,(1−T_a)f⟩=∫f²−∬K_a f(x)f(y)`
`=½∬K_a[f(x)²+f(y)²−2f(x)f(y)]`（第 1 项用行和；第 2 项用对称 ⇒ 列和）
`=½∬(f(x)−f(y))²K_a(x,y)≥0`；`K_a>0` ⇒ 等号 ⟺ `f` 常数 ⇒ `E_a` 在 `Ω⊥` 上严格正（逐点，不主张一致下界）。∎

### 1.3 ③ 两向
`λ:=λ₂∈(0,1]`、`a>0`：`Δ≥δ ⟺ lnλ≤−aδ ⟺ λ≤e^{−aδ} ⟺ g≥1−e^{−aδ}`。
充分：`g≥aδ ⇒ λ≤1−aδ`；`x:=aδ`，`1−x≤e^{−x}`（`add_one_le_exp(−x)`）⇒ `λ≤e^{−aδ} ⇒ Δ≥δ`（`x>1` 时前件假；`x=1` 只 `λ=0` 可满足，`Δ=+∞`）。
反向假：`a=δ=1`、`λ=e^{−1}` ⇒ `g=0.6321205588<1=aδ`，而 `Δ=1=δ`。∎

## §2 Lean 化

**判**：三条都可形式化；但本仓 `lean/` 只有 OS 公理陈述层（`Basic.lean`：`SatisfiesAllOS`/`MassGapStatement`），**无 `T_a`、`λ₂`、核** ⇒ 需新层；Mathlib 有算子零件、无格点对象。下表路径省略前缀 `Mathlib/`。

| 需要 | Mathlib（master，2026-11-27 抓取） | 状态 |
|---|---|---|
| `spectralRadius` | `Analysis/Normed/Algebra/Spectrum.lean:78`；`spectralRadius_eq_of_unital:259`；`spectralRadius_le_eNorm:151`；`ContinuousLinearMap.spectralRadius_eq_eNorm`（`Rayleigh.lean:182`，自伴 ⇒ 谱半径=范数） | 有 |
| `ContinuousLinearMap` / 伴随 | `Analysis/Normed/Operator/ContinuousLinearMap.lean`；`InnerProductSpace/Adjoint.lean`（`adjoint`、`isSelfAdjoint_iff':276`、`IsSelfAdjoint.adjoint_eq:362`、`IsSymmetric.isSelfAdjoint:389`） | 有 |
| 紧自伴谱定理 | `InnerProductSpace/Spectrum.lean:435` `orthogonalComplement_iSup_eigenspaces_eq_bot`；`:455` `finite_dimensional_eigenspace`（均需 `IsCompactOperator`） | 有（紧版） |
| 变分/瑞利 | `InnerProductSpace/Rayleigh.lean:56` `ContinuousLinearMap.rayleighQuotient`；`:120` `norm_eq_iSup_rayleighQuotient`；`:297` `IsSelfAdjoint.hasEigenvector_of_isMinOn`；`:105` `iInf_rayleigh_eq_iInf_rayleigh_sphere` | 有 |
| 有限维 min–max | 同上 `:346` `LinearMap.IsSymmetric.hasEigenvalue_iInf_of_finiteDimensional` | 有（仅有限维） |
| 一般 `inf σ=inf 形式` | 已下载 13 文件内检索 `sInf.*spectrum.*rayleigh`：0 命中（非全库否定） | **无单条**（可拼） |
| 紧算子 min/max 达到 | `Rayleigh.lean:32–36` 明写 **TODO** | **未形式化** |
| CFC / `log` | `Analysis/CStarAlgebra/ContinuousFunctionalCalculus/Basic.lean`；`SpecialFunctions/ContinuousFunctionalCalculus/ExpLog/Basic.lean`；`StarOrder.lean:59` | 有（①不用；`H_a` 用） |
| 对数/指数 | `Real.log_le_iff_le_exp`（`Log/Basic.lean:163`）；`Real.le_log_iff_exp_le:167`；`add_one_le_exp`（`Log/Basic.lean:85`、`Exp.lean:215` 使用） | 有 |
| Fubini | `MeasureTheory/Integral/Prod.lean:481` `integral_integral_swap` | 有 |
| Dirichlet 形 / L² 积分算子 | 全树**路径**检索：`Dirichlet` = 10 文件（全数论）、`IntegralOperator` = 0；`Probability/Kernel` 有 Markov 核 | **路径级无**（需自建核层；未做全库声明检索） |

**代价（估）**：① 80–150 行（子空间 + `sSup/sInf` 换序 + 正元谱非负；字面版 10 行但零信息）；② 150–300 行（`L²(μ)` 核算子、对称、Fubini、代数）；③ 40–80 行（`Real.log 0=0` 是坑：`λ₂=0` 单列或 `WithTop`）。

**AGENTS §1 四要素**：`pwsh`+`Invoke-WebRequest` **OK**（HTTP 200；树 `2.86 MB` + 13 份源文件；5 个猜测路径 404，均未用于判词）；`grep` **OK**（行号正文）；检索式 `spectralRadius`、`IsSelfAdjoint`、`rayleigh`、`sInf.*spectrum.*rayleigh`、`Dirichlet|IntegralOperator`、`log_le_iff_le_exp|add_one_le_exp`、`singleton_sub_eq`；命中见表；`Dirichlet 形/积分算子`=0（路径级；工具 OK，非 TOOL-FAIL）；`singleton_sub_eq` 仅见 2 处使用、声明处未定位（工具 OK，非 TOOL-FAIL）。

## §3 裁决

**能**：三条都无条件 `[T]`，不依赖 L3/GAP。① 只用自伴 + `Ω⊥` 不变；② 只用核表示与行和（核存在=HS，`FFF4-LAM1.md:38`）；③ 只用 `log/exp` 单调性。三者对 `(L,a,β)` 逐点成立；不需一致常数、不需 `λ₂<1` 的一致下界、不需连续极限。

**但**：① 把 `g` 写成 `inf`（试态只给上界 ⇒ `Δ` 的上界）；② 把同一 `inf` 写成 Dirichlet 形（仍 `inf`）；③ 只翻译下界。要 `Δ≥δ` 仍须 `∀f⊥Ω: E_a(f)≥1−e^{−aδ}`，即 Dirichlet 形的**一致正下界**（H1/Poincaré）= 靶。⇒ 本轮是记账/纠错：③ 的 `[T]` 只对指数阈值，`JJJ1-VAR` 的 Q2 `[NG]`（`≥aδ` 版）保留，不被提升；亦不产生隙。

**下一手**：① `lean/YMStatement/Transfer/Var.lean` 先形式化 `Ω⊥` 版（`norm_eq_iSup_rayleighQuotient`+`spectralRadius_eq_eNorm`，绕开缺的 `inf σ=inf 形式`）；② 用 `Prod.lean` 的 Fubini 建核层，先证 `K_a≥0` 版（不证 `K_a>0`，免极值原理）；③ 把 `g≥1−exp(−aδ)↔Δ≥δ` 写成引理，`g≥aδ` 只作 `sufficient`。

*状态：① `[T]`（须 `⊥Ω`；`H` 上版零信息）；② `[T]`（核表示下逐字）；③ `[T]`（指数阈值 ⟺；线性阈值单向+反例）。三条均不依赖 L3/GAP，亦不产生隙下界。零新 MC；未编译；既有文件零改动；未 commit。*

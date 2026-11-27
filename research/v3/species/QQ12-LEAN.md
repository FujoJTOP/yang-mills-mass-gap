# QQ12-LEAN — 三格各一条 `axiom`（类型齐 · 与 `T14`/`T15` 同址）；判：三格全 `[T]` 时 `T16` **不是**只需改名（三参只直吃一格）

> **轮** 244 · 子代理 `agent-89` · `D:\Dev\math`。**写域**：只新建本文件（≤120 行）；既有文件一字未改；未 `git commit`；**零新 MC**；**不编译**（题给）；`docs/position-current.html` 交 `lead`（`AGENTS §6`）。
> **靶（逐字）**：写三条 Lean 声明（A/B/C 各一条 `axiom`），类型签名齐、与 `T14`/`T15` 的消费形式匹配；并判「若三格全 `[T]`，`T16` 是否只需改名」。
> **口径**：三格 = 本轮读法 `A: ϖ≤C`（`KK2:25`）· `B: ϖ≥c₀`（`KK2:27`）· `C: w₁≥w₀`（`KK2:29`）；`T16` 三参 = `IotaL2Seam`/`L3Bound`/`BFloor`（`HY6:101`）。`N` 固定（`KK2:24`）；`w₀=1/(4N⁴)`（`KK1:50`）。
> **标签**：只 `[T]`/`[NG]`；无百分数；每条带「下一手」。直读 15 份（§4）。

## §0 判词

| # | 断言 | 判 | 一行依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | 三条 `axiom` 写得成、类型齐（变量／量词／域三级核过）；缺的只是内容 | **`[T]`**(形式)/**`[NG]`**(内容) | §1；同 `HV5:67` 三级核 | 落 D0 后再过一次 elaboration |
| **0.2** | **A 与 `hA` 槽逐字同型**：`qq12_A_to_L3Bound` 三步（`obtain`＋特化＋`exact`），零新内容 | **`[T]`** | §1；同构于 `KK1:39–40` 的 `L3Bound.weaken` | 台账 `hA` 行写「A 直通」，不写「A 推 L3」 |
| **0.3** | **C 与 `hB` 槽不同型**：`W1Floor`（`w₁`，β-尾，`L≥2`）≠ `BFloor`（`s₁`，`c` 自由）⇒ 必有一条 4 行桥（`c:=w₀`） | **`[T]`** | §1 `qq12_C_to_BFloor`；`KK1:0.4,62–70` | 落桥 body 4 行（`div_le_iff₀`→`le_trans`→`field_simp`） |
| **0.4** | ★ **B 在 `T16` 里无槽**：三参没有 `ϖ≥c₀`；B 在装配里只出现于 `A∧B⇒(Z-last)` | **`[T]`**（结构） | §2；`HY6:101`；`KK2:76` | 二选一：`T16` 加第 4 参／把 B 显式写进 `+OS` 与纸面框 |
| **0.5** | **`hC`（`IotaL2Seam`）不在三格内**：它是第 4 条 `axiom`（`HV5:40`）；「可归约到 `L3`」的 4 行未落 | **`[NG]`**（内容） | §2；`HZX:0.4`、`KK2:56` 均记为下一手 | 写 `L3⇒memLp` 与 `L3⇒strong` 各两行 |
| **0.6** | ★ **判：不是改名**（两读法皆否，§3）：三格闭后仍缺 1 桥＋1 seam＋1 槽＋1 OS；若「三格」指 `T16` 三参本身，则要做的是 `axiom→theorem` 并让 `#print axioms` 归零 | **`[T]`**（判词） | §3(a)(b)(c) | 先做 0.4 的槽决定（改签名会动 §3 全部计数） |
| **0.7** | **可数的量**：`#print axioms hy6_clay_of_three` 现状 **3** 条（`KK1:0.5`）；**只改名 ⇒ 仍 3 条**；三格＋桥＋seam 真闭 ⇒ **0** 条，但 `ClayStatement` 仍相对纸面框（`KK2:6.1` 五项） | **`[T]`**（记账） | §3(c) | 把 3→0 分两段：先 seam，后框 |

**一句话。** 三条签名都不缺类型：A 已是 `T15` 的 `hA`（只差一次特化），C 差一条 4 行桥（`c:=w₀`），B **没有槽**；`hC` 槽的供给是一条三格之外的 seam。故「三格全 `[T]`」既不自动给 `T16` 一个无条件证明，也不只是改名——要改的是**依赖表**（3→0）与**签名**（3 参→4 参），名字反倒可以一字不动。

## §1 三个声明（`axiom` 三条 ＋ 槽定义与靶常数 ＋ 两条桥；**未编译**）

```lean
-- QQ12-LEAN 声明稿（signature-only；`theorem` 不带 body，同 `KK1:§1` 规；**未编译**）
import Mathlib.Data.Real.Basic
import Mathlib.MeasureTheory.Function.LpSpace.Basic
set_option autoImplicit false
namespace MassGap.QQ12
/- D0 变量面：逐字 `HV5:24–28`；`IotaL2Seam` 的权威定义在 `HV5:34–38`（本文件不重述） -/
variable {Ω H Obs : Type*} [MeasurableSpace Ω]
variable [NormedAddCommGroup H] [InnerProductSpace ℂ H] [CompleteSpace H]
variable (μ : Measure Ω) [IsProbabilityMeasure μ]
variable (ι : ℝ → Obs → Ω → H)
variable (Ha : ℝ → H →L[ℂ] H) (Hlim : H →L[ℂ] H)
/- 参数面：`N` 固定（`KK2:24`）；`Par` := `(β,L,b.c.,k)`，`b̂` 取固定物理块（`HZ6:62`） -/
variable (Par : Type*) (varpi : Par → ℝ → ℝ) (a0 : ℝ)
/-- A 格：`∃C>4 ∀p∀a∈(0,a₀], ϖ(p,a)≤C`（逐字 `KK2:25`；`∃C` 在全部参数前 = 全部内容）。 -/
def VarpiUpper (w : Par → ℝ → ℝ) (b : ℝ) : Prop :=
  ∃ C : ℝ, 4 < C ∧ ∀ p : Par, ∀ a : ℝ, a ∈ Set.Ioc 0 b → w p a ≤ C
/-- ★ A 格 axiom：内容 `[NG]`（缺恰一个 `N_B` 的幂，`II1-C4:4.1`）；类型 `[T]`。 -/
axiom qq12_A_varpi_upper (ha0 : 0 < a0) : VarpiUpper varpi a0
/-- `T15`/`T16` 的 `hA` 槽（逐字 `KK1:37`）。 -/
def L3Bound (v : ℝ → ℝ) (b : ℝ) : Prop :=
  ∃ C : ℝ, 4 < C ∧ ∀ a : ℝ, a ∈ Set.Ioc 0 b → v a ≤ C
/-- 桥 A→`hA`：`∃C` 的 `∀p` 特化（三行 body，无新内容）。 -/
theorem qq12_A_to_L3Bound (hA : VarpiUpper varpi a0) (p : Par) : L3Bound (varpi p) a0
/-- B 格：`∃c₀>0 ∀p∀a∈(0,a₀], c₀≤ϖ`（逐字 `KK2:27`，值下界／防退化）。 -/
def VarpiLower (w : Par → ℝ → ℝ) (b : ℝ) : Prop :=
  ∃ c0 : ℝ, 0 < c0 ∧ ∀ p : Par, ∀ a : ℝ, a ∈ Set.Ioc 0 b → c0 ≤ w p a
/-- ★ B 格 axiom：`[T]` 只到 `b̂=2a` 的 `v(β)/4=1.6048e−3`（`NN1:0.4`）；`∀b̂` 缺口 `N_B/k_max=4`（`NN1:0.9`）。 -/
axiom qq12_B_varpi_lower (ha0 : 0 < a0) : VarpiLower varpi a0
/-- `T14` 的 `hB` 槽（逐字 `HY6:100`；**不是** `BFloor`）。 -/
def DeltaFloor (Delta : ℝ → ℝ) (δB : ℝ) : Prop := ∀ a : ℝ, δB ≤ Delta a
/-- 靶常数 `w₀ = 1/(4N⁴)`（`KK1:50`；`N=2` 给 `1/64`）。 -/
noncomputable def w0 (N : ℕ) : ℝ := 1 / (4 * (N : ℝ) ^ 4)
/-- C 格：`∃β₀>0 ∀β≥β₀ ∀L≥2, w₀≤w₁`（逐字 `KK1:54–55` 的 `R5`；`w₁≡ε`，`∀L≡inf_L`）。 -/
def W1Floor (w1 : ℝ → ℕ → ℝ) (N : ℕ) : Prop :=
  ∃ β0 : ℝ, 0 < β0 ∧ ∀ β : ℝ, β0 ≤ β → ∀ L : ℕ, 2 ≤ L → w0 N ≤ w1 β L
/-- ★ C 格 axiom：等价塌缩形 `r := Σ_{n≥2}|E⁰ⁿ|²/|E⁰¹|² ≤ 4N⁴ = 64`（`NN12:0.3`）。 -/
axiom qq12_C_w1_floor (w1 : ℝ → ℕ → ℝ) (N : ℕ) (hN : 2 ≤ N) : W1Floor w1 N
/-- `T15`/`T16` 的 `hB` 槽（逐字 `KK1:58–60`）。 -/
def BFloor (s1 chi M0 : ℕ → ℝ → ℝ) : Prop :=
  ∃ βA : ℝ, ∃ c : ℝ, 0 < c ∧ ∀ L : ℕ, 2 ≤ L → ∀ β : ℝ, βA ≤ β →
    s1 L β ≥ c * (M0 L β / chi L β)
/-- 桥 C→`hB`：`c:=w₀`；前件 `w₁≡f₁²/M₀` ＋ 夹逼左腿 `f₁²≤χs₁`；body 4 行（`KK1:62–70`）。 -/
theorem qq12_C_to_BFloor
    (w1 : ℝ → ℕ → ℝ) (f1sq s1 chi M0 : ℕ → ℝ → ℝ) (N : ℕ) (hN : 2 ≤ N)
    (hχ : ∀ (L : ℕ) (β : ℝ), 0 < chi L β) (hM : ∀ (L : ℕ) (β : ℝ), 0 < M0 L β)
    (hW : ∀ (L : ℕ) (β : ℝ), w1 β L = f1sq L β / M0 L β)
    (hLow : ∀ (L : ℕ) (β : ℝ), f1sq L β ≤ chi L β * s1 L β)
    (hC : W1Floor w1 N) : BFloor s1 chi M0
end MassGap.QQ12
```

**不在三格内的一件**：`hy6_N1_iota_l2_seam : IotaL2Seam μ ι Ha Hlim`（`HV5:40` 逐字）—— `hC` 槽的供给（0.5）。

## §2 槽位匹配（`T16` 吃三参，本轮三格供三题）

| 消费点 | 逐字签名（`HY6:99–101`） | 供给链 | 判 | 下一手 |
|---|---|---|---|---|
| **T14** `hy6_gap_positive` | `(hres : StrongResolvent Ha H) (hB : ∀ a, Delta a ≥ δB) : DeltaCont ≥ δB` | `hres` ← `IotaL2Seam.strong`（`HV5:38`）；`hB` ← **无直供**：须 `C→BFloor`＋`BFloor→(B)→Δ`（`T13`），且 `δB(c)` 未写 | `[T]`(槽位)/`[NG]`(`δB`) | 写 `δB=δB(w₀)` 一行并接 `T13` |
| **T15** `hy6_unique_limit` | `(hC : IotaL2Seam) (hA : L3Bound varpi a0) (hB : BFloor s1 chi M0)` | `hA` ← **A 直通**；`hB` ← C 经桥；`hC` ← 第 4 条 axiom | `[T]`(2/3)/`[NG]`(`hC`) | seam 的 4 行（0.5） |
| **T16** `hy6_clay_of_three` | `(hC) (hA) (hB) : ClayStatement`；`T16 ← {T14,T15}+OS`（`HY6:61`） | 同上；**B 无槽**；OS 在框外 | `[T]`(形式) | §3 |

**计数（可复算）**：三格→`T16` 三参，**直通 1/3**（A↔`hA`）；**桥 1**（C↔`hB`）；**无槽 1**（B）；**框外 1**（OS）；**第 4 条 axiom 1**（seam）。⇒「三格三条 axiom」与「`T16` 三参」不是同一个三元组。

## §3 判：三格全 `[T]` ⇒ `T16` 是否只需改名

**(a) 两套三元组不重合。** 本轮三格 = `{VarpiUpper, VarpiLower, W1Floor}`（`KK2:25–29`）；`T16` 三参 = `{IotaL2Seam, L3Bound, BFloor}`（`HY6:101`）。交集只有「上界」一项——「三格全 `[T]`」在本轮口径下**不**蕴含三参全 `[T]`。

**(b) 最小残余 = 四件，逐件都不是改名。** ①`qq12_C_to_BFloor` 的 4 行 body——真证明（`[T]` 可写），改名替不掉；②seam：`hC` 槽的 `IotaL2Seam`——第 4 条 `axiom`（`HV5:40`），4 行归约未落（`HZX:0.4`、`KK2:56`）；③B 的槽：三参里没有 `ϖ≥c₀`，而它在 `A∧B⇒(Z-last)`（`KK2:76`）——要么给 `T16` 加第 4 参（**改签名**），要么把 B 写进 `+OS`／纸面框（`HY6:61`、`KK2:6.1`）；④OS 本体（`T16←{T14,T15}+OS`）。

**(c) 只改名的后果（可数）。** `#print axioms hy6_clay_of_three` 现状列 3 条（`N1/N2/N4`，`KK1:0.5`）：**只改名 ⇒ 仍列 3 条**——名字变了，依赖表逐字不变。反之，若「三格」指 `T16` 三参本身，则要做的是三条 `axiom` → 三条 `theorem`（body 填入）并让该表归 0；此时 `T16` 的签名与名字都可一字不动，**变的恰是「是不是定理」这件事**。⇒「改名」与「闭」互不替代：改名 0 成本且不产生定理；闭必须付 (b) 的四件。

**(d) 剩余框。** 即便 (b) 全付，`ClayStatement` 仍相对纸面框（`Adm`＋Hyp 1＋FRD＋`(S1)-(S2)`＋KP，`KK2:6.1`）⇒ 无条件 Clay 仍 `[NG]`。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际调用／检索式（逐字） | 命中 | 失败 |
|---|---|---|---|---|---|
| 1 | `pwsh` | **OK**（exit 0，正文非空） | `Get-ChildItem research\v3\species -Filter *.md`；`Get-ChildItem -Recurse -Filter *.lean`（排除 `vendor`）；`(Get-Content …).Count` | 目录全列；`.lean` 全表（截断，spill 落盘）；行数实测 | 0 |
| 2 | `grep` | **OK**（无 provider failure／ENOENT／空响应） | ①`hy6_N1_iota_l2_seam\|hy6_N4_b_floor` @ `D:\Dev\math` ②`T16\|hy6_clay_of_three\|ClayStatement` @ `species/*.md` ③自检 `axiom qq12_\|theorem qq12_\|def (Varpi\|W1Floor\|L3Bound\|BFloor\|DeltaFloor\|w0)` @ 本文件 | ①**18**（承重 `HV5:40,51`、`KK1:16`、`HZX:14`、`QQ8:22,81`）②**27**（承重 `HY6:51,61,101`、`KK1:78,88`、`KK2:11,31`、`HZM:49`）③**12**（3 axiom＋2 theorem＋7 def） | 0 |
| 3 | `read` | **OK**（带行号正文，15/15 非空） | `HV5-LEAN`(97)·`KK1-LEAN`(109)·`HY6-ASSEMBLE`(130)·`KK2-FINAL`(109)·`NN1-BPROOF`(87)·`OO2-ALIGN`(83)·`PP1-LATRP`(110)·`OO6-FORMAL:1–45`·`NN11-KAPPA:1–45`·`II2-OB`(102)·`II1-C4`(94)·`HZY-C0`(94)·`NN12-S`(96)·`IM1-MEASURE:1–60`·`OO9-ALT:1–50` | 15 份 | 0 |
| 4 | 失败标记／`MC`·Lean·`git`·`web_*`·`skill` | **TOOL-FAIL 0；空响应 0；PATH-MISS 0；provider failure 0**；工具栈余项 **0**（题给「不编译」；catalog 无匹配技能） | — | — | 0 |

**直读 vs 转引（`AGENTS §3`）**：直读 = #3 十五份。只经命中行或他件转引、未开正文：`HZV-R5:0.1–0.5`（经 `KK1:14,62`）·`HZX-COMPLETE:0.4`（经 `KK2:56`）·`IAE-M:0.1`（经 `KK2:48`）·`HZF-COV`（经 `NN1`/`II2`/`OO2`）·`HZ6:62`（经 `KK2:24`）⇒ 记转引级。

## §5 边界与下一手

1. **未编译**（题给）：实例面（`H →L[ℂ] H` 的 `Algebra ℂ`、`MemLp` 的 `2:ℝ≥0∞` 强制、`IsProbabilityMeasure`）按 `HV5:93` 同款写法，未过 elaboration；不据此宣称新机检。
2. **`Par` 是抽象参数空间**：`(β,L,b.c.,k)` 的积未展开；`N` 固定、`N`-一致未证（`KK2:5,105`）⇒ 三条 axiom 都是「固定 `N`」形，不是 `N`-一致。
3. **内容全 `[NG]`**：A = 一个 `N_B` 的幂（`II1-C4:4.1`）；B = 块混合亏缺 `N_B/k_max=4`（`NN1:0.9`）；C = 静态比 `r≤4N⁴=64`（`NN12:0.3`）＋渐近起点 `1/t_A=1.28`（`NN12:0.5`）。本轮不产新数，引用数逐条带出处；零新 MC、零改既有文件、未 `commit`、未动 `docs/position-current.html`。
4. **两条下一手（按可执行性）**：① 先做 0.4 的槽决定——`T16` 若加第 4 参，本文件的三格计数与 `#print axioms` 记账一起改；② 写 seam 的 4 行（`L3⇒memLp`、`L3⇒strong`），它是 `hC` 槽的唯一出入口。

*QQ12-LEAN · 只新建本文件；未编译；三格各一条 `axiom`（`VarpiUpper`/`VarpiLower`/`W1Floor`）＋两桥（A→`hA`、C→`hB`）；判：三格全 `[T]` 时 `T16` **不是**只需改名——直通 1/3 槽、桥 1、无槽 1（B）、框外 1（OS）、第 4 条 axiom 1（seam）；只改名 ⇒ `#print axioms` 仍列 3 条。*

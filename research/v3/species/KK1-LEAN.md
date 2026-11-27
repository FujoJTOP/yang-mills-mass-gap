# KK1-LEAN — 两条剩余 `axiom`：柱 1 的 `∃C>4, ϖ≤C` 与柱 2 的 `R5`（类型齐 · 与 `T14`/`T15` 同址）

> **轮** KK1 · 子代理 `kk1-lean` · `D:\Dev\math`。**写域**：只新建本文件（≤110 行）；既有文件一字未改；未 `git commit`；**零 MC**；**未编译 Lean**；`docs/position-current.html` 交 `lead`（`AGENTS §6`）。进展单位 = 轮。
> **靶（逐字）**：把「存在 `C>4` 使 `varpi ≤ C`」与 `R5` 各写成一条 `axiom` —— 类型签名齐、与 `T14`/`T15` 的消费形式匹配。
> **直读（行号级）**：`HY6-ASSEMBLE`(130)·`HV5-LEAN`(97)·`HZV-R5`(88)·`IAE-M`(71)·`HZU-M`(88)·`HZG-BLOCK`(91)·`HZE-KIN`(84)·`HZP-EPS`(77)·`HZS-FINAL2`(103)·`HZX-COMPLETE`(107)·`Reduction.lean:88–132`·`L3B-Assembly.lean`(133，只读格式)。`HZF-COV` 只经命中行（`:18,66,72`）。
> **标签**：只 `[T]`/`[NG]`；无百分数；每结论带「下一手」。

## §0 判词

| # | 断言 | 判 | 一行依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | 两条 `axiom` 写得成、类型齐（变量 / 量词 / 域三级核过）：柱 1 = `L3Bound varpi a0`（`ha0 : 0<a0` 防空窗）；柱 2 = `R5 eps N`（`2≤N`，内层 `2≤L`）。缺的是内容不是类型 | **`[T]`**（形式）/ **`[NG]`**（内容） | §1；`HY6:43,46,49–51`；`HV5:67` | 先落 D0 再一次 elaboration（本文件不编译） |
| **0.2** | `4 < C` 进类型**不改**消费面：`L3Bound.weaken` 逐字给 `HY6:88` 的弱形 `∃C ∀a∈Ioc 0 a0, …≤C`，故 `hA` 槽不动 | **`[T]`** | §1.1；`HZE-KIN:16,64`（L3 只吃 `ϖ≤C`、`a`-一致） | `C` 要具体值 ⇒ 反解 `m₀(C)=2arccoth((C/4)^{1/4})`（`IAE-M:0.1`）；须先有非 log 型率（`HZU-M:60`） |
| **0.3** | 柱 2 的 `axiom` 是 `R5` 本体（`∃β₀∀β≥β₀∀L`）：不许换 `liminf`（见证 `w₀−1/k`）、不许塌成 `∃L₀`；`∀L ≡ inf_L` 可换 | **`[T]`** | `HZV-R5:0.1,0.2,0.4` | 量词面封；下一手只在 `β`-模数（`HZS-FINAL2:83` 三选一） |
| **0.4** | ★ **B 格的 `c` 不再自由**：`R5 ⇒ BFloor`，`c := w₀ = 1/(4N⁴)`；前件 = 夹逼左腿 `f₁²≤χs₁`（`Reduction:91–104`，`0 sorry`）+ 恒等 `ε=f₁²/M₀`（`HZP-EPS:0.1`） | **`[T]`** | §1.5；`HZS-FINAL2:59` 的 `R2` | `w₀` 与 `T9` 的 `(1−ε)c_L` 对齐仍未做（`HV5:95`） |
| **0.5** | **记账更正**：`HV5` 的 `hy6_N4_b_floor`（`s₁`-形、`c` 自由）**退役**（本条经 0.4 导出它；不删原文）；仓内 `axiom` 台账 = **3 条**：`N1`(C，`HV5:40`)·`N2`(A，本条)·`N4`(B，本条 R5 形)。`HZX:0.4` 的「2 条」按此登记 | **`[T]`**（记账） | §2；`HZX:14` 与 `HZX:75` 的 `N1/N4` vs `N1/N2` 两处指代不一致 | `hy6_clay_of_three` 的 `hA`/`hB` 补参（§2 表） |
| **0.6** | 柱 1 的供给路**不被本 `axiom` 选择**：`(V)@2.5 = [NG]`（`u=1.01439>1`）、`ρ_P`（一致 Poincaré = 谱隙，是墙）；`ϖ≤4coth(m/2)^4`、`ϖ≤144c₁²/ρ_P` 是两条**充分**路 | **`[T]`**（记账）/ **`[NG]`**（两条供给） | `HZU-M:16,60`；`HZG-BLOCK:13`；`HZF-COV:72` | 只登记出口；`C` 的存在与 `m>0` 是否可给出**解耦**（`HZU-M:67`） |
| **0.7** | **Clay**：两条 `axiom` = 两处缺件占位，装配仍是蕴含 | **`[NG]`** | §2；`HY6:52` | 见 §4 |

**一句话。** 两条声明都不缺类型：柱 1 缺的**一个数**（`C>4`，`(4,∞)` 开、`C=4` 逐字给 `m≥∞`），柱 2 缺的**一层量词**（`β`-一致模数；`g` 是下确界，点读只给上界）；`c` 被 `R5` 钉成 `w₀`，所以 `N4` 的 `s₁`-形是导出而非第二条假设。

## §1 Lean 声明（`axiom` 两条 + 桥；签名稿，`theorem` 不带 body）

```lean
-- KK1-LEAN 声明稿（signature-only，同 `HY6:§3.2–3.3` 规；**未编译**）
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

set_option autoImplicit false
namespace MassGap.KK1

/-! ## 1. 柱 1（A 格 = `HY6` 的 `N2`）：存在 `C>4` 使 `ϖ≤C` -/
variable (varpi : ℝ → ℝ) (a0 : ℝ)

/-- 靶。`ϖ = Var(q) + O_B`（`HZE-KIN:25`）；`C` 只需 `>4`（`C=4` 逐字给 `m≥∞`，`IAE-M:0.2`）；
    量词逐字同 `HY6:88` 的 `l3_scalar_pair`（`a ∈ (0,a0]`，`a`-一致）。 -/
def L3Bound : Prop := ∃ C : ℝ, 4 < C ∧ ∀ a : ℝ, a ∈ Set.Ioc 0 a0 → varpi a ≤ C
/-- 消费弱形：正是 `HY6:88` 的 `∃C ∀a∈Ioc 0 a0, ρ^ren_a 𝒢^ren_a ≤ C`（`varpi := ρ^ren𝒢^ren`）。 -/
theorem L3Bound.weaken (h : L3Bound varpi a0) :
    ∃ C : ℝ, ∀ a : ℝ, a ∈ Set.Ioc 0 a0 → varpi a ≤ C
variable (Varq OB : ℝ → ℝ)
/-- 靶的裂形（"缺的就一个数"）：`ϖ = Var(q) + O_B` ⇒ `O_B ≤ C − Var(q)`。 -/
theorem L3Bound_iff_split (hsplit : ∀ a : ℝ, varpi a = Varq a + OB a) :
    L3Bound varpi a0 ↔ ∃ C : ℝ, 4 < C ∧ ∀ a : ℝ, a ∈ Set.Ioc 0 a0 → OB a ≤ C - Varq a

/-- ★ **柱 1 的缺数**：`C>4`。`ha0` 防空窗（`a0≤0 ⇒ Ioc 0 a0 = ∅` ⇒ 界 vacuous）。内容 `[NG]`，类型 `[T]`。 -/
axiom kk1_N2_varpi_bound (ha0 : 0 < a0) : L3Bound varpi a0
/-! ## 2. 柱 2（B 格 = `HY6` 的 `N4`）：`R5` 与 `BFloor` -/
/-- 靶常数 `w₀ = r₀/(N²/4) = 1/(4N⁴)`（`HZP-EPS:41`；`N=2` 时 `1/64`）。 -/
noncomputable def w0 (N : ℕ) : ℝ := 1 / (4 * (N : ℝ) ^ 4)

/-- `R5`（`HZS-FINAL2:63` 逐字）：`ε := f₁²/Var(O)`（`HZP-EPS:0.1`）；`∃a₀∀a≥a₀` 在线读法下饱和删除
    （`HZV-R5:0.3`），`∀L` 与 `inf_L` 等价（`HZV-R5:0.1`）。 -/
def R5 (eps : ℝ → ℕ → ℝ) (N : ℕ) : Prop :=
  ∃ β0 : ℝ, ∀ β : ℝ, β0 ≤ β → ∀ L : ℕ, 2 ≤ L → w0 N ≤ eps β L

/-- `B` 格的消费类型（逐字 `HV5:51–55`；`βA`/`c` 装进 `∃` ⇒ `hB` 是无参槽）。 -/
def BFloor (s1 chi M0 : ℕ → ℝ → ℝ) : Prop :=
  ∃ βA : ℝ, ∃ c : ℝ, 0 < c ∧ ∀ L : ℕ, 2 ≤ L → ∀ β : ℝ, βA ≤ β →
    s1 L β ≥ c * (M0 L β / chi L β)

/-- ★ **桥**：`R5 ⇒ BFloor`，`c := w₀`。前件 = 夹逼左腿 + 等同；body 4 行
    （`le_div_iff₀` → `div_le_div_of_nonneg_right` → `le_trans`）。 -/
theorem kk1_R5_implies_BFloor
    (f1sq : ℕ → ℝ → ℝ) (eps : ℝ → ℕ → ℝ) (s1 chi M0 : ℕ → ℝ → ℝ) (N : ℕ)
    (hN : 2 ≤ N) (hχ : ∀ (L : ℕ) (β : ℝ), 0 < chi L β)
    (hM : ∀ (L : ℕ) (β : ℝ), 0 < M0 L β)
    (hLow : ∀ (L : ℕ) (β : ℝ), f1sq L β ≤ chi L β * s1 L β)
    (hEps : ∀ (L : ℕ) (β : ℝ), eps β L = f1sq L β / M0 L β)
    (hR5 : R5 eps N) : BFloor s1 chi M0

/-- ★ **柱 2 的缺量词**：`β`-一致模数（`g` 是下确界，点读只给上界，`HZV-R5:41`）。内容 `[NG]`，类型 `[T]`。 -/
axiom kk1_N4_r5_floor (eps : ℝ → ℕ → ℝ) (N : ℕ) (hN : 2 ≤ N) : R5 eps N

/-! ## 3. 消费面（`HY6:49–51,99–102`；本文件不改它们，只对位）
    T14 `hy6_gap_positive  (hres) (hB : ∀ a, Delta a ≥ δB) : DeltaCont ≥ δB`
    T15 `hy6_unique_limit   (hC) (hA : L3Bound varpi a0) (hB : BFloor s1 chi M0)`
    T16 `hy6_clay_of_three  (hC) (hA : L3Bound varpi a0) (hB : BFloor s1 chi M0) : ClayStatement` -/
end MassGap.KK1
```

## §2 消费形式匹配

| 消费点 | 签名（逐字） | 本文件供给（链） | 判 | 下一手 |
|---|---|---|---|---|
| **T14** `hy6_gap_positive` | `(hres : StrongResolvent Ha H) (hB : ∀ a, Delta a ≥ δB) : DeltaCont ≥ δB` | 不直给：`kk1_N4_r5_floor` → `kk1_R5_implies_BFloor`（`c=w₀`）→ `T3`（`s₁≥cM₀/χ ⟺ factor≥c`，`factor:=s₁χ/M₀`，`HY6:32`）→ `T13`（`(B)`，`HY6:42`）→ `hB` | **`[T]`**（四环全 `[T]`） | `T13` 的 `δB` 与 `w₀` 的显式关系未写（`HY6:42` 只到 `(B)`）⇒ 把 `δB` 写成 `δB(c)` 并代 `c=w₀` |
| **T15** `hy6_unique_limit` | `(hC : IotaL2Seam) (hA : L3Bound) (hB : BFloor)`（`HY6:50`；`§1` 记「C∧A∧B ⇒ 唯一极限」） | `hA` ← `kk1_N2_varpi_bound`；`hB` ← 上链；`hC` 仍由 `HV5:40` 的 `hy6_N1_iota_l2_seam` | **`[T]`**（槽位）/ **`[NG]`**（内容） | `hA`/`hB` 补参（`L3Bound varpi a0`/`BFloor s1 chi M0`）：`HY6` 的无参写法是缩写，不是更强 |
| **T16** `hy6_clay_of_three` | `(hC) (hA) (hB) : ClayStatement` | 同上 | 同上 | Clay 仍 `[NG]` |

## §3 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际调用 / 检索式（逐字） | 命中 |
|---|---|---|---|---|
| 1 | `pwsh` | **OK**（exit 0，正文非空） | `Get-ChildItem research\v3\species\ -Filter "HZ*" -Name`；`Select-String -Path research\v3\species\*.md -Pattern "T14\|T15" -List`；`Select-String … -Pattern "T14"`；`foreach($f in …){ (Get-Content $f \| Measure-Object -Line).Lines }`；`Get-ChildItem -Recurse -Filter *.lean`；`Get-ChildItem research\v3\lean -File`；`Select-String research\v3\lean\*.lean -Pattern "^\s*(axiom\|theorem\|lemma\|def\|structure\|noncomputable def)"` | 13 文件 / 40 行 / 17 行数 / 19962 `.lean`+60 路径 / 6 文件 ≈70 声明 |
| 2 | `read` | **OK**（带行号正文，全返回） | `HY6-ASSEMBLE`(130)·`HV5-LEAN`(97)·`HZV-R5`(88)·`IAE-M`(71)·`HZU-M`(88)·`HZG-BLOCK`(91)·`HZE-KIN`(84)·`HZP-EPS`(77)·`HZS-FINAL2`(103)·`HZX-COMPLETE`(107)·`Reduction.lean:88–132`·`L3B-Assembly.lean`(133) | **12/12** |
| 3 | `grep` | **OK**（有效响应；非 provider failure / 非 ENOENT / 非空响应） | ①`hy6_(N1\|N4\|gap_positive\|unique_limit\|clay_of_three\|L3\|b_floor)` @ `*.lean` ②`L3Bound\|BFloor\|VarpiBound\|hy6_L3\|hy6_A_\|kk1_` @ `*.md` ③`144` @ `species/*.md` | ①**0**（`No matches found`）②**2**（`HY6:101`、`HV5:67`）③**505**（显示 250；承重 `HZF-COV:18,66,72`） |
| 4 | 失败标记 | **TOOL-FAIL 0；空响应 0；PATH-MISS 0；provider failure 0** | ③ 命中溢出（250/505）已按 spill 路径读全，不记为失败 | — |

**直读 vs 转引（`AGENTS §3`）**：直读 = #2 十二份。只经命中行采信、未开正文：`HZF-COV:18,66,72`·`HQ4:19`（`w₀=r₀/(N²/4)`，经 `HZV-R5:5`/`HZP-EPS:41`）·`HS-GAP-8`·`HV9`·`GY1`·`HZ4-S1FLOOR` ⇒ 记转引级。

## §4 边界与下一手

1. **只交声明**：两条 `axiom` + 四个定义（`L3Bound`/`w0`/`R5`/`BFloor`）+ 三条桥的签名；无一 body 过 elaboration，不据此宣称新机检。
2. **窗口**（`AGENTS §2.1`）：`L3Bound` 的 `∀` 只压 `a`（`HZE-KIN:13,16` 的读法：`β/L/N/b̂` 在 `ϖ≤4S(m)` 里被消）；若要把 `∀β∀L` 也压进类型，等式更强，须先判它是否仍等于靶 —— 本文件不裁。
3. **读法**（`AGENTS §2.1`）：`R5` 按线读法（`a₀=a(β)` 饱和，`HZV-R5:0.3`）；乘积读法要 `∃a₀` 在 `∀β` 之外，是更强的 `R5^⊗`，不进本 `axiom`。
4. **两条都是内容 `[NG]`**：柱 1 缺 `C`（一个数，`>4` 即可）；柱 2 缺 `β`-一致模数（结构，不是数）。
5. 未改既有文件、未 `commit`、零 MC、无百分数、无 QED。

*KK1-LEAN · 两条 `axiom`（`N2` = `L3Bound`、`N4` = `R5`）+ 三桥签名；`c := w₀` 已钉 ⇒ `N4` 的 `s₁`-形降为导出。未编译；Clay `[NG]`。*

# Yang–Mills Lean 基础设施：包设计、精确类型与 mathlib 缺口

> 2026-11-22 · 范围：为 4D 纯 Yang–Mills 质量隙（MG）与其余并行线提供可复用的 Lean 4 基础设施。
> 交付物：`research/lean-infra/`（1,570 行 Lean，12 个模块，**0 sorry，0 自定义 axiom**，逐文件编译通过）。
> 2026-11-22 增补：并入 nontriviality 线的 Lemma C / RGCharge / 反例测试（§3.11），以及 group-uniformity 线的群参数化交接引理与无隐藏依赖证书（§3.12）。
> 本文不重复 `research/physics_handoff.md` 的交接控制量分析；它只负责把那些数学对象**变成类型**。

---

## 0. 结论摘要

1. **四个子系统各有明确的 Lean 载体**：格点规范场/轨道测度（`Gauge/`）、OS 重构接口（`Reconstruction/`）、无界自伴算子与二次型谱隙（`Spectral/`）、RG 步进与交接（`RG/`）。全部模块以 mathlib-only 为核心，OSforGFF 仅通过适配层接触（未在本包内导入，理由见 §4.1）。
2. **可证的部分已经证了，不可证的部分变成了结构字段。** 本包编译通过且无 `sorry`：凡是依赖 mathlib 尚未提供的定理（无界谱定理、紧群商测度分解、Wightman 重构）的地方，都写成**命名的假设束**（`SpectralTheoremInput`、`GaugeFixing.integral_eq`、`OSPositiveSemidefinite.cauchy_schwarz`），使下游定理可以"条件于 X"地陈述，并让缺口在类型层面可见——这正是 hexagonmath 的 16 输入模式。
3. **已经真正证出的七件事**（不是接口）：
   - 格点 Wilson 作用的规范不变性（任意群、任意类函数，`Gauge/Lattice.lean:wilsonAction_invariant`）；
   - 规范不变可观测量的期望在"构型空间"与"轨道空间"上相等（`Gauge/Orbits.lean:integral_orbitPushforward`）;
   - 弱→强耦合交接引理：由两侧界 + 重叠条件推出全尺度一致 gap（`RG/Handoff.lean:gap_uniform`、`regions_cover`、`overlap_iff_order`）；
   - 带速率的指数聚类 ⇒ OS4 裸聚类（`Statement/MassGap.lean:exponentialClustering_implies_clustering`），以及转移矩阵衰减 ⇒ 矩阵元指数衰减（`Spectral/Transfer.lean`）；
   - **Lemma C 的充分方向**：Gaussian 理论的高阶累积量全为零（`Statement/NonTriviality.lean:isGaussianTheory_cumulants_vanish`）；
   - **反例测试**：具体相互作用玩具理论非 Gaussian（`toyTheory_not_gaussian`）——证明非平凡性谓词不是空洞的；
   - **群一致性无隐藏依赖证书**：四个常数一致 ⇒ gap 界在群族上一致，且证明不触碰群（`RG/GroupUniform.lean:uniform_over_groups`）。
4. **最大的 mathlib 缺口有三个**：无界自伴算子的谱定理（G-4/G-5/G-6/G-10）、紧（不可数）群作用的商测度/轨道分解（G-2）、以及 OS→Wightman 重构所需的复分析与不定号正交群（G-8/G-9）。前两个直接卡住 MG，第三个只卡住 Clay 陈述的最终形式、**不卡住 MG 本身**（§4.2）。
5. **给并行线的建议**：MG 的陈述层应当以"转移矩阵在真空补上的指数衰减率"为准（`TransferData.HasGap`），而不是谱测度；`Statement/Basic.lean` 里现有的 `ExponentialClustering` 与本包的 `ExponentialClustering` 是同一命题的两个写法，应合并到一处。

---

## 1. 交付物与验证状态

```
research/lean-infra/
├── YangMillsInfra.lean              -- umbrella（全部模块）
├── lakefile.toml                    -- 独立包设计（mathlib v4.33.0-rc1）
├── check.sh                         -- 逐文件编译脚本（用 vendored olean 缓存）
├── REPORT.md / STATUS.md
└── YangMillsInfra/
    ├── Basic.lean                    88 行
    ├── Gauge/Lattice.lean           145 行
    ├── Gauge/Measure.lean           119 行
    ├── Gauge/Orbits.lean            111 行
    ├── Spectral/Unbounded.lean      120 行
    ├── Spectral/Transfer.lean       107 行
    ├── Reconstruction/OS.lean       140 行
    ├── RG/Handoff.lean              160 行
    ├── RG/GroupUniform.lean         156 行
    ├── RG/Step.lean                 130 行
    ├── Statement/MassGap.lean        86 行
    └── Statement/NonTriviality.lean 208 行
```

**验证方式.** 本 session 不做全量 `lake build`（任务要求）。环境里 `vendor/OSforGFF/.lake/packages/mathlib/.lake/build.badcache/` 存有一份完整的 mathlib olean 缓存（8,279 个 `.olean`，与本机 `leanprover/lean4:v4.33.0-rc1` 兼容），`check.sh` 用它设置 `LEAN_PATH` 后直接调用 `lean` 编译每个文件及其依赖，全部 **exit 0**。每个文件顶部只 import mathlib（外加本包自身模块），不 import OSforGFF、BochnerMinlos、GaussianField——因此本包可以脱离 Yang–Mills 主仓库独立构建与复用。

**没有做的事.** 没有证明任何 mathlib 缺口（那正是缺口）；没有全量构建；没有把 OSforGFF 的 `FieldConfiguration`、`SatisfiesAllOS` 直接接进来（见 §4.1 的适配层设计）。

---

## 2. 设计原则

### 原则 1：陈述层与证明层分离

质量隙的形式化最大的风险是**空洞陈述**（DeepMind `formal-conjectures` PR #7636 的教训：OS 字段写成 `True`）。本包的对策是：每个物理概念只定义一个类型或谓词，且该谓词必须**可被反例检验**。

例子：`Spectral.SelfAdjointOperator.HasGapAbove A V Δ` 定义为

```lean
0 < Δ ∧ ∀ ψ ∈ A.domain, ψ ≠ 0 → (∀ v ∈ V, ⟪ψ, v⟫ = 0) →
  Δ * ‖ψ‖^2 ≤ re ⟪ψ, A ψ⟫
```

它有具体的反例语义：对自由有质量场（`A = √(-Δ+m²)`）取 `Δ = m` 成立；对无质量场对任何 `Δ > 0` 都不成立（Rayleigh 商可任意小）。这不是 `True`。

### 原则 2：未证步骤 = 命名结构字段，而不是 `sorry`

`sorry` 会污染 axiom footprint，也让"哪些是假设"变得不可审计。本包把每个未证步骤写成结构字段，且字段名就是它对应的数学定理名：

| 结构 | 字段 | 对应的未证数学 |
|---|---|---|
| `GaugeFixing` | `integral_eq` | 紧群作用的商测度分解（G-2） |
| `WilsonGaugeInvariance` | `haar_invariant` | 紧群 Haar 的双不变性（G-3） |
| `OSPositiveSemidefinite` | `cauchy_schwarz` | OS 形式的极化/Cauchy–Schwarz（可证，但需 ~50 行代数） |
| `CertificateImpliesGap` | `transfer_gap` | 间隙证书 ⇒ 转移矩阵隙 |
| `SpectralTheoremInput`（设计） | `pvm`, `resolves` | 无界谱定理（G-4/G-5） |

于是任何下游定理都可以写成 `theorem foo (h : GaugeFixing μ) : ...`，读者立刻知道它条件于什么。

### 原则 3：核心 mathlib-only，适配层单独放

`YangMillsInfra/Basic.lean` 自己定义了 `SpaceTime`、`TestFunction`、`Configuration`、`EuclideanTheory`，与 `OSforGFF.Spacetime.Basic` 的对应物**名字相同、定义相同**，但不 import。这样：

- 本包可独立编译（已验证）；
- 其它线（os2-rotation、os3-reflection-positivity、nontriviality）可以直接用本包的类型陈述定理，不必拖入 31k 行 OSforGFF；
- 将来写 `YangMillsInfra/Bridge/OSforGFF.lean`，只写 `abbrev` 与实例，把 `OSforGFF.OS.Axioms.SatisfiesAllOS` 接到本包的 `EuclideanTheory` 上。

代价是两边定义需要保持同步；缓解办法是适配层里用 `rfl` 级别的 `example` 断言二者 defeq，任何一侧改动都会立刻编译失败。

---

## 3. 包结构与精确类型

### 3.1 `Basic.lean` — 共享词汇

```lean
abbrev SpaceTime (d : ℕ) : Type := EuclideanSpace ℝ (Fin d)
abbrev TestFunction (d : ℕ) : Type := SchwartzMap (SpaceTime d) ℝ
abbrev TestFunctionℂ (d : ℕ) : Type := SchwartzMap (SpaceTime d) ℂ
abbrev Configuration (d : ℕ) : Type := WeakDual ℝ (TestFunction d)

structure EuclideanTheory (d : ℕ) where
  Z : TestFunction d → ℂ            -- Glimm–Jaffe 生成泛函

structure TranslationAction (d : ℕ) where
  translate : SpaceTime d → TestFunction d → TestFunction d
  translate_zero : ∀ f, translate 0 f = f
  translate_add : ∀ a b f, translate (a + b) f = translate a (translate b f)

def ExponentialClustering (T) (τ) (m : ℝ) : Prop :=
  0 < m ∧ ∀ f g, ∃ C ≥ 0, ∀ a, ‖Z(f + τ_a g) − Z f · Z g‖ ≤ C · exp(−m‖a‖)

def Clustering (T) (τ) : Prop :=            -- OS4 原样
  ∀ f g ε > 0, ∃ R > 0, ∀ a, R < ‖a‖ → ‖Z(f + τ_a g) − Z f · Z g‖ < ε
```

**设计要点**：`EuclideanTheory` 只带 `Z`，不带测度。OS 公理是 `Z` 的性质，不是类型的一部分。这样同一个类型可以表示"只知 OS0–OS3"的理论，也不必在陈述层决定用 `𝒮′` 还是 `𝒟′`。`TranslationAction` 把 `SchwartzMap.translate`（mathlib 没有，OSforGFF 有）抽象成参数。

### 3.2 `Gauge/Lattice.lean` — 有限格点，精确层

组合数据：

```lean
structure LatticeData where
  V L Plaq : Type
  [Fintype V] [DecidableEq V] [Fintype L] [DecidableEq L] [Fintype Plaq]
  src tgt : L → V
  rev : L → L
  rev_src : ∀ l, src (rev l) = tgt l
  rev_rev : ∀ l, rev (rev l) = l
  p1 p2 p3 p4 : Plaq → L
  h12 : ∀ q, src (p2 q) = tgt (p1 q)
  h23 : ∀ q, src (p3 q) = tgt (p2 q)
  h34 : ∀ q, src (p4 q) = tgt (p3 q)
  h41 : ∀ q, src (p1 q) = tgt (p4 q)
```

(plaquette 存成四条有序边 + 四个边界等式，而不是 `Fin 4` 环；这样 `holonomy` 的四重积与规范变换的伸缩积都是纯 `group` 代数，见 §4.1。)

规范作用：

```lean
abbrev Connection (G) : Type := Λ.L → G
abbrev GaugeGroup (G) : Type := Λ.V → G
def gaugeAct (g) (U) (l) : G := g (src l) * U l * (g (tgt l))⁻¹
instance : MulAction (Λ.V → G) (Λ.L → G)

def holonomy (U) (q) : G := U (p1 q) * U (p2 q) * U (p3 q) * U (p4 q)

theorem holonomy_gaugeAct (g U q) :
    holonomy (g • U) q = g (src (p1 q)) * holonomy U q * (g (src (p1 q)))⁻¹

structure CharacterData (G) [Group G] where
  χ : G → ℝ
  χ_one : χ 1 = 1
  χ_conj : ∀ g h, χ (h * g * h⁻¹) = χ g

def wilsonAction (χ) (β) (U) : ℝ := β * ∑ q, (1 - χ (holonomy U q))

theorem wilsonAction_invariant (β g U) : wilsonAction χ β (g • U) = wilsonAction χ β U
```

**为什么不用 `SU(N)` 矩阵**：`CharacterData` 只需要"类函数 + 归一化"，`SU(N)` 的 `Re tr / N` 是一个实例。这样群论部分对所有紧群统一，也避免 mathlib 里 `Matrix.SpecialUnitaryGroup` 的完整 Lie 群结构（mathlib 目前没有 `SU(N)` 作为 Lie 群，只有 `Matrix.UnitaryGroup` 与 `SpecialLinearGroup` 的可定义性）。

### 3.3 `Gauge/Measure.lean` — Wilson 测度

```lean
class CompactGaugeGroup (G) extends Group G, TopologicalSpace G, MeasurableSpace G,
    BorelSpace G, CompactSpace G, T2Space G, IsTopologicalGroup G, SecondCountableTopology G

def haarProb (G) : Measure G := haarMeasure (univPositiveCompacts G)   -- 总质量 1
def linkHaar (Λ G) : Measure (Λ.L → G) := Measure.pi (fun _ => haarProb G)

def wilsonWeight (χ) (β) (U) : ℝ≥0∞ := ofReal (exp (-wilsonAction χ β U))
def partitionFunction (χ) (β) : ℝ≥0∞ := ∫⁻ U, wilsonWeight χ β U ∂(linkHaar Λ G)
def wilsonMeasure (χ) (β) (Z) := (linkHaar Λ G).withDensity (fun U => wilsonWeight χ β U / Z)

theorem wilsonMeasure_univ (hZ : Z = partitionFunction χ β) (h0 : Z ≠ 0) (htop : Z ≠ ∞) :
    wilsonMeasure χ β Z univ = 1

structure WilsonGaugeInvariance (χ) (β) (Z) : Prop where
  haar_invariant : ∀ g, Measure.map (fun U => g • U) (linkHaar Λ G) = linkHaar Λ G
  weight_invariant : ∀ g U, wilsonWeight χ β (g • U) = wilsonWeight χ β U

theorem wilsonWeight_invariant ...      -- 已证（来自 wilsonAction_invariant）
```

**关键点**：`wilsonMeasure` 是**定义**，不是假设。`wilsonMeasure_univ` 证明它是概率测度，只条件于 `Z` 是配分函数且有限非零——这两个是紧性/连续性的例行事实，留给后续。测度的规范不变性只剩 `haar_invariant` 一个字段（G-3）。

### 3.4 `Gauge/Orbits.lean` — 轨道测度

```lean
abbrev gaugeOrbitRel : Setoid (Λ.L → G) := MulAction.orbitRel (Λ.V → G) (Λ.L → G)
abbrev OrbitSpace := Quotient gaugeOrbitRel

def orbitLift (O) (hO : ∀ g U, O (g • U) = O U) : OrbitSpace → E :=
  Quotient.lift O (…)

def orbitPushforward (μ : Measure (Λ.L → G)) : Measure (OrbitSpace) := μ.map Quotient.mk''

theorem integral_orbitPushforward (μ) (O) (hO) (hπ) (hOm) :
    ∫ ω : OrbitSpace, orbitLift O hO ω ∂(orbitPushforward μ) = ∫ U, O U ∂μ

structure GaugeFixing (μ : Measure (Λ.L → G)) where
  slice : Set (Λ.L → G)
  isFundamentalDomain : IsFundamentalDomain (Λ.V → G) slice μ
  integral_eq : ∀ O (hO), ∫ U, O U ∂μ = ∫ ω, orbitLift O hO ω ∂(orbitPushforward μ)

def treeSlice (T : Finset Λ.L) : Set (Λ.L → G) := {U | ∀ l ∈ T, U l = 1}
```

**两个层次要分清**（这是本文件最重要的设计判断）：

- **推前（pushforward）** 是初等的、已证的：规范不变的观测量的期望与"在轨道空间上积分"一致。格点计算实际只用这个。
- **分解（disintegration）** 才是难的：把上层测度写成轨道上的测度积分，等价于"基本域 + 商测度定理"。mathlib 只在**可数**群作用情形有 `QuotientMeasureEqMeasurePreimage`；格点规范群 `Λ.V → G` 是紧 Lie 群、不可数，所以这一条是本包的 `GaugeFixing.integral_eq` 字段（缺口 G-2）。
- **树规范固定**是应该有具体实例：格点上沿极大树把每个构型唯一地规范到"树上全为 1"。`treeSlice` 给出切片定义；补上存在唯一性证明即可填 `GaugeFixing`，这是本包最值得优先攻的具体目标（§8）。

### 3.5 `Spectral/Unbounded.lean` — 无界自伴算子

```lean
structure SelfAdjointOperator (H) [NormedAddCommGroup H] [InnerProductSpace ℂ H] [CompleteSpace H] where
  op : H →ₗ.[ℂ] H
  selfAdjoint : IsSelfAdjoint op

def domain (A) : Submodule ℂ H := A.op.domain
theorem dense_domain (A) : Dense (A.domain : Set H)          -- 来自 mathlib
theorem isClosed (A) : A.op.IsClosed                          -- 来自 mathlib
def rayleigh (A) (ψ : A.domain) : ℝ := (⟪ψ, A ψ⟫).re / ‖ψ‖^2
def HasGapAbove (A) (V : Submodule ℂ H) (Δ : ℝ) : Prop :=
  0 < Δ ∧ ∀ ψ ∈ A.domain, ψ ≠ 0 → (∀ v ∈ V, ⟪ψ, v⟫ = 0) → Δ‖ψ‖² ≤ re ⟪ψ, A ψ⟫
def HasGapAboveVacuum (A) (Ω : H) (Δ) : Prop := A.HasGapAbove (span ℂ {Ω}) Δ
theorem inner_apply_im_zero (A) (ψ) : (⟪ψ, A ψ⟫).im = 0
theorem hasGapAbove_mono ...
```

**核心设计决策**：质量隙用**二次型下界**陈述，不用谱测度。理由：

1. 无界自伴算子的谱测度在 mathlib 里完全没有（G-4/G-5），而且构造它需要 Stieltjes/Herglotz 反演（G-10），是一个独立的大工程；
2. 二次型下界是"一阶"陈述，只需要 `LinearPMap`、内积与自伴性，全部现成；
3. 它与 `spec(A) ∩ (0,Δ) = ∅` 在谱定理下等价——这个等价正是 `SpectralTheoremInput` 要补的，补上之后所有下游定理不用改。

### 3.6 `Spectral/Transfer.lean` — 衰减 ⇒ 聚类

```lean
structure SemigroupGap (T : ℝ → H →L[ℂ] H) (V : Submodule ℂ H) (Δ : ℝ) : Prop where
  Δ_pos : 0 < Δ
  decay : ∀ t ≥ 0, ∀ ψ, (∀ v ∈ V, ⟪ψ, v⟫ = 0) → ‖T t ψ‖ ≤ exp(−Δ t) ‖ψ‖

theorem SemigroupGap.twoPoint_bound (h) (ht) (φ ψ) (hψ) :
    ‖⟪φ, T t ψ⟫‖ ≤ ‖φ‖ · exp(−Δ t) · ‖ψ‖

theorem SemigroupGap.decay_at_nat_multiples (h) (hmul : ∀ s t ≥ 0, T (s+t) = T s * T t)
    (hpres : T t 保持 Vᗮ) (ht) (n) (ψ) (hψ) :
    ‖T (n·t) ψ‖ ≤ exp(−Δ · n · t) ‖ψ‖
```

这是"容易的一半"：从算子的指数衰减到矩阵元与多步衰减。难点（从 Hamiltonian 下界推出 `decay`）在 §3.5 的接口之后。

### 3.7 `Reconstruction/OS.lean` — OS 重构接口

```lean
structure OSVector (d) where            -- ∑ c_i e^{i⟨φ,f_i⟩}，f_i 正时间支撑
  n : ℕ
  f : Fin n → TestFunctionℂ d
  c : Fin n → ℂ

def osForm (Z) (starF) (v w) : ℂ := ∑ i j, conj(v.c i) * w.c j * Z (v.f i - starF (w.f j))

structure OSPositiveSemidefinite (Z) (starF) : Prop where      -- OS3 + 推论
  nonneg : ∀ v, 0 ≤ (osForm Z starF v v).re
  cauchy_schwarz : ∀ v w, ‖osForm Z starF v w‖² ≤ re(osForm v v) * re(osForm w w)

def osNullSet (Z) (starF) (v) : Prop := (osForm Z starF v v).re = 0

structure TransferData where            -- 重构的**输出**
  H : Type
  [NormedAddCommGroup H] [InnerProductSpace ℂ H] [CompleteSpace H]
  vacuum : H
  vacuum_ne : vacuum ≠ 0
  transfer : ℝ → H →L[ℂ] H
  semigroup : ∀ s t ≥ 0, transfer (s+t) = transfer s * transfer t
  selfAdjoint : ∀ t ≥ 0, IsSelfAdjoint (transfer t)
  contraction : ∀ t ≥ 0, ‖transfer t‖ ≤ 1
  vacuum_invariant : ∀ t ≥ 0, transfer t vacuum = vacuum

def TransferData.HasGap (T) (Δ) : Prop :=
  Spectral.SemigroupGap T.transfer (span ℂ {T.vacuum}) Δ
```

**分层理由**：
- 第 0 层（`osForm`）：OS3 直接给出半正定，已由 OSforGFF 对 GFF 证明；
- 第 1 层（商 + 完备化）：数学上是"零空间取商"，需要 Cauchy–Schwarz 才能证明零集是子模；本包把它写成 `OSPositiveSemidefinite.cauchy_schwarz` 字段。补上它之后，用 `InnerProductSpace.Core` + `SeparationQuotient`/`Completion`（mathlib 都有实例，已验证 `UniformSpace.Completion.innerProductSpace`、`SeparationQuotient.instInnerProductSpace`）即可完成构造；
- 第 2 层（`TransferData`）：时间平移给出自伴压缩半群，真空不变。这一层是 MG 的直接载体，**不需要 Wightman 重构**。

### 3.8 `RG/Handoff.lean` — 交接引理（已证）

```lean
structure HandoffData where
  q g2 gap : ℝ → ℝ
  μup gstar q1 δW δS : ℝ
  q_strictAnti : StrictAntiOn q (Set.Ioi 0)
  g2_strictAnti : StrictAntiOn g2 (Set.Ioi 0)
  μup_pos : 0 < μup
  g2_μup : g2 μup = gstar^2
  gstar_pos : 0 < gstar

def WeakSide (D) : Prop := ∀ μ > 0, g2 μ ≤ gstar² → δW ≤ gap μ
def StrongSide (D) : Prop := ∀ μ > 0, q1 ≤ q μ → δS ≤ gap μ
def Overlap (D) : Prop := q1 ≤ q μup

theorem gap_uniform (hW : WeakSide D) (hS : StrongSide D) (hOv : Overlap D) :
    ∀ μ > 0, min δW δS ≤ gap μ

theorem regions_cover (hOv) : Set.Ioi 0 ⊆ weakRegion D ∪ strongRegion D
theorem overlap_iff_order (hqstar : q μstar = q1) : Overlap D ↔ μup ≤ μstar
```

这是 `physics_handoff.md` §3 的逐字形式化。三个定理都真正证了（`StrictAntiOn` + 分情况 + `linarith`，无 `sorry`）。`overlap_iff_order` 把"重叠条件"翻译成"`q` 参数下的序关系"，即交接判据 `q₀ ≥ q₁` 与"`μ_up ≤ μ*`"等价——这正是 handoff 报告强调的"不变量表述"。

### 3.9 `RG/Step.lean` — 步进定理与间隙证书

```lean
def gapCertificate (Z : ℕ → ℝ) (n : ℕ) : ℝ := 4 * Real.log (Z n) - Real.log (Z (2*n))

structure RGStepData where
  scale : ℕ → ℝ
  scale_pos : ∀ n, 0 < scale n
  scale_halves : ∀ n, scale (n+1) = scale n / 2
  error : ℕ → ℝ
  error_nonneg : ∀ n, 0 ≤ error n
  error_summable : Summable error

theorem RGStepData.accumulated_error_le (N) : ∑ n<N error n ≤ ∑' n, error n
theorem RGStepData.accumulated_error_tendsto : Tendsto (∑ n<N error n) atTop (nhds (∑' n, error n))

def VolumeUniform (P : ℕ → Prop) : Prop := ∀ L, P L

structure StepTheorem where
  blockScale : ℕ → ℝ
  blockScale_pos / blockScale_halves
  stepError : ℕ → ℕ → ℝ                     -- n 尺度, L 体积
  stepError_nonneg : ∀ n L, 0 ≤ stepError n L
  uniform_in_volume : ∀ n, ∃ ε ≥ 0, ∀ L, stepError n L ≤ ε
  summable_in_scale : ∀ L, Summable (fun n => stepError n L)

structure UniformGapCertificate (Z : ℕ → ℝ) where
  Z_pos : ∀ k, 1 ≤ Z k
  n n_pos lower lower_pos
  cert_bound : lower ≤ gapCertificate Z n

structure CertificateImpliesGap (T : Reconstruction.TransferData) (Z : ℕ → ℝ) where
  certificate : UniformGapCertificate Z
  transfer_gap : T.HasGap (certificate.lower / (4 * certificate.n))
```

`gapCertificate` 是 OpenAI sigma-model 模板里的 `Δ(n) = 4 log Z(n²) − log Z(2n²)`（以 `L`-盒配分函数记）。`StepTheorem` 的字段就是 Balaban/BBS 分析要输出的定量结论；`CertificateImpliesGap` 是唯一缺失的推理环节，显式列为字段。

### 3.10 `Statement/MassGap.lean`

```lean
structure MassGapStatement (T : EuclideanTheory d) (τ : TranslationAction d)
    (NonTrivial : EuclideanTheory d → Prop) : Prop where
  nonTrivial : NonTrivial T
  gap : ∃ m, ExponentialClustering T τ m

theorem exponentialClustering_implies_clustering {m} (h : ExponentialClustering T τ m) :
    Clustering T τ
```

非平凡性作为**参数**传入，因为它是 nontriviality 线的所有物，陈述层不应重复定义。这里的定理证明了"带速率聚类 ⇒ OS4"，即速率形式确实是 OS4 的加强而非无关条件。

### 3.11 `Statement/NonTriviality.lean` — Gaussianity、累积量与 RGCharge

并入 nontriviality 线的结论（`research/nontriviality/REPORT.md` §1.2）：**纯规范理论中，"某个规范不变连通 4 点 ≠ 0"是空洞见证**——`tr F²` 是 Gaussian 场的 Wick 平方，自由规范场已经通过。唯一同时满足"规范不变"与"自由值 = 0"的见证是重整化耦合 `g_ren`（等价地 `Λ > 0`）。

```lean
def IsGaussianTheory (T : EuclideanTheory d) : Prop :=
  ∃ B : TestFunction d →ₗ[ℝ] TestFunction d →ₗ[ℝ] ℝ,
    ∀ J, T.Z J = Complex.exp (-(1/2 : ℂ) * ((B J J : ℝ) : ℂ))

def NonTrivial (T) : Prop := ¬ IsGaussianTheory T

def cumulantGF (T) (J) : ℝ → ℝ := fun t => (Complex.log (T.Z (t • J))).re
def nthCumulant (T) (J) (n) : ℝ := iteratedDeriv n (cumulantGF T J) 0

theorem isGaussianTheory_cumulants_vanish (h : IsGaussianTheory T) :
    ∀ n ≥ 3, ∀ J, nthCumulant T J n = 0                      -- Lemma C 充分方向，已证

structure CumulantCharacterization (T) : Prop where           -- Lemma C 必要方向（接口）
  analytic_cgf : ∀ J, AnalyticAt ℝ (cumulantGF T J) 0
  gaussian_of_vanishing : (∀ n ≥ 3, ∀ J, nthCumulant T J n = 0) → IsGaussianTheory T

theorem isGaussianTheory_iff_cumulants_vanish (h : CumulantCharacterization T) :
    IsGaussianTheory T ↔ ∀ n ≥ 3, ∀ J, nthCumulant T J n = 0

def toyTheory (L : TestFunction d →ₗ[ℝ] ℝ) (lam : ℝ) : EuclideanTheory d where
  Z J := ((Real.exp (-(1/2) * (L J)^2 - lam * (L J)^4) : ℝ) : ℂ)

theorem toyTheory_not_gaussian (hL : L ≠ 0) (hlam : lam ≠ 0) :
    ¬ IsGaussianTheory (toyTheory L lam)                     -- 反例测试，已证

structure RGCharge (T) where
  value : ℝ
  free_vanishes : IsGaussianTheory T → value = 0

def GaugeNonTrivial (T) : Prop := ∃ c : RGCharge T, c.value ≠ 0
theorem gaugeNonTrivial_implies_not_gaussian : GaugeNonTrivial T → ¬ IsGaussianTheory T
```

**约定说明**：`cumulantGF` 取 `log Z(t • J)` 的实部。对反射正/对称理论 `Z(t•J)` 是实正数，这就是通常的累积量生成函数；取实部避免了 `log` 的分支问题，且不影响 Lemma C 的充分方向。必要方向作为 `gaussian_of_vanishing` 显式接口——它需要 OS0 型解析性，不能由"累积量全灭"单独推出。

**与 `YMStatement/Basic.lean` 的对应**：`IsGaussianTheory` 与 `IsGaussianMeasure` 在测度支撑的理论上相同（同样的 `∃ B` 形式，只是前者说 `Z`、后者说生成泛函）；`RGCharge.free_vanishes` 与新增的 `YMStatement.RGCharge` 字段一致。区别仅在于本包不依赖 OSforGFF。

### 3.12 `RG/GroupUniform.lean` — 群参数化交接与无隐藏依赖

并入 group-uniformity 线的 N1 建议（`research/group-uniformity/REPORT.md` §9.2）：把交接引理的四个群相关常数 `q₀(G), q₁(G), δ_W(G), δ_S(G)` 显式参数化，并证明结论只依赖其中三个。

```lean
structure GroupHandoffData (Grp : Type u) where
  q g2 gap : ℝ → ℝ
  q0 q1 δW δS : ℝ
  q_strictAnti : StrictAntiOn q (Set.Ioi 0)
  g2_strictAnti : StrictAntiOn g2 (Set.Ioi 0)
  μup gstar : ℝ
  μup_pos / gstar_pos / g2_μup / q0_eq : ...
  overlap : q1 ≤ q0

def GroupHandoffData.WeakSide (D) : Prop := ∀ μ > 0, g2 μ ≤ gstar² → δW ≤ gap μ
def GroupHandoffData.StrongSide (D) : Prop := ∀ μ > 0, q1 ≤ q μ → δS ≤ gap μ

theorem GroupHandoffData.gap_uniform (hW) (hS) : ∀ μ > 0, min δW δS ≤ gap μ

theorem uniform_over_groups {ι} {Grp : ι → Type u} (D : ∀ i, GroupHandoffData (Grp i))
    (hW : ∀ i, (D i).WeakSide) (hS : ∀ i, (D i).StrongSide)
    (hq1 : ∀ i, (D i).q1 = q1) (hWc : ∀ i, (D i).δW = δW) (hSc : ∀ i, (D i).δS = δS) :
    ∀ i μ > 0, min δW δS ≤ (D i).gap μ

structure Ledger where            -- 台账 (α, β, γ)：dim G、h∨、dim f 的幂次
  alpha beta gamma : ℕ
def Ledger.intensive / representation / extensive
```

**证书的含义**：`uniform_over_groups` 的证明从不提及 `Grp`——群只通过 `q1, δW, δS` 进入结论（`q0` 只进入 overlap 假设）。这就是"N1：把群一致性局部化为四个数的一致界"的可检查版本。`Ledger` 把 group-uniformity 报告 §6 的 `(α,β,γ)` 台账变成可携带在陈述上的注解。

---

## 4. 四个关键设计决策

### 4.1 为什么不（直接）导入 OSforGFF

OSforGFF 是当前最好的起点（25.5k 行、0 sorry、`#guard_msgs` 冻结 axiom footprint），但它：

- 依赖 `bochner`、`gaussian-field`、`kolmogorov_extension4` 三个外部包，构建链长；
- 它的 `FieldConfiguration`、`SatisfiesAllOS` 是**Gaussian 专用**的接口（`GJGeneratingFunctional` 定义在 `ProbabilityMeasure (WeakDual …)` 上）；
- 其余并行线（如 uniqueness、os3）需要的是陈述层类型，不是测度构造。

因此本包复制了 4 个 `abbrev`（`SpaceTime`、`TestFunction`、`TestFunctionℂ`、`Configuration`）与一个 `structure EuclideanTheory`，其余全部独立。适配层的正确写法是：

```lean
-- YangMillsInfra/Bridge/OSforGFF.lean（未交付，设计）
import OSforGFF.OS.Axioms
import YangMillsInfra.Basic
namespace YangMillsInfra
example : OSforGFF.SpaceTime d = SpaceTime d := rfl      -- 守住同步
example : OSforGFF.FieldConfiguration d = Configuration d := rfl
def EuclideanTheory.ofMeasure (μ : ProbabilityMeasure (Configuration d)) : EuclideanTheory d :=
  ⟨fun J => OSforGFF.GJGeneratingFunctional μ J⟩
theorem satisfiesAllOS_iff ... : OSforGFF.SatisfiesAllOS μ ↔ ... := …
end
```

好处：任何一侧的定义漂移都会在 `example` 处编译失败。

### 4.2 质量隙不需要 Wightman 重构

Clay 陈述要求"至少不弱于 Wightman/OS 公理的公理性质"，所以最终必须做 OS→Wightman。但 **MG 是纯 Euclidean 陈述**：体积一致、标号一致的指数聚类率。给定 `SatisfiesAllOS`（已是 `EuclideanTheory` 的性质），MG 的证明与验证不需要 Minkowski 侧。

因此本包把重构分成两层，并把 Wightman 层明确推迟：

- **Level A（本包）**：OS 形式 → 商 Hilbert 空间 → 转移矩阵 → `TransferData.HasGap`。这是 MG 的全部载体。
- **Level B（不做）**：`TransferData` → Wightman 理论。需要：多元复分析（管域定理、边值分布）、Laplace 变换、不定号正交群 `O(1,3)`。mathlib 全无（G-8/G-9）。

这个切分对项目有直接影响：即使 Level B 永远不做，Level A + 强耦合 + 交接仍然是可发表的形式化工件；反过来，先做 Level B 会把资源投在 MG 不需要的地方。

### 4.3 谱隙用二次型而不是谱测度

见 §3.5。补充一点：用二次型还有一个工程优势——它是**局部**陈述（对每个 `ψ` 验证一个不等式），而谱隙是**全局**陈述（关于集合 `spec(H)`）。前者在 Lean 里由 `linarith`/`nlinarith` 直接处理，后者需要拓扑与测度。

### 4.4 格点优先，连续统只留接口

连续统规范场（主丛、Ehresmann 联络、曲率）在 mathlib 与 PhysLean 都**不存在**（G-1）。4D 构造本身就完全在格点上做、最后取标度极限，所以：

- 本包的 `Gauge/` 全部是有限格点，所有对象有限维、可计算；
- 连续统只出现在 `Basic.Configuration`（`𝒮′`）与 `EuclideanTheory`，作为**极限的输出**，不作为中间对象；
- 这意味着"连续统联络空间上的测度"这个表述在本项目里应当**避免**——它不是必需的，而且无法形式化。

---

## 5. mathlib 缺口清单

按对 MG 的阻塞程度排序。证据列给出缺失声明所在的目录/文件，均为 `mathlib v4.33.0-rc1`（commit `79d0395`）。

| 编号 | 缺口 | 证据 | 对 MG | 建议 |
|---|---|---|---|---|
| **G-4** | 无界自伴算子的谱定理 | `Mathlib/Analysis/InnerProductSpace/LinearPMap.lean` 只有 `adjoint`、`IsSelfAdjoint`、`isClosed`；无 `SpectralTheorem`、无 `Mathlib/Analysis/Spectral` | 阻塞"二次型隙 ↔ 谱隙" | 先做二次型（本包）；谱定理作为独立长线 |
| **G-5** | 投影值测度（PVM）/谱族 / Stone 定理 | 全库 grep `ProjectionValued`、`SpectralMeasure` 无结果 | 阻塞"⇒ `e^{-tH}` 衰减" | 与 G-4 同一工程 |
| **G-6** | 无界算子的可测函数演算 | `Mathlib/Analysis/CStarAlgebra/ContinuousFunctionalCalculus/` 只对 C*-代数元素（有界） | 阻塞从 `HasGapAbove` 推 `SemigroupGap` | 先在有界/解析半群情形做特例 |
| **G-10** | Stieltjes/Herglotz 变换与谱测度构造 | 无 `Stieltjes` 变换、无 Herglotz 定理 | G-4 的技术核心 | 需要先有半平面上的函数论 |
| **G-2** | 紧（不可数）群作用的商测度/轨道分解 | `MeasureTheory/Group/FundamentalDomain.lean` 的 `QuotientMeasureEqMeasurePreimage` 及其非平凡引理要求 `[Countable G]`；`Measure/Haar/Quotient.lean` 只处理 `G ⧸ Γ`（子群） | 阻塞"轨道空间上的测度"严格化 | 优先做**树规范固定**：`IsFundamentalDomain.mk'` 的 `∃! g` 形式不需要可数性 |
| **G-3** | 紧群 Haar 的右不变性/幺模性 | `Haar/Basic.lean` 有 `isMulLeftInvariant_haarMeasure`；全库无 `isMulRightInvariant_haarMeasure`；`unimodular` 只出现在注释 | 阻塞 Wilson 测度的规范不变性（`haar_invariant`） | 约 50–100 行：用 `haarMeasure_unique` + 反转测度 |
| **G-1** | 连续统规范理论几何 | 无主丛、无 Ehresmann 联络；`Mathlib/Geometry/Manifold/VectorBundle/CovariantDerivative` 只是向量丛协变导数；PhysLean 也只有张量记法 | 不阻塞（格点足够） | 不做 |
| **G-8** | OS→Wightman 重构 | 无管域定理、无分布边值、无 Wightman 公理 | 不阻塞 MG，阻塞 Clay 最终陈述 | Level B，长期 |
| **G-9** | 不定号正交群 `O(1,d−1)`、Minkowski 空间 | 只有 `Matrix.orthogonalGroup`（正定）；无 `LorentzGroup` | 同上 | 同上 |
| **G-11** | 标准 Borel 商上的可测选择 | `Probability/Kernel/Disintegration/StandardBorel.lean` 有 `condKernel`，但无"紧群作用商是标准 Borel + 可测截面"定理 | 与 G-2 同源 | 可先对有限格点（紧度量空间）特化 |
| **G-12** | 乘积测度在"逐坐标双射"下的不变性 | `MeasureTheory/Constructions/Pi.lean:385 pi_map_pi` 已有；只差每坐标的 `MeasurePreserving` | 不阻塞 | 依赖 G-3 |
| **G-13** | Sobolev / 加权 Sobolev / 高斯测度在 Sobolev 上 | 无 `Sobolev` 目录（`Mathlib/Analysis/` 下无） | 阻塞 Balaban 型 RG 范数 | 长期；可先用离散 ℓ² 范数替代 |
| **G-14** | 多元 Gaussian、Hilbert 空间值 Gaussian | `GaussianField` 包提供了核化的 Gaussian 构造（见 §6.2），但 mathlib 本身无 | 不阻塞（已有外部包） | 用 `GaussianField` |
| **G-15** | 紧群上连续函数积分的正性/有限性 | 有 `IsHaarMeasure`、`Continuous`、`integral_pos` 相关引理，但拼装成"配分函数 ∈ (0,∞)"需要若干 glue | 不阻塞 | 小工程，可顺手做 |
| **G-16** | Casimir 元素与指标公式 | 全库 grep `Casimir` **零命中**（group-uniformity §8 O7 已核对）；有根系/Weyl 群/Killing 型/Cartan 矩阵 | 阻塞 2D 精确解与 Casimir scaling（L0 的一半、L6） | 需要对半单 Lie 代数的泛包络代数中心；中大型工程 |
| **G-17** | Weyl 积分公式 | grep `Weyl integration` 零命中；Haar 只有抽象测度 | 阻塞单 plaquette 系数（L2）、Weyl 特征标展开 | 可在有 G-16 之后做；需要极大环面与 Weyl 群作用 |
| **G-18** | Weyl 特征标/维数公式 | grep `Weyl character formula` 零命中 | 阻塞 L2 的收敛估计与 `dim R` 权重 | 与 G-16/G-17 同一工程 |
| **G-19** | 紧 Lie 群的 Peter–Weyl 定理 | grep `peter-weyl`/`peterweyl` **零命中** | 阻塞 2D 精确解（L6）与所有表示展开 | 长期；N2 建议先做 `U(1)`、`SU(2)` 的有限维截断绕开 |

**关于 G-2 的补充**：`IsFundamentalDomain.mk'` 的形式是

```lean
theorem IsFundamentalDomain.mk' (h_meas : NullMeasurableSet s μ)
    (h_exists : ∀ x : α, ∃! g : G, g • x ∈ s) : IsFundamentalDomain G s μ
```

它**不需要** `[Countable G]`。所以树规范固定的路线是：构造切片 `treeSlice T`，证明 `∀ U, ∃! g, g • U ∈ treeSlice T`（沿树归纳），再用 `mk'` 得到基本域。真正需要可数性的是 `QuotientMeasureEqMeasurePreimage` 的**商测度**部分；但如果我们只用推前 `orbitPushforward`（§3.4），就不需要商测度定理。**这是本包认为最值得优先完成的具体 Lean 目标。**

---

## 6. 对现有工件的调研

### 6.1 mathlib（v4.33.0-rc1）

可用且已验证（本包实际使用）：

- `LinearPMap`：`domain`、`toFun`、`adjoint`、`Star` 实例、`IsSelfAdjoint`、`IsSelfAdjoint.dense_domain`、`IsSelfAdjoint.isClosed`；
- `MeasureTheory.Measure.pi`（有限乘积）、`pi_map_pi`、`withDensity`、`lintegral_mul_const'`、`ENNReal.mul_inv_cancel`；
- `MeasureTheory.Measure.haarMeasure`（注意 namespace 是 `MeasureTheory.Measure`，且参数是 `TopologicalSpace.PositiveCompacts G`）、`haarMeasure_self`；
- `MeasureTheory.IsFundamentalDomain`、`QuotientMeasureEqMeasurePreimage`（可数群）；
- `UniformSpace.Completion.innerProductSpace`、`SeparationQuotient.instInnerProductSpace`、`InnerProductSpace.ofCore`、`InnerProductSpace.Core.mk`；
- `Mathlib.Probability.Kernel.Disintegration.StandardBorel`（`condKernel`）；
- `StrictAntiOn`、`intermediate_value_Ioo`、`Summable.sum_le_tsum`、`Summable.hasSum.tendsto_sum_nat`。

明确缺失（见 §5）。

### 6.2 OSforGFF（`vendor/OSforGFF`，Douglas–Hoback–Mei–Nissim）

结构（52 文件，31,429 行）：

```
OS/Axioms.lean          OS0–OS4 的谓词（Glimm–Jaffe 形式，𝒮′ 上）
OS/Master.lean          GFF 满足全部 OS 公理（d ≥ 2 通用 + d=2,3,4,5 实例）
OS/OS3_ReflectionPositivity.lean / OS3_MixedRep*.lean
OS/OS4_Clustering.lean  OS4 聚类（Gaussian 分解 + 卷积衰减）
OS/OS4_Ergodicity.lean  多项式聚类 α=6 ⇒ L² 遍历
Measure/GaussianFreeField.lean, Measure/Minlos.lean, Measure/NuclearSpace.lean
Spacetime/{Basic,Euclidean,PositiveTimeTestFunction,TimeTranslation,DiscreteSymmetry}.lean
Guardrails.lean         #guard_msgs 冻结 6 个主定理的陈述与 axiom footprint
```

可直接复用（设计上）：

- `SatisfiesAllOS` 的字段结构 → 本包 `EuclideanTheory` 的适配目标；
- `OS3_ReflectionPositivity` 的 `star` 约定（时间反射 + 复共轭）→ 本包 `osForm` 的 `starF` 参数；
- `PositiveTimeTestFunctionℂ` 的实现（`tsupport ⊆ positiveTimeSet`）→ 重构层的底空间；
- `QFT.E d`（`EuclideanSpace ⋊ O(d)`，有 `Group` 实例与作用引理）→ OS2 的表示论基础；
- `Guardrails` 模式 → 本包与其它线应复制这一做法（用 `#guard_msgs` 冻结陈述）。

**与 Douglas et al. arXiv:2603.15770 的关系**：论文标题 "Formalization of QFT"，摘要只声称 4D 自由玻色场 + Glimm–Jaffe/OS 变体公理，并明确说"不声称证明重构定理或谱理论"。所以本包 §4.2 的 Level A/Level B 切分与论文的范围一致：**公开工件里没有 OS→Wightman 重构**，prior-work 文档 §4.6 的判断（"公开搜索未找到独立确认"）成立。

### 6.3 依赖包（都在 `vendor/OSforGFF/.lake/packages/`）

- **BochnerMinlos**：`minlos_theorem`（已证，非 axiom）：对核 Fréchet 空间 `E` 与连续正定归一化 `Φ`，存在唯一概率测度在 `WeakDual ℝ E` 上实现 `Φ`。这是连续统 Gaussian 参考测度的唯一入口。
- **GaussianField**：`Construction.lean`（`charFun`）、`Density.lean`（`FinLatticeField d N` 上的高斯密度）、`SpectralTheorem.lean`（**紧自伴算子的谱定理**，实 Hilbert 空间，已证）、`Wick.lean`、`Hypercontractive.lean`、`Tightness.lean`。注意 `SpectralTheorem.lean` 是**紧**算子，不能用于 Hamiltonian。
- **KolmogorovExtension4**：`MeasureTheory.projectiveLimit`（对 `IsProjectiveMeasureFamily` 给出 `Π i, α i` 上的测度），是无限体积格点极限的现成工具。

### 6.4 PhysLean

社区物理库（Relativity ~27.5k 行、QFT ~20k、Particles ~19k、Electromagnetism ~9.4k）。关键事实（2026 年状态）：

- 规范理论用**张量记法**，不是主丛/联络；
- 主丛上的 Ehresmann 联络是"规划中"，Bianchi 恒等式 ⇒ Maxwell 尚未打通；
- 没有构造性 QFT 的测度论。

结论：PhysLean 目前**不能**为本项目提供连续统规范几何或 QFT 测度，只能提供指标记法与某些经典场论结果。与 prior-work 文档 §4.2 的判断一致。

---

## 7. 与并行线的接口

| 并行线 | 应使用本包的 | 需要本包补的 |
|---|---|---|
| os2-rotation | `Basic.EuclideanTheory`、`TranslationAction`（可扩展为 `EuclideanAction`） | 一个 `EuclideanAction` 结构（把 `QFT.E d` 的作用包成测试函数作用） |
| os3-reflection-positivity | `Reconstruction.osForm`、`OSPositiveSemidefinite` | `starF` 与 OSforGFF `star` 的桥接引理 |
| nontriviality | `Statement.NonTriviality`：`IsGaussianTheory`、`nthCumulant`、Lemma C、`RGCharge`、`GaugeNonTrivial`、`toyTheory_not_gaussian` | 与 `YMStatement` 的 `IsGaussianMeasure`/`RGCharge` 做 defeq 桥接；补 `CumulantCharacterization.gaussian_of_vanishing` |
| group-uniformity | `RG/GroupUniform`：`GroupHandoffData`、`uniform_over_groups`、`Ledger` | `SU(N)`/例外群的 `C₂`、`h∨`、`dim f` 实例（受阻于 G-16） |
| uniqueness | `Basic.Clustering`、`RG/Handoff` | 体积一致性的表述（`VolumeUniform`）可复用 |
| mass-spectrum | `Spectral/Transfer`、`Reconstruction.TransferData` | 两点函数的 Källén–Lehmann 表示（需要 PVM，G-5） |
| confinement | `Gauge/Lattice`（Wilson 环、holonomy） | `WilsonLoop` 与面积律谓词（可直接在 `LatticeData` 上加） |
| gauge-invariance | `Gauge/Lattice`、`Gauge/Orbits` | 树规范固定的存在唯一性（§5 G-2 建议） |
| group-uniformity | `Gauge/Measure` 的 `CharacterData` | `SU(N)` 的 `Re tr / N` 实例与 `Haar` 的 `SU(N)` 测度 |
| 主 session（交接） | `RG/Handoff`（已证） | 无；把 `HandoffData` 的实例提供出来即可 |

**命名冲突提醒**：`OSforGFF.OS.Axioms.OS4_Clustering` 与本包 `Basic.Clustering` 是同一命题；建议在本包加适配引理 `Clustering_iff_OS4`，并在整合阶段统一到一处。

---

## 8. 路线图与工作量估计

按"可验收工件"排序，估计为单人 + AI 辅助的日历时间（含调试）。

| 优先级 | 任务 | 产出 | 估计 |
|---|---|---|---|
| **P0** | 树规范固定：`∀ U, ∃! g, g • U ∈ treeSlice T` | 填 `GaugeFixing.isFundamentalDomain`（用 `IsFundamentalDomain.mk'`） | 1–2 周 |
| ✅ 已完成 | 群参数化交接 + 无隐藏依赖证书（group-uniformity N1） | `RG/GroupUniform.lean`（`uniform_over_groups`，已证） | — |
| ✅ 已完成 | 非平凡性：Lemma C 充分方向、反例测试、`RGCharge` | `Statement/NonTriviality.lean`（已证） | — |
| **P0** | 紧群 Haar 右不变性 | 填 `WilsonGaugeInvariance.haar_invariant`（G-3） | 1 周 |
| **P0** | OS 形式的 Cauchy–Schwarz 与零空间子模 | 填 `OSPositiveSemidefinite.cauchy_schwarz`，完成 `OSHilbert` 构造 | 2–3 周 |
| **P1** | 连续统/格点适配层 `Bridge/OSforGFF.lean` | `EuclideanTheory.ofMeasure` + `SatisfiesAllOS` 桥接 | 1 周 |
| **P1** | `SU(N)` 实例：`CharacterData`、Haar、`U_p` 的矩阵实现 | `Gauge/Instances/SU_N.lean` | 2–3 周 |
| **P1** | 间隙证书 ⇒ 转移矩阵隙（有限维/紧算子情形） | 填 `CertificateImpliesGap`（用 `GaussianField.SpectralTheorem`） | 3–4 周 |
| **P2** | 有界算子情形的函数演算 ⇒ `SemigroupGap` | `Spectral/CFC.lean`：`H` 有界时的完整推导 | 4–6 周 |
| **P2** | 两点函数的谱表示（Källén–Lehmann） | 需要 PVM；`Spectral/PVM.lean` | 2–4 月 |
| **P3** | 无界谱定理 / Stone 定理 | 上游 mathlib 贡献级 | 6–12 月 |
| **P3** | Wightman 重构 / 管域定理 | `Reconstruction/Wightman.lean` | 1 年+ |

**判断**：P0 三项是"低风险、高价值"，做完之后本包从"设计 + 接口"变成"可用的格点 + OS 基础设施"，足够让强耦合与 2D YM 线独立推进。P1 的间隙证书是 MG 论证链上最或缺的一环。P3 不应作为近期承诺。

---

## 9. 明确不形式化的内容

- **连续统主丛/联络/曲率**（G-1）：格点足够，且 mathlib 完全空白。
- **Wightman 公理与重构定理**（G-8/G-9）：MG 不需要；Clay 最终陈述才需要，属于 Level B。
- **Balaban 多尺度展开的具体范数**（G-13）：本包只固定 `StepTheorem` 的接口；范数空间的选择（大场/小场分解）是 RG 线的核心工作。
- **hexagonmath 的"有聚类无谱隙"反例模型**：形式化它需要构造一个具体的非 Gaussian OS 测度，代价高、收益是"证明速率必要性"。本包用 `ExponentialClustering` 的**定义**（带速率）已经排除了该反例路径；把反例形式化列为可选的后续。
- **数值/格点模拟**：不在本包范围（`numerics/` 独立）。

---

## 10. 复现

环境：`leanprover/lean4:v4.33.0-rc1`（`~/.elan` 已装），mathlib olean 缓存在
`vendor/OSforGFF/.lake/packages/mathlib/.lake/build.badcache/lib/lean`（8,279 个 `.olean`，与该 toolchain 兼容）。

```bash
cd research/lean-infra
./check.sh YangMillsInfra/Gauge/Lattice.lean \
           YangMillsInfra/Gauge/Measure.lean \
           YangMillsInfra/Gauge/Orbits.lean \
           YangMillsInfra/Spectral/Unbounded.lean \
           YangMillsInfra/Spectral/Transfer.lean \
           YangMillsInfra/Reconstruction/OS.lean \
           YangMillsInfra/RG/Handoff.lean \
           YangMillsInfra/RG/GroupUniform.lean \
           YangMillsInfra/RG/Step.lean \
           YangMillsInfra/Statement/MassGap.lean \
           YangMillsInfra/Statement/NonTriviality.lean \
           YangMillsInfra.lean
```

`check.sh` 先编译 `YangMillsInfra/Basic.lean`，再把每个目标文件编译为 `.olean`。全程约 2–4 分钟（主要是 Lean 启动与 import 解析），不需要 `lake build`。

**若将来做全量构建**：`lakefile.toml` 只 require mathlib，`lake build` 会从源码重建 mathlib（本机约数小时、数 GB），不建议在开发循环里做。

---

## 附：文件与声明索引

| 文件 | 主要声明 |
|---|---|
| `Basic.lean` | `SpaceTime`, `TestFunction`, `Configuration`, `EuclideanTheory`, `TranslationAction`, `ExponentialClustering`, `Clustering` |
| `Gauge/Lattice.lean` | `LatticeData`, `gaugeAct`, `MulAction`, `holonomy`, `holonomy_gaugeAct`, `CharacterData`, `wilsonAction`, `wilsonAction_invariant` |
| `Gauge/Measure.lean` | `CompactGaugeGroup`, `haarProb`, `linkHaar`, `wilsonWeight`, `partitionFunction`, `wilsonMeasure`, `wilsonMeasure_univ`, `WilsonGaugeInvariance`, `wilsonWeight_invariant` |
| `Gauge/Orbits.lean` | `gaugeOrbitRel`, `OrbitSpace`, `orbitLift`, `orbitPushforward`, `integral_orbitPushforward`, `GaugeFixing`, `treeSlice` |
| `Spectral/Unbounded.lean` | `SelfAdjointOperator`, `domain`, `dense_domain`, `isClosed`, `rayleigh`, `HasGapAbove`, `HasGapAboveVacuum`, `inner_apply_im_zero`, `hasGapAbove_mono` |
| `Spectral/Transfer.lean` | `SemigroupGap`, `twoPoint_bound`, `decay_at_nat_multiples` |
| `Reconstruction/OS.lean` | `OSVector`, `osForm`, `OSPositiveSemidefinite`, `osNullSet`, `TransferData`, `TransferData.HasGap` |
| `RG/Handoff.lean` | `HandoffData`, `WeakSide`, `StrongSide`, `Overlap`, `weakRegion`, `strongRegion`, `gap_uniform`, `regions_cover`, `overlap_iff_order` |
| `RG/GroupUniform.lean` | `GroupHandoffData`, `WeakSide`, `StrongSide`, `gap_uniform`, `uniform_over_groups`, `Ledger` |
| `RG/Step.lean` | `gapCertificate`, `RGStepData`, `accumulated_error_le`, `accumulated_error_tendsto`, `VolumeUniform`, `StepTheorem`, `UniformGapCertificate`, `CertificateImpliesGap` |
| `Statement/MassGap.lean` | `MassGapStatement`, `exponentialClustering_implies_clustering` |
| `Statement/NonTriviality.lean` | `IsGaussianTheory`, `NonTrivial`, `cumulantGF`, `nthCumulant`, `isGaussianTheory_cumulants_vanish`, `CumulantCharacterization`, `isGaussianTheory_iff_cumulants_vanish`, `toyTheory`, `toyTheory_not_gaussian`, `RGCharge`, `GaugeNonTrivial`, `gaugeNonTrivial_implies_not_gaussian` |

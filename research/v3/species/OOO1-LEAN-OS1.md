# OOO1-LEAN-OS1 — OS1 三块缺件的 Lean 声明与可写性判

> 轮 OOO1 · 子代理 `ooo1-lean-os1` · `D:\Dev\math`。写域：本文件 + `lean/os1/`（新建）。既有文件一字未改；未 `commit`；**零 Lean 编译**；零新 MC（未编译 ⇒ 本轮无一条「机器验证」断言，全部为源码级判读）。
> 口径：只 `[T]`/`[NG]`；无百分数；每结论带「下一手」。`AGENTS §1` 四要素见 §4。
> `[T]` 的判据（本轮口径，因按令不编译）：**使能件全部在 Mathlib（核对到 `文件:行`）+ 无需新定义 + 组装步数为常数**；`[NG]` = 缺件或缺定义。凡「名字待核」处已在代码内标注，它影响能否一次编译过，不影响 `[T]/[NG]` 分类。

---

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | **① 「一般 Polish 弱收敛 ⇒ tight」可写成，且是纯组装** | **`[T]`** | §1.①；六步全有：`Tendsto.isCompact_insert_range`(`Compact.lean:634`)、`IsCompact.isClosed`、`closure_minimal`、`IsCompact.of_isClosed_subset`(`:103`)、`isTightMeasureSet_of_isCompact_closure`(`Prokhorov.lean:609`)、`ProbabilityMeasure.t2Space`(`ProbabilityMeasure.lean:422`) | 把 OS1 骨架里的 `hTight` 假设**删掉**，改由 `isTightMeasureSet_range_of_tendsto` 供 |
| 0.2 | **NNN1 §2 「一般 Polish：`μ_a⇒μ ⇒ tight` 缺」须改一格** | **`[T]`**（更正） | §1.① 步 1–6。Mathlib 缺的是**Prokhorov 方向**（tight ⇒ 相对紧）；OS1 要的方向（收敛 ⇒ tight）= 「收敛 ⇒ 相对紧」+「相对紧 ⇒ tight」，**两步都在** | 台账把该行由「缺」改为「缺引理不缺件，≈6 行组装」 |
| 0.3 | **② 分两层：逐点 `(t,x)↦x+t` 是 `[T]`；测度层「移动试验函数」是 `[NG]`** | **`[T]` / `[NG]`**（分账） | §1.②a：`continuous_add`、`Homeomorph.addRight`；②b：全库 0 命中（§4 检索式 ⑤），Mathlib 只有**固定** `f` 版（`ProbabilityMeasure.lean:639`） | 按 §1.②b 补齐 ~80 行 ε-论证（材料齐：`HeineCantor.lean:47`、`IsCompact.prod`、`Tight.lean:55`） |
| 0.4 | **②b 不需要新定义，只需一次真证明** | **`[NG]`** | §1.②b 注释列出四件用料与余项拆法 | 写成后 `hJoint` 从假设降为定理 |
| 0.5 | **③ `𝒮′(ℝ⁴)` 实例：`[NG]`，且不是「漏写 instance」，是「现行载体不相容」** | **`[NG]`** | §1.③：`TemperedDistribution` = 点态拓扑（`TemperedDistribution.lean:57` + 作者自注「不 barrelled、UBP 失效」）；点态拓扑 = `{s｜Finite s}` 上的 `UniformConvergenceCLM`（`PointwiseConvergenceCLM.lean:60`）⇒ 不可数坐标积 ⇒ 无 `PseudoMetricSpace/SecondCountableTopology`；`CompleteSpace` 亦无 | 换强对偶（Fréchet–Montel 型）载体，或退到「可数试验族 + 等度连续」；两条都要**新定义** |
| 0.6 | **③ 里有两行 `[T]`**：`MeasurableSpace S' := borel S'`、`BorelSpace S' := ⟨rfl⟩` | **`[T]`** | §1.③ 顶部；但**不够**喂块 ①（缺的是度量/可数/完备） | 这两行可立即落地；不要据此宣布「`𝒮′` 实例到手」 |
| 0.7 | **OS1 主定理在抽象 Polish 载体上是 `[T]` 组装，条件 = `hConv`+`hInv`+`hJoint`** | **`[T]`** | §1 末 `os1_map_shift_eq_self`；分离步用 `FiniteMeasure.ext_of_forall_integral_eq`(`FiniteMeasure.lean:379`)，其 `HasOuterApproxClosed` 由 `PseudoMetrizableSpace` 自动给(`HasOuterApproxClosed.lean:218`) | 位置图记「OS1 Lean 化 = 组装 `[T]`，唯一 `[NG]` = ②b」 |

**一句话**：① 与骨架都是 `[T]`；② 半 `[T]` 半 `[NG]`；③ `[NG]` 且卡的是载体不是 instance。**OS1 的 Lean 化缺口收敛为一条**：②b。

---

## §1 三个 Lean 声明（正文见 `lean/os1/*.lean`）

### ① tight_of_weak_convergence —— `[T]`

文件：`lean/os1/OS1Tight.lean`

```lean
theorem isTightMeasureSet_range_of_tendsto
    {E : Type*} [PseudoMetricSpace E] [MeasurableSpace E] [BorelSpace E]
    [SecondCountableTopology E] [CompleteSpace E]
    {μs : ℕ → ProbabilityMeasure E} {μ : ProbabilityMeasure E}
    (h : Tendsto μs atTop (𝓝 μ)) :
    IsTightMeasureSet (Set.range fun n => (μs n : Measure E))
```

- **需要什么**：`E` Polish（`PseudoMetricSpace`+`SecondCountableTopology`+`CompleteSpace`）、`[BorelSpace E]`；`PseudoMetricSpace` 顺带给 `HasOuterApproxClosed`，故 `T2Space (ProbabilityMeasure E)` 到手。
- **Mathlib 有没有**：件全有（六步，逐行标在文件里）；**这条引理本身没有**（0 命中，见 §4 ⑥），但它就是六行组装。
- **缺的怎么补**：不缺。旁边附 `exists_compact_measure_compl_le_of_tendsto`（纯改写 `Tight.lean:60`）供 OS1 以 `∀ε ∃K` 消费。

### ② 联合推移连续 —— ②a `[T]` / ②b `[NG]`

文件：`lean/os1/OS1JointShift.lean`

```lean
-- ②a-1 [T]
theorem continuous_add_pair : Continuous fun p : E × E => p.1 + p.2
-- ②a-2 [T]
def addRightHomeo (t : E) : E ≃ₜ E := Homeomorph.addRight t
-- ②b-0 [T]  固定 f 版（白拿）
theorem tendsto_integral_shift_const_of_tendsto
    (h : Tendsto μs atTop (𝓝 μ)) (f : E →ᵇ ℝ) (y₀ : E) :
    Tendsto (fun n => ∫ x, f (x + y₀) ∂(μs n : Measure E)) atTop (𝓝 (∫ x, f (x + y₀) ∂(μ : Measure E)))
-- ②b [NG]  移动 f 版
theorem tendsto_integral_shift_of_tendsto
    (f : E →ᵇ ℝ) (hTight : IsTightMeasureSet (Set.range fun n => (μs n : Measure E)))
    (hμ : Tendsto μs atTop (𝓝 μ)) (hy : Tendsto y atTop (𝓝 y₀)) :
    Tendsto (fun n => ∫ x, f (x + y n) ∂(μs n : Measure E)) atTop (𝓝 (∫ x, f (x + y₀) ∂(μ : Measure E)))
```

- **需要什么**：②a 要 `[IsUniformAddGroup E]`（含 `ContinuousAdd`）；②b 要 `[UniformSpace E] [IsUniformAddGroup E] [T2Space E] [BorelSpace E] [PseudoMetrizableSpace E]` + `hTight`。
- **Mathlib 有没有**：②a 有（`continuous_add` 为 `continuous_mul` 的 to_additive；`Homeomorph.addRight`，方向由 `Group/Basic.lean:1081` 的 doc 引用确认）。②b **没有**：全库 0 命中；只有固定 `f` 的 `tendsto_map_of_tendsto_of_continuous`（`ProbabilityMeasure.lean:639`）。
- **缺的怎么补**：**不新定义**，补 ~80 行：`hTight` 取 `K`；`K × ({y_n}∪{y₀})` 紧（`IsCompact.prod`，`Compact.lean:1062`）；`(x,y) ↦ f (x+y)` 在 `K×S` 上一致连续（`HeineCantor.lean:47`，只要 `[UniformSpace]`，**不要度量**）；余项 `‖f‖·μ_n(Kᶜ)` 用 `restrict`/`setLIntegral` 拆。另注：这里**不需要**平移不变度量（对 `(x,y_n)` vs `(x,y₀)` 用积一致结构即可）——比 OS1 §1 S3 的表述少一个前提。

### ③ `𝒮′(ℝ⁴)` 实例 —— `[NG]`

文件：`lean/os1/OS1Schwartz.lean`。`E4 := EuclideanSpace ℝ (Fin 4)`，`S' := 𝓢'(E4, ℂ)`。

```lean
-- [T] 已有（Mathlib 现成）：AddCommGroup(UniformConvergenceCLM.lean:167)、Module ℂ(:238)、
--      IsUniformAddGroup(:197)、IsTopologicalAddGroup(:203)、T2Space(PointwiseConvergenceCLM.lean:70)
-- [T] 两行 Borel 记账
noncomputable instance instMeasurableSpaceS' : MeasurableSpace S' := borel S'
instance instBorelSpaceS' : BorelSpace S' := ⟨rfl⟩
-- [NG] 三条目标（不落 sorry：命题会假）
--   instance : PseudoMetricSpace S'      instance : SecondCountableTopology S'      instance : CompleteSpace S'
-- [NG] 平移作用（新定义）
def shift (t : E4) (T : S') : S' := T.comp (SchwartzMap.compSubConstCLM ℂ (-t))
theorem continuous_shift_pair : Continuous fun p : E4 × S' => shift p.1 p.2
theorem os1_at_schwartz (…) : …   -- 真结论 μ.map (shift y₀) = μ
```

- **需要什么**：`PseudoMetricSpace S'` + `SecondCountableTopology S'` + `CompleteSpace S'`（`Prokhorov.lean:545/605` 的变量块），及 `τ_t` 的联合连续 + 可测性。
- **Mathlib 有没有**：拓扑群/一致结构**有**；Polish/Borel 类的度量与可数性**没有**；`compSubConstCLM` 只有单点 CLM（`SchwartzSpace/Basic.lean:1068`），**没有**「对 `a` 联合连续」。
- **缺的怎么补（关键判）**：不是补 instance。点态拓扑是 `{s｜Finite s}` 上的 `UniformConvergenceCLM`（`PointwiseConvergenceCLM.lean:60`），坐标集 `𝓢(ℝ⁴,ℂ)` 不可数 ⇒ 无可数邻域基；`CompleteSpace` 又会碰上作者自注的「不 barrelled、UBP 失效」。要 Polish 只能**换载体**：强对偶（Fréchet–Montel）拓扑，或退到「可数试验族 + 等度连续」。两条都是新定义 ⇒ `[NG]`。故 OS1 主定理本轮只在抽象 Polish 载体上交付。

### 骨架（`[T]`，条件于 ②b）

文件：`lean/os1/OS1Statement.lean`

```lean
theorem probabilityMeasure_ext_of_forall_integral_eq (h : ∀ f : E →ᵇ ℝ, ∫ f ∂μ = ∫ f ∂ν) : μ = ν
theorem os1_map_shift_eq_self
    (hConv : Tendsto μs atTop (𝓝 μ)) (hy : Tendsto y atTop (𝓝 y₀))
    (hInv : ∀ n, (μs n).map (· + y n) = μs n)
    (hJoint : ∀ f : E →ᵇ ℝ, Tendsto (fun n => ∫ x, f (x + y n) ∂(μs n)) atTop (𝓝 (∫ x, f (x + y₀) ∂μ))) :
    μ.map (· + y₀) = μ
```

要点：`hTight` **已不是假设**（① 供）；`hJoint` = ②b 当显式假设（`hMinlos` 体例）；分离步 `FiniteMeasure.ext_of_forall_integral_eq` 白拿。`[T]`。

---

## §2 可写性判

| 块 | 判 | 是「组装」还是「新定义」 | 缺的精确内容 |
|---|---|---|---|
| ① tight | **`[T]`** | 纯组装（6 步） | 无；Mathlib 只是缺这条**六行引理** |
| ②a 逐点 | **`[T]`** | 纯组装（0 步，`exact`） | 无 |
| ②b 测度层 | **`[NG]`** | 组装 + 一次真证明（**不**新定义） | 80 行 ε-论证；四件用料齐但未接 |
| ③ `𝒮′` | **`[NG]`** | **需新定义**（换拓扑 / 新作用） | Polish 三件在**现定义下为假**；`τ_t` 联合连续缺（含试验函数侧） |
| 骨架 OS1 | **`[T]`** | 纯组装（条件于 ②b） | 无（`hJoint` 当假设） |

---

## §3 最小可编译子集（先写就能过，0 `sorry`）

**A. 声明级清单（按依赖序）**

1. `OS1Tight.isTightMeasureSet_range_of_tendsto` —— 0 假设外的新件，0 `sorry`。
2. `OS1Tight.exists_compact_measure_compl_le_of_tendsto` —— 依赖 1 + `Tight.lean:60`。
3. `OS1JointShift.continuous_add_pair` / `addRightHomeo` —— 0 依赖。
4. `OS1JointShift.tendsto_integral_shift_const_of_tendsto` —— 依赖 `ProbabilityMeasure.lean:346`。
5. `OS1Statement.probabilityMeasure_ext_of_forall_integral_eq` —— 依赖 `FiniteMeasure.lean:379`。
6. `OS1Statement.os1_map_shift_eq_self` —— 依赖 3、4、5（**把 ②b 当假设**）。
7. `OS1Schwartz.instMeasurableSpaceS'` / `instBorelSpaceS'` + 五条 `example : … := inferInstance`。

**B. 不进子集（因 `sorry` 或缺件）**：`tendsto_integral_shift_of_tendsto`（②b）、`continuous_shift_pair`（③）、`os1_at_schwartz`（③）。

**C. 构建前置（已核，必须记）**：`lean/` 下**没有** `lakefile.lean`、**没有** `lean-toolchain`（检索命中 0，§4 ⑬）；本机唯一带 mathlib 的 lake 工程是 `vendor/OSforGFF`（`lakefile.lean:17` 要求 mathlib `v4.33.0-rc1`，`lean-toolchain` = `leanprover/lean4:v4.33.0-rc1`）。⇒ 子集 A 目前**无编译目标**。

**下一手**：① 给 `lean/` 建 `lakefile.lean`（require mathlib `v4.33.0-rc1`）或把 `lean/os1/*.lean` 挂进 `vendor/OSforGFF` 的 `lean_lib`；② 编译 A；③ 再攻 ②b（80 行）；④ ③ 只在换载体后才谈。

---

## §4 检索核验（`AGENTS §1` 四要素）

| # | 工具状态 | 实际检索式 | 命中 |
|---|---|---|---|
| ① | `grep` **OK** | `theorem tendsto_iff_forall_integral_tendsto\|IsTightMeasureSet\|def IsTightMeasureSet` @ `Mathlib/MeasureTheory` | **4**（`FiniteMeasure.lean:726`、`ProbabilityMeasure.lean:346`、`Portmanteau.lean:688`、`Tight.lean:55`） |
| ② | `grep` **OK** | `isCompact_insert_range\|insert_range\|isCompact_of_tendsto` @ `Mathlib` | **18**（含 `Compact.lean:612/628/634`） |
| ③ | `grep` **OK** | `isTightMeasureSet_of_isCompact_closure` @ `MeasureTheory/Measure` | **3**（`Prokhorov.lean:609`、`Tight.lean:31`、`Prokhorov.lean:30`） |
| ④ | `grep` **OK** | `tendsto_map\|map_continuous\|Continuous fun μ…\.map` @ `Mathlib/MeasureTheory` | **8**（全为**固定** `f` 版：`ProbabilityMeasure.lean:639`、`FiniteMeasure.lean:957` 等） |
| ⑤ | `grep` **OK** | `map_add_of_tendsto\|tendsto_map_add\|jointly_continuous_map\|continuous_map_add\|tendsto_integral_add` @ 全 `Mathlib` | **0** |
| ⑥ | `grep` **OK** | `PolishSpace .*→SL\|PointwiseConvergenceCLM.*Polish\|MeasurableSpace .*→SLₚₜ\|instance .*TemperedDistribution.*(MeasurableSpace\|BorelSpace\|PseudoMetric\|Polish)` @ 全 `Mathlib` | **0** |
| ⑦ | `grep` **OK** | `compSubConstCLM\|subConst` @ `SchwartzSpace/Basic.lean` | **7**（`:1068,1073,1077,1082`；**无**「对 `a` 连续」） |
| ⑧ | `grep` **OK** | `^instance\|^noncomputable instance` @ `TemperedDistribution.lean` | **19**（全为 Fourier/LineDeriv 类；**无** `MeasurableSpace/BorelSpace/Polish`） |
| ⑨ | `grep` **OK** | `^instance` @ `UniformConvergenceCLM.lean` | **18**（`:167` AddCommGroup、`:197` IsUniformAddGroup、`:203` IsTopologicalAddGroup、`:238` Module） |
| ⑩ | `grep` **OK** | `HasOuterApproxClosed (α\|X\|Ω\|E)\|instance .*HasOuterApproxClosed` @ `Mathlib` | **16**（`HasOuterApproxClosed.lean:218`：`PseudoMetrizableSpace` ⇒ 实例） |
| ⑪ | `grep` **OK** | `uniformContinuousOn_of_continuous` @ `HeineCantor.lean` | **3**（`:47` 主体） |
| ⑫ | `grep` **OK** | `toFiniteMeasure_injective\|structure ProbabilityMeasure\|def ProbabilityMeasure (` @ `Mathlib` | **1**（`ProbabilityMeasure.lean:103`：`def` 包装的 subtype ⇒ 分离步只差一层 `ext`） |
| ⑬ | `glob`/`pwsh` **OK** | `lean/` 下 `lakefile*`、`lean-toolchain`；`vendor/OSforGFF/.lake/packages/mathlib/lean-toolchain` | lakefile@`lean/`：**0**；toolchain：**1**（`v4.33.0-rc1`） |

**失败标记**：`TOOL-FAIL` **0**。所有「缺」均为**已执行检索 0 命中**（⑤⑥）或**实例清单里查无**（⑧），非工具失败。**本轮零编译**，§1/§2 全为静态读源码（含 `文件:行`）。

---

## §5 产出文件

| 文件 | 内容 | 判 |
|---|---|---|
| `lean/os1/OS1Tight.lean` | 块 ① 两条声明（六步逐行标注） | `[T]`，0 `sorry` |
| `lean/os1/OS1JointShift.lean` | ②a-1、②a-2、②b-0（`[T]`）；②b（`[NG]`，1 个 `sorry`） | 半 `[T]` 半 `[NG]` |
| `lean/os1/OS1Schwartz.lean` | ③：已有实例子清单 + Borel 两行（`[T]`）+ Polish/`shift`/联合连续（`[NG]`） | 2 行 `[T]`，其余 `[NG]` |
| `lean/os1/OS1Statement.lean` | 分离引理 + OS1 主定理骨架（`hJoint` 当假设） | `[T]`，0 `sorry` |

**与 `nnn1-os1` 的接口**：骨架三件（`tendsto_iff_forall_integral_tendsto`、`Tight.lean:55`、`HeineCantor.lean:47`、`Homeomorph.addLeft`）**全部核实存在**（§4 ①③⑪、`Group/Basic.lean:1081`）。`nnn1` 说「缺两件 ≈40 + ≈100 行」：本轮的修正是——**tightness 那 40 行不必写**（组装 6 行，`[T]`），`hTight` 假设可直接删；剩下的 100 行全在 ②b（`[NG]`）。`𝒮′` 实体仍 `[NG]`，但理由从「实例不在 Mathlib」精确到「现行点态拓扑与 Polish 不相容」。

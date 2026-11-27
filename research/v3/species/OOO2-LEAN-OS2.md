# OOO2-LEAN-OS2 — OS2 缺的三块：`clipₙ` / `Ω` 的拓扑 / `θ_*μ = μ`

> 靶：`NNN2-OS2.md:51` 的「下一手（任一即升档）」三条（㈠ 截断换成 `clipₙ ∘ A`、㈡ 补 `θ_*μ=μ` 的 Lean 引理、㈢ 密度步）。
> 本轮只做㈠㈡ 的**代码**：`research/v3/species/ooo2/OS2Clip.lean`、`OS2Topology.lean`、`OS2Theta.lean`、`OS2Wiring.lean`（新文件，共 4 个）。
> 纪律：不改既有文件、不 commit、零新 MC、**不编译**；只交 `[T]`/`[NG]`，每条附「下一手」。
> 判据（NNN2-OS2.md §3）：`GL4Transfer.lean:48–66` 的 `rp_transfer` 真覆盖 `C_b` 层 `[T]`；缺三块 = ①`trunc` 不连续（`GL4L2.lean:28–29`）、②`GL4L2.lean:25` 无 `TopologicalSpace`、③仓内 `θ_*μ=μ` 定理 0 条。

---

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | **①** `clipₙ x = max(−n, min(n, x))` 写出，且 `Continuous (clip n)`、`|clipₙ x| ≤ |x|`、`clipₙ x → x`、`|clipₙ x| ≤ n`、`clipₙ x = x`（`|x|≤n`）、`clipₙ ∘ A : Ω →ᵇ ℝ`（**`A` 只要连续，不必有界**）、`clipₙ ∘ A ∈ L²` 全部写出 | **`[T]`** | `ooo2/OS2Clip.lean`；API 逐个核到 mathlib 源（§4 表） | 编译一次即升机检；`trunc` 可整体弃用 |
| 0.2 | **②** `Ω` 的拓扑给出：抽象层 4 个类换齐 GL4 全套实例（`OpensMeasurableSpace`/`HasOuterApproxClosed`/`T2Space (ProbabilityMeasure Ω)`/`CompactSpace (ProbabilityMeasure Ω)`）；具体层 `Ω := κ → [−1,1]`（Wilson 环坐标）、积拓扑/Tychonoff/`Pi.borelSpace`、`θ`=坐标置换且 `Continuous`+`≃ₜ`+`≃ᵐ` | **`[T]`** | `ooo2/OS2Topology.lean`；§2 | 把 `GL4L2.lean:25` 的 `variable` 段按文件末补丁替换 |
| 0.3 | **③** `MeasureTheory.Measure.map θ μ = μ` 的声明写出，并判清它**需要什么**：`θ` 连续 + 逐 `a` 不变 + 弱收敛 + `ProbabilityMeasure Ω` Hausdorff；**不需要** `Ω` 紧、tightness、Prokhorov、`a`-一致矩、`θ` 对合/双射 | **`[T]`** | `ooo2/OS2Theta.lean`；§3 | 把 `GL4L2` 的 `hθ : MeasurePreserving` 降为 `measurePreserving_of_weak_limit` 的推论 |
| 0.4 | **①②③ 并成一条**：`rp_full_cone` —— `A` 连续 + `L²(μ)`、`μ_a` 侧 clip 层 RP ⇒ 极限全锥 RP，且 `hθ` 不再是假设 | **`[T]`**（代码已写）／**`[NG]`**（机检） | `ooo2/OS2Wiring.lean` | 见 §5 编译配方，一次 `lean` 收口 |
| 0.5 | 四文件**机检**（`lean` 退出码 + `#print axioms` 脚印） | **`[NG]`** | 本轮硬规矩「不编译」；除「跑一次编译」外**无其它缺口** | 下一轮执行 §5 的两条命令 |

**一句话判词。** ① 把「有界连续层 ⇒ 截断层」的**类型错配**修掉了（`trunc` 不连续是唯一死因，`clipₙ` 连续且保留 `trunc` 的三条估计）；② 的拓扑只需 4 个类，具体模型就是 REPORT.md §3.1 的 `∏[−1,1]` 坐标空间，`θ`＝坐标置换自然是 `≃ₜ`/`≃ᵐ`；③ 的真内容是**不变性过弱极限**，它严格弱于 GL4 的其余两段，**不消费紧性也不消费 tightness**——所以仓内 `hθ : MeasurePreserving (⇑θ) μ μ` 是把结论当假设（其字段 `map_eq` 字面等于结论）。三块合起来只在签名里剩两条输入：`μ_a` 侧的 clip 层 RP、`A ∈ L²(μ)`。

---

## §1 ① `clipₙ` + 连续性引理（`ooo2/OS2Clip.lean`）

### 1.1 为什么必须换

`GL4L2.lean:28–29` 的 `trunc A n ω = if |A ω| ≤ n then A ω else 0`：

| 需要 | `trunc` | `clipₙ` |
|---|---|---|
| 逐点收敛（`rp_l2_closure` 的 `hpt`） | `trunc_tendsto` :37 ✓ | `clip_tendsto` ✓ |
| `‖·‖ ≤ ‖A ω‖`（`hbound`） | `trunc_abs_le` :31 ✓ | `clip_abs_le` ✓ |
| **有界连续（`rp_transfer` 要 `Ω →ᵇ ℝ`）** | ✗ **不连续** | ✓ `continuous_clip` |
| 有界 | ✗（截断值随 `A` 无界） | ✓ `clip_abs_le_n` |

`trunc` 的缺陷是**单一**的：不连续 ⇒ 进不了 `BoundedContinuousFunction`；三条估计里它已有两条。

### 1.2 定义与引理

```lean
noncomputable def clip (n : ℕ) (x : ℝ) : ℝ := max (-(n : ℝ)) (min (n : ℝ) x)

lemma continuous_clip (n : ℕ) : Continuous (clip n) :=
  continuous_const.max (continuous_const.min continuous_id)
```

配套（同文件）：`clip_abs_le : |clip n x| ≤ |x|`、`clip_norm_le : ‖clip n x‖ ≤ ‖x‖`、
`clip_abs_le_n : |clip n x| ≤ n`、`clip_eq_self : |x| ≤ n → clip n x = x`、`clip_tendsto : clip n x → x`、
`clip_memLp : A ∈ L² ⇒ clipₙ ∘ A ∈ L²`。

### 1.3 送进 `rp_transfer` 的槽位

```lean
noncomputable def clipBCF (n : ℕ) (A : Ω → ℝ) (hA : Continuous A) : Ω →ᵇ ℝ :=
  BoundedContinuousFunction.ofNormedAddCommGroup (fun ω => clip n (A ω))
    ((continuous_clip n).comp hA) (n : ℝ) (…)

@[simp] lemma clipBCF_apply … : clipBCF n A hA ω = clip n (A ω) := rfl
lemma clip_integrand_eq … : (clipBCF n A hA) (θ ω) * (clipBCF n A hA) ω
                          = clip n (A (θ ω)) * clip n (A ω) := rfl
```

**关键增益（写清不吃亏）**：`clipₙ ∘ A` 是 `Ω →ᵇ ℝ` 只要求 `A` **连续**，不要求 `A` 有界 ——
这正是 `C_RP ∩ L²(μ)` 里那些**无界**连续 `A` 需要的形状。`trunc` 做不到：它既不连续，也不有界。

### 1.4 残余（写进签名，不藏）

`clip` 修的是类型/连续性的错配，**不是** `μ_a` 侧的 RP 本身。`OS2Wiring.rp_full_cone` 的假设

```
hRP : ∀ n i, 0 ≤ ∫ ω, clip n (A (θ ω)) * clip n (A ω) ∂(μs i)
```

承载了「`clipₙ ∘ A` 在 `μ_a` 上满足 RP」。两种读法：

- **RP 是对全体 `C_b(Ω)` 的不等式**（NNN2-OS2.md §1 的读法，`C_RP := {μ : ∀A∈C_b, ∫θA·A dμ ≥ 0}`）⇒ `clipₙ ∘ A ∈ C_b`，该假设**免费**；
- **RP 只在环多项式代数上**（非线性截断不在代数里）⇒ 该假设是一条**额外输入**。

两种读法下 `trunc` 都不行：不连续 ⇒ 连 `C_b` 都不属于，而 `clipₙ` 至少落在 `C_b`。

**下一手**：若采第二种读法，需在 `REPORT.md:164` 的接缝处补一句「`C_RP` 的定义域取 `C_b(Ω)` 而非环多项式代数」，否则 §1.4 的残余会以「假设换假设」形式重现。

---

## §2 ② `Ω` 的拓扑（`ooo2/OS2Topology.lean`）

### 2.1 抽象层：四个类

```lean
[TopologicalSpace Ω] [PseudoMetrizableSpace Ω] [MeasurableSpace Ω] [BorelSpace Ω]
```

由此**自动**得到（文件里各写一条 `example … := inferInstance` 钉死，不留口头承诺）：

| 要的实例 | 谁给 | mathlib 源 |
|---|---|---|
| `OpensMeasurableSpace Ω`（`rp_transfer` 要） | `BorelSpace.opensMeasurable` | BorelSpace/Basic.lean:183 |
| `HasOuterApproxClosed Ω`（`ProbabilityMeasure` 的 T2 要） | `PseudoMetrizableSpace` 全局实例 | HasOuterApproxClosed.lean:236 附近 |
| `T2Space (ProbabilityMeasure Ω)`（极限唯一） | `ProbabilityMeasure.t2Space` | ProbabilityMeasure.lean:422 |
| `CompactSpace (ProbabilityMeasure Ω)`（极限存在） | Prokhorov | Prokhorov.lean:168（节内变量含 `[T2Space E] [BorelSpace E]`，:65） |

注意：**只加 `OpensMeasurableSpace` 不够**（`GL4Transfer.lean:39` 现有那段就是它），
T2 与推前连续性都要 `BorelSpace`。这是②相对现状的实质增量。

### 2.2 具体层：REPORT.md §3.1 的公共空间就是现成模型

```lean
abbrev LoopCoord : Type := Set.Icc (-1 : ℝ) 1        -- 归一化 Wilson 环坐标
abbrev LoopSpace (κ : Type*) : Type := κ → LoopCoord  -- Ω := ∏_{i∈κ} [−1,1]

def coordSwap (σ : κ ≃ κ) (x : LoopSpace κ) : LoopSpace κ := fun i => x (σ i)  -- θ
lemma continuous_coordSwap (σ : κ ≃ κ) : Continuous (coordSwap σ) :=
  continuous_pi fun i => continuous_apply (σ i)
noncomputable def coordSwapHomeo (σ : κ ≃ κ) : LoopSpace κ ≃ₜ LoopSpace κ
noncomputable def coordSwapMeasEquiv [Countable κ] (σ : κ ≃ κ) : LoopSpace κ ≃ᵐ LoopSpace κ
```

实例来源：积拓扑（`Pi` 自动）、紧= Tychonoff（`Function.compactSpace` + `compactSpace_Icc`）、
拟可度量化（`pseudoMetrizableSpace_pi`，**不需可数**）、`T2Space`（`Pi` 实例）、
`BorelSpace`（`Pi.borelSpace`，**要 `[Countable κ]`** + 纤维第二可数，`Subtype.secondCountableTopology` 给）。
可数性不是造假：取可数多个 Wilson 环做坐标即可（REPORT.md §3.1 只用到「有限联合」的推前）。

`θ` 的自同胚性与对合性**分开**：`coordSwap` 对任意 `σ : κ ≃ κ` 都连续且是 `≃ₜ`；
反射给的是 `σ = σ⁻¹`，那是物理侧的额外事实，**代码不需要**。

### 2.3 对 `GL4L2.lean:25` 的补丁（不改既有文件）

```
-- 原: variable {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]
variable {Ω : Type*} [TopologicalSpace Ω] [PseudoMetrizableSpace Ω]
  [MeasurableSpace Ω] [BorelSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]
```

**下一手**：打上这段后，`GL4L2` 的 231 行可原地用 `clip` 版替换 `trunc` 版（`rp_l2_closure_trunc` → 见 §5 的 wiring），不必重写证明。

---

## §3 ③ `θ_*μ = μ`：声明 + 它要什么（`ooo2/OS2Theta.lean`）

### 3.1 声明

```lean
-- 主定理
theorem map_eq_of_weak_limit
    (θ : Ω → Ω) (hθ : Continuous θ)
    {ι : Type*} {F : Filter ι} [NeBot F]
    (μs : ι → ProbabilityMeasure Ω) (μ : ProbabilityMeasure Ω)
    (hweak : Tendsto μs F (𝓝 μ))
    (hinv : ∀ i, Measure.map θ (μs i : Measure Ω) = (μs i : Measure Ω)) :
    Measure.map θ (μ : Measure Ω) = (μ : Measure Ω)

-- 同一事实的 MeasurePreserving 形式：直接填 GL4L2.lean:87 的 hθ 槽位
theorem measurePreserving_of_weak_limit
    (θe : Ω ≃ᵐ Ω) (hθc : Continuous (⇑θe)) … :
    MeasurePreserving (⇑θe) (μ : Measure Ω) (μ : Measure Ω)
```

证明三行：`ProbabilityMeasure.tendsto_map_of_tendsto_of_continuous`（推前弱连续）→
`hinv` 把左端识别为 `μs` 本身 → `tendsto_nhds_unique`（Hausdorff）得 `μ.map = μ`，再 `toMeasure_map` 落地。

### 3.2 它需要什么（逐条，含反例口）

| 输入 | 是什么 | 能否去掉 |
|---|---|---|
| `Continuous θ` | 推前 `ν ↦ ν.map θ` 弱连续的唯一结构条件 | **不能**：`θ` 不连续则 `f∘θ ∉ C_b`，弱收敛测试族被破坏 |
| `hweak : μs ⇀ μ` | Mathlib 真弱收敛（`tendsto_iff_forall_integral_tendsto`） | **不能**：这是「过极限」的全部机制 |
| `hinv : ∀ i, θ_*μ_a = μ_a` | 逐 `a` 的格点不变性 = **物理输入** | **不能**：F4-4 的失效口（反射面/分块不对齐）先坏这里；与 os3 的 `(B1)`（`A ∈ W ∩ RWR⁻¹`）同源 |
| `[PseudoMetrizableSpace Ω] [BorelSpace Ω]` | `ProbabilityMeasure Ω` 的 Hausdorff 性 ⇒ 极限唯一 | **不能**（但这是②给的，非额外物理假设） |

**不需要**（这是本条唯一真正的新信息）：`Ω` 紧、tightness、Prokhorov、`sup_a` 任何阶矩、
`θ` 是对合、`θ` 是双射。逐条对照 `NNN2-OS2.md §2` 的「不需要一致矩界」——③ 把那条表述**再削弱一档**：
连**极限存在**都不消费（存在性归 OS0/GL3，不动）。

### 3.3 对 GL4L2 的后果

`GL4L2.lean:74,87,211` 的三处 `hθ : MeasurePreserving (⇑θ) μ μ` 全是**假设**，而
`MeasurePreserving` 的字段 `map_eq` 字面就是 `Measure.map (⇑θ) μ = μ`（`map_eq_of_measurePreserving` 一行证完）。
即：**当前 Lean 链条在这一步是把结论当假设**。③ 把它降级为推论 ⇒ 全锥链条的显式输入从 5 条降到 2 条
（`hRP` 的 clip 层、`MemLp A 2 μ`；原 `hθ`、`hAR` 皆可导出）。

**下一手**：`REPORT.md:215`／`GGG2-COMPACT.md:29` 的「唯一显式假设＝`trunc A n ∈ L²`」应改写为
「**2 条**显式假设：`μ_a` 侧 clip 层 RP、`A ∈ L²(μ)`；`θ` 测度保持是**定理**（`OS2Theta`）」。
另：`GGG3-OSCLOSED.md:19` 的「`C_RP` 弱闭凸锥」标点修（概率测度集不是锥）已在 `NNN2-OS2.md:21` 记，本条不改。

---

## §4 检索核验（AGENTS §1 四要素）

| 项 | 值 |
|---|---|
| 工具状态 | `pwsh` OK · `grep` OK · `glob` OK · `read` OK · `write` OK（**无 `TOOL-FAIL`**） |
| 实际检索式 | ① `MeasurePreserving\|Measure\.map\|map_eq` @ `research/glue-a-limits/*.lean`；② `TopologicalSpace\|PseudoMetrizableSpace\|BorelSpace` @ 同目录；③ `theorem MeasurePreserving\.map_eq`、`tendsto_map_of_tendsto_of_continuous`、`toMeasure_map`、`toMeasure_injective`、`tendsto_nhds_unique`、`abs_le`、`neg_abs_le`、`Continuous.max\|Continuous.min`、`ofNormedAddCommGroup`、`MemLp.of_le\|MemLp.mono`、`Continuous.aestronglyMeasurable`、`Pi.borelSpace`、`Subtype.secondCountableTopology`、`compactSpace_Icc`、`Function.compactSpace`、`pseudoMetrizableSpace_pi`、`continuous_pi\|continuous_apply`、`Homeomorph` 字段、`Homeomorph.toMeasurableEquiv` @ `vendor/OSforGFF/.lake/packages/mathlib/Mathlib/**` |
| 命中 | ①`MeasurePreserving` **3**（`GL4L2.lean:74,87,211`，全为假设行：`Measure.map` 定理 **0**）；②`TopologicalSpace` **1**（`GL4Transfer.lean:39`；`GL4L2.lean` 内 **0** ← 即②的靶）；③上述 mathlib 名逐个命中，行号已写进各 `.lean` 文件头与 §2.1 表 |
| 失败标记 | 无 |

补充（本次新查得的两条关键 API，均为 `[T]` 级引用，行号见文件头）：
`ProbabilityMeasure.tendsto_map_of_tendsto_of_continuous`（ProbabilityMeasure.lean:639）——③ 的全部引擎；
`BoundedContinuousFunction.ofNormedAddCommGroup`（Bounded/Normed.lean:118）——① 的 `Ω →ᵇ ℝ` 构造器。

---

## §5 交付物与编译配方（本轮**未执行**）

| 文件 | 内容 | 行数级 |
|---|---|---|
| `research/v3/species/ooo2/OS2Clip.lean` | ① `clipₙ` + 7 条引理 + `clipBCF` + 交接缝 | ~150 |
| `research/v3/species/ooo2/OS2Topology.lean` | ② 抽象 4 类 + 具体 `κ → [−1,1]` + `θ`（连续/`≃ₜ`/`≃ᵐ`） | ~140 |
| `research/v3/species/ooo2/OS2Theta.lean` | ③ 平凡版 + 主定理 + `MeasurePreserving` 形式 | ~120 |
| `research/v3/species/ooo2/OS2Wiring.lean` | ①②③ 并一条 `rp_full_cone` | ~80 |

编译（下一轮执行；本轮被「不编译」禁止）：

```bash
ROOT="D:/Dev/math/vendor/OSforGFF/.lake/packages"
MP="$ROOT/mathlib/.lake/build.badcache/lib/lean"
export LEAN_PATH="$MP;research/glue-a-limits;research/v3/species/ooo2"
~/.elan/bin/elan run leanprover/lean4:v4.33.0-rc1 lean research/v3/species/ooo2/OS2Wiring.lean
```

（`OS2Wiring.lean` 同时 `import GL4Transfer GL4L2 OS2Clip OS2Theta`，故它一条命令覆盖四个文件的解析；
`#print axioms` 的脚印检查可照 `research/glue-a-limits/check-glue-a.sh:39–43` 的写法加。）

**下一手优先级**：① 跑上面这条命令（`[NG]`→`[T]` 的唯一缺口）；② 若 `OS2Topology` 的
`example … := inferInstance` 有任一不过，按 §2.1 表的 mathlib 源逐条替换为显式 `⟨…⟩` 构造；
③ 通过后把 `GL4L2.lean:25` 的 `variable` 段按 §2.3 补丁升级、弃用 `trunc`。

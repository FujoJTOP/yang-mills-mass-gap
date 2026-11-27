# OOO3-LEAN-FIXED — `CompactTransfer` 与 `Δ_a > 0`（Lean 接口 + 主定理骨架 + 有/无判定）

> 第 OOO3 轮。只新建本文件与 `ooo3/CompactTransfer.lean`、`ooo3/LambdaOne.lean`；`TransferData` 零改动；未 commit；零新 MC；**未编译**（无 `lake build`、无 `#check`、无 `#eval`）。
> 口径：固定 `(L,a)`、Wilson 作用、紧群 `G`、`0 ≤ β < ∞`、`a > 0`；`T_a = transfer a`、`Ω = vacuum`、`λ₀ = 1`、`λ₁ = ‖T_a|_{Ω^⊥}‖`、`Δ_a = −(1/a)·log λ₁`。标签只用 `[T]` / `[NG]`。

---

## §0 判词

| # | 问 | 判 | 下一手 |
|---|---|---|---|
| 1 | `CompactTransfer` 四字段写出（含类型） | **`[T]`**（接口级；未编译，见 §3 的 olean 硬约束） | 字段里 `letI` 的写法要过一遍真编译（`CompactTransfer.lean:110–130`） |
| 2 | 主定理声明 + 证明骨架（四字段 →「顶端单重」→ `λ₁<1`） | **`[T]`**（§2；谱链每一步都点得出 Mathlib 名字与行号） | 把 SORRY 0/1/2/3/4/6/9/10/11/12 十处 A 级胶水补成真证明（≈110 行） |
| 3 | `Δ_a > 0` 的**无条件**闭合（固定 `(L,a)`） | **`[NG]`**：数学上 `[T]`（NNN4 §1 已判），Lean 上缺两条 B 级 | B 级第 1 步：造 `T_a`（`Lattice` 时间分层，800–1500 行） |
| 4 | 用 `IsCompactOperator`、避开 Hilbert–Schmidt | **`[T]`**（本轮复核：Mathlib 8268 文件，`HilbertSchmidt` 命中 **0**；`IsCompactOperator` 命中 **4** 文件） | 不再投人；HS 路线判死 |
| 5 | 四输入不可去 | **`[T]`**（承接 `FFF4-LAM1.md:24` 的 5×5 见证；`quadPos` 在 §2 的 `μ ≥ 0` 一步被用掉，缺它 `λ₁` 可负） | 写进 `CompactTransfer` 的 doc（已完成） |
| 6 | **修正 NNN4 §1 第 7 行** | 四字段对**每个** `a>0` 各给一个 `Δ_a`，**不给** `SemigroupGap`（全体 `t≥0` 共用一 `Δ`） | 要 `SemigroupGap` 需再加「`t ↦ λ₁(t)` 半群插值 + cutoff 一致的 `Δ₀`」——正是真靶，不在四字段内 |
| 7 | 结构形态裁决 | 不 `extends TransferData`，改为 `Prop` 结构 + `T` 作参数 | 否则 `compact` 字段类型找不到 `NormedAddCommGroup T.H`（理由见 §1 末） |
| 8 | `kernel_pos / row_normalized` 的位置 | 留在四字段；**但**"`T_a` 是 `κ` 的积分算子"必须**单列**为 `RepresentsKernel` | 否则四字段只是关于 `κ` 的两句话，与 `transfer` 无连接（§1 末，这是本轮最重要的一条） |

**一句话**：接口与谱链 `[T]`；无条件闭合 `[NG]`，缺的是 `T_a` 的 `L²(μ)` 表示与极值原理，不是谱定理。

---

## §1 `CompactTransfer`（structure，四字段）

文件：`research/v3/species/ooo3/CompactTransfer.lean`

```lean
structure CompactTransfer (T : TransferData) (X : Type*) [MeasurableSpace X]
    (μ : Measure X) [IsProbabilityMeasure μ] (κ : ℝ → X → X → ℝ) : Prop where
  compact : letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ a : ℝ, 0 < a → IsCompactOperator (T.transfer a)
  quadPos : letI := T.instNormed; letI := T.instInner; letI := T.instComplete
    ∀ a : ℝ, 0 ≤ a → ∀ f : T.H, 0 ≤ (inner ℂ f (T.transfer a f)).re
  kernel_pos : ∀ a : ℝ, 0 < a → ∀ x y : X, 0 < κ a x y
  row_normalized : ∀ a : ℝ, 0 < a → ∀ x : X, ∫ y, κ a x y ∂μ = 1
```

字段与 FFF4 见证的对应：

| 字段 | 数学内容 | 去掉之后 | 在证明里哪一步用 |
|---|---|---|---|
| `compact` | `T_a` 紧 | `λ₁ = 1` 可发生 | §2 第 5 步（取本征值） |
| `quadPos` | `0 ≤ T_a`（Loewner） | `λ₁` 可负，`Δ_a` 无定义 | §2 第 6 步（`μ ≥ 0`） |
| `kernel_pos` | `κ_a > 0` | `λ₁ = 1`（常数向量） | §2 第 4 步（极值原理） |
| `row_normalized` | `T_a 1 = 1` | `λ₁ = 1`（概率向量） | §2 第 4 步（极值原理） |

### 三处形态裁决（都会咬人，逐条给理由）

1. **不 `extends`。** `TransferData` 的 `instNormed / instInner / instComplete` 是**普通字段**、不是 instance binder（`OS.lean:99–102`）。若 `CompactTransfer extends TransferData`，`compact` 的字段类型里 `T.transfer a` 需要 `NormedAddCommGroup T.H`，而 `toTransfer.instNormed` 不在 instance 搜索里、`letI` 又没法写在 `extends` 的字段类型前，除非每个字段自己 `letI := toTransfer.instNormed`。改成 `T` 作参数、字段类型自带 `letI`，等价且干净；四字段全是 `Prop`，故整个结构取 `Prop`（这是"假设包"，不是数据）。
2. **`λ₁ := ‖T_a|_{Ω^⊥}‖`（压缩算子范数），不是"第二大本征值"。** 这样定义不需要本征向量基、不需要 HS、不需要 `λ₁` 存在性定理；对正紧自伴算子它**等于** `Ω^⊥` 上谱的上确界（`gapOp_norm_lt_one` 的证明里就是用谱半径取的），所以是同一个量的无坐标形式。
3. **`RepresentsKernel` 必须单列。** `kernel_pos`、`row_normalized` 只说 `κ` 的性质，**没有任何字段把 `κ` 和 `T_a` 连起来**。连起来的那句话是"`T_a` 是 `κ_a` 的积分算子、`Ω` 是常函数 `1`"，它需要 `H ≅ L²(μ)` 的模型——本仓 `T_a` 根本不存在（`Gauge/Lattice.lean:38–42` 只有 `V / L / Plaq`，无时间层）。故单列为

```lean
def RepresentsKernel (T : TransferData) (X : Type*) [MeasurableSpace X] (μ : Measure X)
    (κ : ℝ → X → X → ℝ) : Prop :=
  letI := T.instNormed; letI := T.instInner; letI := T.instComplete
  ∀ a : ℝ, 0 < a → ∀ ψ : T.H, T.transfer a ψ = ψ →
    (∀ v ∈ vacuumLine T.vacuum, inner ℂ ψ v = 0) →
    ∃ u : X → ℝ, Integrable u μ ∧ (∫ x, u x ∂μ = 0) ∧
      (∀ᵐ x ∂μ, u x = ∫ y, κ a x y * u y ∂μ) ∧
      ((∀ᵐ x ∂μ, u x = 0) → ψ = 0)
```

最后一条 `(u = 0 a.e. → ψ = 0)` 是**忠实性**，不可省：没有它，`u ≡ 0` 推不出 `ψ = 0`，极值原理就断在最后一步。

---

## §2 主定理：声明与证明骨架

文件：`research/v3/species/ooo3/LambdaOne.lean`

### 声明链（逐级变弱，每级都可单独编译）

```lean
-- ① 谱侧核心（抽象 Hilbert 空间，与 TransferData 无关）
theorem gapOp_norm_lt_one (A : E →L[ℂ] E) (hAsym : A.IsSymmetric) (hAk : IsCompactOperator A)
    (hA1 : ‖A‖ ≤ 1) (hApos : ∀ x, 0 ≤ (inner ℂ x (A x)).re)
    (Ω : E) (hΩ : A Ω = Ω) (hΩne : Ω ≠ 0)
    (hsimple : ∀ ψ, A ψ = ψ → ψ ∈ (vacuumLine Ω)ᗮ → ψ = 0) :
    ‖gapOp A hAsym Ω hΩ‖ < 1

-- ② 四字段 → 顶端本征值单重
theorem TransferData.top_eigenspace_trivial (CT : CompactTransfer T X μ κ)
    (hrep : RepresentsKernel T X μ κ) (ha : 0 < a) :
    ∀ ψ, T.transfer a ψ = ψ → (∀ v ∈ vacuumLine T.vacuum, inner ℂ ψ v = 0) → ψ = 0

-- ③ λ₁ < 1
theorem TransferData.lambdaOne_lt_one (CT) (hrep) (ha : 0 < a) : T.lambdaOne a ha.le < 1

-- ④ Δ_a > 0
theorem TransferData.Delta_pos (CT) (hrep) (ha : 0 < a) : 0 < T.Delta a ha.le

-- ⑤ 单步 decay（四字段真正给的东西；不是 SemigroupGap）
theorem TransferData.stepGap (CT) (hrep) (ha) : T.StepGap a (T.Delta a ha.le)
```

### 证明骨架（①的八步，逐步骤对 Mathlib）

```
假设 1 ≤ ‖S‖（S := gapOp A …）
  1. ‖S‖ = ‖A‖ 不成立，只知 ≤；用不着            -- 兜底：压缩不增范（SORRY 3）
  2. spectralRadius ℂ S = ‖S‖₊                    -- Rayleigh.lean:182  spectralRadius_eq_nnnorm
  3. ∃ μ ∈ spectrum ℂ S, ‖μ‖₊ = spectralRadius     -- Normed/Algebra/Spectrum.lean:229
                                                    --    （谱紧 :150 + 非空 spectrum.nonempty）
  4. μ ≠ 0                                        -- ‖μ‖ = ‖S‖ ≥ 1
  5. HasEigenvalue S μ                            -- FredholmAlternative.lean:221
  6. 取本征向量 ψ ≠ 0；‖μ‖ ≤ 1 ⇒ ‖μ‖ = 1            -- 本征方程 + 压缩不增范
  7. μ 实（conj_eigenvalue_eq_self）+ μ ≥ 0（A 正）⇒ μ = 1
  8. 抬回 E：A ψ = ψ、ψ ⊥ Ω、ψ ≠ 0 ⇒ 与 hsimple 矛盾   -- 矛盾 → ‖S‖ < 1
```

②的两步（B 级）：

```
ψ 不动 + ψ ⊥ Ω
  -- hrep 给出实值 u：均值 0、满足 u = ∫κ u
  9. 极值原理：正核 + 行归一 ⇒ u 是 a.e. 常值        -- fixed_point_eq_const（SORRY 8）
 10. 常值 + 均值 0 ⇒ u = 0 a.e. ⇒ ψ = 0            -- hrep 的忠实性
```

④的算术（容易踩的一脚）：

```lean
rw [TransferData.Delta]; split_ifs with h0
· norm_num                       -- λ₁ = 0 分支：钉死值 1（因为 Real.log 0 = 0！）
· exact mul_pos_of_neg_of_neg hinv (Real.log_neg hpos hq)
```

`Real.log 0 = 0` 是 Lean 的约定，`−(1/a)·log 0 = 0`，会给出 `Δ_a = 0` 而非 `> 0`。这是本轮唯一一处"数学上无所谓的边界在 Lean 里会直接假"的地方，必须在定义里钉死。

---

## §3 哪几步 Mathlib 有、哪几步要自建

### 有（本地 `vendor/OSforGFF/.lake/packages/mathlib`，逐条名字 + 行号已核）

| 用途 | 名字 | 位置 |
|---|---|---|
| 紧算子定义 | `IsCompactOperator` | `Analysis/Normed/Operator/Compact/Basic.lean:71`（`Prop`，故可进 `Prop` 结构） |
| 子空间压缩仍紧 | `IsCompactOperator.restrict'` / `.restrict` | 同上 `:349` / `:336`（**以 `M →ₗ M` 陈述**，CLM 需来回） |
| 紧自伴谱定理 | `ContinuousLinearMap.orthogonalComplement_iSup_eigenspaces_eq_bot` | `Analysis/InnerProductSpace/Spectrum.lean:443` |
| 紧自伴本征空间有限维 | `ContinuousLinearMap.finite_dimensional_eigenspace` | 同上 `:463` |
| 谱半径 = 范数（自伴） | `ContinuousLinearMap.spectralRadius_eq_nnnorm` | `Analysis/InnerProductSpace/Rayleigh.lean:182` |
| 非零谱点 = 本征值（紧） | `IsCompactOperator.hasEigenvalue_iff_mem_spectrum` | `Compact/FredholmAlternative.lean:221` |
| 谱上取到谱半径 | `spectrum.exists_nnnorm_eq_spectralRadius_of_nonempty` | `Analysis/Normed/Algebra/Spectrum.lean:229` |
| 谱紧 | `spectrum.isCompact` | 同上 `:150` |
| Rayleigh 极值即本征向量 | `IsSelfAdjoint.hasEigenvector_of_isMaxOn` | `Rayleigh.lean:276`（备选路线） |
| 闭子空间完备 | `Submodule.isClosed_orthogonal`、`IsClosed.completeSpace_coe` | `InnerProductSpace/Orthogonal.lean:114` |
| 本征值实（自伴） | `IsSymmetric.conj_eigenvalue_eq_self` | 已在 `Spectrum.lean:106` 用过 |
| CLM 压缩 | `ContinuousLinearMap.restrict` | `Topology/Algebra/Module/ContinuousLinearMap/Restrict.lean:31` |
| Loewner 正定 | `ContinuousLinearMap.IsPositive` | `InnerProductSpace/Positive.lean:266`（本轮未用，`quadPos` 直接给二次型） |

### 无 / 要自建

| 缺什么 | 证据 | 自建代价 |
|---|---|---|
| Hilbert–Schmidt | 三式检索 `HilbertSchmidt` / `Hilbert–Schmidt` / `hilbertSchmidt` @ Mathlib 8268 文件 → 各命中 **0** | 不需要：路线已避开 |
| 积分算子紧性（`integralOperator` 是紧） | `IsCompactOperator` @ `MeasureTheory/**` → 命中 **0** | B 级，500–1200 行（或走 Ascoli，300–800 行） |
| `T_a` 本身（时间分层 + 核 + `T_a 1 = 1`） | `Gauge/Lattice.lean:38–42` 无时间层 | 800–1500 行 |
| 极值原理（正核 + 行归一 ⇒ 常值） | 无现成引理 | 200–400 行（`essSup` + 单调收敛拼装） |
| `IsSelfAdjoint → IsSymmetric`（CLM 版） | 线性映射版 `Adjoint.lean:710` 在，CLM 版名字本机无法 `#check` | 3–10 行 |
| 压缩不增范 `‖S‖ ≤ ‖A‖` | 未查到现成引理 | 10 行 |
| 半群插值 ⇒ 全体 `t` 的 `SemigroupGap` | 无 | **这是真靶**，不在四字段内 |

### 最小可编译子集（本轮实测，不是推断）

环境实测：Mathlib 只有 **2264 / 8268** 个 `.olean`；`Analysis/InnerProductSpace` 只有 **3** 个（`Basic/Defs/Minimal`）；`MeasureTheory` 只有 **6** 个。`YangMillsInfra` 12 个 `.olean` 全在（含 `Reconstruction/OS.olean`、`Spectral/Transfer.olean`）。

| 本文件需要的模块 | olean |
|---|---|
| `Mathlib.Analysis.InnerProductSpace.Basic` | ✅ |
| `Mathlib.Topology.Algebra.Module.ContinuousLinearMap.Restrict` | ✅ |
| `Mathlib.Analysis.Normed.Operator.Basic` | ✅ |
| `Mathlib.Analysis.SpecialFunctions.Log.Basic` | ✅ |
| `YangMillsInfra.Reconstruction.OS` | ✅ |
| `Mathlib.Analysis.InnerProductSpace.Orthogonal` | ❌ |
| `Mathlib.Analysis.InnerProductSpace.Adjoint` | ❌ |
| `Mathlib.Analysis.InnerProductSpace.Symmetric` | ❌ |
| `Mathlib.Analysis.InnerProductSpace.Spectrum` | ❌ |
| `Mathlib.Analysis.InnerProductSpace.Rayleigh` | ❌ |
| `Mathlib.Analysis.Normed.Operator.Compact.Basic` | ❌ |
| `Mathlib.Analysis.Normed.Operator.Compact.FredholmAlternative` | ❌ |
| `Mathlib.Analysis.Normed.Algebra.Spectrum` | ❌ |
| `Mathlib.MeasureTheory.Integral.Bochner.Basic` | ❌ |
| `Mathlib.MeasureTheory.Measure.ProbabilityMeasure` | ❌ |

**结论**：`IsCompactOperator`、`ᗮ`、`∫`、`IsSymmetric`、`spectralRadius` 全都落在未构建的模块里。今天 0 build 能编译的只有**标量骨架**（`vacuumLine` 的 `span`、`Delta` 的 `if/log`、`SemigroupGap` 的消费端）。要编译本文件，最小集 = 上表 **10 个 ❌ 根模块**（其 import 闭包未测，禁止本轮 build）：

```
Mathlib.Analysis.InnerProductSpace.Orthogonal
Mathlib.Analysis.InnerProductSpace.Adjoint
Mathlib.Analysis.InnerProductSpace.Symmetric
Mathlib.Analysis.InnerProductSpace.Spectrum
Mathlib.Analysis.InnerProductSpace.Rayleigh
Mathlib.Analysis.Normed.Operator.Compact.Basic
Mathlib.Analysis.Normed.Operator.Compact.FredholmAlternative
Mathlib.Analysis.Normed.Algebra.Spectrum
Mathlib.MeasureTheory.Integral.Bochner.Basic
Mathlib.MeasureTheory.Measure.ProbabilityMeasure
```

编译入口（下一轮用，本轮不跑）：在 `vendor/OSforGFF` 里 `lake env lean <绝对路径>/CompactTransfer.lean`（先 `lake build` 上述 10 个模块）。

---

## §4 `sorry` 台账（标签 0–12 共 13 个，其中 7 在本轮已消去 ⇒ **现存 12 处**）

| 标签 | 位置 | 内容 | 级别 | 估行 |
|---|---|---|---|---|
| 0 | `CompactTransfer.lean` `gapOp` | `A` 在 `ℂ·Ω` 上 = id | A | 5 |
| 1 | `LambdaOne.lean` §A | `IsCompactOperator S`（CLM/线性映射来回） | A | 10 |
| 2 | §A | `S.IsSymmetric`（`restrict_invariant` CLM 包装） | A | 10 |
| 3 | §A | `‖μ‖ ≤ 1`（压缩不增范 + 本征方程） | A | 10 |
| 4 | §A | `μ = 1`（实 + 非负 + `‖μ‖=1`） | A | 25 |
| 5 | `CompactTransfer.lean` | `IsSelfAdjoint → IsSymmetric`（CLM） | A | 5 |
| 6 | §A | `A ↑ψ = ↑ψ`（`restrict_apply`） | A | 5 |
| 7 | ~~`LambdaOne.lean` §A~~ | ~~`ψ ⊥ Ω`~~ **本轮已消去**：抽象引理的 `hsimple` 直接取成员形式 `ψ ∈ (vacuumLine Ω)ᗮ`，最后一步就是 `ψ.property` | — | 0 |
| 8 | §B | 极值原理 `fixed_point_eq_const` | **B** | 200–400 |
| 9 | §A | `p = ⊥` 时 `‖S‖ = 0` | A | 3 |
| 10 | §D | `T_a` 与 `gapOp` 在 `Ω^⊥` 上一致 | A | 5 |
| 11 | §D | `exp(−Δ_a·a) = λ₁`（`Real.exp_log` + 分支） | A | 15 |
| 12 | §D | `n·a` 的归纳 decay | A | 15 |

A 级合计 ≈ **110 行**（`sorry` 8 除外），B 级 ≈ **200–400 行**，另加 `RepresentsKernel` 的模型侧 800–1500 行。

---

## §5 工具状态（AGENTS §1 四要素）

| 工具 | 状态 | 检索式 | 命中 |
|---|---|---|---|
| `pwsh` | OK（列目录、取行、`Test-Path` 判 olean、`Select-String` 计数） | `Get-ChildItem … -Filter *.lean` 后 `Measure-Object` @ Mathlib | 8268 文件 |
| `pwsh` | OK | `Select-String "HilbertSchmidt"`、`"Hilbert–Schmidt"`、`"hilbertSchmidt"` @ Mathlib（三式各跑） | **0** / **0** / **0** |
| `pwsh` | OK | `Select-String "IsCompactOperator"` @ Mathlib | **4** 文件（`-List` 口径；NNN4 报的 132 是匹配行口径，两者不冲突） |
| `pwsh` | OK | `Select-String "integralOperator"` 与 `"IsCompactOperator"` @ `MeasureTheory/**` | **0** |
| `pwsh` | OK | `Test-Path …\build\lib\lean\<mod>.olean` ×15 | 5 ✅ / 10 ❌（§3 表） |
| `pwsh` | OK | `(Get-ChildItem … -Filter *.olean).Count` @ Mathlib build | **2264** |
| `grep` / `glob` | OK | `**/OS.lean`、`research/v3/species/**/*.lean`、`**/Spectrum.lean` | 2 / 2 / 12 |
| `read` | OK | `OS.lean`（全 140 行）、`Spectral/Transfer.lean`（全 107 行）、`NNN4-FIXED.md`（全 52 行）、`Spectrum.lean:95–165`、`Rayleigh.lean:141–320`、`Compact/Basic.lean:305–360`、`FredholmAlternative.lean:196–229` | — |

无 `TOOL-FAIL`。本文件所有"0 命中"均附检索式与工具状态。**未编译**是任务约束，不是工具失败；§3 的"能不能编译"用 olean 存在性判定，属静态判定，不是编译通过记录。

---

## §6 下一手（三条，按顺序）

1. **build 10 个根模块**（§3 列表），然后 `lake env lean CompactTransfer.lean`：先验 `letI` 写在结构字段类型里能否过（本轮三处形态裁决里唯一有编译风险的一条）。失败则退回"每个字段加一个 `letI` 助手 `def` + 字段类型引用该 `def`"的写法。
2. **补 A 级 6 处胶水**（SORRY 1/2/3/4/6/7，≈110 行）→ `lambdaOne_lt_one` 与 `Delta_pos` 变成 `sorry`-free；此时"固定 `(L,a)` 的 `Δ_a>0`"在 Lean 里是**条件定理**（条件 = `RepresentsKernel` + 四字段）。
3. **B 级两件**：先做 `fixed_point_eq_const`（200–400 行，只依赖 `MeasureTheory` 已有零件）；再做 `T_a` 的时间分层（800–1500 行）——后者做完，`RepresentsKernel` 才可从模型推出，`Δ_a>0` 才从条件定理升为无条件。

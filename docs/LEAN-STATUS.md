# LEAN-STATUS — 条件归约链的机器检验台账

> 2026-11-23 · 机器：Windows 11，x86_64-w64-windows-gnu。
> 本台账**诚实**记录每个命名结果的状态 ∈ {PROVEN（0 sorry）/ STATED（n sorry + 精确缺失内容）}。
> **Clay 质量间隙未证；`H1` 始终是假设。** 本文件检验的是 **"假设 ⇒ 结论"的逻辑有效性**，不是物理结论。

## 0. 环境与构建方式

| 项 | 值 |
|---|---|
| Lean 工具链 | `Lean (version 4.33.0-rc1, x86_64-w64-windows-gnu, commit 62eed1db4d67327ec8120be05f1a1b0847d74561, Release)` |
| 机器上另一套 | `leanprover/lean4:v4.36.0`（**未使用**；vendored Mathlib olean 是 4.33.0-rc1 构建的） |
| Mathlib | 随 `vendor/OSforGFF/.lake/packages/mathlib` 提供的 **预编译 olean 缓存**（`build.badcache`），版本 `v4.33.0-rc1` |
| 构建方式 | `lean/` 目录**不是 lake 工程**（无 `lakefile`）；沿用现有约定：用 `elan run` 直接调 `lean`，`LEAN_PATH` 指向上面的 olean 缓存。项目"build"即下面两个脚本。 |

> 为什么不用 `lake build`：`lean/` 无 `lakefile`；`research/lean-infra/` 虽有 lakefile，但其设计注释明确说明本会话是按 olean 缓存逐模块直接编译（"a full `lake build` was out of scope"）。完整 `lake build` 会触发整份 Mathlib 重编，超出范围且非本任务所需。**Lean/Mathlib 本身是安装好的**（oleans 可用，文件全部编译通过），不是环境缺失。

## 1. 新文件：`lean/YMStatement/Reduction.lean`（namespace `YMStatement.Reduction`）

### (A) 完全证明 · FULLY PROVEN（0 sorry）

| 命名结果 | 状态 | 内容 / 依据 |
|---|---|---|
| `feshbach_bogoliubov_schur` | **PROVEN** | Feshbach–Bogoliubov 双块 Schur 引理：`A ≻ 0`、`Σ = D − Bᵀ A⁻¹ B ⪰ 0` ⇒ `[[A,B],[Bᵀ,D]] ⪰ 0`。经 Mathlib `Matrix.PosDef.fromBlocks₁₁`（两向等价）导出。 |
| `gDouble` / `gDoubleOf` | def | `G″ = E[H] − Var + W` 的曲率分解量。 |
| `gDouble_eq` | **PROVEN**（`rfl`） | 定义性恒等式 `gDouble EH Var W = EH − Var + W`。 |
| `gDouble_eq_S_sub_D` | **PROVEN** | 等价形式 `G″ = S − D`，`S=E[H]`、`D=Var−W`（`ring`）。 |
| `gDouble_nonneg` | **PROVEN** | 胶合：`Var ≤ E[H] ∧ 0 ≤ W ⇒ 0 ≤ G″`。 |
| `gDouble_free` | **PROVEN** | 自由场酸测：`G″(0,0,0)=0`（工具在自由场上失效，非 vacuous）。 |
| `gdouble_lower_bound` | **PROVEN** | `W_geo≥0 ⇒ E[H] − Var ≤ G″`（`G″` 分解的下界一侧，G2.0 等价形式）。 |
| `div_le_of_le_mul` | **PROVEN** | 下侧：`f ≤ χ·lower ⇒ f/χ ≤ lower`。 |
| `le_div_of_mul_le` | **PROVEN** | 上侧：`χ·upper ≤ M₀ ⇒ upper ≤ M₀/χ`。 |
| `susceptibility_squeeze` | **PROVEN** | 两步挤压 `f₁²/χ ≤ inf supp ν ≤ M₀/χ`。 |
| `m_sq_ge_f_lam` | **PROVEN** | `0<f ∧ χ ≤ 1/λ ⇒ f·λ ≤ f/χ`。 |
| `variance_susceptibility_squeeze` | **PROVEN** | 精确陈述：`0<f ∧ f≤χ·lower ∧ χ≤1/λ ⇒ f·λ ≤ m=f/χ ≤ lower`。 |
| `m_sq_ge_of_squeeze` | **PROVEN** | REDUCTION 充分条件版：`r₀>0 ∧ χ≤1/λ ∧ r₀≤f₁² ⇒ r₀·λ ≤ f₁²/χ`。 |
| `transported_residue` | **PROVEN** | `P3 ∧ CP ⇒ 0 < ρ ∧ ρ ≤ f₁²`，`ρ = Z_O² r₀`（乘性残数搬运）。 |

### (C) 归约定理 · THE REDUCTION THEOREM（0 sorry）

| 命名结果 | 状态 | 内容 |
|---|---|---|
| `mass_gap_of_spectral` | **PROVEN** | `H1 ∧ CP ∧ P3 ∧ 无溶解 ⇒ m_gap>0`（经 A 段挤压；不调用任何 sorry 定理）。 |
| **`clay_gap_conditional`** | **PROVEN（0 sorry）** | **`H1 ∧ (CP) ∧ (W_geo≥0) ∧ P3 ∧ 无溶解 ⇒ MassGap`**，其中 `MassGap := 0 < mgap`。这是机器检验的核心：只证**逻辑**，未证输入全是**假设**。`W_geo≥0` 用于给出 `G″ ≥ E[H]−Var` 的下界一侧。 |

### (B) 以 sorry 陈述的未证解析输入 · STATED

| 命名结果（行号） | 状态 | 精确缺失内容 |
|---|---|---|
| `H1_holds`（L209） | **STATED（1 sorry）** | 粗模 Poincaré / 方差上界 `Var_ν(φ) ≤ (1/λ)∫\|∇φ\|²`（等价 `χ≤1/λ`）。缺：`(L,a)` 一致的 `∇_∥`（`k→0`）Poincaré 常数下界 = 质量间隙本身（MASTER §6.1、PROOF-SKELETON G2.5）。 |
| `CP_holds`（L215） | **STATED（1 sorry）** | (CP) 的乘性搬运内核 `f₁² = Z_O²·f₁²_strong`，`Z_O>0`。缺：跨 Savvidy 交叉的 Kato/隔离 + 连续极限（CONTINUUM-CP）。 |
| `W_geo_nonneg_fact`（L222） | **STATED（1 sorry）** | `W_geo = E[Φ″] − Var(Φ′) ≥ 0`。缺：Haar 纤维上有效作用的 **Prékopa 对数凹/凸性**全 β 无条件版（PROOF-JACOBIAN-SAVES 仅 β=0 闭合 + 数值）。 |
| `P3_strong_coupling`（L229） | **STATED（1 sorry）** | `β<1/48 ⇒ 0<r₀ ∧ r₀≤f₁²_strong`。缺：graded-strip KP 引理 + `a₀>0` 的**一致余项**（PROOF-P3-UNCONDITIONAL）；且仅到 `β<1/48`，未搬到 `β=2.5`。 |
| `no_dissolution_kinematic`（L237） | **STATED（1 sorry）** | 无溶解运动学裕度 `m₁²>0 ⇒ m_gap>0`。缺：连续 `a→0` 下把边原子零点与物理谱边**一致地**等同（PROOF-NO-DISSOLUTION-KINEMATIC 只给固定 `a`/固定盒版本）。 |

**`(B)` 的 5 条定理不会被 `(C)` 调用**——`clay_gap_conditional` 直接把它们的结论当作假设。因此 `(C)` 的
公理足迹不含 `sorryAx`（见 §3）。

## 2. 既有条件链（STEP 1：只确认、未改逻辑）

| 文件 | 命名结果（示例） | 状态 |
|---|---|---|
| `lean/YMStatement/Handoff.lean` | `handoff_gap`、`mass_gap_of_certificate` | **PROVEN**（编译通过） |
| `lean/YMStatement/RGHandoff.lean` | `lam_const`、`qfun_strictAnti`、`crossing_exists_unique`、`coverage`、`gap_lower_bound`、`certificate_transfer` 等 | **PROVEN**（编译通过） |
| `lean/YMStatement/Conditional.lean` | `tUp2_le_tUp`、`q0_le_q0TwoLoop`、`ClayInputs.twoLoop_overlap`、`gap_of_inputs`、`clay_of_inputs`、`clustering_rate_of_gap` | **PROVEN**（`#guard_msgs` 冻结公理足迹，编译通过） |
| `research/handoff-proof/HandoffMargin.lean` | `handoff_iff_margin` 等 | **PROVEN**（编译通过） |
| `research/glue-a-limits/GL9Decay.lean` | `GlueA.spectral_decay` | **PROVEN**（公理 `[propext, Classical.choice, Quot.sound]`） |

## 3. 构建命令与输出

### 命令 1（既有链）
```
bash lean/check-conditional.sh
```
输出摘要（exit 0）：
```
=== lean/YMStatement/Handoff.lean (build olean)
=== research/handoff-proof/HandoffMargin.lean (build olean)
=== research/glue-a-limits/GL9Decay.lean (build olean)
'GlueA.spectral_decay' depends on axioms: [propext, Classical.choice, Quot.sound]
=== lean/YMStatement/Conditional.lean
CONDITIONAL COMPILES (axioms frozen in #guard_msgs)
```

### 命令 2（本任务新文件）
```
bash lean/check-reduction.sh
```
输出（exit 0）：
```
=== lean/YMStatement/Reduction.lean
lean/YMStatement/Reduction.lean:209:8: warning: declaration uses `sorry`
lean/YMStatement/Reduction.lean:215:8: warning: declaration uses `sorry`
lean/YMStatement/Reduction.lean:222:8: warning: declaration uses `sorry`
lean/YMStatement/Reduction.lean:229:8: warning: declaration uses `sorry`
lean/YMStatement/Reduction.lean:237:8: warning: declaration uses `sorry`
'YMStatement.Reduction.feshbach_bogoliubov_schur' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.gDouble_nonneg' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.gdouble_lower_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.susceptibility_squeeze' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.variance_susceptibility_squeeze' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.m_sq_ge_of_squeeze' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.transported_residue' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.mass_gap_of_spectral' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.clay_gap_conditional' depends on axioms: [propext, Classical.choice, Quot.sound]
```

**sorry 计数**：`lean/YMStatement/Reduction.lean` 共 **5** 处，全部在 `(B)` 段
（L209 `H1_holds`、L215 `CP_holds`、L222 `W_geo_nonneg_fact`、L229 `P3_strong_coupling`、
L237 `no_dissolution_kinematic`）。`(A)` 与 `(C)` 段 **0 sorry**——由上面的 `#print axioms`
输出确认：不含 `sorryAx`，只有标准公理 `[propext, Classical.choice, Quot.sound]`。

## 4. 诚实判词

- **归约逻辑已机器检验**：`clay_gap_conditional` 严格证明
  **`H1 ∧ (CP) ∧ (W_geo≥0) ∧ P3 ∧ 无溶解 ⇒ m_gap>0`**，0 sorry，只含标准公理。
  即"若这些解析输入成立，则质量间隙成立"这条**逻辑蕴含**是有效的。
- **未证的是输入，不是逻辑**：5 条 `(B)` 输入（H1、(CP) 搬运内核、`W_geo≥0`、P3 下界、
  无溶解连续版）仍是 **sorry 的假设**，其精确缺失内容已逐条命名。
  **`H1` 仍然是假设，绝不称为已证。**
- **Clay 本身**：**未证**。本台账不改变这一事实。

---
*源文件*：[`lean/YMStatement/Reduction.lean`](../lean/YMStatement/Reduction.lean)、
[`lean/check-reduction.sh`](../lean/check-reduction.sh)、既有链
[`lean/YMStatement/Conditional.lean`](../lean/YMStatement/Conditional.lean)。

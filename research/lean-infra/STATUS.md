# STATUS — lean-infra

> 2026-11-22 · Lean 4.33.0-rc1 · mathlib v4.33.0-rc1 (`79d0395`)
> 交付状态：**设计 + 可编译骨架**。12 个模块、1,570 行、0 sorry、0 自定义 axiom，逐文件编译通过。
> 本轮增补（来自并行线）：`Statement/NonTriviality.lean`（nontriviality）与 `RG/GroupUniform.lean`（group-uniformity）。

## 一眼看板

| 子系统 | 文件 | 状态 | 说明 |
|---|---|---|---|
| 共享词汇 | `Basic.lean` | ✅ 编译 | `SpaceTime`, `TestFunction`, `Configuration`, `EuclideanTheory`, `TranslationAction`, `ExponentialClustering`, `Clustering` |
| 格点规范场 | `Gauge/Lattice.lean` | ✅ 编译 + **已证** | `MulAction` 实例、`holonomy_gaugeAct`、`wilsonAction_invariant`（任意群、任意类函数） |
| Wilson 测度 | `Gauge/Measure.lean` | ✅ 编译 | `wilsonMeasure` 是定义；`wilsonMeasure_univ` 已证；规范不变性剩 `haar_invariant`（G-3） |
| 轨道测度 | `Gauge/Orbits.lean` | ✅ 编译 + **已证** | 推前等式 `integral_orbitPushforward`；`GaugeFixing` 记录商测度缺口（G-2） |
| 无界算子 | `Spectral/Unbounded.lean` | ✅ 编译 | 二次型 gap 定义；`inner_apply_im_zero` 已证 |
| 转移衰减 | `Spectral/Transfer.lean` | ✅ 编译 + **已证** | `twoPoint_bound`、`decay_at_nat_multiples` |
| OS 重构接口 | `Reconstruction/OS.lean` | ✅ 编译 | OS 形式、正半定接口、`TransferData`（真空、半群、自伴、压缩） |
| 交接引理 | `RG/Handoff.lean` | ✅ 编译 + **已证** | `gap_uniform`、`regions_cover`、`overlap_iff_order` |
| 群参数化交接 | `RG/GroupUniform.lean` | ✅ 编译 + **已证** | `GroupHandoffData.gap_uniform`、`uniform_over_groups`（无隐藏 G 依赖）、`Ledger` |
| RG 步进 | `RG/Step.lean` | ✅ 编译 | `gapCertificate`、`StepTheorem`、`UniformGapCertificate`、`CertificateImpliesGap` |
| 陈述层 | `Statement/MassGap.lean` | ✅ 编译 + **已证** | 速率聚类 ⇒ OS4 裸聚类 |
| 非平凡性 | `Statement/NonTriviality.lean` | ✅ 编译 + **已证** | `IsGaussianTheory`、`nthCumulant`、Lemma C 充分方向、`toyTheory_not_gaussian`、`RGCharge`、`GaugeNonTrivial` |

## 编译方式（本 session 实际做的）

- **不做**全量 `lake build`（任务要求）。
- 用 `vendor/OSforGFF/.lake/packages/mathlib/.lake/build.badcache/lib/lean`（8,279 个兼容 olean）作为 `LEAN_PATH`，`check.sh` 逐文件编译，全部 exit 0。
- 未导入 OSforGFF / BochnerMinlos / GaussianField；本包 mathlib-only。
- 源码中无 `sorry`、无 `axiom` 声明（已 grep 确认；`axiom` 仅出现在注释里）。

## 已证 vs 假设

**已证（real theorems）**
1. `wilsonAction_invariant` — 格点 Wilson 作用规范不变。
2. `integral_orbitPushforward` — 规范不变观测量的期望 = 轨道空间上的期望。
3. `gap_uniform` / `regions_cover` / `overlap_iff_order` — 弱→强交接引理及其等价形式。
4. `exponentialClustering_implies_clustering` — 速率形式强于 OS4。
5. `SemigroupGap.twoPoint_bound` / `decay_at_nat_multiples` — 衰减传递到矩阵元与多步。
6. `wilsonMeasure_univ` — 归一化 Wilson 测度是概率测度（条件于 Z 有限非零）。
7. `inner_apply_im_zero` — 自伴算子的二次型实值。
8. `isGaussianTheory_cumulants_vanish` — **Lemma C 充分方向**：Gaussian ⇒ `κₙ = 0` (n ≥ 3)。
9. `toyTheory_not_gaussian` — 具体相互作用玩具理论非 Gaussian（**反例测试通过**）。
10. `gaugeNonTrivial_implies_not_gaussian` — 非零 `RGCharge` 排除 Gaussian。
11. `uniform_over_groups` — 四个常数一致 ⇒ gap 界在群族上一致，证明不触碰群。

**假设（以结构字段命名，非 sorry）**
- `GaugeFixing.integral_eq` ← 紧群商测度分解（G-2）
- `WilsonGaugeInvariance.haar_invariant` ← 紧群 Haar 右不变性（G-3）
- `OSPositiveSemidefinite.cauchy_schwarz` ← OS 形式的极化不等式（可证，~50 行）
- `CumulantCharacterization.gaussian_of_vanishing` ← Lemma C 必要方向（需 OS0 解析性）
- `CertificateImpliesGap.transfer_gap` ← 间隙证书 ⇒ 转移矩阵隙
- `StepTheorem` 的全部定量字段 ← Balaban/BBS 分析的输出

## 下一步（按优先级）

1. **P0 树规范固定**：`∀ U, ∃! g, g • U ∈ treeSlice T`，用 `IsFundamentalDomain.mk'`（不需可数性）→ 填 `GaugeFixing`。1–2 周。
2. **P0 紧群 Haar 右不变性**：补 G-3，约 50–100 行。1 周。
3. **P0 OS 形式 Cauchy–Schwarz**：完成 `OSHilbert` 的商 + 完备化构造。2–3 周。
4. **P1 适配层** `Bridge/OSforGFF.lean`：`EuclideanTheory.ofMeasure` + `SatisfiesAllOS` 桥接；以及 `IsGaussianTheory` ↔ `YMStatement.IsGaussianMeasure` 的 defeq 桥接。1 周。
5. **P1 `SU(N)` 实例**：`CharacterData` 的 `Re tr / N`、Haar、配分函数正性。2–3 周。
6. **P1 间隙证书 ⇒ 转移矩阵隙**（有限维/紧算子情形，用 `GaussianField.SpectralTheorem`）。3–4 周。
7. **P2 group-uniformity N2/N3**：2D 精确解 + Casimir scaling 的有限维截断（受阻于 G-16/G-19 的 Peter–Weyl）。

## 已知问题 / 注意

- `Basic.ExponentialClustering` 与 `YMStatement/Basic.lean` 的 `ExponentialClustering` 是同一命题的两个写法；`Basic.IsGaussianTheory`/`RGCharge` 与 `YMStatement.IsGaussianMeasure`/`RGCharge` 同理。整合时应统一，并用 `example` 断言 defeq。
- `OSforGFF.OS.Axioms.OS4_Clustering` 与 `Basic.Clustering` 同义；建议加 `Clustering_iff_OS4` 适配引理。
- `CharacterData` 目前不含 `χ` 的连续性/可测性；`Gauge/Measure` 的配分函数正性证明需要补 `χ_continuous`（或 `χ_measurable`）字段。
- mathlib 的 `haarMeasure` 在 `MeasureTheory.Measure` namespace，参数是 `TopologicalSpace.PositiveCompacts`；`haarMeasure_self` 需 `K₀ = univ` 才给出总质量 1。
- `Complex.exp_injective` 在本 mathlib 版本不存在；用 `Complex.ofReal_exp` + `Real.exp_injective` 替代（`NonTriviality.lean` 已按此写）。
- `Statement/NonTriviality.lean` 的 `cumulantGF` 取 `log Z` 的实部；对反射正/对称理论这就是通常的累积量，见 REPORT §3.11 的约定说明。
- `.olean` 文件（由 `check.sh` 生成）留在目录里以加速增量检查；它们不是源码，删除后 `check.sh` 会重建。
- 缓存里个别 mathlib 模块（如 `Mathlib.LinearAlgebra.LinearIndependent`、`Mathlib.Analysis.SpecialFunctions.Complex.Exp`）缺 olean；避免直接 import 它们。

## 文件

见 `REPORT.md` §1 与 §附录；缺口清单见 `REPORT.md` §5（G-1 … G-19）。

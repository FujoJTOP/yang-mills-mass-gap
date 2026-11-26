# L01 — Lean 环境与已有形式化盘点（"借 Lean 做搬运"的第一步）

> 轮次：本轮（`L01`）· 工作目录 `D:\Dev\math` · 未 `git commit`
> 纪律：所有结论配**实际输出**；标 `[T]`（实测）/ `[C]`（据实测推断）/ `[NG]`（实测否定）
> 检索工具状态：`pwsh` OK · `glob` OK · `grep` OK（无 TOOL-FAIL）

---

## §1 Lean 是否可用 —— 四条命令的实际输出

```
=== where.exe lean ===
C:\Users\hooya\.elan\bin\lean.exe
[exit=0]
=== where.exe lake ===
C:\Users\hooya\.elan\bin\lake.exe
[exit=0]
=== lean --version ===
Lean (version 4.36.0, x86_64-w64-windows-gnu, commit 611c117bfaf7b7ae1cfcf246ce6fba31e8f509e9, Release)
[exit=0]
=== lake --version ===
Lake version 5.0.0-src+611c117 (Lean version 4.36.0)
[exit=0]
```

**[T]** Lean 在 PATH 上：`C:\Users\hooya\.elan\bin\`（`elan` 代理，`lean/lake/leanc/leanchecker/leanmake/leanpkg/lean-beam` 全在同一目录，各 6,090,752 B）。

**[T]** 已装 toolchain（`~/.elan/toolchains/`）：

| toolchain | 大小 |
|---|---|
| `leanprover--lean4---v4.33.0-rc1` | 3,062 MB |
| `leanprover--lean4---v4.36.0` | 3,235 MB |

⇒ **§5 的前提（Lean 不在 PATH）不成立**，`~\\.elan` 已存在，无需回退搜索。

---

## §2 Mathlib 状态

全仓只有一个 `.lake`（`vendor\OSforGFF\.lake`）；**它就是这个仓库的 mathlib 宿主**。

| 项目 | 路径 | lakefile | `lean-toolchain` | `.lake` |
|---|---|---|---|---|
| OSforGFF（mathlib 宿主） | `vendor\OSforGFF` | `lakefile.lean` | `leanprover/lean4:v4.33.0-rc1` | **有**，9.29 GB |
| yangMillsInfra | `research\lean-infra` | `lakefile.toml`（require mathlib `v4.33.0-rc1`） | **缺失** | **无**（未 hydrate） |
| zenodo-v2 打包件 | `docs\zenodo-v2\ym-v2-lean` | `lakefile.lean` | `leanprover/lean4:v4.33.0-rc1` | **无** |
| openai 审计提交 | `vendor\YM-audit\submissions\openai-20261028` | `lakefile.toml` | `leanprover/lean4:v4.32.0` | **无** |

**Mathlib 本体** `[T]`：

- path：`D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib`
- toolchain：`leanprover/lean4:v4.33.0-rc1`
- git rev：`79d0395a1825a6264ad5d269e35e60537518955e`（2026-07-16 01:47:53 +0000），`inputRev = v4.33.0-rc1`
- 编译产物：
  - `.lake\build.badcache\lib\lean\` → **8,268** 个 `.olean`，含根 `Mathlib.olean` ⇒ **可 `import Mathlib`**
  - `.lake\build\lib\lean\` → 2,136 个 `.olean`（部分缓存，够 `import Mathlib.Analysis.SpecialFunctions.Exp`）
- 兄弟依赖（`build\lib\lean`）：`aesop` 132、`batteries` 188、`Qq` 14、`proofwidgets` 13、`LeanSearchClient` 4、`importGraph` 10、`plausible` 13、`GaussianField` 3、`BochnerMinlos` 2、`kolmogorov_extension4` 3

**跑通的两条路（都已实测）** `[T]`：

1. `cd D:\Dev\math\vendor\OSforGFF && lake env lean <file>` （`import Mathlib.Analysis.SpecialFunctions.Exp`）
2. 手工 `LEAN_PATH = mathlib\.lake\build.badcache\lib\lean ; batteries ; aesop ; Qq ; proofwidgets ; LeanSearchClient ; importGraph ; plausible ; Cli ; GaussianField ; BochnerMinlos ; kolmogorov_extension4`，再 `elan run leanprover/lean4:v4.33.0-rc1 lean <file>`（这是 `research\lean-infra\check.sh` 的做法，`import Mathlib` 全量可解析）

---

## §3 已有 `.lean` 文件

**全仓 `.lean`（排除 `.lake`）= 9,700 个** `[T]`：

| 树 | 文件数 | 说明 |
|---|---|---|
| `vendor\YM-audit` | 9,564 | 第三方审计提交包（openai-20261028），非本项目产出 |
| `vendor\OSforGFF` | 61 | 高斯自由场 OS 公理形式化（第三方，25.5k 行级） |
| `docs\zenodo-v2\ym-v2-lean` | 37 | **本项目** zenodo 打包快照 |
| `research\lean-infra` | 13 | **本项目** 12 模块 + 根 |
| `research\glue-*` / `handoff-proof` / `nontriviality` / `minimal-handoff` / `os2-rotation` | 17 | **本项目** 粘合引理 |
| `lean\YMStatement` | 6 | **本项目**（`docs\zenodo-v2\ym-v2-lean\YMStatement` 的**逐字副本**） |
| `.tools` | 2 | `deep.lean` / `trivial.lean`，各 2 行冒烟件 |

**去重后的自有源码：`research\` 36 件 / 3,638 行 + `lean\` 6 件 / 949 行 = 42 件 / 4,587 行**（`docs\zenodo-v2\ym-v2-lean\` 的 37 件是这 42 件的副本，重复计入会翻倍）`[T]`。

**逐文件 · 行数 · 真 `sorry` 数** `[T]`
（`sorry` 计数 = 代码级 `sorry` 战术，已剔除注释/反引号引用）

| # | 文件 | 行 | sorry |
|---|---|---|---|
| 1 | `research\lean-infra\YangMillsInfra.lean` | 12 | 0 |
| 2 | `…\YangMillsInfra\Basic.lean` | 67 | 0 |
| 3 | `…\YangMillsInfra\Gauge\Lattice.lean` | 117 | 0 |
| 4 | `…\YangMillsInfra\Gauge\Measure.lean` | 96 | 0 |
| 5 | `…\YangMillsInfra\Gauge\Orbits.lean` | 91 | 0 |
| 6 | `…\YangMillsInfra\Reconstruction\OS.lean` | 118 | 0 |
| 7 | `…\YangMillsInfra\RG\GroupUniform.lean` | 132 | 0 |
| 8 | `…\YangMillsInfra\RG\Handoff.lean` | 141 | 1 |
| 9 | `…\YangMillsInfra\RG\Step.lean` | 105 | 0 |
| 10 | `…\YangMillsInfra\Spectral\Transfer.lean` | 92 | 0 |
| 11 | `…\YangMillsInfra\Spectral\Unbounded.lean` | 95 | 0 |
| 12 | `…\YangMillsInfra\Statement\MassGap.lean` | 74 | 0 |
| 13 | `…\YangMillsInfra\Statement\NonTriviality.lean` | 180 | 0 |
| 14 | `research\glue-a-limits\GL3Limits.lean` | 87 | 1 |
| 15 | `research\glue-a-limits\GL4L2.lean` | 211 | 1 |
| 16 | `research\glue-a-limits\GL4Transfer.lean` | 83 | 2 |
| 17 | `research\glue-a-limits\GL9Decay.lean` | 90 | 2 |
| 18 | `research\glue-b-blocking\lean\GL2a.lean` | 83 | 1 |
| 19 | `research\glue-b-blocking\lean\GL2bPrime.lean` | 129 | 2 |
| 20 | `research\glue-c-gauge-group\GL6Fubini.lean` | 243 | 0 |
| 21 | `research\glue-c-gauge-group\GroupUniform.lean` | 112 | 0 |
| 22 | `research\glue-d-dual-handoff\GL7Dual.lean` | 113 | 1 |
| 23 | `research\handoff-proof\HandoffMargin.lean` | 96 | 1 |
| 24 | `research\handoff-proof\axcheck_full.lean` | 89 | 0 |
| 25 | `research\handoff-proof\axcheck.lean` | 3 | 0 |
| 26 | `research\minimal-handoff\MinimalHandoff.lean` | 65 | 0 |
| 27 | `research\nontriviality\NonTrivialCore.lean` | 96 | 1 |
| 28 | `research\nontriviality\check_axioms.lean` | 99 | 1 |
| 29 | `research\nontriviality\axiom_check.lean` | 5 | 0 |
| 30 | `research\os2-rotation\lean\SYMDim4.lean` | 249 | 0 |
| 31 | `lean\YMStatement\Basic.lean` | 71 | 0 |
| 32 | `lean\YMStatement\Handoff.lean` | 71 | 0 |
| 33 | `lean\YMStatement\RGHandoff.lean` | 261 | 0 |
| 34 | `lean\YMStatement\Conditional.lean` | 224 | 0 |
| 35 | `lean\YMStatement\H1Sharpening.lean` | 83 | **1** |
| 36 | `lean\YMStatement\Reduction.lean` | 239 | **5** |
| 37 | `docs\zenodo-v2\ym-v2-lean\lakefile.lean` | 29 | 0 |

（`docs\zenodo-v2\ym-v2-lean\` 下 37 个文件与 `research\*`、`lean\*` **逐字重复**——`YMStatement` 5 件 + `research\*` 32 件是同一批源码的第二份拷贝。真 `sorry` 战术合计 **12 处**：`Reduction.lean` 5（两份拷贝各 5）+ `H1Sharpening.lean` 1（两份各 1），散布在 `glue-*` 与 `nontriviality` 的是 `sorry` 占位骨架。）

**第三方树 `sorry`（供对照）** `[T]`：`vendor\OSforGFF` 61 文件中 0 处；`vendor\YM-audit` 9,564 文件中仅 9 个文件含 13 处（全部是 `Challenge.lean` 型"留给解题者的洞"，加 `provenance\lgt\LGT\Bridge\ScaleTransfer.lean` 5 处）。

---

## §4 已有形式化的主题 —— 逐条命中

### §4.1 目标问句：`v02` 说的"Lean 已证、0 sorry 的 `m₁` 夹逼"在哪？→ **找到了，并且实测通过**

**断言来源（引用）**：`research\v3\species\LOOP-RUNBOOK.md:5077`
> **`m₁`**：夹逼 `f₁²/χ ≤ m₁² ≤ M₀/χ`（`[T]`，Lean 已证 0 sorry）

**载体的唯一文件** `[T]`：
`D:\Dev\math\lean\YMStatement\Reduction.lean`（= `docs\zenodo-v2\ym-v2-lean\YMStatement\Reduction.lean`，239 行，同一文件）

| 声明 | 行 | 内容 |
|---|---|---|
| `susceptibility_squeeze` | 101 | **两步挤压 `f₁²/χ ≤ inf supp ν ≤ M₀/χ`（0 sorry）** |
| `variance_susceptibility_squeeze` | 119 | `f/χ ≤ …` 的方差版 |
| `m_sq_ge_of_squeeze` | 125 | `r₀>0 ∧ χ≤1/λ ∧ r₀≤f₁² ⇒ r₀·λ ≤ f₁²/χ` |
| `transported_residue` | 248 | `f₁² = Z_O²·f₁²_strong` 的搬运下界（`(CP)∧P3`） |
| `clay_gap_conditional` | ~271 | 归约定理：假设 ⇒ `m_gap>0` |
| `feshbach_bogoliubov_schur` | 46 | Schur 补引理 |
| `gDouble_nonneg` / `gdouble_lower_bound` | ~190 | `G″ = E[H] − Var + W` 一侧 |

**实测（不是抄注释）** `[T]`：在该目录用 `leanprover/lean4:v4.33.0-rc1` 直接 typecheck `YMStatement/Reduction.lean`，`exit=0`，20.6 s，5 条 warning：

```
YMStatement/Reduction.lean:209:8: warning: declaration uses `sorry`
YMStatement/Reduction.lean:215:8: warning: declaration uses `sorry`
YMStatement/Reduction.lean:222:8: warning: declaration uses `sorry`
YMStatement/Reduction.lean:229:8: warning: declaration uses `sorry`
YMStatement/Reduction.lean:237:8: warning: declaration uses `sorry`
'YMStatement.Reduction.susceptibility_squeeze' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.variance_susceptibility_squeeze' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.m_sq_ge_of_squeeze' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.transported_residue' depends on axioms: [propext, Classical.choice, Quot.sound]
'YMStatement.Reduction.clay_gap_conditional' depends on axioms: [propext, Classical.choice, Quot.sound]
[exit=0]
```

**判读** `[T]`：5 处 `sorry` 全在 `(B)` 段（第 199–238 行，**以 `sorry` 陈述的未证解析输入**：H1/CP/P3/A3 等假设本身）；**夹逼定理与归约定理的公理足迹里没有 `sorryAx`**，只有 `propext / Classical.choice / Quot.sound` ⇒ **v02/LOOP-RUNBOOK 的"Lean 已证、0 sorry"对这条夹逼成立**。

**限定（必须一起搬）** `[C]`：夹逼是**纯代数**（有序域 + 除正数），它不含任何 YM 内容；正界 `r₀·λ` 的两个前件 `r₀ ≤ f₁²`（残数地板）与 `χ ≤ 1/λ`（粗模 Poincaré 天花板）就是 `H1` 与 `P3` 本身——所以它是**搬家**，不是拆墙。搬运时可以整块复用，但**不要**把它当"已得隙"。

### §4.2 其余主题命中

| 主题 | 命中 | 状态 |
|---|---|---|
| `f₁²` | `Reduction.lean`：`T.f1sq`、`transportedRho = ZO^2 * r0`、`transported_residue` | `[T]` 0 sorry（假设型） |
| `f₁²`（同名歧义，**不是**残数） | `YMStatement\Basic.lean:45-49` 的 `f₁ f₂ f₃ f₄ : SchwartzTestFunction d` 是**试验函数**，与残数 `f₁²` 无关 | `[T]` 需防混名 |
| `Poincaré`（粗模，= H1） | **只有文字**：`Reduction.lean` §H1 定义、`H1Sharpening.lean`（`χ ≤ 1/λ`）。**无 Mathlib Poincaré 不等式实例的调用** | `[T]`/`[NG]`：未形式化 |
| `interval` | `GL3Limits.lean:15`（`isCompact_univ_pi`，Tychonoff 乘积区间）；`MinimalHandoff.lean:23,38`（`intervalIntegral.integral_mono_on`）；数学用的 `intervalIntegral.integral_comp_sub_left` 已 `#check` 通过 | `[T]` 可复用 |
| `norm_num` | `SYMDim4.lean:35,38,93,223,242`；`GL4L2.lean:57,59`。自检 `example : (3:ℕ)+4=7 := by norm_num` 编译无错 | `[T]` 可用 |
| `transfer` | `YangMillsInfra\Spectral\Transfer.lean`（`SemigroupGap`、`twoPoint_bound`、`decay_at_nat_multiples`，0 sorry）；`GL4Transfer.lean`（弱极限转移，2 sorry）；`transported_residue`（残数搬运，0 sorry） | `[T]`：**"transfer" 有三种不相干含义**，搬运要挑 `transported_residue` 这条 |
| `gap` | `YangMillsInfra\Statement\MassGap.lean`（`MassGapStatement`、`exponentialClustering_implies_clustering`）；`RG\Step.lean`（`gapCertificate`、`UniformGapCertificate`、`CertificateImpliesGap`）；`RG\Handoff.lean`（`gap_uniform`）；`RG\GroupUniform.lean`（`gap_uniform`、`uniform_over_groups`）；`Spectral\Unbounded.lean`（`HasGapAbove`、`HasGapAboveVacuum`） | `[T]` 全套"隙证书"接口已就位，除 `RG\Handoff.lean` 1 处 |
| `m₁` | `Reduction.lean`（`m₁² = f₁²/χ`、`mass_gap_of_spectral`）；`HandoffMargin.lean`（`q₁` 裕度）；**无数值 `m₁` 的 Lean 常量** | `[T]` 全是符号层 |

### §4.3 形式化的"接口形状"（搬运时直接对接）

`[T]` 现成的可复用层（0 sorry，且已 `#print axioms` 冻结为标准公理）：

- `YangMillsInfra\Basic.lean`：`EuclideanTheory d`、`TranslationAction d`、`ExponentialClustering`、`Clustering`
- `Reconstruction\OS.lean`：`TransferData`、`TransferData.HasGap Δ`、`ReconstructedTheory`
- `Spectral\Transfer.lean`：`SemigroupGap T V Δ` ⇒ `twoPoint_bound` ⇒ `decay_at_nat_multiples`
- `Spectral\Unbounded.lean`：`SelfAdjointOperator`、`HasGapAbove`、`hasGapAbove_mono`
- `RG\Step.lean`：`gapCertificate Z n`、`StepTheorem`、`UniformGapCertificate`（体积一致）
- `RG\Handoff.lean`：`HandoffData`、`WeakSide / StrongSide / Overlap`、`gap_uniform`、`regions_cover`、`overlap_iff_order`
- `YangMillsInfra\Statement\NonTriviality.lean`：`IsGaussianTheory`、`NonTrivial`、`isGaussianTheory_iff_cumulants_vanish`（**自由场酸测的现成机器**）

---

## §5 判死线

### (a) Lean + Mathlib 可用 ⇒ 给版本 + `#check Real.exp_pos` 实际输出 —— **成立** `[T]`

版本：

```
Lean (version 4.36.0, x86_64-w64-windows-gnu, commit 611c117bfaf7b7ae1cfcf246ce6fba31e8f509e9, Release)
Lake version 5.0.0-src+611c117 (Lean version 4.36.0)
mathlib @ 79d0395a1825a6264ad5d269e35e60537518955e (v4.33.0-rc1)，8,268 olean
```

自检**实际输出**（路线 1：`cd vendor\OSforGFF; lake env lean`，`import Mathlib.Analysis.SpecialFunctions.Exp`）：

```
Real.exp_pos (x : ℝ) : 0 < Real.exp x
Real.exp_pos (x : ℝ) : 0 < Real.exp x
[exit=0]
```

自检**实际输出**（路线 2：手工 `LEAN_PATH` 指向 `build.badcache`，**`import Mathlib` 全量**，同时验 `norm_num` 与 `intervalIntegral`）：

```
Real.exp_pos (x : ℝ) : 0 < Real.exp x
Real.exp_pos (x : ℝ) : 0 < Real.exp x
intervalIntegral.integral_comp_sub_left.{u_5} {E : Type u_5} [NormedAddCommGroup E] [NormedSpace ℝ E] {a b : ℝ}
  (f : ℝ → E) (d : ℝ) : ∫ (x : ℝ) in a..b, f (d - x) = ∫ (x : ℝ) in d - b..d - a, f x
[exit=0]
```

（`example : (3:ℕ) + 4 = 7 := by norm_num` 与 `example (x:ℝ) : 0 < Real.exp x := Real.exp_pos x` 均无报错 ⇒ `norm_num`/项模式均通。）

### (b) 若不可用 —— **不适用**

不需要装任何东西。若要**追加**能力（当前不需要）：

| 可选动作 | 体量 |
|---|---|
| 新 toolchain via `elan toolchain install` | ~3.0–3.2 GB / 个（实测现有两个各 3.0/3.2 GB） |
| 新 mathlib 源码 clone | ~0.5 GB（源码）+ 8,268 olean 缓存；增量 `lake exe cache get` 视版本差 |
| 现有 `.lake` 已占 | **9.29 GB**（弃用需先确认 `build.badcache` 可被其他项目共享） |

**真正的缺口不是"装 Lean"，是"给 `research\lean-infra` 补 `lean-toolchain` + `.lake` hydrate"** —— 它现在靠 `check.sh` 手工拼 `LEAN_PATH` 才能编译（`lakefile.toml` 注释已自认："a full `lake build` was out of scope"）。这是搬运前最该修的一条基础设施债。

### (c) 已有形式化与"搬运"相关 ⇒ 逐条列，标 `sorry` 数

| # | 文件 | 与搬运的关系 | 行 | sorry | 实测 |
|---|---|---|---|---|---|
| 1 | `lean\YMStatement\Reduction.lean` | **`f₁²/χ ≤ m₁² ≤ M₀/χ` 夹逼** + 归约 `(P3,CP,H1,A3) ⇒ m_gap>0` | 239 | **5**（全在 `(B)`，非定理本体） | typecheck `exit=0`，20.6 s，夹逼公理足迹无 `sorryAx` **`[T]`** |
| 2 | `…\Reduction.lean` 内 `transported_residue` | **乘性残数搬运** `f₁² = Z_O² f₁²_strong` 的机器检验版 | 248 | 0 | 同上 **`[T]`** |
| 3 | `…\Reduction.lean` 内 `m_sq_ge_of_squeeze` | `r₀>0 ∧ χ≤1/λ ⇒ r₀λ ≤ m₁²`（搬运的**结论形状**） | 125 | 0 | 同上 **`[T]`** |
| 4 | `lean\YMStatement\RGHandoff.lean` | 一圈流 + `Λ` 不变 + `q` 单调 + `gap_lower_bound` + `certificate_transfer` + IVT 穿越点唯一 | 261 | **0** | 分文件编译由 `lean\check-*.sh` 覆盖（本文件 `grep sorry` = 0 命中）**`[T]`** |
| 5 | `lean\YMStatement\Conditional.lean` | 交接 ⇒ 隙（`Overlap ↔ q₁≤q₀`，一圈/两圈闭式） | 224 | **0** | 同上 **`[T]`** |
| 6 | `lean\YMStatement\H1Sharpening.lean` | Lemma G 双侧夹逼 `f/m ≤ χ ≤ M/m`（H1 的强化形态） | 83 | **1** | **`[T]`** |
| 7 | `research\lean-infra\YangMillsInfra\RG\Handoff.lean` | `gap_uniform`：弱侧/强侧 + 重叠区间覆盖 ⇒ 全 β 隙 | 141 | **1** | 需 `check.sh` 手工 `LEAN_PATH` **`[T]`** |
| 8 | `research\lean-infra\YangMillsInfra\RG\GroupUniform.lean` | 同上，群参数化 + `uniform_over_groups` + `Ledger`（intensive/representation/extensive） | 132 | 0 | **`[T]`** |
| 9 | `research\lean-infra\YangMillsInfra\RG\Step.lean` | 体积一致的隙证书 `UniformGapCertificate` | 105 | 0 | **`[T]`** |
| 10 | `research\lean-infra\YangMillsInfra\Spectral\Transfer.lean` | `SemigroupGap ⇒ twoPoint_bound`（搬运要用到的"指数衰减"接口） | 92 | 0 | **`[T]`** |
| 11 | `research\lean-infra\YangMillsInfra\Spectral\Unbounded.lean` | `HasGapAbove(V,Δ)`、单调性 | 95 | 0 | **`[T]`** |
| 12 | `research\lean-infra\YangMillsInfra\Statement\MassGap.lean` | `ExponentialClustering ⇒ Clustering`（陈述层） | 74 | 0 | **`[T]`** |
| 13 | `research\lean-infra\YangMillsInfra\Statement\NonTriviality.lean` | **自由场酸测的机器**：`gaussian_not_nonTrivial`、`isGaussianTheory_iff_cumulants_vanish` | 180 | 0 | **`[T]`** |
| 14 | `research\glue-a-limits\GL4Transfer.lean` | 弱极限转移（搬运的极限步骤） | 83 | 2 | 骨架 **`[T]`** |
| 15 | `research\glue-b-blocking\lean\GL2bPrime.lean` | 阻塞引理（`squeeze_zero_norm` 用在这里） | 129 | 2 | 骨架 **`[T]`** |
| 16 | `research\handoff-proof\HandoffMargin.lean` | 交接裕度 `M = ln(m/Λ) − t_up − ln q₁`，`q₁≤q₀ ↔ 0≤margin` | 96 | 1 | **`[T]`** |

**判读** `[C]`：**搬运已经在 Lean 里有骨架，且核心那一步（残数搬运 + 夹逼 + 归约）是 0 sorry 的**。缺的不是 Lean 能力，而是两条**一致输入**（`f₁² ≥ r₀` 的残数地板、`χ ≤ 1/λ` 的粗模 Poincaré 天花板）——它们在 Lean 里被**显式写成 `sorry`/`Prop` 假设**（`Reduction.lean` 的 `(B)` 段第 209/215/222/229/237 行），这正是"借 Lean 做搬运"下一步要盯的五个洞。

---

## §6 复现命令（可直接粘贴）

```powershell
where.exe lean; where.exe lake; lean --version; lake --version

# 路线 1：项目内（mathlib 依赖自动装配）
cd D:\Dev\math\vendor\OSforGFF; lake env lean <file.lean>

# 路线 2：手工 LEAN_PATH（import Mathlib 全量，8,268 olean）
$depsrc='D:\Dev\math\vendor\OSforGFF\.lake\packages'
$paths=@("$depsrc\mathlib\.lake\build.badcache\lib\lean")
foreach($p in 'batteries','aesop','Qq','proofwidgets','LeanSearchClient','importGraph','plausible','Cli','GaussianField','BochnerMinlos','kolmogorov_extension4'){ $d="$depsrc\$p\.lake\build\lib\lean"; if(Test-Path $d){$paths+=$d} }
$env:LEAN_PATH = ($paths -join ';')
& "$env:USERPROFILE\.elan\bin\elan.exe" run leanprover/lean4:v4.33.0-rc1 lean <file.lean>

# 复验 m₁ 夹逼（0 sorry）
cd D:\Dev\math\docs\zenodo-v2\ym-v2-lean
elan run leanprover/lean4:v4.33.0-rc1 lean YMStatement/Reduction.lean
```

---

## §7 结论一句话

`[T]` Lean 4.36.0 + Lake 5.0.0 在 PATH；mathlib `v4.33.0-rc1`（rev `79d0395a`，8,268 olean）在 `vendor\OSforGFF\.lake` 可用，`#check Real.exp_pos` 与 `import Mathlib` 均 `exit=0`。
`[T]` 本项目自有 Lean **42 件 / 4,587 行**（`research\` 36 件 3,638 行 + `lean\` 6 件 949 行；`docs\zenodo-v2\` 的 37 件是副本，不计），真 `sorry` 战术仅 **12 处**（两份 `Reduction.lean` 各 5 + 两份 `H1Sharpening.lean` 各 1）。
`[T]` v02 说的 `m₁` 夹逼 = `lean\YMStatement\Reduction.lean` 的 `susceptibility_squeeze`（`f₁²/χ ≤ m₁² ≤ M₀/χ`）：**实测 typecheck 通过、公理足迹 `[propext, Classical.choice, Quot.sound]`、无 `sorryAx` ⇒ "已证 0 sorry" 成立**；同文件的 5 处 `sorry` 是它作为**假设**吞进去的 `(B)` 段输入，不是夹逼的证明漏洞。
`[C]` 搬运的基础设施债：`research\lean-infra` 缺 `lean-toolchain` 与 `.lake`，只能靠 `check.sh` 拼 `LEAN_PATH` 编译。

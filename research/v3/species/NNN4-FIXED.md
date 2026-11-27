# NNN4-FIXED — 「固定 `(L,a)`：`Δ_a>0`」闭成 `[T]`（含紧性修正）

> 第 NNN4 轮 · 只新建本文件；既有文件零改动；未 commit；零新 MC；未编译 Lean（无 `lake build`、无 `#check` 运行）。
> 口径：紧群 `G`、Wilson 作用、有限 `L`、`0≤β<∞`、固定 `a>0`；`T_a=e^{−aH}`、`Ω≡1`、`λ₀=1`、`Δ_a=(1/a)ln(λ₀/λ₁)`。标签只用 `[T]`/`[NG]`。

## §0 判词

| 问 | 判 | 下一手 |
|---|---|---|
| 三行整体 | 逐条 `[T]`，需三处修正：真核是跨时间片积分、`f` 复化要拆、缺第四条输入 `0≤T_a≤1` | 见 §1 |
| 固定 `(L,a)` 的 `Δ_a>0` | **无条件 `[T]`**（窗口见上），零新输入 | 记为「地板免费」，靶仍是 `∃Δ₀>0∀(L,a)` |
| 同一命题的 Lean 形 | 仓内现状 `[NG]`：`T_a` 未定义、核未定义、无 HS | §2 的 A/B 两级 |
| `TransferData` | 不动；新增紧 + 正的新结构 | §3 裁决 3 |
| 紧性修正 | FFF2 的修正成立；再补一条：缺的不只是紧，还有**谱非负** | §1 第 5 行 |

## §1 三行逐条

| # | 该条 | 出处（`file:line`） | 判 | 下一手 |
|---|---|---|---|---|
| 1 | `K` 有界 | Lean：`wilsonWeight = ofReal(exp(−S_β))`（`research/lean-infra/YangMillsInfra/Gauge/Measure.lean:66–68`）；`S_β=βΣ_p(1−χ(U_p))`（`…/Gauge/Lattice.lean:130–131`）。`|χ|≤1` ⇒ `|S_β|≤2β·#Plaq` ⇒ `e^{−C}≤K≤e^{C}`。文献：`D02-DIRECT.md:59`、`TT4-WILSON.md:23` | `[T]` | 给 `CharacterData` 补 `χ_le_one`/`χ_ge_neg_one`（现只有 `χ_one`、`χ_conj`，`Lattice.lean:120–128`） |
| 2 | `K>0` 处处 | `Real.exp_pos` 对任意实 `S_β`；真核 `K_a(U,U')=∫dU₀(跨时间片项)`，被积函数严格正且连续 ⇒ 正测度集上积分为正（`os3-reflection-positivity/STATUS.md:20`、`HARNACK-KILLTEST.md:22`：`e^{−S}=∏_t T(U_t,U_{t+1})`） | `[T]` | Lean 需「严格正连续函数积分 `>0`」一条（`MeasureTheory` 有零件，未拼装） |
| 3 | `T_a1=1` | 归一化 Wilson 测度全测度 1（`Gauge/Measure.lean:79–90` 的 `wilsonMeasure_univ`）；`0≤T_a≤1`、`T_a1=1`（`research/glue-b-blocking/REPORT.md:267`） | `[T]` | 这条在 Lean 里是结构缺口：无时间片分解，`T_a` 本身不存在 |
| 4 | 紧性（HS） | `K` 有界 + `μ` 概率 ⇒ `‖K‖_{L²(μ⊗μ)}≤‖K‖_∞` ⇒ `T_a` 是 Hilbert–Schmidt ⇒ 紧（`D02-DIRECT.md:59`；`FFF2-FREE.md:15`）；作用空间无穷维 `L²(SU(2)^{3L³})`（`D02-DIRECT.md:41`） | `[T]`（数学）；Lean `[NG]`（无 HS，§2） | 二选一：自建 HS（Hilbert 基 + 截断）或 Arzelà–Ascoli 走 `C⁰` |
| 5 | **第四条输入 `0≤T_a≤1`** | `T_a=BB*`（slab 因子的 Gram 形）⇒ `⟨f,T_af⟩≥0`（RP；`research/gluing/REPORT.md:350`「`T≥0`」、`FFF2-FREE.md:17` 的 `2′`）；`‖T_a‖≤1` 同处 | `[T]` | 无此条时 `λ₁` 可为负、`Δ_a` 无定义 ⇒ 必须写进 §3 的结构 |
| 6 | 极值原理 ⇒ 顶端单重 | `K(U₀,·)>0` 且 `∫K(U₀,·)dμ=1`：`M=f(U₀)=∫Kf≤M`，等号 ⇒ `f≡M` ⇒ `M=∫f=0` ⇒ `f=0`（`FFF4-LAM1.md:39`） | `[T]`；修正：`f` 复 ⇒ 实/虚部分别做；`M≤0` 由 `∫f=0` 先排除 | 把这两步补进证明文字 |
| 7 | 收尾 `λ₁<1` | 紧自伴 ⇒ 1 孤立、次一点取到；更强形 `q:=‖T_a|_{1^⊥}‖<1` ⇒ `Δ_a≥−(1/a)ln q>0`（`β=0`：`q=0`，`Δ=+∞`） | `[T]` | Lean 友好陈述用 `q<1`，不用「`λ₁`」 |

**修正汇总**：`K` 不是裸 `e^{−βS_W}`（含空间片因子与时间链积分），但 `[T]`；`f` 复化必须拆；`0≤T_a≤1` 是独立第四条，不能由「紧 + `K>0` + 归一」推出（`FFF4-LAM1.md:24` 的 5×5 见证只覆盖前三条）。

## §2 Lean 化

**有**（本地 `vendor/OSforGFF/.lake/packages/mathlib`，8268 个 `.lean`）：`IsCompactOperator`（`Analysis/Normed/Operator/Compact/Basic.lean:71`）；Fredholm 二择一与非零谱 = 本征值（`Compact/FredholmAlternative.lean:164`、`:221`）；紧自伴谱定理（`Analysis/InnerProductSpace/Spectrum.lean:443` 本征空间全体稠、`:463` 本征空间有限维）；Rayleigh 极值即本征向量（`Rayleigh.lean:276`）与 `spectralRadius_eq_nnnorm`（`:182`）；Loewner 星序 `IsPositive`（`Positive.lean:59`、`:266`）与 `IsPositive.spectrumRestricts`（`StarOrder.lean:34`：`0≤T` ⇒ 实谱 `≥0`）；Arzelà–Ascoli（`Topology/UniformSpace/Ascoli.lean`）。
**无**：`HilbertSchmidt` 全库 0 命中；`MeasureTheory` 内 `integralOperator|IsCompactOperator` 0 命中（无「积分算子紧」定理）；开库只见到 PR #41477（限有限维 HS 内积）。

**引理清单**（`A` = 条件式，`B` = 模型级）：
1. `A/χ`：`|χ|≤1` ⇒ `0<K≤e^{2βN_p}`；工具 `IsCompact.exists_bound_of_continuousOn`（仓内先例 `GW1-bridge.lean:275`）。~60 行。
2. `B/T`：定义 `T_a`（`L²(links)` 上的核 + 时间片分解 + `T_a1=1`）；现 `Lattice.lean:38–42` 只有 `V/L/Plaq`，无时间层。~800–1500 行。
3. `B/紧`：`K∈L²(μ⊗μ)` ⇒ `IsCompactOperator`；无 HS API ⇒ 自建（~500–1200 行）或 Ascoli（~300–800 行）。
4. `B/连续+极值`：`T_af∈C⁰`（参数积分连续，需拼 `MeasureTheory.Integral.*`）~200–400 行；极值原理与「等号 ⇒ 常值」~200–400 行。
5. `A/收尾`：`1^⊥` 限制（`IsCompactOperator.restrict'`，`Spectrum.lean:448` 已用）+ `spectralRadius_eq_nnnorm` ⇒ `q<1`。~100–250 行。
6. `A/谱非负`：`0≤T≤1` ⇒ 实谱 ⊂ `[0,1]`（`spectrumRestricts`）。~50 行。

**量级**：A 级（`CompactTransfer` 结构 + 第 5/6 条）~300–600 行 / 3–5 轮；B 级（第 2–4 条，真正无条件）~2500–5000 行 / 10–20 轮。

## §3 裁决

1. **数学裁决**：固定 `(L,a)`、Wilson、紧 `G`、`0≤β<∞` ⇒ `Δ_a>0` 无条件 `[T]`。四条输入不可去：去 `K>0`、去紧、去 `T1=1` 各给 `λ₁=1`（`FFF4-LAM1.md:24`）；去 `0≤T_a` 则 `λ₁` 可负。下一手：本条封为「地板免费」，不再投人。
2. **判别力裁决**：`[T]` 但零判别力——自由 `U(1)`、2D YM 同真（酸测，`GH1-SPECTRAL.md:111`）。靶行只留 `∃Δ₀>0∀(L,a)`。下一手：归 `H1` 栏。
3. **定义裁决**：`TransferData`（`Reconstruction/OS.lean:98–116`）**不动**；新增 `CompactTransfer`，字段 `compact : ∀t>0, IsCompactOperator (transfer t)`、`quadPos : 0≤transfer t`、`kernel_pos`、`row_normalized`。只加紧不够，谱非负不可省。下一手：先做 A 级把 `q<1` 证出来。
4. **Lean 裁决**：现状 `[NG]`（`T_a` 不存在）；A 级后 `[T]`（条件式）；B 级后才无条件。未编译 Lean，故这是接口级判断，不是编译通过记录。下一手：B 级第一步只做 `Lattice` 的时间分层。

**工具状态（AGENTS §1 四要素）**：`pwsh` OK（列目录、取行、文件计数）；`grep` OK——式 `HilbertSchmidt` @ Mathlib 8268 文件 → 命中 0；式 `IsCompactOperator` @ Mathlib → 132；式 `IsCompactOperator` @ 仓内非 vendor `*.lean` 94 文件 → 6（全部位于 `out/mathlib-probe/`）；式 `HilbertSchmidt|Hilbert–Schmidt` @ 仓内非 vendor 2094 文件 → 28（均为文档引述）；`read` OK（`OS.lean`、`Transfer.lean`、`Spectrum.lean:380–470`、`Lattice.lean:100–150`、`Measure.lean:55–95`、`D02-DIRECT.md:40–69`）；`web_search` OK（3 式，返回 mathlib4 `Compact/Basic.lean` 与 PR #41477）。无 `TOOL-FAIL`。未查到（检索式如上，工具 OK，命中 0）：Mathlib 的 Hilbert–Schmidt API、积分算子紧性定理。

# 16 输入审计：交接不等式 `q₀ ≥ q₁` 是否已在其中？

> 2026-11-22 · 来源：Gross 委员会评估 `vendor/YM-audit/docs/assessment-1.pdf`、OpenAI 提交的 Lean 源、hexagonmath 重构 `hexagon:2611.00492v1`。

## 0. 结论（先看这三条）

1. **`q₀ ≥ q₁` 不在 16 个输入里。** 16 个输入全部是 UV / OS / 结构 / 有限格点输入；没有一个给出强耦合控制、IR 唯一性或质量隙。
2. **因此交接不等式是一个新输入，不是隐藏的已有输入。** `handoff-proof` 会话的任务是"补一个新估计"，不是"证明一个已经写好的输入"。
3. **另一个独立风险：** 如果交接接口把 `cert` 塌缩成单个常数（gluing §2.7 的示意签名如此），则 `H3` 变成不参与结论的空假设，整个"两域覆盖"论证消失。正确签名必须保持 `cert : ℝ → ℝ`。

## 1. 来源与引用更正

- **权威来源**是委员会评估（10 页正文 + 检查记录）和 `submissions/openai-20261028/submission/` 的 Lean 源。16 个标签来自 hexagonmath 的**独立重构**，不是委员会评估使用的标签；委员会的 A1–A8 是检查项，不是假设项。
- **`hexagon:2611.00492v1` 不是 arXiv 编号。** 它是 hexagonmath 站内 ID（作者 Chunping Zeng，站点 `https://hexagonmath.org/pdf/2611.00492v1`）。`arXiv:2611.00492` 实际是一篇语言学论文（Dialogic Humor）。引用时必须写 `hexagon:2611.00492v1`，不能写 arXiv。

## 2. 16 个外部输入

按 hexagonmath 的列表（PH4 由 AN1 E2+E4 推出，不计入 16；RG4 由 RG2+RG3 推出，但仍在"已认证包"里保留编号）。

| # | 标签 | 一句话内容 |
|---|---|---|
| 1 | **RG1** | 保对称的规范固定：规范固定对破坏对称的选择取平均，blocking kernel 在格点对称群下等变 |
| 2 | **RG2** | 单步估计，对体积一致：存在 `κ ∈ (0,1)`、`K>0`、`η>0`，`2b₀ log 3 − K/η > 0`，且 block 耦合 `b_k ≥ 3` |
| 3 | **RG3** | 插入读出稳定性：reference level `r(d)` 的归一化系数等于 `v`；从 level `N` 回到 `r(d)` 的 RG 传输线性映射可逆 |
| 4 | **RG4** | 差分估计（由 RG2+RG3 推出）：`|C_{N+1,p}(v;f) − C_{N,p}(v;f)| ≤ A(v,f) κ^N` |
| 5 | **RG5** | 一致光滑性：`sup_{N,p} sup_{x∈K} |∂^α C_{N,p}(v;x)| ≤ M_{α,K,v} < ∞` |
| 6 | **UV1** | 局部有界 + 等连续：对每个标签表、紧集和多指标，存在 `p₀`、`M`、模数 `ω` 使 `|S_p(v;x)|` 一致有界且等连续 |
| 7 | **PH1** | 时间平移共轭：`B_{s,A} = T_s B_{0,A} T_{−s}`，且 `P₀T_s = P₀ = T_sP₀` |
| 8 | **PH2** | 交换性：`A` 交换 |
| 9 | **PH3** | 忠实性：`aΩ = 0 ⟹ a = 0`（`Ω` 是分离向量） |
| 10 | **AN1** | 所选相的关联满足 OS 公理 E0（温性）、E0′（arity 增长）、E1（欧氏不变）、E2（反射正性）、E3（置换）、E4（聚类）、归一化与解析性 |
| 11 | **Y1** | 应力张量唯一性：满足对称、中心化、张量变换、散度为零的元素空间一维 |
| 12 | **Y2** | 两条定律：通量（shell）律与 Ward 恒等式，对每个标签表和不同点成立 |
| 13 | **Y3** | 短距输入：`r→0` 时 `K(r)` 的两点关联展开，含 running coupling `g_r = g(1/r)` 与 `β(g)` |
| 14 | **Y4** | 一圈值：`b₀ = 11 C_g/(48π²)`（规范玻色子 + ghost + measure 的贡献之和） |
| 15 | **B1** | 场内容字典：`F_d ≅ Maxwell.shellSector(G,d)`，且对每个相成立 |
| 16 | **C1** | 与独立格点形式化的桥：有限 Wilson 律在 `β_submission = m β_lgt` 下逐点一致 |

**这些输入覆盖的是**：UV 极限、子列无限体积极限、相分解与聚类、OS 重构、应力张量与 Ward 恒等式、迹反常、一圈 running、OPE、有限格点桥。
**它们不覆盖**：强耦合区的控制、弱→强交接、质量隙的速率、无限体积唯一性、禁闭。

## 3. 审计问题：有没有交接 / 质量隙？

**没有。** 证据链：

1. **委员会评估的检查范围表**（`submissions/openai-20261028/docs/checks.md`）明确写：
   > "The combined theorem does not assert uniqueness of the large-volume limit, a mass gap, confinement, or all-orders perturbative matching."

2. **A5 / `TwoPointScale.MassLimit`** 只给出"正的格点探针"，并且 `checks/TwoPointPhysicalObservable.lean` 的注释写明：
   > "no intrinsic scale or spectral gap is asserted."
   `checks/NonzeroPhysicalObservable.lean` 同样写："No intrinsic mass, mass gap, or IR uniqueness is claimed."

3. **评估 §7.4** 只证明短距离 `g_W(r) > 0` 与 `log(1/r) g_W(r)² → 1/(2b₀)`；并明确写：
   > "The submission does not prove `g_W(r)` is independent of the subsequential infinite-volume limit."
   没有任何关于 `r → ∞`（IR）或强耦合的陈述。

4. **评估 §8 Open Questions** 把"uniform mass gap"和"single-particle states"列为未解决。

5. **AN1 的 E4 是定性聚类**，不带速率；它不蕴含指数衰减，也不蕴含谱隙。PH4（真空一维）由 AN1 E2+E4 推出，仍不含速率。

结论：`q₀ ≥ q₁` 是一个**新的、外部的**输入。把它写成条件定理时，必须显式声明它是新增假设（或新增定理），不能声称它已在 16 个输入中。

## 4. 接口空洞性风险（必须修）

`physics_handoff.md` §3 已给出警示，`gluing/REPORT.md` §2.7 的示意签名仍把它塌缩成常数。两种写法的差别是决定性的：

**错误（H3 空洞）：**
```lean
-- gapBar 是单个常数
theorem handoff (gapBar : ℝ) (hW : δW ≤ gapBar) (hS : δS ≤ gapBar) :
    0 < gapBar → ...
```
这里 `H1`/`H2` 已经直接给出 `gapBar` 的下界，`H3` 不参与结论；再加 `0 < gapBar` 就是同义反复。

**正确（`cert` 依赖 μ）：**
```lean
theorem handoff (cert : ℝ → ℝ) (q g2 : ℝ → ℝ)
    (hW : ∀ μ > 0, g2 μ ≤ gstar^2 → δW ≤ cert μ)
    (hS : ∀ μ > 0, q1 ≤ q μ → δS ≤ cert μ)
    (hOv : q1 ≤ q0) :          -- 等价地 q1 ≤ q μup；因为 q μup = q0
    ∀ μ > 0, min δW δS ≤ cert μ
```
再用证书语义 `∀ μ, cert μ ≤ Δ/Λ` 收尾，得到无量纲的 `Δ/Λ ≥ min(δW,δS) > 0`。

**行动**：gluing §2.7 的示意签名必须改成 `cert : ℝ → ℝ`；否则九变一的论证会在组合时静默失效。

## 5. 对 `handoff-proof` 会话的直接含义

- **不要花时间在 16 个输入里找交接项**——已经查过，没有。
- 任务应表述为：**新增一个输入（或证明它）**，即 `q₀ ≥ q₁`，并保持 `cert : ℝ → ℝ`。
- 如果找不到证明，正确的交付是：把 `q₀ ≥ q₁` 写成一个带证伪测试的**新猜想**，而不是把它伪装成已有输入。
- 审计同时确认：2026 构造给出的短距 `g_W(r)` 非零与一圈 running 是**弱侧**证据，不覆盖强侧；这正是交接缺口所在。

## 6. 审计的边界

- 本文未重建 Lean 源，只读了评估、检查范围表、验收规格和 hexagonmath 重构的假设陈述。
- 16 个标签的"外部 vs 派生"记账有细微不一致：PH4 由 AN1 推出（外部数 17→16），RG4 由 RG2+RG3 推出（但仍在包内保留编号）。这不影响主结论：**没有一个输入是交接不等式**。
- hexagonmath 是 AI 参与度较高的重构，引用时应以委员会评估和 Lean 源为准；本文用它只为**定位假设标签**，不把它当作独立验证。

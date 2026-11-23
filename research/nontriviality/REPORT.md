# 4D 纯 Yang–Mills 连续极限的非平凡性（非自由场）

> 2026-11-22 · 范围：紧单群 `G`（主例 `SU(N)`）纯 Yang–Mills₄。
> 目标：把 Clay 陈述中的"非平凡（non-trivial）"落实成一个**非空洞、可判真假、可形式化**的命题，并给出定理、引理链、障碍与下一步。
> 边界：本文不攻 MG（质量隙）；与主 session 的弱→强 RG 交接（`research/physics_handoff.md`）正交互补，见 §5.3、§8。

---

## 0. 结论摘要（≤300 字）

"非平凡"的严格含义是**非 Gaussian（非准自由）**。对**纯规范**理论有一个致命的陷阱：任何**规范不变**的局域算子是场强的二次及以上组合（`tr F²` 等），而自由规范场的这类组合**本身就不 Gaussian**（`tr F²` 是 Gaussian 变量的 Wick 平方，连通 4 点不为零）。因此"某个规范不变连通 4 点 ≠ 0"是**空洞**见证——自由理论也通过。唯一同时满足"'规范不变'且'自由值 = 0'"的见证量是**重整化耦合 `g_ren(μ)`**（等价地量纲 transmutation 标度 `Λ > 0`）。

结论：**4D YM 的非平凡性 ≡ 连续极限对应的是非平凡的渐近自由轨迹（`Λ > 0`），而非 Gaussian 平凡轨迹（`g₀≡0`）**。这一步是**定义/审计层**（2026 构造是否真的给出 `Λ>0`，还是只给了 Gaussian 定点处的结构性 RG/anomaly 方程），不是独立的 UV 定理。一旦 `Λ>0`，"短程归一化顶角比 ≠ 0"是 UV 侧可控的连续性陈述（2026 的 one-loop coupling 输入应可覆盖），**不需要** MG 的 IR 强耦合控制。φ⁴₄ 平凡性、OPE、连通 4 点都**不能**单独见证 — 原因各异，见 §2–§4。

---

## 1. 问题：把"非平凡"钉死

> **测度来源（A0 引用）**：非平凡性针对**共享格点测度** `μ_{L,a,β} ∝ exp((β/N_c)Σ_p Re tr U_p)∏_e dU_e`，`β=2N_c/g²`、`N_c:=dim f`（逐字引用 `research/gluing/REPORT.md` §1.2 与 `research/glue-c-gauge-group/REPORT.md` (1.1)–(1.2)；显式式与 β/γ 陷阱见 `PROOF.md §1.1`）。九个接口必须引用同一个 `μ_{L,a,β(a)}`。

### 1.1 正确的形式化对象：非 Gaussian / 非准自由

Clay 原文要求"a **non-trivial** quantum Yang–Mills theory"，并在别处把它解释为"不是自由场"。在 OS / Wightman 框架里，自由场的精确对应物是**准自由态（quasi-free state）**，Euclidean 侧即 **Gaussian 测度**：其生成泛函为
```
G(J) = exp( −½ B(J,J) ) ,   B 双线性。
```
这与 `lean/YMStatement/Basic.lean` 的 `IsGaussianMeasure` 完全一致：
```lean
def IsGaussianMeasure (μ : ProbabilityMeasure (FieldConfiguration d)) : Prop :=
  ∃ B : SchwartzTestFunction d →ₗ[ℝ] SchwartzTestFunction d →ₗ[ℝ] ℝ,
    ∀ J, GJGeneratingFunctional μ J = Complex.exp (-(1/2:ℂ) * ((B J J : ℝ):ℂ))
def NonTrivial μ : Prop := ¬ IsGaussianMeasure μ
```
**这是对的抽象层。** 问题在于：把它接到**纯规范**理论上时，见证量选取极易出错（§1.2）。

### 1.2 关键陷阱：规范不变复合算子对自由场也不 Gaussian

设 `A_μ^a` 是规范联络。自由（Gaussian）规范场意味着 `A` 是 Gaussian 的，故场强 `F = dA`（线性）也 Gaussian。但**任何规范不变局域算子都要色 trace，因而至少是 `F` 的二次式**：
```
O_G(x) := tr F_{μν}(x) F^{μν}(x)   （dim 4，"胶子凝聚算子"，简记 O_G）
```
对 Gaussian 的 `F`，`O_G = :F²:` 是 **Wick 平方**，其连通关联**不为零**：
```
⟨:F²:(x) :F²:(0)⟩_c = 2 · (F-传播子)² ≠ 0,
⟨:F²:(x) :F²:(y) :F²:(z)⟩_c ≠ 0,  ...
```
**因此**：

> **命题 1.1（空洞性）。** 若把"非平凡"定义为"存在规范不变连通 4 点函数 ≠ 0"，则该命题被**自由规范场满足**，即见证是空洞的。

这与 `formal-conjectures` 把 OS 字段写成 `True` 是同一类错误（见 `docs/yang-mills-prior-work.md` §4.5）：陈述看似有内容，实则自由理论即可满足。

连接 `A` 本身的连通 4 点 `⟨A A A A⟩_c` **确实**在自由时为零、相互作用时不为零（有三/四点顶角），但 `A` 不是物理可观测（在物理 Hilbert 空间上无局域规范不变的原场）。所以：

> **命题 1.2（选择约束）。** 合格的见证量 `W` 必须同时满足
> (i) **规范不变**（OS 可观测），(ii) **自由值 = 0**（排除 §1.2 的 Wick 平方）。
> 规范不变的复合算子违反 (ii)；联络违反 (i)。二者兼顾的唯一自然对象是**重整化耦合 `g_ren`**（§3.3、§6.3）。

### 1.3 为什么"连通 4 点"在纯标量理论里却可用

对标量 `φ⁴`，`φ` 本身是规范不变的原场；自由 ⟹ Gaussian `φ` ⟹ 由 Wick 定理 `S₄^c(φφφφ)=0`。故 `Basic.lean` 的 `ConnectedFourPoint` 对**标量**理论是合法见证。这正是现有 Lean 层用标量 `φ` 的原因——但**不能直接照搬到纯规范**，必须换成 `g_ren`（或固定规范下的 `A`，接受非物理性）。这是对 `Basic.lean` 的第一条修改建议（§9）。

---

## 2. φ⁴₄ 的平凡性：survey

### 2.1 定理与文献

- **Aizenman (1981/82)**：`d>4` 上临界维数、平均场行为。
- **Fröhlich (1982)**：由红外界 / 关联不等式得到平凡性框架。
- **Aizenman–Duminil-Copin, *Marginal triviality of the scaling limits of critical 4D Ising and λφ⁴ models*, Ann. Math. 194 (2021), arXiv:1912.07973**：`d=4` 临界标度极限是 **Gaussian**。这是里程碑：`φ⁴₄` 的连续极限**只能是自由场**。
- 早先严格 RG：Gawȩdzki–Kupiainen (1985)、Feldman–Magnen–Rivasseau–Sénéor、Brydges–Mitter–Scoppola；现代：Bauerschmidt–Brydges–Slade (LNM 2242, 2019)。

### 2.2 机制：四次耦合 irrelevant，Gaussian 不动点 IR-吸引

`d=4` 时四次耦合 `λ` 在 Gaussian 定点是**边缘 irrelevant**：
```
β_λ(λ) = (3/16π²) λ² + O(λ³) > 0 ,    λ > 0。
```
故 `λ` 沿 IR 方向**增大**、UV 方向**减小**——Gaussian 不动点是 **IR-稳定**。取连续极限（`UV 截断 → ∞`，固定重整化量）时，一致性强迫裸耦合以足够快的速度趋于 0，使**重整化耦合 → 0**。结论：极限是自由场。平凡性 = **Gaussian 不动点 IR-吸引**。

### 2.3 为何不适用于 YM（决定性对比）

| | `φ⁴₄` | 纯 `YM₄` |
|---|---|---|
| `β` 符号 | `β_λ > 0` | `β_g = −b₀ g³ + … < 0`（渐近自由） |
| Gaussian FP | **IR-吸引** | **UV-吸引** |
| 连续极限下 `g_ren` | `→ 0`（被迫） | **非零**（量纲 transmutation） |
| 结论 | 平凡（自由） | **非平凡**：UV 侧自由、物理（IR）尺度相互作用 |

**平凡性论证的引擎是"Gaussian FP 是 IR 端的吸引子"；YM 恰好相反，Gaussian FP 在 UV 端。** 因此平凡性定理**不适用**，YM 逃逸的机制是**相互作用结构在 IR 端被生成**（`Λ_QCD`）。这也预告了 §1.2 与 §6 的核心困难：非平凡性**最终是 IR 陈述**，只是可以用 UV 侧的连续性去**承接**它。

---

## 3. 渐近自由与量纲 transmutation

### 3.1 β 函数与 Balaban

`SU(N)` 一圈：`β(g) = −b₀ g³ + O(g⁵)`，`b₀ = 11N/(48π²) > 0`（与 `research/physics_handoff.md` 约定一致）。前两系数 scheme-无关（Gross–Wilczek、Politzer 1973）。**非微扰**渐近自由 = Balaban 的 UV 稳定性（CMP 1980s），也是 2026 构造 UV 侧的主体。

### 3.2 量纲 transmutation

```
1/g²(μ) = 1/g²(μ₀) + 2b₀ ln(μ/μ₀) ,     g²(μ) = 1/(2b₀ ln(μ/Λ)) ,
Λ = μ₀ e^{−1/(2b₀ g²(μ₀))} .
```
`Λ` 对任何 `g²(μ₀)>0` 都 **> 0**；`Λ` 在单一 scheme 内 RG-不变，跨 scheme 只乘性 rescale（`Λ' = Λ e^{c/(2b₀)}`）。

> **关键**：`Λ>0 ⟺ g²(μ)≠0 ∀ μ>Λ`。这正是"非平凡轨迹"的定量内容——它排除了 `g₀≡0` 的**平凡（Gaussian）轨迹**（那条轨迹 `Λ=0`，极限是自由场）。

### 3.3 "`g_ren≠0`"是 scheme-无关的（对照 `g²=1` 不是）

scheme 变换 `g'² = g² + c g⁴ + O(g⁶)`。**`g=0` 是 scheme 变换的不动点**：`g=0 ⟹ g'=0`。因此谓词
```
"∃ μ, g_ren(μ) ≠ 0"    在 scheme 变换下不变。
```
而主 session 已数值证明 `g²=1` 的穿越点**随 `c` 平移**（`numerics/t2_scheme_shift.txt`：`c=+0.1` 时 `shift=+0.62`，`c=−0.1` 时 `−0.76`）。所以：
> **用 `g_ren≠0` 作非平凡见证是良定义的；用 `g²=1` 作判据则不是。** 这与 handoff §2 的 scheme 分析一致，并强化了"非平凡性 ≠ 交接点位置"。

---

## 4. OPE：survey

### 4.1 YM 的 OPE 词汇

规范不变局域算子：`1`、`O_G = tr F²`（dim 4）、`tr F F̃`（dim 4，CP-odd）、`T_{μν}`（dim 4）、`tr F³`、`tr(DF)²`、`tr F⁴`、……无限多。短距展开 `O_i(x) O_j(0) ~ Σ_k C_{ijk}(x) O_k(0)`，`C ~ |x|^{Δ_k−Δ_i−Δ_j}`（带对数）。

### 4.2 为何 OPE **系数本身**不能见证非平凡

`O_G O_G → O_G` 的 OPE 系数在**自由**理论里就非零（三个 `:F²:` 的自由连通 3 点非零，Wick 收缩给出）。故"某 OPE 系数 ≠ 0"是 §1.2 空洞性的翻版：**规范不变的 OPE 系数自由理论也不为零**。OPE 是 **UV 工具**（组织短距展开、给出 anomalous dimension 与 running coupling 的形式），**不是**自由/相互作用的判别器。

### 4.3 能见证的是什么

真正的判别来自 OPE 里的**非微扰矩阵元**：
```
⟨O_G(x) O_G(0)⟩ = (identity 系数，Wick/微扰部分) + C_{GGG}(x)·⟨O_G⟩ + …
```
自由 ⟹ `⟨O_G⟩ = 0`（无凝聚）；相互作用 YM ⟹ **胶子凝聚 `⟨(g²/4π²) tr F²⟩ ≈ (0.012 GeV⁴)`（scheme-依赖但非零）**。所以：
> "OPE + 非零凝聚 `⟨O_G⟩ ≠ 0`" 可见证非平凡——但**凝聚是 IR 量**，回到 §2.3 / §6 的结论。

---

## 5. 2026 存在性证明的非平凡性主张

> 来源：`docs/yang-mills-prior-work.md` §3.3–3.5、Gross 委员会报告、hexagonmath 2611.00492。本文**未接触原始论文**（在知识截止之后），以下均为"据该项目调研"的转述并标注**需审计**。

### 5.1 声称了什么

委员会认可清单包含：**非平凡性、守恒应力张量、OPE 界、一圈 running coupling、迹反常**。hexagonmath 的 16 输入重构也复现了"迹反常、一圈 running coupling、OPE"。同时 hexagonmath **证明 16 输入推不出质量隙**（给出满足全部 16 输入、有聚类但无谱隙的 OS 模型），并把质量隙归约为 MG。

### 5.2 审计问题（核心）

把 §1–§4 汇合，可精确提出**三个审计问题**，决定"非平凡性"是否真的成立：

1. **结构性 vs 定标性**：2026 证明给出的是**迹反常方程的结构**（`T^μ_μ = (β/2g) tr F²` 作为算子等式），还是**系数非零的定标**（`β(g)≠0` **且** `⟨tr F²⟩≠0` / `g_ren≠0` 在极限中成立）？前者对自由场也形式成立（`β=0` 时右边为 0），只有后者才是非平凡见证。
2. **轨迹非平凡**：连续极限是沿 `Λ>0` 的 AF 轨迹，还是包括了允许塌缩到 `g₀≡0` 的平凡（Gaussian）轨迹？若是后者，**自由理论即满足陈述**——陈述空洞。
3. **子列塌缩**：只有**子列**极限（委员会明确未认可唯一性 / 收敛性）。不同子列可能对应不同 `Λ`；是否存在某子列塌缩为自由场？若未排除，非平凡性须逐子列声明。

### 5.3 判断（与主 session 正交）

- 若 (1)(2) 成立 → 非平凡性是**非平凡 AF 轨迹的 UV/定义层后果**，其承接（短程归一化顶角比 ≠ 0 传到极限）是 **UV 侧连续性**，2026 的 one-loop coupling 输入**应**能覆盖。**不需要** MG。
- 若只是结构性方程 → 非平凡性继承 MG 的 IR 障碍，与 handoff 同类。

> **这正是本文最主要的新结论：非平凡性是否是"易问题"取决于 (1)(2)，而 (1)(2) 是**可审计**的（读 Lean 陈述签名 / 输入列表），这正是 `docs` 里 M0 审计的直接子任务。** 我倾向 (2) 最可能被"物理重整化条件固定物理标度"隐含地满足，但必须核实。

---

## 6. 精确定理与引理链

### 6.1 核心可证引理 C（Gaussian ⟺ 高阶累积为零）— Lean 现在可做

> **引理 C.** 对满足 OS（尤其 OS1 指数正则 + OS0 解析）的测度 `μ`，
> ```
> IsGaussianMeasure μ  ↔  ∀ n ≥ 3, ∀ f₁…f_n, κₙ(μ; f₁,…,f_n) = 0,
> ```
> 其中 `κₙ` 是连通（累积）Schwinger 函数。若再有**偶性**（`φ ↦ −φ` 不变，标量 `φ⁴` 型），奇数 `κₙ` 自动为零，但**每个偶阶 `4,6,8,…` 都仍需为零**。

⚠️ **修正（机器检验驱动）**：本报告早先写"偶性下 `n=4` 一项即可" —— **这是错的**。反例：`{-1,0,1}` 上权重 `(1/6,2/3,1/6)` 的**偶**测度，`κ₄=0` 但 `κ₆=−2/9≠0`（非 Gaussian）。`PROOF.md §2.3`、`NonTrivialCore.lean` 的 `four_point_insufficient`（无公理）已证。**结论：合格见证须控制全部 `n≥3` 阶累积，连通 4 点单独不足。**

证明骨架：`(⇒)` 由 Gaussian 的 log-MGF 二次 ⟹ 高阶累积为零；`(⇐)` 由 OS1 的正则性保证测度被其指数矩唯一确定，累积**全**零 ⟹ log-MGF 为次数 ≤2 多项式 ⟹ Gaussian（Minlos/Bochner 唯一性，OSforGFF 已实现该唯一性）。**这是把 `NonTrivial` 转化为具体累积见证的全部代数内容**。（完整证明见 `PROOF.md §2`。）

### 6.2 定理 N（非平凡性，coupling form，条件性）

> **定理 N.** 设 `G` 紧单，`{μ_a}` 为 Wilson 测度，`a↓0` 沿 AF 轨迹，重整化耦合在某**物理 scheme** 下固定；`μ` 为 OS **子列**极限。假设
> - **(N1)［UV 稳定性 / AF，2026 输入］** 定义 `g_ren(p)` 的规范不变**顶角/传播子比**在紧动量集上一致收敛到连续函数，且短程与一圈公式 `g_ren²(p) = 1/(2b₀ ln(p/Λ)) + O(1/ln²)` 相符（`b₀=11N/(48π²)`）；
> - **(N2)［非平凡轨迹］** `Λ > 0`（不塌缩到 Gaussian 平凡轨迹）。
>
> 则 `∀ p ∈ (Λ,∞): g_ren^μ(p) ≠ 0`，从而 `¬ IsGaussianMeasure μ`。

### 6.3 为何 `g_ren` 是唯一合格见证

`g_ren(p)` 定义为**规范不变的顶角/传播子比**（MOM 型）：规范依赖在比值中**相消**，得规范不变；且自由理论所有 `n≥3` 连通顶角为零 ⟹ `g_ren ≡ 0`。故它同时满足命题 1.2 的 (i) 与 (ii)，是"某个高阶连通关联非零"的**规范不变化身**。

### 6.4 引理链与依赖

| 引理 | 内容 | 依赖 | 状态 |
|---|---|---|---|
| **C** | Gaussian ⟺ 高阶累积全零（偶性消奇阶；但**每偶阶 4,6,… 都需**；4 点单独不足） | 测度论 + OS1 | **Lean 部分已证**：`NonTrivialCore.lean`（`conn4_eq_zero_iff`、`four_point_insufficient`、`radem_K4_ne_zero`，均无公理）；解析方向待接入 Mathlib |
| **C′** | `NonTrivial ↔ ∃ n≥3` 非零累积见证（标量） | C | 解析方向待 Mathlib；代数骨架已 Lean |
| **G** | `g_ren ≡ 0 ⟺ μ` Gaussian（纯规范） | C + 唯一性 | 可证（定义层） |
| **S** | scheme 变换保持 `g=0`；故 `g_ren≠0` scheme-无关 | 代数 | 平凡（已验证） |
| **AF** | `β(g)=−b₀g³+…`，`b₀>0` | 微扰 + Balaban | 2026 输入 |
| **DT** | `Λ = μ e^{−1/(2b₀g²)}`，`Λ>0 ⟺ g_ren≠0 ∀μ>Λ` | AF | 定义层 |
| **N1** | 短程 `g_ren(p)` 一致收敛、匹配一圈 | **2026 UV 输入** | **待审计** |
| **N2** | 极限沿 `Λ>0` 非平凡轨迹 | **轨迹定义 / 子列选择** | **待审计** |
| **T** | 定理 N | N1+N2+DT+G | 由上述拼出 |

**逻辑定性**：C、C′、G、S、DT 是**可证/定义层**（近期 Lean 目标）；AF 是既有输入；**真正的"非平凡性内容"全部压在 N1、N2**——而它们是**审计问题**（§5.2），不是新数学。这是本文最诚实也最有用的结论。

---

## 7. 障碍（按严重度排序）

1. **O1 复合算子陷阱（致命但已诊断）。** §1.2：规范不变连通 4 点对自由场非零。任何"连通 4 点 ≠ 0"型的非平凡性陈述都是**空洞**的。必须改用 `g_ren` 或固定规范下的 `A`。
2. **O2 渐近自由是"敌人"。** UV 局域可观测量随 `a→0` 趋于自由值（`⟨tr U_□⟩→1` 等）；纯 UV 局域见证**检测不到**相互作用。见证必须落在**物理（IR）尺度**或**跨尺度积分量**（`Λ`、凝聚、`g_ren`）。
3. **O3 `Λ>0` 是真正内容。** 非平凡性 ⟺ 量纲 transmutation ⟺ `Λ>0`。这是**全局/IR** 陈述；虽**弱于** MG（只需"存在一个非零尺度"，不需一致指数速率），但**共享 IR 控制需求**。
4. **O4 子列 / 轨迹塌缩 vacuity。** 子列极限 + 未排除 `g₀≡0` 平凡轨迹 ⟹ 陈述可能被自由理论满足。这是 2026 陈述的**首要审计风险**。
5. **O5 迹反常"结构 ≠ 定标"。** 反常**方程**成立 ≠ 系数/凝聚非零；自由场（`β=0`）也满足结构性方程（右边为 0）。
6. **O6 OPE 不判别。** §4.2：OPE 系数自由也不为零。
7. **O7 Lean 层是标量场。** `Basic.lean` 用标量 `φ`；纯规范需复合算子 OS 层 + `g_ren` 见证，尚未搭建。

---

## 8. 下一步（最小可验证动作）

**短期（1–2 周，可独立完成，不碰 MG）**

- **P0（Lean）✅ 核心已做**：`NonTrivialCore.lean`（core Lean，秒级编译，`#print axioms` 无公理）已证：`conn4_eq_zero_iff`；Rademacher `κ₄=−2≠0`（4 点能见证其非 Gaussian）；**`four_point_insufficient`**（偶、`κ₄=0`、`κ₆≠0`，证明 4 点单独不足）。**待补**：接 Mathlib 证解析方向（`κₙ=0 ∀n≥3 ⟹` Gaussian）；把 `NonTrivial` 与"(`n≥3` 全体)非零累积"等价化（注意 §2.3：**不能**收紧到 n=4）；自由 GFF `NonTrivial=false` 的 falsification 测试（需 OSforGFF 构建）。
- **P1（陈述层，纯规范）**：把见证从标量 `φ` 换成 `g_ren`（抽象为"RG 荷 = 规范不变顶角/传播子比"），并**显式证明**自由 ⟹ 荷为零（引理 G）。验收：`isGaussian → rgCharge = 0` 无 sorry。
- **P2（审计，读而非建）**：读 `yang-mills-group/YM` 的陈述签名与 16 输入列表，回答 §5.2 三问：迹反常是结构还是定标？是否显式要求 `Λ>0`／nonzero `g_ren`？子列自由塌缩是否排除？`MassGap` 是否出现？（对应 `docs` 的 M0。）

**中期（与主 session 协同）**

- **P3（定理 N 的 N1）**：把"短程 `g_ren(p)` 一致收敛 + 匹配一圈"从 one-loop 输入推成引理；这是**UV 侧**，期望 2026 UV 机制可覆盖。若成立，则非平凡性**不依赖 IR**，与 MG 解耦。
- **P4（数值 certificate，非严格）**：用 `numerics/` 已有基础设施，在**梯度流耦合** `g²_GF(t)`（handoff 候选 #4，非微扰全尺度）上数值认证"在某个物理尺度 `g²_GF ≠ 0`"沿 AF 轨迹与外推极限稳定。作**非严格但具体**的见证，支撑 P3。

**长期（若 P2 显示只给结构性方程）**

- 非平凡性降级为与 MG 同源的 IR 问题，转攻 `Λ>0` 的严格证明（= 非平凡 AF 不动点的存在性 + 标度标定），与 handoff 的 `q(μ)=m/μ` 判据合并。

---

## 9. 对 `Basic.lean` 的具体修改建议

1. **保留** `IsGaussianMeasure` / `NonTrivial`（抽象层正确）。
2. **新增** `CumulantZeroOfOrderAtLeast n μ`，并**证明**引理 C：`IsGaussianMeasure μ ↔ ∀ n≥3, CumulantZeroOfOrderAtLeast n μ`；偶性下加推论 `↔ CumulantZeroOfOrderAtLeast 4 μ`。
3. **标注** `ConnectedFourPoint` 的适用范围：**仅标量**。对纯规范，**新增** `RGCharge`（抽象为 `SchwartzTestFunction → ℝ` 或 `μ → ℝ`，表示 `g_ren(μ_phys)`），并证 `IsGaussianMeasure μ → RGCharge μ = 0`（引理 G）。**在文档中显式写明 §1.2 陷阱**，防止未来把 `tr F²` 之类的连通 4 点误用作见证（否则谓词空洞）。
4. **新增** falsification 静例（在 `Smoke` 或 `Tests` 文件）：自由 GFF ⟹ `¬NonTrivial`；相互作用 toy ⟹ `NonTrivial`。

（以上为**陈述/验证骨架**，无需建立重型依赖；引理 C 的量测论部分可复用 OSforGFF 的多项式矩与唯一性。）

---

## 10. 参考文献

**平凡性**
- M. Aizenman, *Geometric analysis of φ⁴ fields...*, CMP 86 (1982).
- J. Fröhlich, *On the triviality of λφ⁴₄...*, Nucl. Phys. B 200 (1982).
- M. Aizenman, H. Duminil-Copin, *Marginal triviality of the scaling limits of critical 4D Ising and λφ⁴ models*, Ann. Math. 194 (2021), arXiv:1912.07973.
- R. Bauerschmidt, D. Brydges, G. Slade, *Introduction to a renormalisation group method*, LNM 2242 (2019), arXiv:1403.7424.

**渐近自由 / RG**
- D. Gross, F. Wilczek, PRL 30 (1973); H. D. Politzer, PRL 30 (1973).
- T. Balaban, *Renormalization group approach to lattice gauge field theories I–IV*, CMP 1980s.
- （主 session 的 `research/physics_handoff.md` 有当前的 scheme/交接分析，本文不重复。）

**OPE / 反常**
- K. Wilson, *Non-Lagrangian models of current algebra*, Phys. Rev. 179 (1969).
- J. Collins, *Renormalization*（OPE 与 MOM scheme 的标准参考）。
- M. Shifman, A. Vainshtein, V. Zakharov（胶子凝聚），Nucl. Phys. B 147 (1979); Nielsen (1977).

**OS / 形式化**
- K. Osterwalder, R. Schrader, CMP 31 (1973) 83; 42 (1975) 281.
- M. R. Douglas, S. Hoback, A. Mei, R. Nissim, *Formalization of QFT*, arXiv:2603.15770 (2026); `github.com/mrdouglasny/OSforGFF`.

**2026（据项目调研转述，需审计）**
- Yang–Mills Working Group (Gross 主持), *Validation of the OpenAI Yang–Mills construction*, 2026-11-10.
- OpenAI, `yang-mills-group/OpenAI-YM-20261028`; `yang-mills-group/YM`.
- hexagonmath, *From Four-Dimensional Wilson Lattice Gauge Theory to Quantum Yang–Mills Theory*, 2611.00492v1 (2026-11-15).
- 说明：2026 条目的 arXiv 编号与细节**未经本文独立核实**，仅作线索；非平凡性主张的具体形式（§5.2）以审计为准。

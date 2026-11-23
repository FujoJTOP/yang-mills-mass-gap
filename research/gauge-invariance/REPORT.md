# 规范不变性、协变 RG 与控制量：精确框架、引理链、障碍与下一步

> 2026-11-22 · 范围：4D 纯 Yang–Mills（紧单群 G），Wilson 格点 + 标度极限；目标是把质量间隙（MG）与 RG 控制量写成规范不变/规范协变的形式，并把规范固定构造接到规范不变可观测量上。
>
> **结论边界**：本报告不证明 MG，也不改进 `q₀ ≥ q₁` 交接的解析估计。它给出陈述层与结构层，并把"规范记账"从解析硬点中分离出来。所有断言都带状态标签：`【定理】` 文献已证、`【标准】` 标准但本次未逐条核实、`【本文】` 本报告的推导、`【开放】` 未解、`【待核实】` 需进一步查证。

---

## 0. 摘要（结论先行）

**主论点：规范不变性不是 4D 质量间隙的解析障碍，而是可分离的结构/记账问题。**

1. **测度层**：有限格点 Wilson 测度 `μ_β` 在 `𝒢_Λ = G^{V(Λ)}` 下不变（Haar 双不变 + 迹不变）。它下降到轨道空间 `𝒜/𝒢` 上的概率测度 `μ̄`。物理理论是 `μ̄`，不是任何规范切片上的测度。`【定理/初等】`
2. **RG 层**：存在**精确等变**的块自旋 `Φ`（面平均 + 平行移动到块中心），配块规范同态 `b: 𝒢_Λ → 𝒢_{Λ'}`，满足 `Φ(g·U) = b(g)·Φ(U)`。于是 pushforward 自动 `𝒢_{Λ'}`-不变，RG 在轨道空间上良定义。**Gribov 被绕开，而不是被解决**：全程不需要全局截面。`【本文构造 + 标准测度论】`
3. **控制量层**：`q = m/μ`、谱隙 `Δ`、关联长度 `ξ` 都定义为**不变局部算子**的连通关联，故是 `μ̄` 的函数，规范不变。规范固定只作为**局部**计算工具（Balaban 的小场规范），通过"不变可观测量的规范固定独立性"桥接到不变量的计算。`【定理/标准 + 本文精确化】`
4. **真正的开放点**：(a) Balaban 型块自旋的**定量类保持**（4D SU(N) 的小场/大场控制）；(b) `q₀ ≥ q₁` 交接估计。两者都与规范无关。`【开放】`
5. **陈述层手术**：`lean/YMStatement/Basic.lean` 目前是**标量场**陈述（`FieldConfiguration d = 𝒮′(ℝᵈ)`），不是规范理论陈述；`NonTrivial = ¬Gaussian` 也不对（β=0 的乘积 Haar 测度非 Gaussian 但平凡）。文件里已经定义了 `ConnectedFourPoint` 却没用；改用它即可修复非平凡性。`【本文诊断】`

**一句话**：把"规范固定 + 控制量"换成"轨道空间测度 + 等变块自旋 + 不变控制泛函"，规范不变性变成构造性的，而不是需要额外证明的。

---

## 1. 问题定位

### 1.1 三件必须分开的事

| # | 命题 | 性质 | 状态 |
|---|---|---|---|
| (a) | 测度 `μ_β` 在规范群下不变 | 有限维积分/Haar | 已证、初等 |
| (b) | RG 映射把规范不变理论送到规范不变理论 | 结构 + 定量 | 结构【本文】、定量【开放】 |
| (c) | 控制量 `q, Δ, ξ` 是规范不变的 | 定义在不变代数上 | 由 (a)(b) 推出 |
| (d) | 规范固定构造能计算 (c) | FP/带号公式 | 有限格点已证、连续统形式化 |

把 (a)–(d) 混在一起谈，是"规范不变性很难"这一印象的主要来源。实际上只有 (b) 的**定量**部分是真难，而它不是规范问题。

### 1.2 为什么不能只在规范固定后工作

朴素路线：选一个规范（Landau/Coulomb/轴规范），在切片 `𝒮` 上定义测度 `e^{-S} Δ_FP ∏dA`，然后做 RG、算关联。问题：

- **Singer 定理**：非阿贝尔规范理论在紧化 ℝ⁴（S⁴）上**不存在全局连续截面** `𝒜/𝒢 → 𝒜`。任何"处处良好"的规范固定都不可能。`【定理】`（Singer 1978）
- **Gribov 拷贝**：切片与轨道相交多次；`Δ_FP` 变号；朴素 `|Δ_FP|` 高估/低估；BRST 规范固定的配分函数按 `χ(G)=0` 消失（Neuberger 零点），得到 0/0。`【标准】`
- **结论**：规范固定**不能**用来定义理论。它只能作为计算不变可观测量的工具，并且必须带正确的（带号）轨道权重。

正确的定义顺序是反过来的：先在 `𝒜` 上定义不变测度 → 下降到 `μ̄` → 在 `μ̄` 上定义物理量 → 需要计算时才用局部规范固定。

### 1.3 `Basic.lean` 的三处问题（诊断，不是否定）

当前文件（`D:\Dev\math\lean\YMStatement\Basic.lean`）的结构是：

```lean
def IsGaussianMeasure (μ : ProbabilityMeasure (FieldConfiguration d)) : Prop := ...
def NonTrivial (μ) : Prop := ¬ IsGaussianMeasure μ
def ExponentialClustering (μ) (m) : Prop :=
  m > 0 ∧ ∀ f g, ∃ C, ∀ a, ‖GJ μ (f + g.translate a) - GJ μ f * GJ μ g‖ ≤ C * exp (-m * ‖a‖)
structure MassGapStatement (μ) : Prop where
  os : SatisfiesAllOS μ
  nonTrivial : NonTrivial μ
  gap : ∃ m, ExponentialClustering μ m
```

问题：

1. **配置空间是标量场**。`FieldConfiguration d = 𝒮′(ℝᵈ)` 是标量场的分布空间；Yang–Mills 的配置空间是**联络模规范**（或格点上的 `G^{E}` 模 `G^{V}`）。文件里没有群 `G`、没有链接变量、没有曲率、没有 Wilson 环。所以 `MassGapStatement` 目前不是 YM 陈述。
2. **`NonTrivial = ¬Gaussian` 判据不对**。β=0 的格点规范理论（每条链接独立 Haar）是**乘积测度**：它非 Gaussian，但它是最平凡的理论，且（强耦合下）有质量间隙。反过来，一个有质量的广义自由场是 Gaussian 的、有间隙的、平凡的。正确的判据是"存在非零的连通 n≥3 矩"——它同时排除 Gaussian 与乘积测度。文件里已经定义了 `ConnectedFourPoint`，但 `NonTrivial` 没有用它。`【本文】`
3. **`ExponentialClustering` 的源不是规范不变的**。`GJ μ (f + g.translate a)` 用的是任意测试函数 `f`，对应源 `j·φ`；在规范理论里 `j·A` 依赖规范。物理的源必须是规范不变的（局部复合算子 `tr F²` 的 smear，或 Wilson 环）。而且**质量间隙是关于局部算子的**；Wilson 环的关联由弦张力/面积律支配，不是同一回事（见 §5.3）。

修复方案见 §6。

---

## 2. 精确框架

### 2.1 格点数据

**共享构造（逐字采用 A0 `research/gluing/REPORT.md` §1.1–§1.4）**：`G` 紧、单、单连通（结构群取单连通覆盖 `G̃`）；格点环面 `Λ_{L,a} = (ℤ/Nℤ)⁴`，`N := L/a ∈ ℕ`，即 `T⁴_L`（周期边界），固定 `θ=0`、平凡 't Hooft 通量扇区；每边一个 Haar 概率测度 `dU_e`。为简洁，下文记 `Λ := Λ_{L,a}`，边集 `E(Λ)`，点集 `V(Λ)`。

```
𝒜_Λ = G^{E(Λ)}                    （链接变量 U = (U_e)_{e∈E}）
𝒢_Λ = G^{V(Λ)}                    （规范变换 g = (g_x)_{x∈V}）
(g·U)_{x,y} = g_x U_{x,y} g_y^{-1}          （x,y 为边 e 的端点）
```

Wilson 作用量与测度（`S_β ≡ S_W(β)`，A0 §1.2；**`N_c := dim f` 为定义表示维数**，`SU(N)` 时 `N_c = N`；`tr` 为未归一化定义表示迹，`tr_f 1 = N_c`）：

```
S_β(U) = β Σ_p (1 - (1/N_c) Re tr U_p),      U_p = 绕 plaquette p 的 ordered product
β = 2N_c / g²,      γ := β/N_c = 2/g²          （进入测度的唯一组合是 γ）
μ_{L,a,β} = Z_{L,a}(β)^{-1} e^{-S_β(U)} ∏_{e∈E(Λ)} dU_e      （dU_e = Haar）
```

**RT 声明**：`β(a)` 沿 A0 §1.4 的重整化轨迹；本文与 A0 §1.4 用**同一条 RT、同一个 scheme 常数 `C`**，不引入独立轨迹（不同 `C` 是不同理论；A0 §1.4、GL1）。

**引理 2.1（测度不变性）** `【定理/初等】`
对每个 `g ∈ 𝒢_Λ`，`(g·)_* μ_{L,a,β} = μ_{L,a,β}`。
*证明*：`tr U_p` 在 `U_p ↦ g_x U_p g_x^{-1}` 下不变；`∏ dU_e` 在 `U_e ↦ g_x U_e g_y^{-1}` 下不变（Haar 左/右不变性逐边相乘）。∎

**引理 2.2（标度极限保持不变性）** `【标准】`
若 `μ_{a}` 是一列 `𝒢`-不变测度，弱收敛到 `μ`，则 `μ` 不变（对连续有界不变 `F` 取极限）。格点→连续统时规范群随格点变化，需要一致的连续统规范作用；这一点在 2026 构造里应被显式检查（见 §8.2 审计清单）。

### 2.2 不变可观测量代数

对闭合环 `C ⊂ Λ`：

```
W_C(U) = (1/N) tr ∏_{e∈C} U_e        （Wilson 环；trace 使其与基点无关）
𝔄_Λ = 由 {W_C} 生成的代数
𝔄_Λ^loc = 由有界尺寸的环生成的局部子代数
```

**引理 2.3（Giles 型重构）** `【定理/标准；本次未逐条核实 Giles 1981】`
对 `G = U(N)`（及 SU(N)），`𝔄_Λ = C(𝒜_Λ)^{𝒢_Λ}`：不变连续函数由 Wilson 环生成，且分离轨道。对一般紧 `G`，用所有不可约表示的迹（spin network）生成。

用途：测度 `μ̄` 由它在 `𝔄_Λ` 上的值唯一确定；"物理 = 不变矩"。

### 2.3 轨道空间测度

`π: 𝒜_Λ → 𝒜_Λ/𝒢_Λ` 为商映射。

**引理 2.4（轨道测度）** `【定理/标准】`
`𝒜_Λ` 紧度量、`𝒢_Λ` 紧，故 `𝒜_Λ/𝒢_Λ` 紧度量、`π` 固有。`μ` 不变 ⟹ 存在唯一概率测度 `μ̄` 使 `π_*μ = μ̄`。对不变有界 `F`：

```
∫_{𝒜} F dμ = ∫_{𝒜/𝒢} F̄ dμ̄,      F̄(ω) = F(任意代表元)
```

**可测截面**存在（标准 Borel 选择），但**连续截面不存在**（Singer）。这不影响定义，只影响"用切片计算"。

### 2.4 规范不变的 OS 层

**定理 2.5（Osterwalder–Seiler 1978）** `【定理】`
Wilson 测度对**规范不变可观测量**满足反射正性：时间反射 `Θ` 作用在链接上（并把规范变换一起反射），对正半空间代数的 `F`，

```
⟨ΘF, F⟩_{L²(μ)} ≥ 0.
```

由此得到转移矩阵 `T` 作用在

```
ℋ = L²(𝒜_spatial/𝒢_spatial, μ_spatial)
```

上，`T ≥ 0`，真空 `Ω = 1`，`TΩ = Ω`。`ℋ` 是**规范不变** Hilbert 空间；规范变换在 `ℋ` 上平凡作用。

注：`Menotti–Pelissetto (1988)` 把正性推广到过格点的反射面，对**所有规范不变可观测量、所有时间间隔**成立。`【标准】`

### 2.5 质量间隙的规范不变定义

**定义 2.6（质量间隙）**
在 `ℋ` 上，令 `Δ(μ̄) = inf spec(-log T |_{Ω^⊥})`（格点单位）。等价地：`∃m>0`，使得对**每个局部规范不变算子** `O`（如 smeared `tr F_{μν}F^{μν}`），存在常数 `C_O` 使

```
|⟨O(x) O(y)⟩_c| ≤ C_O e^{-m |x-y|}.
```

`m` 与 `O` 的标号无关（label-uniform）；`C_O` 可以依赖 `O`。`【标准等价；形式化见 §4 L1.1】`

**规范不变性**：`T` 作用在不变子空间，`Δ` 只依赖 `μ̄`。`【本文，由 2.4–2.5】`

**必须强调的三点**（详见 §5.3）：
- 用**局部**算子，不用 Wilson 环（后者是面积律，不是质量间隙）。
- 用**连通**关联（减去真空期望）。
- 有质量的广义自由场也有间隙；所以间隙本身不排除平凡理论，非平凡性要单独陈述（§6.3）。

### 2.6 等变块自旋与 RG 交换图

**定义 2.7（等变块自旋）**
块自旋是映射 `Φ: 𝒜_Λ → 𝒜_{Λ'}`（`Λ' = LΛ`，`L` 为块边长），配群同态 `b: 𝒢_Λ → 𝒢_{Λ'}`，满足

```
(E)   Φ(g·U) = b(g)·Φ(U)      ∀ g ∈ 𝒢_Λ, U ∈ 𝒜_Λ.
```

`b` 称为**块规范变换**。若 (E) 只在误差内成立，记

```
δ(g,U) := Φ(g·U) · (b(g)·Φ(U))^{-1}   （等变缺陷）
```

**构造 2.8（面平均块自旋，精确等变）** `【本文】`
取块中心 `X`（每个块一个），固定块内路径。对相邻块 `X, Y` 的公共面 `F_{XY}`，定义粗链接

```
V_{XY} = Avg_{x∈F_{XY}} [ P_X(x) · U_{x,y} · P_y(Y) ],
```

其中 `y = x + μ̂` 是跨面的细边，`P_X(x)` 是从 `x` 到 `X` 的块内平行移动，`P_y(Y)` 是从 `Y` 到 `y` 的块内平行移动，`Avg` 是**共轭等变**的平均（Karcher 均值，或带共轭不变截断的对数平均）。

**引理 2.9（精确协变）** `【本文】`
在规范变换 `g` 下，块内中间点的 `g` 全部相消，端点剩下 `X` 与 `Y`：

```
V_{XY} ↦ g_X V_{XY} g_Y^{-1}.
```

故取 `b(g)_X = g_X`，`Φ(g·U) = b(g)·Φ(U)` **精确成立**，缺陷 `δ ≡ 1`。`b` 是限制映射 `g ↦ g|_中心`，是满同态。

*证明*：`P_X(x) ↦ g_X P_X(x) g_x^{-1}`，`U_{x,y} ↦ g_x U g_y^{-1}`，`P_y(Y) ↦ g_y P_y(Y) g_Y^{-1}`，乘积中 `g_x, g_y` 相消。平均的共轭等变性给出 `Avg(g T g^{-1}) = g Avg(T) g^{-1}`。∎

**为什么朴素块平均失败** `【标准】`
若令 `Ā_μ(X) = (1/|B|) Σ_{x∈B} A_μ(x)`，`A = log U`，则在 `g` 下

```
A_μ(x) ↦ g_x A_μ(x) g_x^{-1} + g_x ∂_μ g_x^{-1}.
```

第一项的平均**不是** `Ad_{ḡ_X}` 作用在平均上（除非 `g` 在块内常数）；第二项平均后不是任何单个 `ḡ_X` 的 `∂_μ ḡ_X^{-1}`。所以不存在 `ḡ_X` 使 `Ā` 正确变换。**块内路径 + 端点相消**是修复的关键。

**引理 2.10（pushforward 不变）** `【本文】`
若 `μ` 是 `𝒢_Λ`-不变，`ν = Φ_*μ`，则 `ν` 是 `𝒢_{Λ'}`-不变（在 `b` 满射时；限制映射满足）。
*证明*：对 `h ∈ 𝒢_{Λ'}`，取 `g` 使 `b(g)=h`。对连续有界 `f`：

```
∫ f(h·V) dν(V) = ∫ f(h·Φ(U)) dμ(U)
              = ∫ f(Φ(g·U)) dμ(U)
              = ∫ f(Φ(U)) dμ(U)          （μ 不变）
              = ∫ f(V) dν(V).            ∎
```

**引理 2.11（有效作用量不变）** `【本文】`
若 `ν` 关于粗 Haar 测度有密度 `ρ = dν/dHaar`，则 `ρ(h·V) = ρ(V)`，故 `S_eff = -log ρ` 是粗规范不变的。

**引理 2.12（轨道空间 RG）** `【本文】`
`Φ` 等变 ⟹ 诱导良定义的 `Φ̄: 𝒜_Λ/𝒢_Λ → 𝒜_{Λ'}/𝒢_{Λ'}`，且 `ν̄ = Φ̄_* μ̄`。交换图：

```
        Φ
 𝒜_Λ ────────► 𝒜_{Λ'}
  │ π              │ π'
  ▼                ▼
 𝒜_Λ/𝒢_Λ ─────► 𝒜_{Λ'}/𝒢_{Λ'}
        Φ̄
```

测度层同样交换：`μ̄' = Φ̄_* μ̄`。

**注 2.13（定量部分是唯一的难点）** `【开放】`
上述结构（2.9–2.12）是初等的。Balaban 的真正内容是**定量类保持**：证明 `S_eff` 留在一个小场可分析的类里，并给出误差界。`Φ` 的存在与等变性不是问题；`Φ` 的正则性、`log` 分支、大场行为才是问题（见 §3.1、§7）。**规范不变性不进入这一层。**

### 2.7 控制泛函 `q`

**定义 2.14** 令 `μ` 为尺度（物理动量标度），`m(μ)` 为尺度-`μ` 有效理论的关联长度倒数（用不变局部算子的连通关联定义，见 §5.1）。定义

```
q(μ) = m(μ)/μ.
```

等价地，`q = a_k · m_phys`，即第 `k` 层格点单位下的无量纲质量。`q` 是 `μ̄_k` 的泛函，故**规范不变**。`【本文】`

`q` 沿 RG 流演化；它是**坐标**，不是自由参数（见 §5.4 的循环性讨论）。

---

## 3. 五条路线的精确对照

### 3.1 Balaban：规范协变块自旋

**文献**（本次核实）：
- T. Bałaban, *Averaging operations for lattice gauge theories*, Comm. Math. Phys. **98** (1985) 17–51. 摘要明确：研究**格点规范场与规范变换**的平均操作，刻画平均操作正则（解析）的配置类。`【定理】`
- T. Bałaban, *Spaces of regular gauge field configurations on a lattice and gauge fixing conditions*, Comm. Math. Phys. **99** (1985) 75–102. `【定理】`
- T. Bałaban, *Renormalization group approach to lattice gauge field theories I. Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions*, Comm. Math. Phys. **109** (1987) 249–301. `【定理】`
- II. *Cluster expansions*, Comm. Math. Phys. **116** (1988) 1–22. `【定理】`
- *Ultraviolet stability of three-dimensional lattice pure gauge field theories*, Comm. Math. Phys. **102** (1985) 255–275.
- *Large field renormalization I/II*, Comm. Math. Phys. **122** (1989) 175–202 / 355–392.
- 其他同系列：CMP **95** (1984) 17–40; **96** (1984) 223–250; **119** (1988) 243–285.

**结构**（本文的解读）：
- Balaban 的 RG 变换 = 链接与规范变换的**平均操作**（CMP 98），作用在**局部化有效作用量**上（CMP 109），用 cluster expansion 控制（CMP 116）。
- 小场区：在每个块内做**局部规范固定**，使块链接接近单位元；局部规范固定在小场邻域内唯一，**无 Gribov 拷贝**。大场区：单独估计（CMP 122）。
- 关键：有效作用量是**规范不变**的局部泛函之和（Wilson 环展开），RG 保持这个类。

**哪些是定理、哪些不是**：
- 已证：UV 稳定性、小场有效作用量的控制、耦合常数重整化、cluster expansion。
- 未证：IR / 无限体积极限 / OS 公理 / 质量间隙（Balaban 自己也没有声称）。`【开放】`
- **精确定量等变界**（缺陷 `δ` 的界、类保持的常数）在文献中技术性很强，未被独立形式化或逐条审计。`【待核实】`

**与现代重构的关系**：`Dimock, The renormalization group according to Balaban I–III`（Rev. Math. Phys. 25 (2013) 1330010；J. Math. Phys. 54 (2013) 092301；Ann. Henri Poincaré 15 (2014) 2133–2175）是**φ⁴₃**的讲解，**不是规范理论**的讲解。可作为方法参考，不能当作 YM 的现成证明。`【本次核实】`

**本文的定位**：精确等变块自旋是初等的（§2.6）；Balaban 的贡献是**定量类保持**，它可以在"等变"这一结构假设之下独立陈述与攻关。这给了我们一个干净的分解（见 §4 的 L2.5）。

### 3.2 随机量子化在轨道上（CCHS）

**文献**（本次核实）：
- A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Langevin dynamic for the 2D Yang–Mills measure*, Publ. Math. IHÉS **136** (2022) 1–147, arXiv:2006.04987.
- A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Stochastic quantisation of Yang–Mills–Higgs in 3D*, Invent. Math. **237** (2024) 541–696, arXiv:2201.03487.
- I. Chevyrev, H. Shen, *Uniqueness of Gauge Covariant Renormalisation of Stochastic 3D Yang–Mills–Higgs*, Arch. Rational Mech. Anal. **250**(1) (2026). `【标题已核实；arXiv 号未核实】`

**2D 的机制**（精确）：
- 状态空间 = 分布联络的**轨道空间** `𝒜/𝒢`，在自然 Hausdorff 度量下是 Polish 空间。
- Langevin 方程 `dA = -∇S(A) dt + dW` 沿规范方向**退化**（`S` 规范不变 ⟹ `∇S` 与轨道正交）。用**时间依赖规范变换（DeTurck trick）**把方程化为严格抛物，证明解**依分布规范协变**（gauge covariance in law）：初值律不变 ⟹ 解律不变。
- 不变测度 = 2D YM 测度；2D 的 OS 由精确构造（Driver/Lévy/Sengupta）给出。

**3D 的机制**：
- regularity structures 构造重整化流；counterterm 的选择必须使解**依分布规范协变**；CCHS 证明存在唯一这样的选择（Chevyrev–Shen 2026 加强唯一性）。
- 这是"规范协变性 = 重整化方案的约束"的范例。

**4D 的障碍**：**不是规范不变性**，而是无穷多 counterterm / 非超可重整。规范方向的退化在 4D 同样无害；困难在横向方向。`【开放，但与规范无关】`

**对 OS 的地位**：CCHS 给的是**不变测度/过程**，不自动给 OS 正性。2D 因精确构造而有 OS；3D YM–Higgs 的 OS 未由 CCHS 建立。`【标准】`

### 3.3 Osterwalder–Seiler

**文献**：K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Physics **110** (1978) 440–471; E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, LNP 159 (1982). `【本次核实】`

**内容**：
- Wilson 作用的反射正性（规范不变可观测量）。
- 转移矩阵/哈密顿量在规范不变子空间。
- 强耦合**体积一致**的规范不变关联指数聚类，速率 `m(β)>0`；唯一真空；无限体积极限。
- 这正是 **MG 的形状**——只不过只在强耦合。`【定理】`

**意义**：OS 1978 已经示范了"规范不变的 mass gap 陈述"长什么样。我们的 MG 陈述应复刻它的形式，把"强耦合"换成"所有耦合 / 标度极限"。

### 3.4 Gribov 问题

**文献**：V. N. Gribov (1977, SLAC-Trans-0176); I. M. Singer, *Some remarks on the Gribov ambiguity*, Comm. Math. Phys. **60** (1978) 7–12. `【本次核实】`

**内容**：
- Singer：非阿贝尔规范理论在 S⁴ 上商映射**没有全局连续截面**；Gribov 歧义不是 Coulomb 规范特有，而是所有规范都有。
- FP 拷贝：切片与轨道相交多次，`det` 变号；限制到第一 Gribov 区或使用带号求和是标准补救；BRST 规范固定的配分函数按 `χ(G)=0` 消失（Neuberger 零点）。
- **对不变可观测量**：拷贝的贡献在推到轨道空间时坍缩；`O[A]` 不变 ⟹ 与规范/拷贝无关。非不变量（胶子/鬼传播子）依赖规范——它们是规范赝像。

**对 RG 的意义（本文的关键区分）**：
- **UV / 小场**：局部规范固定（Balaban 的块内规范）在小场邻域唯一，**无 Gribov 拷贝**。所以 UV RG 的规范固定是合法的。
- **IR / 大场**：Gribov  horizons 可能出现；这是"大场问题"的一部分，不是 UV 问题。
- **结论**：Gribov 不阻碍 Balaban 的 UV RG；它阻碍的是"用一个全局规范固定定义整个理论"的幻想。`【本文】`

### 3.5 规范固定独立性

**定理 3.1（有限格点）** `【定理/初等；本文给出形式化版本】`
设 `G` 有限，`μ` 在 `𝒢_Λ` 下不变，`F` 不变。对每条轨道 `ω` 取一个代表元 `s(ω)`（任意截面），则

```
∫_{𝒜} F dμ = Σ_{ω ∈ 𝒜/𝒢} μ(ω) · F(s(ω)).
```

特别地，右端与截面 `s` 无关：任意两个规范固定给出相同的不变期望。若写成"切片 + 权重"：

```
∫_{𝒜} F dμ = Σ_{ω} [ μ(ω)/|ω| ] · |ω| · F(s(ω)),
```

其中 `|ω|` 是轨道体积，起 Faddeev–Popov 行列式的作用。

*证明*：`F` 在 `ω` 上常数；`μ|_ω = μ(ω) · Haar_ω`（不变性 ⟹ 轨道上均匀）；求和。∎

**定理 3.2（紧群版本）** `【标准；连续统形式化】`
对紧 `G`，`Σ_ω` 换成 `∫_{𝒜/𝒢}`，`|ω|` 换成 `Δ_FP`（带号），得到

```
∫_{𝒜} F dμ = ∫_{𝒮} F Δ_FP dμ_𝒮      （F 不变，带号 FP 权重）.
```

若切片与轨道相交多次且丢掉符号，等式失败；此时"规范固定测度"不是正测度，非不变量的"期望"是规范赝像。

**推论 3.3（桥）** `【本文】`
控制量 `q`、谱隙 `Δ`、关联长度 `ξ` 都是不变关联的泛函，故**可以在任意规范固定下计算**（带正确 FP 权重），结果相同。这就是"把规范固定构造接到规范不变可观测量"的精确内容。

**推论 3.4（MRG 的地位）** `【本文】`
Magnen–Rivasseau–Sénéor (CMP 155 (1993) 417–430) 的 YM₄（带 IR 截断、规范固定、SU(2)）构造，只有在"计算不变可观测量"的意义上才是物理的；其 Ward 恒等式正是桥的体现。它不证明 OS / 质量间隙，且规范固定是微扰/弱场意义下的（Gribov 未处理）。`【标准】`

### 3.6 现代功能不等式路线

**SZZ（强耦合、任意维）**：Shen–Zhu–Zhu, Comm. Math. Phys. **400** (2023) 805–851, arXiv:2204.12737. `【本次核实】`
- Langevin 动力学；整条链空间上的唯一不变测度、指数遍历、LSI/Poincaré（Bakry–Émery）。
- 一大类可观测量的关联指数衰减 ⟹ 严格正的质量间隙。

**Bauerschmidt–Dagallier（U(1)、4D）**：*Log-Sobolev inequality for U(1) lattice gauge theory in d=4*（2024；`arXiv 号未核实`，见 §9）。`【待核实】`
- 4D 物理维数上第一个规范理论的严格功能不等式；Polchinski 型多尺度 + 单尺度 LSI + telescoping。
- U(1) 阿贝尔，无 Gribov；推广到 SU(N) 是公开的下一步。

**"限制原理"** `【本文，标准论证】`
若 `μ` 在**整条链空间**上满足 LSI(c) 或 Poincaré，则限制到**规范不变子空间**仍成立：不变函数的规范方向梯度为零，Dirichlet 形式只看到横向方向；不变函数构成子空间，故不等式自动限制。于是得到不变子空间的谱隙，**无需规范固定**。
*注意*：整条空间的谱隙可能被规范模式压低甚至为零；物理间隙是**不变子空间**的间隙。这正是"必须在不变子空间上谈 gap"的又一个理由。

**限制**：这些结果都在强耦合/小 `β`（SZZ）或阿贝尔（BBD）。弱耦合 SU(N)、4D 仍开放。`【开放】`

---

## 4. 引理链

记号：`Λ` 细格点，`Λ'` 粗格点，`Φ` 块自旋，`b` 块规范同态，`μ` 测度，`μ̄` 轨道测度，`q` 控制量。难度按 Lean 形式化估计（易/中/难/开放）。

### 4.1 总表

| 编号 | 陈述（简） | 状态 | 难度 | 依赖 |
|---|---|---|---|---|
| L0.1 | 格点数据、Wilson 作用量、测度 | 定义 | — | — |
| L0.2 | `μ_β` 是 `𝒢_Λ`-不变 | 定理/初等 | 易 | Haar |
| L0.3 | 轨道测度 `μ̄` 存在唯一 | 定理/标准 | 中 | 商测度 |
| L0.4 | 不变函数由 Wilson 环生成 | 定理（Giles 型） | 中 | 不变量理论 |
| L0.5 | OS 正性 + 转移矩阵 | 定理（OS 1978） | 难 | 反射正性 |
| L0.6 | OS/不变性过弱极限 | 标准 | 中 | 测度收敛 |
| L1.1 | 质量间隙 ⟺ 不变局部算子指数聚类 | 标准等价 | 难 | 谱定理 |
| L1.2 | `Δ` 只依赖 `μ̄` | 本文 | 易 | L0.3, L1.1 |
| L2.1 | 面平均块自旋精确等变 | 本文 | 易 | 平行移动 |
| L2.2 | pushforward 保持不变 | 本文 | 易 | L2.1 |
| L2.3 | 有效作用量不变 | 本文 | 易 | L2.2 |
| L2.4 | 轨道空间 RG 良定义 | 本文 | 易 | L2.2 |
| L2.5 | Balaban 定量类保持（4D SU(N)） | 开放/文献 | 极难 | Balaban |
| L3.1 | `q` 的定义与不变性 | 本文 | 易 | L1.2 |
| L3.2 | `q` 对耦合 scheme 变换不变 | 本文 + 交接报告 | 中 | 物理质量 |
| L3.3 | `q` 随 `μ` 严格单调，`q=m` 处为 1 | 本文 | 易 | 微积分 |
| L3.4 | 交接引理：`q₀≥q₁ ⟹ gap≥min(δ_W,δ_S)` | 交接报告 | 易 | 实分析 |
| L4.1 | 有限格点规范固定独立性 | 定理/初等 | 易 | 轨道分解 |
| L4.2 | FP/带号公式（紧群） | 标准 | 中 | 带号测度 |
| L4.3 | 桥：`q` 可在任意规范固定下算 | 本文推论 | 易 | L4.1, L3.1 |
| L5.1 | 强耦合 gap（SZZ/OS） | 定理 | 难（引用） | LSI/Poincaré |
| L5.2 | MG 陈述（体积一致、标号一致） | 定义 | — | L1.1 |
| L5.3 | MG ⟹ 质量间隙 | 由 L1.1 | 中 | 谱理论 |
| L5.4 | 证明 MG | 开放 | — | L2.5, L3.4, L5.1 |

### 4.2 关键引理的详细陈述

**L0.2（测度不变）** 已在 §2.1 给出。

**L0.3（轨道测度）** 已在 §2.3 给出。形式化要点：紧群作用的商是标准 Borel 空间；不变测度下降到商；对不变函数用 `π_*` 的积分公式。

**L1.1（质量间隙 ⟺ 指数聚类）** `【标准；形式化难】`
在 OS + 重构下，`T` 在 `Ω^⊥` 上的谱隙 `Δ = -log(λ₂/λ₁)` 与不变局部算子的连通关联衰减率等价。方向"gap ⟹ 聚类"用谱定理（谱测度支撑在 `[Δ,∞)`，Fourier 变换指数衰减）；方向"聚类 ⟹ gap"用谱测度的下界（若无隙，则存在算子的关联不指数衰减）。技术细节：有限体积的 `T` 是有限维矩阵（格点），谱理论是有限维的；空间衰减由 Euclidean 不变性/反射正性给出。`【形式化难度：难，但有限维部分可控】`

**L2.1（面平均精确等变）** 已在 §2.6 给出。形式化要点：`Avg` 的共轭等变性是关键引理；对有限群可用显式求和。

**L2.2（pushforward 不变）** 已在 §2.6 给出。形式化要点：`b` 满射（限制映射）；换元 + `μ` 不变。

**L2.5（Balaban 定量类保持）** `【开放】`
陈述（示意）：存在常数 `c, C` 与"小场类" `𝔖`，使得若 `S_eff` 在 `𝔖` 中且 `β` 足够大，则 `T(S_eff)` 仍在 `𝔖` 中，且 `‖T(S_eff) - S_1-loop‖ ≤ C β^{-2}`（或相应的误差界），并且 `Φ` 在 `𝔖` 上正则（解析）、缺陷 `δ` 可控。
这是 Balaban 的核心内容；它**不涉及规范不变性**（等变性已由 L2.1 保证），只涉及正则性与类保持。

**L3.1（`q` 的不变性）** 已在 §2.7 给出。形式化要点：`q` 是 `μ̄` 的泛函（由 L1.2），`m` 由不变关联定义。

**L3.2（`q` 的 scheme 无关性）** `【本文 + 交接报告】`
- 对**耦合 scheme 变换** `g'² = g² + c g⁴`：`m/μ` 不变（`m` 是物理质量，`μ` 是物理标度）。`【交接报告的论证成立】`
- 对**blocking 方案**：有效质量 `m(μ)` 有 `O((a/ξ)²)` 的修正，故 `q` 只在连续极限严格不变。`【本文的精确化；交接报告未区分】`
建议：陈述里用**物理质量**定义 `q`，把有效质量当作计算代理。`【本文】`

**L3.3（单调性）** 若 `m(μ)` 满足 `μ dm/dμ = -γ_m m`（`γ_m>0`），则

```
dq/dμ = (μ m' - m)/μ² = -(γ_m + 1) m/μ² < 0,
```

故 `q` 随 `μ` 严格递减；`q=1 ⟺ μ=m`。弱侧 `μ≫m ⟹ q≪1`；强侧 `μ≲m ⟹ q≳1`。`【本文】`

**L3.4（交接引理）** 见 `research/physics_handoff.md` §3：若 `q₀:=q(μ_up) ≥ q₁`，则弱域 `[μ_up,∞)` 与强域 `(0,μ*]` 覆盖 `(0,∞)`，故 `∀μ>0, gap(μ) ≥ min(δ_W,δ_S)`。`【交接报告；纯实分析】`

**L4.1（有限格点规范固定独立性）** 已在 §3.5 给出。形式化要点：轨道分解 + `F` 在轨道上常数 + `μ` 在轨道上均匀。这是**最小的可验证 Lean 目标**（见 §8.1 N2）。

**L4.3（桥）** 由 L4.1 + L3.1 直接得到：`q` 是不变泛函，故任意规范固定给出相同值。

**L5.2（MG 陈述）** `【定义】`
```
MG(m_*) :  m_* > 0  ∧
  ∀ (Λ : 体积), ∀ (O₁ O₂ : 局部规范不变算子),
    ∃ C_{O₁O₂} > 0,   （C 不依赖 Λ）
    ∀ x y ∈ Λ,  |⟨O₁(x) O₂(y)⟩_{c,Λ}| ≤ C_{O₁O₂} e^{-m_* |x-y|}.
```
"体积一致" = `C` 与 `m_*` 不依赖 `Λ`；"标号一致" = `m_*` 不依赖 `O`。

**L5.3（MG ⟹ 质量间隙）** 由 L1.1 + L5.2；注意有限体积的 `T` 是有限维，`Δ_Λ ≥ m_*` 对所有 `Λ`，取极限得无限体积的 `Δ ≥ m_*`。

**L5.4（证明 MG）** 开放。路线见 §8。

---

## 5. 控制量的精确化

### 5.1 `q` 的定义

**定义 5.1** 尺度 `μ` 的有效理论是 `μ̄_μ`（轨道测度）。令

```
m(μ) := inf { 衰减率 of ⟨O(x)O(y)⟩_{c,μ̄_μ} : O 局部规范不变 }.
```

则 `m(μ)` 是尺度-`μ` 有效转移矩阵的谱隙（格点单位下 `= Δ_μ`）。定义 `q(μ) = m(μ)/μ`。

**引理 5.2** `q(μ)` 只依赖 `μ̄_μ`，故规范不变。`【由 §2.5、L1.2】`

### 5.2 与 `physics_handoff.md` 的对接

- 交接报告的 `q = m(μ)/μ` 就是本文的 `q`；`m(μ) = 1/ξ(μ)` 是有效关联长度倒数。
- 本文的精确化：
  - `q` 的**规范不变性**由轨道测度给出（报告未讨论规范）。
  - `q` 的**scheme 无关性**：对耦合 scheme 变换成立；对 blocking 有 `O((a/ξ)²)` 修正。建议陈述用物理质量。
  - `q` 的**单调性**与 `q=1` 处的穿越与报告一致。
- 交接报告的 `q₀ ≥ q₁` 判据在本文框架下是**规范不变的不等式**；两侧的 `gap` 都是不变子空间的谱隙。

### 5.3 三个陷阱

**(i) 局部算子 vs Wilson 环。**
质量间隙是**局部**规范不变算子的关联衰减率；Wilson 环 `W_C` 是非局部的，其大环行为由弦张力/面积律 `e^{-σ Area}` 支配，不是 `e^{-m Perimeter}`。用 Wilson 环定义 `m` 会得到错误的对象。正确：用 `tr F_{μν}F^{μν}` 的 smear（小环展开的首项），或一般的局部复合算子。
`【本文强调；交接报告与 M1 计划可能未区分】`

**(ii) 拓扑荷。**
`q_top(x) = tr F_{μν}F̃^{μν}` 的关联在纯 YM（θ=0）预期指数衰减（`0^{-+}` glueball），但其积分（拓扑磁化率）非零；这不矛盾。若理论有凝聚/序参量（纯 YM 没有），连通关联要减去真空期望。陈述用**连通**关联即可。`【标准】`

**(iii) 跑动质量 vs 物理质量。**
`m(μ)`（有效质量）是 scheme/blocking 依赖的计算量；`m_phys`（物理质量）是 `μ→0` 的极限/谱隙。`q` 应用物理质量定义才严格不变；有效质量是代理。`【本文】`

### 5.4 交接判据的规范不变表述

```
弱侧:  g²(μ) ≤ g*²  ⟺  μ ≥ μ_up,      gap(μ) ≥ δ_W
强侧:  q(μ) ≥ q₁    ⟺  μ ≤ μ*,        gap(μ) ≥ δ_S
交接:  q₀ := q(μ_up) ≥ q₁   ⟺   μ_up ≤ μ*
结论:  ∀μ>0, gap(μ) ≥ min(δ_W, δ_S) > 0
```

所有量都是规范不变的；`g²` 只是弱侧的坐标（scheme 依赖），不进入最终不等式。`【交接报告 + 本文】`

**循环性说明**：`q(μ)` 含未知的 `m(μ)`，故 `q` 是坐标而非自由参数；`q₀ ≥ q₁` 的非循环内容是两侧**估计**给出的 `q` 的界，而不是 `q` 的定义。`【交接报告】`

---

## 6. `Basic.lean` 手术方案

### 6.1 诊断（重申）

| 现状 | 问题 | 修复 |
|---|---|---|
| `FieldConfiguration d = 𝒮′(ℝᵈ)` | 标量场，不是联络模规范 | 引入 `G`、链接、`𝒜/𝒢`（或格点优先） |
| `NonTrivial = ¬Gaussian` | β=0 乘积测度非 Gaussian 但平凡 | 用"存在非零连通 n≥3 矩"（`ConnectedFourPoint` 已在文件里） |
| `ExponentialClustering` 用任意源 | 规范理论中 `j·A` 非不变 | 源 = 不变局部算子（或环） |
| `SatisfiesAllOS` | 标量 OS 谓词 | 规范不变 OS（OS 1978 型） |
| 质量间隙用生成泛函 | Wilson 环 ≠ 局部算子 | 局部不变算子的连通关联 |

### 6.2 建议结构（Lean 草图，示意，未编译）

```lean
-- 示意：不可编译；依赖 MeasureTheory、紧群 Haar 测度、有限格点组合。
namespace YMGauge

variable (Λ : Box 4) (G : Type*) [Group G] [TopologicalGroup G]
  [CompactSpace G] [MeasurableSpace G] [BorelSpace G]

abbrev Link := Edge Λ → G
abbrev Site := Vertex Λ → G

def gaugeAct (g : Site Λ G) (U : Link Λ G) : Link Λ G :=
  fun e => g (src e) * U e * (g (tgt e))⁻¹

def GaugeInvariant (μ : Measure (Link Λ G)) : Prop :=
  ∀ g, Measure.map (gaugeAct g) μ = μ

def IsInvariantObservable (F : Link Λ G → ℝ) : Prop :=
  ∀ g U, F (gaugeAct g U) = F U

-- L0.2（有限群/紧群）
theorem wilson_gaugeInvariant
    (β : ℝ) : GaugeInvariant (wilsonMeasure Λ G β) := by
  sorry

-- L2.1 + L2.2：等变块自旋 ⟹ pushforward 不变
structure EquivariantBlockSpin
    (Φ : Link Λ G → Link Λ' G) (b : Site Λ G →* Site Λ' G) : Prop where
  equiv : ∀ g U, Φ (gaugeAct g U) = gaugeAct (b g) (Φ U)

theorem pushforward_invariant
    {Φ : Link Λ G → Link Λ' G} {b : Site Λ G →* Site Λ' G}
    (hΦ : EquivariantBlockSpin Φ b) (hb : Function.Surjective b)
    (hμ : GaugeInvariant μ) :
    GaugeInvariant (Measure.map Φ μ) := by
  intro h
  obtain ⟨g, rfl⟩ := hb h
  calc Measure.map (gaugeAct (b g)) (Measure.map Φ μ)
      = Measure.map (fun U => gaugeAct (b g) (Φ U)) μ := by simp
    _ = Measure.map (fun U => Φ (gaugeAct g U)) μ := by
          simp [hΦ.equiv]
    _ = Measure.map Φ (Measure.map (gaugeAct g) μ) := by simp
    _ = Measure.map Φ μ := by rw [hμ g]

-- L1.1 / L5.2：规范不变的质量间隙
def ConnectedCorrelator
    (O₁ O₂ : Link Λ G → ℝ) (μ : Measure (Link Λ G)) (x y : Vertex Λ) : ℝ :=
  ∫ U, O₁ (translate x U) * O₂ (translate y U) ∂μ
    - (∫ U, O₁ (translate x U) ∂μ) * (∫ U, O₂ (translate y U) ∂μ)

def UniformExpClustering (μ : Measure (Link Λ G)) (m : ℝ) : Prop :=
  m > 0 ∧ ∀ (O₁ O₂ : Link Λ G → ℝ),
    IsInvariantObservable O₁ → IsInvariantObservable O₂ →
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x y : Vertex Λ,
      |ConnectedCorrelator O₁ O₂ μ x y| ≤ C * Real.exp (-m * latticeDist x y)

structure GaugeMassGapStatement (μ : Measure (Link Λ G)) : Prop where
  os : GaugeInvariantOS μ
  nonTrivial : ∃ n ≥ 3, ∃ O : Fin n → (Link Λ G → ℝ),
    (∀ i, IsInvariantObservable (O i)) ∧ cumulant μ O ≠ 0
  gap : ∃ m : ℝ, UniformExpClustering μ m

-- L4.1：有限格点规范固定独立性
-- 关键辅助引理：不变函数在轨道上常数
theorem invariant_const_on_orbit (F : Link Λ G → ℝ)
    (hF : IsInvariantObservable F) {U V : Link Λ G}
    (hUV : OrbitRel U V) : F U = F V := by
  sorry

theorem gaugeFixing_independence [Fintype G] [DecidableEq G]
    (μ : Measure (Link Λ G)) (hμ : GaugeInvariant μ)
    (F : Link Λ G → ℝ) (hF : IsInvariantObservable F)
    (s₁ s₂ : OrbitSection Λ G) :
    ∑ ω : Quotient, μ (orbitOf ω) * F (s₁ ω)
      = ∑ ω : Quotient, μ (orbitOf ω) * F (s₂ ω) := by
  apply Finset.sum_congr rfl
  intro ω _
  rw [invariant_const_on_orbit F hF
        (OrbitSection.same_orbit s₁ s₂ ω)]

end YMGauge
```

说明：
- 这是**示意**，不是可编译代码；`Box`、`Edge`、`orbitOf`、`OrbitSection`、`GaugeInvariantOS`、`cumulant` 等需要定义。
- 最短路径是**有限群 + 有限格点**：`L0.2`、`L2.1`、`L2.2`、`L4.1` 都是有限组合/有限和，mathlib 可支撑。
- `GaugeInvariantOS` 是重活（反射正性），可先作为 `structure` 字段，后续再证。

### 6.3 `NonTrivial` 的正确版本

**定义 6.1（非平凡）**
```
NonTrivial(μ) := ∃ n ≥ 3, ∃ O₁,…,O_n 局部规范不变,  cumulant_μ(O₁,…,O_n) ≠ 0.
```

- 排除 Gaussian（所有 `n≥3` 连通矩为零）。
- 排除 β=0 乘积测度（所有 `n≥2` 连通矩为零）。
- 包含交互 YM（`⟨(tr F²)⁴⟩_c ≠ 0` 预期成立）。
- 是规范不变的（算子在不变代数里）。
`【本文】`

文件里已有 `ConnectedFourPoint`，直接用它（取 `n=4`）即可作为第一版；更稳的是 `n≥3` 的 cumulant。

### 6.4 两条路线

**路线 A（格点优先，推荐形式化）**
在有限格点上定义 `𝒜_Λ, 𝒢_Λ, μ_β, 𝔄_inv, OS, gap`；证明 L0.2、L2.1、L2.2、L4.1；把 `GaugeMassGapStatement` 作为**格点陈述**。连续统 Clay 陈述作为"标度极限 + OS/不变性过极限"的定理（L0.6、L2.2 的极限版）。
优点：有限维、无解析硬点、可编译；与 2026 构造的格点起源一致。
缺点：不是直接的连续统陈述；需要极限定理。

**路线 B（连续统抽象）**
定义"轨道空间上的 OS 数据"结构（测度 + 不变代数 + 反射正性 + 聚类），把 `OSforGFF` 的自由场作为**平凡性测试**（自由 Maxwell 满足 OS 但**不**满足质量间隙——光子无质量；有质量 Proca 不是规范不变，不在范畴内）。
优点：贴近 Clay 陈述。
缺点：抽象、需要大量测度论基础设施；容易变成空洞陈述（参见 `formal-conjectures` 的教训）。

**建议**：A 为主，B 作为接口；`OSforGFF` 只用于"自由场测试"，不直接承载 YM 陈述。

---

## 7. 障碍清单

### 7.1 总表

| # | 障碍 | 性质 | 4D 特有？ | 可能绕法 | 严重性 |
|---|---|---|---|---|---|
| O1 | Balaban 块自旋的定量类保持 | 解析 | 是（4D SU(N)） | 多尺度 + cluster expansion；或随机量子化 | 极高 |
| O2 | 弱→强交接 `q₀≥q₁` | 解析 | 是 | 两侧估计 + 不变量判据 | 极高 |
| O3 | 4D 无穷 counterterm | 解析 | 是 | 渐近自由 + Ward/BRST；regularity structures 不够 | 极高 |
| O4 | Gribov / 无全局截面 | 拓扑 | 否（所有维） | 轨道空间表述；局部规范；带号 FP | 中（UV 可绕） |
| O5 | 商空间奇异/非流形 | 几何 | 否 | 测度论表述；Polish 轨道空间（CCHS 2D） | 中 |
| O6 | OS 重构在规范理论的形式化 | 形式化 | 否 | OS 1978 + 有限维转移矩阵 | 中 |
| O7 | 质量间隙的算子选择（局部 vs 环；拓扑荷） | 陈述 | 否 | 局部不变算子 + 连通关联 | 低（但必须写对） |
| O8 | 陈述层 `Basic.lean` 是标量 | 工程 | 否 | §6 手术 | 低（但阻塞 M1） |
| O9 | 标度极限的规范不变性未审计 | 审计 | 是 | 检查 16 输入 / Lean 陈述 | 中（可能隐藏假设） |
| O10 | 文献中 AI 生成的"YM 证明"噪声 | 文献 | — | 只信同行评审/可复现 | 低（但耗时） |

### 7.2 前五条的说明

**O1（Balaban 定量类保持）**：这是本报告框架下**唯一**真正的规范相关硬点——但注意，它难在**正则性与类保持**，不难在等变性。等变性由 L2.1 精确给出。Balaban 的系列（CMP 95/96/98/99/102/109/116/119/122）是唯一完整的候选；其细节技术性极强，未被独立形式化。`【开放/待核实】`

**O2（交接）**：见 `physics_handoff.md`。本文只保证 `q` 是规范不变的坐标；`q₀≥q₁` 的估计不在本文范围。

**O3（counterterm）**：4D 非超可重整；无穷多 counterterm 必须由 Ward/BRST 组织。规范协变性在这里是**约束**（counterterm 必须规范不变/协变），不是障碍；障碍是重整化本身。CCHS 3D 的"唯一 gauge-covariant 重整化"是范例，4D 无对应。`【开放】`

**O4（Gribov）**：**UV RG 不受影响**（局部小场规范无拷贝）；IR/大场才受影响。不要把 Gribov 当作 UV RG 的障碍。`【本文】`

**O5（商奇异）**：有限格点无问题；连续统需要 CCHS 式的 Polish 轨道空间（2D 已做，4D 未做）。对格点优先的陈述层无影响。

---

## 8. 下一步

### 8.1 最小可验证工件（按优先级）

**N1（本周，文档 + Lean 草图）—— 规范不变陈述层**
写 `GaugeStatement.lean`（新文件，或在 `YMStatement/` 下加 `Gauge.lean`）：有限格点数据、`gaugeAct`、`GaugeInvariant`、`IsInvariantObservable`、`GaugeMassGapStatement`（§6.2 草图）。先证 L0.2（有限群版）。
**验收**：`lake build` 通过；`GaugeMassGapStatement` 的每个字段非空洞（用 β=0 测度测试：满足 `os` 与 `gap`，但**不**满足 `nonTrivial`；用自由 Maxwell 测试：满足 `os`，不满足 `gap`）。
**价值**：修复 `Basic.lean` 的标量错位；这是 M1 的真正前提。

**N2（1–2 周，Lean）—— 有限格点规范固定独立性**
`gaugeFixing_independence`（§6.2 末）：有限群、有限格点、任意截面、不变 `F`。
**验收**：0 sorry；一个具体例子（如 SU(2) 单 plaquette 或小环面）数值验证。
**价值**：这是"把规范固定构造接到不变可观测量"的**第一个可机器检验定理**；也是本子项目唯一完全自足的新 Lean 工件。

**N3（并行，文档）—— 等变 RG 的结构引理**
把 §2.6 的 L2.1–L2.4 写成精确命题（含 `b` 满射、`Avg` 共轭等变、有限群显式公式）；作为 N1 的一部分或独立 note。
**价值**：把"规范不变 RG"从口号变成可引用的引理；主 session 的交接工作可以直接用。

**N4（中期，文献/审计）—— Balaban 的精确陈述**
读 Balaban CMP 98（平均操作）、CMP 99（规范固定条件）、CMP 109（有效作用量），提取：
- 块自旋的精确定义与正则类；
- 有效作用量的类与误差界；
- 是否有显式的等变缺陷界。
**注意**：Dimock I–III 是 φ⁴₃，不是规范理论；只能借方法。
**价值**：O1 的现状评估；决定是否值得投入形式化。

**N5（接口）—— 与主 session 的交接**
把 §5 的 `q` 定义、不变性、单调性交给主 session；明确 `q₀≥q₁` 的不变性论证已就绪，不需要重复推导。**不重复**交接报告的解析估计。

### 8.2 审计清单（针对 2026 4D 构造 / 16 输入）

| 问题 | 为什么重要 |
|---|---|
| 构造是否在**链接**上定义有限体积测度并证明规范不变？ | 若否，整个 OS 陈述可能建立在规范切片上 |
| OS 公理是否陈述在**规范不变代数**上？ | 若是标量代理，mass gap 陈述不对应物理 |
| 标度极限的**规范不变性**是证明的还是假设的？ | 可能是隐藏输入 |
| `16` 个输入中哪个承载"规范不变性"？ | 审计的可定位性 |
| `MassGap` 是否出现在任何主定理中？ | 委员会说没有；确认陈述层的缺口 |
| 非平凡性如何陈述？ | 若是 `¬Gaussian`，需按 §6.3 修正 |

**建议**：这是与 `research/lean-infra` 或审计 session 的接口；本报告不重复代码审计，只给清单。

### 8.3 明确不做的事

- 不证明 MG，不改进 `q₀≥q₁` 的常数。
- 不重建 Balaban 的估计。
- 不形式化 OS 重构的完整测度论（超出本子项目）。
- 不把 Zenodo/vixra 上的"YM 证明/GZ 量子化"预印本当作依据（AI 参与度高、venue 可疑；本次检索中已多次出现）。

---

## 9. 文献（核实状态）

### 9.1 本次核实（2026-11-22，web 检索）

| 文献 | 出处 | 备注 |
|---|---|---|
| Bałaban, *Averaging operations for lattice gauge theories* | CMP **98** (1985) 17–51, DOI 10.1007/BF01211042 | 摘要明确：链接与规范变换的平均、正则类 |
| Bałaban, *Spaces of regular gauge field configurations…gauge fixing conditions* | CMP **99** (1985) 75–102 | 局部规范固定 |
| Bałaban, *RG approach … I. Generation of effective actions…* | CMP **109** (1987) 249–301, DOI 10.1007/BF01215223 | 小场有效作用量、耦合重整化 |
| Bałaban, *… II. Cluster expansions* | CMP **116** (1988) 1–22, DOI 10.1007/BF01239022 | |
| Bałaban, *UV stability of 3D lattice pure gauge…* | CMP **102** (1985) 255–275 | |
| Bałaban, *Large field renormalization I/II* | CMP **122** (1989) 175–202 / 355–392 | |
| Bałaban, *Propagators and renormalization transformations I/II* | CMP **95** (1984) 17–40; **96** (1984) 223–250 | |
| Bałaban, *Convergent renormalization expansions…* | CMP **119** (1988) 243–285 | |
| Dimock, *RG according to Balaban I. Small fields* | Rev. Math. Phys. **25** (2013) 1330010, arXiv:1108.1335 | **φ⁴₃，不是规范** |
| Dimock, *… II. Large fields* | J. Math. Phys. **54** (2013) 092301, arXiv:1212.5562 | 同上 |
| Dimock, *… III. Convergence* | Ann. Henri Poincaré **15** (2014) 2133–2175, arXiv:1304.0705 | 同上 |
| Osterwalder–Seiler, *Gauge field theories on a lattice* | Ann. Physics **110** (1978) 440–471, DOI 10.1016/0003-4916(78)90039-8 | 反射正性、强耦合聚类、质量间隙 |
| Chandra–Chevyrev–Hairer–Shen, *Langevin dynamic for the 2D YM measure* | Publ. Math. IHÉS **136** (2022) 1–147, arXiv:2006.04987 | 轨道空间 Polish；DeTurck；gauge covariance in law |
| Chandra–Chevyrev–Hairer–Shen, *Stochastic quantisation of YMH in 3D* | Invent. Math. **237** (2024) 541–696, arXiv:2201.03487 | 唯一 gauge-covariant 重整化 |
| Chevyrev–Shen, *Uniqueness of Gauge Covariant Renormalisation…* | ARMA **250**(1) (2026) | 标题已核实；arXiv 号未核实 |
| Shen–Zhu–Zhu, *A stochastic analysis approach to lattice YM at strong coupling* | CMP **400** (2023) 805–851, arXiv:2204.12737 | LSI/Poincaré、唯一测度、质量间隙 |
| Singer, *Some remarks on the Gribov ambiguity* | CMP **60** (1978) 7–12, DOI 10.1007/BF01609471 | 无全局连续截面 |
| Gribov (1977) | SLAC-Trans-0176 | Coulomb 规范歧义 |
| Menotti–Pelissetto, *OS positivity for the Wilson action* | Nucl. Phys. B (Proc. Suppl.) **4** (1988) 644 | 过格点反射面 |
| Bauerschmidt–Bodineau–Dagallier, *Stochastic dynamics and the Polchinski equation* | arXiv:2307.07619, Probab. Surveys 21 (2024) | 多尺度 Bakry–Émery |
| Bauerschmidt–Dagallier, *LSI for φ⁴₂/φ⁴₃* | arXiv:2210.13487 | |

### 9.2 项目文档来源（未在本次独立核实）

- M. R. Douglas, *The Yang–Mills Millennium Prize Problem* (Stony Brook 2026); Nature Rev. Phys. 8 (2026) 86–97.
- hexagonmath, *From 4D Wilson Lattice Gauge Theory to Quantum Yang–Mills Theory*, 2611.00492v1 (2026-11-15).
- Yang–Mills Working Group, *Validation of the OpenAI YM construction*, 2026-11-10.
- M. R. Douglas, S. Hoback, A. Mei, N. Nissim, *Formalization of QFT*, arXiv:2603.15770; `github.com/mrdouglasny/OSforGFF`.
- Magnen–Rivasseau–Sénéor, CMP **155** (1993) 417–430.
- Seiler, LNP **159** (1982).
- Lüscher, arXiv:1006.4518; Cao–Park–Sheffield, arXiv:2307.06790; Cao–Nissim–Sheffield, arXiv:2505.16585 / 2509.04688.

### 9.3 未核实 / 需查证

- Giles, *Reconstruction of gauge potentials from Wilson loops*, Phys. Rev. D **24** (1981) 2160（经典，本次未核实）。
- Bauerschmidt–Dagallier, *LSI for U(1) lattice gauge theory in d=4*（2024；标题/作者来自二手来源，**arXiv 号未核实**）。
- Neuberger 零点的精确出处（标准 folklore；本次未定位原始文献）。
- Balaban 各篇的精确定理编号与常数（需读原文）。

---

## 附录 A：记号

| 记号 | 含义 |
|---|---|
| `Λ, Λ'` | 细/粗格点（下文 `Λ := Λ_{L,a}`） |
| `T⁴_L`, `Λ_{L,a}` | 格点环面 `(ℤ/Nℤ)⁴`，`N=L/a`，周期边界，`θ=0`，平凡 't Hooft 通量 |
| `N_c` | `dim f`，定义表示维数（`SU(N)` 时 `=N`）；`β = 2N_c/g²` |
| `C` | A0 §1.4 的 RT scheme 常数（全文同一条 RT） |
| `𝒜_Λ = G^{E}` | 链接配置空间 |
| `𝒢_Λ = G^{V}` | 规范群 |
| `g·U` | 规范作用 `g_x U_{xy} g_y^{-1}` |
| `μ_β, μ̄` | Wilson 测度 / 轨道测度 |
| `π` | 商映射 `𝒜 → 𝒜/𝒢` |
| `𝔄_inv, 𝔄_inv^loc` | 不变代数 / 局部不变代数 |
| `W_C` | Wilson 环 |
| `Θ` | 时间反射 |
| `T, ℋ` | 转移矩阵 / 不变 Hilbert 空间 |
| `Δ, m, ξ` | 谱隙 / 质量 / 关联长度 |
| `Φ, b` | 块自旋 / 块规范同态 |
| `δ` | 等变缺陷 |
| `q = m/μ` | 无量纲控制量 |
| `Δ_FP` | Faddeev–Popov 行列式（带号） |

## 附录 B：一页纸

```
规范不变性 = 三件可分离的事：
  (a) 测度不变        已证（Haar）
  (b) RG 等变         结构初等（面平均 + 平行移动）；定量是 Balaban 的开放内容
  (c) 控制量不变      由 (a)(b) 推出；q, Δ, ξ 定义在不变代数上

框架：
  轨道测度 μ̄  +  等变块自旋 (Φ,b)  +  不变控制泛函 q
  Gribov 被绕开（不需要全局截面）；局部规范固定仅作计算工具

桥：
  不变可观测量的规范固定独立性（有限格点精确；紧群带号 FP）
  ⟹ q 可在任意规范固定下计算，结果相同

陈述层：
  Basic.lean 是标量场，不是 YM；NonTrivial 应为"非零连通 n≥3 矩"
  质量间隙用局部不变算子的连通关联，不用 Wilson 环

真正的开放点（与规范无关）：
  O1 Balaban 定量类保持（4D SU(N)）
  O2 q₀ ≥ q₁ 交接
  O3 4D 无穷 counterterm

下一步：
  N1 规范不变陈述层（Lean）
  N2 有限格点规范固定独立性（Lean，最小可验证）
  N3 等变 RG 结构引理（文档）
  N4 Balaban 精确陈述（审计）
```

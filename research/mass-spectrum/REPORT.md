# 孤立单粒子态（上限隙）与低能胶球谱：定理、引理链、障碍

> 2026-11-22 · 范围：4D 纯 Yang–Mills（G = SU(N), N ≥ 2）。存在性条件化于 16 输入、仅子列极限；质量隙（MG）等价于体积一致、标号一致的指数聚类速率。
> 本报告**不重复**主 session 的弱→强耦合交接控制量分析（`research/physics_handoff.md`），也**不重复**禁闭/面积律（`research/confinement/`，那里已说明面积律管的是**带电弦扇区** `E₀(L) ≥ σL`）。本报告处理**中性胶球扇区**：`Δ = inf spec(H)|_{Ω⊥}` 是否**被一个孤立原子达到**。
> 结论先行：**"有质量隙"与"存在孤立单粒子态"是两个不同的谱命题**。Clay 只要前者；后者严格更强，且 4D YM 中未证。可证的核心是二者之间的精确等价（定理 A）与不等价（定理 B），加上一个约 150 行的抽象定理；真正的墙是**上限隙**，其 O(3) 式证明依赖一个 YM 不具备的 Z₂。

---

## 0. 结论摘要（≤300 字）

质量隙 = `inf spec(H)\{0} > 0`；孤立单粒子态 = 这个下确界是谱测度的一个**孤立原子**。二者不等价：正测度可在正阈值 `Δ` 处有连续密度而没有原子，此时两点函数仍按 `e^{−Δt}` 衰减（带 `1/t` 前因子），质量隙成立却无粒子。**所以 MG 不蕴含胶球存在，Clay 问题也不要求它。**

可证的核心是 KL 表示在时间方向的精确形式：**孤立原子 + 上限隙 ⇔ 两点函数是"尖锐指数衰减" `Z e^{−m₁t} + O(e^{−m₂t})`（定理 A）**。这可在 mathlib 里完整证明（无需谱测度、无需无界算子谱定理）。

真正的缺口是**上限隙**（`(m₁,m₂)` 无谱）：O(3) 用自旋翻转 Z₂ 得到"偶扇区 ≥ 2e"，把上限隙变成对称性事实；SU(N) 的 `0^{++}` 通道没有这样的 Z₂（C 与 P 都不能把单胶球和双胶球分开），该步不能移植。C 奇/P 奇通道是一个可用的类比点，但仍需粒子诠释。

---

## 1. 精确定理

### 1.1 记号与 Källén–Lehmann 表示

**A0 测度来源（GL1 合规）。** 格点层取 `research/gluing/REPORT.md`（A0）§1.2 与 `research/glue-c-gauge-group/REPORT.md` (1.1)–(1.2) 的**共享族**，不另立变体：

```
S_W(β) = β Σ_p (1 − (1/N_c) Re tr U_p),      β = 2N_c/g²,     N_c := dim f,
μ_{L,a,β} ∝ exp(−S_W(β)) ∏_e dU_e  =  exp((β/N_c) Σ_p Re tr U_p) ∏_e dU_e ,      (A0)
```

这是**概率测度**（`Z^{-1}` 归一）。进入测度的唯一组合是 `γ := β/N_c = 2/g²`。

> **写法警告。** 若把权重写成裸 `exp(β' Σ_p Re tr U_p)`，则必须 `β' = γ = 2/g²`（`SU(3), g²=1 → β' = 2`），**不是** `β = 2N_c/g² = 6`。二者差因子 `N_c`（`glue-c` 在 `SU(3), g²=1` 处给出 `6` vs `2` 的具体反例）。本报告不写裸形式。

本报告的 `μ` 是 `μ_{L,a,β}` 在**先 UV 后 IR**（固定物理体积 `L` 先取 `a→0`，再沿子列 `L→∞`；A0 §0.2 第 1 条，子列性见 `docs/yang-mills-prior-work.md` §3.3）下的极限理论。§1 的全部谱陈述都在该理论里陈述；`μ` 不唯一时，谱测度按每个极限理论分别陈述（与 `uniqueness` 的接口）。

**记号表（避免双关）。**

| 符号 | 含义 |
|---|---|
| `μ` | OS 概率测度，即 (A0) 的极限。**只**表示测度 |
| `m`, `m_n`, `m_ρ^{(n)}`, `m₁, m₂` | 质量 / 谱质量 / 通道质量。**不用 `μ` 表示质量** |
| `Δ` | 质量隙 |
| `Δ₊(x; m)` | 质量 `m` 的自由（Euclidean）传播子 |
| `s` | RG 尺度（`physics_handoff` / `confinement` 两报告以另一字母兼作尺度，本报告统一记 `s`） |
| `μν`（下标） | Lorentz 指标，与测度无关 |
| `ν_O` | `O` 的两点谱测度（§1.1 定义） |

设 `μ` 为 Euclidean OS 测度（OSforGFF 的 `SatisfiesAllOS`），`(ℋ, H, Ω)` 为 OS 重构给出的物理 Hilbert 空间、Hamiltonian、唯一真空：`H ≥ 0`，`HΩ = 0`，`ker H = ℂΩ`。固定一个**规范不变**的 smeared 算符 `O`（有界；UV 正规化下成立。连续统中取 OS1 型界保证 smear 后有界），令

```
A := O − ω(O)·1,        ω(O) := ⟨Ω|O|Ω⟩,        A(t) := e^{itH} A e^{−itH}.
```

**谱测度。** `ν_O(B) := ⟨AΩ | P_H(B) AΩ⟩`，`B ∈ 𝔅([0,∞))`。这是 `[0,∞)` 上的**有限正测度**，`ν_O([0,∞)) = ‖AΩ‖²`。两点函数是它的 Fourier–Laplace 变换：

```
W_O(t) := ω(A(t) A(0)) = ∫ e^{−iEt} dν_O(E).                          (1.1)
```

**这就是时间方向的 Källén–Lehmann 表示**，`ν_O` 即 KL 谱测度。它是"谱定理作用于向量 `AΩ`"的全部内容——**不需要**无界自伴算子的谱定理，只需要 `H` 的谱测度在单个向量上的值。

**时空形式（连续统）。** 对 Wightman 场，`W_O(x−y) = ∫₀^∞ dρ_O(m²) Δ₊(x−y; m)`，`Δ₊` 为质量 `m` 的自由传播子，`ρ_O ≥ 0` 为正测度（Usui 2012 在格点上由 Hermiticity + 平移不变 + 反射正性 + 多项式有界证明了这个表示与正性；`arXiv:1201.3415`）。**零空间动量**下，时间方向核为 `e^{−m|t|}/(2m)`，于是

```
C_O(t) := W_O^{p=0}(t) = ∫₀^∞ e^{−m|t|} · dρ_O(m²)/(2m),             (1.2)
```

即 `ν_O(dm) = dρ_O(m²)/(2m)`。**格点上 (1.2) 是精确的有限和**（见 §1.5），连续统中是极限。下文一律用质量变量 `m`；两点函数的能量变量记 `E`。

### 1.2 三个定义

**定义 1（质量隙）。** `MassGap(μ) :⇔ ∃ Δ > 0, spec(H) ∩ (0,Δ) = ∅`（`Basic.lean` 的 `MassGapStatement` 的谱形式；`Basic.lean` 写的是生成泛函的一致指数聚类，二者关系见 §2.2）。

**定义 2（孤立单粒子态 / 上限隙）。** 对通道 `O`：

```
IsolatedOneParticle(O; m₁, m₂) :⇔
   0 < m₁ < m₂
   ∧ ν_O({m₁}) > 0                                  -- 原子（粒子）
   ∧ ν_O((0, m₂) \ {m₁}) = 0.                       -- 下方无谱 + 上方直到 m₂ 无谱（上限隙）
```

等价的谱表述：`inf spec(H |_{ℋ_O}) = m₁` 是特征值（`ℋ_O` 为 `O` 生成的循环子空间；重数可 >1），且 `spec(H |_{ℋ_O}) ∩ (m₁,m₂) = ∅`。**"上限隙"就是 `m₁` 在 `supp ν_O` 中孤立。** 注意这里没有假设 `m₁` 是最轻的总体质量——它只是该通道最轻的。

**定义 3（尖锐指数衰减）。**

```
SharpDecay(O; m₁, Z, m₂, C) :⇔
   0 < Z ∧ m₁ < m₂ ∧
   ∀ t ≥ 0, |W_O(−it) − Z e^{−m₁t}| ≤ C e^{−m₂t}.
```

### 1.3 定理 A（KL 原子 ⇔ 尖锐指数衰减）——**可证**

> **定理 A.** 设 `ν` 是 `[0,∞)` 上的有限正测度，`0 < m₁ < m₂`。则
> ```
> IsolatedAtom(ν; m₁,m₂)  ⇔  ∃ Z C, SharpDecay(ν; m₁,Z,m₂,C).
> ```

**证明（⇐ 是全部内容；⇒ 显然）。** 假设 `|W(t) − Z e^{−m₁t}| ≤ C e^{−m₂t}`，其中 `W(t) = ∫ e^{−Et}dν`，且 `t ≥ 0`。

1. **下方无质量。** 若 `ν([0, m₁−ε]) > 0`，则 `W(t) ≥ ν([0,m₁−ε]) e^{−(m₁−ε)t} ≫ e^{−m₁t}`，故 `W(t)e^{m₁t} → ∞`，与假设 `W(t)e^{m₁t} = Z + O(e^{−(m₂−m₁)t}) → Z` 矛盾。故 `ν([0,m₁)) = 0`。
2. **原子质量恰为 Z。** 现在 `∫ e^{−(E−m₁)t} dν = W(t)e^{m₁t} → Z`。被积函数在 `[m₁,∞)` 上 `≤ 1` 且逐点收敛到 `1_{E=m₁}`，由控制收敛（`ν` 有限）得 `ν({m₁}) = Z > 0`。
3. **`(m₁,m₂)` 无质量。** 记 `ν' := ν − Zδ_{m₁} ≥ 0`。对任意 `E₁ < m₂`，
   ```
   ν'([0,E₁]) · e^{−E₁t} ≤ ∫ e^{−Et} dν' = |W(t) − Z e^{−m₁t}| ≤ C e^{−m₂t},
   ```
   故 `ν'([0,E₁]) ≤ C e^{−(m₂−E₁)t} → 0`。于是 `ν'([0,m₂)) = 0`。∎

**为什么这重要。** 它把 KL 表示（"原子"）、谱理论（"孤立特征值"）、转移矩阵（"单指数 + 指数小余项"）三层**精确**接起来，而且**只用测度论**。这是本子问题真正的"陈述层"，也是唯一可以现在就在 Lean 里 0 sorry 做完的部分。注意方向：`⇐` 需要**正性**（`ν ≥ 0`，即 `ν` 真是某个谱测度）——没有正性，符号测度的 Laplace 变换可以任意快衰减，结论假。

**推论 A′（质量隙 = 最轻胶球质量，精确含义）。**
```
MassGap ∧ (∃ m₁, ν_O({m₁}) > 0 ∧ supp ν_O ⊂ {0} ∪ [m₁,∞) ∧ inf supp(ν_O\{0}) = m₁)
  ⇔  Δ = m₁  且  m₁ 是孤立特征值。
```
即"质量隙等于最轻胶球质量"**本身**就是原子命题，而不是定义。若底部不是原子，`Δ` 只是下确界，没有"质量"可言。

### 1.4 定理 B（质量隙 **不**蕴含孤立粒子）——**可证，显式反例**

> **定理 B.** 存在 `[0,∞)` 上的有限正测度 `ν`（某个自伴算子的谱测度）与 `Δ > 0`，使得
> (i) `supp ν ⊂ [Δ,∞)` 且 `inf supp ν = Δ`（即 `MassGap` 成立，`Δ` 为隙）；
> (ii) `ν` **没有原子**（故 `¬ IsolatedOneParticle`）；
> (iii) `W(t) ≤ (∫dν) e^{−Δt}`，且 `−(1/t) log W(t) → Δ`（指数衰减率**恰为** `Δ`）。

**构造。** 取 `ν := ρ(E)dE`，`ρ` 在 `[Δ,∞)` 连续、`ρ(Δ) > 0`、`∫ρ < ∞`（例如 `ρ(E) = e^{−(E−Δ)}`）。则 `ν` 无原子；`W(t) = ∫_Δ^∞ e^{−Et}ρ(E)dE ≤ e^{−Δt}∫ρ`；且令 `E = Δ + s/t`，
```
W(t)e^{Δt} = ∫₀^∞ e^{−s} ρ(Δ+s/t) ds/t → ρ(Δ)/t,
```
故 `W(t) ~ ρ(Δ)e^{−Δt}/t`：衰减率精确为 `Δ`，**没有**原子。`ν` 可实现为 `L²([Δ,∞),ρ)` 上乘 `E` 算子的谱测度（取循环向量 `1`）。

**含义（本报告最重要的逻辑点）。**
- `ExponentialClustering`（`Basic.lean` 的定义）与 MG 都是"有速率"的陈述，**都允许** `1/t` 型前因子，因此都**不**给出原子。
- 反过来，`IsolatedOneParticle` ⇒ `MassGap`（显然），也 ⇒ 更强的 `SharpDecay`。
- 故严格包含：**`IsolatedOneParticle ⊊ MassGap`**，且 `MG` 落在中间（`IsolatedOneParticle ⇒ MG ⇒ MassGap`）。
- 这与 OpenAI 的 O(3) 论文开头一字不差："A mass gap and an isolated particle mass are different spectral statements. A positive spectral measure may start at a strictly positive threshold and still have no atoms."

**定理 B′（OS 层面的加强）。** 上面的 `ν` 可升级为一个**满足全部 OS 公理**的模型：取**广义自由场**（Gaussian），KL 密度为 `ρ(m²)dm²`，`ρ` 在 `[Δ,∞)` 连续、`ρ(Δ)>0`、`∫ρ(m²)dm²/(1+m²) < ∞`。其协方差 `C(x) = ∫dρ(m²)Δ₊(x;m)`，生成泛函 `Z[J] = exp(−½⟨J,CJ⟩)`：

- 反射正性由 `ρ ≥ 0` 保证（KL 构造），故满足 OS0–OS4（Gaussian 版的 OSforGFF 论证**逐字照搬**）；
- 两点函数按 `e^{−Δt}/t^{…}` 衰减 ⇒ `ExponentialClustering` 与 MG 成立，`Δ` 为质量隙；
- 谱测度无原子 ⇒ **没有孤立单粒子态**。

这是对 `MG ⇏ IsolatedOneParticle` 的**构造性反例**（Gaussian、非平凡性被 `NonTrivial` 排除，但逻辑反例只需放弃该条）。它说明：**"有质量隙"这个陈述本身，在 OS 公理内部，不蕴含任何粒子。**

**极/原子字典。** KL 测度的原子 `ν({m₁}) = Z` 就是传播子 `W(p)` 在 `p² = −m₁²` 处的**实极点**，留数 `Z`。阈值以上的复极点 = 共振；实极点 = 稳定粒子。故"存在孤立单粒子态" = "传播子有一个孤立实极点"，与"粒子"的物理语言完全一致。

### 1.5 定理 C（格点通道分解；`Δ = min_ρ m_ρ`）——**可证（有限维）**

设有限格点 `Λ`，时间切片上的 Hilbert 空间 `ℋ_Λ`（有限维），转移矩阵 `T = e^{−aH}` 为**对称正**算子。由 mathlib 的 `eigenvectorBasis`（`Mathlib/Analysis/InnerProductSpace/Spectrum.lean`，对有限维内积空间上的对称算子给出按特征值降序的幺正特征基），

```
T ψ_n = λ_n ψ_n,   1 = λ₀ > λ₁ ≥ λ₂ ≥ … ≥ λ_{N−1} ≥ 0.
```

定义**质量** `m_n := −a⁻¹ log λ_n`（`λ_n > 0`；强耦合下 `λ_n` 全正）。则

> **定理 C（格点谱分解）。**
> (a) `TransferGap :⇔ λ₁ < 1  ⇔  MassGap`，且 `Δ = m₁`。
> (b) 对时间零切片算符 `A, B`，`⟨Ω|A T^t B|Ω⟩ = Σ_n λ_n^t ⟨Ω|Aψ_n⟩⟨ψ_n|BΩ⟩` — **这就是格点上的 KL 表示**；取 `B = A` 时 `ν = Σ_n |⟨Ω|Aψ_n⟩|² δ_{m_n}` 是正原子测度。
> (c) 若 `T` 与有限对称群 `G_h = O_h × {1,C}`（晶格八面体群 × 电荷共轭）对易，则 `ℋ_Λ = ⊕_ρ ℋ_Λ^ρ` 按 `G_h` 的不可约表示分解，`T` 保持每个 `ℋ_Λ^ρ`。定义通道质量 `m_ρ^{(n)} := −a⁻¹ log λ^{(ρ)}_n`（`λ^{(ρ)}_n` 为 `T|_{ℋ_Λ^ρ}` 的特征值，降序）。令 `m_ρ^{(1)}` 记通道 `ρ` 中**最小正质量**（真空只在平凡通道，是 `m_平凡^{(0)} = 0`）。则
> ```
> Δ = min_ρ m_ρ^{(1)}.
> ```
> (d) **胶球谱定义**：`{ (ρ, n) ↦ m_ρ^{(n)} }`，即每个晶格不可约表示通道的转移矩阵特征值质量。这就是数值格点谱学所测的对象（§4.5）。

**(c) 的物理内容。** 由于物理 Hilbert 空间里每个态都是规范不变的，`spec(H)\{0}` 就是各通道非真空谱的并，故 `Δ` 是各通道最小正质量的**最小值**——这是定理，不是假设。**"最轻的是 `0^{++}`" 才是物理断言**（数值上 `1.71 GeV`，§4.5），不是定理。

**注意 (c) 的两个未形式化点**：(i) mathlib 的有限群表示论需要搭建"特征投影 `P_ρ = (d_ρ/|G|)Σ_g χ_ρ(g)U(g)` + 同型分量分解"；(ii) 连续统 `J^{PC}` 标注需要 `O_h` 表示在 `a→0` 下的 subduction，这是另一个未形式化的极限步骤（§7 O6）。

### 1.6 定理 D（4D YM 的孤立单粒子态；**目标，公开**）

> **定理 D（目标）。** 设 `μ` 为 4D 纯 `SU(N)` YM 的 OS 测度（16 输入 + MG 之下），`(ℋ,H,Ω)` 为其重构。固定 `0^{++}` 通道的规范不变算符 `O`（例如 smeared plaquette 减去 VEV，或 `tr F²` 的正规化版本），`ω(O) = 0`，`OΩ ≠ 0`。则存在 `m₁ > 0`、`Z > 0`、`δ > 0`：
> ```
> ν_O = Z δ_{m₁} + ν_rest,     supp ν_rest ⊂ [(1+δ)m₁, ∞),
> ```
> 且 `m₁ = Δ`（质量隙）。

**超出 MG 的额外输入（三者都未证）：**
- **(I1) 底部可达**：`inf spec(H)|_{Ω⊥}` 是特征值，且 `⟨Ω|O|n⟩ ≠ 0`（算符看见最低态）。
- **(I2) 上限隙**：`spec(H) ∩ (m₁,(1+δ)m₁) = ∅`。
- **(I3) 算符定义良好**：`O` 的重整化（`tr F²` 是 4 维复合算符，需减除与混合；§4.3、§7 O5）与规范不变性。

**定理 A 说明 (I1)+(I2) 恰好等价于 D 的结论**，所以 D = `(I1) ∧ (I2)`（在 MG 与 (I3) 之上）。本报告剩下的全部工作是说明 `(I1)`、`(I2)` 各自有多难（§5、§7）。

### 1.7 定理 E（非空洞性：自由场）——**可证，OSforGFF 已具备**

对 `d ≥ 2` 的有质量自由场 `GFF(d,m)`，OSforGFF 给出 `S₂(f,g) = ∫∫ f(x) C(x−y) g(y)`，`C = Cprofile(‖·‖)`，且

```
GFFPropagator.fourier_eq :  𝓕(Cprofile(‖·‖))(k) = 1/((2π)²‖k‖² + m²).
```

即动量空间**单个极点**，谱测度 `ν = δ_m`（单原子，`Z = 1`，`m₂ = ∞`）。故 `IsolatedOneParticle(GFF; m, ∞)` 成立，且 `SharpDecay` 对任意 `m₂ > m` 成立（余项为 0）。**定理 A、D 的定义在自由场上非空洞**；无质量极限 `m → 0` 失败（`ν = δ₀`，隙不成立）。

---

## 2. 与 Clay 质量隙的关系

### 2.1 严格强于 Clay 的隙

```
IsolatedOneParticle  ⇒  SharpDecay  ⇒  MG  ⇒  MassGap.
```

（在**谱测度层面**：`SharpDecay` 就是两点函数的 `e^{−m₁t}` 衰减，`MG` 是它的一致版本。`MG` 在 `Basic.lean` 里写成生成泛函形式，二者的桥见 §2.2。）

反向全部不成立：`MassGap ⇏ MG`（hexagonmath 2611.00492 的反例，见 `docs/yang-mills-prior-work.md` §1.4）；`MG ⇏ IsolatedOneParticle`（定理 B）。**Clay 问题只要最右端的 `MassGap`。** 所以：

> **本子问题是严格强于 Clay 质量隙的目标。** 证明它免费得到 Clay 的隙；但 Clay 的隙不需要它。它在获奖问题里**不在关键路径上**，不该阻塞主 session 的 MG 攻坚。

### 2.2 `Basic.lean` 的 `ExponentialClustering` 落在这条链的哪里

`Basic.lean` 的 `ExponentialClustering(μ,m)` 是**生成泛函** `Z[J]` 的一致指数聚类（`‖Z[f + T_a g] − Z[f]Z[g]‖ ≤ C e^{−m‖a‖}`）。它与谱测度的关系是：由 OS4/连通两点函数，`Z` 的聚类率控制的是**连通关联**的衰减，而 `ν_O` 的底部控制的是**同一**衰减率的**极限**。二者在"率"的层面一致，但在"是否原子"的层面无关——**定理 A/B 正是在这个界面上**。形式化时需要补一条：`GenerationFunctionalClustering ⇒ 两点函数 `C_O(t)` 的谱表示 (1.2)`，这是 OSforGFF 的 `gaussianFreeField_satisfies_OS4` 的抽象化（自由场已证，互作用场需 `SatisfiesAllOS` + 重构）。

### 2.3 "质量隙 = 最轻胶球质量"

由推论 A′：这个等式**当且仅当**底部是原子。若只证了 `MassGap`，正确的说法是"谱在 `(0,Δ)` 为空"，不能说"有一个质量 `Δ` 的胶球"。物理文献默认两者等价，是因为他们假设了粒子诠释（Haag–Ruelle/LSZ）；**在严格构造里这一步是额外的**。

### 2.4 接口对齐：confinement 与 physics_handoff

三处对接，都是记号/角色的对应，无冲突：

| `confinement/REPORT.md` §2.1 | 本报告 | 关系 |
|---|---|---|
| `E₀(L) = −lim_T (1/T) log⟨W(L×T)⟩` | —（弦扇区） | **不同谱扇区**；与 `Δ` 逻辑独立（他们 §2.1 的结论） |
| `Δ = inf{ E_n(0) : n 非真空 }` | `Δ = min_ρ m_ρ^{(1)}` | 本报告按对称通道 `ρ` **细化**；`E_n(0)` 是各通道的并 |
| `√σ/Λ` certificate，`Q(s)=√σ(s)/s` | `m/Λ` certificate，`q(s)=m(s)/s` | 同一交接判据 `q₀≥q₁`，**不变量不同**（他们 §4.4 已说明同构） |

**关键接口（对主 session，最重要的一条）。** RG 交接（`physics_handoff.md` 的 `q₀≥q₁`）控制的是**衰减率** `m`，即 MG 的"率"那一侧。定理 B / B′ 说明：

> **率 ⇒ 无原子。** 无论 RG 交接多么成功，它单独**不能**产出孤立单粒子态。原子需要独立的**状态计数**输入（近 `m₁` 的态有限）。

所以定理 D 被干净地切成两段、互不重叠：

```
RG 交接（主 session）      →  Δ > 0                 （率）
状态计数（本子问题）        →  Δ 被原子达到 + 上限隙   （原子性）
两者合用                    →  定理 D（孤立单粒子态）
```

> 状态计数本身**不比** RG 交接更难（O(3) 的 S5 是热迹倍增，结构与维数无关），但它需要的**输入**（方盒倍增的配分函数缺陷）正是 RG 输出——所以它是"同一堵墙上的第二根钉子"，而不是一条独立的墙。

---

## 3. 引理链 L0–L10

状态：✅ 现在可形式化；◐ 需搭建基础设施（有限群表示 / 复合算符）；❌ 公开。

| # | 引理 | 内容 | 依据 | 状态 |
|---|---|---|---|---|
| L0 | 转移矩阵 | 有限格点 `T = e^{−aH}` 对称正，`ker(T−1) = ℂΩ` | 反射正性 | ◐（格点层未形式化） |
| L1 | 有限维谱分解 | `T ψ_n = λ_nψ_n`，正交基，`λ` 降序 | mathlib `eigenvectorBasis` | ✅ **已有** |
| L2 | 质量与隙 | `m_n = −a⁻¹log λ_n`；`TransferGap ⇔ MassGap`，`Δ = m₁` | L1 + log 性质 | ✅ |
| L3 | 格点 KL | `⟨A T^t B⟩ = Σ λ_n^t ⟨Aψ⟩⟨ψB⟩`；`ν = Σ c_n δ_{m_n} ≥ 0` | L1 + 谱分解 | ✅ |
| L4 | **定理 A** | 原子 + 上限隙 ⇔ 尖锐指数衰减 | 测度论（控制收敛） | ✅ **~150 行** |
| L5 | **定理 B** | 反例：有隙、无原子、仍指数衰减 | 构造 `ρ(E)dE` | ✅ **~80 行** |
| L6 | 自由场非空洞 | `ν = δ_m`（`fourier_eq`） | OSforGFF | ✅ |
| L7 | 通道分解 | `ℋ = ⊕_ρℋ_ρ`，`Δ = min_ρ m_ρ^{(1)}` | 有限群表示 + L1 | ◐ |
| L8 | 强耦合格点隙 | `∃β₀, ∀β<β₀`：`λ₁ < 1` 一致 | OS78 / SZZ23 | ❌ 未形式化（大工程） |
| L9 | 连续统 `MG ⇒ Δ>0` | 体积一致聚类率 ⇒ 谱隙 | 主 session RG | ❌（条件于 MG） |
| L10 | 4D 原子 + 上限隙 | 定理 D | 需 (I1)+(I2) | ❌ **公开** |

**关键结构**：L1–L7 是"陈述层 + 等价 + 非空洞性"，**全部可证**，合计约 600–900 行，且**不需要**谱测度理论（mathlib 里没有 `spectralMeasure`，见 §8）。L8–L10 是真正的分析内容。

---

## 4. 综述

### 4.1 Källén–Lehmann 表示

KL 是"谱测度存在"的同义反复，其内容是：`W_O` 是 `[0,∞)` 上正测度的 Laplace 变换，因此**完全单调**（`(−1)^n W^{(n)} ≥ 0`）且由 Bochner 定理表征。三种等价形式：

| 形式 | 陈述 | 用途 |
|---|---|---|
| 时间方向 | `W(t) = ∫ e^{−iEt}dν(E)`，`ν ≥ 0` | 定理 A 的出发点；格点上精确 |
| 时空（Wightman） | `W(x) = ∫ dρ(m²) Δ₊(x;m)` | 相对论协变、`m²` 变量 |
| 动量空间 | `W(p) = ∫ dρ(m²)/(p²+m²)`；单粒子 = `δ(p²−m²)` | 极点语言："pole at `p² = −m²`" |

**严格性现状**：格点上由反射正性 + Hermiticity + 平移不变 + 多项式有界 ⇒ 正谱密度（Usui, `arXiv:1201.3415`）；连续统中随 OS 重构一并成立。这个格点版本正是 L3/L4 的直接来源，**比连续统版本更形式化友好**。

### 4.2 谱理论 / OS 重构 / 转移矩阵

三层结构：

```
格点转移矩阵 T（有限维，特征值）
  → 反射正性 ⇒ 物理内积 ⟨A,B⟩ := ∫ θ(A)* B dμ，T 对称正
  → OS 重构 ⇒ (ℋ, H, Ω)，T = e^{−aH}
  → spec(H) = {−a⁻¹ log λ_n}
```

**形式化现实**（已核对本机 vendored mathlib）：mathlib **没有**自伴算子谱定理 / 谱测度（无 `spectralMeasure`、无 `spectralProjection`），但**有**有限维对称算子的 `eigenvectorBasis` 与 `card_filter_eigenvalues_eq`（重数 = 特征空间维数）。因此正确的形式化策略是**两层**：

- **有限维层（L1–L3, L7）**：用 `eigenvectorBasis`，完全可做。
- **抽象测度层（L4–L5）**：用定理 A，绕开算子谱定理。
- **连续统无界算子谱定理：不做，也不需要。** 这是本子问题在形式化上**比预期容易**的原因。

### 4.3 OPE 的位置（UV 尾部，不是谱底）

OPE `O(x)O(0) ~ Σ_k c_k(x) O_k(0)` 的作用是**两件事，都与"是否存在孤立态"无关**：

1. **UV 尾部**：`O(x)O(0)` 的短距展开给出 `ρ_O(m²)` 在大 `m` 的行为（partonic 尾部，对数修正），以及可积性。以 `O = tr F²`（dim 4）为例，领头项是单位算符 `c₁(x)·1`，`c₁ ~ 1/x⁸`，其 Fourier 变换是动量多项式（接触项），对数项给出 `1/s²` 尾部。这决定 `ν` 是否有限、KL 积分是否收敛——**技术必需**。
2. **复合算符重整化**：`tr F²` 与低维算符混合，需 OPE 做系统减除。**没有这一步，胶球算符的谱测度根本没有定义。**

**结论**：OPE 是 (I3) 的工具，**不是** (I1)/(I2) 的工具。把 OPE 当成"证明胶球存在"的路线是范畴错误：OPE 是短距渐近展开，谱底是长距性质。

### 4.4 Fredenhagen–Marcu（带电扇区，与禁闭对偶）

- **原文**：K. Fredenhagen, M. Marcu, *Charged states in ℤ₂ gauge theories*, Comm. Math. Phys. **92** (1983) 81–119。
- **内容**：先证**没有局域算符能产生带电态**；再用"电荷对间距 `L → ∞`、能量一致有界"的规范不变态序列构造带电态。在屏蔽/禁闭区极限落在真空扇区；在自由（Higgs）区极限是**新扇区**（真正的带电态，GNS 给出平移协变表示 + 正 Hamiltonian）。真空分量系数 `c` 是序参量：`c ≠ 0` ⇒ 无带电态（禁闭），`c = 0` ⇒ 存在有限能量带电态。
- **FM 比**：`ρ = ⟨W_{γ₁}⟩⟨W_{γ₂}⟩/⟨W_γ⟩`（开 Wilson 线补成矩形），`ρ → 0` ⇒ 带电态、无禁闭；`ρ ≠ 0` ⇒ 禁闭。
- **与本报告的关系**：FM 完全在**带电/弦扇区**，对偶于禁闭与面积律（见 `research/confinement/`）。它**不**触及中性胶球谱。它的价值是一个方法论提示：**用"一致有界能量 + 间距 → ∞"提取谱信息**，与 O(3) 的 state-count 同源。
- **诚实标注**：我没有核实 FM 是否有第二篇专门题为"mass gap"的论文；本报告只依赖 CMP 92 (1983) 这一篇的内容，并把它标为背景而非引理。

### 4.5 格点胶球谱（数值）

标准方法：构造规范不变算符基（Wilson 环的多种形状 + 模糊化/APE smearing），投影到 `O_h` 的不可约表示（`A₁, A₂, E, T₁, T₂, …`）× `P` × `C`；用变分法（GEVP）从 `C(t) = Σ_n |⟨0|O|n⟩|² e^{−m_n t}` 提取 `m_n`。

Morningstar–Peardon 1999（`hep-lat/9901004`, PRD 60 034509，淬火 SU(3)，各向异性格点）：

| 通道 `J^{PC}` | 质量 (MeV) | 与 `0^{++}` 之比 |
|---|---|---|
| `0^{++}` | 1710 ± 50 ± 80 | 1 |
| `2^{++}` | 2390 ± 30 ± 120 | 1.39(4) |
| `0^{−+}` | 2560 ± 35 ± 120 | 1.50(4) |
| `0*^{++}`（激发） | 2670 ± 180 ± 130 | 1.54(11) |

**严格内容**：这些数字是**有限体积、固定截断的转移矩阵特征值**（§1.5(c)）。"`0^{++}` 最轻"是数值事实，**不是定理**；"连续统 `J^{PC}` 标注"需要 `a→0` 的群论极限（未形式化）。

**对证明的用处**：它们给出 (I2) 的**数值靶**：物理上 `δ ≈ 0.5`（`0^{++}` 激发在 `1.56 m₁`），且 `(m₁,m₂)` 内确实无态。但它们不是证明。

### 4.6 既有的"孤立质量"定理

| 模型 | 结果 | 方法 | 状态 |
|---|---|---|---|
| `P(φ)₂` 弱耦合 | Wightman 公理 + 粒子结构，单粒子质量在双粒子阈值下孤立 | Glimm–Jaffe–Spencer, Ann. Math. 100 (1974) 585 | ✅ 已证（未形式化） |
| `(φ⁴)₃` 弱耦合 | Wightman 公理 + 质量隙 | Feldman–Osterwalder, Ann. Phys. 97 (1976) 80 | ✅ 已证 |
| 2D O(3) 连续统 | 谱测度在最低质量处有正原子，与其余谱分离 | OpenAI, 2026-10-04 | ✅ 声称 + 论证（未独立审计） |
| 4D 纯 YM | **质量隙未证**；孤立单粒子态**未证** | — | ❌ |

**O(3) 论文是本子问题的模板**，其定理 1.1 与本文定理 D 同构（`ρ = Zδ_{m₁²} + ρ_rest`，`supp ρ_rest ⊂ [(m₁+ε)²,∞)`）。§5 逐条分析哪些步骤能移植。

---

## 5. O(3) 模板的移植分析（逐条）

O(3) 的证明链（`_o3-pole.txt`，本目录）：

| 步骤 | O(3) 的机制 | 4D YM 的状态 | 判定 |
|---|---|---|---|
| **输入** | 连续统构造 + 方盒倍增的配分函数缺陷 + 阶乘矩界 | 16 输入给 OS；倍增缺陷在 4D 是 RG 输出 | ◐ 条件 |
| **S1 矩形迹界** | 压力界 + 方盒比较 ⇒ 每个大周长有正隙 | 同构（热迹论证是维数无关的） | ✅ 可移植 |
| **S2 偶扇区 ≥ 2e** | `q → −q` 的 Z₂；偶多项式协方差 ≤ 两点函数² ⇒ 偶非真空能量 ≥ 2e | **无对应 Z₂**：C 与 P 都不把单胶球与双胶球分开（单胶球 `gg` 与双胶球 `(gg)(gg)` 都 C-even、P-even） | ❌ **失败（关键）** |
| **S3 可见性** | FKG + Edwards–Sokal + KSS 穿越定理 ⇒ 最低奇本征函数与自旋场有界重叠 | SU(N) 无 FKG；替代：Ginibre/OS 不等式 + 特征标/polymer 展开，仅在强耦合 | ◐ 部分（仅强耦合） |
| **S4 `liminf e ≥ m`** | S3 + 局域电路 | 需 S3 的定量版本 | ❌ 依赖 S2/S3 |
| **S5 状态计数** | `N(3m/2) ≤ Cw`：热迹倍增 + 反射正性 + 能量截断的正多项式逼近 | 需 4D 的配分函数倍增缺陷（RG） | ◐ 结构可移植，输入缺 |
| **S6 原子** | `J` 个近 `m` 的质量 ⇒ `≥ Jc₀w` 个态 ⇒ 有限 ⇒ 最小点孤立 ⇒ 正测度给原子 | 一般（只用正性 + 计数） | ✅ 给定 S5 |

**核心结论**：S2 是 O(3) 证明的支点，而它**不可移植**。O(3) 里"最低激发是奇态、所有偶态 ≥ 2e"是**对称性 + 正性**的直接推论；YM 的 `0^{++}` 通道没有这个结构，所以"单胶球与双胶球阈值之间无态"必须靠**动力学**（粒子诠释），不能靠对称性。

**一个可用的类比点（本报告的独立观察）。** 若改取 **C 奇**（或 **P 奇**）通道的算符 `O`，则：
- `ν_O` 在所有 C-even 态上为零，特别是 `0^{++}` 胶球；
- 该通道最轻的单粒子态是 `0^{−+}`（`≈2.56 GeV`）；
- **任何 C 奇的（渐近）多粒子态**至少含一个 C 奇粒子，故最轻者是 `0^{−+}` + `0^{++}`，能量 `≥ m_{0^{−+}} + m_{0^{++}}`（`≈4.3 GeV`），远高于 `m_{0^{−+}}`。

于是 S2 的类比成立：**C 奇通道的两粒子阈值被 Z₂ 抬到 `m_{0^{−+}} + m_{0^{++}}`，上限隙是"对称性事实"**——只要粒子诠释（阈值以下由单粒子态穷尽）可用。这把"最难的一步"从 `0^{++}` 通道挪到一个有 Z₂ 保护的通道。代价：`0^{++}` 是物理上最轻的胶球，C 奇通道证明不了它；但定理 D 只要求**存在**一个孤立单粒子态，C 奇通道足够。**建议把 C 奇通道作为第一个攻击目标。**

---

## 6. 已证 / 未证状态表

| 命题 | 状态 |
|---|---|
| 定理 A（原子 ⇔ 尖锐衰减） | ✅ **可证**（测度论，~150 行） |
| 定理 B（隙 ⇏ 原子） | ✅ **可证**（~80 行） |
| 定理 C（格点谱分解；`Δ = min_ρ m_ρ`） | ✅/◐（L1–L3 ✅，L7 需有限群表示） |
| 定理 E（自由场原子） | ✅（OSforGFF `fourier_eq`） |
| 格点强耦合质量隙（`β<β₀`） | ✅ 文献（OS78/SZZ23），未形式化 |
| 连续统 YM 质量隙 | ❌（条件于 MG） |
| 4D 质量隙底部是原子 (I1) | ❌ **公开** |
| 4D 上限隙 (I2) | ❌ **公开** |
| 胶球算符重整化 (I3) | ❌（形式化缺口） |
| 连续统 `J^{PC}` 标注 | ❌（`O_h` subduction 未形式化） |

---

## 7. 障碍（按致命性排序）

**O1 — 原子（底部可达）。** MG 只给衰减率；定理 B 表明"率 `Δ` + 正测度"**不**给原子。要原子必须有**状态计数**（近 `m₁` 的态有限）。这是超出 MG 的**新**数学输入，且不是对称性/正性的免费推论。**最致命。**

**O2 — 上限隙。** `(m₁,m₂)` 无谱。O(3) 的 S2 机制在 YM 失效（无分离单/多胶球的 Z₂）。C 奇通道可绕开（§5），但需粒子诠释。**次致命。**

**O3 — 可见性。** 需规范不变算符 `O` 与最低态有**定量**非零重叠。O(3) 用 FKG/Edwards–Sokal/KSS；SU(N) 无 FKG。替代是特征标/聚合物展开，**只在强耦合**给出下界；弱耦合（连续统）无。**技术硬点。**

**O4 — 4D 状态计数。** S5 需要"方盒倍增的配分函数缺陷一致于截断"，即 RG 输出——与 MG **同一堵墙**。注意：这一步**不**比 MG 更难（结构相同），只是同一障碍的另一个投影。

**O5 — 复合算符重整化。** `tr F²` 是 dim-4 算符，需减除与混合；OSforGFF 只有基本场，**没有复合算符理论**。没有 (I3)，定理 D 的陈述层都是空的（"`ν_{tr F²}`" 未定义）。这是一个**独立的形式化缺口**，不是分析障碍。

**O6 — 连续统 `J^{PC}`。** 格点 `O_h` 不可约表示在 `a→0` 下混合成 `SO(3)` 的 `J`；subduction（`A₁ → 0⁺, T₁ → 1, E⊕T₂ → 2, …`）是群论，但"连续统谱按 `J` 分解"需要重构定理的协变性，未形式化。影响"胶球谱刻画"（定理 C(d) 的连续统版），不影响定理 D。

**O7 — 有限体积原子不持久。** 有限体积的特征值在无限体积可以消失成连续谱。O(3) 论文**显式避免**假设其持久（"does not assume that an eigenvalue at finite circumference persists"），用 state-count 绕开。这是方法上的正确选择，也说明 (I1) 不能用"有限体积特征值 + 紧性"直接得到。

---

## 8. 下一具体步骤

### 8.1 第一步（可立即做，2–4 周）：陈述层 + 等价 + 非空洞性

新建 `lean/YMStatement/Spectrum.lean`（约 600–900 行，目标 0 sorry；依赖 `Mathlib` 与 `OSforGFF`，不依赖 `RGHandoff`）。内容 = L1–L6：

```lean
import Mathlib
import OSforGFF.OS.Axioms
import OSforGFF.Covariance.Propagator
import YMStatement.Basic

noncomputable section
namespace YMStatement.Spectrum

open MeasureTheory
open scoped InnerProductSpace

/-! ## L1–L2: 有限维转移矩阵 -/

variable {𝕜 E : Type*} [RCLike 𝕜] [NormedAddCommGroup E] [InnerProductSpace 𝕜 E]
  [FiniteDimensional 𝕜 E] {n : ℕ}

/-- 格点质量：`m = −a⁻¹ log λ`。-/
noncomputable def massOf (a λ : ℝ) : ℝ := -(Real.log λ) / a

/-- 转移矩阵的谱隙 ⇔ 质量隙。`λ₁` 为第二大特征值（`eigenvalues` 降序）。-/
theorem transferGap_iff_massGap
    (T : E →ₗ[𝕜] E) (hT : T.IsSymmetric) (hn : Module.finrank 𝕜 E = n)
    (a : ℝ) (ha : 0 < a) (hλ : 0 < hT.eigenvalues hn 1) :
    (hT.eigenvalues hn 1 < 1) ↔
      (∃ Δ > 0, ∀ i : Fin n, i ≠ 0 → hT.eigenvalues hn i ≤ Real.exp (-a * Δ)) := by
  sorry

/-! ## L3: 格点 KL 表示 -/

/-- `⟨Ω| A T^t B |Ω⟩ = Σ_n λ_n^t ⟨Ω|Aψ_n⟩⟨ψ_n|BΩ⟩`。-/
theorem twoPoint_spectral
    (T : E →ₗ[𝕜] E) (hT : T.IsSymmetric) (hn : Module.finrank 𝕜 E = n)
    (A B : E →L[𝕜] E) (Ω : E) (hΩ : T Ω = Ω) (t : ℕ) :
    ⟪Ω, A ((T ^ t) (B Ω))⟫_𝕜 =
      ∑ i : Fin n,
        (hT.eigenvalues hn i) ^ t * ⟪Ω, A (hT.eigenvectorBasis hn i)⟫_𝕜 *
          ⟪hT.eigenvectorBasis hn i, B Ω⟫_𝕜 := by
  sorry

/-! ## L4 定理 A: 原子 ⇔ 尖锐指数衰减 -/

/-- 零动量两点函数：`W(t) = ∫ e^{−Et} dν(E)`。-/
noncomputable def correlator (ν : Measure ℝ) (t : ℝ) : ℝ :=
  ∫ E, Real.exp (-E * t) ∂ν

/-- 孤立原子 + 上限隙。-/
def IsolatedAtom (ν : Measure ℝ) (m₁ m₂ : ℝ) : Prop :=
  0 < m₁ ∧ m₁ < m₂ ∧ ν {m₁} ≠ 0 ∧ ν (Set.Ioo 0 m₂ \ {m₁}) = 0

/-- 尖锐指数衰减。-/
def SharpDecay (ν : Measure ℝ) (m₁ Z m₂ C : ℝ) : Prop :=
  0 < Z ∧ m₁ < m₂ ∧ 0 ≤ C ∧
    ∀ t : ℝ, 0 ≤ t → |correlator ν t - Z * Real.exp (-m₁ * t)| ≤ C * Real.exp (-m₂ * t)

/-- **定理 A**：`ν` 有限、`ν([0,∞))` 全质量、`0 < m₁ < m₂` 下二者等价。-/
theorem isolatedAtom_iff_sharpDecay
    (ν : Measure ℝ) (hν : IsFiniteMeasure ν) (hpos : ν (Set.Iic 0) = 0)
    {m₁ m₂ : ℝ} (hm : 0 < m₁) (hm' : m₁ < m₂) :
    IsolatedAtom ν m₁ m₂ ↔ ∃ Z C, SharpDecay ν m₁ Z m₂ C := by
  sorry

/-! ## L5 定理 B: 反例 -/

/-- 存在有隙、无原子、仍按 `e^{−Δt}` 衰减的谱测度。-/
theorem massGap_not_isolated :
    ∃ (ν : Measure ℝ) (Δ : ℝ), 0 < Δ ∧ IsFiniteMeasure ν ∧
      ν (Set.Iic 0) = 0 ∧
      (∀ t : ℝ, 0 ≤ t → correlator ν t ≤ (ν Set.univ).toReal * Real.exp (-Δ * t)) ∧
      (∀ m, ν {m} = 0) := by
  sorry

/-! ## L6 自由场非空洞性 -/

/-- 自由场：`ν = δ_m`（由 `GFFPropagator.fourier_eq`），故对任意 `m₂ > m` 孤立。-/
theorem freeField_isolated (d : ℕ) [Fact (2 ≤ d)] (m : ℝ) [Fact (0 < m)]
    [GFFPropagator d m] (m₂ : ℝ) (hm₂ : m < m₂) :
    IsolatedAtom (Measure.dirac m) m m₂ := by
  sorry

end YMStatement.Spectrum
```

**验收**：`lake build` 通过、0 sorry、`#print axioms` 只剩 `propext, Classical.choice, Quot.sound`。四个定理（`transferGap_iff_massGap`, `twoPoint_spectral`, `isolatedAtom_iff_sharpDecay`, `massGap_not_isolated`）加一个非空洞性测试。

**为什么这是正确第一步**：它是 `docs/yang-mills-prior-work.md` §5.3 的 M1 里程碑在**谱**方向的精确化，且**不需要** mathlib 里不存在的谱测度理论（已核实：mathlib 无 `spectralMeasure`/`spectralProjection`；但有 `eigenvectorBasis`）。它产出的是"非空洞陈述 + 等价 + 反例测试"，正是 §4.5(反面教材) 要求的三件套。

### 8.2 第二步（2–4 周，可证但需搭建）：通道分解

`lean/YMStatement/GlueballChannel.lean`：有限群 `G_h = O_h × {1,C}` 的特征投影、`ℋ = ⊕_ρℋ_ρ`、`Δ = min_ρ m_ρ^{(1)}`（定理 C(c)(d)）。需要：有限群特征标正交性（mathlib 有部分）+ 投影算子自伴 + `T` 保同型分量。与 `confinement` 的 character-expansion 库**共享**（那里 L1 也要 Peter–Weyl）。

### 8.3 第三步（长期，公开）：上限隙

按 §5 的建议，先攻 **C 奇通道**：
1. 定义 C 奇规范不变算符（例如 `tr F² F̃²` 型或 `ε^{μνρσ}tr F_{μν}F_{ρσ}` 的 smeared 版本），需 (I3) 重整化；
2. 给出阈值 `m_{0^{−+}} + m_{0^{++}}` 的严格下界（用粒子诠释/渐近完备，未形式化）；
3. 状态计数（S5，需 RG 输入）；
4. 用定理 A 收口。

**诚实的期望**：第 2 步是**物理输入**而非定理；把它作为假设写出，得到形如
```
(H_particle) + (H_count) + MG  ⇒  ∃ 孤立单粒子态（C 奇通道）
```
的条件定理。这是本子问题在可见未来能达到的最强形式。

---

## 9. 文献

**KL 与谱表示**
- G. Källén, *Helv. Phys. Acta* 25 (1952) 417；H. Lehmann, *Nuovo Cim.* 11 (1954) 342.
- K. Usui, *A Note on Reflection Positivity and the Umezawa–Kamefuchi–Källén–Lehmann Representation of Two Point Correlation Functions*, [arXiv:1201.3415](https://arxiv.org/abs/1201.3415)（格点上由反射正性证正谱密度）。
- R. Haag, B. Schroer, *J. Math. Phys.* 3 (1962) 248（KL 与谱条件）。

**OS、重构、粒子**
- K. Osterwalder, R. Schrader, CMP 31 (1973) 83；42 (1975) 281.
- J. Glimm, A. Jaffe, T. Spencer, *The Wightman axioms and particle structure in the P(φ)₂ quantum field model*, Ann. Math. 100 (1974) 585–632（**孤立单粒子质量的经典定理**）。
- J. Feldman, K. Osterwalder, Ann. Phys. 97 (1976) 80–135.
- R. Haag, B. Swieca, *When does a quantum field theory describe particles?*, CMP 1 (1965) 308–320.
- D. Buchholz, E. H. Wichmann, *Causal independence and the energy-level density of states in local quantum field theory*, CMP 106 (1986) 321–344.

**Fredenhagen–Marcu**
- K. Fredenhagen, M. Marcu, *Charged states in ℤ₂ gauge theories*, CMP 92 (1983) 81–119.
- 后续严格结果：*The phase transition of the Marcu–Fredenhagen ratio in the abelian lattice Higgs model*, Electron. J. Probab. (2024), DOI 10.1214/24-EJP1183。

**格点胶球谱**
- C. Morningstar, M. Peardon, *The glueball spectrum from an anisotropic lattice study*, PRD 60 (1999) 034509, [hep-lat/9901004](https://arxiv.org/abs/hep-lat/9901004).
- C. Morningstar, *What do we know about glueballs and hybrids from lattice QCD*（综述幻灯片）。
- A. Athenodorou, M. Teper, 2020（非淬火标量胶球质量抑制）。

**4D YM 构造与 MG**
- `docs/yang-mills-prior-work.md`（本仓库）。
- OpenAI, *An Isolated Particle Pole for the Two-Dimensional O(3) Spin Field*, 2026-10-04（**本子问题的直接模板**）。
- hexagonmath, *From Four-Dimensional Wilson Lattice Gauge Theory to Quantum Yang–Mills Theory*, 2611.00492（16 输入 + MG 归约 + 反例）。
- `research/physics_handoff.md`（RG 交接控制量）；`research/confinement/`（面积律，带电扇区）。

**Lean**
- `github.com/mrdouglasny/OSforGFF`（本仓库 `vendor/OSforGFF`）；`GFFPropagator.fourier_eq`、`gaussianFreeField_satisfies_OS4_PolynomialClustering`。
- mathlib `Mathlib/Analysis/InnerProductSpace/Spectrum.lean`：`eigenvectorBasis`、`card_filter_eigenvalues_eq`（**无** `spectralMeasure`）。

---

## 附录 A：一页速查

```
层次：IsolatedOneParticle ⇒ SharpDecay ⇒ MG ⇒ MassGap（反向全不成立）
Clay 只要 MassGap；本子问题严格更强，不在获奖关键路径上

可证（现在）：
  定理 A：ν 有原子 m₁ + 上限隙 m₂  ⇔  W(t) = Z e^{−m₁t} + O(e^{−m₂t})   [测度论，~150 行]
  定理 B：存在 ν，有隙 Δ、无原子、仍按 e^{−Δt} 衰减                     [~80 行]
  定理 C：格点 T 对称正 ⇒ m_n = −a⁻¹log λ_n，Δ = min_ρ m_ρ^(1)           [eigenvectorBasis]
  定理 E：自由场 ν = δ_m                                                 [OSforGFF fourier_eq]

未证（分析）：
  (I1) 底部可达 = 原子        ← MG 不给；需状态计数（近 m₁ 的态有限）
  (I2) 上限隙 (m₁,m₂) 无谱    ← O(3) 靠自旋翻转 Z₂，YM 的 0^{++} 无此 Z₂
  (I3) tr F² 重整化           ← 形式化缺口（无复合算符理论）

关键观察：
  C 奇通道的上限隙是 Z₂ 保护的：阈值被抬到 m_{0^{−+}} + m_{0^{++}}
  ⇒ 建议第一个攻击目标是 C 奇通道，而不是 0^{++}

形式化策略（重要）：
  mathlib 无谱测度/无界谱定理；但有有限维 eigenvectorBasis
  ⇒ 有限维转移矩阵层 + 抽象测度层（定理 A），完全绕开无界算子谱定理

数值靶（Morningstar–Peardon 1999, SU(3)）：
  0^{++} 1710, 2^{++} 2390, 0^{−+} 2560 MeV；0^{++} 激发/基态 ≈ 1.56
```

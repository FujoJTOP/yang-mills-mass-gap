# L2C-SD — `G″ = S − D` 里的 `D` 到底是什么：`D ≠ Var(G)`（`D = Var(G+Φ′)`），`D` 与 `1−f₁²` **无函数关系**（类型不同：矩 vs 原子），`D/S≈2%` 是**传播错误**（真值 `1−D/S≈3%`，`D/S≈0.97`）

**日期**：2026-11-25 · **只写** 本文件 · **禁 git commit**（遵守）· 本轮 ≤25 min
**输入（逐字核过，行号在下）**：`docs/paper-v2/ch18-scalar-h1.tex:28–29/33–35/50–70/110–128/139–150`、
`docs/paper-v2/ch05-theorems.tex:233–240`（`eq:split`）、`docs/paper-v4/sec02-scalar.tex:94`（字母版 `S,D` 的**唯一**逐字定义句）、
`research/unconditional/INVENT-VARIANCE-INEQ.md:3/8/10/19/27/30`、`research/v3/species/P01-WHATG.md:14/30–34/57–71/78/110–112`、
`research/v3/species/DS-RATIO.md:12/15–19/21–25/42–51/143–158/233–253/274–297`、
`research/unconditional/PROVE-Gpp-positive.md:7–8/31`、`docs/paper-v2/app.tex:132/232–240/468–469/754`、
`research/v3/species/T01-IMPLY.md:30–36/64–67/96–102`、`research/v3/R-WAVE-SYNTHESIS.md:45`、
`research/unconditional/A4-SPECTRAL-BOOTSTRAP.md:13`、`research/v3/Q6-H1-STATUS.md:35/118/214–216`、
`research/v3/P6-PERIODIC-TO-H1.md:37/235/246`、**`docs/paper-v6/sec07-assembly.tex:70`**（本项目**自己撤回**`G″`-槽与`ν`-槽可比的句子）。
**标签**：`[T]` 真定理/逐字恒等式 · `[C]` 条件定理 · `[V]` 本轮可复算数值 · `[E]` 未过自校验、不作证据 · `[O]` 精确缺失步 · `[J]` 判断。
**无 QED。Clay 0%。**

---

## §0 判词（先给）

| 问题 | 判定 |
|---|---|
| **`S`、`D` 的精确式** | `G″(0) = S − D`，`S = E[H_{vv}+Φ″|0]`，`D = Var(G_v+Φ′|0)`，`H_{vv}=∂_v²S_W`，`G_v=∂_vS_W`，`Φ=−log J`。**`[T]`**（`sec02-scalar.tex:94` 逐字；`ch18:28–29` 同式无字母版）。 |
| **`D` 是不是 `Var(G)`** | **不是。** `D = Var(G)+Var(Φ′)+2Cov(G,Φ′)`；`Var(G)` 只是首项。语料自己的三分为 `Δ_dyn=E[H]−Var(G)`、`W_geo=E[Φ″]−Var(Φ′)`、`X=−2Cov` ⇒ `S−D=Δ_dyn+W_geo+X`。**`[T]`**（`P01:48` 一行代数；`ch18:29`）。任务书前提"`D` 是 `Var(G)`"**逐字为假**，只在 `Φ≡0`（Lebesgue 平测度）时成立。 |
| **`D/S≈2%` 出处** | **不是 `D/S`**。它源自 `PROVE-Gpp-positive.md:8` 的 **`G″(0)=1/Var(a)` 吻合到 ≈2%**（胶子模、弱耦、SU2 `0.274 vs 0.268`；SU3 `0.485 vs 0.45`），经 `INVENT-VARIANCE-INEQ.md:3/30` 改写成"measured `D/S≈2%`"后逐层传播（`new-tools-registry.md:118/266`、`CLAUDE.md:37`、`YANG-MILLS-MASTER.md:72`、`paper/main.tex` 系 body4:41、`ch13-machine.tex:43`）。**真值：`1−D/S≈3%`，`D/S≈0.97`**（`DS-RATIO:233` `[T]`+`[V]`；`app.tex:132/469` 独立给 `0.94–1.00`、5/12 模 `G″<0`）。 |
| **`D` 与 `1−f₁²` 有函数关系吗** | **没有。** 类型不同（**矩 vs 原子**），且存在**精确族**证明 `D` 不是 `f₁²` 的函数：单 plaquette 上 `f₁²` 无定义/退化而 `D=Bβ²` 无界；`β=0` Haar 与 `β=10` 单 plaquette 同 `f₁²=0` 而 `D=0` vs `1.106`。**`[T]`** |
| **判死线** | **(b) 成立 —— L2 的桥不存在**；**(c) 同时成立 —— `2%` 的来源不可靠**（是补、是分裂商、且无归档 witness）。**(a) 不成立**，且 `(a)` 形在量纲上必含第三个未知标度。 |

**一句话**：`D` 是**条件 Gibbs 纤维上力变量的二阶中心矩**（`L²(μ_{a₀})` 对象，`P01:68–71`），`f₁²` 是**谱测度 `ν` 在极点上的原子**（`ν({m₁²})=|⟨Ω|O|𝔾⟩|²`，`LIT-ERA-05:310`/`R-WAVE:41`）。**前者是"矩"，后者是"原子"；所有精确定理都只把原子从**上方**夹（`f₁² ≤ M₀`、Christoffel–Darboux），从无一条把矩从上方用 `1−f₁²` 夹住。**

---

## §1 逐字核 `S`、`D`（任务 1）

### 1.1 字母版定义（`S,D` 的**唯一**逐字句）

`docs/paper-v4/sec02-scalar.tex:94`：

> `$\mathrm{gap}(\theta):=S-D=\Gpp$ with $S=E[H+\Phi''\mid0]$ and $D=\Var(G+\Phi'\mid0)$`

即

```
S := E_μ[ H_vv + Φ″ | a₀=0 ] ,      H_vv = ∂_v² S_W
D := Var_μ( G_v + Φ′ | a₀=0 ) ,     G_v  = ∂_v S_W ,   Φ = −log J
G″(0) = S − D = Δ_dyn + W_geo + X
```

### 1.2 无字母版（同一式，`ch18`）

`ch18-scalar-h1.tex:28–29` 逐字：

> `G''(0)=E[H_{vv}+\Phi''\mid\xi=0]-\operatorname{Var}(G_v+\Phi'\mid\xi=0)=\Delta_{\mathrm{dyn}}+W_{\mathrm{geo}}+X`,
> `with $\Delta_{\mathrm{dyn}}=E[H]-\operatorname{Var}(G)$`

`ch18:33–35`（`eq:scalar-identity`）把 `Var` 显式写成 `Var_μ(⟨v,∇S⟩|0)`。`ch05-theorems.tex:233–237`（`eq:split`）同式，标签在 `:239`。

**核验记录（AGENTS §1）**：`grep` 工具 **OK**；检索式 `S - D|S-D|S−D` 于 `ch18-scalar-h1.tex` → **命中 0**（该章通篇用 `Var/E[H]/W_full`，不用字母 `S,D`；字母版只在 `paper-v4/sec02-scalar.tex:94` 与 `ch05-theorems.tex:548` 的表格行）。**记 0，不记"没有"。**

### 1.3 `S`、`D` 各自装了什么（⇒ 与 `T01:35` 的对应）

```
S = E[H]  +  E[Φ″]                                  ← 两个"平均曲率"（Wilson + Jacobian）
D = Var(G) + Var(Φ′) + 2Cov(G,Φ′)                   ← 三个"噪声"项
S − D = [E[H] − Var(G)]  +  [E[Φ″] − Var(Φ′)]  −  2Cov(G,Φ′)
        =     Δ_dyn       +        W_geo       +       X
```

**⇒ 任务书问"`S` 与 `D` 分别对应哪两个"的精确答**：

| `T01:35` 的项 | 取自 `S` 的部分 | 取自 `D` 的部分 |
|---|---|---|
| `Δ_dyn` | `E[H]`（Wilson 曲率均值） | `Var(G)`（Wilson 力方差，**带负号**） |
| `W_geo` | `E[Φ″]`（Jacobian 曲率均值，`≥N/12` `[T]`） | `Var(Φ′)`（Jacobian 力方差，**带负号**） |
| `X` | **无**（`S` 是线性期望，没有交叉项） | `2Cov(G,Φ′)`（**带负号**，整项来自 `D`） |

即 **`S` 只装"两个均值"，`D` 装"两个方差 + 一个协方差"**；`T01:35` 的三个减号项与 `D` 的三项**一一对应**，`S` 的两项与 `Δ_dyn`、`W_geo` 的正部**一一对应**。`[T]`（`P01:48` 自算一行代数，与 `ch18:29` 一致）。

### 1.4 两处必须钉死的逐字事实

1. **`D` 不是 `Var(G)`。** `Var(G)` 只在 `D` 里占一项；`Φ≡0`（Lebesgue，如 M6）时才有 `D=Var(G)`。**任务书该前提逐字为假。**
2. **`S` 不是作用量。** 字母 `S` 在本项目里撞车：LaTeX 正文里 `S`（=`S_W`）是 Wilson 作用量；`G″=S−D` 里的 `S` 是 `E[H+Φ″]`。**`D` 同样撞车**：与 `DS-RATIO:48` 的 H2 分裂 `S=E[H], D=Var(V′)−E[Φ″]` 不同。

### 1.5 还有**第二对** `(S,D)`：SCBL 分歧版（这决定了 `D/S` 不定值）

`INVENT-VARIANCE-INEQ.md:8` 逐字（`G″=S−D` 的**原始**推导）：

> `G″_0 = S(v) − D(v)`，`S:=E_μ[h_n−⟨K^{−1}c,c⟩|0]`，`D:=Var_μ(f|0)−E_μ[⟨K^{−1}c,c⟩|0]`
> （`h_n=vᵀHess V v`，`f=∂_nV`，`Π₀=I−vvᵀ`，`c=Π₀ Hess V v`，`K=Π₀ Hess V Π₀`）

**⇒ 同一个字母 `D` 在语料里有两个精确式**：

* **`D_A2 := Var(G+Φ′)`**（`sec02:94`、`ch18:29`，**无减法**）——恒有 `D_A2 ≥ 0`；本文件 §1.3 的对应表是按这个版本写的。
* **`D_SCBL := Var(f) − E⟨K⁻¹c,c⟩`**（`INVENT:8`）——**可正可负**，且在 `β=2.5` 的 Savvidy 不定纤维上 `K` **不可逆 ⇒ 无定义**（`INVENT:30`「`det K=0 ⟺ 横向无质量模 ⟺ S=0`」；`DS-RATIO:47/158` 记 `tci_comp.log` 的 `elim? UNDEF`）。

`S_A2 − D_A2 ≡ S_SCBL − D_SCBL ≡ G″` `[T]`（两者只差一个平移 `ψ=E⟨K⁻¹c,c⟩`，`DS-RATIO:143–147`），**但比值 `D/S` 随分裂变**。这不是估计误差，是**定义层的自由度**。

---

## §2 `D/S ≈ 2%` 的出处（任务 2）

### 2.1 传播链（每一跳都给了文件:行）

| 跳 | 内容（逐字/近似逐字） | 出处 | 级别 |
|---|---|---|---|
| **源头** | 「Numerically **`G″(0)=1/Var(a)` to ≈2%**（SU2: `0.274 vs 0.268`；SU3: `0.485 vs 0.45`；quartic coeff ≈0），while the **literal `E[H]−Var` ≈ 0 (+0.06, −0.01)**」 | `research/unconditional/PROVE-Gpp-positive.md:8` | `[V]` 数值比较 |
| **改写** | 「measured `D/S≈2%`」；「Numerically `D/S≈2%` (the measured **Gaussian-exactness** `G″≈1/Var`)」 | `INVENT-VARIANCE-INEQ.md:3`、`:30` | **改写本身无 witness** |
| **传播** | `G″(0)>0`…`=S−D`，`D/S≈2%`；「β=2.5 缺陷 `D/S≈2%`」；「D/S≈2%」 | `docs/new-tools-registry.md:266`、`:118`；`CLAUDE.md:37`；`docs/YANG-MILLS-MASTER.md:72`、`:178`；`docs/paper/body4.tex:41`；`docs/paper-v2/ch13-machine.tex:43`；`docs/paper-v2/ch05-theorems.tex:548` | 全部**二手** |
| **同书反面** | 「the `D/S≈2%` quoted in the main text is a **compactness/Haar number measured on that gauge-dependent mode — not gap evidence** — and is not robust under a susceptibility normalisation (`−0.45` to `+0.60` at strong coupling; the `2%` is a **weak-coupling value**)」；「**`D/S≈2%` is a gauge-dependent number.** It is measured on the gluon mode at weak coupling, **not** on the gauge-invariant `0⁺⁺` at the physical coupling」 | `docs/paper-v2/app.tex:238–240`、`:754`（另 `:132/469` 给 `D/S ∈ 0.94–1.00`、5/12 模 `G″<0`） | `[T]` 级别警告 |
| **本轮裁定** | `2%` = **`1−D/S`**；归档数据重算 `D/S_A2=0.9718`、`D/S_H2=0.9651`、`1−0.97=0.03`；精确单 plaquette 给闭式 `D/S=β²/β*²`（`β*=4.30661`），β=2.5 处 **0.33698**；`0⁺⁺` 上**无**任何 `D/S` 测量 | `DS-RATIO.md:12/21–25/57–68/233–235/274–297` | `[T]`+`[V]` |

### 2.2 判定：**不是实测，是"两个不同量的吻合度"**

- **它不是 `D/S` 的测量**：源头 `PROVE-Gpp-positive:8` 比的是 `G″` 与 `1/Var(a)`，而 `1/Var(a)` 是**乘积恒等式**的量，**不是** `S`（`DS-RATIO:49` 的 `PROD` 行；`:296` 实测 `S_A2·Var(a)=6.32≠1`，边缘远非高斯）。
- **它不是 `β=2.5` 的那个点**：`app.tex:240` 逐字说这是**弱耦值**；SU3 那半个 `0.485 vs 0.45` 在 `β=5.7`。
- **它的对象不是 `0⁺⁺`**：`app.tex:754` 逐字；`0⁺⁺` 的 `D/S` 是**数据缺口**（`DS-RATIO:306`：只有 `Var(q̄)∝V^{-1}`，做不出曲率比值）。
- **它连"单值"都不是**：`DS-RATIO:233` `[T]`——`D/S` 是**分裂商**，同一 `G″` 下 `0.02 / 0.34 / 0.97 / 0/0` 皆可。

**⇒ (c) 成立：`D/S≈2%` 的来源不可靠，应改写为「`1−D/S≈3%`（近抵消），且 `D/S` 是分裂商」。**

---

## §3 ★ 关键：`D` 与 `1−f₁²` 的关系（任务 3）

### 3.1 先把两个对象钉在各自的槽里（逐字）

| | `D`（或 `Var(G)`） | `f₁²` |
|---|---|---|
| 定义 | `Var_μ(G_v+Φ′|a₀=0)`，`G_v=∂_vS_W` | `f₁² = ν({m₁²}) = |⟨Ω|O|𝔾⟩|²`，`O=¼trF²` |
| 概率空间 | **条件 Gibbs 纤维** `(Ω_{a₀}, μ_{a₀})`，`μ_{a₀}∝e^{−S_W}J dX δ(⟨v,X⟩−a₀)` | **Källén–Lehmann 谱测度 `ν`**（真空扇区，有 Hilbert 空间、有能量极点） |
| 类型 | 实值随机变量的**二阶中心矩**（`L²(μ)` 对象，**无时间**、无极点） | 测度的**原子**（极点留数，**需要**转移矩阵 + 孤立本征态） |
| 逐字禁配 | `P01-WHATG.md:71`：「**不是**…不是算符（**无 Hilbert 空间**）」；`:78`：「**不是谱测度 `ν` 上的期望**」；**主文自撤**：`docs/paper-v6/sec07-assembly.tex:70` 逐字「`Δ_0` and `Λ_★` **live on different carriers（`G″` vs `supp ν`）, are of different type (position vs distance), and have different quantifier domains, so they are incomparable**」—「That identification is withdrawn」 | `LIT-ERA-05.md:310`：「`f₁² := ν({m₁²}) = |⟨0|O|G⟩|² > 0`」；`R-WAVE:41`：「the missing step is the **state overlap** `Ω=|⟨G_tgt|O|G_eff⟩|²=1`」 |
| 量纲 | `[D]=[G″]=[a]^{−2}`（`P01:110–112` 逐字：`[E[H]]=[Var(G)]=[a]^{−2}`，`G″_a=a²m_phys²`） | **随归一三变**：`4`（`A-THEOREM.md:35` 未归一 `|⟨0|trF²|m₁⟩|²`）/ `2`（`A1-SUMRULE-A.md:42–43` 即 `f₁²/χ≤m₁²` 的归一）/ `0`（`ch02` 的原子分数） |

**⇒ 第一层答案**：`f₁²` 是**矩阵元平方 = 极点处的原子**；`D` 是**方差 = 构型测度上的二阶矩**。**不是同一类数学对象，也不是同一个随机变量的方差**——连"都是方差"这个退让读法都不成立，因为承载原子的矩是 `M₀=∫dν = Var(¼trF²)`（一个**规范不变**复合算符的方差），而 `D` 是**规范相关力** `G_v=∂_vS_W` 的方差（另加 `Var(Φ′)` 与 `2Cov`）。

### 3.2 为什么"矩 ↔ 原子"只有单向（`[T]`）

- **原子 ≤ 总质量**：`f₁² = ν({m₁²}) ≤ M₀ = ∫dν`；Christoffel–Darboux 更紧：`f₁² ≤ 1/K_n(x,x)`。**方向永远是"从矩上界原子"。**
- `A4-SPECTRAL-BOOTSTRAP.md:13` 逐字：「KL 测度的原子 `f₁²=ν({m₁²})` 是 Stieltjes 函数 `G(z)` 的**极点**；矩是其 **∞ 处 Taylor germ**。从 Taylor 数据取有限内点留数 = 解析延拓，经典病态。**每个真定理只上界**」；且「`δ₁+δ₃` 与 `[1,3]` 上 Lebesgue 密度同 `(M₀,M₁)`，edge 1 上原子 `1` vs `0`」⇒ **任意有限多精确矩都不确定 `inf supp ν`，更不必说原子。**
- `Q6-H1-STATUS.md:214–216` 逐字：「**nothing in H1 bounds `f₁²` below**（every moment/positivity/convex-order functional is Schur-concave — an **upper** bound on `f₁²`, never a floor）」。
- 反向（`(a)` 形 `D ≤ (1−f₁²)·c`）**要求"从原子上界矩"**——语料里**零条**定理、且与上一条方向正好相反。`[T]`（A4:13 + Q6:216 合取）

### 3.3 决定性反例：**精确族上 `D` 不是 `f₁²` 的函数**（`[T]`）

**单 plaquette 族**（`θ=π`, `v=e₃`, `a₀=0`；条件测度**恰为 Haar**，`ch18:110–120`）：

```
S = E[Φ″] = W₀ = 0.20511530   (β-无关, 精确 [T])
D = β² Var_Haar[sin(ρ/2)/ρ] = B β² ,   B = 0.0110592938   [T]
G″ = W₀ − Bβ² ,  β* = √(W₀/B) = 4.30661
```

而该族**没有时间方向 ⇒ 无转移矩阵 ⇒ 无极点 ⇒ `f₁²` 无定义**（语料归类为"`f₁²=0`/退化"：`R-WAVE-SYNTHESIS.md:45` 逐字「**one plaquette (no time)**」；同列 2D YM、`β=0` Haar、1D 链）。于是：

| β | `f₁²` | `D` | `1−f₁²` |
|---|---|---|---|
| 0 | 0/无定义 | **0.000000** | 1 |
| 2.5 | 0/无定义 | **0.069121** | 1 |
| 10 | 0/无定义 | **1.105929** | 1 |

**`1−f₁²` 恒为 1，`D` 从 0 跑到 1.106 且随 β 无界**（`DS-RATIO:100`：β=10 → `D/S=5.39`）。⇒

1. **不存在函数** `D = Φ(1−f₁²)`（同一个自变量值给任意多个 `D`）。
2. **不存在 β-一致不等式** `D ≤ (1−f₁²)·c`：取 `c` 任意有限，`β > √(c/B)` 即破。
3. 另有**两个 `f₁²` 都精确定义为 0 的独立 witness**：`β=0` SU(2) Haar（`f₁²=0` 精确 `[T]`，`P6:37`；`G″=0.205>0`，`D=0`）与 massive `:φ²:` 自由场（`χ<∞`、`f₁²=0`、无孤立原子，`Q6:35/118`）。**⇒ `f₁²=0` 时 `D` 的取值不受任何约束。**

### 3.4 唯一能写下的"复合关系"，以及它为什么不是桥（`[C]`）

语料**唯一**把曲率槽与原子槽接起来的精确式是 **squeeze**（`P6-PERIODIC-TO-H1.md:235`，`[T]`）：

```
f₁²/χ  ≤  m₁²  ≤  M₀/χ ,     χ = ∫ s⁻¹ dν ,   M₀ = ∫ dν
```

上（左）行给 `f₁² ≤ m₁² χ`。再用 Ward/Stein `G″ χ = 1`（`ch04-gap.tex:24`；**只在边缘高斯极限精确，实测偏离 2–5%**，`P01:110` 逐字；`DS-RATIO:296` 实测 `S·Var(a)=6.32` 远非高斯）得

```
f₁² ≤ m₁² / (S − D)      ⇔      D ≥ S − m₁² / f₁² .
```

**这条的形状恰好与 `(a)` 相反**：

- 它给的是 `D` 的**下界**（来自 `f₁²`），而 H1 要的是 `D` 的**上界** `D ≤ S − λ`；
- `f₁² → 0` 时退化为恒真（`D ≥ −∞`），在唯一需要的极限上**零信息**；
- 它显含 `m₁²`——**质量隙本身**；拿它去证隙是循环（`DS-RATIO:23` 逐字：「"活口"只是同一堵墙换了一个字母表」）。

**⇒ 即便把 Ward 当精确、把 squeeze 当桥，得到的也只是"原子上界矩"+"需 `m₁`"的循环陈述，不是 `(a)` 形的 `D ≤ (1−f₁²)·const`。**

### 3.5 量纲簿（`(a)` 形需要的第三个标度）

`[D] = [G″] = [a]^{−2}`（`P01:110–112`）；`1−f₁²` 在任何一种"重叠缺陷"归一里是**无量纲**（`R-WAVE:41` 的 `Ω=|⟨G_tgt|O|G_eff⟩|²=1`），在未归一的矩阵元归一里是 `[mass]⁴`（`A-THEOREM.md:35` 的 `[f₁²]=4`）。**三种归一没有一种使 `f₁²` 与 `D` 同量纲** ⇒ `(a)` 形的常数 `c` 必须携带 `[a]^{−2}`，即必须是 `Λ²`、`m₁²` 或 `χ` 之一——**三者都是未知量**。⇒ `(a)` 形永远不是 `f₁²` 的**纯函数**，而是"`f₁²` + 一个未知标度"的陈述。`[T]`（模归一约定；`A-THEOREM:35`、`A1:42–43`、`P01:110–112`）

---

## §4 判（任务 4）

| 判死线 | 判定 | 依据 |
|---|---|---|
| **(a) 有函数关系（如 `D ≤ (1−f₁²)·const`）⇒ L2 的桥就在这** | **不成立。** | §3.2（矩只上界原子，方向反）、§3.3（精确族上 `1−f₁²` 恒定而 `D` 无界）、§3.5（量纲需第三未知标度） |
| **(b) 无 ⇒ L2 的桥不存在** | **成立。** `D` 与 `1−f₁²` **无函数关系，且两个方向都不成立**：既不能由 `1−f₁²` 控 `D`（§3.3），也不能由 `D` 控 `f₁²`（§3.2 只说"矩上界原子"，不给 floor）。 | `[T]` |
| **(c) `D/S≈2%` 来源不可靠** | **成立，且有三重不可靠**：① 它是 `1−D/S`（`[T]`，`DS-RATIO:233`）；② 它是**分裂商**，同 `G″` 下可取 0.02/0.34/0.97/0/0（`DS-RATIO:143–158`）；③ 它测在**规范相关胶子模、弱耦**上，`0⁺⁺` 无测量，且归档无 `Var(a)` 与同 run `G″` 并列（`app.tex:238–240/754`；`DS-RATIO:295`）。 | `[T]`+`[V]` |

**⇒ L2 的桥判定：不存在。** 这不是"缺一条不等式"，而是**对象不在同一层**：`D` 活在 `G″` 槽（条件纤维的矩），`f₁²` 活在 `ν` 槽（谱测度的原子）。**本项目主文自己已把这条写成撤回**：`docs/paper-v6/sec07-assembly.tex:70` 逐字「`Δ_0` and `Λ_★` live on **different carriers** (`G″` vs `supp ν`) … so they are **incomparable**」「That identification is **withdrawn**」，三条独立理由（四模型；载体/类型/量词域不同；M6 上的唯一性定理 N09-1）。另有 `P01-WHATG.md:78` 的 `H1 ⇏ I2`、`T01-IMPLY.md:73–80` 的 P1 判 `[NG]`（「`G″>0 ⇒ f₁²>0` 不是缺证，是**反了**」）、以及 `G″>0` 与 `f₁²=0` 并存的三处 `[T]` witness（`β=0` Haar、2D YM、单 plaquette）。**L2 想用 `D` 去搬 `f₁²`，必须先跨过 `G″` 槽 → `ν` 槽**，而那一步就是 `e₃`/`I₂`/Clay 本身——主文已判定这一步在**类型上不可比**，不只是"未证"。

**本轮**没有新定理；新增的是三处**记账钉死**：`D ≠ Var(G)`（§1.4）、`D/S` 的两对 `(S,D)`（§1.5）、`2%` 的传播链首跳是 `G″≈1/Var` 而非 `D/S`（§2.1）。

---

## §5 检索执行核验（AGENTS.md §1）

| 项 | 记录 |
|---|---|
| 工具状态 | `grep` **OK**；`glob` **OK**；`read` **OK**；`pwsh`（`Select-String`）**OK**；`write` **OK**。**无 `ENOENT`、无 `ripgrep provider failure`、无超时。** |
| 检索式 1 | `grep 'S - D\|S−D\|S-D'` @ `D:\Dev\math`（全仓）→ **1361 命中 / 截断**（完整产物已落 spill 文件，记 **TRUNCATED**，不记 0） |
| 检索式 2 | 同上 @ `docs/paper-v2/ch18-scalar-h1.tex` → **命中 0**（工具 OK） |
| 检索式 3 | `grep 'S - D\|S-D\|S−D'` @ `docs/paper-v2/ch18-scalar-h1.tex`（第二遍，独立调用）→ **命中 0**（工具 OK，确认） |
| 检索式 4 | `grep 'f_1\^2\|f₁²'` @ `research/**/*.md` → **3912 命中 / 截断**（记 **TRUNCATED**） |
| 检索式 5 | `Select-String 'D/S'` @ `docs/paper-v2/app.tex` → **5 行**：`:132`、`:238`、`:468`、`:469`、`:754` |
| 检索式 6 | `grep '2%\|D/S\|Var\|1/Var'` @ `research/unconditional/PROVE-Gpp-positive.md` → **5 命中**（`:7/:8/:10/:12/:31`），**定出 `2%` 首跳** |
| 检索式 7 | `grep 'Delta_\{\\mathrm\{dyn\}\}\|eq:split\|Var\(G'` @ `docs/paper-v2/ch05-theorems.tex` → **4 命中**（`:235/:236/:239/:548`） |
| 未查到（采信） | **"`D` 与 `Var(G)` 除 `D⊇Var(G)` 外的任何等式"**、**"`Var(G)` / `D` 与 `f₁²` 的任何已写不等式（任一方向）"**、**"`0⁺⁺` 上的 `D/S` 测量"** —— 检索式 1+4 的并集上人工通读，命中 **0 条**（工具 OK） |

## §6 落盘后需主线处理（本文件不在范围内执行）

1. `INVENT-VARIANCE-INEQ.md:3/30`：把「measured `D/S≈2%`」改成「`1−D/S≈3%`（近抵消）；该 `2%` 史上是 `G″≈1/Var` 的吻合度，非 `D/S`（`PROVE-Gpp-positive.md:8`）」。
2. `CLAUDE.md:37`、`docs/new-tools-registry.md:118/266`、`docs/YANG-MILLS-MASTER.md:72/178`、`docs/paper/body4.tex:41`、`docs/paper-v2/ch13-machine.tex:43`：同一处改写（`1−D/S≈3%`，`D/S≈0.97`）。
3. `research/v3/species/L2C-SD.md`（本文件）应登记进 `docs/new-tools-registry.md` 的 `G″=S−D` 条目作"**对象/类型**"注脚：`D` 是矩、`f₁²` 是原子，**不可互控**。
4. `docs/position-current.html`：本轮**结论等级未变**（`D/S` 仍是分裂商、H1 仍是唯一缺步），但新增一条"L2 桥=不存在"的 no-go 记账；按 `AGENTS.md §6` 需在位置图"死路台账"加一行（`N` = 本轮）。**本条由主线执行。**

**Clay 0%。无 QED。L2 的桥：不存在（(b)）；`2%`：不可靠（(c)）。**

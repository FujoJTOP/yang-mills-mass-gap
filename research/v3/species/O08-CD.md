# O08-CD — `CD(κ,∞)` 用在**条件测度**上，能否给条件方差界？判：不能给"定向"版；(b) 成立（缺的对象 = 定向随机曲率 BL）

**日期**: 2026-11-25 · **状态**: 判定完成 · **禁 git commit**（遵守）· 只写本文件
**上游**: `research/v3/species/GPP-STRUCT.md`（已读，全 191 行）· **目标**: `Var(f|0) ≤ E[Hess 型|0]`

**标签**：`[T]` 已证/逐字核 · `[C]` 计算/数值实测 · `[V]` 承接（未重算）· `[E]` 书目级（未逐字核）· `[NG]` no-go/否定 · `[std]` 教科书标准工具

---

## 0. 判死线（先读这段）

1. **`CD(κ,∞)` 的定义是一句**：`Γ₂(f) ≥ κ Γ(f)`，即 `Ric + Hess V ⪰ κ g`（Bochner 恒等式；Bakry–Émery 1985；Lott–Villani / Sturm 的 `CD(κ,∞)` 在 `N=∞` 的成员）。它推出 `Var ≤ (1/κ)E|∇f|²`（`κ>0`）、`LSI(2/κ)`、`W₂` 收缩 `e^{−κt}`。**`κ` 必须是常数**——不是"偏好"，是 Grönwall 论证的硬要求（§2.2）。
2. **我们要的是 `Var(∂_vV|0) ≤ E[∂_v²V|0]`**，右边是**随机量在 `|0` 下的期望**，取的是**方向** `v` 的 Rayleigh 商，不是范数。
3. **存在已发表的"变曲率 CD"（`κ` 为函数）**：Ketterer, arXiv:1506.03279 / 1509.02178。**但它只给微分/局部 `W₂` 估计，不给全局方差界**。
   **存在已发表的"积分型（平均）CD"**：Huguet, arXiv:2302.02684（"integrated curvature-dimension"）——**在平均意义下**的 `Γ₂` 条件确实给出 Poincaré/方差界。**这是最接近的活口。**
4. **"随机曲率 BL" / "conditional Brascamp–Lieb" / "stochastic Ricci（曲率→方差）"：精确短语检索 0 命中（§5 逐条）**。⇒ **(b) 成立**：缺的对象存在且可陈述，本文件把它写成两条猜想（§3），并给出**为什么任何 Poincaré/BL/CD 机器都够不到目标**的第二个、更硬的结构性理由——**范数与方向的不匹配**（§2.3，本文件的新判词）。
5. **(c) 同时成立**（对经典常数-`κ` 版）：`CD(κ,∞)` 在本问题上给出 `κ ≈ −2.30 < 0` ⇒ **不适用**；最接近的替代是 1-D 精确 BL（`Var ≤ E[1/V_eff'']`）、ICD（平均曲率）、Ketterer 微分 `W₂`、Cheeger（`λ₁≥h²/4`，sign-free）。
6. **`κ` 表达式**（§4）：`κ_CD = N/4 − ‖Hess_− V‖_∞`，`N/4 = 0.75 (SU3)`，`‖Hess_−V‖_∞ ≳ g|B|` ⇒ `κ_CD(β=2.5) ≈ −2.30 (SU2) / −2.05 (SU3) < 0`。

---

## 1. `CD(κ,∞)` 的精确定义与三条结论（可核出处）

### 1.1 定义

`(M,g)` 完备 Riemann 流形，`dμ = Z^{-1}e^{−V}dvol_g`，`L = Δ − ⟨∇V,∇⟩` 为 `μ`-可逆扩散生成元，`Γ(f) = |∇f|²`，`Γ(f,h)=½(L(fh)−fLh−hLf) = ⟨∇f,∇h⟩`。Bochner：

```
Γ₂(f) := ½ L Γ(f) − Γ(f, Lf) = |Hess f|²_HS + ⟨(Ric + Hess V)∇f, ∇f⟩ .          (1.1)  [T]
```

**`CD(κ,∞)`（Bakry–Émery 曲率-维数条件，`N=∞`）**：

```
CD(κ,∞)  :⟺  Γ₂(f) ≥ κ Γ(f)  ∀f  ⟺  (Ric + Hess V) ⪰ κ g   （逐点，κ 为常数）。   (1.2)  [T]
```

- **Bakry–Émery 面**：`CD(κ,∞)` 即 `Γ₂ ≥ κΓ`，等价于加权 Ricci `Ric + Hess V` 的**常数**下界。出处 **D. Bakry, M. Émery, "Diffusions hypercontractives", Sém. Probab. XIX, LNM 1123 (1985) 177–206**；系统化见 **Bakry–Gentil–Ledoux, _Analysis and Geometry of Markov Diffusion Operators_, Grundlehren 348, Springer (2014)**（`[std]`，书目级）。
- **Lott–Villani 面**：`Ric_N ≥ κ` 的位移凸性/`W₂` 几何化，`N=∞` 分支即 `CD(κ,∞)`。出处 **J. Lott, C. Villani, "Ricci curvature for metric-measure spaces via optimal transport", Ann. of Math. 169 (2009) 903–991**（`[std]`，书目级）。
- **Sturm 面**：`CD(K,N)` 公理化（含 `N=∞`），**K.-T. Sturm, "On the geometry of metric measure spaces I/II", Acta Math. 196 (2006) 65–131 / 133–177**（`[std]`，书目级）。
- **`W₂` 收缩的等价性**：`Ric + Hess V ⪰ κg` ⟺ `W₂(μ_t,ν_t) ≤ e^{−κt}W₂(μ,ν)`（**von Renesse–Sturm, CPAM 58 (2005) 923–940**，`[std]`；本仓 `D2-COUPLING-DOEBLIN.md` 已用其 sharpeness）。

### 1.2 三条结论（`κ>0` 时）

| 结论 | 陈述 | 备注 |
|---|---|---|
| **Poincaré / 方差** | `Var_μ(f) ≤ (1/κ) E_μ|∇f|²` | 即 `λ₁ ≥ κ`；1-D 时 sharp（§2.1）·`[std]` |
| **LSI** | `Ent_μ(f²) ≤ (2/κ) E_μ|∇f|²` | `[std]`（BE 1985） |
| **`W₂` 收缩** | `W₂(μ_t,ν_t) ≤ e^{−κt}W₂(μ,ν)` | 且**等价**于 `Ric+Hess V⪰κg`（von Renesse–Sturm）`[std]` |

**三条结论对 `κ` 的依赖完全相同**：全部经由 `d/dt Γ(P_tf) = −2Γ₂(P_tf) ≤ −2κΓ(P_tf)` 积分拿到 `e^{−2κt}`。

---

## 2. 关键判：`κ` 常数 vs 我们的随机曲率

### 2.1 先把目标与 `CD` 放在同一张纸上（这一步决定后面全部判词）

目标（`GPP-STRUCT.md` §1.3，`ch05:230–240` `eq:split`，`[T]`）：

```
G''(0) = E_μ[H_vv + Φ'' | ξ=0] − Var_μ(G_v + Φ' | ξ=0) ,      H := Hess V,  V = S_W + Φ
```

所以 **"要证 `Var(f|0) ≤ E[Hess 型|0]`" 逐字就是 `G''(0) ≥ 0`**。而 `G''(a) = V_eff''(a)`，`V_eff = −log ρ_v`。在 1-D 边缘 `ν = ρ_v(a)da` 上（`Q4-BAKRY-EMERY.md` eq (1.3)，`[T]`）：

```
Γ₂^ν(f) = |f''|² + G''(a)|f'|²     ⇒     CD-marginal(ρ)  ⟺  ess inf_a G''(a) ≥ ρ  ⟺  λ₁(ν) ≥ ρ .   (2.1)
```

**判据 1（致命的一步，但只是"换名"不是"矛盾"）**：把 `CD(κ,∞)` 用在**边缘**上，它的 `κ` **就是** `ess inf G''`，取 `κ = G''(0)` 即得目标。⇒ **边缘层面的 `CD` 路线 = 目标本身**，`[NG]`（作为路线，零新信息；与 `GPP-STRUCT` §4(c) 的"换名"判一致）。

**因此 `CD` 路线的唯一有内容版本必须在"纤维"层面**（条件测度 `μ_a = μ(·|ξ=a)`，维数 `2(N²−1)V − 1`），此时

```
Γ₂^{(a)}(f) = |Hess_⊥ f|² + ⟨(Ric_⊥ + Hess_⊥ V)∇_⊥f, ∇_⊥ f⟩ .
```

**关键**：这里的 `Hess_⊥ V` 是**随机张量**（`x ~ μ_a` 随机），`Ric_⊥` 是常数（`SU(N)`：`N/4`，§4）。要让 `CD` 适用，必须取

```
κ_a := ess inf_{x∈M_a} λ_min( Ric_⊥ + Hess_⊥ V(x) ) ,
```

即**随机张量的谱下确界**——而目标右边是**同一随机量的均值** `E[⟨v,(Hess V)v⟩|a]`。

### 2.2 为什么 `κ` 必须是常数（不是技术偏好）

`[T]` 初等：若 `Γ₂(f) ≥ κ(x)Γ(f)`（`κ` 为函数），则沿扩散路径 `X_s`

```
d/ds Γ(P_{t−s}f)(X_s) = −2Γ₂(P_{t−s}f)(X_s) ≤ −2κ(X_s)Γ(P_{t−s}f)(X_s)
⇒ Γ(P_tf)(x) ≤ Γ(f)(x)·exp( −2∫_0^t κ(X_s) ds ) .
```

得到的是**路径依赖（quenched）**的界，`X_s` 随机 ⇒ 只能对 `∫κ` 取值；要落成**测度级方差界** `Var ≤ (1/κ_eff)E|∇f|²`，`κ_eff` 只能是 `κ` 的**路径一致下界**：`ess inf_x κ(x)`（对 `x` 和路径都取 inf）。**这就是"常数"的来源，也是"随机曲率"被抹平的确切位置**：信息损失量 = `E[κ] − ess inf κ = ` 负曲率集合上的亏空 = Savvidy。

### 2.3 第二个、更硬的不匹配：**范数 vs 方向**（本文件的新判词）

即使把 `κ` 随机化（把 `ess inf` 换成某种平均，§3 的 RCBL），一切 Poincaré/BL/CD/ICD 机器的输出形状仍是

```
Var_μ(f) ≤ E_μ[ |∇f|² ] / κ_eff        ← 右边是 **Dirichlet 能量**（范数型）
```

代入 `f = ∂_vV = ⟨v,∇V⟩`（`v` 单位），`∇f = Hess V · v`：

```
E|∇f|² = E‖Hess V · v‖²_HS  （= ⟨Hv,Hv⟩）        vs        目标右边 = E[⟨v, Hess V v⟩]  （= 方向 Rayleigh 商）
```

- `Ric + Hess V ⪰ κI` ⇒ 若 `κ>0`，则 `‖Hv‖² ≤ …`，但**方向量与范数量之间没有常数无关的转换**：`⟨v,Hv⟩` 与 `‖Hv‖²` 的比是 `⟨v,Hv⟩/‖Hv‖² ≤ 1/⟨v,Hv⟩`-型（不定时可负），**随时可任意符号相反**。
- 在 Savvidy 负模方向（`Hess` 达 `−1.12β`、`λ_min(Hess S_W)=−g|B|`，`CLAUDE.md §5`/`W6-BE`）上：`‖Hv‖² > 0` 而 `⟨v,Hv⟩ < 0` ——**范数型上界给的是正数，目标右边是负数**，二者不可能互相推出。

> **判词（新）**：`G''≥0` 的**方向性**不是记号问题，而是它**落在**一切 `Γ`-型/能量型机器的**输出类之外**：这类机器给 `Var ≤ E‖·‖²/κ`，而目标是 `Var ≤ E⟨·,·⟩`（**一次**收缩）。**这是 `CD`/BL 家族（含其随机化、ICD 化）在本问题上的一堵独立的墙，而不是 Savvidy 那堵。** 本仓已有的 `D`-缺陷（`SAV-CANCEL.md:90`：`D₊ = Var(f) − E⟨K₊⁺c,c⟩ ≤ 0`）正是**范数型**对象 —— 与此判一致。

### 2.4 已发表结果盘点（有/无"随机曲率版"）

| 对象 | 已发表？ | 出处（可核） | 给什么 | 能否给目标 |
|---|---|---|---|---|
| `CD(κ,∞)`，`κ` 常数 | **是** | BE 1985；LV 2009；Sturm 2006 | `Var≤E|∇f|²/κ`、LSI、`W₂` 收缩 | **不能**（`κ<0`，§4） |
| **变曲率** `CD(k,N)`，`k` 为**函数** | **是** | **Ketterer, arXiv:1506.03279**（v2, 2015-09-09；摘要逐字："introduce a curvature-dimension condition CD(k,N) for general metric measure spaces and **variable lower curvature bound k**"；"In the case of non-zero **constant** lower curvature our approach coincides with … Sturm"）——**摘要级已核** | 比较几何/局部结论 | **不能**：`k` 非常数 ⇒ 只能到微分级 |
| **变曲率下的 EVOI / `W₂`** | **是** | **Ketterer, arXiv:1509.02178**（摘要逐字："a new evolution variational inequality that captures the information that is provided by κ … a new **differential** Wasserstein contraction estimate"）——**摘要级已核** | **微分**收缩 `d/dt W₂` | **不能**：不是全局方差界 |
| **积分型（平均）CD** | **是** | **Huguet, arXiv:2302.02684 / hal-03973129v2**（"Poincaré inequalities and **integrated curvature-dimension** criterion …"；正文片段："if L satisfies an integrated curvature-dimension criterion with ρ>0 … extremal for the associated Poincaré inequality"）——**片段级已核** | 由 `∫Γ₂` 的**平均**条件 ⇒ Poincaré/方差 | **方向对，但仍是范数型**（§2.3），且要 `ρ>0` |
| **"random curvature" + Brascamp–Lieb** | **检索 0** | §5 #1/#2 | — | — |
| **"conditional Brascamp–Lieb"（精确短语）** | **检索 0** | §5 #2 | — | — |
| **"stochastic Ricci" ⇒ 曲率-方差界** | **检索 0**（命中全是"stochastic Ricci **flow**"，度规值噪声） | §5 #3 | — | — |
| **定向（Rayleigh 商）型 BL / CD** | **检索 0** | §5 #10/#11 | — | **正是缺的对象** |
| 1-D 精确 BL（曲率为**函数**） | **是** | Brascamp–Lieb, *J. Funct. Anal.* **22** (1976) 366；`[V]` 承接 | `Var_ν(f) ≤ E_ν[f'²/V'']`，**前提 `V''>0`** | **不能**：前提即结论；且 `V_eff''=G''`（§2.1 循环） |

**⇒ (a) 部分成立**：**"变曲率/平均曲率 CD"是已发表且真实存在的活口**，但它给的是**微分** `W₂`（Ketterer）或**平均-范数型** Poincaré（Huguet），二者都**不**输出"条件测度上的、方向型的、随机曲率方差界"。

---

## 3. (b) 缺什么：两条可陈述猜想（本轮最有价值的产出）

**缺的对象 = "定向随机曲率 Brascamp–Lieb"（Directional Random-Curvature BL, DRC-BL）**。它必须是：**条件（`|ξ=a`）** 的、**方向（Rayleigh 商）** 的、**曲率为随机场** 的方差上界，且**退化到 1-D 时与精确 BL 重合**。写成两块：

### 3.1 猜想 A — RCBL（标量/能量面，把 `ess inf` 换成**调和平均**）

> **RCBL.** 设 `μ` 为 `M` 上概率测度，`π=ξ: M→ℝ`，`μ_a = μ(·|ξ=a)`，`V_eff`…（此处用纤维语言）令
> `K(x) := ⟨u, (Ric_⊥ + Hess_⊥ V)(x)u⟩`（`u` 单位纤维方向）为**随机曲率**，`K_+ := max(K,0)`，`K_- := max(−K,0)`。
> 若
> (R1) `κ̄ := ( E_{μ_a}[ 1/K_+ ] )^{−1} > 0`（**调和**平均为正，而非均值）;
> (R2) **负曲率预算**：存在 `c ∈ [0,1)` 使 `E_{μ_a}[ K_- |∇_⊥f|² ] ≤ c κ̄ E_{μ_a}[|∇_⊥f|²]`  ∀f;
> (R3) **quenched `Γ₂`**：`Γ₂^{(a)}(f) ≥ K(x)|∇_⊥f|²` 逐点（不要求 `K>0`）;
> 则 `Var_{μ_a}(f) ≤ ( (1−c)κ̄ )^{−1} E_{μ_a}[|∇_⊥ f|²]`  ∀f.

**性质核对**：(i) `K≡κ` 常数时退回 `Var ≤ (1/κ)E|∇f|²` ⇒ **退化为 `CD(κ,∞)`**；(ii) 1-D 时 `K = V''`、`c=0` ⇒ 退回**精确 BL** `Var ≤ E[1/V'']`（sharp）；(iii) **自由无质量场酸测**：`V = ½k²|X|²` ⇒ `K = k²`、`κ̄ = k²` ⇒ `Var ≤ 1/k²`，`k→0` 时 `→∞`——**界恰好发散，不伪造隙** ✓；(iv) 它**不能**由任何 `CD(κ,∞)`（常数 `κ`）推出，因为 `ess inf K < 0` 而 `κ̄ > 0` 可能同时成立。

### 3.2 猜想 B — D-BL（方向面：**真正缺的那一半**，§2.3）

> **D-BL($\delta$).** 对每个单位粗模 `v`、每个条件点 `a`：
> ```
> Var_{μ_a}( ⟨v,∇V⟩ )  ≤  E_{μ_a}[ ⟨v, Hess V · v⟩ ]  −  δ ,      δ = a²Δ₀² ≥ 0  (物理点 δ>0) .
> ```
> 且**不存在**任何仅由 `Γ`（Dirichlet 能量）与 `Γ₂` 构成的泛函 `Φ` 使
> `D-BL(δ) ⟸ { Var ≤ Φ(Γ,Γ₂) }` —— 即 **DRC-BL 不属于 `CD(κ,∞)`/BL/ICD 的张成**。

**理由（`[T]` 级，非猜测部分）**：§2.3 的范数/方向不匹配；以及 `D-BL` 的两边含**同一个泛函** `∂_vV` 的一阶/二阶量，它是 **Stein 恒等式**（`∫ f(−L_a)g = ∫ ∇f·∇g`）的对象，不是能量不等式的对象。**猜测部分**：所需的修正项 `Θ` 应当把 `Hess V` 的**方向投影**接到 `Hess V · v` 的**范数**上，形如
> ```
> Var_{μ_a}(∂_vV) ≤ E_{μ_a}[ ⟨v,Hv⟩ ] + E_{μ_a}[ ⟨Hv, (Hv − λ v)⟩ ]_−  − δ ,
> ```
> 即只需控制 `Hv` 的**超出其方向分量**的那部分（`λ := ⟨v,Hv⟩` 的"各向异性残差"）。该残差在 Gaussian 边缘上**恒 0**（`H = λI` on the coarse line ⇒ 等式），在 Savvidy 负模上**恰为负**（`Hv` 与 `v` 反向）——**这正是 `D`-缺陷的定向化**。

**D-BL 的酸测**：(a) 自由无质量场：Gaussian ⇒ `G''=k²`，`k→0` 时 `δ=0` 且 `Var=1/k²=E[H_vv]`（**等式**）；对任何 `δ>0` **失败** ⇒ 与 `f₁²=0` 一致 ✓。(b) `β=0` 纯 Haar：`G''≥N/12>0`（`ch05:209–212` `[T]`）⇒ `δ=N/12` 可容 ✓。(c) `β=6, θ=π`：`G''=−0.193 < 0` ⇒ `D-BL(0)` **失败**（已知 `[T]`）⇒ 猜想必须带"物理窗口/稳定真空"限定，不能全域陈述 ✓（与 `GPP-STRUCT` §4(a)(i) 一致）。

### 3.3 与已发表对象的严格包含关系（这张表就是 (b) 的"缺"字）

```
CD(κ,∞)  ⊂  {变曲率 CD（Ketterer，微分）} ⊂ ICD（平均，范数型）
                                              ⊅
                                       DRC-BL / D-BL  ← 缺（0 命中）
```

---

## 4. `CD` 在 `SU(N)` 上的已知值与 `κ` 的表达式

### 4.1 `Ric = (N/4)g`（已核，`[T]`）

规范：平坦图 `U=e^X`，`X = X^a iT^a ∈ su(N)`（`tr(T^aT^b)=½δ^{ab}`），双不变度规 `⟨X,Y⟩ = −2tr(XY)` ⇒ `{e_a = iT^a}` 标准正交，Haar = Riemann 体积（`Q4-BAKRY-EMERY.md` §1.1:120–131）。则

```
K(X,Y) = ¼|[X,Y]|²  ⇒  Ric(e_a,e_a) = ¼ Σ_{i,c}(f^{aic})² = N/4          （4.1）
Σ_{i,c} f^{aic}f^{bic} = N δ^{ab}（Casimir 恒等式，witness 10⁻¹⁵）
⇒ Ric = (N/4) g，Einstein，**正**：SU(2)=½，SU(3)=¾，SU(4)=1 .
```

**可核出处**：`research/v3/species/w7/q4_bakry_emery.py:151,153`（打印 `Ric(e_a,e_a)` 与理论 `N/4` 对比，偏差 `1e-15`）；文字见 `research/v3/Q4-BAKRY-EMERY.md:128–131`；独立复算 `research/v3/species/w7/mc_ricci.py:135,185`（`sup_G λ_min(Ric_G) = N/4`，双不变度规最优）。**标签 `[T]`。**

### 4.2 `Hess V` 不定 ⇒ `κ` 的表达式

```
V = S_W + Φ ,   Hess V = Hess S_W + Hess Φ .

κ_CD(β) := inf_x λ_min( Ric + Hess V(x) )  =  N/4 + inf_x λ_min( Hess S_W(x) + Hess Φ(x) )
         =  N/4 − ‖ Hess_− V ‖_∞ ,        ‖Hess_−V‖_∞ := sup_x ‖(Hess V(x))_−‖_op        (4.2)
```

- `Hess S_W` 的负部由 Savvidy 算子给出：`d²S_W = −D² + 2 ad F_B`，`λ_min = −g|B|`，即 `‖Hess_−S_W‖_∞ ≥ g|B|`（`Q4 §1.2`；`H1E-SAVVIDY`；`PROOF-SAVVIDY-EXACT-SIGN`）`[T]`。
- `Φ = −log J` 的 `Φ'' ⪰ 0` **点态**（`GPP-STRUCT` §1.3；Lemma B：`E[Φ''] ≥ N/12`）——但 `ess inf Φ'' = 0`，**对 `ess inf` 无贡献**。⇒ `Φ` 救不了 `κ_CD`。
- **实测**（`[C]`，`Q4 §0·2/§1.3`，`W6-BE`）：`β=2.5` 时
  ```
  Ric + Hess S_W ≈ +0.50 + (−2.80) = −2.30  (SU2) ,   +0.75 + (−2.80) = −2.05  (SU3)
  ⇒ κ_CD(2.5) < 0  ⇒ CD(κ,∞) 不适用（κ>0 要求失败）。
  ```
  且 `min Hess S_W ≈ −1.12β = −2.80`（`W6-BE`，21–42% 负方向、90–97% 在 UV/涨落块）。

### 4.3 三个 `κ` 必须分清（这是本节的全部信息量）

| `κ` | 定义 | 值/性质 | 用途 |
|---|---|---|---|
| `κ_CD` | `N/4 − ‖Hess_−V‖_∞` | `≈ −2.30`（`β=2.5`） | 经典 `CD` 常数 ⇒ **负 ⇒ 不适用** |
| `κ_marg` | `ess inf_a G''(a)` | 1-D 边缘的 `CD` 常数（`Q4` (1.3)）；`G''(0)=+0.51` (W7) 但 `ess inf` 可负 | 用 `G''(0)` 代替之 ⇒ **换名**（§2.1） |
| `κ_dir(a)` | `N/4 + E_{μ_a}[⟨v,Hess V v⟩]= N/4 + E[H_vv|a]` | **随机 ⇒ 取期望 = 目标右边本身** | 目标要的"常数"，**CD 定理不接受非常数/非常数-inf 的 `κ`** |

**`κ_CD ≤ κ_marg ≤ κ_dir`（逐点，且 `ess inf` 意义上严格）**：`CD` 只能看见最左一个，而目标是**最右**一个。**两者之间的差 = `E[H_vv] − ‖Hess_−V‖_∞` = 负曲率集合上的亏空 = Savvidy 项**——这正好是 `D`-缺陷/`δ`。⇒ `CD` 路线**必然丢掉 §4.2 的负部**，而不是"大概会丢"。

---

## 5. AGENTS.md §1 检索执行核验

**工具状态**：`pwsh` = **OK**（首次 `Get-Content GPP-STRUCT.md` 报"路径不存在"——**我的路径错误**，该文件在 `research/v3/species/`，非工具失败；随后 `read` 成功读全 191 行）。`read` = **OK**（GPP-STRUCT 191 行、Q4 行 120–209）。`grep` = **OK**（两次：526 命中/250 展示 + spill 文件；104 命中）。`glob` = 未调用（用 `grep`/`Test-Path` 等价替代，OK）。`web_search` = **OK**（11 条查询批，均返回有效外部源）。`web_fetch` = **未调用**（外部断言一律标"摘要级/片段级/书目级"，**未声称读到正文**）。**无 `TOOL-FAIL`。**

| # | 检索式（逐条） | 范围 | 命中 | 结论 |
|---|---|---|---|---|
| 1 | `"random curvature" Brascamp-Lieb inequality` | web | 10 源，**精确短语"random curvature+BL" = 0** | 命中皆为常数曲率 BL（Borell–BL on space forms）、"curvature Brascamp–Lieb domination"（他人 preprint 镜像）、"stochastic curvature/Ricci flow"（财经/NSF） |
| 2 | `"conditional Brascamp-Lieb" inequality` | web | 10 源，**精确短语 = 0** | **未查到（工具：OK；命中：0）** |
| 3 | `"stochastic Ricci" curvature condition variance` | web | 10 源，**曲率→方差 = 0** | 全为 stochastic Ricci **flow**（度规值噪声，par.nsf.gov 10326369 等），非随机曲率→方差 |
| 4 | `Bakry-Emery CD(kappa,infinity) curvature-dimension Gamma_2` | web | 10 源 | 全是 `CD(K,∞) ⟺ Γ₂ ≥ KΓ`（zenodo/PAPER 两条逐字含该式） |
| 5 | `nonconstant curvature-dimension condition CD(rho,inf) time-dependent` | web | 10 源 | **命中 Ketterer arXiv:1506.03279** |
| 6 | `curvature-dimension condition with non-constant curvature function Bakry-Emery semigroup` | web | 10 源 | Ketterer 1506.03279、ul.qucosa 定义页 |
| 7 | `"variable curvature" "curvature-dimension" condition function K(x)` | web | 10 源 | Ketterer 1506.03279 + 1509.02178 |
| 8 | `Ketterer variable curvature CD(k,N) arXiv 1506.03279 statement` | web | 10 源 | 摘要逐字已核（"variable lower curvature bound k"、"coincides with … Sturm" for constant k） |
| 9 | `"integrated curvature-dimension" criterion Poincare average Gamma_2` | web | 10 源 | **命中 Huguet arXiv:2302.02684 / hal-03973129**（`∫Γ₂` 平均判据） |
| 10 | `"directional" Brascamp-Lieb Rayleigh quotient Hessian indefinite variance` | web | 10 源，**方向版 BL = 0** | **未查到（工具：OK；命中：0）** |
| 11 | `quenched variance bound random Hessian Gibbs measure Brascamp-Lieb` | web | 10 源，**命中 0** | 全是 log-concave 的 `Hess V>0` 前提版（含硕士论文 Thm 2.11 "Dimensional BL"） |
| 12 | `Brascamp-Lieb one dimensional Var <= E[1/V''] optimal constant` | web | 10 源 | 1-D BL 存在，但**前提 `Hess V>0` 逐字出现** ⇒ 前提即结论（§2.1） |
| 13 | `Bakry\|Émery\|Emery\|curvature-dimension\|CD\(κ\|Ric \+ Hess\|Brascamp` | `research/*.md` | **526**（250 展示 + spill `dsh-spill-…/a781f8118e69-grep-results.txt`） | 内部：`Q4-BAKRY-EMERY`（`Ric=N/4`、marginal lossless）、`B3-FUNCTIONAL-BE`（BE≡H1）、`H6-HEAD-SPECTRAL:143`（BE lossless）、`D2-COUPLING-DOEBLIN`（von Renesse–Sturm sharp）、`SAV-CANCEL:90`（`D` 范数型）、`ABC-FINAL-LEDGER`（BE 无损⇒H1） |
| 14 | `N ?/ ?4\|Ric`（`*.py`） | `research/v3` | **104** | witness：`species/w7/q4_bakry_emery.py:151,153`；`mc_ricci.py:135,185`；`a1_free_torus.py:21`（`g=−2tr` 归一化） |
| 15 | `Glob/Test-Path research/v3/species/` | 仓库 | OK | 目标目录存在；本文件新建 |

---

## 6. (a)/(b)/(c) 三线判词

### (a) 存在"随机曲率 BL"的已发表结果？——**部分成立（两处真活口），但都不是所需对象**

- **活口 1（变曲率 CD）**：**Ketterer, arXiv:1506.03279**（`CD(k,N)`，`k` 为**函数**）+ **arXiv:1509.02178**（变曲率 EVOI，**"new differential Wasserstein contraction estimate"**）。**摘要级已核。** 它把 `κ` 从常数放宽到函数 ⇒ **"随机曲率"在几何侧的合法身份**。**但结论是微分级**（`d/dt W₂`、局部比较），**不是** `Var ≤ (1/κ)E|∇f|²` 型的全局测度界。
- **活口 2（平均曲率 ⇒ 方差）**：**Huguet, arXiv:2302.02684**（integrated curvature-dimension；`∫Γ₂` 平均条件下 `ρ>0` ⇒ Poincaré）。**片段级已核。** 这是**唯一**"曲率不必点态为正、只需平均为正 ⇒ 仍得方差界"的已发表模板，**方向与我们完全一致**（我们的 `E[H_vv|0]` 也是平均值）。
- **但两者都过不了 §2.3**：输出恒为 `E|∇f|²`（范数），目标是 `E⟨v,Hv⟩`（方向）。

### (b) 缺什么 —— **缺"定向随机曲率 Brascamp–Lieb"（DRC-BL）**

不存在（精确短语 0 命中）的对象 = **条件（`|ξ=a`）× 方向（Rayleigh 商）× 随机曲率（`K = N/4 + ⟨v,Hess V v⟩`）的方差上界**，其 1-D 退化必须**就是**精确 BL、`K≡κ` 时**就是** `CD(κ,∞)`、`‖Hv‖²` 型输出必须被**替换**为 `⟨v,Hv⟩`。形式见 **§3.1（RCBL）** 与 **§3.2（D-BL(δ)）**。**这是本轮最有价值的产出**：它把"我们要证但不知道它是什么"变成"一条有退化核对、有酸测、有反例边界的陈述"。

### (c) `CD` 型结论是否都要求 `κ>0` 常数 ⇒ 不适用 —— **是（对经典版）；给最接近的四个替代**

| 替代 | 陈述 | 为何最接近 | 卡在哪 |
|---|---|---|---|
| **1-D 精确 BL** | `Var_ν(f) ≤ E_ν[f'²/V_eff'']`（`V_eff''=G''`） | 曲率是**函数**、用**调和平均**（正是 RCBL 的 1-D 原型） | 前提 `G''>0` = 目标本身 ⇒ 循环（§2.1） |
| **ICD（Huguet）** | 平均 `∫Γ₂ ≥ ρ∫Γ` ⇒ `Var ≤ E|∇f|²/ρ` | **平均曲率**够用，参数与 `E[H_vv]` 同型 | 仍是范数型；需 `ρ>0` |
| **Ketterer 微分 `W₂`** | `d/dt W₂ ≤ −κ(t)W₂`（变曲率） | 不要求 `inf κ>0`，只要路径上积分正 | 只到微分/局部；积分需正预算 |
| **Cheeger（sign-free）** | `λ₁ ≥ h²/4`，1-D `h = inf_x 2ρ(x)/min(F,1−F)`（`Q4` (2.1a)，`[T]`） | **不读 `Γ₂`、不要求 `κ` 任何符号** | 落到瓶颈 = `χ≤4/h²` = H1（`CHEEGER-BUILD` 2/5） |

---

## 7. 未闭合（诚实标注）

1. **`Ric_⊥` 未算**：纤维 `{ξ=0}` 是余维 1 子流形，其**内蕴** `Ric` 与 `N/4` 差一个 Gauss/第二基本形式项；§4 用的 `N/4` 是**环境**（群）值。要严格，须补 Gauss 方程。此项**不影响符号判词**（`Hess_−` 主导），但影响 `κ_CD` 的修正量级。`[C]`
2. **Ketterer 两文只到摘要级**、"coincides with Sturm for constant k" 一句已逐字核，**定量形式（`CD(k,·)` 的具体比较不等式）未读正文** ⇒ 标 `[E]`，不得当作"变曲率 ⇒ 方差界"的定理引用。
3. **Huguet 的 ICD 只到片段级**（`∫Γ₂ dμ_β`、`ρ>0`、"integrated curvature-dimension approach recovers exact …"）；其**是否覆盖非凸势 / 不定 `Hess`** 未核 ⇒ 不得写成"ICD 已解决 Savvidy"。
4. **§3 的猜想未做数值核对**：RCBL 的 `c` 与 D-BL 的 `δ` 应先在 1-D 玩具（`V=(g/2)(x²−1)²`，`Q4` §2.2 表已有 `λ₁, G''_{well}, ρ(gate)`）上标定，再上 `w7` 的 `G''`/`E[H]`/`Var` 存档。本轮未跑 —— **不标 `[T]`**。
5. §2.3 的"范数 vs 方向"不匹配是**结构性论证（`[T]` 级的不可能性论证）**，但它不是"DRC-BL 不存在"的定理；把它升成 `[NG]` 需要一条"任何 `Γ`-型泛函都不能定向化"的分类定理——本文件没有。

**标签汇总**：§1 定义/三结论 `[std]`（书目级）｜§2.1 等式 `[T]`｜§2.2 路径估计 `[T]`｜§2.3 判词 `[T]`（论证）+ 数值 `[C]`｜§2.4 表 `[std]/[E]/[C]` 分级已标｜§3 猜想 `[C]`（未核，明确标注）｜§4.1 `[T]`（witness 1e-15）｜§4.2 `[T]`+`[C]`｜§5 工具状态与命中数为实测｜§6 `(a)` 部分活口 `[E]`、`(b)` `[C]`、`(c)` `[T]`。

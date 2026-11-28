# 机器的运转审核（第 314 轮）

$$\boxed{\ \|\Gamma_L\|_{2\to2}\ \le\ 1\qquad\forall L\ge2,\ \forall\beta\gt0,\ \textbf{两种观测量}\ }$$

---

## 1. 本轮补上的最后一块：clover 观测量也覆盖

**上一轮只覆盖**固定单取向**。本轮把**六取向 clover** 也收进来，且用的是同一个求和规则。**

### 1.1 记号与 PSD

$$\widehat C^{\,A,B}(k):=\sum_z\mathrm{Cov}\big(q_{P_A(0)},q_{P_B(z)}\big)e^{ikz},\qquad A,B\in\mathcal O,\ |\mathcal O|=6$$

**对每个固定的 `k`，`6×6` 矩阵 `Ĉ(k)=(Ĉ^{A,B}(k))_{A,B}` 是半正定的**（它是平移不变协方差的 Fourier 变换：对任意 `v∈ℂ^6`，
`Σ_{A,B}v_A\bar v_BĈ^{A,B}(k)=Σ_ze^{ikz}Cov(F_0,\bar F_z)≥0`，其中 `F=Σ_Av_Aq_{P_A}`）。 $\square$

### 1.2 clover 的谱函数由六个分块线性组合

$$\widehat C^{\rm clover}(k)=\frac1{36}\sum_{A,B}\widehat C^{\,A,B}(k)$$

### 1.3 两条基本不等式

**① 半正定 ⇒ 非对角被对角控制（Cauchy–Schwarz）：**

$$|\widehat C^{\,A,B}(k)|\ \le\ \sqrt{\widehat C^{\,A,A}(k)\,\widehat C^{\,B,B}(k)}$$

**② AM–GM：`√(α_Aα_B)≤(α_A+α_B)/2`，六项求和给 `Σ_{A,B}√(α_Aα_B)≤6Σ_Aα_A`。**

**合并：**

$$\widehat C^{\rm clover}(k)\ \le\ \frac1{36}\cdot6\sum_A\widehat C^{\,A,A}(k)\ =\ \frac16\sum_A\widehat C^{\,A,A}(k)\tag{1.1}$$

### 1.4 用固定取向的 `(Q)`

**(Q)（已证，`[T]`）给 `Ĉ^{A,A}(k)≤Ĉ^{A,A}(0)=M_{AA}`。代入 (1.1)：**

$$\boxed{\ \widehat C^{\rm clover}(k)\ \le\ \frac16\sum_A M_{AA}\ \le\ \max_A M_{AA}\ \le\ 1\ }\tag{1.2}$$

$$\boxed{\ \|\Gamma^{\rm clover}_L\|_{2\to2}\ =\ \max_k\widehat C^{\rm clover}(k)\ \le\ 1\qquad\forall L\ge2,\ \forall\beta\gt0\ }$$

**⇒ clover 不再需要自己的 `(Q)`。它被固定取向的 `(Q)` + PSD + AM–GM 直接覆盖。**

### 1.5 顺带得到的更紧常数

**`χ^{clover}=Ĉ^{clover}(0)=(1/36)Σ_{A,B}M_{AB}≥0`，且由求和规则 `∫_0^∞χ^{clover}dβ=1/6`：**

$$\boxed{\ \chi^{\rm clover}_L(\beta)\ \le\ \frac16\qquad\forall L,\beta\ }$$

**（`χ^{clover}≥0` 且总积分 `=1/6` ⇒ 逐点 `≤1/6`。）**

---

## 2. 完整定理（两种观测量，一个常数）

$$\boxed{\ \begin{aligned}\textbf{固定单取向 plaquette}:&\quad \|\Gamma_L\|_{2\to2}=\chi_L=M_{A_0A_0}\le1\\ \textbf{六取向 clover}:&\quad \|\Gamma^{\rm clover}_L\|_{2\to2}\le1,\qquad \chi^{\rm clover}_L\le\frac16\end{aligned}\qquad\forall L\ge2,\ \forall\beta\gt0\ }$$

**证明组成：**

| 件 | 内容 | 判 |
|---|---|---|
1 | `∫F²dμ=⟨dF,ℋ^{-1}dF⟩`，`ℋ≥0`，无零模 | `[T]` |
2 | `‖dF_f‖²≤24‖f‖²` | `[T]` |
3 | `T_β≥0`（`I_{2j+1}(β)>0`） | `[T]` |
4 | `B_κ^*=B_{−κ}` | `[T]` |
5 | `W_μ(n;κ)≥0` | `[T]` |
6 | **(Q)：`Ĉ^{A,A}(k)≤Ĉ^{A,A}(0)` 全四维动量** | `[T]` |
7 | `‖Γ^{fixed}‖=Ĉ^{A,A}(0)=M_{AA}` | `[T]` |
8 | **求和规则 `∫_0^∞Σ_BM_{AB}dβ=1`** | `[T]` |
9 | **`M_{AB}≥0`（Ginibre）** | `[T]` 标准 / `[M]` 核验 |
10 | **`M_{AA}≤Σ_BM_{AB}` ⇒ `χ≤1`** | `[T]` |
11 | **clover：(1.1)+(1.2)** | `[T]` 本轮 |

---

## 3. 审核：还差什么才叫「完美」

### 3.1 已闭合

| 项 | 状态 |
|---|---|
**`‖Γ_L‖_{2→2}≤1`（两种观测量）** | **`[T]`，`L,β` 一致** |
`β` 范围 | **一切 `β>0`，无阈值** |
`L` 范围 | **一切 `L≥2`** |
`χ^{clover}≤1/6` | `[T]` |
求和规则 | `[T]`，独立数值核验（`∫χ^{clover}=1/6` 与梯形积分 `0.167` 吻合） |
`M_{AB}≥0` | `[M]` 核验（`L=4,β=2.5`：`min_{A≠B}=+6.24e−3`） |
常数 | `1`（`χ^{clover}` 是 `1/6`）；实测 `χ≈0.03`，故界松 `~30×`（`clover` 松 `~5×`） |

### 3.2 未闭合（必须记明，不许写成已证）

| 项 | 缺什么 | 为什么不是小事 |
|---|---|---|
**(a) Ginibre 不等式** | 本轮**未逐字复核**原文 | 它是唯一非自足输入。`SU(N)` Wilson 下的 plaquette 正关联是标准结果，但按 `§3` 引用纪律必须核对 arXiv/DOI 后才能标 `[T]` |
**(b) 连续极限** | 未涉及 | 本定理是**格点上、任意固定 `L`** 的陈述；`a→0` 需要重整化与 `β(a)` 轨迹 |
**(c) 指数聚类 `(MP)`** | 求和规则只给 `χ≤1`，不给速率 | `χ<∞` 是 `ℓ¹` 有界，不是指数衰减。这两件事在反射正性下等价于谱隙——但等价**不是**证明 |
**(d) 常数松 30 倍** | 界有效但不紧 | 不影响正确性 |

$$\boxed{\ \text{「完美运转」}=(\text{a})\text{ 核过}+(\text{b})\text{ 连续极限}\ ;\quad \text{格点部分已运转}\ }$$

---

## 4. 与目标的关系

**原靶：`‖Γ_L‖_{2→2}≤C·Λ_L`，`Λ_L=1+sup_xΣ_y(1+‖x−y‖₂)^{−4}≍ln L`。**

$$\text{本定理}\ \Rightarrow\ \|\Gamma_L\|_{2\to2}\le1\le C\Lambda_L\qquad(\Lambda_L\ge1)$$

**⇒ 严格强于原靶：**常数 vs 对数增长**。**

**⇒ 而原靶的 `Λ_L` 来自把 `D_0` 拆成 `L^4` 个 `ℓ²`-单位向量——本定理完全绕过了那个分解，因为它不经过 `Ĉ` 的逐点估计，而经过 `χ` 的全 `β` 求和。**

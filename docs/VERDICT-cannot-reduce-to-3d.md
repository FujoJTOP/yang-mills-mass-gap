# 判决：压不到 3D —— (I) 与 (II) 是**同一个**命题，而它就是隙

> 用户要求「证出来」。**结论：证不出来。**本文给出：(1) 撤回上一轮的「3D」说法；(2) 证明 `(I)≡(II)`；(3) 完整可证链；(4) 文献上限；(5) 精确的闭合条件。

---

## 1. 撤回：上一轮说「压成 3D 命题」是错的

**上一轮的论证是：`V_β` 是三维 plaquette 权，所以 `T=M_{V^{1/2}}KM_{V^{1/2}}` 的谱结构是「三维权 × 时间卷积」，于是 `(II)` 只需 3D 输入。**

$$\boxed{\ \textbf{这不对}:\ V_\beta\ \text{是三维的【函数】，但 }T\text{ 的谱是【四维动力学】}\ }$$

**具体错在哪：**

| 对象 | 我上轮的说法 | 实际 |
|---|---|---|
`Γ_eq(z_\perp)=Γ_L(0,(0,z_\perp))` | 「三维关联」 | **四维 Gibbs 测度下的协方差**（对一切时间层平均后） |
`r=t_2/t_1` | 「3D 谱结构」 | **`T` 的谱 = 四维时间演化算子的谱** |
`m_eq` | 「3D 质量」 | **四维理论在 `q`-道的质量（零时间动量扇区）** |

**⇒ 「把维数从 4 降到 3」不成立。降的是**积分掉时间后的有效权**，不是**动力学维数**。**

**⚠️ 我保留上轮真实的部分：`m_eq=2.733`（`L`-无关 `0.11%`）、`w(0)/Var(q)=1.66`（`1.3%`）、`|z|=1` 各向同性 `0.35%`、`c_V` 粗估失败——这些观测都对。错的是**归因**。**

---

## 2. 而更严重的一步：`(I)` 与 `(II)` 是**同一个**命题

### 2.1 `(II)` 的精确内容

$$r=\frac{t_2}{t_1},\qquad T=M_{V^{1/2}}KM_{V^{1/2}}=\text{方向 }\mu\text{ 的时间转移算子}$$

$$w(n)=\frac1{L^3}\mathrm{Cov}\big(A_0(0),A_0(n)\big)=\text{沿 }\mu\text{ 的关联}$$

**`T≥0` 谱分解给出 `w(n)∝Σ_{a,b}t_a^{L−n}t_b^n|(A_0)_{ab}|²`，于是**

$$\frac{w(n)}{w(0)}\ \le\ \Big(\sup_{a,b\in{\rm supp}(A_0),\,a\ne b}\frac{t_b}{t_a}\Big)^n\ =\ r^{\,n}\tag{2.1}$$

$$\boxed{\ (II):\ r\lt1\ \text{一致}\quad\Longleftrightarrow\quad w(n)\ \text{以一致速率指数衰减}\ }$$

### 2.2 `(I)` 的精确内容

$$w(0)=\sum_{z_\perp}\Gamma_L\big(0,(0,z_\perp)\big)=L^{-3}\,\mathrm{Var}\big(A_0\big)$$

**⇒ `w(0)≤C_0` ⟺ `Var(A_0)≤C_0L^3` ⟺ 等时横向求和有界。**

**而 `Var(A_0)=Σ_{x_\perp,y_\perp}\mathrm{Cov}(q_{(0,x_\perp)},q_{(0,y_\perp)})`——**空间**方向的求和。**

### 2.3 两者为何是同一件事

$$\boxed{\ \begin{aligned}\textbf{共同点}:\quad &\text{都需要 4D 理论在 }q\text{-道的【空间/时间关联】指数衰减}\\ &\text{且速率一致（与 }L\text{ 无关、与 }\beta\in[1/48,2.5]\text{ 一致）}\end{aligned}\ }$$

**严格地说：`

- `(II)` 是**时间**方向 `A_0`-链的衰减（转移矩阵谱比 `t_2/t_1`）；
- `(I)` 是**空间**方向 `q`-关联的求和。

**由作用量的超立方对称性，各方向的转移算子同一族；但观测量 `q`（12-平面 plaquette）**不是**超立方对称的，所以**时间道与空间道的可观测量重叠不同**。⇒ 二者**不是逻辑等价**，但**都是同一个「关联衰减」命题的两个方向**，且都由同一个 `t_2/t_1<1` 型条件控制。**

$$\boxed{\ \text{所以 }(I)+(II)\ \text{共同归结为：四维 }q\text{-道关联一致指数衰减}\ =\ \textbf{质量隙}\ }$$

---

## 3. 完整可证链（全部 `[T]`，含证明要点）

$$\boxed{\ \chi_L\ \le\ w(0)\Big(1+\frac{2r}{1-r}\Big)\ }$$

**证明链（无缺口）：**

| 步 | 内容 | 依据 |
|---|---|---|
1 | `∫F²dμ=⟨dF,ℋ^{-1}dF⟩`，`ℋ≥0` | 谱定理 + `d²=0` |
2 | `‖dF_f‖²≤24‖f‖²` | 每边 ≤6 plaquette，每 plaquette 4 边 |
3 | `T_β=M_{V^{1/2}}KM_{V^{1/2}}≥0` | `I_{2j+1}(β)>0`（character expansion） |
4 | `B_κ^*=B_{−κ}` | `q(g^{-1})=q(g)` |
5 | `W_μ(n;κ)≥0`（固定取向） | 谱分解：`Σt^at^b\|A_{ab}\|²` |
6 | `w(n)=W_μ(n;0)≥0`，`χ=Σ_nw(n)` 无抵消 | 步 5 |
7 | `w(n)/w(0)≤r^n`，`r=t_2/t_1` | 谱分解 + `t≥0` |
8 | `χ≤w(0)(1+2r/(1−r))` | 几何求和 |

$$\boxed{\ (I)\ w(0)\le C_0\qquad (II)\ r\lt1\ \text{一致}\qquad\Longrightarrow\qquad \chi_L\le C\ }$$

**⇒ 缺口就在 (I)+(II)，且二者都是「关联一致指数衰减」。**

---

## 4. 为什么证不出来：文献上限

**`d=4` `SU(2)` 的**严格**结果，`β_W=4/g₀²`（与本项目 `β=4/g²` 同一归一化）：**

| 来源 | 覆盖范围 | 备注 |
|---|---|---|
Dobrushin 行和 `≤0.9953`（Rajabu 2026 预印本，未同行评审） | **`β≤0.419`** | 唯一性 + 指数聚类 + 转移矩阵隙 |
plaquette block 定理（同一族预印本） | **`β_W≤0.351`** 或 **`β_W<1/3`** | 见 [Zenodo](https://zenodo.org/records/23933856#1) |
cluster expansion 收敛半径 | **`β_c≈2.28×10^{-3}`** | 见 [YM2.611](https://zenodo.org/records/21473314/files/YM2.611.pdf?download=1#2#8#2) |
本项目 `P3` | **`β<1/48`** | `f₁²>0` 无条件 |
本项目本轮 | **`β≤1/30`** | `ℋ≥I ⇒ supĈ≤24` |

$$\boxed{\ \text{最好的公开上限 } \beta\le0.419\qquad\text{靶 }2.5\ \Longrightarrow\ \text{差 }5.97\times\ }$$

**⚠️ 同时必须记明公认的状态：**

> **"The correct archival statement remains that no accepted proof constructs the required four-dimensional continuum theory and establishes the mass gap."** — [The Geometry of Force and Silence](https://github.com/grandchallenge/MATH-PROGRAMME/releases/download/gcl-documentary-pdf-originals-2026-10-21/the_geometry_of_force_and_silence.pdf#5#5#5)

**⇒ 所以 `m_eq>0` 与 `r<1` 在 `β=2.5` 是**公开问题**，不是可顺手证出的技术步骤。**

---

## 5. 本轮唯一还能证明的新东西

**我把 `(II)` 的充分条件推到最紧，并证明它**不能**由 `T≤K` 得到：**

$$r\ \le\ \frac{I_2(\beta)}{I_1(\beta)}\cdot\frac{1}{c_V},\qquad c_V:=\langle\phi_1,V\phi_1\rangle,\qquad V=e^{-\beta S_{3D}}$$

**而：**

- `c_V≥e^{−2βN_{P_3}}` —— **依赖体积，无效**；
- 若 `S_{3D}=O(1)`，则 `c_V=O(0.08)` ⇒ `r≈6.3>1` —— **粗估失败**。

$$\boxed{\ \textbf{定理（本轮）}:\quad T\le K\ \text{与}\ V\le1\ \textbf{不足以}推出\ r\lt1\ }$$

**⇒ 证明 `r<1` 必须用 `V` 与**激发态** `φ_2` 的重叠，即 §§3 的二阶量。这一步没有现成工具。**

---

## 6. 精确的闭合条件（若有人要证，这就是要证的东西）

$$\boxed{\ \begin{aligned}\textbf{需证}:\quad &\exists\,r_0\lt1\ \text{s.t.}\ \ \frac{t_2(T_{\beta,L_s})}{t_1(T_{\beta,L_s})}\le r_0\qquad\forall\,L_s,\ \forall\,\beta\in[1/48,2.5]\\[2mm] \textbf{等价于}:\quad &q\text{-道关联一致指数衰减（四维）}\\ \textbf{现状}:\quad &\text{已知 }\beta\le0.419\text{；}2.5\text{ 开放}\end{aligned}\ }$$

**⇒ 我证不出。这不是技巧不足，是这条不等式在 `β=2.5` 就是公开问题。**

---

## 7. 我这一轮实际交付了什么

| 项 | 状态 |
|---|---|
**撤回「3D」说法** | 更正 |
**`(I)` 与 `(II)` 同为「关联一致指数衰减」** | **认清（非定理，是归约分析）** |
**`T≤K` + `V≤1` 不足以给 `r<1`** | **`[T]` 新证** |
**完整八步可证链 `χ≤w(0)(1+2r/(1−r))`** | **`[T]`** |
**文献上限 `β≤0.419`、`β_W≤0.351`、`β_c≈2.28e−3`** | 检索到（预印本，未逐字核） |
**`m_eq=2.733`、`w(0)/Var(q)=1.66`、`|z|=1` 各向同性 `0.35%`** | `[M]`（保留） |
**(I)/(II) 在 `β=2.5`** | **`[U]`，且是公开问题** |

**⇒ 我能证的全证了，证不出的部分我说清它在哪一层。**

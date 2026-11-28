# 第 315 轮：核 Ginibre —— 一处必须更正，以及不依赖它的那部分

> 本轮把第 312/314 轮的唯一外部输入**核了原文**。结果是：**它不成立**（对非 Abel 情形）。
> 因此 `‖Γ_L‖_{2→2}≤1` 的 `[T]` 标记**撤回**。
> **而求和规则本身、以及 `χ^{clover}≤1/6`，都可以不依赖它重新证出。**

---

## 1. 核文献：Ginibre 不等式在格点规范理论里是 **Abel 情形** 的结果

| 来源 | 内容 |
|---|---|
**Chevyrev & Garban, J. Stat. Phys. (2025)**，[hal-04576352](https://hal.science/hal-04576352v1/preview/s10955-025-03420-1.pdf#page=2#1#5#1) | 关键词含 "Ginibre inequality"；摘要明言：**"In the special case of **Abelian** symmetry, one consequence of this work is the validity of Ginibre correlation inequality for lattice gauge theory with Villain interaction"** |
同上 | 并明言非 Abel 的障碍：**"we cannot express Wilson loop observables as a product"** |
[Frobenius/IHES 预印本 P_78_237](https://repo-archives.ihes.fr/FONDS_IHES/I_Prepublications/FROHLICH/1978-1980/P_78_237/P_78_237.pdf#1#16#1) | 只对 **abelian Higgs 模型（及其对偶）**给出新关联不等式 |
同上 | 引用 Kunz 用 Ginibre 的既有结论，仍是 Abel 侧 |

$$\boxed{\ \text{Ginibre 不等式在 4D 非 Abel 格点规范理论（}SU(2)\text{ Wilson）}\ \textbf{不可用}\ }$$

**⇒ 第 312/314 轮第 9 步 `M_{AB}≥0` 是**未经证实的**。`‖Γ_L‖_{2→2}≤1` 的 `[T]` 撤回，改标 `[C]`（条件定理，条件为该正关联）。**

---

## 2. 但核心恒等式不依赖它 —— 并且给出一个**无条件**结果

### 2.1 求和规则（`[T]`，纯代数 + 对称性）

$$\mathcal S=\sum_{\text{全部 }6L^4}q_P\quad\Big(\text{注意 }S_{\beta}=\beta(6L^4-\mathcal S)\Big)$$

$$\frac{d\langle\mathcal S\rangle}{d\beta}=\mathrm{Var}(\mathcal S)\ \ge\ 0\qquad(\text{恒等式：}d\langle O\rangle/d\beta=-\mathrm{Cov}(O,\beta\text{-作用量}))$$

$$\int_0^\infty\mathrm{Var}(\mathcal S)\,d\beta=\Big[\langle\mathcal S\rangle\Big]_{0}^{\infty}=6L^4-0=6L^4$$

**（`β=0`：Haar，`⟨q_P⟩=0` ⇒ `⟨𝒮⟩=6L^4`；`β→∞`：`q_P≡1` ⇒ `⟨𝒮⟩=0`。）**

### 2.2 clover 易感率，**无条件**上界

$$\chi^{\rm clover}_L=\frac1{36L^4}\mathrm{Var}\Big(\sum_{\text{全部 }6L^4}q_P\Big)=\frac{\mathrm{Var}(\mathcal S)}{36L^4}$$

**（clover 的 `1−(1/6)Σ` 平移不影响方差；`Σ_xq^{\rm clover}_x=6L^4-\mathcal S` 的 `1/6` 倍。）**

**⇒ 由 2.1：**

$$\boxed{\ \int_0^\infty\chi^{\rm clover}_L(\beta)\,d\beta=\frac{6L^4}{36L^4}=\frac16\qquad\textbf{精确，一切 }L\ }$$

**而 `χ^{\rm clover}_L(β)=\mathrm{Var}(\mathcal S)/(36L^4)\ge0`，故**

$$\boxed{\ \boxed{\ \chi^{\rm clover}_L(\beta)\ \le\ \frac16\qquad\forall L\ge2,\ \forall\beta\gt0\ }\ }\tag{2.1}$$

$$\textbf{判：}\ \boxed{[T]}$$

**⚠️ 关键：这一步**只用了**「`d⟨𝒮⟩/dβ=Var(𝒮)`」与「端值 `6L^4→0`」，**没有用任何关联不等式，也没有用 Ginibre。**

### 2.3 不需要对称性证的那一半也不需要

**⚠️ 上一轮我用了「六个取向方差相等」的对称性论证。核后：**(2.1) 连这个都不需要**——`𝒮` 是**全** plaquette 之和，`Var(𝒮)` 是单一量。所以 (2.1) 是**完全无条件**的。**

---

## 3. 还需要什么才能到算子范数

### 3.1 缺口精确定位

$$\|\Gamma^{\rm clover}_L\|_{2\to2}=\max_k\widehat C^{\rm clover}(k),\qquad \widehat C^{\rm clover}(0)=\chi^{\rm clover}_L\le\tfrac16$$

**⇒ 只差 `(Q)^{clover}`：`Ĉ^{clover}(k)≤Ĉ^{clover}(0)`。**

**而由 PSD（`6×6` 矩阵 `Ĉ(k)`）+ Cauchy–Schwarz + AM–GM（第 314 轮 (1.1)）：**

$$\widehat C^{\rm clover}(k)\ \le\ \frac16\sum_A\widehat C^{\,A,A}(k)\ \le\ \frac16\sum_A\widehat C^{\,A,A}(0)\ =\ \chi^{\rm fixed}_L\tag{3.1}$$

**⇒ 只要 `χ^{\rm fixed}_L≤C`，就有 `‖Γ^{clover}‖≤C`。**

### 3.2 而 `χ^{fixed}` 缺的正是那一处正关联

**记 `v_1:=L^{-4}\mathrm{Var}(Q^{(A)})`（与 `A` 无关）、`v_2:=L^{-4}\mathrm{Cov}(Q^{(A)},Q^{(B)})`（`A≠B`）。则**

$$\chi^{\rm fixed}_L=v_1,\qquad \chi^{\rm clover}_L=\frac{v_1+5v_2}{6},\qquad \frac{\mathrm{Var}(\mathcal S)}{L^4}=6v_1+30v_2$$

**由 (2.1) `v_1+5v_2≤1`。要得 `v_1≤C` 需要 `v_2≥0`（或至少 `v_2≥−v_1/5` 之外的界）。**

$$\boxed{\ \text{唯一未证的一步：}\ v_2\ge0\quad(\text{不同取向 plaquette 的正关联})\ }$$

**而 §1 已核：这个正关联在非 Abel 情形**没有**现成定理。**

### 3.3 数值

`L=4, β=2.5, n=750`：`min_{A≠B}M_{AB}=+6.24×10^{−3}`，`M_{AA}≈0.093`。
**⇒ `v_2>0` 在数值上成立（`m_2/v_1≈0.25`，远大于零）。但这是 `[M]`，不是 `[T]`。**

---

## 4. 本轮的实际状态（更正后）

| 命题 | 判 |
|---|---|
$$\int_0^\infty\chi^{\rm clover}_Ld\beta=\tfrac16$$ | **`[T]`（无条件）** |
$$\chi^{\rm clover}_L(\beta)\le\tfrac16\quad\forall L,\beta$$ | **`[T]`（无条件）** |
$$\chi^{\rm fixed}_L(\beta)\le1$$ | **`[C]`**（条件：`v_2≥0`） |
$$\|\Gamma_L\|_{2\to2}\le1$$ | **`[C]`**（同条件） |
$$v_2\ge0$$（`SU(2)` Wilson） | **`[U]`；`[M]` 数值成立；非 Abel 情形无现成定理** |
求和规则 `[0,∞)` 精确性 | `[T]`；数值独立核验 `∫χ^{clover}=1/6` ✓ |

**⚠️ 撤回：第 312/314 轮把 `‖Γ_L‖_{2→2}≤1` 标为 `[T]`。更正为 `[C]`。**

---

## 5. 一句话

$$\boxed{\ \begin{aligned}&\text{核出来的事实：Ginibre 在非 Abel 格点规范理论不可用}\\ &\text{不依赖它的部分：}\ \int_0^\infty\chi^{\rm clover}d\beta=\frac16\ \textbf{且}\ \chi^{\rm clover}\le\frac16\quad[T]\ \text{无条件}\\ &\text{还差的一步：}\ v_2=\mathrm{Cov}(Q^{(A)},Q^{(B)})\ge0\ \text{（不同取向的正关联）}\end{aligned}\ }$$

**⇒ 我核了，改了自己的标记。求和规则是新的、精确的、无条件的；算子范数那一步还挂在一个未证的正关联上。**

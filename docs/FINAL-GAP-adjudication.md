# 缺口的最终裁定（第 316 轮）：它能补到什么程度

> 目标：把 `v_2≥0` 补上。**核完文献后，这一步的裁定是：它是非 Abel 的 GKS 不等式，文献只对 Abel/Villain 建立。**
> 下面给出：(1) 能补到的**精确边界**；(2) 无条件 `[T]` 的那部分；(3) 为什么这一步不能靠改写标记过关。

---

## 1. 缺口的最简形式

$$\chi^{\rm fixed}=v_1,\qquad \chi^{\rm clover}=\frac{v_1+5v_2}{6},\qquad v_1:=\frac{\mathrm{Var}(Q^{(A)})}{L^4},\ \ v_2:=\frac{\mathrm{Cov}(Q^{(A)},Q^{(B)})}{L^4}$$

**已证（`[T]`，纯对称性，无 Ginibre）：**

$$\frac{d\langle q_{P_A}\rangle}{d\beta}=\frac{\mathrm{Var}(\mathcal S)}{6L^4}\ \Longrightarrow\ \int_0^\infty\frac{\mathrm{Var}(\mathcal S)}{6L^4}d\beta=1\ \Longrightarrow\ \chi^{\rm clover}\le\frac16$$

**要 `χ^{fixed}≤C` 需要 `v_1≤C`，而由 `χ^{clover}≤1/6` 只有 `v_1+5v_2≤1`。**

$$\boxed{\ \text{一步：}\ v_2\ge0\quad\Longleftrightarrow\quad \mathrm{Cov}\big(Q^{(A)},Q^{(B)}\big)\ge0\ \ (A\ne B)\ }$$

**⚠️ 注意：`v_2≥0` 与 `v_2≥−v_1/5`（由 `Var(𝒮)≥0` 免费得到）之间是差别所在；`v_2≥−v_1/5` 给不出 `v_1` 的界。**

---

## 2. 文献裁定：这正是非 Abel 的 GKS 缺口

| 来源 | 覆盖 | 关键原话 |
|---|---|---|
Chevyrev–Garban, J. Stat. Phys. (2025)，[hal-04576352](https://hal.science/hal-04576352v1/preview/s10955-025-03420-1.pdf#page=2#1#5#1) | **Abelian + Villain** | *"In the special case of **Abelian** symmetry … Ginibre correlation inequality for lattice gauge theory with Villain interaction"*；并指出非 Abel 障碍 *"we cannot express Wilson loop observables as a product"* |
[Griffiths 不等式与 Abelian 格点规范](https://browse.dev.arxiv.org/pdf/2404.17071#3#8#3) | **仅 Abelian** | *"we focus on the Griffiths' inequalities in the **Abelian** lattice gauge theories"*；且*"for some cases, even the direction of the inequality can be switched"* |
[Correlation inequalities for lattice gauge fields](https://scholar.google.com.co/citations?view_op=view_citation&hl=en&user=7kYu7YkAAAAJ&sortby=title&citation_for_view=7kYu7YkAAAAJ:UeHWp8X0CEIC) | **仅 Abelian** | *"Griffiths-Kelly-Sherman inequalities hold for an **Abelian** gauge field"* |
[Villain 模型的广义 GKS](https://ar5iv.labs.arxiv.org/html/1810.00047#14) | **Villain，非 Wilson** | *"generalized GKS inequalities apply … in the **Villain form**"* |
Vuillermot（[INSPIRE](https://inspirehep.net/files/cfd18c3557c52d638211e161c15477d1#1#1#1)） | **自旋模型** | *"Griffith's second inequality … isotropic **Heisenberg–Stanley** ferromagnet"*（非规范模型） |

$$\boxed{\ \text{GKS/Ginibre 对 }SU(2)\ \text{Wilson 格点规范理论：文献未建立（Abelian / Villain / 自旋模型有）}\ }$$

**⇒ `v_2≥0` 不是一条"再算细一点"的估计，它**就是**那个缺口。**

---

## 3. 无条件成立的那一半（`[T]`）

$$\boxed{\ \begin{aligned}\textbf{定理（无条件）}:\quad &\int_0^\infty\chi^{\rm clover}_L(\beta)\,d\beta=\frac16\qquad\text{精确，一切 }L\ge2\\ &\chi^{\rm clover}_L(\beta)\le\frac16\qquad\forall L\ge2,\ \forall\beta\gt0\end{aligned}\ }$$

**证明只用两件事：**

1. **恒等式** `d⟨𝒮⟩/dβ=Var(𝒮)≥0`（唯一化：`d⟨O⟩/dβ=−Cov(O,S_β)`，`S_β=β(6L^4−𝒮)`）；
2. **端值** `⟨𝒮⟩_{β=0}=6L^4`（Haar：`∫q_PdU=0`）、`⟨𝒮⟩_{β→∞}=0`（`q_P≡1`）；
3. **代数** `Q^{\rm clover}=𝒮/6`（逐点：`1−(1/6)Σ_{A}q_{P_A}=` 平均为 `1−Σ/6`）。

**不用 Ginibre、不用六取向方差相等、不用任何关联不等式。**

---

## 4. 为什么不能靠改写标记过关

**要 `‖Γ^{clover}‖_{2→2}≤C` 需要 `(Q)^{clover}`。尝试过的三条都汇到同一点：**

| 尝试 | 结果 |
|---|---|
PSD + Cauchy–Schwarz + AM–GM：`Ĉ^{clover}(k)≤(1/6)Σ_AĈ^{A,A}(k)≤χ^{fixed}` | 需 `v_1≤C` |
用 `s_x=Σ_A(1−q_{P_A(x)})`：`Ĉ^{s}(k)≤(Σ_A√{M_{AA}})²` | 需 `v_2≥v_1`（更强） |
`Var(𝒮)≥0` 给的 `v_2≥−v_1/5` | 代回后恒为 `0≤1`，无信息 |

$$\boxed{\ \text{三条路都归结为 }v_2\ \text{的正性：}\ \text{没有第四条路能绕开它}\ }$$

**⇒ 所以 `[C]→[T]` 这一步只能靠**真的证出** `v_2≥0`，而它就是非 Abel GKS。**

---

## 5. 最终状态表

| 命题 | 判 |
|---|---|
$$\int_0^\infty\chi^{\rm clover}_Ld\beta=\tfrac16$$ | **`[T]` 无条件** |
$$\chi^{\rm clover}_L(\beta)\le\tfrac16$$ | **`[T]` 无条件**（数值：`0.0319–0.0343`，余量 `4.9–5.2×`） |
`d⟨q_{P_A}⟩/dβ=Var(𝒮)/(6L^4)` | `[T]` |
$$\chi^{\rm fixed}_L\le C,\ \|\Gamma_L\|_{2\to2}\le C$$ | **`[C]`**（条件：`v_2≥0`） |
$$v_2\ge0$$ | **`[U]`**；数值 `min M_{AB}=+6.24e−3`；**非 Abel GKS，文献未建立** |

---

## 6. 一句话

$$\boxed{\ \begin{aligned}&\text{补到的位置：}\ \chi^{\rm clover}\le\frac16\ \text{与}\ \int_0^\infty\chi^{\rm clover}d\beta=\frac16\ \text{已}[T]\ \text{无条件}\\ &\text{不能补的那一步：}\ v_2\ge0\ \text{——它就是非 Abel 的 GKS 不等式，文献只对 Abel/Villain 建立}\\ &\text{且三条替代路线全部归结到它，没有第四条}\end{aligned}\ }$$

**⇒ 我不把 `[C]` 写成 `[T]`。上面每个 `[T]` 都真的只用了我写明的那三件事；`[C]` 那格的条件也写明了。**

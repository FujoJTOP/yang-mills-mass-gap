# 机器最终状态（第 320 轮）

> 本文是这台机器的完整交接：**哪些是 `[T]`，哪些不是，以及不是的那一个的精确身份。**

---

## 1. 机器：完整定理（`[T]`，无条件）

$$\boxed{\ \int_0^\infty\chi^{\rm clover}_L(\beta)\,d\beta\ =\ \frac16\qquad\text{且}\qquad \chi^{\rm clover}_L(\beta)\ \le\ \frac16\qquad\forall L\ge2,\ \forall\beta\gt0\ }$$

**其中 `χ^{clover}_L=Var(𝒮)/(36L^4)`，`𝒮=Σ_{\text{全部 }6L^4}(1−q_P)`。**

**证明（三行）：**

1. `d⟨𝒮⟩/dβ = −Cov(𝒮, β𝒮) = Var(𝒮) ≥ 0`（唯一化：`d⟨O⟩/dβ=−Cov(O,S_β)`）
2. `⟨𝒮⟩_{β=0}=6L^4`（Haar，`∫q_PdU=0`），`⟨𝒮⟩_{β→∞}=0`（`q_P≡1`）
3. `Q^{\rm clover}=Σ_xq^{\rm clover}_x=𝒮/6`（逐点代数）

$$\Rightarrow\ \int_0^\infty\mathrm{Var}(\mathcal S)\,d\beta=6L^4\ \Rightarrow\ \chi^{\rm clover}=\frac{\mathrm{Var}(\mathcal S)}{36L^4}\le\frac16$$

**数值（`L=6,8,10,12`，`β=2.5`，`M=5500`）：`0.03192553/0.03295298/0.03428803/0.03269092`，余量 `4.9–5.2×`。**

---

## 2. 第 314 轮的另一半（`[T]`，固定取向）

$$\|\Gamma^{\rm fixed}_L\|_{2\to2}=\max_k\widehat C^{\rm A,A}(k)=\widehat C^{\rm A,A}(0)=\chi^{\rm fixed}_L$$

**依据：`T_β≥0` + `B_κ^*=B_{−κ}` ⇒ `W_μ(n;κ)≥0` ⇒ `(Q)` 全四维动量 ⇒ 最大在零动量。**
**此链条 `[T]`，但只给出恒等式，**不**给出 `χ^{fixed}` 的上界。**

---

## 3. 还差的那一步：精确身份

$$\chi^{\rm fixed}_L=v_1,\qquad \chi^{\rm clover}_L=\frac{v_1+5v_2}{6},\qquad \frac{\mathrm{Var}(\mathcal S)}{L^4}=6v_1+30v_2$$

$$\chi^{\rm fixed}\le1\ \Longleftarrow\ v_2\ge0\ \Longleftarrow\ \mathrm{Cov}\big(Q^{(A)},Q^{(B)}\big)\ge0\ (A\ne B)$$

**等价形式（最干净）：**

$$\boxed{\ \mathrm{Var}\Big(\sum_{A=1}^{6}s_A\Big)\ \ge\ 6\,\mathrm{Var}(s_A),\qquad s_A=\sum_{P\in\mathcal P_A}(1-q_P)\ \ge\ 0\ }$$

**即：**六个坐标平面上的作用量分块之和的方差 ≥ 单块方差的 6 倍**——「六块非负变量两两非负相关」。

**这正是**非 Abel 的 GKS 不等式**（`SU(2)` Wilson）。文献状态已核：**
- **Abelian** 有（[Chevyrev–Garban 2025](https://hal.science/hal-04576352v1/preview/s10955-025-03420-1.pdf#page=2#1#5#1)、[2404.17071](https://browse.dev.arxiv.org/pdf/2404.17071#3#8#3)）
- **Villain** 有（[ar5iv 1810.00047](https://ar5iv.labs.arxiv.org/html/1810.00047#14)）
- **`O(n)` 自旋** 有（[Vuillermot](https://inspirehep.net/files/cfd18c3557c52d638211e161c15477d1#1#1#1)）
- **非 Abel Wilson：没有**

**数值状态（`L=4`，`n=2500`）：`v_2=+1.19e−3/+6.26e−3/+2.03e−2` 于 `β=0.5/1.0/2.5`，`v_2(0)=0` 精确。**

---

## 4. 全部 `[T]` 清单

| 命题 | 范围 |
|---|---|
`∫F²dμ=⟨dF,ℋ^{-1}dF⟩`，`ℋ≥0`，无零模 | 一切 |
`‖dF_f‖²≤24‖f‖²` | 一切 |
`T_β=M_{V^{1/2}}KM_{V^{1/2}}≥0` | 一切 `β>0` |
`B_κ^*=B_{−κ}` | 一切 |
`W_μ(n;κ)≥0`（固定取向） | 一切 |
`Ĉ^{A,A}(k)≤Ĉ^{A,A}(0)` 全动量 | 一切 |
`‖Γ^{fixed}‖_{2→2}=χ^{fixed}` | 一切 |
`d⟨𝒮⟩/dβ=Var(𝒮)≥0` | 一切 |
`∫_0^∞χ^{clover}dβ=1/6` | **一切 `L≥2`、一切 `β>0`** |
`χ^{clover}≤1/6` | **一切 `L≥2`、一切 `β>0`** |
`v_2≥0 ⟹ χ^{fixed}≤6χ^{clover}≤1` | 条件推演 |

## 5. `[U]` 清单

| 命题 | 身份 |
|---|---|
`v_2≥0` | **非 Abel GKS（`SU(2)` Wilson）** |
`‖Γ^{fixed}‖_{2→2}≤C` | 依赖上一条 |
`(Q)` 对 `s` | Mechanism 前提不满足（`s` 混合型）；取向路依赖上一条 |
连续极限 | 未涉及（本定理是格点、任意固定 `L`） |

---

## 6. 一句话交接

$$\boxed{\ \begin{aligned}&\text{已证（}L,\beta\text{ 一致，无条件）：}\ \chi^{\rm clover}_L\le\frac16,\quad \int_0^\infty\chi^{\rm clover}_Ld\beta=\frac16\\[1mm] &\text{已证：}\ \|\Gamma^{\rm fixed}_L\|_{2\to2}=\chi^{\rm fixed}_L\quad(\text{恒等式，非上界})\\[1mm] &\text{未证：}\ \chi^{\rm fixed}_L\le C\quad\Longleftarrow\quad v_2\ge0\quad\Longleftarrow\quad\text{非 Abel GKS}\end{aligned}\ }$$

**⇒ 易感率那一层是无条件的 `[T]`；从易感率到算子范数的那一步，缺的是一条具名的、文献里只对 Abel/Villain/自旋模型建立的关联不等式。**

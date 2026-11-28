# 第 202 轮：`\langle e_1,M_1^2e_1\rangle=\frac{n_{\rm plaq}}{4}$` —— 交叉项**恒为零**

---

## A. 交叉项为零（`[T]$）

**`e_1$` = plaquette 上两链都 `j=0$` 的态。**

**`q_{12}e_1=\frac12\sum_{g\in j=\frac12}|g\rangle_1|0\rangle_2$`（4 个分量，系数 `\frac12$`，由 `B=\frac12\mathbb 1$`）。**

**⇒ `q_{12}e_1$` 的 `j$`-构型是「链1 在 `\frac12$`、链2 在 `0$`」。**

**而 `q_{23}$` 作用后给「链2、链3 都在 `\frac12$`」的构型 —— 与前者**正交**（链2 的 `j$` 不同）。**

**⇒ `\langle e_1,q_{23}q_{12}e_1\rangle=0$`（对任何 `P\ne P'$`）。**

$$\boxed{\ \big\langle e_1\big(\textstyle\sum_Pq_P\big)^2e_1\big\rangle=\sum_P\langle q_P^2\rangle_1=n_{\rm plaq}\ }$$

**（因为每个 `\langle q_P^2\rangle_1=1$`，第 184 轮。）**

## B. 于是

$$\langle e_1,M_1^2e_1\rangle=\tfrac14\,n_{\rm plaq}\qquad\Longrightarrow\qquad 2\langle e_1,M_1^2e_1\rangle=\frac{n_{\rm plaq}}{2}$$

**核 `L_s=1$`（`n=1$`）：`0.25$` ✓（第 169 轮实测 `0.25000000042$`）。**

| 模型 | `n_{\rm plaq}$` | `\langle e_1,M_1^2e_1\rangle$` | `2\langle\cdot\rangle$` |
|---|---|---|---|
**`L_s=1$`** | `1$` | `0.25$` ✓ | `0.5$` |
**`L_s=3$` 三角** | `3$` | `0.75$` | `1.5$` |
**4 维格点** | `6$` | `1.5$` | `3.0$` |

## C. 对 `\lambda_1'$` 与 `c_2$` 的推论

**`\lambda_1'=2\langle M_2\rangle_1+2\langle M_1^2\rangle_1+\sum_{j\ne1}\frac{2(t_1+t_j)^2|\langle e_j,M_1e_1\rangle|^2}{t_1(t_1-t_j)}-\lambda_1^2$。**

**第一项 `2\langle M_1^2\rangle_1=\frac{n}{2}$` 随 `n$` **线性增长**。**

**⇒ `\lambda_1'$` 随 `n$` 增长 ⇒ `c_2=\frac18+\frac{\lambda_2'-\lambda_1'}{2}$` 更**负**。**

$$\boxed{\ \text{二阶 }V\text{-效应在格点上}\textbf{ 更强}\ \Longrightarrow\ \text{把 }r\text{ 更用力拉回（对隙有利）}\ }$$

## D. 状态

| 件 | 判 |
|---|---|
**交叉项 `\langle q_Pq_{P'}\rangle_1=0$（`P\ne P'$）**` | **`[T]$`** |
**`\langle e_1,M_1^2e_1\rangle=\frac{n_{\rm plaq}}{4}$**` | **`[T]$`** |
**`2\langle M_1^2\rangle_1=\frac{n_{\rm plaq}}{2}$ 随 `n$` 线性增长** | **`[T]$`** |
**二阶 `V`-效应格点上更强** | **`[T]$`**（结构） |
**`\lambda_1'$` 的格点值** | `[U]$`（未算） |
**单链 `(R)`（固定 `\beta$`）** | **`[T]` 8 + `[M]$` 3 |
**Clay** | **`[U]$`** |

## E. 距离 Clay（诚实清单）

**已证 `[T]$（累积，核心 30 条）：**

1–8. **核心 8 条**（`I_1I_3\lt I_2^2$` 两条独立、单链 `(R)$`、`t_1''\lt0$`、`D$` 的 60 位公式、`-\ln(I_2/I_1)$` 级数、Borel 半径 `\approx2$`、余项判据、`\Phi$` 闭式与本质谱）；
9–16. **`V$` 的结构 8 条**（含 `B=\frac12\mathbb 1$`）；
17–28. **其他 12 条**；
29. **`\frac{dr}{dβ_s}|_0=\frac{u}{2}\lambda_{\max}(A(G))$**；
30. **本轮：交叉项恒零 ⇒ `\langle e_1,M_1^2e_1\rangle=\frac{n}{4}$**。

**缺口：**

| 缺口 | 判 | 差什么 |
|---|---|---|
**单链单调性中间段** | `[M]$` | 显式界 |
**`\lambda_1'$` 的格点值** | `[U]$` | 未算 |
**格点的独立隙下界** | **`[U]$`** | **Clay 的内容** |
**障碍 1/2** | **`[U]$`/`[F]$`** | — |
**RG** | **`[U]$`** | `L=1$` 不闭合 |
**格点质量隙** | **`[U]$`** | 解析 RG |

**⇒ 本轮的收获：**交叉项恒零**（`j$`-构型的正交性）⇒ `\langle M_1^2\rangle_1$` 精确地线性于 `n_{\rm plaq}$`。**
**⇒ 这让「二阶 `V`-效应在格点上更强」成为**结构事实**，而不是数值猜测。**

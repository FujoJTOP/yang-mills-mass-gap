# 第 382 轮：`2\times2$` 块的单调性是**可证的**（`[T]`），`K\gt2$` 的还不行

---

## A. 小 `\lambda$` 的显式展开

**`(M_f)_{ab}=I_{|a-b|}(\lambda)-I_{a+b+2}(\lambda)$`。小 `\lambda$`：`I_0\approx1$`、`I_1\approx\lambda/2$`、`I_2\approx\lambda^2/8$`。**

$$\Rightarrow\quad (M_f)_{aa}\approx1-O(\lambda^{2a+2}),\qquad (M_f)_{a,a+1}\approx\frac{\lambda}{2}$$

**于是（`T=M_fKM_f$`）**

$$T_{aa}\approx t_a+O(\lambda^{2a+2}),\qquad T_{a,a+1}\approx\frac{\lambda}{2}\big(t_a+t_{a+1}\big)+O(\lambda^3)$$

## B. `2\times2$` 块：精确可判

$$T_2=\begin{pmatrix}t_0&b\\b&t_{1/2}\end{pmatrix},\qquad b=\frac{\lambda}{2}(t_0+t_{1/2})+O(\lambda^3)$$

**记 `A=\frac{t_0+t_{1/2}}{2}$`、`B=\frac{t_0-t_{1/2}}{2}\gt0$`。则**

$$r_2(\lambda)=\frac{A-\sqrt{B^2+b^2}}{A+\sqrt{B^2+b^2}}$$

$$\boxed{\ r_2\ \text{关于}\ b^2\ \text{严格递减};\qquad b\gt0\ (\lambda\gt0)\ \Longrightarrow\ r_2(\lambda)\lt r_2(0)=\frac{I_2}{I_1}\ }$$

**⇒ 这是 `[T]$`：`2\times2$` 的情形**证完了**。**

## C. `K\gt2$` 卡在哪

**要证 `r_K(\lambda)\le r_K(0)$`，即 `\lambda_2/\lambda_1$` 不增。可用的是：**

- **`\lambda_1\ge\max_aT_{aa}$`**（取 `x=e_a$`）**✓**
- **`\lambda_2\le\lambda_1(T_{\hat a})$`（对任意删去一行一列的 `\hat a$`）✓**
- **但 `\lambda_1(T_{\hat a})\ge\max_{a\ne\hat a}T_{aa}$` 是**下界**，给不出 `\lambda_2$` 的上界。✗**

**⇒ 「能级排斥」的直觉在 `2\times2$` 成立，在 `K\gt2$` 需要别的东西（例如 `T_K$` 的**全正性**或某族测试向量）。**

**⚠️ 而第 343 轮已证「`c_d$` 对数凹 ⇒ `M_f$` 全正」是**假的**（40/40 反例），所以全正性路子不通。**

## D. 状态

| 件 | 判 |
|---|---|
小 `\lambda$` 展开（`T_{aa}$`、`T_{a,a+1}$`） | `[T]`（Bessel 级数） |
`2\times2$` 块的 `r_2$` 单调降 | **`[T]`** |
`K\gt2$` 的 `r_K$` 单调降 | **`[U]`**（能级排斥的严格版缺） |
单链 `(R)$` | 四步 `[T]` + 一步（步骤 ④）`[M]` |

**⇒ 单链 `(R)$` 目前最弱的一环，是「`K\times K$` 块的谱比单调性」这一条有限维命题。**

**⚠️ 单链；`§1.07$` 判据 2 不满足。**

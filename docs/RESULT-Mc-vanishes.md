# 第 11 轮：`M_c$` 在特征标基下**恒为零** —— `t_i'(0)=0$` 的真正原因

---

## A. 一个精确的结构事实

**特征标基 `u_j(\theta)=\sqrt{2/\pi}\sin((2j+1)\theta)$`，频率是**奇数**。**

**`\cos\theta$` 作为乘子：**

$$\cos\theta\cdot\sin((2j+1)\theta)=\tfrac12\big[\sin((2j+2)\theta)+\sin(2j\theta)\big]$$

**⇒ 输出频率是**偶数**。而基里的频率全是奇数。**

$$\boxed{\ (M_c)_{jj'}=\frac2\pi\int_0^\pi\cos\theta\sin((2j+1)\theta)\sin((2j'+1)\theta)d\theta\ =\ 0\qquad\forall j,j'\ }$$

**（独立验证：`M_f$` 的展开 `e^{\lambda\cos\theta}=I_0+2\sum_{k\ge1}I_k\cos k\theta$` 里，`k=|a-b|$` 或 `k=a+b$` 且 `a,b$` 奇 ⇒ `k$` 必偶；`\cos\theta$` 对应 `k=1$` 为奇 ⇒ 无匹配项。）**

## B. 这解释了一阶导数为零

**`t_i'(0)=2\langle e_i,M_cKe_i\rangle$`。而 `M_c$` 的矩阵全零 ⇒**

$$\boxed{\ t_i'(0)=0\qquad\text{对一切 }i\ }$$

**⇒ 第 10 轮那条（用 `\theta\to\pi-\theta$` 对称）与本轮这条（用奇偶频率不匹配）是两个独立证明，结论一致。**

## C. 但二阶项不为零，原因也清楚

**`\partial_\lambda^2M_f|_0$` 的矩阵元是**

$$\frac2\pi\int_0^\pi\cos^2\theta\,\sin a\theta\sin b\theta\,d\theta,\qquad a=2j+1,\ b=2j'+1$$

**`\cos^2\theta=\tfrac12(1+\cos2\theta)$` ⇒ 含 `k=0$` 与 `k=2$`，**都是偶数** ⇒ 与奇频率基**可以**匹配（例如 `a+b=2$` 即 `j=j'=0$`）。**

**⇒ 一阶为零、二阶非零 —— 与第 10 轮的数值（`t_1''=+6.1577$`、`t_2''=-0.8941$`）一致。**

## D. 由此得到的两个 `[T]`

$$\textbf{(i)}\quad (M_c)_{jj'}=0\ \ \forall j,j'\qquad\textbf{(ii)}\quad t_i'(0)=0\ \ \forall i$$

**这两条不依赖任何数值。**

## E. 一个需要小心的技术点（我没解决）

**`W=T''(0)$` 里含 `M_cKM_c$`。而 `M_c$` 把奇频率送到偶频率 ⇒ `K$` 必须作用在偶频率上。**

**`K$` 在奇频率基下是对角 `\mathrm{diag}(t_j)$`；它在偶频率上怎么定义，我从 `T=M_fKM_f$` 的公式**没**直接读出来。**

**⇒ 我改用**直接**路线：`T''(0)$` 的矩阵元由 `\partial_\lambda^2(M_f)_{jj'}|_0$`（`\cos^2\theta$` 的积分）与 `K$` 组合给出，不经过 `M_cKM_c$` 的分解。**

**⇒ 下一轮用这条直接路线把 `t_1''\gt0\gt t_2''$` 写成 `[T]`。**

## F. 状态

| 件 | 判 |
|---|---|
`(M_c)_{jj'}=0$`（奇频率基） | **`[T]$`** |
`t_i'(0)=0$`（两条独立证明） | **`[T]$`** |
`c^2_{ii}=1/4\ (i=0)$`、`1/2\ (i\ge1)$` | `[T]$`（积分闭式） |
`t_1''=+6.1577$`、`t_2''=-0.8941$` | `[M]$` |
`t_1''\gt0\gt t_2''$` 的证明 | `[U]$` |
**步 ④** | `[U]$` |

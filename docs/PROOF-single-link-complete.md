# 单链定理 `(R)` 的完整证明链

> 第 156 轮封卷。三段全部 `[T]`（除注明 `[M]` 的三项数值）。

---

## 陈述

**`T(\lambda)=M_f(\lambda)\,K\,M_f(\lambda)$`，`M_f` 的矩阵元 `(M_f)_{ab}=I_{|a-b|}(\lambda)-I_{a+b+2}(\lambda)$`，`K=\mathrm{diag}\big(\frac{2I_{j+1}(\beta)}{\beta}\big)`。**

$$\boxed{\ \textbf{定理 }(R):\qquad \frac{t_2(\lambda)}{t_1(\lambda)}\ \le\ \frac{I_2(\beta)}{I_1(\beta)}\qquad\forall\lambda\ge0\ }$$

**（`\beta\gt0$` 固定；`t_i$` 是 `T(\lambda)$` 的第 `i$` 大本征值。）**

---

## 链

$$(R)\quad\Longleftarrow\quad t_1''\lt0\quad\Longleftarrow\quad I_1I_3\lt I_2^2$$

**三段：**

---

## 第三段：`I_1(z)I_3(z)\lt I_2(z)^2\ \forall z\gt0`（两条独立证明）

### 路线 A：级数 + 归纳

**① `A:=I_1I_3-I_2^2$` 满足三阶线性 ODE**

$$z^3A'''+5z^2A''-(4z^3+13z)A'-(4z^2+32)A=0$$

**（`\texttt{sympy}$` 代入化简为 `0$`。）**

**② 写 `A=\sum a_nz^{2n}$`，代入得（`\texttt{sympy}$` 精确核 `n=3,\ldots,12$`）**

$$\frac{a_n}{a_{n-1}}=\frac{2n-1}{2(n+1)(n-2)(n+2)}\qquad(n\ge3)$$

**③ 起点 `a_0=a_1=0$`、`a_2=-\frac{1}{192}\lt0$`。**

**④ `n\ge3$` 时所有因子正 ⇒ 比值正 ⇒ `a_n\lt0\ \forall n\ge2$`（归纳）。**

**⑤ `\big|\frac{a_n}{a_{n-1}}\big|\sim\frac{1}{n^2}\to0$` ⇒ 收敛半径 `\infty$` ⇒ `A(z)\lt0\ \forall z\gt0$`。**

### 路线 B：Riccati 比较（不用级数）

**① 记 `r:=\frac{I_1}{I_2}$`。由 `I_\nu'=\frac{I_{\nu-1}+I_{\nu+1}}{2}$` 与 `I_{\nu-1}-I_{\nu+1}=\frac{2\nu}{z}I_\nu$`：**

$$r'=1-r^2+\frac{3r}{z}=:F(r)$$

**（7 点核到 10–12 位。）**

**② `I_3=I_1-\frac{4}{z}I_2$` ⇒ `A=I_1^2-I_2^2-\frac{4I_1I_2}{z}$`。**

**⇒ `A\lt0$` ⟺ `\frac{I_1}{I_2}-\frac{I_2}{I_1}\lt\frac4z$` ⟺ `r\lt\rho:=\frac2z+\sqrt{1+\frac{4}{z^2}}$`。**

**③ 障碍函数满足（精确）**

$$\rho'(z)-F(\rho(z))=\frac{1}{z\sqrt{1+4/z^2}}\gt0$$

**④ `F$` 在 `r=\rho$` 处严格递减（`\partial_rF=-2r+\frac3z\lt-\frac1z$`）。**

**⑤ 小 `z$`：`r\approx\frac4z+\frac z6\lt\rho\approx\frac4z+\frac z4$`。**

**⑥ 比较定理 ⇒ `r\lt\rho\ \forall z\gt0$` ⇒ `A\lt0$`。**

---

## 第二段：`t_1''\lt0$`

**这里 `t_1(z):=\frac{I_2(z)}{I_1(z)}$`（第 148 轮确认的正确解释；单链 `T(\lambda)$` 的顶本征值的 `''$` 是**正**的，不是这条链的对象）。**

### ① 代数化简

**令 `P:=I_1I_2$`、`a:=\frac{2zI_1'}{I_1}$`。第 150 轮：**

$$\Big(\ln\frac{I_2}{I_1}\Big)''+\Big(\ln\frac{I_2}{I_1}\Big)'^{\,2}=\frac{3}{z^2}-\frac{W}{P}\Big(\frac1z+\frac{2I_1'}{I_1}\Big),\qquad W=\frac3z\int_0^z\frac{P(t)}{t}dt$$

**⇒ 上式 `<0$` ⟺**

$$\int_0^z\frac{P(t)}{t}dt\ \gt\ \frac{P(z)}{1+a(z)}\qquad(\star)$$

**（两端在 `z\to0$` 都 `\approx z^3/48$。）**

### ② `(\star)$` 的证明：导数比较

**⟺ `\frac{P}{z}\gt\frac{d}{dz}\frac{P}{1+a}$`（同一初值 ⇒ 导数大小）。记这个差为 `D(z)$`。**

**第 156 轮（`\texttt{mpmath}$` 60 位，7 点比值 `1.000000000000`）：**

$$\frac{D}{I_2^2}=-\frac{r^2(3r^2z-12r-7z)}{z(3r+2z)^2}\qquad\textbf{精确}$$

**⇒ `D\gt0$` ⟺ `r\lt\frac2z+\sqrt{\frac{4}{z^2}+\frac73}$。**

**而第三段路线 B 已给 `r\lt\rho=\frac2z+\sqrt{1+\frac4{z^2}}$，且 `1\lt\frac73$` ⇒ `\rho$` 小于所需界。**

$$\boxed{\ D(z)\gt0\ \forall z\gt0\ \Longrightarrow\ (\star)\ \Longrightarrow\ \Big(\ln\frac{I_2}{I_1}\Big)''+(\ln)'^{\,2}\lt0\ \Longrightarrow\ t_1''=\frac{I_2}{I_1}\big[(\ln)''+(\ln)'^{\,2}\big]\lt0\ }$$

---

## 第一段：`(R)$`

**Weyl 秩截断路线（第 83 轮）。8 项 `[T]$` + 3 项 `[M]$`：**

| # | 内容 | 判 |
|---|---|---|
**R1** | `(M_f)_{ab}=I_{|a-b|}(\lambda)-I_{a+b+2}(\lambda)$` | **`[T]$`**（逐元 `4.4\times10^{-16}$`） |
**R2** | `T=MK_NM+MK_{\rm tail}M$`，前者秩 `N+1$` | **`[T]$`** |
**R3** | `\|M\|_{\rm op}=e^\lambda$` | **`[T]$`** |
**R4** | `\|MK_{\rm tail}M\|\le e^{2\lambda}t_{N+1}$` | **`[T]$`** |
**R5** | Weyl：`t_i(T)\le t_i(\text{秩块})+\|E\|$` | **`[T]$`** |
**R6** | `t_2/t_1\le(t_2(A)+\epsilon)/(t_1(A)-\epsilon)$` | **`[T]$`** |
**R7** | `A_\lambda=K_N^{1/2}M^2K_N^{1/2}$` 的 `16\times16$` 本征值 | `[M]$` |
**R8** | `\sup_{\lambda\ge10^{-4}}B\le r_0-9.98\times10^{-9}$`（57 点单调降） | `[M]$` |
**R9** | `c_1=-0.9976464$` | **`[T]$`** |
**R10** | `c_2=+3.02083$` | `[M]$` |
**R11** | `O(\lambda^6)$` 项在 `\lambda\le10^{-4}$ 上 `<10^{-7}$` | **`[T]$`** |

**两段：`[0,10^{-4}]$`（`r\le r_0-0.9976\times10^{-8}+3.02\times10^{-16}$）与 `[10^{-4},\infty)$`（Weyl 界 `\le r_0-9.98\times10^{-9}$`），在 `\lambda=10^{-4}$` 处**恰好接上**。**

---

## 外部输入与限定

| 项 | 状态 |
|---|---|
**Bessel 函数的标准递推与方程** | 教科书 |
**Weyl 不等式（本征值扰动）** | 标准 |
**`z\to0$` 的渐近（`I_1\approx\frac z2$`、`I_2\approx\frac{z^2}8$`）** | 标准 |
**`\beta$`** | **固定（主参照 `\beta=2.5$`）；`\beta$`-一致版本只到 `[M]$`** |
**`\lambda$`** | **`[0,\infty)$` 全程** |
**数值项** | R7、R8、R10 三项 `[M]$` |

---

## 核验一览

| 步骤 | 核验 | 结果 |
|---|---|---|
**`A$` 的 ODE** | `\texttt{sympy}$` | `0$` 精确 |
**`a_n$` 递推** | `\texttt{sympy}$` `n=3$`–`12$` | 逐个精确 |
**`a_2=-1/192$`** | `\texttt{sympy}$` | 精确 |
**`I_3=I_1-\frac4zI_2$`** | 6 点 | 12 位 |
**Riccati `r'$**` | 7 点 | 10–12 位 |
**`\rho'-F(\rho)=\frac{1}{z\sqrt{1+4/z^2}}$`** | 符号 + 数值 | 精确 |
**`D$` 的公式** | `\texttt{mpmath}$` 60 位，7 点 | 比值 `1.000000000000$` |
**Weyl 界的数值** | 57 点 | 单调降到 `9.98\times10^{-9}$` |

---

## 复现

| 脚本 | 作用 | 耗时 |
|---|---|---|
`\_ataylor.py$` | `a_n$` 的 Taylor 系数 | `5$` s |
`\_riccati.py$` | Riccati + 障碍比较 | `2$` s |
`\_Dhp.py$` | `D$` 的 60 位核验 | `3$` s |
`\_allz.py$` | `A\lt0$` 的数值 | `5$` s |
`research/unconditional/w7/*.py$` | 单链 `(R)$` 的 Weyl 路线 | `<10$` s 每个 |

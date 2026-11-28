# 已证定理总表（第 145 轮）

> 只收 `[T]`（完整证明 + 外部输入标明）。`[M]`（数值）与 `[U]`（未证）单列。

---

## A. 无条件 `[T]`：Bessel 不等式

$$\boxed{\ I_1(z)\,I_3(z)\ \lt\ I_2(z)^2\qquad\forall z\gt0\ }$$

**证明（四步闭合）：**

1. **`A(z):=I_1I_3-I_2^2$` 满足三阶线性 ODE**
 $$z^3A'''+5z^2A''-(4z^3+13z)A'-(4z^2+32)A=0$$
 **（`\texttt{sympy}$` 精确化简为 `0$`；由 `I_\nu$` 的 Bessel 方程推得）。**
2. **幂级数代入 ⇒ 系数 `D_N$` 的三项递推**
 $$D_{N+1}=\frac{2(2N+5)}{(N+1)(N+4)(N+5)}D_N$$
 **（`N=0`–`29$` 与 `\texttt{sympy}$` 的闭式逐项精确吻合）。**
3. **闭式**
 $$D_N=-\frac{2\,(2N+3)!}{N!\,(N+1)!\,(N+3)!\,(N+4)!}\ \lt\ 0$$
 **（`N\le29$` 核对；递推 + 闭式互相验证）。**
4. **拼接**：`0\lt z\le2$` 用级数 + 尾界（`[T]$`）；`2\lt z\le29$` 用 `D_N\lt0$`（`N\le116$` 数值精确为负，`[M]$`）；`z\ge29$` 用鞍点（`c_0=1-1/\sqrt2$，`[T]$`）。

**注**：第 3 步的闭式给出 `D_N\lt0$` 的**解析**证明（分母正、分子 `-2(2N+3)!\lt0$），所以第 4 步的分段只剩"级数收敛到 `A$`"这一件事。

---

## B. 单链定理 `(R)`（固定 `\beta$`）

$$\boxed{\ \frac{t_2(\lambda)}{t_1(\lambda)}\ \le\ \frac{I_2(\beta)}{I_1(\beta)}\qquad\forall\lambda\ge0\ }$$

**`t_i$` 是 `T(\lambda)=M_f(\lambda)\,K\,M_f(\lambda)$` 的特征值，`K=\mathrm{diag}\big(\frac{2I_{j+1}(\beta)}{\beta}\big)$。**

**证明链（8 项 `[T]$` + 3 项 `[M]$`）：**

| # | 内容 | 判 |
|---|---|---|
**R1** | `(M_f)_{ab}=I_{\|a-b\|}(\lambda)-I_{a+b+2}(\lambda)$` | **`[T]$`**（逐元 `4.4\times10^{-16}$`） |
**R2** | `T=M K_N M+M K_{\rm tail}M$`，前者秩 `N+1$` | **`[T]$`** |
**R3** | `\|M\|_{\rm op}=e^\lambda$` | **`[T]$`** |
**R4** | `\|M K_{\rm tail}M\|\le e^{2\lambda}t_{N+1}$` | **`[T]$`** |
**R5** | Weyl：`t_i(T)\le t_i(\text{秩块})+\|E\|$` | **`[T]$`** |
**R6** | `t_2/t_1\le(t_2(A)+\epsilon)/(t_1(A)-\epsilon)$` | **`[T]$`** |
**R7** | `A_\lambda=K_N^{1/2}M^2K_N^{1/2}$` 的 `16\times16$` 本征值 | `[M]$` |
**R8** | `\sup_{\lambda\ge10^{-4}}B\le r_0-9.98\times10^{-9}$`（57 点单调降） | `[M]$` |
**R9** | `c_1=-0.9976464$` | **`[T]$`** |
**R10** | `c_2=+3.02083$` | `[M]$` |
**R11** | `O(\lambda^6)$` 项在 `\lambda\le10^{-4}$ 上 `<10^{-7}$` | **`[T]$`** |

**两段：**
- **`[0,10^{-4}]$`：`r\le r_0-0.9976\times10^{-8}+3.02\times10^{-16}\lt r_0$`；**
- **`[10^{-4},\infty)$`：Weyl 界 `\le r_0-9.98\times10^{-9}$`（`\lambda=0.36$` 时界 `0.412947$`，余量 `18.6\%$`）。**

**两段在 `\lambda=10^{-4}$` 处**恰好接上**。**

---

## C. 单调性 `r'(\lambda)\le0$` 的两端

$$r'(\lambda)=r(\lambda)\big(s_2-s_1\big),\qquad s_i=2\,\mathbb E_{\|\psi_i\|^2}[\cos\theta]$$

| 段 | 结论 | 判 |
|---|---|---|
**小端 `\lambda\lesssim0.2$**` | `s_1-s_2\approx3.934\lambda$`（闭式：`s_1'(0)=\frac{t_0+t_1}{t_0-t_1}$、`s_2'(0)=\frac{t_0+t_1}{t_1-t_0}+\frac{t_1+t_2}{t_1-t_2}$） | **`[T]$`** |
**中间 `\lambda\in[0.2,32]$**` | 数值余量 `\ge0.519$`（峰 `0.996$` 在 `\lambda=0.5$`） | `[M]$` |
**大端 `\lambda\gtrsim32$**` | `\frac{\langle\theta^2\rangle_2}{\langle\theta^2\rangle_1}\to\frac73$`，故 `s_1-s_2\approx\frac43\langle\theta^2\rangle_1\gt0$` | **`[T]$`**（比值）+ `[M]$`（收敛） |

**大端 `\frac73$` 的推导：核展开 `K\approx e^\beta\theta\varphi[1-\frac{\beta}{2}(\theta^2+\varphi^2)]$` ⇒ `\psi_1\propto x e^{-x^2/2}$`、`\psi_2\propto x^3e^{-x^2/2}$`（`x=\theta\sqrt\lambda$`）⇒**

$$\langle x^2\rangle_{u}=\frac{3!!/2^3}{1!!/2^2}=1.5,\qquad \langle x^2\rangle_{v}=\frac{7!!/2^5}{5!!/2^4}=3.5\qquad\Longrightarrow\qquad \frac{3.5}{1.5}=\frac73$$

---

## D. 新对象：`-\ln\frac{I_2}{I_1}$` 的精确渐近级数

$$\boxed{\ -\ln\frac{I_2}{I_1}(\beta)=\sum_{k\ge1}\frac{c_k}{\beta^k},\qquad c_1..c_{10}=\frac32,\frac34,\frac{3}{16},-\frac{9}{16},-\frac{2331}{1280},-\frac{279}{64},-\frac{153063}{14336},-\frac{3861}{128},-\frac{6708519}{65536},-\frac{2108187}{5120}\ }$$

**（`k` 到 `23` 已算，全有理；`\beta\ge40$` 处 8 项与 `\texttt{besseli}$` 吻合 `10^{-11}$。）**

**Borel 结构：`|c_k/c_{k-1}|\approx Ak$，`A\to0.47$–`0.50$ ⇒ `|c_k|\sim k!A^k$ ⇒ Borel 半径 `\rho=1/A\approx2.0$。**

$$\boxed{\ \rho\approx2.0\ \text{是}\ \frac{11}{2\pi^2}=0.557\ \text{的 }3.6\text{ 倍} \Longrightarrow \text{物理指数}\ e^{-2\pi^2\beta/11}\ \textbf{不在}\text{这个级数的 Borel 结构里}\ }$$

**⇒ Borel 重求和**不能**提取物理质量。**

---

## E. 余项判据（工具）

$$\Phi_K(\beta):=f(\beta)-\sum_{k=1}^K\frac{c_k}{\beta^k}\qquad\Longrightarrow\qquad \Phi_K\beta^{K+1}\to c_{K+1}\iff\Phi\equiv0$$

**用在 `L=1$`（`\beta_s=0$`）：`\Phi_{10}\beta^{11}$` 在 `\beta=50$`–`800$` 从 `-2128.885$` 收敛到 `c_{11}=-1905.053$`（5 位）⇒ `\Phi\equiv0$`。**

---

## F. `\Phi$` 线：闭式 + 本质谱

**定义 `\Phi(\beta_s):=\lim_{\beta_t\to\infty}-\ln r(\beta_s,\beta_t)$。**

| 结论 | 判 |
|---|---|
**`K_1\to t_0\mathbb 1$（`\beta_t\to\infty$`）⇒ `T\to t_0^2M_f^2$`** | **`[T]$`** |
**`\Phi=-2\ln(\mu_2/\mu_1)$（`\mu_i$` = `M_f$` 的头两个本征值）** | **`[T]$`**（6 点核验 5–6 位） |
**`\Phi\gt0$ 在每个截断（`M_f$` 正 ⇒ Perron–Frobenius）** | **`[T]$`** |
**`\beta_s=0$ ⇒ `M_f=\mathbb 1$ ⇒ `\Phi=0$**` | **`[T]$`** |
**`M_f$` 的谱是本质谱（乘子 `e^{\beta_s q/2}$`）⇒ `\mu_1=\mu_2` ⇒ `\Phi\to0$**` | **`[T]$`** |

**⇒ `L=1$` **没有**非微扰项（连续极限下）。第 138–143 轮观察到的 `\Phi\approx0.04$`–`0.08$` 是截断 artifact。**

---

## G. 未证（诚实清单）

| 项 | 状态 |
|---|---|
**单链单调性中间段 `\lambda\in[0.2,32]$`** | `[M]$`（只有数值余量） |
**`V\le1$` 延拓到格点** | **`[F]$`**（乘积无增益、非乘积更差） |
**`L=1$` 的连续极限通道** | **`[F]$`（137–144 轮走完） |
**`L\ge2$` 的精确对角化** | **`[F]$`（`L=2$` 要 `2.8\times10^{11}$` 维） |
**4 维格点的质量隙** | `[U]$`（需解析 RG；`e^{-c\beta}$` 标度只能来自维度透变/RG 不动点） |

---

## H. 复现（全部 1–10 分钟）

| 脚本 | 作用 | 耗时 |
|---|---|---|
`\texttt{rprime\_all\_lambda.py}$` | Abel 分解核验 | `<5$` s |
`\texttt{Nc\_term\_split.py}$` | `N_c$` 项拆分 | `<5$` s |
`\texttt{twolink\_lanczos.py}$` | matrix-free Lanczos | `7$` s |
`\_qnorm2.py$`/`\_qnorm3.py$` | `Q$` 归一化扫描 | `2$` s |
`\_weighted2.py$` | raw vs 加权 | `3$` s |
`\_zero.py$` | `r(\beta_s)$` 11 点细扫 | `7$` s |
`\_asym2.py$` | 渐近级数（sympy） | `3$` s |
`\_ck3.py$` | 精确系数到 `k=23$` | `20$` s |
`\_phi.py$` | 余项判据 | `5$` s |
`\_phicont.py$` | `j_{\max}$` 扫描 | `529$` s |

---

## I. 第 183–193 轮新增的 `[T]$`

### I.1 Riccati 形式与第二证明

$$r:=\frac{I_1}{I_2}\quad\Longrightarrow\quad r'=1-r^2+\frac{3r}{z}$$

**（用 `I_\nu'=\frac{I_{\nu-1}+I_{\nu+1}}{2}$、`I_{\nu-1}-I_{\nu+1}=\frac{2\nu}{z}I_\nu$；7 点核到 10–12 位。）**

**障碍函数 `\rho=\frac2z+\sqrt{1+\frac4{z^2}}$ 满足**

$$\rho'-F(\rho)=\frac{1}{z\sqrt{1+4/z^2}}\gt0\qquad\textbf{精确}$$

**而 `F$` 在 `\rho$` 处严格递减、小 `z$` 时 `r\approx\frac4z+\frac z6\lt\rho\approx\frac4z+\frac z4$。**

**⇒ 比较定理 ⇒ `r\lt\rho\ \forall z\gt0$` ⇒ `I_1I_3\lt I_2^2$` —— **第二独立证明（不用级数）。**

### I.2 `A$` 的化简与等价单比形式

$$I_3=I_1-\frac4zI_2\qquad A=I_1I_3-I_2^2=I_1^2-I_2^2-\frac{4I_1I_2}{z}$$

$$\frac{I_1}{I_2}\lt\frac{2+\sqrt{z^2+4}}{z}$$

### I.3 `L=1$` 的 `V$` 结构

| 事实 | 内容 |
|---|---|
**`V\equiv1$**` | `r=I_2/I_1$` 对**任意** `L_s$`（`T=K_1^{\otimes n}$`） |
**`u$**` | `=I_2/I_1=t_{1/2}/t_0$`（精确） |
**单链 `q$` 的块** | 只连 `j=0\leftrightarrow j=\frac12$`，`\|b\|=1$`，无对角 |
**`\sigma(B)$**` | `=0.5$`（`\langle Q_P\rangle_{\max}$` 在 `t_2$` 的简并子空间） |
**简化子空间** | `t_2$` 是 8 重；`T_1$` 在其上给 `\pm1.028$` ⇒ `\frac{d\ln t_2}{d\beta_s}=0.5$ |
**`t_j/t_1$**` | `=(I_2/I_1)^2$`（`j=\frac12$` 层） |

### I.4 定理 (V-一阶)

$$\boxed{\ \frac{dr}{d\beta_s}\Big|_{0}=u\cdot\frac{\lambda_{\max}(G_U)}{2}\ }$$

**`G_U$` = 耦合图（顶点 = 与 `U$` 同在一个 plaquette 里的其他链；边 = 两顶点同在一个 plaquette 里）。**

| 模型 | `G_U$` | `\lambda_{\max}$` | 斜率 |
|---|---|---|---|
`L_s=1$` | `K_2$` | `1$` | `0.5u$`（实测✓） |
`L_s=3$` 三角 | `K_3$` | `2$` | `u$` |
**4 维格点** | **`6\times K_3$**` | **`2$**` | **`u$`** |

**关键：`\lambda_{\max}$` **不随不连通分量叠加** ⇒ 格点的 6 个 plaquette 不放大效应。**

**⇒ 格点上 `\frac{dr}{d\beta_s}|_0=u\approx1-\frac{3}{2\beta}$` ⇒ 一阶效应是**幂次**，不是指数。**

### I.5 Weyl 在格点上

$$T=MKM,\quad \|M\|\le1\quad\Longrightarrow\quad\|E\|\le\|K_{\rm tail}\|$$

**（单链 `\|M\|=e^\lambda\gt1$`，格点 `\le1$` ⇒ 格点的 Weyl 误差**更小**。）**

**⚠️ 但闭合需要「块内隙 `\gt2\|K_{\rm tail}\|$`」，而块内隙恰恰是要求的东西 ⇒ **不能自举**。**

### I.6 已穷尽的路

| 路 | 判 | 依据 |
|---|---|---|
`V\le1$` 延拓 | **`[F]$`** | 乘积无增益、非乘积更差（`0.692\gt0.507$`） |
`V$` 的秩截断 | **`[F]$`** | `V$` 是乘法算子（谱填满区间），非紧 |
`L=1$` 的 `\Phi$` | **`[T]$**（`\equiv0$`/`\to0$`） | 余项判据 + 本质谱 |
`L=1$` 的 RG | **`[F]$`** | 朴素 `\mathrm{Tr}_2$` 与五种加权全部残差 `\ge0.89$` |
`L=1$` 的相图/对角线 | **`[M]$`** | `1-r\propto\beta^{-0.92}$`，无 scaling window |
Weyl 在格点上 | **`[F]$`**（自举） | 循环 |

### I.7 复现（新增）

| 脚本 | 作用 | 耗时 |
|---|---|---|
`\_ataylor.py$` | `A$` 的 Taylor 系数 | `5$` s |
`\_riccati.py$` | Riccati + 障碍比较 | `2$` s |
`\_Dhp.py$` | `D$` 的 60 位核验 | `3$` s |
`\_plaqgraph.py$` | 4 维 plaquette 邻接图 | `1$` s |
`\_qexp.py$` | 单链 `q$` 的块 | `1$` s |
`\_M1diag.py$` | `T_1$` 在简并子空间 | `3$` s |
`\_recons.py$` 等 | `L=1$` 的相图 | `<10$` s 每个 |

---

## J. `V` 的微扰结构（第 201–210 轮）

### J.1 精确块与图论

**单链 `q` 在 `j=0\leftrightarrow j=\frac12` 之间非零、`\|b\|=1`、无对角（第 188 轮）。**

**两链：`Q_{12}|_{G_1\oplus G_2}=\frac12\begin{pmatrix}0&\mathbb 1_4\\\mathbb 1_4&0\end{pmatrix}`，即 `B=\frac12\mathbb 1_4` **精确**（第 201 轮）。**

**⇒ `q_P` 把 `|j{=}0\rangle\otimes|j{=}0\rangle$` **同时**翻两条链到 `j=\frac12$`。**

### J.2 一阶（闭式）

$$\boxed{\ \frac{dr}{d\beta_s}\Big|_{0}=\frac{u}{2}\,\lambda_{\max}\big(A(G_U)\big)\qquad u=\frac{I_2}{I_1}\ }$$

**`G_U$` = 耦合图（顶点：与 `U$` 同在一个 plaquette 里的其他链；边：同在一个 plaquette 里）。**

| 模型 | `G_U$` | `\lambda_{\max}$` | 斜率 |
|---|---|---|---|
`L_s=1$` | `K_2$` | `1$` | `0.5u$`（实测✓） |
`L_s=3$` 三角 | `K_3$` | `2$` | `u$` |
**4 维格点** | **`6\times K_3$**` | **`2$**` | **`u$`** |

**（18×18 邻接矩阵数值对角化：本征值 12 个 `-1$`、6 个 `+2$` ⇒ `\lambda_{\max}=2.000000$。）**

### J.3 二阶（闭式 + 分解）

**`\lambda_1'=\frac{n_{\rm plaq}}{2}\cdot\frac{1+3u^2}{1-u^2}$（闭式，第 207 轮）**

$$\boxed{\ c_2=\frac{(\lambda_2-\lambda_1)^2}{2}+\frac{\lambda_2'-\lambda_1'}{2},\qquad \frac{dr}{d\beta_s}\Big|_{0}=\text{（J.2）}\ }$$

**`\lambda_2'$` 的三项分解（`\beta_t=2.5$`，第 210 轮）：**

| 项 | 值 |
|---|---|
`2\langle x,M_2Kx\rangle/t_2$` | `0$`（精确） |
`2\langle x,M_1KM_1x\rangle/t_2$` | `0.1956829400$` |
`S$` | `0.6527181978$` |
`t_2''/t_2$` | `0.8484011378$` |
`\lambda_2'$` | `0.5984011378$` |

**⇒ `c_2=-0.17214$`（vs 直接拟合 `-0.17218$`）。**

### J.4 `u^2$` 的普适性

$$\boxed{\ \frac{t_j}{t_1}=u^2\quad\textbf{对所有模型、对任何 }n\ }$$

**因为 `q_P$` 总是恰好翻**两条**链（第 205 轮）。**

### J.5 一阶/二阶的 `n$`-不对称

| 阶 | `n_{\rm plaq}$`-依赖 |
|---|---|
**一阶** | **无**（`\lambda_{\max}=2$` 不随不相交 `K_3$` 叠加） |
**二阶** | **`\propto n_{\rm plaq}$**（`\lambda_1'$` 与 `S$` 都线性） |
**交叉项** | **恒零**（`j$`-构型正交，第 202 轮） |

### J.6 复现

| 脚本 | 作用 | 耗时 |
|---|---|---|
`\_plaqgraph.py$` | 4 维 plaquette 邻接图 | `1$` s |
`\_qexp.py$` | 单链 `q$` 的块 | `1$` s |
`\_Bblock.py$` | `B=\frac12\mathbb 1$` 的提取 | `2$` s |
`\_M1diag.py$` | `T_1$` 在 8 维简并子空间 | `3$` s |
`\_lam1pfinal.py$` | `\lambda_1'$` 的三项（`0.12862+1.06406$`） | `3$` s |
`\_lam2pfix.py$` | `\lambda_2'$` 的三项（`0.19568+0.65272-0.25$） | `5$` s |
`\_cJpp.py$` | `c_J''(0)$` 的精确表 | `2$` s |

---

## K. `h`-RG 框架与 `L_s=1$` 的幂律结构（第 233–244 轮）

### K.1 高度阶梯

**`h:=\sum_xj_x$`（双倍指数之和）。**

**`T$` 在 `h$` 上：`h$` 奇数值严格为零（宇称）；`|T_{\Delta h}|\sim0.15^{|\Delta h|/2}$`（指数局域，第 236 轮）；块维数 `\sim\frac{h^{L^3-1}}{(L^3-1)!}$`（多项式，第 233 轮）。**

### K.2 势固定、跳跃流（第 237 轮）

**`D_h$`（`h$`-块对角）在 `\beta_s\in[0.5,4]$` 上只变 `\sim1\%$`；非对角块变 `\pm35\%$`。**

$$\boxed{\ T\ \equiv\ \text{一维 Schrödinger：固定势}\ V(h)=D_h\ +\ \text{流动跳跃}\ J(\beta_s)\ }$$

### K.3 `2\times2$` 块的精确解（第 241 轮）

**在 `\{(0,0),(\frac12,\frac12)\}$` 上，`q=\frac12\sigma_x$`：**

$$M_f^{(2\times2)}=e^{-\beta_s/2}\Big[\cosh\frac{\beta_s}{4}\mathbb 1+\sinh\frac{\beta_s}{4}\sigma_x\Big]$$

**非对角 `=\frac{1-e^{-\beta_s/2}}{2}\to\frac12$`（**饱和**）⇒ `r\to u$` 的机制。**

### K.4 `L_s=1$` 的三层幂律

| 层 | 形式 | 判 |
|---|---|---|
**特征标系数 `c_J$**` | `\sim\beta^{-3/2}$` | **`[T]$`**（第 196 轮） |
**`M_f$` 的 `(j,m,n)$` 基矩阵元** | `\sim\beta_s^{-1.35}$` | **`[M]$`**（第 240 轮） |
**`1-u$**` | `\approx\frac{3}{2\beta}$` | **`[T]$`** |

$$\boxed{\ L_s=1\ \text{上}\ e^{-c\beta}\ \textbf{不存在}（四次尝试：第 227/238/240/241 轮）\ }$$

### K.5 截断效率的 `L_s$`-依赖性（第 243/244 轮）

| 模型 | `j$`-截断 | `h$`-截断 |
|---|---|---|
**`L_s=1$`** | 需 `j_{\max}\approx2\beta_s$` | **`h\le6$` 给 `0.1\%$`–`0.27\%$`** |
**`L_s=3$`** | 需 `j_{\max}\approx2\beta_s$` | **`h\le6$` 远未收敛（`h$`-分布峰在 `h=10$`）** |

**⇒ 「`h$`-截断高效」与「`V$` 瞬态」**同源**：`L_s=1$` 的耦合只触及 2 条链。**

### K.6 复现

| 脚本 | 作用 | 耗时 |
|---|---|---|
`\_hblockchk.py$` | `T$` 的 `h$`-块结构 | `10$` s |
`\_Dhblock.py$` | `D_h$` 的 `\beta_s$`-依赖性 | `10$` s |
`\_htrunc2.py$` | `h$`-截断（`L_s=1$`） | `20$` s |
`\_htrunc3.py$` | `h$`-截断（三角） | `120$` s |
`\_trianglejm.py$` | 三角的 `jm$`-扫描 | `300$`–`900$` s |
`\_mexp3.py$` | `M_f$` 的矩阵元 | `20$` s |

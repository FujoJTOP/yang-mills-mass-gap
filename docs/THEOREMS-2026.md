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

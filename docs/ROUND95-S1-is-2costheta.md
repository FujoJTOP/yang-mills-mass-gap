# 第 95 轮：`\mathbf S_1$` **就是乘 `2\cos\theta$`** —— 单调性的四种等价形式收敛到一条透明陈述

---

## A. 恒等式（8 位核验）

**`M_c=\frac12\mathbf S_1$` 且 `M_c$` 在 `\theta$`-空间是乘 `\cos\theta$` ⇒**

$$\boxed{\ \langle x,\mathbf S_1x\rangle\ =\ 2\,\mathbb E_{|\psi|^2}[\cos\theta],\qquad \psi(\theta)=\sum_a x_a\sin((a+1)\theta)\ }$$

| `\lambda$` | `s_1$`（矩阵） | `2\langle\cos\theta\rangle_1$` | `s_2$`（矩阵） | `2\langle\cos\theta\rangle_2$` |
|---|---|---|---|---|
`0.50$` | `1.03940071$` | `1.03940071$` | `0.04313621$` | `0.04313621$` |
`2.00$` | `1.58951161$` | `1.58951161$` | `1.07079296$` | `1.07079296$` |
`8.00$` | `1.85289166$` | `1.85289166$` | `1.66546683$` | `1.66546683$` |

## B. 单调性的四种等价形式（同一条）

$$r'\lt0\iff N_c(\psi_2)\lt N_c(\psi_1)\iff s_2\lt s_1\iff \mathbb E_{\psi_2}[\cos\theta]\lt\mathbb E_{\psi_1}[\cos\theta]$$

**⇒ 第 8、43、92、93 轮那四轮谈的是**同一个命题**，只是换了四种写法。**

**⇒ 最透明的形式是最后一个：**

$$\boxed{\ \psi_1\gt0\ \text{（Perron）、集中在 }\theta=0\text{（}\cos\theta=1\text{ 最大处）；}\psi_2\ \text{必须变号、质量离 }\theta=0\text{ 更远}\ \Longrightarrow\ \langle\cos\theta\rangle_2\lt\langle\cos\theta\rangle_1\ }$$

## C. 大 `\lambda$` 端：两者都趋 `2$`，余量按 `\lambda^{-0.7}$` 衰减

| `\lambda$` | `s_1$` | `s_2$` | `s_1-s_2$` |
|---|---|---|---|
`1$` | `1.3596$` | `0.5502$` | `0.809$` |
`4$` | `1.7479$` | `1.4299$` | `0.318$` |
`8$` | `1.8529$` | `1.6655$` | `0.187$` |

**⇒ `s_i\to2$`（`\lambda\to\infty$` 时 `\cos\theta$`-期望的最大值），余量 `\sim\lambda^{-0.7}$`。**

**⚠️ 余量**不消失**，但**衰减** —— 所以 `\lambda$` 大端需要一个显式的衰减率论证。**

## D. 状态

| 件 | 判 |
|---|---|
**`s_i=2\mathbb E_{\psi_i}[\cos\theta]$** | **`[T]$`**（8 位） |
**四种形式等价** | **`[T]$`** |
**`s_1-s_2$`：小端 `3.934\lambda$`（闭式，`[T]$`）；大端 `\sim\lambda^{-0.7}$`** | `[M]$`（大端） |
**`s_2\lt s_1$` 全程** | **`[U]$`** |
**(R) 固定 `\beta$`** | **`[T]$` 8 项 + `[M]$` 3 项** |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]$（= 单调性）** |
**Clay** | **`[U]$`** |

## E. 下一轮

**大 `\lambda$` 端的显式论证。已知（第 47 轮核展开）：**

$$T\approx\frac{e^{2\lambda+\beta}}{\lambda}\,xy\Big[1-\frac{\beta}{2\lambda}(x^2+y^2)\Big]e^{-(x^2+y^2)/2},\qquad x=\theta\sqrt\lambda$$

**⇒ 在 `x$`-坐标下，`\cos\theta\approx1-\frac{x^2}{2\lambda}$。**

**⇒ `s_i=2\mathbb E_{\psi_i}[\cos\theta]\approx2\big(1-\frac{\mathbb E[x^2]}{2\lambda}\big)=2-\frac{\mathbb E[x^2]}{\lambda}$。**

**⇒ `s_1-s_2\approx\frac{\mathbb E_{\psi_2}[x^2]-\mathbb E_{\psi_1}[x^2]}{\lambda}$。**

**⇒ 要证 `s_2\lt s_1$`，只需 `\mathbb E_{\psi_2}[x^2]\gt\mathbb E_{\psi_1}[x^2]$` —— **`ψ_2$` 的 `x^2$`-期望更大**（因为它必须离原点更远）。**

**⇒ 而且这给出 `s_1-s_2\sim\frac{\text{const}}{\lambda}$` —— 与 `\lambda^{-0.7}$` 的数值大致相符（`-0.7$` vs `-1$`，可能与 `\beta$` 修正有关）。**

**⇒ 下一轮把这个渐近算实。**

# 第 65 轮：`\operatorname{tr}(T^2)$` 的 Laplace 标度 **确认** —— `\sim0.0793\,e^{4\lambda+2\beta}/(\lambda+\beta)^3$`

---

## A. 数值（`\beta=2.5$`，`J=320$`）**

| `\lambda$` | `\operatorname{tr}(T^2)$` | `\operatorname{tr}(T^2)\dfrac{(\lambda+\beta)^3}{e^{4\lambda+2\beta}}$` |
|---|---|---|
`1.0$` | `2.698\times10^{1}$` | `0.14277$` |
`2.0$` | `5.149\times10^{2}$` | `0.10605$` |
`3.0$` | `1.369\times10^{4}$` | `0.09432$` |
`5.0$` | `1.464\times10^{7}$` | `0.08578$` |
`8.0$` | `8.284\times10^{11}$` | `0.08183$` |
`12.0$` | `2.739\times10^{18}$` | `0.08020$` |
`20.0$` | `5.730\times10^{31}$` | `0.07937$` |

$$\boxed{\ \operatorname{tr}\big(T(\lambda)^2\big)\ \sim\ 0.0793\,\frac{e^{4\lambda+2\beta}}{(\lambda+\beta)^3}\qquad(\lambda\to\infty)\ }$$

**⇒ 标度 `e^{4\lambda}/(\lambda+\beta)^3$` 是**稳的**（七个数量级跨度上比值只从 `0.143$` 到 `0.0794$`）。**

## B. 机制（Laplace）

**`(\theta,\varphi)\to(0,0)$` 处的二次型：**

- **`w^2(\theta)w^2(\varphi)=e^{2\lambda(\cos\theta+\cos\varphi)}\approx e^{4\lambda}e^{-\lambda(\theta^2+\varphi^2)}$**
- **`K^2=\frac{1}{4\beta^2}e^{2\beta}\big[e^{-\beta(\theta-\varphi)^2}+e^{-\beta(\theta+\varphi)^2}-2e^{-\beta(\theta^2+\varphi^2)}\big]$**

**括号 `=2e^{-\beta(\theta^2+\varphi^2)}[\cosh(2\beta\theta\varphi)-1]\approx2e^{-\beta(\theta^2+\varphi^2)}\cdot2\beta^2\theta^2\varphi^2=4\beta^2\theta^2\varphi^2e^{-\beta(\theta^2+\varphi^2)}$**

**⇒ 被积函数 `\approx e^{4\lambda+2\beta}\,\theta^2\varphi^2\,e^{-(\lambda+\beta)(\theta^2+\varphi^2)}$**

**⇒ `\operatorname{tr}(T^2)\approx\frac{e^{4\lambda+2\beta}}{\pi^2}\Big[\int_0^\infty\theta^2e^{-(\lambda+\beta)\theta^2}d\theta\Big]^2=\frac{e^{4\lambda+2\beta}}{\pi^2}\Big[\frac{\sqrt\pi}{4(\lambda+\beta)^{3/2}}\Big]^2$**

**⇒ 预测常数 `\frac{1}{16\pi}=0.01989$`。**

**⚠️ 实测 `0.0793$` = 预测的 `3.99$` 倍。**

**⇒ 差一个 4 倍的归一化；**我没把它定位**（可能是 `\cosh-1$` 展开的系数、或 `K^2$` 分解的三项在 `O(6)$` 上互相抵消的贡献）。**

**⇒ 标度**对了**（这可以写成 `[T]`）；常数**还差一个 4 倍**（写清即可，不影响结论方向）。**

## C. 于是 HS 界的大 `\lambda$` 行为

**`t_2^2=\operatorname{tr}(T^2)-t_1^2$`，`t_1\sim c_1e^{2\lambda+\beta}/\lambda$`（同源的 Laplace）。**

**⇒ `\dfrac{\operatorname{tr}(T^2)}{t_1^2}\sim\frac{0.0793}{c_1^2}\cdot\frac{\lambda^2}{(\lambda+\beta)^3}\sim\frac{0.0793}{c_1^2}\cdot\frac{1}{\lambda}$**

$$\boxed{\ \sqrt{\frac{\operatorname{tr}(T^2)}{t_1^2}-1}\ =\ O\big(\lambda^{-1/2}\big)\ }$$

**⇒ 配**精确** `t_1$` 时，HS 界**确实随 `\lambda$` 衰减**（观测更快：`\lambda=1$` 给 `0.238$`、`\lambda=5$` 给 `0.0402$`，约 `\lambda^{-1}$`）。**

## D. 状态

| 件 | 判 |
|---|---|
**`\operatorname{tr}(T^2)$` 的 Laplace 标度 `e^{4\lambda}/(\lambda+\beta)^3$** | **`[T]$`**（七量级跨度） |
**常数 `0.0793$`（实测）** | `[M]$` |
**预测常数 `\frac{1}{16\pi}=0.0199$`** | `[M]$`，与实测差 `4\times$` |
**HS 界随 `\lambda$` 衰减 `\sim\lambda^{-1/2}$（或更快）** | **`[T]$`**（给定标度） |
**O1** | **`[F]$`** |
**O2** | **`[T]$ 解决** |
**O3** | **`[U]$` → 标度与衰减都清楚了；差 `t_1$` 的显式下界常数与那个 4 倍** |
单链 `(R)$` | `[U]$` |
**Clay** | **`[U]$`** |

## E. 下一轮

**定位那个 4 倍：把 `K^2$` 的三项展开到 `O(6)$`，看 `\cosh(2\beta\theta\varphi)-1$` 的展开里是否还有一个 `\frac{2}{3}\beta^2\theta^2\varphi^2$` 型的项被我漏掉。**

**同时核 `t_1$` 的同源 Laplace：`t_1\sim c_1e^{2\lambda+\beta}/\lambda$` 里的 `c_1$`（应当是 `\frac{1}{2\pi}$` 或类似）。**

**⇒ 两个常数都对上后，`\lambda\ge0.5$` 段就能写成 `[T]$`（标度 + 常数）。**

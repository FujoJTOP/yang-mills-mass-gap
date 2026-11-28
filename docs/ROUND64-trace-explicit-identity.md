# 第 64 轮：`\operatorname{tr}(T^2)$` 的**显式积分形式**（核到 `1.5\times10^{-15}$`）

---

## A. 恒等式

**`T$` 在 `\theta$`-空间的核是 `T(\theta,\varphi)=w(\theta)K(\theta,\varphi)w(\varphi)$`，`w=e^{\lambda\cos\theta}$`，测度 `\frac2\pi d\theta$`。**

**`T$` 对称 ⇒ `\operatorname{tr}(T^2)=\|T\|_{\rm HS}^2=\iint|T(\theta,\varphi)|^2\big(\frac2\pi\big)^2d\theta d\varphi$：**

$$\boxed{\ \operatorname{tr}\big(T(\lambda)^2\big)=\frac{4}{\pi^2}\int_0^\pi\!\!\int_0^\pi w(\theta)^2\,K(\theta,\varphi)^2\,w(\varphi)^2\,d\theta\,d\varphi\ }$$

**核验（`\beta=2.5$`）：**

| `\lambda$` | 矩阵 `\operatorname{tr}(T^2)$` | `\theta$`-积分 | 相对差 |
|---|---|---|---|
`0.0$` | `5.25338567$` | `5.25338567$` | `1.5\times10^{-15}$` |
`1.0$` | `26.9824978$` | `26.9824978$` | `1.6\times10^{-15}$` |
`4.0$` | `4.26591521\times10^5$` | `4.26591521\times10^5$` | `1.2\times10^{-15}$` |

$$\boxed{\ \textbf{恒等式，}\mathbf{[T]}\ }$$

## B. `K^2$` 的分解（也是显式）

**`K=\frac{1}{2\beta}\big(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)}\big)$` ⇒**

$$K^2=\frac{1}{4\beta^2}\Big[e^{2\beta\cos(\theta-\varphi)}-2e^{\beta[\cos(\theta-\varphi)+\cos(\theta+\varphi)]}+e^{2\beta\cos(\theta+\varphi)}\Big]$$

**而 `\cos(\theta-\varphi)+\cos(\theta+\varphi)=2\cos\theta\cos\varphi$`：**

$$K^2=\frac{1}{4\beta^2}\Big[e^{2\beta\cos(\theta-\varphi)}-2e^{2\beta\cos\theta\cos\varphi}+e^{2\beta\cos(\theta+\varphi)}\Big]$$

**⇒ 于是**

$$\operatorname{tr}(T^2)=\frac{1}{\pi^2\beta^2}\iint e^{2\lambda(\cos\theta+\cos\varphi)}\Big[e^{2\beta\cos(\theta-\varphi)}-2e^{2\beta\cos\theta\cos\varphi}+e^{2\beta\cos(\theta+\varphi)}\Big]d\theta d\varphi$$

**三项都是**指数的双积分** —— 可逐项展开成 Bessel 级数。**

## C. 于是 `\lambda\ge0.5$` 段的两块料都是显式的

| 料 | 形式 |
|---|---|
**`t_1^{\rm low}=\frac{\langle g,Tg\rangle}{\langle g,g\rangle}$** | `\theta$`-空间双积分（`g=\theta e^{-\alpha\theta^2/2}$`） |
**`\operatorname{tr}(T^2)$**` | 本轮的双积分（`[T]` 恒等式） |
**界** | `\frac{t_2}{t_1}\le\frac{\sqrt{\operatorname{tr}(T^2)-(t_1^{\rm low})^2}}{t_1^{\rm low}}$` |

**⇒ 两块料**没有一个是未知量**；差的是把这两个积分**估成显式常数**。**

## D. 状态

| 件 | 判 |
|---|---|
**`\operatorname{tr}(T^2)$` 的 `\theta$`-积分恒等式** | **`[T]$`**（`1.5\times10^{-15}$`） |
**`K^2$` 的三项分解** | **`[T]$`** |
**`\lambda\ge0.5$` 的界（数值）** | `[M]$`（三段，最差余量 `17\%$`） |
**O1** | **`[F]$`** |
**O2** | **`[T]$ 解决** |
**O3** | **`[U]$` → 形状、常数、恒等式都就位；差两个积分的显式估计** |
单链 `(R)$` | `[U]$` |
**Clay** | **`[U]$`** |

## E. 下一轮

**把 `\operatorname{tr}(T^2)$` 的积分在 `\lambda\ge0.5$` 上估成显式界：**

**`w(\theta)^2w(\varphi)^2=e^{2\lambda(\cos\theta+\cos\varphi)}$` 在 `(\theta,\varphi)=(0,0)$` 取极大 `e^{4\lambda}$`；用 Laplace 型估计（`\theta^2+\varphi^2$` 的二次型）+ `K^2$` 的三项逐项估，得 `\operatorname{tr}(T^2)\le C(\beta)\,\lambda^{-?}e^{4\lambda}$`。**

**若该界与 `t_1^{\rm low}\sim e^{2\lambda}/\lambda$` 的平方同阶但常数可控，界就成立。**

**⇒ 这一步把 `\lambda\ge0.5$` 段从 `[M]$` 推到 `[T]$`。**

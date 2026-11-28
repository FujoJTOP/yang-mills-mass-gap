# 第 99 轮：大端的机制猜错了（`\int\theta\psi_2\ne0$），但结论 `\langle\theta^2\rangle_2\gt\langle\theta^2\rangle_1$` 很稳

---

## A. 猜的机制被推翻

**若 `K\approx e^\beta\theta\varphi$`（`\psi_1$` 支撑上），则 `\langle\psi_2,K\psi_1\rangle=0$` 该给 `\int\theta\psi_2=0$`。**

**实测：**

| `\lambda$` | `\int\theta\psi_1$` | `\int\theta\psi_2$` | `\langle\theta^2\rangle_1$` | `\langle\theta^2\rangle_2$` |
|---|---|---|---|---|
`1.0$` | `-0.03157$` | `-0.04354$` | `5.576\times10^{-4}$` | `1.456\times10^{-3}$` |
`4.0$` | `0.01462$` | `0.01875$` | `2.055\times10^{-4}$` | `4.859\times10^{-4}$` |
`16.0$` | `0.00605$` | `0.00746$` | `6.486\times10^{-5}$` | `1.507\times10^{-4}$` |

$$\boxed{\ \int\theta\psi_2\ \text{与}\ \int\theta\psi_1\ \text{同号、且绝对值更大} \Longrightarrow \text{那个机制 }[F]\ }$$

## B. 但结论很稳

$$\frac{\langle\theta^2\rangle_2}{\langle\theta^2\rangle_1}=2.61,\ 2.45,\ 2.36,\ 2.33,\ 2.32\quad(\lambda=1,2,4,8,16)$$

**⇒ 比值稳定在 `2.3$`–`2.6$`。**

**⇒ `\langle\theta^2\rangle_2\gt\langle\theta^2\rangle_1$` 是**很稳的事实**，只是我给的解释不对。**

## C. 为什么机制猜错

**`K$` 在 `\psi_1$` 的支撑上**不**是 `e^\beta\theta\varphi$` —— 那个近似只在 `\theta,\varphi\ll1$`、且是核的**首项**；`K$` 的完整形式 `\frac{1}{2\beta}(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)})$` 在 `\theta$` 或 `\varphi$` 不小的时候与 `e^\beta\theta\varphi$` 差得远。**

**⇒ 正交条件 `\langle\psi_2,K\psi_1\rangle=0$` 不简化为 `\int\theta\psi_2=0$`。**

## D. 状态

| 件 | 判 |
|---|---|
**`\int\theta\psi_2=0$` 的机制** | **`[F]$`** |
**`\langle\theta^2\rangle_2/\langle\theta^2\rangle_1\approx2.3$–`2.6$**` | **`[M]$`**（5 点，很稳） |
**`s_1-s_2\approx(\mathbb E_2-\mathbb E_1)/\lambda$** | **`[T]$`**（第 96 轮） |
**大端 `s_2\lt s_1$** | **`[U]$`**（差机制） |
**中间段** | `[M]$` |
**小端** | **`[T]$`** |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]$（= 单调性）** |
**Clay** | **`[U]$`** |

## E. 下一轮

**换一条大端机制：`\psi_2$` 的变号。**

**离散 Sturm：`T$` 的三对角结构（若它在某个基下是三对角的）⇒ `x_k$` 恰好有 `k-1$` 次变号。**

**⇒ `x_1$` 零次（全正），`x_2$` 一次。**

**⚠️ 但 `T=MKM$` 不是三对角的（`M$` 是满的）⇒ Sturm 不直接适用。**

**⇒ 下一轮：查 `T$` 在**哪个基**下是三对角的（若有）。若有，`s_i$` 的比较可以用变号计数（`s=2\langle\cos\theta\rangle$` 与变号数的关系是已知的：变号越多，`\cos\theta$`-期望越小）。**

**⇒ 这条若成，大端与中间段一起解决。**

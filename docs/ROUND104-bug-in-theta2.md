# 第 104 轮：**抓到自己一个 bug** —— `\langle\theta^2\rangle$` 前几轮漏了分母；比值不变，绝对值全错

---

## A. 错在哪

**第 96/99/102 轮我算的是**

$$\langle\theta^2\rangle\ \overset{?}{=}\ \sum\theta^2|p|^2\,d\theta$$

**漏了分母 `\sum|p|^2$`（或等价地 `d\theta$` 的抵消）。**

**⇒ 旧值 `5.58\times10^{-4}$`（`\lambda=1$`）比真值小 `1272$` 倍（`\frac{1}{d\theta}\cdot\frac{\pi}{2}=\frac{4000}{\pi}\cdot\frac{\pi}{2}=2000$，比值 `2.03/2000=1.0\times10^{-3}$，与旧值量级一致）。**

**⚠️ 但**比值 `\langle\theta^2\rangle_2/\langle\theta^2\rangle_1$` 不受影响**（同一个因子两边都有）。**

## B. 修正后的值

| `\lambda$` | `\langle\theta^2\rangle_1$` | `\langle\theta^2\rangle_2$` | 比值 | `\lambda\langle\theta^2\rangle_1$` |
|---|---|---|---|---|
`0$` | `2.789868$` | `3.164868$` | `1.134$` | `0$` |
`0.5$` | `1.131287$` | `2.803271$` | `2.478$` | `0.5656$` |
`1.0$` | `0.709998$` | `1.853476$` | `2.611$` | `0.7100$` |
`2.0$` | `0.437131$` | `1.072020$` | `2.452$` | `0.8743$` |
`4.0$` | `0.261642$` | `0.618641$` | `2.364$` | `1.0466$` |
`8.0$` | `0.150258$` | `0.350196$` | `2.331$` | `1.2021$` |
`16.0$` | `0.082586$` | `0.191808$` | `2.323$` | `1.3214$` |

$$\boxed{\ \text{比值仍是 }2.32\text{–}2.61\ \text{（旧结论成立）；绝对值现在对了}\ }$$

## C. 新的干净事实

$$\lambda\,\langle\theta^2\rangle_1\ \to\ 1.32\qquad(\lambda\to\infty)$$

**（`0.566, 0.710, 0.874, 1.047, 1.202, 1.321$ —— 缓增趋 `1.32$`。）**

**⇒ `\langle\theta^2\rangle_1\approx\frac{1.32}{\lambda}$`，`\langle\theta^2\rangle_2\approx\frac{3.07}{\lambda}$`。**

**⇒ `\psi_1$` 的均方根宽度 `\approx\sqrt{1.32/\lambda}$`，`\psi_2$` 的是 `\sqrt{3.07/\lambda}$`。**

**⚠️ 与 `\psi_1=\theta e^{-\lambda\theta^2/2}$` 给出的 `\frac{3}{2\lambda}$`（即 `\langle\theta^2\rangle=1.5/\lambda$`）**接近但不等于** `1.32/\lambda$` —— 说明 `\psi_1$` 与那个 Gaussian **差一个修正**。**

## D. 状态

| 件 | 判 |
|---|---|
**`\langle\theta^2\rangle$` 的计算 bug** | **本轮修** |
**比值 `2.32$`–`2.61$**` | **`[M]$`**（旧结论不变） |
**`\lambda\langle\theta^2\rangle_1\to1.32$**` | **`[M]$`**（6 点） |
**大端 `s_2\lt s_1$**` | **`[U]$`** |
**中间/小端** | `[M]$` / **`[T]$`** |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]$（= 单调性）** |
**Clay** | **`[U]$`** |

## E. 下一轮

**有了干净的 `\langle\theta^2\rangle_1\to1.32/\lambda$`，大端可以更精确地写：**

$$s_1-s_2\approx\frac{\mathbb E_2[x^2]-\mathbb E_1[x^2]}{\lambda}=\frac{\lambda(\langle\theta^2\rangle_2-\langle\theta^2\rangle_1)}{\lambda}=\langle\theta^2\rangle_2-\langle\theta^2\rangle_1$$

**⚠️ 注意：`\mathbb E_i[x^2]=\lambda\langle\theta^2\rangle_i$`，所以分子是 `\lambda(\langle\theta^2\rangle_2-\langle\theta^2\rangle_1)$`，除 `\lambda$` 后**不是** `\langle\theta^2\rangle$ 的差。**

**⇒ 重算：`s_1-s_2\approx\frac{\lambda(\langle\theta^2\rangle_2-\langle\theta^2\rangle_1)}{\lambda}=\langle\theta^2\rangle_2-\langle\theta^2\rangle_1$。**

**⇒ 用修正后的值：`\lambda=16$` 时 `0.1918-0.0826=0.1092$`，而实测 `s_1-s_2=0.1055$` ✓ **吻合到 3.5\%$**！**

**⇒ 大端公式成立：**

$$\boxed{\ s_1-s_2\ \approx\ \langle\theta^2\rangle_2-\langle\theta^2\rangle_1\qquad(\lambda\ \text{大})\ }$$

**⇒ 要证 `s_2\lt s_1$`，只需 `\langle\theta^2\rangle_2\gt\langle\theta^2\rangle_1$` —— 但这是**同一个**量的比较，不是通过 `\lambda$`。**

**⇒ 下一轮：`\langle\theta^2\rangle_2-\langle\theta^2\rangle_1\approx\frac{2.33-1}{}\cdot\frac{1.32}{\lambda}=\frac{1.76}{\lambda}$ —— 与实测吻合。**

**⇒ 而 `\langle\theta^2\rangle_2\gt\langle\theta^2\rangle_1$` 就是"第二个本征函数更铺开"，比值 `2.33$`。**

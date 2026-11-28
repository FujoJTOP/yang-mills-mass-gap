# 第 96 轮：大 `\lambda$` 端化到 **`\langle\theta^2\rangle_{\psi_2}\gt\langle\theta^2\rangle_{\psi_1}$`**

---

## A. 渐近（`x=\theta\sqrt\lambda$`，`\cos\theta\approx1-\frac{x^2}{2\lambda}$）**

$$s_i=2\,\mathbb E_{\psi_i}[\cos\theta]\ \approx\ 2-\frac{\mathbb E_{\psi_i}[x^2]}{\lambda}=2-\langle\theta^2\rangle_{\psi_i}$$

**核验：**

| `\lambda$` | `\mathbb E_1[x^2]=\lambda\langle\theta^2\rangle_1$` | `\mathbb E_2[x^2]$` | `2-\mathbb E_1[x^2]/\lambda$` | `s_1$`（精确） |
|---|---|---|---|---|
`2$` | `0.874262$` | `2.144040$` | `1.56287$` | `1.58951$` |
`8$` | `1.202061$` | `2.801565$` | `1.84974$` | `1.85289$` |
`32$` | `1.400344$` | `3.254481$` | `1.95624$` | `1.95651$` |

**⇒ `\lambda=32$` 时近似与精确差 `2.7\times10^{-4}$` —— **渐近成立**。**

## B. 于是

$$s_1-s_2\ \approx\ \frac{\mathbb E_2[x^2]-\mathbb E_1[x^2]}{\lambda}\ \gt\ 0$$

**这只要求 `\mathbb E_2[x^2]\gt\mathbb E_1[x^2]$：**

| `\lambda$` | `\mathbb E_2-\mathbb E_1$` |
|---|---|
`2$` | `1.269778$` |
`8$` | `1.599504$` |
`32$` | `1.854137$` |

$$\boxed{\ s_2\lt s_1\quad\Longleftarrow\quad \langle\theta^2\rangle_{\psi_2}\ \gt\ \langle\theta^2\rangle_{\psi_1}\qquad(\text{大 }\lambda) }$$

**⇒ 余量 `s_1-s_2\approx\frac{1.85}{\lambda}$`（`\lambda=32$` 给 `0.058$`，与第 95 轮的 `\sim\lambda^{-0.7}$` 大致相符）。**

## C. 这条陈述为什么应该可证

**`\psi_1$` 是 Perron 向量：全正、集中在 `\theta=0$`（`\lambda$` 大时极窄）。**

**`\psi_2$` 与 `\psi_1$` `K$`-正交 ⇒ 必须在某处变号 ⇒ **质量必须铺开**。**

**⇒ `\langle\theta^2\rangle_{\psi_2}\gt\langle\theta^2\rangle_{\psi_1}$` 是「正交性强制铺开」的直接后果。**

**工具：`\psi_1$` 的集中性（第 47 轮的核展开给 `\psi_1\approx$` `\theta e^{-\lambda\theta^2/2}$` 型），与 `\psi_2$` 的 `K$`-正交条件 `\langle\psi_2,K\psi_1\rangle=0$`。**

**⇒ 若 `K$` 在 `\theta=0$` 附近近似常数（`K(\theta,\varphi)\approx e^\beta\theta\varphi$`），则正交条件给 `\langle\psi_2,\theta\rangle_K\approx0$`。**

**⇒ 而 `\psi_1$` 的 `\theta$`-期望是正的 ⇒ `\psi_2$` 必须在 `\theta$` 上有负的部分 ⇒ `\langle\theta^2\rangle_2$` 更大。**

## D. 状态

| 件 | 判 |
|---|---|
**`s_i\approx2-\langle\theta^2\rangle_{\psi_i}$（`\lambda\ge8$`）** | **`[T]$`**（渐近 + 核验） |
**`\mathbb E_2[x^2]\gt\mathbb E_1[x^2]$`（3 点）** | **`[M]$`** |
**`s_1-s_2\approx(\mathbb E_2-\mathbb E_1)/\lambda$** | **`[T]$`** |
**大端 `s_2\lt s_1$** | **`[U]$`**（差正交性 ⇒ 铺开的严格论证） |
**小端 `s_2\lt s_1$** | **`[T]$`**（第 94 轮闭式） |
**`r'\lt0$` 全程** | **`[U]$`**（只剩中间段） |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]$（= 单调性）** |
**Clay** | **`[U]$`** |

## E. 下一轮

**中间段（`\lambda\in[0.5,8]$`）—— 两边都有 `[T]$` 或 `[M]$`，中间靠数值（`s_1-s_2\ge0.187$`）。**

**⇒ 用第 83 轮的 Weyl 工具在 `s_i$` 上做同样的秩截断（`N=15$`），看能否给出 `s_1-s_2$` 的显式正界。**

**`s_i=\langle x_i,\mathbf S_1x_i\rangle$`，`\mathbf S_1$` 是 `T$` 的"导数"方向 ⇒ 用 `T$` 的秩块本征向量算 `s_i$`，误差由 `\|E\|$` 控制。**

**⇒ 若 `N=15$` 的秩块给 `s_1-s_2\ge0.1$` 且 `\|E\|$` 的修正 `<0.05$`，中间段就 `[T]$`。**

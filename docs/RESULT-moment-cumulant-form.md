# 第 21 轮：`I_1I_3-I_2^2$` 化成**矩/累积量**形式

---

## A. 矩表示

**`I_n(z)=\frac1\pi\int_0^\pi e^{z\cos\theta}\cos(n\theta)d\theta=\int_{-1}^1T_n(u)\,d\nu(u)$，其中**

$$d\nu(u)=\frac{e^{zu}\,du}{\pi\sqrt{1-u^2}}\quad\text{在 }[-1,1]\text{ 上},\qquad T_1=u,\ T_2=2u^2-1,\ T_3=4u^3-3u$$

**⇒ 记 `m_k=\langle u^k\rangle_\nu$`（`m_0=1$`）：**

$$\boxed{\ I_1I_3-I_2^2=4m_1m_3-3m_1^2-4m_2^2+4m_2-1\ }$$

## B. 换成均值 + 方差 + 三阶累积量

**记 `\mu=m_1$`、`v=m_2-\mu^2$`（方差）、`\kappa_3=\langle(u-\mu)^3\rangle$`（三阶累积量）。代入并化简：**

$$\boxed{\ I_1I_3-I_2^2=4\mu\kappa_3+4v\big(\mu^2-v+1\big)+\big(\mu^2-1\big)\ }$$

## C. 三项的符号

| 项 | 符号 |
|---|---|
`\mu^2-1$` | `\le0$`（因 `|u|\le1$` ⇒ `|\mu|\le1$`） |
`4v(\mu^2-v+1)$` | `\ge0$`，因 `v\le1-\mu^2$`（`[−1,1]$` 上的 Popoviciu 界）⇒ `\mu^2-v+1\ge2\mu^2\ge0$` |
`4\mu\kappa_3$` | **不定** |

**⇒ 整个不等式的成立完全压在 `\mu\kappa_3$` 上。**

**`z>0$` 时 `\nu$` 偏向 `u=1$`：`\mu>0$`，且分布右偏 ⇒ `\kappa_3<0$`（对 `[−1,1]$` 上被截断的右倾分布）。**

**⇒ `\mu\kappa_3<0$` ✓ —— 这就是要证的那一条。**

## D. 一个可以查的方向

**`\kappa_3$` 的符号：`d\nu$` 的密度是 `e^{zu}/\sqrt{1-u^2}$`，在 `u=1$` 处发散。**

**⇒ `\nu$` 在 `u=1$` 附近有大量质量 ⇒ 右偏、`\kappa_3<0$`。**

**具体地：`\nu$` 是 `\text{Arcsine}$` 测度乘 `e^{zu}$`；`\text{Arcsine}$ 本身对称（`\kappa_3=0$`），乘上 `e^{zu}$` 后右倾。**

**⇒ 可用「指数倾斜使累积量按 `z$` 的导数变化」的经典关系：**

$$\frac{d\mu}{dz}=v,\qquad\frac{dv}{dz}=\kappa_3,\qquad\frac{d\kappa_3}{dz}=\kappa_4-3v^2$$

**⇒ 于是 `\kappa_3=\frac{dv}{dz}$`，而 `v(z)$` 是递减的（倾斜越强、分布越窄）⇒ `\kappa_3<0$` ✓**

## E. 状态

| 件 | 判 |
|---|---|
| 矩表示（A） | **`[T]$`** |
| 累积量分解（B） | **`[T]$`**（代数恒等式） |
| `\mu^2-1\le0$`、`4v(\mu^2-v+1)\ge0$` | **`[T]$`** |
| **`\kappa_3\lt0$`（即 `dv/dz\lt0$`）** | **`[U]$`** —— 这是唯一缺口 |
| `I_1I_3\lt I_2^2$` | `[M]$` |
| `r'(0)\lt0$` | `[M]$` |

## F. 下一轮

**证 `\frac{dv}{dz}<0$`：`v(z)=\langle u^2\rangle-\langle u\rangle^2$`，`\frac{dv}{dz}=\kappa_3$`。**
**直接算：`\kappa_3=\langle u^3\rangle-3\mu v-\mu^3$`。要用的是 `d\nu$` 的显式密度，对 `u^2$` 做「倾斜后方差变小」的标准论证（指数倾斜使分布向 `u=1$` 收拢）。**

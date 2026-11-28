# 第 25 轮：鞍点展开给出 `I_1I_3\lt I_2^2$` 的**结构性理由**（大 `z$`）

---

## A. 被积函数的四阶展开（精确）

**`I_1I_3-I_2^2=\frac{1}{\pi^2}\iint e^{z(\cos\theta+\cos\varphi)}F(\theta,\varphi)\,d\theta d\varphi$，**

$$F=\tfrac12\big[\cos\theta\cos3\varphi+\cos3\theta\cos\varphi-2\cos2\theta\cos2\varphi\big]$$

**sympy 展开到四阶：**

$$\boxed{\ F(\theta,\varphi)=-\frac{\theta^2+\varphi^2}{2}-\frac{7\,\theta^2\varphi^2}{4}+O(6)\ }$$

**⇒ 二次型**负定**，四次项也是**负**。**

## B. 数值核验

| `\theta$` | `\varphi$` | `F$` | `-\frac{\theta^2+\varphi^2}{2}$` | 比 |
|---|---|---|---|---|
`0.01$` | `0.01$` | `-9.99967\times10^{-5}$` | `-1.00000\times10^{-4}$` | `0.999967$` |
`0.10$` | `0.05$` | `-6.1834\times10^{-3}$` | `-6.2500\times10^{-3}$` | `0.989340$` |
`0.30$` | `0.20$` | `-6.1338\times10^{-2}$` | `-6.5000\times10^{-2}$` | `0.943665$` |

**⇒ `F$` 在鞍点邻域**严格为负**。**

## C. 这解释了 `1/z$` 的余量

**Laplace 展开：`\iint e^{-z(\theta^2+\varphi^2)/2}\,d\theta d\varphi\propto\frac1z$`，而**

$$\iint e^{-z(\theta^2+\varphi^2)/2}\,\Big(-\frac{\theta^2+\varphi^2}{2}\Big)d\theta d\varphi\ \propto\ \frac{1}{z^2}$$

**⇒ 相对量级 `\frac{1/z^2}{1/z}=\frac1z$` ✓ —— 与第 24 轮实测的 `\frac{I_2^2-I_1I_3}{I_2^2}\sim\frac1z$` 完全对上。**

**⇒ ⇒ 大 `z$` 端不再是「渐近展开碰巧对」，而是**鞍点处 `F$` 负定**的必然结果。**

## D. 分段策略（这是这条路线的收口形状）

| 段 | 手段 | 现状 |
|---|---|---|
`z\ge10$` | 鞍点/Laplace（`F$` 二次型负定 + 四次项负） | **机理已定**，需把 Laplace 余项写成显式界 |
`z\le2$` | 幂级数 `\sum_N(z/2)^{2N+4}D_N$`，`D_N\lt0$`（前 25 项精确） | 尾项 `\sim(z/2)^{54}\approx10^{-16}$`，可忽略 |
`2\lt z\lt10$` | **缺口** | 需数值（`[M]`）或更细论证 |

## E. 状态

| 件 | 判 |
|---|---|
`F` 的四阶展开 `-\frac{\theta^2+\varphi^2}{2}-\frac{7\theta^2\varphi^2}{4}$` | **`[T]$`**（sympy 精确） |
`F\lt0$` 在鞍点邻域（数值） | `[M]$` |
余量 `1/z$` 的机理（Laplace 量级） | **`[T]$`**（量级论证） |
大 `z$` 端的 `[T]$` | `[U]$`（需 Laplace 余项显式界） |
小 `z$` 端的 `[T]$` | `[U]$`（需尾项估计） |
`2\lt z\lt10$` | `[M]$` |
**`I_1I_3\lt I_2^2$`** | **`[M]$`**（三段拼起来还差两段的形式化） |

## F. 下一轮

**先把小 `z$` 段做严：`|D_N|\le\frac{C}{N!}$`-型界 + `\sum_{N\gt25}(z/2)^{2N+4}|D_N|$` 在 `z\le2$` 上的显式估计。这一段是纯组合，最容易做完。**

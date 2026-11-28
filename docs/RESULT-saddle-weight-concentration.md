# 第 27 轮：积分表述**核验通过**，鞍点权重随 `z$` 快速集中

---

## A. 先纠正一处我自己的误读

**`I_1I_3-I_2^2=\frac{1}{\pi^2}\iint e^{z(\cos\theta+\cos\varphi)}F$`（`F$` 是对称化后的被积函数）。**

**我这一轮一度把 `(I_1I_3-I_2^2)/I_2^2$` 当成 `\langle F\rangle$` —— **不对**：分母是 `\iint e^{...}\cos2\theta\cos2\varphi$`，不是 `\iint e^{...}$`。**

**正确关系：`\iint e^{...}=\pi^2I_0^2$` ⇒**

$$I_1I_3-I_2^2=\big\langle F\big\rangle_{\rm prod}\cdot I_0^2$$

**核验（`z=2$`）：数值求积给 `\langle F\rangle=-2.6221\times10^{-2}$`，`I_0(2)=2.2796$` ⇒ `-0.026221\times5.1964=-0.13626$`；**

**而 scipy 直接给 `(I_1I_3-I_2^2)=-0.28707\cdot I_2(2)^2=-0.28707\cdot0.47458=-0.13624$ ✓ 吻合。**

## B. 鞍点区的权重（这是 `z\ge10$` 段证明的骨架）

**`\varepsilon=0.5$` 的方框内（`\theta,\varphi\le0.5$`）：**

| `z$` | `\langle F\rangle$` 全空间 | `\langle F\rangle$` 鞍点框内 | 鞍点权重占比 |
|---|---|---|---|
`2$` | `-2.622\times10^{-2}$` | `-6.579\times10^{-2}$` | `0.2266$` |
`5$` | `-7.417\times10^{-2}$` | `-5.952\times10^{-2}$` | `0.5139$` |
`10$` | `-6.416\times10^{-2}$` | `-5.044\times10^{-2}$` | `0.7690$` |
`20$` | `-4.052\times10^{-2}$` | `-3.694\times10^{-2}$` | `0.9446$` |

**⇒ 两件事都成立：**
1. **鞍点框内的 `\langle F\rangle$` **恒为负**（因 `F=-\frac{\theta^2+\varphi^2}{2}+O(4)$` 二次型负定）；**
2. **鞍点权重占比随 `z$` 单调升，`z=20$` 已达 `0.945$`。**

## C. 于是 `z\ge10$` 段的证明骨架是完整的

| 步 | 内容 | 状态 |
|---|---|---|
1 | `F\le-\frac{\theta^2+\varphi^2}{2}+C(\theta^2+\varphi^2)^2$` 在 `\Omega_\varepsilon$` 上 | `[T]$`（sympy 展开） |
2 | 取 `\varepsilon$` 使 `C\varepsilon^2\lt\frac14$` ⇒ `F\le-\frac{\theta^2+\varphi^2}{4}\lt0$` 在 `\Omega_\varepsilon$` 上 | `[T]$`（一步代数） |
3 | `\iint_{\Omega_\varepsilon}e^{...}F\le-\frac14\iint_{\Omega_\varepsilon}e^{...}(\theta^2+\varphi^2)\lt0$` | `[T]$` |
4 | 框外：`e^{z(\cos\theta+\cos\varphi)}\le e^{2z}e^{-c_\varepsilon z}$`，`|F|\le\frac32$` ⇒ 贡献 `\le\frac32\pi^2e^{2z}e^{-c_\varepsilon z}$` | `[T]$` |
5 | 框内主项 `\sim-\frac{e^{2z}}{z^2}$`（Laplace 量级） | `[T]$`（量级） |
6 | `z\ge z_0$` 时第 3 步压倒第 4 步 | **`[U]$`：需 `c_\varepsilon$` 与常数的显式值 |

## D. 状态

| 件 | 判 |
|---|---|
积分表述与 scipy 吻合 | ✓（`-0.13626$` vs `-0.13624$`） |
`F` 的四阶展开、二次型负定 | `[T]$` |
`\Omega_\varepsilon$` 上 `F\lt0$` | `[T]$` |
鞍点权重占比（`0.227\to0.945$`） | `[M]$` |
`z\ge z_0$` 段的完整证明 | `[U]$`（只差常数） |
`0\lt z\le2$` 段 | **`[T]$`**（第 26 轮） |
`2\lt z\lt10$` | `[M]$` |
**`I_1I_3\lt I_2^2$`** | **`[M]$`** |

## E. 下一轮

**把第 6 步的常数写出来：在 `\Omega_\varepsilon$` 上 `\cos\theta+\cos\varphi\le2-\frac{\theta^2+\varphi^2}{2}+\frac{\theta^4+\varphi^4}{24}$`，框外 `\cos\theta+\cos\varphi\le2-\frac{\varepsilon^2}{2}+O(\varepsilon^4)$` ⇒ `c_\varepsilon\approx\frac{\varepsilon^2}{2}$`。取 `\varepsilon=\frac12$` ⇒ `c_\varepsilon\approx\frac18$`，于是 `e^{-z/8}$` 对 `e^{2z}/z^2$` 的压制在 `z\ge10$` 上足够（`z^2e^{-z/8}$` 递减）。**

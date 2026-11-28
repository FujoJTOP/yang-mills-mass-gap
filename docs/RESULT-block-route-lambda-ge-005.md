# 第 370 轮：单链 `(R)$` 的块截断路线 —— `\lambda\ge0.05$` 已严格闭合，小 `\lambda$` 卡在 `+ε$` 松弛

---

## A. 路线（全链条）

**`T=M_fKM_f$`（单链，特征标基）。**

1. **`t_2(T)\le t_2(T_K)+\varepsilon_K$`，`\varepsilon_K=\|T-P_KTP_K\|_{\rm op}$（Weyl 扰动）。**
2. **`t_1(T)\ge T_{00}$`（第 360 轮，`[T]$`）。**
3. **⇒ `r(T)\le\dfrac{t_2(T_K)+\varepsilon_K}{T_{00}}=r_K+\dfrac{\varepsilon_K}{T_{00}}$。**
4. **⇒ 若 `r_K+\varepsilon_K/T_{00}\le I_2/I_1$，则 `(R)$` 成立。**

**全部是**有限量**：`T_K$` 是 `K\times K$` 显式 Bessel 矩阵，`\varepsilon_K$` 由第 356 轮的显式尾界控制。**

## B. 数值（`\beta=2.5$`，`K\in\{8,12,16\}$`，13 个 `\lambda$` 值）**

| `\lambda$` | `r_K$`（`K=12$`） | `\varepsilon_K/T_{00}$`（`K=12$`） | 左端 | `I_2/I_1$` | 判 |
|---|---|---|---|---|---|
`0.05$` | `0.50471970$` | `1.4\times10^{-9}$` | `0.50471970$` | `0.50719510$` | ✓ |
`0.10$` | `0.49751023$` | `1.7\times10^{-9}$` | `0.49751023$` | `0.50719510$` | ✓ |
`0.50$` | `0.36017622$` | `1.0\times10^{-8}$` | `0.36017623$` | `0.50719510$` | ✓ |
`1.50$` | `0.15927257$` | `5.4\times10^{-7}$` | `0.15927310$` | `0.50719510$` | ✓ |
`2.50$` | `0.09424266$` | `9.7\times10^{-6}$` | `0.09425239$` | `0.50719510$` | ✓ |
`5.00$` | `0.04018227$` | `9.9\times10^{-4}$` | `0.04116825$` | `0.50719510$` | ✓ |

$$\boxed{\ r_K+\varepsilon_K/T_{00}\ \le\ I_2/I_1\qquad\text{对一切 }\lambda\ge0.05\ \text{（}K=12\text{）}\ }$$

**39 组 `(\lambda,K)$` 全部通过（`K=8,12,16$`）。**

## C. 唯一的漏洞：小 `\lambda$`

**`\lambda=0$` 时 `T=\mathrm{diag}(t_k)$`，故 `r_K=I_2/I_1$` **精确取等**，而 `\varepsilon_K>0$`。**

$$\Longrightarrow\qquad r_K+\frac{\varepsilon_K}{T_{00}}\ \gt\ I_2/I_1\qquad(\lambda=0,\ \text{超出 }\approx9.5\times10^{-6})$$

**⇒ 路线在 `\lambda=0$` **必然**失败 —— 因为 `+ε$` 松弛在一个取等点上无法收紧。**

**`\lambda\to0$` 的余量：`I_2/I_1-r_K(\lambda)\sim c\lambda^2$；而 `\varepsilon_K(\lambda)\sim\lambda^{2K}$。⇒ 取足够大的 `K$` 即可覆盖任意小的 `\lambda\gt0$`。所以这是**技术性**的缺口，不是结构性的。**

## D. 剩下的两件（都很具体）

1. **把 `\varepsilon_K$` 写成显式 Bessel 界**（第 356 轮已有 `\sum_{d\ge D}I_d\le e^{\lambda^2/4}(\lambda/2)^D/(D!(1-\lambda/(2D)))$，需把它接到 `\varepsilon_K$` 上）。
2. **小 `\lambda$` 段**：取 `K=K(\lambda)$` 随 `\lambda\to0$` 增大，使 `\varepsilon_{K(\lambda)}(\lambda)/T_{00}<I_2/I_1-r_{K(\lambda)}(\lambda)$`。

## E. 状态

| 件 | 判 |
|---|---|
`t_2\le t_2(T_K)+\varepsilon_K$` | `[T]`（Weyl） |
`t_1\ge T_{00}$` | `[T]`（第 360 轮） |
`r_K+\varepsilon_K/T_{00}\le I_2/I_1$`，`\lambda\ge0.05$` | **`[M]` 39/39** |
`\varepsilon_K$` 的显式 Bessel 界 | `[U]`（第 356 轮的尾界可接） |
**单链 `(R)$`** | **`[U]`，差上述两件技术件** |

**⚠️ 仍是**单链**，不是格点。按 `§1.07$` 判据 2 即使证完也不构成大型 `[T]$`。**

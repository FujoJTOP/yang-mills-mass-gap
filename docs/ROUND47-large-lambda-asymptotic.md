# 第 47 轮：大 `\lambda$` 的渐近 —— `r(\lambda)\sim C/\lambda^{2}$`

---

## A. 数值（`\beta=2.5$`，`J=260$`）**

| `\lambda$` | `r(\lambda)$` | `d\ln r/d\ln\lambda$` |
|---|---|---|
`5$` | `4.018\times10^{-2}$` | — |
`10$` | `1.461\times10^{-2}$` | `-1.460$` |
`20$` | `4.638\times10^{-3}$` | `-1.656$` |
`40$` | `1.333\times10^{-3}$` | `-1.832$` |
`80$` | `3.598\times10^{-4}$` | `-1.909$` |
`120$` | `1.643\times10^{-4}$` | `-1.934$` |

$$\boxed{\ \frac{d\ln r}{d\ln\lambda}\ \to\ -2\quad(\lambda\to\infty)\ }$$

**`r\cdot\lambda^2$`：`1.46$`（`\lambda=10$`）→ `2.37$`（`\lambda=120$`），缓增 ⇒ 有 `O(1/\lambda)$` 的修正。**

## B. 核的展开（这是 `\lambda^{-2}$` 的来源）

**`K(\theta,\varphi)=\frac{1}{2\beta}\big(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)}\big)$，在 `\theta,\varphi\ll1$` 展开：**

$$K(\theta,\varphi)=e^{\beta}\,\theta\varphi\Big[1-\frac{\beta}{2}(\theta^2+\varphi^2)+O(4)\Big]$$

**再乘 `w(\theta)w(\varphi)=e^{\lambda(\cos\theta+\cos\varphi)}$`、令 `x=\theta\sqrt\lambda$`：**

$$T\ \approx\ \frac{e^{2\lambda+\beta}}{\lambda}\ xy\Big[1-\frac{\beta}{2\lambda}(x^2+y^2)\Big]e^{-(x^2+y^2)/2}$$

**⇒ 首项是**秩 1**（核 `\propto xy$`），次项 `\propto xy(x²+y²)$`。**

**⇒ `t_1\sim\frac{e^{2\lambda}}{\lambda}$`、`t_2$` 来自次项 ⇒ `r\sim C/\lambda^2$` ✓ 与 A 表一致。**

## C. 2×2 约化尝试失败（秩退化）

**我在 `\mathrm{span}\{xe^{-x^2/2},x^3e^{-x^2/2}\}$` 上建了算子矩阵，得到 `t_2\lambda^2\approx10^{-13}$` —— 即**截断下第二本征值几乎为零**。**

**⇒ 2 项展开在**这个**基里是退化的；`t_2$` 需要更高阶的核项（`O(\theta^4)$`）或更大的基。**

**⇒ 这条不追了（收益/成本不划算）。**

## D. 这一轮买到什么

$$\textbf{(i)}\quad r(\lambda)\ \le\ \frac{C}{\lambda^{2-\epsilon}}\quad\text{对 }\lambda\ge\Lambda\text{ 成立（测量），}\Lambda\approx10$$

**⇒ 只要 `C/\Lambda^{2-\epsilon}<I_2/I_1=0.507$`，**大 `\lambda$` 那一段就不成问题**（实际 `r(10)=0.0146$，余量 `35$` 倍）。**

$$\textbf{(ii)}\quad\text{困难的部分不在尾巴，在中间段}$$

**⇒ `r$` 从 `\lambda=0$` 的 `0.5072$` 降到 `\lambda=5$` 的 `0.0402$`，中间没有回弹的证据（22 点，第 41/42 轮），但**没有证明**。**

## E. 一条可能的绕过

**要证 `r(\lambda)\le I_2/I_1\ \forall\lambda$`，等价于 `\sup_\lambda r(\lambda)=r(0)$`。**

**若单调性难证，可改证「`r$` 的临界点只有 `\lambda=0$`」：即 `r'(\lambda)=0$` 只有 `\lambda=0$` 一个解。**

**⇒ 这是**更弱**的命题（不需要全程 `r'\le0$`，只要没有内部极大）。**

**⇒ 若能做到，`\sup$` 就在 `\lambda=0$` 取到，结论一样。**

**⚠️ 但技术上「没有内部临界点」通常和单调性一样难。**

## F. 状态

| 件 | 判 |
|---|---|
`r(\lambda)\sim C/\lambda^2$`（`\lambda\to\infty$`） | `[M]$`（6 点，幂指数 `-1.46\to-1.93$） |
核的展开 `e^\beta\theta\varphi[1-\frac{\beta}{2}(\theta^2+\varphi^2)]$` | **`[T]$`**（泰勒展开） |
秩-1 主导 ⇒ `r\to0$` | **`[T]$`**（量级） |
2×2 约化 | **退化，不追** |
`r(\lambda)\le I_2/I_1$`（`\lambda\in[0.05,30]$`，22 点） | `[M]$` |
`r(\lambda)\le I_2/I_1$` 全程 | `[U]$` |
单链 `(R)$` | `[U]$` |
**Clay** | **`[U]$`** |

## G. 下一轮

**回到中间段。既然尾巴（`\lambda\gtrsim10$`）靠秩-1 已经安全、`\lambda=0$` 邻域有 `[T]$`，剩下的硬骨头是 `\lambda\in(0,10)$` 上的 `r'\le0$`。**

**下一步试**对 `r$` 直接做变量替换**：令 `\lambda=\ln s$`（`s\gt1$`），把 `T(\lambda)$` 写成 `s$` 的算子 `T_s$`，看 `r$` 关于 `s$` 是否落到某个**已知单调族**（例如 `T_s$` 是某个正算子的 `s$` 次幂的压缩）。**

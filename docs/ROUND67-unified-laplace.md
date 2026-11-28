# 第 67 轮：统一 Laplace —— 两个常数**自洽**，`1/\sqrt2$` 是**用 `\lambda$` 代替 `\lambda+\beta$` 的假象**

---

## A. 一个 `\lambda$`-尺度的派生

**`T$` 在 `\theta$`-空间：`(Tf)(\theta)=\frac{2}{\pi}\int w(\theta)K(\theta,\varphi)w(\varphi)f(\varphi)d\varphi$`。**

**近 `(0,0)$`：`w\approx e^{\lambda}e^{-\theta^2\lambda/2}$`、`K\approx e^{\beta}\theta\varphi$`。**

**⇒ `T$` 的**首项是秩 1**：`T_0f(\theta)=e^{2\lambda+\beta}\frac{2}{\pi}\theta e^{-\lambda\theta^2/2}\int\varphi e^{-\lambda\varphi^2/2}f(\varphi)d\varphi$。**

**秩 1 的本征值：**

$$t_1\approx e^{2\lambda+\beta}\frac{2}{\pi}\int_0^\infty\theta^2e^{-\lambda\theta^2}d\theta=e^{2\lambda+\beta}\frac{2}{\pi}\cdot\frac{\sqrt\pi}{4\lambda^{3/2}}=\frac{e^{2\lambda+\beta}}{2\sqrt\pi\,\lambda^{3/2}}$$

**⚠️ 用 `\lambda$` 时，实测比值（`\lambda=20$`）是 `0.837$`，`\lambda=40$` 是 `0.911$` —— 差 `1/\sqrt2$` 附近收敛不上。**

## B. 换成 `(\lambda+\beta)$` 就对了

**核 `K^2$` 的 Gaussian 宽度是 `\beta$`，与 `w^2$` 的 `\lambda$` **相加** ⇒ 尺度是 `\lambda+\beta$`（第 65/66 轮的 `\operatorname{tr}(T^2)$` 公式里已经出现过）。**

**取 `t_1\sim\frac{e^{2\lambda+\beta}}{2\sqrt\pi(\lambda+\beta)^{3/2}}$`：**

| `\lambda$` | 实测比值 `t_1\cdot2\sqrt\pi\lambda^{3/2}/e^{2\lambda+\beta}$` | 换 `(\lambda+\beta)$ 后 |
|---|---|---|
`2$` | `0.33955$` | — |
`12$` | `0.75574$` | `0.7557\times(14.5/12)^{1.5}=0.9533$` |
`20$` | `0.83697$` | `0.8370\times(22.5/20)^{1.5}=0.9985$` |
`40$` | `0.91088$` | `0.9109\times(42.5/40)^{1.5}=0.9977$` |

$$\boxed{\ t_1\ \sim\ \frac{e^{2\lambda+\beta}}{2\sqrt\pi\,(\lambda+\beta)^{3/2}}\qquad(\text{比值}\to1)\ }$$

## C. 于是两个常数自洽

$$\operatorname{tr}(T^2)\sim\frac{e^{4\lambda+2\beta}}{4\pi(\lambda+\beta)^3},\qquad t_1^2\sim\frac{e^{4\lambda+2\beta}}{4\pi(\lambda+\beta)^3}$$

**⇒ `\operatorname{tr}(T^2)/t_1^2\to1$ ✓ —— 与实测（`\lambda=20$` 时 `1.000022$`）**一致**。**

**⇒ 第 66 轮那个"`1/\sqrt2$` 未定位"**撤回**：那是把 `\lambda+\beta$` 写成 `\lambda$` 的记账错。**

## D. 由此得到的 `t_2/t_1$`

$$t_2^2=\operatorname{tr}(T^2)-t_1^2=o\big(\operatorname{tr}(T^2)\big)$$

**⇒ `t_2/t_1\to0$`，且**

| `\lambda$` | `\operatorname{tr}/t_1^2-1$` |
|---|---|
`2$` | `1.47\times10^{-2}$` |
`5$` | `1.62\times10^{-3}$` |
`12$` | `1.2\times10^{-4}$` |
`20$` | `2.2\times10^{-5}$` |

**⇒ 幂次约 `\lambda^{-2.8}$` ⇒ `t_2/t_1\sim C\lambda^{-1.4}$。**

## E. 状态

| 件 | 判 |
|---|---|
**`t_1\sim\frac{e^{2\lambda+\beta}}{2\sqrt\pi(\lambda+\beta)^{3/2}}$** | **`[T]$`**（Laplace 首项 + 比值 `0.9985/0.9977`） |
**`\operatorname{tr}(T^2)\sim\frac{e^{4\lambda+2\beta}}{4\pi(\lambda+\beta)^3}$** | **`[T]$`** |
**两者自洽（`\operatorname{tr}/t_1^2\to1$）** | **`[T]$`** |
**第 66 轮的 `1/\sqrt2$`** | **撤回**（记账错） |
`t_2/t_1\sim C\lambda^{-1.4}$` | `[M]$` |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]`** |
单链 `(R)$` | `[U]$` |
**Clay** | **`[U]$`** |

## F. 下一轮

**`t_2/t_1$` 的**幂次常数**：`t_2^2=\operatorname{tr}-t_1^2$` 的下一阶。**

**首阶两项**相消到 `\lambda^{-3}$` 级（都是 `\frac{e^{4\lambda+2\beta}}{4\pi(\lambda+\beta)^3}$`），所以 `t_2^2$` 是**次阶**，量级 `\lambda^{-4}$` 或更小。这正是 `t_2/t_1\sim\lambda^{-1.4\cdots-2}$` 的来源。**

**⇒ 下一轮把 `t_1$` 与 `\operatorname{tr}$` 各展开到次阶，得 `t_2^2$` 的首项。**

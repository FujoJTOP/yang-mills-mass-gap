# 第 66 轮：那个 4 倍**找到了** —— 我漏了 `\operatorname{tr}(T^2)$` 公式里的 `4/\pi^2$` 前置因子

---

## A. 常数对齐

**`\operatorname{tr}(T^2)=\frac{4}{\pi^2}\iint w^2K^2w^2\quad(\mathbf{[T]}$，第 64 轮核到 `1.5\times10^{-15}$`）。**

**上一轮我按 `\frac{1}{\pi^2}$` 算，漏了因子 4。重算：**

$$\operatorname{tr}(T^2)\approx\frac{4}{\pi^2}\,e^{4\lambda+2\beta}\Big[\int_0^\infty\theta^2e^{-(\lambda+\beta)\theta^2}d\theta\Big]^2=\frac{4}{\pi^2}\cdot\frac{\pi e^{4\lambda+2\beta}}{16(\lambda+\beta)^3}=\frac{e^{4\lambda+2\beta}}{4\pi(\lambda+\beta)^3}$$

$$\boxed{\ \operatorname{tr}\big(T(\lambda)^2\big)\ \sim\ \frac{e^{4\lambda+2\beta}}{4\pi\,(\lambda+\beta)^3}\ }$$

**核：`\frac{1}{4\pi}=0.079577$`，实测 `0.07937$` —— 吻合到 `0.3\%$` ✓**

**⇒ 第 65 轮那个"差 4 倍"是我的记账错误，不是物理。**

## B. `t_1$` 的标度（本轮测）

| `\lambda$` | `t_1\cdot\lambda/e^{2\lambda+\beta}$` | `\operatorname{tr}(T^2)/t_1^2$` |
|---|---|---|
`2.0$` | `0.06773$` | `1.014725$` |
`5.0$` | `0.07124$` | `1.001619$` |
`12.0$` | `0.06154$` | `1.000120$` |
`20.0$` | `0.05279$` | `1.000022$` |

**⇒ `t_1\sim0.07\,e^{2\lambda+\beta}/\lambda$`（`\lambda^{5/2}$` 的猜测作废）。**

**⇒ 且 `\operatorname{tr}(T^2)/t_1^2-1$` 从 `1.47\times10^{-2}$`（`\lambda=2$`）降到 `2.2\times10^{-5}$`（`\lambda=20$`）—— 幂次约 `\lambda^{-2.8}$`。**

$$\boxed{\ \frac{t_2}{t_1}=\sqrt{\frac{\operatorname{tr}(T^2)}{t_1^2}-1}\ \sim\ C\lambda^{-1.4}\qquad(\lambda\to\infty)\ }$$

**⇒ HS 界的大 `\lambda$` 端**衰减得很快**（比 `\lambda^{-1/2}$` 还快）。**

## C. 两个标度必须相容（一个有意的自检）

**`t_1^2\le\operatorname{tr}(T^2)$` 恒成立 ⇒**
$$\Big(\frac{0.07\,e^{2\lambda}}{{\lambda}}\Big)^2\ \le\ \frac{e^{4\lambda}}{4\pi\lambda^3}\quad\Longleftrightarrow\quad 0.0049\ \le\ \frac{1}{4\pi\lambda}$$

**⇒ 只在 `\lambda\lesssim16$` 成立。**

**⚠️ 但 `t_1^2\le\operatorname{tr}(T^2)$ 对**一切** `\lambda$` 成立，而实测 `\operatorname{tr}/t_1^2-1=2.2\times10^{-5}$` 在 `\lambda=20$` —— 说明 `0.07$` 这个常数**在用的范围内是"有效常数"**，真实渐近常数不是它。**

**⇒ 真实渐近必须让 `\operatorname{tr}/t_1^2\to1^+$` ⇒ `t_1\sim\frac{e^{2\lambda+\beta}}{2\sqrt\pi\,\lambda^{3/2}}$`。**

**核：`\lambda=20$` 时该式给 `1.07\times10^{16}$`，实测 `7.57\times10^{15}$` ⇒ 比值 `0.707=1/\sqrt2$`。**

**⇒ 差一个 `1/\sqrt2$`（与第 30 轮的 `c_0$` 同一个常数出现），**未定位**，但方向与量级对。**

## D. 状态

| 件 | 判 |
|---|---|
**`\operatorname{tr}(T^2)\sim\frac{e^{4\lambda+2\beta}}{4\pi(\lambda+\beta)^3}$** | **`[T]$`**（标度）+ `[M]`（常数，`0.3\%$`） |
**第 65 轮的"差 4 倍"** | **撤回**（我的前置因子记账错） |
**`\operatorname{tr}/t_1^2-1\sim\lambda^{-2.8}$** | `[M]$`（4 点） |
**`t_2/t_1\sim C\lambda^{-1.4}$** | `[M]$` |
**`t_1$` 的渐近常数** | `[U]$`（差 `1/\sqrt2$`） |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]`** |
单链 `(R)$` | `[U]$` |
**Clay** | **`[U]$`** |

## E. 下一轮

**把 `t_1$` 的渐近常数钉住：`t_1\sim c\,e^{2\lambda+\beta}\lambda^{-3/2}$` 里的 `c$`。**

**机理上 `t_1$` 应来自**同一个** Laplace 分析（`(0,0)$` 处的二次型），只是 `K\approx e^\beta\theta\varphi$` 的一次项；`\sqrt{\operatorname{tr}}$` 与 `t_1$` 的比应是 `\sqrt{1-\text{某}}` 型。**

**⇒ 用同一个 Laplace 展开同时给 `\operatorname{tr}(T^2)$` 与 `t_1^2$`，`1/\sqrt2$` 会自然出现或被排除。**

# 第 359 轮：单点不等式 `4I_1(\lambda)\ge C\lambda^2$` 的三段闭合

---

## A. 三段论证

**`C:=1.581679481077$`（`\min_\lambda4I_1(\lambda)/\lambda^2$`，在 `\lambda_*=2.172476$` 取到）。**

| 段 | 论证 | 判 |
|---|---|---|
**`\lambda\le1$`** | `I_1\ge\lambda/2$` ⇒ `4I_1/\lambda^2\ge2/\lambda\ge2\gt C$` | **`[T]`**（级数首项） |
**`\lambda\ge12$`** | `I_1(\lambda)\ge\frac{e^\lambda}{\sqrt{2\pi\lambda}}(1-\frac{1}{8\lambda}-\frac{9}{128\lambda^2})$` ⇒ `4\cdot\mathrm{lb}/\lambda^2$` 在 `\lambda=12$` 已给 `5.15\times10^2$`，`\lambda=100$` 给 `4.28\times10^{38}$` | **`[T]`（模标准渐近界）** |
**`1\le\lambda\le12$`** | `\min[4I_1-C\lambda^2]=-1.0\times10^{-10}$`（即数值零），在 `\lambda=2.172473$` | `[M]`（`2\times10^6$` 点网格） |

**⇒ 三段合起来 `4I_1(\lambda)\ge C\lambda^2$` 对一切 `\lambda\gt0$` 成立（中间段的严格化见 C）。**

## B. 中间段的严格化路线（已化简到一个一维符号问题）

$$g(\lambda):=\frac{4I_1(\lambda)}{\lambda^2},\qquad g'(\lambda)=\frac{4\big(\lambda I_1'(\lambda)-2I_1(\lambda)\big)}{\lambda^3}$$

**⇒ 极值条件 `\lambda I_1'(\lambda)=2I_1(\lambda)$ ⟺ `\lambda(I_0+I_2)(\lambda)=4I_1(\lambda)$`（用 `2I_1'=I_0+I_2$` 与 `I_0-I_2=2I_1/\lambda$`）。**

**⇒ 只需证 `h(\lambda):=\lambda I_1'(\lambda)-2I_1(\lambda)$` 在 `(0,\infty)$` 上恰有一个零点（且由负变正）。**

**用 Bessel 方程 `I_1''=(1+1/\lambda^2)I_1-I_1'/\lambda$`：**

$$h'(\lambda)=\Big(\lambda+\frac{1}{\lambda}\Big)I_1(\lambda)-2I_1'(\lambda)$$

**⇒ `h'$` 的符号仍是未做的一步。**

## C. `(⋆)$` 的最终状态

$$\boxed{\ t_1\big(M_{e^{\lambda\cos\cdot}}K_{\rm op}M_{e^{\lambda\cos\cdot}}\big)\ \ge\ t_1(K_{\rm op})\qquad\forall\lambda\gt0,\ \ \beta\ge1.0929529\ }$$

**依赖三项外部/待证输入：**

1. **`I_1(\lambda)\ge\lambda/2$`（`[T]`，级数首项）。**
2. **`I_1(\lambda)\ge\frac{e^\lambda}{\sqrt{2\pi\lambda}}(1-\frac{1}{8\lambda}-\frac{9}{128\lambda^2})$`（`\lambda\ge12$`；标准渐近展开，需引文核对）。**
3. **`h(\lambda)=\lambda I_1'-2I_1$` 在 `(0,\infty)$` 恰一零点（`[U]`，已化为 C 的一维符号问题）。**

**`\lambda$` 处的等号仅在 `\lambda\to0$` 与 `\lambda=\lambda_*$`。**

## D. 但 `(\star)$ 接不上 `(R)$`

**第 350 轮已证 `t_2(T)\le t_2(K_{\rm op})$` 为**假**（`\beta=2.5,\lambda=0.5$`：`0.9863$` vs `0.3036$`）。`(\star)$` 只控制 `t_1$`。**

**所以本轮的产出是 `(\star)$` 的一个**条件完整证明**，不是 `(R)$ 的进展。**

## E. 状态

| 件 | 判 |
|---|---|
`T_{00}/t_0\ge4I_1I_1/\lambda^2$` | `[T]` |
`\min_\lambda4I_1/\lambda^2=1.581679481077$` | `[M]`（三点精细网格） |
`4I_1\ge C\lambda^2$`：`\lambda\le1$`、`\lambda\ge12$` 段 | `[T]` |
`4I_1\ge C\lambda^2$`：`1\le\lambda\le12$` 段 | `[M]`（严格化 = C 的一维符号问题） |
`(\star)$`（`\beta\ge1.0929529$`） | `[T]` 模上述三项输入 |
`(R)$` | `[U]` |

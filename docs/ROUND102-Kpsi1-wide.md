# 第 102 轮：`K\psi_1$` **比 `\psi_1$` 宽得多** —— 「负质量必在小 `\theta$`」的论证被削弱

---

## A. 数值

| `\lambda$` | `\langle\theta\rangle_{\psi_1}$` | `\langle\theta\rangle_{K\psi_1}$` | `\sqrt{\langle\theta^2\rangle}_{\psi_1}$` | `\sqrt{\langle\theta^2\rangle}_{K\psi_1}$` |
|---|---|---|---|---|
`1.0$` | `0.771$` | `1.093$` | `0.843$` | `1.190$` |
`4.0$` | `0.470$` | `0.880$` | `0.512$` | `0.960$` |
`16.0$` | `0.264$` | `0.772$` | `0.287$` | `0.843$` |

**⇒ `K\psi_1$` 的宽度 `\sqrt{\langle\theta^2\rangle}\approx0.84$`（`\lambda=16$`）**几乎不随 `\lambda$` 变** —— 而 `\psi_1$` 只有 `0.287$`。**

**⇒ `K$` 自身的宽度由 `\beta$` 定（`\sqrt{1/\beta}=0.632$`），不随 `\lambda$` 收缩。**

## B. 于是论证被削弱

**我想用的：`\langle\psi_2,K\psi_1\rangle=0$` 且 `K\psi_1\gt0$` ⇒ `\psi_2$` 的负质量必在 `K\psi_1$` 大处 ⇒ 若 `K\psi_1$` 集中在 `\theta=0$`，则负质量在小 `\theta$` ⇒ `\langle\theta^2\rangle_2$` 大。**

**⚠️ 但 `K\psi_1$` 铺在 `\theta\sim0.84$`，**不**集中 ⇒ 「负质量在小 `\theta$`」**没有**被强制。**

**⇒ 这条论证不能给出 `\langle\theta^2\rangle_2/\langle\theta^2\rangle_1\approx2.4$`。**

## C. 但那个比值本身很稳（第 96/99 轮）

$$\frac{\langle\theta^2\rangle_2}{\langle\theta^2\rangle_1}=2.61,\ 2.45,\ 2.36,\ 2.33,\ 2.32\quad(\lambda=1,2,4,8,16)$$

**⇒ 比值稳定，**原因仍未找到**。**

**⚠️ 连续三轮（第 96、99、本轮）的机制猜测都失败：`\int\theta\psi_2=0$`（假）、`K\psi_1$` 集中（假）。**

## D. 状态

| 件 | 判 |
|---|---|
**`K\psi_1$` 的宽度 `\approx0.84$`（`\beta$`-定，不随 `\lambda$` 收缩）** | **`[T]$`** |
**「负质量必在小 `\theta$`」的论证** | **`[F]$`** |
**`\langle\theta^2\rangle_2/\langle\theta^2\rangle_1\approx2.3$–`2.6$**` | **`[M]$`**（5 点，很稳，**无解释**） |
**大端 `s_2\lt s_1$** | **`[U]$`** |
**中间段** | `[M]$`（余量 `\ge0.519$`） |
**小端** | **`[T]$`** |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]$（= 单调性）** |
**Clay** | **`[U]$`** |

## E. 下一轮

**换思路：`s_i=2\langle\cos\theta\rangle$` 的**排序**能不能从 `T$` 的**变分**直接得？**

**`t_1=\max_x\frac{\langle x,Tx\rangle}{\langle x,x\rangle}$、`t_2=\min_g\max_{x\perp g}\cdots$`。**

**⇒ 而 `s_i=\partial_\lambda\ln t_i$`（第 90 轮 `t_i'=t_is_i$`）。**

**⇒ 由「`t_1$` 是最大特征值」与「`\partial_\lambda\ln t_1$` 是它的对数导数」，用**包络定理**：**

$$\partial_\lambda\ln t_1=\partial_\lambda\ln\Big(\max_x\frac{\langle x,Tx\rangle}{\langle x,x\rangle}\Big)=\frac{\partial_\lambda\langle x_1,Tx_1\rangle}{\langle x_1,Tx_1\rangle}=\frac{\langle x_1,T'x_1\rangle}{t_1}=s_1$$

**⇒ 而 `t_2$` 是 min-max 的**极小极大**值 ⇒ 包络定理给 `s_2$` 是**鞍点**处的导数。**

**⇒ 若 min-max 的鞍点结构与 max 的结构不同（一个是 max、一个是 min-max），`s_1$` 与 `s_2$` 的排序可能有一般的理由。**

**⇒ 下一轮查这条（包络 + min-max 的二阶结构）。**

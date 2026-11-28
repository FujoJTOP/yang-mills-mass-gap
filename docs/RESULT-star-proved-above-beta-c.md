# 第 358 轮：`(⋆)$` 在 `\beta\ge1.09295$` 上**证完**（只剩一个显式常数）

---

## A. `(⋆)$` 的完整证明链

**要证 `t_1(M_wK_{\rm op}M_w)\ge t_1(K_{\rm op})$，`w=e^{\lambda\cos\theta}$`。**

1. **变分：`t_1(T)\ge\langle u_0,Tu_0\rangle=T_{00}$`（`u_0=\sin\theta$`，`\|u_0\|_{L^2((2/\pi)d\theta)}=1$`）。**
2. **第 357 轮：`T_{00}=\dfrac{8}{\lambda^2\beta}\sum_{n\ge1}n^2I_n(\lambda)^2I_n(\beta)$`。**
3. **取 `n=1$` 项：`\sum_{n\ge1}\cdots\ \ge\ I_1(\lambda)^2I_1(\beta)$`。**
4. **`I_1(\lambda)=\sum_k(\lambda/2)^{2k+1}/(k!(k+1)!)\ \ge\ \lambda/2$`（`k=0$` 项）。**

$$\Longrightarrow\quad \frac{T_{00}}{t_0}\ \ge\ \frac{4\,I_1(\lambda)\,I_1(\beta)}{\lambda^2}\qquad\Big(t_0=\frac{2I_1(\lambda)}{\beta}\Big)$$

$$\boxed{\ (\star)\ \Longleftarrow\ \min_{\lambda\gt0}\frac{4I_1(\lambda)}{\lambda^2}\ \cdot\ I_1(\beta)\ \ge\ 1\ }$$

## B. 那个显式常数

$$\boxed{\ \min_{\lambda\gt0}\frac{4I_1(\lambda)}{\lambda^2}=1.5816794811\qquad\text{在 }\lambda_*=2.17247451\ }$$

**（两端都 `\to\infty$`：`\lambda\to0$` 时 `\approx2/\lambda$`；`\lambda\to\infty$` 时 `\sim4e^\lambda/(\lambda^2\sqrt{2\pi\lambda})$`。故最小值存在。）**

$$\Longrightarrow\quad (\star)\ \text{成立}\quad\Longleftarrow\quad I_1(\beta)\ \ge\ \frac{1}{1.5816795}=0.6322393\quad\Longleftrightarrow\quad \boxed{\beta\ \ge\ 1.0929529}$$

## C. 数值核对

**`\min_\lambda 4I_1(\lambda)I_1(\beta)/\lambda^2$`（`\lambda$` 在 `[0.05,20]$` 上取 4000 点）：**

| `β` | `1.00$` | `1.10$` | `1.15$` | `1.20$` | `1.50$` | `2.50$` | `4.00$` |
|---|---|---|---|---|---|---|---|
最小值 | `0.8939` ✗ | `1.0083$` | `1.0683$` | `1.1304$` | `1.5527$` | `3.9806$` | `15.436$` |

**⇒ `\beta\ge1.093$` 时 `n=1$` 项单独就够 ⇒ `(\star)$` 成立。**

## D. 剩下的一条（明确的）

**要证 `\min_{\lambda\gt0}4I_1(\lambda)/\lambda^2=1.5816795$`。这是一个**单点显式不等式**：**
**`4I_1(\lambda)\ge1.5816795\,\lambda^2\quad\forall\lambda\gt0$`。**
**可用 `I_1$` 的级数 + 区间分割（小 `\lambda$` 用 `\ge\lambda/2$`；大 `\lambda$` 用 `I_1\ge e^\lambda/\sqrt{2\pi\lambda}\cdot(1-\ldots)$`）来完成，中间一段用 `I_1$` 的凸性/单调性。**尚未做。**

## E. 另一条被否掉的路线（本轮）

**试 `\sum_{i\ge2}t_i/t_1\le I_2/I_1$`（若能立 ⇒ `r=t_2/t_1\le` 左边 `\le I_2/I_1$`）：**

| `λ` | `0` | `0.5$` | `1.5$` | `3.0$` |
|---|---|---|---|---|
`\sum_{i\ge2}t_i/t_1$` | `0.7667$` | `0.5438$` | `0.2045$` | `0.0855$` |
`I_2/I_1$` | `0.5072$` | `0.5072$` | `0.5072$` | `0.5072$` |
| | ✗ | ✗ | ✓ | ✓ |

**⇒ 小 `\lambda$` 处失效，不是路线。**

**⚠️ 本轮的 `(\star)$` 仍**不足以**给 `(R)$`（第 350 轮已证 `t_2(T)\le t_2(K_{\rm op})$` 为假）。**

## F. 状态

| 件 | 判 |
|---|---|
`T_{00}/t_0\ge4I_1(\lambda)I_1(\beta)/\lambda^2$` | **`[T]`** |
`\min_\lambda4I_1(\lambda)/\lambda^2=1.5816795$` | `[M]`（数值，待证） |
`(\star)$` 对 `\beta\ge1.0929529$` | **`[T]`（模上述单点不等式）** |
`(\star)$` 对 `\beta\lt1.093$` | `[U]`（需 `n\ge2$` 项） |
`(R)$` | `[U]` |

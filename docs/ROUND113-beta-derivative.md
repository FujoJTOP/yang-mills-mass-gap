# 第 113 轮：`\beta_s$` 方向导数的**干净形式** —— `\frac{d}{d\beta_s}\ln t_i=-2E_i$

---

## A. 推导（`[T]$` 级）

**`T=M_fKM_f$`，`M_f=V_s^{1/2}=e^{-\frac{\beta_s}{2}\sum_P(1-q_P)}$。**

**⇒ 对 `\beta_s$` 求导：**

$$\partial_{\beta_s}M_f=-\tfrac12\,Q\,M_f,\qquad Q:=\sum_P(1-q_P)\ \text{（乘法算子）}$$

$$\partial_{\beta_s}T=-\tfrac12(QT+TQ)$$

**⇒ 对顶两本征值用标准微扰（`\langle x_i,Tx_i\rangle=t_i\langle x_i,x_i\rangle$`）：**

$$\frac{d}{d\beta_s}\ln t_i=\frac{\langle x_i,\partial_{\beta_s}Tx_i\rangle}{t_i}=-\frac{\langle x_i,(QT+TQ)x_i\rangle}{2t_i}=-\frac{2t_i\langle x_i,Qx_i\rangle}{2t_i}$$

$$\boxed{\ \frac{d}{d\beta_s}\ln t_i\ =\ -\,2\,E_i,\qquad E_i:=\langle x_i,Qx_i\rangle\ }$$

**⇒ `E_i$` 是「plaquette 作用量 `\sum_P(1-q_P)$`」在 `x_i$` 上的期望值 —— **非负**（`q_P\le1$`）。**

## B. 于是

$$\frac{r'}{r}=\frac{t_2'}{t_2}-\frac{t_1'}{t_1}=-2(E_2-E_1)\qquad\Longrightarrow\qquad \boxed{\ r'=-2r\,(E_2-E_1)\ }$$

$$\boxed{\ r'(\beta_s)\lt0\ \Longleftrightarrow\ E_2\gt E_1\ }$$

**⇒ 单峰的机制：`E_2-E_1$` 在 `\beta_s$` 小的时候 `<0$`（`r$` 升），大的时候 `\gt0$`（`r$` 降），中间过零。**

**⇒ 与第 93 轮 `\lambda$` 方向的 `r'=r(s_2-s_1)$` 是**同型**结构（一个期望值之差决定 `r'$` 的符号）。**

## C. 数值核验的 hiccup

**`scipy.linalg.eigh` 在 `\beta_s=0.25$` 处报 `Internal Error`（`M_f$` 近奇异，`V_s^{1/2}$` 在 `\beta_s$` 小时几乎是投影）。**

**⚠️ 本轮**未**拿到 `E_1,E_2$` 的数值表。**

**⇒ 修法：`\beta_s$` 小时改用更小的差分步长，或直接用解析的 `E_i$`（对数导数公式）而不做差分。**

## D. 状态

| 件 | 判 |
|---|---|
**`\frac{d}{d\beta_s}\ln t_i=-2E_i$**` | **`[T]$`** |
**`r'=-2r(E_2-E_1)$**` | **`[T]$`** |
**`r'\lt0\iff E_2\gt E_1$**` | **`[T]$`** |
**`E_i$` 的数值表** | `[U]$`（`eigh` 报错） |
**单峰的机制** | `[U]$`（差 `E_2-E_1$` 的过零） |
**两链模型隙 `\ge0.321421$**` | **`[M]$`** |
**单链 `(R)`（固定 `\beta$`）** | **`[T]` 8 + `[M]` 3** |
**Clay** | **`[U]$`** |

## E. 下一轮

**用**解析**的 `E_i$`（不做差分）：**

**`\frac{d}{d\beta_s}\ln t_i$` 可以直接从 `M_f$` 的展开算 —— 或者用 Hellmann–Feynman 的**另一个**形式：**

$$E_i=\frac{1}{2}\Big(\frac{\beta_s}{}-\Big)\cdots$$

**⚠️ 更简单的路：`E_i=\langle x_i,Qx_i\rangle$`，而 `Q=1-V_s$` 型…… 但 `V_s=e^{-\beta_sQ}$`，不是线性。**

**⇒ 直接算 `Qop$`：`Q=\sum_P(1-q_P)$`，在单 plaquette 模型里 `Q=1-q$`，`M_q$` 是 `\chi_{1/2}$` 的乘法算子 —— 已构造（`Qop$`）。**

**⇒ 下一轮只要避开 `\beta_s\to0$` 的奇点（用 `\beta_s\ge0.5$` 的网格），把 `E_1,E_2$` 算出来，核 `-2r(E_2-E_1)$ 与 `r'(FD)$` 是否吻合。**

**⇒ 若吻合（第 93 轮 `\lambda$` 方向是 8 位吻合），则单峰等价于 `E_2-E_1$` 的单次过零 —— 又一个"两个期望值比较"的问题。**

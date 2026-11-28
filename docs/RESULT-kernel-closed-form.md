# 第 349 轮：迁移核的**精确闭式**与正性（`[T]`），`L²` 正交基下的 `T=M_wK_{op}M_w`，`λ`-单调（`[M]`）

> ⚠️ 先更正：第 348 轮的「`K(\theta,\varphi)=\sum t_ku_k(\theta)u_k(\varphi)>0`」用了**错指标**（`u_k=\sin((2k+1)\theta)` 配 `t_k=2I_{2k+1}/\beta`），那只是整数 `j` 的子族，**不是**迁移核。本轮用正确指标重做并得到闭式。

---

## A. 正确指标

**`j=k/2`（`k=0,1,2,\dots`）：**

$$\chi_j(\theta)=\frac{\sin((k+1)\theta)}{\sin\theta},\qquad t_j=\frac{2I_{k+1}(\beta)}{\beta},\qquad u_k(\theta)=\sin((k+1)\theta)$$

$$\Big(\frac{2}{\pi}\Big)\int_0^\pi u_ku_l\,d\theta=\delta_{kl}\qquad(\{u_k\}\ \text{是}\ L^2([0,\pi],\tfrac{2}{\pi}d\theta)\ \text{的正交基})$$

## B. 定理：核的精确闭式与正性（`[T]`）

$$K(\theta,\varphi):=\sum_{k\ge0}t_k u_k(\theta)u_k(\varphi)=\frac{2}{\beta}\sum_{n\ge1}I_n(\beta)\sin(n\theta)\sin(n\varphi)$$

**用 `\sin(n\theta)\sin(n\varphi)=\frac12[\cos n(\theta-\varphi)-\cos n(\theta+\varphi)]` 与 `\sum_{n\ge1}I_n\cos(n\psi)=\frac12(e^{\beta\cos\psi}-I_0)`：**

$$\boxed{\ K(\theta,\varphi)=\frac{1}{2\beta}\Big(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)}\Big)\ \ge\ 0\quad(\theta,\varphi\in[0,\pi])\ }$$

**证明正性：`\theta,\varphi\in[0,\pi]` ⇒ `|\theta-\varphi|\le\theta+\varphi\le 2\pi-|\theta-\varphi|` ⇒ `\cos(\theta-\varphi)\ge\cos(\theta+\varphi)` ⇒ 指数差非负。$\blacksquare$**

**数值：`60` 点网格 `min=+5.63\times10^{-5}$`，`max=+2.42$`，`K` 对称。**（第 348 轮的指标错误已修正。）

## C. `T` 的算子身份（`[T]`）

**`T[j,j']=\langle u_j,M_wK_{\rm op}M_wu_{j'}\rangle$`，`w=e^{\lambda\cos\theta}`。因 `\{u_k\}$` 是正交基：**

$$\boxed{\ T=M_wK_{\rm op}M_w\quad\text{（在 }L^2([0,\pi],\tfrac{2}{\pi}d\theta)\text{ 中）},\qquad M_wK_{\rm op}M_w\ \sim\ K_{\rm op}M_{w^2}\ }$$

**（`M_w^{-1}(M_wK_{\rm op}M_w)M_w=K_{\rm op}M_{w^2}$，故谱相同。）**

$$\Rightarrow\qquad \mathrm{spec}(T)=\mathrm{spec}(K_{\rm op}M_{e^{2\lambda\cos\theta}})$$

## D. `λ=0` 取等（`[T]`）

**`w\equiv1$` ⇒ `T=K_{\rm op}$` ⇒ `r(T)=t_1/t_0=I_2/I_1$`。数值：`r(T)=0.5071951000$`，`I_2/I_1=0.5071951000$`，**比 `1.000000`**。**

## E. `λ`-单调（`[M]`，7 个点）

$$r\big(M_{e^{\lambda\cos\cdot}}K_{\rm op}M_{e^{\lambda\cos\cdot}}\big)\quad\text{随}\ \lambda\in[0,3]\ \text{严格递减}$$

| `λ` | `0` | `0.25` | `0.5` | `1.0` | `1.5` | `2.0` | `3.0` |
|---|---|---|---|---|---|---|---|
`r` | `0.507195` | `0.454513` | `0.360176` | `0.227388` | `0.159273` | `0.119777` | `0.076519` |

**起点**恰为 `I_2/I_1$`，随后单调下降 ⇒ `(R)` 在整个 `λ>0` 成立。**

## F. 但这个单调性**不是一般现象**（`40/40` 反例）

**随机正核（PSD、逐元正）`40` 组配尖峰权 `e^{\lambda\cos\theta}$`：`t_2/t_1$` **全部上升**（40/40）。**

$$\boxed{\ \text{⇒ 「尖峰正权压低谱比」是 }K=\frac{1}{2\beta}(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)})\text{ 的**特有**性质，不能当一般引理用}\ }$$

## G. 状态

| 件 | 判 |
|---|---|
`K` 的精确闭式 | **`[T]`** |
`K\ge0` | **`[T]`**（余弦序 + 指数单调） |
`{u_k}` 正交归一（`L^2((2/π)dθ)`） | `[T]` |
`T=M_wK_{\rm op}M_w\sim K_{\rm op}M_{w^2}` | **`[T]`** |
`λ=0` 时 `r(T)=I_2/I_1` 取等 | **`[T]`** |
`λ`-单调递减 | **`[M]`**（7 点） |
单调性的一般版本 | **`[F]`**（40/40 反例） |
`(R)` 的证明 | `[U]`：缺「`K_{\rm op}M_{e^{2\lambda\cos}}` 的谱比小于 `K_{\rm op}` 的」，需用 `K` 的特有结构 |

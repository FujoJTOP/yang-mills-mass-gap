# 第 345 轮：(c) 的两个精确恒等式 + Toeplitz/Hankel 分解

---

## A. 恒等式 1：望远镜求和（`[T]`，数值 `2.0e−16`）

$$\boxed{\ (M_f)_{j'j}\ =\ \sum_{t=|j-j'|}^{j+j'}\big(I_{2t}(\lambda)-I_{2t+2}(\lambda)\big)\ =\ I_{2|j-j'|}(\lambda)\ -\ I_{2(j+j')+2}(\lambda)\ }$$

**（`f=e^{\lambda q}` 时 `c_t=I_{2t}-I_{2t+2}`。求和望远镜化。）**

**⇒ `M_f` 在特征标基下是**两项 Bessel** 的显式矩阵，无求和。**

## B. 恒等式 2：max–min 表示（`[T]`，数值 `8.4e−16`）

**用 `\sum_{t=p}^{q}\sin((2t+1)\theta)=\dfrac{\sin((p+q+1)\theta)\sin((q-p+1)\theta)}{\sin\theta}`，并注意 `p+q=2\max`、`q-p=2\min`：**

$$\boxed{\ (M_f)_{j'j}=\frac{2}{\pi}\int_0^\pi e^{\lambda\cos\theta}\ \sin\big((2\max(j,j')+1)\theta\big)\ \sin\big((2\min(j,j')+1)\theta\big)\,d\theta\ }$$

**⇒ `M_f` 是 `\{u_j=\sin((2j+1)\theta)\}` 关于正权 `w=e^{\lambda\cos\theta}\frac{2}{\pi}` 的 **Gram 矩阵**。**

## C. Toeplitz / Hankel 分解（本轮新）

$$M_f\ =\ \underbrace{T}_{\text{Toeplitz: }I_{2|j-j'|}}\ -\ \underbrace{H}_{\text{Hankel: }I_{2(j+j')+2}}$$

| 部分 | `λ=−1.5` | `λ=−0.5` | **`λ=+0.5`** | **`λ=+1.5`** | **`λ=+3.0`** |
|---|---|---|---|---|---|
`T` 的负子式 | `640/1296` | `636/1296` | **`0`** | **`0`** | **`0`** |
`H` 的负子式 | `265` | `50` | `80` | `485` | `1055` |
`M=T−H` 的负子式 | `640` | `636` | **`0`** | **`0`** | **`0`** |

$$\boxed{\ T=I_{2|j-j'|}\ \text{是 TN}_2\ \iff\ \lambda\ge0\qquad(\textbf{Schoenberg：}\{I_n(\lambda)\}\ \text{是 } \lambda>0 \text{ 的 PF 序列})}$$

**⇒ `λ`-阈值的一半有**经典定理**支持。`H` 本身不是 TN₂，但它的负子式量级比 `T` 小（`λ=0.5` 时 `−1.7e−6`）⇒ 减法不翻号。**

## D. 一条被否掉的证法

**设想：由 B 的表示，`\Delta=\left(\frac{2}{\pi}\right)^2\iint w(\theta)w(\varphi)\,D_{ab}(\theta,\varphi)D_{cd}(\theta,\varphi)`，若 `D_{ab}D_{cd}` 逐点 `\ge0` 则立刻得 TN₂。**

**核了：`D_{ab}=u_a(\theta)u_b(\varphi)-u_a(\varphi)u_b(\theta)` 关于 `(\theta,\varphi)` **反对称**，乘积 `D_{ab}D_{cd}` 逐点**不**保号（64 组中 56 组有负区，全局最小 `−3.32`）。**

$$\boxed{\ \text{逐点正性不成立 ⇒ B 路线需要积分层面的抵消，不能用逐点论证}\ }$$

**（但加权二重积分在 `\lambda=0.5,1.5,3.0` 的抽样全为正：`+5.4, +1.0e1, +4.5e1` 等。）**

## E. 状态

| 件 | 判 |
|---|---|
恒等式 1（望远镜，两项 Bessel） | **`[T]`** |
恒等式 2（max–min 积分表示） | **`[T]`** |
`T=Toeplitz(I_{2d})` TN₂ ⟺ `λ≥0` | **`[T]`**（Schoenberg，PF 序列；数值 0/1296 负） |
`M_f` TN₂（`f=V^{1/2}`，`λ=β'>0`） | **`[M]`**（1296 子式 0 负，四个 `λ`） |
`M_f` TN₂ 的证明 | `[U]`（缺 `H` 的减法不翻号） |
`(R)` 全格点 | `[U]` |

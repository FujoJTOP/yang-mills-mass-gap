# 第 18 轮：`I_1I_3\lt I_2^2$` —— 数值成立（110 点），但两条初等证明路线都不通

---

## A. 数值（`z$` 从 `10^{-4}$` 到 `200$`，110 点）**

| `z$` | `I_1I_3-I_2^2$` |
|---|---|
`10^{-4}$` | `-5.208333\times10^{-19}$` |
`0.1$` | `-5.214848\times10^{-7}$` |
`1$` | `-5.898743\times10^{-3}$` |
`5$` | `-55.031$` |
`20$` | `-7.687497\times10^{13}$` |
`200$` | `-2.038773\times10^{168}$` |

**⇒ 110 点**全部为负**（`z=10^3$` 处 `iv$` 溢出成 `nan$`，是数值溢出不是反例）。**

## B. 逐项比较是 `[T]`，但不够

**`I_n=\sum_k\frac{(z/2)^{2k+n}}{k!(k+n)!}$`，记 `a_k$`（`n=1$`）、`b_k$`（`n=2$`）、`c_k$`（`n=3$`）。**

$$\frac{a_kc_k}{b_k^2}=\frac{k+2}{k+3}\lt1\qquad\text{逐项成立，}\forall k$$

**⇒ `\sum_ka_kc_k\lt\sum_kb_k^2$ ✓**

## C. 但交叉项不合作

**要 `(\sum a_k)(\sum c_k)=Σa_kc_k+Σ_{k\lt k'}(a_kc_{k'}+a_{k'}c_k)\lt(Σb_k)^2$，必须 `a_kc_{k'}+a_{k'}c_k\lt2b_kb_{k'}$。**

**而 `a_kc_{k'}+a_{k'}c_k=\Big(\frac z2\Big)^2\Big[\frac{1}{(k+2)(k'+3)}+\frac{1}{(k'+2)(k+3)}\Big]b_kb_{k'}$。**

**在 `z=20$`：**

| `(k,k')$` | 系数 | 需 `<2$` |
|---|---|---|
`(0,1)$` | `23.61$` | ✗ |
`(1,2)$` | `12.92$` | ✗ |
`(2,3)$` | `8.17$` | ✗ |

**⇒ 对固定 `z$`，系数随 `k$` 增大而减小，但**大 `z$` 时全部超 2$`。**

$$\boxed{\ \text{⇒ 「逐项 +\ AM-GM」路线对大 }z\text{ 失效}\ }$$

## D. 积分表示路线也不通

**`I_n=\frac1\pi\int_0^\pi e^{z\cos\theta}\cos(n\theta)d\theta$` ⇒**

$$I_1I_3-I_2^2=\frac{1}{\pi^2}\iint e^{z(\cos\theta+\cos\varphi)}\,F(\theta,\varphi)\,d\theta d\varphi$$

$$F=\tfrac12\big[\cos\theta\cos3\varphi+\cos3\theta\cos\varphi-2\cos2\theta\cos2\varphi\big]$$

**若 `F\le0$` 处处成立，结论立得。**但 `F$` 变号**：**

- **`\theta=\varphi=0$`：`F=0$`（`=c^2-1$`，`c=1$`）**
- **`\theta=0$`、`\varphi=\pi/3$`：`F=+0.25\gt0$` ✗**

**⇒ 「被积函数逐点非正」不成立。积分为负是因为权重 `e^{z(\cos\theta+\cos\varphi)}$` 集中在 `F\lt0$` 的区域，不是逐点。**

## E. 两端渐近可算（这是可写的部分）

**小 `z$`：`I_1\approx z/2$`、`I_2\approx z^2/8$`、`I_3\approx z^3/48$` ⇒**

$$I_1I_3-I_2^2\approx\frac{z^4}{96}-\frac{z^4}{64}=-\frac{z^4}{192}\lt0$$

**大 `z$`：`I_n(z)\approx\frac{e^z}{\sqrt{2\pi z}}\Big(1-\frac{4n^2-1}{8z}+\cdots\Big)$ ⇒**

$$I_1I_3-I_2^2\approx\frac{e^{2z}}{2\pi z}\Big[-\frac{(4-1)+(36-1)-(2\cdot15)}{8z}\Big]=\frac{e^{2z}}{2\pi z}\Big(-\frac{1}{z}\Big)\lt0$$

**⇒ 两端都严格负。**

## F. 状态

| 件 | 判 |
|---|---|
`I_1I_3\lt I_2^2$`（110 点，`z\in[10^{-4},200]$`） | `[M]$` |
逐项 `a_kc_k/b_k^2=(k+2)/(k+3)\lt1$` | **`[T]$`** |
逐项+AM-GM 路线 | **`[F]$`**（`z=20$` 时系数 `8$`–`24$`） |
积分被积函数逐点非正 | **`[F]$`**（`\theta=0,\varphi=\pi/3$`：`F=+0.25$`） |
两端渐近（小 `z$`、大 `z$`） | **`[T]$`** |
**`I_1I_3\lt I_2^2$` 的证明** | **`[U]$`**，需 Hankel 全正性型论证 |
`r'(0)\lt0$` | `[M]$` |

## G. 下一轮

**查 Hankel 行列式 `\det\begin{pmatrix}I_2&I_3\\I_1&I_2\end{pmatrix}$` 的已知结果：`I_n(z)$` 作为 `n$` 的序列是**全正序列**（`Toeplitz+Hankel`），而 `\det\gt0$` 正是二阶全正性。若文献有，直接引用并核；没有，用 `I_n=\frac1\pi\int e^{z\cos\theta}\cos n\theta$` 的 `\int\cos` 正交性做 Gram 矩阵论证。**

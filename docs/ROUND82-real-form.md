# 第 82 轮：`M_f(i\mu)$` 的**实形式** `\widetilde M$`（非对称，但恰好实）

---

## A. 辫子相似

**取 `D=\mathrm{diag}(i^a)$`，`\widetilde M:=D^{-1}M_f(i\mu)D$。**

**因 `a-b$` 与 `a+b$` 同奇偶，逐项化简：**

$$\widetilde M_{ab}=i^{\,b-a+|a-b|}J_{|a-b|}(\mu)-i^{\,2b+2}J_{a+b+2}(\mu)$$

- **`b\ge a$`：`b-a+|a-b|=2(b-a)$ 偶 ⇒ `i^{2(b-a)}=(-1)^{b-a}$`（实）；**
- **`b\lt a$`：`b-a+|a-b|=0$` ⇒ `i^0=1$`（实）；**
- **`2b+2$` 偶 ⇒ `i^{2b+2}=(-1)^{b+1}$`（实）。**

$$\boxed{\ \widetilde M_{ab}=(-1)^{\max(b-a,\,0)}J_{|a-b|}(\mu)-(-1)^{b+1}J_{a+b+2}(\mu)\quad\text{（全实）}\ }$$

**核：`\max|\mathrm{Im}\,\widetilde M|=0.000\times10^{0}$。**

**（`\mu=0.3$` 的 `5\times5$` 块：对角 `0.988792$`、`0.977605$`；`\widetilde M_{01}=-0.148878$`、`\widetilde M_{10}=+0.148878$` ⇒ 非对角**反对称**。）**

## B. 谱对得上

**`\mathrm{spec}\big(M_fKM_f\big)=\mathrm{spec}\big(\widetilde MK\widetilde M\big)$`，`\max$` 差 `6.7\times10^{-16}$。**

**⇒ `T(i\mu)$` 相似于一个**实矩阵**。**

## C. 但它不是对称的，且 `\widetilde M$` 的本征值复数

| 量 | 值 |
|---|---|
`\max\|\widetilde M-\widetilde M^{\mathsf T}\|$` | `0.2978$` |
`\max\|\mathrm{Im}\,\mathrm{eig}(\widetilde M)\|$` | `0.2857$` |
**`\mathrm{spec}(T(i\mu))$** | **实（`\mu\lt\mu^*$` 时，第 71 轮）** |

**⇒ `\widetilde M$` 是**实、非对称、有复本征值**，但 **`\widetilde MK\widetilde M$` 的谱是实的**。**

**⇒ 这不是"对称化"的路线，是"实矩阵 + 实谱"的路线。**

## D. 结构（从矩阵看出来的）

**`\widetilde M = \mathrm{sym} + \mathrm{anti}$`：**

- **非对角元反对称（`\widetilde M_{a,a+1}=-\widetilde M_{a+1,a}$`）；**
- **对角 `\widetilde M_{aa}=J_0-(-1)^{a+1}J_{2a+2}$`（随 `a$` 变号）。**

**⇒ `\widetilde M$` 是 Hankel + 反对称 Toeplitz 的组合形状。**

## E. 状态

| 件 | 判 |
|---|---|
**`\widetilde M$` 的实形式（`\max|\mathrm{Im}|=0$）** | **`[T]$`** |
**`\mathrm{spec}(MKM)=\mathrm{spec}(\widetilde MK\widetilde M)$（`6.7\times10^{-16}$）** | **`[T]$`** |
**`\widetilde M$` 非对称、有复本征值** | **`[T]$`**（数值） |
**`\widetilde MK\widetilde M$` 的谱实（`\mu\lt\mu^*$）** | `[M]$`（第 71 轮） |
**`\tilde r(-u)$` 的完全单调性** | **`[U]$`** |
**`\lambda\ge0.36$ 段** | **`[T]$`-级** |
**O1 / O2 / O3** | **`[F]` / `[T]` 解决 / `[U]$（单链不需要）** |
单链 `(R)$` | `[U]$` |
**Clay** | **`[U]$`** |

## F. 下一轮

**`\tilde r(-u)$` 的完全单调性。有了 `\widetilde M$` 的实形式，可试：**

$$\tilde r(-u)=\frac{\lambda_2(u)}{\lambda_1(u)},\qquad \lambda_i\ \text{是}\ \widetilde M(\mu)K\widetilde M(\mu)\ \text{的本征值},\ \mu=\sqrt u$$

**若 `\widetilde M(\mu)$ 作为 `\mu$` 的函数有显式单调结构（如 `\widetilde M = A(\mu)$` 且 `A'` 半正定），则 `\lambda_2/\lambda_1$` 的单调性可由本征值的单调性得。**

**⇒ 下一轮算 `\frac{d}{d\mu}\widetilde M$` 并查它的符号结构。**

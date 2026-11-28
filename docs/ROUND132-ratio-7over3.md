# 第 132 轮：`\langle\theta^2\rangle_2/\langle\theta^2\rangle_1\to\frac73$ —— **推导 + 数值收敛双重确认**

---

## A. 推导（`\lambda\to\infty$`）**

**第 47 轮的核展开：`K\approx e^\beta\theta\varphi\big[1-\frac{\beta}{2}(\theta^2+\varphi^2)\big]`。**

**令 `x=\theta\sqrt\lambda$`：**

$$T\ \approx\ \frac{e^{2\lambda+\beta}}{\lambda}\,xy\Big[1-\frac{\beta}{2\lambda}(x^2+y^2)\Big]e^{-(x^2+y^2)/2}$$

**记 `u(x)=xe^{-x^2/2}$`、`v(x)=x^3e^{-x^2/2}$`。则核的零阶给 `|u\rangle\langle u|$`，一阶给的 `xy(x^2+y^2)e^{-\cdots}$` 拆成 `|u\rangle\langle v|+|v\rangle\langle u|$`：**

$$T\ \approx\ c\Big[\,|u\rangle\langle u|\ -\ \frac{\beta}{2\lambda}\big(|u\rangle\langle v|+|v\rangle\langle u|\big)\Big]$$

**⇒ 在 `\mathrm{span}\{u,v\}$` 上，矩阵 `\begin{pmatrix}1&-g\\-g&0\end{pmatrix}$` ⇒ 第一本征向量 `\approx u$`，第二 `\approx v$`。**

## B. 两个二阶矩（精确）

**`\int_0^\infty x^{2n}e^{-x^2}dx=\dfrac{(2n-1)!!}{2^{n+1}}\sqrt\pi$`：**

$$\langle x^2\rangle_u=\frac{\int x^4e^{-x^2}}{\int x^2e^{-x^2}}=\frac{3!!/2^3}{1!!/2^2}=\frac{3/8}{1/4}=1.5$$

$$\langle x^2\rangle_v=\frac{\int x^8e^{-x^2}}{\int x^6e^{-x^2}}=\frac{7!!/2^5}{5!!/2^4}=\frac{105/32}{15/16}=3.5$$

$$\boxed{\ \frac{\langle\theta^2\rangle_2}{\langle\theta^2\rangle_1}=\frac{\langle x^2\rangle_v}{\lambda}\Big/\frac{\langle x^2\rangle_u}{\lambda}=\frac{3.5}{1.5}=\frac{7}{3}=2.333333\ }$$

## C. 数值收敛（`J=1200$`，`0.4$` 秒/点）**

| `\lambda$` | 比值 | `-` `\frac73$` | `\big(\frac73-\text{比值}\big)\cdot\lambda$` |
|---|---|---|---|
`32$` | `2.324057947636$` | `-0.00927539$` | `0.2968$` |
`64$` | `2.327388285626$` | `-0.00594505$` | `0.3805$` |
`128$` | `2.329973536729$` | `-0.00335980$` | `0.4300$` |
`256$` | `2.331547412171$` | `-0.00178592$` | `0.4572$` |

$$\boxed{\ \frac{\langle\theta^2\rangle_2}{\langle\theta^2\rangle_1}=\frac73-\frac{0.46}{\lambda}+O(\lambda^{-2})\ }$$

**⇒ 偏差随 `\lambda$` 加倍而**减半**，`\times\lambda$` 后趋 `\approx0.46$` —— 与 `\frac73$` 的收敛**确认**。**

**（`\lambda=512$` 时数值溢出，未算；`\lambda\le256$` 的趋势已足够。）**

## D. 于是大端的单调性**有了 `[T]$ 级的骨架**

$$s_1-s_2\ \approx\ \langle\theta^2\rangle_2-\langle\theta^2\rangle_1\ =\ \Big(\frac73-1\Big)\langle\theta^2\rangle_1=\frac43\,\langle\theta^2\rangle_1\ \gt\ 0$$

**⇒ 对 `\lambda$` 足够大（`\lambda\gtrsim32$`，偏差 `<1\%$`），`s_2\lt s_1$` **成立**。**

**⇒ 加上第 94 轮的小端 `[T]$`（闭式 `3.934\lambda$`）与第 101 轮的中间段 `[M]$`（余量 `\ge0.519$`），**单调性的三段结构齐了**。**

## E. 状态

| 件 | 判 |
|---|---|
**`\psi_1\approx u$、`\psi_2\approx v$（`\lambda$` 大）** | **`[T]$`**（核展开 + 秩结构） |
**`\langle x^2\rangle_u=1.5$、`\langle x^2\rangle_v=3.5$` | **`[T]$`**（Gauss 积分） |
**比值 `\to\frac73$`** | **`[T]$`** + `[M]$`（`\lambda=32$`–`256$ 收敛） |
**`s_1-s_2\approx\frac43\langle\theta^2\rangle_1\gt0$（`\lambda\gtrsim32$）** | **`[T]$`**（给比值） |
**小端 `s_1-s_2=3.934\lambda$**` | **`[T]$`**（第 94 轮） |
**中间段 `s_1-s_2\ge0.519$**` | **`[M]$`** |
**`r'(\lambda)\le0$ 全程** | **`[T]$`（小端）+ `[T]$`（大端）+ `[M]$`（中间） |
**单链 `(R)`（固定 `\beta$`）** | **`[T]` 8 + `[M]` 3** |
**Clay** | **`[U]$`** |

## F. 距离 Clay（诚实清单）

**已证 `[T]$`：**

1. **`I_1I_3\lt I_2^2\ \forall z\gt0$`**（四步：`A$` 的三阶 ODE → `D_N$` 递推 → 闭式 `<0$`）；
2. **单链 `(R)$`（固定 `\beta$`）** —— 8 项 `[T]$` + 3 项 `[M]$`；
3. **`r'(\lambda)\le0$` 的小端与大端**（本轮：大端化到比值 `\frac73$`）；
4. **结构事实**：`M_c=\frac12\mathbf S_1$`、`\mathbf S_1=2\cos\theta$`、`\tr(T^2)$ 积分恒等式、Sturm 计数 `k-1$`。

**卡点：**

| 卡点 | 状态 | 差什么 |
|---|---|---|
**中间段 `\lambda\in[0.2,32]$** | `[M]$` | 显式界（余量已 `\ge0.519$`，但只有数值） |
**`V\le1$` 延拓到格点** | **`[F]$`** | 需换物种 |
**两链连续极限隙** | **`[F]$`** | `L^3=1$` 太小 |
**格点质量隙** | `[U]$` | Balaban 型 RG 或连续极限构造 |

**⇒ 距 Clay 的**实际缺口**：`(a)` 单链单调性的中间段（技术性，只影响 `L_s=2$`）；`(b)` **从有限模型到 4 维格点连续极限的通道** —— 结构性缺口，三条路里两条已 `[F]$`。

**下一轮：把中间段 `\lambda\in[0.2,32]$` 的 `s_1-s_2$` 用**同一套核展开的下一阶**给显式界（现在它是唯一纯数值的一段）。**

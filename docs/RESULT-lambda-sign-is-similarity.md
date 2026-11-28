# 第 346 轮：`λ → −λ` 是符号相似变换 ⇒ TN₂ 是基依赖量，不是 `(R)` 的正确不变量

---

## A. 精确的对偶关系

**`j=a/2`、`j'=b/2`（`a,b` 整数）。望远镜式给出**

$$(M_f)_{j'j}=I_{|a-b|}(\lambda)-I_{a+b+2}(\lambda)$$

**且 `a+b\equiv|a-b|\pmod 2`。用 `I_n(-\lambda)=(-1)^nI_n(\lambda)`：**

$$\boxed{\ (M_f)_{\,j'j}(-\lambda)\ =\ (-1)^{|a-b|}\,(M_f)_{\,j'j}(\lambda)\ }$$

$$\Longleftrightarrow\qquad M_f(-\lambda)\ =\ D\,M_f(\lambda)\,D,\qquad D=\mathrm{diag}\big((-1)^{a}\big)_{a}$$

**（`D` 是行列式的对角符号矩阵，`D²=1`。）**

## B. 直接推论：谱不变

**`K` 在特征标基下是对角的（本征值 `t_j`），故 `DKD=K`。于是**

$$T(-\lambda)=M_f(-\lambda)KM_f(-\lambda)=D\,M_f(\lambda)KM_f(\lambda)\,D=D\,T(\lambda)\,D$$

$$\boxed{\ \mathrm{spec}\,T(-\lambda)=\mathrm{spec}\,T(\lambda)\qquad\Longrightarrow\qquad r(T)\ \text{与}\ \lambda\ \text{的符号无关}\ }$$

**⇒ 第 342 轮 `r(T)/(I_2/I_1)\in[0.14,0.85]` 的结论不受符号影响 ✓**

## C. 而 TN₂ 不是这样

**数值（1296 个 `2×2` 子式，`j\in\{0,\frac12,\dots,4\}`）：**

| `λ` | `−3.0` | `−1.5` | `−0.5` | **`+0.5`** | **`+1.5`** | **`+3.0`** |
|---|---|---|---|---|---|---|
负子式 | `640` | `640` | `626` | **`0`** | **`0`** | **`0`** |

**⇒ 谱相同（B），TN₂ 不同（C）⇒ **TN₂ 是基依赖量**。**

$$\boxed{\ \textbf{TN}_2\ \text{不是 }(R)\ \text{的正确不变量 ⇒ 第 344–345 轮把它当作证明目标是路线错误}\ }$$

## D. 这个发现的价值

**（1）它解释了为什么 `λ<0` 时 TN₂ 失败而 `r(T)` 仍 `≤I_2/I_1`：两者是不同层次的对象。**
**（2）它排除了「用 TN₂ + 变差缩减推 `(R)`」这条路线——不是因为 TN₂ 证不出来，而是因为即便证出来也不够。**
**（3）它把 `(R)` 的正确对象锁定为**谱**：`t_2(M_fKM_f)/t_1(M_fKM_f)` 与 `t_2(K)/t_1(K)` 的比较。**

## E. 仍然站得住的

| 件 | 判 |
|---|---|
`(M_f)_{j'j}=I_{|a-b|}(\lambda)-I_{a+b+2}(\lambda)` | `[T]`（望远镜；数值 `2.0e−16`） |
`(M_f)_{j'j}(-\lambda)=(-1)^{|a-b|}(M_f)_{j'j}(\lambda)` | `[T]` |
`M_f(-\lambda)=DM_f(\lambda)D`，`spec T(-\lambda)=spec T(\lambda)` | `[T]` |
`r(T)\le I_2/I_1`（单链，`\lambda>0`） | `[M]` 12/12，比值 `0.14–0.85` |
`(R)` 的证明 | `[U]`；**正确对象是谱，不是 TN₂** |

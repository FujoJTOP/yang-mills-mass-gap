# 第 348 轮：两条正性（`K` 核逐点正、`M_f^{-1}` 棋盘符号稳）

---

## A. 迁移矩阵核逐点为正（`[M]`，`β=2.5`，`KMAX=40`）

$$K(\theta,\varphi):=\sum_{k\ge0}t_k\,u_k(\theta)u_k(\varphi),\qquad t_k=\frac{2I_{2k+1}(\beta)}{\beta},\quad u_k=\sin((2k+1)\cdot)$$

**`401×401` 网格（`40401` 点）：**

$$\boxed{\ \min K=+6.13\times10^{-6}\ ,\quad \max K=+2.42\ ,\quad \text{负点数}=0\ }$$

**⇒ `K(\theta,\varphi)\ \gt\ 0$ 逐点成立。**

**⚠️ 它不是 `e^{\beta\cos(\theta-\varphi)}`（那个在 `(0.5,0.5)` 处是 `12.18`，而 `K` 是 `0.85`）。正确读法：由 `SU(2)` **加法定理**，`K` 是 `e^{\beta\cos\omega}` 对相对方位角平均后的量，平均不改变正性。**

## B. `M_f^{-1}` 的棋盘符号（`[M]`，18/18）

$$\boxed{\ D\,M_f^{-1}\,D\ \gt\ 0\ \text{逐元}\qquad\forall\lambda\in\{0.2,0.5,1,2,4,8\},\ J_2\in\{8,14,20\}\ }$$

**18 组全部 0 负元（最小 `+1.1\times10^{-25}$` … `+1.0$`）。**

**⇒ `M_f^{-1}` 的符号恒为 `(-1)^{b+j}` ⇒ `M_f` 属**振荡矩阵**（Gantmacher–Krein）类。**

## C. 由此得到的表示

**用 max–min 恒等式（第 345 轮）：`(M_f)_{xy}=\frac{2}{\pi}\int w(\theta)u_xu_y\,d\theta`，`w=e^{\lambda\cos\theta}`。**

$$T[j,j']=\Big(\frac{2}{\pi}\Big)^2\iint w(\theta)w(\varphi)\,u_j(\theta)\,K(\theta,\varphi)\,u_{j'}(\varphi)\,d\theta\,d\varphi$$

**即 `T=A_w^{\mathsf T}K_{\rm op}A_w`，`A_w[j,\theta]=w(\theta)u_j(\theta)`。**

**⇒ `T` 是**正核 `K_{\rm op}` 在函数系 `{w\,u_j}` 下的矩阵**。**

## D. 这条表示解释了 `(R)` 在 `λ=0` 的边缘情形

**`λ=0` 时 `w≡1`、`M_f=I`，`T=K` 精确，`r(T)=t_2/t_1=I_2/I_1` **取等**（与第 333 轮 `L_s=2` 取等的结论一致）。**

**`λ>0` 时数值上 `r(T)<I_2/I_1`（12/12，比值 `0.14–0.85`）。**

$$\boxed{\ (R)\iff\text{「乘上峰值在 }\theta=0\text{ 的正权 }w\text{ 会**压低**谱比」}\ }$$

**这一步仍是 `[U]`：`A_w` 的列不正交，`T` 不是 `K_{\rm op}` 的压缩，没有现成的交错定理直接给出。**

## E. 状态

| 件 | 判 |
|---|---|
`K(\theta,\varphi)>0` 逐点 | `[M]`（`40401` 点，0 负） |
`D M_f^{-1}D>0` | `[M]`（18/18） |
`T=A_w^TK_{\rm op}A_w` 表示 | `[T]`（由 max–min） |
`r(T)\le I_2/I_1` | `[M]` 12/12；证明 `[U]` |
| 缺的那一步 | `A_w` 列不正交，无交错定理可用 |

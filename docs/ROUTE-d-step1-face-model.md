# (d) 第一步：`SU(2)$` 面模型的显式权重写出来了，但**权重不定**

---

## A. 显式形式

**`e^{\beta q_P}=\sum_jt_j\chi_j(U_P)$，`q_P=\tfrac12\mathrm{tr}U_P=\chi_{1/2}$，**

$$t_j=\frac{2(2j+1)\,I_{2j+1}(\beta)}{\beta}\ \gt0$$

**要证（`\chi_j(I)=2j+1$` 归一，`\int\chi_j\chi_{j'}=\delta_{jj'}$`）：`\sum_jt_j(2j+1)=e^{\beta}$`。**

**每条 link 属 6 个 plaquette，逐 link 积分给一个 `6j$` 符号 ⇒**

$$\boxed{\ Z=\sum_{\{j_P\}}\ \prod_{P}t_{j_P}\ \prod_{\text{links}}\ \Big\{\begin{matrix}j_1&j_2&j_3\\j_4&j_5&j_6\end{matrix}\Big\}\ }$$

**这就是 `d=4$` `SU(2)$` 对偶面模型：面 = plaquette，自旋 = `j_P$`，顶点权重 = `6j$`。**

## B. 但 `6j$` 符号**半数取负**

**扫描 `j\in\{0,\frac12,1,\frac32,2,\frac52,3\}$` 的全部六元组（`7^6=117649$` 组）：**

| 量 | 值 |
|---|---|
非零 `6j$` | `3408$` |
**负值** | **`1670$`（`49.00\%$`）** |
零 | `11488$` |

**样例如下：**

$$\Big\{\begin{matrix}0&0&0\\ \frac12&\frac12&\frac12\end{matrix}\Big\}=-0.7071,\qquad \Big\{\begin{matrix}0&0&0\\\frac32&\frac32&\frac32\end{matrix}\Big\}=-0.5000,\qquad \Big\{\begin{matrix}0&1&1\\1&1&1\end{matrix}\Big\}=-0.3333$$

## C. 这直接阻止规格里的 Peierls 计划

**Peierls 下界要求「面权重 `w_\beta(\Sigma)\ge0$`」，于是「去掉一部分面 ⇒ `Z$` 只减不增」。**

$$\boxed{\ \text{`SU(2)$` 的 }6j\text{ 有 }49\%\text{ 取负 ⇒ 面权重不定 ⇒ 「子集求和单调」不成立 ⇒ Peierls 不能直接用}\ }$$

**⇒ `docs/ROUTES-ade-specifications.md` 里 (d) 的 D4 第 2 步（Peierls 型论证）的**前提不满足**。**

## D. 还剩下的两条（这是 (d) 的真实入口）

| 入口 | 说明 | 已知 |
|---|---|---|
**反射正性 (OS)** | 转移矩阵 `T\ge0$` ⇒ `\langle W(C)\rangle$` 的 RP 结构 ⇒ 弦张力 | `SU(2)$` 的 `T\ge0$` 已由第 312 轮的机器给出（`\chi_L\le1$` 那条链用到 `T_\beta\ge0$`） |
**Pauli / Biedenharn–Elliott 重排** | 把 `6j$` 的符号用正交关系消掉，写成「平方和」 | 未核 |

**⇒ (d) 的下一步应改成：用 RP 而不是 Peierls；即从 `T\ge0$` 出发，构造 Wilson 环的 RP 分解，把弦张力写成 `-\lim\frac{1}{\mathrm{Area}}\ln\langle W\rangle$` 并求下界。**

## E. 状态

| 件 | 判 |
|---|---|
`t_j=\frac{2(2j+1)I_{2j+1}}{\beta}\gt0$` | `[T]$`（`I_\nu$` 正） |
面模型显式形式（面/自旋/`6j$`） | `[T]$`（标准特征标展开） |
`6j$` 半数取负（`49.00\%$`，`3408$` 个非零样本） | **`[T]$`**（`sympy.wigner_6j` 精确） |
**Peierls 直接用** | **`[F]$`** |
RP 路线 | 未开工 |

## F. 下一轮

**(d) 第 2 步改走 RP：从第 312 轮已证的 `T_\beta\ge0$` 出发，写 Wilson 环的 RP 分解，求弦张力下界。**

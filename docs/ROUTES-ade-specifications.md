# (a) / (d) / (e) 三条路线的完整写出

> 用户指示：先把 (a)(d)(e) 三条全部写出来，(b) 已派出，(a)(d)(e) 一个一个做，计数器重置为 0/60。
> 本文件是三条路线的**完整规格**：精确命题、它能给什么、已知记录、攻击步骤、判定标准。每条都可以独立开工。

---

# 路线 (a)：非 Abel 高斯支配 / 横截 Hessian 控制

## A1 要证的精确命题

**记 `\Lambda_L$` 为 `L^4$` 环面格点，`U=(U_e)$` 为链变量，`d\mu_\beta$` 为 Wilson 测度（`S=\beta\sum_P(1-q_P)$`，`q_P=\frac12\operatorname{tr}U_P$）。**

**「高斯支配」（GD）的标准形式：对任意 `f=(f_\gamma)_{\gamma\in\mathcal P^{or}}$（有向 plaquette 上的 `\mathfrak{su}(2)$`-值函数），**

$$\boxed{\ \Big\langle \exp\Big(\sum_\gamma f_\gamma(U)\Big)\Big\rangle_{\Lambda_L}\ \le\ \exp\Big(\frac{1}{4\beta}\sum_\gamma\|f_\gamma\|^2\Big)\ }$$

**其中 `\|f_\gamma\|^2=\sup_x\sum_{a=1}^{3}\|f_\gamma(x)\|_{L^2(dU)}^2$`-型。**（精确范数取法与文献一致即可，关键是 `1/\beta$` 系数与「平方和」右端。）

**「横截 Hessian 控制」是它的充分形式：**

$$\boxed{\ \lambda_{\min}\big(\mathscr H\big|_{\text{横截面 transverse}}\big)\ \ge\ 0\ }$$

**其中 `\mathscr H=\nabla^*\nabla+\mathrm{Ric}+\mathrm{Hess}\,S$`，横截面 = 正交于规范轨道。**

## A2 它能给什么

**(GD) ⇒ Fröhlich–Simon–Spencer 红外界 ⇒ 对 `|\mathrm{Re}\,k|\lt m$`，**

$$\big|\hat C(k)\big|\ \le\ \frac{\text{const}}{\beta}\cdot\frac{1}{m^2-|\mathrm{Re}\,k|^2}\qquad\Longrightarrow\qquad m\ \ge\ m_0(\beta)$$

**⇒ 逐点衰减（正是 `§1.07$` 要的形状）。**

## A3 已知记录（本轮次实测）

**作用量 Hessian 在 Savvidy 构型（`q_P=0$`、`|a_P|=1$`）的单 plaquette 谱精确为**

$$\{-(1+\sqrt2)\beta\times2,\ -(2-\sqrt2)\beta\times2,\ 0\times4,\ +(2-\sqrt2)\beta\times2,\ +(1+\sqrt2)\beta\times2\}$$

$$\boxed{\ \|\mathrm{Hess}\,S\|\ =\ 6(1+\sqrt2)\,\beta\ =\ 14.4853\,\beta\quad(\text{每边属 6 个 plaquette，Savvidy 点取等})\ }$$

**⇒ `\mathscr H\succeq(2-14.4853\beta)I$`，在 `\beta\gt0.138$` 不定。**
**且**规范方向精确落在核里**（数值 `0.000\times10^0`），负模是**物理横向模**。**
**文献：`(GD)$` 对 `U(1)$`、Villain、`O(n)$` 自旋模型有证明。**

## A4 攻击步骤（按此顺序）

1. **把 `\mathscr H\succeq(2-14.4853\beta)I$` 替换成「横截面上的」版本**：扣掉 Savvidy 的 4 个负方向后，剩下 5 个方向是否正？给出精确的横截面谱。
2. **若横截面仍不定**，改问：负模方向上的**测度权重**是否指数小（Savvidy 构型作用量 `=\beta$` 每 plaquette）？⇒ 走「条件 (GD)」。
3. **条件 (GD) 的坏集估计**：第 336 轮已做并集界 `4L^4e^{-6\beta}$`，在 `L\gtrsim35$` 失效。需换**非并集**（块/RG）版本。
4. **判定**：得到 `\beta`-一致的 `m_0\gt0$` 即成功；但必须过渐近自由酸测（`m_0(\beta)\to0$`）。

## A5 判定标准

| 结果 | 判 |
|---|---|
横截面谱精确给出且 `\ge0$` | `[T]$`（新的结构事实） |
横截面谱仍不定 | `(a)$` 关闭，转条件版 |
条件版在 `L$`-一致意义下立 | **大型 `[T]$`** |

---

# 路线 (d)：`d=4$` `SU(2)$` 对偶面模型的表面张力下界

## D1 要证的精确命题

**`SU(2)$` `d=4$` 规范理论有对偶表示：配分函数写成对「面」（世界面）的求和**

$$Z=\sum_{\partial\Sigma=\emptyset}w_\beta(\Sigma),\qquad w_\beta(\Sigma)\ \ge\ 0$$

**「弦张力」`\tau(\beta)$` 由 Wilson 环的周长律定义：**

$$\big\langle W(C)\big\rangle\ \sim\ \exp\big(-\tau(\beta)\,\mathrm{Area}(C)\big)\qquad(|C|\to\infty)$$

$$\boxed{\ \text{要证：}\ \exists\,\tau_0(\beta)\gt0\ \text{ 使 } \tau(\beta)\ge\tau_0(\beta)\ \text{ 对一切 }\beta\ }$$

**（或至少：`\tau(\beta)\gt0$` 对一切 `\beta\gt0$`。）**

## D2 它能给什么

**「质量 = 表面张力」：`m\ge\tau$`。⇒ 逐点衰减，`§1.07$` 要的形状。**

## D3 已知记录

**`d=4$` 非 Abel 情形该命题未见证明。`d=3$` `Z_N$` 有证明（Göpfert–Mack 型）。**
**`SU(2)$` 的面模型权重 `w_\beta$` 来自特征标展开 `e^{\beta q}=\sum_jt_j\chi_j$`，`t_j=2I_{2j+1}(\beta)/\beta\gt0$`（本项目已核）。**

## D4 攻击步骤

1. **写出 `SU(2)$` 面模型的显式权重**（用 `t_j$` 与 `SU(2)$` 的 `6j$`/`9j$` 符号；`SU(2)$` 的 `3j$`/`6j$` 组合有闭合公式）。
2. **找「面张力」的下界形式**：对「边界固定为 `C$`」的面求和，用**面积-周长不等式**或 **Peierls 型论证**给出 `\ge e^{-\tau_0\mathrm{Area}}$`。
3. **关键难点**：`d=4$` 面可以「鼓包」而不增大边界 ⇒ 需要**等周**控制 `\mathrm{Area}\le c\,\mathrm{diam}^2$`-型。
4. **判定**：拿到与 `\beta$` 无关（或 `\beta$`-显式）的 `\tau_0\gt0$` 即成功。

## D5 判定标准

| 结果 | 判 |
|---|---|
`SU(2)$` 面模型权重的显式形式 | `[T]$`（表示论） |
`\tau(\beta)\gt0$` 对一切 `\beta$` | **大型 `[T]$`** |
只在 `\beta$` 大或小成立 | 落入已知窗口，不算 |

---

# 路线 (e)：第七种机制

## E1 要证的精确命题

**构造一条新机制 `I$`，满足：**

**(E-1) 形状（逐点）：存在可观测量 `O$` 与显式 `F$`，使**

$$\big\langle O(x)O(0)\big\rangle_c\ \le\ F(|x|),\qquad F(n)\le F(0)\,e^{-m_0\,n}$$

**且 `m_0$` 由 `I$` 的输入显式给出。**

**(E-2) 四项排除（不得预设）**

| 不得预设 | 已核理由 |
|---|---|
`\mathrm{Hess}(-\log d\mu)\ge\kappa I$` | `\|\mathrm{Hess}\,S\|=6(1+\sqrt2)\beta$`（路线 (a) 的 A3） |
`\langle\phi_2,V\phi_2\rangle\le\langle\phi_1,V\phi_1\rangle$` | 两链已判 `[F]$`（`r=0.5891,\ 0.6743\gt I_2/I_1=0.5072$`） |
小参数 `\beta\ll1$` | 簇展开只到 `\beta\lesssim0.419$` |
`(GD)$` 型 | 与 (E-2) 第一条同 |

**(E-3) 酸测（`CLAUDE.md §1.1`）**

- **`m_0(\beta)\to0$` 当 `\beta\to\infty$`（渐近自由）。** 若给出 `\beta$`-一致的 `m_0\gt0$` ⇒ 与渐近自由矛盾 ⇒ 推倒。
- **自由场上必须失效**（否则 vacuous）。

## E2 已排除的八条承载结构

| # | 结构 | 排除理由 |
|---|---|---|
1 | 求和规则 `\chi` + 逐点非负 | 只给 `\le\chi$`，无衰减 |
2 | FKG / 正关联 | 给**下界**，方向相反 |
3 | Bakry–Émery | 同 (E-2) 第一条 |
4 | LSI / 输运不等式 | 需 `(GD)$` 型输入 |
5 | RP + 转移矩阵 | 折叠到两链比较（已 `[F]$`） |
6 | 特征标正定性 | 只给正性 |
7 | 离散极值原理 | 等价于作用量-作用量关联 ⇒ `\chi$`-类 |
8 | 求和规则 + `\beta$`-单调 | **实测 `\chi$` 非单调**（`L=4$`、6 点） |

## E3 第 9 条候选（尚未检查）——本轮次留下的唯一入口

**本项目唯一证明过的**逐点**陈述是 `t_1$`-增长（`T1$`），其证明方式值得注意：**

$$t_1(T)\ \ge\ T_{00}\quad\text{（变分，取**显式测试向量** }u_0=\sin\theta\text{）}$$

**⇒ 下界来自「显式测试对象」，不是来自泛函不等式。**

**⇒ (e) 的第 9 条攻击方向：对逐点量构造显式测试对象，把下界算出来。**

**具体化**：要逐点上界 `\langle q(x)q(0)\rangle_c$`，等价于要 `\chi$`-测度的「低谱质量」控制。已知 `\nu_0([0,s])\le s/6$`（由 `\chi\le1/6$`）。**这条已经是逐点型（在 `s$` 上逐点）。**

**⇒ 若能把 `\nu_0$` 与 `\langle q(x)q(0)\rangle_c$` 的**空间衰减**联起来（例如通过某个显式测试函数族），就得到 (e)。这是唯一没试过的入口。**

## E4 攻击步骤

1. **把 `\nu_0([0,s])\le s/6$` 与 `\chi=\sum_x\langle q(x)q(0)\rangle_c$` 的关系写出来**（已有）。
2. **构造一族显式测试函数 `\psi_k$`，使 `\langle\psi_k,T\psi_k\rangle$` 给出 `\hat C(k)$` 的下界** ⇒ 得到 `\lambda_{?}$` 的信息。
3. **检查这条是否又落回 `\chi$`-类**（很可能是）。
4. **判定**：若 (E-1)～(E-3) 全过 ⇒ **大型 `[T]$`**。

## E5 判定标准

| 结果 | 判 |
|---|---|
显式测试对象给出 `\hat C(k)$` 的**上**界且带衰减 | **大型 `[T]$**` |
只给出 `\chi$`-类 | 落回 (b) 的死路 |
无候选 | (e) 关闭 |

---

## 三条路线的优先级（按「可开工性」排序，不含难度评价）

| 路线 | 第一步是可执行的有限计算吗 | 第一步 |
|---|---|---|
**(a)** | 是 | 算 Savvidy 点 Hessian 的**横截面谱**（扣掉 4 个负方向后剩下 5 个方向），看是否 `\ge0$` |
**(d)** | 是 | 写出 `SU(2)$` 面模型权重的显式形式（用 `t_j$` 与 `6j$`/`9j$` 符号） |
**(e)** | 是 | 用 `\nu_0([0,s])\le s/6$` 构造显式测试函数族，看能否给出 `\hat C(k)$` 的上界 |

**按用户指示：(a) → (d) → (e)，一个一个做。**

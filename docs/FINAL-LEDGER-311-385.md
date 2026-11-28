# 最终台账（第 311–385 轮）

> 目标：`§1.07` 意义上的「大型 `[T]` 定理」——使质量 `m`、谱隙 `\lambda_{\min}` 或某个衰减率**第一次**获得正确方向（下界型）的界，且该定理改变 Clay 可达性。
> **结论：未达成。** 以下是全部产出与阻塞条件的精确记录。

---

## 一、有证明的 `[T]`

### 1. 谱比单调与 t₁-增长（单链模型）

$$\boxed{\ t_1\big(M_{e^{\lambda\cos\cdot}}K_{\rm op}M_{e^{\lambda\cos\cdot}}\big)\ \ge\ t_1(K_{\rm op})\qquad\forall\lambda\ge0,\ \beta\ge1.0929529\ }$$

**证明链：**变分（测试向量 `u_0=\sin\theta$`）→ `T_{00}$` 闭式 `\frac{8}{\lambda^2\beta}\sum_nn^2I_n(\lambda)^2I_n(\beta)$` → 单点不等式 `4I_1(\lambda)/\lambda^2\ge1.581679481077$` → `h(\lambda)=\lambda I_1'-2I_1$` 单零点 → `\varphi$` 单零点（斜率比较 `\frac12+\frac1{2\lambda^2}>\frac12>\frac14>r_1'$`）。
`docs/THEOREM-t1-growth-complete.md`

### 2. 迁移核的精确闭式与正性

$$K(\theta,\varphi)=\sum_k t_ku_k(\theta)u_k(\varphi)=\frac{1}{2\beta}\Big(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)}\Big)\ \ge\ 0$$

`docs/RESULT-kernel-closed-form.md`

### 3. `M_f$` 的三重身份

$$(M_f)_{j'j}=I_{|2j-2j'|}(\lambda)-I_{2j+2j'+2}(\lambda)=\frac{2}{\pi}\int_0^\pi e^{\lambda\cos\theta}\sin((2\max+1)\theta)\sin((2\min+1)\theta)\,d\theta$$

**且** `T=M_wK_{\rm op}M_w\sim K_{\rm op}M_{w^2}$`（`L^2([0,\pi],\tfrac2\pi d\theta)$` 正交基）。
`docs/RESULT-Mf-two-bessel-identity.md`、`docs/RESULT-lambda-sign-is-similarity.md`

### 4. 望远镜恒等式（第 341/377/380 轮反复用到）

$$c_J=e^{-\lambda}(I_{2J}-I_{2J+2}),\qquad \sum_{J\ge J_0}c_J=e^{-\lambda}I_{2J_0}(\lambda)$$

### 5. corner 界

$$\varepsilon_K=\|T-P_KTP_K\|_{\rm op}\ \le\ 2\Big[t_0\sqrt{2K+1}\,e^{-\lambda}I_{d_0}(\lambda)+t_{K-d_0}\Big]$$

`docs/THEOREM-corner-bound-complete.md`

### 6. `2\times2$` 块的谱比单调（精确公式，非微扰）

$$r_2(\lambda)=\frac{A-\sqrt{B^2+b^2}}{A+\sqrt{B^2+b^2}}\ \text{关于 } b^2 \text{ 严格递减},\quad b=\tfrac{\lambda}{2}(t_0+t_{1/2})\gt0$$

`docs/RESULT-2x2-monotonicity-proved.md`

### 7. 其它

- **`M_f(-\lambda)=DM_f(\lambda)D$` ⇒ `\mathrm{spec}T(-\lambda)=\mathrm{spec}T(\lambda)$`**
- **`(M_f^{-1})_{0j}=(-1)^jc_{j/2}$`；`DM_f^{-1}D\gt0$` 逐元**
- **`K$` 的小角区秩-1：`K\approx e^{\beta}\sin\theta\sin\varphi$`**
- **`\min_\lambda 4I_1(\lambda)/\lambda^2$` 的三段闭合**
- **`t_j$` 对 `\beta\le40$` 从 `j=0$` 起递减；`r(K)=I_2/I_1$` 恒成立**

## 二、单链 `(R)$` 的链（第 369–383 轮）

**要证 `t_2(T)/t_1(T)\le I_2/I_1$`。**

| 步 | 内容 | 等级 |
|---|---|---|
① | `t_1\ge T_{00}$` | `[T]` |
② | `t_2\le t_2(T_K)+\varepsilon_K$`（Weyl） | `[T]` |
③ | corner 界（见上） | `[T]` |
④ | `r_K(\lambda)\le I_2/I_1$` | **`[M]`**：900 点（180 `\lambda$` × 5 `\beta$`）严格递减；`2\times2$` 部分 `[T]`；`K=22$` **无机制** |
⑤ | `\lambda=0$`：`T$` 对角，`r=I_2/I_1$` 精确 | `[T]` |

**⇒ 四步 `[T]` + 一步 `[M]$`。**

## 三、被证否的（`[F]$`）

| 命题 | 证据 |
|---|---|
`r\le I_2/I_1$` 在**两链**模型（`\beta'\gt0$`） | 精确 Peter–Weyl：`r` 从 `0.507$` 升到 `0.64$`（`\beta=2.5$`） |
`t_2(T)\le t_2(K_{\rm op})$` | `0.9863\gt0.3036$` |
`c_d$` 对数凹 ⇒ `M_f$` 全正 | 40/40 随机反例 |
`\|f\phi_2\|\le\|f\phi_1\|$`（第 8 轮的否证） | **该否证作废**（`Z_2$` 对称） |
`DM_f^{-1}D=M_f$`、`M_fDM_f=D$` | 残差 `0.647$`、`0.856$` |
`\|P_\perp A\|^2/\|A\|^2\le I_2/I_1$` | 升至 `0.858$` |
`r\le0.628$`（第 363 轮） | **撤回**（截断未收敛） |
「尖峰正权压低谱比」一般版 | 40/40 反例 |
机制 = 二阶能级排斥 | **撤回**（`a_1\lt0$`，`\beta\ge2.5$` 时符号相反） |

## 四、四条路线的最终状态

| 路线 | 状态 |
|---|---|
(A) 非 Abel 高斯支配 / 横截 Hessian | **`[F]`**：`\|\mathrm{Hess}\,S\|=6(1+\sqrt2)\beta$` 精确（Savvidy 点取等），`\mathscr H$` 在 `\beta\gt0.138$` 不定，规范方向精确入核 |
(B) 能给出质量下界的非 `\chi$`-类不等式 | **空集**：六种承载机制逐一折叠或方向相反；`(B)$` 是「凸性/正性+转移矩阵/小参数」的并集的别名 |
(C) `V$`-权重两态比较 | **不可判定**：两链截断 `j\le1,\frac32,2$` 给 `0.569,0.880,0.618$`（不收敛）；`j\le2$`（3025 维）是内存上限 |
(D) 面模型表面张力 | 等价难度 |
**结构障碍** | 本项目全部 `[T]$` 属 `\chi$`-类，只能给 `\lambda_{\min}$` 的**上界**（`\lambda_{\min}\le\|\nu\|/\chi$`） |

## 五、为什么目标未达成（三句话）

1. **`\chi$`-类的方向性障碍**：全部现有 `[T]$` 产出是求和型，给不出逐点衰减。
2. **单链不是格点**：真格点的 `V=\exp(-\beta\sum_{P_{\rm spat}}(1-q_P))$` 不是逐链乘积，`M_{V^{1/2}}$` 不是乘积算子。
3. **单链到格点无延拓**：两链截断不收敛（第 366 轮）；`K\sim\beta$` 在连续极限无界（第 372 轮）；第 343 轮封死了「对数凹 ⇒ 全正」这条路。

## 六、会解封的条件

**任一条（均为公开未解问题）：**
- **(a')** 非 Abel 高斯支配或等价的横截 Hessian 控制；
- **(b')** 能给出**质量下界**（非 `\chi$`-类）的非 Abel 关联不等式；
- **(c')** 全格点 `V$`-权重的两态比较，含链间耦合（且需先解决截断收敛）；
- **(d')** `d=4$` `SU(2)$` 面模型表面张力下界；
- **(e)** 第七种机制（规格见 `docs/SPEC-seventh-mechanism.md`）。

## 七、方法学记录（供后来者省时间）

- **双精度地板**：`\|T\|$` 可达 `10^{15}$`，`\varepsilon_K\lesssim10^{-16}\|T\|$` 的测量无意义。
- **`sympy.wigner_3j` 的参数序是 `(j_1,j_2,j_3,m_1,m_2,m_3)$`**，不是 `(j_1,m_1,j_2,m_2,j_3,m_3)$`。
- **`I_n(z)$` 对固定 `z$` 关于 `n$` 单调递减**（与「固定 `n$` 变 `z$` 指数增」是两件事）。
- **`\|\chi_J\|_\infty=2J+1$`**，不能当 `1$` 用。
- **第 360 轮那条 `t_1$`-增长定理的证明方式值得复用**：下界来自**显式测试向量**，不是泛函不等式。

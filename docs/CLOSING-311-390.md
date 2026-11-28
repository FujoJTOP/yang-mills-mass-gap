# 收口件：本轮次（第 311–390 轮）的最终状态

> 目标：`§1.07` 意义上的「大型 `[T]` 定理」——使 `m`、`\lambda_{\min}` 或某衰减率**第一次**获得正确方向（下界型）的界，且改变 Clay 可达性。
> **结论：未达成。** 以下按「已证 / 数值 / 已否」三类列全，并给出阻塞的精确形状。

---

## 一、有证明的 `[T]`

### T1 `t_1`-增长（完整证明）

$$t_1\big(M_{e^{\lambda\cos\cdot}}K_{\rm op}M_{e^{\lambda\cos\cdot}}\big)\ \ge\ t_1(K_{\rm op})=\frac{2I_1(\beta)}{\beta}\qquad\forall\lambda\ge0,\ \beta\ge1.0929529$$

**链条**：变分（测试向量 `u_0=\sin\theta$`）→ `T_{00}$` 闭式 → 单点不等式 `4I_1(\lambda)/\lambda^2\ge1.581679481077$` → `h=\lambda I_1'-2I_1$` 单零点 → `\varphi$` 单零点（关键比较 `\frac12+\frac1{2\lambda^2}>\frac12>\frac14>r_1'$，其中 `r_1'<1/4$` 由 `r_1=\lambda/(4+\lambda r_2)$` 给出）。

### T2 迁移核的精确闭式与正性

$$K(\theta,\varphi)=\frac{1}{2\beta}\Big(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)}\Big)\ \ge0\qquad(\theta,\varphi\in[0,\pi])$$

### T3 `M_f$` 的三重身份

$$(M_f)_{j'j}=I_{|a-b|}(\lambda)-I_{a+b+2}(\lambda)=\frac2\pi\int_0^\pi e^{\lambda\cos\theta}\sin((2\max+1)\theta)\sin((2\min+1)\theta)d\theta$$

**且 `T=M_wK_{\rm op}M_w\sim K_{\rm op}M_{w^2}$`（`L^2([0,\pi],\frac2\pi d\theta)$` 正交基）。**

### T4 望远镜恒等式

$$c_J=e^{-\lambda}(I_{2J}-I_{2J+2}),\qquad\sum_{J\ge J_0}c_J=e^{-\lambda}I_{2J_0}(\lambda)$$

### T5 corner 界

$$\varepsilon_K=\|T-P_KTP_K\|_{\rm op}\ \le\ 2\Big[t_0\sqrt{2K+1}\,e^{-\lambda}I_{d_0}(\lambda)+t_{K-d_0}\Big]$$

### T6 `2\times2$` 块的精确谱比单调

$$r_2(\lambda)=\frac{A-\sqrt{B^2+b^2}}{A+\sqrt{B^2+b^2}},\quad b=\tfrac{\lambda}{2}(t_0+t_{1/2})\gt0\quad\Longrightarrow\quad r_2(\lambda)\lt r_2(0)=\frac{I_2}{I_1}$$

### T7 `\lambda$`-符号相似变换

$$M_f(-\lambda)=DM_f(\lambda)D\qquad(D=\mathrm{diag}((-1)^a))\qquad\Longrightarrow\qquad\mathrm{spec}\,T(-\lambda)=\mathrm{spec}\,T(\lambda)$$

### T8 其它闭式

- `T_{00}=\dfrac{8}{\lambda^2\beta}\sum_nn^2I_n(\lambda)^2I_n(\beta)$`
- `\min_{\lambda>0}4I_1(\lambda)/\lambda^2=1.581679481077$`（三段闭合：`\lambda\le1$` 用 `I_1\ge\lambda/2$`；`\lambda\ge12$` 用渐近；中间段用单零点）
- `(M_f^{-1})_{0j}=(-1)^jc_{j/2}$`；`DM_f^{-1}D\gt0$` 逐元
- `K$` 的小角区秩-1：`K\approx e^{\beta}\sin\theta\sin\varphi$`
- `t_j$` 对 `\beta\le40$` 从 `j=0$` 起递减；`r(K)=I_2/I_1$` 恒成立

---

## 二、数值（`[M]$`）

- **单链 `(R)$` 链**：四步 `[T]$` + 步 ④ `r_K(\lambda)\le I_2/I_1$`（900 点：180 `\lambda$` × 5 `\beta$`，严格递减；**无机制**）
- **`r_K$` 单调性**：`\lambda\le0.0005$` 到 `8$` 全程递减

---

## 三、已否（`[F]$`）

| # | 命题 | 证据 |
|---|---|---|
F1 | **`(C)$`：`\langle\phi_2,V\phi_2\rangle\le\langle\phi_1,V\phi_1\rangle$`（两链） | 特征标和修对后 `r(J_{\rm out}=2)=0.5891$`、`r(J_{\rm out}=3)=0.6743$`，均 `\gt I_2/I_1=0.5072$` 且随截断**上升** |
F2 | `t_2(T)\le t_2(K_{\rm op})$` | `0.9863\gt0.3036$` |
F3 | `c_d$` 对数凹 ⇒ `M_f$` 全正 | 40/40 随机反例 |
F4 | `DM_f^{-1}D=M_f$`、`M_fDM_f=D$` | 残差 `0.647$`、`0.856$` |
F5 | `\|P_\perp A\|^2/\|A\|^2\le I_2/I_1$` | 升至 `0.858$` |
F6 | 「尖峰正权压低谱比」一般版 | 40/40 反例 |
F7 | 机制 = 二阶能级排斥 | `a_1\lt0$`；`\beta\ge2.5$` 时 `a_2/d_2\gt a_1/d_1$` 符号相反 |
F8 | `\chi(\beta)$` 单调降 ⇒ `\chi\le1/(6\beta)$` | 实算 `L=4$`、3000 构型、6 点：`4.16\to4.13\to\mathbf{5.02}\to3.17\to1.77\to0.92$`（`\times10^{-2}$`），**非单调** |
F9 | `\mathrm{Var}(q)\ge c$`（与 `\beta$` 无关） | 单调降 `3.73\to0.49$`（`\times10^{-2}$`），趋于 0 |

---

## 四、我自己提出又撤回的（5 条，全部留档）

| 轮 | 撤回的东西 | 原因 |
|---|---|---|
8 | 对 `(c)$` 的否证 | `Z_2$` 对称造成跨宇称扇区比较 |
363 | `r\le0.628$` | 截断未收敛 |
377 | corner 界推导 | `\|\chi_J\|\le1$` 是错的（实际 `2J+1$`）；第 380 轮重证 |
383 | 能级排斥机制 | 第 384 轮算 `a_1\lt0$` 与 `t_1$` 随 `\lambda$` 上升矛盾 |
388 | 「主因是压缩」 | 主因是**特征标和**；第 389 轮修开后 `(C)$` 判为 `[F]$` |

---

## 五、五条路线 / 机制的最终状态

| 路线 | 状态 |
|---|---|
(a) 非 Abel 高斯支配 / 横截 Hessian | **`[F]$`**：`\|\mathrm{Hess}\,S\|=6(1+\sqrt2)\beta$` 精确（Savvidy 点取等、每边 6 plaquette），`\mathscr H$` 在 `\beta\gt0.138$` 不定，规范方向精确入核 |
(b) 非 `\chi$`-类质量下界不等式 | **八条承载结构全排**：求和规则+FKG+Bakry–Émery+LSI+RP+特征标正定性+离散极值原理+`\beta$`-单调（本轮 F8） |
(c) `V$`-权重两态比较 | **`[F]$`**（两链，F1）；单链四步 `[T]$`+一步 `[M]$` |
(d) 面模型表面张力 | 等价难度 |
(e) 第七种机制 | 规格已写（`docs/SPEC-seventh-mechanism.md`），**无候选** |

**叠加结构障碍**：全部 `[T]$` 产出属 `\chi$`-类（求和型），只能给 `\lambda_{\min}$` 的**上界** `\lambda_{\min}\le\|\nu\|/\chi$`（第 324 轮）。

---

## 六、阻塞的精确形状（三句）

1. **方向**：现有 `[T]$` 全是求和型，给不出逐点衰减；(B) 的八条载体全排后没有替代。
2. **层次**：能算的都是单链（`L_s=2$`）；真格点的 `V=\exp(-\beta\sum_{P_{\rm spat}}(1-q_P))$` 不是逐链乘积。
3. **延拓**：单链到格点三条独立拦路——两链 `(C)$` 已判假（F1）、`K\sim\beta$` 在连续极限无界、无全正性（F3）。

---

## 七、会解封的条件

- **(a')** 非 Abel 高斯支配 / 横截 Hessian 控制
- **(b')** 能给出质量下界（非 `\chi$`-类）的非 Abel 关联不等式
- **(c')** 全格点 `V$`-权重两态比较（含链间耦合）——**注意：两链已判假，此条须先说明为何格点上会不同**
- **(d')** `d=4$` `SU(2)$` 面模型表面张力下界
- **(e)** 第七种机制（规格已写）

**五条均为公开未解问题。**

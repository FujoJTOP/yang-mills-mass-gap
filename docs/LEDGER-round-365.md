# 第 365 轮：完整台账与资源边界

---

## A. 数值路线的资源边界（本轮核清）

**两链截断 Peter–Weyl 的规模：每链维数 `d=\sum_{2j\le J}(2j+1)`，两链 `d^2$`。`M_f$` 由 `\sum_{J'}(2J'+1)^2$` 个 `d^2\times d^2$` 的 Kronecker 积累加而成。**

| `j_{\max}$` | 每链 `d$` | 两链 `d^2$` | 单个 Kronecker 积内存 | 块数 | 是否可行 |
|---|---|---|---|---|---|
`1$` | `14$` | `196$` | `0.3$` MB | `9$` | 可行 |
`\frac32$` | `30$` | `900$` | `6.5$` MB | `16$` | 可行（本轮已跑） |
`2$` | `55$` | `3025$` | `72$` MB | `25$` | **累积约 1.8 GB，勉强** |
`\frac52$` | `91$` | `8281$` | `540$` MB | `36$` | **不可行** |

$$\boxed{\ \text{收敛检验的最高可及截断是 }j\le2\ (3025\times3025)\text{，且需按块累积而非一次性 Kronecker}\ }$$

**第 364 轮已证明 `j\le1$` 与 `j\le\frac32$` 之间差 `0.31$` ⇒ 未收敛。`j\le2$` 能否收敛**未知**；`j\le2$` 之上不可及。**

## B. 35 轮累计台账

### B1. 完整的 `[T]` 定理（有证明、外部输入已标）

| 定理 | 位置 |
|---|---|
**`t_1(M_{e^{\lambda\cos\cdot}}K_{\rm op}M_{e^{\lambda\cos\cdot}})\ge t_1(K_{\rm op})$**，`\forall\lambda\ge0$`、`\beta\ge1.0929529$` | `THEOREM-t1-growth-complete.md` |
**`K(\theta,\varphi)=\frac{1}{2\beta}(e^{\beta\cos(\theta-\varphi)}-e^{\beta\cos(\theta+\varphi)})\ge0$** | `RESULT-kernel-closed-form.md` |
**`(M_f)_{j'j}=I_{|2j-2j'|}(\lambda)-I_{2j+2j'+2}(\lambda)$**（望远镜） | `RESULT-Mf-two-bessel-identity.md` |
**`(M_f)_{j'j}=\frac{2}{\pi}\int e^{\lambda\cos\theta}\sin((2\max+1)\theta)\sin((2\min+1)\theta)d\theta$** | 同上 |
**`M_f(-\lambda)=DM_f(\lambda)D$` ⇒ `\mathrm{spec}T(-\lambda)=\mathrm{spec}T(\lambda)$** | `RESULT-lambda-sign-is-similarity.md` |
**`T_{00}=\frac{8}{\lambda^2\beta}\sum_nn^2I_n(\lambda)^2I_n(\beta)$** | `RESULT-T00-closed-form-bessel-reduction.md` |
**`\min_\lambda4I_1(\lambda)/\lambda^2=1.581679481077$`**（三段闭合） | `RESULT-single-point-inequality-3-regions.md` |
**`(M_f^{-1})_{0j}=(-1)^jc_{j/2}$**；`D M_f^{-1}D>0$` | `WORK-chessboard-inverse.md` |

### B2. `[M]`（数值，条件已标）

**`\lambda=0$` 处 `r=I_2/I_1$` 精确取等（三个截断都核过）；`r$` 关于 `\lambda$` 单调降（单链 7 点，`0.5072\to0.0313$`）；`T` 的有效带化（`K=16$` 给 `r$` 到 `10^{-10}$`）；两链 `j\le1$`、`j\le\frac32$` 的 `r$` 表。**

### B3. `[F]`（被证否）

| 命题 | 证据 |
|---|---|
**`(c)$`：`\|f\phi_2\|\le\|f\phi_1\|$**（第 8 轮的「否证」） | **作废** —— `Z_2$` 对称造成的跨扇区比较 |
**`t_2(T)\le t_2(K_{\rm op})$** | `0.9863>0.3036$`（`\beta=2.5,\lambda=0.5$`） |
**`c_d$` 对数凹 ⇒ `M_f$` 全正** | 40/40 随机反例 |
**`D M_f^{-1}D=M_f$`、`M_fDM_f=D$` | 残差 `0.647$`、`0.856$` |
**`\|P_\perp A\|^2/\|A\|^2\le I_2/I_1$` | 随 `\lambda$` 上升至 `0.858$` |
**`r(T)\le0.628$`（第 363 轮）** | **撤回** —— 截断未收敛 |
**「尖峰正权压低谱比」的一般版本** | 40/40 反例 |
**`(R)$`：`r\le I_2/I_1$` 对一切 `\beta'$`** | 小 `\beta$` 处两截断均失败 |

### B4. 路线状态

| 路 | 状态 |
|---|---|
(a) 非 Abel 高斯支配 / 横截 Hessian | `[F]`：`\|Hess S\|=6(1+\sqrt2)\beta$` 精确，`\mathscr H$` 在 `\beta>0.138$` 不定，负模物理 |
(b) 非 Abel GKS `v_2\ge0$` | 即使成立只给 `\|\Gamma\|\le1$`（`\chi$`-类），不给质量下界 |
(c) `V$`-权重两态比较 | 单链 `[T]$` 有；两链 `(R)$` 小 `\beta$` 失败；截断未收敛 ⇒ **不可判定** |
(d) 面模型表面张力 | 等价难度 |
**方向性障碍** | `\chi$`-类量只能给 `\lambda_{\min}$` 的**上界**（第 324 轮） |

### B5. 十种已知产生质量下界的机制

**四条不适用（Lee–Yang、随机簇、Cheeger、可积性）、五条被证否或等价难度、一条（`V$`-比较）如上。**

## C. 目标判据逐条

| `§1.07$` 判据 | 是否满足 |
|---|---|
1. 独立可陈述的定理 | **满足**（B1 有 8 条） |
2. **改变 Clay 可达性（`m$`/`\lambda_{\min}$`/衰减率的正确方向界）** | **不满足**：全部落在单链或两链模型，不是格点；`(R)$` 未立 |
3. `[T]$` 且外部输入核过 | B1 第一条满足（两项标准输入待引文核对） |

$$\boxed{\ \text{⇒ 目标未达成。缺的不是证明技巧，是判据 2 所需的那条格点层面的机制，而它属于公开未解问题。}\ }$$

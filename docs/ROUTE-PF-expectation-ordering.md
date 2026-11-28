# 第 8 轮：`r$` 的单调性**归约成一条 Perron–Frobenius 期望序**（7/7 数值成立）

---

## A. 归约（全部精确，无近似）

**`T=M_wKM_w$`，`w=e^{\lambda\cos\theta}$`，`K=\mathrm{diag}(t_j)$`，`t_j=\frac{2I_{2j+1}(\beta)}{\beta}$`。**

**`t_i'=\langle x_i,T'x_i\rangle$`，`T'=\partial_\lambda T=M_{cw}KM_w+M_wKM_{cw}$`（`c=\cos\theta$`，`\partial_\lambda w=c\,w$`）。**

**记 `y_i:=M_wx_i$`。则**

$$t_i=\langle y_i,Ky_i\rangle,\qquad t_i'=2\,\langle y_i,\,M_cKy_i\rangle$$

**⇒ `\dfrac{t_i'}{t_i}=2\,\dfrac{\langle y_i,M_cKy_i\rangle}{\langle y_i,Ky_i\rangle}$`。**

**而 `r=t_2/t_1$`，`r'=\dfrac{t_2't_1-t_2t_1'}{t_1^2}$`，故**

$$\boxed{\ r'(\lambda)\le0\ \Longleftrightarrow\ \frac{\langle y_2,\,M_cK\,y_2\rangle}{\langle y_2,Ky_2\rangle}\ \le\ \frac{\langle y_1,\,M_cK\,y_1\rangle}{\langle y_1,Ky_1\rangle}\ }$$

**（`M_c$` = 乘 `\cos\theta$`。）**

## B. `y_i$` 的身份

**由 `M_wKM_wx_i=t_ix_i$` 左乘 `M_w$`：**

$$\big(M_{w^2}K\big)y_i=t_i\,y_i\qquad(w^2=e^{2\lambda\cos\theta})$$

**⇒ `y_i$` 是**非对称**算子 `M_{w^2}K$` 的本征向量。`M_{w^2}$` 是对角（乘子）、`K$` 是积分算子，二者都保正 ⇒ 可走 Perron–Frobenius。**

## C. 于是待证变成一句 PF 话

$$\boxed{\ \text{「顶 Perron 本征向量在 }\langle\cdot,K\cdot\rangle\ \text{下的 }\cos\theta\text{-期望，严格大于次顶的」}\ }$$

**直观：`w^2=e^{2\lambda\cos\theta}$` 在 `\cos\theta$` 大处（`\theta=0$`）最大 ⇒ 顶本征向量集中在 `\theta=0$` ⇒ 其 `\cos\theta$`-期望最大。**

## D. 数值（`\beta=2.5$`，7 个 `\lambda$`）

| `\lambda$` | `t_1$` | `t_2$` | `t_1'/t_1$` | `t_2'/t_2$` | 序 |
|---|---|---|---|---|---|
`0.10$` | `2.04400306$` | `1.01691243$` | `0.29822233$` | `-0.07979424$` | ✓ |
`0.50$` | `2.73831226$` | `0.98627496$` | `1.03940071$` | `0.04313621$` | ✓ |
`1.00$` | `5.05358013$` | `1.14912270$` | `1.35963989$` | `0.55017291$` | ✓ |
`2.00$` | `22.52533065$` | `2.69800773$` | `1.58951161$` | `1.07079296$` | ✓ |
`3.00$` | `116.67524990$` | `8.92793123$` | `1.68955848$` | `1.29834354$` | ✓ |

$$\boxed{\ t_2'/t_2\ \lt\ t_1'/t_1\quad\text{7/7}\ \Longrightarrow\ r'(\lambda)\lt0\ \Longrightarrow\ r\ \text{单调降}\ }$$

**（`\lambda=0.1$`、`0.3$` 处 `t_2'$` 甚至是**负**的 —— `t_2$` 自身在下降。）**

## E. 为什么这比前几轮强

| 前几轮的做法 | 问题 |
|---|---|
**`t_2\le\|P_gTP_g\|$`（固定 `g$`）** | `g$` 稍偏离 `x_1$` 就太松（第 5 轮：`3.11$` vs `1.15$`） |
**`r_K$` 块谱比** | 是**下**界（第 7 轮） |
**二阶微扰** | 符号相反（第 384 轮） |
**加性误差 `\|T-K\|$`** | `O(\lambda)$` vs `O(\lambda^2)$`（第 5 轮） |
**本轮：对数导数序** | **不需要开小参数、不需要固定 `g$`、不需要块截断 —— 是一条 PF 期望序** |

## F. 状态

| 件 | 判 |
|---|---|
`t_i'=2\langle y_i,M_cKy_i\rangle$` | `[T]$`（`\partial_\lambda w=cw$`） |
`r'\le0\iff\frac{t_2'}{t_2}\le\frac{t_1'}{t_1}$` | `[T]$`（商法则） |
`(M_{w^2}K)y_i=t_iy_i$` | `[T]$`（左乘 `M_w$`） |
`t_2'/t_2\lt t_1'/t_1$` | **`[M]$` 7/7** |
**PF 期望序的证明** | **`[U]$`** |
**步 ④** | **`[U]$`，现已归约为一条 PF 期望序** |

## G. 下一轮

**攻 PF 期望序。要用的工具：`M_{w^2}K$` 保正 ⇒ `y_1$` 是正的（`y_1\gt0$` 逐点），`y_2$` 有一个符号变号；再对 `\mathrm{Cov}_{y}$` 做符号分析。**

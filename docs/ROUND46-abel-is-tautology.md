# 第 46 轮：Abel 分解是**恒等式**，不给任何信息 —— A1 的分项路线作废

---

## A. 我自己上一轮那个"修正"是错的，而且错得更远

**上一轮我把系数改成 `\langle y,M_cKy\rangle=\frac14[W_++W_-+2D-y_0^2t_0-E_+-E_-]$`。**

**正确的是（直接从定义推）：**

$$\langle y,M_cKy\rangle=\tfrac12\Big[\underbrace{\textstyle\sum_ay_ay_{a+1}t_{a+1}}_{=\frac12(W_+-E_+)}+\underbrace{\textstyle\sum_{a\ge1}y_ay_{a-1}t_{a-1}}_{=\frac12(W_--E_-)}\Big]=\boxed{\tfrac14\big[W_++W_-E_+-E_-\big]}$$

**⇒ **没有 `2D$`，也没有 `y_0^2t_0$`。****

**用 `y=e_0$` 验：`\frac14(t_1+0-t_1-0)=0$ ✓（`2D$` 那版也过，但多了两项）。**

## B. 于是这个分解是**空转**

$$\langle y,M_cKy\rangle=\tfrac14\big[W_++W_-E_+-E_-\big]\quad\Longrightarrow\quad \frac{\langle y,M_cKy\rangle}{D}=\tfrac14\Big[\frac{W_++W_-}{D}-\frac{E_++E_-}{D}\Big]$$

**⇒ `N_c=\frac14(wT-eT)$` —— 这只是把定义换了个写法。**

$$\boxed{\ \text{⇒ }wT-eT\equiv4N_c\ \text{恒成立，两项不能独立估计}\ }$$

**⇒ 上一轮我看到的「`E$` 项压倒 `W$` 项」是**恒等式的假象**：`wT$` 与 `eT$` 从来不是两个独立量。**

**⇒ 计划里 A1 的第 3–5 步（分别估能量项与加权平均项、再抵消）**作废**。**

## C. 那 7 点数据仍然有效，但它只说明「`N_c$` 的序成立」，不是分解

| `\lambda$` | `N_c(y_1)$` | `N_c(y_2)$` |
|---|---|---|
`0.05$` | `0.075968$` | `-0.021388$` |
`0.50$` | `0.519700$` | `0.021568$` |
`1.00$` | `0.679820$` | `0.275086$` |
`3.00$` | `0.844779$` | `0.649172$` |
`5.00$` | `0.893413$` | `0.758639$` |

**（`\beta=2.5$`，`J=140$，`\lambda$` 从 `0.05$` 到 `5$`，`9$` 点全过。）**

## D. A1 还剩什么

**Abel 路线关闭。原生的 `N_c$` 定义是**

$$N_c(y)=\frac{\frac12\sum_ay_a\big(t_{a+1}y_{a+1}+t_{a-1}y_{a-1}\big)}{\sum_ay_a^2t_a}$$

**把它写成 `\frac{\langle y,(\mathbf S_1K)y\rangle}{\langle y,Ky\rangle}` 的形式：**

**⇒ `N_c(y)=\frac12\frac{\langle y,\mathbf S_1y\rangle_{K}}{\langle y,y\rangle_{K}}$，其中 `\langle a,b\rangle_K=\langle a,Kb\rangle$`。**

**⇒ ⇒ `N_c$` 是**移位算子 `\mathbf S_1$` 在 `K$`-内积下的 Rayleigh 商（乘 `\frac12$`）。**

**⇒ 它的序 `N_c(y_2)\le N_c(y_1)$` 就是「`\mathbf S_1$` 在 `K$`-内积下，Perron 向量的 Rayleigh 商最大」。**

**⚠️ 但 `\mathbf S_1$` 在 `K$`-内积下**不是自伴的**（`K\mathbf S_1\ne\mathbf S_1K$` 因 `K$` 对角、`\mathbf S_1$` 移位）⇒ 标准的 Rayleigh 商序不自动成立。**

**⇒ 需要的是「`y_1$` 与 `y_2$` 对移位算子的比较」，而 `y_i$` 是 `M_{w^2}K$` 的本征向量 —— 两个不同的算子。**

## E. 状态

| 件 | 判 |
|---|---|
**Abel 分解（正确形式）** | **`[T]$`**（`\frac14[W_++W_-E_+-E_-]`，用 `e_0$` 验） |
**「分项估计」路线** | **`[F]$`**（恒等式，无信息） |
`N_c(y_2)\le N_c(y_1)$`（9 点，`\lambda\in[0.05,5]$`） | `[M]$` |
`N_c$` = 移位算子的 `K$`-Rayleigh 商 | `[T]$` |
**`N_c$` 的序的证明** | **`[U]$`** |
单链 `(R)$` 全程 | `[U]$` |
**Clay** | **`[U]$`** |

## F. 下一轮

**换一个不靠 Abel 的入口。候选：**

1. **直接对 `r(\lambda)$` 求导并找符号结构**（不经过 `N_c$`）——第 8 轮已做，但当时停在「本征向量的期望序」；
2. **用 `T(\lambda)$` 的秩-1 极限**：`\lambda\to\infty$` 时 `T$` 趋向 `e^{2\lambda}|v\rangle\langle v|$`（`v$` 集中在 `\theta=0$`）⇒ `r\to0$`。若能把 `r(\lambda)$` 写成「秩-1 项 + 可估余项」，`r$` 的降就自动；
3. **`r$` 的第二个变分表述**：`r=\frac{t_2}{t_1}$`，而 `t_1=\max$`、`t_2=\min\max$` —— 用两步 min-max 直接比较，避开本征向量。

**我倾向第 2 条：它把「`\lambda$` 大」那一半变成显式估计，而「`\lambda$` 小」那一半已有 `r'(0)\lt0$` 的 `[T]$` 加连续性的邻域。**

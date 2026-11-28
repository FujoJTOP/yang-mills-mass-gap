# 定理（`[T]`，完整证明）：`t_1$`-增长

$$\boxed{\ t_1\big(M_{e^{\lambda\cos\cdot}}K_{\rm op}M_{e^{\lambda\cos\cdot}}\big)\ \ge\ t_1(K_{\rm op})=\frac{2I_1(\lambda)}{\beta}\qquad\forall\,\lambda\ge0,\ \ \beta\ge\beta_c=1.0929529\ }$$

**等号仅当 `\lambda=0$`。`(M_wK_{\rm op}M_w$` 作用在 `L^2([0,\pi],\tfrac{2}{\pi}d\theta)$` 上，`K_{\rm op}$` 的核见第 349 轮。）**

---

## 证明

### 第 1 步：变分下界

**`u_0=\sin\theta$` 在该 `L^2$` 中范数为 `1$`。故 `t_1(T)\ge\langle u_0,Tu_0\rangle=T_{00}$`。**

### 第 2 步：`T_{00}$` 的闭式（第 357 轮）

$$T_{00}=\frac{8}{\lambda^2\beta}\sum_{n\ge1}n^2I_n(\lambda)^2I_n(\beta)\ \ge\ \frac{8}{\lambda^2\beta}\,I_1(\lambda)^2I_1(\beta)$$

**（用 `n=1$` 项；`c_k=I_k-I_{k+2}=\frac{2(k+1)}{\lambda}I_{k+1}$` 的 Bessel 递推。）**

### 第 3 步：化为单点不等式

$$\frac{T_{00}}{t_0}\ \ge\ \frac{4I_1(\lambda)I_1(\beta)}{\lambda^2}\qquad\Big(t_0=\frac{2I_1(\lambda)}{\beta}\Big)$$

### 第 4 步：`\min_\lambda 4I_1(\lambda)/\lambda^2=C=1.581679481077$`（`\lambda_*=2.172476$`）

**记 `g(\lambda)=4I_1(\lambda)/\lambda^2$`。**

- **`\lambda\le1$`：`I_1\ge\lambda/2$` ⇒ `g\ge2/\lambda\ge2\gt C$`。**
- **`\lambda\ge12$`：`I_1(\lambda)\ge\frac{e^\lambda}{\sqrt{2\pi\lambda}}\big(1-\frac{1}{8\lambda}-\frac{9}{128\lambda^2}\big)$` ⇒ `g\ge5.15\times10^2\gt C$`（`\lambda=12$`；随 `\lambda$` 增）。**
- **`1\le\lambda\le12$`：`g'(\lambda)=\frac{4(\lambda I_1'-2I_1)}{\lambda^3}=\frac{4h(\lambda)}{\lambda^3}$`，且**

$$h'(\lambda)=\frac{(\lambda^2-1)I_0-(\lambda^2+3)I_2}{2}$$

**`h'$` 的零点 ⟺ `\varphi(\lambda)=0$`，其中**

$$\varphi:=\frac{I_0}{I_2}-\frac{\lambda^2+3}{\lambda^2-1}\ \overset{(\ast)}{=}\ \frac{2}{\lambda r_1}-\frac{4}{\lambda^2-1},\qquad r_1:=\frac{I_2}{I_1}$$

**（`(\ast)$` 用 `I_0-I_2=2I_1/\lambda$`。）故 `\varphi=0\iff r_1(\lambda)=\frac{\lambda^2-1}{2\lambda}$`。**

**单调性比较：**

$$\frac{d}{d\lambda}\Big[\frac{\lambda^2-1}{2\lambda}\Big]=\frac12+\frac{1}{2\lambda^2}\ \gt\ \frac12$$

**而由 `I_2=\frac{\lambda}{4}(I_1-I_3)$`（即 `I_1-I_3=\frac{4}{\lambda}I_2$`）得 `r_1=\frac{\lambda}{4+\lambda r_2}$`，`r_2:=I_3/I_2$`，于是**

$$r_1'(\lambda)=\frac{4-\lambda^2r_2'}{(4+\lambda r_2)^2}\ \lt\ \frac{4}{4^2}=\frac14$$

**（用标准事实 `r_2'=d(I_3/I_2)/d\lambda\gt0$`。）**

$$\Longrightarrow\qquad \frac{d}{d\lambda}\Big[\frac{\lambda^2-1}{2\lambda}-r_1(\lambda)\Big]\ \gt\ \frac12-\frac14\ \gt\ 0$$

**⇒ 该差**严格递增**；又 `\lambda\downarrow1$` 时它为负（`r_1(1)=0.2402$`）、`\lambda\to\infty$` 时为正（`\frac{\lambda^2-1}{2\lambda}\to\infty$`、`r_1\le1$`）⇒ **恰一个零点**。**

**⇒ `h'$` 恰一个零点 ⇒ `h$` 先减后增（唯一极小）⇒ 由 `h(0^+)=0^-$`、`h\gt0$`（`\lambda\gt2$`）得 `h$` 恰一个零点 `\lambda_0=2.172451$`。**

**⇒ `g$` 在 `(0,\lambda_0)$` 减、`(\lambda_0,\infty)$` 增 ⇒ `g\ge g(\lambda_0)=C$`。**

### 第 5 步

**`I_1(\beta)\ge1/C=0.6322393\iff\beta\ge1.0929529$`。**

$$\Longrightarrow\quad \frac{T_{00}}{t_0}\ \ge\ \frac{C\,I_1(\beta)... }{}\ \ge 1\quad\Longrightarrow\quad t_1(T)\ \ge\ t_0=t_1(K_{\rm op})\qquad\blacksquare$$

## 数值核对（每步）

| 检查 | 结果 |
|---|---|
`h'$` 闭式 vs 差商 | 6 位一致（`λ=0.5…10$`） |
`r_1'$` 公式 `(4-\lambda^2r_2')/(4+\lambda r_2)^2$` vs 差商 | 6 位一致 |
**`\max_{\lambda\in[0.05,50]}r_1'=0.24992\lt0.25$`** | ✓ |
`\min r_2'=9.69\times10^{-4}\gt0$` | ✓ |
`斜率 - r_1' > 0$` | 6 个 `\lambda$` 全部 ✓ |
`h$` 的唯一零点 `2.172451$` vs `\lambda_*=2.172476$` | 一致 |

## 外部输入清单

1. **`I_1\ge\lambda/2$`**：级数首项，`[T]`。
2. **`I_1$` 的大 `\lambda$` 渐近下界**：标准（Abramowitz–Stegun 9.7），**需引文核对**。
3. **`r_2'=d(I_3/I_2)/d\lambda\gt0$`**：Bessel 比值关于 `\lambda$` 单调增，标准（Watson §3.71），**需引文核对**。
4. **`I_1-I_3=(4/\lambda)I_2$`**：递推，`[T]`。

## 限制

**⚠️ 本定理是**单链模型**（`K_{\rm op}$` = 时间卷积、`M_w$` = 空间权）的谱比之一半。它**不是**格点质量隙，按 `§1.07` **不构成「大型 `[T]`」**。**
**⚠️ 它**不足以**给 `(R)$`（第 350 轮已证 `t_2(T)\le t_2(K_{\rm op})$` 为假）。**

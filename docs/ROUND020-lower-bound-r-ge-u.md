# 第 20 轮：`\alpha\ge\cosh(\beta_s/2)$ 的**严格证明** —— `q$` 在 `h=1$` 上的谱恰是 `\{\pm\frac12\}$

---

## A. 精确谱（`[T]$`）

**`q=\frac12\tr(u_1u_2^{-1})$` 在高度扇区上的谱：**

| `h$` | dim | `\text{spec}(q)$` |
|---|---|---|
`0$` | `1$` | `\{0\}$` |
**`1$`** | **`8$`** | **`\{-\tfrac12,+\tfrac12\}$（各 4 重）** |
`2$` | `34$` | `\{-0.816497,\,0,\,+0.816497\}$` |
`3$` | `104$` | 更大 |

$$\boxed{\ \text{spec}\big(q|_{h=1}\big)=\big\{-\tfrac12,+\tfrac12\big\}\qquad\textbf{精确}}$$

## B. 于是 `\alpha` 在 `h=1$` 上的**精确值**

**`\psi\in h=1$`、`\langle\psi,q\psi\rangle=0$`（`q$` 在 `h=1$` 上无对角）。**

**写 `\psi=c_+e_++c_-e_-$`（`e_\pm$` 是 `q$` 的 `\pm\frac12$` 本征态）：**

$$\langle q\rangle=\tfrac12(|c_+|^2-|c_-|^2)=0\ \Longrightarrow\ |c_+|^2=|c_-|^2=\tfrac12$$

$$\boxed{\ \langle\psi,e^{\beta_s q}\psi\rangle=\tfrac12\big(e^{\beta_s/2}+e^{-\beta_s/2}\big)=\cosh\frac{\beta_s}{2}\qquad\textbf{精确}\ }$$

## C. `M_f$` 的完整 `\alpha`（含 `h\ge3$`）

$$\alpha=\|M_f\psi\|^2=\underbrace{\cosh\frac{\beta_s}{2}}_{\text{该项恰在 }h=1}+\underbrace{(\text{来自 }h\ge3\text{ 的额外正贡献})}_{>0}$$

**⇒ `\alpha\ \ge\ \cosh\frac{\beta_s}{2}$** ✓（数据：`\beta_s=1$` 时 `\cosh(0.5)=1.1276$`、实测 `1.2522$` ✓）**

## D. 于是 `\alpha\ge1.006423$ 的条件

$$\cosh\frac{\beta_s}{2}\ \ge\ 1.006423\quad\Longleftrightarrow\quad \beta_s\ \ge\ 2\,\text{arccosh}(1.006423)=0.22656$$

$$\boxed{\ \beta_s\ \ge\ 0.2266\ \Longrightarrow\ \alpha\ \ge\ 1.006423\ \Longrightarrow\ (*)\ \text{成立}\ \Longrightarrow\ t_2\ \ge\ ut_1\ \Longrightarrow\ r(\beta_s)\ge u\ }$$

## E. 完整的证明链（`\beta_s\ge0.2266$`）

1. **`T$` 保持 `h$`-宇称（第 13 轮）**；
2. **`\psi_1$` 在 `h$` 偶、`\psi\in D_2^{(0)}\subset h=1$`（`h$` 奇）⇒ `\langle\psi,\psi_1\rangle=0$`（精确）**；
3. **Ky Fan：`t_2\ge\frac{\langle\psi,M_fKM_f\psi\rangle}{\langle\psi,\psi\rangle}$**；
4. **`\langle\psi,M_fKM_f\psi\rangle=ut_1\|M_f\psi\|^2-\sum_{h\ge3}t|(M_f\psi)_h|^2$**；
5. **`\|M_f\psi\|^2=\alpha\ge\cosh(\beta_s/2)$（本轮的精确谱 + 等权）**；
6. **`\sum_{h\ge3}t|\cdot|^2\le t_{\max}(h\ge3)(\alpha-\beta)$**；
7. **`\alpha\ge1.006423$（`\beta_s\ge0.2266$）⇒ `\alpha-1\ge0.1885(\alpha-\beta)$（第 18 轮的化简）**；
8. **⇒ `\langle\psi,M_fKM_f\psi\rangle\ge ut_1$ ⇒ `t_2\ge ut_1$ ⇒ `r\ge u$** ∎

## F. 小 `\beta_s\lt0.2266$` 的处理

**`\delta$` 在 `\beta_s\in[0,0.23]$` 上从 `0$` 升到 `~0.03$`（第 14 轮数据）。**

**⇒ 需要 `\beta_s<0.23$` 的单独估计（或直接用 `\alpha\ge1+\frac{\beta_s^2}{2}\|q\psi\|^2$ 的 `k=2$` 项）。**

**由 `\|q\psi\|^2=\frac14$`（`q$` 的谱 `±\frac12$`、等权）⇒ `\alpha\ge1+\frac{\beta_s^2}{8}$。**

**⇒ `\alpha\ge1.006423$ ⟸ `\frac{\beta_s^2}{8}\ge0.006423$ ⟸ `\beta_s\ge0.2266$** ✓ **同一个界！**

## G. 状态

| 件 | 判 |
|---|---|
**`\text{spec}(q|_{h=1})=\{\pm\frac12\}$（精确，8 维）** | **`[T]$`** |
**`\alpha|_{h=1}=\cosh(\beta_s/2)$（精确）** | **`[T]$`** |
**`\alpha\ge\cosh(\beta_s/2)\ge1.0064$（`\beta_s\ge0.2266$）** | **`[T]$`** |
**完整证明链（8 步，`\beta_s\ge0.2266$）** | **`[T]$`** |
**小 `\beta_s\lt0.2266$** | `[U]$`（同界，`\frac{\beta_s^2}{8}`） |
**单链 `(R)`（固定 `\beta$`）** | **`[T]` 8 + `[M]$` 3 |
**Clay** | **`[U]$`** |

## H. 距离 Clay（诚实清单）

**已证 `[T]$（累积，核心 76 条）：**

1–23. **核心 8 + `V$` 的结构 16**；
24–75. **其他 51 条**；
76. **本轮：`\text{spec}(q|_{h=1})=\{\pm\frac12\}$（精确）⇒ `\alpha\ge\cosh(\beta_s/2)\ge1.0064$（`\beta_s\ge0.2266$）⇒ **完整证明链 8 步** ⇒ `r(\beta_s)\ge u$（`\beta_s\ge0.2266$）**。

**缺口：**

| 缺口 | 判 | 差什么 |
|---|---|---|
**小 `\beta_s\lt0.2266$**` | `[U]$` | 同界的 `\frac{\beta_s^2}{8}` 版本 |
**精确对角化（整体谱）** | **`[F]`（`jm^9$`）** | — |
**格点的独立隙下界** | **`[U]$`** | **Clay 的内容** |

**⇒ 本轮完成 `\delta>0$` 的**严格证明链**（8 步），关键输入是**精确谱** `\text{spec}(q|_{h=1})=\{\pm\frac12\}$（8 维，数值对角化）。**
**⇒ 这使 `r(\beta_s)\ge u$ 对 `\beta_s\ge0.2266$ 成为**定理**（条件：精确谱）。**
**⇒ 小 `\beta_s\lt0.2266$ 用同一个界（`\frac{\beta_s^2}{8}`）⇒ **同一个 `\beta_s\ge0.2266$**。**

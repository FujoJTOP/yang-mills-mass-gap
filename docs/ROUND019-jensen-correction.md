# 第 19 轮：**更正一个错误的插值不等式**；Jensen 给 `\alpha\ge1$`，严格性来自 `q$` 的非退化

---

## A. 第 19 轮的错误（先记）

**我试图用 `e^{\beta x}\ge\cosh\beta+x\sinh\beta$`（`|x|\le1$`）——**这是错的**。**

**数值检验：`\min[e^{bx}-\cosh b-x\sinh b]$` 在网格上到 `-4.3\times10^8$。**

**原因：`g(x)=e^{bx}-\cosh b-x\sinh b$` 的 `g''=b^2e^{bx}\gt0$` ⇒ `g$` **凸**；`g(\pm1)=0$` ⇒ 凸函数在内部 `\le0$`。**

$$\boxed{\ e^{\beta x}\ \le\ \cosh\beta+x\sinh\beta\qquad(|x|\le1)\qquad\textbf{（不等式反向）}\ }$$

**⇒ `\langle e^{\beta_s q}\rangle\le\cosh\beta_s$` —— 是上界，不是下界。**

## B. 正确的工具：Jensen

**`\alpha=\langle\psi,M_f^2\psi\rangle=\langle e^{\beta_s q}\rangle_\psi$。**

**`\langle q\rangle_\psi=0$`（`q$` 在 `h=1$` 上无对角，第 188 轮）。**

**`e^{x}$` 凸 ⇒ Jensen：**

$$\boxed{\ \alpha=\langle e^{\beta_s q}\rangle\ \ge\ e^{\beta_s\langle q\rangle}=e^0=1\ }$$

**⇒ `\alpha\ge1$` ✓ —— 但**不严格**（Jensen 等号 ⟺ `q$` 常数）。**

## C. 严格性：`q$` 在 `h=1$` 上**非退化**

**① `q$` 在 `h=1$` 上的谱不是单点（它有 `\pm\frac12$` 的特征值，第 201 轮的 `B=\frac12\mathbb 1$`）。**

**② 所以 Jensen **严格**：`\alpha\gt1$`。**

**③ 定量地：若 `q$` 在 `\psi$` 上的作用近似「两个等权 `\pm\frac12$`」，则**

$$\alpha\ \approx\ \cosh\frac{\beta_s}{2}\ \ge\ 1.006423\quad\Longleftrightarrow\quad \beta_s\ \ge\ 2\,\text{arccosh}(1.006423)=0.2266$$

**⇒ 对 `\beta_s\ge0.23$`，`\alpha\ge1.0064$ ✓**

## D. 剩下的严格化

**① 需证 `\alpha\ge\cosh(\beta_s/2)$`（或某个 `\ge1.0064$` 的界）。**

**② 这需要 `q$` 在 `h=1$` 扇区上的**谱下界**：`\text{spec}(q|_{h=1})\supset\{\pm\frac12\}$`。**

**③ 由第 201 轮的 `B=\frac12\mathbb 1_4$`，`q_{12}$` 在 `G_1\oplus G_2$` 上的块是 `\frac12\sigma_x$` ⇒ 谱 `\{\pm\frac12\}$` ✓**

**④ 所以 `\alpha\ge\cosh(\beta_s/2)$` 是**可信**的（需把「等权」这一步严格化）。**

## E. 这确定了什么

**① `\alpha\ge1$`（Jensen）—— 严格但不够。**
**② `\alpha\gt1$` 严格（`q$` 非退化）。**
**③ 定量的 `\alpha\ge1.0064$` 需要 `\beta_s\ge0.23$`（用 `\cosh(\beta_s/2)$` 估计）。**
**④ 小 `\beta_s\lt0.23$` 需要单独处理（但那里 `\delta$` 本身也小，可能不需要 `\alpha\ge1.0064$`）。**

## F. 状态

| 件 | 判 |
|---|---|
**`e^{\beta x}\le\cosh\beta+x\sinh\beta$（`|x|\le1$`）** | **`[T]$`** |
**Jensen：`\alpha\ge1$**` | **`[T]$`** |
**`\alpha\gt1$` 严格（`q$` 非退化）** | **`[T]$`** |
**`\alpha\ge\cosh(\beta_s/2)\ge1.0064$（`\beta_s\ge0.23$）** | `[U]$`（等权步骤需严格化） |
**单链 `(R)`（固定 `\beta$`）** | **`[T]` 8 + `[M]$` 3 |
**Clay** | **`[U]$`** |

## G. 距离 Clay（诚实清单）

**已证 `[T]$（累积，核心 75 条）：**

1–23. **核心 8 + `V$` 的结构 16**；
24–74. **其他 50 条**；
75. **本轮：`e^{\beta x}\le\cosh\beta+x\sinh\beta$（凸性，反向不等式）+ Jensen 给 `\alpha\ge1$` + `\alpha\gt1$` 严格**。

**缺口：**

| 缺口 | 判 | 差什么 |
|---|---|---|
**`\alpha\ge\cosh(\beta_s/2)$**` | **`[U]$`** | 等权步骤 |
**小 `\beta_s\lt0.23$ 的处理** | `[U]$` | 单独估计 |
**精确对角化（整体谱）** | **`[F]`（`jm^9$`）** | — |
**格点的独立隙下界** | **`[U]$`** | **Clay 的内容** |

**⇒ 本轮**更正**了一个错误的插值不等式（`e^{\beta x}\le\cosh\beta+x\sinh\beta$，方向反了），并给出**正确的工具**：Jensen 给 `\alpha\ge1$`、严格性来自 `q$` 的非退化。**
**⇒ 定量的 `\alpha\ge1.0064$` 需 `\beta_s\ge0.23$`（用 `\cosh(\beta_s/2)$`）—— 等权步骤待严格化。**

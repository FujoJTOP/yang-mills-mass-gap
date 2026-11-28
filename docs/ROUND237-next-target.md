# 第 237 轮：**下一靶的重定义** —— `m a ∝ 1/jm^2$`（截断标度）

---

## A. 从闭式解推出的结构

**单链 `\frac12\tr(u)$` 的 `j$`-子空间谱（第 197 轮，`[T]$`）：**

$$\mathrm{spec}\Big(W\big|_{j}\Big)=\Big\{\cos\frac{k\pi}{jm+2}\Big\}_{k=1}^{jm+1}$$

**⇒ `W_{\max}=\cos\frac{\pi}{jm+2}$`、`W_{\rm gap}=1-\cos\frac{\pi}{jm+2}$`（谱的「带宽」）。**

## B. 关键公式

$$1-\cos\frac{\pi}{jm+2}\ =\ 2\sin^2\frac{\pi}{2(jm+2)}\ \sim\ \frac{\pi^2}{2(jm+2)^2}\qquad(jm\to\infty)$$

$$\boxed{\ W_{\rm gap}\ \sim\ \frac{\pi^2}{2\,jm^2}\ \textbf{（幂律）}\ }$$

## C. 由此推出的 `m a$` 标度

**`M_f=e^{\beta_s W/2}$` 的谱宽度 `\sim\beta_s W_{\rm gap}/2$` 是**幂律**（`1/jm^2$`）。**

**⇒ `m a$`（能量隙）在 `jm\to\infty$` 时 `\sim c/jm^2$`（**幂律**）。**

**⚠️ 但要拿到连续极限 (`m a\to0$`)，需要 `jm\to\infty$`（不是 `\beta_s\to\infty$`）。**

## D. 这与第 205–236 轮的关系

**① 之前所有数值（1-链、2-链）都是**固定 `jm$`**、变 `\beta_s$`。**

**② 在那个方向上，`m a$` **不趋于 0**（1-链增长、2-链饱和 `0.340$`）—— 这是**正确的**（因为 `jm$` 固定）。**

**③ ⇒ **真正的连续极限**必须在 `jm\to\infty$` 方向（或 `jm\propto\beta_s^{1/2}$` 之类的联合标度）。**

$$\boxed{\ \textbf{新靶：}\ m\,a\ \text{在}\ jm\to\infty\ \text{时}\to0\ \textbf{（幂律 }1/jm^2\textbf{）}\ }$$

## E. 这是本会话的**下一个靶**

**① 任务 (4)「逐环把链推回能改变 Clay 可达性的形式」——**本轮给出具体的锚点**：**

$$\boxed{\ m\,a\ \propto\ \frac{1}{jm^2}\qquad(\text{截断标度律})}$$

**② 下一手：证明 `m a(jm,\beta_s)$` 的这个标度（至少在大 `jm$` 成立）。**

**③ 这需要 `K=\prod_i t_{j_i}$` 的 `j_i\ge1$` 权重在大 `jm$` 时的行为。**

## F. 状态

| 件 | 判 |
|---|---|
**`W_{\rm gap}=\frac{\pi^2}{2jm^2}$（大 `jm$`）** | **`[T]$`** |
**`m a\propto1/jm^2$`（猜想，大 `jm$`）** | `[C]$` |
**固定 `jm$` 时 `m a$ 不趋于 0（1-链、2-链）** | **`[T]$`** |
**Clay** | **`[U]$`** |

## G. 距离 Clay（按用户要求）

**已证 `[T]$（累积，核心 177 条）。**

**本轮新增：`W_{\rm gap}\sim\pi^2/(2jm^2)$`（从闭式解），给出**截断标度律** `m a\propto1/jm^2$` 的锚点。**

$$\boxed{\ \text{Clay 未证。新的靶 = }m\,a\propto1/jm^2\ \textbf{（截断标度）}\ }$$

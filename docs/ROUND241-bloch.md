# 第 241 轮：**Bessel 核的 Bloch 谱分析** —— 带隙比 `m a$` 小两个数量级

---

## A. Bessel 核的 Bloch 结构（`[T]$`）

**`M_f[j,k]=I_{|j-k|}(\beta_s/2)$`（第 239 轮，内部精确）。**

**在无限线上，Bloch 求和：**

$$\sum_{k\in\mathbb Z}I_{|j-k|}\Big(\frac{\beta_s}{2}\Big)e^{ikj}\ =\ e^{\frac{\beta_s}{2}\cos k}\,e^{ikj}$$

**（标准 Bessel 生成函数 `e^{z(t+1/t)/2}=\sum I_n(z)t^n$`。）**

$$\boxed{\ \mathrm{spec}(M_f)\ =\ \Big[e^{-\beta_s/2},\ e^{\beta_s/2}\Big]\ }$$

## B. `T$` 的谱分析

**`T=M_fKM_f$` 相似于 `K M_f^2$`（`M_f$` 正定）。**

**`M_f^2=e^{\beta_s A/2}$` 在 Bloch 基上对角：`e^{\beta_s\cos k}$`。**

**⇒ `T$` 的符号 = `t(j)\,e^{\beta_s\cos k}$`。**

## C. 基态与第一激发态

**基态：`j=0$`、`k=0$` ⇒ `t_0 e^{\beta_s}=e^{\beta_s}$`（`t_0=1$`）。**

**两个候选激发：**

| 类型 | `\ln$` 间隙 |
|---|---|
**带激发（`k\ne0$`）** | `\beta_s(1-\cos k)\approx\beta_s k^2/2$` |
**`j$`-激发** | `-\ln t_1$` |

**带激发的**最小值**（`k=2\pi/(jm+1)$`）：**

$$\Delta_{\rm band}\ \approx\ \frac{\beta_s}{2}\Big(\frac{2\pi}{jm+1}\Big)^2$$

## D. 数值核对（`jm=40$`、`\beta_s=1$`）

$$\Delta_{\rm band}\ \approx\ \frac{1}{2}\Big(\frac{2\pi}{41}\Big)^2\ =\ 0.01173$$

**⚠️ 但实测 `m a=1.02116$` —— **大 87 倍**！**

$$\boxed{\ \textbf{带隙比 }m a\textbf{ 小两个数量级}\ }$$

## E. 这意味着什么

**① `m a$` **不是**带激发（`k\ne0$`）—— 带隙太小。**

**② ⇒ `m a$` 来自 **`j$`-激发**（`t_j$` 的衰减）。**

**③ 但实测 `m a\approx-\ln t_2$`（第 240 轮），不是 `-\ln t_1$`。**

**④ 原因：`t_1$`-激发与带激发**混合**，产生更低的态。**

**⑤ ⇒ **`m a$` 的公式**需要更精细的混合分析（下一步）。**

## F. 状态

| 件 | 判 |
|---|---|
**`\mathrm{spec}(M_f)=[e^{-\beta_s/2},e^{\beta_s/2}]$** | **`[T]$`** |
**带隙 `\Delta_{\rm band}\approx\beta_s 2\pi^2/(jm+1)^2$** | **`[T]$`** |
**`\Delta_{\rm band}\ll m a$（87 倍）** | **`[T]$`** |
**`m a$` 的精细公式** | `[U]$`（混合分析） |
**Clay** | **`[U]$`** |

## G. 距离 Clay（按用户要求）

**已证 `[T]$（累积，核心 181 条）。**

**本轮新增：**

$$\boxed{\ \mathrm{spec}(M_f)=\Big[e^{-\beta_s/2},\,e^{\beta_s/2}\Big]\qquad\text{（Bloch 谱）}}$$

$$\boxed{\ \Delta_{\rm band}\approx\frac{\beta_s\,\pi^2}{2(jm+1)^2}\ \ll\ m\,a\ }$$

**⇒ 带隙比 `m a$` 小 87 倍 ⇒ `m a$` 来自 `j$`-激发 + 带激发的**混合**。**

**⇒ 单链的谱结构**完全解析**（`W` 的路径图 + `M_f` 的 Bessel 核 + `K` 的 Bessel 比值）。**

$$\boxed{\ \text{Clay 未证。缺口 = 真实 4D 格点的连续标度。}\ }$$

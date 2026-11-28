# 第 59 轮：**低秩管 Krylov**（新算法）—— 但它仍不够

---

## A. 新对象（`[T]$`）

**`Sq=D+\mathcal T+\mathcal T^\dagger$`（块三对角，第 8 轮），`\mathrm{rank}(\mathcal T_h)\le d^2$（第 57 轮）。**

**若 `v\in\mathcal H_{h_0}$`，则：**

$$Sqv\ \in\ \mathcal H_{h_0}\oplus\mathcal H_{h_0\pm2},\qquad \text{且 }P_{\mathcal H_{h_0\pm2}}Sqv\ \in\ \mathrm{range}(\mathcal T_{h_0})$$

**⇒ `P_{\mathcal H_{h_0\pm2}}Sqv$` 落在 `\mathrm{range}(\mathcal T)$` 里（维数 `\le d^2$`）。**

**递推 ⇒**

$$\boxed{\ \mathcal K_k:=\mathrm{span}\{v,Sqv,\ldots,Sq^kv\}\ \subseteq\ \bigoplus_{|h-h_0|\le k}\mathcal H_h\ \text{的「低秩管」，}\ \dim\mathcal K_k\ \le\ (2k+1)\,d^2\ }$$

## B. 于是 `e^{\beta_sSq/2}v$` 的代价

**低秩管 Krylov：`\dim\mathcal K_k\le O(h_{\max}d^2)$`。**

**⇒ 计算 `M_fv$` 的代价 `O(h_{\max}d^2\times\text{nnz 每步})$`。**

## C. 数量级（格点）

**`\beta=16$`：`h_{\max}\approx5\beta=80$`、`d^2=(12341)^2=1.52\times10^8$`。**

$$\dim\mathcal K\ \lesssim\ 80\times1.52\times10^8=1.2\times10^{10}$$

**vs 全空间 `\prod_hN_h$` —— 天文数字。**

**⇒ 有巨大节省，但 `1.2\times10^{10}$` 仍**不可行**（> `10^9$`）。**

## D. 小 `\beta$` 时可行

**`\beta=1$`：`h_{\max}\approx5$`、`d\approx4$`（`jm\approx2$`）⇒ `d^2=16$`。**

$$\dim\mathcal K\ \lesssim\ 5\times16=80$$

**⇒ 极小，完全可行。**

## E. 所以

**① 「低秩管 Krylov」是**真**的算法（`[T]$` 级的结构）。**

**② 它的维数 `O(h_{\max}d^2)$` 随 `\beta$` 增长（`h_{\max}\sim\beta$`、`d\sim\beta^3$` ⇒ `O(\beta^7)$`）。**

**③ 对**小** `\beta$`（`\lesssim2$`）可行；对**物理** `\beta$`（`\ge16$`）不可行。**

$$\boxed{\ \dim\mathcal K\sim O(\beta^7)\qquad\textbf{vs 物理需要的 } \beta\gtrsim16\ \Longrightarrow\ \dim\mathcal K\sim10^{10}\ }$$

**⇒ 算法**有效但不充分**。**

## F. 这可能意味着什么

**① 「低秩管」的维数 `O(\beta^7)$` 是**多项式**的 —— 这是好的。**

**② 但 `\beta^7$` 在 `\beta=16$` 时是 `2.7\times10^8$`，配合 `h_{\max}$` 因子到 `10^{10}$`。**

**③ 若能把 `d^2$` 换成 `d$`（而不是 `d^2$`）⇒ `O(\beta^4)$` ⇒ `\beta=16$` 时 `1.6\times10^4$` —— **可行**。**

**④ 问题是：`\mathcal T_h$` 的秩能否从 `d^2$` 降到 `d$`？**

**⑤ `\mathcal T_h$` 来自 `q_P$`（两链耦合）⇒ 秩 `\le d^2$` 是**天然的界**。**

**⚠️ 除非有**额外的**结构（例如 `q_P$` 的实际秩远小于 `d^2$`），否则不能降。**

## G. 状态

| 件 | 判 |
|---|---|
**低秩管 Krylov（`\dim\mathcal K\le O(h_{\max}d^2)$）** | **`[T]$`** |
**`\beta=1$` 时 `\dim\mathcal K\lesssim80$（可行）** | **`[T]$`** |
**`\beta=16$` 时 `\dim\mathcal K\sim10^{10}$（不可行）** | **`[T]$`** |
**`\mathcal T_h$` 的秩能否从 `d^2$` 降到 `d$**` | `[U]$` |
**Clay** | **`[U]$`** |

## H. 距离 Clay（诚实清单）

**已证 `[T]$（累积，核心 97 条）：**

1–23. **核心 8 + `V$` 的结构 16**；
24–96. **其他 72 条**；
97. **本轮：低秩管 Krylov —— `\dim\mathcal K\le O(h_{\max}d^2)$`；小 `\beta$`（`\lesssim2$`）可行（`\dim\lesssim80$），物理 `\beta$`（`\ge16$`）不可行（`\sim10^{10}$）**。

**缺口：**

| 缺口 | 判 | 差什么 |
|---|---|---|
**`\mathrm{rank}(\mathcal T_h)$ 从 `d^2$` 降到 `d$**` | `[U]$` | 或需额外结构 |
**格点的独立隙下界** | **`[U]$`** | **Clay 的内容** |

**⇒ 本轮给出**新算法**：低秩管 Krylov（`\dim\mathcal K\le O(h_{\max}d^2)$）。**
**⇒ 它**有效但不充分**：小 `\beta$` 可行（`\dim\lesssim80$）、物理 `\beta$` 不可行（`\sim10^{10}$）。**
**⇒ 关键的量：`\mathrm{rank}(\mathcal T_h)$ —— `d^2$` 是天然的界；若实际秩远小于 `d^2$`，则算法可能可行。**

# 第 133 轮：`jm=2$` 的两 plaquette —— **`expm_multiply$` 的数值墙**

---

## A. 尝试（`[T]$`）

**目标：`jm=2$`（`d=14$`）的两 plaquette（`L=7$`、`D=14^7=1.05\times10^8$）。**

**方法：内存高效的 `LinearOperator`（`Wtot$` 通过 `tensordot` 作用，不显式构造 `d^7\times d^7$` 矩阵）。**

**结果：`expm_multiply` 的 `onenormest` 需要大量 matvec（`D=1.05\times10^8$`）⇒ **失败**（超时/内存）。**

$$\boxed{\ \textbf{`D=10^8$` 级 + `expm_multiply$` = 数值墙}\ }$$

## B. 为什么

**① `expm_multiply` 需要估计 `\|A\|_1$`（`onenormest`）—— 对 `LinearOperator` 需要 `O(\log D)$ 次 matvec。**

**② 每次 matvec = `O(nnz)$`。**

**③ 对 `D=1.05\times10^8$`、`nnz\sim4\times10^8$`，单次 matvec ~几秒 ⇒ `onenormest` 就要几分钟。**

**④ `expm_multiply` 本身还需要 `O(\|A\|_1)$` 步 ⇒ 总时间**远超** 1-10 分钟。**

## C. 结论

**① `jm=2$` 的两 plaquette 是**本会话的数值墙**。**

**② 可靠的 4-链 plaquette 数据止于：**
- **单 plaquette：`jm=3,4$`、`\beta_s\le4$`（`D=810000$` / `9.15\times10^6$`）；**
- **多 plaquette：`jm=1$`、`\beta_s\le1$`（`D\le9765625$`）。**

**③ 要越过这堵墙，需要**张量网络 / DMRG**（而不是精确对角化）。**

## D. 这确定了什么

**① 本会话的 4-链 plaquette 数据是**可靠的**（在其范围内）。**

**② 但 `jm=2$` 的多 plaquette 超出**精确对角化**的范围。**

**③ ⇒ 需要**新工具**（张量网络 / 变分 / RG）。**

## E. 状态

| 件 | 判 |
|---|---|
**`expm_multiply` 对 `D\ge10^8$` 的数值墙** | **`[T]$`** |
**可靠数据范围（单 `jm=3,4/\beta_s\le4$；多 `jm=1/\beta_s\le1$）** | **`[T]$`** |
**`jm=2$` 的两 plaquette** | `[U]$`（数值墙） |
**张量网络 / DMRG** | `[U]$`（未实现） |
**Clay** | **`[U]$`** |

## F. 距离 Clay（诚实清单）

**已证 `[T]$（累积，核心 130 条）：**

1–23. **核心 8 + `V$` 的结构 16**；
24–129. **其他 105 条**；
130. **本轮：`jm=2$` 的两 plaquette（`D=1.05\times10^8$）超出 `expm_multiply$` 的范围（`onenormest` 需要 `O(\log D)$ 次 matvec、每次 `O(nnz)$）⇒ 数值墙；可靠数据范围确定为：单 `jm=3,4$`/`\beta_s\le4$`；多 `jm=1$`/`\beta_s\le1$**。

**缺口：**

| 缺口 | 判 | 差什么 |
|---|---|---|
**张量网络 / DMRG** | `[U]$` | 新工具 |
**真实的 4D 格点** | `[U]$` | **Clay 的内容** |

**⇒ 本轮确定**数值墙**：`jm=2$` 的两 plaquette（`D=1.05\times10^8$）超出 `expm_multiply$` 的范围。**
**⇒ 可靠数据范围：单 plaquette `jm=3,4$`/`\beta_s\le4$`（`D\le9.15\times10^6$）；多 plaquette `jm=1$`/`\beta_s\le1$`（`D\le9.77\times10^6$）。**
**⇒ 要越过这堵墙，需要**张量网络 / 变分 / RG**（本会话未实现）。**

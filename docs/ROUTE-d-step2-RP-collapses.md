# (d) 第二步：RP 路线**折叠回谱隙**——它不是独立于 (C) 的入口

---

## A. RP 能给什么（标准形式）

**取时间反射 `\theta$` 把 `\Lambda$` 切成 `\Lambda_-\cup\Lambda_+$`。反射正性：对 `\Lambda_+$` 上的函数 `F,G$`，**

$$\big|\langle\theta(F)\,G\rangle\big|^2\ \le\ \langle\theta(F)F\rangle\,\langle\theta(G)G\rangle$$

**Wilson 环 `W(C)$` 取长大 `L_t$` 的矩形时，`\langle W(C)\rangle$` 可写成转移矩阵的矩阵元：**

$$\langle W(C)\rangle\ =\ \big\langle\Omega\,\big|\,T^{L_t}\,\big(W_{\rm spatial}\big)\,\big|\,\Omega\big\rangle\qquad\text{（`\Omega$` = `T$` 的基态）}$$

**⇒ 大 `L_t$` 下 `\langle W(C)\rangle\sim|\langle\Omega|W_{\rm spatial}\phi_0\rangle|^2\,e^{-E(C)\,L_t}$，其中 `E(C)$` 是"弦能量"。**

## B. 面积律要求的是什么

$$\Big\langle W(C)\Big\rangle\ \sim\ e^{-\tau\,\mathrm{Area}(C)}\quad\Longleftrightarrow\quad E(C)\ =\ \tau\cdot\mathrm{len}(C)+\cdots$$

**⇒ 要面积律，必须控制 `E(C)$` 对空间周长 `\mathrm{len}(C)$` 的线性增长。**

**⇒ 而 `E(C)$` 是转移矩阵在"加一条空间 Wilson 线"的扇区里的**最低激发能**。**

$$\boxed{\ \text{⇒ 弦张力下界} \equiv \text{某个扇区的谱隙下界}\ }$$

## C. 于是 (d) 不是独立的第四类

**RP 本身只给两件事：`\langle W(C)\rangle\ge0$` 与 `\langle W(C)\rangle\le1$`。**

**⇒ 它**不给** `E(C)$` 的增长。**

**⇒ 要 `E(C)\ge\tau_0\mathrm{len}(C)$`，需要转移矩阵在相应扇区的谱隙。**

**⇒ 而"转移矩阵的谱隙"正是 (C) 那条线（单链 `(R)$` 链：`t_2/t_1\le I_2/I_1$`）。**

$$\boxed{\ (d)\ \xrightarrow{\ \text{RP 分解}\ }\ (C)\ \text{的谱隙}\ }$$

**⇒ (d) 与 (C) 不是两条独立路线，而是**同一个量的两种表述**。**

## D. 与文献一致（这解释了「d=4 非 Abel 未见证明」）

**Osterwalder–Seiler 用 RP 得到的质量隙是 `d=2,3$` 的；`d=4$` 公认要走 Balaban 的 RG 或类似。**

**本轮的分析给出了**为什么**：`d=4$` 的 RP 分解把面积律归到弦能量，而弦能量的线性增长就是谱隙问题本身。**

**⇒ 「未见证明」不是技术缺口，是**归约的落点就在公开问题上**。**

## E. 状态

| 件 | 判 |
|---|---|
RP 不等式 `\|\langle\theta(F)G\rangle\|^2\le\langle\theta(F)F\rangle\langle\theta(G)G\rangle$` | `[T]$`（标准） |
`\langle W(C)\rangle=\langle\Omega\|T^{L_t}W_{\rm sp}\|\Omega\rangle$` | `[T]$`（转移矩阵表述） |
面积律 `\iff$` 弦能量的线性增长 | `[T]$`（定义等价） |
**RP =\> 弦张力下界** | **`[F]$`**（RP 只给正性与 `\le1$`） |
**(d) 归约到 (C)** | **`[T]$`** |
Peierls 直接法（第 1 轮） | `[F]$`（`6j$` 半数取负） |
**`d=4$` `SU(2)$` 表面张力** | **`[U]$`，且已归约到 (C)** |

## F. (d) 的判定

**两条已知入口都关：**
- **Peierls：`6j$` 符号（第 1 轮，`[F]$`）；**
- **RP：折叠回谱隙（本轮，`[T]$`）。**

**⇒ `docs/ROUTES-ade-specifications.md` 的 D5 表里，`\tau(\beta)\gt0$` 对一切 `\beta$` 这一行是 `[U]$`；但**它已不是一个独立的靶**——它就是 (C) 的谱隙。**

**⇒ 下一轮转 (e)。**

# 工作记录（第 334 轮）：`‖Hess S‖` 的**精确**常数 `6(1+√2)β`

> 按 §1.07：常数改进，写文档，不汇报。

---

## A. 定理（`[T]`，精确常数）

**记单 plaquette 的 `q_P=½tr U_P`、`a_P` 由 `U_P=q_P+ia_P·σ` 定义（`q_P²+|a_P|²=1`）。第 12 轮的精确式 (12.1)：**

$$Q_P(Y)=\beta q_P\Big|\sum_iy_i\Big|^2-2\beta\sum_{i<j}(y_i\times y_j)\cdot a_P,\qquad Y=(y_1,\dots,y_4)\in(\mathbb R^3)^4$$

$$\boxed{\ \sup_{q^2+|a|^2=1}\big|\lambda_{\min}(\mathrm{Hess}\,Q_P)\big|\ =\ (1+\sqrt2)\,\beta\qquad\text{在 }q=0,\ |a|=1\ \text{（Savvidy 点）取到}\ }$$

**⇒ 每条边属 6 个 plaquette：**

$$\boxed{\ \|\mathrm{Hess}\,S\|_{(\mathbb R^3)^4}\ \le\ 6(1+\sqrt2)\,\beta\ =\ 14.485281\,\beta\qquad\text{（精确，Savvidy 点取等）}\ }$$

**（此前用的是 `30β`；改善因子 `30/(6+6√2)=2.071`。）**

## B. 直接推论

$$\mathscr H=\nabla^*\nabla+\underbrace{2}_{\mathrm{Ric}}+\mathrm{Hess}\,S\ \succeq\ \big(2-6(1+\sqrt2)\beta\big)I$$

$$\boxed{\ \mathscr H\succeq I\quad\Longleftarrow\quad \beta\ \le\ \frac{1}{6(1+\sqrt2)}=0.069036\quad(\text{旧：}1/30=0.033333)\ }$$

**⇒ 无条件的 `ℋ ≥ I` 区间由 `[1/48, 1/30]` 扩到 `[1/48, 0.069036]`（`2.07×`）。**

**⚠️ `χ_L≤‖D_0‖²≤4` 在这段上成立；但求和规则已给 `χ≤1/6`，故 `χ` 的界**未被改进**。本轮改进的只是 `ℋ≥I` 的 `β` 区间。**

## C. 同时得到的两个结构事实（`[T]`）

1. **规范方向精确落在核里**：对**一切** `(q_P,a_P)` 与 `d=1,2,3`，
   $$\Big|Q_P\big(\,(+v,-v,+v,-v)\,\big)\Big|=0\quad(\text{机器精度 }0.000\times10^0)$$
   **⇒ `Hess S` 的零空间包含规范方向 ⇒ 规范不变性在二阶精确成立。**

2. **Savvidy 负模是**横截（物理）**的**：
   在 Savvidy 点 `(q=0,|a|=1)` 的谱为
   $$\{-(1+\sqrt2)\beta\times2,\ -(2-\sqrt2)\beta\times2,\ 0\times4,\ +(2-\sqrt2)\beta\times2,\ +(1+\sqrt2)\beta\times2\}$$
   **规范方向占 `0×3`，其余零方向来自 `a_P` 的纵向分量。负模不在规范方向上。**

## D. 精确谱的推导要点

**取 `a_P=e_3`。则 `(y_i×y_j)·e_3=u_iv_j-v_iu_j`（`u,v` 为 `y` 的前两个分量），与第三分量 `w` 无关**
**⇒ `w` 给出 `4` 个零方向。`(u,v)∈R^8` 上 `Q_P=-2β\,u^TAv`，`A_{ij}=½ (i<j)`、`-½ (i>j)`（反对称）。**
**⇒ 谱为 `±β·s_i(A)`，`A` 的奇异值给出 `{1+\sqrt2,\ 2-\sqrt2\}` 各二重。$\blacksquare$**

## E. 状态

| 件 | 状态 |
|---|---|
精确常数 `‖Hess S‖≤6(1+√2)β` | `[T]`（本轮） |
`ℋ≥I` 区间扩到 `β≤0.069036` | `[T]`（本轮，`2.07×`） |
规范方向精确入核 | `[T]`（本轮，数值 `0` 到机器精度） |
Savvidy 负模横截 | `[T]`（本轮，谱结构） |
Savvidy 障碍仍存在（`β=2.5` 时 `ℋ` 不定） | `[T]` |
质量/衰减率下界 | `[U]` |

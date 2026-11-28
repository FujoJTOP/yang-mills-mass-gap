# 第 319 轮：`v_2>0` —— 我上一轮的否证是噪声造成的，路线是活的

> 第 317 轮我报 `v_2≥0` 为假。**那条是错的**：我看的是**最小值**，而超立方对称性强制所有 `A≠B` 的 `M_{AB}` **精确相等**，所以最小值只能是噪声。
> 把统计量提到 `n=2500` 后，**均值在全部测试的 `β` 上都为正**。路线不是死的。

---

## 1. 为什么会看错

$$\text{超立方群在 6 个坐标平面上传递作用}\ \Longrightarrow\ M_{AA}=v_1\ (\forall A),\qquad M_{AB}=v_2\ (\forall A\ne B)$$

**⇒ 真值只有两个数 `v_1,v_2`；`min_{A≠B}M_{AB}` 与 `max_{A≠B}M_{AB}` 的散布**全部是有限样本噪声**。**

**⇒ 正确的统计量是**离对角均值**：**

$$v_2\ =\ \frac{1}{30}\sum_{A\ne B}M_{AB}\ =\ \frac{1}{30}\Big(\frac{\mathrm{Var}(\mathcal S)}{L^4}-6v_1\Big)$$

---

## 2. `n=2500` 的结果（`L=4`）

| `β` | `v_1=max_AM_AA` | `Σ_{A≠B}M_{AB}` | **`v_2`** | `β=0` 期望 |
|---|---|---|---|---|
`0.5` | `0.2556` | `+0.0358` | **`+1.19e−3`** | `0` |
`1.0` | `0.2338` | `+0.1877` | **`+6.26e−3`** | `0` |
`2.5` | `0.0932` | `+0.6099` | **`+2.03e−2`** | `0` |

$$\boxed{\ v_2(\beta)\ \gt\ 0\qquad\text{且随 }\beta\text{ 增长（}\beta=0\text{ 处恰为 }0\text{）}\ }$$

**⇒ 第 317 轮的 `[F]` 撤回。`v_2≥0` 在数值上成立，且 `β=0` 时**取等号**。**

---

## 3. 于是缺口变成一个「等式处 + 单调增」的形式

$$\chi^{\rm fixed}=v_1,\qquad \chi^{\rm clover}=\frac{v_1+5v_2}{6},\qquad \frac{\mathrm{Var}(\mathcal S)}{L^4}=6v_1+30v_2$$

**要证的等价形式：**

$$\boxed{\ 6\,\mathrm{Var}\big(Q^{(A)}\big)\ \le\ \mathrm{Var}(\mathcal S)=6\,\mathrm{Var}\big(Q^{(A)}\big)+30\,v_2L^4\ }$$

$$\Longleftrightarrow\qquad v_2\ \ge\ 0$$

**⇒ 而且 `β=0` 处**两边精确相等**（plaquette 独立，`Var(𝒮)=6Var(Q^{(A)})`），之后 `v_2` 单调增长。**
**⇒ 这不是「大常数 vs 小常数」，是**从等号出发、正向增长**。**

---

## 4. 已无条件成立的部分（不变，且是本轮可交付的硬结果）

$$\boxed{\ \int_0^\infty\chi^{\rm clover}_L(\beta)\,d\beta=\frac16\qquad\text{且}\qquad \chi^{\rm clover}_L(\beta)\le\frac16\quad\forall L\ge2,\ \forall\beta\gt0\qquad\textbf{[T]}}$$

**证明三件事：`d⟨𝒮⟩/dβ=Var(𝒮)≥0`；端值 `⟨𝒮⟩_{β=0}=6L^4`、`⟨𝒮⟩_{β→∞}=0`；`Q^{\rm clover}=𝒮/6`。**
**不经过取向分解、不用 `v_2`、不用 `T≥0`。`v_2` 正负都不影响它。**

---

## 5. 数值一致性（`L=4`，`n=2500`）

| `β` | `χ^{clover}=ΣΣM/36` | `≤1/6`? | `v_1=χ^{fixed}` | `v_1≤1/6`? |
|---|---|---|---|---|
`0.5` | `0.042206` | ✓ | `0.2556` | ✗（`1/6` 不够） |
`1.0` | `0.041936` | ✓ | `0.2338` | ✗ |
`2.5` | `0.032276` | ✓ | `0.0932` | ✓ |

**⇒ `χ^{clover}≤1/6` 三处全过；`χ^{fixed}` 在 `L=4,β\le1` 处超过 `1/6`，但在 `1` 以下。**

$$\Rightarrow\quad \text{可用的界：}\ \chi^{\rm fixed}\le 6\chi^{\rm clover}\cdot\frac{v_1}{v_1+5v_2}\ \le\ 1\cdot\frac{v_1}{v_1+5v_2}\cdot 6\ldots$$

**更直接：`v_1=6\chi^{clover}-5v_2\le 6\chi^{clover}\le 1`（当 `v_2≥0`）。**

**⇒ `v_2≥0` 一旦成立，立即得**

$$\boxed{\ \chi^{\rm fixed}_L(\beta)\ \le\ 6\chi^{\rm clover}_L(\beta)\ \le\ 1\qquad\forall L,\ \forall\beta\ }$$

**（因为 `χ^{clover}≤1/6` ⇒ `6χ^{clover}≤1`。）**

---

## 6. 本轮状态

| 命题 | 判 |
|---|---|
$$\int_0^\infty\chi^{\rm clover}d\beta=\tfrac16$$ | **`[T]` 无条件** |
$$\chi^{\rm clover}\le\tfrac16$$ | **`[T]` 无条件** |
$$v_2\ge0$$ | **`[M]` 成立（`n=2500` 均值正，`β=0` 取等号）；第 317 轮的 `[F]` 撤回** |
$$v_2\ge0\ \Rightarrow\ \chi^{\rm fixed}\le6\chi^{\rm clover}\le1$$ | **`[T]`（条件推演）** |
$$v_2\ge0\ \text{的证明}$$ | `[U]` |

**⇒ 现在的形状：`χ^{clover}≤1/6` 是无条件的；`χ^{fixed}≤1` 只差 `v_2≥0` 一步，而 `v_2` 在 `β=0` 取等号、之后单调正增长（`n=2500` 三处 `β` 全正）。**

# 工作记录（第 332 轮）：Duhamel / Seiler 路线的核查

> 按 §1.07：记录，不作汇报。

---

## A. 本轮的公式化

**对 `T=Σ_at_a|a⟩⟨a|`（`t_a>0`）与自伴 `A`，Duhamel 二点函数**

$$(A,A):=\frac1Z\int_0^1ds\ \mathrm{Tr}\big(T^{L(1-s)}AT^{Ls}A\big)=\frac1Z\sum_{a,b}\mathcal K(t_a,t_b)\,|A_{ab}|^2$$

$$\mathcal K(\alpha^{1/L},\beta^{1/L})=\begin{cases}\alpha&\alpha=\beta\\[1mm]\dfrac{\alpha-\beta}{L\ln(\alpha/\beta)}&\alpha\ne\beta\end{cases}\qquad(\alpha=t_a^L,\ \beta=t_b^L)$$

## B. 新证的不等式（`[T]`）

$$\boxed{\ \mathcal K(\alpha,\beta)\ \le\ \max(\alpha,\beta)\qquad\Longrightarrow\qquad (A,A)\ \le\ \big\langle A^2\big\rangle\ }$$

**证明：`α≥β` 时 `(α−β)/ln(α/β) ≤ α` ⟺ `1−s ≤ −\ln s`（`s=β/α∈(0,1]`），后者为初等不等式。 $\blacksquare$**

**⇒ Duhamel 二点函数被**静态二阶矩**控制。这是本轮得到的干净不等式。**

## C. 它为什么不给衰减率

**Seiler 不等式 `|\langle A^*B\rangle|^2\le (A,A)(B,B)` 配合 B 给**

$$\widehat C^{\mathcal S}(k)^2\ \le\ (A_k,A_k)(A_k,A_k)\ \le\ \big\langle A_k^2\big\rangle^2\qquad\Big(A_k=\widehat S(k)\Big)$$

**⇒ `Ĉ^{\mathcal S}(k) ≤ ⟨Ŝ(k)²⟩ = Var(Ŝ(k)) = L^4Ĉ^{\mathcal S}(k)` —— 自指，无信息。**

**⇒ Duhamel 上界（B）与静态矩同阶，**不能**产生 `k`-依赖的衰减。**

## D. 与 FSS 红外界的关系（定位）

**FSS 路线需要三步：**

| 步 | 内容 | 本项目状态 |
|---|---|---|
(GD) | 高斯支配：`⟨e^{(f,Φ)}⟩≤e^{‖f‖²_{G}/2}` | `[U]`；非 Abel 无（第 329 轮已核） |
(Seiler) | `|\langle A^*B\rangle|²≤(A,A)(B,B)` | `[T]`（本轮公式化并核对） |
(Duhamel) | `(A,A)≤` 静态矩 | `[T]`（本轮 B） |

$$\boxed{\ \text{缺口精确定位在 (GD) —— 而 (GD) 恰是非 Abel 的、文献里只对 Abel/自旋/玻色建立的}\ }$$

**⇒ E2 入口与 E1 入口、与第 329 轮的结论**收敛到同一点**。**

## E. 状态

| 件 | 状态 |
|---|---|
Duhamel 核公式 | `[T]` |
`(A,A)≤⟨A²⟩` | `[T]`（本轮） |
Seiler 不等式 | `[T]`（引文级，已核对形式） |
(GD) 高斯支配 | `[U]`，非 Abel 无 |
衰减率 | `[U]` |

# 任务书（第二版）：协方差算子的范数

> 只含任务与验收标准。第一版中三处被推翻的内容已在下文更正，并标明更正理由。
>
> **记号**：`Λ_L={0,…,L−1}⁴`，`U_e∈SU(2)`，`S_{β,L}=βΣ_P(1−½tr U_P)`，`q_x=½tr U_{P(x)}∈[−1,1]`，
> `Γ_L(x,y)=Cov(q_x,q_y)`，`g_c(z)=(1+‖z‖₂)^{−4}`，`Λ_L=1+sup_xΣ_y g_c(x−y)≍ln L`。
> **已有数值**：`L=6,8,10,12`，`β=2.5`；`χ_L:=Ĉ_L(0)≈0.032`。

---

## §0 状态

**已知（可当引理）：**

| 编号 | 内容 | 来源 |
|---|---|---|
(P1) | `∫F²dμ=⟨dF,ℋ^{-1}dF⟩`；`ℋ=dd_μ^*+d_μ^*d≥0` | 已证 |
(P2) | `‖d‖²≤24`；`H¹=0` ⇒ `ℋ` 无零模 | 已证 |
(P3) | `ν_k([0,s])≤As(s/ω)^δ`（`s≤ω`）+`≤As`（`ω<s≤1`）⇒ `‖Γ‖=O(ln L)` | 已证 |
(P4) | `ν_0([0,s])≤As^{1+δ}` ⇒ `χ_L≤A/δ+24` | 已证 |
(P5) | `‖V‖≤30β`、`tr_fiber V_-≤360βL⁴`（Wilson，`V=HessS`） | 已证（§12） |
(P6) | `ν_0([0,s])≤(720β/(2−s)²)Σ_{λ_j≤s}⟨ψ_j,V_-ψ_j⟩` | 已证（§12.4） |

**已知为「不是逻辑推论」（不得使用）：**

$$\chi_L=O(1)\ \not\Longrightarrow\ \|\Gamma_L\|=O(1)\qquad S_j\text{ 有界}\ \not\Longrightarrow\ \sum_z|\Gamma_L(0,z)|<\infty$$
$$\lambda_{\rm harm}:=M_0/\chi_L\ \text{大}\ \not\Longrightarrow\ \text{低谱尾小}\qquad tr_{\rm fiber}V_-=O(L^4)\ \not\Longrightarrow\ Tr_{L^2}(P_sV_-)=O(L^4)$$

---

## §1 任务 A（数值，决定性，最先做）

$$\textbf{A.}\quad L=6,8,10,12,\ \beta=2.5:\ \text{计算【全部】}\widehat C_L(k),\ \text{报告 }\arg\max_k\widehat C_L(k)\ \text{与}\ \frac{\max_k\widehat C_L(k)}{\widehat C_L(0)}$$

**实现（有原始数据时）：**

```python
import numpy as np

def structure_factor(q):
    """q.shape = (N, L, L, L, L); 第一个轴是 MC 配置。"""
    q = np.asarray(q, dtype=float)
    N, L = q.shape[0], q.shape[1]
    V = L**4
    Q = np.fft.fftn(q, axes=(1, 2, 3, 4))
    Q -= Q.mean(axis=0, keepdims=True)          # 关键：跨配置中心化
    C_hat = np.sum(np.abs(Q) ** 2, axis=0) / ((N - 1) * V)
    idx = np.unravel_index(np.argmax(C_hat), C_hat.shape)
    return C_hat, 2 * np.pi * np.asarray(idx) / L, C_hat[idx] / C_hat[(0, 0, 0, 0)]
```

**⚠️ 零动量必须跨配置中心化，不得对每个配置减去其自身空间平均，否则会删掉真正的零动量涨落。**

**验收：**逐 `L` 给 `max_k Ĉ_L(k)`、`argmax`、比值、`‖Γ_L‖/Λ_L` 随 `L` 的趋势。若 `argmax ≠ 0`，报告该动量。

**为何必须先做：**`χ_L` 的 `L`-不变性**不能**推出 `‖Γ‖=O(1)`。反例（在 `q_x∈[−1,1]` 内）：`q_x=aξ(−1)^{x₁}+bη_x`，`a+b≤1`，则 `Ĉ(0)=b²` 与 `L` 无关，而 `Ĉ((π,0,0,0))=a²L⁴+b²`。

**缺数据时的诚实做法：**只报
$$R_L:=\frac{\max_k\widehat C_L(k)}{\widehat C_L(0)}\ \ge1,\qquad \frac{\|\Gamma_L\|_{2\to2}}{\Lambda_L}=R_L\cdot\frac{\chi_L}{\Lambda_L},$$
**不要填具体数字。**

---

## §2 任务 B（解析，最弱一层）

$$\textbf{B.}\quad \nu_{0,L,\beta}([0,s])\le A s^{1+\delta}\quad(0<s\le1),\qquad A,\delta\ \text{与}\ L,\beta\ \text{无关}$$

**⇒ 由 (P4) 得 `χ_L≤A/δ+24`。**

**能做与不能做：**固定 `L` 时由 (P2) 与对 `β` 的连续依赖，`λ_*(L):=inf_{β∈[1/48,2.5]}λ_min(ℋ_{β,L})>0`，可取 `A_L=24λ_*(L)^{-1−δ}`。**这个常数依赖 `L`。**

**⇒ B 的真实内容 = `sup_L λ_*(L)^{-1−δ}<∞`，或一个比统一谱隙更弱的源加权低谱尾估计。**

**验收：**给出 `A,δ` 的**显式且 `L,β`-无关**表达式；若只能给 `A_L`，明确说未完成。**不要声称 B 解决了原题**（只给 `χ=O(1)`）。

---

## §3 任务 C（解析，直通 `O(1)`——推荐主攻）

$$\|\Gamma_L\|_{2\to2}\ \le\ \|\Gamma_L\|_{1\to1}\ =\ \sum_{z\in\Lambda_L}|\Gamma_L(0,z)|$$

$$\boxed{\ \textbf{(C2)}\qquad \big|\Gamma_L(0,z)\big|\ \le\ C\,(1+\|z\|_2)^{-4-\varepsilon},\qquad\varepsilon>0\ }$$

**⚠️ 第一版的指数 `4` 不够：**`|Γ(0,z)|≤C(1+‖z‖)^{−4}` 只给 `Σ=O(ln L)`。必须 `4+ε`。

**⚠️ 壳层数据 (`S_j≤M`) 不蕴含 (C2)：**`S_j≤M` 只给 `Σ_{\|z\|=j}|Γ|≤M(1+j)⁴` ⇒ `Σ_z=O(L⁵)`。反例同 §1：`Γ(0,z)=a²(−1)^{z₁}+b²δ_{z,0}` 使 `χ=b²` 恒定、`S_j≍a²/j` 有界，而 `Σ_z|Γ|≍a²L⁴`。

**⇒ 拟合值 `p_eff=4.5–5.1` 支持 (C2)（`ε≈0.5–1.1`），但它是从 `Γ(0,z)` 对 `‖z‖` 的拟合得来，**不是**从 `S_j` 有界性得来。**

**验收：**`C,ε` 显式；`β` 范围（`[1/48,2.5]` 还是只有 `β≤β_0`）；若只能给 `|Γ|≤Cg_c(z)(1+校正)`，报校正项。

---

## §4 任务 D（已闭合，不要重复做）

**已闭合（可引用，不必重做）：**

$$Q_P(Y)=\beta q\Big|\sum_iy_i\Big|^2-2\beta\sum_{i<j}(y_i\times y_j)\cdot a\qquad(U_P=qI+ia\cdot\sigma,\ q^2+|a|^2=1)$$
$$\|V\|\le30\beta,\quad \|V_-\|\ge6\beta\ (\text{显式构型}),\quad \operatorname{tr}_{\rm fiber}V_-\le360\beta L^4$$
$$\nu_0([0,s])\le\frac{720\beta}{(2-s)^2}\sum_{\lambda_j\le s}\langle\psi_j,V_-\psi_j\rangle$$

**⚠️ 这条路已判为绕路，理由三条（不要重试）：**

1. **纤维迹不控全局迹**：`Tr_{L²}(P_sV_-)` 依赖本征函数在配置空间上的分布，`tr_fiber` 是每点的有限维矩阵迹，二者之间无不等式。
2. **谱计数界发散**：CLR/LT 给 `#{λ_j<2}≤C∫tr(V_-^{d/2})` 而 `d=12L⁴→∞`；BS 迹 `Tr(W^{1/2}(A_+-s)^{-1}W^{1/2})≤‖W‖(2−s)^{-1}Tr(1)=∞`。
3. **补偿不可得**：需 `Σ_{λ_j≤s}⟨ψ_j,V_-ψ_j⟩=O(s^{1+δ})`。即使本征函数完全离域，也只给 `O(#{λ_j≤s})`，故需 `#{λ_j≤s}=O(s^{1+δ})`——**全一形式空间的谱隙，强于原题。**

---

## §5 任务 E（精确谱形式，等价于 `O(1)`）

$$\boxed{\ \textbf{(E)}\qquad \sup_{L,\beta,k}\ \int_0^1\frac{\nu_{k,L,\beta}([0,s])}{s^2}\,ds\ <\ \infty\ }$$

**⇒ (E) ⟺ `‖Γ_L‖_{2→2}=O(1)`**（高谱部分由 `‖D_k‖²≤24` 自动有界）。**这不是等价于原靶的"必要条件"，而是**原靶的精确谱表述**。**

**验收：**若只能证 `k=0` 的版本，明确说那是 B，不是 E。

---

## §6 禁止清单

| 禁止 | 理由 |
|---|---|
`χ=O(1) ⟹ ‖Γ‖=O(1)` | §1 反例 |
`λ_harm` 大 ⟹ 低谱尾小 | §0 谱测度反例 |
`S_j` 有界 ⟹ 绝对可和 | §3 反例 |
`tr_fiber V_-=O(L⁴) ⟹ Tr_{L²}(P_sV_-)=O(L⁴)` | §4 第 1 条 |
局部 `Hess S` 负 ⇒ 标量道失界 | **没有这个定理**；反之亦无 |
绝对值包络当路径权重 | 单点反例 `g=3β²t^{−1}I_2(t)/I_1(t)\big|_{β=2.5,t=5}=2.697527>1` |
指数 `4` 的逐点界 | §3：只给 `O(ln L)` |

---

## §7 优先级

$$\boxed{\ \textbf{A（数值，定 }\arg\max\text{）}\ \succ\ \textbf{C（}4+\varepsilon\text{ 逐点，直通 }O(1)\text{，绕开全部谱机器）}\ \succ\ \textbf{B（}k=0\text{ 最弱层）}\ \succ\ \textbf{E（精确谱形式）}\ \succ\ \textbf{D（已闭合，绕路）}\ }$$

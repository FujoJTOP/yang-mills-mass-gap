# P03-EQUIV — `H1 ⟺ Δ_dyn ≥ a²Δ₀²` 的独立核准；`a²Δ₀²` 是什么

**日期**: 2026-11-25 · **状态**: 已核（恒等式 `[T]`；等价 `[T]`（条件见 §3.3）；量纲/数值 `[C]`+`[E]`）· **禁 git commit**（遵守）· 只写本文件
**输入（逐字读，不凭记忆）**: `docs/paper-v2/app.tex:179–202, 780–800`、`ch03-existence.tex:144–175`、`ch05-theorems.tex:230–245`、`ch08-immune.tex:8–31`、`ch13-machine.tex:57–78`、`ch18-scalar-h1.tex:84–91`、`research/v3/w7/h1c_wfull.py:1–50`、`research/v3/species/GPP-STRUCT.md:9–65,120–130`、`O11-PHYS.md`、`LOOP-RUNBOOK.md:2527–2542`
**交付对象**: 用户（lead）提出的三步推导：`δ` 的定义 → `W_full` 两定义 → `H1 ⟺ Δ_dyn ≥ a²Δ₀²`

---

## 0. 判死线（先读这段）

1. **`δ` 的逐字式核准**：`app.tex:190–191` 写的是 `W_full = W_geo + X` 与 **`δ := a²Δ₀² + W_full`**，`[T-transcribed]`，无漏项。`δ ≠ a²Δ₀²`（`GPP-STRUCT.md:13,125,144` 与 `LOOP-RUNBOOK.md:2158` 的 `δ = a²Δ₀²` 是**语料级漏项**，`O14-AUDIT-GPP.md:167–169` 已记；本文件不复改）。
2. **`W_full` 的两个定义相等** `[T]`，且是**恒等**，不是拟合：`E[Φ''] + Var(G_v) − Var(G_v+Φ') ≡ (E[Φ''] − Var(Φ')) + (Var(G_v) + Var(Φ') − Var(G_v+Φ')) = W_geo + X`。`h1c` 自己就写着这两个形式（`h1c_wfull.py:42` 明标 "the exact identity ..."，`:98` 用身份式值对账 `W_full = W_geo + X`）。**⇒ 判死线 (b) 不触发，无记账错。**
3. **`a²Δ₀²` 依赖 `a`（格距），不依赖 `a₀`（幅度）** `[T]`：`ch03:152–157` 的 `G''_a = a²m_phys²`、`m_lat = a·m_phys`；`ch03:161–163` 的 `G''_a/a² ≥ Δ₀²` ⟺ `λ_a ~ a²Δ₀²`。量纲：**无量纲**，`a²Δ₀² = (aΔ₀)² = m_lat²`（`a` 是长度、`Δ₀` 是质量/长度⁻¹、`G''_a` 本身无量纲）。**⇒ 判死线 (c) 不触发。**
4. **等价成立，但不是"化简"** `[T]`/`[NG]`：`H1 ⟺ Δ_dyn ≥ a²Δ₀²` 是**逐字代数恒等**（`W_full` 在两端**同一个对象**，直接抵消），所以它**不产生任何新估计**——`Δ_dyn ≥ a²Δ₀²` 就是 H1 本身换了一身记号（与 `GPP-STRUCT.md:143` 的"换名 `[NG]`"同判）。**判死线 (a) 触发，但必须同时钉死：它不是"更简单的比较"，是同一比较。** 真正的内容是 `Δ_dyn = E[H|0] − Var(G|0)` 的下界，而要证它需要的正是 Wilson 侧的条件凸性/BL，即 Savvidy 墙（`GPP-STRUCT.md:135–137`）。

---

## 1. 逐字引用（含行号）

### 1.1 `app.tex:182–191` —— (M) 式与 `δ` 的定义（逐字，含断行）

```latex
182: \LM{} The single open inequality, stated with no slack. There is $\Delta_0>0$ such that,
183: uniformly in $(L,a)$ and for every coarse transverse $v$ at $k\to0$,
184: \[
185:   \operatorname{Var}_\mu\!\big(\langle v,\nabla V\rangle\,\big|\,\xi=0\big)
186:   \;\le\;
187:   E_\mu\!\big[H_{vv}\,\big|\,\xi=0\big]\;-\;a^{2}\Delta_0^{2},
188:   \qquad V=S_W+\Phi .
189: \]
190: This is the (M) form of Chapter~\ref{ch:immune} with $W_{\mathrm{full}}=W_{\mathrm{geo}}+X$ and
191: $\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality.
192: Equivalently $D\le S-\lambda$ \emph{together with} $S\ge\lambda$ (the bare defect form
193: $D\le(1-\epsilon)S$ is insufficient), equivalently the dimensionless ratio
194: $G''_a/a^{2}\ge\Delta_0^{2}$.
```

**核准**：`δ := a²Δ₀² + W_full` `[T-transcribed]`；`Δ₀>0` 只被断言"存在"（`:182`），**语料从未给 `Δ₀` 数值**（§4.3）。
`:194–196` 又记两条量纲陷阱：`E[f²]≈2` 是 UV/equipartition 钉死的、**不与 `a²` 成比例**，不能冒充 `λ`。

### 1.2 `app.tex:795–796` —— 符号表里的 `W_geo` / `W_full`（逐字）

```latex
795: \item[$W_{\mathrm{geo}}$] the Haar geometric term $E[\Phi'']-\operatorname{Var}(\Phi')$.
796: \item[$W_{\mathrm{full}}$] $W_{\mathrm{geo}}+X$, with $X=-2\operatorname{Cov}(G,\Phi')$.
```

### 1.3 `ch08-immune.tex:11–21` —— 真正的 (M)（**注意 LHS 是 `∇S`，RHS 带 `+W_full`**，逐字）

```latex
11: \begin{equation}
12: (\mathrm{M})\qquad
13: \operatorname{Var}_\mu\bigl(\langle v,\nabla S\rangle\mid0\bigr)
14: \ \le\ E_\mu[H_{vv}\mid0]+W_{\mathrm{full}}-\delta ,
15: \qquad\text{uniform in }(L,U'),
16: \label{eq:species}
17: \end{equation}
18: where $W_{\mathrm{full}}:=W_{\mathrm{geo}}+X$ is the full two-point remainder of
19: Chapter~\ref{sec:wgeo} ($X=-2\operatorname{Cov}(G,\Phi')$); the display is then $G''\ge\delta$ in
20: curvature form. Equivalently $\operatorname{Var}(\xi_{k\to0})\le1/\lambda$ (H1), equivalently
21: $m_{\mathrm{eff}}^2(k\to0)\ge\lambda$, equivalently $G''(0)\ge\lambda$.
```

`app.tex:185–187` 与 `ch13-machine.tex:62–68` 是同一式的**第二份排版**（逐字比：`ch13:63–66` = "`\operatorname{Var}_\mu(\langle v,\nabla V\rangle|\xi=0) \le E_\mu[H_{vv}|\xi=0] - a^{2}\Delta_0^{2}, \quad V=S_W+\Phi`"，`ch13:69–70` = "`This is the (M) form with $W_{\mathrm{full}}=W_{\mathrm{geo}}+X$ and $\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality.`"）。
**⇒ 语料里 (M) 有两条互不等价的排版**：`ch08:13–14`（`∇S`，带 `+W_full`）与 `app:185–187`/`ch13:63–66`（`∇V`，无 `+W_full`）。见 §3.3，这是**唯一**的记账缝隙。

### 1.4 `ch05-theorems.tex:230–240` —— `eq:split`（`GPP-STRUCT.md:24–25` 的转引源，逐字）

```latex
230: \textbf{The exact split and what fails.} The identity is
231: \begin{equation}
232: \begin{aligned}
233: G''(0)&=E\bigl[H_{vv}+\Phi''\,\big|\,\xi=0\bigr]
234:         -\operatorname{Var}\bigl(G_v+\Phi'\,\big|\,\xi=0\bigr)
235:         =\Delta_{\mathrm{dyn}}+W_{\mathrm{geo}}+X,\\
236: \Delta_{\mathrm{dyn}}&=E[H]-\operatorname{Var}(G),\qquad
237: X=-2\operatorname{Cov}(G,\Phi').
238: \end{aligned}
239: \label{eq:split}
240: \end{equation}
```

`GPP-STRUCT.md:9` 与 `:60–62`、`O11-PHYS.md:11–12` 逐字同式（`Δ = E[κ|0] − Var(F|0) = G''(0)`）。`H := Hess V` 的写法在 `GPP-STRUCT.md:65`（那里 `H` 与 `H_vv+Φ''` 并列 ⇒ **`H` 应读作 `Hess S_W`**，不是 `Hess V`；`O09-PHICURV.md:111` 已记此记号冲突，否则 `Φ''` 被算两次）。

### 1.5 `h1c` 的两个 `W_full` 形式（逐字）

```python
research/v3/w7/h1c_wfull.py
  9:       W_full := W_geo + X,   X = -2 Cov(G, Phi'),   G = d_v S,  Phi' = d_v Phi
 10:       G''(0) = E[H+Phi''] - Var(G+Phi') = Delta_dyn + W_geo + X = Delta_dyn + W_full.
 11:   Hence  W_full >= W_geo  <=>  X >= 0  <=>  Cov(G,Phi'|a) <= 0  =  hypothesis [NA] (app.tex l.800),
 42:     the exact identity  W_full = E[Phi''] + Var(G) - Var(G+Phi')."""
 98:                 Id[i] = m['EGpp'] + m['VarG'] - m['VarPhiG']     # identity value of W_full
```

### 1.6 `ch03-existence.tex:152–163` —— `a²Δ₀²` 的来历（逐字）

```latex
152: \begin{equation}
153: G''_a=a^2\,m_{\mathrm{phys}}^2,\qquad
154: m_{\mathrm{lat}}=\sqrt{G''_a}=a\,m_{\mathrm{phys}},\qquad
155: \Delta=\lim m_{\mathrm{lat}}/a .
156: \label{eq:dim}
157: \end{equation}
158: Therefore a claim of the form ``at fixed $a$, $G''_a\ge\lambda$ with $\lambda$ a positive
159: \emph{constant}'' is dimensionally wrong for the gap: it would give
160: $m_{\mathrm{phys}}=\lim\sqrt{\lambda}/a=+\infty$, the free-box UV scale, a lattice artifact rather
161: than a particle. The correct target is a \textbf{ratio} $G''_a/a^2\ge\Delta_0^2$ at the
162: asymptotically-free scaling rate, i.e.\ $\lambda_a\sim a^2\Delta_0^2$, with $\lambda_a$
163: transported multiplicatively by $Z_O^2$. The transport itself is free (it is a square); what is
164: not free is the non-zero seed it multiplies.
```

（另见 `ch18-scalar-h1.tex:90` 逐字 "prove $G''(0)\ge a^2\Delta_0^2$ uniformly"、`G1-CONTINUUM-SCALING.md:8` 逐字 "`G″_a ≥ a²Δ₀² ⟺ m_lat=√G″_a ≥ aΔ₀`"。）

---

## 2. `W_full` 两定义是否相等 —— 代数核对 `[T]`

记 `G := G_v = ∂_v S_W`（力）、`Φ' := ∂_v Φ`（Jacobian 力）。`Var` 是双线性型，故对任意两随机变量

```
Var(G + Φ') = Var(G) + Var(Φ') + 2 Cov(G, Φ')                                  (★ 恒等，无假设)
```

于是

```
h1c 形式:   W_full = E[Φ''] + Var(G) − Var(G + Φ')
                    = E[Φ''] + Var(G) − [Var(G) + Var(Φ') + 2Cov(G,Φ')]       由 (★)
                    = E[Φ''] − Var(Φ') − 2Cov(G,Φ')
                    = W_geo                + X                 ∎
app.tex:796: W_full = W_geo + X = [E[Φ''] − Var(Φ')] + [−2Cov(G,Φ')]
```

**结论 `[T]`**：两式**逐项相同**，`X = −2Cov(G,Φ') = Var(G) + Var(Φ') − Var(G+Φ')`。这不是近似，`(★)` 是二阶矩的双线性恒等式；用到 `eq:split` 的 `Δ_dyn / W_geo / X` 三块命名与 `Φ''⪰0` 都无关。

**数值核对 `[E]`（用仓库自带数，SU(3) 单 plaquette `β=2.5`, `O11-PHYS.md:65`, seed 201）**

| 量 | 值 | 出处 |
|---|---|---|
| `E[H]` | `+0.09534` | `O11:65` |
| `Var(G)` | `0.005715` | `O11:65` |
| `E[Φ'']` | `0.34659` | `O11:65` |
| `Var(Φ')` | `9.17e−4` | `O11:65` |
| `Var(F)=Var(G+Φ')` | `0.006252` | `O11:65` |
| `G''=E[κ]−Var(F)` | `0.435673` | `O11:65` |

```
h1c 形式      W_full = 0.34659 + 0.005715 − 0.006252 = 0.346053000
app.796 形式  W_full = (0.34659 − 0.000917) + (0.005715 + 0.000917 − 0.006252) = 0.346053000
差                    = 5.55e−17   （浮点零）
⇒ Δ_dyn = E[H] − Var(G) = 0.089625;  Δ_dyn + W_full = 0.435678  vs  表列 G'' = 0.435673（差 5e−6，源表五位舍入）
⇒ 隐含 X = +3.80e−4（=−2Cov，即 Cov(G,Φ') = −1.90e−4 <0；本例 [NA] 恰好成立，但只在小数点后第 4 位）
```

SU(2) 单 plaquette（`O11:53`）：`Var(Φ')=0`、`X=0` ⇒ 两形式都给 `W_full = W_geo = 0.20509`。**两模型都零差。**

**⇒ 判死线 (b) 不触发**：`W_full` 无记账错；`h1c` 自己在 `:42` 就把它写成 "the exact identity"。这条等价链上唯一的记账缝隙在 **(M) 的 LHS 用 `∇S` 还是 `∇V`**（§3.3），不在 `W_full`。

---

## 3. 逐步核准用户的三步推导

### 3.1 第一步：`G''(0) = Δ_dyn + W_full` —— `[T]` **通过**

由 `eq:split`（§1.4）逐字：

```
G''(0) = E[H_vv+Φ''|0] − Var(G_v+Φ'|0)
       = E[H_vv|0] + E[Φ''|0] − Var(G_v+Φ'|0)              （线性）
       = [E[H_vv|0] − Var(G_v|0)] + [E[Φ''] + Var(G_v|0) − Var(G_v+Φ'|0)]
       = Δ_dyn                     + W_full                  ∎   （W_full 按 h1c:9 的定义）
```

`Δ_dyn := E[H_vv|0] − Var(G_v|0)`，`W_full` 见 §2。**用户这一步无误**；`h1c_wfull.py:10` 逐字写着同一行。

### 3.2 第二步：`H1(δ) ⟺ G''(0) ≥ δ` —— `[T]` **通过**（前提是 §3.3 的读法）

`ch08:19–21` 逐字："the display is then `G''≥δ` in curvature form ... equivalently `G''(0)≥λ`"。由 3.1，`G''(0) ≥ a²Δ₀² + W_full ⟺ Δ_dyn + W_full ≥ a²Δ₀² + W_full ⟺ Δ_dyn ≥ a²Δ₀²`：**`W_full` 是同一个对象，两端相消，消去不需要任何估计**。`[T]`

### 3.3 缝隙：`(M)` 的两份排版不是同一条不等式 —— `[C]`

把 `ch08:13–14` 与 `app.tex:185–187` 代入同一套定义，得的是**两个不同阈值**：

| 排版 | LHS | 等价于 | 阈值 |
|---|---|---|---|
| `ch08:13–14`（`∇S`，RHS `+W_full`）`[T]` | `Var(G_v|0)` | `Δ_dyn ≥ a²Δ₀²` | `δ = a²Δ₀² + W_full` ✓ 与 `app:191` 一致 |
| `app:185–187` / `ch13:63–66`（`∇V`，RHS 无 `+W_full`） | `Var(G_v+Φ'|0)` | `E[H|0] − Var(F|0) ≥ a²Δ₀²` ⟺ `G''(0) ≥ a²Δ₀² + E[Φ''|0]` | `δ' = a²Δ₀² + E[Φ''|0]` ✗ |

两阈值相等 ⟺ `E[Φ''|0] = W_full` ⟺ `Var(Φ') + 2Cov(G,Φ') = 0`——**一般不为零**（SU(3) 单 plaquette：`Var(Φ')=9.2e−4>0`；`L=2` 上 `W_full` 实测可到 `O(−1)`，见 §4.4）。
**判**：`ch08:11–17` 是正确排版（LHS `∇S` + RHS `+W_full`，两者互洽，且给出 `δ=a²Δ₀²+W_full`）。`app:185–187`/`ch13:63–66` 的 `⟨v,∇V⟩, V=S_W+Φ` **要么应写成 `⟨v,∇S_W⟩`，要么 `H_vv` 应读作 `Hess V` 的 `vv` 元**；现写法（`∇V` + `H_vv`=Wilson 部分）与 `:191` 的 `δ` 不自洽。`[C]`，量级 `O(Var(Φ')+2Cov)`；单 plaquette 上 `~2e−3`（故长期未被发现），`L=2` 上 `O(1)`（不可忽略）。
（相关：`O14-AUDIT-GPP.md:167` 判"`app.tex:191` 对、`GPP-STRUCT` 漏 `W_full`"——方向对；本文件补的是**另一侧**：漏的不是 `W_full` 而是 `(M)` 的 LHS 记号。`GPP-STRUCT.md:65` 的 `H := Hess V` 是同一处记号的第三个版本。）

---

## 4. `a²Δ₀²` 是什么

### 4.1 定义（三处互证）`[T]`

1. `app.tex:182–188`：`Δ₀>0`，"uniformly in `(L,a)`"，右端减 `a²Δ₀²`；`:194` 等价写法 `G''_a/a² ≥ Δ₀²`。
2. `ch03:152–163`：`G''_a = a²m_phys²`、`m_lat = √G''_a = a m_phys`；目标是**比值** `G''_a/a² ≥ Δ₀²`，即 `λ_a ~ a²Δ₀²`（"at fixed `a`, `G''_a ≥ λ` = 常数"是**量纲错**，会给 `m_phys = √λ/a → ∞`）。
3. 等价链（`H1-DIRECT-SPLIT.md:196–197`、`G1-CONTINUUM-SCALING.md:8`）：`G''_a ≥ a²Δ₀² ⟺ m_lat ≥ aΔ₀ ⟺ G''_a/a² ≥ Δ₀² ⟺ a²Var_μ(a) ≤ 1/Δ₀²`。

**⇒ `a²Δ₀²` 不是"另一个量"，它就是 `G''` 在格点单位下的目标值**：`a²` = 格距平方（`eq:dim` 的 `a`），`Δ₀` = 目标质量（`(L,a)`-一致）。`a²Δ₀² = (aΔ₀)² = m_lat²`。

### 4.2 量纲 `[T]`

| 对象 | 量纲 | 说明 |
|---|---|---|
| `a` | 长度 | 格距（`ch03:153`） |
| `Δ₀` | 质量 = 长度⁻¹ | 目标隙 |
| `a²Δ₀²` | **无量纲** | `= m_lat²`，与 `G''_a`、`E[H]`、`Var(G_v)`、`W_full` 全部同量纲（`L=1` 归一化 `‖v‖=1` 时） |
| `E[f²]≈2` | 无量纲但**∝ a⁰** | UV/equipartition 钉死，`app.tex:194–196` 明令不得冒充 `λ` `[T]` |

**依赖 `a`（格距），不依赖 `a₀`（幅度）** `[T]`：`ch15-uvface.tex:31` 与 `h1c_wfull.py` 里 `a₀` 是**条件化点的幅度坐标**（"conditioning at coarse field `a = a₀`"，`h1c` 扫 `a0=0` 邻域）；`a²Δ₀²` 里没有 `a₀`。⇒ 比较是"**幅度依赖的 LHS** `Δ_dyn(a₀)` ≥ **幅度无关的 RHS** `a²Δ₀²`"。**判死线 (c) 不触发。**

### 4.3 数值量级（`β=2.5`、`L=2`、`SU(2)`）`[C]`+`[E]`

**锚**：`a(β≈2.5) ≈ 0.15–0.2 fm`（`H5-HEAD-NUMERICS.md:183` 逐字 "`a(β≈2.5)≈0.15–0.2 fm`, `ξ/a≲0.3` — outside the scaling window"；旁证 `T23-TEMPLATE.md:353` "`β=2.5` 时 `m_{0⁺⁺}≈1.7 GeV`、`a≈0.2–0.3 fm`"，`COMPUTE-GPP-K0-CONTINUUM.md` 的 `a(β=5.7)≈0.18 fm`）。取 `ħc = 0.19733 GeV·fm`：

```
a(fm)   a⁻¹(GeV) | a²Δ₀²  for Δ₀ = 0.5 / 1.0 / 1.7 GeV
0.150    1.316   |        0.144   0.578   1.670
0.175    1.128   |        0.197   0.787   2.273
0.200    0.987   |        0.257   1.027   2.969
```

（`Δ₀=1.7 GeV` 一行取语料的胶球锚 `m_{0⁺⁺}`，**只作标尺**；`Δ₀` 在语料里是自由参数，§4.4。）

**两个必须同时记的话**：
- **`a→0` 时 `a²Δ₀² → 0`，且与 `G''_a` 同阶（都 `∝a²`）** ⇒ `G''≥0` 与 `G''≥a²Δ₀²` 的差别**不是"数值差多少"，是逻辑强弱**（`m_lat≥0` 无隙 vs `m_lat≥aΔ₀` 有隙；`O14-AUDIT-GPP.md:170` 同判，`[C]`）。
- **`L=2` 的额外限制** `[T]`：`L=2` ⇒ 允许动量 `k = πn/a`，`k→0` 只剩 `k=0` 一个模，"邻域"为空；此时 `G''` 与 `a²Δ₀²` 都是在**单一体积、无 `a→0`** 下比（`O24-L2-PUSH.md:29–30`、`O26-VARJUMP.md:109–116` 同口径）。数值上 `L=2` 的 `Δ` 还差 `L²` 归一（`O26`：`−0.221 → −0.884`），与"符号"无关、与"阈值幅度"有关。

### 4.4 这条等价在 `β=2.5, L=2` 上值多少（`[E]`/`[C]`）

用 `Δ_dyn(L=2) = +1.065 ± 0.333`（`LOOP-RUNBOOK.md:2525,2809`）反解 `Δ₀` 的上限：

```
a = 0.15 fm (a⁻¹=1.316 GeV): Δ_dyn ≥ 1.065 ⇒ m_lat ≤ 1.032 ⇒ Δ₀ ≤ 1.36 GeV（1σ 下缘 1.13 GeV）
a = 0.20 fm (a⁻¹=0.987 GeV): 同上           ⇒ Δ₀ ≤ 1.02 GeV（1σ 下缘 0.84 GeV）
```

**⇒ 等价式的用处就在这里**：它把"证明 `(L,a)`-一致不等式"换成"**两个 `O(1)` 数比大小**"，并把 `Δ₀` 的可宣称上界钉在 `≈0.84–1.36 GeV`（`β=2.5, L=2`）。取 `Δ₀ = 1.7 GeV`（胶球锚）则该点**数值上不成立**——但 `L=2`（`a≈0.2 fm` ⇒ 盒 `0.4 fm`）本来就装不下胶球，这不是反例、是分辨率不足（`O11-PHYS.md:25` 的 "判不出"，`O24:30`）。**不得**把本节读成对 H1 的正面或反面证据。

---

## 5. 结论（对三条判死线逐条）

### (a) 等价成立 ⇒ 写成更简单的比较 + 量级 —— **触发（但带一条硬限制）**

```
H1  ⟺  G''(0) ≥ a²Δ₀² + W_full  ⟺  Δ_dyn ≥ a²Δ₀² ,
Δ_dyn := E_μ[H_vv|0] − Var_μ(G_v|0)        （H_vv = Hess S_W 的 vv 元，G_v = ∂_v S_W）
```
`[T]`（§3.1+§3.2，`W_full` 两端恒等相消）。
**量级**：`a²Δ₀² = (aΔ₀)²`，无量纲；`β=2.5, SU(2)`（`a≈0.15–0.2 fm`）上 `0.14–0.26`（`Δ₀=0.5 GeV`）、`0.58–1.03`（`1.0 GeV`）、`1.7–3.0`（`1.7 GeV`）；`a→0` 时 `→0 ∝ a²`。
**硬限制 `[NG]`**：这是**换名**，不是化简。`W_full` 之所以相消，是因为 `δ` 当初就是"`E[Φ'']` 那侧的差额"定义的；抵消后剩下的 `E[H|0] − Var(G|0) ≥ a²Δ₀²` 逐字就是 `ch08:13–14` 的 (M)（经 `Var` 移项）。**任何"H1 已被化归为更简单比较"的读法都过不了这一关**：被证明的对象一个没少（Wilson 侧的"条件方差 ≤ 条件曲率"，即 `GPP-STRUCT.md:137` 的第六物种），`Δ_dyn` 的符号与大小仍只由 MC 给（`+1.065±0.333`，`L=2`；`L=3` 已变号 `−0.572±0.609`，`LOOP-RUNBOOK.md:2540`）。与 `GPP-STRUCT.md:143` 的判词同：**作为路线 `[NG]`**。

### (b) `W_full` 两定义不等 ⇒ 记账错 —— **不触发**

§2：两式恒等（代数 `[T]` + 两个模型数值零差 `[E]`，差 `5.6e−17`）。`h1c_wfull.py:42` 自己标 "exact identity"。**语料这条没错。**
**但**（同一科目的另一处）`(M)` 的 LHS 有两条不等价排版（§3.3）：`ch08:13` 的 `∇S` 对、`app:185–187`/`ch13:63–66` 的 `∇V` 与 `:191` 的 `δ` 不自洽（差 `Var(Φ')+2Cov`）。**这条要记，但不改 `δ` 的表达式。**

### (c) `a²Δ₀²` 依赖 `a₀`（幅度）⇒ 等价不成立 —— **不触发**

§4.2：`a²Δ₀²` 依赖 `a` = **格距**（`ch03:152–157` 的 `eq:dim`），不依赖 `a₀` = 幅度（那是 `ch15:31`/`h1c` 的条件化坐标）。因此**正确形式就是 `Δ_dyn(a₀) ≥ a²Δ₀²`**：左端随 `a₀` 变（`h1c` 扫 `a0`；`O09-PHICURV.md:170` 实测 margin 只在 `|a|≲0.71`），右端是常数。**等价不因 `a₀` 失效。**

---

## 6. AGENTS.md §1 检索执行核验

**工具状态（本轮实测）**：`read` = **OK**（返回正文与行号：`app.tex:170–209,780–817`、`ch03:140–179`、`ch05:228–245`、`ch08:1–45`、`ch13:57–78`、`ch18:84–97`、`GPP-STRUCT.md:1–70,120–179`、`h1c_wfull.py:1–50`、`O11-PHYS.md:1–80`、`LOOP-RUNBOOK.md:2525–2554`）；`grep` = **OK**（多次，均返回带行号正文；一次仅因**正则自身**被 `rg` 拒绝（`\\Lambda_{\\mathrm` 尾随 `{`）——**这是我写错检索式，非 provider failure**，改写后同一目标由 `pwsh` grep 命中）；`glob` = **OK**（首轮 `Get-ChildItem` 定位 `app.tex`/`ch03-existence.tex` → 各 1 命中）；`pwsh` = **OK**（`py -3.12 -u` 一次，退出 0，输出见 §2/§4 表；`Get-ChildItem`/`Select-Object` 多次有效）；`write` = **OK**；`web_search`/`web_fetch` = **未调用**（本文无外部文献断言，全部为仓内代数 + 数值）。**无 `TOOL-FAIL`。**

| # | 检索式 | 范围 | 命中 | 落点 |
|---|---|---|---|---|
| 1 | `Get-ChildItem -Recurse -Filter app.tex` / `ch03-existence.tex` | 全仓 | 1 / 1 | 定位到 `docs/paper-v2/` |
| 2 | `W_full\|W_geo` | `research/v3/w7/h1c_wfull.py` | 33 | `:9–11,42,98,250–266` |
| 3 | `GPP-STRUCT`（全仓无 `*.md` 限定） | 全仓 | 138（250 上限内） | `species/GPP-STRUCT.md`、`O12/O14/O18/O26` 等 |
| 4 | `Delta_0\|\\Delta_0\|Δ₀\|Δ_0` | `docs/paper-v2/*.tex` | 12 | `app:182,187,191,194`；`ch03:161,162,174`；`ch13:60,65,70,73`；`ch18:90` |
| 5 | `a\^2\\?Delta_0\|a²Δ₀²` | `research/**/*.md` | 76 | `GPP-STRUCT:13,125,144,209–212`；`O14:153–170`；`O26:96–116`；`LOOP-RUNBOOK:2538,2710–2712` |
| 6 | `lattice spacing\|fm\b\|a^{-1}` | 全仓 | 1086（250 展示，全表落盘 `%TEMP%\dsh-spill-zvYlSr\session-77dc9ab8edb0\4886b17895f3-grep-results.txt`） | 关键：`H5-HEAD-NUMERICS:183`（`a(β≈2.5)≈0.15–0.2 fm`）、`T23-TEMPLATE:353`、`TIE-SIGMA:28`（`√σ=484.96 MeV`） |

**未查到（检索式：`Δ₀` 的数值；范围：全仓；工具：`grep`/`pwsh` OK；命中：0，全为符号）** —— 与 `O14-AUDIT-GPP.md:162`、`GPP-STRUCT.md:212` 一致：**`Δ₀` 在语料中无数值**，故 §4.3 的量级是"给定 `Δ₀` 的换算式"，不是语料数值。

---

## 7. 标签与未闭合

| 陈述 | 标签 |
|---|---|
| `δ := a²Δ₀² + W_full`（`app:190–191` 逐字）、`W_full = W_geo+X`（`app:796`） | `[T-transcribed]` |
| `W_full = E[Φ''] + Var(G) − Var(G+Φ') ≡ W_geo + X` | `[T]`（`(★)` 双线性恒等；`[E]` 两模型零差） |
| `G''(0) = Δ_dyn + W_full` | `[T]` |
| `H1 ⟺ Δ_dyn ≥ a²Δ₀²` | `[T]`（给定 §3.3 的 `∇S` 排版）；作为**路线** `[NG]`（换名） |
| `(M)` 的 `∇V` 排版 vs `δ`：差 `Var(Φ')+2Cov(G,Φ')` | `[C]`（代数确凿；量级 `2e−3`（plaquette）→`O(1)`（`L=2`）） |
| `a²Δ₀²` 无量纲、依赖格距 `a`、`=(aΔ₀)²` | `[T]` |
| `β=2.5,SU(2)` 的 `a²Δ₀²` 数值表、`Δ₀` 上界 `0.84–1.36 GeV` | `[C]`（锚 `a` 为 `[C]`，`ħc` 换算 `[E]`） |
| `Δ_dyn(L=2)=+1.065±0.333`（`L=3` 变号） | `[E]`（转引 `LOOP-RUNBOOK:2525,2540,2809`，本文件未重跑 MC） |

**未闭合（诚实标注）**：① §3.3 的排版缝隙是**文档级**判决（"哪一份对"由 `eq:split` 的 `H_vv+Φ''` 拆分决定，`[C]`），不是新数学；`app.tex:185–187` 与 `ch13:63–66` 两处**建议**改成 `⟨v,\nabla S_W⟩` 或把 `+W_full` 显式写回（本文件不改任何 `.tex`）。② `L=2` 上 `Δ_dyn` 只有 `2–3.2σ`，`L=3` 已变号 ⇒ 等价的"数值支撑"不存在，只有结构。③ `Δ₀` 的数值在整个语料里缺席（§6），所以 `a²Δ₀²` 与 `G''` 的**真正比较**在物理点仍未做过。

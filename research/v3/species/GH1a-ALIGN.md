# GH1a-ALIGN — `H1` 的四个版本对齐：它们**不是**逐字同一个命题

> 2026-11-25 · 子代理 `gh1a-align` · **只写本文件** · **未 git commit** · **无 MC** · No QED · **Clay 0%**
> 任务：把散在 `sec04-survey.tex` / `Q6-H1-STATUS.md` / `H1-FINAL.md` / `N11-FALLBACK.md` / `N04-2D.md` / `GW5-ACID.md` / `INEQ-R5` 里的四个"`H1` 版本"**逐字抄出**，写清量词次序，逐对判等价，并检查**载体错配**。
> 标签：`[T]` 真定理/恒等 · `[C]` 条件 · `[NG]` 具名类 no-go 或错判 · `[M]` 口径 · `[E]` 数值/精确算术 · `[O]` 精确缺失步 · `[J]` 判断。
> 前置声明：本轮**不产生新数学**。全部判决由**仓内恒等式**与**仓内反例**完成；凡未给出证明的方向，一律不标 `[T]`。

---

## §0 判词（一页）

**一句话**：四版里**没有任何两版是逐字同一个命题**。`(1)⟺(2)` 有一支是**恒等换名**（`[T]`，但**零信息**），其余六个方向全是 `[NG]` / `[O]`。语料反复指向的**不是** `H1`，而是 `O=¼trF²` 的 Källén–Lehmann 测度上的**原子 + 隔离**这一对 `(f₁², Λ★−m₁²)`；`H1` 只是它上方的**天花板**。

| # | 问 | 裁定 | 级 | 依据（file:line） |
|---|---|---|---|---|
| 1 | 四版是同一个命题吗 | **否。0/4 对等价**：`(1)⟺(2)` 只有"换名"这一支为 `[T]`（不携带信息）；`(2)⟺(3)`、`(1)⟺(3)`、`(3)⟺(4)` 均在**至少一个方向** `[NG]` | `[NG]` | §2 判决表 |
| 2 | `(1)⟺(2)` 的真相 | **(1b)（`(M)` 形）与 `(2)` 的 `Δ_dyn` 半边是同一句的两种写法**（移项恒等，零新估计）；但 `(2)` 还把 `G″` 写成**希尔伯特空间算子** `−D²+2adF_B`，那一半与 `Δ_dyn` **不等价** | `[T]`+`[NG]` | `app.tex:190–192`；`ch13-machine.tex:69–73`；`Q02-H1I2.md:44`；对照 `N04-2D.md:14,16` |
| 3 | `(2)⟺(3)` 成立吗 | **否，且语料自己已有 `[T]` 判词**：`H1 ⇎ I2` 两方向皆假（`Q02`），五名链已撤为三名（`RUNBOOK:9823–9842`） | `[NG]` | `Q02-H1I2.md:17–19`；`LOOP-RUNBOOK.md:9823,9841`；`:10430` |
| 4 | `(1)⟺(3)` 成立吗 | **否**：方向 A 有**带隙、无原子**的精确见证（`χ≤1/λ` 且 `inf supp ν=m₁>0` 且 `f₁²=0`）；方向 B 有 `I2` 真而 `χ=∞` 的见证 | `[NG]` | `R6-ATOM-ADVERSARIAL.md:254–260`；`Q02-H1I2.md:17` |
| 5 | `(3)⟺(4)` 成立吗 | **否，而且不是"未证"**：`(4)` 已被 `[T]` 判死（`lsc` 类的正 barrier 一个都不存在），`(3)` 仍是 `[O]`；把 `(3)` 读成 `(4)` 是把一个**开**命题换成一个**已假**命题 | `[NG]` | `GW5-ACID.md:125,127,130`；`INEQ-R5-donsker-varadhan.md:9,15` |
| 6 | 真身是哪一个 | **都不是 `H1`**。真身 = `(O=¼trF², f₁²≥r₀ ∧ I₂)` 这一对（原子地板 + 隔离）；`H1`（无论哪一版）是它的**必要天花板**，Q6 已判 `H1` **严格弱于**间隙 | `[T]` | `Q6-H1-STATUS.md:67–68,96,237–239` |
| 7 | "收敛到 `H1`"这句话 | **必须改**：改成"收敛到 `H*`，其中 `H1` 只占天花板那一半"（给出可直接粘贴的替换文本，§4） | `[J]` | §4 |
| 8 | 对齐可做吗 / 先验 | **可做，本轮已做完**。先验（读原文前，按语料自述"same object in several vocabularies"，`ch13-machine.tex:34`）：四版应逐字同一；后验：**0/4 对等价**，唯一 `[T]` 等价是恒等换名 | `[J]` | §1.6、§5 |
| 9 | 这是台账第几条同类错 | **1 条主账 + 2 条关联**：主账 = `Symbol polysemy (7/7)`（`H1` 一名多物，`sec05-ledger.tex:130–148` 该类现为 `(1/6)…(6/6)`）；关联① = `quantifier scope` 第 2 例（第 1 例 = 台账 #29，`sec05-ledger.tex:233–236`）；关联② = `Scope error` 第 6 条（该类现 `(1/5)…(5/5)`，`sec05-ledger.tex:176–190`，先例 #15）。**并触发台账 #5 的未执行补救**（"every `G''` statement must name its carrier"，`sec05-ledger.tex:140–141`） | `[J]` | §3、§4.3 |
| 10 | Clay | **0%** | — | — |

**必须先钉死的一条时间线（否则本节会被读成"新发现"）**：语料**已经**撤过这条等价链 ——
`LOOP-RUNBOOK.md:9841`（逐字）："| `H1 ≡ I2 ≡ (E-I2) ≡ H* ≡ Δ₀>0` | **`Λ★下界 ≡ I2 ≡ (E-I2)`**（三项） |"；
`:10430` 逐字："*The correct equivalence class is three, not five: `{Λ_★ lower bound ≡ I2 ≡ (E-I2)}`; `H*` and `Δ₀` do not belong to it.*"；
`H01.md:158` 把它记为 `[NG]`："**已撤**：`:9823`"；
`H01.md:170` 逐字："论文**两章各自写了一个 `H1`**：`ch02:89` 的 `χ≤1/λ` 与 `ch18:43` 的 `G″(0)≥δ`。二者**不是同一个命题**"。
⇒ 本轮**不是**发现新错，而是：**(a)** 首次把**四版原文 + 量词表**并排对齐（`H01` 只对齐了 `A/B/C` 三版，未含 `(3)`/`(4)`）；**(b)** 首次指出 `(3)` 的 `⟺` 在 `⟹` 方向**为假**（§2.2，新）；**(c)** 首次指出 `(2)` 内部"`Δ_dyn` 地板 ⟺ `G″` 谱隙"**为假**（§2.1，新）。

---

## §1 四版原文与量词表

### 1.1 `(H1-1)` 方差版（`sec04` 表第 17 行）

**表内逐字**（`docs/paper-final/sec04-survey.tex:109`）：

```
17 & $H1$ & a one-sided variance bound at typical field and $k\to0$ & \LN & $H1\not\Rightarrow m>0$: the thin atom witnesses a gap-free state satisfying $H1$\\
```

**该行的源**（`sec04-survey.tex:137`，逐字）：

```
17: \code{Q6-H1-STATUS.md}:208 (thin atom); \code{R6-ATOM-ADVERSARIAL.md}:254 (gapped, atomless); \code{LOOP-RUNBOOK.md}:19048--19053.
```

**操作形一**（`research/v3/Q6-H1-STATUS.md:76`，逐字）：

```
H1:    χ := Var_μ(ξ_{k→0}) = ∫_0^∞ s^{-1} dν(s) ≤ 1/λ ,   uniformly in (L,a),  β=2.5.     (1.1)
```

**操作形二**（`docs/paper-v2/ch02-squeeze.tex:89–90`，逐字）：

```
89: \mathrm{H1}:\qquad \operatorname{Var}_\mu(\xi_{k\to0})=\chi\ \le\ 1/\lambda,
90: \qquad\text{uniform in }(L,a).
```

**"at typical field" 的落点**（`docs/paper-v2/ch13-machine.tex:51–53`，逐字）：

```
51: the ``typical field $\times$ high momentum'' cell of the $(\text{field}\times\text{momentum})$
52: plane, requiring the ceiling $\operatorname{Var}(\xi_{k\to0})\le1/\lambda$ (curvature form
53: $\operatorname{Var}_\mu(\langle v,\nabla S\rangle|0)\le E_\mu[H_{vv}|0]+W_{\mathrm{full}}-\delta$;
```

> **注意（本轮第一处逐字冲突）**：`ch13:52–53` 把 `Var(ξ_{k→0})` 与 `Var_μ(⟨v,∇S⟩|0)` 写在**同一句的括号里**当作同一个 ceiling；而 `Q6-H1-STATUS.md:186–188` 逐字说 "force ≠ Hessian ≠ `χ`"。这是 §3 的 (a) 条。
> **另注**：`H1-FINAL.md:14`（row 1）给出的是 `H1(δ) ⟺ G″(0) ≥ δ`，即**曲率**形，不是方差形 —— 所以 `H1-FINAL` 的 row 1 不属于 `(1)`，而属于 `(2)` 的载体。任务书把 `H1-FINAL.md` 列为 `(1)` 的源，此处**只对"`H1 ⇏ m>0`"这一句成立**（`H1-FINAL.md:18`，thin atom）。

### 1.2 `(H1-2)` 谱隙版（`sec04` 表第 11 行）

**表内逐字**（`sec04-survey.tex:103`）：

```
11 & minimality & $\lambda_1(\nu)\ge\operatorname{ess\,inf}G''$ & \LN & reduces to $\operatorname{ess\,inf}G''=H1$, and the Savvidy tail \cite{savvidy,no,leutwyler} makes $\operatorname{ess\,inf}G''<0$\\
```

**"缺的那一步"的逐字**（`research/v3/species/H1-FINAL.md:20`）：

```
| — | 缺的那一步（唯一） | **`Δ_dyn` 的 `(L,a)`-一致正下界**（等价：粗模 `G″` 的 `(L,a)`-一致谱隙），即 `ch13-machine.tex:49–50` 那一格 | `[O]` | `N11-FALLBACK.md:§6`、`M05-RATIO.md:128–129`（"新靶的已证部分 = ∅"） |
```

**"换名"的逐字**（`research/v3/species/N11-FALLBACK.md:162`）：

```
- 唯一**降级**的建议：`Δ₀>0` 从"唯一缺口"降为"H1 的换名"（§0 ①）；`app.tex:182–194` 的等价链是 `[T]`，故"缺口"与 `H1` 是**同一件账**，不应在台账里占两格。
```

**载体恒等式的逐字**（`research/v3/species/N04-2D.md:14` 与 `:49`）：

```
14: | 1 | 2D 上 `G″` 是多少 | **不是一个数**：`G″ = W_geo + Δ_dyn`，`W_geo = −∂_φ²log J = 2/sin²φ ≥ 2`（class angle，**与面积无关**），`Δ_dyn = −∂_φ²log K_t = O(e^{−2t})`，`t=g²A/2`。**T03 的 `G″=0` 与 H1A:132 的 `+2.0101` 是同一模型、两个粗变量坐标、`t→∞` 的两个值，差一整块 `W_full`** | `[T]`+`[E]` |
49: ### [A] 逐点 `[E]`：恒等 `G″ = Δ_dyn + W_full` 在 2D 上精确成立（残差 `3.7e−40`）
```

**算子形**（`docs/paper-final/sec01-reduction.tex:310`，逐字片段）：`the Hessian of the action at a background, $G''=-D^{2}+2adF_B$, and its indefinite direc…`；其不定性见 `docs/paper-v2/ch18-scalar-h1.tex:184–186`："the Hessian $-D^2+2adF_B$ of the Savvidy cell is indefinite (a negative transverse magnetic mode surviving the gauge quotient; the measured negative fraction is $7.1\%$)"。

**双轴一致的逐字出处**（`docs/paper-v2/app.tex:182–188`）：

```
182: \LM{} The single open inequality, stated with no slack. There is $\Delta_0>0$ such that,
183: uniformly in $(L,a)$ and for every coarse transverse $v$ at $k\to0$,
184: \[
185:   \operatorname{Var}_\mu\!\big(\langle v,\nabla V\rangle\,\big|\,\xi=0\big)
186:   \;\le\;
187:   E_\mu\!\big[H_{vv}\,\big|\,\xi=0\big]\;-\;a^{2}\Delta_0^{2},
188:   \qquad V=S_W+\Phi .
189: \]
190: This is the (M) form of Chapter~\ref{ch:immune} with $W_{\mathrm{full}}=W_{\mathrm{geo}}+X$ and
191: $\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality
192: (the $W_{\mathrm{full}}$ cancels on both sides, since $G''_a=\Delta_{\mathrm{dyn}}+W_{\mathrm{full}}$).
```

### 1.3 `(H1-3)` 隔离版（`sec04` §9.4 末段）

`main.toc:39` 逐字 `\contentsline {subsection}{\numberline {9.4}What survives}{31}{subsection.9.4}%` ⇒ §9.4 = `sec04-survey.tex:221` 的 `\subsection{What survives}`。

**逐字**（`sec04-survey.tex:246–251`）：

```
246: The threshold form of that inequality is
247: \[
248: w_1\ge w_0 \iff \frac{\delta_1}{\varepsilon}\le 1-w_0,
249: \]
250: which is precisely the isolation of the second support point, i.e.\ $H1$. So the FKG route
251: does not deliver a floor either: it delivers $H1$, and $H1$ does not imply $m>0$ (row 17).
```

**源文件的逐字**（`research/v3/species/U09.md:160`）：

```
2. **唯一还开着的形状 = 把单调性问题换成谱矩坐标**：`w₁ ≥ ϑ = 1 − δ₁/ε` 已是 `[T]` 恒等式的推论；要 `w₁ ≥ w₀` ⟺ `δ₁/ε ≤ 1−w₀` ⟺ **第一矩不能贴近第二支撑点** ⟺ 第二能级隔离 `s₂>s₁`（`I2`）。**这与 `H1` 是同一个数**（`QI-BRACKET` 全篇；`GKS-GAUGE:24` 同判"amplitude→edge is `Λ`=isolation=H1"）。
```

**该等价所依赖的恒等式（逐字）**（`research/v3/species/QI-BRACKET.md:18–24`）：

```
18: 1. **桥的精确内容被一条恒等式钉死 `[T]`。** 令 `ν_w = Σ w_nδ_{s_n}`（概率测度）、`s₁ = inf supp ν_w`、
19:    `s₂ = inf(supp ν_w ∩ (s₁,∞))`（**第二支撑点**，T3 约定）、`Ē = ∫s dν_w`、`ε := s₂−s₁`、`δ₁ := Ē−s₁ ≥ 0`。则
20:    ```
21:       ϑ := (s₂−Ē)/(s₂−s₁) = 1 − δ₁/ε = w₁ − Δ/ε ,      Δ := ∫_{[s₂,∞)}(s−s₂)dν_w ≥ 0 .      (QI-1)
22:    ```
23:    **桥 = 二能级近似；它的 slack 恰是 `[s₂,∞)` 上的"超额一阶矩"`Δ` 除以能级间距。**
24:    数值（500 随机测度）：`max|ϑ − (w₁−Δ/ε)| = 3.6×10⁻¹⁵`，`max(ϑ−w₁) = 0`（不越界）。
```

### 1.4 `(H1-4)` 熵率前件（`GW5-ACID.md` §2）

**逐字**（`research/v3/species/GW5-ACID.md:111–130`，节选不影响逐字性）：

```
111: 设定：`Ψ` 是任意 **lsc** 泛函，路线希望 `Ψ` 给残数一个正的下界。`f₁²(ν) := ν({m₁²})`。
123: 要"用熵率给残数正下界"，逐字等于断言
124: ```
125:    ∃ Ψ ∈ lsc 类：  (i) Ψ(ν) ≤ f₁²(ν)  ∀ν      且   (ii) Ψ(ν_root) > 0 .                    (P)
126: ```
127: `(P)` 的否定（`¬(P)`）逐字是 `r₀ ≤ 0`：`INEQ-R5:15` "Positing a positive sup-barrier is verbatim `H*(r₀)`"。**注意：`(i)` 是单侧不等式，`(ii)` 是正性；两者缺一不可。**
130: 取任意 `ν*`，把它在 `m₁²` 处的边缘原子替换为宽 `1/n` 的 bump：`ν_n ⇐ ν*`（弱收敛），`ν_n({m₁²}) = 0`。则 `Ψ(ν_n) ≤ f₁²(ν_n) = 0`（用 `(i)`）；lsc 给 `Ψ(ν*) ≤ liminf_n Ψ(ν_n) ≤ 0`。`ν*` 任意 ⇒ **`Ψ ≤ 0`**（`INEQ-R5:9` 逐字）。⇒ `(ii)` 不可满足 ⇒ `(P)` 假 ⇒ **`r₀ ≤ 0`**。
```

**源**（`research/unconditional/INEQ-R5-donsker-varadhan.md:7,9,15`，逐字）：

```
7: **2. The target is an atom — wrong class.** `f₁²=ν({m₁²})` is the mass of a point. … hence `f₁²` is not a sup of continuous functionals — *a fortiori* not a DV rate. DV and the atom are dual cones.
9: **3. Forcing lemma (uniform over the whole sup class).** … **`Ψ≤0`** — no sup-representable barrier is positive, i.e. `r₀≤0`.
15: … Positing a positive sup-barrier is verbatim `H*(r₀)`.
```

### 1.5 量词表（本任务的核心）

记号：`x := (L,a)`；`v` = 粗横向模；`U′` = 规范组态；`Typ` = "at typical field"（`ξ=0` 条件化）。

| 版 | 载体（**到底是什么对象**） | 界的方向 | 一致的轴 | 量词次序（**符号**） | `(L,a)` 双轴？ | 逐点还是全类 |
|---|---|---|---|---|---|---|
| `(1a)` `H1_χ` | `χ = Var_μ(ξ_{k→0}) = ∫s⁻¹dν`（测度的一阶矩，`O=¼trF²`） | **上界（天花板）** | `(L,a)`；`β=2.5` 钉死 | `∃λ>0 ∀(L,a)` | **是** | 全类（`Var_μ` 已把场积掉） |
| `(1b)` `(M)` 形 | `Var_μ(⟨v,∇V⟩|ξ=0)` 对 `E_μ[H_vv|ξ=0]`（**力方差 vs Hessian，条件在典型场**） | **下界**（对 `Δ_dyn` 而言） | `(L,a)` **+ `v`** | `∃Δ₀>0 ∀(L,a) ∀v` | **是** | 典型场（`ξ=0` 条件化） |
| `(2)` `Δ_dyn` 地板 | `Δ_dyn := E[H_vv|0] − Var(G_v|0)`（曲率，标量逐模） | **下界** | `(L,a)` **+ `v`** | `∃δ>0 ∀(L,a) ∀v` | **是** | 逐模（"mode by mode"） |
| `(2′)` `G″` 谱隙 | **希尔伯特空间算子** `−D²+2adF_B`（`(μ,a)` 混合指标） | 谱底 | `(L,a)` **+ 背景场** | `∃δ ∀(L,a) ∀(背景)`：`inf spec ≥ δ` | 名义上是 | 算子谱（全谱底） |
| `(3)` 隔离判据 | 测度 `ν` 的**三点谱数据** `(s₁,s₂,Ē)` | **等价式**（不是不等式） | **无轴**（单条测度） | `∀w₀∈[0,1]`：`(w₁≥w₀ ⟺ δ₁/ε ≤ 1−w₀)` | **否**（轴必须外加） | 逐测度 |
| `(4)` `H*(r₀)` | 测度空间上的**泛函类** `lsc` 中的 `Ψ` | **下界**（`Ψ ≤ f₁²`）+ 单点正性 | **无轴**（测度空间自身） | `∃Ψ∈lsc ∀ν`，且 `Ψ(ν_root)>0` | **否** | 全测度类 |

**三条由该表直接读出的事实**：

1. **只有 `(1b)`、`(2)` 把量词写成 `∃δ ∀L∀a∀v` 的"好次序"**，而且这是语料的**明文**（`app.tex:182–183`："There is `Δ₀>0` such that, **uniformly in `(L,a)` and for every coarse transverse `v`**"，同一句在 `ch13-machine.tex:60–61` 重复）。这是仓内**唯一**把 `∃δ∀(L,a)` 写全的地方。
2. **`(3)`、`(4)` 根本没有 `(L,a)` 量词**。`sec04:248`、`U09:160` 印出来的是**单条测度上的等价式**；"一致版"必须外加 `∃w₀>0 ∀(L,a)`（或 `∀(L,a)` 的 `Ψ`），而**这一步在原文里没有**。`sec04:250` 的 `i.e. $H1$` 就是在这个缺口上搭的桥：它把**无轴命题**叫成 `H1`，而 `H1` 的两版都是**有轴**的。
3. **裸 `H1` 与均匀 `H1` 是两个命题**，`ch18` 自己把它们分开写：`ch18-scalar-h1.tex:43–44` 逐字 `H1(\delta)\iff G''(0)\ge\delta \qquad\text{(pointwise, mode by mode, for every }(\beta,U')\text{)}`（逐点 `∀x∃δ(x)`）；而 `:78–79` 逐字 "the open content is now (i) the sign and size of `G''(0)` at the physical point and (ii) **its `(L,a)`-uniformity**, and (ii) is the part that is not closed"（`∃δ∀x`）。⇒ **语料明文承认：`∀∃` 升到 `∃∀` 就是那个未闭合的缺口**。这正是本任务要的"`∀L∃δ` 与 `∃δ∀L` 是两个命题"。

### 1.6 `(2)` 是**两条**命题被 `⟺` 粘起来

`H1-FINAL.md:20` 的原句把两件事写成等价：`Δ_dyn` 的一致正下界 **⟺** 粗模 `G″` 的一致谱隙。拆开：

- `(2Δ)`：`∃δ>0 ∀(L,a)∀v: Δ_dyn ≥ δ`；
- `(2G)`：`∃δ>0 ∀(L,a): inf spec(−D²+2adF_B) ≥ δ`（或 `ess inf G″ ≥ δ`，两种读法在各文件里混用）。

由 `G″ = Δ_dyn + W_full`（`[T]` 恒等，`N04-2D.md:14,26,49`；机检残差 `6.7e−16`，`H1-FINAL.md:48`）可知：`(2Δ) ⇒ (2G)` 需 `W_full ≥ −c`（`c<δ`）；`(2G) ⇒ (2Δ)` 需 `W_full ≤ C`（上界）。**两者仓内都没有 `[T]`**，而且 `W_full` 本身是**两个对象同名**（台账 #10，`sec05-ledger.tex:159–161`："`W_full` and `W_geo` are two objects under one name, used as one object. → Prékopa's `W_geo ≥ 0` does not transfer to `W_full`"）。⇒ `(2)` 的 `⟺` 是**载荷所在**，见 §2.1 的 `[NG]`。

---

## §2 逐对等价判决

**判据口径（先写死，避免"未证标 `[T]`"）**：`[T]` 只给**代数恒等**或**仓内已证定理**；`[NG]` 要给**见证**（模型内或测度级）；`[O]` 要写清**缺哪一步**。

### 2.1 `(1)⟺(2)`

| 方向 | 判 | 证明 / 反例 / 缺失步 |
|---|---|---|
| `(1b) ⇒ (2Δ)` | `[T]` | 纯移项，**零估计**：`G″ = E[H|0] − Var(G|0) + W_full`（`ch18-scalar-h1.tex:30–36`；`N04-2D.md:26` `eq:split`）。`(M)`：`Var ≤ E[H] − a²Δ₀²`。两边同加 `W_full` 再移项 ⇒ `E[H] − Var + W_full ≥ a²Δ₀² + W_full` ⇒ `G″ ≥ δ`，`δ:=a²Δ₀²+W_full` ⇒ `Δ_dyn ≥ a²Δ₀²`。语料同判：`app.tex:190–192`、`ch13-machine.tex:69–73`、`Q02-H1I2.md:44`（"零新估计，只是换名"）、`P03-EQUIV.md:253–258` |
| `(2Δ) ⇒ (1b)` | `[T]` | 同上（同一句，反方向移项） |
| **净判** | **`[T]` 但是"换名"，不是等价** | 两条式子**共用同一组未知量**（`E[H]`、`Var`、`W_full` 一个都没被单独界住）；`ch18-scalar-h1.tex:59–63` 逐字："`Var`, `E[H]` and `W_full` are **not separately constrained** by `H1` beyond their appearance in (eq:scalar-identity); … A proposal that bounds `Var` without controlling the sign of `G''(0)` is a claim about `G''(0)` in disguise" ⇒ **信息量 = 0** |
| `(2Δ) ⇒ (1a)` | `[O]` | 需两步仓内都没有的件：①`W_full ≥ −c`（`c<δ`）的一致下界 —— 台账 #10 已判 Prékopa 的 `W_geo≥0` **不传给** `W_full`（`sec05-ledger.tex:159–161`）；且 `W_geo ≥ (N/12)‖v‖²` 在归一化 `‖v‖²=1/L` 下地板 `(N/12)/L→0`（`H1-FINAL.md:61`；`ch05-theorems.tex:224`）；`W_full ≥ N/12` 是 `[V]`/`[M]` 且**机检 144/144 点 `W_full < N/12`，min `−0.157`**（`H1-FINAL.md:61`）。②需 `ess inf` 读法（`Q6-H1-STATUS.md:145–146` `[T]`："`ess inf G″≥λ ⇒ H1`"）；**典型场/点值读法已被判 `[NG]`**（`Q6-H1-STATUS.md:147–150`："A **point** value `G″(a₀*)≥λ` is not an `ess inf`"） |
| `(1a) ⇒ (2Δ)` | `[NG]` | 见证 W1（薄原子）：`ν = wδ_ε+(1−w)δ_{s₂}`，`w=ε/2λ` ⇒ **`χ=1/λ` 精确**（`Q6-H1-STATUS.md:114`；`H1-FINAL.md:121` 逐字"解出 `s₂=2λ(1−w)`、`χ=1/λ`"），而 `(1a)` 对每个 `ε>0` 都成立、**同一个 `λ`**；同时取平坦坐标 `W_full ≡ 0`（`N04-2D.md:15` 逐字"平坦坐标里恒为 **0**"）⇒ `Δ_dyn = G″ =` 该点曲率 `= ε`（曲率=该点稳定性，`ch18-scalar-h1.tex:182–183` 逐字 "``the curvature at the conditioned point'' *is* ``the stability of that point''"）⇒ 沿 `ε↓0` **无 `δ>0`** |
| `(2G)`（算子/`ess inf` 读法） | `[NG]` | ①**坐标依赖**：同一 2D 模型给出 `G″=0` 或 `+2.0101`，只取决于粗变量坐标（`N04-2D.md:14` 逐字"同一模型、两个粗变量坐标"）；②**`ess inf G″<0`**：`sec04-survey.tex:103` 逐字 "the Savvidy tail … makes `ess inf G″<0`"，数值 `K_S(2.5)=−119`（`GAP2.md:72`）；`ch18-scalar-h1.tex:184–186` 逐字 7.1% 负模；③**写法的良定性**：`−D²+2adF_B` 作色矩阵**不自伴**（`F03-2DTEST.md:26–29` 逐字"`ad F` 反对称，本征值 `±iθ_α`，`−D²+2adF` 不自伴（"谱"是复数）。自伴对象是**带指标混合**的 `H_μν^{ab}`"） |

**结论**：`(1)⟺(2)` **作为"两个不同命题"为 `[NG]`**；**作为"同一个式子的两种写法"为 `[T]` 但空转**。两种读法都**不支持**语料的"等价"用法（后者不携带信息，前者假的）。

### 2.2 `(2)⟺(3)`

| 方向 | 判 | 证明 / 反例 / 缺失步 |
|---|---|---|
| `(2Δ) ⇒ (3)` | `[NG]` | 见证 W2（`:φ²:`，自由大质量复合场）：`ν` 是 `[4m²,∞)` 上的**割测度**（`f₁²=0`，`Q6-H1-STATUS.md:117–119` 逐字"gapped (`m₁²=4m²>0`), but `ν` is a **cut** on `[4m²,∞)`, `f₁²=0`"），高斯边缘 ⇒ `G″=1/χ>0`（`Q6-H1-STATUS.md:154` 逐字"Under `η=0`, `G″=1/χ` exactly"），故 `(2Δ)` 真；而 `(3)` 在此处**无定义**：`ε = s₂−s₁ = 0`，`δ₁/ε` 是 `0/0` 不定（`QI-BRACKET.md:26–31` 逐字"**`ε=0` 处三元组 `(s₁,s₂,Ē)` 对 `w₁` 零信息**…**0/0 是真实的不定性，不是可去奇点**"）⇒ `(3)` 恰在"无原子"处失效，**不能**是 `(2)` 的等价物 |
| `(3) ⇒ (2Δ)` | `[NG]` | 见证 W3（测度级，两原子）：`ν = ½δ₀ + ½δ₂` ⇒ `s₁=0`、`ε=2`、`Ē=1`、`δ₁=1` ⇒ `δ₁/ε = ½ = 1−w₀`，**取 `w₀=½` 时 `(3)` 成立**；而 `s₁=0` ⇒ 平坦坐标下曲率 `=0` ⇒ `Δ_dyn=0`，**无正地板**。同一见证也否定 `(3)⇒(1a)`（`χ=∫s⁻¹dν=∞`）。**模型内同型见证**：2D YM 的 `w₁=∫…`（`O₂=½trU` 给精确 `f₁²=1`，`ATOM2.md:80–81`）而 `Δ_dyn=O(e^{−2t})` 变号且 `sup|Δ_dyn|=24e^{−2t}→0`（`N04-2D.md:16`） |
| **`(3)` 自身的 `⟺`（`⟹` 方向）** | **`[NG]`（本轮新）** | `(3)` 的原文是"`w₁ ≥ w₀ ⟺ δ₁/ε ≤ 1−w₀`"。由 QI-1（`[T]` 恒等）：`w₁ = ϑ + Δ/ε`，`Δ ≥ 0` ⇒ **`⟸` 真**（`δ₁/ε ≤ 1−w₀ ⇒ w₁ ≥ ϑ ≥ w₀`）；**`⟹` 假**：只要 `Δ>0` 就能用尾矩补足，`ϑ` 可以远低于 `w₀`。**仓内数值见证**（`QI-BRACKET.md:36–38` 逐字）：`L=8: ϑ = −1.74  w₁ = 1.05e−1`；取 `w₀=0.10` ⇒ `w₁ ≥ w₀` **真**而 `δ₁/ε = 1−ϑ = 2.74 ≤ 0.90` **假**。⇒ 正确的写法是**单侧**：`δ₁/ε ≤ 1−w₀ ⟹ w₁ ≥ w₀`，`⟺` 只在 `Δ=0`（两点测度）时成立 |
| **语料自己的判词** | `[T]`（非等价） | `Q02-H1I2.md:17–19`：`I2 ⇒ H1` **假**（见证 A `ν=δ₀+1_{[δ,δ+1]}`；见证 B Savvidy 双井 `I2: 2.76→27.92` 而 `χ→∞`）；`H1 ⇒ I2` **假**（Lean 反模型 `f=0` 支）；"**两个方向都不成立 ⇒ 不可比**" `[T]` |

### 2.3 `(1)⟺(3)`

| 方向 | 判 | 证明 / 反例 / 缺失步 |
|---|---|---|
| `(1a) ⇒ (3)` | `[NG]` | 见证 W4 = R6-a（`[T]` 定理，`R6-ATOM-ADVERSARIAL.md:254–255` 逐字）："For every `λ>0`, `M₀>0`, `m₁>0` there is a positive measure `ν` with `M₀(ν)=M₀`, `inf supp ν = m₁` (**gapped, fixed**), `f₁²=ν({m₁²})=0` (**atomless**), and `χ=∫s^{−1}dν ≤ 1/λ`." ⇒ `(1a)` **真**、间隙**真**、原子**零** ⇒ 对任意 `w₀>0`，`w₁ ≥ w₀` **假** |
| `(3) ⇒ (1a)` | `[NG]` | 见证 = `Q02-H1I2.md:17` 见证 A：`ν = δ₀ + 1_{[δ,δ+1]}`：`I2` 真（`Λ★−m₁² = δ`）而 `χ = ∫s⁻¹dν = ∞`（`δ₀` 处）⇒ `(1a)` 假。简化版即 §2.2 的 `ν=½δ₀+½δ₂` |
| `(1b) ⇒ (3)` | `[NG]` | 同 §2.1 的 W1：`ε↓0` 时 `w₁ = ε/2λ → 0`，而 `(M)` 的松弛 `a²Δ₀²` 可在 `ε` 的同一极限下保持不动（`Δ₀` 是**存在量词常数**，`H1-FINAL.md:87` 逐字"`Δ₀` 是一个**存在量词常数**"）⇒ 两个极限方向相反 |

### 2.4 `(3)⟺(4)`

| 方向 | 判 | 证明 / 反例 / 缺失步 |
|---|---|---|
| `(3) ⇒ (4)` | `[NG]` | `(4)` 已被 `[T]` 判死：任意 `lsc Ψ ≤ f₁²` 必 `Ψ≤0`（`GW5-ACID.md:130` L5；`INEQ-R5-donsker-varadhan.md:9` 强制引理）⇒ `(4)` 的**类为空**。故 `(3)` 即便处处为真，也交不出 `(4)`；"原子地板 + 隔离 ⇒ 存在 `lsc` barrier"这一步**不是缺估计，是类错配**（`usc` vs `lsc`，`INEQ-R5:7`） |
| `(4) ⇒ (3)` | 空真、**零内容** | `(4)` 假，故蕴含空真；且 `(4)` 的直接后果只有 `f₁²(ν_root) ≥ Ψ(ν_root) > 0`（正原子质量），**不含隔离 `ε>0`** ⇒ 即使 `(4)` 真也只够 `(3)` 的一半 |
| **净判** | `[NG]` | 把 `(3)` 标成 `(4)` 的战果，等价于把一个**开**命题记成一个**已证否**命题的推论 |

### 2.5 判决表（`§0` 与回答所贴版本）

| 对 | 方向 1 | 方向 2 | 净判 |
|---|---|---|---|
| `(1)⟺(2)` | `(1b)⇒(2Δ)` `[T]`（移项恒等，零估计） | `(2Δ)⇒(1b)` `[T]`（同句）；`(2Δ)⇒(1a)` `[O]`；`(1a)⇒(2Δ)` `[NG]`（W1 薄原子） | **恒等换名 `[T]`，但不是等价**；按"两条命题"读 ⇒ `[NG]` |
| `(2)⟺(3)` | `(2Δ)⇒(3)` `[NG]`（W2 `:φ²:`，`ε=0` 处无定义） | `(3)⇒(2Δ)` `[NG]`（W3 `½δ₀+½δ₂`；2D 同型） | `[NG]`；语料自身 `H1⇎I2` `[T]` |
| `(1)⟺(3)` | `(1a)⇒(3)` `[NG]`（W4 R6-a：带隙+无原子+`χ≤1/λ`） | `(3)⇒(1a)` `[NG]`（Q02 见证 A：`I2` 真、`χ=∞`） | `[NG]` |
| `(3)⟺(4)` | `(3)⇒(4)` `[NG]`（`(4)` 被 L5 判死，类为空） | `(4)⇒(3)` 空真、无隔离 | `[NG]` |
| **`(3)` 内部** | `RHS⇒LHS` `[T]`（QI-1 + `Δ≥0`） | `LHS⇒RHS` **`[NG]`（新）**（QI-BRACKET:36–38 的 `L=8` 数） | **原文的 `⟺` 不成立** |

### 2.6 为什么"语料把它们叫等价"是**同名不同物**，不是笔误

三条互相独立的证据：

1. **载体**：`(1a)` 是 `χ`（测度的一阶矩）；`(1b)` 是**力方差 vs Hessian**（`Var(⟨v,∇V⟩|0)` 与 `E[H_vv|0]`）；`(2)` 是**曲率** `Δ_dyn`；`(3)` 是**支撑点间距** `s₂−s₁`；`(4)` 是**泛函类**里的元。五个载体，五个量纲/类型。
2. **语料自己的载体警告**（`Q6-H1-STATUS.md:58–64` 逐字）："**Carrier hazard** … on the gauge-invariant composite `tr F²`, `Var(q̄)=χ/V→0`, so H1 is **vacuously true** and carries no gap content; on the gluon … Both branches are reported … Hence **"H1 ⟺ gap" presupposes a gauge-invariant, non-self-averaging carrier whose linear existence is exactly the structural theorem that does not hold.**" ⇒ 语料**已明说**该等价预设了一个它没有的载体。
3. **先例**：台账 #15（`sec05-ledger.tex:179–180` 逐字）"$\Phi$ is attached to the wrong carrier: $\Delta_0$ sits on $G''$, not on the measure $\nu$. & The name $\Delta_0>0$ is withdrawn"，以及 `sec05-ledger.tex:70–73` 的撤回原文（"the name sat on a carrier that is not the carrier of the statement"）。⇒ 本任务是**同一类**：名字（`H1`）坐在四个不同载体上。

---

## §3 载体错配检查（七处，逐处带 `file:line`）

| # | 错配 | 逐字证据 | 后果 / 级 |
|---|---|---|---|
| (a) | **三个"方差"**：`χ` / 力方差 / `E[H_vv]` | `ch13-machine.tex:52–53` 把 `Var(ξ_{k→0})` 与 `Var_μ(⟨v,∇S⟩|0)` 写在同一括号里；`Q6-H1-STATUS.md:186–188` 逐字 "**force ≠ Hessian ≠ `χ`**" | 三个量纲相同、内容不同的对象被当成一个 ceiling。`(1)` 的"variance"在三处指尖之间滑动 ⇒ `[NG]`（作为同一命题） |
| (b) | **`G″`：数 vs 算子** | `ch18:43`（数：`G″(0)≥δ`，逐模）vs `sec01-reduction.tex:310`（算子：`G''=-D^{2}+2adF_B`）vs `ch18:184–186`（该算子不定，7.1% 负） | `(2)` 的两半永远不可能同真同假 ⇒ `[NG]`（§2.1） |
| (c) | **`G″ = 1/χ` 与 `G″ = (1−η)/χ`** | `Q6-H1-STATUS.md:83` 逐字 "**The four faces (all one object).** `G″(k=0) = V_eff″(0) = Γ″(0) = 1/χ = 1/C̃(0)`" **对照** `:151–155` 逐字 "`G″=(1−η)/χ`；`χ≤1/λ` gives `G″≥(1−η)λ` **only** … Exact equality `G″·χ=1` needs `η=0` (Gaussian)" | **同一文件内**两处冲突；`η` 实测 `2–5%`（`Q6:152`）⇒ "all one object" 只在 `η=0` 成立。`[NG]`（作为恒等） |
| (d) | **`W_geo` vs `W_full`** | 台账 #10（`sec05-ledger.tex:159–161`）；`N04-2D.md:15` 逐字 "`W_full` … 是**运动学量**…平坦坐标里恒为 `0`" | Prékopa 的 `W_geo≥0` 不传给 `W_full` ⇒ `(2)` 的 `⟺` 直接断在这里。`[NG]` |
| (e) | **两个 `H*`** | `GW5-ACID.md:125` 的 `H*(r₀)`（`lsc` barrier，已判死）vs `Q6-H1-STATUS.md:67–68,96` 的 `H*`（充分对 `(f₁²≥r₀ ∧ I₂) ∧ (χ≤1/λ)`） | 同名两物；`(4)` 用的是前者，`(3)` 指向后者 ⇒ `(3)⟺(4)` 的"等价"从一开始就不是同一个符号 |
| (f) | **两个 `Δ_dyn`** | 台账 #7（`sec05-ledger.tex:146–148`）逐字 "$\Delta_{\mathrm{dyn}}$ carries two definitions, one of them $\Delta_{\mathrm{dyn}}=-W_{\mathrm{full}}$ on 2D Yang--Mills"；`H1-FINAL.md:96` 的 2D 行写 "`= 0`（取等）"，`N04-2D.md:14,16` 写 `O(e^{−2t})` 且**变号**（`inf_φ = −24e^{−2t}/a²`） | 2D 上 `Δ_dyn` 到底是 `0` 还是 `O(e^{−2t})` 变号 —— 两个定义都在用。`[M]`+`[E]` 口径差 |
| (g) | **`G″` 的坐标依赖** | `N04-2D.md:14` 逐字 "**T03 的 `G″=0` 与 H1A:132 的 `+2.0101` 是同一模型、两个粗变量坐标、`t→∞` 的两个值，差一整块 `W_full`**"；`:43` "`G″` 本身**依赖粗变量坐标**" | `(2G)` 的"谱隙"不是模型的函数，是**模型+坐标**的函数 ⇒ 不能作为 `(2)` 的等价半边 |

**载体错配的净结论**：`(3)` 说两个支撑点、`(2)` 说算子谱、`(1)` 说方差 —— 三者**没有共同的载体**，也没有共同的量纲。语料把它们称为"等价"，机制与台账 #15 完全相同：**名字代替了载体**（`RUNBOOK:9848` 自己的话："'形状相同'不等于'同一命题' —— 这是 AGENTS §2.3「算式对 ≠ 有内容」的第二次犯"）。

---

## §4 对"收敛到 `H1`"这句话的修正

### 4.1 该句在仓内的三处出现

| 位置 | 逐字 | 问题 |
|---|---|---|
| `sec04-survey.tex:250–251` | "which is precisely the isolation of the second support point, i.e.\ `$H1$`. So the FKG route does not deliver a floor either: it delivers `$H1$`" | 把**无轴**的隔离判据叫成 `H1`（`(L,a)` 量词缺失，§1.5） |
| `ch13-machine.tex:49–50` | "These are not six problems; they are one --- the volume-uniform lower bound on `m_{0^{++}}` along the path, which is the mass gap itself." | 六项里含 `(1)`、`(2)`、`(3)` 型对象，且 `:34` 说 "These five conditions define the same object in several vocabularies" |
| `H1-FINAL.md:23` | "`H1` 与 `Δ₀>0` **是同一句话**" | 与 `:157`（"五名链 ⇒ 改为三名 `{Λ★下界, I2, (E-I2)}` + 两名另列"）自相矛盾；该项**未落盘** |

### 4.2 替换文本（逐字可粘贴）

> **旧**："所有路线都收敛到 `H1`。"
>
> **新**："所有路线都收敛到同一个**对**：`H* = (f₁² ≥ r₀ ∧ I₂) ∧ (χ ≤ 1/λ)`（`Q6-H1-STATUS.md:67–68,96`）。其中**硬的一半**是 `O=¼trF²` 的 KL 测度在 `m₁²` 处的**原子 + 隔离** `(f₁², Λ★−m₁²)`；**软的一半**是天花板 `χ ≤ 1/λ`。`H1` 只是**软的那一半**，且 `H1 ⇏ m>0` 已证（`sec04:109`、`Q6:26–30`、Lean `H1_does_not_imply_atom` 0 sorry）。因此**不许**再写"`H1` = 间隙"、也不许把 `(3)`/`(4)` 记作 `H1` 的战果。"

### 4.3 三条必须跟着改的次级句子

1. `ch13-machine.tex:34` 的 "These five conditions **define the same object** in several vocabularies" ⇒ 加限定："same object **as a target**（同一靶），not the same proposition（非同句）"，或按 `H01.md:174–177` 的 A/B/C + 主版写法。
2. `sec04-survey.tex:248` 的 `\iff` ⇒ 改成 `\Longleftarrow`（`⟸`），并把 `Δ` 的 slack 写出来（`QI-BRACKET.md:21–23`）。
3. `H1-FINAL.md:20` 的 `⟺`（`Δ_dyn` 地板 ⟺ `G″` 谱隙）⇒ 改成"`Δ_dyn` 地板 ⇒ `G″` 地板（需 `W_full≥−c`）"，并注明 `W_full` 的一致下界是 `[O]`（台账 #10）。

### 4.4 台账落账（按 `sec05-ledger.tex:95–97` 的计数规则）

计数规则逐字："*One entry* is one minimal unit of debt such that copying it literally produces a wrong conclusion. A polysemous symbol counts once per symbol, with the number of senses in brackets"。

- **主账 1 条**：`Symbol polysemy (7/7)` —— `H1` 有 **5** 个义（`H1_χ` / `(M)` 力形 / `Δ_dyn` 地板 / `G″` 算子谱 / 隔离判据），照抄任一处都会产出错结论。该类现记 `(1/6)…(6/6)`（`sec05-ledger.tex:130–148`），加本条须改 `(7/7)`。
- **关联 1**：`quantifier scope` 第 **2** 例（第 1 例 = 台账 #29/item 38，`sec05-ledger.tex:233–236`）；本条记的是 `∀∃` vs `∃∀` 的**次序**，与 #29 的"`O` 固定与否"不同项。
- **关联 2**：`Scope error` 第 **6** 条（该类现 `(1/5)…(5/5)`，`sec05-ledger.tex:176–190`），先例 = #15（`:179–180`）。
- **不新增第三类**：`Symbol polysemy` 与 `Scope error` 的界线由 #5/#10 与 #15 划好（前者"一名多物"，后者"名字坐到错载体"），本条两者都沾；按"一条 debt"取**主账只记一次**，避免把一次错记成三次。
- **台账 #5 未闭环**：`sec05-ledger.tex:140–141` 的补救逐字是 "every `G''` statement must name its carrier" —— 本轮 `(2)`（`H1-FINAL.md:20`）与 `(1b)`（`app.tex:185–187`）都没有给 `G″` 命名载体（数？算子？哪个坐标？）⇒ #5 的补救**从未执行**。

---

## §5 工具台账（`AGENTS.md §1` 四项）

| 项 | 记录 |
|---|---|
| **工具状态** | `pwsh`（`Get-ChildItem` / `Get-Content` / `Select-String`）**OK**，全部返回有效输出 · `read` **OK** · `grep`（ripgrep）**部分 FAIL，见下** · `write` **OK**。**无 `ENOENT`、无空响应。** |
| **`TOOL-FAIL` #1** | `grep` 工具，检索式 `H1\|variance\|spectral gap\|Delta_{\\mathrm{dyn}}\|\\Delta_\{?dyn\|two support\|support point\|What survives`，`path=docs\paper-final\sec04-survey.tex` ⇒ **`rg: regex parse error: repetition quantifier expects a valid decimal`**（`\{?` 不被接受）。**记 `TOOL-FAIL`，不记"未查到"。** 改用无转义式重试成功（16 命中）。 |
| **`TOOL-FAIL` #2** | `pwsh`，检索式 `Select-String -Pattern 'H1 ≡'`，`scope='.'`（全仓递归含 `vendor/`）⇒ **超过 120000 ms 未完成，被移入后台 `pwsh-25baef8693aa4c79`，由本代理 `job_kill` 终止**。**记 `TOOL-FAIL`（超时），不据此下任何"0 命中"结论。** 改为限定目录重跑（见下 S9/S10）。 |
| **检索式 S1** | `read docs\paper-final\sec04-survey.tex` 全 258 行（逐行）⇒ 表 22 行、§9.4、`\src` 全读 |
| **检索式 S2** | `read research\v3\species\H1-FINAL.md` 全 198 行；`read research\v3\Q6-H1-STATUS.md:1–70,95–164,185–244`；`read GW5-ACID.md:100–149`；`read QI-BRACKET.md:11–40`；`read U09.md:118–162`；`read Q02-H1I2.md:1–45`；`read N04-2D.md:1–60`；`read H01.md:1–60,149–190`；`read INEQ-R5-donsker-varadhan.md` 全 17 行 |
| **检索式 S3** | `grep 'H1\|variance\|spectral\|dyn\|support\|survives' @ sec04-survey.tex` ⇒ **16 命中** |
| **检索式 S4** | `grep '^## \|^### \|G''\|G″\|Delta_\{?dyn\|Δ_dyn\|spectral\|谱隙' @ N11-FALLBACK.md` ⇒ **23 命中** |
| **检索式 S5** | `grep 'H\*\|H\*\(r\|Ψ\|Psi\|f₁²\|Donsker\|Varadhan\|H1' @ GW5-ACID.md` ⇒ **28 命中** |
| **检索式 S6** | `grep '1-w_0\|w₀\|δ₁\|delta_1\|support point\|支撑点\|隔离' @ U09.md` ⇒ **5 命中** |
| **检索式 S7** | `grep 'object mismatch\|same-name\|ledger\|\\Phi' @ sec05-ledger.tex` ⇒ **23 命中** |
| **检索式 S8** | `grep 'adF_\{?B\|2adF' @ paper-final+paper-v2+species` ⇒ **`docs\paper-final` 1 / `docs\paper-v2` 1 / `research\v3\species` 27** |
| **检索式 S9** | `Select-String '∃δ∀L'` / `'∀L∃δ'`，scope = `docs\paper-final`,`docs\paper-v2`,`research\v3\species`,`research\unconditional`（`*.md;*.tex`）⇒ **0 / 0 / 0 / 0** 与 **0 / 0 / 0 / 0**。**采信（工具 OK；命中 0）**：语料**从不**把量词次序写成符号，只用散文写（`app.tex:183` "uniformly in `(L,a)` and for every coarse transverse `v`"）。 |
| **检索式 S10** | `Select-String '∀(L,a)'` 同 scope ⇒ `species` **1 命中**（`H1-FINAL.md:86` 的 `∃Δ₀>0 ∀(L,a)∀v`），其余 **0** |
| **检索式 S11** | `Select-String 'w₁ ≥ w₀'` ⇒ `species` **39**；`'δ₁/ε'` ⇒ `species` **11**；`'Δ_dyn ≥ δ'` ⇒ `species` **3**（`species` 为 `research\v3` 子集，行数按目录分别计） |
| **检索式 S12** | `Select-String 'H1 ≡'` ⇒ `species` **6**、`unconditional` **1**；`'H1≡'` ⇒ `species` **0**、`unconditional` **2**；`'define the same object'` ⇒ `docs\paper-v2` **1**（`ch13-machine.tex:34`）、`species` **1**（`H1B-FIVE.md:90` 转引）；`'五名链'` ⇒ `species` **10**；`'两名另列'` ⇒ **0**（全 4 个 scope） |
| **未查到（采信）** | **"四个 `H1` 版本逐字同一"的证明**：S3–S7 的 **95** 命中逐条读过 + `H01.md` 全读，**0 处**（工具 OK）。**"`(3)` 的 `⟹` 方向的证明"**：S6 的 5 命中 + `QI-BRACKET.md` 全读，**0 处**（工具 OK；QI-1 反而给出反例）。 |
| **外部文献** | **0 条**。本文件不引任何 arXiv/DOI；所有断言只用仓内 `file:line`，故 `AGENTS.md §3` 的引用要求**空满足**。 |
| **未做** | 未跑任何 MC / 采样 / 数值作业（`无 MC`）；未改任何既有文件（只新建本文件）；未 `git commit`/`push`。 |

---

## §6 诚实边界

1. **本轮不证任何新定理。** `[T]` 只用在两类东西上：代数恒等（`G″=Δ_dyn+W_full`、QI-1、`(M)⟺Δ_dyn≥a²Δ₀²`）与仓内已证定理（R6-a、L5 强制引理、BL 单向、Lean `H1_does_not_imply_atom`）。**No QED。Clay 0%。**
2. **`(1a)⇒(2Δ)` 的 `[NG]` 依赖一条建模前提**：把"典型场曲率"读成"较低井的曲率"，并令二次井的极点质量 = 该井曲率。前者由 `ch18:182–183`（"the curvature at the conditioned point **is** the stability of that point"）支持，后者是标准对应但**仓内未逐字立过**；若改读"较重井"，该见证失效，但也不产生 `[T]`（它会退化成 §3(a) 的载体歧义）。**故本条按最保守规则记 `[NG]` + 本限定，不记 `[T]`。**
3. **`(3)⇒(2Δ)`/`(3)⇒(1a)` 的 2D 见证跨了载体**：`f₁²=1` 用的是 `O₂=½trU`（`ATOM2.md:78–81`），`Δ_dyn` 的数值是 class-angle 坐标下的（`N04-2D.md:40`）。**这正是 §3(g) 的错配本身**，不是论证漏洞；但因此 2D 只作**同型见证**，主见证用测度级的 `ν=½δ₀+½δ₂`（无模型假设）。
4. **两处 `t=0.5` 的 2D 读数不同**：`H1-FINAL.md:96` 写 `−26.478`，`N04-2D.md:53` 在 `φ=1.5` 写 `−21.8474`。前者**未给 `φ`**，无法逐点核对 ⇒ 记 `[E]`/`[M]` 口径差，**不改任何判决**。
5. **`(3)` 的 `⟹` 反例是数值的（`[E]`）**：`QI-BRACKET.md:36–38` 的 `L=8, ϑ=−1.74, w₁=1.05e−1`。解析反例我给了 `ν = ½δ₀+½δ₂` 一类的构造，但**没有**写成"对所有 `ν` 的解析命题"并机检；若要升 `[T]`，需补一条 `Δ>0 ⇒ ϑ<w₀≤w₁` 的显式三点测度族。
6. **台账编号是建议值**，不是落盘事实：`sec05-ledger.tex` 现为 38 条，本轮**未改**该文件（任务纪律：不改既有文件）。`(7/7)`、`quantifier scope (2/2)`、`Scope error (6/6)` 三处都需**父代理**决定是否落账与如何重排 `(n/N)`。
7. **位置图未更新**：`AGENTS.md §6` 要求结论变化时更新 `docs/position-current.html`；本代理被限定"不改既有文件"，故**未动**。若父代理采纳 §4.2 的替换文本，该文件与本文件的结论都属于**触发条件 3（靶的定义变了）**，须同步落盘。
8. **一条不成立的替代读法**：有人会争辩"四版只是同一对象的四个**面**，等价指'靶相同'而非'命题相同'"。本轮**部分承认**这一点：`ch13:34` 的"same object"在"同一靶"读法下可成立。但 `sec04:248` 用的符号是 `\iff`（逻辑等价），`sec04:250` 用的是 "i.e."（同一物），`H1-FINAL:20` 用的是"等价" —— **三处都是命题层的断言，不是靶层的**。故该争辩不改变本轮的 `[NG]`。

---

*判词：四版不是逐字同一命题；`0/4` 对等价；唯一 `[T]` 是恒等换名（零信息）。真身是 `O=¼trF²` 的 KL 测度上的原子 + 隔离，`H1` 只是它的天花板。`(3)` 的 `⟺` 与 `(2)` 内部的 `⟺` 各有一支为假（本轮新）。"收敛到 `H1`"须改为"收敛到 `H*`，`H1` 只占天花板一半"。No QED；Clay 0%。*

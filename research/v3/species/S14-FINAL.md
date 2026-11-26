# S14-FINAL — `(★)` 的防撞名终审：它与 `T01:35` / `H1` / `l2a` / `N09-1` 的四条关系

> 2026-11-25 · 子代理 `s14-fixture` · **只写本文件** · **未 git commit** · 未动 `docs/`（含 `position-current.html`：受"只写一个文件"约束）
> 被审对象（`(★)`，逐字）：
> `Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E_{μ(·|U′)}⟨PᵀH_jP v,v⟩ + W_j(v) − δ_W`，uniform in `(j,U′,a)` ⟺ `G″_j ≥ δ_W > 0`。
> 逐字读过（`read`/`Get-Content` 带行号）：`docs/paper-v2/ch08-immune.tex:6–30` · `ch18-scalar-h1.tex:28–67,75–88` · `app.tex:179–200` · `ch05-theorems.tex` 之 `:230–240`（经 `N05:24–37` 逐字转引并核 label）
> · `species/N05-SPLIT.md` 全文 · `Q02-H1I2.md` 全文 · `H1-FINAL.md` 全文 · `L2A-BRIDGE.md:1–120` · `L3A-STAR.md` 全文 · `L3C-STARSTRONG.md` 全文 · `N09-RESTATE.md:44–55,137` · `LOOP-RUNBOOK.md:11251–11315`。
> 标签：`[T]` 定理/代数恒等 · `[C]` 条件 · `[NG]` no-go/错判 · `[O]` 精确缺失步 · `[E]` 机检/数值 · `[M]` 未做的一步 · `[J]` 判断。**No QED；Clay 0%。**

---

## §0 判词（四条关系 + 三条判死线，先给）

| # | 问 | **裁定** | 级 |
|---|---|---|---|
| 1 | `(★)` 与 `T01:35` 的 `G″ = Δ_dyn + W_geo + X` 是同一个东西吗 | **是，且不止"同一个"：`(★)` 就是该恒等式的移项版**——`(★)` 与 `ch08-immune.tex:11–17` 的 `(M)` 逐字同式（仅 `(L,U′)`→`(j,U′,a)`），而 `ch18:30–37` 逐字把它写成 `Var − E[H] − W_full = −G″(0)` 并自判 "an identity with no content beyond that display"。⇒ `(★)` 的 `⟺` 是**一次移项**，零内容 | `[T]` |
| 2 | `(★)` 与 `H1` 是同一个吗 | **是同一族的具体成员**：`H1(δ) ⟺ G″(0) ≥ δ`（`ch18:43–44`，量词 "pointwise, mode by mode, for every `(β,U′)`"）；`(★)` 加 "uniform in `(j,U′,a)`" = 取**一致 `δ`** 的那个成员，即 `app.tex:182–192` 的 `∃Δ₀>0 ∀(L,a)∀v` 形（`δ := a²Δ₀²+W_full`，逐字"the two displays are the same inequality"）。⇒ `(★)` = `H1` 的 `δ`-一致形 = `Δ₀>0`；而 `H1 ⇏ m>0` | `[T]` |
| 3 | `(★)` 与 `l2a` 的乘法式是同一个吗 | **不是**。`l2a` 是**双因子**乘法桥 `f₁²≥c ∧ G″≥λ ⟹ m₁²≥cλ`（`L2A:104–113`，`BL-1` `f₁²/χ ≤ inf supp ν`、`BL-3` `χ≤1/λ ⟺ G″≥λ`）；`(★)` 只占其中**第二个合取项**（`λ`/天花板槽），**不含残数 `f₁²`**。⇒ `(★)` 是 `l2a` 的输入槽，不是 `l2a` | `[T]` |
| 4 | `(★)` 撞 `N09-1` 吗 | **撞。双角都撞**：① 若 `(★)` 的载体（`G″`，`∀v` 读法）过四模型 ⇒ 按 N09-1 必是 `m` 的函数，而 `ch18:81` 逐字 `G″_a = a²m_phys²` ⇒ 它就是靶的**量纲重标**；② 若不过 ⇒ 它就是被 2D（必要性）+ thin atom（充分性）逐出的那一族（`H1-FINAL §3.3`）⇒ 不是路 | `[T]`（条件于 M6 规格 + `L→∞`，`N09:137`）/`[J]` |
| — | `(★)` 的新增内容 | **∅**（§5）。唯一可称"新"的候选是"**单独约束 `Var`**"，而 `ch18:59–63` 逐字预先判死："`Var`, `E[H]` and `W_full` are **not separately constrained** by `H1` … **A proposal that bounds `Var` without controlling the sign of `G″(0)` is a claim about `G″(0)` in disguise, not an independent route.**" | `[T]` |
| (a) | 是 `T01:35`/`H1` 的换名 ⇒ 总攻点为空 | **触发（主判）。逐字见 §6(a)。** | `[T]` |
| (b) | 若"不是" ⇒ 指出新增内容 | **不触发**：`(★)` 是上述两个的换名，无可指认的新增内容；`(b)` 无可填 | `[T]` |
| (c) | 撞 `N09-1` ⇒ 明说 | **命中，且已明说**（§4）。 | `[T]` |

**一句话**：`(★)` **不是**第三个物种，它是 `ch08` 的 `(M)`（`eq:species`）——本项目给"一切诚实路线最终撞上的那一条不等式"起的名字——移项后的写法；论文自己在 `ch18:38` 写了 "Substituting the curvature form … and moving terms again yields …"。

---

## §1 关系 1：`(★)` ≡ `T01:35`（同一个东西，且只差一次移项）

### 1.1 逐字对照（`[T]`）

`docs/paper-v2/ch08-immune.tex:6–21`（`read` 逐字）：

```
 6: \subsection{The species}
 8: Every honest route terminates at the same inequality. Written in its most general form (the
 9: $(\blacklozenge)$ of \code{PROOF-JACOBIAN-SAVES.md}), the missing input is a
10: \textbf{single-sided upper bound on a conditioned variance}:
11: \begin{equation}
12: (\mathrm{M})\qquad
13: \operatorname{Var}_\mu\bigl(\langle v,\nabla S\rangle\mid0\bigr)
14: \ \le\ E_\mu[H_{vv}\mid0]+W_{\mathrm{full}}-\delta ,
15: \qquad\text{uniform in }(L,U'),
16: \label{eq:species}
17: \end{equation}
18: where $W_{\mathrm{full}}:=W_{\mathrm{geo}}+X$ is the full two-point remainder of
19: Chapter~\ref{sec:wgeo}; the display is then $G''\ge\delta$ in
20: curvature form. Equivalently $\operatorname{Var}(\xi_{k\to0})\le1/\lambda$ (H1), equivalently
21: $m_{\mathrm{eff}}^2(k\to0)\ge\lambda$, equivalently $G''(0)\ge\lambda$.
```

⇒ `(★)` 与 `(M)` **同式**：`⟨v,∇S⟩↔⟨Pᵀ∇S_j,v⟩`、`H_vv↔⟨PᵀH_jPv,v⟩`、`δ↔δ_W`、`(L,U′)↔(j,U′,a)`；而 `ch08:19` 自己就写了 "the display is then `G″ ≥ δ`"——`(★)` 的 `⟺` 是**引用**，不是推导。

### 1.2 它是恒等式的移项（`[T]`，两行双线性）

`ch18-scalar-h1.tex:28–37`（逐字）：

```
28: G''(0)=E[H_{vv}+\Phi''\mid\xi=0]-Var(G_v+\Phi'\mid\xi=0)
29:  =\Delta_{\mathrm{dyn}}+W_{\mathrm{geo}}+X , with \Delta_{\mathrm{dyn}}=E[H]-Var(G)
30: and X=-2Cov(G,\Phi'). Subtracting E[H]+W_{\mathrm{full}} from both sides,
31: W_{\mathrm{full}}:=W_{\mathrm{geo}}+X, gives an identity with no content beyond that display:
33:  Var_\mu(⟨v,∇S⟩|0) − E_\mu[H_{vv}|0] − W_full = − G''(0).
38: It uses no model, no uniformity, and no hypothesis that a gap exists.
40: Substituting the curvature form ... Var ≤ E[H]+W_full−δ --- and moving terms again yields the
43: H1(δ) ⟺ G''(0) ≥ δ
```

⇒ **`(★)` = `:38` 的恒等式代入自身、再移一次项**。机器核准：`max |G″ − (Δ_dyn+W_geo+X)| = 6.7e−16`（`N05 §3`，`w7/n05_split_audit.py`）——它是**定义**，不是不等式。

### 1.3 与 `T01` 判词的接口

`N05-SPLIT.md:52` 早已把 `(★)` 一字不差地记为 `ch08` 的 `(M)`：
> `(M)` | `Var_μ(⟨v,∇S⟩|0) ≤ E_μ[H_vv|0] + W_full − δ`，`(L,U′)`-一致 ⇒ `ch08:19` "the display is then `G″ ≥ δ`" | `ch08:11–19`（`eq:species`）

⇒ **判**：`(★)` 与 `T01:35` **是同一个东西**（`(★)` 是它的移项形），故 `(★)` 亦为**恒等式、零内容**。`[T]`

---

## §2 关系 2：`(★)` ≡ `H1` 的 `δ`-一致形（= `Δ₀>0`），而 `H1 ⇏ m>0`

### 2.1 `(★)` 落在 `H1` 的哪一格（`[T]`）

```
ch18:43–44   H1(δ) ⟺ G″(0) ≥ δ   （pointwise, mode by mode, for every (β,U′)）
ch18:55–58   族 {H1(δ)} 全序；"H1" without δ 就是 G″(0)>0
app:182–187  There is Δ₀>0 such that, uniformly in (L,a) and for every coarse transverse v at k→0,
             Var_μ(⟨v,∇V⟩|ξ=0) ≤ E_μ[H_vv|ξ=0] − a²Δ₀²
app:190–192  This is the (M) form ... with δ := a²Δ₀²+W_full; the two displays are the same inequality
```

⇒ `(★)`（"uniform in `(j,U′,a)`"的 `G″_j ≥ δ_W`）**就是** `app:182–187` 的 `∃δ∀` 形 = `H1(δ)` 的**一致版** = `Δ₀>0`。它不是 `H1` 之外的新命题，而是 `H1` 族里**唯一有物理名字的那一个成员**（`H1-FINAL §1.1`：`Δ_dyn ≥ a²Δ₀²` 不是第二定义，`δ_phys = a²Δ₀²+W_full`）。

### 2.2 因此继承 `H1` 的那条 `[T]` 否证（`[T]`，逐点）

`H1 ⇏ m>0`：thin atom `ν = wδ_ε + (1−w)δ_{s₂}`、`w = ε/2λ` ⇒ 用 `χ = ∫s⁻¹dν` 解出 `s₂ = 2λ(1−w)`、`χ = 1/λ` ⇒ **`H1` 真（取等）**，而 `inf supp ν = ε → 0`（**无隙**）。Lean `H1_does_not_imply_atom`，0 `sorry`；精确分数机核见 `Q02 §3.2`（`λ=1,ε=1e−6`：`χ=1.0000000000`、余量 `1.999998`）。

**限定（`AGENTS §2.1`，必须随结论引）**：该见证是**逐测度/逐 `(L,a)` 点**的；把 `(★)` 的**一致版**整句打掉还需"把反模型常值嵌入 `(L,a,β)` 族"，语料未做 ⇒ 一致版的反驳是 `[M]`，逐点版是 `[T]`（`Q02 §4.2(iv)`、`H1-FINAL §4`）。

---

## §3 关系 3：不是 `l2a` 的乘法式，只占它的第二个因子

`L2A-BRIDGE.md:104–118`（逐字，`[T]`）：

```
(BL-1)  f₁²/χ ≤ inf supp ν ≤ M₀/χ                 （ch02-squeeze.tex:51）
(BL-2)  (f₁²≥r₀) ∧ (χ≤1/λ)  ⟹  m₁² ≥ f₁²/χ ≥ r₀λ > 0   （ch02:199–201）
(BL-3)  χ ≤ 1/λ  ⟺  G″(0) ≥ λ                     （ch18 eq:scalar-rule）
⇒ 唯一成立的桥 = f₁² ≥ c ∧ G″ ≥ λ ⟹ m₁² ≥ cλ        （乘法桥；作为"残数⇒曲率缺陷"的蕴含不成立）
```

- `l2a` 的乘法式要**两件**：残数 `f₁² ≥ c`（谱/测度槽，`I1`）+ 曲率 `G″ ≥ λ`（`G″` 槽）。
- `(★)` 只给**第二件**，且给法是恒等式移项（§1）。它**一处也没提 `f₁²`**。
- ⇒ `(★)` **不是** `l2a`，而是 `l2a` 的**输入槽**。**唯一非空的邻接**是：`(★) + I1 ⇒ m₁² ≥ r₀·λ > 0`（`BL-2`，位置 = 重量 + 天花板，`Q02 §5.3`），而 `I1` 自身的窗只有 `β < 0.018394`（`[T] mod [KL]`）⇒ 这条路**不缺 `(★)`，缺 `I1` 的窗**（`l2b`：`(D/S)/(1−f₁²) = Θ(u^{−(T+1)}) → ∞`，方向反 + `β=0` Haar 上 `G″=+0.2051` 而 `f₁²=0` ⇒ 残数连必要性都没有）。

---

## §4 关系 4：撞 `N09-1`（双角都撞）

`N09-RESTATE.md:46–53`（`[T]`，条件于 M6 规格 + `L→∞`）：M6 上唯一趋于 `0` 的尺度是谱底 `m₆ := min(w,2π²/L²)` ⇒ "任何**过四模型**的量必是 `m` 的（重标）函数"。

`(★)` 的载体 `G″`（`∀v` 读法）两种命运，都使它不成为新对象：

| 角 | 假如 | 后果 | 级 |
|---|---|---|---|
| **①** | `G″` **过**四模型 | 四模型签名与 `m` 同：自由场 `0` / 2D `+2.0000` / M6 `inf_v = 0`（`Φ≡0 ⇒ W_full=0 ⇒ G″=Δ_dyn`，`b`-模 `H_vv=0`、`Var=0`）/ 单 plaquette 无 transfer 向 | ⇒ 按 `N09-1` 必是 `m` 的函数；且 `ch18:81` 逐字 `G″_a = a²m_phys²` ⇒ **靶的量纲重标**，不是新物种 | `[T]` |
| **②** | `G″` **不过**四模型 | 它正是被逐出的那一族 | 2D 假阴性（`m=3g²/8>0` 而 `Δ_dyn=0`）+ thin atom（`H1` 真而隙假）⇒ **不是路**，是反例载体（`H1-FINAL §3.2–3.3`） | `[T]`+`[E]` |

⇒ **`(★)` 撞 `N09-1`：明说。** 它既不是"非 `m` 的载体"（角 ①），也不能靠不过筛来救（角 ②）。这正是 `Q03` 的结论（靶 = `m>0`；`G″>0`/`H1`/`Δ₀>0` 不是书写）。

---

## §5 `(★)` 的新增内容：**∅**（唯一候选被论文逐字预先判死）

三项可能的"新"逐条排除：

| 候选新增 | 判 | 逐字/机制 | 级 |
|---|---|---|---|
| ① **单独约束 `Var`** | **否** | `ch18:59–63`（claimbox 第 2 条）："`Var`, `E[H]` and `W_full` are **not separately constrained** by `H1` beyond their appearance in `eq:scalar-identity`; the sign lives in the curvature. **A proposal that bounds `Var` without controlling the sign of `G″(0)` is a claim about `G″(0)` in disguise, not an independent route.**" `[T]` | |
| ② **"uniform in `(j,U′,a)`" 的量词升级** | **否** | 一致量词正是 `app:182–187` 的 `∃Δ₀∀(L,a)∀v`（= `H1(δ_phys)` = `Δ₀>0`），已判"是靶的换名、四模型 2/4、缺步 = `Δ_dyn` 的 `(L,a)`-一致正下界"（`H1-FINAL §3,§6(a)`；`ch13-machine.tex:49–50`）。**升级量词 ≠ 升级内容**：`G″` 的 `(L,a)`-一致性正是 `ch18:78–80` 逐字点名的那一格（"the open content is now (i) the sign and size of `G″(0)` … and (ii) its `(L,a)`-uniformity, and (ii) is the part that is not closed"） | `[T]`/`[O]` |
| ③ **`j`-指标（逐 species）** | **否** | `ch08` 的该式就叫 `eq:species`、章名 "The species"，其命题是"**每一条**诚实路线都终止于**同一条**不等式"（`ch08:8–9`）⇒ 逐 `j` 不产生新物种，只产生同一个 species 的多份拷贝 | `[T]` |

**⇒ `(★)` 的新增内容 = ∅。** 它精确等于"已入账的**唯一那一格**"（`Δ_dyn`/`G″` 的 `(L,a)`-一致正下界）的**改写**——把缺件重述一遍不是进展（`AGENTS §2.3`）。

### 5.1 旁证（不承重，但同向，两条独立）

- **`l3b`（机器）**：`I2 ∧ m₀ ≥ 0 ⇒ m ≥ δ > 0` 的 Lean 装配 `sorry=0`、三行、`g = δ`，而 **linter 直接报 `hstar` 未被引用** ⇒ 该装配里"`(★)`"这一格（能量谱写法 `∀y∈spec, y≠m₀ → y≥c`）**没被用到**（`LOOP-RUNBOOK:11262–11271`）。
- **`L3A`/`L3C`（谱写法）**：谱写法的 `(★)` 已定稿为 `inf_c Δ_c ≥ δΛ ≡ (a0)∧(b′)`，且 `(a0)` = **目标在 `0⁺⁺` 道的投影**、与目标"同一台机器"⇒ `L3` 是**记账不是归约**（`L3A §5(b)`、`L3C §5(a)`）；`(★-low)` 领头阶假（`w95c`）。

⇒ **两种写法、三条独立路径，同一结论**：`(★)` 不承重。

---

## §6 判死线（逐字）

### (a) **触发（主判）——"总攻点是空的"，逐字写清**

```
(★) 与 ch08-immune.tex:11–17 的 (M)（eq:species）逐字同式（仅 (L,U')→(j,U',a)）；
ch08:19 自己写 "the display is then G'' ≥ δ"；
ch18:30–37 把同一式写成 Var − E[H] − W_full = −G''(0)，并自判
          "an identity with no content beyond that display"（机检残余 6.7e−16）；
ch18:38–43 "Substituting the curvature form ... moving terms again yields H1(δ) ⟺ G''(0) ≥ δ"。
⇒ (★) 的 "⟺" 不是等价式的两端，是一个恒等式移了一次项。
⇒ (★) 是 T01:35 的移项形 + H1(δ) 的一致形（= Δ₀>0）。
⇒ **"全军的总攻点"是空的：它不是新物种，是同一 species 的第三种写法。**
      且它作为路线已按两条独立理由判死：H1 ⇏ m>0（thin atom，[T] 逐点）、
      Δ₀>0 是 m 的量纲重标（G″_a = a²m_phys²，[T]）。
```

### (b) **不触发**

`(★)` **是** `T01:35`/`H1` 的换名，无可指认的新增内容（§5 三项候选全被逐字排除）。若强行要"新增"，只能写"它把缺件重述了一遍"——那不是内容。

### (c) **命中 `N09-1`（明说）**

`(★)` 的载体 `G″` 要么过四模型（⇒ `N09-1` 判它是 `m` 的函数，且 `ch18:81` 给等式 `G″_a = a²m_phys²`），要么不过（⇒ 被 2D + thin atom 逐出）。两角都不产生新对象。`[T]`/`[J]`

### (d) 附带记账（不属判死线）

- **两个 `(★)` 不是同一个式子**：语料里"`(★)`"同时指 ① 本文审的微分式（`ch08` `(M)` = `eq:species`），② 谱式 `inf_c Δ_c ≥ δΛ`（`L3A §2`）。二者**不等价**（`Q02`：`H1` 与 `I2` 不可比），但**同槽**（"一切诚实路线撞上的那一条"）。**建议此后禁用裸 `(★)`**，改称 `(M)/eq:species`（微分）与 `(★-gap)`（谱）。
- 若后续仍要攻这一格，正确名字只有一个：**`Δ_dyn`（= `G″ − W_full`）的 `(L,a)`-一致正下界**，即 `ch13-machine.tex:49–50`；障碍已定位为 `N07` 的 `Faa` 变号（`Δ_dyn = −β·Faa + O(β²)`，`θ=π,a₀≳0.4` 处变号）与 `N03` 的 `(N/12)/L → 0`（`H1-FINAL §6(a)`）。

---

## §7 检索执行核验（`AGENTS.md §1` 四项）

| 项 | 记录 |
|---|---|
| 工具状态 | `grep` **OK**（2 次，均返回带行号正文，非 `ripgrep provider failure`）· `pwsh`/`Get-Content` **OK**（4 次，逐行带号）· `read` **OK**（8 次全文/大段）· `glob` **OK**（1 次）· `write` **OK**。**无 `ENOENT`、无超时、无 `TOOL-FAIL`。** |
| 检索式 1 | `T01:35\|N09-1\|l2a\|l3b` @ `research/v3` → **74 命中**（取 `LOOP-RUNBOOK:9630,9867–9927,10323–10336,11100,11251–11315,12048`；`H1-FINAL:143,145`；`L2C-SD:58–75`；`N09-RESTATE:14,46–55,137`；`P02-M6DYN:13,23`；`P03-M6REL:13,36`；`Q03-TARGET:20,82–91,168`；`FROZEN-POSITION:17,31,34`） |
| 检索式 2 | `species/S1*.md` → **0 命中**（`S14-FINAL.md` 为新建，未覆盖既有文件；工具 OK） |
| 逐字核（开卷，非转引） | `ch08-immune.tex:6–30`（`(M)`/`eq:species`/`eq:neqassoc2` 全文）· `ch18-scalar-h1.tex:28–67`（`eq:scalar-identity`、`eq:scalar-rule`、claimbox 三条）· `ch18:75–88`（`G″_a = a²m_phys²`、`(L,a)`-一致性未闭）· `app.tex:179–200`（`Δ₀` 块、"the two displays are the same inequality"、`D≤S−λ together with S≥λ`） |
| 未查到（采信） | **"`(★)` 作为独立于 `G″` 的不等式被证明/被使用"**：检索式 1 + 上述 4 处开卷 + `N05`/`Q02`/`H1-FINAL`/`L2A`/`L3A`/`L3C` 通读 ⇒ **0 处**（工具 OK）。反向证据 3 处：`ch18:31` 自判 "identity with no content"、`ch18:61–63` 预先判死"bounds `Var` ⇒ a claim about `G″(0)` in disguise"、`l3b` 的 linter 报 `hstar` 未用。 |

---

## §8 限定与自否证

- **限定 1（窗口）**：本文件的 `[T]` 全是代数/恒等层（无窗）；凡涉及 `m>0`、`I1`、`I2` 的判词只在 `β < 0.018394`（`[T] mod [KL]`，有限 `L`）与逐点口径成立，`β≈2.5` 只有单点 `[E]`（`H1-FINAL §2.3`）。
- **限定 2（`∀v` 口径）**：四模型表按 `app.tex:183` 逐字 `∀v` 读；取单模 `v=∂_a` 会在 M6 上造成"假阳性"错判（已撤回，`P02-M6DYN`、`H1-FINAL §3.2`）。
- **限定 3（`N09-1` 依赖）**：`N09-1` 的 `[T]` 条件于 M6 的模型规格 + `L→∞` 取值（`N09:137`）；M6 规格若改，§4 角 ① 需重算。
- **自否证**：若存在一条**不用 `G″` 符号**就能给 `Var` 单向上界（并因此定出 `G″` 的符号）的引理 ⇒ `ch18:61–63` 的判词被推翻，`(★)` 成为独立路线；该引理必须过自由场酸测（自由场 `f₁²=0`、`G″=0` 上失效）。语料内唯一候选是"负关联 / anti-FKG"`Cov(G,Φ′|0) ≤ 0`（`ch08:22–30`，`(eq:neqassoc2)`），`ch08:29–30` 逐字 "absent from every known dependence class" ⇒ 仍未过（`X` 符号不定，机检 144/144 点 `W_full < N/12`，min `−0.157`，`N05 §3.1`）。

---

*Status：`(★)` **不是**新物种。它 = `ch08` 的 `(M)`（`eq:species`）= `T01:35` 恒等式的一次移项 = `H1(δ)` 的 `δ`-一致形 = `app.tex:182–192` 的 `Δ₀>0`；`H1 ⇏ m>0`（thin atom，`[T]` 逐点）、`Δ₀>0` 四模型 2/4、且 `G″_a = a²m_phys²` ⇒ 撞 `N09-1`（双角）。`l2a` 是双因子乘法桥，`(★)` 只占其第二因子（无 `f₁²`）。**新增内容 = ∅**（`ch18:61–63` 逐字预先判死"bounds `Var`"这一路）。**判死线：(a) 触发（主判）——总攻点为空；(b) 不触发；(c) 命中 `N09-1`。** 唯一合法名字仍是"`Δ_dyn` 的 `(L,a)`-一致正下界"（`ch13-machine.tex:49–50`）。No QED；Clay 0%。*

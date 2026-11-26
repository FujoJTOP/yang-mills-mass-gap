# N08 · `app.tex` 中 `Δ₀` 的出处与地位

- 目标文件：`docs/paper-v2/app.tex`（817 行，52364 B，LastWriteTime 2026-11-25 00:58）
- 检索式：`grep -n "Delta_0|Δ₀|delta_0|\\Delta_\{?0\}?|Delta0"`，工具 OK，命中 **4** 行
- 检索式（上游）：`grep -rn` 于 `research/unconditional/` 同式，工具 OK，命中 53 行（见 §4）
- 标签：`[T]` 定理 / `[C]` 条件定理 / `[M]` 缺失步（`\LM`）/ `[NG]` no-go / `[E]` 数值

---

## 1. `Δ₀` 的全部出现处（逐字 + 行号）

| 行 | 逐字原文 | 作用 |
|---|---|---|
| **182** | `\LM{} The single open inequality, stated with no slack. There is $\Delta_0>0$ such that,` | **首次出现 = 存在量词引入**，带 `\LM`（Missing step）标签 |
| **187** | `E_\mu\!\big[H_{vv}\,\big|\,\xi=0\big]\;-\;a^{2}\Delta_0^{2},` | 不等式右端，`Δ₀` 出现在 `a²Δ₀²` |
| **191** | `$\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality.` | 定义 `δ` 时**消费** `Δ₀`（`δ` 才是被定义者，`Δ₀` 不是） |
| **194** | `$G''_a/a^{2}\ge\Delta_0^{2}$. Two dimension traps are recorded: the unconditioned budget` | 等价的**无量纲比形式** |

全部集中在 `§ app:sharp`（179–202，标题 `\subsection{The sharpest form of H1 and the exact missing estimate}`）。**无第 5 处**；`\section` 其它位置（含 verification map、notation key）**不含** `Δ₀`。

## 2. 上游：定义？假设？待证？

**无定义**：全文检索 `\Delta_0\s*:=` / `\Delta_0=` 命中 **0**（工具 OK）。`:=` 在 191 行只定义 `\delta`、`W_{\mathrm{full}}`。

**无"设"**：182 行不是假设句式（`设/Suppose/Assume`），是 `There is $\Delta_0>0$ such that ...` —— 即 **存在量化**：
```
∃ Δ₀>0  ∀(L,a)  ∀ coarse transverse v at k→0:
    Var_μ(⟨v,∇V⟩ | ξ=0)  ≤  E_μ[H_vv | ξ=0] − a²Δ₀²
```
`Δ₀` 是该量词式里的**约束变量**（bound variable），**不是**被假定的常数，**也**不是被构造出的常数（无 `:=`、无 `Δ₀=`、无上/下界估计、无"取 Δ₀ = ⋯"）。

标签判据（不是推断）：`main.tex:40` `\newcommand{\LM}{\labb{cM}{M}}`；`figures-extra.tex:90` `11.4/cM/\LM/{Missing step}/{one open inequality, no slack hidden}` ⇒ **`\LM` = Missing step = 待证**。同一节的结论句 `app.tex:246` `No claim in this appendix closes the gap. The reduction is machine-checked; the single missing estimate is stated with no slack` 与之自洽。

上游同式（独立命中，非转述）：`research/unconditional/H1-DIRECT-SPLIT.md:131` `∃Δ₀>0: ∀ U′, ∀ coarse transverse v ⊥ orbit, uniformly in (L,a):` → `:134` `⟺ G″_a(v)/a² ≥ Δ₀² (the mass gap in curvature language).`；同文 `:93` 表行 `| ratio (continuum) | G″_a/a² ≥ Δ₀², any Δ₀>0 | ... | **the target** |` ⇒ 上游同样把 `Δ₀>0` 标为 **the target（靶）**，而非输入。`H1-DIRECT-SPLIT.md:144` 还给出酸测：自由场 `0 ≤ 0 − a²Δ₀²` 不可满足 ⇒ 该式**非空**、正确排除自由场 `[C]`（条件式；在自由场上被证伪是设计目的，不是定理）。

## 3. 下游：谁依赖 `Δ₀>0`

1. `app.tex:191` `\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}` —— (M) 形式的 `δ` 只是同一不等式的改名，**依赖但无新增** `[M]`。
2. `app.tex:194` `G''_a/a^{2}\ge\Delta_0^{2}` —— 无量纲比形式，等价改写 `[M]`。
3. `app.tex:192-193` `D\le S-\lambda` **together with** `S\ge\lambda` —— 等价链，明确否掉 `D\le(1-\epsilon)S`（后者不足）`[M]`。
4. `app.tex:334` verification map：`Savvidy surviving margin $=$ H1 & \LM & SAVVIDY-MULTICELL-RG.md` —— 依赖方仍是 **`\LM`**，不是 `\LT`/`\LC`。
5. **条件定理链**：`app.tex:21` `clay_gap_conditional (T) (hH1) (hCP) (hWgeo) (hP3) (hno) : MassGap T` `[C]`；`app.tex:23` 明列 "The five analytic facts --- H1, (CP), W_geo≥0, P3, no-dissolution --- are"（假设）。此处吃的是 **H1**（`app.tex:790` `[H1] the ceiling $\chi\le1/\lambda$, uniform in $(L,a)$`），不是 `Δ₀`；`Δ₀` 只在 `H1-DIRECT-SPLIT.md:113` 以 `λ ≥ a²Δ₀²` 与 H1 挂钩。
6. **无** `[T]`、**无** `[LC]`、**无** `[LE]` 结论以 `Δ₀>0` 为条件：182–202 全段是"陈述 + 等价改写 + 两个量纲陷阱"，不含任何由它推出的定理。

## 4. 判：(b)，并附 (c) 的悬空成分

- **(b) 待证命题 —— 成立。** 精确陈述（`app.tex:182–189`，标签 `\LM`）：
  > `∃ Δ₀>0`，使对一切 `(L,a)` 一致、对每个 `k→0` 粗横向 `v`：
  > `Var_μ(⟨v,∇V⟩ | ξ=0) ≤ E_μ[H_vv | ξ=0] − a²Δ₀²`，`V = S_W + Φ`。
  这是 **H1 的最锐形式 / 待证命题**，其证明**不在本仓库**（`app.tex:248-250`：要证它的那台机器"does not exist in the present toolbox"）。
- **(a) 不成立（对 `Δ₀` 本身）。** `182` 行**不是**"设 `Δ₀>0`"；没有把结论塞进前提。**但 (a) 在另一层成立、必须记明**：条件定理 `clay_gap_conditional`（`:21`）的前提 `hH1` 就是 `χ ≤ 1/λ` 一致成立（`:790`），而上游 `H1-DIRECT-SPLIT.md:191`/`:113` 显示 `H1 ⟺ D ≤ S − a²Δ₀² ⟺ G″_a/a² ≥ Δ₀²` —— 即 **`H1` 这个前提与"存在质量隙"的曲率表述等价**。故：**逐字上 `Δ₀>0` 是待证命题；`H1` 作为定理前提才是那个"把结论当假设"的位置。** `[C]` 而非 `[T]`。
- **(c) 不成立（字面）。** `Δ₀` 并非"只在 190–191 出现一次"：实测 **4 次，182 / 187 / 191 / 194**，且 191、194 两处是**等价的改写**（有推导链：`\delta:=a²Δ₀²+W_full`；`G″_a/a²≥Δ₀²`），不是孤立悬空。
- **(c) 的悬空成分（须记）**：`Δ₀` **无定义式、无构造、无估计、无单位指定**（"dimensionless ratio" 的角色由 `a²Δ₀²` 承担）；它在被写下的那一刻起就是一个**未赋值的存在量词变量**，任何"取 `Δ₀ = f(β)`"式的引用都在原文里**没有依据**。这一点是 `Δ₀` 与 `λ`（`H1` 里的常数）之间的差别：`λ` 在 `:790` 同样未被赋值。`[M]`

## 5. 一句话

`Δ₀` 不是从任何推导里"来"的 —— 它是 `app.tex:182` 那个 `\LM` 待证不等式里**被存在量化出来的一个名字**；出现 4 次（182/187/191/194），2 次是等价改写；`Δ₀>0` 是**待证命题（靶）本身**，不是假设；真正"把结论当假设"的是条件定理 `:21` 的 `hH1`（`:790`），其内容与靶等价。

---
路径：`docs/paper-v2/app.tex`（L182–194 主体，L21/L334/L790 交叉引用）；上游 `research/unconditional/H1-DIRECT-SPLIT.md`（L93、L113、L131–134、L144、L191）。

# IIII1-WDEF · `w` 是什么、有没有上界

> 靶：`A = c₁ = a₁ ν̃σ([0,∞))`，`a₁ = 1`，`A = ∫ w dνσ`（νσ 概率）；门 `A ≤ 11.797261`（`HHHH1-AVAL.md:60`）。
> 纪律：只 `[T]`/`[NG]`；每结论带「下一手」；出处全部为实际 `read`/`grep` 命中；0 命中明写。
> `fax` = `research/v3/species/_ht7/faizal.txt`（32615 行）。

## §1 `w` 的定义（逐字）

**1.1 原文——同一句，两处**

- `fax:22293`（(6.21)）：`Ck = ∫₀^∞ e^{−t∆Ah} dν̃σ(t)`；`fax:22295–22297`：「for a finite positive Borel measure `dν̃σ` supported in a compact subset of (0, ∞) that differs from `dνσ` only by a **smooth nonnegative weight** arising from the quadratic part of the effective action. A dyadic partition of unity produces a decomposition with uniform exponential off-range bounds.」
- `fax:9463`（(C.1)）：`Cσ = ∫₀^∞ e^{−t∆Ah} dν̃σ(t)`；`fax:9466–9468`：「for a finite positive Borel measure `dν̃σ` supported in a compact subinterval of (0, ∞); this differs from the raw projector measure only by a **smooth nonnegative weight** coming from the quadratic effective action.」

⇒ 同一句在 (6.21) 与 (C.1) 各一次。⇒ `HHHH1-AVAL.md:13/65` 的「`smooth nonnegative weight` 仅在 `22296` 出现一次」**不成立**（第二处 `9467`）。`[T]`

**1.2 钉死**：`w` **不是符号、不是定义式**，是短语；语义 = `ν̃σ` 相对 `νσ` 的 **Radon–Nikodym 密度比**（「differs from `dνσ` only by … weight」，`22296`/`9467`），支在热时间 `(0,∞)`。fax 给 `w` 的全部信息 = 「光滑 + 非负 + 来自二次有效作用」。`[T]`

**1.3 不是指示函数、不是单位分解成员**：`{η_j}` 是热时间上的 FRD 分解（`fax:22322–22325`；`Σ_{j≥0}η_j ≡ 1` @`fax:22333`；`Σ_{j≤N}η_j ↑ 1_(0,∞)` @`fax:10159–10162`），作用在**测度**上（(6.24) `fax:22337`、(C.63) `fax:10085`）。`w` 与 `{η_j}` 无交集。`[T]`

**1.4 双重符号碰撞（会误导检索）**
- `ν̃σ` 有两支语义：§6/(C.1) 带权（`22293–22296`）；附录 A `(2.15)` 纯重标度 `dν̃σ(τ) = dν(σ²τ)`（`fax:30736`，`ν` 为 probability，`fax:30713–30714`）。`[T]`
- `w` 另有完全无关的含义：`w ∈ N` 板厚半宽，(8.20) `Col^(k)_w = {x : −w a_k ≤ x₀ ≤ w a_k}`（`fax:23038–23041`；`(2w+1)` 片 @`fax:23161`）。`[T]`

## §2 上界

**2.1 0 公式.** `w` 无显式表达式：`quadratic part` 于 fax 命中 2（`9468`/`22296`），均无公式。`[T]`

**2.2 0 上界.** fax 用到 `ν̃σ` 的地方只取「有限」：
- `fax:10144`：「holds with `c1 := a1 ν̃σ([0,∞))` … **the finiteness of c1 is part of the hypotheses**」
- `fax:10153`：「⟨ψ,Cσ^(j)ψ⟩ ≤ ∫η_j dν̃σ ≤ **ν̃σ[0,∞)**」
- `fax:10157`：「`c3 := ν̃σ([0,∞))`，independent of j」；`fax:22346`：「finiteness of `dν̃σ`」；`fax:10087`：`ν̃σ` 只被称「finite positive measure」

⇒ `w ≤ 1`、`w = O(1)`、`w` 无界：fax **三者皆未写**（`[T]` 关于「未写」）。检索 `w ≤|w≤|0 ≤ w` 于 fax 命中 4（`16992`/`18194` surface weight、`23041`/`23161` collar），**无一条是 `w` 的上界**。`[T]`

**2.3 唯一自动通道 = 概率归一，不是单位分解.** `fax:479`：「The bound `0 ≤ f ≤ 1` is **automatic if µf is a probability measure**; otherwise we may rescale µf to ensure ‖Π‖ ≤ 1.」而 `ν̃σ` 在三处**就是概率测度**：`fax:30713–30714`+`30736`（Def 2.1/(2.15)）、`fax:28255`+`28342`（Def 8.1/(8.16)，原句 "using that `ν̃σ` is a probability measure"）、`fax:31202`（Def 4.1）。`[T]`

## §3 判定：`A ≤ 11.797261`

- 恒等：`A = c₁ = a₁ ν̃σ([0,∞))`（`fax:10144`）= `a₁` × (ν̃σ 总质量)。`[T]`
- **支 (i) §6 未归一口径**（只写 finite positive：`22295`/`9466`/`10087`）：质量无值 ⇒ `A` 无上界 ⇒ `A ≤ 11.797` = **`[NG]`**。下一手：走 `α ≥ 2.327964` 支。
- **支 (ii) Def 2.1/8.1/4.1 口径**（`ν̃σ` 概率 ⇒ 质量 = 1）：`A = a₁`；且 `a₁ ≤ 1` 可取（(A.3) `fax:32061` 无前因子 ⇒ 前缀 1；`(C.68) fax:10146–10153` 的 semi-group contractivity 直接给 `0 ≤ Cσ^(j) ≤ 1`，亦不需 `a₁`）⇒ `A ≤ 1 ≤ 11.797261`，**余 11.797×** ⇒ **`[T]`**（在归一约定下）。下一手：核 §6 的 `dν̃σ` 是否即 `ν(σ²·)` 类推前；若是，把「`mass(ν̃σ)=1 ⇒ A=a₁≤1`」入约定表。
- 裁决：**门由归一放行，不由 `w` 的上界放行**；「`w ≤ 1`」提法作废，正确提法是 `⟨w⟩_{νσ} = 1`（`mass(ν̃σ)=1`）。`[T]`
- 警告：`a₁` 亦符号碰撞——`fax:10121–10124` 的 `a₁` 是 Combes–Thomas/Davies–Gaffney 仅存性常数，`a₁ = 1` 只在 `fax:30300`（另一章 `aR(β)` 语境）。`[T]`；`a₁` 在 §6 口径的具体数值 = `[NG]`（fax 未给数）。

## §4 检索核验（AGENTS §1）

| # | 检索式 / 工具 | 结果 | 命中 |
|---|---|---|---|
| 1 | `read` 22270–22329 / 9455–9494 / 10078–10172 / 22330–22355 / 30710–30745 / 28248–28349 / 470–499 / 23036–23047 | OK | 8/8 有效 |
| 2 | `nonnegative weight\|quadratic part\|Radon` @fax | OK | 6（`5487/9467/18045/18521/22296/27920`） |
| 3 | `partition of unity` @fax | OK | 8（`249/865/9477/10082/10187/22297/22322/22332`） |
| 4 | `11\.797\|11\.80` @fax | OK | **0 命中** |
| 5 | `ν̃_σ`（下划线写法）@fax | OK | **0 命中**（OCR 为波浪号+σ，非下划线） |
| 6 | `probability measure\|total mass\|mass of ν` @fax | OK | 53（含 `479/15479/19052/28255/28342/30713/31202`） |
| 7 | `ν̃σ [0, ∞)` @fax | OK | 2（`10141/10153`） |
| 8 | `w ≤\|w≤\|0 ≤ w\|w < 1\|weight bound` @fax | OK | 4（均与 `w` 上界无关：`16992/18194/23041/23161`） |
| 9 | `a1 =\|a₁ =\|exist constants a1` @fax | OK | 2（`10121`/`30300`） |

无 `TOOL-FAIL`；无百分比；未改既有文件；未 `commit`；零 MC。

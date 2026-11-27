# IIII2-PARTITION — `w` 的归类，与 `w ≤ 1` / `A ≤ 1` 裁决

> 二值等级 `[T]` / `[NG]`。`fax` = `research/v3/species/_ht7/faizal.txt`。
> 本轮不改任何既有文件、不 commit、零新 MC。

## 检索执行核验（AGENTS §1）

| 项 | 内容 |
|---|---|
| 工具状态 | `pwsh`(Get-Content 定行读) = OK；`Select-String` = OK；**无 TOOL-FAIL** |
| 检索式 1 | `Get-Content fax \| Select-Object -Skip 22284 -First 22` → 命中 22 行（22285–22306） |
| 检索式 2 | `... -Skip 21904 -First 36` → 命中 36 行（21905–21940） |
| 检索式 3 | `Select-String 'smooth nonnegative weight\|quadratic part\|dν̃\|w dν'` → 命中 33 行 |
| 检索式 4 | `Select-String 'partition of unity\|χ\|dyadic'` → 命中 40 行（截断显示） |
| 检索式 5 | `Select-String 'w\s*=\s*e\|w\(τ\)\|w\(t\)\s*=\|w ≤ 1\|0 ≤ w\|weight w'` → **13 行命中，其中 w 的显式公式 0 命中** |
| 检索式 6 | `Select-String 'normalized Gaussian measure\|total mass\|probability measure'` → 命中 20 行 |

## §0 判词

`[T]` 靶里被当成同一个 `w` 的其实是**三个互不相同的对象**，行号逐字如下：

1. `fax:22295–22297`（(6.21) 下方）：「for a finite positive Borel measure dν̃σ supported in a compact subset of (0, ∞) that differs from dνσ only by a smooth nonnegative weight arising from the quadratic part of the effective action. A dyadic partition of unity produces a decomposition with uniform exponential off-range bounds.」
   ⇒ `w` = **dν̃σ 相对 dνσ 的密度**，来源 = 有效作用量的**二次部分**；只给了 *smooth + nonnegative*。
2. `fax:21910–21913`：「where j∗ is the largest index with Rj∗ ≤ diam(B). … define the truncated (cumulant) expectation with respect to the normalized Gaussian measure with covariance Ckint by」 ⇒ 这里没有任何 `w`；只有**协方差 Ckint 的截断高斯**与累积量 `T⟨F1;…;Fn⟩`。
3. `fax:21922–21930`：`δp(U) = ½‖1 − Up‖²F ∈ [0, 2N]` (5.22)；`φ : [0,∞) → [0,∞)`，`φ(s)=s`（`0≤s≤1`）、`φ(s)=s²`（`s≥2`）、`φ″≥0`；`Rk(U) = exp(−λ Σ_{p⊂Λk} φ(δp(U)))` (5.23)。`fax:21935`：「multiplying dµσ by Rk preserves reflection positivity」。
   ⇒ `Rk` 是**乘在 dµσ 上的大场调节子**（plaquette 泛函），不是热时间测度的密度。

`[NG]` 因此「`w` 就是 `R_k`」这一归属不成立：`Rk = Rk(U)` 活在 plaquette 变量上，`dνσ` 活在热时间 `s ∈ [0,∞)`（`fax:15478–15481`、`19049–19054`）——类型不同。

## §1 两类候选（实际是三类）

| 候选 | 逐字证据 | 该对象自己的性质 | 是否是 `w` |
|---|---|---|---|
| (a) `w = R_k = exp(−λΣφ(δ_p))` | `fax:21930`；`φ ≥ 0`、`λ > 0` | `[T]` `0 < R_k ≤ 1`，等号 iff 全 `δ_p = 0` | `[NG]` 反例方式：`R_k` 乘 `dµσ`（`21935`），非 `dνσ` 的密度 |
| (b) `w` = 单位分解权 `χ_j` | `fax:865`「smooth nonnegative partition of unity {ψj}j≥0 on (0, ∞) with the usual dyadic supports」；`fax:9477`「a smooth partition of unity {ηj} on time」；`fax:22337` `Ck = Σ ηj(t) e^{−t∆Ah} dν̃σ(t)` | `[T]` `Σ_j η_j = 1`、`η_j ≥ 0` ⇒ 每个 `η_j ≤ 1` | `[NG]` `η_j` 是**积分内的时间局部化因子**，对算子和分解；`dν̃σ` 本身未变 |
| (c) `w = dν̃σ/dνσ`（密度/推前权） | `fax:22296`；`fax:9466–9468`「differs from the raw projector measure only by a smooth nonnegative weight coming from the quadratic effective action」；显式形态仅 `fax:30736` `dν̃σ(τ) = dν(σ²τ)`、`fax:32056` `ν̃σ(·) := ν(σ²·)` | `[T]` 只有 *smooth nonnegative*；`w` 的**显式公式 0 命中**（检索式 5） | `[T]` 这是唯一与 22296 逐字吻合的归属 |

附带逐字事实（裁决要用）：`fax:15479`「dνσ(s) on [0, ∞) with **total mass χσ(0) = 1**」；`fax:19052`「with a finite positive Borel measure νσ of **total mass fσ(0) = 1**」。
⇒ `νσ` 是**概率测度**（质量恒为 1）。

## §2 `w ≤ 1` 是否自动：逐条判

1. `[T]` `R_k ≤ 1` 本身**自动**：`φ: [0,∞)→[0,∞)`、`λ > 0`、`Σφ(δ_p) ≥ 0` ⇒ 指数 ≤ 0 ⇒ `R_k ≤ 1`（`fax:21927–21930`）。**但这条不给 `w` 任何信息**（对象不同，见 §1）。
2. `[T]` 单位分解权 `≤ 1` **自动**：`Σ_j η_j = 1`、`η_j ≥ 0` ⇒ 逐点 `η_j ≤ 1`（`fax:865`、`9477`）。**同样不给 `w`**：`η_j` 在 `fax:22337` 是求积因子，`dν̃σ` 不变。
3. `[NG]` **`w ≤ 1` 不自动**：`fax:22296` / `fax:9467` 对 `w` 只断言 *smooth nonnegative*；全仓检索 `w ≤ 1` / `0 ≤ w ≤ 1` / `w = e^{…}` 关于该权 **0 命中**（检索式 5，工具 OK）。
4. `[NG]` 若 `w` 取 `fax:30736` 的推前形态 `dν̃σ(τ) = dν(σ²τ)`，则 `w(τ) = σ²ρ(σ²τ)/ρ(τ)` 型雅可比，**没有 ≤ 1 的机制**（`σ` 任取 admissible，`fax:30731`）。

判：`w ≤ 1` = `[NG]`（不自动）。**但注意**：`A ≤ 1` 不需要逐点 `w ≤ 1`。

## §3 判定（`A ≤ 1` 能否成立）

1. `[T]` **精确等价**：`∫ dνσ = 1`（`fax:15479`、`19052`）⇒
   `A = ∫ w dνσ = E_{νσ}[w]`，故 `A ≤ 1 ⟺ E_{νσ}[w] ≤ 1`。这是**平均界**，不是逐点界。
2. `[T]` 若 (a) 或 (b) 成立（`w ≤ 1` 逐点），则 `A = E[w] ≤ 1`，门那一格过——**但 (a)(b) 的归属都被 §1 判 `[NG]`，此路不通**。
3. `[NG]` 若 `w` 是 (c) 的有效作用量二次部分权：fax **未给 `w` 的显式式子**（0 命中），故 `E_{νσ}[w] ≤ 1` **未证**，`A ≤ 1` 不能裁决为真。
4. `[T]`（退化支路，须标明）若额外取 `fax:30736` 的推前 `ν̃σ = ν(σ²·)`，推前保质量 ⇒ `A = ν̃σ((0,∞)) = ν((0,∞)) = 1` **恒等** ⇒ `A ≤ 1` 成立，但 `A ≡ 1` 使该格不含信息（不是对 `w ≤ 1` 的检验）。

**裁决**：`w ≤ 1` = `[NG]`（不自动）；`A ≤ 1` = `[NG]`（在原判据下未成立），唯一成立的形态是退化恒等 `A ≡ 1`。
下一手：①拿 `w` 的显式式子（`fax:22296` 只命名来源，未给公式）钉死 `E_{νσ}[w]`；②把 `A ≤ 1` 改成**平均**靶 `E_{νσ}[w] ≤ 1`，用 `νσ` 的 CM/Laplace 表示（`fax:15476–15481`）做分部积分/单调查表；③若走退化支路，须先说明 `A` 已被归一化吸收、门那一格应换量。

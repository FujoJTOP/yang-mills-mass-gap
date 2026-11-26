# Z11 · 裁一处：`δ ≥ κ/3` 的假是否落进 `I1`（洞 B 的范围）

> 2026-11-25 · 子代理 `z11-hole-b-scope` · 判官位。**只裁一处**：`DRESSED:94` 的 `δ≥κ/3`（`W05:66` 判假）是否落进 `I1`。
> 冲突双方：`Z01-P3ASSEMBLY.md:14`（`I1=[T] mod [KL]`，说不落）· `Z05-AUDITI1.md:63,89`（洞 B，说落）。
> 实读：`research/unconditional/P3-CYLINDER-RESIDUE.md:100–186`（含 `:151`、`:168` 全句）· `P3-DRESSED-MASS.md:80–174`（含 `:86–127` 全段）· `PROOF-P3-UNCONDITIONAL.md`（grep）· `research/v3/species/X01-L1C.md:40–80,113–118`（grep 逐字）· `Z01-P3ASSEMBLY.md`（grep，含 `:14,:148,:203–204,:271–279,:304`）· `Z05-AUDITI1.md:36–95` · `W05-TDECAY.md`（grep，含 `:66,:113,:123–125`）· `X02-THEOREMA.md`（grep，含 `:67–77,:110–113,:124`）· `W02-L1.md`（grep）· `LOOP-RUNBOOK.md`（grep，含 `:5500,:5707–5708,:5800–5807`）。
> 标签 `[T]/[C]/[NG]/[M]`。No QED。

---

## 0 · 判决（先给结论）

**`[T]`** —— **`D_∞` 不依赖 `δ`**（§2.4 三路逐字）；故 **洞 B 不落 `I1`，`z01` 对：`I1 = [T] mod [KL]`**。判死线 **(a)**。

**但 `z01` 的**依据**错位**：`CYL-RES:151/:168` 谈的是**谱隙** `δ_L=E_{2,L}−E_{1,L}`（`I2` 的 `δ`），不是 `DRESSED:94` 的**余项衰减率** `δ`；不能逐字援引它把 `DRESSED:94` 的 `δ` 摘掉。真正承重的依据是 `X01:47/:67`（`D_∞` 是**级数和**，"与 `δ` 无函数关系"）+ `DRESSED:116`（(4.1) 只用绝对收敛）。结论不动，依据换。

---

## 1 · `I1` 的逐字陈述：里面有没有 `δ`

`PROOF-P3-UNCONDITIONAL.md:15`（逐字）：`f₁²_L ≥ r₀(β) > 0`, **uniformly in `L`**（`t→∞` 在固定 `L` 先取；`L→∞` 不取）。
`:31`（逐字）：`f₁²_L ≥ r₀(β)=½c₁u^p>0`。
`Z01:203–204`：`f₁²_L = ν_L({m₁^{dr}}) = a₀D_∞ ≥ a₀(1−Cu) ≥ ½a₀ ≥ ½(2N)^{-2} > 0`；`inf_L f₁²_L ≥ ¼(2N)^{-2} =: r₀(β) > 0`。

**`δ` **不在** `I1` 的陈述里。** 逐项拆：

| `I1` 的成分 | 定义处 | 含 `δ`？ |
|---|---|---|
| `a₀ = (2N)^{-p}·s`，`s∈ℤ_{≥1}` | `DRESSED:36`、`:119` | 否：只含 Weyl/Haar 重耦合多重数与 `c_f>0` |
| `D_∞ := exp(Σ_{α≥0}(1−α)A(α))` | `X01:42` | 否：`A(α)` 是**附着跨度** `α` 的 Ursell 簇和，与任何率无关 |
| `D_∞ = 1+O(u) ≥ ½` | `DRESSED:116,:123–124`；`Z01:148` | 否：只用 `|log D_fin(T)| ≤ Cu` |
| `r₀` | `¼(2N)^{-2}` 或 `½c₁u^p`（`[M]` 两口径） | 否 |

⇒ `f₁² = a₀D_∞` 里**没有** `δ`；`D_∞` 的定义（`X01:42`）里**没有** `δ`。`δ` 在整条链里只出现在 (3.3) 的余项与 `η` 的界上（`X02:66–77`、`LOOP-RUNBOOK:5707–5708`）。

---

## 2 · 四问逐一

### 2.1 `DRESSED:94` 的 `δ≥κ/3` 用在哪一步、证什么

逐字（`P3-DRESSED-MASS.md:94–95`）：
> `C_L(t) = A^{dressed}e^{−m₁^{dressed}t}(1+O(u)) + O(e^{−m₁^{dressed}t}e^{−δt})`, `δ ≥ ½κ−O(u/a) ≥ ⅓κ` `(3.3)`
> **uniformly in `L`.** I.e. `e^{m₁^{dressed}t}C_L(t) → A^{dressed}` **with an exponentially small error at rate `δ≥κ/3`.**

它在文件里被消费于**两处**，都是"率/间距"，不是"值"：
1. `:95` 的**误差率**：`e^{m t}C_L(t) → A` 的**接近速率**；
2. `:107` 末句：`m₁^{dressed}(L)` 是谱底，"the next level at `≥1.5κ−O(u/a) > m₁^{dressed}`" —— **次级分离**。

而 `I1` 取的是 (3.1) 的**精确等式**（`DRESSED:86`：`Σ_{P∈𝔓_tube}w(P)=a₀u^{4T}Z_dec = A^{dressed}e^{−m₁^{dressed}t}D_fin(T)`，**无 `δ`**）+ Lemma 3.4 的**极限存在**（`DRESSED:99–101`，只要 `C(t)e^{mt}→A`）+ `D_fin→D_∞`。
⇒ **那个结论（定量余项率 `κ/3`）不进 `I1`；`I1` 至多消费"余项相对领头项可忽略"这一陈述（`∃δ>0`，分割固定后）。** 这正是 `Z01:14`、`:254`(S7)、`:304`(R5) 的做法：**撤 `κ/3`、换 `∃δ>0`**；`W05:113`(E1) 的勘误建议逐字同此（"撤 `δ≥κ/3`"、改 `|η|≤Cu+Be^{−δt}`）。

**洞 B 即使成立，方向对 `I1` 有利、不利两说均无负项**：`Z05:63` 说的同速率族（横向位移四壁管 `|P|=4(t/a)+4d`）按 `DRESSED:44`/`:69`、`FOURWALL:79,86`（"finite positive transverse sum"）**并进 `a₀` 的横向求和**——是**正因子**，对下界只增不减。故 `f₁²_L ≥ r₀` 这一**不等式**对洞 B 鲁棒（价值 `a₀` 的重标定，不是 `D_∞` 的污染）。

### 2.2 `CYL-RES:151/:168` 的逐字，与"系数是值、隙另算"

`P3-CYLINDER-RESIDUE.md:151`（逐字，全句）：
> `PROOF-P1P2` §4 already flags: the `O(u)` remainder bound `e^{m₁t}R_L=O(u)` is *"a `t`-uniform multiplicative bound, satisfied by any `λ₂/λ₁∈[1−Cu,1)` — hence it certifies the residue I1, not the isolation I2."* This is **consistent** with the lemma above and sharpens it: the `O(u)` bound fixes the **coefficient** `f₁²_L` (the value, clause (d)) but is **compatible with** `δ_L=E_{2,L}−E_{1,L}→0` as `L→∞`. That is precisely why (b) uses only `δ_L>0` *at finite `L`* (not a uniform `δ`): **the coefficient is a value, the gap is a separate object.** **The residue identification does not secretly smuggle in a gap.** `[T]`

`:168`（逐字）：
> **NOT needed** for the **finite-`L` residue floor** `f₁²_L≥r₀(β)`: at each finite `L`, compactness makes `δ_L>0` automatically, and the *value* `f₁²_L` is bounded by the polymer amplitude (clause (d)), **independent of how small `δ_L` is**. **The residue theorem (I1) stands without `I2`.** `[T]/[C]`

**意思**：`O(u)` 余项界定的是**系数**（一个值），它对 `δ_L→0` 也成立 —— 系数与间距是**两个对象**，故残数识别不暗带一个隙。

**能不能援引它把 `δ` 从 `I1` 里摘掉？**
- **能援引的**：`I1` 不需要**谱隙**（`I2`）这一层 —— `:147` 逐字："The cost grading is a statement about polymer **amplitudes**, not transfer-matrix **eigenvalue spacing** … it cannot and need not produce a uniform gap `δ` for `f₁²_L≥r₀`."
- **不能援引的**：把 `DRESSED:94` 的 `δ`（**余项衰减率**）说成"被 `:151/:168` 摘掉"。`:151/:168` 的 `δ_L` 是 `E_{2,L}−E_{1,L}`；`X02:124` 已把 `δ` 分名（`δ_L1` 质量率 / `δ_I2` 质量平方隙 / `δ` 余项）。**同名异物，逐字不对号。**
⇒ `z01` 的**依据错位**（应改用 `X01:47/:67` + `DRESSED:116`），**结论仍对**。

### 2.3 ★ `D_∞` 是否依赖 `δ` —— 三路，全部指向"不依赖"

**(i) 定义路（`[T]`）** `X01:42`（逐字）：`φ_dec := Σ_{α≥0}A(α)`, `log D_∞ := Σ_{α≥0}(1−α)A(α)`。`X01:47`（逐字）："(L1-c.1) 同时**证明 `D_∞` 存在**（右端是绝对收敛级数，不先取极限）"。⇒ `D_∞` 由 `A(α)` 级数**直接定义**，不经极限、不经任何率。

**(ii) 逐字否认路（决定性）** `X01:67`（逐字）：
> 若把 `−(−δT)` 写进定义，等价于定义 `D_∞:=\lim e^{δt}D_fin`，而 `D_fin` 只**收敛**到 `D_∞`、不指数增长（`(L1-c.1)` 里 `D_∞` 是级数和，**与 `δ` 无函数关系**）⇒ 该写法在 `T=0` 就不成立。

**(iii) 界路（`[T]`）** `DRESSED:116`（逐字）：`|\log D_fin(T)| ≤ C₁u+O(1)·u = Cu … uniformly in T,L`，"**Crucial:** (4.1) uses only **absolute convergence**, *not* termwise nonnegativity"；`:117` 同。⇒ `D_∞=1+O(u)≥½` 只吃绝对收敛（`DRESSED:123–124`；`Z01:148`）。

`δ` 在 `D_∞` 附近唯一出现处是 `ε_T` 的**界**：`X01:80` `|log D_fin(T)−log D_∞| ≤ Bū^T`，`δ:=ln(1/ū)/a`；且这个 `δ` 是 `X01` **自给**的（来自跨度⇒成本 + 绝对收敛），**不是** `DRESSED:94` 的 `κ/3`。`LOOP-RUNBOOK:5707–5708` 同构：`η=(D_fin/D_∞−1)+(R_hi/f₁²)e^{m₁^{dr}t}` ⇒ `|η|≤2Be^{−δt}` "**只需 `D_∞≥½`**" —— `δ` 进 `η` 的界，不进 `D_∞`。

⇒ **`D_∞` 不依赖 `δ`。** 洞 B 的假 `δ≥κ/3` 不进 `D_∞`，故不进 `f₁²_L=a₀D_∞`，故不落 `I1`。

### 2.4 判死线

**(a) 成立。** `D_∞` 不依赖 `δ` ⇒ **洞 B 不落 `I1`，`z01` 对：`I1 = [T] mod [KL]`**（`[KL]` 单一名义输入，`Y01-KL:118` 判 `[T]`）。(b) 不触发。(c) 不触发——语料对此**明确**（`X01:67` 逐字否认），不是"不明"。

---

## 3 · 影响面（不越界，只记一句）

- **不动的**：`I1` 的等级（`[T] mod [KL]`）、`r₀` 的下界链、`D_∞≥½`、`X01` 的 `(L1)`。
- **要动的（措辞级，非本轮裁项）**：`DRESSED:94–95` 的 `δ≥κ/3` 与 `:107` 末"next level `≥1.5κ`"仍按 `Z01:304`(R5)/`W05:113`(E1) 撤换为 `∃δ>0`（分割固定后）——它们属 `I2` 形状与 `E_{1,L}` 的账面，不进 `I1`。
- **`z01` 的引用要改**：`CYL-RES:151/:168` 只能支持"`I1` 不需谱隙"，不能支持"摘掉 `DRESSED:94` 的余项率 `δ`"。主依据换成 `X01:47/:67` + `DRESSED:116`。
- **`Z05:89`（洞 B 行"对 `I1` 的后果"）** 应改写为"对 `I2`/`E_{1,L}` 账面的后果"，不列 `I1`。

---

## 4 · 检索核验（AGENTS §1）

| # | 工具 | 检索式 | 状态 | 命中 |
|---|---|---|---|---|
| 1 | `glob` | `research/v3/**/*.md` / `research/v3/species/*.md` | OK | 362 / 250（截断列表，已列前 100） |
| 2 | `glob` | `research/v3/species/Z1*.md` | OK | **0**（本文件此前不存在） |
| 3 | `grep` | `CYL-RES` @ `research` | OK | 15 |
| 4 | `grep` | `DRESSED` @ `research`（`*.md`） | OK | 208 |
| 5 | `grep` | `D_∞\|D_\\infty\|D_\{?\\?infty` @ `research` | OK | 81（关键：`X01:41–42,47,67,80`；`X02:66–77,110`；`Z01:145–148,203–204`；`LOOP-RUNBOOK:5500,5707,5800–5807`） |
| 6 | `grep` | `I1\|residue` @ `PROOF-P3-UNCONDITIONAL.md` | OK | 8（`:15,:31`） |
| 7 | `read` | `P3-CYLINDER-RESIDUE.md:100–186` | OK | 87 行（`:151`、`:168` 逐字已录） |
| 8 | `read` | `P3-DRESSED-MASS.md:80–174` | OK | 95 行（`:94–95`、`:107`、`:115–125` 逐字已录） |
| 9 | `read` | `Z05-AUDITI1.md:36–95` | OK | 60 行（`:63` 洞 B、`:89` 定位表） |

**无 TOOL-FAIL。** 未逐字核：`X01-L1C.md`、`W05-TDECAY.md`、`X02-THEOREMA.md`、`W02-L1.md`、`LOOP-RUNBOOK.md` 的**全文**（只经 `grep` 命中行取逐字，均含行号）；`P3-AUDIT-CORRECTED.md`、`P3-FOURWALL-DERIVATION.md`、`P3-FULLREP-KP.md` 未直读（本轮不裁它们）。

*Refs.* 内部：`research/unconditional/P3-CYLINDER-RESIDUE.md:147,151,168` · `P3-DRESSED-MASS.md:44,69,86,94–95,99–101,107,116–117,123–125` · `PROOF-P3-UNCONDITIONAL.md:15,31` · `research/v3/species/X01-L1C.md:41–47,67,80` · `Z01-P3ASSEMBLY.md:14,148,203–204,254,271–279,304` · `Z05-AUDITI1.md:60–63,88–89` · `W05-TDECAY.md:66,113,123–125` · `X02-THEOREMA.md:66–77,110–113,124` · `W02-L1.md:22–25,53,78` · `LOOP-RUNBOOK.md:5500,5707–5708,5800–5807,6307`。外部（仅转引，未新核）：Osterwalder–Seiler Ann. Phys. **110** (1978) 440；Kotecký–Preiss CMP **103** (1986) 491；Seiler LNP **159** (1982)。

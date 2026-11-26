# Y06-P3DOWN — `MASTER:42` 到底降了什么：降的**正是** `f₁² ≥ r₀`（= `I1` 本身）；⇒ 按总册当下口径 `I1` = `[C]`，`I2` 链第一环**不在** `[T]`（判 (b)）；但三条降级理由**全部已有补证文件**（判 (c) 同时触发）⇒ 回填后 `I1` 可回 `[T] mod [KL]`（`[KL]=[T]`）

> 2026-11-25 · 任务：核 `y02` 报的「`MASTER:42` 已把 `f₁² ≥ r₀ > 0` 降为 `[条件/待重证]`」。
> 只写本文件；未 git commit。标签 `[T]/[C]/[NG]/[M]/[O]/[J]`。

---

## §0 检索执行记录（`AGENTS.md §1` 四要素）

| # | 检索式 | 工具 | 工具状态 | 命中 |
|---|---|---|---|---|
| 1 | `Get-Content docs/YANG-MILLS-MASTER.md` 取 `[35..55]` 带行号 | `pwsh` | OK（有效响应） | 文件 241 行；`:42` 逐字取出 |
| 2 | `Select-String -Path docs/YANG-MILLS-MASTER.md -Pattern 'P3\|残数'` | `pwsh` | OK | 11（`:35` 小节标题、`:41`、`:42`、`:152`、`:218` 等） |
| 3 | `grep "f₁² ≥ r₀"` @ `D:\Dev\math` `*.md` | `grep` | OK（250/298，全量落盘 spill） | 298 |
| 4 | `grep "P3-AUDIT-CORRECTED\|待重证\|MASTER:42\|不再是 \[T\]"` @ `research/v3/species/*.md` | `grep` | OK | 12 |
| 5 | `read` `D:\Dev\math\research\unconditional\P3-AUDIT-CORRECTED.md` `:60–133` | `pwsh` | OK | 133 行全文关键段 |
| 6 | `read` `research/v3/species/U07-CPSTRONG.md` 全 148 行；`X02-THEOREMA.md:1–60`；`Y01-KL.md:1–60`；`X01-L1C.md:1–40`；`Y05-P3GEN.md:1–30` | `pwsh` | OK | — |
| 7 | `grep "\[R\]\|\[C\]\|uniform\|1−O\(u\)\|resum"` @ `research/unconditional/P3-DRESSED-MASS.md` | `grep` | OK | 23 |

**无 TOOL-FAIL。** 未逐字核：`P3-FULLREP-KP.md` 证明主体（只经 `U07 §2` 与 `P3-AUDIT-CORRECTED:57–58` 转引 `K≤10³N`、`Δ=20`）。

---

## §1 `MASTER:42` 逐字（`docs/YANG-MILLS-MASTER.md:42`）

**这是 SZZ 给不了的一半（SZZ 只给隙的相关系数上界；P3 想给残数下界）。** ⚠ **2026-11-24 审计更正**：`PROOF-P3-UNCONDITIONAL.md` 曾称 KP 引理"完整写出并闭合"——**不成立**：最小管道几何应为 `m₁=4ln(1/u)/a`（四壁管道，非"每片 1 个"；factor-4，Gupta §14.2.2 / Münster）、全表示 KP 未控、dressed-mass 重求和未证。**残数结论 `f₁²>0` 存活**（因子 4 只重标定质量），但证明**不完整** → P3 标 **[条件/待重证]**，不再是 [T]。见 `P3-AUDIT-CORRECTED.md`。

同文件 `:35` 小节标题逐字：`### 2.1 P3 强耦合残数定理（[条件/待重证]，2026-11-24 审计）`。

### §1.1 它降的是**哪一个断言** —— 逐项排除

| 候选 | 是否被 `:42` 降？ | 依据（`MASTER:42` 逐字） |
|---|---|---|
| **`f₁² ≥ r₀ > 0`（= `I1`，P3 残数定理本身）** | **是** | "`PROOF-P3-UNCONDITIONAL.md` 曾称 KP 引理'完整写出并闭合'——**不成立**" + "P3 标 **[条件/待重证]**，不再是 [T]" |
| `(CP)`（`= I1 ∧ I2`） | **否**（`:42` 未提）；`(CP)` 本来就不是 `[T]` | `(CP)` 的 `I2` 半在 `β<1/48` 上从未证（`U07 §5(a)`；`IMP-ARROWS.md:128`） |
| `r₀>0` 这个**量**的存在性 | **否** | "因子 4 只重标定质量"、"残数结论 `f₁²>0` **存活**" ⇒ `r₀>0` 未被否证 |
| `⟨Ω|O|𝔾⟩ ≠ 0`（幅度地板） | **否**（间接：属被降证明的 §4 一环） | `:42` 只点名 "dressed-mass 重求和未证" |
| `ν ≥ 0`（`[KL]`） | **否**（`:42` 不涉） | `[KL]` 是 OS'78 标准定理，`y01 §2` 判 `[T]` |

**判词**：降的是**证明完整性标签**，被降对象**逐字就是 P3 = `f₁²>0`**，亦即任务里的 `I1`。
"存活"指的是**断言未被反驳**，不是"仍是 `[T]`"——`:42` 明写 "**不再是 [T]**"。

### §1.2 降级的三条理由（`:42` 逐字列出）

1. 最小管道几何应为 `m₁=4ln(1/u)/a`（四壁管道，非"每片 1 个"；factor-4，Gupta §14.2.2 / Münster）；
2. 全表示 KP 未控；
3. dressed-mass 重求和未证。

`P3-AUDIT-CORRECTED.md:131`（同源判词，逐字）："**Proof as written … NOT valid** … Label **[条件 / 待重证]** (matching the 2026-11-24 header already added). `[J]`"；

`:131` 末："What does **not** survive is the label 'unconditional/closed': the correct status is **[C] conditional on the dressed-mass resummation (§1.3) + the standard full-rep KP (§1.2)**."

⇒ `MASTER:42` = `P3-AUDIT-CORRECTED` 的摘要，两者**一致，无冲突**。

---

## §2 降级理由 vs `u07` / `x02`：是否一致

### §2.1 `u07` 的判词（`U07-CPSTRONG.md`）

- `:9` 判 **(a)**：「`f₁²>0` 在 `β<1/48`（固定 `a`、柱面、`L`-一致、有限 `L`）**已经 `[T]`**」。
- `:71` §3：「`r₀>0` 是 `[T]` **mod OS/KL**」，并"严格等级（**不得升格**）"列出：算术部分（`c₁>0`,`u>0`）`[T]` 恒等；**有内容的是 `f₁² ≥ r₀` 这条不等式**，`[T] mod OS/KL`。
- `:101` 链图把 `[KP]`、`[重求和]` 两格都标 `[T]`。

### §2.2 `x02` 的判词（`X02-THEOREMA.md`）

- `:9`：「链的**唯一**载重物是 `(L1)`（`[C]`，其 L1-c 是 `[O]`）；Theorem A 贡献**零动力学**，只做'率 ⇒ 隔离'的翻译。」
- `:36` 载重行：**"`⇐` 是全部内容，用 `ν≥0`"**；`w02` 用的正是 `⇐`（先证率，再得隔离），"方向正确"。
- `:42`：`ch05:435` 是等价改写，可作桥**只在 `⇐` 方向且 SharpDecay 已被独立证明**时。

### §2.3 判定：**两者不冲突，但 `u07` 与 `MASTER:42` 冲突**

| 比较 | 结论 | 依据 |
|---|---|---|
| `u07` vs `x02` | **一致（不冲突）** | 二者谈**不同环**：`u07` 谈 `I1`（残数地板 `f₁²≥r₀`）；`x02` 谈 `I2` 的翻译（率⇒隔离）。`x02` 明说 Theorem A "贡献**零动力学**"、不生产 `I1`；`u07 §4` 明说 `f₁²>0` **不需要** `(CP)`、不需要 `Z_O`、不需要 `I2`。⇒ `x02` 的 `[C]` **不落在 `I1` 上**。 |
| `u07` vs `MASTER:42` | **冲突** | `u07` 判 `[T] mod OS/KL`；`MASTER:42` 判 `[条件/待重证]`、明写"不再是 [T]"。`u07` 全文只引 `P3-AUDIT-CORRECTED` 的 `:51`、`:87`（SZZ 引用与 `p=0` 记账），**未处理审计的 §1.3 `[R]`**，也未引 `P3-DRESSED-MASS`。⇒ `u07` 的 `[T]` 是**越过审计更正的标签**，不能覆盖 `MASTER:42`。 |
| `x02` 的 `(L1)` `[C]` vs `x01` | `x01` **已修** | `X01-L1C.md:5`：「`(L1)` ⇒ `[T]`（在 KP 盘内），`I2` 的路径因此闭合（余 `[KL]` 与两处记号勘误）」。⇒ `x02` 的 `(L1)[C]` 已被 `x01` 补成 `[T]`。 |

---

## §3 关键：降的是**哪一步** —— 逐条判 `[T]` 等级

任务列出的三个候选步，逐个判（**域统一为 `β<β₁`、固定 `a`、柱面 `𝒞_L`、不取 `L→∞`**）：

| 步 | 断言 | 等级（当下，按总册口径） | 依据 |
|---|---|---|---|
| **①** | `r₀ = ½c₁u^p > 0`（组合/代数） | **`[T]`** | 纯算术：`c₁>0`（SU(N) 正群常数）、`u(β)>0`、`p` 有限。`U07:73` 逐字："`r₀>0` 的**算术**部分（`c₁>0`, `u>0`）是 `[T]` 恒等，平凡"。`p` 的取值只移前因子（`ch05:28`），`P3-AUDIT-CORRECTED:87` 取 `p=0` 于归一 A1g plaquette——**记账级不齐，不改正性**。 |
| **②** | `⟨Ω\|O\|𝔾⟩ ≠ 0` / 等价地幅度地板 `a₀>0` 与 `inf_L A^{dressed} ≥ ¼(2N)^{-2}>0` | **`[T]` mod `[KL]`，但挂在被降的证明上 ⇒ 当下 `[C]`** | `P3-DRESSED-MASS.md:16` 逐字："Decorations are `O(u)`-small … so `A = a₀(1±Cu) ≥ (1−Cu)a₀ > 0`; `a₀>0` uses only `c_f>0` and the singlet recoupling `f⊗f̄→1` (multiplicity `1`)"，标 `[T]`；`:170`："(d) uniform positive amplitude floor — PASS `[T]`"。**但**这正是不在 `PROOF-P3-UNCONDITIONAL.md §4` 内、被审计判 `[R]` 的那一环（`P3-AUDIT-CORRECTED §1.3`:60 逐字："**Proving** this uniformly (items 1–4) is the load-bearing missing step; the file's §4 does not contain it. `[R]`"）。⇒ 证明在**伴文**里成立，在**主证明**里未写入 ⇒ `MASTER:42` 仍判不完整。 |
| **③** | `ν ≥ 0`（`[KL]`：柱面 KL/OS 表示、`dμ_O≥0`） | **`[T]`** | `Y01-KL.md:56` 判："既不是'框架假设（不可证）'，也不是'数学缺口（未解）'，而是'**标准定理、当时未在 P3 文件内写出**'。`[T]`"。`P3-DRESSED-MASS.md:18` 同判："The single remaining **analytic input** is (c)'s KL/OS … is standard `[T]` (OS'78) but not written in the P3 file." |

**哪个不是被降的对象**：`ν≥0`（`[KL]`）**完全不在 `:42` 的清单里**；`r₀>0` 的算术部分也**不在**。`:42` 降的清单是**三条**：四壁几何、全表示 KP、dressed-mass 重求和 —— 都属 `②` 那一步（`f₁² ≥ r₀` 的**证明**），不是 `①`、不是 `③`。

### §3.1 三条理由的补证现状（这是 (c) 的核心）

| 降级理由 | 伴文 | 判 | 逐字证据 |
|---|---|---|---|
| 四壁几何 `m₁=4ln(1/u)/a` | `P3-FOURWALL-DERIVATION.md §1.1–1.5` | **`[T]` 已补** | `:18`："The **four-wall minimal-surface geometry is now a theorem [T]**"；`P3-AUDIT-CORRECTED §2` Lemma 2.1 亦标 `[T]` |
| 全表示 KP 未控 | `P3-FULLREP-KP.md` | **`[T]`（标准，Seiler）** | `P3-AUDIT-CORRECTED §2`:57："**Full-rep KP estimate** (factorial domination) — *standard, Seiler LNP 159* `[T]`"；`U07:60` 转引 `K≤10³N`（`β≤1/48`）、`Δ=20` |
| dressed-mass 重求和未证 | `P3-DRESSED-MASS.md` | **(a)(b)(d) `[T]`；(c) ⇒ 单一命名输入 `[KL]`=`[T]`** | `:13,14,16` 三项 PASS `[T]`；`:18`："closed for (a),(b),(d)"；`:170` 全 PASS |

⇒ **三条理由全部已有可引文件**；剩下的唯一命名输入 `[KL]` 由 `y01` 判 `[T]`（标准 OS'78）。
⇒ 但这**不等于** `MASTER:42` 已被否证：`:42` 说的"证明不完整"是就**主证明文件** `PROOF-P3-UNCONDITIONAL.md` 而言，而该文件**没有**回填这些伴文（`P3-DRESSED-MASS:3` 明写 "This file does **not** edit the proof or the paper"）。总册口径 = 主证明口径。

---

## §4 判死线

### (b) **触发** —— 降的正是 `f₁² ≥ r₀`，`I1` 不在 `[T]`

逐字对照，无歧义：`MASTER:42` 的 `[条件/待重证]` 落在 **P3** 上，而 **P3 就是 `f₁²>0`**（`MASTER:35` 标题逐字："### 2.1 P3 强耦合残数定理（[条件/待重证]，2026-11-24 审计）"）。任务把 `I1` 定义为 `f₁² ≥ r₀`，`U07 §4` 表第 0 行逐字："`f₁² ≥ r₀(β)` —— **就是 P3 的结论本身**（不是一步推理）"。

```
(I1)[C]  ←  MASTER:42 + P3-AUDIT-CORRECTED:131（[条件/待重证]）
   ↓
(L1)[T]  ←  X01-L1C:5（L1-c 已写成可核一页，[C]→[T]）
   ↓
X02 代数（log D → η）        [T]（X02 §2.1，初等）
   ↓
Thm A 的 ⇐（用 ν≥0）         [T]（X02 §1.1，REPORT.md:95–98,:100,:110）
   ↓
I2                           [C] —— 载重全在 (L1) 与 I1；I2 本身在 β<1/48 仍 ✗
```

**⇒ 第一环 `I1` = `[C]`**；`I2` 链因此**不在 `[T]`**。按任务 (b) 的措辞明说：
**「`I2` 路径在 `[C]` 级闭合」**——不是 `[T]` 级。`Clay` 进度 **0%** 不变（`U07 §6`、`X02:7` 同判）。

### (c) **同时触发** —— 重叠处已定位，且**已被那几步补上（除回填这一行政步骤）**

| 重叠处 | 与哪一步重叠 | 是否已补 | 判定 |
|---|---|---|---|
| dressed-mass 重求和 `[R]` | `P3-DRESSED-MASS`（(a)(b)(d) `[T]`）+ `y01` 的 `[KL]`（`[T]`） | **已补（数学层）**，未回填主证明/总册 | `I1` 可回 `[T] mod [KL]`，而 `[KL]=[T]` |
| 四壁几何 factor-4 | `P3-FOURWALL-DERIVATION`（`[T]`） | 已补 | 不再威胁结论（与 `:42`"残数结论存活"一致） |
| 全表示 KP | `P3-FULLREP-KP` + Seiler（`[T]`） | 已补 | 同上 |
| `x02` 的 `(L1)` `[C]`、`L1-c` `[O]` | `X01-L1C`（`[C]`→`[T]`） | 已补 | `I2` 侧的第一环恢复 `[T]` |
| `ν≥0`（`[KL]`） | `y01`（`[T]`，标准 OS'78） | 已补 | 唯一命名输入不构成缺口 |

**⇒ 精确结论（两档，必须分开写，不得合并）**：
1. **总册口径（当下）**：`I1` = `[条件/待重证]` = `[C]`。`I2` 链第一环不在 `[T]`。**判 (b)。**
2. **证据层（可回填）**：三条降级理由各有 `[T]` 伴文，`[KL]` = `[T]` ⇒ 回填后 `I1` = `[T] mod [KL]` = `[T]`（域：`β<1/48`、固定 `a`、柱面、`L`-一致，**不取 `L→∞`、不取 `a→0`**）。**判 (c) 且"已补"。**
3. **不许做的一步**：在 `MASTER:42` 与 `PROOF-P3-UNCONDITIONAL.md` 未回填前，把 `I1` 写回 `[T]` —— 那正是 `U07` 的越级（它未处理审计 §1.3，也未引 `P3-DRESSED-MASS`）。`[J]`

---

## §5 路径（下一个靶）

1. **P0 行政**：把 `P3-DRESSED-MASS` 的 (a)(b)(d) 与 `P3-FOURWALL-DERIVATION §2` 的 Lemma 2.1/Cor 2.2–2.3 写回 `PROOF-P3-UNCONDITIONAL.md §3–§4`（替换被审计判错的 1-per-slice 几何与"asserted" 重求和），并把 `[KL]` 显式进口。→ 之后 `P3` 可正式改标 `[T] mod [KL]`，`MASTER:42` 同轮回填。
2. **不要**再证 `β<1/48` 上的 `f₁²>0`（`U07 §6` 同判：零新地盘；该域的隙早由 OS'78/SZZ'23 独立证过）。
3. **真靶 = `I2`（隔离 `Λ★>m₁²`）在 `β<1/48`**：这是 `(CP)` 唯一另一半，且 `x02` 已判 Theorem A **不生产** `I2`（只翻译率）。障碍逐字在 `ch05:74–76`：KP/OS 余项是 `t`-一致 `O(u)`，**永远不是** `O(e^{−δt})` ⇒ 需第二指数尺度（换物种，不是加常数）。
4. 平台记法：`T⁴_L` 应回填为柱面 `𝒞_L`（`U07:73`、`X02:67` 双记账，勘误级）。

---

*本文件只写 `research/v3/species/Y06-P3DOWN.md`；未 git commit。No QED。Clay 0%。*

# Y01-KL — `[KL]` 逐字定位、载体、可证性判定

> 2026-11-25 · 任务：把 `I2` 路径的最后一项 `[KL]` 逐字定死并判能否证。
> 上游：`X01-L1C.md`（`:5`、`:107–118`）、`LOOP-RUNBOOK.md`（`:5824`、`:5834`、`:5838`、`:5842`）。
> 标签：`[T]` 定理 · `[C]` 条件定理 · `[R]` 缺步 · `[M]` 记号/勘误 · `[J]` 判词。

---

## §0 检索记录（AGENTS §1 四项）

| # | 检索式 | 工具 | 工具状态 | 命中 |
|---|---|---|---|---|
| 1 | `grep -r "\[KL\]\|KL\b" research docs --include=*.md` | grep (ripgrep) | OK（有效响应，非 provider failure） | **558**（显示截断于 250，全量落盘 `dsh-spill-…/98ff96084c12-grep-results.txt`） |
| 2 | 同上 `--include=*.tex` | grep (ripgrep) | OK | **5** |
| 3 | `grep -r "\[KL\]" research docs --include=*.md`（只取**标签式**，排除全称 Källén–Lehmann 与 bianchi 的 KL-divergence） | grep (ripgrep) | OK | **56** |
| 4 | 同上 `--include=*.tex` | grep (ripgrep) | OK | **0** |
| 5 | `glob **/*CYLINDER*` | glob | OK | 1（`research/unconditional/P3-CYLINDER-RESIDUE.md`） |

**说明**：检索 1 的 558 命中里绝大多数是**裸 `KL`**（= Källén–Lehmann 全称、以及 `research/bianchi/NEW-ARCHITECTURE.md:140,154` 的 KL-divergence，与本题无关）。`[KL]` 作为**记号**只在 `.md` 出现（`.tex` 0 命中 ⇒ 论文正文里没有这个标签）。

**带方括号标签的 56 处分布**（全部列出前 20 处，其余为同一用法的重复）：

- **定义处（唯一）**：`research/unconditional/P3-DRESSED-MASS.md:135–137`
- **转述处（唯一）**：`research/v3/species/W03-L2.md:149`
- **把 `[KL]` 指到具体载体**：`research/v3/T2-TUBE-OVERLAP.md:124`（“`P3-CYLINDER-RESIDUE` for the KL/OS spectral identification `[KL]`”）、`:471`；`X01-L1C.md:118`（“`[KL]`（KL 表示，`DRESSED-MASS §5`）”）
- **`(mod [KL])` 用法（不改陈述，仅降级标记）**：`T2-TUBE-OVERLAP.md:226,238,386`；`U3-LARGE-N-BINDING.md:196,356`；`W1B-P3-AMPLITUDE.md:194`；`W1-STRONG-ALIGN.md:17,85,248,261`；`W2-MASTER-FIELD.md:47,210,255,355,427`；`W05-TDECAY.md:101`；`X02-THEOREMA.md:112,162`；`P3-DRESSED-MASS.md:15,18,107,125,167,170`；`W03-L2.md:28,138,149,162,203`
- **`I2` 路径清单处**：`LOOP-RUNBOOK.md:5641,5643,5695,5728,5739,5780,5796,5824,5834,5838,5842,5865,5881`

---

## §1 `[KL]` 的逐字陈述

### 1.1 定义（`research/unconditional/P3-DRESSED-MASS.md:135–137`）— 逐字

> **`[KL]` (named input). Källén–Lehmann / OS transfer-matrix representation on `T³_L ×_a ℤ` (companion `P3-CYLINDER-RESIDUE`).** Reflection positivity of the Wilson measure (Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; project: PROOF-GAUSS-INVARIANT-PF) yields a self-adjoint transfer matrix `T=e^{−aH}`, `H≥0`, on the gauge-invariant physical Hilbert space, and for the Hermitian `O`,
> `⟨Ω|O(t)O(0)|Ω⟩_{c} = ∫_{[0,∞)} e^{−Et}\,dμ_O^{(L)}(E)`, `dμ_O^{(L)}≥0`, `μ_O^{(L)}=Σ_n |⟨Ω|O|n⟩|² δ_{E_n} + (continuum)`.
> This is the classical Källén (1952) / Lehmann (1954) representation, here for the lattice transfer matrix (Reed–Simon IV, §XIII).

同页 `:133` 逐字给出它的**地位**：

> Everything in §1–§4 is `[T]` **except one line**: the Källén–Lehmann representation (3.4) on the cylinder together with `H≥0` self-adjoint and `Ω` the unique vacuum. This is a standard theorem, but it is **not written in the P3 file** and is the load-bearing spectral input:

### 1.2 转述（`research/v3/species/W03-L2.md:149`）— 逐字

> | `[KL]` | 柱面 Källén–Lehmann 表示 `C_L(t)=∫e^{−Et}dμ_O^{(L)}`，`μ≥0` | `[KL]`（单一解析输入，`P3-DRESSED-MASS §3.3`） |

### 1.3 与 `X01` 的对应

`X01-L1C.md:118` 逐字：「`[KL]`（KL 表示，`DRESSED-MASS §5`）是唯一余下的命名输入。」⇒ `X01` 的 `[KL]` **就是** §1.1 这一条，不是别的东西。

---

## §2 框架假设 还是 数学缺口

**判：既不是"框架假设（不可证）"，也不是"数学缺口（未解）"，而是"标准定理、当时未在 P3 文件内写出"。** `[T]`

`x01` 说「框架性，不是数学缺口」（`LOOP-RUNBOOK:5643`、`:5838`）——**半对**：

| 方面 | 核 | 标签 |
|---|---|---|
| 不是数学缺口 | 对。它**不是开放问题**：它有标准证明（OS'78 + 有限 `L` 转移矩阵紧性 + 谱定理），无任何 YM 动力学内容。 | `[T]` |
| 「框架项」的措辞会误导 | 它**不是"假定"**：它没有留下可调的自由参数，也不是只能"假定存在"。`P3-CYLINDER-RESIDUE.md` 已把它**证出**。 | `[T]` |
| 「未完成」（`:5838` 第 2 条） | **过期判断**。载体文件 `P3-CYLINDER-RESIDUE.md`（2026-11-24，晚于 `PROOF-P3-UNCONDITIONAL`）已逐条做完 §1–§3，且全部标 `[T]`。 | `[M]` |

---

## §3 它出现在 `I2` 路径的哪一环 —— **精确定位：Theorem A 的 `⇐`，作为其"正性"前提**

路径（`X01-L1C.md:111–117`）：

```
L1-a [T] → (L1-c.1) [T] → (L1-c.2) [T] → (L1) [T] → X02 §2 恒等代数 [T] → Theorem A ⇐ → I2
```

- **不在 `(L1)`**：`(L1)`（`|log D_fin−log D_∞|≤B e^{−δt}`）是 KP/Mayer 的组合内容，**零谱论**（`X01:125` 逐字：「定理 A 是纯测度论…L1-c 是 KP/Mayer 形式下的'边界条件缺陷'引理…**定理 A 不知道 `Z_dec`，L1-c 不知道 `ν`**」）。
- **不在 `X02` 代数**：`X02` 给的是 `log D→η` 的恒等代数（`X01:114`），把组合量翻成率。
- **在 Theorem A 的 `⇐`**：`REPORT.md:110` 逐字：「注意方向：`⇐` 需要**正性**（`ν ≥ 0`，即 `ν` 真是某个谱测度）——没有正性，符号测度的 Laplace 变换可以任意快衰减，结论假。」⇒ **`[KL]` 供给的正是这个 `ν≥0`**：它断言 `C_L(t)=∫e^{−Et}dμ_O^{(L)}(E)`、`dμ_O^{(L)}≥0`。没有它，`e^{−Et}` 的 Laplace 反演只是形式运算，`Z=f₁²=a₀D_∞` 无法被称为"原子重量"。

**一句话**：`[KL]` = `(L1)+X02` 的**输出（一个数 `η`/`Z`）**被读成**谱测度原子**的那一步的许可证；缺它，`I2` 的 `Λ★>m₁²` 只是"两个指数率不同"，不是"谱上有一个孤立原子"。

---

## §4 判：能否证？—— **(a)**

### 4.1 `[KL]` **可证，且已在语料内证出**（载体 `P3-CYLINDER-RESIDUE.md`，逐字）

| 步 | 内容 | 出处逐字 | 标签 |
|---|---|---|---|
| 1 | 转移矩阵存在、自伴收缩、`H≥0`、唯一真空 | 「By Osterwalder–Schrader reflection positivity of `(S,Haar)`, there is a **self-adjoint contraction** `T_L = e^{−aH_L}` on the **physical (gauge-invariant) Hilbert space** … with `0 ≤ T_L ≤ 1`, `H_L ≥ 0`, and (Perron–PF, positive kernel) a **unique** top eigenvector `Ω_L`」 | `:46–50` `[T]` |
| 2 | 有限 `L` 谱**纯离散** ⇒ 每个谱点都是原子 | 「**at finite spatial volume the spectrum is purely discrete** — there is *no continuum*, every spectral point is an **atom**」 | `:52–57` `[T]` |
| 3 | **KL 表示本身** | 「`C_L(t) = ∫_0^∞ e^{−s t} dν_L(s) = Σ_{n≥1} f_{n,L}² e^{−E_{n,L} t}`, `ν_L := Σ f_{n,L}² δ_{E_{n,L}}`, `ν_L ≥ 0`. This is the **Källén–Lehmann (spectral) representation** on the lattice. `f_{n,L}²≥0` is the spectral theorem … `ν_L≥0` is Osterwalder–Schrader reflection positivity.」 | `:67–69` **`[T]`** |
| 4 | 原子抽取公式 | 「`f₁²_L = lim_{t→∞} e^{E_{1,L}t} C_L(t)` … `f₁²_L = ν_L({E_{1,L}})` … Both are legitimate **because the finite-`L` spectrum is discrete**」 | `:78–83` `[T]` |

⇒ **`[KL]` 的逐字陈述 = `P3-CYLINDER-RESIDUE.md:67–69` 那一行，标签 `[T]`。** 它是 `[T]`（OS'78 + Seiler LNP 159 + 紧算子谱定理），**不是**假定、**不是**开放问题。

### 4.2 判死线登记：**`(a)` 触发**

> **`I2` 的路径在 `[KL]` 这一环上闭合。** `[KL]` = 「柱面 `T³_L×ℝ` 上 `C_L(t)` 有正谱测度表示 `dμ_O^{(L)}≥0`」= `P3-CYLINDER-RESIDUE.md §1.2–1.3` = **`[T]`**。
> 因此 `X01` 的链可以逐字改写为：
> ```
> L1-a[T] → L1-c[T] → (L1)[T] → X02 代数[T] → Theorem A ⇐ [T, 由 P3-CYLINDER-RESIDUE §1.3 供 ν≥0] → I2 [T]
> ```
> **余项：无。**（`x01` 的 `[KL]` 从"命名输入"降为"已证引理"。）

### 4.3 必须并列的三条限定（不许略）

1. **道的范围**：这仍只是 `0⁺⁺` 一道的 `I2`（`LOOP-RUNBOOK:5837`、`:5843`）。`m>0` 要**每一道**的隔离，且与 `I2` 正交。`[T]`
2. **`[KL]` 的载体是**有限 `L`**：`ν_L≥0` 在有限空间体积成立，`inf_L f₁²_L` 只作**一致性轴**（`P3-CYLINDER-RESIDUE §1.4`）。若把 `[KL]` 用在 `L=∞`（"原子不是连续谱阈、 survives `L→∞`"），那一层**恰好就是 `I2` 本身**（`§3` 逐字：`[R]`）⇒ 会循环。`X01` 的链只吃有限 `L` 的 `ν≥0`，**不循环**。`[T]/[R]`
3. **值口径的旁注（不属 `[KL]`）**：`P3-CYLINDER-RESIDUE §VERDICT Gap 2` 把"领头管 = 领头原子**系数**"的识别降为 `(R1)` 打扮质量重求和 `[R]` + `(R2)` 全表示 KP `[T]`。`X01` 的 `Z=f₁²=a₀D_∞` 走的是 `P3-DRESSED-MASS §4.3`（绝对收敛栈，不走 `(R1)`）——但若有人要求把 `a₀D_∞` 与 `|⟨G|O|Ω⟩|²` **逐字等同**，那里仍挂一条 `(R1)` `[R]`，与 `[KL]` 无关。`[T]/[R]`

---

## §5 结论（供上层引用）

| 项 | 结论 | 标签 |
|---|---|---|
| `[KL]` 逐字 | 柱面 `T³_L×_aℤ` 上 `⟨Ω|O(t)O(0)|Ω⟩_c=∫_{[0,∞)}e^{−Et}dμ_O^{(L)}(E)`，`dμ_O^{(L)}≥0` | `[T]` |
| 出处 | `research/unconditional/P3-DRESSED-MASS.md:135–137`（命名）；证于 `research/unconditional/P3-CYLINDER-RESIDUE.md:67–69` | `[T]` |
| 在路径哪一环 | **Theorem A 的 `⇐`，作为 `ν≥0` 前提**；不在 `(L1)`、不在 `X02` 代数 | `[T]` |
| 框架假设 / 数学缺口 | 两者都不是：**标准定理**，`x01` 说"框架性、非数学缺口"方向对，但「未完成」已过期 | `[M]` |
| 能否证 | **能，且已证**（OS'78 + 有限-`L` 紧性 + 谱定理） | `[T]` |
| 判死线 | **(a)** ⇒ `I2` 路径真正闭合，余项无 | `[T]` |

**No QED on Clay。** `I2` 闭合 ≠ `m>0`（一道 ≠ 所有道）。Clay 0%。

# M3B-OSWIN — 第二关核验：`β₀^{OS} > 0.018394` 吗？

> 靶（`L3C-STARSTRONG.md` P1，第 154 行逐字）："查 OS'78（Ann. Phys. 110 (1978) 440）的 `β₀` 常数，写 `β₀^{OS} > 0.018394`（或证明 `0.018394 < min(β₀^{OS}, 1/48)`）"。
> 本文件只做一件事：**给数与行号，或按 `AGENTS §1` 明记"未读到"**。禁止把 `[O]` 写成"不存在"。

---

## §0 结论（一句话）

**`β₀^{OS}` 的显式数值：未读到 `[O]`。** 语料 0 命中、OS'78 正文在 Elsevier 付费墙后（只读到摘要），二次文献读到的都是符号 `β < β₀` 或 "β small enough"。⇒ **`β₀^{OS} > 0.018394` 判不出（判死线 (c)）**，**不是**"不成立"，也**不是**"成立"。
**可替代的那条路成立**：SZZ（arXiv:2204.12737，摘要逐字）给 `|β| < 1/(16(d−1))`，`d=4 ⇒ 1/48 = 0.020833 > 0.018394` ⇒ **窗被一个已发表的全观测量簇聚结果覆盖，但那不是 OS，且其量词是 "a large class of observables"** ⇒ `[C]`，不足以把 P1 升 `[T]`。

---

## §1 逐字核语料：`β₀^{OS}` 有没有数？

**工具状态**：`pwsh` + `Select-String` OK（非 `TOOL-FAIL`；`grep` 工具本会话未用，改用 `Select-String` 递归，返回有效行）。

| # | 检索式（对 `research/**`、`docs/**` 的 `*.md`、`*.tex`） | 命中 | 结果 |
|---|---|---|---|
| 1 | `Osterwalder\|OS'78\|Seiler\|β₀\|beta_0` | **40+（显示 40）** | 全部是 `β<β₀` 符号式、或一环系数 `β₀=11N/48π²`（命名撞车，非收敛半径） |
| 2 | `β₀\^\{?OS\|β_0\^\{?OS\|beta_0\^\{?OS\|0\.018394\|0\.073576` | **0**（对 `β₀^{OS}` 这个记号）；`0.018394/0.073576` 命中在 `docs/paper-v2/ch20-unconditional-i1.tex:40,169,426` | 记号 `β₀^{OS}` **语料里一次都没出现过**（除 `L3C-STARSTRONG` 自己） |
| 3 | `1/48\|0\.0208333` 且同行含 `OS\|Osterwalder\|Seiler\|半径\|radius` | **0** | 语料里 `1/48` 的归属**全部**是 SZZ / P3，**没有一处**归给 OS'78 |
| 4 | `OS'78\|OS-78\|Osterwalder.Seiler\|osterwalder-seiler` | **25（显示 25）** | 逐条见 §1.1：全是"簇聚/公理化/面积律"，**无一个数** |

### §1.1 语料里 OS'78 的每一条说法（逐字，带文件:行）

- `research/confinement/REPORT.md:167`："收敛半径 `β₀>0`" —— **只给正性，不给值**。
- `research/confinement/REPORT.md:218`："格点强耦合面积律（全 N，`β<β₀`）… Osterwalder–Seiler 1978"。
- `research/confinement/REPORT.md:153` / `STATUS.md:37,39`："`Δ(β)>0`（Osterwalder–Seiler）"、"强耦合质量隙 ✅ OS78"。
- `docs/paper-v2/ch20-unconditional-i1.tex:104`："`β<β₁` the positivity of the mass is Osterwalder--Seiler's cluster expansion"（符号 `β₁`，非 OS 的 `β₀`）。
- `docs/paper-v2/ch11b-related.tex:32`："Osterwalder--Seiler … the source of P3's rigorous polymer form"（无论文数值）。
- `docs/paper-v6/sec07-assembly.tex:91`："在 `β<β₁` 上 `m>0` 由 OS'78 的全观测量聚类"（符号）。
- `docs/yang-mills-prior-work.md:91,152`、`docs/conditional-theorem.md:19`：仅"公理化/已发表"。
- **唯一带具体半径的数**：`research/v3/species/LIT-CHECK.md:94,111` 逐字 `|β| < 1/(16(d−1))` for `SU(N)`，`d=4 ⇒ 1/48` —— **归属 SZZ 2023，不是 OS**；`LIT-CHECK.md:23,113` 还专门记录了 `1/24` 在源文**不存在**的更正。
- **自证**：`L3C-STARSTRONG.md:191` 逐字："**`β₀^{OS}` 的数值**：**未查到**（检索式 `β₀\|β_0`，工具 OK，命中全是 `β₀=11N/48π²` 的一环系数，**非** OS'78 的收敛半径）⇒ P1 是真空缺，不是"没找"。" —— 本轮**独立复现**了这条（§1 表 1–4，含新增的 `β₀^{OS}` 记号检索与 `1/48` 归属检索两条它没做的）。

**⇒ 语料裁决：`β₀^{OS}` 数值 0 命中（检索式 1–4，工具 OK）。属 `[O]`，不属 `[NG]`。**

---

## §2 查文献：OS 1978 的强耦合隙窗是多少？

### §2.1 一手源（OS'78）

`K. Osterwalder, E. Seiler, "Gauge field theories on a lattice", Ann. Phys. **110** (1978) 440–471`，**DOI `10.1016/0003-4916(78)90039-8`**（语料 `docs/paper-v2/refs.tex:23`、`docs/paper/body4.tex:283` 同）。

- **读到**：摘要/著录（Google Scholar + INSPIRE 著录页）："prove existence and analyticity of the infinite volume limit of strongly coupled Yang-Mills theories on the lattice … verify Wilson's confinement bound … rigorous treatment of the Higgs mechanism"。
- **未读到**：**正文**（Elsevier 付费墙；1978 年文，无 arXiv 版）。**故摘要层面不含任何显式 `β₀` 数值。** ⇒ `[O]`，不得记 `[NG]`（`AGENTS §1`、`§3`）。

### §2.2 二次源（能读到的都读了）

| 源 | 读到什么 | 是否给 `β₀^{OS}` 的数 |
|---|---|---|
| SZZ, arXiv:**2204.12737**（CMP 400 (2023) 805，DOI `10.1007/s00220-022-04609-1`）| **摘要逐字**：`|β| < (N−2)/(32(d−1)N)`（SO(N)）、`|β| < 1/(16(d−1))`（SU(N)）；"correlations of **a large class** of observables decay exponentially … strictly positive mass gap" | **不给 OS 的数**（这是 SZZ 自己的窗） |
| Faria da Veiga–O'Carroll, arXiv:**2509.03513**（v2，2026-01）| **摘要 + 正文前段**（HTML）：只写 `0 < β ≡ 1/g² ≪ 1`、"for small β" | **不给** |
| Suzuki 强耦合综述 slides（Yukawa QFT2008）| "strong coupling expansion has a finite radius of convergence (Osterwalder–Seiler, Kotecký–Preiss)" | **只给"有限半径存在"**，无数 |
| Seiler, LNP **159** (1982)（专著，语料 `body4.tex:286`）| 只读到著录/被引 | **未读正文，不给数** |
| 若干 2025–2026 preprint / zenodo / vixra（如 zenodo 23354163 写 `Δ_OS=(1/a)ln(1/(12β))`）| 出现 `1/12` 型数字 | **不采信**：非同行评审、无定理号可核、与 OS'78 原文无法对齐 ⇒ 记 `[U]`（unverified source），**不得进账** |

**检索式（`web_search`，工具 OK，共 6 组 13 条 query）**：
`"Osterwalder Seiler 1978 Gauge field theories on a lattice strong coupling convergence radius beta < 1/"`、
`"Osterwalder-Seiler strong coupling threshold beta < 1/ SU(N) 4d polymer expansion"`、
`"Seiler 1982 LNP 159 strong coupling radius beta_0 value"`、
`"\"Osterwalder\" \"Seiler\" \"β < 1/12\""`、
`"lattice gauge \"strong coupling\" \"|β| <\" \"(2d−1)\" Osterwalder Seiler polymer"`、
`"\"1/(12β)\" mass gap lattice strong coupling"`。
**命中：每条返回 10 源；其中含 `β₀^{OS}` 显式数值的：0。** 失败标记：无 `TOOL-FAIL`（`web_fetch` 一次 404 是 `2204.12737v3` 版本号不存在，改 v1/abs 后 200，已排除）。

---

## §3 判：`β₀^{OS} > 0.018394` 吗？

**判不出。触发判死线 (c)。** 逐字：

- `β₀^{OS}` = **未读到**（`[O]`；检索式 §1 表 1–4 + §2.2 六组，工具 OK，命中 0 个显式数值）。
- ⇒ 命题 `β₀^{OS} > 0.018394` **本轮既未证成、也未证伪**；`(a)` 不能引用，`(b)` 的区间也无法点名（区间 `(β₀^{OS}, 0.018394]` 的**左端点未知**）。**不许**把"未读到"写成"不存在"，也**不许**把它写成"成立"。
- **P1 状态维持 `[O]`**，与 `L3C-STARSTRONG.md:154` 的记账一致（该行本来就写 `[O]→[T]`，P1 未做）。

### §3.1 顺带钉死的两条算术（是 `[T]`，且不依赖 OS）

1. `1/48 = 0.0208333…`；`0.018394 = 0.88291×(1/48)`（`docs/paper-v2/ch20-unconditional-i1.tex:169,426` 逐字）。⇒ `0.018394 < 1/48`，比值 `1.1327`，绝对差 `2.439e-3`。
2. **SZZ 的窗覆盖我们的窗**：`|β| < 1/(16(d−1))`，`d=4` ⇒ `1/48 = 0.020833 > 0.018394` ✔（摘要逐字，上表）。⇒ **若把 P1 的"OS"换成"SZZ"，窗对齐这一步在数值上成立**。

### §3.2 但 SZZ 替换**不能**直接把 P1 升 `[T]`（这条必须留痕）

P1 要 OS 的东西是**逐道、`β`-一致**的 `m>0`（`L3C-STARSTRONG.md:88–89` 逐字："`m>0` 在 `β<β₀` 由 OS'78 逐道给出"）。而 SZZ 摘要逐字是 "correlations of **a large class of** observables decay exponentially" —— **量词不是"所有道"**，且其质量隙结论来自 LSI/Poincaré + 指数簇聚，**给出的是"存在 δ>0"，不是"每道一个可比 δ"**。⇒ SZZ 替换 ⇒ `[C]`（窗覆盖 `[T]`，逐道一致 δ 仍 `[O]`），**不封 P1**。

---

## §4 记账 / 下一个靶

| 项 | 状态 |
|---|---|
| `β₀^{OS}` 数值 | **`[O]`（未读到；语料 0 命中 + OS'78 正文付费墙；检索式与命中数见 §1–§2）** |
| `β₀^{OS} > 0.018394`？ | **不可判（判死线 (c)）** |
| `0.018394 < 1/48` | `[T]`（§3.1） |
| 窗被 SZZ 的 `1/48` 覆盖 | `[T]`（摘要逐字） |
| 逐道一致 `δ`（窗内所有道） | `[O]`（SZZ 量词是 "a large class"；OS'78 未读到正文） |
| **P1** | **仍 `[O]`**；`(★-gap)` 在 `β<0.018394` 上**未封口** |
| 下一个靶（≤1 步） | 取 OS'78 正文（馆际/馆藏扫描或 Seiler LNP 159 第 3 章）读出定理号后的 `β₀(d,N)` 显式式；**若读不到 ⇒ 换 P2+P3（道盲余项 + `C ≤ 2.691`）**，不在这条上耗轮 |

**限定**：`β` 为 Wilson 逆耦合（`β=g⁻²`），`d=4`，`SU(N)`，固定 `a`、有限 `L`；`β≥1/48` 一律不谈；`β≈2.5` 侧无关。
**No QED；Clay 0%。** 本文件新增的只有：`β₀^{OS}` 的**空账被独立复核**（含两条语料未做的检索式）+ `0.018394 < 1/48` 的可替代路径与它的量词上限。

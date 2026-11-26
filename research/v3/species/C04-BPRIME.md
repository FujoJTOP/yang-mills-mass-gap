# C04 — `(b′)` 是什么：§14.2 第五项的定名、出处与结算

> 本文件只回答一件事：`LOOP-RUNBOOK §14.2` 限定⑤ `(b')` 的**准确含义**、**出处**、**是否已解决**。
> 不改任何源文件，不 commit。

---

## 0 · 判（一句话）

**`(b′)` = IN-11 的弱形式（成员性前提）：`δ` 之上不存在落在 `{0⁺⁺,2⁺⁺}` 之外的 bulk 态**（`0⁻⁺`、奇 C、自旋 ≥3 均 `≥δ`）。
它是 **IN-b 的 `(b)`**（"最轻 bulk 态在 `tr F²` 道，即 `0⁺⁺`"）的**弱化**，**不是** `A02 §4(b)` 那条判死线分支的变体。
**在本文件任务的窗口（`β<β₁`）内它已解决 `[T]`**（`PROOF-IN11.md:19`）；窗口外（`β≈2.5`）仍是 `[C]` 物理输入。**它不是第 6 个未关项。**
**但 §14.2 的"从未定义"是真的**：`LOOP-RUNBOOK.md:7018` 写 `⑤ (b') —— 见下`，**§14.2 之后并无定义**，前向引用悬空。

---

## 1 · 检索执行核验（`AGENTS.md §1` 四项）

**工具状态**：`grep` = **OK**（5 次调用，全部返回带行号正文，**非** `ripgrep provider failure`）；`read` = OK（5 次：`IN-B.md` 全 33 行、`IN-11-attack.md` 全 37 行、`PROOF-IN11.md` 1–45、`A02-2PPI2.md` 30–109、`LOOP-RUNBOOK.md` 6990–7034、`B-CHECK.md` 1–36）；`glob` = OK（定位 `A02-2PPI2.md`）；`pwsh` = OK（列目录）。**无 `TOOL-FAIL`。**
**检索式与命中数**（路径 `D:\Dev\math`）：

| # | 检索式 | include | 命中 |
|---|---|---|---|
| 1 | `b′` | `*.md` | **117** |
| 2 | `\(b′\)`（带括号的标号，去噪 `b′(0)=G″` 曲线族） | `*.md` | **23** |
| 3 | `\(b'\)`（ASCII 撇） | `*.md` | **5** |
| 4 | `b-prime` | 全部 | **0** |
| 5 | `b′` | `*.tex` | **0** |
| 6 | `b'` | `*.tex` | **0** |
| 7 | `见下`（查前向引用是否落地） | `*.md` | **51**（`LOOP-RUNBOOK` 内 2：`:99`、`:7018`） |

⇒ `(b′)` **不是误记**（判死线 (b) 不成立）：语料里有 23 处带括号的显式标号，其中 9 处集中定义同一个东西。

---

## 2 · 出处（逐字）

| 位置 | 逐字 | 标签 |
|---|---|---|
| `research/defect-bootstrap/IN-B.md:12` | **「**(b′) the lightest bulk state lies in `0⁺⁺` **or** `2⁺⁺` — i.e. **no bulk state below `δ` carries other quantum numbers**」 | 定义句（原始出处） |
| `IN-B.md:17` | 「This is a conditional theorem with the **weaker** premise `(b′)`. **Caveat:** `(b′)` must exclude sub-`δ` states in *every other* channel (`0⁻⁺`, odd C, higher spin), since the bootstrap bounds only the two C-even, P-even channels it sees.」 | 承重范围 |
| `docs/conditional-theorem.md:207` | 「等价弱形式 `(b′)`：在 `δ` 之上**不存在**带其他量子数的 bulk 态（`0⁻⁺`、奇 C、更高自旋均不低于 `δ`）」 | 等价重述 |
| `research/unconditional/IN-11-attack.md:4` | 「Weak form `(b′)`: no bulk state below the certificate `δ` lies outside the C-even, P-even scalar/tensor sectors」 | 同 |
| `research/unconditional/PROOF-IN11.md:19` | 「This **proves (b′)** at strong coupling: no `0⁻⁺`, odd-C, or spin-≥3 state is the lightest.」 | **已解决的判据** |
| `research/unconditional/C2-bypass-INb.md:17,23` | 「Only `K=0` transfers the bound, and `K=0` *is* `(b′)` = IN-11」；「The weakest monotone repair is `Δ_{c₀} ≥ cert`, i.e. IN-11 weak form `(b′)`」 | 不可再弱 |
| `research/defect-bootstrap/CONDITIONAL-SUMMARY.md:45` | 「弱形式 `(b′)` 要求 `δ` 之上无其他量子数的 bulk 态……**物理输入，不可由 bootstrap 导出**」 | `[C]` 侧 |
| `research/v3/species/W95C-ORDER.md:89` | 「**成员性 `(b′)`** … 强耦合 `[T]`（`PROOF-IN11.md:19`）；`β≈2.5` 侧 `[C]` 物理输入」 | 两栏 |
| `research/v3/species/A02-2PPI2.md:102` | 「4. **`(b′)` 成员性仍承重**：`0⁻⁺/1⁻⁻` 不需要各自证书，但需要"`δ` 之上无 bulk 态"」 | **`a02` 自己就定义了它** |

⇒ `A02` **不是**"引用了未定义的 `(b′)`"：它在 `§4(a)` 的四条限定第 4 条里写清了。`§14.2`（`LOOP-RUNBOOK:7018`）的 `⑤ (b') —— 见下` 才是那个"从未定义"的洞——**"见下"无下文**（§14.2 在 `:7024` 结束，`:7020–7024` 只讲"剥掉 `0⁺⁺` 严格最轻"，无 `(b′)` 定义）。

---

## 3 · `(b)` 是什么，`(b′)` 是不是它的变体

**IN-b 的 `(b)`**（`B-CHECK.md:1,19` 逐字）：
> 「(b) **the lightest bulk state couples to `tr F²`**」/「(b) needs that the **lightest** bulk state is spin-0. The trace only ties `tr F²` to the spin-0 *sector*, **not** to the *lightest* state.」

**`(b′)` 是它的弱化，不是同名重写**（`B-CHECK.md:24–29` 逐字）：
> 「## Weakening (not promotion) … one could bound the spin-0 **and** spin-2 channels and take the min. Then (b) weakens to "the lightest state couples to the bulk-defect OPE (`D·D`)"」

⇒ **结论**：`(b′)` **是 `(b)` 的变体（弱化）**，但变体对象是 **IN-b 的 `(b)`**（`0⁺⁺` 单独最轻），**不是** `A02 §4(b)`「不触发：阈值不同 ≠ 链要改」那条分支。两者只是同名的不同枚举：
- `A02 §4` 的 `(a)/(b)/(c)` = **判死线三分支**（触发/不触发/方向相反）；
- `IN-B` 的 `(a0)/(a2)/(b′)` = **前提清单**；`(b′)` 取自这里。

`LOOP-RUNBOOK:7015–7018` 的限定①②③④⑤里，前面四项都是 a02 的行内标签，唯独⑤照抄了 IN-B 的标号却没带定义 ⇒ 混淆的**生成点**在这里。

---

## 4 · 最可能的读法之辨：是不是"`O(1)` 余项符号未定"

**不是。两者正交，且互补。**

| | `O(1)` 符号未定 | `(b′)` 成员性 |
|---|---|---|
| 讲什么 | `0⁺⁺` 与 `2⁺⁺` **谁更轻 + 是否被隔离**（对内的排序） | **`{0⁺⁺,2⁺⁺}` 之外**（`0⁻⁺`/奇 C/自旋≥3）在 `δ` 之下有没有态（对外的排除） |
| 出处 | `PROOF-IN11.md:3,23`「`a_{0⁺⁺}=a_{2⁺⁺}=1` … the ordering and isolation of `0⁺⁺` live entirely in the `O(1)` subleading splitting」；`A02:84`「这落在 `O(1)` 次领头，正是 `W95C §2` 说的"符号未定"的那一层」`[M]/[C]` | `IN-B.md:12` 起，见 §2 表 |
| 若它不成立 | 无碍：`min(Δ₀,Δ₂)` 照样 `≥δ` ⇒ 结论**不动** | **致命**：出现更轻的他道态 ⇒ `inf spec` 不在 `{Δ₀,Δ₂}` ⇒ 链断 |
| 状态 | 未做（"not carried out to the needed uniformity"，`PROOF-IN11:23`） | 强耦合 `[T]`；`β≈2.5` `[C]` |

**正是 `(b′)` 把"`O(1)` 符号未定"变成无害**：`A02:96` 逐字「**`min(Δ₀,Δ₂)` 的两个因子都有证书，`m>0` 在强耦合窗不再依赖"`0⁺⁺` 严格最轻"**」——因为只要两道都 `≥δ` 且无他道更轻，谁在 `O(1)` 上赢**不进结论**。
⇒ 把 `⑤ (b')` 读成"`O(1)` 余项符号未定"会把一个**已证 `[T]`** 的项误记成一个**未做**的项，方向相反（前者是"减承重"，后者是"加承重"）。

---

## 5 · 结算：它解决了吗

**(a) 有明确出处与含义 ⇒ 采信。已解决的部分与未解决的部分必须分开记：**

| 窗口 | `(b′)` 状态 | 依据 |
|---|---|---|
| `0<β<β₁`（**= `A02`/§14.2 自己的窗口**） | **已解决 `[T]`** | `PROOF-IN11.md:15–19`：OS/SZZ 收敛 + `O_h` 立方群引理 ⇒ `a_c≥2 (c∉{0⁺⁺,2⁺⁺})` ⇒ `m_c − m_{0⁺⁺} ≥ (a_c−1)ln(1/u) → +∞` ⇒ `inf_c Δ_c = min(Δ_{0⁺⁺},Δ_{2⁺⁺})` |
| `β≈2.5`（交叉点侧） | **未解决 `[C]`**，且**不可导出** | `PROOF-IN11.md:27` 解析性墙；`IN-11-attack.md:32`「**IN-11 is definitively physical** — neither derivable nor conditionally derivable」；翻转条件 `IN-11-attack:35` 也不存在 |
| 带物质（轻夸克） | **假** | `IN-11-attack.md:20–21`：`π`（`0⁻⁺`）最轻；"是动力学，非对称性保护" |

⇒ 对 `A02` 的强耦合命题：`(b′)` **已关**，**不构成第 6 项**。`A02:100` 限定②「不跨 `β≈2.5`」已经把它挡在窗外——两项不自相矛盾，⑤ 是**记账用**的诚实标签，不是新债。

**两处必须一并记的文书缺陷（都不是新数学洞）：**

1. **前向引用悬空**（成因）：`LOOP-RUNBOOK.md:7018` `⑤ (b') —— 见下`，§14.2 后无定义。RUNBOOK 全文 ASCII `(b'` 仅此 1 处，Unicode `b′` 0 处 ⇒ **RUNBOOK 从未定义过 `(b′)`**。修法：`:7018` 改写为 `⑤ (b′)_{IN-11}：δ 之上无 {0⁺⁺,2⁺⁺} 之外的 bulk 态（IN-B.md:12；强耦合已证 PROOF-IN11.md:19）`。
2. **窗口记账差一行**：`PROOF-IN11` 的 `β₁`（OS/SZZ 收敛半径）与 `A02:96` 的 `β₁ = min(1/48, c/N)` 必须再取 `min` 才保证**同一窗**。`A02` 未显式并入 ⇒ 建议写成 `β₁ = min(1/48, c/N, β₁^{IN-11})`。一行记账，`[T]` 不改判。

**（判死线归属：`(a)`）** —— 不是 (b)（命中 ≠ 0，`b-prime` 的 0 命中不构成否证，定义用的是 Unicode `′`）；也不是 (c)（它不指向一个未关的实质限定，在 `A02` 窗口内已 `[T]`）。

---

## 6 · 附：`(b′)` 在仓里是**多义标号**（引它必须带限定）

| 位置 | 该处的 `(b′)` 指 |
|---|---|
| `IN-B.md:12`（**本义**，`A02`/§14.2 引用的是这个） | IN-11 弱形式：`δ` 之上无 `{0⁺⁺,2⁺⁺}` 之外的 bulk 态 |
| `research/unconditional/OS-P-INTERFACE.md:105,192` | `coarse⊕fiber` Feshbach–Schur 归约在 `H_phys` 上成立（`[真定理]`，靠正交性） |
| `research/unconditional/TOOL-FREE-PROBABILITY.md:117,123` | `ρ` 的 hard support gap ⇒ `spec(H)|_{0⁺⁺} ⊆ {0}∪[Δ,∞)` |
| `research/unconditional/BS-KERNEL-SIGN.md:15` | 某 `(b)` 的攻击分支（符号/阶梯控制） |
| `research/v3/species/R2-H2GLOBAL.md:79` | §1.3 `U6` 的 `B_perp`（勘误 E1：该对象在语料里不存在） |

⇒ 建议全仓把 §14.2 那一项写成 **`(b′)_{IN-11}`**，与 `(b′)_{Schur}`、`(b′)_{supp}` 区分。这是本次"从未定义"事故的结构性原因：**同一个符号在至少 5 个文件里各指一件事**。

---

*Refs.* 内部：`research/defect-bootstrap/IN-B.md:3,10–19,23–26`；`research/defect-bootstrap/B-CHECK.md:1,17–29,34`；`research/defect-bootstrap/GAP-DICTIONARY.md:20,25,32`；`research/defect-bootstrap/CONDITIONAL-SUMMARY.md:45`；`docs/conditional-theorem.md:207–208,389`；`research/unconditional/IN-11-attack.md:4,20–21,32,35`；`research/unconditional/PROOF-IN11.md:3,11–19,23,27`；`research/unconditional/C2-bypass-INb.md:17,23`；`research/v3/species/W95C-ORDER.md:85,89,94,121`；`research/v3/species/A02-2PPI2.md:48,84,96,98–102,131,161`；`research/v3/species/LOOP-RUNBOOK.md:7000–7024`；`research/v3/species/R2-H2GLOBAL.md:79`；`research/unconditional/OS-P-INTERFACE.md:105,192`；`research/unconditional/TOOL-FREE-PROBABILITY.md:117–124`；`research/unconditional/BS-KERNEL-SIGN.md:15`。
外部（转引，本轮未新核）：Morningstar–Peardon, `hep-lat/9901004` / PRD **60** (1999) 034509（`0⁺⁺` 1710 < `2⁺⁺` 2390 < `0⁻⁺` 2560 MeV）；Osterwalder–Seiler, Ann. Phys. **110** (1978) 440；Shen–Zhu–Zhu, CMP **400** (2023) 805 / arXiv:2204.12737。
No QED；Clay 0%。

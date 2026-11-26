# U02-CP — `(CP)` 是什么：一个符号、**两个**命题；义 A 可证伪（判 (a)），义 B 是速率一致化；**`(CP) ⇒ H1`，故 `t02` 的"两个目标"在合取上冗余**

> 2026-11-25 · 任务：把 `(CP)` 写成可证的命题。判词：**(a) 义 A 有明确、可证伪的陈述**（`I1∧I2`，逐字有出处）；但语料同一符号 `(CP)` 承载**两个不同命题**（义 A 谱/隔离、义 B 速率一致化 `d5/(CP)=p=1`），从未并列辨析 ⇒ `t02` 的"两个目标"一半是**记号冲突**。义 A 严格强于 `f₁²>0`（多出 `I2` 与窗口一致化），且 `[CP-b]⇒[CP-a]`、`I2` 一致化 ⟺ `m₁≥μ` = gap ⇒ **`(CP) ⟹ H1`；`H1∧(CP) ≡ (CP)`**，即 `t02` 的两目标是**包含关系，不是两条独立路**（判 (c) 的实质成立，见 §4）。

---

## §1 检索记录（`AGENTS.md §1` 四要素）

| # | 工具状态 | 检索式 | 命中数 |
|---|---|---|---|
| 1 | `grep` OK（ripgrep 正常返回，非 ENOENT/非 provider failure） | `\(CP\)` in `research/**/*.md` | **165** |
| 2 | 同上 OK | `\(CP\)` in `docs/**/*.md` | **15** |
| 3 | 同上 OK | `\(CP\)` in `**/*.tex`（全仓） | **8** |
| 4 | 同上 OK | `\(CP\)` in `research/unconditional/INVENT-TRANSPORT-INEQ.md` | **2** |
| 5 | 同上 OK | `\(CP\).{0,80}(L=1|N=∞)|(L=1|N=∞).{0,80}\(CP\)` in 全仓 | **0** |

**假阳性剔除**：165 命中中至少 1 条不是本对象 —— `research/defect-bootstrap/ANOMALY.md:27` 的 `(CP)` 是 **charge-parity**（`θ=π` CP），与谱/隔离无关。⇒ 真命中 ≈164。**无 TOOL-FAIL。**

---

## §2 逐字出处（两义）

### 义 A — 谱/隔离义（搬运定理的条件）

**A1** `docs/YANG-MILLS-MASTER.md:119`：

> `2. **(CP)**：单 T=0 禁闭相（无体相变）+ 0⁺⁺ 边 Kato-保持 + 一致隔离。`

**A2** `docs/paper-v2/ch13-machine.tex:40-41`：

> `\item \textbf{(CP):} a single $T=0$ confining phase (no bulk transition) $+$ a $0^{++}$ edge that`
> `  is Kato-preserved $+$ uniform isolation.`

**A3** `research/unconditional/INVENT-TRANSPORT-INEQ.md:17`（原始出处，带量词的版本）：

> `Assume **(CP)**: β_eff<1/48 and β=2.5 lie on one constant-physics line — a single confining phase, lightest 0⁺⁺ Kato-isolated (Λ(β)−m₁²≥κ>0) and carrying the edge throughout (no bulk transition).`

**A4** `research/unconditional/CONTINUUM-CP.md:1` / `:11`（拆分与方向）：

> `(CP) = [CP-a] ∧ [CP-b] = I1 ∧ I2`，`[CP-b] ⇒ [CP-a]`
> `[CP-b]=I1 (residue f₁²≥r₀) ∧ I2 (isolation Λ−m₁²≥κ)；... I2 is the upper gap, strictly above m_g. Only I1 closes, at β<1/48 (P3)`

**A5** `research/unconditional/G2-RESIDUE-ATOM.md:1`（方向判词）：

> `isolation I2 is the gap one rung up, and (CP) is downstream of H1`

### 义 B — 速率一致化义（与 A 不同轴）

**B1** `docs/paper-v4/sec04-five.tex:89`：`$=$ the corpus object $d5$\,/\,$(CP)$`；`:244` 与 `sec07-lit.tex:93`、`sec01-position.tex:89`：`rate uniformisation = d5/(CP)=p=1`。

**B2** `research/v3/species/RATE-UNIFORM.md:129`（把 B 写成量词式）：

> `d₅ / (CP)  ≡  ∃r>0 : inf_{β∈W} c_route(β) ≥ r ，  W ∋ 2.5`

**B3** `research/v3/species/RATE-UNIFORM.md:182`（明文承认两义不同）：

> `(CP) = [CP-a] ∧ [CP-b] 谈的是 ν 的支撑/隔离，与两者的轴都不同`

⇒ **[T]** `(CP)` 在本仓是**一个符号两个命题**，且 B 的作者自己声明与 A 不同轴。`t02` 的"两个目标"引的是 A；`paper-v4/H1B-FIVE` 引的是 B。

---

## §3 量词形式：它是关于**谱底部在整个 β 窗口上的一致结构**，不是关于极限、也不是关于残数单点

规范化陈述（把 A1–A3 补全量词；**[C]**，因语料只给散文）：

固定 `N`、`T=0`、开边界无关。记 `ν_β` 为 `β` 处 `0⁺⁺` 通道 KL 谱测度（`V→∞` 后），`m₁(β)=inf supp ν_β`，`Λ(β)` = 同通道次一能级，`f₁²(β)=ν_β({m₁(β)²})`。

> **(CP)** ≡ `∃κ>0 ∃r₀>0 : ∀β∈[β_eff,2.5] : (i) f₁²(β) ≥ r₀` ∧ `(ii) Λ(β) − m₁(β)² ≥ κ`
> 且 `(iii)` 同窗口 `L`-一致、`a`-一致（常数与 `V,L,a` 无关）。

- **对象**：`(i)` 是**谱底部成原子**（原子性/残数），`(ii)` 是**该原子与次一级的间隔**（隔离）。**不是** `a→0` 极限陈述，**不是**某个残数的单点值陈述，**不是**关于 `f=−(1/V)log Z` 的解析性（那是 `[CP-a]`，是**推论**）。
- **量词型**：`∃κ,r₀ ∀β ∀(L,a)` —— **"全体 β × 一致"**。既不是 `∃β` 单点，也不是 `∀n` 序列。
- **强度顺序** `[T]`：`(CP) ⇔ [CP-b] ⇔ I1∧I2`；`[CP-b] ⇒ [CP-a]`（`CONTINUUM-CP §2`）；`I1⟺I2` given moments（`GAP-G2-edge`）。
- 用 `Λ≈(1.5m₁)²` 的 grading：`Λ−m₁²≈1.25m₁²`，故 **`(ii)` 一致 ⟺ `m₁≥μ>0` 一致** ⇒ **`(CP) ⟹ H1`**（`OTHER-F1SQ-TRANSPORT:14`、`PROOF-CP-COMBINE:13`）。**[T]**

---

## §4 与 `f₁²>0` 的关系：**严格更强**，不是等价；且因此 `H1∧(CP)` 冗余

| 关系 | 判 | 依据 |
|---|---|---|
| `(CP) ⟹ f₁²>0`（在窗口上） | **[T]** | `(CP)` 的 `(i)` 就是它；`[CP-a]⇏residue` 不构成反方向 |
| `f₁²>0 ⟹ (CP)` | **✘** | 单相/有原子可**无隔离**：单相可无隙、可连续溶解（`CONTINUUM-CP:13`、`PROOF-NO-MELT-PHASE`）；`[CP-a]⇏[CP-b]` **`[T]`** |
| `(CP) ⟺ f₁²>0` | **✘**（差一个 `I2` + 一个 `∀β`） | `T02-F1SQ.md:130` 逐字：`(CP)`（`=I1∧I2`，**严格更强**）`[T]` |
| `(CP) ⟺ H1` | 方向 **`(CP)⟹H1`**；`H1⟹(CP)` **未证、很可能假** | `G2-RESIDUE-ATOM`（`(CP)` downstream of H1）；`IMP-ARROWS §6.3`：把 `(CP)` 当前件即"隔离⟹隙"，**上一级蕴含下一级** |

**⇒ 对 `t02` 的冲突点（判 (c) 的实质）**。`t02` 逐字（`LOOP-RUNBOOK.md:4590`）：

> `**目标不是 1 个，是 **2 个**：H1 与 (CP)（t02：H1 ⇏ f₁²>0，反例 β=0 给 G″=0.2688>0 而 f₁²=0）。**`

这句要分两半判：
- "`H1 ⇏ f₁²>0`" —— **[T]**，反例成立（`β=0`）。**但那是 `H1 ⇏ I1`，不是 `H1 ⇏ (CP)`。**
- "两个目标" —— **[M]**：既然 `(CP) ⟹ H1`（§3），则 `H1 ∧ (CP) ≡ (CP)`。**两个合取项里一个是另一个的推论 ⇒ 不是两条独立路，是一条更强的 + 它的推论。** 正确写法是 **`(CP)` 单一目标**（比 H1 强一级），而不是 `H1 + (CP)` 两项并列。

---

## §5 已知成立的窗口（判 4）

| 条件 | `(CP)` 状态 |
|---|---|
| `β<1/48` | **只有 `I1`**：`f₁²≥r₀(β)>0`，`V,L`-一致，**无条件**（P3，`PROOF-P3-UNCONDITIONAL`）**`[T]`**；`[CP-a]` 亦成立（单 gapped 相，OS78/SZZ23 `[T]`）。**`I2` 未证** ⇒ **`(CP)` 在 `β<1/48` 也不成立**（不是"窗口内已证"）。 |
| `β≤12/17≈0.71` | 面积律 `[T]`（Borgs–Seiler），**≠ 单相**（`CONTINUUM-CP:5`）。 |
| `β≈2.5` | **未证**；`f₁²` 与 `I2` 都开。 |
| `L=1` | **无记录**（检索式 5，命中 **0**，工具 OK）⇒ 语料从未把 `(CP)` 与 `L=1` 挂钩。 |
| `N=∞` | **无记录**（同上，命中 **0**）；large-`N` 只到强耦合，`N→∞→3` 无定理（`PROOF-LARGE-N-RIGOROUS`）。 |
| 小圆（`r→0`） | **反例方向**：`(CP)` 读作 adiabatic continuity 时对纯 YM **假**（GPY 相变），`QCD(adj)` 有四相（`SMALLCIRCLE-ADIA §5`）。**`[T]` 这是 `(CP)` 最接近的可证伪邻域。** |

---

## §6 可证伪性 + 今天能攻的形式（判 (a) 的交付）

**可证伪**（三条，任一命中即 `(CP)` 假）：
1. **K2 触零**：网格上 `f₁²(β_t)` 命中 `0` 后反弹（`MS-ANALYTIC-NOZERO` 的 K2）⇒ `(i)` 假。
2. **隔离塌**：测出 `Λ(β)−m₁(β)² ≤ 0`（次一级并入最轻 `0⁺⁺`）⇒ `(ii)` 假。
3. **相变**：`[1/48,2.5]` 上出现体相变/Fisher 零点压到实轴 ⇒ `[CP-a]` 假 ⇒ `(CP)` 假（`CONTINUUM-CP §2/§3`）。

**今天能攻的形式**（把 `∀β` 收成单点，成本最低）：
> **(CP)@2.5** ≡ `f₁²(2.5) ≥ r₀>0` ∧ `Λ(2.5) − m₁(2.5)² ≥ κ>0`，单点、无需跨 `β` 一致。
> 依据：`COV-T.md:233` 逐字——"**只需要 `β=2.5` 处的单点**、不需要跨 `β` 的一致性。这是比 `LOOP-RUNBOOK §8.2` 的 `L`-一致性死点窄得多的靶。"

⇒ 三个可执行的攻击面：(α) `Λ(2.5)−m₁²` 的 GEVP 直接测（`PROOF-EDGE-E2` 说现有工具全是**底部泛函**、对 `λ₁−λ₂` 盲，故这是新仪表需求，不是重跑）；(β) 用 Fredholm 行列式看 `z₁` 的**简单性**（Kato 域 = `I2`，`SYNTHESIS-NEWCONSTRUCTS`）；(γ) 把 `(i)+(ii)` 合成 `ν` 的**总≤原子** Schur-凸强制（`SIGNED-FOURIER-CORE §99` 已判其 Schur-凹，方向反 ⇒ 此路 `[NG]`）。

---

## §7 记账

| 项 | 判 |
|---|---|
| `(CP)` 有精确可证伪陈述（义 A） | **`[T]` 有** —— A1–A4 逐字，量词补全见 §3 ⇒ **判 (a)** |
| `(CP)` 只是口号 | **否**（但"两个目标"那句是） |
| `(CP)` 一符号两义，语料未辨析 | **`[T]`**（B3 自认） ⇒ `t02` 的记账需分义 |
| `(CP) ⟺ f₁²>0` | **`[T]` 否**，严格更强 |
| `H1 ∧ (CP)` 是否两项 | **`[M]` 否**：`(CP)⟹H1` ⇒ 合取冗余，单一目标是 `(CP)` |
| `(CP)` 在 `β<1/48` 是否已证 | **否**：只有 `I1`，`I2` 全开 |
| `L=1` / `N=∞` 是否给 `(CP)` | **无记录**（命中 0，工具 OK） |
| Clay | **0%** |

*Refs（内部）*：`INVENT-TRANSPORT-INEQ:17`、`CONTINUUM-CP:1/5/11/13`、`G2-RESIDUE-ATOM:1/7/14`、`OTHER-F1SQ-TRANSPORT:14/16`、`PROOF-CP-COMBINE:13`、`MS-ANALYTIC-NOZERO:27`、`RATE-UNIFORM:129/182`、`T02-F1SQ:24/130`、`IMP-ARROWS §6.3`、`LOOP-RUNBOOK:4590`、`COV-T:233`、`SMALLCIRCLE-ADIA §5`、`SIGNED-FOURIER-CORE:99`、`docs/YANG-MILLS-MASTER.md:119`、`docs/paper-v2/ch13-machine.tex:40-41`、`docs/paper-v4/sec04-five.tex:89`。*外部*：Osterwalder–Seiler Ann. Phys. **110** (1978) 440；Shen–Zhu–Zhu CMP **400** (2023) 805（arXiv:2204.12737）；Borgs–Seiler CMP **91** (1983) 329；Cossu–D'Elia JHEP 2009（arXiv:0904.1353）——均为**摘要级核**，未逐字核正文。

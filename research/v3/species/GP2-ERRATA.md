# GP2-ERRATA — 全项目勘误与引用错位台账（只记录，不改既有文件）

> 第 GP2 轮 · 子代理 `gp2-errata` · 工作目录 `D:\Dev\math`。
> **纪律**：只新建本文件；**未改任何既有文件**；未 `git commit`；**无 MC**；`.lean` 只读未编译。数值只用 `py -3.12` 做一次确定性算术（`exit 0`）。
> **等级**（`CLAUDE.md §2`）：`[T]` 真定理/恒等式 · `[C]` 条件 · `[NG]` 具名 no-go/缺步 · `[M]` 口径/记账 · `[E]` 可复跑数值 · `[O]` 开放 · `[J]` 判断。
> **进展单位 = 轮**（`AGENTS §5`）。**No QED。Clay = 0%。**

---

## §0 判词

```
欠账总数：18 条 = 13 条（任务给定，逐条核实）+ 5 条（查漏新增）
已核：13 / 13（每条给 文件:行 + 原文摘录 + 修正文本 + 等级）
其中：必须是错的（改正数/改符号/换证书）   = 8 条  E1 E2 E3 E5 E6 E8 E9 E12
      等级/口径需重标（结论大致可留）       = 3 条  E4 E7 E10(未落盘) E13 → 4 条
      纯记账（位置图/注册表）               = 2 条  E11 E13
查漏新增 5 条（全部 `[M]`/`[E]` 级，无一条推翻任何既有 `[T]`）：
      L1 `U3` 双指（UNIF3.md vs U3-LARGE-N-BINDING.md）——11 处引用按行号只能指后者
      L2 GA3 §1.6「悬空引注」是**路径查错**（查 species/ 而文件在 v3/）⇒ 该勘误本身要撤
      L3 `Q04-BAKRY-EMERY.md`（B01-RP:134）真悬空（Q4 vs Q04 之误）
      L5 `\bar u` 槽跨论文不一致（paper-v6 63.8=错值；paper-final sec05 用 16.99=\bar C=1 槽）
      L6 `662×` 已判口径错但仍以现读数残留在 B14/GAP2 四处（已记录、未回填）
必动论文（paper-v6 / paper-final）：E9 · E12 · E10  → 3 条
只动工作文件（research/）：E1 E2 E3 E4 E5 E6 E7 E8 + L1~L3
只动位置图：E11  只动注册表：E13
Clay = 0%；本文件不产生任何地板，也不把任何 `[C]` 升 `[T]`。
```

**一句话**：**13 条给定欠账全部成立、全部可定位到行；其中 8 条是"数是错的/式子是退化的/证书是坏的"，3 条是"等级或口径标错"，2 条是记账。查漏另抓 5 条，最重要的一条不是新错，而是"一条旧勘误本身查错了目录"（L2）——按 `AGENTS §1`，`0 命中` 若不同时报"工具真的执行过、且查的是哪个目录"，就会生产出这种反噬。**

---

## §1 十三条核实表（逐条：位置 / 原文摘录 / 修正 / 等级）

> 摘录一律**逐字**、≤30 字；修正文本以 §4 的"可粘行"为准（两处一致）。

| # | 位置（`文件:行`） | 原文摘录（逐字，≤30 字） | 判：问题 | 修正文本（要点） | 等级 |
|---|---|---|---|---|---|
| **E1** | `research/v3/species/GP9-RB.md:24` | ``d(π²−4)/a²`（`d=3,a=10⁻¹`：`586.96`）` | **错数**：`586.96` = `d=1` 的值；`d=3,a=10⁻¹` 应为 `1760.88`。`py -3.12`：`3(π²−4)/1e−2 = 1760.8813203`、`1(π²−4)/1e−2 = 586.9604401` | `（`d=3,a=10⁻¹`：`1760.88`；`586.96` 是 `d=1`）` | **[E]**（数错）· 结论不动 |
| **E2** | `docs/FINDING-02-TWO-HALVES.md:83` | `反例 liminf Δ_a = 1/2 < 1 = Δ_cont [T] REPORT:281` | **坏证书**：`REPORT:281` 的例 `‖H_a−H‖=ε_a↛0` 却声称强预解成立 ⇒ 算术自相矛盾（`GN2:47–49`，`[T]`）。`[T]` 级标签不成立 | 换 `GO2:66–67` 族：`Δ_a=ε_a→0`、`liminf Δ_a=0<1=Δ_cont`，且 **A4 强预解收敛 `[T]`**、A5 非范数收敛、A6 Hausdorff `[T]`。该证书已有独立定理：**SEP-1**（`GR1-SPECRES.md:144` 定理 / `:131` 见证 W-A / `:146` 证） | **[T]（换证书后仍成立）**+**[M]**（原证书记账） |
| **E3** | `research/v3/species/GN2-VARIATIONAL.md:41` | `与靶**等价**（量词写清：**对任意** `ι_a`…` | **只证了 `⟹`**：§3.2 命题 GN2-ADD 给 `(iii) ⟹ Δ_a ≥ Δ_cont − O(η_a)`；`⟸`（靶 ⟹ `(iii)`）全文无证明，且**同文件 `:135` 逐字写"二者不可互推"** —— 与 `:41` 的"等价"自相矛盾 | `与靶**不弱于**，但只证了 `⟹`（充分条件）；`⟸` 未证（§2.4 第 2 点逐字「二者不可互推」）⇒ 等级 `[C]`，**不是**等价` | **[C]**（由"等价/`[C]`"降为"充分条件 `[C]`"） |
| **E4** | `research/v3/species/GV2-KAPPAC-ALGO.md:29` | ``BianchiLemma.md:224` 的 `c₀>κ` 已假定` | **引用错位**：`BianchiLemma.md` 全文 **0 命中** `c₀>κ`；`:224` 实为"自编号：BL-1…T1–T4"，与该式无关（工具 `grep` **OK**；检索式见 §5.1 #2；命中 **0**） | `BianchiLemma.md:111`（式 `c₀ ≥ 0.5+γ_ent+ln3`，`:23` 同式）与 `:116`（`γ_ent=1 ⇒ 2.60` 数表行）——由二者立得 `c₀>κ`（`κ=0.5`，`:113`） | **[M]**（引号→行号） |
| **E5** | `research/v3/species/GV2-KAPPAC-ALGO.md:91` | `（`BianchiLemma.md:117`：`γ_ent=1 ⇒ c₀≥2.60`）` | **行号漂一格**：`:117` 是 `γ_ent=1.5 ⇒ 3.10`；`γ_ent=1 ⇒ 2.60` 在 **`:116`**，式在 **`:111`** | `(`BianchiLemma.md:111,116`)`（`:117` 是 `γ_ent=1.5⇒3.10`） | **[M]** |
| **E6** | `research/v3/species/GV2-KAPPAC-ALGO.md:75` | `e^{−c₀q|C|} ≤ |w(C)| ≤ e^{−c₀q|C|}` | **退化**：左右两式**逐字同一** ⇒ 等价于等式 `|w(C)|=e^{−c₀q|C|}`，而命题 AS 只需"下界在 `C₀`、上界对所有 `C`"（`:83` 证明只用 `A₀ ≥ e^{−c₀q|C₀|}` 与 `A_j ≤ e^{γj}e^{−c₀q(|C₀|+j)}`） | 拆成两条：`|w(C₀)| ≥ e^{−c₀q|C₀|}` ∧ `|w(C)| ≤ e^{−c₀q|C|}`；若坚持一式，须带两个 O(1) 常数 `c₁e^{−c₀q|C|} ≤ |w(C)| ≤ c₂e^{−c₀q|C|}` | **[T]/[M]**（命题 AS 的证明不动，仅重述前提） |
| **E7** | `docs/paper-v2/ch20-unconditional-i1.tex:298` vs `research/unconditional/P3-CYLINDER-RESIDUE.md:118`（同 `:110`） | ch20：`The floor is \LT{} with the Ursell` ／ P3：``[C]` — conditional on **(R1)**… **(R2)**` | **等级错位**：`\LT`=`labb{cT}{T}`（`paper-v2/main.tex:38`）⇒ ch20 判 `[T]`；P3 判 `[C]`。本轮裁决（`GU1-I1T.md:5` 逐字"`ch20:298` 的 `\LT` 正确；`P3:118` 的 `[C]` 过保守"）与 `:25`（`(R1)` 已由 `P3-DRESSED-MASS.md` §2–4 交付、`(R2)` 早已 `[T]`）⇒ **`\LT` 正确、`[C]` 过保守** | P3 `:110` 的 `(d) **[C] …` 与 `:118` 的 `**[C]** — conditional on (R1)…(R2)` → **`[T]`**，括注"(R1) 见 `P3-DRESSED-MASS.md` §2–4；(R2) `[T]`（Seiler + `M02-CBAR.lean` sorry=0）；定级依据 `GU1-I1T.md` §0.6"。**ch20 侧不改**（已是 `\LT`） | **[T]**（P3 由 `[C]` 升）· 定级条件已由 `gu1-i1t` 满足 |
| **E8** | `research/v3/species/GK5-CONTINUUM.md:61`、`:277`；`docs/FINDING-02-TWO-HALVES.md:133` | `U3:150–155`（段一贡献 = 0） | **指向不存在的文件**：全仓 **0 个** `U3.md`（`Get-ChildItem -Recurse -Filter U3.md`；工具 **OK**；命中 **0 文件**）。真身 `research/v3/species/UNIF3.md`（247 行），逐字句在 **`:151`**（"任何**固定 β（含 2.5）**上的结论对 `{a→0}` 的子列都不构成约束"），整块 `:149–154` | `UNIF3.md:149–154`（逐字句 `:151`）。**注**：该勘误**已**被 `GN1-NODIP.md:248` 逐字记过一次（"真身是 `research/v3/species/UNIF3.md:149–154`"）⇒ 本条为**复核 + 待落盘**，不是新发现 | **[M]**（引号→文件） |
| **E9** | `docs/paper-v6/sec01-position.tex:179`、`sec06-kl.tex:114`、`sec07-assembly.tex:97` | `$\bar u=e\Delta Nu\approx63.8>1$`（三处同句） | **头条数错**：`63.8` 对不上论文自己的窗（`0.018394` ⇒ `2.5/0.018394=135.914`），也≠`\bar u` 的正确值。仓内裁：**`\bar u(2.5)=33.9785`**（`\bar C=N` 现行口径）、**`Ĉ·\bar u=135.914`**（`Ĉ=4`，且 `=eΔ·β_W` 精确恒等）；`63.8 = eΔ·N·[2I₂/I₁](β_W/2)=63.92`（多一个 `N` ＋ Bessel 辐角减半两处偏离叠加）（`GU1-UBAR.md:1,13,17,111,120`） | `$\bar u=e\Delta Nu=33.98>1$（$\hat C\bar u=135.91$，且 $\hat C\bar u=e\Delta\beta_W$ 精确）` | **[E]**（算术不符）+**[M]**（口径未标）· `>1` 的结论不动 |
| **E10** | `research/v3/species/V32-TODO.md:51–58`（R5） | `§8.4 的「体积依赖 8.3%」必须撤换` | **已记、未落**：登记本身正确（`GH2-LTART.md` 复核 `−6.1%±2.4%` 符号相反、`L_t=20` 给 `−2.2%±1.7%`），但 `.tex` **一处未改**：`sec03-measurement.tex:191–192` 仍写 `8.3\%`＋"does not dilute"，`:203–206`（`\src`）仍以 `RUN1:59` 背书 `8.3\%`，`:304`（warnbox 第 3 条）与 `:362`（v3box）仍以 `8.3\%` 当误差棒 | 贴 `GH2-LTART.md:287–298` 的 `.tex` 就绪句；同步 `:203–206`、`:304–306`、`:361–363` 三处（`GH2:303–307` 已给逐行替换）；`tab:m-results`（`:167–172`）**不改数**，只在 `:174–176` 加表注 | **[M]**（未落盘）· 替换文本 **[E]**（已测） |
| **E11** | `docs/position-current.html` | 标题 `第 236 轮`／页脚 `第 238 轮` | **自标轮次不一致 + 编年缺行**：`:6,:103` 写 236、`:540` 写 238、`:545` 已有 238 卡片、`:567` 已有 237 卡片；编年段 `:506–529` 实有轮次 = 178,179,216–224,227–234,225,236 ⇒ **缺 226、235**，且 225 排在 234 之后（顺序错乱）；`第 239 轮`/`r239` 全仓 **0 命中**（工具 **OK**），`docs/position-r239.html` **不存在**（最新 `position-r238.html`）⇒ **r239 未落** | 见 §3-P2 的清单（标题/页脚统一到同一轮次；补 226/235 或显式写"该轮无结论变化"；新快照 `position-r239.html` + 覆盖 `current`） | **[M]**（记账） |
| **E12** | `docs/paper-final/sec04-survey.tex:103` | `reduces to $\operatorname{ess\,inf}G''=H1$` | **恒等读法是假的**：`ess inf G″`（`U2` 读法，测度上的 `ess inf`）与 `H1`（`U1` 读法，`E_μ[·│ξ=0]` 的量）**类型不同**，且 **`NA2:228` 活见证**：同一行 `⟨Σ│a₀=0⟩=+0.060>0` 而 `P(Σ<0)=0.288>0` ⇒ 平均正、`ess inf` 负。判：作为恒等式 `[T]`（假）；作为"等难度"`[O]`（`GA3-MINVEXP.md:146–151,335`） | `reduces to the same open cell as $H1$ (\code{[O]}: an equal-difficulty statement, \emph{not} an identity --- \code{GA3-MINVEXP.md} \S2.2)` | **[T]**（恒等读法为假）+**[O]**（等难度）· 行 11 的 kill 仍合法杀 `U2`，但不得转嫁 `U1` |
| **E13** | `research/v3/species/GR4-SPECIES7.md:17`（§0 第 4 条） | `过 4/4 的有原则通过者 = `m` 的重标` | **未入册**：该判词（四模型只值 2 bit；过 4/4 ⇒ 只能是 `m` 的重标）是全项目**复用率最高**的一条筛子判据，却只出现在 GR4 正文（源：`B03.md:145`、`B05`），`docs/new-tools-registry.md` 无条目 | 写入 `docs/new-tools-registry.md`（一行，见 §4-E13）；并给"任何候选量第一次过四模型时必须同时报「它是否是 `m` 的单调函数」" | **[T]**（判词本身）+**[M]**（登记缺失） |

---

## §2 查漏结果（先给 12 条指定检索式的四项，再给新增欠账）

### 2.1 指定检索式（工具：`pwsh` + `Select-String`；扫描面 `*.md,*.tex,*.html`，**2527 文件**；状态全部 **OK**）

| 检索式（逐字） | 命中行 | 命中文件 | 说明 / 落点 |
|---|---:|---:|---|
| `662×|662x` | **84** | 38 | 现读数已改 `165.5×`（`GW1-FLOW-TARGET.md:22` 判 `[M]` 口径错）；残留见 **L6** |
| `165\.5` | **46** | 28 | 现行同口径读数（`FINDING-01:92` 逐字"距离 165.5×（不是 662×）"） |
| `16\.99` | **48** | 11 | `\bar C=1` 槽的合法格（`GU1-UBAR.md:14`），与 `β₁=0.147152` 成对自洽 |
| `33\.98` | **98** | 27 | 现行正确值 `33.9785`（`PAPER-CPB-CLOSURE*.md:63` 等） |
| `63\.8` | **25** | 13 | 唯一公式性出现 = **E9 那三处**；其余为 `main.tex:690` 的 GEVP 本征值 `63.82`（**不同对象，勿误改**） |
| `4\.84` | **99** | 54 | 已修（`main.tex:243`、`REVISION-v3.1.md:30`）；`ch20:204` 给 `0.018269↔4.84` 正确配对 |
| `4\.87` | **104** | 50 | 已修（`0.018394↔4.87`），`GU2-FIVERINGS.md:158` 另给 `=1/g(a*)` 的互证 |
| `17\.97` | **36** | 21 | 未见自相矛盾（`GC2`/`GX2` 均以其为比值读数） |
| `135\.9` | **91** | 31 | `Ĉ·\bar u=135.914=eΔβ_W`（精确恒等） |
| `U3\.md` | **19** | 11 | **真指 `U3.md` 的 0 行**：4 行在 `GN1-NODIP.md:100,104,238,248`（**记录其不存在**），其余 15 行是 `*-SU3.md` 子串（`LEAD-NOTE-NA-SU3.md` 等）⇒ 见 E8 与 **L1** |
| `BAKRY-EMERY` | **39** | 19 | 两类：`Q4-BAKRY-EMERY.md`（**存在**，`research/v3/`）与 `Q04-BAKRY-EMERY.md`（**不存在**）⇒ 见 **L3** |
| `\[C\] on` | **13** | 6 | 显式条件式共 6 处实指（`H07-RETRACT:20,36,81`·`S5-N-UNIFORM:114`·`G4-GENERAL-BOUND:182`·`S3-UV-FESHBACH:62`），与 `GT2-GRADES.md:103` 的计数一致；其余 7 行是**元讨论**行（`GT2:52,103,115,244,254`、`YANG-MILLS-CONDITIONAL-REDUCTION:112,661`） |

**补两条本文件自己用的检索式（四项齐）**

| # | 检索式 | 工具 | 状态 | 命中 | 判 |
|---|---|---|---|---:|---|
| a | `c₀\s*>\s*κ` @ `research/bianchi/BianchiLemma.md` | `grep` | **OK** | **0** | `BianchiLemma.md` 全文无 `c₀>κ` ⇒ E4 的引用错位成立 |
| b | `Get-ChildItem -Recurse -File -Filter U3.md` | `pwsh` | **OK** | **0 文件** | `U3.md` 不存在 ⇒ E8 成立 |
| c | `(?<![A-Za-z])U3:\d` @ `research/v3/species/*.md` | `pwsh` | **OK** | **27 行 / 5 文件** | 见 L1（先跑的 `U3:\d` 给 32 行，被 `SU3:` 子串污染 ⇒ 加边界后重跑，**旧的 32 不采信**） |
| d | `第 239 轮|r239` @ `research/v3/species/*.md` | `pwsh` | **OK** | **0** | r239 无正文落点 ⇒ E11 |

**`TOOL-FAIL`：0 次。** 一次**过滤器缺陷**（c 的先跑版本，`SU3:` 污染）已记录并重跑；按 `AGENTS §1`，该次结果**作废**，不作为任何结论。

### 2.2 查漏新增欠账（5 条）

| # | 位置 | 事实 | 修正 | 等级 |
|---|---|---|---|---|
| **L1** | `research/v3/species/GR2-EDGEATOM.md:56,117,138,211,224,246,247`；`GU1-I1T.md:22,108,140`；`GT1B-REDUCTION.md:183` | **`U3` 双指**：在 `GK5` 里 `U3`≡`UNIF3.md`（247 行），在 GR2/GU1/GT1B 里 `U3`≡`research/v3/U3-LARGE-N-BINDING.md`（**398 行**）。硬证据：GR2 引 `U3:313,322,344,353` 与 `U3:305–344` —— **全部越出 UNIF3 的 247 行**；而 `U3-LARGE-N-BINDING.md:144` 逐字正是 "`tr F²=tr E²=H` … conserved"，`:313` 正是 `(A1) trF^2 residues = [0.75, 0.0,…]`；`U3-LARGE-N-BINDING:305–344` 正是 §8 witness。**同一行 `GR2:246` 同时出现 `U3:305–344` 与 `UNIF3:145–162`** ⇒ 同文件内两 token 指两文件 | 全仓把 `U3` 的唯一指代写成**全名**：段一贡献 → `UNIF3.md`；2D `f₁²≡0`/`w₁=0`/witness → `U3-LARGE-N-BINDING.md`。**禁止**把 E8 的修正**反向套用**（若照 E8 把 `U3:144` 也改成 `UNIF3.md:144`，会指到"## 4. 判：技术还是本质"，**内容完全错**） | **[M]**（命名碰撞，11 处） |
| **L2** | `research/v3/species/GA3-MINVEXP.md:126–132,335`（§1.6） | **一条旧勘误查错目录**：GA3 判 "`ess inf G″=H1` 的唯一逐字书证指向不存在的文件"，其核验是 `Test-Path research/v3/species/Q4-BAKRY-EMERY.md = False`。**实际文件在 `research/v3/Q4-BAKRY-EMERY.md`（384 行）**，且 `:50–58` 逐字含 "the marginal BE constant is **literally** `ess inf G″`"（`:53–54`）⇒ `N09-RESTATE.md:76` 的引用**是有效的**。GA3 的"全 `species/*.md` 检索 `literally` → 3 命中"也因此**扫错了目录** | **撤回 GA3 §1.6 的"悬空引注"判定**（`[M]`→作废）；保留其 §2.2 的实质判词（恒等读法 `[T]`（假）、等难度 `[O]`），因该判词的独立见证是 `NA2:228`，不依赖这条书证 | **[M]**（元勘误：撤一条勘误） |
| **L3** | `research/v3/species/B01-RP.md:134` | 引 `Q04-BAKRY-EMERY.md`（**0** 与 **4** 之误）：`Get-ChildItem Q04*` 只回 `Q04-WEYL.md` ⇒ `Q04-BAKRY-EMERY.md` **不存在**（工具 **OK**；命中 **0 文件**） | `Q4-BAKRY-EMERY.md`（`research/v3/`，384 行） | **[M]** |
| **L5** | `docs/paper-v6/sec01:179/sec06:114/sec07:97` vs `docs/paper-final/sec05-ledger.tex:220` | `\bar u(2.5)` 的**槽**跨论文不一致：paper-v6 用错值 `63.8`；paper-final 用 `16.99`（`\bar C=1` 槽，与 `β₁=0.147152` 成对，**自洽**）。两处都"`>1`"，但**不是同一个数**，评审对读会当成矛盾 | 两文各标槽：paper-v6 → `\bar u(2.5)=33.98`（`\bar C=N`）＋ `Ĉ\bar u=135.91`；paper-final → 保留 `16.99` 但写"`\bar C=1` 槽"。**同一张表内不得混槽** | **[M]** |
| **L6** | `B14.md:106,154`；`GAP2.md:54,105` | `662×` 已由 `GW1-FLOW-TARGET.md:22,172` 判 `[M]`（分子 `β_W`、分母 `β_A` 混口径，同口径 `165.5×`），但四处**仍以现读数出现**（`B14:106` 列为表值、`B14:154` 列为出处、`GAP2:54` 列为判据、`GAP2:105` 作为"120×/30×/662×"三兄弟之一） | 已在 `GW1` 记录，但**未回填**：这四处应加"（已改 `165.5×`，`GW1:22`）"或删；`GAP2:105` 的"不产生定理"结论不受影响（`GW1:138` 逐字"不改变 `D04/U13` 任何判词"） | **[M]**（已记录·未回填） |

### 2.3 查了但**不构成**欠账（防止把正常写法当错，逐条留痕）

| # | 位置 | 为什么不是欠账 |
|---|---|---|
| N1 | `GV2:20` 引 ``BianchiLemma.md:24`` 支撑"判据 `κ_c ≥ 0.81`" | `:24` 逐字是 `CPB-W` 需 `κ_c ≲ 0.80`。`0.81` 是它的**反写阈值**（判据取反），不是错引；建议加"（`:24` 的 `≲0.80` 取反）"以免下一轮再被当错（`[M]` 级建议） |
| N2 | `docs/paper-v3/main.tex:690` 的 `63.82` | 是 GEVP 谱值 `228.17,119.79,63.82,33.40` 中的一项，**与 `\bar u` 无关**；`63\.8` 检索的 25 行里它最容易被误改 |
| N3 | `FINDING-02:84` 引 `REPORT:283` 的 `T_a=diag(1,e^{−c₁a},e^{−c₂a})` | 该族是**极限不编码 `Δ`** 的反例（`GL-B:283`），与 E2 要换的是 `:281`（**严格性**例）不同条，**保留** |

---

## §3 清账优先级（三档 + 一档待裁）

### P0 · 论文级 —— 必须动 `docs/paper-v6/` 或 `docs/paper-final/`（3 条）

| 顺序 | 条 | 文件:行 | 为什么是 P0 | 落盘形态 |
|---|---|---|---|---|
| 1 | **E9** | `paper-v6/sec01-position.tex:179`、`sec06-kl.tex:114`、`sec07-assembly.tex:97` | 头条**错数**，且被摘要级句子（"只用来显示窗的边界"）承载；错法可诊断（多一个 `N` ＋ Bessel 辐角减半），改后 `>1` 结论不变 ⇒ **零风险、高收益** | §4-E9 三行；改后重跑 `pdflatex`（`main.pdf` 1557747 B 需重建） |
| 2 | **E12** | `paper-final/sec04-survey.tex:103`（表第 11 行）；同类表述另见 `:131` | 把"等难度"写成恒等式，是**类型错**；`NA2:228` 一行数据即可反证。审稿人一旦拿 `NA2` 对读，第 11 行的 kill 会被读成整行作废 | §4-E12 一行；**只改措辞，不改 kill**（`ess inf G″<0` 杀 `U2` 仍成立） |
| 3 | **E10** | `paper-final/sec03-measurement.tex:191–192` + 同步 `:203–206`、`:304–306`、`:361–363` | 已登记未落盘 ⇒ 论文与复核报告**直接矛盾**（`−6.1%±2.4%` 符号相反 vs `+8.3%`）；这是全台账里唯一"读者拿两份文件一对就炸"的一条 | `GH2-LTART.md:287–298` 的 `.tex` 就绪句 + `:303–307` 的三行同步表；`tab:m-results` 不改数，只加表注 |

### P1 · 工作文件级 —— 只动 `research/`（8 条 + L1~L3）

| 批 | 条 | 落点 | 要点 |
|---|---|---|---|
| B1（数） | **E1** | `GP9-RB.md:24` | `586.96`→`1760.88`，并注明 `586.96` 是 `d=1`；`§0.10` 是"关键数（全可核）"行，错一个数整行信用受损 |
| B2（证书） | **E2** | `FINDING-02:83` | 换 `GO2:66–67`/SEP-1；**同轮**把 `:84` 的 `REPORT:283` 留着（不同条） |
| B3（等级） | **E3** | `GN2:41` | `等价`→`只证 `⟹``，等级 `[C]`；引用同文件 `:135` 的"不可互推"作内部证据 |
| B4（引号） | **E4·E5·E6** | `GV2:29`、`:91`、`:75` | 同一文件三处：两处行号、一处退化不等式；**一起改**，否则 `:75` 的 `(ii)` 与 `:83` 的证明仍对不上 |
| B5（文件号） | **E8** | `GK5:61`、`:277`、`FINDING-02:133` | `U3:150–155`→`UNIF3.md:149–154`；**改前先读 L1**，不要把 `U3:144` 一起改 |
| B6（等级·已裁） | **E7** | `P3-CYLINDER-RESIDUE.md:110`、`:118` | `[C]`→`[T]`，括注定级依据 `GU1-I1T` §0.6；**ch20 不改**（`\LT` 正确，且 `paper-v2` 为存档线） |
| B7（新发现） | **L1·L2·L3** | `GR2`/`GU1-I1T`/`GT1B` 11 处；`GA3:126–132,335`；`B01-RP:134` | L2 是**撤勘误**，须排在 L1 之前做（先把"引文可信"恢复，再改引文写法） |

### P2 · 位置图级 —— 只动 `docs/position-current.html`（1 条 = E11，附清单）

给 `lead` 的**逐项清单**（`AGENTS §6` 要求每次结论变化同时留快照）：

1. **轮次自标统一**：`:6`（`<title>`）与 `:103`（`<h1>`）写"第 236 轮"，`:540`（`footer`）写"第 238 轮"，`:362` 另有"第 236 轮坐标审计"块 ⇒ 三处必须统一到**当前轮**。
2. **编年段补行**（`:506–529`）：现序列 = 178,179,216→224,**227**→234,**225**,236。缺 **226**、**235**；且 `225` 落在 `234` 之后（位置错，应插在 `224` 后）。补法二选一：补一条，或写"第 226/235 轮：无结论变化（只做审计，等级未动）"。
3. **r236 的位置**：`:528` 的 236 条在 `225` 之后，属**排序错**（不是缺）——本轮任务书写"编年段缺 `r236`"，实测是**位置错 + 缺 226/235**，以本行为准。
4. **r239 落盘**：全仓 `第 239 轮|r239` **0 命中**、`position-r239.html` **不存在** ⇒ 若第 239 轮有结论（`git log` 顶端三条 `a734e1f`/`2481fc2`/`190d1e2` 均属论文线），须**新建快照** `docs/position-r239.html` 再覆盖 `current`；否则显式记"第 239 轮：仅论文出版线，位置不变"。

### P3 · 登记级 —— 只动 `docs/new-tools-registry.md`（1 条 = E13）

- 加一行（§4-E13）；触发条件是 `AGENTS §6` 条件 ①（断言等级变化）——四模型判别力上界本就是 `[T]`，故这里只是**入册**，不动位置图。

### P4 · 待裁 —— 先裁决再动（不属本轮可落盘）

- **L1** 的"`U3` 唯一指代"要用一条**命名纪律**收口（例如：`U3` 前缀只许 `U3-LARGE-N-BINDING.md`，`UNIF3` 一律写全名），否则下一轮还会再裂一次。
- **E7** 已由 `gu1-i1t` 定级，可直接落；**但**若 `lead` 让 `paper-v2` 复活为投稿线，则 `ch20:298` 与 P3 的等级必须**同时**落，不能单边改。

---

## §4 可直接粘的修正文本（每条一行；`文件:行` + 旧 → 新）

```
E1  research/v3/species/GP9-RB.md:24
    旧：`d(π²−4)/a²`（`d=3,a=10⁻¹`：`586.96`）
    新：`d(π²−4)/a²`（`d=3,a=10⁻¹`：`1760.88`；`586.96` 是 `d=1`）

E2  docs/FINDING-02-TWO-HALVES.md:83
    旧：反例   liminf Δ_a = 1/2 < 1 = Δ_cont    [T]  REPORT:281（抽象族，非 YM）
    新：反例   liminf Δ_a = 0   < 1 = Δ_cont    [T]  GO2:66–67（定理 SEP-1，GR1-SPECRES.md:131,144；ε_a→0、N_a→∞，A4 强预解收敛 [T]；REPORT:281 的例算术自相矛盾，GN2:47–49，不再引用）

E3  research/v3/species/GN2-VARIATIONAL.md:41
    旧：与靶**等价**（量词写清：**对任意** `ι_a` 在 `η_a→0` 时 (iii) 失效 ⟺ 污染 ⇒ 不弱于靶）
    新：**不弱于**靶，但只证了 `⟹`（充分条件 `[C]`）：`(iii) ⟹ Δ_a ≥ Δ_cont − O(η_a)`（§3.2）；反向 `⟸` 未证，且本文件 §2.4 第 2 点逐字「二者不可互推」

E4  research/v3/species/GV2-KAPPAC-ALGO.md:29
    旧：（这恰好是"平凡 B2（`κ_c=1`）成立"的同一条件，`BianchiLemma.md:224` 的 `c₀>κ` 已假定）
    新：（这恰好是"平凡 B2（`κ_c=1`）成立"的同一条件；`c₀>κ` 由 `BianchiLemma.md:111`（`c₀ ≥ 0.5+γ_ent+ln3`，`:23` 同式）与 `:116`（`γ_ent=1 ⇒ 2.60`）立得，`κ=0.5` 见 `:113`）

E5  research/v3/species/GV2-KAPPAC-ALGO.md:91
    旧：（`BianchiLemma.md:117`：`γ_ent=1 ⇒ c₀≥2.60`）
    新：（`BianchiLemma.md:111` 的式 + `:116` 的 `γ_ent=1 ⇒ 2.60` 行；`:117` 是 `γ_ent=1.5 ⇒ 3.10`）

E6  research/v3/species/GV2-KAPPAC-ALGO.md:75
    旧：(ii) 双侧权标度：`e^{−c₀q|C|} ≤ |w(C)| ≤ e^{−c₀q|C|}` 的首项（即最小曲面 `C₀` 自身的质量 `A₀ ≥ e^{−c₀q|C₀|}`）与上界；
    新：(ii) 双侧权标度（**上下界不同常数，原文左右同一 ⇒ 已退化为等式**）：`|w(C₀)| ≥ e^{−c₀q|C₀|}` ∧ `|w(C)| ≤ e^{−c₀q|C|}`（§1.4 的证明只用这两条）；若写成一式，须写成 `c₁e^{−c₀q|C|} ≤ |w(C)| ≤ c₂e^{−c₀q|C|}`，`c₁,c₂∈(0,∞)` 为 `O(1)`。

E7  research/unconditional/P3-CYLINDER-RESIDUE.md:110 + :118
    旧：:110  `(d) **[C] Polymer identification (P3).**`
         :118  `**[C]** — conditional on **(R1)** … and **(R2)** …`
    新：:110  `(d) **[T] Polymer identification (P3).**`
         :118  `**[T]** — (R1) the dressed-mass/amplitude resummation, delivered at `[T]` in `P3-DRESSED-MASS.md` §2–4; (R2) the full-rep KP estimate, `[T]` (Seiler; `research/v3/lean/M02-CBAR.lean` 25 theorems, `sorry=0`); grading per `GU1-I1T.md` §0.6 (the earlier `[C]` was over-conservative; `ch20:298`'s `\LT` was right).`
    （`docs/paper-v2/ch20-unconditional-i1.tex:298` **不改**：`\LT` = `labb{cT}{T}`，本就 `[T]`）

E8  research/v3/species/GK5-CONTINUUM.md:61 与 :277；docs/FINDING-02-TWO-HALVES.md:133
    旧：`U3:150–155`（段一贡献 = 0）
    新：`UNIF3.md:149–154`（逐字句在 `:151`：任何**固定 β（含 2.5）**上的结论对 `{a→0}` 的子列都不构成约束）
    ⚠ 不得把 `U3:144`／`U3:313`／`U3:305–344` 一并改：那些指 `research/v3/U3-LARGE-N-BINDING.md`（见 L1）

E9  docs/paper-v6/sec01-position.tex:179
    旧：故 $\beta=2.5$ 处的任何读数（$\bar u\approx63.8>1$）只用来显示窗的边界
    新：故 $\beta=2.5$ 处的任何读数（$\bar u(2.5)=e\Delta Nu=33.98>1$，$\hat C\bar u=135.91=e\Delta\beta_W$）只用来显示窗的边界

    docs/paper-v6/sec06-kl.tex:114
    旧：$\beta=2.5$ 处 $\bar u=e\Delta Nu\approx63.8>1$，更无从谈起
    新：$\beta=2.5$ 处 $\bar u=e\Delta Nu=33.98>1$（$\hat C\bar u=135.91$），更无从谈起

    docs/paper-v6/sec07-assembly.tex:97
    旧：$\beta=2.5$ 处 $\bar u=e\Delta Nu\approx63.8>1$，KP 窗外的任何一步本文都没有
    新：$\beta=2.5$ 处 $\bar u=e\Delta Nu=33.98>1$（$\hat C\bar u=135.91=e\Delta\beta_W$），KP 窗外的任何一步本文都没有
    （三处均只改数，`>1` 的判词不变；`GU1-UBAR.md:1,13,17,111,120` 为唯一权威读数）

E10 docs/paper-final/sec03-measurement.tex:191–192 → 贴 `GH2-LTART.md:287–298`（`.tex` 就绪句，逐字）：
    \item \emph{Volume dependence: not resolved, and the earlier $8.3\%$ was statistics.} The pair (0.373,0.342) is a $1.3\sigma$ difference … $w_1(6,12)=0.3354\pm0.0057$ / $w_1(8,12)=0.3558\pm0.0057$ … $-6.1\%\pm2.4\%$, \emph{opposite in sign} … at $L_t=20$ $w_1(6,20)=0.3491\pm0.0036$ / $w_1(8,20)=0.3568\pm0.0045$, i.e.\ $-2.2\%\pm1.7\%$, so that we quote no volume trend and only $|w_1(L_s{=}6)-w_1(L_s{=}8)|<5.4\%$ at $2\sigma$.
    同步三处（`GH2-LTART.md:303–307` 给逐行替换）：
      :203–206 \src{}      → 旧的 $8.3\%$ 归 RUN1:59；新对归 GH2-LTART.md §2–3（w7/gh2_ltA.json、gh2_ltB.json、gh2_pool.log）
      :304–306 warnbox 3   → `is $8.3\%$` 改为 `is not resolved ($|{\cdot}|<5.4\%$ at $2\sigma$, \code{GH2-LTART.md})`
      :361–363 v3box       → `the dependence … is $8.3\%$` 改为 `the volume dependence … is not resolved: $|w_1(6)-w_1(8)|<5.4\%$ at $2\sigma$ after the re-measurement of \code{GH2-LTART.md}`
    `tab:m-results`（:167–172）不改数，只在 :174–176 加表注指向 `GH2-LTART.md`

E11 docs/position-current.html
    :6/:103 的「第 236 轮」与 :540 的「第 238 轮」统一为当前轮；:506–529 编年段补 226、235（或写"无结论变化"），并把 225 从 234 之后移回 224 之后；新建 docs/position-r239.html 快照后再覆盖 current

E12 docs/paper-final/sec04-survey.tex:103
    旧：11 & minimality & $\lambda_1(\nu)\ge\operatorname{ess\,inf}G''$ & \LN & reduces to $\operatorname{ess\,inf}G''=H1$, and the Savvidy tail \cite{savvidy,no,leutwyler} makes $\operatorname{ess\,inf}G''<0$\\
    新：11 & minimality & $\lambda_1(\nu)\ge\operatorname{ess\,inf}G''$ & \LN & reduces to the same open cell as $H1$ (\code{[O]}: an equal-difficulty statement, \emph{not} an identity --- \code{GA3-MINVEXP.md} \S2.2; the Savvidy tail \cite{savvidy,no,leutwyler} makes $\operatorname{ess\,inf}G''<0$, which kills this cell, not $H1$)\\

E13 docs/new-tools-registry.md（新增一行）
    | 四模型判别力上界（**筛子**，不是工具） | 过 4/4 的"有原则通过者" = `m` 的重标；四模型只值 2 bit | `[T]`（`GR4-SPECIES7.md` §0 第 4 条；`B03.md:145`；`B05`） | 用途：任何候选量首次过四模型时，必须同报「它是否为 `m` 的单调函数」；若是 ⇒ 不得当作"第 7 物种/新范式" |

L1  （命名纪律，非逐行改）全仓：`UNIF3.md` 只许写全名；`U3-LARGE-N-BINDING.md` 只许写全名或 `U3-LNB`。
    GR2-EDGEATOM.md:56,117,138,211,224,246,247 → `U3:` 改 `U3-LARGE-N-BINDING.md:`（同 :246 行内的 `UNIF3:145–162` 保持不变）
    GU1-I1T.md:22,108,140 → 同改；GT1B-REDUCTION.md:183 → 同改

L2  research/v3/species/GA3-MINVEXP.md:126–132 与 :335 第 4 条
    旧：**"`ess inf G″` = H1"这条等式的唯一"逐字"书证，指向一个不存在的文件。**
    新：**（**撤回**，本行作废）该书证有效**：`Q4-BAKRY-EMERY.md` 实存于 `research/v3/`（384 行），其 `:50–58` 逐字含 "the marginal BE constant is literally `ess inf G″`"（`:53–54`）；原核验只查了 `research/v3/species/`（路径错，非文件缺失）。**实质判词不变**（`[T]` 假 / `[O]` 等难度），独立见证 `NA2:228`。

L3  research/v3/species/B01-RP.md:134
    旧：`Q04-BAKRY-EMERY.md`
    新：`Q4-BAKRY-EMERY.md`（`research/v3/Q4-BAKRY-EMERY.md`）

L6  B14.md:106,154 · GAP2.md:54,105
    旧：`662×`（作为现读数）
    新：`165.5×`（同口径；`662×` 是 `β_W/β_A` 混口径读法，`GW1-FLOW-TARGET.md:22` 已判 `[M]`）
    ⚠ 两文件的**结论**不动（`GW1:138`：不改变任何判词）
```

---

## §5 台账（`AGENTS §1` 四要素：工具状态 / 实际检索式 / 命中数 / 失败标记）

### 5.1 工具状态表

| # | 工具 | 状态 | 实际调用（逐条） | 命中 |
|---|---|---|---|---|
| 1 | `read` | **OK**（全部返回正文） | `GP9-RB.md`(335 行，全文) · `GN2-VARIATIONAL.md`(316，全文) · `GV2-KAPPAC-ALGO.md`(298，全文) · `GK5-CONTINUUM.md`(381，全文) · `V32-TODO.md`(59，全文) · `BianchiLemma.md:100–139,205–235` · `FINDING-02-TWO-HALVES.md:75–90,125–140` · `GR4-SPECIES7.md:1–45` · `GU1-I1T.md:1–30` · `GA3-MINVEXP.md:126–155,330–340` · `Q4-BAKRY-EMERY.md:46–60` · `GH2-LTART.md:275–312` | 全出 |
| 2 | `grep`（ripgrep） | **OK**（有效响应，非 provider failure） | `SEP-1|SEP1|SEP–1` @ `*.md` ⇒ **4**（全在 `GR1-SPECRES.md:73,131,144,431`）；`c₀\s*>\s*κ` @ `BianchiLemma.md` ⇒ **No matches（0）** | 4 / 0 |
| 3 | `pwsh` + `Get-Content` 定段 | **OK** | `ch20:290–310` · `P3-CYLINDER:110–125` · `paper-v6/sec04-delta.tex:103` · `paper-final/sec04-survey.tex:103` · `GO2-SOFTMODE.md:55–150` · `UNIF3.md:140–160` · `U3-LARGE-N-BINDING.md:1–12,140–160,305–325` · `sec03-measurement.tex:188–208` · `position-current.html:495–587` | 全出 |
| 4 | `pwsh` + `Select-String`（12 条指定式 + 2 条自用式） | **OK** | 逐条见 §2.1；扫描面 `*.md,*.tex,*.html` = **2527 文件** | 见 §2.1 |
| 5 | `pwsh` + `Get-ChildItem`（存在性核验） | **OK** | `-Filter U3.md` ⇒ **0 文件**；`-Filter Q04*` ⇒ **1**（`Q04-WEYL.md`）；`-Filter Q4*` ⇒ **3**（`research/v3/Q4-BAKRY-EMERY.md`、`species/Q4-LOGCONVEX.md`、`w7/q4_bakry_emery.py`）；`position-r*.html` ⇒ **51**（最大 `position-r238.html`，**无 r239**） | 0/1/3/51 |
| 6 | `pwsh` + `git log --oneline -12` | **OK** | 见 §2 尾注；顶端 `a734e1f / 2481fc2 / 190d1e2`（论文线）、`e35ac6f`(r238)、`d2900e3`(r237) | 12 |
| 7 | `py -3.12` | **OK**（`exit 0`，确定性，无随机、无 MC，<1 s） | `3(π²−4)/1e−2` = `1760.8813203268073`；`1(π²−4)/1e−2` = `586.9604401089358`；`1−4/π²` = `0.5947152654306489`（与 `GP9-RB:87` 逐位一致 ⇒ 该文件其余常数可信） | — |
| 8 | `web_search`/`web_fetch` | **未调用** | — | 本文件无文献断言；`Seiler`/`Reed–Simon` 一律转引仓内既有判词，按仓内口径记 `[R]`/`[T]` |
| 9 | `.lean` | **未编译**（只读纪律） | 仅转引 `research/v3/lean/M02-CBAR.lean` 的 `sorry=0`（经 `GU1-I1T:25`），**未复跑** `check-*.sh` | — |

### 5.2 `TOOL-FAIL` 与过滤器缺陷

- **`TOOL-FAIL`：0 次。**
- **过滤器缺陷 1 次（已作废并重跑）**：`U3:\d` @ `species/*.md` 首跑 **32 行 / 6 文件**，被 `SU3:` 子串污染（`LEAD-NOTE-NA-SU3.md`、`COMPUTE-GPP-SU3.md` 等）。加负向后视 `(?<![A-Za-z])U3:` 重跑 ⇒ **27 行 / 5 文件**。**32 那次不采信**（`AGENTS §1`：过滤器错与"没有"是两件事）。
- **同类反例（本文件自己也踩）**：`U3\.md` 的 19 行里，**15 行**是 `*-SU3.md` 子串 —— 若只报"19 命中"会得出"`U3.md` 到处被引"的相反结论。这正是 L2 的机制（扫错目录 ⇒ 报"不存在"）。

### 5.3 诚实边界（8 条）

1. **本文件只记录**：未改任何既有文件、未 commit、未动 `position-current.html`（`AGENTS §6` 的位置更新**留给 `lead`**，清单见 §3-P2）。
2. **E7 的等级是"引用已裁结论"**：`[C]`→`[T]` 的依据是 `GU1-I1T.md:5,25` 的裁决，不是我重证了 `(R1)`；`(R1)` 的 `[T]` 依赖 `P3-DRESSED-MASS.md` §2–4，本轮**未逐行重读**该文件。
3. **E2 的新证书只核了"够用"，未核"最强"**：`GO2:66–67` 的 A4/A7 我按 `GO2` 自陈 `[T]` 采信（其证明是两行对角矩阵算术，我逐字核过 `‖D_aψ‖² ≤ C_z(|ψ_{N_a}|²+Σ_{n>N_a}|ψ_n|²)` 这一步成立）；是否还有更强的族（例如同时满足 `‖·‖` 有界）**未查**。
4. **E9 的修正是"槽内正确"**：`33.98` 是 `\bar C=N` 现行口径的值；若论文改用 `\bar C=1`，同一位置应为 `16.99`（`GU1-UBAR.md:14`）。**改前须确认论文自己在 §2 用的是哪一格**（`paper-v6/sec02-i1.tex:104–130` 已裁 `\bar C≡N`）。
5. **L1 的行号判据是"越界"**：我用"`U3:313` 越出 `UNIF3.md` 的 247 行"作硬证据；这是**充分**证据（行号不可能落在文件外），但不是"每一条 `U3:` 都已逐条对号"——11 处里我逐处看了 7 处（`GR2:56,117,138,211,224,246,247`、`GU1:22,108,140`、`GT1B:183` 的引文语义），**行号未逐一到目标行逐字比对**。
6. **L6 是"已记录·未回填"**，不是新欠账；把它列进 §2.2 是为了让 `lead` 能一次清完，而**不是**指控 `B14`/`GAP2` 有新错。
7. **`17.97` 与 `135.9` 两条检索式**：命中数很大（36 / 91），我只做了"是否自相矛盾"的**抽样**（`GC2`/`GX2`/`GU1-UBAR`/`PAPER-CPB-CLOSURE` 四处），**没有**逐行判 127 行。故"未见矛盾"是**定向面上的 0**，不是全称命题。
8. **No QED；Clay = 0%。** 本文件的全部产出是 18 条记账与 18 条可直接粘的修正行；**没有一条**把任何断言升为 `[T]`（唯一"升级"E7 是引用他轮已裁的等级变更），也**没有一条**触及质量隙本身。

---

*Status：13 条给定欠账**全部成立**（8 条是错数/错式/坏证书，3 条是等级或口径，2 条是记账），查漏另得 5 条（`[M]` 级）。最硬的三条：E9（论文头条 `63.8`→`33.98`，多付一个 `N`）、E12（`ess inf G″=H1` 恒等读法为假，`NA2:228` 活见证）、E2（`REPORT:281` 证书坏，换 SEP-1）。最值得记的一条是 **L2**：一条旧勘误自己查错了目录（`species/` vs `v3/`），把一条有效引文判成"指向不存在的文件"——`AGENTS §1` 的"四项齐"正是为这件事写的。**No QED；Clay 0%。***

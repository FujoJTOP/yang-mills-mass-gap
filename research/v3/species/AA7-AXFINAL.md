# AA7-AXFINAL — AX 的最终形态与它到 Clay 的账

> 轮 AA7 · `agent-114` · `D:\Dev\math`。写域仅本文件；既有文件一字未改；未 `commit`；零新 MC；零 Lean。
> **编号**：按题面（OS2=RP、OS3=`O(4)`）；仓内 `docs/conditional-theorem.md:45` 反号（OS2=欧氏旋转、OS3=RP）⇒ 跨文件必带号。
> **口径**：`Clay ⟸ AX ∧ NT ∧ GAP`（`research/unconditional/B2-bypass-gap.md:6` 逐字：`Clay(μ*) = AX(OS0–4) ∧ NT(Λ>0, non-Gaussian) ∧ GAP(Δ>0) ∧ [Wilson construction]`）；重构只兑 `AX`。**只出 `[T]`/`[NG]`**。
> **窗**：OS0 的 `L` **固定**（`UVS-L` 字面不可满足，`AA1-OS0.md:1.7`）；OS3 在 `a→0`、`L` 固定；OS4/W1 的定量路线在物理窗 `β≥2.5`（`AA3-OS4.md:53`）。跨窗不搬结论。

## §0 判词（四条，各带下一手）

1. **[T] 产出 = (a) Wightman 场 + 三限定**。五行产出全部核到 `file:line`（§1）；但第 4 行（Poincaré 表示/谱条件）的等级是**标准定理级**，不是机器级——`lean-infra/REPORT.md:486` 逐字是"Level B（不做）"。→ **下一手**：AX 的接口签名分两列「机器已证 / 标准定理」，禁混写。
2. **[NG]「AX = OS0∧OS1∧OS2∧OS3∧OS4∧重构 是 AX 的最小充分分解」——假。** 开义务上的最小形是 **`{L3} ∧ {GAP}` + 三条零代价件**（§2）；`SYM-4D` 是独立**陈述**但同**物种**。→ **下一手**：公理表压成两行 `L3 | G2`，五条只留"哪条承重"注。
3. **[T] OS0 不可省，只能换层级**：它是全表唯一 `∀n` 的**增长/温性**条，也是 Minlos 的唯一入口；其余四条（不变/正性/衰减）没有一条能产 `(U1)` 的 `a`-一致常数。省它的合法方式只有把靶改到 Minkowski 侧、让 `W0/W1` 承担温性（`docs/yang-mills-prior-work.md:46–47`）。→ **下一手**：删掉"OS0 不必要"的任何读法，统一写「OS0 = UVS = L3 = Minlos 面」。
4. **[T] OS4 在 `∧GAP` 下其聚类半冗余、其唯一真空半不冗余**；两向要分别写死（`GAP⇒聚类` 是定理，反向不是）。→ **下一手**：台账禁止在「AX-OS4」与「GAP」两处重复记同一个物种（§3 边 E3）。

## §1 ① 逐字核：重构的产出（5 行）+ 三限定（3 行）

| # | 产出 | 逐字 | `file:line` | 判 |
|---|---|---|---|---|
| 1 | `ℋ_phys` | 「`(ℋ, H, Ω)` 为 OS 重构给出的物理 Hilbert 空间、Hamiltonian、唯一真空」 | `research/mass-spectrum/REPORT.md:48` | `[T]` |
| 2 | `H ≥ 0`、`HΩ=0`、`ker H=ℂΩ` | 「`H ≥ 0`，`HΩ = 0`，`ker H = ℂΩ`」 | 同 `:48` | `[T]` |
| 3 | Wightman 场 | 「OS 重构定理：OS0–OS4 ⇒ 存在 Wightman 理论」（把 Euclidean 测度翻成 Minkowski QFT 的标准桥） | `docs/yang-mills-prior-work.md:64` | `[T]`（文献级） |
| 4 | Poincaré 表示 + 谱条件 + 局域性 + 场为算子值分布 | 「W0：可分 Hilbert 空间 ℋ 上有幺正的 Poincaré 表示；谱条件 `spec(P) ⊂ V₊`；存在唯一的 Poincaré 不变真空 Ω」；重构侧盒逐字「I-RECON: OS0–OS4 ⇒ Wightman + 谱条件 + 唯一真空 + Δ>0」 | `yang-mills-prior-work.md:46`；`research/gluing/REPORT.md:1348–1350` | `[T]`（文献/标准）；**非机器** |
| 5 | 机器化的真实边界 | 「Level B（不做）：`TransferData` → Wightman 理论。需要：多元复分析（管域定理、边值分布）、Laplace 变换、不定号正交群 `O(1,3)`。mathlib 全无（G-8/G-9）」 | `research/lean-infra/REPORT.md:486` | `[T]`（缺口记录） |

**三限定（逐字）**：① **所选相** —— AN1 =「所选相的关联满足 OS 公理 E0…E4…」（`research/handoff-audit/16-input-audit.md:31`）；② **子列** —— 「每一 `μ^{(k)}` 满足 OS0–OS4（含 E4 聚类）并重构出 Wightman 理论」，只给子列的唯一原因是"没有极限点唯一的论证"（`research/uniqueness/REPORT.md:78,80–85`）；③ **规范不变场** —— 「continuum OS reconstruction needs a gauge-invariant algebra; gauge-fixed limits do not supply it」（`research/unconditional/T3-construction.md:11`；另 `mass-spectrum/REPORT.md:48` 固定"**规范不变**的 smeared 算符 `O`"）。

**两条核出的错位（本文件新增）**：
- **(A) `Δ>0` 被塞进了重构的产出。** `gluing:1348–1350` 的盒把 `Δ>0` 写在 `I-RECON` **输出**里；正确形态是 `docs/conditional-theorem.md:301` 逐字「OS0–OS4 **+ 隙** ⇒ Wightman 理论 + 谱隙」⇒ 重构不含隙，`Δ>0` 是并列合取项。反例：无质量自由标量过 OS0–4 + 唯一真空 + 聚类而 `Δ=0`（`B2-bypass-gap.md:10`）。判 `[NG]`（"OS0–OS4 ⇒ Δ>0"这一读法）。
- **(B) 第 4 行的 `[T]` 与 `lean-infra:486` 不同级。** 该行逐字是"Level B 不做"⇒ 引用它是**诚实的缺口记录**，不是证明；`docs/conditional-theorem.md:317` 同记「S12 是散文 + 条件」。⇒ 第 4 行只能记"标准定理 + 未机器化"，不得写成"已 Lean 化的产出"。

## §2 ② AX 的最小充分分解

**“可省”只指从开义务表删，不指从陈述删**——从陈述删任何一条都换靶（每行给见证）。

| 件 | 形态 | 开义务 | 删它 ⇒ 换靶？ | `file:line` |
|---|---|---|---|---|
| **OS0** | `UVS(L)`：`(U1)` a-一致 `L^p` 界 + `(U3)` `n!` 增长 + Minlos | **`[NG]`** = L3 | **否（不可删）**：唯一 `∀n` 条件；`(U1)` 供给 0；`Σe^kδ_k` 证明弱收敛只保"是分布"不保温性 | `AA1-OS0.md:8,19,24–27`；`glue-a-limits/REPORT.md:58–62,88–92,109–119`；`GGG3-OSCLOSED.md:17` |
| **OS1** | 平移（群像 `aℤ⁴` 稠密） | **`[T]`**（零代价） | 是（删则无 `E(4)`、Poincaré 表示失据） | `GGG3:18` |
| **OS2** | 反射正性（弱闭凸锥，**不需一致矩界**，Lean `rp_transfer`/`rp_l2_closure_trunc`） | **`[T]`**（零代价） | 是（删则 `ℋ_phys` 无内积正性） | `GGG3:19`；`glue-a-limits:145–177` |
| **OS3** | `O(4)` = `KIN-2[T]` + `SYM-4D[NG]` | **`[NG]`** | **否**：`ℓ=4` 各向异性 Gaussian 过其余四条、过隙、仍破 `O(4)` | `AA2-OS3.md:11,20,50`；`G1-OS2-INFVOL.md:9`；`G1-ASSEMBLE-100.md:52–54` |
| **OS4** | 聚类 + 唯一真空（`ker H=ℂΩ`） | 聚类半 **`[T]`⟸GAP**；唯一真空 **`[T]`**（VAC）；极限唯一 = W1 **不在本条** | 是（删则 `H` 的真空无名、`Ω` 不唯一） | `conditional-theorem.md:294,301`；`GGG3:21,30–33`；`uniqueness:66,78,85` |
| **重构** | `Re`：OS0–OS4 ⇒ `(ℋ,H,Ω)` + Wightman | **`[T]`**（标准，无额外假设） | 是 | §1 行 1–4 |

**判（三条，逐条带下一手）**：
1. **[T] 逻辑上五条合取最小（每条不可删），义务上不最小。** 开义务压缩为 `AX = {L3} ∧ {SYM-4D} ∧ {GAP}` + 三条零代价件（OS1、OS2、`Re`）；而 `SYM-4D` 的内容 = `(a,L)`-一致矩界（=L3 物种）+ 盒标度面的聚类型输入（=`G2`）——逐字：「OS2/SYM-4D is a genuinely separate open **statement** … yet its content is the L3 species, plus a G2/MG-adjacent clustering input」（`G1-ASSEMBLE-100.md:143`，另 `:57,81`）。⇒ **AX 的最小开义务集合 = `{L3, GAP}`**（`:132,145`）。→ **下一手**：公理表两行，`SYM-4D` 单列「独立陈述、同一物种」。
2. **[NG]「OS0 可省略」；[T]「OS0 可换层级」。** 类型论证：OS0 是唯一对 `n` 与对 Schwartz 半范数量化的**增长**条件；其余四条都是固定 `(f,g)`/固定锥上的不变性·正性·衰减，不能合成 `(U1)` 的 `a`-一致常数（`AA1:19` 判 `(U3)` 只是 `(C-form)` 的推论、`(U1)` 供给 0）。换层级的唯一合法形式：靶改到 Minkowski 侧，温性由「场是算子值**温**分布」（`yang-mills-prior-work.md:47` W1）承担——义务不减少，只换名字（Minlos 那一步消失，场分布性那一步出现）。→ **下一手**：把 OS0 写成「Minlos 面」，并在接口注"层级 A（Euclidean）/层级 B（Minkowski）二选一"。
3. **[NG]「OS4 可省」但 [T]「OS4 的聚类半在 `∧GAP` 下冗余」。** `GAP`（`∃δ>0 ∀L,a` 的转移矩阵界）⇒ 指数聚类已 Lean（`GlueA.spectral_decay`，`conditional-theorem.md:294`；`glue-a-limits:15,278–285`），方向**不可逆**（2D 临界 Ising：唯一、聚类、幂律、无隙）。⇒ 在合取 `AX ∧ GAP` 里 OS4-聚类是被蕴含项：**记账上不能两处都算**。而"唯一的理论"要 `W1`（`W1 ⟺ UNI`），`W1` 率任意、不被 `GAP` 蕴含（`MG ⇒ W1`，`W1 ⇏ MG`，`GGG3:33`；`AA5-W1.md:10–11`）。→ **下一手**：`AX` 的 OS4 格只留"VAC 免费 + 聚类由 GAP 供"，W1 单列唯一性表。

## §3 ③ 三件接口表（哪件的输出进哪件的输入）

| # | 边 | 输出（源） | 输入（汇） | 内容 / 为什么 | 判 | `file:line` |
|---|---|---|---|---|---|---|
| E1 | AX → NT | `μ`（极限测度） | NT 的见证 | 只需 `μ` 与短程；见证取 `g_ren≠0` / `Λ>0`（自由理论 `trF²` 的连通 4 点≠0 是**空洞**见证） | `[T]` | `nontriviality/REPORT.md:11,13,23`；`GGG5-RECON.md:25` |
| E2 | AX(OS1∧OS2=RP) → GAP | `ℋ_gi`、`Ω`、`T=e^{−aH}` | GAP 的**陈述** | GAP 逐字 `(GAP) ∃δ>0 ∀L,a: x⊥Ω ⇒ ‖Tx‖≤e^{−aδ}‖x‖` ⇒ 只吃**算子框架**，不吃 `μ` 的隙 | `[T]` | `glue-a-limits:250,253–257` |
| E3 | GAP → AX(OS4) | 一致隙 | 指数聚类 | `GlueA.spectral_decay`；速率 `δ/2`（坐标方向技巧），常数与 `L` 无关只对有界类免费 | `[T]` | `conditional-theorem.md:294`；`glue-a-limits:32,278–285` |
| E4 | GAP ↛ NT | — | — | Banks–Zaks 见证：`Λ>0 ∧ Δ=0`；NT 也可由 UV 侧 `E2-g` 供给、**不经** MG/禁闭 | `[T]` | `B2-bypass-gap.md:11`；`gluing/REPORT.md:1353` |
| E5 | NT ↛ GAP | — | — | 无质量自由标量过 OS0–4 + 唯一真空 + 聚类而 `Δ=0`；Gaussian 亦过 OS0–4 | `[T]` | `B2-bypass-gap.md:10`；`GGG5-RECON.md:25` |
| E6 | AX → `Re` → Wightman；GAP 必走**同一** `Re` | OS0–OS4 | `H≥0`、Poincaré rep、谱条件、`spec(H)∩(0,Δ)=∅` | `Δ` 只有经同一 `(相, 子列, 规范不变场)` 三元组的重构才成为 `H` 的谱性质（Källén–Lehmann）；三元组不同则 `Δ` 与 `H` 不是同一对象 | `[T]` | `conditional-theorem.md:301`；`mass-spectrum:48,60–68`；§1 三限定 |
| E7 | AX(OS3) ⊥ GAP | — | — | 有隙不必 `O(4)`（有质量 `ℓ=4` Gauss 有隙、`W`-不变、非 `O(4)`）；`O(4)` 也不逼出隙 | `[T]` | `AA2-OS3.md:50` |
| E8 | W1/UNI ← ？ | — | — | `MG ⇒ W1`、`W1 ⇏ MG`；`W1` 的率**任意** ⇒ 不是 GAP 的输出；且 `W1 ⟺ UNI` ⇒ 不比唯一性便宜 | `[T]` | `GGG3:33`；`AA5-W1.md:10–11,15` |

**接口一句话**：`AX` 给 `NT` 与 `GAP` 的是**同一个对象 `μ`**，但给的东西不同——给 `NT` 的是"短程耦合非零"的载体，给 `GAP` 的是"`T=e^{−aH}` 这个算子框架"；`GAP` 反过来只回给 `AX` 一样东西：**OS4 的聚类半**（E3），此外不供给任何公理。三件的汇合点只有一个：**E6 的同一重构**。

## §4 记账 / 限定 / 下一手

- **新增判**：(i) `Δ>0` 不在 `I-RECON` 输出内（§1-A）；(ii) 第 4 行产出是标准定理级（§1-B）；(iii) 五条合取义务上可压成 `{L3} ∧ {GAP}`（§2-1）；(iv) OS0 只能换层级、不能省（§2-2）；(v) OS4-聚类在 `∧GAP` 下冗余、唯一真空半不冗余（§2-3）；(vi) 三件接口 8 条边、汇合点唯一（§3）。
- **未做（明示）**：`SYM-4D` 的 `(a,L)`-一致 `a²` 系数；`L3` 的 `(★)`；`(U1)` 的大场侧；`W1` 的弱耦合 `∀ξ` 机制；`hMinlos`（Mathlib 无）。**"`L3=H1=G2` 一个物种"按 `[NG]`（作为定理）记**——仓内自陈 "No claim that solving L3 solves G2 … same functional form, not literally the same statement"（`G1-ASSEMBLE-100.md:154`）；作为**分类记账**引 `:145`。
- **不动**：`docs/`、`research/unconditional/*`、`docs/position-current.html`（AGENTS §6 由 `lead` 回填）。零新 MC、未 commit、未改既有文件。

## §附 工具状态（AGENTS §1 四要素）

| 工具 | 状态 | 检索式 / 命令（逐字） | 命中 |
|---|---|---|---|
| `grep` | **OK**（有效响应，无 `provider failure`/`ENOENT`/空响应） | ①`OS0\|OS1\|OS2\|OS3\|OS4\|UVS\|IN-2\|GL9` @ `*.md` ②`AX\|重构\|reconstruct\|Wightman\|ℋ_phys\|Poincar` @ `research/v3/species/*.md` | ①**1061**（显示 250，全表落盘）②**535**（显示 250） |
| `read` | **OK**（带行号、非空、无截断） | `GGG5-RECON`(47)·`AA1`(74)·`AA2`(75)·`AA3`(81)·`AA4`(69)·`AA5`(89)·`GGG3`(53)·`G1-ASSEMBLE-100`(157)·`mass-spectrum/REPORT:30–74`·`uniqueness/REPORT:50–119`·`16-input-audit:20–79`·`nontriviality/REPORT:5–34`·`lean-infra/REPORT:472–511`·`T3-construction`(17)·`B2-bypass-gap`(31)·`conditional-theorem:40–61,294–321`·`gluing/REPORT:1335–1379,1612–1629`·`yang-mills-prior-work:36–79` | 18 文件 |
| `pwsh` | **OK**（exit 0；首次输出长、已落盘） | `Get-ChildItem -Force D:\Dev\math\research\v3\species`；`Get-ChildItem ...\AA*.md`；`Test-Path ...\AA7-AXFINAL.md` | `AA*.md` = **5**（AA1–AA5）+ AAA1–AAA5；`AA7-AXFINAL.md` = **False**（新建不覆盖） |
| `glob` | 未用（路径发现由 `pwsh Get-ChildItem` 覆盖） | — | — |

**失败标记：0**（无 `TOOL-FAIL`；无 `0 命中` 结论）。**直读 vs 转引**：§1 五行、§2 六行、§3 八边全部逐字直读上列文件；`G1-OS2-INFVOL:9`、`SZ/SZZ` 等均经 `AA2/AA3` 与 `G1-ASSEMBLE-100` 转引，未另开原件。**限定**：本文件只做记账与接口，不含任何新分析估计；"`L3` 真而未证"这一档不在本文件重判（`G1-ASSEMBLE-100:104,152`）。

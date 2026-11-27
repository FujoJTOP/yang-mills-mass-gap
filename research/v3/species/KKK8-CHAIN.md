# KKK8-CHAIN — Clay 的最终最小链（逐环核）

> 轮 KKK8 · 子代理 `kkk8-chain` · `D:\Dev\math`。**写域仅本文件**；既有文件一字未改；未 `commit`；零新 MC；零 Lean。
> **口径（逐字）**：`Clay(μ*) = AX(OS0–4) ∧ NT(Λ>0, non-Gaussian) ∧ GAP(Δ>0) ∧ [Wilson construction]`（`research/unconditional/B2-bypass-gap.md:6`）。
> **编号**：按题面（OS2=RP、OS3=`O(4)`）；仓内 `docs/conditional-theorem.md` 反号（OS2=欧氏旋转、OS3=RP）⇒ 跨文件必带号。
> **窗**：`L3` 在固定 `L`、`a→0`、`g_j→0`；`GAP` 必须 `(L,a)`-**一致**（`glue-a-limits/REPORT.md:341` 逐字"必须关于 `L` 一致"）；`NT` 在连续极限 `μ` 上。跨窗不搬（`AGENTS §2.1`）。

## §0 判词（三条，各带「下一手」）

1. **[T] 最小链 = 三节点、两原子：`{L3} ∧ {GAP} ∧ {NT}`。** 节点 3（陈述账：`B2:6` 三个合取项各自承重），原子 2（物种账：`L3`=UV 面、`G2`=IR 面，`G1-ASSEMBLE-100.md:145,151`）。→ **下一手**：位置图分两列记「节点 3 / 原子 2」，禁把两数混写成一个。
2. **[NG]「NT 可从最小链上删」——假。** `AX ∧ GAP ⇏ NT`（有质量 Gaussian：`Δ=m>0` 且 Gaussian，`AA6-NT.md:52`），`NT ⇏ GAP`（Banks–Zaks：`Λ>0 ∧ Δ=0`，`B2-bypass-gap.md:11`）⇒ `NT` 是独立合取项，删它即换靶。→ **下一手**：`conditional-theorem.md:298` 的 S9（`μ` 非高斯）单列第三行，禁并入 S5/S6。
3. **[T]「NT 不添第三原子」——真。** `NT` 的缺件只有一个 `N2`（`Λ>0`/不塌缩到 `g₀≡0`），经审计归 `G2`（`docs/IN-AUDIT.md:17`，IN-7→G2）；`N1` 近 `IN-0` 且 16 输入 15/16 已兑现（`IN-AUDIT.md:10`）；条件形 `N1∧N2 ⇒ ¬Gaussian` 已 `[T]`（`research/nontriviality/PROOF.md:145–151`）。→ **下一手**：闭 `G2` 的同一台机器须**另起一行**推 `N2`，不得记账为已付。

## §1 最小链（8 步：`{L3} ∧ {GAP} (∧ NT)` → Clay）

| # | 步 | 判 | `file:line` |
|---|---|---|---|
| 1 | `L3` = 矩的 `a`-一致界 + Minlos 阶乘；开放，**唯一** `∀n` 增长条、Minlos 唯一入口 ⇒ 不可省、只能换层级 | `[NG]`（开件）/ 不可省 `[T]` | `HQ1-L3.md:16`；`G1-L3.md:3`；`AA7-AXFINAL.md:12` |
| 2 | `L3 ⇒ UVS(L) ⇒ GL3`（相对紧+对角子列+Minlos）⇒ `OS0`；`OS1` 平移 | `[T]`（给定 `L3`） | `G1-EXISTENCE-THEOREM.md:31,33`；`GGG3-OSCLOSED.md:18` |
| 3 | `OS2`=RP：弱闭锥、不需一致矩界（Lean `rp_transfer`/`rp_l2_closure_trunc`） | `[T]` 零代价 | `GGG3-OSCLOSED.md:19`；`glue-a-limits:145–177` |
| 4 | `OS3`=`O(4)`=`KIN-2[T]`+`SYM-4D[NG]`；`SYM-4D` 内容 = `(a,L)`-一致矩界（`L3` 物种）+盒标度（`G2`）⇒ 不添原子 | `[T]`（分类级） | `AA2-OS3.md:20,26`；`G1-ASSEMBLE-100.md:54,143` |
| 5 | `OS4`：聚类 ⟸`GAP`（`δ/2`，Lean `spectral_decay`）；VAC `[T]`；**极限唯一 = `W1`**（`W1 ⟺ UNI`）不被 `GAP` 蕴含 ⇒ 不入链 | 聚类 `[T]`；`W1` `[NG]` | `conditional-theorem.md:294,297`；`G1-ASSEMBLE-100.md:59`；`AA5-W1.md:11`；`AA7:41,47` |
| 6 | `GAP`：`∃δ>0 ∀L,a: x⊥Ω ⇒ ‖Tx‖≤e^{−aδ}‖x‖` ⇒ `Δ≥δ`；经**同一** `Re` 成 `H` 的谱性质 `spec(H)∩(0,Δ)=∅` | `[NG]`（开件）/ 搬运 `[T]` | `glue-a-limits:250,252`；`AA7:54,58`；`conditional-theorem.md:301` |
| 7 | `Re`：`OS0–OS4 + 隙 ⇒ Wightman + 谱隙`（标准定理，未机器化） | `[T]` 零代价 | `conditional-theorem.md:301`；`AA7:22,29`；`lean-infra/REPORT.md:486` |
| 8 | `NT`：定理 N `N1∧N2 ⇒ ¬IsGaussianMeasure μ`；`N2` 缺 ⇒ 无条件形 | `[NG]`（条件形 `[T]`） | `PROOF.md:145–151`；`AA6-NT.md:16,57`；`IN-AUDIT.md:17` |
| ⇒ | `Clay = AX ∧ NT ∧ GAP`：第 2–5 步兑 `AX`（原子 `L3\|G2`）、第 6 步兑 `GAP`、第 8 步兑 `NT`；汇合点唯一 = 第 7 步**同一** `Re` | `[T]`（结构） | `B2:6`；`AA7:62` |

**两向不可互推（记账铁律）**：`GAP ⇏ NT`（有质量 Gaussian）、`NT ⇏ GAP`（BZ）；`L3 = G2` **只作分类记账** `[T]`（分类级），**作定理** `[NG]`（`G1-ASSEMBLE-100.md:154` 逐字 "same functional form, not literally the same statement"；`AA7-AXFINAL.md:67`）。

## §2 NT 是否可省——裁决

**裁决：三节点、两原子。链是三件，账是两件；任何单写都必须注明是哪一账。**

- **`[T]` 不省节点**：`Clay` 是存在型、但合取项逐字含 `NT`（`B2:6`）；而 `AX∧GAP` 放过有质量 Gaussian（`AA6:52`）⇒ 链上必须有第 8 步（仓内装配链的 S9，`conditional-theorem.md:298`）。
- **`[T]` 不添原子**：`N2 ∈ G2` 原子（`IN-AUDIT.md:17`），`N1` 已由 16 输入兑现 15/16（`IN-AUDIT.md:10`）⇒ `NT` 不构成第三次"新物种"。
- **`[NG]`「同原子 ⇒ `GAP` 买 `NT`」**：同原子 ≠ 同陈述（`G1-ASSEMBLE-100.md:154`）；有质量 Gaussian 是 `GAP ∧ ¬N2` 的见证 ⇒ `N2` 须由闭 `G2` 的同一台机器**另起一行**推出。
- **下一手**：S9 输入栏写 `N1 | N2(G2 原子)`、输出 `¬Gaussian`；位置图 `NT` 行只记「条件于 `N2`，与 `G2` 同原子」，不记「先闭合」。

## §3 接口表（三件：输入 / 输出 / 谁喂谁）

| 件 | 输入 | 输出 | 喂谁 / 不喂谁 |
|---|---|---|---|
| **`L3`** | `(★)` 逐壳单侧方差**上**界 `Var(Pᵀ∇S_j)≤E[PᵀH_jP]+W_j−δ_W` + `(★dec)`，跨 `J∼log(1/a)` 塔；单块 `G″≥δ_W`；`L2` BPHZ | `UVS(L)`=矩的一致界+Minlos 阶乘 ⇒ `μ_L`、`OS0`；`SYM-4D` 的 `(a,L)`-一致面 | → `Re`（`OS0`）；→ `NT` 的载体 `μ` 与 `N1` 短程侧；**⊣ `GAP`**（无蕴含，仅同物种） |
| **`GAP`** | 算子框架 `(ℋ_gi, Ω, T=e^{−aH})`；`(GAP) ∃δ>0 ∀L,a` | `Δ≥δ`；体积一致指数聚类 `δ/2`（GL9 Lean）；经同一 `Re` 的 `spec(H)∩(0,Δ)=∅` | → `AX` 的 `OS4`-聚类半；→ `Re` 的谱隙；**⊣ `NT`**（不给 `N2`） |
| **`NT`** | `μ`（`N1`：短程一致收敛+一圈匹配）；`N2`：极限不塌缩到 `g₀≡0`（`Λ>0`） | `¬IsGaussianMeasure μ`（定理 N） | → `Clay` 第三合取项；**⊣ `GAP`**（BZ 见证） |

**一句话**：`AX` 给 `NT` 与 `GAP` 的是**同一个对象 `μ`**，但给的东西不同（`NT` 取短程耦合载体，`GAP` 取 `T=e^{−aH}` 算子框架）；`GAP` 只回给 `AX` 一样东西（`OS4` 聚类半）；三件的汇合点只有一个 = `Re` 的**同一**三元组（相、子列、规范不变场）。

## §4 工具状态（`AGENTS §1` 四要素）

| 工具 | 状态 | 检索式/命令（逐字） | 命中 |
|---|---|---|---|
| `grep` | **OK**（有效响应，无 `provider failure`/`ENOENT`/空响应） | ①`AA7-AXFINAL\|AA6-NT` @ `research` ②`OS0` @ `research/**/*.md` ③`GAP\|δ>0\|x\^\\perp\|e\^\{-aδ\}` @ `research/glue-a-limits/REPORT.md` | ①**4** ②**167** ③**14** |
| `glob` | **OK** | `research/v3/species/*.md`；`research/v3/**/*.md`；`**/IN-AUDIT.md`；`research/v3/species/KKK*.md` | ①**983**（截 100 落盘）②**1097**（截 100）③**1**=`docs\IN-AUDIT.md`（`research/handoff-audit/IN-AUDIT.md` **不存在**）④**1**=`KKK4-UTILE.md`（`KKK8-CHAIN.md` 先前不存在 ⇒ 新建不覆盖） |
| `read` | **OK**（带行号、非空、无截断） | `AA7-AXFINAL`(79)·`AA6-NT`(70)·`B2:1–20`·`docs/IN-AUDIT`(36)·`nontriviality/PROOF:140–169`·`conditional-theorem:288–309`·`G1-ASSEMBLE-100:140–157`·`HQ1-L3:1–16`·`AA5-W1:1–18` | 9 文件 |
| `pwsh` | 未用（路径发现由 `glob` 覆盖） | — | — |

**失败标记：0**（无 `TOOL-FAIL`）。**直读 vs 转引**：`AA7/AA6/B2/IN-AUDIT/PROOF/conditional-theorem/G1-ASSEMBLE-100/HQ1/AA5` 直读；`GGG3:18,19`、`AA2:20,26`、`glue-a-limits:145–177,250,252,341`、`G1-EXISTENCE:31,33`、`lean-infra:486` 经 `grep` 行命中或所属文件引文，标**转引**。**未改既有文件、未 `commit`、零新 MC。**

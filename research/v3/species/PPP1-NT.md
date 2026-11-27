# PPP1-NT — 把 `NT` 闭成 `[T]`：直路（引理 C + `κ₃≠0`）逐字核 ⇒ `[NG]`；`NT` 无条件 `[NG]`，条件形 `[T]`

> 轮 PPP1 · 子代理 `ppp1-nt` · `D:\Dev\math`。**写域**：新建本文件 + `research/v3/species/ppp1/NT.lean`、`ppp1/NT-C3.lean`；既有文件一字未动；未 `commit`；零新 MC；**未编译 Lean**（只写代码 + 逐行标依赖）；`docs/position-current.html` 交 `lead`（`AGENTS §6`）。
> 标签二值（只 `[T]`/`[NG]`）；无百分数；每条判词带「下一手」。口径承 `AA6-NT:4`：`Clay ⟸ AX ∧ NT ∧ GAP`，`NT = ¬Gaussian`；窗口承 `AA6-NT:6`：`κ₃` 读数在 `SU(2)`、`a=1`、`β=2.5`、有限 `L`，与连续极限 `μ` **不同窗**（`AGENTS §2.1`）。

## §0 判词

| # | 断言 | 判 | 依据（`file:line`） | 下一手 |
|---|---|---|---|---|
| **0.1** | **直路自身即断（不是缺 N2）**：`κ₃(q̄)≠0 ⇒ ¬Gaussian` 用在 `q` 上**假**。`q=¼tr F²` 是 Wick 平方，其**自由**连通三点非零 ⇒ `κ₃(q̄)≠0` 被自由场满足 | **`[NG]`** | `nontriviality/REPORT.md:42–45`（`⟨:F²:(x):F²:(y):F²:(z)⟩_c ≠ 0` 逐字，Gaussian `F`）；`REPORT:55–57`（命题 1.2：见证须 free-vanishing）；`NT-C3.lean:composite_third_cumulant_not_a_witness` | 直路从 NT 路线表删除；`κ₃` 只留诊断栏 |
| **0.2** | **引理 C 真，`n=3` 逆否是 `[T]`**——但只辖**线性**泛函 `⟨φ,J⟩`；`q̄` 不在辖域内。且直路所需的**推理免费**（充分方向取逆否，无需解析性） | **`[T]`**（蕴含） | `PROOF:52–53,57–59`；`lean-infra/…/NonTriviality.lean:78–80`（已机检，无公理）、`:112–123`；`PPP1-NT §1` | 引理 C 只对线性见证引用；对复合算子写「不适用」 |
| **0.3** | `¬Gaussian` **恰好**是 Clay 的 `NT`，不缺更强的物；`GaugeNonTrivial`（=`g_ren≠0`）是**见证形**（更强，唯一合法供货形态） | **`[T]`** | `lean/YMStatement/Basic.lean:52–53,84–87`；`B2:6`；`REPORT:11,55–57` | 位置图 NT 行不改陈述，只改见证栏 |
| **0.4** | `NT` 能否 `[T]`：**无条件 `[NG]`**。唯一完整物 `N1∧N2 ⇒ NT`（定理 N，**归约**） | **`[NG]`**（条件形 `[T]`） | `PROOF:145–151,203`；`KKK5:0.5`；`BREAK-IN-7:18` | NT 留缺件栏；成本记 `IN-3 ∧ Y3` |
| **0.5** | ★「重整化后 `κ₃ = 8.62e−9`，有限且非零」**不成立为见证**：既非 free-vanishing（0.1），也无 `(a,L)` 族（§2a）。它是**诊断读数** | **`[NG]`** | `OO7:0.6`；`QQ2:0.4,51–53`；`AA6-NT:45` | 报 `κ₃` 必带 `(L,β,N,JK)` + 「非见证」 |
| **0.6** | 勘误（对 `AA6-NT:0.3,:42`）：其「蕴含 `[T]`／合格见证类型（自由值 `=0`）」与 `REPORT:11,55–57` 自相矛盾。判词 `[NG]` 不变，理由从「供给不足」升级为「**见证类型不合法**」 | **`[T]`**（勘误） | `AA6-NT:42` vs `REPORT:45,57`；`PROOF:109` | `AA6-NT` 该行加 errata 指针，不改其 `[NG]` |

## §1 引理 C：逐字、`file:line`、辖域

**逐字（`research/nontriviality/PROOF.md:52–53`）**
> **引理 C.** 设 `μ` 是概率测度且 `G_μ` 在 `0` 邻域解析、被矩唯一确定。则
> `IsGaussianMeasure μ  ↔  ∀ n ≥ 3, ∀ f₁,…,f_n : κ_n(f₁,…,f_n) = 0.`

证明 `:57–59`（`⇒`：`log G_μ` 是次数 ≤2 的多项式；`⇐`：解析 + 唯一性延拓 + Minlos/Bochner）。注 `:61`（不能收紧到 `n=4`；反例 `:63–83`）。
Lean：`isGaussianTheory_cumulants_vanish`（`NonTriviality.lean:78–80`，**充分方向，已证，无公理**）、`CumulantCharacterization`（`:112–117`，必要方向=接口）、`isGaussianTheory_iff_cumulants_vanish`（`:120–123`，条件等价）。机检状态与 `PROOF:185–188` 一致。

**它是否真的 `κ₃≠0 ⇒ ¬Gaussian`？** 是——但**只在 `κ₃` 是线性泛函 `⟨φ,f⟩` 的累积量时**（`κ_n(f₁,…,f_n)` 的 `f_i ∈ 𝒮` 就定义了这一点）。两点推论：

1. 这一步**不需要** OS0/OS1：取充分方向的逆否即可（`NT.lean:nthCumulant_three_ne_zero_implies_NT`，`sorry` 为零）。⇒ 直路缺的**不是推理**。
2. `q = 1−½Re tr U_□ = c·a⁴trF²+O(a⁶)` 是 `F` 的**二次**泛函。`REPORT.md:45` 逐字：对 Gaussian `F`，`⟨:F²:(x):F²:(y):F²:(z)⟩_c ≠ 0` ⇒ **自由场的 `κ₃(q̄)` 本来就非零**。故「`κ₃(q̄)≠0 ⇒ ¬Gaussian`」不是引理 C 的实例，而是**另一条假推理**（`NT-C3.lean` 给 `ℚ` 上精确反例：偶源 `κ₃=0`、`Φ=X²` 的 `κ₃=81/32`）。

## §2 直路的三个缺口

**(a) 一致性：要，而且要的是 `(a,L)` 族，不是「再测一个点」。** 已有六格（`QQ2-ALOCAL:51–53`）：`L=6 iso +8.620e−9±3.1e−9`、`L=4 iso +1.037e−7±2.4e−8`、`6³×12 +1.772e−9`、`6³×20 +5.878e−10`、`8³×12 +1.076e−9`、八副本合并 `+6.272e−10±6.8e−11`。**全正**，但级差 `≈176×`（`1.037e−7/5.878e−10`）、无标度律、`a` 全为 `1` ⇒ 无 `a→0` 族。叠 `a¹²` 幂次（`AA6-NT:45`）⇒ 裸读数**不可搬**；要的具名输入 =「重整化三阶累积的 `(a,L)`-一致非零下界」（`NT-C3.lean:C3Uniform`）。

**(b) 号：不只在别处会零——已证它在窗内必过零至少两次。** (i) `χ_q(0)=1/24` 精确（Haar）+ `∫_0^{2.5}χ_q=0.108822>2.5/24` ⇒ `χ_q` 在 `[0,2.5]` 内不可能单调 ⇒ 号必变（`VV2:0.7`；`OO7:0.4`；`QQ2:1(2)` 的 `κ=6` 修正不改结论）。(ii) 精确锚：`χ_q'(0)=0`、`χ_q''(0)=−1/96<0`（`β=0` 局部极大，`DDD1-CHI24:19`）⇒ 由恒等式 `∂_β𝒮_q=−κV²κ₃(q̄)`（`QQ2:0.1(2)`，`κ=6`）得 **`κ₃(q̄)(β=0)=0` 精确**；再与内部峰（`DDD1:25`：`β=0` 峰 + 内部峰夹一极小）合起来 ⇒ `κ₃(q̄)` 在 `(0,2.5)` 内变号 ≥2 次。⇒ `κ₃(q̄)≠0` 是**点值**陈述，零点位置已知非空，无法升成极限陈述。(iii) 对象不可混：`probe-cumulant` 的 `κ₃(V)` 是 `S` 的、`κ₃(q̄)` 是缺陷的（`QQ2:59` 记同款符号坑）⇒ 不互搬（`AGENTS §2.1`）。

**(c) `¬Gaussian` 够不够 Clay 的 `NT`：够，且这是判据不是缺口。** `lean/YMStatement/Basic.lean:52–53` 定义 `NonTrivial μ := ¬IsGaussianMeasure μ`，`MassGapStatement`（`:84–87`）的 `nonTrivial` 字段就是它；`B2:6` 的 `NT(Λ>0, non-Gaussian)` 同一谓词。`g_ren` 是**见证**（`REPORT:55–57` 要求 (i) 规范不变 (ii) 自由值 `=0`），不是更强的 **NT**。⇒ 缺口在「哪个量能当见证」，不在「NT 要求多强」。

## §3 `NT` 能否 `[T]`

**无条件 `[NG]`。** 三条全断：直路 `[NG]`（0.1，蕴含层就假）；`AF+DT ⇒ Λ>0` `[NG]`（`AA6-NT:0.2`；`BREAK-IN-7:3,12`）；`A-abs` 反推 `[NG]`（`AA6-NT:0.4`）。**唯一完整的命题是条件形**：

> **条件陈述（quote-ready）**：`NT` 于 `N1 ∧ N2` 下为 `[T]`（定理 N，`PROOF:145–151`；此为**归约**，`PROOF:203`）；`N2 = IN-3 ∧ Y3`（`BREAK-IN-7:18`）；**无条件形 `[NG]`**。`κ₃(q̄)` 直路**不**是 NT 的证据——`q` 是 Wick 平方，其自由连通三点非零（`REPORT:45`），违反 `REPORT:57` 的见证约束 (ii)；且无 `(a,L)`-一致下界（§2a）。**下一手**：NT 的候选栏只剩 `g_ren`（`GaugeNonTrivial`）；`κ₃` 移诊断栏并附「非见证」标签。

## §4 Lean（本轮**未编译**；每定义/定理的作用与依赖）

`research/v3/species/ppp1/NT.lean`
- `def NT (T) : Prop := ¬ IsGaussianTheory T` — 与 `YMStatement.NonTrivial = ¬IsGaussianMeasure` 同形（`lean/YMStatement/Basic.lean:52–53`；其 `:84–87` 的字段即此谓词）。依赖：`YangMillsInfra.Statement.NonTriviality`（`EuclideanTheory`/`TestFunction` 经 `lean-infra/…/Basic.lean:31,52`）。
- `theorem nt_iff_not_gaussian` — `Iff.rfl`，无依赖额外。
- `theorem lemmaC (h : CumulantCharacterization T)` — 引理 C 全陈述；依赖 `isGaussianTheory_iff_cumulants_vanish`（`NonTriviality.lean:120`）+ 其必要方向接口（`:112`）。
- `theorem nthCumulant_three_ne_zero_implies_NT` — **主定理**：`κ₃≠0 ⇒ NT`；依赖 `isGaussianTheory_cumulants_vanish`（`:78`，已证，无公理）。`sorry` **0**。
- `theorem nt_of_gaugeNonTrivial` — `g_ren≠0 ⇒ NT`；依赖 `gaugeNonTrivial_implies_not_gaussian`（`:202`，已证）。

`research/v3/species/ppp1/NT-C3.lean`（依赖 `Mathlib` + 同上模块）
- `def κ3` + `even_source_kappa3`/`source_m1,m2,m3`/`source_kappa3_zero`/`composite_moments`/`composite_kappa3` — `ℚ` 精确算术，`norm_num`；`sorry` 0。
- `theorem composite_third_cumulant_not_a_witness` — **空洞性**：源 `κ₃=0` 而 `Φ=X²` 的 `κ₃=81/32≠0`；`sorry` 0。
- `structure C3Uniform (T) : Prop`（`c`、`c_pos`、`witness`）— 具名输入的**类型化**。
- `theorem c3_uniform_implies_NT` — 已证（`abs`+`linarith`）；`sorry` 0。
- `theorem hole_bare_to_uniform (T) : C3Uniform T` — **唯一 `sorry`**，docstring 明写「不主张为真，仅作缺口位置」；`#print axioms` 会显示 `sorryAx`，不作公理使用。
- 核验命令（未执行）：`cd research/lean-infra && lake build YangMillsInfra.Statement.NonTriviality`（mathlib `v4.33.0-rc1`）；两个新文件需先入包（或加 `LEAN_PATH`）方可 `lake env lean`。

## §5 台账（`AGENTS §1` 四要素）

| 工具 | 状态 | 检索式（逐字） | 命中／失败 |
|---|---|---|---|
| `grep` | **OK**（带行号，非 `provider failure`） | ①`引理 C\|Lemma C\|κ₃\|kappa_3` @ `D:\Dev\math` ②`IsGaussianMeasure` @ 同 ③`绕过 N2\|引理 C 直\|直证 NT` @ `species` ④`自由值\|Wick 平方\|8σ` @ `research` ⑤`BREAK-IN-7` @ `research` ⑥`κ₃\|6.272\|副本` @ `QQ2-ALOCAL.md` | ①**355**（显示 250，落盘）②**42** ③**1** ④**250/285**（显示上限）⑤**18** ⑥**11**；失败 **0** |
| `read` | **OK**（带行号、非空） | `PROOF.md:40–206`、`NonTriviality.lean`(208)、`Basic.lean`(89)、`NonTrivialCore.lean`(126)、`REPORT.md:30–64`、`OO7-Q3.md`(106)、`AA6-NT.md`(70)、`KKK5-N2.md`(72)、`BREAK-IN-7.md`(18)、`B2-bypass-gap.md`(31)、`lean-infra/Basic.lean`(88)、`lakefile.toml`(13) | 12/12；`TOOL-FAIL` **0** |
| `glob` | **OK** | `research/v3/species/*KKK5*`、`{AA6,OO7,QQ2,VV2,BREAK-IN-7,IN-AUDIT}*`、`**/{BREAK-IN-7,B2-bypass-gap}*`、`research/lean-infra/**/*.{lean,toml}` | 1 / 5 / 2 / 17 |
| `pwsh` | **OK**（exit 0） | `Get-ChildItem … -Filter 'QQ2*'`；`(Get-ChildItem 'species\*.md').Count`；`Select-String`（VV2 行） | **1012**；`QQ2-ALOCAL.md` 13316 B；**1 次 `PATH-MISS`**（`QQ2-ALOC.md` 不存在，真名为 `QQ2-ALOCAL.md`）⇒ 记路径错，非检索失败 |
| `web_*`/MC/Lean 编译/`git`/`skill` | **0**（未调用） | 本靶 = 仓内逐字核 + 写代码；catalog 无匹配技能 | — |

## §6 诚实边界（5 条）

1. **未编译**（任务书明令）：代码按 mathlib `v4.33.0-rc1` 的 API 与仓内既有 `NonTriviality.lean` 的词表写；可能的类型错（如 `obtain` 解构深度、`norm_num` 对 `ℚ` 幂）**未经机器确认**，故 §4 只称「声明与骨架」。
2. **`β=0` 的精确零是格点模型内陈述**（`χ_q'(0)=0` 走 `DDD1:19` 的 Haar 计算，只有经 `QQ2:0.1(2)` 的恒等式才转成 `κ₃(q̄)(0)=0`），**不搬到连续极限**；它的作用只是证明「`κ₃(q̄)≠0` 不是可搬量」。
3. **`composite_third_cumulant_not_a_witness` 是 `ℚ` 上有限反例**，不是 Gaussian 测度的形式化；Gaussian 侧的权威依据是 `REPORT:42–45` 的逐字式。
4. **勘误只动措辞不动判词**：`AA6-NT:0.3/:42` 的 `[NG]` 保留；改的是「蕴含 `[T]`、自由值 `=0`」这一句。
5. 零新 MC、零 `web`、未 `commit`、未动 `docs/position-current.html`；本文件只写 `research/v3/species/PPP1-NT.md` 与 `ppp1/{NT,NT-C3}.lean`。

*PPP1-NT · `[T]` 三条（引理 C 辖域、`¬Gaussian = NT`、引理 C `n=3` 逆否）+ `[NG]` 三条（直路、无条件 NT、`κ₃` 见证）+ 勘误一条。出口：**`NT` 无条件 `[NG]`；条件形 `N1∧N2 ⇒ NT` 为 `[T]`；`κ₃(q̄)` 直路在蕴含层即断**。*

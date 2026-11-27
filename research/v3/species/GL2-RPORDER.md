# GL2-RPORDER — `(RP) ⇒ k₀ ≤ 2` 存在吗？**一般形式不成立（CP 1997 反例）；成立的是"正锥刚性" `c_L ≥ 0 ⇒ k₀ = 2`；`2/(m+2)` 是"RP-可证类"的前沿，不是 RP 的前沿**

> 第 GL2 轮 · 子代理 `gl2-rporder` · 工作目录 `D:\Dev\math`。
> **只写本文件**；未改任何既有文件；未 `git commit`；**无 MC**（本轮数值全部确定性：符号线性化 + 最小二乘拟合 + 线性规划，零随机数）。
> 上游实读：`research/v3/species/B01-RP.md`（全 179 行，R1–R5 五条）· `research/v3/species/X01-L1C.md`（全 157 行）·
> `research/v3/species/GK3-SYMANZIK.md`（全 361 行，dim-6 ⇒ 正性丢失的唯一出处）· `research/v3/species/GK1-FREQSPLIT.md`（`m` 三义）·
> `research/gluing/REPORT.md:60–99,1485–1524`（`:79` 固定 Wilson、`:1504` R16）· `research/new-paths/SURVEY.md:355–394`（`:374`）。
> 外部（**全部只到摘要/书目级或二手正文片段，未逐字核正文**）：Lüscher–Weisz NPB **240** (1984) 349 · Lüscher CMP **54** (1977) 283 · Caracciolo–Pelissetto PLB **402** (1997) 335 · Osterwalder–Seiler 1978 · Menotti–Pelissetto CMP **113** (1987) 369 · Necco 博士论文 hep-lat/0306005 App. B · hep-lat/0208052。DOI/arXiv 见 §5.4。
> 标签：`[T]` 已证 · `[C]` 条件 · `[NG]` 具名 no-go / 不可作定理 · `[M]` 口径/记账 · `[E]` 可复跑算术 · `[R]` 外部（未逐字核正文）· `[O]` 精确缺失步。
> 纪律：每数带 `file:line`；未证不标 `[T]`；**无 QED**；**Clay = 0%**。

---

## §0 · 判词（先行）

| # | 问 | 答 | 等级 |
|---|---|---|---|
| **0.1** | 仓内/文献里存在定理「RP ⇒ `k₀ ≤ 2`」吗 | **不存在。** 仓内**0** 先例（§1.3）；文献里**没有**这一条（§2.1 检索式与命中，Q1–Q9）。这不是"没查到"，是"该命题的一般形式假"（0.2） | `[NG]` |
| **0.2** | 那这个命题本身成立吗 | **一般形式不成立。** Caracciolo–Pelissetto（PLB 402 (1997) 335, arXiv:hep-lat/9610015）在 2D σ 模型上**显式构造**了「`O(a²)` on-shell 改进 **且** 满足反射正性」的作用量 ⇒ `(RP) ∧ (k₀^{os}=4)` 可共存 ⇒ **RP 不强制 `k₀ ≤ 2`** | **`[R]`**（摘要逐字，正文未核） |
| **0.3** | 那么"RP 强制 `k₀ ≤ 2`"这条 no-go 拿不到，是不是整条思路作废 | **不是，但要换主语。** 真正成立的是**正锥刚性**：`c_L ≥ 0`（全部回路系数非负，= Osterwalder–Seiler/Lüscher 正性证明唯一适用的类）⇒ **RP 成立且 `k₀^{shell} = 2`**；反过来 `k₀^{shell} ≥ 4` ⇒ **必须**有负回路系数 ⇒ 掉出正锥 ⇒ RP 只能另证。⇒ **强制 `k₀ ≤ 2` 的是"正锥"，不是"RP"。** | **`[T]`+`[E]`** |
| **0.4** | 所以 `2/(m+2)` 是什么 | **它是"RP-可证类"的前沿，不是 RP 的前沿。** 逐字：`k ≤ k₀/(m+2)` 是 `[T]`+`[E]`（KLMN/T3，`GR1-SPECRES.md:27,257–259`）；在正锥里 `k₀ = 2` 是 `[T]`；**在正锥外 `k₀` 无界，而 RP 在正锥外既无定理也无反例（4D gauge）**。⇒ `2/(m+2)` 的"RP-保护"身份只能是 `[C]`：**条件 = "作用量留在正锥"**，而这是本项目唯一能证 RP 的类。 | **`[C]`**（不是 `[T]`） |
| **0.5** | 那 4D gauge 里能不能"特制"一个保住 RP 的改进（像 2D σ 那样） | **开问题（`[O]`）。** 无定理、无反例。三条结构性理由说明"2D 能"不能直接搬（§3.4）：① 2D σ 的换号通道是**约束恒等式** `φ·∂φ=0`，规范回路没有；② 4D gauge 里唯一同类通道是多迹/EOM-trivial 的 `tr(D_μF_{μν})²`，它**在壳上为零** ⇒ 只能抬 `k₀^{shell}` 不能抬 `k₀^{os}`，且是场的重定义（RP 是否保持未知）；③ 改进点在 (plaq,rect) 族里是**孤立点** `c₁=−1/12<0`，每多消一个 `p⁴` 不变式就要多一个负系数（N4：≥2 个）。 | `[O]` |
| **0.6** | `k₀` 是三义（本轮结构发现，仿 `GK1` 的 `m` 三义） | **`k₀^{shell}`（方向一致的相对小壳界，T3 乘子真正消费的那个）· `k₀^{axis}`（只沿轴/on-shell 横模，`GK3 §2.3` 的 `α₂(c₁)`）· `k₀^{os}`（Symanzik on-shell 幂）。**三者可分裂，且确实分裂**：树级 LW 作用量 `k₀^{axis}=4` 但 **`k₀^{shell}=2`**（N2，本轮新）。⇒ 讨论"`k₀=4`"必须点名读法，否则不可判。 | `[T]`+`[E]` |
| **0.7** | 本轮的实质增量 | 三件可复跑的：**(i)** 正锥刚性（**16** 个形状 × 3–6 方向共 ≈124 点的 `κ_L>0` + LP 不可行 + 多形状解 ≥2 负系数）；**(ii)** 树级 LW 的 `k₀^{shell}` **不**从 2 升到 4（`k₀` 三义分裂的定量见证）；**(iii)** 命题形式化与三读法判定表（§3.5）。**没有一条把 Clay 往前推。** | `[T]`+`[E]`/`[M]` |
| **0.8** | 对用户原假设的裁决 | 用户的逻辑链「dim-6 破坏 RP（真）＋ RP 强制 `k₀≤2`（**假**）⇒ `k ≤ 2/(m+2)` 是 RP 的推论」——**第二环断了**。可救的版本：「dim-6 破坏 RP（真）＋ **正锥**强制 `k₀ = 2`（真）⇒ **留在正锥 ⇒ `k ≤ 2/(m+2)`**；**离开正锥 ⇒ 必须先付出一份 4D gauge 的 RP 证明**，40 年无人付过」。⇒ 指数路线从"开放缺口"变成"**有名字的、定价的**条件"，但仍**不是** no-go 定理。 | — |

**一句话**：**`(RP) ⇒ k₀ ≤ 2` 不是定理，它的一般形式被 2D σ 的 CP 构造否掉；但"留在正锥"确实强制 `k₀ = 2`，而正锥是唯一能证 RP 的类 —— 所以 `2/(m+2)` 是"RP-可证类"的硬前沿，等级 `[C]`，不是 `[T]` 的 no-go。Clay 0%。**

---

## §1 · 前科核（`AGENTS §1` 四要素齐；判：仓内有没有人提过这个方向）

### 1.1 工具状态 / 实际检索式 / 命中数 / 失败标记

| # | 工具 | 状态 | 实际检索式（逐条） | 命中 |
|---|---|---|---|---|
| 1 | `glob` | **OK**（返回路径清单） | `research/v3/**/*.md`（733，截断显 100，全表落盘）；`**/B01-RP.md`；`**/gluing/REPORT.md`；`**/new-paths/SURVEY.md`；`**/*X01*` | 733 / 1 / 1 / 1 / 2 |
| 2 | `grep`（内置） | **OK**（返回行号 + 正文，非 provider failure） | ① `k_0\|k₀\|2/\(m\+2\)\|k_0/\(m\+2\)` @ `research/**/*.md` → **63 行 / 11 文件**（全表落盘）<br>② `RP.{0,12}(保护\|强制\|封顶\|禁令)\|k₀\s*≤\s*2\|k_0\s*≤\s*2\|Caracciolo\|Pelissetto\|perfect action\|Hasenfratz` @ 同路径 → **75 行**<br>③ `RP\|反射正性\|reflection\|正性\|k₀\|改进\|Symanzik\|Lüscher` @ `X01-L1C.md` → **1**（无关：`REPORT.md:110` 的"定理 A 用正性"）<br>④ 同上 @ `GK1-FREQSPLIT.md` → **1**（`GK1:168`，讲 `h_S`，与 RP 无关） | 63 / 75 / 1 / 1 |
| 3 | `read` | **OK**（全文返回） | `B01-RP.md`(179) · `X01-L1C.md`(157) · `GK3-SYMANZIK.md`(361) · `gluing/REPORT.md:60–99,1485–1524` · `new-paths/SURVEY.md:355–394` | — |
| 4 | `web_search` | **OK**（有效摘要 + 书目，无超时/空响应） | 9 条查询，见 §2.1 | 见 §2.1 |
| 5 | `web_fetch` | **OK**（HTTP 200） | `export.arxiv.org/abs/hep-lat/9610015v2`（CP 摘要全文）；`ar5iv.labs.arxiv.org/html/hep-lat/0306005`（Necco 论文目录 + App. B.2.1 标题 + 摘要） | 2 |
| 6 | `pwsh` | **OK**（`exit 0` 三次；一次 `exit 1` = **脚本自身断言**，斜边回路不是格点 link，**非环境失败**） | `py -3.12` + numpy 2.4.6 / scipy 1.18.1 | 3 次成功 / 1 次脚本 bug |
| 7 | `git commit` / MC / 采样 | **未执行 / 0 次 / 0 次** | — | — |

**`TOOL-FAIL`：0 次。空响应：0 次。`ENOENT` / `ripgrep provider failure`：0 次。**

### 1.2 前科清单（`file:line` 逐条）

| # | 前科 | 逐字给什么 | 与本轮的关系 | `file:line` | 级 |
|---|---|---|---|---|---|
| P1 | **R1–R5 五条**（RP 的定量内容） | R1 `T=e^{−aH}`、`0≤T≤1`、`H≥0`；R2 `spec T ⊂[0,1]`（砍负半轴）；R3 `ν≥0`、`C(t)` 完全单调；R4 `λ₁=1` 单重 + 正核；R5 强耦合下 `‖T|_{Ω^⊥}‖≤e^{−c(β)}` | **R1 是 `Δ_a` 链的唯一正性来源**（`GK3:229–234` 逐条对上）。本轮不新增 R；只问"R1–R5 是否隐含 `k₀≤2`"⇒ **不隐含**（R1–R4 全是符号/端点，R5 的定量输入是活动度不是 RP，`B01-RP:32–34`） | `B01-RP.md:24–34` | `[T]`（R1–R4）/`[C]`（R5） |
| P2 | **RP 对谱隙的盲性命题** | 「RP 约束的是权重的**符号**，`m` 是权重的**指数率**；两者是正权重空间上的独立坐标。RP 不蕴含 `m>0`」 | **本轮的命题是同一个"盲性"的第二层**：RP 连**作用量符号的阶**（`k₀`）也不约束 —— 因为 `k₀` 是**符号的 Taylor 阶**，不是测度的性质。⇒ 盲性从"符号 vs 大小"扩到"符号 vs 阶" | `B01-RP.md:97–105`（§3.4） | `[T]` |
| P3 | **固定 Wilson 的理由** | 「作用量取 Wilson…理由：Wilson 是**已知反射正性**的最简选择（Osterwalder–Seiler 1978）…改进作用量（Symanzik）虽有更好的 `a²` 行为，但**其 RP 需要单独证明**」 | **本轮把"需要单独证明"翻译成可判形式**：不是记账缺口，是**正锥外的整类**（§3.2 命题 R 的反向） | `gluing/REPORT.md:79` | `[T]` |
| P4 | **R16 风险行** | 「Wilson vs **改进作用量：某 session 用 Symanzik** … 低 ｜ **RP 未证** ｜ §1.2：固定 Wilson」 | 同 P3；说明这条已被当"风险"记过，未被当"定理"想过 | `gluing/REPORT.md:1504` | `[M]` |
| P5 | **换正则化的收益判** | 「**换正则化**（Symanzik/fixed-point/不同分块）——**不是清楚的正收益**…**RP 需重证**；项目固定 Wilson 正是为了 RP」 | 同 P3/P4；`fixed-point` 也在列 ⇒ FP 作用量的 RP 同样被记为"需重证" | `new-paths/SURVEY.md:374` | `[M]` |
| P6 | **GL2 接口**（RP vs 分块） | 「RP vs 分块的规范协变性…二者互不蕴含；`os3` 给出判据 `(B1) ⟺ A ∈ W ∩ R W R⁻¹`」 | 说明 RP 在仓内是**结构化**对象（反射面 ↔ 分块对齐），不是"一个不等式"；本轮的正锥也是结构化对象（系数符号） | `gluing/REPORT.md:1489` | `[T]` |
| P7 | **`m` 三义（方法论前例）** | 「`m` 是三义，任务书把其中两个并列了：`m_phys` / `m_KL` / …」 | **本轮 `k₀` 三义直接仿此**：`k₀^{shell}`/`k₀^{axis}`/`k₀^{os}`，且分裂是**可算的**（N2） | `GK1-FREQSPLIT.md:6,26,210,213` | `[T]`+`[M]` |
| P8 | **dim-6 ⇒ 正性丢失** | Lüscher–Weisz 1984 摘要逐字（见 §2.2）；⇒ `H_a≥0`（`B01-RP:26` R1）丢失，`Δ_a` 可能不良定 | 本轮**不重开**这条；它只提供"正锥外 RP 不可得"的**已观测**一侧 | `GK3-SYMANZIK.md:32,219–241` | `[R]` |
| P9 | **CP 已被登记为"特制改进存在"** | 「**不是"改进必然破坏 RP"的定理**。Lüscher–Weisz 的措辞是 "in general"：存在**特制**改进作用量保住 RP。**2D `N`-矢量模型**里就有显式构造」 | **本轮它就是反例本体**（§3.3）。GK3 只把它当"防过度声称"的脚注（`GK3:239`、`GK3:330` 第 10 条），**没把它当成"`(RP)⇒k₀≤2` 的反例"来判** —— 这是本轮与 GK3 的分界 | `GK3-SYMANZIK.md:239,330,343` | `[R]` |
| P10 | **`k ≤ 2/(m+2)` 的出处** | T3/KLMN：`\|𝔞^{(3)}_a(ψ)\| ≤ ε𝔮(ψ) + C_ε a²/ε·‖ψ‖²` ⇒ `ε=a^α` ⇒ 指数 `min(α, 2−(m+1)α)` ⇒ `k ≤ 2/(m+2) ≤ 1` | 本轮**不改**这个 `[T]`；只问它的 `k₀` 是不是"RP 强制"的 | `GR1-SPECRES.md:27,257–259`；`GP9-RB.md:120` | `[T]`+`[E]` |
| P11 | **`k₀=4` 的三道封顶** | ① KLMN `4/(m+2)≤2`；② T2 角点 `O(1)`（改进只减 22.7%）；③ 正性丢失 | 本轮**新增第四道**：④ 正锥刚性（`k₀^{shell}=4` 需 ≥2 个负回路系数）。④ 与 ①② 一样**与正性无关地**成立 | `GK3-SYMANZIK.md:33,243–268` | `[T]`+`[E]` |
| P12 | **`X01-L1C` 与本方向无关** | 全 157 行 = 端点簇分解簿记（`L1-c.1`、`\|ε_T\|≤Bū^T`、`δ=ln(1/ū)/a`）；**零 RP、零 `k₀`、零改进** | 任务书要求核它 —— 结论：**不是前科**（检索式 1.1 第 2③ 行，命中 1 条为无关引用） | `X01-L1C.md:1–157` | `[T]`（"无"） |
| P13 | **全仓 `2/(m+2)` 只落在 2 个文件** | `GR1`（定义 + C 层封顶）· `GT1B`（登记） | 说明"把 `2/(m+2)` 的出身改成 RP"这件事**从未被提出** | `GT1B-REDUCTION.md:44,103` | `[M]` |
| P14 | **位置图的四模型期望** | 自由场 `0` · 2D YM `>0` · M6 `0` · 单 plaquette 无定义 | 本轮 §3.6 用它做酸测；结论：本命题**进不了四模型筛**（它不给 `m` 一个数） | `GK3-SYMANZIK.md:15,299–306` | `[M]` |

### 1.3 判：有人提过吗 —— 三个问题分开答

1. **「`(RP) ⇒ k₀ ≤ 2` 这个方向」：0 命中，0 先例。** 全仓 `k₀` 只出现在 `GK3`（`GK3:28,33,184,249–252,266,322,330,357`）与 `GR1/GT1B` 的 KLMN 表里，全部把 `k₀` 当**给定输入**（Wilson ⇒ 2、all-orders 改进 ⇒ 4），**没有任何一处**问"`k₀` 的上界能不能由 RP 推出"。检索式见 §1.1 第 2①② 行；工具 OK；命中 63 行已逐行看，**0 行**涉及该内涵。
2. **「RP 是否限制改进」：有人提过，但停在"需重证/未证"（P3–P5 三处独立记载）**，没有一次走到"能否**证明**它被禁止"。
3. **「CP 1997 是反例」：`GK3` 已读到该摘要并逐字引用（P9），但只用作"防过度声称"，未用作"命题为假的见证"。** ⇒ 本轮不是发现新文献，是**把已登记的脚注升级为判定**（并补上正锥定理，使"能救的那一半"有 `[T]`）。

---

## §2 · 外部检索结果（必须给 DOI/arXiv，不许编造）

### 2.1 检索执行（`AGENTS §1` 四要素）

| # | 查询（逐条） | 工具状态 | 命中/结果 |
|---|---|---|---|
| Q1 | `Lüscher Weisz "Definition and general properties of the transfer matrix" "continuum limit improved" Nuclear Physics B 240 1984 349` | **OK**，有效摘要 | ScienceDirect `pii/0550321384902700`（摘要逐字）+ Scopus `0000806223`（240(3) 349–361）+ MathSciNet MR0762195（卷期页一致） |
| Q2 | `Lüscher "Construction of a selfadjoint, strictly positive transfer matrix" Comm Math Phys 54 1977 283` | **OK** | CMP 54(3) 283–292（1977-10）；MathSciNet MR1552616；**作者 = M. Lüscher，不是 K. Wilson**（勘误见 §2.5） |
| Q3 | `Caracciolo Pelissetto "Comparing different improvement programs for the N-vector model" reflection positivity 1997` | **OK** | ADS `1997PhLB..402..335C` + arXiv abs（**v2 摘要全文已 fetch**）+ DOI `10.1016/S0370-2693(97)00495-4` |
| Q4 | `Menotti Pelissetto positivity transfer matrix improved lattice gauge action reflection positivity class of actions` | **OK**，但**未见**"改进作用量的 RP 定理" | 只命中 Menotti–Pelissetto 1987/1988（**Wilson 作用量**的 site-plane 推广）；zbMATH 01573464 摘要逐字 "Osterwalder–Seiler's proof … is extended to the case in which the reflection is performed with respect to planes containing sites" ⇒ **仍限 Wilson** |
| Q5 | `Hasenfratz Niedermayer perfect action reflection positivity 4D SU(3) positivity transfer matrix` | **OK**，**未见**"FP 作用量的 RP 定理" | 命中 HN 1994（`hep-lat/9308004`，方法示例在 2D σ）、DeGrand–Hasenfratz–Hasenfratz–Niedermayer SU(3) FP（`hep-lat/9506030`，"instanton 解尺度不变、谱到 1-loop 精确、无 `a^n`/`g²a^n` cutoff 效应"）；**RP 一项：0 命中** |
| Q6 | `improved gauge action reflection positivity violation Iwasaki DBW2 transfer matrix negative eigenvalues` | **OK** | `hep-lat/0306005`（Necco, DESY-THESIS-2003-016）App. B.2.1「**Unphysical poles in the propagator**」；正文逐字 "for the Iwasaki and DBW2 action, due to **the violation of physical positivity** described in appendix B, one expects also negative contributions"；`hep-lat/0208052`（Iwasaki/DBW2 的 `c_j` = `0.331` / `1.4088`，**符号约定见 §2.3 第 1 行注**）+ LW 1984 的 "complex eigenvalues may occur" 转述；`0710.0425`（DBW2 下 domain-wall 转移矩阵出现**负本征模** ⇒ 振荡） |
| Q7 | `"fixed point action" SU(3) gauge loop coefficients negative signs positivity problem` | **OK**，**未见** RP 结论 | `1004.3462` 综述片段："The proposal of Hasenfratz and Niedermayer (1994) is to simulate the action `βA^{FP}` for β large"；**无 RP 判定** |
| Q8 | `arXiv 1004.3462 improved actions review Hasenfratz Niedermayer perfect action positivity` | **OK** | 同上；仅书目 |
| Q9 | 上一批的补充（Symanzik σ 模型 nnn 负耦合 / `beta_2` 与 RP） | **OK**，**0 命中** | 无一条给出"2D σ 的 nnn 耦合必须 `≥0` 才能 RP"的**可引定理**；故本轮**不主张**该句，只引 CP 的构造 |

**`TOOL-FAIL`：0。空响应：0。超时：0。**

### 2.2 决定性的两条（逐字，`[R]`）

**X1 — Lüscher–Weisz, NPB 240 (1984) 349–361**, DOI `10.1016/0550-3213(84)90270-0`（摘要逐字，两处独立索引一致：ScienceDirect + MathSciNet MR0762195）：
> "When operators of **dimension 6** are added to the standard Wilson action in lattice gauge theories, **physical positivity is lost in general**. We show that a transfer matrix can nevertheless be defined. Its properties are, however, unusual: **complex eigenvalues may occur** …, and **there are always contributions in the spectral decomposition of two-point functions that come with a negative weight.**"

**X2 — Caracciolo–Pelissetto, PLB 402 (1997) 335–340**, arXiv:`hep-lat/9610015`（v2 摘要全文已 fetch，逐字）：
> "We discuss the connection between various types of improved actions in the context of the **two-dimensional sigma-model**. We also discuss spectrum-improved actions showing that these actions **do not have any improved behaviour**. **An `O(a²)` on-shell improved action with all couplings defined on a plaquette and satisfying reflection positivity is also explicitly constructed.**"

**X2 的三层读法（必须分开，不许合并）**：
1. `k₀^{os}=4` 与 **RP** 在该构造上**并存** ⇒ **`(P-os)` 假**，`[R]`。
2. 它**不**直接说 `k₀^{shell}=4`（摘要只说 on-shell 改进；"all couplings defined on a plaquette" 是其 RP 的来源）⇒ **`(P-shell)` 在 2D σ 上是否假，未定**（正文未核）。
3. 它**不**涉及规范回路基 ⇒ 与 §3.2 的正锥命题**不冲突**（§3.3）。

### 2.3 已观测的 4D gauge 一侧（不是"未证"，是"违例被记录"）

| # | 事实 | 来源 | 级 |
|---|---|---|---|
| 1 | 已知 4D 改进规范作用量的 `1×2` 系数：LW `c₁=−1/12`（`GK3` N6，本轮 N2 独立复现）；Iwasaki `0.331`、DBW2 `1.4088`（`hep-lat/0208052` 表的 `c_j` 列，**该表的符号约定本文件未逐字核**；按标准约定即 `c₁=−0.331`、`−1.4088`）⇒ **全部落在正锥外** | `hep-lat/0208052`（片段） | `[R]` |
| 2 | 负系数在自由传播子上就是**非物理极点**（T 矩阵有复/负模） | `hep-lat/0306005` App. B.2.1「Unphysical poles in the propagator」（目录逐字） | `[R]` |
| 3 | "due to **the violation of physical positivity** … one expects also negative contributions"（Iwasaki/DBW2） | `hep-lat/0306005` 正文片段 | `[R]`（片段） |
| 4 | DBW2 下 domain-wall 转移矩阵出现**负本征模** ⇒ 关联函数振荡 | `0710.0425` 摘要片段 | `[R]` |
| 5 | Wilson 的 RP 有定理；**改进作用量的 RP/正转移矩阵：无证明** | Lüscher 1977（`[R]`，摘要逐字见下）+ Osterwalder–Seiler 1978 + Menotti–Pelissetto 1987（仓内 `os3-reflection-positivity/REPORT.md:13,148,624`；`R1-BRIDGE.md:53`） | `[R]`/`[T-标准]`（Wilson 侧） |

**X3 — Lüscher, CMP 54 (1977) 283–292**, DOI `10.1007/BF01614090`（摘要逐字，Springer + Project Euclid `cmp/1103900872`）：
> "It is shown that **physical positivity holds in Wilson's lattice gauge theories** …"

⇒ **RP 属于 Wilson（正锥）；`k₀=4` 属于正锥外。两句话在文献里是同一个二分，只是从没被写成一条定理。**

### 2.4 Hasenfratz–Niedermayer（完美作用量）与 RP：**检索到底，无结论**

- FP 作用量的**定义性质**：classically perfect（自由谱精确）+ 1-loop perfect ⇒ `k₀^{os} = ∞`（"no `a^n` nor `g²a^n` cut-off effects"，`hep-lat/9506030` 摘要片段 `[R]`）。
- **RP：0 命中。** 检索式 Q5/Q7（§2.1），工具 OK，命中 0 条给出"FP 作用量满足/违反 RP"的陈述。仓内同样：`FORGE5-FPACTION.md:19` 只谈"positivity = BPS positivity"，`B03-BLOCK.md:59` 只谈构造，**都不谈 RP**。
- **但本轮有一条结构性推论（可算，不需要文献）**：若 FP 的**自由符号**是精确的（`h(p) ≡ |p|²`，这正是 classically perfect 的含义），则 `κ = 0`；而正锥里 `κ = Σ c_L κ_L > 0`（`c_L ≥ 0`，N1）⇒ **FP 作用量必不在正锥内**（必有负回路系数，或系数和不可绝对收敛）。⇒ `gluing/REPORT.md:79` 那句"其 RP 需要单独证明"**不是记账，是结构必然**。`[T]`（模 N1 的 16 形状扫描；一般 `L` 的 `κ_L>0` 为 `[O]`）。

### 2.5 勘误（只记账，不改原文件）

| # | 位置 | 原文 | 事实 | 影响 | 级 |
|---|---|---|---|---|---|
| **E-1** | 任务书 §2 第 2 条 | 「**Wilson 1977**（CMP 54 283）：RP 对 Wilson 作用量的证明」 | **CMP 54 (1977) 283–292 的作者是 M. Lüscher**，标题 *Construction of a selfadjoint, strictly positive transfer matrix for Euclidean lattice gauge theories*；它证的是 **Wilson 的格点规范理论**（摘要逐字 "physical positivity holds in **Wilson's** lattice gauge theories"）。**K. Wilson 1977 没有 CMP 54 283 这篇。** | 引用归属错；**结论不变**（RP 属 Wilson 作用量） | `[T]`（书目，三处独立索引） |
| **E-2** | 任务书 §2 第 4 条 | 「Caracciolo–Pelissetto … 2D σ 模型特制改进保住 RP ⇒ 这是反例」 | **对，但只对 `k₀^{os}` 成立**：摘要只说 "`O(a²)` **on-shell** improved"。说它否掉 `(RP)⇒k₀^{shell}≤2` 需要"其自由符号的 `p⁴` 系数也为零"——**摘要不含此信息，正文未核**。 | 反例的作用域必须写清（§3.3 第 2 层） | `[R]`/`[M]` |
| **E-3** | `GK3-SYMANZIK.md:184,249–252` | 隐含「`k₀=4` ⇒ KLMN 乘子从 `a²` 换 `a⁴`」 | **在 `k₀^{shell}` 读法下，连树级 LW 都不成立**：本轮 N2 给出 `c₁=−1/12` 时轴向 `p⁴` 系数为 0，但 2-平面方向 `+0.0278`、单向对角 `+0.0833` ⇒ 方向一致的相对界仍是 `O(a²Λ²)` ⇒ `k₀^{shell}=2`。要 `k₀^{shell}=4` 必须**同时**消掉 `Σp⁴` 与 `(Σp²)²` 两个不变式（N4：需 ≥2 个负系数，量级 `−4.5, −71.5`）。 | `GK3` 的"`k:2→4`"**只在 `k₀^{axis}`/`k₀^{os}` 读法下成立**；不该外推到 KLMN 乘子 | `[T]`+`[E]` |

---

## §3 · 命题与判定

### 3.1 命题的确切形式（三个版本，缺一不可）

**类 `𝔄`（唯一被考虑的类）**：4D 超立方格 `Λ_{L,a}=(ℤ/Nℤ)^4`，紧单群 `G`，作用量
```
S = Σ_L c_L · Σ_{γ ∈ shape L} (1 − (1/N_c) Re tr U_γ),      c_L ∈ ℝ
```
- **形状集**：闭合格点回路（单迹），对超立方群与时间反射**封闭**；
- **局域性：允许非局域。** 只要求 `Σ_L |c_L| < ∞`（系数可求和）⇒ FP 类作用量**在被考虑的范围内**；`1 − Re tr U_γ` 逐项定义。
- **允许负系数：是**（这正是全部内容所在）。
- **测度**：`μ_{L,a} = Z^{-1} e^{−S} ∏ dU`，Haar；`Θ` = 时间反射（过链平面 OS'78 / 过格点 Menotti–Pelissetto'87；本轮的回路线性化与反射面选择无关）。

**`k₀` 三义（本轮结构发现；仿 `GK1-FREQSPLIT.md:26` 的 `m` 三义）**

| 读法 | 定义（`h` = 自由横模符号，`Q_{μν}(p) = Σ_L c_L w_{Lμ}w^*_{Lν}` 的横本征值） | 谁消费 | Wilson | 树级 LW (`c₁=−1/12`) |
|---|---|---|---|---|
| **`k₀^{shell}`** | 最大 `k`：`\|h(p) − \|p\|²\| ≤ C(aΛ)^k \|p\|²`，`∀\|ap\|≤Λ`，**方向一致** | **`GR1:257` 的 KLMN 乘子 `C_ε a^{k₀}/ε`** ⇒ `k ≤ k₀/(m+2)` | **2** | **2**（N2，本轮新） |
| **`k₀^{axis}`** | 同上但只沿轴 / on-shell 横模（`p=(iE,p₁,0,0)`） | `GK3 §2.3` 的 `α₂(c₁)=−(1/12+c₁)`、"`k:2→4`" | **2** | **4** |
| **`k₀^{os}`** | on-shell 物理量（质量/色散/顶点）的 `a` 幂 | Symanzik 纲领；CP 1997 的 "on-shell improved" | 2 | 4 |

**三个命题**
- **(P-shell)** `μ_{L,a}` 反射正性 ⇒ `k₀^{shell} ≤ 2`；
- **(P-axis)** 同上 ⇒ `k₀^{axis} ≤ 2`；
- **(P-os)** 同上 ⇒ `k₀^{os} ≤ 2`。

> **为什么要三个**：T3 乘子消费的是 `k₀^{shell}`（**唯一**能改 `2/(m+2)` 的那个）；GK3 算的是 `k₀^{axis}`；CP 造的是 `k₀^{os}`。**把三者当一个数讨论，就是本轮要修的那个错**（同 `GK1` 修的 `m` 并列）。

### 3.2 正锥定理（本轮的主结果，`𝒞₊` 内成立）

**定义（正锥）** `𝒞₊ := { c : c_L ≥ 0 ∀L }`。

> **命题 R.** 设 `S` 由 `c ∈ 𝒞₊` 给出（形状集反射封闭、`Σc_L<∞`）。则
> **(R-a)** `μ` 反射正性成立；
> **(R-b)** `k₀^{shell}[S] = 2` 且 `k₀^{axis}[S] = 2`（**等号**，不是 `≤`）；
> **(R-c)** 反之 `k₀^{shell}[S] ≥ 4` ⇒ `c ∉ 𝒞₊`（至少一个 `c_L < 0`；定量版见 N4：要消 4 个方向的 `p⁴` 需 **≥2** 个负系数）。

**证明骨架（4 步，每步可核）**

1. **线性化 ⇒ 秩一 PSD 之和。** 对闭合回路 `γ`，`1 − (1/N)Re tr U_γ = (g²/4)|∮_γ A|² + O(A³)`，其中
   `∮_γ A = Σ_{k} s_k a A_{μ_k}(p) e^{−ip·x_k} ≡ a w_γ(p)·A(p)`。
   ⇒ `Q_{μν}(p) = (g²/4) Σ_L c_L w_{Lμ}w^*_{Lν}`，每个 `w_L w_L^†` 秩一 PSD。`[T]`
2. **闭路 ⇒ `w_L(0)=0`。** `Σ_k s_k e_{μ_k} = 0`（向量闭合）⇒ `w` 在 `p=0` 无常数项 ⇒ `h_L := |w_{L,T}|² = A_L p² + B_L p⁴ + O(p⁶)`，`A_L > 0`（回路张成时空）。`[T]`
3. **符号刚性（N1，`[E]`）。** **16** 个形状 × 3–6 个方向（共 ≈124 个 (形状,方向) 点），全部 `B_L(u) < 0`（等价 `κ_L := −B_L/A_L > 0`，实测区间 `[0.0833, 3.2745]`，最小 = 单 plaquette 的 `1/12`）。含平面矩形 `1×1…4×4`、`5×1`、`6×2`、L 形、zigzag、3D corner（chair）、8-link twist。⇒ **每一个正回路系数都往 `p⁴` 上添一个负项**。`[E]`（一般 `L` 未证，`[O]`）
4. **锥上求和。** 按 `Σ c_L A_L = 1` 归一后 `h = Σ c_L h_L = |p|² − κ(p̂)|p|⁴ + O(a⁶|p|⁴)`，`κ(p̂) = Σ c_L κ_L(p̂) > 0`（正项和，`c≠0`）⇒ `k₀^{shell} = 2`。反向：`k₀^{shell} ≥ 4` 要 `Σ c_L B_L(u) = 0 ∀u`；由第 3 步 `B_L(u)<0` ⇒ 必须 `c ∉ 𝒞₊`。`[T]`（模 N1）

**(R-a) 的依据**：`e^{−βc_L(1−(1/N)Re tr U_L)} = e^{−βc_L}·e^{(βc_L/N)Re tr U_L}`；`c_L>0` 时 `e^{x Re tr U}` 在特征标基下的系数非负（SU(2) 显式：`2I_{2r+1}(x)/x > 0`；一般紧群：由 `e^{x(χ_f+χ̄_f)/2}` 展开 + Littlewood–Richardson 系数非负），乘积的正系数由 LR 保持 ⇒ Osterwalder–Seiler/Lüscher 的正核论证逐项适用。`[T-标准]`（仓内依据：`os3-reflection-positivity/REPORT.md:13`；`R1-BRIDGE.md:53`；`gluing/REPORT.md:79` 的"Osterwalder–Seiler 1978"）。

**数值证人（`[E]`，确定性，脚本 §5.2）**

| # | 计算 | 结果 |
|---|---|---|
| **N1** | 单形状 `κ_L = −B_L/A_L`，**16** 形状 × 3–6 方向（≈124 点） | 全部 `> 0`；最小 `+0.08333`（plaq），最大 `+3.2745`（`6×2`，单向对角）；**0 个反例** |
| **N2** | (plaq,rect) 族 `c₀+8c₁=1`，轴向 `B/A` 与 `c₁` | `A ≡ 3.00000000`（`c₁`-无关 ✓ 归一化）；`B/A = −(1/12+c₁)`，逐点差 `≤5e−7` ⇒ **复现 `GK3` N6**；唯一零点 `c₁=−1/12`；**但在 2-平面方向 `c₁=−1/12` 时 `B/A=+0.0278`、单向对角 `+0.0833` ≠ 0 ⇒ `k₀^{shell}=2`** |
| **N3** | LP：`c ≥ 0, Σc_LA_L=3, Σc_LB_L(u)=0` | **不可行**（HiGHS status 2），**每个方向单独也不可行** |
| **N4** | 允许负系数、5 形状解 4 方向条件 | 6 个形状组合的解**全部带 ≥2 个负系数**：如 `[−71.52, +27.87, −4.50, +4.00, −1.11]` |

⇒ **`𝒞₊` 里 "RP" 与 "`k₀ ≥ 4`" 互斥；`2/(m+2)` 是正锥的硬前沿。**

### 3.3 一般形式不成立：CP 1997 是反例（且它不推翻 §3.2）

**判：`(P-os)` 假，`[R]`。** CP 的构造同时给出 "`O(a²)` on-shell improved" 与 "satisfying reflection positivity"（§2.2 X2 逐字）⇒ **RP 不强制 `k₀^{os} ≤ 2`**。

**三层作用域（必须分开）：**

| 层 | 命题 | CP 是否否掉 | 依据 |
|---|---|---|---|
| L1 | `(RP) ⇒ k₀^{os} ≤ 2` | **否掉** | 摘要逐字 |
| L2 | `(RP) ⇒ k₀^{shell} ≤ 2` | **未定**（摘要无"自由符号 `p⁴` 为零"的信息；正文未核） | `[O]` |
| L3 | `c ∈ 𝒞₊ ⇒ k₀^{shell}=2` | **不涉及**（CP 的动作不是规范回路基；见 §3.4 理由 1） | §3.2 |

**为什么 L1 被否掉却不推翻 L3 —— 机制（不是借口）**：2D σ 的二次型可以由**约束** `φ²=1` 产生的恒等式重写：`(∂φ)² = −φ·∂²φ`，于是 `O(a²)` 修正的负 `p⁴` 项可以被换成**正**的 `p⁴` 项（`(φ·∂²φ)²` 型），而**正性不破**（它是平方）。规范回路没有这条恒等式：`|∮_γ A|²` 的 `p⁴` 系数符号由第 3 步（N1）刚性锁死。**这不是维数差异，是"有约束 vs 无约束"的差异**（下面 §3.4 第 1 条）。

### 3.4 4D gauge 里能否"特制"改进？—— 判定：**开问题**，附三条结构性理由

**判：无定理、无反例 ⇒ `[O]`。** 换句话说，**本轮的 no-go 只在 `𝒞₊` 内成立**；"4D gauge 存在 RP-保住的 `k₀=4` 作用量"这句话，既没被证，也没被否。

三条理由（解释"2D 能、4D gauge 未必能"）：

| # | 机制 | 2D σ（CP 能） | 4D gauge（？） | 级 |
|---|---|---|---|---|
| 1 | **约束恒等式** | `φ·∂_μφ = 0` ⇒ 负 `p⁴` 可换成正的平方项，**在正锥内** | **无此恒等式**；`|∮A|²` 的 `p⁴` 符号刚性（N1） | `[T]`+`[E]` |
| 2 | **多迹 / EOM-trivial 算符** | 约束把多迹项化为单迹 | 唯一同类：`tr(D_μF_{μν})²`（平方 ⇒ 二次型 `(∂F)² ~ p⁴` 为正，可抵负 `p⁴`），但**在壳上为零** ⇒ 只抬 `k₀^{shell}`、不抬 `k₀^{os}`；且它是场的重定义（Jacobian = 一个局域 dim-6 插入）⇒ **RP 是否保持未知**。**本轮未逐项算它对两个 `p⁴` 不变式的贡献** | `[M]`（未算）+ `[O]`（RP） |
| 3 | **改进点的余维** | 改进条件 codim 1，耦合空间 ≥2 维 ⇒ 留一条**线**，RP 是开条件 ⇒ 有可能 | (plaq,rect) 族里改进点是**孤立点** `c₁=−1/12<0`（`GK3` N6 + 本轮 N2）；每多消一个 `p⁴` 不变式就多一个负系数（N4：≥2 个，量级 `−71.5`）⇒ 正锥外解是**高余维**的 | `[T]`+`[E]` |

**必须同时记住的"已观测"一侧**（比"未证"更硬）：LW/Iwasaki/DBW2 的 `c₁<0`，而自由传播子上就出现**非物理极点**、DBW2 出现**负本征模**（§2.3 第 2–4 行）；LW 1984 的 "in general" 明确留了口子。⇒ **既不能写成"已证不可能"，也不能写成"只是没证"。**

### 3.5 判定表（三读法 × 三等级）

| # | 命题 | 一般形式 | 限 `𝒞₊` | 等级 | 依据 |
|---|---|---|---|---|---|
| 1 | `(RP) ⇒ k₀^{shell} ≤ 2` | **不可作定理**（无反例、无定理；正锥外无结论） | **成立，且 `=2`** | `[NG]` / `[T]`+`[E]` | §3.2、§3.3 L2 |
| 2 | `(RP) ⇒ k₀^{axis} ≤ 2` | 同上 | **成立，且 `=2`** | `[NG]` / `[T]`+`[E]` | 同上 |
| 3 | `(RP) ⇒ k₀^{os} ≤ 2` | **假**（CP 1997） | 成立（`[T]`） | **`[R]`** | §2.2 X2 |
| 4 | `c ∈ 𝒞₊ ⇒ k₀^{shell}=2`（正锥刚性） | **成立（等号）** | — | **`[T]`+`[E]`** | 命题 R |
| 5 | 4D gauge `𝒞₊` 外的 RP 保住改进 | **开** | — | `[O]` | §3.4 |
| 6 | `2/(m+2)` 是"RP-保护前沿"吗 | **不是**（RP 不强制）；**是"RP-可证类"前沿**（`𝒞₊` 内强制） | — | **`[C]`** | 汇总 1–5 |
| 7 | 用户的 no-go 链是否成立 | **不成立**（第二环 `RP⇒k₀≤2` 缺）；可救版本见 §0.8 | — | `[NG]` | — |

### 3.6 酸测（自由场 / 2D（CP 构造）/ M6 / 单 plaquette）

| # | 模型 | 期望（`GK3:299–306` 转引 `position-current.html:395`） | 本轮结果 | 判 |
|---|---|---|---|---|
| 1 | **自由场** | `0` | 本命题**不给 `m` 一个数**：`k₀` 只依赖符号 `h`，与 `m`、`L` 无关（`[T]`）⇒ 无法在自由场上"取到 0 还是非 0" | **N/A（`[M]`）** |
| 2 | **2D（CP 构造）** | 2D 必须不能凭空制造 `a`-依赖；2D 精确零点要保住 | **一般形式 `(P-os)` 在此被一票否决**（CP：RP + on-shell 改进并存）⇒ **命题的一般形式 FAIL**；`𝒞₊` 形式 PASS（2D σ 不在规范回路基内） | **FAIL（一般）/ PASS（`𝒞₊`）** |
| 3 | **M6** | `0` | 无 `a`、无 `β`、无格点回路 ⇒ `k₀` **不可陈述**（同 `GK3 §5` 第 3 行） | **N/A（`[M]`）** |
| 4 | **单 plaquette** | 无定义 | 无 transfer 方向 ⇒ 回路基退化 ⇒ `k₀` 无定义 | **N/A（`[M]`）** |

> **酸测的诊断价值（必须写清）**：本命题**进不了四模型筛**——它不是"给 `m` 一个数"的工具，而是**关于"哪些作用量类里 RP 与 `k₀` 可共存"的元命题**。⇒ 它的正确用法是**改写路线状态**（§0.8），不是当新工具去跑试金石。**任何把它当"新物种工具"申报的做法都该被驳**。

---

## §4 · 对 Clay 的影响（一句话 + 百分比）

**一句话**：本命题**不改变任何 Clay 量**——`k` 是定量 `SPEC-RES` 的 C 层（`CERT-C1` 的 `θ_a ≤ Ca^k`）的技术输入（`GR1:336,350`），不是隙本身；即使 `k ≤ 2/(m+2)` 升成 `[T]` 的 no-go，也既不证 `Δ_cont>0`、也不证 `(NEW)`。**唯一可说的净效应**：把"指数路线为什么卡在 `2/(m+2)`"从"工具缺"改写成"**作用量类的代价**"，并给出天花板 —— 即便有人把 4D gauge 的 RP 改进造出来（`k₀=4`），KLMN 硬顶也只给 `4/(m+2) ≤ 2`（`m=0` 取到 2，`GK3:28,252`），即**最好情况是把 `k` 从 `≤1` 抬到 `≤2`，永不到 4**。

**Clay = 0%。**（本轮：无 QED；未证之物未标 `[T]`；未 commit；无 MC。）

---

## §5 · 台账

### 5.1 本轮改动清单

| # | 内容 | 级 |
|---|---|---|
| 1 | 把「`(RP) ⇒ k₀ ≤ 2`」形式化为三读法（shell/axis/os）× 三命题 | `[T]`（定义） |
| 2 | 正锥刚性 `c_L ≥ 0 ⇒ k₀^{shell}=k₀^{axis}=2`，反向 `k₀^{shell}≥4 ⇒ c∉𝒞₊` | `[T]`+`[E]` |
| 3 | **16** 形状 `κ_L>0` 扫描（≈124 点）；LP 不可行；4 方向 shell 改进需 ≥2 负系数 | `[E]` |
| 4 | `k₀` 三义分裂的定量见证：树级 LW 的 `k₀^{axis}=4` 但 `k₀^{shell}=2` | `[T]`+`[E]` |
| 5 | 一般形式被 CP 1997 否掉（`(P-os)` 假）；作用域三层 | `[R]` |
| 6 | 4D gauge 特制改进 = 开问题 + 三条结构性理由 | `[O]` |
| 7 | 三处勘误 E-1（CMP 54 283 作者）/E-2（CP 作用域）/E-3（`GK3` 的 `k₀` 读法外推） | `[T]`/`[R]`/`[M]` |

### 5.2 数值与脚本（确定性，零随机，零 MC）

| # | 脚本（`$env:TEMP`，**未写入仓库**） | 内容 | 状态 |
|---|---|---|---|
| 1 | `gl2_rporder.py` | 8 形状 `κ_L`；(plaq,rect) 族 `B/A(c₁)` 三方向；根定位 | `exit 0` |
| 2 | `gl2_shapes.py` | 形状族 `κ_L` 扫描（矩形族 1×1…4×4、5×1、6×2、L、zigzag、chair、twist） | `exit 1`（**脚本自身断言**：`parallelogram` 含斜边、不是格点 link，中断于此）。**已定义 22 个形状，中断前 16 行结果已打印**（含全部矩形族与两个非平面回路），**非环境失败** |
| 3 | `gl2_lp.py` | 6 形状 × 6 方向 `B_L`；LP 可行性（全方向 + 逐方向） | `exit 0`，HiGHS status 2（不可行） |
| 4 | `gl2_shell.py` | 5/6 形状解 `A·c=3, B·c=0`（4 方向），报负系数个数 | `exit 0` |

复现：`py -3.12 %TEMP%\gl2_rporder.py`（numpy 2.4.6 / scipy 1.18.1）。

### 5.3 诚实边界（必读）

1. **`k₀^{shell}` 的一般 `L` 命题（"每个闭合回路 `κ_L>0`"）本轮只做了 16 个形状 × 3–6 方向的扫描，不是定理。** 它的失败方式很明确：某个回路若 `κ_L<0`（`|w|²=(Σs cos)²+(Σs sin)²` 的 `O(p⁴)` 项为 `−⅓XY+¼Z²`，符号不先验固定），正锥刚性就有反例。**这是本文件最该被攻的一格**，`[O]`。
2. **`𝒞₊ ⊆ {RP}` 是包含，不是相等。** 负系数作用量**可能**仍然 RP（Lüscher–Weisz 只说 "in general"）。所以"正锥外无 RP"**不是**本文件的断言；本文件只断言"正锥外**没有** RP 证明"。
3. **`(R-a)`（`c_L≥0 ⇒ RP`）依据 OS/Lüscher 的标准论证**，本轮**未重证**（仓内依据 `os3` §2.1、`R1-BRIDGE:53`）；`[T-标准]`，不是本轮新证。
4. **CP 1997 只读到摘要（v2 全文摘要已 fetch），正文未核。** 因此 `(P-shell)` 在 2D σ 上是否也假，**未定**（§3.3 L2）。
5. **FP 作用量的"`κ=0` ⇒ 负系数"是推论不是文献结论**（§2.4）；若 FP 的自由符号并非逐点等于 `|p|²`（例如只在 on-shell 上精确），该推论不成立 —— 这正是"classically perfect"的精确含义我**没有**逐字核的地方。
6. **`E-3` 是对 `GK3` 的读法收窄，不是推翻 `GK3`**：`GK3` 的 `k=4` 判 `[NG]` 依赖三道封顶，其中 ①② 与 `k₀` 读法无关；本轮只指出"`k₀=4` ⇒ 乘子 `a⁴`"这一步在 `k₀^{shell}` 读法下**连树级 LW 都不给**。未改 `GK3`。
7. **未做**：任何 4D `SU(N)` 的 `β=2.5` 数值（= 墙）；FP 作用量系数符号的逐项核；2D σ 的 `p⁴` 符号（若做，就能判 §3.3 L2）。
8. **不改 `docs/position-current.html`**：`AGENTS §6` 的触发条件 2（"某条路线判死或复活"）**在本轮触发**——「指数路线因 `k₀` 封顶而不可抬」由"开放缺口"变为"**正锥外**才有机会、且正锥外无 RP 证明"。按任务约束无 `docs/` 写权，**留给 Lead**。建议死路台账新增一条：**GL2-RPORDER = `[C]`（正锥前沿 `k₀=2`，`2/(m+2)`）+ `[O]`（4D gauge 正锥外）**，不要记成 `[NG]`。
9. **No QED；Clay = 0%。**
10. **不许说的话**：①"已证 RP ⇒ `k₀ ≤ 2`"（假）；②"`2/(m+2)` 已是 `[T]` 级 no-go"（只是 `[C]`）；③"4D gauge 不可能有 RP 改进"（开问题）；④"CP 是 `(P-shell)` 的反例"（只到 `(P-os)`）；⑤"本轮造了一个新工具"（它是元命题，进不了四模型筛）。

### 5.4 外部引用（DOI / arXiv，全部**摘要或书目级**，未逐字核正文；标注读入级别）

| # | 引用 | 标识 | 读到哪 | 本轮用途 |
|---|---|---|---|---|
| X1 | **M. Lüscher, P. Weisz, *Definition and general properties of the transfer matrix in continuum limit improved lattice gauge theories*, NPB 240 (1984) 349–361** | DOI `10.1016/0550-3213(84)90270-0`；MathSciNet MR0762195 | **摘要逐字**（ScienceDirect + MathSciNet 卷期页互证） | dim-6 ⇒ physical positivity lost / 复本征值 / 负权（§2.2） |
| X2 | **S. Caracciolo, A. Pelissetto, *Comparing different improvement programs for the N-vector model*, PLB 402 (1997) 335–340** | DOI `10.1016/S0370-2693(97)00495-4`；arXiv:`hep-lat/9610015` | **v2 摘要全文 fetch**（`export.arxiv.org/abs/hep-lat/9610015v2`） | **反例本体**（§3.3） |
| X3 | M. Lüscher, *Construction of a selfadjoint, strictly positive transfer matrix for Euclidean lattice gauge theories*, CMP 54 (1977) 283–292 | DOI `10.1007/BF01614090`；MR1552616 | 摘要逐字（"physical positivity holds in **Wilson's** lattice gauge theories"） | RP 属 Wilson/正锥；勘误 E-1 |
| X4 | K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. 110 (1978) 440–471 | DOI `10.1016/0003-4916(78)90039-8` | **仅仓内转引**（`gluing/REPORT.md:79`；`B01-RP:149` 标注"摘要级"） | Wilson RP 的标准引用（不独立承重） |
| X5 | P. Menotti, A. Pelissetto, *General proof of Osterwalder–Schrader positivity for the Wilson action*, CMP 113 (1987) 369–373 | DOI `10.1007/BF01221251`；zbMATH 01573464 | zbMATH 摘要片段（"extended to … planes containing sites"） | 过格点平面 RP，**仍限 Wilson** |
| X6 | S. Necco, *The static quark potential and scaling behavior of SU(3) lattice Yang-Mills theory*, DESY-THESIS-2003-016 | arXiv:`hep-lat/0306005` | 目录 + 摘要 + 正文片段（ar5iv） | 4D gauge 改进作用量的**非物理极点 / 正性违例**（§2.3） |
| X7 | *Universality and RG-improved gauge actions* | arXiv:`hep-lat/0208052` | 片段（`c_j` 表：Iwasaki 0.331、DBW2 1.4088；LW 1984 转述）；**符号约定未逐字核**（§2.3 第 1 行） | 已知改进作用量的 `1×2` 系数在正锥外 |
| X8 | P. Hasenfratz, F. Niedermayer, *Perfect lattice action for asymptotically free theories*, NPB 414 (1994) 785 | arXiv:`hep-lat/9308004` | 摘要级 / 二手片段 | FP 纲领（**RP：0 命中**） |
| X9 | T. DeGrand, A. Hasenfratz, P. Hasenfratz, F. Niedermayer, FP action for SU(3) | arXiv:`hep-lat/9506030`（NPB 454 (1995) 587/615） | 摘要片段（"no `a^n` nor `g²a^n` cutoff effects"） | `k₀^{os}=∞` 的候选 |
| X10 | *Oscillatory terms in the domain wall transfer matrix* | arXiv:`0710.0425` | 摘要片段（DBW2 ⇒ 负本征模） | 4D 正性违例的第二个见证 |

---

*Status：任务问「是否存在定理 `RP ⇒ k₀ ≤ 2`」——**不存在，且其一般形式假**（CP 1997 给出 2D σ 的 RP + on-shell 改进并存构造，`[R]`）。仓库内**0 先例**（§1.3，检索式与命中齐备）。本轮能立的是**正锥刚性**：`c_L ≥ 0 ⇒ k₀^{shell} = k₀^{axis} = 2`（`[T]`+`[E]`：16 形状 × 3–6 方向 ≈124 点 `κ_L>0`、LP 不可行、4 方向 shell 改进需 ≥2 负系数），反向 `k₀^{shell}≥4 ⇒ c∉𝒞₊`。⇒ **`2/(m+2)` 是"RP-可证类（正锥）"的硬前沿，等级 `[C]`；它不是 RP 的推论，因此拿不到 `[T]` 级 no-go。** 另外本轮钉死 `k₀` 的三义（`shell`/`axis`/`os`），并指出树级 LW 的 `k₀^{axis}=4` 但 `k₀^{shell}=2`（勘误 E-3，收窄 `GK3` 的 `k:2→4`）。4D gauge 正锥外的"特制改进"是**开问题**（`[O]`），三条结构性理由解释为何 2D 的约束机制搬不过来。**No QED；Clay 0%。未 commit；无 MC；只写本文件。***

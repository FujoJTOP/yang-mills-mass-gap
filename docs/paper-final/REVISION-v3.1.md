# YANG-MILLS-MERGED-V3 → v3.1：优化记录

对象：`docs/paper-final/`（`main.tex` + `sec01`–`sec06` + `refs.tex`），输出 `main.pdf`。
编译：`pdflatex main.tex` × 3（`thebibliography`，无需 bibtex）。结果：**46 页，0 error，0 overfull，0 undefined citation/reference，0 hyperref token warning**（v3：41 页，7 underfull，30 hyperref token warning）。

四条硬指标（本轮全部为"改前 → 改后"）：

| 项 | v3 | v3.1 |
|---|---|---|
| 正文 `\cite` | **0** | **76**（书目 76 条全部接通，无孤儿） |
| 重复章节号 | 有两套 §3/§4/§5/§6 | 连续 1–11 + 附录 A |
| 无编号表格 | 6 张（裸 `center`） | 11 张全部编号 + 题注 |
| 编译告警 | 30 hyperref token + 7 underfull | 0 + 5 |

---

## 一、加入引用（0 → 76）

v3 的书目是完整的、但**正文一条没用**。本轮按"哪里是继承、哪里是自己证"来放引用，而不是撒胡椒面：

1. **新增 §6 "Where the reduction sits in the literature"**（Part I 末）。五段分类：框架（`jaffe`,`sw`,`kl`,`os1973`,`osterwalder-seiler`）／强耦合质量间隙不属本文（`seiler`,`bgdf`,`gjs`,`fo`,`kp`,`szz`）／RG 路线（`balaban`,`balaban-conv`,`mrs`）／分析工具（`prekopa`,`bl`,`kato`,`reed-simon`,`vnw`,`temple`,`feller`,`js`,`sinclair`,`bbs`,`chatterjee`）／其余通往质量标度的路线（格点 `wilson,creutz,luscher-stepscaling,luscher-flow,morningstar,aoki`；大 N 与圈方程 `thooft-planar,eguchi-kawai,makeenko-migdal`；弦 `polchinski-strominger,luscher-weisz`；全息 `maldacena,gkp,witten,witten-thermal,polchinski-strassler,sakai-sugimoto`；bootstrap `elshowk,prv`；resurgence `abs`；β 函数与 a-定理 `caswell,vrvl,banks-zaks,polchinski,komargodski-schwimmer,dksst,zamolodchikov,osborn,hofman-maldacena,hofman-li`；红外方案 `gribov,singer,zwanziger,vsmekal,fischer-alkofer,huber,cornwall`；受限模型中的禁闭 `gopfert`）。
2. **正文就地引用**（不改结论，只标注归属）：§1.1（KL/反射正性）、§1.2（变分上界）、§1.3(i)（Efron–Stein 方差、输运–熵判据）、§1.3(ii)（Lüscher 体积依赖）、§1.3(iii)（Savvidy/Nielsen–Olesen/Leutwyler 负模）、§1.5（格点谱学的算符依赖）、§7 五环（KP、谱定理、KL/O–S）、§7.2（von Neumann–Wigner 水平交叉）、§7.6（Griffiths/Bochner/Ginibre/Kelly–Sherman）、§7.7（SZZ、O–S）、§8.2（Creutz 的 SU(2) 参考值）、§8.4（Grade 提升为 T）、§8.5（Morningstar 的变分谱学）、§9.1（Clay 声明）、§9.2 表内 4 行 + 新增"表后不等式的出处"一段（12 条）、§10.1（Prékopa/Brascamp–Lieb）。

> 说明：正文里 `\code{U01.md}, §3` 这类"文件 + 行号"**不是**文献引用，是本文的 provenance 规则，一律保留。

## 二、准确性核验（抓到 8 处数/文不一致 + 3 处书目硬错）

内部算术复核（逐条重算，非目测）：

1. **过报倍数配对错**。正文/摘要/§11.3 把 `0.018394` 配 `4.84×`；实测 `0.018394/0.003777 = 4.87`，`4.84` 属更早的 `0.018269`（`0.018269/0.003777 = 4.837`）。摘要、§7.1.1 的两行比值、§11.3 条目 1、修订说明统一改为**按各自分母**印刷。
2. **`0.25%` 与 `0.24%` 打架**。§9.2 第 7 行写 `0.24%`（`(0.5700−0.5686159)/0.5700 = 0.2428%`，比值 `0.9976` 亦为 0.24%），§10 表第 11 行写 `0.25%`。统一为 **0.24%**。
3. **"四个反例"实为 2 反例 + 2 算术证书**。§7.5 正文自己写明第 3、4 行"N 证书算术"，而摘要与节标题都叫"four counterexamples that kill four natural implications"。标题改 **"Four machine-verified lemmas: two counterexamples, two certificates"**，摘要同步。
4. **"isolation 已在 Lean 4 验证"过头**。形式上验证的是**有限支撑影子**（Proposition C，`I2 ⟺ SecExists`），§7.2 的两条限制自己写明了。摘要改为 "a finite-support shadow of the isolation, verified as a Lean 4 proposition with no sorry"。
5. **`TV` 归一化**。§1.2 表列 `‖μ_h−μ_0‖_TV = 0.600`，按标准 TV 应为 `0.300`（表中用的是 `L¹ = 2·TV`）。不修改数字，改**标注归一化**，题注写明 `∫|dμ−dν|` 是通常 TV 的两倍。
6. **`m₁` 与 `m₁²` 记号冲突**。§1 令 `m₁² := inf supp ν`（ν 载在 `m²` 上），§8.1 却同屏写 `C(t)=∫e^{−mt}dν(m)` 与 `ν({m₁²})`。§1.1 增一段记号约定（`ν` 载 `m²`，`C(t)=∫e^{−mt}dν(m²)=∫e^{−√s t}dν(s)`），`eq:m-defs` 跟着改。
7. **`#35/#36` 是 item 号，不是表行号**。§11.4 写 "#35/#36"，而表里 35/36 行是 class-level residual，真正的行是 **26/27**。改为 "corpus item 35/36（表 9 第 26/27 行）"。
8. **`N09-1` 等级与其自身表述矛盾**（见下节 L3）。

书目核验（写 `docs/paper-final/BIBCHECK-REPORT.md`，工具状态逐条记录；12 条抽检）：

| 键 | 判定 | 处置 |
|---|---|---|
| `mrs` | **DOI 错**：`10.1007/BF02079397` 指向一篇**地下水**论文（Fluid Dynamics 30 (1995) 758） | 改为 `10.1007/BF02097397`（数字换位 `797↔973`）。卷/页/年/标题原样 |
| `gopfert` | 年份错：应为 **1982**（CMP 82(4) 卷头 Jan 1982） | `1981 → 1982`，补全标题与 DOI |
| `balaban` | `119 (1988) 243–285` 被塞进 "Large field renormalization I–II" 的括号，实为**独立论文** | 拆出新键 `balaban-conv`；I/II、Large field I/II 卷页经 NII-REO + Crossref 双证无误 |
| `szz` | 号对（arXiv:2204.12737 = Shen/Zhu/Zhu, CMP 400 (2023) 805–851），但正文表述**把假设当常数** | 见下 |
| `sw` | 不是可引书目（"as cited in the Clay statement"） | 换成 Streater–Wightman 书（Benjamin 1964 / Princeton Landmarks 2000） |
| `dksst` / `hofman-li` / `fo` / `bgdf` / `kp` / `osterwalder-seiler` / `luscher77` | 卷页年全对 | 补标题/终止页/DOI |

**`1/48` 的表述更正（实质）**：SZZ 证明的是"当 `|β| < 1/48`（`d=4`, `SU(N)`）时 log-Sobolev 与 Poincaré 成立"，`1/48` 是**窗口边界**，不是"log-Sobolev 半径"这个常数。v3 的 "the radius 1/48 is Shen–Zhu–Zhu's log-Sobolev radius" 已改为 "the edge 1/48 of the strong-coupling window is the threshold of SZZ's log-Sobolev and Poincaré inequalities"。SZZ 是否另给显式 LSI 常数 `c e^{−κβ}`：**UNVERIFIED**（未取全文）。

## 三、逻辑漏洞（6 处，其中 3 处是实质性的）

**L1｜leg (i) 写的是"点态正"，靶要的是"下确界正"。** v3 的 `Clay ⟺ (i)+(ii)+(iii)` 里 (i) 为 `∀(L,a) w₁>0`，而 §9.1 的靶是 `inf_{L,t} w₁ > 0`；§1.2 全节论证的正是"点态正 ⇏ 一致地板"。(i) 已改写为 `inf_{L,t} w₁(L,2.5,t)>0`（`a=a(β)`），点态形式单独标注为更弱、且是测量唯一能给的（表 1 第 (i) 行）。

**L2｜"Clay ⟺ 三腿"与 §9.3 第三层直接冲突。** §9.3 说 `w₁` 与 `m` 在 2D YM、M6 上**方向相反**、完全解耦；若如此，`w₁` 侧的地板**不可能**等价于 Clay。v3 把该等价当无条件陈述印出。处置：左端改为 **`Clay`（on the atom route）**，并在框内写明三条限定 (a) 路线 / (b) 算符 / (c) 地板而非符号，指向 §9.3；(a) 明说"在该类之外等价为假，且不能靠加强 (i)–(iii) 得到"。摘要同步改为 "exact *on the atom route* … outside that class it is false"。

**L3｜`N09-1` 标 `[T]`，但它自己的陈述带印刷假设。** 标题页的等级定义是"带显式假设的定理 = `[C]`"。`N09-1` 明写 "conditional on the M6 specification and the L→∞ reading"，却打了 `[T]`。**这是本轮唯一的等级下调：`[T] → [C]`**，并在 §7.3 增一段说明（"与窗口过报同种性质，就地更正而非悄悄改"）；§7.4.1 里"`N09-1` 是 `[T]` 所以等级不转移"的论证改为"两者载体不同、各自带对方没有的假设，故等级不转移"。

**L4｜§7 导语的计数与它的六小节不符。** "two are theorems in the strict sense / one [C] / one [M] / two carry corrections" 漏掉了后加的 block-Loewner 小节。改为：三 `[T]`（I1、I2、block-Loewner）／两 `[C]`（N09-1、SG1-1）／一 `[M]`／两处更正。

**L5｜`\cos²θ` 见证在 §8.4 标 `[M]`、在 §1.5 标 `[T]`。** 同一命题两个等级。§8.4 改为 `[T]`（证明在 §5），数据层面另标 `[LE]`。

**L6｜§9.2 第 15 行（`Δ₀`）归在 "Circular"，但其"死因"写的是两个反例。** 该行死因补上真正的循环理由（"floor 就是靶换名，见 §3.3(iii)"），并说明为何仍留在该组。

另：§11.1 里 leg (i) 一句"as written, is `[E]` … `[O]`"两个等级相邻出现，已改为"拆分：点态 `[E]`、一致形式 `[O]`"。

## 四、排版与可读性

- **章节号**：源文件在 Part III 重置 section 计数，导致**两套 §3–§6**（旧 TOC 里 "6 What is established" 与 "6 The verdict" 并存）。改为全篇连续 1–11 + 附录 A；正文所有 `§N` 硬编码改成 `\ref`。
- **表格**：6 张裸 `center` 表（§1.1 legs、§1.2 紧族、§1.5 算符表、§7.1 五洞、§10 两处核验、§11 二十路分类）全部改为 `table` 浮动体 + 编号 + 题注（现 Table 1–11）；两个 `longtable` 补 `\caption`（Table 8、9）。
- **PDF 书签**：`\LT/\LC/\fcolorbox` 与行内公式进标题导致 30 条 hyperref token warning，全部用 `\texorpdfstring` 包裹（显示 chip、书签用 `[T]`/`SG1-1`）。现为 **0**。
- **`\src` 出处行**：改为 `\noindent` + 局部 `\sloppy`，消掉缩进与两处 overfull（`19.6pt`、`87.5pt`）。
- **版式**：`\part` 加"小体字母 Part N + 全宽细线"，Part II/III 补 `\clearpage`（v3 里 Part II 挤在页面下半）；标题页加分隔线、`v3.1` 标识、2020 MSC 与关键词；页眉 `(v3) → (v3.1)`。
- **题注位置**统一在表下；`tocdepth=3`（收 subsubsection、不收 paragraph 噪声）；超链接加书签编号/展开。
- **修订说明**：TOC 后新增 "Note on this revision (v3.1)" 框（R1–R10），把上面每一类更正写在读者能看到的地方 —— 与 §10 台账同一伦理。

## 五、未解决 / 明确留白

- `szz` 是否给出显式 LSI 常数：**UNVERIFIED**（未取全文；条目本身号、刊、卷、页已双源确认）。
- Balaban I 的副标题：只见到检索片段，条目**只写主标题 + `I`**，不补副标题。
- Streater–Wightman 页数两源不一致（x+181 / viii+181）：条目不写页数。
- 台账 9 条 class-level residual（行 30–38）仍是残项，本轮未伪造身份填补。
- v3 的 `p8`、`V6`、`F07`、`ADJ6`、`Q13` 五个短标签在 §10 路径约定里**明写为不可解析**（改为"读者请 grep"），不编造路径。

## 六、改动文件

| 文件 | 内容 |
|---|---|
| `main.tex` | preamble（`\src`、`\part`、题注、tocdepth、hyperref 元数据）、标题页、摘要、TOC 后新增 v3.1 修订说明 |
| `sec01-reduction.tex` | 记号约定、(i) 改写、`Clay on the atom route` 三条限定、表 1/2 成表、Savvidy 段、新增 §6 文献句 |
| `sec02-established.tex` | 导语计数、五环引用、`N09-1` 降 `[C]`、表 4、§7.5 标题与摘要对齐、`\code{szz}`→`\cite` |
| `sec03-measurement.tex` | `eq:m-defs` 记号、Grade 提为 `[T]`、Creutz/Morningstar 引用、去掉 section 计数重置 |
| `sec04-survey.tex` | 表 8 题注、第 2/4/11/12/15 行引用与死因、表后"不等式出处"段、§9.3 加路线限定 |
| `sec05-ledger.tex` | 路径约定补五标签、`0.24%`、表 9/10/11 成表、Prékopa 引用 |
| `sec06-verdict.tex` | 全部 `\S N` → `\ref`、`#35/#36`、`4.87`、`cM` 框改中性、附录标注 |
| `refs.tex` | 12 条记录修正/补全（含 `mrs` DOI、`gopfert` 年份、`balaban` 拆分、`sw` 换书目）、顶部说明改写 |
| `BIBCHECK-REPORT.md` | 书目核验原始报告（含工具状态、检索式、命中数） |
| `../position-current.html` + `position-r179.html` | 按仓库 §6 记录本轮：`N09-1 [T]→[C]`、`4.84→4.87`、`0.25%→0.24%`、引用 0→76、三处书目硬错 |

**结论层面**：本轮**没有**新增结论，`Clay` 仍为 **0%**；唯一等级变动是 `N09-1` 由 `[T]` 下调为 `[C]`（更诚实，不是更强）。

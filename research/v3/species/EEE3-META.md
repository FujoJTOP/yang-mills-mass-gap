# EEE3-META — 一个质量隙证明需要哪些成分；4D 缺哪一条

> 结构调查（不改既有文件、不 commit）。只交 `[T]`/`[NG]`。每条引用给 arXiv/DOI + 【读摘要】/【读正文】。
> 行数 ≤90。仓库既有同题工作：`research/v3/species/T23-TEMPLATE.md`（C1–C4 承重步分类）——本文件以它为输入，不重复其数表。

## §0 判词

`[T]` **一个质量隙证明 = 三件套**：① **正性载体**（反射正性/OS 重构，或正 Gibbs 权）；② **一个承重步**，把载体变成与 `(L,a)` 无关的指数衰减率；③ **一个 regime**，使②的前提成立。已证样本的②只有三形态：**C1** 转移算子=已知算子（可积/热核）、**C2** 稀薄缺陷气体+屏蔽、**C3** 小活动度团簇展开。

`[T]` **4D 纯 SU(N) 缺的不是 (i)–(vii) 中的某条工具，而是 (vii) 的跨 regime 延拓**。C3 在 4D **已成立**（SZZ：任意 `d>1`，`SU(N)` 需 `|β|<1/(16(d−1))`，`d=4` ⇒ `β<1/48`，摘要逐字 "the infinite volume measure has a strictly positive mass gap"）。缺的是把 `β` 从 `0.0417` 搬到物理点 `λ=g²N=3.2`。**口径（AGENTS §2.1）**：`77×` 是 `λ` 口径（`T23:137`），活动度口径为 `u(2.5)=0.507` vs `0.031`（16×）；两个数不通用。

`[NG]` **(ii) 单极凝聚在 4D 纯 YM 不是"未做到"，是对象不存在**：`π₂(SU(N))=0`，无有限作用量单极（无 Higgs 无 't Hooft–Polyakov 解）；abelian projection 的单极依赖规范固定，且中心相位对 N-ality-0 局域算子精确为 `1`（CV-1）。以 (ii) 为承重步的 4D 路线结构性关闭，不是待办。

`[NG]` **引用纠错**：arXiv:2207.05443 **不是** 2D 纯 YM 的隙证明。它是 *Gauge field marginal of an Abelian Higgs model*（2D 格点 Abelian Higgs、Villain 作用；loop expansion + UV 稳定性），CMP **405**, 143 (2024)【读摘要】。2D YM 可解性的可核引用见 §1。

`[NG]` **项目清单无需新增 4D 特有成分**：两条 4D 特有结构——(a) 中性盲（CV-1）、(b) `d=4` marginal `[g]=0` 无 gauge-covariant BPHZ——都已在册（§2.3）。

**下一手**：把 77× 缺口写成一条**具名引理**（对活动度 `u(β)` 的解析延拓，或对 `G″(k→0)` 的 `β`-一致下界），而不是再找第八条候选成分。

## §1 各模型的隙证明用了什么成分

| 模型 | 承重成分 | 其余候选的角色 | 引用与级 |
|---|---|---|---|
| **2D YM**（含 1+1D，二者同一模型） | **(iii) 可积性**：场空间有限维 ⇒ 转移算子=紧群热核，Peter–Weyl 对角化 ⇒ Casimir 谱 | (i) 面积律是**输出**；(ii)(vi) 平凡（`F∧F≡0`）；(iv) 非必需 | Gross–Witten, PRD **21** (1980) 446, DOI 10.1103/PhysRevD.21.446【读摘要】；Witten, CMP **141** (1991) 153, DOI 10.1007/BF02100009【读摘要】`[T]` |
| **N=∞ 的 2D** | **(iii)+(iv)**：`N=∞` 使 2D 配分函数用 steepest-descent 显式求值，得确定性极限 + `λ=2` 三阶相变 | 隙仍由 (iii) 的 Casimir 谱给；(iv) 不加新机制 | 同上（Gross–Witten 摘要逐字 "explicitly evaluated for all fixed λ=g²N … third-order phase transition, at λ=2"）`[T]` |
| **3D compact U(1)** | **(ii)+(vi)**：`π₂(U(1))=Z` ⇒ 有限作用量单极 ⇒ Boltzmann 压制 `ζ≪1` ⇒ Debye 屏蔽；严格版：Villain 作用对偶到 3D `Z`-铁磁体/Coulomb 气体 | (i) 弦张力是**输出**（"α nonzero for all values of g²"）；(iii)(iv) 不用 | Polyakov, NPB **120** (1977) 429, DOI 10.1016/0550-3213(77)90086-4【读摘要】`[T]`；Göpfert–Mack, CMP **82** (1982) 545, DOI 10.1007/BF01961240【读摘要】`[T]`（连续统极限=质量 `m_D` 的自由标量 ⇒ 隙严格） |
| **3D SU(N)** | 只有 **(vii) C3**（强耦合，小活动度）；连续统/弱耦合无 | (ii) 无单极对象；(iii) 无解 | SZZ, arXiv:2204.12737（`d=3` 在定理内）【读摘要】`[T]`；3D 弱耦合 `[NG]`（无证明） |
| **4D SU(N)，β<1/48** | **(vii) C3**：正性载体（RP）+ 小活动度聚合物/团簇展开 + LSI/Poincaré | (i) 面积律是输出；(iv) 大 `N` 因子化在同一 regime 内成立（摘要逐字 "factorized correlations"） | Osterwalder–Seiler, Ann. Phys. **110** (1978) 440, DOI 10.1016/0003-4916(78)90039-8【读摘要】`[T]`（摘要给正性/转移矩阵/禁闭界；**正文的 gap 条款未逐字核**）；SZZ arXiv:2204.12737【读摘要】`[T]` |
| **4D N=2 SU(2)（Seiberg–Witten）** | **(iii)′ 精确可解性**（SW 解 + 超对称保护）；(ii) 单极凝聚是**结论**，不是承重步 | (i)(iv)(v)(vi) 均非承重 | Seiberg–Witten, NPB **426** (1994) 19, DOI 10.1016/0550-3213(94)90124-4 / arXiv:hep-th/9407087【读摘要】`[T]`（摘要逐字 "with a suitable perturbation confinement is described by monopole condensation"） |
| **4D 纯 SU(N) 弱耦合** | **无**。`(ii)` 只有 't Hooft 的拓扑算符/相结构论证：3+1 维四相含"无 Higgs 而永久禁闭" | 该论证**不给隙**（摘要四相含"有质量粒子的中间相"）；单极图像属其 abelian-projection 后续，本报告未核 | 't Hooft, NPB **138** (1978) 1, DOI 10.1016/0550-3213(78)90153-0【读摘要】`[NG]`（作为隙证明） |

**脚注（不进主表）**：2D `CP^{N−1}` 是"真局部自由度+真隙"的 2D 样本，承重步为 (ii) 分数瞬子气体稀薄（`ζ∝e^{−cN}`）+(iv) 大 `N`；外部原文**未核**，仅仓内转述 ⇒ `[NG]`（作为本报告引用）。

## §2 4D 纯 SU(N) 缺哪一条：逐条判

| 候选 | 4D 可用？ | 承重性 | 判 |
|---|---|---|---|
| **(i) 面积律** | 可用（`β<1/48` 严格；弱耦合未证） | **不承重**：面积律给 1-形式弦通道 `√σ`；2D 有精确面积律而 `f₁²(trF²)≡0`（`B03.md:49`、`T23:95`） | `[NG]`（充分性）；下一手：需 `0⁺⁺`–弦 overlap 引理，未证 |
| **(ii) 单极凝聚** | **对象不存在** | — | `[NG]`；下一手：若复活须先造规范不变单极算符（存在性未证） |
| **(iii) 可积性** | 无已知结构 | 2D 的可积性靠"有限维场空间"；4D 转移算子含非二次磁项、逐点无穷维 | `[NG]`；下一手：无（作为承重步关闭） |
| **(iv) 大 N 因子化** | 强耦合**有**（SZZ 摘要）；连续统无 | **不承重**：因子化只把同一 regime 的隙重述；4D Makeenko–Migdal 不闭合 | `[NG]`；下一手：4D `N=∞` master field 闭合 = 未证 |
| **(v) 质量形变** | 无规范不变相关质量 | `gap(m₀>0)` 只给 `inf` 的上界，不能 floor 未形变隙（`DEFORM-DECOUPLE-4D:41`）；adjoint Higgs 去禁闭 | `[NG]`；下一手：无 |
| **(vi) 拓扑荷紧致性** | 有 `π₃(SU(N))=Z`（瞬子） | **不承重**：瞬子是有限尺寸事件、无长程场 ⇒ 无稀薄等离子体屏蔽；拓扑量在 2D gapped 模型恒 0（`B04-INDEX.md:112` `[T]` 二分法） | `[NG]`；下一手：无（拓扑路线已整体判死） |
| **(vii) 其他 = C3** | **有，且 `d=4` 在定理内** | **缺 regime**：`u(2.5)=0.507` vs `0.031` = 77×；展开是渐近级数，半径外无更强估计 | `[T]`（定理存在）/`[NG]`（跨 `β`） |

**§2 结论**：4D 缺的是 **(vii) 的跨 `β` 延拓**；`(i)–(vi)` 全部**不承重**，其中 `(ii)` 与 `(vi)` 是对象不存在，`(iii)` 无对象，`(v)` 改理论，`(i)`/`(iv)` 只给别的通道或同一 regime。

**§2.3 Q3：有没有 4D 特有、清单未覆盖的成分？** `[T]` **无**。(a) 中性盲 `CV-1`（中心相位对 N-ality-0 局域算子 `≡1`）与 (b) `d=4` marginal（`[g]=0`，无 gauge-covariant BPHZ）都已在册。七条候选可完整映射到 C1–C4：`(iii)=C1`、`(ii)+(vi)=C2`、`(vii)+(i)=C3`、`(iv)⊂C4 的标签侧`；`(v)` 在任何纯 YM 样本里都不是承重步——它空的原因是"改理论"，不是"没做到"。

**下一手**：唯一未被否证的新方向不是"新成分"，而是"**中性 `0⁺⁺` 的规范不变局域算符**"（CV-1 只排除中心相位，未排除全部中性算符）；这是造算符的活，不是列清单的活。

## §3 台账（AGENTS §1 四要素）

**① 工具状态**：`pwsh` OK（3 次，含目录列举）；`glob` OK；`grep` OK（4 次，均返回有效计数，末次 18 命中）；`read` OK；`web_fetch` 有效 8 次，`TOOL-FAIL` 1 次（见④）。

**② 实际检索式**

- `grep "Migdal|2207\.05443|Chandra|Chevyrev"` @ `research/**/*.md`
- `grep "EEE[0-9]|EEE-"` @ `research/v3/**/*.md`
- `grep "瞬子|instanton|π_?3\(|pi_3"` @ `research/v3/**/*.md`
- `grep "Curci|质量形变|mass deformation|topologically massive|Deser"` @ `research/**/*.md`
- `read research/v3/species/T23-TEMPLATE.md:88–217`
- `web_fetch https://arxiv.org/abs/{2207.05443, 2204.12737, hep-th/9407087}`
- `web_fetch https://inspirehep.net/api/literature?...doi:(5 个 DOI)`
- `web_fetch https://inspirehep.net/api/literature?...doi:(2 个 DOI)`
- `web_fetch https://inspirehep.net/api/literature?q=t "quark confinement and topology of gauge theories"`
- `web_fetch https://inspirehep.net/api/literature?q=t "confinement of static quarks" and t three-dimensional`
- `web_fetch https://api.crossref.org/works?query.bibliographic=Driver…`

**③ 命中数**：449（第 1 式，显示截断于 250）；**0**（第 2 式；grep 时刻无 `EEE` 文件，随后并发代理写入 `EEE1-P3CHAN`/`EEE2-F1SQ` ⇒ 本文件非系列首作，且与二者不冲突：它们在 `2⁺⁺` 道/算符层，本文件在成分层）；53；18；T23 段 130 行；arXiv 3/3 摘要；INSPIRE 5/5、2/2（其一为无关核物理论文，按 DOI 排除）、0、1；crossref 失败。

**④ 失败标记与未核清单**

- `TOOL-FAIL`：crossref 查 Driver 1989 → `Error: unsupported content type "unknown"`；**该文未引**。
- INSPIRE 两式为 **0 命中但工具 OK**（非"不存在"）：`t "quark confinement and topology of gauge theories"`；`t "three-dimensional compact electrodynamics"`。
- **未逐字核正文**：Migdal 1975（旧刊无 DOI，未引其文本）；Osterwalder–Seiler 的 gap 条款（摘要未含 "gap" 字样）；'t Hooft 1981 abelian projection；Witten 1979 `CP^{N−1}`；Driver/Sengupta/Lévy 的 2D 测度严格化（均未引）。
- **行数**：本文件 90 行内（末行即此）。

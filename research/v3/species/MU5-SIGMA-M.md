# MU5-SIGMA-M — `σ > 0 ⇒ m > 0` 的严格性裁决

**任务**：判 `σ>0`（弦张力）是否蕴含 `m>0`（质量隙），以决定 `σ = 2g²N·Λ²` `[C]` 能否升级为 Clay 级结果。
**裁决日**：2026-11-25 · **标签**：`[A]` 已逐字核对并成立 · `[B]` 领域公认但本次未逐字核原文 · `[C]` 条件性/本战役结果 · `[X]` 未查到（附检索式）

---

## 0. 判死线（第一段）

### **(b) —— 无已证蕴含；`σ` 只是 confinement 侧，与 Clay 差一层。且比 (b) 更坏：文献中已证的箭头是反的。**

三条硬结论：

1. **文献里没有"已证 `σ>0 ⇒ m>0`"这个定理。** 检索 Osterwalder–Seiler 1978 的原文摘要，其已证清单是：反射正性 ⇒ 正自伴转移矩阵；强耦合无限体积极限的存在性与解析性；**Wilson confinement bound（面积律）**；格点 Higgs 机制的严格处理。**摘要中不含"mass gap"字样** `[A]`（逐字见 §1.1）。所谓"OS 证明了质量隙"这一常见转述，在摘要层面无法核对——**这是一个必须停止传播的引用误记**。
2. **已证的是反方向：`m>0 ⇒ σ>0`。** Chatterjee（Comm. Math. Phys. 2021）严格证明：**指数关联衰减（+ 非平凡中心）⇒ 中心对称不破缺 ⇒ 该理论是 confining 的** `[A]`（逐字见 §1.2）。即：**谱隙是面积律的充分条件输入，不是它的推论。** 对我们的链条，这是"倒扣一分"而非"白送一分"。
3. **`σ = 2g²N·Λ²` 的 `[C]` 链到此终止。** 它是一条 confinement 侧陈述（线性势 / 面积律系数），**不可**升级为 `m_{0++}>0`。要接 Clay，必须另造"面积律 → 局部算符谱隙"这一层机器——而这一层目前**不存在于文献**（§2 列出它需要什么）。

> **一句话**：`σ>0` 与 `m>0` 之间不是"差一步推导"，是**缺一个已证的定理，且现有的定理朝相反方向走**。把 `σ = 2g²N·Λ²` 说成"Clay 级"= 把 confinement 指标冒充 mass-gap 定理。

---

## 1. 文献：逐字引用（全部经原文摘要页核对）

### 1.1 Osterwalder & Seiler 1978 — 面积律已证，隙未见

- K. Osterwalder, E. Seiler, **"Gauge Field Theories on the Lattice"**, *Annals of Physics* **110** (1978) 440–471.
- **DOI**: `10.1016/0003-4916(78)90039-8` · INSPIRE record `120117` · 核对途径：INSPIRE HEP API 返回的 Elsevier 官方摘要 `[A]`
- **逐字引用**（Elsevier abstract，全文照录关键句）：

  > "We discuss the lattice approximation of Yang-Mills fields and of fermion fields in the Euclidean setup and we verify physical positivity for the Schwinger functions of these approximations. **This implies the existence of a positive self-adjoint transfer matrix.** We then **prove existence and analyticity of the infinite volume limit of strongly coupled Yang-Mills theories** on the lattice and **we verify Wilson's confinement bound**. Finally we present a rigorous treatment of the Higgs mechanism in lattice gauge theories."

- **判读**：
  - "confinement bound" = Wilson 面积律下界 ⇒ **`σ>0` 在强耦合已证** `[A]`。
  - 摘要**只到转移矩阵的存在性**，未声明谱隙；"positive self-adjoint transfer matrix" 是**谱隙的必要前提**，不是谱隙本身。
  - 因此：**"OS 1978 给出 `σ>0 ⇒ m>0`" 这一说法在可核文本中查无实据** `[X]`（检索式：`"Gauge field theories on a lattice" Osterwalder Seiler abstract mass gap`；`string tension mass gap implication`）。

### 1.2 Chatterjee 2021 — 已证的箭头是 `隙 ⇒ confinement`

- S. Chatterjee, **"A probabilistic mechanism for quark confinement"**, *Comm. Math. Phys.* (2021).
- **arXiv**: `2006.16229`（v5, 2021-03-26）· **DOI**: `10.1007/s00220-021-04086-y` · 核对途径：arXiv abs 页 `[A]`
- **逐字引用**（abstract）：

  > "This article gives mathematical confirmation of this heuristic, by rigorously defining of center symmetry in lattice gauge theories and **proving that a theory is confining when center symmetry is unbroken**. Furthermore, **a sufficient condition for unbroken center symmetry is given: It is shown that if the center of the gauge group is nontrivial, and correlations decay exponentially under arbitrary boundary conditions, then center symmetry does not break.**"

- **判读**：**"correlations decay exponentially"（即隙型输入）⇒ 不破缺 ⇒ confining**。这是文献中唯一严格闭合的蕴含，**方向为 `m>0 ⇒ σ>0`**，与我们的需求相反 `[A]`。
- **附带**：该文正文明确指出面积律在弱耦合仍是开放问题（"Proving that the area law holds at weak coupling remains a largely open problem above dimension two"）`[B]`（PDF 文本片段，未做页面级逐字定位）。

### 1.3 Chatterjee 2018 — 4D 质量隙仍是"believed"

- S. Chatterjee, **"Yang-Mills for probabilists"**, arXiv `1803.01950` · DOI `10.48550/arXiv.1803.01950` `[A]`（arXiv abs 页逐字）
- **逐字引用**（abstract）：

  > "The rigorous construction of quantum Yang-Mills theories, especially in dimension four, is one of the central open problems of mathematical physics."

- **判读**：2018 年概率学派权威综述将 4D 构造（含隙）列为**未解**。与 1.1/1.2 合并 ⇒ **不存在把 `σ>0` 兑换成 `m>0` 的现成定理** `[A]`。

### 1.4 Seiler 1982 专著 — 隙的讨论所在，但本次未取得可逐字核对的公开文本

- E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, **Lecture Notes in Physics 159**, Springer, 1982, iv+192 pp. `[B]`（书目信息经 Springer/图书馆编目核对：LNP 159，1982）
- **状态**：该书是"强耦合下面积律 + 隙由同一簇展开给出"的标准出处 `[B]`，但**本次未取得可公开逐字引用的原文段落** `[X]`（检索式：`Seiler "Gauge Theories as a Problem of Constructive Quantum Field Theory" mass gap chapter`；`Seiler Lecture Notes in Physics 159 Springer DOI 1982`）。
- **纪律声明**：在拿到该书的逐字段落之前，**禁止**以"Seiler 已证 `σ>0 ⇒ m>0`"作为引用。目前可安全陈述的只有：**强耦合区面积律与隙由同一收敛簇展开分别给出**（两条平行结论，非蕴含）`[B]`。

---

## 2. 逻辑：从面积律到谱隙，中间缺什么（清单）

面积律 `W(C) ~ e^{−σA}` 是**大圈渐近**；谱隙是**转移矩阵本征值间隙**。二者之间需要下列**全部**环节，缺一不可：

| # | 缺的环节 | 为什么必需 | 现状 |
|---|---|---|---|
| L1 | **反射正性 / OS 正性** | 才能定义 `T = e^{−H}` 且 `H ≥ 0` 自伴 | OS1978 已给"positive self-adjoint transfer matrix" `[A]` |
| L2 | **转移矩阵谱定理 + 真空唯一性** | 隙 = `inf spec(H) \ {0}`；真空简并则"隙"无定义（Perron–Frobenius / 不可约性需要） | 强耦合有 `[B]`；一般情形无 `[X]` |
| L3 | **`L→∞` 与 `T→∞` 的极限交换次序** | 面积律是 `A→∞`（双向）渐近；隙是先在有限盒取 `T→∞`、再 `L→∞` 的 `inf`。**顺序不可交换**，交换性本身是一个定理级缺口 | **未见任何文献证明此交换** `[X]` |
| L4 | **从"静态源扇区"跨到"局部算符扇区"** | 面积律只测 `W(C)` 的**基态**能量 `E₀(R) = σR + O(1/R)`；`m_{0++}` 由**局部规范不变算符** `O(x)` 的连通关联 `⟨O(x)O(0)⟩_c` 的衰减率给出。**两个扇区不同**，面积律对后者无发言权 | **无蕴含定理** `[X]`（本节核心缺口） |
| L5 | **0^{++} 通道的算符重叠 `⟨0|O|1⟩ ≠ 0`** | 即使有隙，也需证明该通道**真的**被 `O` 看到（否则只是别的通道的隙） | 数值上成立 `[B]`；严格无 `[X]` |
| L6 | **OS 重构回到 Wightman/相对论理论** | 才有 Clay 意义上的"量子场论 + 质量隙" | OS1978 只做 Euclidean 侧 `[A]` |

**L4 是致命的，且可以量化**：面积律自身的修正项就是反证——有效弦理论给出

`V(R) = σR − π(D−2)/(12R) + O(1/R³)`（**Lüscher term**；Lüscher, Symanzik, Weisz, *NPB* **173** (1980) 365 `[B]`，本次经二手引文核对，未逐字核原文 `[X]`）

即：**产生面积律的那个扇区，其激发能级间距是 `O(1/R) → 0`**（弦的振动态）。**面积律在 `R→∞` 极限下自带一个无隙方向。** 因此"面积律 ⇒ 隙"在结构上不可能由面积律单独提供——隙必须来自**另一个扇区**（glueball/局部算符），而面积律恰恰对该扇区沉默。

> **精确缺失的一步**（写死）：**L4** —— "面积律（静态源扇区基态能量的线性增长）⇒ 局部规范不变算符连通关联的指数衰减"。这一步**文献中不存在**，且**不能**由 L1–L3 拼出。

---

## 3. 反例检查

### 3.1 "有面积律但无隙" —— 未查到明确反例，但**未被任何定理排除**

- **检索结果**：`[X]`（检索式：`"area law" "does not imply" "mass gap" lattice gauge`；`confinement without mass gap counterexample`）。**未查到**一个明确构造：反射正性格点规范理论中 `σ>0` 而 `spec(H)` 在真空上方连续到 0。
- **但注意**：查不到反例 ≠ 蕴含成立。当前状态是**"蕴含既未证、也未反证"**，即 §0 的 (b)。若有人造出该反例，则立即降为 (c)。
- **最接近反例的结构**：§2-L4 的弦激发谱 `O(1/R)`——**面积律所在的扇区本身无隙**。这已是"`σ>0` 不保证整个谱有隙"的结构性证据，只是它不构成"整个理论无隙"的反例（glueball 扇区可能仍有隙）`[A]`（结构性论证，非文献定理）。

### 3.2 "有隙但无面积律" —— **已证存在**（这是文献中真实存在的一侧）

- **基本表示物质 / Higgs 区**：`Fradkin–Shenker` 型解析延拓——规范理论 + 基本表示 Higgs 场，**有质量隙但 Wilson 圈为周长律**（无面积律）`[B]`。OS1978 摘要逐字含 "we present a rigorous treatment of the Higgs mechanism in lattice gauge theories" `[A]`，为该侧提供严格基础。
- **4D U(1) Coulomb 相**：周长律 + 无质量光子（无隙无面积律）`[B]`。
- **有限温度去禁闭相**：时间方向周长律；空间 Wilson 圈仍有面积律（磁禁闭）——说明**面积律与隙可以解耦地出现在不同扇区** `[B]`。
- **含义**：`σ` 与 `m` **不是同一枚硬币的两面**。存在 `m>0, σ=0` 的严格例子 ⇒ **`m` 不蕴含 `σ`**（与 §1.2 的 Chatterjee 定理不矛盾：后者要求非平凡中心 + 任意边界条件下的指数衰减；基本表示物质**显式破坏中心对称**，前提失效）。**两侧互不蕴含的图景是一致的。**

---

## 4. 量纲：`σ = m²` 是普遍关系吗？

**不是。** `σ` 与 `m²` 同量纲（质量²）**只说明可以相比，不说明相等**。

| 体系 | `σ` | `m` | `m²/σ` | 标签 |
|---|---|---|---|---|
| 4D SU(3) 格点 | `√σ ≈ 440 MeV` | `m_{0++} ≈ 1710 MeV` | `≈ 15` | `[B]`（标准格点值，本次未逐字核 Morningstar–Peardon 原文 `[X]`） |
| 4D SU(∞) | 同上标度 | `m_{0++}/√σ` 随 `1/N²` 缓变 | `≈ 11–12` | `[B]` |
| 2D YM | `σ = g²N/2`（Migdal/Witten 精确） | 谱 ∝ `g²N` | 纯数（模型依赖） | `[B]`，**未逐字核** `[X]`（检索式：`two-dimensional Yang-Mills exact string tension g^2 N/2 glueball mass spectrum`） |

**结论**：
- **`σ = m²` 在任何模型里都不成立**；4D 里差一个 `O(10)` 的因子。
- 唯一"两者同源"的情形是 2D YM：`σ` 与 `m` **都正比于 `g²N`**，故比值是纯数——这**不是** `σ = m²`，而是"同一耦合定标"。
- 4D 中 `σ` 与 `m` 仅通过 `Λ_QCD` 联系，比值 `m_{0++}/√σ` 是**模型依赖的无量纲数**，无普遍关系 `[B]`。
- **对本战役的直接后果**：即便 `σ = 2g²N·Λ² > 0` 成立，也**不能**由此读出 `m² = 2g²NΛ²`。任何"由 `σ` 直接给 `m`"的步骤都需要额外输入（§2-L4 或格点定标），不可白拿。

---

## 5. 对 `σ = 2g²N·Λ²` `[C]` 的处置

- `σ = 2g²NΛ²`（`[C]`，条件于 `ε=1/Λ`；5:1 票）**保持原标签不动**。本次裁决**不支持**任何升级。
- 它现在的位置：**confinement 侧、强耦合可比、量纲正确、但通向 `m_{0++}` 的路未铺**。
- **禁止的写法**（列入纪律）：
  - ✗ "由 `σ = 2g²NΛ² > 0` 得 `m > 0`"
  - ✗ "OS 1978 已证面积律蕴含质量隙"
  - ✗ "`σ = m²`，故 `m = √(2g²N)Λ`"
- **允许的写法**：
  - ✓ "`σ = 2g²NΛ² > 0`（`[C]`）⇒ 强耦合面积律 / 线性势系数为正（confinement 指标）"
  - ✓ "文献已证反向蕴含：指数关联衰减 ⇒ confining（Chatterjee 2021, CMP）"

---

## 6. 若要接 Clay，唯一可走的路（路径）

不是修补 `σ`，而是**造 §2-L4**——即本宪章第 1 条要求的"新机器"：

1. **不能**从 `W(C)` 出发。面积律已被证明对局部算符扇区沉默（L4 + Lüscher 项）。
2. **必须**造一个新单调量 `M`，满足：
   - `M` 定义在**局部规范不变算符**的连通关联上（不经 Wilson 圈）；
   - 对 YM **正**，对自由无质量场 `f₁²=0` **失效**（过酸测）；
   - 能把强耦合区的 `m(β)>0` 与 `σ(β)>0` **各自独立**搬运过 Savvidy 交叉点；
   - 搬运过程**不依赖**"面积律蕴含隙"这一未证命题。
3. **候选交接点**：本战役已有的 `f₁² = Z_O² f₁²` 乘性残数搬运——注意它搬的是**残数**，不是**圈渐近**；这正是它与 L4 兼容的原因（残数是局部量）。**这是唯一不自相矛盾的现有线索。**
4. **验收标准**：新机器必须同时给出 `m_{0++} > 0` **与** `σ > 0`，且**不**通过二者互推。

---

## 7. 检索记录与"未查到"清单（纪律留痕）

**实际使用的检索式**：`string tension mass gap implication` · `area law implies mass gap` · `confinement mass gap rigorous Osterwalder Seiler` · `Osterwalder Seiler "Gauge field theories on a lattice" 1978 Annals of Physics mass gap` · `Seiler "Gauge Theories as a Problem of Constructive Quantum Field Theory" mass gap chapter` · `"area law" "does not imply" "mass gap" lattice gauge` · `confinement without mass gap counterexample` · `string tension glueball mass ratio 2D Yang-Mills exact` · `two-dimensional Yang-Mills exact string tension g^2 N/2 glueball mass spectrum Gross Witten` · `Luscher term static potential 1/R string excitation spectrum effective string` · 以及 INSPIRE HEP API 直查（`a Osterwalder and a Seiler`；`t "Gauge field theories on a lattice"`）。

**核对到的**：OS1978 摘要（INSPIRE/Elsevier 全文照录）`[A]`；Chatterjee 2021 摘要 `[A]`；Chatterjee 2018 摘要 `[A]`；Chatterjee 2021 DOI `10.1007/s00220-021-04086-y` `[A]`。

**未查到 / 未逐字核**（`[X]`，逐条列出）：
1. 任何文献中**已证**的 `σ>0 ⇒ m>0`（一般情形或 4D SU(N)）— **查无**。
2. Seiler 1982 (LNP 159) 关于隙的**原文段落**（书未公开取全文）。
3. Osterwalder–Seiler 1978 正文中是否含"mass gap"定理（仅核到摘要；摘要不含）— **需查原文 PDF 才能定论**，在此之前不得引用。
4. "有面积律但无隙"的**明确反例** — 查无；状态为"未证未反证"。
5. Lüscher–Symanzik–Weisz 1980 (NPB 173, 365) 原文逐字（经二手引文核对 `[B]`）。
6. 2D YM 精确 `σ = g²N/2` 与 glueball 谱的原始出处逐字（`[B]`）。
7. 4D SU(3) 格点 `√σ≈440 MeV`、`m_{0++}≈1710 MeV` 的原始格点论文逐字（`[B]`）。

**检索负面信号（值得记录）**：`string tension mass gap implication` 等检索式的头部结果被大量 Zenodo/Preprints.org 的"YM mass gap 已解"自出版稿件占据（含声称"由面积律得隙"的稿件）。**这些稿件互相引用、无一经过同行评审，且其"面积律 ⇒ 隙"步骤均落在 §2-L4 的空缺上而不自知。本文件一律不引用。**

---

**裁决：`(b)`。`σ>0` 不蕴含 `m>0`（无已证定理、无已知反例、缺 §2-L4）；且文献已证的箭头为 `m>0 ⇒ σ>0`（Chatterjee 2021）。`σ = 2g²NΛ²` 保持 `[C]`，属 confinement 侧，距 Clay 尚缺 L4 这一整层机器。**

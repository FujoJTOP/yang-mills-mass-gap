# A5-LIT-UNIFORM — 有没有"从有限 `a` 的谱隙推出 `a`-一致谱隙"的已发表定理？

> 任务：只回答一个问题，并逐个核引用。检索工具：arXiv API（export.arxiv.org）、Crossref REST、web_search/web_fetch。
> 阅读深度标注：**[元数据]** = 只核到 Crossref/arXiv 的题名-作者-年-卷-页-DOI；**[摘要]** = 读到官方摘要；
> **[正文片段]** = 读到正文某段（非全文精读）。**凡未标注者，均未读到全文**。

---

## 0. 判死线（先给裁决）

**主判：(b)。** 已发表数学里**没有**一条定理从"有限格距 `a` 的谱隙"推出"`a`-一致、`L`-一致的谱隙"，也没有一条在 `a→0` 之后给出 `m_{0^{++}} > 0`。
连续极限这一步**在已发表文献中是空白**——不是我们搜漏了，是领域边界。

**形式层面落在 (c)。** 文献里唯一存在的"一致"是**两类不同的东西**，都不等于我们要的：

| 文献中的"一致" | 它一致的是 | 它是否给谱隙 |
|---|---|---|
| Balaban 的 UV 稳定性估计（`CMP` 102/109/116） | **有效作用量**在 RG 步数（即 `a`）上一致 | **否**：控制的是配分函数/有效作用的解析性，谱隙从未被推出 |
| Dimock–Hurd / Dimock 的 Balaban 型 RG 分析 | **`a→0` 的奇异控制** | **否**：明文"只做第一步 / 不控制余项" |
| Battle–Federbush、Federbush 团簇展开 | **体积 `L`**（热力学极限、团簇展开收敛） | 部分（对 P(φ)₂ 类模型），**与 `a` 无关** |
| Osterwalder–Seiler | **边界/反射正性 → OS 重构** | 提供"谱条件 ⟺ 指数聚类"的桥，**不含任何 `a` 命题** |

即：**"`a`-一致"这个词在文献里存在，但它挂在"有效作用量估计"上，不挂在"谱隙"上。** 把两者接起来的那一步，没人在已发表文献里做过。

**若只有 `Λ`-依赖形式 ⇒ 与 `t3` 的 `ma ≥ c·aΛ` 同型（见 §3）**，只能得到"随 `a` 趋零的界"，得不到 `a`-一致正下界。这与 (c) 完全吻合。

---

## 1. 四条（实为六条）候选，逐条核过

### C1. Osterwalder–Seiler, 1978 — 轴：**边界/反射正性 + OS 重构**
- K. Osterwalder, E. Seiler, "Gauge field theories on a lattice", *Annals of Physics* **110** (1978) 440–471.
- DOI: `10.1016/0003-4916(78)90039-8` **[元数据，Crossref 核过]**
- 统一/提供的是：格点 Wilson 理论的反射正性 ⇒ OS 正性 ⇒ 可重构 Wightman/谱数据。
- **不含 `a`-一致谱隙**：它是"若有隙则得指数聚类"的桥，不是"隙存在且一致"的定理。

### C2. Balaban, 1985 — 轴：**格距 `a`（UV cutoff），但对象是有效作用量**
- T. Balaban, "Ultraviolet stability of three-dimensional lattice pure gauge field theories", *Commun. Math. Phys.* **102** (1985) 255–275.
- DOI: `10.1007/BF01229380` **[元数据，Crossref 核过]**
- 一致的是：3D 纯规范理论配分函数/有效作用量在**多重 RG 步（= 逐层缩小的 `a`）**下的稳定性。
- **没有谱隙结论**；且是 3D。

### C3. Balaban, 1988 — 轴：**格距 `a`（4D，小场近似 + 耦合常数重整化）**
- T. Balaban, "Renormalization group approach to lattice gauge field theories"（I: generation of effective actions in a small field approximation and a coupling constant renormalization）, *Commun. Math. Phys.* **116** (1988) 1–22.
- DOI: `10.1007/BF01239022` **[元数据，Crossref 核过]**
- 同一系列另有 `CMP` **95** (1984) 17–40（`10.1007/BF01215753`）、`CMP` **99** (1985) 389–434（`10.1007/BF01240355`）、`CMP` **102** (1985) 277–309（`10.1007/BF01229381`）——**均核到元数据**。
- 一致的是：**RG 步数**上的估计（这就是"uniform in the lattice spacing"在文献中的实际含义）。
- **谱隙仍缺**：Balaban 程序止于 UV 稳定性，未接 OS 重构，未给 `m>0`。

### C4. Magnen–Rivasseau–Sénéor, 1993 — 轴：**UV cutoff 移除（固定 IR cutoff，规范固定 SU(2)）**
- J. Magnen, V. Rivasseau, R. Sénéor, "Construction of YM₄ with an infrared cutoff", *Commun. Math. Phys.* **155** (1993) 325–383.
- DOI: `10.1007/BF02097397` **[元数据，Crossref 核过；正文片段：projecteuclid PDF 首页——原文自述"we never try to lift"红外 cutoff，且以 Slavnov 恒等式作为构造正确性的检验]**
- 统一的是：**去掉 UV cutoff**（`a→0` 方向！）。
- **但不是我们要的东西**：(i) 红外 cutoff 永不解除；(ii) 规范固定，未构造规范不变 OS 理论；(iii) 不给规范不变 `0^{++}` 谱隙。

### C5. Battle–Federbush, 1982 — 轴：**体积 `L`（phase-cell 团簇展开）**
- G. A. Battle III, P. Federbush, "A phase cell cluster expansion for Euclidean field theories", *Annals of Physics* **142** (1982) 95–139.
- DOI: `10.1016/0003-4916(82)90230-5` **[元数据，Crossref 核过]**
- 同组：Battle–Federbush, *CMP* **88** (1983) 263–293（`10.1007/BF01209480`，hierarchical φ⁴₃）；Battle–Federbush, *CMP* **109** (1987) 417–419（`10.1007/BF01206144`）。
- 统一的是：**体积/热力学极限**上的展开收敛。`a` 不是它的轴。

### C6. Dimock（Balaban 型 RG，QED₃）— 轴：**`a→0` 的 UV 问题，但明言未完成**
- J. Dimock, "Quantum electrodynamics on the 3-torus — I", arXiv:`math-ph/0210020` **[摘要：读到官方摘要]**
  > "We start with the lattice gauge theory on a toroidal lattice and seek to control the singularities **as the lattice spacing is taken to zero** … In this paper we explore the **first step** in detail."
- J. Dimock, "Multiscale block averaging for QED in d=3", arXiv:`1712.10029` **[摘要]**；"QED on the 3-torus II — the RG flow", arXiv:`math-ph/0407063` **[摘要]**（"we **do not control remainders**"）。
- 轴 = `a`，方向对；**结论 = 没有谱隙，且 3D，且 QED**。这是最接近"统一 `a`"的已发表程序，但它自己声明未完成。

### 附：Wilson 1974 — 轴：**`a` 与 `Λ` 的标度关系（dimensional transmutation）**
- K. G. Wilson, "Confinement of quarks", *Phys. Rev. D* **10** (1974) 2445–2459. DOI: `10.1103/PhysRevD.10.2445` **[元数据，Crossref 核过]**
- 提供 `a = a(β)`、`Λ ~ a^{-1} f(g)` 的标度关系；**不是谱隙定理**，但是 §3 形式判定的出处。

---

## 2. 判定：有没有一条统一 `a` 的？

**没有。** 明确记录：

- **未查到**任何以"谱隙在格距 `a` 上一致"或"`a→0` 极限后 `m_{0^{++}} > 0`"为**结论**的已发表定理。
- 检索式（可复现）：
  1. arXiv API: `abs:"mass gap" AND abs:"lattice gauge" AND abs:uniform` → 命中 **1** 条：arXiv:`1507.06624`（Milsted, *PRD* 93, 085012 (2016)，MPS **数值**，非定理）。
  2. arXiv API: `au:"Dimock J" AND abs:gauge` → 命中 **6** 条，逐条看摘要，**无**谱隙定理（见 C6）。
  3. Crossref: `query.author=Balaban & bibliographic=lattice gauge renormalization ultraviolet stability` → 命中 8 条，全部为 UV 稳定性/RG 有效作用量（见 C2/C3）。
  4. Crossref: `query.author=Federbush & bibliographic=phase cell cluster expansion gauge mass gap` → 命中 8 条，全部体积/团簇展开（见 C5）。
  5. web_search: `"mass gap" "independent of the lattice spacing" Yang-Mills theorem`；`"uniform in the lattice spacing" spectral gap Yang-Mills Balaban` → 只出到两类来源：(i) Balaban 的 UV 稳定性被称为 "uniform in the lattice spacing"（**作用量**层面）；(ii) Zenodo/GitHub 上的**未发表**自印稿（"conditional proof"、"Phase Closure"、"YM Bridge Theorem" 等），**非同行评议，不引用**。
- 因此：**"`a`-一致谱隙定理"在我们的检索范围内不存在**；`Dimock–Hurd` 的格点规范 RG 分析在文献中被反复引作 "in preparation"/未发表（本次检索 Crossref 只核到二人**已发表**的红外 QED 与 KT 工作：*J. Math. Phys.* **33** (1992) 814–821，`10.1063/1.529761`；*CMP* **137** (1991) 263–287，`10.1007/BF02431881`）——**此点标注为"未核到正式出版记录"，不作定理引用**。

---

## 3. `Λ`（UV cutoff）与 `m` 的两种形式：文献实际写哪一种？

**文献实际写的是"无量纲、随 `a` 趋零"的那一种，即 `t3` 的 `ma ≥ c·aΛ` 同型。** 逐条说明：

1. **格点文献里不存在"`a`-无关的 `m ≥ c(g)`"。** 一旦 `m` 带量纲而 `a` 是唯一长度单位，所有界都以 `a^{-1}` 为单位给出，写成无量纲形式 `a m ≥ c·f(g)`，`f(g) ~ exp(-c/g²)`（弱耦合）。这正是 `t3` 的 `ma ≥ c·aΛ`（取 `Λ=1/a`，`f(g)=aΛ'`）。
2. **沿临界线 `a→0` 时 `f(g)→0`** ⇒ 界**随 `a` 趋零**，不给 `a`-一致正下界。物理上这是必然的：`a m_phys → 0`，任何"`a m ≥ c>0` 一致"的断言都会与 `m_phys` 有限矛盾。
3. **要得到"`m_{0^{++}} > 0`"，必须先证 `m ≍ Λ_QCD`（`Λ_QCD` 由 `a(β)` 的标度关系动力学生成）**，即 `m ≥ c·Λ_QCD`——这才是 `a`-一致的形式。它需要的正是 dimensional transmutation 的严格版本。出处（标度关系本身）：Wilson, *PRD* **10** (1974) 2445，DOI `10.1103/PhysRevD.10.2445` **[元数据]**。
4. **结论（对齐判死线 (c)）**：文献中可得的只有"随 `a` 趋零的界"。`m ≥ c·Λ·f(g)` 这种"`Λ` 挂在界里"的写法，在 `Λ = 1/a` 时**恒等于**无量纲 `a m ≥ c f(g)`；把它误读成"物理下界"就是 `a→0` 时的符号/量纲错位。**没有出处支持 `m ≥ c(g)`（`a`-无关）**——本次检索未查到任何一篇这样写的已发表文献（检索式同 §2）。

---

## 4. 一句话结论（供总册引用）

> 已发表数学的边界在 **UV 稳定性**：Balaban（`CMP` 102/109/116）、MRS（`CMP` 155）、Battle–Federbush（`Ann. Phys.` 142）分别在 `a`（有效作用量）、UV cutoff（固定 IR）、`L` 上给出了一致性；**"谱隙在 `a` 上一致"这一步没有任何定理**，`Dimock` 的 Balaban 型程序止于"第一步/不控制余项"。故 `A5` 的连续极限缺口属**领域边界 (b)**，其形式表现是 `t3` 的 `ma ≥ c·aΛ`（**界随 `a` 趋零，(c)**）。

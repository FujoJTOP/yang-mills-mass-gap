# LIT-CHECK — 文献核查 + 新颖性裁决：`(NA)` 有无定理、五机制占位、4D 最强半径、`F3∧F4` 可计算性对照

> 标签：**[R] 已核外部文献** / **[J] 判断（非定理）** / **[C] 条件** / **[NG] no-go** / **[X] 抓错更正**
> 纪律：本条任务**只**做核查，不做 QED。凡未查到证据的，一律写"未查到（检索式：…）"，**不补**。
> 全部检索于 2026-11-24 由 `web_search` / `web_fetch` 实做；每条给 URL 或 DOI/arXiv。**未读原文的，明写"未读原文"。**
> 上游依赖：`SIXTH-SPECIES-CHARTER.md` §1/§2/§4/§6、`docs/new-tools-registry.md`（F1 待审路线）、`research/unconditional/SPECIES-EXISTENCE.md` §2–§3、`docs/paper-v2/ch11b-related.tex` + `refs.tex`。

---

## §1 四个问题的判定表（本文件的正文核心）

| # | 问题 | 判定 | 一句话依据 | 锚点 |
|---|---|---|---|---|
| **Q1** | Haar×Wilson 型测度上，`Cov_μ(G,Φ′\|0) ≤ 0` 有无**已发表**负相依/负相关定理？ | **不存在** | 负相依整套机器（Pemantle / BBL / Joag-Dev–Proschan / Hargé）全部建立在**格偏序 / 多项式 / log-concave / 高斯**四类结构上；Haar×Wilson 一样都不是。**不是漏找**：Pemantle 本人说 negative dependence "not nearly as robust" | [Pemantle 2000](https://doi.org/10.1063/1.533200)；[BBL 2009](https://arxiv.org/abs/0707.2340) |
| **Q1′** | 该协方差的**一般符号定理**存在吗？ | **存在，但给出相反结论** | 分部积分恒等式把该协方差**精确等于期望 Hessian**：`Cov(∂_iS, ∂_jS) = E[Hess S]_{ij}` ⇒ 无普适符号律；而"`Hess S ⪰ ρ ⇒ gap ≥ ρ`"就是 **Bakry–Émery Γ₂**。YM 的 Savvidy Hessian **不定** ⇒ 判据失效（这正是 SZZ 必须退到小 β 的原因） | [Stein/Efron–Stein](https://doi.org/10.1214/aos/1176345513)；SZZ [arXiv:2204.12737](https://arxiv.org/abs/2204.12737) |
| **Q1″** | Weingarten / 6j 符号导致"正相依失效"，有无**定理**？ | **现象部分存在；定理不存在** | Weingarten 函数不定号（sign problem）是公认事实且有综述；**把它转成负相关/反 FKG 定理的论文：未查到** | [Collins–Matsumoto–Novak, Notices AMS 69(5) 734–745 (2022)](https://doi.org/10.1090/noti2474) |
| **Q2-1** | Harish–Chandra `c`-函数 / 非紧对称空间 Plancherel / 球变换，用于格点规范或谱测度原子？ | **不存在** | `c`-函数确实在非紧调和分析里是主力（Gindikin 50 年综述、Helgason），也出现在**别的**模型的谱权里，但**格点规范/YM 质量间隙语境：未查到** | [Gindikin, Ann. Fac. Sci. Toulouse 25(2-3) 385–402 (2016)](https://doi.org/10.5802/afst.1498) |
| **Q2-2** | 伪谱 / 非正规 / Kreiss 常数 / resolvent 谱包含，用于转移算子谱隙**下界**？ | **部分存在（非 YM）** | 伪谱 + Kreiss 矩阵定理是**瞬态增长**的经典工具；非可逆链有 "pseudo-spectral gap" `γ_ps`，但用于**混合时间**，且由可逆对称化 `γ((M*)^kM^k)/k` 构造。**YM / Gribov 转移算子：未查到** | [Trefethen–Embree 2005](https://press.princeton.edu/books/hardcover/9780691119427/spectra-and-pseudospectra)；[Paulin, EJP 20 (2015), arXiv:1212.2015](https://arxiv.org/abs/1212.2015) |
| **Q2-3** | max-plus / tropical / Maslov 去量子化，用于统计力学转移矩阵的**严格界**？ | **部分存在（非 gap）** | 有成体系的 **tropical 热力学形式化**（Bousch 算子 = Ruelle 算子的 tropical 对应）与 Maslov 去量子化；但都是 `β→∞ / h→0` 的**极限陈述**，给的是**领头指数**，对"第二特征值间隙"没有 tropical 对应物。**严格 gap 界：未查到** | [Li–Sun, arXiv:2408.10169](https://arxiv.org/abs/2408.10169)；[Litvinov, arXiv:1005.1247](https://arxiv.org/abs/1005.1247) |
| **Q2-4** | Lorentzian 多项式 / M-凸 / stretched LR 对数凹，用于格点规范特征标展开？ | **不存在** | Lorentzian 理论住在 `{0,1}`-支撑 / M-凸支撑上；Wilson 特征标展开系数**变号**（Weingarten 符号）、底空间是紧群对偶、不是 M-凸集。**并且** stretched LR 对数凹本身在一般情形**仍未证** | [Brändén–Huh, Ann. of Math. 192(3) 821–891 (2020)](https://doi.org/10.4007/annals.2020.192.3.4)；[arXiv:2607.22301 (2026)](https://arxiv.org/abs/2607.22301) |
| **Q2-5** | 小球概率 / 倍增 / Nazarov–Turán，用于把测度方差与原子（谱下确界）挂钩？ | **不存在** | Nazarov–Turán 是**局部**指数多项式估计（不确定性/唯一延拓）；小球概率是高斯过程的 `L₂`-小球渐近。**方向也不对**：小球给"质量不能集中"的**上界**，而 gap 要"质量确实集中在 `m₁²>0`"——符号反向 | [Nazarov, Algebra i Analiz 5 (1993)](https://www.mathnet.ru/eng/aa416) |
| **Q3** | 2023–2026 是否已有 4D 纯 YM 的**非微扰、(L,a)-一致、控号下界**？ | **不存在** | 最强仍是**强耦合 / 't Hooft regime**：SZZ `\|β\|<1/48`（d=4, SU(N)）；CNS 2025 只做**面积律**且自认 d=4 大 λ / 小 N 开；Nissim 2025 仍是 't Hooft regime；Chatterjee 2019 强耦合大 N；Balaban 只到 UV。**无一条 (L,a)-一致、无一条过 β≈2.5** | SZZ [arXiv:2204.12737](https://arxiv.org/abs/2204.12737)；CNS [arXiv:2505.16585](https://arxiv.org/abs/2505.16585) |
| **Q3′** | 数字核对：`β<1/48` ✔ / `λ_tH<1/24` ✘ | **`1/24` 在源文中不存在** | SZZ 条件是 `\|β\|<1/(16(d−1))`，d=4 ⇒ `1/48`；且该条件**与 N 无关**，写成 't Hooft 耦合是 `λ=βN < N/48`，**不是** `1/24`。语料自己的 `SOTA-VERIFIED.md` 已标 UNVERIFIED，本文件**独立复核确认** | [SZZ abstract](https://arxiv.org/abs/2204.12737) |
| **Q4** | 元定理 `F3∧F4=∅` 的可计算性层面，文献有等价结果吗？ | **部分存在（且反向修正一处）** | (a) `T↦σ(T)` 只有**下半可计算**、`σ(T)` 是 r.e.-紧但**非递归紧**；任一 r.e.-紧集都是某可计算正规算子的谱 ⇒ "谱可计算"一般**假**。(b) **但**紧自伴算子的谱**完全可计算** ⇒ "原子不可计算"作为**无条件断言是假的**，必须改成"**族上一致**不可计算"。(c) 有效 Tauberian：Littlewood Tauberian 的**有穷化**已发表，且明说存在"**没有可计算收敛速率**"的情形 —— 这正是 B6 的证明论对应物 | [Brattka–Dillhage, J. UCS 11(12) 1884–1900 (2005)](https://zbmath.org/?q=an:1120.47063)；(c) [APAL, S0168007222001464](https://doi.org/10.1016/j.apal.2022.103117) |

---

## §2 逐条证据（含"我核到的是…"更正记录）

### §2.1 Q1 — 负相依/负相关：机器在，但不是这台

**(a) 已发表的负相依理论（四类结构，全部不覆盖 Haar×Wilson）。**

1. **Pemantle, "Towards a theory of negative dependence"**, *J. Math. Phys.* **41**(3) 1371–1390 (2000), DOI [10.1063/1.533200](https://doi.org/10.1063/1.533200), arXiv:[math/0404095](https://arxiv.org/abs/math/0404095)。原文（检索到的正文片段）明写：*"Negative dependence, by contrast, is not nearly as robust… since a random variable is always positively correlated with itself, one cannot expect all monotone functions to be negatively correlated."* 定义域是 `B_n = {0,1}^n`（布尔格）。
2. **Borcea–Brändén–Liggett, "Negative dependence and the geometry of polynomials"**, *J. Amer. Math. Soc.* **22** (2009) 521–567, arXiv:[0707.2340](https://arxiv.org/abs/0707.2340)。核心定理链：**strongly Rayleigh ⇒ NA ⇒ pairwise negative correlation**；承载物是**多项式**（实稳定 / half-plane property）与**行列式测度**、生成树测度、对称排空过程。
3. **Joag-Dev & Proschan, "Negative association of random variables with applications"**, *Ann. Statist.* **11** (1983) 286–295。连续情形**唯一**的一般判据是高斯：`NA ⟺ 两两 `cov ≤ 0``（见 zbMATH 转述：[zbmath 1218.62059](https://zbmath.org/?q=ci%3A0482.62046)）。
4. **Hargé (2004)**, 凸/对数凹相关不等式（偶凸 与 偶对数凹 **负相关**），*Probab. Theory Related Fields*；转述见 [Bulletin UOC 64](https://bulletin.math.uoc.gr/vol/64/64-1-31.pdf) Theorem 19 与 [ALEA 21-25](https://alea.impa.br/articles/v21/21-25.pdf) Theorem 10.1。**卷页未核到**（只核到年份+刊名）。
5. **Scott–Sokal** 完全单调性/排斥格气：见 [Pemantle review](https://www2.math.upenn.edu/~pemantle/neg-dep-midwest221021.pdf) 的 SR 层级图与 [MPP 讲义](https://indico.mpp.mpg.de/event/10191/contributions/38073/attachments/20667/26201/17.10.pdf) 的"Szegő ⇒ Scott–Sokal"。
6. **Brascamp–Lieb**（方差上界，非负相关）：*J. Funct. Anal.* **22** (1976) 366–389。

**[X] 更正 1（抓错）**：`research/unconditional/SPECIES-EXISTENCE.md` 参考文献列表写 "Brascamp–Lieb, *Adv. Math.* **22** (1976) 366"。**我核到的是** *J. Funct. Anal.* **22** (1976) 366–389（Brunn–Minkowski/Prékopa–Leindler 那篇）；`docs/paper-v2/refs.tex:51` 的写法是对的，SPECIES-EXISTENCE 那条需改。

**[X] 更正 2（抓错）**：任务书提到的 "Borodin–Diaconis–Fulman negative dependence zoo"。**我核到的是** BDF 名下可核论文为 *"On adding a list of numbers (and other one-dependent determinantal processes)"*, arXiv:[0904.3740](https://arxiv.org/abs/0904.3740)（zbMATH 引文片段转述）。**题名含 "zoo of negative dependence" 的 BDF 论文：未查到**（检索式：`Borodin Diaconis Fulman "zoo of negative dependence"`）。"zoo / 层级图"的通行载体是 Pemantle 的综述讲义，不是 BDF 论文。

**(b) 结论：不存在。** 四条理由都是结构性的，不是检索运气：

- **无格结构**。NA 的定义（`f,g` 在**不交坐标集**上单调）需要坐标偏序；`SU(N)` 上 Haar 测度的底空间没有这种偏序。这一条与语料 `C8-SPECIES-NOGO` 的"unstatable"独立吻合（内部 [R]）。
- **非 log-concave**。`μ ∝ e^{−βΣ tr(1−U_p)}` 的对数密度在紧群上**周期**，周期凹=常数（语料 B1 周期 no-go）；Hargé/BL 的整套假设（log-concavity / `Hess ⪰ 0`）全部落空。
- **NA 不封闭**。条件化、边缘化、混合都破坏 NA；`Cov_μ(G,Φ′|0)` 恰恰是**条件化**后的量，即便整体有 NA 也搬不过来。
- **变号是本质的**。Wilson 特征标展开的 Weingarten 系数不定号（见 (c)），任何"正相依"或"负相依"的**单调性**论证都必须先有一个符号确定的权重。

**(c) Weingarten / 6j 符号：现象已发表，定理未发表。**
- Weingarten 计算的权威综述：**Collins–Matsumoto–Novak, "The Weingarten calculus"**, *Notices Amer. Math. Soc.* **69**(5) 734–745 (2022), DOI [10.1090/noti2474](https://doi.org/10.1090/noti2474)；原始积分公式 **Collins–Śniady**, *Comm. Math. Phys.* **264** (2006) 773–795, DOI [10.1007/s00220-006-1554-3](https://doi.org/10.1007/s00220-006-1554-3)；中心化/符号尊重公式的新进展 **Collins–Nagatsu**, *Algebraic Combinatorics* **9**(2) 465–479 (2026), DOI [10.5802/alco.480](https://doi.org/10.5802/alco.480)。
- **"正相依失效 ⇒ 负相依成立"这一步，未查到任何已发表定理**（检索式：`Weingarten 6j symbol negative correlation lattice gauge`, `"negative dependence" fails SU(N) lattice gauge sign problem`, `anti-FKG lattice gauge`）。判定：**部分存在（现象）/ 不存在（定理）**。
- [J] 这不是缺口而是**逻辑事实**：不定号只否定 FKG 的**前提**，不给出相反的**结论**。要得到 `≤0` 必须另造不等式——与章程 §5 的 (α)/(β) 分类一致。

**(d) Stein / score–Hessian：一般定理存在，答案与 (NA) 相反方向。**
对 `μ ∝ e^{−S}`，分部积分（Stein 恒等式）给

```
E_μ[Φ · ∂_i S] = E_μ[∂_i Φ]        ⇒   Cov_μ(Φ, ∂_iS) = E_μ[∂_iΦ]
取 Φ = ∂_j S:   Cov_μ(∂_iS, ∂_jS) = E_μ[Hess S]_{ij}
```

来源：Stein (1981) *Ann. Statist.* **9** 1135–1151；Efron–Stein *Ann. Statist.* **9** (1981) 586（`refs.tex:88` 已列）。

**这三个推论直接决定 `(NA)` 的命运：**
1. **`(NA)` 不是"负相依"，而是一条 Hessian 符号陈述**。`Cov(G,Φ′|0)` 经上述恒等式等价于某个**条件期望 Hessian**；要它 `≤0`，必须证明该 Hessian 的某个方向**负定**。而语料 A2 的恒等式 `G″ = E[h] − Var(f) + W_geo` 说明 `G″` 是**两个正项减一个**——`(NA)` 正是"`Var` 项赢"的陈述。**这正是 H1，不是 H1 的新物种。**
2. **唯一已发表的一般"符号定理"是 Bakry–Émery Γ₂**：`Hess S ⪰ ρ > 0 ⇒ gap ≥ ρ`（Bakry–Émery 1985, *Sém. Probab. XIX*, LNM **1123**, 177–206）。它是**充分**条件，方向是 `Hess ⪰ 0` ⇒ gap。YM 的 Savvidy Hessian 是 `−D² + 2ad F_B`，**不定**（负模横向磁），所以判据**按定义**失效。SZZ 明确用 Bakry–Émery，且正因为 U(N) 的 Ricci 曲率**不一致正**而失效，才要另造随机环境分解（[arXiv:2510.22788](https://arxiv.org/abs/2510.22788) 摘要原文：*"this approach fails because its Ricci curvature is not uniformly positive"*）。
3. **因此 `(NA)` 在文献里的"近亲"恰恰是它想绕开的东西**：要么 `Hess ⪰ 0`（假），要么 `Var − E[h] ≤ 0`（= H1）。这不是文献解决的，是文献**证明这条路必须回到 H1** 的旁证。

---

### §2.2 Q2 — 五机制的占位情况（逐条）

| 机制 | 占位判定 | 已发表的用法 | 与本战役 gap 语境的差 |
|---|---|---|---|
| 1. Harish–Chandra `c`-函数 / 非紧 Plancherel / 球变换 | **未被占用（不存在）** | 非紧对称空间谱理论主力：[Gindikin 2016](https://doi.org/10.5802/afst.1498)；Helgason 的 `c`-函数 ↔ Plancherel 密度；在 **SO₀(L;N)** 型模型的"谱权 = `c`-函数"里出现过（[math-ph/0601049](https://arxiv.org/abs/math-ph/0601049)） | YM 的底空间是**紧**群上的 Wilson 测度，不是非紧 `G/K`；`c`-函数编码的是**沿 Weyl 腔的渐近**，而我们要的是**紧群上的离散谱原子**。**无桥**（[J]） |
| 2. 伪谱 / 非正规 / Kreiss / resolvent 谱包含 | **部分占用，非 gap** | 伪谱+Kreiss 矩阵定理 = 瞬态增长经典（[Trefethen–Embree 2005](https://press.princeton.edu/books/hardcover/9780691119427/spectra-and-pseudospectra)）；非可逆链的 `γ_ps = max_k γ((M*)^kM^k)/k` 用于**混合时间**（[Paulin, EJP 20 (2015)](https://arxiv.org/abs/1212.2015)） | `γ_ps` 是 `1/t_mix` 的**上界**工具（≡ 谱隙的**下界**只在可逆对称化后成立），且**逐链**、非 `(L,a)` 一致；YM/Gribov 转移算子非正规性来自符号问题，`γ_ps` 的对称化恰好把符号问题**抹平**（= 换了测度）。**未查到 YM 用法** |
| 3. max-plus / tropical / Maslov 去量子化 | **部分占用，非 gap** | [Li–Sun, "Tropical thermodynamic formalism", arXiv:2408.10169](https://arxiv.org/abs/2408.10169)：Bousch 算子 `𝓛_A` = Ruelle 算子的 tropical 对应；[Litvinov, arXiv:1005.1247](https://arxiv.org/abs/1005.1247) 的 Maslov 去量子化；BBS/超离散化 | 去量子化是 `β→∞`（或 `h→0`）的**极限**，给 `log Z` 的**领头指数**（= 基态能量 / 最大特征值），用 Collatz–Wielandt 得**单侧**界。**tropical 谱论没有稳定的"第二特征值"概念** ⇒ 对**间隙**结构性无输出（[J]，但这是硬结构事实） |
| 4. Lorentzian / M-凸 / stretched LR 对数凹 | **未被占用（不存在）** | [Brändén–Huh, Ann. of Math. 192(3) 821–891 (2020)](https://doi.org/10.4007/annals.2020.192.3.4)：Lorentzian ⇒ 支撑是 M-凸集、系数对数凹、"Independence and negative dependence" 一节 | Wilson 特征标展开的系数**变号**、指标集是 `SU(N)` 的表示环（非 M-凸集），且展开**不是**多项式而是无穷级数。**并且** stretched LR 对数凹一般情形**仍未证**（[arXiv:2607.22301 (2026)](https://arxiv.org/abs/2607.22301) 明列为 FrontierMath 未解问题）⇒ 连"借用"的原料都没有 |
| 5. 小球概率 / 倍增 / Nazarov–Turán | **未被占用（不存在）** | Nazarov–Turán：局部指数多项式估计，用于不确定性/唯一延拓（[Nazarov 1993](https://www.mathnet.ru/eng/aa416)，亦见 [Tao 讲义 Lemma 13](https://math.berkeley.edu/~ztao/paper2.pdf)）；小球概率：高斯过程 `L₂`-小球的**精确渐近**（[Nazarov–Nikitin 综述](http://www.pdmi.ras.ru/EIMI/2005/sd/talk/nazarov.pdf)） | **方向反**：小球给"质量不能压进小集合"的**上界**，即**下界**在展布上；gap 要的是"质量确实压在 `m₁²>0`"。要让两者相遇必须先把 `Var` 换成 `inf supp`，而语料的 squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` 已证明 `Var`（=`χ`）只给 `m₁²` 的**上界**。**机制 5 撞的是同一条 squeeze** |

**Q2 汇总判决**：五条里 **0 条被 gap 语境占用**；**2 条（2、3）在别的语境被占用且其占用手法恰好抹掉我们要的信号**；**3 条（1、4、5）不但未被占用，其原料还缺（4 的 stretched LR 本身未证）或方向反（5）**。

---

### §2.3 Q3 — 4D 纯 YM 的最强已知半径（截止 2026-11）

#### 已核到的、**同行评议或 arXiv 正刊**的严谨结果（按强度排序）

| 结果 | 精确陈述 | 是否 (L,a)-一致 | 是否过 β≈2.5 | 引用 |
|---|---|---|---|---|
| **SZZ 2023** | `d>1`、't Hooft 标度 `βN`；**`\|β\| < (N−2)/(32(d−1)N)`（SO(N)）或 `\|β\| < 1/(16(d−1))`（SU(N)）** 下：Langevin 不变测度唯一、指数遍历、**无穷体积极限唯一**、**LSI + Poincaré** 成立、关联指数衰减（**严格正质量间隙**） | **体积一致 ✔；格距不一致 ✘**（间隙在格点单位，`a→0` 时 `m_phys=m_lat/a` 反而增大，不是连续极限） | **✘**（强耦合，`d=4` ⇒ `β<1/48`） | *Comm. Math. Phys.* **400** (2023) 805；arXiv:[2204.12737](https://arxiv.org/abs/2204.12737)；DOI [10.1007/s00220-022-04609-1](https://doi.org/10.1007/s00220-022-04609-1)；Zbl 1522.81127 |
| **CNS 2025a（面积律扩张）** | "extend the parameter regimes for which **area law** is proven for pure `U(N)` … in particular when N is large"；改进 Osterwalder–Seiler 1978；方法：master loop equation 当线性非齐次方程 + a priori 界 + truncated model | 未声称 `(L,a)` 一致 | **✘**（这是**面积律**，不是间隙；且自认 d=4 大 λ / 小 N 开） | arXiv:[2505.16585](https://arxiv.org/abs/2505.16585)（v3, 2025-09-26） |
| **CNS 2025b（动力学面积律）** | 把 SZZ 的动力學方法用到**'t Hooft regime** 证 Wilson 面积律，验证 [DF80] 的 mass gap 条件；`G ∈ {U(N), SU(N), SO(2N)}` | 未声称 | **✘** | arXiv:[2509.04688](https://arxiv.org/abs/2509.04688)（v2, 2025-09-28） |
| **Nissim 2025** | `U(N)` 在 't Hooft regime 的**质量间隙 + 唯一无穷体积极限 + 大 N 极限新证明**；技巧：`U(N)` 重写为随机环境 `SU(N)`（`U(1)` 场），cluster expansion + Langevin | 体积一致；格距不一致 | **✘**（'t Hooft regime；原文：U(N) 的 **Ricci 曲率不一致正**，Bakry–Émery 直接失效） | arXiv:[2510.22788](https://arxiv.org/abs/2510.22788) |
| **Chatterjee 2019** | 强耦合 `SO(N)` 大 `N` 的 Wilson 环**严格解** | — | **✘** | *Comm. Math. Phys.* **366** 203–268 (2019)，DOI [10.1007/s00220-019-03353-3](https://doi.org/10.1007/s00220-019-03353-3)；arXiv:[1502.07719](https://arxiv.org/abs/1502.07719) |
| **Osterwalder–Seiler 1978** | 格点规范场论的构造 + 强耦合质量间隙（P3 的技术源头） | — | **✘** | *Ann. Phys.* **110** 440–471 (1978)，DOI [10.1016/0003-4916(78)90039-8](https://doi.org/10.1016/0003-4916(78)90039-8) |
| **Balaban 系列** | 格点规范 RG 的**小场**分析、大场重整化；只到 **UV**，IR 未闭合 | — | **✘** | *Comm. Math. Phys.* **109** 249–301 (1987)，DOI [10.1007/BF01215223](https://doi.org/10.1007/BF01215223)（经 Crossref 参考文献表核对） |
| **Magnen–Rivasseau–Sénéor 1993** | 纯 SU(2) `YM₄` 在**固定红外截断**下的构造（无 UV 截断） | ✘ | ✘ | *Comm. Math. Phys.* **155** 325–384 (1993)，DOI [10.1007/BF02097397](https://doi.org/10.1007/BF02097397) |
| **Chandra–Chevyrev–Hairer–Shen 2024** | 3D `Yang–Mills–Higgs` 随机量化 | 3D | ✘ | *Invent. Math.* **237** 541–696 (2024)，DOI [10.1007/s00222-024-01264-2](https://doi.org/10.1007/s00222-024-01264-2) |
| **Cao–Chatterjee 2024** | 3D 欧氏 YM 的**状态空间** | 3D | ✘ | *Comm. Math. Phys.* **405** 3 (2024)，DOI [10.1007/s00220-023-04870-y](https://doi.org/10.1007/s00220-023-04870-y) |

#### Douglas 2026 综述（Nature Rev. Phys.）
- **已核到的书目记录**（Crossref 全记录）：Michael R. Douglas, **"The Yang–Mills Millennium problem"**, *Nature Reviews Physics* **8**(2), **86–97** (2026)；DOI [10.1038/s42254-025-00909-2](https://doi.org/10.1038/s42254-025-00909-2)；online **2026-01-12**，accepted **2025-11-14**；单一作者，Harvard；74 条参考文献（含 SZZ、CNS、Chatterjee、Balaban、MRS、Meurice 的 tensor RG 综述）。
- **未读原文**：`nature.com` 正文页与 PDF 均被 `idp.nature.com` 重定向拦下（`web_fetch` 报 cross-origin redirect），我**没有**读到它的论断原文。因此**本文件不引用它的结论句**；`new-tools-registry.md` F1 里"明说 no clear path known"这句，我**无法独立核实**，标 [R-未核原文]。二手转述（[Meurice 报告摘要](https://math.uiowa.edu/event/36209/0)）引 Douglas 说 *"Conceivably, the tensor RG [that we developed] will be adequate…"* —— 这是**乐观**语气，与"no clear path"方向相反。**[X] 更正 3**：引用 Douglas 时应改写为"综述梳理了包括 tensor RG 在内的路线"，不要写"明说无路"。

#### 数字核对（任务的"例如 `β<1/48`、`λ_tH<1/24`"）

- **`β < 1/48`：✔ 正确，来源确定。** SZZ 摘要原文 `|β| < 1/(16(d−1))` for `SU(N)`；`d=4` ⇒ `1/(16·3) = 1/48`。语料 P3 的 `β<1/48` 与 SZZ 阈值**同源**（这一点应在 P3 文档里显式写明，否则读者会以为是独立巧合）。
  - **重要限定**：SZZ 的 `β` 是"`βN` 作逆耦合"约定下的量，条件**与 N 无关**。若换成 't Hooft 耦合 `λ = βN`，则条件是 `|λ| < N/48`，**不是常数**。
- **`λ_tH < 1/24`：✘ 源文中不存在。** 我把三篇（SZZ / CNS-2505.16585 / Nissim-2510.22788）的摘要与可得正文片段全查过，**没有任何一篇出现 `1/24` 或 `1/(8(d−1))`**（检索式：`"1/24" 't Hooft coupling lattice Yang-Mills mass gap threshold theorem`）。
  - 语料 `research/unconditional/SOTA-VERIFIED.md:27` 已自行标注 **UNVERIFIED**，`TRANSPORT-RIGOROUS-NUMERICS.md:89` 进一步说 CNS 的真实阈值是 `β ≤ β₀(d) = 10^{-10d}/d`（均匀于 N、未归一化迹）。**我无法独立复核该式**（CNS PDF 不能经 `web_fetch` 取文本，`application/pdf` 不支持），故该式标 **[R-未核原文]**；但"`1/24` 不在源文中"这一否定结论我**独立确认**。
  - **必须清理的位置**（我在仓库里 grep 到 8 处仍在用 `1/24`）：
    - `docs/paper-v2/ch07-nogos.tex:235`
    - `docs/paper-v2/ch08-immune.tex:120`
    - `docs/paper-v2/ch13-machine.tex:110, 123`
    - `docs/paper-v2/app.tex:563`
    - `docs/paper/body3.tex:252, 395`；`docs/paper/body4.tex:73`
    - `docs/YANG-MILLS-MASTER.md:100, 145`；`docs/YANG-MILLS-CONDITIONAL-REDUCTION.md:518, 578, 635`；`docs/v2.txt:2511, 2640, 4004, 4014, 5287`
    - `research/unconditional/BROWNIAN-LOOP-SOUP.md:9, 30`；`PROOF-CONFINEMENT-RANDOM-SURFACE.md:7, 13`；`PROOF-MULTISCALE-4D.md:3, 9`；`research/lit-v3/LIT-ERA-03.md:181`
    → 建议统一改为 `β<1/48（SZZ, d=4, SU(N)）`，并**删掉"~75×"**（该因子由 `1/24` 推出，随之失效）。

#### 灰文献噪声（必须记录，但**不计入**证据）
2025–2026 年出现大量 Zenodo / PhilArchive / vixra 的"4D YM 质量间隙证明"（例：[zenodo 20114011](https://zenodo.org/records/20114011)、[zenodo 20300895](https://zenodo.org/records/20300895)、[zenodo 18210079](https://zenodo.org/records/18210079)、[philarchive NIETEO-12v3](https://philarchive.org/archive/NIETEO-12v3)）。**全部未经同行评议**，按 deep-research IRON RULE（Tier 3，不可核 = FAIL）**不作为文献证据**。这是 Q3 判定"不存在"的一部分理由：**不是没人声称，是没人证成**。

#### Q3 判定
**不存在**。截止 2026-11，**没有任何** 4D 纯 YM 的"非微扰 + `(L,a)` 一致 + 控号下界"结果。最强者是 SZZ 的 `β<1/48`（体积一致、强耦合、格点单位），**离 β≈2.5 的物理交叉点差 ~2 个数量级在 β 上，且完全不给 `a→0`**。语料 `CLAUDE.md §5` 的"这台机器尚不存在"在文献侧**成立**。

---

### §2.4 Q4 — `F3∧F4=∅` 的可计算性对照

**(a) 已发表且**精确对应**：谱/原子的可计算性程度。**
- **Brattka & Dillhage, "Computability of the spectrum of self-adjoint operators"**, *J. UCS* **11**(12) 1884–1900 (2005), Zbl [1120.47063](https://zbmath.org/?q=an:1120.47063)。主结果（原文片段）：
  1. (Uniform) `T ↦ σ(T)` **下半可计算**（只有正信息）；
  2. (Non-uniform) 任一可计算 `T` 的 `σ(T)` 是 **r.e.-紧**集；
  3. **不能**强化为**递归紧**：**任一 r.e.-紧集都可实现为某可计算正规算子的谱**；
  4. **但**：`T` 若**紧**自伴，则谱**完全可计算**（距离函数可计算）。
- 后续：Brattka & Dillhage (2008), *ENTCS*, "…and the Computable Operational Calculus" —— 紧自伴算子的**特征值可计算**（[ACM DL 摘要](https://dl.acm.org/doi/abs/10.1016/j.entcs.2008.03.026)）。
- **Pour-El & Richards**：*Computability in Analysis and Physics* (Springer, 1989)；Second Main Theorem（可计算自伴算子其谱不可计算）——Brattka–Dillhage 明言其非一致版本"also follows from the Second Main Theorem of Pour-El and Richards"。

**(b) 这直接修正语料的一处表述。** [X] **更正 4**：
`SPECIES-EXISTENCE.md` §5 与 `C8-SPECIES-NOGO` 若把 `F3∧F4=∅` 的根据写成"**原子不可计算**"，则**作为无条件断言是假的**：对**紧**自伴算子（固定 cutoff 下的转移算子 `T_O = e^{−aH}` 恰是紧的、且数据可计算），谱**完全可计算**。**真实的障碍不是"不可计算"，而是"不可一致（uniform）"**——即
```
   {T_{L,a}} 的谱没有可从 (L,a,数据) 一致计算的表示
```
这才是 `F3∧F4=∅` 在可计算性层面唯一站得住的形态，且它是**族上**陈述，恰好对应章程 §2 的"**均匀**"要求。这条修正**加强**而非削弱元定理：它把"computable ⇒ floor-blind"从一条可疑的一般断言，收紧成一条**关于族的一致可计算性**的精确断言。

**(c) 侧证：原子在"点"层面的真实程度。**
- [arXiv:1907.07815](https://arxiv.org/abs/1907.07815)（"Degrees of Randomized Computability"）：*"If `P` is a left-c.e. semi-measure and `X` is an atom of `P`, then `X` is computable."* ⇒ 原子**位置**往往可计算，难的是**质量**与**族上一致**。
- Hoyrup–Rojas：原子less 可计算概率空间的同构定理（转述见 [arXiv:1203.5535](https://arxiv.org/abs/1203.5535) Thm 7.16）；原子定义 `μ({X})>0` 见 [PSU 论文](https://etda.libraries.psu.edu/files/final_submissions/21992)。
- **"atom of a measure 的 Weihrauch 度"：未查到**（检索式：`Weihrauch degree "atom of a measure"`, `"atoms of measures are not computable"`）。最接近的是上两条 + Brattka–Dillhage 的谱度。

**(d) 有效 Tauberian：存在，且正是 B6 的证明论对应物。**
- **"A finitization of Littlewood's Tauberian theorem and an application in Tauberian remainder theory"**, *Ann. Pure Appl. Logic*（[DOI 10.1016/j.apal.2022.103117](https://doi.org/10.1016/j.apal.2022.103117)）。原文片段：*"…the case where **no computable rate of convergence exists**"* —— 即 Tauberian 余项**可以没有可计算速率**。
- **"A note on the finitization of Abelian and Tauberian theorems"**, arXiv:[1910.07926](https://arxiv.org/abs/1910.07926)。
- 判定：**这部分 = 语料 B6（全单调 / 解析延拓 no-go，`t→∞` 速率）的已发表证明论版本**。B6 说"唯一能印 floor 的是 `t→∞` 速率 = 循环"；文献说"该速率**可以没有可计算界**"。两者讲同一件事的两个面。**这是 Q4 里"已被占用"最实的一条**：B6 的 **no-go 强度**已在文献中有独立支撑，但**没有人把它接到 KL 测度/质量间隙上**（那一步仍是本战役的）。

**(e) 一般原理层。** 可计算分析的标准事实：**可计算 ⇒ 在表示上连续**（Weihrauch, *Computable Analysis: An Introduction*, Springer 2000；[Brattka 的 discontinuity 层级讲义](https://ar5iv.labs.arxiv.org/html/cs/0607114)）。若"检出 `inf supp ν > 0`"在所用表示下**不连续**，则任何 `F4`-函数**必不能**检出。**这是 `F3∧F4=∅` 的抽象形式**，而它是否成立**完全取决于 `ν ↦ inf supp ν` 的不连续点集**——这正是 §3 的酸测 2（薄原子 `ν=wδ_ε+(1−w)δ_S`，`ε→0`）在做的事。**文献给了一般原理，没给这个实例**。

#### Q4 判定
**部分存在**。(i) 谱/原子的**可计算性程度**有精确已发表定理；(ii) **无条件"原子不可计算"是错的**，正确形态是"**族上不可一致计算**"；(iii) Tauberian 速率不可计算有已发表证明论结果（= B6 的镜像）；(iv) 把 (i)–(iii) 接成 `F3∧F4=∅` 这一步**未查到**，仍属本战役所有。

---

## §3 按章程七节骨架重排（本节为交付骨架，§1/§4/§6 已按任务书改写）

### §3.1 = 任务 §1（四问判定表）
见上文 **§1 判定表**（12 行）。核心三句：
1. `(NA)`：**无已发表定理**；其一般近亲（Stein 恒等式 / Bakry–Émery Γ₂）给出的是**相反方向的充分条件**，YM 按定义不满足。
2. 五机制：**0/5 被 gap 语境占用**；机制 2、3 在别处被占用且占用手法抹掉信号；机制 4 原料本身未证；机制 5 方向反。
3. 4D 最强：**SZZ `β<1/48`**（体积一致、强耦合、非 `(L,a)` 一致）；`λ_tH<1/24` **在源文中不存在**。

### §3.2 = 控号定理（本文件能证的 vs 缺的）
- **[R] 已发表控号定理（可引，不需重证）**：SZZ（`β<1/48` 的 LSI+Poincaré+指数簇聚）；Bakry–Émery Γ₂；Brascamp–Lieb；Osterwalder–Seiler 强耦合。
- **[T-here] 本文件能当场证的两条平凡但有用的推论**（一行）：
  1. **`(NA)` ⟺ Hessian 符号陈述**：由 Stein 恒等式，`Cov_μ(∂_iS,∂_jS) = E_μ[Hess S]_{ij}`。故 `(NA)` 若成立，必是"某方向期望 Hessian ≤ 0"；反之 `Hess ⪰ ρ>0`（Bakry–Émery）给出 `≥0`。**两方向不可能同时无条件下成立** ⇒ `(NA)` 与 `Hess ⪰ 0` **互斥**，而后者是 YM 的**假**命题。**这不是 no-go**（`(NA)` 只要求**条件化后某一特定方向**，不要求全域），但它把 `(NA)` 精确钉成"Savvidy 负模方向上的条件 Hessian 不等式"= 语料已经知道的那个对象。
  2. **`1/24` 不可救**：即便把 CNS 的阈值取到最优读法，`λ_tH = βN` 与 `β<1/48` 的换算给 `λ_tH < N/48`（N-依赖），**任何常数阈值都写不成 `1/24`**。⇒ 该数字**必须删除**，不是"待核"。
- **[O] 缺的步**：`(NA)` 的条件化 Hessian 不等式（= 语料 A2 的 `D<S`）；五机制的桥；族上一致可计算的表示；4D 过 β≈2.5 的 `a→0`。

### §3.3 = 酸测五点记分（针对"用文献结果能否过酸测"）
| # | 测试 | 若拿文献结果来顶，结果 | 理由 |
|---|---|---|---|
| 1 | 自由场失效 | **PASS（但空洞）** | Bakry–Émery / BL 在自由场 `Hess=0` 上给 `ρ=0` ⇒ 报零，形式上过酸测；但正因为 `ρ=0` 它也**证不出**任何 YM 的 `ρ>0` ⇒ 空洞 |
| 2 | 薄原子 | **FAIL** | BL/`γ_ps` 全是**方差/二阶**工具，`ν=wδ_ε+(1−w)δ_S` 上 `Var` 有界而 `m₁²=ε→0`（语料 squeeze 已证）⇒ 会给薄原子发证书 |
| 3 | 深双井 | **FAIL** | `γ_ps` 与 BL 都无 barrier-uniform 常数 |
| 4 | 归一化/真空性 | **N/A** | 文献工具作用在测度上，不涉及 `Var(q)` vs `Var(q̄)` |
| 5 | 量纲 | **FAIL** | 全部只给格点单位，无 `a` 依赖 |

### §3.4 = **新颖性裁决**（任务 §4）
> 裁决口径：**"新颖"= 本战役声称的新意，在文献里未被占用**。分四档：**已占用 / 部分占用 / 未被占用但原料缺 / 未被占用且原料齐**。

| 战役声称 | 裁决 | 说明 |
|---|---|---|
| **(NA) 是"variance-free 面"，是唯一的非正性逃逸口** | **未被占用** | 文献无任何 `Haar×Wilson` 负相依定理；但 **Stein 恒等式把它降级为 Hessian 不等式**，所以它**不是新物种**，是 H1 的一个坐标。裁决：**新颖 = 无（已被自己的 A2 恒等式吸收）** |
| **五机制（Harish-Chandra / 伪谱 / tropical / Lorentzian / 小球）** | 1、4、5 **未被占用**；2、3 **部分占用** | 但 1 无桥、4 原料缺、5 方向反 ⇒ **可探索性 ≠ 可用性**。裁决：机制 4 **最不该走**（原料本身是未解问题）；机制 1 **最值得再想**（唯一"零占位 + 零方向冲突"） |
| **`F3∧F4=∅` 是元定理级 no-go** | **部分占用，且需改口径** | 一般原理（可计算⇒连续）与谱度（Brattka–Dillhage）已发表；**但**"原子不可计算"作为无条件断言**假**；正确形态"族上不可一致计算"**未查到** ⇒ **这一形态是本战役可主张的新意，但必须改写现有表述** |
| **B6（全单调/Tauberian 速率不可计算）** | **已占用（证明论侧）** | APAL 2022 的 Tauberian 有穷化明说"无可计算速率"。裁决：B6 作为 **no-go 有独立文献支撑**（这增强它），但**不能**当成新贡献主张 |
| **P3 `β<1/48`** | **已占用** | 阈值与 SZZ 完全同源（`1/(16(d-1))`, d=4）。裁决：P3 的**新意只能在残数 `f₁²` 的一致性下界**，**不能**在 `1/48` 这个数上。**必须显式引用 SZZ**，否则是重复发现 |
| **SZZ/CNS 只到强耦合、`1/24`** | **`1/24` 是自造** | 必须删。裁决：不是新颖性问题是**错误**问题 |
| **"无均匀可计算证书、一般证书 = gap 本身"（Theorem 3 等价）** | **未被占用（就我所查）** | 未查到等价的"certificate class = 目标"形式化（检索式：`"uniform computable certificate" gap equivalence class no-go`, `computable functional cannot detect support infimum`）。裁决：**可主张**，但 (a) 需把"不可计算"改成"不可一致可计算"；(b) 需承认它是 `E2-FLIP-NOGO` 四角的重述 |

**总裁决**：本轮的**新意不在 `(NA)`、不在五机制、不在 `1/48`**；唯一可主张的新意是
> **"`F3∧F4=∅` 的正确形态是【族上不可一致可计算】，而非【原子不可计算】；且它与 Tauberian 速率不可计算（APAL 2022）同源。"**

这是一条**收紧**，不是一条进展。

### §3.5 = 第五数据模 / 新坐标归属（对 §5 α/β/γ）
- 本文件**不产出** α/β/γ 候选（任务书只要核查）。
- 但要记一条约束：Q4(b) 的修正意味着 **α（第五数据模）必须同时是"可一致计算的"**，而不是"可计算的"。任何只在固定 cutoff 下可计算的数据模（谱、矩、`γ_ps`）都**自动**落在 `C0`，无论多新。**这条约束会淘汰掉大部分 α 候选。**

### §3.6 = **若文献已解决，指到具体定理**（任务 §6）
| 战役的哪一步 | 文献是否已解决 | 精确指针 |
|---|---|---|
| `(NA)` `Cov ≤ 0` on Haar×Wilson | **否** | 无。最近的是 Bakry–Émery Γ₂（充分、方向相反）：`Hess S ⪰ ρ ⇒ gap ≥ ρ` |
| `Cov(∇S, ∇Φ)` 的符号 | **是（作为恒等式，不是不等式）** | Stein 恒等式 `Cov(Φ,∂_iS)=E[∂_iΦ]`；取 `Φ=∂_jS` 得 `Cov(∂_iS,∂_jS)=E[Hess S]_{ij}` |
| 强耦合质量间隙 | **是（但只在 `β<1/48`）** | SZZ, CMP 400 (2023) 805, arXiv:2204.12737 — 主定理（LSI + Poincaré + 指数簇聚 + 唯一无穷体积极限） |
| 强耦合面积律 | **是（扩张 regime）** | CNS, arXiv:2505.16585；CNS, arXiv:2509.04688；Nissim, arXiv:2510.22788 |
| `U(N)` 的 Bakry–Émery 失效 | **是** | Nissim, arXiv:2510.22788 摘要（Ricci 曲率不一致正） |
| 4D 非微扰 `(L,a)` 一致控号下界 | **否** | 无。Douglas, Nat. Rev. Phys. 8(2) 86–97 (2026) 是综述（**我未读原文**） |
| "谱/原子可计算" | **是，且**对紧算子**可计算** | Brattka–Dillhage, J. UCS 11(12) 1884–1900 (2005)；紧自伴 ⇒ 完全可计算 |
| Tauberian 余项无可计算速率 | **是** | APAL, DOI 10.1016/j.apal.2022.103117；arXiv:1910.07926 |

### §3.7 = 可执行的下一步
1. **清 `1/24`**：按 §2.3 的 8 处清单全局替换为 `β<1/48`，并删"~75×"。**这是本轮最高优先级**——它是**错误**，不是缺证据。
2. **改 `F3∧F4` 口径**：`SPECIES-EXISTENCE.md` §5 与 `C8-SPECIES-NOGO` 中凡以"原子不可计算"为据处，改为"**族上不可一致可计算**"，并引 Brattka–Dillhage 的 (1)(2)(3) 与"紧自伴可计算"作对照——**后者是反例，必须写进去**，否则审稿人一句话打掉。
3. **P3 显式挂 SZZ**：`PROOF-P3-UNCONDITIONAL.md` 开头加一行"阈值 `β<1/48` 与 SZZ (CMP 400, 2023) 一致；本文新意仅在残数 `f₁²` 的一致性下界"。
4. **Douglas 引用改写**：不再写"明说 no clear path"，改为"综述梳理 tensor RG 等路线"（我未读原文，见 §2.3）。
5. **修 Brascamp–Lieb 出处**：`SPECIES-EXISTENCE.md` 参考文献的 `Adv. Math. 22 (1976) 366` → `J. Funct. Anal. 22 (1976) 366–389`。
6. **若要继续找机制**：按 §3.4，**唯一"零占位且零方向冲突"的是机制 1（Harish–Chandra `c`-函数）**，但先要回答"紧群 vs 非紧 `G/K`"的桥问题——建议派一个 200 行以内的探针：在 `SU(2)` 单 plaquette 上算 `c`-函数的**离散类比**（`c(λ)` 的紧群对应物 = Weingarten 系数作为 `λ` 的函数），看它是否**非零且可解析**。若无，机制 1 当场死（成本极低，收益明确）。

---

## §4 附：核实/未核实清单（可审计）

**已核（书目 + 关键论断原文片段）**
SZZ arXiv:2204.12737（摘要全文，含 `1/(16(d−1))`）· CNS arXiv:2505.16585（摘要全文）· CNS arXiv:2509.04688（摘要全文）· Nissim arXiv:2510.22788（摘要全文）· Douglas Crossref 全记录（卷/期/页/日期/DOI/74 refs）· Chatterjee CMP 366 203–268 + DOI · Osterwalder–Seiler Ann. Phys. 110 440–471 + DOI · MRS CMP 155 325–384 + DOI · Balaban CMP 109 249–301 + DOI · Chandra et al. Invent. Math. 237 541–696 · Cao–Chatterjee CMP 405 3 · Pemantle JMP 41 1371–1390 + DOI + arXiv:math/0404095 · BBL JAMS 22 521–567 + arXiv:0707.2340 · Brändén–Huh Ann. Math. 192(3) 821–891 · Collins–Matsumoto–Novak Notices AMS 69(5) 734–745 · Collins–Śniady CMP 264 773–795 · Collins–Nagatsu Algebr. Comb. 9(2) 465–479 + DOI · Gindikin Ann. Fac. Sci. Toulouse 25(2-3) 385–402 + DOI · Brattka–Dillhage J. UCS 11(12) 1884–1900 + Zbl 1120.47063（含四条主结果原文片段）· Paulin arXiv:1212.2015 · Li–Sun arXiv:2408.10169 · Litvinov arXiv:1005.1247 · APAL Tauberian 有穷化 + DOI（含"no computable rate"原文片段）· arXiv:1910.07926 · arXiv:1907.07815（原子⇒可计算片段）· arXiv:2607.22301（stretched LR 仍开）

**未核到（明写检索式，不补）**
- "BDF zoo of negative dependence" 论文 —— 检索式：`Borodin Diaconis Fulman "zoo of negative dependence"`
- Haar×Wilson 的负相依/负相关**定理** —— 检索式：`negative dependence Haar measure compact group lattice gauge theory negative correlation`；`"negative correlation" lattice gauge theory Wilson action plaquette correlation inequality`；`anti-FKG lattice gauge`；`"negative dependence" compact group theorem`
- 把 Weingarten/6j 符号转成反 FKG 结论的定理 —— 检索式：`"negative dependence" fails SU(N) lattice gauge sign problem Weingarten 6j symbol`
- "atom of a measure 的 Weihrauch 度" —— 检索式：`Weihrauch degree "atom of a measure"`；`"atoms of measures are not computable"`
- Harish–Chandra `c`-函数用于格点规范/YM —— 检索式：`Harish-Chandra c-function lattice gauge theory mass gap`；`spherical transform spectral measure atom mass gap`
- 伪谱/Kreiss 用于 YM/Gribov 转移算子 —— 检索式：`pseudospectrum transfer operator spectral gap lower bound`；`Kreiss constant resolvent spectral inclusion gap`
- tropical/max-plus 用于统计力学转移矩阵**严格 gap 界** —— 检索式：`max-plus tropical dequantization transfer matrix statistical mechanics rigorous bound`；`tropical Perron-Frobenius rigorous bound partition function`
- Lorentzian/M-凸 用于格点规范特征标展开 —— 检索式：`Lorentzian polynomials lattice gauge character expansion`；`log-concavity character expansion SU(N) Wilson loop`
- Nazarov–Turán/小球 用于把 `Var` 与 `inf supp` 挂钩 —— 检索式：`small ball probability Nazarov-Turan spectral gap`；`small ball mass gap lower bound`
- SZZ 的正式勘误 —— 检索式：`"Erratum" "A stochastic analysis approach to lattice Yang-Mills at strong coupling"`；zbMATH Zbl 1522.81127 无 erratum 记录。**注意**：Zenodo 灰文献里有"the Shen-Zhu-Zhu erratum"的说法（[zenodo 22099425](https://www.zenodo.org/records/22099425)），**未获正式记录支持，不采信**。
- Douglas 2026 综述**正文** —— `nature.com` 与 PDF 均被 `idp.nature.com` 重定向拦下；只有 Crossref 书目记录。
- CNS 的 `β₀(d)=10^{-10d}/d` 精确式 —— arXiv PDF 不支持文本抓取（`unsupported content type "application/pdf"`）；只有语料 `TRANSPORT-RIGOROUS-NUMERICS.md:89` 的转述，标 **[R-未核原文]**。

**抓错记录（4 条，已标 [X]）**
1. `SPECIES-EXISTENCE.md` 的 "Brascamp–Lieb, *Adv. Math.* 22 (1976) 366" → 应为 *J. Funct. Anal.* **22** (1976) 366–389（`refs.tex:51` 正确）。
2. 任务书/语料的 "Borodin–Diaconis–Fulman negative dependence zoo" → 无此论文；BDF 名下可核为 arXiv:0904.3740。
3. `new-tools-registry.md` F1 的 "Douglas 明说 no clear path known" → 我**未读原文**；二手转述反呈**乐观**语气（Meurice 报告）；引用须改写。
4. `(NA)`/`F3∧F4` 若以"原子不可计算"为据 → **对紧自伴算子为假**；正确形态是"族上不可一致可计算"。

---

## §5 一句话

**`(NA)` 在文献里没有定理，而它的一般近亲（Stein 恒等式 + Bakry–Émery Γ₂）恰好把它翻译成"Hessian 符号"——也就是 H1 的另一个坐标，不是新物种；五机制 0/5 被 gap 语境占用，但三条因原料缺或方向反而不可用；4D 最强半径是 SZZ 的 `β<1/48`（体积一致、非 `(L,a)` 一致），而语料里流通的 `λ_tH<1/24` 在源文中根本不存在，必须全库删除；`F3∧F4=∅` 的可计算性内核已发表（Brattka–Dillhage 的谱度 + APAL 的 Tauberian 无可计算速率），但"原子不可计算"这句话是错的，正确形态是"族上不可一致可计算"——这是本轮唯一可主张的新意，而它是一条收紧，不是进展。**

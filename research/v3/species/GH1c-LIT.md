# GH1c-LIT — 跨年代文献侦察：哪些技术**能生产** H1 要的那种下半连续 / 谱边下界

> **承接**：`SG1`（定理 `SG1-1`：`{G1-协变, F ≤ s₁, 有限矩可算} = ∅`，`SG1.md:62`）· `INEQ-R5-donsker-varadhan.md`（DV 率是 `sup` ⇒ lsc 但类不相交，`:5,7,9,15`）· `GW5-ACID.md:16`（五行杀）· `docs/paper-final/sec01-reduction.tex:220–228`（「七个 uniformity 工具」含 `Γ`-收敛）· `sec01:211–218`（`a=0`/`L=∞` 不是族成员 ⇒ 无紧性）。
> **任务口径**：**不找「谁证过 H1」**，找「哪些技术能生产 lsc / 谱边下界」，并给每条线的**具体断点**。
> **标签**：`[T]` 真定理/恒等式 · `[C]` 条件 · `[NG]` 具名类 no-go · `[M]` 混合 · `[E]` 可复跑数值 · `[J]` 判断 · `[O]` 缺失步/未核实物。
> **阅读等级**：`R2` 读到完整摘要页 · `R1` 读到含关键句的检索片段 · `R3` 仅书目记录（未读摘要）· `R0` = `[O]` 未核。
> **纪律**：AGENTS §1 四要素（§4）· 每个数带 `file:line` · 未证不标 `[T]` · 无 QED · **Clay = 0%** · 不改既有文件、不 commit、只写本文件。
> **Skill**：已扫 skill 目录；`deep-research` 的 lit-review 模式**部分匹配**，但其 13-agent 流水线自带报告骨架与覆盖门，与本文件固定的 §0–§5 骨架、以及「首要失败模式是幻觉 ID」这一风险点不匹配（后者只能靠逐条台账控制，不能靠代理数量控制）⇒ **不加载**，全手写核验。R0 级条目不补、不猜。

---

## §0 判词

### §0.1 侦察的第一结论：**「类外」不等于「lsc 泛函」——两条路都已被夹死**

任务书的直觉是「`SG1-1` 杀的是弱-* 连续泛函 ⇒ 逃到 lsc」。**侦察结果否掉这个读法**，理由是仓内已有一条比 `SG1-1` 更强的刀：

- `INEQ-R5:9` 的 forcing lemma：设 `Ψ` lsc 且 `Ψ(ν) ≤ f₁²(ν)` **逐点**，则把边原子换成 `1/n`-宽 bump 得 `Ψ ≤ 0`。⇒ **任何 lsc 的、逐点从下方压住原子的泛函都恒 `≤ 0`**。`[T]`（该证明只用 portmanteau + lsc 的定义）。
- 所以 lsc **不是**逃逸方向，它是 `SG1-1` 之后**第二个封闭类**。`[T]`

于是真正的问题变成：**有没有一种东西，它不是「测度上的泛函」，或者在薄原子族上根本不被 `f₁²` 压住？** 侦察给出判据如下（本文件新增，`[T]`；**注意它不能写成「不连续」**，见下面的自更正）：

> **C1（类外判据，两条，任一即可）**
> **(a) 不是泛函**：`Q` 的定义域不是 `𝒫(ℝ)`，而需要 `𝒫` 之外的额外结构（算子、场、族参数）⇒ `SG1-1`（「由 `S_W` 直接算出的 `F`」）与 forcing lemma（「`Ψ ≤ f₁²` 逐点」）**两者都不可陈述**。
> **(b) 不是 minorant**：即便 `Q` 是 `𝒫` 上的泛函，只要在薄原子族 `ν_ε` 上 `Q(ν_ε) ↛ 0` 而 `f₁²(ν_ε) = ε → 0`，则 `Q ≰ f₁²` ⇒ 刀砍不到。
> **证明**：(a) 是定义域陈述；(b) 是逐点不等式在一个具体族上的否证。∎
>
> **[X] 自更正（本文件第一稿写错，此处留痕）**：C1 **不能**写成「`Q` 非弱-* 连续 ⇒ 类外」。`lsc` 而非连续的量遍地都是，而 `lsc` 恰是 `INEQ-R5:9` 砍掉的那一类 ⇒ 「不连续」**既不充分**（lsc 也在类内）**也不必要**（(b) 与连续性无关）。任务书「逃到 lsc」的读法在这一行被否掉；能逃的只有 (a) 与 (b)。

**C1 的第二个推论（对 H1 是好消息，也是坏消息）**：仓库自己反复用的**薄原子见证** `ν_ε = ε δ_{m²} + (1−ε) δ_Λ`（`Q6-H1-STATUS.md:208–215` 转述于 `sec01:240–242`）**同时**证明：`f₁² = ε → 0` 而**谱隙 `Δ = Λ − m²` 不变**。⇒ 「隙」不是 `f₁²` 的 minorant，反过来 `f₁²` 也不是「隙」的 minorant。**两者在薄原子族上是正交的两个坐标**。`[T]`（初等；同一族，反方向使用）。
这条同时说明：**能生产「隙」的技术不必能生产「原子」**。本文件因此把两条目标分开记账，不替仓内裁决 `CLAUDE.md §5`（视为同一对象）与 `INEQ-R5:13`（`f₁²` 严格细于 `Δ`）的冲突（见 §5-B2）。

### §0.2 最有希望的两台（各一句断点）

| 名次 | 技术 | 为什么它**在类外**（过 C1） | 它交付什么 | **4D 非阿贝尔中等耦合处的具体断点（一句话）** |
|---|---|---|---|---|
| **①** | **Kato 谱半连续 + Birman–Schwinger 判据**（`[R-09]`,`[R-10]`） | **C1(a)+(b) 双过**：(a) `Q = gap(σ(H))` 要算子结构，不是 `𝒫` 上的泛函；(b) 在薄原子族上 `gap ≡ Λ − m² > 0` 而 `f₁² = ε → 0` ⇒ `gap ≰ f₁²`（见 §2.1(a) 的 `H_ε` 构造） | **「`(c, c+δ)` 内无谱」型下界**（BS：`‖|V|^{1/2}(H₀−λ)^{-1}V^{1/2}‖ < 1 ⇒ λ ∉ σ(H)`）；以及**「隙的搬运」**：norm-resolvent 收敛 ⇒ 有限重孤立本征值连续 | **它要的 norm-resolvent 收敛恰恰是 H1 的另一半**：`a→0` 时粗模算子族活在不同 Hilbert 空间上（`sec01:211–218` 的「族在 `a=0` 无成员」），只可能拿到强-resolvent；而 `min–max` 单独只给 `λ₁` 的 **usc**（`λ₁ = inf_ψ⟨ψ,Aψ⟩` 是连续函数的下确界，`[T]`）⇒ 地板必须**作为输入**进来，循环。 |
| **②** | **红外界（IRB / Gaussian domination，`[R-01]`,`[R-02]`,`[R-37]`,`[R-38]`）** | **C1(a) 过**：它作用在**场**（两点函数）上，是**逐模**不等式 `⟨|σ̂(p)|²⟩ ≤ 1/(2βĴ(p))`，不经过 `ν ↦ ℝ` 的积分泛函。**C1(b) 亦过**（`[J]`）：IRB 是两点函数的上界，与「原子权重」无逐点可比关系 | **方差上界 ⇒ 谱边下界**：susceptibility `χ = Σ_p ⟨|σ̂(p)|²⟩ ≤ C(β) ⇒ m² ≥ 1/C(β)`；这正是 H1「一致方差界」的**同形**陈述（`[T]` 作为蕴含，`[R]` 作为 IRB 本身） | **IRB 的非阿贝尔版不存在，而且 4D 连阿贝尔版都不给地板**：IRB 的推导靠 diamagnetic inequality + 标量/矢量序参量的 RP 比较（`[R-01]` 摘要自述 "strong diamagnetic inequalities; infrared bounds"），非阿贝尔无 diamagnetic 类比；更强的否决是 **Guth 1980（`[R-32]`）**：4D U(1) 有**退禁闭 Coulomb 相** ⇒ 4D 阿贝尔**没有** `β`-一致地板 ⇒ 唯一的「聚合物机给地板」先例（3D U(1)，`[R-03]`）是**维数锁定**，不是群锁定。 |

**一句话合判**：**H1 要的下半连续不在「更聪明的泛函」里，在「算子的谱」和「场的两点函数」这两个不走弱-* 拓扑的层里**；而这两层各自缺的那一样东西，都是同一件——**跨 `a` 的一致性**（Kato 缺 norm-resolvent 一致的算子族；IRB 缺非阿贝尔的 diamagnetic/Gaussian-domination 步）。

### §0.3 其余各条的一句话处置

| 线 | 处置 | 一句 |
|---|---|---|
| 1970s 构造性 QFT（GJS/GRS/Federbush） | **技术源头，非工具** | 高温展开给出的是 `β`-小处的**指数会聚**，不是地板；`[R-05]`,`[R-06]`,`[R-07]` |
| 1980s 聚合物（KP/Seiler） | **只有一台给地板，且维数锁定** | KP 要**绝对可和的活动度**（小参数）；`[R-04]` |
| 1990s–2000s 格点 RG（Balaban/MRS/GK） | **UV 侧** | 全部在**小场/小耦合**；GK 的块自旋是唯一 4D 边缘耦合的严格控制，但控的是**临界（无隙）**理论 `[R-29]` |
| 2000s 熵/曲率（BÉ/HS/Zegarlinski/SZZ） | **方向对、窗口窄** | 曲率判据 `Hess ⪰ ρ > 0` 是充分条件，Savvidy Hessian 不定 ⇒ 按定义失效；SZZ 的窗口 `|β| < 1/48` `[R-22]` |
| 2010s 输运-熵（vRS/Kuwada/Erbar–Maas） | **是等价，不是工具** | `LSI ⟺ T₂` 是**双向**陈述：不产生新信息，只把「隙」翻译成「输运-熵」；见 §1.5 |
| 2020s 概率化格点 YM（Chatterjee/Cao/SZZ 系） | **同窗口** | 全部住在强耦合/'t Hooft；`[R-22]`,`[R-23]`,`[R-24]`,`[R-25]`,`[R-26]` |
| 大偏差下界（FW/DZ） | **另一回事，但同一侧门** | 下界 LDP 是**开集**陈述（`liminf log P(X∈O) ≥ −inf_O I`）；点不是开集 ⇒ 它**证成**而非破除抹平对手（`INEQ-R5:11` 同结论，本文件补上「逻辑形式不同」这一层） |
| `Γ`-收敛（Attouch/Dal Maso） | **论文归类错、判词对** | `Γ-liminf` 由定义即 lsc（§1.3 给一行证明）⇒ 它不是「promotion」，而是**本来就在 lsc 类里** ⇒ 死于 forcing lemma，不是死于「只是搬运」 |
| 自旋玻璃 typical vs uniform（Guerra/Parisi/浓度） | **机制不可搬** | 「典型→一致」靠**高斯浓度** + RSB 插值；格点规范无高斯无序，浓度给的是**聚合量**（自由能）的涨落上界 ⇒ 与 `INEQ-R5` 同一堵墙 |

---

## §1 十条年代线

### §1.1 主表（技术层）

| # | 年代线 | 技术 | 交付什么 | 要什么假设 | 能否跨到中等耦合 `β≈2.5`、4D 非阿贝尔 | 能否给 **lsc** 或**谱边下界** | 文献 | 等级 |
|---|---|---|---|---|---|---|---|---|
| 1 | 1970s 构造性 QFT | Glimm–Jaffe–Spencer 高温展开 | 弱耦合 `P(φ)₂` 的粒子结构、谱条件 | `β`（或耦合）**小**；`d=2` | **否**（小参数展开，`β≈2.5` 处无小参数） | 否（给的是谱条件的**存在**，不是下界常数） | `[R-05]` `[R-06]` | R3 |
| 1′ | 同上 | **红外界 / Gaussian domination**（FSS；BFS 的规范版） | `⟨|σ̂(p)|²⟩ ≤ 1/(2βĴ(p))`，**方差上界** | 反射正性 + 高斯支配 + **diamagnetic**（阿贝尔）+ 标量序参量 | **否**（非阿贝尔无 diamagnetic；4D 阿贝尔被 Guth 否决） | **是**：`χ ≤ C ⇒ m² ≥ 1/C` `[T]`（蕴含） | `[R-02]` `[R-01]` `[R-37]` `[R-38]` | R1 |
| 1″ | 同上 | Federbush 聚合物表述 | 把高温展开写成聚合物气 | 活动度可和 | 否 | 否 | `[R-07]` | R2 |
| 2 | 1980s 团簇/聚合物 | Kotecký–Preiss 判据 | 抽象聚合物模型的会聚准则 | `Σ_{γ∋x}|φ(γ)|e^{a|γ|} ≤ a`（**绝对可和 / 小活动度**） | **否**（`β≈2.5` 时 Wilson 活动度不小；非阿贝尔还有符号问题） | 否（它给会聚，不给地板） | `[R-04]` | R1 |
| 2′ | 同上 | **Göpfert–Mack：3D U(1) 全耦合禁闭** | **`β`-一致的质量隙地板**（唯一已知「聚合物机给地板」先例） | 3D + 对偶到**正** Coulomb 气 + Debye 屏蔽 + Villain 作用 | **否**，且**双重**否：①4D U(1) 本身无一致地板（`[R-32]`）②非阿贝尔对偶非正气体 | **是**（`m ≥ c(β) > 0`, `∀β`） | `[R-03]` `[R-32]` | R1 |
| 2″ | 同上 | Seiler 的规范场构造性综述 | 把上述工具搬到规范场 | 同上 | 否 | 否 | `[R-08]` | R3 |
| 3 | 1990s–2000s 格点 RG | Balaban 多尺度 | UV 稳定性（3D 已证；4D 部分） | 小场/大场分解 + 小耦合 | **否**（只到 UV；`sec01` 的「机器尚不存在」在文献侧成立） | 否（UV 侧无 IR 地板） | `[R-27]` | R1 |
| 3′ | 同上 | Magnen–Rivasseau–Sénéor：YM₄ 固定 IR cutoff | `SU(2)` Schwinger 函数（轴向规范、无 UV cutoff） | 固定 IR cutoff ⇒ **不求隙** | 否 | 否（明确不涉及隙） | `[R-28]` | R1 |
| 3″ | 同上 | Gawedzki–Kupiainen 严格块自旋 | 4D `λ(∇φ)⁴` **临界**理论的构造 | 高斯不动点 + 小 `λ` | 否（控的是**无隙**临界点 ⇒ 与目标反向） | 否（且是「无隙」的证明） | `[R-29]` | R1 |
| 3‴ | 同上 | Wetterich / Polchinski 泛函 RG（严格版） | 精确流方程（恒等式） | 截断 + 小耦合 | **否**（本轮未查到严格版给出 4D YM 一致隙；见 §3-G4 灰文献） | 否（是**方程**，不是不等式） | `[R-30]` `[R-31]` | R1 |
| 4 | 2000s 熵/曲率 | Bakry–Émery `Γ₂` | `Hess S ⪰ ρ > 0 ⇒ gap ≥ ρ` | 曲率**一致正** | **否**：Savvidy `Hess = −D² + 2adF_B` **不定**（仓内实测 7.1% 负）⇒ 假设假 | **是（形式上是谱边下界）**，但前件假 ⇒ 空洞 | `[R-16]` | R1 |
| 4′ | 同上 | Holley–Stroock / Zegarlinski | 格气 LSI 的**逐格**判据 | Dobrushin 唯一性型条件（高温） | 否（高温侧） | 是（隙下界），窗口在高温 | `[R-17]` `[R-18]` | R1 |
| 4″ | 同上 | **SZZ 的 LSI + Poincaré + 指数簇聚** | 无穷体积测度的 LSI/Poincaré，**质量隙严格正** | `|β| < 1/(16(d−1))`（`SU(N)`，与 `N` 无关）`⇒ d=4: β<1/48` | **否**：强耦合；且是**格点单位**（`a→0` 时 `m_phys = m_lat/a` 反而增大） | **是**（文献里最强的「一致隙」）；窗口 `β<1/48` | `[R-22]` | R2（abs 全文已读） |
| 5 | 2010s 输运-熵 | von Renesse–Sturm `T₂ ⟺ LSI` | 曲率 ⇔ 输运-熵 ⇔ LSI 的**等价链** | 等价，非充分条件（循环风险） | 否（等价不产生常数） | **不产生新界**（把「隙」翻译成「输运-熵」，两者同难） | `[R-19]` `[R-20]` `[R-21]` | R1 |
| 6 | 2020s 概率化 YM | Chatterjee；Cao–Chatterjee；Cao–Shen–Zhu；Nissim | 强耦合/'t Hooft regime 的构造、面积律、唯一无穷体积测度 | 小 `β` / 大 `N` / 't Hooft 标度 | 否（同 §1-4″ 窗口） | 是（在窗口内），窗口不动 | `[R-25]` `[R-26]` `[R-24]` `[R-23]` | R2/R3 |
| 7 | 谱理论/扰动论 | **Kato 谱半连续**；min–max；**Birman–Schwinger**；Agmon 隧穿下界 | ①norm-resolvent 收敛下 `σ` 的半连续 + 有限重孤立本征值连续；②`‖K(λ)‖<1 ⇒ λ∉σ(H)`；③本征函数的几何指数下界 | ①norm-resolvent + 固定 Hilbert 空间；②`H=H₀+V` 且 `V` 相对紧、耦合可控；③解析性 + 几何 | **部分**：`a→0` 时族活在不同空间 ⇒ 只有强-resolvent；`β≈2.5` 处无小 `V` | **是（本层唯一原生给出「谱边下界」的）**，但前件见 §2.1 | `[R-09]` `[R-10]` `[R-11]` | R3/R2/R1 |
| 8 | 大偏差下界 | Freidlin–Wentzell / Dembo–Zeitouni 的**下界半** | `liminf log P(X∈O) ≥ −inf_O I`（**开集**） | 开集；速率函数 `I` 的 lsc + 紧水平集 | 否（同 `INEQ-R5` 的结论，只是逻辑形式不同） | **否**：点不是开集 ⇒ 不产生原子质量；**反而证成**抹平对手 | `[R-14]` `[R-15]` | R3/R0 |
| 9 | `Γ`-收敛 | `Γ-liminf` / `Γ-lim`、Dal Maso 基本定理 | `Γ-liminf` **恒 lsc**（一行证明见 §1.3）；等 coercivity 下 `min` 收敛 | coercivity（`Γ`-lim 的 `min` 收敛侧）；本条**不需要** coercivity 就已是 lsc | 否（lsc 类 ⇒ forcing lemma） | **「是 lsc」但方向相反地无用**：它给的是 `inf` 的**下界**，而原子是 usc | `[R-12]` `[R-13]` | R3 |
| 10 | 自旋玻璃/高维统计 | Guerra 插值、Parisi 公式、浓度不等式 | 自由能的**变分上界** + 匹配下界（RSB） | **高斯无序**（插值/浓度都靠 Gaussian） | 否（格点规范无高斯无序；Haar 上无 Gaussian concentration） | 否（浓度是**聚合量**上界，不是原子下界） | `[R-33]` `[R-34]` `[R-35]` `[R-36]` | R1 |

### §1.2 引用明细（`[R-xx]` → ID/DOI + 题名 + 年份 + 阅读等级）

> 规则：**本表只写检索里实际出现的号**。凡 DOI/卷页未在命中中出现，一律不写、改注「未核」。`R2` 读到完整摘要页 · `R1` 读到含关键句的检索片段 · `R3` 仅书目记录 · `R0` = `[O]` 未核。

| Key | 文献 | 年 | 等级 |
|---|---|---|---|
| `[R-01]` | Brydges, Fröhlich, Seiler, **On the construction of quantized gauge fields. I. General results**, *Ann. Physics* **121** 227–284；DOI `10.1016/0003-4916(79)90098-8` | 1979 | R1 |
| `[R-02]` | Fröhlich, Simon, Spencer, **Infrared bounds, phase transitions and continuous symmetry breaking**, *Comm. Math. Phys.* **50** 79；DOI `10.1007/BF01608557`（单源：ADS） | 1976 | R1 |
| `[R-03]` | Göpfert, Mack, **Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant**, *Comm. Math. Phys.* **82** 545–606；Project Euclid `euclid.cmp/1103920658`（**DOI 未核**） | 1982 | R1 |
| `[R-04]` | Kotecký, Preiss, **Cluster expansion for abstract polymer models**, *Comm. Math. Phys.* **103** 491–498；DOI `10.1007/BF01211762` | 1986 | R1 |
| `[R-05]` | Glimm, Jaffe, Spencer, **The particle structure of the weakly coupled P(φ)₂ model and other applications of high temperature expansions, Part I**，in *Constructive Quantum Field Theory*（Velo–Wightman）；CERN record 413610 | 1973 | R3 |
| `[R-06]` | Guerra, Rosen, Simon, **The P(φ)₂ Euclidean quantum field theory as classical statistical mechanics**, *Ann. of Math.*（卷页未核） | 1975 | R3 |
| `[R-07]` | Federbush, **For the Quantum Heisenberg Ferromagnet, a Polymer Expansion and its High T Convergence**, `arXiv:math-ph/0108002` | 2001 | R2 |
| `[R-08]` | Seiler, **Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics**, *Lecture Notes in Physics* **159**, Springer | 1982 | R3 |
| `[R-09]` | Kato, **Perturbation Theory for Linear Operators**, 2nd ed., Grundlehren **132**, Springer；DOI `10.1007/978-3-642-66282-9` | 1995 | R3（**未逐字核**） |
| `[R-09b]` | *Math. Comp.* **36** (1981) 153 的片段：「see, for example, **Kato [6, p. 210]** … lower semicontinuity is true under weaker conditions」 | 1981 | R1（**半连续的二手源**） |
| `[R-10]` | **Operator theoretic methods for the eigenvalue counting function in spectral gaps**, `arXiv:0809.4365`（Birman–Schwinger 型判据的现代处理） | — | R2 |
| `[R-11]` | Agmon, **Lectures on Exponential Decay of Solutions of Second-Order Elliptic Equations…**, Princeton Math. Notes **29** | 1982 | R1 |
| `[R-12]` | Dal Maso, **An Introduction to Γ-Convergence**, Birkhäuser；Zbl `0816.49001` | 1993 | R3 |
| `[R-13]` | Attouch, **Variational Convergence for Functions and Operators**, Pitman | 1984 | R3 |
| `[R-14]` | Dembo, Zeitouni, **Large Deviations Techniques and Applications**, 2nd ed., Springer；DOI `10.1007/978-1-4612-5320-4` | 1998 | R3 |
| `[R-15]` | Freidlin, Wentzell, **Random Perturbations of Dynamical Systems**（版次/DOI **未核**） | — | **R0** |
| `[R-16]` | Bakry, Émery, **Diffusions hypercontractives**, *Sém. Probab. XIX*, LNM **1123**, 177–206 | 1985 | R1 |
| `[R-17]` | Holley, Stroock, **Logarithmic Sobolev inequalities and stochastic Ising models**, *J. Stat. Phys.* **46** 1159–1194 | 1987 | R1 |
| `[R-18]` | Zegarlinski, **Dobrushin uniqueness theorem and logarithmic Sobolev inequalities**, *J. Funct. Anal.*；DOI `10.1016/0022-1236(92)90073-R` | 1992 | R1 |
| `[R-19]` | von Renesse, Sturm, **Transport inequalities, gradient estimates, entropy and Ricci curvature**, *Comm. Pure Appl. Math.*（年/卷页未核） | — | R1 |
| `[R-20]` | Kuwada, **Duality on gradient estimates and Wasserstein controls**（刊/年未核） | — | R1 |
| `[R-21]` | Erbar, Maas, **Ricci curvature of finite Markov chains via convexity of the entropy**（刊/年未核） | — | R1 |
| `[R-22]` | Shen, Zhu, Zhu, **A stochastic analysis approach to lattice Yang–Mills at strong coupling**, `arXiv:2204.12737`；*Comm. Math. Phys.* **400** 805–851；DOI `10.1007/s00220-022-04609-1` | 2022/2023 | **R2（摘要全文已读）** |
| `[R-23]` | Nissim, **U(N) lattice Yang–Mills in the 't Hooft regime**, `arXiv:2510.22788` | 2025 | R2 |
| `[R-24]` | Cao, Shen, Zhu, **Expanded regimes of area law for lattice Yang–Mills theories**, `arXiv:2505.16585` | 2025 | R2 |
| `[R-25]` | Chatterjee, **Yang–Mills for probabilists**, `arXiv:1803.01950` | 2018 | R2 |
| `[R-26]` | Cao, Chatterjee, **A state space for 3D Euclidean Yang–Mills theories**, `arXiv:2111.12813`；*Comm. Math. Phys.* **405** 3（DOI 见 `LIT-CHECK.md:103`） | 2021/2024 | R3 |
| `[R-27]` | Balaban, **Ultraviolet stability of three-dimensional lattice pure gauge field theories**, *Comm. Math. Phys.* **102** 255–275；Euclid `cmp/1104114382`（**与仓内所记 109/1987 不是同一篇**，见 §3-G12） | 1985 | R1 |
| `[R-28]` | Magnen, Rivasseau, Sénéor, **Construction of YM₄ with an infrared cutoff**, *Comm. Math. Phys.* **155** 325–383（NII）/ 325–384（`LIT-CHECK`）；DOI `10.1007/BF02097397` | 1993 | R1 |
| `[R-29]` | Gawedzki, Kupiainen, **A rigorous block spin approach to massless lattice theories**, *Comm. Math. Phys.*（年/卷页未核） | — | R1 |
| `[R-30]` | Wetterich, **Exact evolution equation for the effective potential**, *Phys. Lett. B* **301** 90–94；DOI `10.1016/0370-2693(93)90726-X` | 1993 | R1 |
| `[R-31]` | Polchinski, **Renormalization and effective Lagrangians**, *Nucl. Phys. B* **231** 269–295；DOI `10.1016/0550-3213(84)90287-6` | 1984 | R1 |
| `[R-32]` | Guth, **Existence proof of a nonconfining phase in four-dimensional U(1) lattice gauge theory**, *Phys. Rev. D*（**卷页未核**） | 1980 | R1（两源） |
| `[R-33]` | Guerra, **Broken replica symmetry bounds in the mean field spin glass model**, *Comm. Math. Phys.*（卷页未核） | 2003 | R1 |
| `[R-34]` | Talagrand, **The Parisi formula**, *Ann. of Math.*（卷页未核） | 2006 | R1 |
| `[R-35]` | Aizenman, Sims, Starr, **Extended variational principle for the Sherrington–Kirkpatrick spin-glass model**, *Phys. Rev. B* **68** 214403 | 2003 | R1 |
| `[R-36]` | Panchenko, **The Parisi formula for mixed p-spin models**, `arXiv:1112.4409` | 2011 | R2 |
| `[R-37]` | Bauerschmidt, **Lent 2022 lecture notes**（含 §3.6 "Reflection positivity and proof of IRB"），`cims.nyu.edu/~bauerschmidt/teaching/lm2022/notes.pdf` | 2022 | R1（**正文未逐字核**） |
| `[R-38]` | jstage 期刊片段（`soken` 70(2)）：**Theorem 4.2.4 (Infrared bound or Gaussian domination)** — 「the two-point function is bounded by that of the free or Gaussian field」；**篇名/作者未核** | — | R1 |
| `[R-39]` | Fröhlich, **Massless phases and symmetry restoration in abelian gauge theories and spin systems**, IHES preprint `P/81/40` | 1981 | R1 |

### §1.3 `Γ`-收敛的归类裁决（任务 §9 直问）

**判**：`sec01:220–228` 把 `Γ`-收敛列为「transport 或 promotion（pointwise → uniform, local → global, weak → strong, `L→L+1`）」——**归类不正确**；**但结论（它救不了）正确，理由是另一个**。

**(a) 为什么归类不正确（`[T]`，一行证明）。** 按定义
```
(Γ-liminf_h F_h)(x) := sup_{U ∋ x, U open} liminf_h inf_{y∈U} F_h(y).
```
固定 `U`，函数 `x ↦ liminf_h inf_{y∈U} F_h(y)` 在 `U` 上**是常数**，故在 `U` 上连续（特别地 lsc）；`lsc` 函数族对上确界封闭 ⇒ `Γ-liminf F_h` 是 lsc。**即 `Γ`-收敛的第一半 `Γ-liminf` 不是「把点态升成一致」，它天生就在 lsc 类里**。`Γ`-收敛也不是 promotion，而是**泛函空间上的一种收敛**（Attouch/Dal Maso 的定义对象是「函数列 → 函数」，不是「点态信息 → 一致信息」）。

**(b) 为什么结论仍对（`[T]`）。** 正因为 `Γ-liminf` 是 lsc，它就是 `INEQ-R5:9` 的 forcing lemma 的适用对象：若它还被要求从下方压住原子，则它 `≤ 0`。⇒ `Γ`-收敛属于**已被 `INEQ-R5` 关掉的那个类**，与「搬运」无关。

**(c) 可修的一句。** 把 `sec01:220–228` 的名单改成两组更诚实：**「promotion / transport」六个**（Banach–Steinhaus、Ascoli–Arzelà、Prokhorov、Baire、uniform LDP、`L→L+1` 归纳）＋ **「lsc 型、因而被 `INEQ-R5` 关闭」一个**（`Γ`-收敛）。这不是文字问题：它决定「下一步该往哪找」——不该往**任何**泛函收敛里找。

### §1.4 下界 LDP 与 DV 是不是同一件事（任务 §8 直问）

**判**：**逻辑形式不同，失效机制相同；且下界 LDP 是「证成抹平对手」的那一半。** `[T]`

- DV 率是 `H(Q|μ) = sup_F{E_Q[F] − log E_μ[e^F]}` ⇒ **表现为 lsc 泛函**（`INEQ-R5:5`）。
- LDP 的下界半（Dembo–Zeitouni 型）：`liminf_n (1/n) log P(X_n ∈ O) ≥ −inf_O I`，`O` **开**。它是**对概率的下界**，不是「某泛函可表示为 sup」。
- 关键：**点不是开集**。开集下界对 `O = (m²−δ, m²+δ)` 给的是**邻域质量**的下界；而抹平对手 `ν_n`（边原子变 `1/n`-bump）恰恰使邻域质量**大**、原子质量 `= 0`（`INEQ-R5:11`）。⇒ 下界 LDP 与「原子 `>0`」**相容地共存于抹平族**，因此它**不能**排除 `f₁² = 0`；它甚至**保证**邻域不空。
- 因此答案：**是另一回事**（形式不同：概率界的下半 vs 泛函的 sup 表示），**但不是另一条路**（同一侧门）。

### §1.5 输运-熵是不是「lsc 方向」（任务 §5 直问）

**判**：**都不是**。`[T]`（对「不产生信息」这一条）

`von Renesse–Sturm` 型结果是**等价链**：`Ric ≥ ρ ⟺ T₂(ρ) ⟺ LSI(ρ) ⟺ 熵沿 `W₂` 测地线 `ρ`-凸`。等价意味着：任何一端都**不能**比另一端更容易拿到。仓内的诊断已经写在 `sec01:259–261`：「唯一一条不同类型的下界机制（`vrs` 的输运-熵判据）需要 `β=2.5` 的格点测度不提供的曲率输入」。侦察补充两点：

1. **最优常数是 `inf` 型的比值**，如 `C_LSI = inf_f Ent(f²)/Var(f)`。**下确界**（对连续泛函族）是 **usc**，不是 lsc ⇒ 「隙的下界」不是这类判据的免费产物；这正是 `INEQ-R5` 的对偶面（`[T]`）。
2. `Kuwada` 的对偶（`L^p`-Wasserstein 控制 ⟺ `L^q`-梯度估计）与 `Erbar–Maas` 的「熵凸性 = Ricci 下界」都是**把几何翻译成熵**，不是**造**几何。

⇒ 它们属于「**改写靶**」类，不属于「**生产下界**」类。

---

## §2 深挖前两名

### §2.1 第一名：Kato 谱半连续 + Birman–Schwinger（算子层）

#### (a) 它为什么**在类外**（过 C1，`[T]`）

用 §0.1 的薄原子族。取两点谱测度 `ν_ε = ε δ_{m²} + (1−ε) δ_Λ`，`Λ > m²` 固定，取与之相容的算子

```
H_ε = m² P_ε + Λ (1 − P_ε),      rank P_ε = ε·dim ,
```
于是 `σ(H_ε) = {m², Λ}`、`gap(H_ε) = Λ − m²`（**与 `ε` 无关**），而 `f₁²(ν_ε) = ε → 0`、`ν_ε ⇀ δ_Λ`。

- **C1(a)**：`gap` 要「算子」这个 `𝒫` 外的结构才定义，`F = F(矩)` 的 `SG1-1` 论证对它不可陈述。`[T]`
- **C1(b)**：`gap ≡ Λ − m² > 0` 而 `f₁² → 0` ⇒ `gap ≰ f₁²` ⇒ **forcing lemma 前提不成立**，刀砍不到。`[T]`

这条不是文字游戏：它是仓内「`f₁²` 严格细于 `Δ`」（`INEQ-R5:13`）的**反向使用**——正因为细，粗的那个才躲得开针对细的那个的刀。

#### (b) 它交付什么（三件，按可用性排序）

1. **「无谱」型下界（Birman–Schwinger，`[R-10]`）**：`H = H₀ + V`，若 `‖ |V|^{1/2}(H₀ − λ)^{-1} V^{1/2} ‖ < 1` 则 `λ ∉ σ(H)`。这是**逐点、算子层、非变分**的：不需要 maximizer，不需要凹性，直接排除谱。**这是本文件找到的唯一一台「原生输出谱边下界且不经过测度泛函」的机器。** `[T]`（作为原则）
2. **隙的搬运（Kato 半连续，`[R-09]`）**：norm-resolvent 收敛 + 有限重孤立本征值 ⇒ 本征值连续。这正好是「把 `β₁` 的隙搬到 `β₂`」所需要的**局部化版本**，也是 `sec01:236–267` 所说的「搬过 Savvidy 交叉点」的形状。`[R-09]` p.210 经 `[R-09b]` 转引；**未逐字核原文** ⇒ R3/R1 混。
3. **Agmon 型隧穿下界（`[R-11]`）**：本征函数的几何指数下界 ⇒ 控制「薄原子」的**空间**形状（不是权重）。**用途有限**：它给的是本征函数不消失，不是谱隙。

#### (c) 断点（三条，逐条带数/带出处）

- **断点 1（Hilbert 空间不一致）**：Kato 的半连续是**同一 Hilbert 空间上**的算子列陈述；`a→0` 时 `T_{L,a}` 的载体维数 `~(L/a)^4 → ∞`（`sec01:214`），族在 `a=0` **无成员**（`sec01:212`）。⇒ 只有「强-resolvent 嵌入」可用，而**强**-resolvent 只保住 `inf σ` 的一侧（`[T]`：`inf σ` 在强-resolvent 下**上**半连续），不保地板。
- **断点 2（min–max 的方向）**：`λ₁(A) = inf_{‖ψ‖=1}⟨ψ,Aψ⟩`，是连续（线性）函数族的**下确界** ⇒ `λ₁` 是 `A` 的**凹**泛函 ⇒ **usc**。`[T]`（初等）⇒ **min–max 本身永不可能给出「`λ₁ ≥ c` 一致」**；地板必须以**强制性的二次型下界** `⟨ψ,Hψ⟩ ≥ c‖ψ‖²` 形式**输入**，而那正是 H1。这条把「用 min–max 造地板」的整条支线判死（`[NG]`，不依赖任何模型细节）。
- **断点 3（BS 的 `V` 不够小）**：BS 需要 `H = H₀ + V` 且 `V` 相对 `H₀` 可控。4D 非阿贝尔在 `β≈2.5` 处：`g² = 4/β ≈ 1.6`（`GW5:15` 的约定），`V` 与 `H₀` 同阶 ⇒ `‖K(λ)‖ < 1` 在目标窗口**不可望**；而在 `β<1/48` 处 BS 用得上，但那已是 SZZ 的地盘（`[R-22]`）。

#### (d) 一句话

**Kato/BS 是「把隙搬过去」的正确工具，但它搬不动一个不存在的隙**：它要求的 norm-resolvent 一致性与 `a→0` 的载体变化不相容，而 min–max 只给 usc ⇒ 地板只能作为输入进来。

### §2.2 第二名：红外界 IRB（场层）——为何它「形状对、维数错、群不对」

#### (a) 它交付什么（`[R-02]`,`[R-37]`,`[R-38]`）

标准陈述（`[R-38]` 的 `Theorem 4.2.4` 逐字：「Infrared bound or Gaussian domination … the two-point function is bounded by that of the free or Gaussian field」）：
```
⟨ |σ̂(p)|² ⟩ ≤ 1 / (c · β · Ĵ(p)),   Ĵ(p) = Σ_x J(x)(1 − cos p·x) > 0  (p ≠ 0)
```
**常数 `c` 依归一化而定，本文件不固定它**（`[R-02]`/`[R-38]` 的正文未逐字核 ⇒ 不自造常数）。由 sum rule `χ = Σ_p ⟨|σ̂(p)|²⟩` 得 `χ ≤ C(β)`，再由 `χ ≍ 1/m²` 得

```
m² ≥ 1/C(β)        （谱边下界）      [T]（作为蕴含；IRB 本身是 [R]）
```
**这就是 H1 的同形物**：H1 的「粗糙模算子的一致谱隙」在测度侧读出来正是「一致方差上界」，而 IRB 是文献里**最成体系**的「方差上界 ⇒ 质量下界」机器。它的推导**不经过**任何 `𝒫 → ℝ` 的积分泛函（它是逐 `p` 的不等式），因此**过 C1(a)**。

#### (b) 假设（三条，缺一不可）

1. **反射正性**（RP）：把 `⟨·⟩` 分解到半空间并比较；
2. **高斯支配**：`⟨e^{f}⟩ ≤ e^{(1/2)⟨f,f⟩_J}` 型不等式（`[R-37]` §3.6 的证明路线）；
3. **阿贝尔性的 diamagnetic 步**：`[R-01]` 摘要自述其成果含 "strong **diamagnetic** inequalities; **infrared bounds**" ⇒ BFS 的规范版 IRB 是**挂着 diamagnetic 的**，而 diamagnetic inequality 是 `|e^{iθ}|` 型（阿贝尔）工具。`[J]`（由摘要措辞推断，**未读正文**）

#### (c) 为何「唯一的聚合物地板先例」不可搬（任务 §2 直问，`[J]`+`[R-03]`+`[R-32]`）

Göpfert–Mack 的 3D U(1) 全耦合禁闭是**唯一**已知「聚合物机给地板」的先例。它的三条支柱**全部是 3D 或 U(1) 专属**：

1. **对偶到正 Coulomb 气**：3D U(1) Villain 模型对偶为整数流/单极子气，活动度**正**；KP（`[R-04]`）才可逐字使用。4D 非阿贝尔的对偶是**面/膜**模型，权重**不定号**（Weingarten 符号，`LIT-CHECK.md:53–55` 已核）⇒ 无绝对可和活动度 ⇒ **KP 前提直接失败**（不是「常数不好」，是**前提假**）。`[T]`（作为 KP 前提的否证）
2. **Debye 屏蔽给可算的 `m(β) > 0`，一致于 `β`**：3D 库仑势对数增长 ⇒ 屏蔽长度对一切 `β` 有限。4D 库仑势 `1/r` **不**给这种一致性。
3. **维数锁定（最强的一条）**：4D 阿贝尔 U(1) **本身没有** `β`-一致地板——Guth 1980（`[R-32]`）证明 4D U(1) 存在**退禁闭相**（`[R-39]` 的综述同结论）。⇒ **「4D 的阿贝尔类比都不成立」** ⇒ 把 3D 结果失败归因于「非阿贝尔」是**错误归因**；正确的归因是**维数**。这一条是本轮侦察最有信息量的否定结果，因为它把「聚合物路线在 4D 失效」从「群的非交换性」改写成「`d=4` 的相图」，而后者的难度是已知的（4D U(1) 相变本身是定理）。

#### (d) 4D 非阿贝尔的具体断点（一句 + 一句）

- **一句**：IRB 的「高斯支配」在测度侧等价于「指数函数的二次型在 RP 内积下**正**」；而仓内的墙恰恰是 Savvidy `Hess = −D² + 2adF_B` **不定**（`CLAUDE.md §5`，`W4` 实测 7.1% 负模）⇒ **IRB 与负模在同一方向上对立**。`[J]`（未证蕴含，属假设层对立，不作定理）
- **再一句**：即便 IRB 成立，它控的是**联络**的两点函数（规范依赖）；物理谱隙要的是**规范不变**态上的转移矩阵，两者之间的桥（轴向规范 + 重构）在 `[R-28]` 里被明确**不**用来求隙。

---

## §3 幻觉引用剔除表（本轮实际命中 vs 任务书/仓内转述）

| # | 待查字符串 | 本轮检索结果 | 处置 |
|---|---|---|---|
| **G1** | 「**Brydges–Fröhlich–Seiler–Lehmann**」（任务书 §1 的四作者组合） | 检索式 `"Brydges" "Frohlich" "Seiler" "Lehmann"` → **0 条四作者文献**；可见命中全为文献表里的并列人名。BFS 三作者的 I/II/III 系列**确实存在**（`[R-01]`，Ann. Phys. 121 (1979) 227–284，DOI 已有） | **不写 Lehmann**。★ **本任务首要失败模式命中一次**：该四作者组合无法核实。凡引用一律写 `Brydges–Fröhlich–Seiler`（BFS） |
| **G2** | 「**Alexandrov–Ovchinnikov**」（任务书 §7） | 检索式 `"Alexandrov" "Ovchinnikov" spectral theory paper` → 命中 2 条，均为**气溶胶遥感**（Kassianov/Ovchinnikov/Alexandrov, AMT 2010 等）⇒ **数学文献 0 条** | `[O] 未核`。最近的真实同名线索：`A. Ovchinnikov`（与 Dejak–Lushnikov–Yu–Sigal 合作，非自伴算子谱分析，见 `arXiv:1911.10884` 的引文）与 `A. B. Aleksandrov`（算子函数光滑性）——**两者都不是任务书所指对象，不替代** |
| **G3** | 「**红外界**」 | 检索 `"infrared bound" lattice gauge theory statement Gaussian domination proof` → 命中 IRB 的标准陈述（`[R-38]` Theorem 4.2.4、`[R-37]` §3.6）。**「red external field / 红外界」这个术语在检索里 0 命中** | 判为**术语误写**，按 **infrared bound（红外界）** 处理；不保留「红外界」字样作为术语（只在本行作为剔除记录保留） |
| **G4** | Wetterich/Polchinski「严格版给 4D YM 一致隙」 | 检索 `Wetterich exact renormalization group proof mass gap` → 有效同行评议结果 **0**；命中含 **Zenodo 灰文献**（"Topological Origin of the Mass Gap…", zenodo 17910093） | 方程本身可引（`[R-30]`,`[R-31]`）；「严格版给一致隙」判 **`[O]` 未核**，且灰文献**不计入证据**（沿用 `LIT-CHECK.md:126` 的 Tier-3 口径） |
| **G5** | Göpfert–Mack 的 DOI | 只核到 CMP **82 (1982) 545–606** + Project Euclid `euclid.cmp/1103920658`；**DOI 字符串未在任何命中里出现** | **不写 DOI**，只写卷页 + Euclid id |
| **G6** | Guth 1980 的卷页 | 两次独立命中均只给「Phys. Rev. D, 1980」+ 题名，**未给卷/页** | **不写 `21 (1980) 2291`**（记忆值），只写刊年 |
| **G7** | Guerra 2003 / Talagrand 2006 / vRS 的卷页 | 命中只给刊名+年（Guerra：CMP 2003；Talagrand：Ann. of Math. 2006；vRS：CPAM） | **不写卷页**，只写刊年；标 R1 |
| **G8** | MRS 1993 的页范围 | NII 记 `325–383`，仓内 `LIT-CHECK.md:101` 记 `325–384` | **两记并列，不裁决**（未逐字核） |
| **G9** | Kato 的「lsc 定理号」 | 只核到二手引述：MCom 36 (1981) 153 明写「see, for example, **Kato [6, p. 210]** … lower semicontinuity is true under weaker conditions」 | 写 **p.210 经转引**，**不写定理号**；Kato 原书 `[R-09]` 未逐字核 |
| **G10** | Fröhlich–Simon–Spencer 的 DOI 末两位 | ADS 记录给 `10.1007/BF01608557`（注意：与常见的 `…56` 不同） | 写 **`BF01608557`**（以 ADS 记录为准），并标注「单源」 |
| **G11** | 「Temple/Kato 下界」 | 检索 `Temple inequality lower bound eigenvalue Kato` → 有效书目 **0**（命中为 Zenodo 灰文献） | `[O] 未核`，只在 §2.1(b) 作提示，不进引用表 |
| **G12** | 仓内 `LIT-CHECK.md:100` 的「Balaban, *CMP* **109** 249–301 (1987)」 | 本轮独立命中到的是**另一篇**：Balaban, **Ultraviolet stability of three-dimensional lattice pure gauge field theories**, *CMP* **102** 255–275 (1985)，Euclid `cmp/1104114382`；两篇**都存在**（不同篇） | 两记**并列**，不互相替代。凡引「Balaban 只到 UV」时**不指定**篇号即可；要指定时必须先核 109/1987 的原文（`[O]`） |

---

## §4 检索台账（AGENTS §1 四要素）

**工具限制声明（不是 `TOOL-FAIL`，但必须标）**：`web_search` 返回的是**合并去重后的前 10 条源**（工具自报 "Showing the first 10 sources"），**逐条检索式的独立命中数不可观测**。因此下表「可见命中」一栏报的是**该批次合并后可见源条数**，并标注哪些子式有明确的对应片段。凡「0」均指**该批次内无对应片段**，不等同于「不存在」（见 §5-B1）。

| 批 | 检索式（逐条） | 工具 | 状态 | 可见命中 | 备注 |
|---|---|---|---|---|---|
| S1 | ①`Brydges Fröhlich Seiler infrared bound lattice gauge theory` ②`Göpfert Mack … 1982` ③`Kotecký Preiss cluster expansion abstract polymer models 1986` ④`Fröhlich Simon Spencer infrared bounds … 1976` ⑤`Glimm Jaffe Spencer phase transitions P(phi)2` ⑥`Federbush polymer expansion` ⑦`Guerra Rosen Simon boundary conditions` ⑧`Seiler gauge theories constructive … lecture notes` | `web_search` | **OK** | 10（截断） | ①→BFS I 题名+摘要片段（含 diamagnetic/IRB）②→OSTI+Springer 目录（82, 545–606）③→Scholar（103, 491–498）④→ADS（DOI）⑤→Cambridge 参考文献 ⑥→arXiv abs ⑦→JSTOR 转述 ⑧→NII 书目（LNP 159） |
| S2 | ①`Kato perturbation theory lower semicontinuity of the spectrum norm resolvent convergence` ②`Kato 1958 lower semicontinuity …` ③`Birman Schwinger principle` ④`Alexandrov Ovchinnikov …` ⑤`Agmon Lectures on exponential decay …` ⑥`Courant-Fischer min-max …` | `web_search` | **OK** | 10（截断） | ①/②→MCom 36 (1981) 153 的转引片段（**唯一**可用的 lsc 二手陈述）③→arXiv:0809.4365 **④→0 条数学文献（全是气溶胶）**⑤→HAL 引文（Princeton Math Notes 29, 1982）⑥→CMU 讲义 |
| S3 | ①`arXiv 2204.12737 …` ②`Chatterjee Yang-Mills for probabilists` ③`Cao Chatterjee stochastic quantisation` ④`Shen Zhu Zhu log-Sobolev window 2510` ⑤`Balaban … CMP` ⑥`Gawedzki Kupiainen block spin` ⑦`Wetterich exact RG proof mass gap` ⑧`Polchinski exact RG Fermi surface` | `web_search` | **OK** | 10（截断） | ①→Springer DOI ②→arXiv:1803.01950 ③→引文表给出 **arXiv:2111.12813** ④→（无直接对应；见 S5）⑤→NII（3D UV 稳定性）⑥→NII（标题+摘要首句）⑦→**Zenodo 灰文献 1 条 + 有效 0** ⑧→hep-th/0108018 与多份引文（B231, 269） |
| S4 | ①`Dal Maso Gamma-convergence …` ②`Attouch variational convergence …` ③`Freidlin Wentzell lower bound` ④`Dembo Zeitouni … lower bound open sets` ⑤`von Renesse Sturm transport inequalities` ⑥`Kuwada …` ⑦`Erbar Maas …` ⑧`Guerra interpolation Parisi` | `web_search` | **OK** | 10（截断） | ①→zbMATH Zbl 0816.49001 ②→NII 书目 ③→arXiv:2610.36815 转述 ④→CWI 论文引 Lemma 4.1.6 ⑤→Scholar 摘要片段 ⑥→zbMATH 摘录 ⑦→Erbar 作者页 ⑧→ar5iv 1112.4409 |
| S5 | ①`"Göpfert" … doi` ②`Kotecký Preiss doi` ③`Fröhlich Simon Spencer doi` ④`Kato book doi` ⑤`Dal Maso doi` ⑥`Attouch 1984` ⑦`Dembo Zeitouni doi` ⑧`Alexandrov Ovchinnikov spectral` | `web_search` | **OK** | 10（截断） | ①→Project Euclid 首页（**DOI 未出现**）②→Scopus/Mendeley 双源 `10.1007/BF01211762` ③→ADS `10.1007/BF01608557` ④→arXiv:2404.01422 引文给出 `10.1007/978-3-642-66282-9` ⑤→zbMATH ⑥→NII ⑦→dblp `10.1007/978-1-4612-5320-4` **⑧→0 条数学（气溶胶×2）** |
| S6 | ①`Bakry Émery … LNM 1123` ②`Holley Stroock … 1987` ③`Zegarlinski …` ④`Guerra broken replica … 2003` ⑤`Talagrand Parisi formula 2006` ⑥`Aizenman Sims Starr` ⑦`Cao Chatterjee state space arXiv` ⑧`Gawedzki Kupiainen …` | `web_search` | **OK** | 10（截断） | ①→LNM 1123, 177–206 双源 ②→JSP 46(5/6) 1159–1194 ③→JFA, DOI `10.1016/0022-1236(92)90073-R` ④→Scholar ⑤→Scholar ⑥→PRB 68, 214403 ⑦→引文表 `arXiv:2111.12813` ⑧→NII |
| S7 | ①`"Brydges" "Frohlich" "Seiler" "Lehmann"` ②`Guth 1980 … PRD` ③`Fröhlich Spencer 4D abelian phase transition` ④`Nissim 2510.22788` ⑤`Temple inequality …` ⑥`"infrared bound" lattice gauge statement` | `web_search` | **OK** | 10（截断） | ①→**0 条四作者文献**（文献表并列人名）②→Google Patents+Scholar **两源**（题名+年份）③→IHES P/81/40（4D U(1) 退禁闭）④→arXiv abs 题名 ⑤→**0 条有效** ⑥→**IRB 标准陈述**：jstage Theorem 4.2.4 + Bauerschmidt 讲义 §3.6 |
| S8 | ①`GJS 1973` ②`GRS 1975` ③`Balaban CMP 1987 109` ④`MRS 1993 CMP 155` ⑤`Shen Zhu Zhu master loop` ⑥`Federbush polymer Heisenberg` | `web_search` | **OK** | 10（截断） | ①→CERN 记录 413610（1973）②→Scholar（Ann. of Math. 1975）③→Project Euclid **CMP 102 (1985) 255–275**（与仓内记的 109/1987 **不同篇**，需注意）④→NII（155(2) 325–383）⑤→arXiv:2604.16252 引 `[SSZ25a]` ⑥→arXiv abs |
| F1 | `https://arxiv.org/abs/2204.12737` | `web_fetch` | **OK（HTTP 200）** | 1 | **完整摘要已读**（R2）：含 `|β| < (N−2)/(32(d−1)N)`（SO(N)）与 `|β| < 1/(16(d−1))`（SU(N)）；正式 DOI `10.1007/s00220-022-04609-1` |

**TOOL-FAIL**：**0 条**。（`web_search` 8 批全部有效响应；`web_fetch` 1 次 HTTP 200。）

---

## §5 诚实边界

- **B1（工具分辨率的边界）**：`web_search` 合并截断（前 10 源）⇒ §4 的「0」**只是本批次内无对应片段**，**不是**「文献不存在」的证明。任何 `[O] 未核` 都可能在别处存在；本文件只断言「本轮检索式下未命中」。
- **B2（不裁决 `f₁²` vs `Δ`）**：仓内两种读法并存——`CLAUDE.md §5` 把「隙 / `f₁²>0` / H1 / (CP)」视为同一对象的多个名字，`INEQ-R5:13` 说 `f₁²` **严格细于** `Δ`。本文件 §0.1 的 C1 与 §2.1 的薄原子用法**只**依赖「两者在薄原子族上可分离」这一弱事实，**不裁决**二者是否等价。若等价，则 §2.1 的「类外」论证需改口径（`f₁²` 也会被 C1 归入类外，因为薄原子族同时破坏 `Δ` 与 `f₁²` 的因子化）。`[J]`
- **B3（Kato 的引用层级）**：`[R-09]` 是书，**本轮未逐字核**；「p.210 的半连续陈述」来自 MCom 36 (1981) 153 的转引片段（R1）。任何写进论文的定理号必须另行核原文。
- **B4（IRB 的规范版未核）**：`[R-01]` 只读到摘要片段（"strong diamagnetic inequalities; infrared bounds"）。BFS 的 IRB **具体陈述、适用群、是否覆盖非阿贝尔**——**未核**。§2.2(b)(3) 的「diamagnetic ⇒ 阿贝尔」是由**措辞推断**的 `[J]`，不是定理。
- **B5（Göpfert–Mack 的失效归因）**：§2.2(c) 把「3D 先例不可搬」归因为**三条**（符号/活动度、库仑势、维数）。其中**第 3 条**（4D U(1) 无一致地板，`[R-32]`）是**硬**的；第 1、2 条是标准图像，**未逐字核 GM 原文**，标 `[J]`。
- **B6（未做）**：①没有查任何「谱流 / 指标定理」线（`Fredholm`/`Atiyah–Singer` 方向，仓内 `LOC-DETERMINANT.md` 已另册）；②没有查「非自伴 / 伪谱」线（`LIT-CHECK.md:18,79` 已判部分占用）；③没有核任何正文 PDF（本会话 `web_fetch` 未试 `application/pdf`）；④没有对 `[R-22]`/`[R-23]`/`[R-24]` 的**正文**做逐字核，只读摘要。
- **B7（Clay）**：**0%**。本文件是侦察台账，不产生定理、不改任何等级、不声称 H1 有解。§0 的「两台最有希望」是**阅读判断**（`[J]`），不是可主张的新意；其中唯一 `[T]` 级的新东西是 **C1（类外判据）** 与 §1.3 的 **`Γ-liminf` 一行 lsc 证明**，两者都是**收紧**，不是进展。
- **B8（本文件自身的错，留痕）**：C1 的第一稿写成「非弱-* 连续 ⇒ 类外」，**是错的**（lsc 非连续者遍地，而 lsc 正是被砍的类）。§0.1 已用 `[X]` 自更正为 (a) 定义域 + (b) 非 minorant 两条。凡读者引 C1，**只能引更正版**。
- **B9（`Γ`-收敛的裁决强度）**：§1.3 只裁「`Γ-liminf` 是 lsc」这一条（定义级，`[T]`）与「因而落在被 `INEQ-R5:9` 关闭的类里」（`[T]`）。**没有**裁决 `Γ`-收敛在别的问题里是否好用；也没有裁决 `sec01` 其余六个工具的归类（未核，`[O]`）。

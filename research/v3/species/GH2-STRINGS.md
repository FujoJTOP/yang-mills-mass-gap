# GH2-STRINGS — `σ>0 ⇒ m ≥ c√σ` 的等级、`σ>0` 是否比 `m>0` 更可达、以及这条路是不是换名

> 第 GH2 轮 · 子代理 `gh2-strings` · 工作目录 `D:\Dev\math` · **只写本文件**（未改任何既有文件、未 `git commit`、未跑 MC）
> 标签：`[T]` 真定理 · `[C]` 条件定理 · `[O]` 精确缺失步 · `[NG]` 对具名类的 no-go · `[M]` 缺件/外推 · `[E]` 数值或实测 · `[R]` 标准文献（未开正文）· `[X]` 记账/口径错 · `[J]` 判词
> **Clay = 0%。无 QED。未证之物一律不标 `[T]`。**

---

## §0 判词（先给）

### 0.1 四问四答

| 问 | 答 | 标签 | 依据 |
|---|---|---|---|
| `m ≥ c√σ` 是**定理**还是**启发式** | **都不是**。它是「**窗口内的 `[C]`** + **窗外的 `[O]`**」。仓内枚举 **11 条命题**，**0 条把蕴含本身标 `[T]`**（详见 §1.1 的口径）；外部**无**同行评议定理（只有非评议预印本） | `[C]`/`[O]` | `ROOF-CONFINEMENT-GAP.md:68–74`、`IMP-ARROWS.md:56`、`BLAST-STRING-TENSION.md:24`、`SS-SECTOR.md:53,360` |
| `3.405(21)` 是什么 | **纯测量**（AT2021 Tab.34 的连续极限比值）。**不是**任何定理的常数；窗内可证的常数是 `4√ln(1/u)`，在窗边 `≈9.76(N=2)/10.40(N=3)`，是它的 **2.87–3.06 倍** ⇒ 两者之差**就是**缺的那次搬运 | `[E]`+`[M]` | `LIT-VERDICT.md:75`、`L2-GLUEBALL-LAMBDA.md:36,48`、`ROOF-CONFINEMENT-GAP.md:68`、§1.3 |
| `σ>0` 是否比 `m>0` **更可达** | **否**。两者由**同一台 KP 引擎、同一半径**给出（`uK<1`）；唯一被引作"更宽"的 `β<0.706` 在仓内有**两处相反记录**（说它是非零温 `β_c` 界、不含 `m` 界）⇒ 那句"宽 34×"不能引用 | `[NG]` | `IMP-ARROWS.md:67`、`SS-SECTOR.md:303–307`、`L5-CLASSICAL.md:56,80`、`LIT-VERDICT.md:47`、§2.4 |
| 这条路是不是只是把 `m>0` 换成 `σ>0` | **是：换名 + 一格空洞**。`σ` 是带电(1-形式)数据，对中性(0-形式)质量**无泛函桥**（定理 SS `[T]`）；强耦合唯一那一格与 `σ>0` 的证明**同域锁定**，窗内 `m>0` 已由 OS/SZZ 独立给出 ⇒ **净增益 0** | `[NG]` | `SS-SECTOR.md:353–356`、`IMP-ARROWS.md:197,216–217,230`、§3 |

### 0.2 五条判死线

1. **等级**：仓内**没有任何文件**把 `σ>0 ⇒ m ≥ c√σ` 标 `[T]`。最强的是 `[C]`（`ROOF-CONFINEMENT-GAP.md:3,68–74`、`IMP-ARROWS.md:56`），且它**只在一个半径里**成立：`β*(N) ≤ 8.9×10⁻³`（`ROOF:88`），比目标 `β≈2.47` 低 **278×**。
2. **数值不是定理**：`3.405(21)` 是测量（`LIT-VERDICT.md:75`；`L2-GLUEBALL-LAMBDA.md:36` 逐字「Tab.34 … `M(J^{PC})/√σ` continuum limit」）。仓内已用它反推过 `c₀=1/3.405²=0.086251`（`GD2-WSTAR.md:131`）——**把 3.405 当"定理的 c"是循环**。
3. **可达性**：`σ>0` 与 `m>0` 的已证窗口**同源同半径**（同一个 `uK<1`：`IMP-ARROWS.md:67`；`SS-SECTOR.md:305–307`）。"`σ` 的窗更宽"这个说法在仓内**被两条记录顶掉**（§2.4）⇒ `σ>0` **不**更可达。
4. **换名**：`σ` 是 1-形式带电扇区的率，`m` 是 0-形式中性扇区的谱底；**在装饰封闭类上不存在任何 `F` 使 `F(带电数据) ≤ m`**（定理 SS `[T]`，`SS-SECTOR.md:353`）。所以这条路**不能**是"更弱的充分条件"，只能是"同一个 `Λ` 的两次读数"（`ROOF:74`）。
5. **净增益为 0 的算术**：唯一已证的动力桥的半径 `β* ≤ 8.9e−3` **包含在** `m>0` 已独立证明的区域 `β<1/48` 之内（`IMP-ARROWS.md:56,197,216`；`IMP-ARROWS.md:135` 逐字「那儿的隙早被 OS/SZZ 独立证过 ⇒ **零新地盘**」）⇒ **{桥成立} ⊆ {`m>0` 已知}** ⇒ 任何用这条桥的证明都在已证区域里打转。这就是 `IMP-ARROWS.md:217` 的「有效 **0** 条」。

### 0.3 一句话

**这条路不是新格子：它把一个带电(1-形式)的量当成中性(0-形式)质量的"替代前提"，但 (i) 这两个前提的已证窗口完全重合、(ii) 唯一的动力桥只在隙已被独立证明的地方收敛、(iii) 窗内证出的常数（≥4，窗边 ≈9.8）不是物理的 3.405 且随 `β→0` 发散 —— 所以它既没把 `m>0` 变简单，也没把 `m>0` 换掉，只是给同一个 `u` 加了一次读数。判 `[NG]`（作为绕道），保留 `[C]`（作为强耦合恒等式的读法）。**

---

## §1 `m ≥ c√σ` 的等级：`[C]`（窗内）/ `[O]`（窗外）

### 1.1 逐字核仓内：11 条命题，**0 条把蕴含本身标 `[T]`**

> 口径：`grep c√σ|c′√σ|…|m_\{0⁺⁺\}≥` 原始命中 **57** 行；把"注入/其它对象"的命中剔掉、再逐条去重后 = **11 条命题**（表内 #1–#11）＋ **2 条"形状对但空"**（表后）。#1、#2 的 `[T]` 是**两条腿**（面积律系数、`m₁`）各自的身份，**不是**蕴含本身的身份；#3 才是蕴含本身，它是 `[C]`。

| # | `file:line` | 逐字命题 | 标签 | 证了什么 / 没证什么 |
|---|---|---|---|---|
| 1 | `unconditional/ROOF-CONFINEMENT-GAP.md:61–62` | 「Area law `[T]` …：`σ a² = −ln u + O(u) > 0` for small `β`」 | `[T]` | 面积律**系数**的形式（OS/Wilson 级数） |
| 2 | `unconditional/ROOF-CONFINEMENT-GAP.md:63–64` | 「four-walled tube … `m₁ = 4 ln(1/u)/a` 且 `f₁² ≥ r₀(β)>0`，`L`-uniform」 | `[T]` mod OS/KL | 质量腿（= P3 的 `I1`，不是 `I2`、不是隙） |
| 3 | `unconditional/ROOF-CONFINEMENT-GAP.md:66–74` | 「`m₁/√σ = 4√(ln(1/u)) ≥ 4`」⇒「**`[C]`** `m_{0⁺⁺} ≥ 4√σ·√(ln(1/u)) ≥ 4√σ > 0`」 | **`[C]`** | 窗口内的不等式；窗口由 KP 半径定 |
| 4 | `unconditional/ROOF-CONFINEMENT-GAP.md:88` | 「`β*(N) ≤ 8.9×10⁻³`（优化分级；粗 KP `3.36×10⁻³`；benchmark `1/48 ≈ 2.08×10⁻²`）」 | `[C]` | 桥的半径 |
| 5 | `v3/species/IMP-ARROWS.md:56`（行 4） | 「面积律 **+ `β<β*≤8.9×10⁻³`** ⇒ `m_{0⁺⁺} ≥ 4√σ`」；状态「**`[C]` 已证·同域锁定**」 | **`[C]`** | 与 `σ>0` 消费**同一个** `uK<1` |
| 6 | `unconditional/BLAST-STRING-TENSION.md:24` | 「Bridge `[C]`, β<1/48 … `m₁ ≥ 4√σ`」；`:44` 该界在**无原子**配置上饱和 | **`[C]`** | 只给 `H1`（倒易率），**永远推不到原子 `f₁²`** |
| 7 | `v3/species/SS-SECTOR.md:53,360` | 「不存在 `(L,a)`-一致、`β` 穿过 2.5 的 `m_{0⁺⁺} ≥ c√σ` 定理」（作为**可证性命题**被引用）；「语料行 4 `m_{0⁺⁺}≥4√σ`，`β<8.9e-3`」→ **`[C]`** | `[C]`+`[O]` | 一致版**只能**证"链 I 必须借动力"，**不能**用来补链（`:57–59`） |
| 8 | `unconditional/AREALAW-GAP.md:14,18,27` | 「`m_{0⁺⁺}=c√σ`（`c≈3.5–4`）是**物理恒等式，不是定理**」；缺失步 = 中性道下界 | `[O]` | 面积律在上、`G″` 在下、中间无定理 |
| 9 | `unconditional/TRANSMUTE-TENSION-GAP.md:3,21` | 闭弦 `inf_{L≥L_c}E=0`（`L_c√σ=1.447`）⇒ `inf_L E_closed ≥ c′√σ` **为假** | `[NG]` | 弦图像自己给不出下界 |
| 10 | `unconditional/PROOF-STRING-TENSION-STABILIZER.md:3` | 「the bridging ratio `m_{0⁺⁺}/√σ=const` is a *measured* transmutation regularity, not the area law」 | `[E]`+`[O]` | 比值 = 测量 |
| 11 | `docs/YANG-MILLS-CONDITIONAL-REDUCTION.md:589,597` | 死路台账 #18：`σ > 0 ⇒ m_{0⁺⁺} ≥ c√σ` →「sector … and scale mismatch; **the available inequality points the wrong way**」 | `[NG]` | 已登记 |

补充两处（**形状对但空**）：`FORGE5-AREALAW-GAP-LINK.md:16` 逐字「A bound `m≥c√σ` degenerates correctly (`σ=0` for free field), so the *form* passes the acid test — the route's one genuine attraction. **But no such theorem exists, so it is empty.**」；`LOOP-RUNBOOK.md:40` 逐字「唯一已证的"面积律⟹隙"，`m_{0⁺⁺}≥4√σ` 只在 `β<8.9e−3`，与 `σ>0` **同域锁定**、同一 KP 引擎」。

> **核验判读**：11 条命题里——**把蕴含本身标 `[T]` 的：0 条**；`[C]` **5 条**（#3–#7）；`[O]` **3 条**（#7 与 `[C]` 重合、#8、#10）；`[NG]` **2 条**（#9、#11）；`[E]` **1 条**（#10）；两条腿各自的 `[T]` 只在 #1/#2。仓内**没有任何文件**把蕴含写成定理；`SS-SECTOR.md:53` 那一句是**可证性命题的引用**（"若一致版成立，则链 I 被关闭为不可能"），**不是**存在性断言 —— 引用时必须连这句一起引，否则会把 `[O]` 读成 `[NG]`。

### 1.2 唯一有推导的地方：三行骨架

```
(1)  σ a² = −ln u + O(u)                [T] 窗内（OS'78 / Wilson 强耦合；ROOF-CONFINEMENT-GAP.md:61）
(2)  m₁   = 4 ln(1/u)/a                 [T] mod OS/KL（PROOF-P3-UNCONDITIONAL.md:11–16,22）
(3)  ⇒ m₁/√σ = 4√(ln(1/u)) ≥ 4          [T] 算术（u ≤ e^{−1}）∘ 同一个 u
```
⇒ 结论行本身是 **`[C]`**，因为 `(1)+(2)` 的合并吃掉一条**额外假设**：**`u`-同一性**——"中性道的聚合物类由四壁管主导，且其权重与带电道是同一个 `u`"。这条假设**不是装饰不变的**（`SS-SECTOR.md:305–307` 逐字：「这些箭头的关键假设是"中性道的聚合物类由四壁管主导（且其权重与带电道同一个 `u`）"，而这条假设恰恰不是装饰不变的」），而它只在收敛半径里被验证过 ⇒ 这正是它不能升 `[T]` 的**精确理由**，不是"还没做"。

### 1.3 窗内常数 vs 实测 `3.405(21)`：2.9–3.1 倍，且随 `β→0` 发散

`u = β/(2N²)`（`PROOF-P3-UNCONDITIONAL.md:22` 的 SU(2) 校验式），`m₁/√σ = 4√ln(1/u)`：

| `β` | `u`（N=2） | `m₁/√σ` | `u`（N=3） | `m₁/√σ` | 对 `3.405` 的倍数（N=2 / N=3） |
|---|---|---|---|---|---|
| `1/48 = 2.083e−2` | `2.604e−3` | **9.758** | `1.157e−3` | **10.401** | `2.87` / `3.06` |
| `8.9e−3`（桥边） | `1.113e−3` | 10.432 | `4.944e−4` | 11.036 | `3.06` / `3.24` |
| `1e−3` | `1.250e−4` | 11.992 | `5.556e−5` | 12.521 | `3.52` / `3.68` |

三个后果（全部 `[E]`/`[M]`，不是新定理）：
1. 窗内可证的**常数下界**是 `4`，而**窗边的实际值**是 `9.8–10.4` —— **不是 3.405**。所以"用强耦合证 `m ≥ c√σ` 然后把 `c` 当物理常数"是把两个不同 window 的数混读（AGENTS §2.1）。
2. 它随 `β→0` **发散**（`ln(1/u)` 增长而 `√(…)` 只取半次幂）⇒ 它**不可能**是"同一个物理常数"的低端截断；要把 `4` 连到 `3.405` 必须穿过 `[8.9e−3, 2.47]`（**278×**），这一步就是 `[M]`。
3. 因此 `3.405(21)` 的地位是：**它是"这条路的靶值"，不是"这条路的常数"**。它唯一被仓内用过的合法身份是**归一化锚**（反推 `c₀`，`GD2-WSTAR.md:131`）与**否证阈值**（`DT-DIRICHLET` 的 `m₁ ≥ π√σ = 3.1416` 与 `3.3–3.9` 只差 5–25%，`docs/yang-mills-frontier.md:363`）。

### 1.4 窗外：`[O]`，且缺件已被具名

- `SS-SECTOR.md:53` 把一致版写成**可证性命题**；`SS-SECTOR.md:324` 逐字给出缺件名字：「`e₃` ⟺ 在 `β=2.5` 证明纯 YM 的**带电完备性** ⟺ 证明中性谱无'带电致盲'分量 —— **这就是 Clay 的核心内容**」。
- `AREALAW-GAP.md:27` 的缺件形态：「A **proven neutral-channel bound** `m_{0⁺⁺}≥c√σ` (`c>0`), or an **OS-positive cluster/surface expansion crossing Savvidy (β≈2.47) yielding area law ∧ 0⁺⁺ gap simultaneously**」。
- `TRANSMUTE-TENSION-GAP.md:21–23` 的两面：「(R) scale: `σ_phys>0`（连续统）— 只在 `β<1/48`/`0.706` 证过 `σ_lat>0`；(B) sector: `inf_L E_closed(L) ≥ c′√σ` — **FALSE**」。
⇒ **窗外不是"没证"，是"已知一端为假 + 另一端循环"**（闭弦可缩 ⇒ `inf=0`；要补就把胶球尺寸当输入 = 隙）。

### 1.5 外部检索：无定理

检索式（`web_search`，3 组，见 §5.1）：
- `rigorous theorem glueball mass lower bound string tension lattice Yang-Mills m >= c sqrt sigma` → **0 条同行评议定理**；命中的两条是 **Zenodo 预印本**（其一逐字「We identify `m≥c√σ_str` for an explicit group theoretic constant `c>0`」），**非同行评议、未核正文**，标 `[X-外]`，**不得引用**。
- `"mass gap" "string tension" rigorous inequality lattice gauge theory theorem proven` → 命中的是 (i) 「Toward the Yang–Mills Mass Gap」自述"partial results, failed proof strategies, open problems"；(ii) 2+1 维变分**上界**（`Exact ground state … rigorous upper bounds of the mass gap`）—— 方向与所需相反。
- 结论：**不存在** `m ≥ c√σ` 的定理。仓内对"为什么不存在"的判据是 `BS-BLIND.md:31` 逐字：「A rigorous lower bound `m_{0⁺⁺} ≥ c′√σ` would need a **rigorous flux-tube excitation-gap bound**. The effective-string (Lüscher-term) expansion is asymptotic order-by-order, **not a bound**」。

**⇒ §1 判：`m ≥ c√σ` 的等级 = `[C]`（`β<β*`，`c=4√ln(1/u)`，窗边 ≈9.8）/ `[O]`（窗外）/ `[E]`（`3.405(21)` 是测量）。它不是定理，也不只是启发式：它是一个"在窗口内可证、窗外一端为假、且常数不是物理常数"的条件恒等式。**

---

## §2 `σ>0` 是否比 `m>0` 更可达：**否**

### 2.0 窗表（同一坐标、同一约定 `β=2N/g²`、定 `a`）

| 量 | 已证窗口 | 出处 | 到 `β≈2.47` |
|---|---|---|---|
| `σ_lat>0`（面积律） | `β<1/48 ≈ 2.08e−2`（SZZ，体积一致） | `PROOF-STRING-TENSION-STABILIZER.md:18` | **119×** |
| `m_lat>0`（`0⁺⁺` 隙） | `β<1/48`（SZZ）/ `β<β_OS`（OS'78，常数按群与维数） | 同上 `:6`；`lit-v3/LIT-ERA-04.md:130` | **119×** |
| 桥 `m₁ ≥ 4√σ` | `β ≤ 8.9e−3`（⊂ `1/48`） | `ROOF-CONFINEMENT-GAP.md:88`；`IMP-ARROWS.md:56` | **278×** |
| `σ>0` @ large-`N`/'t Hooft | `λ<λ₀`：`λ₀=1/(8(d−1))=1/24`（Thm 1.6）**或** ~`1e−41`（另一读法，`AREALAW-GAP.md:8` 的 11-24 更正） | `PROOF-CONFINEMENT-RANDOM-SURFACE.md:5,7` | 物理点 `g²N≈3.2` 是阈值的 **~75×**；CNS 自认 open |
| `σ>0` @ `β<12/17≈0.706` | **争议**（§2.4） | 记录 A vs 记录 B | 3.5× |
| 靶 | `β≈2.47`（SU(3) crossover） | `PROOF-STRING-TENSION-STABILIZER.md:6` | 1 |

**读法**：`σ` 与 `m` 的两行**完全同宽**（都 `1/48`），桥更窄（`8.9e−3`）。**唯一**看起来"σ 更宽"的格子是 `0.706` 那一行，它在 §2.4 被顶掉。

### 2.1 (a) 面积律的已知证明：**小 `β`，不是大 `β`；`d≥3` 只买到 `U(1)`；'t Hooft 代数只给二分法**

- **逐字核 `PROOF-P3-UNCONDITIONAL.md`：那里没有 `σ`。** 它的定理是残数 `f₁²_L ≥ r₀(β)>0`（`:11–15`），并在 `:32` 逐字划界：「**Does NOT** give: isolation `I2`, **the strong-coupling gap itself** (that is Osterwalder–Seiler '78 / SZZ '23, not P3), the continuum `a→0`, or the crossover `β→2.5` (the analyticity wall `1/48`, a *method* limit)」。⇒ 把 P3 当"面积律定理"是转引错；`σ>0` 的源是 OS/SZZ，不是 P3。
- **逐字核 `research/mass-spectrum/REPORT.md`：那里没有 `σ` 的正面内容。** `:4` 逐字：「本报告**不重复**禁闭/面积律（`research/confinement/`，那里已说明面积律管的是**带电弦扇区** `E₀(L) ≥ σL`）。本报告处理**中性胶球扇区**」；`:228` 逐字：「`√σ/Λ` certificate，`Q(s)=√σ(s)/s`」对「`m/Λ` certificate，`q(s)=m(s)/s`」——「同一交接判据 `q₀≥q₁`，**不变量不同**」；`:313` 逐字：FM「完全在**带电/弦扇区**，对偶于禁闭与面积律 … 它**不**触及中性胶球谱」。⇒ **中性扇区的权威报告明确把 `σ` 排除在外**，这与 `FORGE5-AREALAW-GAP-LINK.md:3–7` 的电荷分级判据一致。
- **"`β` 大到某个值"方向是反的。** 大 `β`（弱耦合）正是面积律**最难**、且对 abelian **为假**的地方：4D `U(1)` 在大 `β` 是 Coulomb 相（`Guth, PRD 21 (1980) 2291`，`[R]`，经 `IMP-ARROWS.md:63` 与 `PROOF-STRING-TENSION-STABILIZER.md:9` 转引）；CNS 的定理住在 `λ<λ₀` **小耦合**，而 `β=2.5` 对应 `g²N≈3.2`，被 CNS 自己列为 open（`PROOF-CONFINEMENT-RANDOM-SURFACE.md:7` 逐字引 CNS-a Remark 1.3(5)/1.4："it remains a major open problem to show area law for large `λ` and `d=4` for arbitrary `N≥1`"）。
- **`d≥3` 单独不给任何东西。** 唯一的**全耦合**面积律定理是 3D 紧致 `U(1)`（Polyakov 1977；Göpfert–Mack CMP **82** (1982) 545），`IMP-ARROWS.md:63` 标 `[T]` 并注明「**维数+群锁定**：4D 非阿贝尔无对应物」。3D `Z₂` ≡ 对偶 3D Ising 同理（离散群）。
- **'t Hooft 代数 = 二分法，不是选择。** `PROOF-CENTER-VORTEX.md:3` 逐字：「't Hooft's commutation gives only the area/perimeter **dichotomy**, not the selection; no 1-form anomaly forces `σ>0`」；Tomboulis–Yaffe 是**充分条件**且其前提（涡旋凝聚）**就是禁闭重述**（`DEFECT-PEIERLS.md:6,20`）。
- **正权重的缺件**：CPS 的曲面和是 **signed Weingarten**，正性重组是 CNS 自己点名的缺件（`PROOF-CONFINEMENT-RANDOM-SURFACE.md:9` 逐字引 CNS-a Remark 1.4："the positive and negative weights very nearly cancel"）。

### 2.2 (b) 格点单调性：FSS 是**上界、方向错**；SU(N) 无 FKG ⇒ 无 `σ(β)` 单调定理

- **Fröhlich–Simon–Spencer 的高斯支配是上界**：`BASE-INFINITE-VOL.md:41–42` 逐字「**(B) Gaussian domination (Fröhlich–Simon–Spencer) `[T]` — upper bound, wrong direction** … For a *lower* floor on `R_vac` it gives nothing」；`Q1-INFRARED-BOUND.md:401` 把 `{FSS, chessboard, compact-T, heat-kernel, Dobrushin}` 整类判 `[NG]`（到不了 `β=2.5`），且 FSS 在 `k=0` **空转**（`:27–32`，自由无质量参考 `Ĉ(0)=∞`）。
- **没有 FKG ⇒ 没有 `σ(β)` 单调性**：`PROOF-STRING-TENSION-STABILIZER.md:14–15` 逐字「SU(N) has no FKG (PROOF-RANDOM-CLUSTER; PROOF-TRANSPORT-MONOTONE: **no FKG ⇒ no `σ(β)` provably monotone**)」；`BASE-INFINITE-VOL.md:45` 逐字「**no monotonicity of a KL atom weight in `β` is proved** for nonabelian lattice gauge」（那是 `f₁²`，但对 `σ` 同样只有 `[NG]`）。
- **`m` 侧同病**：`B04.md:167` 逐字「`m·a` 的严格单调性未证（`[NG]`）… 文献与本仓都**没有**这个工具」；且即使有，方向**递减** ⇒ 左端点锁的是 **ceiling 不是 floor**（`B04.md:151–156`）。
⇒ 单调性这条路**对两个量都不通**，因此它不构成 `σ` 相对 `m` 的任何优势。

### 2.3 (c) Osterwalder–Seiler / 对偶：`β<1/48` 是 `[T]`，但它是**收敛半径**，不能给整个 `β`

- **是 `[T]`**：OS'78 给「reflection positivity, positive transfer matrix, unique analytic infinite-volume limit, Wilson's area law, exponential clustering with a positive mass gap」（`lit-v3/LIT-ERA-04.md:130`）⇒ 面积律在窗内是定理，**答案：是**。
- **不能给整个 `β`**：`lit-v3/LIT-ERA-04.md:132` 逐字「The **true radius of convergence** of the strong-coupling series is not known … The rigorous `1/48` is therefore a **method radius**, not a phase boundary and not the physical coupling」；`MONOTONE-COUPLING.md:15` 逐字「Analytic `β`-continuation (PROOF-NO-MELT-PHASE: mass is `β`-analytic) carries the **value**, not a **sign-controlled uniform bound**」。⇒ 解析延拓只能搬"数"，不能搬"符号受控的界"，所以把 `σ_lat(β)>0` 从 `1/48` 延到全 `β` 的这一步**在方法上不存在**。
- **对偶只在有正测度的地方存在**：3D 紧致 `U(1)`（单极气体，Polyakov/Göpfert–Mack，全耦合 `[T]`）、3D `Z₂` ≡ Ising（`[T]`）；4D 非阿贝尔**没有**正的对偶面测度（`AREALAW-GAP.md:11` 逐字：「For non-abelian SU(N) the inverse (area-law) series is **Weingarten**, whose local weights are **signed** … There is no FKG …, no order parameter, **no 4D duality**」）。4D `U(1)` 有 Coulomb 相（`Guth`）⇒ 连 abelian 4D 都"不是全 `β`"。⇒ **3D `U(1)` 那条 `[T]` 无法给 4D `SU(N)` 的 `σ>0` 在任意 `β`。**
- 顺带：连续统 `σ_phys>0` 处处 open（`PROOF-STRING-TENSION-STABILIZER.md:20` 逐字「**`σ_phys>0` (continuum):** open everywhere (`a→0` transport)」）——**与 `m` 的连续统是同一堵墙**。

### 2.4 冲突裁决：`β<0.706` 是"面积律"还是"`β_c`"？——**降为 `[X]/[M]`**

| 记录 | 逐字 | 标 |
|---|---|---|
| **A（面积律）** | `PROOF-STRING-TENSION-STABILIZER.md:6,18`：「Widest explicit published range: Borgs–Seiler `β<12/17≈0.706`」；`IMP-ARROWS.md:53`：「`β<12/17≈0.706` Borgs–Seiler」；`CONTINUUM-CP.md:5`：「area law to `β<12/17≈0.71`」；`v3/species/U02-CP.md:103`：「`β≤12/17≈0.71` 面积律 `[T]`」 | `[T]`（转引级） |
| **B（`β_c`、非零温、无 `m` 界）** | `v3/species/L5-CLASSICAL.md:56`：「Borgs–Seiler = 无（是 `β_c` 界）」；`:80`：「Borgs–Seiler 的红外界是**关于 `β_c` 的界**、且在**非零温** … 没有任何 `m` 的界」；`v3/species/LIT-VERDICT.md:47`：「Borgs–Seiler, CMP 91 (1983)｜只给 `β_c` 下界的红外界｜无（且非零温）」 | `[T]`（转引级） |

**外部核**：ADS 记录 CMP **91** (1983) 329–380 的题名是 *“Lattice Yang–Mills theory **at nonzero temperature** and the confinement problem”*，DOI `10.1007/BF01208780`（[ADS](https://ui.adsabs.harvard.edu/abs/1983CMaPh..91..329B/exportcitation)）。**正文在付费墙后，未逐字核**（`[R]`，单源）。
**裁决**：题名与记录 B 一致 ⇒ **0.706 这一行不能当"T=0 面积律半径"引用**；它至多是"非零温的 `β_c` 界"。于是：
- 记录 A 的 4 处引用需加限定（`[X]` 记账）；
- **"σ 的窗比 m 宽 34×"这句话撤回**（`[X]`）；
- **判不变**：不管谁对，2.47 那堵墙没动（0.706 也仍是 3.5× 之外），而 §2.0 的"同宽"结论只依赖 `1/48` 那一行，它两边都无争议。

### 2.5 §2 判

**`σ>0` 不比 `m>0` 更可达。** 三条逐一：
- (a) 面积律只在**小 `β`** 证（`[T]`），`d≥3` 只买到 abelian/离散，'t Hooft 代数只给二分法 ⇒ **没有第二个入口**；
- (b) 没有任何 `β`-单调工具（FSS 上界、无 FKG、无 KL 权重单调）⇒ **不能从窗内推到窗外**；
- (c) `β<1/48` 的 `[T]` 是**方法半径**，解析延拓搬不动符号受控界；唯一的全耦合 `[T]` 住在 3D `U(1)` ⇒ **不能延到整个 `β`**。
**而且它们同窗**：`σ`、`m`、桥三者由**同一个 `uK<1`** 控制（`IMP-ARROWS.md:67`、`SS-SECTOR.md:305–307`）⇒ 可达性**不是** `σ` 的优势，是**同一个限制**。

---

## §3 ★ 这条路是否只是把 `m>0` 换成 `σ>0` 而没有变简单：**是（换名 + 一格空洞）**

### 3.1 五个检验

| 检验 | 结果 | 标签 | 依据 |
|---|---|---|---|
| **① 逻辑类型** | `σ` 是 1-形式带电扇区的率；`m` 是 0-形式中性扇区的谱底。**不同类型** ⇒ 不能互为"更弱前提" | `[T]` | `IMP-ARROWS.md:37–42`（"跨类型时是病态比较"）、`ROOF-CONFINEMENT-GAP.md:18–35` |
| **② 泛函依赖** | **不存在**任何 `F`（连可测性都不要）使 `F(带电数据) ≤ m_n` 在装饰封闭类上一致成立。⇒ **`σ` 的函数（含 `c√σ`）不可能给中性下界** | **`[T]`**（定理 SS） | `SS-SECTOR.md:40,353`（定理 SS）；`:289–296`（推论一：补 `e₃` 必须消费非装饰不变假设，三型穷尽） |
| **③ 来源** | 唯一已证的动力实现是 **`u`-同一性**（同一个聚合物权重在两扇区各读一次率）。这是**恒等式**，不是蕴含 | `[C]` | `ROOF:74`；`IMP-ARROWS.md:118–122`；`T01-IMPLY.md:87` |
| **④ 窗口包含** | `{桥成立}` `=β≤8.9e−3` **⊂** `{m>0 已知}` `=β<1/48` ⇒ 增量**空集** | `[T]`（算术）+`[E]` | `ROOF:88`；`IMP-ARROWS.md:56,197,216`；`IMP-ARROWS.md:135` |
| **⑤ 加宽可组合吗** | `σ>0` 在 large-`N`/'t Hooft 有 `[T]`（CNS），但**那里没有桥**（CNS 的隙独立地骑在 Durhuus–Fröhlich 质量隙条件上）⇒ 前提加宽与箭头**不组合** | `[T]`/`[E]` | `AREALAW-GAP.md:8`（逐字"the gap there rides on the Durhuus–Fröhlich mass-gap condition — again an *independent* input"）；`IMP-ARROWS.md:14`（"CNS 的隙**本来就不来自面积律**"） |

### 3.2 集合论一句（这是"有没有新格子"的全部答案）

```
{β : σ>0 已证}            ⊇ {β≤8.9e−3}          （桥）
{β : m>0 已独立证}         ⊇ β<1/48 ⊇ {β≤8.9e−3}
⇒ {桥成立} ∖ {m>0 已知}    = ∅
```
⇒ **不存在**"`σ>0` 已证、桥已证、而 `m>0` 尚未知"的格子。任何用这条桥的论证都落在"隙已被 OS/SZZ 证明"的区域里 ⇒ **净增益 0**。这正是 `IMP-ARROWS.md:216–217` 的逐字结论：「满足「更弱 + 已证 + 结论 = Clay + 非同域锁定」的：**0 条**」与 `:230` 的四分支之 (iii)。

### 3.3 "同域锁定"是结构性的，不是运气

`SS-SECTOR.md:303–307`（推论三，`[T]`+`[E]`）逐字：「这些箭头的关键假设是"中性道的聚合物类由四壁管主导（且其权重与带电道同一个 `u`）"，而这条假设恰恰不是装饰不变的 … 所以"同域锁定"与"证明必须用动力输入"是**同一件事的两种说法** … **这不是巧合，是定理 SS 的必然**。」
⇒ 强耦合那一格的 `[C]` **不可能**被改造成"更宽"的桥：它之所以成立，靠的正是"两扇区共用一个 `u`"，而这个共用只在收敛半径里被验证 ⇒ 桥的半径 = `σ` 的半径 = `m` 的半径。

### 3.4 唯一可能的新格子（如果一定要找）

按推论一，补 `e₃` 的假设只有三型：**(i)** 中性扇区数据（循环）、**(ii)** 带电×中性的**联合/状态层**约束（动力）、**(iii)** 显式解。`u`-同一性属 (iii) 的强耦合形态。⇒ 真正的"新格子"必须是 **(ii)** 的一种**非 `u`** 形态，具体化即 `SS-SECTOR.md:330` 的 P2：「`H1` 是两点方差（`(d2)` 级），若升级成"带电×中性**联合方差**上界"，它就不再是扇区数据，可能越过装饰障碍」。
**这不是"把 `m>0` 换成 `σ>0`"，而是"造一个联合对象"**；且它一旦成立，`SS-SECTOR.md:324` 说它等价于在 `β=2.5` 证明**带电完备性** = Clay 的核心 ⇒ 那是正面攻坚，不是绕道。

### 3.5 §3 判

**这条路 = 换名 + 一格空洞 ⇒ `[NG]`（作为绕道）。**
- **换名**：把"证 `m>0`"换成"证 `σ>0` 再乘一个常数"；但 `σ` 的函数被定理 SS 禁止给中性下界，乘常数不改变类型。
- **空洞**：唯一已证的乘法（`m₁ ≥ 4√σ`）只在隙已被独立证明的地方收敛；窗外一端为假（闭弦 `inf E=0`）、一端循环（胶球尺寸 = 隙）。
- **它唯一真实的正面内容**（必须记账，不许抹掉）：**强耦合窗内 `σ>0 ⇒ m_{0⁺⁺}≥4√σ` 是可证的真不等式 `[C]`**，且它把"为什么两扇区物理上跟随"解释成"同一个 `u` 的两次读数"。这是**机制**，不是**路线**。
- 与既有台账一致：`docs/YANG-MILLS-CONDITIONAL-REDUCTION.md:589`（#18，`[NG]`）、`LOOP-RUNBOOK.md:40`、`IMP-ARROWS.md:56,197`、`SS-SECTOR.md:53,360`。**本轮没有复活任何一格；本轮新增的是把"能不能更简单"这个问题化成一个集合包含（§3.2）与一条定理（定理 SS 的引用），并撤回一句记账（§2.4 的 34×）。**

---

## §4 酸测（自由场 / 2D / M6 / 单 plaquette）

| 模型 | `σ` | `m` | `σ>0 ⇒ m≥c√σ` 的判定 | 标签 |
|---|---|---|---|---|
| **自由无质量场** | `0`（周长/Coulomb 律） | `0` | `0 ≥ 0`，**空真**：不捏造隙（作为判别器 PASS），但**零信息** | `[T]` |
| **2D 纯 YM** | `(g²/2)C₂(f) > 0` **精确** | `m_n=(2N²/(N²−1))m_ch ≥ 2m_ch > 0` | **成立**（两个都真）——但是**恒等格**，不是证据；且 **`m/√σ ∝ L√σ` ⇒ 无体积一致 `c`** | `[T]`+`[T]*` |
| **M6**（`0+1`：`a` 谐振子 ⊗ `b` 自由环） | **无定义**（无规范场、无 Wilson 环） | `m_L=min(ω,½(2π/L)²)→0` | **上诉词空真 ⇒ 盲**；按 `N02-FOURMODEL.md:21` 的口径"失明 ⇒ 靶不合法" | `[T]` |
| **单 plaquette** | 无定义（无 transfer） | 无定义 | 不可判 | `[T]` |

### 4.1 自由场：过，但过法正是它空的理由

`σ=0`（周长/Coulomb）而 `m=0` ⇒ 蕴含空真。`FORGE5-AREALAW-GAP-LINK.md:16` 已把这一点逐字写成"唯一的真实吸引力"与"空"的同一句话。**注意方向**：`CLAUDE.md §1` 要求新工具"对自由场失效"——这里它确实失效（不产生假隙），但**失效的方式是"没有前件"**，所以它同时也**没有判别力**。

### 4.2 2D：唯一"两个都真"的格，但它把 `c` 变成体积的函数

- **`σ` 是精确的**：2D 面积律精确、`σ_R ∝ C₂(R)`（`IMP-ARROWS.md:57,62`；`SS-SECTOR.md:357`）。**口径修正（必写）**：语料四模型表里的「2D YM `m=3g²/8`」**其实是 `σ` 不是 `m`**（`GN4-REDTEAM.md:242` 逐字引 `D03-ML:33,101–102`：「`3g²/8` 是**弦张力 `σ`（单位长度）**」；外部同向：标准讲义亦给 `σ = 3g²/8 − O(1/g²)` 为强耦合弦张力）。有限环上的隙是 `gap(L)=L·3g²/8`（`ATOM2.md:70`：`3/8, 3/4, 3/2`；`:71` 逐字「`L` 倍是**弦张力×周长**，不是同一读数」）。
- **中性隙也正**：`m_n=(2N²/(N²−1))m_ch ≥ 2m_ch>0`（`SS-SECTOR.md:357`，`[T]`+`[V]`）⇒ **2D 不构成 `σ>0 ∧ m=0` 的反例**。
- **但比值不是常数**：`m_ch = σL`（`ATOM2.md:70–71` 的 `gap(L)=Lσ` 就是带电弦的能量；`:71` 逐字"弦张力×周长"）⇒ `m_n/√σ = (2N²/(N²−1))·L√σ`（SU(2)：`(8/3)L√σ`）。⇒ **`L→0` 时 `m/√σ→0`** ⇒ **不存在体积一致的 `c`**。
- **更狠的一刀**：2D 里靶的**原子**是**空测度**——`ν_{trF²}` 为空（`SS-SECTOR.md:358`，`[T]`+`[V]`；`IMP-ARROWS.md:57` 逐字「2D 精确反例 `σ_fund=g²C₂/2>0` 而 `f₁²(trF²)=0`」）⇒ **`σ>0` 在 2D 甚至推不出"`tr F²` 道非空"**。所以 2D 只支持最弱的那一句（"都有隙"），**不支持**任何通往 Clay 的形式。
- **教训（可搬的部分）**：比值需要一个**长度尺度**（`L_min ≥ κ/√σ`，即胶球尺寸）。这与 `TRANSMUTE-TENSION-GAP.md:21–22` 的缺件 **(B)** 是同一个数 ⇒ 2D 不是反例，而是把缺件**显示出来**。

### 4.3 M6：不是"假"，是"盲"

M6 没有规范场 ⇒ 没有 Wilson 环 ⇒ `σ` **无定义**；而 `m_L = min(ω,½(2π/L)²) → 0`（`U09-M6.md:11,28–29`）。于是命题在 M6 上**既不能真也不能假**（前件不存在）。
按本仓四模型口径，靶必须在 M6 上**为假**（`N02-FOURMODEL.md:21` 逐字：「逐模（`v` = 谐振方向 `a`）`Δ_dyn = 2ω > 0` 而 `m=0` ⇒ **失明，靶不合法**」）；这里对 `σ>0⇒m>0` 发生的是**同一种失明**。⇒ **M6 格不过。**

### 4.4 单 plaquette：无 transfer ⇒ `σ`、`m` 双无定义

（`B04-INDEX.md:62` 同口径；`U09-M6.md` 的 `m` 需 transfer。）

### 4.5 酸测判词

- **自由场 PASS（空真）/ 2D PASS（恒等）/ M6 FAIL（盲）/ 单片 undefined** ⇒ 三个可判格里 **2 过 1 不过**，**但两个"过"都不携带信息**（一个前件不存在、一个是同一个 `u` 的两次读数）。
- **决定性的一条**：**仓内不存在 `σ>0 ∧ m=0` 的模型**——`ROOF-CONFINEMENT-GAP.md:102–117`（§5）穷举后逐字：「we found **no** rigorous model with area law but a gapless `0⁺⁺`; the implication, where testable, holds. But (per §1) it is **not kinematic**, so its 4D truth is dynamical = the wall, not a theorem.」
- ⇒ **这条蕴含在四模型上不可否证**。不可否证的**不是判据**，是"同一个 `Λ` 的恒等式"。**"过酸测"在这里是负面证据，不是正面证据**（与 `FORGE5-AREALAW-GAP-LINK.md:16` 同判）。

---

## §5 台账

### 5.1 工具四要素（AGENTS §1）

| # | 工具 | 状态 | 检索式 / 调用 | 命中 |
|---|---|---|---|---|
| 1 | `pwsh` | **OK** | `Get-ChildItem D:\Dev\math`、`research\{v3,v3\species,unconditional,mass-spectrum}` 枚举 | 有效（`species/` ~700 `*.md`，`unconditional/` ~600 `*.md`，后者输出被截断并落 spill） |
| 2 | `grep` | **OK ×9 / TOOL-FAIL ×1** | ①`3\.405` ②`sqrt\{?σ\|√σ\|string tension\|弦张力\|area law\|面积律` ③`Fröhlich.Simon.Spencer\|FSS` ④`B04\|B04\.md` ⑤`σ\|面积律\|area law\|string tension\|Fröhlich\|monoton` @`mass-spectrum` ⑥`c√σ\|c′√σ\|c\\sqrt\{\\sigma\}\|≥ *c.?√σ\|m_\{0⁺⁺\}≥` ⑦`同域锁定` ⑧`β_OS\|beta_OS\|Borgs–Seiler\|Borgs.Seiler` ⑨`M_\{\-2\}` @`unconditional` | ①22 ②250/1575（截断，全表落 spill）③113 ④170 ⑤11 ⑥57 ⑦17 ⑧44 ⑨11 |
| 2b | `grep` **TOOL-FAIL** | **失败 ⇒ 已换式重试** | `M_\{-2\}\|M_{-2}\|M₋₂` → `rg: regex parse error … repetition quantifier expects a valid decimal`（`{-2}` 被读成区间） | **记 TOOL-FAIL，不记 0 命中**；改写为 `M_\{\-2\}` → 11（见 §5.1 #9） |
| 3 | `read` | **OK**（逐行） | 全文/定段：`ROOF-CONFINEMENT-GAP`(138)、`IMP-ARROWS:1–150`、`SS-SECTOR:40–84,288–367`、`AREALAW-GAP`(31)、`TRANSMUTE-TENSION-GAP`(26)、`PROOF-STRING-TENSION-STABILIZER`(24)、`FORGE5-AREALAW-GAP-LINK`(22)、`DEFECT-PEIERLS`(37)、`PROOF-CONFINEMENT-RANDOM-SURFACE`(15)、`BLAST-STRING-TENSION`(56)、`PROOF-P3-UNCONDITIONAL`(36)、`INEQ-R10`(24)、`BASE-INFINITE-VOL:36–53`、`B04:80–213`、`N02-FOURMODEL:1–70`、`U09-M6:1–60`、`GN4-REDTEAM:225–284`、`TIE-SIGMA`(105)、`mass-spectrum/REPORT:1–30,222–235`、`YANG-MILLS-CONDITIONAL-REDUCTION:575–602`、`GH2-LTART:1–40` | 逐行有效 |
| 4 | `web_search` | **OK**（4+3+3 组） | 见 §5.3 | 10/10/14 条（工具自报截断） |
| 5 | `web_fetch` | 未用 | 摘要级足够；付费墙条目未取正文 | — |
| 6 | `pwsh` + `py -3.12 -I -c` | **OK**，`exit 0` | §1.3 的比值表（纯标量算术，**无 MC**） | 6 行 |

**无其它 TOOL-FAIL。未 `git commit`。本轮只写：`research/v3/species/GH2-STRINGS.md`。**

### 5.2 标签汇总

- **`[T]`**：定理 SS（`SS-SECTOR.md:353`）、推论一/三（`:354,356`）、电荷分级运动学无桥（`ROOF:33`、`FORGE5`#1）、`σ a²=−ln u+O(u)`/`m₁=4ln(1/u)/a` 两条腿（窗内）、2D `σ_R∝C₂(R)` 与 `ν_{trF²}` 空、`{桥}⊆{m>0}` 的算术、四模型四条、`M_{-2}` 给**上界**（`SLOT-BOOTSTRAP.md:95`）。
- **`[C]`**：`m_{0⁺⁺} ≥ 4√σ`（`β≤8.9e−3`）、`m₁ ≥ 4√σ`（`β<1/48`）、`BLAST` 的 `σ⇒H1`（需 `λ≤16σ/M`）。
- **`[O]`**：窗外 `m≥c√σ`；缺件 = `e₃`/带电完备性（`SS-SECTOR:324`）；缺件 = `σ_phys>0` 的 `a→0`（`PROOF-STRING-TENSION-STABILIZER:20`）。
- **`[NG]`**：作为绕道（§3）；闭弦统一下界（`TRANSMUTE:21`）；FSS/RP/单调/FKG 四类（`BASE-INFINITE-VOL:36–50`、`Q1-INFRARED-BOUND:401`）。
- **`[M]`**：把 `c=4` 连到 `3.405` 的那次搬运（278×）；`β*(N)` 的精确值未定（候选 `8.9e−3/3.36e−3/1/48`）。
- **`[E]`**：`3.405(21)`；`9.758/10.401/10.432/11.036/11.992/12.521`；2.47 的 119×/278×/3.5×。
- **`[X]`**：`β<12/17` 的"面积律"读法（§2.4）；四模型表「2D `m=3g²/8`」应为 `σ`（`GN4:242`）。

### 5.3 外部引用（DOI/arXiv + 读到什么程度）

| 文献 | 标识 | 读到 |
|---|---|---|
| Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440–471 | DOI `10.1016/0003-4916(78)90039-8`（Scopus 记录；仓内 `MON2.md:23` 同值） | **摘要/题录级**；正文未开 |
| Shen–Zhu–Zhu, *CMP* **400**(2) (2023) 805–851 | arXiv:2204.12737；DOI `10.1007/s00220-022-04609-1`（Wiley 引文 + arXiv 引文双源） | 摘要 + 仓内逐字核（`LIT-ERA-04:130`） |
| Borgs–Seiler, *CMP* **91** (1983) 329–380 | DOI `10.1007/BF01208780`（[ADS](https://ui.adsabs.harvard.edu/abs/1983CMaPh..91..329B/exportcitation)）；题名含 "**at nonzero temperature**" | **仅题录**；正文付费墙，**未逐字核**（§2.4 的裁决就架在这条限定上） |
| Cao–Nissim–Sheffield | arXiv:2505.16585（Thm 1.2）；arXiv:2509.04688（Thm 1.6） | 摘要 + 仓内逐字核（`PROOF-CONFINEMENT-RANDOM-SURFACE:5–7`） |
| Cao–Park–Sheffield | arXiv:2307.06790 | 摘要 + 仓内 |
| Fröhlich–Simon–Spencer, *CMP* **50** (1976) 79–95 | DOI `10.1007/BF01608557`（仓内 `GH1c-LIT.md:250,269–275` 单源核过） | 标准工具，未开正文 |
| Durhuus–Fröhlich, *CMP* **75** (1980) 103 | 题名见 [Semantic Scholar](https://www.semanticscholar.org/paper/c3bc3e8506da6992ac5a2c193fd9823d0f113cf3) | 题录级 |
| Gaiotto–Kapustin–Seiberg–Willett | arXiv:1412.5148（JHEP 02 (2015) 172） | 题录级 |
| Athenodorou–Teper | arXiv:2106.00364，JHEP **12** (2021) 082（Tab.34 p.82；Tab.15 p.44） | 仓内 `L2-GLUEBALL-LAMBDA.md:36,48` 已逐字核表头 |
| Morningstar–Peardon | *PRD* **60** (1999) 034509 | 题录级（`c≈3.5–4`） |
| Caron-Huot–Pokraka–Zahraee | arXiv:2309.04472，JHEP **01** (2024) 195 | 仓内 `INEQ-R10:12` 逐字引："positivity alone forces the minimal residue to zero" |
| Guth, *PRD* **21** (1980) 2291；Polyakov *NPB* **120** (1977) 429；Göpfert–Mack *CMP* **82** (1982) 545；'t Hooft *NPB* **138** (1978) 1 / **153** (1979) 141；Tomboulis–Yaffe *CMP* **100** (1985) 313；Kotecký–Preiss *CMP* **103** (1986) 491 | 仓内转引（`IMP-ARROWS.md:53–67,369`、`DEFECT-PEIERLS:37`） | `[R]` 未开正文 |

**非同行评议、不得引用（`[X-外]`）**：Zenodo 上的「5.4 Group Theoretic Lower Bound」（自称 `m≥c√σ_str`）、「Toward the Yang–Mills Mass Gap」（自述 partial results / failed strategies）、「A Kinetic-Half-Life Correspondence and the Yang–Mills Mass Gap」（用 `√σ=7/16` GeV 凑 `m=1.75` GeV）。三者**均无同行评议、均未核正文**；列在此处**只为记账"检索确实执行过且命中了什么"**。

### 5.4 未做 / 限定

1. **未**读 Borgs–Seiler 正文（付费墙）⇒ §2.4 的裁决是"题名 + 两条仓内记录"级 `[R]`，不是逐字定理级。
2. **未**核 CNS 的 `λ₀` 真值：仓内两读并存（`1/24` vs `~1e−41`，`PROOF-CONFINEMENT-RANDOM-SURFACE.md:7` vs `AREALAW-GAP.md:8`）。**两者都远在物理点之外**，判不变；但"75×"这个数只在第一种读法下成立 ⇒ 引用时必须带读法。
3. **未**跑任何 MC、**未**改任何既有文件、**未** `git commit`。
4. §1.3 的比值表是**首阶级数**（`u=β/(2N²)`、`σ a²=−ln u+O(u)` 丢 `O(u)`）⇒ 标 `[E]`，不当定理。
5. 2D 的 `gap(L)=Lσ` 采仓内 `[T]*` 口径（`ATOM2.md:70–71`、`GN4-REDTEAM.md:242`）；2D 超可重整、`L` 是物理周长，**其小体积极限不可直接搬到 4D**；能搬的只有"比值需要长度尺度"这一条（§4.2）。

### 5.5 复现（§1.3 表）

```
py -3.12 -I -c "import math;print([(b,N,4*math.sqrt(math.log(2*N*N/b))) for N in (2,3) for b in (1/48,8.9e-3,1e-3)])"
→ (0.0208333,2,9.7576) (0.0208333,3,10.4012) (0.0089,2,10.4316) (0.0089,3,11.0360) (0.001,2,11.9915) (0.001,3,12.5208)
```

### 5.6 Clay

**Clay = 0%。** 本轮没有证任何新定理；新增的是一条集合包含（§3.2）＋一次记账撤回（§2.4 的 34×）＋把"能不能更简单"化成定理 SS 的一次引用。**No QED。**

---

**一句话**：`m ≥ c√σ` 是**窗内的 `[C]`**（`c=4√ln(1/u)`，窗边 ≈9.8–10.4，`β≤8.9e−3`）**加窗外的 `[O]`**，`3.405(21)` 是**测量**；`σ>0` 与 `m>0` **同窗同引擎**（同一 `uK<1`），"σ 更宽 34×"是记账错（0.706 是**非零温 `β_c`**）；唯一已证的桥**只收敛在隙已被独立证明的地方**（`{桥}⊆{m>0}`）⇒ 这条路是**换名 + 一格空洞**，判 `[NG]`；酸测上它"过"的方式（自由场空真、2D 恒等、M6 盲）正是它**没有判别力**的证据。

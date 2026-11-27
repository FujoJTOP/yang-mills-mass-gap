# GY1-NOGAP-CLUSTER — 蕴含链（聚类/多项式衰减/χ-界 ⇒ 无光模）+ 强耦合窗文献映射

> **轮**：GY1 · 子代理 `gy1-nogap-cluster` · 工作目录 `D:\Dev\math`（Windows Git Bash）。
> **写域**：只新建本文件 `research/v3/species/GY1-NOGAP-CLUSTER.md`。**未 `git commit`；未改任何既有文件**（未碰 `docs/position-current.html`、`LOOP-RUNBOOK.md`、`HN1-H1.md`、`Q6-H1-STATUS.md` 等）；零 MC、零 Lean 编译。临时下载件（`~/gy1_szz.html`、`~/gy1_szz.pdf`，在仓库外）已删。
> **靶（逐字，出处 `HN1-H1.md:253–256` §3.3 (c-3)）**：
> ```
> 4D 纯 YM 无光模定理：不存在 L→∞ 上满足 m₁²(L) ≲ L^{−2} 且 f₁² 有正下界（f₁² ≥ r₀ > 0）的族。
> ```
> 载体与量（逐字，`Q6-H1-STATUS.md:76` + `HN1-H1.md:41–47` §1.1 (F1)）：`ν_{L,β}` = 4D SU(N) Wilson 测度上 `O = ¼tr F²` 的 Källén–Lehmann 测度（`C(x)=∫e^{−s|x|}dν(s)` 型分解）；`m₁²(L) = inf supp ν_{L,β}`；`f₁²(L) = ν_{L,β}({m₁²(L)})`；`χ = ∫s⁻¹dν`。族 = `{μ_{L,β}}_{L→∞}`，β 固定（(★) 窗 β∈[1/48, 2.5]，先钉单 β）。
> **酸测 F1（硬要求）**：定理必须在自由无质量场失效（`χ_L ≍ L²`、箱隙 `4π²/L²`、原子 `O(1)` 残数，`Q6-H1-STATUS.md:183` 逐字）；不用 YM 特有输入的"证明"自动作废。
> **标签**：只交 `[T]`/`[C]`/`[O]`/`[NG]`；`[R-corpus]` = 仅仓内转述、本轮未读原文。**No QED；Clay 本轮记 0 项推进。**

---

## §0 一页判词表

| # | 断言 | 判 | 一行理由 / 数 | 出处 |
|---|---|---|---|---|
| 0.1 | **引理 A**：`L`-一致指数聚类（速率 `m_c>0` 任意小）⇒ `m₁²(L) ≥ m_c²` ⇒ 无光模（**不需要 f₁²**） | **`[T]`** | 一行：`C̄(t) ≥ ν̄([s₀,s])e^{−√s t} ⇒ limsup −(1/t)ln C̄ ≤ √s`，而假设给 `liminf ≥ m_c`；令 `s↓s₀` 得 `inf supp ν̄ ≥ m_c²` | 本文件 §1.1；谱表示承重 = `PROOF-P3-UNCONDITIONAL.md:31–33`（离散谱 + `C_L(t)=Σf²e^{−Et}`） |
| 0.2 | **引理 B**：`L`-一致多项式衰减 `\|C(x)\| ≤ A r^{−(2+ε)}`（任 `ε>0`）⇒ `χ_L = o(L²)` ⇒（squeeze + `f₁²≥r₀`）`m₁²(L) ≳ r₀L^{ε−2} ≫ L^{−2}` | **`[T]`** | d=4 壳计数 `(2r+1)⁴−(2r−1)⁴ = 32r³+8r` ⇒ `Σ_x r^{−α} = Θ(L^{4−α})`（α<4）⇒ `o(L²) ⟺ α>2`。**阈值是 α>2，不是 α>3**（α>3 是把"逐点可和需 α>4"混入）；数值核：α=2 时 `Σ/L²→16.13`（饱和=自由场），α=2.5/3 → 0 | 本文件 §1.2；squeeze = `Q6-H1-STATUS.md:90` |
| 0.3 | **引理 C**：`χ_L = o(L²)`（特例 `χ_L ≤ C` 一致）∧ `f₁² ≥ r₀` ⇒ `L²·m₁²(L) ≥ r₀L²/χ_L → ∞` ⇒ 无光模 | **`[T]`** | 一行：`χ_L = ∫s⁻¹dν ≥ f₁²/m₁²(L) ≥ r₀/m₁²(L)`。**f₁² 承重不可去**：thin atom `wδ_ε+(1−w)δ_{s₂}`，`w=ε/2λ`，`χ=1/λ=O(1)=o(L²)` 而 `m₁²=ε→0` | 本文件 §1.3；thin atom = `Q6-H1-STATUS.md:114–116`；Schur 盲 = `HN1-H1.md:266–274` |
| 0.4 | **推论 D（本轮交付的窗口化 T）**：d=4、SU(N)、固定 a、全 L、β_A∈(0,1/48)（= β_W∈(0, N²/48)；SU(2) 1/12、SU(3) 3/16）上，**无光模定理成立**，载子 `O=¼tr F²` | **`[T]`**（D2 逐 L 路模 OS/KL；D1 外部路见 0.5） | **D2（仓内路）**：P3 残数定理 `f₁²_L ≥ r₀(β)>0` 逐 L（`β<β₁=min(1/48,c/N)`）∧ SZZ-路 `χ ≤ C/K_S` 全 L（`Q6-H1-STATUS.md:176`）⇒ 引理 C ⇒ `m₁²(L) ≥ r₀K_S/C > 0` 逐 L | `PROOF-P3-UNCONDITIONAL.md:20–31`（定理逐字）+ `Q6-H1-STATUS.md:176` + §1.3 |
| 0.5 | **D1（独立外部路，同窗）**：SZZ Corollary 4.11（逐字读到）给 `C∞_cyl` 观测量的体积无关指数聚类 ⇒ 引理 A ⇒ `inf supp ν̄ ≥ c_N²` | **`[T]`**（无穷体积测度）/ **`[C]`**（逐有限 L 版：同一 GZ 论证 + 体积一致 Poincaré，但推论原文只对紧极限陈述） | SZZ Cor 4.11 逐字：`Cov(f,g) ≤ c₁d(g)e^{−c_N d(Λf,Λg)}(‖‖f‖‖∞‖‖g‖‖∞+‖f‖₂‖g‖₂)`，`c_N` 依赖 `K_S,N,d`、与 L 无关；`O=¼tr F²` ∈ `C∞_cyl`（12 条 link 的有界光滑函数） | SZZ arXiv:2204.12737v1 §4.3（本轮读 PDF 正文，逐字）；`Q6-H1-STATUS.md:176` |
| 0.6 | **文献映射主题**：(ii) `L`-一致指数聚类在 4D SU(N) 强耦合窗**已证**（OS'78 框架 + SZZ 显式半径 1/48_{β_A}）；(iii) "质量隙"在 **SZZ 的定义（=指数聚类）下已证**、在"0⁺⁺ 原子 + 正残数"意义下**未证**（= atom/H*，正是 P3 的残数定理那一半） | **`[T]`**（映射本身） | SZZ 摘要逐字 "correlations of a large class of observables decay exponentially, namely the infinite volume measure has a strictly positive mass gap"；SZZ §4.3 标题即 "Application II: mass gap"，其隙 = 聚类 | §2 表；SZZ 摘要+§4.3（读到） |
| 0.7 | **(★) 窗覆盖判定**：`(★)=[1/48,2.5]` 的归一化是 P0 混写（`BETA-CONVENTION.md:176` 已判）。β_W 读法：已证子窗 `[1/48, N²/48) ⊂ (★)` 非空（SU(2) `[1/48,1/12)`、SU(3) `[1/48,3/16)`）；β_A 读法：`(★)∩已证窗=∅`（只触 excluded 端点 1/48） | **`[T]`**（两种读法各自算术） | `β_W = N²β_A`（`BETA-CONVENTION.md:64` 逐字核出）；SZZ 窗 `|β_A|<1/48` N 无关（Assumption 1.1：`K_S=(N+2)/2−1−24N|β_A|>0`） | §2.5；`BETA-CONVENTION.md:27,64,115–123,176` |
| 0.8 | **反方向自检**：无任何文献证明 4D SU(N) 全 β 窗（含 β_W=2.5）的指数聚类或质量隙 | **`[T]`**（0 命中） | arXiv "mass gap"+"lattice gauge" 最近 25 篇 sweep：0 篇 4D SU(N) 严格隙；`L4-PREPRINTS.md:177–179` 已杀三条自称（2606.19362 公设 / 2307.00788 标签 / 2506.00284 撤回）；CNS/Nissim 't Hooft 窗仍是强耦合邻域 | §3；检索核验表 R10–R11 |
| 0.9 | **酸测 F1**：三条引理在自由无质量 4D 场**全部正确失效**：(H-a) 失败（无指数聚类）；(H-b) 恰在 α=2 边界失败（`C(x)=1/(4π²x²)`）；(H-c) 失败（`χ_L≍L²`）。在 HN1 否证族 `ν_L=½δ_{1/L²}+½δ₄` 上也全部失效（`χ_L=L²/2+1/8≍L²`，`C̄(t)≍½e^{−t/L}`） | **`[T]`** | 引理不使用任何"对自由场也成立"的输入：聚类/多项式衰减/`χ=o(L²)` 都是自由场**没有**的性质 ⇒ 非空转 | §1.4；`Q6-H1-STATUS.md:183,298–300` |
| 0.10 | 引用更正三条：Göpfert–Mack DOI = `10.1007/bf01961240`（任务猜的 BF01265111 与检索摘要的 BF01212772 均错——后者解析到肝胆外科论文）；Kotecký–Preiss DOI = `10.1007/BF01211762`（任务猜的 BF01240679 解析到 A/D 转换器论文）；"Münster–Weisz《On the convergence of the strong coupling expansion》"**未查到**（INSPIRE 合著仅 1 篇 = 流管厚度 NPB 180 (1981) 1）；Münster 名是 **G.** Münster（仓内 `P3-AUDIT-CORRECTED.md:33` 写 F. 系误） | **`[T]`** | Crossref API 逐条解析核实；INSPIRE 合著检索 total=1 | §2 表 + §4 检索核验表 |
| 0.11 | **Clay**：本轮**不推进**。窗口化 T 在 β_W<N²/48，距靶点 β_W=2.5 约 30×（SU(2)），且即便搬过去也只给"无光模"，不给 0⁺⁺ 原子+间距（I₂） | — | `HN1-H1.md:31`（"H1 即使全证也只到天花板"）同构：无光模是靶的**必要半** | `BETA-CONVENTION.md:122`；`HN1-H1.md:31` |

---

## §1 蕴含链：三条引理（引理级陈述 + 完整一行代数）

### §1.0 场地修补（必须先钉）

靶的量 `m₁²(L)=inf supp ν_{L,β}`、`f₁²=ν({m₁²})` 需要时间轴 KL 分解 `C̄_L(t)=∫e^{−√s|t|}dν̄_L(s)`，而**周期时间 `T_L⁴` 上 `t→∞` 极限不存在**（周期传播子有界）——仓内已钉（`PROOF-P3-UNCONDITIONAL.md:31–32` 逐字 "`T⁴_L` is ill-posed (periodic propagator bounded ⇒ limit does not exist) — hence the cylinder"）。故本文件全部引理在**圆柱 `𝒞_L = T³_L × ℝ`**（空间周期、时间无穷）上陈述；`m₁²(L) = E_{1,L}²`（O-道最低能量的平方）、`ν̄_L = Σ_n f_{n,L}²δ_{E_{n,L}²}`（`f²=|⟨n|O_c|Ω⟩|²≥0`，OS 正性 ⇒ `C̄_L(t) ≥ 0`；`PROOF-P3-UNCONDITIONAL.md:31–33`）。自由场箱隙机制在圆柱上原样存活（`m₁²=(2π/L)²≍L^{−2}`）⇒ 换场地不碰酸测。`χ_L = Σ_{x∈T_L⁴}C_L(x)` 在周期环面上本来就有定义（`Q6-H1-STATUS.md:78–79`）。

### §1.1 引理 A（指数聚类 ⇒ 无光模；不需要 f₁²）`[T]`

> **假设 (H-a)**：存在 `A<∞`、`m_c>0` 与 `L` 无关，使 `0 ≤ C̄_L(t) ≤ A e^{−m_c|t|}` 对一切 `t∈ℤ`、一切 `L`。
> **结论**：`m₁²(L) = inf supp ν̄_L ≥ m_c²` 对一切 `L`。从而 `m₁²(L) ≥ m_c² ≫ L^{−2}`（`L→∞`），无光模定理成立（`f₁²` 不需要）。

**证明（完整代数）**：记 `s₀ := inf supp ν̄_L`。对任意 `s ≥ s₀` 使 `W_s := ν̄_L([s₀, s]) > 0`（这种 `s` 存在且可取 `s↓s₀`，因 `s₀` 在闭支撑里）：
```
C̄_L(t) = ∫ e^{−√u|t|} dν̄_L(u) ≥ ∫_{[s₀,s]} e^{−√u|t|} dν̄_L(u) ≥ W_s · e^{−√s|t|}
  ⇒  limsup_{t→∞} −(1/t)ln C̄_L(t) ≤ √s .                     (A.1)
(H-a)  ⇒  C̄_L(t) ≤ A e^{−m_c t}  ⇒  liminf_{t→∞} −(1/t)ln C̄_L(t) ≥ m_c .   (A.2)
(A.1)∧(A.2) ⇒ m_c ≤ √s 对一切上述 s ⇒ 令 s↓s₀ ⇒ s₀ = inf supp ν̄_L ≥ m_c² . ∎
```
**窗口标注**：(H-a) 的输入与结论都逐 L、固定 β、固定 a；`m_c` 任意小皆可——这正是"比隙弱"（隙定理要 `m_c` 有物理下界，靶只要 `>0`）。**适用场地**：`𝒞_L`（§1.0）。

### §1.2 引理 B（多项式衰减 ⇒ `χ_L=o(L²)` ⇒ 无光模；阈值 α>2 不是 α>3）`[T]`

> **假设 (H-b)**：`|C_L(x)| ≤ A|x|^{−(2+ε)}` 对一切 `x ≠ 0`、一切 `L`（`|x|` = 提升到 `ℤ⁴` 的距离），某 `ε>0`。
> **结论**：`χ_L = o(L²)`，且连同 `f₁² ≥ r₀`：`m₁²(L) ≳ r₀ L^{ε−2} ≫ L^{−2}`。

**Step 1（d=4 求和标度，精确）**：max-范数壳计数 `#{x∈ℤ⁴ : max|x_i|=r} = (2r+1)⁴−(2r−1)⁴ = 32r³+8r`，故（`R := ⌊L/2⌋`）
```
Σ_{x∈T_L⁴\{0}} |x|^{−α} = Σ_{r=1}^{R} (32r³+8r) r^{−α} = 32Σ_{r≤R} r^{1−α} + 8Σ_{r≤R} r^{−1−α}
   =  32·R^{2−α}/(2−α)·(1+o(1)) + O(1)      (α<4)  ⇒ Θ(L^{4−α})
   =  Θ(ln L)                                (α=4)  =  O(1)        (α>4)
⇒ Σ = o(L²)  ⟺  α > 2   （α=2 给 Θ(L²)：恰是自由场饱和）
```
数值核（`py -3.12`，本轮，max-范数精确求和）：`Σ/L²` 在 α=2 时 `21.0→18.3→17.1→16.5→16.3→16.1`（L=4→128，收敛 16.13，**不**是 o(L²)）；α=2.5 时 `13.1→8.4→5.6→3.9→2.7→1.9`（`≍L^{−1/2}→0`）；α=3 时 `8.7→4.2→2.0→1.0→0.50→0.25`（`≍L^{−1}→0`）。**阈值钉死：α>2；α>3 系混淆**（"逐点绝对可和"才要 α>4，那是另一条陈述，这里不需要）。
**Step 2（squeeze 接力）**：
```
χ_L = Σ_x C_L(x) ≤ C_L(0) + Σ_{x≠0}|C_L(x)| ≤ C₀ + A·Θ(L^{2−ε}) =: C₀ + A' L^{2−ε}   (0<ε<2)
squeeze 下 rung（Q6-H1-STATUS.md:90 逐字 f₁²/χ ≤ m₁²，本轮精确复算 5 值全中）：
m₁²(L) ≥ f₁²/χ_L ≥ r₀/(C₀+A'L^{2−ε}) ≍ r₀ L^{ε−2} ≫ L^{−2}   (ε>0)  ∎
```
**窗口标注**：(H-b) 是对**每 fixed β、fixed a、全 L** 的假设；`ε>0` 任意小皆可。`ε=2` 给 `Θ(ln L)`、`ε>2` 给 `O(1)`，二者更过。

### §1.3 引理 C（`χ_L=o(L²)` + `f₁²≥r₀` ⇒ 无光模）`[T]`

> **假设 (H-c)**：`f₁²(L) ≥ r₀ > 0` 且 `χ_L = o(L²)`（特例：`χ_L ≤ χ* < ∞` 与 L 无关——强耦合窗给的就是这个）。
> **结论**：`L² m₁²(L) → ∞`，即 `m₁²(L) ≫ L^{−2}`。

**一行代数**：
```
χ_L = ∫ s⁻¹ dν_L(s) ≥ ν_L({m₁²(L)})·(m₁²(L))⁻¹ = f₁²/m₁²(L) ≥ r₀/m₁²(L)
  ⇒  m₁²(L) ≥ r₀/χ_L   ⇒   L²·m₁²(L) ≥ r₀·L²/χ_L → ∞   (χ_L=o(L²))  ∎
```
**f₁² 承重不可去（精确反例边界）**：thin atom `ν = wδ_ε + (1−w)δ_{s₂}`，`w = ε/(2λ)`：`χ = 1/λ = O(1) = o(L²)` 而 `m₁² = ε → 0`（`Q6-H1-STATUS.md:114–116`；Lean `H1_does_not_imply_atom`）。纯 `χ` 信息（矩型泛函）按 Schur-凹性看不见支撑底（`HN1-H1.md:266–274`）——**这就是为什么 (H-c) 必须"连同 `f₁²≥r₀` 一起用"**，也正是本路线必须接 P3 残数定理的原因。

### §1.4 推论 D（窗口化 T——本轮交付物）与酸测

```
D2（仓内路，逐 L）：  β < β₁ = min(1/48, c/N)（β_A 口径）：
   f₁²_L ≥ r₀(β) > 0  逐 L    [T mod OS/KL]   PROOF-P3-UNCONDITIONAL.md:20–31（定理逐字）
   χ_L ≤ C/K_S        全 L    [T]             Q6-H1-STATUS.md:176（SZZ-BE 路，体积一致）
   ⇒ 引理 C ⇒  m₁²(L) ≥ r₀(β)·K_S/C > 0  逐 L  ⇒ 无光模定理（β_A∈(0,1/48), 固定 a, 全 L）。

D1（外部路，独立）：  β_A < 1/48（Assumption 1.1, SU(N), d=4）：
   SZZ Corollary 4.11（逐字）⇒ (H-a) 以 m_c = c_N(β,N,d) > 0（体积无关常数）
   ⇒ 引理 A ⇒ inf supp ν̄ ≥ c_N²   [T]（无穷体积测度）；逐有限 L 版 = [C]（同一 GZ 论证 +
   SZZ §4.1 体积一致 Poincaré；推论原文只对 tight limit 陈述——缺的不是机制是逐字行）。
```
**酸测 F1（三条引理全部正确失效）**`[T]`：
- 自由无质量 4D 场：`C(x)=1/(4π²x²)` 恰 α=2 ⇒ (H-b) 对一切 ε>0 失败；无指数聚类 ⇒ (H-a) 失败；`χ_L≍L²` ⇒ (H-c) 失败（`Q6-H1-STATUS.md:183`）。
- HN1 否证族 `ν_L = ½δ_{1/L²}+½δ₄`（`HN1-H1.md:229–231`）：`C̄_L(t) = ½e^{−t/L}+½e^{−2t}` 在 `t≤L` 为 O(1) ⇒ (H-a)(H-b) 失败；`χ_L = L²/2+1/8 ≍ L²` ⇒ (H-c) 失败。
⇒ 三条引理的假设都**不是**自由场具有的性质，没有空转。

---

## §2 文献映射（主任务）：4D SU(N) Wilson 的 L-一致聚类/隙，谁在哪个 β 窗证明了什么

### §2.1 三格判据表（(i) 面积律 / (ii) 指数聚类 / (iii) 质量隙——三者差几格）

| 文献 | (i) 面积律 | (ii) 指数聚类 | (iii) 质量隙/转移矩阵隙 | β 窗（含归一化） | 读到 |
|---|---|---|---|---|---|
| **Osterwalder–Seiler '78**, Ann. Phys. **110** 440–471, DOI `10.1016/0003-4916(78)90039-8`（`LIT-ERA-04.md:197` VERIFIED） | ✓（"verify Wilson's confinement bound"，摘要逐字） | ✓（SZZ 正文转述逐字："as explained in [OS78], the existence of a mass gap (exponential clustering) follows from convergence of the cluster expansion"；观测量类 = "bounded functions of finitely many edge variables"，SZZ §1 转述） | **摘要不写隙**；转移矩阵存在性在摘要（"existence of a positive self-adjoint transfer matrix"）。"uniform strong-coupling gap"的归属（`LIT-CHECK.md:99`、`P3-CYLINDER-RESIDUE.md:186`）靠 SZZ 的已发表转述 + 原文正文（付费墙未读） | 强耦合某邻域 `β < β₀^{OS}`，**显式数值 [O]**（`M3B-OSWIN.md:76` 独立复核：语料+web 六组检索 0 命中显式数）；OS 自己不带 't Hooft 缩放，"βN small"（SZZ §1 逐字） | 摘要（INSPIRE 逐字）；正文未读 |
| **Seiler, LNP 159 (Springer 1982)**, DOI `10.1007/3-540-11559-5`（Crossref 核，本轮） | 转述级 | 转述级（"the book [Sei82] for these expansion techniques and results"，SZZ §1 逐字） | 同上 | 同 OS'78 | 只读到著录；`Z02-SEILER.md:27` 记 Springer 重定向 TOOL-FAIL |
| **Shen–Zhu–Zhu '23**, CMP **400** 805, arXiv:[2204.12737](https://arxiv.org/abs/2204.12737), DOI `10.1007/s00220-022-04609-1` | —（不做面积律） | **✓ 显式定理**：Corollary 4.11 逐字 `Cov(f,g) ≤ c₁d(g)e^{−c_N d(Λf,Λg)}(‖‖f‖‖∞‖‖g‖‖∞+‖f‖₂‖g‖₂)`，`f,g∈C∞_cyl`，常数体积无关；`c_N` 依赖 `K_S,N,d` | **✓ 在"隙=聚类"义**（摘要逐字 "strictly positive mass gap"；§4.3 标题 "Application II: mass gap"）；**✗ 不给 0⁺⁺ 原子/残数/间距** | **`|β_A| < 1/48`**（SU(N), d=4, N 无关；Assumption 1.1 ⟺ (1.3) 逐字）；β_W 下 = `N²/48`（SU(2) 1/12、SU(3) 3/16） | **摘要 + 正文 §1.1–1.2、§4.1、§4.3（本轮读 arXiv v1 PDF，43 页，逐字）** |
| **Kotecký–Preiss '86**, CMP **103** 491–498, DOI `10.1007/BF01211762`（Crossref 核） | 工具，不是定理 | 工具（KP 判据供聚合展开收敛） | 工具 | 通用抽象聚合模型 | 著录（Crossref；`LIT-ERA-04.md:259` VERIFIED） |
| **Göpfert–Mack '82**, CMP **82** 545–606, DOI `10.1007/bf01961240`（Crossref 核） | ✓（3D U(1)） | — | "absence of massless photons"（摘要逐字）——**3D U(1) + 3D CP¹，非 4D SU(N)** | `β < β₁ = 1.668`（3D U(1) 自有口径） | 摘要 |
| **Borgs–Seiler '83**, CMP **91** 329–380, DOI `10.1007/BF01208780`（INSPIRE 作者核：Borgs, Christian + Seiler, Erhard） | `β < 12/17 ≈ 0.706` = **`[R-corpus]`**（`cpb-attack/09-area-law/REPORT.md:8,35`、`CONTINUUM-CP.md:5`、`DEFECT-PEIERLS.md:6` 三处一致转述；摘要无此数；**β 归一化未核**——原文正文未读） | — | — | 12/17 的口径未核 ⇒ 不进判词 | 摘要（INSPIRE："string of inequalities between the corresponding string tensions"） |
| **Münster（G.）**，"The glueball mass in strong coupling lattice gauge theories" PLB **95** (1980) 59；"Physical strong coupling expansion parameters and glueball mass ratios" PLB **121** (1983) 53（INSPIRE 核）；NPB **190** (1981) 439（`LIT-ERA-04.md:250` 标 UNVERIFIED） | — | — | **级数** `m = −4ln u + Σ_{k≥1}m_k u^k`（`P3-AUDIT-CORRECTED.md:33`）——`[C]/[V]`，非定理 | 强耦合级数，收敛性未严格化 | 著录（INSPIRE）；正文未读 |
| **"Münster–Weisz《On the convergence of the strong coupling expansion》"**（`REPORT.md:289`） | **未查到** | — | — | — | INSPIRE 合著检索 total=1（仅 "How Thick Are Chromoelectric Flux Tubes?" NPB **180** (1981) 1，DOI `10.1016/0550-3213(81)90151-6`）；`ti "strong coupling expansion"` 全库 137 条无此题。`REPORT.md:289` 自标"未逐条核实" ⇒ **该引用名不可用** |
| **Faria da Veiga–O'Carroll**, arXiv:[2509.03513](https://arxiv.org/abs/2509.03513)（**最新 v2，2026-01-24**，arXiv API 核） | — | — | **否**：最新版摘要逐字 "in the gluon field expansion of the action, we show there are no local terms, such as a mass term. However, the action associated with the exponent of the exponentiated Haar measure density has a positive local mass term ∝ β=g⁻² and the SU(N) quadratic Casimir" + "at least for small β, a multiplicity-two one-particle glu[on state is expected]"——**仍是 expected，无质量下界** ⇒ 与 `B05.md:128`、`SURVEY-RECENT-ARXIV.md:37` 的仓内判定一致，无需改账 | "small β"（无显式半径） | 摘要（arXiv API 逐字） |
| **Cao–Nissim–Sheffield**（arXiv:2509.04688 / 2505.16585）、**Nissim**（arXiv:2510.22788 "U(N) … 't Hooft regime"：摘要逐字 "We establish a mass gap…in the 't Hooft regime"） | ✓（动态法验证 DF80 前提） | ✓ | ✓（'t Hooft 窗） | 仍是强耦合邻域（`|β| < 1/(16(d−1))` 型）；不越 (★) 上端（`L4-PREPRINTS.md:104–116`） | 摘要（仓内已核）+ 本轮未重复 |
| **Balaban** 系列（CMP **109** (1987) 249–301，DOI `10.1007/BF01215223`；`RG2-MAP.md:86–90` 逐条核） | — | — | **UV 稳定性，非 IR 隙**（"由此完成 4D 纯规范场紫外稳定性的证明"——UV 端） | UV | 著录+摘要（仓内核） |

### §2.2 核心检索题的答案（(iii) 4D SU(N) 强耦合质量隙到底有没有人证过？）

**答案分两半（这是本路最重要的分辨）**：
1. **"隙 = 指数聚类"的意义：已证。** OS'78 建立强耦合框架（转移矩阵 + 全观测量聚合展开体积一致收敛，观测量 = "bounded functions of finitely many edge variables"，SZZ §1 逐字转述），并按 SZZ 的已发表转述给出 mass gap = exponential clustering；**SZZ'23 给出显式窗** `|β_A|<1/48`（SU(N)、d=4、N 无关）+ 体积无关常数 + LSI/Poincaré + 无穷体积唯一性。CNS/Nissim 把同型结果扩到 U(N)/'t Hooft regime，窗仍是强耦合邻域。
2. **"隙 = 0⁺⁺ 原子 + 正残数（+I₂ 间距）"的意义：未证。** SZZ 的聚类对 `C∞_cyl` 观测量是**协方差衰减**，不给出 `f₁² = ν({m₁²}) ≥ r₀` 的原子——那是仓内 P3 残数定理（`PROOF-P3-UNCONDITIONAL.md:20–31`，`[T] mod OS/KL`，逐 L、圆柱、`β<β₁`）补的另一半。**"转移矩阵谱隙 + 谱原子"的合取式，文献里没有任何 4D SU(N) 版本。**

### §2.3 (ii) 聚类与 (iii) 隙差几格——逐格写明

```
(i) 面积律（Wilson 圈）      ：一个观测量类的上界；⇏ 聚类（Durhuus–Fröhlich 条件才桥接，
                                cpb-attack/09 "Area law ⟹ gap: conditional, not automatic"）。
(ii) 指数聚类（协方差）      ：全部圆柱观测量的协方差指数衰减 ⇒ 引理 A ⇒ m₁² ≥ m_c²（不需要残数）。
(iii) 质量隙（SZZ 义）       ：= (ii)（SZZ 自己等同二者，§4.3 + §1 逐字）。
(iii′) 质量隙（原子义）      ：(ii) + f₁² ≥ r₀ + I₂ —— 文献无；仓内 P3 给 f₁² 半（强耦合窗）。
```

### §2.4 β=2.5 端缺什么（逐格）

1. **无任何聚类定理**：`β_W ≈ 2.5` 处连面积律都没有严格版（最好 `[R-corpus]` 12/17≈0.706，且归一化未核；crossover SU(2) β_c=2.2986（`BETA-CONVENTION.md:102`））；KP/LSI 机器的收敛判据（`K_S>0`、KP `ū<1`）在 crossover 外全部失效。
2. **缺的那一格 = 跨 crossover 的非微扰搬运**：把强耦合聚合/随机分析工具从 `β_A<1/48` 搬到 `β_W≈2.5`，需要控制一个非微扰 analyticity/single-phase 陈述——即 Clay 墙本身（`CONTINUUM-CP.md:5`："Not proved, any β≳0.71: nothing — no single-phase theorem"）。
3. **即便搬到**，得到的也只是"无光模"（m₁² ≥ 常数 ≫ L^{−2}），不是原子+间距——`I₂`/连续极限仍开（`PROOF-P3-UNCONDITIONAL.md:56–58`）。

### §2.5 (★) 窗覆盖判定（两种读法的精确算术）`[T]`

`(★)=[1/48, 2.5]` 的写法已被 `BETA-CONVENTION.md:176`（C-1）判为 P0 混写（`β_A` 端点与 `β_W` 端点并列）。两种自洽读法：

| 读法 | (★) 的含义 | 已证强耦合窗 | 交集 | 判 |
|---|---|---|---|---|
| β_A | β_A∈[1/48, 2.5]（上端 β_W=N²·2.5，SU(2)=10，越过 β_c=2.30） | β_A∈(0,1/48) | **∅**（只触 excluded 端点） | (★) 内无已证子窗 |
| β_W | β_W∈[1/48, 2.5]（靶点=2.5_W，语料口径） | β_W∈(0, N²/48) | **[1/48, N²/48)** 非空：SU(2) `[1/48, 1/12)`、SU(3) `[1/48, 3/16)` | (★) 的强耦合子窗**已证** |

换算恒等式 `β_W = N²β_A`（`BETA-CONVENTION.md:64`，由 SZZ 作用量 (1.2) `S=Nβ_A Re ΣTr` 与标准 Wilson 的逐 plaquette 系数匹配精确导出）。**语料的 like-for-like 距离 = 30×**（SU(2)：2.5/(1/12)=30.0；精确头 2.2986/(1/12)=27.6；SU(3)：5.6924/(3/16)=30.4；`BETA-CONVENTION.md:122`）。

---

## §3 反方向自检（若有人已证全窗 ⇒ 升级 T）

**检索执行与结果**：
1. arXiv API（`export.arxiv.org/api/query`，工具 OK）："mass gap" ∧ "lattice gauge"，按提交日降序 25 篇：**0 篇**是 4D SU(N) 严格隙（命中全是 U(1)/Schwinger/2+1 维/数值/其他模型；最近的相关严格结果仍是 't Hooft regime 强耦合族）。
2. 仓内全库扫描（`L4-PREPRINTS.md:148–179`，先轮已做、本轮复核）：三条自称"4D YM 质量隙已证"的预印本全部有具名缺步——`2606.19362`（uniform-in-a 靠 Hypothesis 1 公设）、`2307.00788`（隙=Casimir 标签，无标度）、`2506.00284`（arXiv 管理员撤回）；诚实标本 `2503.15539` 自承 "uniform physical mass gap remain unproved"。
3. 结论：**无升级**。不存在任何已证的 4D SU(N) 全 β 窗（含 β_W=2.5）指数聚类或质量隙。（四项核验：①工具 = arXiv API + 仓内文件，OK；②检索式见上两条；③命中 = 0；④无 TOOL-FAIL。）

---

## §4 检索核验表（AGENTS §1：四项缺一不可）

| # | 工具 | 状态 | 实际检索式/调用 | 命中 | 失败标记 |
|---|---|---|---|---|---|
| R1 | Read（语料） | OK | `HN1-H1.md`（385 行全）· `Q6-H1-STATUS.md`（407 全）· `BETA-CONVENTION.md`（227 全）· `M3B-OSWIN.md`（105 全）· `L4-PREPRINTS.md:100–179`· `cpb-attack/09-area-law/REPORT.md`（36 全）· `PROOF-P3-UNCONDITIONAL.md:1–60` | 各返回正文 | 0 |
| R2 | bash grep（仓内） | OK | ①`"Göpfert\|Gopfert\|Münster\|Muenster\|Munster"` @ research/ docs/ `-l` → 20 文件 ②`"2509.03513\|Faria da Veiga\|O'Carroll"` → 6 文件 ③`"Weisz"` → 15 行 ④`"10.1007\|10.1016…"` @ LIT-CHECK/RG2-MAP/LIT-ERA-04 → 40 行 ⑤`"Borgs"` → 10 行 ⑥`"r₀\|r_0"` @ PROOF-P3 → 8 行 | 见左 | 0 |
| R3 | curl + Crossref API | OK | `works/10.1007/bf01961240`→Göpfert–Mack CMP 82:545 ✓；`works/10.1007/BF01212772`→**肝胆外科论文**（检索摘要的 DOI 错）；`works/10.1007/BF01211762`→Kotecký–Preiss ✓；`works/10.1007/BF01240679`→**A/D 转换器论文**（任务猜的 DOI 错）；`works?query.bibliographic=…` ×3（KP、Göpfert–Mack、Borgs–Seiler、Seiler LNP159→`10.1007/3-540-11559-5` ✓） | 8 次调用全 200 | 0（1 次 429 后重试成功） |
| R4 | curl + INSPIRE API | OK | ①`doi:(10.1016/0003-4916(78)90039-8)`→OS'78 摘要逐字 ✓ ②`doi:(10.1007/BF01208780)`→Borgs–Seiler 题名+作者+摘要 ✓ ③`a "Munster, G." and a "Weisz, P."`→**total=1**（flux tube NPB 180）④`a "Munster, G." and ti strong coupling`→1（glueball mass）⑤`ti "glueball mass in strong coupling"`→1 ⑥`ti "strong coupling expansion"`→**137**（含 Münster PLB 121 (1983) 53；无 Münster–Weisz 收敛题）⑦`a "Borgs, C." and a "Seiler, E."`→1（NPB 215 (1983) 125） | 见左 | 0 |
| R5 | curl + Semantic Scholar API | OK（摘要被出版方抹除） | `DOI:10.1016/0003-4916(78)90039-8`、`DOI:10.1007/bf01208780`（abstract: null——"elided by publisher"） | 摘要 2/2 为 null | 0（非检索失败，改走 INSPIRE 得到摘要） |
| R6 | curl + arXiv API | OK | `id_list=2509.03513`→v2 更新 2026-01-24 + 摘要逐字；`search_query=all:"mass gap" AND all:"lattice gauge"` 降序 25 | 1 + 25 | 0 |
| R7 | curl 下载 SZZ v1 PDF + pymupdf 抽取 | OK | 43 页 104k 字符；关键词 `mass gap`×7、`large class of observables`×2、`Theorem`×43；定位并逐字抄出 §1 转述 OS78 段、§4.3 引理/推论段、Corollary 4.11 全文 | ✓ | 0 |
| R8 | WebFetch | 6 成 5 败 | 成：arxiv abs 2204.12737、html 2204.12737v1（截到 §1.2）。败：ar5iv（空返回）×1；link.springer.com（idp 重定向墙）×2；doi.org→springer（302 跨域）×1；projecteuclid（Incapsula 拦截）×1；arxiv pdf（content-type 不支持）×1 | — | **TOOL-FAIL ×5**（均有替代路：改 curl 直取 PDF/API，未据失败下任何"未查到"结论） |
| R9 | WebSearch | 3 成 2 败 | 成：Göpfert–Mack（给出 DOI **BF01212772，后被 R3 证伪**——搜索引擎 DOI 不可采信的实测案例）、Kotecký–Preiss（给出 **BF01212562，后被 R3 证伪**）、Borgs 12/17（返回**乱题名**的不可靠综述，仅作"存在该数值的转述生态"旁证，不采信内容）。败：超时 ×2（Münster–Weisz query、Borgs query 首试） | — | **TOOL-FAIL ×2（超时）**，均改 API 路完成 |
| R10 | py -3.12（精确计算） | OK（exit 0） | 壳计数恒等式 `(2r+1)⁴−(2r−1)⁴=32r³+8r` + `Σ_x r^{−α}` 在 α∈{2,2.5,3}、L∈{4,…,128} 的精确求和 + 积分估计 `2π²∫₁^L r^{3−α}dr` + squeeze 在 HN1 族上 5 个 `Fraction` 精确值 | 全部通过（数值见 §1.2） | 0 |
| R11 | 临时文件 | OK | `~/gy1_szz.html`（62,718 B）、`~/gy1_szz.pdf`（525,950 B）下载后已删；仓库内零新增文件（本文件除外） | — | 0 |

**TOOL-FAIL 汇总**：WebFetch ×5（墙/拦截/类型），WebSearch 超时 ×2。**每一次都换用了替代工具并完成了检索；没有任何"未查到"结论建立在失败的工具上。** 检索引擎返回的 DOI 三分之二与 Crossref 解析不符——**本轮所有 DOI 一律以 Crossref API 解析为准**。

---

## §5 诚实边界（8 条）

1. **本文件没有证 4D YM 的无光模定理在 β_W≈2.5 成立**。交付的是：三条引理 `[T]` + 强耦合窗上的窗口化 `[T]`（推论 D）。窗口化 T 是合法交付物（任务书 §2.1 明认）。
2. **D2 的两条输入各带窗口**：P3 残数定理 `[T] mod OS/KL`、圆柱场地、`β<β₁=min(1/48,c/N)`（`PROOF-P3-UNCONDITIONAL.md:20–31`）；χ 界 `[T]` 的窗口是"全 N 全 L 固定 a"（`Q6-H1-STATUS.md:176`）——两者都在 β_A 口径。合成结论的窗口是两者的交。
3. **D1 的逐有限-L 版是 `[C]`**：SZZ Corollary 4.11 原文对 tight limit（无穷体积）陈述；有限 L 版同一机制（体积一致 Poincaré + GZ 论证）但**不是原文逐字行**。无穷体积版（`inf supp ν̄_∞ ≥ c_N²`）是 `[T]`（读到推论原文）。
4. **`O=¼tr F² ∈ C∞_cyl` 的判定**：格点上 `O(x)` 是 12 条 link 的光滑有界函数（compact 群 ⇒ 有界），落在 SZZ 的 `C∞_cyl` 类里——这是按 SZZ §4.3 的类定义逐条核对得出的 `[T]`；但 SZZ 未逐点列举例子的观测量清单，故记一条诚实的 caveat：**类判据核对过、清单未列**。
5. **Borgs–Seiler 12/17 未核**：数值只存在于仓内三文件转述，原文付费墙，β 归一化不明 ⇒ 全文件不依赖它；它只出现在 §2 表的 `[R-corpus]` 行。
6. **(★) 的 1/48 归一化歧义是继承债**（`BETA-CONVENTION.md:176` 已判 P0）：本文件在 §2.5 把两种读法都算清，不在歧义上取巧。
7. **引理 A/B/C 是测度论/谱论层面的蕴含**，其非平凡输入（OS 正性 ⇒ 谱表示、体积一致 Poincaré、`f₁²≥r₀`）全部来自已引文献与仓内已证定理；本文件没有发明任何 YM 特有估计，因此**它不越过酸测 F1 的空转线**——空转线由假设 (H-a)/(H-b)/(H-c) 承担，而它们在自由场上恰全部失败（§1.4）。
8. **Clay 0%**：本轮交付不改变 `docs/position-current.html` 的任何百分比；无 QED。

---

## §6 记账

- **本轮新增**：三条引理 + 阈值钉死（α>2）+ 双路窗口化 T（β_A<1/48）+ 文献映射表（10 行，含 3 条 DOI/署名更正与 1 条引用名证伪）+ 反方向自检（0 命中）。
- **累计**：语料侧无文件被修改；全部更正建议（`F.→G. Münster`、`REPORT.md:289` 引用名、DOI 表）登记在本文件，留给 lead 决定是否回写。
- **下一个靶**：(★) 窗内最近未证点 = β_W∈[1/12, 2.5)（SU(2)）的第一格——需要跨 crossover 的任一非微扰陈述（单相/解析/聚类任一）；按 `CONTINUUM-CP.md:5`，β≳0.71（最宽 [R-corpus] 面积律）之外一切未证。

*复现：本轮唯一计算 = §1.2 的壳计数求和与 §1.3 的 squeeze 分数（`py -3.12`，确定性，无随机数）；外部检索全部走 Crossref/INSPIRE/arXiv 公开 API 与 arXiv 公开 PDF。只新建本文件；未 `git commit`；未动既有文件；未跑 MC；未编译 Lean。*

**判决：T（窗口：d=4 SU(N)、固定 a、任意 N、全 L、β_A∈(0,1/48) ⟺ β_W∈(0,N²/48)（SU(2) 1/12、SU(3) 3/16）、载子 O=¼tr F²（有界圆柱类）——无光模定理成立；承重引理：C（`χ_L≤C/K_S` 全 L + P3 `f₁²_L≥r₀(β)` 逐 L ⇒ `m₁²(L)≥r₀K_S/C>0`）与 A（SZZ Corollary 4.11 指数聚类 ⇒ `inf supp ν̄ ≥ c_N²`，无穷体积 `[T]`/逐 L `[C]`）；文献：Shen–Zhu–Zhu CMP 400 (2023) 805 / arXiv:2204.12737（摘要+§4.3 逐字读到）+ Osterwalder–Seiler Ann. Phys. 110 (1978) 440（摘要+SZZ 转述）+ PROOF-P3-UNCONDITIONAL.md:20–31；(★)=[1/48,2.5] 的 [1/48,β₀] 子窗在 β_W 读法下已证（β₀=N²/48），β_A 读法下 (★)∩已证窗=∅（1/48|2.5 并列系 BETA-CONVENTION 已判的 P0 混写）；β_W=2.5 端断——缺的那格 = 跨 crossover 的非微扰搬运（Clay 墙本身），且即便补上也只到"无光模"，原子+间距（I₂）仍另开。**）

# 条件定理：4D 纯 Yang–Mills 的存在性与质量隙

> 2026-11-23 · 在具名输入 `IN-0 … IN-11` 下，Clay 陈述的**条件**证明骨架。
> 依赖：`research/gluing/REPORT.md`（共享构造、接口、GL1–GL9、风险）、`research/handoff-proof/`、`research/os2-rotation/`、`research/glue-{a,b,c,d}/`、`research/defect-bootstrap/`。
> **这不是 Clay 问题的解。** 它把"解"精确地写成"`IN-0 … IN-11` ⇒ 结论"，并逐条标注每个输入的状态。其中 **`IN-10`（`hcert`）** 与 **`IN-11`（`IN-b`）** 是新增的两条显式前提：前者是 Lean 合成链写进签名的量级桥（§8 旧记"第五条声明"），后者是缺陷-bootstrap 路线唯一无法内证、必须以物理事实命名的输入。
> **缺陷-bootstrap 路线的终点是条件定理 + 部分界，不是无条件证明**（§3bis、§3ter）：`(iii)`（4D 纯规范 CFT ⇒ 自由）**无法**被 bootstrap 关闭——4D 无 CFT 唯一性定理。

---

## 0. 这句话的准确含义

Clay 问题要求：**任意**紧单规范群 `G`，存在 ℝ⁴ 上非平凡的量子 Yang–Mills 理论，含质量隙。

本文证明：**若**下列 `IN-0 … IN-11` 全部成立，**则** Clay 陈述成立。十二个输入中：

| 状态 | 数量 | 是哪些 |
|---|---|---|
| 已机器验证或已证 | 2 | `IN-1` 的有界类（但空洞，常数=1）、`IN-8` 的组装引理 |
| 方法层面有已发表结果 | — | 强耦合（Osterwalder–Seiler 1978；SZZ）；但 `IN-4` 的 4D 移植未证 |
| **开放** | **7** | `IN-2`、`IN-3`、`IN-4`（4D 部分）、`IN-5`、`IN-6`、`IN-7` 的定量常数、`IN-9`。**若 `(Q-CONS)` 成立，`IN-2` 并入 `IN-3/IN-4`，独立输入降为 6** |
| **显式假设（量级桥）** | **1** | `IN-10`（`hcert`）——`cert ≤ Δ/Λ` 的语义桥，Lean 签名里点名 |
| **物理输入（非可导出）** | **1** | `IN-11`（`IN-b`）——最轻 bulk 态在 `0⁺⁺/2⁺⁺`，缺陷 bootstrap 证不了 |

**没有一个输入是结论的同义反复**（§4 逐条给出分离模型）。

**粘合引理 9/9 封闭**：九条 GL 全部已证明或已归约为具名输入（§5）——这是"九变一"的全部内容。其中 6 条无条件 ✅、6 条 ✅-蕴含（有两条同时属于两类，按拆分的部分计）。

---

## 1. 记号

全部沿用 `research/gluing/REPORT.md` §1.8，不另立变体。关键几个：

| 符号 | 含义 |
|---|---|
| `G` | 紧、连通、单连通单李群；结构群取单连通覆盖 |
| `Λ_{L,a}` | 格点环面 `(ℤ/Nℤ)⁴`，`N = L/a` |
| `β` | `2N_c/g²`，`N_c := dim f`（定义表示维数） |
| `γ` | `β/N_c = 2/g²`（唯一有物理意义的权重系数） |
| `b₀` | `11 h^∨/(48π²)`；`B := 2b₀ = 11h^∨/(24π²)` |
| `λ_G` | `h^∨ g²`（RT 的不变量） |
| `q(μ)` | `m(μ)/μ`（交接坐标） |
| `cert μ` | 尺度 `μ` 分析给出的**下界证书**（函数） |
| `Δbar` | `Δ/Λ` 或 `√σ/Λ`（**单个**被估的数） |
| OS 编号 | OS0 正则性 · OS1 平移 · **OS2 欧氏旋转** · **OS3 反射正性** · OS4 聚类 |

极限顺序固定为**先 `a→0`（固定 `L`）、后 `L→∞`**，不是联合极限（gluing §1.6）。

---

## 2. 定理

**定理（条件）。** 设 `G` 为任意紧、连通、单连通单李群。假设 `IN-0 … IN-11`（§3）成立——即十个构造/物理输入，外加量级桥 `hcert = IN-10` 与物理前提 `IN-b = IN-11`。则存在 ℝ⁴ 上规范不变的测度 `μ`，使：

1. **公理**：`μ` 满足 OS0–OS4，即正则性、平移不变、**O(4) 不变**、反射正性、聚类；
2. **非平凡**：`μ` 非高斯，且 `g_ren/Λ > 0`；
3. **质量隙**：`Δ/Λ ≥ min(δ_W, δ_S) > 0`，且存在 `m ≥ c·Δ`（`c>0` 绝对常数）与 `C>0`，使对所有规范不变局部可观测 `f, g` 与所有 `a ∈ ℝ⁴`，
   ```
   |S(f τ_a g) − S(f) S(g)| ≤ C e^{−m‖a‖} ;
   ```
4. **唯一真空**：`μ` 的 `L→∞` 极限唯一，GNS 真空一维；
5. **禁闭**：`μ` 满足面积律，弦张量 `σ = c_σ Λ² > 0`；
6. **群一致**：所有常数 `m, σ, C, δ_W, δ_S` 只依赖 `(h^∨, dim G, dim f)`。

**结构**：`μ = lim_{L→∞} μ_L`，其中 `μ_L = lim_{a→0} μ_{L,a,β(a)}` 沿同一条重整化轨迹（`IN-0` 的单条 RT，固定 scheme 常数 `C`）。

**转移矩阵形式（定理的精确内核）。** 设 `T_a = e^{−H_a}` 为格点规范不变 Hilbert 空间 `ℋ_a` 上的转移矩阵，`Ω` 为真空。结论 (1)+(3) 等价于：`H_a` 在真空扇区 `1^⊥ = Ω^⊥` 上的谱有**关于 `a` 一致的正隙**

```
∃ δ := min(δ_W, δ_S) > 0 :   inf spec( H_a|_{1^⊥} ) ≥ δ · Λ   (所有足够小的 a) ,
```

且连续统极限（沿 `IN-0` 的 RT、`a→0`）的转移矩阵 `T` 满足 `inf spec H|_{1^⊥} = Δ ≥ δ·Λ`，次大本征值 `e^{−Δ} < 1`（即 `T` 的最大本征值 `1` 与其余谱严格分离），从而给出 (3) 的指数聚类。把格点隙搬到连续统需要 `IN-9`（`SPEC-RES`）；把 bulk-channel bootstrap 界搬到 `Δ` 需要 `IN-10`（`hcert`）**与** `IN-11`（`IN-b`）——这正是本节条件性的全部来源。

---

## 3. 具名输入 `IN-0 … IN-11`

每个输入给出：陈述、状态、来源、**证伪测试**（具体失败模式）。

### IN-0 基础构造（2026 构造的 16 个输入，去掉 `sym` 字段）

- **陈述**：16 个外��输入成立，给出：应力张量（对称、守恒、Ward 恒等式）、迹反常、一圈 running `b₀ = 11h^∨/(48π²)`、子列无穷体积极限、OS 重构框架、有限格点桥。
- **状态**：已认证外部输入（Gross 委员会，2026-11）。
- **关键限定**：**去掉 `Selection.Certificate` 的 `sym` 字段**（`vendor/YM-audit/submissions/openai-20261028/submission/YM/Phase/Mixture.lean:136–148`）。该字段是旋转不变性的**假设**，由 `IN-2` 取代。`checks/SourceScalar.lean:34` 只证"字段 ⇒ E1"，不证字段本身。
- **不在其中**：质量隙、唯一性、禁闭、强耦合控制（`research/handoff-audit/16-input-audit.md`）。

### IN-1 UV 极限（固定体积）

- **陈述**：对每个 `L`，沿 RT 存在子列 `a_n→0` 使 `μ_{L,a_n} ⇒ μ_L`；`μ_L` 满足 OS0/OS1、OS3（反射正性）、规范不变。
- **状态**：**GL3 = ✅-蕴含**（消费 `UVS ∧ UVS-L`；Banach–Alaoglu + 对角 + Minlos），有界类（常数=1）标为**空洞、非结果**。**GL4 = ✅ 无条件（作为引理）**，已机器验证；但作为组合步骤消费弱极限存在（= `UVS`）。反射正性经 GL2a + GL4，`(U_a)` 必需（Gaussian 反例给出锐性）。
- **证伪测试（F3-5）**：无质量 GFF 在 `T⁴_L` 上取 `f=1_{B(0,L/4)}/|B|` 有 `S^{(2)}_L(f,f) ~ L²` ⟹ `UVS-L` 不自动；或 `(U_a)` 在某 `a` 序列失效。

### IN-2 O(4) 不变性

- **陈述**：每个 `μ_L` 在 `O(4)` 下不变，定义为角度荷湮灭真空：
  ```
  (NG)   Q_g Ω = 0      对所有 g ∈ so(4)，
  Q_g = ∫d³x (x_i T_{4j} − x_j T_{4i}) .
  ```
- **状态**：**开放；若 `(Q-CONS)` 成立则不是独立输入**。这是 2026 构造里被假设的 `sym` 字段（R34）。见下面的条件反转。
- **为什么不能从 IN-0 直接推出**：各向异性高斯满足平移、`W`-不变、RP、Ward 却破坏 `O(4)`（os2 `ROUTE-D.md`）。
- **精确结构**（os2 `ROUTE-D.md` §13）：
  ```
  定理 S： Q 反自伴、[H,Q]=0 于 Ω 的定义域、ker H = ℂΩ  ⟹  QΩ = 0.   （三行，不用 Goldstone）
  但局部 Ward 只给：
      [H,Q_g] = i B_g ,    B_g = lim_{R→∞} ∮_{|x|=R} x_{[i} T_{j]k} n_k .
  ```
  于是有一条**条件反转**：
  ```
  MG（⇒ 聚类）+ 守恒对称应力张量 + (Q-CONS: B_g 在稠密定义域上为零)
      ⟹  (NG)  ⟹  U2 .
  ```
  若 `(Q-CONS)` 成立，`U2` 就是质量隙的推论，`IN-2` 不再是独立输入。
- **最锐的检验（已做：`REVERSAL-OPEN`）**：有质量 ℓ=4 高斯
  ```
  Ĉ(p) = 1/(p² + m² + (ε/Λ²)Σ_i p_i⁴)
  ```
  它有隙、`W`-不变、RP、聚类，但破坏 `U2`。计算结果（`research/os2-rotation/check-d14.py`，四个一般动量）：
  ```
  T_{μν} = ∂_μφ∂_νφ − λ∂_νφ∂_μ³φ − δ_{μν}L
  反对称部分  A(p) = λ p_μ p_ν (p_μ² − p_ν²) ≠ 0   ⟹ T 不对称；
  改进 B（deg 3，parity-even，反对称于 (ρ,μ)）的基：
      δ_{ρν}p_μp² − δ_{μν}p_ρp²  与  ε_{ρμνσ}p_σp²
  两者经 p_ρ 收缩恒为零  ⟹  改进像 = {0}  ⟹  无对称守恒局部 T。
  ```
  所以这个高斯**不满足 `Conserved.symmetric`**，**不是反例**，反转仍可能。
  **收窄**：任何**自由**的非旋转不变标量理论都没有对称守恒应力张量（角动量不守恒）；所以反例必须是**有相互作用的**理论——而按谱论证，那类理论需要无隙模式，与质量隙矛盾。反转因此可信，但未证。
- **ℓ=4 形式的地位**：`(NG-ℓ4)`（连通两点函数的 ℓ=4 投影为零）是**必要条件，不是充分条件**——`W`-不变函数在 `S³` 上有 `ℓ=4,6,8,…` 多个通道。最锐的中间形式是 `(NG-ℓ4^T)`（应力张量两点函数），因为 `2⊗2 = 0⊕1⊕2⊕3⊕4`，`T–T` 没有 `ℓ≥6`。
- **证伪测试**：见"最锐的检验"。若该高斯同时满足 `Conserved.symmetric + StressLaws`，则 `MG + QA ⇏ U2`。

### IN-3 弱侧类保持（`CPB-W`，含 `RB2/NPB`）

- **陈述**：存在 `λ_w > 0` 与 `δ_W > 0`，使
  ```
  ∀μ>0,  λ_G(μ) ≤ λ_w  ⟹  cert μ ≥ δ_W .
  ```
  交接所需的达到范围是 `λ_w ≳ 7.3`（`q₀=1` 处；`q₀(λ=6)=0.571`）。
  定量核心是两圈界
  ```
  NPB:  |β(λ) − β_{≤2}(λ)| ≤ C₃ λ⁴ ,   0 < λ ≤ 7.3 ,
  ```
  对 UV 截断与体积一致（`A = 1/C_RG`，`B = 51/(121C_RG²)`，`C_RG = 24π²/11`）。
- **状态**：**开放——唯一的物理墙**。Balaban 型 RG 按构造停在 `λ ~ 1`；中间耦合区无严格方法。
- **已证部分**：两圈闭式 `t_up^{2L}(λ) = C_RG/λ − (51/121)ln(1+121C_RG/(51λ))`；`λ=6` 时 `t=2.639`、`q₀=0.571`；`q₀=1` 在 `λ≈7.28`。
- **与 `IN-4` 的共同机制**：类保持的 Gaussian 计权路线需要 `C_cum ≲ 1.6`（标准 Kotecký–Preiss 在 `g²~2` 给 `≥2`）。唯一能压低它的是 **Bianchi 带号抵消 `κ_c`**（见 `IN-4`）。**`κ_c ≈ 1` ⟹ 两条路线同时关闭。**
- **判定（BL-6，`research/bianchi/BianchiLemma.md` §8.5）：现有机制下无已知路线。** 三条方向一致的证据：(i) 强耦合 surface/character 表象中 `κ_c=1` 可精确证（SU(2) 系数全正 ⟹ 各项非负 ⟹ signed=abs）；(ii) 权重局部因子化 ⟹ 符号不随面积反相；(iii) `dF=0` 不提供可加缺陷荷（4D 非阿贝尔纯 YM 无 2D O(N) 式约束/指标相消）。**这是关于路线的负结果，不是 `CPB-W` 为假的证明。** 翻转只有两条路：`q≈1` 实测 `κ_c<1`，或出现新的可加符号/指标结构。
- **证伪测试**：三圈/格点 `t_up` 在 `g²≈2` 处超出允许的 `Δ`；或 `β` 在 `λ<7.3` 有零点或拐点。

### IN-4 强侧 signed gas（`CPB-S` / `CPB-S^σ`）

- **陈述**：存在 `q_c` 与 `δ_S > 0`，使 `q(μ) ≥ q₁ ⟹ cert μ ≥ δ_S`；对弦张量 `σ` 有对应形式。
- **状态**：2D O(N) 有定理；**4D 移植未证**（规范约束 + Bianchi 替代 O(N) 硬约束）。
- **更锐的子命题（Proposition B，未证）**：若缺陷是闭 2-曲面（`dF=0`）、张力 `σ(q) ~ c₀q²`、熵 `N(n) ≤ e^{γn}`，则 `z(q) ≤ C e^{−cq²}`，`q_c ≈ 0.835/√c₀`。**收益完全取决于 `c₀`**：`c₀=1 → λ*≈6.8`；若张力是 `Λ²` 而 `m²`，则 `c₀≈0.02–0.06`，`q_c` 到 3–4，**反而更差**。反例分支（域壁，张力 `~q`）给 `z ~ 1/q`，命题不成立。
- **证伪测试**：signed-gas 活度 `z(q_c) > 1/2`；或缺陷作用量 `S_def/|C|` 线性于 `q`（而非二次）。
- **与 `IN-3` 的共同机制（关键）**：`q²` 缩放在 `q≈1` 处**中性**（`q²≈q`，BL-5 已证为定理），所以 Proposition B 不能救交叉点。唯一有希望的是**固定边界上的带号抵消**
  ```
  |Σ_{∂C=∂C₀} w(C)| ≤ e^{−c₀q|C₀|} κ_c^{|C₀|} ,   κ_c<1
  ```
  `IN-3` 需 `κ_c ≤ 0.8`；`IN-4` 需更强，故 **`IN-4` 过则 `IN-3` 自动过**。若 `κ_c ≈ 1`，两者同时关闭。

### IN-5 重叠（`H3` × 2）

- **陈述**：`q₁ ≤ q₀ = q(μ_up)`，且两个 certificate（`m` 与 `σ`）**各自**成立（GL7）。数值上等价于 `λ* ≳ 7.3`。
- **状态**：由 `IN-3` 的达到范围与 `IN-4` 的阈值共同决定；两个 `H3` 要各证一次（常数与 `μ*` 不同）。
- **证伪测试**：`q₁ > q₀`（两域不重叠，`cert` 可为零）。

### IN-6 IR 极限与唯一真空

- **陈述**：`{μ_L}` 有唯一弱极限 `μ`；`μ` 满足 OS4。破环方式：MG 表述为**有限体积**族的一致性（gluing §5.1），不是关于极限的命题。
- **状态**：条件于 `W1_0`（中性扇区边界屏蔽）与 `W1_1`（弦扇区）。
- **证伪测试**：存在两个不同的极限点（例如不同 't Hooft 扇区）。

### IN-7 非平凡性

- **陈述**：`μ` 非高斯；`g_ren/Λ > 0`。
- **状态**：由迹反常 + `IN-0` 给出（I-NT）；须用六点或 `g_ren`，**不能用四点**（`κ₄=0` 但 `κ₆≠0`，机器核验）。
- **证伪测试**：四点与六点同时为零。

### IN-8 群一致常数

- **陈述**：`IN-1 – IN-7` 的常数只依赖 `(h^∨, dim G, dim f)`；`min(δ_W,δ_S)` 对每个 `G` 严格正。
- **状态**：组装引理**已证**（GL8 ✅\*，`GroupUniform.lean`）；逐群四常数**未供给**（属 `IN-3`/`IN-4`）。
- **已测**：固定 `g²=1` 时所有群在一圈下失败（`q₀`：SU(2) `1.69e−4`、SU(3) `6.11e−3`、G2 `0.037`）；固定 `λ*=10.4` 时全通过，但需 `g*² = 5.2/3.47/2.60`。
- **证伪测试**：某 `G` 的 `q₀(G) < 1` 且无法由更强耦合补偿。
- **GL8 的拆分**：组装（无条件 ✅）+ 供给（✅-蕴含，消费 `IN-3 ∧ IN-4`）+ **类一致 `GL8-U`（条件于 `IN-8`，不由前两者推出）**。反例 `1/(n+1)`：逐成员正不蕴含一致正；Clay 的"任意 G"要覆盖无穷族（`SU(N)` 全部 `N`），所以要么逐 `N` 证明，要么给出一致界。

### IN-9 谱预解收敛（`SPEC-RES`）

- **陈述**：转移矩阵 `T_a` 在 `1^⊥` 上的谱集在 Hausdorff 意义下收敛（强预解收敛），从而
  ```
  Δ_cont ≥ liminf_{a→0} Δ_a .
  ```
- **状态**：**开放**。这是 GL2b′ 的具名输入。原 GL2b 的朴素方向（`Δ_coarse ≥ δ ⟹ Δ_fine ≥ δ`）**为假**（Cauchy 交错只给上界方向），已由 GL2b′ 替换。
- **为什么不能从 GL3 推出**：glue-b 给出了必要性反例。
- **证伪测试**：格点隙在 `a→0` 时消失、而连续统隙不消失（或反之）。

### IN-10 量级桥（`hcert`）

- **陈述**：`hcert : ∀μ, cert μ ≤ Δ/Λ`。`cert μ` 是尺度分析给出的谱隙**下界**证书（§1 记号）；`mass_gap_of_certificate` 消费它得到 `Δ/Λ ≥ min(δ_W,δ_S)`。无它时 `IN-3/IN-4` 只证"证书为正"，不证"真实谱隙为正"。
- **状态**：**显式假设**（`research/minimal-handoff/MINIMAL.md` §6）。它写进 Lean 签名 `Handoff.lean:77`，不在 `IN-0…IN-9`、也不在 §5 的四件声明里，因此单列为 `IN-10`（旧 §8 记为"第五条声明 / 或 `IN-10`"）。
- **与缺陷 bootstrap 的关系**：bootstrap 供给的是 `0⁺⁺/2⁺⁺` 两个 **bulk channel** 的谱维下界；`cert` 取这两个 channel 界的 `min`。从"channel 下界"到"真实 `Δ` 下界"的对齐，正是 `hcert` 那一步（配合 `IN-b`）。
- **证伪测试**：存在一处 certificate 值严格大于真实 `Δ/Λ`（证书不是下界），或某理论 `cert>0` 而 `Δ=0`。

### IN-11 最轻 bulk 胶球在 `0⁺⁺/2⁺⁺`（`IN-b`）

- **陈述**：最轻 bulk 态落在 `0⁺⁺`（标量，由 `tr F²` 创建）或 `2⁺⁺`（张量，由 `T_{μν}` 创建）这两个 C-even、P-even 扇区；质量隙由这两个 channel 饱和。等价弱形式 `(b′)`：在 `δ` 之上**不存在**带其他量子数的 bulk 态（`0⁻⁺`、奇 C、更高自旋均不低于 `δ`）。
- **状态**：**物理输入，不可由 bootstrap 导出**。格点胶球谱支持：Morningstar–Peardon 1999（`hep-lat/9901004`，淬火 SU(3)）得 `0⁺⁺` 1710 MeV 最轻、`2⁺⁺` 2390 MeV 次之（与 `0⁺⁺` 之比 1.39）。**非对称性保护**：迹反常只把 `tr F²` 绑到自旋-0 **扇区**（`T^μ_μ = (β/2g)·tr F² + Euler`，Euler 为全导数、不产生态），但**不**保证最轻态是自旋-0；无 Ward 恒等式/anomaly 固定最轻 channel。缺陷 bootstrap 只 bound 它看得见的两个 C-even、P-even 扇区（`0⁺⁺` = `tr F²`、`2⁺⁺` = `T_{μν}`），对其他通道无信息——所以 `(b′)` 必须作为独立前提。
- **证伪测试**：格点或连续统出现 `δ` 以下、属于非 `0⁺⁺/2⁺⁺` 的 bulk 态（更轻的 `0⁻⁺`、奇 C、或更高自旋）；或某纯规范 CFT 的最轻态落在缺陷不可见通道，使 bootstrap 之 bound 无法转移。

---

## 3bis. 缺陷-bootstrap 机器（诚实摘要）

> 缺陷-bootstrap（DB）路线（`research/defect-bootstrap/ROUTE.md`）的设想：1-形式 `Z_N` 对称性不作用在局域算符上（`tr F²` 中性），只作用在线缺陷；故要对含 `Z_N`-荷线缺陷的谱做水平 bootstrap，才能排除"相互作用的纯规范 4D CFT"。本路线为 `IN-3/IN-4` 的证书提供一条非 RG 的候选来源。以下为机器摘要，**全部未跑通代码**，边界随条目标注。

**(1) 分块形状（BLOCK-CORRECTED 已定）.** 线在 4D 为 codim-3，保留群 `SO(2,1)×SO(3)`。bulk–defect 块形如

```
𝔽_j(u, θ) = f(u) · P_j(cos θ) ,   P_j = SO(3) 横向分波 ,
```

含**一个** `K_z`-不变 H² 交比 `u` 加 `cos θ`。共形 Casimir 是 `SO(2,1)` 二次 Casimir，即 `u` 上的**超几何型**算子（在 `u=0,1` 正则奇异），**不是**朴素 H² Laplacian `(z²−1)f''+zf'`。`Δ_φ` 只经 `u→0`（bulk-OPE）归一化与动力学系数 `a_O²` 进入；缺陷 channel 形状与 `Δ_φ` 无关。**缺陷道块：已解出**（`BLOCK-FINAL.md` 裁決，两路推导 + 独立数值验证）——`f(u) = (1−u)^{Δ_j}·₂F₁(Δ_j−Δ_φ, Δ_j−Δ_φ; 2Δ_j; 1−u)/N`，`N=Γ(2Δ_j)Γ(2Δ_φ)/Γ(Δ_j+Δ_φ)²`，满足其控制方程 `u(u−1)²f''+(1−2Δ_φ)(u−1)²f'+Δ_φ²(u−1)f=Δ_j(Δ_j−1)f`；指标 `{Δ_j,1−Δ_j}`（u=1，primary+shadow）、`{0,2Δ_φ}`（u=0，identity+double-trace）直接编码两端物理；数值验证 `f(0)→1`、`f/(1−u)^{Δ_j}→1/N` 常数、ODE 残差 ~1e−7。先前"单块无法同满足两端"只在**自由 KG 调和拉普拉斯**（错算子）上成立。**bulk 道块**仍为 Dolan–Osborn import（缺陷运动学归约未验）。

**(2) crossing / junction 的 SDP 形式（SDPB-FORM、SDPB-LAYER）.** 关键更正：**未扭曲**的 1-形式反常**不是新求解器**——它是**标准 SDPB + 每算符相位旋转 + 一条线性等式**。反常由 `ω(q_i,q_j)=e^{2πi k q_i q_j/N}` 为每个缺陷算符固定一个相位 `φ_i`：

```
max  Λ(z)
s.t. Re(e^{−2iφ_i} Λ(F_i)) ≥ 0       对每个允许的 (Δ_i, j_i)   [ray positivity，非新正锥]
     Λ(G_O) ≥ 0                       对每个 bulk O              [实]
     Re Λ(crossing) = 0 ,  Im Λ(crossing) = 0   每个 Z_N 扇区 q
     Λ(F_identity) = 1
```

实现为预处理（`F_i ↦ e^{−2iφ_i}F_i` 的列相位旋转 + 多一条 `Im Λ(cross)=0`）。代价 `~N ×` 标准 SDPB（每个 `Z_N` 扇区一条 crossing）；SU(3) 约 `3×`。**真正新的**是**扭曲（变线）扇区**——两线/域壁 SDP，**在标准 SDPB 之外**。

**(3) SDPB 输入组装（SDPB-INPUT）.** `SO(3)` 精确（`P_j`），把 crossing 投到每个横向分波 `j`，问题拆成 `(j,q)` 上独立的 1D 问题。节点 `N_u = 60`、`w = log((1−u)/u)` 上的映射-Chebyshev，`u ∈ [10⁻⁴, 1−10⁻⁴]`；`j = 0…J`（`J≈10`）；`q ∈ Z_N`（SU(3): 3 扇区）。缺陷块 `F_{Δ,j}(u_k)=f_Δ(u_k)`（一列一个允许 `(Δ,j)`）；bulk 块 `G_{O,j}(u_k)` = Dolan–Osborn 块的 `j` 分波（`O` = identity、`T_{μν}`(ℓ=2)、`tr F²`(ℓ=0)）。逐 `(q,j)` 的 crossing `Σ_{i:q_i=q} a_i² F_i = Σ_O C_O² G_O`，junction 相位在投影**前**乘 `a_i² → e^{2πi k q_i q_j/N}a_i²`。目标：对固定试探块 `max Λ(F_test)`（标准排除测试）。**未验证弱点**：bulk 分波 `G_{O,j}` 依赖缺陷运动学约定（`z,z̄↔u` 映射是弱点）；`φ_i` 依赖 level `k`/BF 归一化；`u`-节点须分辨 `u=1` 附近的 `(1−u)^{Δ_j}`。**未跑代码。**

**(4) 更底层的缺口（RESOLUTION）.** 真正的 bulk–defect 块有 **3 个交比**（codim-3），而现有文献里闭式超几何块（2501.06900、2208.11715、1806.01862）是**缺陷 4-点**块（一个交比，**不同的对象**）。非 SUSY 纯规范 `Z_N` 设置**尚未建立**。故即便形状已定，3-交比 bulk–defect 块本身**缺失**。缺陷 4-点 cross 里含一个 bulk channel，可间接约束 bulk 谱，但更弱、更间接。

**净结论**：DB 的机器**部分建成**——形状已定；**缺陷道显式块已解出并数值验证**（`BLOCK-FINAL.md`）；SDP 形式确认可复用标准 SDPB（未扭曲 = 标准 SDPB + 相位旋转 + 一条线性等式）；SDPB 管线已在 1D GFF 玩具上端到端验证（`SDPB-DERISK.md`，`Δ_min=2Δ_φ` 复原）。**仍缺**：`Z_N` 纯规范设置、3-交比 bulk–defect 块、真 PSD/SOS 求解器及实际应用。它至多给出 `0⁺⁺/2⁺⁺` 两个 bulk channel 的**下界**，方向恰是 `IN-11`（`IN-b`）所需的那一步。

---

## 3ter. `(iii)` 的精确降级：bootstrap 关不掉唯一性

> 这是本路线诚实端点的核心（`III-UNIQUENESS.md`）。**判定：NO。**

DB 路线原想证 `(iii)`：*纯规范、无物质、无 SUSY 的 4D CFT ⇒ 自由*。**它证不了。** 理由分三点：

1. **不存在 4D CFT 唯一性定理。** 4D bootstrap 产出的是**岛**（`N=4`、应力张量 bootstrap）与**解析界**（extremal functional），但**没有任何 4D CFT 被以这种方式证成唯一**。离散 `Z_N` 对称性是**选择定则**（修剪谱），**不是**唯一性引擎。先例：无。
2. **SDPB 只给逐点界，不给全局排除。** 一个 SDPB 岛是对谱空间每个采样点的**严格界**（区间算术），**不是**唯一性陈述。排除**所有**相互作用 CFT 需要下列之一：(i) **解析 extremal functional**——被自由解饱和、被每个相互作用解违反的闭式不等式（罕见，仅存在于解析/large-`N` 极限）；(ii) **先验紧致性归约**——用 gap 界 + 截断把谱压进有限维紧集，再做严格区间扫描（归约本身是解析输入）；(iii) **bootstrap+定理混合**——证明岛边界塌缩到自由点。在采样 `(Δ_i,j_i)` 上的有限扫描只给**逐点不可行**，**永不**给连续统上的全局排除。**没有"多撒 SDPB 点"通往唯一性的路。**
3. **条件于 `IN-b` 时，方向是反的。** 若已有 bulk gap（`IN-b`），bootstrap 反而去 bound **缺陷谱**——这是证 gap 的**错误方向**。

**能诚实声称的最强陈述（部分排除）：**

> 不存在具有纯规范 `Z_N` 结构、且 channel `X` 中最轻算符低于某个**可计算** `δ_bootstrap` 的相互作用 CFT。

这**严格弱于** `(iii)`：`δ_bootstrap` **之上**仍有**相互作用窗口存活**。因此本路线的产出是

> **条件定理（消费 `IN-0 … IN-11`）+ 部分排除界 `δ_bootstrap`。**

**它绝不是 `(iii)` 的无条件证明，也绝不是质量隙的无条件证明。** 翻转条件：出现一个被自由解饱和、被所有相互作用解违反的解析 extremal functional；或一个已证的紧致性归约 + 严格扫描；或一条 4D CFT 唯一性定理。三者目前都**不存在**。

---

## 4. 非空洞性：每个输入都不可删

| 删掉 | 会发生什么 | 分离模型 |
|---|---|---|
| `IN-2` | `μ` 可能各向异性 | 各向异性高斯：满足其余 OS 条件，破坏 `O(4)` |
| `IN-3` | 弱侧无证书，`q₀` 不可控 | 强耦合只覆盖 `λ` 大；`λ≈2` 处无信息 |
| `IN-4` | 强侧无证书 | 弱侧只覆盖 `λ≲1`；`q₁` 无法达到 |
| `IN-5` | 两域不交，`cert` 可为零 | `q₁ > q₀` 时覆盖论证失效 |
| `IN-6` | 极限不唯一 | 两个不同 't Hooft 扇区 |
| `IN-7` | 结论平凡化 | 自由高斯：满足 OS0–OS4 且有隙，但不是 YM |
| `IN-8` | 只在单个 `G` 成立 | Clay 要求"任意紧单群"（无穷族） |
| `IN-9` | 格点隙无法搬到连续统 | 谱集在 `a→0` 时不收敛，或极限谱更密 |
| `IN-10`（`hcert`） | `cert` 不再绑定真实谱隙，`cert ≤ Δ/Λ` 无依据 | 某理论里 `cert > 0` 但 `Δ = 0`：合成链给不出隙 |
| `IN-11`（`IN-b`） | bootstrap 只 bound `0⁺⁺/2⁺⁺`，gap 若在其他 channel 则界不转移 | bulk gap 落在缺陷不可见 / 非 C-even-P-even 的通道 |

**结论没有被假设**：定理的 (1)–(6) 没有一条出现在 `IN-0 – IN-9` 的陈述里；`IN-3`/`IN-4` 只给 `cert` 的下界，`cert ≤ Δbar` 由语义桥 `hcert`（`IN-10`）单独承担（`mass_gap_of_certificate`）；从 bulk-channel 下界到真实 `Δ` 的转移由 `IN-11`（`IN-b`）承担（§3bis）。

---

## 5. 证明骨架

| 步 | 内容 | 依据（GL 分类） | 状态 |
|---|---|---|---|
| S1 | 每个 `L` 存在 `μ_L`，规范不变，OS0/OS1/OS3 | **GL3 ✅-蕴含**（`UVS∧UVS-L`；Lean 化到 Minlos 前，`Minlos` 为显式假设）、**GL4 ✅**（引理；组合处消费弱极限）、**GL6 ✅**（有限层机器验证）、GL2a ✅ | 有界类空洞；无界类消费 `IN-1` |
| S2 | `μ_L` 的 `O(4)` 不变性 | **GL5-KIN ✅**（dim-4 无各向异性）+ **GL5-DYN ✅-蕴含**（`IN-2`） | 条件于 `IN-2` |
| S3 | `∀μ, cert μ ≥ min(δ_W,δ_S)` | `handoff_gap`（Lean） | **已证**（条件于 `IN-3/4/5`） |
| S4 | `Δ/Λ ≥ min(δ_W,δ_S)` | `mass_gap_of_certificate`（Lean） | **已证**（条件于 `hcert`） |
| S5 | 一致隙 ⇒ 关于 `L` 一致指数聚类 | **GL9 ✅-蕴含**（`CPB-W`），`GlueA.spectral_decay`（Lean） | **已证（作为蕴含）** |
| S6 | 隙跨正则化传递：`Δ_cont ≥ liminf Δ_a` | **GL2b′ ✅-蕴含**（`IN-9`），4 行证明 | 条件于 `IN-9` |
| S7 | `{μ_L}` 有唯一极限 `μ`，继承 OS0/1/2/3 | `IN-6`、GL4 | 条件 |
| S8 | `μ` 满足 OS4 且隙保持 | S5 + S6 + GL4 | 条件于 `IN-6` |
| S9 | `μ` 非高斯 | `IN-7`、迹反常 | 条件 |
| S10 | 面积律 `σ>0` | `IN-4`（σ 版）、**GL7 ✅-蕴含**（`CPB-W ∧ CPB-S^σ`），`GL7Dual.lean` | 条件 |
| S11 | 群一致常数对每个 `G` 为正 | **GL8 ✅**（组装）+ **GL8-S ✅-蕴含**（`IN-3∧IN-4`）；`GL8-U` 条件于 `IN-8` | 条件 |
| S12 | OS0–OS4 + 隙 ⇒ Wightman 理论 + 谱隙 | I-RECON | 标准 |

**缺陷-bootstrap 在合成链中的位置**：S3 的 `cert ≥ min(δ_W,δ_S)` 可由 DB 在 `0⁺⁺`（`tr F²`）与 `2⁺⁺`（`T_{μν}`）两个 **bulk channel** 的谱维下界来供给（§3bis，作为 `IN-3/IN-4` 的非 RG 候选来源）。但把"channel 下界"转成"真实 `Δ` 下界"**恰好**需要 `IN-11`（`IN-b`）：若最轻 bulk 态不在这两个 channel，界就不转移。随后 `IN-10`（`hcert`）再把 `cert` 接到 `Δ/Λ`。**故 DB 给的是 `IN-3/IN-4` 的证书候选，不是它们本身；它无法替代 `IN-10/IN-11`。**

**GL 架构：9/9 封闭** —— 每条要么**无条件 ✅**（GL1、GL2a、GL4、GL6、GL5-KIN、GL8-A），要么**✅-蕴含**（GL3、GL2b′、GL5-DYN、GL8-S、GL9、GL7），消费的输入在陈述里点名。**9/9 的含义是"每条粘合引理已证明或已归约为具名输入"，不是"九个问题都解决了"。**

**机器化状态（11 块）**：**9 块已 Lean 化**（GL2a、GL2b′、GL3、GL4、GL5-KIN、GL6、GL7、GL8、GL9 —— 其中 GL3/GL4/GL6 各带一个显式声明）；**GL5-DYN 无可形式化内容**（`IN-2 ⟹ OS2` 是改名/记账）；**GL1 是审计**（跨文档合规，非定理）。

**四件"显式声明而非机器证明"的事项**（都在签名或文件头里点名）：
1. **Minlos**（GL3）：`hMinlos` 写进签名，Mathlib 缺失；
2. **截断的 Borel 可测性**（GL4）：`trunc A n ∈ L²` 作为显式假设；
3. **变分隙 ↔ 谱隙的等同**（GL2b′）：声明性 scope；
4. **Reed–Simon VIII.24b**（GL2b′）：Mathlib 缺该引理。

**"全机器"的准确含义**：每个**可形式化的**陈述都已 Lean 化、无 `sorry`、公理冻结。四件声明是外部假设，不是证明——它们被推到类型层面，正是为了让读者能看见。

**机器验证的边界**：S3、S4、S5、S10、S11 的核心已 Lean 化（§6）；S6 是 4 行散文证明；S1、S2、S7–S9、S12 是散文 + 条件。

---

## 6. 已机器验证的部件

```
bash lean/check-handoff.sh
  → Handoff.lean            handoff_gap, mass_gap_of_certificate
  → RGHandoff.lean          certificate_transfer
  → HandoffMargin.lean      handoff_iff_margin, required_lambda
  → ALL THREE COMPILE       公理冻结为 [propext, Classical.choice, Quot.sound]

bash research/glue-a-limits/check-glue-a.sh
  → GlueA.spectral_decay    （GL9：隙 ⇒ 聚类）
  → GL4Transfer.lean        （GL4：RP 与聚类界的弱极限传递）
  → GL4L2.lean              （GL4：L² 截断闭包，全锥机器验证）
  → GL3Limits.lean          （GL3：紧性 + 可数对角；Minlos 为显式假设）
  → ALL GLUE-A CHECKS PASS

bash research/glue-d-dual-handoff/check-gl7.sh
  → GL7Dual.lean            （GL7：两个 hcert 同轨塌缩）
  → 全部通过

bash research/glue-b-blocking/check.sh
  → GL2a.lean, GL2bPrime.lean （RP 保持分块、隙传输蕴含）
  → GL2a + GL2bPrime COMPILE (no sorry, standard axioms frozen)

bash lean/check-conditional.sh
  → Conditional.lean        （clay_of_inputs：NPB → margin → Δ/Λ > 0）
  → 公理冻结为 [propext, Classical.choice, Quot.sound]

bash research/glue-c-gauge-group/check-glue-c.sh
  → GroupUniform.lean       （GL8 组装引理）
  → GL6Fubini.lean          （GL6：有限轨道 Fubini + S₃ 精确枚举，机器验证）
  → GL1 合规：9/9 接口同一测度、同一拓扑、同一 RT
  → ALL GLUE-C CHECKS PASS

bash research/os2-rotation/check-symdim4.sh
  → SYMDIM4 OK               W₄ 平均 = (tr A/4)·1（dim-4 无各向异性）
```

另有可复现脚本：`verify.py`（GL2a 判据）、`gl6_check.py`（GL6 精确枚举）、`gl1_check.py`（测度合规）、`gl4_gaussian_counterexample.py`（`(U_a)` 的锐性）。

---

## 7. 与 2026 构造的关系

- 2026 构造证明了：**存在性**（子列极限）+ OS 重构 + 应力张量/Ward + 迹反常 + 一圈 running。
- 它**没有**：质量隙、唯一性、禁闭、强耦合控制。
- 它把 `O(4)` 放在 `Selection.Certificate` 的 `sym` 字段里（**假设**）。
- 本定理做的事：把 `sym` 提出来记为 `IN-2`；把其余缺口记为 `IN-3 … IN-9`；把量级桥与物理前提分别记为 `IN-10`（`hcert`）与 `IN-11`（`IN-b`）；证明"输入 ⇒ Clay"。

**因此**：2026 构造 + `IN-2 + … + IN-9` + `IN-10`（`hcert`）+ `IN-11`（`IN-b`）= Clay。这就是本项目的全部内容。**注意**：`IN-11`（`IN-b`）来自缺陷-bootstrap 路线，是该路线唯一无法内证、必须以物理事实命名的输入；该路线本身给出的是条件定理 + 部分界 `δ_bootstrap`（§3bis、§3ter），**不是**无条件证明。

---

## 8. 诚实清单

1. **`IN-3` 是唯一的物理墙。** 它要求中间耦合区（`λ ∈ [1, 7.3]`）的非微扰 `β` 界；Balaban 型 RG 停在 `λ~1`，强耦合展开要求 `λ` 大，中间区无严格方法。**多开并行会话不会改变这一点。**
2. **`IN-2` 是独立输入**，不是 OS 公理的推论（各向异性高斯是分离模型）。
3. **`IN-4` 的 4D 部分未证**；2D 有定理。
4. **`IN-6` 条件于 `W1_0/W1_1`**，两者都未证。
5. **GL2b′ 是 ✅-蕴含**（`SPEC-RES ⟹ Δ_cont ≥ liminf Δ_a`），`SPEC-RES`（`IN-9`）本身开放；`IN-1` 的无界类条件于 `UVS ∧ UVS-L`，`GL4` 的**组合**用法条件于弱极限存在。
6. **S1、S2、S7–S9、S12 是散文**，未形式化；S3、S4、S5、S10、S11 已形式化（S6 是 4 行散文证明）。
7. **群一致是纪律不是定理**：逐群常数未供给。
8bis. **`IN-3` / `IN-4` 在现有机制下无已知路线（BL-6）。** Bianchi 不足；B1 在 `q≈1` 中性（定理）；B2 默认 `κ_c≈1`（判定）。这否定的是"用 Bianchi / 符号抵消证明它们"，**不是**否定它们本身。翻转条件写在 `IN-3` 条目里。

9. **`IN-2` 可能不是独立输入。** 若 `(Q-CONS)`（边界算符 `B_g` 在真空上为零）成立，则 `MG ⇒ U2`，`IN-2` 降为 `IN-3/IN-4` 的推论。有质量 ℓ=4 高斯的检验**已做，结果 `REVERSAL-OPEN`**：它没有对称守恒应力张量，因此不是反例；反例必须来自**有相互作用**的理论。这是"U2 与 MG 逻辑独立"这条结论可能被反转的唯一位置。

10. **`hcert` 现记为 `IN-10`（第五条显式假设，未列入"四件声明"）。** `Handoff.lean` 的 `mass_gap_of_certificate` 以 `hcert : ∀ μ, cert μ ≤ Δ/Λ` 为假设（[Handoff.lean:77](../lean/YMStatement/Handoff.lean:77)）；它不在 `IN-0…IN-9`，也不在 §5 的四件声明里。最小化版本与残余缺口的完整记账见 `research/minimal-handoff/MINIMAL.md` §6。**本文将其正式记为 `IN-10`。**

11. **`IN-11`（`IN-b`）是物理输入，不可由 bootstrap 导出。** 最轻 bulk 态在 `0⁺⁺/2⁺⁺` 只在格点上是事实（Morningstar–Peardon：`0⁺⁺` 1710 MeV 最轻、`2⁺⁺` 2390 MeV），非对称性保护（迹反常不保证最轻态是自旋-0）。DB 只 bound 这两个 C-even、P-even 扇区，对其他通道无信息——`(b′)` 必须独立假设（`IN-B.md`、`GAP-DICTIONARY.md`、`B-CHECK.md`）。

12. **缺陷-bootstrap 证不了 `(iii)`（精确降级，§3ter）。** 4D **无** CFT 唯一性定理；SDPB 逐点给严格界、不给全局排除，**没有"多撒点"通往唯一性的路**。最强诚实陈述是**部分排除**："不存在 channel `X` 中低于可计算 `δ_bootstrap` 的纯规范 `Z_N` 相互作用 CFT"——`δ_bootstrap` **之上仍有相互作用窗口存活**。条件于 `IN-b` 时 bootstrap 反而去 bound 缺陷谱（**错误方向**）。**本路线的产出是条件定理（消费 `IN-0 … IN-11`）+ 部分界 `δ_bootstrap`，不是 `(iii)` 或质量隙的无条件证明。** 详见 `research/defect-bootstrap/III-UNIQUENESS.md` 与 `CONDITIONAL-SUMMARY.md`。

---

## 9. 一句话

> **本项目的产出是：Clay 问题 = 12 个具名输入（`IN-0 … IN-11`）+ 一条已机器验证的合成链（9 条粘合引理全部封闭）。** 其中 7 个构造输入开放、`IN-3` 是已知困难的那一个；`IN-10`（`hcert`）是量级桥、`IN-11`（`IN-b`）是物理输入。缺陷-bootstrap 路线给出的是**条件定理 + 部分排除界 `δ_bootstrap`**，它**不能**关闭 `(iii)`（4D 无 CFT 唯一性定理——§3ter）。

---

## 附录：输入来源索引

| 输入 | 文档 |
|---|---|
| `IN-0` | `docs/yang-mills-prior-work.md` §4；`research/handoff-audit/16-input-audit.md` |
| `IN-1` | `research/glue-a-limits/REPORT.md`；`research/glue-c-gauge-group/REPORT.md` §3 |
| `IN-2` | `research/os2-rotation/NPSYM.md`、`ROUTE-D.md`；`research/os2-rotation/lean/SYMDim4.lean` |
| `IN-3` | `research/handoff-proof/REPORT.md` §0.5、§2.7.1 |
| `IN-4` | `research/handoff-proof/REPORT.md`；`research/confinement/REPORT.md` |
| `IN-5` | `research/handoff-proof/REPORT.md`；`research/glue-d-dual-handoff/REPORT.md` |
| `IN-6` | `research/uniqueness/REPORT.md`、`R32.md` |
| `IN-7` | `research/nontriviality/REPORT.md` |
| `IN-8` | `research/group-uniformity/REPORT.md`；`research/glue-c-gauge-group/GroupUniform.lean` |
| `IN-9` | `research/glue-b-blocking/REPORT.md` §3.8（`SPEC-RES`） |
| `IN-10`（`hcert`） | `lean/YMStatement/Handoff.lean:77`；`research/minimal-handoff/MINIMAL.md` §6 |
| `IN-11`（`IN-b`） | `research/defect-bootstrap/IN-B.md`、`GAP-DICTIONARY.md`、`B-CHECK.md`；lattice `research/mass-spectrum/REPORT.md` §4.5（Morningstar–Peardon `hep-lat/9901004`） |
| DB 机器 | `research/defect-bootstrap/`：`ROUTE.md`、`BLOCK-CORRECTED.md`、`SDPB-FORM.md`、`SDPB-INPUT.md`、`SDPB-LAYER.md`、`RESOLUTION.md`、`III-UNIQUENESS.md`、`CONDITIONAL-SUMMARY.md` |
| 合成 | `research/gluing/REPORT.md`（§1 构造、§2 接口、§3 GL、§3.5 台账） |

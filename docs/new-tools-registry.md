# 新工具登记册 · 战役造出的一切

> 2026-11-23 · 这轮"造工具 / 硬算"总攻造出的**新对象、新 no-go、结构定理**总登记，免得忘。
> 铁律：每个构造都过**自由场酸测**（自由无质量场 `f₁²=0`、`G″=0` 时必须失效，否则 vacuous）。

---

## A. 可用的新构造（过酸测 / 真能碰原子）

### A1. 乘性残数搬运（multiplicative residue transport）
**公式**：`f₁²(μ) = Z_O(μ,μ₀)² · f₁²(μ₀)`，`Z_O = (β/2g)(μ₀)/(β/2g)(μ) > 0`，`γ_O = −d ln(β/2g)/d ln μ = 2b₀g²+…`
- 依据：`O_S=¼tr F²` **乘性重整**（恒等混合对角、`F F̃` 为 P-奇、BRST-恰当导数在物理态为零）。
- 过酸测：`0→0`（乘性）；自由场 `:φ²:` 由 Z₂ 宇称 `a≡0`；常数 `Z_O²>0` 凭"是平方"、`O(1)` in `V`。
- 定理：`f₁²(2.5) ≥ Z_O² r₀(β_eff)`，**条件于 (CP)**（单 T=0 禁闭相 + 0⁺⁺ 边 Kato-保持）。
- 残留：**(CP)**（= 单相 + 隔离）。文件 `INVENT-TRANSPORT-INEQ.md`。

### A2. 精确方差分解（Savvidy-corrected Brascamp–Lieb）
**公式**：`G″(0) = S(v) − D(v)`，
`S = E_μ[h_n − ⟨K⁻¹c,c⟩ | 0]`（高斯–Schur 曲率，`K=Π₀ Hess V Π₀`，`c=Π₀ Hess V v`，`h_n=vᵀ Hess V v`）；
`D = Var_μ(∂_nV | 0) − E_μ[⟨K⁻¹c,c⟩ | 0]`（条件化纤维的 BL 缺陷）。
符号律：**`D≤0 ⇔ K⪰0`（真 BL）；`D=0 ⇔` 高斯；`D>0 ⇔` 纤维 Savvidy**。`c_Φ=N/12`。
- 过酸测（非平凡）：自由场高斯 `D=0`、`k→0` 模 `Var(a)=∞⇒S=0`，`G″=0−0=0`（不造伪正性）。
- 定理：**β=0**：`V=Φ` 凸 ⇒ `D≤0` ⇒ `G″≥S≥N/12>0`。**β=2.5**：`D/S≈2%`（边缘基本高斯）⇒ `G″>0` 数值稳健、裕度 ~98%。
- 残留：一致界住 2% 纤维缺陷 `H_ff`（= H1）。文件 `INVENT-VARIANCE-INEQ.md`。

### A3. Fredholm 行列式恒等式（KL 谱行列式 / Selberg–Gutzwiller 型）
**公式**：`f₁² = −λ₁ Res_{z=z₁} G(z) = −(1/z₁)·N_O(z₁)/D_O′(z₁)`，
`D_O(z)=det(1−z T_O)`，`log D_O = −Σ_k z^k Tr(T_O^k)/k`，`z₁=e^{am₁}`（最小模零点），
**"周期轨道"** `Tr T_O^k = Z_{S¹_k×T³}`（闭合时间通量管）。
- 过酸测（严格）：自由场纯割无极 ⇒ 留数 0；YM 有自己双胶球割下的孤立零点（`m₁²<4m₁²`）⇒ `f₁²>0`。
- **看得见原子**：极点与割低矩全同，唯**零点集**不同——线性盲、行列式可见；且**供给 `m₁²`**（零点位置）。
- 残留：`a→0` 求零点（连续统）= 边际 BPHZ 墙。文件 `INVENT-NEW-IDENTITY.md`。

### A4. 规范不变 IR 质量（= `tr F²` 谱边）
**公式**：`m_GI := inf supp ρ_s = lim_{|x|→∞} −∂_{|x|} ln|G_s(x)| = E_1^{(s)} − E_0`，
在 `H_phys=L²(𝒜/𝒢, μ̄)` 上取 `s=(β/2g)tr F_{μν}F^{μν}`，`ρ_s≥0`（反射正性）。
- **移除伪障碍**：`A∉H_phys`、`ρ_A` 符号不定——都只来自用势 `A`；用 `F²` 则 `ρ_s≥0`。**规范不变性从来不是墙。**
- 过酸测（定义）：自由 Maxwell `m_GI=0`（阈值 0 两光子割），YM `=m_{0⁺⁺}`。
- 残留：边原子 `f₁²>0`。文件 `INVENT-GAUGE-INVARIANT-MASS.md`。

### A5. Gauge W-熵（Perelman 移植）
**公式**：`F_W[μ_t] = ∫_𝒜 (S_W(A)+|∇_A f_t(A)|²_G) e^{−f_t(A)} dvol_G`，沿共轭流
`∂_tA=−∇_A S_W`、`∂_tf=−Δ_A f+|∇_A f|²−S_W`。保隙单调 ⟺ `Γ₂=Ric_𝒜+Hess_𝒜 S_W ≥ ρ`。
- 酸测（镜像过）：自由 `Γ₂≥0` 单调但只认证自由格点隙（非空）；YM `Γ₂≈−2.0<0` —— **在无隙侧对、有隙侧阻**。
- 判词：Perelman 的 W 要"几何=流"同对象；平空间 YM 拆开二者。**要越界的不是新泛函，是那个缺失的非微扰流本身。** 文件 `INVENT-NEW-MONOTONE.md`。

---

## B. 新 no-go 定理（"为什么这类不行"的新理由）

### B1. 周期 no-go（convexification 不可能）
`μ` 的对数密度**周期**（紧群），**周期凹=常数** ⇒ 无对数凹替身等于 `μ`；**Savvidy 负性=周期性在回绕 `|x|>π/2` 的代价**。（热磨光只给半凹；Prékopa 提升被 G″/H_ff 不定禁。）`INVENT-CONVEXIFICATION.md`。

### B2. 嵌套 no-go（第二锥不可能）
每个变号/Karlin 锥 `K₀⊂K₁⊂…` **都含正锥** ⇒ 在每层吞下严格正的自由场 ⇒ vacuous。能过的锥必须切到**正性之下**（符号不定），而那里 Schur-凸**预设了要证的原子**。（另：plaquette 积张量化非符号正则。）`INVENT-SECOND-CONE.md`。

### B3. 刚性 no-go（指标不可能）
任何**离散整数**、且在保对称变形下**局部常数**的指标，沿连通 β 族（`β=∞` 达自由高斯）**必为常数** ⇒ 等于自由值 ⇒ **YM 全为 0**。**隙必须是 β-非解析跳变，不是任何离散不变量。** 统一解释 NCG/localization/BPS/Witten 旋量为何全在自由点返 0。（唯一过酸测的 Gribov-FP 指标靠数规范副本。）`INVENT-NEW-INDEX.md`。

### B4. Feshbach/Schur + 图角 no-go（消去负块只会帮幸存者）
块分解 `H=[[A,B],[Bᵀ,C]]`，`C≺0`（LLL 负块）⇒ `Σ_P=A−BC⁻¹Bᵀ=A+B|C|⁻¹Bᵀ⪰A`：修正**只加不减**（Haynsworth `inertia(H)=inertia(C)+inertia(Σ)`，验证 `λ_min(Σ−A)=−6·10⁻¹⁶`）。故"小角度保 `Σ` 正"是符号错误：`Σ>0⇔A=PᵀHP>0`，`A` 即隙。图角唯一通过重叠 `ε=sin²θ` 进入，需 `ε<½`（`θ<45°`）才有 `A>0`：**冻结常背景上 `θ=0` 但 vacuous（`F≡0` 同 0、`g`-依赖）**，**物理 `k=0` 上 `ε=O(1)`（实测 `PᵀHP` min `−0.3`）、无动量分离**。整条路 `⇔ ε<½ ⇔ H1`。`FORGE-FESHBACH-ANGLE.md`。

### B5. Krein 空间 / Pontryagin no-go（计数≠底，UV≠IR；分叉死结）
Savvidy Hessian `−D²+2adF_B` 是**不定**算子，Krein 空间谱论（Pontryagin Π_κ、Langer definitizable、谱函数）确是其精确语言，且能把 UV 塔锥空严格证成形（`λ_j>k_⋆` 壳 `κ=0`、`Tr H_−≲∫|F|³`）——**真收获（UV 半定理）**。但一切 Krein 不变量是**计数/符号**（负指标 `r_−`、符号差）或 **UV 尺度**（J-正底 `=Θ(g|B)`，Landau 能级），而隙是 `k→0` 底部的 **min**、经量纲透射的 **IR** 尺度。两路互斥封死：`J=sign H`（使地板物理）⇒ `κ=dim H_+=∞` 且循环（需正谱）；先验 BRST/Gupta–Bleuler 度规 ⇒ `κ=O(1)` 但自由场**逐字相同**（酸测空洞）。signed 负指标是整数 ⇒ B3 刚性 ⇒ 恒等自由值 `0`。**definitizing 劈裂恰在 `k=0`（Savvidy 带底、粗模与最负模重合）失去一致裕度——Krein 临界点正是间隙要印数的那一点。** 余下缺失不等式 `inf spec(H|_物理侧)≥κ_{IR}>0`（`(L,a)`-一致）= `(H)` = H1。`FORGE-KREIN-SPECTRAL.md`。

### B6. 全单调 / 解析延拓 no-go（原子是整函数，floor 只在 `t→∞` 速率里）
`C(t)=∫e^{−st}dν(s)` 全单调（Bernstein–Widder）**⟺ `ν≥0`**，与 KL/RP 正性同一件事，**无超出正性的刚性**。故（B1 型线性数据）任何有限 `t` 标量都是 `m₁²` 的**递减上界**：`r(t):=−(1/t)log(C(t)/M₀)` 满足 `m₁²≤r(t)≤M₁/M₀` 且 `r↓m₁²`（`log C` 凸，Cauchy–Schwarz）。**方向反。** 解析孪生**结构性瞎**：边原子贡献 `f₁²e^{−m₁²t}` 是**整函数**（无奇点），而收敛半平面 `Re t>0` 由 UV 延伸（`supp` 到 `+∞`）固定，**与 `m₁²` 无关** ⇒ 半径/延拓/奇异点看不见 IR floor；`t=0` 导数增长=矩比（又上界）且连续统 UV 发散；转移变量 Padé 才把原子看成极点，代价 `M≳2/(am₁)→∞`（`PADE-ROUCHE-DEEP`，relocate 到 `W3`）。Wiener–Ikehara 的"主项"= 领头指数 `m₁²/f₁²`，其正则性假设需**先减去答案**⇒循环。**唯一能印 floor 的是 `t→∞` 速率 = 一致指数簇聚 = 隙（循环）；缺失步仍是搬运（H\*）。** `TOOL-COMPLETELY-MONOTONE.md`。

### B7. 梯度流 no-go（B6 实例：流能 = 完全单调）
Lüscher 梯度流 `∂_tU=−D_νG_{νμ}` 保反射正性；流向能量 `E(t)=⟨tr G_{μν}G_{μν}⟩(t)` 是 Laplace 变换 `∫h(u)e^{−2tu}du`、`h≥0` ⇒ **完全单调 = B6 = `ν≥0`**，无超出正性的动力学内容。流时"有效质量" `m_flow²(t)=−½∂_t ln E(t)` 是 tilt 谱均值，`m₁² ≤ m_flow²(t)` 且 `↓ m₁²` —— **只上界（= susceptibility / H1）**，`t→∞` 才等式（循环）。**酸测败**：自由无质量 `t²E(t) ≡ 1/(64π²) > 0` ⇒ 任何由 `E, tE, t²E, t²∂_tE` 造的下界都被无隙场满足（vacuous）。`k→0` 粗曲率 `G″=1/χ` 是流不变量。缺失步 `H*(r₀)`。`TRANSPORT-GRADIENT-FLOW.md`。

---

## C. 战役中证出的结构定理（真结论，非 no-go）

| 定理 | 内容 | 文件 |
|---|---|---|
| **P3 强耦合残数定理** | `f₁² ≥ c₁(β)u^p(1−O(u)) > 0`，`β<1/48`、体积一致、显式 `β₁=O(1/N)`（模一个标准 graded-strip KP 引理） | `PROOF-P3/P1P2` |
| **无溶解结构定理** | `0⁺⁺` 字面并进 `(gg)` 割**谱上不可能**（`m₁<2m₁`、`0⁺⁺` vs `0⁺⁺*` von Neumann–Wigner 避免交叉） | `PROOF-NO-DISSOLUTION-KINEMATIC` |
| **W_geo≥0 引理** | Haar 几何项非负：`Φ=−logJ=Σ g(λ_i−λ_j)` 凸（`g″≥1/6`）⇒ Prékopa ⇒ `W_geo=E[Φ″]−Var(Φ′)≥0`，全模全 k | `PROOF-JACOBIAN-SAVES` |
| **规范不变 PF 定理** | 规范不变转移矩阵是 Perron–Frobenius（满正 Haar ⇒ 严格正核，无 Gribov/FP 符号）；但只钉真空（与胶球同 `0⁺⁺` 道） | `PROOF-GAUSS-INVARIANT-PF` |
| **Lemma SC** | 强耦合端粗模 Poincaré 成立：SZZ 一致 LSI 过 1-Lipschitz 块映射 `Var_ν(φ)≤(C_blk²/λ_μ)∫|∇φ|²` | `PROOF-COARSE-POINCARÉ` |
| **-0.165 是假阴** | 该值漏 `+W`；同模精确 `G″=+0.51`（SU3 2.5），`G″` 在可测处为正 | `PROOF-JACOBIAN-SAVES` |

---

## D. 残留对象词典（所有路最终都落到这里）

1. **H1**：粗模 Poincaré `Var_ν(φ)≤(1/λ)∫|∇φ|²`，`λ>0` 一致 ⟺ `m_eff²(k→0)≥λ`。（`∇∥` 半边。）
2. **(CP)**：单 T=0 禁闭相（无体相变）+ 0⁺⁺ 边 Kato-保持 + 一致隔离。（搬运定理所缺。）
3. **边原子 `f₁²>0`**：`ν({m₁²})>0`（usc、非对角、需 `E₂`）。
4. **`G″(0)>0` at `k→0`**：粗模曲率，`=S−D`，`D/S≈2%`。
5. **行列式零点 in 连续统**：`D_O(z)=det(1−zT_O)` 的零点在 `a→0` 的求值。
6. **Savvidy 符号**：非阿贝尔/Weingarten/`6j`/`H_ff` 不定——一个物的五个面。
7. **一切路 = 隙**：H1、NEO、I2、束缚能、`(CP)`、边原子……皆与质量间隙同构。

---

## F. 待审历史独立路线（未充分审查，勿当已解 / 勿当已死）

### F1. Battle–Federbush 相位单元 / Federbush chunk 展开（MZP 候选）
**来源（已逐一核实，2026-11-24，宪章 §2）**：
- Battle–Federbush, "A phase cell cluster expansion for Euclidean field theories", *Ann. Phys.* **142**(1) 95–139 (1982), DOI 10.1016/0003-4916(82)90230-5 —— **mass-zero** 簇展开、相位单元系数为节点、关联函数=有限耦合模式子集、`d=2,3` 有 vacuum-stability 替代、**`d=4` 明列为未来目标**。
- Federbush, "Phase cell approach to YM: IV. The choice of variables", *CMP* **114** 317–343 (1988), DOI 10.1007/BF01225039 —— 三类变量：pure modes / averaging-correction / **chunks**（大场、处理 `π₃(G)` 拓扑障碍、点状规范奇点）。
- Federbush, "V. Analysis of a chunk", *CMP* **127** 433–457 (1990), DOI 10.1007/BF02104497。
- 近亲：Magnen–Rivasseau–Sénéor, "Construction of YM₄ with an IR cutoff", *CMP* **155**(2) 325–383 (1993), DOI 10.1007/BF02097397（纯 SU(2)、**固定红外截断**无 UV 截断、轴向规范大场正性）。
- 背景：Douglas, "The Yang–Mills Millennium problem", *Nature Rev. Phys.* **8**(2) 86–97 (2026), DOI 10.1038/s42254-025-00909-2（明说"no clear path known"）。

**结构（与已死物种的真实不同处）**：簇展开在**模空间**（相位单元系数）而非时空区域；**mass-zero**（不预插质量传播子喂进展开）；chunks 作大场**真积分变量**、不被 Bernoulli 坏块概率排除；小模对大场背景经 **`F` 而非 `A`** 耦合、纯规范背景解耦、带 `(ℓ/L)²` 因子。

**为何是盲区**：免疫图有 Balaban、MRS，却漏了这条 —— 值得正式在册。

**已知弱点（预判，子代理正在证/伪）**：`(ℓ/L)²` 是尺度分离工具 —— 在 UV→IR（`aΛ→0`，即我们的 G1-L3）有效；但质量在 `1/Λ` 处生成，那里软模 `ℓ~L~1/Λ`，`(ℓ/L)²→1`，**抑制恰在生成处消失**。故只有两条路：要么模空间森林靠**组合学**在 `ℓ/L→1` 硬控（=Federbush 当年 `d=4` 未证成处），要么仍需一个新动力学输入（=H1）。`(II) chunk-induced IR screening` 作者自认**未证、不预设符号**。

**状态**：**[已审 · 2026-11-24] `[NG]`（G2 死）+ `[T]`（G1-UV 真定理）** —— 见 `research/unconditional/TRANSPORT-PHASE-CELL.md`（酸测 1–4 + 与 H1 逐项对齐，witness `w7/transport_phase_cell.py`）。判词：`(ℓ/L)²` 是尺度分离因子，在 `ℓ~L` **→O(1)**；质量由**尺度无关且不定**的 Savvidy 二次耦合给出；chunk 积分是**正测度（只给上界）**；森林活性在交叉点为边际 `g²=O(1)`。**(II) = H1/I2 的描述重写；MZP-1 = L3/H1 同一物种；无新动力学输入。** 原文本作者从未发表 d=4 YM 收敛证明（系列止于 VI，VII 未发表），4D 仅 φ³₄（Williamson 1987，超可重整）。**G1 真定理保留：VI 的 Thm 1–2（格点–连续对偶）、纯规范退耦（精确）。**

---

## E. 一句话

**这轮造出的最好的三件（乘性搬运、`G″=S−D`、Fredholm 零点），都真正"碰上"了原子，但都把墙留在"连续统 + 搬到 2.5"那一侧——即 P3 的边际解析墙。** 三条新 no-go（周期、嵌套、刚性）把"为什么老/新工具那一类都不行"又钉死三层。**墙的坐标、名字、和逃逸口，现在全在册。**

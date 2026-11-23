# GL3 / GL4 / GL9 — 极限级粘合引理

> 2026-11-22 · glue-a · 4D 纯 Yang–Mills，`G` 紧单连通单群。
> **记号不重定义**：`μ_{L,a,β(a)}`、`Obs`、`S^{(n)}`、`θ`、`C_RP`、`Λ_{L,a}`、`β=2N_c/g²`、极限阶梯（阶段 1 `a→0` 固定 `L`，阶段 2 `L→∞`）**逐字引自** `research/gluing/REPORT.md` §1。本文件只增加引理层。
> 四件交付（C16）：每个引理给 ①精确陈述（含量词与常数） ②证明**或具名缺口** ③依赖 ④证伪测试。

---

## 0. 摘要（状态先行）

| 引理 | 交付结论 | 状态 |
|---|---|---|
| **GL3** | **✅-蕴含**（消费 `UVS(L) ∧ UVS-L`）：二者 ⇒ `{S^{(n)}_{L,a}}` 相对紧、子列弱-* 收敛、界关于 a 一致；`UVS-L` 下阶段 2 子列紧。**Lean 模块化**：紧性 + 可数对角已机器检验（`GL3Limits.lean`，`exists_subseq_tendsto_pi`）；**Minlos 作为显式假设**写在签名里（`gl3_imp_with_minlos` 的 `hMinlos`），`P` = 极限正性（由 GL4 供给）。**有界可观测类常数 = 1，是空洞的，不算成果。** | ✅-蕴含（Lean 模块 + 显式 Minlos） |
| **GL4** | **✅（无开放输入）**。**全锥机器检验**：`rp_transfer`（有界连续、Mathlib **真实弱收敛**）、`clustering_transfer`（`(U_a)` 是显式假设 `hU`）、`rp_l2_closure` + `rp_l2_closure_trunc`（**L² 截断闭包**：`A_R = A·1_{|A|≤R}`、`A_R → A` 于 `L²`、Hölder/Cauchy–Schwarz 收尾）。唯一显式假设是 `trunc A n ∈ L²`（由 `|trunc| ≤ |A|` 显然）。 | ✅（全锥 Lean） |
| **GL9** | **✅-蕴含**（消费 `GAP = CPB-W`）：一致隙 ⇒ 体积一致指数聚类。核心谱不等式已 **Lean 机器检验**（`GL9Decay.lean`，`spectral_decay`，公理 = `propext, Classical.choice, Quot.sound`）。隙本身是 `I-HANDOFF[m]` 的公开输入，属**假设**而非本引理内容。 | ✅-蕴含（核心 Lean） |

**分类（按本轮要求，不含糊）**：

```
GL4（作为引理）：无条件 ✅（全锥已验证）
GL4（作为组合步骤）：消费 GL3 的 UVS（弱极限存在）
GL3 = ✅-蕴含（点名消费 UVS(L) ∧ UVS-L）
GL9 = ✅-蕴含（点名消费 GAP = CPB-W）
```

板上不得把 GL4 读成"组合中也无条件"。

**三条关键发现**（本 session 的增量）：

1. **GL3 的"紧性"对有界可观测是免费的。** `{μ_{L,a}}` 推前到 `∏[-1,1]` 后自动相对紧（Tychonoff/Prokhorov），真正的数学内容是 (i) 无界可观测（smeared `tr F²`）的 a-一致界 = `UVS`，(ii) 极限的物理识别 = `I-RECON`/唯一性。**有界类是完全空洞的（常数 = 1），不得记作成果**；GL3 的全部内容在 `UVS ∧ UVS-L`。
2. **GL4 的 RP 部分白送；聚类部分不白送。** 弱极限**不**保持速率：massive Gaussian `mₙ→0` ⇒ 弱极限 = massless GFF，无指数聚类，而 RP 与弱收敛都成立。故 `MG` 必须写成**体积一致**——这是定理前提，不是文风。→ 新风险 **R28**。
3. **GL9 的常数与 `L` 无关是概率测度的免费推论**（有界可观测的 ℋ-范数 `≤‖·‖_∞`）。Lieb–Robinson 只在两处进入：无界 smeared 算子、以及"Hamiltonian 隙 ⇒ 转移矩阵隙"。且一般欧氏方向的速率是 **`δ/2`**（坐标方向技巧），不是 `δ`。→ 新风险 **R29**。

---

## 1. 共享假设（不自带变体，C1）

只列出本文件真正使用的字段。

```
(A0)  μ_{L,a,β(a)} 概率测度 = Z^{-1} e^{-S_W} ∏dU_e，S_W 为 Wilson 作用量，β=2N_c/g²   (§1.2)
(RT)  a ↦ β(a) 是 §1.4 的单一轨迹，scheme 常数 C 与其它八接口相同                         (§1.4, GL1)
(OB)  Obs = 规范不变局部算子 *-代数；S^{(n)} 为 Schwinger 函数；|tr U(C)|/N_c ≤ 1           (§1.3)
(RP)  有限 (L,a)：对 A∈C_RP，∫ \overline{θA}A dμ_{L,a} ≥ 0（Osterwalder–Seiler 1978）        (§2.2, GL2a)
(阶梯) 阶段 1：a→0 固定 L；阶段 2：L→∞。禁止联合极限                                       (§1.6)
```

`θ` 是时间反射 `(t,x)↦(−t,x)`（空间边固定、时间边取逆并共轭）；`C_RP` 是 `t>0` 半空间可测的 Obs 元素锥。本文件对**每个坐标方向** `μ=0,1,2,3` 使用同一套 `θ_μ, C_RP^μ`（超立方不变 + OS'78 对每个坐标方向的 link-plane 反射给出 RP）。

---

## 2. GL3 — 关于 a 的一致紧性/收敛

### 2.1 陈述（可记账形式：一条具名蕴含 + 一个显式空洞性说明）

**分类约定**：`✅` = 不消费开放输入；`✅-蕴含` = 消费**点名**输入并在陈述中写出。GL3 属后者。

**定义 `UVS(L)`（无界可观测的 a-一致界）。** 存在 `a₀(L)>0` 与常数 `C_{n,k}(L)<∞`，使对所有 `a ≤ a₀(L)`、`n ≤ n₀`、`f_i ∈ C_c^k(B(0,L/4))`：
```
|S^{(n)}_{L,a,β(a)}(f_1,…,f_n)| ≤ C_{n,k}(L) ∏_i ‖f_i‖_{C^k} ,                       (3.U1)
```
且 `Z_{L,a}(f) := Σ_n i^n/n! S^{(n)}_{L,a}(f^{⊗n})` 在 `0` 处**等度连续**（a-一致）。三个子句见 §2.2。

**定义 `UVS-L`（L-一致）。** 对每个 `n,k`，`sup_{L≥L₀} C_{n,k}(L) < ∞`。

**GL3-IMP（定理，✅-蕴含；消费 `UVS(L) ∧ UVS-L`）。** 在 `UVS(L) ∧ UVS-L` 下：
```
(i)   ∀ n,k：{S^{(n)}_{L,a} : a ∈ (0,a₀(L)]} 在 (C_c^k)′ 中相对紧，
      且界 (3.U1) 关于 a 一致；
(ii)  对每个 a_j→0 存在子列 a_{j_l} 与 S^{(n)}_L ∈ (C_c^k)′，使对**所有 n 同时**
      S^{(n)}_{L,a_{j_l}} → S^{(n)}_L （weak-*），且极限满足同一界；
(iii) 若 UVS 的等度连续子句成立，则 S^{(n)}_L 是某个概率测度 μ_L 的矩（Minlos）；
(iv)  在 UVS-L 下，{S^{(n)}_L} 关于 L 相对紧：阶段 2 的子列极限存在。
```
**这条蕴含的证明只是 `UVS/UVS-L` 的直接推论**（§2.2），没有新分析。**全收敛（非子列）仍需唯一性（W1，R19），不属 GL3。**

**GL3-B（有界可观测；常数 = 1）：空洞，不算成果。** 归一化 Wilson 环满足 `|Φ_j| ≤ 1`，故一切矩 `∈ [-1,1]`，"紧性"是 Tychonoff/Prokhorov 的**免费**推论，与 `a,L` 无关。它只说明"有界类没有障碍"，**不构成 GL3 的任何数学内容**；**不得记作 ✅**。（一维 Bolzano–Weierstrass 步骤已 Lean 化：`exists_convergent_subseq_of_abs_le_one`；多数列对角化是散文。）

**GL3-U / GL3-L 不再是独立陈述**：它们就是 `GL3-IMP` 的两个假设。

### 2.2 证明（GL3-IMP）

1. **Banach–Alaoglu。** `C_c^k(B(0,L/4))` 是赋范空间，其对偶的范数球弱-* 紧。（`S^{(n)}` 是 n-重线性泛函；经射影张量积 `⊗̂^n_π C^k` 视为线性泛函，Banach–Alaoglu 适用。）`(3.U1)` 说明每个 `S^{(n)}_{L,a}` 落在半径为 `C_{n,k}(L)` 的球内，故 `{S^{(n)}_{L,a}}` 相对紧。得 (i)。
2. **可数对角化。** 元组 `(n, f_1,…,f_n)`（取可数稠密集）可数；沿 `a_j→0` 对每个元组取子列，再对角线，得一个子列使**所有** `n` 的矩同时收敛。弱-* 极限满足同一界（范数的下半连续性）。得 (ii)。对一般 `f_i∈C_c^k` 用稠密族延拓。
3. **Minlos / 正性。** `Z_{L,a}` 在 `0` 等度连续 ⇒ 极限 `Z_L` 在 `0` 连续；每个 `Z_{L,a}` 是概率测度 `μ_{L,a}` 的生成泛函（正、归一），逐点极限仍正、归一；Minlos（等价地：柱代数上的正归一矩泛函 + Riesz–Markov）给出 `μ_L`。得 (iii)。正性这一步与 `GL4-RP` 同源（§3）。
4. **(iv)** 同一论证，常数换成 `sup_L C_{n,k}(L)`。

**`UVS` 的三个子句（具名缺口的内容，不是本 session 的定理）**：
- `(U1)` a-一致的 `L^p` 界：对 plaquette/场强的 smeared 多项式，`∫|P(F)|^p dμ_{L,a} ≤ C_{p,R}`，`a`-一致。Balaban 的小场/大场分解给出。
- `(U2)` 重正化减除的正确性：`tr F²` 有 `a^{-4}` 加性发散，必须按**同一** scheme 常数 `C`（§1.4、R26）减除；否则 `(3.U1)` 的 `C_{n,k}(L)` 不是 a-一致。这是 RT 进入 GL3 的唯一位置。
- `(U3)` `0` 处的等度连续性（即 `n!` 增长控制 `C_{n,k}(L) ≤ A^n n! ⋅ (L-多项式)`），Minlos 所需，强于 `(U1)`。
> 项目状态：`(U1)(U2)` 是 Balaban 框架对固定 `L` 盒子的已知结论，但**未被本项目重证**；`(U3)` 需要显式的 factorial 生长控制。故 GL3 整体是 **✅-蕴含**。

**空洞性自检（必须保留）**：若把 `UVS/UVS-L` 换成 `|Φ|≤1`，上述证明在 (1) 仍然"成功"（常数 = 1），这正说明**有界类的紧性是免费的**；GL3 的全部内容在 `UVS` 的 `(U1)-(U3)`。

### 2.5 Lean 化（已验证）：紧性 + 可数对角，Minlos 作为显式假设

模型：`k : ℕ` 是格距指标（`a_k → 0`），`j : ℕ` 枚举可数矩坐标 `(n, f₁,…,f_n)`；矩族为 `S : ℕ → (ℕ → ℝ)`，`|S k j| ≤ C j` 即 UVS 的 a-一致界。

```lean
theorem exists_subseq_tendsto_pi (S : ℕ → (ℕ → ℝ)) (C : ℕ → ℝ)
    (hS : ∀ k j, |S k j| ≤ C j) :
    ∃ (L : ℕ → ℝ) (φ : ℕ → ℕ), StrictMono φ ∧
      ∀ j, Tendsto (fun k => S (φ k) j) atTop (𝓝 (L j))

def HasRealization (L : ℕ → ℝ) : Prop :=
  ∃ μ : Measure (ℕ → ℝ), IsProbabilityMeasure μ ∧ ∀ j, ∫ x, x j ∂μ = L j

theorem gl3_imp_with_minlos (S : ℕ → (ℕ → ℝ)) (C : ℕ → ℝ)
    (hS : ∀ k j, |S k j| ≤ C j) (P : (ℕ → ℝ) → Prop)
    (hMinlos : ∀ L, (∀ j, |L j| ≤ C j) → P L → HasRealization L) :
    ∃ (L : ℕ → ℝ) (φ : ℕ → ℕ), StrictMono φ ∧
      (∀ j, Tendsto (fun k => S (φ k) j) atTop (𝓝 (L j))) ∧
      (∀ _hP : P L, HasRealization L)
```

- **步骤 1（紧性）依赖 Mathlib 的哪条**：具体矩编码用 `isCompact_univ_pi`（**Tychonoff**：紧区间之积紧）；抽象 Banach–Alaoglu 由 `WeakDual.isCompact_closedBall` 直接给出（本文件的 `banach_alaoglu_compact` 是它的包装）。
- **步骤 2（可数对角）**：`IsCompact.tendsto_subseq`（紧 + 第一可数 ⇒ 子列收敛）+ `tendsto_pi_nhds`（积收敛 = 坐标收敛）。全部机器检验。
- **步骤 3（Minlos）是显式外部假设**：`hMinlos` 就是 Minlos 定理（等价的 Prokhorov 核形式），**Mathlib 没有**，本文件**不证明**它，只把它逐字写进签名。`P` 抽象极限矩泛函的**正性/归一化**，在应用中由 `GL4-RP` 供给。
- 公理：`#guard_msgs` 冻结为 `[propext, Classical.choice, Quot.sound]`；`bash research/glue-a-limits/check-glue-a.sh` 复现。

### 2.3 依赖

| 方向 | 内容 |
|---|---|
| **消费** | `A0`/`GL1`（同一测度族与同一 RT）、`GL4-RP`（把极限矩泛函识别为 RP 测度的正性）、**`UVS(L) ∧ UVS-L`**（具名输入）、`MG`/`GL9`（可选：给 n=2 的 L-一致界） |
| **被消费** | `I-OS2` 的 A2（`S^{(n)}` 分布收敛）、`I-OS3` 的极限（B3 需先有极限）、`I-UNIQ` 的 TIGHT/子列、`GL4`（极限存在） |

> GL3-B 的有界类**不在**消费清单里（空洞）。Lean：`GL3Limits.lean` 的 `exists_subseq_tendsto_pi` 机器检验紧性 + 可数对角；`gl3_imp_with_minlos` 把 Minlos 作为显式假设；`GL4Transfer.lean` 另有 1D Bolzano–Weierstrass 备用。

### 2.4 证伪测试

- **F3-1（RT 必需，可执行）**：取**固定 β**（不走 RT）的序列。此时 `ξ_lat` 固定、`L` 固定 ⇒ `ξ_phys = a ξ_lat → 0`；未减除的 2 点函数 `~a^{-4}`。沿 RT 测 plaquette 与重正化 2 点函数：若 `a→0` 时无稳定 scaling window（`ξ/a` 不增长而 `ξ_phys` 塌缩），则 `(3.U1)` 的常数**不是** a-一致 ⇒ RT/UVS 假设不成立。这是标准 scaling test。
- **F3-2（因子 2，流程级）**：故意用 `48π²` 而非 `24π²`（R26/C15）。SU(3) `g²=1` 处的 `t=24π²/(11·3)=7.178` 是唯一真值；若数值匹配错的那个，`UVS` 的 scheme 常数 `C` 不一致，a-一致界在受检量上差 `e^{3.59}≈36`。检查：所有 `t_up/λ/β` 表项由 `C_RG:=24π²/11` 导出。
- **F3-3（量词反例，检查证明）**：若某"证明"只用 `|Φ|≤1` 而不碰 a-一致界，它就会对**任意**概率测度族"证明"存在连续极限。把它应用到"固定 β 的强耦合序列"：它声称存在连续极限，而该序列物理关联长度 `→0`。故正确证明必须在大场处消费 `UVS`；grep 证明中是否出现 `(U1)-(U3)`。
- **F3-4（紧性空洞性，自检）**：把 `J` 取成有限集，`(i)(ii)` 平凡。若证明没用"可数对角化 + 稠密族延拓"，它只证了有限维紧性，不足以给 `S^{(n)}` 全族的分布收敛。
- **F3-5（`UVS-L` 的证伪，可执行）**：把"L-一致"应用到 `T⁴_L` 上的 **massless GFF**（`λ₂=e^{-2π/L}→1`，无隙）：取 `f = 1_{B(0,L/4)}/|B(0,L/4)|`，则 `S^{(2)}_L(f,f) = ∫∫ G(x−y)f(x)f(y) ~ L²`，任何 L-一致界失败。这不证伪 `GL3-IMP`（`UVS-L` 是它的假设），而是证伪"`UVS-L` 无条件成立"；与 F9-2 一致。**任何声称 `UVS-L` 无条件的证明必须排除此例**（即必须使用隙或 RT）。

---

## 3. GL4 — 弱极限下 RP 与聚类的传递

### 3.1 陈述

**GL4-RP（RP 锥在弱极限下闭）。**

- 设定：`Ω` 紧度量；`θ:Ω→Ω` **连续**对合；`μ_a ⇒ μ` 概率测度；`(θ)_*μ_a = μ_a`；每个 `μ_a` 满足 RP：`∀A∈C_RP^b(Ω)`（有界连续、`t>0` 半空间可测），`∫\overline{θA} A dμ_a ≥ 0`。
  - **公共空间**：把每个格点测度经有界可观测族（归一化 Wilson 环的有限联合）推前到公共紧空间 `Ω ⊂ ∏[-1,1]`；`θ` 在这些坐标下是**连续的坐标置换**（反射把每条 Wilson 环映到另一条，环的集合按反射置换）。这是 `(a)` 的具体内容；若分块/反射面不对齐，`θ` 仍连续但 `(θ)_*μ_a=μ_a` 会失效（F4-4）。
- 结论：`(θ)_*μ = μ`，且 `∀A ∈ C_RP ∩ L²(μ)`，`∫\overline{θA} A dμ ≥ 0`。
- **常数**：无。**不需要**一致矩界。

**GL4-C（聚类界的传递）。**

- 设定：同 GL4-RP；`O,P` 有界规范不变局部算子，`‖O‖_∞,‖P‖_∞ ≤ M`；`c_a(x) := ⟨O P(x)⟩_{μ_a,c}`。
- 假设 `(U_a)`：`∀a∈(0,a₀]`、`∀|x|≥r₀`：`|c_a(x)| ≤ C e^{-δ|x|}`，`C,δ` **与 a 无关**。
- 假设 `(CONV)`：`c_a(x) → c(x)` 对每个 `x`（由 GL3 的 n=2 矩收敛给出）。
- 结论：`|c(x)| ≤ C e^{-δ|x|}`。常数 `C,δ` 精确继承。

### 3.2 证明

**GL4-RP 证明（完整）。**
1. `θ_*μ = μ`：`θ` 连续 ⇒ `∀f∈C_b(Ω)`，`∫f dθ_*μ = ∫ f∘θ dμ = lim_a ∫ f∘θ dμ_a = lim_a ∫ f dμ_a = ∫f dμ`。
2. **有界 `A`**：`f_A := \overline{θA}A ∈ C_b(Ω)`（`θ` 连续、`A` 连续、实/复乘连续）。弱收敛 ⇒ `∫f_A dμ_a → ∫f_A dμ`；每项 `≥0` ⇒ 极限 `≥0`。
3. **无界 `A ∈ C_RP ∩ L²(μ)`**：令 `A_R := A·1_{|A|≤R}`。支撑条件与代数闭包在截断下保持，故 `A_R ∈ C_RP`，且 `A_R → A` in `L²(μ)`。由第 2 步 `∫\overline{θA_R}A_R dμ ≥ 0`。Cauchy–Schwarz：
   ```
   |∫(\overline{θA}A − \overline{θA_R}A_R) dμ|
     ≤ ‖θA‖₂‖A−A_R‖₂ + ‖A_R‖₂‖θA−θA_R‖₂  → 0 ,
   ```
   故 `∫\overline{θA}A dμ = lim_R ∫\overline{θA_R}A_R dμ ≥ 0`。∎

**关于"RP 是对所有测试函数的不等式"（任务要求明确声明的连续性）**：
- 需要 `(a) θ` 连续（保证 `f_A ∈ C_b`）；`(b) C_RP` 在 `L²(μ)` 中闭（支撑条件 + Obs 的代数闭包）；`(c)` 极限测度存在（GL3）。
- **不需要** `sup_a ‖A‖_{L²(μ_a)} < ∞`：第 3 步只在**极限**测度里做截断。这正是 RP 与聚类的分界（聚类需要一致界，RP 不需要）。
- 只控制**有限个** `A` 的论证是错的：`C_RP` 不是有限生成的，必须走 L² 截断闭包。见 F4-3。

**GL4-C 证明（两行）。** `|c(x)| = lim_a |c_a(x)|`（绝对值的连续性与 `(CONV)`）`≤ C e^{-δ|x|}`（`(U_a)`）。∎
> 全部内容在 `(U_a)`。GL4-C 本身平凡；它是**正确陈述的位置**：聚类是"一致界 + 逐点收敛"，不是"弱收敛"的推论。

### 3.2.1 Lean 化（已验证）

`GL4Transfer.lean` 机器检验两条（公理审计见 §9）：

```lean
theorem rp_transfer
    (θ : Ω → Ω) (hθ : Continuous θ)
    {ι : Type*} {F : Filter ι} [NeBot F] (μs : ι → ProbabilityMeasure Ω)
    (μ : ProbabilityMeasure Ω) (hweak : Tendsto μs F (𝓝 μ))
    (A : Ω →ᵇ ℝ)
    (hRP : ∀ i, 0 ≤ ∫ ω, A (θ ω) * A ω ∂(μs i : Measure Ω)) :
    0 ≤ ∫ ω, A (θ ω) * A ω ∂(μ : Measure Ω)

theorem clustering_transfer
    {ι : Type*} {F : Filter ι} [NeBot F] (c : ι → ℝ → ℝ) (clim : ℝ → ℝ)
    (hconv : ∀ x, Tendsto (fun i => c i x) F (𝓝 (clim x)))
    (C δ : ℝ) (hU : ∀ i x, |c i x| ≤ C * Real.exp (-δ * x)) :
    ∀ x, |clim x| ≤ C * Real.exp (-δ * x)
```

- `hweak` 是 Mathlib 在 `ProbabilityMeasure Ω` 弱拓扑上的**真实弱收敛**（`ProbabilityMeasure.tendsto_iff_forall_integral_tendsto` 对一切 `Ω →ᵇ ℝ` 给出积分收敛），不是抽象的 `hconv`；`hθ : Continuous θ` 就是"声明的连续性"。
- `hU` 就是 `(U_a)`，**显式出现在签名里**。`clustering_transfer` 完整覆盖 GL4-C。

**L² 截断闭包（`GL4L2.lean`，已机器检验）**：

```lean
theorem rp_l2_closure_trunc
    (θ : Ω ≃ᵐ Ω) (hθ : MeasurePreserving (⇑θ) μ μ)
    (A : Ω → ℝ) (hA : MemLp A 2 μ)
    (hAR : ∀ n, MemLp (trunc A n) 2 μ)
    (hRP : ∀ n, 0 ≤ ∫ ω, trunc A n (θ ω) * trunc A n ω ∂μ) :
    0 ≤ ∫ ω, A (θ ω) * A ω ∂μ
```

- `trunc A n ω = A ω · 1_{|A ω| ≤ n}`；`trunc_tendsto` 给出 `A_R → A` 逐点，`trunc_abs_le` 给出 `|A_R| ≤ |A|`。
- `rp_l2_closure` 内部：控制收敛 ⇒ `∫ (A_R − A)² → 0`（**这就是 `A_R → A` 于 `L²`**）；`holder2`（Hölder，p=q=2）给出 `∫ ‖f‖‖g‖ ≤ ‖f‖₂‖g‖₂`；`nrm_comp`（`θ` 测度保持）给出 `‖f∘θ‖₂ = ‖f‖₂`；三角不等式收尾得 `|∫P_R − ∫P| ≤ 2‖A‖₂‖A_R−A‖₂ → 0`，再由 `hRP` 取极限。
- **范围（不含糊）**：GL4 现在**全锥机器检验**；唯一的显式假设是 `trunc A n ∈ L²`（写在签名里，由 `|trunc| ≤ |A|` 显然成立）。
- 复现：`bash research/glue-a-limits/check-glue-a.sh`；`GL4L2.lean` 用 `#guard_msgs` 冻结 `[propext, Classical.choice, Quot.sound]`，其余两文件由 `#print axioms` 输出并检查。

### 3.3 反例：为什么一致性不能去掉（GL4 的 sharpness）

ℝ⁴ 上的 **massive Gaussian 族**，质量 `m_n → 0`：
- 每个测度都 RP（Gaussian ⇒ RP），并满足 OS0–OS4 的相应版本；`c_n(x) = C_n e^{-m_n|x|}`。
- 协方差核 `G_{m_n}(x) → G_0(x) = 1/(4π²|x|²)`（`x≠0` 逐点；4D 中 `1/|x|²` 局部可积，故作为分布收敛），测度**弱收敛**于 massless GFF。
- 极限的 2 点函数 `~1/|x|²`，**无指数聚类**。

故 `(U_a)` 以任何固定 `δ` 失败，而 RP 与弱收敛均成立。**结论：`MG` 必须"体积一致"；弱极限不制造速率。**（R28）

### 3.4 依赖

| 方向 | 内容 |
|---|---|
| **消费** | `GL3`（极限存在 + 矩收敛）、`GL2a`/OS'78（有限 a 的 RP），GL4-C 另需**一个 a-一致的聚类界**（若格点隙在 a 上一致，由 `GL9` 提供；或直接由 `MG` 在阶段 2 提供） |
| **被消费** | `I-OS3` 的 B3（弱闭）、`I-UNIQ`（把 MG 搬到极限）、`I-OS2`（极限不变性的连续性） |

### 3.5 证伪测试

- **F4-1（一致界必须）**：Gaussian `m_n→0` 反例（可执行：`m=1/n` 的 4D massive Green 函数，验证 `G_m(x)→1/(4π²|x|²)` 且 `e^{-m x}→1`）。若去掉 `(U_a)` 后"定理"仍成立，则它必错。
- **F4-2（RP ⇏ 聚类）**：2D Ising 在 `β_c`：RP 成立、唯一 Gibbs 态，但幂律、无隙。故任何**只用** RP 的 GL4 证明是错的。审计：证明是否消费 `(U_a)` 或 `GL9`。
- **F4-3（有限测试族必错）**：fat-tail 序列 `μ_n = (1−1/n)δ_0 + (1/n)δ_{n²}`（`Ω=ℝ, θ=id`，平凡 RP）：弱极限 `δ_0`，`A(x)=x` 的 `L²(μ_n)` 矩 `n³→∞`。RP 本身不失败，但"用格点不等式直接对无界 `A` 取极限"在没有 L² 控制时**不成立**；正确证明必须出现"L² 截断 + 支撑闭包"。审计：grep 该步骤。
- **F4-4（θ-协变）**：若极限的 `θ` 不是连续对合（例如反射面与分块不对齐，违反 §1.7 约束 1），`θ_*μ=μ` 失败。检查与 os3 的 `(B1)` 判据（`A∈W∩RWR⁻¹`）一致。

---

## 4. GL9 — 一致隙 ⇒ 体积一致指数聚类

### 4.1 陈述

- 设定：`L,a`，`N=L/a∈2ℕ`，`μ_{L,a}` 为 Wilson 测度。固定坐标方向 `μ`，令 `T=T_{L,a}` 为该方向的转移矩阵，作用在 `ℋ = L²(Ω_slice)`（物理上取其规范不变中性扇区 `ℋ_gi ⊂ ℋ`），满足
  ```
  (T1)  T 自伴、正、‖T‖=1，TΩ=Ω，‖Ω‖=1                     （GL2a/OS3）
  (GAP) ∃δ>0（与 L,a 无关）：∀x∈ℋ_gi，x⊥Ω ⇒ ‖T x‖ ≤ e^{-aδ}‖x‖
  ```
  （等价：`λ₂(T|_{ℋ_gi}) ≤ e^{-aδ}`；物理隙 `Δ ≥ δ`。这是 `I-HANDOFF[m]` 的输出，公开输入 `CPB-W`。）
- 结论：存在 `C=C(δ,ρ,M)<∞`（**与 `L,a` 无关**），使对所有规范不变局部算子 `O,P`（`‖O‖_∞,‖P‖_∞≤M`，可由时间厚度 `≤ρ` 的切片表示）以及所有 `0<|x|≤L/2`：
  ```
  |⟨O(0)P(x)⟩_{L,a,c}| ≤ C e^{-(δ/2)|x|} .                                             (4.1)
  ```
  更精确：时间方向 `e^{-δ'|x⁰|}`（`0<δ'<δ`）；任一坐标方向 `e^{-δ'|x^μ|}`；合并得 `δ'/2`。

### 4.2 证明

**引理 A（Ω⊥ 在 T 下不变 + 范数衰减）。** 若 `x⊥Ω`，则 `T x⊥Ω`，且 `‖Tⁿ x‖ ≤ e^{-aδn}‖x‖`。
> `T x⊥Ω ⟺ ⟪Ω,Tx⟫=0`。由 `T` 自伴：`⟪Ω,Tx⟫=⟪TΩ,x⟫=⟪Ω,x⟫=0`。范数由 `(GAP)` 对 `Tⁿx` 归纳。

**引理 B（沿真空方向分解）。** `Tⁿ B = ⟪Ω,B⟫•Ω + Tⁿ(B − ⟪Ω,B⟫•Ω)`，且 `‖B − ⟪Ω,B⟫Ω‖ ≤ ‖B‖`。
> 对 `n` 归纳：`n=0` 是 `P` 的定义；归纳步用 `TΩ=Ω` 与线性。

**主证明。**
1. 时向：令 `t=na`。由引理 B 与 `⟪A,Ω⟫=⟪Ω,A⟫`，
   ```
   ⟨O|Tⁿ|P⟩_c = ⟪O, Tⁿ P⟫ − ⟪Ω,P⟫⟪O,Ω⟫ = ⟪O, Tⁿ(P − ⟪Ω,P⟫Ω)⟫ ,
   ```
   故 `|⟨O|Tⁿ|P⟩_c| ≤ ‖O‖ ‖Tⁿ(P−…)‖ ≤ e^{-δt}‖O‖‖P − ⟪Ω,P⟫Ω‖ ≤ e^{-δt}‖O‖‖P‖`。
2. 角向：每个坐标方向都有转移矩阵（超立方不变 + OS'78 逐方向 link-plane RP）；平移不变（周期测度）⇒ 把 `O` 移到 `x^μ=0`。取 `|x^μ|` 最大的坐标方向，`|x| ≤ 2 max_μ|x^μ|`，于是该方向上的 `e^{-δ'|x^μ|}` 给出 `e^{-(δ'/2)|x|}`。
3. 常数与 `L` 无关：`‖O‖ = ‖O‖_{L²(μ)} ≤ ‖O‖_∞`，因为 `μ` 是**概率**测度且 `|O|≤M`。时间厚度 `ρ` 的切片表示引入一个只依赖 `δρ` 的因子，吸收进 `C`。
4. 周期性：`t ≤ L/2` 时镜像项 `≤ e^{-δt}`（距离 `L−t ≥ t`）；`|x|>L/2` 不在陈述内（缠绕）。
5. 扇区：`(GAP)` 必须取中性/规范不变扇区；Wilson/Polyakov 环的隙是另一本账（`I-HANDOFF[σ]`）。

**Lean 化（已验证）**：`GL9Decay.lean` 机器检验引理 A+B 的算子形式
```
|⟪A, Tⁿ B⟫ − ⟪Ω,B⟫⟪A,Ω⟫| ≤ λⁿ ‖A‖‖B‖ ,      -- GlueA.spectral_decay
```
条件：`T` 线性（`E →ₗ[ℝ] E`）、`TΩ=Ω`、`T` 保持 `Ω⊥`（`PreservesPerp`）、且在 `Ω⊥` 上 `‖Tx‖≤λ‖x‖`（`GapOnPerp`，`0≤λ`）。**不需要有限维**——论证是纯代数的，故对无限体积转移矩阵同样有效。这是 GL9 唯一的非平凡代数步骤（`T` 自伴 ⇒ `T` 保持 `Ω⊥`；`offVac Ω x = x − ⟪Ω,x⟫Ω` 是去掉真空方向）。其余（`λ=e^{-aδ}`、`t=na`、常数、方向、扇区）是记账。
- 复现：`bash research/glue-a-limits/check-glue-a.sh`；末行 `ALL GLUE-A CHECKS PASS`。
- 公理审计：`#print axioms GlueA.spectral_decay` → `[propext, Classical.choice, Quot.sound]`，无 `sorryAx`。
- 证明要点（对应 `REPORT.md` 引理 A/B）：`inner_offVac`（正交性）、`norm_offVac_le`（投影不增范数）、`iterate_decomp`（沿真空分解）、`perp_iterate`（`Ω⊥` 不变）、`iterate_norm_le`（几何衰减）。

### 4.3 Lieb–Robinson 在哪里进入（以及哪里不需要）

- **不需要**：有界局部算符的 `ℋ`-范数 `≤‖·‖_∞`——因为 `μ` 是概率测度。故 §4.2 第 3 步**不含** LR。
- **需要**（具名 **`LR-G`**，标准）：(i) 无界 smeared 场强算子的局部 `ℋ`-范数界（`‖|O_f⟩‖ ≤ C_f`，与 `L` 无关）；(ii) 若隙是 Kogut–Sitaram Hamiltonian 的隙，用 Duhamel + LR 把它搬到欧氏转移矩阵；(iii) 若要在**无限体积**直接陈述，用 LR 把有限体积条件期望与无限体积态比较。
- `LR-G` 的常数只依赖局部作用与维数，不依赖体积（Nachtergaele–Sims 型；对格点规范模型需在规范不变代数上证）。**这与 `UVS`/`CPB-W` 不同：`LR-G` 是标准可证引理，不是开放问题。**

### 4.4 依赖

| 方向 | 内容 |
|---|---|
| **消费** | `GL2a`/`OS3`（`T≥0`）、OS'78（逐方向 RP）、`I-HANDOFF[m]`（`GAP`，公开输入 `CPB-W`）、可选 `LR-G`（无界类） |
| **被消费** | `I-UNIQ` 的 `W1_0`（经 `MG`）、`I-OS2`、`I-SPEC`（谱表示）、`GL3-L`/`GL4-C`（一致聚类界） |

### 4.5 证伪测试

- **F9-1（一致隙不可弱化）**：2D Ising 在 `T_c` 以下、周期 bc：`Δ_L ~ e^{-σL} → 0`；假设失败，结论也失败（`±` 两态）。故若某证明只用"`Δ_L>0 ∀L`"，它在 Ising 上会"证明"体积一致聚类——错。检查：证明是否在某处用到了 `δ` 与 `L` 无关。
- **F9-2（自由无质量场）**：`T⁴_L` 上的 massless GFF：`λ₂=e^{-2π/L}→1`；假设失败；幂律。同上。
- **F9-3（数值，可执行）**：SU(3) `β≈2.6`，`L=16,20,24` 的 `0^{++}` 关联。计算 `m_eff(t) = −log(c(t+1)/c(t))` 与拟合前置因子 `C(L)`。GL9 预言：平台与 `L` 无关；`C(L)` 常数。**证伪**：平台随 `L` 下降（隙闭合）或 `C(L)` 随 `L` 增长（体积相关范数/LR 常数误用）。同时检查对角方向速率 `≥δ/2`（不要声称 `δ`）。
- **F9-4（扇区量词）**：用 Wilson 环（带电）检验中性隙不控制它：把 GL9 用于弦扇区是**量词反例**（陈述用了错误的隙）；弦的速率应由 `I-HANDOFF[σ]` 给出。
- **F9-5（LR 常数）**：若 `C(L)` 随 `L` 增长，检查是否把 `LR-G` 的常数误用成了体积相关。

---

## 5. 组合：这三条引理在依赖图中的位置

```
GL3-B（有界紧性，✅）──┐
GL4-RP（RP 闭，✅）  ──┼──► 阶段 1：μ_L 存在 + RP（I-OS3_fin → I-OS3_lim）
GL3-U（UVS，条件）   ──┘
        │
I-HANDOFF[m]（GAP = CPB-W，公开输入）
        │
        ▼
GL9（✅ 蕴含；核心 Lean 化） ──► 格点/有限体积：体积一致指数聚类
        │                                   │
        │                     GL4-C（✅，需一致界）
        │                                   ▼
        └──────────────────────────────► MG（体积一致）
                                            │
                                 GL3-L（有界类免费）
                                            ▼
                                 阶段 2：{μ_L} 子列紧 + 子列极限
                                            │
                                 I-UNIQ（W1_0）+ 唯一性
                                            ▼
                                          EXIST
```

**两个不能省的东西**：
- `GL3-U` 的 `UVS`：否则无界可观测（smeared `tr F²`、OS4 的 Schwartz 测试函数）的连续极限没有定义；
- `I-HANDOFF[m]` 的 `GAP`（`CPB-W`）：否则 GL9 没有输入。

**GL4 与 GL9 的非输入部分都是真的**：RP 的弱闭性是免费的；一致隙 ⇒ 聚类的谱论证是初等的（且已 Lean 化）。**本 session 没有证明质量隙**——那是 `CPB-W`。

**一个顺序警告**：GL9 的 `GAP` 必须关于 `L` 一致。若只在每个 `L` 有某个正隙而 `inf_L Δ_L = 0`，则 §4.2 的界变成 `e^{-δ_L t}`，常数 `C` 会随 `L` 爆炸，与 F9-1 的 Ising 情形相同。这就是 `MG` 必须"体积一致"的另一处体现。

---

## 6. 新增风险（接 R0–R27）

| # | 风险 | 位置 | 级别 | 后果 | 缓解 |
|---|---|---|---|---|---|
| **R28** | **弱极限不保持聚类速率**：massive Gaussian `mₙ→0` 的弱极限是 massless GFF（RP、OS 公理都满足），无指数聚类 | GL4 ↔ MG/I-UNIQ | 高（陈述层） | 任何"RP + 弱收敛 ⇒ 聚类"的论证错；`MG` 若写成"逐 L 有隙"而非"体积一致"，则不可用 | GL4-C 显式消费一致界；GL4 的证据包含 Gaussian 反例；`MG` 的陈述强制体积一致（§0.2 判断 4） |
| **R29** | **GL9 的速率是 `δ/2` 而非 `δ`；常数与 L 无关只对有界类免费** | GL9 ↔ I-UNIQ/I-CONF | 中 | 下游若用 `e^{-δ|x|}` 作面积律输入会高估 2 倍；若把 LR 当承重部分，会掩盖真正输入（GAP） | 陈述写 `δ'/2`（`δ'<δ`）；明确 LR 只在无界类进入；F9-3 数值检验 |
| **R30** | **GL3 的"紧性"对有界可观测是空洞的** | GL3 ↔ §3.4 追踪板 | 中（记账） | 把 GL3 记 ✅ 时若只说"μ_{L,a} 紧"，会掩盖 `UVS` 缺口 | GL3 的 ✅ 只能记在 GL3-B（软）；GL3-U 必须写 `UVS` |

---

## 7. 诚实的边界

1. **GL4 的机器证据 = 三条 Lean 定理 + 一个反例脚本。** `rp_transfer`（有界连续、真实弱收敛）、`clustering_transfer`（`(U_a)` 显式）、`rp_l2_closure_trunc`（L² 截断闭包）已编译，公理标准；反例脚本验证 `G_m(r)→1/(4π²r²)`（`m=1/64` 时比值 `0.99–0.9994`），即"去掉一致性 ⇒ 聚类丢失"。**GL4 不再有散文步骤**；唯一显式假设是 `trunc A n ∈ L²`（明示在签名里，由 `|trunc| ≤ |A|` 成立）。
2. **GL9 作为蕴含已证**，核心谱不等式已 Lean 化（`GL9Decay.lean`，`spectral_decay`，公理仅标准三条，`check-glue-a.sh` 复现）。但 `GAP` 是 `CPB-W`——**公开输入**，不是 GL9 的内容。GL9 不能证明质量隙。
3. **GL3 是 ✅-蕴含**：`UVS(L) ∧ UVS-L ⟹ 紧性/子列/一致界`。有界类常数 = 1，**空洞**，不算成果；内容全在 `UVS` 的 `(U1)-(U3)`。Lean 覆盖紧性与可数对角（`exists_subseq_tendsto_pi`）以及带显式 Minlos 假设的模块（`gl3_imp_with_minlos`）。**Minlos 本身是显式外部假设**（`hMinlos`），不是本文件证明；Mathlib 无此定理。
4. **我没有验证** `UVS` 的具体常数，也没有验证 `LR-G` 对 Wilson 格点规范模型的逐字适用性（引用的是标准框架，未逐页核对）。
5. **GL9 的 `δ/2` 不是最优**；若需要更好的角向速率，应显式构造对角方向的转移矩阵（若能），而不是靠坐标方向技巧。
6. **GL9 的 Lean 化是"核心谱不等式"，不是"GL9 整条"**：从 `λ=e^{-aδ}`、`t=na` 到 `e^{-δt}` 的代入、四方向合并、以及"范数 ≤ ‖·‖_∞"都在散文里。它们初等，但**没有逐行机器检验**。

---

## 8. 与 §3.0 供给矩阵的一致性

- `GL3 → I-OS2/I-OS3/I-UNIQ`：一致（§3.0 行"GL3 紧性"）。
- `GL4 → I-OS3(B3)/I-UNIQ`：一致（§3.0 行"GL4 弱极限"）。
- `GL9 → I-UNIQ（经 MG）`：一致（§3.0 行"GL9"）。
- 新增的两条**内部**依赖（不在 §3.0 的表里，但由本报告建立）：
  - `GL4-C` 消费 `GL9`（a-一致聚类界）；
  - `GL3-L` 消费 `MG`（L-一致界，仅无界类）。
  这两条应补进 §3.0 的下一版。

---

## 9. 复现

```bash
bash research/glue-a-limits/check-glue-a.sh
```

依次执行：
1. 用与 `lean/check-handoff.sh` 相同的 pin（`v4.33.0-rc1` + vendored Mathlib olean cache）编译三个 Lean 文件，并强制公理足迹为 `[propext, Classical.choice, Quot.sound]`、无 `sorryAx`：
   - `GL9Decay.lean`：`#print axioms GlueA.spectral_decay`（GL9 核心）；
   - `GL4Transfer.lean`：`#print axioms` 覆盖 `GlueA.rp_transfer`、`GlueA.clustering_transfer`、`GlueA.exists_convergent_subseq_of_abs_le_one`（GL4 有界层 + GL3 一维紧性）；
   - `GL4L2.lean`：`#guard_msgs` 冻结 `GlueA.rp_l2_closure_trunc`、`rp_l2_closure`、`holder2` 的公理足迹（GL4 全锥）；
   - `GL3Limits.lean`：`#guard_msgs` 冻结 `GlueA.exists_subseq_tendsto_pi`、`gl3_imp_with_minlos`（GL3 紧性 + 对角 + 显式 Minlos）。
2. 运行 `gl4_gaussian_counterexample.py`（`python -I`，仅标准库）：4D massive Green 函数 `G_m(r)=mK₁(mr)/(4π²r)` 经精确积分表示数值求值，验证 `m→0` 时逐点收敛到 `1/(4π²r²)`。

预期末行：`ALL GLUE-A CHECKS PASS`。

文件清单：

| 文件 | 内容 |
|---|---|
| `REPORT.md` | 本文件（三引理的四件交付 + 组合 + 风险） |
| `STATUS.md` | 一页状态：证明/条件/缺口/下一步 |
| `GL9Decay.lean` | GL9 核心谱不等式的机器检验（`GlueA.spectral_decay`） |
| `GL4Transfer.lean` | GL4-RP（真实弱收敛）、GL4-C（`(U_a)` 显式）、GL3 一维紧性步骤的机器检验 |
| `GL4L2.lean` | GL4 全锥：L² 截断闭包（`trunc`、`trunc_tendsto`、`holder2`、`nrm_comp`、`rp_l2_closure_trunc`） |
| `GL3Limits.lean` | GL3 紧性 + 可数对角（`exists_subseq_tendsto_pi`）与带显式 Minlos 假设的模块（`gl3_imp_with_minlos`） |
| `gl4_gaussian_counterexample.py` | GL4-C sharpness 的可复现数值反例 |
| `check-glue-a.sh` | 上述三项的一键复现 |

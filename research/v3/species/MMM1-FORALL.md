# MMM1-FORALL — GAP 的 `∀L` 轴：有限体积隙 vs 无穷体积隙（方向 / 定理 / 是否免费）

> 第 MMM1 轮 · 子代理 `mmm1-fl` · `D:\Dev\math`。**写域**：只新建本文件（≤85 行）；既有文件一字未改；未 `commit`；零新 MC、零 Lean；未动 `docs/position-current.html`（交 `lead`，`AGENTS §6`）。标签只 `[T]`/`[NG]`；无百分数；每条带「下一手」。
> **靶（逐字转引）**：`GAP: ∃δ>0 ∀L,a`，真空扇区（中性、局域）最轻态质量 `≥δ`（`AA4-CARRIER:45`、`KKK1-MLAM:4`）。kkk1 判（转引 `KKK1-MLAM:13,26`）：`sup_a[ρ𝒢]<∞` 是 `L3` 的 `a`-内容（**固定 `L`**），而 `χ₂>0 ⇒ C₂(L)≳cχ₂L⁴ ⇒ sup_L=∞`。
> **口径**：`m(L) := −ln(λ₂/λ₁) = inf spec(H_L\{E₀})`（`D03-ML:43`）；`d=4`、`SU(N)`、Wilson、**周期**、`T³_L`。`L` 窗与 `a` 窗不互借（`AGENTS §2.1`）。

## §0 判词（三行，各带下一手）

| # | 断言 | 判 | 一句（具名障碍） | 下一手 |
|---|---|---|---|---|
| **0.1** | 方向：`m(L) ≥ m(∞)`？ | **`[NG]`**（类内） | **"体积增强"是边界条件的产物，不是体积的**：周期自由场 `m(L)≡m` **恒等**（`D03-ML:142`）；能抬隙的是**扭曲** bc（`√(m²+π²/L²)` 递减）。相互作用的 `L`-修正指数小、符号不定：仓内唯一 `0⁺⁺` 公式给**负**（`Δm∝−g²e^{−√3mL}/(mL)`，`CONSTRUCT-LUSCHER:102`），TFIM `g>1` 数值给**正**（`D03-ML:130`） | 取同 β 同 `a` 的 `m(L)` 序列（`L=6,8,10,12`），只读 `k=0∩0⁺⁺`，判**符号**（§3 末） |
| **0.2** | 定理：有严格结果吗 | **`[NG]`**（torus）/ **`[T]`**（仅 `S¹` 几何） | **min–max 是唯一合法形状，且只在 trial 态类嵌套时成立**：`ℝ³×S¹_L` 上 `m(S¹_L) ≥ m(ℝ⁴)` **`[T]`**（`FORGE5-SMALL-S1:14` 逐字，且明写它只证 `m(ℝ⁴) ≤ m(L)`）；`T³_L` 上两个配置空间不嵌套（周期识别、无边界可钉、电能项不保周期子流形）⇒ 同一论证缺席 | 若要用 min–max，须先造 `T³_L ↪ T³_{L'}` 的**能量不增扩张**（一条真引理，仓内 0 条） |
| **0.3** | 若 `m(L) ≥ m(∞)`，`∀L` 免费？ | **`[NG]`**（不免费：量词转移，剩余合项 = 靶） | `GAP ∧ 单调 ⟺ 单调 ∧ [m_∞≥δ]`；而 `m_∞` 与 `Δ_cont:=inf spec(H_∞\|_{Ω^⊥})` 是**同一数**（§3.2）⇒ 剩余合项就是 Clay 的 `GAP(Δ>0)`；反方向 `GAP ⟹ Δ_cont>0` 便宜，`⟸` 恰是单调性 | 记 `∀L` 为**未兑**；若符号为负，`∀L` **严格强于** `Δ_cont>0` ⇒ GAP 须写成 `inf_L m(L)>0`（= 一致界） |

**一句判**：`m(L) ≥ m(∞)` 在**周期**几何里没有免费来源（自由场恒等），符号未定；定理只在 `S¹` 几何成立且给的是目标的**上**界；`∀L` 不是折扣，是**把两轴并成一轴、留下的那一轴等于靶**。★ 是对准了错锁的钥匙。

## §1 方向

**1.1 三个机制：一真两伪（`[T]`，模型内）**

| 机制 | 给什么 | 在周期 `T³_L` 上 |
|---|---|---|
| 动量量子化 / 域单调 | Dirichlet Laplacian：`Ω₁⊂Ω₂ ⇒ λ_k(Ω₁)≥λ_k(Ω₂)`；**Neumann 不成立**（逐字 "domain monotonicity does not hold for the Neumann Laplacian"）【读摘要·片段；`web_fetch` 取 PDF 失败，见 §4#5】 | **不适用**：周期 bc 无边界，常数模 `k=0` 不动 ⇒ 自由场 `m(L)≡m` **精确**（`D03-ML:142`） |
| 边界 / 拓扑 | 扭曲 bc：`√(m²+π²/L²)` **严格递减** | 只在扭曲 bc；`D03-ML:143` 逐字"自由场里 `m(L)` 的全部非平凡 `L`-依赖来自拓扑/边界条件，**不是"体积"**" |
| 自能云（盒内求和 > 积分） | `Δm ∝ −g_{0⁺⁺0⁺⁺0⁺⁺}²·e^{−√3mL}/(mL)`，比例带**负号**（`CONSTRUCT-LUSCHER:102`） | 纯 YM `0⁺⁺` 里唯一写下的公式；按 `D03-ML:52` 的 `Δm:=M−m` 读 ⇒ 有限 **<** 无穷。**符号约定未逐字核 ⇒ 只作证据，不作定理** |

**1.2 精确模型（转引 `D03-ML:147–155`，`[T]`/`[T-数值]`）**：周期自由标量 = **常数**；扭曲自由标量 = 递减到 `m∞`；TFIM `g>1` = `m_∞+正修正` 从上方递减；M6 = `min(ω,½(2π/L)²)` 递减；2D YM 的"递增"是 torelon `σL→∞`（`m/L=0.4038776782` 逐位恒定），**不是隙**。⇒ **无一例给"递增 + 有限正极限"**。
**1.3 反向见证（`[T-数值]`，`D03-ML:23–25`）**：TFIM `g<1`（有序相）`m(L): 3.55e−2→7.25e−5` 而真隙 `=1` ⇒ 真空近简并时 `m(L)` **测的是隧穿劈裂，不是隙**。⇒ 方向问题不只是"未证"，是**对象可能不是隙**。

## §2 定理

**2.1 唯一合法形状 = min–max（变分原理）**：`m = inf spec(H)` 是**下确界** ⇒ "加 trial 态只会压低、不会抬高"。小体积 trial 类若能源不增地嵌入大体积，则 `m(小) ≥ m(大)`。`FORGE5-SMALL-S1:14` 在 `S¹` 几何给出这条 **`[T]`**，并明写其代价：正的小-`L` 界只封 `m(ℝ⁴) ≤ m(L)`，是目标的**上界**。
**2.2 `T³_L` 上嵌套缺席（`[T]`，结构事实）**：`H_L`、`H_{L'}` 作用于**不同**配置空间；`L`-周期子流形在电能项（`∂²`）下不闭合 ⇒ 无能量不增扩张，min–max 无从施加。
**2.3 PF 给什么、不给什么（`[T]`）**：`T_a>0`、`λ₁=1` 单重、`λ₂<λ₁`（`B01-RP:11` 逐字"**没有 `λ₂` 的任何界**"）⇒ **逐 `L` 成立，不给 `L` 间比较**。`SU(N)` Wilson 非铁磁（全表示展开号不定）⇒ **无 FKG ⇒ 无体积单调机器**（`G02-LINF:102`、`GH3-MONOTONE:104`）。
**2.4 文献（`[T-文献·摘要级]` / `[NG]`）**：Lüscher CMP **104**(1986)177–206 只给 `Δm→0` 的**指数率**、**不给符号**（`D03-ML:52`）；"`m(L)` 对 `L` 单调"的定理检索 **0 条**（`D03-ML:55`、`BN-RESIDUAL:148`，工具 OK）。

## §3 `∀L` 是否免费

**3.1 形式化（`[T]` 逻辑）**：记 `U := [∀L: m(L)≥m(∞)]`。则 `GAP ∧ U ⟺ U ∧ [∃δ∀a: m_∞(a)≥δ]`。⇒ `U` **只做量词消去，不产生数**。
**3.2 `L=∞` 那一端 = Clay 的 `Δ_cont>0`（`[T]` 同一对象 + `[NG]` 存在性）**：`m_∞ := lim_L m(L) = inf spec(H_∞|_{Ω^⊥}) = Δ_cont`（定义同一，`D03-ML:43`）；而 `GAP(Δ>0)` 逐字是 `Clay = AX ∧ NT ∧ GAP(Δ>0)` 的合项（`B2-bypass-gap:6`）。⇒ 用 `L=∞` 端替换 `∀L`，恰是**补上那条未证的 `U`**：零折扣。另注：`m_∞≥δ` 还需唯一真空（`GAP⇏NT` 是另一轴，`AA6-NT:52`），不计在此。
**3.3 若符号为负，`∀L` 严格强于靶（`[NG]`）**：`m(L) ≤ m(∞)` ⇒ `inf` 在有限 `L` ⇒ `∀L` 变成 `inf_L` 一致界 = 本仓的墙（`RATE-UNIFORM`、`CTL-CLOSEST`；`D03-ML:31` 逐字"那不是新入口，那是墙"）。⇒ **两支都不免费**：上支靠未证的 `U`，下支靠未证的一致界。
**3.4 ★ 的唯一真贡献 = 拆掉一条伪障碍（`[T]` 对象判别）**：`χ₂L⁴` 讲的是**矩**的 `sup_L C_{2,k}(L)=∞`（`HS-L3-5:48–54`），GAP 要的是**隙**的 `inf_L m(L)`；前者是 `C₂` 的 `L¹` 质量、后者是谱底（`D03-ML:43`）⇒ 不同对象。故 kkk1 的"`L`-轴挡 GAP"是**类型错**，应从 GAP 的障碍栏删除，只留 `a` 轴（`L3`）。

**下一手（零新数学，可执行）**：用已有 `GEVP-4D`/`WG` 数据取 `m(L)` @ 同 β 同 `a`、`L=6,8,10,12`；先过 torelon 诊断（`m/L≈const` ⇒ 绕行；`m→const` ⇒ 体态；`m∝1/L²` ⇒ 无质量塔，`D03-ML:171–173`），再只读 `k=0∩0⁺⁺∩N-ality 0`（`HZZ-M2:51`：不钉扇区则 `m(L)~½(2π/L)²→0`，读到的是塔）。判 `m(L)−m(∞)` 的**符号**：符号定 ⇒ `U` 真值定 ⇒ `∀L` 免费与否才可判。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 检索式 / 调用 | 命中 |
|---|---|---|---|---|
| 1 | `grep` | **OK** | ①`kkk1\|KKK1` @`research/v3` ②`χ₂\|chi_2\|cχ₂L⁴` @`research/v3` ③`m\(L\)\|m\(∞\)\|m_∞\|体积.*单调` @`research` | ①**6** ②**66** ③**250/276**（spill） |
| 2 | `glob` | **OK** | `research/v3/species/MMM*` · `KKK*` | MMM **1**（`MMM2-RGINV`；`MMM1` 本文件首建）· KKK **7** |
| 3 | `read` | **OK** | `KKK1-MLAM`(64 全)·`D03-ML:1–70`·`FORGE5-SMALL-S1`(22 全)·`CONSTRUCT-LUSCHER:80–124` | 4 份 |
| 4 | `web_search` | **OK** | ①`finite volume mass gap monotonicity in volume lattice gauge theory` ②`Luscher volume dependence of the energy spectrum stable particle states mass shift sign` ③`correlation length increases with volume FKG monotonicity ferromagnetic theorem` ④`domain monotonicity eigenvalues Dirichlet Laplacian inclusion theorem Rayleigh Faber Krahn` ⑤`finite volume gap exponentially small tunneling degenerate vacua symmetric boundary condition` | 各 10 源；可用 4 条（Dirichlet 域单调 + **Neumann 反例**；Lüscher I/II；隧穿塌缩）。"`m(L)` 单调定理"仍 **0** |
| 5 | `web_fetch` | **`TOOL-FAIL`** | `math.ucdavis.edu/~saito/tutorials/ssp14tutorial.pdf` | `unsupported content type "application/pdf"` ⇒ 域单调/Neumann 反例只到**检索片段级**，已按【读摘要·片段】标，不升 `[T]` |
| 6 | MC / Lean / `git` | **0 调用** | — | — |

**`TOOL-FAIL` 1（#5，已降级标注）；空响应 0；`PATH-MISS` 0；MC/Lean/`commit` 0。** 直读：`KKK1-MLAM`、`D03-ML:1–70`、`FORGE5-SMALL-S1`、`CONSTRUCT-LUSCHER:80–124`；其余为命中行**转引**（`B01-RP:11`、`BN-RESIDUAL:148`、`GH3-MONOTONE:104`、`G02-LINF:102`、`HZZ-M2:51`、`HS-L3-5:48–54`、`B2-bypass-gap:6`、`AA6-NT:52`、`AA4-CARRIER:45`）。

*MMM1-FORALL · 只新建本文件；出口：方向 `[NG]`（周期几何无体积增强，自由场恒等）· 定理 `[T]` 仅 `S¹` 几何 / `[NG]` torus · `∀L` 不免费（剩余合项 = `Δ_cont>0`）；下一手 = §3 末的符号判决。*

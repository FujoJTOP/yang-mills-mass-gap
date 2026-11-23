# A4-SPECTRAL-BOOTSTRAP — SHAPE A round 4：平面谱密度的 bootstrap 不闭合 #5 —— moment problem 的**可行性非空且隙盲（自由场是实现），严格弱于 H1**；而"逼出隙"必须加谱边数据 `Λ` = H1 换名

> 2026-11-23 · 判词 (verdict). **不存在**任何有限闭集，由这五条精确约束（SD sum rule `E[H]=E[f²]`、Ward `G″·Var(a)=1`、范数正性、crossing/RP、AF UV 尾）在**不假设**的前提下逼出 `m₁²≥λ`。精确的失败是**方向**：每条约束都是**矩（∞ 处 Taylor germ）或 Schur-凹正性（上界机）**；隙是**支撑下确界（edge）**，是有限矩永不见、且可行集里**被 0 实现**的量。**关键裁决**：moment problem 的**可行性不归约到 H1——它严格更弱**：可行集非空，自由场（两光子 cut，`f₁²=0`）就是解，故"可行 = 不排除"对无隙点也成立，零隙信息。**恰恰"逼出隙"（让无隙点不可行）才归约到 H1**——须加被排除的谱边数据（第二能级隙 `Λ`）= H1 换名。**1 名义通过（非微扰）；#2/#3/#5 败；#4 空洞（"败于 f₁²=0"正是它排不掉自由场之同一事实）。**

## (1) 五条约束全是"上界机"（方向错）
- **SD sum rule** `E[H]=E[f²]`：per-direction 的**总权重/一阶矩**，钉总权、不钉低-q 分布（ABC2-DUALITY §1）。无条件等分 `E[H]−Var≈0`，条件化即隙（ABC2-SYNTHESIS §4）。
- **Ward** `G″·Var(a)=1`：导出二次式 `x²−E[H]x+E[H]≥0`，判别式 `Δ=E[H](E[H]−4)`。实测 `E[H]≈2.0`（SU3 β=5.7: 2.03/1.99；SU2: 2.46）⇒ `Δ=−4<0`，对**一切** `x` 成立，钉不住任何东西（ABC3 §1）。只在 `E[H]+W≥2√E[H]`（`E=2` 需 `W≥0.83`，已证地板 `N/12=0.25`）才咬——达不到。
- **范数正性 / RP**：PSD ⇒ Cauchy–Schwarz/Hölder，**Schur-凹**，只给上界（SCHUR-CONVEX；INEQ-T5；GOLDS-POS-MASS）。文献明文："positivity alone forces the minimal residue to zero"（Caron-Huot–Pokraka–Zahraee JHEP 2024）。
- **Crossing**：纯 YM 非 CFT（`T^μ_μ≠0`，GLUEBALL-BOOTSTRAP §0）；且**自由 CFT 恰是 crossing 精确解**，永远排不掉（FORGE5-BOOTSTRAP-EXCLUDE Lemma）。
- **AF UV 尾**：钉 `s→∞` 渐近/总权，对 `x∼1/m₁` 的长距 pole 盲（ABC3 §3；PROOF-DISSOLUTION §3）。

## (2) 核心定理：矩是 ∞ 处 Taylor germ，edge 是有限处 pole
KL 测度的原子 `f₁²=ν({m₁²})` 是 Stieltjes 函数 `G(z)` 的**极点**；矩是其**∞ 处 Taylor germ**。从 Taylor 数据取有限内点留数 = 解析延拓，经典病态。**每个真定理只上界**：Christoffel–Darboux `f₁²≤1/K_n(x,x)`；Hankel PSD、Nevanlinna–Pick、Stieltjes 反演同向（INEQ-T1）。矩锥的下端**恰为 0 且被实现**：`δ₁+δ₃` 与 `[1,3]` 上 Lebesgue 密度同 `(M₀,M₁)=(2,4)`，edge 1 上原子 `1` vs `0`；`δ₀+δ₂` 与 `2δ₁` 同 `M₀=M₁=2`，edge `0` vs `1`（GAP-G2-edge）。故**任意有限多精确矩都不确定 `inf supp`**——"钉够矩以逼隙"失败的精确原因。

## (3) 可行性归约到 H1 吗？—— 否，严格更弱；逼隙才归约到 H1
Bootstrap 的 **primal feasible = 不排除**（要排除无隙点须 **infeasible**）。而约束系统**恒可行**：自由场（两光子 cut，绝对连续，`f₁²=0`，`M₀=1`）实现全部被护矩且无隙（SCHUR-CONVEX 主障碍）。故**可行性被一个无隙测度回答"是" ⇒ 零隙信息，严格弱于 H1**（H1⇒可行，反之不然）。要**排除一切无隙测度**，须加**非矩数据**：第二能级隙 `Λ>m₁²`（residue–level-gap 引理 `f₁²≥(ΛM₀−M₁)/(Λ−m₁²)`，GAP-G2-edge / INEQ-T1）。`Λ` 正是粗模 Poincaré 质量 = **H1 换名**。故：**可行 = 自明真（空洞）；逼隙 = H1。**

## (4) 具体闭合再加一层空
SD+Ward 的"自洽闭合"就是 `x≥0`（§1，`Δ<0`）。额外 `W_geo≥0`（Prékopa）只给 `N/12` 漂移，不吸收典型场 `Var(G|0)`（SUB1：`D≤cW_geo` 缺 ~10×）。故即便五条全上，闭合仍是套套逻辑。

## (5) 与 conformal bootstrap 的类比三处断裂
(i) YM 无 IR CFT（`T^μ_μ≠0`），无 conformal crossing 可解；(ii) 自由点 crossing 精确，排不掉；(iii) bootstrap 返回 **Δ 坐标的单侧 island/上界**，而 #2/#5 是 `k→0` 曲率坐标的**下界**——坐标与方向皆错（GLUEBALL-BOOTSTRAP §b；FORGE5 §B）；S-matrix bootstrap 更须先假定隙以定义渐近态。

## 五点记分
1. **非微扰 — PASS。** 精确谱表示 + SD + Ward，全 β。
2. **控号 `m_eff²(k→0)≥λ` — FAIL。** 全为上界机/矩；SD+Ward 二次式 `Δ<0` 空洞。
3. **非循环 — FAIL at closing。** 可行自明真（非循环，但空洞）；逼隙须加 `Λ` = H1，零松弛换名。
4. **非-SUSY + 败自由场 — PASS 但空洞。** 自由场 `f₁²=0` 恰是可行 crossing/正性解 ⇒ "败于 `f₁²=0`"与"排不掉无隙点"是**同一事实**（FORGE5）。
5. **均匀 `m_eff²(k→0)≥λ` — FAIL。** 矩不定 `inf supp`；`1/N²` 非均匀（ABC2；PROOF-LARGE-N-RIGOROUS）。

**最锐一句.** Bootstrap 是正性/矩机，把谱边原子**从上方**夹（`f₁²≤1/K_n`）；隙要的是支撑下确界的**下方**界。钉任意多精确矩，`inf supp≥0` 仍自由（自由场把它压到 0，`δ`/Lebesgue 纤维证明其可实现）。故 moment problem 的**可行性非空且隙盲**——严格弱于 H1、被无隙点满足；而**逼出隙**（排除无隙测度）恰是那条被排除的谱边数据 `Λ` = H1。Bootstrap 不闭合：它把"低-q 权重须载隙"换成"支撑边须被排除"，逐字是 H1 的矩语言。

*Refs.* ABC3-SPECTRAL-SUMRULE（`Δ=E[H](E[H]−4)<0`、`E[H]≈2.0`、`E+W≥2√E` 阈值）；ABC2-DUALITY（§1 谱表示=平面隙、sum rule 只钉总权）；ABC2-SYNTHESIS（§4 条件化=隙）；PROOF-IR-H1（H1=隙、rename）；INEQ-T1-hankel-stieltjes（矩=∞ 处 Taylor、Christoffel 上界、原子 0 vs CD 最大）；INEQ-T5 / INEQ-R10（正性只给平凡 `f₁²≥0`、正权 ⇒ 上界符号引理）；GAP-G2-edge（矩不定 `inf supp`、residue–level-gap 需 `Λ`）；SCHUR-CONVEX（自由场 `f₁²=0, M₀=1` 边界、coercivity=`H*`）；GOLDS-POS-MASS（正性只钉 `m=0`）；FORGE5-BOOTSTRAP-EXCLUDE（自由 CFT crossing 精确、class-P）；GLUEBALL-BOOTSTRAP（YM 非 CFT、SDPB 绑 `Δ` 非质量）；PROOF-DISSOLUTION-SPECTRAL-DENSITY（OPE 对 `s→∞` 盲）；k0-CONSTRAINT-EFFPOT / PROOF-JACOBIAN-SAVES（`G″=Δ_dyn+W_geo+X`、`W_geo≥N/12`）；SUB1-COV-BOUND（`D≤cW_geo` 缺 10×）；COMPUTE-GPP-SU3 / W7-Gpp0（`E[H]≈2.0`）；PROOF-LARGE-N-RIGOROUS（`1/N²` 不传输）；Källén 1952；Lehmann 1954；Caron-Huot–Pokraka–Zahraee, JHEP 01 (2024) 195 [arXiv:2309.04472]。

# C4-COV-LASTMILE — certifying `|X|<W_geo` without a variance bound is the `Var(G|0)` ceiling = H1; the rearrangement is a re-description, and the only variance-free route (a negative-correlation theorem) is absent *and* contradicted per-mode

> **判词 (verdict).** B2 把 β=2.5 的符号压到 `G″=W_geo+X+Δ_dyn`，`X=−2Cov(G,Φ′|0)`。**逐模实测**（COMPUTE-GPP-SU3/w7 log）`X` 取 `−0.09…+0.16`，均值≈0，`|X|≤0.16≪W_geo≈0.3`；而唯一严格界 `|X|≤2√(VarG·VarΦ′)≈0.40>W_geo` 短 `O(0.1)`。**裁决：这 `O(0.1)` 是 Cauchy–Schwarz 的松，不是物理的差——`X` 本就极小；而证明它小，逐字是 `Var(G|0)≤W_geo²/4VarΦ′≈0.8–1.3` vs 实测 `VarG≈2`，即 H1。** rearrangement/Chebyshev **不闭合**：SU(N)+Wilson 无 FKG（且 FKG 指向反号），fiber 无标量充分统计量，而不知 copula 的多元 rearrangement **恰退化为** CS；唯一免方差的逻辑是**符号**定理 `Cov≤0`（则 `X≥0`、`G″≥W_geo+Δ_dyn≥N/12`），但它要一条不存在的**负**关联不等式，且被逐模**变号**的数据否证。不是 margin 问题，是 certification-without-variance 问题；=(◆)=H1。

## §1 精确短差，与它的松
`G″=Δ_dyn+W_geo+X`，`W_geo=E[Φ″|0]−Var(Φ′|0)≥N/12`（Prékopa，rigorous）。开格上 `Δ_dyn≈0`（等分，`E[H]−Var≈0`），故 `G″≈W_geo+X`。CS：`|X|≤2√(VarG·VarΦ′)≈0.40`（`VarG≈2,VarΦ′≈0.02`），而 `W_geo≈0.32` ⇒ 下界 `−0.08<0`。但 log 逐模值（`−0.088,+0.158,−0.018,−0.059,+0.073,+0.047`）表明 `X` 是**零均值、按模变号的小量**：诊断 `corr(a,∇S)≈0`——作用力几乎与粗幅 `a` 正交，故与 `Φ′≈0.85a+噪声` 的正交投影几乎不协变。算术：若 `Φ′≈ca`，则 `E[aG]=1−cVar(a)`（Stein），`Cov(G,Φ′)≈c(1−cVar(a))≈0` 恰在**Jacobian 主导刚度**（`cVar(a)≈1`，即"平坦 Haar 供曲率"）时成立——`X` 小与 `W_geo` 主导是**同一件事**。**小的是 `X`，松的是界**。

## §1′ 对题面的直接回答
"rearrangement 闭合 `O(0.1)`，还是紧 `|X|` 即 `VarG`=H1？"——**两者都不是好结局，且是同一结局**：rearrangement **不**闭合（§2–§4），而任何**从方差出发**的紧 `|X|` 逐字就是 `VarG`=H1（§5）。唯一真正独立于 H1 的目标是**符号** `Cov≤0`（§3），它若成立确实直穿 H1，但既无不等式支撑、又逐模变号——所以题目要证的 "`Cov<0`" 不是一个自洽的普遍命题，而是一个**测量 artifact**：`X` 小是因为零均值/正交，不是因为恒负。

## §2 为何 Chebyshev/rearrangement 三重失效
1. **无 FKG，且方向反。** Harris/Chebyshev 求和不等式需格序+log-超模密度；`e^{−S_W}` 是 plaquette 相位的非单调余弦积，在平坦图卡上非 log-超模。更致命：FKG 给**正**关联 `Cov(G,Φ′)≥0`（⇒`X≤0`，减号），恰是要证的反面。要 `Cov≤0` 需**负**关联（D-相关/确定性点过程类），SU(N) Haar 无此已知结构。
2. **无标量充分统计量。** 一维 Chebyshev（增函数×减函数 ⇒ 负协变）只在两者是**同一标量**的函数时成立。"`|F|` vs 趋近中心元"的分裂只是**启发**：fiber 是 `2(N²−1)V` 维，`G`、`Φ′` 无公共全序。没有单变量约化。
3. **多元 rearrangement=CS。** Fréchet–Hoeffding：仅给边际，`Cov` 的极值由反单调耦合给出，模长即 `√(VarG·VarΦ′)`。想更紧必须输入 **copula**（测度本身）——即 H1。

## §3 符号路：唯一免方差的路，但其"机制"是图卡局部的
若 `Cov(G,Φ′|0)≤0` 则 `X≥0`，`G″≥W_geo+Δ_dyn≥N/12`，**不需任何方差界**，直穿 H1。所需是 SU(N)+Wilson 的负关联定理——§2.1 已述其不存在（无 FKG，无负关联）。所谓"log-cut 反局部化"机制（`Φ′` 在中心元附近大而 `F≈0`）**正是 PROOF-SAVVIDY-WGEO-DOMINATION 证明的图卡局部现象**：`Φ` 的超额凸性 `g″→∞` 只在 `d→2π`（`U≈z∈Z_N`，`U_□≈1`，`F≈0`）出现，`|F|`-无关且**与 Savvidy 驱动 `∝|F|` 反相关**。这不是测度的负关联性质，而是**单一坐标奇性的几何**；把它提升为协方差符号定理需它在**所有**高 `|F|` 区控场——而 §1′/log 显示它逐模变号（最高模 `X=+0.16`，最低 `−0.09`），β=5.7 均值 `Cov>0`（`X<0`），仅 β=2.5 Part-B 的 `k→0`-like 模 `X≈+0.22`。故"`Cov<0` 普遍成立"**按题面即假**；稳健事实只是 `|X|≪W_geo`（正交性），非定号。

## §4 二阶 Brascamp–Lieb 横项 —— 循环
BL 型协方差界 `Cov(G,Φ′)≤∫⟨B∇G,∇Φ′⟩` 欲优于 CS，须取最优矩阵 `B`；而 `B` 的正定性/量级由联合 Hessian `Hess V` 决定——β=2.5 处即 Savvidy 不定 `−D²+2adF`（`λ_min=−g|B|`）。取 `B` 就是求那台不定块之逆，逐字 H1（SCBL `G″=S−D` 的 `D>0` 即此）。若 `Hess V⪰0`（低场/β<1/48）BL 已给 `D≤0`、`G″≥S≥N/12`——**那是已闭格**，非本开格。二阶 BL 只在 `Hess V` **不定**处才 "new"，而不定处它退化为求逆问题。故二阶横项界**不含**一条独立于 `Var(G|0)` 的信息。

## §4′ 一个真正新物种长什么样（speculative，不在手）
要直穿 H1 又不用方差，唯一形态是**确定性、无耦合**的负相关性：对 Haar–Wilson 建立 `Cov(G,Φ′|0)≤0` 或 `≤ε·W_geo`（`ε<1`）而**不**经 `Var(G)`。类比的候选是 log-次模测度的 D-相关（Borcea–Brändén–Liggett）在 `SU(N)` 上的推广，或把 `(G,Φ′)` 映到**同一个**单调标量（如逐 link 的主相位 `θ_ℓ`）后用一维 Chebyshev——但 §2.2 已证无此标量，且 §3 证机制是图卡局部。这要求新几何不变量，不是重排现有分裂。

## §5 紧 `|X|` 逐字是 `VarG≤C`
`2√(VarG·VarΦ′)<W_geo ⟺ VarG<W_geo²/(4VarΦ′)`。取 `W_geo≈0.25–0.32`：阈 `0.78–1.28`，而 `VarG≈2`——差 `1.5–2.5×`。而 `Var_μ(G|0)≤C` 均匀正是 (◆)=`Var(a)≤1/λ`=H1（PROOF-JACOBIAN-SAVES §2）。**故"紧 `|X|` 界"不是比 H1 小的一块，就是 H1。**

## 五点裁决
1. **非微扰 — PASS（条件）。** 关联/rearrangement 是测度不等式，无圈展开。
2. **控号下界 — 条件 PASS / 严格 FAIL。** 符号路给 `G″≥N/12`（真下界，但需不存在的负关联）；方差路 `G″≥N/12−|X|` 不定号。
3. **非循环 — FAIL。** 紧 `|X|`⇔`Var(G|0)≤C`⇔H1；§4 循环。
4. **非-SUSY + 败自由场 — PASS。** 自由场 `F≡0⇒Φ′=0⇒X=0`，CS 与 rearrangement 都产 0，不造伪正性。
5. **均匀下界 — FAIL。** `|X|≤0.40>W_geo`；要均匀只余符号路（假）。

## Exact failing estimate
```
G″ ≥ W_geo − |X|,   |X| ≤ 2√(VarG·VarΦ′) ≈ 0.40 > W_geo ≈ 0.32  ⇒  G″ ≳ −0.08.
Closing needs  Var_μ(G|0) < W_geo²/(4VarΦ′) ≈ 0.8–1.3  (vs measured ≈2)
   ⟺ (◆) Var(G|0)−E[H|0]<W_full  ⟺ Var(a_{k→0})≤1/λ  =  H1.
Only escape:  Cov(G,Φ′|0)≤0  ⇒ X≥0 ⇒ G″≥N/12  — needs a NEGATIVE-association
   inequality for SU(N)+Wilson (no FKG; FKG gives Cov≥0, wrong sign);
   contradicted per-mode (X∈[−0.09,+0.16], mean≈0).
```

**Sharpest line.** The cross term is genuinely tiny — but because it is **near-orthogonal** (`corr(a,∇S)≈0`, `X` zero-mean per-mode), not because it is negatively correlated; so the task's "`Cov<0`" premise is neither the mechanism nor true per-mode. The `O(0.1)` is Cauchy–Schwarz's looseness (`VarG≈2`), and tightening it is exactly the `Var(G|0)` ceiling — H1 in curvature language. The one variance-free path, a sign theorem `Cov≤0`, would close it cleanly but requires a negative-association inequality SU(N) does not have and the data contradict. **The rearrangement does not close the last mile; a tight `|X|` bound *is* `VarG`=H1.** What would be genuinely new: a deterministic, coupling-free negative-correlation / log-submodular comparison for the Haar–Wilson measure — not a rearrangement of the existing split.

*Refs.* COMPUTE-GPP-SU3 (`W=+0.28/+0.31`, `X=−0.04/−0.01`, `G″=+0.42/+0.33`; diagnostics `corr(a,∇S)≈0`, `corr(a,Φ′)≈0.85`); w7/su3_L4.log (per-mode `X∈[−0.09,+0.16]`); ABC2-CONSTRUCTIVE (§4 `|X|≤2√VarG·VarΦ′≈0.40>0.32`); PROOF-JACOBIAN-SAVES (§1–2, (◆)=H1); OTHER-H1-POINCARE (`G″=Δ_dyn+W_geo+X`); OTHER-GPP-KZERO (measured +0.51, `Δ_dyn≈0`, (◆)); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`, cross-`B`=Hess); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`); PROOF-SAVVIDY-WGEO-DOMINATION (log-cut anti-correlated, `ρ≳10`); SUB1-COV-BOUND (`D≤cW_geo` false); Fortuin–Kasteleyn–Ginibre 1971; Harris 1960; Chebyshev sum ineq.; Fréchet 1951/Hoeffding 1940; Rüschendorf 1983; Brascamp–Lieb 1976; Prékopa 1971; Savvidy PLB 71 (1977) 133.

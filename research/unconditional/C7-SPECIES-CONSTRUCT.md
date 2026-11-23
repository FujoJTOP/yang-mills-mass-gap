# C7-SPECIES-CONSTRUCT — SHAPE C round 6: five mechanisms outside the 9-class table are five dialects of one wall — the variance route is dead by equipartition, the sign route needs a negative-association/orthogonality theorem none supplies, and the observed near-orthogonality is the Ward/gap relation in disguise

> **判词 (verdict).** (V) 方差上界 与 (N) 负号 **不对称**。(V) 被**等分律**封死：`Var(G|0)≈E[H|0]≈2`（B4：L4 2.01→L6 2.05，体积稳），则 `G″≥E[H]−VarG≈0`，**永不给 `≥λ>0`**——上界路线是**恒等式形状**的死路，非缺一条不等式。(N) 是唯一活路（`X=−2Cov≥0⇒G″≥W_geo+Δ_dyn≥N/12`，免方差直穿），但 (a)–(e) **皆不在** 9 类（NA/NR-ND/strong-Rayleigh/log-submod/DPP/BL/FKG/Simon–Lieb/Chebyshev）内，却**无一**造出 `e^{−S_W}` 上的负关联/正交性定理。**(a)** "IR 条件测度更 log-凹"逐字是 SCBL 缺陷 `D_k≤0`；β=2.5 的 `D>0` 即不定纤维 `K_k`（W6-BE 90–97% 负向在被积 band 内；`PᵀHP|ran P` −0.3），"条件化整合负 band"**正是** D 的生成，唯 `k→0` 处干净、那里 `D_k≤0`=`G″(0)≥λ`=H1。**(b)** 可证号仅 `β<1/48`（P3），那是残数**直接**正性（或 FKG **正**号=错号），非负关联；`f₁²` 跨 `[1/48,2.5]` 稳定性=β-非解析跳变 no-go，且 NA 在 2.5 **逐模变号**（`X∈[−0.09,+0.16]`）**不存在**，无极限可取。**(c)** HWI `W_2²≤(2/λ)·H`，相对熵 `H(μ_{2.5}|μ_0)∼βV` **广延**：欲得小 `W_2` 需 `λ=O(1)`=隙，否则 `W_2∼√V`；扰动 `Hess(βS_W)` 为 Savvidy 不定（`λ_min=−g|B|`≫`N/12`）；给正曲率的传输不等式常数即 `1/λ`=LSI=隙，反单调耦合=Fréchet–Hoeffding=C-S `0.40`。**(d)** 模量界是**逐点值**乘子= B3 coercion 类，`ρ` 不变（零乘子即隙），规范商去其杠杆（Savvidy 负模物理、过商仍在；Landau FP 零模/Gribov），标量模量违反 C4 §2.2"无公共充分统计量"。**(e)** `corr(a,∇S)≈0` 经 Ward `E[af]=1`、`Φ′≈ca` **等价于** `E[aG]=1−E[aΦ′]≈0 ⟺ cVar(a)≈1`，即 `G″Var(a)=1` 一族=H1；Elitzur 只杀单个不变算子、不杀叉协方差，规范不变性给 Ward **恒等式**（重述而非界定）。**裁决：五种皆同一堵墙的方言；方差死于等分，符号无人供给。五点：五者皆 1 过、2 败、3 败、4 过、5 败。**

## The two routes are not symmetric
`G″=Δ_dyn+W_geo+X`, `Δ_dyn=E[H]−Var(G)`, `W_geo=E[Φ″]−Var(Φ′)≥N/12` (Prékopa), `X=−2Cov(G,Φ′|0)`. An **upper** bound `Var(G|0)≤C` acts through `Δ_dyn=E[H]−Var≥E[H]−C`; equipartition pins `Var(G|0)≈E[H|0]≈2`, so `Δ_dyn≈0` and *no* `C` yields `G″≥λ>0`. (V) is dead by identity shape, not a missing estimate. Only (N) — a **negative-sign** theorem (`Cov≤0`, or `|Cov|≤ε·W_geo`, `ε<1`) — adds to `W_geo` and crosses without variance. The 9-class table gave no (N); (a)–(e) were the remaining non-member candidates, and they fail the same way.

## (a) RG / conditioned log-concavity *is* the SCBL defect
"More log-concave at the IR scale" is exactly `D_k≤0` in `G″=S_k−D_k` (INVENT-VARIANCE-INEQ), `D_k=Var(f|0)−E[⟨K_k^{−1}c,c⟩|0]`, `K_k` the Hessian of the band integrated to `k`. But `D>0⟺K_k` indefinite, and W6-BE places 90–97% of negatives *in that band* (`PᵀHP|ran P` −0.3, OTHER-H1): conditioning **over** the negative band creates `D>0`, it does not remove it. The only cut where the conditional is provably log-concave is `k→0`, where `D_k≤0` is verbatim `G″(0)≥λ`. Bankable: *"IR convexification ⟺ SCBL defect ≤0"*; missing: a contraction of `D_k/S_k` under blocking — itself the emergence of IR stiffness.

## (b) Lattice association + stability limit
The only rigorous positive sign is P3 (`f₁²>0`, `β<1/48`, cluster expansion) — *direct residue positivity*, not negative association; strong-coupling FKG-type inequalities, where they hold, give the **opposite** sign. The needed step — `f₁²>0` stable across `[1/48,2.5]` — is the documented β-nonanalytic-jump no-go, and the NA it would carry is false at 2.5 (per-mode `X` changes sign): there is **no limit object**. Passes only where already banked.

## (c) Transport coupling — only C-S or the LSI constant
HWI gives `W_2²≤(2/λ)·H`, with relative entropy `H(μ_β|μ_0)∼βV` **extensive**: a small, uniform `W_2` perturbation requires `λ=O(1)` — i.e. the gap — else `W_2∼√V`; the curvature perturbation is `Hess(βS_W)`, the Savvidy indefinite form exceeding `N/12`. Any transport/HWI inequality that *yields* positive curvature has constant `1/λ`=LSI=H1; an anti-monotone (Fréchet–Hoeffding) coupling is C-S `0.40`, short. Kendall/Mallows are permutation-space metrics, inapplicable to continuous gauge fields.

## (d) Gauge-fixed modulus bound — a re-skin of B3
A modulus `M(U)` is a *pointwise value*; bounding it rescales `Var(G)`, `Var(Φ′)` together, leaving `ρ=Cov/√(VarG·VarΦ′)` — and `G″` — invariant ("the null multiplier is the gap"). Gauge-fixing does not help: Savvidy negative fiber modes are physical (`δF≠0`, survive the quotient, PROOF-SAVVIDY-WGEO-DOMINATION); Landau FP zero modes/Gribov copies defeat any `M`-lower-bound exactly at the soft mode; a gauge-invariant scalar modulus is delocalized/non-monotone — no common scalar sufficient statistic (C4 §2.2).

## (e) Near-orthogonality via Elitzur *is* the Ward relation
`f=G+Φ′`, Ward `E[af]=1` ⇒ `E[aG]=1−E[aΦ′]`; with `Φ′≈ca` (`corr(a,Φ′)≈0.85`), `Cov(G,Φ′)=c(1−cVar(a))`. Hence `corr(a,∇S)≈0 ⟺ E[aG]≈0 ⟺ cVar(a)≈1` = the Ward/gap relation `G″Var(a)=1`. Elitzur zeroes *single* non-invariant operators (`⟨a⟩=⟨G⟩=0`), not cross-covariances; gauge invariance gives Ward *identities* that **recast** `Cov(G,Φ′)` via the gap, never bound it. The observed orthogonality is the gap in disguise, not an independent structural fact a gauge argument can certify.

## Five-point scorecard
1. **Nonperturbative — PASS** (all five measure-level; no loops/`1/N`/`g`).
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL** (V dead by equipartition; N absent; per-mode sign-changing).
3. **Non-circular — FAIL** ((a) `D_k≤0`≡H1 at `k→0`; (b) β-stability≡no-go; (c) const=`1/λ`; (d) value multiplier; (e) Ward≡`G″Var(a)=1`).
4. **Non-SUSY + fails free field — PASS** (free `G=0,Φ′=0⇒Cov=0,G″=S=0`; (c)/(d) vacuous there).
5. **Uniform `≥λ` — FAIL** (none produces a `λ`).

## Exact failing estimate
```
(V) Var(G|0)≤C<E[H|0]≈2 : impossible — equipartition Var(G|0)≈E[H|0] (L4 2.01, L6 2.05)
    ⇒ G″ ≥ E[H]−VarG ≈ 0, never ≥λ.        [variance route dead by identity shape]
(N) Cov(G,Φ′|0)≤0.11 (⟺ X≥−0.22): C-S ⇒ ≤0.20; four association hypotheses absent
    (continuous, frustrated cos-sums, non-attractive ⇒ neither real-stable nor ordered).
(a) IR-convex ⟺ D_k≤0 ; D>0 = indefinite K_k (90–97% negatives in band; PᵀHP|ranP −0.3).
(b) provable sign only β<1/48 (P3); β-stability across [1/48,2.5] = β-nonanalytic no-go; NA false per-mode.
(c) HWI W_2²≤(2/λ)H, H(μ_β|μ_0)∼βV extensive ⇒ small W_2 needs λ=O(1)=gap; Hess-pert λ_min=−g|B|>N/12; positive-curvature transport const=1/λ=H1.
(d) modulus bound = value multiplier, ρ-invariant ("null multiplier = gap"); Gribov breaks it at soft mode.
(e) E[aG]=1−E[aΦ′]≈0 ⟺ cVar(a)≈1 ⟺ G″Var(a)=1; Elitzur kills operators, not Cov(G,Φ′).
```

**Sharpest line.** No new species: the five are five write-ups of the one exact missing step (◆), `Var_μ(G|0)` uniform, equivalently `G″Var(a)=1` at `k→0` = H1. The variance route is dead by *identity shape* (equipartition `Var≈E[H]`, so `Δ_dyn≈0`; no upper bound can be positive), and the sign route — the only survivor — needs a deterministic, coupling-free **negative-correlation/orthogonality** theorem for `e^{−S_W}`, a measure that is continuous, frustrated (cos-sums), non-attractive, carrying no order, no real-stable polynomial, no DPP kernel, no chart-local common scalar. Each mechanism controls *size* (`VarG`, `|X|` via C-S), never the *angle* `ρ` that carries the sign. Two results are worth banking: **(a)** "IR log-convexification ⟺ SCBL defect `D_k≤0`" (a conditional identity; numerically `D/S≈2%`, but no RG monotone), and **(e)** the proof that `corr(a,∇S)≈0` is *provably* the Ward/gap saturation `cVar(a)≈1` — the orthogonality is not a free structural gift but the gap re-expressed. The last step remains: a negative-association species for the continuous frustrated non-attractive SU(N) Wilson measure, which does not yet exist.

*Refs.* C4-COV-LASTMILE, C5-NEG-ASSOC, C6-DEPENDENCE-CLASS, B4-VARIANCE-MAG, B3-TYPICAL-COERCIVITY; ABC-NEW-INEQ (Ward `G″Var(a)=1`); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`, `D/S≈2%`); PROOF-SAVVIDY-WGEO-DOMINATION (`λ_min=−g|B|`, negatives physical); OTHER-H1-POINCARE (`PᵀHP|ran P` −0.3, blocking keeps k→0); W6-BE-at-2.5; COMPUTE-GPP-SU3/w7 (per-mode `X`, `corr(a,∇S)≈0`, `corr(a,Φ′)≈0.85`). External: Elitzur 1975; Ward/Schwinger–Dyson; Fréchet 1951/Hoeffding 1940; Villani 2009; Mallows 1972; Kendall 1938; Bakry–Émery; Fortuin–Kasteleyn–Ginibre 1971; Prékopa 1971; Brascamp–Lieb 1976; Savvidy PLB 71 (1977) 133.

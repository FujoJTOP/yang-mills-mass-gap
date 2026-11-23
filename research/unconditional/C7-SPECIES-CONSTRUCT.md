# C7-SPECIES-CONSTRUCT — SHAPE C round 6: five mechanisms outside the 9-class table are five dialects of one wall — the variance route is dead by equipartition, the sign route needs a negative-association/orthogonality theorem none supplies, and the observed near-orthogonality is the Ward/gap relation in disguise

> **判词 (verdict).** 题设要的两样——(V) 方差上界 `Var(G|0)≤C` 或 (N) 负号 `Cov(G,Φ′)≤0`——**不对称**。(V) 被**等分律**结构性封死：`Var(G|0)≈E[H|0]≈2`（B4：L4 2.01→L6 2.05，体积稳），于是 `G″≥E[H]−VarG≈0`，**永远给不出 `≥λ>0`**——上界路线是恒等式形状的死路，不是缺一条不等式。(N) 是唯一活路（`X=−2Cov≥0⇒G″≥W_geo+Δ_dyn≥N/12`，免方差直穿），五个新机制 (a)–(e) **都不在** 9 类（NA/NR-ND/strong-Rayleigh/log-submod/DPP/BL/FKG/Simon–Lieb/Chebyshev）内，却**无一**造出 `e^{−S_W}` 上的负关联或正交性定理。**(a) RG 凸化 = SCBL 缺陷 `D_k≤0`**：逐字同一方程；β=2.5 的 `D>0` 就是不定纤维 `K_k`（W6-BE 把 90–97% 负向放进被积掉的band；`PᵀHP|ran P` min≈−0.3），条件化"整合高动量负向"**正是** `D` 的生成；只在 `k→0` 处干净，那里 `D_k≤0` 即 `G″(0)≥λ`=H1。**(b) 格点关联+极限**：可证区仅 `β<1/48`（P3 已闭），能得到的是**正** FKG 号（错号）或残数直接正性，非所需负关联；`f₁²` 跨 `[1/48,2.5]` 的 β-稳定性 = "β-非解析跳变"no-go，且 NA 假设在 2.5 **逐模变号（`X∈[−0.09,+0.16]`）根本不存在**，无物可取极限。**(c) 传输**：`W_2(μ_{2.5},μ_0)≳√(βV)`（广延相对熵），扰动 Hessian 为 Savvidy 不定（`λ_min=−g|B|`≫`N/12`），任何给出正曲率的传输不等式常数为 `1/λ`=LSI=隙，反单调耦合退化为 Fréchet–Hoeffding=C-S `0.40`。**(d) 规范固定模量**：模量界是**逐点值**乘子=B3 coercion 类，`ρ` 不变（"零乘子即隙"），规范不变性/商去其杠杆（Savvidy 负模为物理、过商仍在；Landau FP 零模/Gribov），标量模量违反 C4 §2.2"无公共标量充分统计量"。**(e) 近正交性**：`corr(a,∇S)≈0` 经 Ward `E[af]=1` 与 `Φ′≈ca` **恰等价于** `E[aG]=1−E[aΦ′]≈0 ⟺ cVar(a)≈1`，即 `G″Var(a)=1` 一族=H1；Elitzur 只杀单个非不变算子、不杀叉协方差，规范不变性给的是 Ward **恒等式**（重述而非界定）。**裁决：五种新机制皆同一堵墙的方言；方差死于等分，符号无人供给。五点：五者皆 1 过、2 败、3 败、4 过、5 败。**

## The two routes are not symmetric
`G″(0)=Δ_dyn+W_geo+X`, `Δ_dyn=E[H]−Var(G)`, `W_geo=E[Φ″]−Var(Φ′)≥N/12` (Prékopa), `X=−2Cov(G,Φ′|0)`. A variance **upper** bound `Var(G|0)≤C` enters through `Δ_dyn=E[H]−Var≥E[H]−C`; since equipartition pins `Var(G|0)≈E[H|0]≈2`, `Δ_dyn≈0` and no `C` gives `G″≥λ>0`. So (V) is dead *by identity shape*, not by a missing estimate. Only (N), a **negative-sign** theorem (`Cov≤0` or `|Cov|≤ε·W_geo`, `ε<1`), adds to `W_geo` and crosses without variance. The 9-class table already showed none of its members supplies (N); (a)–(e) were the remaining non-member candidates. They fail the same way.

## (a) RG / conditioned log-concavity — this *is* the SCBL defect
"Make the conditioned coarse measure more log-concave at the IR scale" is exactly `D_k≤0` in `G″=S_k−D_k` (INVENT-VARIANCE-INEQ), where `D_k=Var(f|0)−E[⟨K_k^{−1}c,c⟩|0]` and `K_k` is the Hessian of the band integrated out down to `k`. But `D>0 ⟺ K_k` indefinite, and W6-BE puts 90–97% of negatives *in that band*, while `PᵀHP|ran P` min≈−0.3 (OTHER-H1): conditioning **over** the negative band is what *creates* `D>0`, it does not remove it. The only cut where the conditional is provably log-concave is `k→0`, where `D_k≤0` is verbatim `G″(0)≥λ`. Bankable structural fact: *"IR convexification ⟺ SCBL defect ≤0"*; missing: a contraction of `D_k/S_k` under blocking — which would itself be the emergence of IR stiffness.

## (b) Lattice association + stability limit — the provable set is β<1/48
The only rigorous positive sign is P3 (`f₁²>0`, `β<1/48`, cluster expansion). That is *direct residue positivity*, not the negative association needed; strong-coupling FKG-type inequalities, where they hold, are the **opposite** sign. The needed step — `f₁²>0` stable across `[1/48,2.5]` — is exactly the documented β-nonanalytic-jump no-go (no crossing lemma), and the NA hypothesis it would carry is false at 2.5 (per-mode `X` changes sign), so there is **no limit object** to take. Passes only in the region already banked.

## (c) Transport coupling — only C-S or the LSI constant
`W_2(μ_β,μ_0)` is controlled by relative entropy `β(E[S_W]−min)∼βV`, hence `≳√V`; the curvature perturbation is `Hess(βS_W)`, the Savvidy indefinite form (`λ_min=−g|B|`) exceeding the `N/12` floor. Any transport/HWI inequality that *yields* positive curvature has constant `1/λ` = LSI = H1; an anti-monotone (Fréchet–Hoeffding) coupling is C-S `0.40`, short. Kendall/Mallows metrics are for permutation spaces, inapplicable to continuous gauge fields.

## (d) Gauge-fixed modulus bound — a re-skin of B3
A modulus `M(U)` is a *pointwise value*; bounding it multiplies `Var(G)` and `Var(Φ′)` together, leaving `ρ=Cov/√(VarG·VarΦ′)` — hence `G″` — invariant ("the null multiplier is the gap", B3). Gauge-fixing does not help: Savvidy negative fiber modes are physical (`δF≠0`, survive the quotient, PROOF-SAVVIDY-WGEO-DOMINATION), and Landau-gauge FP zero modes/Gribov copies defeat any `M`-lower-bound exactly at the soft mode. A gauge-invariant scalar modulus is delocalized and non-monotone — no common scalar sufficient statistic (C4 §2.2).

## (e) Near-orthogonality via Elitzur — it *is* the Ward relation
With `f=G+Φ′`, Ward `E[a f]=1` (periodic, no boundary) gives `E[aG]=1−E[aΦ′]`; since `Φ′≈ca` (`corr(a,Φ′)≈0.85`, C4), `Cov(G,Φ′)=c·E[aG]=c(1−cVar(a))`. Hence `corr(a,∇S)≈0 ⟺ E[aG]≈0 ⟺ cVar(a)≈1`, i.e. exactly the Ward/gap relation `G″Var(a)=1`. Elitzur's theorem zeroes *single* non-invariant operators (`⟨a⟩=⟨G⟩=0`), not their cross-covariance; gauge invariance supplies Ward *identities* that **recast** `Cov(G,Φ′)` in terms of the gap, never bound it. The measured orthogonality is the gap in disguise, not an independent structural fact a gauge argument could certify.

## Five-point scorecard
1. **Nonperturbative — PASS** (all five are measure-level; no loops/`1/N`/`g`).
2. **Sign-controlling lower bound `m_eff²(k→0)≥λ` — FAIL** (V dead by equipartition; N absent in all five; per-mode sign-changing).
3. **Non-circular — FAIL** ((a) `D_k≤0`≡H1 at `k→0`; (b) β-stability≡no-go; (c) constant=`1/λ`; (d) value multiplier; (e) Ward≡`G″Var(a)=1`).
4. **Non-SUSY + fails free field — PASS** (free `G=0,Φ′=0⇒Cov=0,G″=S=0`; (c)/(d) vacuous there).
5. **Uniform `≥λ` — FAIL** (none produces a `λ`).

## Exact failing estimate
```
(V) Var(G|0)≤C<E[H|0]≈2 : impossible — equipartition Var(G|0)≈E[H|0] (L4 2.01, L6 2.05)
    ⇒ G″ ≥ E[H]−VarG ≈ 0, never ≥λ.          [variance route dead by identity shape]
(N) Cov(G,Φ′|0)≤0.11 (⟺ X≥−0.22): C-S ⇒ ≤0.20; the four association hypotheses absent
    (continuous, frustrated cos-sums, non-attractive ⇒ neither real-stable nor ordered).
(a) IR-convex ⟺ D_k≤0 ; D>0 = indefinite K_k (90–97% negatives in band W6-BE; PᵀHP|ranP −0.3).
(b) provable sign only β<1/48 (P3); β-stability across [1/48,2.5] = β-nonanalytic no-go; NA false per-mode.
(c) W_2(μ_β,μ_0)≳√(βV); Hess-perturbation λ_min=−g|B|>N/12; positive-curvature transport const=1/λ=H1.
(d) modulus bound = value multiplier, ρ-invariant ("null multiplier = gap"); Gribov breaks it at the soft mode.
(e) E[aG]=1−E[aΦ′]≈0 ⟺ cVar(a)≈1 ⟺ G″Var(a)=1; Elitzur kills operators, not Cov(G,Φ′).
```

**Sharpest line.** No new species: the five are five write-ups of the one exact missing step (◆), `Var_μ(G|0)` uniform, equivalently `G″Var(a)=1` at `k→0`, = H1. The variance route is dead by *identity shape* (equipartition `Var≈E[H]`, so `Δ_dyn≈0` and no upper bound can be positive), and the sign route — the only survivor — needs a deterministic, coupling-free **negative-correlation / orthogonality** theorem for `e^{−S_W}`, which is a continuous, frustrated (cos-sum), non-attractive measure carrying no order, no real-stable polynomial, no DPP kernel and no chart-local common scalar. Each mechanism controls *size* (`VarG`, `|X|` via C-S), not the *angle* `ρ` that carries the sign. Two results are worth banking: **(a)** the equivalence "IR log-convexification ⟺ SCBL defect `D_k≤0`" (a conditional/structural identity, numerically `D/S≈2%` but with no RG monotone), and **(e)** the proof that the observed `corr(a,∇S)≈0` is *provably* the Ward/gap saturation `cVar(a)≈1` — so the orthogonality is not a free structural gift but the gap itself, re-expressed. The last step remains: a negative-association species for the continuous frustrated non-attractive SU(N) Wilson measure, which does not yet exist.

*Refs.* C4-COV-LASTMILE, C5-NEG-ASSOC, C6-DEPENDENCE-CLASS, B4-VARIANCE-MAG, B3-TYPICAL-COERCIVITY; ABC-NEW-INEQ (Ward `G″Var(a)=1`); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`, `D/S≈2%`); PROOF-SAVVIDY-WGEO-DOMINATION (Savvidy negatives physical, `λ_min=−g|B|`); OTHER-H1-POINCARE (`PᵀHP|ran P` min −0.3, blocking keeps k→0); W6-BE-at-2.5 (90–97% negatives in UV fiber); COMPUTE-GPP-SU3 / w7 log (per-mode `X∈[−0.09,+0.16]`, `corr(a,∇S)≈0`, `corr(a,Φ′)≈0.85`). External: Elitzur 1975; Ward/Schwinger–Dyson; Fréchet 1951/Hoeffding 1940; Villani, *Optimal Transport* (2009); Mallows 1972; Kendall 1938; Bakry–Émery; Fortuin–Kasteleyn–Ginibre 1971; Prékopa 1971; Brascamp–Lieb 1976; Savvidy PLB 71 (1977) 133.

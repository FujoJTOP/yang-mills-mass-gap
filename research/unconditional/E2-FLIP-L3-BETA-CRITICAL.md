# E2-FLIP-L3-BETA-CRITICAL — can a β-critical **non-analytic** functional lower-bound `m₁` *without* the nonperturbative crossover transport (construct it), or is **L-3 ≡ transport≡confinement** (prove it)?

> 2026-11-24 · Direct attack on the **L-3 loophole** of `E2-FLIP-NOGO §5` (and the same object flagged in `CLAUDE.md §5` / `BRIEF §4`): *"a functional that is β-critical — NON-ANALYTIC at the Savvidy crossover β≈2.5 — evades the deformation-rigidity no-go (horn β / B3), because a functional that JUMPS across β is not deformation-invariant."* The task: **(a)** CONSTRUCT such a `Φ` that is acid-clean, `(L,a)`-uniform, lower-bounds `m₁`, and is genuinely **not** "transport≡confinement" (a β-non-analytic order parameter / condensation onset / anomaly-jump / critical-point quantity that lower-controls `m₁` without solving the crossover); **OR (b)** PROVE `L-3 = transport≡confinement`. This note tries (a) hard, item by item, and proves (b).
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[V-here]** fresh computation (`research/unconditional/e2_flip_l3_witness.py`). External facts are **[R]** (standard/as-reported, not re-verified this session). **No QED on anything unproven. This does NOT prove Clay and does NOT prove no tool can exist** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`; `CONVERGENCE-redteam §d`).

---

## Verdict up front

> ### **=TRANSPORT≡CONFINEMENT. [NG] for the class `L-3` (β-critical / non-analytic functionals). (b) is proved; no (a) constructed.** The escape's premise is, up to one genuine sub-fact, a **category error**: L-3 fuses **two different "β-non-analyticities"**, and neither is a handle that yields a lower bound without the transport machine.
>
> 1. **Non-analytic *at `g=0`*** (the true one): `m₁(g,a) = c·Λ(g,a)`, `Λ = a^{-1}f(g)`, `f(g) ∼ e^{-c'/g²}` **flat at `g=0`** (essential singularity, β=∞). This is **real** — it is dimensional transmutation — but it is the **SCALE `Λ`**, carrying **no information about the pure ratio `c = m₁/Λ`**. `[T]` dimensional.
> 2. **Non-analytic at an intermediate *crossover `β_c≈2.5`*** (the one that would evade rigidity / horn β): **does not exist in the T=0 theory.** `PROOF-NO-MELT-PHASE` `[T-cond.]`: on a single-phase T=0 interval the gap is **real-analytic** in β; `β_c≈2.47` is the *analytic* running-coupling crossover (`λ≈7.3↔β_eff=18/λ`) or the **finite-`N_t` deconfinement** (a genuine `T>0` transition, in the **center / 1-form** channel, not the `0⁺⁺` residue). A "jump at `β_c≈2.5`" has no T=0 jump to lock onto; the finite-T one is the wrong channel and points the wrong way (§5.4).
>
> **The decisive structural reason (why β-criticality is the wrong lever):** *non-analyticity evades only the **rigidity** horn, which was never the binding one.* The binding obstruction to a **lower** bound is the **Schur / direction** wall, and the moment lemma (`GAP-G2-edge`, `SCHUR-CONVEX`) is **analyticity-blind**: a functional built from **aggregates** (moments, moments-of-moments, any `∫φ dν`, the free energy, a condensate, an anomaly) is **Schur-concave**, hence bounds the edge **from above**, *no matter how wildly non-analytic in β it is*. (`E2-FLIP-DYNAMICAL` already exhibited a fully coupling-dependent, non-invariant functional — the Savvidy Hessian — that **evades** horn β and still fails direction + uniformity.) So β-criticality moves `Φ` across the *rigidity* sign, not the *Schur* sign.
>
> **Two clean logical lemmas close it.** **(Far-side lemma [T-logic])** a *discontinuous* minorant must be `≤ m₁` on the **physical side**; its admissible jump-amplitude is exactly `inf_{phys} m₁`, i.e. the far-side gap — the seam supplies no calibration on either side (§3, `[V-here]` W3). **(Scale-vs-ratio [T]+[J])** the β-critical content is at most `Λ`; a lower bound needs the pure ratio `c` *(or the isolation margin `(ρ²−1)m₁²`)*, both `= ` the gap up to a universal factor; controlling either across `β_c` **is** the transport machine (§4).
>
> **Exact stuck step:** *the only β-critical Schur-**convex** minorant of `m₁` is `m₁`'s own singularity (the atom / the dissolution `Λ−m₁²→0`); and no such object is simultaneously **β-monotone/protected** (`PROOF-TRANSPORT-MONOTONE`: atom-coercive ⟂ β-monotone — opposite sides of the Schur order). Turning the β-critical datum into `m₁ ≥ Δ` therefore requires the crossover transport `= confinement`.*

---

## §0 The object, the class `L-3`, and what a flip must deliver

**Setup** (`E2-FLIP-NOGO §0`, `BRIEF §1`). Lattice Wilson measure `μ_β`, `β=2N/g²`, transfer operator `T_β=e^{-aH_β}`, gap `m₁(β)=E₁−E₀`, KL spectral measure `ν` of `O=¼tr F²`. Physical theory = **continuum limit** (`β→∞` *with* `a(β)→0`, `m_phys=c·Λ` fixed). Strong-coupling `β<1/48`: SZZ/P3 gap `[T]`. Free massless = Gaussian at `g=0` (β=∞, *fixed a*), gapless `[T]`.

**Class `L-3` (stated precisely — the target of this note).** A functional `Φ` on the β-family `{μ_β}` that
- **(LB)** `Φ(μ_β) ≤ m₁(μ_β)` for every `β` on the physical trajectory, uniformly in `(L,a)`, with `inf Φ ≥ Δ>0`;
- **(acid)** `Φ(μ_free)=0` (or `≤0`);
- **(β-crit)** is **non-analytic at `β_c`** (this is the proposed escape of horn β / rigidity B3);
- **(non-transport)** is computable **without** the nonperturbative RG transport that carries a strong-coupling fact to `β_c`/continuum.

By `E2-FLIP-NOGO` **Theorem 0 [T]**, (LB)+(acid) is *exactly* "a uniformly-positive **minorant** of the gap `β↦m₁(β)`, computable from accessible data." So the entire question is: **can a β-critical datum certify such a minorant without transport?** (a) = yes; (b) = no.

**What "transport≡confinement" means here** (`BRIEF §4`, `SMALLCIRCLE-ADIA §3`, `MONOTONE-COUPLING`, `PROOF-TRANSPORT-MONOTONE`). The controllable end (strong coupling / the more-gapped end in lattice units) is separated from the physical point by the crossover; *every* monotone comparison points toward the more-gapped end, so proven values are **upper** bounds (`MONOTONE-COUPLING §2`); local transport is a theorem **given** uniform isolation, and uniform isolation is `H*(r₀)` one rung above the gap (`PROOF-TRANSPORT-MONOTONE`). "The nonperturbative machine that carries the strong-coupling gap *continuously* across the Savvidy crossover and takes the continuum limit" is the missing object (`CLAUDE.md §5`).

---

## §1 The disambiguation that decides L-3 — *two* β-non-analyticities, neither usable

This is the core of the verdict and the sharpening over `E2-FLIP-NOGO §5`'s one-paragraph dismissal.

**(A) Non-analyticity at `g=0` — real, but it is the SCALE, not the gap.** `m₁(g,a) = c·Λ(g,a)` with `Λ a = f(g) ∼ C\,g^{...}\exp(-1/(2b₀g²))` — a **flat** function of `g²` at `0` (all derivatives `0`, essential singularity at `β=∞`). This is the genuine "β-non-analytic" content of the gap: **dimensional transmutation.** It is `[T]` (the standard form of asymptotic freedom / transmutation; `W3-adjacent-4D`, `MS-GAP-EQUATION`). But it yields **`Λ`, a dimensionful scale, up to scheme** — and a lower bound `m₁≥Δ` requires the **dimensionless ratio `c=m₁/Λ`** (`=1.7` for SU(3) glueball/`Λ`, the same pure number by single-scale covariance `G2D-ISOLATION-I2 §2`). The non-analytic datum gives the exponent `c'` (the one-loop coefficient) and the location `Λ`; it gives **nothing** about `c`. §4.

**(B) Non-analyticity at an intermediate `β_c≈2.5` — the rigidity-evading one — is absent at T=0.** `PROOF-NO-MELT-PHASE` `[T-cond. on single-phase]`: on a single connected T=0 confined phase, every correlator and `m₁` are **real-analytic** in `β` (finite `V`: entire kernel; infinite `V`: Kato + no Lee–Yang/Fisher zeros). `β_c≈2.47` is *(i)* the λ≈7.3 running-coupling crossover (`β_eff=18/λ`, an *analytic* finite-slope crossover) and/or *(ii)* the SU(2) **finite-`N_t`** deconfinement (`β_c(4)≈2.30…β_c(8)≈2.51`) — a real Svetitsky–Yaffe transition, but on a **finite temporal circle** (`T>0`), gone as `N_t→∞`. A *continuous* dissolution (`f₁²→0`, `m₁` rising into its threshold `Λ−m₁²→0`) is **analytic** and allowed (`PROOF-NO-MELT-PHASE §3`, `PROOF-TRANSPORT-CONTINUITY §3`) — **not a jump.**

> **The tension the task flags is thereby dissolved.** The phrase *"the gap is a β-non-analytic jump"* (from the rigidity no-go, `INVENT-NEW-INDEX §3`) is a statement about what an **index** would need (`B3`): an index evading rigidity would have to jump at finite β, but jumps at regular couplings are non-index, *built from the gap data*. It is **not** a literal claim that `m₁(β)` jumps at `β_c`. `PROOF-NO-MELT-PHASE` is the sharper truth: `m₁` is analytic at T=0 and non-analytic only at `g=0` — i.e. the "natural home" for the gap's non-analyticity is **(A) the transmutation scale `Λ`**, which is the **dimensionful** datum (wrong thing, §4), not (B) a T=0 jump (which does not exist). L-3 assumed the home was (B); it is (A).

> **Label.** (A) `[T]` dimensional; (B)-absent `[T-cond.]` (single-phase; transition-freeness of all of `[1/48,2.5]` is *unproved* — `PROOF-NO-MELT-PHASE`, `PROOF-TRANSPORT-CONTINUITY §5` — but even the hypothetical transition is handled in §3/§5.4 and lands on transport). The *disambiguation* itself is `[T]`: the two non-analyticities are at different points (`g=0` vs `β_c`) and carry different data (scale vs jump).

---

## §2 Lemma A — β-criticality is orthogonal to the **binding** (Schur / direction) wall  **[T]**

**Why horn β was never the load-bearing horn.** The rigidity no-go (`B3`) says: a datum *constant on continuous deformations* is constant on the connected β-family, hence `= ` its free value `≤0`. A β-critical functional is **not** deformation-invariant, so it **evades** horn β — but that only removes a constraint on a *particular proposed mechanism* (topological indices). The **binding** obstruction to a **lower** bound is the **direction/Schur** wall (moment lemma `GAP-G2-edge`; `SCHUR-CONVEX`; `EDGE-f1sq`), and β-criticality does not touch it:

> **Lemma A (analyticity-blindness of the direction wall). [T]** Let `Φ(ν)=F(∫φ₁dν,…,∫φ_kdν)` be any functional of **finitely many linear statistics** of the spectral measure (moments, generalised moments, susceptibilities, free energy, any `⟨q⟩`). Then `Φ` cannot lower-bound `inf supp ν`, **regardless of its dependence on β** — smooth, discontinuous, or wildly non-analytic. *Proof.* The linear-data lemma / `δ₀+δ₂` vs `2δ₁` witness (`E2-FLIP-NOGO §2-α`; `[V-here] W1`): measures with identical `(M₀,M₁)` (and, by the LP witness, identical `M₀..M₃`) and floors `0` vs `1`. `Φ` takes the **same** value on both; it cannot be `≤` both floors and positive. The argument uses only the values `∫φ_i dν`, not the β-smoothness of `F`. ∎
>
> Equivalently: an aggregate is **Schur-concave** (smallest on concentrated measures), so it is an **upper** bound on the edge; a **lower** bound needs a **Schur-convex coercivity** — a *total bounded above by the atom* (`EDGE-f1sq §1`, `PROOF-TRANSPORT-MONOTONE`) — and **no aggregate, analytic or not, supplies one.**

> **Remark (functionals of the *whole* β-profile, not finitely many moments).** The same conclusion holds with the finite list replaced by the entire profile `β↦(M_k(β))_k`. The banked **linear-data / Jaffe** witness (`E2-FLIP-NOGO §2-α`; `LEVER §0B`; Jaffe 1967) gives, at *each* β, `ν_ε = wδ_ε + (1−w)δ_S` whose **every** moment is finite and equals the gapped reference to `O(ε)` while the edge `=ε→0`. Hence the whole moment-profile can carry a sharp **β-critical onset** (β where `M_k` crosses over) *and* an edge that vanishes — exactly the decoupling `[V-here]` W2 exhibits at the level of `m₁(β)` itself. A functional of the profile still cannot lower-bound the floor, and its β-critical onset is calibrated to `Λ`, not to the gap (§4).

**Consequence.** L-3's claim "β-critical ⇒ escapes the direction wall" is **false**. The escape it actually achieves (evading horn β) was already achieved by *any* coupling-dependent functional, including the one (`Φ_Sav`, the Savvidy Hessian) that `E2-FLIP-DYNAMICAL` proved **fails direction + uniformity**. β-criticality is a property of the **wrong horn**.

---

## §3 Lemma B — the far-side-knowledge lemma: a discontinuous minorant must know the far side  **[T]** (logic) + `[V-here]`

The dream of L-3 is a functional that is `0` on the strong side and **jumps up** to `Δ>0` on the physical side, landing below `m₁` there "for free," using the seam rather than transport. This is impossible as a matter of the logic of the inequality:

> **Lemma B (far-side knowledge). [T]** Let `Φ` be any functional with `Φ(μ_β) ≤ m₁(μ_β)` for all `β` in the physical regime `P`, and suppose `Φ` is **discontinuous** at an interior `β_c` (so it "jumps"). Then the **largest value `Φ` may take on `P` while still satisfying `Φ≤m₁` there is `inf_{β∈P} m₁`**: to certify `Φ ≥ Δ` on `P` one must have `Δ ≤ inf_P m₁`, i.e. a positive lower bound on the gap over **all** of `P`. The "jump value" is admissible **iff** it is `≤` the far-side infimum; the seam calibration **is** the far-side gap.
>
> *Proof.* `Φ≤m₁` pointwise on `P` plus `Φ≥Δ` on `P` ⇒ `Δ ≤ inf_P m₁`. The discontinuity of `Φ` at `β_c` changes the *shape* of `Φ`, not the pointwise constraint it must satisfy on `P`. Whether `Φ`'s positivity on `P` "comes from" a jump or from anything else, the statement `Δ ≤ m₁` **on all of `P`** is the gap on `P`. ∎

This is tautological *as logic* and that is precisely its force: **the seam gives no calibration on either side.** To use a jump you must already place its upper plateau `≤ inf_P m₁` — you must know the gap on the physical side. `[V-here] W3` makes it numerical: sweeping the admissible constant `Δ` of `Φ=Δ·𝟙[β>β_c]` over a physical-side edge profile `m₁(β)↓g_phys`, the largest valid `Δ` is **exactly `inf_{β>β_c} m₁ = g_phys`** (to machine precision) — the far-side gap, for every `g_phys∈{1,0.5,0.1,…,0}`.

> **Why this is not just "Theorem 0 rephrased."** Theorem 0 says the *value* on the free endpoint is forced `≤0`. Lemma B is about an **interior seam**: the value on the **whole physical plateau** is forced `≤ inf_{P}m₁`, i.e. the seam-functional is *calibrated to the extremal gap on the far side*. The rigidity escape (jump at `β_c`) and the acid test (free endpoint) act at **different points**; Lemma B shows the interior jump cannot bypass the far-side extremum, independent of where `β_c` sits. This is the L-3-specific content.

---

## §4 Lemma C — the scale-vs-ratio wall: β-critical data give `Λ`, not `m₁`  **[T]** dimensional + **[J]** single-scale

Even granting a *Sharp* β-critical datum that is **not** an aggregate (so Lemma A does not apply) — e.g. the true non-analyticity (A) at `g=0`, or a genuine critical point — it supplies a **dimensionful scale**, never the gap's pure ratio:

> **Lemma C. [T]+[J]** Pure 4D YM is classically scale-free; on an `(L,a)`-uniform gauge-invariant local functional the only scales are `{Λ_QCD,√σ,g|B|,1/a,1/L}` (`E2-FLIP-NOGO §3, Theorem S`). Uniformity kills `1/a,1/L`; `√σ` is 1-form/cross-channel (`BS-THRESHOLD §4`); `g|B|` is UV/static-curvature (`TEST-C`). Leaves `Λ_QCD` — **and `Λ_QCD` is dimensional transmutation, i.e. exactly the non-analytic datum (A).** So any β-critical dimensionful `Φ = Ψ(Λ_QCD)` satisfies `(LB)` only if `Ψ(Λ) ≤ c·Λ` with `c=m₁/Λ` — **requiring the pure ratio `c`** (`[J]` single-scale: `c` β-independent). A dimensionless β-critical datum (a shape, an isolation margin) bounds `f₁²`, `Λ−m₁²=(ρ²−1)m₁²`, or `ρ=m_{0⁺⁺*}/m₁`, each `= ` the gap up to a universal factor (`L2` + SSC; `G2D-ISOLATION-I2`; `RATIO-RHO`). **Either way the β-critical content is `Λ` or a gap-ratio — the dimensionless nonperturbative number, not the gap's floor.**

> **One line.** *A β-critical functional is a scale-detector; the theory's only scale is the gap in disguise; and the gap is `c·Λ` with `c` — the pure number — the genuinely nonperturbative content. The non-analyticity pins the exponent (⟹ the scale); it leaves the ratio free.*

This is exactly `MS-MONOTONE-INVARIANT`'s "invariant = covariance `Z_O²` (`0↦0`)" lifted to the β-critical setting: RG-invariance / a transmutation datum **transports a known residue across scales; it never lifts an unknown one off zero.**

---

## §5 The construction attempts — each run, in its strongest form, to its exact death

I tried every named shape from the task. None survives; each dies on a *distinct, named* wall — not one repeated excuse.

### 5.1 Savvidy tachyon-condensation **onset**  → direction + wrong-sign  `[NG]`
The natural "β-critical onset": at the Savvidy crossover the perturbative vacuum is unstable (`λ_min=−g|B|`), so a condensate "turns on." **This is `E2-FLIP-DYNAMICAL`, fully worked:** the Savvidy Hessian depth is **anti-correlated** with the transfer gap (solvable witness: depth `−g→−∞` while the true gap `→0` as `e^{−S_inst}`), is extensive/UV (`∝V,∝β`, fails uniformity), and the direction-correct version (`Φ=e^{−S_inst}`, the tunnelling action) **is the transport machine `=H1`.** The onset is the wrong object (§2 direction) and the wrong direction (§4). **Dead.**

### 5.2 Condensation onset / β-non-analytic **order parameter**  → Schur-concave ⇒ UPPER  `[NG]`
Any order parameter is an expectation `⟨q⟩=∫q\,dμ` — an **aggregate** (Lemma A). `MS-MONOTONE-INVARIANT` `[NG]`: the gluon condensate `⟨¼trF²⟩`, `χ_t`, the anomaly coefficient are **totals**, the moment lemma turns them into **upper** bounds on the edge; `χ_t` even runs the wrong way (`m₁²≤χ_t/f₁²`, blows up as `f₁²→0`) and is the `P`-odd `0⁻⁺` channel. The "same marginal operator" saves nothing: `f₁²=Z_O²f₁²` is a **covariance** (`0↦0`) transporting a *known* residue. **Dead** (direction, before any β question).

### 5.3 **Anomaly-jump**  → free-field-false + protected-total ⇒ UPPER + does not jump  `[NG]`
Three independent kills (`MS-MONOTONE-INVARIANT §4`): (i) the conformal `a`/`c` charges are **nonzero for the free gauge multiplet** while `f₁²=0` ⇒ `f₁²≥c·a` is **false on the free field** (fails acid test in the *wrong direction*); (ii) `a` is a **protected monotone total** (`a`-theorem `[R]`), so it does not *jump* (it flows monotonically) and gives an upper bound; (iii) the `β²` log coefficient of `⟨θθ⟩` **vanishes as `g→0`** ⇒ acid-vacuous and a UV/high moment. **Dead** on three counts.

### 5.4 **Critical-point quantity** — finite-T deconfinement  → wrong channel + wrong direction  `[NG]`
The deconfinement transition (`Svetitsky–Yaffe`) is a **genuine** β-critical non-analyticity — and the *only* real one in the problem. But: (i) it lives at `T=1/(N_t a)>0`, gone at `T=0` (`N_t→∞`, `PROOF-NO-MELT-PHASE`); (ii) it is the **center `Z_N` / 1-form** transition — its order parameter is the nonlocal Polyakov loop, **orthogonal** to the `0⁺⁺` residue (`PROOF-1FORM-ANOMALY`: `trF²` is N-ality-0, blind to 1-form center); (iii) the **spatial screening mass** `m_scr(T) ≥ m_scr(0)=m_gap` (`SMALLCIRCLE-ADIA §3`, `PIROGOV-SINAI`) is an **upper** bound on the T=0 gap, and `T_c` itself is a fixed multiple of the single scale `Λ` (SSC ⟹ needs the ratio, §4). **Dead** (wrong channel, wrong direction, scale-vs-ratio).

### 5.5 **Dissolution-onset / isolation-endpoint `β_d`**  → circular = `H*`/`H1`  `[NG]`
The sharpest β-critical shape: let `β_d = sup{β : Δ(β)=Λ−m₁² ≥ κ}` be the **endpoint of the isolation region** — manifestly β-critical (a boundary). Claim: for `β<β_d`, `m₁ ≥ √(κ/(ρ²−1))` — a lower bound "from isolation." **Dies at definition:** `β_d` is defined via `Δ=Λ−m₁²`, a **gap object**; computing `β_d` *is* computing the gap's isolation margin. Uniformizing the margin is verbatim `ρ≥1+c = H*` (`RATIO-RHO §4`); and the dissolution that sets `β_d` is exactly the closure failure of `PROOF-TRANSPORT-CONTINUITY` (`Λ−m₁²→0`, `m_{0⁺⁺}→0`). **Dead** (circular = transport).

### 5.6 Strong-coupling **analytic continuation** to `β_c`  → value ≠ signed bound; crosses the singularity  `[NG]`
The strongest honest version: `m₁` is real-analytic (`PROOF-NO-MELT-PHASE`), so continue the strong-coupling series (radius `β_*≈O(1/N)≈1/48`) to `β=2.5`. Dies twice: (i) analytic continuation carries the **value**, not a **sign-controlled bound** (`MONOTONE-COUPLING §2`, `PROOF-NO-MELT-PHASE`); we know only `m₁≥Δ₀` on `[0,1/48]` — an inequality, not the function, so infinitely many analytic continuations are consistent. (ii) `2.5` is `~10²×` beyond the convergent disk; the continuation crosses **unknown singularity structure** — which *is* the crossover (`PROOF-TRANSPORT-CONTINUITY`). **Dead** (= transport).

**Net:** six shapes, six distinct deaths; **none** is "we didn't try hard enough" — three die on direction (Schur/aggregate, before any β question), one on scale-vs-ratio, one on wrong channel, one on circularity. β-criticality never moves any of them across the binding wall.

---

## §6 The dichotomy theorem for `L-3`, assembled  **[NG]**

> **Theorem (L-3 dichotomy).** For every `Φ ∈ L-3` (§0) satisfying (LB)+(acid), at least one holds:
> **(a)** `Φ` factors through **aggregates** of the spectral measure ⇒ Schur-concave ⇒ **upper** bound, not a lower bound — *unless acid-vacuous* (§2, `[T]`);
> **(b)** `Φ` is **dimensionful** and β-critical ⇒ lands on `Λ_QCD` (the transmutation datum) ⇒ `(LB)` needs the pure ratio `c=m₁/Λ` — the gap (§4, `[T]`+`[J]`);
> **(c)** `Φ` is a **dissolution / isolation / atom** object ⇒ it **is** the gap's own β-critical singularity (`Λ−m₁²`, `f₁²`, `ρ`), and to be a uniform lower bound it must be `β`-monotone/protected — which **no** atom-coercive object is (`PROOF-TRANSPORT-MONOTONE`, opposite sides of the Schur order) ⇒ needs the transport (§4, §5.5, `[NG]`);
> **(d)** `Φ` uses **field-space** structure beyond the spectral measure ⇒ that is `L-2`, and `L-2 ≡ H1` (`E2-FLIP-NOGO §5-L2`).
> **There is no fifth, non-circular, uniform lower-bound route.** In particular, "β-critical" (non-analyticity) only affects the **rigidity** horn — which is not binding — and the **far-side lemma** (§3) shows a discontinuous minorant must be calibrated to the far-side gap. ∎
>
> **Label: `[NG]` for the class `L-3`, uniform statement.** Assembled from `[T]` (Lemma A moment lemma; Lemma B logic; Lemma C dimensional analysis + `PROOF-NO-MELT-PHASE` for the T=0 absence; `PROOF-TRANSPORT-MONOTONE` atom-coercive ⟂ β-monotone) and `[NG]` per-species (§5). **Not a global barrier** (`∀class ≠ ∀possible`).

**Exact stuck step (L-3, one box):**
```
need:   a beta-critical minorant  Phi <= m1  (uniform, acid-clean) computable WITHOUT transport.
have:   the only beta-critical content of m1 is  (A) Lambda  at g=0  [the SCALE, not the ratio c],
        or (B) a T=0 beta_c JUMP  [DOES NOT EXIST: PROOF-NO-MELT-PHASE; the real beta_c is finite-T,
        wrong channel], or the dissolution Lambda-m1^2->0  [= the gap's own singularity].
gap:    to get  m1 >= Delta  one needs  c = m1/Lambda  (or  rho=1+c),  a PURE nonperturbative
        number; and a discontinuous minorant must be <= inf_phys m1 (far-side lemma).
        controlling c / the isolation / the far side  ACROSS beta_c  IS the nonperturbative
        transport = confinement.                                            ==> [O]
```
Equivalently, the missing object is `PROOF-TRANSPORT-MONOTONE`'s **"β-monotone **and** atom-coercive functional"** = uniform isolation `Δ(β)≥δ` on `[1/48,2.5]` = `H*(r₀)` — one rung above the gap — under its β-critical name.

---

## §7 `[V-here]` computational witness  (`research/unconditional/e2_flip_l3_witness.py`)

Pure `numpy`, single process, ~0.3 s (`py -3.12`). It tests the **logic** of L-3 exactly (no claim about the real SU(N) spectrum):

- **W1 (moment lemma).** `δ₀+δ₂` and `2δ₁` have identical `(M₀,M₁)=(2,2)` and floors `0` vs `1`; any `Φ(M₀,M₁)` takes the same value on both ⇒ cannot lower-bound the floor. Direction check: `M₁/M₀ ≥` floor, i.e. moments are an **upper** bound. (Order-3 upgrade = banked linear-data lemma / `E2-FLIP-NOGO` W1.)
- **W2 (β-criticality ⟂ gap).** A crossover edge with sharpness `κ`: `max|dm₁/dβ|` grows without bound (`0.25→2.5→25→250→2500` for `κ=1,10,10²,10³,10⁴`) — approaching a true **jump** (fully "β-critical") — while the physical-side gap `inf_{β>β_c}m₁=g_phys` is set **independently** (`g_phys∈{1,0.1,10⁻³,0}`). **Sharpness carries no lower-bound information about `m₁`.**
- **W3 (far-side lemma).** The largest admissible constant of `Φ=Δ·𝟙[β>β_c]` is **exactly** `inf_{β>β_c}m₁ = g_phys` (machine precision, every `g_phys`). A discontinuous minorant is calibrated to the **far-side gap**.

**Reading.** W1 kills the aggregate class at direction; W2 shows β-criticality is decoupled from the gap (no free lunch from sharpness); W3 is Lemma B numerically. These **witness already-banked lemmas** (`GAP-G2-edge` moment lemma; `SCHUR-CONVEX`; Theorem 0) plus the L-3-specific **decoupling W2**; they are not new evidence about the world.

---

## §8 Honest boundary, flip condition, status

- **Does NOT prove Clay; does NOT prove no tool can exist** (`CONVERGENCE-redteam §d`, `KEY-EXISTS §5`). It is `[NG]` for the named class `L-3` (β-critical / non-analytic functionals), assembled from `[T]` lemmas + per-species `[NG]`, not a global barrier.
- **(B)-absent is `[T-cond.]`** (single-phase): transition-freeness of all of `[1/48,2.5]` is unproved. But the hypothetical T=0 transition is covered: even granting a genuine `β_c` jump, the **far-side lemma (§3)** forces the jump amplitude `≤ inf_{phys}m₁` (the gap), so the escape still requires transport. The **scale-vs-ratio** and **direction** halves are `[T]` independent of the phase question.
- **Do not under-sell what is real:** the gap *is* β-non-analytic — **at `g=0`** (transmutation, `[T]`). That is a genuine, central, nonperturbative fact; it is simply the **scale `Λ`**, and the gap is `c·Λ`, not `Λ`. The non-analyticity is real; it is the **ratio** `c` it does not give.
- **Flip condition (exactly what would make L-3 FLIP).** Exhibit `Φ` that is (i) a uniform `(L,a)` lower bound on `m₁`, (ii) acid-clean, (iii) **β-critical** yet **not** an aggregate (Lemma A), **not** dimensionful-onto-`Λ` (Lemma C), **not** the gap's own isolation/atom (case c), **not** field-space `=H1` (case d), and (iv) whose far-side plateau is certified **without** the far-side gap (Lemma B). By §1 such a `Φ` must lock onto a **T=0 β_c jump that does not exist** (or a finite-T datum in the wrong channel), **or** carry the pure ratio `c` nonperturbatively — i.e. it must **be** the transport machine. **None is known; finding one is the task.** (Per `E2-FLIP-NOGO §6`, a genuine 6th species would live in `L-2`/`H1`, *not* here.)
- **Relation to the corpus.** This note does not overturn `E2-FLIP-NOGO §5-L3`; it **sharpens** it from "a functional that jumps must detect confinement" to a five-part dichotomy with the **`g=0`-vs-`β_c` disambiguation** (new), the **far-side lemma** (new, `[T]`-logic + `[V-here]`), and the **analyticity-blindness of the direction wall** (new framing), and it corrects the loose "the gap is a β_c jump" reading using `PROOF-NO-MELT-PHASE`.

---

## §9 References

**Internal (absolute paths; the banked verdicts this note assembles and sharpens):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall / missing machine); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five conditions; §4 transport≡confinement); `D:\Dev\math\docs\new-tools-registry.md` (§A1–A5 tools, §B1–B6 no-gos incl. **B3 rigidity**, §D residual dictionary D.1–D.7).
- `D:\Dev\math\research\unconditional\E2-FLIP-NOGO.md` (**§5-L3 the target**; §2-α moment/linear-data lemma; §3 Theorem S scale dichotomy; §5-L2; §6 flip condition); `E2-FLIP-DYNAMICAL.md` (Savvidy lever anti-correlated; `e^{−S_inst}`=H1); `E2-FLIP-KINEMATIC.md` (threshold anchor; L1/SSC); `E2-REGRESS-CIRCULARITY.md` (§0.4 SSC; §2 vNW; §6 flip condition).
- **`PROOF-NO-MELT-PHASE.md`** (T=0 real-analyticity; `β_c≈2.47` = analytic crossover / finite-`N_t` deconfinement; **the disambiguation §1**); `PROOF-TRANSPORT-CONTINUITY.md` (clopen; closure fails by dissolution `Λ−m₁²→0`); `PROOF-TRANSPORT-MONOTONE.md` (**atom-coercive ⟂ β-monotone**; totals→upper); `MONOTONE-COUPLING.md` (direction wall: proven end is more-gapped; analytic continuation ≠ signed bound); `MS-MONOTONE-INVARIANT.md` (RG-invariant = covariance `Z_O²` `0↦0`; aggregates→upper); `MS-GAP-EQUATION.md` (transmutation `m²∼e^{−c/g²}`; DSE branch collapses to gapless at 2.5); `RATIO-RHO.md` (`ρ`=`H*`; LLR is free/vacuous).
- `INVENT-NEW-INDEX.md` (**B3 rigidity**; "gap must be β-non-analytic-jump" as an index statement); `G2D-ISOLATION-I2.md` (§2 SSC; §3 vNW pointwise; §5 `ρ`=`H*`); `GAP-G2-edge.md` (moment lemma); `SCHUR-CONVEX.md`; `EDGE-f1sq` / `PROOF-TRANSPORT-MONOTONE` (Schur-convex coercivity absent); `SUPPORT-NOT-ATOM.md`; `BS-THRESHOLD.md` (`√σ` 1-form); `TEST-C-GLOBAL-SCHUR.md` (`~0.6β` ≠ transfer gap); `SMALLCIRCLE-ADIA.md` (§3 screening `≥m_gap`; backward inequality); `W3-adjacent-4D.md` (scale-free ⇒ `Λ`); `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_min=−g|B|`; `N_neg∝V`); `PROOF-GAUSS-INVARIANT-PF.md`, `PROOF-P3-UNCONDITIONAL.md` (`β<1/48` gap `[T]`); `PIROGOV-SINAI.md` (screening = susceptibility).
- **`D:\Dev\math\research\unconditional\e2_flip_l3_witness.py`** (`[V-here]` W1–W3; run `py -3.12 e2_flip_l3_witness.py`).

**External — `[R]` standard (not re-verified this session):** Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737; Kato, *Perturbation Theory for Linear Operators* (2nd ed., Ch. VII); Svetitsky–Yaffe, *Nucl. Phys. B* **210** (1982) 423; Gross–Pisarski–Yaffe, *RMP* **53** (1981) 43 (running coupling / finite-T); Gross–Wilczek, *PRL* **30** (1973) 1343 & Politzer, *PRL* **30** (1973) 1346 (asymptotic freedom / transmutation); Savvidy, *PLB* **71** (1977) 133; Komargodski–Schwimmer, *JHEP* **07** (2011) 018 (`a`-theorem). No fabricated bound; no QED.

---

*Status: **=TRANSPORT≡CONFINEMENT — `[NG]` for the class `L-3` (β-critical / non-analytic functionals), uniform statement.** No (a) constructed; (b) proved. Sharpened by the **`g=0`-vs-`β_c` disambiguation** (the true non-analyticity is transmutation at `g=0` = the *scale* `Λ`, not a T=0 `β_c` jump — `PROOF-NO-MELT-PHASE`), the **analyticity-blindness of the direction/Schur wall** (`[T]` moment lemma — β-criticality evades only the non-binding rigidity horn), the **far-side-knowledge lemma** (`[T]`-logic + `[V-here]` W3 — a discontinuous minorant is calibrated to the far-side gap), and the **scale-vs-ratio wall** (`[T]`+`[J]` — β-critical data give `Λ`, the lower bound needs the pure ratio `c=m₁/Λ`). Exact stuck step: the only β-critical Schur-convex minorant is the gap's own isolation/atom singularity, and no such object is simultaneously β-monotone/protected (`PROOF-TRANSPORT-MONOTONE`); controlling `c`/isolation/far-side across `β_c` **is** the nonperturbative transport machine. NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). **Clay stays open.***

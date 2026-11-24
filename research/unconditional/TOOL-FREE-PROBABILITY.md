# TOOL-FREE-PROBABILITY — does the N=∞ **spectral curve** (Makeenko–Migdal / Douglas–Kazakov / GWW / Ginsparg–Zinn-Justin) give a NON-CIRCULAR **support-gap = mass-gap** handle for the 4D `0⁺⁺` sector?

> 2026-11-24 · Scouting report on the **free-probability / RMT / large-N master-field spectral-curve** species for the 4D pure SU(N) mass gap (Clay G2). Body of the question: at `N=∞` the master field is governed by the Makeenko–Migdal loop equation, whose *reduced* solutions are spectral curves with a support (cuts) and possibly gaps; does a **gap in the spectral curve** (a support gap, the needed *lower*-bound shape for `spec(H)⊆{0}∪[Δ,∞)`) survive as a rigorous handle, or is it planar-only / strong-coupling / circular?
> **Labels (charter §2): [T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断. **[V-here]** = fresh reasoning/witness in this session (no repo state mutated, no heavy job run). External facts are **[R]/[standard]** unless marked **[V]**; I re-ran no external source. **No QED is claimed on anything unproven, and this document does NOT prove Clay, nor that no tool can exist.**
> **Prior art in-repo (must not be re-listed as new):** `E3-LOOPSPACE` (loop-space generator ⇒ Giles ⇒ field generator / Savvidy `Γ₂`), `DYSON-GAS` (transfer not a Vandermonde log-gas; 2D GWW/DK edge = free energy), `E4-RMT-EDGE` (RMT edge = random-Schrödinger tachyon edge, anti-gap), `k0-FREEPROB-AVERAGE` (free CLT pins the `1/V` volume average, not the `k→0` gap), `H1-T8-free-prob` (freeness is the right *language* for N=∞, not a gap route), `PROOF-LARGE-N-RIGOROUS` (rigorous large-N = Chatterjee strong-coupling master field; no continuum gap), `ABC-DUALITY` (MM/1-form exact structures = 跨道/错向/错区). This report is additive only where marked **[V-here]** and in the *classification* (§0, §2, §4).

---

## Verdict up front

> ### **[NG] for the spectral-curve class as a non-circular `0⁺⁺`-gap handle + [O] one exact missing step. NO genuine leverage. This species dies of a wall that is *mostly not* the campaign's "wrong-direction" one — that is the refinement this report contributes.**

**Four-line answer to the four posed questions.**
1. **Does the N=∞ spectral curve have a `0⁺⁺` gap, rigorously at `N=∞`?** — *The premise is ill-posed in 4D.* In genuine 4D the loop equation is a **nonlinear functional PDE on loop space with no finite Riemann-surface solution**; a "spectral curve" exists only after a *finite-dimensional reduction* (0D semicircle; 2D GWW/DK; Eguchi–Kawai). Those reductions carry **no `J^{PC}`/`0⁺⁺` label at all** — the curve is the eigenvalue density of a *holonomy/single-matrix*, all channels mixed. So "a gap in the `0⁺⁺` sector of the spectral curve" is a **category error**, not a hard theorem. Where the curve *is* rigorous (2D), the "gap" is a large-N **free-energy / eigenvalue-support transition** (Douglas–Kazakov, third order), i.e. H1, **not** the mass gap (§1, §2).
2. **Does it survive `N→∞→3`?** — No theorem. The curve *is* the N=∞ object; at N=3 there is no curve (finite matrix, discrete spectrum), and the only bridge is the `1/N²` genus series, **asymptotic with no N=3 bound** (`PROOF-LARGE-N-RIGOROUS §c`). A uniform N=3 shadow requires the support gap to exceed the `O(1/N²)` corrections — which *assumes* the gap is uniform (circular). §5.
3. **Free probability / free convolution / R-transform — a gap criterion?** — A **rigorous gap *calculus* exists** (Voiculescu R-transform linearizes ⊞; Belinschi–Bercovici support/atom theory; Biane free-Brownian-motion support evolution). But it is a **calculus that needs its input measures**, and the inputs ("single-link / master-field measures") are the *unbuilt 4D master field*. Moreover it is a *calculus of the auxiliary eigenvalue measure*, not of the physical Källén–Lehmann (KL) measure. §3.
4. **Genuine handle, or planar-only / strong-coupling / circular?** — **All three dials read *bad*, but in a specific order:** **(a) planar-only** — yes, the curve is by construction the `N=∞` limit and there is no controlled `1/N` backterm; **(b) strong-coupling/2D** — yes, the *rigorous* curves live exactly at fixed small `β` (Chatterjee strong coupling), in 2D, or at 0D; the 4D physical point `β→∞, a→0` is unbuilt; **(c) circular** — yes, and this is the *decisive* dial: the **support edge is pinned by the equilibrium-measure potential, and that potential IS the effective action / Hamiltonian** (`DYSON-GAS §3`). §4.

**What is genuinely *not* the old objection (the refinement).** The campaign note "rigorous large-N gives `S_N ≥ S_∞` (wrong direction)" and "large-N only strong-coupling / wrong side" are **free-energy/value** objections. A **support-gap** statement is a *topological/support* statement (is an interval empty?), **not a value inequality** — so in *form* it is **not** the `S_N ≥ S_∞` wrong-direction species, and (crucially) the curve-gap **detector passes the free-field acid test non-vacuously**: the free Gaussian matrix model is **one-cut** (no gap, correct: free massless is gapless), whereas a confining two-cut potential *would* show a gap. **So the acid test is not what kills it.** What kills it is: **(i) measure mismatch** (`supp ρ_curve` = auxiliary eigenvalue density ≠ `supp ν` = physical KL support; the bridge is OS + transfer, and the transfer's gap is a *different functional*), **(ii) circularity** (support pinned by effective action = H1), **(iii) regime** (rigorous only strong-coupling/2D/0D), **(iv) transfer** (no `N=3`). This is the honest classification: **not wrong-direction; it is H1 disguised as a Riemann surface.** The rest of this report substantiates each.

---

## §0 The object, stated precisely (so the verdict is checkable)

**The physical target** (`BRIEF §1`, `SUPPORT-NOT-ATOM`): gap ⟺ support floor `supp ν ⊂ [m₁²,∞)`, `m₁>0`, where `ν ≥ 0` is the **Källén–Lehmann spectral measure of the physical `0⁺⁺` two-point function** `⟨0|O(t)O(0)|0⟩`, `O=¼tr F²`. Clay G2 needs only the **support** statement `spec(H)⊆{0}∪[Δ,∞)` (it does *not* need the atom `f₁²>0`; `E2-REGRESS §0.1`). So a *support gap* is **direction-correct in shape** — this is the one thing the spectral-curve species has going for it.

**The spectral curve** (`Makeenko–Migdal PLB 88 (1979) 135`; `Douglas–Kazakov PLB 319 (1993) 219`; `Gross–Witten PRD 21 (1980) 446` / `Wadia PLB 93 (1980) 403`): in a *reduced* large-N model the loop-equation/partition function becomes a matrix integral `Z=∫dM e^{−N tr V(M)}`. At `N=∞` the eigenvalue density `ρ` is the **equilibrium measure** of a Coulomb-gas variational problem,
```
ρ = argmin_{ρ≥0, ∫ρ=1} [ ∬ log|x−y|^{-1} dρ(x)dρ(y) + ∫ V(x) dρ(x) ],
```
whose **support** is a finite union of cuts `∪[a_i,b_i]`; the "spectral curve" is the (hyper)elliptic curve carrying the resolvent `ω(z)=∫dρ/(z−x)`. The **gaps** are the intervals `(b_i,a_{i+1})` not in the support. Classical cases: 0D ⇒ Wigner semicircle (single cut, **no gap**); 2D YM ⇒ Gross–Witten–Wadia unitary matrix model (single **arc** on the circle; a **gap/non-uniform arc** opens at the GWW critical coupling); the `(2,2m−1)`-type multicritical double-scaling limits ⇒ `Ginsparg–Zinn-Justin` 1990–93 (**[R]**, exact volume not re-verified) — all of these are **equilibrium-measure cuts of an auxiliary matrix**, and "criticality"/"gap-opening" is exactly where the **free energy is non-analytic** (third order for GWW/DK).

**Two structural facts carried from the banked corpus and re-emphasized [V-here in placement, not content]:**
- **F1 (no 4D curve).** Makeenko–Migdal is exact as a tower at finite N (linear, never closes) and closes at `N=∞` only by **factorization** into a *quadratic functional* equation on loop space (`PROOF-LARGE-N-RIGOROUS §b`). It is **not** a finite Riemann–Hilbert problem in 4D; no spectral curve exists except in the solvable *reductions*. (The word "spectral curve" is imported from the reductions.)
- **F2 (measure mismatch).** `ρ_curve` is the **equal-time / static eigenvalue density** of a holonomy or reduced matrix. The mass is the **Euclidean-time decay rate** (atom of `ν`), read through the transfer/KMS generator. `ρ_curve` and `ν` are **different measures**; the only bridge is OS reconstruction + the transfer operator, whose gap is the `k→0` generator's gap — a **different functional** of `ρ_curve` than its support (`H1-T8 §2`, `k0-FREEPROB §3`, `E3-LOOPSPACE §2` via Giles: `spec(L_loop)=spec(L_field|_{gauge-inv})`).

---

## §1 Does the N=∞ spectral curve have a `0⁺⁺` gap, rigorously?

**NG as stated (category error); [O] for the only well-posed version.** By **F1**, in 4D there is no curve to inspect. By construction, the curves that *do* exist (0D/2D/EK) are eigenvalue densities of a matrix with **no `J^{PC}` decomposition**: the GWW curve is a density on the eigenvalue circle `θ∈(−π,π]` of a *unitary holonomy*, and "gap" = an **arc with no eigenvalues** — a statement about **holonomy eigenvalue clustering**, which is the **center/1-form deconfinement order parameter**, not the `0⁺⁺` glueball. There is no map from "arc gap on `U(1)`/`SU(N)` holonomy" to "`0⁺⁺` support floor" — the former is 1-form/center data, the latter 0-form/`tr F²` data, and the cross-sector barrier is banked (`BS-THRESHOLD §4`, `E2-REGRESS §4 item 2`, `ABC-DUALITY §2`, `PROOF-1FORM-ANOMALY`).

**The one well-posed reading:** "in a *reduced model whose eigen-observable is `O=¼tr F²`*, does the equilibrium measure have a cut gap that lower-bounds the physical `0⁺⁺` edge?" — This is **[O]**: it presupposes (i) a 4D physical-coupling master field (F1, unbuilt) and (ii) the identification `ρ_curve = ν` (F2, unproven). No theorem supplies either. **Verdict §1: [NG] as posed (no `0⁺⁺`-labeled curve in 4D); [O] as the only sensible reading.**

---

## §2 Support-gap ≠ physical floor: the 2D GWW/DK witness [V-here: the witness framing]

This is the sharpest single point, and it is a **genuine acid test of the *curve-gap detector*** — the exact analogue of "must fail on the free field," applied to the one place the curve is exactly solvable.

**The witness.** 2D pure YM is **exactly solvable** at `N=∞` as the Gross–Witten–Wadia unitary matrix model; its eigenvalue support on the holonomy circle undergoes a **third-order (Douglas–Kazakov) large-N transition** at the critical area `βA = π` (conventions aside): one side the eigenvalues fill the circle (**no gap**), the other they collapse to an arc leaving a **gap**. Both the spectral curve and the gap are **rigorous at N=∞** ([T]-level for the model). **Yet this "gap" is NOT the 2D mass gap.** It is:
- a **free-energy singularity** (`−log Z` non-analytic, third order) ⇒ the **effective action** ⇒ H1; and
- the **center/1-form confinement-deconfinement order parameter** (eigenvalue clustering on the holonomy) ⇒ **cross-sector** to the local `0⁺⁺` `tr F²` channel.

So we have a **rigorous, exactly-computable spectral-curve support gap that provably does not correspond to a `0⁺⁺` mass gap.** (2D pure YM has no propagating glueball at all; the GWW gap tracks the *area-law/center* structure, not any `0⁺⁺` pole.) **A support gap in the spectral curve is therefore demonstrably insufficient** — the identification "curve gap ⇒ physical mass gap" fails in the *only* dimension where all terms are known exactly. This is the decisive structural objection, and it is stronger than "large-N is strong coupling": it is a **counterexample to the implication**, not a regime excuse.

**Morphology (subsidiary, [J]).** Even in a Hermitian one-matrix model the natural gap is a **symmetric two-cut** `[−2σ,−Δ]∪[Δ,2σ]` around the potential's double-well minimum: a gap **centered at 0** (or at the potential's symmetry point), with **soft (`~√(x−edge)`) edges**, not a floor `[Δ,∞)` above a vacuum at `0`. The physical spectrum has the **vacuum pinned at 0** by Perron–Frobenius (`PROOF-GAUSS-INVARIANT-PF`) and excitations above; the curve has **no vacuum**. So the *geometric meaning* of the curve's gap (empty interval between two symmetric bands) is not the physical gap (nothing between vacuum and the lightest excitation). **This matters only if one tried to lean on the curve's "hardness"; it does not rescue F2** — the measure mismatch already decides it.

**Summary §2:** the curve-gap detector is **acid-clean in form** (free Gaussian = one-cut = no gap; a confining potential *would* be two-cut), but a **support gap of the auxiliary measure is neither necessary nor sufficient for the physical `0⁺⁺` floor** — the 2D witness shows it is not sufficient, and F2 shows the auxiliary measure is the wrong object. **[NG].**

---

## §3 Free probability / free convolution / R-transform — the gap criterion that exists but does not bite

**What is rigorous (the calculus is real).** [standard/R] Voiculescu: free additive convolution `⊞` and the **R-transform linearizes it**, `R_{μ⊞ν}=R_μ+R_ν`; free CLT ⇒ semicircle; the N=∞ limit of *independent-in-color* matrix ensembles is **free**. **Belinschi–Bercovici**: a complete **support and atom theory** for free convolutions — `supp(μ⊞ν)` is read off the **subordination functions** (analytic, on the physical sheet), **gaps open/close at critical parameters of the semigroup** `μ^{⊞t}` (free Brownian motion, **Biane**), and **atoms of free convolutions are special** (they occur only under specific subordination identities). Free convolution with a semicircle is **absolutely continuous with a real-analytic density** on the interior of its support (**Bercovici–Voiculescu**) ⇒ **no atoms**, and gaps are *support* features of an `L¹` density, not spectral floors. **So yes — there is a rigorous gap criterion for free convolutions.** Sub-question (3)'s answer is *yes to the existence of a criterion*.

**Why it does not bite the mass gap (`H1-T8 §2-3`, `k0-FREEPROB`, re-derived [V-here] as a three-failure chain).**
1. **Wrong measure.** `⊞`/R describe the **eigenvalue measure of sums of freely independent *ensembles*** — i.e. the **static master-field eigenvalue density** (F2). The physical target is the **time-ordered KL measure** `ν` of `⟨O(t)O(0)⟩`. Freeness is a property of the **tracial algebra / Wilson-loop factorization**, *not* of `ν`. No theorem identifies `ν` with a free convolution.
2. **Unbuilt inputs.** A gap criterion for `μ₁⊞μ₂` needs `μ₁,μ₂`. The "elementary measures" here are the single-link / single-site master-field measures = the **unbuilt 4D master field**. A calculus is not an oracle: it needs its arguments.
3. **Atom-free output.** Even granting (1)+(2), the free-convolution support is that of an **absolutely continuous** density (Bercovici–Voiculescu smoothing) — a **soft-edged** gap, the *opposite* morphology to the mass-consistent KL measure whose bottom is an **atom/floor**. The `⊞` machine smooths atoms away; the mass lives exactly at an atom. (`DYSON-GAS §3`'s "soft Airy edge = no atom = anti-mass," here in free-convolution dress.)

**Verdict §3: [NG] as a `0⁺⁺`-gap handle (wrong measure + unbuilt inputs + atom-smoothing); [O] as "prove `ν` is a free convolution of computable elementary measures" — no such identification exists, and it is not the free-field-limited step, it is a structural mismatch.** The R-transform is a genuine *language* for the N=∞ limit (this is `H1-T8`'s honest credit), not a gap oracle.

---

## §4 The "wrong-direction (`S_N ≥ S_∞`) / strong-coupling" objections — does the support gap avoid them?

**This is where the report adds precision over the banked death-map line "large-N = strong-coupling / wrong-direction."**

- **It avoids the *letter* of `S_N ≥ S_∞`.** That objection is about **value/entropy inequalities** on the free energy (the rigorous large-N objects are typically one-sided in a way that bounds a gap from *above*, or a sector min–max `gap(sector) ≥ gap(full)`, `ABC-DUALITY §1(iv)`). A **support-gap** statement ("this interval is empty") is **not** a value inequality; in form it is the right *kind*. **Credit where due: the spectral-curve support gap is not killed by the wrong-direction argument per se.**
- **But it lands on the *equilibrium-measure variational problem*, whose potential IS the effective action.** `ρ_curve = argmin` with potential `V`; the **support edge `a₁` is the smallest `x` with `U(x)+V(x)=ℓ`** (Robin constant `ℓ`), i.e. the **free boundary is pinned by `V`**. In the physical reading `V = aH` (or `V` = the effective action coupling/area). Then `E₀ = inf supp ρ_curve = inf spec H = m₁` = **the gap**. `DYSON-GAS §3`: "the edge location = potential `W=aH` = `inf spec H` = the gap ⇒ circular." So the *support-gap* formulation is **H1 relocated onto a free-boundary condition**: lower-bounding the edge `⇐` lower-bounding `U+V` `⇐` lower-bounding the **action at the `0⁺⁺` configuration = the gap**. This is the *same* structure as `E2-REGRESS L1` (§0.2): a variational problem pins an edge only through a **hard spectral input**, and here the input is `V`, which is the answer.
- **Strong-coupling / wrong side — YES, the curves that are rigorous are exactly the wrong regime.** Chatterjee's master field is at **fixed small `β`** (strong coupling, large 't Hooft `λ`); GWW/DK are **2D**; Eguchi–Kawai is a quenched/twisted **single-site** reduction (center-unstable off the patch). The physical point `β→∞, a→0` (small `λ`) is **not covered** (`PROOF-LARGE-N-RIGOROUS §a`, `ABC-DUALITY §1`, `H1-T8 §1`).

**Verdict §4: the support gap *does* dodge the literal `S_N ≥ S_∞` value-direction objection — that is the genuine refinement — but it re-faces the wall as (i) free-boundary pinning by the effective action = H1, and (ii) the strong-coupling/2D regime. [NG]/[O], not a clean escape.**

---

## §5 Transfer `N=∞ → 3` [NG]

- The spectral curve **is** the `N=∞` object. At `N=3` there is **no curve**: a finite matrix has finitely many eigenvalues, a discrete spectrum, no equilibrium density. "The spectral curve at `N=3`" does not denote anything.
- The only bridge is the **`1/N` genus expansion** `m(N)=m_∞+c/N²+…`, which is **asymptotic**: for each `k` a remainder `O(N^{−k−1})`, **no bound at `N=3`** (`PROOF-LARGE-N-RIGOROUS §c`, `H1-T8 §4`, `C1`, `GAP-IN11`). The correction `c/9` is uncontrolled.
- Even granting a `1/N²`-controlled regime in principle: to transport an `N=∞` **support gap** to `N=3` you need the gap to **dominate the `O(1/N²)` eigenvalue fluctuations** near the edge. That is precisely the statement "the gap is uniform / bounded away from `0`" — **the answer, assumed**. A support gap whose width could close as `O(1/N²)` transports to **nothing**.

**Verdict §5: [NG] — no theorem, and the only bridge (genus series) is asymptotic and circular at the needed step.**

---

## §6 Five-point scorecard (`BRIEF §2`)

1. **Nonperturbative — PASS in principle, unbuilt at the physical point.** Freeness / free CLT / equilibrium measures are exact `N=∞` statements (not `1/N` series); *but* the 4D physical-coupling master field is the unbuilt G1 object, and the rigorous curves are 0D/2D/strong-coupling. (`H1-T8 §1`.)
2. **Sign-controlling LOWER bound — FAIL as a *gap* bound.** A **support gap** is the right *shape* for G2, and the detector is free-field-clean (§2) — **but** it bounds the support of the **wrong measure** (F2), and the support edge is pinned by the effective action (§4). It cannot produce `m₁²(k→0) ≥ λ` physically. **Pass in letter, fail in object.**
3. **Non-circular — FAIL.** Support edge = equilibrium free boundary = potential `V` = `aH` / effective action = `inf spec H` = the gap (`DYSON-GAS §3`); or, to go `ρ_curve → ν`, OS + transfer, whose gap *is* the answer (`H1-T8 §2`, `E3-LOOPSPACE §2`). Circular by two independent roads.
4. **Non-SUSY + fails free field — PASS (narrow, non-vacuous in form).** Non-SUSY ✓; free Gaussian matrix model is **one-cut, no gap** ⇒ correctly returns gapless ✓. This is the **genuine credit** (the detector *would* distinguish one-cut free from two-cut confining) — the failure is not the acid test but §§2–5. (Same credit pattern as `ABC-DUALITY`'s constant-mode/Haar branch and `DYSON-GAS`'s vacuous pass — but here it is *non*-vacuous because a two-cut potential is a genuine gap in the auxiliary measure.)
5. **Uniform `(L,a)` lower bound — FAIL.** Edge location is ensemble-set with `N^{−2/3}`/Tracy–Widom fluctuations and soft-edge density `ρ~√(x−a₁)`; `1/N²` transfer uncontrolled; rigorous curves only at fixed small `β`/2D/0D — no `λ` uniform over `β∈[1/48,2.5]`, `V→∞`, `a→0`.

**Score: #1 pass-in-principle; #4 pass (non-vacuous, the real credit); #2, #3, #5 fail, and #2/#3/#5 are ONE object — "the support edge is pinned by the effective action = H1."**

---

## §7 The exact missing step, precisely stated [O]

```
(O-FP)  ∃ a 4D physical-coupling (β→∞, λ(1/a)→0, a→0) large-N master field whose
        eigenvalue measure ρ on a gauge-invariant O=¼tr F² reduction
        (a) has a hard support gap  (b) that pins the 0⁺⁺ channel,
        AND a theorem
        (b')  hard support gap of ρ  ⇒  spec(H)|_{0⁺⁺} ⊆ {0} ∪ [Δ,∞),  Δ>0,
        proved WITHOUT using the transfer operator's own gap (the effective action),
        uniform in (L,a), and stable at N=3.

No such (O-FP) exists. (a),(b),(a') are the unbuilt 4D master field (F1) + the
unproven ρ=ν identification (F2) + the 0⁺⁺-sector map (cross-sector, §1);
(b') is the OS+transfer "trace→dynamical" bridge — whose content IS H1.
The 2D GWW/DK witness (§2) shows (a)⇒(b') is FALSE with all terms known:
a rigorous curve gap coexisting with no 0⁺⁺ mass gap.
```

This is the same wall in a third coordinate. **"A gap in the spectral curve = a support gap"** is **true only after** the circular identification `ρ_curve = ν` and the cross-sector map to `0⁺⁺`, and **only if** "support gap" is upgraded from the auxiliary measure's soft-edged cut to a physical floor — and the one exactly-solvable case shows the upgrade fails.

---

## §8 What would flip it (and why it is the same wall)

**Flip condition.** Exhibit **one** of:
(A) a 4D physical-coupling master field, constructed, whose `¼tr F²`-reduced equilibrium measure is **provably hard-gapped in the `0⁺⁺` sector** *and* whose gap is established **independently of the effective action**; or
(B) a theorem `hard support gap of ρ_curve ⇒ exp. clustering / KL floor`, **non-circularly**, i.e. without the transfer gap in the hypothesis.

**(A)** is the unbuilt weak-coupling 4D master field (the G1 object) **plus** the gap's independence from the action — but the action *is* the potential pinning the edge (§4), so (A) is self-defeating as an *independent* handle. **(B)** is the OS+transfer bridge = H1 in spectral-curve language. **Neither can be built without the missing machine**, exactly the charter's "the machine does not exist." There is a *third* conceivable lever — **the concentration/factorization of the master field as a rigidity input** — but `k0-FREEPROB`/`E4-RMT-EDGE §3` already place it: concentration pins the **aggregate** (`1/V`) and the **edge's fluctuation**, and the conditioning that would fix a sign/edge is a **large-deviation tilt whose rate is the effective action = H1**. So the apparent "third lever" is H1 again. **No flip is available; record the class as [NG] + [O].**

**One honest partial to keep.** The **support-gap detector is the cleanest acid-test-clean *shape* in the whole free-probability family**: it is the only member that (i) points at a **support** statement (G2's letter), and (ii) **fails the free field non-vacuously** (one-cut vs two-cut). *If* the measure-identification and regime problems were ever solved, this is the sub-family one would return to. That is a real, narrow [J]-level note — not a handle.

---

## §9 Honest boundary — what this document does NOT say

- It does **not** prove Clay and does **not** prove no tool can exist (`∀class ≠ ∀possible`; `CONVERGENCE-redteam §d`, `KEY-EXISTS §5`). It is **[NG] for the spectral-curve/free-convolution class + [O] the exact missing step + [J]** the refinement of *which* wall.
- It does **not** claim the wrong-direction (`S_N ≥ S_∞`) objection applies — it explicitly **does not** (§4). The spectral-curve species is *not* primarily a wrong-direction species; it is a **measure-mismatch + free-boundary-circularity** species. Conflating it with the free-energy large-N tools would be the imprecision this report corrects.
- The **2D GWW/DK witness** (§2) is the load-bearing fresh point: a rigorous spectral-curve support gap that provably is **not** a `0⁺⁺` mass gap. Its external-model facts are **[R]/[standard]** (exact-solvability, third-order transition); the *use* of it as a curve-gap-detector acid test is **[V-here]**.
- Bibliographic details of the free-probability theorems are **[R]/[standard]**, not re-verified this session: Voiculescu (R-transform), Voiculescu–Dykema–Nica, Nica–Speicher, Biane (free Brownian motion), Bercovici–Voiculescu (smoothing / unbounded support), Belinschi–Bercovici (support/atoms of free convolutions), Ben-Arous–Guionnet (large deviations). Makeenko–Migdal, Douglas–Kazakov, Gross–Witten, Wadia, Ginsparg–Zinn-Justin are cited by the eras/volumes already in-repo (Ginsparg–Zinn-Justin exact volume **[R]**, not re-verified). No fabricated bound; no QED.

---

## §10 References

**Internal (absolute paths; the banked verdicts this report classifies against):**
- `D:\Dev\math\CLAUDE.md`; `D:\Dev\math\docs\BRIEF-FOR-MODEL.md`; `D:\Dev\math\research\unconditional\HISTORICAL-TOOL-SWEEP.md` (§2 item 9 RMT/free-prob; §3 corner); `D:\Dev\math\research\unconditional\E2-REGRESS-CIRCULARITY.md` (§0.1 support-not-atom; §0.2 L1; §4 item 2 area-law⇏gap, §4 item 7 large-N).
- `E3-LOOPSPACE.md` (Giles ⇒ loop generator = field generator; `Γ₂`=Savvidy); `DYSON-GAS.md` (§1 Vandermonde on holonomy not energies; §2 2D GWW/DK = free energy; §3 edge location = potential = gap; §4 LD tilt = H1); `E4-RMT-EDGE.md` (RMT edge = random-Schrödinger tachyon edge = anti-gap; conditioning = H1).
- `k0-FREEPROB-AVERAGE.md` (free CLT pins the `1/V` aggregate, not the `k→0` gap; normalization fence); `H1-T8-free-prob.md` (freeness = N=∞ language; trace static vs gap dynamical; transfer no theorem); `PROOF-LARGE-N-RIGOROUS.md` (Chatterjee strong-coupling master field; `N=∞→3` asymptotic no bound).
- `ABC-DUALITY.md` (MM/1-form exact structures = 跨道/错向/错区; min–max sector = upper); `k0-SPECTRAL-REDUCTION.md` (tail's gap = the gap); `PROOF-DISSOLUTION-SPECTRAL-DENSITY.md`; `PROOF-COARSE-POINCARE.md` (H1 = coarse Poincaré = bare gap); `PROOF-1FORM-ANOMALY.md` (1-form blind to local `tr F²`); `BS-THRESHOLD.md` (area-law⇏gap; `√σ` 1-form); `SUPPORT-NOT-ATOM.md`; `PROOF-GAUSS-INVARIANT-PF.md` (vacuum pinned at 0); `GAP-G2-edge.md` (moments upper-only, edge = non-moment datum).

**External — [standard]/[R] (textbook; not re-verified this session):**
- Makeenko–Migdal, *Phys. Lett. B* **88** (1979) 135; Migdal, *Sov. Phys. JETP* **42** (1976) 413; 't Hooft, *Nucl. Phys. B* **72** (1974) 461.
- Gross–Witten, *Phys. Rev. D* **21** (1980) 446; Wadia, *Phys. Lett. B* **93** (1980) 403; Douglas–Kazakov, *Phys. Lett. B* **319** (1993) 219; Ginsparg–Zinn-Justin, 1990–93 matrix-model / double-scaling work **[R]** (exact volume not re-verified); Eguchi–Kawai, *Phys. Rev. Lett.* **48** (1982) 1063; Gopakumar–Gross, hep-th/9411021 (2D YM master field).
- Voiculescu, *J. Funct. Anal.* **66** (1986) 323 (R-transform / free convolution); Voiculescu–Dykema–Nica, *Free Random Variables* (AMS 1992); Nica–Speicher, *Lectures on the Combinatorics of Free Probability* (2006); Biane, *Indiana Univ. Math. J.* **46** (1997) 1 (free Brownian motion / support evolution); Bercovici–Voiculescu, *Indiana Univ. Math. J.* **42** (1993) 733 (free convolution with unbounded support / smoothing); Belinschi–Bercovici, support and atoms of free convolutions **[R]**; Ben-Arous–Guionnet, *Comm. Math. Phys.* **218** (2001) 177 (large deviations / LD tilt).
- Erdős–Yau–Yin, arXiv:1212.0833 (rigidity / edge universality); Tracy–Widom, *Comm. Math. Phys.* **159** (1994) 151; Bleher–Its, *Ann. Math.* **150** (1999) 185 (soft edge).
- Chatterjee, *Comm. Math. Phys.* **366** (2019) 203, arXiv:1502.07719; Jafarov, arXiv:1610.03821; Chatterjee–Jafarov, arXiv:1604.04777; Chatterjee, arXiv:2006.16229; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Shen–Zhu–Zhu, *Comm. Math. Phys.* **400** (2023) 805 (P3 strong-coupling).

---

*Status: **[NG] for the spectral-curve / free-convolution class as a non-circular `0⁺⁺`-gap handle + [O] the exact missing step (O-FP: a physical-coupling 4D master field whose eigenvalue measure is hard-gapped in the `0⁺⁺` channel and whose gap is independent of the effective action) + [J] the refinement (this species dies of measure-mismatch + free-boundary circularity, NOT primarily of the `S_N ≥ S_∞` wrong-direction objection; its support-gap detector is the family's cleanest acid-test-clean shape but identifies the wrong measure).** Load-bearing fresh point: the 2D GWW/DK witness — a rigorous spectral-curve support gap that is demonstrably not a `0⁺⁺` mass gap. NOT a proof of Clay; NOT a global impossibility theorem. Clay stays open.*

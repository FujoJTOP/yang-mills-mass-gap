# SMALLCIRCLE-ADIA — R³×S¹ small circle + adiabatic continuity: a softer corridor to the 4D pure-SU(N) mass gap, or the same wall on the radius axis?

> 2026-11-23 · Assessment of the small-circle / deformed-YM / magnetic-bion route as a *distinct door* to Clay, audited against the project's iron rules (label 真定理/条件定理/no-go/精确缺失步; free-massless acid test). Consolidates and sharpens `R1-small-S1-bions`, `R2-R4-limit`, `FORGE5-SMALL-S1`, `FORGE5-DUAL-SUPERCONDUCTOR`, `EGUCHI-KAWAI`, `P4-RP-uniform-gap`, `MONOTONE-COUPLING`, `W2-transmuted-gap`.

## Verdict (up front)

**SAME WALL, relocated to the circle radius — and structurally the *worst* direction: an upper-bound machine.** There is a genuinely softer *near end* (a real small parameter `e^{−S₀}` and controlled semiclassics), but the hard step `r→∞` is the mass gap itself, and **the only free inequality is`m(S¹_r) ≥ m(ℝ⁴)`, which runs backward for a floor.** Concretely:

- **What is rigorously proven at small r for 4D pure SU(N): essentially nothing about a gap.** Pure YM on a small periodic circle is **center-broken / deconfined** (GPY one-loop), so its small-r "confined phase" does not exist without a deformation. The center-symmetric small-r theory is **deformed YM / QCD(adj) / 't Hooft-twisted** — a *different theory*. The bion-induced dual-photon mass is **controlled semiclassics, not a theorem**: no rigorous continuum SU(N) measure on ℝ³×S¹ exists at all, and the only rigorous "monopole gas ⇒ area law" is **abelian** (Göpfert–Mack).
- **The continuity step.** `Δ(r) := gap(ℝ³×S¹_r)` is **non-increasing** in `r` (min–max / finite-volume upper bounds: `Δ(ℝ⁴)=inf_r Δ(r)`). Hence a small-r lower bound `Δ(small) ≥ δ` certifies only `Δ(ℝ⁴) ≤ Δ(small)` — an **upper** bound. To floor `Δ(ℝ⁴) ≥ δ` you must prove `Δ(r) ≥ δ` for **all** r (small *and* large); the large-r end is the gap. This is exactly the campaign's own L\* / transport wall (`P4 §5`, `MONOTONE-COUPLING §2`), transported to the r-axis. **精确缺失步.**
- **Acid test: PASSES (not vacuous).** The free massless field on ℝ³×S¹ has no monopoles/bions ⇒ `m_γ=0`, and the zero KK mode is a massless 3D photon ⇒ genuine gap 0. The mechanism correctly fabricates **no** gap. (Guard: the KK spacing `1/r` is a finite-volume artifact, not the mass gap; the acid test is applied to the zero-KK-mode mass.)
- **Not the same theory; (CP) does not help.** Bridging deformed/twisted small-r to pure YM on ℝ⁴ requires taking (deformation→0)∘(r→∞) — both of which *are* the gap; and (CP) ("no phase transition as r grows") is precisely adiabatic continuity, which is **false** for pure YM (GPY transition) and, in the closest analog QCD(adj), lattice finds **four phases** (Cossu–D'Elia), i.e. no smooth connection.

**One-line answer:** *Not a softer corridor. A softer near-end grafted onto the same backward inequality; the floor at large r is the mass gap. The different variable (r) buys real semiclassical control but points the inequality the wrong way.*

---

## §1 Problem posed precisely

Clay = **G1** (existence of 4D pure SU(N) YM) + **G2** (`spec(H)⊆{0}∪[Δ,∞)`, `Δ>0`) on the physical Hilbert space. Let `Δ(r)` be the mass gap of the theory on `ℝ³×S¹_r` (one direction compactified, circumference `2πr`… radius `r`; the lightest excitation above the vacuum). Two ends:

- **Small r** (`NΛr ≪ 1`): weak coupling, semiclassical. GPY one-loop, monopole-instantons, magnetic bions, dual photon, area law — Ünsal–Yaffe programme.
- **Large r** (`r→∞`): `ℝ⁴`, strong coupling — the physical theory, the gap.

"Adiabatic continuity" (Ünsal–Yaffe, `0803.0344`) is the **conjecture** that `Δ(r)` connects smoothly (no phase transition) from the controlled small-r phase to `ℝ⁴`. The route proposes: *prove the gap at small r, then transport it to `r=∞`.* The assessment is whether that transport is (a) done, (b) doable, and (c) pointing the right way.

---

## §2 What is ACTUALLY proven at small r (label every item)

**真定理 / rigorous operator-level statements:**
- **Osterwalder–Seiler** (`Ann. Phys. 110 (1978) 440`): the lattice Wilson measure is reflection-positive; a self-adjoint transfer matrix `T=e^{−aH}`, `0≤T≤1`, `H≥0`, and a Källén–Lehmann measure exist. `RP` fixes *signs*, **not** the size of `1−λ₁`. This is a transfer operator — but on the **lattice at finite spacing**, and it yields no gap. (It is the same machinery the whole campaign used; it does not specialise to small r.)
- **Center symmetry, the Polyakov/'t Hooft loop algebra, and the classification of twisted sectors / monopole-instantons** are rigorous operator/lattice statements (Elitzur; 't Hooft flux). What is *not* rigorous is any dynamical consequence.

**Controlled one-loop (算得，不是定理):**
- **GPY (`RMP 53 (1981) 43`) one-loop effective potential for the Polyakov loop.** Its load-bearing, **honest** content here is the *opposite* of the task framing: for **pure YM with periodic (= thermal) boundary conditions, the one-loop potential favours the center-broken (deconfined) configuration.** Small periodic circle of pure YM is **deconfined**, `⟨tr Ω⟩ ≠ 0`. Center symmetry at small r is restored only by **(i)** periodic adjoint fermions (`QCD(adj)`), **(ii)** a **double-trace deformation** (Shifman–Ünsal `0802.1232`; Ünsal–Yaffe `0803.0344`), or **(iii)** an 't Hooft twist. *(This corrects the task's parenthetical "pure YM: center symmetric" — for pure YM it is center-broken; only the deformed/twisted/adjoint theories are center-symmetric small-r.)* [真定理 for the one-loop computation; its dynamical extrapolation is not.]

**半经典 / semiclassical (the bion programme — a small parameter, NOT a proof):**
- **Monopole-instantons** (`Affleck–Harvey–Witten, NPB 206 (1982) 413`; Lee–Yi, van Baal) and **magnetic bions** (`Ünsal, PRD 80 (2009) 065001, arXiv:0709.3269`) generate a **dual-photon mass** `m_γ ~ e^{−S₀}`, `S₀ = 8π²/(g²N) ≫ 1`, a 3D area law, and a gap. Diluteness is `e^{−S₀}`; determinants and moduli are **computed, not bounded**; the bion is a **marginal `1/r`** object (a threshold molecule, not a clean separated particle). `真定理: no.`

**What is rigorously absent:** an SU(N) continuum measure on ℝ³×S¹ at small r; a convergent cluster expansion over the interacting nonabelian monopole–bion gas giving area law + exponential clustering. The single rigorous "monopole gas ⇒ area law" is **abelian** (`Göpfert–Mack, CMP 82 (1982) 545`, 3D U(1) lattice). The UV machine that could assemble it (Balaban) stops at small fields / small volume. **结论：small-r gap for pure SU(N) is 半经典/heuristic, not 真定理.**

---

## §3 The hard step: the small-r → ℝ⁴ continuity (and why it is the same wall, backward)

**Lemma (decompactification monotonicity of the gap).** For `r₂>r₁`: `Δ(r₂) ≤ Δ(r₁)`. *Sketch (min–max).* The gap is the infimum of the Rayleigh quotient over unit vectors `⊥ Ω`; a larger circle admits a larger space of trial states (more KK/Fourier modes), so the minimum over the larger space is `≤` the minimum over the smaller. Equivalently (`P4 §2`, exact): `Δ(ℝ⁴) = inf_r Δ(r)`, and every finite-r gap is an **upper bound** on `Δ(ℝ⁴)`. ∎

**Consequence (the backward inequality).** Suppose — granting the unproven small-r gap — that `Δ(r_small) ≥ δ` is proven. By the Lemma this certifies only
```
Δ(ℝ⁴) ≤ Δ(r_small)          (an UPPER bound),
```
not `Δ(ℝ⁴) ≥ δ`. To floor `Δ(ℝ⁴)` one needs `Δ(r) ≥ δ` **uniformly in r, including `r→∞`** — which is `Δ(ℝ⁴) ≥ δ`, the mass gap. **This is exactly the campaign's residual object L\* / the transport wall**, written on the r-axis instead of the β-axis.

**Corroboration via the screening mass (temporal-circle language).** At `T=1/r`, the spatial screening mass satisfies `m_scr(T) ≥ m_scr(T{=}0) = m_gap` and is **non-increasing in r**. A *lower* bound on `m_scr(small r, high T)` is therefore a value *above* `m_gap` and cannot floor it. Both formulations agree; the operative object is the one the semiclassical expansion computes (the finite-r gap / dual-photon mass), and it sits at the **more-gapped end**.

**Is this the SAME wall as the main corridor?** Yes — structurally identical, and the corpus already named the isomorphism: the controllable end is the **more-gapped** end, and every monotone comparison points *toward* it, so the proven value is an upper bound (`MONOTONE-COUPLING §2`: gap decreasing in β / increasing in deformation ⇒ proven references are upper bounds; `PROOF-TRANSPORT-MONOTONE`: transport ≡ confinement). Replacing β by r does not change the sign of any inequality. **Verdict on Q2: same wall, not genuinely different; the different *variable* buys semiclassical control at one end only.**

**Could it be made rigorous — construct on ℝ³×S¹, gap uniform as r→∞?** Two independent gaps, both labeled:
1. **Construction (`G1`-type).** No rigorous reflection-positive continuum SU(N) measure on ℝ³×S¹ at small r exists (`R1`). The abelian Göpfert–Mack gas + Balaban UV control are ingredients, never assembled for nonabelian YM. **精确缺失步 (construction).**
2. **Uniform floor (`G2`).** Even granting (1) and the small-r gap, a **uniform-in-r lower bound** `Δ(r) ≥ δ` is required, and by the Lemma its large-r end *is* `Δ(ℝ⁴) ≥ δ`. Completing this route proves the **upper** bound `Δ(ℝ⁴) ≤ Δ(small)` unless it *simultaneously* proves the large-r floor — i.e. the gap. **精确缺失步 (= Clay G2).**

---

## §4 Acid test: free massless gauge field on ℝ³×S¹

Take the free (abelian, massless) gauge field on ℝ³×S¹. **Correct output: gap = 0, at every r.**

- **Zero-KK-mode route:** the massless 3D photon is gapless (continuum down to 0); there is no mass gap.
- **Bion route:** the free/`U(1)` field with no monopole condensation has **no monopole-instantons and no bions**; the dilute-gas dual-photon mass term is **absent**, so `m_γ = 0`, `f₁² = 0`. The resummation correctly returns **no mass** — it does not fabricate a gap.

**Outcome: PASS — the route is not vacuous on the free massless field.** This is the one genuine pass (`R1`, `FORGE5-SMALL-S1 #4`, `FORGE5-DUAL-SUPERCONDUCTOR #4`, `PIROGOV-SINAI §3`). Two honest guards:

- **Spurious-`1/r` trap.** One could mistakenly read "the KK gap `~1/r`" as a mass gap. That is finite-volume quantisation of the *compact* momentum, not a dynamical mass: it →0 as `r→∞` and is present for the free field too. The acid test must be applied to the **mass gap of the zero-KK-mode sector**, where the free field is correctly gapless.
- **Pass is necessary, not sufficient.** What generates the gap (bion dual-photon mass) lives in the **magnetic / dual** channel (`√σ`, 1-form center), whereas the Clay residue `f₁²=|⟨0|tr F²|G⟩|²` is the **neutral `0⁺⁺`** channel. `tr F²` is N-ality-0 and provably **blind** to the 1-form center symmetry (`PROOF-1FORM-ANOMALY`). So the object that passes the acid test is the wrong channel for the residue (`FORGE5-DUAL-SUPERCONDUCTOR §5`, `H1-T10-dual`).

---

## §5 Is the small-circle confined phase the *same* theory as ℝ⁴ pure YM? Can (CP) help?

**No — not the same theory.** On ℝ³×S¹:
- **Exact compactification, periodic BC, pure YM** = pure YM, but at small r it is **center-broken / deconfined** (GPY). So the **confined** small-r phase is *not* pure YM.
- The **center-symmetric confined small-r phase** requires a **different theory**: **deformed YM** (double-trace, Shifman–Ünsal; Ünsal–Yaffe), **QCD(adj)** (periodic adjoint fermions), or an **'t Hooft twist**. Each is *not* ℝ⁴ pure YM. The deformation is, by construction, "irrelevant at large L" **only under the adiabatic-continuity conjecture** — i.e. taking deformation→0 returns pure YM *only if the continuity holds*.

**The bridge is a composition of two infima, both = the gap.** `deformed small-r → pure ℝ⁴` requires `(deformation→0) ∘ (r→∞)`. Deformation→0 is an infimum over the deformation strength (every `m₀>0` reference is an upper bound; `MONOTONE-COUPLING §4`), and `r→∞` is the mass gap. So identifying the phases *is* the missing step, not a step around it. **条件定理 at best; the identification is circular.**

**(CP)-type arguments.** "(CP)" in the campaign = *single T=0 confining phase (no bulk transition) + `0⁺⁺` edge Kato-preservation + uniform isolation* — one of the equivalent residual objects, i.e. = the gap. In this route, the (CP)-content is exactly "no phase transition as r grows" = **adiabatic continuity itself.** Two honest obstacles:
- It is **false for pure YM**: the GPY one-loop + lattice give a genuine **deconfinement transition at `L_c=1/T_c`**; deleting it requires the deformation (`FORGE5 §3`).
- It is **not robust even for the analog**: in `QCD(adj)` lattice finds **four phases** with two intermediate center-broken phases between the small- and large-`L` confined phases (`Cossu–D'Elia, JHEP 2009, arXiv:0904.1353`) — i.e. **no** smooth connection in the closest controlled analog. (Trace-deformed YM "reconfinement" — `Bonati–Cardinali–D'Elia, PRD 98 054508 (2018); PRD 101 034508 (2020); Bonati et al., PRD 103 034506 (2021)` — shows the *deformation* can restore center symmetry above `T_c`, which is a statement about the *deformed* theory, not continuity of pure YM.) So (CP) cannot be invoked; it **is** the conjectured bridge and has lattice evidence against it.

---

## §6 Exact missing step, and citable rigorous results

**精确缺失步 (the whole route, once):**
> Construct a reflection-positive continuum SU(N) measure on `ℝ³×S¹` at small r *with center stabilization*, prove a convergent cluster expansion over the nonabelian monopole–bion gas giving `Δ(r) ≥ δ` at small r, **and** prove `Δ(r) ≥ δ` uniformly as `r→∞`. By the Lemma (§3) the uniform statement *is* `Δ(ℝ⁴) ≥ δ`. **The transport step is indistinguishable from the mass gap** — the same object as L\*, and it points backward (small-r bound ⇒ upper bound on `Δ(ℝ⁴)`).

**Rigorous results one may cite (each with its honest limitation):**
- **Göpfert–Mack, CMP 82 (1982) 545** — rigorous monopole-gas ⇒ area law, **3D abelian U(1) only**.
- **Osterwalder–Seiler, Ann. Phys. 110 (1978) 440** — rigorous OS transfer matrix + KL measure; **no gap**, lattice/finite-spacing.
- **Shen–Zhu–Zhu, CMP 400 (2023) 805, arXiv:2204.12737** — the only unconditional 4D pure-YM lattice gap, **strong coupling `β<1/48`** (wrong regime; no continuum, no semiclassical).
- **GPY, RMP 53 (1981) 43** — one-loop Polyakov potential; shows pure YM small-r **center breaking**.
- **2D `CP^{N-1}` / `O(N)` "continuity": Sulejmanpasic, PRL 118, 011601 (2017), arXiv:1610.04009.** *Physics-level* confirmation (not a rigorous math theorem) that with **symmetry-twisted** boundary conditions these 2D sigma models avoid the Affleck transition; reduces to a smooth effective QM. **Does not transfer:** (i) 2D-specific (dimensional reduction is drastic), (ii) the **twist is essential** — thermal BC has the Affleck transition, (iii) no rigorous 4D analog; in 4D the twisted-sector/anomaly data constrain the *vacuum structure*, not the `0⁺⁺` **gap** (cross-channel: `PROOF-1FORM-ANOMALY`).
- **Lattice adiabatic-continuity tests (numerical, not proof):** `Cossu–D'Elia 2009` (four phases, *against* smooth connection); `Bonati et al.` (trace-deformed reconfinement). *(Note on the task's "Davies et al": I could not verify any "Davies et al" paper on small-circle adiabatic continuity; the relevant lattice works are Cossu–D'Elia / Bonati–Cardinali–D'Elia. Flagged per the citation-checking rule; not propagated.)*

---

## §7 Five-point scorecard

| # | Criterion | Verdict |
|---|---|---|
| 1 | **Nonperturbative** | **PASS (method), FAIL (theorem).** Bion gas is genuinely nonperturbative (`e^{−S₀}`); but *no rigorous small-r SU(N) gap* exists — determines are computed, bion is marginal, the only rigorous gas is abelian. Nothing to transport. |
| 2 | **Sign-controlling LOWER bound `m²(k→0)≥λ`** | **FAIL.** The sign fixed at small r is the 3D dual-photon / monopole-gas sign — a confinement-scale object; it never touches the 4D Savvidy-cone Hessian `−D²+2ad F_B` (`λ_min=−g|B|<0`) which *is* the wall. Wrong object. |
| 3 | **Non-circular** | **FAIL.** Undeformed pure YM small-r is center-broken (GPY); the controlled theory is deformed/twisted = *not* pure YM; identifying its large-r limit with pure YM is adiabatic continuity = the gap. Circular. |
| 4 | **Non-SUSY + fails free field** | **PASS (sole genuine pass).** Non-SUSY monopoles; free/abelian field has no bion ⇒ `m_γ=0`, `f₁²=0` — correctly no gap. Caveat: it is the dual/screening mass, not the neutral `0⁺⁺` residue (wrong channel). |
| 5 | **Uniform lower bound `≥λ` in `(L,a)` / all r** | **FAIL (structurally).** See §3: `Δ(r)` non-increasing, `Δ(ℝ⁴)=inf_r Δ(r)`; small-r bound ⇒ **upper** bound on `Δ(ℝ⁴)`. The uniform floor is L\* = the gap. |

**Score: 1 genuine pass (#4), #1 pass-method/fail-theorem, #2/#3/#5 are one object — the backward decompactification inequality.**

**Same wall or softer corridor?** Same wall. The near end is genuinely softer (a true small parameter, controlled semiclassics — an *advantage* over the lattice strong-coupling end), but "soft at one end" is not a soft corridor: the hard step is `r→∞`, it is the mass gap, and **no inequality reverses its direction**. The variable change β↔r relabels the axis without flipping the sign.

---

## §8 Sharpest line

> **Small S¹ is an upper-bound machine:** min–max forces `Δ(S¹_r) ≥ Δ(ℝ⁴)`, non-increasing in r — the better you control small r, the more you **cap** the ℝ⁴ gap, never **floor** it. No RP, no adiabatic-continuity gesture, and no 't Hooft-anomaly argument reverses the direction; the floor is `Δ(ℝ⁴) ≥ δ` at large r, which is the mass gap. The route's one real asset — semiclassical control at small r — is real, but it sits on the **more-gapped** end of a monotone quantity that points the wrong way.

---

## §9 References (verified 2026-11-23)

**External:**
- D. J. Gross, R. D. Pisarski, L. G. Yaffe, *QCD and instantons at finite temperature*, Rev. Mod. Phys. **53** (1981) 43. [GPY]
- I. Affleck, J. Harvey, E. Witten, *Instantons and (super-)symmetry breaking in (2+1) dimensions*, Nucl. Phys. **B206** (1982) 413–439. [monopole-instanton / Polyakov mechanism]
- M. Ünsal, *Magnetic bion condensation: a new mechanism of confinement and mass gap in four dimensions*, Phys. Rev. D **80** (2009) 065001, arXiv:0709.3269.
- M. Ünsal, *Abelian duality, confinement, and chiral symmetry breaking in QCD(adj)*, Phys. Rev. Lett. **100** (2008) 032005, arXiv:0708.1772.
- M. Shifman, M. Ünsal, *QCD-like theories on R³×S¹: a smooth journey from small to large r(S¹) with double-trace deformations*, Phys. Rev. D **78** (2008) 065004, arXiv:0802.1232.
- M. Ünsal, L. G. Yaffe, *Center-stabilized Yang–Mills theory: confinement and large-N volume independence*, Phys. Rev. D **78** (2008) 065035, arXiv:0803.0344.
- G. 't Hooft, *A property of electric and magnetic flux in nonabelian gauge theories*, Nucl. Phys. **B153** (1979) 141. [twisted sectors]
- K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. **110** (1978) 440. [OS transfer matrix]
- C. Göpfert, G. Mack, *Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory…*, Commun. Math. Phys. **82** (1982) 545. [rigorous abelian monopole gas]
- T. Sulejmanpasic, *Global symmetries, volume independence, and continuity in quantum field theories*, Phys. Rev. Lett. **118** (2017) 011601, arXiv:1610.04009. [2D CP^{N-1}/O(N), symmetry-twisted — physics-level]
- G. Cossu, M. D'Elia, *Finite size phase transitions in QCD with adjoint fermions*, JHEP **0907** (2009) 048 / arXiv:0904.1353. [four phases — against smooth continuity]
- C. Bonati, A. Cardinali, M. D'Elia, *θ-dependence in trace-deformed Yang–Mills*, Phys. Rev. D **98** (2018) 054508; Phys. Rev. D **101** (2020) 034508; Bonati et al., Phys. Rev. D **103** (2021) 034506. [trace-deformed reconfinement lattice]
- Shen, Zhu, Zhu, *…mass gap for lattice Yang–Mills*, Commun. Math. Phys. **400** (2023) 805, arXiv:2204.12737. [strong-coupling 4D gap]
- E. Poppitz, *Notes on the confined phase / deformed YM*, review, arXiv:2111.10423. [programme's own "conjecture" status]

**Internal (this project):** `R1-small-S1-bions` (rigor status), `R2-R4-limit` (ℝ⁴ limit open), `FORGE5-SMALL-S1` (RP/no-crossing continuity — backward inequality), `FORGE5-DUAL-SUPERCONDUCTOR` (channel/circularity), `EGUCHI-KAWAI` (same center-breaking wall in volume form), `P4-RP-uniform-gap §2/§5` (`Δ_{ℝ⁴}=inf_V Δ_V`, L\*), `MONOTONE-COUPLING §2/§4` (backward direction; deformation infimum), `PIROGOV-SINAI §3` (screening = susceptibility; free-field pass), `W2-transmuted-gap (iv)` (bions on ℝ³×S¹ only), `PROOF-1FORM-ANOMALY` (`tr F²` blind to 1-form center), `PROOF-SAVVIDY-EXACT-SIGN` (`λ_min=−g|B|`), `CONTINUUM-3D-TO-4D` (construction wall), `H1-T10-dual` (Meissner cross-channel).

---

## §10 Bottom line for the campaign

This route is an **audit-negative with a genuine structural asset**, not a new door: it reproduces the campaign's own wall (`P4`/`MONOTONE-COUPLING`/`PROOF-TRANSPORT-MONOTONE`) on the **r-axis**, with the same wrong-direction signature (control the *more-gapped* end; the floor is the *less-gapped* physical end). Its distinguishing merit — a real semiclassical small parameter `e^{−S₀}` at small r — does **not** soften the hard step, because the step is one-directional: `r→∞` uniformity = `Δ(ℝ⁴) ≥ δ` = Clay G2. The single open sub-question worth one further probe (not a corridor) is whether a **4D analog of the Sulejmanpasic twisted-BC mechanism** could yield a *gap* (not merely vacuum/anomaly data) without the backward inequality; the honest expectation is **fatal** — 2D-specific reduction, twist-essential, and cross-channel for the `0⁺⁺` residue. **The route does not move the wall; it sharpens where the wall sits on the radius axis.**

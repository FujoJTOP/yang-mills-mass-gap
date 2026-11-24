# DEFORM-DECOUPLE-4D — can a "construct-then-decouple" deformation deliver the 4D pure-SU(N) mass gap?

> 2026-11-24 · Strategic research scout on the one constructive route: find a **4D constructible deformation** of pure Yang–Mills whose **decoupling/undeforming limit is pure YM with the same IR mass gap**, then transfer the construction's gap.
>
> **Labels (charter §2).** **[T]** 真定理 / true theorem (unconditional within its stated domain) · **[C]** 条件定理 / conditional theorem · **[NG]** no-go **for a named class** (not for the problem) · **[O]** 精确缺失步 / exact missing step / open · **[J]** 判断 / judgment · **[R]** external/library fact — *re-verified this session* where marked, otherwise standard/as-reported · **[V-here]** fresh computation this session.
>
> **Scope / honesty.** This document does **NOT** prove Clay, and does **NOT** prove that no route can exist (`∀class ≠ ∀possible`, KEY-EXISTS §5). It is a **[NG] for the enumerated deformation classes**, an **[O] map of the exact missing steps**, and a **[J] strategic verdict**. No QED on anything unproven.

---

## VERDICT (up front)

**No deformation in the enumerated classes satisfies (i)∧(ii)∧(iii).** Precisely: there is **no** 4D deformation that is simultaneously **(i)** rigorously constructible, **(ii)** has a *proven dynamical gap of the deformed nonabelian theory*, and **(iii)** has a *proven gap-preserving limit to pure YM*. The candidates split into two conjugate classes, each dying for its own reason — and the two deaths are the same obstruction seen from opposite ends.

- **Class A — explicit-scale / weak-coupling deformations** (mass-deformed, softly-broken, Gaussian 4D YM–Higgs). These **do** have a gap, but the gap is the **deformation's own mass** (Proca / Higgs-mechanism / gaugino), not the nonabelian dynamics. The pure-YM limit **removes the very scale that produced the gap**, so the gap closes; and they **fail the free-field acid test in the fatal direction** — they give a gap to a *free abelian* field too, hence are not detecting the YM gap. **[NG] for the class.**
- **Class B — "same-theory" deformations** (center-stabilized / twisted small-`S¹`, Eguchi–Kawai volume reduction). These are genuinely nonabelian and **pass** the free-field acid test, but they clear **(i),(ii),(iii) not at all** rigorously, and the **only free inequality is an upper-bound machine**: the controlled end is the *more-gapped* end, and min–max points the floor the wrong way. **[NG] for the class.**
- **Single unifying obstruction.** *The deformation → pure-YM limit is a min–max infimum over the deformation parameter (undeform `m₀→0`, decompactify `r→∞`, `N→∞`, restore the interaction as `g` is un-sent to 0), and every controlled/räsoniert statement about the deformed theory is an **upper bound** on the limit gap. A gap-preserving limit therefore requires a **uniform floor at the physical, less-gapped end** — which is exactly `m_{0⁺⁺}>0` = H1 = (CP) = the mass gap.* Equivalently: **every constructible deformation's pure-YM limit crosses the Savvidy / strong-coupling crossover at fixed nonzero coupling, and "adiabatic continuity" is precisely the unproven (CP).** The three textbook phrases — *min–max goes the wrong way*, *adiabatic continuity is conjectural*, *crosses the crossover = H1* — are one statement. **[T]/[J]** (min–max lemma **[T]**; the identification with H1 **[J]**, banked in `SMALLCIRCLE-ADIA`, `MONOTONE-COUPLING`, `E2-FLIP-SECOND-SCALE`).
- **Acid-test framing (Task 4).** A clean construct-then-decouple *must* fail on the free massless field (free is gapless). Class B passes — but being sign-correct is not enough: it is an upper-bound machine. Class A **fails in the wrong direction** (gap for a free field ⇒ wrong object). **The nonabelian interaction enters only at the Savvidy/strong-coupling crossover — exactly where no construction crosses.**

**One-line answer.** *No constructible 4D deformation has a proven gap-preserving pure-YM limit. The most promising constructible candidate is Chatterjee's 4D YM–Higgs scaling limit; its exact missing link is a non-Gaussian fixed-`g` construction plus a gap-preserving Higgs-decoupling limit, which is the 4D marginal/Savvidy wall in disguise. The unifying obstruction is that decoupling is a min–max infimum whose gap-preserving limit is the mass gap.*

---

## §1 Candidate table

Legend: **(i)** constructible? · **(ii)** gap in deformed theory? · **(iii)** decoupling → pure YM proven gap-preserving?

| deformation | (i) constructible? | (ii) gap? | (iii) decoupling proven? | obstruction | citation |
|---|---|---|---|---|---|
| **(a) 4D YM + fundamental Higgs** (fixed length, degenerate potential) | **Partial [R]**: a 4D scaling limit *exists* — but only at `g=O(ε^{50d})→0`, `αg=cε`, Higgs length `α→∞`, unitary gauge | **Trivial [R]**: the limit is a massive **Gaussian (Proca)** field; the mass is the **Higgs-mechanism** mass, not YM dynamics | **No [O]**: paper explicitly leaves non-Gaussian limits open; pure YM is not obtained | explicit Higgs mass; `g→0` kills the nonabelian interaction; fixed-`g` non-Gaussian = 4D marginal/Savvidy wall | Chatterjee, arXiv:2401.10507, *Prob. Math. Phys.* **7** (2026) 339–381 **[R, verified]** |
| **(a′) 4D YM + adjoint Higgs** | **No [O]**: no 4D construction; the 3D YM–Higgs construction is local-in-time only (no measure, no gap) | would be the Higgs mechanism; no theorem | **No**; and a **genuine phase boundary can exist** for adjoint Higgs (Fradkin–Shenker), so `m_H→∞` crosses a transition | construction wall **+** adjoint-Higgs transition ('t Hooft–Polyakov monopoles) | CCHS, arXiv:2201.03487, *Invent. Math.* **237** (2024) 541 **[R, verified]**; Fradkin–Shenker, *PRD* **19** (1979) 3682 **[R, verified]** |
| **(b) mass-deformed / softly-broken YM** | **No [O]**: a Proca mass breaks gauge invariance (no OS/rigorous 4D); soft-SUSY has no 4D construction | **Trivial [T]**: `gap ≥ m₀`; but this is the deformation's mass | **No [NG]**: `gap(m₀)` ↑ in `m₀` ⇒ `gap(YM)=inf_{m₀}`; every `m₀>0` reference is an **upper** bound; `m₀→0` removes the gap | the mass term *is* the gap; decoupling is an infimum; acid test fails (free abelian gets a gap too) | `MONOTONE-COUPLING` §2/§4; `E2-FLIP-SECOND-SCALE` §(iii) **[T-conditional]/[NG]** |
| **(c) twisted / center-stabilized small-`S¹`** (Ünsal–Yaffe, Shifman–Ünsal) | **No [O]**: no rigorous continuum SU(N) measure on `ℝ³×S¹`; bions are *computed* semiclassics, not bounded; the only rigorous gas⇒area-law is **abelian** (Göpfert–Mack, 3D) | **Semiclassical only [R]/[O]**: dual-photon mass `m_γ~e^{−S₀}`; and it lives in the **1-form / `√σ` channel**, not the neutral `0⁺⁺` residue | **No [NG]**: pure YM small-`S¹` is **center-broken** (GPY); the confined small-`r` phase is *deformed / twisted / QCD(adj)*; lattice finds **four phases**, two center-broken intermediates; `Δ(r)` non-increasing, `Δ(ℝ⁴)=inf_r Δ(r)` | deformation = different theory; center breaking of pure YM; backward min–max inequality; cross-channel | Ünsal–Yaffe, arXiv:0803.0344, *PRD* **78** (2008) 065035 **[R, verified]**; Shifman–Ünsal, arXiv:0802.1232, *PRD* **78** (2008) 065004 **[R, verified]**; Cossu–D'Elia, arXiv:0904.1353, *JHEP* 0907:048 **[R, verified]**; `SMALLCIRCLE-ADIA` **[T]/[NG]** |
| **(d) SUSY-adjacent** (N=2 SW; N=1; N=1\*; N=4) | **No [O]**: no 4D *constructive* construction; SW is exact semi-classical/duality, not an OS measure | **Wrong channel / no gap [R]**: N=2 Coulomb branch gapless; BPS atoms are **charged** (`Z≠0`), while `tr F²` has `Z=0`; N=4 conformal gapless; N=1 gaugino-condensation "gap" is EFT-level, unproven | **No [NG]**: N=2→0 soft breaking destroys protection and crosses the SW point; N=1→0 `m_gluino→∞` is physics, not a theorem | charge-sector orthogonality (`tr F²` blind to `Z`); soft breaking; no proof | Seiberg–Witten, hep-th/9407087, *NPB* **426** (1994) 19 **[R, verified]**; `BPS-SUSY-PROTECTION` **[T]/[NG]** |
| **(e) spatial compactification / Eguchi–Kawai** | **No [O]**: volume independence is an `N=∞` statement; fails when center breaks (Bhanot–Heller–Neuberger); no theorem at finite `N` | **Not a spectral gap [R]/[NG]**: the one-site model is a **0D matrix integral** — no Euclidean time, no Källén–Lehmann measure; `N=∞→3` has no theorem | **No [NG]**: exact reduction ⇒ Gaussian slice `D=0` ⇒ `gap_red ≥ gap_true` (upper bound); center breaking = same small-`S¹` wall | volume independence unproven; no time direction; `N=∞→3`; min–max backward | Eguchi–Kawai, *PRL* **48** (1982) 1063 **[R, verified]**; Kovtun–Ünsal–Yaffe, hep-th/0702021, *JHEP* 06 (2007) 019 **[R, verified]**; `EGUCHI-KAWAI` **[NG]** |
| **(ref) 4D lattice strong coupling** (SZZ) | **Yes [T] — but not a decoupling route** | **Yes [T]**: exponential decay / mass gap for `β<1/48` | **No**: wrong regime (opposite to asymptotic freedom); `a→0`, `β→∞` is the crossover wall | the only unconditional 4D gap lives at strong coupling; no continuum limit | Shen–Zhu–Zhu, arXiv:2204.12737, *CMP* **400** (2023) 805–851 **[R, verified]** |

**Notes on the table (each a claim).**
- **[R] (a) construction is Gaussian, not interacting.** Chatterjee's theorem is a *scaling limit* of lattice SU(2) YM + fundamental Higgs; the limit field is a **massive Gaussian** ("converges to a massive Gaussian field"), and the paper says "the question of constructing a non-Gaussian scaling limit remains open." So (i) is cleared only for a **free** limit, and (ii) is a **free Proca mass**, not a dynamical YM gap. **[R, verified]**
- **[R] (a′) 3D is not a template.** The 3D YM–Higgs construction (CCHS) is **local-in-time**, with no invariant measure, no OS positivity, no gap. Even in 3D the strategy does not already possess a gap to transfer. **[R, verified; `CONTINUUM-3D-TO-4D`]**
- **[T] (b) direction is decisive, independent of existence.** `gap(YM) = lim_{m₀→0} gap(m₀)` is an **infimum**; every proven `gap(m₀>0)` is an upper bound (`MONOTONE-COUPLING §4`). A mass deformation therefore cannot floor the undeformed gap.
- **[T] (c) min–max lemma.** For circumference `r`: `Δ(r)` is non-increasing in `r` and `Δ(ℝ⁴)=inf_r Δ(r)` (finite-volume gap ≥ infinite-volume gap). Hence a small-`r` lower bound certifies only `Δ(ℝ⁴) ≤ Δ(small)` — an **upper** bound. Completing the route needs `Δ(r)≥δ` for **all** `r`, including `r→∞`, which *is* the gap. (`SMALLCIRCLE-ADIA §3`.)
- **[R] (c) adiabatic continuity is not a theorem.** Ünsal–Yaffe present center stabilization as making volume independence hold and small-`L` dynamics "tractable analytically," with **no order parameter separating** the small- and large-radius regimes — a continuity *assumption*. Poppitz's review presents "evidence" against phase transitions, not a proof. Lattice QCD(adj) finds **four phases** with two center-broken intermediates (Cossu–D'Elia), i.e. **no smooth connection** in the closest controlled analogue. **[R, verified]**
- **[NG] (e) reduction is the wrong direction.** The ledger records the general fact: any *exact* reduction gives a Gaussian slice (`D=0`), hence `gap_red ≥ gap_true` — an upper bound (`ABC-FINAL-LEDGER` §A). Eguchi–Kawai adds the loss of a time direction and the unproven `N=∞→3` transfer.

---

## §2 The most promising candidate and its exact missing link

### §2A Most promising: **4D YM + fundamental Higgs** (Chatterjee) — the only class with a real 4D construction

**Why it is the (comparative) best.** It is the **only** deformation in the table that clears (i) *at all*: there is a published, peer-reviewed 4D scaling limit of nonabelian lattice gauge theory coupled to a Higgs (Chatterjee 2026). It also has the structural advantage that, for a **fundamental** Higgs, the Higgs and confinement regimes are the **same phase** (Fradkin–Shenker complementarity) — so, unlike the adjoint case (a′), there is **no known phase boundary** that the decoupling limit must cross. **[R, verified]**

**What is actually constructed. [R, verified]** `SU(2)` lattice Yang–Mills + a **fixed-length fundamental** Higgs (degenerate potential), in unitary gauge, with `ε→0`, `g=O(ε^{50d})→0`, `α→∞`, `αg=cε`; a stereographic projection of the gauge field **converges to a massive Gaussian field**. The paper states it is "the first rigorous proof that the Higgs mechanism generates mass," and explicitly: "the question of constructing a non-Gaussian scaling limit remains open."

**Exact missing link (the whole route, once). [O]**
> Construct the 4D YM–Higgs theory **non-Gaussian, at fixed nonzero gauge coupling `g`**, and prove that its gap **survives the Higgs-decoupling limit** (`α→0` / `m_H→∞`) to pure YM.

Two independent gaps, both named:
1. **Construction [O].** Fixed-`g` 4D YM–Higgs is **marginal**: `[g²]=0`, the scale-sum is harmonic, and the tree-level degree condition of regularity structures fails. The only solved marginal SPDE (4D `φ⁴`) is **trivial** (Aizenman–Duminil-Copin), so it is no template. **[R, verified]**
2. **Gap survival [O].** The controlled limit has `g→0`; the "gap" is the **Gaussian Proca mass**, and the nonabelian interaction is gone. Proving the gap at fixed `g`, then decoupling the Higgs, is **not** the same problem — and the transverse gauge sector is unchanged by the scalar: the Savvidy Hessian `−D²+2 ad F_B` remains **indefinite** on transverse modes (`λ_min=−g|B|<0`). A fundamental Higgs adds a **convex radial direction**, but **not** to the transverse gauge block where the negative mode lives. **The Savvidy/marginal wall is untouched by the Higgs.** **[T] (Savvidy sign) / [J] (that the Higgs does not remove it).**

**Acid-test verdict for (a).** The constructed theory is **free in the limit**; a free Proca field **has a mass**, so the mechanism does **not** fail on the free massless field in the required direction — it produces a gap for a free field. That is precisely the signature that the "gap" is the **explicit Higgs/Proca mass**, not the nonabelian gap. **[NG]** 

### §2B Physically closest: **center-stabilized small-`S¹`** (Ünsal–Yaffe) — fixed coupling, real nonabelian candidate, but zero rigorous legs

**Why one might prefer it.** At **fixed 't Hooft coupling** it is a genuine nonabelian theory with a real small parameter `e^{−S₀}`, and it **passes** the free-field acid test (the free/abelian field has no monopoles/bions, so `m_γ=0`; the resummation correctly fabricates no gap). It is the only class whose controlled regime is **not** a `g→0` free limit.

**Exact missing link (the whole route, once). [O]**
> Construct a reflection-positive continuum SU(N) measure on `ℝ³×S¹` at small `r` **with center stabilization**, prove a convergent cluster expansion over the nonabelian monopole–bion gas giving `Δ(r)≥δ` at small `r`, **and** prove `Δ(r)≥δ` uniformly as `r→∞`.

By the min–max lemma, the uniform statement **is** `Δ(ℝ⁴)≥δ` — the mass gap. Three independent deaths:
1. **Construction [O]:** no rigorous continuum SU(N) measure on `ℝ³×S¹`; the only rigorous "gas ⇒ area law" is **3D abelian** (Göpfert–Mack). **[R]**
2. **Different theory [NG]:** pure YM on a small periodic circle is **center-broken/deconfined** (GPY one-loop); the confined small-`r` phase requires a **deformation** (double-trace / adjoint fermions / twist) — *not* pure YM. Recovering pure YM is `(deformation→0)∘(r→∞)`, and the deformation→0 limit is an **infimum** (upper bound). **[T]/[NG]**
3. **Continuity [NG]:** adiabatic continuity is conjectural, **false** for pure YM (GPY transition), and unsupported in the closest analogue (four phases, Cossu–D'Elia). **[R]/[NG]**

**Acid-test verdict for (c).** Passes (free field correctly gapless) — but **sign-wrong**: the controlled end is the *more-gapped* end, so every proven statement is an **upper** bound on `Δ(ℝ⁴)`. Being non-vacuous is not being a floor. **[NG]**

---

## §3 The verdict (the decisive question, answered)

**Question.** Is there **any** 4D deformation that is (i) rigorously constructible, (ii) has a *proven* gap in the deformed theory, **and** (iii) a *proven* limit back to pure YM with the gap preserved?

**Answer: No — not in the enumerated classes.** The closest is §2A: it clears (i) only for a **Gaussian** limit, clears (ii) only for a **Proca/Higgs mass**, and fails (iii) outright — the constructed limit is **free at `g→0`**, and the paper itself leaves the non-Gaussian fixed-`g` construction open. Every genuinely nonabelian, acid-test-passing candidate (§2B and Class B generally) clears none of (i),(ii),(iii) rigorously and dies on the **direction** of the only free inequality.

**The single unifying obstruction (stated once, sharply). [T]/[J]**
> **The deformation → pure-YM limit is a min–max infimum over the deformation parameter, and every controlled reference sits on the more-gapped side of a monotone quantity (`β`, `r`, `N`, `m₀`, `g`). The only free inequality therefore runs *backward*: a gap-preserving limit requires a uniform floor at the physical, less-gapped end — which is `m_{0⁺⁺}>0`, i.e. H1 = (CP) = the mass gap. In deformation language: every constructible deformation's pure-YM limit crosses the Savvidy / strong-coupling crossover at fixed nonzero coupling, and "adiabatic continuity" *is* the unproven (CP).**

This is not three obstructions; it is one. "Min–max goes the wrong way," "adiabatic continuity is conjectural," and "crosses the crossover = H1" are the **same sentence** written on the `r`-axis, the deformation axis, and the `β`-axis. The charter's wall reappears here with a different subscript: the **deformation-axis twin of the transport wall**.

**What would flip it (the flip condition, stated honestly). [O]** A genuine path would be one of:
- a **non-Gaussian 4D construction at fixed nonzero `g`** whose gap is **not** an explicit deformation mass (kills §2A's "free Proca" failure), **together with** a proof that the gap survives undeforming at fixed `g`; **or**
- a **proven adiabatic continuity** with a **uniform-in-`r` lower bound** `Δ(r)≥δ` on the *nonabelian* small-`S¹` theory (kills §2B's direction failure) — but by min–max this *is* `Δ(ℝ⁴)≥δ`.

Both flip conditions are **H1/(CP) in disguise**: the missing object is the nonperturbative, uniform, sign-controlled floor at the physical end. There is no fourth option in the enumerated set. This is a **[NG] for the named classes**, not a global impossibility theorem (`∀class ≠ ∀possible`).

**Explicitly honest boundaries.**
- **[T] true statements we do bank:** the min–max lemma (finite-volume gap ≥ infinite-volume gap; `Δ(ℝ⁴)=inf_r Δ(r)`); the only unconditional 4D pure-YM gap is lattice strong coupling `β<1/48` (SZZ); the 3D YM–Higgs construction exists but is local-in-time with no measure/gap; Chatterjee's 4D YM–Higgs scaling limit exists but is massive Gaussian at `g→0`.
- **[R] external facts we do *not* overstate:** center stabilization and volume independence (Ünsal–Yaffe / Kovtun–Ünsal–Yaffe) are **physics-level**, and adiabatic continuity is **conjectural** (Poppitz review: "evidence"); Eguchi–Kawai is an `N=∞` reduction with no finite-`N` theorem; N=2 Seiberg–Witten is exact but its protected atoms are in the **charged** sector.
- **No QED.** The mass gap remains open. This document localizes the failure of construct-then-decouple; it does not close Clay.

---

## §4 Acid test (Task 4)

A correct construct-then-decouple must **fail on the free massless field** (free is gapless); the mechanism must be silent on the free theory and positive only on the nonabelian one.

| class | behavior on free massless | verdict |
|---|---|---|
| **A** (mass-deformed; Gaussian YM–Higgs/Proca; soft-SUSY) | Gives a **nonzero mass to a free/abelian field** (Proca / Higgs mass). The free field is supposed to be gapless; the mechanism is not silent. | **Fails acid test in the fatal direction** — the "gap" is the deformation's explicit mass, not the nonabelian gap. **[NG]** |
| **B** (center-stabilized small-`S¹`; Eguchi–Kawai) | Free/abelian field has **no monopoles/bions** ⇒ `m_γ=0`, correctly gapless; the KK spacing `1/r` is a finite-volume artifact, not a dynamical mass. | **Passes** (non-vacuous) — but is an **upper-bound machine**: the controlled end is the more-gapped end. **[NG]** |

**Where the nonabelian interaction enters — and why that is the death.** The interaction is *not* in the deformation mass (Class A) and *not* in the controlled small-parameter end (Class B). It enters **only** through the Savvidy/strong-coupling crossover at **fixed nonzero coupling**: the transverse Hessian `−D²+2 ad F_B` is indefinite (`λ_min=−g|B|<0`), and this is the single non-abelian sign behind the Weingarten/`6j` sign, the gluon condensate, and the residue `f₁²`. **No construction in the enumerated classes crosses that crossover.** A deformation that passes the acid test (Class B) points its inequality *backward* across it; a deformation that has a proven gap (Class A) puts the gap in *by hand* on the wrong side of it. **The acid test therefore does not select a route; it selects the wall.**

---

## §5 References (verification status 2026-11-24)

**Verified this session [R, verified] (title / authors / venue / arXiv from publisher or arXiv abstract pages):**
- A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Stochastic quantisation of Yang–Mills–Higgs in 3D*, *Invent. Math.* **237** (2024) 541–696, arXiv:2201.03487. — 3D only; local-in-time; no measure, no gap.
- S. Chatterjee, *A scaling limit of SU(2) lattice Yang–Mills–Higgs theory*, arXiv:2401.10507, *Prob. Math. Phys.* **7** (2026) 339–381. — massive **Gaussian** limit; `g=O(ε^{50d})→0`; non-Gaussian open.
- J. Magnen, V. Rivasseau, R. Sénéor, *Construction of YM₄ with an infrared cutoff*, *CMP* **155** (1993) 325–383. — UV-renormalised SU(2) Schwinger functions, **IR cutoff retained**; no OS positivity, no gap.
- H. Shen, R. Zhu, X. Zhu, *A stochastic analysis approach to lattice Yang–Mills at strong coupling*, *CMP* **400** (2023) 805–851, arXiv:2204.12737. — unconditional 4D lattice gap for `|β|<1/[16(d−1)]` (SU(N), `=1/48` in `d=4`); strong coupling only.
- M. Ünsal, L. G. Yaffe, *Center-stabilized Yang–Mills theory: confinement and large-N volume independence*, *PRD* **78** (2008) 065035, arXiv:0803.0344. — continuity **assumed**, not proven; no order parameter separating regimes.
- M. Shifman, M. Ünsal, *QCD-like theories on R³×S¹: a smooth journey from small to large r(S¹) with double-trace deformations*, *PRD* **78** (2008) 065004, arXiv:0802.1232. — double-trace deformation; continuity argued, not proven.
- P. Kovtun, M. Ünsal, L. G. Yaffe, *Volume independence in large N_c QCD-like gauge theories*, *JHEP* **06** (2007) 019, hep-th/0702021. — volume independence as orbifold equivalence; symmetry-realization conditions.
- G. Cossu, M. D'Elia, *Finite size phase transitions in QCD with adjoint fermions*, *JHEP* **0907** (2009) 048, arXiv:0904.1353. — **four phases**, two center-broken intermediates (no smooth continuity).
- T. Eguchi, H. Kawai, *Reduction of dynamical degrees of freedom in the large-N gauge theory*, *PRL* **48** (1982) 1063. — `N=∞` one-site reduction; no finite-`N` theorem.
- E. Fradkin, S. H. Shenker, *Phase diagrams of lattice gauge theories with Higgs fields*, *PRD* **19** (1979) 3682–3697. — fundamental Higgs: Higgs/confinement same phase; adjoint Higgs: a real boundary **can** exist.
- N. Seiberg, E. Witten, *Monopole condensation, and confinement in N=2 supersymmetric Yang–Mills theory*, *NPB* **426** (1994) 19–52, hep-th/9407087. — exact BPS spectrum; atoms in the charged sector.
- M. Aizenman, H. Duminil-Copin, *Marginal triviality of the scaling limits of critical 4D Ising and λφ₄⁴ models*, *Ann. Math.* **194** (2021) 163–235, arXiv:1912.07973. — the one solved marginal SPDE is **trivial**; no gauge template.
- E. Poppitz, *Notes on confinement on R³×S¹: from Yang–Mills, super-Yang–Mills, and QCD(adj) to QCD(F)*, arXiv:2111.10423 (review; related DOI 10.3390/sym14010180). — "evidence" against phase transitions, not proof.

**Standard / repo-verified (not re-checked this session) [R]:** Gross–Pisarski–Yaffe, *RMP* **53** (1981) 43 (GPY one-loop center breaking); C. Göpfert, G. Mack, *CMP* **82** (1982) 545 (rigorous 3D abelian monopole gas); G. 't Hooft, *NPB* **138** (1978) 1 & **153** (1979) 141 (twisted sectors); T. Sulejmanpasic, *PRL* **118** (2017) 011601, arXiv:1610.04009 (2D twisted-BC continuity; 2D-specific); Bhanot–Heller–Neuberger, *PLB* **113** (1982) 47 (EK center breaking); Bauerschmidt–Brydges–Slade, *LNM* **2242** (2019) (φ⁴₄ marginal template); Shen–Zhu–Zhu / Borgs–Seiler, *CMP* **91** (1983) 329.

**Correction log / citation hygiene.**
- "Davies et al." small-circle adiabatic continuity: **no such paper found**; the relevant lattice works are Cossu–D'Elia (2009) and Bonati–Cardinali–D'Elia (2018–2021). Flagged, not propagated (also flagged in `SMALLCIRCLE-ADIA §6`).
- The task's parenthetical "pure YM small-`S¹` is center symmetric" is **false** for pure YM (GPY: center-broken); only **deformed / twisted / QCD(adj)** is center-symmetric small-`r`. Corrected in the table.
- SZZ: the published bound is `|β|<1/[16(d−1)]`; the often-quoted `1/48` is its `d=4`, SU(N) specialisation.

**Internal cross-references.** `SMALLCIRCLE-ADIA.md` (§3 min–max; §5 deformed≠pure YM; §7 scorecard); `EGUCHI-KAWAI.md`; `CONSTRUCTION-METHODS.md` (a)–(f); `CONTINUUM-3D-TO-4D.md` (marginal vs Savvidy wall); `MONOTONE-COUPLING.md` (§2/§4 direction, deformation = infimum); `E2-FLIP-SECOND-SCALE.md` §(iii) (deformation scale [NG]); `BPS-SUSY-PROTECTION.md` (charge-sector orthogonality); `HIGGS-DIM-TRANS.md` (dilaton ≠ mechanism); `FORGE5-CONFINING-MONOPOLE.md` (rigor ⟺ wrong phase); `ABC-FINAL-LEDGER.md` §A (exact reduction ⇒ upper bound); `PROOF-P3-UNCONDITIONAL.md` / SZZ `β<1/48`; `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_min=−g|B|<0`); `PROOF-1FORM-ANOMALY.md` (`tr F²` blind to 1-form center).

---

## §6 Sharpest line

> **Construct-then-decouple is the right constructive instinct and it dies at exactly one place: the decoupling limit is a min–max infimum, and the gap-preserving infimum is the mass gap.** The single candidate with a real 4D construction (Chatterjee's YM–Higgs) buys its "gap" by sending `g→0`, so the gap is a free Proca mass — the nonabelian interaction is precisely what is decoupled away. The single genuinely nonabelian candidate (center-stabilized small-`S¹`) buys its control at the **more-gapped** end, so every theorem it could give is an **upper** bound. The interaction enters only at the Savvidy/strong-coupling crossover; **no deformation in this enumeration crosses it**, and *crossing it is the mass gap*. This is the deformation-axis twin of H1 — not a new door, the same wall with a new coordinate.

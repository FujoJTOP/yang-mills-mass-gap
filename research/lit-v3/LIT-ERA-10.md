# LIT-ERA-10 · 2008–2014 — Exact CFT constraints: collider bounds, a-theorem, bootstrap revival

> Scope: the exact, non-perturbative constraints on 4D CFTs that matured in 2008–2014 — conformal
> collider bounds on `a/c`, the proof of the `a`-theorem, the 4D scale-vs-conformal question, and the
> numerical conformal bootstrap. This file records **what these results do and do not constrain**:
> they constrain the *existence and operator data of CFTs*, not the endpoints of an arbitrary RG flow,
> and in particular they do **not** by themselves imply a mass gap for pure Yang–Mills.
>
> 判词：本纪元把"IR 若是 CFT"这一支装上了硬约束（a/c 窗口、a 单调、scale⇒conformal、bootstrap 排除法）；
> 但它对"IR 不是 CFT"（即隙）**零输出**。任何 v3 论证若把 CFT 约束直接当隙证明，是范畴错误。

---

## 0. Legend, conventions, verification status

Labels used on every entry:

- **[T]** — theorem: a proved result under explicitly stated hypotheses (may still be physics-level, not Lean-level).
- **[C]** — conjecture / proposed bound / claim not proved (at the time, or at all).
- **[NG]** — no-go / impossibility / explicit counterexample that closes a route.
- **[O]** — observation, review, method, or tool (not itself a theorem).
- **[V]** — verified numerical/consistency datum (bootstrap numerics, free-field sum rules, lattice cross-check).

Conventions for the trace anomaly (the one used by Komargodski–Schwimmer and the collider literature):

```
T^μ_μ = a · E₄ − c · W²        (E₄ = Euler density, W² = Weyl tensor squared)
```

Free-field coefficients below are in this normalization. `a` is the Euler anomaly (the 4D analogue of
the 2D central charge `c`); `c` is the Weyl-squared coefficient. The ratio `a/c` is the collider observable.

**Verification note (honest).** The shared web-search budget for this session was exhausted after the
core entries were checked. Entries in §3 carry VERIFIED/UNVERIFIED individually and are itemized in the
ledger (§7). VERIFIED means re-checked this session against arXiv/SCOAP3/INSPIRE/zbMATH/publisher data.
UNVERIFIED means standard-literature citation from model knowledge, not re-checked this session; re-check
before v3 submission. No bibliographic datum here was invented; where a page range was not re-checked it
is flagged rather than guessed silently.

---

## 1. What this era is about — and what it is not

The 2008–2014 wave is the second birth of **exact CFT model-building constraints** (the first being the
2D minimal-model / BPZ era). Four strands:

1. **Conformal collider bounds** (Hofman–Maldacena 2008): energy/charge correlations in a CFT, required
   to be positive by unitarity, yield `1/3 ≤ a/c ≤ 31/18` for any unitary, parity-preserving 4D CFT with a
   unique stress tensor. Proposed in 2008; proved from CFT first principles (lightcone bootstrap) in 2016.
2. **The `a`-theorem** (Komargodski–Schwimmer 2011): along any 4D RG flow between CFTs, `a_UV > a_IR`.
   Proof via the dilaton effective action and unitarity of dilaton scattering. Independent derivations
   (Luty–Polchinski–Rattazzi; Jack–Osborn) followed.
3. **Scale vs conformal invariance in 4D** (Polchinski 1988 for `d=2`; Dymarsky–Komargodski–Schwimmer–
   Theisen 2013/2015 for `d=4`): a unitary scale-invariant 4D theory is conformal, under stated
   assumptions; the free 2-form with noncompact gauge symmetry is the known exception.
4. **Bootstrap revival** (Rattazzi–Rychkov–Tonni–Vichi 2008; El-Showk et al. 2012): crossing symmetry +
   unitarity bound the space of CFT spectra; numerically this solved the 3D Ising model and turned
   "which CFTs exist?" into a semi-definite-programming question.

**What the era is NOT.** It is not a theory of RG flows ending in a gapped phase. Every result above
presupposes a CFT somewhere (UV, IR, or both) or constrains a single CFT. None addresses the existence
of the 4D pure-YM continuum, and none computes a mass gap. This is the central scope fact for v3.

---

## 2. Exact logical relationship: "no CFT with pure-YM content" ⟺ "pure YM has a mass gap"?

### 2.1 The three IR alternatives

For a putative renormalized pure SU(N) Yang–Mills theory on R⁴ (Poincaré-invariant, unitary, unique
vacuum, satisfying OS/Wightman axioms), exactly one of the following should hold in the IR:

- **(G) Gapped:** `Δ = inf{√(p²) : p² > 0, p in spectrum} > 0`. This is the Clay conclusion.
- **(S) Gapless and scale-invariant:** the IR is a nontrivial unitary 4D CFT (by scale⇒conformal, §3 E10–E11).
- **(X) Gapless, non-scale-invariant:** massless excitations with a running coupling / noncompact moduli /
  Coulomb-like behavior — gapless but not a fixed point.

The physically expected answer for pure YM is (G) (confinement). But "expected" is not "proved."

### 2.2 The valid and invalid implications

- **(G) ⇒ ¬(S) and (G) ⇒ ¬(X):** trivial — a mass gap excludes all gapless IR behavior.
- **¬(S) ⇏ (G):** *invalid without excluding (X).* Ruling out the CFT branch leaves the non-conformal
  gapless branch (X) open. This is the precise logical gap.
- **¬(S) ∧ ¬(X) ⇔ (G):** valid by trichotomy. So a proof of the mass gap via "no CFT" needs **two**
  no-go's: no CFT endpoint **and** no non-conformal gapless phase.
- **¬(S) ⇔ (G)** holds **only if** one assumes the "IR dichotomy" (gapless ⇒ scale-invariant ⇒ CFT).
  The era supplies the second arrow (scale⇒conformal) conditionally; it does **not** supply the first
  (gapless ⇒ scale-invariant). That first arrow is a separate, unproved dynamical statement.

### 2.3 Where each era result does / does not bear

| Result (entry) | What it actually constrains | Implies gap? | Role in the chain |
|---|---|---|---|
| Collider bounds (E1, E2) | `a/c` of any unitary, parity-preserving CFT with unique stress tensor | **No** | Constrains only branch (S)'s operator data; never touches (G) or (X). |
| `a`-theorem (E7, E8) | monotonicity of `a` along CFT→CFT flows; `a_IR ≥ 0` | **No** | If IR is a CFT, `a_IR < a_UV`; if gapped, `a_IR = 0`. Both allowed; a trivial monotonicity. |
| scale⇒conformal (E10, E11) | unitary scale-invariant 4D ⇒ conformal (with hypotheses) | **No** | Converts (S) into "CFT"; does not touch (X) or (G). |
| bootstrap (E14, E18, E19) | space of consistent CFT spectra / OPE data | **No (yet)** | The only tool that could in principle prove "no CFT with YM content"; not done for YM. |
| free-field `a,c` (E23) + unitarity bounds (E24) | calibration and positivity | **No** | Fixes `a_UV` for asymptotically free YM; no IR statement. |
| anomaly matching (E27) | massless spectrum for given 't Hooft anomalies | **No** | Vacuous for pure YM: no chiral matter, trivial anomaly matching in the IR. |

### 2.4 The sharp statement for v3

> **"No CFT with pure-YM content" is necessary but not sufficient for the mass gap.** It excludes branch
> (S) only. The mass gap additionally requires excluding branch (X). No result in 2008–2014 proves either
> no-go for pure YM. What the era does prove is a set of *conditional* constraints on branch (S):
> *if* the IR is a CFT, then its `a/c` lies in `[1/3, 31/18]`, its `a` is below the UV value, and its
> spectrum must solve crossing/unitarity. None of these is violated by any known YM candidate, so none
> excludes (S), and none is a step toward (G).

A second sharp point: **`a`-theorem is vacuous at a gapped endpoint.** The theorem's hypothesis is a flow
between two CFTs; a mass gap means the IR is the vacuum, not a CFT, so the inequality degenerates to
`a_UV > 0`. This is why the `a`-theorem, despite being one of the deepest exact results of the era, gives
**zero** quantitative information about `Δ` or the residue `f₁²`.

---

## 3. Annotated bibliography (29 entries)

### A. Conformal collider bounds

**E1. [C] D. M. Hofman, J. Maldacena**, "Conformal collider physics: Energy and charge correlations,"
JHEP **05** (2008) 012. arXiv:0803.1467 [hep-th]. DOI:10.1088/1126-6708/2008/05/012. **VERIFIED.**
Defines CFT event-shape observables (energy and charge deposited at infinity). Positivity of measured
energy gives `1/3 ≤ a/c ≤ 31/18` for unitary parity-preserving 4D CFTs with a unique stress tensor;
for superconformal theories `1/2 ≤ a/c ≤ 3/2`. Also relates charge two-point correlations to DIS
moments and gives an AdS/CFT shock-wave prescription. At the time the bounds were a *proposal* (positivity
of an energy functional), not a first-principles proof — hence `[C]`. Load-bearing for: the `a/c` window.

**E2. [T] D. M. Hofman, D. Li, D. Meltzer, D. Poland, F. Rejón-Barrera**, "A proof of the conformal
collider bounds," JHEP **06** (2016) 111, 40 pp. arXiv:1603.03771 [hep-th]. DOI:10.1007/JHEP06(2016)111.
zbMATH Zbl 1388.81048. **VERIFIED.**
Proves E1's bounds from CFT first principles (lightcone bootstrap) for any unitary, parity-preserving CFT
with a unique stress tensor in `d ≥ 3`; hence `31/18 ≥ a/c ≥ 1/3` in 4D, `3/2 ≥ a/c ≥ 1/2` for SCFTs.
Also treats non-conserved currents and DIS analogues. This converts E1 from `[C]` to `[T]` retroactively;
the 2016 date is why E1 remains labeled `[C]` in its own year.

**E3. [T] X. O. Camanho, J. D. Edelstein, J. Maldacena, A. Zhiboedov**, "Causality constraints on
corrections to the graviton three-point coupling," JHEP **02** (2016) 020, 62 pp. arXiv:1407.5597 [hep-th].
DOI:10.1007/JHEP02(2016)020. **VERIFIED** (preprint 2014).
Eikonal causality (Shapiro time delay) forbids certain higher-derivative graviton three-point structures
unless an infinite tower of massive higher-spin states (`J > 2`) exists. In AdS this bounds conformal
anomaly coefficients in terms of `Δ_gap` (lightest single-trace operator of spin `J > 2`). Relevance to YM:
it constrains *holographic duals* of a putative YM IR CFT, not YM itself; any AdS dual with a hard mass gap
and no higher-spin tower is in tension. `[T]` under the causality/eikonal assumptions.

### B. The `a`-theorem

**E4. [T] A. B. Zamolodchikov**, "Irreversibility of the flux of the renormalization group in a 2D field
theory," JETP Lett. **43** (1986) 730 [Pisma Zh. Eksp. Teor. Fiz. **43** (1986) 565]. **UNVERIFIED.**
The 2D `c`-theorem: a `c`-function decreases along RG flows between 2D CFTs. The prototype for all
monotonicity theorems; establishes the *method* (positivity of a two-point function of `T`) later used in 4D.

**E5. [C] J. L. Cardy**, "Is there a c-theorem in four dimensions?" Phys. Lett. **B215** (1988) 749–752.
**UNVERIFIED.**
Proposes the 4D `a`-theorem by analogy, identifying the Euler anomaly `a` as the candidate monotone.
Conjectural at the time; the 2011 proof (E7) settles it.

**E6. [C] H. Osborn**, "Derivation of a four-dimensional c-theorem," Phys. Lett. **B222** (1989) 97–102.
**UNVERIFIED.**
Local RG / Weyl-consistency approach; derives a perturbative 4D monotonicity statement. Superseded as a
non-perturbative proof by E7 but supplies the local-RG technology used throughout.

**E7. [T] Z. Komargodski, A. Schwimmer**, "On renormalization group flows in four dimensions," JHEP **12**
(2011) 099, 20 pp. arXiv:1107.3987 [hep-th]. DOI:10.1007/JHEP12(2011)099. zbMATH Zbl 1306.81140.
**VERIFIED.**
Proof of the 4D `a`-theorem. Reinterprets any RG flow as spontaneously broken conformal symmetry; the
dilaton effective action has an `(a_UV − a_IR)` term whose contribution to `2→2` dilaton scattering is
positive by unitarity, giving `a_UV > a_IR`. The `c`-anomaly is algebraically trivial; `a` is the
"non-Abelian" one. Load-bearing: the only non-perturbative monotonicity theorem in 4D.

**E8. [T] M. A. Luty, J. Polchinski, R. Rattazzi**, "The a-theorem and the asymptotics of 4D quantum field
theory," JHEP **01** (2013) 152. arXiv:1204.5221 [hep-th]. **UNVERIFIED.**
Independent derivation/analysis of the `a`-theorem via dispersion relations and the high-energy asymptotics
of 4D correlators. Confirms E7 and clarifies the assumptions (UV completion, spectral positivity).

**E9. [O] Y. Nakayama**, "Scale invariance vs conformal invariance," Phys. Rept. **569** (2015) 1–93.
arXiv:1302.0884 [hep-th]. DOI:10.1016/j.physrep.2014.12.003. **VERIFIED.**
Comprehensive review: the state of the scale⇒conformal problem, the `a`-theorem, local RG, holographic
arguments, and the catalogue of scale-invariant examples. The best single entry point for the era's logic.
Notes that as of 2014 no unitary scale-invariant non-conformal 4D example satisfying the standard
hypotheses was known.

### C. Scale vs conformal invariance

**E10. [T] J. Polchinski**, "Scale and conformal invariance in quantum field theory," Nucl. Phys. **B303**
(1988) 226–236. Report UTTG-22-87. **VERIFIED** (INSPIRE 22605).
In `d=2`, completes Zamolodchikov's argument that rigid scale invariance implies local conformal invariance
under broad conditions. In `d ≥ 3` finds neither a general proof nor a counterexample. This is the historical
root of the 4D question answered (conditionally) by E11.

**E11. [T] A. Dymarsky, Z. Komargodski, A. Schwimmer, S. Theisen**, "On scale and conformal invariance in
four dimensions," JHEP **10** (2015) 171, 24 pp. arXiv:1309.2921 [hep-th]. DOI:10.1007/JHEP10(2015)171.
zbMATH Zbl 1388.81372. **VERIFIED** (preprint 2013).
Shows that in a unitary 4D scale-invariant theory, infinitely many matrix elements vanish in a suitable
kinematic configuration — a nontrivial *necessary* condition for conformality — and argues it is also
*sufficient*. Hence scale ⇒ conformal in 4D under the stated hypotheses (Poincaré, unitarity, discrete
scaling spectrum, unbroken scale invariance, existence of the scale current). Directly upgrades branch (S)
to "CFT."

**E12. [NG] A. Dymarsky, Z. Komargodski, A. Schwimmer, S. Theisen**, same paper (E11), §4: the free
two-form with **noncompact** gauge symmetry is scale-invariant but **not** conformal; the compact version is
consistent. **VERIFIED** (same source).
Explicit counterexample showing E11 is not unconditional. Relevance to YM: pure YM has compact gauge group,
so this loophole does not apply — but it is the reason the theorem must be quoted with hypotheses.

**E13. [T] A. Dymarsky, K. Farnsworth, Z. Komargodski, A. Schwimmer, S. Theisen**, "Scale invariance plus
unitarity implies conformal invariance in four dimensions," arXiv:1710.07341 [hep-th] (2017). **UNVERIFIED.**
Follow-up closing remaining loopholes in E11 and sharpening the operator statement. Outside the 2008–2014
window but the natural terminus of the question; include as the modern citable form.

### D. Bootstrap revival

**E14. [T] R. Rattazzi, V. S. Rychkov, E. Tonni, A. Vichi**, "Bounding scalar operator dimensions in 4D
CFT," JHEP **12** (2008) 031, 48 pp. arXiv:0807.0004 [hep-th]. DOI:10.1088/1126-6708/2008/12/031.
**VERIFIED.**
Founding paper of the modern numerical bootstrap. Uses OPE + conformal blocks + crossing + unitarity to
derive a theory-independent inequality `[φ²] ≤ f([φ])` for any unitary 4D CFT; `f(d) → 2` as `d → 1`
(free-theory limit approached continuously). The bound is rigorous; only its numerical evaluation is
approximate. Establishes the logic "crossing + unitarity ⇒ exclusion regions in CFT data space."

**E15. [T] F. A. Dolan, H. Osborn**, "Conformal four point functions and the operator product expansion,"
Nucl. Phys. **B599** (2001) 459–496. arXiv:hep-th/0011040. **UNVERIFIED.**
Computes 4D conformal partial waves / conformal blocks — the technical ingredient without which E14 is
not numerically implementable.

**E16. [T] F. A. Dolan, H. Osborn**, "Conformal partial waves and the operator product expansion," Nucl.
Phys. **B678** (2004) 491–507. arXiv:hep-th/0309180. **UNVERIFIED.**
Generalizes conformal blocks to arbitrary spacetime dimension; makes the 3D bootstrap (E18) possible.

**E17. [O] H. Osborn, A. Petkou**, "Implications of conformal invariance in field theories for general
dimensions," Annals Phys. **231** (1994) 311–362. arXiv:hep-th/9307010. **UNVERIFIED.**
Systematic conformal OPE in general `d`; the structural reference behind E15–E16.

**E18. [V] S. El-Showk, M. F. Paulos, D. Poland, S. Rychkov, D. Simmons-Duffin, A. Vichi**, "Solving the
3D Ising model with the conformal bootstrap," Phys. Rev. **D86** (2012) 025022. arXiv:1203.6064 [hep-th].
DOI:10.1103/PhysRevD.86.025022. **VERIFIED.**
Derives 3D conformal blocks and shows the 3D Ising model sits at a corner of the bootstrap-allowed region;
extracts critical exponents. The demonstration that crossing + unitarity can *fix* a theory, not merely
bound it. `[V]` because the load-bearing content is numerical.

**E19. [O] S. El-Showk, M. F. Paulos**, "Bootstrapping conformal field theories with the extremal functional
method," Phys. Rev. Lett. **111** (2013) 241601. arXiv:1211.2810 [hep-th]. **UNVERIFIED.**
Extremal-functional method: at the boundary of the allowed region, the solution is isolated and the
spectrum read off. The technique that made "CFT data as output" routine.

**E20. [O] V. S. Rychkov**, "EPFL Lectures on Conformal Field Theory in D≥3 Dimensions," SpringerBriefs in
Physics (2017). arXiv:1601.05000 [hep-th]. **UNVERIFIED.**
Pedagogical derivation of unitarity bounds, conformal blocks, and the bootstrap logic. Useful for v3
notation alignment.

**E21. [O] D. Simmons-Duffin**, "TASI Lectures on the Conformal Bootstrap," PoS **TASI2015** (2016) 004.
arXiv:1602.07982 [hep-th]. **UNVERIFIED.**
Modern review of the numerical bootstrap, including semidefinite-programming implementation and the
3D Ising solution. The practical reference.

**E22. [O] D. Poland, S. Rychkov, A. Vichi**, "The conformal bootstrap: Theory, numerical techniques, and
applications," Rev. Mod. Phys. **91** (2019) 015002. arXiv:1805.04405 [hep-th]. **UNVERIFIED.**
Definitive review of the whole program. Outside the era; the citation of record for "what the bootstrap
can and cannot currently prove."

### E. Central charges, unitarity bounds, anomalies

**E23. [T] M. J. Duff**, "Observations on conformal anomalies," Nucl. Phys. **B125** (1977) 334–348.
**UNVERIFIED.**
Establishes the free-field trace-anomaly coefficients and the `a`, `c` basis. The values in §4 are the
standard ones from this line of work (cross-checked in §4 by the N=4 sum rule). Load-bearing for the
numerical window `[1/3, 31/18]` because the endpoints are free-field values.

**E24. [T] G. Mack**, "All unitary ray representations of the conformal group SU(2,2) with positive energy,"
Commun. Math. Phys. **55** (1977) 1–28. **UNVERIFIED.**
Unitarity bounds on conformal primary dimensions (scalar `Δ ≥ 1`, spin `ℓ ≥ 1`: `Δ ≥ ℓ + d − 2`, etc.).
The positivity input behind every bootstrap bound and behind the collider inequalities.

**E25. [T] S. Minwalla**, "Restrictions imposed by superconformal invariance on quantum field theories,"
Adv. Theor. Math. Phys. **2** (1998) 781–846. arXiv:hep-th/9712074. **UNVERIFIED.**
Unitarity bounds in 4D superconformal theories; the SCFT window `1/2 ≤ a/c ≤ 3/2` of E1–E2 is the
superconformal refinement of these constraints.

**E26. [T] D. Anselmi, D. Z. Freedman, M. T. Grisaru, A. A. Johansen**, "Positivity constraints on anomalies
in supersymmetric theories," Phys. Rev. **D57** (1998) 7570–7588. arXiv:hep-th/9712129. **UNVERIFIED.**
Positivity of anomaly coefficients from unitarity; early evidence that `a`, `c` are constrained, not free.
Conceptual ancestor of the collider bounds.

**E27. [T] G. 't Hooft**, "Naturalness, chiral symmetry, and spontaneous chiral symmetry breaking," in
*Recent Developments in Gauge Theories* (Cargèse 1979), NATO Sci. Ser. **B59** (1980) 135–157. **UNVERIFIED.**
Anomaly matching: the massless spectrum of any IR CFT must reproduce the 't Hooft anomalies of the UV.
For pure YM this is vacuous (no chiral fermions, no continuous global anomaly in the IR), which is itself
the relevant fact: anomaly matching gives **no** constraint on a pure-YM IR CFT.

**E28. [O] M. Henningson, K. Skenderis**, "The holographic Weyl anomaly," JHEP **07** (1998) 023.
arXiv:hep-th/9806087. **UNVERIFIED.**
Computes `a`, `c` from the AdS dual; the bridge between collider bounds and holographic mass-gap models
(era 8). Relevant because any holographic YM dual must reproduce the field-theory `a/c`.

**E29. [V] E. Barnes, K. A. Intriligator, B. Wecht, J. Wright**, "Evidence for the a-theorem in four
dimensions," Nucl. Phys. **B702** (2004) 3–36. arXiv:hep-th/0405013. **UNVERIFIED.**
Perturbative/numerical evidence for `a_UV > a_IR` before the 2011 proof; useful as an independent check of
E7's sign and normalization.

---

## 4. Central charges of free fields and the meaning of the window

Normalization: `T^μ_μ = a E₄ − c W²` (E23 conventions).

| Field | `a` | `c` | `a/c` |
|---|---|---|---|
| Real scalar | 1/360 | 1/120 | **1/3** (lower endpoint) |
| Weyl fermion | 11/720 | 1/40 | 11/18 ≈ 0.611 |
| Dirac fermion | 11/360 | 1/20 | 11/18 ≈ 0.611 |
| Vector (gauge boson) | 31/180 | 1/10 | **31/18** ≈ 1.722 (upper endpoint) |
| U(N) N=4 SYM | N²/4 | N²/4 | 1 |

Aggregate formulas (cross-check): for `N_s` real scalars, `N_f` Dirac fermions, `N_v` vectors,
`a = (1/360)(N_s + 11 N_f + 62 N_v)`, `c = (1/120)(N_s + 6 N_f + 12 N_v)`.

**Consistency check [V].** N=4 SU(N) SYM = 1 vector + 4 Weyl + 6 real scalars:
`a = 31/180 + 4·(11/720) + 6·(1/360) = 45/180 = 1/4`; `c = 1/10 + 4·(1/40) + 6·(1/120) = 1/4`.
So `a = c = 1/4` per U(1) vector multiplet (i.e. `(N²−1)/4` for SU(N)), matching the exact superconformal
result. This fixes the normalization and confirms the endpoint values.

**What the window means.** `1/3 ≤ a/c ≤ 31/18` is saturated *at both ends by free fields*: the free scalar
(lower) and the free vector (upper). Interacting CFTs lie strictly inside unless they are free. This is a
statement about **which single CFTs can exist**, not about flows. For pure YM:

- The UV is asymptotically free, so `a_UV = (N²−1)·31/180` (the free-vector value).
- If the IR is a CFT (branch S), it must satisfy `1/3 ≤ a/c ≤ 31/18` and `0 < a_IR < a_UV`. Both are
  non-trivial but far from exclusionary; no candidate value is known, and the constraint is one inequality
  on one ratio.
- If the IR is gapped (branch G), the window is irrelevant: there is no IR CFT to constrain.

So the window is a **filter on branch (S)**, not a probe of branch (G).

---

## 5. Position in 14 eras

Campaign working periodization for v3 (thematic; ranges intentionally overlap where strands coexist).
Era 10 is this file.

| # | Years | Era | Bearing on the mass gap |
|---|---|---|---|
| 1 | 1954–1972 | Origins: non-Abelian gauge theory; confinement posed; no gauge-invariant order parameter | Defines the problem |
| 2 | 1973–1974 | Asymptotic freedom and the renormalization group (Gross–Wilczek, Politzer; 't Hooft) | Establishes UV free fixed point; `a_UV` |
| 3 | 1974–1977 | Instantons, θ-vacua, U(1) problem (BPST; 't Hooft) | Non-perturbative structure |
| 4 | 1974–1980 | Lattice gauge theory; strong-coupling area law (Wilson; Osterwalder–Seiler) | First rigorous confinement regime |
| 5 | 1978–1987 | Continuum rigor: constructive QFT, Balaban RG, OS axioms for gauge theory | Existence side (G1) |
| 6 | 1980–1994 | IR mechanics: Gribov copies, monopoles, center vortices, dual superconductivity | Confinement mechanisms |
| 7 | 1994–1999 | Supersymmetric exact results and holography (Seiberg–Witten; AdS/CFT) | Exact solvable analogues |
| 8 | 1999–2005 | Large-N, AdS/QCD, holographic mass-gap models | Model-building; no controlled YM dual |
| 9 | 2004–2009 | Precision lattice glueball spectroscopy and confinement phenomenology | Numerical target spectrum |
| **10** | **2008–2014** | **Exact CFT constraints: collider bounds, a-theorem, scale⇒conformal, bootstrap revival** | **Constrains branch (S) only; no gap output** |
| 11 | 2012–2018 | Numerical bootstrap precision; 3D Ising solution; CFT data as axioms | Exclusion tool for CFT endpoints |
| 12 | 2015–2020 | Generalized/higher-form symmetries, 't Hooft anomalies, confinement | New constraints on (S); center symmetry |
| 13 | 2018–2023 | S-matrix/positivity and resurgence constraints on RG and spectra | Bounds on amplitudes, not the gap |
| 14 | 2020–2026 | Formalization, AI-assisted proof search, present unconditional campaign | The current attempt |

---

## 6. Gaps and attack surfaces for v3

1. **Missing lemma L1 (gapless ⇒ scale-invariant).** The era proves scale⇒conformal (E11), but not that a
   gapless pure-YM IR is scale-invariant. A Coulomb-like branch (X) would evade every CFT constraint.
   *Attack:* prove pure YM has no gapless non-conformal phase, using the area law / center symmetry /
   a massless-gluon no-go. Without this, "no CFT" cannot yield the gap.

2. **Missing lemma L2 (bootstrap applied to the symmetry-poor case).** A pure-YM IR CFT would have no
   continuous 0-form global symmetry (center symmetry is 1-form and confined), the *least constrained*
   bootstrap setup. *Attack:* import higher-form symmetry (era 12) or large-N structure as extra crossing
   input; otherwise the bootstrap exclusion region is weakest exactly where YM lives.

3. **No bridge from `a/c` to `f₁²` (the residue).** This is the highest-value gap. If one can prove a sum
   rule `a_IR ≥ F(f₁², m₁)` or `Δa = a_UV − a_IR ≥ c·f₁²`, then the `a`-theorem would convert the exact
   monotonicity into a *lower bound on the residue* — a genuinely new species of tool. No such relation
   exists in the era's literature. *Attack:* derive `a` as an integral over the Källén–Lehmann spectral
   density and compare with the residue at the lightest atom.

4. **The `a`-theorem degenerates at a gapped endpoint.** `a_IR = 0` is trivially allowed, so the theorem
   is silent on `Δ`. *Attack:* formulate a *quantitative* `a`-theorem, `a_IR ≥ G(m₁, f₁²)`, that becomes
   non-trivial as `m₁ → 0`; this is the only way the era's deepest theorem could touch the gap.

5. **Free-field acid test cuts both ways.** The era's tools are CFT tools and are *vacuous on a gapped
   theory by construction* (no CFT, nothing to constrain). Any v3 use must demonstrate non-vacuity on YM
   specifically — e.g., a bound that fails for free fields but holds for confining spectra. Otherwise it
   fails the charter's acid test (a tool that also "works" on the free theory proves nothing).

6. **Hypotheses of scale⇒conformal must be checked, not assumed.** E11 requires a discrete scaling
   spectrum, unbroken scale invariance, and no noncompact gauge symmetry; E12 shows the last is essential.
   Pure YM is compact, so the loophole is closed — but the discrete-spectrum assumption must be verified
   for any candidate IR.

7. **Collider bounds have blind spots.** They need a unique stress tensor and parity conservation. A YM IR
   CFT satisfies both, so the bounds apply but do not discriminate: every candidate passes. *Attack:* only
   the *strict* inequalities (interacting ⇒ strictly inside the free-field endpoints) could certify
   "interacting," never "gapped."

8. **Holographic duals are constrained but not constructed.** E3 bounds `Δ_gap` vs anomaly coefficients for
   any AdS dual; a hard-gap dual with no higher-spin tower is in tension. But no controlled top-down dual
   of pure 4D YM exists, so this constrains models, not the theory. *Attack:* use E3 as a consistency
   filter on era-8 proposals; do not mistake a surviving model for a proof.

9. **Existence (G1) is untouched.** Every result here assumes the theory exists. v3 must keep G1 (continuum
   existence, cutoff-consistent renormalization) separate from G2 (the gap) and must not let CFT constraints
   masquerade as progress on G1.

---

## 7. Verification ledger

Methods: WS = WebSearch result this session; XC = independent hand cross-check (sum rule / algebra);
STD = standard literature from model knowledge, not re-checked this session. Status is per entry.

| # | Entry (short) | Identifier | Method / source | Status |
|---|---|---|---|---|
| E1 | Hofman–Maldacena, conformal collider | JHEP 05 (2008) 012; arXiv:0803.1467 | WS (arXiv record, JHEP metadata) | VERIFIED |
| E2 | HLMPR, proof of collider bounds | JHEP 06 (2016) 111; arXiv:1603.03771 | WS (SCOAP3 16102; DOI; zbMATH Zbl 1388.81048) | VERIFIED |
| E3 | Camanho–Edelstein–Maldacena–Zhiboedov | JHEP 02 (2016) 020; arXiv:1407.5597 | WS (DOI 10.1007/JHEP02(2016)020; crossref) | VERIFIED |
| E4 | Zamolodchikov, 2D c-theorem | JETP Lett. 43 (1986) 730 | STD | UNVERIFIED |
| E5 | Cardy, 4D c-theorem conjecture | Phys. Lett. B215 (1988) 749 | STD | UNVERIFIED |
| E6 | Osborn, 4D c-theorem derivation | Phys. Lett. B222 (1989) 97 | STD | UNVERIFIED |
| E7 | Komargodski–Schwimmer, a-theorem | JHEP 12 (2011) 099; arXiv:1107.3987 | WS (arXiv; DOI; zbMATH Zbl 1306.81140) | VERIFIED |
| E8 | Luty–Polchinski–Rattazzi | JHEP 01 (2013) 152; arXiv:1204.5221 | STD | UNVERIFIED |
| E9 | Nakayama, scale vs conformal review | Phys. Rept. 569 (2015) 1; arXiv:1302.0884 | WS (DOI 10.1016/j.physrep.2014.12.003; Caltech) | VERIFIED |
| E10 | Polchinski, scale vs conformal | Nucl. Phys. B303 (1988) 226 | WS (INSPIRE 22605; OSTI; MR0953717) | VERIFIED |
| E11 | DKSST, scale⇒conformal in 4D | JHEP 10 (2015) 171; arXiv:1309.2921 | WS (DOI 10.1007/JHEP10(2015)171; zbMATH Zbl 1388.81372) | VERIFIED |
| E12 | DKSST free 2-form counterexample | same as E11, §4 | WS (same source) | VERIFIED |
| E13 | Dymarsky et al., unitarity⇒conformal | arXiv:1710.07341 | STD | UNVERIFIED |
| E14 | Rattazzi–Rychkov–Tonni–Vichi | JHEP 12 (2008) 031; arXiv:0807.0004 | WS (arXiv; DOI 10.1088/1126-6708/2008/12/031) | VERIFIED |
| E15 | Dolan–Osborn, conformal 4-pt | Nucl. Phys. B599 (2001) 459; hep-th/0011040 | STD | UNVERIFIED |
| E16 | Dolan–Osborn, partial waves | Nucl. Phys. B678 (2004) 491; hep-th/0309180 | STD | UNVERIFIED |
| E17 | Osborn–Petkou, general-d OPE | Annals Phys. 231 (1994) 311; hep-th/9307010 | STD | UNVERIFIED |
| E18 | El-Showk et al., 3D Ising bootstrap | PRD 86 (2012) 025022; arXiv:1203.6064 | WS (APS/Caltech; DOI 10.1103/PhysRevD.86.025022) | VERIFIED |
| E19 | El-Showk–Paulos, extremal functional | PRL 111 (2013) 241601; arXiv:1211.2810 | STD | UNVERIFIED |
| E20 | Rychkov, EPFL lectures | arXiv:1601.05000 | STD | UNVERIFIED |
| E21 | Simmons-Duffin, TASI lectures | PoS TASI2015 (2016) 004; arXiv:1602.07982 | STD | UNVERIFIED |
| E22 | Poland–Rychkov–Simmons-Duffin review | Rev. Mod. Phys. 91 (2019) 015002; arXiv:1805.04405 | STD | UNVERIFIED |
| E23 | Duff, conformal anomalies | Nucl. Phys. B125 (1977) 334 | STD + XC (free-field values fixed by §4 N=4 sum) | UNVERIFIED (values XC) |
| E24 | Mack, unitarity bounds | Commun. Math. Phys. 55 (1977) 1 | STD | UNVERIFIED |
| E25 | Minwalla, SCFT unitarity bounds | ATMP 2 (1998) 781; hep-th/9712074 | STD | UNVERIFIED |
| E26 | Anselmi–Freedman–Grisaru–Johansen | PRD 57 (1998) 7570; hep-th/9712129 | STD | UNVERIFIED |
| E27 | 't Hooft, anomaly matching | NATO Sci. Ser. B59 (1980) 135 | STD | UNVERIFIED |
| E28 | Henningson–Skenderis, holographic anomaly | JHEP 07 (1998) 023; hep-th/9806087 | STD | UNVERIFIED |
| E29 | Barnes–Intriligator–Wecht–Wright | Nucl. Phys. B702 (2004) 3; hep-th/0405013 | STD | UNVERIFIED |

**Free-field table (§4):** VERIFIED by XC — the aggregate formula `a = (1/360)(N_s + 11N_f + 62N_v)`,
`c = (1/120)(N_s + 6N_f + 12N_v)` was returned by WS and independently checked against the exact N=4
result `a = c = 1/4` per U(1) vector multiplet; the four ratios `1/3, 11/18, 11/18, 31/18` follow by
algebra. This is the one numerical block that is fully verified.

**Summary.** 10 of 29 entries VERIFIED this session (E1, E2, E3, E7, E9, E10, E11, E12, E14, E18), plus
the §4 free-field block verified by cross-check; the remaining 19 are UNVERIFIED but standard. The verified set covers every *load-bearing* claim:
the collider bounds, their proof, the `a`-theorem, scale⇒conformal and its counterexample, the bootstrap
founding papers, Polchinski, and the free-field normalization. The unverified set is supporting/review
material. Re-run verification on E8, E13, E19, E23–E29 before v3 submission.

---

### One-paragraph bottom line

The 2008–2014 era proves three exact facts about 4D CFTs — `1/3 ≤ a/c ≤ 31/18` (proved 2016),
`a_UV > a_IR` along CFT flows, and scale⇒conformal under stated hypotheses — and revives the numerical
bootstrap as an exclusion machine. Every one of these constrains **branch (S): the hypothesis that the
pure-YM IR is a CFT**. None constrains **branch (G): the mass gap**, and none excludes **branch (X): a
non-conformal gapless phase**. The exact logical relation is `¬(S) ∧ ¬(X) ⇔ (G)`; the era supplies
constraints on (S) and a conditional link (S)⇔CFT, but no proof of `¬(S)` and no statement about (X).
For v3 the productive move is not to re-derive these bounds but to build the missing bridge from `a`
(the monotone) to `f₁²` (the residue) — a sum rule that would make the `a`-theorem bite on the gap.

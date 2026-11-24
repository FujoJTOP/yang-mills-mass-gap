# INVENT-SYNTHESIS — 16 routes (8 existing theory + 8 bare-handed species): the renormalon is AVOIDABLE by changing the observable (two renormalon-free reformulations), but every route lands on the same strut

> **VERDICT (no QED).** Sixteen routes were dispatched to build the missing strut (a renormalon-free,
> phase-sensitive, non-order-parameter characterization of the `0⁺⁺` anomaly). **None lays the roof.** But
> two genuinely new structural facts emerged. **(1) The `u=2` renormalon is avoidable by changing the
> observable**: two independent routes — the **large-N master field** (`m_{0⁺⁺}(λ)` as a spectral pole)
> and **dual superconductivity** (the dual-Higgs mass from a dim-2 condensate) — write the strut *without*
> the `u=2` pole, moving the obstruction from "renormalon ambiguity" to "**existence / condensation**" (a
> different, plausibly more attackable problem). **(2) The whole known landscape is exhausted with a
> multi-dimensional no-go atlas**, each species dying of its own exact theorem. One bare-handed species
> (**level-set percolation on `tr F²`**) survives as a genuine `P1`-type object, blocked at a single
> inequality. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric · **[V-here]** new this round. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The 16-route scorecard

### Wave A — along EXISTING theory

| route (file) | verdict | what it gives / exact step |
|---|---|---|
| **Large-N master** `ROOF-LARGEN-MASTER` | `[T]` renormalon-free + `1/N`; `[NG]` | **first** strut with NO `u=2` renormalon AND a small parameter (`m/√σ` N-independent: 3.79→3.30, SU3 +7.7%); blocked at **4D master-field existence / closing the loop equation** (Singer did 2D only) |
| **Dual superconductor** `ROOF-DUALSUPERCONDUCTOR` | `[T]` structural; `[NG]` | **renormalon-free by construction** (dual Higgs mass, dim-2 `v` ⇒ pole at `u=1` not `u=2`); cost: **monopole condensation = confinement = H1**; gap ≠ KL atom (fails `:φ²:`); ratio short |
| **Resurgence** `ROOF-RESURGENCE` | `[T]`/`[NG]` | ambiguity = **Stokes constant (computable)**; condensate = **free transseries parameter** (no perturbation fixes it); rewrites the wall, doesn't cross it |
| **OPE** `ROOF-OPE` | `[NG]` | `u=2` pole; Ward transports not removes; positivity ambiguous; **only `ρ(s)` renormalon-free**, but bounding the atom from it = `H1` (= bootstrap `Δ_min=0 ∀K`) |
| **Holography** `ROOF-HOLOGRAPHY` | `[C]`/`[NG]` | definite ratios (`m_{0⁺⁺}z₀=j_{2,1}`, `0⁺⁺*/0⁺⁺=1.64`), but **ambiguity relocated into IR boundary** (`z₀,c,R`); assumes, not derives, `H1` |
| **Stochastic quantization** `ROOF-STOCHASTIC` | `[T]`/`[NG]` | works `d≤3` (CCHS); **`d=4` critical**; counterterms = `u=2` renormalon; FP gap is a tautology |
| **Tensor-RG** `ROOF-TENSOR-RG` | `[NG]` | renormalon-free but **upper bounds only**; entanglement gap ≠ physical gap (acid-dirty); no 3+1D theorem; missing lower half = `H1` |
| **Defect CFT** `ROOF-DEFECT-CFT` | `[T]`/`[NG]` | protected **1-form** data (displacement `Δ_D=4`), but **wrong regime** (no IR CFT) + **channel** (1-form vs 0⁺⁺) |

### Wave B — BARE-HANDED new species

| species (file) | verdict | mechanism of death / yield |
|---|---|---|
| **Lee–Yang / Fisher zeros** `INVENT-LEEYANG` | type 4/5, content 0/5 | new TYPE exists (analytic-structure functional); content dies on **kinematic `β=0` zeros** (acid-dirty: fires on free field); exact step = **IR/kinematic zero separation = UV/IR = H1** |
| **Percolation / random-cluster** `INVENT-PERCOLATION` | `[NG]` + yield | genuine representation (FK identity `Q_p=⟨σ₀σ_x⟩`, FK–FKG); dies on channel + B2 + no-q-state; **yield: level-set percolation on `tr F²` = a real `P1` species**, missing "percolation ⇒ floored χ" |
| **Non-equilibrium housekeeping heat** `INVENT-HOUSEKEEPING` | `[T]`/`[NG]` | outside Theorem BE, but **reversibility dichotomy**: intrinsic ⇒ `σ≡0` (re-derives BE); driven ⇒ controller-dependent (`σ=2γ²/m`) + acid-dirty (max at `m→0`) |
| **Center-flux twisted gap** `INVENT-FLUX` | `[T]` | **center-neutrality**: `tr F²` is a singlet ⇒ `Γ(𝔷)=Γ(1)`; the species is the 't Hooft loop in spectral clothing |
| **Modular Hamiltonian** `INVENT-MODULAR` | `[NG]` | provably **one of four**; KMS ⇒ modular gap `=β·`(mass gap) (circular); acid fails on `:φ²:`; the non-spectral Connes type is `III₁` universally (gap-independent) |
| **Thermodynamic curvature** `INVENT-CURVATURE` | `[T]`/`[NG]` | **Theorem BE in Riemannian dress** (`g_hh=χ`); single-parameter ⇒ flat; Gaussian ⇒ `R≡−1` (acid-vacuous); not a `μ`-functional (needs an embedding choice) |
| **Gradient-flow response** `INVENT-FLOWRESP` | `[NG]` + yield | renormalon-free but **atom-blind** (`E(t)`=Laplace of spectrum; rate set by continuum edge); response is a **moment** (homogeneous); **yield: renormalon-free `0⁺⁺` scale measurement** |
| **Extreme-value tail dependence** `INVENT-TAILDEP` | `[NG]` | reduces to **κ relabelled** (`θ_EV` monotone in correlation = homogeneous); tail index floor-blind; atom-blind; **closes the EVT/large-deviation class** (3rd B2 casualty) |

**0 / 16 lay the roof. 16 / 16 converge on the strut.**

---

## §2 · The two recoveries: renormalon is AVOIDABLE by changing the observable

The decisive finding of the round. `ROOF-CONDENSATE-FLOOR` had proved the **gluon condensate** is
irreducibly renormalon-ambiguous (ambiguity `~Λ⁴` = its own size). But two routes **write the strut without
the `u=2` pole** by not flooring the condensate at all:

| route | the observable that replaces `⟨tr F²⟩` | why no `u=2` | new obstruction |
|---|---|---|---|
| **Large-N master** | the spectral pole `m_{0⁺⁺}(λ)` (a physical mass) | a mass ratio `m/√σ` has no renormalon | **4D master field** / close the Makeenko–Migdal hierarchy |
| **Dual superconductor** | the dual-Higgs curvature `m_H=√(2λ)v` | `v` is dim-2 ⇒ Borel pole at `u=1`, UV-finite on dual lattice | **monopole condensation** = confinement = `H1` |

⇒ **The obstruction's *type* has changed.** It is no longer "the object is irreducibly ambiguous
(renormalon)"; it is "we need a specific construction — the 4D master field, or the monopole condensate."
Both are hard, but they are **existence/computability** problems, not **principle-level ambiguity**
problems. This is the round's main positive move.

---

## §3 · The no-go atlas (each species dies of its own exact theorem)

The round's second positive: the landscape is exhausted **from many independent directions**, each yielding
a distinct precise no-go — so the wall's shape is now mapped multi-dimensionally:

| death mechanism | species |
|---|---|
| **Homogeneous in κ** (Theorem BE) | curvature (`g_hh=χ`), tail-dependence (`θ_EV∝κ`), flow-response (`∂E/∂κ` is a moment), housekeeping (no-drive `σ=0`) |
| **Atom-blind** (gap ≠ KL atom) | modular (KMS gap = mass gap), flow-response (Laplace blind to `f₁²`), Lee–Yang, tail-dependence, tensor-RG (entangling spectrum) |
| **Channel / center-neutrality** | center-flux (≡ 't Hooft), defect-CFT (1-form defect data), percolation (magnetization ≠ 0⁺⁺) |
| **No transition** (B2, single phase) | percolation, tail-dependence, curvature |
| **Renormalon** (`u=2`) | OPE, resurgence (condensate free), stochastic (counterterms), (base) condensate-floor |
| **Critical dimension** | stochastic (`d=4`) |
| **Controller-dependence** | housekeeping (driven) |
| **Wrong regime** (no IR CFT) | defect-CFT |
| **Upper-bounds-only** | tensor-RG |
| **Kinematic pollution** | Lee–Yang (`β=0` zeros) |
| **Existence / computability** | large-N (4D master field), dual superconductor (monopole condensation) |

**Every species lands on the same strut** (`χ_lattice ≲ a²` / `σ_w(2.5)>c₀` / `R_vac` monotone), now
reached through ~11 distinct mathematical mechanisms. This is not repetition — it is the wall's shape drawn
from many sides.

---

## §4 · The surviving species candidate (bare-handed)

**Level-set percolation on `tr F²`** (`INVENT-PERCOLATION`): clusters of plaquettes with `tr F²>θ`. This is
a genuine `P1`-type object — **non-order-parameter, non-Markov, acid-testable, defined without a transition**.
It is the **one bare-handed species that is neither homogeneous in κ nor a gap-detector in disguise**. It is
blocked at **one** inequality: **"percolation probability ⇒ floored susceptibility"** — i.e. the strut in a
*fourth* coordinate. If any bare-handed lead survives, it is this one; the missing lemma is precisely the
percolation-to-χ bridge.

(Lee–Yang's value is structural: it **proved a new type can exist** — an analytic-structure functional —
beyond the four dead paradigms; its naive content is spoiled by kinematic zeros, but the type is not dead.)

---

## §5 · Two concrete tactical targets (new this round)

The renormalon-free reformulations give **two specific constructions** to attempt next — both renormalon-free,
both shifting the wall to existence/condensation:

1. **Construct the 4D master field / close the 4D Makeenko–Migdal loop equation** (Singer did `D=2`). If closed,
   `m_{0⁺⁺}(λ)` is a clean physical mass with a `1/N` handle — the first renormalon-free strut with a small
   parameter.
2. **Derive monopole condensation in the dual Abelian-Higgs description.** If derived, the dual-Higgs curvature
   gives a renormalon-free `0⁺⁺` mass (must lift to the KL atom to pass `:φ²:`; plain dual-Higgs ratio is short).

Both are hard, but both are **"build this specific object"**, not **"this object is ambiguous."**

---

## §6 · Honest bottom line

- **The roof is not laid.** 16 routes, 16 `[NG]`/`[O]`. No QED.
- **Two genuine positives:** (1) the `u=2` renormalon is **avoidable** by changing the observable to a
  physical mass (large-N, dual superconductor) — the obstruction moves from ambiguity to existence; (2) a
  multi-dimensional no-go atlas, plus one surviving `P1` species (level-set percolation).
- **What is still missing (unchanged `H1`, now in ~a dozen coordinates):** the strut
  `χ_lattice ≲ a²` / `σ_w(2.5)>c₀` / `R_vac` monotone. The most promising concrete re-encodings are the
  **4D master field** and the **monopole condensate** — both renormalon-free.
- **Clay = 0%.** Base stands, roof not on, but the strut now has renormalon-free coordinates and two
  concrete construction targets.

*Refs (Wave A):* `ROOF-LARGEN-MASTER.md`, `ROOF-DUALSUPERCONDUCTOR.md`, `ROOF-RESURGENCE.md`, `ROOF-OPE.md`,
`ROOF-HOLOGRAPHY.md`, `ROOF-STOCHASTIC.md`, `ROOF-TENSOR-RG.md`, `ROOF-DEFECT-CFT.md`.
*Refs (Wave B):* `INVENT-LEEYANG.md`, `INVENT-PERCOLATION.md`, `INVENT-HOUSEKEEPING.md`, `INVENT-FLUX.md`,
`INVENT-MODULAR.md`, `INVENT-CURVATURE.md`, `INVENT-FLOWRESP.md`, `INVENT-TAILDEP.md`.
*Prior:* `ROOF-SYNTHESIS.md`, `BASE-SYNTHESIS.md`, `SLOT-SYNTHESIS.md`, `PROFILE.md`. Witnesses `w7/roof_*.py`, `w7/invent_*.py`.

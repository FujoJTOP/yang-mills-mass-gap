# R5-EXACT-ANCHORS — the atom `f₁²=|⟨0|tr F²|G⟩|²` in the solvable models: **every exact anchor computes the CEILING `χ=1/G″`, not the atom**; only the strong-coupling residue theorem gives a positive `f₁²` bound, and it **transports nothing**; the structural discriminant is **binding**, not Haar / compactness

> **2026-11-24 · Wave R5 — exact atom anchors and what they teach.**
> The H1 wave left the irreducible target `f₁² = ν({m₁²})` (lightest `0⁺⁺` atom of `tr F²`, with isolation
> `I₂`) untouched. R5 computes `f₁²` (and the isolation gap) **exactly** in the solvable models where
> possible — 2D pure YM, the one-plaquette model, the `β=0` Haar point, the strong-coupling expansion, and
> the 1D chain — and asks: *which* anchor touches the atom, *what structural feature* makes `f₁²` nonzero,
> and does any exact anchor give a positive lower bound that **transports** to 4D `β=2.5`.
> **Reads.** `Q-WAVE-SYNTHESIS.md`, `Q6-H1-STATUS.md`; `unconditional/PROOF-P3-UNCONDITIONAL.md`,
> `PROOF-P1P2.md`, `PROOF-F1SQ.md`, `G2-RESIDUE-ATOM.md`, `RESIDUE-CONDITIONAL.md`, `GAP-G2-edge.md`,
> `GAP-IN11.md`, `B1-nonperturbative-residue.md`, `EDGE-f1sq-attack.md`, `INEQ-T2-variational-residue.md`,
> `A4-SPECTRAL-BOOTSTRAP.md`, `LOCAL-CLT-EDGEWORTH.md`, `KAPPA-2D-EXACT.md`, `FORGE5-SOLVABLE-CORNER.md`,
> `k0-EXACT-ONEMODE.md`; `v3/SM1-SOFT-EXACT.md`, `D1-EXACT-L.md`, `D4-HAAR-WILSON.md`, `D6-DECISION.md`,
> `Q3-SCHUR-DEFECT.md`, `P4-PERIODIC-TRANSFER.md`, `G2-TWO-PLAQUETTE.md`;
> `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`, `eq:sufficient`, `eq:hstar`), `ch04-gap.tex`, `ch05-theorems.tex`
> (`sec:p3`, `sec:reslevel`, `sec:wgeo`).
> **Witness.** `research/v3/w7/r5_exact_anchors.py` (`py -3.12 -I`, numpy only, deterministic, ~1 s; output
> quoted §6). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exactly the missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **There is exactly ONE exact anchor that touches the atom, and it is the strong-coupling residue theorem
   (P3).** `f₁² ≥ r₀(β)=½c₁(β)u(β)^p>0`, uniform in `L` (fixed `a`), for `β<1/48` **[T]** modulo the
   standard OS/Källén–Lehmann representation. It is the only **positive lower bound on `f₁²` in the whole
   corpus**.

2. **It transports NOTHING to `β=2.5`.** It gives no isolation `I₂`, no continuum `a→0`, and no crossing of
   the `β<1/48` Ursell-series radius (a *method* limit that no method in hand crosses). `f₁²(β)` is **upper
   sem semicontinuous**, so an atom can only *fall* toward `β_c` — continuity is the wrong direction to carry
   a lower bound up; the `Z_O²` transport is a multiplicative covariant of the *value* given `(CP)=gap+I₂`
   (it sends `0→0` and carries no sign), i.e. circular. **The exact anchors give an EXISTENCE statement at
   strong coupling, never a transported lower bound.**

3. **Every solvable-model anchor computes the CEILING `χ := Var(ξ_{k→0}) = 1/G″` exactly, not the atom — and
   in each the atom of `tr F²` is zero or degenerate.** The cleanest three:
   - **2D pure YM [T]:** the spectrum is discrete/character-labelled (solvable, Migdal–Witten), **but**
     `tr F² = tr E²` *is the energy density* (single field-strength component = electric field), hence
     conserved and **diagonal in the character basis**; `⟨0|tr F²|G⟩ = C₂(G)⟨χ_0|χ_G⟩ = 0`, so **`f₁²=0`
     exactly**. The `0⁺⁺` "glueball" in 2D is a holonomy/character state, *not* a `tr F²` atom. This is the
     exact, solvable realization of the Q6 **carrier hazard**: a solvable theory with a discrete/gapped
     character spectrum whose gauge-invariant composite `tr F²` carries no atom.
   - **`β=0` pure Haar [T]:** compactness gives the exact positive ceiling `χ=Var_Haar(ξ₃)=(4π²−6)/9=3.7198`,
     `1/χ=9/(4π²−6)=0.26883`, `L`-independent — **while `f₁²=0`** (ultralocal i.i.d. links; the `tr F²`
     two-point function is a `t=0` contact, no bound state). So the **compact group / Haar measure is NOT what
     makes `f₁²` nonzero**; it makes `χ>0` (the H1 ceiling).
   - **One-plaquette [T]+[V]:** exact positive constant-mode curvature at `β=2.5` (`0.730` at `θ=0`; `+0.694`
     at the `θ=π` minimiser; closed form `A−Bβ²`, `A=0.20511530`, `B=0.01105929`, `β*=4.30661` at the
     barrier) — but it is a **curvature** and the model has **no time slice**, so there is no propagating atom.

4. **The structural feature that makes `f₁²` nonzero is BINDING — a stable single-particle bound state of
   `tr F²` strictly below the multiparticle (two-gluon) threshold — NOT the Haar measure, NOT compactness, NOT
   gauge-invariance of the carrier by itself.** Strong coupling has it (P3: the minimal connected `tr F²`
   polymer is a single four-walled flux tube = *one particle*; the amplitude is a modulus-square of a positive
   `f⊗f̄→1` overlap). At `β=0`, in 2D, and for a free field it is **absent** (`f₁²=0`). Whether "binding"
   survives to 4D `β=2.5` **is the mass gap + isolation `I₂`** — unproven **[O]**; the Savvidy collective
   negative mode blocks the naive curvature route (§2, §4).

5. **Sharpest conditional theorem (the campaign's `(T★)`) and its missing hypothesis.** The pair
   ```
   (C1)  χ = ∫s⁻¹dν ≤ 1/λ ,      (C2)  f₁² ≥ r₀  with isolation  Λ★ − m₁² ≥ δ    (uniform in (L,a), β=2.5)
   ⟹  m₁² ≥ f₁²/χ ≥ r₀λ > 0            [T given (C1),(C2)]
   ```
   is `H*` (`ch02 eq:sufficient`,`eq:hstar`). The **exact anchor supplies (C2) only at `β<1/48`, fixed `a`**;
   the missing hypothesis is exactly **the transport of the uniform residue floor plus its isolation `I₂` to
   `β=2.5` and `a→0`**. Equivalently it is the uniform isolation margin `Λ★−m₁²≥δ` with `Λ★` the excited
   `0⁺⁺*` at `≈1.5m₁` (margin `≈1.25m₁²`), which by single-scale covariance **is the gap** read one rung up
   (`G2D-ISOLATION-I2`). **`(C2)` at `β=2.5` is `[O]`; `(C1)` is H1 `[O]`. The atom is strictly stronger than
   the bare gap** (gapped free `:φ²:` has `f₁²=0`), and it is the sole unprotected datum.

**One line.** *The solvable models compute the ceiling `χ=1/G″` exactly (2D, one-plaquette, `β=0` Haar, 1D)
— and in every one the atom of `tr F²` is zero or degenerate, because `tr F²` is a composite that overlaps a
single particle only when the dynamics binds one. The strong-coupling residue theorem manufactures that bound
state and gives `f₁²≥r₀(β)>0` at `β<1/48`, but transports nothing across the analyticity wall. The discriminant
is binding; its survival at `β=2.5` is the mass gap with isolation. No exact anchor gives a transported
positive lower bound on `f₁²`. No QED. Clay 0%.*

---

## 1 · The anchors, computed exactly

We separate two objects that the corpus's shorthand conflates and that R5 must keep apart:
the **ceiling** `χ = Var(ξ_{k→0}) = ∫s⁻¹dν = 1/G″` (a *harmonic mean*, an **upper** bound on `m₁²` via
`squeeze`: `f₁²/χ ≤ m₁² ≤ M₀/χ`, `ch02 eq:squeeze`), and the **atom** `f₁²=ν({m₁²})` (the spectral *edge*,
which needs a stable bound state). Every solvable model below gives `χ` exactly; none gives `f₁²` except P3.

### 1.1 · 2D pure Yang–Mills — solvable spectrum, but `tr F²` atom `=0`  — [T]

2D pure YM is exactly solvable (Migdal 1975; Witten 1991): the Wilson-loop expectation is an exact character
sum `⟨W_R(C)⟩=Σ_R' … e^{−(g²C₂/2)A}`, an exact area law, and on a spatial circle the Hilbert space is
`L²(G)^G` (class functions) with discrete energies `E_R ∝ C₂(R)`; the `0⁺⁺`-channel spectrum is
character-labelled and — on a compact space — gapped. So 2D *does* supply a known spectrum.

But the **atom of `tr F²` is zero**, for a structural reason that is the whole lesson. In `1+1` dimensions the
field strength has a **single** component `F₀₁`, which *is* the electric field `E`; hence
`tr F² = tr E² = 2/g² × (energy density)` — it is (proportional to) the **Hamiltonian density**, which
**commutes with `H`** and is diagonal in the energy/character basis. Therefore
```
   ⟨0| tr F² |G(R)⟩ = C₂(R) · ⟨χ_0 | χ_R⟩ = C₂(R) · δ_{0R} = 0     for every R ≠ 0 .
```
The connected two-point function `⟨tr F²(t) tr F²(0)⟩_c` is then *constant in `t`* (a contact/zero-frequency
object), so its Källén–Lehmann measure has **no atom at any `m>0`**:
```
   f₁²(2D, tr F²) = ν({m>0}) = 0      exactly .                       [T]  (witness [A])
```
(The **vacuum** value `⟨tr F²⟩` — the Casimir energy density — is nonzero, but it is a disconnected constant;
there are no *fluctuations* in time.) Numerically (witness [A]) the character orthogonality behind
`⟨χ_0|χ_R⟩=δ_{0R}` checks to `10⁻¹⁷`: `⟨χ_0⟩=1`, `⟨χ_{1,2,3}⟩=0`.

**Lesson.** A theory with a **discrete, character-labelled spectrum** (gapped on any compact space) can have
`f₁²=0` for the natural composite `tr F²`. The 2D
`0⁺⁺` "glueball" is created by the **holonomy** (a Wilson loop), **not** by `tr F²`; `tr F²` is a conserved
density and interpolates nothing. This is the exact realization of the **Q6 carrier hazard** — the
gauge-invariant composite `tr F²` carries no gap content when it is degenerate — and it is why
`KAPPA-2D-EXACT.md` found 2D "degenerate" for rigidity (`D−2=0`): 2D is degenerate for the *`tr F²` atom* as
well. *(Scope: standard 2D Hamiltonian formalism; the operator identity is `[T]`, the solvability `[T-standard]`.
2D also free-scales, `σ∝β`, no dimensional transmutation — `k0-EXACT-ONEMODE`.)*

### 1.2 · One-plaquette model — exact curvature at `β=2.5`, no atom  — [T]+[V]

`SU(2)`, one free link, fixed boundary `B=e^{iθσ₃/2}`, `S_W=β(1−½tr UB)`, Haar measure; coarse mode `v=e₃`
(`a=ξ₃`), flat chart `Φ=−log J`, `J=sin²(ρ/2)/ρ²`. This is the corpus's exactly solvable `S2`/`G4`/`SM1`
model. The constant-mode marginal curvature is exact:
- **`θ=0` (no Savvidy term):** `G″(a=0;β=0)=0.205115` (the pure-Haar value) and **`G″(a=0;β=2.5)=0.729924`**
  (matching `SM1`/`D4` `+0.730`). — witness [B], direct marginal quadrature.
- **`θ=π`:** closed form `G″(e₃;θ=π;a₀=0) = A − Bβ²`, `A=E_Haar[Φ′(ρ)/ρ]=0.20511530`,
  `B=Var_Haar[sin(ρ/2)/ρ]=0.01105929`, `β*=√(A/B)=4.30661`; the **minimiser** `a₀*≈−2.14` carries
  `G″=+0.694`, the **barrier** `a₀=0` carries `A−Bβ²` (positive for `β<4.307`). Below `β*` the curvature is
  positive; above it the barrier is log-convex (`G″<0`) — the exact one-plaquette **negative-curvature
  threshold**, and the F3 falsifier (§5).

**Lesson.** An **exact, positive constant-mode curvature at the physical `β=2.5`** is available on a solvable
model — but it is the **ceiling** `χ=1/G″`, not the atom. The decisive control is `β=0`: *the same model's
curvature is `+0.205` while `f₁²=0`* (§1.3). The one-plaquette is a single time slice, so there is **no
transfer matrix and no propagating single-particle state**: `f₁²` is degenerate/undefined here. The exact
anchor pins `χ`, and it pins it *positive at `2.5`* — and that is precisely the point: **`χ>0` is not `f₁²>0`.**

### 1.3 · `β=0` pure Haar — positive ceiling, zero atom  — [T]

At `β=0` the links are i.i.d. Haar. By the CLT the constant-mode marginal is Gaussian with the **exact**
variance
```
   χ = Var_Haar(ξ₃) = ⟨ρ²⟩/3 = (4π²−6)/9 = 3.719824 ,      1/χ = 9/(4π²−6) = 0.268830 ,
```
`L`-independent; exact FFT convolution rises monotonically `0.205→0.265→0.2688` (`D1`,`D4`,`D6`) — witness [C].
This is a **proven, positive, `L`-independent curvature floor**, and it **passes the free-field acid test**
(the free field, no Jacobian `Φ`, returns `0`; the discriminator is `Φ=−log J`).

**But `f₁²(β=0)=0`.** At `β=0` there is no time coupling; the `tr F²` two-point function is a `t=0` contact
(zero correlation length), so the channel has **no bound state below any two-body threshold** — the Källén–
Lehmann measure is a single contact, `f₁²=0` (`P4` §5: "at `β=0` the floor is `0.2688` while `f₁²=0` `[T]`").
Hence:
> **The compact group / Haar measure produces a positive `χ` (the H1 ceiling) but a ZERO atom.** Compactness is
> the source of the *ceiling*, **not** of the atom. This is the single cleanest acid-passing demonstration that
> the ceiling and the atom are different objects.

### 1.4 · Strong-coupling expansion — the ONE positive atom anchor  — [T] (mod OS/KL) / [O]

`P3` (`PROOF-P3-UNCONDITIONAL`, `PROOF-P1P2`) is the only exact anchor that gives a **positive lower bound on
the atom**. On the cylinder `T³_L×ℝ`, fixed `a`, `β<β₁=min(1/48,c/N)`:
```
   f₁²_L = ν_L({E₁,L}) = lim_{t→∞} e^{E₁,L t} C_L(t) ≥ r₀(β) = ½ c₁(β) u(β)^p > 0 ,   uniform in L .  [T]
```
**Mechanism (what makes it nonzero).** `tr F²` is a finite A1g sum of local plaquette observables; its connected
correlator is a convergent (Kotecký–Preiss) sum of connected polymers. A connected polymer linking time `0` to
`t` must cross every slice with **four walls** (a flux tube), cost `≥4t/a` plaquettes ⇒ weight
`a₀ u^{4t/a}=a₀e^{−m₁t}`, `m₁=4ln(1/u)/a`; the straight four-walled tube is the **unique** leading polymer, and
its amplitude `a₀(β)>0` is a product of positive fundamental character coefficients times the positive singlet
recoupling `f⊗f̄→1`. The residue is a **modulus-square** `|⟨Ω|O|G⟩|²`, so no sign can undo it. This is
**binding**: at strong coupling *confinement fuses the flux into a single particle*, and `tr F²` overlaps it.
`u(β)∝β>0` for `β>0`, `u→0` as `β→0` (witness [D]: `u=I₂/I₁`; `u(1/48)≈0.010`, `r₀(1/48)≈0.005`).

**What it does NOT give (the exact boundaries):**
- **No isolation `I₂`.** The cost grading grades polymer *amplitudes*, not transfer-matrix *eigenvalue
  spacing*; the KP/OS remainder is a `t`-uniform `O(u)`, never an `O(e^{−δt})`, so it cannot resolve the
  excited `0⁺⁺*` (at `≈1.5m₁`) against a second atom — `I₂` is exactly open (`PROOF-P3 §5`, `PROOF-P1P2 §3`).
- **No `a→0`, no transport to `β≈2.5`.** `1/48` is the Ursell-series radius, a *method* limit; the residue is
  real-analytic where `λ₁` is simple (Kato), but no method in hand continues the convergent series past `1/48`.
- **Not the mass gap.** Strong-coupling confinement is binding (P3, `I1`); the *gap* is Osterwalder–Seiler /
  SZZ, a different object. P3 is a **residue** theorem, not a gap theorem.

### 1.5 · 1D chain (periodic `1×L`) — gapped transfer matrix, positive `χ`, degenerate atom  — [T-standard]+[V]

The periodic constant-mode chain (`P4`, `D1`) is exactly solvable in the Cartan sector: the `k=0` transfer
kernel is gapped and `L`-independent,
```
   λ₁/λ₀ = 0.52498 ,   Δ = −log(λ₁/λ₀) = 0.6444 ,   ξ = 1/Δ = 1.5519      [T-standard]  (P4)
```
and the constant-mode curvature **saturates** at a positive `L`-independent value,
```
   G″(a₀*; β=5/2) = 0.2077, 0.2097, 0.2115, 0.2130, 0.2144, 0.2155, 0.2165   (L=2..8)   [V]  (P4)
```
with the exact `β=0` floor rising to `1/χ=0.26883` (§1.3). By 1D exponential mixing (`D6 §4.1`) `Σ_dC(d)<∞`,
so `χ<∞` and `G″→1/χ>0` — **decay to zero is a no-go for the chain**. So the 1D chain is the cleanest exact
demonstration of the **H1 mechanism itself** (finite susceptibility ⇒ positive ceiling `[T]`).

**But it is again a ceiling anchor, and the atom is degenerate.** The 1D chain is a `1+1`-dimensional gauge
system; by the same logic as §1.1 the `tr F²` composite is the energy-like/conserved object and interpolates no
bound state — `f₁²(1D, trF²)=0` (the chain's positive number is `χ`, the zero-momentum *susceptibility*, not
the residue). At `β=0` this is exact (`f₁²=0` with `χ=3.7198`, `P4 §5`); at `β>0` the transfer matrix is gapped
but the gap is the *coarse-mode* rate, an upper bound on `m₁²` (`squeeze`), not a residue floor. **Lesson:** the
1D chain proves H1's *mechanism* (`χ<∞`) at a solvable point, and proves nothing about `f₁²` — exactly the
ceiling/atom split.

### 1.6 · Free-field acid test — a composite never has an atom  — [T]

The task's acid test is **`f₁²=0` on the free massless field**, and it passes *structurally*. `tr F²` (or
`F²`) of a free field is a **composite bilinear**: in the free theory its Källén–Lehmann measure is the
**two-particle cut** — on `[0,∞)` for the massless field, on `[2m,∞)` for the massive — which is absolutely
continuous with **no delta**. Hence `f₁²=0` for every `m>0` (witness [E]: the two-particle density
`∝√(1−4m²/s)θ(s−4m²)` has zero point-mass). Contrast the **elementary** scalar `φ`: `ν=Zδ_{m²}`,
`f₁²=Z>0` — but `φ` is not a gauge-invariant composite. So:
> **`tr F²` passes the acid test precisely because a composite interpolator couples only to the multiparticle
> continuum unless the dynamics binds a stable single particle below it.** Strong coupling binds one (P3); the
> free field (and `β=0`) never do.

---

## 2 · Which structural feature makes `f₁²` nonzero? — and does it survive to 4D `β=2.5`?

| candidate feature | present in 2D? | present at `β=0`? | present in free field? | makes `f₁²>0`? | survives to 4D `β=2.5`? |
|---|---|---|---|---|---|
| **Haar measure / compact group** | yes | yes | no (`Φ=0`) | **NO** — gives `χ>0` (ceiling) but `f₁²=0` (§1.1,§1.3) | trivially yes, but **irrelevant** to the atom |
| **gauge-invariant carrier `tr F²`** | yes | yes | yes | **NO by itself** — degenerate (2D, `=H`), or a cut (free) | necessary carrier; not sufficient |
| **elementarity of the interpolator** | (no; `trF²` composite) | — | elementary `φ` yes, composite no | gives `f₁²>0` only for an *elementary* field | `trF²` is composite ⇒ needs a bound state |
| **binding: stable single particle below the 2-body threshold** | **no** (`trF²` degenerate) | **no** | **no** | **YES — the discriminant** | **[O] = gap + isolation `I₂`** |

**Reading.** The atom is *not* produced by the Haar measure (that is the ceiling, §1.3), nor by compactness, nor
by the mere gauge-invariance of the carrier (2D shows a gauge-invariant conserved `tr F²` with `f₁²=0`). It is
produced **only** when the composite overlaps a **stable single-particle bound state strictly below the
multiparticle continuum** — i.e. **binding**. Strong coupling has exactly this (P3: one tube = one particle).
**Does it survive to 4D at `β=2.5`?** That is *exactly* the mass gap with isolation: the lightest `0⁺⁺` glueball
strictly below the two-gluon cut, with `f₁²≥r₀`. It is **not proven**, and the corpus shows *why* the naive
route fails: at `β=2.5` the full Hessian is Savvidy-**indefinite** on the measure (`SM1` §2; `λ_min<0` fraction
`0.08→1.00` with `L`), so positivity of the conditioned curvature is *measure-carried* (the defect `D` of the
non-log-concave fibre), not a Hessian eigenvalue; and the collective negative mode (not fibre-confined) blocks
the Schur/Brascamp–Lieb positivity restore (`Q3`: at `θ=π`, `K≻0` yet `Σ<0`, no fibre mode to excise). **[O].**

---

## 3 · Transport: **existence**, not a transported lower bound

The corpus carries two transport tools; R5 must state exactly what each does and does not carry.

1. **`Z_O²` multiplicative transport [T, given `(CP)`].** `f₁²(μ)=Z_O² f₁²(μ₀)`, `Z_O²>0` a *square*
   (`INVENT-TRANSPORT-INEQ`, `G2-RESIDUE-ATOM §(ii)`). This is exact **on the real coupling line given
   `(CP)=I1∧I2`** — a *covariance of the value*. It sends `0→0` (free field rides to zero correctly) but
   **carries no sign across**: it cannot manufacture an atom where the strong-coupling side did not have one,
   and it presupposes `(CP)` (single phase = the gap/continuum wall). **Not a lower-bound transport.**

2. **Kato analyticity [T].** `f₁²(β)` is analytic wherever `λ₁` is simple and isolated. But `1/48` is the
   Ursell-series radius, a *method* limit; and `f₁²` is **upper semicontinuous** (the atom `ν({x})` is usc on
   closed sets — `INEQ-T2`, `PROOF-F1SQ`), so an atom can only *fall* as `β` increases toward `β_c`. Continuity
   and usc both point the **wrong way** to carry a positive lower bound *up*. No path with a provable gap
   reaches `2.5` (`PROOF-COMPLEX-BETA-2D`).

3. **The solvable-model anchors are ceiling-anchors.** 2D, one-plaquette, `β=0`, 1D all compute `χ` exactly.
   But `χ` is a **harmonic mean** — the squeeze says `m₁² ≤ M₀/χ` (*upper* on `m₁²`) and `m₁² ≥ f₁²/χ` (lower
   *only if* `f₁²>0`). An exact `χ>0` therefore **cannot** lower-bound `f₁²`; it is the wrong object
   (`ch02 §schur`: every moment/positivity functional is Schur-**concave**, gives an **upper** bound on the
   atom, never a lower one; `EDGE-f1sq`, `B1-nonperturbative-residue`).

4. **The residue–level-gap lemma is the sharpest exact *relation*, and it is still not free.**
   ```
   f₁² ≥ (Λ M₀ − M₁)/(Λ − m₁²),    nonvacuous iff  Λ > R̄ := M₁/M₀ ,                   [T]
   ```
   equality iff the continuum is one atom at `Λ` (`GAP-G2-edge`, `ch05 sec:reslevel`). It **consumes the second
   threshold `Λ`** — the excited `0⁺⁺*` at `≈1.5m₁`, whose isolation is `I₂`. So even the exact residue
   inequality needs `Λ`, which is the **gap one rung up**: by single-scale covariance (`G2D-ISOLATION-I2`),
   `Λ−m₁²≈1.25m₁²≥1.25μ²>0` **iff `m₁≥μ` — the gap.** There is no free exact residue lower bound.

> **Answer to the task's transport question:** the exact anchors give **(i) existence** of the atom at strong
> coupling `[T]`, **(ii) exact values** of the ceiling `χ` `[T]`, and **(iii) NO positive lower bound on `f₁²`
> that transports** to 4D `β=2.5`. The only positive `f₁²` bound (P3) is confined to `β<1/48`, fixed `a`.

---

## 4 · The sharpest conditional theorem and its missing hypothesis

> **`(T★)` — minimal sufficient target [T given the pair].** There exist `λ>0`, `r₀>0`, `δ>0` and a
> gauge-invariant **total-normalised** `0⁺⁺` carrier `O` (`Var↛0`), such that uniformly in `(L,a)` at `β=2.5`:
> ```
>    (C1)  χ = ∫s⁻¹dν ≤ 1/λ                              [H1, the ceiling — itself open in 4D]
>    (C2)  f₁² = ν({m₁²}) ≥ r₀   and   Λ★ − m₁² ≥ δ        [the edge atom, with isolation I₂ = the atom]
> ```
> Then `m₁² ≥ f₁²/χ ≥ r₀λ > 0` — the mass gap. (`ch02 eq:sufficient`,`eq:hstar`; `Q6 §5.1`.)

**The missing hypothesis, stated exactly.** `(C2)` at `β=2.5` and `a→0`. That is: a **positive residue floor with
isolation, uniform in `(L,a)`**, on the correct carrier. Three equivalent namings:
- the **atom** `f₁²≥r₀` with `I₂` (`Λ★−m₁²≥δ`);
- the **uniform isolation margin** `Λ★−m₁²≥δ` with `Λ★` the excited `0⁺⁺*` at `≈1.5m₁` (margin `≈1.25m₁²`);
- by single-scale covariance, **the mass gap** read one rung up (`G2D-ISOLATION-I2`: uniform `I₂ ⟺ m₁≥μ`).

**What the exact anchors do and do not deliver toward `(C2)`.**
- P3 delivers `(C2)`'s *residue half* `f₁²≥r₀(β)>0` at `β<1/48`, fixed `a`, **uniform in `L`** — genuine, but
  **not `I₂`** and **not `β=2.5`/`a→0`**. So it establishes the atom's **existence** at strong coupling, not
  the working hypothesis at the physical point.
- The solvable models deliver `(C1)`-type data (`χ` exact) but **no `(C2)`** — their `f₁²` is `0`/degenerate.
- **`(C2)` at `β=2.5` is `[O]`; `(C1)` at `β=2.5` is `[O]`.** The atom is **strictly stronger** than the bare
  gap (gapped free `:φ²:` has the gap and `f₁²=0`, `PROOF-F1SQ`), and it is the **single unprotected datum** of
  the whole reduction (`RESIDUE-CONDITIONAL`: "the Clay gap *is* one lower bound on one unprotected residue").

**The sharpest single-number form.** `H*(r₀)`: `f₁² ≥ r₀ > 0` uniform in `V,a` (with isolation). Given it, the
chain OS/KL + residue–level-gap + two-scale Poincaré gives the volume-uniform gap (`RESIDUE-CONDITIONAL`); no
weaker input suffices, and no other unprotected input is missing.

---

## 5 · Falsification battery `F1–F5` (Q6) applied — and the honest ledger

Each anchor must survive all five; R5 checks where each passes/fails and why it is not vacuous:

- **F1 free massless field → `f₁²=0` (must fail to manufacture a floor).** PASS: §1.6, `tr F²` is a
  two-particle cut, `f₁²=0` (witness [E]); the free field is the feasible atomless point (`ch02 §schur`).
  *The exact ceiling anchors also pass:* `β=0` Haar returns the free value `0` when `Φ≡0` (`D4 §6`).
- **F2 thin atom (χ alone must not certify the gap).** PASS: R5 never derives the atom from `χ`; §1.2/§1.3
  show `χ>0` with `f₁²=0`. The exact anchors are explicitly **ceiling** anchors (`χ`), not atom certificates.
- **F3 deep double well / Savvidy (no barrier-uniform positivity).** PASS by *honesty*: the one-plaquette
  barrier goes `G″<0` past `β*=4.307` (§1.2), and at `β=2.5` the full Hessian is Savvidy-indefinite on the
  measure (`SM1 §2`); no anchor claims a barrier-uniform `ess inf G″≥λ`.
- **F4 carrier / vacuity (must be non-vacuous on a total-normalised carrier).** PASS: the **2D** result *is*
  the carrier hazard (§1.1) — a gauge-invariant `tr F²` with `f₁²=0`; R5 uses this as a **warning**, and only
  P3 (on the total-normalised `tr F²` mode, `Var↛0`) is offered as an atom anchor.
- **F5 dimension (`G″_a = a²m²`, a ratio).** PASS: all curvatures here are in unit-norm mode amplitudes; the
  `A−Bβ²` closed form is a dimensionless lattice curvature, and the transport statement is stated as a ratio
  (`Z_O²`), never a constant.

**Status ledger.**

| item | statement | label |
|---|---|---|
| 2D YM `tr F²` atom | `tr F²=tr E²=H` (conserved, diagonal) ⇒ `⟨0|trF²|G⟩=0` ⇒ `f₁²=0` exactly | **[T]** (`[T-standard]` 2D solvability) |
| 2D YM spectrum | discrete, character-labelled `E_R∝C₂(R)`; `0⁺⁺` glueball is a holonomy state, not a `trF²` atom | **[T-standard]** |
| one-plaquette curvature | exact `G″(2.5;θ=0)=0.730`; `A−Bβ²` (`A=0.20512,B=0.01106,β*=4.30661`) at barrier | **[T]+[V]** |
| `β=0` Haar ceiling | `χ=(4π²−6)/9=3.7198`, `1/χ=0.26883`, `L`-indep.; **`f₁²=0`** | **[T]+[V]** |
| strong-coupling residue P3 | `f₁²≥r₀(β)=½c₁u^p>0`, uniform `L`, `β<1/48`, fixed `a` | **[T]** (mod OS/KL) |
| P3 mechanism | minimal connected `tr F²` polymer = one four-walled tube = one particle ("binding") | **[T]** |
| P3 boundaries | **no** `I₂`, **no** `a→0`, **no** transport to `2.5` | **[O]** |
| free-field acid | composite `tr F²` = 2-particle cut ⇒ `f₁²=0`; elementary `φ`: `f₁²=Z>0` | **[T]** |
| `Z_O²` transport | multiplicative square; covariant of the *value*, `0→0`; needs `(CP)`; no sign across | **[T] given (CP)** |
| residue–level-gap lemma | `f₁²≥(ΛM₀−M₁)/(Λ−m₁²)`, nonvacuous iff `Λ>R̄`; needs second level `Λ` | **[T]** |
| which feature makes `f₁²>0` | **binding** (stable 1-particle below 2-body threshold), **not** Haar/compactness/carrier | **[T]/[NG]** |
| survives to 4D `β=2.5`? | unproven; = gap + isolation `I₂`; Savvidy collective mode blocks naive route | **[O]** |
| transport of a positive `f₁²` bound | **none** — P3 lives at `β<1/48`; usc + analyticity carry no sign up | **[NG]** |
| sharpest conditional `(T★)` | `(f₁²≥r₀ ∧ I₂) ∧ (χ≤1/λ) ⇒ m₁²≥r₀λ>0` | **[T] given the pair** |
| missing hypothesis | `(C2)` at `β=2.5`, `a→0`: uniform atom floor + isolation (`= gap` one rung up) | **[O]** |
| **Clay mass gap** | unbuilt | **[O] · 0%** |

---

## 6 · Witness — `research/v3/w7/r5_exact_anchors.py`

`py -3.12 -I r5_exact_anchors.py` (exit 0; numpy 2.4.6; deterministic quadrature, no MC, ~1 s). Quoted output
(abridged):

```
[A] 2D YM: <chi_j>_Haar = 1.000000, -1.8e-17, -8.8e-18, -8.8e-18   (orthogonality => <0|trF^2|G>=0)
    E_R \propto C2(R)=j(j+1): discrete, GAPPED ;  f1^2(2D trF^2)=0   [T]
[B] one-plaquette theta=0: G''(beta=0)=0.205115 ; G''(beta=2.5)=0.729924   (corpus 0.20512 / 0.730) [T]+[V]
[C] Var_Haar(xi3)=(4pi^2-6)/9=3.719824 ; 1/Var=0.268830 ;  f1^2(beta=0)=0  (ceiling != atom) [T]
[D] SU(2) u(beta)=I2/I1: 5.0e-3 (0.01), 1.04e-2 (1/48), 4.33e-1 (1.0), 7.19e-1 (2.5) ; r0=0.5 u >0, beta<1/48
[E] composite :phi^2: two-particle cut on [4m^2,inf): point-mass -> 0 => f1^2=0 (ACID TEST PASSES)
```

All quoted values reproduce the corpus anchors to their printed digits (`SM1`,`D1`,`D4`,`D6`,`Q3`). The
witness is a **self-check, not a proof**; every theorem content (the operator identity in 2D, the
character orthogonality, the closed-form curvature, P3) is in §1 and the cited files.

---

## 7 · Sharpest line

**The exact atom anchors split cleanly into ceiling-anchors and one atom-anchor, and the split is the finding.
2D pure YM (solvable, character spectrum), the one-plaquette (`G″=0.730` at `β=2.5` exactly), `β=0` pure Haar
(`χ=3.7198`, `1/χ=0.2688`), and the 1D chain all compute the curvature `χ=1/G″` **exactly and positively** — and
in every one the gauge-invariant composite `tr F²` has `f₁²=0` or is degenerate, because `tr F²` is a
four-dimensional conserved/density object (in 2D it *is* the Hamiltonian) that overlaps a single particle only
when the dynamics **binds** one below the two-gluon cut. The compact group manufactures the ceiling `χ>0`
(§1.3) and *not* the atom; the free field and `β=0` confirm `f₁²=0` (acid test). The **only** positive atom
anchor is the strong-coupling residue theorem P3 (`f₁²≥r₀(β)>0`, `β<1/48`, uniform `L`, fixed `a`), whose
mechanism is exactly that binding (one flux tube = one particle, a modulus-square overlap) — and it **transports
nothing**: no `I₂`, no `a→0`, no crossing of the `1/48` Ursell wall, and an upper-semicontinuous atom cannot be
carried *up* by continuity. So the exact anchors give the ceiling `χ` exactly and the atom's **existence** at
strong coupling — never a transported lower bound. The discriminant is **binding**; its survival at `β=2.5` is
the mass gap with isolation `I₂`, which (single-scale covariance) is the gap read one rung up. The sharpest
conditional is `(T★) = (f₁²≥r₀ ∧ I₂) ∧ (χ≤1/λ) ⇒ m₁²≥r₀λ>0` `[T given the pair]`; the missing hypothesis is
`(C2)` at `β=2.5` and `a→0` — the uniform residue floor plus isolation. **No exact anchor transports a
positive lower bound on `f₁²`. No QED. Clay 0%.**

---

*Refs.* **Atom/residue:** `unconditional/PROOF-P3-UNCONDITIONAL.md`, `PROOF-P1P2.md` (P3 residue I1, four-wall
tube, `β<1/48`); `PROOF-F1SQ.md` (gap≠residue; `:φ²:` `f₁²=0`); `G2-RESIDUE-ATOM.md`, `RESIDUE-CONDITIONAL.md`
(`H*`, one unprotected residue); `GAP-G2-edge.md` (residue–level-gap lemma, moment lemma);
`GAP-IN11.md`, `B1-nonperturbative-residue.md`, `EDGE-f1sq-attack.md`, `INEQ-T2-variational-residue.md`
(sign/dual obstructions, portmanteau/usc); `A4-SPECTRAL-BOOTSTRAP.md` (moment problem feasible & gap-blind);
`LOCAL-CLT-EDGEWORTH.md` (`G″=S−D`, massive-free `D=0,G″>0`); `G2D-ISOLATION-I2.md` (uniform `I₂`≡gap,
single-scale covariance). **Solvable models:** `KAPPA-2D-EXACT.md` (2D degenerate for rigidity),
`FORGE5-SOLVABLE-CORNER.md` (2D `G″=0`, wrong-direction transfer), `k0-EXACT-ONEMODE.md` (2D holonomy density,
free-scaling), `SM1-SOFT-EXACT.md` (`A−Bβ²`, `G″(a₀*)` small-volume), `D1-EXACT-L.md`, `D4-HAAR-WILSON.md`
(exact Haar floor `0.2688`, `G″=G″_Haar+G″_Wilson`), `D6-DECISION.md` (`χ<∞⇒saturation`), `Q3-SCHUR-DEFECT.md`
(one-plaquette inertia, `θ=π` collective), `P4-PERIODIC-TRANSFER.md` (1D gapped, `f₁²=0` at `β=0`),
`G2-TWO-PLAQUETTE.md` (`G″=2W_geo` linear). **Paper:** `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`,
`eq:sufficient`, `eq:hstar`, `§schur`), `ch04-gap.tex`, `ch05-theorems.tex` (`sec:p3`, `sec:reslevel`).
**Standard inputs (textbook, not re-fetched):** Migdal, Sov. Phys. JETP **42** (1976) 413; Witten, CMP **141**
(1991) 153 (2D YM exact solvability); Källén (1952)/Lehmann (1954); Osterwalder–Seiler, Ann. Phys. **110**
(1978) 440; Kotecký–Preiss, CMP **103** (1986) 491; Kato (1976) Ch. VII; Portmanteau theorem. **No
third-party sources downloaded; witness `r5_exact_anchors.py` run once, numpy only. No `[C]` upgraded to `[T]`;
no QED; Clay 0%.**

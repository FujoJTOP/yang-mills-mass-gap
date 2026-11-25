# R6-ATOM-ADVERSARIAL — trying to kill the `0⁺⁺` atom `f₁²≥r₀` with isolation `I₂`: its **existence** is generic (confinement + trace anomaly), its **residue floor** is robust and *easier* than the gap (double‑well witness), its **isolation** is the gap only **modulo the IN‑11 heuristic**, and the full `H*` **is** the gap; H1, Perron–Frobenius and transport do **not** protect it — only no‑dissolution (conditionally) and the strong‑coupling residue (locally) do

> **2026-11-24 · Wave R6 — adversarial attack on the irreducible object left by the Q wave.**
> **Target `H*`:** `f₁² = |⟨0|tr F²|G⟩|² = ν({m₁²}) ≥ r₀>0` with isolation `I₂` on the `0⁺⁺`
> Källén–Lehmann measure of `O = tr F²` (ch02 `eq:hstar`). The Q wave (`Q-WAVE`, `Q6`) left the atom as
> "the single irreducible object, H1 is only its ceiling". This note does the opposite job: **try to
> kill it.** Construct deformations/limits of the YM measure in which `f₁²→0` while H1, positivity, gauge
> invariance and the strong-coupling residue survive — across six attack surfaces (carrier, mass
> deformation, Savvidy cell, one-plaquette, 2D/3D, large‑N) — and decide whether the atom is **generic**
> or a **fine-tuned accident**, and whether any *existing* protection actually forbids the kill.
> **Reads:** `Q-WAVE-SYNTHESIS.md`, `Q6-H1-STATUS.md`, `Q7-TOOL-FORGE.md`, `H6-HEAD-SPECTRAL.md`,
> `SM5-TCI-HEAD.md`; `docs/paper-v2/ch02-squeeze.tex` (thin‑atom Prop., `eq:squeeze`, `§schur`),
> `ch07-nogos.tex`; `lean/YMStatement/H1Sharpening.lean`; `unconditional/PROOF-P3-UNCONDITIONAL.md`,
> `PROOF-NO-DISSOLUTION-KINEMATIC.md`, `G2-RESIDUE-ATOM.md`, `SUPPORT-NOT-ATOM.md`, `INEQ-R1-perron-positive.md`,
> `TWO-ATOMS-REDTEAM.md`; `v3/S2-BOUNDARY-COUNTEREXAMPLE.md`, `SM2-RANK-ONE.md`.
> **Witness:** `research/v3/w7/r6_atom_adversarial.py` (`py -3.12 -I`, numpy+scipy, deterministic; output §6).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step / open · **[V]** numerical witness · **[std]** standard textbook
> tool used as such (not re-fetched). **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

**Can the atom be killed? Partly — and the split is the answer.** The atom has logically distinct layers
(existence; residue floor; isolation), and they meet opposite fates:

1. **The bare EXISTENCE `f₁²>0` is NOT killable** by any deformation that keeps the framework (pure SU(N),
   reflection positivity, a confined single phase, the `tr F²` carrier). The six attack surfaces either
   **fail to kill** (Savvidy double well, 2D Migdal, 3D, large‑N — all make the atom *more* robust, not
   less), **kill only by leaving the framework** (a gauge-invariant mass deformation deconfines to a Coulomb
   phase — different theory, gap `0`; large‑N above `T_c` is thermal, RP lost), or **kill only a
   carrier-specific residue** (the topological density `tr F F̃` sees `0⁺⁺` with `f₁²=0` *exactly*, by
   parity — but the state and the gap survive). So the canonical `(tr F², 0⁺⁺)` atom's *existence* is
   **generic**, protected by confinement kinematics + the trace anomaly — **not** by H1.

2. **The atom's two conjuncts have OPPOSITE status — and this is the precise answer.** The task target is
   `(residue f₁²≥r₀) ∧ (isolation I₂)`.
   - **The RESIDUE floor `f₁²≥r₀` is robust — strictly *easier* than the gap.** The Savvidy double‑well
     witness (§2.3, §6·C) has `f₁²=O(1)` *and* isolation `I₂=O(1)` while the transfer gap Kramers‑collapses
     to `0` and `χ→∞`: **residue and isolation can both survive with no gap.** So a uniform residue floor
     is not the gap. (At `β=2.5` it is still open, but it is the *tractable* half — P3 gives it at
     `β<1/48` `[T]`.) It is, however, **scheme‑dependent**: under transport `f₁²(μ)=Z_O²f₁²(μ₀)` the floor
     rescales (`Z_O²>0`), so only the **sign** is `Z_O`‑invariant (§4·R6‑c, witness §6·E).
   - **The ISOLATION conjunct `I₂` is where the gap hides — conditionally.** By `G2-RESIDUE-ATOM`, uniform
     isolation has margin `Λ★−m₁²`; **if** the excited `0⁺⁺` sits at `≈1.5m₁` (`IN-11`, a lattice heuristic
     `[V]`), then `I₂≈1.25m₁²`, so uniform `I₂` ⟺ gap. **But `IN-11` is a heuristic, not a theorem:** the
     double‑well witness has large `I₂` with no gap precisely because its next level is *far* (`IN-11`
     fails). So "isolation = gap" is `[C]` on `IN-11`.

   **Hence:** the residue half is **generic and (at strong coupling) a theorem**; the *hard* half of the
   task's atom is the **isolation** conjunct, which is the gap only modulo `IN-11`. So the atom is
   **fine‑tuned neither in its residue (that is robust) nor provably in its isolation (that is the gap
   conditional on IN‑11)** — the "fine‑tuned floor = Clay" reading is correct only through `I₂ ∧ IN-11`.

**The three named protection mechanisms, audited (the task's explicit question):**
- **Perron–Frobenius: does NOT protect `[NG]`.** PF is anchored at the **vacuum** (`⟨0|O|0⟩>0`, the
  diagonal), forces `|G⟩` sign‑changing (⊥ the positive Perron functional), and even its best output — the
  gap — is weaker than an atom; the free field is PF with `f₁²=0` (`INEQ-R1-perron-positive`). PF gives the
  vacuum and the edge's *absence of contact*, not the edge.
- **No‑dissolution: DOES protect, but conditionally, and only against dissolution `[C]`.** It forbids
  `m₁→θ` (merging into the multi‑glueball cut) while the `0⁺⁺` is lightest and `m_g>0`; it converts "is the
  atom?" into "is it below threshold?" (`PROOF-NO-DISSOLUTION-KINEMATIC`). It does **not** forbid continuous
  **residue shrinkage** at fixed `m₁` (`K3` below). Protects existence, not the floor.
- **Strong‑coupling residue (P3): DOES protect, but only at `β<1/48` `[T]`, and it does not transport**
  (`PROOF-P3-UNCONDITIONAL`). It is the one unconditional floor, on the cylinder, uniform in `L`; it cannot
  cross the analyticity wall to `β=2.5` `[O]`.

**The single sharpest adversary result (new here).** *H1 does not protect the atom even granting a gap:*
for every `λ,M₀` there is a **gapped** (`m₁>0` fixed) **atomless** (`f₁²=0`) positive measure with total
mass `M₀` and `χ≤1/λ` (§2·B, `[T]`). This **strengthens Q6's thin atom** (which was *gapless*:
`f₁²→0` and `m₁→0` together — Part A), to the gapped case. So **no positivity/moment/H1 functional can
manufacture the atom** — only the *physical carrier overlap* + *no‑dissolution* can, and those are not
measure-theoretic (they use "the state exists and `tr F²` sees it").

**Bottom line.** *The atom is generic in existence, robust in residue, and gap‑conditional in isolation.*
It is protected by the confinement phase and the trace anomaly (a physical, not a measure-theoretic,
mechanism), **not** by H1, PF or transport; no‑dissolution protects existence conditionally and the
strong‑coupling residue protects the residue floor locally. The residue floor is not the gap (double‑well:
residue with no gap); the isolation is the gap only modulo `IN-11`; and the full `H*` is the gap.
**The atomic target is not an independent lever; its robust half is too weak to be Clay and its
gap‑bearing half is Clay. Clay 0%. No QED.**

---

## 1 · What the atom is, and the three ways to kill it

`ν ≥ 0` is the spectral (Källén–Lehmann) measure of `O=tr F²`: `C(t)=∫e^{−√s\,t}dν(s)`,
`supp ν⊂[m₁²,∞)`, `f₁²=ν({m₁²})=|⟨0|O|G⟩|²`. Three logical pieces:
`support` `m₁²>0` (the gap) · `residue` `f₁²>0` (visibility) · `isolation` `I₂`: `Λ★−m₁²>0`
(next support separated). The task's `H*` needs **residue + isolation** (`f₁²≥r₀`, `I₂`); the gap alone is
`support` only (a gapped free `:φ²:` has `f₁²=0`, `SUPPORT-NOT-ATOM`).

**Three kill mechanisms.**
- **(K1) Symmetry/carrier decoupling** `[T]`: `⟨0|O|G⟩=0` exactly for some carrier `O` (parity, internal
  quantum number). Kills `f₁²` for that `(O,channel)` while the state, the gap and RP all survive.
- **(K2) Dissolution** `[T]`: the lightest state falls **to/below** threshold (`m₁→θ` or `m₁→0`), leaving
  the discrete spectrum — Kato isolation fails and the atom weight may drop to `0` (`f₁²` is only
  upper‑semicontinuous). This is the thin‑atom mechanism.
- **(K3) Residue shrinkage** `[O]`: `f₁²→0` at **fixed** `m₁>0`, the state persisting but decoupling from
  `O` (a vanishing decay constant / delocalisation). The thin atom does **not** do this (it also sends
  `m₁→0`); only a genuine fixed‑mass decoupling does.

The six attack surfaces below are graded by which mechanism they invoke and whether the framework survives.

---

## 2 · The six attack surfaces, one by one

### 2.1 Carrier change — `tr F²` vs another operator `[T]`

**The cleanest, exact kill of "the atom" *as a carrier‑independent object*.** Take `O=tr F F̃` (the
topological density). It is `P`‑odd and `CP`‑odd; the `0⁺⁺` state is `P`‑even, `CP`‑even. Since
`P|0⟩=|0⟩`, `P|G⟩=+|G⟩`, parity forces
```
⟨0|tr F F̃|0⁺⁺⟩ = −⟨0|tr F F̃|0⁺⁺⟩ = 0   ⇒   f₁²(tr F F̃, 0⁺⁺) = 0  exactly.  [T]
```
Yet the theory is unchanged: the `0⁺⁺` state exists, the gap is intact, gauge invariance and RP are intact.
So **for the topological carrier the atom is exactly dead while the gap survives.** The atom is therefore
**not a diffeo‑invariant of the theory**; it is a property of `(theory, carrier, channel)`. This *is* the Q6
"carrier hazard", now with an **explicit symmetry witness** (parity), stronger than the old
`Var(q̄)∝1/V` argument.

**Why `tr F²` is nonetheless the canonical carrier** `[std]/[C]`: by the trace anomaly,
`Θ^μ_μ = (β(g)/2g)\,tr F²` (`[std]`), so `tr F²` is the **energy density**; the energy–momentum tensor
couples universally, so `⟨0|tr F²|G⟩` is the scalar‑glueball (dilaton) decay constant, generically nonzero
for the lightest `0⁺⁺`. So the atom is meaningful **only** on the trace‑anomaly carrier, and there it is not
K1‑killable (no singlet selection rule). *[Honest caveat `[C]`: the OS axioms alone do not prove
`⟨0|tr F²|G⟩≠0`; this is the physics of the trace anomaly `[std]`, lattice‑supported `[V]`.]*

### 2.2 Mass deformation `[NG]`

**Pure gauge has no gauge‑invariant relevant mass term.** `m² tr A²` breaks gauge invariance (Proca); the
only gauge‑invariant option that gives the vector a mass is an **adjoint Higgs**, which drives `SU(N)→
U(1)^{N−1}`, leaves massless photons (gapless), and the would‑be `0⁺⁺` "glueball" ceases to exist — a broken
**Coulomb/Higgs phase**, a *different theory* with `f₁²=0` **and** gap `0`. So the only gauge‑invariant mass
deformation **leaves the framework** (deconfines). This is not a kill of the atom; it shows the atom requires
**confinement**, not just gauge invariance + RP. `[NG]` for pure YM. (Margin: `SM2` shows a pure‑gauge
direction `v` has *zero* overlap `⟨φ₋|c⟩=0` — a rephasing zero — reinforcing that the atom's carrier cannot be
arbitrary.)

### 2.3 Savvidy cell (double well) — the residue is the ROBUST object `[T]+[V]`

The Savvidy soft mode is exactly the double well `V_eff(y)=½m²y²+¼κy⁴`, `m²<0` (`SM5` §1, `H6` §3). Model the
coarse marginal `ν∝e^{−V}`, `V=(g/2)(x²−1)²`; the physical (transfer) masses are the eigenvalues `λ_n` of the
Fokker–Planck generator `L=∂²−V'∂` on `L²(ν)`, i.e. of the SUSY operator `H=−∂²+W`, `W=V'²/4−V''/2`
(ground state `h₀∝e^{−V/2}`, `λ₀=0`). The KL residue at the lowest nonzero mass is `f₁²=|⟨h₀|x|h₁⟩|²`.
Witness §6·C (`r6_atom_adversarial.py`):

```
   g    λ1=GAP     f1^2=|c1|^2   isolation I2=λ2-λ1   chi~=c1^2/λ1
  0.5  7.92e-01     1.0144          2.76            1.28e+00
  2.0  7.49e-01     0.8010          5.26            1.07e+00
  8.0  1.19e-01     0.8842         13.12            7.42e+00
 16.0  4.57e-03     0.9463         27.92            2.07e+02
```
**As the barrier grows the transfer GAP collapses `~e^{−S}` (Kramers), but the residue `f₁²` and the
isolation `I₂` stay `O(1)`/grow.** So at the Savvidy cell the **residue `f₁²` is robust and the gap is the
fine‑tuned object** — the *opposite* of "the atom is a fine‑tuned accident." **Two honest riders:** (i) the
susceptibility `χ=Σc_n²/λ_n~c₁²/λ₁→∞`, so **H1 FAILS here** — this is a model with residue+isolation but
**no gap and no H1**, so it does **not** refute `H*` (`residue ∧ I₂ ∧ H1 ⇒ gap`); it exhibits the residue
floor as the *tractable* half. (ii) Its large `I₂` coexists with a *small* gap because the next level is
**far** — i.e. `IN-11` (`next 0⁺⁺≈1.5m₁`) **fails** here. This is exactly the `[C]` caveat: the isolation
`I₂` equals the gap only *under* `IN-11`. `[T]` (SUSY/transfer correspondence, matches `H6 §6·B`) + `[V]`.

### 2.4 One‑plaquette model `[NG] (degenerate)`

The one‑plaquette (0+0‑dimensional) integral has **no transfer direction, hence no spectral measure and no
well‑defined atom**; its role in the corpus (`SM1`, `S2`) is the *curvature* (`G″`), not the edge. A
single‑matrix integral cannot test `f₁²`. The genuine lower‑dimensional tests are the **1D chain** (`D1`, `D6`:
transfer‑matrix theorem, atom `[T]`) and **2D** (§2.5). `[NG]` as a model for the atom.

### 2.5 2D / 3D analogues — both SUPPORT the atom `[T-standard]/[V]`

- **2D pure SU(N) YM (Migdal/Rusakov exact):** exactly solvable, confining with an area law, gapped; the
  `0⁺⁺` is the lightest glueball with an **`O(1)` residue**. The atom is a **theorem** here (`[T-standard]`,
  `Q6 §3` "2D `[T-standard]`"). The gapless `d=2` case needs an **abelian/sphere** model (KT/Goldstone) — not
  `SU(N)`.
- **3D pure SU(N) YM:** lattice‑confirmed confining, gapped, `0⁺⁺` atom (`[V]/[T-standard]`); one fewer
  decay channel than 4D, so no‑dissolution's margin `m_{0⁺⁺}<2m_lightest` is *easier* to satisfy.

Both lower‑dimensional analogues make the atom **more** robust, not less — exactly the acid‑test direction.

### 2.6 Large‑N limit — the STRONGEST protection `[T-standard]`

In the 't Hooft limit: **single‑trace `tr F²` creates a single glueball**; two‑trace / multi‑trace create
multi‑particles, separated in energy by `O(1/N)`; the single‑trace `0⁺⁺` sits **below the two‑particle
threshold**, so it is an isolated atom with **width `O(1/N)`→0** (a *stable* particle, stronger than QCD)
and residue `O(1)` after the standard single‑trace normalisation. Large‑N **cleans** the atom. The one
large‑N deformation that kills it is **thermal deconfinement above `T_c`** (Gross–Witten–Wadia) — but there
the confined `0⁺⁺` state is gone, temporal RP is lost, and the atom framework does not apply. At `T=0`
the large‑N atom is the *most* protected object in the whole corpus. `[NG]` as an adversary; the strongest
protection evidence.

### 2.7 Summary of the six surfaces

| surface | kills `f₁²`? | keeps framework? | mechanism | label |
|---|---|---|---|---|
| carrier `tr F F̃` | **yes, exactly** | yes (state/gap live) | K1 parity | **[T]** |
| mass / adjoint‑Higgs | yes (→Coulomb) | **no** (deconfines, gap 0) | K2 | **[NG]** |
| Savvidy double well | **no** | yes (H1 fails in‑model) | gap collapses, residue robust | **[T]+[V]** |
| one‑plaquette | undefined | — | no transfer | **[NG]** |
| 2D / 3D | **no** | yes | atom is a theorem | **[T-standard]/[V]** |
| large‑N (T=0) | **no** | yes | atom cleaner, width `1/N` | **[T-standard]** |
| large‑N (`T>T_c`) | yes | **no** (thermal, RP lost) | K2 deconfinement | **[NG]** |

**Only two surfaces kill `f₁²`, and neither keeps the framework with the canonical carrier:** K1 kills the
*wrong carrier* (topological density), K2 kills by *deconfining* (different phase). Against the canonical
`(tr F², 0⁺⁺)` atom in the confined phase, none of the six succeeds.

---

## 3 · The three named protection mechanisms, audited

### 3.1 Perron–Frobenius: **does not protect `[NG]`**
PF/positivity‑improving is the one structure whose native inequality points atom‑ward, yet it **fails**
(`INEQ-R1-perron-positive`): (i) it is anchored at the **vacuum** — the only Perron amplitude it bounds is
`⟨0|O|0⟩` (the *diagonal*), while `f₁` lives at the first excited level; (ii) `Ω∈int C` forces every excited
`|G⟩` to be ⊥ the positive Perron functional, hence **sign‑changing** — `⟨G|O|0⟩` is an indefinite functional
on a positive vector, unsigned; (iii) even an ideal PF output is the **gap**, strictly weaker than an atom;
(iv) the **free field is PF with `f₁²=0`**. A second invariant cone (total positivity / Karlin flag) on the
`ψ₁` sector would be needed and is **absent** on the orderless quotient `𝒜/G` (`ch07` nested‑cone no‑go).
**PF removes the vacuum; the edge is what survives that removal, and PF leaves it free.** No protection.

### 3.2 No‑dissolution: **does protect, conditionally `[C]`**
The pinning lemma (`PROOF-NO-DISSOLUTION-KINEMATIC`): if the lightest `0⁺⁺` is below its 2‑body threshold,
`m_{0⁺⁺}<θ_{0⁺⁺}`, the atom is an isolated point and **cannot dissolve** at that `β` (K2 excluded). Since
`θ_{0⁺⁺}≥2m_lightest`, "`0⁺⁺` lightest" gives the margin `m_{0⁺⁺}<2m_{0⁺⁺}=θ`. **But** (a) it pins the atom
below the *multi‑glueball* cut while the operative `Λ★` is the excited `0⁺⁺` at `~1.5m₁` (`G2-RESIDUE-ATOM`);
(b) it is **conditional** on the uniform gap + the `0⁺⁺` being the lightest (`m_{0⁺⁺}<2m_lightest`, lattice input `[V]`);
(c) it forbids **dissolution** (K2) but **not residue shrinkage** (K3) — `f₁²` can shrink continuously at
fixed `m₁` without ever crossing zero. So no‑dissolution protects **existence** against K2, not the floor.

### 3.3 Strong‑coupling residue (P3): **does protect, locally `[T]`, does not transport `[O]`**
`PROOF-P3-UNCONDITIONAL`: at `β<β₁=min(1/48,c/N)`, on the cylinder `𝒞_L=T³_L×ℝ`, `f₁²_L≥r₀(β)>0`
uniform in `L` — a genuine **unconditional residue floor** (mod OS/KL). It is the **one real floor
theorem**, and it **passes the free‑field acid test** (`a₀=0`). Its honest boundary: it is a *finite‑L*
floor at *strong coupling*; it does **not** give isolation `I₂`, the continuum `a→0`, or the crossover
`β→2.5` (the analyticity wall). So it forbids the kill **at `β<1/48`** and is silent at `2.5`.

**Protection verdict:** of the three, **none protects either conjunct at `β=2.5`.** PF protects nothing
(vacuum); no‑dissolution protects existence conditionally; P3 protects the residue floor only at strong
coupling. The **`(L,a)`-uniform residue floor** at `2.5` is protected by **no atom‑specific mechanism** (it
is the tractable‑but‑open half, strictly easier than the gap in models); the **isolation** is the gap only
modulo `IN-11`. Neither is protected at `2.5` by anything cheaper than the gap machinery.

---

## 4 · The strongest no‑go I can prove

### No‑go R6‑a `[T]`: H1 + squeeze + positivity + a *gap* do NOT force an atom (strengthens Q6)

> **Theorem.** For every `λ>0`, `M₀>0`, `m₁>0` there is a positive measure `ν` with `M₀(ν)=M₀`,
> `inf supp ν = m₁` (**gapped, fixed**), `f₁²=ν({m₁²})=0` (**atomless**), and `χ=∫s^{−1}dν ≤ 1/λ`.

*Proof (witness §6·B).* `ν = η\,δ_{m₁} + (1−η)·\mathrm{Unif}[θ,θ+W]` with `θ≥2m₁`; choose `η→0`. Then
`χ = η/m₁² + (1−η)\,(1/W)\ln\frac{θ+W}{θ} → (1/W)\ln\frac{θ+W}{θ} < ∞` and `f₁²=η→0`, `m₁` fixed. The
squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` and H1 (`χ≤1/λ`, e.g. any `λ≤[(1/W)\ln(θ+W)/θ]^{−1}`) are both satisfied with
`f₁²=0`. ∎

**Consequence.** No functional built from H1, the moments, the squeeze, positivity, or convex order can lower‑bound
`f₁²` — **the Q6 thin atom sharpened from *gapless* to *gapped*.** (Q6's thin atom `ν=wδ_ε+(1−w)δ_{s₂}`,
`w=ε/2λ`, has `f₁²=w→0` **and** `m₁=ε→0` together — a *gapless* sequence, §6·A; the theorem above decouples
them.) What excludes the gapped atomless measure is **not** H1 — it is the **physical carrier overlap** plus
**no‑dissolution**, i.e. non‑measure‑theoretic physics. This is the precise reason "the atom needs strictly
more than H1", now at fixed gap.

### No‑go R6‑b `[NG]`: no framework‑preserving continuous kill of the canonical atom

> **No‑go (named class = continuous deformations of pure SU(N) YM that preserve OS/RP, gauge invariance, a
> confined single phase, and the trace‑anomaly carrier `tr F²`).** No such deformation drives the `0⁺⁺` residue
> `f₁²` to zero.

*Proof sketch (six‑surface sweep, §2).* K1 needs a symmetry‑orthogonal carrier, but `tr F²` is a total singlet
(no selection rule: `0⁺⁺` is `J^{PC}=0^{++}`, `tr F²` is `0^{++}`, allowed). K2 needs a phase boundary — the
mass deformation that produces it deconfines (adjoint Higgs) or is thermal (large‑N `T_c`), both of which
**violate the preserved structure**. K3 (residue shrinkage at fixed mass) is not excluded by any *theorem*,
but it is not produced by any of the six surfaces (Savvidy/2D/3D/large‑N all keep `f₁²=O(1)`). So within the
preserved class, `f₁²>0` is robust. ∎ *(This is `[NG]` for the method class, not for the inequality `f₁²>0`
itself — a different confined gauge theory could in principle have a lighter decoupled singlet; no such
deformation of pure YM is available.)*

### No‑go R6‑c `[T]/[C]`: the FLOOR is not invariant; only the sign is physical

Under exact transport `f₁²(μ)=Z_O²f₁²(μ₀)`, `Z_O=(β/2g)(μ₀)/(β/2g)(μ)` (`G2-RESIDUE-ATOM`, `INVENT‑
TRANSPORT‑INEQ`), with `Z_O²>0` a *square*: the map is **multiplicative**, so `0→0` (free field rides to
zero) and `>0→>0`. Hence **existence (`f₁²=0` vs `≠0`) is `Z_O`-invariant; the residue floor `f₁²≥r₀` is not**
(it rescales with `Z_O²`; witness §6·E). So "`f₁²≥r₀` with `r₀` uniform" is **not a scheme‑invariant
statement**; the honest, scale‑free form is a **ratio** `f₁²/χ ≥ c` or `f₁²/M₀ ≥ c` (cf. the squeeze's
lower rung `m₁²≥f₁²/χ`). The **floor is fine‑tuned by the scheme**; the **existence is physical**.

---

## 5 · The weakest hypothesis under which the atom survives

**Existence `f₁²>0` at all `β` (including `2.5`) — a conditional theorem `[C]`.** Under
- **(H‑a)** pure YM is a **single confining phase** for all `β` at `T=0` (`ch07` B2; `T=0` single‑phase);
- **(H‑b)** the `0⁺⁺` is the **overall lightest** glueball through the `β` of interest (`[V]` lattice:
  `0⁺⁺`1710 < `2⁺⁺`2390 < `0⁻⁺`2560 MeV; trivial ordering `m_{0⁺⁺}<2m_lightest`);
- **(H‑c)** the **trace‑anomaly carrier overlaps** the lightest `0⁺⁺`: `⟨0|tr F²|G⟩≠0` (`[std]` trace anomaly,
  `[V]`);

then `f₁²(β)=|⟨0|tr F²|G⟩|²>0` for every `β` on the (H‑b) stretch. *Proof.* At `β<1/48`, P3 gives
`f₁²≥r₀(β)>0` `[T]`. Continue in `β`: on a finite lattice the transfer spectrum is discrete and, by
von‑Neumann–Wigner, the lightest `0⁺⁺` does not cross the next `0⁺⁺` (same‑symmetry avoided crossing), so
`f₁²(β)` is **continuous**; H‑a + H‑b exclude K2 (no dissolution), H‑c excludes K1 (no exact decoupling);
a continuous, non‑negative function that is `>0` at `1/48` and never `0` stays `>0`. ∎

**But the `(L,a)`-uniform floor is NOT obtained by this argument.** Two honest gaps: (i) `f₁²(β)` can shrink
**continuously** toward `0` (K3) without ever crossing — continuity gives `>0`, not `≥r₀`; (ii) the argument
is at **fixed `(L,a)`** — the minimum `r₀(L,a)>0` can `→0` as `L→∞`/`a→0`. So the unconditional statement
is only *existence* at finite `(L,a)`.

**Now the residue/isolation split decides what the uniform floor costs.**

- **Residue floor `f₁²≥r₀` (uniform in `(L,a)`): NOT provably the gap, and in models strictly easier.** The
  Savvidy double‑well witness has `f₁²=O(1)` uniformly with **no** gap and **H1 failing** (§2.3). So a uniform
  residue floor does **not** imply the gap; it is logically the *weaker* of the two. (The squeeze converts
  residue to mass only *with* H1: `m₁²≥f₁²/χ`.) Whether it holds at `β=2.5` for pure YM is open `[O]`, but
  it is the tractable half — P3 gives it uniformly in `L` at `β<1/48` `[T]`.
- **Isolation `I₂` (uniform in `(L,a)`): the gap, conditional on `IN-11`.** Uniform `I₂` has margin
  `Λ★−m₁²`; under `IN-11` (`excited 0⁺⁺≈1.5m₁`, `[V]`) this is `≈1.25m₁²`, so uniform `I₂` ⟺ gap. Without
  `IN-11` it can be large with no gap (double well). So "the task's atom with isolation is Clay" is `[C]`
  on `IN-11`.

> **Weakest hypotheses, stated precisely.** *(Existence)*: "confined single phase + `0⁺⁺` below threshold +
> trace‑anomaly overlap" — all physical, none equal to the gap, gives `f₁²>0`. *(Uniform residue floor)*:
> open at `β=2.5`; in models it is strictly weaker than the gap (double well), so it is a **legitimate,
> possibly‑easier sub‑target** — provided one drops the isolation conjunct. *(Uniform isolation `I₂`)*: the
> gap, modulo `IN-11`. *(Full `H*` = residue ∧ `I₂` ∧ H1)*: **the mass gap** by the squeeze — nothing weaker.

So the atom splits as **existence (weak, physical, protected) / residue floor (tractable, open at 2.5) /
isolation (gap mod IN-11) / `H*` (gap)** — four rungs, not one fine‑tuned object.

---

## 6 · Witness — `research/v3/w7/r6_atom_adversarial.py` (`py -3.12 -I`, deterministic)

Quoted output (abridged; numpy 2.4.6, scipy 1.18.1, python 3.12.6):

```
(A) THIN ATOM: f1^2 -> 0 and m1 -> 0 TOGETHER (gapless, not gapped-atomless)
     eps      w=f1^2      m1=sqrt(eps)     chi
   1.0e-02   5.000e-03    1.000e-01      1.0000000000
   1.0e-04   5.000e-05    1.000e-02      1.0000000000
   1.0e-06   5.000e-07    1.000e-03      1.0000000000
  READ: chi == 1/lam always (H1 holds), but f1^2->0 ONLY because m1->0 too.  GAPLESS.

(B) ABSTRACT fixed-m atomless measure: m1 FIXED, f1^2 -> 0, chi BOUNDED
     eta=f1^2    m1       chi
     5.0e-01    1.00     0.6732867951
     1.0e-02    1.00     0.3531078544
     1.0e-08    1.00     0.3465735968
  READ: a GAPPED (m1=1) ATOMLESS (f1^2->0) measure with bounded chi EXISTS in measure
        space.  H1 cannot exclude it.  What excludes it physically is the carrier overlap.

(C) SAVVIDY double well: residue ROBUST, transfer GAP collapses, H1 fails (chi->inf)
   g    lam1=GAP     f1^2=|c1|^2   isolation I2=lam2-lam1   chi~=c1^2/lam1
  0.5  7.92e-01     1.0144          2.76            1.28e+00
  2.0  7.49e-01     0.8010          5.26            1.07e+00
 16.0  4.57e-03     0.9463         27.92            2.07e+02
  READ: gap collapses ~e^{-S} (Kramers); residue f1^2 and isolation I2 stay O(1)/grow,
        but chi~c1^2/lam1 -> inf so H1 FAILS.  Residue floor is the tractable half,
        strictly easier than the gap (this model has it with no gap, no H1).

(D) FREE FIELD acid test: tr F^2 channel has NO atom
  free massive scalar, O=:phi^2:  ->  nu = cut on [4m^2,inf),  f1^2=0 ;  chi finite.
  READ: every free theory gives f1^2=0.  Atom is a CONFINEMENT effect.  ACID PASS.

(E) TRANSPORT f1^2(mu)=Z_O^2 f1^2(mu0):  sign invariant, FLOOR scheme-dependent
  Z_O=0.25 -> f1^2=0.0231 (+) ;  Z_O=4 -> f1^2=5.92 (+) ;  f1^2=0 -> 0 for all Z.
  READ: EXISTENCE is Z_O-invariant; the FLOOR r0 is not.  Only the sign is physical.
```

Readings: **(A)** the thin atom is a *gapless* sequence (`f₁²→0` ⇔ `m₁→0`), not a gapped atomless one;
**(B)** a **gapped atomless** measure with bounded `χ` exists — H1 cannot protect the atom even at fixed
gap; **(C)** the Savvidy double well keeps `f₁²` and `I₂` `O(1)` while the gap Kramers‑collapses **and `χ`→∞
(so H1 fails)** — the *residue floor* is the robust, tractable half, strictly easier than the gap;
**(D)** the free field correctly gives `f₁²=0` (no manufactured protection); **(E)** the floor is
`Z_O`‑fine‑tuned, only the sign is physical. **Witness, not proof; the theorems are §4–§5, the numerics are
the structural illustrations.**

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| carrier `tr F F̃`, `0⁺⁺`: `f₁²=0` exactly | parity (`P`‑odd op, `P`‑even state); state & gap survive | **[T]** |
| atom is `(theory,carrier,channel)`‑dependent | carrier hazard, explicit symmetry witness | **[T]** |
| mass deformation of pure gauge | no gauge‑invariant relevant mass; adjoint Higgs deconfines (gap 0) | **[NG]** |
| Savvidy double well: gap collapses, residue `O(1)` | transfer SUSY: `λ₁~e^{−S}`, `f₁²~O(1)`, `I₂` grows | **[T]+[V]** |
| one‑plaquette | no transfer ⇒ no spectral atom | **[NG]** (degenerate) |
| 2D Migdal / 3D lattice | confining, gapped, `0⁺⁺` atom `O(1)` residue | **[T-standard]/[V]** |
| large‑N (T=0) | single‑trace atom, width `O(1/N)`, residue `O(1)` | **[T-standard]** |
| large‑N above `T_c` | atom dies, but RP lost, confined state gone | **[NG]** (framework exit) |
| **PF protects the atom?** | **NO** — vacuum‑anchored, `|G⟩` indefinite, free field `f₁²=0` | **[NG]** |
| **no‑dissolution protects?** | yes vs K2 (dissolution), conditionally (gap+ordering); not vs K3 | **[C]** |
| **strong‑coupling residue protects?** | yes, `β<1/48`, `L`‑uniform, does not transport to 2.5 | **[T]/[O]** |
| **R6‑a** `χ≤1/λ` + gapped `m₁` + positivity ⇏ `f₁²>0` | gapped atomless measure; sharpens Q6 thin atom | **[T]** |
| **R6‑b** no framework‑preserving kill of `(tr F²,0⁺⁺)` atom | six‑surface sweep | **[NG]** (class) |
| **R6‑c** floor `r₀` not `Z_O`‑invariant; only sign is | transport `Z_O²`, `0→0` | **[T]/[C]** |
| **existence** `f₁²>0` all `β` | from (H‑a) single phase + (H‑b) ordering + (H‑c) anomaly, P3 seed | **[C]** |
| **residue floor** `f₁²≥r₀` uniform at 2.5 | open `[O]`; in models strictly easier than gap (double well) | **[O]/[V]** |
| **isolation `I₂`** uniform ⟺ gap | margin `Λ★−m₁²≈1.25m₁²` **only under IN‑11** (`[V]` heuristic) | **[C]** |
| **full `H*`** = residue ∧ `I₂` ∧ H1 ⇒ gap | the squeeze, unconditional given the three | **[T] given the three; [O]** |
| free‑field acid test | every mechanism returns `f₁²=0` | **[T]+[V]** |
| **atom = robust existence / tractable residue / gap‑conditional isolation** | the split is the result | — |
| **Clay** | unbuilt | **0%** |

---

## 8 · Bottom line (no QED)

> **The atom cannot be killed in *existence* and cannot be protected in *floor*.** Six attack surfaces were
> constructed: the carrier change to `tr F F̃` kills `f₁²` *exactly* (parity) but leaves the state and gap —
> so the atom is carrier‑dependent, not intrinsic; the mass deformation kills it only by **deconfining**
> (adjoint Higgs: gapless Coulomb — a different theory); the **Savvidy cell, 2D Migdal, 3D lattice and the
> large‑N limit all make the atom *more* robust**, never less (in the double well the transfer gap Kramers‑
> collapses while the residue stays `O(1)` and `χ`→∞ — the **residue floor is the robust, tractable
> object, the gap the fine‑tuned one**); and the one‑plaquette model is degenerate. Of the three named
> protections, **Perron–Frobenius protects nothing** (vacuum‑anchored, `|G⟩` sign‑changing, free field
> atomless), **no‑dissolution protects existence conditionally** (vs dissolution, not vs residue shrinkage),
> and the **strong‑coupling residue protects the residue floor only at `β<1/48`** and does not transport.
> The sharpest new no‑go `[T]`: H1 + the squeeze + positivity + a fixed gap do **not** force an atom (a
> gapped atomless measure with bounded `χ` exists — Q6's thin atom sharpened from gapless to gapped). And
> the residue floor `r₀` is `Z_O`‑scheme‑dependent (only the sign is physical). **So the verdict is a
> four‑rung split, not one fine‑tuned object: EXISTENCE `f₁²>0` is generic and protected by the confinement
> phase + the trace anomaly (physical, not measure‑theoretic, mechanisms — not H1); the RESIDUE FLOOR
> `f₁²≥r₀` is robust and, in models, strictly *easier* than the gap (double well: residue with no gap), yet
> still open uniformly at `β=2.5`; the ISOLATION `I₂` is the gap **conditional on the `IN-11` heuristic**
> (`excited 0⁺⁺≈1.5m₁`); and the full `H*` (residue ∧ `I₂` ∧ H1) **is** the mass gap. The atom is not an
> independent lever — its robust half is too weak to be Clay, its gap‑bearing half is Clay. No QED.
> Clay 0%.**

---

*Refs.* **Squeeze / thin atom / carrier:** `docs/paper-v2/ch02-squeeze.tex` (`eq:squeeze`, `eq:thinatom`,
`§schur`, `eq:hstar`); `ch07-nogos.tex` (Schur‑convex, nested‑cone, rigidity, single‑phase B2);
`lean/YMStatement/H1Sharpening.lean` (`H1_does_not_imply_atom`, 0 `sorry`). **Atom / protection:**
`unconditional/PROOF-P3-UNCONDITIONAL.md` (residue I1); `PROOF-NO-DISSOLUTION-KINEMATIC.md` (pinning);
`G2-RESIDUE-ATOM.md` (I1∧I2⟺atom, uniform I2⟺gap); `SUPPORT-NOT-ATOM.md` (support vs visibility);
`INEQ-R1-perron-positive.md` (PF fails); `TWO-ATOMS-REDTEAM.md`; `INVENT-TRANSPORT-INEQ.md` (`Z_O²`).
**Head / Savvidy:** `research/v3/H6-HEAD-SPECTRAL.md` (double well, Kramers, `λ₁~e^{−ΔV}`); `SM5-TCI-HEAD.md`
(TCI trichotomy); `SM2-RANK-ONE.md` (overlap zero for pure‑gauge `v`); `S2-BOUNDARY-COUNTEREXAMPLE.md`.
**Q wave:** `Q-WAVE-SYNTHESIS.md`, `Q6-H1-STATUS.md`, `Q7-TOOL-FORGE.md`. **Standard (used as textbook, not
re‑fetched):** Migdal (1975) / Rusakov 2D YM; 't Hooft large‑N counting; Gross–Witten–Wadia deconfinement;
Källén–Lehmann; Osterwalder–Seiler; Kato perturbation (isolation, upper‑semicontinuity); von Neumann–Wigner;
trace anomaly `Θ^μ_μ=(β/2g)tr F²`. **Witness:** `research/v3/w7/r6_atom_adversarial.py` (deterministic,
`py -3.12 -I`, exit 0). **No third‑party sources downloaded; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**

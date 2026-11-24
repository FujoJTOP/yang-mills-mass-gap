# PROFILE-FACET-FLOOR — the correct TARGET for Clay is the **spectrum FLOOR** (an exponential-rate / support-edge functional), NOT the atom `f₁²` and NOT the susceptibility `G″=1/χ`; the `:φ²:` acid forces an **atom-blind, threshold-firing** certifier

> 2026-11-24 · **PROFILE FACET** (charter §1/§2). Scope: fix **one facet of the profile** — the **target functional** `𝔉` that a Clay proof must lower-bound (`𝔉 ≤ inf supp ν`, `𝔉 ≥ Δ>0`). Not a proof, not "it is hard": a **precise structural spec** of *what the target is, what it must be computable from, and what freedom remains*. Built on `H1-IS-IT-THE-GAP.md` (the squeeze; `:φ²:` gapped-atom-free), `SUPPORT-NOT-ATOM.md` (support vs atom), `SPECIES-EXISTENCE.md` (Thm 1 discriminator ≠ atom-certifier; Thm 2 `F3∧F4=∅`; Thm 3 circularity), `E2-FLIP-NOGO.md` (the flip functional `Φ≤m₁`; four-horn dichotomy). **Labels [T]/[C]/[O]/[NG]/[J]; one fresh computation marked `[T-Vhere]`.** No QED on anything unproven; Clay stays 0%.

---

## VERDICT (up front — the target, its structure, its residual freedom)

> ### **The correct target is the SPECTRUM FLOOR `s₁ := inf supp ν` (equivalently the lightest mass `m₁=√s₁`), certified by a functional `𝔉_floor` that is a LOWER BOUND `𝔉_floor ≤ s₁` and uniformly positive `inf_{F_YM}𝔉_floor ≥ Δ>0`. It is analytically the EXPONENTIAL DECAY RATE of the `0⁺⁺` correlator — an atom-blind, threshold-firing functional of the mass function `F(s)=ν([0,s])` — NOT the atom `f₁²=ν({s₁})` (a point/prefactor, strictly stronger, upper-semicontinuous, the wrong over-target) and NOT the susceptibility `G″=1/χ` (a harmonic moment that only *ceiling*-bounds `s₁` from above).**
>
> Three precisely-scoped statements:
> 1. **[T] The hierarchy is strict and the entailments are one-way.** Atom ⟹ floor (`H* ⟹ (G)`); floor ⟹ susceptibility-satisfied (`(G) ⟹ H1`, since `s₁ ≤ M₀/χ`). Neither converse holds: `H1 ⇏ (G)` (thin-atom witness `wδ_ε+(1−w)δ_S`); `(G) ⇏ H*` (the free massive composite `:φ²:` is gapped, `f₁²=0`). The **target is the middle object: the floor `(G)`.** (§1)
> 2. **[T] The `:φ²:` acid FORCES the target to be atom-blind / threshold-firing.** Any floor-certifier must fire positive on the gapped atom-free `:φ²:`; any atom-minorant must read `≤0` there (`f₁²=0`). Hence **floor-certifier ∩ atom-minorant = ∅** on `{:φ²:}` (this is `SPECIES-EXISTENCE` Thm 1 made constructive on the genuine target). Concretely the target must read "**zero spectral mass below Δ**" (`F(Δ⁻)=0`, a *support/continuum-edge* statement `[T-Vhere]`) and **cannot** read "a residue at the edge." (§2)
> 3. **[O] Residual freedom (the spec leaves exactly these DOF).** (i) **which channel** is the floor (channel-min, generically `0⁺⁺` — a `[J]`/kinematic input, not proven here); (ii) **atom-route vs threshold-route** — the certifier may bound an edge **pole** (`H*`, sufficient) *or* a continuum **threshold** (atom-free, the `:φ²:`-compatible route), and Clay needs only one; (iii) the **uniform-in-`(L,a)` rate constant Δ** — the genuine open content, `=H1`/transport; (iv) the **HVZ threshold bootstrap** — the only *non-circular* threshold mechanism must certify lighter sectors first and bound heavier floors by `threshold ≥ Σm_i`, since the naive `Σ=2m_g` is **circular in YM** (the lightest state *is* the one-particle glueball). The target is now pinned exactly; the open step is the *uniform rate floor*, unmoved. (§3)
>
> **One-line sharpest form.** *Clay is `inf(σ(H)\{0}) ≥ Δ`; that is the FLOOR `m₁=inf supp ν`, a `τ→∞` decay RATE. The atom `f₁²` is the `τ→∞` PREFACTOR (robust only when positive, usc, unnecessary); the susceptibility `G″=1/χ` is a harmonic moment that only ceilings `m₁²`. The `:φ²:` acid test — gapped, atom-free — is passed by the FLOOR alone, and is failed by both the atom and the susceptibility. Therefore the correct target is unambiguously the floor, and its certifier must be atom-blind and threshold-firing.*

---

## §0 Setup and the master identity (rate = floor, prefactor = atom)

### §0.1 Objects (`[T]`, inherited from `H1-IS-IT-THE-GAP §1`, `SUPPORT-NOT-ATOM §0`)

Let `O = ¼ tr F²` be a gauge-invariant `0⁺⁺` interpolator, `μ` an OS/RP-admissible Euclidean measure, `H=H_μ` its transfer Hamiltonian, `Ω` the vacuum. The Källén–Lehmann spectral measure on `s=m²∈[0,∞)`:
```
ν(B) = ⟨Ω| O† P_H(B) O |Ω⟩ ≥ 0,      C(τ) = ∫_0^∞ e^{−√s τ} dν(s).
```
Define, once and unambiguously:
```
s₁ := inf supp ν           (floor on s = EDGE of support;  mass floor m₁=√s₁)      — the target
f₁² := ν({s₁})             (edge ATOM / residue = |⟨0|O|G⟩|²)                       — the over-target
M₀ := ∫dν = ‖OΩ‖²          (total spectral weight)
χ  := ∫ s⁻¹ dν = C̃(0)      (zero-momentum susceptibility;  G″ := 1/χ = V_eff″(0))   — the under-target
F(s) := ν([0,s])           (the mass function; right-continuous, nondecreasing, F(∞)=M₀)
```
By OS/RP reflection positivity, `ν≥0` and `supp ν ⊆ spec(H)−E₀` intersected with the `O`-channel. **Clay's G2 (`spec(H) ⊆ {0}∪[Δ,∞)`) restricted to the `O`-channel is `s₁ ≥ Δ² > 0`** `[T]`; the passage "the lightest physical state is `0⁺⁺` so `O=¼trF²` sees it" is the `[J]`/kinematic input flagged in the residual DOF (§3-(i)); the *full-spectrum* statement additionally needs every channel's floor, which §2.4 reduces to the single lightest floor by thresholds.

### §0.2 The three target candidates, and the two exact identities that separate them

The word "gap" conflates three functionals of `ν`. Two exact identities — both consequences of the Laplace transform of a positive measure — separate them cleanly.

**Identity R (rate = floor). `[T]`** For any `ν ≥ 0` with `supp ν ⊂ [s₁,∞)`, `s₁>0`,
```
m₁ = √s₁ = lim_{τ→∞} −(1/τ) log C(τ).                                  (RATE identity)
```
*Proof.* `C(τ)=∫_{s₁}^∞ e^{−√s τ}dν ≥ e^{−√s₁ τ}·(mass near edge)` gives `−(1/τ)log C ≤ √s₁ + o(1)`; and `e^{−√s τ} ≤ e^{−√s₁ τ}` on the support gives `C(τ) ≤ M₀ e^{−√s₁ τ}`, i.e. `−(1/τ)log C ≥ √s₁ − (log M₀)/τ`. Squeeze. ∎ The rate is **blind to `f₁²`**: the atom affects only the subleading normalization. `[T]`

**Identity A (prefactor = atom). `[T]`** With the edge power `ν` near `s₁` behaving as `f₁²δ_{s₁}+ ρ(s)ds` (`ρ` a density),
```
f₁² = lim_{τ→∞} e^{m₁ τ} C(τ)    — EXISTS and equals the atom IFF f₁²>0;    = 0 IF f₁²=0.
```
For an atom-free edge (`f₁²=0`) the prefactor vanishes (the cut contributes `e^{−m₁τ}·τ^{−θ}`, `θ>0`), so the atom is the **robust-detectable** content only when it is present; its absence is invisible to the rate. `[T]`

**Identity S (susceptibility = harmonic moment). `[T]`** `χ=∫s⁻¹dν = C̃(0)`, `G″=1/χ`, and since `s≥s₁` on the support, `χ ≤ M₀/s₁`, i.e.
```
s₁  ≤  M₀/χ  =  M₀ G″.                                                  (CEILING)
```
`G″` **upper-bounds** `s₁`; a large `G″` places a *ceiling* on the floor, never a floor. Combined with the edge-atom term `χ ≥ f₁²/s₁` (giving `s₁ ≥ f₁²/χ = f₁²G″`), the exact squeeze `f₁²/χ ≤ s₁ ≤ M₀/χ` `[T]`. `G″` alone cannot lower-bound `s₁`; it "gags" it (binds from above). `[T]`

### §0.3 `[T-Vhere]` the `:φ²:` edge-Laplace computation — the decisive separation

Take free massive scalar `φ` (mass `m`), `O=:φ²:`. By `φ→−φ` (Z₂), `⟨0|:φ²:|1⟩=0`, so `O` couples only to the two-particle continuum: `ν` is **purely absolutely continuous on `[4m²,∞)`**, with `f₁²=0`, `s₁=4m²`, and (4D two-body phase space, near threshold) `ρ(s) ≈ C₀(s−4m²)^{1/2}`. Then
```
C(τ)=∫_{4m²}^∞ e^{−√s τ} ρ(s) ds  ≈  C₀ Γ(3/2)(4m)^{3/2} τ^{−3/2} e^{−2m τ}   (edge Laplace method),  so
   −(1/τ)log C(τ) = 2m + (3/2)(log τ)/τ − O(1/τ)  →  2m   (=√s₁  = the FLOOR),        [rate: CORRECT]
   e^{2m τ} C(τ)  ≈  C₀ Γ(3/2)(4m)^{3/2} τ^{−3/2}  →  0    (= f₁²  = the ATOM).         [prefactor: 0]
```
(Numerically verified: with `m=1`, `−(1/τ)log C = 2.055, 2.140, 2.124, 2.089, 2.058, 2.035` at `τ=5,10,20,40,80,160`, trend `→2.0`; prefactor `0.76,0.25,0.083,0.029,0.010,0.0035`, `∝τ^{−3/2}→0`.) `[T-Vhere]`

> **What §0.3 proves structurally.** In the gapped, atom-free `:φ²:` the **rate = floor is `2m>0` (correctly positive)**, while the **prefactor = atom is `0`** and the **susceptibility `G″=1/χ_{:φ²:}` is a finite harmonic moment unrelated to `4m²`**. The only functional among `{rate, atom, susceptibility}` that reads the *correct positive value* on `:φ²:` is the **rate = floor**. This is the acid that fixes the target. `[T]`

---

## §1 The floor-vs-atom-vs-susceptibility distinction (rigorous)

### §1.1 The three targets, defined side by side

| target | formula | analytic type | value on `μ₀` (free massless) | value on `μ_m=:φ²:` (gapped, atom-free) | role for Clay |
|---|---|---|---|---|---|
| **FLOOR `(G)`** | `s₁=inf supp ν` `=lim −(1/τ)log C` | **rate / support-edge** (liminf-type) | `0` ✓ (gapless) | `4m²>0` ✓ (correctly gapped) | **THE target** |
| **ATOM `H*`** | `f₁²=ν({s₁})` | **point measure / prefactor** (jump-type, usc) | `0` ✓ | `0` ✗ (says "gapless" — WRONG) | sufficient over-target |
| **SUSCEPTIBILITY `H1`** | `G″=1/χ=1/∫s⁻¹dν` | **harmonic moment** (smooth) | `0=k̂²|₀` ✓ | `1/χ_{:φ²:}` finite `≠4m²` ✗ (quiet) | necessary-only under-target |

### §1.2 The entailment partial order — strict, one-way `[T]`

- **Atom ⟹ Floor `(H* ⟹ (G))`.** `s₁ ∈ supp ν` by definition; if `f₁²>0` *with isolation* `s₂−s₁≥δ`, then `s₁` is an isolated eigenvalue above `0`, hence a gap. `[T]`
- **Floor ⟹ Susceptibility-satisfied `((G) ⟹ H1)`.** `s₁≥Δ²` with `s≥s₁` on support ⇒ `χ=∫s⁻¹dν ≤ M₀/Δ²` ⇒ `G″=1/χ ≥ Δ²/M₀`. `[T]` (`H1-IS-IT-THE-GAP §3(c)`).
- **Converse 1 fails: `H1 ⇏ (G)`.** Thin-atom witness `ν = wδ_ε + (1−w)δ_S`: `χ ≤ 1/λ` (H1 holds) but `s₁=ε→0`. `[T]` (`REDUCTION-AUDIT` (b)).
- **Converse 2 fails: `(G) ⇏ H*`.** `:φ²:`: floor `4m²>0`, atom `0`. `[T]` (`PROOF-F1SQ`, `§0.3`).

> **Partial order (as certifiers of Clay):** `H1 ⊊ (G) ⊊ H*` in *strength of hypothesis*: `H*` implies the gap but is unnecessary; `H1` is necessary but insufficient; `(G)` **is** the gap. **The target is `(G)`, the floor — the unique element that is both necessary-and-sufficient *and* passable by the atom-free witness.**

### §1.3 The `:φ²:` acid is the decisive discriminator `[T]`

Of the three, only the floor survives `:φ²:`:
- **susceptibility fails**: `G″_{:φ²:}=1/χ_{:φ²:}` is finite and equal to neither `4m²` nor `0`; it certifies nothing about the gap's *value* and is compatible with any `s₁` by the squeeze (`0 ≤ 4m² ≤ ∞`). `[T]`
- **atom fails**: `f₁²=0` reads "no atom," i.e. the atom functional *mis-classifies* the gapped `:φ²:` as gapless. `[T]`
- **floor passes**: `s₁=4m²>0` correctly, atom-free. `[T]`

Hence the `:φ²:` test — *gapped, atomless* — **uniquely selects the floor**. This is the constructive form of `SPECIES-EXISTENCE` Theorem 1 on the *genuine target*: a gap-*discriminator* must fire on `μ_m`, an atom-*minorant* must read `0` on `μ_m`, so the correct target's certifier and the atom-minorant are **disjoint**.

### §1.4 Why the susceptibility "gags" (the ceiling, made exact) `[T]`

The task's phrasing "`G″` only gags" is the ceiling `s₁ ≤ M₀ G″` (Identity S, §0.2): a large curvature places an **upper** bound on the floor. The naive "`1/χ ≤ s₁`" is a **sign error** — it would require harmonic-mean ≤ min, false (`H1-IS-IT-THE-GAP §3(a)`). The only *lower* bound through `χ` is `s₁ ≥ f₁²/χ = f₁² G″`, which needs the atom `f₁²≥r₀` and is then `H*`, not `H1`. **So `H1` is a necessary face and a ceiling, never a floor.** `[T]`

### §1.5 Why the atom is the wrong over-target (three independent reasons) `[T]/[J]`

1. **Unnecessary** (Clay needs support, not visibility): `(G) ⇏ H*` via `:φ²:`. `[T]`
2. **Unstable / upper-semicontinuous**: `f₁²` is a jump of `F` at `s₁`; under weak convergence an atom can smear to `0` (`SUPPORT-NOT-ATOM §4`, `PROOF-TRANSPORT-ANALYTIC K1`) — it can *drop to 0 at the Savvidy crossing*. The floor `s₁` (a rate) has no such drop; it is a `liminf`-type datum. `[T]/[J]`
3. **Harder**: `H*` = atom + isolation `I2` = "one rung above the gap" (`GAP-G2-edge`), strictly more than the floor. `[T]`

**Conclusion of §1.** The correct target is `s₁ = inf supp ν`, the **floor**, certified bottom-up (`𝔉 ≤ s₁`, uniformly positive). It is the **rate**, not the residue; the **support edge**, not the harmonic moment; the **middle** of the strict chain `H1 ⊊ (G) ⊊ H*`.

---

## §2 The threshold-edge structure the target must certify

### §2.1 What "certify the floor" means — the support/mass-function form `[T]`

Since `s₁ = inf{s : F(s)>0} = inf supp ν`, the statement `s₁ ≥ Δ²` is **exactly**
```
sup{ σ : ν([0,σ)) = 0 } ≥ Δ²     ⟺     F(Δ²⁻) = 0     ("no spectral mass below Δ²").     (SUPPORT FLOOR)
```
This is a statement about the **mass function `F` on the window `[0,Δ²)`**, not about the point density `ρ` nor the point atom `f₁²`. It is an **evaluation of `ν` on sets** (the low window), i.e. a **support/continuum-edge** readout. A certifier `𝔉_floor` must be a functional of `F` (equivalently of "no sub-Δ mass"), with `𝔉_floor(μ) ≤ s₁(μ)` and uniform positivity. `[T]`

### §2.2 The threshold requirement is FORCED by the atom-free acid `[T]`

**Theorem F (floor-certifier ⊥ atom-minorant on `:φ²:`). `[T]` `[V-here-assembly of SPECIES-EXISTENCE Thm 1].`** Let `𝔉` be a functional with `𝔉(μ) ≤ s₁(μ)` (genuine floor lower bound) that is uniformly positive on the gapped family, so it must read `𝔉(μ_m) > 0` on the atom-free gapped `:φ²:`. Any atom-minorant `𝔊 ≤ f₁²` reads `𝔊(μ_m) ≤ f₁²(μ_m) = 0`. Hence on `{:φ²:}`, `𝔉(μ_m)>0 > 𝔊(μ_m)`, so **`𝔉 ≠ 𝔊`; no single functional is both a floor-certifier and an atom-minorant.** ∎

> **Consequence.** The target **must be a `:φ²:`-gapped-detector**: it must fire on the **continuum edge** `[4m²,∞)` where there is **no atom**. In operator language it must certify `inf(σ(H)\{0}) ≥ Δ` **including the case where this infimum is a two-particle (or N-particle) threshold, not a bound-state pole.** The target is therefore a **threshold / continuum-edge functional**, and this is **not a choice** — it is forced by the one acid the atom cannot pass.

### §2.3 The threshold / continuum-edge structure (pole ⊕ cut) `[T]/[O]`

Decompose the floor `s₁ = min(s₁^pole, s₁^cut)` where
- **pole edge** `s₁^pole`: an isolated eigenvalue below the continuum (a genuine `0⁺⁺` bound state). Here `f₁²>0` (the atom route `H*` applies — *sufficient*), and the floor is the pole location.
- **cut edge / threshold** `s₁^cut`: the bottom of the continuum (a multi-particle threshold), where `f₁²=0` possibly (`:φ²:` is purely this). By **HVZ**, `σ_ess(H) = [Σ,∞)` with `Σ` the infimum of **thresholds built from the spectra of sub-clusters** — for the two-body case `Σ = 2 m_{lightest}` (sum of the lightest one-particle masses admitted by the channel).

The **target functional must be a lower bound on this `min`** — hence on the *edge of the essential spectrum* when there is no pole, and on the pole when there is one. A certifier that handles **only** poles is *not* a floor-certifier (fails `:φ²:`); a certifier that handles the threshold/cut is. **What this changes, precisely:**

1. **The readout shifts from a point to a set-boundary.** Atom = `ν` at a point (jump of `F`); floor = left endpoint of `{F>0}` (support boundary). The certifier is a **mass-function** functional, not a density/residue functional. `[T]`
2. **The mechanism shifts from residue to rate.** Because the atom is unavailable on the cut, the only atom-blind certificate is a **lower bound on the exponential decay rate** (Identity R), i.e. an **anti-concentration / clustering-rate** bound of the coarse modes — which is exactly the `H1` shape (`SUPPORT-NOT-ATOM §3`, `E2-FLIP-NOGO (L-2)`). So the floor route **reduces the target to a uniform rate floor**, not to a residue. `[T]/[O]`
3. **Cuts need edge handling.** A pole is a delta; a cut is `ρ(s)~(s−s₁)^θ ds`. The rate (R) is insensitive to the edge exponent `θ` (`§0.3`: rate `→2m` regardless of `τ^{−3/2}`), so the *rate* certifier is robust — but any certifier built instead from the density/residue near the edge inherits the edge-exponent analytic difficulty. This is why the **rate** is the correct analytic object. `[T-Vhere]`

### §2.4 The threshold route's circularity — HVZ is self-referential in YM `[O]` (inherited, `SUPPORT-NOT-ATOM §4`)

The natural threshold bound `Σ = 2m_g` **presupposes `m_g>0`** — and in YM the lightest one-particle state (the `0⁺⁺` glueball) **is** the gap itself, so `Σ = 2m_g = 2·(the answer)` — **circular**. HVZ also does **not** lower-bound the lightest *bound* state: a discrete eigenvalue may sit at `ε→0` below the threshold (exactly the case a pole below the cut). Hence: **the threshold structure makes the target atom-free *in principle* (good), but its canonical HVZ mechanism is circular for YM**, and the only non-circular floor certificate must be a **direct rate/anti-concentration bound** that does **not** route through `2m_g` — or a **mass-ordered bootstrap** (§3-(iv)) in which lighter sectors are certified first and a heavier channel's floor is then bounded by `threshold ≥ Σ(lighter masses)`. `[O]`

---

## §3 The requirement spec (what this facet imposes; the residual freedom)

### §3.1 The target functional — precise spec

> **TARGET SPEC (floor).** A Clay proof must exhibit a functional `𝔉_floor : F_adm → ℝ` with:
> - **(T1) genuine lower bound** `[T]`: `𝔉_floor(μ) ≤ s₁(μ) = inf supp ν_μ` for all `μ ∈ F_adm` (equivalently `𝔉_floor ≤ −lim_{τ→∞}(1/τ)\log C_μ(τ)`).
> - **(T2) uniform positivity** `[T]`: `inf_{μ ∈ F_YM} 𝔉_floor(μ) ≥ Δ² > 0`, uniform in `(L,a)`, at the crossover `β≈2.5`.
> - **(T3) atom-blind / threshold-firing** `[T]`: `𝔉_floor(μ_m) > 0` on the gapped atom-free `:φ²:` (it reads the continuum/threshold edge, `F(Δ²⁻)=0`), and it is a functional of the **mass function `F(s)=ν([0,s])`**, not of the point atom.
> - **(T4) acid** `[T]` (implied by T1, `E2-FLIP-NOGO Thm 0`): `𝔉_floor(μ₀) ≤ 0` on free massless (`s₁=0`).
>
> The object to be *bounded* is `s₁ = inf supp ν` = the **exponential decay rate** = the **lightest `0⁺⁺` mass** = `inf(σ(H)\{0})` in the channel. The certificate is `𝔉_floor ≤ s₁`, uniformly positive. **Nothing about the residue `f₁²` or the curvature `G″=1/χ` is part of the target.**

**This is now a complete, unambiguous target.** The facet's job is to *fix the target*, and it is fixed: **the floor, certified atom-blind, as a rate.**

### §3.2 The residual degree of freedom (what the spec deliberately leaves open)

The target is pinned; the following are the **genuine DOF** a proof may exercise — this is where freedom lives, not in the target itself:

1. **Channel selection `[J]`.** `s₁` here is the `O=¼trF²` (`0⁺⁺`) floor; the *full-spectrum* Clay needs the minimum over all channels. The reduction "the lightest glueball is `0⁺⁺`" is a `[J]`/kinematic input (empirical; not proven here). **Freedom:** which interpolator/channel carries the global floor; the target still applies channel-by-channel, but the *global* `inf` requires the kinematic identification.
2. **Atom-route vs threshold-route `[T]/[C]`.** The certifier may prove the floor by either (a) an **isolated edge pole** below threshold (`H*` — sufficient, atom route), or (b) a **continuum edge / rate lower bound** (atom-free, the `:φ²:`-compatible route). **Clay needs only one.** For the lightest YM glueball (a bound state `m_g < 2m_g`), the atom route is *available* and sufficient; for an atom-free edge the rate route is *required*. **Freedom:** pick the cheaper certificate per channel — one is always a valid floor bound.
3. **The uniform rate constant `[O]`.** The *only* genuinely open content is T2: a **uniform-in-`(L,a)` lower bound on the exponential rate** across the Savvidy crossover. This is precisely the residual object `H1` / the transport machine (`SUPPORT-NOT-ATOM §3`, `E2-FLIP-NOGO (L-2)`): *the target is pinned, the uniform rate floor is not*. **This is the exact missing step — not a target ambiguity.** `[O]`
4. **HVZ bootstrap ordering `[O]`.** To make the threshold route non-circular, certify **lighter sectors first** (their floors are proven bottom-up), then bound a heavier channel by `threshold ≥ Σ(lighter one-particle masses)`. **Freedom/obligation:** choose a sector ordering such that each threshold is built from *already-certified* lighter masses, terminating at the single lightest floor (which must be certified directly as a rate, not via its own threshold). Whether such an ordering closes *all* channels from one direct lightest-floor certificate is the open structural question. `[O]`
5. **Crossing/transport stability `[O]`.** The floor (a rate/liminf) is **more stable than the atom** (a jump/usc) under the Savvidy crossing, but still needs a crossing argument: a uniform rate floor must be transported across `β≈2.5`. This is the same transport machine (`CLAUDE.md §5`), now cast on the **rate** (robust) instead of the atom (fragile) — a genuine *simplification of the target's fragility*, not of its difficulty. `[O]`

### §3.3 What the facet does and does not claim (`[J]`, no QED)

- **Does** (this facet): fix the **target** unambiguously — `s₁ = inf supp ν`, the floor, certified atom-blind & threshold-firing as a **rate**; prove the **strict hierarchy** `H1 ⊊ (G) ⊊ H*` with witnesses; prove **floor-certifier ⊥ atom-minorant** on `:φ²:` (Theorem F); exhibit the **rate=floor / prefactor=atom** separation by exact edge-Laplace computation (`§0.3`); and enumerate the **residual DOF** (§3.2).
- **Does not**: prove Clay; prove the uniform rate floor (T2 = the open `[O]`); prove the channel-kinematics `[J]`; or claim the transport machine exists. **The target is now precise; the open object (uniform rate ≥ Δ) is unmoved.** `[J]`

---

## References

**Internal (banked; labels carried):**
- `D:\Dev\math\research\unconditional\H1-IS-IT-THE-GAP.md` — the squeeze `f₁²/χ ≤ s₁ ≤ M₀/χ`; `G″=1/χ=V_eff″`; `:φ²:` gapped atom-free (§2.3-iii); sign of the harmonic mean (§3(a)); H1 necessary-not-sufficient (§3).
- `D:\Dev\math\research\unconditional\SUPPORT-NOT-ATOM.md` — support vs atom vs visibility (§0); Bushell/Krein–Rutman rate `=G″=H1` (§2–3); HVZ `Σ=2m_g` circular, no lightest-bound-state bound, Weyl blind (§4).
- `D:\Dev\math\research\unconditional\SPECIES-EXISTENCE.md` — Thm 1 discriminator ≠ atom-certifier (§1, the `:φ²:` witness); Thm 2 `F3∧F4=∅`; Thm 3 circularity; §0.1 support (G) vs atom (H*).
- `D:\Dev\math\research\unconditional\E2-FLIP-NOGO.md` — the flip functional `Φ≤m₁` (≡ floor certifier); Thm 0 acid implied by (LB); four-horn dichotomy; (L-2) field-space inequality ≡ H1.
- `D:\Dev\math\research\unconditional\GAP-G2-edge.md` — residue–level-gap lemma; edge = one rung above the gap.
- `D:\Dev\math\research\unconditional\PROOF-F1SQ.md` — `ν`, `M_k`, `:φ²:` `f₁²=0`.
- `D:\Dev\math\docs\new-tools-registry.md` — §A4 gauge-invariant IR mass = `trF²` spectral edge (free Maxwell: threshold-0 two-photon cut); §D residual dictionary D.1–D.7.
- `D:\Dev\math\docs\PROOF-SKELETON.md` — dimension `G″_a=a²m_phys²`; target = ratio `G″_a/a²≥Δ₀²`.

**External `[R]` (standard; not re-verified this session):** Källén (1952)/Lehmann (1954) — KL representation; Osterwalder–Schrader (1973)/Osterwalder–Seiler, *Ann. Phys.* **110** (1978) — OS/RP; Jaffe, *Phys. Rev.* **158** (1967) — generalized free fields; Bernstein–Widder — completely monotone ⟺ `ν≥0`; Hunziker (1968)/van Winter (1964)/Zhislin (1960) — HVZ essential-spectrum theorem; Widder *Laplace Transform* (1941) — edge Laplace asymptotics.

---

*Status: **PROFILE FACET — TARGET SPEC.** The correct Clay target is the **spectrum floor** `s₁=inf supp ν` (an exponential-rate / support-edge functional, atom-blind, threshold-firing), certified by `𝔉_floor≤s₁` uniformly positive. It is rigorously separated from the atom `f₁²` (point/prefactor, usc, strictly stronger, wrong over-target; fails `:φ²:`) and the susceptibility `G″=1/χ` (harmonic moment, only ceilings `s₁≤M₀/χ`; fails `:φ²:`). The strict chain `H1⊊(G)⊊H*` is proved with witnesses; floor-certifier ⊥ atom-minorant is proved on `:φ²:` (Theorem F); the rate=floor / prefactor=atom split is computed exactly by edge-Laplace asymptotics (`[T-Vhere]`, numerically confirmed). The target is fixed; the residual DOF are channel-kinematics `[J]`, atom-vs-threshold route, the uniform rate floor `[O]=H1`, and the HVZ bootstrap ordering `[O]`. Clay NOT proved; the uniform rate floor is unmoved and open; 0%. No QED on anything unproven.*

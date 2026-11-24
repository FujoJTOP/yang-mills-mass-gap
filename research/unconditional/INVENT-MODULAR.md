# INVENT-MODULAR — the Tomita–Takesaki modular spectrum is NOT a new species: it is the transfer/thermal spectrum in disguise, or trivial, or the entanglement spectrum (a gap-detector, not an atom-certifier)

> **VERDICT (no QED).** The candidate "new algebraic species" — the `0⁺⁺`-channel **modular** spectral
> quantity (smallest modular eigenvalue above the vacuum / "modular gap") of the Tomita–Takesaki operator
> `Δ_μ = e^{-K}` — **fails**. It is provably one of four things, none of which is a new handle on the
> `0⁺⁺` **atom**: **(T1) trivial** (`Δ=1` for the full-algebra pure vacuum); **(T2) the transfer/thermal
> spectrum** (`KMS ⇒ K=βH`, so the *modular gap = β·(mass gap) verbatim — circular*); **(T3) the kinematic
> Bisognano–Wichmann boost** (subregion wedge; state-independent continuous spectrum, no mass gap); or
> **(T4) the entanglement spectrum** — which is a **gap-detector** (positive for *any* gapped measure,
> including the atomless `:φ²:` where `f₁²=0`), so it **fails the discriminating free-field acid test**
> and is circular for the gap. The deep reason: the one gap-robust modular invariant, Connes' **type
> `III₁`**, is **universal** for continuum thermal states — **gap-agnostic**. The modular operator is the
> *wrong invariant*. **`[NG]` for the new-species class; `[O]` exact failing step = the modular gap is a
> gap-detector, not an atom-certifier.** Clay 0%.

---

## §0 · Labels, iron laws, acid test (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**IL-i** non-positivity (must fail on free massless: `f₁²=0`, `G″=0`; and discriminatingly on `:φ²:`: `f₁²=0`, gap`>0`). **IL-ii** nonabelian. **IL-iii** no small parameter (`α_s~O(1)`).

---

## §1 · The definition (what the species would be)

Fix the gauge-invariant algebra `M` (the algebra of observables) and a state `ω` (the physical measure `μ`,
represented by a cyclic separating vector `Ω` in the GNS representation). Tomita–Takesaki gives the
**modular operator** `Δ_ω` (with `Δ_ω = S_ω*S_ω`, `S_ω` the closure of `AΩ ↦ A*Ω`), the **modular
conjugation** `J_ω`, and the **modular Hamiltonian** `K_ω = −log Δ_ω`, generating the modular flow
`σ_t = Ad Δ_ω^{it}`.

The proposed `0⁺⁺`-channel species is the **modular gap**
```
Δ_mod(ω) :=  (smallest modular eigenvalue of K_ω above the vacuum)  =  the gap of the 0⁺⁺ block of K_ω,
```
or, for a subregion `A`, the corresponding **entanglement gap** of the modular Hamiltonian `K_A = −log ρ_A`.
The hope: `Δ_mod` is *inhomogeneous* in the coarse mass (`κ = F(κ) + 𝒮` with `𝒮` a genuine source), or at
least *phase-sensitive* at `β=2.5`, and bounds `χ_lattice ≲ a²` / `σ_w(2.5)>c₀` / `R_vac`.

**It is none of these.** Four exact faces exhaust it.

---

## §2 · The four faces (a trichotomy-plus-one theorem)

**(T1) Full-algebra pure vacuum ⇒ trivial `[T]`.** On the full algebra `B(H)` with a pure state
`ω(A)=⟨Ω|A|Ω⟩` (`Ω` cyclic), the GNS modular data is `J = `(Ω-basis conjugation), **`Δ = 1`**. Hence
`K_ω = 0`: the modular spectrum is `{0}`, no modular gap, no information. (The zero-temperature lattice
gauge vacuum is exactly this pure state on the full algebra.)

**(T2) KMS/thermal state ⇒ the transfer spectrum `[T]`.** For a KMS state at inverse temperature `β_th`
(equivalently any Gibbs state `ρ = e^{−βH}/Z` on `B(H)`, the lattice's finite-`T` / transfer-matrix setting),
Tomita–Takesaki gives `Δ = e^{−βH}` (up to a central constant), i.e.
```
K_ω = −log ρ = β H + const·1   ⇒   modular spectrum = { β E_n }  = the Boltzmann spectrum,
⇒   Δ_mod = β·(E₁ − E₀) = β·(mass gap)  —  the mass gap VERBATIM, rescaled.
```
**The modular gap is a function of the mass gap: circular.** (Witness PART A/B: `Δ_mod = β·Δ` at every
`β`; the full modular spectrum equals `β×(Hamiltonian spectrum)`, no new invariant.)

**(T3) Subregion wedge ⇒ kinematic boost `[T]`.** For the vacuum of a Lorentz-invariant theory and the
algebra of a Rindler wedge `W`, **Bisognano–Wichmann** gives `Δ_W = e^{−2π K_boost}`: the modular Hamiltonian
is the *boost generator*, whose spectrum is **continuous** (no mass gap) and **state-independent** (depends
only on the wedge and Lorentz invariance). No gap, no phase information.

**(T4) General subregion ⇒ entanglement spectrum (gap-detector) `[T]/[V-num]`.** For a general region `A` in
the vacuum, `K_A = −log ρ_A` is the **entanglement Hamiltonian**; its gap is the **entanglement gap**. For a
gapped system `Δ_ent ≳ 1/ξ ~ mass gap > 0`; for a gapless system the spectrum is dense and `Δ_ent → 0`. It
**tracks the gap** — circular — and, decisively, it is **positive for any gapped measure including the
atomless `:φ²:` (`f₁²=0`, gap`>0`)**, so it **fails the discriminating acid test** (§3).

Together: the modular spectrum is *trivial* (full-algebra vacuum), *the gap rescaled* (KMS/thermal), or *a
gap-detector* (subregion/entanglement). **No face carries the atom.**

---

## §3 · Acid test (the decisive failure)

Two-adversary acid (`IL-i`): must vanish on the **free massless** field (`f₁²=0`, gap`=0`) **and** be `≤0`
(ideally vanish) on the **gapped-atomless `:φ²:`** (`f₁²=0`, gap`>0`).

- Free massless: `Δ_mod = β·0 = 0` ✓ (the acid *returns* 0 — necessary).
- **Gapped atomless `:φ²:`: `Δ_mod = β·Δ > 0` ✗.** The modular gap **fires where `f₁² = 0`**.

Witness (PART A): at `β=2.5`, the atomless sector (`gap=0.6`, `f₁²=0`) gives `Δ_mod = 1.5 > 0`, while the
free massless sector gives `0`. So `Δ_mod` is a **gap-detector**, not an **atom-certifier** — precisely the
already-established "gap ⇏ atom" failure (`:φ²:` is gapped, atomless). The modular gap cannot separate them.

---

## §4 · Phase-sensitivity (ii) and the bound (iii)

- **(ii) Phase-sensitive at `β≈2.5`?** By (T2)/(T4) the modular gap is a *function of the mass gap*, hence a
  smooth positive function through the crossover. It is "phase-sensitive" only in the sense that *it is the
  gap* — circular; it adds no independent phase information beyond the object it is supposed to certify.
- **(iii) Bounds `χ_lattice ≲ a²` / `σ_w(2.5)>c₀` / `R_vac`?** No intrinsic relation. The modular gap is a
  spectral quantity; any inequality to `χ`/`σ`/`R_vac` would be a *new derivation*, not an identity, and
  since `Δ_mod` is itself the gap in disguise, such a bound is circular. There is no free lunch.

---

## §5 · The deep reason: Connes type is gap-agnostic

The one modular invariant that is *not* just the spectrum is the **type** of the von Neumann algebra
(Connes' classification). In algebraic QFT, the local algebra of any thermal (KMS) state in `≥3` spacetime
dimensions is **type `III₁`** — *universally*, independent of whether the theory is gapped or gapless, and
independent of the coupling. On the finite lattice the algebra is **type `I`** (finite-dimensional) and `Δ`
is literally the reduced density matrix (the entanglement spectrum). In **both** cases the modular *type*
**does not distinguish gapped from gapless**. And the modular *spectrum* is the Boltzmann/transfer spectrum
(T2) or the entanglement spectrum (T4). So the modular structure is **the wrong invariant** for the atom:
its only non-spectral invariant (type) is universal, and its spectral content is the gap itself.

---

## §6 · The one honest "new" feature — and why it does not break the wall

The modular operator is genuinely **state-dependent** (unlike the Hamiltonian), and the **relative modular
operator / relative entropy** `S(ω‖φ) = −⟨log Δ_{ω|φ}⟩` is a bona fide *state-discriminator*: it is finite
for "close" states and distinguishes phases. **But**: (a) relative entropy is a **positivity-clean
(lower-semicontinuous) functional** — it cannot *print* an atom (the project's "atom is usc, positivity is
lsc" obstruction applies verbatim); and (b) extracting a rigorous *uniform* mass-gap bound from
`S(ω_true‖Gaussian)` is the same crossover/transport wall in yet another coordinate (it is an `H1`-type
quantitative statement). So the state-sensitivity is real but *vacuous for the atom* and *circular for the
gap* — the same two ways every previous route died.

---

## §7 · Verdict and exact failing step

- **Not a new species.** The modular spectrum is (T1) trivial, (T2) the transfer/thermal spectrum (`K=βH`;
  modular gap `=β·Δ`, circular), (T3) the kinematic Bisognano–Wichmann boost, or (T4) the entanglement
  spectrum (gap-detector). None carries the `0⁺⁺` atom.
- **Exact failing step `[O]`:** the modular gap `Δ_mod = β·(mass gap)` is a **gap-detector, not an
  atom-certifier** — it is *positive for the gapped-atomless `:φ²:` (`f₁²=0`)*, so it fails the
  discriminating acid test, and for the KMS state it *is* the mass gap exactly (circular). The one
  non-spectral modular invariant (Connes type `III₁`) is universal and gap-agnostic.
- **Reduction:** the modular operator reduces to the transfer matrix (KMS), to the identity (pure vacuum),
  or to the entanglement Hamiltonian (subregion) — i.e. **to the gap itself, circularly**. It is **not** a
  sixth species.

**Clay 0%; no QED.** The modular/Algebraic route is closed; the search for the strut (a renormalon-free
`0⁺⁺` anomaly / dimensional transmutation, `=H1`) continues along *other* invented species.

---

*Refs:* `ROOF-SYNTHESIS.md` (§4: the strut must be a non-order-parameter, non-Markov, phase-sensitive
functional on `μ`), `BASE-SYNTHESIS.md` (P1/P4), `SLOT-SYNTHESIS.md` (atom usc vs lsc), `PROFILE.md`.
Relevant standard theory: Tomita–Takesaki (1970s), Bisognano–Wichmann (1975–76), Connes' type
classification (`III₁` universality), Casini–Huerta–Myers (entanglement Hamiltonian ↔ modular Hamiltonian).
Witness: `w7/invent_modular.py` (runs clean, `py -3.12 -I`, exact, `<0.1 s`).

# WI-THEOREM — weight-invariance of `𝒜=π/|Im τ_zero|` characterises an *isolated* atom (I2), not merely its weight (I1)

> **VERDICT (no QED).** For `Z(τ)=w₀+Σ_i w_i e^{−E_i τ}` (`w₀>0`, `w_i>0`, `0<E₁≤E₂≤…`) and
> `𝒜:=π/min_k|Im τ_k|` over zeros, the map **weights ↦ 𝒜 is constant iff the spectral support above
> the VEV is a *single (isolated) atom*** `E₁` — and then `𝒜=E₁` exactly. A second level or a continuum
> makes `𝒜` weight-dependent, with sensitivity a **linear functional of the cut density**. In KL language
> this is the statement that **the WI test detects isolation (I2 / `H*`), not merely the atom weight
> (`f₁²>0`, I1)** — the harder half of the gap. **Clay 0%.**

---

## §0 · Setting and labels

`Z(τ)=w₀+Σ_i w_i e^{−E_i τ}`, `w₀=|⟨0|O|0⟩|²>0`, `w_i=|⟨i|O|0⟩|²>0`, energies `0<E₁≤E₂≤…` fixed;
zeros `{τ_k}` of `Z` in `ℂ`; `𝒜:=π/min_k|Im τ_k|` (the `PNT-IR-ANALYTIC` functional). **Reweighting** =
`(w₀,{w_i})↦(w₀',{w_i'})` at *fixed* energies (an operator choice `O→O_λ` changes overlaps `w_i`, not the
Hamiltonian spectrum). Physical reweightings `O→λ₁O+λ₂O'` give `w_i=|λ₁a_i+λ₂b_i|²`, a subset of arbitrary
positive reweighting, so theorems for all positive reweightings apply.

Labels **[T]** exact · **[C]** conditional/numeric · **[V-num]** verified numerically · **[O]** open.

---

## §1 · Theorem A (sufficiency, `[T]`): a single atom ⇒ `𝒜` invariant, `=E₁`

If the support is one atom, `Z(τ)=w₀+w₁e^{−E₁τ}`. Then `Z(τ)=0 ⟺ e^{−E₁τ}=−w₀/w₁`. Writing `τ=σ+iθ`:
`|e^{−E₁τ}|=e^{−E₁σ}=w₀/w₁ ⇒ σ=ln(w₁/w₀)/E₁`, and `e^{−iE₁θ}=−1 ⇒ E₁θ=(2k+1)π`. Hence
```
τ_k = ( ln(w₁/w₀) − i(2k+1)π ) / E₁ ,     min_k |Im τ_k| = π/E₁ ,
𝒜 = π/(π/E₁) = E₁ .
```
**The weight `w₁` and the VEV `w₀` only move `Re τ = ln(w₁/w₀)/E₁`; the modulus `w₀/w₁` sits on the real
axis and never touches `Im τ`.** So `𝒜=E₁` for every `(w₀,w₁)>0`. **Verified** to machine precision
(`w₁∈{0.1,0.3,1,3,5}`, `w₀∈{0.5,1,2}`): `𝒜−E₁=0.0`.

---

## §2 · Theorem B (non-invariance, `[T]`/`[V-num]`): a second level / a continuum makes `𝒜` move

**First-order (generic).** At the atom-only zero `τ_a=(ln(w₁/w₀)−iπ)/E₁` with `Z′_atom(τ_a)=E₁w₀`, adding a
level `w₂e^{−E₂τ}` shifts the zero by `δτ=−w₂e^{−E₂τ_a}/(E₁w₀)`, so
```
δ(Im τ) = w₂ · e^{−E₂ Re τ_a} · sin(E₂π/E₁) / (E₁ w₀) .
```
This vanishes only if `E₂/E₁∈ℤ`; generically `≠0`, so `𝒜` moves with `w₂`. In the integer-ratio case the
first-order term cancels but `Z` is a polynomial in `u=e^{−E₁τ}` whose roots (hence `arg u`, hence `Im τ`)
still depend on the weights — **non-invariant for every `N≥2`**.

**Verified (`E₁=1, E₂=√2`):** `𝒜 = 1.000307, 1.003062, 1.029840, 1.121243, 1.185802, 1.281846` for
`w₂=10⁻³,10⁻²,0.1,0.5,1,3` — a spread of `0.2815` at fixed energies. **A continuum** (discretised `ρ` on
`[E₁,∞)`) likewise moves `𝒜` with the weight scale (and the zero set *crowds* — the `count=21` vs sparse
`count=5` the numeric agent observed). **→ `𝒜` is NOT reweighting-invariant without isolation.**

---

## §3 · Proposition C (sensitivity, `[C]`/`[V-num]`): the derivative is a linear functional of the cut density

Add a gapped continuum `ε·ρ(E)` above `E₁+Δ`. The nearest zero shifts by `δτ=−ε p(τ_a)/(E₁w₀)` with
`p(τ_a)=∫ρ(E)e^{−Eτ_a}dE`, hence
```
d𝒜/dε = E₁ · Im[p(τ_a)] / (π w₀) ,     p(τ_a)=∫_Δ^∞ ρ(E) e^{−E τ_a} dE ,
```
a **linear functional of the cut density `ρ`**. **Verified** (`Δ=0.5`, `ρ=1`, `ε≤0.008`): finite-difference
`d𝒜/dε=+0.016243` vs predicted `+0.016324` (relative agreement `5×10⁻³`). So the *sensitivity* of `𝒜` to
reweighting is exactly the continuum content — the diagnostic side of the WI test.

---

## §4 · The characterisation, and why it is I2 not I1

Combining §1–§3:
```
𝒜 is reweighting-invariant  ⟺  supp(spectral measure above w₀) is a SINGLE point E₁  (𝒜=E₁).
```
- **A single atom** (support `{E₁}`) ⇒ invariant (`𝒜=E₁`), §1.
- **Any second level or continuum** (even with a **gap** `Δ>0`) ⇒ non-invariant, §2. **A gap alone does
  not restore invariance** — the second level's `sin(E₂π/E₁)` shift is nonzero regardless of `Δ`.

In KL language for the `0⁺⁺` two-point function (`w₀=|⟨tr F²⟩|²`, atom `f₁²` at `m₁²`), "support is a
single isolated atom" means **`f₁²>0` *and* a gap to the continuum** — i.e. exactly the **isolation `I2`**
(`H*`), the harder half of the gap. **The atom weight alone (`f₁²>0`, I1) is *not* enough**: an atom with a
*gapless* continuum at the floor is non-invariant. This is the precise content of the numeric agent's
"weight-invariance is a finite-L atom signature": it is an **I2 (isolation) signature**, sharper than I1.

---

## §5 · Scope, and the honest limit

- **What is proved:** a clean characterisation of *when* the `𝒜` functional is operator-reweighting-stable
  (Theorem A exact; Theorem B by first-order + integer-ratio arguments + witness; Proposition C numeric).
  This is a statement about the **Laplace-transform zero structure**, finite or infinite volume, for any
  positive spectral weights.
- **What is NOT proved / `[O]`:** that the physical YM spectrum is a single isolated atom — that *is* the
  gap (`H1`/`I2`). At finite `L` the raw `𝒜` fires on adversaries (`construct_ir_analytic_num`); the
  reweighting-invariance is the finite-`L`-meaningful invariant, but proving it **uniformly in `(L,a)` at
  `β≈2.5`** (with reweighting that preserves the coupling content) is `H1`. The WI test is a **handle for
  I2, not a proof of it.**
- **Caveat:** requiring `w_i>0` strictly; a zero weight just shrinks the effective support. Physical
  reweightings are a subset of arbitrary positive reweightings, so the results apply a fortiori.

**Clay 0%; no QED.** *Files:* `w7/wi_theorem_witness.py` (verified: §1 `𝒜−E₁=0`, §2 spread `0.2815`,
§3 agreement `5×10⁻³`).

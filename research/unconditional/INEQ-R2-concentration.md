# INEQ-R2 — is there an INDEPENDENT spectral-concentration inequality for the tail, or does it collapse to the level gap E₂?

> 2026-11-23 · Attack on `f₁²=|⟨0|O|G⟩|²` via the variance decomposition. **Verdict: NO independent concentration inequality exists. `f₁²=Var−tail` is an identity, so a *useful* lower bound `f₁²≥Var−tail` requires `tail≤(1−c)Var` with `c>0`; that tail bound is provably FALSE from positivity/moments/operator-norm alone, and is *equivalent* to the residue–level-gap lemma with `Λ=E₂`. The route collapses to GAP-G2-edge.** No fabrication.

**Setup.** RP gives `ν≥0`, `supp ν⊂[m₁²,∞)`, with `Var_μ(O)=M₀=∫dν=Σ_{n≥1}|⟨0|O|n⟩|²`, atom `f₁²=ν({m₁²})`, `tail=M₀−f₁²=Σ_{n≥2}|⟨0|O|n⟩|²`. "Spectral concentration" = `f₁²≥c·M₀`, i.e. the lightest carries a definite fraction `c>0`. This is *precisely* "atom relative to total," the object of INEQ-T1.

**Positivity + moments cannot force a fraction (explicit).** For any `ε∈(0,1)` take `ν_ε=εM₀δ_{m₁²}+(1−ε)ρ`, `ρ` a density on `[m₁²+δ,∞)`. Then `f₁²=εM₀`, `tail=(1−ε)M₀`, all positivity intact, support law intact. As `δ→0` the two "levels" merge; `inf_ν f₁²=0` over the moment fiber (INEQ-T1: identical `(M₀,M₁)` realised by atom `1` and atom `0`). Operator-norm gives only `Var≤‖O−⟨O⟩‖²` — an **upper** bound on the total, never a lower bound on one summand. Matrix-valued/cross-correlator positivity is still PSD ⇒ same upper-only sign (EDGE-f₁² sign obstruction). Dimensionless `c` would need a moment *ratio*; moments give only the ratio lemma `M_a/M_b` (upper). Free field sits at `f₁²=0` (continuum, no atom). So **`c>0` is unprovable; any such claim must inject a second scale.**

**The tail bound *is* the level-gap lemma (equivalence).** Put the rest on `[E₂,∞)`, `E₂>m₁²` the second squared-mass threshold. On the tail `s≥E₂`:
`M₁=f₁²m₁²+∫_tail s dν ≥ f₁²m₁²+E₂·tail`. With `tail=M₀−f₁²`:
`M₁ ≥ E₂M₀+f₁²(m₁²−E₂)` ⟹ **`f₁² ≥ (E₂M₀−M₁)/(E₂−m₁²)`** — exactly GAP-G2-edge's `(ΛM₀−M₁)/(Λ−m₁²)`, `Λ=E₂`. A Chebyshev `tail≤M₁/E₂` gives the same at first order. So the concentration input is nothing but `Λ=E₂`.

**What it costs.** (i) a *strict* second level `E₂>m₁²` (the two-glueball threshold — a binding/dynamical fact); (ii) a bound on the unprotected total `M₁=⟨0|(O−⟨O⟩)H(O−⟨O⟩)|0⟩`. Both are absent unconditionally. The apparent "second-moment structure" is just `M₁`; it is unprotected and useless without `E₂`.

**Sharpest line.** `f₁²=Var−tail` is a tautology. Turning it into a bound needs mass *beneath* the tail, and "no mass between `m₁²` and `E₂`" is the gap one rung up. There is no concentration independent of `E₂`: concentration **is** `E₂`. Collapse to GAP-G2-edge; the edge theorem rewritten in variance variables.

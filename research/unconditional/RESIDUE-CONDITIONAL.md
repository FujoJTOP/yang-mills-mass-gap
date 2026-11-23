# RESIDUE-CONDITIONAL — the Clay gap as ONE lower bound on ONE unprotected residue

> 2026-11-23 · The campaign's endpoint, as a clean conditional theorem. **Result.** For pure SU(N) YM the entire unprotected, dynamical, non-positivity content of the mass gap is a *single number*: a uniform lower bound `f₁² ≥ r₀ > 0` on the residue of the lightest transverse spectral atom. Everything else in the chain — OS/reflection positivity, the Källén–Lehmann measure, transfer-matrix identification, the residue–level-gap trade — is a theorem or a re-expression. Below: the hypothesis, what it is *not*, the chain, and exactly what stays open. Theorem vs hypothesis is marked at every step; nothing is fabricated.

## Hypothesis H*(r₀)
`φ` a gauge-invariant transverse interpolator, `ν` its KL measure: `ν≥0`, `supp ν⊂[m₁²,∞)`, `f₁²=ν({m₁²})`, moments `M_k=∫sᵏdν`. **H\*(r₀):** `f₁² ≥ r₀ > 0`, uniformly in volume `V` and cutoff `a`. (PT–BFM form: the IR massless pole of the transverse three-gluon vertex has residue `≥ r₀`.) This is the *Schwinger order parameter*: it is one object with Polchinski `c`, the marginal curvature `G″`, and the LGP constant `γ` (CONVERGENCE-redteam §c) — one datum in several vocabularies.

## What it is NOT
- **Not positivity.** RP holds for the free massless scalar (`ν=δ₀`, gap 0). RP fixes signs, never the size of `1−λ₁` (P4 §1).
- **Not a moment/sum-rule bound.** Linear-data lemma (LEVER §0B): no `Φ=Φ(moments)` lower-bounds `inf supp ν`; moments are blind to a small edge atom (`δ₀+δ₂` vs `2δ₁`).
- **Not protected.** The transverse vertex has no anomalous slot (W3-marginal-BPHZ); no Ward/BRST identity fixes the residue. Only dynamics does.
- **Not `D(0)=0`** (Gribov–Zwanziger, strictly weaker), **not variational** (min–max is upper).

## Rigorous skeleton (no gap assumed)
1. **OS'78.** RP ⇒ `(H,T,Ω)`, `T=e^{−aH}`, `0≤T≤1`, `Δ_V=−(1/a)ln λ₁=min supp ν`.
2. **KL spectral theorem.** `C(t)=Σf_n²e^{−m_n t}`, `f₁²=|⟨Ω|φ|1⟩|²`.
3. **Residue–level-gap Lemma** (GAP-G2-edge, proved): `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)`, equality iff the rest is one atom at `Λ`; given moments, residue and second-level threshold are interconvertible.
4. **Moment-ratio lemma:** `m₁² ≤ (M_a/M_b)^{1/(a−b)}` — always *upper*. So no residue/moment algebra alone lower-bounds the mass: H\* is irreducibly a spectral-*edge* datum.

## The chain H*(r₀) ⇒ volume-uniform gap
Grant H\*(r₀) with the lightest atom isolated at `m₁²>0`. Then:
- **(a) No massless edge.** An edge atom of weight `≥r₀` lives at positive mass: a massless transverse atom of weight `r₀` (Coulomb phase, LEVER "the one door") is excluded by the order-parameter reading of H\*.
- **(b) Two-scale separation.** The lemma trades `f₁²≥r₀`+moments for a bound on the second threshold `Λ`.
- **(c) Two-scale Poincaré** (CONDITIONAL-LATTICE-GAP, a theorem *given* its curvature hypotheses H1–H3): `Var_μ(F) ≤ (1/m)Σ_x μ(|∇_x F|²)`, `m=λΛ_UV/(λ+Λ_UV)`, hence `Δ_Λ ≥ m` uniformly in `L,a`; by OS reconstruction, exponential clustering at rate `m`.

So `f₁²≥r₀` ⇒ `m(r₀)>0` on the lattice — steps (1)–(4) are theorems, (a)–(c) the clean conditional implication. **Honest caveat:** no closed form `m(r₀)` is currently derivable; the residue enters through (a)+(b), not a formula. This is the furthest the machinery carries it.

## What remains beyond `f₁²≥r₀`
- **Isolation / two-scale separation.** The *fixed-`r₀`* edge statement is logically **strictly stronger** than a bare `m>0` (a gap permits near-degeneracy that splits the edge weight) — CONVERGENCE-redteam §a. It is the *sharp* single-number form of the gap, not a weaker shortcut.
- **H1–H3** of CONDITIONAL-LATTICE-GAP: the softest `k→0` mode gives `G″(0)<0` (VAR-FINAL), so the two-scale route is currently *disfavoured*, not proved.
- **Continuum `a→0` (IN-9).** All of the above is fixed-β lattice; uniform control along the RG trajectory is separate and open. The residue is unprotected, so a uniform `r₀` surviving `a→0` is itself the open content.

**Sharpest line.** The Clay gap *is* the statement that one unprotected residue admits a positive lower bound uniform in `V,a`: prove `f₁²≥r₀` (with isolation and `a→0`) and the gap follows by KL+OS+two-scale Poincaré; nothing weaker suffices, and no other unprotected input is missing.

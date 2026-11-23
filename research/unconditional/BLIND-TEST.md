# BLIND TEST — 4D pure SU(N) Yang–Mills: Is There a Mass Gap?

*Independent, first-principles analysis. No project reduction, framework, or named conjecture is assumed.*

## 1. Definition

Work on Λ = (aℤ/Lℤ)⁴ with β = 2N/g₀² and the Wilson action
S = β Σ_p (1 − N⁻¹ Re tr U_p),
Haar measure on links. Let ⟨·⟩ be the Gibbs state and 𝒜 the gauge-invariant local algebra (traces of Wilson loops and of field-strength bilinears).

Reflection positivity of (S, Haar) gives a self-adjoint transfer matrix T = e^{−aH}, 0 ≤ T ≤ 1, on the physical Hilbert space H_gi (spatial links modulo Gauss law), with Perron–Frobenius vacuum Ω and λ₁ = 1. Define the lattice gap
Δ(a,β,L) := −log(λ₂/λ₁) ∈ (0,∞].

**Continuum mass gap.** Along the asymptotic-freedom curve a(β)Λ_lat → 0, the theory has a mass gap iff there exist aₙ → 0, Lₙ → ∞ with
m := lim Δ(aₙ,β(aₙ),Lₙ)/aₙ ∈ (0,∞),
and the limit is nontrivial (interacting, confining). The uniform form to be proved is: ∃ c > 0, a₀ with Δ(a,β(a),L) ≥ c·a for all a < a₀ and L > L₀(a). By the spectral theorem Δ is the exponential decay rate in Euclidean time of connected gauge-invariant correlators; locality upgrades this to spatial decay. The matching upper bound Δ = O(a) belongs to *existence*, not to the mass-gap content.

## 2. Three obstructions

(i) **Two scales.** For β < β₀ a convergent expansion works. At β → ∞ the continuum point is the Gaussian UV fixed point while the IR is strongly coupled; no single expansion reaches both ends.

(ii) **The answer is nonperturbative.** At every finite order in g₀² the invariant-mass spectrum of the tr F² two-point function starts at zero (two massless gluons). Any bound that is analytic in g₀² and vanishes at g₀ = 0 is compatible with perturbation theory and therefore cannot be the gap. A proof must produce m ≍ Λ ∼ a⁻¹ e^{−8π²/(β₀g₀²)}, non-analytic at g₀ = 0. This rules out, as a matter of principle, every "curvature/convexity in the coupling" scheme.

(iii) **Gribov.** The Hessian of S is indefinite at large field (Savvidy-type negative modes), the Faddeev–Popov operator on a gauge slice can reach zero (the Gribov horizon), and the orbit space has singularities. Hence no uniform log-Sobolev/Poincaré inequality and no convexity argument can yield a β-uniform gap. No known theorem even excludes a β-transition.

## 3. What I can prove

**A (strong coupling).** For fixed N there is β₀ > 0 such that for β ≤ β₀, connected correlators of 𝒜 decay exponentially, uniformly in L, with rate m(β) > 0. Sketch: expand e^{β Re tr U} = Σ_λ c_λ(β)χ_λ(U), |c_λ| ≤ (Cβ)^{|λ|}; the resulting surface gas satisfies the Kotecký–Preiss criterion (Σ_{γ∋x}|w(γ)|e^{α|γ|} ≤ α) for β < β₀, giving clustering. This is genuine, but β₀ is a fixed small number — not the continuum.

**B (perturbative no-go).** By (ii), every finite-order spectral density is gapless; this is a rigorous statement about what a proof cannot be.

**C (equivalence).** Reflection positivity plus the spectral theorem identifies Δ with the clustering rate; once uniform clustering holds, L → ∞ and finite-volume corrections are exponentially small and harmless (localization). The a-uniform rate is the only hard part.

## 4. Attack: RG transport and the exact missing lemma

Block-spin gauge-invariantly with scale factor 2: effective action A_b, coupling g_b, block spacing a_b = 2^b a. The physical mass is RG-invariant; let m_b := Δ(a_b)/a_b be its dimensionless version at scale b. At b = B with a_B ≈ 1/Λ the running coupling is O(1); *if* the block theory lies in the strong-coupling class of A, then m_B ≥ c > 0 in block units, hence Δ ≥ c/a_B ≥ c′Λ > 0, uniformly in a, and the proof is complete.

The minimal statement that settles everything:

**(RG-transport).** There exist a gauge-invariant block map and a class 𝒞 of effective actions containing the strong-coupling Wilson actions such that (a) 𝒞 is invariant up to errors η_b with Σ_{b≤B(a)} η_b ≤ ½c for all a; (b) the map carries the small-g₀ Wilson action into 𝒞 at scale a_B; (c) the mass functional changes by at most η_b per step.

I cannot prove this. Parts (a)+(b) are scale-uniform large-field control of A_b — Gribov copies, topological sectors, locality, and control of the marginal operator tr F² — over B(a) ∼ 1/g₀² steps. Part (c) requires the character-expansion radius to dominate the O(1) coupling at the confinement scale; this is not known. For SU(2) all character coefficients are positive, so Griffiths-type domination is conceivable; for SU(N≥3) they alternate in sign, which blocks the standard ferromagnetic comparison. No exact 4D duality (center-vortex or random-surface) is available to replace (c).

## 5. Honest conclusion

(a) **Proved:** strong-coupling mass gap, uniform in volume; reflection-positivity/transfer-matrix equivalence; perturbative gaplessness; Gribov obstruction to convexity/curvature strategies.

(b) **Conditional:** (RG-transport) ⇒ continuum mass gap.

(c) **Open (Clay).** I am stuck exactly at scale-uniform effective-action control from the Gaussian UV fixed point down to coupling O(1) at scale 1/Λ, with a summable mass error — equivalently, at proving the strong-coupling radius dominates the dynamically generated coupling. No sixth-species shortcut is visible to me: the mass gap *is* the statement that this RG transport is stable, and that stability is precisely what is missing.

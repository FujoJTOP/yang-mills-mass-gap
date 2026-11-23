# G2H — Independent (Blind) Attack on the 4D Yang–Mills Mass Gap

**Granted, nothing else.** Pure SU(N) YM in 4D exists as a Wightman QFT: Hilbert space ℋ, unique Poincaré-invariant vacuum Ω, a self-adjoint Hamiltonian H ≥ 0 with HΩ = 0, and Källén–Lehmann (KL) spectral measures; the theory is the scaling limit of the Wilson lattice theory. I assume no prior reduction.

## 1. What "mass gap" means, precisely

Let ℋ_phys be the closure of the gauge-invariant (singlet) polynomial algebra applied to Ω (physical/neutral sector), and H_phys := H restricted to ℋ_phys ⊖ ℂΩ. Define Δ := inf spec(H_phys). **Mass gap** means Δ > 0.

Equivalently, for the translation spectral measure ν of the neutral scalar sector — or, for any gauge-invariant scalar O with ⟨Ω|O|Ω⟩ = 0, its KL measure ν_O — we need **inf supp ν ≥ Δ² > 0** (convention s = m²), i.e. no spectral weight in (0, Δ). Uniformity in the limit: there is a single Δ > 0 such that for every lattice spacing a, **after** the thermodynamic limit L→∞, the measure ν_a has no support in (0, Δ), with Δ independent of a as a→0. A finite volume L has a spurious O(L⁻²) gap from level discreteness; the content is uniformity in volume and spacing.

## 2. WARNING — the susceptibility trap (do not fall in)

The tempting move is to prove a **susceptibility bound** χ = ∫ dx ⟨O(x)O(0)⟩_c = ∫ (1/s) dν(s) < ∞ and conclude a gap. **This is false.** Bounded χ controls only a first weighted moment; a small-weight near-zero continuum keeps χ finite while inf supp ν = 0.

Explicit admissible counterexample. Take dν(s) = s^{1/2} 1_{[0,1]}(s) ds (a plausible CFT-like spectral density; for a 4D CFT primary of dimension Δ_O, ρ(s) ∼ s^{Δ_O−2}, any Δ_O > 2 gives the same structure). Then inf supp ν = 0 yet χ = ∫₀¹ s^{−1/2} ds = 2 < ∞. Indeed ⟨O(x)O(0)⟩ ∼ |x|^{−2Δ_O} decays only polynomially, χ converges if 2Δ_O > 4, and there is **no gap**. So no finite set of moment/susceptibility/energy-moment bounds can yield "inf supp > 0": the support of a positive measure is not determined by moments (indeterminate moment problem), and the thin low-mass tail is exactly invisible to low moments. Any argument that uses only χ and not YM-specific dynamics is invalid.

## 3. What is genuinely provable (PROVED)

- **(P1) KL equivalence.** Gap ⇔ exponential clustering of connected gauge-invariant correlators. For spacelike x: ν([0,ε]) ≥ e^{√ε|x|}⟨O(x)O(0)⟩_c. Hence decay ⟨O(x)O(0)⟩_c ≤ C e^{−γ|x|} forces ν([0,ε]) = 0 for every ε < γ²; conversely inf supp ν = Δ² gives decay ≤ (total mass)·e^{−Δ|x|}. Thus a uniform decay rate γ for all spacelike x, uniform in a and L, is *equivalent* to a gap ≥ γ². (Standard KL; I can prove the implication cleanly.)
- **(P2) Moment/susceptibility bounds do NOT imply a gap.** Section 2; rigorous counterexample.
- **(P3) Strong-coupling lattice gap.** For β < β₀ the Wilson lattice theory has m_lat(β) ≥ c > 0 in infinite volume (polymer/cluster expansion, Osterwalder–Seiler). Known theorem; only covers small β.
- **(P4) Continuum map.** By asymptotic freedom the critical curve is β = β(a) ∼ (11N/48π²) log(1/(aΛ)); the physical mass is m_phys = m_lat(β(a))/a. Standard.

## 4. The conditional skeleton (PROVED implication)

**(C) Uniform physical lower bound ⇒ gap.** Assume there are Δ > 0, C < ∞ with, for all a and (after L→∞),
⟨O(x)O(0)⟩_c^{(a)} ≤ C e^{−Δ|x|}, x spacelike.
Since the scaling limit gives pointwise W_a(x) → W(x), the bound passes to the limit: |W(x)| ≤ C e^{−Δ|x|}. By (P1) the limiting ν has no weight in (0, Δ²). **This step is airtight.** The whole problem is the hypothesis — *uniformity of the decay rate* — not the limit passage.

## 5. Where I get stuck — the minimal missing statement

The granted Wightman existence is **insufficient by itself**: the free massless field satisfies every Wightman axiom with inf supp ν = 0. So the proof must use the specific construction. The one statement that settles everything:

> **(MS)** On the critical curve β = β(a), the lattice mass gap in physical units is bounded below uniformly: there is Δ > 0 with m_lat(β(a))/a ≥ Δ for all a ∈ (0, a₀), uniformly in volume.

Equivalently: a volume- and coupling-uniform lower bound on the inverse correlation length of *every* gauge-invariant connected correlator, persisting from β<β₀ to β→∞.

I cannot prove (MS) from first principles. Why each route fails:
- **Monotonicity in β?** Unknown; m_lat(β) is not proven monotone, and even monotone positivity would need a uniform floor along the critical line.
- **Weak-coupling cluster expansion?** The gauge-fixed quadratic form is indefinite (Savvidy-type negative modes), the perturbative series is asymptotic with renormalon/IR ambiguities, so it cannot certify a *nonperturbative* mass.
- **Positivity / Perron-type bounds on the transfer matrix T = e^{−H}?** Reflection positivity gives 0 ≤ T ≤ 1, TΩ = Ω, but T is non-compact; no uniform spectral-gap bound at strong magnetic coupling follows.
- **Moment / sum-rule methods?** Blocked by (P2): they see moments, not support.

## 6. Triage

- **Proved:** (P1) gap ⇔ uniform exponential clustering; (P2) susceptibility/moment bounds are insufficient; (C) uniform physical mass ⇒ gap; (P3) strong-coupling gap; (P4) AF continuum map.
- **Conditional:** [[the gap holds]] **iff** (MS) holds; gap in the limit is exactly equivalent to uniform-in-(a,L) exponential decay (P1 + C).
- **Open (the missing machine):** (MS) — the nonperturbative, cutoff- and volume-uniform transport of the strong-coupling gap across the whole coupling range to the scaling limit. This is the precise gap in the argument; everything else is bookkeeping.

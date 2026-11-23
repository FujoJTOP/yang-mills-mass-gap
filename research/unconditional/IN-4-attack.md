# IN-4 attack — CPB-S strong-side signed gas

## 1. Statement (IN-4 / Proposition B)
IN-4 (CPB-S): ∃ q_c, δ_S>0 with q(μ)≥q₁ ⟹ cert μ ≥ δ_S, plus a σ-version. Its sharp kernel, **Proposition B**: if defects are closed 2-surfaces (dF=0), tension σ(q)≥c₀q², entropy N(n)≤e^{γn}, then z(q)≤C e^{−cq²}, with q_c≈0.835/√c₀ (the z=½ crossing, ≈√(ln2/c₀)).

## 2. True status
Proposition B is **proven as a conditional lemma** (geometric series; `handoff-proof` §4.4). The unproven content is H1–H2, and neither follows:
- dF=0 is a *flux* cocycle, not *support*-closure (BL-1). BL-3: the x₄=0 slice is a closed 2-chain of codim-1 support with q-independent per-block cost ⟹ H1 ⇏ H2.
- q² is **neutral at q≈1** (BL-5, q²=q). At the crossover Prop B gives zero gain; only the absolute σ(1) or κ_c enters.

2D O(N) is a theorem; the 4D port (no hard constraint, no additive index) is unproven. No rigorous 4D signed-gas result exists.

## 3. Advance attempt
Not provable with current mathematics; three independent obstructions:
1. **Sector mismatch** (NP-7): the only geometric candidate, the center vortex, carries charge in the *charged* sector; the mass gap lives in the *neutral* sector. BL-6 default κ_c≈1.
2. **Neutrality at q≈1** (BL-5) makes Prop B the wrong object at the crossover — the requirement is σ(1), a dynamical weight, not a scaling exponent.
3. **c₀ is dynamical** (BL-3). Physically c₀=T₂/m² (membrane tension over m²). If T₂~m², c₀~O(1): c₀=1 ⟹ q_c=0.83, λ*≈6.8. If T₂~Λ², c₀=(Λ/m)²≈0.02–0.06 (m/Λ=4–8) ⟹ q_c≈3.4–5.9, λ*≈14–30 — the wall worsens by up to an order of magnitude.

Real advance: c₀ is a *physical ratio*, in principle measurable (multi-level IMC, center projection), whereas κ_c has no direct route. The σ(1) lattice proxy failed (Gaussian-tail artifact: it measures field width, not action). No established value exists.

## 4. Flip condition
IN-4 flips only if one of:
(a) prove T₂≥c₀m² with c₀≳2.6 in the neutral sector (confinement-scale dynamics — essentially the target);
(b) construct a block-invariant, θ=0, gauge-visible additive Z₂ charge in the *neutral* sector (NP-7; NP-T7b no-go pending) — no candidate exists;
(c) exhibit κ_c≤0.20 on fixed boundaries;
(d) exclude codim-1 closed 2-chains (BL-3) as dominant defects — dF=0 alone cannot.

Absent these, IN-4 remains open and Proposition B is vacuous exactly where it is needed.

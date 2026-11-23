# W4 — The Polchinski trapping-region target

**Verdict: clear, non-excluded. The log log(1/a) is not a hidden wall: it is the marginal coupling, resummed as in the solved marginal templates. The risk is technical — a gauge-covariant, all-orders, large-field uniform bound at the one dimension current methods stop.**

## 1. The lemma
Gauge-covariant P(𝒜/𝒢) (CCHS), Polchinski flow ∂_j L_j = ½(δL/δA)·Ċ·(δL/δA) − δ²L/δA²:Ċ, j=0..J, e^{−J}=a. Write L_j = g_j^{−2}∫Tr F² + R_j (θ=Tr F F̃: total derivative, drops); ‖·‖_X = BKAR/Kotecký–Preiss norm on gauge-invariant dim>4 functionals.

**Trapping lemma.** ∃ g*,r,c>0, bare g₀(g*) with, ∀J:
(T1) 0<g_j≤g*, ‖R_j‖_X ≤ r j^{−θ} ∀j≤J;
(T2) ‖R_{j+1}‖_X ≤ e^{−c}‖R_j‖_X + C g_j⁴, uniform in cutoff and volume;
(T3) g_{j+1}=g_j−½β₀g_j³+O(g_j⁵), remainder uniform.
Then sup_j‖R_j‖_X<∞; j→∞ limit ⇒ continuum measure (G1). Since Σg_j⁴<∞ although Σg_j²~(1/β₀)log J diverges, the log is absorbed into g_j²≈1/(2β₀j).

## 2. Footholds
**(i) Templates.** GN₂ (FMRS, CMP 103 (1986) 67): AF makes the UV a *convergent* OS construction — but fermionic, massive by input, not gauge. KT (Fröhlich–Spencer, CMP 81 (1981) 527): trapping box around the marginal fixed line ⇒ exponential decay; structurally exact but Abelian/Gaussian. d=4 marginal Polchinski (Bauerschmidt–Brydges–Slade, LNM 2242, 2019) *proves* a trapping region with log resummation at exact marginality in d=4 — scalar, weak-coupling. Gauge covariance and non-Gaussian marginals missing.

**(ii) Beta coefficients.** b₀ (Gross–Wilczek/Politzer 1973) and 4-loop (van Ritbergen–Vermaseren–Larin, PLB 400 (1997) 379) are computed, not proved nonperturbatively; none rigorous for 4D gauge. Coefficients fix asymptotics, not uniformity.

**(iii) Gauge partials.** Balaban (CMP 109 (1987) 249): UV stability, gauge-fixed, non-uniform. Müller–Schiemann (CMP 110 (1987) 261): rigorous 4D SU(2) continuum limit — but hierarchical toy. CCHS (2201.03487), Chevyrev–Shen (2503.03060): gauge-covariant but d=3 (super-renormalizable; Hairer's degree condition fails at d=4).

## 3. First sub-lemma (provable now)
**L0 (free contraction).** The linearized free Polchinski semigroup on gauge-invariant dim>4 functionals contracts, ‖P_{>4}e^{−tL₀}‖≤Ce^{−ct}, uniformly in cutoff and volume — the half needing no β, only lattice-uniform Gaussian estimates and gauge-invariant projection. **L1** (next): the one-loop box (T3) with uniform O(g⁵).

## 4. Honest assessment
Hard but clear: the log log is resummable, not a wall. The obstruction is one place — a gauge-covariant, all-orders, Savvidy-safe large-field bound at d=4, where regularity structures (degree condition), Balaban (uniformity), BBS (non-Gaussian marginal coupling) stop. The lemma lives in the UV weak-coupling regime where AF sits, so G1 is plausible; G2 still needs the strong-coupling IR. L0 is easy; physics is in (T1)+(T3).

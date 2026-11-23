# W6 — Bakry–Émery at β = 2.5: attempt and exact failure

**Theorem attempted.** Gauge quotient M = G^E/G_site, bi-invariant metric, μ_Λ = Z⁻¹e^{−S_W}dU. Claim CD(m,∞) uniform in volume:
`Ric_M(v,v) + Hess S_W(v,v) ≥ m|v|², ∀v∈TM, m(β)>0,` (★)
⇒ gap ≥ m, LSI constant ≤ 2/m [Otto 2001; Sturm 2006; Lott–Villani 2009; von Renesse–Sturm 2005] ⇒ volume-uniform mass gap [SZZ, arXiv:2204.12737, CMP 400 (2023) 805]. S_W gauge-invariant ⇒ orbit Hessian 0, (★) descends [Sturm 2006]. Is m > 0 at β = 2.5?

**Verdict: FALSE — pointwise curvature negative; not a proof.**

**1.** Ric of the group factor is a β-independent positive constant (SU(2): ½g; SU(3): ¾g). Hess S_W = β·Hess(−Re Tr/N) ∝ β, bounded but reaching ≈ −1.1β. Numeric (`w6_check.py`): SU(2) L=4, β=2.5, full fine Hessian min = −2.81 (−1.12β), 21.7% of 3072 directions negative; SU(3) L=3 min = −2.77 (−1.11β), 42% negative. True Ric ⇒ min(Ric+Hess S) ≈ −2.3 / −2.0; even generous Ric = +1.5 gives −1.3. So m < 0 for every m > 0. (O'Neill +¾|A|² ≥ 0 only helps; still fails.)

**2.** On the *same* configs the plane-wave collective subspace (the 11-slc-blocked object) is positive — SU(2) min +0.40β, 0% negative, reproducing "H<0 fraction 0.000" — but the most-negative eigenvector has only 3.4% (SU(2) L=4) / 9.7% (SU(3)) weight there: negatives are 90–97% orthogonal. The probe tests a low-momentum slice; BE quantifies over all of TM. At SU(3) β=2.5 the broader L=3 plane-wave subspace is itself 27% negative — positivity is construction/volume-dependent.

**3.** Negatives are high-momentum (~90–95% in the fluctuation block, MARGINAL.md), yet "UV" is irrelevant pointwise: blocking acts on the measure, not microscopic curvature, so every tangent direction enters Γ₂. At k→0, where blocked/induced mass would enter, the restoring p²→0 and the marginal proxy G″ ≤ E[PᵀHP] soft value ≈ +0.043 − 0.208 = −0.165 < 0 (VAR-FINAL.md); the induced mass would be the gap itself — circular.

**4.** Hess S ∝ β with Ric = O(1) ⇒ BE only for β ≲ Ric/1.1 ≈ 0.5–1.4, SZZ's strong-coupling window (rigorous β < 1/48 in d=4). β = 2.5 is far outside; blocking does not move the pointwise curvature.

**Exact break:** BE curvature in the Savvidy-negative directions [Savvidy 1977; Nielsen–Olesen 1978; Leutwyler 1981], orthogonal (≤10%) to the collective plane-wave subspace, Hess S ≤ −1.1β exceeding any Ric; and the k→0 mode with G″ ≈ −0.165 < 0. A volume-uniform BE/LSI at β = 2.5 is **not** obtainable this way.

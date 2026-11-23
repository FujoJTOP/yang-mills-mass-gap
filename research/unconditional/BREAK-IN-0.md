# BREAK-IN-0 — Y3: subsequence-independence of g_W(r)

**Setup.** A volume subsequence yields a limit state μ on the gauge-invariant algebra 𝒜_gi, with ergodic (Choquet) decomposition μ = ∫_X ω ν(dω), ν a probability on pure phases ω (each clustering, OS-reconstructible; §6.1/Thm 6.1). Let c_r(μ) be the normalized curvature-squared coefficient in the scale-r blocked effective action; `g_W(r) = c_r(μ)^{−1/2}` (Def. 9.17).

**Affine data survive (the proven model).** The trace amplitude is a *linear* functional, D = μ(φ) = ∫_X ω(φ) ν(dω). By Thm 7.1 the constant `3 dim_R g/4` is the same for *every* phase ω, so ω(φ) = d_* for all ω and D = d_* for every ν. Linearity ⇒ phase-mixing invariance; this is why the leading trace law (17) is subsequence-independent.

**Why nonlinearity breaks it.** Blocking is deterministic, so the coarse measure is the ν-barycenter: `dμ_coarse/dHaar = ∫_X e^{−S_ω} ν(dω)`, with S_ω the scale-r action of phase ω. Hence
`S_eff(μ) = −log ∫_X e^{−S_ω} ν(dω)`,
and c_r(μ) is the quadratic coefficient of this log-sum-exp — *not* linear in ν:
`c_r(∫ μ_ω ν(dω)) ≠ ∫ c_r(μ_ω) ν(dω)`.
Schematically `c_r(mix) = c_r(ω) + δ(ν)`, where δ(ν) is a variance/entropy correction, nonzero whenever phases differ in their *subleading* blocked action. Thus even if every phase shares the same quadratic coefficient c_*, the mixture — and different subsequences with different ν — carry a ν-dependent shift. **Exact first failure: coefficient extraction does not commute with the phase average; the variance term δ(ν) is the obstruction.** Thm 7.2 pins only the leading c_*, not δ(ν).

**Does clustering / uniqueness close?** Clustering: **no** — E4 holds per component, not for the mixed μ (§6.2); a clustering limit pins each subsequence only to *some* pure phase, so different ones may pick different phases (two-phase Ising). Uniqueness: **yes**, trivially.

**Clean lemma (closes).** If the net {μ_Λ} has a *unique* limit point μ, then g_W(r) is subsequence-independent. *Proof:* `g_W(r) = c_r(μ)` is a function of μ. ∎

**Verdict / honest correction.** Y3 closes only under (i) UNI-lim [open = IN-6], or (ii) *scale-r phase-coincidence*: all phases across all limit decompositions share the same scale-r blocked action, forcing δ(ν)=0. (ii) is strictly weaker than UNI-lim, stronger than the proven leading coefficient, and not implied by clustering. Also, Y3 (even ∀r) is a single running function and does **not** imply UNI-lim: the correct logic is `UNI-lim ⇒ Y3`, not "Y3 = unique limit state." Absent (i)/(ii), the failure is precisely `δ(ν) ≠ 0`.

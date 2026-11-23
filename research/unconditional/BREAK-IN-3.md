# BREAK-IN-3 — T1: κ_c in the SU(2) surface representation

> 2026-11-23 · scripts `research/unconditional/t1_kappa.py`, `t1_beta6.py` · numpy/scipy, reuses `lattice-cpb/su3._sample_w0`.

**Question.** Is the sole IN-3 escape real — `κ_c = (z_signed/z_abs)^{1/|C₀|} ≤ 0.80` at q≈1 (CPB-W)?

**What is computable.** Only the *microscopic* SU(2) Wilson surface/character weights can be evaluated on a small lattice. I ran, exactly (no MC error), at β=2.0–2.6 (intermediate-coupling crossover) and β=6.0–6.3 (the SU(2) point with ξ≈1, i.e. q≈1 per `lattice-cpb/STATUS.md`):

- **Character coefficients** `c_r(β)=(2(2r+1)/β)I_{2r+1}(β)`. Closed form verified to 1e-15. **c_r>0 for all r, all β** — e.g. β=2.2: (1.740, 0.812, 0.074, 0.003…); β=6.0: (20.45, 30.15, 13.28, 2.93, 0.39, 0.034, 0.002).
- **Link multiplicities** `M({r_i})=dim[⊗r_i]₀`: non-negative integers (verified 1,1,2,3,5; M(½,½,½,½)=2).
- **Exact signed-vs-absolute enumeration** of a 12-plaquette 4D defect 2-chain carrying six 4-valent links (so M=2 enters), all 3¹² assignments: **0 negative terms**, `z_signed/z_abs = 1.000000000000`, hence **κ_c = 1.000000000000** at every β (2.0, 2.2, 2.4, 2.6, 6.0, 6.3).
- **Checkerboard SU(2) heat bath** (validated vs Metropolis and an exact 1D integral; ⟨P⟩ matches known SU(2)), L=6: ⟨P⟩ = 0.501(1), 0.570(1), 0.632(1), 0.671(1) at β=2.0/2.2/2.4/2.6; 0.870(1) at β=6.0.

**Result (i).** In the surface/character representation **κ_c = 1 exactly**, at the crossover and at q≈1. The surface-representation T1 escape (κ_c≤0.80) is **closed**.

**Result (ii) — analytic bound κ_c≤0.8: not established, and false here.** Since c_r>0 and every multiplicity ≥0, every surface weight w(C)≥0, so z_signed=z_abs and κ_c=1. The general inequality is κ_c≤1 (|Σw|≤Σ|w|); SU(2) saturates it. The IN-3§2 "high-rep alternation" crack does **not** exist for SU(2): c_r>0 for *all* r (at β=6 even c₁>c₀). No argument reaches κ_c<1.

**Honest boundary.** The CPB object is `w_q(C)`, the *cluster* weight of the RG-blocked effective action at block spacing a_k (q=m a_k≈1). That is not a bare surface weight. Blocking a positive measure gives a positive Boltzmann factor, but the Ursell/connected expansion of log carries alternating signs, so κ_c<1 is not excluded by positivity there — and a small-lattice MC cannot reach w_q (`STATUS.md`: "κ_c has no direct numerical route"). So T1 on the *blocked* weights remains unmeasured; the *surface-representation* T1 measured here returns κ_c=1, confirming BL-6's default. Net: one escape route removed; the other still unmeasured, with no known mechanism. This does not prove CPB-W false.

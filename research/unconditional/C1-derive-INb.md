# C1 — Deriving IN-b: five angles and the cross-channel obstruction

> 2026-11-23 · Dedicated pass over the five prescribed routes. **Verdict: none is a theorem; IN-b is irreducibly physical.** The five collapse to one object. §Verdict isolates a sharp reason that is new here: the ordering relates two *orthogonal* symmetry sectors, and at θ=0 every cross-channel positivity inequality is identically trivial.

**Setup.** ν_c = Källén–Lehmann measure of a gauge-invariant interpolator in channel c; supp ν_c ⊂ [m_c²,∞); f₁² = ν_c({m_c²}). IN-b needs m_{0⁻⁺} ≥ m_{0⁺⁺} (weak form: every channel outside 0⁺⁺/2⁺⁺ above δ).

**(a) Lattice–continuum ratio theorem.** None exists. Measured ratios (Morningstar–Peardon) are finite-`a` Monte Carlo plus an uncontrolled `a→0` extrapolation. A proven continuum ratio presupposes uniform convergence of both edges — a gap statement per channel; with any positive 0⁺⁺ certificate it already yields the 0⁻⁺ gap, so it is the target, not a shortcut. Rigorous control exists only at β<β₀ (Osterwalder–Seiler), where the gap is already a theorem, and never reaches scaling.

**(b) Large-N with controlled 1/N.** No. Both channels are single-trace, m∼N⁰; the 4D master field is unconstructed; the 1/N series is asymptotic with no error bound. An exact N=∞ ordering does not transfer to N=3 without control. The one exact 0⁻⁺ relation, Witten–Veneziano, needs dynamical quarks and makes the η′ *light* (O(1/N)); pure glue is the singular N_f→0 limit of *full* QCD, not quenched, so WV does not continue.

**(c) Bootstrap bound on m_{0⁻⁺}.** No. The defect bootstrap sees 0⁺⁺ (tr F², trace anomaly) and 2⁺⁺ (T_{μν}, normalization fixed by its role as translation generator) because those are protected. q = tr F F̃ has no protected source; and q×q is P-even at leading order (two P-odd factors), so 0⁻⁺ is no s-channel pole of a local four-point bootstrap. The S-matrix bootstrap takes external masses as input and bounds couplings.

**(d) New symmetry / non-invertible structure on local operators.** No. At θ=0 every exact symmetry of the local gauge-invariant algebra is Poincaré×P×C; all are linearly realized and non-anomalous in pure YM, so they label J^{PC} but cannot order it. 1-form Z_N, 2-groups and non-invertible defects act on lines and surfaces; local glueballs are colour singlets, N-ality 0. A 4D non-invertible 0-form duality is known only with extended SUSY. At θ≠0, P/CP mixing couples the sectors, but the mixing — hence all ordering information — vanishes as θ→0.

**(e) Weingarten / QCD-inequality lower bound.** No, and here the obstruction is sharpest. Weingarten inequalities need a positive *same-channel* comparison and bound the light state from **above**. A cross-channel comparison needs a nonzero mixed correlator; at θ=0
$$\langle q(x)\,\mathrm{tr}\,F^2(0)\rangle \equiv 0 \qquad (P:\ q\ \text{odd},\ \mathrm{tr}\,F^2\ \text{even}),$$
so Cauchy–Schwarz degenerates to `0 ≤ C_q C_S`. A nonzero cross-channel kernel needs θ≠0, excluded. A lower bound on one correlator's decay rate is verbatim confinement.

**Verdict — why irreducibly physical.**
(i) *Cross-channel triviality.* `m_a ≥ m_b` relates two superselection sectors. Every rigorous ordering tool is a positivity inequality needing a positive cross-channel kernel; at θ=0, P and C make every two-point mixed kernel vanish identically, and a mixed-parity superposition has no definite edge. No two-point comparison can order the channels.
(ii) *Single-sector residue.* The only remaining route is a lower bound on the 0⁻⁺ edge residue `f₁²>0`. The lowest-dimension 0⁻⁺ interpolator q=∂·K is a total divergence; its residue is unprotected (subtraction ambiguous), and the free-theory no-go kills every g-continuous protected functional (SCHUR-CONVEX). No anomaly or Ward identity supplies `f₁²>0`.
Thus a cross-channel kernel, a channel-blind edge certificate, and `f₁²>0` are one object, and each is gap-strength. The same obstruction recurs for every remaining channel in the weak form (odd C, spin ≥ 3).

The only regime where an ordering is even controllable is β<β₀ strong coupling, which never reaches the continuum. **No theorem; IN-b stays physical.** *No fabricated proof.*

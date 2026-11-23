# IN-3 attack — κ_c at the crossover

## 1. Restatement
IN-3 (CPB-W, weak side): ∃`λ_w,δ_W>0` with `λ_G(μ)≤λ_w ⟹ cert μ≥δ_W`. Reach needed: `λ_w≳7.3` (`q₀=1` at `λ≈7.28`; `β≈2.47`, SU(3)) — the intermediate-coupling crossover. The Gaussian-weighted class route reduces to a local 3rd-cumulant bound with constant `C_cum`; passage needs `C_cum≲1.6`, while standard Kotecký–Preiss at `g²≈2` gives `≥2`. BL-5 (theorem): additive `q²` is neutral at `q=1`. The only mechanism that lowers `C_cum` by an O(1) factor is **signed cancellation** in the fixed-boundary surface sum
`|Σ_{∂C=∂C₀} w_q(C)| ≤ e^{−c₀q|C₀|} κ_c^{|C₀|}`, `κ_c<1`.
`κ_c` is a per-block **weight** ratio `(z_signed/z_abs)^{1/|C₀|}`; CPB-W needs `κ_c≤0.80` (`δ≥ln1.25`), CPB-S `≲0.20`. One `κ_c` serves both.

## 2. BL-6 — defended
BL-6 ("no known route; default `κ_c≈1`") is **robust**. I widen its evidence into an exhaustive search for a non-trivial local sign cocycle: strong-coupling positivity `c_r(β)>0` ⇒ `κ_c=1` exactly; instanton/Chern–Pontryagin θ-phase inert and non-per-block at `θ=0`; CS surface phase needs the θ-term (absent); 1-form center `Z_N^(1)` charges **loops**, not closed surfaces, fixed boundary ⇒ fixed charge, no relative sign, and no 2-form symmetry exists; global/Witten anomaly is fermionic; self-intersection/genus is embedding data, not additive per block; Griffiths positivity forces *same*-sign ⇒ `κ_c=1` (the opposite). All fail.
**One crack (steelman):** Evidence 1's high-rep alternation for SU(N>2) is unproven subdominant at `β≈2.47`; an O(1) rep-correlated flip could lower `κ_c`. Untested, but those coefficients are exponentially suppressed vs. leading reps — implausible, not excluded.

## 3. Advancing κ_c — all blocked
- **A. Analytic at q≈1:** needs the blocked weights `w_q(C)`, unknown at the crossover (that *is* the problem). Blocked.
- **B. New additive index:** none survives (§2). Confirms BL-6.
- **C. Other expansions** (Wasserstein–CD, thimble/CL, resurgence): bypass `κ_c` but die elsewhere (Savvidy concavity; real-β has nothing to complexify; 4D renormalons). No `κ_c<1`.
- **D. Empirics** (T1, `z_signed/z_abs` at `q=1`): not yet measured. Hessian proxies (W4 positive blocked band; 11-slc `G''_soft≈−0.165`, DEAD-lean) are *consistent* with `κ_c≈1` but are not the ratio.

## 4. Flip condition
Any one flips IN-3: (i) lattice T1 measures `κ_c≤0.80` at `q=1`; (ii) a proven additive surface cocycle (a new `Z₂` index with non-trivial holonomy on surface-space); (iii) escape without `κ_c` via larger single-block weight `σ(1)≥0.5+γ_ent+ln3` (T2). Absent these, the no-route verdict stands and IN-3 remains the physical wall.

*Status: §2 exhaustive-search is an argument, not a theorem; BL-6 is a negative result about routes, not a proof CPB-W is false.*

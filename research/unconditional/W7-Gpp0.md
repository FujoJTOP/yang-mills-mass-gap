# W7 — Exact sign of the marginal curvature G″(0)

**Verdict.** The physical (Haar-measure) marginal curvature G″(0) is **positive in every mode tested**: SU(2) β=2.5 +0.27 (k=π/2) and +0.32 (k=π/3); SU(3) +0.49 (β=5.7) and +0.86 (β=5.7 soft blocked); SU(3) β=2.5 soft blocked +0.51. The *naive flat* energy-variance object E[H]−Var(∇S) is ≈0 (range −0.03…+0.04): the force variance saturates the curvature. So SLC / a k→0 mass-gap direction **survives** (G″>0 physically); the earlier cond_var −0.165 was a ~5× over-subtraction of a true ≈0 margin, not a real negative.

**Method (exact, no Gaussian proxy).** Fix one collective mode v — k=2π/L transverse plane waves, and the b=2 block-uniform softest eigenvector of Ĥ=E[PᵀHP] — and let a=vᵀ log U. The equilibrium marginal P(a) gives V_eff(a)=−log P(a) and V_eff″(0)=E[vᵀHv|a=0]−Var(vᵀ∇S|a=0)=G″(0) exactly (constrained free-energy identity). H and vᵀ∇S come from the Wilson action along U→e^{isv}U (finite differences, log-free). κ from an even-polynomial fit of −log P, bootstrap SE (Gaussian control recovers 1/σ² to 5%). SU(2): L=4 N=4000, L=6 N=800 (48 modes); SU(3): L=4 N=1500 (36 modes).

| case | modes | V_eff″(0)=κ (β-units) | 1/Var(a) | E[H] | Var(∇S) | flat E[H]−Var |
|---|---|---|---|---|---|---|
| SU(2) β=2.5 L=4 (k=π/2) | 48/48>0, +0.22..+0.34 | **+0.274±0.004** (+0.11) | 0.27 | 2.46 | 2.44 | +0.02 |
| SU(2) β=2.5 L=6 (k=π/3) | 48/48>0, +0.20..+0.52 | **+0.318±0.009** (+0.13) | 0.27 | 2.45 | 2.46 | −0.01 |
| SU(3) β=5.7 L=4 (k=π/2) | 36/36>0, +0.22..+0.63 | **+0.485±0.012** (+0.085) | 0.45 | 3.18 | 3.19 | +0.04 |

Part B — the cond_var object. SU(3) β=2.5 (cond_var's deep-SC regime; Ĥ min=+0.137) block-uniform softest blocked mode: exact G″=**+0.51±0.08** (E[H]=0.22, Var=0.24, flat E−Var=−0.03, branch-trimmed +0.56). [SU(2) β=2.5 same object +0.270±0.034; SU(3) β=5.7 +0.86 (N=300).] cond_var's −0.165 came from a 3-config Ĥ eigen-direction and a mixed-coordinate regression (absolute-log A′ vs left-variation ∇S, R²≈0) that subtracted the *full* force variance → over-subtraction.

**Reading.** Two objects. (i) Physical V_eff″>0 — the transverse inverse propagator; part of it (≈+0.25) is the compact-group Haar volume-element curvature, geometric not action. (ii) Flat E[H]−Var≈0: ∇S is nearly independent of a, so conditioning barely lowers Var(∇S) and the action curvature is saturated by force variance (bare E[H]≈99% cancelled). SLC holds physically, but the *dynamical* (action) margin is ≈0 — near-critical, not robust.

**Limits.** L=4/6 reach only k=π/2, π/3 — NOT k→0; L=4≠L=6 volume, so not a clean k-scan. A positive low-k curvature is consistent with, not a proof of, a k→0 mass gap. SU(3) principal-log branch cuts (every β=2.5/5.7 config has a near-±π eigenvalue) give ~10–20% systematic. Scripts: `research/unconditional/w7/`.

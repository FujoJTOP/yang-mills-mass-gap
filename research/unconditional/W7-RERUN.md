# W7-RERUN — G″(0) with valid configs and the correct +W definition

**Verdict.** Re-run on valid configs with the full identity `G″(0)=E[PᵀHP|A′]−Var(Pᵀ∇S|A′)+W`:
**G″(0) > 0 in all four cases and volume-robust.** The positive curvature is the **Haar-Jacobian
term W**, *not* `E[H]−Var`: the action bracket is ≈0 everywhere while the Jacobian curvature
`E[Φ″]` is the large, stable positive piece. W7's *number* (κ≈0.27) was right; its *claimed
identity* `V_eff″=E[H]−Var` was wrong and omitted +W (VERIFY correct).

**Method (flat exponential chart).** With `X=Σxₐtᵃ`, `a=vᵀlogU` (v a transverse plane mode),
displace flat `X→X+sV`. Then exactly
`V_eff″(0)=E[H_vv+Φ″|a=0]−Var(G_v+Φ′|a=0)`, `Φ=−logJ`, `H_vv`,`G_v` from `S_W`, `Φ′`,`Φ″`
analytic in the eigenangles `φᵢ` (`Φ=Σ_{i<j}g(φᵢ−φⱼ)`, `g(δ)=−2log[sin(δ/2)/(δ/2)]`).
`G″` is taken as the KDE curvature of −logP(a) (cross-checked by `1/Var`); `W≡E[Φ″]−Var(Φ′|0)−2Cov`.
Gaussian-kernel conditioning at a=0 is essential (`Var(Φ′)` is dominated by its `a`-dependence,
corr(a,Φ′)=0.91; unconditioned `E[H]−Var` is meaningless here).

**Branch fix (SU3).** W7's "0/1500 branch-safe" used `max|eigangle|<π−0.15` — a max over ~3000
eigenangles, combinatorially impossible. The real defect: principal SU3 eigenangles sum to ±2π
(a Z₃ centre twist), so `exp(iX)=e^{±i2π/3}U` (wrong action). Fix: traceless Hermitian log with the
nearest-0 angle shifted by ∓2π, giving `exp(iX)=U` **exactly for every config**. All configs are now
valid for the action derivatives; the residual twist (2.49% of links) is quantified, not discarded.

| case | N | valid | G″(KDE) | 1/Var | E[H]−Var | E[Φ″] | W | (E[H]−Var)+W |
|---|---|---|---|---|---|---|---|---|
| SU2 L=4 β2.5 | 3000 | 3000/3000 | **+0.281±0.044** | 0.272 | −0.02 | +0.34 | +0.22 | +0.20 |
| SU2 L=6 β2.5 | 1200 | 1200/1200 | **+0.276±0.071** | 0.272 | +0.07 | +0.33 | +0.24 | +0.31 |
| SU3 L=4 β5.7 | 2000 | 2000/2000 | **+0.425±0.083** | 0.404 | +0.02 | +0.46 | +0.28 | +0.30 |
| SU3 L=6 β5.7 | 500 | 500/500 | **+0.329±0.151** | 0.380 | −0.06 | +0.46 | +0.31 | +0.25 |

Errors: bootstrap SE (pooled over 6 modes; mode-to-mode scatter included).

**Resolutions.**
1. **Volume-robust?** Yes. SU2: 0.281(L4)→0.276(L6), `1/Var`=0.272 at both. SU3: 0.425→0.329
   (L6 500 cfgs, ±0.15), `1/Var`=0.40→0.38. G″ stays positive and O(0.3–0.4); the soft-mode
   `1/Var(k)` is unchanged with volume — no collapse toward 0.
2. **W or E[H]−Var?** W. `E[H]−Var`∈[−0.06,+0.07]≈0 in every case; `E[Φ″]`=+0.34(SU2)/+0.46(SU3)
   is volume- and mode-independent to ±0.01; `W`=+0.22–0.31 carries it. Pure-Jacobian control
   (β=0, S=0): SU2 `V_eff″`=0.248 vs `E[Φ″]−Var(Φ′|0)`=0.246 (exact closure); SU3 `V_eff″`≈0.48–0.55
   with `E[Φ″]`=0.45. The geometric term alone is positive — this is the decisive mechanism.
3. **Branch-safety fixed?** Yes — corrected traceless log makes `exp(iX)=U` exact for 100% of
   configs (all four runs show N/N valid); twist systematic quantified at 2.49% of links.

**Caveats (honest).** (i) The identity `(E[H]−Var)+W` closes to within ~0.06–0.12 of the KDE `G″`
(≈1–2σ of the component errors). It is exact at β=0 (SU2: 0.248 vs 0.246); at β=2.5/5.7 it closes
to ~1σ for SU2, with a larger residual for SU3 whose log amplitude carries the chart-cut — the
*sign* and the W-dominance are robust.
(ii) `E[Φ″]`/`W` use the flat (exponential-chart) action Hessian (`E[H]`≈1.52/2.03), not W7's
group-geodesic ray Hessian (2.46); the flat one is the one in the marginal-curvature identity.
(iii) `G″` is amplitude-definition dependent beyond quadratic order (log vs sin differ), as any
lattice gluon-mode curvature; we use the standard `A∝logU`.

**Files.** `w7/rerun_core.py` (log, Jacobian, KDE/conditioning), `rerun_run.py` (driver),
`rerun_analyze.py` (pooled table); data `w7/rerun_su{2,3}_L{4,6}.npz`; logs `w7/su3_L{4,6}.log`.

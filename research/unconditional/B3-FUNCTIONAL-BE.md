# B3-FUNCTIONAL-BE — Bakry–Émery Γ₂ on the conditioned coarse measure: the criterion is *exactly* H1 (lossless on flat space), and the Wilson cone is **not** compensated by Φ's log-cut convexity — anti-localized (floor `N/12` at `F≈0`, surge at `d→2π` where `F≈0`) plus a coarse-slice variance it cannot absorb

> **判词 (verdict).** Γ₂ on `ν=P_*μ` (full `V=S_W+Φ`) invents no route: on the flat coarse space `Γ₂(f,f)=‖Hess f‖²+Hess S_eff(∇f,∇f)`, and `Hess S_eff=G″` **exactly** (k0-CONSTRAINT), so `Γ₂≥ρ ⟺ G″⪰ρ ⟺ m_eff²(k→0)≥ρ ⟺ H1` — sharp, zero slack. The real question — is the **Wilson cone-negative exactly compensated by Φ's convexity** `Φ″⪰1/6` on the coarse slice? — answers **NO**, two faces of one object (`G″Var(a)=1`). *Microscopic/fiber*: `Φ″`'s usable floor is `N/12`, the cone `2g|F|~1.1β`; Loewner `−H⪯Φ″+k²` is **false**, spectral ratio `ρ≈3–16`, ~10× short (true only for `β<N/(12·1.1)≈0.15`). *Coarse slice*: bare `PᵀHP` min `−0.3` is *comparable* to `W_geo≈0.32` — **not** a magnitude gap — but the conditioning variance `Var(V′|0)` subtracts, tipping it (`slack W_geo−|X|≈−0.08`). `Γ₂≥ρ>0` does **not** close; it dies at the `k→0` cone mode for exactly H1.

**§1 The criterion is lossless.** The coarse generator reversible w.r.t. `ν` is `L_ν=Δ_c−∇S_eff·∇_c` on the flat pushforward of `P` (orthogonal onto a linear subspace of `X=log U`). Take `f` affine (`Hess f≡0`, killing the always-positive `‖Hess f‖²` term): `Γ₂(v,v)=Hess S_eff(v,v)=G″(v,v)`. So Poincaré(`ρ`) on `ν` is **equivalent** (no loss) to `G″⪰ρ` — flatness removes BE's only possible advantage over the identity `Hess V_eff=G″`. Distinct from W6's *full-space* `Ric+Hess S_W` (min `−2.0`, a genuine no-go): that quantifies all of `TM`; here the *marginal* curvature `G″`.

**§2 The full-V curvature, two faces.** `G″(0)=E[V″|0]−Var(V′|0)=E[H|0]+E[Φ″|0]−Var(G+Φ′|0)=Δ_dyn+W_geo+X`, `H=Hess S_W`, `W_geo≥0` (Prékopa), `E[Φ″|0]≥N/12`. Unconditioned at the `P`-invariant vacuum `E[H]=E[f²]` (equipartition, `Δ_dyn=0`), so `G″` **is** the conditioning defect — *conditioning is 100% of the gap*. The faces connect: high-momentum fiber negatives (90–97% of W6) Schur-integrate into `Var(G|0)`, then Ward `G″·Var(a)=1` turns variance into curvature.

**§3 The key question — compensation on the coarse slice.**
- **(i) Shallow cell `2g|F|≤N/12`:** YES. `K^V=k²−2g|F|+Φ″≥k²−2g|F|+N/12≥0`; BL fires, `Γ₂≥0`, all modes incl `k→0`, all `L`. A real theorem — the provable curvature.
- **(ii) Open cell `N/12<2g|F|≲μ²`, high momentum:** NO. Pointwise Loewner `−H⪯Φ″+k²` fails: uniform `c<1` is **false** (WGEO-DOMINATION), `ρ≈3–16` at β=2.5. Magnitude `N/12≈0.17–0.25` vs `2g|F|~O(1)`.
- **(iii) Localization — the decisive reason.** `Φ″`'s **floor** `1/6` is at `d→0` (eigenvalue collision, small field); the **surge** `g″→∞` is only at the log-cut `d→2π` (links→center `z`), where `U_□=zzz⁻¹z⁻¹=1` so **`F≈0`**. The cone `∝|F|` lives at per-link phases `O(π)<2π`, **away** from the cut, where `Φ″=O(1)=N/12`. The compensating convexity is **anti-correlated** with `|F|` — large where the cone is small. Structurally the wrong corner.
- **Conditioning makes it worse.** `−Var(V′|0)` is a *negative-definite* addition absent from pointwise `Hess V`; compactness bounds it (finite) but never flips its sign. The fiber's own Lichnerowicz `Γ₂` bounds the **fast** fiber relaxation — the wrong object; the coarse stiffness is transverse (PROOF-COARSE-POINCARE (b)).

**§4 Where `Γ₂≥ρ>0` dies — the cone negative.**
```
∃ρ>0: ∀U′, ∀ coarse v⊥orbit,  Γ₂^ν(v,v)=G″(0)[v,v]
   = E[H+Φ″|0][v,v] − Var(G+Φ′|0)[v,v] ≥ ρ‖v‖²,  uniform (L,a);
⟺ G″⪰ρ ⟺ m_eff²(k→0)≥ρ ⟺ Var_ν(φ)≤(1/ρ)ℰ_ν(φ) = H1.
Dies at v₀ = the k→0 Savvidy-negative coarse mode (PᵀHP min≈−0.3):
   E[Φ″|0]≈N/12 (floor; log-cut surplus measure-suppressed, F≈0 there)
   < −E[H|0] + Var(V′|0) + ρ.   [fiber ρ_spectral≈3–16; coarse slack −0.08]
```
The missing piece is the same conditioned-variance bound (◆): `Var(V′|0) ≤ E[V″|0]−ρ`.

**§5 Five-point verdict.**
1. **Nonperturbative — PASS.** Exact pushforward criterion, Prékopa/BL, Savvidy spectrum; no loop.
2. **Sign-controlling = LOWER bound `m_eff²(k→0)≥λ` — FAIL at the cone.** Shallow gives `Γ₂≥0` (true); open/`k→0` indefinite (`E[Φ″]≈N/12` vs `−E[H]≈0.3`+`Var`).
3. **Non-circular — FAIL at closing.** `Γ₂^ν=G″=H1` verbatim (route B); route A (push LSI through `P`) needs the 2.5 gap; reweighting a verified no-go. A faithful rename, zero slack.
4. **Non-SUSY + fails free field — PASS (caveat).** Continuum free field `Φ=0`, zero mode flat ⇒ `G″=Γ₂=0`, `λ=0`. Caveat: the lattice floor `N/12` is present at β=0 (`0.248`) — a compactness artifact absent in the continuum — so positivity must be attributed to the **conditioning** (`E[H|0]=+0.22`, `W_full≈0.54`), not the floor, or it certifies a kinematic mass as the gap.
5. **Uniform lower bound `≥λ` — FAIL.** Exactly H1.

**Sharpest line.** BE `Γ₂` on `ν` is *lossless* (flat coarse space) hence *exactly* `G″⪰ρ=H1` — not a new functional, a faithful dialect. The key question answers **no**: Φ's convexity does not compensate the Wilson cone on the coarse slice, for a *structural* reason (anti-localization: floor `N/12` at `F≈0`, surge at the log-cut `d→2π` where `F≈0`, cone `∝|F|` at `O(π)` sees only `N/12`) and a *measure* reason (conditioning subtracts `Var`, tipping the coarse slice by `−0.08`). The compact fiber's positivity bounds the fast fiber relaxation, not the slow coarse stiffness; the cone negative is in the dynamical direction where compactness gives nothing. Residual: one conditioned-variance bound = H1, unmoved.

*Refs.* OTHER-H1-POINCARE (§2B BE≡H1, two-route no-go); PROOF-SAVVIDY-WGEO-DOMINATION (`c<1` false, `ρ≈3–16`, log-cut anti-correlation, β<0.15); OTHER-SAVVIDY-SIGN (coarse = variance ceiling, not magnitude); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`, (◆)); ABC3-NTRANSFER §0 (`G″Var=1`); ABC-NEW-INEQ (Ward, `E[H]=E[f²]`); ABC2-SYNTHESIS (open cell = H1; 10×); PROOF-COARSE-POINCARE (PᵀHP −0.3; fast≠slow); W6-BE-at-2.5 (full-space min −2.0, 90–97% fiber); k0-CONSTRAINT-EFFPOT (`V_eff″=G″`, β=0 0.248/0.246); FORGE5-CONE-VARIANCE (conditioning = gap); PROOF-SAVVIDY-EXACT-SIGN (`λ_min=−g|B|`); Bakry–Émery 1985; Prékopa 1971; Brascamp–Lieb 1976; Holley–Stroock 1987; Savvidy PLB 71 (1977) 133.

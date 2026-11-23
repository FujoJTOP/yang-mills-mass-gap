# MS-GAP-EQUATION — the Ward/DSE gap equation is exact and acid-clean; its only *provable* fixed point is gapless, and the positive branch's bifurcation needs a k→0 Savvidy tachyon that at β=2.5 lives at high momentum

> **判词 (verdict).** Write the gap as the fixed point `κ=G″(0)=m₁²` of the exact Ward–Stein/DSE map `κ=F(κ)`, `F(κ)=(1−η)/Var(a;κ)`, `Var(a;κ)=∫wC_κ`, `C_κ=[k̂²+Π]^{-1}`, `Π=loop[C_κ]`. The gapless branch `κ=0` is **always** a solution and is the **contraction attractor** (free field `F≡0`, acid test clean). A second, positive, `(L,a)`-uniform fixed point exists **iff `F'(0)>1`** — a pitchfork, i.e. a linear IR instability of the gapless vacuum. That criterion is **not met at the coarse mode**: the Savvidy tachyon that would supply it is transverse-Landau-gapped and, at β=2.5, its negatives are 90–97% high-momentum while the low-`k` projection is positive. The two candidate generators — dynamical tachyon condensation (absent at `k→0`) and the non-dynamical Haar floor `W_geo≥N/12` (real at β=0; at 2.5 folds into `W>−Δ_dyn`=(◆)) — both terminate in H1. It banks β<1/48 (P3) and β=0 (Haar); **at β=2.5 the provable branch collapses to gapless.** Scorecard: #1,#4 pass; #2(2.5),#3,#5 fail. **Not 5/5; Clay 0%.**

## §1 The exact identity and the equation
Flat chart `X=log U`, mode `a=⟨v,log U⟩`, `|v|=1`, `μ∝e^{−V}dX`, `V=S_W+Φ`, `f=∂_aV`. Integration by parts on the compact fiber (no boundary) gives the **Ward–Stein identity**
```
E_μ[f a] = 1   (all β, nonperturbative).
```
With marginal `p(a)`, `m(a)=E[f|a]`, `∂_a log p=−m`, `m′(a)=E[H|a]−Var(f|a)=G″(a)`, this is `∫am(a)p da=1`; reflection symmetry (`E[a]=0`) gives
```
G″(0)·Var(a) = 1  ⟺  marginal exactly Gaussian;
   generally  κ·Var(a) = 1−η,  η = non-Gaussianity defect (0 ⇔ Gaussian).
```
So `κ:=G″(0)=m₁²` is the inverse coarse variance plus defect, and the **self-consistent gap equation** is the fixed-point problem
```
κ = F(κ) := (1−η(κ)) / Var(a;κ),   Var(a;κ)=Σ_k |w(k)|²/[k̂²+Π(k;κ)].
```
`f₁²` is the atom weight in the mode's spectral measure (`Var(a)∋f₁²/κ`); the positive fixed point **is** the atom (`κ=m₁²>0`, `f₁²>0`); the gapless branch carries none (`f₁²=0`).

## §2 Fixed-point structure
- **Gapless `κ=0` always.** Free field `Π=0`, `Var(a)~L²→∞` (IR zero mode), `F(0)=0`. On the torus `κ~1/L²→0` — the trivial geometric "gap," not `L`-uniform.
- **Positive `κ>0` iff `Var(a)` bounded uniformly** `⟺ Π(0)≥λ` — **verbatim H1** (`G″≈1/Var` to 2%, PROVE-Gpp §5). The equation does not manufacture it; it is equivalent to it.
- **The Ward moment does not force it.** `∫amp=1` is a *global* normalization, compatible with `m′(0)<0`: odd `m=−εa+a³` gives `∫amp=E[a⁴]−εE[a²]=1>0` while `m′(0)=−ε<0`. The sign `κ=m′(0)` lives in the **IR curvature of the force response**, unconstrained by the moment — the gap-equation face of "moment ≠ atom."

## §3 Steelman — spontaneous generation — and the bifurcation criterion
The classical mechanism is **tachyon condensation**: Savvidy's `V_1(B)` has negative curvature at `B=0` (`λ_min=−g|B|`), the perturbative vacuum is unstable, and a condensate `⟨F²⟩>0` yields `m₁~√(2g|B|)>0`. In gap-equation language this is a **pitchfork**: `F'(0)>1` bifurcates a positive `κ*` — the sharp, falsifiable criterion. **It fails on two counts:**
1. **The tachyon is not a `k→0` mode.** Spectrum `λ_{n,s}=k_∥²+(2n+1)g|B|+2s gB` — lowest Landau level: `k_∥=0` but **transverse-Landau-gapped** at `√(g|B|)`; it modulates transversally, breaking translation. At β=2.5 the coarse plane-wave subspace is **positive** (+0.40β, 0% negative) and negatives are 90–97% high-momentum (W6-BE). `F'(0)>1` at the coarse mode is **not present**; a fast-mode gap cannot gap the slow mode (PROOF-IR-H1 §4).
2. **The measured positivity is Haar, not dynamical.** `G″≈1/Var` with `E[H]−Var≈0` (equipartition), so the surviving curvature is `W_geo+X`, the Jacobian term (PROVE-Gpp; CAVITY-MEAN-FIELD). Real (`W_geo≥N/12`, Prékopa) but **non-dynamical**; at 2.5 its survival needs `W_full>−Δ_dyn`, i.e. `Var(G|0)≤E[H]+W−λ` — **(◆)=H1**.

Even a proven `κ>0` from the DSE is gauge-fixed, sign-indefinite (`∫ρ_A=0`; gluon ∉ BRST cohomology): dynamical mass, not the physical gap (PROOF-FRG-DSE-GAP-EQUATION).

## §4 Acid test + scorecard
Free massless: `Π=0⇒F≡0⇒κ=0` unique; no manufactured gap. ✓
1. **Nonperturbative — PASS** (exact Ward/Stein; no expansion).
2. **Sign `κ=m₁²≥λ` at 2.5 — FAIL.** Provable only gapless.
3. **Non-circular — PASS as identity, FAIL at closing** (`κ>0 ⟺ Var(a)≤1/λ`=H1; vertex pole = order parameter).
4. **Non-SUSY + fails free field — PASS (cleanest).**
5. **Uniform `(L,a)` — FAIL** at 2.5; only β=0 `N/12`, β<1/48 (P3).

**Sharpest line.** The gap equation is exact with a clean free-field degeneracy (`F≡0`) and genuine positive fixed points at strong coupling (β<1/48 P3; β=0 Haar `N/12`). At β=2.5 spontaneous generation needs `F'(0)>1` — a `k→0` Savvidy tachyon — but the tachyon is transverse-Landau-gapped, its negative weight in the high-momentum block the coarse mode does not couple to; the Haar-floor generator is non-dynamical, its survival (◆). The Ward moment `∫amp=1` is a global normalization that cannot force IR curvature. **Provable outcome: collapse to the gapless branch; the positive branch (measured +0.51, W7) is H1. 5/5 not achieved; Clay 0%.**

*Refs.* ABC-NEW-INEQ §1 (Ward `G″·Var(a)=1` ⟺ Gaussian; defect); FORGE5-CONE-VARIANCE (SD `E[H]=E[f²]`, equipartition); PROOF-FRG-DSE-GAP-EQUATION (DSE closure, transfer failure, Oehme–Zimmermann `∫ρ_A=0`); ATOM-FRG-mass (transmutation `m²~Λ²e^{−c/g²}`, three-gluon vertex order parameter); BANACH-FIXED-POINT (gapless contraction attractor, zero mode neutral, bifurcation=order parameter); PROOF-SAVVIDY-EXACT-SIGN (`λ_{n,s}` ladder, `λ_min=−g|B|`); W6-BE-at-2.5 (negatives 90–97% high-momentum, coarse plane-wave +0.40β 0% neg); PROOF-IR-H1 §4 (fast-mode cannot gap slow-mode); PROVE-Gpp-positive §5 (`G″≈1/Var`; H1⟺`Var(a)≤1/λ`); CAVITY-MEAN-FIELD (positivity is Haar, `Δ_dyn≈0`); k0-CONSTRAINT-EFFPOT (`V_eff″=G″`, β=0 `N/12`); PROOF-JACOBIAN-SAVES (`W_geo≥N/12`, (◆)); GAP-G2-mass (edge/localization, moment≠atom); PROOF-COARSE-POINCARE (H1); PROOF-P3-UNCONDITIONAL (β<1/48). Savvidy PLB **71** (1977) 133; Prékopa 1971; Shen–Zhu–Zhu, CMP **400** (2023) 805.

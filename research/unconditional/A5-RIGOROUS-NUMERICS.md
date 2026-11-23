# A5-RIGOROUS-NUMERICS — a certified finite-region `G″>0`: exact architecture, positive point-value margin, three walls

> **判词 (verdict).** There is an **exact finite-region certificate architecture** and — using the *proven* pointwise floor `E[Φ″|0] ≥ N/12 = ¼` (Lemma B) instead of the weaker `W_geo≥0` — it carries a **positive point-value margin** on the existing SU(3) β=5.7 data: **+0.092 (L=4) / +0.044 (L=6)** for the panel-averaged softest mode at `k=2π/L`. So "measured positivity + an a-priori floor" *does* assemble into a positive lower bound at point values. But the certificate is **not a theorem**, for three independent reasons. **(1) Error bars ≠ enclosures:** the MC block-bootstrap σ is *statistical*; it never yields a rigorous interval, so the three scalars the certificate needs — `bracket`, `Var(Φ′|0)`, `X` — are estimated, not bounded. **(2) The a-priori route is vacuous:** the only deterministic ceiling on `|X|` is Cauchy–Schwarz `2√(VarG·VarΦ′) = 1.06`, **24× the whole 0.044 budget**, and the range/perturbation bounds are astronomically worse. **(3) The validated computation is infeasible:** the enclosure it needs is a ratio of integrals over `G^{4L⁴}` (**8192 real dims** at L=4), nonsmooth (branch cuts, eigenvalue collision where `g′=−cot(d/2)+2/d→∞`), sign-indefinite (Savvidy) — beyond sparse-grid / TT-quadrature / certified transfer-matrix. **And even if closed, it would not be a gap result:** at the *finite* mode `k=2π/L` the **free massless field already has `G″=4sin²(k/2)=2 (L4) / 1 (L6) > 0`**, so finite-k positivity is *not* a YM-vs-free discriminator — it does not engage the charter's free-field acid test (which lives at `k→0`, where the free field gives `G″=0`). The honest object this route proves is a **finite-(β,L,k>0) lattice-curvature statement**, conditional on one validated-numerics input that is exactly (◆)/H1 in Cauchy–Schwarz clothing.

## §1 The exact finite-region object and the certificate identity

Fix β=5.7, L∈{4,6}, and a transverse plane mode `v` at `k=2π/L` (unit-normalized in the real algebra over `4L⁴` links). Let `a=⟨v,log U⟩` (branch-safe log, rerun_core `su3_log_herm`). Then
`G″(β,L,k) := d²/da²[−log P(a)]|_{a=0}` is a **finite, well-defined** number. In the flat chart the exact identity (rerun_core; PROOF-JACOBIAN-SAVES §1, W7-Gpp0) reads
```
G″(0) = bracket + W_geo + X ,
bracket = E[H_vv|0] − Var(G_v|0) ,   W_geo = E[Φ″|0] − Var(Φ′|0) ,   X = −2Cov(G_v,Φ′|0) ,
```
with `H_vv=d²S/ds²`, `G_v=dS/ds` along `X→X+sV`, `Φ=−log J` the Haar log-Jacobian. **Proven (Lemma B, PROOF-JACOBIAN-SAVES §2):** the pair function is convex, `g″(d)=½csc²(d/2)−2/d²≥1/6`, giving the pointwise floor `Φ″[v]≥(N/12)‖v‖²`, hence under **any** measure (conditioned or not)
```
E[Φ″|0] ≥ N/12 = 1/4    (SU(3), unit mode) .                                          (FLOOR)
```
The certificate is then a three-scalar inequality:

> **Theorem A5′ (conditional certificate).** Let `b* ≤ bracket`, `u* ≥ Var(Φ′|0)`, `x* ≤ X` be rigorous bounds. Then
> `G″(β,L,k) ≥ b* + ¼ − u* + x* =: λ .`   If `λ>0`, the panel-averaged softest-mode curvature at fixed (β,L,k=2π/L) is **rigorously positive**.

**This is exact arithmetic**; the whole question is whether `b*,u*,x*` can be made rigorous. Measured point values (pool `_pool_su3_L{4,6}.npz`; `a5_certificate_margins.py`) give `λ_point = +0.092 (L4) / +0.044 (L6) > 0` — so the architecture has a **thin positive margin**, and the certified object must be the **6-mode panel average** (per-mode `λ_point` ranges −0.48…+0.72 at L6, so individual-mode certification is not even point-value-supported).

## §2 The three legs — what is proven vs measured

| β | L | N | bracket | E[Φ″] | Var(Φ′\|0) | W_geo | X | G″(KDE) | λ_point (¼ floor) |
|---|---|---|---|---|---|---|---|---|---|
| 5.7 | 4 | 2000 | +0.019 | 0.458 | 0.139 | +0.318 | −0.038 | **+0.425** | **+0.092** |
| 5.7 | 6 | 500 | −0.058 | 0.458 | 0.139 | +0.319 | −0.009 | **+0.329** | **+0.044** |

- **Proven (any measure, any k, any L):** `E[Φ″|0]≥¼` (FLOOR) and `W_geo≥0` (Prékopa/Φ-convex). These are the *only* rigorous positive numbers. Note the discriminator: FLOOR `=¼` is **sharper than** `W_geo≥0` because it does not pre-subtract `Var(Φ′)`; using it lifts the point-value certificate from PROOF-VALIDATED-GPP-SIGN's `−0.13/−0.27` to `+0.09/+0.04`. That lift is the **one new tool** this memo contributes: *reframe the certificate around the pointwise floor, not the marginal Prékopa sign.*
- **Measured, not bounded:** `bracket`, `Var(Φ′|0)`, `X`, and `E[Φ″]`'s exact value 0.458 (vs floor 0.25). All four are conditional moments at `a=0`, i.e. ratios of finite integrals. The margin `λ_point` is the **δ-budget**: the three rigorous enclosures `u*,−b*,−x*` may together consume at most **0.092 (L4) / 0.044 (L6)** and still certify `G″>0`. The measured widths are small (`Var(Φ′)=0.14`, `|bracket|=0.02–0.06`, `|X|=0.01–0.04`), so a *tight-but-rigorous* enclosure would close — **but none exists** (§3–§4).
- **Chart systematic:** the identity `bracket+W+X` closes to `G″(KDE)` only to `−0.13/−0.08` pooled (COMPUTE-GPP-SU3; per-mode closure −0.59…+0.49). This is an *estimator* artifact (KDE vs kernel-conditioned vs histogram give slightly different splits), not a failure of the identity — a rigorous implementation uses one fixed chart definition, so the systematic would not enter a true enclosure — but it flags that these scalars are delicate and the 0.04 L6 margin sits *below* the estimator noise.

## §3 Why "rigorous error bars" cannot come from the MC samples

Block bootstrap (block=20, 1500 resamples, PROOF-VALIDATED-GPP-SIGN) excludes 0 at 5.1σ (L4) / 2.2σ (L6) conservative — a genuine, strong *statistical* statement. But resampling returns an **estimate and a standard error, never an interval that contains the true conditional moment with certainty**. To pass from σ to a bound one needs one of: (i) a certified sampler (a proven mixing-time / spectral-gap bound → concentration), nonexistent for lattice gauge at β=5.7 on these volumes (that *is* the gap, PROOF-FSS / H1-T3); (ii) an a-priori deterministic ceiling on the integrand's fluctuation (a.s. range), which gives `Var(G)≤(range²)/4` with range `O(β·√(4L⁴))≈O(10²–10³)` → bounds `O(10⁴)`, vacuous; (iii) `e^{sup S}` perturbation of a Haar Poincaré inequality, with `sup S` up to `O(10³)` → constant `e^{10³}`, vacuous. **No route turns the σ into an enclosure.**

## §4 The certified computation: design and feasibility

The rigorous input Theorem A5′ needs is a validated enclosure of three conditional moments at `a=0`. Each is a ratio of finite integrals over the compact group `G^{4L⁴}` (`dU` Haar):
```
bracket, Var(Φ′|0), X  =  [∫ (…) e^{−S−Φ} δ(a) dU] / [∫ e^{−S−Φ} δ(a) dU] .
```
At L=4 this is a **8192-real-dimensional** integral. Enumeration of the design space of validated-numerics:

- **Validated (interval) quadrature / sparse grids.** Need smoothness and bounded mixed derivatives. The integrand is **nonsmooth**: `log U` has branch cuts (center twist, handled but discontinuous across it), `Φ` has eigenvalue-collision log-divergences (`g′=−cot(d/2)+2/d→∞` at `d→2π`). Sparse grids and tensor-train cross-interpolation both fail (no rigorous TT-rank bound for 4D plaquette couplings). **Infeasible.**
- **Certified transfer matrix / Peter–Weyl Feshbach (PROOF-COMPUTER-ASSISTED).** This is the project's one *rigorous* finite-lattice engine, but (a) it certifies a finite-volume **gap** (an *upper* bound `Δ_V≥Δ_∞`), not a marginal curvature, and (b) it already exceeds certified dense linear algebra at **24 links** (2³), while L=4 has 1024 links. **Wrong object + wall.** FORGE5-CERTIFIED-TN confirms certified PEPS contraction gives no operator-norm lower bound (`ε(D)` volume-exponential; sign-blind; FORGE5-SL-UV-CONE empty).
- **Deterministic a-priori bounds.** Cauchy–Schwarz on `X`: `|X|≤2√(Var(G|0)·Var(Φ′|0)) = 2√(2.01·0.139) = 1.06`, against a total budget `0.044–0.092` — **24× too loose**, and that is the *best* determinant bound (range/perturbation are far worse, §3). This is B4-VARIANCE-MAG's finding in the floor language: the sign lives in `corr(G,Φ′)≈0` near-orthogonality, which every magnitude bound discards.
- **The one place it closes rigorously — and it is not the target.** At **β=0** the measure is a *product* of one-link Haar: the marginal of `a` is the `4L⁴`-fold convolution of the single-link (8-dim) projection density, so `G″(0)` is **rigorously computable by 1-link validated quadrature + Berry–Esseen `O(1/√n)`**, giving `G″≥1/(n·σ_link²)−Ce^{−cn}>0` unconditionally. This is a real certified computation (light, feasible) — but β=0 is the infinite-coupling limit, and the small-β window is already the landed strong-coupling theorem (PROOF-P3-UNCONDITIONAL, `f₁²>0` at β<1/48). **It does not reach β=5.7.**

**Feasibility verdict:** the finite-region bound at β=5.7, L∈{4,6} is **not certifiable with current validated-numerics technology**. The exact missing input — a rigorous enclosure of `Var(Φ′|0)`, `bracket`, `X` good enough to beat the 0.044–0.092 budget — is, via `G″≈1/Var(a)` (PROOF-JACOBIAN-SAVES §2), **identical to the variance bound (◆): `Var(G|0)−E[H|0]<W_full`, uniform in (L,U′)`, i.e. H1 in curvature language.** The certificate does not dodge H1; it *is* H1, wearing a finite-dimensional costume.

## §5 Exactly what this proves — and what it does not

**Proven today (unconditional):**
1. `E[Φ″|0] ≥ ¼` for every mode, every k, every L (Lemma B) — a rigorous *positive* lattice number.
2. `W_geo ≥ 0` for every mode including k→0 (Prékopa).
3. The **conditional** Theorem A5′: *if* `b*,u*,x*` are supplied rigorously, then `G″(β,L,2π/L) ≥ b*+¼−u*+x*`.
4. A **β=0** (and small-β via P3) *unconditional* `G″>0` by 1-link validated quadrature — non-physical coupling.

**NOT proven — the exact boundaries:**
- **No `G″>0` at β=5.7 is certified.** The point-value margin (+0.09/+0.04) is positive but the enclosures are statistical; the a-priori and validated-numerics routes are vacuous/infeasible. To certify it would need the (◆)/H1 variance bound.
- **No `k→0`.** The mode is fixed at `k=2π/L` (`π/2, π/3`); nothing is said about the `k→0` limit, which is where the free field gives `G″→0` and the gap actually lives.
- **No `L→∞`, no continuum, no mass gap.** Two fixed volumes; no volume-uniformity; the continuum `a→0` along the RG trajectory (IN-9) is untouched. `Δ_V≥Δ_∞` makes any finite-V lower bound an *upper* bound on the continuum gap.
- **Fails the free-field acid test as a gap tool.** At finite `k=2π/L`, `G″^free(k)=4sin²(k/2)=2 (L4) / 1 (L6) > 0` — *the same sign, a larger value than the interacting 0.42/0.33*. A certified finite-k positivity therefore also holds for the free massless field and is **not** a YM-vs-free discriminator; it does not touch the `k→0` object (`G″^free(0)=0`) the charter's acid test protects. So even a fully rigorous finite-region bound would be a true statement *about lattice curvature*, **not a step toward the gap in the acid-test sense.**

## Sharpest line

The floor `E[Φ″|0]≥¼` is genuinely rigorous and turns measured positivity into a **positive point-value certificate `G″ ≥ +0.09 (L4) / +0.04 (L6)`** — a thin margin whose δ-budget (enclose three scalars to total width `<0.05`) is real. But three independent walls block it from becoming a theorem: the bootstrap σ is not an enclosure; the only a-priori ceiling (`|X|≤1.06` by Cauchy–Schwarz) is 24× the budget; and the validated 8192-dim integral is beyond interval/TT/transfer-matrix methods. The missing input is exactly the H1 variance bound (◆), and even if supplied the result would be a **finite-k** statement the **free field also satisfies** (`G″^free=2/1>0`), so it is a partial lattice theorem, *not* a gap discriminator. The honest gains are: the pointwise-floor certificate (sharper than `W_geo≥0`), the explicit δ-budget, and the clean statement that rigorous finite-*k* positivity is the wrong acid test.

---

*Scripts.* `research/unconditional/a5_certificate_margins.py` (reproduces §2 and the CS shortfall from `w7/_pool_su3_L{4,6}.npz`); underlying `w7/rerun_core.py` (`kernel_cond`, `phi_derivs_analytic`), `w7/kzero_analyze.py`. *Internal refs.* COMPUTE-GPP-SU3 (+0.425/+0.329, decomposition); B4-VARIANCE-MAG (CS 0.40>budget, sign-blind); PROOF-VALIDATED-GPP-SIGN (bootstrap is statistical; `W_geo≥0`-only gives −0.13/−0.27); PROOF-JACOBIAN-SAVES (Lemma B floor `N/12`; (◆)=H1); CONDITIONAL-LATTICE-GAP (H1); OUT-certified-computation, PROOF-COMPUTER-ASSISTED (8192-dim / 24-link wall, *upper*-bound direction); FORGE5-CERTIFIED-TN (no operator-norm lower bound); PROOF-P3-UNCONDITIONAL (β<1/48 `f₁²>0`). *External.* Prékopa, Acta Sci. Math. 34 (1973) 335; Brascamp–Lieb, J. Funct. Anal. 22 (1976) 366; Wilson, PRD 10 (1974) 2445.

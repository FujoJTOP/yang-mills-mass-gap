# PROOF-VALIDATED-GPP-SIGN — the bootstrap excludes 0; a *rigorous* finite-(β,L) bound does not follow from the samples

> **Verdict.** The coarse-mode `G″>0` on SU(3) β=5.7 survives as a **statistical** statement of exceptional strength, not as a proof. Block bootstrap (circular blocks = 20, 1500 resamples, 6 transverse modes at k=2π/L) excludes 0 at **5.1σ one-sided at L=4** (conservative single-mode SE) and **2.2σ at L=6** conservative — **6.3–12.7σ** when the SE is that of the pooled mean. That is *statistical*: resampling returns estimates and standard errors, never enclosures. The exact identity `G″=(E[H]−Var(G))+W_geo+X` is genuine and one term has a **proven** sign (`W_geo≥0` by Prékopa; `Φ″≥(N/12)‖v‖²` pointwise, Lemma B). But the positivity rides on the *measured* `W_geo≈+0.32`, whose rigorous content is only `≥0`: interval arithmetic that keeps **only** the proven part and treats `bracket,X` statistically yields the bound **−0.13 (L4) / −0.27 (L6)** — it does **not** certify `G″>0`. A certified finite-(β,L) `G″>δ` is **not achievable from these samples**.

## (1) Statistical bound — STATISTICAL, not proof
Reproduced exactly from `_pool_su3_L{4,6}.npz` (G″=+0.4248 / +0.3287):

| β | L | N | G″ | SE (joint / conservative) | z (joint / cons.) | 95% 1-sided LB (joint / cons.) | modes >0 |
|---|---|---|---|---|---|---|---|
| 5.7 | 4 | 2000 | +0.425 | 0.033 / 0.084 | 12.7 / 5.1 | +0.37 / +0.29 | 6/6 |
| 5.7 | 6 | 500 | +0.329 | 0.053 / 0.149 | 6.3 / 2.2 | +0.24 / +0.08 | 6/6 |

One-sided Gaussian p: `2×10⁻⁷` (L4, conservative), `1.4×10⁻²` (L6, conservative); all 1500 pooled replicates >0 on both volumes. KDE-bandwidth systematic ±0.006 (L4) / ±0.026 (L6); chart/branch closure `(bracket+W)−G″ = −0.13/−0.08` (~25%) is the dominant systematic and does **not** move the sign. This justifies "the sign is not a single-volume fluctuation" — and nothing more.

## (2) Validated numerics on the exact identity
Decomposition (pooled means): `bracket=E[H]−Var(G)=+0.019 (L4)/−0.058 (L6)` — a near-cancellation of two O(2) forces (`E[H]≈2.03`, `Var(G)≈2.01`, equipartition), sign-ambiguous; `W_geo=E[Φ″]−Var(Φ′)=+0.318/+0.319` (bootstrap CI +0.312…+0.325, stable to conditioning bw 0.20–0.40); `X=−2Cov(G,Φ′)=−0.038/−0.009` (small, uncertified). **Rigorous:** `W_geo≥0` (all modes, all k, uniform L); `E[Φ″]≥(N/12)‖v‖²` (≈0.25 floor; measured 0.458). **Not rigorous:** the *values* of `W_geo`, `bracket`, `X`. Interval arithmetic trusting only `W_geo≥0` + statistical CIs for `bracket,X`: lower bound **−0.13 / −0.27** — 0 not excluded. Trusting instead `W_geo`'s measured CI gives +0.18/+0.03, but that number is bootstrap-statistical, not certified. The sign lives in `W_geo`'s *magnitude*, which the theorem bounds from below by 0 only.

## (3) The honest limit — not certifiable from MC
Certifying `G″>δ` at fixed (β,L) needs one of: **(i)** exact/validated evaluation of the finite integral over `G^{4L⁴}` — 1024 links × 8 = **8192 real dimensions** at L=4, far beyond rigorous quadrature; **(ii)** a certified transfer-matrix enclosure (Peter–Weyl + Feshbach, PROOF-COMPUTER-ASSISTED) — the 2³ (24-link) threshold is already past dense linear algebra, so 1024 links is hopeless, *and* that engine certifies a finite-volume gap (an **upper** bound), not the marginal curvature; **(iii)** a rigorous MCMC spectral-gap/mixing bound plus concentration — nonexistent for lattice gauge at β=5.7 on these volumes; a Hoeffding step on `G` (sd 1.4, data range 11) already half-widths ~0.14 non-rigorously and the geometric a.s. bound is O(β), vacuous. The single missing inequality is exactly **(◆)** `Var(G|0)−E[H|0] < W_full`, uniform in (L,U′) ⟺ `Var(a)≤1/δ` ⟺ `m_eff²(k→0)≥δ` — H1.

**Sharpest line.** The bootstrap turns "+0.43±0.08, 6/6 modes" into a 5σ *statistical* exclusion of 0 — real, and still only statistical. The rigorous content of the measured `W_geo=+0.32` is the single inequality `W_geo≥0`; what survives on top of it is `bracket+X` — point value **−0.02 (L4) / −0.07 (L6)**, and −0.13/−0.27 at the 95% lower end — so no finite-(β,L) `G″>δ` is certified from MC alone. The gap is one quantitative variance bound — H1 — not a better dataset.

*Refs.* COMPUTE-GPP-SU3 (+0.425/+0.329, decomposition); W7-Gpp0 (exact marginal identity); PROOF-JACOBIAN-SAVES (`W_geo≥0`, Lemma B, (◆)); PROOF-COMPUTER-ASSISTED (certified spectral engine, finite-volume upper bound, link-count wall); Prékopa 1971; Brascamp–Lieb 1976; scripts `research/unconditional/w7/` (`_pool_su3_L{4,6}.npz`, bootstrap reproduced block=20/1500).

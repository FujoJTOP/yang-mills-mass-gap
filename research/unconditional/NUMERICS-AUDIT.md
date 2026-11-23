# NUMERICS-AUDIT — the L=8 curvature sign probe is a measure-dominated lattice number, not evidence for a continuum gap

**Verdict (blunt).** The campaign's reading — "measured `G″(k=2π/L)>0`, `W_geo`-dominated, therefore evidence FOR a mass gap" — does not survive audit. What the runs actually establish: at **one β, one ~1.4 fm box, and one finite momentum** per gauge group, the *amplitude marginal* `−log P(a)` is locally concave at `a=0` (`G″>0`), and that concavity is carried (94–107%) by the Haar-measure term `W_geo` — a β-independent, proven-`≥0` log-concavity of the gauge measure that the **free theory also satisfies**. The dynamical term that would carry a gap (the Wilson-action bracket) is measured to be statistically **zero**. **No rigorous positivity is established** (honest 95% one-sided lower bound −0.098 / −0.035). It is a *sign probe of a lattice curvature*, nothing about the continuum gap — and at finite k it does not even separate Yang–Mills from the free field, which is positive **and larger**.

## (a) Finite-k positivity is not a gap signal — it is the wrong object
For the free massless transverse mode `G″^free(k)=4sin²(k/2)`: positive at every finite k, →0 only as k→0. At the measured modes it is **larger than the interacting value in every case**:

| mode | `G″^free` | SU(3) β=5.7 | SU(2) β=2.4 |
|---|---|---|---|
| k=π/4 | 0.586 | +0.307 | +0.192 |
| k=π/2 | 2.00 | +0.396 | +0.262 |
| (L4 π/2 / L6 π/3) | 2.0 / 1.0 | +0.425 / +0.329 | — |

Two consequences. (i) Free and interacting are **both positive** — no discrimination. (ii) The interacting curvature is **below free** (ratios 0.20–0.52). If `G″` were an inverse propagator `k̂²+m²`, a gap would make it *larger* than free; it is *smaller*. So `G″` here is not a mass proxy — it is the flat-chart amplitude curvature, whose sign the measure sets. The charter's acid test lives at **k→0** (`G″^free→0`); at finite k the free field passes trivially, so this measurement **never engages** the acid test.

## (b) Zero rigorous positivity
Only two statements are theorems: `W_geo≥0` (Prékopa) and `E[Φ″|0]≥N/12` (Lemma B) — both **value-free**. The sign-carrying magnitude (+0.32) is *measured*. Keeping only the proven pieces, the 95% one-sided lower bound is **−0.098 (SU3) / −0.035 (SU2)**: 0 is not excluded. A5's "point-value margin" (+0.09/+0.04) needs three rigorous enclosures it lacks: bootstrap σ is not an enclosure; the Cauchy–Schwarz ceiling `|X|≤1.06` is 24× the budget; the validated 8192-dim integral is infeasible. The pooled +7.8σ is **statistical**, and worse — per mode the identity `bracket+W_geo+X` misses the gold-standard `G″` by a **mean |0.19| (SU3) / 0.16 (SU2)** — i.e. by ~the size of `G″` itself; and per-mode `bracket` ranges −0.36…+0.28 (SU3). So "94–107% `W_geo`" and "bracket≈0 / equipartition" are *pooled-average* attributions that hide large opposing per-mode errors, not pointwise facts.

## (c) Not a window, not a continuum
Single β per group (one lattice spacing, `a≈0.18 fm`, **no continuum extrapolation**), `L≤8` (box ≈1.4 fm — a volume *scan*, not a thermodynamic limit), softest mode `k=π/4` (never `k→0`). That is **one point in (β,L,k)**; nothing is uniform in L, at a→0, or at k→0. `G″` at k=π/4 has no continuum-mass meaning; the only bridge to `m_eff²(k→0)` is exactly H1 (`Var(a)≤1/δ`), which these runs do not touch.

## (d) What a rigorous continuum-safe gap measurement requires — and feasibility
It requires a spectral gap with a rigorous **lower** bound *uniform in volume* (the Clay condition `Δ_V≥Δ>0` for all V), taken at a→0 along the RG trajectory and V→∞. Numerically that is the standard lattice-glueball extraction (gauge-invariant correlators, several spacings, continuum + finite-size extrapolation → ~1.7 GeV SU(3) 0⁺⁺) — that is **evidence, never proof**, and it is *not* the amplitude curvature measured here. Rigorous certification of a continuum gap **is the Clay problem**: infeasible; the project's own certified transfer matrix only gives the finite-V gap in the *upper*-bound direction (`Δ_V≥Δ_∞`) and caps at ~24 links versus 1024 at L=4. Feasible → no; numerical → yes, and long done elsewhere, but not what these runs did.

**Sharpest line.** The runs establish a **sign probe of a Haar-measure log-concavity** on one L=8 box: `G″>0`, `W_geo`-dominated, action-bracket ≈0, statistically real but rigorously inconclusive (LB ≈ −0.10), *below* the free-field value, at finite k. That is a lattice-curvature **consistency check — not evidence for a continuum mass gap**: the dynamical part that would carry any gap is measured flat, and the surviving sign belongs to a term the free measure satisfies too.

---
*Audited against:* `COMPUTE-GPP-K0-CONTINUUM.md` (+0.307/+0.192; LB −0.098/−0.035), `COMPUTE-GPP-SU3.md` (+0.425/+0.329; "the flat Haar measure, not the Wilson action, supplies the curvature"), `A5-RIGOROUS-NUMERICS.md` (finite-k free value 2/1; CS 24× budget), `B4-VARIANCE-MAG.md` (equipartition bracket≈0, CS sign-blind). Raw data re-read: `w7/k0_su{3,2}_*_results.json` (pooled `W_frac_of_Gpp`=0.935/1.07; per-mode closure spread). Free values recomputed: `4sin²(k/2)`.

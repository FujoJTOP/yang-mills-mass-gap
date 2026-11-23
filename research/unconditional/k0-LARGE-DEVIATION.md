# k0-LARGE-DEVIATION — the LDP / Laplace rate of the volume-averaged field: does its curvature at the minimum give `Var(a_0)≤1/λ`?

> **Verdict (判词).** Switching the target from the atom `f₁²` to the **average** `a_0` is *correct and does* dodge both objections that killed the DV route on the residue (INEQ-R5): the average is **lsc-friendly** (a rate function of an average, being `Λ*`, is *automatically* convex and lsc) and its LDP bounds open neighbourhoods on **both** sides — so "usc atom vs lsc sup" and "LDP runs uphill" do **not** apply here. But the LDP does **not** deliver #5. The rate curvature at its minimum is an **identity**, not an inequality: `I''(0)=1/Var_μ(a_0)=G″_∞=m_eff²` (Gärtner–Ellis / Varadhan, second order). Hence "curvature ≥ λ" is *verbatim* "`Var≤1/λ` `⟺` `G″≥λ`" = **H1**: the rate function does not *bound* the curvature, it *names* it. Convexity is free but yields only `I''≥0` (the trivial `Var<∞`/`Var≥0`), never the uniform modulus `λ`. **1 equivocal, 2 free-but-trivial, 3 circular at closing, 4 compact-false-positive, 5 = H1.**

## The object, and the identity that is an equality
`a_0=⟨v_∞,X⟩`, `v_∞` the `k→0` mode; the volume average is `\bar X=a_0/√V`, so `Var(\bar X)=Var(a_0)/V`. The scaled pressure `Λ(θ)=lim_V V⁻¹log E_μ e^{V⟨θ,\bar X⟩}` is convex (a limit of log-Laplaces); `I=Λ*` is the convex, lsc rate. When `Λ` is twice differentiable at `0`, Gärtner–Ellis gives the exact dictionary
```
Λ''(0)=Var(a_0)=1/G″_∞ ,   I''(0)=1/Λ''(0)=G″_∞=m_eff² .
```
So the brief's identification — *rate curvature at the minimum = gap, positivity = variance bound* — is **exactly true**: `I''(0)=m_eff²`, `I''(0)≥λ ⟺ Var(a_0)≤1/λ`. But an equality cannot bound. To *use* `I''≥λ` you must supply `I''≥λ`, which is `G″≥λ`.

## Why the free structure cannot supply it
**(2) Convexity is real but sign-trivial.** `I` convex is robust to Savvidy (unlike Bakry–Émery) — a genuine gain — but it gives only `I''(0)≥0` (`Var<∞`). Strictness with modulus `λ` is *not* implied: convex `I` is compatible with a flat minimum (`I''=0`, massless) and a kink (source first-order transition); distinguishing strict from flat *is* the gap.
**(Direction.)** The moderate-deviation bound `P(⟨θ,√V a_0⟩>t)≲e^{−V t²/(2Λ''_θθ)}` carries `Λ''(0)=Var(a_0)` on its **right side** — encoding "`Var≤1/λ`" presupposes the variance bound; the free half gives only the variance **lower** bound.
**(4) Compact false positive.** On the torus `|X|≤π`, `Var(a_0)≤π²` by compactness *regardless of gap*: `I''(0)≥1/π²>0` even for a massless compact theory — the acid test fails on the lattice where β=2.5 lives (as FORGE5-SL-MARTINGALE #4). Massive-`ℝ^V` gives `I''=m²`, massless `0`: correct only off-compact.
**(5)/(3) Uniformity and data.** Even granting `I''(0)>0` at fixed `(L,a)`, `λ` independent of `(L,a)` is the volume-uniform LSI/Poincaré — the continuum/single-phase wall `(CP)` (PROOF-CONCENTRATION §4). And existence of `I''(0)` (`Var<∞`) *is* the gap. At β=2.5 the coarse proxy is `G″≈−0.165<0` (VAR-FINAL): the putative minimum is *concave*.

## Five-point verdict (honest)
1. **Nonperturbative — EQUIVOCAL.** `Λ(s)`, the uniform-source pressure, exists (finite-range Gibbs, subadditive). But *twice*-differentiability at `0` is not free; it is `Var(a_0)<∞` = the gap.
2. **Sign control — FREE, but sign-trivial.** "`I` is convex" (`I''≥0`) is value-free and Savvidy-robust — the one real gain. It does not give `I''≥λ`.
3. **Non-circular — FAIL at closing.** The dictionary is exact, but the closing input `I''≥λ` is `G″≥λ`; the MDP tail bound assumes the variance bound.
4. **Non-SUSY / acid test — FAIL on compact.** Free massive returns `m²`, massless returns `0` ✓ off-compact; compactness forces the false positive `Var≤π²`.
5. **Uniform `m_eff²≥λ` — FAIL.** Coefficient is `G″=H1`; uniformity is `(CP)`; data `G″≈−0.165`.

**Sharpest line.** The average is the right semicontinuity class — this genuinely escapes INEQ-R5 — but the Legendre machinery it buys computes the gap's *value* by an equality (`I''(0)=1/Var=G″`), not a *bound*; and the only free inequality, convexity, is the trivial half (`Var≥0`), while the needed half (`Var≤1/λ`) is exactly the uniform strict-convexity coefficient = the volume-uniform LSI. **#5 = H1, renamed as rate curvature.**

*Refs.* Donsker–Varadhan, *Comm. Pure Appl. Math.* **28** (1975) 1; Gärtner, *Th. Prob. Appl.* **22** (1977) 24; Ellis; Varadhan; Bryc (Legendre duality); Dupuis–Ellis (1997). Internal: INEQ-R5-donsker-varadhan (atom usc vs lsc sup — evaded by the average); PROOF-STOCHASTIC-LOCALIZATION (`A_0` reciprocal `=1/Var`); INVENT-VARIANCE-INEQ (`G″=S−D`); PROOF-JACOBIAN-SAVES (`G″≈1/Var`); PROOF-COARSE-POINCARE (H1=`m_eff²`); PROOF-CONCENTRATION (§4 volume-uniform LSI=gap); IDEA-finite-size (d, LD rate curvature = `G″`); FORGE5-SL-MARTINGALE (#4 compact `Var≤π²` false positive); VAR-FINAL (`G″≈−0.165`); INVENT-TRANSPORT-INEQ (`(CP)`).

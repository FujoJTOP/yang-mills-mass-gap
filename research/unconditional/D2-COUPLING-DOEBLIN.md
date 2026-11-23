# D2-COUPLING-DOEBLIN — the Doeblin / martingale coupling on the coarse slice: a genuine dynamical species whose contraction rate is *identically* the marginal curvature `G″`, hence sharp-equivalent to H1, with the Savvidy tachyon the exact missing drift

> **2026-11-23 · NEW SPECIES D2 (dynamical / coupling).** Forge a two-copy coupling of the equilibrium chain and ask whether the **coarse amplitude** `a=⟨v,X⟩` contracts by `(1−λ)` per transfer step — a Wasserstein/martingale **drift** (Doeblin / Hairer–Mattingly / Meyn–Tweedie) giving Poincaré `λ`. **Verdict.** The coupling is a genuine non-static species and it **passes the acid test** — yet its contraction rate is *identically* the coarse marginal curvature `b′(0)=G″` (von Renesse–Sturm / Bakry–Émery: `W₂`-contraction(λ) ⟺ `Hess V_eff⪰λ`, **sharp, zero slack**). The Savvidy cone **blocks the drift**: the anti-self-dual tachyon is exactly the `D>0` in `G″=S−D`, it sits at `k→0` (`PᵀHP` min ≈−0.3), and at β=2.5 it is **typical** (equipartition `Var(G|0)≈E[H|0]≈2`), not a rare set a coupling can marginalize — so the *pointwise* drift is expanding on the cone and only the *measure-averaged* `G″` survives, whose positivity **is** H1. Full-chain Doeblin adds an `e^{−cβ}` metastability loss. **1 PASS; 2,3,5 FAIL; 4 PASS (fails free field cleanly).**

## 1. The two coupled copies, the drift, and the exact rate
Equilibrium chain `P` for `μ∝e^{−S_W}dU` on `𝒜=G^{E}` (heat-bath = electric heat kernel `>0` × `e^{−a·magnetic}`, PROOF-GAUSS-INVARIANT-PF). Coarse map `a(x)=⟨v,X⟩`, `|v|=1` (the `k→0` block-uniform mode). Two copies `X,X′` under a synchronous (common-noise) or reflection coupling; coarse separation `δ_n=a(X_n)−a(X′_n)`. One step, `b(a):=E_μ[−⟨v,∇S⟩|a]` the **conditional mean force** (`b′<0` = restoring):
> `E[δ_{n+1}²]=δ_n²(1−2λΔt)+O(Δt)`,  **`λ=−b′(0)=G″_0`.**

`b′(0)=−G″` because `G″` is the second derivative of the pushed free energy `V_eff`: `G″=E[h_n|0]−Var(f|0)` (INVENT-VARIANCE). Pushed-forward coarse chain `ν=P_*μ`. **von Renesse–Sturm:** `W₂(νP_t,ν′P_t)≤e^{−λt}W₂` ⟺ `Hess V_eff⪰λ` — an *equivalence*, so **no coupling constant can beat `G″`**: the drift route *is* route (B) of OTHER-H1, and reflection/shift (Eberle, Lindvall–Rogers) couplings reduce to the same monotonicity `b′≤−λ`.

## 2. Acid test — mandatory
**Free massless.** The `k=0` mode is a **flat direction**: `⟨v,∇S⟩≡0`, so `b≡0`, `P` is a pure Brownian in `a` (Var grows `∝t`), synchronous coupling gives `Eδ²=const` ⇒ `λ=0`, `G″=0`. **No contraction — acid test PASSED.** Massive free: `b=−m²a`, `λ=m²`. **Pure Haar (β=0).** `S_W=0`, `V=Φ` convex, `Φ″≥N/12` ⇒ `b′(0)≤−N/12`, `λ≥N/12` — here the coupling **is a theorem** (`D≤0, S≥N/12`), exactly the PROOF-JACOBIAN-SAVES / -STOCHASTIC-LOCALIZATION calibration. All three exact ⇒ non-vacuous.

## 3. Does the Savvidy cone block the drift? — yes, exactly
`b′(0)=−G″=−(S−D)`; `D>0` is the Brascamp–Lieb defect of the conditioned fiber, non-zero **iff** `K=Π₀Hess VΠ₀` is indefinite — the anti-self-dual tachyon. On the coarse slice the **bare** Hessian `PᵀHP` has min ≈−0.3 (blocking is a low-pass filter; the Savvidy instability is *at* `k→0`, so the coarse slice is the worst one; MARGINAL). Thus the **pointwise** drift is *expanding* on the tachyon — no pointwise `(1−λ)`, `λ>0`. The averaged drift can be positive only through `S>D`. Crucially, at β=2.5 the Savvidy region is the **bulk** of the measure (equipartition `Var(G|0)≈E[H|0]≈2`; W7 action margin ≈0): a Doeblin "small set" cannot be the bad set, and the barrier the coupling must clear is `O(β)` *everywhere*. The measured marginal `G″=+0.51` (W7) is positive, but it is carried by the conditioning `E[H|0]=+0.22` — a **measure** fact, i.e. control of `Var(G|0)`, which is H1.

## 4. Two-route exhaustion (dynamical mirror of OTHER-H1)
**(A)** push the full LSI through `P`: `Var_ν(a)≤(1/λ_μ)E_μ|∇a|²` ⇒ needs the 2.5 gap — **circular**. **(B)** Bakry–Émery on `ν`: `Hess V_eff=G″⪰λ` — **H1**. The coupling *is* (B) in two-copy language; there is no third route.

## 5. Five-point scorecard
1. **Nonperturbative — PASS.** Pathwise drift, no loop expansion, exact transfer kernel.
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL@2.5.** Proven pointwise drift `≥0` only; `λ≥N/12` genuine at β=0.
3. **Non-circular — FAIL.** Drift rate `=G″`; von Renesse–Sturm is sharp, so proving `λ>0` *is* H1, zero slack.
4. **Non-SUSY + fails free field — PASS (cleanest).** Massless: driftless `k=0` Brownian, `λ=0`; massive `λ=m²`; β=0 `λ≥N/12`.
5. **Uniform lower bound — FAIL.** Doeblin minorization pays `e^{−cβ}`; `G″≥λ` uniform in `(L,a)` is H1.

## Exact failing estimate
```
∃λ>0: ∀x,x′,  E|a(X₁)−a(X′₁)| ≤ (1−λ)|a(x)−a(x′)|        (coarse W₂-contraction)
  ⟺  Hess V_eff ⪰ λ                          [von Renesse–Sturm, SHARP ⇒ no better λ]
  ⟺  G″(0)=S(v)−D(v) ≥ λ                      [drift b′(0)=−G″]
pointwise: bare PᵀHP|coarse min ≈ −0.3  ⇒ drift expanding on the tachyon;
  only measure-averaged S>D (D/S≈2%) survives ⇒ drift = G″ = H1.
full-chain Doeblin: minorization ε ~ e^{−cβ} ⇒ rate e^{−cβ}, non-uniform, ≪ true gap.
```

**Sharpest line.** The two-copy coupling is a real dynamical species: it never asks for static log-concavity and cleanly reproduces `λ=0` for the massless free field (a driftless `k=0` Brownian) and `λ≥N/12` for pure Haar. But von Renesse–Sturm makes "`W₂` contracts at rate `λ`" and "`Hess V_eff⪰λ`" the **same statement, sharply**, so the coupling is the dynamical face of H1 with no slack to exploit. The Savvidy cone blocks it exactly as it blocks Bakry–Émery, because the anti-self-dual tachyon *is* `D>0`, it sits at `k→0`, and at β=2.5 it is typical (equipartition) — no Doeblin small-set can push it into a rare event. The proven drift is `≥0`; the measured `+0.51` is H1 in coupling clothing. **#2,3,5 fail (one object: the drift = `G″`); #1,4 pass.**

*Refs.* Doeblin, Bull. Sci. Math. France **62** (1938); Meyn–Tweedie, *Markov Chains and Stochastic Stability* (2009); Hairer–Mattingly, Ann. Probab. **36** (2008) 1943; Hairer–Mattingly, *Yet another look at Harris' ergodic theorem* (2011); von Renesse–Sturm, CPAM **58** (2005) 923; Bakry–Émery 1985; Lindvall–Rogers, Ann. Probab. **14** (1986) 860; Eberle, PTRF **166** (2016) 851. Internal: PROOF-GAUSS-INVARIANT-PF (positive transfer kernel); PROOF-STOCHASTIC-LOCALIZATION (`Hess S_t=Hess S+tI`, covariance process); FORGE5-SL-MARTINGALE (QV identity, tight for Gaussian — the one-copy precedent); OTHER-H1-POINCARE (two-route exhaustion A/B, `Hess V_eff=G″`); PROOF-COARSE-POINCARE (`PᵀHP` indefinite at `k→0`); INVENT-VARIANCE-INEQ (`G″=S−D`, `D/S≈2%`); PROOF-JACOBIAN-SAVES (`Φ″≥N/12`, `G″=+0.51`); W6-BE-at-2.5 (min −1.12β, 90–97% fiber); MARGINAL (min −0.3); W7-Gpp0. No fabricated bound; this is a confirmed reduce-to-H1, derived fresh in the coupling vocabulary.

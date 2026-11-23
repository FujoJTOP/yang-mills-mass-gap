# GAP-IN12 — can a variational upper bound force per-channel sharpness?

> 2026-11-23. Attacks the matching bound `Δ_X ≤ K_X·cert_X` (the per-channel form of `Δ ≤ K·cert`, `X ∈ {0⁺⁺, 2⁺⁺}`) through a min–max trial state, completing `BREAK-IN-10`. Verdict: the variational route gives the correct, **non-circular architecture** — but closing it is *strictly stronger* than the gap, so it neither follows from IN-0…IN-11 nor bypasses the wall.

## (1) The sandwich is sound, not circular

Take a trial state `ψ_X = O_X|Ω⟩` on `Ω⊥`, with `O_X = tr F²` (0⁺⁺) or `T_{μν}` (2⁺⁺). Its spectral measure `ρ_ψ` has `inf supp ρ_ψ ≥ Δ_X`, so the Rayleigh quotient is an upper bound on the infimum:
```
Δ_X = inf supp ρ_X  ≤  (∫s dρ_ψ)/(∫dρ_ψ) =: R(ψ_X) .
```
Certificate validity gives `cert_X ≤ Δ_X`. Hence `cert_X ≤ Δ_X ≤ R(ψ_X)`, **both endpoints computable from the axioms** — no knowledge of `Δ_X` needed. To get `Δ_X ≤ K_X cert_X` it suffices to exhibit one explicit `ψ_X` with

> **(V)  R(ψ_X) ≤ K_X·cert_X .**

If (V) holds, `cert_X ≤ Δ_X ≤ R(ψ_X) ≤ K_X cert_X`. This is the *right* shape and proves the two-sided bound without solving the gap — the route is logically valid.

## (2) Why (V) does not follow, and is gap-strength

Combining (V) with validity forces `Δ_X ≤ R(ψ_X) ≤ K_X Δ_X`, i.e. `ψ_X` is a **near-optimal trial glueball** (factor-`K` window, not the exact mass). That is strictly weaker than *knowing* `Δ_X` but strictly stronger than the gap (below). Concretely it is unavailable:

- `R(ψ_X)` is the zero-momentum ratio of connected 4- to 2-point functions of `tr F²`/`T_{μν}` — UV-divergent, renormalizable only along the RT (consumes IN-0/IN-9 data), and computable to the needed precision only by controlling the IR spectral weight = the gap.
- `cert_X` is a feasibility threshold of a *different* constraint set (SDP / RG scale); nothing forces that floor within O(1) of an attainable trial energy — loose bootstrap bounds are generic. Min–max supplies the upper endpoint `R`, never the comparison `R ≤ K cert`.

## (3) Sharpness is gap-strength and then some

`cert_X ≤ Δ_X ≤ K_X cert_X` with `cert_X>0` **implies** `Δ_X>0`; the converse fails (BREAK-IN-10). So IN-12 (`hsharp`):

- is **not derivable** from IN-0…IN-11;
- is **not a way around** the wall — any proof must already prove the mass gap, plus construct a near-optimal trial state;
- is at best reachable at finite `(a,L)` (`cert_a ≤ Δ_a ≤ R_a(ψ)`), but the sandwich neither runs at the crossover coupling (IN-3/IN-4) nor survives `a→0` without IN-9.

**Cleanest outcome.** "cert ≡ gap" stays false as a theorem. Min–max upgrades it at most to
`cert ≍ gap ⟺ (gap ∧ ∃ψ_X: R(ψ_X) ≤ K_X cert_X)`,
a conditional whose second conjunct is a constructive glueball-mass bound **strictly stronger than the gap**. IN-12 remains a fresh axiom; the missing lemma is a renormalized trial-state energy bound `R(ψ_X) ≤ K_X cert_X`, which is the honest content of sharpness and is currently unavailable in 4D.

# BREAK-IN-10 — can the matching upper bound `Δ/Λ ≤ K·cert` be proved?

> 2026-11-23. Attempts the reverse direction that would restore `cert ≍ gap`. Verdict: **not provable from IN-0…IN-11**; the upper bound is *logically independent* of the current inputs and *strictly stronger* — it is exactly per-channel **sharpness** of the certificate.

## (1) Exact reduction

Let `Δ_X` be the channel gap, `cert μ = min_X cert_X(μ)`, `X ∈ {0⁺⁺,2⁺⁺}`. IN-11 gives `Δ = min_X Δ_X`. The one-sided `cert ≤ Δ` uses only `cert_X ≤ Δ_X`. The wanted `Δ ≤ K·cert` is, with `X*` the argmin, `Δ_{X*} ≤ K·cert_{X*}`:
> **the gap-carrying channel's certificate must be within a constant `K` of the truth**, i.e. `cert_X ≥ Δ_X/K_X` for each visible `X`.

This is a *new object* — optimality-up-to-constant (sharpness) of the bootstrap bound — not its validity.

## (2) The variational route cannot close it

Min–max gives only `Δ_X ≤ R(ψ)` for a trial state (an upper bound on the *truth*). Combined with `cert_X ≤ Δ_X ≤ R(ψ)` this **brackets** `Δ_X ∈ [cert_X, R(ψ)]`. Collapsing to `Δ_X ≤ K·cert_X` requires `R(ψ) ≤ K·cert_X`, i.e. `cert_X ≈ Δ_X` — sharpness again. The variational principle controls the bracket's *top*; it cannot make a valid lower bound tight. It yields a *bracketing* certificate `[cert_low, cert_high]`, not a one-sided bound against `cert` alone.

## (3) Exact failure (separation)

IN-3/IN-4 supply only positivity: `cert ≥ δ_W` or `δ_S > 0`, with **no relation between `δ` and `Δ`**. Every stated inequality (`δ ≤ cert ≤ Δ`, channel lower bounds, IN-11) holds for each `Δ/δ ∈ [1,∞)`; letting `δ/Δ → 0` makes `Δ/cert` unbounded, so no admissible `K` exists. The upper bound is therefore not implied — it is a fresh quantitative axiom (bootstrap completeness/sharpness). It is also not realizable now: the 4D DB is unconverged (only a 1D GFF toy verified end-to-end), so `cert_X ≪ Δ_X` and `Δ ≤ K·cert` fails at any small `K`. (At deep UV `cert ~ 1/a ≫ Δ`, so the upper bound is UV-trivial; its entire content lives at the IR/crossover wall — it adds a requirement there, it does not bypass it.)

## (4) Effect on the equivalence claim

`cert>0 ⇔ Δ>0` is exactly `(IN-10 ∧ hsharp)`: forward uses `cert ≤ Δ`, backward uses `hsharp: cert ≥ Δ/K`. Adding `hsharp` makes the biconditional *true but conditional* — no net gain over one-sided sufficiency. "cert ≡ gap" stays false as a theorem; it is restorable only by *axiom*.

## (5) Lean-level change

In `Handoff.lean`, `mass_gap_of_certificate` (line 77) consumes `hcert : cert μ ≤ Δ_over_Λ`. Add
```
(hsharp : ∀ μ, 0<μ → Δ_over_Λ ≤ K * cert μ)
```
and prove `∀ μ, 0<μ → (0 < cert μ ↔ 0 < Δ_over_Λ)` (→ via `hcert`, ← via `hsharp`). Label it **IN-12 (`hsharp`)** — a *new axiom*, not derived.

# IN-10 (`hcert`) — attack record

> 2026-11-23. Companion to `W4-cert-vs-gap.md`.

## (1) Precise statement

`hcert : ∀ μ>0, cert μ ≤ Δ/Λ` (`Handoff.lean:77`). `Δ/Λ` is **one fixed** physical number; `cert μ` is the scale-`μ` certificate, defined as `min(cert_{0⁺⁺}(μ), cert_{2⁺⁺}(μ))` — the min of the two bulk-channel lower bounds (§3bis). So IN-10 asserts: *the certificate never exceeds the true gap*. It is strictly **one-sided** — no lower bound on `Δ/Λ`, no tightness, no `cert>0` (that is IN-3/IN-4). It **cannot be tightened to two-sided by logic**: a two-sided form needs a genuine *upper*-bound certificate `Δ/Λ ≤ K·cert μ`, i.e. tightness of `cert` up to a constant — a new object, not a consequence of the existing one.

## (2) True status

Explicit assumption, not proved. Docs name it the fifth assumption / `M6` (`MINIMAL.md` §6); the Lean theorem consumes it as a hypothesis. It is **nonvacuous**: `MINIMAL.md` §8.3 shows that at deep UV `cert ~ 1/a ≫ Δ/Λ`, violating it — so it is meaningful only at the IR/crossover scale.

## (3) Can it be advanced?

Partially, yes. If each channel certificate is a **true lower bound**, `cert_X(μ) ≤ Δ_X`, then
`cert(μ) = min_X cert_X ≤ min(Δ_{0⁺⁺}, Δ_{2⁺⁺}) = Δ` **by IN-11**.
So IN-10 is **not an independent fifth axiom**: it is *implied* by IN-11 plus the lower-bound direction of the two channel certificates. The residual content is exactly that direction (`cert_X ≤ Δ_X`), which is the supplier obligation of IN-3/IN-4/DB — plus the §8.3 scale restriction. It is neither pure definition nor irreducible bridge; it is a corollary of (IN-11 ∧ "certificates are lower bounds").

## (4) What W4-3 changes

W4-3 kills `gap ⇒ cert>0` (necessity), **not** `cert ≤ Δ` (the only direction Lean uses). Hence the conditional theorem's **proof is untouched** — it stays valid *sufficiency*: `cert>0 ∧ IN-10 ∧ IN-11 ⇒ gap`. What is retracted is the *meta*-claim "cert ≡ gap": certificate-free routes (RP/transfer-matrix, measure-geometric) are not excluded. The theorem is "inputs ⇒ gap," never "gap iff inputs."

## (5) Flip conditions

- **Equivalence returns** if one proves an upper-bound certificate `Δ/Λ ≤ K·cert μ` at a physical scale (cert becomes a surrogate for the gap).
- **IN-10 fails** if a bulk state below `δ` appears outside `0⁺⁺/2⁺⁺` (falsifies IN-11) or if `cert>0` at `Δ=0`.
- **The `∀μ` form fails as written** at deep UV; it must be read at IR/crossover (§8.3).

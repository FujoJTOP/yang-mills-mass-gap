# G1-LARGE-N — does large N give a controlled (G1) ("pure YM has no IR fixed point")?

**Result: `[T]` 2-loop β<0 for pure YM at all λ, but `[NG]` large-N does NOT give a controlled nonperturbative (G1) — it reduces (G1) to large-N confinement, the same strong-coupling wall.**

## The 2-loop 't Hooft β (large-N leading) `[T]`

With `λ=g²N`, `N_f=x·N`, the N factors cancel and

```
β_λ(λ,x) = −(2/3)[ (11−2x) λ²/(16π²) + (34−12x) λ³/(16π²)² ].
```

- **Pure YM (`x=0`)**: `(11−2x)=11>0`, `(34−12x)=34>0` ⇒ **both terms negative ⇒ `β_λ(λ)<0` for every `λ>0`** — **no IR fixed point** (witness: β_λ = −4.7e−4 … −1.4e+3 across λ=0.1…100, no zero).
- **Conformal window**: an IR fixed point `λ* = 16π²(11−2x)/(12x−34) > 0` requires **`34/12 < x < 11/2`**, i.e. `N_f/N ∈ (2.833, 5.5)` (upper edge `11/2` = Banks–Zaks `β₀=0`; lower edge `34/12` = `β₁=0`). Witness: `x=4 ⇒ λ*=33.84` (β crosses zero), `x=0 ⇒ none`.
- **Window-boundary scaling**: in the 't Hooft variable the window edges are **`O(1)` constants** (`11/2`, `34/12`), so pure YM (`N_f/N=0`) sits **far below** the lower edge — "outside the window" is **robust** at large N (not marginal).

## Why this is NOT a controlled (G1)

1. **Perturbative only `[C]/[O]`**: the β beyond 2 loops is **scheme-dependent**, and the fixed-point *existence* is a nonperturbative question. The 2-loop "no zero at x=0" does not exclude a strong-coupling IR fixed point (see `STAR-NO-BZ`, `STAR-COUNTEREXAMPLE`: 2-loop β<0 ≠ no nonperturbative zero).
2. **The small parameter at large N is `1/N`, for FLUCTUATIONS — not for the fixed-point question.** Large N controls the master field / meson spectrum, not the existence of a fixed point; the coupling `λ` is `O(1)` at the confining scale, exactly where (G1) lives.
3. **The 4D master field is unconstructed** (Singer 2D only) — large-N YM's "solution" in 4D is itself open.
4. **"Pure YM has no IR fixed point at N_f=0" ≡ large-N confinement** — a rigorous nonperturbative statement does not exist (lattice strongly supports it; no proof).

## Verdict

`[T]` β_λ<0 ∀λ at x=0 and the window `(34/12, 11/2)` (2-loop). But **large N does NOT give a controlled (G1)**: the window boundary is an `O(1)` statement only perturbatively, the nonperturbative "no fixed point at N_f=0" is large-N confinement (unproven), and `1/N` controls fluctuations, not the fixed point. (G1) at large N **reduces to the same strong-coupling/confinement problem = `H1` in the large-N coordinate**. Clay 0%; no QED.

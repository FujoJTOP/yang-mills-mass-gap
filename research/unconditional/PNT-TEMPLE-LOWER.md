# PNT-TEMPLE-LOWER — Temple/Weinstein do NOT supply the missing lower bound on `E₁`: the bound is dominated by the Rayleigh quotient and collapses to `E₀`

> **VERDICT (no QED).** Temple's and Weinstein's bounds give *valid*, renormalon-free lower bounds
> on the first excited energy `E₁` from variational data (`ρ = ⟨H⟩`, variance `σ²`, an upper bound `W`
> on `E₁`). But **Temple's bound is rigorously `≤ ρ`** — dominated by the Rayleigh quotient — and `ρ`
> lower-bounds `E₁` only under the *unverifiable* hypothesis `ρ < E₁`. For **every tractable ground-state
> trial** (`ρ ≈ E₀`) the delivered gap lower bound sits **at or below `E₀`** (exactly `E₀` for a two-level
> mixture; numerically `≤ E₀` over 20 000 random trials). Hence Temple/Weinstein **do not close the
> tensor-RG sandwich**: the missing lower half is *not* a Temple certificate, it is the same
> nonperturbative input — a near-exact `0⁺⁺` wavefunctional with controlled variance, or an independent
> spectral-separation certificate. **Clay 0%.**

---

## §0 · The bounds and their validity conditions

- **Temple (1928):** `E₁ ≥ ρ − σ²/(W − ρ)`, valid when `ρ < E₁` (and `ρ < W`), with `W ≥ E₁` an
  *upper* bound on `E₁`, `ρ = ⟨H⟩`, `σ² = ⟨H²⟩ − ρ²` on a normalized trial.
- **Weinstein (1934):** `E₁ ≥ ρ − σ`, valid when `ρ < E₂`.
- Both are **renormalon-free**: only Rayleigh quotients and variances (finite expectations) enter — no
  `u=2` Borel pole. First requirement (a finite, scheme-independent lower-bound *type*) is met `[T]`.

## §1 · Rigorous dominance (the decisive fact)

`σ² ≥ 0` and (validity) `W > ρ` give, identically,

> **Temple bound `= ρ − σ²/(W − ρ) ≤ ρ`.**  `[T]`

So Temple never beats the **Rayleigh quotient** `ρ`. But `ρ` lower-bounds `E₁` *only if* `ρ < E₁` — which
is the unknown the bound is supposed to certify. Either (a) `ρ < E₁` is known, then `ρ` itself is the
better bound and Temple is redundant; or (b) `ρ < E₁` is unknown, and Temple certifies nothing beyond
`ρ − σ²/(W−ρ) ≤ ρ` using the checkable `W`. **In no case does Temple add a gap lower bound beyond `ρ`.**

## §2 · Calibration (exact solvable models, `w7/pnt_temple_lower.py`)

- **Harmonic oscillator** (`E₀=0.5, E₁=1.5`): Gaussian trials give valid bounds **≤ ρ**, peaking at
  `0.5 = E₀` at the optimal (exact-ground) width; never above `E₀`.
- **Random symmetric `8×8` matrix** (exact gap `Δ=0.482`): trial `√(1−ε²)v₀ + ε v₁`. **Two-level law:**
  `ρ = E₀ + ε²Δ`, `σ² = ε²(1−ε²)Δ²`, so `σ²/(W−ρ) = ε²Δ` and **Temple `= E₀` exactly for all `ε`**
  (0 % of the gap). Random near-ground scan (20 000 trials): `max(bound − ρ) = 0` (dominance exact);
  `max bound = −3.44 < E₀ = −3.21` — **never exceeds `E₀`**.

## §3 · 3+1D assessment — the exact missing input

To turn Temple into a positive *gap* certificate one needs `σ² ≪ (W−ρ)(ρ−E₀)`, i.e. a trial with
**variance much smaller than `gap²`** and `ρ < E₁` certified. In 3+1D nonabelian gauge theory that trial is
the **0⁺⁺ vacuum wavefunctional with controlled variance** (plus a tight `W` and a gap-separation
certificate). Constructing it nonperturbatively is the theory itself; computing `σ² = ⟨H²⟩` needs a 4D
nonperturbative expectation. The only escape from the dominance (`ρ) is an **independent spectral-
separation input** — which is the same wall.

**Labels:** valid-bound-type `[T]`; dominance `≤ ρ` `[T]`; collapse-to-`E₀` law `[T]`/`[V-num]`;
3+1D gap certificate `[NG]` for the Temple class; exact missing input `[O]` = controlled-variance `0⁺⁺`
wavefunctional / spectral separation. **No QED.**

*Witness:* `research/unconditional/w7/pnt_temple_lower.py` (numpy only, `<1 s`).

# INEQ-R3-hylleraas — two-state / Hylleraas–Undheim on the TRANSITION amplitude

> 2026-11-23 · Tests whether the *off-diagonal* residue `f₁²=|⟨Ω|O|1⟩|²` (RESIDUE-CONDITIONAL's `H*(r₀)`) is more variationally accessible than the diagonal edge bounds of LEVER-nonconvex-edge. This is the **Hilbert-space / Ritz** route, complementary to the measure-level dual of INEQ-T2-variational-residue. **Verdict: not more accessible.** A transition amplitude is *stationary, not extremal*; every certified *lower* bound (Eckart–Löwdin–Shull, Miller, Temple) consumes a **lower** bound on the next exact level `E_{k+1}` — the same missing spectral datum, in a direction Ritz cannot supply. Naming: HUM = Hylleraas–Undheim–**MacDonald** (interlacing, *upper*); **Miller** (1969) = the stationary *lower*-bound method, which needs `H²` moments.

## (0) Direct inequalities run the wrong way
Rayleigh–Ritz / HUM on nested trial subspaces gives `ε_k ≥ E_k` for *every* level (interlacing, min–max): upper bounds. There is no extremum principle for a single matrix element. Cauchy–Schwarz yields `|t|² ≤ ⟨Ω|O†O|Ω⟩·⟨1|O†O|1⟩` — *upper*. The sum rule `Σₙ|tₙ|²=⟨Ω|O²|Ω⟩` is a fixed linear datum, blind to a lone edge atom (linear-data lemma, LEVER §0B). So the readily available inequalities bound `|t|` from **above**.

## (1) What the two-state / HUM construction delivers
In the trial space `{φ₀,φ₁}≈{Ω,ψ₁}`, diagonalising `H` gives Ritz values `ε₀≤ε₁` with `ε₁≥E₁` (HUM, *upper*). The projected `⟨ε₀|O|ε₁⟩` is an **estimate** of `t`: because `t` is *stationary* under state variation, its error is `O(δ²)` when the trial vectors are first-order accurate — accuracy without a fixed sign. Concretely `⟨φ₀|O|φ₁⟩ = c₀c₁ t + (leakage)`, no one-sided inequality.

## (2) Turning the estimate into a bound needs `E_{k+1}`
To get `|t|≥…` one needs `|c_k|` bounded **below**. The Eckart–Löwdin–Shull overlap-accuracy inequality `|c_k|² ≥ 1−(ε_k−E_k)/(E_{k+1}−E_k)` requires a **lower bound on the next exact level `E_{k+1}`** — exactly what Ritz (upper) cannot give and the campaign cannot produce (one-rung regress, LEVER §1). Miller's stationary lower-bound method needs `H²` matrix elements (second spectral moments); Temple lower-bounds `E₀` from an upper `E₁` but lower-bounds `E₁` only from an upper `E₂`. The residue–level-gap lemma `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)` has this same shape: a lower bound on `f₁²` *given* the second threshold `Λ=E₂`.

## (3) The dual factorisation is a re-expression, not a gift
If `O=[H,A]`, hypervirial gives `t=(E₁−E₀)⟨Ω|A|1⟩`; Feynman–Hellmann-style eigenvector differentiation gives `⟨Ω|∂_λ|1⟩ = ⟨Ω|∂_λH|1⟩/(E₁−E₀)`. Both move the unknown into a *gap denominator* times another unprotected eigenvector element — no cheaper object.

**Constructed statement (honest).** A two-state/HUM scheme yields a *consistency pair* `ε₁≥E₁` plus `|t|≈|⟨ε₀|O|ε₁⟩|` to `O(δ²)`. A certified lower bound on `|t|` requires, as input, a **lower bound on the gap / next level** — i.e. `H*(r₀)`'s spectral edge. The transition amplitude is thus *less* accessible to certified bounds than the eigenvalue, which at least possesses min–max. The off-diagonal gives no new lever; `f₁²≥r₀` remains the irreducible unprotected datum.

*Refs.* INEQ-T2-variational-residue (measure-level dual, infimum/usc); EDGE-f1sq-attack (sign obstruction, §3 residue≡`Λ`); RESIDUE-CONDITIONAL (residue–level-gap lemma, `H*`); LEVER-nonconvex-edge §0B,§1; GAP-G2-edge; Eckart (1930), Shull–Löwdin, Cohen–Feldmann–McEachran (1971), Miller, J. Chem. Phys. 50, 2758 (1969); Temple–Weinstein lower bounds.

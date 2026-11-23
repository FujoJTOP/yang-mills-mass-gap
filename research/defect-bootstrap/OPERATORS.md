# Pure-gauge `Z_N`: bulk–defect operator content and crossing

> CPB-attack probe (sub-task 8). Verdict: the operator content is clean and `Z_N`-structured;
> the spectrum `{q_i, Δ_i}` is what the bootstrap must supply.

## Bulk

All local operators are **`Z_N`-neutral** (the 1-form symmetry acts trivially on local
operators): `tr F²` (`Δ=4`), `tr F F̃` (`Δ=4`), `tr F⁴` (`Δ=8`), … — all `q=0`.

## Defect (on the Wilson line `W_R`, charge `q` = N-ality of `R`)

Defect operators `O_i(τ)` are **graded by `q_i ∈ Z_N`** (`H_def = ⊕_{q_i} H^{(q_i)}`).
Lightest: displacement `D^a` (adjoint, `Δ=2`), then transverse-gluon / field-strength
primaries `Δ ≥ 2`. The `q_i` of each is dynamical defect CFT data.

## OPE coefficients (allowed by `Z_N`)

- Bulk–defect: `a_i ≠ 0` only when `q_bulk + q_i = q_line` ⟹ **`q_i = q_line`** (bulk `q=0`).
- Defect–defect: charge adds mod `N`; the `a_i²` are `Z_N`-graded.

## Crossing (neutral bulk `φ`)

```
⟨φφ⟩_def(q) = Σ_O C^{φφ}_O 𝒢_O(ξ) = Σ_{q_i} a_i² 𝔽_i(ξ) ,
```

the defect sum over `Z_N`-graded sectors with relative phases fixed by the anomaly 2-cocycle
`ω(q_i,q_j) = e^{2πi k q_i q_j/N}` (`k` = level, θ-related), an **irremovable projective rep**.
The bulk channel has no such phase.

## Honest boundary

The exact phase/level (`k`, θ-normalisation) is BF-convention-dependent; the `q_i` of `D` and
the defect spectrum are data the bootstrap must supply; non-SUSY pure-gauge defect CFT data
are unbuilt.

Sources: GKSW 1412.5148; Cordova–Tachikawa; 2501.06900; 1601.02883.

# Casimir equation for the defect block — codim-2 caveat

> CPB-attack probe (sub-task 3, core). Produced the **codim-2** (1-variable) hypergeometric
> block; whether this is the 4D-line block is **disputed** (resolve via arXiv:2501.06900).

## Setup (as briefed: `SO(2,1)×SO(2)`)

Defect along `τ=x³`; group `SO(2,1)_τ × SO(2)_⊥`; a defect primary `O_j` has dimension
`Δ_j` and transverse spin `s_j`. Two bulk scalars of dimension `Δ_φ` at `(τ, y)`. The
bulk–defect block `𝒢_j(z)` is a function of **one** cross-ratio `z`.

## Casimir

With the `SO(2,1)` Casimir `L₀² − ½{L₁,L₋₁}`, the eigenvalue on the exchanged defect primary
is

```
C_j = Δ_j(Δ_j−1) + s_j² .
```

The ODE:

```
z(1−z) 𝒢'' + [c_j − (a_j+b_j+1) z] 𝒢' − a_j b_j 𝒢 = 0 ,
a_j = b_j = (Δ_j+s_j)/2 ,   c_j = Δ_j+s_j .
```

## Solution

```
𝒢_j(z) = z^a ₂F₁(a, a; 2a; z) ,   a = (Δ_j+s_j)/2 .
```

## Boundary conditions

(i) OPE `z→0`: `𝒢_j ~ z^a`, subleading `z^{a+1}`. (ii) `z→1`: branch point; the crossing
kernel carries the anomalous phase `e^{iθ/N}` fixed by the `Z_N^{(1)}`–θ anomaly. (iii)
single-valuedness / positivity on the physical sheet fixes the contour.

## Honest boundary

The 1D global block `z^Δ ₂F₁(Δ,Δ;2Δ;z)` is known; the defect block **with transverse spin** is
the object being constructed. The codimension-dependent shift (full `SO(3)→SO(2)` reduction,
non-SUSY) and the closed form are **not in the literature**. Established here: the eigenvalue
`C_j = Δ_j(Δ_j−1)+s_j²` and the hypergeometric structure; the precise constants need the full
derivation.

## Key caveat

This is the **codim-2** (`SO(2)`, 1 cross-ratio) block. If the 4D line is **codim-3**
(`SO(3)`, 3 cross-ratios — per `KINEMATICS.md`), this is the **wrong object**. Resolution
pending (arXiv:2501.06900).

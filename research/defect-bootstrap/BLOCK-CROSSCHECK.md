# Block cross-check — independent derivation of the codim-3 bulk–defect two-point block

> Independent route (defect-OPE → 1D conformal block + a direct differential-equation solve in the
> H² cross-ratio `u`), **not** by copying the `₂F₁` guess of `BLOCK-CORRECTED.md`.
> Result: the physical block is a *single* `₂F₁` with a `(1−u)^{Δ_j}` prefactor. The two old
> candidates `f₁=P_{Δ_j−1}(cosh η)` and `f₂=Q_{Δ_j−1}(cosh η)` are **harmonic (H²-Laplacian) partial
> waves**, not conformal blocks — that is why neither, and no linear combination of them, can hold
> both OPE limits. Both requirements, the exact ODE and numerics are below.

---

## 0. Result

Physical reduced block (normalised `f(0)=1`):

```
                     Δ_j                                ___
  f(u) = (1 − u)        ·  ₂F₁( Δ_j − Δ_φ , Δ_j − Δ_φ ; 2Δ_j ; 1 − u )   /   N ,

  N = ₂F₁( Δ_j − Δ_φ , Δ_j − Δ_φ ; 2Δ_j ; 1 )  =  Γ(2Δ_j) Γ(2Δ_φ) / Γ(Δ_j + Δ_φ)²
```

with `u = (cosh η − 1)/(cosh η + 1) = d²/(4ρ₁ρ₂ + d²)`, `d² = z₁₂² + (ρ₁−ρ₂)²` (θ=0 slice).
The full block is `f(u)·P_j(cos θ)`.

**Governing equation (defect-channel conformal Casimir, eigenvalue `Δ_j(Δ_j−1)`):**

```
  u (u−1)² f''  +  (1 − 2Δ_φ) (u−1)² f'  +  Δ_φ² (u−1) f   =   Δ_j(Δ_j−1) f .        (⋆)
```

Indicial exponents of (⋆):

| point | exponents | physics |
|---|---|---|
| `u = 1` (defect OPE, a point on the line) | `{Δ_j , 1 − Δ_j}` | exchanged primary and its **shadow** |
| `u = 0` (bulk OPE, `x₁→x₂`) | `{0 , 2Δ_φ}` | identity and the first bulk operator |

So the two OPE endpoints are `(1−u)^{Δ_j}` at `u→1` and `f→1` at `u→0` **by construction** — they are
the two physical indicial branches, and `f` is the unique branch regular at `u=0` normalised to 1.

## 1. Why `f₁` and `f₂` fail (the real reason)

`f₁ = ₂F₁(Δ_j,1−Δ_j;1;−u/(1−u)) = P_{Δ_j−1}(cosh η)` and `f₂ = Q_{Δ_j−1}(cosh η)` are the two
solutions of the **H² Laplacian ("harmonic") equation**

```
  u g'' + g' − Δ_j(Δ_j−1)/(1−u)² g = 0     (the free / mass-Δ_j KG eq.; = the Legendre eq. in cosh η)
```

i.e. the *free* bulk-to-bulk propagator in the transverse H² (cf. `FREE-CODIM3.md`). That is **not**
the defect Casimir (⋆). Numerical confirmation (θ=0, Δ_j=2, defect-OPE wants `~(1−u)^2`):

| u | `f₁ = P₁(coshη)` | `f₂ = Q₁(coshη)` | `(1−u)^2` |
|---|---|---|---|
| 0.02 | 1.04 | 1.04 | 9.6e−1 |
| 0.20 | 1.50 | 2.07e−1 | 6.4e−1 |
| 0.80 | 9.00 | 4.15e−3 | 4.0e−2 |
| 0.95 | 39.0 | 2.19e−4 | 2.5e−3 |
| 0.99 | 199  | 8.4e−6 | 1.0e−4 |

`f₁` grows (wrong defect OPE); `f₂` decays correctly at `u→1` but **log-diverges** at `u→0`
(`Q₀ ∼ −½ln u` etc.), so it fails the bulk OPE. Because they are the two harmonic solutions, **no
linear combination** `A f₁ + B f₂` fixes both ends (finiteness at `u→0` forces `B=0`). The fix is not
a different combination — it is a **different equation**: the conformal block (⋆).

## 2. Independent derivation

**(a) Defect OPE reduces the radial part to a 1D conformal block.** For a codim-3 line the symmetry is
`SO(2,1)×SO(3)`. The `SO(3)` part gives `P_j(cos θ)`; the `SO(2,1)` part acts on the H² point
`w = z + iρ` (Möbius; verified: the special-conformal `(z,ρ)` map is `w ↦ w/(1−bw)`). A bulk scalar of
dimension `Δ_φ` is, under `SO(2,1)`, a weight-`Δ_φ` field on H²; its defect-OPE partner (the defect
identity) has dimension `0`. So the radial block is the standard **1D conformal four-point block** with
external dimensions `(Δ_φ, 0, 0, Δ_φ)`, internal `Δ_j`, in the cross-ratio `z = 1 − u`:

```
  g(z) = z^{Δ_j} ₂F₁( Δ_j + Δ₂ − Δ₁ , Δ_j + Δ₃ − Δ₄ ; 2Δ_j ; z ) ,   (Δ₁,Δ₄)=Δ_φ , (Δ₂,Δ₃)=0
       = (1−u)^{Δ_j} ₂F₁( Δ_j − Δ_φ , Δ_j − Δ_φ ; 2Δ_j ; 1 − u ) .
```

The variable `z = 1−u` is forced by the limits: `z→0` must be the defect OPE (`u→1`) and `z→1` the bulk
coincidence (`u→0`).

**(b) Equivalent ODE.** Eliminating `₂F₁` gives exactly (⋆). Its indicial exponents (computed with sympy)
are `{Δ_j,1−Δ_j}` at `u=1` and `{0,2Δ_φ}` at `u=0` — the physical OPE data, an independent check that
(⋆) is the defect Casimir and not the harmonic equation.

**(c) Both limits fixed.** `f(0)=1` (normalise by `N`); `f(u→1) → (1−u)^{Δ_j}·₂F₁(Δ_j−Δ_φ,Δ_j−Δ_φ;2Δ_j;0)/N = (1−u)^{Δ_j}/N`.
The bulk-OPE correction is `O(u^{2Δ_φ})` (double-trace twist), and the defect-OPE tower is the
integer-spaced `(1−u)^{Δ_j+n}` — both inherited from `₂F₁`.

## 3. The `Δ_φ` normalisation

The reduced block above is `Δ_φ`-normalised by stripping the external factor. The full two-point block is

```
  ⟨φ(x₁) φ(x₂)⟩_{O_j}  =  a_j · (4 ρ₁ ρ₂ u)^{−Δ_φ} · (1 − u)^{Δ_j}
                               · ₂F₁( Δ_j−Δ_φ , Δ_j−Δ_φ ; 2Δ_j ; 1−u ) · P_j(cos θ) / N .
```

Equivalently, with the standard defect two-point form `⟨φφ⟩ = (ρ₁ρ₂)^{−Δ_φ} 𝔽(u)` one has
`𝔽(u) = (4u)^{−Δ_φ} f(u)`. `Δ_φ` enters in **two** places only: (i) the external factor `(4ρ₁ρ₂u)^{−Δ_φ}`
(= `ρ^{−Δ_φ}` per point × `(4u)^{−Δ_φ}`), which produces the bulk-OPE singularity
`⟨φφ⟩ ~ (4ρ²u)^{−Δ_φ}` as `u→0`; and (ii) the `₂F₁` parameters `Δ_j−Δ_φ`. With this prefactor the
defect-OPE limit is `⟨φφ⟩ ~ (1−u)^{Δ_j−Δ_φ}` (since `ρ₁ ∝ 1−u`), as it must be.

## 4. Integer / half-integer `Δ_j` (log case)

At `u=1` the exponents `{Δ_j, 1−Δ_j}` differ by `2Δ_j−1`. When `Δ_j ∈ ½ℤ` (in particular all integer
`Δ_j ≥ 1`) this is an integer, so the **second (shadow) solution** of (⋆) carries a `(1−u)^{1−Δ_j}·log(1−u)`
term. The **physical block** `f` is the `Δ_j` branch — the analytic series
`(1−u)^{Δ_j}(1 + c₁(1−u) + …)`, which has **no logarithm** for any `Δ_j>0` (verified numerically at
`Δ_j=1,2,3` below: `f/(1−u)^{Δ_j}` converges smoothly). The log lives entirely in the discarded shadow
branch (and in the harmonic `Q_{Δ_j−1}`). So there is no separate `f`-log to handle; `f` is the
`(1−u)^{Δ_j}` branch of (⋆), and the second solution of (⋆) is the shadow `(1−u)^{1−Δ_j}` branch, which
acquires a `log(1−u)` (standard integer-exponent hypergeometric behaviour) precisely when `2Δ_j−1∈ℤ`.

## 5. Numerical verification

Closed form evaluated with `scipy.special.hyp2f1` (Python 3.12, scipy 1.18.1). `Δ_φ = 2` (representative;
`Δ_φ` only shifts the `₂F₁` parameters and the prefactor). `f` normalised so `f(0)=1`.
Two OPE tests: **bulk** `f(u→0)→1`; **defect** `f/(1−u)^{Δ_j} → const` as `u→1`.

| `Δ_j` | `u` | `f(u)` | `f/(1−u)^{Δ_j}` |
|---|---|---|---|
| 1.0 | 0.05 | 9.341667e−01 | 9.833333e−01 |
| 1.0 | 0.20 | 7.466667e−01 | 9.333333e−01 |
| 1.0 | 0.50 | 4.166667e−01 | 8.333333e−01 |
| 1.0 | 0.80 | 1.466667e−01 | 7.333333e−01 |
| 1.0 | 0.95 | 3.416667e−02 | 6.833333e−01 |
| 1.5 | 0.05 | 9.220990e−01 | 9.958459e−01 |
| 1.5 | 0.20 | 7.037503e−01 | 9.835210e−01 |
| 1.5 | 0.50 | 3.391953e−01 | 9.593893e−01 |
| 1.5 | 0.80 | 8.370293e−02 | 9.358272e−01 |
| 1.5 | 0.95 | 1.033320e−02 | 9.242295e−01 |
| 2.0 | 0.05 | 9.025000e−01 | 1.000000e+00 |
| 2.0 | 0.20 | 6.400000e−01 | 1.000000e+00 |
| 2.0 | 0.50 | 2.500000e−01 | 1.000000e+00 |
| 2.0 | 0.80 | 4.000000e−02 | 1.000000e+00 |
| 2.0 | 0.95 | 2.500000e−03 | 1.000000e+00 |
| 3.0 | 0.05 | 8.437233e−01 | 9.840774e−01 |
| 3.0 | 0.20 | 4.827611e−01 | 9.428927e−01 |
| 3.0 | 0.50 | 1.098138e−01 | 8.785108e−01 |
| 3.0 | 0.80 | 6.626529e−03 | 8.283161e−01 |
| 3.0 | 0.95 | 1.008455e−04 | 8.067637e−01 |

**Bulk-OPE column** `f → 1`: e.g. `f(1e−5) = 0.999987, 0.999984, 0.999980, 0.999967` at
`Δ_j = 1, 1.5, 2, 3`. **Defect-OPE column**: `f/(1−u)^{Δ_j}` is flat/converging as `u→1`
(`→ 0.6668, 0.9204, 1.0000, 0.8000` at `u=0.9997`). Both limits hold for every `Δ_j`.

Cross-checks (non-circular):
- **ODE residual** of the closed form in (⋆): `max|resid| ∈ [1.5e−7, 1.4e−6]` over `u∈(0.2,0.8)` for
  `Δ_j=1,1.5,2,3` (finite-difference floor). The closed form solves the defect Casimir.
- **Independent series coefficient**: the ODE (⋆) alone gives `f'(0) = (Δ_j(Δ_j−1)+Δ_φ²)/(1−2Δ_φ)`;
  this matches the closed form's `f'(0)` to 6 digits (`−1.333333, −1.583333, −2.000000, −3.333333`).
- **Indicial exponents** of (⋆) reproduce `{Δ_j,1−Δ_j}` (u=1) and `{0,2Δ_φ}` (u=0).
- **Full-block normalisation check**: `⟨φφ⟩=(4ρ₁ρ₂u)^{−Δ_φ}f` satisfies
  `⟨φφ⟩/(1−u)^{Δ_j−Δ_φ} → const` as `ρ₁→0` (`Δ_j=3,Δ_φ=2`: `0.2095→0.2012→0.2024→0.2012` as `ρ₁` halves),
  confirming the `Δ_φ` prefactor and the defect-OPE exponent `Δ_j−Δ_φ`.

## 6. Honest boundary

- The shape `f(u)P_j(cosθ)` and one `SO(2,1)` invariant `u` are assumed from the parent docs; I did not
  re-derive the two-invariant structure.
- The cleanest closed derivation I could complete is **(a)+(b)**: defect-OPE ⇒ standard 1D block with
  external `(Δ_φ,0,0,Δ_φ)` ⇒ the `₂F₁` ⇒ its ODE (⋆), whose indicial exponents independently reproduce
  the OPE data. I could **not** re-derive the operator (⋆) directly from a bulk `(z,ρ)` Casimir by brute
  force: the na"ive total/single-point `SO(2,1)` Casimirs *annihilate* any `ρ`-dressed `u`-function (I
  verified `Cas[…]=0` in sympy), so the defect-channel "reduced" operator must be imported from the 1D
  block rather than obtained by acting the group Casimir on the two-point function. This is the one step
  I am flagging as not closed from first principles here.
- `Δ_φ = 2` is representative; any `Δ_φ` gives the same structure with `Δ_j−Δ_φ` in the `₂F₁` and the
  `(4ρ₁ρ₂u)^{−Δ_φ}` prefactor. No claim about the dynamical OPE coefficients `a_j` / large-`N` here.

## 7. Bottom line vs. the old question

`BLOCK-CORRECTED.md` framed the answer as a combination of `f₁=P` and `f₂=Q`. That premise is the
problem: those are harmonic H² partial waves solving the free KG equation, not defect conformal blocks.
The correct block is the single `₂F₁` above, solving the defect Casimir (⋆) with both OPE endpoints
`f(0)=1` and `f ∼ (1−u)^{Δ_j}` built in.

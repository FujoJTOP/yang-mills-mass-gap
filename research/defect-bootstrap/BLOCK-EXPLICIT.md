# Explicit defect blocks: closed forms, the impossibility, and the honest boundary

> Supersedes the "one closed form fixes both OPE ends" hope in `BLOCK-CORRECTED.md` (and the raw `f`
> of `U-ODE.md`). Convention is the one shared by all three notes:
> `u = tanh²(η/2)`, `z = cosh d = (1+u)/(1−u)`, `u ∈ (0,1)`; **`u→0` = bulk OPE**, **`u→1` = defect OPE**.
> All numerics: `verify_block.py` → `verify_block_OUT.txt`, mpmath at 40 digits, scipy cross-check
> (`eval_legendre`, `lqmn`) for integer degree.

---

## 1. Closed forms, and why no single function works

### PROVEN: the two independent H² solutions

On H², a function of the geodesic distance `d` satisfies `f'' + coth(d) f' = λ f`. The `SO(2,1)`
quadratic Casimir of the exchanged defect operator of dimension `Δ_j` acts with eigenvalue
`λ = Δ_j(Δ_j−1)` (`verify_block_OUT.txt` (c): machine-precision agreement, rel. err. ~1e−40). Writing
`λ = ν(ν+1)` fixes

```
ν = Δ_j − 1 ,       z = cosh d = (1+u)/(1−u) ,
```

and the two independent solutions are the Legendre functions of degree `ν`:

**Bulk-channel block (`P`-type, first kind) — regular at `u→0`:**
```
P_{Δ_j−1}(z) = ₂F₁(Δ_j, 1−Δ_j; 1; −u/(1−u))            [DLMF 14.3.6]
             = ₂F₁(1−Δ_j, Δ_j; 1; (1−z)/2) .
```
This is exactly `U-ODE.md`'s `f`, now correctly identified: `f` **is** `P_{Δ_j−1}`.

**Defect-channel block (`Q`-type, second kind) — decays at `u→1`:**
```
Q_{Δ_j−1}(z) = (√π Γ(Δ_j) / (2^{Δ_j} Γ(Δ_j+½))) · z^{−Δ_j} · ₂F₁((Δ_j+1)/2, Δ_j/2; Δ_j+½; 1/z²)
```
[DLMF 14.3.7]. The `z^{−Δ_j}` factor is the manifest `(1−u)^{Δ_j}` defect decay.

### PROVEN: the OPE-end behaviour (numerics, `verify_block_OUT.txt` (a),(b))

| block | `u→0` (bulk OPE, need `→1`) | `u→1` (defect OPE, need `(1−u)^{Δ_j}`) |
|---|---|---|
| `P_{Δ_j−1}` | `P → 1` ✓ | `P ~ z^{Δ_j−1} ~ (1−u)^{1−Δ_j}` ✗ **grows** |
| `Q_{Δ_j−1}` | `Q ~ −½ ln u + O(1)` ✗ **log-diverges** | `Q ~ (1−u)^{Δ_j}` ✓ |

Precise leading constants (verified; `γ` Euler, `ψ` digamma):
```
P_{Δ_j−1}·(1−u)^{Δ_j−1}  →  2^{2(Δ_j−1)} Γ(Δ_j−½)/(√π Γ(Δ_j))     (= {1, 1.27324, 2, 6} for Δ_j={1,1.5,2,3})
Q_{Δ_j−1}/(1−u)^{Δ_j}     →  √π Γ(Δ_j)/(2^{2Δ_j} Γ(Δ_j+½))         (= {0.5, 0.19635, 0.083333, 0.016667})
Q_{Δ_j−1} + ½ ln u        →  −γ − ψ(Δ_j)                            (= {0, −0.613706, −1, −1.5})
```
So `P` carries a **nonzero** `(1−u)^{1−Δ_j}` growing branch (`→ 2^{2(Δ_j−1)}Γ(Δ_j−½)/(√πΓ(Δ_j)) ≠ 0`),
while `Q` carries **none** (`Q·(1−u)^{Δ_j−1} → 0`).

### PROVEN: the exact connection between the two solutions (Ferrers)

Because `(z+1)/(z−1) = 1/u` **exactly**,
```
Q_{Δ_j−1}(z) = −½ ln(u) · P_{Δ_j−1}(z) − W_{Δ_j−2}(z) ,     W analytic,  W_{Δ_j−2}(1) = γ + ψ(Δ_j).
```
(Verified: `Q + ½ ln u · P → −(γ+ψ(Δ_j))` smoothly for all three test `Δ_j`.) This is the precise sense
in which the two blocks differ by a **non-analytic** `log(u)·P` term at the bulk end: they cannot be
related by any finite constant rescaling or combination that removes the log.

### PROVEN: impossibility of a single `α P + β Q` meeting both ends

Two facts, each numerically exhibited:

1. **Bulk finiteness (u→0) forces `β = 0`.** `P → 1` (bounded), `Q ~ −½ ln u → +∞`. Since
   `½ ln(1/u) = −½ ln u` exactly, any `β ≠ 0` gives `βQ → ±∞` that a bounded `αP` cannot cancel.
   (Numerics, Δ_j=2: `Q = {0.407, 1.349, 2.461, 3.606, 4.757}` at `u = {1e−1 … 1e−5}`, unbounded; `P→1`.)
2. **Pure defect decay (u→1) forces `α = 0`.** The growing `(1−u)^{1−Δ_j}` branch is present only in
   `P`, with nonzero coefficient `2^{2(Δ_j−1)}Γ(Δ_j−½)/(√πΓ(Δ_j))`; `Q` has none. A combination that is
   *pure* `(1−u)^{Δ_j}` must have `α·(nonzero) = 0`, i.e. `α = 0`.
   (Numerics, `u=0.9999`: `P·(1−u)^{Δ_j−1} = {1.273, 2.000, 5.999}`, `Q·(1−u)^{Δ_j−1} = {2e−9, 8e−14, 2e−22}`.)

`β=0` *and* `α=0` ⟹ the only combination is `f ≡ 0`. **The 2-dimensional solution space cannot meet two
incompatible boundary conditions.** There is no both-limits single function; none is fabricated here.

**Physical reading (not a defect in the maths).** The defect-OPE block is the `Q`-type solution
(unique: regular, pure `(1−u)^{Δ_j}`). The bulk-OPE block is the `P`-type solution (`→1`). They are
*different blocks of two different OPEs*, connected by the Ferrers relation above. The full correlator
is a **defect-OPE sum** `Σ_j a_{Δ_j}² Q_{Δ_j−1}(z)·P_j(cosθ)`; its bulk limit is controlled by the
**bulk-OPE sum**. "Both ends fixed by one closed form" is impossible; "both ends fixed by the two
blocks + their connection + the coefficient data" is the correct statement.

---

## 2. The `Δ_φ` normalization

Write the correlator as
```
⟨φ(x₁)φ(x₂)⟩_def = (ρ₁ρ₂)^{−Δ_φ} · F(u, cos θ) ,   F = Σ_j a_{Δ_j}² · Q_{Δ_j−1}(z) · P_j(cos θ) (+ bulk-channel terms).
```

- **PROVEN — the `u`-shape is `Δ_φ`-independent.** The H² radial ODE `f''+coth(d)f' = Δ_j(Δ_j−1)f`
  contains only `Δ_j` (the exchanged defect dimension); `Δ_φ` does not appear. Hence the functional
  form `P_{Δ_j−1}(z)` or `Q_{Δ_j−1}(z)` is fixed by `Δ_j` alone.
- **`Δ_φ` enters through (i) the prefactor `(ρ₁ρ₂)^{−Δ_φ}`** (engineering dimension / the `ρ^{−Δ_φ}`
  scaling at each point), and **(ii) the dynamical coefficients `a_{Δ_j}²`** (defect OPE
  coefficients, tied to `C^{φφ}_{O_j}`). These two are where all `Δ_φ`-dependence lives.
- **Bulk-OPE "→1" is not a property of any single block.** It is recovered only (a) in the **full
  `Δ_j`-sum** (the defect blocks must conspire, using the bulk-OPE data `(Δ_φ, C^{φφ}_O)`, to cancel the
  `Q`-type logs and reproduce the free propagator), equivalently (b) **after multiplying by
  `x₁₂^{2Δ_φ}`** (dividing out the bulk two-point). A single `Q_{Δ_j−1}` keeps its `−½ ln u`
  divergence; a single `P_{Δ_j−1}` keeps its `(1−u)^{1−Δ_j}` growth. Neither alone `→1`.

Mark: the ODE's `Δ_φ`-independence is **PROVEN**; the exact normalization/coefficients that make the
*sum* reproduce `→1` are **CONVENTION / not closed here** (they are the dynamical bootstrap data, not
kinematics).

---

## 3. Numerical verification table (`Δ_j × u`)

From `verify_block_OUT.txt` (mpmath, 40 digits; identical to scipy `eval_legendre`/`lqmn` for integer
`Δ_j` to ≤1e−16). Columns: `P`, and its defect-end growth probe `P(1−u)^{Δ_j−1}`; `Q`, its defect-end
decay probe `Q/(1−u)^{Δ_j}`, and its bulk-end log probe `Q+½ln u`. The **limits** (last row) are what
the probes converge to; the finite-`u` values illustrate the approach.

| `Δ_j` | `u` | `P` | `P(1−u)^{Δ_j−1}` | `Q` | `Q/(1−u)^{Δ_j}` | `Q+½ ln u` |
|---|---|---|---|---|---|---|
| 1 | 0.05 | 1.00000 | 1.0000 | 1.49787 | 1.5767 | 0.0000 |
| 1 | 0.95 | 1.00000 | 1.0000 | 2.5647e−2 | 0.5129 | 0.0000 |
| 1.5 | 0.05 | 1.03884 | 1.0125 | 0.95704 | 1.0336 | −0.5408 |
| 1.5 | 0.5 | 1.59739 | 1.1295 | 0.11289 | 0.3193 | −0.2337 |
| 1.5 | 0.95 | 5.62443 | 1.2577 | 2.2809e−3 | 0.2040 | −0.0234 |
| 2 | 0.05 | 1.10526 | 1.0500 | 0.65554 | 0.7264 | −0.8423 |
| 2 | 0.5 | 3.00000 | 1.5000 | 3.9721e−2 | 0.1589 | −0.3069 |
| 2 | 0.95 | 39.000 | 1.9500 | 2.1924e−4 | 0.0877 | −0.0254 |
| 3 | 0.05 | 1.33241 | 1.2025 | 0.33788 | 0.3941 | −1.1600 |
| 3 | 0.5 | 13.000 | 3.2500 | 5.4567e−3 | 0.04365 | −0.3411 |
| 3 | 0.95 | 2281.0 | 5.7025 | 2.2490e−6 | 0.01799 | −0.0256 |
| **limit** | `u→0` | **`→1`** ✓ | — | **`→+∞`** ✗ | — | **`→−γ−ψ(Δ_j)`** |
| **limit** | `u→1` | **`→∞`** ✗ | **`→2^{2(Δ_j−1)}Γ(Δ_j−½)/(√πΓ(Δ_j))`** | — | **`→√πΓ(Δ_j)/(2^{2Δ_j}Γ(Δ_j+½))`** | — |

Both "wrong-end" failures are visible: `P` **grows** at `u→1` (`P(1−u)^{Δ_j−1}` stays finite ≠ 0),
and `Q` **log-diverges** at `u→0` (`Q+½ln u` stays finite while `Q→+∞`). The small-`u` convergence of
`Q+½ln u` to `−γ−ψ(Δ_j)` is slow (log corrections) — pinned at `u=1e−5` in
`verify_block_OUT.txt` (e.g. Δ_j=2: −0.999885 vs −1). By contrast the defect-end probes converge fast.

---

## 4. Free-limit calibration vs `FREE-CODIM3.md`

`FREE-CODIM3.md` (free scalar, `Δ_φ=1`, codim-3) expands the free propagator in transverse partial
waves `f_l ∝ (2l+1) Q_l(cosh η) P_l(cos θ)` with integer degree `l`, and notes `l=0` gives
`Q₀ = −½ ln u` ("matching the 1D `Δ=1` log"). Its closed form is exactly the `ν=l` case of my `Q` above
(`(√π Γ(l+1)/(2^{l+1}Γ(l+3/2)))z^{−l−1}₂F₁((l+2)/2,(l+1)/2;l+3/2;1/z²)` equals `Q_l(z)`).

- **PROVEN:** the degree is `ν = Δ_j − 1`; at the free point the partial-wave label is `l = ν`, so the
  natural identification is **`Δ_j = l + 1`** (equivalently `l = Δ_j − 1`). With `Δ_j=l+1` my
  `Q_{Δ_j−1} = Q_l`, reproducing `FREE-CODIM3`'s `f_l`. The `l=0 ↔ Δ_j=1` anchor (the `−½ ln u` log) is
  consistent with this.
- **CONVENTION / ambiguous:** `FREE-CODIM3`'s `l` is the **SO(3)⊥ partial-wave label of the free
  propagator**, not independently an OPE dimension; `Δ_j = l+1` is the natural continuation, but the
  free two-point is the *trivial* defect, so it does not by itself fix which operator (unit `Δ=l+1`,
  or a double-trace `[φφ]_l` of dimension `2+l`) sits at each `l`. For **non-integer** `Δ_j` there is no
  free anchor at all: the offset `Δ_j = l+1` is a continuation choice. The unambiguous, PROVEN content
  is only `degree ν = Δ_j − 1`.

---

## Status ledger

| claim | status |
|---|---|
| Casimir eigenvalue `Δ_j(Δ_j−1)`, degree `ν=Δ_j−1` on H² | **PROVEN** (ODE check, 1e−40) |
| `P_{Δ_j−1} = ₂F₁(Δ_j,1−Δ_j;1;−u/(1−u))`; `Q` closed form | **PROVEN** (DLMF + 40-digit match) |
| `P→1` (u→0), `P~(1−u)^{1−Δ_j}` (u→1); `Q~(1−u)^{Δ_j}` (u→1), `Q~−½ln u` (u→0) | **PROVEN** (numerics) |
| No `αP+βQ` (nonzero) meets both ends; only `f=0` | **PROVEN** (2-fact argument + numerics) |
| Ferrers connection `Q = ½ln(1/u)P − W`, `W(1)=γ+ψ(Δ_j)` | **PROVEN** (numeric) |
| `u`-shape `Δ_φ`-independent; `Δ_φ` in prefactor + `a_{Δ_j}²` | **PROVEN** (ODE) / coefficients CONVENTION |
| exact normalization making the defect-OPE *sum* `→1` (bulk end) | **OPEN / CONVENTION** (dynamical data) |
| free offset `Δ_j = l+1` | **CONVENTION** (natural; operator id & non-integer continuation ambiguous) |

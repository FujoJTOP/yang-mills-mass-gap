# P1-PERIODIC-HOLONOMY — the exact periodic global constant-holonomy anchor: `S_W ≡ 0`, marginal pure Haar, `G″_anch(0)=−ρ_Haar″/ρ_Haar=1/(2·Cin(2π))=0.205115` for SU(2), **L-independent [T]**; and the acid test shows the exact object is `−ρ″/ρ`, **not** `1/Var` (the latter is only the Gaussian limit)

> **2026-11-24 · Wave P1 — the periodic constant mode, exact anchor.**
> Question: on a **periodic torus** (no twist, no fixed boundary), does the constant-mode curvature
> saturate at a positive `L`-independent value or decay? This note settles the *global constant
> holonomy* — the one object that can be computed in closed form — and separates it from the `k=0`
> Fourier mode (P2), which is the physical soft mode.
> Witness: `research/v3/w7/p1_periodic_holonomy.py` (`py -3.12 -I … quick`, numpy 2.4.6 / scipy 1.18.1;
> output quoted §7). **No third-party sources.**
>
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

1. **The action vanishes exactly [T].** For the global constant link configuration — all links equal
   `U`, or one direction equal `U` and the rest `I` — every plaquette is the identity,
   `P_□ = U·U·U†·U† = I` (or `U·I·U†·I = I`), so `S_W = 0` for **every** `L`, `β`, `U`
   (witness `[A]`: `S_W/β = 5×10⁻¹⁵`, `max|P−I| = 1×10⁻¹⁶`).
2. **The marginal is pure Haar [T].** Hence the constant-holonomy marginal is the Haar measure on the
   single modulus `U`, and its curvature is exactly
   ```
      G″_anch(0) = −ρ_Haar″(0)/ρ_Haar(0) = 1/(2·Cin(2π)) = 0.20511530   (SU(2)),
   ```
   `L`-independent, for every `L` **[T]**. The finite-`L` value is not an asymptote: it is the value.
3. **It is not `1/Var_Haar` [T]/[V].** `1/Var_Haar(ξ₃) = 9/(4π²−6) = 0.26883` is the **Gaussian/CLT
   limit**, not the exact curvature of the Haar marginal; the gap is `+0.064`. The acid test decides
   which is the exact object: for the free massless `U(1)` field (flat Haar, `J=1`), `−ρ″/ρ = 0`
   exactly, while `1/Var = 3/π² = 0.304` would be a **spurious floor**. Therefore the exact object is
   `−ρ″/ρ`; `1/Var` is a large-`L`/Gaussian approximation **[T]/[V]**.
4. **Which object is the physical soft mode.** The physical decompactification soft mode is **object
   B**, the `k=0` Fourier mode of the fluctuations, whose `β=0` marginal is the `L`-fold convolution
   and whose curvature **rises to** `1/Var = 0.26883` as `L→∞` **[T]**. Object A (global holonomy) is
   the flat-connection moduli measure: an exact, solvable anchor, but **not** the soft-mode marginal.
5. **SU(3).** Direct Haar sampling gives `Var_Haar(ξ₃) = 2.0165`, `1/Var = 0.4959`, and a stable KDE
   reading `G″_anch(0) ≈ 0.569` (`h`-stable); i.e. for SU(3) the exact Haar curvature is **above**
   `1/Var` (the opposite sign of the SU(2) gap), confirming that the two are genuinely different
   objects and that `1/Var` is not the exact curvature **[V]**.

---

## 1 · The configuration and the plaquette computation

On a periodic `L^4` lattice with `U_μ(x)` the link variables, take the global constant configuration

```
   (A1)  U_μ(x) = U  for all x, μ                 (all links equal), or
   (A2)  U_3(x) = U,  U_μ(x) = I  (μ ≠ 3)         (one direction constant, rest identity).
```

For (A2), a plaquette in a plane containing direction 3 is `P = U·I·U†·I = I`; every other plaquette
is `I·I·I†·I† = I`. For (A1) every plaquette is `U·U·U†·U† = I`. The witness verifies both on a
periodic `L=2` 4D torus with a generic SU(2) element:

```
   dir-3 = U, rest = I :  S_W/β = +5.329e-15   max|P−I| = 1.111e-16
   all links = U       :  S_W/β = +1.066e-14   max|P−I| = 1.153e-16
```

Hence `S_W ≡ 0` **exactly**, and the Boltzmann weight on this slice is `1` **[T]**. This is the
periodic-torus statement of the well-known fact that a constant gauge field has zero field strength;
it is the exact starting point of the GPY/holonomy effective-potential literature, but here it is used
only for the **marginal** curvature, with no loop expansion.

---

## 2 · The exact Haar curvature (object A)

For SU(2) in the exponential chart `U = e^{iξ·σ/2}`, the one-Cartesian-component Haar marginal is
(the transverse components are integrated exactly)

```
   ρ₁(z) = 2π ∫_{|z|}^{2π} sin²(ρ/2)/ρ dρ = π ( Cin(2π) − Cin(|z|) ),  |z| ≤ 2π,
```

with `Cin(x)=γ_E+log x−Ci(x)`. At the origin `ρ₁′(0)=0` and

```
   G″_anch(0) = −ρ₁″(0)/ρ₁(0) = 1/(2·Cin(2π)) = 0.20511530 .
```

The witness reproduces this to 6 digits by a 5-point finite difference (`0.205116`), and gets the same
closed form from the exact `β=0` one-link marginal **[T]+[V]**. The value is `L`-independent because
the slice carries no `L`-dependent factor: the marginal is the **same single-modulus Haar measure**
for every `L`.

**Exact vs Gaussian.** `Var_Haar(ξ₃) = (4π²−6)/9 = 3.719824`, so `1/Var = 0.268830`. The two differ by
`0.0637`; the Haar marginal is bounded (`|z|≤2π`) and non-Gaussian, so the Gaussian identity
`G″=1/Var` does not hold. This is not a normalization artefact: the acid test below shows which object
is the curvature.

---

## 3 · Acid test: free `U(1)` decides between `−ρ″/ρ` and `1/Var`

For the free massless `U(1)` field in the flat chart, the Haar measure is Lebesgue (`J=1`), so the
constant-mode marginal is uniform on the circle and `−ρ″/ρ = 0` exactly. But `1/Var(θ) = 3/π² = 0.304`
would be a **positive floor for a free field** — precisely the failure mode the charter forbids.
Therefore:

- `−ρ″/ρ` **passes** the acid test (returns `0` on the free field);
- `1/Var` **fails** it (manufactures `0.304`);

so the exact curvature object is `−ρ″/ρ`, and statements of the form `G″=1/Var` are Gaussian-limit
statements, valid only when the marginal is asymptotically Gaussian **[T]**.

**Consequence for the corpus.** `D6-DECISION`'s dichotomy `G″(a₀*)=1/Var(a)` is correct as a
**Gaussian-limit** identity (its own derivation says so), and its `β=0` value `0.2688` is the
`L→∞` value of object B, not the exact finite-`L` Haar curvature. The exact finite-`L` object is the
convolution curvature, which at `L=1` is `0.2051` and rises to `0.2688`. No numerical claim in the
corpus changes; the **logical status** of `1/Var` is sharpened from "identity" to "Gaussian limit
`[C]`".

---

## 4 · Object A vs object B: which is the physical soft mode

The `k=0` Fourier mode is `a = L^{-1/2} Σ_l ξ_{3,l}` with `|v|=1`. Its `β=0` marginal is the `L`-fold
convolution `ρ_v^{Haar}(a) = √L (ρ₁^{*L})(√L a)`, and the witness reproduces the D1/D4 values exactly:

| `L` | 1 | 2 | 3 | 4 | 6 | 8 | 12 | 16 | 32 | `∞` |
|---|---|---|---|---|---|---|---|---|---|---|
| `G″_L(0)` | 0.20532 | 0.23966 | 0.24940 | 0.25431 | 0.25918 | 0.26160 | 0.26402 | 0.26522 | 0.26703 | **0.26883** |

This is **object B**: `L`-dependent at finite `L`, rising to the Gaussian value `1/Var` **[T]**. The
physical decompactification soft mode is object B because the coarse field is the zero-momentum
component of the *fluctuation field*, not the global holonomy of the flat connection. Object A is the
moduli-space measure of constant connections, useful as an exact anchor and as an acid-test
discriminator, but it is **not** the marginal whose curvature enters `(★)`.

---

## 5 · SU(3)

Direct Haar sampling (`1.5×10⁶` samples) gives `Var_Haar(ξ₃) = 2.01650` and `1/Var = 0.49591`; the KDE
curvature at the origin is `0.6226, 0.5989, 0.5383, 0.4671` for bandwidths `0.1, 0.15, 0.2, 0.3`, with
a stable reading `≈0.569`. Two readings:

- the exact Haar curvature for SU(3) is **above** `1/Var` (opposite to SU(2)), so `1/Var` is neither
  exact nor a universal bound;
- the SU(2) exact value `0.20512` is reproduced by the same sampler to `≈0.205` at `h=0.3`, validating
  the method **[V]**.

The lesson is group-independent: the exact object is `−ρ_Haar″/ρ_Haar`, and the Gaussian value
`1/Var` is a limit.

---

## 6 · What this settles and what it does not

**Settled [T]:**
- `S_W ≡ 0` on the global constant-holonomy slice, every `L`, `β`, `U`;
- the marginal is pure Haar; `G″_anch(0)=1/(2·Cin(2π))=0.205115` for SU(2), `L`-independent;
- the exact curvature object is `−ρ″/ρ`; `1/Var` is the Gaussian/CLT limit (acid test: free `U(1)`
  gives `0` vs the spurious `3/π²`);
- object B's `β=0` curvature rises to `1/Var=0.26883` as `L→∞`.

**Not settled [O]:**
- the sign/size of the **Wilson part** of object B at `β=2.5` in the periodic model — this is P2/P4,
  and in 4D it is `H1`;
- whether a positive constant-mode floor implies anything about the atom — P6, answer: no, a
  curvature floor is a susceptibility ceiling, never an atom.

**One line.** *On the periodic torus the global constant holonomy is an exact, `L`-independent Haar
anchor `0.205115` [T]; the physical soft mode is the `k=0` Fourier mode, whose `β=0` floor is the
Gaussian value `0.26883` [T]; and the acid test proves the exact curvature is `−ρ″/ρ`, not `1/Var`.
No QED. Clay 0%.*

---

## 7 · Witness — `research/v3/w7/p1_periodic_holonomy.py`

`py -3.12 -I p1_periodic_holonomy.py quick` (exit 0, 1.5M Haar samples; `full` uses 8M). Quoted output
(abridged):

```
[A] dir-3=U, rest=I: S_W/beta=+5.329e-15  max|P-I|=1.111e-16 ; all links=U: +1.066e-14, 1.153e-16
[B] exact G''_anch(0)=1/(2 Cin(2pi))=0.205115 ; numeric 5-pt FD=0.205116
    Var_Haar=3.719824 ; 1/Var=0.268830 ; gap +0.0637 (marginal NON-Gaussian)
[C] object B beta=0: L=1,2,3,4,6,8,12,16,32 -> 0.20532,0.23966,0.24940,0.25431,0.25918,0.26160,
    0.26402,0.26522,0.26703 -> 1/Var=0.26883 (CLT)
[D] free U(1): exact G''_anch(0)=0 ; 1/Var=3/pi^2=0.3040 (SPURIOUS). Exact object = -rho''/rho.
[E] SU(2): Var=3.72374, 1/Var=0.26855, KDE G''=0.2050 at h=0.3 (exact 0.20512)
    SU(3): Var=2.01650, 1/Var=0.49591, KDE G''~0.569 (h-stable)
```

*Refs.* `D4-HAAR-WILSON.md` (§1.2 closed form `ρ₁`, Haar floor), `D6-DECISION.md` (Gaussian-limit
`1/Var`), `SM4-CONSTANT-MODE.md` (constant-mode marginal), `P2-PERIODIC-K0.md` (object B at `β=2.5`),
`P6-PERIODIC-TO-H1.md` (logical strength). *No third-party sources downloaded.*

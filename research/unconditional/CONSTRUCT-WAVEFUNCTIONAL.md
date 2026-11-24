# CONSTRUCT-WAVEFUNCTIONAL — a controlled-variance `0⁺⁺` wavefunctional **does** close the Temple sandwich, but its construction *is* the eigenstate: variance-control ⟺ solving the theory (= `H1`)

> **VERDICT (no QED).** PNT-TEMPLE-LOWER showed Temple is dominated by the Rayleigh quotient and named
> the missing input: a *controlled-variance* `0⁺⁺` wavefunctional (variance `≪` gap²) plus spectral
> separation. **Built and tested on a genuinely interacting solvable model** (1D anharmonic oscillator
> `H = −½∂² + ½ω²x² + λx⁴`, `λ=0.5`): the Temple sandwich **closes** — a trial with `Var(H)≈3.2×10⁻³`
> certifies **92%** of the true gap (`E₁−E₀ = 1.628`, certificate `1.50`). **But the same trial is a small
> perturbation off the exact first-excited eigenvector.** The sharp dichotomy (witness):
> **generic Gaussian families carry a variance floor `4.64×10⁻²` and certify NO gap (Temple ratio `5.47`)**,
> while **variance → 0 ⟺ trial → exact eigenstate ⟺ ratio → 0**. So the "controlled-variance `0⁺⁺`
> wavefunctional" *exists* but is *the solution itself*. In 3+1D that is `H1`.

---

## §0 · Labels
**[T]** 真定理 · **[V-num]** new numeric · **[C]** conditional (given `E₀^L, W`) · **[NG]** no-go for the named class · **[O]** exact missing step. No QED.

---

## §1 · The Temple sandwich and its ratio (recap)

Temple (1928): for normalized `ψ`, `ρ=⟨H⟩`, `σ²=Var H=⟨H²⟩−ρ²`, `E₀≤ρ<E₁`, and `W≥E₁` an upper bound on `E₁`,
```
E₁ ≥ ρ − σ²/(W − ρ).
```
It certifies a **gap** `E₁ ≥ E₀^L` iff the **Temple ratio**
```
R := σ² / ((W − ρ)(ρ − E₀^L))  <  1 .
```
So a positive lower bound on the `0⁺⁺` gap needs three inputs: **(i)** a certified lower bound `E₀^L ≤ E₀`
(hard); **(ii)** an upper bound `W ≥ E₁` (Rayleigh on an orthogonal trial — easy); **(iii)** a trial `ψ`
with `E₀<ρ<E₁` and **small variance** `σ²`. The variance is the bottleneck.

---

## §2 · The construction on the anharmonic oscillator (interacting, solvable)

Exact spectrum by dense ED in the harmonic-oscillator basis (no factorial overflow); trial functions
evaluated in spectral (FFT) position space (`w=1`).

| model | `E₀` | `E₁` | gap |
|---|---|---|---|
| **free** `λ=0` (acid) | 0.500000 | 1.500000 | 1.000000 |
| **interacting** `λ=0.5` | 0.696176 | 2.324406 | 1.628231 |

### (a) Generic Gaussian family — the variance floor blocks the sandwich
Best Gaussian `0⁺⁺`-type trial (`even` for `ρ`, `odd` for `W`), scanned over width `s∈[0.45,1.9]`:

| model | generic residual `Var(H)` floor | best generic Temple bound | gap certified | ratio |
|---|---|---|---|---|
| **free** `λ=0` | `5.6×10⁻¹⁷` (≈0) | — (trivial) | — | — |
| **interacting** `λ=0.5` | `4.64×10⁻²` | `0.6588 < E₀^L=0.6962` | **NO** | `5.47` |

**Free field:** the exact ground state *is* a Gaussian ⇒ variance floor `≈0` ⇒ controllable variance is
**trivial and vacuous** (the acid baseline: no interaction, no mechanism).
**Interacting field:** **no Gaussian is an eigenstate**, so the Gaussian family carries a hard variance
floor `4.64×10⁻²`; the best generic Temple bound falls *below* `E₀` ⇒ **certifies no gap** (ratio `5.47`).

### (b) Controlled variance — perturbing the exact eigenvector closes the sandwich
Take `ψ = normalize(U₁ + εU_k)` (small admixture off the exact first-excited eigenvector):

| `ε` | `Var(H)` | Temple bound | gap certified | ratio |
|---|---|---|---|---|
| 0.30 | `2.83×10⁻¹` | 0.2336 | **no** (`<E₀`) | 1.31 |
| 0.10 | `3.54×10⁻²` | 1.3832 | 0.687 | 0.574 |
| **0.03** | `3.24×10⁻³` | **2.1947** | **1.50** | **0.079** |

At `ε=0.03` the sandwich closes and certifies `E₁ ≥ E₀ + 1.50` against the true `E₁−E₀ = 1.628` — **92% of
the gap**, rigorously, given `(E₀^L, W)`. As `ε→0`, `Var→0`, `ratio→0`, and the certificate `→ E₁` exactly.

### (c) The dichotomy (the round's content)
```
generic Gaussian family:   Var floor 4.6e-2  ->  ratio 5.47  ->  NO gap
controlled variance:       Var -> 0          ->  ratio -> 0  ->  gap -> exact
```
**A controlled-variance `0⁺⁺` wavefunctional exists and closes the Temple sandwich — but every such
wavefunctional is (a small perturbation of) the exact eigenstate.** The "controlled variance" datum `≡`
knowledge of the eigenvector.

---

## §3 · Honest caveats
- The calibration used `E₀^L = E₀` (exact) as a stand-in for a *certified* lower bound. A fully rigorous
  certificate needs `E₀^L < E₀` from an independent method (ground-state Temple/Weinstein or interval
  arithmetic) — itself a lower-bound problem of the same kind.
- `W` came from the odd-trial variational energy (`2.348` vs true `E₁=2.324`, a small overshoot) — a genuine
  upper bound, as required.
- The model is `0+0`-dimensional (single mode). The *mechanics* of the Temple sandwich and the
  variance bottleneck are dimension-independent; the gauge-specific content enters only at §4.

---

## §4 · The 3+1D assessment
The anharmonic oscillator is the cleanest *interacting* solvable instance, and it already exhibits the
relevant structure. In 3+1D pure `SU(N)` YM the analogous object is the gauge-invariant `0⁺⁺`
wavefunctional on the physical Hilbert space `H_phys = ker Q_BRST`. The controlled-variance requirement
`⟨H²⟩ − ⟨H⟩² ≪ (gap)²` becomes:

- a **variational `0⁺⁺` trial** with a *certified small* `Var(H)`, i.e. a trial whose overlap with the
  exact `0⁺⁺` eigenstate is `≈1`, together with a nonperturbative evaluation of `⟨H²⟩` on `H_phys`;
- **plus** a certified lower bound `E₀^L` on the vacuum energy;
- **plus** an upper bound `W ≥ E₁`.

But "a trial within `Var ≪ gap²` of the `0⁺⁺` eigenstate" *is* the eigenstate to that accuracy, and the
nonperturbative `⟨H²⟩` on `H_phys` is the transfer-matrix second moment — the theory itself. **So the
controlled-variance `0⁺⁺` wavefunctional exists (existence: it is the eigenvector) but its construction is
exactly the `0⁺⁺` eigenvalue problem = the mass gap = `H1`.** The Temple route does not reduce difficulty; it
*relocates* the missing input from "a lower bound on `E₁`" to "a `0⁺⁺` wavefunctional with certified
sub-gap variance," which is the same wall in the wavefunctional coordinate.

---

## §5 · Bottom line
- **Temple does not close the sandwich with generic trial families** (`[NG]` for Gaussian/ansatz families):
  the variance floor `4.6×10⁻²` (interacting) ⇒ ratio `5.47` ⇒ no gap; free field is trivially exact (`[T]`).
- **A controlled-variance `0⁺⁺` wavefunctional DOES close it** (`[V-num]`/`[C]`): at `Var=3.2×10⁻³` it
  certifies 92% of the gap — **but only as a near-exact eigenvector**.
- **The controlled-variance datum ⟺ eigenstate knowledge ⟺ solving the model**; in 3+1D this is `H1` (`[O]`).
  Existence is not the issue; construction is.
- **Clay 0%; no QED.** The Temple/wavefunctional coordinate is another face of the same wall.

*Ref:* prev `PNT-TEMPLE-LOWER.md`, `ROOF-TENSOR-RG.md`; witness `w7/construct_wavefunctional.py`
(runs `py -3.12 -I`, ~<5 s, numpy only).

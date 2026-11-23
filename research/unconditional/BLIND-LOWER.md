# BLIND-LOWER — A genuine lower bound on the continuum Yang–Mills gap?

*Independent derivation, first principles. No project files consulted.*

We are granted: 4D pure SU(N) YM exists as a Wightman QFT, with transfer matrix
`H ≥ 0`, vacuum `Ω` at energy 0, and the gap `Δ = inf(spec H ∖ {0})`. Strong coupling
(`β` small) is proven gapped. The question: a **lower** bound `Δ ≥ m > 0` on the
physical theory, ideally uniform across the Savvidy crossover.

## 0. Why the direction is the whole problem
Spectral data `G(t) = ∫ e^{-λt} dμ(λ)`. If `G(t) ≤ A e^{-mt}` for all t, then the
bottom of `supp μ` lies above `m`, i.e. **an upper bound on the correlator is the
lower bound on the gap**. Reflection positivity, checkerboard/cs monotonicity, and
all "fat" correlation inequalities instead give `G(t) ≥ c e^{-Δt}` — **lower**
bounds on `G`, hence **upper** bounds on `Δ`. They cannot yield `Δ ≥ m`. The only
positivity-type statement that yields a lower gap is a *massive-Gaussian /
Dirichlet-form domination* of the interacting measure: an **upper** bound of
log-Sobolev / Poincaré type. Keep that distinction fixed; it decides every route.

## 1. Temple / Weinstein — circular
Temple gives `λ₁ ≥ ρ − (⟨H²⟩−ρ²)/(β − ρ)` where `ρ` is a trial Rayleigh quotient and
`β ≤ λ₂` a **lower bound on the NEXT eigenvalue**. To bound `Δ = λ₁` we must already
possess a lower bound on `λ₂` (the "gap of the gap"), and so on ad infinitum. In a
continuum QFT the ladder has no anchor: there is no `λ_n` at infinity to bootstrap
from (Weyl governs large `λ`, the wrong end). Every variational lower-bound engine is
therefore either **conditioned on the next eigenvalue** or reduces to the same gap one
level up. Verdict: circular.

## 2. Comparison / coupling — local in β
GKS/FKG coupling needs a total order or ferromagnetic (positive) couplings; non-abelian
lattice gauge has neither, so the clean comparison theorems do not apply. Polymer /
cluster expansions (Brydges–Federbush–Kennedy; Balaban's small/large-field RG) do prove
`Δ_lat(β) ≥ c(β) > 0`, but `c(β)` degrades monotonically and the proof requires sign-
definite (positive) plaquette action — exactly the quadratic positivity that **flips at
Savvidy**. The comparison is genuine but **only at strong coupling**; it does not cross.

## 3. Weyl / spectral-support — wrong regime and no scale
Weyl's law `N(λ) ~ C λ^{d/2}` counts the **UV** spectrum (large `λ`). The gap is an **IR**
statement about the bottom of the support, which is dimensionally transmuted: `Δ = κ·Λ`
with `κ` a pure number. No scale-free spectral asymptotic fixes `κ` from below; a Weyl
argument constrains the density of large eigenvalues, never the existence or size of a
lowest one. "No eigenvalue below `m`" is not a Weyl statement.

## 4. RG monotones — flow constraints, not IR magnitude
The a-theorem (`a_UV > a_IR`, Komargodski–Schwimmer) and QNEC/ANEC constrain how
central charges and entanglement fall along a flow. They rule flows *in or out*; they do
not bound the IR mass scale `Δ/Λ` from below. A monotone that "wants" a gapped IR still
allows the gap to be arbitrarily small relative to the cutoff. Wrong observable.

## 5. The one honest lower-bound engine — and its certificate
Log-Sobolev / Bakry–Émery is a real lower-bound engine: if `Ent_μ(f²) ≤ 2C·E_μ(f,f)` for
the Euclidean gauge measure `μ` (Dirichlet form `E_μ`), then OS reconstruction gives a
mass gap `Δ ≥ 1/(2C)`. This is an **upper** bound on `C`, i.e. the right direction. The
standard certificate for `C` is a **curvature lower bound**: `Ric + Hess S ≥ K > 0`
gives `C ≤ 1/K`. Here `S` is the gauge action and its second variation is
`½∫ tr ξ(−D² + 2\,ad F_B)ξ`. At the **Savvidy crossover** the transverse magnetic mode
makes this Hessian **indefinite** (`−D²+2adF_B` has negative eigenvalues), so the natural
curvature certificate `K` is not positive and not uniform in `β`, cutoff, or background;
`C(β) → ∞`. The engine's output `Δ ≥ 1/(2C)` collapses to `0` precisely at the wall.
This is not a defect of log-Sobolev; it is the honest content: the only lower-bound engine
we have certifies a gap only where a **positive-definite comparison form** exists.

## 6. Reduction to ONE inequality
Everything above funnels into a single statement. Prove:

> **(H)** There is `κ > 0` such that, uniformly in cutoff `Λ` and `β`, the
> *physical* gauge-reduced second variation of the effective action — the Hessian
> `−D² + 2\,ad F_B` **plus the one-loop/measure back-reaction** that resums the negative
> magnetic mode into a positive-definite weight on gauge-invariant transverse modes —
> satisfies `⟨ξ, (−D²+2\,adF_B)_{phys} ξ⟩ ≥ κ‖ξ‖²` for every unit transverse mode `ξ`.

If (H) holds, Bakry–Émery upgrades it to a uniform LSI constant `C ≤ 1/κ`, whence
`Δ ≥ κ/2` uniformly across all `β`. The obstruction is therefore one inequality: the
**gauge-reduced effective Hessian is positive-definite with a β-uniform lower bound.**
The naive Hessian is `~7%` negative; the entire difficulty is whether the measure/determinant
back-reaction provably lifts it to positive-definite **without** the lift depending on the
scale that sets `Δ` (the circularity trap: a lift sized `∝ κ` proves nothing).

## Verdict
- **PROVED:** strong-coupling lattice gap `Δ_lat(β) ≥ c(β) > 0` (`β < β₀`); 2D YM exactly
  gapped; a-theorem, LSI/Poincaré inequalities as *conditional* engines. Genuine, but local.
- **CONDITIONAL:** (H) ⇒ uniform LSI ⇒ `Δ ≥ κ/2` across the Savvidy crossover. This is the
  whole Clay-adjacent gap.
- **NO-GO (for the naive tools):** Temple/Weinstein (circular — needs `λ₂`); correlator
  positivity (upper bound on `Δ`, wrong way); Weyl (UV, dimensionally blind); RG monotones
  (flow, not scale); cluster expansion/comparison (sign-definite β, fails at Savvidy).
- **OPEN:** a genuine `Δ ≥ m` for the continuum theory across the crossover. Absent (H),
  every available "lower-bound" engine is either **local in β** (works only strong-coupling)
  or **circular** (consumes the gap it produces). The single precise missing inequality is (H):
  β-uniform positivity of the gauge-reduced effective Hessian.

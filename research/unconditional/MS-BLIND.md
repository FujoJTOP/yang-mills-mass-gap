# MS-BLIND — Transporting the strong-coupling gap across the Savvidy crossover

*Independent first-principles attempt; no other project files consulted.*

**Granted.** 4D pure SU(N) YM exists as a Wightman QFT in the a→0, L→∞ limit of the
Wilson lattice (transfer matrix T = e^{-H}, H≥0; Källén–Lehmann measure ν≥0). At strong
coupling β<β₀ a convergent cluster expansion gives a gap m(β) ≥ m₀ > 0 and a positive
residue f₁²(β) > 0. The physical point is β→∞ on the asymptotic-freedom curve,
a(β) ≈ Λ⁻¹·(…)·e^{-cβ}. There the one-loop Hessian −D² + 2ad F_B is **indefinite**:
in a constant chromomagnetic background the lowest transverse mode sits at −gB < 0
(Savvidy), and the negative part is *extensive* — per unit volume Σ_{λ<0}|λ| ≍ (gB)².
Target: show liminf_{β→∞} m_lat(β)/a(β) = m_phys > 0 uniformly.

## 1. Five candidate mechanisms and the exact death of each

1. **RG flow of the gap.** m is dimensionful: m(μ) = m(μ₀)(μ/μ₀)^{γ_m}. Transport needs
   γ_m bounded below along the whole path and integrable; γ_m is known only perturbatively
   at the two ends, and the crossover is nonperturbative. No transport.

2. **Gap equation with a positive fixed point.** The exact closed systems (Makeenko–Migdal
   loop equations, Schwinger–Dyson) are infinite hierarchies. Any truncation admitting a
   fixed point m\* > 0 controls the true m only to the truncation error, which is O(1) at
   the crossover where the condensate back-reacts on the propagator. m\* is a property of
   the truncation, not a theorem about the theory. Conditional at best.

3. **Analyticity + zero-free region.** m(β) is analytic only on the strong-coupling disk
   |β<β₀| (radius set by nearest complex singularity). The Savvidy point lies outside; no
   zero-free region is known across it, and even the *absence of a phase transition* for
   SU(N≥2) in 4D is unproven (numerically smooth, not theorem). Cannot cross.

4. **Variational / monotone.** Rayleigh quotients give **upper** bounds on individual
   eigenvalues; the gap λ₁−λ₀ admits no variational lower bound without a prior gap
   (Temple–Kato are circular). Any FKG-type monotonicity in β would yield
   m(β) ≤ m(β_small): again an upper bound, wrong direction.

5. **The cluster expansion itself.** Its engine is a stable quadratic minorant:
   H″(A) ≥ ε(−∂²) − C with ε>0 *uniform*. At Savvidy ε cannot be chosen >0 — the quadratic
   form has an extensive negative spectrum — so no sign-definite polymer weight exists and
   the expansion cannot even begin. This is the sharpest form of the obstruction.

## 2. What *does* transport: the sign of residues

A residue is a squared matrix element with multiplicative renormalization,
f₁²(μ) = Z_O(μ)² f₁²(μ₀). For a single (non-mixing) operator Z_O² > 0, so **sgn f₁² is
RG-invariant along the entire flow, crossover included**, provided γ_O is integrable (no
pole) and mixing with lower-dimension operators is triangular with positive diagonal.
This is the only robust transport available — and it is conditional on single-operator
covariance. It transports positivity, not size.

## 3. Why the residue is not the gap

f₁² > 0 does **not** imply m > 0: a massless excitation carries f₁² > 0 at a pole s=0.
The residue is the *intensity* of spectral weight; the gap is its *lower support edge*.
Multiplicative positivity protects the former and is structurally blind to the latter —
hence this report's name.

## 4. The precise obstruction

Spectral positivity (ν≥0, supp ⊂ [m²,∞)) gives moment ratios M_{k+1}/M_k ≥ m² — an
**upper** bound on m². Variational methods give upper bounds. The cluster expansion is the
only source of a *lower* bound, and it dies at Savvidy (§1.5). Every robust tool is thus
either **one-sided** (upper bounds) or **local** (a β-disk). The transport needs the exact
opposite: a *global lower* bound. Stated cleanly, the missing inequality is

> **(SGL)**  m⁴ ≥ κ·⟨¼ tr F²⟩,  with κ>0 universal and ⟨¼ tr F²⟩>0 for all β,

i.e. the gap is bounded below by the fourth root of the gluon condensate. Positivity gives
only the reverse-sign family (m² ≤ moment ratio via OPE/sum rules). Nothing known implies
SGL. Consequently transporting the proven condensate *and* residue positivity yields a
positive spectral weight but **not** a positive gap: the pole may sit at zero.

## 5. Verdict — proved / conditional / open

- **PROVEN.** Strong-coupling gap m≥m₀>0 (granted); one-sidedness of positivity and
  variational methods (upper bounds on m); multiplicativity of Z_O² (sign preservation of
  residues).
- **CONDITIONAL.** Sign-transport of f₁² across the crossover (needs integrable γ_O and
  triangular/positive mixing); existence of a gap-equation fixed point m\*>0.
- **OPEN — the genuine obstruction.** (SGL) / the uniform lower bound
  liminf m_lat/a > 0. The crossover is not merely "hard": it is the precise point where the
  only lower-bound engine (a stable quadratic minorant) fails because the Savvidy negative
  mode is extensive, while every remaining tool — positivity, moments, variational,
  monotonicity, multiplicative RG — is structurally one-sided or local and can certify
  upper bounds and sign preservation only. The gap is the *location* of a pole, which no
  positivity principle protects.

**Net.** A nonperturbative transport exists for the *residue's sign*, not for the *gap*.
Turning residue positivity into a gap positivity is exactly inequality (SGL); that single
inequality is the crossover obstruction.

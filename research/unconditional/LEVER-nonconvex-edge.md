# LEVER-nonconvex-edge — can a NON-CONVEX functional give a genuine LOWER bound on the spectral edge?

> 2026-11-23 · The "missed lever" of CONVERGENCE-redteam §(b)3, attacked head-on. **Verdict: it does not exist from present provable data.** The four non-convex routes collapse by four *distinct* mechanism-tags — direction, hypothesis, and insensitivity — and two of them (Temple, Combes–Thomas) do not merely fail: their hypotheses *are* the gap. No construction survives; the exact open door is named. No fabricated proof.

**The target object.** `Φ` on positive measures `ν` (KL measure, `ν≥0`, `supp ν⊂[m²,∞)`), built from protected data + provable bounds, with `Φ ≤ inf supp ν` on every admissible theory, and **not** implied by finitely many moments.

## (0) Two model-independent walls
**(A) Free-field adversary.** Any `Φ` whose value is fixed by the OS/reflection-positivity axioms alone is `≤ 0`: the **free massless scalar** (`ν=δ₀`, gap 0, propagator `1/4π²x²`) satisfies every OS axiom. So `Φ>0` requires injecting YM dynamics. The only proven YM-specific data — trace-anomaly coefficient, stress-tensor normalisation, asymptotic freedom, Ward/SD identities — are UV/total (moments, ⇒ upper, GAP-G2-edge) or unclosed (SD). The lever cannot live in "protected data."

**(B) Linear-data lemma (rigorous).** Fix any finite set `∫φ_i dν = C_i` (moments, generalised moments, f-sum, susceptibilities). For any `ν`, any `ε`, and any neighbourhood `U∋0` there is `ν'` with an atom in `U` and the same `C_i` up to `O(ε)`. Hence no `Φ=Φ(C_i)` lower-bounds `inf supp ν` uniformly. This sharpens GAP-G2-edge's `δ₀+δ₂` vs `2δ₁` from moments to **all fixed linear spectral data**.

## (1) Temple / Kato–Weinstein — *direction incompatibility*
Temple gives a true lower bound **only given `ρ ≤ E₂`**. The **only** variational output is the opposite: `E₂ ≤ RQ(ψ₂)`. The two directions are incompatible — substituting a variational second state into Temple is **unsound** (`ρ ≥ E₂` violates the hypothesis; it can return a false, too-high bound). There is no self-contained min–max lower bound for a bottom eigenvalue: Courant–Fischer is a *min* (upper); Weinstein's bound — the sharpest — also needs a spectral separation. In the mass-operator frame this is a **one-rung regress**: a lower bound on `E_{n+1}` needs `E_{n+2}`; the ladder has no top except a proven continuum threshold (two-scale separation, GAP-G2-edge). Min–max carries a lower bound *upward from a known base*; there is no base at the top.

## (2) Determinant / Rayleigh on a restricted class — *monotone the wrong way*
Restricting the class `𝒞→𝒞'⊂𝒞` only **raises** `min_{φ∈𝒞'}RQ`, i.e. worsens the upper bound; it can never flip sign. `RQ(φ) ≥ m²` for every `φ`, so any min over subclasses stays an upper bound. Determinants are total (trace-type) functionals; every positive-operator determinant inequality (Hadamard, Fischer, Lidskii) bounds `det` from above given diagonals, or needs the other eigenvalues from below — the same regress. Infinite-volume `det` → Fredholm/zeta-regularised = moments.

## (3) Double commutator / Combes–Thomas / ξ — *hypothesis = gap*
The f-sum gives `M₁=½⟨[[H,A],A†]⟩`, `M₀=⟨{A†,A}⟩_conn`, so `gap ≤ M₁/M₀`: again **upper**. A lower bound needs the **Bogoliubov inequality** `gap·|⟨[A,B]⟩|² ≤ double-commutator`, which needs a **local order parameter** `⟨[A,B]⟩≠0`. Pure YM has none (confinement is not SSB of a local symmetry; Polyakov loop is non-local, finite-T). Combes–Thomas localisation is an *inside-a-gap resolvent* statement: it yields rate `ξ` **given** `ξ<∞=1/gap`. A massless scalar obeys the same local commutator bounds with `ξ=∞`, so "ξ finite controlled independently" is the gap, renamed.

## (4) Sign-definite sum rule — *insensitivity*
A sum rule is a fixed linear functional of `ν`; by (B) it cannot pin the edge. A proven-negative-weight rule is still one linear constraint. The fully *nonlinear* spectral consequence of RP — `D(p²)=∫dν/(p²+s)` is Stieltjes / completely monotone — **admits** a massless atom `c/p²`, `c>0` (Widder). Complete monotonicity does not forbid the edge atom. No proven-sign YM nonlinear inequality of the needed direction exists.

## The one door, and why it closes on the gap
A lever would need a **proven-sign, YM-specific nonlinear inequality** forbidding the transverse massless mode *and* a **lower** bound on its residue `f₁²`. But the transverse-vertex residue is symmetry-unprotected — no anomalous slot (GAP-G2-mass §a) — so only dynamics fixes its sign, and "`f₁²>0`" *is* the Schwinger order parameter, i.e. the gap. Even the strongest conceivable edge statement, `D(0)=0`, is strictly weaker (Gribov–Zwanziger: `D(0)=0 ⇏ gap`) and still unproven in continuum.

**Sharpest line.** Every non-convex attempt collapses to the gap by **direction** (Temple, f-sum, determinants: upper), by **hypothesis** (CT `ξ`, Bogoliubov order parameter, Temple `E₂`: gap assumed), or by **insensitivity** (linear data, Stieltjes: massless atom admissible). No `Φ ≤ inf supp ν` exists from present provable data; constructing one is proving a sign-definite YM spectral inequality — the gap.

**Still true:** Temple/Weinstein/f-sum/Bogoliubov/CT are correct theorems; the conditional residue-`⇒`-gap implication is clean; strong-coupling lattice gaps are real. None yields a continuum lower bound.

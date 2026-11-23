# INVENT-GAUGE-INVARIANT-MASS — the 0⁺⁺ spectral edge of the field-strength correlator: right object, unprovable sign

> **Verdict. A genuine gauge-invariant mass is constructible — the spectral edge `m_GI = inf supp ρ_s` of the OS-positive two-point function of `s = tr F²` — and it is non-vacuous: `0` for free Maxwell, `= m_{0⁺⁺}` for Yang–Mills.** It is manifestly gauge-invariant, BRST-closed, positive-norm, which **removes both** obstructions that killed every gauge-fixed route: the sign-indefinite `ρ_A` (Oehme–Zimmermann) and `A ∉ H_phys` — because `s ∈ ker Q_BRST`. The transfer to `m_GI > 0` still **fails**, now for one reason: the **edge atom**. OS positivity is Schur-concave, controls only moments of `ρ_s`, never its edge; free Maxwell satisfies every axiom with `m_GI = 0`. The construction converts "the mass is not gauge-invariant" into "the gap *is* the edge" — a strict and final reduction, no fabricated theorem.

## 1. Why the prompted combination fails
`⟨A_min²⟩ := ⟨inf_g‖A^g‖²⟩` is gauge-invariant (orbit functional; Singer bypassed) but a **scale** (`∼Λ²`), not a spectral edge: the minimiser leaves the transverse modes, so `⟨A_min²⟩ > 0` for the **massless** field too (acid test fails), and it carries no Källén–Lehmann weight. `γ_G` is tied to the Landau horizon — gauge-fixing-dependent (the `A^h` completion is a proposal, not a theorem). `u(0)` is a covariant-gauge **confinement witness**, not an edge — and the mass-bearing decoupling branch has `u(0) ≈ −0.7…−0.8 ≠ −1` (Fischer). Mass and confinement parameter sit on *opposite* branches; the triple does not cohere.

## 2. The construction (the missing object)
On `H_phys = ker Q_BRST / im Q_BRST = L²(𝒜/𝒢, μ̄)`, let `s(x) = (β(g)/2g) tr F_{μν}F^{μν}(x)` — the trace-anomaly density: gauge-invariant, BRST-closed, dim 4, and in the Wilson-loop algebra (small-loop expansion), whose transfer matrix is `T` on `L²(𝒜/𝒢, μ̄)`.

```
DEFINITION (gauge-invariant IR mass)
   G_s(x) = ⟨s(x) s(0)⟩_c = ∫_0^∞ ρ_s(μ²) Δ(x; μ²) dμ² ,   ρ_s ≥ 0   (OS)
   m_GI   = inf supp ρ_s = lim_{|x|→∞} −∂_{|x|} ln|G_s(x)|
```

`ρ_s ≥ 0` since `s` is Hermitian and in `H_phys`; `m_GI` is the prefactor-independent decay rate.

**Transfer (the theorem).** `T = e^{−aH}` is self-adjoint on `H_phys`, Perron vector `Ω` (PROOF-GAUSS-INVARIANT-PF). The spectral theorem in the cyclic subspace of `s` gives the Källén–Lehmann form with `m_GI = E_1^{(s)} − E_0`. **Conditional theorem:** *if* `⟨s(x)s(0)⟩_c ≤ C e^{−m|x|}`, *then* `m_GI = m > 0`, the lightest `0⁺⁺` mass coupled to `F²` — the full gap iff the lightest bulk state is `0⁺⁺`-coupled (GAP-DICTIONARY (b); lattice-true, not symmetry-protected).

## 3. Acid-test result
- **Free Maxwell (`g→0`):** `s` is quadratic → two-photon states; `ρ_s` continuous on `[0,∞)` from threshold `0` (soft/collinear photons, total momentum → 0). `m_GI = 0`. ✓ (Finite-volume `T` is *still* gapped — Krein–Rutman is **not** the order parameter; only the thermodynamic edge is.)
- **Yang–Mills:** gapped `0⁺⁺` ⇒ `m_GI = m_{0⁺⁺} > 0`. ✓

Nonzero for YM, zero for free: the object is **non-vacuous**.

## 4. Exactly where it fails — the edge atom
Both gauge-fixed obstructions are gone: `ρ_s ≥ 0` (no Oehme–Zimmermann flip — that struck `ρ_A` only because `A ∉ H_phys`) and `s ∈ ker Q`. What remains is purely spectral — and it is the wall. Every Cauchy–Schwarz/Hölder/moment/OS inequality is Schur-concave and sees only moments and total mass. `inf supp ρ_s` is non-convex edge data: `M_0 = M_1 = 2` is realised by both `δ_0 + δ_2` (edge 0) and `2δ_1` (edge 1), so moments never pin the edge (GAP-G2-edge), and free Maxwell realises the same axioms with a continuum from `0` (GOLDS-POS-MASS). Perron–Frobenius anchors the *vacuum* — itself `0⁺⁺` — so the glueball lies in `Ω^⊥`: PF seals the channel top, never the edge. Hence **no theorem `m_GI > 0` exists**; proving it is the 4D dynamics itself.

**Sharpest line.** The gauge-invariant mass is the spectral edge of the OS-positive `tr F²` correlator — the correct, acid-test-valid order parameter, with the `ρ_A`-sign and BRST-cohomology objections *removed*; but positivity is Schur-concave and cannot force a positive-mass atom, so the transfer stops at the edge, and free Maxwell is the exact witness. Gauge invariance was never the wall; the edge is.

*Refs.* Osterwalder–Seiler Ann. Phys. 110 (1978) 440; Lüscher CMP 54 (1977) 283; Källén 1952; Lehmann 1954; Oehme–Zimmermann PRD 21 (1980) 1661; Kugo–Ojima PTPS 66 (1979) 1; Gribov NPB 139 (1978) 1; Zwanziger NPB 323 (1989) 513; Singer CMP 60 (1978) 7; Kondo PRD 58 (1998) 105019 (hep-th/9805153, `A^h`); Fischer arXiv:0812.2745; Morningstar–Peardon PRD 60 (1999) 034509. Internal: PROOF-GAUSS-INVARIANT-PF; GOLDS-POS-MASS; GAP-G2-mass; GAP-G2-edge; EDGE-f1sq-attack; PROOF-LANDAU-GAUGE-IR; PROOF-FRG-DSE-GAP-EQUATION; GLUEBALL-BOOTSTRAP; GAP-DICTIONARY; gauge-invariance/REPORT.

# FINITE-VOL-CERTIFY — a certified finite-volume gap at 4D SU(2) β=2.5: the exact barrier

> 2026-11-24 · Witness `research/unconditional/w7/finite_vol_certify.py` (`py -3.12`, single job).
> Weights: **W4 / PROOF-COMPUTER-ASSISTED / FORGE5-CERTIFIED-TN / BS-NUMERICAL-RIGOROUS / PROOF-FSS-RIGOROUS / OUT-certified-computation.**

Labels: **[T]** rigorous theorem · **[C]** conditional theorem · **[O]** open / no-go ·
**[V-here]** verified numerically in this work (code above) · **[R]** prior registry result.

---

## VERDICT

**There is no certified finite-volume 4D SU(2) mass-gap VALUE at any physically
meaningful volume (L ≥ 2) at β = 2.5. The barrier is exact and quantified, not a
lack of effort.** The only spatial volume small enough to certify by dense linear
algebra is the **degenerate single-site L = 1** (no glueball can live there).

The *certified computation that does close rigorously here* is the **free-field box
gap** through the full certificate pipeline (a-posteriori residual enclosure,
error ≈ 10⁻³⁰) — which is exactly the **acid sanity**: a fixed-`(a,L)` free massless
field already has a positive box gap, so a fixed-`(a,L)` gap is **not** the Clay gap.
The value-certification machinery is genuine; the 4D-interacting object it would
certify is unreachable at β = 2.5 for L ≥ 2.

| object | largest certifiable volume | certified value | status |
|---|---|---|---|
| free massless scalar box gap (acid test) | L = 4 (trivially any L) | `1.41421356237 ± 10⁻³⁰` (L=4) | **[T]** + **[V-here]** |
| 4D SU(2) β=2.5 gap, L=1 (degenerate) | L = 1³ (dim = 8470 at j≤3) | dense-feasible, *not assembled* (recoupling build) | **[O]** |
| 4D SU(2) β=2.5 gap, L=2 | 2³: dim(j≤1/2)=**13,251,039** | **INFEASIBLE** | **exact barrier [V-here]** |
| 4D SU(2) β=2.5 gap, L=3,4 | — | **INFEASIBLE** (dim ≫ 10¹⁰) | **exact barrier** |

**No QED.** A certified finite-volume gap is a real result, but it is **not** the
Clay mass gap: (i) at fixed `(a,L)` the free massless field is also gapped (proved
below), and (ii) `Δ_V ≥ Δ_∞ = inf_L Δ_V`, so a *lower* bound on a finite-volume gap
is an **upper** bound on the continuum gap — the wrong direction. **The Clay gap
remains open; this deliverable is a frontier map, not a step across it.**

---

## §1 Setup

**Theory.** 4D pure SU(2) Wilson lattice gauge on the torus `T⁴_L = L_t × L³`
(spatial `L³`), spacing `a = 1`, Wilson action `S = β Σ_p (1 − ½ Re tr U_p)` with
the standard normalisation `β = 2N/g² = 4/g²` (N=2). At the physical crossover
**β = 2.5 ⇒ g² = 1.6, ½g² = 0.8, 1/g² = 0.625** (SU(2) bulk crossover is a
*analytic* crossover near β≈2.2–2.5, not a phase transition).

**Transfer matrix / gap.** `T = e^{−aH_V}`, with `H_V` the Kogut–Susskind
Hamiltonian on the spatial slice,
```
H_V = (g²/2) Σ_{links} C₂(ℓ)  +  (1/g²) Σ_{plaq} (N − Re tr U_p),     C₂(j)=j(j+1),  N=2.
```
`T` acts on the gauge-invariant Hilbert space `ℋ = L²(G^{E_V})^{G^{V}}` (Haar), which
is **infinite-dimensional even at L = 1** (`G = SU(2)`; Peter–Weyl `⊕_j V_j ⊗ V_j^*`).
The **mass gap in lattice units** is `Δ_V = E₁ − E₀ = −ln(λ₁/λ₀)`, the lightest
gauge-invariant excitation (the 0⁺⁺ channel is the lightest in pure YM; we take the
overall gap, which is 0⁺⁺). A finite-volume `Δ_V` is one number; the continuum
transfer gap is `Δ_∞ = inf_L Δ_V`.

**Certified-gap definition used here.** A *certified lower bound*
`Δ_V ≥ Δ_cert > 0` requires a **rigorous enclosure** `E₀ ∈ [a,b]`, `E₁ ∈ [c,d]`
with `c > b`, giving `Δ_cert = c − b`. Dense diagonalisation alone is *not* a
certificate; you also need (i) interval/verified eigenvalue bounds on a **finite**
matrix and (ii) a **rigorous tail bound** relating the finite matrix's spectrum to
the true infinite-dimensional `T`. Both are mandatory; the second is the wall.

---

## §2 The certificate — what closes, with numbers

### §2.1 Rigorous, certified, and non-trivial: the free-field box gap **[T][V-here]**

The free massless scalar single-particle Hamiltonian `M = −Δ_lattice` on `L³`
(periodic, 7-point stencil) is an **exact integer matrix** (position space,
**non-diagonal**, banded). We compute its two lowest eigenvalues with the verified
eigensolver `mpmath.eigsy` at 60 decimals and certify each by the a-posteriori
residual bound: for a unit eigenvector `v`, `ρ = ⟨v,Mv⟩`, the spectral theorem gives
an eigenvalue within `‖Mv − ρv‖₂` of `ρ`. The residuals are `≈10⁻⁶⁰`; the two
enclosures are disjoint (`ρ₁ − ρ₀ > r₀ + r₁`), so the box gap is **certified**:

| L | dim | `μ₀` enclosure | `μ₁` enclosure | certified box gap `√μ₁−√μ₀` | exact `2 sin(π/L)` |
|---|---|---|---|---|---|
| 2 | 8  | `[−1.6e-60, 5.5e-61]` | `[4,4]` | **`2.0000000000 ± 7×10⁻³¹`** | 2.0 |
| 3 | 27 | `[−8.3e-61, 3.0e-60]` | `[3,3]` | **`1.7320508076 ± 1×10⁻³⁰`** | 1.732050808 |
| 4 | 64 | `[−5.3e-60, 1.0e-60]` | `[2,2]` | **`1.4142135624 ± 1×10⁻³⁰`** | 1.414213562 |

This is the **first-of-its-kind certified box-gap value** in this program — with an
explicit controlled error (`~10⁻³⁰`) and a fully verifiable pipeline. Its content is
the acid test (§4): the **free massless field is gapped in a box**. **[T]**

### §2.2 The 4D SU(2) β=2.5 gap — why no value is quoted at L ≥ 2

Dense diagonalisation of the *actual* gauge-invariant Kogut–Susskind matrix needs
the finite truncation `P_Λ` (per-link rep ceiling `C₂ ≤ Λ`). Two independent walls
(§3) forbid a **certificate** at L ≥ 2. At **L = 1** the matrix *is* small (dim = 15
at j≤½, 421 at j≤3/2, **8470 at j≤3** — §3) and fits in 0.6 GB, but
(i) the volume is a **single site** (three commutator-plaquettes `½Re tr[U_i,U_j]`,
no glueball), and (ii) assembling the full SU(2) electric-basis Hamiltonian requires
the plaquette recoupling (6j) engine; that is a genuine but separate engineering
task, and its output would be a physically meaningless number. **We do not
manufacture it; the honest deliverable is the barrier below.** **[O]**

---

## §3 Feasibility frontier — the exact barrier

We compute two decisive quantities. All numbers reproduced by the witness.

### §3.1 Gauge-invariant Hilbert-space dimension `dim_GI`

For a spatial lattice with per-link rep ceiling `J` (link space `⊕_{j≤J} V_j`), the
gauge-invariant dimension is the lattice-gauge integral
`dim_GI = ∫ ∏_s dg_s ∏_{links} Σ_{j≤J} χ_j(g_s)χ_j(g_t)`.

* **L = 1 (EXACT, Gauss–Legendre, [T][V-here])**:  `dim_GI = 15 (j≤½), 99 (j≤1),`
  `421 (j≤3/2), 1359 (j≤2), 8470 (j≤3), 34497 (j≤4), 107745 (j≤5)`.
* **L = 2, j≤½ (EXACT enumeration of all 2²⁴ link assignments, [T][V-here])**:
  `dim_GI = 13,251,039` (**1.3×10⁷**) — even at the *minimal* truncation.
* **L = 3, j≤½ (MC, order-of-magnitude [V-here])**: `dim_GI ~ 10¹⁴` (MC variance is
  O(mean); the rigorous content is `≫ 10⁴`, far past the dense ceiling).
* **L = 4**: ≫ 10²⁰.

**Dense-certified ceiling.** A dense `dim×dim` matrix needs `dim²·8` bytes; a 16 GB
workstation caps `dim ≲ 3×10⁴` for eigh, and **verified/interval** dense LA (the
certificate needs interval arithmetic + a validated eigensolver, 4–8× overhead) caps
`dim ≲ 10⁴`. Hence:

| lattice | `dim_GI` (min trunc.) | dense bytes | vs certified ceiling (~10⁴) |
|---|---|---|---|
| **L=2** (24 links) | **13,251,039** | **1.4 PB** | **×10³ infeasible** |
| L=3 (81 links) | ≳10¹⁴ | 10²⁸ B | hopeless |
| L=4 (192 links) | ≫10²⁰ | — | hopeless |

So certified dense diagonalisation is **infeasible at L = 2 by a factor ≳10³** in
dimension (and the sparse-Lanczos route gives eigenvalues but **not a certificate**,
since it supplies neither intervals nor a rigorous tail bound).

### §3.2 Electric-tail (Feshbach) certificate — vacuous at β=2.5

The only rigorous route to a genuine operator-norm lower bound
(PROOF-COMPUTER-ASSISTED; Weyl/Kato): keep reps `C₂ ≤ Λ`, use the electric tail
floor on the discarded sector `γ(J) = (g²/2)·C₂(J+½) = 0.8·(J+½)(J+3/2)`, and the
Schur/Weyl correction `η = ‖P V_mag Q‖²/(γ − E_max)`. The certificate is non-vacuous
only if `η < Δ_V ∼ O(1)`. The magnetic boundary coupling is
`‖P V Q‖² ≃ (1/g²)² · N_plaq · (2J+1)` (one fundamental per plaquette link, reduced
matrix element grows as `√(2J+1)`). Solving `η<1` at β=2.5:

| L | `n_plaq`=3L³ | `J_min` for `η<1` | `dim_GI` at that `J_min` |
|---|---|---|---|
| **1** | 3 | **3.0** | **8470** (dense-feasible) |
| **2** | 24 | **22.5** | ≳10¹³ (min trunc. already 1.3×10⁷) |
| 3 | 81 | 78.0 | ≫10²⁰ |
| 4 | 192 | 186.5 | ≫10²⁰ |

**Exact barrier.** At L ≥ 2 the electric-floor certificate demands a rep ceiling
`J_min ∼ O(10–10²)`, at which `dim_GI` is ≳10¹³ — astronomically beyond the 10⁴
certified ceiling. (The `J_min` figures use the standard Kogut–Susskind normalisation
above; they shift by an O(1) if the trace/½ convention changes, but the barrier does
**not** depend on them: the exact `dim_GI(L=2, j≤½)=1.3×10⁷` alone already exceeds
the ceiling by ≳10³, independent of any Hamiltonian normalisation.) At the *minimal*
truncation j≤½ (dim 1.3×10⁷, already infeasible) the floor is only `γ=1.6` while
`‖V_mag‖ ∼ N_plaq/g² ∼ 15`, so `η ≫ Δ_V` — the certificate is **vacuous even if the
matrix were assembled**. Both walls are the same object: the electric-energy floor
and the magnetic coupling compete at fixed β=2.5, and the gauge-invariant dimension
grows as `~exp(c L³ log d)`.

### §3.3 What could change the frontier (honest)

* **Sparse/structured + a genuine operator-norm tail bound**: blocked — the only
  operator-norm lower-bound device (Weyl on a PEPS/TN truncation) has a
  **volume-exponential** `ε(D)` in 4D and is sign-blind (FORGE5-CERTIFIED-TN [R]);
  no d≥2 rigorous TN error theorem exists.
* **Larger `J` to rescue `γ`**: makes `dim_GI` worse (exponential).
* **L = 1 with the recoupling engine**: feasible (dim 8470) but degenerate — not a
  gap of the 4D theory.

---

## §4 Acid sanity — the method is *faithful*, and it is why this is not Clay

**[T][V-here].** On the *same* `(a,L)` the free **massless** field is run through the
*same* certified pipeline and returns a **positive** box gap: `2.0 ± 10⁻³¹ (L=2)`,
`1.732… ± 10⁻³⁰ (L=3)`, `1.414… ± 10⁻³⁰ (L=4)`, matching the exact `2 sin(π/L)` to
30 digits. Therefore:

1. **A positive certified finite-volume gap does not discriminate YM from a free
   field** — any correct faithful finite-volume pipeline must return one for both.
   So a fixed-`(a,L)` gap certificate, however rigorous, **cannot** be the Clay
   statement, and a pipeline that *failed* to return it for the free field would be
   broken (unfaithful). This is the "free-field acid test" of the Charter, passed.
2. The gaps **decrease with L** (`2.0 → 1.73 → 1.41 → …→ 2 sin(π/L) → 0`): the box
   gap is a *box artefact* that vanishes as the box grows — precisely why
   `Δ_V ≥ Δ_∞ = inf_L Δ_V` and a lower bound on `Δ_V` is an **upper** bound on the
   continuum answer (PROOF-FSS §1, BS-RIGOROUS [R]).
3. **Faithfulness**: the pipeline returns the exact analytic spectrum where known;
   the residual is `~10⁻⁶⁰`, so the certificate is tight, not slack.

---

## §5 What this does and does not prove

**Does prove ([T]/[V-here]):**
* The **certified eigenvalue-inclusion pipeline** (build exact finite matrix →
  a-posteriori interval enclosure → certified gap + explicit error) is correct and
  tight: it reproduces the exact free-field box gap to `±10⁻³⁰` on non-diagonal
  integer matrices. **[T]**
* **Exact frontier numbers**: `dim_GI(L=1)=15…107745` (j≤½…5, exact);
  `dim_GI(L=2, j≤½)=13,251,039` (exact); and the electric-floor `J_min` vs dimension
  table. These **quantify** the barrier. **[V-here]**
* A faithful finite-volume method returns a **positive** box gap for the free
  massless field — the acid sanity. **[T][V-here]**

**Does NOT prove / does NOT deliver:**
* **No** certified finite-volume gap value for 4D SU(2) at β=2.5 for any **L ≥ 2** —
  certified dense LA is infeasible (L=2 minimal dim 1.3×10⁷ ≫ certified ceiling 10⁴;
  electric-floor certificate vacuous at tractable `J`). This **confirms** the prior
  structural no-go (PROOF-COMPUTER-ASSISTED, FORGE5, BS-RIGOROUS, PROOF-FSS) with
  concrete numbers.
* **Not Clay.** Even a perfect fixed-`(a,L)` certificate is `Δ_V ≥ Δ_∞` — an *upper*
  bound on the continuum gap. The free field is gapped at the same `(a,L)`.
* **L=1** is the only dense-feasible volume and is degenerate (single site); its full
  SU(2) Hamiltonian was **not** assembled (a separate recoupling task), and it would
  not be a gap of the 4D theory.

**Sharpest line.** The certificate machinery is real and closed (free-field box gap
`±10⁻³⁰`, acid test passed). But at β=2.5 the only object it can *faithfully*
certify at a physical volume — the interacting gauge gap — sits behind a
**quantified exact barrier**: the minimal L=2 gauge-invariant space is already
`1.3×10⁷`-dimensional (≈10³× the certified-dense ceiling), and the electric-tail
Feshbach bound is vacuous at every tractable truncation. `Δ_V > 0` is automatic and
proves nothing about `Δ_∞`; the step across the barrier is uniform-in-`(V,a)` control
= the mass gap itself.

---

## Reproduce

```
py -3.12 research/unconditional/w7/finite_vol_certify.py     # single job, ~10 s (incl. L=2 exact enum)
```
Output: free-field box gaps (§2.1); certified enclosures (§2.1/§4); exact `dim_GI`
L=1,2 + MC L=3 (§3.1); electric-floor `J_min` table (§3.2); verdict.

*External references.* Wilson, PRD **10** (1974) 2445; Kogut–Susskind, PRD **11**
(1975) 395 (Hamiltonian form, β=2N/g²); Creutz, *Quarks, Gluons and Lattices* (1983);
Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Temple, Proc. R. Soc. A **119**
(1928) 276; Kato, *Perturbation Theory*; Davis–Kahan, SIAM J. Numer. Anal. **7**
(1970) 1; Rump (IntLab) / Neumaier (interval methods).
*Internal.* `PROOF-COMPUTER-ASSISTED` (the engine + the 24-link wall);
`FORGE5-CERTIFIED-TN` (no operator-norm TN lower bound in 4D); `BS-NUMERICAL-RIGOROUS`
(lower-bound direction is a structural no-go); `PROOF-FSS-RIGOROUS` (`Δ_V ≥ Δ_∞`,
wrong sign); `OUT-certified-computation` (finite-V = upper bound); `A5-RIGOROUS-NUMERICS`.

# Crossover Preservation and the Weak–Strong Handoff for Four-Dimensional Pure Yang–Mills: A Reduction and a Map of Obstructions

**Yuxuan Jiang**

*Preprint · 22 November 2026 · prepared for Zenodo*

---

## Abstract

We reduce the existence and mass-gap problem for four-dimensional pure Yang–Mills to a **weak–strong handoff inequality** `q₀ ≥ q₁`, formulate the two-sided *crossover preservation* statements `CPB-W` and `CPB-S` that would discharge it, and prove the handoff algebra in a **machine-checked** form. We then report, with proofs and falsification tests, that every mechanism currently available for `CPB-W`/`CPB-S` fails at the crossover `q ≈ 1`: Bianchi/sign cancellation is absent (`κ_c = 1` provable in the strong-coupling representation; no additive defect charge in 4D non-abelian gauge theory), Gaussian-weighted class preservation is blocked by a cluster constant (`C_cum ≥ 2` against a required `≤ 1.6`), and the numerical proxies that appeared to suggest otherwise are Gaussian tails of a blocked field that is itself Gaussian.

The outcome is not a proof of the mass gap. It is a precise statement of the obstruction: at `q ≈ 1` the single-block activity must be small, either because the defect action is large or because signs cancel, and no known structure supplies either. We also record the surrounding **conditional theorem**: ten named inputs, nine glue lemmas, nine machine-checked blocks, and four explicit external declarations.

**Keywords:** Yang–Mills, mass gap, lattice gauge theory, renormalization group, crossover, Bianchi identity, Lean 4, numerical falsification.

---

## 1. Introduction

The existence and mass gap of four-dimensional pure Yang–Mills is one of the Clay Millennium Problems. A recent construction (submitted in 2026 and assessed by an independent committee) yields a **subsequential infinite-volume limit** and OS reconstruction, but asserts neither a mass gap, nor uniqueness, nor confinement. The starting point of this paper is to isolate the **connective part** between existence and the gap — the *handoff* — and to write it as an inequality whose inputs can be checked one by one.

This paper **does not** prove the mass gap. It does three things:

1. **Reduction.** The mass gap is reduced to an explicit, dimensionally consistent handoff inequality (Section 2), whose algebraic core is proved and machine-checked.
2. **Precise two-sided statements.** The two crossover-preservation statements `CPB-W` (weak side) and `CPB-S` (strong side) are written in falsifiable form (Section 3), with a proof skeleton and its load-bearing lemma (Section 4).
3. **Map of obstructions.** Every known mechanism is shown or judged to fail at the crossover (Sections 5–6), and the remaining paths are listed (Section 7).

Section 8 records the surrounding **conditional theorem** and its machine-verification status, so that the reader can judge the step "if the inputs hold, then the statement holds".

---

## 2. The handoff inequality

### 2.1 Setting and constants

Fix a compact, connected, simply connected simple Lie group `G`; let `h^∨` be its dual Coxeter number, `dim f` the dimension of the defining representation, and `λ_G := h^∨ g²` the RG-invariant coupling. In the physical normalization `T(f) = 1/2`,

```
b₀ = 11 h^∨ / (48π²) ,      B := 2b₀ = 11 h^∨ / (24π²) ,
1/g²(μ) = B ln(μ/Λ) + O(g²) ,      C_RG := 24π²/11 .
```

**Notational discipline** (enforced throughout): the measure weight carries `γ := β/N_c = 2/g²` with `N_c := dim f`; the running coupling uses `λ_G = h^∨ g²`. The two coincide only for `SU(N)`, where `h^∨ = dim f = N`. Mixing them changes `q₀` by a factor `e^{3.59} ≈ 36` at `g² = 1`.

### 2.2 Theorem 1 (equivalent forms of the handoff)

Let `q(μ) := m/μ` (with `m` the physical mass gap, RG-invariant), let `μ_up` be defined by `g²(μ_up) = g*²`, and set `q₀ := q(μ_up)`; let `q₁` be the strong-side threshold. Then the following are equivalent:

```
(i)   q₁ ≤ q₀ ;
(ii)  ln(m/Λ) ≥ t_up(λ*) + ln q₁ ;
(iii) M := ln(m/Λ) − t_up(λ*) − ln q₁ ≥ 0 ;
(iv)  t_λ ≤ t_S   (the two domains meet in the (t, ln q) plane).
```

Here `t_up(λ*) = 24π²/(11λ*)` at one loop. **Geometric reading:** the weak side covers `t ≥ t_λ`, the strong side covers `t ≤ t_S`; full coverage holds iff `t_λ ≤ t_S`. This is pure real analysis and has been formalized in Lean 4 (`handoff_iff_margin`, `required_lambda`; axioms `propext`, `Classical.choice`, `Quot.sound`).

### 2.3 Theorem 2 (reach ⇒ explicit mass lower bound)

Suppose weak-side class preservation holds for `λ ≤ λ_c` (reach `W`) and the strong-side signed gas holds for `q ≥ q_c` (reach `S`). Then the handoff holds if and only if

```
m/Λ ≥ q_c · exp(t_up(λ_c)) .
```

**This is the central reduction.** The mass gap is split into two independently attackable quantities: the *reach of two methods*, and the *physical mass* `m/Λ`. Section 3 states the two reach hypotheses; Sections 4–6 explain why they are hard to obtain at `q ≈ 1`.

### 2.4 Two-loop closed form and numbers

The exact solution of the truncated two-loop flow is

```
t_up^{2L}(λ) = C_RG/λ − (51/121) ln(1 + 121 C_RG/(51 λ)) .
```

For `SU(3)` at `λ = 6` (`g² = 2`): `t_up^{2L} = 2.639`; with `m/Λ = 8` this gives `q₀ = 0.571`; `q₀ = 1` at `λ ≈ 7.28` (`g² ≈ 2.43`). One loop gives `3.589` (`q₀ = 0.221`); three-loop/numerical values give about `2.36/2.0` (scheme-dependent and **not** used as rigorous input).

**Consequence: `CPB-W` must supply at least a two-loop rigorous bound.** With one-loop control only, `q₁` cannot exceed `0.22` and the handoff cannot close.

---

## 3. Two-sided crossover preservation

### 3.1 `CPB-W` (weak side)

**Conjecture.** There exist `λ_c ≈ 7.3` (two-loop; `≈ 12` in the conservative one-loop version) and a weak-side class `C_W(ε₀)` — a massive Gaussian reference of mass `m` plus a bounded interaction, with relative error `ε₀ ≤ 1/2` — such that the effective action stays in `C_W` for all `λ ≤ λ_c`, and `C_W` is invariant under one blocking step; and such that a rigorous two-loop bound `RB2` on `t_up(λ_c)` holds.

The precise content of `RB2` is

```
|β(λ) − β_{≤2}(λ)| ≤ C₃ λ⁴ ,      0 < λ ≤ 7.3 ,
```

with `C₃` explicit and uniform in the UV cutoff and the volume, or the one-sided sign `β ≤ −Aλ² − Bλ³`.

**Status: open.** Its quantitative core `NPB` is a non-perturbative beta-function bound; rigorous RG of Balaban type stops at `λ ~ 1` (Section 4.3).

### 3.2 `CPB-S` (strong side)

Write the effective action as `S_k = S_k^G + V_k` (massive Gaussian plus interaction) and define the **crossover activity**

```
z_k(q) := sup_x Σ_{C∋x} |w_q(C)| · e^{+κ q |C|} ,     κ = 1/2 ,
```

where `w_q(C)` is the signed weight of a connected cluster `C` and `q = m a_k` is the dimensionless mass in block units. **Conjecture.** There are absolute constants `q_c ∈ [1, 3/2]` and `z_c = 1/2` such that `q_k ≥ q_c` implies `z_k(q_k) ≤ z_c`, and the class is preserved under one blocking step.

**Key structure:** `q_{k+1} = 2 q_k` is exact kinematics, and `z` is non-increasing in `q` (provided `c₀ > κ`); hence the worst point is the smallest `q` in the window, namely `q_c`. Therefore `CPB-S` **only needs to be proved once, at `q_c`**.

### 3.3 Why the standard routes fail at the crossover

- **Balaban.** The small-field expansion is around the trivial connection, with the coupling as expansion parameter; by construction it stops at `λ ~ 1`. At `λ ~ 6` the large-field region is not suppressed by `β`.
- **Strong-coupling expansion.** Its radius is governed by `β/(2 dim f)`; `q₁ ~ 1` needs `β ~ 2.2`, far beyond the thresholds of rigorous cluster expansions.
- **Bakry–Émery.** The Wilson Hessian contains `−[F,F]`-type negative directions; only a mass term can control them, and the mass is unknown.

---

## 4. Proof skeleton and the load-bearing wall

### 4.1 The ten-lemma skeleton

The class-preservation part of `CPB-W` decomposes as follows (`[B]` = Balaban/standard, `[new]`, `[open]`):

```
1  Blocking and covariance splitting      [B]        G = B^T G' B + G_⊥
2  Gaussian conditional expectation/Wick  standard
3  Cumulant expansion                     [B]+[new]  remainder O(ε³L^{3/2}), needs εL^{1/2} ≤ 1/2
4  Fluctuation sum (one-loop log)         [new]      L(q) ≤ C₃ ln(1/q) + C₄
5  Lipschitz contraction                  [B]+[new]  sup does not contract; C¹/Lip does
6  Large-field Lipschitz/C² bound         [B] small g² only; q~1 [open]   ← the wall
7  Gauge independence / equivariant blocks [B]+[new]
8  Volume/UV uniformity                   [B] partial
9  One-step class preservation            [new]      conditional in U(1)
10 Iteration + RB2/NPB                    [open]
```

### 4.2 Lemma 6: the large-field bound fails at `q ≈ 1`

Controlling `V` in the `C²` norm (sup + Lip + Hessian) on the ball `‖A‖ ≤ R`, non-degeneracy requires

```
S'' ≥ (c₀ m² − Cε)·Id > 0 ,     ε < c₀ m² / C .
```

The non-abelian commutator gives `‖V''‖ ≲ R³/g² + R²/g`; at `q ≈ 1`, `g² ≈ 2`, and typical field amplitude `R ~ g/q ≈ 1.4`, this is `≈ 2.8`, **exceeding the admissible `c₀ m²/C`**. Balaban's large-field estimates are of the same class but are proved only for small `g` (UV); no `q ~ 1` version exists.

**This is the core of `CPB-W`'s class preservation; it is not `RB2/NPB`.**

### 4.3 The Gaussian-weighted class and the cluster constant

A natural remedy is the local `L²` class

```
C_W^G(ε) = { S = S^G + V : sup_x ∫ |V_x|² dμ_G ≤ ε² } ,
```

with `μ_G` the massive Gaussian of covariance `g²(−Δ+q²)^{-1}`. It moves the expansion parameter from `ε ≈ 2` (a definite failure) to

```
εL^{1/2} = 0.24 + C_cum·R₃ ,      R₃ = 0.16–0.46 ,
```

which passes only if `εL^{1/2} ≤ 1/2`, i.e. `C_cum ≲ 0.6–1.6`. The standard Kotecký–Preiss bound at `g² ~ 2` gives **`C_cum ≥ 2`**:

```
C_cum = 1 :  0.40–0.70  (passes only at the low end)
C_cum = 2 :  0.56–1.16  (fails)
C_cum = 3 :  0.72–1.62  (fails)
```

**Hence Gaussian weighting does not bypass Lemma 6; it moves the wall to `C_cum`, and the standard constant is already out.**

---

## 5. The Bianchi mechanism and its failure

A recurrent hope is that the Bianchi identity `dF = 0` produces sign cancellations that suppress the activity. We settle this point by point.

### 5.1 BL-1: `dF = 0` is a topological statement

`dF = 0` is exactly the flux 2-cocycle condition (the product of plaquettes around every cube is the identity). Only in the abelian case with quantized flux does it imply that defect 2-chains are closed. It constrains flux conservation; it does **not** constrain the sign or the size of single-block weights.

### 5.2 BL-3: a closed 2-chain can be codimension one

The slice `x₄ = 0` is a closed 2-chain, but it is a **codimension-one slab**, whose per-block cost is independent of `q`. Hence "defects are closed 2-surfaces" does **not** imply an area-type action `σ(q) ~ q²`.

### 5.3 BL-4: the improvement is in entropy, not in action density

`dF = 0` acts on the **counting rate** of connected 2-chains, `γ_all → γ_cyc`, a **`q`-independent pure number**; it does not act on `|w_q(C)|`. In the additive form `|w| ≤ exp(−c₀q|C| − c₁q²|C|)`, one can prove **`c₁ = 0`**.

### 5.4 BL-5: quadratic scaling is empty at `q ≈ 1`

At `q = 1`, any `q^α` with `α > 1` equals `q`. **Any "quadratic scaling" improvement has exactly zero gain at the point where it is needed.** The threshold is the absolute efficiency `σ(1)` at `q = 1`, not the scaling.

### 5.5 BL-6: `κ_c = 1` in the strong-coupling regime

The only mechanism that could serve both sides is **signed cancellation at fixed boundary**,

```
|Σ_{∂C=∂C₀} w_q(C)| ≤ e^{−c₀q|C₀|} · κ_c^{|C₀|} ,      κ_c < 1 ,
```

where `κ_c` is the per-block cancellation factor (`κ_c = e^{−δ}`). `CPB-W` needs `κ_c ≤ 0.80`; `CPB-S` needs more.

**Verdict (BL-6, a judgment rather than a theorem): the default is `κ_c ≈ 1`, and no known mechanism gives `κ_c < 1`.** Three converging lines of evidence:

1. **In the strong-coupling surface/character representation, `κ_c = 1` is provable exactly:** the `SU(2)` coefficients `2I_{2r+1}(β)/β > 0` are all positive, so all surface-sum terms are non-negative and the signed sum equals the absolute sum.
2. **Weights factorize locally**, so signs do not alternate with area; `κ_c < 1` would require a per-block `e^{iπ · additive charge}`.
3. **`dF = 0` supplies no additive defect charge:** 2D `O(N)` uses exact constraints plus a topological index, 4D `U(1)` uses duality/Poisson, and **4D non-abelian pure Yang–Mills has neither** (the instanton density is not a defect charge; `θ` only produces a global phase).

`U(1)` provides a clean counterexample: `F = dA` makes `dF = 0` automatic, yet there is no non-trivial signed cancellation. **Thus `dF = 0` and `κ_c < 1` are logically independent.**

> **This is a negative result about routes, not a disproof of `CPB-W`/`CPB-S`.** What is ruled out is the route that proves them via Bianchi / signed cancellation.

---

## 6. Numerical tests

### 6.1 Apparatus and validation

`SU(3)` Wilson action, Cabibbo–Marinari heat bath (exact `SU(2)` subgroup sampler), cross-validated against Metropolis, strong-coupling series, and the known `β = 6` plaquette. Lattices `L = 4, 6, 8, 10`; primary grid `β = 2.6, 2.8, 3.0, 3.3` (`g² = 1.8–2.3`), 8000 measurements each.

### 6.2 `R₄` and `κ_eff`

`R₄` (the connected four-point ratio of the blocked plaquette) satisfies `|R₄| ≲ 0.13` for all `β`, `L`, `b` (`≲ 0.01` for `b ≥ 2`); no droplet structure is seen.

`κ_eff` (lowest eigenvalue of the projected Hessian) is negative throughout `β = 1.7 … 6`, **and remains `−1.74` at the weak-coupling control `β = 40`**. Its scaling `|κ_eff| ≈ (0.85–1.6)·β(1−⟨P⟩)` identifies it as a Savvidy/Nielsen–Olesen negative mode of the **microscopic** Wilson action — present for any non-convex action, including continuum Yang–Mills. **It is therefore not a discriminator for `CPB-W`**; it is reported but not used.

### 6.3 All character coefficients are positive

Character coefficients of the plaquette distribution, `c_R = ⟨χ_R(U_p)⟩`:

| `β` | `c₃` | `c₆` | `c₈` | `c₁₀` | `r₁` |
|---|---|---|---|---|---|
| 2.6 | +0.52107(19) | +0.13495(18) | +0.23809(26) | +0.02323(16) | 1.0000 |
| 3.0 | +0.61476(19) | +0.18739(20) | +0.32620(27) | +0.03800(17) | 1.0000 |
| 3.3 | +0.68861(19) | +0.23549(20) | +0.40444(29) | +0.05317(18) | 1.0000 |

`L = 6` and `L = 8` agree. **No `c_R < 0` occurs**, and the single-block signed/absolute ratio is `r₁ = 1.0000`: there is no sign cancellation at the `c₁` level. This is consistent with BL-6's `κ_c ≈ 1` (a **confirmation**, not a failure), but `r₁` is a `c₁` proxy, not a measurement of `κ_c`.

### 6.4 Defect fractions are Gaussian tails

A seemingly powerful proxy estimates the defect action as `c₀ = −ln(fraction)` from threshold defect fractions. **This proxy is falsified**: the blocked field is itself **Gaussian**, so every threshold fraction is merely a Gaussian tail.

```
β=2.6, b=2, L=8:  μ=0.17355, σ=0.02805
Gaussian-tail prediction:  f(⟨P⟩/2)=0.00099 ,  f(0.75⟨P⟩)=0.0610
Measured:                  f(⟨P⟩/2)=0.00089 ,  f(0.75⟨P⟩)=0.0602
```

The variances at `L = 6` and `L = 8` agree to `<1%` (`0.000781` vs `0.000787`), so this is not a finite-volume effect. Changing the threshold or the coupling moves `c₀` between `2.4` and `10.1`: **it is a function of the threshold and the distribution width, not an independent action.**

### 6.5 `ξ/a ≈ 7`: feasibility, and the evidence/certificate distinction

The σ certificate requires `q₁^σ ≤ 0.143`, i.e. `ξ/a ≈ 7`. Measured: `ξ/a(β=6) = 0.93(8)`; at `β ≤ 3.3`, `ξ/a ≲ 0.3` (the correlator dies within one lattice spacing). Reaching `ξ/a = 7` needs `L ≈ 28–34`.

**Numerics do not constitute the certificate.** `q₁^σ ≤ 0.143` is a rigorous inequality (a theorem or an input); neither a point value nor a confidence interval can prove it. Numerics can only supply feasibility/consistency **evidence**, and if used as evidence the correct statistical object is a **2σ lower bound**: at 5% precision one must measure `ξ/a ≈ 8.5` (since `8.5 − 2×0.425 = 7.65 > 7`), and one needs `L/ξ ≥ 4`. At `ξ ≈ 7`, `L = 12` gives `L/ξ = 1.7`, so an `L = 12` scan determines `ξ/a` only up to `≈ 3–4` and **cannot support `ξ/a ≥ 7`**.

**Cost.** With pure NumPy, cost scales as `L⁴`: at `L = 32`, about `38 s/measurement` (about `105 h` for `10⁴` measurements). A C/OpenMP kernel (16 cores, ≥50×) plus an anisotropic volume reduces a `32⁴` run at 5% precision to about **74 h** with a GEVP variational operator (3 days on 16 independent chains). This is a feasibility limit, not numerical noise.

---

## 7. Remaining paths

Given the verdicts above, `CPB-W`/`CPB-S` have no known route under existing mechanisms. The remaining candidates are:

1. **Inverse Monte Carlo (IMC).** Fit the couplings of the blocked effective action to obtain `w_q(C)`, hence `c₀` and `κ_c` directly. This is the **only remaining numerical handle**, and it is a project of its own.
2. **A new additive sign/index structure.** The object that `κ_c < 1` requires. 2D `O(N)` and 4D `U(1)` each have their own version; whether 4D non-abelian admits one is unknown.
3. **A different handoff architecture.** A different weak/strong matching point, a different RG scheme, or abandoning the handoff for an all-coupling construction.
4. **A different certificate.** The string tension `σ` instead of the mass `m`, or an observable with a cleaner signal; whether the `ξ/a ≈ 7` threshold can be avoided is unknown.

---

## 8. The conditional theorem and machine verification

The reduction above is embedded in a larger **conditional theorem**: ten named inputs `IN-0 … IN-9` (the certified inputs of the 2026 construction, the UV limit, `O(4)`, `CPB-W`, `CPB-S`, the overlap, IR uniqueness, non-triviality, group uniformity, and spectral-resolution convergence) imply the Clay statement.

Its "nine-to-one" part consists of nine glue lemmas `GL1–GL9`, of which:

```
Unconditional ✅ (6 blocks): GL1 (audit) · GL2a · GL4 · GL6 · GL5-KIN · GL8-A
✅-implication (6 blocks)  : GL3 · GL2b′ · GL5-DYN · GL8-S · GL9 · GL7
No formal content          : GL5-DYN (IN-2 ⟹ OS2 is bookkeeping)
Audit (not a theorem)      : GL1 (cross-document compliance)
```

**Machine-verification status:** 9 of the 11 blocks are formalized in Lean 4, with no `sorry`, and the axiom footprint is frozen to `[propext, Classical.choice, Quot.sound]`. Four items are **explicit declarations rather than machine proofs**, all pushed into the **type level**:

1. Minlos (the `hMinlos` hypothesis of GL3; absent from Mathlib);
2. Borel measurability of the truncation (an explicit hypothesis in GL4; Mathlib lacks `comp_stronglyMeasurable`);
3. the identification of the variational gap with the spectral gap (a declarative scope in GL2b′);
4. Reed–Simon VIII.24b (GL2b′; absent from Mathlib).

**These four are external assumptions, not proofs.** They appear in theorem signatures or file headers precisely so that a reader can see them, rather than having them hidden in prose.

---

## 9. Conclusion

We have reduced the mass-gap problem for four-dimensional pure Yang–Mills to a handoff inequality, and then shown that every known mechanism capable of suppressing the single-block activity at `q ≈ 1` fails at exactly the point where it is needed.

- Bianchi is topological, acting on entropy rather than weights; quadratic scaling is empty at `q = 1`; and `κ_c = 1` is provable in the strong-coupling regime.
- Gaussian weighting moves the wall to the cluster constant `C_cum`, where the standard bound is already out.
- Numerically, `κ_c` is not directly measurable, the σ defect-fraction proxy is a Gaussian tail, and `ξ/a ≈ 7` is a feasibility boundary.

**This is not a failure report but a map:** it records which routes have been walked, why each is closed, and under what conditions it would reopen (a measured `κ_c < 1` at `q ≈ 1`, or a new additive sign/index structure).

---

## Appendix A. Constants and symbols

| Symbol | Meaning | Do not confuse with |
|---|---|---|
| `λ_G` | `h^∨ g²` (RT invariant) | not `dim f · g²` |
| `b₀` | `11h^∨/(48π²)` | not `11 dim f/(48π²)` (they coincide only for `SU(N)`) |
| `B` | `2b₀ = 11h^∨/(24π²)` | not `b₀` |
| `γ` | `β/N_c = 2/g²` (measure weight) | a different variable from `λ_G` |
| `κ` | `1/2` (Kotecký–Preiss factor) | not `κ_c` |
| `κ_c` | signed-cancellation factor `<1` | not `κ` |
| `σ(q)` | defect action density | not the string tension `σ` |

## Appendix B. Reproduction

```
bash lean/check-handoff.sh            # handoff algebra, axioms frozen
bash lean/check-conditional.sh        # conditional-theorem skeleton
bash research/glue-a-limits/check-glue-a.sh
bash research/glue-b-blocking/check.sh
bash research/glue-c-gauge-group/check-glue-c.sh
bash research/glue-d-dual-handoff/check-gl7.sh
bash research/os2-rotation/check-symdim4.sh
python research/os2-rotation/check-d14.py
```

Numerics: `research/lattice-cpb/` (heat bath, analysis, `t_up.csv`, `block_summary.txt`).

## References

1. K. Wilson, *Confinement of quarks*, Phys. Rev. D 10 (1974).
2. K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. 110 (1978).
3. T. Balaban, *Renormalization group approach to lattice gauge field theories*, Comm. Math. Phys. 102 (1985), 109 (1986), 122 (1989). *(volume numbers to be verified)*
4. R. Kotecký, D. Preiss, *Cluster expansion for abstract polymer models*, Comm. Math. Phys. 103 (1986).
5. S. Cao, M. Park, S. Sheffield, *Random surfaces and lattice Yang–Mills* (2024). *(to be verified)*
6. A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Stochastic quantisation of Yang–Mills–Higgs in 3D*, Invent. Math. 237 (2024).
7. The 2026 four-dimensional pure Yang–Mills construction (technical report, assessed by an independent committee; unpublished).
8. The hexagonmath reconstruction `hexagon:2611.00492v1` (a site identifier, **not** an arXiv number; an AI-assisted preprint).
9. Constructive results for 2D `O(N)`. *(exact source to be verified)*

---

## Acknowledgments

This work was carried out with substantial AI assistance. **Deepseek v4.1 flash participated deeply in the research process**, including literature synthesis, derivation checks, and numerical analysis; AI assistants were also used for drafting and code generation. All mathematical claims, and any errors, remain the responsibility of the author.

---

> **Disclaimer.** This paper does not claim a proof of Yang–Mills existence or of the mass gap. Every "open" input is labelled as such in the text; every negative result is accompanied by a falsification test. Where a judgment (rather than a theorem) is asserted, it is marked as a judgment. Readers who find an error should use the reproducible commands and numerical data referenced above.

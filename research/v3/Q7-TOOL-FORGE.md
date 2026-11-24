# Q7-TOOL-FORGE — searching for a **sixth species** for H1: the census of seven candidate families
# collapses into two sharp, computed no-gos (a *direction* theorem and a *frozen-outlier/fine-tuning*
# theorem); both leave exactly one named escape hatch, and neither is a live wire to the gap

> 2026-11-24 · **H1** (the last open input): `χ = Var_μ(ξ_{k→0}) ≤ 1/λ` uniformly in `(L,a)` at
> `β=2.5`, equivalently `G″_{k=0} ≥ λ>0` on the coarse marginal. Five paradigms
> (positivity / convexity / spectral / Perron / RG) are exhausted; the charter (§1) demands a
> **sixth species**, not a repackaging. This file proposes and *first-tests* seven candidate
> families, then drills into the two that survive the census.
> Reads: `docs/new-tools-registry.md` (A1–A5, B1–B8), `P-WAVE-SYNTHESIS.md`, `STAR2-CLASSIFICATION.md`,
> `H6-HEAD-SPECTRAL.md`, `unconditional/PROOF-JACOBIAN-SAVES.md`, `SM2-RANK-ONE.md`, `SM5-TCI-HEAD.md`,
> `docs/paper-v3/main.tex` §6.
> Witnesses: `research/v3/w7/q7_direction.py` (Q1), `research/v3/w7/q7_rmt.py` (Q2); `py -3.12 -I`, numpy only.
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step / open · **[V]** numerical witness here · **[std]** standard
> textbook tool used as such. **No QED. Clay 0%.**

---

## 0 · Verdict (no QED)

**The seven candidate families do not contain a live sixth species that proves H1.  Five of the seven
collapse — provably — onto objects the campaign has already banked (`Γ₂`, `H1`-by-definition, the
frozen Savvidy edge, a phase constraint, or the spectral no-go).  The two that do not collapse reduce to
two sharp, *computed* no-gos, each with exactly one escape hatch:**

1. **[Q1] Direction theorem `[T]`/`[NG]`.** *Analytic-duality species* (Fourier/uncertainty +
   Fisher–Rao + Paley–Wiener/zero-set).  **Unconditional theorem:** `E_ρ[G″] = I(ρ) ≥ 1/Var_ρ(a)`
   (Cramér–Rao; sharp for the Gaussian).  This is a **lower bound on average curvature = a lower bound
   on the variance** — the *logical opposite* of H1 (`Var ≤ 1/λ`).  The entire class of
   "shape/width/analyticity/information" functionals of the marginal is structurally a **variance-floor
   (curvature-ceiling)** machine: Bernstein gives `G″ ≤ τ²` (a *ceiling*), never a floor, and a
   band-limited density `ρ_ε = 1+ε cos(τa)` has `G″ = ετ²/(1+ε) → 0` at *fixed bandwidth* as `ε→0`
   (the free/flat field).  **Verified** (witness Q1).  This closes the analytic half of Picard–Lefschetz
   (the zero-set/strip shadow) and all of Fisher–Rao.

2. **[Q2] Frozen-outlier / fine-tuning theorem `[T]`/`[NG]`.** *Random-matrix / free-probability
   species.*  **Two exact facts:** (a) the Savvidy soft mode is an **exact deterministic eigenvalue**
   `−g|B|` below the positive Landau edge `+g|B|` (separation `2g|B|`) — a fixed outlier, **not** a
   BBP/BaW fluctuation, so the frozen-Hessian RMT *returns the negativity* `G″<0` and is structurally
   **blind to the measure rescue** `D_+ ≤ 0`; (b) on the measure-averaged object `G″ = S − D` with
   measured `D/S ≈ 0.97`, the gap is the thin residue `G″/D ≈ 3%`, while every available RMT/free-prob
   *relative* error on `D` is `≥10%` — so the tool **cannot resolve `sign(G″)`** at the head.
   **Verified** (witness Q2).

**The one escape hatch of each.**  Q1: a *converse-direction* inequality `Var ≤ 1/λ` that does **not**
come from a shape/width duality — which is H1 itself, and which for a non-log-concave measure is exactly
the missing uniform strong-concavity.  Q2: *relative*-error control `< G″/D ≈ 3%` on the force variance
`D` — i.e. the **variance of the variance** / a sub-leading free cumulant, which no leading-order RMT
result supplies.  Both escape hatches are, word for word, H1.

**Clay = 0%.  The wall did not move.**  What is new is a *map*: the sixth-species search is now
**classified**, and the reason it fails is a single one — every candidate is a **duality/fluctuation
functional of `ρ`** (or of the frozen operator), and every such functional is a *variance floor*
(wrong direction) or *leading-order blind to a 3% fine-tuning*.  The gap lives in the *shape region*
`{ρ bounded away from uniform}` of the marginal, whose quantitative strength **is** `ess inf G″ ≥ λ`.

---

## 1 · What the sixth species must beat (the bar)

The five exhausted paradigms all answer one of two questions about the **scalar coarse potential**
`V_eff(a) = −log ρ(a)` (`ν` = coarse marginal, `a ≡ ξ_{k→0}`):

| paradigm | object | why it dies (banked) |
|---|---|---|
| positivity | `ν ≥ 0`, KL/RP, complete monotonicity | sign only; **B6**: `r(t)` is an *upper* bound on `m₁²`, direction reversed |
| convexity | Brascamp–Lieb, `Φ` convex, Prékopa, `W_geo ≥ 0` | sign of `W_geo` is `[T]` but dominance over Savvidy = `Var(G)` bound = H1 |
| spectral | PF/Krein–Rutman/Bushell, Bakry–Émery, Cheeger, LSI, block transfer | **H6 `[NG]`**: `λ₁(ν) ≥ ess inf G″`, but `ess inf` is negative (Savvidy tail) |
| Perron | gauge-invariant PF kernel | pins the vacuum (same `0⁺⁺` channel), no modulus |
| RG | one-step descent, Feshbach/Schur | transports a floor `4π²/L² → 0`; fails on the soft mode (`η_j ≥ 1`) |

**The single shared obstruction** (`H6 §3`): **typical-point positivity `G″(a₀*) > 0` vs essential
infimum `ess inf G″ < 0`**, the Savvidy double well.  And (B6/B8-IV): every *moment / Laplace / sum-rule*
readout is an **upper** bound on `m₁²` or `G″` whenever the source is a positive density.

> **Bar.** A genuine sixth species must produce `G″_{k=0} ≥ λ > 0` (or `Var ≤ 1/λ`) **uniformly**,
> using a functional of `ν` that is *neither* (a) the pointwise sign of `V_eff''` / its essential
> infimum, *nor* (b) a single-operator spectral gap, *nor* (c) a Perron eigenstructure, *nor* (d) an
> RG contraction.  It must also **fail the free-field acid test** — i.e. return `0` on the massless
> field *non-vacuously*.

---

## 2 · The census of seven candidate families (+ extras)

Columns: **mechanism** (what the species actually measures); **what it would prove**; **free-field
acid test** (must return `0`/`≤0`); **verdict** and *where it collapses*.

| # | species | mechanism | would prove | free-field acid test | verdict |
|---|---|---|---|---|---|
| **(i)** | **OT / Otto / Ricci** on the gauge group | entropic gradient flow in `W₂`; synthetic `Ric_ν ≥ ρ ⇒` contraction `λ₁(ν) ≥ ρ` | `λ₁(ν) ≥ ρ` | free Gaussian: `W₂`-rate `=k²→0` ✓ | **[NG]** — in 1-D `RCD(ρ,∞) ⇔ V_eff'' ≥ ρ` a.e. (1-D characterization of synthetic lower Ricci); so "synthetic" = classical `Γ₂`, banked **negative** (A5: `Γ₂ ≈ −2.0`; W4: negative mode survives gauge quotient, 7.1%). *Escape:* a `W₂`-contraction **not** via `Γ₂` — forbidden by the 1-D collapse. |
| **(ii)** | **Picard–Lefschetz thimbles** / positive decomposition | complexify lattice field; real cycle = Σ thimbles `w_σ G″_σ`, `w_σ` = Lefschetz measures `≥0` | weighted-avg `G″ ≥ λ` if dominant thimble log-concave | free Gaussian: single real thimble, `G″=k²→0` ✓ | **[O]** (as saddle-sum) + **[NG]** (as analyticity). Real weight is positive on the real contour (no sign problem to cure); the Savvidy real saddle is a Lefschetz *saddle* — the real cycle crosses a Stokes line, so it is **not** one thimble; controlling `Σ w_σ` *is* the measure control H1 needs. Its *analytic shadow* (zero-set/strip) is **Q1** → ceiling only. |
| **(iii)** | **Information geometry** (Fisher–Rao, Wasserstein-RG) | `E_ρ[G″]=I(ρ)` exact; Cramér–Rao `I·Var ≥ 1`; RG as `W₂` flow | `Var ≥ 1/I` (**a floor**) | free flat `ρ`: `I=0, Var=∞` ✓ | **[NG]** — Cramér–Rao is a *variance floor* (= curvature ceiling), wrong direction; `W₂`-RG rate `=λ₁=G″` (1-D collapse). This is **Q1**'s Fisher half. *Escape:* a de-Bruijn identity giving variance *decay* along RG — absent for non-log-concave measures. |
| **(iv)** | **Random matrix / free probability** | spectral edge + outliers; free cumulants for `D=Var(force)` | distribution of `λ_min`; `D` to leading cumulant | free field: no spike, edge `0` ✓ | **[NG]** — **Q2**: soft mode is an *exact* deterministic outlier `−gB` (frozen RMT → `G″<0`, the **negativity**; blind to `D_+`); on `S−D`, RMT rel-error `≥10% ≫ 3%` gap. *Escape:* `<3%` relative control of `D`. |
| **(v)** | **Tensor network / MPS** transfer, rigorous gap certificate | MPS ground state + tangent-space / injective-MPS gap bound | rigorous gap lower bound *if* volume-uniform certificate exists | gapless free MPS: vanishing gap ✓ | **[NG]** — subsumed by **H6** (block transfer); no rigorous volume-uniform MPS gap theorem exists in `d≥2` (all rigorous MPS gap results are 1-D); the coarse transfer gap *is* `H1` by H6-B. *Escape:* a rigorous 2-D/3-D gap theorem — open in its own right. |
| **(vi)** | **Topological / θ-vacuum / 't Hooft anomaly** | anomaly matching + 1-form center: gapped IR ⇒ TQFT or center SSB; θ-periodicity | *existence* of a gap / phase, **not** a magnitude | free `U(1)`: trivial θ, no anomaly → no constraint ✓ | **[NG]** for H1 + **[T]** for the phase. It is the **wrong logical type**: a yes/no constraint on massless states, not a *uniform lower bound* `λ`. (Consistent with P3's unbroken center at θ=0.) *Escape:* none — phase ≠ magnitude. |
| **(vii)** | **others** | see below | — | — | all **[NG]** |

**§(vii) extras.**
- **Kubo / KMS / f-sum rule.** `D = ∫⟨ff⟩_c dt`; gap = threshold of the spectral function; the f-sum
  rule bounds the *weight below the threshold*, not the threshold. **[NG]** = spectral/positivity
  rename (cf. B8-IV: superconvergence gives only upper bounds).
- **Modular / algebraic QFT (Tomita–Takesaki).** Mass gap = bottom of the modular Hamiltonian
  spectrum; but KMS = transfer matrix, and it is not lattice-computable at `β=2.5`. **[NG]** formal.
- **Lee–Yang / complex-parameter zeros.** Gap = distance of the nearest partition-function zero to the
  real axis; zero-free regions ⇐ cluster expansion ⇐ small `β`. **[NG]** = RG/cluster rename; no
  zero-free region at `β=2.5` (that *is* the gap).
- **Integrability / Bethe ansatz.** Solvable points are `β=0` (Haar) and `β=∞` (free) — both
  **off-tower**; finite-`β` lattice YM is not integrable. **[NG]** = dual of the rigidity no-go **B3**
  (the gap is a `β`-non-analytic jump; solvable points return the free value).
- **Holographic / AdS-CFT (soft wall).** `m ~ Λ`; a model, not a rigorous tool. **[O]** (not a tool).
- **Mori–Zwanzig / projection operators.** Effective mass from fast-mode elimination; the sign is
  convention (B) Feshbach (`Σ ⪰ A`), already banked. **[NG]** = Feshbach species (SM5 §1.2).

> **Census result.** Of the seven task-named families, (i),(iii),(iv),(v),(vi) are **[NG]** and
> provably collapse onto banked objects; (ii) splits into an [O] saddle-sum (circular) and an [NG]
> analyticity shadow.  The only *unbanked* content is the two directional/fluctuation theorems below.
> **There is no seventh live wire.**

---

## 3 · [Q1] Direction theorem — analytic duality is a **variance floor** `[T]`/`[NG]`

### 3.1 Statement and proof

Let `ρ` be the coarse marginal density of `a ≡ ξ_{k→0}` on its support (an interval; for the physical
marginal `ρ → 0` at the chart edge `|a| = 2π`, so boundary terms drop).  Define
`G″ = −(log ρ)″ = −ρ″/ρ + (log ρ)′²` (the `V_eff″` object of H6; equals `−ρ″/ρ` at a max `a₀*`).

> **Theorem Q1 [T].** (a) *Exact:* `E_ρ[G″] = I(ρ) := E_ρ[(∂_a log ρ)²]` (integration by parts).
> (b) *Cramér–Rao:* `I(ρ)·Var_ρ(a) ≥ 1`, i.e. **`E_ρ[G″] ≥ 1/Var_ρ(a)`**, with equality iff Gaussian.
> (c) *Paley–Wiener/Bernstein:* if `ρ` has exponential type `τ` (band-limited to `[−τ,τ]`), then
> `G″ = −ρ″/ρ ≤ τ²` pointwise — a **ceiling**.

**Proof.** (a) `E[−(log ρ)″] = −∫ρ(ρ″/ρ − (log ρ)′²) = −∫ρ″ + ∫(ρ′)²/ρ = 0 + I`.  (b) Cauchy–Schwarz:
`Var(a)·I = ∫(a−ā)²ρ ·∫(ρ′/ρ)²ρ ≥ (∫(a−ā)ρ′ da)² = (∫ρ da)² = 1` (boundary term `[(a−ā)ρ]` vanishes).
(c) Bernstein `|ρ′|≤τρ` ⇒ `|ρ″|≤τ²ρ`. ∎

**Why this is a no-go for H1.**  H1 is `Var ≤ 1/λ` (equivalently `G″ ≥ λ`, *uniformly*, and at the
essential infimum, not just in average).  Theorem Q1 gives the **reverse**: `Var ≥ 1/E[G″]` and
`G″ ≤ τ²`.  It (i) *permits* arbitrarily large `Var` — a spread marginal satisfies it with a weak lower
bound `E[G″] ≥ 1/Var ≪ λ` — and (ii) caps `G″` at `τ²`.  Neither can *forbid* the gapless scenario,
and (c) shows the *only* sign the analytic data (bandwidth `τ`) controls is the **upper** one.  The
free field (`ρ` flat) is the extremal case: `I=0`, `Var=∞`, `G″=0` — **saturating (b) and (c) from the
gapless side**.

### 3.2 The exact-anchor reading (why P1/P2 are ceilings)

`P-WAVE §2` proved the exact global-holonomy curvature is `−ρ″/ρ = 1/(2·Cin(2π)) = 0.2051`, and the
periodic `k=0` value `c = 0.222`.  Q1 says both are **values strictly below the Bernstein ceiling**
`τ²` of the exact marginal and are on the "peak" side of the duality.  That is *why* the corpus keeps
finding them as **floors of a tractable model** but they never transport to 4D: the object `−ρ″/ρ` is,
by Theorem Q1, a quantity whose *analytic data* bounds it only from above.  **The anchor is a ceiling
statement, not a floor statement** — the logical status that `P-WAVE §2,§5` already flagged
(`−ρ″/ρ` vs `1/Var`) is now explained structurally: every duality functional is one-sided the wrong way.

### 3.3 Witness (`w7/q7_direction.py`, quoted)

```
(A) BERNSTEIN CEILING vs the eps->0 (free-field) limit  [bandwidth tau fixed]
   tau =  2.0
      eps= 0.01:  -rho''/rho(peak)=   0.03960  grid max=   0.03960  (ceiling tau^2=4.0)
      eps= 0.50:  -rho''/rho(peak)=   1.33333  grid max=   1.33333  (ceiling tau^2=4.0)
      eps= 0.99:  -rho''/rho(peak)=   1.98995  grid max=   1.98995  (ceiling tau^2=4.0)
      => eps->0  gives  G''->0  at tau=2;  the ceiling is NEVER a floor.
(B) CRAMER-RAO:  I(rho)*Var(a) >= 1  is a VARIANCE FLOOR (wrong direction for H1)
   sigma=  1.0:  G''(peak)= 1.000000   I= 1.000000   Var=     1.000   I*Var=  1.0000
   sigma= 10.0:  G''(peak)= 0.010000   I= 0.010000   Var=   100.000   I*Var=  1.0000
   sigma= 40.0:  G''(peak)= 0.000625   I= 0.000625   Var=  1600.000   I*Var=  1.0000   (free: Var->inf, G''->0)
      two-Gaussian mixture: I= 0.89479  Var= 2.20250  I*Var= 1.9708 (>1)
(C) EXACT FISHER IDENTITY:  E_rho[-(log rho)''] = I(rho)
   sigma=1.0:  E[-(log rho)'']= 1.00000   I= 1.00000   |diff|= 7.31e-14
   sigma=5.0:  E[-(log rho)'']= 0.04000   I= 0.04000   |diff|= 3.43e-15
```
Readings: (A) fixed bandwidth, `G″→0` as the marginal flattens — **the ceiling `τ²` is never a floor**;
(B) every entry satisfies `Var ≥ 1/I`, none certifies `Var ≤ 1/λ`; the free field is the `σ→∞` row;
(C) the exact identity `E[G″]=I` to `1e-14`.

### 3.4 Exact missing step (Q1)

> **`[O]` (Q1).** An *upper* bound `Var_ν(a) ≤ 1/λ`, equivalently a **lower** bound on the curvature
> `ess inf G″ ≥ λ`, obtained from a functional that is **not** a shape/width duality.  By Theorem Q1
> every duality functional (uncertainty, band-limit, Fisher, Cramér–Rao, zero-set/strip) is a *variance
> floor*; the needed ingredient must certify the marginal is **bounded away from uniform** with a
> *quantitative modulus* — which is, verbatim, the uniform strong-concavity H1.

**The positive residue `[T]`.** Theorem Q1(b) is unconditional and gives `Var_ν(ξ_{k→0}) ≥ 1/I(ρ)` — a
*lower* bound on the susceptibility.  It says the coarse mode is **not over-localized**; combined with
the free field (`Var→∞`) and the periodic plateau (`Var≈4.45`, P-wave) it is consistent.  It is a true
statement in the *wrong* direction — precisely the content of the no-go.

---

## 4 · [Q2] Frozen-outlier / fine-tuning theorem — RMT is blind to the measure rescue `[T]`/`[NG]`

### 4.1 Two exact facts

> **Fact Q2-a [T].** The Savvidy soft mode of `Q = −D² + 2\,ad F_B` is an **exact, deterministic**
> eigenvalue `λ_- = −g|B|` (`SM5` eq.(1.1)); the positive rungs sit at `~(2n+1)g|B|`.  The separation
> is `2g|B|`.  This is a **fixed outlier** of the operator given the background `B` — **not** a Baik–Ben
> Arous–Péché random fluctuation of a Wigner bulk (the tower is determined by `B`, not sampled).

> **Fact Q2-b [T]/[V].** The measure-averaged curvature is `G″ = S − D`, `D = D_+ + D_{\rm Sav}`,
> `D_+ ≤ 0` (Brascamp–Lieb sign).  Measured: `D/S ≈ 0.97` (`SM2`, β=2.4, L=8), so the gap is the thin
> residue `G″ = O(3\%\,S)`.

**Consequence (the no-go).**  (i) A frozen-Hessian RMT computes `λ_min = −g|B| < 0` — it returns the
**negativity** and is *structurally blind* to `D_+` (the Brascamp–Lieb buy-in, a *measure* property,
not an operator-spectrum property).  (ii) On the measure-averaged object, an RMT/free-probability
estimate of `D` with **relative** error `δ` yields `G″` in `[G″ − δD, G″ + δD]`; this interval contains
`0` unless `δ < G″/D ≈ 3\%`.  Every available relative error is `≥10\%`:

| source of `δ` (from the *same* rank-one machinery, `SM2`) | `δ` | resolves `sign(G″)`? (need `δ<3%`) |
|---|---|---|
| finite-`c` vs `2/λ₄` (exact quadrature correction) | `10\%` | **no** |
| LLL overlap model (`SM2 §1.2`, `O(1)` lossy) | `O(1)` | **no** |
| free-prob leading cumulant (model) | `30\%` | **no** |

So **RMT/free-probability cannot resolve `sign(G″)` at the head even in principle at leading order** —
the gap is a `3\%` fine-tuning, below the method's intrinsic relative accuracy.

### 4.2 Witness (`w7/q7_rmt.py`, quoted)

```
(A) SAVVIDY OUTLIER IS EXACT (deterministic spectrum, not a Wigner fluctuation)
   gB= 1.0: lambda_min =  -1.00   positive edge =   1.00   outlier separation =   2.00 (=2 gB)
   gB= 4.0: lambda_min =  -4.00   positive edge =   4.00   outlier separation =   8.00 (=2 gB)
   => lambda_min<0 EXACTLY (outlier escapes; no BBP absorption). Frozen RMT -> G''<0.
(B) FINE-TUNING: gap G'' = S - D is a THIN residue of two O(S) quantities
   D/S=0.970:  S=1.400  D=1.358  G''= 0.0420   sign resolved only if (rel err on D) <  3.09%
   needed for sign resolution: delta < G''/D = 3.09%
      finite-c vs 2/lam4 (exact quadrature)        : delta ~  10.0%  -> DOES NOT RESOLVE
      LLL overlap model (SM2 Sec.1.2, O(1) lossy)  : delta ~ 100.0%  -> DOES NOT RESOLVE
   delta=10.0%:  G'' in [-0.0938,  0.1778]   sign fixed? NO (straddles 0)
(D) FREE-FIELD ACID TEST:  F=0 => no outlier; edge -> 0. RMT returns lambda_min = 0 (gapless). PASSES.
```
Readings: (A) the outlier is exact and deterministic (`−gB`, separation `2gB`) — frozen RMT returns the
*negativity*; (B) the `3.09%` fine-tuning is below every achievable relative error, so `sign(G″)` is
unresolved; (D) free-field acid test passes (returns `0`), confirming the tool is acid-clean and
therefore silent on the gap.

### 4.3 Exact missing step (Q2)

> **`[O]` (Q2).** *Relative*-error control `< G″/D ≈ 3\%` on the force variance
> `D = Var_μ(\partial_a V|a=0)` — i.e. the **variance of the variance** (a sub-leading free cumulant /
> the fourth cumulant), with a sign-definite contribution from `|D_+|` (the Brascamp–Lieb magnitude).
> No leading-order RMT / free-probability result provides this; it is the fine-tuning that *is* H1.

---

## 5 · Why these two are "the two most promising"

They are the only census entries that are (a) genuinely outside the five paradigms, (b) yield **sharp,
*computed*** obstructions (not hand-waves), and (c) each concentrate the entire failure into **one named
escape hatch**:

- **Q1 (direction / duality).** Closes the **entire analytic-information family** — Picard–Lefschetz's
  zero-set/strip shadow, Paley–Wiener, Fourier-uncertainty, and all of Fisher–Rao — with a *single
  sign*: every duality functional floors the variance.  This is a *unifying reason*, new to the corpus,
  and it *explains* the standing puzzle that the exact anchors (`0.2051`, `0.222`) are positive floors of
  tractable models but never transport: they are **ceiling-side** values.  Escape hatch = `Var ≤ 1/λ`.
- **Q2 (frozen outlier / fine-tuning).** Closes **RMT / free-probability**, and does so at the exact
  point the corpus's own `SM2` data points: `D/S ≈ 0.97`.  The `3\%` threshold is quantitative and
  *data-anchored* (`SM2`, `β=2.4`, `L=8`); the failure is a *precision* barrier, not a modeling gap.
  Escape hatch = `3\%`-accurate `D`, equivalently a lower bound on `|D_+|`.

Neither is a live wire to a proof; both are **the sharpest available statements of where any future
sixth-species attempt must break the symmetry of the problem**.  Together with `B6` (the moment/Laplace
family is one-sided by construction), they give a *three-way* classification of the entire non-five space:
**dualities floor the variance (Q1), moments floor the mass (B6), and precision cannot resolve the
fine-tuned sign (Q2).**

---

## 6 · Ledger

| item | statement | label |
|---|---|---|
| **Q1(a) Fisher identity** | `E_ρ[G″] = I(ρ) = E[(∂ log ρ)²]` | **[T]** (exact; witness `1e-14`) |
| **Q1(b) Cramér–Rao** | `I·Var ≥ 1` ⇒ `E[G″] ≥ 1/Var`; Gaussian equality | **[T]** (witness `I·Var=1.0000`) |
| **Q1(c) Bernstein** | band-limited `ρ` ⇒ `G″ = −ρ″/ρ ≤ τ²` (ceiling) | **[T]** (witness grid max ≤ τ²) |
| **Q1 directional no-go** | analytic/information family is a **variance floor**; wrong direction for H1 | **[NG]** (proof + witness) |
| **Q1 flat-limit** | fixed `τ`, `ρ_ε=1+εcos(τa)`: `G″=ετ²/(1+ε)→0` as `ε→0` (free field) | **[T]** (witness) |
| **Q1 exact-anchor reading** | `0.2051`/`0.222` are below the Bernstein ceiling ⇒ *ceiling-side* values | **[T]**/**[V]** |
| **Q2-a exact outlier** | Savvidy `λ_- = −g|B|`, deterministic, separation `2g|B|` (not BBP) | **[T]** (witness) |
| **Q2-b fine-tuning** | `G″=S−D`, `D/S≈0.97` ⇒ `G″/D≈3%` (measured `SM2`) | **[T]**/**[V]** |
| **Q2 precision barrier** | RMT rel-error `δ≥10% ≫ 3%` ⇒ `sign(G″)` unresolved | **[NG]** (witness table) |
| **(i) OT/Ricci** | 1-D `RCD ⇔ V''≥ρ` ⇒ = `Γ₂`, banked `−2.0` | **[NG]** |
| **(ii) PL thimbles** | real saddle is a Lefschetz saddle; sum-weights = measure control; analytic shadow = Q1 | **[O]**+**[NG]** |
| **(iii) information geom.** | Cramér–Rao variance floor; `W₂`-RG = `λ₁ = G″` | **[NG]** |
| **(v) MPS gap cert.** | subsumed by H6; no `d≥2` rigorous certificate | **[NG]** |
| **(vi) anomaly / θ** | constrains *phase* (existence), not magnitude | **[NG]** / **[T]** (phase) |
| **(vii) Kubo/modular/Lee–Yang/integrability/holo/Mori** | each = banked rename (spectral/RG/B3/Feshbach) or non-rigorous model | **[NG]** |
| **free-field acid tests (all)** | every family returns `0` on the massless field, non-vacuously | **[T]**+**[V]** |
| **H1** | `Var ≤ 1/λ` / `ess inf G″ ≥ λ` uniformly | **[O]** — unmoved |
| **Clay** | — | **0%** |

---

## 7 · Bottom line (no QED)

> **The sixth-species search is classified.**  Seven candidate families were proposed and first-tested;
> five are provably not new (they reduce to the banked `Γ₂`, `H1`-by-definition, the frozen Savvidy
> edge, a phase constraint, or the spectral no-go).  The two that survive are **two sharp, computed
> no-gos** with data-anchored obstructions:
> **(Q1)** analytic-duality functionals (Fourier/Paley–Wiener, Fisher–Rao, zero-set/strip — hence
> Picard–Lefschetz's analytic shadow) are **variance floors** — the wrong direction, saturated by the
> free field; **(Q2)** RMT/free-probability sees the **exact** deterministic Savvidy outlier `−g|B|`
> (returns the *negativity*, blind to the BL rescue) and cannot resolve the `3\%` fine-tuning
> `G″=S−D` with its `≥10\%` relative error.  Each has **one** escape hatch, and each hatch is, word for
> word, H1 (`Var≤1/λ`; `3\%`-accurate `D`).  **The gap lives in the shape region
> `{ρ bounded away from uniform}`; every sixth-species functional is one-sided by construction —
> dualities floor the variance, moments floor the mass, precision cannot resolve the fine-tuned sign.**
> **Clay = 0%. No QED.**

*Refs:* `docs/new-tools-registry.md` (A1–A5, B1–B8); `P-WAVE-SYNTHESIS.md`;
`STAR2-CLASSIFICATION.md`; `H6-HEAD-SPECTRAL.md`; `unconditional/PROOF-JACOBIAN-SAVES.md`;
`SM2-RANK-ONE.md`; `SM5-TCI-HEAD.md`; `H1-RG-DESCENT.md`; `docs/paper-v3/main.tex` §6.
Witnesses: `research/v3/w7/q7_direction.py`, `research/v3/w7/q7_rmt.py`.
Standard (textbook, not re-fetched): Bernstein/Paley–Wiener; Cramér–Rao; Prékopa (1971);
Brascamp–Lieb (1976); Bakry–Émery (1985); Lott–Villani–Sturm RCD; Picard–Lefschetz/Stokes;
Baik–Ben Arous–Péché (2005); Savvidy (1977). No third-party sources downloaded.

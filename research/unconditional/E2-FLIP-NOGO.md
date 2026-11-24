# E2-FLIP-NOGO — can *any* coupling-dependent dynamical "flip functional" break the E₂-regress, or is the regress irreducible for a precisely-stated natural class (with the surviving loophole named)?

> 2026-11-24 · **ADVERSARY pass** on the flip condition of `E2-REGRESS-CIRCULARITY §6`: *"the only thing that flips it = exhibit one non-Schur, non-regress, acid-test-clean trial functional `Φ`."* This note tries to **prove such `Φ` cannot exist in a named class** `C0`, by assembling a **four-horn dichotomy** from the banked `[T]`/`[NG]` lemmas (`LEVER-nonconvex-edge §0B`, `SCHUR-CONVEX`, `TOOL-COMPLETELY-MONOTONE`, `INVENT-NEW-INDEX` B3, `G2D-ISOLATION-I2` SSC) plus one fresh logical spine (Theorem 0 + the scale dichotomy). It is as much a test of the *adversary* as of the lever: a real loophole is reported as loudly as a no-go.
> **Labels (charter §2):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 (not theorem). External facts carried from the repo are **[R]**; the one fresh computation is **[V-here]** (`research/unconditional/w7/e2_flip_witness.py`). **No QED on anything unproven; this does NOT prove Clay and does NOT prove no proof can exist** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`, `CONVERGENCE-redteam §d`).

---

## Verdict up front

> ### **NO-GO PROVED [NG] for the class `C0`** = *functionals whose value is determined by the Källén–Lehmann spectral measure `ν` (equivalently the local Euclidean 2-pt data) at accessible scales*, accessed by any of four modes: **finite-moment (α), deformation-invariant (β), smooth/analytic/Tauberian (γ), quantitative-nonlinear-extremal (δ)**. For every `Φ ∈ C0` the dichotomy holds: **`Φ` either (a) cannot lower-bound `inf supp ν` at all, or (b) is circular (needs the gap / `E₂` / the mass scale / a certified eigenvector), or (c) reduces to `H1` (a uniform anti-concentration / coarse-mode variance bound), which is the open equivalent target — not an escape.**
>
> ### **LOOPHOLE THAT SURVIVES — named precisely: the complement of `C0`.** Two concrete, non-vacuous members:
> **(L-1) the vNW pointwise isolation** — same-channel von Neumann–Wigner *does* give a non-circular, acid-clean, coupling-dependent functional that isolates the lightest `0⁺⁺` below the excited cut **without assuming the gap**; it is **[T]-shaped and genuinely survives the three horns**, but it is **sub-uniform** (it gives a strict pointwise inequality, never a uniform `κ`; uniformizing it = SSC / the ratio `ρ` = `H*`).
> **(L-2) a field-space (non-spectral) sign-definite nonlinear inequality.** Any such `Φ` that is positive on YM and vanishes on free massless must forbid a sub-mass transverse mode, i.e. *is* an anti-concentration/propagation-floor statement — which is **`H1` by content**. So **(L-2) is not a break of the wall; it is the wall** (`H1`), and a genuine 6th species would have to be exactly this object (open `[O]`).

**So the honest score:** the E₂-regress is **irreducible for the UNIFORM statement within `C0`** (a real `[NG]`, assembled horn-by-horn), and the *only* known loopholes are (L-1) a real but **pointwise/sub-uniform** break, and (L-2) the `H1` reduction (the open target, not a shortcut). **No loophole survives that gives a uniform, non-circular, acid-clean lower bound and is *not* `H1`.**

---

## §0 The object, the target, and why "acid test" is *implied* (the logical spine)

**Setup (`BRIEF-FOR-MODEL §1`, `E2 §0.1`).** `μ` = reflection-positive Euclidean measure with transfer operator `T` and unique vacuum `Ω` (Perron–Frobenius, `PROOF-GAUSS-INVARIANT-PF`). For a gauge-invariant `O=¼tr F²`, the KL spectral measure is `ν(B)=⟨Ω|O†P_H(B)O|Ω⟩`, `supp ν⊂[m₁²,∞)`, `m₁=E₁−E₀`, atom `f₁²=ν({m₁²})`. Clay G2 = a uniform-in-`(L,a)` lower bound `m₁≥Δ>0` (a **support** statement; it does *not* need the atom, `SUPPORT-NOT-ATOM`).

**Definition (flip functional).** `Φ` is a *flip functional on an admissible family `F={μ}`* iff (LB) `Φ(μ)≤m₁(μ)` for all `μ∈F`, and `Φ` is **uniformly positive**: `inf_{μ∈F} Φ(μ) ≥ Δ>0`. The five task-constraints are: (i) `Φ` not Schur-concave; (ii) non-regress (no `E₂`/upper-bound input); (iii) coupling-dependent/dynamical (not deformation-invariant); (iv) acid-test-clean (`Φ(μ_free)=0`); (v) non-circular.

**Theorem 0 (the acid test is implied by (LB), not an extra filter). [T]**
Let `{μ_n}` be any admissible sequence with `m₁(μ_n)→0` (e.g. the free massive scalar as `m→0`; the `g→0`/`β→∞` limit of lattice YM is a gapless endpoint). Any `Φ` with `Φ(μ)≤m₁(μ)` satisfies `limsup_n Φ(μ_n) ≤ limsup_n m₁(μ_n) = 0`. In particular `Φ(μ_free)≤0` is **forced** by (LB); it is not an independent constraint.
*Proof.* Immediate from `Φ≤m₁` pointwise. ∎

> **Consequence (the correct framing of the whole problem) [T].** A flip functional on `F` is **exactly a uniformly-positive minorant of the gap function `β↦m₁(β)` on `F`, computable from accessible data.** "Construct a flip functional" ≡ "prove `m₁≥Δ` on `F`." The only non-trivial content is therefore **computability from accessible data** — i.e. *which data*, and whether a minorant of `m₁` can be squeezed out of it *without already containing `m₁`*. This is why the regress is not a curiosity: it is the statement that a lower bound below the true value is, by definition, a proof of the gap. The dichotomy below is a classification of the *data* by which such a minorant could be certified.

---

## §1 The class `C0`, stated precisely

**`C0` := the functionals `Φ` whose value is a function of the KL spectral measure `ν` of some gauge-invariant interpolator (equivalently of the local Euclidean 2-/multi-point data at accessible scales), and whose evaluation uses `ν` through one of four access modes:**

| mode | `Φ` is a function of … | paradigm |
|---|---|---|
| **(α) finite data** | finitely many `∫φ_i dν` (moments, generalised moments, f-sum, susceptibilities) or any finite list of smooth spectral statistics | moment problem, Stieltjes/Christoffel |
| **(β) invariant data** | a datum constant on continuous deformations of `μ` (homotopy type, index, signed count, `E_d`-type, derived/TQFT invariant) | Morse/index/spectral flow |
| **(γ) smooth/analytic data** | the full `ν` through a smooth functional (`C(t)`, its log-derivatives `r(t)`, analytic continuation, singularity/growth data) | completely-monotone / Tauberian |
| **(δ) quantitative-extremal data** | `ν` through an extremal/quantitative functional (capacity, Cheeger, min-support/quantile, large-deviation rate, transport/Wasserstein) | potential theory / optimal transport |

`C0` is the *union of the four modes*. It is not a contrivance: it is "every way to read a scalar off a positive spectral measure by finite data, by an invariant, by analytic structure, or by a variational/quantitative extremum." Every tool in the `~110`-route death map is an element of `C0`.

**Remark (what `C0` excludes — the loophole boundary, stated up front).** `C0` excludes functionals that use **structure beyond the spectral measure** — a *sign-definite nonlinear inequality on the field-space measure itself* (a correlation/anti-concentration inequality not reducible to spectral data), a *certified near-eigenvector*, or `E₂`/the isolation margin. §5 shows these are exactly `H1` (open target) or `H*`/`I2` (the wall) or the vNW pointwise break.

---

## §2 The four-horn dichotomy — each horn killed with an exact statement

### Horn (α) — finite data ⇒ cannot lower-bound the edge  **[T]**
**Theorem α (finite-data blindness).** For any finite set `φ_1,…,φ_k` and any moment target `C=(C_i)`, the solution set `{ν≥0 : ∫φ_i dν=C_i}` — when non-singleton — contains measures whose support floor is *arbitrarily far below* that of any particular member. Hence **no `Φ=F(∫φ_i dν)` lower-bounds `inf supp ν`** (uniformly over the fiber).
*Proof.* The classical `δ₀+δ₂` vs `2δ₁` witness (`E2 §0.5`, `LEVER §0B`): identical `(M₀,M₁)`, floors `0` vs `1`. The `[V-here]` **W1** extends this to **order-3 moments exactly**: a measure with `inf supp=0` (mass `w₀=0.0275` at `s=0`) shares `M₀,M₁,M₂,M₃=(1, 1.995, 4.4105, 10.544)` with a gapped reference of floor `1.0` (LP residual `2.7×10⁻¹⁴`). The general perturbative form is the **linear-data lemma**: for any finite exact data and any `ε` there is `ν'` with an atom in any neighbourhood of `0` and the same data up to `O(ε)`. ∎ (The one protected negative-power total that *would* help, `χ_t=M_{−1}`, is in the `0⁻⁺` channel and **unprotected/upper-type**, `SCHUR-CONVEX §1`.)

### Horn (β) — deformation-invariant data ⇒ frozen at the free value  **[NG]**
**Theorem β (rigidity blindness).** Let `μ ↦ Φ(μ)` depend only on a datum constant on continuous deformations. Then `Φ` is locally constant on the connected `β`-family `(0,∞)`, hence `Φ(μ_β)=Φ(μ_free)` (the free `β→∞` / `g=0` endpoint). By (LB) and Theorem 0, `Φ(μ_free)≤0`; therefore `Φ≤0` on the whole family ⇒ **not a uniformly positive minorant.**
*Proof.* `INVENT-NEW-INDEX §3` rigidity (B3): a discrete/homotopical invariant on a connected parameter family is constant; the free theory lies in the closure; `m₁(free)=0` forces `Φ(free)≤0`. ∎ This is the exact content of `TOOL-FACTORIZATION-ALGEBRA` (two-tier blindness), `TOOL-MORSE-ATIYAH-BOTT` (index strata), `FORGE-SPECTRAL-FLOW`, `PROOF-MORSE-FLOER`. **Verdict: `[NG]` — a deformation-invariant functional cannot be a positive lower bound.**

### Horn (γ) — smooth/analytic spectral functionals ⇒ upper bounds or the answer  **[T]/[NG]**
**Theorem γ.** For `C(t)=∫e^{−st}dν ≥0` (completely monotone ⟺ `ν≥0`, Bernstein–Widder):
1. *(finite-t rate is an upper bound)* `r(t):=−(1/t)\log(C(t)/C(0))` satisfies `m₁² ≤ r(t) ≤ M₁/M₀` and `r(t)↓m₁²` as `t→∞`. So every finite-`t` smooth scalar yields `m₁² ≤ r(t)` — an **upper** bound, never a lower bound.
2. *(the limit is the answer)* `m₁² = lim_{t→∞} r(t)` is an **exact Tauberian identity**; using it requires `C(t)` controlled uniformly on `t≳1/m₁` = exponential clustering = the gap (**circular**).
3. *(analytic continuation is floor-blind)* the edge atom contributes the **entire** function `f₁²e^{−m₁²t}` (no singularity); the half-plane of convergence is `{Re t>0}` for **every** floor (fixed by the UV extent). Radius/singularity/growth data is **independent of `m₁²`**.
4. *(no extra rigidity)* complete monotonicity `⟺ ν≥0`; there is no side condition beyond positivity.
*Proof.* `TOOL-COMPLETELY-MONOTONE` §1–§3 (direction theorem, atom-entire, CM≡ν≥0); `[V-here]` **W2** confirms (1)+(acid): on the massless-type `C(t)=(1+t)^{−2}` (floor `0`), `r(t)=1.62,1.39,0.92,0.48,0.23,0.09` at `t=0.5,1,3,10,30,100` — **positive at every finite `t`**, so the finite-`t` rate is positive exactly where it must be `0` (vacuous), and only `t→∞` returns `0` (=the answer). ∎ **Verdict: `[NG]` Family-I.**

### Horn (δ) — quantitative nonlinear extremal ⇒ reduces to `H1` or to the effective action  **[NG]/[O]**
**Theorem δ (each quantitative extremal lands on `H1` or on the answer).** The nonlinear extremal functionals below each **either** reduce to the coarse-mode anti-concentration identity `H1`, **or** pin the support via a variational potential that *is* the effective action:
- **Capacity / Wiener / Cheeger:** the `L`-uniform conductance floor `Φ_eq≥c` *is* `H1` (`POTENTIAL-THEORY-CAPACITY`, `HISTORICAL #4`; `∫Γ₂ dμ=‖Lf‖²=λ₁` identity).
- **Large-deviation tilt:** the conditioning that would fix an edge sign is an LD tilt whose **rate function is the effective action** = `H1` (`DYSON-GAS §4`, `E4-RMT-EDGE §3`).
- **Transport / Wasserstein Γ₂ (Otto–Bakry–Émery):** `∫Γ₂dμ=‖Lf‖²=λ₁` — *exactly* the gap by an identity = `H1` (`BRIEF §3` row "最优传输").
- **Min-support / quantile of `ν`:** computing it needs the edge itself (the answer); it is `LB` only because it *is* `m₁` renamed.
*Proof / attribution.* Each is a banked `[NG]`/`[O]` (`HISTORICAL-TOOL-SWEEP` §2 items 4,9,11,13; `BRIEF §3`). It is one object: `H1` (residual-dictionary `D.1`) or `H*`. ∎ **Verdict: `[NG]` → `H1` [O].**

---

## §3 The scale dichotomy — why a *dimensionful* flip functional is necessarily circular  **[NG-for-class | T-conditional on SSC]**

This is the sharpest *new* adverserial spine (the one the task flags: *"any `Φ` vanishing on free massless must scale with the same coupling-dependent quantity that the gap does"*). It closes the dimensionful horn that §2's modes leave slightly implicit.

**Theorem S (dimension dichotomy).** Let `Φ` be a flip functional (§0) uniform in `(L,a)`. Then:
- **(D) dimensionful (`Φ` is an energy):** pure 4D YM is classically scale-free; the only intrinsic scales on an `(L,a)`-uniform gauge-invariant local functional are `{Λ_QCD, √σ, g|B|, 1/a, 1/L}` (cutoff/volume/Savvidy). Uniformity kills `1/a,1/L`; `g|B|` is UV/static-curvature (`TEST-C-GLOBAL-SCHUR`: the `~0.6β` mass is *not* the transfer gap); `√σ` is **1-form/cross-sector** (`BS-THRESHOLD §4`). Leaves **`Λ_QCD`** (dimensional transmutation, the *only* scale). So `Φ=Ψ(Λ_QCD)`, and (LB) requires `Ψ(Λ_QCD) ≤ m₁ = c·Λ_QCD` with `c=m₁/Λ_QCD`. **Bounding `m₁` therefore requires knowing the dimensionless ratio `c=m₁/Λ_QCD` — which is the gap** (`circular`).
- **(d) dimensionless (`Φ` bounds a shape datum):** then `Φ` bounds the atom `f₁²`, the isolation margin `Λ−m₁²`, or the ratio `ρ=m_{0⁺⁺*}/m₁`. By **L2** (`GAP-G2-edge`) the atom ⟺ isolation `E₂≥E₁+δ`; by **SSC** (`G2D-ISOLATION-I2 §2`, `E₂=ρE₁`) isolation ⟺ the gap up to the universal `ρ`; and bounding `ρ` is `H*` (`G2D §5`). **`[O]`.**

*Proof.* The scale inventory (D) is `W3-adjacent-4D` (marginal/scale-free: must generate its own `Λ`) + `BS-THRESHOLD` (the only other scale `√σ` is 1-form) + `TEST-C-GLOBAL-SCHUR` (`g|B|` is the wrong object); the ratio step is definitional. The `(d)` half is `L2`+`SSC`. The single-scale hypothesis (`c` β-independent) is `[J]`-empirical, so (D)'s *circularity* conclusion is `[T]`-conditional on SSC; the *dimensional analysis* `[T]`. ∎

> **Sharpest line (the adversary's core).** *A flip functional that vanishes on the free massless field and is uniform in `(L,a)` is a **scale-detector**; the acid test says it must vanish exactly where the theory's scale vanishes (`Λ=0`); but the theory's only scale `Λ` is the gap in disguise (`m₁=c·Λ`), and detecting "a mass is present" from accessible local data *non-circularly* means showing the local measure has no weight below that mass — which is `H1`.* The acid test does not just *permit* the same coupling-dependent quantity as the gap — it **forces it**.

---

## §4 The dichotomy theorem for `C0`, assembled  **[NG]**

**Theorem (four-horn dichotomy).** For every `Φ ∈ C0` satisfying (LB) and (iv):
> `Φ` is **(a)** not a lower bound (§2-α finite data, §2-γ finite-`t` smooth data; both are *upper*-bound/inconsistent), **or (b) circular** (§2-β rigid ⇒ `≤0`; §3-D dimensionful ⇒ needs `m₁/Λ`), **or (c) reduces to `H1`** (§2-δ capacity/Cheeger/LD/transport; §3-d atom/isolation ⇒ `L2`/`SSC`, or `H1` via the coarse-mode identity).
> There is **no fourth, non-circular, uniform, lower-bound output.**

*Assembly.* `C0 = α∪β∪γ∪δ` by definition. α, γ: `[T]` (Theorems α, γ). β: `[NG]` (Theorem β, rigidity B3). δ: `[NG]/[O]` (Theorem δ, banked per-species). D/d: §3. ∎

> **Label: `[NG]` — no-go for the class `C0`, uniform statement.** This is assembled from `[T]` lemmas on the mathematical horns (α, γ are theorems about operators/measures, tool-independent) and `[NG]` on the rigidity/quantitative horns. It is **not** a global barrier: `C0 ≠ all possible` (§5, §6).

**Why this is not merely a re-list of the death map.** Two things are new here relative to `E2-REGRESS-CIRCULARITY`:
1. **Theorem 0** recasts "acid-test-clean + lower-bounding" as *a uniformly-positive minorant of `m₁` computable from data* — so the regress is the tautological-but-decisive statement that a minorant below the true gap is a proof of the gap, and the only question is the **data mode**.
2. **The scale dichotomy (§3)** closes the one horn (`dimensionful`) that the four "positivity-family" horns leave ambiguous: a uniform dimensionful `Φ` is forced onto `Λ_QCD`, and `Λ_QCD` is the gap. This is `HISTORICAL §3`'s "sector/threshold mismatch" question answered on the scale axis: there *is no second scale* to mismatch onto, so `D` inevitably lands on `c=m₁/Λ`.

---

## §5 The loopholes that survive — named precisely (the real path)

A genuine loophole is more valuable than a fake no-go. Two survive, and one class-boundary must be stated honestly.

### (L-1) The vNW pointwise isolation — **a real, non-circular, sub-uniform break**  **[T-shape]**
Same-channel von Neumann–Wigner forbids crossing of two `0⁺⁺` levels; strong-coupling cost-grading puts `m₁<m_{0⁺⁺*}` and the ordering is sign-constant on the connected `β`-line (confinement is a crossover). Hence `m₁<m_{0⁺⁺*}` **pointwise wherever both are bound**, *without presupposing `m₁>0`* — a genuine, acid-clean (no free-field tower to order), coupling-dependent functional that the horns α–δ **do not kill**, because it is not a functional of `ν` alone but of the *ordering of two levels* (`G2D-ISOLATION-I2 §3`, `PROOF-NO-DISSOLUTION-KINEMATIC`). **It survives all five constraints at the POINTWISE level.** Why it is not Clay: (i) the excited `0⁺⁺*` **dissolves into the continuum at large `β`** (Kato isolation fails there); (ii) it gives strict inequality, never a **uniform `κ`**; (iii) uniformizing — `ρ−1` bounded below — is **verbatim `H*`** (`G2D §5`). **Verdict: `[T]` pointwise survivor; not uniform.** This is the *one* place the E₂-regress is genuinely broken, and the honest adversary must yield it.

### (L-2) Beyond the spectral measure: the field-space nonlinear inequality  **[O] ≡ `H1`**
`C0`'s complement = functionals using the **field-space measure**, not just `ν`. Any such `Φ` that is (iv) acid-clean and (LB) positive must forbid a transverse mode propagating below the mass — i.e. must assert the local measure has **no spectral weight at scales `<m₁`**, which is **uniform anti-concentration of the coarse modes = `H1`** (residual-dictionary `D.1`; `PROOF-COARSE-POINCARÉ`: `Var_ν(φ)≤(1/λ)∫|∇φ|²`, `λ>0` uniformly ⟺ `m_eff²(k→0)≥λ`). **So (L-2) is not a shortcut — it *is* the open target.** A "6th species" (charter `§1`) is precisely a proven-sign, YM-specific, free-field-failing field-space inequality of this kind; none is known, and constructing one is the missing machine. **[O]/[J], prior negative.**

### (L-3) The "continuous except at the crossover" class — a β-critical non-analytic functional  **[NG]/[O → transport≡confinement]**
The task names "any `Φ` continuous in the coupling except at the crossover." Such a `Φ` evades horn **β** (it is *not* deformation-invariant — it jumps at `β_c`). Is it a loophole? **No, in the examined forms:** a functional that jumps only at the confinement crossover must *detect* the onset of confinement (an area-law coefficient, percolation-type order parameter, IR mass generation). Detecting that nonperturbatively across the crossover **is** `transport ≡ confinement` (`BRIEF §4`, `SMALLCIRCLE-ADIA`, `CHATGPT-METHOD-RGENTRY`): either it is circular (uses the path gap to move the gap) or it *silently proves confinement*. Moreover a sharp `β_c` exists only in the thermodynamic limit, so an `(L,a)`-uniform `Φ` cannot lean on it. **Verdict: this class is real as a *shape* but reduces to `transport≡confinement` `[NG]` — no escape.** (It is, however, the *only* non-invariant way to evade horn β, so it is the precise object to watch: a β-critical functional that is neither rigid nor transport-circular would be a genuine break; none is known.)

---

## §6 Honest boundary — what this does NOT prove

- It does **not** prove Clay and does **not** prove no proof exists (`∀class ≠ ∀possible`; `KEY-EXISTS §5`, `CONVERGENCE-redteam §d`). It proves a **`[NG]` for the class `C0`** (four horns + scale dichotomy) and certifies that `C0`'s complement contains **no uniform, non-circular, acid-clean escape that is not `H1`** — but `C0`'s complement is **not proven empty**; a not-yet-invented 6th species could live there, and by (L-2) it would *be* `H1`.
- The **vNW pointwise (L-1)** is **not** under-sold: it is a real break of the regress at the pointwise level (`[T]`-shaped). It simply cannot furnish the **uniform** `(L,a)` margin Clay needs; uniformizing it is `SSC`/`ρ` = `H*`.
- **Theorem S (D)**'s *circularity* conclusion is `[T]`-conditional on the Single-Scale hypothesis (`ρ` β-independent, `[J]`, empirical `ρ≈1.5–1.6`); the *dimensional analysis* (few scales ⇒ `Λ_QCD`) is `[T]`.
- The fresh computation (`W1–W3`, `[V-here]`) is a two-minute check of already-banked lemmas (linear-data lemma, direction theorem, moment-ratio direction); the substantive numerics/partials (`Temple 0/3474`, vNW, `ρ`) are inherited `[R-internal]`.
- **Flip condition (unchanged, now sharpened).** This verdict flips to **BROKEN** iff someone exhibits a `Φ` that is **(i)** a uniform `(L,a)` lower bound, **(ii)** acid-clean, **(iii)** *not* in `C0` (uses structure beyond the spectral measure), **(iv)** *not* `H1`, and **(v)** *not* the vNW pointwise statement nor transport-circular. By §3 such a `Φ` must be a **dimensionless, non-`H1`, non-`L2`/`SSC` field-space invariant of pure YM** — the name of the 6th species. None is known; finding it is the task.

---

## §7 References

**Internal (absolute paths; the banked lemmas this note assembles):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall/machine); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five constraints; §3 death map; §4 transport≡confinement); `D:\Dev\math\docs\new-tools-registry.md` (§D residual dictionary D.1–D.7; B3 rigidity).
- `D:\Dev\math\research\unconditional\E2-REGRESS-CIRCULARITY.md` (§0.2 L1; §0.4 SSC; §0.6 two horns; §6 flip condition — the target of this note); `HISTORICAL-TOOL-SWEEP.md` (#12 Tauberian; §2 rows 4,9,11,13; §3 sector/threshold); `KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 missed lever; §7 calibration); `CONVERGENCE-redteam.md` (§b3 non-convex edge; §d `∀class≠∀possible`).
- `LEVER-nonconvex-edge.md` (§0A free-field adversary; §0B linear-data lemma — Theorem α); `SCHUR-CONVEX.md` (master Schur no-go; `χ_t` unprotected — Theorem α caveat); `TOOL-COMPLETELY-MONOTONE.md` (§1–3 direction theorem, atom-entire, CM≡ν≥0 — Theorem γ); `G2D-ISOLATION-I2.md` (§2 SSC; §3 vNW pointwise; §5 ratio=`H*` — §3, L-1).
- `INVENT-NEW-INDEX.md` (B3 rigidity — Theorem β); `TOOL-FACTORIZATION-ALGEBRA.md` (two-tier blindness); `TOOL-MORSE-ATIYAH-BOTT.md` (index/gap-above strata); `FORGE-SPECTRAL-FLOW.md`, `PROOF-MORSE-FLOER.md` (rigidity instances).
- `GAP-G2-edge.md` (residue–level-gap lemma L2; moments upper); `SUPPORT-NOT-ATOM.md` (support vs atom); `PROOF-COARSE-POINCARÉ.md` (H1 = coarse Poincaré = `m_eff²`); `POTENTIAL-THEORY-CAPACITY.md`, `B3-ISOPERIMETRIC.md` (capacity ⇒ H1); `DYSON-GAS.md` (§3 edge=potential=gap; §4 LD tilt=H1); `E4-RMT-EDGE.md` (LD conditioning = H1 — Theorem δ); `BS-THRESHOLD.md` (`√σ` 1-form — Theorem S); `TEST-C-GLOBAL-SCHUR.md` (`g|B|` ≠ transfer gap); `W3-adjacent-4D.md` (scale-free ⇒ generates `Λ`); `PROOF-NO-DISSOLUTION-KINEMATIC.md` (vNW ordering — L-1); `SMALLCIRCLE-ADIA.md`, `CHATGPT-METHOD-RGENTRY.md` (transport≡confinement — L-3); `PROOF-GAUSS-INVARIANT-PF.md` (vacuum pinned at 0).
- `D:\Dev\math\research\unconditional\w7\e2_flip_witness.py` (`[V-here]` W1–W3; run `py -3.12 e2_flip_witness.py`).

**External — [standard]/[R] (textbook; not re-verified this session):** Bernstein (1928/29), Widder *The Laplace Transform* (1941) — completely monotone ⟺ `ν≥0`; Temple (1928)/Weinstein (1932)/Kato (1949); von Neumann–Wigner (1929); Markov–Krein–Nudelman (moment problem); Bakry–Émery / Otto (transport `Γ₂`); Wiener (1924)/Kakutani (1944) capacity; Ben-Arous–Guionnet (large deviations). No fabricated bound; no QED.

---

*Status: **[NG] for the class `C0` (four-horn dichotomy: finite-moment / deformation-invariant / smooth-analytic / quantitative-extremal, plus the §3 scale dichotomy), uniform statement — assembled from `[T]` lemmas (α, γ, Theorem 0, scale inventory) and `[NG]` (β rigidity, δ capacity/LD/transport). LOOPHOLE THAT SURVIVES: (L-1) vNW pointwise isolation (`[T]`-shaped, genuine but sub-uniform; uniformizing = SSC/`ρ`=`H*`); (L-2) the field-space nonlinear inequality `≡ H1` (the open target, not a shortcut); (L-3) β-critical functionals `→ transport≡confinement` `[NG]`. The E₂-regress is irreducible for the UNIFORM statement within `C0`; the only non-`H1` loophole is the pointwise vNW break, which does not close Clay. NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay stays open.***

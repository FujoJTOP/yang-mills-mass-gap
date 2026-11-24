# E2-REGRESS-CIRCULARITY — can the "reverse tool's input = the answer" regress be broken by a genuinely NON-CIRCULAR anchor, or is it the irreducible wall?

> 2026-11-23 · Attack on the **input-circularity** obstruction of Family II (Temple–Kato–Weinstein / Birman–Schwinger / spectral flow / capacity), the deep wall named in `HISTORICAL-TOOL-SWEEP §2–3` and `KEY-EXISTS-OR-IMPOSSIBLE §5`. This is the structural question *"the reverse tool's input is the very answer."* Scope is narrow and honest: does any anchor give a **uniform lower bound on `E₁` (the lightest above-vacuum energy) without presupposing `E₂`/the gap**?
> **Labels (charter §2):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 (judgment, not theorem). External facts carried from the repo are **[R]** (as-reported; I re-ran no external source this session). The one fresh computation here is the two-line witness in §0.5, labelled **[V-here]**. **No QED on anything unproven; this document does NOT prove Clay, and does NOT prove that no tool can exist.**

---

## Verdict up front

> ### **IRREDUCIBLE for the uniform lower bound (Clay's G2) — with ONE genuine non-collapsing corner (same-channel von Neumann–Wigner *pointwise* isolation) that breaks the circularity pointwise but cannot supply the uniform margin.**

**BROKEN?** No — not for what Clay needs. **The precise obstruction, in one paragraph:**

> For the *uniform* statement Clay needs, Temple/Weinstein/Kato reduce a lower bound on a bottom eigenvalue to a lower bound on the **next** level (`ρ ≤ E₂`), i.e. to *isolation* away from the lightest state. In pure Yang–Mills, a **one-scale** theory, the second `0⁺⁺` level is not an independent object: it is the first radial excitation `0⁺⁺*`, and `m₂ = ρ·m₁` with a universal (β-,`L`-independent) ratio `ρ` (`Single-Scale Covariance`, §0.4). Hence *"uniform isolation `m₂²−m₁² ≥ κ`"* and *"uniform gap `m₁ ≥ μ`"* are **one inequality up to the fixed number `ρ`** — supplying `E₂` independently is impossible because `E₂` is not independent. The only non-circular edge datum that non-perturbatively forces a bottom eigenvalue is **the atom itself** (positivity forces only `m=0`; Goldstone uniqueness — `GOLDS-POS-MASS`). Therefore every Family-II input is one of exactly three things: **(i) the answer one rung up** (`E₂ = ρE₁`), **(ii) the answer as a certified near-eigenvector** (Temple's sharpness needs `Var_⊥` small = the glueball), or **(iii) an upper bound** (Schur-concave positivity, `SCHUR-CONVEX`). The regress is not an infinite ladder; in a one-scale theory it **terminates at a one-step identity `E₂ = ρE₁`**, and that identity is the circularity.

**What is *not* claimed:** this is a **[NG] for the identified classes** (Family-I Schur-concave; Family-II regress/rigidity) plus the **exact missing step [O]** for the uniform statement, with a **[T]-conditional structural reason** (Single-Scale Covariance). It is **not** a global impossibility theorem: `∀class ≠ ∀possible` (`CONVERGENCE-redteam §d`, `KEY-EXISTS §5`). The verdict is a *precise negative*, not a barrier theorem.

**The one genuine partial** (do not under-sell it): same-channel **von Neumann–Wigner** level repulsion (§2.4, `G2D-ISOLATION-I2 §3`) gives *pointwise* isolation of the lightest `0⁺⁺` below the excited cut **without assuming the gap** — a real break of the circularity at the pointwise level. It is a [T]-shape partial; it does not reach a uniform `κ`, so it does not close Clay.

---

## §0 The exact logical schema of the obstruction

### 0.1 The target object
PF closes the vacuum *exactly* (`PROOF-GAUSS-INVARIANT-PF`: Perron eigenvector = vacuum, `E₀` = free energy), so `m₁ = E₁ − E₀` reduces **exactly** to a lower bound on `E₁` on `(vac)^⊥` (`TEMPLE-WEINSTEIN §2`). The KL measure `ν ≥ 0`, `supp ν ⊂ [m₁², ∞)`, atom `f₁² = ν({m₁²})`. Clay's G2 needs `spec(H) ∩ (0, Δ) = ∅`, i.e. a uniform `m₁ ≥ Δ > 0` — a **support** statement (it does *not* need the atom `f₁²>0`; `SUPPORT-NOT-ATOM §1`, the massive free `:φ²:` is gapped with `f₁²=0`).

### 0.2 Lemma L1 — mean + variance cannot pin inf-support  **[T]**
*(`TEMPLE-WEINSTEIN §3`, re-stated abstractly.)* For any trial, the only data a pure variational calculation exposes about the spectrum below are `R = ⟨H⟩_ψ` and `Var_ψ(H)`. **No function `Φ(R, Var)` lower-bounds `inf supp ν`.** Fix `R`, `Var`; put mass `ε = 1/N²` at `−N` and `1−ε` at `s = −εt/(1−ε)`. Then `R` and `Var` are unchanged but `inf supp ν → −∞` (and can be made `0` in the physical non-negative frame). **Consequence:** *every* valid trial lower bound must smuggle a **hard spectral edge** `E₂` (Temple/Weinstein) or a **certified near-eigenvector** (small `Var_⊥`); neither is a function of `(R, Var)` alone. This is the energy-space face of `GAP-G2-edge`'s "moments cannot pin the edge."

### 0.3 Lemma L2 — residue ⟺ level gap  **[T]**
*(`GAP-G2-edge`, residue–level–gap lemma.)* If `ν = f₁²δ_{m₁²} + ν_rest`, `supp ν_rest ⊂ [Λ, ∞)`, `Λ > m₁²`, then
`f₁² ≥ (ΛM₀ − M₁)/(Λ − m₁²)`, equality iff the rest is one atom at `Λ`.
**Consequence:** a lower bound on the atom (`I1`: `f₁² ≥ r₀`) is, given the total moments, **equivalent** to a lower bound on the second level (`I2`: `Λ ≥ m₁² + κ`). Edge control at level `n` needs level `n+1`. The atom, the isolation, and the gap are one wall seen three ways.

### 0.4 Lemma L3 — Single-Scale Covariance (SSC): `E₂ = ρ·E₁`  **[T]-conditional + [J]**
*(`G2D-ISOLATION-I2 §2`.)* Pure YM is marginal and scale-free (`W3-adjacent-4D`: "must generate its own `Λ`"). Its spectrum has a **single dynamical scale**, so every glueball mass is `c·Λ_QCD` and the ratio `ρ = m_{0⁺⁺*}/m₁` is a fixed number, β-,`L`-,`a`-independent. Then, with `m₂² = ρ²m₁²`,
> `m₂² − m₁² = (ρ² − 1)·m₁²`.
Hence for fixed `ρ > 1`:
> **uniform isolation  `m₂² − m₁² ≥ κ`  ⟺  uniform gap  `m₁ ≥ √(κ/(ρ²−1))`.**  *(both directions)*
The **`ρ>1`** half is the von Neumann–Wigner ordering §2.4 (**[T]** given both states bound + strong-coupling ordering); the **`ρ` bounded above / β-independent** half is the one-scale hypothesis (**[J]**, empirically `ρ≈1.5–1.6`, `PROOF-ISOLATION §4`). *The structural content:* in a one-scale theory, "isolated by `κ`" and "gapped by `μ`" are **the same sentence up to a universal constant**. Supplying `E₂` from a "different, provable" source is impossible because `E₂` is not a different object — it is `ρ` times the answer.

### 0.5 Fresh witness of L1 + the kinematic-`ρ` overclaim  **[V-here]**
Pure-Python, this session (no repo state touched; banked numerics such as Temple `0/3474`, commutator `5000/5000`, naive overclaim are **inherited [R-internal]** from `TEMPLE-WEINSTEIN`/`TEMPLE-RIGOROUS-H`):
```
L1: nu1=delta0+delta2 and nu2=2*delta1  ->  M0=2.000, M1=2.000 (identical), inf_supp 0 vs 1.  ✔
Temple (E0=0,E1=1,E2=1.5; trial R=1.0100, Var=0.0049, true E1=1.0):
   rho = E2 = 1.5 (VALID)        -> RHS=1.0000 <= E1  (sound)
   rho = 2m1 = 2.0 (>E2, INVALID)-> RHS=1.0051 >  E1  (OVERCLAIMS)
   rho = 3.0      (INVALID)      -> RHS=1.0075 >  E1  (OVERCLAIMS)
```
The second/third rows are the concrete failure of the "kinematic `E₂`" escape: a threshold set **above** the true dynamical second level is an **unsound** `ρ` and returns a false too-high "lower" bound. (§2 makes this the closure of the `(d)` anchor.)

### 0.6 The Temple schema — two horns, no third
Temple on `H|_{(vac)^⊥}` (bottom `E₁=m₁`, second `E₂=m₂`): if `ρ ≤ m₂` and `m₁ ≤ R < ρ`, then `m₁ ≥ R − Var/(ρ−R)`.
- **Horn (T-i):** valid only with a **lower** bound `ρ ≤ m₂` = *isolation* `I2`. In one-scale YM, `m₂ = ρ·m₁`, so by L3 this **is** the gap.
- **Horn (T-ii):** sharp only if `Var_⊥(ψ₁) ≤ (R_⊥ − E₀ − λ)²` — i.e. `ψ₁` is a **certified near-eigenvector**, a certified glueball trial (`f₁² ≥ z₀²`, the `(M1)` overlap). That is the object being bounded (`PROOF-F1SQ`-class circularity).
- **No third horn.** L1 is a theorem about operators, not tools: `(R, Var)` alone pins no edge. Weinstein/Kato are strictly *more* demanding and land on the same `I2` (`TEMPLE-RIGOROUS-H §4`).

---

## §1 Scorecard of candidate non-circular anchors

| anchor | (i) gives **LOWER** bound on `E₁`? | (ii) free-field acid test | (iii) relocates circularity? | verdict |
|---|---|---|---|---|
| **(a)** variational upper + "no continuum below" | **No** — min–max is a *min* = upper (`Courant–Fischer`) | n/a (upper) | "no continuum below" **is** the gap | **[NG] / tautology** |
| **(b)** reflection positivity / OS | **No** — PSD ⇒ Cauchy–Schwarz ⇒ upper; RP is scale-covariant | free massless is RP with `E₁=0` ⇒ no floor (correct) **but no floor on YM either** (vacuous) | best lower-shaped output = Bushell cone rate `=G″=H1`, *insufficient* (`SUPPORT-NOT-ATOM §3`) | **[NG] (vacuous) / relocates to H1** |
| **(c)** trace anomaly / OPE / sum rules | **No** — anomaly is a *moment* `∫ρ/s²∝Λ⁴`, upper-type | free Maxwell: `a_free≠0`, `f₁²=0` ⇒ 4-pt vertex `V=∞`, bound vacuous | `V` diverges exactly where `f₁²` vanishes ⇒ no `f₁²` isolation (`INEQ-R6`) | **[NG] / [O]** |
| **(d)** kinematic `E₂` / multi-particle cut | **No** — collapses to `I2` (§2) | free: no isolated constituent, threshold `0`, returns nothing (correct, vacuous) | two independent closures: next level *dynamical* (§2.2) + SSC `E₂=ρE₁` (§2.3); **but** a genuine *pointwise* corner survives (§2.4) | **[NG uniform] + real pointwise partial** |
| **(e)** symmetry / 't Hooft anomaly / a-theorem | **No** — anomaly matching is a *consistency* condition (forbids *symmetric trivially-gapped*, i.e. forces gapless/topo order), **opposite direction**; a-theorem is upper on a dimensionless number | 1-form anomaly blind to local `tr F²` (a point carries no 1-form charge); θ=0 anomaly trivial ⇒ no force | θ=π forces at most a TQFT (still gapped), never a local mass | **[NG]** |
| **(f)** two-sided / bootstrap / exclusion | **No** — conformal bootstrap: YM not a CFT, no crossing; S-matrix bootstrap *needs* asymptotic atoms = `I1∧I2` = `H*` (**stronger** than gap) | free Maxwell is crossing-exact ⇒ any lower-bound exclusion would wrongly kill the free field (fails acid test) | "exclude the gapless point" = exclude everything or nothing; feasible ⇒ not excluded; forcing = add `Λ` = `H1` renamed | **[NG] / relocates** |
| **(g)** any unconsidered anchor | **No** (§4) | each either Schur-concave (upper) or needs a second scale pure YM lacks | strongest untried = ratio `ρ` = `O(u)` = `H*` | **[O] / [J]** |

**Read the table as the precise negative:** every row is **direction-correct at most in form**, **acid-test-vacuous or acid-test-incompatible**, and **relocates the circularity onto `H1` (edge/curvature), `I2` (isolation), or `f₁²>0` (atom)** — the same wall in different coordinates.

---

## §2 Deep dive on (d): the kinematic-`E₂` escape — closed, with one real corner

This is the escape §3 of `HISTORICAL-TOOL-SWEEP` flagged as the *only* conceivable non-collapsing shape ("a sector/threshold mismatch giving a genuinely kinematic `E₂`"). We re-examine it in full.

### 2.1 The steelman — why it *should* work (Bazley–Kato)
For a genuine two-body problem with **asymptotic** constituents of mass `m_g` and a *single* bound state, the level above the ground state is the **kinematic** ionization continuum at `s_th = (2m_g)²`. Then Temple's `ρ = (2m_g)²` is *free* (kinematic), the bound-state energy lies strictly below it, and one gets a rigorous binding-energy bound. This is exactly how rigorous He/helium bounds work (`BS-VARIATIONAL-LOWER §2`, Bazley–Kato 1961). If pure YM behaved like this, `E₂` would be kinematic and the circularity would break.

### 2.2 Closure α — the next `0⁺⁺` level is DYNAMICAL, not kinematic
In the physical `0⁺⁺` channel the level above the lightest glueball is the **first radial excitation `0⁺⁺*`** (≈`1.5–1.6 m₁`; `PROOF-ISOLATION §4.4`, `G2D-ISOLATION-I2 §1`), which sits **below** the two-glueball cut `(2m₁)²`. It is a dynamical Bethe–Salpeter bound state, not a continuum threshold. So `E₂ = m²_{0⁺⁺*} = I2`. The constituent "two-gluon threshold" `4m_g²` is (i) **not physical** (single gluons are confined, not asymptotic), (ii) sits **above** the true `E₂` when `m_{0⁺⁺*} < 2m_g` ⇒ **invalid `ρ`** ⇒ **overclaims** (§0.5 witness; `BS-VARIATIONAL-LOWER §3,§5`), and (iii) even granting `m_g`, the clean sub-threshold window `[m₁², 4m_g²)` is **empty** on lattice/DSE numbers. **Verdict α: `E₂` is dynamical (`=I2`); the kinematic steelman does not instantiate.** [T about the level structure, given both states bound]

### 2.3 Closure β — Single-Scale Covariance makes any threshold `∝ m₁`
Even if a genuinely kinematic cut existed in some channel, L3 (§0.4) says any threshold is a fixed multiple of the single scale: a multi-particle threshold is `(n·m_lightest)²` with `m_lightest` the target itself or a fixed multiple; a "different channel" mass is `c·Λ_QCD = (c/m₁)·m₁²` — again `∝ m₁²`. **To escape SSC one needs a threshold set by a *different* scale** (topological, order-parameter) — and pure YM has none (`G2D-ISOLATION-I2 §2`; confinement is *not* SSB of a local symmetry, so no local order parameter; the only other scale candidate, the string tension, is **1-form** — §4). So the sector/threshold mismatch cannot manufacture an independent `E₂`: it either reproduces `∝m₁²` or imports a scale the theory lacks. **Verdict β: SSC is a dimensional closure, not a missing theorem.** [T-conditional on one-scale + [J]]

### 2.4 The genuine non-collapsing corner — vNW *pointwise* isolation
**This is the one place the circularity is genuinely broken, and it deserves the credit it gets** (`G2D-ISOLATION-I2 §3`, `PROOF-ISOLATION §4.4`). Ground and excited `0⁺⁺` are both `A1g` (same symmetry sector). **von Neumann–Wigner** forbids the crossing of two same-sector levels. Strong-coupling cost grading puts `m₁ < m_{0⁺⁺*}` and the level difference is sign-constant on the connected β-line (no zero-`T` transition; confinement is a crossover). Hence `m₁ < m_{0⁺⁺*}` **pointwise at every β where both are bound**, and the lightest `0⁺⁺` **is** isolated below the excited cut — **without presupposing `m₁>0`**, i.e. non-circularly (it presupposes only the strong-coupling *ordering*, a different input).
**Why it is not Clay:** (i) the excited `0⁺⁺*` **dissolves into the continuum** at large β (Kato isolation fails *there*), so vNW protects only the strong-coupling side; (ii) it gives strict inequality (`ρ>1`), never a **uniform `κ`**; (iii) the margin `ρ−1` is set by the unresolved `O(u)` decoration (`SUB2`) — the same estimate as `H*`. The only scale-**invariant** escape is the **ratio `ρ = m_{0⁺⁺*}/m₁ ≥ 1+c`** itself (§4), unproven. **Verdict: pointwise break = real partial; uniform break = none.**

### 2.5 Bottom line on (d)
The `LEVER-nonconvex-edge` collapse is **truly closed for the uniform statement**, and for a sharper reason than "our tools are weak": the next level is dynamical (`α`) *and* tied to the answer by one scale (`β`). The genuinely non-collapsing corner (`vNW pointwise`) is real, is [T]-shaped, but is **sub-uniform**, so it does not close Clay. This is the precise sense in which "the reverse tool's input = the answer" survives.

---

## §3 Why "irreducible" is structural, not merely a toolbox gap — and what it does NOT prove

- **The closure is dimensional, not a missing lemma.** SSC (§0.4) says: in a one-scale theory, the isolation margin is a fixed multiple of the mass it is supposed to bound. This is not "we haven't found the theorem"; it is "there is no independent scale for `E₂` to live on." The corresponding statement for *tools* — every Family-I functional is Schur-concave, every Family-II input is `E₂`/`I2`/a certified trial — is the [NG]-for-classes content of `SCHUR-CONVEX` and this document.
- **But it is NOT `∀possible`.** Per `CONVERGENCE-redteam §d` and `KEY-EXISTS §5`, what is proven is that every **named** positive-measure/moment tool is one-sided and every **named** reverse tool inputs the answer. A not-yet-invented **non-positivity, non-regress, sign-controlling** inequality could evade both — it would have to be acid-test-clean (fail the free massless field) and would have to break either **SSC** (supply a genuinely second scale — impossible in pure YM without a new structure) or **L1** (pin an inf-support from trial data — impossible as a *theorem about operators*, unless the "trial" data is richer than `(R, Var)`, e.g. a certified near-eigenvector = the answer). So the irreducible verdict is **robust within the current architecture**, and the escape would have to be genuinely new (charter §1's "6th species"), **not** a relabeling.
- **The precise honest boundary:** "IRREDUCIBLE" here means *"no non-circular anchor exists in the identified space, and the reason is structural (one-scale + atom-only)."* It does **not** mean *"no proof of Clay exists"* (`KEY-EXISTS §1,§7`: proof in principle `~85–90% [J]`). The two must not be conflated.

---

## §4 The new-anchor sweep (g) — strongest untried shapes, and why none breaks the regress

The charter (§1) forbids stopping at "it's hard," so here is the honest attempt to *build* the missing anchor. Each candidate is judged by the acid test (must give **no** lower bound on the free massless field) and by whether it merely **relocates** the circularity.

1. **Scale-invariant ratio `ρ = m_{0⁺⁺*}/m₁` (the closest to real).** β-,`L`-invariant, acid-test-clean (the free field has no `0⁺⁺` tower to order, so the tool is vacuous there — correct). If one could prove **`ρ ≥ 1+c`** (e.g. **Landau-ladder rigidity** `m_n²=(2n+1)m₁²` ⇒ `ρ=√3`, a *symmetry-fixed* ratio), then + the proven strong-coupling gap (`SZZ`) gives uniform `I2` **without** touching `a→0`. **But:** the cost grading fixes only a *weight* ratio `(p'−p)ln(1/u)/a`; resolving it against the `O(u)` self-energy is **verbatim** the `SUB2` failure (`Cu` vs `Be^{−δt}`) and the `I1` residue-sign wall `H*` (`G2D-ISOLATION-I2 §5`). **Label [C]/conjecture, not theorem.** This is the one candidate that would *flip* the verdict if proven; it is recorded as the sharpest open sub-target, not as an anchor that already works.
2. **A threshold from a different (topological / order-parameter) scale.** Pure YM has **no local order parameter** (confinement is not SSB of a local symmetry; the Polyakov loop is non-local and vanishes at `T=0`), so no second scale is available (`G2D §2`). The only other scale, `√σ` (string tension), is **1-form**; `m_{0⁺⁺}=c√σ` is a *measured* transmutation, not a theorem, and `area law ⇏ gap` — cross-sector (`BS-THRESHOLD §4`). **[NG] for this shape.**
3. **'t Hooft-anomaly / symmetry forcing.** 1-form `Z_N` anomaly is blind to the local `tr F²` residue (a point carries no 1-form charge); anomaly matching *forbids* symmetric trivially-gapped vacua, the **opposite** of forcing a gap, and at θ=0 the electric anomaly is trivial (no force); θ=π forces at most a gapped TQFT (`PROOF-1FORM-ANOMALY`). 4D `a`-theorem is a one-sided *upper* bound on a dimensionless monotone, admits both gapped (`a_IR=0`) and gapless CFT, so it cannot exclude the CFT (`A-THEOREM`). **[NG].**
4. **Bootstrap / crossing / S-matrix.** Conformal bootstrap: YM is not a CFT (`T^μ_μ≠0`), free Maxwell is crossing-exact and unexcludable (`FORGE5-BOOTSTRAP-EXCLUDE` Lemma), and the output is a one-sided `Δ`-island (wrong coordinate/direction) (`GLUEBALL-BOOTSTRAP`, `A4-SPECTRAL-BOOTSTRAP`). S-matrix bootstrap **needs asymptotic atoms** = `I1∧I2` = `H*`, strictly **stronger** than the gap; Froissart/Regge give only `α(0)≤1` (critical Pomeron `m₀²=0` is consistent) (`PROOF-REGGE-FROISSART`). **[NG].**
5. **Krein/Pontryagin spectral function on the indefinite Savvidy Hessian.** The right language for the non-Hilbert operator, but its invariants are **counts** (negative index, signature) or **UV scales** (`g|B|` Landau bottom); the gap is the `k→0` **min**, and the definitizing split loses its uniform margin exactly at `k=0` (`FORGE-KREIN-SPECTRAL`). **[NG].**
6. **Index / spectral flow.** A discrete invariant is rigid on the connected β-line ⇒ equals its free value `0` (`INVENT-NEW-INDEX` rigidity no-go); the only acid-test-passing index counts **gauge copies**, not physical states. **[NG].**
7. **Entanglement / area law; large-N; holography; `θ`-vacuum convexity; small-circle/interface.** Area law is one-directional in `≥2D` (gapless area-law states exist) and is 1-form/line data; large-N/holography are not rigorous 4D anchors; `d²E/dθ²=χ_t>0` caps `m_{0⁻⁺}` **from above**; the small-circle deformation's only free inequality `m(S¹_r) ≥ m(ℝ⁴)` is an **upper** bound (wrong way) (`SMALLCIRCLE-ADIA`, `W3-adjacent-4D`). **[NG].**

**Verdict (g): no unconsidered anchor breaks the regress; the ratio `ρ` is the single sharpest untried, and it relocates onto `O(u)=H*`.** [O]/[J]

---

## §5 The exact missing step(s), precisely stated

Everything collapses to **one** of two *equivalent* statements (by L2), uniform in `(β>β_c, L, a)`:

```
(I2, isolation)   ∃δ>0:  spec(H)|_{0⁺⁺} ∩ (m₁, m₁+δ) = ∅ ,   i.e.  E₂ ≥ E₁ + δ
                  — equivalently, with SSC,  m₂² − m₁² ≥ (ρ²−1)m₁²  = uniform m₁ ≥ μ.

(T-ii, trial)     ∃ certified ψ₁:  Var_⊥(ψ₁) ≤ (R_⊥ − E₀ − λ)²     (a known glueball)
                  — equivalently  f₁² ≥ z₀²  (the (M1) overlap), i.e. the atom.

(By L2 both are equivalent to  f₁² ≥ r₀  (I1)  and to "the gap.")
```

So **support (`E₁≥Δ`) / isolation (`I2`) / residue (`f₁²≥r₀`) are one wall** (`SUPPORT-NOT-ATOM`, `GAP-G2-edge`, `PROOF-ISOLATION`). The genuinely *new* (non-species) object, if one insists on naming a direction the corpus has not closed, is the **scale-invariant ratio**
```
ρ := m_{0⁺⁺*}/m₁ ≥ 1 + c   (β,L-uniform, c>0),
```
which is **not** a restatement of the gap (it is `gap ∧ [ratio]`, strictly stronger, `G2D-ISOLATION-I2 §2`), is acid-test-clean (free field: no tower), and — **if proven** — would close uniform `I2` from the *proven* strong-coupling gap without `a→0`. **But proving `ρ ≥ 1+c` rigorously is the same `O(u)`-decoration sign estimate that is `H*`.** [O]

---

## §6 Honest boundary — what this document does NOT say

- It does **not** prove Clay, and does **not** prove that no proof exists. "IRREDUCIBLE" is scoped to *"no non-circular anchor in the identified Family-I/Family-II space, with a structural (one-scale, atom-only) reason"* — a **[NG]-for-classes + [O]-exact-missing-step + [J]**, never a global barrier (`∀class ≠ ∀possible`).
- It does **not** under-sell the pointwise **vNW** corner: it is a genuine, non-circular break of the circularity *at the pointwise level* ([T]-shape). It simply does not furnish the **uniform** margin Clay requires.
- The **SSC `E₂=ρE₁`** equivalence is **[T]-conditional**: the `ρ>1` half is vNW ([T]); the `ρ` bounded/β-independent half is the one-scale hypothesis ([J], empirical `ρ≈1.5–1.6`). If `ρ` were *unbounded*, uniform isolation would not imply the gap — but pure YM gives no mechanism to unbounded `ρ` (single scale; the only scale-invariant escape is the ratio itself, §5).
- The fresh computation (§0.5) is a two-line witness of *already-banked* lemmas; the substantive numerics (`Temple 0/3474`, commutator `5000/5000`, naive overclaim `5/5`, impossibility lemma exact) are **inherited [R-internal]**. I re-ran **no** external source and **no** heavy lattice job this session.
- **Flip condition:** exhibit **one** functional `Φ` on positive KL measures that is (i) non-Schur-concave (not an upper bound), (ii) computed from trial data richer than `(R,Var)` *without* a certified near-eigenvector, and (iii) fails the free massless field — then this verdict flips to **BROKEN**. None is known; building one is the "6th species" of `CLAUDE.md §1`.

---

## §7 References

**Internal (absolute paths; the banked verdicts this document classifies against):**
- `D:\Dev\math\CLAUDE.md`; `D:\Dev\math\docs\BRIEF-FOR-MODEL.md`; `D:\Dev\math\research\unconditional\HISTORICAL-TOOL-SWEEP.md` (§2–3, the two walls); `D:\Dev\math\research\unconditional\KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 the missed lever, §7 calibration); `D:\Dev\math\research\unconditional\CONVERGENCE-redteam.md` (§b3 non-convex edge, §d `∀class≠∀possible`).
- `LEVER-nonconvex-edge.md` (Temple/determinant/Combes–Thomas/sum-rule all collapse); `TEMPLE-WEINSTEIN.md` (naive `E₀`-shift false; L1 impossibility lemma; `E₁⟺E₂=I2`); `TEMPLE-RIGOROUS-H.md` (rigorous `H` gap-agnostic; kinematic-edge steelman killed); `PROOF-EDGE-E2.md` (regress one rung up; bottom-functional blindness); `BS-VARIATIONAL-LOWER.md` (next `0⁺⁺` level = radial excitation, not kinematic); `BS-THRESHOLD.md` (two thresholds; string tension 1-form, cross-sector); `BS-RP-STRUCTURE.md` (composite RP still upper-bound; Källén threshold `(2m_g)²` conditional).
- `GAP-G2-edge.md` (moment ratio upper; residue–level-gap lemma); `GAP-G2-mass.md` (edge/localization; transverse vertex unprotected); `PROOF-ISOLATION.md` (L1=I1∧I2; operative `Λ≈1.5m_g`); `G2D-ISOLATION-I2.md` (Single-Scale Covariance; vNW pointwise; ratio `ρ`); `G2F-ATOM-PLUS-ISOLATION.md` (determinant zero = atom+isolation, edge-locked); `SUPPORT-NOT-ATOM.md` (support vs atom; Bushell rate `=G″=H1`); `TRANSPORT-MACHINE.md` (Stage 3 = H1/I2, one scale).
- `SCHUR-CONVEX.md` (master Schur no-go); `GOLDS-POS-MASS.md` (positivity pins only `m=0`); `A-THEOREM.md` (4D `a`-theorem upper, CFT not excluded); `PROOF-1FORM-ANOMALY.md` (1-form anomaly blind to local `tr F²`); `FORGE-KREIN-SPECTRAL.md` (Krein count/UV, not `k→0` min); `INVENT-NEW-INDEX.md` (rigidity no-go); `FORGE5-BOOTSTRAP-EXCLUDE.md`, `GLUEBALL-BOOTSTRAP.md`, `A4-SPECTRAL-BOOTSTRAP.md` (bootstrap upper/feasible, free CFT crossing-exact); `PROOF-REGGE-FROISSART.md` (S-matrix needs atoms, upper/free-mass input); `PROOF-HAMILTONIAN-VARIATIONAL.md` (threshold from above); `P4-RP-uniform-gap.md` (RP no gap; finite-block criterion needs uniform block gap); `PROOF-SPECTRAL-INDEPENDENCE.md` (sign re-imported at the gap step); `W3-adjacent-4D.md` (marginal, scale-free, must transmute).
- `INEQ-R6-4pt-ward.md` (4-pt protected contact homogeneous `f₁²·c`/disconnected); `INEQ-T4`, `INEQ-T5`, `INEQ-T1`, `INEQ-R10`, `GAP-IN11`, `INEQ-R3` (referenced within the above).

**External (carried from the repo; **[R]**, not re-verified this session):**
- Temple, *Proc. London Math. Soc.* (1928); Weinstein, *Phys. Rev.* **40** (1932) 797; Kato, *J. Phys. Soc. Japan* **4** (1949) 334; Courant–Fischer; Bazley–Kato, *J. Math. Mech.* **10** (1961) 291; von Neumann–Wigner (1929); Reed–Simon IV Thm XIII.5.
- Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Osterwalder–Schrader, *CMP* **31** (1973) 83; Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737.
- Gaiotto–Kapustin–Seiberg–Willett, arXiv:1412.5148; Cordova–Ohmori, arXiv:1910.04962; Cordova–Dumitrescu–Intriligator, arXiv:1802.04790; Komargodski–Schwimmer, arXiv:1107.3987.
- Froissart, *Phys. Rev.* **123** (1961) 1053; A. Martin (1963/1969); Chew–Frautschi, *PRL* **7** (1961) 394; Mandelstam (1958–59); Steinmann, LNP **11** (1971); Lehmann–Symanzik–Zimmermann (1955); Haag–Ruelle (1958/1962).
- Morningstar–Peardon, *PRD* **60** (1999) 034509; Källén (1952); Lehmann (1954); Caron-Huot–Pokraka–Zahraee, *JHEP* 01 (2024) 195, arXiv:2309.04472; Bushell, *Arch. Rational Mech. Anal.* **52** (1973) 330; Birkhoff, *Trans. AMS* **85** (1957) 219; von Renesse–Sturm, *CPAM* **58** (2005) 923; Anari et al. (spectral independence), referenced within `PROOF-SPECTRAL-INDEPENDENCE`.

---

*Status: **[NG] for the named Family-I (Schur-concave) and Family-II (regress/rigidity) classes + [O] exact missing step (uniform I2 / certified trial / scale-invariant ratio `ρ`) + [T]-conditional structural reason (Single-Scale Covariance) + [J] (no third anchor in the identified space).** Verdict: **IRREDUCIBLE for the uniform bound, with one genuine non-collapsing corner (vNW pointwise isolation).** NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay stays open.*

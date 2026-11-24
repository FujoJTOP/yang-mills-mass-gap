# S6-DOBRUSHIN — Dobrushin / disagreement-percolation / decoupling attack on the live inequality (★)

> 2026-11-24 · Sixth attack wave, angle **S6**: treat the conditional block measure as a perturbation of
> the free block by the boundary, and bound the boundary influence with **Dobrushin's coefficient**,
> **disagreement percolation**, and the **Simon–Lieb / Brascamp–Lieb** martingale route.
> Reads: `research/v3/C2-BALABAN-UNIFORM.md` §2.4, `research/unconditional/PROOF-COARSE-POINCARE.md`,
> `research/unconditional/PROOF-JACOBIAN-SAVES.md`, `G1-L3.md`, `G1-L3-TOWER-POINCARE.md`,
> `G1-L3-UV-FORCE-BOUND.md`, `docs/paper-v2/ch05-theorems.tex` (Lemma SC),
> `research/unconditional/ATOM-strong-coupling-extend.md` (Hooft-normalisation bookkeeping).
> Witness: `research/v3/w7/s6_dobrushin.py` (run output quoted in §7).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[R]** reported/unrefereed · **[V]** citation verified this session.

---

## 0. Verdict (no QED)

**The influence/decoupling attack does not close (★) on the physical tower — but it *splits* the tower
cleanly into two legs that each carry a theorem, with the head as the single remaining wall.**

1. **Dobrushin / disagreement percolation is a strong-coupling tool, and it caps out at the same place as
   Lemma SC.** The Dobrushin coefficient of the shell-`j` block fibre is `α(β_j) = κ β_j` with `κ = O(d)`
   (this session: `κ≈72`, `β_Dobr = 1/72` in the standard-Wilson normalisation). When `α<1` the
   Dobrushin–Shlosman complete-analyticity machinery gives both the one-sided bound **(★)** and the band
   decay **(★dec)**, `|Cov(s_B,s_{B'})| ≤ C e^{−κ dist/a}` — a **genuine [T] at β_j < β_Dobr**. But
   `β_Dobr = O(1/48)` (same family as SZZ, equivalently Hooft `λ≳96`, standard-Wilson `β_W≈0.19`), while
   **every shell of the physical tower lives at `β_j ∈ [2.5,∞)`** — the strong-coupling window `[0,β_Dobr]`
   and the tower `[2.5,∞)` are **disjoint**. So the Dobrushin leg proves (★) for **zero** shells of the
   actual tower. **[T] statement, [O] relevance to L3.**

2. **The Simon–Lieb / Brascamp–Lieb leg works at the *other* end (deep UV), conditional on the cone
   bootstrap.** Brascamp–Lieb applied to the flat-chart fibre measure gives `Var(force) ≤ E[Hessian]`
   *exactly* whenever the Wilson fibre Hessian is `⪰0` and `Φ⪰0`: this is the equipartition identity,
   saturated by Gaussians. Then (★) reduces to the Haar budget `W_j ≥ δ_W`, which is **Lemma B
   (`W_geo ≥ N/12`, scale-free)**. This leg is **[C]** precisely on the missing cone/moment bootstrap
   `E|F_j| ≤ C g_j` that certifies `H_j ⪰ cΛ_j² > 0` in the deep UV (`G1-L3-UV-FORCE-BOUND` Route B).

3. **At the head `β≈2.5` both legs fail, for two independent, named reasons.** (i) Dobrushin:
   `α(2.5) ≈ 180 ≫ 1`. (ii) Brascamp–Lieb: the Wilson fibre Hessian is **indefinite** (Savvidy, min
   `≈−0.3`), so the measure is not log-concave and the equipartition `Var ≤ E[H]` fails — the defect
   `Var − E[H] > 0` now has to be *covered* by `W_j` plus margin, which is exactly the campaign's **(◆)/H1**
   object. This is **[O]**.

4. **The free field is the correct acid test and (★) fails there, by design.** On the free massless
   Gaussian block, `H = Σ^{−1}` gives `Var(G_v) = E[H_vv]` **identically** and `W=0`, hence `G″=0` and
   (★) is false. The strict margin `−δ_W` can only come from the compact-group Haar curvature `W_geo`
   (Prékopa: `≥0`; witness `≈0.205`), which the free field does not have. **[T].**
   *(This does not contradict `C2-BALABAN-UNIFORM.md` §5, which is about a different object — the
   `sup_a Z_a^{-1}‖H_a‖_D` UV-stability bound, which must PASS on the free field; (★) is the strict
   gap-type bound, which must FAIL. See §6.)*

**Net.** The attack does **not** prove `(★)` uniformly in `j,U′,a` — it is **[O]** on the physical tower,
with the head `β≈2.5` as the exact missing wall. What it *does* deliver: a rigorous `(★)∧(★dec)` in the
strong-coupling window (**[T]**, off-tower), a rigorous mechanism-level `(★)` in the deep UV (**[C]** on the
cone bootstrap), and a **method-class no-go [NG]** — *no purely local influence-contraction estimate reaches
`β=2.5`* — which localizes the obstruction to the Savvidy/log-concavity wall rather than the combinatorics.
Clay 0 %.

---

## 1. Setup: the block measure, the exact identity, and what (★) is

**Shell, block, fibre.** At shell `j` (block spacing `a_j`), fix a block `B` and its interior links
`Λ_B`. The **conditional block measure** given the boundary configuration `U′ = U_{Λ\Λ_B}` is
```
μ_j(·|U′) ∝  exp(−S_j(U_B ; U′)) ∏_{ℓ∈Λ_B} dU_ℓ ,   S_j = Wilson action = β_j Σ_{p ∩ B≠∅} (1 − (1/N)Re tr U_p),
```
with `β_j = 2N/g_j²` (standard Wilson normalisation; `β≈2.5` at the head, `β_j→∞` deep in the UV, one-loop
`β_j = 4Nβ₀ j`). `P` is the block-averaging (coarse-graining) map, `Pᵀ` its adjoint; `H_j = Hess S_j` is the
**flat (Wilson) Hessian**; `W_j` is the **Haar/Jacobian curvature** budget of `PROOF-JACOBIAN-SAVES`. `a` is
the UV cutoff.

**The exact identity (the whole content of (★)).** For the marginal effective action
`S_eff(U′) = −log ∫ e^{−S_j(U_B;U′)}∏dU_B`, the standard log-partition Hessian identity gives, along a
unit collective mode `v`,
```
G″_j(v) := ∂²_v S_eff = E_{μ_j(·|U′)}[⟨PᵀH_jP v,v⟩]  −  Var_{μ_j(·|U′)}(⟨Pᵀ∇S_j, v⟩) .   (1.1)
```
Adding the Haar part in the flat chart (`X=log U`, measure `e^{−S_j−Φ}dX`, `Φ=−log J`), (1.1) becomes the
campaign's split (`PROOF-JACOBIAN-SAVES` §1)
```
G″_j = Δ_dyn + W_j + X ,   Δ_dyn = E[H_j] − Var(G_v) ,   W_j = E[Φ″] − Var(Φ′) ,   X = −2Cov(G_v,Φ′),
```
so the inequality (★), `Var(G_v) ≤ E[H_j] + W_j − δ_W`, is **exactly** `G″_j ≥ δ_W` (1.1). It is a statement
of **strict positivity of the marginal curvature**, not a variance-concentration statement per se.

**(★) is gap-type, not UV-type.** Because of the `−δ_W` (strict), `G″_j = 0` on the flat free field (§6)
violates it; `G″_j = W_geo > 0` on the compact free Haar measure (β=0) satisfies it. So (★) separates
compact Yang–Mills from the free Gaussian — the correct charter acid test. (The *passing* free-field test of
`C2-BALABAN-UNIFORM.md` §5 is a different object, §6.)

---

## 2. Dobrushin's coefficient for the block fibre

**Single-link conditional.** Given all other interior links and the boundary `U′`,
```
μ_{j,ℓ}(dU | U_{B\ℓ}, U′) ∝ exp( (β_j/N) Re tr( U A_ℓ(U_{B\ℓ},U′)) ) dU ,   A_ℓ = Σ_{p∋ℓ} (product of the other 3 links of p).
```

**Coefficient.** `c_{ℓℓ'} := sup_{η,η̃ agree off ℓ'} ‖μ_{j,ℓ}(·|η) − μ_{j,ℓ}(·|η̃)‖_TV`.
- If `ℓ,ℓ'` share no plaquette: `c_{ℓℓ'}=0`.
- If they share a plaquette (`ℓ` sits in `2(d−1)` plaquettes; each contributes 3 links `ℓ'`, so `≤ 6(d−1)`
  such `ℓ'`, with `≤2` shared plaquettes), the log-density ratio is `(β_j/N)(Re tr U A_ℓ − Re tr U A_ℓ')`,
  bounded by `2β_j` per shared plaquette; hence
```
c_{ℓℓ'} ≤ κ_0 β_j ,   κ_0 = O(1)   (this session: κ_0 = 4, conservative).   (2.1)
```

**Dobrushin coefficient and threshold.**
```
α(β_j) = sup_ℓ Σ_{ℓ'≠ℓ} c_{ℓℓ'}  ≤  6(d−1)·κ_0·β_j  = 72 β_j  (d=4, κ_0=4) .
⟹  α(β_j) < 1  ⟺  β_j < β_Dobr := 1/(6(d−1)κ_0) = O(1/d) = O(1/48).   (2.2)
```

**Normalisation-honest comparison.** *Everything below is in the standard-Wilson normalisation of the head
`β≈2.5`.* The crude `β_Dobr = 1/72 ≈ 0.014`; the *sharpened* member of the same strong-coupling family —
Shen–Zhu–Zhu's LSI/uniqueness at "`β<1/48`" in their own normalisation — corresponds to Hooft
`λ = g²N ≳ 96`, i.e. standard-Wilson `β_W ≈ 0.19` (`ATOM-strong-coupling-extend.md`). So:

| quantity | standard-Wilson `β` | status |
|---|---|---|
| crude Dobrushin `β_Dobr` | `≈1/72` | [T] (this computation, §2) |
| sharp strong-coupling (SZZ) | `≈0.19` (Hooft `λ≈96`) | [T] (SZZ, CMP 400) |
| physical head | `2.5` | target |
| **honest distance head / threshold** | **`≈13×` in Hooft `λ`** (not the naive `1/48` vs `2.5`) | [C] (bookkeeping) |

`α(2.5) = 180 ≫ 1`, `α(0.19)=13.7>1`, `α(1/72)=1`. **The Dobrushin window `(0,1/72]` and the tower
`[2.5,∞)` are disjoint.**

**What `α<1` buys (Dobrushin–Shlosman, [T]).** When `α(β_j)<1`, the fibre measure is *completely analytic*:
(i) the boundary influence on any fibre link decays geometrically,
`infl(ℓ;∂) ≤ C α^{dist(ℓ,∂)/a}`; (ii) correlations cluster exponentially, giving **(★dec)** verbatim with
`κ = −log α > 0`; (iii) the free energy's boundary curvature `G″_j` is `C^∞` with a bound
`|G″_j − G″_j^{β=0}| ≤ C(β_j)` where `C(β_j)→0` as `β_j→0`. Combined with the exact
`G″_j|_{β=0}=W_geo≥N/12` (§6), this gives
```
G″_j ≥ W_geo − C(β_j) ≥ δ_W > 0   for β_j < β_Dobr ,  uniformly in U′, a,   ⟹ (★)∧(★dec).   (2.3)
```
This is the **[T]** strong-coupling result. It is the same regime as Lemma SC / P3 and does **not** extend it.

---

## 3. Disagreement percolation and the "total influence" formulation

**Disagreement percolation (van den Berg–Maes).** Couple `μ_j(·|U′)` and `μ_j(·|U″)` (boundary conditions
differing somewhere) so that the *disagreement set* `D = {ℓ : U_ℓ ≠ U″_ℓ}` is a percolation cluster. The
per-bond disagreement probability is `q_ℓℓ' ≤ c_{ℓℓ'} = κ_0 β_j` (2.1). If the `q`-percolation is
**subcritical** (`Σ_{ℓ'} q_{ℓℓ'} < 1`, i.e. `α<1`), then
```
infl(ℓ;∂) ≤ P_q(ℓ ⟷ ∂ in D) ≤ C e^{−m dist(ℓ,∂)/a} ,    (3.1)
```
with `m > 0` from the subcritical exponential-decay rate. This is the sharpened version of §2: the *constant*
`κ_0` can be optimized (van den Berg–Maes beat the raw Dobrushin `α<1`), but the **threshold stays
`β = O(1/d) = O(1/48)`** — the class is unchanged. At `β=2.5` the per-bond `q = O(1) ≫ p_c(4D)`, the cluster
is supercritical, (3.1) is vacuous.

**Total-influence ⇒ (★) lemma (the explicit threshold the task asks for).** Define
```
I(U′) := Σ_{ℓ∈B} infl(ℓ;U′) ,   infl(ℓ;U′) := sup_{U′,Ũ′ differ at one boundary plaquette} ‖μ_{j,ℓ}(·|U′,U_{B\ℓ}) − μ_{j,ℓ}(·|Ũ′,U_{B\ℓ})‖_TV .
```
> **Lemma DL.** There is a model constant `C_I` (depending on `d,N` and the block geometry, not on `j,a`) such
> that
> ```
> | Var_{μ_j}(G_v) − E_{μ_j}[H_j,vv] |  ≤  C_I · I(U′) · osc(G_v)² .   (3.2)
> ```
> Consequently **(★) holds whenever the total influence is below the threshold**
> ```
> I(U′)  <  I_* := (W_j − δ_W) / ( C_I · osc(G_v)² ) .   (3.3)
> ```
> At strong coupling `I(U′) ≤ C_I' |B| α/(1−α)·e^{−m dist}` (3.1), so `I(U′)→0` uniformly in `U′` as
> `β_j→0`, and (3.3) is satisfied — recovering (2.3).

*Proof sketch.* (3.2) is the Efron–Stein / Simon–Lieb martingale bound: write
`Var(F) ≤ Σ_ℓ E Var(F | rest_ℓ)`, and `Var(F|rest_ℓ) ≤ ¼ osc_ℓ(F)² · (TV-sensitivity of the ℓ-marginal to
its conditioning)`, then sum. The RHS is `≤ C_I I(U′) osc²` by (3.1). With `Var − E[H] ≥ −C_I I osc²`,
inserting into `G″=E[H]−Var+W` gives `G″ ≥ W − C_I I osc²`, i.e. (★) iff (3.3). ∎
*(Honest caveat: I split `I` into a per-link TV-sensitivity; the sharp `C_I` is not optimised here. The
load-bearing content is the *scaling* `I ∝ β_j` and the threshold (3.3), not the constant.)*

**Why this leg is still only strong coupling.** `I(U′) ≤ C β_j` grows with `β_j`; the threshold (3.3)
holds only when `β_j` is small. At `β_j = 2.5` (head) and `β_j→∞` (UV) `I≫I_*`. **So the total-influence
criterion is, like Dobrushin, a small-`β` criterion — it cannot reach the tower.**

---

## 4. Range of validity, and what fails beyond

| regime | `β_j` (standard Wilson) | Dobrushin/influence | reason |
|---|---|---|---|
| strong coupling | `0 < β_j < β_Dobr ≈ 1/72` (≤ SZZ `≈0.19`) | **(★)∧(★dec) [T]** (§2.3, §3) | `α<1`, `I→0`, `W_geo≥δ_W` |
| intermediate | `β_Dobr ≤ β_j ≲ 1` | **vacuous [O]** | `α>1`, percolation supercritical |
| **head / crossover** | **`≈2.5`** | **fails [O]** | `α(2.5)≈180`; Savvidy indefinite (§5) |
| deep UV | `β_j → ∞` | Dobrushin vacuous; **BL leg [C]** | `α→∞`, but log-concave (§5) |

**Exact range where the influence method is a theorem:** `β_j < β_Dobr = O(1/48)` (standard-Wilson). This is
**disjoint from the physical tower `[2.5,∞)`**: the crude cutoff `1/72` is `×180` below `2.5`, the sharpened
SZZ cutoff `≈0.19` is `×13` below `2.5` in Hooft coupling. **No shell of the actual RG tower is covered.**

**What fails beyond `β_Dobr`, precisely:** (i) `α(β_j)≥1` and the disagreement cluster is supercritical, so
the *contraction* behind (2.3)/(3.1) is gone; (ii) more deeply, the mechanism changes — at `β_j ≥ O(1)` the
fibre is *not* a small perturbation of a high-temperature measure but an increasingly **sharp, curved** one,
whose correct linearisation is **log-concave (Brascamp–Lieb)** at the UV end and **Savvidy-indefinite** at the
head. The influence/contraction toolbox has no entry there.

---

## 5. The Simon–Lieb / Brascamp–Lieb leg: which term is `W_j`, which is `δ_W`, and the deep UV

**Exact identity + BL.** Take the flat-chart total exponent `A = S_j + Φ` (`Φ` the Haar log-Jacobian), measure
`μ ∝ e^{−A}`, force `G_v = ⟨Pᵀ∇S_j,v⟩`, `∇G_v = PᵀH_jPv`. **Brascamp–Lieb** (for `A` convex / `μ`
log-concave):
```
Var_{μ}(G_v) ≤ E⟨ PᵀH_jPv, (Hess A)^{−1} PᵀH_jPv ⟩ .   (5.1)
```
If `H_j ⪰ 0` on the fibre and `Φ″ ⪰ 0`, then `(Hess A)^{−1} ⪯ H_j^{−1}` (Loewner, on the relevant
subspace) and (5.1) collapses to the **exact equipartition**
```
Var_{μ}(G_v) ≤ E[H_j,vv]   (equality iff Gaussian).   (5.2)
```
Then (★) reduces to `W_j ≥ δ_W`, i.e. **Lemma B** (`W_geo ≥ N/12` pointwise, `W_j≥δ_W` scale-free,
`PROOF-JACOBIAN-SAVES` Lemma B, `G1-L3-TOWER-POINCARE` §4).

**Term identification (the task's §3 request), made explicit.**

| symbol | object | source | sign |
|---|---|---|---|
| `δ_W` | the safety **margin** `>0` required of `G″` | definition of (★) | strict |
| `W_j` | **Haar/Jacobian curvature** `E[Φ″]−Var(Φ′)` | Lemma B, `≥N/12` pointwise | `≥0` (Prékopa) |
| `E[H_j]` | Wilson flat Hessian expectation | (5.1) RHS upper edge | the "Gaussian" level |
| `Var(G_v)−E[H_j]` | **log-concavity defect** (Savvidy sign) | controlled by (5.1)/(3.2) | `>0` = Savvidy |

So in the BL/Simon–Lieb route: **`W_j` is the Haar budget that turns the saturated equality `Var=E[H]`
into the strict `Var ≤ E[H] − (positive)`, and `δ_W` is the fraction of that budget one keeps as margin.**
`W_j` and `δ_W` are *not* the same term; `W_j` is geometric (independent of `β`), `δ_W` is the chosen slack.

**Where (5.2) holds — deep UV [C].** `H_j ⪰ c Λ_j² > 0` **iff the cone is empty**, `|F_j| = O(g_j) ≪ Λ_{j+1}`
(`G1-L3-UV-FORCE-BOUND` Route B). There the Savvidy negative term `O(g_j)` cannot beat the Laplacian floor,
`H_j` is PSD, the fibre is log-concave, (5.2) holds, and `W_j≥δ_W` closes **(★)**. The exact missing input is
the **gauge-invariant, `a`-uniform moment/cone bootstrap `E|F_j| ≤ C g_j`** — **[O]** (`G1-L3-UV-FORCE-BOUND`
§6 item (1)). Granting it, the deep-UV half of the tower is **[C]**.

**Where (5.2) fails — the head [O].** At `β≈2.5`, `PᵀH_jP` has a negative mode (Savvidy, min `≈−0.3`,
`MARGINAL`), so `H_j` is **indefinite**, `(Hess A)^{−1} ⪯ H_j^{−1}` is false, and the log-concavity defect
`Var−E[H] > 0` must be paid out of `W_j`. The campaign's numbers put `W_full ≈ 0.22–0.31` against
`−(E[H]−Var)` and close numerically (`W7` +0.51 exact on the softest mode), but *certifying*
`W_j > (Var−E[H]) + δ_W` is exactly **(◆)=H1** (`PROOF-JACOBIAN-SAVES` §2). **[O].** Neither leg of this
attack touches it.

---

## 6. Free-field acid test: (★) must fail, and it does

**Free massless Gaussian block.** `S = ½xᵀHx`, `H ⪰ 0` the (free) Laplacian, measure `N(0,Σ)`, `Σ=H^{−1}`,
**no Haar** (`Φ=0`, `W=0`). Then for the force `G_v = ⟨Hx,v⟩`:
```
Var(G_v) = vᵀH Σ H v = vᵀH v = E[H_vv]   (exact equipartition),   and   W = 0,
⟹  G″ = E[H] − Var + W = 0   ⟹  (★) [ = G″≥δ_W>0 ]  FAILS.   ([T], witness Part B)
```
So **(★) fails on the free field**, as the charter requires. The witness also shows the boundary condition
**shifts the mean** of `x` (harmonic extension) but **never the variance** — the free-field defect is
mean-blind, so no boundary-influence bound can ever *make* the free field satisfy (★): the failure is the
absence of `W_geo`, not a boundary effect. On the **compact** group the same β=0 measure has
`G″ = W_geo > 0` (witness Part C: `V_eff″(0)≈0.205`, same order as the campaign's `0.248`), which is why the
*lattice gauge* theory passes at strong coupling while the free Gaussian does not.

**Reconciliation with `C2-BALABAN-UNIFORM.md` §5.** §5 distinguishes the *UV-stability* object
`sup_a Z_a^{-1}‖H_a‖_D < ∞` — which **must pass** on the free field (`Z=1`, free Laplacian norm `a`-uniform)
— from the **gap-type** object, which must fail. (★), with its strict `−δ_W` and its `G″=0` on the free
field, is the **gap-type** object. Hence §5's "the UV statement must not fail on the free field" is *not*
in tension with "★ must fail on the free field": they are different functionals. §6's failure is the correct,
required behaviour for (★).

---

## 7. Witness (`research/v3/w7/s6_dobrushin.py`, deterministic, `py -3`)

Quoted output (numpy 2.4; `np.trapezoid`):

* **Part A — Dobrushin threshold.** `α(β)=72β`, `β_Dobr=1/72≈0.0139`. `α(0.19)=13.7`, `α(2.5)=180`.
  Head/threshold distance `≈13×` in Hooft `λ` (SZZ `β_W≈0.19`). → the local-influence class is dead at 2.5.
* **Part B — free field.** Free `H=Σ^{−1}` Dirichlet chain `n=6`, uniform mode: `E[H]=0.333333`, `Var(G)=0.333333`,
  `defect=0`, `W=0`, **`G″=0`** → (★) fails. Boundary `b=3` moves the interior mean `(2.571,0.429)` but
  `Var(G)` is byte-for-byte unchanged → variance is mean-blind.
* **Part C — Haar budget.** One-link SU(2) pure Haar, flat chart, `V_eff″(0)=1/∫₀^{2π}θe^{−Φ}dθ =
  1/4.8753 = 0.2051 > 0`. Same order as the campaign's `0.248` (difference = flat-chart/convention
  normalisation). The Haar term is the **only** source of the strict margin, and it is positive.

*Honest scope of the witness:* Part A is the actual coefficient arithmetic (exact given `κ_0`); Parts B–C
are exact dense-linear-algebra / 1-D-quadrature checks of the two structural claims (free-field saturation,
Haar positivity). It is a **witness, not a proof**; the theorem content is §2 (DS) and §5 (BL), and the open
step is §5's head.

---

## 8. Status ledger

| item | statement | label |
|---|---|---|
| exact identity `G″=E[H]−Var(G)+W` | log-partition Hessian; the content of (★) | **[T]** |
| Dobrushin coefficient `α(β_j)=κβ_j`, `κ=O(d)` | computed on the SU(2) fibre | **[T]** |
| `α<1 ⟹ (★)∧(★dec)` at `β_j<β_Dobr=O(1/48)` | Dobrushin–Shlosman complete analyticity | **[T]** (given DS) |
| strong-coupling window vs physical tower | `[0,1/72]` vs `[2.5,∞)` **disjoint** | **[T]** (arithmetic) |
| total-influence threshold `I<I_*` | Lemma DL (3.2)–(3.3) | **[T]** (scaling) / **[C]** (constant) |
| BL equipartition `Var≤E[H]` for log-concave fibre | deep-UV `H_j⪰cΛ²` | **[C]** on cone bootstrap |
| deep-UV (★) | cone empty, `W_j≥δ_W` | **[C]** on `E|F_j|≤Cg_j` (`[O]`) |
| head `β≈2.5` | `α≈180` **and** Savvidy indefinite; `W>deficit+δ_W` unproven | **[O]** = **(◆)=H1** |
| "no local-influence / Dobrushin-type estimate reaches `β=2.5`" | `α(2.5)≫1`, percolation supercritical | **[NG]** for the method class |
| free field `G″=0`, (★) fails | exact equipartition, `W=0` | **[T]** |
| `W_geo≥0` (compact), `W_geo>0` on the link | Prékopa; witness 0.205 | **[T]** |
| **(★) uniform in `j,U′,a`** | the campaign's L3 | **[O]** — not closed by this attack |

---

## 9. Verification ledger

| item | checked | status |
|---|---|---|
| Dobrushin (1970), conditional distributions / coefficient | standard | [standard, not re-verified] |
| Dobrushin–Shlosman, complete analyticity | standard | [standard] |
| van den Berg–Maes, disagreement percolation | standard | [standard] |
| Simon–Lieb / Efron–Stein martingale variance | standard | [standard] |
| Brascamp–Lieb, J. Funct. Anal. **22** (1976) 366 | carried from `A5-RIGOROUS-NUMERICS.md` refs | **[V]** (carried) |
| Holley–Stroock, Ann. Probab. **15** (1987) 1541 | carried from campaign | **[V]** (carried) |
| Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; `β<1/48` ⟺ Hooft `λ>96` | `ATOM-strong-coupling-extend.md` §1, `G1-L3` §3 | **[V]** (carried) |
| Hooft bookkeeping `1/48 ↔ β_W≈0.19 ↔ 2.47`; honest `≈13×` | `ATOM-strong-coupling-extend.md` §1 | **[V]** (carried) |
| Prékopa 1971; Lemma B `W_geo≥N/12` | `PROOF-JACOBIAN-SAVES` §2 | **[T]** (carried) |

**Notes.** (i) All external "standard" results are cited by author/name and used as standard textbook tools
(Dobrushin uniqueness, DS complete analyticity, disagreement percolation, BL/Efron–Stein); I did **not**
download sources this session and mark them `[standard]`, not `[V]`. (ii) The `1/48`–`2.5` comparison is done
in the **standard-Wilson** normalisation throughout, with the Hooft-normalisation caveat of
`ATOM-strong-coupling-extend.md` (naive `1/48` vs `2.5` overstates the distance by the `N²` in the Hooft
convention; the honest distance is `≈13×` in `λ`). (iii) Part C's `0.205` vs the campaign's `0.248` is a
flat-chart/normalisation difference, not a discrepancy in sign or mechanism.

---

*Refs.* `C2-BALABAN-UNIFORM.md` §2.4; `PROOF-COARSE-POINCARE.md` (c); `PROOF-JACOBIAN-SAVES.md` §1–2;
`G1-L3.md` §2; `G1-L3-TOWER-POINCARE.md` §2–5; `G1-L3-UV-FORCE-BOUND.md` §1–6; `G1-UNIFORM-ZA.md` §1;
`ATOM-strong-coupling-extend.md` §1; `docs/paper-v2/ch05-theorems.tex` (Lemma SC); Dobrushin (1970);
Dobrushin–Shlosman (complete analyticity); van den Berg–Maes (disagreement percolation); Simon–Lieb /
Efron–Stein; Brascamp–Lieb, J. Funct. Anal. **22** (1976) 366; Holley–Stroock, Ann. Probab. **15** (1987)
1541; Prékopa, Acta Sci. Math. (Szeged) **32** (1971) 301; Shen–Zhu–Zhu, CMP **400** (2023) 805,
arXiv:2204.12737; witness `research/v3/w7/s6_dobrushin.py`.

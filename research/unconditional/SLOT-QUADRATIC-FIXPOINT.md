# SLOT-QUADRATIC-FIXPOINT — the "sign-definite quadratic minorant" as a **FIXPOINT**: `Q = 𝒬[Hess, μ]`

> 2026-11-24 · Direct constructive attempt on the ONE open slot (`PROFILE.md §3`, `PROFILE-FACET-ANCHOR §3.2`): a `k→0` sign-definite gauge-invariant quadratic minorant `Q` (`Q ≤ Hess`, `Q ⪰ 0`) that is **supplied by the interaction** rather than **assumed**. The route: make `Q` the **solution of an extremal relation** `Q = 𝒬[Hess, μ]` over the measure, so `Q=0` is a *fixed point* to be tested, not a datum to be granted.
>
> **Verdict up front.** The fixpoint can be written and **exists**, but it **does not fill the slot.** Three tiers, exact:
> **(1) `[T]` Pointwise** (`Q⪰0`, `Q ≤ Hess` as forms): **empty** on the Savvidy shell — `Q⪰0` & `Q≤H` ⇒ `H⪰0`, so a PSD minorant exists only where the Hessian is already PSD = where the gap is assumed. *(the trap)*
> **(2) `[T]` Measure-averaged** (`Q ⪰ 0`, `Q ≤ E_μ[Hess]`): **exists** and is `E_μ[Hess]=E_μ[∇V⊗∇V]⪰0` (Stein), sign-definite *without* presupposing the gap — **but gap-blind**: the free massless field gives `Q=0` exactly. This is the only construction that is genuinely non-circular, and it certifies nothing at `k→0`.
> **(3) `[O]` Conditioned / effective-potential** (`Q := Hess S_eff = G″(k→0)`): the object where the gap lives. `Q=0` is admissible (free field), there is **no source at `Q=0`**, and `Q>0` uniformly = `H1 = I2`. The Savvidy indefiniteness is **compensated only at the gap-blind expectation level (Stein)**; at the conditioned level it is **adversarial** (`Q*≤0`), and overcoming it needs a positive source = the gap.
>
> **The single exact failing step (`[O]`):** the fixpoint relation is **source-free at `Q=0`** — `Q=0 ⟺ ∂_vV≡0 μ-a.s.`, which the free/gapless field satisfies — and pure YM is classically scale-free, so the only explicit scale-breaker (trace anomaly) enters **homogeneously and cancels**. A nonzero fixpoint therefore *is* spontaneous scale breaking (dimensional transmutation) *is* the mass gap. **The fixpoint does not exclude the gapless point; it labels it `Q=0`.** `Q=0` admissible, no source, uniform attainment = transport. **Clay remains 0%. No QED.**
>
> Labels (charter §2): **[T]** true theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[V-here]** this round · **[R]** inherited/external. Witness: `research/unconditional/w7/slot_quadratic_fixpoint.py` (`py -3.12`, exact linear algebra + 1D/2D quadrature only, no lattice MC).

---

## §0 · The object, and the three faces of "quadratic minorant"

**Setup (`H1-HESSIAN-DIRECT §1`, `PROOF-SAVVIDY-EXACT-SIGN §1`, `PROOF-JACOBIAN-SAVES §1`).** Lattice SU(N) Wilson, `μ=Z^{-1}e^{-S_W}dU`, `S_W=(1/2g²)Σ_□tr F²`. Flat chart `X_ℓ=log U_ℓ∈su(N)`, `dU=J dX`, `Φ=−log J`, so `μ∝e^{-V}dX`, `V=S_W+Φ`. The **Hessian** on the gauge-invariant soft `k→0` mode `v` (coarse block-uniform / zero-momentum transverse, `|v|=1`, `a=⟨v,X⟩`) is

```
   Hess V = −D² + 2 ad F_B + Φ″ ,       H_vv = vᵀ(Hess V)v ,
```
the three summands being (i) covariant Laplacian `−D²≥0`, (ii) spin–field `2adF_B` **indefinite**, Landau bottom `λ_min=−g|B|` `[T]` (Savvidy, nonabelian-only, physical after gauge quotient: `7.1%` neg, `min −0.82β`, W4), (iii) Haar Jacobian `Φ″⪰N/12>0` pointwise `[T]`.

**The "Dirichlet form" and the three faces.** The task's "quadratic minorant `Q ≤ Hess`, `Q>0`" is the **Poincaré form**: the coarse Dirichlet form `∫|∇φ|²dν` dominates the variance, `Var_ν(φ) ≤ (1/λ)∫|∇φ|²dν`, iff the effective-action Hessian `Hess S_eff ⪰ λ` on coarse modes (`PROOF-COARSE-POINCARÉ §0`). Equivalently `Q := G″(k→0) ≥ λ`, where `G″(0)=−(log p₀)″(0)` is the curvature of the coarse marginal `p₀` of `a_{k→0}`. So a "sign-definite quadratic minorant at `k→0`" is exactly `λ·g ≤ Hess S_eff`, and its top eigenvalue `λ` **is** the gap in curvature language (`H1-IS-IT-THE-GAP`, `KEY-DIRECT-DERIVATION §0`).

The fixpoint route must define `Q` **from the measure `μ` and the Hessian**, *not* from `λ=max{m₁…}`. There are exactly three ways to read "`Q≤Hess`", and they behave completely differently:

| reading of `Q ≤ Hess` | what `Q` is | sign-definite when? | verdict |
|---|---|---|---|
| **(A) pointwise / form order** | `Q⪰0`, `Hess V − Q ⪰0` ∀ fields | only if `Hess V ⪰0` everywhere | **[T]** empty on Savvidy |
| **(B) measure-averaged** | `Q⪰0`, `E_μ[Hess V] − Q ⪰0` | always (Stein) | **[T]** exists, gap-blind |
| **(C) conditioned / effective** | `Q := Hess S_eff = G″` | iff `D≤S` | **[O]** = `H1` |

The whole content of this note is that **(A) is impossible, (B) is vacuous at `k→0`, (C) is the gap.** The "fixpoint" only exists non-trivially as (B), and (B) does not floor.

---

## §1 · The fixpoint `Q = 𝒬[Hess, μ]`, and the maximal-minorant theorem

**Definition (the extremal relation).** For the coarse mode `v`, define the **self-consistent minorant**

```
       Q(μ)  :=  max { q ⪰ 0 :  q ≤ E_μ[Hess V] }         (maximal PSD form under the AVERAGED Hessian)
```
(the pointwise version — replace `E_μ` by the essential infimum `μ-essinf Hess V` — is (A); the conditioned version — replace `E_μ` by the Hessian of `−log p₀` — is (C)). The **fixpoint** is `(Q,μ)` solving `Q = Q(μ)` together with the equilibrium consistency `μ = equilib(β, Q)`.

### Theorem Q1 (maximal-PSD-minorant dichotomy). `[T]`

Let `H` be a symmetric form. A PSD form `Q` with `Q ≤ H` (form order, i.e. `H−Q⪰0`) exists **iff** `H ⪰ 0`; and then `Q=H` is the maximal one.

*Proof.* `Q ⪰ 0` and `H−Q ⪰ 0` give `H = Q + (H−Q) ⪰ 0`. Conversely `H⪰0` admits `Q=H`. ∎

**Corollary (pointwise face is empty on Savvidy).** At any field with `λ_min(Hess V)=−g|B|<0` (the Savvidy shell — measured `21.7%` negative directions at SU2 β=2.5, W6; `7.1%` after the exact gauge quotient, W4), **no** PSD form `Q` satisfies `Q ≤ Hess V`. *Witness Part 1:* the indefinite `H=[[2,0.6],[0.6,−1]]` (`λ_min=−1.12`) admits no PSD `Q≤H` (exhaustive grid search finds none — as forced). So a **pointwise** quadratic minorant is sign-definite only where the Hessian is *already* PSD, i.e. where the gap is assumed. This is the circularity trap `PROFILE` flagged: `[NG]` for the pointwise class.

### Theorem Q2 (expected-Hessian minorant = force covariance). `[T]`

For every measure `μ` on the compact fiber (periodic, no boundary), the **full** expected Hessian is PSD, by the exact Stein / Schwinger–Dyson identity:
```
   E_μ[Hess V]  =  E_μ[∇V ⊗ ∇V]  ⪰  0 .          (SD/Stein — exact, all β, no convexity)
```
Hence the maximal PSD form `≤ E_μ[Hess V]` is `E_μ[Hess V]` itself, and on the coarse mode
```
   Q(μ)_v  =  E_μ[(∂_vV)²]  =  E_μ[f_v²]  ⪰  0 ,     f_v := ∂_vV .
```
*Proof.* IBP on the periodic fiber: `0=∫∂_v(f_v e^{-V})=∫(H_{vv}−f_v²)e^{-V}` ⇒ `E[H_{vv}]=E[f_v²]`; PSD because `f_v²≥0`. ∎ *(Witness Part 2(a): `E[H_xx]=E[f_x²]=1` to machine precision for all couplings of the bounded Gaussian; H1-HESSIAN-DIRECT §1.1 banked as `[T]`.)*

**This is the only non-circular construction.** `Q(μ)=E[f_v²]` is sign-definite **without presupposing the gap**: it is the force second moment, defined from `μ` and `V` alone. **But** it is **gap-blind** (Corollary below). Note the *budget reading* of Q2: `E[(2adF_B)_vv]=E[f²]−E[(−D²)_{vv}]−E[Φ″] ≥ −(E[−D²]+E[Φ″])` — the Savvidy term provably never eats the budget `[T]`. This is the honest sense in which "the Savvidy indefiniteness is compensated by the measure" — **in expectation** (§2).

### Theorem Q3 (the `Q=0` fixed point, and its indifference to the gap). `[T]`

For the expected-Hessian fixpoint `Q(μ)_v = E_μ[f_v²]`:
```
   Q(μ)_v = 0   ⟺   f_v = ∂_vV = 0   μ-a.s.   ⟺   μ is supported on the critical set {∂_vV=0} .
```
Such `μ` include: (i) the **free massless field** — on the null mode `n` of the free inverse propagator `A` (`An=0`), `∂_nV ≡ 0` **identically**, so `Q_free(k→0)=0`; (ii) any **gapless interacting** theory with a flat coarse mode (a CFT with `m₁=0`). Hence `0 ∈ {Q : Q=Q(μ)}`, i.e. **the gapless point is a fixed point of the averaged relation**, and the relation *cannot exclude it* — it outputs `0` there. (For the conditioned relation (C) the same holds: `G″_free = S−D = 0−0 = 0`, `INVENT-VARIANCE-INEQ §4`.)

Theorem Q3 is the averaged face of the corpus's **gapless-solution obstruction** (Schema C, `MS-GAP-EQUATION §2.3`; Theorem SA, `PROFILE-FACET-ANCHOR §2.6`): `F(0)=0` for the free field. **`Q=0` is admissible. There is no source at `Q=0`** — see §3 for the dimensional reason.

---

## §2 · Can the `k→0` Savvidy indefiniteness (`λ_min=−g|B|`) be compensated self-consistently?

The task's central question. The answer splits exactly by the three faces:

**(2a) Pointwise — NO. `[T]`** By Theorem Q1, where `λ_min=−g|B|<0` no PSD minorant exists. The Savvidy negative mode is pointwise physical (not gauge: it survives the transverse gauge quotient, `7.1%` neg, W4; every negative mode is transverse-magnetic with `δF≠0`). So pointwise compensation is impossible — you cannot subtract an indefinite direction and remain `⪰0`. *(Witness Part 1.)*

**(2b) In expectation — YES, provably, but gap-blind. `[T]`** The Stein identity (Q2) shows the pointwise Savvidy negativity is *exactly* compensated by the force second moment in expectation: `E[f_v²] ≥ 0` for every measure. This is a genuine true theorem (H1-HESSIAN-DIRECT's "expectation-positive though pointwise indefinite", `[T]`). **But** it holds for the **free massless field too** (`E[H]=E[f²]=0`, `Var(a)=∞`): it is met by the very theory `Q` must *fail*. It is therefore **vacuous in the separating sense** — it does not fire. The compensation is real but lives at the *unconditioned* level (L3 of H1-HESSIAN-DIRECT), which is **on the free-field side of the acid test**.

**(2c) Self-consistently at the level of the effective mass — NO, it is ADVERSARIAL. `[V-here]`** Make the compensation self-consistent: let the coarse curvature `Q` be the curvature of the effective potential `S_eff`, and let the effective flat mass be induced by the interaction. In the exact marginal toy (`KEY-DIRECT-DERIVATION §D`, re-run Witness Part 3),
```
   G″(0)  =  m²_bare − g² Var(y²|0) ,        Var(y²|0)>0 ,
```
the Savvidy-type term enters with a **minus**: the nonabelian interaction **reduces** `G″` (drives it toward/through zero), it does not raise it. Closing the loop (`k:=Q`, induced mass) gives the fixpoint `Q = b − g²Var(y²|0;k=Q)`; with `b=0` (pure YM, no bare mass) the converged fixed point is **`Q*<0`** (`−0.39` at `g=0.6`, `−1.09` at `g=1.0`, Witness Part 3b) — the Savvidy sign. So at the conditioned/self-consistent level the Savvidy indefiniteness is **not compensated but reinforced**: it is the *adversary* of a positive `Q` (matches `KEY-DIRECT-DERIVATION §2` "nonabelian `K` is adversarial"; `PROOF-COARSE-POINCARÉ §a` `PᵀHP` min `−0.3`). The conditioned compensation is the SCBL inequality `D<S` (`G″=S−D`), and `D>0` **is** the indefinite fiber: proving `D<S` uniformly is `H1` (`INVENT-VARIANCE-INEQ §5`).

**Net [T]/[O].** "Savvidy compensated self-consistently" is true **only** at the measure-expectation level, where it is provably **gap-blind**; at the conditioned level that would floor the spectrum it is **adversarial**, and the compensation inequality is `H1`.

---

## §3 · The gapless-fixed-point test: is there a source at `Q=0`?

### 3.1 `Q=0` is admissible — yes, and it is the free field `[T]`

By Theorem Q3, `Q(μ)_v=0 ⟺ ∂_vV≡0 μ-a.s.`. The free massless field lies on this set exactly (null mode `n` of `A`). Hence the fixpoint relation `R(Q,μ): Q=Q(μ)` is **satisfied by `Q=0`** (the free/gapless point), both in the averaged face (Q2) and in the conditioned face (SCBL: `S=D=0`). Witness Part 3: free massless `m²=0,g=0` gives `G″=0.000000`. **`Q=0` admissible — the relation does not exclude the gapless point.**

### 3.2 Is there a source at `Q=0`? — No. The no-source theorem `[T-structure]` + `[V-here]`

A "source" would be a `Q`-independent, additive, **positive** term `s>0` in the relation, so that `R(0)=s≠0` and `Q=0` is *not* a solution. There are only three candidate sources, and each fails:

**(S1) Bare mass.** A `+m²_bare` term is a source — it is exactly what makes the toy's `Q*>0` for `b>0` (Witness Part 3b). **Pure YM has no bare mass** (`m²_bare=0` by gauge invariance; the lattice has no mass parameter). Remove it and the self-consistent curvature is `Q*≤0` (Witness Part 3b, `b=0` rows). So the only demonstrated source is one pure YM lacks.

**(S2) The Haar/Jacobian term `Φ″⪰N/12`.** This is the one pointwise-positive, measure-independent, β-independent addend (`[T]`, Lemma B). It is a real candidate source — *but* (i) it is **gauge-dependent** (it is the curvature of the flat-chart Jacobian `Φ=−log J`, and `X=log U` is chart-dependent): the *gauge-invariant* `Q` cannot equal `Φ″` alone; it enters only through the gauge-fixed marginal `W_geo=E[Φ″]−Var(Φ′)≥0` (Prékopa `[T]`). (ii) In that gauge-fixed marginal it is **absorbed by the Stein budget**: `E[2adF_B] ≥ −(E[−D²]+E[Φ″])`, so the negative Savvidy term can consume all of `E[Φ″]`; the *total* `G″=Δ_dyn+W_geo+X` is what is physical, and at the free field `Δ_dyn+X` cancels `W_geo` to give `G″=0` (SCBL acid test). So `Φ″` is a *positive addend*, not a *net source*: it does not prevent `Q=0`. (`PROOF-JACOBIAN-SAVES`: `W_geo≥0` is provable all k, but `G″>0` needs `W>Var−E[H]` = the conditioned variance bound = `H1`.)

**(S3) The trace anomaly.** YM is classically scale-free; the **only** explicit scale-breaker is the anomaly `T^μ_μ=(β/2g)tr F²`, `β<0`, `b₀=11N/48π²>0`. Could it source the scale? **No**: in the dilation low-energy theorem the anomaly relation is `⟨0|O|G⟩ = f_G m_G² g/(2β)` with `O=¼trF²` — it is **homogeneous in the residue `f_G` and the mass `m_G`** (`KEY-DIRECT-DERIVATION §4 Attack IV`, sums rules are moments ⇒ *upper* bounds; the `0⁻⁺` partner `trF\tilde F=∂_μK^μ` is a decoupled total derivative, so there is **no signed superconvergence**). Witness Part 4: the relation fixes `m_G` in terms of `f_G` but carries **no `Q`-independent constant** — it cannot *floor* either. So the anomaly does not provide a source at `Q=0`; it provides a *homogeneous* relation.

**Theorem Q-S (no source ⇒ nonzero fixpoint = the gap). `[T-structure]`** In pure YM (no bare mass, only the dimensionless coupling and the regulator `a`), any fixpoint `Q=𝒬[Hess,μ]` of mass dimension 2 that (i) is built from `V` and `μ`, (ii) is not given a scale by hand, is **scale-covariant**: under `x→sx` (`a→a/s`), `Hess V→s² Hess V`, so `Q→s²Q` and `Q=0` is invariant. A **nonzero** solution `Q>0` therefore **spontaneously breaks scale invariance** — it *is* dimensional transmutation, *is* `Λ_YM`, *is* the mass gap. There is no `Q`-independent constant to lift `Q=0` off the relation (S1 absent, S2 absorbed, S3 homogeneous). **Corollary:** the quadratic fixpoint **cannot supply** the scale; it can only **have** it if the gap exists. `Q=0` admissible, **no source** — the exclusion of `Q=0` is exactly `H1`.

*This is the sharp, new result of this round (`[V-here]`): the "6th species as fixpoint" is pinned not to a technical gap but to a **structural no-source property** — a sign-definite quadratic minorant defined self-consistently over the measure has the gapless point as a fixed point because pure YM is classically scale-free and the only scale-breaker is homogeneous.*

---

## §4 · Uniform attainment at `β≈2.5` (and why it is `I2`)

The fixpoint must be **attained uniformly in `(L,a)`** across the crossover `β≈2.47`. Three levels, matching §0:

**(4a) Measure-averaged (B): attained trivially, not a floor. `[T]`** `Q(μ)=E[f_v²]` exists for every `(L,a)` by construction (Q2). But it is **gap-blind**: it is `≥0` for every measure and `=0` on the free field, and it does not imply the conditioned `G″≥λ` (a bimodal `p₀` with a dip at `0` has `G″<0` while `E[H]=E[f²]≥0`; H1-HESSIAN-DIRECT §3). So uniform attainment *of (B)* is vacuous.

**(4b) Conditioned (C): measured positive, uniformity open. `[T-measured]/[O-uniform]`.** The conditioned fixpoint `Q=G″(k→0)` is **measured positive** at the crossover: `G″=+0.51±0.08` (SU3 β=2.5, exact marginal on the softest blocked mode, W7 Part B), `+0.27/+0.32` (SU2 β=2.5 / SU3 β=5.7 finite-`k` plane waves, W7-RERUN), with `1/Var` volume-independent (0.272 at L=4 *and* L=6). But the *uniform* statement `G″≥λ` ⟺ `Var_μ(f|0) ≤ E_μ[H|0]+W−λ` ⟺ `D<S` ⟺ `Var_ν(φ)≤(1/λ)∫|∇φ|²dν` is **not proved** and **is `H1`** (`PROOF-COARSE-POINCARÉ §exact`, `H1-HESSIAN-DIRECT §7`, `INVENT-VARIANCE-INEQ §5`). The `a→0` uniformity — turning "positive at each `(L,a)`" into "positive in the continuum limit, uniformly" — is the zero-finding/continuum `I2` wall (`META-CROSSDOMAIN §1.4`, `LOC-DETERMINANT §5`).

**(4c) Constructive attainment (iterating to the fixpoint): fails uniformly. `[T][NG]`** The physical Riccati/FRG map that would *iterate* to `Q` has `ρ(DT)→1` at `k→0` (**neutral zero mode**) and `ρ(DT)=∞` on the Savvidy shell; the `k→0` mode is **neutral**, so no uniform contraction, and the only contractively-certified fixed point near `κ=0` is the **gapless** one; the nonzero branch is a **bifurcation** whose insertion is the order parameter (`BANACH-FIXED-POINT §2–3,§5`, `[T]`). So *constructive* uniform attainment is `[NG]`.

**Uniform-attainment verdict.** (B) is attained but vacuous; (C) is measured positive and its uniform attainment **is** `H1 = I2`; (C-constructive) is `[NG]` (neutral mode). The fixpoint is **not** uniformly attained across `β≈2.5` by any route that does not already contain the gap.

---

## §5 · Acid test, the three iron laws, and the profile constraints P1–P4

### 5.1 Acid test `[T]` — passed only in the vacuous sense

A gap-certifying `Q` must **fail** on the free massless field (give `≤0`), and **fail** on the gapped-atomless `:φ²:` if targeting the atom. The fixpoint `Q`:
- free massless `k→0`: `Q=0` (Witness Part 3). **Meets `≤0`** — acid-clean, non-spurious (not `>0` on the free field).
- massive free composite `:φ²:`: `G″=1/Var(a)>0` (massive ⇒ finite `Var`), but `f₁²=0`. So the **curvature** fixpoint is `>0` on a gapped atom-free theory. It is therefore a **support** object (`G″`), not an **atom** minorant: `Q>0 ⇏ f₁²>0` (the squeeze `f₁²/χ≤m₁²≤M₀/χ` only *upper*-bounds `m₁²` by `G″=1/χ`; `H1-IS-IT-THE-GAP`). So the fixpoint passes the acid test only as a **support/curvature** quantity; as an **atom** certifier it **falsely fires on `:φ²:`** (`[T]`, Theorem 1 of `SPECIES-EXISTENCE`). **It is not an atom minorant.**

### 5.2 Three iron laws

| law | the fixpoint `Q` | label |
|---|---|---|
| **(IL-i) non-positivity** (must not be satisfied by free field *and* must fail there) | (B) is sat by free field (`0`) ⇒ **vacuous**; (A) impossible; (C) is the variance bound = `H1` | (B) ✗ / (C) `[O]` |
| **(IL-ii) uses nonabelian interaction** | yes where it matters: Savvidy `2adF_B` is nonabelian-only (`[T]`), but it enters **adversarially** (drives `Q` down) | ✓ but adversarial |
| **(IL-iii) no small parameter** | (B) exact all β; (C) unproven at `α_s~O(1)`; the constructive iteration needs contraction (fails) | (B) ✓ / (C) `[O]` |

### 5.3 Profile constraints P1–P4 (`PROFILE §2`)

- **(P1) a functional on the field measure `μ`** — ✓ the fixpoint is `Q=𝒬[Hess,μ]`, a measure functional (not a spectral statistic of `ν`).
- **(P2) `𝔉=𝔉_*+J`, all content in the usc floor spike `J`** — the fixpoint gives `Q` in `(B)` as an **lsc/continuous** measure functional (`E[f²]`), i.e. the `𝔉_*` part (`≤0`-vacuous on the gapless class); it does **not** supply the floor spike `J`. In `(C)` the spike is the positivity of `G″`, which is `H1`. So P2 is **not** met by the fixpoint — it reproduces the `𝔉_*` half only.
- **(P3) `J=inf_λ φ_λ`, computable upper bounds, atom-blind, threshold-firing on the two-particle cut** — the fixpoint is atom-blind ✓ but does **not** fire on the `:φ²:` cut (it is `>0` there); it has no `inf_λ` structure. ✗.
- **(P4) fixed point of an intrinsic relation `R` excluding the gapless solution by a `k→0` sign-definite mechanism** — the fixpoint **is** intrinsic (controller-free) ✓ and **is** `k→0` ✓, but it **does not exclude** the gapless solution: `Q=0 ∈{Q:R}` (Thm Q3), **no source** (Thm Q-S). ✗ — and this is the exact slot.

**Scorecard: P1 ✓, P4-partial (intrinsic & k→0 but `Q=0` admissible), P2/P3 ✗.** The fixpoint is the intrinsic `k→0` object but it is the **`𝔉_*`** (vacuous-positivity) half, not the floor spike `J`.

---

## §6 · Verdict — constructed or exact failing step

**Not constructed.** The quadratic-minorant fixpoint exists in exactly two forms and neither fills the slot:

1. **(A) pointwise** — **empty** on the Savvidy shell (`[T]`, Thm Q1: `Q⪰0`&`Q≤Hess`⇒`Hess⪰0`). A pointwise sign-definite minorant presupposes the gap.
2. **(B) measure-averaged** — **exists** (`Q=E[f²]`, Stein `[T]`), non-circular, sign-definite — **but gap-blind and vacuous at `k→0`**: `Q_free=0`, it never fires, and it does not imply the conditioned curvature. It is the `𝔉_*` half.
3. **(C) conditioned / effective-potential** — the object where the gap lives, `Q=G″(k→0)=S−D`; **`Q=0` admissible**, **no source at `Q=0`**, uniform `Q≥λ` is `H1=I2` (`[O]`).

**The exact failing step (`[O]`, sharpest form, `[V-here]`):**
> The fixpoint relation is **source-free at `Q=0`**: `Q=0 ⟺ ∂_vV≡0 μ-a.s.` — which the free/gapless field satisfies — and **pure YM is classically scale-free**, so the only explicit scale-breaker (the trace anomaly) enters **homogeneously and cancels** (no `Q`-independent constant). The Haar Jacobian's `Φ″⪰N/12` is a real positive addend but is gauge-dependent and absorbed by the Stein budget (free field gives `G″=0`); the only demonstrated source is a **bare mass**, which pure YM lacks. Hence `Q=0` is a fixed point and **a nonzero fixed point `= spontaneous scale breaking = dimensional transmutation = the mass gap.** The Savvidy indefiniteness is compensated **only** at the gap-blind expectation level (Stein); at the conditioned level it is **adversarial** (`Q*≤0`).

**This is `H1`, in the fixpoint coordinate.** It is the same wall as `PROFILE-FACET-ANCHOR`'s "exclude the gapless solution by a `k→0` sign-definite mechanism," `MS-GAP-EQUATION`'s `F(0)=0`, and `SPECIES-EXISTENCE`'s `F3∧F4=∅` — now stated as a **no-source / scale-homogeneity** property of self-consistent quadratic minorants: **the relation cannot supply the scale it would need to exclude its own zero.**

**Flip condition (`[T]`-criterion).** The slot is filled iff someone exhibits a sign-definite quadratic minorant fixpoint `Q=𝒬[Hess,μ]` with a **`Q`-independent source** `s>0` making `Q=0` inadmissible — i.e. a positive, gauge-invariant, interaction-generated constant of mass dimension 2 in pure YM. By Thm Q-S that constant **is** the mass gap; by `KEY-DIRECT-DERIVATION` no such object is known (the Savvidy scale `√(2gB)→0` at the physical vacuum; the anomaly is homogeneous; every gapped reference is on the wrong side of the min–max order). So the flip condition is unsatisfiable without `H1`.

**Honest boundary.** `∀class ≠ ∀possible`: this is `[NG]` for the *quadratic-minorant-fixpoint class* (three faces A/B/C, the natural extremal definitions `[T]`/`[O]`). A non-quadratic, non-extremal, non-source-free object could in principle escape — but by `P2`/`PROOF-FACET-ANCHOR §3.3` any such object that is floor-sensitive and self-anchoring is the Fredholm-zero / transport object, which is `[T-shape]`/`[O]` and stalls at `a→0`. The quadratic fixpoint is **not** it. **Clay 0%; no QED.**

---

## §7 · References

**Witness (this round):** `research/unconditional/w7/slot_quadratic_fixpoint.py` (`py -3.12`; Part 1 pointwise trap exact linear algebra; Part 2 bounded-Gaussian Stein `E[H]=E[f²]` to machine precision + truncation = conditioning; Part 3 exact `G″=m²−g²Var` marginal; Part 3b self-consistent closure `Q=b−g²Var(k=Q)`, `b=0`⇒`Q*<0`; Part 4 anomaly homogeneous). No lattice MC; one small job.

**Internal (absolute paths; banked, labels carried):** `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall = Savvidy concavity / key = `G″(0)`/H1); `research/unconditional/PROFILE.md` (§3 the one open slot; P1–P4); `PROFILE-FACET-ANCHOR.md` (§1.2 Thm C; §1.3 Thm B; §2.6 Thm SA; §3.2 residual DOF; §3.3 flip); `KEY-DIRECT-DERIVATION.md` (§1 Savvidy attack; §2 `G″=S−D` cancellation, `G″=m²−g²Var` toy; §4 no signed channel; §7 exact step); `H1-HESSIAN-DIRECT.md` (§1.1 Stein `E[H]=E[f²]`; §3 the gap is L3→L4 conditioning; §4 four levels; §6–7 exact stuck step); `INVENT-VARIANCE-INEQ.md` (SCBL `G″=S−D`; acid test free `S=D=0`; `D<S`=`H1`); `PROOF-JACOBIAN-SAVES.md` (`W_geo≥0` Prékopa; `G″=Δ_dyn+W_geo+X`; `W>Var−E[H]`=`H1`); `PROOF-COARSE-POINCARÉ.md` (`Hess S_eff⪰λ`⟺Poincaré⟺`H1`; `PᵀHP` min `−0.3`); `PROOF-SAVVIDY-EXACT-SIGN.md` (`−D²+2adF_B`, `λ_min=−g|B|`, count/depth, nonabelian-only, `Tr(2adF)=0`); `H1-IS-IT-THE-GAP.md` (squeeze; `G″=1/χ`; `H1≠H*`); `MS-GAP-EQUATION.md` (§2.3 Schema C `κ=F(κ)`, `F(0)=0`); `BANACH-FIXED-POINT.md` (neutral `k→0`; Savvidy expansive; bifurcation); `SPECIES-EXISTENCE.md` (Thm 1 discriminator≠certifier; Thm 2 `F3∧F4=∅`; Thm 3 equivalence); `META-GAPFORCING-IDENTITY.md` (READ≠FORCE); `W7-Gpp0`/`W7-RERUN` (`G″` +0.27/+0.32/+0.51; `E[H]≈Var≈2.45`; `1/Var≈0.27`); `MARGINAL` (`PᵀHP −0.3`); `W4`/`W6-BE-at-2.5` (Savvidy bottom `−1.12β`, 7.1%/21.7% neg).

**External `[R]`:** Savvidy, *Phys. Lett. B* **71** (1977) 133; Nielsen–Olesen, *Nucl. Phys. B* **144** (1978) 376; Leutwyler, *Nucl. Phys. B* **179** (1981) 129; Brascamp–Lieb, *Adv. Math.* **22** (1976) 366; Prékopa, *Acta Sci. Math.* **32** (1971) 301; Stein's lemma / Schwinger–Dyson IBP; Efron–Stein (1981); Courant–Fischer (Reed–Simon IV); Shen–Zhu–Zhu, *CMP* **400** (2023) 805 — strong-coupling gap; Cornwall, *Phys. Rev. D* **26** (1982) 1453 (effective mass `m_eff²=m₀²+Π(0)`); Coleman–Weinberg (1973) dimensional transmutation.

---

*Status: **Quadratic-minorant-fixpoint constructive attempt — NOT constructed.** `[T]` pointwise minorant empty on Savvidy (Thm Q1); `[T]` expected-Hessian minorant exists = `E[f²]` but gap-blind, `Q_free=0` (Thm Q2/Q3); `[O]` conditioned effective-potential fixpoint `G″` has `Q=0` admissible and **no source** (Thm Q-S: pure YM scale-free, bare mass absent, Jacobian absorbed, anomaly homogeneous) ⇒ nonzero fixpoint = dimensional transmutation = the gap; uniform attainment = `H1=I2`; constructive iteration `[NG]` (neutral `k→0`). The Savvidy indefiniteness is compensated only at the gap-blind expectation level and is adversarial at the conditioned level. **`Q=0` admissible; no source; = `H1`.** Clay 0%; no QED. `∀class ≠ ∀possible`.*

# PROFILE-FACET-SEMICONT — the SEMICONTINUITY facet of the missing species `𝔉`: what upper-semicontinuity structure a gap-certifier must have, and why lsc/continuous functionals are provably blind

> 2026-11-24 · **Facet note** of the "profile" (wanted-species specification, `CLAUDE.md §1/§5`). Not a slogan; a **precise structural spec** for `𝔉`. Builds on `SPECIES-EXISTENCE.md` (`ℭ`, F1–F4, `F3∧F4=∅`), `META-GAPFORCING-IDENTITY.md` (READ ≠ FORCE), `SUPPORT-NOT-ATOM.md` (Clay = support floor), `H1-IS-IT-THE-GAP.md` (support target (G) vs atom target (H\*)), `PROOF-NEVANLINNA-POLE.md` (atom usc, germ fiber convex), `TOOL-COMPLETELY-MONOTONE.md` (finite-`t` rate = upper bound; only the `t→∞` limit is clean), `LEVER-nonconvex-edge.md §0B` (linear-data lemma).
>
> Labels (charter §2): **[T]** true theorem (proved here or an elementary consequence of a banked `[T]`) · **[C]** conditional theorem · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[J]** judgment · **[V-here]** new this round · **[R]** inherited/standard. **No QED on anything unproven. Clay stays 0%. `∀class ≠ ∀possible`.**
>
> **Scope of this facet.** The spec has many facets (lagrangian/field-space, transport/Fredholm, carrier, uniformity). THIS note isolates **one**: the *semicontinuity* structure. It answers: *what type of function of the measure must `𝔉` be, so that it can see a **usc** support-floor without being (a) a positivity/lsc bound (vacuous) or (b) an exact identity (READ ≠ FORCE)?* The answer is a hard constraint plus a residual DOF — derived below.

---

## VERDICT (up front)

> ### The species must be **upper semicontinuous and NOT lower semicontinuous** on the space of Källén–Lehmann measures; equivalently it must carry a **strict upward jump ("floor spike") whose lower-semicontinuous envelope is `≤ 0` everywhere**. Precisely, every admissible `𝔉` decomposes `𝔉 = 𝔉_* + J` with `𝔉_*` lsc and `≤0` on the dense gapless class (the **vacuous positivity part**) and `J = 𝔉−𝔉_* ≥ 0` a usc functional supported on the (closed, nowhere-dense) gapped set (the **entire species**). No lsc/continuous `𝔉`, and no functional that factors through a data map whose fibers contain a vaguely-dense gapless set, can satisfy F1 ∧ F3.
>
> **Three exactly-scoped statements, one per question:**
> 1. **[T] THE OBSTRUCTION.** On `(𝔐, vague)` the floor `m₁ = inf supp ν` and the edge atom `f₁² = ν({m₁})` are **upper semicontinuous**; the gapped sets `G_δ = {m₁ ≥ δ}` are **closed with empty interior**; the gapless measures `{m₁ = 0}` are **vaguely (and weakly, and TV-) dense**. Hence any **lsc** `𝔉` with `𝔉 ≤ 0` on the gapless class is `≤ 0` on *all* of `𝔐` (§1, Theorem 1). Continuous `𝔉` is a special case (§1, Cor 1). This *is* the obstruction, stated topologically; it subsumes the four-mode enumeration of `SPECIES-EXISTENCE` Theorem 2 and sharpens it (`F4` ⟹ continuity/lsc ⟹ blind, §1.5).
> 2. **[T]/[NG] WHICH lsc CLASSES SURVIVE.** The class "lsc **and** support-floor-sensitive **and** acid-clean" is **EMPTY** (Theorem 1 with the acid clause). The class that survives is exactly the **usc-with-a-jump** class; its canonical members are `m₁`, `f₁²`, the closed-set indicators `1_{G_δ}`, and the Fredholm minimal-zero `z₁` (§2, catalog). Each survives the semicontinuity test **structurally** but is either *the answer* or *reduces to the transport/zero-finding problem* — i.e. fails **F4**, not F3 (§2.3).
> 3. **[O] THE REQUIREMENT + RESIDUAL FREEDOM.** `𝔉` must (R1) be usc-not-lsc, (R2) have lsc envelope `≤0` ("vacuous-positivity separation"), (R3) have its spike supported on `{m₁ ≥ δ_𝔉}` with **spike location = `inf supp ν`**, (R4) attain the limit **uniformly in `(L,a)` across the crossover**. The **residual DOF** is the *parametrization of the spike*, not its value: choice of (i) data fiber/carrier operator `O`, (ii) threshold `δ`, (iii) the **sharp family of computable upper bounds** `{φ_λ}` whose *infimum* realizes the spike (Tauberian `t`, Christoffel `k`, Fredholm `z`, continuum `a`). All members evaluate to `m₁` in the limit — the content is the *location* of the spike, which **is** `inf supp ν` = the gap (§3). The facet therefore **localizes** the missing step to "uniform attainment of the semicontinuous limit" and **forbids** any continuous/lsc substitute.

---

## §1 · The obstruction, stated rigorously

### §1.1 The space and the topology

**Space.** `𝔐` := finite positive Borel measures `ν` on `[0,∞)` (`s = m²`), the Källén–Lehmann spectral measures `ν_μ(B) = ⟨Ω|O† P_H(B) O|Ω⟩`, `O = ¼ tr F²`. `supp ν ⊂ [m₁²,∞)`, `m₁ = inf supp ν`, edge atom `f₁² = ν({m₁²})`, `M₀ = ∫dν = ‖O‖² < ∞`, `χ = ∫s⁻¹dν`. Physical theories need `ν` **renormalized (point-split)** to make `M₀ < ∞` (the UV tail is a fixed upper end; it does not affect the IR floor — `TOOL-COMPLETELY-MONOTONE §2`). The sub-class of interest is the **OS/RP-admissible** `ν`; it contains every generalized free field (Jaffe: any `ν ≥ 0` gives an OS-admissible GFF two-point function, `[R]`), in particular the free massless field (`ν` a.c., `m₁ = 0`, `f₁² = 0`) and the corpus witnesses `wδ_ε+(1−w)δ_S` (`REDUCTION-AUDIT (b)`).

**Topology.** 
- **Vague** = weak-* against `C_c((0,∞))` (resp. `C_0` for finite `ν`): `ν_n → ν` iff `∫φ dν_n → ∫φ dν` for all bounded continuous `φ`. 
- **Weak** = vague + tightness/mass-preservation (Prokhorov; on `[0,∞)` Polish, `ν_n ⇒ ν` iff vague and `ν_n([0,∞)) → ν([0,∞))`). 
- **TV** (strong) = total variation.

The physically relevant compactness (the corpus's `GL3`, Prokhorov) is **vague/weak**. I prove below that the obstruction is **topology-robust** (holds in all three), so the choice does not rescue anything; but the *canonical* statement is vague/weak. Throughout, "nbd" means vague-nbd.

### §1.2 Two structural facts (both `[T]`)

Recall `supp ν = {s : ν(N) > 0 for every nbd N of s}`, `m₁ = inf supp ν` (with `m₁ = +∞` if `ν = 0`).

**Proposition 1 (floor and edge atom are usc). [T]** If `ν_n → ν` vaguely (hence weakly, hence in TV), then
```
limsup_n m₁(ν_n) ≤ m₁(ν),        and        limsup_n f₁²(ν_n) ≤ f₁²(ν).
```
*Proof.* *(Floor.)* Suppose `limsup m₁(ν_n) > m₁(ν) =: m`. Fix `m < η < limsup m₁(ν_n)`. Along a subsequence `m₁(ν_n) ≥ η`, i.e. `ν_n((0,η)) = 0`. Two cases. If `ν((m,η)) > 0`: since `(m,η)` is open, Portmanteau's lower bound gives `liminf ν_n((m,η)) ≥ ν((m,η)) > 0`, contradicting `ν_n((m,η)) = 0`. If `ν((m,η)) = 0`: then `m ∈ supp ν` forces `ν({m}) > 0`; take `φ ∈ C_c`, `φ ≡ 1` on `[0,m]`, `φ` decreasing to `0` on `[m,η]`, `φ ≡ 0` on `[η,∞)`; then `∫φ dν_n = 0` (support `≥ η`, `φ=0` there) for all `n`, so `∫φ dν = lim ∫φ dν_n = 0`, i.e. `ν([0,m]) = 0`, contradicting `m ∈ supp ν`. Both impossible. ∎ *(Edge atom.)* Same `φ`-test with `φ ≡ 1` on a nbd of the (moving) floor and the fact that the total mass of `ν_n` on `[0,η]` is `≤ ν([0,η']) + o(1)` for `η' > η`; the atom `ν_n({m₁(ν_n)}) ≤ ν_n([0,η])` climbs to at most `ν([0,η])`, and letting `η ↓ m₁(ν)` yields `limsup f₁²(ν_n) ≤ f₁²(ν)`. ∎

**Proposition 2 (gapless measures are dense; gapped sets are closed-nowhere-dense). [T]** For every `δ > 0`:
```
(2a)  {ν : m₁(ν) = 0}  is vaguely-, weakly-, and TV-dense in 𝔐;
(2b)  G_δ := {ν : m₁(ν) ≥ δ}  is vaguely-, weakly-, and TV-closed, with EMPTY interior.
```
*Proof.* (2a) For any `ν ∈ 𝔐` set `ν_ε := (1−ε)ν + ε δ_0`. Then `ν_ε({0}) = ε > 0` so `m₁(ν_ε) = 0`, and `ν_ε → ν` vaguely (`∫φ dν_ε = (1−ε)∫φ dν + εφ(0) → ∫φ dν`), weakly (`‖ν_ε‖ = ‖ν‖`), and in TV (`‖ν_ε − ν‖_TV = ε‖δ_0 − ν‖_TV → 0`). (2b) Closed = Prop 1 (usc ⟹ superlevel sets closed). Empty interior: for `ν ∈ G_δ`, the `ν_ε` of (2a) lie in `G_δ^c` and `→ ν`; hence every `ν ∈ G_δ` is a limit of `G_δ^c`, so `int G_δ = ∅`. ∎

> **This is the whole story in two lines.** The gap is a **closed** property (stable under limits — the good news), but its certificate set is **nowhere dense** and the *complement is dense* (the obstruction). Any functional continuous/lsc enough to be pinned on the dense gapless set is pinned `≤0` everywhere.

### §1.3 The obstruction theorem (the core)

Let `𝔉 : 𝔐 → ℝ` and recall F1/F3 from `SPECIES-EXISTENCE §1`: **F1 (acid)** `𝔉(ν) ≤ 0` for all gapless `ν` (contains the free massless field); **F3 (detect/ certify)** `∃ ν ∈ G_δ : 𝔉(ν) ≥ c > 0` for some `δ, c > 0` (uniform over a family).

**Theorem 1 (lsc no-go). [T]** If `𝔉` is lower semicontinuous on `(𝔐, τ)`, `τ ∈ {vague, weak, TV}`, then
```
F1  ⟹  𝔉 ≤ 0  on ALL of 𝔐,   hence  F3 FAILS.
```
*Proof.* Fix `ν ∈ 𝔐`. By Prop 2a choose gapless `ν_n → ν` (`τ`). By lsc, `𝔉(ν) ≤ liminf_n 𝔉(ν_n) ≤ liminf_n 0 = 0`. As `ν` was arbitrary, `𝔉 ≤ 0` on `𝔐`. ∎

**Corollary 1 (continuous no-go). [T]** Continuous `𝔉` is lsc; apply Thm 1. In particular every `𝔉 = Φ ∘ D` with `D` continuous in `(𝔐,τ)` and `Φ` arbitrary *with `Φ∘D` continuous* fails F3 under F1. This includes all of the corpus's **C0** modes: finite moments `∫φ_i dν` (continuous), the finite-`t` CM rate `r(t) = −(1/t)log(C(t)/M₀)` (continuous, `TOOL-COMPLETELY-MONOTONE §2`), deformation invariants (locally constant), and the continuous quantitative extremals.

**Theorem 1′ (fiber no-go; does not assume continuity). [T]** Let `𝔉 = Φ ∘ D` for a data map `D : 𝔐 → Y`. If for some `ν ∈ G_δ` the fiber `D^{-1}(D(ν))` contains a gapless `ν'` (i.e. `m₁(ν') = 0` and `D(ν') = D(ν)`), then F1 ⟹ `𝔉(ν) = 𝔉(ν') ≤ 0`, so F3 fails at `ν`. 
*Proof.* `𝔉(ν) = Φ(D(ν)) = Φ(D(ν')) = 𝔉(ν')`. ∎

Thm 1′ is the `SPECIES-EXISTENCE` Theorem 2 mechanism in one line, but it is *stronger in scope*: it kills **discontinuous** `Φ` too, because the filler `ν'` has *the same data*. So a survivor must avoid both Thm 1 (lsc) and Thm 1′ (gapless-dense fibers). The only way: `D` must **separate** `G_δ` from the gapless class **within every fiber** — i.e. `D` must *see the edge*. That is the entry point of §3.

### §1.4 Why an lsc lower bound cannot separate a floor at `δ` from `δ' → 0`

Take two measures in one admissible family: `ν^{(δ)}` with `m₁ = δ` and `ν^{(δ')}` obtained by pushing the same total/computable data onto a thin atom at `δ' → 0` (the corpus witnesses: `ρ = f₁²/χ`-preserving `wδ_{δ'} + (1−w)δ_S`, `H1-IS-IT-THE-GAP §3(b)`; and the exact three-moment witness below). Then `m₁(ν^{(δ)}) = δ ≠ δ' = m₁(ν^{(δ')})`. A certifier must output `≥ c` on the first, `≤ 0` on the second. But:
- **(lsc route).** If `𝔉` is lsc and `≤0` on the (dense) gapless class, then `𝔉(ν^{(δ)}) ≤ liminf` over gapless approximants `= 0` (Thm 1). The value at the gapped point is *capped by the values at arbitrarily close gapless points*, which are `≤0`. So `𝔉` cannot be `≥ c` where the floor is `δ` — it reads `δ` and `δ'→0` **identically as `≤0`**. 
- **(fiber route).** If the two share computable data (`D(ν^{(δ)}) = D(ν^{(δ')})`), `Φ∘D` is literally equal on them. To tell `δ` from `0` it must evaluate `δ` — i.e. contain `m₁` (F4 fails). 

Linear data is *insufficient* for this, provably: finitely many `∫φ_i dν` cannot separate floors across the fiber (`LEVER §0B`). **Exact witness [V-here]:** the gapped `ν_g = ½δ_1 + ½δ_3` (floor 1) and the gapless `ν_w = ⅕δ_0 + ⅘δ_{5/2}` (floor 0) satisfy `(M₀,M₁,M₂)(ν_g) = (M₀,M₁,M₂)(ν_w) = (1,2,5)` **exactly** (verified in rational arithmetic; they split only at `M₃`: `14` vs `25/2`). So the first three moments admit both a gap-1 and a gap-0 measure; any `Φ(M₀,M₁,M₂)` is blind, and *lsc* functionals of even the full germ are blind by Thm 1.

### §1.5 Sharpening of `F3 ∧ F4 = ∅` `[V-here]`

`SPECIES-EXISTENCE §2` proved `F3 ∧ F4 = ∅` **mode-by-mode** (α/β/γ/δ). The semicontinuity facet gives a **single topological proof** independent of the four-mode enumeration:
```
F4 (𝔉 reads continuous spectral data D, hence 𝔉 = Φ∘D is lsc/continuous on the data topology)
   ⟹ (Thm 1 / Cor 1)  F1 ⟹ 𝔉 ≤ 0  ⟹  ¬F3.
```
So the incompatibility is not an accident of four specific data modes — it is the **semicontinuity direction**: *positivity/lsc is the wrong semicontinuity class to certify a usc floor.* This is the load-bearing identification of this facet. (Caveat, honest: a "computable" `Φ` need not be continuous if `D` is discontinuous; but Thm 1′ then applies unless `D` separates the edge. Either way the survivor must break continuity/see the edge — `[T]` for Thm 1/1′, `[J]` for the identification "F4 ⟹ continuous **for constructed/analytic `𝔉`**".)

---

## §2 · Candidate lsc-but-floor-sensitive classes: which survive

The task asks: *is there a functional that is lsc but still support-floor-sensitive?* The answer separates cleanly.

### §2.1 The lsc-and-floor-sensitive class is EMPTY (given the acid test). [T]/[NG]

By Theorem 1, an lsc `𝔉` that is `≤0` on the dense gapless class is `≤0` everywhere, hence **cannot be `≥ c > 0` on any gapped point** — it is *floor-insensitive in the certifying direction* (it may take different *negative* values, but it can never certify `m₁ ≥ δ`). With the acid/F1 clause, **`{lsc} ∩ {floor-sensitive} ∩ {acid-clean} = ∅`.** This is a clean **[NG]** for the naturally-hoped class. (If one drops F1, the class is nonempty — e.g. `−m₁` is lsc and floor-sensitive — but every member is *sign-wrong* or *vacuous*, which is the F1 failure. `[T]`)

### §2.2 The lsc-envelope decomposition — the precise characterization (new). [T] [V-here]

Let `𝔉_*` be the **lower-semicontinuous envelope** of `𝔉`: `𝔉_*(ν) := sup{ φ(ν) : φ ≤ 𝔉, φ lsc } = min( 𝔉(ν), liminf_{ν'→ν} 𝔉(ν') )`. It is the largest lsc function below `𝔉`.

**Lemma 2 (vacuous-positivity / spike decomposition). [T]** Let `𝔉` satisfy **F1** (`𝔉 ≤ 0` on the dense gapless class). Then:
```
(i)   𝔉_*(ν) ≤ 0   for ALL ν ∈ 𝔐;                       ("the lsc part is vacuous")
(ii)  for every gapped ν with 𝔉(ν) > 0:   𝔉(ν) > 𝔉_*(ν);  ("strict upward jump")
(iii) 𝔉 = 𝔉_* + J,   J := 𝔉 − 𝔉_* ≥ 0   is usc,   J(gapless) = 0,   J(gapped) = 𝔉 − 𝔉_* > 0.
```
*Proof.* (i) Fix `ν`. By F1 every gapless `ν'` has `𝔉(ν') ≤ 0`; gapless is dense (Prop 2a), so along a gapless net `ν' → ν`, `liminf_{ν'→ν} 𝔉(ν') ≤ liminf over the gapless subset ≤ 0`; and `𝔉(ν)` may be `>0`, so `𝔉_*(ν) = min(𝔉(ν), liminf) ≤ 0`. (ii) `𝔉(ν) > 0 ≥ 𝔉_*(ν)`. (iii) `𝔉_*` lsc ⟹ `−𝔉_*` usc; `J = 𝔉 − 𝔉_*` with `𝔉` arbitrary — but at any point where `𝔉 > 𝔉_*`, `J > 0`, and `J = 0` wherever `𝔉` is itself lsc (in particular on the gapless class where `𝔉_* = 𝔉 = ≤0` on the nose if `𝔉` is lsc there). *[The precise statement: `J ≥ 0`, `J ≡ 0` on the lsc-continuity set of `𝔉`, and `J > 0` exactly at the upward-jump points.]* ∎

> **Interpretation — the facet's central formula.** **`𝔉 = (𝔉_* ≤ 0, lsc, vacuous) + (J ≥ 0, usc spike, the entire species).`** The "positivity/lsc bound" of part (a) of the task is **exactly `𝔉_*`**, and by (i) it is `≤0` everywhere — it carries **zero** certification. **All** of the gap content is in the **upward jump `J`**, which is a usc functional vanishing on the dense gapless set and positive on the (closed-nowhere-dense) gapped set — i.e. `J` is itself a **pure floor detector**. This is a *reduction*: WLOG one may take `𝔉 = J` (usc, `𝔉_* ≡ 0`). The task's requirements (a) "not a positivity/lsc bound" and (b) "not an exact identity" become: **`J ≠ 0` (not the vacuous `𝔉_*`), and `J` is not *exactly computable* (not the READ object).** A genuine FORCE is a usc spike `J > 0` uniformly on `G_δ` that is **not** the identity `m₁·1_{G_δ}` nor a Fredholm-zero readout.

### §2.3 Catalog — the candidate classes, and which survive

| candidate class | semicontinuity on `(𝔐,vague)` | floor-sensitive? | acid (F1)? | computable (F4)? | **survives as certifier?** | label |
|---|---|---|---|---|---|---|
| **C0 modes** — moments, invariant, finite-`t` rate `r(t)`, quantitative extremals | continuous (or lsc) | **no** (Thm 1/2) | achievable | yes | **NO** — blind | [T]/[NG] |
| **positivity/Schur/CM bounds** (`M_k`, `f₁²/χ`… harmonic mean, Christoffel) | lsc/continuous | no (only *upper* bounds on `m₁`) | yes | yes | **NO** | [T] |
| **`−m₁`, `−log…`** | lsc | yes | **fails** (sign) | yes | NO (sign) | [T] |
| **`m₁ = inf supp ν`** (support floor) | **usc, NOT lsc** (Prop 1 + witness) | yes | yes (`=0` gapless) | **no** (= answer) | **YES structurally; fails F4** | [O] |
| **`f₁² = ν({m₁})`** (edge atom) | **usc, NOT lsc** | yes | yes | **no** (= answer) | YES but = answer | [O] |
| **`1_{G_δ}`** (closed-set indicator, support target only) | **usc** (`G_δ` closed) | yes (binary) | yes | **no** (need `G_δ`) | YES for target (**G**); non-computable | [O] |
| **Fredholm minimal zero `z₁`** (`D_O(z₁)=0`, pole-below-cut) | **usc** (defined only on closed gapped set) | yes | yes (no isolated zero ⟹ `f₁²=0`) | **no** (= transport) | YES structurally; = transport | [O] |
| **pointwise limit of continuous `φ_k ↓`** (Baire-1; `inf_k φ_k`) | **usc** (inf of continuous) | yes **only in the limit** | only the limit is clean | limit = answer | YES iff the limit is taken | [T]/[O] |
| **lsc envelope `𝔉_*` of any candidate** | lsc | **no** (≤0 everywhere, Lemma 2i) | yes | yes | **NO** (vacuous) | [T] |

**Reading the catalog.**
1. **Every survivor is usc-not-lsc** (§2.1: the lsc class is empty). The dichotomy is *exhaustive and directional*: blind-and-lsc vs sensitive-but-non-computable.
2. **The four natural usc survivors are all "the answer in different coordinates."** `m₁`, `f₁²`, `1_{G_δ}`, `z₁` are mutually inter-derivable on the gapped class (given the Fredholm data, `z₁` gives `m₁` and the residue gives `f₁²`; `1_{G_δ}` is the thresholding of `m₁`). None is *computable* uniformly across the crossover: `z₁` evaluation = the `a→0` zero-finding = transport (`LOC-DETERMINANT §5`, `META-GAPFORCING §5`); `1_{G_δ}` needs `G_δ`; `m₁, f₁²` are the gap. So they survive **F3** but fail **F4** — exactly `SPECIES-EXISTENCE` Theorem 3 (equivalence/circularity), reached topologically.

### §2.4 The "limit of bounded-continuous approximants" — the universal shape. [T]

By the **Baire/Choquet representation** (every usc function on a perfectly normal space is an infimum of continuous functions; `𝔐` with vague topology is metrizable on bounded-mass sets, hence perfectly normal), `[R]`:
```
m₁ = inf_{λ} φ_λ ,     φ_λ continuous,   φ_λ ≥ m₁   (each φ_λ is an UPPER approximant).
```
So the usc floor is **literally `inf` of continuous upper bounds**. This is *exactly* the corpus's `r(t) ↓ m₁` (`TOOL-COMPLETELY-MONOTONE §2`), the Christoffel sequence, `M₀/χ` (squeeze upper side, `H1-IS-IT-THE-GAP §2`), etc. Consequences, all `[T]`:
- **Each finite `φ_λ` fails the acid test**: `φ_λ ≥ m₁`, so `φ_λ(free massless) ≥ m₁(free)=0`, and for the sharp families `φ_λ(free) > 0` strictly at finite `λ` (`r(t)` is positive at every finite `t` on the free field — `TOOL-COMPLETELY-MONOTONE §5`). The upper approximants are *positive where they must be 0*.
- **Only the infimum `inf_λ φ_λ` is acid-clean** (it returns `0` on the free field), and taking it is the **`λ→∞` limit** = `m₁` = the answer (Thm 1: the limit is not lsc; it is the usc floor).
- Therefore: **the entire "semicontinuity defect" `J` of a candidate is the gap between `inf_λ φ_λ` and its finite-level values; it lives purely in the limit.** The residual freedom (§3) is *which* family `{φ_λ}` — but the *value at the limit is universal* (`= m₁`).

This **unifies** the facet with the corpus's "only the `t→∞`/`k→∞`/`z→edge`/`a→0` limit is clean" observations: they are all instances of "`m₁` = inf of a sharp family of continuous *upper* bounds, and the infimum is a genuine usc upward spike attained only in the limit."

---

## §3 · The requirement spec (what this facet imposes) — and residual freedom

### §3.1 The hard constraint (R1–R4)

> **R1 (semicontinuity — non-negotiable).** `𝔉` must be **upper semicontinuous and not lower semicontinuous** on `(𝔐, τ)` (`τ ∈ {vague, weak, TV}`): its hypograph `{(ν,t): t ≤ 𝔉(ν)}` is closed, but its epigraph is not; equivalently `𝔉_* ≠ 𝔉` strictly on `G_δ`. — **[T]** (Lemmas 2, Thm 1).
>
> **R2 (vacuous-positivity separation).** `𝔉_* ≤ 0` everywhere (Lemma 2(i)); the certifier does **not** contain any "positivity part" that is itself `>0`. All certification is in the jump `J = 𝔉 − 𝔉_* ≥ 0`. (This is the formal content of "not a positivity/lsc bound.") — **[T]**.
>
> **R3 (spike location = the floor).** `supp J = {ν : 𝔉(ν) > 0} ⊆ G_{δ_𝔉} = {m₁ ≥ δ_𝔉}`, and `δ_𝔉` is **not** an independently-supplied constant: the *location* of the spike is `inf supp ν`. Equivalently, `𝔉` must **evaluate a singularity/edge of an analytic object, not a value of a continuous one** (`inf supp`, the Fredholm minimal zero `z₁`, the `t→∞` Tauberian rate, the support-boundary functional). This is the "singularity vs value" distinction, made precise: a *value* of a continuous functional is lsc/continuous (blind, Thm 1); a *singularity location* (edge/zero/pole of a transform) is usc (Lemma 2catalog). — **[T]/[O]**.
>
> **R4 (uniform attainment).** The usc spike must be attained **uniformly in `(L,a)` across the crossover `β ≈ 2.5`**: `inf_{F_YM} J ≥ c > 0` as `a→0, L→∞`. By §2.4 this is a statement about the *infimum over the approximating family* being `≥ c` — i.e. no member `φ_λ` of the sharp family dips below `c` at the crossover. — **[O] = the transport/H1 wall** (`TRANSPORT-FIVE-WAY-VERDICT §5`).

**Forbidden (what `𝔉` CANNOT be — this facet's no-gos):**
- **[NG]** any **continuous** functional of computable spectral data (moments, finite-`t` rate, invariants, continuous extremals) — Thm 1/Cor 1. (Kills all of C0 **at once**, not mode-by-mode.)
- **[NG]** any **lsc** functional satisfying F1 — Thm 1.
- **[NG]** any `Φ∘D` whose fibers contain a gapless element with the same data — Thm 1′ (kills discontinuous `Φ` too).
- **[NG]** any **finite** member `φ_λ` of an upper-bound family, used alone — acid-fails (§2.4).

### §3.2 Residual degrees of freedom (what the facet leaves open)

The facet pins the *semicontinuity class and the location of the spike*; it does **not** pin:
1. **The data fiber / carrier operator `O`** — which gauge-invariant interpolator and channel (`0⁺⁺` `tr F²` is the canonical one; the argument only needs some `O` with an edge atom). Equivalent spikes arise for any `O` with `f₁²(O) > 0`.
2. **The threshold `δ_𝔉`** — the level at which the spike is triggered; any `δ < m₁` works, and `δ_𝔉 → m₁` in the sharp limit. (Fixes only the *support* of `J`, not its value on it.)
3. **The sharp family `{φ_λ}` realizing the spike as `inf_λ φ_λ`** — Tauberian (`t→∞`), Christoffel/moment (`k→∞`), Fredholm-zero (`z → z₁`), continuum (`a → 0`), Baire/indicator approximants. **All give the same limit value `m₁`**; the freedom is the *parametrization of the limit*, not its content.
4. **The quantitative flavor** — `m₁` (floor) vs `f₁²` (atom residue) vs `1_{G_δ}` (binary support). The support target (G, = Clay) needs only `1_{G_δ}`'s content (the closed-set membership); the atom target (H\*) needs `f₁²`. Both are the same `J` read at different coarseness (`SUPPORT-NOT-ATOM`).

> **The content vs parametrization split (the facet's sharpest sentence).** The facet forces `𝔉` to live entirely in a **usc upward spike whose location is `inf supp ν`**. The *value* of the spike is forced (in the sharp limit it is the floor/atom). The *residual freedom is parametrization only* — which family, which carrier, which threshold — and **every parametrization's limit is the same floor.** Hence this facet **cannot manufacture** the species by clever choice of data or regularization: it reduces the entire freedom to "pick a representation of `inf supp ν`." The only remaining content is **R4**: attaining that usc limit *uniformly across the crossover*, which is verbatim the transport machine (`CLAUDE.md §5`, "尚不存在").

### §3.3 The exact missing step this facet isolates [O]

```
MISSING STEP (semicontinuity coordinate):  a usc functional 𝔉 with lsc envelope ≤0
   (Lemma 2) whose upward spike J is (i) supported on {m₁ ≥ δ} with spike-LOCATION = inf supp ν
   (R3), and (ii) evaluable UNIFORMLY in (L,a) across β≈2.5 without the spike-location itself
   being supplied (R4).           Equivalently:  a way to ATTAIN the usc limit inf_λ φ_λ
   uniformly, when each φ_λ is computable but acid-failing and only the limit is clean.
```
This is **not** "find a smarter continuous functional" (Thm 1 forbids) and **not** "find an exact identity" (`READ ≠ FORCE`; a `J` equal to `m₁·1_{G_δ}` is an identity-read = F4 fails). It is **"make the semicontinuous limit uniform"** — the transport/zero-finding problem in its semicontinuity coordinate. The facet thus **certifies a clean no-go** (continuous/lsc/gapless-fiber functionals are blind, `[T]/[NG]`) and **localizes** the open content to R4 (`[O]`), without claiming Clay and without a QED.

### §3.4 Honest boundaries. 

- **`∀class ≠ ∀possible`.** Thm 1/1′ and the catalog are exact for the named classes (continuous/lsc/gapless-fiber functionals of the renormalized KL measure). They do **not** prove no `𝔉` whatsoever exists: a functional on a data fiber where the gapless class is **not** dense could escape Thm 1′ *without* being lsc. **The candidate such fiber is the full field-space measure** (not just `ν`): the filler `(1−ε)ν + εδ_0` is an admissible GFF, but it is **not** a limit point of interacting YM in field-space, so field-space data fibers may not be gapless-dense. That escape is precisely the loophole **L-2 ≡ H1** of `E2-FLIP-NOGO §5` (a field-space nonlinear anti-concentration inequality). By R1, even a field-space `𝔉` would have to be **usc-not-lsc** — so the semicontinuity requirement is **necessary for every data mode**; what changes is only whether the *density* (hence vacuity vs non-vacuity) holds. Proving a usc floor survives on the field-space fiber **is** H1. — **[T] for the barrier + [O]/[J] for the boundary.**
- **The "computable ⟹ continuous" step** is `[J]` for arbitrary computable `Φ` (a discontinuous computable `Φ` exists); the *rigorous* killers are Thm 1 (lsc) and Thm 1′ (fiber). The *identification* "constructed/analytic `𝔉` ⟹ continuous ⟹ blind" is `[T]`; that "all C0 modes are continuous/lsc" is `[T]` (each mode checked). 
- **No new floor is manufactured here.** This facet's contribution is the **decomposition `𝔉 = 𝔉_* + J`,** the **topological proof of `F3∧F4=∅`** (`F4`⟹continuity⟹`Thm 1`), and the **identification of the residual freedom with the parametrization of `inf supp ν`** — a sharpening and unification, not a solution.

---

## References

**External — [standard]/[R] (textbook; not re-verified this session):** 
- P. Billingsley, *Convergence of Probability Measures* (2nd ed., Wiley 1999) — Portmanteau theorem, vague/weak convergence, Prokhorov. 
- G. Choquet, *Topologie* / *Lecture on Analysis*; R. Engelking, *General Topology* — usc/lsc, closed hypograph ⟺ usc, usc = inf of continuous on perfectly normal spaces (Baire class 1/2). 
- A. Källén (1952)/H. Lehmann (1954); K. Osterwalder, R. Schrader, *CMP* **31** (1973) 83; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440 — KL measure, OS/RP. 
- A. Jaffe, *Phys. Rev.* **158** (1967) 1454 — generalized free fields (`ν ≥ 0 ⇒` OS-admissible). 
- S. Bernstein (1928/29); D. V. Widder, *The Laplace Transform* (1941) — completely monotone ⟺ `ν ≥ 0`. 

**Internal (absolute paths; the banked results this facet assembles):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall/key/machine "尚不存在").
- `SPECIES-EXISTENCE.md` — class `ℭ` / F1–F4, `F3∧F4=∅` (Theorem 2), equivalence/circularity (Theorem 3), `F4=C0`. **This facet's §1.5 is a topological proof thereof.**
- `META-GAPFORCING-IDENTITY.md` — READ (Fredholm/B6) ≠ FORCE; every exact identity admits the gapless free field. → R3 "singularity vs value".
- `SUPPORT-NOT-ATOM.md` — Clay = support floor `inf supp ν`, atom-free gapped `:φ²:`; the cone is `H1`-or-`f₁²`.
- `H1-IS-IT-THE-GAP.md` — squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ`; `M₀/χ` and `r(t)` are *upper* bounds (the `φ_λ` of §2.4); thin-atom witness `wδ_ε+(1−w)δ_S` (§3(b)) = the gapless filler.
- `PROOF-NEVANLINNA-POLE.md` — atom **upper-semicontinuous (portmanteau)**; germ fiber convex with residue ranging `[0,λ_N]`; Mechanism 2 "continuous germ-functional cannot lower-bound a `0↔λ_N` jump". → Prop 1 / Thm 1.
- `TOOL-COMPLETELY-MONOTONE.md` — `CM ≡ ν≥0`; `r(t) ↓ m₁` is a decreasing *upper* bound, positive on the free massless field at every finite `t`; only `t→∞` is clean. → §2.4.
- `LEVER-nonconvex-edge.md §0B` — linear-data lemma (finite moments cannot pin the edge). → §1.4.
- `E2-FLIP-NOGO.md §2,§5` — four-horn dichotomy; loophole **L-2 ≡ H1** (field-space anti-concentration). → §3.4 boundary.
- `REDUCTION-AUDIT.md (b)`, `GAP-G2-edge.md` — thin-atom witness; moment-ratio upper bound; residue–level-gap lemma.
- `LOC-DETERMINANT.md §5`, `INVENT-NEW-IDENTITY.md` — Fredholm minimal zero = READ; `a→0` zero-finding = transport. → catalog `z₁`.
- `TRANSPORT-FIVE-WAY-VERDICT.md §5` — the transport machine "尚不存在"; uniform crossing = H1. → R4.
- `PROOF-P3-UNCONDITIONAL.md` — computable `f₁²>0` at `β<1/48` (cluster expansion; diverges before crossover — the sub-regime member of the usc class, F4 holds *only* there).

*Status: **[T]/[NG] THE OBSTRUCTION (floor & atom usc, gapped closed-nowhere-dense, gapless dense; Thm 1 lsc no-go, Thm 1′ fiber no-go, Cor 1 continuous no-go; topology-robust vague/weak/TV) + [T] the DECOMPOSITION `𝔉=𝔉_*+J` (Lemma 2, vacuous-positivity part `𝔉_*≤0` + usc spike `J`=species) + [V-here] topological proof of `F3∧F4=∅` and exact two-atom moment witness `½δ_1+½δ_3` vs `⅕δ_0+⅘δ_{5/2}` sharing `(M₀,M₁,M₂)=(1,2,5)` + [O] the exact missing step = uniform attainment of the usc limit across the crossover (= transport/H1). REQUIREMENT: `𝔉` usc-not-lsc, lsc-envelope ≤0, spike located at `inf supp ν` (R1–R4). RESIDUAL DOF: parametrization of the spike (fiber/carrier, threshold, sharp upper-bound family, quantitative flavor) — all limits = `m₁`; content = spike location = the gap. Named lsc-and-floor-sensitive class is EMPTY. NOT a proof of Clay; NOT a global impossibility theorem (`∀class≠∀possible`; field-space escape = L-2 ≡ H1). Clay stays 0%. No QED.***

# PROFILE-FACET-ANCHOR — the **SELF-ANCHORING** facet of the wanted-species profile `𝔉`: what structural property lets a lower bound on `m₁`/`f₁²` be produced with **no controller** and **without running backward across the crossover**

> 2026-11-24 · SPEC facet (charter §1/§2). One facet of the "profile" (wanted-species specification); sibling facets are nonlinearity / non-positivity / IR-support-sensitivity / gapped-vs-gapless discrimination / `(L,a)`-uniformity (`META-GAPFORCING-IDENTITY §4`, `SPECIES-EXISTENCE §1`). **This note is a PRECISE STRUCTURAL SPEC, not a slogan and not a promise.** It does **NOT** prove Clay; it does **NOT** prove no self-anchoring tool can exist (`∀class ≠ ∀possible`). Labels: **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for a *named class* · **[J]** 判断 · **[V-here]** new this round · **[R]** inherited/banked (label carried).

---

## VERDICT (the structural requirement + the residual DOF)

> **A functional `𝔉` is SELF-ANCHORING iff its lower bound `Δ(𝔉,H)` is the *solution of an intrinsic relation* `R(𝔉,H) ≥ 0` that (a) uses **no controller** — no fixed continuum edge, no separation oracle `E₂`, no quenched ensemble, no rigid saddle — and (b) **is violated by the gapless (free) field**, so that `R`'s own solution set lies in `[Δ₀,∞)` with `Δ₀>0`.**
>
> **Formally, self-anchoring is always available** — the intrinsic schemas (contraction fixed point, order/cone fixed point, gap equation, eigenvalue inequality, bootstrap feasibility) each **escape controller-necessity by construction** (§2, all `[R]` banked). **But every one of them is defeated by the SAME self-anchoring-specific obstruction that is NOT controller-necessity:** the intrinsic relation has a **gapless solution** (free field / neutral fixed point / order-reversing map / trivial feasible set), and the only way to remove it is to make the relation **nonlinear, IR-support-sensitive, and `k→0` sign-definite** — which is the gap (`H1`) in yet another coordinate.
>
> **Structural requirement on `𝔉` (the spec): `𝔉` must be the fixed point / minimum-norm solution of a *nonlinear, IR-support-sensitive, controller-free* relation `R`, and that relation must *intrinsically exclude* the gapless solution by a `k→0` sign-definite mechanism.** Equivalently: `𝔉` must be **self-anchoring** (no external edge) **and floor-sensitive** (sees `inf supp ν`) **simultaneously**, and the two are jointly satisfiable only by a `k→0` sign-definite self-consistency. **[NG] for every structural relation; [O] for the unique floor-sensitive survivor (Fredholm zero).**
>
> **Residual degree of freedom: the `k→0` sign-definite self-consistency mechanism** — the single slot Theorem SA (§2.6) leaves open. The *representational* choices (which nonlinear transform, which coordinate `G″`/`f₁²`/`E₂−m₁²`/Fredholm-zero, which term decomposition `G″=S−D`/QV-sum) are free; the *content* of the exclusion is invariant and equals `H1 = I2 = H*`. Its two faces — **exclusion of the gapless fixed point** and **`a→0` zero-finding uniformity** — are one object, the transport wall. Clay stays open; 0%.

---

## §0 · Setting, and the exact dichotomy the facet must respect

**Object (`BRIEF §1`, `SPECIES-EXISTENCE §0.1`).** `H` = self-adjoint transfer Hamiltonian on `H_phys=ker Gauss`, unique vacuum `Ω`, `HΩ=0`; on `(vac)^⊥`, bottom `E₁=m₁`, second level `E₂`. `ν≥0` = KL spectral measure, `supp ν ⊂ [m₁²,∞)`, atom `f₁²=ν({m₁²})`; `χ=∫s⁻¹dν`; squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` `[T]`. `F_YM` = lattice family `(β,L,a)`, `β` crossing `~2.47`. Two targets: **(G)** support floor `m₁≥Δ` (Clay proper); **(H\*)** atom `f₁²≥c` + isolation `E₂−m₁²≥δ` (strictly stronger; `:φ²:` is gapped-atom-free).

**Control parameter `λ`.** Family `{H_λ}` with `λ ∈ I`; `λ_ctrl` = controlled (more-gapped) end; `λ_phys` = physical/crossover end. Four instantiations: volume (`λ=1/L`; `λ_phys=0`), coupling (`λ=β`; `λ_phys=2.5`), deformation (`λ=m₀`; `λ_phys=0`), compact radius (`λ=1/r`; `λ_phys=0`).

**The two failure modes this facet must navigate (`SPECIES-EXISTENCE §5`, exact):**
- **Controller-necessity:** a lower bound that compares to a reference/controller; in YM the controller *is* the gap (`H1`).
- **Min–max backward:** a monotone/order lower bound whose controlled end sits on the *more-gapped* side ⇒ every proven inequality caps, not floors.

**Self-anchoring** is the conjunction of two negations: **no controller** (§1) **and** **not backward** (§1.3). The facet asks whether a functional can supply a floor with *both* negations met — i.e. with **no comparison at all** to an external anchor and **no monotone transport** from a gapped end.

---

## §1 · Controller-necessity and the backward-Lemma as theorems

### 1.1 The two objects being formalized

> **Definition (controller).** Let `D(H)` be the **accessible data** of `H` (Ward/moment/DSE data, finite-`t` correlators, deformation invariants, quantitative extremals — the class `C0` of `E2-FLIP-NOGO §1`). A **control datum** for a bound `Φ` is `C(H)` with a **validity condition** `(VC)`, such that `Φ=Φ(D,C)` is sound (`Φ≤m₁`) *only under* `(VC)`, and `(VC)` is **not** a consequence of `D(H)` alone: varying `C` at fixed `D` can flip the sign of `Φ−m₁`. The bound is **controller-free** iff it is sound with `C≡∅` (no `(VC)` beyond `D(H)`).

> **Definition (self-anchoring).** `𝔉` is **self-anchoring on `F`** iff there is a closed **intrinsic relation** `R(𝔉,H)` on accessible data such that **(SA-sound)** `R(𝔉,H) ⟹ 𝔉(H) ≤ m₁(H)` (or `≤f₁²`) for every `H∈F`; **(SA-nonvacuous)** `R` holds with `𝔉>0` only for interacting/gapped `H` and is *empty or `≤0`* at the free/gapless `H`; **(SA-floor)** the solution set `{Δ : R(Δ,H)} ⊂ [Δ₀(H),∞)` with `Δ₀≥δ>0` uniform on `F`, computed from `R` and non-edge data.

### 1.2 Theorem C — order/controller duality (which schemas force a reference). `[T]`

> **Theorem C.** Let `A,H` be self-adjoint, bounded below, `A ≤ H` in quadratic-form sense on a common form domain. Then `E_k(A) ≤ E_k(H)` for all `k` (Courant–Fischer / min–max monotonicity `[R]`; `E2-FLIP-MINORANT §0`). Hence:
> **(i) lower-bound-by-comparison** `E₁(H) ≥ E₁(A)` is sound **iff the order oracle `A ≤ H` holds** *and* `E₁(A)` is independently computable;
> **(ii)** if only `H ≤ A` holds (the reference is *more gapped*), min–max gives `E₁(H) ≤ E₁(A)` — an **upper** bound;
> **(iii)** if neither order is provable, comparison yields **nothing**.
> *Proof:* `E_k` is monotone non-decreasing under form order (same trial space, pointwise `⟨f,Af⟩≤⟨f,Hf⟩`); the assertions are that monotonicity read in the three cases. ∎

**Corollary C′ (why every comparison reference is on the wrong side of YM). `[T]`**
A comparison lower bound requires the reference to be a **form-minorant** of `H` — the *floppier* (less-gapped) side — with a computable bottom. In YM:
- the **free field** is the only universally valid minorant candidate, but (a) `H_free ≤ H_YM` fails wherever the interaction is sign-indefinite — the Savvidy fiber `Hess = −D²+2ad F_B`, `λ_min=−g|B|<0` `[T]` (`PROOF-SAVVIDY-EXACT-SIGN`), so **no order** exists there; (b) even if `V≥0`, `E₁(H_free)=0`, giving the vacuum bound `m₁≥0` and failing the acid test (certifies `0`). So the free reference is either **no order** or **`≤0`**.
- every **proven-gapped reference** (strong coupling, small volume, mass deformation) is on the **more-gapped** side: `H_YM ≤ H_ref`, not `H_ref ≤ H_YM`; Theorem C(ii) returns an **upper** bound.

### 1.3 Theorem B — backward monotone transport (why the free inequality is always an upper bound). `[T]`

> **Theorem B (backward monotone transport).** Let `m(λ)=E₁(H_λ)` on `(vac)^⊥`, and suppose the **certified monotonicity** holds: `m(λ_phys) ≤ m(λ_ctrl)` (the gap is non-increasing toward the physical/crossover end). Then any certified lower bound at the controlled end, `m(λ_ctrl) ≥ δ>0`, implies only
> ```
>      m(λ_phys) ≤ m(λ_ctrl)           (an UPPER bound on the physical gap);
> ```
> a lower bound `m(λ_phys) ≥ δ` is **not implied** by monotonicity alone.
> *Proof:* `m(λ_phys) ≤ m(λ_ctrl)` is the hypothesis; a lower bound at `λ_phys` is logically independent of it. ∎

> **Lemma B″ (the variational root cause — the "free inequality" points up). `[T]`**
> `E₁(H) = inf_{ψ⊥Ω,‖ψ‖=1}⟨ψ,Hψ⟩` is an **infimum over a variational class**. Every controlled anchor is obtained by **restricting** that class (compactify / fewer trial states), **up-shifting the operator** (`+m₀²`, strong coupling — small lattice β — localizing the measure), or **imposing more boundary constraint**. By Courant–Fischer both *restriction* and *up-shift* **raise** the inf, so `E₁(controlled) ≥ E₁(physical)`; the inequality an order/min–max argument can prove runs **backward**. **The direction is forced by the variational characterization, not by a choice of reference.** The three named monotonicities (each `[T]` in its family):
> - **Volume** (`λ=1/L`, min–max over a larger trial space): `Δ(ℝ⁴)=inf_L Δ(L)`, `Δ(L₂)≤Δ(L₁)` for `L₂>L₁` `[T]` (`SMALLCIRCLE-ADIA §3` Lemma; `P4: Δ_{ℝ⁴}=inf_V Δ_V`).
> - **Coupling** (`λ=β`): lattice gap `m(β)↓` in `β` (AF `√σ a²=−log(β/2N²)`), strong coupling more-gapped `[T]` (`MONOTONE-COUPLING §2`).
> - **Deformation** (`λ=m₀`): `gap(m₀)↑` in `m₀`, `gap(YM)=inf_{m₀}gap(m₀)` `[T]` (`MONOTONE-COUPLING §4`).
> In **every** case the controlled end is the more-constrained / more-gapped end, so **the free inequality is an upper bound**. ∎

> **Corollary B‴ (no scalar ODE shortcut). `[T]`** The derivative form fails identically: `m₁(β)=−(1/a)ln λ₁(β)` is a spectral-edge datum, `dλ₁/dβ=⟨1|∂_βT_β|1⟩` needs the eigenvector and the full matrix (no `F(λ₁,β)` exists); the only rigorous β-Grönwall bounds (Weyl/Lipschitz, Holley–Stroock) have **extensive** constants `≍V` and `→0` as `V→∞,a→0` (`PROOF-GAP-DIFFERENTIAL-INEQ §1,§3,§5`). The backward transport cannot be closed by an ODE either.

### 1.4 Theorem ENCL — controller circularity at the crossover (`SSC`). `[T + C-cond]`

> **Theorem ENCL (minorant ⇒ separation ⇒ `E₂`). `[T]`** For normalized `ψ`, `R=⟨ψ,Hψ⟩`, `Var=‖(H−R)ψ‖²`, the Temple bound `m₁ ≥ R − Var/(ρ−R)` is sound iff `ρ ≤ E₂`, and is strictly increasing in `ρ`; **no function of `(R,Var)` alone lower-bounds `m₁`** (impossibility lemma: fix `R,Var`, move an atom to `−N`) `[R]` (`E2-FLIP-MINORANT §1(i)`; `TEMPLE-WEINSTEIN §3`).
> **SSC (`Single-Scale Covariance`). [C — physical single-scale assumption]** pure YM has one dynamical scale, so `E₂ = ρ_r·m₁` with `ρ_r=m_{0⁺⁺*}/m₁` a fixed `β,L,a`-independent number (`G2D-ISOLATION-I2 §2`; `E2-REGRESS §0.4`). Thus `ρ ≤ E₂ = ρ_r m₁`: **the separation oracle is a fixed multiple of the very `m₁` being bounded** ⇒ the minorant *is* the gap, one rung up. **`[NG]` for the enclosure class; the exact input `ρ≤E₂` is `I2=H*`.**

**§1 verdict (`[T]`+`[NG]` for named classes).** Four theorem-schemas **force** a controller:
1. **Comparison/order** ⇒ order oracle `A≤H` with `A` on the floppier side (Theorem C) — in YM the only floppier reference is free (`≤0`, acid-vacuous), all gapped references are on the wrong side (Corollary C′).
2. **Min–max enclosed lower bounds** (Temple/Kato/Weinstein/Gershgorin/Brauer) ⇒ separation oracle `ρ≤E₂` = `ρ_r m₁` (Theorem ENCL, SSC) — circular.
3. **Monotone parameter transport** (β, r, m₀, L) ⇒ a controlled more-gapped end (Theorem B, Lemma B″) — backward.
4. **Semiclassical splitting / FMM** ⇒ rigid saddle (no size modulus) / quenched ensemble (`META-CROSSDOMAIN-PROTOTYPE §3,§2`) — at the crossover this is an uncontrolled `O(1)` liquid, `=H1`.
The **common core:** a comparison lower bound needs a reference on the **dominated** side with a computable bottom; YM's only dominated reference is the free field, and every gapped reference is on the **more-gapped** side. **Controller-necessity is the statement that a comparison floor must import an anchor that, in YM, is the gap.**

**Honest boundary `[J]`.** Theorem B and Corollary C′ are **not** a proof that no lower bound exists — they are a proof that **no comparison/order/min–max/transport lower bound** can cross the crossover (`∀class ≠ ∀possible`). A floor must therefore be **controller-free** (§2).

---

## §2 · The class of self-anchoring schemas — which escape controller-necessity, and which survive

### 2.0 The class, and the escape criterion

By §1, escaping controller-necessity means using **no** `A≤H` oracle, **no** `ρ≤E₂`, **no** external ensemble/saddle — i.e. the bound is the **fixed point / solution of an intrinsic relation** `R(𝔉,H)`. Five canonical schemas exist in the corpus. Each is **controller-free by construction** (hence escapes Theorem C/ENCL/B). The question is whether each satisfies **SA-floor** (`Δ₀>0`) and **SA-nonvacuous** (fails on free) — i.e. whether `R` **excludes the gapless solution**.

### 2.1 Schema A — contraction fixed point (Banach). `[NG]` for the application; escapes controller-necessity

`𝔉` = unique fixed point of an explicit contraction `T` (`ρ(T)<1` uniformly); computable by iteration, no external reference. **SA-sound holds vacuously** (no anchor needed).
**Failure (`BANACH-FIXED-POINT §2–3,§5`, `[T]`).** The physical map (constrained curvature Riccati / FRG) has `ρ(DT)→1` at `k→0` (**zero mode neutral**) and `ρ(DT)=∞` on the Savvidy shell (`(k∂_kR)/(R−2g|F_B|)²` blows up as `R_k↓2g|F_B|`) ⇒ **no uniform contraction**, and `{κ≥λ}` is **not `T`-invariant**. The only contractively-certified fixed point near `κ=0` is the **gapless** one (`T(0)=0`); the nonzero branch is a **bifurcation** (eigenvalue of `DT` crossing `1`) whose insertion *is* the order parameter `=H1` (circular).
**Verdict:** escapes controller-necessity; defeated by **neutral/expansive direction** — the gap lives on the *neutral* `k→0` mode of `T`, unreachable by contraction; Savvidy makes `T` expansive exactly where the gap must be read. `[NG]`.

### 2.2 Schema B — monotone / order fixed point with invariant cone (Schauder–Krasnosel'kii). `[T]` at β=0; `[NG]` at the crossover

`{𝔉≥Δ₀}` an invariant set of a monotone map `T` (order-preserving in a cone), nonempty ⇒ fixed point `≥Δ₀`; controller-free (order-theoretic, no external reference).
**Survivor (`BANACH-FIXED-POINT §4`, `[T]`).** At β=0, `W=Φ″⪰N/12` is convex, Prékopa makes the marginal log-concave, `{κ≥N/12}` **is** `T`-invariant, Schauder/Krasnosel'kii give `κ*≥N/12`, uniform in `L` — a genuine self-anchored bound, **β=0 only**.
**Failure at β=2.5 (`§3`, `[T]`).** `{κ≥λ}` is **not `T`-invariant** (a `κ≥λ` profile maps to `κ<0` on the typical-field tail; log-concavity fails); no invariant order cone for SU(N) (torsion, `C8`); `Cov(G,Φ′)>0` is order-**reversing**; `W` indefinite on the cone.
**Verdict:** escapes controller-necessity; defeated by **order-reversal / non-convexity** (the indefinite Savvidy Hessian destroys the invariant cone). `[NG]` at the crossover.

### 2.3 Schema C — gap equation / one-sided fixed point (`Δ ≤ F(Δ)` forcing `Δ ≥ Δ₀`). `[T]` at β<1/48; `[NG]` at the crossover

This is the schema the task names. Two equivalent forms:
- **Equality (gap equation).** `Δ = F(Δ)` with `F` intrinsic (Ward–Stein/DSE closure). **[T] (`MS-GAP-EQUATION §1–2`)** the exact self-consistent equation `κ = F(κ) := (1−η(κ))/Var(a;κ)` holds for every measure; the **gapless branch `κ=0` is always a solution and is the contraction attractor** (`F(0)=0`, free field). A positive fixed point exists **iff `F′(0)>1`** (a pitchfork = linear IR instability).
- **Inequality (forcing).** `∀Δ:[Δ ≤ F(Δ)] ⟹ Δ ≥ Δ₀` is the containment `{Δ : Δ ≤ F(Δ)} ⊆ [Δ₀,∞)`, which needs **two** facts: **(f-i)** the origin is **excluded**, `Δ=0 ∉ {Δ : Δ ≤ F(Δ)}`, i.e. `F(0) < 0` (otherwise `0 ≤ F(0)` puts `Δ=0` in the set and forces `Δ₀=0`); **(f-ii)** an **upward crossing** at `Δ₀` (`F(Δ)<Δ` just below, `F(Δ₀)≥Δ₀`). Fact (f-i) is the sharp one: since the physical/free relation has `F(0)=0` (gapless), the forcing **collapses to `Δ₀=0`** unless `F(0)<0` — i.e. unless the relation is *violated at the gapless point*. So both forms demand the **same** thing through different faces: the equality form's `F′(0)>1` (pitchfork) and the inequality form's `F(0)<0` (exclude the origin) are both **"`R` must fail at the gapless solution."**

**Failure at the crossover (`MS-GAP-EQUATION §3`, `[T]`).** Both faces require the gapless point to be removed: `F′(0)>1` (equality pitchfork) or `F(0)<0` (inequality exclusion). Both need a **`k→0` Savvidy tachyon** to make the coarse self-energy negative at `k=0`. It is **absent**: the tachyon is **transverse-Landau-gapped** (`λ_{n,s}=k_∥²+(2n+1)g|B|+2s gB`, lowest level `k_∥=0` but gapped at `√(g|B|)`); at β=2.5 the coarse `k→0` plane-wave block is **positive** (`+0.40β`, 0% negative), and 90–97% of the negatives are **high-momentum** (`W6-BE-at-2.5`); a fast-mode gap cannot gap the slow mode (`PROOF-IR-H1 §4`). Hence `F′(0)≤1` **and** `F(0)=0` (the free/gapless relation) at the coarse mode, and the **provable branch collapses to the gapless one**. The Ward moment `∫a m p=1` is a *global normalization* compatible with `m′(0)<0` — **moment ≠ atom**; the sign of `κ` lives in the IR curvature and is unconstrained by the intrinsic moment data.
**Verdict:** escapes controller-necessity (genuine self-consistency); defeated by the **gapless equality/fixed-point solution** — the free field sits on `R` (`F(0)=0`), so SA-floor fails. The escape (a `k→0` sign-definite instability forcing `F′(0)>1` / `F(0)<0`) is exactly `H1`. `[NG]` at the crossover; the genuine positive branch exists only at β<1/48 (P3) and β=0 (Haar `N/12`).

### 2.4 Schema D — spectral gap as solution of an eigenvalue-inequality (Courant–Fischer-as-definition). `[O]` / tautology

`m₁ = sup{λ : H ≥ λ P_⊥}` (Rayleigh characterization). This is a *definition*; used as a lower bound, it requires proving `H ≥ λP_⊥` for explicit `λ>0` from non-edge data. But `H ≥ λP_⊥` **is verbatim the gap statement** (Theorem POS, `E2-FLIP-MINORANT §1(v)`; `PROOF-COARSE-POINCARE`). 
**Verdict:** the schema "escapes controller-necessity" only by **being the gap**; it restates the target. Its SA-sound relation `R=(H≥λP_⊥)` has SA-floor iff `R` is *provable*, which is the gap. **No independent content; the inequality's provability is the whole problem.** `[O]`/tautology.

### 2.5 Schema E — bootstrap / moment-feasibility (island from consistency). `[NG]`

`Δ` must lie in the intersection of a set of exact constraints (OS positivity, Ward, sum rules, crossing). **`[NG]` (`A4-SPECTRAL-BOOTSTRAP §3`)** the feasible set is **nonempty and contains the free two-photon-cut field** (`f₁²=0`, gapless, implements every protected moment — `SCHUR-CONVEX`) ⇒ the island **contains `Δ=0`** ⇒ SA-floor fails. Forcing `Δ₀>0` (excluding every gapless measure) requires adding the **excluded spectral-edge datum `Λ`** (the residue–level-gap lemma `f₁²≥(ΛM₀−M₁)/(Λ−m₁²)`) — which is the coarse-mode Poincaré mass `=H1` renamed (`A4 §3`, `GAP-G2-edge`).
**Verdict:** escapes controller-necessity (consistency, no anchor); defeated by the **gapless-feasible-point** — the free field is admissible, so the feasible set does not floor. Same shape as Schema C, in the moment coordinate.

### 2.6 Theorem SA — self-anchoring rigidity: the **gapless-solution obstruction**. `[T-schema]` + `[NG]` for the structural class + `[O]` for the survivor

> **Theorem SA.** Let `R` be an **intrinsic relation** defining a self-anchoring functional (SA-sound). Suppose `R` is **closed under the free/gapless field** — i.e. the free field `ν_free` (and a gapless interacting `ν_CFT`) satisfies `R` with `𝔉≤0`. (This holds for *every* structural relation, `[T]`: Banach `T(0)=0`; Schauder β=0 cone `∋0`; gap equation `κ=0`; bootstrap feasible `∋ν_free`; eigenvalue inequality at `λ=0`.) Then `{Δ:R}` is **not** contained in `[Δ₀,∞)` for any `Δ₀>0` — **SA-floor fails**. To meet SA-floor, `R` must **exclude the gapless field** (`R` violated by some `inf supp ν=0`, `ν∈𝔐`). But `R` built from structural data (Ward/IBP, Stein, KL moments, OS positivity, DSE, eigenvalue identities, OPE, Bianchi) is satisfied by the free/gapless field — **`META-GAPFORCING-IDENTITY §3.1 Theorem`**, the "identity-forcing is impossible" no-go. Hence SA-floor can be met **only** if `R` is a **non-structural, nonlinear, IR-support-sensitive** relation that the gapless field *violates*; and by `SPECIES-EXISTENCE Theorem 2` (`F3∧F4=∅`) and `META-GAPFORCING §4`, the only such object is the **Fredholm zero / edge**, whose evaluation is `m₁=inf supp ν` `=H1`.
> ∎ (`[NG]` for the structural-relation class; `[O]` for the unique survivor.)

**Reading.** Theorem SA is the **self-anchoring face** of the corpus's `computable ⇒ floor-blind / floor-sensitive ⇒ circular` dichotomy (`SPECIES-EXISTENCE §5`). The facet's contribution: it isolates *why* self-anchoring fails for structural schemas — not because they need a controller (they don't), but because **the free field is always a solution of a structural relation**, and removing it is the gap. **The "anchor" a controller would have supplied is replaced, in the self-anchoring schema, by the *exclusion of the gapless solution* — and that exclusion is `H1`.**

### 2.7 The unique floor-sensitive self-anchoring object — and why it is still `H1`

The **Fredholm determinant** `D_O(z)=det(1−zT_O)`, atom `f₁²=−(1/z₁)N_O(z₁)/D_O′(z₁)`, `z₁=e^{am₁}` the **smallest zero** (`I7`; `docs/new-tools-registry.md §A3`; `META-CROSS-DOMAIN §1`), is the **unique** object that is simultaneously:
- **self-anchoring** — `z₁` is a root of an *intrinsic* spectral equation `D_O(z₁)=0`, no external controller (`INVENT-NEW-IDENTITY`);
- **nonlinear in `ν`** — not a moment (linear data lemma, `LEVER-nonconvex-edge §0B`);
- **IR-support-sensitive** — gapped: `z₁` a pole below the cut; gapless: the zero merges into the cut (`META-GAPFORCING §4.5`, `LOC-DETERMINANT §5`);
- **acid-clean** — free field: pure cut, no isolated zero, `f₁²=0` (`META-GAPFORCING §3.2 W3`).
So it is the **only known self-anchoring + floor-sensitive relation** — the facet's goldilocks object. **But it is `READ`, not `FORCE`** (`META-GAPFORCING §3.5`): to *lower-bound* `f₁²` one must (i) locate `z₁` (bracket the zero `=know m₁`) and (ii) upper-bound the kernel derivative `D_O′(z₁)` — both need control at `α_s~O(1)`, and making them `(L,a)`-uniform as `a→0` is the **marginal BPHZ zero-finding wall** `=transport=H1` (`META-CROSSDOMAIN §1.4`; `LOC-DETERMINANT §5`). Its consistency condition is the isolation `S₂=o(λ₁)` `=I2`. **Verdict: `[T-shape]` (self-anchoring + floor-sensitive in form) but `[O]` (uniform evaluation = transport).**

### 2.8 §2 synthesis table

| schema | controller-free? | gapless solution of `R`? | survives `(L,a)` crossover? | label |
|---|---|---|---|---|
| **A** contraction FP (Banach) | **yes** | yes (`T(0)=0`, contraction near 0 certifies the gapless FP) | **no** — neutral `k→0` mode; Savvidy expansive | `[NG]` |
| **B** order/cone FP (Schauder) | **yes** | yes (β=0 cone `∋0`; and no invariant cone at 2.5) | **no** — order-reversal; `Cov>0`; torsion | `[T]` β=0 / `[NG]` 2.5 |
| **C** gap equation `κ=F(κ)` | **yes** | yes (`κ=0` always; free `F≡0`) | **no** — needs `F′(0)>1` = `k→0` tachyon (absent) | `[T]` β<1/48 / `[NG]` 2.5 |
| **D** eigenvalue inequality | **yes (tautology)** | =the gap | — | `[O]`/tautology |
| **E** bootstrap feasibility | **yes** | yes (`ν_free` feasible) | **no** — forcing needs `Λ`=H1 | `[NG]` |
| **Fredholm zero** | **yes** | **no** (unique floor-sensitive) | **no** — `a→0` zero-finding = transport | `[T-shape]`/`[O]` |

**Synthesis [T]/[NG]/[O].** *Every* self-anchoring schema escapes controller-necessity; *none* escapes the **gapless-solution obstruction** except the Fredholm zero; and the Fredholm zero's uniform evaluation is the transport wall. So the facet's exact result: **self-anchoring is necessary for a crossover floor but not sufficient; the binding residual is the `k→0` sign-definite exclusion of the gapless solution, which is `H1`.**

---

## §3 · The requirement spec on `𝔉` (and the residual freedom)

### 3.1 The structural requirement (what `𝔉` must be)

A self-anchoring `𝔉` must satisfy **SA-sound + SA-nonvacuous + SA-floor** (§1.1) with an **intrinsic, controller-free** relation `R`. Theorem SA pins the only viable shape:

> **(SA-1) Controller-free.** `R(𝔉,H)` uses only data of `H` — **no** `E₂`, fixed continuum edge, quenched ensemble, or rigid saddle (`C≡∅`). [escapes §1]
> **(SA-2) Nonlinear in `ν`.** `R` is **not** a linear/moment/positivity functional (linear data are floor-blind; `LEVER-nonconvex-edge §0B`, `SCHUR-CONVEX`). [necessary: structural relations are gapless-closed, Theorem SA]
> **(SA-3) IR-support-sensitive.** `R` **distinguishes gapped from gapless** — it is violated by the free/gapless field (acid-clean), not merely "satisfied differently" (`META-GAPFORCING §3.3b`, `IR-RESPONSE-TEST`).
> **(SA-4) Gapless-solution excluded intrinsically.** The solution set `{Δ:R}` must **not** contain `Δ=0`; the exclusion must be done by `R` itself (not by an added edge datum `Λ`/`E₂`, which would reinstate a controller and land on Schema E/C).
> **(SA-5) `k→0` sign-definite mechanism.** The exclusion in (SA-4) must be provided by a **`k→0` sign-definite structure** — not by a UV/Landau/Savvidy object. (The Savvidy sector provably fails: it is *expansive* for `T`, *wrong-sign* for BS, *high-momentum* for the tachyon; `W6-BE-at-2.5`, `PROOF-SAVVIDY-EXACT-SIGN`, `BANACH-FIXED-POINT §3`.)

Equivalently, in one sentence: **`𝔉`'s bound must be the fixed point of a nonlinear IR-support-sensitive relation that the free field does not satisfy — and the failure of the free field must be intrinsic (a `k→0` sign statement), not imported (an `E₂`/edge anchor).**

### 3.2 The residual degree of freedom

Theorem SA leaves exactly one slot free: **the `k→0` sign-definite self-consistency mechanism** — equivalently, the **content of the exclusion in (SA-5)**. All else is constrained. Precisely:

> **Residual DOF = the `k→0` sign structure.** The prover is free to choose **representations**:
> - **(i) the intrinsic nonlinear transform `R`** — Fredholm determinant (`I7/A3`), spectral shift, bootstrap constraint system, DSE/Ward closure, Riccati/FRG map. *All structural choices provably fail (Theorem SA); the Fredholm zero is the unique floor-sensitive one.*
> - **(ii) the coordinate of exclusion** — curvature `G″`, residue `f₁²`, isolation `E₂−m₁²` (=SSC ratio), or Fredholm-zero location. *These are `H1`-equivalent under the squeeze/SSC (`H1-IS-IT-THE-GAP`; `SUPPORT-NOT-ATOM`; `G2D-ISOLATION-I2`).*
> - **(iii) the term decomposition** making the `k→0` bound controllable — e.g. `G″=S−D` (variance split) or the QV martingale sum-vs-product (`INVENT-VARIANCE-INEQ`; `SELFOPT-5OF5`). *A representational gain (the QV martingale removes the `e^{Θ(N_□)}` backward cost) that exposes the same residual `G″_∞≥λ`.{=H1}*
>
> These are **presentation** freedoms. The **invariant content** they all reach is the **same `k→0` sign-definite object**, with two faces that are one:
> - **exclusion of the gapless fixed point** (Schemas A/B/C/E all bottom out here), and
> - **`a→0` zero-finding / isolation uniformity** (Schema Fredholm; `S₂=o(λ₁)=I2`).
> Both are `H1 = I2 = H*`. **There is no second free slot: the `k→0` sign structure is the residual DOF, and its content is the gap.**

Equivalent one-line statement of the DOF: *the prover may choose any intrinsic nonlinear IR-support-sensitive relation and any coordinate/decomposition; whatever they choose, the only unproven step is the `k→0` sign-definite exclusion of the gapless solution — which no structural relation provides and which is `H1`.*

### 3.3 Flip condition (what would make the facet satisfiable) — `[T]`-criterion

The self-anchoring facet becomes satisfiable iff someone exhibits an intrinsic relation `R` (SA-1..SA-5) with a **`k→0` sign-definite structure that excludes the gapless solution without an external anchor**. By §2.7 that is exactly:
```
   an intrinsic nonlinear transform whose solution set at k→0 has no gapless point —
   i.e. a stable, sign-definite, gauge-invariant k→0 quadratic structure (the "6th species")
   that the free field violates.  None is known; it IS H1 = the transport machine.
```
Equivalently: a **provable `H ≥ λP_⊥` at `k→0`** from non-edge data (Schema D realized), or a **`k→0` tachyon / invariant positive cone / `F′(0)>1`** (Schemas A/B/C realized at the coarse mode) — each absent at β=2.5 by a *computed* obstruction, not by lack of effort.

### 3.4 Honest boundary

- **Not a proof of Clay; not a global impossibility theorem** (`∀class ≠ ∀possible`; `KEY-EXISTS §5`). Theorem SA is `[NG]` for the *structural-relation class* + `[O]` for the survivor; the survivor (Fredholm zero) is real and self-anchoring, it simply stalls at `a→0`.
- **What is genuinely new here (`[V-here]` framing):** the identification that **self-anchoring is the "READ" side of `READ ≠ FORCE`** — the intrinsic relation reads the atom when given the spectrum but cannot *floor* it without excluding the gapless solution, and that exclusion is `H1`. This is the facet-level statement of `SPECIES-EXISTENCE`'s `F3∧F4=∅` and `META-GAPFORCING`'s identity-forcing no-go, specialized to "no controller."
- **The facet's payload:** it converts "we cannot bound without an anchor" (a negative) into a **positive structural profile** — `𝔉` must be *self-anchoring ∧ floor-sensitive ∧ `k→0` sign-definite*, and the residual freedom is the `k→0` sign structure. This tells the next construction **where to aim** (the `k→0` sector, intrinsically), not merely that it is hard.
- **No QED.** Clay remains open; the mass gap is not proved. This document is a PRECISE STRUCTURAL SPEC.

---

## §4 · References

**Internal (absolute paths; the banked lemmas this facet assembles — labels carried):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall = Savvidy concavity / key = `G″(0)`/`H1`).
- `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§1 atom `f₁²`; §2 five bars; §3 death map; §4 transport≡confinement); `D:\Dev\math\docs\new-tools-registry.md` (§A3 Fredholm determinant identity).
- `D:\Dev\math\research\unconditional\SPECIES-EXISTENCE.md` (class `ℭ`, `F1–F4`, `F3∧F4=∅` Theorem 2 `[NG]`, equivalence Theorem 3 `[T]`, vNW `[T-shape]`); `META-CROSSDOMAIN-PROTOTYPE.md` (controller-necessity; BS topology inverted; four-prototype table); `META-GAPFORCING-IDENTITY.md` (identity-forcing no-go §3.1; `READ ≠ FORCE` §3.5; requirement portrait §4); `E2-FLIP-MINORANT.md` (Theorem ENCL, DET, CMP, POS; `Φ_sc`).
- `D:\Dev\math\research\unconditional\SMALLCIRCLE-ADIA.md` (decompactification min–max Lemma §3, backward inequality); `MONOTONE-COUPLING.md` (§2 direction; §4 deformation infimum; free inequality backward); `DEFORM-DECOUPLE-4D.md` (deformation→pure-YM = min–max infimum; unifying obstruction).
- `D:\Dev\math\research\unconditional\BANACH-FIXED-POINT.md` (Schema A contraction; neutral `k→0`; Savvidy expansive; β=0 Schauder invariant cone); `MS-GAP-EQUATION.md` (Schema C gap equation; `κ=0` always; `F′(0)>1` pitchfork; tachyon UV); `PROOF-GAP-DIFFERENTIAL-INEQ.md` (no scalar β-ODE; Grönwall extensive); `SELFOPT-5OF5.md` (QV martingale; residual `G″_∞≥λ`=H1); `A4-SPECTRAL-BOOTSTRAP.md` (Schema E; feasible set contains free field; forcing needs `Λ`=H1).
- `D:\Dev\math\research\unconditional\G2D-ISOLATION-I2.md` (SSC `E₂=ρ_r m₁`; vNW pointwise; uniform `I2`=gap); `E2-REGRESS-CIRCULARITY.md` (§0.4 SSC); `SUPPORT-NOT-ATOM.md`; `H1-IS-IT-THE-GAP.md` (squeeze; `H1<H*`); `PROOF-COARSE-POINCARÉ.md` (`H1`=Poincaré); `SCHUR-CONVEX.md` (free-field boundary `f₁²=0`); `INVENT-VARIANCE-INEQ.md` (`G″=S−D`); `INVENT-NEW-IDENTITY.md` / `LOC-DETERMINANT.md` (Fredholm zero; `S₂=o(λ₁)=I2`; `a→0` zero-finding); `PROOF-IR-H1.md` (§4 fast-mode cannot gap slow-mode); `W6-BE-at-2.5.md` (negatives 90–97% high-momentum; coarse plane-wave +0.40β 0% neg); `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_min=−g|B|<0`); `IR-RESPONSE-TEST.md` (trilemma).
- `D:\Dev\math\research\unconditional\PROOF-P3-UNCONDITIONAL.md` (`f₁²>0` at β<1/48 — the β<1/48 survivor of Schema C); `REDTEAM-SPEC.md` (spec self-consistency; (2)≡(5); non-positivity necessity).

**External `[R]` (standard; not re-verified this session):** Courant–Fischer / Rayleigh–Ritz (Reed–Simon IV, Thm XIII.1); Temple (1928), Weinstein (1932), Kato (1949/1966); Banach (1922), Schauder (1930), Leray–Schauder (1934), Krasnosel'kii (1964); Reid, *Riccati Differential Equations* (1972); Källén (1952)/Lehmann (1954); Osterwalder–Schrader (1973), Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Jaffe, *Phys. Rev.* **158** (1967) 1454; Bernstein–Widder (CM ⟺ `ν≥0`); Savvidy, *Phys. Lett. B* **71** (1977); Shen–Zhu–Zhu, *CMP* **400** (2023) 805 (strong-coupling gap); Weyl/Hoffman–Wielandt; Holley–Stroock, *J. Stat. Phys.* **46** (1987) 1159; Caron-Huot–Pokraka–Zahraee, *JHEP* **01** (2024) 195.

---

*Status: **SELF-ANCHORING facet SPEC.** Structural requirement: `𝔉`'s bound must be the fixed point of a nonlinear, IR-support-sensitive, controller-free relation `R` that *intrinsically* excludes the gapless solution by a `k→0` sign-definite mechanism (SA-1..SA-5). **Theorem C + B** formalize controller-necessity / min–max backward (comparison ⇒ order oracle on the floppier side; controlled anchors more-gapped ⇒ free inequality is an upper bound). **Theorem SA** proves every *structural* self-anchoring relation is defeated by the gapless-solution obstruction, and the unique floor-sensitive survivor (Fredholm zero) is `READ`, not `FORCE`, stalling at `a→0` zero-finding. **Residual DOF = the `k→0` sign-definite self-consistency** (representational choices free; invariant content = `H1=I2=H*`). NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay 0%; no QED.*

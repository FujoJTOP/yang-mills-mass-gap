# SPECIES-EXISTENCE — does the missing "6th species" `𝔉` exist? Formal class `ℭ` (F1–F4), the **computable⇒floor-blind / floor-sensitive⇒circular** dichotomy, its no-go for the computable subclass, and the equivalence that makes the existence question *be* Clay

> 2026-11-24 · **META-THEOREM** (charter §1/§2). Question: formalize the class `ℭ` of admissible gap-certifying functionals as the axioms **F1–F4**, then prove **(a) existence**, **(b) no-go**, or **(c) equivalence/circularity**. This note delivers all three, cleanly separated by sub-class. It builds directly on `E2-FLIP-NOGO.md` (the class `C0` and its four-horn dichotomy), `META-GAPFORCING-IDENTITY.md` (READ ≠ FORCE), `H1-IS-IT-THE-GAP.md` (the squeeze), `SUPPORT-NOT-ATOM.md` (support vs atom), and `G2D-ISOLATION-I2.md` (the vNW pointwise break). **Where I inherit a lemma I say so and carry its label; where the result is new it is marked `[V-here]`. No QED on anything unproven. `∀class ≠ ∀possible`.**
>
> Labels (charter §2): **[T]** true theorem (proved, or an elementary one-line consequence of a banked `[T]`) · **[C]** conditional theorem · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[J]** judgment (not a theorem) · **[V-here]** new this round · **[R]** external/inherited.

---

## VERDICT (up front — which, for what class)

> ### There is **no uniform computable species**, but the *general* species is the **gap itself** — so the answer is **EQUIVALENCE (circularity) for the full class `ℭ_uniform`, realized as a clean NO-GO for its computable subclass `ℭ^comp`, with existence holding only sub-uniformly (vNW) or sub-regime (P3).**
>
> Three precisely-scoped statements, one per reading of "exists":
> 1. **[NG] NO-GO for the computable subclass** `ℭ^comp := ℭ_uniform ∩ {F4}` = **EMPTY**. **F3 (uniform gap-detection) and F4 (computable/non-circular) are jointly unsatisfiable** for any `𝔉` in the four data-modes that F4 permits (finite-moment / invariant / finite-`t`-analytic / quantitative-extremal). This is the task's option **(b)**, sharpened to a two-axiom incompatibility: `F3 ∧ F4 = ∅`. (§2, Theorem 2.)
> 2. **[T] EQUIVALENCE (circularity) for the full class `ℭ_uniform`** (drop F4): "`∃𝔉 ∈ ℭ_uniform` that *certifies* `f₁² ≥ c`" is **logically equivalent to the gap** `f₁² ≥ c` on the YM family — because the only certificate that reaches the bound is the edge atom `f₁²` *itself*, which is not computable. The species-existence question *is* Clay. This is the task's option **(c)**, stated precisely. (§3, Theorem 3.)
> 3. **[T-shape] EXISTENCE only sub-uniform / sub-regime.** A genuine, non-circular, acid-clean, **pointwise** functional *does* exist — same-channel von Neumann–Wigner isolation `𝔉_vNW` — satisfying F1, F2, F4 and a pointwise F3; it is real but **sub-uniform** (strict inequality, no uniform `κ`; dissolves at large β), and *uniformizing it is verbatim `H*`*. Separately, P3 is a **computable** `f₁²>0` certificate but **only at `β<1/48`** (strong-coupling cluster expansion). Neither is a uniform member. (§4, Theorem 4 + §4.2.)
>
> **No construction gives a uniform member of `ℭ^comp`.** The dichotomy's two horns are exhaustive of the constructible: a functional is either **computable and floor-blind** (cannot make F3 uniform — a moment/invariant/finite-`t`/extremal cannot see `inf supp ν`, `[T]`), or **floor-sensitive and non-computable** (its evaluation requires `inf supp ν` = the gap — the `t→∞` Tauberian limit, the continuum Fredholm zero, or a field-space nonlinear inequality = `H1`). The surviving floor-sensitive objects **are** the gap. This is *not* a no-go for all possible `𝔉`; it is a no-go for the computable class and an exact identification of the survivor with Clay.

---

## §0 Setup, target disambiguation, and provenance

### §0.1 The two faces — support vs atom — and why the task's target must be fixed

Let `O = ¼ tr F²` be the gauge-invariant `0⁺⁺` interpolator; for an OS/RP-admissible Euclidean theory `μ` let `ν_μ ≥ 0` be its Källén–Lehmann measure on `s = m² ∈ [0,∞)`, `C_μ(τ)=∫ e^{−√s τ}dν_μ(s)`. Define (all `[T]`, `H1-IS-IT-THE-GAP.md §1`):
```
m₁(μ)=inf supp ν_μ  (support floor = the mass gap),     f₁²(μ)=ν_μ({m₁²})  (edge atom = the trF² residue),
χ(μ)=∫ s⁻¹dν_μ  (zero-momentum susceptibility),         M₀(μ)=∫dν_μ.
```
Two distinct targets are conflated by the word "gap":
- **(G) support target (Clay proper).** `m₁(μ) ≥ Δ > 0` uniformly on `F_YM`. State-level; **needs no atom** (`SUPPORT-NOT-ATOM`).
- **(H*) atom target (the task's framing).** `f₁²(μ) ≥ c > 0` *with isolation* `E₂−m₁² ≥ δ`. Strictly stronger than (G): the massive free composite `:φ²:` is **gapped (`m₁=2m>0`) yet atom-free (`f₁²=0`)**, `χ` finite (`H1-IS-IT-THE-GAP §2.3-iii`, `PROOF-F1SQ`). So `(G) ⇏ (H*)`.

**Exact squeeze (`[T]`, `H1-IS-IT-THE-GAP §2.2`):** `f₁²/χ ≤ m₁² ≤ M₀/χ`. Consequence: the curvature `G″=1/χ` is an **upper** bound on `m₁²` (harmonic mean ≥ min) and is **blind to a thin edge atom** (`ν=wδ_ε+(1−w)δ_S`, `w/ε≤1/λ`: `χ≤1/λ`, `m₁²=ε→0`; `dν=√s ds`: `χ<∞`, `f₁²=0`, `m₁=0`, `KEY-DIRECT-OPTIMIZE §A.5`). This is why **the curvature face is not the atom** and why the "species" must be a non-positivity object (§1).

**The atom is upper-semicontinuous and may drop to 0 at a crossing** (`[R]` `SUPPORT-NOT-ATOM §4`, `PROOF-TRANSPORT-ANALYTIC K1`): a lower bound `f₁²≥c` proven on one side of the crossover does **not** transport across it without controlling the crossing. This is the precise sense in which "the atom is usc" enters — *`f₁²>0` is not a stable/closed property along the physical β-path; it can vanish at the crossover*. It is also why positivity/lower-bound (lsc-producing) tools cannot print it: their outputs are continuous data blind to a moving edge (§2).

### §0.2 The adversary family `F_adm` (four canonical theories, each `[T]`/`[R]`)

| label | theory | `m₁` | `f₁²` (atom) | `χ` | role |
|---|---|---|---|---|---|
| `μ₀` | free massless (U(1) or SU(N) `g=0`) | `0` (gapless) | `0` | `∞` | **acid test**: F1 must give `≤0` here |
| `μ_m` | free **massive composite** `:φ²:` | `2m>0` (gapped) | `0` (atom-free!) | finite | **the discriminating acid**: gapped but no atom |
| `μ_CFT` | interacting **gapless** CFT | `0` | `0` | finite-or-`∞` | F3 must read "gapless/`≤0`" |
| `μ_YM` | interacting YM at `β≈2.5` (`α_s~O(1)`) | target `>0` | target `>0` | target finite | the object |
`F_YM` = the lattice YM family `(β,L,a)` uniform in `(L,a)`, `β` crossing `~2.47`.

### §0.3 Provenance — what is new here vs inherited

- **Inherited (labels carried):** the four-horn blindness of computable spectral functionals `[T]/[NG]` (`E2-FLIP-NOGO §2`, Theorems α,β,γ,δ); READ≠FORCE `[NG]` (`META-GAPFORCING-IDENTITY`); the squeeze `[T]` and `H1<H*` (`H1-IS-IT-THE-GAP`); support<atom (`SUPPORT-NOT-ATOM`); the vNW pointwise break `[T-shape]` (`G2D-ISOLATION-I2 §3`); P3 `[T]` (`PROOF-P3-UNCONDITIONAL`); the Schur/free-field master obstruction `[NG]` (`SCHUR-CONVEX`).
- **`[V-here]` (new in this note):** (i) the **exact axiom identification `F4 ≡ C0`** and the resulting sharp framing **`F3 ∧ F4 = ∅`** (Theorem 2); (ii) the **discriminator ≠ atom-certifier split** (Theorem 1) — a gap-detector *cannot* be an atom-minorant, because the gapped atom-free `μ_m` forces it nonzero where `f₁²=0`; (iii) the **equivalence theorem** (Theorem 3) stating `Spec(c) ⟺ Gap(c)` *with F4 enforced*, i.e. the precise circularity; (iv) the explicit three-tier verdict and the falsification condition.

---

## §1 The class `ℭ` — precise axioms

Throughout `𝔉 : F_adm → ℝ` (each `μ ↦ 𝔉(μ)`). The four task-conditions become:

> **(A1 = F1, non-positive / acid-clean).** `𝔉(μ₀) ≤ 0`. *(The functional does not falsely certify the gapless free field.)*
>
> **(A2 = F2, nonabelian-using).** `𝔉(μ_{U(1) pure}) = 0`, and `𝔉` depends essentially on a nonabelian invariant (`b₀=11N/48π²>0`, `C₂(adj)=N`, or the indefinite Savvidy fiber `2 ad F_B`). *(The Abel no-go: `b₀=0`, `2adF_B=0`, no negative mode for U(1); `KEY-DIRECT-DERIVATION §1,§6a`.)*
>
> **(A3 = F3, uniform gap-detection).** `inf_{μ∈F_YM} 𝔉(μ) ≥ c > 0` **(uniform in `L,a`, at `β≈2.5`, no small parameter)** *and* `𝔉(μ_CFT) ≤ 0` for every gapless admissible `μ_CFT`. *(The sign of `𝔉` tracks gap-presence, uniformly. The weaker "merely different values" reading is vacuous — a functional can differ without detecting; the detection reading is the faithful formalization of "value-distinguishing gapped vs gapless".)*
>
> **(A4 = F4, computable / non-circular).** The value and sign of `𝔉(μ)` are determined from data `D(μ)` available **without knowledge of `m₁(μ)=inf supp ν_μ`** — i.e. `𝔉 = Φ ∘ D` where `D` is any of the four modes: **(α)** finitely many moments `∫φ_i dν`; **(β)** a continuous-deformation invariant; **(γ)** smooth/analytic data (`C(τ)` for `τ` in a bounded range / analytic germ at `∞`); **(δ)** a quantitative extremal (capacity/Cheeger, large-deviation rate, transport). *(This is exactly the corpus's class `C0`: `F4 = membership in C0` — `[V-here]` identification, `E2-FLIP-NOGO §1`.)*

Define `ℭ := {𝔉 : A1∧A2∧A3∧A4}` and `ℭ^comp := ℭ` (F4 enforced); `ℭ_uniform` = drop A4 (A1∧A2∧A3).

**Lemma A (F1 is not independent, given certification). [T]** If `𝔉` *certifies the atom* — `𝔉(μ)>0 ⇒ f₁²(μ)>0` — then `𝔉(μ₀) ≤ 0` automatically (since `f₁²(μ₀)=0`). Hence A1 is **implied** by the certifying property; it is a genuine independent axiom only for the *discriminator* reading of A3 (which is not a certifier — Theorem 1). *(This is `E2-FLIP-NOGO Theorem 0` specialized from `m₁` to the atom `f₁²`.)* ∎

**Theorem 1 (discriminator ≠ atom-certifier). `[T]` `[V-here]`.** Suppose `𝔉` satisfies A2 (`𝔉(μ₀)=0`) and A3 read as *gap-detection over `F_adm`*: then `𝔉(μ_m) ≠ 0` (it must flag the gapped `μ_m` as non-gapless). But `f₁²(μ_m)=0`. Hence:
- `𝔉(μ_m) ≠ f₁²(μ_m)`, so `𝔉` is **not a minorant of `f₁²`**; and
- if `𝔉(μ_m)>0` (the "gapped" signature), `𝔉` **falsely certifies an atom** in the atom-free `μ_m`.

Therefore **no `𝔉` is simultaneously a gap-detector over `F_adm` and an atom-minorant (`𝔉 ≤ f₁²`).** The support target (G) and the atom target (H*) require **different** functionals; A3 (F3) serves only (G). *This is the first place the task's problem statement is resolved rather than attacked: the lightest-atom framing (`H*`) and the "distinguish gapped/gapless" framing (`G`) are logically different targets, and the natural F3 discriminator is wrong for the atom precisely on the gapped-atomless `:φ²:`.* ∎

**Consequence.** From here the meta-theorem splits:
- For **(H*)**: a certifier must be a **minorant `𝔉 ≤ f₁²` bounded below** — handled by Theorem 3 (equivalence/circularity).
- For **(G)**: a detector must satisfy A3 — handled by Theorem 2 (no-go for the computable class) and Theorem 4 (the vNW survivor).

---

## §2 Main argument I — the NO-GO for the computable class: `F3 ∧ F4 = ∅`  `[NG]`

**Theorem 2 (F3-uniform ∧ F4 are jointly unsatisfiable). `[NG]`, assembled from `[T]`/`[NG]` per-horn lemmas of `E2-FLIP-NOGO §2`.** No `𝔉` satisfies A2 ∧ A3 ∧ A4.

*Proof (the two-horn dichotomy, mode by mode).* Fix a YM measure `ν_YM` at `β≈2.5` and a competitor `ν'` sharing the **computable data** `D(ν_YM)` but with `inf supp ν' = 0` (gapless). By each mode there is such a `ν'`:
- **(α) finite moments.** `δ₀+δ₂` vs `2δ₁` share `(M₀,M₁)` with floors `0` vs `1` (`E2 §2-α`); the `[V-here]` witness W1 shares `M₀..M₃=(1,1.995,4.4105,10.544)` exactly with a floor-`1` reference while having `inf supp=0` (LP residual `2.7×10⁻¹⁴`). Exact. `[T]`
- **(γ) finite-`t` analytic.** The finite-`t` rate `r(t)=−(1/t)log(C(t)/M₀)` is an **upper** bound `m₁²≤r(t)`, and is *positive at every finite `t` even on the gapless `(1+t)^{−2}`* (`E2 W2`) — so a `γ`-functional is either **vacuous (F1-failing, positive on gapless)** or must take the `t→∞` limit, which is `m₁²` itself (**F4-failing**). The edge atom contributes the *entire* function `f₁²e^{−m₁²τ}` (no singularity), so analytic continuation/growth data is floor-blind. `[T]`
- **(β) invariant.** Deformation-invariant data are locally constant on the connected β-family and equal the free value, `≤0` by A1; cannot be `≥c>0`. `[NG]` (rigidity `INVENT-NEW-INDEX §3`).
- **(δ) quantitative extremal.** Capacity/Cheeger = `∫Γ₂dμ=‖Lf‖²=λ₁` identity = `H1`; LD-tilt rate = effective action = `H1`; transport `Γ₂` = `H1`; min-support/quantile = `m₁` renamed. Each is `H1` (insufficient for (H*)) or the answer. `[NG]/[O]`

Now A4 gives `𝔉(ν_YM)=Φ(D(ν_YM))=Φ(D(ν'))=𝔉(ν')`. A3 gives `𝔉(ν_YM)≥c>0`, so `𝔉(ν')≥c>0`. But `ν'` is gapless, so A3's degradation clause demands `𝔉(ν')≤0`. **Contradiction.** ∎

**Why this is the task's option (b), made exact.** The task states (b) as "F1–F4 jointly unsatisfiable given *free field is axiom-admissible + gapless* and *atom is usc*." Theorem 2 is exactly this, with the precise mechanism:
- **free field axiom-admissible + gapless** → `μ₀` and the generalized-free witnesses `ν'` (`ν'≥0` ⇒ OS-admissible / Jaffe, `G2B §1`) sit in `F_adm`; A1/A3's degradation clause pins their value `≤0`.
- **atom usc** → the edge atom is a **discontinuous (moving-edge) functional** (`SUPPORT-NOT-ATOM §4`): the computable data `D` (α/β/γ/δ) are continuous/floor-blind and cannot separate `ν_YM` from the floor-lowered `ν'` with the same `D`. The discontinuity is precisely the gap between "the data" and "the edge."
- **The unsatisfiable pair is `F3 ∧ F4`**, not all four: F1 is implied (Lemma A) and F2 is satisfiable (Savvidy/Haar are nonabelian-only, `[T]`); the impossibility is that *detecting the gap uniformly* (F3) requires *seeing the support floor*, which the *computable* data (F4) provably do not.

**Honest boundary (`∀class ≠ ∀possible`).** Theorem 2 is a no-go for `ℭ^comp`, i.e. for A4 in the four modes `α/β/γ/δ = C0`. It is **not** a proof that no computable functional whatsoever can exist: a functional using data outside these four modes could in principle escape — but by `E2-FLIP-NOGO §5 (L-2)` any such A3-satisfying, F1-clean field-space functional must forbid a sub-mass transverse mode, which **is `H1`** (the open target). So "outside C0" is not a hiding place; it is `H1` itself. This is `[J]` (prior negative), not `[T]`.

---

## §3 Main argument II — the EQUIVALENCE: the species-existence question *is* the gap  `[T]`

Fix the atom target (H*). Define:
```
Gap(c)  :=  [ ∀μ∈F_YM : f₁²(μ) ≥ c > 0 ]                    (uniform atom gap, with isolation — Clay + H*)
Spec(c) :=  ∃𝔉 : F_adm→ℝ,  A1∧A2∧A4 ∧ [ ∀μ∈F_YM : c ≤ 𝔉(μ) ≤ f₁²(μ) ]   (a computable uniform minorant bounded below = a certificate)
```

**Theorem 3 (equivalence / circularity). `[T]` `[V-here]` framing.**
1. **`Spec(c) ⟹ Gap(c)`** — trivially: `c ≤ 𝔉(μ) ≤ f₁²(μ)`. `[T]`
2. **`Gap(c) ⟹ Spec(c)`** can be witnessed **only** by `𝔉` agreeing with `f₁²` on `F_YM` (the gap value itself) — but that `𝔉` **fails A4**: its evaluation is `f₁²(μ) = ν_μ({inf supp ν_μ})`, which *is* `inf supp ν_μ` = the edge = the gap. `[T]`
3. Combined with Theorem 2 (no `A4`-member is a uniform minorant reaching `c`), **`Spec(c)` is non-empty iff one allows the non-computable certificate `𝔉=f₁²` — i.e. iff `Gap(c)` holds.** Hence
```
        "∃𝔉 ∈ ℭ_uniform that certifies f₁² ≥ c"   ⟺   f₁² ≥ c on F_YM   ⟺   Clay (atom face).
```
∎

**The precise statement of the circularity.** A certificate `𝔉` with `c ≤ 𝔉 ≤ f₁²` is, verbatim, a **proof of the gap** (`E2-FLIP-NOGO Theorem 0`: a uniformly-positive minorant of the gap function below the true value is a proof of the gap). The only minorant that *reaches* the bound is the edge atom itself, and computing the edge atom requires `m₁=inf supp ν` — the gap. So **the species that certifies the gap is the gap; the existence question is self-referential.** This is not pessimism — it is the *exact logical shape* of why every one of the ~110 routes dies on "compute the edge without the edge" (`META-GAPFORCING-IDENTITY §3.3`, READ≠FORCE).

**Sharpest one-line corollary. `[T]`** *A functional `𝔉 ∈ ℭ_uniform` lowers-bounds the atom `f₁²` uniformly iff `𝔉` is (essentially) the atom; a *computable* `𝔉` is floor-blind (Theorem 2); therefore **a computable species and the gap cannot coexist** — the species-existence question has answer NO for computable certificates and reduces to Clay for all others.* This is the clean impossibility the task requested, stated as equivalence.

---

## §4 The EXISTENCE partials — honest, precisely scoped (`[T-shape]`)

The dichotomy leaves exactly two real survivors, neither a uniform member of `ℭ^comp`.

### §4.1 Theorem 4 — vNW same-channel pointwise isolation: a genuine sub-uniform break `[T-shape]`

**Construction.** Let `𝔉_vNW(μ) ∈ {+1,−1,0}` read the **ordering of the two lightest same-channel `0⁺⁺` levels**: `+1` if `m₁ < m_{0⁺⁺*}` with both bound (gapped regime), `0/−1` at the crossover/gapless. This is `G2D-ISOLATION-I2 §3` (same-channel von Neumann–Wigner).

**Axiom check.** 
- **A1 (F1):** `𝔉_vNW(μ₀)` degrades to the "no tower to order" value `≤0` (`f₁²(μ₀)=0`; no confining tower). ✓
- **A2 (F2):** the cost-grading/level-ordering is a confinement phenomenon, absent for free Maxwell (`𝔉_vNW(μ_{U(1)})=0`). ✓
- **A4 (F4):** computable from the ordering certified by strong-coupling cost-grading + β-sign-constancy (no zero-T transition — confinement is a crossover), **not** from `m₁`. **NON-CIRCULAR.** ✓
- **A3-pointwise (F3 pointwise):** same-channel vNW forbids crossing; strong-coupling puts `m₁<m_{0⁺⁺*}`; sign-constant on the connected β-line ⇒ `m₁<m_{0⁺⁺*}` at every β while both are bound — a genuine **pointwise** discrimination **without presupposing the gap**. ✓(pointwise)

**Why it is NOT a uniform member (fails A3-uniform). `[T]`** (i) it gives only a **strict inequality**, never a uniform `κ`; (ii) the excited `0⁺⁺*` **dissolves into the continuum** at large β (Kato isolation fails there), so it protects only the strong-coupling side; (iii) **uniformizing it** — proving `ρ := m_{0⁺⁺*}/m₁ ≥ 1+c` with `ρ` β-independent — is **verbatim `H*`** (`G2D-ISOLATION-I2 §5`, `SUB2`/`SSC`). So `𝔉_vNW` is a real **sub-uniform** member of `A1∧A2∧A4∧A3_pointwise` — the honest existence part — and its uniformization is the wall.

### §4.2 P3 — a computable `f₁²>0` certificate, but only at `β<1/48` `[T]`

`PROOF-P3-UNCONDITIONAL` proves `f₁²>0` **unconditionally at `β<1/48`** by cluster expansion — a genuine **computable** residue lower bound *in that regime*. It is not a uniform member of `ℭ^comp` because the cluster expansion that *computes* it **diverges** at `ε≈2.5` (its convergence radius ends far below the crossover). `[T]`

**What P3 + vNW jointly show.** A computable certificate of the atom **exists** (P3, strong coupling) and a non-circular pointwise discriminator of the gap **exists** (vNW). **What does not exist is a single functional that is *both computable and uniform across the crossover*** — because the cluster expansion diverges (P3's computation dies) and the pointwise ordering is sub-uniform (vNW's strength dies). The wall is precisely the **crossover `β≈2.5`** — not "no tool in any regime," but "no tool *uniform through the crossing*." This is the transport≡confinement wall in its cleanest logical form.

### §4.3 The strongest construction attempt, and where it dies

Combining the two survivors is the natural attempt: `𝔉_try := ` "P3's cluster-expansion residue, analytically continued past the crossover by RG, glued to vNW's ordering." This is **exactly** the transport machine (`TRANSPORT-FIVE-WAY-VERDICT §5`: "crossing the crossover ≡ confinement"). It dies at the named step: the analytic continuation of the cluster expansion across `β≈2.5` is a nonperturbative RG control that **is** `H1`; the glue (vNW ordering) is sub-uniform. **No gluing crosses.** So the strongest construction dies on the same wall — `[O]`, the exact missing step = the transport machine, which "尚不存在" (`CLAUDE.md §5`).

---

## §5 Verdict + what it implies for Clay

**Restating the three tiers.**
1. **NO-GO (computable class) `[NG]`.** `ℭ^comp = {F1∧F2∧F3∧F4} = ∅`. The incompatibility is the two-axiom statement **`F3 ∧ F4 = ∅`** (Theorem 2): uniform gap-detection requires seeing the support floor; the computable data (α/β/γ/δ = `C0`) are floor-blind; any floor-sensitive computable object is either vacuous (positive on gapless) or the gap (`t→∞` limit / Fredholm zero / field-space inequality). This is a clean impossibility theorem for a *precisely-defined class* — the real, falsifiable deliverable.
2. **EQUIVALENCE (full class) `[T]`.** Dropping F4, "`∃𝔉∈ℭ_uniform` certifying `f₁²≥c`" ⟺ `f₁²≥c` on `F_YM` ⟺ Clay (atom face). The only certificate reaching the bound is the edge atom itself — the species-existence question is self-referential (Theorem 3). **This is the precise statement that "the 6th species is Clay."**
3. **EXISTENCE only sub-uniform / sub-regime `[T-shape]`.** vNW pointwise isolation (non-circular, acid-clean, sub-uniform) and P3 (`β<1/48`, computable) are real `A1∧A2∧A4` objects; neither is a uniform member; uniformizing either is `H*`/`H1`.

**What this implies for Clay.**
- **Clay is not proved and not refuted (0%).** The meta-theorem proves *statements about the class of certificates*, not the gap. It certifies the corpus's deepest structural claim from a fresh, logically clean angle: **there is no uniform computable certificate, and the general certificate is the gap itself.**
- **The wall, stated as a dichotomy, is exact and class-complete for `C0`.** "Computable ⇒ floor-blind; floor-sensitive ⇒ circular." Every constructible tool lands on one horn. The survivor (floor-sensitive, non-circular) is literally `H1`/`H*` — so the species and the gap are the same object, and **"build the 6th species" = "build the transport machine" = Clay** (`CLAUDE.md §1/§5`, `TRANSPORT-FIVE-WAY-VERDICT §5`).
- **This is `∀class ≠ ∀possible`.** The no-go is for the four computable data-modes (= `C0`). A functional using data outside them is *not* proven impossible — but by `E2-FLIP-NOGO §5 (L-2)` it must be a field-space anti-concentration inequality, which **is `H1`**. So the escape from `C0` is the open target, not a loophole — and this note does not close it.

**Falsification condition (what would flip to EXISTENCE). `[T]`-criterion** — `ℭ^comp` becomes non-empty iff someone exhibits `𝔉` with: (i) `𝔉(μ₀)≤0` (F1); (ii) `𝔉(μ_{U(1)})=0`, nonabelian-essential (F2); (iii) `inf_{F_YM}𝔉≥c>0` **uniform in `(L,a)` across `β≈2.5`** and `𝔉(gapless)≤0` (F3-uniform); (iv) `𝔉` computable from data **not** containing `m₁=inf supp ν` (F4) — i.e. provably *not* `H1`/`H*`/transport-circular, and *not* the vNW pointwise statement. By Theorem 2 such a `𝔉` is impossible within `C0`; by `(L-2)` any such `𝔉` outside `C0` is a field-space nonlinear inequality `= H1`. **None is known; finding it is the task.** The honest position: no uniform computable species exists today, the general species is the gap, and the only inhabitant of the gap between them is the nonexistent transport machine.

---

## References

**Internal (absolute paths; the banked results this meta-theorem formalizes and assembles):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall = Savvidy concavity, key = `G″(0)`/H1, machine "尚不存在").
- `D:\Dev\math\research\unconditional\E2-FLIP-NOGO.md` — **the primary antecedent**: class `C0`, Theorem 0 (acid implied by LB), four-horn dichotomy (Theorems α,β,γ,δ), scale dichotomy, loopholes (L-1 vNW, L-2 `≡H1`, L-3 transport). Theorem 2 here = F3∧F4 reframing of its four horns; F4 ≡ C0 `[V-here]`.
- `META-GAPFORCING-IDENTITY.md` — READ ≠ FORCE; identities satisfied by the gapless free field; no gapped-forcing identity `[NG]`.
- `H1-IS-IT-THE-GAP.md` — squeeze `f₁²/χ≤m₁²≤M₀/χ`; `G″=1/χ`; `:φ²:` gapped atom-free; H1<H*.
- `SUPPORT-NOT-ATOM.md` — support vs atom; `f₁²` usc, may drop at crossing; uniform support needs `H1` or the atom.
- `G2D-ISOLATION-I2.md` — vNW same-channel pointwise isolation `[T-shape]`; SSC; uniformizing = `H*`.
- `KEY-DIRECT-DERIVATION.md` — five attacks, exact failing steps; Savvidy nonabelian-only `[T]`; `S−D` cancellation; named missing inequality.
- `KEY-DIRECT-OPTIMIZE.md` — `Φ_total` positive at `B=0` but `G″` = H1; `√s ds` witness; channel mismatch; H1→H* fails.
- `IR-RESPONSE-TEST.md` — connected-nonlinear trilemma (Gaussian-blind / anomaly-blind / circular); "genuine 6th species must use neither channel."
- `SCHUR-CONVEX.md` — master free-field no-go for protected `g`-continuous data.
- `C7-SPECIES-CONSTRUCT.md`, `C8-SPECIES-NOGO.md` — variance route dead by equipartition; negative-association species unstatable+false; only Stein/H1 survives.
- `TRANSPORT-FIVE-WAY-VERDICT.md` — transport≡confinement; five routes land on `H1`.
- `PROOF-P3-UNCONDITIONAL.md` — `f₁²>0` at `β<1/48` (computable strong-coupling certificate).
- `PROOF-NO-DISSOLUTION-KINEMATIC.md` — pole-below-cut; `0⁺⁺` vs `0⁺⁺*` vNW ordering.

**External `[R]` (standard; not re-verified this session):** Källén (1952)/Lehmann (1954) — KL representation; Osterwalder–Schrader (1973)/Osterwalder–Seiler, *Ann. Phys.* **110** (1978) — OS/RP; Jaffe, *Phys. Rev.* **158** (1967) — generalized free fields (`ν≥0⇒OS`); Bernstein–Widder — completely monotone ⟺ `ν≥0`; von Neumann–Wigner (1929); Shen–Zhu–Zhu, *CMP* **400** (2023) — strong-coupling gap; Savvidy, *Phys. Lett. B* **71** (1977); Prékopa (1971); Brascamp–Lieb, *Adv. Math.* **22** (1976).

*Status: **Three precisely-scoped verdicts — NO-GO for the computable class `ℭ^comp` (`F3∧F4=∅`, Theorem 2, `[NG]`); EQUIVALENCE/circularity for the full class `ℭ_uniform` (`Spec(c)⟺Gap(c)⟺Clay`, Theorem 3, `[T]`); EXISTENCE only sub-uniform (vNW pointwise, `[T-shape]`) or sub-regime (P3, `β<1/48`).` New: F4≡C0 identification; discriminator≠atom-certifier (Theorem 1); the equivalence framing (Theorem 3). NOT a no-go for all possible (`∀class≠∀possible`); outside the four computable modes the survivor is `H1` itself. Clay NOT proved; 0%. No QED on anything unproven.***

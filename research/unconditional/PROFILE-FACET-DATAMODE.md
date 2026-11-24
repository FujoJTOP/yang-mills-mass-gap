# PROFILE-FACET-DATAMODE — the DATA-MODE facet of the species profile: what `C0` is, why all four of its modes are floor-blind (exact witnesses), the *necessary structure* of a floor-carrying datum, and the verdict on the one residual degree of freedom — **no spectral data mode `D*` escapes (CLOSED); the only exit leaves the spectral category for the field-space inequality `≡ H1`**

> 2026-11-24 · **CRUX FACET** (charter §1/§2). This is the data-mode facet of the "profile" of the object that would prove the 4D pure SU(N) mass gap. It owns the **one residual degree of freedom** left by the meta-theorem `F3 ∧ F4 = ∅` (`SPECIES-EXISTENCE.md §2`): *the data mode* `D` through which a gap-detector reads the theory. Question: is there a data mode `D*` that **(i) carries the floor `inf supp ν`, (ii) is NOT in `C0`, (iii) is computable non-circularly**? If yes, name it; if no, the facet is closed and the profile has **no free DOF in the spectral category**.
>
> Labels (charter §2): **[T]** true theorem (proved, or one-line from a banked `[T]`) · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[J]** judgment (not a theorem) · **[V-here]** new this round · **[R]** external/standard (as-reported, not re-verified this session).
>
> **`∀class ≠ ∀possible`. This note closes the DATA-MODE facet for the *spectral-data category* (all functionals of `ν`, finite or infinite order). It does NOT prove Clay, and it does NOT close the beyond-`ν` (field-space) exit — which is exactly `H1` and is named as the surviving opening. No QED on anything unproven.**

---

## VERDICT (up front)

> ### **NO SPECTRAL DATA MODE `D*` — the DATA-MODE facet is CLOSED in the spectral category.**
>
> Within the category *"data = a function of the Källén–Lehmann measure `ν`"* — which contains `C0` (the four finite/bounded modes `α/β/γ/δ`) **and** every natural infinite-order / analytic extension (the `t→∞` Tauberian rate, the Fredholm determinant zero `z₁`, a radius of convergence / exponential type, a complex-`ω` singularity, the full moment sequence) — **there is no `D*`.** Every floor-carrying member is **circular**: its certified evaluation *is* the floor (`m₁` renamed, `H*`, or exponential clustering / isolation `I2` / marginal gauge BPHZ). Every non-circular / computable member is **floor-blind**: it is vague-continuous (finite data), reads the *top* of the support not the bottom (radius/type), or is positive where it must be zero (finite-`t` rate). The two horns are **exhaustive of the spectral category**:
>
> ```
>   (spectral datum D)   ⇒   [ floor-blind  (D in C0 or a top/UV-reader) ]   XOR   [ circular (D is the floor) ].
> ```
>
> **This is the sharpening of `F3∧F4=∅` that the facet was chartered to determine.** The meta-theorem proved `F3∧F4=∅` *for `C0`* (the four bounded modes). This note proves the incompatibility is **not an artifact of the four-mode restriction**: lifting `C0` to *all* spectral data — adding exactly the floor-**carrying** infinite-order objects the task lists as candidates (a)–(d) — does **not** produce an escape, because every such object lands on the **circular** horn. The floor is a *support-boundary / real-axis-infinity* datum; it is **not** an analytic / radius / type / zero-location datum of the Laplace transform, and its only computable carriers are the rate/zero/pole objects that *are* `m₁`.
>
> **What survives — and it is NOT a data mode.** The one opening left by the whole profile is the exit **out of the spectral category**: a **field-space, sign-definite nonlinear inequality on the measure `μ` that is not a function of `ν`** — proven (`E2-FLIP-L2-COMPLEMENT`, Lemma L0+L2) to be **content-identical to `H1`** (coarse-mode anti-concentration / `m_eff²(k→0)≥λ`). That is a *structural inequality*, not a "data mode," so it does **not** inhabit this facet. **Consequence for the profile: the "one residual DOF" (the data mode) is not a genuine continuum of choices — it is a rigid binary {blind, circular}, and the only real freedom is the species change `ν → μ` into `H1`.** The data-mode facet contributes **no opening**; the profile's opening is forced entirely into the beyond-`ν` facet `= H1 ≡ Clay`.

---

## §0 Setup, the target, and the precise object the facet governs

**Setup (inherited, `SPECIES-EXISTENCE §0`, `H1-IS-IT-THE-GAP §1`).** `μ` = reflection-positive (OS) Euclidean measure, transfer operator `T=e^{−aH}`, unique vacuum `Ω`. For a gauge-invariant `0⁺⁺` interpolator `O=¼tr F²`, the KL spectral measure is `ν(B)=⟨Ω|O†P_H(B)O|Ω⟩` on `s=m²≥0`, `supp ν ∋ m₁²=:inf supp ν` (the floor = the mass gap), edge atom `f₁²=ν({m₁²})`. The correlator is `C(τ)=∫_0^∞ e^{−√s τ}dν(s)=∫ e^{−sτ}dν_s` (equivalently on `s`). Clay (G) = `m₁≥Δ>0` **uniform in `(L,a)`** — a **support** statement; it needs no atom (`SUPPORT-NOT-ATOM`). The atom target `(H*)` = `f₁²≥c>0` with isolation `E₂−m₁²≥δ` is strictly stronger.

**The floor as a functional.** `inf supp ν = sup{c≥0 : ν([0,c))=0}` — it is determined by the **tail masses** `F(c):=ν([0,c))` for all `c`, i.e. by the support *boundary* of `ν`. This is the load-bearing structural fact of the whole facet:

> **[T] Fact 0 (the floor is a support-boundary datum).** `m₁(ν)` depends only on the *cumulative* `c↦ν([0,c))` in a right-neighbourhood of the edge. It is **not** a value of any specific transform of `ν`; it is the location where the distribution function first leaves `0`.

**Two topologies, the dichotomy's engine.** Let `M₊[0,∞)` carry the **vague** topology (`ν_n→ν` iff `∫φ dν_n→∫φ dν` for all bounded continuous `φ` with compact support).

> **[T] Fact 1 (vague-continuity of bounded data).** For every bounded continuous `φ`, `ν ↦ ∫φ dν` is vague-continuous. Hence any functional of finitely many such integrals — moments `M_k=∫s^k dν` on a UV-cut lattice, generalized moments `∫φ_i dν`, susceptibilities, any bounded spectral statistic — is **vague-continuous**.
>
> **[T] Fact 2 (the floor is only upper-semicontinuous, not continuous — `usc`).** `{ν : inf supp ν ≥ c}` is vague-**closed**. *Proof.* If `ν_n→ν` vaguely and `inf supp ν_n ≥ c` ∀n, then for `φ∈C_c` with `supp φ⊂[0,c)`, `∫φ dν=lim∫φ dν_n=0`; so `ν=0` on `[0,c)`, i.e. `inf supp ν≥c`. ∎ Equivalently `ν↦ inf supp ν` is **upper**-semicontinuous. It has **jump discontinuities**: a floor can drop under a vague limit (thin atom pushed to `0`). **[R→T]** (`PROOF-NEVANLINNA` Mechanism 2; `LEVER-nonconvex-edge §0B`.)

**The logical gap the facet must adjudicate.** A **vague-continuous** functional that is **constant on a data-fiber** cannot distinguish two members of that fiber with different floors (Fact 1 vs Fact 2: continuity vs usc). A **discontinuous** (infinite-order / boundary) functional *could* carry the floor — but only if it can be **evaluated with uniform error control that does not presuppose the floor** (else, by Fact 2, a floor-lowered competitor fits the same data within any tolerance). That "uniform control without presupposing the floor" is the exact object this facet interrogates. §1 shows the four `C0` modes are all on the continuous/floor-blind side; §2 shows the four floor-carrying *infinite-order* candidates are all on the discontinuous/**circular** side; §3 states the necessary structure and the closure.

---

## §1 `C0` and its blindness — the four modes, and exact witnesses

**`C0` := functionals `Φ` whose value is a function of the KL measure `ν` of a gauge-invariant interpolator, accessed through one of four modes** (`E2-FLIP-NOGO §1`; `SPECIES-EXISTENCE §1 A4`; this **is** `F4`, `[V-here]` identification `F4 ≡ membership in C0`):

| mode | `Φ` is a function of … | paradigm |
|---|---|---|
| **(α) finite data** | finitely many `∫φ_i dν` (moments `M_k`, generalized moments, f-sum, susceptibilities) or any finite list of smooth spectral statistics | Stieltjes/Christoffel moment problem |
| **(β) invariant data** | a datum constant on continuous deformations of `μ` (homotopy type, index, signed count, `E_d`-type, derived/TQFT invariant) | Morse / index / spectral flow |
| **(γ) smooth/analytic data** | the full `ν` through a **bounded-range** smooth functional (`C(τ)` for `τ` in a bounded interval, its log-derivatives `r(τ)`, the **analytic germ at `∞`** — but *not* a limit taken) | completely-monotone / Tauberian (finite part) |
| **(δ) quantitative-extremal data** | `ν` through an extremal functional (capacity, Cheeger, min-support/quantile, large-deviation rate, transport/Wasserstein) | potential theory / optimal transport |

**Why each mode is floor-blind — with the exact witness.**

### §1.α Finite data — `[T]`
**Theorem α (finite-data blindness).** For any finite `φ_1,…,φ_k` and any moment target `C=(∫φ_i dν)`, if the fiber `{ν≥0 : ∫φ_i dν=C_i}` is non-singleton, it contains members whose support floor is **arbitrarily far below** any prescribed member. Hence no `Φ=F(∫φ_i dν)` lower-bounds `inf supp ν` uniformly over the fiber.
*Exact witnesses.*
- *(order-1, textbook).* `ν_1=δ₀+δ₂` vs `ν_2=2δ₁`: both have `(M₀,M₁)=(2,2)`, floors `inf supp=0` vs `1`. **[T]**
- *(order-3, exact, `[V-here]` in `w7/e2_flip_witness.py` W1).* A gapped reference on `[1,3]` with moments `M₀..M₃=(1,\,1.995,\,4.4105,\,10.544)` is matched **exactly** (LP residual `2.7×10⁻¹⁴`) by a measure that additionally carries mass `w₀=0.0275` at `s=0`, so `inf supp=0`. **Four exact moments do not pin the edge.** **[T]**
- *(general).* The **linear-data lemma**: for any finite exact data and any `ε>0` there is `ν'` with an atom within `ε` of `0` and the same data up to `O(ε)` — because the fiber is infinite-dimensional (a finite set of bounded-continuous constraints on an infinite-dimensional cone). **[T]**
*Mechanism.* Fact 1 (moments are vague-continuous) + Fact 2 (floor is usc): the moments see only the *weak average*, never the *support boundary* `F(c)` near `0`. The floor is not a weak-continuous functional, so a weak-continuous functional fixes it only through the fiber-free choice.

### §1.β Deformation-invariant — `[NG]`
**Theorem β (rigidity).** A datum constant on continuous deformations is locally constant on the connected `β`-family `(0,∞)`, hence equals its free-theory (`β→∞`) value; by `(LB)` and the acid test the free value is `≤0`, so `Φ≤0` everywhere ⇒ not a uniformly-positive floor bound.
*Witness.* A discrete/homotopical invariant (index, signed count) on a connected parameter family is constant (`INVENT-NEW-INDEX §3`, B3). **[NG]**

### §1.γ Bounded-range smooth / analytic germ — `[T]/[NG]`
**Theorem γ.** For `C(τ)=∫e^{−sτ}dν` (completely monotone ⟺ `ν≥0`, Bernstein–Widder):
1. *(finite-`t` rate is an UPPER bound.)* `r(t):=−(1/t)\log(C(t)/C(0))` satisfies `m₁²≤r(t)≤M₁/M₀` and `r(t)↓m₁²` as `t→∞`. So every **finite-`t`** scalar yields `m₁²≤r(t)` — an upper bound, never a lower bound. **[T]** (`TOOL-COMPLETELY-MONOTONE §2`, direction theorem; `[V-here]` W2.)
2. *(the analytic germ at `∞` is floor-blind — atom is ENTIRE.)* An atom at `m₁²` contributes `f₁²e^{−m₁²t}`, **entire in `t`, no singularity**. The a.c. part `∫_{edge}^∞g(s)e^{−st}ds` is analytic exactly on `{Re t>0}`, fixed by the **UV extent** `s→+∞`, *independent of `m₁²`*. So the **radius of convergence** of `C` about any real `t₀>0` is `t₀` (distance to `Re t=0`), the **same for every floor** — gapped atom, gapped continuum, free massless alike. **[T]** (`TOOL-COMPLETELY-MONOTONE §3`; `MS-ANALYTIC-NOZERO`.)
3. *(acid test fails in the computable range.)* On the massless-type `C(t)=(1+t)^{−2}` (floor `0`), `r(t)` is **strictly positive at every finite `t`** (`r=1.62,1.39,0.92,0.48,0.23,0.09` at `t=0.5,1,3,10,30,100`, W2) — positive exactly where it must be `0`, so a finite-`t` functional is **vacuous** (F1-failing). Only the `t→∞` limit returns `0` — and taking it is candidate (c), §2. **[T]**
*Mechanism.* Bounded-range `C(t)` = finitely many (compact family of) Laplace samples = the time-domain face of the linear-data lemma (§1.α); the atom being entire removes any singularity to read.

### §1.δ Quantitative extremal — `[NG]/[O]`
**Theorem δ.** Each quantitative extremal is either the coarse-mode identity `H1` or `m₁` renamed:
- **Capacity / Wiener / Cheeger:** the `L`-uniform conductance floor **is** `H1` (`∫Γ₂dμ=‖Lf‖²=λ₁` identity; `POTENTIAL-THEORY-CAPACITY`).
- **Large-deviation tilt:** the conditioning that fixes an edge sign has rate function **= the effective action = `H1`** (`DYSON-GAS §4`, `E4-RMT-EDGE §3`).
- **Transport / Wasserstein `Γ₂` (Otto–Bakry–Émery):** `∫Γ₂dμ=‖Lf‖²=λ₁` — the gap by an identity = `H1`.
- **Min-support / quantile of `ν`:** computing it needs the edge itself — it is `LB` only because it **is `m₁` renamed**.
**[NG]/[O]** (`E2-FLIP-NOGO §2-δ`.)

### §1.scale The dimensionful closure (why `C0` is not merely four-positivities) — `[NG]`
A uniform dimensionful `Φ` is forced onto the single scale `Λ_QCD` (pure YM is classically scale-free; the only scales are `Λ_QCD,√σ,g|B|,1/a,1/L`; uniformity kills the regulators, `g|B|` is UV/static ≠ transfer gap, `√σ` is 1-form/cross-sector). Then `Φ=Ψ(Λ)` with `(LB)` giving `Ψ(Λ)≤m₁=c·Λ`, so bounding `m₁` requires the dimensionless ratio `c=m₁/Λ` — **the gap** (`E2-FLIP-NOGO §3`, `E2-FLIP-SECOND-SCALE`). **[NG]**

> **§1 verdict.** `C0`'s four modes are floor-blind for **one** structural reason (`[T]`): each reads `ν` through a **vague-continuous / bounded / bounded-range / upper-direction** functional, while the floor (Fact 0) is a **support-boundary** datum that is only **usc** (Fact 2) — so a `C0` functional is constant on a fiber whose floors are unbounded below, and cannot certify a floor from below. This is `F3∧F4=∅` (`SPECIES-EXISTENCE Thm 2`), re-derived here structurally.

---

## §2 The floor-carrying candidates, and which escape `C0` — precise per-candidate verdicts

The facet's real question is whether *leaving* `C0` — using data that is **not** a bounded/vague-continuous `C0` functional — gives a `D*`. The task names four; each is genuinely floor-**carrying** (unlike `C0`), so each **does** "escape" the trivial `F3∧F4=∅` block on the `C0` horns (it is not floor-blind). The decisive question for each is then **(iii): is it computable non-circularly?** For all four, the answer is **no** — but for four *distinct*, exactly-named reasons.

### Candidate (a) — the smallest zero of the Fredholm determinant `D_O(z)=det(1−zT_O)`, `z₁=e^{am₁}` `[T]→[NG]`

*Carries the floor?* **YES, at fixed `a`.** `D_O(z)=∏_{f_n≠0}(1−z e^{−aE_n})`, zeros `z_n=e^{aE_n}`; `m₁=(1/a)\log z₁`. The determinant is a **perfect probe per box**: the smallest zero *is* the atom position, the residue `N_O(z₁)/D_O′(z₁)` *is* the atom (exact fixed-`a` identity, registry A3 = Birman–Schwinger). So (a) is floor-carrying and **not** a bounded-continuous `C0` functional (it is an infinite-order entire object) — it genuinely "escapes `C0`" *as a datum type*.

*But is the zero location "computable data" outside `C0`, and is it non-circular?* **NO, on three independent counts — the location is `m₁` renamed, its extraction is isolation, and its continuum limit is BPHZ.**
1. **`m₁` renamed (δ-horn in determinant coordinates).** `z₁=e^{am₁}⟺m₁=(1/a)\log z₁` is an *identity*. Reading the zero **is** reading the floor — no lower bound is produced; the "tool" outputs the answer. This is exactly the `δ`-horn "min-support = `m₁` renamed" (`§1.δ`), dressed as an analytic zero.
2. **Extraction is isolation `I2`.** To *locate* the smallest zero you must separate `z₁` from the accumulating `z_{n≥2}`; Rouché on the factored determinant closes **iff `S₂=Σ_{n≥2}λ_n=o(λ₁)`** = strict one-state dominance = the isolation `Λ−m₁²≥κ` with multiplicity control = `I2` (`LOC-DETERMINANT §2–3`). The residue cannot supply this (residues are values; the obstruction is *multiplicities*; a tiny-weight atom is still a zero). So extracting the zero from `D_O` **is** `H*`.
3. **Continuum limit is marginal gauge BPHZ.** `z_n=e^{aE_n}→1` as `a→0`; the isolated smallest zero **coalesces into the KL pole** (the edge of a cut), and `Tr T_O^k=Z(S¹_k×T³)` are closed 4D partition functions whose `a→0` renormalization **is** the marginal 4D gauge BPHZ flow (`Σ_j g_j²∼\log\log(1/a)`) — the wall (`CONTINUUM-FREDHOLM`). The zero is not `(L,a)`-uniform; it is a `z`-coordinate artifact.

*Verdict (a).* **Floor-carrying ✓; escapes the *trivial* `C0` block ✓; computable non-circularly ✗** — it is the floor/atom (`m₁` renamed), its extraction presupposes isolation `I2`, and it does not survive `a→0` without BPHZ. **Not `D*`.** It is the `(δ)+(γ)`-horn combination relocated to the transfer variable `z` — the meta-theorem's "the zero is the gap" (`E2-FLIP-MINORANT` Theorem DET: *the determinant coefficients ARE the moments*; `META-CROSSDOMAIN-PROTOTYPE §1`: A3 = Birman–Schwinger, stalls at `a→0` zero-finding).

### Candidate (b) — an infinite-order / non-perturbative datum (an `n→∞` rate, a radius of convergence, a Tauberian boundary) `[T]/[NG]`

*Decompose by what the infinite-order datum actually reads.*
- **Radius of convergence / exponential type.** PER §1.γ-2, the ROC about `t₀>0` is `t₀` (distance to the `Re t=0` UV line), **independent of `m₁²`** for every floor. Growth/type reads the **top of the support (`sup supp`, the UV extent)**, not the bottom. **Floor-blind** ⇒ not even floor-carrying. `[T]` (`TOOL-COMPLETELY-MONOTONE §3`.)
- **Moment growth `lim M_n^{1/n}`.** `=sup supp` (the **top**), and `m₁²≤M_n/M_{n−1}` (moment-ratio) is an **upper** bound on the floor. Floor-blind to the *bottom* and wrong-signed for a lower bound. `[T]` (`GAP-G2-edge`; `TOOL-COMPLETELY-MONOTONE §3`.)
- **The full moment sequence `{M_k}`** *does* carry the floor (Hamburger/Stieltjes determinacy ⟹ `ν` ⟹ `inf supp`). But (i) every finite truncation is floor-blind (§1.α); (ii) `{M_k}↦inf supp` is **discontinuous (usc)** — no *finite* error control converts moment approximations into a *certified lower bound* on the floor (the floor-lowered competitor has moments `O(ε)`-close, W1/W2 of `E2-FLIP-L2-COMPLEMENT`); (iii) for `trF²` the moments are **UV-divergent** in the continuum, and on the lattice to obtain them uniformly to order `M≳2/(am₁)→∞` is exactly the marginal-BPHZ / transport log-sum (`PADE-ROUCHE-DEEP`; `W3`). **Circular / non-uniform.**
- **The Tauberian boundary** = the `t→∞` rate `r(t)→m₁²` = candidate (c), below.

*Verdict (b).* The infinite-order data **split**: those that are *computable / complex-analytic* (radius, type, moment growth) are **floor-blind (read the top)**; the one that *carries the bottom* (the moment sequence, the Tauberian boundary) is **circular / non-uniform** (usc + BPHZ). **No infinite-order datum threads the needle.** Not `D*`.

### Candidate (c) — the regularized support functional `lim_{t→∞}−(1/t)\log C(t)` (= the floor) `[T] → circular`

*Carries the floor?* **YES, exactly** — the Tauberian identity `−(1/t)\log C(t)→inf supp ν` for `ν≥0` (Bernstein–Widder makes `ν≥0` the exact Tauberian side-condition). `[T]` (`TOOL-COMPLETELY-MONOTONE §1`.) It is **not** a bounded-range `C0` datum — it is a limit, hence genuinely outside `C0`; it is **acid-clean in the limit** (returns `0` on free massless).

*Is its computability circular?* **YES — this is the canonical circularity, and it is exact.** Two facets:
1. **The limit *is* the answer.** `m₁² = inf_{t>0} r(t)`; every finite `t` gives `r(t)≥m₁²` (an upper bound, §1.γ-1). A *certified lower bound* needs the **infimum over all `t`**, i.e. the `t→∞` rate controlled uniformly on the scale `t≳1/m₁`. But `C(t)≤M₀e^{−Δt}` uniformly for large `t` **is exponential clustering = the gap** (`m₁²≥Δ ⟺` such `C`). So evaluating the limit to a certified bound **presupposes** the gap. **[T]** (`TOOL-COMPLETELY-MONOTONE §7 (i)`.)
2. **The competitor obstruction (quantitative).** Finitely many exact `C(t_i)` confine the floor only to `[0, −(1/t)\log v]` — an **upper** endpoint, **no lower endpoint**: a floor-lowered atom `+` tail fits the same finite data (the single-time witness, `TOOL-COMPLETELY-MONOTONE §2`). So error-controlled evaluation cannot certify a floor without already knowing it.

*Verdict (c).* **Floor-carrying ✓, outside `C0` ✓, computable non-circularly ✗ — circular in the exact sense:** the limit is well-defined and *equals* `m₁²` (a theorem), but any **certified** evaluation is exponential clustering `= H*`/gap. This is the meta-theorem's "floor-sensitive ⇒ circular" horn in its purest form. **Not `D*`.**

### Candidate (d) — a complex-`ω` analytic-continuation datum (the lightest complex-`ω` singularity of the retarded `0⁺⁺` correlator) `[T]→[NG]/[J]`

*Carries the floor?* **YES in principle** — the gap is the lightest singularity (pole/branch point) of the retarded Green's function `G(ω)` on the real `ω`-axis in the `0⁺⁺` channel; the pole is at `ω=m₁`.

*Is it computable non-circularly?* **NO — it is (i) the resolvent pole = the atom, and (ii) an ill-posed inverse problem.**
1. **The singularity is the resolvent pole = the atom.** On the Euclidean side the edge atom is **entire** (§1.γ-2): there is *no* singularity to continue; the pole at `m₁` appears **only** after analytic continuation, which requires the full `ν` = `circular` (`MS-ANALYTIC-NOZERO`; `TOOL-COMPLETELY-MONOTONE §3`). The complex-`ω` pole *is* the resolvent singularity — the atom `f₁²` at `m₁` — read in another variable.
2. **Ill-posedness (no certified bound).** Analytic continuation from finite-precision Euclidean data to real `ω` is exponentially ill-conditioned (Backus–Gilbert / maximum-entropy / Nevanlinna); it yields an **estimate** of the lowest peak, not a **certified lower bound**. A certified floor from `G(ω)` needs the continuation error control = the same clustering/BPHZ wall. **[J]/[O]**
3. **Positivity lives in energy, not `ω`; the `ω`-plane zero/singularity is a *numerator* (residue) object, not a positivity object** — the `β`/`ω` sign-protection is exactly what `MS-ANALYTIC-NOZERO` shows is absent (`a(β)` is sign-indefinite off the Perron cone).

*Verdict (d).* **Floor-carrying ✓, outside `C0` as a datum type ✓, computable non-circularly ✗** — the `ω`-singularity **is** the resolvent pole/atom (circular), and as an inversion it is ill-posed (not certified). **Not `D*`.**

---

## §3 The verdict — necessary structure, the closure, and the precise opening

### §3.1 The necessary structure a `D*` would have to have `[T within spectral class]`

From §1–§2 the requirements for a floor-carrying computable datum `D*` are **exact**:

> **Structure Theorem (necessary conditions). `[T]`** Any data mode `D` that (i) *carries* `inf supp ν` and (ii) admits a **certified lower bound** `inf_{F_YM}Φ ≥ c>0` (uniform in `(L,a)`) must satisfy:
> - **(S1) Discontinuity.** `D` cannot be a function of finitely many vague-continuous (bounded / bounded-range) data (Fact 1 vs Fact 2). It must encode a **support-boundary** functional (Fact 0): the cumulative `c↦ν([0,c))` near the edge, or an object **isomorphic** to it — the `t→∞` exponential rate, the smallest Fredholm zero, the resolvent pole.
> - **(S2) Error control.** The evaluation of that boundary-functional must come with a **uniform-in-`(L,a)` error bound** — otherwise the floor-lowered competitor (`ν↦ν'` with the same finite/shadow data, W1/W2 of `E2-FLIP-L2-COMPLEMENT`) fits within tolerance and no lower bound is certified.
> - **(S3) Non-circularity.** That error control must **not presuppose** `inf supp ν` (= the gap).

`(S1)` is *achievable* and is exactly why candidates (a)–(d) are not instantly dead: each is an isomorphic carrier of the support boundary. The facet's verdict is that **`(S2)∧(S3)` is unsatisfiable in the spectral category.**

### §3.2 The No-Go `[NG]` — every spectral `D*` is circular

> **No-Go Theorem (no spectral `D*`). `[NG]` for the class `S` = `C0` ∪ {support-boundary-carrying infinite-order spectral data: `t→∞` rate, Fredholm zero, full moment sequence, complex-`ω` singularity}.** For every `D∈S`, `(S2)∧(S3)` cannot hold: the uniform error control on the support-boundary functional `(S2)` **is**, in each of the three isomorphic coordinates, a statement about the gap itself:
> - in the **time** coordinate, `(S2)` ⟺ `C(t)≤M₀e^{−Δt}` uniformly = **exponential clustering = the gap** (`(c)`, `TOOL-COMPLETELY-MONOTONE §7 (i)`);
> - in the **transfer** coordinate, `(S2)` ⟺ `S₂=o(λ₁)` = **isolation `I2`** (`(a)`, `LOC-DETERMINANT §2–3`);
> - in the **resolvent/continuum** coordinate, `(S2)` ⟺ the **marginal gauge BPHZ** renormalization of `det(1−zT_O)` = the transport wall (`(a)`, `CONTINUUM-FREDHOLM §3–4`; `(b)`, `PADE-ROUCHE-DEEP`).
> Hence `(S2)` presupposes the gap in every coordinate: `(S3)` fails. **`S` contains no non-circular floor-carrying computable `D*`.** ∎

### §3.3 The verdict and the precise opening `[NG]` + `[O]`/`[J]`

1. **The DATA-MODE facet is CLOSED in the spectral category `[NG]`.** The "one residual DOF" left by `F3∧F4=∅` is **not** a free degree of freedom: within the spectral category, the floor is a support-boundary datum and every computable carrier of it (rate / zero / pole / full moments) is **`m₁` renamed and/or needs isolation/BPHZ** — the circular horn. Every non-carrier is **floor-blind** — the `C0` horn. The dichotomy `floor-blind ⊕ circular` is **exhaustive of spectral data**; adding the infinite-order candidates (a)–(d) to `C0` does **not** open a third horn. The facet therefore contributes **no opening** to the profile.
2. **The only opening leaves the spectral category — into field space `[O]`/`[J]`.** The meta-theorem's `L-2` (`E2-FLIP-NOGO §5`) and its proof (`E2-FLIP-L2-COMPLEMENT`, Lemmas L0+L2) establish the content-identification: any functional that *raises* the certified floor must certify "the softest (`k→0`) mode has no spectral weight below `Δ`," which **is** `H1` (coarse-mode anti-concentration, `Var_μ(φ)≤(1/λ)ℰ_μ(φ)`, `m_eff²(k→0)≥λ`). This is **not a data mode of `ν`** — it is a **sign-definite nonlinear inequality on the field-space measure `μ` itself**. So it does not inhabit the data-mode facet; it is the **species change** the charter §1 asks for. Its construction (a proven-sign, YM-specific, free-field-failing field-space inequality whose positivity lands at the gap cell `{typical field × k→0}`) is **the transport machine / the missing sixth species = Clay** (`CLAUDE.md §5`; `TRANSPORT-FIVE-WAY-VERDICT §5`).
3. **The profile's opening, stated at this facet.**
   ```
     PROFILE opening (data-mode facet)  =  ∅              [no spectral D*]
     PROFILE opening (whole profile)    =  {field-space μ-inequality ≡ H1}   [only via species change]
   ```
   The data-mode facet **forces** the opening out of the spectral category: it is the precise sense in which *any* repaired profile must stop asking "which data of `ν`?" and start asking "which sign-definite inequality on `μ`?" — the latter being exactly the open problem.

### §3.4 Why this is the sharpest possible statement, and its honest boundary

- **This is stronger than `F3∧F4=∅` for `C0`, and precisely scoped.** `SPECIES-EXISTENCE Thm 2` proved the incompatibility for the four *bounded* modes. This facet proves it is **robust under lifting to all spectral data** — the natural response ("use an infinite-order / analytic / rate datum") is *not* an escape. That is a genuine strengthening: the incompatibility is a property of the **category** (functions of `ν`), not of the four-mode restriction. **[NG]** for the spectral category.
- **It is `∀class ≠ ∀possible`.** The no-go is for the **spectral-data category `S`**. It is *not* a statement that no proof of the gap exists: the beyond-`ν` field-space inequality (`H1`) is explicitly **not** in `S` and is left **open** — in fact the facet *derives* that it is the only survivor, agreeing with `E2-FLIP-L2-COMPLEMENT`. Proving `S` contains no `D*` does **not** prove Clay and does **not** prove `H1` is unprovable or provable.
- **The acid test is satisfied by the *verdict's objects*, not evaded.** Every floor-carrying candidate (a)–(d) is checked on the free massless field: (c) gives `0` **only in the circular limit** (positive at all finite `t`, W2); (a) is acid-clean *per box* but its extraction needs `I2`; (d) has a cut from `ω=0` on the free field (no pole) but continuing to read it is the unresolvable inversion. None is both acid-clean **and** non-circular — which is exactly the failure.
- **Falsification condition (what would flip this facet to OPEN).** A `D*` exists iff one exhibits a **spectral** datum `D(ν)` that (i) carries `inf supp ν`, (ii) is **not** a bounded/vague-continuous `C0` functional, (iii) yields a certified lower bound uniform in `(L,a)` across `β≈2.5`, **without** exponential clustering / isolation `I2` / marginal BPHZ — i.e. a support-boundary carrier with error control **and** without presupposing the floor. By §3.2 every such carrier's control is one of those three equivalent objects; a genuinely fourth coordination would break the no-go. None is known. By `(L-2)`, any beyond-`S` functional that qualifies is `H1`.

---

## §4 Summary table — the facet in one view

| datum `D` | carries floor? | in `C0`? | computable non-circularly? | horn | verdict |
|---|---|---|---|---|---|
| `M_k=∫s^k dν` (finite) | ✗ (W1: 4 moments, floors 0 vs 1) | ✓ (α) | — (blind) | blind | `[T]` |
| deformation invariant | ✗ (rigid, free value `≤0`) | ✓ (β) | — (blind) | blind | `[NG]` |
| `C(t)`, `r(t)` finite-`t` | ✗ (upper; positive on free massless) | ✓ (γ) | — / vacuous | blind | `[T]` |
| capacity / LD / transport | ✗ (= `H1` / `m₁` renamed) | ✓ (δ) | reduces to `H1` | blind→H1 | `[NG]/[O]` |
| ROC / exp. type / moment growth | ✗ (reads *top* / UV extent) | γ-ext | — (blind) | blind | `[T]` |
| **Fredholm zero `z₁=e^{am₁}`** (a) | ✓ | outside-`C0` type | ✗ (`m₁` renamed; extraction `=I2`; `a→0`=BPHZ) | **circular** | `[NG]` |
| **`t→∞` rate** (c) | ✓ | outside-`C0` type | ✗ (limit = the answer; clustering = gap) | **circular** | `[T]/[NG]` |
| **full moment sequence** (b) | ✓ | outside-`C0` type | ✗ (usc: no finite error control; BPHZ to order `2/(am₁)`) | **circular** | `[NG]` |
| **complex-`ω` singularity** (d) | ✓ | outside-`C0` type | ✗ (pole = resolvent atom; ill-posed inversion) | **circular** | `[NG]/[J]` |
| **field-space `μ`-inequality** (L-2) | ✓ (as joint constraint) | **not a `ν`-datum** | open — **is `H1`** | **species change** | `[O] ≡ H1` |

**One line.** `C0` (four modes) = floor-blind; every floor-carrying **spectral** datum = circular; the only survivor (field-space inequality) is **`H1`**, *not a data mode*. The data-mode facet is **closed**; the profile's opening is forced out of the spectral category into `H1`.

---

## §5 References

**Internal (absolute paths; the banked results this facet assembles and sharpens):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop + acid test; §2 labels; §5 wall = Savvidy concavity, key = `G″(0)`/`H1`, machine "尚不存在").
- `SPECIES-EXISTENCE.md` — **the meta-theorem this facet sharpens**: `F3∧F4=∅` (Thm 2), the equivalence/circularity (Thm 3), `F4≡C0` `[V-here]`, the `∀class≠∀possible` boundary.
- `E2-FLIP-NOGO.md` — **the `C0` definition** (§1: `α/β/γ/δ`), the four-horn dichotomy (§2 Theorems α,β,γ,δ), the scale dichotomy (§3), loopholes L-1/L-2/L-3 (§5).
- `E2-FLIP-L2-COMPLEMENT.md` — **the beyond-`ν` closure** (`L-2`): Lemma L0 (`m₁` is a 2-pt-edge fact; beyond-`ν` location-orthogonal) + Lemma L2 (content of any floor-raiser **is** `H1`); the four logical levels; `[V-here]` W1–W3 witnesses.
- `TOOL-COMPLETELY-MONOTONE.md` — **candidate (c)/(b)**: Tauberian identity; direction theorem `r(t)↓m₁²` (finite-`t` = upper); atom-**entire** / ROC-floor-blind (§3); Wiener–Ikehara circularity; the `t→∞`∨Padé∨`H1` exact stuck step (§7).
- `CONTINUUM-FREDHOLM.md` — **candidate (a)**: `z₁=e^{am₁}` coalesces into the KL pole; renormalizing `D_O` = marginal gauge BPHZ; zero not `(L,a)`-uniform.
- `LOC-DETERMINANT.md` — **candidate (a) extraction**: factored Rouché closes iff `S₂=o(λ₁)` = isolation `I2`; residues hit the wrong variable (values vs multiplicities); Padé–Rouché relocates to transport.
- `MS-ANALYTIC-NOZERO.md` — **candidate (d)**: positivity lives in energy not `β`/`ω`; min/max-modulus & Rouché need zero-freeness/isolation already; touching zero = `(CP)=H*`.
- `H1-IS-IT-THE-GAP.md` — the squeeze `f₁²/χ≤m₁²≤M₀/χ`; `G″=1/χ` susceptibility not pole; `:φ²:` gapped-atomless; `H1` necessary-not-sufficient; `H1≡H*` is a conjecture.
- `SUPPORT-NOT-ATOM.md` — support vs atom; `:φ²:` gapped `f₁²=0`; uniform support needs `H1` or a sharp cone.
- `E2-FLIP-MINORANT.md` (Theorem DET: determinant coefficients ARE the moments), `META-CROSSDOMAIN-PROTOTYPE.md` (§1: A3 = Birman–Schwinger; `a→0` zero-finding = BPHZ wall), `E2-FLIP-SECOND-SCALE.md` (dimensional rigidity), `E2-FLIP-DYNAMICAL.md` (Savvidy-Hessian anti-correlation; gap = conditioned variance = `H1`), `PROOF-NEVANLINNA-POLE.md` (floor/atom usc; germ fiber convex), `PROOF-COARSE-POINCARÉ.md` (`H1` = coarse Poincaré = `m_eff²(k→0)`).
- `research/unconditional/w7/e2_flip_witness.py` (W1 order-3 moment fiber; W2 finite-`t` rate positive on massless), `w7/l2_complement_witness.py` (W1 thin-atom `H1`-blind; W2 edge-blind bounded functionals; W3 2pt/4pt orthogonality).

**External — `[R]` standard / textbook (not re-verified this session):** Bernstein (1928/29), Widder *The Laplace Transform* (1941) — completely monotone ⟺ `ν≥0`; Hardy–Littlewood / Wiener (1932) / Ikehara (1931) — Tauberian/Wiener–Ikehara; Paley–Wiener (support ↔ exponential type); Stieltjes (1894), de Montessus — continued fraction/Padé; Krein–Nudelman — Markov moment problem; Fredholm–Hadamard factorization; Seeley–DeWitt / Kato (1976) VII; Reed–Simon IV §XIII.1; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Jaffe, *Phys. Rev.* **158** (1967) 1454; Shen–Zhu–Zhu, *CMP* **400** (2023) 805 (arXiv:2204.12737).

---

*Status:* ***DATA-MODE facet CLOSED in the spectral category.* `C0` (four modes α/β/γ/δ) is floor-blind `[T]` — exact witnesses: `δ₀+δ₂` vs `2δ₁`; W1 order-3 LP-exact match with floors 0 vs 1; W2 finite-`t` rate positive on gapless. The necessary structure of a floor-carrying datum is `(S1)` support-boundary carrier `∧ (S2)` uniform error control `∧ (S3)` non-circular `[T]`; the No-Go `[NG]` is `(S2)∧(S3)` unsatisfiable for every spectral carrier — in time `(S2)=` exponential clustering, in transfer `=I2`, in continuum `=` marginal BPHZ, all `=` the gap. Candidates: (a) Fredholm zero = `m₁` renamed + extraction `=I2` + `a→0` BPHZ; (b) ROC/type/moment-growth = floor-blind (reads top), full moment sequence = usc + BPHZ; (c) `t→∞` rate = the answer / clustering; (d) complex-`ω` = resolvent atom + ill-posed inversion. **NO spectral `D*`.** The one opening is the category exit to the field-space `μ`-inequality `≡ H1` (`E2-FLIP-L2-COMPLEMENT`, Lemma L2) — **not a data mode**, `[O] = Clay`. `∀class ≠ ∀possible`: NOT a proof of Clay; does NOT close the beyond-`ν` exit. Clay stays open, 0%. No QED.*

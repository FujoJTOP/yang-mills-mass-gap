# TOOL-FACTORIZATION-ALGEBRA — can the Costello–Gwilliam factorization-algebra / homotopical (derived) QFT framework, or Atiyah–Segal / cobordism extended TQFT, supply a NON-CIRCULAR, LOWER-bound-shaped criterion for the 4D pure SU(N) mass gap — or does it only restate/presuppose it?

> 2026-11-24 · Scouting a specific area for the Clay gap (`CLAUDE.md` §1: build, don't stop at "it's hard"). Area = **factorization algebras / Costello–Gwilliam homotopical-and-derived QFT**, and **algebraic/topological QFT** (Atiyah–Segal, cobordism, extended TQFT). The steelman: a theory *is* a factorization algebra on ℝ⁴; a **gapped** theory is one whose **IR / observable structure is trivial/short-range** (no propagating IR mode), a **gapless** one flows to a nontrivial CFT. Questions: (1) is "trivial IR / OPE algebra has no nontrivial IR module" a *structural criterion* that could establish the gap, and is it **lower-bound-shaped** or does it **restate** the gap; (2) does the **derived/homotopical** structure (a vanishing of a Hochschild/higher structure) give a **non-circular** cohomological obstruction to gaplessness the analytic routes lack; (3) the 4D CG statements are **perturbative (formal ħ)** — does the framework already **presuppose** the gap (perturbative around a massive vacuum) or can it prove it?
>
> **Labels (`CLAUDE.md` §2 / `E2-REGRESS-CIRCULARITY`):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[R]** as-reported (not re-derived here) / **[V]** citation verified this session / **[V-here]** structural lemma proven-here at the level stated. **No QED on anything unproven; this document does NOT prove Clay and does NOT prove that no homotopical tool can exist** (`∀class ≠ ∀possible`).

---

## Verdict up front

> ### **RESTATES / PRESUPPOSES — no genuinely new (non-analytic) obstruction.** Factorization algebras are a **framework, not a mechanism.** They refine the AQFT/operator-algebraic language homotopically; they do not add a lower-bound engine. Three precise reasons, each mapped onto an already-dead species:

1. **Mass-blindness at fixed scale — a two-horn lemma, no third horn [T-here + maps to rigidity B3 / AQFT].** The free **massive** and free **massless** scalar factorization algebras are built by the **same functor**, are the **same structural kind** (same E_d / P₀ / BV type), and differ **only in the analytic propagator** `P_m=(Δ+m²)^{−1}` vs `P_0`; the family `m ↦ F_m` is continuous (indeed real-analytic) through `m=0`. Therefore any **deformation-invariant** homotopical invariant returns the **same value at `m=0` and at `m=ε>0`** and *cannot* certify a gap. Equivalently: the gap is a **β-non-analytic jump** (`INVENT-NEW-INDEX` rigidity no-go, B3), while homotopical invariants (E_d-type, Hochschild/factorization-homology dimension, derived Euler characteristic) are **locally constant on the connected parameter family** ⇒ frozen at the free value.
2. **The only invariant that jumps is the IR / large-scale limit — and that is an *equivalence*, not a bound.** "The IR factorization algebra is locally constant / trivial ⟺ the theory is gapped" is the homotopical face of "gap ⟺ exponential clustering" — an **identity** (the Tauberian equivalence, `HISTORICAL-TOOL-SWEEP` #12), not a lower bound. And *computing the IR limit nonperturbatively* for pure YM **is** the missing machine (`CLAUDE.md` §5). **Two horns, no third** — verbatim the Family-II shape the project already closed (`E2-REGRESS-CIRCULARITY` §0.6).
3. **Perturbatively it produces the MASSLESS sector, not the gapped one.** Costello–Gwilliam (CG), built on Costello 2011, construct the **perturbative** observables as a **formal power series in ħ** around a classical vacuum. Around the pure-YM vacuum `A=0` the linearized kinetic operator is `−Δ` ⇒ the perturbative factorization algebra is that of a **massless** theory (massless gluons); the gap is nonperturbative (`~e^{−c/g²}`), **outside** the formal-ħ construction. So the framework neither presupposes a gap nor can prove one: the **mass is always an *input*** (a chosen quadratic term) **or absent** (pure YM). [T-here about what CG constructs; [R] on the perturbative/ħ content.]

**TQFT / Atiyah–Segal / cobordism** is a **consistency-and-rigidity** condition (functoriality, gluing, duality; cobordism hypothesis = classification by a fully-dualizable object), **dimensionless/topological**, hence **blind to the mass scale Λ** — the same wrong direction as the **'t Hooft anomaly** (`[NG]`, `PROOF-1FORM-ANOMALY`, `E2 §1 row e`). Pure YM is **not** a TQFT; a gapped YM's IR is the **trivial (vacuum) theory**, so "IR = TQFT" ⟺ "gapped" is again a **restatement**, and the classification is deformation-invariant (B3).

**Derived/"cohomological obstruction to gaplessness"?** No known one, and structurally obstructed in the deformation-invariant class (§2). The would-be invariants (Hochschild/factorization-homology trace, an "IR defect" `δ_ℓ`) are **either** signed/Euler-characteristic counts ⇒ rigidity-frozen at the free value (B3), **or** functions that *contain* `m` as an input (the spectral/partition determinant — a homotopical cousin of `INVENT-NEW-IDENTITY`/A3) **or** decay at rate `=m` (restates the gap). The candidate `δ_ℓ` built and acid-tested in §3 passes the free-field test only **vacuously** and fails as a new mechanism.

**One honest asset (do not oversell):** the framework is the **sharpest available *language*** to state the gapped/gapless dichotomy *exactly* ("IR observable algebra locally constant / trivial" vs "nontrivial CFT") and to **diagnose** that the perturbative CG theory is the gapless sector. That is a *dictionary*, the same role AQFT played in `AQFT-CAUSAL` ("correct language, same wall"), not a lever.

---

## §1 What the framework actually constructs — and where the mass enters

**Factorization algebra (CG Vol 1 §2–3).** A prefactorization algebra on a manifold `M` assigns a cochain complex `F(U)` to each open `U`, with multiplication maps `F(U₁)⊗…⊗F(Uₙ)→F(V)` for pairwise-disjoint `Uᵢ⊂V`, and is *factorization* if it satisfies local-to-global (Weiss/Čech) gluing. **[V]** CG, *Factorization Algebras in Quantum Field Theory*: Vol 1, New Math. Monographs **31**, CUP 2017 (ISBN 9781107163102); Vol 2, New Math. Monographs **41**, CUP 2021 (ISBN 9781107163157).

**The CG main theorem (perturbative).** Any QFT in the sense of Costello 2011 yields a factorization algebra `Obs^q` of quantum observables on `M` (CG Vol 2), built from the **Batalin–Vilkovisky** complex of the effective action via renormalization (heat-kernel/effective-action RG). **[V]** K. Costello, *Renormalization and Effective Field Theory*, AMS Math. Surveys and Monographs **170** (2011) — rigorously **perturbative**, everything in **formal power series in ħ**; no convergence, no nonperturbative completion. **[R]** the CG construction is likewise formal in ħ (their own framing: "any QFT in the sense of [Costello 2011]").

**Where the mass lives.** For a **free** theory the factorization algebra is built from `(E, P, ·)` — the equations of motion `E`, a propagator `P` (a distribution, symmetric, with the correct wavefront set), and the commutative product. The propagator enters the **differential** (BV operator) and the **OPE/multiplication**. For a scalar of mass `m`: `P_m=(−Δ+m²)^{−1}`; for `m=0`: `P_0=(−Δ)^{−1}` (well-defined on ℝ^d, d≥3, decaying like `|x|^{−(d−2)}`, up to the usual constants). **The mass is a parameter of `P` (an input), not an output.** There is no place in the axioms where "gapped" is *produced*; it can only be *read off* from `P`'s decay.

**The "gapped ⇒ IR trivial" folklore, stated honestly.** Massive ⇒ no long-wavelength modes ⇒ large-scale observables frozen (the IR is the vacuum); massless ⇒ long-wavelength modes ⇒ nontrivial IR (CFT). The most recent technical treatment in this exact framework is **[V]** M. Kawahira & T. Shigemura, *States and IR divergences in factorization algebras*, arXiv:2412.08183 (2024; v3 2025): it defines the compactification state, gives a method for **handling IR divergences in the massless theory**, and proves the augmentation, compactification, and Schwartz states **agree in both massive and massless cases**. Two honest readings: (i) the framework **notices** the massive/massless distinction; (ii) it does so **through the analytic structure** (presence/absence of long-wavelength modes; *removal* of IR divergences as an **operation**), and the mass is **taken as input** (the propagator). It claims **no** mass-gap criterion. This is the cleanest evidence for the verdict: the framework *encodes* the gap faithfully and *derives* nothing.

**The derived/∞-structure.** A factorization algebra "locally constant" on ℝ^d is an `E_d`-algebra (Ayala–Francis; **[V]** *Factorization homology of topological manifolds*, J. Topology **8** (2015) 1045–1084, arXiv:1206.5522; locally-constant ⟺ `E_d`, from Lurie's framework). Factorization homology `∫_M F` computes observables on `M`; the `S¹`-compactification is a Hochschild-type trace. **All of these are homotopy/deformation invariants of `F`.**

**Extended TQFT (Atiyah–Segal, cobordism).** An `n`-dim TQFT is a symmetric-monoidal functor `Bord_n→Vect` (or `∞`-categorical target); the **cobordism hypothesis** ([V] Lurie, arXiv:0905.0465, 2009; after Baez–Dolan) classifies fully extended TQFTs by a **fully dualizable object** (the value on a point). **[R]** This is a **classification/rigidity** statement: a TQFT is determined by *finite, dimensionless, topology-only* data.

---

## §2 The load-bearing structural finding — the two-tier blindness lemma ("two horns, no third")

This is the substantive content of the memo. It is what makes the verdict structural rather than "we didn't find it."

**Lemma (two-tier blindness). [V-here / T-here, at the level stated]**
Let `F_m` be the factorization algebra of the free scalar of mass `m` on ℝ^d (d≥3), and let `I` be any invariant of factorization algebras that is **invariant under continuous deformation of `m`** (a homotopy invariant of the `E_d`/`P₀`/BV structure — e.g. quasi-isomorphism type, `E_d`-algebra homotopy type, the dimension of Hochschild/factorization homology, the derived Euler characteristic). Then `I(F_0)=I(F_m)` for all `m`, because `m↦F_m` is continuous through `m=0`. Hence `I` **cannot** be `0` on the gapless (`m=0`) theory and nonzero on the gapped (`m>0`) theory: it is **mass-blind**. The **only** invariant that jumps at `m=0` is the **IR / large-scale limit** `F_∞:=lim_{ℓ→∞} F_m|_{scales ≫ ℓ}` — but `F_∞` is the **vacuum/trivial** algebra iff `m>0` and a **nontrivial CFT** iff `m=0`, i.e. it is a **restatement** of "gapped," and its computation **is** the nonperturbative IR problem. **There is no third horn.**

*Why this is the right sharpening of the project's wall.* It is **verbatim** the project's two established obstructions, in homotopical dress:
- Horn (i) = **rigidity no-go B3** (`INVENT-NEW-INDEX`): "any discrete/homotopical invariant locally constant on the connected β-family is frozen to its free value." The homotopical structure is *by construction* deformation-invariant, so horn (i) bites **harder** here than for indices: the entire derived package is a deformation invariant.
- Horn (ii) = **Tauberian equivalence** (`HISTORICAL-TOOL-SWEEP` #12): "gap ⟺ decay" is an **identity** for `ν≥0` (Bernstein/Laplace), never a lower-bound engine; "IR trivial ⟺ gapped" is the factorization-algebra face of the same identity.
- The fact that it is a **two-horn dichotomy with no third** is the same logical shape as `E2-REGRESS-CIRCULARITY` §0.6 (Temple schema, "two horns, no third"): Family-II inputs are *either* the answer one rung up *or* the answer itself. Here "the answer" is the IR limit.

**Honest boundary of the lemma.** It is a statement about **deformation-invariant** invariants. A *genuinely non-invariant* homotopical object (something that jumps at the gap **without** reading `P`) would evade horn (i). But (a) such an object is not a homotopy invariant of `F` — it would have to be a *scale-sensitive* functional, i.e. built from the analytic kernel; and (b) the only scale-sensitive functionals available are the propagator's decay/determinant, which are the gap (§3). So the escape from horn (i) lands in horn (ii). This is the **exact stuck step** (§5), stated precisely, not hand-waved.

---

## §3 The honest "build it" attempt — the IR-defect functional `δ_ℓ`, and the acid test

The charter (`CLAUDE.md` §1) forbids stopping at "framework, not mechanism" without *trying* to build the missing object. Here is the strongest candidate the framework suggests, defined, acid-tested, and judged.

**Definition.** For the factorization algebra `F` and a scale `ℓ`, let `ι_ℓ : F(B_ℓ) → F(B_{2ℓ})` be the inclusion of concentric balls, and let
```
δ_ℓ(F)  :=  ‖ Proj_{coker H^•(ι_ℓ)} ‖_{OS}        (norm of the failure of H^•(ι_ℓ) to be an iso, in the OS inner product)
```
the "IR defect at scale `ℓ`": how much the observables still **change** when you double the ball. Claim to test: *"`δ_ℓ ≤ C e^{−λℓ}` uniformly ⇒ gap `≥λ`"* and *"does `δ_ℓ` separate free massive from free massless without already containing the answer?"*

**Acid test (free fields), run structurally [V-here]:**
- **Free massless scalar.** Long-wavelength mode survives every scale ⇒ `H^•(ι_ℓ)` is not an iso at any `ℓ` ⇒ `δ_ℓ ↛ 0`. **Correctly reports no gap** ⇒ the tool does **not fabricate** a gap on the free massless field. ✅ (letter of the acid test.)
- **Free massive scalar.** `δ_ℓ ~ e^{−mℓ}` ⇒ the decay rate is exactly `m`. **Reports the gap — but `m` was the input.** So the tool's *output rate is the answer it was asked to find*.
- **Interacting pure YM.** `δ_ℓ` is **defined**, but *bounding it uniformly* (a `λ>0` independent of `ℓ,(L,a)`) **is** "uniform exponential clustering" **is** the gap. `δ_ℓ` is `H1` in factorization-algebra dress.

**Verdict on `δ_ℓ`.** It **restates** the gap (it is the exponential-clustering identity, `HISTORICAL #12`) and is **circular** (the uniform rate `λ` is the gap; the free value is `λ=m` as input). It **passes the free-field acid test only vacuously** — it does not fabricate a gap, but it also derives nothing; it is the homotopical face of the Tauberian identity. **It fails as a new mechanism.** [T-here]

**Two close cousins, both already dead in the repo — do not double-count.**
- **Factorization-homology / Hochschild trace = homotopical A3.** `∫_{S¹}F` (a trace/determinant-like object) *is* a function of `m` and **can see the atom** — exactly the role of the **Fredholm determinant identity** `INVENT-NEW-IDENTITY`/A3 (`f₁² = −(1/z₁)N_O(z₁)/D_O'(z₁)`, "看得见原子" via a **zero**). But A3's residual is the **continuum `a→0`** limit; the homotopical trace's residual is the **IR/scaling limit**. Same residual, homotopical dress.
- **Derived Euler characteristic / signed count = rigidity 0.** The `E_d`-structure's Euler-characteristic / Koszul-sign count is a **signed integer**, deformation-invariant ⇒ frozen at the free value `0` (B3). A signed count of *degree* is the wrong kind of sign — cohomological degree, not the physical positivity that would certify a **lower** bound on the bottom of the spectrum.

So the three natural "derived" candidates fall into the three already-closed bins: **restatement (Tauberian `δ_ℓ`)**, **homotopical determinant (A3)**, **rigidity (B3)**.

---

## §4 Answers to the three posed questions, point by point

**(1) Is "trivial IR / OPE algebra has no nontrivial IR module" a structural criterion, and is it lower-bound-shaped?**
It is a **faithful characterization**, *not* a lower-bound-shaped criterion (§2). Formally "gapped ⟺ IR factorization algebra locally constant/trivial" is **[T]-as-an-equivalence** (it is the homotopical form of gap ⟺ clustering) — but as a *tool* it is horn (ii): a restatement, whose evaluation needs the nonperturbative IR limit. It does **not** give `spec(H)⊆{0}∪[Δ,∞)` with `Δ>0`; it *is* that statement in another alphabet. **[restates]**

**(2) Does the derived/homotopical structure give a NON-CIRCULAR cohomological obstruction to gaplessness?**
**No known non-circular one, and none in the deformation-invariant class** (§2, §3). The candidate invariants are either deformation-invariant (blind, B3), or read the propagator (circular), or are determinant-traces (A3, continuum residual). The one structural reason this is not just "we didn't find it": *Hochschild homology, factorization homology, `E_d`-type, and the cobordism classification are all **deformation invariants** — they are **designed** to be invariant under continuous change of the theory, and the gap is a **non-analytic jump**. A designed-invariant object cannot see a designed-non-invariant feature. **[NG]-for-the-class (deformation-invariant homotopical invariants) + [O] (the exact step below).** This is a *precise negative*, not a barrier theorem (`∀class ≠ ∀possible`).

**(3) Does the framework presuppose the gap (perturbative around a massive vacuum), or can it prove it?**
**Neither presuppose nor prove — it produces the gapless perturbative sector and leaves the gap outside.** Around the pure-YM vacuum `A=0` the CG factorization algebra is that of **massless** gluons (linearized operator `−Δ`); the gap is `~e^{−c/g²}` (dimensional transmutation), **nonperturbative**, not in the formal-ħ CG object. The mass is always an **input parameter** of the quadratic operator (for a free or Higgsed/massive theory) or **absent** (pure YM). So the framework does not smuggle in a gap; it also has no route to one. (If one perturbs around a *Savvidy* background ⟨F⟩≠0, the tree-level "mass" is the **wrong** object — the UV/static-curvature mass `~g|B|`/`~g²⟨F⟩` that `TEST-C-GLOBAL-SCHUR` measured (`~0.6β`, weak-coupling `→0`) and explicitly **excluded** as "not the transfer-matrix gap.") **[T-here]**

---

## §5 The exact stuck step [O] — stated precisely

The one object whose existence would flip this verdict:

> A **homotopical (derived/coherent) invariant `I` of the factorization algebra** such that
> **(i)** `I` is **not** a deformation invariant (it jumps when the theory becomes gapped, without being handed `P` or the IR limit),
> **(ii)** `I` **vanishes iff gapped** — an obstruction to gaplessness,
> **(iii)** `I` is **computable from the UV algebra** (the CG micro-local/UV data) **without** solving the IR/scaling limit, and
> **(iv)** it **fails the free massless field** (returns "no obstruction," non-vacuously).

No such invariant is known. The structural argument of §2 says **(i) and the "homotopical" requirement are in tension**: a derived invariant is *built* to be deformation-invariant, so demanding it jump at the gap forces it to read the analytic kernel — which is **(iii) failing** (the kernel *is* the IR data). This is the **factorization-algebra form of the E₂-regress**: the homotopical invariant is *either* frozen (horn i) *or* it is the IR limit/the propagator (horn ii = the answer). **The regress is factored, not broken.**

**What genuinely new thing this area would need (the honest corner, if any).** The only shape not yet foreclosed is a **non-deformation-invariant derived structure** — something like a *scale-graded* or *filtered* homotopical invariant (a "derived distance to locally-constant") that is **exact/cohomological**, not analytic, yet still detects the decay rate. The attempts in §3 (`δ_ℓ`) show the naive version is just the clustering identity. A real attempt would need a **filtered/weighted Hochschild-type invariant** whose *differential* (not its cohomology) carries the scale — but every such construction we can name either reduces to a determinant (A3) or to the Tauberian identity. **This is [J]-level intuition, not a theorem, and the prior — by the rigidity argument — is negative.**

---

## §6 Charter five-point scorecard (treating the framework as a would-be gap tool)

| # | constraint (`BRIEF-FOR-MODEL` §2) | verdict |
|---|---|---|
| 1 | Nonperturbative | **No** — CG is formal in ħ (Costello 2011); the nonperturbative completion is open and is where the gap lives. |
| 2 | Sign-controlling **lower** bound (`m_eff²(k→0)≥λ` or `f₁²≥r₀`) | **No** — the only jumps are the IR limit (restatement) and the determinant (A3); no lower bound. |
| 3 | Non-circular (no presupposed gap/condensate/`m_g>0`) | **No** — the gap enters as the propagator input or as the IR limit; `δ_ℓ`'s uniform rate = the gap. |
| 4 | Non-SUSY + **fails free massless** | **No (vacuous)** — the framework **admits** the free massless field as a factorization algebra; the massive/massless distinction is *data* (`P`), so any criterion is satisfied at `m=0` and is vacuous, exactly the Family-I pattern (`SCHUR-CONVEX`); AQFT analogue: free massless net is a **counterexample-shaped** bystander (`AQFT-CAUSAL` §3). |
| 5 | `(L,a)`-uniform `≥λ` | **No** — uniform clustering is exactly `H1`/the gap; the IR limit is not `(L,a)`-uniform by construction. |

**Death-map placement (do not invent a 13th mechanism):** this route is **not a new species.** It is the **homotopical refinement of three already-classified rows**:
- **Rigidity / index (Family II, B3 `INVENT-NEW-INDEX`)**: derived invariants are deformation-invariant ⇒ frozen at the free value. *Primary.*
- **Operator algebra / AQFT (`AQFT-CAUSAL`; death-map row "算子代数")**: same mass-blindness — the E_d/factorization structure is the homotopical face of the net/type-III₁ structure, and the free massless theory satisfies it. *Structural sibling.*
- **'t Hooft anomaly / symmetry (`PROOF-1FORM-ANOMALY`, `E2 §1 row e`)**: TQFT/Atiyah–Segal/cobordism is a **consistency/rigidity** condition — *opposite direction* to forcing a gap (anomaly matching *forbids* symmetric trivially-gapped vacua). *TQFT face.*
- Plus the **Tauberian equivalence (`HISTORICAL #12`)** for the IR-limit restatement.

---

## §7 Verdict, flip condition, and what is (honestly) gained

**Verdict: RESTATES / PRESUPPOSES. No genuinely new non-analytic obstruction.** [NG]-for-the-class (deformation-invariant homotopical invariants) + [O] exact stuck step (§5) + [T-here] two-tier blindness lemma (§2) + [J] (no third horn, prior negative). **Not** a proof of Clay; **not** a global impossibility theorem.

**Flip condition.** Produce a **derived/homotopical invariant** satisfying §5(i)–(iv). The rigidity argument (§2) says this must be **non-deformation-invariant** and **UV-computable** at once — a combination no known construction achieves. Exhibit one and this verdict flips to **BROKEN**.

**What is gained (honest, small).**
1. **A precise dictionary.** "gapped ⟺ IR observable algebra locally constant/trivial; gapless ⟺ nontrivial CFT" is now stated *exactly*, and the **two-horn blindness lemma** (§2) shows *why* this dictionary can never be an engine: it is the homotopical form of the E₂-regress / rigidity no-go. This is the same *negative* value AQFT-CAUSAL delivered ("correct language, same wall").
2. **A precise diagnosis of the perturbative CG theory** (§4 Q3): it is the **massless** (gapless) sector; the gap is nonperturbative and outside it. This rules out the hope that "CG already proves the perturbative gap ⇒ done"; it doesn't — the perturbative theory is the *gapless* one, and the nonperturbative gap is exactly the missing machine.
3. **A candidate tool built and honestly killed** (`δ_ℓ`, §3): it passes the acid test only vacuously (does not fabricate a gap) and reduces to the Tauberian identity / `H1` under the uniform-rate requirement. Recorded so no future session re-discovers it as "new."

**Sharpest line.** Factorization algebras / derived QFT are the **cleanest possible statement** that "gapped = the IR observable algebra is trivial" — and that statement is a **translation**, not a **proof**: the massive and massless free fields are the *same kind* of factorization algebra differing only in the analytic propagator, every **deformation-invariant** derived invariant is therefore **mass-blind** (frozen at the free value, the rigidity no-go B3 in homotopical dress), and the **only** invariant that jumps — the IR limit — *is* the gap, whose nonperturbative computation is the machine that does not yet exist. The homotopical structure **factors** the E₂-regress (fixed-scale invariant *or* IR limit; **two horns, no third**) but does not break it. **The wall is relabeled, not crossed; no fabrication, Clay open.**

---

## §8 References

**Internal (absolute paths; the verdicts this memo classifies against):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 the wall and the missing machine); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five constraints; §3 death map); `D:\Dev\math\research\unconditional\HISTORICAL-TOOL-SWEEP.md` (#12 Tauberian equivalence; Family I/II split).
- `E2-REGRESS-CIRCULARITY.md` (§0.6 two horns no third; §5 exact missing steps) — the regress this memo re-derives in homotopical form.
- `AQFT-CAUSAL.md` (mass-blindness of modular/algebraic structure; free massless net as counterexample) — the operator-algebraic sibling.
- `new-tools-registry.md`: A3 `INVENT-NEW-IDENTITY` (Fredholm determinant "sees the atom," continuum residual); B3 `INVENT-NEW-INDEX` (rigidity: discrete invariants frozen at the free value) — the two bins the derived candidates fall into.
- `INVENT-NEW-IDENTITY.md`, `INVENT-NEW-INDEX.md`, `SCHUR-CONVEX.md`, `PROOF-1FORM-ANOMALY.md`, `A-THEOREM.md`, `TEST-C-GLOBAL-SCHUR.md` (§4: static-curvature mass `~0.6β` is *not* the transfer gap), `GAP-G2-edge.md` (moments cannot pin the edge), `W3-adjacent-4D.md` (marginal, must transmute its own scale — why a topological/TQFT structure cannot supply `Λ`).
- `DARBOUX-FACTORIZATION.md` (the *other* "factorization" in the repo — Darboux/isospectral, unrelated; listed to avoid name-collision).

**External — [V] verified this session (bibliographic):**
- K. Costello, *Renormalization and Effective Field Theory*, AMS Math. Surveys and Monographs **170** (2011), ISBN 978-0-8218-5288-0 — perturbative, formal ħ.
- K. Costello & O. Gwilliam, *Factorization Algebras in Quantum Field Theory*, Vol 1 (New Math. Monographs **31**, CUP 2017, ISBN 9781107163102), Vol 2 (New Math. Monographs **41**, CUP 2021, ISBN 9781107163157).
- O. Gwilliam & K. Rejzner, *The observables of a perturbative algebraic quantum field theory form a factorization algebra*, arXiv:2212.08175 [math-ph] (2022, v2 2023) — pAQFT ⇒ factorization algebra (perturbative).
- O. Gwilliam & K. Rejzner, *Relating nets and factorization algebras of observables: free field theories*, arXiv:1711.06674 [math-ph]; *CMP* **373** (2020) 107–174 — pAQFT net ↔ CG factorization algebra for free theories (the mass/propagator is input).
- M. Kawahira & T. Shigemura, *States and IR divergences in factorization algebras*, arXiv:2412.08183 [hep-th] (2024; v3 2025) — IR-divergence handling in the massless theory; augmentation/compactification/Schwartz states agree in **both** massive and massless cases (mass as input; no gap criterion).
- D. Ayala & J. Francis, *Factorization homology of topological manifolds*, *J. Topology* **8** (2015) 1045–1084, arXiv:1206.5522 — `E_d`-algebras ⟺ locally constant factorization algebras; factorization homology.
- J. Lurie, *On the Classification of Topological Field Theories*, arXiv:0905.0465 (2009) — cobordism hypothesis (classification by a fully dualizable object; deformation-invariant).

**External — [R] standard / as-reported (not re-derived here):**
- M. Atiyah, *Topological quantum field theories*, *Publ. Math. IHÉS* **68** (1988) 175–186; G. Segal, *The definition of conformal field theory* (LMS Lecture Note Ser. **308**, 2004); J. Baez & J. Dolan, *Higher-dimensional algebra and topological quantum field theory*, *J. Math. Phys.* **36** (1995) 6073 — the Atiyah–Segal/cobordism axioms.
- "gapped QFT → TQFT at long distances; gapless → CFT" framing: standard folklore (e.g. GGI/lecture notes), *not* a theorem for pure YM.

---

*Status: **[NG] for the class of deformation-invariant homotopical/derived/TQFT invariants** (two-tier blindness lemma, §2: fixed-scale ⇒ mass-blind/rigidity; IR limit ⇒ restatement) + **[O] exact stuck step** (§5: a non-deformation-invariant, UV-computable, free-field-failing derived invariant — none known, structurally obstructed) + **[T-here] diagnostic** (the perturbative CG theory around `A=0` is the *gapless* massless sector; §4 Q3). Verdict: **RESTATES / PRESUPPOSES — framework, not mechanism.** NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay stays open.*

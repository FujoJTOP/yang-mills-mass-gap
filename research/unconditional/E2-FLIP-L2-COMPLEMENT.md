# E2-FLIP-L2-COMPLEMENT — attack the complement of `C0`: is there a field-space (beyond-ν) nonlinear functional that lower-bounds `m₁` and is **not** H1 — or is the complement **= H1**?

> 2026-11-24 · Adversary/constructor pass on **loophole L-2** of `E2-FLIP-NOGO §5`: *"the complement of `C0` = field-space nonlinear inequalities on the gauge field `A` not determined by `ν` alone; the no-go asserts any such acid-clean positive functional IS anti-concentration `= H1`."* This note (a) **tries to build a non-H1 field-space functional**, and (b) **proves the precise sense in which the complement is `H1`** — and pinpoints the logical level at which that is a theorem versus the open Clay core. It is a test of the *escape*, not a re-listing of the `C0` death map.
> **Labels (charter §2):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[V-here]** fresh computation (`research/unconditional/w7/l2_complement_witness.py`, run this session). External facts carried from the repo are **[R]**. **No QED on anything unproven. This does NOT prove Clay; it does not prove no non-H1 functional can exist in an unenumerated class (`∀class ≠ ∀possible`, `KEY-EXISTS §5`).**

---

## Verdict up front

> ### **NO FLIP.** No field-space (beyond-ν) nonlinear functional is constructed that lower-bounds `m₁` and is genuinely non-H1; and constructing one **would be a proof of the gap** (so it is not a free structural win). All *nameable* beyond-ν mechanisms collapse: they either (i) are `H1` in operative content, (ii) fail the free-massless acid test, (iii) are upper-bound/one-sided, (iv) are deformation-rigid, or (v) are **UV/IR-blind** — see the exhaustion table §2.
>
> ### **COMPLEMENT = H1** is *proved at the content level* **[T]** and *demonstrated at the mechanism level* **[NG]-for-all-named-mechanisms + [T/C] anti-localization**, and its *sufficiency / ∀possible* reading is **[O] = the Clay core**. Precisely, there are four logical levels, and honesty requires stating which one "=H1" lives on:
> - **Level I (implication):** *every* functional that certifies `m₁ ≥ Δ` certifies `H1` (the gap ⟹ its own coarse-Poincaré consequence). **[T] — trivial/tautological; rules nothing out**, because *everything* that yields the gap yields `H1`. (This is why the phrase "any acid-clean functional IS H1 by content" is *vacuous* when read at Level I.)
> - **Level II (content-identification):** the *datum* a beyond-ν functional must ultimately deliver to raise the certified `m₁` is exactly the **soft-mode no-low-support (= coarse Poincaré `H1`)** — because `m₁` is a **2-pt spectral-edge** fact and beyond-ν structure is *orthogonal to its location*. **[T] — the real content of "L-2 = H1"**, and it is a genuine (if modest) theorem, not a tautology.
> - **Level III (mechanism):** no field-space functional supplies the gap through an *operative estimate* other than the soft-mode variance bound; the structural reason is **anti-localization of positivity / support mismatch**: the nonlinear/higher-moment/higher-point positivity is largest exactly where the gap is *not* (high `k`, high field, log-cut), and the only datum at the gap cell (typical field × `k→0`) is the conditioned coarse variance `= H1`. **[NG]-for-named-mechanisms** assembled with **[T/C] structural** (conditional on the corpus's measured support facts).
> - **Level IV (∀possible / sufficiency):** "no *conceivable* non-H1 functional, and `H1` suffices" — **[O]**; proving it is (equivalent to) proving the sufficiency `H1 (∧ the provable fiber half H2-global) ⟹ gap` (`PROOF-IR-H1`), i.e. the Clay hard core. **Not claimed.**

**The honest headline.** `L-2` is **not a loophole — it is the wall** (`H1`), but the *reason* is sharper than the prior note's assertion: the **content** (`[T]`, §1) plus the **support-mismatch/anti-localization structure** (`[T/C]`, §3) plus the **candidate exhaustion** (`[NG]`, §2). The strong reading "`H1` is the *whole* complement" is *not* a free standalone theorem — it is, at its strongest defensible form, the banked `[NG]-for-class` + the open `H1 ⟺ gap` `[O]`. So the honest score is unchanged in direction from `E2-FLIP-NOGO` but moves the `L-2` verdict from *asserted* `[O]/[J]` to *partly proved* (`[T]` content + `[NG]` mechanisms), with the irreducible core named exactly.

---

## §0 The objects, stated once, unambiguously

**Setup (`E2-FLIP-NOGO §0–1`, `E2-REGRESS §0.1`).** `μ` = reflection-positive (OS) Euclidean measure with transfer operator `T=e^{-aH}`, unique vacuum `Ω` (PF, `PROOF-GAUSS-INVARIANT-PF`). For a gauge-invariant interpolator `O`, the KL spectral measure is `ν_O(B)=⟨Ω|O†P_H(B)O|Ω⟩` on `s=m²`, `supp ν_O ⊂ [m₁²,∞)`, edge atom `f₁²=ν_O({m₁²})`. Clay G2 = a uniform-in-`(L,a)` **support** bound
```
   m₁ := min over gauge-invariant O of ( inf supp ν_O ) ≥ Δ > 0 .      (the gap; a SUPPORT statement, SUPPORT-NOT-ATOM)
```

**`C0` (the no-go class, restated).** Functionals `Φ` whose value is a function of the KL spectral measure `ν_O` (equivalently the local 2-pt data), accessed by any of four modes: **(α)** finite moments / spectral statistics; **(β)** deformation-invariant (topological/index) data; **(γ)** smooth/analytic/Tauberian data; **(δ)** quantitative-extremal (capacity/Cheeger/large-deviation/transport) data. Four-horn dichotomy `[NG]` (`E2-FLIP-NOGO §2–4`): every `Φ∈C0` is (a) not a lower bound, (b) circular, or (c) `=H1`.

**`L-2` = the complement of `C0`.** Functionals `Ψ(μ)` that depend on the **field-space measure `μ` through structure not captured by any 2-pt spectral measure** `ν_O`. Concretely, "field-space nonlinear, not determined by `ν`" means `Ψ` uses at least one of:
1. **Higher connected correlators** — the connected 3-/4-/…-point functions. (In pure YM at `θ=0`, the 3-point of `tr F²` vanishes by CP; the connected **4-point** is the first nontrivial beyond-2pt datum.)
2. **Nonlinear / non-polynomial functionals** — entropy, relative entropy, exponential moments, large-deviation rates, concentration functions (measure of a set in field space), participation ratios, quasi-entropies.
3. **Sign / correlation inequalities** — FKG / anti-FKG, reflection positivity *beyond* its 2-pt consequences, Gaussian domination.
4. **Topological / geometric measure data** — instanton density, Gribov-region geometry, vortex / Polyakov-loop percolation, Dirac spectral flow.
5. **Global "shape"** — support, multi-modality, **depletion** / phase separation of the order-parameter distribution.

**Two readings of `H1` — must not be conflated (this is the whole subtlety).**
- **H1-coarse-Poincaré** *(the task's object: "anti-concentration/variance bound on the soft mode")*: for block-averaged `φ`, `Var_ν(φ) ≤ (1/λ) ℰ_ν(φ)` with `ℰ_ν(φ)=∫|∇φ|²dν`, uniformly `λ>0`; **⟺ `m_eff²(k→0) ≥ λ`** (`new-tools-registry D.1`, `PROOF-COARSE-POINCARÉ`). By `PROOF-IR-H1 §3` and `OTHER-H1-POINCARÉ`, `H1 ∧ H2-global ⟹ gap`, and `H2-global` is the (conditionally provable) **fiber/UV** half; so **H1 is the IR hard core of the gap**. *This is the relevant closure target.*
- **H1-curvature** *(weaker)*: `G″(k=0)=1/χ=1/∫s^{-1}dν ≥ λ` (`H1-IS-IT-THE-GAP`). This is the **inverse zero-momentum susceptibility**; it is *necessary but not sufficient* (thin-atom blindness, witnessed **W1** below). `=H1` at this reading is **doubly insufficient**, so the only closure that could matter is `=H1-coarse-Poincaré`.

**Flip functional (the target).** `Ψ` with `Ψ(μ) ≤ m₁(μ)` on the Clay family `F`, `inf_F Ψ ≥ Δ>0`, acid-clean (`Ψ(μ_free)=0`, forced by `Ψ≤m₁` and `m₁(free)=0`; `E2-FLIP-NOGO` Theorem 0). By `E2-FLIP-NOGO §0`, **existence of a flip functional ≡ proof of the gap on `F`**; the only non-trivial content is **computability from field-space data without already containing `m₁`**.

---

## §1 Two lemmas: the implication is trivial; the content-identification is the theorem

### Lemma L0 (`gap is a 2-pt-edge fact; beyond-ν is location-orthogonal`) — **[T]**
`m₁` is a function of the collection of 2-pt spectral measures `{ν_O}_O` alone:
```
   m₁ = min_O ( inf supp ν_O ).
```
*Consequences.* **(i)** Two admissible measures `μ,μ'` with the **same 2-pt functions** (hence same all `ν_O`) have the **same `m₁`**, however their higher-point functions differ. **(ii)** A functional depending on `μ` only through structure *orthogonal to the 2-pt sector* (e.g. only the connected 4-pt functions, with the 2-pt sector fixed as a parameter) **cannot distinguish two theories of equal gap** — it is gap-**blind**.
*Proof.* Definition of `m₁` (min over interpolators of the edge of the 2-pt KL measure) + OS reconstruction (the 2-pt functions determine the `ν_O` and vice-versa). No higher-point datum enters `inf supp ν_O`. ∎
*Witness.* **W3** [V-here]: three unit-variance distributions have identical 2-pt (variance `=1`, same gap location in a single-mode toy) but connected 4th cumulants `0, −2, +0.25`; a `κ₄ ≥ c` functional is therefore orthogonal to the edge. (And a gapless interacting fixed point has `κ₄≠0`, so `κ₄ ≥ c` also **fails the acid test**.)

> **Lemma L0 is the precise meaning of "field-space nonlinear, not determined by `ν`":** such `Ψ` is *orthogonal to the gap's location*. It can act only as a **joint constraint / certificate** that (together with the 2-pt sector) forces the edge up — it cannot *be* the edge.

### Lemma L1 (`Level-I implication is a tautology`) — **[T]**
For **every** functional `Ψ` (in `C0` *or* the complement) with `Ψ ≤ m₁` and `inf_F Ψ ≥ Δ`: `Ψ` certifies `m₁ ≥ Δ`, hence the **full** Poincaré constant `= the gap` holds, hence the **coarse** Poincaré `H1` holds (coarse test functions are a subclass).
*Proof.* `gap ⟹ full PI ⟹ coarse PI`: one line, spectral theory of the generator + restricting the test-function class. ∎
> **So "every gap-functional implies H1" is `[T]` but vacuous** — *every* functional, in `C0` or not, implies `H1` once it yields the gap. **The `E2-FLIP-NOGO` sentence "any acid-clean functional IS H1 by content," read at Level I, rules out nothing.** This is the first honest correction: the `L-2` claim has force **only** at Level II/III.

### Lemma L2 (`content-identification: the datum a beyond-ν functional must deliver IS H1`) — **[T]**
To *raise the certified lower bound on `m₁`* is to certify that the **softest (`k→0`) coarse mode has no spectral support below `Δ`**, i.e. `m_eff²(k→0) ≥ Δ`, which **is** `H1`-coarse-Poincaré. Hence any beyond-ν functional that successfully helps certify `m₁ ≥ Δ` must, in its operative step, certify exactly `H1`.
*Proof.* By L0, `m₁` is read off the 2-pt spectral edge, and the minimum is attained at the softest mode `k→0` (`GAP-G2-mass`; `PROOF-COARSE-POINCARÉ`). "No spectral weight below `Δ` for the softest mode" ⟺ `m_eff²(k→0) ≥ Δ` ⟺ `H1` (`new-tools-registry D.1`). A beyond-ν constraint that *raises the certified edge* must force the 2-pt spectral measure to have no low support — that *is* the `H1` datum. ∎
> **This is the genuine "L-2 = H1" `[T]`.** It is **not** a tautology: it asserts that the *target datum* is not merely *implied by* the gap (Level I) but is **identical to** what any beyond-ν certificate must establish. The higher-point *origin* of a constraint does not change its *content*: the content is still "the 2-pt soft mode is massive," which is `H1`.

**Sharpest restatement.** `m₁` lives *only* in the 2-pt sector (L0); beyond-ν structure helps only by constraining the 2-pt sector (L0-ii); a constraint that constrains the 2-pt sector's **edge upward** is an **anti-concentration** statement (L2), because "no support near `0`" is "the measure does not extend down to `0`" is "the soft mode is massive" = `H1`. There is no third thing for a beyond-ν functional to do.

---

## §2 The construction attempt — exhaustive candidate sweep (each tested on the 5 charter conditions)

I genuinely tried to *build* a non-H1 field-space functional (§ charter §1 "build, don't stop"). Each row gives the verdict on the binding conditions: **acid** (fails on free-massless, non-vacuously), **uniform `(L,a)`**, **LB (`≤m₁`)**, **not-H1**. Sources are banked.

| # | beyond-ν candidate `Ψ` | acid-clean? | uniform? | lower-bound? | not-H1? | verdict |
|---|---|---|---|---|---|---|
| 1 | **FSS infrared bound** (RP + Gaussian domination; `⟨φ(x)φ(y)⟩≤` Green's fn) | ✔ (free is Gaussian-exact, mass `0` allowed) | ✘ at `β≈2.5` (domination unproven; gauge not scalar-`φ⁴`) | ✔ (2-pt upper ⇒ mass lower) | **✘** — output **is** an upper bound on the 2-pt `= ` anti-concentration `=H1` | **=H1** |
| 2 | **Cluster / polymer expansion** (polymer-smallness) | ✔ (no scale ⇒ percolation ⇒ no decay) | ✘ at `2.5` (expansion parameter not small; Balaban stops at crossover) | ✔ (in convergence regime) | **✘** — convergence ⇔ massive-Gaussian domination; conclusion (2-pt decay upper ⇔ mass lower) **is H1** | **=H1**, regime-bound |
| 3 | **Higher-point (connected 3-/4-pt)** | ✔ (free: connected 4-pt `=0` for Gaussian) | contested (UV-dominated) | ✘ directly (a decay *rate ≥m₁* is an upper bound on propagation ⇒ needs a clustering **upper** bound `=H1`) | **✘** — by L0-ii it acts only as a joint constraint; the constraint that raises the edge **is H1** | **=H1** |
| 4 | **Entropy / mutual info / area law** | ambiguous | — | **✘** one-directional (gapless area-law states exist, `E2-REGRESS §4.7`); 1-form/line data | (upper-type) | **[NG]** upper-type |
| 5 | **Topological / index / spectral flow** | (free value `0`) | ✘ | ✘ | — | **[NG]** rigid = free value `0` (`INVENT-NEW-INDEX` B3) |
| 6 | **κ₄ / depletion / anti-FKG / signed** | **✘** (gapless interacting has `κ₄≠0`; no FKG in non-abelian gauge ⇒ no sign control; positivity is Schur-**concave** ⇒ upper) | — | ✘ on the sign-carrying cell | — | **[NG]** (fails acid; Schur-upper, `SCHUR-CONVEX`) |
| 7 | **Savvidy / action-Hessian (`−g|B|`)** | ✔ letter (vacuous) | **✘** (`∝V,∝β`) | **✘** anti-correlated (`Hess` depth `↓` while gap `↑`) | — | **[NG]** (`E2-FLIP-DYNAMICAL`, exact 1D witness) |
| 8 | **Spectral geometry / Krein signature** | ✔ (counts `=0` free) | ✘ (counts `∝V`, UV `g|B|`) | ✘ | — | **[NG]** count≠bottom, UV≠IR (`FORGE-KREIN-SPECTRAL`) |
| 9 | **Multi-ν (= all interpolators' 2-pt)** | (still 2-pt) | — | ✔ in principle | **✘** it is still 2-pt data ⇒ in `C0`'s spirit; and (L0) gives nothing beyond `m₁` once the optimal `O` is known | **[NG]** ∈ `C0`-closure |
| 10 | **String tension `σ` (Wilson area coefficient)** | ✔ | (1-form, cross-sector) | **✘** `m₁=c√σ` is *measured* transmutation, not theorem; ratio `m₁/√σ` is `H*` (`BS-THRESHOLD`) | — | **[NG]** ratio=`H*` |

**Reading.** **No candidate both (a) passes acid + uniformity and (b) supplies the lower bound through a non-anti-concentration mechanism.** Every "passes acid" mechanism either fails uniformity/direction (§2.1,2.2,2.3 ⇒ `=H1` anyway) or fails acid outright (§6, `κ₄`/signed). Every "gives a number" mechanism is a count / UV scale / upper bound / rigid. This is the **`[NG]-for-named-mechanisms`** content: ten independent beyond-ν shapes, the same collapse. It is **not** `∀possible`.

**Why the "=H1" column is not a coincidence.** It is forced by L0+L2 (§1): a beyond-ν mechanism that *helps* must raise the 2-pt edge, which is `H1`. The candidate table is the *empirical shadow* of the L2 theorem on the 10 bankable shapes.

---

## §3 The structural theorem: anti-localization = support mismatch — **[T/C]**

Why does *every* beyond-ν positive structure fail to reach the gap? The corpus's deepest structural fact answers this, and it is the mechanism-level proof behind L2.

**Anti-localization of positivity (`TRANSPORT-CORE §4`, `B5-CORRELATION-CONTROL`, sharpened here).** In the exact decomposition of the relevant quadratic/bilinear structure into a **positive (compensation) part** and the **gap-defining part**, every positive term the corpus can *prove* is **largest on the complement of the gap's cell** and **vanishes at** the cell where the gap is decided:
```
   positive term           large where                vanishes at the gap cell (typical field × k→0)
   ‖Hess f‖² (Γ₂)          high momentum k ≫ √(g/β)     k → 0
   Φ″ surplus (Prékopa)    log-cut d→2π  (F≈0)          per-link phase O(π)  (Savvidy on)
   W_geo floor N/12        |F|-free;  survives everywhere  but is only the Haar floor, not the sign
   WSTC (suppresses only)  high field ‖F‖ > κ_typ       typical field ‖F‖ ≈ κ_typ
```
Consequently, at the gap cell **the only surviving datum is the conditioned coarse force variance** `Var_μ(⟨v,∇S⟩|a=0) ⟺ Var(a_{k→0}) ⟺ H1` (`TRANSPORT-CORE` stuck step; `OTHER-H1-POINCARÉ`; `PROOF-JACOBIAN-SAVES`).

**Theorem (support mismatch) — [T/C].** A field-space functional `Ψ` can contribute to a *lower* bound on `m₁` only by supplying **positivity at the gap cell** (typical field × `k→0`), because that is where the `inf` over modes is attained (L0). But every *beyond-ν / nonlinear / higher-moment / higher-point* positive structure is **anti-localized off** the gap cell (the table) — its positivity is largest where the mode is *not* the softest. The only object that lives **on** the gap cell is the conditioned coarse variance `=H1`. Hence a beyond-ν functional offers **no positivity at the gap cell beyond H1**: `L-2 = H1`, structurally.
*Label.* **[T]** for the *mode-attainment* logic (L0) and the *exact* `Γ₂`/Schmidt decomposition identities; **[C]** for the anti-localization *support facts* (conditional on the corpus's measured facts: W4 7.1% negative-mode fraction, W7 `G″=+0.51` carried by conditioning `E[H|0]=+0.22`, `Var(G|0)`; `D/S≈2%`). It is **not** `∀possible`: it is a statement about *the structures SU(N)+Wilson actually possesses*, not a barrier.

> **The synthesis.** `L-2`'s "nonlinear field-space inequality" is exactly the class of **positive structure anti-localized away from the gap cell**. That is *why* it is `H1`: the gap cell contains only `H1`, and the nonlinear structure is, by the RG/`Γ₂` decomposition, *elsewhere*. The "nonlinear axis" of the wall is the **support-mismatch face** of the same wall that the `k`-axis (`TRANSPORT-CORE`) and the action axis (`E2-FLIP-DYNAMICAL`) exhibit. Three faces, one wall.

**Cross-check with `<κ₄`/depletion (§2 row 6).** The alert reader may object: a *depletion* of the soft-mode distribution is a beyond-ν statement that *does* live at `k→0` (it is about the order-parameter value at the softest mode). Does it evade anti-localization? **No:** depletion = a **non-convex effective potential** at the soft mode = the **Savvidy negative direction** (`E2-FLIP-DYNAMICAL`; `OTHER-H1-POINCARÉ`), i.e. the *wall itself*, not a positivity at the gap cell. A depletion of the small-field region is the bimodality whose transfer gap is the **tunnelling** action `~e^{-S_inst}` (`H1`, the missing machine). So depletion is `H1` in yet another costume — consistent with row 6.

---

## §4 The verdict, at each logical level (the precise statement of "complement = H1")

| level | precise claim | label | is it a "strong theorem"? |
|---|---|---|---|
| **I** | any `Ψ` certifying `m₁≥Δ` certifies `H1` | **[T]** | **No — tautological** (everything implies `H1`); rules out nothing. |
| **II** | the **datum** any beyond-ν certificate must deliver to raise `m₁` is the soft-mode no-low-support `=H1` (`m₁` is 2-pt-edge; beyond-ν is location-orthogonal) | **[T]** (L0+L2) | **Yes, modest but genuine** — a content-identification, not a tautology. |
| **III** | no *nameable* beyond-ν mechanism supplies the gap through a non-`H1` operative estimate (anti-localization/support-mismatch) | **[NG]**-for-all-named-mechanisms + **[T/C]** structural | **Yes, for the class** — 10 shapes → one collapse; *not* `∀possible`. |
| **IV** | no *conceivable* non-`H1` functional; and `H1` *suffices* (`H1 ∧ H2-global ⟹ gap`) | **[O]** | **No — this IS the Clay core.** Proving it = proving the gap (via `H1∧H2-global`, `PROOF-IR-H1`; `H2-global` = the provable fiber half). |

**The honest adjudication of the task's two options.**
- **(a) FLIP** — a constructed non-H1 field-space functional — **would be a proof of the gap** (`E2-FLIP-NOGO §0`). None is constructed; all 10 nameable shapes collapse (§2–3). **No flip.**
- **(b) COMPLEMENT = H1** — **proved at Levels II–III** (`[T]` content + `[NG]` mechanisms + `[T/C]` anti-localization), **open at Level IV** (`[O]` = the gap). **"`H1` is the whole complement," read as a standalone *strong* theorem, is NOT a free result:** at Level I it is a tautology, at Level IV it is the Clay core, and only at Levels II–III is it a genuine (achievable) theorem — which is what this note proves.

**This corrects and sharpens `E2-FLIP-NOGO`'s `L-2`.** That note labeled `L-2` "`[O] ≡ H1` — the open target, not a shortcut," as an *assertion*. This note: **(i)** exposes the Level-I vacuity (so the *sentence* "any such `Φ` is `H1` by content" must be read at Levels II–III or it is empty); **(ii)** *proves* the Level-II content-identification `[T]`; **(iii)** *grounds* the Level-III reduction in the anti-localization/support-mismatch structure `[T/C]` rather than assertion; **(iv)** exhausts 10 beyond-ν shapes `[NG]`; **(v)** isolates the residual as exactly `[O] = H1 ⟺ gap`. The direction is unchanged; the *status* is upgraded from `[O]/[J]` to `[T]+[NG]+[T/C]+[O]`.

---

## §5 `[V-here]` numeric witnesses (`l2_complement_witness.py`, single process, numpy)

Three two-minute checks, each an illustration of a banked lemma (not a Clay computation):
```
W1  thin-atom H1-blindness   nu = eps*delta_eps + (1-eps)*delta_S:
    eps      m1=sqrt(eps)     chi        G''=1/chi    H1(G''>=0.05)?
    1e-1     3.16e-1          1.2250     0.8163       True
    1e-2     1.00e-1          1.2475     0.8016       True
    1e-4     1.00e-2          1.2500     0.8000       True
    1e-6     1.00e-3          1.2500     0.8000       True
  -> H1-curvature bounded below (G''->0.8) while m1->0 : H1-curvature is NOT sufficient
     (thin atom).  [necessity gap=>chi<=M0/lam is fine; sufficiency fails without the atom]

W2  nonlinear-but-bounded spectral functionals edge-blind to the thin atom:
    eps      m1        Int e^-s dnu    participation   entropy
    1e-2     1.0e-1    0.0701          1.0202          0.0560
    1e-4     1.0e-2    0.0609          1.0002          0.0010
    1e-6     1.0e-3    0.0608          1.0000          0.00002
  -> every bounded/continuous nonlinear functional moves by O(eps) while m1->0;
     only a SHARP 'mass below delta' functional sees the edge -- and its uniform,
     quantitative version (effective mass m_eff^2(k->0) >= lam) is precisely H1.

W3  2-pt vs 4-pt orthogonality (var=1 all three, kappa4 differs):
    Gaussian N(0,1):        var=0.9997  kappa4= 0.0045
    Rademacher +-1:         var=1.0000  kappa4=-2.0000
    two-scale a=2,b=.5,p=.2: var=1.0013  kappa4= 0.2478
  -> identical 2-pt (gap location) but different 4-pt : beyond-2pt is location-orthogonal.
```
These witness, respectively: `H1`-curvature insufficiency (§0 the two readings), the edge-blindness of bounded nonlinear spectral functionals (W2 = Lemma L2's "only the sharp/H1 datum sees the edge"), and Lemma L0's 2-pt/4-pt orthogonality. They are *illustrations of banked algebra*, **not** SU(N) numbers.

---

## §6 Honest boundary — what this does and does not say

- **Does NOT prove Clay; does NOT prove no non-H1 functional can exist** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`, `CONVERGENCE-redteam §d`). Levels II–III are a `[T]` content-theorem and an `[NG]`-for-mechanisms + `[T/C]` structural theorem; Level IV is `[O]` = the gap. **The complement is not proven *globally* empty**; a not-yet-invented 6th species could live in an unenumerated beyond-ν class — but by L2 (L0), if it certifies the gap it certifies `H1`, and by §3 its positivity must land *on* the gap cell, which is exactly what every nameable structure fails to do.
- **The task's framing "H1 is the complement would be a strong theorem" is itself the honest crux:** it is *strong* exactly to the extent it is *not* a tautology — and the non-tautological part (Levels II–III) is what is proved, while the *strongest* form (Level IV) is the open problem, **not** a free corollary. Claiming a ∀possible standalone theorem would be exactly the `∀class ⇒ ∀possible` error the corpus forbids.
- **The `[T/C]` label on §3 is load-bearing:** the exact decomposition identities (`Γ₂`, Schmidt, orthogonal split) are `[T]`; the *anti-localization support facts* are `[C]`, conditional on the corpus's measured support data (W4, W7, `Var(G|0)`, `D/S≈2%`). If a future measurement relocated a positive term *onto* the gap cell, §3's conclusion would weaken to `[O]` — that is the honest falsifier.
- **The FLIP condition (unchanged in spirit, sharpened).** The verdict flips to **BROKEN** iff someone exhibits a `Ψ` that is **(i)** a uniform `(L,a)` lower bound on `m₁`, **(ii)** acid-clean (no gap on free massless), **(iii)** uses beyond-ν field-space structure, and **(iv)** whose *operative positivity is at the gap cell* (typical field × `k→0`) — i.e. a beyond-ν positive term **not** anti-localized off `{typical field × k→0}`. By §3, that object would *by construction* be the missing 6th species and would **prove the gap**. None is known; producing one is the task, and it is not a relabeling.
- **No fabricated flip, no QED. Clay stays open.**

---

## §7 References

**Internal (absolute paths; the banked lemmas this note assembles/extends):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop + acid-test; §2 labels; §5 wall/machine); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five conditions; §3 death map; §7 "positivity proves disorder/upper, atom needs order/lower"); `D:\Dev\math\docs\new-tools-registry.md` (§D residual dictionary D.1 `H1`; A2 `G″=S−D`; A5 `Γ₂`/W-熵; B3/Krein; Savvidy sign).
- `research/unconditional/E2-FLIP-NOGO.md` (**§1 `C0`; §5 L-2 — the target of this note; Theorem 0; §3 scale dichotomy**); `E2-REGRESS-CIRCULARITY.md` (§0.2 L1; §0.4 SSC; §6 flip condition); `E2-FLIP-DYNAMICAL.md` (**Savvidy anti-correlation `[T-in-model]`; the wall is the conditioned variance**); `KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 missed lever; §7 calibration; `∀class≠∀possible`); `CONVERGENCE-redteam.md` (§b3 non-convex edge; §d).
- `research/unconditional/PROOF-COARSE-POINCARE.md` (H1 = coarse Poincaré = `m_eff²(k→0)`; Lemma SC at `β<1/48`); `OTHER-H1-POINCARÉ.md` (direct coarse-Dirichlet classification: (A) circular / (B) H1); `PROOF-IR-H1.md` (`H1` strictly stronger than `Δ≥m`; `H1∧H2-global⟹gap`; exact orthogonal split; fiber=provable half); `H1-IS-IT-THE-GAP.md` (**two readings; `G″=1/χ` susceptibility; squeeze `f₁²/χ≤m₁²≤M₀/χ`; thin-atom**); `TRANSPORT-CORE.md` (**§4 anti-localization of positivity — the structural heart of §3 here**); `B5-CORRELATION-CONTROL.md` (Elitzur/Schmidt close, wall unmoved).
- `research/unconditional/SCHUR-CONVEX.md` (master Schur no-go; positivity ⇒ upper; free field `1≤0`); `GAP-G2-edge.md` (residue–level lemma); `PROOF-F1SQ.md` / `SUPPORT-NOT-ATOM.md` (support vs atom; massive `:φ²:` `f₁²=0`); `INVENT-NEW-INDEX.md` (B3 rigidity); `FORGE-KREIN-SPECTRAL.md` (count≠bottom, UV≠IR); `BS-THRESHOLD.md` (`√σ` 1-form, cross-sector); `TOOL-COMPLETELY-MONOTONE.md` (§1–3 direction/atom/CM).
- `research/unconditional/w7/l2_complement_witness.py` (**`[V-here]` W1–W3**; run `py -3.12 l2_complement_witness.py`).

**External — `[R]` standard / textbook (not re-verified this session):** Fröhlich–Simon–Spencer, *CMP* **50** (1976) 79 (infrared bound / Gaussian domination); Brydges–Fröhlich–Sokal / cluster expansions (Glimm–Jaffe; Brydges–Dimock–Hurd); Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440 (lattice YM OS positivity); Bakry–Émery / Bakry–Gentil–Ledoux (Poincaré of a diffusion generator); Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737 (`β<1/48` gap); Savvidy, *PLB* **71** (1977) 133; Nielsen–Olesen, *NPB* **144** (1978) 376; Elitzur, *PRD* **12** (1975) 3978. No fabricated bound; no QED.

---

*Status: **NO FLIP.** `L-2` (the complement of `C0`) is not a loophole but the wall: **COMPLEMENT = H1** is proved at the **content level `[T]`** (Lemmas L0+L2: `m₁` is a 2-pt-edge datum; beyond-ν structure is location-orthogonal; any certificate that raises `m₁` certifies the soft-mode no-low-support `=H1`), demonstrated at the **mechanism level `[NG]`-for-10-named-shapes + `[T/C]` anti-localization/support-mismatch** (nonlinear/higher-point positivity is anti-localized off the gap cell typical-field × `k→0`; only the conditioned coarse variance `=H1` lives there), and **open at the sufficiency/∀possible level `[O] = the Clay core`**. The Level-I reading ("any acid-clean functional implies H1") is a `[T]` tautology carrying no information — so "H1 is the whole complement" is a *strong* statement only at Levels II–III (proved here) and is the gap itself at Level IV (not a free theorem). Corrects `E2-FLIP-NOGO`'s asserted `L-2` to a partly-proved result with the residual named exactly. NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay stays open.***

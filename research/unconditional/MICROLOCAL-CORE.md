# MICROLOCAL-CORE — 11th species: analytic microlocal analysis (Hörmander's `WF_A`, Sato–Kawai–Kashiwara, propagation of analytic singularities) on the Wightman functions

> 2026-11-23 · CORE attack on the candidate **11th species**: the *analytic* wavefront set `WF_A`, the analytic singular spectrum `SS` of hyperfunctions, and the **propagation of analytic singularities**. The hope is the one the task names: a Hörmander propagation *lower bound* on the distance from the real axis to the nearest singularity of the two-point function — equivalently the tube width / spacelike-decay rate `e^{−m|x|}` — giving the mass `m` **without any measure-positivity statement**.
>
> **判词 (verdict, up front): NO NEW LEVER — score 2/5, and it is the *analytic principal-symbol twin* of `HORMANDER.md`, landing on the same wall by the same mechanism.** It is a genuinely different *category* from measure positivity (distribution theory / complex analysis, not variance; not modular either — cf. the dead 10th species), and it is the first species whose **free-*mass* tracking is exact and non-vacuous on *both* sides** (free massless → distance `0`; free massive → width `m`; both in the continuum). But it fails three beats for a structural reason: **(a)** the propagation theorem is a **principal-symbol / characteristic-variety** statement — `Char(□_YM) = {ξ²=0}` — which is *identical for `m=0` and `m>0`*, the mass `m²` being a **lower-order (subprincipal)** datum that never enters the bicharacteristic flow; **(b)** the quantity that *does* carry the mass — the tube width at spacelike infinity — is `inf supp ν` by a Paley–Wiener identity `= A4's m_GI = the KL edge = H1`, so a lower bound on it *is* the gap restated (circular); **(c)** worse, it is **atom-blind**: the massive free scalar `O=:φ²:` has tube width `m>0` yet `f₁²=0` (`REDUCTION-AUDIT`, `PROOF-TRANSPORT-ANALYTIC`), so even a magic non-circular *distance* bound certifies the **gap**, not the Clay **atom** (`gap ⇏ atom`). The Euclidean theory is **elliptic** (no real bicharacteristics ⇒ only regularity, the `k→∞` side), the Minkowski theory is **nonlinear hyperbolic with no 4D solution theory** (the `G1` existence wall), and the anti-localization principle survives in *phase space*: singularities propagate on the light cone (`ξ²=0`, high frequency), the mass is read at spacelike infinity (`ξ→0`, off-cone). **Score 2/5** (passes #1 nonperturbative, #4 free-field; fails #2 lower-bound, #3 non-circular, #5 uniform) — identical profile to `HORMANDER.md` (2/5), `AQFT-CAUSAL` (1/5), `TRANSPORT-CORE`/`CHEEGER` (2/5): *right language/object, same wall.*

---

## §0. STEP-0 honesty discharge (what the corpus already banks; the one untested limb)

**Microlocal analysis is partly in the corpus — but not this limb.** The grep (`microlocal|wavefront|Hörmander|Sato|Kawai|Kashiwara|analytic singular|hyperfunction|singular support|analytic wavefront`) returns exactly one live file plus incidental `Kawai` (Eguchi–Kawai) and `Wightman/KL` hits:

- **`HORMANDER.md`** (2/5) — the **subellipticity / hypoellipticity** species: the *real-axis* Hörmander theorem (`{X_a}` bracket-generating ⇒ `‖u‖²_{H^{1/r}} ≲ Σ‖X_a u‖² + ‖u‖²`). Its no-go: the gauge frame is a *global* frame ⇒ bracket condition holds **elliptically** (`r=1`, the free Laplacian's principal part), subellipticity is a `k→∞` **regularity** gain, the gap is the `k→0` **spectral bottom**, and `ΣX_a²` is sign-blind. **This is the `C^∞` / real-principal-symbol sibling of the species I test here.** The analytic `WF_A` limb (complexified singularity location, Sato–Kawai–Kashiwara, propagation along bicharacteristics, Paley–Wiener tube width) is **not** in `HORMANDER.md` and **not** anywhere in the corpus — `Sato|Kawai|Kashiwara` matches only Eguchi–Kawai; `analytic wavefront`/`hyperfunction` match nothing. §1–§5 build only on this untested limb, and reconcile with `HORMANDER.md` explicitly (§2, §4).

- **`PROOF-TRANSPORT-ANALYTIC.md`** (banked) — already proves the (Kato) analyticity of `f₁²(β)` is **sign-blind**: it transports the *value*, not the *sign*; the massive free scalar has a simple, isolated lightest `0⁺⁺` and still `f₁²(:φ²:)=0`. "Analyticity + gap ⇏ atom." This is the same defect in a different coordinate (§2, §3).

- **`REDUCTION-AUDIT.md`** (banked) — the load-bearing logical fact: **`gap ⇏ atom`**. A gapped generalized free field, and the gapped free scalar with `O=:φ²:` (KL measure a cut on `[4m²,∞)`), have `m>0` yet `f₁²=0`. The chain `gap ⟺ H1` is FALSE as a biconditional; only `atom ∧ isolation ⇒ gap` (one way) and `gap ⇒ χ≤M₀/λ` (a susceptibility bound) survive.

- **`A4 / INVENT-GAUGE-INVARIANT-MASS.md`** (banked tool) — defines *exactly the object this species would bound*: `m_GI := inf supp ρ_s = lim_{|x|→∞} −∂_{|x|} ln|G_s(x)|`. The "analytic distance / spacelike-decay rate" the 11th species hopes to lower-bound **is** `m_GI = inf supp ρ_s`. A4's stated residual: "the edge atom `f₁²>0`." So the *object* is already registry A4; the open content is the atom.

- **`AQFT-CAUSAL.md`, `H1-T2-AQFT.md`** (banked) — the dictionary `gap ⟺ exponential clustering ⟺ mass shell ⟺ KL atom`, and the conclusion that the **only** mass-sensitive datum is the clustering rate `inf supp ν = m`, which is H1. The 11th species would re-derive the *same* identity with Hörmander's machinery (§2(b)).

- **`GOLDS-POS-MASS.md`** (banked) — "positivity + analyticity alone cannot force an atom at any mass; the free massless field satisfies every OS/Wightman axiom with `ν` purely continuous, `f₁²=0`." Directly caps any analytic-only lever.

**Therefore:** the *microlocal WF_A / Sato–Kawai–Kashiwara propagation limb is genuinely untested*, and §1–§5 attack it. Do **not** re-derive `HORMANDER.md`'s subelliptic no-go; instead show the analytic limb is its **principal-symbol twin** and fails by the same principal-vs-lower-order mechanism, plus the new atom-blindness and phase-space anti-localization wrinkles.

---

## §1. The analytic-wavefront / propagation setup, and what it *would* give

**The analytic wavefront set.** For a distribution `u ∈ 𝒟'(Ω)`, `x₀ ∉ singsupp_A u` if `u` is real-analytic near `x₀`; `(x₀,ξ₀) ∉ WF_A(u)` (`ξ₀≠0`) if there is a neighbourhood in which the **FBI transform** `Tu(x,ξ,λ) = ∫ e^{−λ(x−y)²/2} e^{−iλ ξ·(y−x)} u(y) χ(y) dy` decays like `e^{−λδ}` uniformly for `ξ` in a conic neighbourhood of `ξ₀` (Hörmander, *CPAM* **24** (1971) 671–704; *Analysis of Linear PDO I*, Grundlehren 256, §8.4). Equivalently `WF_A(u)` is the support of the analytic singular spectrum `SS(u)` of the hyperfunction `u` (Sato–Kawai–Kashiwara, LNM 287 (1973) 265–529; Kashiwara–Schapira, *Sheaves on Manifolds*, Grundlehren 292 (1990), Ch. IX). The distinction from `WF` (the `C^∞` wavefront set) is that `WF_A` is defined by *analytic* (not smooth) test regularity.

**The propagation theorem (the lever the species wants).** For `P` of **real principal type** with principal symbol `p = σ(P)`, the **bicharacteristic flow** `H_p` acts on `T*Ω`, and (`Hörmander` Vol I Thm 8.5.1 for `WF`; the analytic version, Hörmander *CPAM* 1971 / Vol I §9.6):

```
   WF_A(u) \ WF_A(Pu)  ⊂  Char(P) = {p = 0},   and WF_A(u) is invariant
   under the bicharacteristic flow of p on Char(P).                              (PROP)
```

For **nonlinear** `P` the analogue holds modulo `C^∞`/analytic smoothing by the Bony paradifferential calculus (Bony, *Ann. Sci. ENS* **14** (1981) 209–246). The structural content: **the *location* of singularities (which `(x,ξ)`), once the characteristic set and the flow are fixed, is governed entirely by the principal symbol `p`.**

**What it would give — the imagined lower bound.** The two-point Wightman function of a gauge-invariant operator `O` has the Källén–Lehmann form (Refs: Källén 1952, Lehmann 1954; `mass-spectrum/REPORT.md`; `G1-CONTINUUM-OS`)

```
   W_O(x−y) = ⟨Ω|O(x)O(y)|Ω⟩ = ∫₀^∞ dρ_O(m²) Δ₊(x−y; m),   ρ_O ≥ 0,   supp ρ_O ⊂ [m₁²,∞).   (KL)
```

At spacelike separation `r = √(−x²)`, `Δ₊(x;m) = (m/4π²r) K₁(mr) ~ C·e^{−m r}/r^{3/2}` — **exponential decay at rate `m`.** The *distance from the real `r`-axis to the nearest singularity of `W_O`* (equivalently the width of the tube of holomorphy, equivalently the abscissa of the Laplace transform of `W_O`) is therefore

```
   dist_an(W_O)  =  inf supp ρ_O  =  m₁   =  A4's  m_GI = lim_{|x|→∞} −∂_{|x|} ln|G_s(x)|.   (PW)
```

**This is the property the species would certify.** If a microlocal argument could bound `dist_an(W_O) ≥ m₀ > 0` *from below, non-circularly, uniformly in `(L,a)`*, that would be a mass gap. The whole question of this file is whether `(PROP)` supplies such a bound. It does not — for three independent, structural reasons (§2).

---

## §2. The crux — elliptic / nonlinear obstruction, and the KL-edge circularity

### (a) The propagation input is structurally unavailable, *and* mass-blind where available

**(a-Euclidean) Elliptic ⇒ no bicharacteristics ⇒ only regularity.** The Euclidean YM system `D_μ F^{μν}=0` (on a gauge slice) and the YM Hessian `Hess S_W = −D² + 2 ad F_B` have a **principal symbol that is elliptic transversally** (`σ_ξ = ξ²·Π`, `Π` the transverse gauge projection): `Char(σ) ∩ (T*ℝ⁴_ℝ) = ∅` over the real covector bundle. **An elliptic operator has no real bicharacteristics** — `(PROP)` is vacuous (there is nothing to propagate along). What ellipticity *does* give is **analytic hypoellipticity**: `WF_A(u) ⊂ WF_A(Pu)` (with real-analytic coefficients), i.e. a **regularity / upper bound on singular support**, *not* a lower bound on any distance. This is precisely the direction `HORMANDER.md` identified: **regularity (`k→∞`) ≠ the spectral bottom (`k→0`).** The gap lives at the soft mode `k→0`, i.e. `ξ→0` — where elliptic propagation says *nothing* (it is a statement about `ξ→∞`, large covectors).

**(a-Minkowski) Nonlinear hyperbolic ⇒ no solution theory.** The reconstructed Minkowski field equation `∂_μF^{μν} + [A_μ,F^{μν}] = 0` is a **nonlinear** hyperbolic system. `(PROP)` for nonlinear `P` (Bony) requires (i) a (weak) solution to exist and (ii) a priori paradifferential regularity to run the calculus. **Neither exists for 4D quantum YM** — that is the `G1` existence wall itself (`PROOF-SKELETON` §2). So the standard propagation input is *unavailable* for the very theory whose gap is at issue.

**(a-principal) Even granting a Minkowski solution, `(PROP)` is principal-symbol — hence mass-blind.** For the linearized operator `□ + m²` (the free UV limit), the principal symbol is `p = ξ²`, `Char = {ξ²=0}` (the light cone), **identical for `m=0` and `m>0`**; the mass `m²` is a **lower-order (subprincipal)** term that enters only the *transport equation* (the amplitude of the singularity along the bicharacteristic), *never* the flow `H_p` that `(PROP)` controls. So the *geometry* of the analytic wavefront set — where the singularities are, and along which curves they propagate — is **mass-independent**. The mass can only show up in the *amplitude* of the singularity, which is a functional of the *boundary/initial data* — i.e. of the measure `ρ`, i.e. of `inf supp ρ`. Which is §2(b).

**Result of (a):** On Euclidean, `(PROP)` degenerates to regularity (upper bound, wrong direction, `k→∞`). On Minkowski, `(PROP)` has no solution to act on (G1), and even with one it is a principal-symbol statement blind to `m²`. **The lever is either empty (Euclidean) or unavailable+blind (Minkowski).** This is the exact analytic twin of `HORMANDER.md`'s "subellipticity is `k→∞` regularity; the gap is the `k→0` bottom."

### (b) The quantity that *does* carry the mass is the KL edge — circular

The only mass-sensitive datum in `(KL)` is the **edge** `inf supp ρ_O`; by `(PW)` that edge *is* the tube width / `dist_an`. Two faces, same identity:

**(b-i) "Lower-bound the distance" = "lower-bound the KL edge" = H1.** To prove `dist_an(W_O) ≥ m₀` is, by `(PW)`, to prove `inf supp ρ_O ≥ m₀` — a **lower bound on the support of the spectral measure at its bottom edge**, i.e. *exactly* the statement that the lightest state in the channel has mass `≥m₀`. That is the mass gap. In the coarse-mode variable, `inf supp ρ` is tied to `G″(k→0)=1/Var(a)` by the Ward–Stein identity (`ABC-NEW-INEQ`; `OTHER-H1-POINCARÉ`), i.e. to `H1`. The lower bound is not a *consequence* of an independent microlocal input; it **is** the target.

**(b-ii) The steelman "analytic singularity of `ν` ⇒ exponential decay of the FT ⇒ gap" is the same circle.** A natural escape: use the *non-stationary-phase / analytic-regularity ⟹ decay* theorem (Hörmander Vol I, Paley–Wiener–Schwartz, Thms 7.3.8/7.4.3): if the spectral measure `ν` is real-analytic (no singular support) on the interval beyond `m₀`, then the Fourier transform `W_O` decays like `e^{−m₀|x|}`. Does *that* give a lower bound? **No — it is the identical identity read backwards.** The exponential rate one obtains **equals** the distance from the real axis to the nearest singularity of `ν`, which is `inf supp ν = m₁`. To conclude a *positive* rate you must *already know* `ν` has no singular support below `m₀` — which is the gap. And the free massive scalar shows the analyticity hypothesis is satisfied by a gapped theory with a continuum `ρ` on `[4m²,∞)` (analytic off the cut) — so analyticity of `ν` is *compatible with* `f₁² = 0`. This is `PROOF-TRANSPORT-ANALYTIC`'s "analyticity transports the value, not the sign," and `GOLDS-POS-MASS`'s "positivity + analyticity cannot force an atom," in Hörmander coordinates. It is also, verbatim, `AQFT-CAUSAL` §5 / `MODULAR-BLIND` §2: the "decay rate = clustering rate = `inf supp ν` = H1."

**(b-iii) The atom-blindness — even a magic distance bound is the wrong object.** Suppose — *contra* (b-i)/(b-ii)* —* a genuine non-circular microlocal lower bound `dist_an(W_O) ≥ m₀` existed. It would certify a **mass** `m_O > 0`. But the Clay object is the **atom** `f₁² = ν({m₁²}) > 0` for `O = tr F²`, and `gap ⇏ atom` (`REDUCTION-AUDIT`): the massive free scalar has `dist_an(W_φ) = m > 0` yet `f₁²(:φ²:) = 0`. So the species, *even in its best case*, targets the **necessary** condition (a gap) and not the **sufficient** one (the residue `I1 ∧ I2 = H*`, strictly stronger than `m>0`; `G2-RESIDUE-ATOM`, `RESIDUE-CONDITIONAL`). It loses the residue content — the exact loss `REDUCTION-AUDIT` diagnoses in the whole "gap ⟺ atom" chain.

**Verdict of §2.** The microlocal lever is **not genuine as a gap-forcer**: `(PROP)` is empty on the elliptic (Euclidean) side and unavailable/blind on the nonlinear (Minkowski) side; the mass-bearing quantity is the KL edge `= H1` (circular in two independent ways); and even in the best case it certifies the gap, not the Clay atom (atom-blind). This is the **analytic principal-symbol twin** of `HORMANDER.md`: there the principal symbol `ΣX_a²` was blind to the lower-order Savvidy curvature; here the principal symbol `ξ²` is blind to the lower-order mass `m²`. **Same mechanism, complexified.**

---

## §3. The two acid tests

### ACID TEST (1) — avoid measure-positivity *and* anti-localization; genuine LOWER bound?

Three sub-questions; the species wins one, loses two.

- **Avoids measure-positivity — YES (in form).** `WF_A`, `SS`, `(PROP)` are distribution-theoretic / complex-analytic; no variance, no LSI, no FKG, no `Γ₂`, no RP. This is a genuinely **different category** from the six dead positivity species and from the modular species (10th). Like Cheeger (curvature-free) and AQFT (operator-algebraic), the *syntax* is new. Real, but — as everyone from `CHEEGER-CORE` to `MODULAR-BLIND` found — a new syntax is not a new lever.
- **Avoids anti-localization — NO (§4).** The anti-localization principle (`TRANSPORT-CORE` §4) survives *in phase space*: `(PROP)` propagates singularities on the **characteristic variety** (`ξ²=0`, the light cone — a *high-frequency* locus), while the mass edge is read at **spacelike infinity / `k→0`** (`ξ→0`, off-cone). The machinery lives precisely where the mass is *not* read. The only proven positive inputs (high-field WSTC starvation, Haar `W_geo` floor) are, as always, anti-localized relative to the "typical-field × spacelike-decay" cell where the edge is set. (§4.)
- **Genuine LOWER bound — NO (circular).** By §2(b), "lower bound on `dist_an`" **is** "lower bound on `inf supp ν`" **is** the gap `= H1`. There is no independent perimeter/distance to bound. The output *anatomy* is right (a positive width, mass-like — better than the variance-lower-bound / decay-upper-bound of the dead positivity species), but the *value* is the open problem. (This is the same "right anatomy, circular value" verdict as `TRANSPORT-CORE`'s `Γ₂`.)

**Result: (1) FAIL** — new category (form) but no anti-localization escape and no non-circular lower bound (content).

### ACID TEST (2) — free fields (the species' genuine, and only, strength)

- **Free MASSLESS (must give `m=0`, no spurious gap).** `W(x)` of the massless field: `ρ(m²)=δ(m²)` at `m=0`; `Δ₊(x;0) = 1/(4π²(−x²))` — **singular exactly on the light cone `x²=0`**, `WF_A(Δ₊) = {(x,ξ): x²=0, ξ ‖ x}`; the spacelike decay rate is `0` (power law `1/r²`, no exponential). So `dist_an = inf supp ρ = 0`. **Returns `m=0`, no spurious gap — correctly failed, non-vacuously, in the continuum** (not via compactness). ✔
- **Free MASSIVE (must give tube width `m`).** `ρ(m²)=δ(m²−m²)`; `Δ₊(x;m) ~ e^{−m r}/r^{3/2}` for spacelike `r→∞`; `dist_an = m > 0`. **Returns `m`, correct scaling.** ✔
- **Free field `f₁²=0` razor (the kill for sufficiency).** Take the **massive** free scalar with `O = :φ²:`. The theory is gapped (`m>0`), and `W_φ` has tube width `m>0` (so the tool *would* report a positive analytic distance). Yet **`f₁²(:φ²:) = 0`** — its KL measure is the two-particle cut `[4m²,∞)`, no atom (`REDUCTION-AUDIT`; `PROOF-TRANSPORT-ANALYTIC` §3; `BS-RP-STRUCTURE`). **The analytic-distance tool is therefore atom-blind: it certifies a gap on a theory whose relevant residue is exactly zero.** This is *the* decisive acid-test finding: the species passes the required "fail the massless field" test cleanly, but its output does not distinguish the atom from the continuum — the precise `gap ⇏ atom` failure.

**Result: (2) PASS on the literal criterion** (free massless → `0`, non-vacuous, in the continuum; free massive → `m`, correct) — **genuinely stronger than `HORMANDER.md` (anti-correlated) and `TRANSPORT-CORE` (lattice `ρ≥N/12` vacuous)** — **but atom-blind** (massive free `:φ²:`), which is why it still fails #2/#3.

---

## §4. Anti-localization — the phase-space mismatch

The anti-localization principle (`TRANSPORT-CORE` §4; sharpened for modular data in `MODULAR-BLIND` §4.5) is: *every positive compensation is largest on the complement of the gap cell* (typical field × `k→0`). For analytic microlocal analysis it reappears **in phase space**, in its cleanest form:

| object | lives where | gap cell | mismatch |
|---|---|---|---|
| `(PROP)` bicharacteristic flow | `Char(□)= {ξ²=0}`, high `\|ξ\|` (light cone) | `ξ→0`, spacelike (`ξ²<0`) | disjoint covector regions |
| `WF_A` / `SS` singularities | conic near `ξ₀≠0` (FBI, `λ→∞`) | spacelike-infinity decay (`λ, ξ→0`) | disjoint scales |
| high-field WSTC / `W_geo` Haar floor | `‖F‖≫κ_typ` / `\|F\|`-free | typical field `‖F‖≈κ_typ` | (the `TRANSPORT-CORE` wall, unchanged) |
| analytic hypoellipticity (Euclidean) | `ξ→∞` regularity | `k→0` bottom | opposite ends of the spectrum |

Two independent mismatches per row: in `ξ`-magnitude (the propagation is a large-`ξ` statement; the mass is a `ξ→0` datum) and in `ξ`-direction (propagation is *on* the light cone `ξ²=0`; the mass is read *off* it, in the spacelike region `ξ²<0`). The one *proven* positive input that could in principle reach a lower bound — the quartic coercivity / WSTC high-field starvation — is (as in every species) anti-localized to the **cone**, while the analytic-singularity edge at spacelike infinity is set at **typical field**. So the coercivity cannot be inserted into a propagation argument to reach the mass-edge cell: it is largest exactly where the edge is not.

**This is the central honest point.** The species does **not** dodge `TRANSPORT-CORE` §4; it reproduces it in a new coordinate (covector phase space instead of field-magnitude). The "positive" structures a microlocal proof would lean on — ellipticity/analytic-hypoellipticity, coercivity, the characteristic-set geometry — are all anti-localized relative to the `ξ→0`, spacelike, typical-field cell where `inf supp ν` is decided. Combined with §2(a), the picture is exact: **the high-`ξ` machinery (propagation) and the low-`ξ` datum (mass) are disjoint, and the only bridge between them — `(PW)` — is the circular KL-edge identity.**

---

## §5. Scorecard (0–5), sharpest line, exact stuck step

### Five-condition grade: **2/5**

| # | condition | grade | one-line basis |
|---|---|---|---|
| 1 | **nonperturbative** | **PASS** | `WF_A`/`SS`/`(PROP)` exact distribution theory; no loop, no `1/N`, all `β`. |
| 2 | **sign-controlling LOWER bound `m_eff²(k→0)≥λ`** | **FAIL** | output *anatomy* right (a positive width, mass-like), but `(PROP)` is principal-symbol (mass-blind); the only lower bound on `dist_an = inf supp ν` **is** the gap `=H1` (§2b). |
| 3 | **non-circular** | **FAIL** | `dist_an = inf supp ν = m_GI =` KL edge `=H1` (`(PW)`, A4); non-stationary-phase "analytic ⇒ decay" is the same circle (§2b-ii); atom-blind (`gap⇏atom`, §2b-iii). |
| 4 | **non-SUSY + fails free massless field** | **PASS (clean, non-vacuous; atom-blind caveat)** | pure YM non-SUSY; free massless → `0` (continuum, non-vacuous), free massive → `m` (correct) — **stronger than `HORMANDER.md`/`TRANSPORT-CORE` on this beat**; but massive free `:φ²:` has width `m>0` & `f₁²=0` (blind to the atom). |
| 5 | **`(L,a)`-uniform `m_eff²(k→0)≥λ`** | **FAIL** | a uniform lower bound on `dist_an` **is** the uniform mass gap `=H1`; finite-`L` is vacuous; no non-circular uniform bound exists (§2b, §4). |

**Score: 2/5** (passes #1, #4; fails #2, #3, #5). Same profile as `HORMANDER.md` (2/5), `CHEEGER` (2/5), `TRANSPORT-CORE` (2/5); one beat better than the refuted modular species (1/5). It is the **most honest 2 in the corpus on the free-field axis** — it is the *only* species that tracks the free *mass* exactly on both the massless and massive sides in the continuum — and that very cleanliness is what exposes the failure: **correct free-mass tracking is not a discriminator of the atom.**

### Sharpest line

**Analytic microlocal analysis is the complexified twin of Hörmander's subellipticity, and dies by the same principal-symbol/blindness mechanism.** `(PROP)` governs the *location* of singularities through the principal symbol `p=σ(P)`; for the YM field equation `Char(□)={ξ²=0}` is *identical for `m=0` and `m>0`* — the mass `m²` is a lower-order (subprincipal) datum that never enters the bicharacteristic flow. The one quantity that carries `m` — the holomorphy-tube width at spacelike infinity, `e^{−m|x|}` — is `inf supp ν` by Paley–Wiener, i.e. `A4's m_GI`, i.e. the KL edge `=H1`; bounding it below *is* the gap (circular). On Euclidean the operator is elliptic (no bicharacteristics ⇒ only `k→∞` regularity, `HORMANDER.md`'s wrong end); on Minkowski it is nonlinear with no 4D solution theory (`G1`). And even a magic non-circular distance bound would certify the **gap**, not the **atom**: the massive free scalar `:φ²:` has width `m>0` and `f₁²=0` (`REDUCTION-AUDIT`). The species is the first with exact, non-vacuous free-mass tracking on both sides — and that is exactly why it shows the analytic distance is *the right number for the wrong question*: it measures `inf supp ν`, the gap, not `ν({m₁²})`, the atom.

### EXACT stuck step (缺失步)

```
╔════════════════════════════════════════════════════════════════════════════╗
║ target:   dist_an(W_O) = inf supp ρ_O ≥ m₀>0,  uniform (L,a)                 ║
║           ⟺  m_GI ≥ m₀  ⟺  the mass gap  =  H1.                              ║
║                                                                            ║
║ propagation lever:  WF_A(u)\WF_A(Pu) ⊂ Char(P), flow of σ(P)   (PROP)        ║
║   • Euclidean YM: elliptic, Char∩T*ℝ⁴_ℝ = ∅ ⇒ no bicharacteristics,         ║
║     only analytic-hypoellipticity = REGULARITY (upper, k→∞) — wrong end.     ║
║   • Minkowski YM: nonlinear hyperbolic, no 4D solution theory (=G1) ⇒        ║
║     (PROP) has nothing to act on; and Char□={ξ²=0} is m-INDEPENDENT          ║
║     (m² is subprincipal) ⇒ flow is mass-blind.                               ║
║                                                                            ║
║ mass carrier:  dist_an(W_O) = inf supp ρ_O  (Paley–Wiener / KL / A4.m_GI)    ║
║   → a lower bound on it IS the KL edge IS the gap IS H1  (circular, ×2:      ║
║     direct, and via non-stationary-phase "analytic⇒decay").                  ║
║                                                                            ║
║ atom:  even dist_an≥m₀>0 certifies only the GAP; the Clay object is           ║
║   f₁²=ν({m₁²})>0 for tr F².  Massive free :φ²: has dist_an=m>0, f₁²=0.        ║
║   gap ⇏ atom  (REDUCTION-AUDIT).  Missing:  H* = I1∧I2 (residue+isolation).  ║
║                                                                            ║
║ missing:  a microlocal estimate on the SUBPRINCIPAL / Paley–Wiener WIDTH     ║
║   that bounds inf supp ν from below using dynamics — which is precisely the   ║
║   clustering-rate = inf supp ν = H1.  No propagation theorem bounds a width. ║
╚════════════════════════════════════════════════════════════════════════════╝
```

**The new no-go to bank (this species' genuine contribution):** *Propagation of analytic singularities is a **principal-symbol (characteristic-variety) theorem**: it fixes **where** singularities are and **along which bicharacteristics** they move, governed by `σ(P)`; the mass is a **lower-order / subprincipal** datum read off the **width of the holomorphy tube at spacelike infinity** (= `inf supp ν`, Paley–Wiener), which no propagation theorem bounds from below. Hence the analytic microlocal species is structurally incapable of a non-circular gap lower bound — the complexified `HORMANDER.md`.* Secondary no-go: *analytic distance is atom-blind (`gap ⇏ atom`, massive free `:φ²:` is the witness), so even a non-circular distance bound would target the necessary condition, not the Clay residue.*

### What would change the verdict (honest)

Not a sharper propagation estimate — a **category the species does not possess**: (i) a *nonlinear* analytic-microlocal theorem that bounds the **transport amplitude** (the subprincipal reduction) by the **dynamics** rather than by the initial KL measure — but the amplitude *is* a functional of the measure's edge, so this is H1 in microlocal dress; or (ii) a same-channel coercion `dist_an ≳ ... ` that yields `f₁²` (not merely `inf supp ν`) — i.e. `H* = I1∧I2`, the atom, which the analytic distance cannot see. Neither is a missing estimate; both are the open problem. Absent one, the 11th species is closed as a **faithful dialect of the KL edge = H1** — the analytic principal-symbol twin of `HORMANDER.md`, with one genuine (atom-blind) free-field win.

---

*Refs.* **Internal (read/cited as ground, not invented):** `HORMANDER.md` (subellipticity = `k→∞` regularity, principal part `ΣX_a²` = free Laplacian's, sign-blind, bracket condition elliptic — **the real-axis sibling; this file is its analytic/complexified twin**); `PROOF-TRANSPORT-ANALYTIC.md` (analyticity transports the *value* not the *sign*; massive free scalar simple isolated `0⁺⁺` but `f₁²(:φ²:)=0`); `REDUCTION-AUDIT.md` (`gap ⇏ atom`; gapped free `:φ²:` has `f₁²=0`; `H* = I1∧I2` strictly stronger than `m>0`; `χ≤M₀/λ` the only one-way theorem); `INVENT-GAUGE-INVARIANT-MASS.md` A4 (`m_GI = inf supp ρ_s = lim −∂_{|x|}ln|G_s|` = the analytic distance; residual = edge atom); `AQFT-CAUSAL.md` / `H1-T2-AQFT.md` (gap ⟺ clustering ⟺ mass shell ⟺ KL atom; only mass-sensitive datum = clustering rate = `inf supp ν` = H1; free massless net's Wightman axioms); `GOLDS-POS-MASS.md` (positivity + analyticity cannot force an atom); `G2H-BLIND-CLEAN.md` (free massless satisfies every Wightman axiom, `inf supp ν=0`; granted Wightman existence insufficient); `BS-RP-STRUCTURE.md` (KL threshold `(2m_g)²` conditional on `m_g>0`=H1; RP ⇒ upper only); `TRANSPORT-CORE.md` §4 (anti-localization of every positive compensation; the gap cell = typical field × `k→0`); `MODULAR-BLIND.md` (the dead 10th species; "different category" + "anti-localization" already-banked pattern; the twin structure this file re-uses); `CHEEGER-CORE.md` (the 9th-species house style: curvature-free syntax ≠ lever); `PROOF-SKELETON.md` §0,§2 (single species = variance upper bound, UV=L3/IR=H1; no 4D YM net = G1); `ABC-NEW-INEQ` / `OTHER-H1-POINCARE` (Ward `G″·Var(a)=1`, coarse `G″=1/χ`); `mass-spectrum/REPORT.md` (KL representation `W_O(x−y)=∫dρ_O Δ₊(x−y;m)`, zero-spatial-momentum kernel `e^{−m|t|}/(2m)`); `G1-CONTINUUM-OS` (KL existence from OS0/1/3); `lean-infra/REPORT.md` (tube-domain theorem / edge-of-wedge = mathlib G-8/G-9 gap; MG needs no Wightman reconstruction — a pure-Euclidean statement). **External (verified):** Hörmander, "Uniqueness theorems and wave front sets for solutions of linear differential equations with analytic coefficients," *Comm. Pure Appl. Math.* **24** (1971) 671–704, DOI 10.1002/cpa.3160240506 (analytic wavefront set `WF_A`, FBI transform, propagation); Hörmander, *The Analysis of Linear Partial Differential Operators I*, Grundlehren 256 (Springer, 1983), §8.4 (def. `WF_A`), §8.5 (propagation, Thm 8.5.1), §9.6 (analytic version), Thms 7.3.8/7.4.3 (Paley–Wiener–Schwartz); Sato–Kawai–Kashiwara, "Microfunctions and pseudo-differential equations," in *Hyperfunctions and Pseudo-Differential Equations*, Lecture Notes in Math. **287** (Springer, 1973) 265–529 (hyperfunctions, singular spectrum `SS`); Kashiwara–Schapira, *Sheaves on Manifolds*, Grundlehren 292 (Springer, 1990), Ch. IX (microlocal / singular support); Bony, "Calcul symbolique et propagation des singularités pour les équations aux dérivées partielles non linéaires," *Ann. Sci. ENS* **14** (1981) 209–246 (nonlinear propagation, paradifferential); Treves, *Introduction to Pseudodifferential and Fourier Integral Operators*, Plenum (1980); Källén, *Helv. Phys. Acta* **25** (1952) 417; Lehmann, *Nuovo Cim.* **11** (1954) 342; Streater–Wightman, *PCT, Spin and Statistics, and All That* (1964) (tube/edge-of-the-wedge analyticity); Osterwalder–Schrader, *CMP* **31** (1973) 83 / **42** (1975) 281; Savvidy, *Phys. Lett. B* **71** (1977) 133; Nielsen–Olesen, *Nucl. Phys. B* **144** (1978) 376; Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737. *(Citation honesty: Sato–Kawai–Kashiwara LNM **287** (1973) 265–529 confirmed via live MathSciNet/Springer listings (MR0420735). The Hörmander *CPAM* **24** (1971) 671–704 reference is the standard, consistently-cited one (DOI 10.1002/cpa.3160240506); the live search backend errored on that query, so it is asserted from canonical knowledge, not live-verified — flag for the next citation-audit pass. All other external refs are corpus-pre-vetted classics.)*

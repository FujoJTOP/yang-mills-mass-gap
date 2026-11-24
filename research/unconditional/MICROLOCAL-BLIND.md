# MICROLOCAL-BLIND — independent blind red-team on the "analytic microlocal analysis / wavefront set (Hörmander, Sato–Kawai–Kashiwara)" species for the 4D pure SU(N) mass gap

> 2026-11-23 · **判词：SPECIES DEAD — and essentially ALREADY-BANKED. Survival score 2/5 (generous; true proof-content = 1/5, the `HORMANDER` row).** The proposal: use the **analytic** wavefront set `WF_A` and **propagation of analytic singularities** (Hörmander; Sato–Kawai–Kashiwara) to *lower-bound the distance from the real axis to the nearest analytic singularity of the two-point function* — claimed `= m`, the gap — *without measure positivity*. It dies because that distance **is** the Källén–Lehmann spectral edge, and **lower-bounding it is the gap** (KILL 1 lands hardest, as a two-horn fork: keep `ρ≥0` and it is the KL identity; drop `ρ≥0` and "distance = mass" is meaningless — the advertised "no positivity" destroys the route's own object). It dies again because microlocal **propagation is a principal-symbol / characteristic-variety statement**, while `m` is the **lower-order** mass-shell datum — the same wall that banked `HORMANDER` (subellipticity is regularity at `k→∞`, the gap is the sign of the indefinite lower-order Savvidy curvature at `k→0`) — and because the machinery has **no input** (Euclidean YM is elliptic; Minkowski YM has no rigorous nonlinear solution theory — constructing it *is* Clay — and the nonlinear interaction terms are distributions of the *same* singular order, not the analytic coefficients propagation needs; meanwhile Osterwalder–Schrader already encodes the analyticity, and OS is *insufficient*: the free massless field satisfies every Wightman axiom with `inf supp ν = 0`, `G2H` §5). It dies a third time on **object mismatch / anti-localization**: the sign lives on the measure-theoretic cell (typical field × `k→0`, `TRANSPORT-CORE` §4), while `WF_A` is a phase-space support datum of a *single* correlator — it cannot see the coarse-mode variance where the sign is decided. KILL 4 **does not** land on its literal clause — the tool returns distance `0` on the free massless field (no `>0` leak, a clean acid test, real credit) — but it lands as **vacuity**: the functional is a *read-out* of the spectral edge (it returns `m>0` equally for the free **massive** field, the same family that fails `f₁²(:φ²:)=0`), so it manufactures no atom. Exact stuck step: a quantitative strip-width / analytic-singularity lower bound `δ ≥ m` *off the characteristic set* — which is `inf supp ν ≥ m²` = the KL edge = H1/H* itself.

---

## §0. The tool, stated precisely (kill the real thing, not a strawman) — and the prior version (Step 0)

For a gauge-invariant interpolator `O`, the two-point Wightman / Schwinger function `W(x)=⟨Ω|O(x)O(0)|Ω⟩_c` is (by Källén–Lehmann, `ρ≥0`) the boundary value of a function analytic off the spectral cuts, with

```
W(x) = ∫₀^∞ dν(s) · Δ_+(x; √s),     supp ν ⊂ [Δ², ∞),  Δ² = inf supp ν,
```

and the **nearest singularity** of the analytic continuation (in `|x|`, or in the complexified invariant `p²`) sits at `√s = Δ = m`. The proposal (stated at full strength):

```
Claim: m = dist( real axis , nearest analytic singularity of W );
       propagation-of-analytic-singularities (WF_A)
       lower-bounds that distance, hence m, WITHOUT using ρ≥0 / reflection positivity.
```

This must deliver (the five standing constraints): (1) nonperturbative; (2) a sign-controlling **lower** bound `m_eff²(k→0)≥λ`; (3) non-circular (not assume the gap); (4) non-SUSY **and fails the free massless field** (`f₁²=0`, `G″=0`); (5) `(L,a)`-uniform, reaches `β≈2.5`.

**Prior version (Step 0, `Grep` of the corpus).** Two in-repo faces are this species already:

- **`HORMANDER.md`** — the **C^∞ / subelliptic** sibling: microlocal *regularity* (Hörmander bracket / hypoellipticity, `H^{1/r}` gain). Its kill: the gauge frame is a global frame ⇒ elliptic `r=1` ⇒ **same principal part as the free Laplacian, which is maximally hypoelliptic and gapless**; subellipticity is `k→∞` regularity, the gap is `k→0` spectral bottom. The **analytic** `WF_A` version is a genuinely *different* object (analytic, not smooth, singularities; Sato–Kawai–Kashiwara, not Hörmander brackets) — so I attack it on its own terms below, using `HORMANDER`'s structural no-gos as corroboration, not as a substitute.
- **`INVENT-NEW-IDENTITY.md`** — the **nearest-singularity face already constructed and banked**: "`f₁²` is the residue at the smallest-magnitude zero of the KL spectral (Fredholm) determinant `D_O(z)=det(1−zT_O)` … free field is a cut, no zero ⇒ residue `0`; YM has an isolated zero *below* its cut ⇒ residue `>0`." That file's own verdict: it **relocates** the wall ("from lower-bound a usc atom to evaluate a determinant zero in the continuum") — the continuum evaluation *is* the analytic-continuation / QFT-existence wall. The analytic-microlocal proposal is `INVENT-NEW-IDENTITY`'s object ("nearest singularity = `m`") dressed in `WF_A` language, with the determinant-evaluation replaced by a propagation claim.

So this is the union of the two banked faces. §1–§4 are the four kills; each maps to a hard constraint.

---

## §1. KILL 1 — KL-edge circularity, as a **two-horn fork**: it is the spectral edge, and the advertised "no positivity" destroys the mass interpretation

**Lands — decisively, and it is the primary kill.**

The identity the species rests on — *"distance from the real axis to the nearest analytic singularity of `W` `= m`"* — **is** the definition of the Källén–Lehmann spectral edge. `Δ² = inf supp ν` is, by construction, exactly the point where the analytic continuation of `W` first fails to be analytic (the lower endpoint of the spectral measure's support is where the branch point of `∫dν(s)/(p²+s)` sits). This is already banked: `G2H-BLIND-CLEAN` §1 ("for its KL measure `ν_O` we need `inf supp ν ≥ Δ² > 0`"), `REDUCTION-AUDIT` (the whole "mass or susceptibility" section). **Lower-bounding the distance-to-nearest-singularity is lower-bounding `inf supp ν`, which is the gap.** The microlocal vocabulary renames it; it does not create slack.

The proposal's distinguishing advertising — *"without measure positivity"* — turns this into a **fork** in which both horns are dead:

- **Horn A (keep `ρ≥0`).** Then `W = ∫dν/(p²+s)` with `ν≥0`, the nearest singularity is at `s=Δ²`, and "lower-bound the distance" is literally "prove `Δ²>0`" = the gap. **Zero independent leverage.** And the one quantitative object propagation could conceivably supply — a strip of real-analyticity of half-width `δ` — has `δ = ` the spectral edge by the KL representation (`W(x)~∫e^{−√s|x|}dν(s)` is analytic only for `|x|` where the exponential converges, i.e. up to the light of the lowest mass). Lower-bounding `δ` **is** lower-bounding `m`. The propagation theorem does not bound `δ`; the KL integral does, and its width is the answer.

- **Horn B (drop `ρ≥0`, as advertised).** Then the Källén–Lehmann representation **fails**: there is no positive spectral measure, `inf supp ν` is not defined as a physical threshold, and the "nearest analytic singularity of `W`" carries **no mass interpretation at all** — an arbitrary complex distribution has analytic singularities wherever it likes, decoupled from any Hilbert-space spectral edge. The route then computes a quantity that **is no longer the mass**. The species' entire selling point (escaping positivity) destroys the very identification (distance = mass) that made it a candidate. **Self-defeating.**

Either horn → dead. **The claimed gain is the claimed bug.** There is no third horn: the "distance = mass" identity is *conditional on* KL positivity, and "lower-bound the distance" *conditional on* KL positivity *is* the gap.

**And the propagation cannot rescue it, for a principal-symbol reason (the `HORMANDER` inheritance).** Propagation of analytic singularities (Hörmander / Duistermaat–Hörmander; Sato–Kawai–Kashiwara) is governed by the **principal symbol / characteristic variety** of the operator `W` solves. For the wave / Klein–Gordon operator `□+m²`, the principal symbol is `ξ²` — the **light cone**, *identical* to the free massless operator — while `m²` sits in the **lower-order** ("mass") term. So the characteristic singularities that propagate sit at **distance `0`** (the cone), universally, free field included; the **mass-shell** analytic singularity (distance `m`) is **off-characteristic** and is *not* a bicharacteristic-flow datum — it is fixed by the global solution (the spectral weight), not by the principal-symbol propagation. This is exactly `HORMANDER.md` §5 ("the gap is not a principal-symbol datum but the sign of the indefinite lower-order Savvidy curvature; a sum-of-squares / principal-symbol tool is blind to it") in analytic dress. **Verdict on Kill 1: LANDS.** Beats #2 and #3 FAIL.

---

## §2. KILL 2 — no propagation input: the equation has no analytic coefficients, the QFT does not exist, and OS already grants (and is insufficient for) the analyticity

**Lands.**

**(a) Euclidean YM is elliptic — the wrong symbol class.** The rigorous YM object in 4D is the **Euclidean** (Schwinger / lattice) theory. Analytic propagation of singularities is a **hyperbolic / real-principal-type** theorem (bicharacteristic flow at finite speed on a Lorentzian characteristic variety); the Euclidean operator `−D²+…` is **elliptic** (no real characteristic variety to flow along). So on the only rigorous side (Euclidean) the machinery has no bicharacteristic to propagate. (This is `HORMANDER` §1–§2 again: elliptic, step 1, same as free.) One must go to **Minkowski** to get a real characteristic variety — and there the theory does not exist (next point).

**(b) Minkowski interacting YM has no rigorous solution theory — constructing it *is* Clay.** To propagate `WF_A(W)` you need `W` to *solve* a hyperbolic equation in the analytic / distributional category. The interacting 4D YM Wightman QFT — an operator-valued distribution satisfying `∂_μ F^{μν}=…` at the needed regularity — **is precisely what the mass-gap problem asks you to construct** (Clay; `PROOF-SKELETON` §0, Clay `0%`; `G2H-BLIND-CLEAN` frame: existence is *granted* as a hypothesis, never proved nonperturbatively). So the object whose wavefront set you would propagate **is the unconstructed QFT**: the input is unavailable, circular at the existence level.

**(c) The interaction terms are not analytic coefficients — the propagation hypothesis fails at the only place it could bite.** Even granting a Minkowski solution formally, the equation of motion for the correlator is
```
⟨(D_μF^{μν})^a(x) O(y)O(0)⟩ = (contact terms) + ⟨(nonlinear interaction terms) O(y)O(0)⟩.
```
Hörmander's propagation theorem requires the **lower-order coefficients to be `C^∞`/real-analytic**, so that the *principal symbol alone* controls the singularity flow. But the "interaction terms" are **products of the same field operators** — distributions of the **same singular order** as the unknown `W`, not analytic background. So the nonlinearity sits exactly where control is lost: you cannot treat it as an analytic lower-order potential, and the propagation estimate has no hypothesis to stand on. The one linear piece — the free `□` acting on `W` — has, as in §1, the *free* (distance-`0`) characteristic variety.

**(d) OS reconstruction already encodes the analyticity — so microlocal adds no axiom, and the OS-granted content is insufficient.** Osterwalder–Schrader reconstruction (the spectral condition `H≥0`) already yields the tube-domain analyticity of the Wightman functions ("primitive analyticity"): `W` is the boundary value of a function analytic in the forward tube. That *is* the analytic-singularity structure microlocal analysis would discuss. So `WF_A` supplies **no new analyticity input** — it is a restatement of the OS spectral condition, which is *assumed*. And the OS-granted analyticity is **provably insufficient**: the free massless field satisfies **every** Wightman axiom with `inf supp ν = 0` (`G2H-BLIND-CLEAN` §5; `ATOM-adversarial`: "RP 是标度协变的（自由 Maxwell 即 RP 且无隙）"); the spectral condition `H≥0` permits singularities as close to the real axis as `0`. So the analyticity that OS grants cannot force `m>0`, and microlocal propagation, being a consequence of that same analyticity, cannot either. **Verdict on Kill 2: LANDS.** Beat #1 nibbled (analyticity is real but already banked by OS), Beats #3/#5 FAIL on input.

---

## §3. KILL 3 — anti-localization / object mismatch: `WF_A` is a single-correlator phase-space datum, the sign lives on a measure-theoretic field×momentum cell

**Lands.**

The gap's **sign** is carried, per `TRANSPORT-CORE` §4 and `PROOF-SKELETON` G2.5, on the **cell "typical field × `k→0`"** — a cell described by the *measure* `μ`: "typical field" means `‖F‖≈κ_typ`, a typical-amplitude configuration of the Gibbs/Wilson measure, and `k→0` is the coarse zero mode `a_{k→0}`. The sign is the positivity of `G″ = E[H|0] − Var(∇S|0) + W_geo` at that cell (`G2.5` open, `0%`). Every positive compensation is **anti-localized** to the complement of that cell (`TRANSPORT-CORE` §4 table: `‖Hess f‖²` at high `k`, `Φ″` surplus at the log-cut, WSTC at high field), and `HORMANDER` §4 records the same support mismatch for the regularity face.

Now, `WF_A(W)` is a **closed conic subset of `T*M∖0`** attached to a **single deterministic distribution** `W`. It reports *whether a point is an analytic-singularity and in which `ξ` direction* — a qualitative (yes/no + direction) support statement about **one correlator**. It is silent on:

- the **field-amplitude** coordinate `‖F‖≈κ_typ` (`WF_A` lives in spacetime × momentum, **not** in configuration-amplitude space);
- the **measure-theoretic** variance `Var(∇S|0)` where the sign is decided (`WF_A` sees no `dμ`);
- the **quantitative** value of `G″` (the scale of the gap — `WF_A` is qualitative and cannot supply a constant `λ>0`, only "singular vs analytic").

So the objects live in disjoint spaces: the sign is a **measure-theoretic field×IR-momentum** datum; `WF_A` is a **single-distribution spacetime×momentum** datum. Propagation is moreover a **global** flow (along the full bicharacteristic), so even its phase-space locality does not map onto the local sign-carrier cell. This is the **anti-localization obstruction** that banked the 8th species (`TRANSPORT-CORE` §4: "all three [positive terms] are largest where the gap is *not*, and all three vanish on the *one* cell … where the sign is decided"), now with sharper force: `WF_A` is not merely anti-localized, it belongs to a **different object category** and cannot address the sign cell at all. **Verdict on Kill 3: LANDS.** Beat #2 FAIL (again), on support/object grounds.

---

## §4. KILL 4 — free-field / vacuity: **no `>0` leak** (its one genuine pass), but a **read-out** of the spectral edge, not a lower bound on it

**Partially lands — I report honestly which sub-clause fires and which does not.**

**Does *not* fire: the literal free-massless acid-test clause.** On the free massless field (`O=tr F²`, two massless gluons), `ν` is purely absolutely continuous with support from `s=0` (`ρ(s)∝s²`, `INVENT-NEW-IDENTITY`; `ATOM-adversarial`), so the nearest analytic singularity is at `s=0` — the **real axis itself**, distance `0`. The tool returns **distance `0`**, i.e. it does **not** leak `>0` on the free massless field. By the charter's acid rule ("自由无质量场若不失效 ⇒ vacuous") this is a **clean non-vacuous discriminator**, and it is the species' one genuine asset (shared with the `Cheeger` / `capacity` family, `CHEEGER-BLIND` §3). Credit where due.

**Does fire: the vacuity clause ("does not merely restate the free spectral edge").** The functional "nearest analytic singularity distance" is a **read-out** of whatever spectral edge is already in `W`:

- On the free **massive** field it returns `m>0` — but that is the mass one *put in*; it manufactures nothing. This is verbatim `TRANSPORT-BLIND` §3(b)'s free-field kill of Bakry–Émery: "it returns `m²>0` on the free **massive** field … it **cannot distinguish** 'free massive gap' from 'YM glueball atom' — it reports whatever dynamical mass is already in the measure and does not manufacture the interacting atom."
- On the **composite** `tr F²` in even a gaped theory the nearest singularity is the **two-particle cut threshold** and the residue (atom weight) can still be `0` — the massive free scalar with `O=:φ²:` has `ν` a pure cut on `[4m²,∞)`, `f₁²(:φ²:)=0` (`REDUCTION-AUDIT` link A; `INVENT-NEW-IDENTITY` acid test). So "nearest singularity at distance `m`" does **not** imply the Clay atom `f₁²>0`, nor the lightest-mass ordering (`REDUCTION-AUDIT` link A: "`gap ⇏ atom`; the atom + isolation ⇒ gap *in that channel*, becoming the mass gap only if that channel carries the lightest state overall — an unproved ordering input"). The distance-to-singularity is a *necessary* datum (like H1), **not sufficient** for `H*` (`RESIDUE-CONDITIONAL`: atom strictly stronger than gap).

So KILL 4 lands **as vacuity**: the tool neither leaks `>0` (good) **nor** supplies an independent lower bound (fatal) — it **restates the spectral edge it is supposed to bound**, exactly as `INVENT-NEW-IDENTITY` admits ("relocates the wall … does not dissolve the continuum problem"). **Verdict on Kill 4: LANDS on the vacuity clause; the leak clause does not fire.** Beat #4 PASS (literal acid test) but it buys only a read-out.

---

## §5. Five-point verdict (honest)

| # | Constraint | Verdict | Exact failure |
|---|---|---|---|
| 1 | Nonperturbative | **PASS** | Analytic `WF_A` / propagation is exact distribution theory, no loop, no `g`-expansion, all `β`. (The analyticity itself is OS-banked, `G2H` — real but not new.) |
| 2 | Sign-controlling **lower** bound `m_eff²(k→0)≥λ` | **FAIL** | "distance to nearest singularity" `= inf supp ν =` the gap (KL identity, §1); propagation targets the **principal symbol** (characteristic cone, distance `0`), the mass `m` is **lower-order / off-characteristic** (§1, `HORMANDER` §5); no quantitative strip-width bound off the characteristic set. |
| 3 | Non-circular | **FAIL** | KL-edge identity (§1, Horn A = tautology); "no positivity" fork makes "distance = mass" meaningless (§1, Horn B); need the Minkowski QFT to apply propagation = the unconstructed object (§2); `WF_A` re-states OS analyticity, which is *insufficient* (free massless, §2d). |
| 4 | Non-SUSY + **fails free massless field** | **PASS (read-out)** | Returns distance `0` on free massless — **no `>0` leak**, clean acid test. **But** returns `m>0` on free **massive** and on any gaped theory (cannot manufacture the atom); vacuous as a *proof* (§4). |
| 5 | Uniform `≥ρ` reach `β≈2.5` | **FAIL** | No `β`-dependence; Minkowski QFT does not exist at `β=2.5` to propagate (§2b); uniformity in `(L,a)` is H1 verbatim; continuum evaluation of the analytic singularity = `INVENT-NEW-IDENTITY`'s wall. |

**Score = 2/5** (conditions 1, 4 pass). It is a **generous 2**: the #4 pass is a clean but **vacuous read-out** (the identical functional certifies the free *massive* "gap" it did not prove), so the species' true *proof-content* is the `HORMANDER` `1/5` row (nonperturbative only) crossed with `INVENT-NEW-IDENTITY`'s banked "nearest singularity = atom residue, wall relocated to the continuum."

---

## §6. Why this is a re-tread, not a new species (decisive prior art)

The proposal is the **analytic face** of two already-banked objects, in no new class:

- **`HORMANDER.md`** — the C^∞ sibling: "**Wrong object; acid-test anti-correlated.** … subellipticity is `k→∞` regularity, the gap is `k→0` spectral bottom … the gap is the sign of the indefinite lower-order Savvidy curvature, which a nonnegative sum-of-squares is blind to." The analytic version **inherits the same principal-symbol blindness** (§1, §2) minus the anti-correlation (it correctly returns `0` on the free massless field, hence a slightly higher score than `HORMANDER`'s 1/5, not a new mechanism).
- **`INVENT-NEW-IDENTITY.md`** — the **exact same object** ("`f₁²` is the residue at the smallest zero of the KL spectral determinant; free = cut/no-pole ⇒ `0`; YM = isolated pole below its cut ⇒ `>0`") already **constructed**, with its own honest verdict that it **relocates the wall** to "evaluate a determinant zero in the continuum" — i.e. to analytic continuation = the QFT-existence wall the microlocal route also hits (§2b).
- **`REDUCTION-AUDIT.md` / `G2H-BLIND-CLEAN.md`** — the "ridge" identity: gap ⟺ `inf supp ν ≥ Δ²`; OS grants the analyticity and the free massless field satisfies **all** axioms with `inf supp = 0`. So "the analytic singularity sits at `m`" is the OS-granted edge, **insufficient** by itself.
- **`TRANSPORT-CORE.md` §4 / `PROOF-SKELETON` G2.5** — the sign lives on the measure-theoretic typical-field × `k→0` cell; `WF_A` is a single-correlator phase-space datum; **object mismatch** (§3).

The proposal's only "addition" over these — *"analytic (not C^∞) singularities; propagation, not regularity"* — is a genuine distinction from `HORMANDER` (`WF_A` ⊋ `WF`) but is **non-load-bearing**: analytic propagation is still **principal-symbol-governed** (so blind to the lower-order mass, §1) and still **needs the solution theory** (so has no input at `β=2.5`, §2). Same lock, a more refined key.

---

## §7. What, if anything, survives (the honest residue)

- **Genuine, already-banked, not new:** (i) the analytic-singularity classification of `W` (a real, exact fact — but it is the OS spectral condition, `G2H`); (ii) the **clean free-massless discriminator** (distance `0` on the free field — non-vacuous, shared in spirit with `Cheeger`/capacity); (iii) the `INVENT-NEW-IDENTITY` **nearest-singularity = atom residue** (pole below cut) identity, exact at fixed `a`, every `L`.
- **A faithful dialect, not a closure:** "distance to the nearest analytic singularity `= m`" is a **restatement of the KL edge** in microlocal vocabulary; it is **lossless**, so it cannot weaken H1 (`H1 ⇒ gap, H1 unreduced`). It lets one cite Hörmander / Sato–Kawai–Kashiwara propagation machinery; the machinery is real, the leverage is nil.
- **The exact missing OBJECT** is a **quantitative, off-characteristic analytic-singularity lower bound**: a real-analyticity strip half-width `δ ≥ m` for `W` (equivalently `inf supp ν ≥ m²`, uniform in `(L,a)`) that does *not* already presuppose the KL width — which is H1, and in the residue/atom form H* = `I1∧I2` (`REDUCTION-AUDIT`; `RESIDUE-CONDITIONAL`). The microlocal route does not supply it; it re-expresses the demand for it.

---

## §8. Sharpest line

**The analytic-microlocal species dies on the identity it was built around: "the distance from the real axis to the nearest analytic singularity of `W` `= m`" *is* the Källén–Lehmann spectral edge `√(inf supp ν)`, so lower-bounding it is the gap — and the advertised escape ("without measure positivity") is self-defeating, because the "distance = mass" identification *requires* KL positivity (`ρ≥0`); drop it and the singularity carries no mass at all (KILL 1, two-horn fork). Propagation of analytic singularities cannot rescue this: it is governed by the **principal symbol / characteristic variety** (the light cone, distance `0`, free field included), while the mass `m` is the **lower-order / off-characteristic** mass-shell datum — the exact `HORMANDER` blindness (subellipticity: regularity at `k→∞`, blind to the lower-order Savvidy curvature at `k→0`) in analytic dress. The machinery has **no input**: Euclidean YM is elliptic, Minkowski YM has no rigorous nonlinear solution theory (constructing it *is* Clay), the interaction terms are distributions of the same singular order (not the analytic coefficients propagation requires), and Osterwalder–Schrader already encodes — and the free massless field proves insufficient — the analyticity. And it is **anti-localized to the wrong object**: the sign lives on the measure-theoretic typical-field × `k→0` cell (`TRANSPORT-CORE` §4), while `WF_A` is a single-correlator phase-space support datum that cannot see the coarse-mode variance. It is a faithful, lossless **restatement** of the KL edge crossed with `INVENT-NEW-IDENTITY`'s nearest-singularity = atom-residue identity, whose wall was already relocated to the continuum; microlocal propagation does not cross it.**

**Exact stuck step.** A quantitative, **off-characteristic analytic-singularity lower bound** for the gauge-invariant two-point function:
```
   W(x) real-analytic for  |x| < 1/δ  with  δ ≥ m  (uniform L, a),
   i.e.  inf supp ν ≥ m²  with m > 0  —  derived WITHOUT assuming the KL width —
```
equivalently the residue/isolation form `H* = I1 ∧ I2` (`f₁² ≥ r₀` + `I2`). This **is** the Källén–Lehmann edge `= H1`/`H*`; the principal-symbol propagation of `WF_A` bounds only the characteristic cone (distance `0`), never the off-cone mass-shell singularity. Producing the strip bound *is* the gap; analytic microlocal analysis does not produce it — it restates the demand for it.

---

*Refs.* **Internal (read/cited as ground, not invented):** `HORMANDER` (C^∞ subelliptic sibling: principal-symbol blindness; gap = lower-order Savvidy curvature sign at `k→0`; elliptic `r=1`; free field gapless; §3 `+‖u‖²` is the mass; §5 sign-blindness); `INVENT-NEW-IDENTITY` (nearest singularity = residue of KL Fredholm determinant; free = cut/no-pole ⇒ `0`; YM = isolated pole below cut ⇒ `>0`; wall **relocated** to continuum evaluation); `REDUCTION-AUDIT` (gap ⇏ atom; atom+isolation ⇒ gap in channel; lightest-mass ordering unproved; `inf supp ν` = the object; susceptibility ≠ mass); `G2H-BLIND-CLEAN` §1,§3,§5 (KL equivalence; free massless satisfies every Wightman axiom with `inf supp ν = 0`; moments insufficient); `ATOM-adversarial` (KL `D(0)=0 ⇏ gap`; RP/center neither prove nor disprove the atom — free Maxwell is RP and gapless); `TRANSPORT-CORE` §3.2,§4 (sign cell = typical field × `k→0`; every positive term anti-localized to its complement; lattice free-field `N/12` artifact; free continuum ⇒ `ρ=0`); `TRANSPORT-BLIND` §3(b) (tool reports whatever mass is already in the measure; cannot distinguish free-massive from interacting atom); `RESIDUE-CONDITIONAL` / `G2-RESIDUE-ATOM` (atom strictly stronger than gap; `H* = I1∧I2`); `PROOF-SKELETON` §0,§1 (Clay `0%`; G2.5 open `0%` = typical field × high momentum; dimension `G″_a=a²m²`); `docs/YANG-MILLS-MASTER.md` (the single missing species = a variance upper bound). **External (canonical; verified to the extent cited, no fabricated identifiers):** Hörmander, *The Analysis of Linear Partial Differential Operators I*, Grundlehren **256** (Springer, 1983), §8 (wave front sets; §8.6 the **analytic** wave-front set `WF_A`), and *LPDO III–IV* for propagation of singularities; Duistermaat–Hörmander, *Acta Math.* **128** (1972) 183 (propagation of singularities, real principal type); Sato–Kawai–Kashiwara, "Microfunctions and pseudo-differential equations," *Lecture Notes in Math.* **287** (1973) 265 (analytic microlocal analysis; singular spectrum); Treves, *Introduction to Pseudodifferential and Fourier Integral Operators* (Plenum, 1980); Baouendi–Goulaouic, *Comm. PDE* **3** (1978) (analytic hypoellipticity); Källén, *Helv. Phys. Acta* **25** (1952); Lehmann, *Nuovo Cim.* **11** (1954); Osterwalder–Schrader, *Comm. Math. Phys.* **31** (1973) 83 and **42** (1975) 281 (reconstruction ⇒ Wightman analyticity / tube domains); Streater–Wightman, *PCT, Spin and Statistics, and All That* (1964); Glimm–Jaffe, *Quantum Physics: A Functional Integral Point of View* (Springer, 1987); Jaffe–Witten, *Quantum Yang–Mills Theory* (Clay, 2000); Savvidy, *Phys. Lett. B* **71** (1977) 133; Shen–Zhu–Zhu, *Comm. Math. Phys.* **400** (2023) 805, arXiv:2204.12737.

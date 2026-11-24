# TOOL-COMPLETELY-MONOTONE — can the completely-monotone / Bernstein–Widder / Wiener–Ikehara structure of the KL correlator `C(t)=∫e^{−st}dν(s)` give a NONPERTURBATIVE **LOWER** bound on the support floor `m₁²=inf supp ν`?

> 2026-11-24 · Scouting a specific AREA (completely-monotone function theory + Abelian/Tauberian asymptotics) for the one direction that keeps failing: a **lower** bound on a spectral support edge. Judged by the charter: **(a)** lower-vs-upper direction, **(b)** free-field acid test (must *fail* on the free massless field), **(c)** already Family I (Schur-concave / a transform functional) or genuinely escaping, **(d)** does it relocate the `E₂`-regress / circularity. Labels: **[T]** 真定理 / **[C]** 条件 / **[NG]** no-go (named class) / **[O]** 精确缺失步 / **[J]** 判词. My own elementary checks are **[V-here]**; external facts are **[standard]/[R]** (textbook — I re-verified no external source this session and ran no lattice job). **No QED on anything; this document does not prove Clay and does not prove no tool can exist.**

---

## Verdict up front

> ### **REDUCES TO FAMILY I — and, in the only limit that could bite, to the `E₂`/transport circularity. NOT a new corner.**
>
> The completely-monotone structure is **exactly equivalent to `ν≥0`** (Bernstein–Widder) — i.e. to Källén–Lehmann / reflection positivity itself — so **every** finite-data functional it generates is a decreasing **UPPER** bound on `m₁²`, never a lower bound. There is **no rigidity hidden in "all derivatives alternate"** beyond positivity: the alternation is *equivalent* to `ν≥0`, and positivity only ever says "the mass is *below* this number." The one object that could give a floor — the `t→∞` exponential rate `−(1/t)log C(t) → m₁²` — is an **exact Tauberian identity**, not a lever: to *use* it you must evaluate `C(t)` on the scale `t ≳ 1/m₁`, whose uniform control **is** exponential clustering = the gap (circular), and the resulting floor datum is the unprotected residue `f₁²` = `H*`. The analytic-continuation / radius-of-convergence / Stieltjes-fraction / Padé twin of the question is **floor-blind by a structural reason**: the edge **atom contributes an ENTIRE function** `f₁²e^{−m₁²t}` and so creates *no singularity*, while the half-plane of convergence `Re t>0` is fixed by the UV extent (support to `+∞`), independent of `m₁²`. Wiener–Ikehara's "main term" is the leading exponential `m₁²/f₁²`; stating its regularity hypothesis requires subtracting the answer. **Nothing here escapes; the area's contribution is a sharper *negative* (temporal localization of the obstruction), not a lever.**

**The three posed sub-questions, answered directly:**
1. *Sharper than elementary Laplace Tauberian?* **No [T] + [NG-class].** `CM ⟺ ν≥0` *is exactly* the hypothesis under which the Laplace Tauberian is an **identity** (`−(1/t)log C(t) → m₁²`); there is nothing to sharpen because the identity is already exact.
2. *Bernstein–Widder / ROC / singularity / Stieltjes-fraction / Padé / derivative-growth lower bound, non-circular?* **No [NG].** Atom is entire ⇒ analytic structure floor-blind; derivative growth at `t=0` is the moment ratio (**upper**) and UV-divergent in the continuum; the transfer-variable Stieltjes/Padé sees the atom as a *pole* but at cost `O(1/(am₁)) → ∞` = `PADE-ROUCHE-DEEP` (relocates to transport).
3. *Wiener–Ikehara analogue for the support floor?* **No [NG].** The "main term" is the leading exponential = the answer; the regularity hypothesis presupposes subtracting `m₁²`.

---

## §0 Setup — the exact object, and one definitional warning

`ν≥0` finite on `[0,∞)` is the KL spectral measure of `O=tr F²` (or any scalar interpolating operator): `ν(B)=⟨Ω|O†P_H(B)O|Ω⟩`, `supp ν ⊂ [m₁²,∞)`, edge atom `f₁²=ν({m₁²})`, moments `M_k=∫s^k dν`. The correlator is
> `C(t)=∫_0^∞ e^{−st} dν(s) = ⟨Ω| O† e^{−tH} O |Ω⟩`,  `t≥0`.

**Facts used.** `C` is **completely monotone (CM)**: `(−1)^n C^{(n)}(t) ≥ 0` for all `n≥0`, `t>0`. By **Bernstein's theorem** (Bernstein 1928/29; Widder, *The Laplace Transform* 1941) the converse holds: `C` is CM on `(0,∞)` with `C(t)<∞` **iff** `C(t)=∫e^{−st}dν(s)` for a unique `ν≥0` **[standard]**. All-derivatives-alternate is therefore **equivalent to** `ν≥0`, i.e. to the KL/RP positivity the theory *already* has — not extra dynamical content (same statement as `PROOF-POLYA-FREQUENCY` (i), in the Bernstein rather than Schoenberg form).

**Definitional warning (direction of the word "abscissa").** The mass² is the **exponential decay rate**, `m₁² = −lim_{t→∞}(1/t)log C(t)` (exact for CM/`ν≥0`; §1). It is **not** the "abscissa of convergence" of the one-sided Laplace integral `∫_0^∞ e^{−st}dν(s)` in the transform variable `t`: that abscissa is `Re t = 0` **regardless of `m₁²`** — the integral converges for `Re t>0` (the `e^{−st}` damping at `s=+∞`) and diverges on `Re t≤0` (the support reaches `+∞`), never sees the bottom edge. (Why this matters is §3: it is exactly why radius/analytic-continuation methods cannot locate `m₁²`.) Also recall the statement split from `SUPPORT-NOT-ATOM`: Clay's `Δ>0` is **support** (`m₁²≥Δ²`), and the massive free `:φ²:` realizes it with `f₁²=0` — so the target is the *floor*, not the atom.

---

## §1 Q1 — CM structure is not sharper than the Laplace Tauberian (it *is* the Tauberian hypothesis)

**Theorem (Tauberian identity for CM)  [T/standard].** For `C(t)=∫e^{−st}dν(s)`, `ν≥0`:
> `− (1/t) log C(t) → inf supp ν = m₁²`  as `t→∞`.

This is the classical Abelian–Tauberian statement made exact by `ν≥0` (the `ν≥0` "one-sidedness" *is* the Tauberian side condition; `HISTORICAL-TOOL-SWEEP` #12 already records that it is an exact equivalence, not an improvable bound). The whole **content** of the theorem is supplied by `ν≥0`.

**Consequence  [T].** Bernstein's theorem says `ν≥0 ⟺ C` CM. Therefore the completely-monotone hypothesis adds **nothing** to the hypothesis "`ν≥0`" under which the Tauberian identity already holds. There is **no sharper support-floor criterion derivable from the CM structure**, because the exact identity is the sharpest possible statement (`m₁²` is literally the limit). Any *stronger* criterion could only come from a **side condition beyond CM** — a known moment, a growth/type bound, or a UV bound — and all of those are either (i) finite linear data (Family I, `LEVER §0B`), or (ii) the uniform UV/scheme control (= `IN-9`/transport), i.e. the banked wall, not a new object.

> **Label: [T] that the identity is exact; [NG] for the class "CM-only support-floor criteria sharper than Tauberian" — none exist, since CM = `ν≥0` = the Tauberian side condition.**

---

## §2 The direction theorem — every finite-`t` CM functional is a decreasing UPPER bound

The CM toolkit does produce canonical scalar functionals of `C` (the "local rate," the "tilted mean," the log-derivatives). Here is where they point.

**Theorem (rate is an upper bound, monotone to the floor)  [T; proof [V-here]].** Let `ν≥0`, `supp ν⊂[m₁²,∞)`, `M₀=C(0)∈(0,∞)`, `S∼ν/M₀`. Define `r(t):=−(1/t)log(C(t)/M₀)=−(1/t)log E[e^{−tS}]`. Then for every `t>0`:
> **`m₁² ≤ r(t) ≤ M₁/M₀ = E[S]`**,  and `t↦r(t)` is non-increasing, `r(t)↓m₁²` as `t→∞`.

*Proof.* (i) Lower bound: `s≥m₁² ⇒ e^{−ts}≤e^{−tm₁²}` pointwise, so `C(t)≤M₀e^{−tm₁²}`, hence `r(t)≥m₁²`. (ii) Upper bound: `x↦e^{−tx}` convex ⇒ Jensen `E[e^{−tS}]≥e^{−tE[S]}`, hence `r(t)≤E[S]=M₁/M₀`. (iii) Monotone: `log C` is convex, because `C''C−(C')² = [∫s²w][∫w]−[∫s w]² ≥ 0` (Cauchy–Schwarz, `w=e^{−ts}dν≥0`) so `(log C)''≥0`; `r` is the (sign-flipped) secant slope of the convex `log C` from `t=0`, decreasing to the asymptotic slope `−(log C)'(∞)=m₁²`. ∎

**So every canonical finite-`t` scalar of the CM machine — `r(t)`, the tilted mean `−C'(t)/C(t)=E_t[S]∈[m₁²,M₁/M₀]` (the "Lévy drift" when `C` is log-CM), the tilted moments `(−1)^n C^{(n)}(t)/C(t)` — is an UPPER bound on `m₁²`, tightening downward to it.** A lower bound on `m₁²` (what a gap proof needs: `m₁²≥Δ²`) is the *opposite* inequality; none is produced.

**The data floor-blindness lemma (exponential-moment form of `LEVER §0B`)  [T; witness [V-here]].**
*Single-time witness.* Fix a time `t` and an observed value `v=C(t)/M₀=E[e^{−tS}]∈(0,1)`. For **any** floor `m*∈[0, −(1/t)log v]` there is a probability measure on `[m*,∞)` with `E[e^{−tS}]=v` (atom of mass `p=v\,e^{tm*}` at `m*`, remainder pushed to `+∞`). Hence the single datum `C(t)` confines the floor only to `[0, −(1/t)log v]` — an **upper** endpoint, **no lower endpoint.** (Numerically: `t=1`, target `E[e^{−S}]=0.3`: floors `0` (atom `0.3`+tail), `0.6`, and `1.204` (single atom) all give the same `C(1)`.) *General statement:* finitely many times `0<t_1<…<t_k` impose `k` real linear constraints on the infinite-dimensional convex set of positive measures; the solution set is infinite-dimensional, so `inf supp` is unbounded below across it. This is the **Laplace/exponential-test-function** face of the master **linear-data lemma** and of the semicontinuity fact (`PROOF-NEVANLINNA` Mechanism 2: the floor/atom is only *upper*-semicontinuous, while `ν↦∫e^{−ts}dν` is a bounded-continuous functional).

> **Label: [T].** No functional of finitely many exact values/derivatives of `C` lower-bounds `m₁²`; the complete alternation adds no constraint (it equals `ν≥0`), so it cannot beat the lemma. **This is `LEVER §0B` / `SCHUR-CONVEX` in the time-domain coordinates — i.e. Family I.**

**UV remark [R/internal].** For `O=trF²`, `M₀=⟨O²⟩` and all `M_k` (`k≥1`) are UV-divergent; use the point-split/normal-ordered correlator (finite small-`t`), which changes no direction — `r` and `E_t[S]` remain upper bounds, and the floor sits in the `t→∞` rate.

---

## §3 Q2 — analytic continuation / ROC / singularities / Stieltjes fraction / Padé / derivative growth

The natural "lower-bound-shaped" classical statement is *minimum-modulus / growth ⇒ first singularity* (`HISTORICAL-TOOL-SWEEP` Corner A). Translated to `C(t)`, it would read "the growth of the analytic continuation of `C` lower-bounds `m₁²`." Here is the structural reason it cannot, in exactly the `ν≥0`/CM setting.

**Structural theorem (atom is entire; analyticity is floor-blind)  [T, elementary].**
(a) An **atom** at `m₁²` contributes `f₁² e^{−m₁² t}` to `C(t)`. This is **entire in `t`** — it creates **no singularity anywhere** in the complex `t`-plane. (b) The a.c. part `∫_{edge}^∞ g(s)e^{−st}ds` (physical `g(s)~s²`) converges absolutely iff `Re t>0` **and is analytic throughout `Re t>0`**, for *every* value of the support floor; its only obstruction is the line `Re t=0` (the UV extent `s→+∞`), *independent of `m₁²`*. Hence:
> the **region of analyticity of `C` in the `t`-plane is `{Re t>0}` regardless of `m₁²`**; the **Taylor radius of convergence of `C` about any real point `t₀>0` is exactly `t₀`** (distance to the line `Re t=0`), **also independent of `m₁²`**; and the bottom atom is **analytically invisible** (entire).

**Consequences for the four candidate mechanisms:**
- **Radius-of-convergence / singularity-continuation.** Blind to the floor: `ROC=t₀` for a gapped atom, a gapped continuum, **and** the free massless field alike (all have `Re t>0` as the convergence half-plane, since all reach UV). Growth-of-continuation sees only the **UV extent**, never the IR floor. (This is the `MS-ANALYTIC-NOZERO` / `PROOF-NEVANLINNA` "min-modulus needs zero-free / pole-vs-cut" wall, sharpened: the atom is not merely hard to force, it is **entire**.)
- **Derivative growth at `t=0`.** `(−1)^n C^{(n)}(0+)=M_n=∫s^n dν`. The full growth gives `lim M_n^{1/n}=\sup supp` (the **top**) and the **moment-ratio lemma** `m₁² ≤ M_n/M_{n−1}` (`GAP-G2-edge`) — again **upper**. In the continuum `M_n=∞` for `trF²` (UV), so the germ at `0` does not even exist. Tilted version: `C^{(n)}(t)/C(t)` are moments of the **tilted** measure `e^{−ts}dν`, which has the **same** support/floor — so the "derivative growth at `t₀`" is the same Family-I data about a same-floor measure; no new lower bound.
- **Stieltjes continued fraction / Padé/Christoffel in the continuum variable `t`.** Same fate: continuum moments diverge; and an S-fraction/Padé of the Laplace transform would be reconstructing from the (nonexistent) Taylor germ, with the atom an entire contribution it cannot localize. The classical **Christoffel** bound at a point is `f₁² ≤ 1/K_n` — an **upper** bound on the atom (`INEQ-T1`).
- **Stieltjes fraction / Padé in the lattice *transfer* variable** `λ=e^{−aE}`, where `C(ak)=∫λ^k dν(λ)` **are** finite power moments and the atom **is** a pole `1/(1−zλ₁)`. Here the machinery *does* see the atom — but it needs `M≳2/(am₁)` exact uniform moments and **collapses as `a→0`** (`PADE-ROUCHE-DEEP`, the `λ`-gap width `≈½am₁→0`). This is the **only** regime where the atom is a visible pole, and its cost *is* the marginal-BPHZ/transport log-sum = the banked wall. **Relocates, does not cross.**

> **Label: [NG].** Every `ν≥0`-native analytic/singularity/continued-fraction criterion is either **floor-blind** (atom entire; half-plane fixed by UV) or **upper-bound** (moment ratio / Christoffel), and the one atom-visible form (transfer-variable Padé) **relocates** to transport (§ `PADE-ROUCHE-DEEP`).

---

## §4 Q3 — Wiener–Ikehara analogue for the support floor

**Wiener–Ikehara  [standard/R].** For a Dirichlet series `A(s)=Σ a_n n^{−s}` with `a_n≥0` (or `∫a(x)x^{−s}dx`), if `A(s)−c/(s−1)` extends continuously/across `Re s=1`, then the partial sums obey the main term `Σ_{n≤x}a_n ~ c\,x`. The engine: **positivity + a boundary-regularity hypothesis ⇒ an asymptotic main term.**

**Translate to the floor.** The analogue of the "main term" is the **leading exponential** of `C(t)`: for an isolated edge, `C(t)=f₁² e^{−m₁²t}(1+o(1))`, so `−(1/t)log C(t)=m₁²+o(1)` and the "amplitude" is `f₁²`. An "Ikehara-for-the-floor" would read: *if `e^{m^2 t}C(t)` (some subtraction) has a regular boundary/short-time behavior, then `m₁²` equals the growth rate.* Two closures:
- **The regularity hypothesis is the answer in disguise.** To state "`C(t)/e^{−\bar m t}` is regular / has a limit," you must already possess `\bar m=m₁²` to subtract it; checking the hypothesis *presupposes* the floor. (This is the Abelian/Tauberian direction that is already free: `supp⊂[m₁²,∞) ⇒ C(t)≤M₀e^{−m₁²t}`.)
- **It outputs the answer, not a bound.** Ikehara produces the **main term** = the floor itself; it is an *equivalence* under `ν≥0`, i.e. verbatim the §1 Tauberian identity. An equivalence is not a lower-bound engine (same status as `HISTORICAL-TOOL-SWEEP` #12).

In the CM/one-sided setting there is no external "regular variation" side condition available to make Ikehara bite harder: the physical side condition is `ν≥0` **already fully used**, and it upgrades the Abelian direction to the exact identity. There is no remaining slack for an Ikehara-type *inequality* in the floor direction.

> **Label: [NG].** The Wiener–Ikehara analogue either presupposes `m₁²` (circular) or restates the exact Tauberian identity (no bound). It does not lower-bound the floor non-circularly.

---

## §5 Acid test + five-point scorecard (best candidate = the "CM rate / Tauberian-at-finite-`t`" tool)

**Free-field acid test — the decisive failure [V-here, numerics].** For the **free massless** field (ν a.c. on `[0,∞)`, floor `0`), `r(t)=−(1/t)log E[e^{−tS}]` is **strictly positive at every finite `t`** and only decreases to `0` as `t→∞`. Explicit (γ-type density on `[0,∞)`): `r(t)=3.65,\,3.13,\,2.20,\,1.30,\,0.61,\,0.27` at `t=0.1,0.3,1,3,10,30`. So the finite-`t` CM rate is **positive on the free massless field** — it can never be a gap certificate (it does not "fail"; it is *positive where it must be 0*). Only the **`t→∞` limit** returns `0`, and taking it *is* knowing the answer. (Contrast the requirement: a valid tool must be **0 on free massless, `>0` on YM** at the *same finite data*. The CM rate is the reverse on the free field, and in the limit is non-constructive.)

| criterion | (1) nonperturbative | (2) sign → **LOWER** floor | (3) non-circular | (4) acid test (fails free field) | (5) `(L,a)`-uniform | verdict |
|---|---|---|---|---|---|---|
| `r(t)` / Tauberian rate at finite `t` | form ✓ | **✗ UPPER** (`r≥m₁²`) | ✓ at finite `t` but then gives no floor | **✗ positive on free massless** | ✗ | **1/5** |
| Tauberian rate *in the limit* | ✓ | ✓ (exact) | **✗ limit = the answer** | ✓ (→0 free) | **✗** | **circular** |
| continuum Padé/S-fraction | ✓ | ✗ (moments `=∞`) | — | — | ✗ | **0/5** |
| transfer-variable Padé (atom = pole) | ✓ | ✗ (upper/absence) | ✗ (`M∼1/(am₁)`) | ✓ | ✗ | **2/5** (`PADE-ROUCHE-DEEP`) |
| Wiener–Ikehara analogue | ✓ | ✗ (main term = answer) | ✗ (subtract answer) | — | ✗ | **0–1/5** |

---

## §6 The one honest refinement — a sharper *negative*, not a lever

This area's real contribution is to **localize the obstruction in the time variable**, which is sharper than the generic "positivity is Schur-concave":

1. **The entire CM toolkit lives at finite `t`, and every finite-`t` scalar is a decreasing upper bound on `m₁²`** (`r(t)↓m₁²`, tilted mean `↓m₁²`). So the wall is not "the CM functionals are the wrong cone" only — it is **temporal**: the needed lower bound lives purely in the **`t→∞` rate**, and at every finite time the data points the *other* way.
2. **The edge atom is ENTIRE, so analytic structure is structurally floor-blind** (ROC `=t₀` independent of `m₁²`; half-plane fixed by UV). This unifies `MS-ANALYTIC-NOZERO`, `PROOF-NEVANLINNA` (pole-vs-cut), and `PADE-ROUCHE-DEEP` (atom visible only in the collapse-prone transfer variable) under one CM/singularity statement.
3. **The "sharper-than-Tauberian" question is closed exactly**: `CM ⟺ ν≥0` is *the* hypothesis that makes the Laplace Tauberian an identity, so no CM-only criterion can be sharper.

This is **negative knowledge, honestly certified** — it does **not** manufacture a floor. It refines *where* the wall sits (in `t→∞`), and why the analytic twin (Padé/continuation) cannot reach it.

---

## §7 The exact stuck step [O]

To extract a floor from `C` you must do one of:
- **(i) take the `t→∞` rate** — requires `C(t)` controlled on the scale `t ≳ 1/m₁`, uniformly in `(L,a)`, i.e. **`C(t)≤M₀e^{−Δt}` = exponential clustering = the gap** (`⇔ m₁²≥Δ`, since `m₁² = inf_t r(t)`); **circular**; or
- **(ii) continue the transfer-variable Stieltjes/Padé to resolve the atom's pole** — requires `M≳2/(am₁)→∞` exact uniform moments = the marginal-BPHZ/transport log-sum `Σ_j g_j²` (`W3`), **not uniform as `a→0`**; or
- **(iii) supply a uniform anti-concentration / stability inequality** turning finite data into a floor bound — which is the quantitative form of the atom's upper-semicontinuity and **is** the uniform UV/scheme control `IN-9`/transport.

All three are the **banked wall**, not a new sub-target:
> **(i) ≡ uniform `I2`/gap; (ii) ≡ `W3` transport; (iii) ≡ `H1`+UV control.** The CM area produces **no** independent fourth option. The exact missing step is unchanged: the **unprotected edge residue / uniform decay rate carried across the Savvidy crossing** (`H*`), i.e. the transport machine that does not yet exist (`CLAUDE.md §5`).

**Flip condition (for honesty).** This verdict flips only if one exhibits a **non-linear functional `Φ(C)`** that (i) lower-bounds `m₁²`, (ii) is NOT the `t→∞` limit (uses only data at some finite window), and (iii) is `0` on the free massless field and `>0` on YM. By §2 (`LEVER §0B`) no functional of finitely many values/derivatives qualifies; by §3 no analytic/continuation functional qualifies (atom entire). Such a `Φ` would have to exploit a **provable YM-specific nonlinear constraint beyond `ν≥0`** — the "6th species" — and the only imaginable one (strict log-convexity / non-atomness) is quantitative-empty (§2: strictness of `r>m₁²` can be arbitrarily small). None is known.

---

## §8 Honest boundary

- **Does not prove Clay; does not prove no tool can exist** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`, `CONVERGENCE-redteam §d`). It proves: within the **completely-monotone / analytic-continuation / Tauberian** class, every criterion is either an upper bound (finite data), floor-blind (atom entire), or circular (the `t→∞` rate). That is a **[NG] for this named class** + **[O]** exact missing step, not a barrier theorem.
- **Not a re-discovery of a dead file dressed anew**: `PROOF-POLYA-FREQUENCY` used the *Schoenberg/total-positivity* face (≡ `ν≥0`); this document uses the *Bernstein–Widder + singularity/entire-atom + Wiener–Ikehara* faces and adds the **direction theorem** (`r(t)↓m₁²`, finite-`t` upper-only), the **atom-is-entire / ROC-floor-blind** structural reason, and the **Wiener–Ikehara circularity** — none of which is stated as such in the corpus. But it lands on the **same wall** (`SCHUR-CONVEX` / `LEVER §0B` / `PADE-ROUCHE` / `PROOF-NEVANLINNA`).
- **Honest bottom line: reduces to Family I / the `E₂`/transport circularity. No genuine untried lower-bound tool in this area.** The candidate that exists (`r(t)`) is an **upper** bound and is **positive on the free massless field at finite `t`** — it cannot certify a gap.

---

## References

**External — [standard]/[R] (textbook; not re-verified this session):**
- S. Bernstein (1928/29); D. V. Widder, *The Laplace Transform* (Princeton, 1941) — completely monotone ⟺ Laplace transform of a positive measure; Widder's representation theorems.
- G. H. Hardy & J. E. Littlewood, Tauberian theorems; N. Wiener (1932), *Tauberian theorems*; S. Ikehara (1931) — the Wiener–Ikehara theorem (positivity + boundary regularity ⇒ main term).
- I. J. Schoenberg (1951), Pólya-frequency / positive-definite additive kernels (the TP face; `PROOF-POLYA-FREQUENCY`).
- T. J. Stieltjes (1894); de Montessus de Ballore — continued-fraction / Padé for Stieltjes transforms (`PADE-ROUCHE-DEEP`).
- Paley–Wiener (support ↔ exponential type / analyticity; the half-line unbounded case = decay-rate, **not** a `t`-plane singularity).
- M. G. Krein & A. A. Nudelman, *The Markov Moment Problem and Extremal Problems* (moment mass interval; `PROOF-NEVANLINNA`).

**Internal (absolute paths; the banked verdicts this document classifies against):**
- `D:\Dev\math\CLAUDE.md`; `D:\Dev\math\docs\BRIEF-FOR-MODEL.md`; `D:\Dev\math\research\unconditional\HISTORICAL-TOOL-SWEEP.md` (§2 #12 Tauberian = equivalence; Corner A Hadamard/Jensen = Family I in disguise); `D:\Dev\math\research\unconditional\E2-REGRESS-CIRCULARITY.md` (§0.4 SSC `E₂=ρE₁`; §5 exact missing step).
- `LEVER-nonconvex-edge.md` — **§0B linear-data lemma**; **(4)** "complete monotonicity does not forbid the massless atom (Widder)" — the direct predecessor of this file.
- `PROOF-POLYA-FREQUENCY.md` — CM ⟺ `ν≥0` via Schoenberg/Karlin (total-positivity face; free-field no separator).
- `PROOF-NEVANLINNA-POLE.md` — pole-vs-cut; atom upper-semicontinuous; germ fiber convex; forcing subclass = residue–level-gap lemma (`Λ` needed).
- `INEQ-T1-hankel-stieltjes.md` — Christoffel/ Hankel/Nevalinna are **upper**-bound machines; lower end of the fiber is `0`.
- `PADE-ROUCHE-DEEP.md` — transfer-variable Stieltjes/Padé sees the atom as a pole but needs `M≳2/(am₁)→∞` (relocates to `W3`/transport).
- `MS-ANALYTIC-NOZERO.md`, `PADE-ROUCHE-DEEP.md` — minimum-modulus/analyticity circular; UV tail masks IR atom.
- `G2A-INFSUPP-DIRECT.md` — exact squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ`; floor datum = bottom-residue.
- `SUPPORT-NOT-ATOM.md` — gap is support, not atom; free massive `:φ²:` gapped with `f₁²=0`.
- `SCHUR-CONVEX.md` — master Schur/positivity [NG] class (the cone this area sits in).

*Status: **[NG] for the completely-monotone / analytic-continuation / Tauberian class + [O] exact missing step (uniform `t→∞` rate = I2, or `W3` transport, or `H1`+UV control) + [T] (CM ≡ `ν≥0`; `r(t)↓m₁²`; atom entire).** Verdict: **reduces to Family I / the `E₂`–transport circularity; no new lower-bound tool in this area.** NOT a proof of Clay; NOT a global impossibility theorem. Clay stays open.*

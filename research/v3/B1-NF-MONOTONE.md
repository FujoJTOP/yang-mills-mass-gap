# B1-NF-MONOTONE — the N_f route is **blocked**: the algebraic (two-loop) monotone is exact but perturbative and partly scheme-dependent; the *nonperturbative* monotone is a **no-go — strict monotonicity of any universal order parameter is logically incompatible with a nonempty open conformal window**, and non-strict monotonicity is insufficient to put N_f=0 on the gapped side

> **2026-11-24 · ATTACK SURFACE B.** 判词 (verdict), no QED, Clay 0%.
> **Target:** a theorem "some universal IR observable is monotone in N_f", which together with a rigorous
> Banks–Zaks (BZ) IRFP at large N_f would put N_f=0 on the gapped side **without** computing the
> strong-coupling β.
> **Result.** **(i) `[T]` algebraic.** The two-loop facts are all exactly true and are proved below:
> β₀,β₁ strictly decreasing in N_f; sign changes at N_f = 306/38 (=153/19≈8.053) and 33/2; **no positive
> zero for N_f ≤ 8**; the fixed-point coupling² `g*² = −β₀/β₁` strictly decreasing on the window with
> `d g*²/dN_f = −642/(38N_f−306)²`; the two-loop β pointwise **strictly increasing** in N_f. **(ii) `[NG]`
> nonperturbative.** *No* universal strictly-monotone order parameter can exist across the phase boundary.
> The intended logic needs **strict** monotonicity (non-strict gives only `O(0)≥0`, vacuous); but any
> universal order parameter **vanishes identically on the conformal window**, which is an **open set** —
> so it is constant (=0) on an interval and therefore **not strictly monotone** there. The two ingredients
> of the route (a strictly monotone order parameter; an open rigorous BZ window) are **mutually exclusive**.
> **(iii) `[O]`.** The hypothesis "rigorous BZ IRFP at large N_f" is itself unavailable (only perturbative,
> `[C]`); so the route fails twice over. **Strongest true monotonicity statement:** *within the two-loop
> window, the scheme-independent IRFP datum `γ_m* = 3 − Δ_{ψ̄ψ}` is monotone in N_f at leading order
> (`[T]` at O(a), `[C]` all orders), running from `γ_m*→0` at the Gaussian edge to `γ_m*∼O(1)` at the lower
> edge* — but it is **defined only where a fixed point exists**, so it cannot reach N_f=0. **Exact missing
> lemma: a strictly monotone, universal order parameter defined across the phase boundary — proved
> impossible.** The N_f route is not hard; it is **structurally closed.**

---

## §0 · Labels

`[T]` theorem · `[C]` conditional · `[NG]` no-go for a named class · `[O]` exact missing step ·
`[V]` numerical evidence · `[LIT]` literature-map statement.

Conventions (campaign / `G1-ACID`, `G1-WINDOW-BOUNDARY`): `β(g) = −β₀ g³ − β₁ g⁵` with `16π²` absorbed in
`g`, SU(3),
```
β₀(N_f) = 11 − (2/3)N_f ,     β₁(N_f) = 102 − (38/3)N_f .
```
This matches the standard MS-bar polynomials `β₀=(33−2N_f)/3`, `β₁ = (34/3)C_A² − (20/3)C_A T_F N_f − 4C_F T_F N_f`
for `C_A=3, C_F=4/3, T_F=1/2` (van Ritbergen–Vermaseren–Larin). **Notation for the small parameter.** The
fixed point of `β = −β₀g³−β₁g⁵` sits at `g*²_{abs} = −β₀/β₁` directly (the absorbed coupling); the physical
MS-bar coupling is `g*²_{MSbar} = 16π²·g*²_{abs}`. Define `a := −β₀/β₁ = g*²_{MSbar}/(16π²)` — the small
parameter. Then the one-loop mass anomalous dimension is `γ_m* = 6C_F·a + O(a²)`, and `√a = g*_{MSbar}/(4π)`.
Throughout this document `g*²` abbreviates `a` (house convention, matching `G1-ACID`).

---

## §1 · Candidate observables — what is universal and what is not

A monotonicity theorem about a **scheme-dependent** quantity is worthless. The decisive sorting:

| Observable | Universal? (scheme-indep.) | Defined at N_f=0? | Defined on the conformal window? | Monotone in N_f? | Verdict |
|---|---|---|---|---|---|
| mass gap `Δ` (as a *number*) | no (needs a scale) | yes (`Δ>0`) | yes (`Δ=0`) | **no theorem**; and `Δ` is dimensionful | fails |
| `Δ/√σ`, `m_{0++}/√σ` | **yes** (dimensionless ratio) | yes | **0/0 — undefined** (σ=0, gap=0) | no theorem; ratio would be non-monotone/divergent at the edge | **domain fails** |
| chiral condensate `⟨ψ̄ψ⟩` / `⟨ψ̄ψ⟩/F³` | value at fixed μ: no; "≠0" yes | n/a (no fermions) | ≡0 (conformal) | **NO — lattice: enhanced** at N_f=6 over N_f=2 | **counterexample** |
| `γ_m` at an IRFP → `γ_m*` | **yes** (a CFT datum: `3−Δ_{ψ̄ψ}`) | **undefined** (no FP) | yes | yes at leading order; `[C]` all orders | **domain fails** (`[O]`) |
| `g*²` (fixed-point coupling²) | **NO** (scheme-dependent value) | undefined | yes | `[T]` exact (below) | **worthless** (scheme) |
| # IR-relevant operators | yes (integer) | yes | yes | integer step-function; **constant on open sub-intervals** | **not strictly monotone** |
| `a` (conformal anomaly) | yes | yes (`a_IR=0` gapped) | yes (`0<a_IR<a_UV`) | `a_UV` increases in N_f | **phase-blind** (`A-THEOREM`) |
| β-function value `β(g)` at fixed `g` | **no** (scheme-dependent) | yes | yes | `[T]` exact pointwise (below) | scheme + **truncated only** |

Three immediate observations, each fatal to at least one candidate:
1. The **only** dimensionless universal observables that compare theories across N_f are **ratios of
   physical scales** (`m/√σ`) and **CFT data at a fixed point** (`γ_m*`, `Δ_{ψ̄ψ}`). The first is `0/0` on
   the conformal window; the second is *undefined* at N_f=0. **The two ends of any candidate monotone live
   in disjoint domains of definition** (§4).
2. The natural order parameter `⟨ψ̄ψ⟩` is **empirically non-monotone** (LSD: `⟨ψ̄ψ⟩/F³` enhanced ≥50% at
   N_f=6 vs N_f=2) — an explicit counterexample, not just an absence of proof (§5).
3. The quantities that *are* exactly monotone (`β(g)`, `g*²`) are **scheme-dependent** (and `β(g)` is the
   truncated polynomial); by the task's own criterion they cannot carry a physical uniformity theorem.

---

## §2 · The algebraic part — proved exactly `[T]`

All statements here are exact and verified by `research/v3/w7/b1_nf_monotone.py` (rational arithmetic).

**Proposition A1 (monotonicity of the coefficients).** `β₀, β₁ : ℝ → ℝ` are affine with slopes `−2/3` and
`−38/3 < 0`; both are **strictly decreasing in N_f**. ∎

**Proposition A2 (sign changes).** `β₀ = 0 ⟺ N_f = 33/2 = 16.5` (asymptotic-freedom upper edge);
`β₁ = 0 ⟺ N_f = 306/38 = 153/19 ≈ 8.0526` (Banks–Zaks lower edge, two-loop). Since `β₁(8)=2/3>0` and
`β₁(9)=−12<0`, the sign change falls at `153/19 ∈ (8,9)`. ∎

**Proposition A3 (no perturbative IRFP for N_f ≤ 8).** For N_f ≤ 8, `β₀ ≥ β₀(8) = 17/3 > 0` and
`β₁ ≥ β₁(8) = 2/3 > 0`, so
`β(g) = −g³(β₀ + β₁ g²) < 0` for every `g > 0`. The two-loop β has **no positive zero** (and only the
trivial zero at `g=0`). In particular N_f=0: `β₀=11, β₁=102`, maximally asymptotically free. ∎
(The boundary of the no-zero set is exactly `N_f ≤ 153/19`, where `β₁ ≥ 0`.)

**Proposition A4 (strict monotonicity of the fixed-point coupling on the window).** On
`153/19 < N_f < 33/2`, `β₀>0>β₁`, and the unique positive zero is
```
a(N_f) = g*²_{abs} = −β₀/β₁ = (2N_f − 33)/(306 − 38 N_f) ,
d/dN_f [a] = −642/(38 N_f − 306)² < 0 .
```
So `a` (and hence `g*²_{MSbar}=16π²a`) is **strictly decreasing**, `a → +∞` as `N_f ↓ 153/19` and `a → 0`
as `N_f ↑ 33/2` (checked against a central finite difference in the witness: e.g. at N_f=12, analytic slope
−0.028533 vs finite-diff −0.028534; `a=0.06`, `√a=0.245=` `g*_{MSbar}/(4π)`). ∎ **But `a`, being a coupling
value, is scheme-dependent — this monotonicity, though exact, is not universal and cannot be the physical
theorem.**

**Proposition A5 (pointwise strict increase of β in N_f, two-loop).** For fixed `g>0`,
```
β_{N_f+2}(g) − β_{N_f}(g) = (β₀(N_f)−β₀(N_f+2)) g³ + (β₁(N_f)−β₁(N_f+2)) g⁵
                          = (4/3) g³ + (76/3) g⁵ > 0 .
```
So at every fixed coupling the two-loop β **strictly increases with N_f** (fermion screening). This is the
sharpest exact N_f-monotonicity. It is a statement about the **truncated polynomial**, not about the
nonperturbative zero structure — the `G1-ACID` lesson: "AF + single coupling" does not forbid a zero, so a
monotone truncated β neither creates nor forbids an IRFP nonperturbatively. ∎

**Proposition A6 (leading-order monotonicity of the universal datum, `[T]`@O(a) / `[C]` all orders).** The
leading IRFP mass anomalous dimension is `γ_m* = 6C_F·a + O(a²)`, `a = −β₀/β₁` (= `g*²_{MSbar}/(16π²)`). Since `a(N_f)`
is strictly decreasing on the window (A4) and `6C_F>0`, `γ_m*` is **strictly decreasing in N_f at leading
order**, `γ_m*→0` at the Gaussian edge. Equivalently the bilinear scaling dimension `Δ_{ψ̄ψ}=3−γ_m*`
**increases** toward 3 (free) as N_f increases (§6 table: `γ_m*` = 3.33, 1.41, 0.48, 0.027 at N_f = 9,10,
12,16). `γ_m*` **is** scheme-independent (a CFT eigenvalue); the *all-orders monotonicity* is not proven —
near the lower edge `a=O(1)` and the expansion fails (`[C]`). ∎

**What §2 delivers and does not.** It delivers every algebraic fact the route needs *on the perturbative
side*: a monotone `g*²`, a monotone leading `γ_m*`, a pointwise-monotone β, exact sign changes, and the
N_f≤8 no-zero theorem. It delivers **no** bridge to N_f=0: every monotone object here is either
scheme-dependent (`g*²`, `β`) or defined **only on the window** (`γ_m*`).

---

## §3 · The nonperturbative attempt — the logical skeleton of the intended route

Write `O(N_f)` for a candidate **universal order parameter**, meaning
```
(A)  O(N_f) > 0  ⟺  the theory at N_f is gapped ;      O(N_f) = 0 on the conformal (gapless) window W .
```
The task's route is the implication
```
   [ O strictly monotone in N_f ]  ∧  [ rigorous BZ IRFP at some N_f*∈W ]   ⟹   O(0) > 0   ⟹   N_f=0 gapped.
```
The logic of the implication is **valid**, and it pins the exact requirement:

> **Lemma (monotone-order-parameter step).** If `O` is **strictly decreasing** and `O(N_f*)=0` for some
> `N_f*>0`, then `O(0) > O(N_f*) = 0`.
>
> *Proof.* Strict decreasing: `N_f < N_f' ⟹ O(N_f) > O(N_f')`; take `N_f=0 < N_f*`. ∎

So **strict** monotonicity is exactly what converts the anchor `O(N_f*)=0` into the conclusion `O(0)>0`.
Two things are now required, and each fails:

**(L1) Strictness is mandatory.** **Non-strict** (non-increasing) monotonicity gives only
`O(0) ≥ O(N_f*) = 0` — i.e. `O(0) ≥ 0`, which is **vacuous** (compatible with N_f=0 gapless). No weaker
monotonicity notion suffices.

**(L2) Strict monotonicity is impossible for an order parameter.** This is the no-go of §4.

---

## §4 · The no-go — strict monotonicity vs a nonempty open window `[NG]`

**Hypothesis of the rigorous BZ input.** "A rigorous Banks–Zaks IRFP at large N_f" means an **open interval**
`W = (N_f^c, 11N/2)` (or any open subinterval — e.g. `(16.4,16.5)`, where `g*` is small and the two-loop IR
analysis is controlled) on which the massless theory is IR-conformal. "Open" is structural: conformality at
one `g*` persists in a neighbourhood of `N_f`, and the perturbative window is open near the Gaussian edge.

**Proposition (no-go).** Let `O` satisfy (A) — in particular `O ≡ 0` on `W`. If `W` contains two distinct
points `N_f1 < N_f2`, then **`O` is not strictly monotone on any interval containing both**. Hence `O` is not
strictly monotone on `[0, N_f*]` for any `N_f* ∈ W`.

*Proof.* `O(N_f1) = O(N_f2) = 0`, violating strict monotonicity (`0 > 0` false). Any interval containing
`N_f1, N_f2` contains two equal values. ∎

**Corollary (the route is self-contradictory).** The anchor needed in §3, `O(N_f*)=0`, is supplied by the
**openness** of `W`; and that same openness forces `O` to be zero on at least two points, **destroying the
strict monotonicity the step requires**. The two ingredients of the intended argument are **mutually
exclusive**. No smoothing/redefining of `O` helps: making `O` positive inside `W` (to restore strictness)
removes (A) — it stops being an order parameter and no longer reports the gap at `N_f*`.

**Why integers do not rescue it.** An integer order parameter (e.g. # relevant operators) is a *step
function*: constant on open sub-intervals of `W`, hence **not strictly monotone** there (strict monotonicity
forbids any plateau). `W` contains many integers (e.g. 9…16), so at least two equal values.

**Why a "conformality measure" does not rescue it.** If instead `O` is a *conformality* measure (`O=0`
gapped, `O>0` conformal, **increasing** in N_f), then `O(N_f*)>0` and monotonicity give only
`O(0) < O(N_f*)` — i.e. `O(0) ∈ [0, O(N_f*))`, which does **not** force `O(0)=0` (gapped). Same failure in
the opposite direction.

**Domain-bifurcation form of the same obstruction (sharpest statement).** The two halves of any monotone are
observables of **different kinds** living in **disjoint domains**:
- on `W` (conformal): the universal data are CFT eigenvalues `γ_m*`, `Δ_{ψ̄ψ}`, defined **only** there;
- at `N_f=0` (gapped, and across the confining side): the universal data are scale ratios `m/√σ`, `m/Λ`,
  defined **only** there.
There is **no single scheme-independent observable with a domain covering both** `N_f=0` and `W`. The monotone
breaks **exactly at the phase boundary it is supposed to locate** — the route is circular in the same way
`(G1)` is.

**Status: `[NG]`.** No general N_f-monotonicity theorem can put N_f=0 on the gapped side. This is not a
limitation of technique; it is a logical obstruction.

**Independent `[O]`.** The hypothesis "**rigorous** BZ IRFP at large N_f" is itself **not available**: the
BZ fixed point is established only in the two-loop / ε expansion (`[C]`, Banks–Zaks; `G1-ACID`; `LIT-ERA-09`
E09-21). There is no nonperturbative, rigorous large-N_f IRFP theorem. So even granting the no-go were
evaded, the route's *other* leg is `[O]`. Both legs fail.

---

## §5 · Counterexample search — does any candidate survive?

We looked for a candidate that (i) is universal, (ii) is strictly monotone in N_f, (iii) is an order
parameter. **None survives; one is an explicit empirical counterexample.**

**(C1) Chiral condensate `⟨ψ̄ψ⟩`, `⟨ψ̄ψ⟩/F³` — explicit non-monotonic counterexample `[V]`.**
LSD (`arXiv:0910.2224`, PRL 104 (2010) 071601) finds `⟨ψ̄ψ⟩/F³` agrees with QCD at N_f=2 and is **enhanced
by ≥50% at N_f=6** — i.e. the condensate *increases* from N_f=2 to N_f=6 before it must *vanish* in the
conformal window. So `⟨ψ̄ψ⟩` is **not monotone in N_f**: it rises then falls. (Its *value* is also
μ-/scheme-dependent; only "≠0" is universal, and "≠0" is the gap itself.) This single lattice fact refutes
the most natural monotone.

**(C2) `m_{0++}/√σ` — structurally disqualified `[NG]`.** Both `m_{0++}` and `√σ` vanish in the conformal
window, so the ratio is `0/0` there (undefined); on the gapped side there is no monotonicity theorem, and
the ratio must diverge or misbehave as `√σ→0` at the edge. It cannot be the cross-boundary monotone.

**(C3) `γ_m*` — universal but window-only `[O]`/`[C]`.** This is the *best* candidate: `γ_m* = 3−Δ_{ψ̄ψ}` is
a genuine scheme-independent CFT datum, provably monotone at leading order (A6), plausibly monotone all
orders. It is the Miransky/Appelquist "critical anomalous dimension" language: χSB sets in when `γ_m*` grows
past `γ_c≈1` (`[C]`, ladder/BKT). **But `γ_m*` exists only inside `W`.** It can locate the *lower edge from
above*; it cannot be evaluated at `N_f=0`. It fails the domain test (the `[O]` of §4). Worse, "`γ_m*>γ_c ⟹
χSB`" is a gap-equation **criterion**, not a theorem — it presupposes the very dynamics in question.

**(C4) # relevant operators — step function `[NG]`.** Integer count, constant on open subintervals of `W`,
not strictly monotone (§4); and it is essentially a restatement of "does a fixed point exist."

**(C5) `a`-anomaly — phase-blind `[NG]`.** `a_UV` increases monotonically with N_f (more UV matter), but
`a` is phase-blind: both gapped (`a_IR=0`) and gapless (`0<a_IR<a_UV`) IR endpoints are allowed
(`A-THEOREM.md`). It cannot be the discriminator.

**(C6) The unrefereed "monotonicity theorems" in the wild are the algebraic part, circular `[LIT]`.**
Web search (2026-11-24) surfaced only non-peer-reviewed Zenodo/Academia preprints. Their "theorems" posit a
"drift constant" `c0(N_f) = κC_adj − N_f σ` that is **affine and strictly decreasing in N_f** — which is
literally the two-loop `β₀` structure of §2 restated, with the window edge defined as its zero
(`N_f^crit = κC_adj/σ`). That is **Proposition A1–A3 dressed as physics**: it *assumes* the linear form and
reports its zero. It is not a nonperturbative theorem, and it does not evade the §4 no-go (it simply declares
`W` and does not attempt to reach `N_f=0` by strict monotonicity). Mainstream sources treat N_f-monotonicity
of the gap as approximate/conjectural. **No peer-reviewed monotonicity theorem exists.**

---

## §6 · Acid test — N_f=0 (gapped) vs N_f=12 (conformal), fixed N=3

A monotone that "distinguishes N_f=0 from a conformal N_f" must output *gapped* at 0 and *gapless* at 12.

| N_f | β₀ | β₁ | two-loop IRFP `√a=g*_{MSbar}/4π` | `a=−β₀/β₁` | `γ_m*≈6C_F a` | phase (evidence) |
|---|---|---|---|---|---|---|
| **0** | +11 | +102 | **none** (`β<0 ∀g`) | — | **undefined** | **gapped** (SZZ `β<1/48` strong coupling; lattice) |
| 8 | +17/3 | +2/3 | none | — | — | confining/walking (lattice) |
| **12** | +3 | −50 | 0.245 | 0.060 | 0.48 | **IR-conformal** (lattice, *contested*: LSD vs Hasenfratz–Peterson) |
| 16 | +1/3 | −302/3 | 0.058 | 0.0033 | 0.027 | perturbative BZ (controlled) |

**The acid test verdict.** Every *continuous* candidate that could carry the signal is either:
- an **order parameter** → `≡0` on the window (N_f=12 and 16 both report `0`) → by §4 **not strictly
  monotone**, so it cannot connect to N_f=0; or
- **scheme-dependent** (`g*²`, `β(g)`) → worthless as a uniform statement; or
- **undefined on one of the two points** (`γ_m*` at N_f=0; `m/√σ` at N_f=12).

The **only** scheme-independent object that correctly reports the two phases is the **phase qualifier itself**
("gap present / absent") — which is the answer the route was supposed to *derive*, not assume. The acid test
therefore **fails for the route**: the "monotone" is either tautological (it *is* the gap) or vacuous
(scheme-dependent / undefined on one side).

**Free-field acid test (charter §1).** Take the free (massless) gauge field as the gapless conformal
reference. Any proposed monotone that "passes" by reporting `0` for the free field and `>0` for YM is merely
reporting the mass gap — the object whose positivity is the theorem. It self-cancels: **vacuous pass**, in
the charter's sense (the tool does no independent work; it presupposes the gap).

---

## §7 · Verdict

- **`[T]` (algebraic, exact).** All two-loop facts — coefficient monotonicity, sign changes at 306/38 and
  33/2, no-zero for N_f≤8, strict monotonicity of `g*²` (scheme-dependent) and of leading-order `γ_m*`
  (universal, window-only), and pointwise strict monotonicity of `β(g)` in N_f.
- **`[NG]` (nonperturbative).** **No universal strictly-monotone order parameter can exist across the phase
  boundary**: (i) strict monotonicity is *required* to conclude `O(0)>0` (non-strict is vacuous); (ii) any
  universal order parameter is `≡0` on the **open** conformal window, hence not strictly monotone there;
  (iii) integer/"conformality" variants fail identically; (iv) the two observables' **domains of definition
  are disjoint** (CFT data at N_f>0; scale ratios at N_f=0). The route's two ingredients are mutually
  exclusive. **This is the precise obstruction, not a technical gap.**
- **`[O]` (independent).** The "rigorous BZ IRFP at large N_f" leg is itself unproven (only perturbative
  `[C]`); and even a *theorem* "no IRFP at N_f=0" would give only "outside the window," which is not yet
  "gapped" (the `(iii)` gapless-non-conformal loophole, `LIT-ERA-09 §7.4`).
- **`[V]` (counterexample).** The natural monotone `⟨ψ̄ψ⟩/F³` is **empirically non-monotone** (LSD: enhanced
  ≥50% at N_f=6 over N_f=2).

**Net.** **The N_f route is blocked.** It cannot put N_f=0 on the gapped side without computing the
strong-coupling β, because the only cross-boundary monotone would have to be a strictly-monotone universal
order parameter, and that object is provably impossible (§4). What the route *can* prove is a
**within-window** monotonicity of the universal IRFP datum `γ_m*` — a real, citable `[T]`@O(a) — but it
lives entirely on the conformal side and cannot reach N_f=0.

---

## §8 · The exact missing lemma, and the precise obstruction

> **Missing Lemma (needed, and impossible).** *There exists a scheme-independent observable `O(N_f)`,
> defined for all `N_f ∈ [0, 11N/2]`, that is (a) an order parameter for the gapped/conformal transition
> (`O>0` gapped, `O≡0` on the conformal window) and (b) **strictly monotone** in `N_f`.*

**This lemma is false** whenever the conformal window is a nonempty open set. Proof: (b) alone would give
`O(0)>0` from the BZ anchor; but (a) makes `O≡0` on an open interval (≥2 points), contradicting strict
monotonicity. Dropping (a) destroys the anchor; dropping strictness makes the conclusion vacuous. **There is
no consistent choice.** The obstruction is the **openness of the gapless phase**: an order parameter is
constant (zero) across an *open* set of N_f, and "strict monotone" and "constant on an open set" cannot
coexist.

**What a *different* route would have to do instead.** Since no cross-boundary monotone exists, the only way
to use N_f information is one of:
1. A **non-monotone but sign-definite** route: prove directly that some *scheme-independent* quantity that is
   strictly positive in the window (e.g. `γ_m*`) **blows up / crosses a critical value** at a finite N_f =
   N_f^c(N) and that χSB is *equivalent* to that crossing — i.e. a genuine nonperturbative criterion, not a
   ladder guess. This is the `γ_m` route, and its gap is exactly "γ_m > γ_c ⟹ χSB, rigorously" — an `[O]`
   of gap strength (a dynamical statement about the chiral transition).
2. A **direct no-IRFP theorem at strong coupling** (a nonperturbative sign theorem for β at N_f=0) — the
   `(G1)` wall, `STAR-NO-BZ` `[O]`.
3. A **transport machine** carrying the SZZ strong-coupling gap across the crossover — the campaign's
   missing machine, `[O]`.

None of these is an N_f-monotonicity theorem; each is gap-strength. The N_f route is a **reduction to (G1)
in a different coordinate, not a bypass of it.**

---

## §9 · References (verified 2026-11-24)

| Ref | Data | Status |
|---|---|---|
| Banks–Zaks | T. Banks, A. Zaks, *On the phase structure of vector-like gauge theories with massless fermions*, NPB **196** (1982) 189–204, DOI 10.1016/0550-3213(82)90035-9 | **VERIFIED** (ADS) |
| 4-loop β | T. van Ritbergen, J. A. M. Vermaseren, S. A. Larin, *The four-loop beta-function in QCD*, PLB **400** (1997) 379–384, arXiv:hep-ph/9701390 (confirms β₀,β₁ polynomials) | **VERIFIED** (arXiv, search) |
| γ_m one-loop | R. Tarrach, *The pole mass in perturbative QCD*, NPB **183** (1981) 384–396, DOI 10.1016/0550-3213(81)90140-1; `γ_m^{(1)}=6C_F` standard (also Vermaseren–Larin–van Ritbergen PLB 405 (1997) 327; Chetyrkin PLB 404 (1997) 161) | **VERIFIED** (search; coefficient cross-checked against ladder/HQET usage) |
| LSD condensate | LSD Collab. (Appelquist et al.), *Toward TeV conformality*, PRL **104** (2010) 071601, arXiv:0910.2224 — `⟨ψ̄ψ⟩/F³` enhanced ≥50% at N_f=6 vs N_f=2 | **VERIFIED** (arXiv abs + search) |
| SZZ strong-coupling gap | H. Shen, R. Zhu, X. Zhu, CMP **400** (2023) 805–851, arXiv:2204.12737; SU(N) d=4 gap for β<1/48 | **VERIFIED** (repo `LIT-ERA-13` 13.01; Zbl 1522.81127) |
| Critical γ_m≈1 | V. A. Miransky, K. Yamawaki, *Conformal phase transition in gauge theories*, PRD **55** (1997) 5051, Erratum PRD 56 (1997) 3768 (γ_m≈1 ladder/BKT lower-edge criterion); see also Appelquist–Lane–Mahanta 1986, Cohen–Georgi 1989 | **UNVERIFIED** here (standard citation, per repo E09-24; the γ_m≈1 *criterion* is `[C]`) |

House: `G1-WINDOW-BOUNDARY` (§1 two-loop window), `G1-ACID` (AF ⇏ gapped; single coupling can have a zero),
`STAR-NO-BZ` ((G1) is dynamical `[O]`), `A-THEOREM` (a is phase-blind, upper-type), `LIT-ERA-09` (BZ `[C]`,
N_f=12 contested, §7.4 "outside window ⇏ gapped"), `LIT-ERA-13` (13.01 SZZ, 13.16/13.17 contested IRFP),
`V3-LIT-SYNTHESIS §5B` (attack surface B: "cleanest new target", no such theorem known).
Witness: `research/v3/w7/b1_nf_monotone.py` (exact rational algebra + toy no-go; runs clean).

*Status ledger:* **`[T]`** A1–A6 (algebraic); **`[NG]`** nonperturbative monotonicity (§4); **`[O]`** rigorous
BZ leg + missing lemma + (iii)-loophole; **`[V]`** LSD condensate counterexample; **`[C]`** γ_m* all-orders
monotonicity, γ_m>γ_c criterion. **Clay 0%. No QED.**

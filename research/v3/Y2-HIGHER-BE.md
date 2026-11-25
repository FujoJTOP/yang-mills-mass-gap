# Y2-HIGHER-BE — the higher-order / nonlocal Bakry–Émery rescue of the Savvidy soft mode **fails**: higher critical-point jets do not determine the gap (fixed-jet double well), the regularized curvature is not a valid BE constant (4-convex but gapless; over-claims on the exact one-plaquette), and the only nonlocal criterion that does control `λ₁` **is** the gate density `ρ(gate)` = H1

> 2026-11-24 · Attack **Y2** on the Savvidy wall. Standard BE/LSI dies because its one-sided
> constant is `ess inf G″` and the Savvidy tail makes `ess inf G″<0` (H6-B, Q4). Task: can a
> criterion with **higher covariant derivatives** (iterated `Γₖ`, `k`-convexity, regularized
> curvature `G″+αG⁗`) or a **nonlocal functional** restore a *uniform log-concavity ⇒ Poincaré*
> bound? Tested on the **exact one-plaquette** (`G″=A−Bβ²`, `β*=4.30661`) and the **Savvidy cell**
> (the double-well normal form). Reads: `Q4-BAKRY-EMERY.md`, `H6-HEAD-SPECTRAL.md`,
> `T-WAVE-SYNTHESIS.md`; model anchors `S2`/`G4`/`H3`/`H4`/`Q3`/`R5` (the `A−Bβ²` closed form);
> `B3-FUNCTIONAL-BE`, `OTHER-H1-POINCARÉ`, `CHEEGER-CORE`.
> Witness: `research/v3/w7/y2_higher_be.py` (`py -3.12 -I`, numpy+scipy, deterministic, 1-D jobs;
> output §5).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## §0 · Verdict (no QED, Clay 0%)

**No.** A criterion built from higher covariant derivatives, or from a nonlocal functional, does
**not** restore a uniform log-concavity bound past the Savvidy negative mode. Three independent
obstructions, each a named no-go, each checked on the exact models:

1. **Higher critical-point jets do not determine the gap `[T]`+`[V]`.** The fixed-jet double well
   `V_w(y)=(w²/2)((y/w)²−1)²` has, at **every critical point**, `V′(cp)=0`, barrier curvature
   `V″(0)=−2`, well curvature `V″(±w)=+4` — **independent of `w`** — and all higher derivatives
   `V^{(k)}(cp)→0` as `w→∞`; yet the barrier height `w²/2→∞` and the spectral gap **Kramers-collapses**
   `λ₁=e^{-w²/2}·(pref)`: `λ₁=9.0×10⁻¹, 3.2×10⁻¹, 1.1×10⁻¹, 9.1×10⁻³, 2.9×10⁻⁴` for
   `w=1,1.5,2,3,4` (§5·A). **Any criterion that is a function of the jet `{V^{(k)}(x)}` read at the
   critical points (or its essential infimum) is blind to `w`, hence to the gap.** The gap is the
   **nonlocal barrier action** `∫√(2(V−E))`, which no finite-order local datum contains.

2. **The regularized curvature is not a valid BE constant `[NG]`+`[V]`.** The Savvidy
   double well `V_g=(g/2)(x²−1)²` is **4-convex** (`V‴=0`, `V⁗=6g>0` everywhere) and its
   regularized curvature `κ_α:=V″+αV⁗` at the barrier is `(6α−1)g>0` for `α>1/6` (`κ₁=5g>0`); yet
   `λ₁→0` (`0.79→4.6×10⁻³` for `g=0.5→16`, §5·B). So "higher-order log-concavity" (`V⁗≥0`) and the
   regularized curvature are **satisfied on a gapless double well** — a criterion positive there
   cannot imply a Poincaré bound. On the **exact one-plaquette** the same `κ₁=G″+G⁗` **over-claims**
   (`κ₁=0.62>λ₁=0.40` at `β=10`), so it is not even a lower bound (§5·D).

3. **The only nonlocal functional that controls `λ₁` is the gate capacity `=H1` `[T]`+`[NG]`.**
   The exact 1-D Cheeger constant `h(ν)=inf_x ρ(x)/min(F,1−F)` gives the rigorous `λ₁≥h²/4`
   (verified on the one-plaquette: `h²/4=0.094,0.197,0.402 ≤ λ₁=0.340,0.414,0.404` at
   `β=2.5,5,10`). But `h` is a **barrier-action** functional: on the double well `h~ρ(gate)~e^{-ΔV}`
   (`g=0.5→16`: `ρ(gate)↓492×`, `h↓492×`, `λ₁↓173×`, §5·F). `h≥h₀ ⟺ ρ(gate)≥c ⟺` no interior
   barrier `⟺` **H1**. The nonlocal rescue is H1 renamed; there is no nonlocal object that certifies
   the gap without being a gate-density bound.

The acid tests pass and are **non-vacuous**: free massless field `λ₁=k²→0` exactly — every higher
jet is `0`, so every higher-order criterion returns `0` (§5·C); the 1-D `SU(2)` Haar marginal is a
single log-concave well, `λ₁=0.2645`, `χ≤3.78` (`β=2.5` floor `0.222⇒χ≤4.5`) — plain BE already
suffices, higher order adds nothing (§5·E).

**The three kill criteria all fire.** (a) *Circular*: the natural repair "require `ess inf G″≥0`"
demands the defect absent, which is exactly the strong form — false at `β>4.30661` in the
one-plaquette and at the 4-D Savvidy (`ess inf<0`); it is H1. (b) *Reduces to H1*: every nonlocal
criterion that actually controls `λ₁` is the Cheeger/gate density. (c) *Orderless quotient*: higher
covariant derivatives `D_A^k` are gauge-covariant and stateable on `𝒜/𝒢`, but the working nonlocal
functional (1-D Cheeger) uses the **total order** of `ℝ` (half-lines are the optimal sets); `𝒜/𝒢`
carries no order (`T-WAVE` T3), so on the quotient only the tautology `λ₁≥h²/4` survives. No QED.

---

## §1 · What "higher-order / nonlocal Bakry–Émery" can mean — and the formal candidates

Let `dμ=Z^{-1}e^{-V}dVol`, reversible generator `L=Δ−⟨∇V,∇⟩`, `Γ(f)=|∇f|²`,
`Γ₂(f)=|Hess f|²+(Ric+Hess V)(∇f,∇f)`. The `k=0` marginal is 1-D, `V=V_eff`, and (H6-B)
`λ₁(ν)≥ess inf V_eff″ = ess inf G″`.

**Candidate (I): global iterated `Γₖ` / `CD(ρ,∞)`, `CD(ρ,n)`.** `CD(ρ,∞)⟺Ric+Hess V⪰ρg`;
`CD(ρ,n)` adds `+(1/n)(Lf)²`. At the Savvidy barrier `Ric+Hess S_W≈−2.3<0` (Q4 §0·2/§5·G) on the
negative mode, so `Γ₂(f,f)<0` for `f` linear there and **no** `ρ>0` version holds; `Γ₃,Γ₄` inherit
the negative leading curvature (they differentiate the same `Ric+Hess V`). **[NG]** for the class.

**Candidate (II): `k`-convexity / higher-order log-concavity.** `V` is `k`-convex if `V^{(k)}≥0`.
`2`-convexity (`V″≥0`, log-concavity) **does** give a Poincaré inequality (Brascamp–Lieb, Bakry–Émery,
standard `[T]`). The higher orders `k≥4` do **not** — §5·A/B are exact counterexamples. **[NG]**.

**Candidate (III): regularized curvature `κ_α=V_eff″+αV_eff⁗` (or `+β(V‴)²/…`).** Designed to turn
`G″<0` at the barrier into `κ_α>0` using the positive fourth derivative. Shown invalid both ways:
satisfied on a gapless well (§5·B) and over-claiming on a gapped one (§5·D). **[NG]**.

**Candidate (IV): nonlocal Dirichlet form / fractional LSI / energy functionals.** The LSI constant
for a nonlocal form is controlled by the "convexity in the small" of the interaction; on a double-well
energy this is again `~e^{-ΔV}` and the operative quantity is the **gate density** (§5·F). Every such
criterion that can certify the gap is a gate-density bound. **[NG]** for "new"; **`[O]=H1`** for the
positive content.

**Candidate (V): an integrated local functional.** Any criterion that *integrates* the jet over `ℝ`
(e.g. `∫(V″)_+ dx`, the barrier width) is nonlocal and sees `w`; but what it computes is the barrier
action / gate, i.e. Candidate (IV). So there is no third category.

---

## §2 · Exact tests

### 2.1 The fixed-jet double well — the decisive local-datum failure `[T]`+`[V]` (§5·A)

`V_w(y)=(w²/2)((y/w)²−1)² = y⁴/(2w²)−y²+w²/2`. Exact jets:

```
V′(y)=2y³/w² − 2y ,   V″(y)=6y²/w² − 2 ,   V‴(y)=12y/w² ,   V⁗(y)=12/w² ,   V^{(k)}=0 (k≥6).
critical points: y=0 (barrier, V″=−2) ,  y=±w (wells, V″=+4)   — independent of w.
as w→∞:  V^{(k)}(cp)→0  for k≥3 ;  barrier V(0)−V(±w)=w²/2→∞ ;  λ₁(ν_w) ~ e^{-w²/2}→0 .
```

The **whole critical-point jet** (`V′,V″,V‴,V⁗,…` at `0,±w`) converges to a fixed limit while the
gap collapses. Hence **no** function `f` of the critical-point jets — of any finite or fixed infinite
order — can lower-bound `λ₁`; `f` is constant along the family while `λ₁→0`. This is the sharpest
statement of the failure: the essential-infimum curvature `−2` and the well curvature `+4` are
**fixed positive-looking data** and still give no gap. (The `[NG]` is *not* about `ess inf V″` being
negative — here it is `−2` and fixed — it is that the gap is a nonlocal barrier action.)

### 2.2 The `4`-convex family — "higher-order log-concavity" is satisfied and useless `[NG]`+`[V]` (§5·B)

`V_g=(g/2)(x²−1)²`, `V⁗=6g>0` for **all** `x` (4-convex, `k`-convex for every even `k`), and
`κ_α(bar)=V″(0)+αV⁗(0)=(6α−1)g>0` for `α>1/6`. Yet the gap Kramers-collapses:
`λ₁=0.79,0.85,0.75,0.45,0.12,4.6×10⁻³` for `g=0.5,1,2,4,8,16`. The **4-convexity and the
regularized curvature grow without bound** while `λ₁→0`. So `{V⁗≥0}` (or `{κ_α≥ρ>0}`) is **not** a
sufficient condition for a Poincaré inequality, and cannot be the "higher-order log-concavity" the
proposal wants. This is consistent with the known sharp result that only `2`-convexity suffices.

### 2.3 The exact one-plaquette — `G″=A−Bβ²` and the regularized curvature `[T]`+`[V]` (§5·D)

Exact `θ=π`, `e₃` slice marginal (Gauss–Legendre, S2/G4): `A=0.20511530`, `B=0.01105929`,
`β*=√(A/B)=4.30661`; the witness reproduces `G″(0)=+0.13599,−0.07137,−0.90081` at `β=2.5,5,10`,
matching `A−Bβ²` to all digits.

```
 β    G″(0)      G⁗(0)     κ₁=G″+G⁗    λ₁(ν)    h²/4     gate
 2.5  +0.13599   +0.0964   +0.2324     0.3398   0.0944   −1.92
 5.0  −0.07137   +0.3312   +0.2598     0.4139   0.1970   −2.50
10.0  −0.90081   +1.5205   +0.6197     0.4042   0.4022   −2.82
```

Two facts. (i) **Below `β*`** `G″(0)>0`: log-concave, plain BE works. **Above `β*`** `G″(0)<0` but
`G⁗(0)>0`, so a higher-order regularization could *formally* keep `κ₁>0`. (ii) But the **exact
marginal gap stays `O(0.4)`** across `β=2.5→10` — because the one-plaquette marginal is a **compact
tilted single well** and the `G″<0` point is an **off-mass shoulder**, not a bottleneck (Q3's
closure: "`G″<0` is a log-convex shoulder, not a symmetric double well"). So the one-plaquette
**does not exhibit a Kramers collapse**; compactness keeps `χ` finite. Consequently (a) the plain
BE lower bound `ess inf G″` goes **vacuous** at `β>β*` (`−0.90<0`) yet the true gap is `0.40` — BE is
*lossy*, not wrong; and (b) the regularized `κ₁=0.62>λ₁=0.40` at `β=10` **over-claims**, so `κ₁` is
**not a valid lower bound**. The one-plaquette is therefore a *control*: it shows the negative
`G″(0)` set is real but not dynamically binding, and it **falsifies** the naive `λ₁≥G″+αG⁗` criterion.
The genuine threat is the **Savvidy cell** of §2.1, where the barrier genuinely separates two wells.

### 2.4 Savvidy cell `[T]`+`[V]`

The Savvidy operator `−D²+2ad F_B` has the negative mode `−g|B|`; its 1-D normal form is the double
well `V_eff(y)=½m²y²+¼κy⁴`, `m²<0`, wells `±√(−m²/κ)`, barrier `ΔV=m⁴/(4κ)` (S1 §4.2). Its
critical-point jet is `V″(well)=−2m²>0`, `V″(barrier)=m²<0`, `V⁗=6κ>0` — *exactly* the §2.2
structure. §2.1 is the cleanest possible Savvidy cell (barrier curvature held fixed, barrier grown by
width); §2.2 is the `κ`-driven cell (barrier grown by depth). Both give `λ₁→0` at bounded
curvature-type data. **The Savvidy cell is a Kramers double well; no jet criterion sees its gap.**

### 2.5 Acid tests `[T]`+`[V]` (§5·C, §5·E)

- **Free massless field.** `V=½k²x²`: `V^{(k)}=0` for all `k≥3`, `λ₁=k²→0`. Every higher-order
  criterion that is a jet function returns `0` (its value is `f(k²,0,0,…)→0`). **ACID PASS** — no
  false floor.
- **1-D `SU(2)`.** Exact one-link Haar marginal `ρ₁(z)=π(Cin2π−Cin|z|)` is log-concave:
  `G″(0)=0.205115`, `Var=3.7198`, `λ₁=0.264451`, `χ≤3.781`; `β=2.5` periodic floor `G″→0.222`,
  `Var≈4.45`, `χ≤4.5`. **Finite `χ` reproduced** — but by *plain* BE (single well, `ess inf G″>0`);
  the higher-order machinery is unnecessary and gives no different number. This is the honest sense
  in which the 1-D test "passes": the criterion is non-vacuous where BE already is, and silent where
  BE fails.

---

## §3 · The three kill criteria, explicitly

**(a) "Requires the defect absent" ⟹ circular.** The natural higher-order repair is to demand
`ess inf G″≥ρ>0` everywhere (no negative-curvature point). This is exactly the strong form killed in
V3/G4: the one-plaquette has `G″(0)=A−Bβ²<0` for `β>4.30661`, and the 4-D Savvidy has
`ess inf G″<0` at `β=2.5`. So the hypothesis is unavailable exactly where it is needed; demanding it
*is* demanding the 4-D gap. **Fires.** [NG]

**(b) "Reduces to H1" ⟹ no new object.** Any nonlocal functional that actually certifies
`λ₁≥ρ` on the double well must read the barrier action; the minimal such object is the Cheeger gate
`h(ν)`, and `h≥h₀⟺ρ(gate)≥c⟺H1` (§5·F). Integrated local jet functionals compute the same barrier
action (Candidate V). **Fires.** [NG] for "new", [O] for the content.

**(c) "Cannot be stated on the orderless quotient" ⟹ partial.** The higher **covariant** derivatives
`D_A^k` are gauge-covariant and transport to `𝒜/𝒢`; so the *local* higher-order criterion **can** be
stated there (this is not the obstruction — its failure is §2.1/§2.2). But the only functional that
*works* is the 1-D Cheeger constant, whose exactness uses the **total order** of `ℝ` (half-lines are
the optimal isoperimetric sets). On `𝒜/𝒢` there is no order (`T-WAVE` T3: `TP₂` cannot be stated on
the orderless quotient), so the exact gate functional has no canonical lift; only the tautology
`λ₁≥h²/4` survives, with `h` un-defined on the quotient. **Fires** for the nonlocal branch. [NG]

So the proposal fails every branch: the local/higher-jet branch is **false** (§2.1–2.2), the
nonlocal branch is **H1** or **unstatable** (§3b–c), and the "repair the sign" branch is **circular**
(§3a).

---

## §4 · The exact missing step — unchanged `[O]`

Nothing in this wave moves the wall. What is missing for 4-D at `β=2.5` is the same single object as
H6/Q4, now with one more equivalent face:

```
∃ c>0, p<1/2 :  ν({G″<0}) ≤ p   AND   ρ(gate) ≥ c   AND   {G″≥ρ_good} connected ,   uniformly (L,a).
                                                                                    [O] = H1
```

Equivalently: a **uniform lower bound on the `k=0` marginal density at the Savvidy gate**, i.e. an
upper bound on the barrier free-energy `V_eff(barrier)` relative to the well. Y2 adds the precise
reason no *local* higher-order criterion can substitute: the barrier action is nonlocal (§2.1), and
the regularized curvature is not a valid BE constant (§2.2–2.3). The higher-order/nonlocal program
therefore does not lower the wall; it confirms the wall is a **gate-density** statement, not a
curvature statement.

---

## §5 · Witness — `research/v3/w7/y2_higher_be.py` (`py -3.12 -I`, deterministic, numpy+scipy)

Quoted output (abridged; full in `y2_higher_be.out`):

```
(A) FIXED-JET DOUBLE WELL V_w=(w^2/2)((y/w)^2-1)^2 : critical-point jets fixed, lambda_1 -> 0
      w    V''(barrier)  V''(well)  V''''(barrier)  barrier  lambda_1(FEM)  e^{-barrier}
     1.0     -2.0000     +4.0064     1.2000e+01      0.500    8.997e-01      6.07e-01
     1.5     -2.0000     +3.9936     5.3333e+00      1.125    3.188e-01      3.25e-01
     2.0     -2.0000     +4.0064     3.0000e+00      2.000    1.137e-01      1.35e-01
     3.0     -2.0000     +3.9936     1.3333e+00      4.500    9.089e-03      1.11e-02
     4.0     -2.0000     +4.0064     7.5000e-01      8.000    2.859e-04      3.35e-04
(B) 4-CONVEX V_g=(g/2)(x^2-1)^2 : V''''=6g>0, kappa_a(bar)=(6a-1)g>0 for a>1/6, but lambda_1->0
       g    V''''=6g  V''(barrier)  kappa_1/3  kappa_1   lambda_1
     0.5     3.000      -0.500       +0.500    +2.500   7.921e-01
     2.0    12.000      -2.000       +2.000   +10.000   7.487e-01
     8.0    48.000      -8.000       +8.000   +40.000   1.191e-01
    16.0    96.000     -16.000      +16.000   +80.000   4.575e-03
(C) FREE FIELD V=1/2 k^2 x^2 : V''''=0, lambda_1=k^2 -> 0   (ACID PASS)
     k=1.00 lambda_1=0.99999 | 0.50 0.25000 | 0.10 0.01000 | 0.02 0.00040
(D) EXACT ONE-PLAQUETTE theta=pi, G''(0)=A-B beta^2 (A=0.20511530 B=0.01105929 beta*=4.30661)
     beta   G''(0)     A-B b^2    G''''(0)   kappa1=G''+G''''  lambda_1   h^2/4   gate
     2.50  +0.13599   +0.13599    +0.0964     +0.2324        0.3398   0.0944  -1.92
     5.00  -0.07137   -0.07137    +0.3312     +0.2598        0.4139   0.1970  -2.50
    10.00  -0.90081   -0.90081    +1.5205     +0.6197        0.4042   0.4022  -2.82
     => kappa1=+0.62 > lambda_1=0.40 above beta*: not a lower bound [NG].  Marginal gap stays O(0.4)
        (tilted single well, G''<0 is an off-mass shoulder, Q3) — one-plaquette has NO Kramers collapse.
(E) 1D SU(2) exact Haar: G''(0)=0.205115  Var=3.7198  lambda_1=0.264451  => chi<=3.781
     (beta=2.5 periodic floor G''->0.222, Var~4.45 => chi<=4.5)  [single well: plain BE suffices]
(F) NONLOCAL Cheeger gate h(nu): h>=h0 <=> rho(gate)>=c <=> H1
       g    rho(gate)    Cheeger h    h^2/4      lambda_1
     0.5   2.561e-01    5.122e-01   6.56e-02   7.921e-01
     2.0   1.864e-01    3.728e-01   3.47e-02   7.487e-01
     8.0   1.932e-02    3.865e-02   3.73e-04   1.191e-01
    16.0   5.207e-04    1.041e-03   2.71e-07   4.575e-03
```

Readings: **(A)** the full critical-point jet is `w`-independent while `λ₁≈e^{-barrier}` collapses —
no local jet criterion can give a uniform floor; **(B)** 4-convexity and `κ₁=5g>0` coexist with
`λ₁→0` — higher-order log-concavity is not sufficient; **(C)** free field returns `0` (acid pass);
**(D)** exact `A−Bβ²` reproduced, the regularized curvature over-claims (`κ₁>λ₁`), and the
one-plaquette marginal is a gapped single well (no Kramers collapse) — the wall is the Savvidy cell,
not the one-plaquette; **(E)** 1-D `SU(2)` finite `χ` via plain BE; **(F)** the only working nonlocal
functional is the gate density `=H1`.

---

## §6 · Status ledger

| item | statement | label |
|---|---|---|
| iterated `Γₖ` / `CD(ρ,n)` at the Savvidy barrier | `Ric+Hess V≈−2.3<0` on the negative mode ⇒ `Γ₂<0`; no `ρ>0` version | **[NG]** (Q4/H6) |
| `2`-convexity ⇒ PI | Brascamp–Lieb / Bakry–Émery (log-concave) | **[T]**, standard |
| `4`-convexity ⇒ PI? | **false**: `V_g` has `V⁗=6g>0`, `κ₁=5g>0`, `λ₁→0` | **[NG]**+**[V]** (§5·B) |
| fixed-jet double well | critical-point jet `w`-independent, `λ₁~e^{-w²/2}→0` | **[T]**+**[V]** (§5·A) |
| regularized curvature `κ_α=V″+αV⁗` as BE constant | not a lower bound (over-claims `κ₁>λ₁` on exact one-plaquette) | **[NG]**+**[V]** (§5·D) |
| exact one-plaquette `G″=A−Bβ²` | `A=0.20511530`, `B=0.01105929`, `β*=4.30661` reproduced | **[T]**+**[V]** (§5·D) |
| one-plaquette marginal above `β*` | `G″(0)<0` but compact tilted single well, `λ₁≈0.34–0.41` (no Kramers collapse) | **[T]**+**[V]** (§5·D) |
| Cheeger `λ₁≥h²/4`, `h=inf ρ/min(F,1−F)` | verified `h²/4≤λ₁` on one-plaquette at all `β` | **[T]**+**[V]** |
| nonlocal gate `h~ρ(gate)~e^{-ΔV}` | `h≥h₀⟺ρ(gate)≥c⟺` no interior barrier `⟺H1` | **[T]**/**[NG]** (§5·F) |
| free massless field | every higher jet `0`, `λ₁=k²→0` | **[T]**+**[V]** (§5·C) |
| 1-D `SU(2)` Haar | log-concave, `λ₁=0.2645`, `χ≤3.78`; `β=2.5`: `χ≤4.5` | **[T]**+**[V]** (§5·E) |
| integral jet functional | computes the barrier action ⇒ same as nonlocal gate ⇒ H1 | **[T]** (§1·V) |
| orderless quotient `𝒜/𝒢` | `D_A^k` stateable, but 1-D Cheeger needs the total order; only tautology survives | **[NG]** (§3c) |
| **uniform `ρ(gate)≥c` / `ess inf G″≥λ`** | the surviving open step | **[O] = H1** |
| **higher-order / nonlocal BE closes the gap** | **no** (all branches fail / circular / H1) | **[NG]**, 0% |

---

## §7 · Sharpest line

**Higher-order and nonlocal Bakry–Émery do not rescue the Savvidy soft mode. The gap of a double well
is a nonlocal barrier action `e^{-ΔV}`, invisible to every finite-order jet: the exact family
`V_w=(w²/2)((y/w)²−1)²` holds its full critical-point jet fixed (`V″=−2` at the barrier, `+4` at the
wells) while `λ₁→0`; the `4`-convex `V_g` keeps `V⁗=6g>0` and `κ₁=5g>0` at a collapsing `λ₁`; and on
the exact one-plaquette (`G″=A−Bβ²`, `β*=4.30661`) the regularized curvature `G″+G⁗` *over-claims*
(`0.62>λ₁=0.40`), so it is not a BE constant. The only criterion that does control `λ₁` is the
Cheeger/gate density `ρ(gate)`, which *is* H1; the local-jet branch is false, the nonlocal branch is
H1 or unstatable on the orderless `𝒜/𝒢`, and the "keep the curvature positive" branch is circular
(it demands the defect absent). Free field returns `0` (acid pass); 1-D `SU(2)` reproduces finite `χ`
under *plain* BE. What is missing is unchanged — one number, a uniform lower bound on the marginal
density at the Savvidy gate. No QED. Clay 0%.**

---

*Refs (internal).* `Q4-BAKRY-EMERY.md` (§0, §1.3 marginal BE, §2 weighted/capacity); `H6-HEAD-SPECTRAL.md`
(§1 H6-A/B, §3 typical-point vs essential-infimum, §4 method-class `[NG]`); `T-WAVE-SYNTHESIS.md` (§1
predicates, T3 orderless quotient/`TP₂`); `S2-BOUNDARY-COUNTEREXAMPLE.md`/`G4-GENERAL-BOUND.md`
(`A−Bβ²`); `H3-AVERAGED-HEAD.md`/`H4-HEAD-COUNTEREXAMPLE.md` (exact one-plaquette, `β*=4.30661`);
`Q3-SCHUR-DEFECT.md` (one-plaquette `G″<0` is a log-convex shoulder, not a double well);
`R5-EXACT-ANCHORS.md`; `S1-CONDITIONAL-POINCARE.md` §4.2 (Savvidy normal form); `B3-FUNCTIONAL-BE`,
`OTHER-H1-POINCARÉ`, `CHEEGER-CORE` (marginal BE lossless, Cheeger tautology).
*External (standard, textbook tools, cited from corpus, not re-fetched):* Bakry–Émery 1985; Brascamp–Lieb,
Adv. Math. **22** (1976) 366; Cheeger, *Problems in Analysis* (1970) 195; Bobkov–Ledoux; Bovier–Eckhoff–
Gayrard–Klein (Kramers/metastability); Savvidy, Phys. Lett. B **71** (1977) 133. Witness:
`research/v3/w7/y2_higher_be.py`.
*No third-party sources downloaded this session; all external inputs are standard and cited from the
corpus, not re-fetched.*

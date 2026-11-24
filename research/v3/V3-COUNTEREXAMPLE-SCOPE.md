# V3-COUNTEREXAMPLE-SCOPE — the `G″<0` counterexample to the coarse-curvature inequality (★) is made rigorous, is **in scope for the true coarse marginal curvature (C)** and **out of scope for the literal C2 variance inequality (U)**, and it forces one specific correction: evaluate `G″` at the coarse minimiser `a₀*`, not uniformly in the coarse field `a₀`

> 2026-11-24 · Decides the scope of the G4 counterexample to the single live inequality `(★)`.
> Reads: `G4-GENERAL-BOUND.md` (the counterexample `G″=0.2057−0.0123β²`), `S2-BOUNDARY-COUNTEREXAMPLE.md`,
> `G1-BLOCK-IDENTITY.md`, `G3-PROJECTION-TERM.md`, `STAR-WAVE-SYNTHESIS.md`,
> `docs/paper-v2/ch05-theorems.tex` (`eq:split`, `eq:vardecomp`, `sec:wgeo`).
> Witness: `research/v3/w7/v3_counterexample_scope.py` (`py -3.12`, deterministic Gauss–Legendre,
> output quoted verbatim below). *(No `RECONCILIATION-STAR.md` exists; the closest is
> `research/v3/w7/r1_model_curvature.py`.)*
>
> **Labels (charter §2).** **[T]** true theorem (proved here) · **[C]** conditional · **[V]** numerical
> (deterministic quadrature; residual / n-independence quoted) · **[NG]** no-go for a named class ·
> **[O]** exact missing step.

---

## 0 · Verdict (no QED; Clay 0%)

**The counterexample is rigorous, and its scope is exactly the true coarse marginal curvature (C).**
Concretely, in the exactly solvable SU(2) one-plaquette / one-free-link model of S2:

1. **[T]+[V] Exact value of the true `G″` at the barrier.** At boundary `B→−I` (`θ=π`), longitudinal
   mode `v=e₃`, coarse field `a₀=0`, the marginal curvature is *exactly*
   ```
     G″(e₃; θ=π, a₀=0) = A − B β² ,     A = E_Haar[Φ′(r)/r] = 0.20511530 ,  B = Var_Haar[sin(r/2)/r] = 0.01105929 ,
     ⟹ β* = √(A/B) = 4.30661 ,   G″(5) = −0.0714 ,  G″(10) = −0.9008 ,  G″(20) = −4.2186 .           (0.1)
   ```
   Validated three independent ways to `≤10⁻⁶` (analytic conditioned moments, finite-difference of
   `−d²/da² log ρ`, and the conditioned decomposition `Δ_dyn+W_geo+X`). This **corrects** G4's quoted
   `A=0.2057, B=0.0123, β*=4.09`: those came from a missing Gauss–Legendre weight in
   `g4_general_bound.py`'s `Gpp_e3_analytic` (`w=exp(−h)·r` instead of `exp(−h)·r·w_GL`). The qualitative
   conclusion `G″<0` for `β>β*` is unchanged; the correct threshold is `β*=4.307`, not `4.09`.
2. **[V] The counterexample is longitudinal-only.** At the same `θ=π, a₀=0`: transverse `v=e₁` gives
   `G″=0.214, 0.362, 0.612, 1.118` at `β=0.5,2.5,5,10` (**positive, growing**); at `θ=0` (B=I) both modes give
   `G″=0.291,0.730,1.342,2.587` (**positive, manifold**). The Savvidy-style trouble is on the boundary
   axis, not the transverse plane — the opposite of the S2 target.
3. **[T]+[V] The sign is a coarse-field effect, not a boundary-only one.** On `θ=π, e₃`, the coarse
   effective action `S_eff(a₀)=−log ρ_{e₃}(a₀)` has its **single well (coarse minimiser) at `a₀*≈−2.1,
   −2.6, −2.9`** for `β=2.5,5,10`, and there `G″(mode)=+0.69, +1.34, +2.61>0` (grows with `β`). The
   negative-curvature set is the **concave shoulder** `a₀∈(−0.2, 3.8)` at `β=5`, which *contains the fixed
   point `a₀=0`* but *not the mode*. So `G″` is positive where the marginal concentrates and negative
   on an off-minimum shoulder.
4. **[NG] The strong form is dead; the weakest survivor is "at the coarse minimiser".** "`G″(v;a₀)≥δ_W>0`
   uniformly in all modes `v`, all boundaries `U′`, all coarse fields `a₀`, all `β` in the tower" is
   **false** (0.1). The weakest, physically-correct correction that excludes (0.1) is: **restate (★) as
   convexity of the coarse effective action at its own minimum**, `G″(v; a₀*(U′,v))≥δ_W(U′,v)>0`, i.e.
   evaluate at the coarse minimiser, not uniformly in `a₀`. This survives the exact model `[V]`, its
   *non-degeneracy* `G″(a₀*)≥0` is `[T]` (trivial), and its *uniform* lower bound `δ_W>0` is exactly the
   known wall **H1** `[O]`.
5. **[T] The dichotomy that decides "in scope".** (★) admits two inequivalent readings:
   - **`(★_U)`, the literal C2/S2 variance inequality** `Var(G)≤E[H]+W_j−δ` `⟺ G″_U:=E[H]−Var(G)+W_geo≥δ`.
     By the exact flat-chart identity this is `=2W_geo`, and S2/G1/G3 prove it positive. **The counterexample
     does NOT refute `(★_U)`** — it is a statement about a different object.
   - **`(★_C)`, the true coarse marginal curvature** `G″(v;a₀)=−d²/da² log ρ_v(a)|_{a₀}`, which is what
     `ch05` `eq:split`/`eq:vardecomp` **actually evaluate** (the `|0` conditioning made explicit) and what
     the coarse effective action / gap needs. **The counterexample refutes the strong form of `(★_C)`.**
   So (★) is *either* true-but-irrelevant (`(★_U)`) *or* false-in-the-strong-form (`(★_C)`). That is the
   rigorous finding; it does not close the program, it fixes the object and the evaluation point.

**Acid test passes.** Free massless Gaussian (`Φ=0`, `V=κx²/2`): `G″=κ→0` as `κ→0` (no manufactured floor);
`β=0` pure Haar: `G″=A=0.2051>N/12=0.1667`, matching `ch05`'s β=0 theorem `G″≥N/12` `[T]+[V]`.

---

## 1 · The model and the exact object — **[T]**

SU(2), one plaquette `P=UB`, one free link `U=exp(ξ)`, `ξ∈ℝ³` with `ρ=|ξ|∈[0,2π]`, boundary
`B=exp(iθσ₃/2)`:
```
  S(U;B)=β(1−½ tr(UB)) = β(1 − cos(ρ/2)cos(θ/2) + sin(ρ/2)sin(θ/2) u) ,   u := ξ₃/ρ ,
  dU = J d³ξ ,  J ∝ sin²(ρ/2)/ρ² ,  Φ := −log J = 2log(ρ/2) − 2log sin(ρ/2) ,  dμ(·|B) ∝ e^{−S−Φ} d³ξ .
```
(This is S2's model. `ρ=2π` is the cut locus where `e^{−Φ}=J→0`.)

**The object.** For a unit coarse direction `v` and coarse coordinate `a=⟨v,ξ⟩`, the coarse marginal
`ρ_v(a)=∫δ(⟨v,ξ⟩−a)e^{−V}d³ξ` (`V=S+Φ`) has curvature, by Lemma G1 (`[T]`, pure differentiation):
```
  G″(v;a₀) := −d²/da² log ρ_v(a)|_{a₀} = E_{μ}[∂_v²V | a₀] − Var_{μ}(∂_vV | a₀)      (marginal identity)
            = Δ_dyn(v;a₀) + W_geo(v;a₀) + X(v;a₀)   (ch05 eq:split, conditioned at a₀)          (1.1)
  Δ_dyn=E[H]−Var(G),  W_geo=E[Φ″]−Var(Φ′),  X=−2Cov(G,Φ′),  G=∂_vS, H=∂_v²S, Φ′=∂_vΦ, Φ″=∂_v²Φ.
```
All expectations/Var are **conditioned** on `a=a₀` (this is the `|0` of `ch05 eq:split` made explicit;
it is what makes the Var non-zero and the object the coarse curvature rather than a deterministic
fine-field Hessian). This `(C)` is **not** the unconditioned `(U):=E[H]−Var(G)+W_geo` that the flat-chart
Stein identity turns into `2W_geo`. The two coincide only when `Var(∂_vV|a₀)` and the `E[H|a₀]` compensation
match, which they do *not* at a barrier.

---

## 2 · Exact computation of the true `G″` — **[T]** formula, **[V]** constants

### 2.1 The barrier, longitudinal `e₃`, `θ=π`, `a₀=0` — exact formula

At `θ=π`: `cos(θ/2)=0, sin(θ/2)=1`, so `S=β(1+sin(ρ/2)u)`. On the slice `ξ₃=a₀=0`: `u=0`, `ρ=r=|(ξ₁,ξ₂)|`, hence
`V=β+Φ(r)` and the slice measure is **β-independent pure Haar** `∝ sin²(r/2)/r dr`. The Cartesian
derivatives along `e₃` are
```
  ∂_{e₃}V|₀ = g(r)/r = β sin(r/2)/r ,      ∂²_{e₃}V|₀ = h′(r)/r = Φ′(r)/r      (h′=Φ′ since the β-part of h is constant at θ=π).
```
Therefore
```
  G″(e₃;θ=π,a₀=0) = E_w[Φ′(r)/r] − β² Var_w[sin(r/2)/r] ,   w(r) ∝ sin²(r/2)/r ,  r∈(0,2π).       (2.1)
```
This is *exact*; `A := E_w[Φ′(r)/r] = 0.20511530`, `B := Var_w[sin(r/2)/r] = 0.01105929`, both converged to
10 digits for `n=800,1500,3000,6000,12000` GL nodes. **Decomposition:** `X=0` exactly on this slice
(`∂_{e₃}Φ=Φ′(ρ)u∝u=0`), `W_geo=E[Φ′/r]=A>0`, `Δ_dyn=−β²B<0` (the Savvidy force variance). Witness table:

```
[1] theta=pi, v=e3, a0=0 :  G'' = A - B beta^2
      beta   A-Bβ²      fd(−d²logρ)  decomp direct    Δ_dyn      W_geo       X
      0.00   0.205115     0.205115     0.205115      −0.00000   0.20512   −0.0e+00
      1.00   0.194056     0.194056     0.194056      −0.01106   0.20512   −0.0e+00
      2.50   0.135995     0.135995     0.135995      −0.06912   0.20512   −0.0e+00
      4.00   0.028167     0.028167     0.028167      −0.17695   0.20512   −0.0e+00
      4.31   0.000001     0.000001     0.000001      −0.20511   0.20512   −0.0e+00
      5.00  −0.071367    −0.071367    −0.071367      −0.27648   0.20512   −0.0e+00
     10.00  −0.900814    −0.900814    −0.900814      −1.10593   0.20512   −0.0e+00
     20.00  −4.218602    −4.218602    −4.218602      −4.42372   0.20512   −0.0e+00
```
Three methods (exact formula, 5-point FD of `−log ρ`, conditioned decomposition) agree to `≤10⁻⁶`. This is
`[T]` for the formula (2.1) and `[V]` for `A,B`. **Mechanism (the whole point of "including X"):** in the
*unconditioned* measure `E[H]` and `Var(G)` are both positive and cancel to give `2W_geo>0` (S2);
*conditioning at `a₀=0` makes `E[H|0]=∂²_{e₃}S|₀ = 0`* (the action is flat in the `e₃`-Hessian at the well
barrier at `θ=π`), removing the compensator and exposing `−Var(G|0)` — so the unconditioned positivity `2W_geo`
is precisely the wrong bookkeeping at a barrier. `X=0` here; at `a₀≠0` all three terms are nonzero and the
direct moment `E[∂²V]−Var(∂V]` (= (1.1)) includes `X` automatically.

### 2.2 Data correction to G4 — **[V]**

G4's quoted `A=0.2057, B=0.0123, β*=4.09` do **not** reproduce with the correct quadrature weight.
`g4_general_bound.py`'s `Gpp_e3_analytic` sets `w = np.exp(−h)*r` and omits the Gauss–Legendre weight
`WRHO`. Replicating that (buggy) sum reproduces `0.2057 / 0.0123 / 4.090` exactly; the weight-included,
n-converged values are `0.20511530 / 0.01105929 / 4.30661`. Neither affects any qualitative conclusion;
both give `G″<0` for `β≳4.1–4.3`. **Scope verdict uses the corrected values.**

### 2.3 Transverse `e₁` and `θ=0` — **[V]** (both positive)

```
[2] theta=pi, v=e1 (TRANSVERSE), a0=0 :    Gpp(e1) = fd    (Var(d1V|0)=0, analytic=FD to 1e-6)
      beta=0.50 0.213594 ;  1.00 0.237435 ;  2.50 0.362312 ;  5.00 0.612218 ;  10.00 1.118117   (all > 0)
[3] theta=0 (B=I), a0=0 :  e3=e1 by isotropy;  Gpp = 0.290877, 0.729924, 1.342360, 2.587494  (β=.5,2.5,5,10)
```
At `θ=0`, `S=β(1−cos(ρ/2))` is radial, `g≡0`, `∂_vV|₀=0`, and `G″=E_w[h′(r)/r]` with
`h′=β/2·sin(r/2)+Φ′(r)>0` on `(0,2π)` — **manifestly positive for every β** `[T]`. So the longitudinal
tachyonic direction exists **only at maximal holonomy `θ=π`**, and only off the coarse minimum.

### 2.4 Well anatomy: mode vs shoulder — **[V]**

At `θ=π, e₃`, `S_eff(a₀)=−log ρ_{e₃}(a₀)` has a single well at `a₀*≈−2.1,−2.6,−2.9` (β=2.5,5,10) and a
concave shoulder where `G″<0`:
```
[4]  beta=2.5: mode a0*=−2.125  G''(mode)=+0.6903   G''(a0=0)=+0.1360
     beta=5.0: mode a0*=−2.625  G''(mode)=+1.3414   G''(a0=0)=−0.0714
     beta=10.0 mode a0*=−2.875  G''(mode)=+2.6057   G''(a0=0)=−0.9008
```
At `β=5` the negative set is `a₀∈(−0.2, 3.8)` (G″ crosses zero at `a₀≈−0.2` and again at `≈3.8`). Thus
`a₀=0` (the point `ch05` writes as `G″(0)`) sits in the concave shoulder for `β>4.31`, while the coarse
minimiser `a₀*` — where the marginal concentrates — always has `G″>0`, growing with `β`.

---

## 3 · Scope analysis — is the counterexample in the admissible class of (★)? — **[T]

The admissible class of `(★_C)` is: "for every shell `j`, every boundary `U′`, every coarse direction `v`,
and every coarse field `a₀`." We test the three sub-questions.

1. **Does the projection `P` admit longitudinal modes? — YES.** `P` is a coarse-graining whose image
   spans the adjoint algebra `su(N)`; `w=Pv=e₃` (the boundary/holonomy axis) is a perfectly admissible
   *fixed* fine direction. Nothing in the definition of `P` or in `(★_C)` excludes the longitudinal mode.
   (The S2/G1/G3 result "`P` harmless, `G″=2W_geo`" is about object `(U)` in the *unconditioned* measure,
   not `(C)`; it is not a restriction on which `v` `(★_C)` quantifies over.) **In scope.**

2. **Is the barrier / shoulder configuration in the support of the conditional measure? — YES.** The
   conditional density `e^{−V}d³ξ/ρ_v` is strictly positive on the whole open fibre (the `ρ=2π` cut locus
   is the only edge, where it vanishes), so *every* interior coarse field `a₀∈(−2π,2π)` — including
   `a₀=0` and the shoulder — has positive probability. `a₀=0` is admissible (it is not the mode at `θ=π`,
   but it is a point of full support). **In scope.**

3. **Is the boundary `B→−I` (`θ=π`) admissible? — YES (at the edge of the holonomy range).** The boundary
   links are *conditioned*, so `(★_C)` is defined for every fixed `U′`; and in the RG application `U′` is
   drawn from the equilibrium boundary measure, which on compact `SU(2)` is absolutely continuous w.r.t.
   Haar and has full support, so `θ=π` is in support (an endpoint of `[0,π]`, not excluded). **In scope for
   the uniform-in-`U′` reading** (which STAR-WAVE already flags as too strong) ; `θ=π` is *not* in the
   small-field set `𝒮_t`, so `(★_sf)` excludes it.

**Verdict on scope.** The counterexample is **in the admissible class of `(★_C)` as stated (uniform in
`U′,v,a₀`)** and **rigorously refutes that strong form** `[NG]`. It does **not** touch the literal C2
variance inequality `(★_U)` (a different object). The *weakest restriction that excludes it* is one of,
in increasing strength:
```
  (i)   restrict the coarse field to the coarse minimiser a₀* (NOT uniform in a₀)   <— weakest, physical;
  (ii)  restrict to the small-field boundary set 𝒮_t  (excludes θ=π);  or transverse modes only;
  (iii) allow boundary- AND field-dependent δ_W(U′,v,a₀) (alone insufficient: at a₀=0, δ must go negative).
```
(i) is the weakest and the physically meaningful one: convexity of the coarse effective action is needed
*at its own minimum*, not at every point of the fibre. (ii) is Balaban's excision; (iii) by itself does not
rescue the uniform-in-`a₀` form.

---

## 4 · The no-go and the weakest corrected inequality

### 4.1 No-go — **[NG]**

> **[NG] "Strong-form `(★_C)`."** There is **no** `δ_W>0` with `G″(v;a₀)≥δ_W` for all coarse fields
> `a₀` (equivalently: the coarse effective action is **not** uniformly `δ_W`-strongly convex), uniformly in
> `U′` and `v`. Refutation: one plaquette, `B→−I` (`θ=π`), `v=e₃` (longitudinal), `a₀=0`, `β>β*=4.30661`,
> where `G″=A−Bβ²<0` exactly (2.1), while the unconditioned S2 object is `2W_geo=2A=0.410>0`. **[T]** (model)
> + **[V]** (constants). By `θ=π` being full-support and `a₀=0` being full-support, this point is admissible.

The failure is *the coarse image of the Savvidy/S1 double-well*: conditioning the coarse field at a shoulder
sets `E[H|a₀]=0`, exposing the negative Wilson-force variance `−Var(G|a₀)=−Bβ²`, while the Haar curvature
`W_geo=A` stays fixed; their crossing at `β*` is the marginal-curvature analogue of the Savvidy negative mode.
It is **the same wall as S1's `[NG]` ("uniform log-concavity fails on chromomagnetic boundaries")**, now seen
on the true coarse object.

### 4.2 Weakest surviving inequality — **[T]** non-degeneracy, **[V]** positivity, **[O] uniform `δ`**

> **(★_loc) — convexity at the coarse minimiser.** For every boundary `U′`, mode `v`, let
> `a₀*(U′,v) = argmax ρ_v` (the coarse minimiser / mode of the marginal). Then
> `G″(v; a₀*) = −d²/da² log ρ_v|_{a₀*} ≥ 0`, with strict positivity iff the coarse well is non-degenerate,
> and one asks for a uniform lower bound `G″(v;a₀*) ≥ δ_W(U′,v) > 0`.

- **`G″(a₀*)≥0` is `[T]` and trivial:** `a₀*` is a local maximum of the smooth density `ρ_v`, so
  `(log ρ_v)″≤0`, i.e. `G″=−(log ρ_v)″≥0`. There is no positivity to manufacture — a genuine coarse well has
  positive curvature at its bottom.
- **`G″(a₀*)>0` strictly** holds on the exact model scan: `+0.69, +1.34, +2.61` at `β=2.5,5,10` (§2.4), and
  `[V]` the value grows with `β` (the Wilson well sharpens); at `θ=0` it is `>0` for all `β` `[T]`. So
  **(★_loc) survives the counterexample.**
- **A uniform `δ_W>0` in `(★_loc)`** (uniform in `U′,v`, per shell) is **exactly H1 / the coarse spectral
  gap** `[O]`: `G″(a₀*)` is the inverse local variance of the coarse field (`=1/Var_ρ(a)` in the Gaussian
  approximation), so a uniform lower bound is a uniform coarse-mass bound. This is the known wall, not a new
  one. Equivalently it is the `S≥δ` factor of S3/STAR-WAVE `[O]`.

So the counterexample *does not force a new machine*; it forces the corrected **object** (`(★_C)`) and the
corrected **evaluation point** (`a₀*`, not uniform in `a₀`), at which the surviving statement is `[T]`/`[V]`
and its hard part is the pre-existing H1 wall. What it *does* prove is that the sharp "uniform-in-`a₀`" strong
form is dead and that `2W_geo` (object `(U)`) is not a certificate of coarse convexity.

---

## 5 · Acid tests — **[T]+[V]**

1. **Free massless field.** `Φ=0`, `V=κx²/2`: `G″=κ→0` as `κ→0` (witness `[5]`). The tool manufactures no
   floor on the free theory; a positive `δ_W` requires `β>0` (the Wilson interaction). **Required, holds.**
2. **`β=0` pure Haar.** `G″(e₃;θ=π,a₀=0)=A=0.205115 > N/12=0.166667` (`[6]`), matching `ch05`'s β=0 theorem
   `G″≥N/12` `[T]`, while the unconditioned `2W_geo→0` (S2). No contradiction: the β=0 spectrum is governed
   by the Haar floor `A`, not by `2W_geo`.
3. **Continuity / sign law.** `G″(e₃;θ=π,a₀=0)=A−Bβ²` is continuous in `β`, equals the Haar floor at `β=0`,
   and crosses zero at `β*=4.30661`; the positivity at small `β` is the Haar geometric curvature, destroyed at
   large `β` by the Wilson-force variance `−Bβ²`. `[T]`+`[V]`.

---

## 6 · Status ledger

| item | statement | label |
|---|---|---|
| marginal identity `G″(v;a₀)=−d²log ρ_v/da²\|_{a₀}=E[∂²V\|a₀]−Var(∂V\|a₀)` | Lemma G1, any block, exact | **[T]** |
| exact barrier formula `G″(e₃;θ=π,a₀=0)=A−Bβ²`, `A=0.20511530`, `B=0.01105929` | (2.1); 3 methods to `≤1e-6`; n-independent | **[T]**+**[V]** |
| `β*=√(A/B)=4.30661`; `G″(5)=−0.0714`, `G″(10)=−0.9008`, `G″(20)=−4.2186` | exact + FD + decomposition | **[T]**+**[V]** |
| G4's `A=0.2057, B=0.0123, β*=4.09` | missing-GL-weight artifact (`w=exp(−h)r`); corrected above | **[V]** (correction) |
| mechanism: conditioning `a₀=0` makes `E[H\|a₀]=0`, exposes `−Var(G\|a₀)`; unconditioned `E[H]` cancels it | (2.1)/(§2.1) | **[T]** |
| transverse `e₁` at `θ=π`, and both modes at `θ=0`, stay positive for all `β` | witness tables [2],[3]; `θ=0` manifest `[T]` | **[V]**+**[T]** |
| coarse well at `a₀*≈−2.1/−2.6/−2.9`; `G″(a₀*)>0` and grows with `β`; negative set is shoulder `(−0.2,3.8)` at `β=5` | witness [4], a0-scan | **[V]** |
| longitudinal `e₃` is admissible under `P`; `a₀=0` and `B→−I` are full-support | scope §3.1–3.3 | **[T]** |
| **strong-form** `G″(v;a₀)≥δ_W` uniform in `(U′,v,a₀,β)` | **false** (counterexample at `θ=π,e₃,a₀=0,β>4.31`) | **[NG]** |
| literal C2 variance inequality `(★_U): Var(G)≤E[H]+W_j−δ ⟺ 2W_geo≥δ` | true (S2/G1/G3) but **wrong object** for coarse curvature; not refuted here | **[T]** (but insufficient) |
| **(★_loc)**: `G″(v;a₀*)≥0` at the coarse minimiser | `[T]` trivial (max of `ρ`); strict `[V]`; survives counterexample | **[T]**+**[V]** |
| uniform `δ_W>0` in `(★_loc)` (uniform in `U′,v`) | `[O]` = H1 / coarse spectral gap; the known wall | **[O]** |
| free-field (`κ→0`) and `β=0` acid tests | `G″=κ→0`; `β=0: G″=A>N/12` | **[T]**+**[V]** |

---

## 7 · Witness — `research/v3/w7/v3_counterexample_scope.py`

Deterministic Gauss–Legendre (`n=6000` slice nodes, `240` angle nodes for `e₁`), `py -3.12`, exit 0. Key
features: cancellation-free `Φ_ρ=1/x−cot x`, `Φ_ρρ=½(−x⁻²+csc²x)`; the **correct** quadrature weight
`e^{−V}·r·w_GL` (fixing G4's omission); three-method agreement `≤10⁻⁶` quoted in §2. The pointwise Hessian
`∂²_{e₃}V` was separately verified against a finite difference at random interior points to `≤7×10⁻⁸`.
Output quoted verbatim in §2.1–2.4 and §5.

---

## 8 · Sharpest line

**Made rigorous, the counterexample is exact: `G″(e₃;θ=π,a₀=0)=A−Bβ²=0.20511530−0.01105929 β²`, negative for
`β>4.30661` (correcting G4's `β*=4.09`, whose constant `B` came from a missing quadrature weight). It is
genuinely in the admissible class of the *true coarse marginal curvature* `(C)` — longitudinal modes are in
`im P`, `a₀=0` and `B→−I` both have full support — so it kills the uniform-in-`a₀` strong form `[NG]`. But it
does not touch the literal C2 variance inequality `(U)=2W_geo` (a different object), and it does not touch the
physically-correct statement at the coarse minimiser `a₀*`, where `G″>0` trivially `[T]` and grows with `β`
`[V]`. The surviving inequality is `(★_loc)`: convexity of the coarse effective action at its own minimum,
whose uniform lower bound is the pre-existing wall H1 `[O]`.** The counterexample is a real, rigorous no-go for
the strong form; it redirects the campaign to `(C)`-at-the-minimiser, and does not close it. **No QED.
Clay 0%.**

---

*Refs.* `G4-GENERAL-BOUND.md` §0,§1,§3(iv) (the counterexample, object (C), the `(U)/(C)` split);
`S2-BOUNDARY-COUNTEREXAMPLE.md` §1–2 (model, flat-chart Stein identity, `(U)=2W_geo`);
`G1-BLOCK-IDENTITY.md` §3 (constant-`v` identity is about `(U)`); `G3-PROJECTION-TERM.md` §2 (linear `P` is
per-direction, about `(U)`); `STAR-WAVE-SYNTHESIS.md` §2 (`(★_sf)`, `(★_avg)`, `(★_{δ(U′)})`, the
uniform-in-`U′` `[NG]`); `docs/paper-v2/ch05-theorems.tex` `eq:wgeo`, `eq:split`, `eq:vardecomp`, Lemma B;
`S3-UV-FESHBACH.md` (Schur split, `S≥δ`=H1); `research/v3/w7/r1_model_curvature.py` (the `(U)`–`(C)`
reconciliation). Witness `research/v3/w7/v3_counterexample_scope.py`. No third-party sources downloaded.

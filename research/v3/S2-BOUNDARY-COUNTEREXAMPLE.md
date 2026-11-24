# S2-BOUNDARY-COUNTEREXAMPLE — falsification attack on (★) via the conditioning boundary U′

> 2026-11-24 · Attack angle S2 on the single live inequality of the campaign.
> Reads: `C2-BALABAN-UNIFORM.md` §2.4, `docs/paper-v2/ch04-gap.tex` (Savvidy cell),
> `docs/paper-v2/ch05-theorems.tex` (variance decomposition), `V3-WAVE-SYNTHESIS.md`,
> `research/unconditional/G1-L3-UV-FORCE-BOUND.md` §1.
> Witness: `research/v3/w7/s2_plaquette.py` (deterministic quadrature + rejection-MC; `py -3.12`).
>
> **Labels (charter §2).** **[T]** true theorem (and machine-checked here) · **[C]** conditional /
> proof-sketch · **[V]** numerical (deterministic quadrature, identity residual quoted) · **[NG]** no-go ·
> **[O]** exact missing step · **[R]** reported/unrefereed.

---

## 0. Verdict (no QED)

**The falsification attempt FAILS. In the exactly solvable one-plaquette model conditioned on the boundary
U′=B, (★) holds strictly, uniformly in U′ (and in the direction v), for every β>0.** No boundary
alignment — including the one aligned with the transverse (Savvidy) negative-mode direction at maximal
holonomy — makes the conditioned variance exceed E[H]+W. Precisely:

1. **[T] Exact reduction.** In the flat (exponential) chart the single-link conditioned model satisfies an
   exact Stein identity `Var_μ(G_v+Φ′_v) = E_μ[H_vv+Φ″_vv]` (density vanishes at the chart edge). With
   `E_μ[∂_v(S+Φ)]=0` this gives, **algebraically and exactly**,
   ```
     Cov_μ(G_v,Φ′_v) = W_geo(v),      and therefore   G″_model(v) := E[H_vv]−Var(G_v)+W_geo(v) = 2 W_geo(v).
   ```
   So **in this model (★) is *equivalent* to `W_geo(v|U′) ≥ δ_W/2 > 0`.** The reduction is
   direction-by-direction, so it survives any coarse projection `P` (take the fixed direction `w=Pv`). It is
   machine-verified to residual `≤5×10⁻⁷` (§3).
2. **[V] The scan finds no counterexample.** `W_geo(v|B) > 0` for all `β>0`, all boundary holonomies
   `θ∈[0,π]`, both inequivalent directions `v=e₃` (longitudinal, along the boundary axis) and `v=e₁`
   (transverse, the Savvidy-like direction). The **infimum over the boundary is attained at `θ=π`**
   (boundary `B→−I`, the far pole, maximal background field) with the **transverse direction** — exactly the
   configuration the attack targeted — and it is still **strictly positive**: at the physical head
   `β=2.5`, `min_θ W_geo = 0.0659`, giving `δ_W = 2 min = 0.132`. It increases along the UV tower.
3. **[T] Acid test (required failure).** The flat / free group (`U(1)`, `J=1`, `Φ≡0`) gives `W_geo=0` and
   `Var=E[H]`, hence `G″_model=0`: **(★) fails there, as it must.** Equivalently, SU(2) pure Haar
   (`β=0`) also gives `W_geo=0` (by `ξ→−ξ` symmetry, `E_0[Φ′]=0`). The positive floor `δ_W` is generated
   *entirely by the Wilson interaction* `β>0`; it vanishes at `β=0` and grows to the Haar floor `≈1/6`.
4. **[C] Correct restriction of (★).** `δ_W = δ_W(β_j)` is uniform in `U′` and `v` at fixed shell, and
   *increases* along the UV tower (`β_j = 4Nβ₀ j → ∞`, where `W_geo → 1/6`). It is **not** a β-independent
   constant and **cannot** be: it →0 as `β→0`. This is not a defect — it is exactly the acid test, and it is
   compatible with the campaign's per-shell reading (each shell carries its own `β_j`). Clay 0%.
5. **[O] What the model does and does not test.** It faithfully tests the boundary-conditioned fine-fiber
   inequality (fixed boundary = fixed `U′`; integrate the block interior). The one scenario it does **not**
   cover is an *additional* coarse-field constraint inside the fiber (pinning the block average, not just its
   boundary); such a constraint truncates the fine integral and leaves a genuine boundary term. That is the
   only place a counterexample could still live (§6).

---

## 1. The model (kept deliberately minimal and exact)

SU(2) lattice gauge theory, **one plaquette, one free link**, three boundary links fixed to `B∈SU(2)`:

* plaquette `P = U·B`, free link `U∈SU(2)`, action `S(U;B) = β(1 − ½ tr(UB))` (Wilson);
* conditional measure `dμ(·|B) = Z_B^{−1} e^{−S(U;B)} dU`, `dU` = Haar.  This **is** `μ(·|U′)` with `U′=B`.
* gauge invariance: set `B = exp(iθ σ₃/2)`; `μ` depends only on `θ` and the angle of `v` to the 3-axis, so
  the two inequivalent directions are `v=e₃` (along the boundary axis) and `v=e₁` (transverse).

**Flat (exponential) chart.** `U=exp(ξ)`, `ξ∈su(2)≅ℝ³`, `ρ=|ξ|`. Haar: `dU = J d³ξ`,
`J ∝ sin²(ρ/2)/ρ²`, and the **geometric potential** is `Φ := −log J = 2 log(ρ/2) − 2 log sin(ρ/2)`
(up to an additive constant). Thus `dμ = Z^{−1} e^{−𝒱} d³ξ` with the full flat-chart potential
`𝒱 = S + Φ`. The chart region is the ball `ρ∈[0,2π]` (the group `SU(2)≅S³`), and `e^{−𝒱}` **vanishes** at
its only singular edge `ρ=2π` (`sin²(ρ/2)→0`), with a regular interior point at `ρ=0`.

The objects in (★) are, for a fixed unit direction `v`:
```
  G_v = ∂_v S,   H_vv = ∂_v² S,   Φ′_v = ∂_v Φ,   Φ″_vv = ∂_v² Φ,
  W_geo(v) := E_μ[Φ″_vv] − Var_μ(Φ′_v),      G″_model(v) := E_μ[H_vv] − Var_μ(G_v) + W_geo(v).
```
Here `∂_v` is the **Cartesian** directional derivative in `ξ∈ℝ³`; `W_geo` is exactly the campaign's
Haar geometric curvature term of `ch05-theorems.tex` eq. (eq:wgeo). *(All ∂ above are along the one fixed
direction `v`; this is what appears in `⟨Pᵀ∇S,v⟩` after setting `w=Pv`.)*

---

## 2. The exact identity and the reduction to `W_geo` — **[T]**

**Lemma 1 (flat-chart Stein identity).** For every fixed unit `v` and every `β,θ`,
```
   Var_μ(G_v + Φ′_v) = E_μ[ H_vv + Φ″_vv ].
```
*Proof.* `G_v+Φ′_v=∂_v𝒱`, `H_vv+Φ″_vv=∂_v²𝒱`. Then
`∫ ∂_v(∂_v𝒱 · e^{−𝒱}) d³ξ = ∫(∂_v²𝒱) e^{−𝒱} − ∫(∂_v𝒱)² e^{−𝒱} = 0` because the integrand is a total
derivative and `(∂_v𝒱)e^{−𝒱}→0` on `ρ=2π`. Rearranged: `E[∂_v²𝒱]=E[(∂_v𝒱)²]`, and since
`E[∂_v𝒱]=−∫∂_v e^{−𝒱}=0` (Lemma 2 below) the RHS is `Var(∂_v𝒱)`. ∎

**Lemma 2.** `E_μ[∂_v𝒱]=0` (total derivative, vanishing boundary density).

**Corollary (the reduction).** `Cov_μ(G_v,Φ′_v)=W_geo(v)` and `G″_model(v)=2W_geo(v)`.
*Proof.* `E[G_v]=−E[Φ′_v]` (Lemma 2). With the cross-Stein identity
`E[∂_v𝒱·Φ′_v]=E[∂_vΦ′_v]=E[Φ″_vv]` (Lemma 1's IBP with `g=Φ′_v`),
`E[G_vΦ′_v]=E[(∂_v𝒱−Φ′_v)Φ′_v]=E[Φ″_vv]−E[(Φ′_v)²]`, hence
`Cov(G_v,Φ′_v)=E[Φ″_vv]−E[(Φ′_v)²]+(E[Φ′_v])²=E[Φ″_vv]−Var(Φ′_v)=W_geo(v)`.
Lemma 1 gives `E[H_vv]−Var(G_v)=2Cov−W_geo=W_geo`; add `W_geo` once more: `G″_model=2W_geo`. ∎

**Consequence.** In this model, (★) `⟺ 2W_geo(v|U′) ≥ δ_W`. The whole content of the UV inequality is
therefore the **strict positivity, and a uniform lower bound, of the interacting Haar geometric curvature**
`W_geo=E[Φ″]−Var(Φ′)`. The campaign proves `W_geo≥0` at `β=0` (pure Haar); the interacting case `β>0`
— which is what (★) needs per shell — is exactly what §3 decides in this model.

**Why the projection `P` is harmless.** (★) uses `⟨Pᵀ∇S,v⟩=∂_wS` with `w=Pv`, a *single fixed* direction.
Lemmas 1–2 and the corollary hold for **any fixed** `w`, so `G″_model(w)=2W_geo(w)` for every `P`. The
reduction is not an artefact of `P=I`.

---

## 3. The scan: no counterexample — **[V]**

`research/v3/w7/s2_plaquette.py`. Deterministic 2-D Gauss–Legendre quadrature in `(ρ,u=cosφ)` (900×900),
with the `ψ`-integral done exactly and `Φ_ρ, Φ_ρρ` in cancellation-free form (naive
`−2/ρ²+½csc²(ρ/2)` suffers catastrophic cancellation near `ρ=0`; that artefact produced spurious
*negative* `W_geo` in a first pass and is removed here). The built-in consistency test
`|E[H]−Var(G) − W_geo|` is **≤ 5×10⁻⁷** at every reported point, confirming both the identity and the
quadrature.

**Identity + reduction check** (excerpt; full table in the witness output):
```
 β     θ      dir   E[H]−VarG      W_geo        G″=2W_geo
 2.50  0.000  e3    0.1319564   0.1319566   0.2639130
 2.50  1.571  e1    0.1231482   0.1231483   0.2462964
 2.50  3.142  e3    0.1309120   0.1309120   0.2618240
 2.50  3.142  e1    0.0658653   0.0658653   0.1317305   ← worst case at the head
```

**`W_geo(v|B)` vs boundary holonomy `θ∈[0,π]`** (min over `θ`; `δ_W=2·min`; `argmin` in the scan is
`θ=π` unless noted):

| `β` | dir | `W(θ=0)` | `W(θ=π)` | `min_θ W_geo` | `δ_W=2 min` |
|---|---|---|---|---|---|
| 0.1 | e₃ | 0.01006 | 0.00035 | 3.53e-4 | 7.1e-4 |
| 0.1 | e₁ | 0.01006 | 0.00015 | **1.55e-4** | **3.1e-4** |
| 0.5 | e₃ | 0.04552 | 0.00860 | 8.60e-3 | 1.72e-2 |
| 0.5 | e₁ | 0.04552 | 0.00380 | 3.80e-3 | 7.59e-3 |
| 1.0 | e₃ | 0.07964 | 0.03195 | 3.20e-2 | 6.39e-2 |
| 1.0 | e₁ | 0.07964 | 0.01440 | 1.44e-2 | 2.88e-2 |
| **2.5** | e₃ | 0.13196 | 0.13091 | 0.1309 | 0.2618 |
| **2.5** | e₁ | 0.13196 | 0.06587 | **0.06587** | **0.1317** |
| 5.0 | e₃ | 0.15301 | 0.22387 | 0.1530 | 0.3060 |
| 5.0 | e₁ | 0.15301 | 0.13182 | 0.1318 | 0.2636 |
| 10.0 | e₁ | 0.16058 | 0.17225 | 0.1606 | 0.3212 |

**Read-off.**
* `W_geo(v|B)>0` everywhere for `β>0`: **(★) is never violated**, in any direction, at any boundary.
* **The worst case is exactly the targeted one:** boundary at maximal holonomy `θ=π` (`B→−I`, largest
  background field, the Savvidy regime) with the **transverse** direction `e₁` (the negative-mode plane).
  Even there `δ_W>0` (≈0.13 at the head `β=2.5`).
* `W(θ=0)` is direction-independent (axial symmetry) and **saturates to `1/6`** as `β→∞`
  (`β=10`: 0.1606 → 1/6=0.1667), the SU(2) Haar floor of Lemma B. For small `β`, `inf_θ W_geo ∝ β²`
  (`β`: 0.005→9.2e-7, 0.01→3.6e-6, 0.02→1.4e-5, all `θ=π,e₁`).
* Same order as the campaign's measured `W_geo=+0.319` at `β=2.5` (`ch05`, numerics ch.) — the model
  is in the right ballpark.

**Fine small-β scan (sign):** `inf_θ W_geo` stays **positive** down to `β=0.005` and approaches `0` *from
above* as `β↓0` (`∝β²`); it is never negative. The `θ=π` worst case is always the argmin.

---

## 4. Acid test — the free / flat limit **must** fail, and it does — **[T]**

* Flat / free space: the group is `U(1)` (or `ℝ`), `J≡1`, `Φ≡0`. Then `W_geo=0`, and `S` is a Gaussian so
  `Var(∂_vS)=E[∂_v²S]` exactly; hence `G″_model = E[H]−Var(G)+W = 0`. **(★) fails (`δ_W=0`).** Required.
* SU(2) pure Haar (`β=0`, no action): `G_v=0=H_vv`; and by the exact pure-`Φ` identity combined with
  `ξ→−ξ` symmetry, `W_geo=Var part = E_0[Φ′]² − Var_0(Φ′) ... = E_0[Φ′]² = 0`. So `G″_model=0` at
  `β=0` too, and `→0` as `β↓0` (table §3). The positive floor is **generated by `β>0`**, i.e. by the
  Wilson interaction tilting Haar and breaking the reflection that forces `E[Φ′]=0`.

This is the crucial consistency: the model does **not** manufacture positivity where the free theory has none.
The gap-causing ingredient is the same one the campaign identifies — the Haar curvature `W`, made strictly
positive and pinned near `1/6` by the interaction; on the free field it is absent and (★) collapses to
`0 ≥ δ_W`, which is false for any `δ_W>0`.

*Normalisation caveat (honest, convention-only).* The **absolute** value of `W_geo` depends on the mode
basis / Weyl-factor convention. In the clean chart used here, pure SU(2) Haar gives `W_geo=0` exactly
(`E_0[Φ′]=0` by symmetry), whereas `ch05-theorems.tex` reports a nonzero value (`0.246`) "at β=0" in its
gauge-mode probe — a different per-mode normalisation of `Φ` and of `v`. **The falsification verdict does not
depend on this:** what is convention-independent is (i) the exact reduction `G″_model=2W_geo` (Lemma 1–Cor.,
residual `≤5×10⁻⁷`), and (ii) the **sign** and **boundary-uniform positivity** of the interacting `W_geo`,
both of which the scan settles. Only conclusions (i)–(ii) are used.

---

## 5. Why the worst case is controlled (a lower bound on `δ_W(U′)`) — **[C]**

The scan is deterministic and the identity is exact, so the positivity of `W_geo` in the scanned box is
**[V]**. A mechanism/proof-sketch that it is *not* a numerical accident:

* **Gaussian limit (exact, [T]-in-limit).** Near the action minimum `ξ₀=log B†`, approximate
  `S ≈ ½ a|ξ−ξ₀|²` (action curvature `a=O(β)`, `a↑` with `β`) and `Φ ≈ ½ k|ξ|²` with `k=1/6`
  (`Φ″⪰1/6`, Lemma B). With `G=a(ξ−ξ₀)`, `Φ′=kξ`, `ξ` Gaussian of variance `(a+k)^{−1}`,
  `W_geo = Cov(G,Φ′) = ak/(a+k) > 0`. This reproduces the observed limits: `→k=1/6` as `a→∞`
  (`β→∞`), `∝β` as `β→0`. On the Gaussian approximation the bound is **uniform in the boundary** `ξ₀`
  (`ak/(a+k)` is `ξ₀`-independent): boundary alignment cannot hurt.
* **Worst case = maximal geodesic displacement.** Varying `θ` moves the measure's mean along the geodesic
  from `0` to `−I` (distance `π`). `θ=π` maximizes the spread of `Φ′` and hence `Var(Φ′)` — the only term
  that subtracts in `W_geo` — so it is the natural argmin, and the transverse `e₁` (where `Φ″` is smallest in
  the spread directions) is the second argmin. Both are captured by the scan; neither crosses `0`.
* **Uniform lower bound.** For each fixed `β>0`, `inf_{U′,v} W_geo = W_geo(θ=π, e₁)` is a finite positive
  constant (table §3), so `δ_W(β) := 2·inf_{U′,v} W_geo > 0` is **uniform in the boundary and direction**,
  and increases along the tower. The bound is *not* uniform in `β` and cannot be: `δ_W(β)→0` as `β→0`
  (the acid test). Along the campaign tower `β_j=4Nβ₀ j → ∞`, so `δ_W(β_j)` *grows*; the worst shell is the
  IR head `β≈2.5`, where `δ_W ≈ 0.13`.

So the correct restriction of (★) is: **`δ_W = δ_W(β_j) > 0`, uniform in `U′` and `v` (and in `a`, the
model has no `a`), but not β-independent.** Stated uniformly in `(j,U′,a)` with a *single* constant it is
false only in the trivial free-field limit, exactly as the acid test demands.

---

## 6. Scope: where a counterexample could still hide — **[O]**

The model integrates the **entire** fine fibre at fixed boundary `U′=B`. Its reduction `G″=2W_geo` used only
`e^{−𝒱}` vanishing on the fibre boundary. The one ingredient it does **not** model is a **coarse-field
constraint inside the fibre** — pinning the block *average* (or a coarse block-spin) to a saddle, not merely
fixing its boundary. Such an interior constraint truncates the fine integral over a sub-fibre with a
non-vanishing edge, so Lemma 1 acquires a **genuine boundary term** — precisely the
"boundary/divergence defect" named in `G1-L3-UV-FORCE-BOUND.md` §1. That term is the only object left that
could push `Var−E[H]−W` positive. The present result sharpens the target:

> **In the boundary-conditioned fibre (no interior coarse constraint), (★) reduces exactly to
> `2W_geo(v|U′) ≥ δ_W`, and this holds with a positive, boundary-uniform floor generated by `β>0`.**
> The surviving risk is *not* the boundary holonomy or the negative-mode direction (both scanned, both
> safe); it is the **interior coarse-projection boundary term**, i.e. the `(★dec)`/clustering half.

This is consistent with `C2-BALABAN-UNIFORM.md` §2.4, which already located the residual in the small-field
conditional analytic bound rather than the (now-excluded) large-field entropy channel.

---

## 7. Status ledger

| item | statement | label |
|---|---|---|
| flat-chart Stein identity `Var(G_v+Φ′_v)=E[H_vv+Φ″_vv]` in the single-link model | IBP, vanishing edge density; residual ≤5e-7 | **[T]** (+[V]) |
| reduction `G″_model(v)=Cov(G_v,Φ′_v)=2W_geo(v)`; survives any projection `P` (per-direction) | §2 | **[T]** |
| scan `W_geo(v|B)>0` for all `β>0`, `θ∈[0,π]`, `v∈{e₃,e₁}`; worst case `θ=π`, transverse | witness, deterministic quadrature | **[V]** |
| `δ_W(β)=2 inf_{U′,v} W_geo > 0` uniform in `U′,v`; `δ_W(β)→0` as `β→0`; grows along tower | §3, §5 | **[C]** (+[V]) |
| free/flat (`U(1)`, `Φ=0`) and SU(2) Haar `β=0`: `G″_model=0`, (★) fails | acid test | **[T]** |
| (★) with a single β-independent `δ_W` uniform down to `β→0` | **false** — must exclude the free field | **[NG]** (of that reading) |
| counterexample to (★) via boundary / negative-mode alignment | **not found**; (★) survives in the model | **[NG] refuted** |
| interior coarse-projection boundary term / `(★dec)` | untested here; only remaining hiding place | **[O]** |

**Sharpest line.** A boundary aligned with the Savvidy negative mode is the *worst case* of the exactly
solvable conditioned one-plaquette model — and it still gives `δ_W(β≈2.5)≈0.13>0`. The attempt to falsify
(★) on the boundary fails: in the boundary-conditioned fibre, (★) is *exactly* `2W_geo(v|U′)≥δ_W`, the
interacting Haar curvature stays strictly positive, and its floor is generated, as it must be, by the
Wilson interaction and not by free-field kinematics. The one door left open is not the boundary but the
interior coarse projection. Clay 0%.

---

## 8. Verification ledger

| item | checked | how |
|---|---|---|
| model action/measure, chart, `Φ=−log J` | exact; `S` matches `½tr(UB)=cos(ρ/2)cos(θ/2)−sin(ρ/2)sin(θ/2)u` | analytic + witness |
| flat-chart identity residual `E[H]−Var(G)−W_geo` | `≤5×10⁻⁷` at all reported points | `s2_plaquette.py` |
| cancellation-free `Φ_ρ=1/x−cot x`, `Φ_ρρ=½(−x⁻²+csc²x)`, `x=ρ/2` | Taylor-matched at `x=0.3`; naive form gives spurious negatives (discarded) | witness note |
| `W_geo` sign/scan over `β∈[0.005,10]`, `θ∈[0,π]`, `v∈{e₃,e₁}` | all positive; `min` at `θ=π`, transverse; `∝β²` small-`β`, `→1/6` large-`β` | witness |
| `W_geo(β=2.5)∈[0.066,0.275]` | same order as campaign measured `+0.319` | consistency |
| free-field / flat-group failure `G″=0` | `Φ≡0 ⇒ W=0`, Gaussian `Var=E[H]` | analytic |
| SU(2) Haar `β=0`, `W_geo=0` | pure-`Φ` Stein + `ξ→−ξ` symmetry | analytic + `β↓0` limit |

*Refs.* `C2-BALABAN-UNIFORM.md` §2.4; `ch04-gap.tex` §open cell / Savvidy; `ch05-theorems.tex`
eqs. (eq:wgeo, eq:split, eq:vardecomp); `V3-WAVE-SYNTHESIS.md` §2–3; `G1-L3-UV-FORCE-BOUND.md` §1; witness
`research/v3/w7/s2_plaquette.py`.

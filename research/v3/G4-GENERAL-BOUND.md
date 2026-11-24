# G4-GENERAL-BOUND — the general-block fate of (★): the S2 reduction `(★)⇔2W_geo` is invalid, the true coarse curvature `G″` is *not* `2W_geo`, it admits `G″<0`, and `W_geo>0` at β>0 is a real open lemma but is *not* the gate

> 2026-11-24 · Attack angle G4 on the single live inequality `(★)` of the campaign, generalizing the
> exactly solvable S2 one-plaquette result to a **general block and general boundary**.
> Reads: `S2-BOUNDARY-COUNTEREXAMPLE.md` (§2 reduction, §3 scan), `S1-CONDITIONAL-POINCARE.md`,
> `S3-UV-FESHBACH.md`, `C2-BALABAN-UNIFORM.md` §2.4, `docs/paper-v2/ch05-theorems.tex`
> (`eq:split`, `eq:vardecomp`, `sec:wgeo`), `V3-WAVE-SYNTHESIS.md`, `STAR-WAVE-SYNTHESIS.md`.
> Witness: `research/v3/w7/g4_general_bound.py` (`py -3.12`, deterministic Gauss–Legendre quadrature,
> output quoted below).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional / proof-sketch · **[V]** numerical
> (deterministic quadrature; identity residual quoted) · **[NG]** no-go for a named class · **[O]** exact
> missing step · **[R]** reported/unrefereed.

---

## 0 · Verdict (no QED; Clay 0%)

**The task's premise is the finding.** The candidate reduction left by S2, `(★) ⇔ 2W_geo(Pw) ≥ δ_W`, is
**not valid** — it compares the wrong functional — and no general bound `G″ ≥ 2W_geo`, `G″ ≥ c·W_geo`, or
`G″ ≥ 2W_geo − C` is true. Concretely:

1. **[T] S2 computed the *unconditional* `Δ_dyn+W_geo`, not the coarse curvature.** The genuine object in
   `(★)` is the **coarse-mode-conditioned marginal curvature** `G″(v;a₀) = −d²/da² log ρ_v(a)|_{a₀}`
   (`ch05` `eq:split`/`eq:vardecomp`, the `|0` made explicit). S2's `G″_model := E[H]−Var(G)+W_geo`
   omits the cross term `X = −2Cov(G,Φ′)` and uses **unconditioned** moments. By the exact
   Stein + cross-Stein identities (Lemma G2 below), `Δ_dyn := E[H]−Var(G) = W_geo = Cov(G,Φ′)` holds for
   **any** flat-chart block, so S2's `G″_model = 2W_geo` is a correct *algebraic identity* — but for the
   wrong object. The literal `(★)` G″ (`Δ_dyn+W_geo+X`, with `X=−2W_geo`) is **identically 0**.
2. **[T]/[V] The true `G″` is NOT `2W_geo` and changes sign.** The correct coarse curvature is
   `G″ = E[Σ]−D` (S3). In the S2 model it equals `2W_geo` at no boundary except accidentally: at `θ=0`
   it is `≈4×` larger (`β=5`: `G″=1.36` vs `2W_geo=0.31`); at maximal holonomy, longitudinal mode, and the
   **barrier of the coarse double well**, it goes **negative** (`β=5`: `G″=−0.10`; `β=10`: `−1.02`;
   `β=20`: `−4.71`) while `2W_geo>0` throughout. So `(★)` is *mis-certified* by the reduction exactly
   where it actually fails.
3. **[NG] There is no universal lower bound of the forms (i)–(iii).** Counterexample (iv) exists:
   `G″<0` in the one-plaquette model at `β>β*=4.09`, boundary `B→−I`, coarse field at the double-well
   barrier. This is the coarse-graining image of the Savvidy/S1 double well, and (mechanism) persists for
   general `SU(N)` (`[C]`).
4. **[T] The clean β=0 contradiction pins the defect.** At β=0 (pure Haar) `W_geo=0` while `ch05`
   *proves* `G″ ≥ N/12 > 0`; the reduction `(★)⇔2W_geo≥δ_W` would predict `(★)` **false** at β=0, where it
   is a **theorem**. No further refutation is needed.
5. **[O] `W_geo>0` at β>0 is real but not the gate.** `Prekopa`/Brascamp–Lieb give `W_geo≥0` only for the
   log-concave `β=0` Haar marginal (where in the centered chart it is `=0`). At β>0 the interacting measure
   is not log-concave (S1 `Prop U′-NG`), so `W_geo>0` is **not** a Prekopa corollary; the exact missing
   lemma is a BL-type variance bound `Var_{μ_{U′}}(Φ′) ≤ E_{μ_{U′}}[Φ″]−c` under the interacting measure.
   **But `W_geo>0` is not what `(★)` needs**, since `(★) ≠ 2W_geo`; the gate is the coarse marginal
   positivity `E[Σ]≥δ_W` (= the coarse spectral gap, H1) at the RG-relevant coarse field.

**Acid test passed.** Free/flat field (`Φ≡0`): `W_geo=0`, Gaussian `Var=E[H]`, hence `G″=0` — no false
positivity; `(★)` fails for any `δ_W>0`, as required. β=0 Haar: `G″→N/12>0` but `W_geo→0`.

---

## 1 · Three distinct objects — and which one is `(★)`

Fix a block `B` and boundary `U′`. Work in the flat product chart `X=(X_ℓ)_{ℓ∈B}∈⊕_ℓ𝔰𝔲(N)`, Haar
`dU_ℓ=e^{−Φ_0(X_ℓ)}dX_ℓ`; the fiber measure is
```
  μ_{U′}(dX) = Z(U′)^{−1} e^{−V(X)} dX ,   V = S_W(X;U′) + Φ(X) ,   Φ = Σ_{ℓ∈B} Φ_0(X_ℓ),     (1.1)
```
with `Φ_0` the single-link Haar potential (`Φ_0''⪰(N/12)I`, `ch05` Lemma B). Because `Φ_0→+∞` on each
link's chart edge, `e^{−V}` **vanishes on the chart boundary** — the only structural input any of the
identities below need.

For a unit vector `v` in the fiber tangent space put `G_v := ∂_vS_W`, `H_vv := ∂_v²S_W`,
`Φ′_v := ∂_vΦ`, `Φ″_vv := ∂_v²Φ`, and the coarse coordinate `a := ⟨v,X⟩`. There are **three** different
second-variation objects in the corpus; they must not be conflated:

| symbol | definition | what it is |
|---|---|---|
| **(L) literal `(★)`** | `G″_lit := E_{μ_{U′}}[H_vv] − Var_{μ_{U′}}(G_v) + W_j` , `W_j=W_geo+X` | C2's written form, **no** coarse-mode conditioning |
| **(U) S2's model** | `G″_model := E_{μ_{U′}}[H_vv] − Var_{μ_{U′}}(G_v) + W_geo` | S2's `G″_model`, **no** `X`, **no** conditioning |
| **(C) true coarse curvature** | `G″(v;a₀) := −d²/da² log ρ_v(a)|_{a₀}` , `ρ_v` = marginal of `a` | `ch05` `eq:split`/`eq:vardecomp`, conditioning `|0` made explicit |

where `W_geo := E[Φ″_vv]−Var(Φ′_v)` and `X := −2Cov(G_v,Φ′_v)`.

**Claim (C is `(★)`).** The quantity that controls the **coarse** measure (the effective action
`S_eff(a)=−log ρ_v(a)` and its Hessian `G″(a)=S_eff''(a)`) — the object `(★)` is *about* — is (C). C2's
`G″_j(v)=E[H_j]−Var_j+W_j≥δ_W` is a shorthand for (C) with the `|0` suppressed (compare `ch05`
`eq:split`: `G″(0)=E[H+Φ″|ξ=0]−Var(G+Φ′|ξ=0)=Δ_dyn+W_geo+X`). **[T]**

---

## 2 · The exact identities for a general block — **[T]**

> **Lemma G1 (marginal-curvature identity; no IBP, any smooth `V`).** For the measure (1.1), every unit
> `v`, and `ρ_v(a)=∫_{⟨v,X⟩=a} e^{−V}`,
> ```
>    G″(v;a₀) = −d²/da² log ρ_v(a)|_{a₀} = E_{μ_{U′}}[∂_v²V | a=a₀] − Var_{μ_{U′}}(∂_vV | a=a₀).     (2.1)
> ```
> *Proof.* `ρ′(a)=−ρ(a)E[∂_vV|a]`, `ρ″(a)=ρ(a)E[(∂_vV)²−∂_v²V|a]`; combine in
> `−ρ″/ρ+(ρ′/ρ)²`. Pure differentiation under the integral; valid because `ρ_v>0`. ∎

> **Lemma G2 (unconditional Stein + cross-Stein; the vanishing-boundary input).** Over the **full** block
> measure (~no coarse conditioning), the vanishing of `e^{−V}` on the chart edge gives, for every `v`,
> ```
>    E[∂_vV] = 0 ,   E[∂_v²V] = E[(∂_vV)²]              (Stein, g=∂_vV),
>    E[∂_vV·Φ′_v] = E[Φ″_vv]                              (cross-Stein, g=Φ′_v).                      (2.2)
> ```
> Consequently, with `H=∂_v²S_W`, `G=∂_vS_W`:
> ```
>    Δ_dyn := E[H]−Var(G) = E[Φ″]−Var(Φ′) = Cov(G,Φ′) =: W_geo,                                     (2.3)
>    Δ_dyn + W_geo = 2 W_geo ,     G″_lit = Δ_dyn + W_geo + X = 2Cov − 2Cov = 0 .                   (2.4)
> ```
> *Proof.* (2.2): total derivatives `∫∂_v[(∂_vV)e^{−V}]=0`, `∫∂_v[(Φ′_v)e^{−V}]=0`; boundary term is
> `(∂_vV)e^{−V}` resp. `Φ′_v e^{−V}` on the edge, which vanish since `e^{−V}→0` there (and `∂_vV,Φ′` grow
> at most like `1/(edge-distance)`, beaten by `e^{−V}~(\text{edge-dist})²`). Then (2.3) is S2's §2
> computation, algebraic from (2.2) and `E[G]=−E[Φ′]`; (2.4) follows with `X=−2Cov`. ∎ Verified
> numerically: `Δ_dyn=W_geo` to `≤10⁻⁶` at every `(β,θ,v)` in witness table (A).

**Three immediate consequences, all for a general block:**

- **(2.3) is a clean general identity: the force-variance defect *equals* the Haar geometric curvature**
  (`= the force–Jacobian covariance`). It is *not* special to one plaquette; it holds for any block, any
  boundary, any `β`. On the free field `Φ≡0 ⇒ W_geo=0 ⇒ E[H]=Var(G)` (Gaussian equipartition) — the acid
  test is built into the identity.
- **(2.4): the literal `(★)` has `G″_lit ≡ 0`.** Read without the coarse-mode conditioning, `(★)` says
  `0 ≥ δ_W`, i.e. **false for every `δ_W>0`, in every model** (including the S2 model). This is not a
  counterexample to the *intended* `(★)`; it shows the conditioning `|a=a₀` is **load-bearing** and cannot
  be dropped. The well-posed `(★)` is (C), not (L).
- **S2's `G″_model = 2W_geo` is correct but is `(U)`, not `(C)`.** So `(★) ⇔ 2W_geo ≥ δ_W` fails under
  **either** reading: read as (L), `(★)` is `0≥δ_W` while S2 compares `2W_geo` (irrelevant); read as (C),
  `G″(v;a₀) ≠ 2W_geo` (different functional). This is the decisive point of G4.

> **Lemma G3 (S3/Schur form of the true coarse curvature).** With `Π₀=I−vvᵀ`,
> `K=Π₀ Hess V Π₀`, `c=Π₀ Hess V v`, `Σ=vᵀ Hess V v − ⟨K⁻¹c,c⟩` (Hess in the full fiber),
> ```
>    G″(v;a₀) = E_{μ_{U′}}[Σ | a₀] − D ,   D := Var(∂_vV|a₀) − E[⟨K⁻¹c,c⟩|a₀] ,                    (2.5)
> ```
> and `K⪰0` (log-concave fiber) `⇒ D≤0 ⇒ G″ ≥ E[Σ]`. `[T]` (S3 Lemma S; Gaussian Schur = marginal
> curvature, witness `s3_feshbach_witness.py`).

---

## 3 · The attempted inequalities (i)–(iii) and the counterexample (iv) — **[T]/[V]/[NG]**

All numbers are the witness output. "2W_geo" is S2's unconditional quantity (table A); "G″" is the true
coarse curvature (C) (table B), computed two ways (analytic slice integral and finite-difference of
`−log ρ`), agreeing to the quadrature.

**(i) `G″ ≥ 2W_geo` for all blocks/boundaries — FALSE.**
- At `θ=0` (no boundary field) `G″ ≫ 2W_geo`: `β=2.5,e₃`: `G″=0.741` vs `2W_geo=0.264` (×2.8);
  `β=5`: `1.358` vs `0.306` (×4.4); `β=10`: `2.605` vs `0.321`.
- At `θ=π` (max holonomy), longitudinal `e₃`: `β=5`: `G″=−0.102 < 0 < 2W_geo=0.448`. So `G″≥2W_geo`
  fails in the strongest possible way (opposite sign).

**(ii) `G″ ≥ c·W_geo`, `c>0` — FALSE.** Same point: `G″<0` while `W_geo>0`. No positive `c` can hold
across the interface. (On the log-concave small-field set the true floor is `c·(N/12)` from `Φ″⪰N/12`,
**not** `c·W_geo` — see (iii)/§4.)

**(iii) `G″ ≥ 2W_geo − C(U′)` with controlled `C` — not meaningful / FALSE.** The `2W_geo` baseline is the
wrong object, so no such inequality is a theorem. The *correct* remainder form is (2.5),
`G″ = E[Σ] − D`: `G″ ≥ E[Σ] − D⁺`, where `D⁺=(D)⁺` is the BL defect, `=0` on log-concave fibers and
**unbounded** on Savvidy boundaries. So a "`−C(U′)`" bound holds only with `C=D⁺`, which is *not uniformly
controlled* — it is exactly the S1 `[NG]`.

**(iv) Counterexample `G″ < 2W_geo`, indeed `G″ < 0` — FOUND. [T] in the model, [V] numbers.**

At boundary `B=e^{iπσ₃/2}` (maximal holonomy, `θ=π`), **longitudinal** mode `v=e₃`, the marginal curvature
factorises exactly (`V=h(ρ)+g(ρ)u`, `u=ξ₃/ρ`; on the slice `ξ₃=a₀=0` the measure is `β`-independent Haar,
`V_3=g(r)/r`, `V_33=h′(r)/r`):
```
  G″(e₃; θ=π) = E_{Haar}[Φ′(r)/r] − β² · Var_{Haar}[sin(r/2)/r] = 0.2057 − 0.0123 β² .            (3.1)
```
So `G″(e₃;θ=π) = 0` at `β* ≈ 4.09` and `G″<0` for `β>β*`: `β=5: −0.102`, `β=10: −1.024`, `β=20: −4.71`,
while the unconditional `2W_geo>0` at every one (`β=5:0.448`, `β=10:0.532`). **The reduction would
certify `(★)` (since `2W_geo>δ`) precisely where the true `G″<0` and `(★)` fails.**

**It is the *longitudinal* mode, not the transverse.** S2 targeted the transverse `e₁` as the worst case;
for the true `G″`, `e₁` stays **positive** (`β=10:1.118`, `β=20:2.131`), while `e₃` turns negative. The
tachyonic direction at maximal holonomy is the one **along the holonomy axis**.

**Double-well anatomy.** `S_eff(a₀)=−log ρ_{e₃}(a₀)` at `θ=π` is a double well (the classical
Savvidy/first-order picture of S1 §4.2): at `β=5` the density peaks at `a₀≈−2` (the coarse **minimum**,
`G″=+1.44` at `a₀=−3`) and has a **barrier** at `a₀≳0` where `G″=−0.40` (`a₀=1.5`) and `−0.19`
(`a₀=3`). The sign of `G″` is exactly the sign of `S_eff''`: positive at the well, negative at the
barrier. So `(★)` holds at the RG minimum and **fails at the barrier** — i.e. the coarse effective action
is **not convex** at strong boundary field. This is the coarse image of the S1 `[NG]` and of
STAR-WAVE's "uniform in all `U′` is false."

**What *is* provable (the honest bound). [T] given log-concavity, [C] on the good set.** On the fiber with
`K⪰0` (small-field boundary, S1 Lemma S1/Cor), Prekopa+Brascamp–Lieb give `D≤0`, hence (2.5)
`G″ ≥ E[Σ] ≥ E[Φ″_n] − E[⟨K⁻¹c,c⟩] ≥ N/12 − \text{Schur loss}`. So the true universal lower bound, *on the
log-concave set*, is the **Haar floor `N/12` minus the Schur loss** — **not** `2W_geo`, **not** `W_geo`.
On the full boundary range even this fails (the counterexample has `K⋡0`). This matches S3 and S1 exactly.

---

## 4 · `W_geo>0` at β>0: what Prekopa/BL can and cannot give — **[T]/[NG]/[O]**

Write `W_geo := E_{μ_{U′}}[Φ″] − Var_{μ_{U′}}(Φ′)` under the **interacting** measure (the object `ch05`
measures as `+0.319` at β=2.5 and that S2's `2W_geo` refers to). By Lemma G2, `W_geo = Cov_{μ_{U′}}(G,Φ′)`.

**(a) What Prekopa gives — and only that. [T]** For the **pure-Haar** measure `e^{−Φ}dX` (`β=0`), `Φ` is
convex (`Φ″⪰(N/12)I`), so the measure is log-concave and Prekopa gives `W_geo ≥ 0`; in the centered chart
`ξ→−ξ` is a symmetry, so `E[Φ′]=0` and Stein `Var(Φ′)=E[Φ″]`, i.e. **`W_geo=0` exactly** at β=0
(witness table A: `β=10⁻⁴`, `W_geo=3.5×10⁻⁸`). So `Prekopa`'s `W_geo≥0` is **saturated** on the free Haar
measure; it carries no positivity for β>0.

**(b) Why Prekopa/BL do NOT transfer to β>0. [NG] for the method.** At β>0 the interacting measure
`e^{−S_W−Φ}` is **not** log-concave: `Hess(S_W)` is Savvidy-indefinite (`S1 Prop U′-NG`; W6 `λ_min` to
`−8` for chromomagnetic boundary fields). Brascamp–Lieb (`D≤0`) needs `K⪰0`; it fails. So there is **no
Prekopa/BL route** to `W_geo>0` at β>0.

**(c) The exact missing lemma. [O]** `W_geo = E[Φ″]−Var(Φ′)`; the deficit is purely `Var(Φ′)`:
```
   (Missing Lemma W+)    Var_{μ_{U′}}(∂_vΦ)  ≤  E_{μ_{U′}}[∂_v²Φ] − c ,     c = c(β)>0 ,              (4.1)
```
i.e. a Brascamp–Lieb-type **variance bound for the Haar force under the interacting measure**. Equivalently
`Cov_{μ_{U′}}(G,Φ′) ≥ c`. This is the same operator inequality `K⪰0` / effective log-concavity that is the
Savvidy wall; it can be certified only on the small-field set (where it is `[C]` via S1's Corollary) and is
**open** at the head. Witness table A shows `W_geo>0` at all tested `(β≤10,θ∈{0,π/2,π},v)` — numerically
robust, but the non-log-concave measure gives no proof.

**(d) `W_geo>0` is NOT the gate for `(★)`. [T].** Since `(★) ≠ 2W_geo` (§2), positivity of `W_geo` is a
statement about the *unconditional* `Δ_dyn+W_geo` object, not about the coarse curvature `G″=E[Σ]−D`. The
gate for `(★)` is `E[Σ]≥δ_W` at the RG-relevant coarse field — the **coarse spectral gap**, i.e. H1 (the
gap) — exactly as S3/S1 conclude. `W_geo>0` is a genuine and interesting open lemma about `ch05`'s sign
theorem and the measured decomposition; it is **not** the decisive step. (The task's conjecture that it is
"the real gate" is correct only under S2's reduction, which §2–3 refute.)

---

## 5 · What would count as proving / disproving `(★)`, given the S2 reduction is void

**The S2 reduction `(★)⇔2W_geo≥δ_W` cannot be used as the definition of the task.** Replace it by (C)/(2.5):

> **Proving `(★)`** = exhibiting `δ_W>0` such that, uniformly in shell `j`, boundary `U′`, and lattice
> spacing `a`, the **coarse marginal curvature at the RG-relevant coarse field** satisfies
> ```
>    G″_j(v;a₀(U′)) = −d²/da² log ρ^{(j)}_{v}(a)|_{a₀} = E_{μ_{U′}}[Σ_j|a₀] − D_j ≥ δ_W        (5.1)
> ```
> for every coarse mode `v` (equivalently: the coarse effective action `S_eff,j` is `δ_W`-strongly convex
> along every `v`). By (2.5) this is the conjunction `[fiber log-concavity enough that D≤E[Σ]−δ_W] ∧
> [E[Σ]≥δ_W]`; the second factor is the coarse gap `H1`, the first is the small-field/`(★_sf)` input.
> *Reading:* uniform in `U′` throughout is **false** (§3); the tenable statements are `(★_sf)` (uniform on
> the small-field set `𝒮_t`, S1 Cor), `(★_avg)` (averaged over `U′`), `(★_{δ(U′)})` (boundary-dependent
> margin) — exactly STAR-WAVE §2.

> **Disproving `(★)`** = a single `(block, boundary U′, coarse field a₀, mode v)` with
> `G″(v;a₀) < δ_W` (in particular `<0`). **This G4 supplies one** (§3(iv)): one plaquette, `B→−I`,
> `v=e₃` longitudinal, `a₀` at the double-well barrier, `β>β*=4.09` (`β=5: G″=−0.10`), while the S2
> functional `2W_geo=0.45>0`. So the **uniform-in-all-`U′`** reading of `(★)` is **disproved**, matching
> S1 `[NG]`; the **small-field / averaged** readings survive.

**Crucial bookkeeping correction.** `2W_geo` is *not* a lower bound for `G″` and *not* an upper bound: the
true `G″` is above `2W_geo` on part of the domain (`θ=0`, small `β`) and below it (negative) elsewhere
(§3). Any proof or disproof phrased through `2W_geo` is answering the wrong question.

---

## 6 · Acid tests (the free field must fail) — **[T]/[V]**

1. **Flat / free field (`Φ≡0`, `U(1)`/abelian):** `W_geo=0`; `S` Gaussian ⇒ `Var(∂_vV)=E[∂_v²V]` ⇒
   `G″=0` exactly. `(★)` fails for any `δ_W>0`. **Required, and holds.** (Lemma G2 with `Φ=0`:
   `E[H]=Var(G)`, `W_geo=0`.) **[T]**
2. **β=0 pure Haar on `SU(2)`:** `W_geo=0` (saturated Prekopa; `ξ→−ξ` symmetry), but the **true** coarse
   curvature `G″(e₃;θ=π)=0.206 > N/12=0.167` (witness; consistent with `ch05`'s `β=0` theorem
   `G″≥N/12`). So `(★)` **holds** at β=0 (for `δ_W<0.206`) while `2W_geo=0<δ_W`. **The S2 reduction
   predicts failure where there is a theorem** — the cleanest one-line refutation of §2. **[T]+[V]**
3. **Continuity:** `G″(e₃;θ=π) = 0.2057 − 0.0123 β²` is continuous in `β`, `=N/12`-floored at β=0, and
   crosses zero at `β*=4.09` — the positivity is generated by the Haar floor at small `β` and destroyed by
   the Wilson double-well at large `β`. **[V]**

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| marginal-curvature identity `G″(v;a₀)=−d²/da²logρ|a₀=E[∂_v²V|a₀]−Var(∂_vV|a₀)` | Lemma G1, any block, no IBP | **[T]** |
| unconditional Stein+cross-Stein; `Δ_dyn=W_geo=Cov(G,Φ′)`; `Δ_dyn+W_geo=2W_geo`; `G″_lit≡0` | Lemma G2, any block (`Δ_dyn=W_geo` to `≤1e-6`) | **[T]** |
| S2 reduction `(★)⇔2W_geo≥δ_W` | **invalid**: compares (U) not (L)/(C); mis-certifies at the counterexample | **[NG]** (of the reduction) |
| `β=0`: `(★)` is a theorem (`G″≥N/12`) with `2W_geo=0` | direct contradiction of the reduction | **[T]** (+**[V]**) |
| `G″ = E[Σ]−D`; `K⪰0⇒G″≥E[Σ]` | S3 Lemma S; the correct coarse-curvature form | **[T]** |
| (i) `G″≥2W_geo` universally | **false**; `G″≫2W_geo` at `θ=0`, `G″<0` at barrier | **[NG]**+**[V]** |
| (ii) `G″≥c·W_geo`, `c>0` | **false**; `G″<0` where `W_geo>0` | **[NG]**+**[V]** |
| (iii) `G″≥2W_geo−C(U′)` bounded `C` | **false**; correct remainder `−D⁺` unbounded on Savvidy | **[NG]** |
| (iv) counterexample `G″<2W_geo`, `G″<0` | `e₃`, `θ=π`, `β>4.09`: `G″=0.206−0.0123β²` (`β=5:−0.10`) | **[T]**(model)+**[V]** |
| provable bound on log-concave fiber | `G″≥E[Σ]≥N/12−`Schur-loss (not `2W_geo`) | **[T]** given `K⪰0`, **[C]** on `𝒮_t` |
| universal (all `U′`) positivity of `G″` | **false** (double well at strong field; S1 `[NG]` restated) | **[NG]** |
| `W_geo>0` at β>0 via Prekopa/BL | **no**: interacting measure not log-concave; β=0 case `=0` | **[NG]** (method) |
| exact missing lemma for `W_geo>0` | `Var_{μ_{U′}}(Φ′)≤E[Φ″]−c` (BL for Haar force, interacting μ) | **[O]** |
| `W_geo>0` is the gate for `(★)` | **no** — `(★)≠2W_geo`; gate is `E[Σ]≥δ_W` = coarse gap `H1` | **[T]** |
| counterexample generalizes to `SU(N)` | Cartan double well along maximal-holonomy axis | **[C]** |
| free field / β=0 acid tests | `Φ=0⇒G″=0`; β=0 `G″>0, W_geo=0` (reduction fails) | **[T]**+**[V]** |

---

## 8 · Witness — `research/v3/w7/g4_general_bound.py` (deterministic, `py -3.12`)

Quoted excerpt:

```
(A) UNCONDITIONAL:  Delta_dyn = W_geo,  Delta_dyn+W_geo = 2 W_geo   (residual <=1e-6)
    beta  th      dir    Delta_dyn      W_geo       2W_geo
    2.50  0.000   e3      0.131956    0.131957    0.263913
    2.50  3.142   e3      0.130912    0.130912    0.261824
    10.0  3.142   e3      0.265782    0.265782    0.531564

(B) TRUE coarse curvature G''(v)=-d^2/da^2 log rho|_0   vs 2W_geo  vs 0
    beta  th      dir   Gpp(analytic)  Gpp(fd e3)   2W_geo       Gpp<0?
    0.10  0.000   e3      0.221681      0.221062    0.020112     False
    2.50  0.000   e3      0.740664      0.729924    0.263913     False
    5.00  0.000   e3      1.357639      1.342360    0.306023     False
    2.50  3.142   e1      0.362312        --        0.131731     False
    5.00  3.142   e3     -0.101786     -0.071367    0.447746     True   <-- counterexample
    10.0  3.142   e3     -1.024306     -0.900814    0.531564     True
    20.0  3.142   e3     -4.714384        --             --      True

(C) acid tests:  beta=1e-4 pure Haar:  G''(e3)=0.2057 (>N/12)  while  W_geo=3.5e-8 (~0)

(D) beta-scan at theta=pi, e3:  G''=0.206-0.0123 beta^2
    beta:  0   1      2      3      4      5      6      8      10
    G'': .206 .193  .157  .095  .0089 -.102  -.237  -.581  -1.024
```

The finite-difference column (`−d²/da² log ρ` by 5-point stencil) confirms the analytic slice integral;
the model and its `2W_geo` reproduce S2 (`s2_plaquette.py`) to quadrature.

---

## 9 · Sharpest line

**S2 proved a true identity — `Δ_dyn+W_geo=2W_geo` for every flat-chart block — but for the wrong object:
the `(★)` curvature is the coarse marginal curvature `G″(v;a₀)=−d²/da²logρ|_{a₀}=E[Σ]−D`, which is not
`2W_geo`, is `≈4×` larger at `θ=0`, and turns *negative* at the double-well barrier of a maximal-holonomy
boundary (`G″=0.206−0.0123β²`, `β>4.09`).** Hence the candidate reduction `(★)⇔2W_geo≥δ_W` is void: it
would certify `(★)` exactly where `G″<0` and mis-predict failure at β=0 where `(★)` is a theorem.
There is therefore **no general bound `G″≥2W_geo`, `≥cW_geo`, or `≥2W_geo−C`**, and a **counterexample with
`G″<0` exists**; the only provable lower bound is the Haar floor `N/12` minus the Schur loss, available
exactly on the log-concave small-field set. `W_geo>0` at β>0 is a real open lemma (Prekopa/BL give only
the saturated `β=0` case; the missing step is a BL variance bound under the non-log-concave interacting
measure), but it is **not** the gate: the gate remains the coarse gap `E[Σ]≥δ_W` = H1. **No QED. Clay 0%.**

---

*Refs.* `S2-BOUNDARY-COUNTEREXAMPLE.md` §1–4; `S1-CONDITIONAL-POINCARE.md` §1.3, §4.2, §6;
`S3-UV-FESHBACH.md` §1–4 (Lemma S, Cor. S, `G″=E[Σ]−D`); `C2-BALABAN-UNIFORM.md` §2.4;
`docs/paper-v2/ch05-theorems.tex` `eq:split`, `eq:vardecomp`, `sec:wgeo` (Lemma B `Φ″⪰N/12`);
`STAR-WAVE-SYNTHESIS.md` §2; Brascamp–Lieb (1976); Prékopa (1971). Witness `research/v3/w7/g4_general_bound.py`;
model cross-checked against `research/v3/w7/s2_plaquette.py`. No third-party sources downloaded.

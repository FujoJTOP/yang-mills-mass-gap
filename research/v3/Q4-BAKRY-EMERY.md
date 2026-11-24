# Q4-BAKRY-EMERY — a local/weighted Bakry–Émery criterion on the gauge group + a defect term for H1: `Ric_{SU(N)}=(N/4)g` is positive but `Ric+Hess S_W<0` (global CD dead = H6), and the *weighted* version is a **capacity** (bottleneck) statement, not a measure statement; the missing 4D estimate is the Savvidy **gate density** `ρ(barrier)>0` (single-well), which is H1 renamed

> 2026-11-24 · Attack surface **Q4** on `H1`: `χ = Var_μ(ξ_{k→0}) ≤ 1/λ` uniformly in `(L,a)` at
> `β=2.5`, equivalently `G″_{k=0} ≥ λ > 0`. Task: a **local/weighted Bakry–Émery criterion** on the
> gauge group plus a defect term — compute the `Ric` of `SU(N)` (heat-kernel/Haar metric) and the
> Hessian of the Wilson action, write the `Γ₂` criterion for the `k=0` marginal; since the global
> criterion fails (H6), formulate a weighted/local criterion on the good set with the defect
> controlled by its measure or by the measured `D/S≈2%`, try to prove a weighted Poincaré/LSI with
> constant `C_good(1+CD)`, and convert it into `χ≤1/λ`; test on 1D `SU(2)` (finite `χ`, ideally
> `≈0.22`) and free `U(1)` (must fail); state exactly what is missing for 4D at `β=2.5`.
> Reads: `H6-HEAD-SPECTRAL.md`; `P-WAVE-SYNTHESIS.md`, `P1`/`P2`/`P4`; `PROOF-JACOBIAN-SAVES.md`;
> `SM3-ANISOTROPIC-RG.md`; `unconditional/W6-BE-at-2.5.md`, `S1-CONDITIONAL-POINCARE.md`,
> `INVENT-VARIANCE-INEQ.md` (SCBL `D/S≈2%`); `CHEEGER-CORE`/`CHEEGER-BUILD`; `docs/paper-v2/ch02-squeeze.tex`;
> `docs/paper-v3/main.tex` §4,§6.
> Witness: `research/v3/w7/q4_bakry_emery.py` (`py -3.12 -I`, numpy+scipy, deterministic; output §5).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## §0 · Verdict (no QED, Clay 0%)

The Bakry–Émery angle splits cleanly into a **true geometry part** (Ricci is positive, the `Γ₂`
criterion is a theorem), a **true failure** (the global weighted Ricci `Ric+Hess S_W` is negative at
`β=2.5`, exactly H6), and a **true but conditional rescue** (the *weighted* criterion) whose defect
term is a **capacity / gate density**, not a measure — and the measured `D/S≈2%` is an averaged
variance ratio that does **not** control that capacity. Precisely:

1. **Ricci of the gauge group `[T]`.** With the bi-invariant (heat-kernel / Haar) metric
   `⟨X,Y⟩=−2tr(XY)`, `SU(N)` is Einstein with
   ```
   Ric_{SU(N)} = (N/4)·g ,   i.e.  Ric(v,v)/|v|² = N/4 > 0   pointwise   (SU2: ½, SU3: ¾).
   ```
   (witness §5·R: from the structure constants, `Σ_{ab,cd} f^{acd}f^{bcd}=Nδ^{ab}`, so
   `Ric(e_a,e_a)=¼Σ_{i,c}(f^{aic})²=N/4`; `SU(2)=½` ↔ round `S³` of radius 2, `κ=1/2`.) A compact
   semisimple group has **positive** Ricci — this is the entire positive input.

2. **Wilson Hessian and the global `Γ₂` failure `[T-standard]`+`[NG]`.** The quadratic Wilson Hessian
   on link fluctuations in a background `F` is the lattice **Savvidy operator** `−D²+2ad F_B`, whose
   most negative eigenvalue is `≈−g|B|≈−1.1β≈−2.8` at `β=2.5` (`S4 §3.2`, `SM3`, `W6`). The
   `Γ₂` criterion for the gauge-measure `dμ∝e^{−S_W}dHaar` is `CD(ρ,∞)`: `Ric+Hess S_W ⪰ ρg`. At
   `β=2.5`
   ```
   Ric + Hess S_W ≈ +0.50 + (−2.80) = −2.30  (SU2) ,   +0.75 + (−2.80) = −2.05  (SU3)   [NG]  (§5·G)
   ```
   — **negative**, reproducing `W6-BE-at-2.5` (`−2.3/−2.0`). The global criterion **fails**; this is
   H6's *essential-infimum* obstruction in Ricci coordinates (`ess inf(Ric+Hess S_W)<0`).

3. **The `k=0` marginal: BE is lossless, so the group `Ric` is silent on `G″` `[T]`.** Pushing the
   full `Γ₂` down the projection `P` onto the coarse mode is *lossless* on flat coarse space
   (`OTHER-H1-POINCARE`, `B3-FUNCTIONAL-BE`, H6 reading (ii)): the marginal generator
   `L_ν=∂_a²−G″∂_a` has `Γ₂^ν(f)=|f″|²+G″f′²`, so the **marginal BE constant is literally
   `ess inf G″`** (H6-B). Hence *no* property of `Ric_{SU(N)}` or of the fast Hessian spectrum can
   bound the `k=0` `χ` beyond `ess inf G″`: the group curvature is integrated out into `ν` and
   reappears only as the marginal curvature `G″`. **The gauge-group BE cannot manufacture the coarse
   gap** — consistent with the "fast fiber is silent on `G″`" two-coordinate obstruction
   (`B3-ISOPERIMETRIC`, H6 reading (iii)).

4. **The weighted criterion — exact form, and why the defect is a CAPACITY not a MEASURE `[T]`+`[C]`.**
   The global `ess inf` is negative only on the Savvidy tail; a *weighted* criterion asks whether `ν`
   restricted to the good set `𝒢={G″≥ρ_good}` controls `λ₁(ν)`. In 1-D (the `k=0` marginal **is** a
   1-D measure) the half-line is the optimal Cheeger set, so the Cheeger constant is **exact**:
   ```
   h(ν) = inf_x  2ρ(x)/min(F(x),1−F(x)) ,   F = CDF ,     and    λ₁(ν) ≥ h(ν)²/4 .          (4.1a)
   ```
   The bottleneck is the **gate** `x₀` where `ψ(x)=2ρ(x)/min(F,1−F)` is smallest. Two exact regimes:
   ```
   single well (no gate):  ψ ≈ 2ρ(median)/½ ≈ 4ρ(median) ~ √ρ_good  ⇒  λ₁ ≥ h²/4 ~ ρ_good   [T]
   interior barrier:       ψ(gate) = 2ρ(gate)/min(ν_𝒢,ν_𝒟) ,  and the exact Kramers asymptotic
                           λ₁(ν) ≍ ρ(gate)·√(V″_min V″_sad)/…  ∝ ρ(gate)   (double well)   [T-std]  (4.1b)
   ```
   where `ρ(gate)=ρ(barrier)=e^{−V_eff(barrier)}/Z` is the **marginal density at the Savvidy barrier**
   (the Kramers factor `e^{−ΔV}`). **`λ₁` is controlled by the gate density `ρ(gate)`, not by the tail
   measure `ν(𝒟)`.** The witness (§5·W) makes it sharp: for `V=(g/2)(x²−1)²`, as `g:0.5→16`, the
   defect **measure** `ν(G″<0)` *shrinks* `0.31→0.0076` (41×), yet the exact Cheeger `h` *collapses*
   `1.02→0.0021` (485×), `λ₁` collapses `0.79→0.0046` (170×), `λ₁/G″_{well}` collapses `0.79→10⁻⁴`
   (matching H6 §3's Kramers ratio `→0`), and `CD:=ρ_good/h²` grows `~1→7×10⁶`. **Defect-measure
   is the wrong control; the gate density `ρ(gate)` (which sets `h`) is the right one.** The requested
   `C_good(1+CD)` form is the **rigorous** Cheeger bound `λ₁ ≥ h²/4 = ρ_good/(4CD)` (`CD=ρ_good/h²`),
   good in the single well (`CD=O(1)`) and degrading exactly when a gate appears; the *sharper* single-
   well statement `λ₁≥cρ_good` is H6-B, while at a barrier only `λ₁≥h²/4` and the Kramers scaling
   `λ₁∝ρ(gate)` survive (the criterion is then a scaling, not a proof).

5. **The measured `D/S≈2%` does not transfer `[V]`.** The `D/S≈2%` (`INVENT-VARIANCE-INEQ`, SCBL with
   `S=E_μ[H]` the *mean Hessian* and `D=Var_μ(⟨∇S,v⟩|0)` the *force variance*) is an **averaged**
   variance-deficit: `E_ρ[G″]=E[S]−E[D]≈0.98 E[S]>0`. It bounds a *mean*, not the **gate density**
   `ρ(barrier)` nor the **mass in the barrier neighborhood** `ν(𝒟)`. It is consistent with — but does
   not prove — the single-well hypothesis that the weighted criterion needs. (It is a different `D`
   from H2's rank-one `D/S≈0.97`; do not conflate.)

6. **1D `SU(2)` — the weighted criterion reproduces a finite `χ` `[T]`+`[V]`.** The exact one-link
   Haar marginal `ρ₁(z)=π(Cin2π−Cin|z|)` is **log-concave**: `ess inf G″=+0.205>0`, `G″(0)=0.2051`,
   `Var=3.7198`, and the generator's spectral gap `λ₁(ν)=0.2645` (witness §5·S), so
   `χ≤1/λ₁=3.78` (measured `Var=3.72≤3.78`). At `β=2.5` the periodic `k=0` floor is `G″→0.222`,
   `Var≈4.45`, so `χ≤1/0.222≈4.5` — reproducing P-wave's saturation `c=0.222`, `χ≈4.45`. The `≈0.22`
   in the brief is the curvature/`λ` floor; `χ≈4.5`. The 1D object is a **single well**, which is
   exactly why plain/weighted BE is non-vacuous there.

7. **Free `U(1)` — the criterion returns `0` `[T]`+`[V]`.** `Ric_{U(1)}=0` (abelian), `Φ=0` (flat
   Haar); the massless `k→0` marginal is flat, `G″=(2π/L)²→0`, `λ₁→0` (witness §5·U). **Acid test
   passed**: no floor is manufactured on a free field.

8. **What is missing for 4D at `β=2.5` `[O]`.** The weighted criterion `(2.1a)` gives a finite,
   uniform `χ` **iff** the `k=0` coarse marginal is a **single well** with a **gate density bounded
   below**, `ρ(gate)≥c>0`, uniformly in `(L,a)` — equivalently, the Savvidy negative directions do
   *not* form an interior bottleneck between two `ν`-massive wells. This is **exactly H1**, now in
   geometric form: *no interior barrier below the wells, and `ρ(barrier)` not exponentially small.*
   The measured `G″(a₀*)>0`, `D/S≈2%`, and the 1D single-well exact result are **evidence** for it,
   not proof. **The missing estimate is one number: a uniform lower bound on the marginal density at
   the Savvidy gate (equivalently an upper bound on the barrier free-energy `V_eff(barrier)` relative
   to the well).** Clay 0%.

---

## §1 · Geometry: `Ric` of `SU(N)`, the Wilson Hessian, and the `Γ₂` criterion

### 1.1 Ricci of the gauge group `[T]`

Use the flat chart `U=e^{X}`, `X=X^a iT^a∈su(N)` (anti-Hermitian generators `e_a=iT^a`,
`tr(T^aT^b)=½δ^{ab}`), and the bi-invariant metric
```
⟨X,Y⟩ = −2 tr(XY)  ⇒  {e_a} orthonormal (|e_a|²=1),  Haar = Riemannian volume.
```
For a bi-invariant metric the sectional curvature on orthonormal `X,Y` is `K(X,Y)=¼|[X,Y]|²`
(the normalisation fixed by `SU(2)`=round `S³`), so
```
Ric(e_a,e_a) = Σ_i K(e_a,e_i) = ¼ Σ_{i,c} (f^{aic})² = ¼·(Casimir sum) = N/4 ,            (1.1)
```
using the identity `Σ_{i,c}f^{aic}f^{bic}=Nδ^{ab}` (verified to `10⁻¹⁵`, witness §5·R). Hence
`Ric=(N/4)g`, Einstein, **positive** — `SU(2):½`, `SU(3):¾`. (This agrees with `W6`.)

### 1.2 The Wilson Hessian = the Savvidy operator `[T-standard]`

`S_W=βΣ_□(1−N^{-1}Re trU_□)`. Quadratic expansion about a background with field strength `F_B` gives
```
d²S_W[ω,ω] = β Σ_□ |d_A ω|² + 2⟨F_B,[ω,ω]⟩ + O(ω³)   ⇒   −D² + 2 ad F_B   (per link, soft limit),
```
the lattice **Savvidy operator** (`S1 §4.2`, `S4 §3.2`, `PROVE-Gpp`). Its spectrum is `[−g|B|, Λ_UV]`:
the **negative mode** `−g|B|≈−1.1β` is the transverse magnetic mode, gauge-physical after the quotient
(`W4`, `CLAY.md §5`). This is "the wall": the one place a positive `Ric` can be overwhelmed.

### 1.3 The `Γ₂` criterion, global and marginal `[T]`

The measure `dμ=Z^{-1}e^{−S_W}dHaar` is reversible for `L=Δ_{𝒜}−⟨∇S_W,∇⟩`. Bakry–Émery:
```
Γ₂(f,f) = |Hess f|² + (Ric_{𝒜} + Hess S_W)(∇f,∇f) ;   CD(ρ,∞): Γ₂ ≥ ρ Γ  ⟺  Ric_{𝒜}+Hess S_W ⪰ ρg.   (1.2)
```
**Global at `β=2.5`:** `Ric+Hess S_W ≈ −2.30 (SU2) / −2.05 (SU3) < 0` (§0·2, §5·G) — fails on the
Savvidy direction, `[NG]`.

**Marginal (`k=0`):** project along the coarse mode `a=⟨v,X⟩`. On flat coarse space the pushforward is
*exact* (`ℰ_ν(φ)=ℰ_μ(φ∘P)`, `OTHER-H1-POINCARE`), and the `1`-D marginal generator
`L_ν=∂_a²−G″(a)∂_a` has `Γ₂^ν(f)=|f″|²+G″|f′|²`. Therefore
```
CD-marginal(ρ):  G″(a) ≥ ρ  a.e.   ⟺   λ₁(ν) ≥ ρ .      (BE on 1-D is lossless; H6-B.)            (1.3)
```
So the **only** positive constant Bakry–Émery can give for the `k=0` `χ` is `ess inf_a G″(a)` — the
group `Ric` has been absorbed into `ν` and is invisible here. This is the precise sense in which
"the gauge-group `Ric` is silent on `G″`": the marginal BE constant *is* `ess inf G″` (§0·3).

---

## §2 · The weighted / local criterion and the capacity defect

### 2.1 The exact 1-D Cheeger form, and the two regimes `[T]`

The `k=0` marginal `ν` is a probability measure on `ℝ` (density `ρ`, `G″=−(log ρ)″`). In 1-D the
optimal isoperimetric set is a half-line, so the Cheeger constant is **exact** and the Cheeger
inequality is a **rigorous lower bound** on the spectral gap:
```
h(ν) = inf_x ψ(x) ,   ψ(x) := 2ρ(x)/min(F(x),1−F(x)) ,    λ₁(ν) ≥ h(ν)²/4 .                 (2.1a)
```
There are two exact regimes for `h` (and for `λ₁`):

- **Single well (no interior gate).** The infimum of `ψ` is attained in the bulk, `ψ(median)≈2ρ(med)/½
  =4ρ(med)~√ρ_good`, so `h~√ρ_good` and `λ₁ ≥ h²/4 ~ ρ_good`. If in fact `G″≥ρ_good>0` **everywhere**
  (`ess inf G″>0`, globally log-concave) then the sharp Bakry–Émery bound `λ₁≥ρ_good` holds directly
  (H6-B). **Criterion closes.** `[T]`
- **Interior barrier (Savvidy gate).** If `ρ` has an interior minimum `ρ(gate)` between two
  `ν`-massive regions, the gate upper-bounds `h`: `h ≤ ψ(gate)=2ρ(gate)/min(ν_𝒢,ν_𝒟)`. More sharply,
  the transition test function gives the exact **Kramers** asymptotic `λ₁(ν) ≍ ρ(gate)·(prefactor)`,
  i.e. `λ₁` is **proportional to the gate density**, while `G″_{well}=ρ_good` stays positive.
  **Criterion fails** as `ρ(gate)→0`. `[T-std]`+`[V]`

The brief's `C_good(1+CD)` form is realized **rigorously** by `CD := ρ_good/h(ν)²`:
```
λ₁(ν) ≥ h(ν)²/4 = ρ_good/(4CD)   ⇒   χ ≤ 4CD/ρ_good .                                       (2.1b)
```
`CD=ρ_good/h²` is `O(1)` in the single well (`h²~ρ_good`) and blows up at a barrier (`h` collapses);
it is the **curvature-to-Cheeger ratio**, and `ρ_good/ρ(gate)` is its physically suggestive proxy
(`h≤ψ(gate)~ρ(gate)/mass` at a gate). `(2.1b)` is a **rigorous** lower bound (Cheeger); the sharper
`λ₁≥cρ_good` in the single well is H6-B.

### 2.2 Why the defect MEASURE is the wrong control `[T]`+`[V]`

The intuition "the negative set is small ⇒ small defect" is **false** for the spectral constant: a
*single-point* interior barrier (measure zero) already gives `λ₁~e^{−ΔV}` (Kramers; H6 §3). The witness
(§5·W) quantifies it on `V=(g/2)(x²−1)²` (FEM `λ₁`, exact Cheeger `h`, exact gate data):

| `g` | `λ₁`(FEM) | Kramers | `G″_{well}=2g` | `λ₁/G″_{well}` | `ν(G″<0)` | `ρ(gate)` | Cheeger `h` | `h²/4` |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 7.92e-1 | 1.75e-1 | 1.0 | 0.79 | 3.1e-1 | 2.6e-1 | 1.02 | 0.26 |
| 2.0 | 7.49e-1 | 3.31e-1 | 4.0 | 0.19 | 2.7e-1 | 1.9e-1 | 0.75 | 0.14 |
| 8.0 | 1.19e-1 | 6.60e-2 | 16 | 7.4e-3 | 6.3e-2 | 1.9e-2 | 7.7e-2 | 1.5e-3 |
| 16.0 | 4.58e-3 | 2.42e-3 | 32 | 1.4e-4 | 7.6e-3 | 5.2e-4 | 2.1e-3 | 1.1e-6 |

As the barrier grows `g:0.5→16`, the defect **measure** `ν(G″<0)` *shrinks* `0.31→0.0076` (41×), yet
the exact Cheeger `h` *collapses* `1.02→0.0021` (485×), `λ₁` collapses `170×`, and `λ₁/G″_{well}`
collapses to `10⁻⁴` — exactly H6's Kramers ratio `→0`; `CD=2g/h²` grows `~1→7×10⁶`. The rigorous
Cheeger bound `λ₁≥h²/4` holds throughout (e.g. `1.1×10⁻⁶ ≤ 4.6×10⁻³` at `g=16`). A single Gaussian well
(`ν(bad)=0`, `h²≈ρ_good`) has `λ₁=ρ=G″` exactly (BE sharp). **So the control is the gate density
`ρ(gate)` (via `h`), never the defect measure `ν(𝒟)`.**

### 2.3 Conversion to `χ≤1/λ` `[T]`

For any probability measure and any unit coarse mode `a` (`|∇a|=1`, so `ℰ_ν(a)=∫|∇a|²dν=1`), the
Rayleigh quotient gives `λ₁(ν) ≤ ℰ_ν(a)/Var_ν(a)=1/χ`. Hence
```
λ₁(ν) ≥ λ   ⇒   χ = Var_ν(a) ≤ 1/λ .                                                       (2.2)
```
This is the clean, exact direction (D3; paper-v2 §"`G″` versus `1/χ`"): a lower bound on the spectral
gap is a **ceiling** on `χ` (`=H1`). In the single-well regime `(2.1a)` and `(2.2)` give
`χ ≤ (1+CD)/ρ_good`; all the content is in `λ₁(ν)≥λ`, and `(2.2)` is one line.

---

## §3 · Tests

### 3.1 1D `SU(2)` — finite `χ≈4.5`, curvature floor `≈0.22` `[T]`+`[V]`

**(a) Exact anchor (β=0, one link).** The exact marginal `ρ₁(z)=π(Cin2π−Cin|z|)` on `(−2π,2π)`:
`G″(0)=−ρ₁″(0)/ρ₁(0)=1/(2Cin2π)=0.205115` exactly; `Var=3.719824`; `1/Var=0.268830`;
`G″·Var=0.763` (non-Gaussian). The FEM generator gap is `λ₁(ν)=0.264451`, so `χ≤1/λ₁=3.78`
(measured `Var=3.72≤3.78`). The marginal is **log-concave** (`ess inf G″=+0.205>0` at `z=0`), so BE is
non-vacuous directly — the 1D object is a single well (§5·S).

**(b) Periodic `k=0` saturation (β=0, exact convolution).** The `L`-fold convolution `Z₀(k)=[φ(k/√L)]^L`
reproduces P2's exact `β=0` sumrule `G″(0)=0.2396,0.2542,0.2616,0.2652` for `L=2,4,8,16` → `1/Var(Haar)
=0.2688` (CLT), `Var→3.7198` (§5·P). At `β=2.5` the periodic `k=0` floor is `c=0.222` with `Var≈4.45`,
`G″·Var=0.963` (`P2`/`P4`) — same single well, lightly perturbed by the bounded-negative Wilson part
(`−0.045`). The weighted criterion therefore gives
```
β=2.5:  λ₁≈G″→0.222   ⇒   χ ≤ 1/0.222 = 4.50   (measured χ≈4.45).      [T]+[V]
```
The brief's "`≈0.22`" is the curvature/`λ` floor; the susceptibility ceiling it implies is `χ≈4.5`.

### 3.2 Free `U(1)` — the criterion must fail; it does `[T]`+`[V]`

`Ric_{U(1)}=0` (abelian), `Φ=0` (flat Haar ⇒ `W_geo=0`). The massless `k→0` marginal is flat Gaussian:
`G″=(2π/L)²`, `λ₁=G″`; as `L→∞`, `G″=2.467,0.154,0.0096,0.0006→0` for `L=4,16,64,256` (§5·U). So
`λ₁→0`, `χ→∞`. **No floor is manufactured** — the acid test passes: the weighted criterion returns
`0` on the free field, exactly as it must, because every positive input (`Ric`, `Φ″`) vanishes.

---

## §4 · What is missing for 4D at `β=2.5` — the exact statement `[O]`

The weighted criterion, once the defect is correctly a **capacity**, reduces H1 to one geometric
question about the **`k=0` coarse marginal** `ν` at `β=2.5`:

> **(Q4-◆) Single-well + bounded gate `[O]`.** The `k=0` marginal `ν` is (i) **single-well on its
> bulk** — `G″≥ρ_good>0` on the connected `𝒢` that carries mass `≥1−p` — and (ii) **gated** — the
> density at its worst interior dip / outer gate satisfies `ρ(gate)≥c>0` uniformly in `(L,a)` (no
> Savvidy barrier with `e^{−ΔV}`-small density between two `ν`-massive wells). Then the exact Cheeger
> `(2.1a)` gives `λ₁(ν)≥h²/4` with `h` bounded below, and `(2.2)` gives `χ≤4/h²` uniformly.
> **(Q4-◆) = H1.**

**Nothing in the bank proves (Q4-◆) at `β=2.5`.** What exists:
- **1D (`β=0`, exact):** `(Q4-◆)` holds — the marginal is log-concave, single well (witness §5·S).
- **1D (`β=2.5`, P-wave):** holds numerically — `G″→0.222>0` saturating, gapped transfer (`Δ=0.6444`),
  single well; `[V]`, not `[T]`.
- **4D (`β=2.5`):** the measured `G″(a₀*)>0` (typical point), `D/S≈2%` (small *mean* variance deficit),
  and `W6`'s fact that the negative Hessian directions are `90–97%` **orthogonal to the collective
  plane-wave subspace** are **evidence** for single-wellness — but the Savvidy double-well possibility
  (`H6 §3`: `m²<0` at the crossover, stabilised to two wells by compactness) is exactly the
  un-controlled case. The averaged `D/S≈2%` bounds `E[D]`, which is **orthogonal** to both
  `ρ(barrier)` (a *minimum* of the density) and `ν(𝒟)` (the barrier mass).

**The single missing estimate, stated once:**
```
∃ c>0, p<1/2 :  ν_β=2.5({G″<0}) ≤ p   AND   ρ(gate) ≥ c   AND   𝒢 = {G″≥ρ_good} is connected,
uniformly in (L,a).     [O]  = H1 in bottleneck coordinates.
```
Equivalently: **a uniform lower bound on the `k=0` marginal density at the Savvidy gate**, i.e. an
upper bound on the barrier free-energy `V_eff(barrier)` relative to the well depth — the one number
the weighted Bakry–Émery criterion needs and the campaign does not yet have. This is the same object as
H6's `ess inf G″≥λ` on the typical bulk, the Cheeger bottleneck `h(ν)≥h₀` (`CHEEGER-BUILD`,
`LOC-SUFFICIENCY`), and the Kramers barrier `ΔV`; the weighted-BE coordinate makes it a **density
bound** rather than a curvature bound.

**Honest scope.** The Bakry–Émery *method class* was already named a `[NG]` in H6; this note does not
overturn that — it **explains the no-go in Ricci language and isolates the one conditional rescue**.
The rescue is real but demands (Q4-◆) = H1. So the geometry (Ricci, `Γ₂`) is a theorem, the global
failure is a theorem, and the weighted criterion is a *conditional* theorem whose hypothesis is
exactly the gap. **No QED; Clay 0%.**

---

## §5 · Witness — `research/v3/w7/q4_bakry_emery.py` (`py -3.12 -I`, deterministic, numpy+scipy)

Quoted output (abridged):
```
(R) RICCI SU(N), metric <X,Y>=-2tr(XY)
  SU(2): dim= 3  Ric(e_a,e_a)=0.500000 (theory N/4=0.5000)  |sum ff - N I|=8.9e-16
  SU(3): dim= 8  Ric(e_a,e_a)=0.750000 (theory N/4=0.7500)  |sum ff - N I|=1.3e-15
(G) GLOBAL BE beta=2.5
  SU(2): Ric=0.50  Hess S_W(min)~-2.80  =>  -2.30 <0   (W6: -2.3)
  SU(3): Ric=0.75  Hess S_W(min)~-2.80  =>  -2.05 <0   (W6: -2.0)
(S) 1D SU(2) exact one-link Haar marginal
  G''(0)=0.205115 (=1/(2Cin2pi))  Var=3.719824  1/Var=0.268830  G''*Var=0.7630
  ess inf G''=+0.2051 at z=0  lambda_1(nu)=0.264451
  exact 1-D Cheeger h=inf 2rho/min(F,1-F)=0.7759 (gate z=0, F=0.5);  h^2/4=0.1505 <= lambda_1=0.2645 (OK)
  => LOG-CONCAVE (essinf G''=+0.205>0);  G''(0)=0.205 <= lambda_1=0.2645 <= 1/Var=0.2688;  chi<=3.78.
(P) periodic beta=0 k=0 convolution: G''(0)=0.2396,0.2542,0.2616,0.2652 (L=2,4,8,16) -> 1/Var=0.2688
     beta=2.5 periodic (P2/P4): G''=0.222, Var=4.45, G''*Var=0.963
(W) double well V=(g/2)(x^2-1)^2   (FEM lam1, exact Cheeger h, exact gate data)
       g  lam1(FEM)  Kramers  G''w=2g  lam1/G''w  nu(G''<0)  rho_gate   Cheeger h   h^2/4
      0.5  7.921e-01  1.75e-01   1.000   0.7921   3.105e-01 2.56e-01  1.024e+00 2.62e-01
      2.0  7.487e-01  3.31e-01   4.000   0.1872   2.662e-01 1.86e-01  7.456e-01 1.39e-01
      8.0  1.191e-01  6.60e-02  16.000   0.0074   6.326e-02 1.93e-02  7.729e-02 1.49e-03
     16.0  4.575e-03  2.42e-03  32.000   0.0001   7.607e-03 5.21e-04  2.083e-03 1.08e-06
     => nu(bad) shrinks 41x, Cheeger h collapses 485x, lam1 collapses 170x, lam1/G''well->0 (Kramers).
(U) free U(1): G''=lambda_1=(2pi/L)^2 -> 0 for L=4,16,64,256  (ACID PASS)
```
Readings: **(R)** Ricci `=(N/4)g>0` exact to `10⁻¹⁵`; **(G)** global weighted Ricci `≈−2.3/−2.05<0`,
reproducing `W6`; **(S)** the exact 1D `SU(2)` marginal is log-concave (Cheeger `h=0.776`, `λ₁=0.2645≥h²/4`),
`χ≤3.78`; **(P)** the `β=0` convolution reproduces P2 to 4 digits, `1/Var→0.2688` CLT; **(W)** a 41×
shrink in defect *measure* coexists with a 485× Cheeger collapse and a 170× `λ₁` collapse (Kramers) —
measure ≠ capacity; **(U)** the free field gives `0`.

---

## §6 · Status ledger

| item | statement | label |
|---|---|---|
| Ricci of `SU(N)`, bi-invariant metric | `Ric=(N/4)g>0` (SU2 ½, SU3 ¾); from `Σ f f = Nδ` | **[T]** (§5·R) |
| Wilson Hessian = Savvidy operator | `−D²+2ad F_B`, min `≈−1.1β` | **[T-standard]** |
| global `Γ₂`: `Ric+Hess S_W` | `≈−2.30 (SU2)/−2.05 (SU3)` at `β=2.5` | **[NG]** (= H6; §5·G) |
| marginal BE is lossless | `CD-marginal(ρ) ⟺ G″≥ρ ⟺ λ₁≥ρ`; group `Ric` silent on `G″` | **[T]** (H6-B) |
| weighted PI, capacity defect | `λ₁≥h²/4` exact 1-D Cheeger `h=inf 2ρ/min(F,1−F)`; `CD=ρ_good/h²` | **[T]** (Cheeger) |
| defect-MEASURE vs capacity | `ν(G″<0)↓41×` while Cheeger `h↓485×`, `λ₁↓170×` (Kramers) | **[T]**+**[V]** (§5·W) |
| conversion to `χ` | `λ₁≥λ ⟹ χ≤1/λ` (test fn `a`) | **[T]** (2.2) |
| 1D `SU(2)` exact | log-concave, `G″(0)=0.2051`, `λ₁=0.2645`, `χ≤3.78` | **[T]**+**[V]** |
| 1D `SU(2)` `β=2.5` floor | `G″→0.222`, `Var≈4.45`, `χ≤4.5` | **[V]** (P2/P4) |
| free `U(1)` | `Ric=0`, `Φ=0`, `G″=λ₁→0`, `χ→∞` | **[T]**+**[V]** (§5·U) |
| measured `D/S≈2%` (SCBL) | averaged variance-deficit `E[D]/E[S]`, `G″≈0.98E[S]>0` | **[V]**; does **not** bound `ρ(gate)` |
| **(Q4-◆) single-well + `ρ(gate)≥c`** | uniform, `(L,a)`; `=H1` bottleneck form | **[O]** — missing |
| **H1 / mass gap / Clay** | unbuilt; this note renames, not closes | **[O]**, 0% |

---

## §7 · Sharpest line

**The gauge group has positive Ricci (a theorem), but at `β=2.5` the Wilson Hessian's Savvidy negative
mode overwhelms it (`Ric+Hess S_W≈−2.3<0`), so the *global* Bakry–Émery constant is dead — exactly H6's
essential-infimum obstruction. Projecting onto the `k=0` mode is lossless, so the marginal BE constant
*is* `ess inf G″`: the group curvature is silent on `χ`. The *weighted* criterion can close H1, but only
with a **capacity** defect — the exact 1-D Cheeger `h(ν)`, i.e. the density at the Savvidy gate
`ρ(gate)=e^{−ΔV}`, not the defect *measure*: a 41× smaller negative set coexists with a 485× smaller
Cheeger `h` and 170× smaller `λ₁` (Kramers). The measured
`D/S≈2%` bounds a mean, not that gate density. In 1D `SU(2)` the marginal is a single log-concave well,
so the criterion reproduces `χ≤3.8` (`β=2.5`: `χ≈4.5`); on the free `U(1)` it returns `0`. What is
missing for 4D is one geometric number — a uniform lower bound on the marginal density at the Savvidy
gate, i.e. no interior barrier below the wells — which is H1 itself. No QED. Clay 0%.**

---

*Refs (internal).* `H6-HEAD-SPECTRAL.md` (§1 H6-A/B/C, §3 double well/Kramers); `P-WAVE-SYNTHESIS.md`,
`P1-PERIODIC-HOLONOMY.md`, `P2-PERIODIC-K0.md`, `P4-PERIODIC-TRANSFER.md`; `PROOF-JACOBIAN-SAVES.md`
(`Φ″≥N/12`, Prékopa `W_geo≥0`); `SM3-ANISOTROPIC-RG.md` (hard-fibre `γ_hard`, good-set floor);
`unconditional/W6-BE-at-2.5.md` (global `Ric+Hess S≈−2.3`); `S1-CONDITIONAL-POINCARE.md` §4.2 (double
well), `INVENT-VARIANCE-INEQ.md` (SCBL `D/S≈2%`); `CHEEGER-CORE`/`CHEEGER-BUILD`/`CHEEGER-BLIND`
(bottleneck, Kramers); `B3-ISOPERIMETRIC`, `OTHER-H1-POINCARE`, `B3-FUNCTIONAL-BE` (lossless marginal
BE); `D3-SUSCEPTIBILITY-VARIATIONAL` (`λ₁≤1/Var`); `docs/paper-v2/ch02-squeeze.tex` (`G″` vs `1/χ`);
`docs/paper-v3/main.tex` §4 (good set), §6 (`ess inf`, spectral no-go, periodic floor).
*External (standard, used as textbook tools, not re-fetched):* Bakry–Émery 1985; Brascamp–Lieb, Adv.
Math. **22** (1976) 366; Cheeger, *Problems in Analysis* (1970) 195; Bobkov–Ledoux (log-concave
isoperimetry); Bovier–Eckhoff–Gayrard–Klein (Kramers/metastability); Savvidy, Phys. Lett. B **71**
(1977) 133; Nielsen–Olesen, Nucl. Phys. B **136** (1978) 45; Prékopa, Acta Sci. Math. (Szeged) **32**
(1971) 301. Witness: `research/v3/w7/q4_bakry_emery.py`.
*No third-party sources downloaded this session; all external inputs are standard and cited from the
corpus, not re-fetched.*

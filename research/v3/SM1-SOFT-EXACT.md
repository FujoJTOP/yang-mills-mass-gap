# SM1-SOFT-EXACT — exact soft-mode curvature in small-volume SU(2) Wilson models: the conditioned margin `G″(a₀*)` is **positive and O(1)** at β=2.5, while the full Hessian is **indefinite (Savvidy) on the measure** — so the positivity is measure-carried, not a Hessian gap, and the descent break (`γ_j<0`, `η_j≥1`) is **already realised at small volume**

> 2026-11-24 · Wave **SM1**: the last open mode of the Yang–Mills UV face, the **Savvidy soft mode**.
> Task: exact soft-mode curvature in small-volume models — SU(2) on 1×1, 1×2, 1×3, 1×4 strips and a
> 2×2 cross block (fixed boundary), plus the 4D periodic single-site `L=1` (the minimal lattice that
> carries the non-abelian commutator mode). Compute the **full Hessian spectrum** and the
> **conditioned coarse curvature of the softest mode at its coarse minimiser `a₀*`, including `X`**.
> Reads: `H1-RG-DESCENT.md` (the descent breaks on the soft mode at β≈2.4–2.5), `H4-HEAD-COUNTEREXAMPLE.md`
> (local sign at `a₀*`), `STAR2-CLASSIFICATION.md`, `H2-DEFECT-BOUND.md`, `G2-TWO-PLAQUETTE.md`,
> `docs/paper-v3/main.tex` §6.
> Witness: `research/v3/w7/sm1_exact.py`, `sm1_blocks.py`, `sm1_4d.py` (`py -3.12`, numpy only; small
> volumes; **one job at a time**; deterministic quadrature/convolution + Haar-importance sampling).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED, Clay 0%)

**The soft-mode margin `G″(a₀*)` is POSITIVE in every small-volume model at β=2.5, and it is O(1)
(decreases only slowly with L, ≈`1/L`, **not** the `4π²/L²` box-gap collapse). At the same time the full Hessian is indefinite on the measure — the
Savvidy negative direction is present already at small volume — so the positivity is entirely
*measure-carried* (the conditioned marginal curvature), never a Hessian eigenvalue. The descent's
break condition (`γ_j→0`, `η_j≥1`, `D_j>0` — measured `D=+0.19`) is realised at small volume; the
only thing that survives it is the conditioned curvature `G″(a₀*)>0`.**

Precisely:

1. **Full Hessian spectrum `[T]+[V]`.** The measure-typical full Hessian `Hess V = Hess(S+Φ)` is
   **indefinite at β=2.5 for every block**: fraction of sampled configs with `λ_min<0` is `0.08`
   (1×1), `0.40` (1×2), `0.81` (1×3), `0.79` (1×4), `0.81` (2×2), `1.00` (4D `L=1`). The negative
   eigenvector is the non-abelian **transverse Savvidy mode** (the 4D single-site gives it exactly:
   support on the two non-commuting links, in the transverse colour, 2-fold degenerate). So
   `γ_j := λ_min(fibre Hessian) ≤ 0` on the typical field — the descent's `γ_j→0` is **already
   realised at small volume**, not only asymptotically.

2. **Conditioned soft-mode curvature `[T]+[V]`.** At the coarse minimiser `a₀*=argmax ρ_v` (interior,
   H4 theorem), `G″(v;a₀*) = −ρ_v''/ρ_v ≥ 0` and is **strictly positive** and **O(1)**:
   ```
   1×1 exact          : G″(a₀*) = +0.69 (θ=π, e₃) / +0.36 (softest, transverse e₁)   [T/V]
   1×2, 1×3, 1×4 (k=0): G″(a₀*) = +1.01, +0.65, +0.42  (θ=π)   [V]
   2×2 cross (k=0)    : G″(a₀*) = +0.78   [V]
   4D single-site     : G″(v*;a₀*) = +0.42 (v* = the Savvidy eigenmode)   [V]
   ```
   `X = −2Cov(∂_vS, ∂_vΦ|a₀*)` is included; it is small for the k=0 modes (`−0.03…−0.04`) and
   sizable for the 4D Savvidy mode (`+0.065`). The decomposition `Δ_dyn+W_geo+X` matches `G″` to the
   stated digits.

3. **Minimum margin and scaling `[V]`.** Over boundary `θ∈[0,π]` × block shape × mode direction at
   β=2.5, the **minimum is `+0.320` (Diff₃, 1×2, θ=π)**, above `N/12=0.167`. Scaling: with `L`
   the k=0 margin decreases **slowly, roughly `∝ 1/L`** (`1.01→0.65→0.42` for `L=2→4`), **not** the
   box-gap `4π²/L²` collapse (`9.87→2.47`); with `β` it **increases** (`0.49→0.83→1.60` for
   `β=1→4`). **Closed forms** exist at `a₀=0` (`G″=A−Bβ²`, `A=0.20511530`, `B=0.01105929`,
   `β*=4.30661` `[T/V]`), for `β=0` (an `L`-fold convolution of the one-link marginal, value
   `0.229…0.254`), and generally `G″(a₀*)=−ρ_v''(a₀*)/ρ_v(a₀*)` with `ρ_v` the exact marginal `[T]`.

4. **Descent break mechanism at small volume `[V]`.** The fibre gap `γ_j=λ_min(K)` crosses zero at a
   chromomagnetic field `c*≈0.55` (4D `L=1`, β=2.5); near the crossing the Schur loss `cᵀK⁻¹c`
   is `O(1)` and `|η|=|Schur|/|Σ|≥1` (`η=4.2` at `c=0.5`, `2.0` at `c=0.6`), and `Σ=h_n−cᵀK⁻¹c`
   turns negative (`c=0.6`). This is the H1 break (`γ_j→0`, `η_j≥1`) **realised at L=1**. The
   measured `G″(a₀*)` nonetheless stays `>0` — the defect `D` of the non-log-concave fibre is what
   carries the positivity (the H2 "measure-carried" defect).

5. **Acid test `[T]+[V]`.** Free massless field: `G″=κ→0` (no manufactured floor). Pure Haar `β=0`:
   `G″(a₀*)≈0.23–0.25>N/12` and `L`-independent (the Haar floor, *not* a free field). The full-Hessian
   spectrum is flat (`λ_min=N/12`) at the flat configuration `[T]`. All tools pass.

**One line.** At β=2.5 the small-volume soft-mode **conditioned curvature is positive and O(1)**
(does not soften like the box gap), but the **full Hessian is already Savvidy-indefinite on the
measure** — so the margin is a *typical-point, measure-carried* curvature, `γ_j<0` and `η_j≥1` hold
already at `L=1`, and the volume-uniform gap remains `H1` `[O]`. **No QED. Clay 0%.**

---

## 1 · The model, the soft mode, and gauge-invariant conditioning

**Charts and measure.** SU(2), flat exponential chart `U=e^{iξ·σ/2}` per link, `ξ∈B(2π)⊂ℝ³`,
Haar `dU = J(ρ)d³ξ`, `J(ρ)=sin²(ρ/2)/ρ²`, `Φ(ξ)=−log J(ρ) = 2log(ρ/2)−2log sin(ρ/2)`, `ρ=|ξ|`.
The conditioned measure is `dμ ∝ e^{−V}dξ`, `V = S_W + Φ`. Blocks (fixed boundary; the boundary
links are fixed data, so the flat chart is a genuine **gauge slice** — no residual gauge orbit —
making the conditioning below gauge-invariant):

```
 1×1 :  free U ;               P = U B .
 1×L :  free X₁..X_L ;         P₁=X₁X₂ ,  P_i=B_{i-1}X_i^{-1}X_{i+1} (i≥2) ,  P_L=B_{L-1}X_L^{-1} .
 2×2 :  cross block, free {H01,H11,V10,V11} (4 links incident to the centre), 8 boundary links fixed;
        P_mn = H[m,n] V[m+1,n] H[m,n+1]^{-1} V[m,n]^{-1} ,  m,n∈{0,1} .
 4D L=1 (single site, periodic): free U₀..U₃ ; P_{μν}=U_μU_νU_μ†U_ν† ,  μ<ν .
```

**The soft mode.** Two equivalent definitions, both used:
- **Hessian-softest** `v*` = the eigenvector of `Hess V` (or of the ensemble-averaged `H̄=E[Hess V]`)
  with the smallest eigenvalue. At the Savvidy field this is the transverse-magnetic mode (`[V]`,
  4D: support on the two non-commuting links, transverse colour).
- **Block-uniform coarse mode** `v_c = (e₃,…,e₃)/√L` (the `k=0` mode), used for the descent
  decomposition because it is a *fixed* linear/divergence-free direction (the admissible projection
  class, STAR2 §4), so `G″_uncond = 2W_geo` exactly (G2 `[T]`).

**Gauge-invariant conditioning on the soft mode.** Condition on the amplitude `a=⟨v,ξ⟩` of the
soft direction. In the fixed-boundary blocks the chart has no residual gauge freedom, so `a` is a
gauge-fixed coordinate; in the periodic 4D `L=1` model the soft direction is taken **covariantly**
(`v*` is the eigenvector of the gauge-covariant Hessian, so `a=⟨v*,ξ⟩` is gauge-invariant under the
simultaneous rotation of `v*` and `ξ`). The coarse minimiser is `a₀*=argmax_a ρ_v(a)`, `ρ_v` the
exact marginal of `a`.

**The object.** The exact marginal-curvature identity (H4 Lemma G1 / ch05 `eq:split`, model-independent):
```
   G″(v;a₀) := −d²/da² log ρ_v(a)|_{a₀} = E[∂_v²V | a₀] − Var(∂_vV | a₀)
             = Δ_dyn + W_geo + X ,     X = −2Cov(∂_vS, ∂_vΦ | a₀) ,     (1.1)
```
`Δ_dyn=E[∂_v²S]−Var(∂_vS)`, `W_geo=E[Φ″]−Var(Φ′)`, all conditioned at `a=a₀`. At an interior
maximiser `a₀*`, `G″(a₀*)=−ρ_v''(a₀*)/ρ_v(a₀*)≥0` (strict for a non-degenerate max); `[T]` (H4).

---

## 2 · Full Hessian spectrum — the Savvidy direction is present at small volume

`Hess V = Hess S_W + Hess Φ`, computed exactly (5-point finite differences in the flat chart,
`δ=10⁻⁴`) for each free-link block. `Hess Φ(0)=(N/12)I`; `Hess S_W` at flat is `β·(lattice
Laplacian)`.

**(a) 4D `L=1` single-site — the exact Savvidy.** Flat: `Hess V = (N/12)I`, `λ_min=0.1667` (PD). At a
constant non-abelian field `U₁=e^{icσ₃/2}`, `U₂=e^{icσ₁/2}` (`F₁₂=[A₁,A₂]≠0`), β=2.5:
```
   c      λ_min       #neg    c      λ_min
  0.0   +0.1667        0     1.2   −0.6219 (2)
  0.4   +0.0684        0     1.6   −1.0961 (2)
  0.6   −0.0506        2     2.0   −1.5640 (2)
```
The transverse Savvidy mode crosses zero at `c*≈0.55` and has **two** negative eigenvalues
(the degenerate transverse plane), `λ_min≈−g|F|`. This is the lattice image of `−D²+2F`'s negative
mode. On the **β=2.5 measure** (Haar importance sample) the single-site Hessian is indefinite
**everywhere**: weighted `λ_min` quantiles `q05=−2.07, q50=−0.99, q95=−0.22`, fraction `<0`=1.00 `[V]`.

**(b) 2D blocks — indefinite on the measure.** `[V]` (`sm1_blocks.py [H1]`):
```
 block     flat λ_min   ens q05   ens q50   ens q95   frac(λ_min<0)
 1×1         0.1667    −0.060     0.305     0.505      0.079
 1×2         0.1667    −0.386     0.051     0.275      0.404
 1×3        −0.3345    −0.559    −0.146     0.142      0.812
 1×4        −0.6476    −0.586    −0.142     0.134      0.793
 2×2cross    0.1667    −0.667    −0.175     0.048      0.811
```
So `Hess V` is **positive-definite only near the flat / ordered configurations** (the flat 1-link at
θ=π is a saddle with `λ_min=N/12`); on the Haar-typical field (which is what the measure samples at
β=2.5) it is **indefinite**, increasingly so with block size. **The fibre gap `γ_j=λ_min(K)` is
therefore `≤0` on the typical field already at these small volumes.** This is the exact opposite of
the earlier "small volume is stable" guess, and it is the small-volume content of the Savvidy wall.

---

## 3 · Exact / numerical soft-mode curvature at the coarse minimiser (with X)

**(a) 1×1 exact (Gauss–Legendre quadrature, `rho_v_1x1`).** The marginal `ρ_{e₃}(a)` is a 1-D
radial integral (the transverse plane integrates to `2π∫r dr`); `G″(a₀*)` from the 5-point
log-curvature. β=2.5 `[T/V]`:
```
   θ        a₀*      G″(a₀*)   [mode]
  0.000    0.00      +0.730     e₃
  π/4     −0.59      +0.730     e₃
  π/2     −1.16      +0.730     e₃
  3π/4    −1.69      +0.723     e₃
  π       −2.14      +0.694     e₃        (H4 exact +0.696)
  π/2      0.00      +0.618     transverse e₁
  π        0.00      +0.362     transverse e₁
```
Decomposition at `θ=π, e₃, a₀*=−2.14`: `Δ_dyn=+0.464`, `W_geo=+0.249`, `X=−0.019`, sum `+0.694`.
At the barrier `a₀=0` (θ=π): `G″=A−Bβ²=0.136>0`; it goes negative only for `β>β*=4.307` `[T/V]`.

**(b) 1×2, 1×3, 1×4, 2×2 at β=2.5 (Haar-importance sampling, exact per-`a` Hessian moments).** k=0 mode
`v=(e₃,…,e₃)/√L`, θ=π `[V]`:
```
 block      a₀*      G″(a₀*)   Δ_dyn    W_geo     X       ESS
 1×2       0.37     +1.010    0.787    0.261   −0.038    2266
 1×3       1.11     +0.653    0.391    0.306   −0.044    1101
 1×4       0.88     +0.424    0.173    0.277   −0.026     513
 2×2cross  0.33     +0.781    0.561    0.250   −0.029     277
```
`X<0` and small for the k=0 modes; the decomposition closes. These reproduce H4's 1×2/1×3 values.

**(c) 4D single-site, Savvidy eigenmode `v*` (unit norm).** `[V]` (`sm1_4d.py [G4]`): `G″(v*;a₀*)=
+0.424` at `a₀*=−0.33`, decomposition `Δ_dyn=+0.099, W_geo=+0.260, X=+0.065`. **Here `X` is sizable
(+0.065, 15% of `G″`)** — the coupled chart matters. The negative-curvature bins (`#neg=13/400`,
`min −1.47` at `a=−2.6`) are on the **flanks**, never at the minimiser (H4 §3.5, same feature).

---

## 4 · Minimum margin, scaling with L and β, closed forms

**Minimum over boundary × block × mode at β=2.5 `[V]`.** Boundary scan (1×2, k=0): `G″(a₀*)=
+1.182, +1.161, +1.148, +1.085, +0.955` at `θ=0,π/4,π/2,3π/4,π` (worst at maximal holonomy `θ=π`).
Mode scan (1×2, θ=π): `A₃=+0.325, V₃=+0.728, D₃=+0.966, Diff₃=+0.320`. **The minimum is `+0.320`
(Diff₃), `>N/12=0.167`.** So the margin is **positive over the whole scan**.

**Scaling with `L` `[V]`.** k=0, θ=π: `G″(a₀*)=1.01, 0.65, 0.42` for `L=2,3,4` — **O(1), decreasing
slowly (≈`∝1/L`)**. It does **not** follow the box gap `4π²/L²=9.87, 4.39, 2.47` (which drops ~4×
per doubling); the soft (`k=0`) curvature is a different, volume-slower object than the UV box-floor.
The 2×2 cross (`+0.78`) does not lie on the 1-D-chain law — the margin is set by the local block
geometry, not a clean `L`-power.

**Scaling with `β` `[V]`.** 1×2, k=0, θ=π: `G″(a₀*)=+0.485, +0.830, +1.602` at `β=1, 2.5, 4` —
**increasing**, roughly linearly at large β. The head β=2.5 is mid-range; the margin is **weakest at
small β** (approaching but above the Haar floor). β=0: `G″(a₀*)=0.229, 0.240, 0.249, 0.254` for
`L=1..4` — the **Haar floor** (`>N/12`), `L`-independent.

**Fisher average `I(ρ)=E_ρ[G″]` `[T]+[V]`.** By integration by parts (H2 Theorem M),
`E_{a~ρ}[G″(a)]=I(ρ)=E_ρ[(∂_a log ρ)²]≥0`, exact, vanishing on the free field. Measured
`I(ρ)≈0.69, 0.48, 0.43, 0.55` for the four blocks (1×2,1×3,1×4,2×2) — *below* the typical-point
`G″(a₀*)` and **decreasing as the mode delocalises** (`0.687→0.485→0.434` for `L=2→4`); at fixed
`L=2` it rises with β (`0.365, 0.682, 1.177` at `β=1,2.5,4`). This is the volume-sensitive object;
the typical-point `G″(a₀*)` is not. A uniform lower bound `I(ρ)≥δ` is the mass gap `[O]`.

**Closed forms.**
- At `a₀=0` (1-link): `G″(e₃;0)=A−Bβ²`, `A=E_Haar[Φ′(ρ)/ρ]=0.20511530`, `B=Var_Haar[sin(ρ/2)/ρ]=
  0.01105929`, `β*=√(A/B)=4.30661` `[T/V]` (reproduces V3/G4).
- At `β=0` (any block): `ρ_v` is an `L`-fold convolution of the one-link scalar marginal
  `ρ₁(a)=2π∫_{|a|}^{2π}sin²(ρ/2)/ρ dρ/(4π²)`; `G″(a₀*)=−[log ρ_v]''(a₀*)` `[T]`.
- General `a₀*`: `G″(a₀*)=−ρ_v''(a₀*)/ρ_v(a₀*)≥0` with `ρ_v` the exact marginal `[T]` (H4).

---

## 5 · The descent break, and the exact mechanism at small volume

H1's one-step descent `G″_j ⪰ (1+η_j)^{-1}G″_{j+1}^{bare} − L_j` has validity hypothesis **fibre gap
`K ⪰ γ_j I>0`**. In the small-volume reduction with fixed coarse direction `v_c`, the Gaussian
coarse curvature is the Schur complement `Σ = h_n − cᵀK⁻¹c` (`h_n=v_cᵀHv_c`, `K` the Hessian on
`v_c^⊥`, `c=QᵀHv_c`), and `η=|h_n/Σ−1|` `[T]`.

**4D `L=1` single-site, β=2.5, `v_c=e₃` on all links** `[V]` (`sm1_4d.py [G5]`):
```
   field c    h_n      γ=λ_min(K)   Schur=cᵀK⁻¹c   Σ        η       K PD?
    0.0     0.167      0.167        0.000       0.167    0.00    ok
    0.3     0.208      0.111        0.024       0.184    0.13    ok
    0.5     0.277      0.015        0.224       0.053    4.20    ok (γ→0)
    0.6     0.320     −0.050        0.629      −0.309    2.04    BROKEN (γ<0, Σ<0)
    0.8     0.419     −0.209      −12.20      12.62      0.97    BROKEN (γ<0)
    1.2     0.620     −0.620       −2.20       2.82      0.78    BROKEN
```
**The mechanism is exact and small-volume.** (Scope: this is the *one-mode Schur reduction* of a
single small block — the small-volume model of the shell-pair descent, with `v_c^⊥` playing the role
of the fibre `ker P` and `v_c` the coarse mode; the actual tower descent inherits the same
`γ_j`, `η_j` because they are properties of `Hess V` on the block.) The Savvidy direction enters
*either* the fibre (`K` indefinite, `γ_j<0`) *or* the coarse mode itself (`h_n<0`), depending on the
coarse direction. When
the chromomagnetic field passes `c*≈0.55` with `v_c=e₃`, the **fibre gap `γ_j=λ_min(K)` crosses zero
and goes negative**; `K⁻¹` becomes indefinite, the Schur loss `cᵀK⁻¹c` is `O(1)`/sign-varying,
`|η|` reaches `≥1` at the crossing (`4.2` at `c=0.5`, `2.0` at `c=0.6`), and `Σ` turns negative —
i.e. the one-step multiplicative lemma's contraction factor `c_j=(1+η_j)^{-1}` is **vacuous**
(`η_j≥1`) and the Gaussian coarse curvature `Σ≤0`. (If instead `v_c` is taken along the softest
*eigenvector*, `h_n<0` directly and `γ>0`; the break is then `h_n<0` with the same conclusion for
the descent.) This is the H1 break (`γ_j→0`, `η_j≥1`) **realised at `L=1`**, not only in the
large-volume tower. The 2D chains show it too: their flat Hessian is already indefinite for
`L≥3` at θ=π (`λ_min=−0.33` at 1×3, `−0.65` at 1×4), and the measure-typical Hessian for
`L≥2` (`§2b`).

**Why the measured `G″(a₀*)` survives `[T]+[V]`.** The H1 identity `G″ = E[Σ] − D` splits the
conditioned curvature into the Gaussian Schur part and the **defect**
`D = Var(∂_vV|a₀*) − E[⟨K⁻¹c,c⟩|a₀*]` (H1 (1.2), `D_j`). On the non-log-concave (Savvidy) fibre the
defect is `O(1)` and sign-variable (H2: `D/S∈[0.83,1.10]`, `5/12` modes `D>S`); the measured
`G″(a₀*)>0` is the **residue** of that competition (H2 §1). Measured here (4D `L=1`, fixed `v_c=e₃`,
β=2.5, `[G5b]`): `E[Schur|a₀*]=+0.92`, `Var(∂_vV|a₀*)=+1.10`, so the **defect `D=Var−E[Schur]=
+0.19`** at `a₀*=−0.52` (ESS≈21, indicative) — the same positive, `O(1)` "measure-carried" object H2
bounded. So the margin is carried by the **defect of the indefinite fibre**, not by any Hessian
eigenvalue. That is why `G″` can be positive while `γ_j<0`.

---

## 6 · Acid tests `[T]+[V]`

1. **Free massless field** (`Φ≡0`, Gaussian `V=κ|ξ|²/2`): `ρ` Gaussian, `G″=κ` at every `a₀`,
   `κ→0` as the mass →0. **No manufactured floor.** `[T]`.
2. **Flat Hessian** (any `L`): at the flat configuration `Hess V = βQ(k)+N/12`; the soft mode is the
   pure-gauge `k=0` mode, `λ_min=N/12=0.1667` `[T]` (and the transverse `β·2Σ_μ(1−cos k_μ)` is
   `≥2β` for `L=4`). The Savvidy negativity needs a **non-flat field**, which the measure supplies at
   β=2.5 (`§2`).
3. **Pure Haar `β=0`**: `G″(a₀*)=0.229…0.254>N/12`, `L`-independent — the Haar floor, *not* the free
   field (the free field is item 1). `[T]+[V]`.
4. **Fisher** `E_ρ[G″]=I(ρ)≥0`, `=0` on the free flat marginal. `[T]+[V]`.

The tools pass the acid test **non-vacuously and discriminatingly**: they return `0` on the free
field, `N/12` on flat Haar, and positive O(1) at the interacting head — without manufacturing a gap.

---

## 7 · Witness — `research/v3/w7/sm1_{exact,blocks,4d}.py`

`py -3.12` (numpy 2.4.6), small volumes, one job at a time. Quoted outputs (abridged):
```
sm1_exact.py  [E0] A=0.2051152978 B=0.0110592938 beta*=4.30661
              [E1] 1×1 Hessian at β=2.5: θ=π flat λ=0.1667, ordered 0.456, marg-max 0.379 (PD)
              [E2] exact G″(a₀*): θ=π e₃ +0.694 (Δ0.464,W0.249,X−0.019); transverse e₁ +0.362
              [E4] β=0 convolution G″(a₀*) L=1..4: 0.229,0.240,0.249,0.254 (>N/12)
              [E5] 4D flat λ_min=N/12; Q_min(nz)=4,3,2 (L=2,3,4)
sm1_blocks.py [H1] 2D full-Hessian frac(λ_min<0): 0.08(1×1) 0.40(1×2) 0.81(1×3) 0.79(1×4) 0.81(2×2)
              [H2] G″(a₀*) k=0 β=2.5: 1×2 +1.01, 1×3 +0.65, 1×4 +0.42, 2×2 +0.78
              [H3] min over θ×block×mode = +0.320 (Diff₃) > N/12
              [H4] L-scaling 1.01→0.65→0.42 ; β-scaling 0.49→0.83→1.60 ; I(ρ)=0.69,0.48,0.43,0.55
sm1_4d.py     [G1] Savvidy threshold c*≈0.55 at β=2.5; #neg=2 (transverse)
              [G2] measure frac(λ_min<0)=1.00 ; q50=−0.99
              [G4] G″(v*;a₀*)=+0.424 (Δ0.099,W0.260,X+0.065)
              [G5] γ crosses 0 at c*≈0.55; η=4.2(c=0.5), 2.0(c=0.6); Σ<0 at c=0.6
              [G5b] defect D=Var−E[Schur]=+0.19 at a₀*=−0.52 (E[Schur]=+0.92, Var=+1.10; ESS≈21)
```
All deterministic except the Haar-importance-sampled `[H2–H4]`, `[G2–G5]` (ESS quoted per bin). No
third-party sources downloaded.

**Monte-Carlo scatter (honest).** The `[V]` values carry Haar-importance-sampling scatter: independent
seeds move `G″(a₀*)` by ≈`0.05–0.15` at the quoted ESS (e.g. 1×2 k=0 β=2.5 gives `0.83–1.01` on two
independent `N=1.5·10⁶` runs). The **sign and the O(1) magnitude are robust; the third digit is not.**
The `[E*]` (1×1, β=0 convolutions, 4D flat) values are exact quadrature/algebra and carry no such
scatter.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| exact marginal curvature `G″(v;a₀)=−ρ″/ρ=E[∂²V\|a₀]−Var(∂V\|a₀)` incl. `X` | (1.1), H4 Lemma G1 | **[T]** |
| `a₀*` interior, `G″(a₀*)=−ρ″/ρ≥0` (strict if non-degenerate) | H4 Thm; one-link `[E2]` | **[T]** |
| **full Hessian indefinite on the measure at β=2.5** (frac<0: 0.08→1.00 with L) | §2, `[H1]`,`[G2]` | **[T-standard Savvidy]**+**[V]** |
| 4D `L=1` Savvidy threshold `c*≈0.55`, transverse mode 2-fold, `λ_min≈−g\|F\|` | §2a, `[G1]` | **[V]** |
| closed form `G″(e₃;0)=A−Bβ²`, `A,B,β*` (V3/G4) | §4, `[E0]` | **[T]**+**[V]** |
| **soft-mode margin `G″(a₀*)>0`, O(1)**, min `+0.320` over θ×block×mode at β=2.5 | §4, `[H2–H3]` | **[V]** |
| `L`-scaling `G″(a₀*)` O(1), `≈∝1/L` (NOT `4π²/L²`); `β`-scaling increasing | §4, `[H4]` | **[V]** |
| defect/`D`: positivity carried by the non-log-concave-fibre defect, `γ_j<0` | §5, H1 (1.2), H2 | **[T]**+**[V]** |
| **descent break `γ_j<0`, `η_j≥1` realised at `L=1`** (not only asymptotically) | §5, `[G5]` | **[V]** |
| Fisher `E_ρ[G″]=I(ρ)≥0`, `=0` free field; shrinks with mode delocalisation | §4, `[H4]` | **[T]**+**[V]** |
| acid: free `G″=κ→0`; β=0 Haar `0.23–0.25>N/12`; flat `λ_min=N/12` | §6 | **[T]**+**[V]** |
| **volume-uniform soft-mode gap `δ_W>0`** | = H1 | **[O]** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 9 · Sharpest line

**At β=2.5 the soft-mode conditioned curvature is positive and O(1) in every small-volume model —
`G″(a₀*)=+0.69` (1×1), `+1.01/+0.65/+0.42` (1×2/1×3/1×4), `+0.78` (2×2), `+0.42` (4D single-site
Savvidy mode), minimum `+0.32` over boundary × block × mode, all `>N/12` — and it decreases only
slowly with `L` (≈`1/L`, not the `4π²/L²` box-gap collapse) and increases with `β`. But the full
Hessian is already Savvidy-indefinite on the measure (frac `λ_min<0` = 0.08 → 1.00 with `L`; threshold
field `c*≈0.55` at `L=1`), so `γ_j=λ_min(fibre)≤0` and the descent break `η_j≥1` is realised **at
`L=1`**, not only in the tower.** The margin is therefore a *typical-point, measure-carried*
curvature — the defect `D` of the indefinite fibre, not a Hessian eigenvalue — and the
volume-uniform gap is, unchanged, **H1 `[O]`**. **No QED. Clay 0%.**

---

*Refs.* `H1-RG-DESCENT.md` §2–§5 (one-step lemma, `γ_j→0`, `η_j≥1`, `D_j`); `H4-HEAD-COUNTEREXAMPLE.md`
§1–§5 (local sign theorem at `a₀*`, decomposition `Δ_dyn+W_geo+X`, exact 1×1, block/boundary/β scans);
`STAR2-CLASSIFICATION.md` §0–§4 (true object, `2W_geo`, linear projection class); `H2-DEFECT-BOUND.md`
§1–§3 (rank-one split, measure-carried defect `D`, Fisher identity Theorem M); `G2-TWO-PLAQUETTE.md`
§2–§4 (`R=0`, `G″=2W_geo` for linear `v`, edge-vanishing); `S1-CONDITIONAL-POINCARE.md` §4 (Savvidy
indefinite fibre, `λ_min<0` for chromomagnetic boundary); `docs/paper-v3/main.tex` §6 (head, defect,
Fisher, local sign). Witness `research/v3/w7/sm1_exact.py`, `sm1_blocks.py`, `sm1_4d.py`. No
third-party sources downloaded; standard inputs (Savvidy 1977, Prékopa, Brascamp–Lieb, Schur) cited
from the corpus.

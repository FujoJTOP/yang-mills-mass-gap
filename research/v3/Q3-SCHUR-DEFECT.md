# Q3-SCHUR-DEFECT — exact elimination of the unstable mode and the Schur-complement bound: the inertia lemma `Σ>0 ⇔ n₋(H)=n₋(K)`, the one-plaquette dichotomy (**fibre-confined ⇒ excisable**, **collective ⇒ irreducible**), the strip's positive constant-mode curvature, and the rank-one defect lemma `χ ≤ 1/ρ_rest` whose binding hypothesis is a **uniform fibre gap**, not the measured `D/S`

> 2026-11-24 · Attack angle **Q3** on `H1` — the last open input: `χ = Var_μ(ξ_{k→0}) ≤ 1/λ`
> uniformly in `(L,a)` at `β=2.5`, equivalently `G″_{k=0} ≥ λ>0`. The obstruction is the Savvidy
> **negative mode**: `Hess V` is indefinite, `ess inf G″<0`, which kills global log-concavity; the
> unstable mode is rank-one / low-rank. Task: integrate the unstable mode out **exactly**, compute
> the Schur complement on the remaining modes, and decide whether a rank-one defect lemma
> `χ ≤ (1+η)/ρ_rest` with `η` from the measured defect `D/S≈2%` can hold.
> Reads: `SM2-RANK-ONE.md` (margin identity `G″=(S−D_Sav)+|D₊|`, ceiling `m_soft≤G″`, overlap lemma),
> `V2-TRUE-BOUND.md` (Schur monotonicity `Σ⪰ρ`, Lemma V2.1, `G″=E[Σ]−D`), `H1-RG-DESCENT.md`
> (Schur loss `⟨K⁻¹c,c⟩=g²/γ_j→∞`, `η_j≥1`), `H6-HEAD-SPECTRAL.md` (typical point vs essential
> infimum), `S3-UV-FESHBACH.md` (`G″=E[Σ]−D`), `P-WAVE-SYNTHESIS.md` (periodic floor), `G4-GENERAL-BOUND.md`
> (exact `G″=A−Bβ²` counterexample), `docs/paper-v3/main.tex` §6.5–6.7.
> Witness: `research/v3/w7/q3_schur_defect.py` (`py -3.12 -I`, numpy only; exact quadrature / exact
> linear algebra; output quoted §5).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[O]** exact
> missing step · **[NG]** no-go for a named class · **[V]** this-round numerics. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

The Schur complement is the right object, and it splits the problem into two inequivalent regimes
that the corpus had conflated. Five exact results decide the task:

1. **Inertia lemma `[T]`.** For `H = [[α,c^ᵀ],[c,K]]` with `K` invertible, `Σ := α − c^ᵀ K⁻¹c`
   satisfies `det H = det K · Σ` and **Haynsworth inertia additivity**
   ```
   In(H) = In(K) ⊕ In(Σ) ,   i.e.   n₋(H) = n₋(K) + n₋(Σ).
   ```
   Hence `Σ>0` (the coarse curvature is positive after eliminating the fibre) **iff** every negative
   mode of `H` already sits in the fibre block `K`. This is the exact answer to "does eliminating the
   unstable mode restore positivity": **it does iff the unstable mode is fibre-confined.** When `K≻0`
   yet `Σ<0`, the negativity is **irreducible** — it lives in the coarse direction, and removing
   fibre modes cannot touch it.

2. **One-plaquette dichotomy `[T]`+`[V]`.** In the exactly solvable `SU(2)` one-plaquette (`S2`/`G4`,
   longitudinal `v=e₃`, exact `G″=A−Bβ²`, `A=0.20511530`, `B=0.01105929`, negative for `β>4.30661`):
   - **`θ=0` (genuine fibre Savvidy).** The fibre `K` becomes indefinite at large `r`, but the coarse
     coupling is **exactly zero** (`c_r=0`, the `e₃` direction is decoupled from the transverse
     negative mode by axial symmetry), so `Σ=α>0` pointwise and `n₋(H)=n₋(K)`. **Excision works;
     the Schur complement is positive.**
   - **`θ=π` (the barrier, `G″<0` past `β=4.30661`).** The fibre is **positive definite everywhere**
     (`n₋(K)=0`), yet `H` is indefinite and `Σ<0` (at `β=6`: `Σ→−1.06`, `n₋(H)=n₋(Σ)=1`). The
     negative mode is a **collective `(radial, e₃)` mix**, not a fibre mode; and `e₃` is not
     orthogonal to it. **Excision of any fibre mode cannot restore positivity** — there is no fibre
     negative mode to excise, and the negativity is the `e₃`-direction force variance `D=Var(g/r)`
     itself (`∝β²`) overtaking the Haar stiffness `S=E[Φ′/r]=O(1)`.

3. **Strip (1×L) `[T]`+`[V]`.** The constant mode has `S_W≡0` **exactly** (uniform shift symmetry), so
   its Hessian block is Haar-only: the Gaussian flat-point Schur complement is `Σ₀=Φ″(0)=1/6`, and the
   exact constant-mode curvature (all non-constant modes integrated out by transfer matrix) is
   `G″(0)=0.2077→0.2165` (`L=2..8`), saturating at `c≈0.223`. **Positive, `L`-uniform, no Savvidy
   instability on the constant mode** (the negatives live in the gapped `4π²/L²` non-constant fibre).

4. **Rank-one defect lemma `[T]`.** Exact variance identity `G″ = S − D` with `S=E[α]` (bare
   stiffness) and `D=Var(f)` (force variance); `η=D/S`, and `1/G″ = 1/(S(1−η))` is its Gaussian-local
   value. Brascamp–Lieb on a log-concave fibre (`K≻0`) gives `G″(a₀) ≥ E[Σ|a₀] =: ρ_rest(a₀)`, hence
   ```
   χ = Var_μ(a) ≤ E[1/G″] ≤ 1/ρ_* ≤ (1+η)/ρ_*     whenever ρ_* := ess inf_{a₀} ρ_rest(a₀) > 0.
   ```
   The lemma is **true**, but the sharp form is `χ ≤ 1/ρ_*` — the `(1+η)` is cosmetic slack, and
   the **binding hypothesis is `ρ_* > 0`, i.e. `E⟨K⁻¹c,c⟩ < S` uniformly in `a₀`, which is NOT implied by
   the small defect `η=D/S`** (BL gives the wrong-direction `E⟨K⁻¹c,c⟩ ≥ D`; the *slack* `E⟨K⁻¹c,c⟩−D`, the
   non-Gaussianity, can dominate). §4 exhibits `G″>0, ρ_rest<0` at a *small* `η=0.26`.

5. **Acid test `[T]`.** Free massless field: `G″=ρ_rest=k²→0` as `k→0`, `χ→∞`; no positive floor is
   manufactured, and the positivity hypothesis `ρ_rest>0` fails exactly at `k=0` as required.

**One line.** *The inertia lemma `Σ>0 ⇔ n₋(H)=n₋(K)` is the exact surgery result: rank-one excision
restores positivity iff the Savvidy mode is fibre-confined (θ=0), and fails at the longitudinal
barrier (θ=π) where the unstable mode is collective. The rank-one defect lemma `χ≤(1+η)/ρ_rest` is
true but reduces to the sharp `χ≤1/ρ_*`, whose hypothesis `ρ_*:=ess inf_{a₀}ρ_rest>0` is a **uniform
fibre gap** `K⪰γI` — the exact missing 4D estimate, not the measured `D/S`.* **No QED. Clay 0%.**

---

## 1 · The inertia lemma `[T]`: what "eliminate the unstable mode" can and cannot do

Let `V` be the block potential, `v` a unit coarse direction, `a=⟨v,X⟩`; write the `v`-block Hessian
as
```
H(v) = [[α, c^ᵀ],[c, K]],   α = ∂_v²V,  c = Π_{v^⊥}∇(∂_vV),  K = Π_{v^⊥}(Hess V)Π_{v^⊥},
```
and the **Schur complement on the coarse direction**
```
Σ = α − c^ᵀ K⁻¹ c        (the Gaussian curvature on `v` after eliminating the fibre `K`).
```

> **Lemma Q3.1 (Schur signature; Haynsworth). [T]** If `K` is invertible,
> `det H = det K · Σ`, and the inertia is additive:
> `In(H) = In(K) ⊕ In(Σ)`, so `n₋(H) = n₋(K) + n₋(Σ)`. Consequently
> ```
> Σ > 0   ⇔   n₋(H) = n₋(K)      ("all negative modes are fibre-confined"),
> Σ indefinite   ⇔   n₋(H) > n₋(K)  ("part of the negativity is in the coarse direction").
> ```

*Proof.* Standard block determinant/inertia (Haynsworth 1968; Schur). Verified to `10⁻⁹` over random
blocks in witness (1). ∎

**What this says about the task.** "Integrate out the unstable mode exactly and compute the Schur
complement on the remaining modes" has a clean, dimension-free answer:
- If the negative mode(s) of `H` are **confined to the fibre** (`n₋(H)=n₋(K)`), then the fibre
  Schur complement `Σ` is **positive** on the coarse direction. Eliminating the fibre — exactly (with
  its double-well / nonlinear marginal) or in the Gaussian approximation — leaves a positive coarse
  curvature. **The surgery works.**
- If `H` has a negative mode **mixing the coarse direction with the fibre** (`n₋(H)>n₋(K)`), then `Σ`
  is indefinite. Eliminating *fibre* modes cannot remove it: the negative eigenvalue of `Σ` is the
  coarse-direction curvature itself. Removing the single negative eigenvector `φ₋` of `H` leaves a
  positive subspace, but `e₃ ∦ φ₋`, so this does **not** give a curvature "along `e₃`" — you would be
  computing the curvature of a *different* mode. **The surgery is the wrong operation for `H1`.**

The companion identity `G″ = S − D` (`S=E[α]`, `D=Var(f)`, the exact marginal curvature) completes the
picture: the Schur complement is the **Gaussian estimate** of `G″`, and Brascamp–Lieb makes it a lower
bound when the fibre is log-concave (§4).

---

## 2 · The one-plaquette model: the two regimes are exactly realized `[T]`+`[V]`

**Model (`S2`/`G4`, exact).** Flat chart `U=e^{iξ·σ/2}`, fixed boundary `B=e^{iθσ₃/2}`,
```
V(ξ)=h(ρ)+g(ρ)u,  u=ξ₃/ρ,  ρ=|ξ|,  h=β(1−cos(ρ/2)cos(θ/2))+Φ(ρ),  g=β sin(ρ/2)sin(θ/2),
Φ(ρ)=2log(ρ/2)−2log sin(ρ/2)   (Haar, ball ρ<2π),
```
coarse mode `v=e₃` (`a=ξ₃`), fibre = the transverse plane. At the slice `ξ₃=a₀=0`, in the basis
`(radial ξ̂, tangential ξ̂⊥, e₃)` the Hessian is
```
H = [[K_rad, 0, c_r],[0, K_tang, 0],[c_r, 0, α]],   α=h′(r)/r,  c_r=∂_r(g/r),  K_rad=h″(r),  K_tang=h′(r)/r.
```
The exact marginal curvature is G4's
```
G″(e₃;θ=π;a₀=0) = A − B β²,   A = −ρ_Haar″/ρ_Haar = 0.20511530,  B = 0.01105929,   (β* = 4.30661).
```
Witness (2) reproduces it to all printed digits by two routes (`G″ = S−D` analytic and a 5-point FD
of `log ρ(a)`), and the reference values are:

| `β` | `θ` | `S=E[α]` | `D=Var(f)` | `G″=S−D` | `E[Σ]=ρ_rest` | `η=D/S` |
|---|---|---|---|---|---|---|
| 2.5 | 0 | 0.72992 | 0.00000 | +0.72992 | +0.72992 | 0.000 |
| 2.5 | π | 0.20512 | 0.06912 | **+0.13599** | **+0.05549** | 0.337 |
| 6.0 | 0 | 1.59064 | 0.00000 | +1.59064 | +1.59064 | 0.000 |
| 6.0 | π | 0.20512 | 0.39813 | **−0.19302** | **−0.65673** | 1.941 |

**`θ=0`: the negative mode is fibre-confined.** Because `g∝sin(θ/2)=0`, the coarse–fibre coupling is
`c_r=0` **exactly**; the coarse direction is decoupled from the transverse negative mode. The fibre
`K` does go indefinite at large `r` (`K_rad=h″<0` for `r` past the Savvidy crossover; witness:
121/600 `r`-values at `β=6`), and `H` is indefinite precisely there, but `n₋(H)=n₋(K)` and
`Σ=α=Φ′/r>0`. **Excision (or simply noting the decoupling) makes `Σ` positive; the coarse curvature
is `G″=S>0`.** This is the clean `rank-one-fibre` regime `SM2` presupposes.

**`θ=π`: the negative mode is collective.** The fibre is positive definite **everywhere**
(`K_rad=Φ″>0`, `K_tang=Φ′/r>0`; witness: 0/600 `r`-values indefinite at `β=6`), yet `H` is indefinite
(411/600 `r`-values) and `Σ<0` (e.g. `Σ(1)=−0.17`, `Σ(2)=−0.81`, `Σ(3)=−1.06`), with `n₋(H)=n₋(Σ)=1`.
The single negative eigenvector is the collective `(radial,e₃)` mix, and `e₃` has O(1) overlap with
it. There is **no fibre negative mode to excise**; the negativity is carried by the `e₃`-direction
force variance `D=Var(g/r)=β²Var(sin(r)/2r)∝β²` overtaking the Haar stiffness `S=E[Φ′/r]=O(1)`. This
is exactly `SM2`'s **overlap** `⟨φ₋|c⟩=⟨H_cf φ₋|v⟩` at its **unsuppressed (worst-case) value**: `c_r`
*is* that overlap, and at `θ=π` it is maximal (no helicity/gauge suppression), which is why the rank-one
margin cannot rescue `G″` here. `SM2`'s ceiling `m_soft≤G″` presumes `φ₋` is a **fibre** eigenvector;
at the longitudinal barrier it is not, so the ceiling does not even bound the barrier negativity.

**Honest correction to the "double well" reading `[V]`.** At `θ=π` the marginal `ρ(a)` is **not**
symmetric: the Savvidy term `g(ρ)u=g(ρ)ξ₃/ρ` breaks `a→−a`. Witness (2) confirms `ρ(a)` is monotone
decreasing from `a=0` (no symmetric trough), and `G″<0` is a **log-convex shoulder** (`log ρ` curved
upward), not a symmetric well. The negativity is a **force-variance defect** `D>S` on a tilted
marginal, not a barrier of a symmetric double well. (The symmetric double-well normal form of `H6`/`S1`
is the *transverse* conditioned soft mode; it is not this longitudinal marginal.)

---

## 3 · The 1×L strip: positive, uniform, and instability-free `[T]`+`[V]`

On the periodic Cartan chain (`P2`/`P4` model: `S_W=−βΣcos(Δz/2)`, per-link exact Haar weight `ρ₁`),
the constant mode is the mean `a=L^{-1/2}Σz_i`. Its defining property is **`S_W≡0` on it** (uniform
shift leaves every plaquette invariant — the exact flat direction `[T]`, `P-WAVE-SYNTHESIS §2`). Hence:
- the constant-mode Hessian block of `S_W` is `0`, and of `Φ=−logJ` is `Φ″⪰1/6`; the coarse–nonconstant
  coupling from `S_W` is `0` (shift symmetry), so at the flat point `c=0` and the Gaussian Schur
  complement is `Σ₀=Φ″(0)=1/6`;
- the **exact** constant-mode Schur-complement curvature (all non-constant modes integrated out
  exactly, via `Z(k)=Tr[M(k)^L]`) is, at `β=2.5`,

| `L` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| `G″(0)` | 0.2077 | 0.2097 | 0.2115 | 0.2130 | 0.2144 | 0.2155 | 0.2165 |

saturating (`c≈0.223`, `q≈0.87`) — the `P2` periodic floor, reproduced here.

**Verdict for (ii).** The strip's Schur-complement curvature is **positive and `L`-uniform**; the
exact value (`≈0.21–0.22`) exceeds the Gaussian flat-point Schur complement (`1/6≈0.167`) because the
measure (non-Gaussian) increases the curvature relative to the quadratic estimate. There is **no
Savvidy instability on the constant mode**: `S_W` is exactly flat there, and the negative modes live in
the gapped (`4π²/L²`) non-constant fibre. The strip is therefore the *good* regime of the one-plaquette
dichotomy — `Σ>0` — in the thermodynamic-transfer-matrix setting.

---

## 4 · The rank-one defect lemma: true, sharp form, and the binding hypothesis `[T]`

Let `α=∂_v²V`, `f=∂_vV`, `K`, `c` as in §1, and let `μ(·|a₀)` be the conditioned fibre measure.

> **Lemma Q3.2 (exact variance identity + Schur/BL bound). [T]**
> (a) `G″(v;a₀) = E[α|a₀] − Var(f|a₀) =: S − D` (differentiation under the integral).
> (b) If `K ⪰ γ I` with `γ>0` on the fibre (log-concave conditioned measure), then Brascamp–Lieb gives
> `D = Var(f|a₀) ≤ E[⟨K⁻¹c,c⟩|a₀]`, hence `G″(v;a₀) = S − D ≥ S − E[⟨K⁻¹c,c⟩|a₀] = E[Σ|a₀] =: ρ_rest(a₀)`.
> (c) Let `ρ_* := ess inf_{a₀} ρ_rest(a₀)`. If `ρ_*>0`, then `V_eff″ = G″ ≥ ρ_*` a.e. in `a₀`, the coarse
> marginal is `ρ_*`-log-concave, and Brascamp–Lieb gives the susceptibility ceiling
> ```
> χ := Var_μ(a) ≤ E[1/G″] ≤ 1/ρ_* ≤ (1+η)/ρ_* ,      η = D/S .
> ```

*Proof.* (a) standard. (b) BL variance bound on the log-concave fibre with Hessian `K`; the Schur
monotonicity `Σ⪰essinf λ_min` is `V2` Lemma V2.1/Cor V2.2. (c) `G″≥ρ_*` a.e. ⇒ BL `Var(a) ≤ E[1/G″] ≤ 1/ρ_*`
(the BL/`1/Var` ceiling of `P-WAVE §2` — the exact object is `−ρ″/ρ`, and `χ=1/G″` is only its Gaussian
equality, here measured to `≈4%`), and `ρ_rest ≤ S−D=G″` because `E⟨K⁻¹c,c⟩ ≥ D`. ∎

*Note.* The essential infimum `ρ_*` (over `a₀`), not the single conditioned value, is what the ceiling
needs — this is exactly the `H6` typical-point-vs-essential-infimum object. At the barrier `G″(a₀=0)<0`,
so `ρ_*≤G″(0)<0` and the hypothesis fails outright.

**What is and is not true.**
- **True:** `χ ≤ E[1/G″] ≤ 1/ρ_* ≤ (1+η)/ρ_*` **whenever the fibre is log-concave and `ρ_*>0`.** This is
  a genuine theorem. Witness (4) verifies `G″≥ρ_rest` in every row of a Gaussian+quartic toy and on the
  one-plaquette numbers (`β=2.5,θ=π`: `ρ_rest(0)=+0.055>0 ⇒ 1/G″(0)=7.35 ≤ 1/ρ_rest(0)=18.0 ≤ (1+η)/ρ_rest(0)=24.1`;
  these are the `Gaussian-local` `1/G″` values, `χ` itself obeying `χ≤E[1/G″]`).
- **The `(1+η)` is cosmetic.** The sharp statement is `χ ≤ 1/ρ_rest`; the `(1+η)` factor weakens it and
  uses the *wrong* defect. The binding hypothesis is `ρ_rest>0`, equivalently `E⟨K⁻¹c,c⟩ < S` (Schur
  loss below bare stiffness).
- **Small `η=D/S` does NOT imply `ρ_rest>0`.** `E⟨K⁻¹c,c⟩` is the Schur loss; BL gives only the
  wrong-direction `E⟨K⁻¹c,c⟩ ≥ D`, so the **non-Gaussianity slack** `E⟨K⁻¹c,c⟩−D` can exceed `S−D` and
  drive `ρ_rest<0` while `G″>0`. Witness (4): `A=0.90, K=0.15, λ=8` gives `η=D/S=0.26` (**small**) yet
  `G″=+0.67>0` while `ρ_rest=−0.12<0` — the Schur bound fails though the curvature is positive. The
  slack is controlled by the **fibre gap** `K(y)⪰γ` (soft fibre ⇒ large slack), not by `D/S`.
- **At the barrier the hypothesis fails outright.** `β=6, θ=π`: `ρ_rest(a₀=0)=−0.657<0`, hence
  `ρ_* = ess inf_{a₀} ρ_rest ≤ −0.657 < 0`; the lemma does not apply and `G″(0)=−0.193<0`.

> **The exact missing estimate (4D) `[O]`.** A **uniform fibre gap**
> ```
> K ⪰ γ I ,  γ = γ_* > 0   uniformly in (L,a) and in the coarse mode,   ⇔   E⟨K⁻¹c,c⟩ ≤ (1−δ) S ,
> ```
> equivalently the lower bound `ρ_rest ≥ ρ_* > 0` in Lemma Q3.2(b). By Lemma Q3.1 this is exactly the
> statement that the negative modes of `H` are **fibre-confined** (`n₋(H)=n₋(K)`) at every coarse mode
> — i.e. the Savvidy mode never mixes into the `k→0` direction. This is `H1-RG-DESCENT`'s `(H-C)` and the
> Feshbach `Σ≥μ` (TCI); it fails precisely at the head where `γ_j→0` (`⟨K⁻¹c,c⟩=g²/γ_j→∞`). The measured
> defect `D/S≈2%` controls the force variance but **not** this gap — the two are related only through the
> BL slack in the wrong direction. Controlling the gap (or the slack) uniformly is the missing step.

---

## 5 · Witness — `research/v3/w7/q3_schur_defect.py` (`py -3.12 -I`, deterministic)

Quoted output (abridged):

```
(1) HAYNSWORTH:  det H = det K * Sigma  (all trials): True ;  n_-(H)=n_-(K)+n_-(Sigma): True
    K PD & Sigma<0 => one negative in H none in K (IRREDUCIBLE); K indef & n-(H)=n-(K) => Sigma>0 (EXCISABLE)
(2) ONE-PLAQUETTE (th=pi): G''=A-B b^2 reproduced (A=0.20511530, B=0.01105929)
    beta=2.5 th=pi: S=0.20512 D=0.06912 G''=+0.13599 E[Sigma]=+0.05549 eta=0.337
    beta=6.0 th=pi: S=0.20512 D=0.39813 G''=-0.19302 E[Sigma]=-0.65673 eta=1.941
    theta=0 : c_r=0 => Sigma=alpha>0, n-(H)=n-(K) (EXCISABLE)
    theta=pi: K>0 everywhere (0/600 indef) but H indef 411/600, Sigma<0 (COLLECTIVE, no fibre mode)
(3) STRIP beta=2.5: G''(0)=0.2077,0.2097,0.2115,0.2130,0.2144,0.2155,0.2165 (L=2..8); fit c=0.223
    Gaussian flat-point Schur Sigma0=Phi''(0)=1/6 ; exact > 1/6 ; positive & L-uniform
(4) DEFECT LEMMA: G''=S-D; Gaussian-local 1/G''=1/(S(1-eta)); G''(a0)>=E[Sigma|a0]=rho_rest
    rigorous ceiling: Var(a)<=E[1/G'']<=1/rho_*<=(1+eta)/rho_*  (rho_*=essinf_a0 rho_rest)
    one-plaquette beta=2.5: rho_rest(0)=+0.055>0 => 1/G''(0)=7.35 <= 1/rho_rest(0)=18.0 <= 24.1
    toy A=0.90,K=0.15,lam=8: eta=0.26 SMALL yet G''=+0.67, rho_rest=-0.12<0 (Schur bound fails)
(5) ACID: free massless G''=rho_rest=k^2 -> 0 (k=1,0.3,0.1,0.02,0.005 -> 1,.09,.01,4e-4,2.5e-5); no floor
```

Evidence quality: (1) exact linear algebra; (2) exact Gauss–Legendre marginal quadrature (two routes
agree; `G″=A−Bβ²` to all digits) + exact per-`r` inertia; (3) exact transfer-matrix `Tr[M^L]` and
5-point `log ρ`; (4) exact `S,D,E[Σ]` on the one-plaquette and deterministic quadrature in the toy;
(5) exact. It is a **witness, not a proof**; the theorem content is Lemmas Q3.1–Q3.2 and the exact
`G″=A−Bβ²`; the open step is the uniform fibre gap of §4.

---

## 6 · Ledger, acid tests, reading

| item | statement | label |
|---|---|---|
| inertia lemma | `det H=det K·Σ`, `In(H)=In(K)⊕In(Σ)`, `Σ>0⇔n₋(H)=n₋(K)` | **[T]** |
| excision criterion | rank-one fibre excision restores positivity **iff** unstable mode is fibre-confined | **[T]** |
| one-plaquette exact `G″=A−Bβ²` (`θ=π`,`v=e₃`) | `A=0.20511530`, `B=0.01105929`, `β*=4.30661` | **[T]** (S2/G4) |
| `θ=0`: `c_r=0`, `Σ=α>0`, `G″=S>0` | fibre-confined; excision works | **[T]**+**[V]** |
| `θ=π`, `β=6`: `K≻0`, `Σ<0`, `G″<0`, collective mode | excision cannot help; no fibre negative mode | **[T]**+**[V]** |
| `θ=π` marginal is tilted (not a symmetric double well) | `g u` breaks `a→−a`; `G″<0` = log-convex shoulder | **[V]** (correction) |
| strip constant mode: `S_W≡0`, `Σ₀=1/6`, `G″=0.2077..0.2165`, saturating | positive, `L`-uniform, no instability | **[T]**+**[V]** |
| exact variance identity `G″=S−D`; Gaussian-local `1/G″=1/(S(1−η))` | differentiation under the integral | **[T]** |
| Schur/BL bound `G″(a₀)≥E[Σ|a₀]=ρ_rest(a₀)` (fibre `K≻0`) | `D≤E⟨K⁻¹c,c⟩`; `ρ_rest≤G″` | **[T]** (BL, V2) |
| defect lemma `χ≤E[1/G″]≤1/ρ_*≤(1+η)/ρ_*` **if `ρ_*>0`** | sharp form `χ≤1/ρ_*`; `(1+η)` cosmetic | **[T]** |
| small `η=D/S` ⇒ `ρ_*>0` | **false**: BL slack `E⟨K⁻¹c,c⟩−D` can dominate (`η=0.26, ρ_rest<0, G″>0`) | **[NG]** (of that implication) |
| barrier `β=6`: `ρ_rest(0)<0` | lemma inapplicable; `G″(0)<0` | **[V]** |
| **uniform fibre gap `K⪰γI` (`ρ_rest≥ρ_*>0`), ⇔ `n₋(H)=n₋(K)` all modes** | the missing 4D estimate (= `(H-C)`/TCI = `H1`) | **[O]** |
| free-field `G″=ρ_rest=k²→0` | acid test PASS | **[T]**+**[V]** |
| **Clay** | unbuilt | **0%** |

**Acid tests.**
- **Free massless `[T]`+`[V]`.** `V=½k²x²`: `α=k²`, `c=0`, `K=k²`, `Σ=k²`, `D=0`, `G″=S−D=k²`. All
  `→0` as `k→0`; `χ→+∞`. The tool manufactures **no** positive floor and the hypothesis `ρ_rest>0`
  correctly fails at `k=0`. (Contrast: the one-plaquette `β=0` Haar floor `G″=0.2051` is a compact-
  measure curvature, not a free-field floor — the "`1/Var` vs `−ρ″/ρ`" distinction of `P-WAVE §2`.)
- **Discriminating `[T]`.** In the one-plaquette `θ=0` the free-`D` (Gaussian, `D=0`) limit gives
  `G″=S=α>0` from the *Haar* measure — a genuine compact-measure curvature, not a manufactured gap;
  it does not certify an atom (`H1_does_not_imply_atom`, `P6`).

**Reading.** The Schur complement answers the task cleanly and negatively for the hard direction:
*rank-one elimination restores positivity exactly when the unstable mode is fibre-confined, and the
longitudinal barrier is not that case.* The rank-one defect lemma is true, but its content is the sharp
`χ≤1/ρ_*`, whose hypothesis `ρ_*>0` is a **uniform fibre gap** — the same object as `(H-C)`, the
Feshbach `Σ≥μ`, and `H1`. The measured defect `D/S≈2%` is the wrong proxy: the Schur loss is bounded
by the fibre gap, related to the force variance only through the non-Gaussianity slack, and at the
barrier it is the `e₃`-direction force variance that is large. What is missing for 4D is exactly the
uniform lower bound `K⪰γI` on the coarse fibre — equivalently `ess inf G″≥λ` on the typical bulk with
the collective Savvidy mode excluded — which is `H1` itself, not a corollary of the defect data.

*No QED. Clay 0%.*

---

## 7 · References (internal, this session)

- `SM2-RANK-ONE.md` §1–§3 (overlap lemma `⟨φ₋|c⟩=⟨H_cfφ₋|v⟩`, margin identity `G″=(S−D_Sav)+|D₊|`,
  ceiling `m_soft≤G″`, refined `D_Sav≤g²Var`).
- `V2-TRUE-BOUND.md` §2 (Lemma V2.1 marginal strong convexity; Cor V2.2 exact Schur, no loss;
  `G″=E[Σ]−D`).
- `H1-RG-DESCENT.md` §2, §4 (Lemma H1.1 `G″_j⪰(1+η_j)⁻¹G″_{j+1}`, Schur loss `⟨K⁻¹c,c⟩=g²/γ_j→∞`,
  `(H-C) K⪰γ_jI`); `S3-UV-FESHBACH.md` (`G″=E[Σ]−D`).
- `H6-HEAD-SPECTRAL.md` §3 (typical point vs essential infimum); `G4-GENERAL-BOUND.md` §2–§3
  (`G″=A−Bβ²`, `β*=4.09` vs one-plaquette `β*=4.30661`); `P-WAVE-SYNTHESIS.md` §2–§3 (periodic floor,
  `−ρ″/ρ` object); `docs/paper-v3/main.tex` §6.5–6.7.
- Standard: Haynsworth (1968) inertia additivity / Schur complement; Brascamp–Lieb (1976) variance
  bound; Prékopa (1971); Savvidy (1977) / Nielsen–Olesen (1978) negative mode. All cited from the
  corpus; no third-party sources downloaded this session.
- Witness `research/v3/w7/q3_schur_defect.py` (deterministic; output quoted §5).

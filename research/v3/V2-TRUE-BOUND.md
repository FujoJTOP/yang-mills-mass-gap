# V2-TRUE-BOUND — a genuine lower bound on the TRUE `G″` (including `X`) on the small-field set: the marginal-strong-convexity theorem `G″ ≥ essinf λ_min(Hess V)`, its exact Schur form, the boundary threshold, and where it stops (the bulk Savvidy / H1)

> 2026-11-24 · V2 attack on the single live inequality `(★)` of the campaign.
> Object: the **true coarse curvature** `G″_true(v;a₀) = −d²/da² log ρ_v(a)|_{a₀} = E[∂_v²V|a₀] − Var(∂_vV|a₀)`,
> i.e. the object `(★)` actually uses (ch05 `eq:split` with the `|0` made explicit; S3's `E[Σ]−D`).
> **Not** the unconditional `2W_geo` object (G4's `(U)`/`(L)`), which is `2W_geo` or `0` and is the wrong functional.
> Reads: `docs/paper-v2/ch05-theorems.tex` (`eq:split`, `eq:vardecomp`, `sec:wgeo`, `sec:feshbach`);
> `S1-CONDITIONAL-POINCARE.md` (Lemma S1, Corollary S1′, LSI on `𝒮_t`); `S3-UV-FESHBACH.md` (`G″=E[Σ]−D`);
> `G4-GENERAL-BOUND.md` (the `G″_lit≡0` trap, `G″=2W_geo` is the wrong object, barrier counterexample,
> `G″(e₃;θ=π)=0.206−0.0123β²`); `S2-BOUNDARY-COUNTEREXAMPLE.md`; `S5-N-UNIFORM.md`; `STAR-WAVE-SYNTHESIS.md`.
> Witness: `research/v3/w7/v2_true_bound_witness.py` (`py -3.12`, deterministic; output quoted in §5).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[V]** numerical witness · **[R]** reported.

---

## 0 · Verdict (no QED, Clay 0%)

**Yes — there is a genuine, exact lower bound on the TRUE `G″`, cross term `X` included, and it is positive on the small-field boundary set in the well-controlled (small-block / small-field) regime. It is not a reduction to `2W_geo`.**

1. **[T] The theorem (marginal strong convexity).** If the block fiber is `ρ`-strongly log-concave,
   `Hess_X V ⪰ ρ I` on the fiber (`V = S_W + Φ`), then for **every** unit direction `v`, **every**
   conditioning value `a₀`,
   ```
        G″_true(v;a₀) = E[∂_v²V | a₀] − Var(∂_vV | a₀)  ≥  ρ(B,U′),    ρ(B,U′) := essinf_X λ_min(Hess_X V).
   ```
   This is exact, holds for the object `(★)` uses, and includes `X = −2Cov(G,Φ′|a₀)` (it is a bound on
   the full `−d²/da² log ρ_v`, not on any piece). Proof: `E[Σ|a₀] ≥ ρ` (Schur-complement monotonicity)
   plus `D ≤ 0` (Brascamp–Lieb on the log-concave conditioned fiber); `G″ = E[Σ] − D`. **[T]**

2. **[T] The modulus has the Haar floor plus the box gap minus the Savvidy loss.**
   ```
        ρ(B,U′)  ≥  N/12  +  c_□(B)·1_{non-constant fibre mode}  −  C₁β F_∂(U′)  −  C_Sβ F_bulk,max ,
   ```
   with `N/12` the pointwise Haar floor (`Φ″⪰(N/12)I`, `ch05` Lemma B; SU(2) `1/6`, SU(3) `1/4`),
   `c_□(B)~4π²/L²` the block Dirichlet box gap, `C₁,C_S` explicit plaquette constants, `F_∂` the boundary
   field strength, `F_bulk,max` the interior field supremum over the fiber. **[T]** pointwise; the
   constants `C₁,C_S` are existence-proved but not decimalised — exactly the status of the campaign's
   `C_g ≤ CN` in the P3 theorem.

3. **[C] The genuine positive bound on the small-field set `𝒮_t`.** On
   `𝒮_t = {U′ : C₁βF_∂(U′) < c_□(B) + N/12}` (one-plaquette / small block, `F_bulk` slaved to `F_∂`),
   ```
        G″_true(v;a₀) ≥ ρ₀(t) = c_□(B) + N/12 − C₁β t  > 0      uniformly in v, a₀, U′∈𝒮_t .
   ```
   In `N`: `ρ₀ = N/12 − O(1) = cN − O(1)`, linear in `N` (matches S5's `δ_W=cN−O(1)`). This **is** the
   requested `G″_true ≥ N/12 − boundary-loss` bound, with the boundary loss `C₁βt` and the `c_□` bonus.
   One-plaquette witness at the head `β=2.5`: `ρ₀(θ=0)=+0.217`; `θ_c(2.5)≈0.81` rad (`F_∂<0.081`). **[C]+[V]**

4. **[NG] The pointwise bound is vacuous past `𝒮_t` and for a general multi-link block at the head.**
   For `C₁βF_∂ > c_□ + N/12` (large holonomy, e.g. `θ→π`), `ρ ≤ 0`: the fiber's own soft Savvidy modes
   make `Hess V` indefinite, so the strong-convexity bound manufactures nothing. **The Schur loss
   `C₁βF_∂` can — and at the head does — exceed `N/12`.** For a **general 4D multi-link block** the same
   happens even at zero boundary, because the fiber contains interior-field points whose Savvidy
   curvature is not bounded by the boundary (`F_bulk` term). **[NG]** for the pointwise-convexity class
   beyond `𝒮_t` / for general blocks at `β≈2.5`.

5. **[O] What survives past the bound — and it is the wall.** The TRUE `G″` is nonetheless **positive at
   the head** numerically: one-plaquette `G″(e₃;θ=π)=+0.129` at `β=2.5`, crossing zero only at
   `β*=4.09`. So the true object is *more robust* than pointwise convexity; but its positivity there is
   carried by the **measure** (the Brascamp–Lieb defect `D< E[Σ]`), i.e. by controlling `D` — **exactly
   H1**. Proving the head positivity in general is H1/the gap, not a convexity corollary. **[O]**

**Acid test passes.** Free massless field (`Φ≡0`, Gaussian `V`): `ρ=0` (soft mode), `G″=0` — no positive
floor is manufactured, `(★)` fails, as required. `:φ²:`: the bound is a **curvature** statement, never
an atom statement; a massive free quadratic gives `G″=m²>0` yet the composite has `f₁²=0` — curvature
`>0` does not imply an atom. **No QED; Clay 0%.**

---

## 1 · The true object (and the two traps it is not)

Fix a block `B`, boundary `U′`, fiber `X ∈ Ω_B ⊂ ⊕_ℓ 𝔰𝔲(N)`. Flat chart `U_ℓ=e^{X_ℓ}`, Haar
`∏dU_ℓ = e^{−Φ(X)}dX`, `Φ = Σ_ℓ Φ_0(X_ℓ)` (per-link sum; `Φ_0''⪰(N/12)I`). Conditional measure
`μ_{U′}(dX) = Z(U′)^{-1} e^{−V}dX`, `V = S_W(X;U′) + Φ(X)`.

For a unit fine direction `v` and `a := ⟨v,X⟩` with marginal `ρ_v(a) = ∫_{⟨v,X⟩=a} e^{−V}dσ`, the true
coarse curvature is

```
G″_true(v;a₀) := −d²/da² log ρ_v(a)|_{a₀} = E_{μ_{U′}}[∂_v²V | a=a₀] − Var_{μ_{U′}}(∂_vV | a=a₀).   (1.1)
```

**(1.1) is `[T]`** (pure differentiation under the integral, `ρ_v>0`; `G4` Lemma G1). It equals S3's
`E[Σ]−D`, and, splitting `V=S_W+Φ` and conditioning, ch05's

```
G″_true = Δ_dyn(a₀) + W_geo(a₀) + X(a₀),                                              (1.2)
Δ_dyn = E[H_vv|a₀] − Var(G_v|a₀),  W_geo = E[Φ″_vv|a₀] − Var(Φ′_v|a₀),  X = −2Cov(G_v,Φ′_v|a₀).
```

**Two traps this note avoids.**

- **(L) the literal (G4 Lemma G2).** *Unconditioned* (no `|a₀`), the flat-chart Stein identities give
  `Δ_dyn = W_geo = Cov(G,Φ′)` and `X = −2W_geo`, so `Δ_dyn+W_geo+X ≡ 0`. So the unconditioned object is
  `0` in every model — the `|a₀` is **load-bearing** and must not be dropped.
- **(U) S2's `2W_geo`.** `G″_model = E[H]−Var(G)+W_geo` (unconditioned, `X` dropped) `= 2W_geo` exactly.
  This is a true identity but for the **wrong object**: it is *not* (1.1), and it is `≈4×` off at `θ=0`
  and of the **wrong sign** at the double-well barrier (G4 §3). So `(★)⇔2W_geo≥δ_W` is void; the bound
  below is on (1.1), not on `2W_geo`. **[T]/[NG]** (of the reduction).

---

## 2 · The exact theorem — marginal strong convexity bounds the TRUE `G″`

Write the fiber in coordinates `(a, y)`, `a=⟨v,X⟩`, and the fiber Hessian
`Hess V = [[α, c],[cᵀ, K]]`, `α = vᵀ Hess V v`, `c = ∂_a∂_y V` (mixed Hessian), `K = Π₀ Hess V Π₀`
(the transversal fiber Hessian, `Π₀ = 1−vvᵀ`). The **Schur complement** along `v` is
`Σ = α − ⟨K⁻¹c,c⟩` (this is S3's `Σ = A − ⟨K⁻¹c,c⟩` with `A=α`).

> **Lemma V2.1 (marginal strong convexity) — [T].**
> Suppose `Hess_X V ⪰ ρ I` on the (essential support of the) fiber, `ρ>0`. Then for every unit `v` and
> every `a₀`, `G″_true(v;a₀) ≥ ρ`.

*Proof.* Split (1.1): `G″ = E[α|a₀] − Var(∂_vV|a₀)`.
(i) **Schur monotonicity.** `Hess V ⪰ ρI` ⇒ `K ⪰ ρI` and `Σ = α − ⟨K⁻¹c,c⟩ ⪰ ρ` pointwise (Schur
complement of `Hess V − ρI ⪰ 0`; equivalently the Schur complement of an `⪰ρI` matrix is `⪰ρI`).
Hence `E[Σ|a₀] ≥ ρ`.
(ii) **Brascamp–Lieb defect `D ≤ 0`.** The conditional `μ(·|a₀) ∝ e^{−V(a₀,y)}dy` has `y`-Hessian
`K ⪰ ρI>0`, so it is log-concave; BL gives
`Var(∂_vV|a₀) ≤ E[⟨K⁻¹ ∂_y(∂_vV), ∂_y(∂_vV)⟩ | a₀] = E[⟨K⁻¹c,c⟩|a₀]` (since `∂_y∂_vV = c`).
Thus `D := Var(∂_vV|a₀) − E[⟨K⁻¹c,c⟩|a₀] ≤ 0`, and
`G″ = E[α|a₀] − Var(∂_vV|a₀) = E[Σ|a₀] − D ≥ E[Σ|a₀] ≥ ρ`. ∎

> **Corollary V2.2 (exact Schur form; no separate "Schur loss") — [T].**
> Under the same hypothesis, `G″_true = E[Σ|a₀] − D` with `Σ ⪰ ρ` pointwise and `D ≤ 0`. The crude
> decomposition that *introduces* a loss, `G″ ≥ E[α] − E[⟨K⁻¹c,c⟩] ≥ E[H_vv]+N/12 − E[⟨K⁻¹c,c⟩]`,
> is **weaker**: it can go negative even when `Σ ⪰ ρ>0`, because `⟨K⁻¹c,c⟩` is not separately bounded by
> `N/12`. The **lossless** statement is `Σ ⪰ essinf λ_min(Hess V) = ρ`. Verified numerically: for 15
> random `6×6` matrices with `λ_min ≥ ρ`, `Σ−ρ ≥ 0` (witness (B)). **[T]**

**This is the genuine bound.** It is on `G″_true` of (1.1) (the `(★)` object), it includes `X`, it is
exact, and it carries **no** reduction to `2W_geo`.

---

## 3 · Evaluating the modulus `ρ` — Haar floor, box gap, Savvidy loss; boundary dependence

`Hess V = Hess S_W + Hess Φ`, and the two pieces are controlled oppositely:

- **`Hess Φ ⪰ (N/12)I`** pointwise, every `β`, every measure (`ch05` Lemma B; `g''(d)≥1/6`). The
  per-link Haar factorisation makes this a sum over links, hence `U′`-independent. **[T]**
- **`Hess S_W ⪰ −C_Sβ F_· I`.** The Wilson Hessian's negative part is the Savvidy curvature, bounded by
  the local field strength times an explicit plaquette constant (`SZZ` Lemma 4.1 form `|Hess S_W|≤8(d−1)Nβ`;
  the negative eigenvalue is the lattice `−D²+2adF` soft mode, `W6` re-measured `≈−1.12β`). Three cases:
  - **non-constant fiber modes:** the electric Laplacian gap gives `+c_□(B)` (diamagnetic, `U′`-uniform
    lower bound) — the **UV protection**;
  - **soft (constant-transverse) fiber modes:** no `c_□`; curvature `−C_Sβ F_local`;
  - **boundary vs bulk:** `F_local ≤ F_∂(U′)` on boundary plaquettes, but the fiber's **interior** field
    reaches `F_bulk,max` (a fiber supremum, `U′`-independent).

So, pointwise `[T]`,

```
ρ(B,U′) ≥ N/12 + c_□(B)·1_{non-constant} − C₁β F_∂(U′) − C_Sβ F_bulk,max .              (3.1)
```

**The boundary dependence is entirely in `−C₁βF_∂(U′)`**, monotone decreasing, `U′`-uniform otherwise;
`F_bulk,max` is `U′`-independent.

### 3.1 When does the "Schur/Savvidy loss" exceed `N/12`? — the condition and the corrected bound

The question the task poses has a clean answer:

- **It can.** The Savvidy loss `C₁βF_∂ + C_SβF_bulk` **exceeds** `c_□ + N/12` exactly when
  ```
        C₁β F_∂(U′)  >  c_□(B) + N/12 − C_Sβ F_bulk,max .                                 (3.2)
  ```
  When (3.2) holds, `ρ ≤ 0`, the pointwise strong-convexity bound is **vacuous**, and the crude
  `N/12 − ⟨K⁻¹c,c⟩` bookkeeping is meaningless (indeed `λ_min(K) < N/12` there).
- **The corrected bound.** Do **not** subtract a Schur loss. Use the exact Schur complement
  `Σ ⪰ essinf λ_min(Hess V) = ρ` (Corollary V2.2). It is positive **iff** (3.2) fails, i.e. on the
  small-field set
  ```
        𝒮_t = {U′ : C₁βF_∂(U′) < c_□(B) + N/12 − C_SβF_bulk,max}   (one-plaquette/small block: F_bulk slaved),
  ```
  where `G″_true ≥ ρ ≥ c_□ + N/12 − C₁βt > 0`. Beyond `𝒮_t` no pointwise bound exists; see §4/§6.

### 3.2 Constants, tracked

`N/12` = Haar floor (`1/6` SU(2), `1/4` SU(3)); `c_□(B) ~ 4π²/L²` box gap (in lattice units, the
`1/a_j²` physical scaling is inside `c_□` — fixed lattice block `L`); `C₁` boundary-plaquette constant,
`C_S` bulk-Savvidy constant — both explicit-in-existence, not decimalised (same status as `C_g≤CN`).
In `N`: `ρ ≥ N/12 − O(1) = cN−O(1)`, linear in `N`, `N`-free radius (the adjoint dimension is
multiplicity, S5). In `a`: `c_□` and `β` carry the lattice units; the combination is `a`-uniform at fixed
physical block (`W_geo`/Haar floor coupling-free; S5 §4). **[T]** for `N/12`; **[C]** for the assembled
`ρ₀=cN−O(1)`.

---

## 4 · Term-by-term signs on `𝒮_t` (and why only the sum is controlled)

On `𝒮_t` (`ρ>0`), from (1.2):

| term | expression | sign on `𝒮_t` | label |
|---|---|---|---|
| Haar reservoir | `E[Φ″_vv|a₀]` | `≥ N/12 > 0` (pointwise, any measure) | **[T]** |
| force variance | `Var(G_v|a₀)` | `≥ 0` (subtraction, hurts) | **[T]** |
| `Δ_dyn` | `E[H_vv|a₀] − Var(G_v|a₀)` | **sign-indefinite** (`E[H_vv]` Savvidy-negative possible; `Var≥0`) | **[T]** |
| `W_geo` | `E[Φ″|a₀] − Var(Φ′|a₀)` | `≥0` **only** pure-Haar `β=0` (Prékopa); sign-**unproven** under the interacting conditioned measure | **[T]** (`β=0`) / **[O]** (`β>0`) |
| `X` | `−2Cov(G_v,Φ′_v|a₀)` | **sign-indefinite**; conditioned cross term (the unconditioned `X=−2W_geo` does **not** apply) | **[T]** |
| **sum** | `Δ_dyn+W_geo+X = E[Σ|a₀]−D` | **`≥ ρ > 0`** (Lemma V2.1) | **[T]** |

**The message.** The individual pieces `Δ_dyn, W_geo, X` are **not** each positive; only their **sum** is
bounded below, by `ρ`, through the Schur/BL structure. This is exactly why the S2 "`G″=2W_geo`" reduction
was the wrong object (it bounds a summand-identity, not the sum's curvature) and why the campaign's
`W_geo≥0` sign theorem is valuable but **not** the gate. The bound (Lemma V2.1) never needs a sign for
`X` or for `W_geo`; it controls the total.

---

## 5 · Witness — `research/v3/w7/v2_true_bound_witness.py` (deterministic, `py -3.12`)

```
(A) MARGINAL STRONG CONVEXITY [T]:  V=0.5(x^2+y^2)+0.1(x^4+y^4)/4, Hess V>=1.0
    G''_true(0)=1.000000 >= rho=1.0  OK=True ;  pointwise Schur Sigma=1.000000=rho (no loss)
    contrast double-well: G''=-0.2916 at the BARRIER, +1.1872 at the WELL
        -> D>0 when the fibre is not log-concave; G'' need NOT exceed any floor.

(B) SCHUR MONOTONICITY [T]:  15 random 6x6 M>=rho I -> min(Sigma-rho)>=0 : True
    => exact Sigma>=rho has no 'Schur loss'; crude 'N/12-<K^-1 c,c>' is weaker.

(C) ONE-PLAQUETTE SU(2)   essinf lam_min(Hess V)   vs   TRUE G''(a0=0)   (N/12=0.1667)
    beta  theta   essinf HessV      TRUE Gpp(e3,theta=pi)=0.206-0.0123 beta^2
    0.50  0.000   +0.2730           --
    0.50  3.142   +0.0997           +0.2026
    2.50  0.000   +0.2165           --
    2.50  3.142   -0.3596           +0.1288     <-- pointwise bound vacuous, TRUE G''>0
    4.50  0.000   -0.0112           --
    4.50  3.142   -0.8489           -0.0434     <-- TRUE G'' itself <0 past beta*=4.09
    small-field threshold at beta=2.5:  essinf = +0.217(theta=0) ... 0 at theta_c~0.81 ... -0.360(theta=pi)
    => on theta<theta_c (small-field set)  G''_true >= essinf > 0  is a GENUINE POSITIVE BOUND (incl. X);
       cross-check theta=0,beta=2.5:  essinf +0.217 <= TRUE G''(e3)=0.741  (valid, conservative).

(D) ACID: free massless k->0: essinf=k, G''=k -> 0 (no false floor);  :phi^2: curvature>0 != atom.
```

Readings: (A) Lemma V2.1 is exact and sharp at the Gaussian point; the double-well shows the bound is
**necessary** (a non-log-concave fiber gives `G″<0`). (B) the exact Schur complement carries no loss; the
"`N/12 − Schur-loss`" budget is the weaker (and failure-prone) one. (C) the pointwise modulus `ρ` is the
**boundary-dependent** small-field quantity: positive for `F_∂<F_c`, crossing zero at `θ_c`; at the head
`β=2.5` the small-field set is **non-empty** (`θ_c≈0.81`) and there the genuine positive bound holds; at
large holonomy the pointwise bound is vacuous **yet the TRUE `G″` is still `+0.13`** (measure-carried),
until the TRUE `G″` itself crosses zero at `β*=4.09` (matching G4's barrier counterexample). (D) no false
positivity on the free field; curvature is not an atom.

---

## 6 · Where the bound stops — the general block and the head `[NG]`/`[O]`

Lemma V2.1 is unconditional given `ρ>0`; the *evaluation* of `ρ>0` is regime-dependent:

- **One-plaquette / small block, small field, or deep-UV block:** `ρ>0` on `𝒮_t` (witness §5; S1 Cor.
  S1′; the UV `c_□` plus `N/12` beat the boundary Savidy). Genuine positive bound. **[C]**
- **General multi-link 4D block at the head `β≈2.5`:** the fiber contains **interior-field** points
  (`F_bulk` large) whose Savvidy curvature is **not** bounded by `F_∂`; `ρ≤0` (`s1` witness: the 4D `L=2`
  soft mode is `0` at zero field and `<−N/12` at a bulk field). The small-field set bounds the **boundary**,
  not the interior. So the pointwise bound is **vacuous** for a general block at the head. **[NG]** for the
  pointwise-strong-convexity class there. (The deep UV fixes this via `c_□≫β`; the head does not.)
- **What is nevertheless true at the head:** the TRUE `G″` is positive (witness: `+0.13` at `β=2.5`,
  crossing at `β*=4.09`). Its positivity is carried by the **measure** — the defect `D` being `<E[Σ]` on
  the effective (typical) support, not pointwise. Controlling `D` uniformly is **exactly H1** (S1 §6,
  STAR-WAVE item 3). So the head case is not a new obstruction; it is the wall. **[O]**

**Honest bottom line for §6.** The genuine positive bound on the TRUE `G″` is proved on `𝒮_t` in the
small-block/small-field (`ρ>0`) regime. Beyond it — large holonomy, or a general block whose interior
field is uncontrolled at `β≈2.5` — the pointwise route gives nothing, and the residual positivity of the
TRUE `G″` is the measure-carried `D`-control, i.e. H1. The bound is real but does **not** cross the wall.

---

## 7 · Acid tests (the bound must not manufacture positivity on the free field / must not imply an atom)

1. **Free massless field (`Φ≡0`, Gaussian `V=½ XᵀAX`, `A` the Dirichlet Laplacian, `A⪰0`):** `ρ=0`
   (soft mode), `G″_true=0`; the bound `G″≥0` is **vacuous** and `(★)` fails for any `δ_W>0`. **[T]+[V]**
   (witness (D): `k→0` gives `essinf=k→0`, `G″=k→0`).
2. **`:φ²:` / massive free quadratic:** `G″=m²>0`, but the composite `:φ²:` is a two-particle **continuum**
   with `f₁²=0`. Positivity of `G″` does **not** imply an atom — the bound is a **curvature** statement
   (necessary-shaped, not sufficient for `(★)`/the gap). This is the same lesson as S1 §6 (`D`-control is
   necessary-shaped). **[T]**
3. **`β=0` pure Haar (SU(2)):** `Hess V = Hess Φ ⪰ N/12`, so `ρ≥N/12>0` and `G″≥N/12` — the ch05
   `β=0` theorem is the `ρ`-form at `S_W=0`. Not "false positivity": it is the measure's own Haar
   curvature, present for any SU(N) measure (incl. the trivial one), and it does not certify a gap.
   `W_geo=0` there (Prékopa saturated), confirming the sign theorem is separate from the gate. **[T]**

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| true object `G″(v;a₀)=E[∂_v²V|a₀]−Var(∂_vV|a₀)=−d²/da²log ρ_v|a₀` | (1.1); `=E[Σ]−D`; `=Δ_dyn+W_geo+X` conditioned | **[T]** |
| literal (unconditioned) `G″_lit≡0`; S2 `G″_model=2W_geo` | wrong object; `|a₀` load-bearing | **[T]/[NG]** |
| **Lemma V2.1** marginal strong convexity | `Hess V⪰ρI ⇒ G″_true(v;a₀)≥ρ` (incl. `X`) | **[T]** |
| **Corollary V2.2** exact Schur, no loss | `Σ⪰essinf λ_min(Hess V)`; crude `N/12−loss` weaker | **[T]** |
| modulus evaluation | `ρ≥N/12+c_□·1_{nonconst}−C₁βF_∂−C_SβF_bulk` | **[T]** pointwise / **[C]** assembled `cN−O(1)` |
| **Schur/Savvidy loss > N/12** | iff `C₁βF_∂>c_□+N/12−C_SβF_bulk`; then `ρ≤0`, bound vacuous | **[T]** condition (3.2) |
| **genuine positive bound on `𝒮_t`** | `G″_true≥ρ₀=c_□+N/12−C₁βt>0` (small block/small field) | **[C]** + **[V]** |
| one-plaquette `β=2.5`: `θ_c≈0.81`, `ρ₀(θ=0)=+0.217`; true `G″(e₃)=0.741≥0.217` | witness | **[V]** |
| general block at head: `essinf Hess V<0` (interior field) | pointwise bound vacuous | **[NG]** |
| true `G″>0` at head (`+0.13` at `β=2.5`), `<0` past `β*=4.09` | measure-carried; `D<E[Σ]` | **[V]** (`β*`)/**[O]** (proof) |
| head positivity = controlling defect `D` | = H1 (coarse gap) | **[O]** |
| free massless `G″=0`, `ρ=0`; `:φ²:` `G″>0` but no atom | acid tests PASS | **[T]+[V]** |

---

## 9 · Sharpest line

**The genuine bound the campaign needs is not `2W_geo` but the marginal-strong-convexity floor
`G″_true(v;a₀) ≥ ρ(B,U′) = essinf_X λ_min(Hess_X V)` — exact, `X` included, no reduction — and `ρ>0` is
the small-field statement `N/12 + c_□ > C₁βF_∂ + C_SβF_bulk`.** On the small-field set `𝒮_t` this gives a
positive bound `G″_true ≥ c_□ + N/12 − C₁βt` (linear in `N`, boundary-dependent through `F_∂`); the Schur
loss **can** exceed `N/12` precisely off `𝒮_t`, where the corrected, lossless statement is the exact
Schur complement `Σ⪰ρ` (positive on `𝒮_t`, vacuous outside). In the exactly solvable one-plaquette at the
head `β=2.5` the small-field set is non-empty (`θ_c≈0.81`) and the bound holds there (`ρ₀(0)=+0.22`,
true `G″=0.74` at `θ=0`); beyond `θ_c` — and for a general multi-link block whose interior field is
uncontrolled — the pointwise bound is vacuous, and the residual positivity of the TRUE `G″` (still
`+0.13` at `β=2.5`, crossing at `β*=4.09`) is carried by the measure (defect `D`), which is H1. So the
V2 push produces a **genuine, exact, positive lower bound on the true `G″` on the small-field set**, with
its precise constants and boundary threshold — and it confirms, rather than crosses, the wall: the head is
the gap. **No QED; Clay 0%.**

---

*Refs.* `docs/paper-v2/ch05-theorems.tex` `eq:split`, `eq:vardecomp`, `sec:wgeo` (Lemma B `Φ″⪰N/12`),
`sec:feshbach` (two-block, `G″=S−D`); `S1-CONDITIONAL-POINCARE.md` §1.2–§3 (Lemma S1, Cor. S1′), §6;
`S3-UV-FESHBACH.md` §1–§2 (Lemma S, Cor S, `G″=E[Σ]−D`); `G4-GENERAL-BOUND.md` §2 (Lemma G1/G2/G3), §3
(barrier counterexample `0.206−0.0123β²`, `β*=4.09`); `S2-BOUNDARY-COUNTEREXAMPLE.md`; `S5-N-UNIFORM.md`
§3–§4; `STAR-WAVE-SYNTHESIS.md` §4. Brascamp–Lieb (1976); Prékopa (1971); Schur-complement monotonicity.
Witness `research/v3/w7/v2_true_bound_witness.py` (numpy only, deterministic). No third-party sources
downloaded this session; all external inputs (BL, Prékopa, Schur/Haynsworth, Savvidy/Nielsen–Olesen) are
standard and cited from the corpus files.

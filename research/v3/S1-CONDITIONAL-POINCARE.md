# S1-CONDITIONAL-POINCARE — conditional Poincaré / LSI for one block with arbitrary boundary, and the boundary-dependence of (★)

> 2026-11-24 · Attack wave on the single live inequality `(★)` of the campaign, angle S1:
> the **conditional (fiber) measure of one block with arbitrary boundary `U′`**.
> Reads: `C2-BALABAN-UNIFORM.md` §2.4, `V3-WAVE-SYNTHESIS.md` §2–3, `docs/paper-v2/ch05-theorems.tex`
> (Feshbach + variance decomposition + `W_geo`), `research/unconditional/PROOF-COARSE-POINCARE.md`,
> `PROOF-L3.md`, `G1-L3-TOWER-POINCARE.md`, `G1-L3-UV-FORCE-BOUND.md`, `INVENT-VARIANCE-INEQ.md`,
> `CONDITIONAL-LATTICE-GAP.md`, `W6-BE-at-2.5.md`, `W7-Gpp0.md`.
> Witness: `research/v3/w7/s1_conditional_poincare.py` (run output quoted in §5, §6).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step / open · **[V]** numerical witness · **[R]** reported/unrefereed.
>
> **`PROOF-JACOBSON-SAVES.md` is a dangling citation in `ch05` §`sec:wgeo`.** The `W_geo≥0` statement is
> nevertheless recorded verbatim in `ch05` (lines 159–202) and its one input, `Φ''≥N/12` / `g''≥1/6`, is
> re-verified here in §5(C). The missing file does not affect this note; flagged for the citation ledger.

---

## 0 · Verdict (no QED)

**The conditional-Poincaré angle does NOT prove `(★)`. It does three things, all useful and all honest.**

1. **It localizes the entire `U′`-dependence of `(★)` into ONE matrix inequality.**
   The exact conditional (fiber) object is the Savvidy-corrected Brascamp–Lieb (SCBL) split
   `G″ = S − D` of `INVENT-VARIANCE-INEQ`, in which the fiber log-concavity defect obeys
   `D ≤ 0 ⟺ K := Π₀ Hess V Π₀ ⪰ 0` on the block. So *every* conditional route to `(★)`
   (Bakry–Émery, Brascamp–Lieb, Prékopa, log-Sobolev, Poincaré) is governed by the sign of one
   fiber-Hessian matrix `K`. **`[T]`** (identity) + **`[T]`** (BL sign law).

2. **Uniformity in `U′` is IMPOSSIBLE for that route — a clean method no-go.**
   The boundary `U′` enters `K` only through the **boundary plaquettes' Savvidy (curvature) term**,
   which on the soft (constant-transverse / `k→0`) mode is `R₀(U′) ≈ −C₁β f_∂(U′)`, with `f_∂` the
   boundary field strength. The Haar Jacobian contributes the `U′`-**independent** floor `Φ'' ⪰ N/12`
   (Prékopa, `W_geo≥0`); the block covariant Laplacian contributes a positive box gap `c_□`. Hence
   `K ⪰ 0 ⟺ c_□ + N/12 ≳ C₁β f_∂(U′)`. This **fails for every sufficiently strong boundary field**
   (`f_∂ > f_c`); the Savvidy instability is present for *any* nonzero chromomagnetic field on the
   constant transverse mode and is re-measured here (`λ_min: 0 → −1.43 → −3.07 → …` as `f: 0→0.12→0.46`).
   **`[NG]` for the named class "prove `(★)` by `U′`-uniform log-concavity / BE / BL / PI / LSI / Prékopa."**
   **`[T-standard]`** (Savvidy/Nielsen–Olesen) + **`[V]`** (witness §5B).

3. **The corrected form, and why it is exactly what the campaign already needs.**
   The `U′`-uniform statement survives on the **small-field boundary set**
   `𝒮_t := {U′ : C₁β f_∂(U′) ≤ c_□ + N/12 − ρ₀}` with an explicit modulus `ρ₀(t) > 0`
   (**Lemma S1, `[T]` given the plaquette constant `C₁`**), i.e. on precisely the boundary fields that
   **survive Lemma B's large-field removal**. The three admissible corrected targets are
   `(★_sf)` (small-field boundary), `(★_avg)` (averaged over `U′ ~ ν`), and `(★_{δ(U′)})`
   (boundary-dependent margin). **`[C]`**

**The precise obstruction, stated once.** A conditional Poincaré/LSI controls the **fiber** variance by the
**mixed** Hessian `c = Π₀ Hess V v`, i.e. it controls only the SCBL defect `D`. It does **not** control the
Gaussian–Schur part `S` — the **coarse** curvature. The free zero mode has `K ⪰ 0` (log-concave) yet
`S = 1/Var(a) → 0`, so `G″ = S − D → 0`: **log-concavity is necessary-shaped, not sufficient, for `(★)`.**
`(★)` needs `S ≥ δ` in addition, and `S≥δ` is the **coarse spectral gap = H1**. This cleanly separates the
two halves of `(★)`: *boundary-uniformity* (this file: provable on `𝒮_t`) and *coarse coercivity*
(`H1`: the gap). **`[T]` + `[V]` (witness §5D).** Clay = 0%.

---

## 1 · The conditional block measure, its Hessian, and the exact identities

### 1.1 The conditional measure

Fix a block `B` (the set of links of one coarse cell, side `≈ 2a_j` at shell `j`) and a boundary
configuration `U′` (the links adjacent to `B` but not in it, plus the coarse/IR links). The **regular
conditional law** of the block links is
```
μ_{U′}(dU_B) = Z(U′)^{-1} exp( −S_W(U_B; U′) ) ∏_{ℓ∈B} dU_ℓ ,                       (1.1)
S_W(U_B;U′) = β Σ_{p : p∩B≠∅} (1 − (1/N) Re tr U_p) ,
```
where the plaquette sum ranges over **interior** plaquettes (`p ⊂ B ∪ int`, all links in `B`) and
**boundary** plaquettes (`p` uses ≥1 link in `B` and ≥1 link in `U′`). `β = 2N/g_j²` at shell `j`.
This is the exact object that `(★)` conditions on; `Z(U′)` is the fiber partition function.

**Flat chart.** Put `U_ℓ = e^{X_ℓ}`, `X_ℓ ∈ 𝔰𝔲(N)`, on the block only (boundary held fixed). Haar is
`dU_ℓ = J(X_ℓ) dX_ℓ = e^{−Φ(X_ℓ)} dX_ℓ`, with
```
Φ(X) = Σ_{ℓ∈B} Φ_0(X_ℓ) ,   Φ_0 = Σ_{i<j} g(λ_i − λ_j) ,   g(d) = −2 log( sin(d/2)/(d/2) ) .   (1.2)
```
Crucially **`Φ` depends only on the block links, NOT on `U′`** — this is the structural fact that makes the
Haar curvature `U′`-independent. Define the total block potential
```
V(X_B; U′) := S_W(X_B; U′) + Φ(X_B) ,   μ_{U′} ∝ e^{−V} dX_B .                          (1.3)
```

### 1.2 The Hessian and where the boundary enters

The second variation splits into three pieces:
```
Hess_{X_B} V = L_{U′}(X_B)  +  R_{U′}(X_B)  +  Φ''(X_B) ,                              (1.4)
```
- **`Φ''`** — Haar Jacobian curvature; **`U′`-independent**, `Φ'' ⪰ (N/12) I` per link mode
  (`g''(d) ≥ 1/6` on `|d|<2π`; SU(2)=1/6, SU(3)=1/4). **`[T]`**, re-verified §5(C).
- **`L_{U′}`** — the lattice **covariant Laplacian** (the "curl/electric" part of the Wilson Hessian),
  self-adjoint, `⪰ 0`; by the lattice **diamagnetic inequality** its spectrum dominates the
  connection-free Dirichlet Laplacian of the block, so `L_{U′} ⪰ c_□ Π_{>0}` on the non-constant modes,
  with `c_□ = c_□(B)` the box/Dirichlet gap (in physical units `≈ c/a_j²` for the shell-`j` block).
  `U′` enters `L_{U′}` through the covariant difference at the boundary, but **only upward**
  (diamagnetism), so the *lower* bound `c_□ Π_{>0}` is `U′`-independent.
- **`R_{U′}`** — the **Savvidy (curvature)** part: the commutator cross-term of the Wilson Hessian,
  `R_{U′} ∝ β · (background field strength)`. **This is the only `U′`-dependent piece that can be
  negative.** `U′` enters `R_{U′}` through the **boundary plaquettes' field strength**; interior
  plaquettes carry the block's own (variable) field. On the constant-transverse mode,
  `R₀(U′) ≈ −C₁ β f_∂(U′)` (an explicit plaquette constant `C₁`), the classic lattice Savvidy operator
  `−D² + 2F` evaluated at `k=0`. **`[T-standard]` Savvidy 1977 / Nielsen–Olesen 1978; re-measured §5(B).**

So the only way the block measure fails to be log-concave is a boundary field strong enough to make
`R_{U′}` beat `c_□ + Φ''`. This is the whole `U′`-dependence of `(★)`, in one line.

### 1.3 The exact identities — and the two readings of `G″`

`(★)` is a statement about the **coarse** second variation. There are two equivalent bookkeepings; both
reduce to the same fiber condition `K ⪰ 0`.

**(A) RG-block reading (marginal free energy).** `ν = P_* μ`, `S_eff(U′) = −log(dν/dU′)`. For a coarse
mode `v` (lifted `Pv` to the fine block),
```
G″(U′)[v] = E_{μ_{U′}}[ ⟨Pᵀ H_S P v , v⟩ ] − Var_{μ_{U′}}( ⟨Pᵀ ∇S , v⟩ ) .               (1.5)
```
Since `Φ` is `U′`-independent, the Haar term does not appear explicitly: the free energy's Hessian is the
action curvature **minus** the fiber force variance. For the free field (`S` quadratic in the fiber with
linear coupling to the boundary) this is **exactly `0`** (translation family / Stein) — the acid test.

**(B) Single-mode SCBL reading (`INVENT-VARIANCE-INEQ`).** Fix one mode `v` (`|v|=1`), `a=⟨v,X⟩`,
condition on `a=0`; `h_n = vᵀ Hess V v = H_vv + Φ''_vv`, `f = ∂_nV = G_v + Φ'_v`. Then
```
G″(0) = E[h_n | a=0] − Var(f | a=0) = S(v) − D(v) ,                                      (1.6)
S = E[ h_n − ⟨K^{-1}c, c⟩ | 0 ] ,   D = Var(f|0) − E[⟨K^{-1}c, c⟩ |0] ,
K = Π₀ Hess V Π₀ ,   c = Π₀ Hess V v ,   Π₀ = 1 − vvᵀ  (transverse within the block) .
```
`S` is the **Gaussian–Schur / coarse** curvature (`S = 1/Var_μ(a)` for a Gaussian);
`D` is the Brascamp–Lieb defect of the conditioned fiber. The **BL sign law** is
```
D ≤ 0  ⟺  K ⪰ 0  (fiber log-concave) ;  D = 0 ⟺ Gaussian ;  D > 0 ⟺ K indefinite (Savvidy).   (1.7)
```
Here `Φ''` enters additively and, being `U′`-independent, contributes the `U′`-uniform floor
`S^Φ = 1/((Φ'')^{-1})_{nn} ≥ N/12` (`ch05`, `INVENT` §3). **`[T]`.**

**Reconciliation.** `(A)` and `(B)` are the same identity for "fiber = transverse-to-`v` within the block."
`(B)` is the form in which `(★)` is normally written (`+W_j`, the Haar curvature); `(A)` is the bare
marginal identity. In `(B)`, the `+W_j ≥ 0` term is a **genuine additive component** of `G″` (up to the
cross term `X = −2Cov(G,Φ')` that the SCBL split folds into `S`; see `ch05` `eq:split`). In `(A)` there is
no explicit `W`. **Both readings are attacked below, and the fiber condition `K ⪰ 0` is identical in both.**

**The structural point (used in §0 and §6).** In `(B)`, a conditional Poincaré/LSI/BL statement controls
**`D` only** (it needs `K ⪰ 0` and bounds the fiber variance of `f` by the *mixed* gradient `c`). It says
**nothing** about `S`, the coarse curvature. `(★)` needs `S − D ≥ δ`, i.e. **both** `D ≤ 0` and `S ≥ δ`.
`S ≥ δ` is the **coarse spectral gap** — the marginal strong log-concavity, i.e. **H1**. So *no
conditional-Poincaré argument can by itself prove `(★)`*; it can only supply the `D ≤ 0` half.

---

## 2 · The toolbox, with the boundary terms displayed

Each tool, applied to the conditional measure `(1.1)`, and exactly what `U′` does to it.

**(a) Bakry–Émery (curvature–dimension / `Γ₂`).** If `Hess V ⪰ ρ I` on the block at every fiber point
(plus `Ric ≥ ρ_G`; SU(2) `½`, SU(3) `¾`, `W6`), then `μ_{U′}` satisfies `LSI(ρ)` and `PI(ρ)`. The modulus is
`ρ(U′;X_B) = λ_min(L_{U′} + R_{U′} + Φ'')`. `U′` enters `R_{U′}`; **BE can be `U′`-uniform only if
`L_{U′}+R_{U′}+Φ'' ⪰ ρ₀ I` uniformly, which is exactly `K ⪰ ρ₀ I`.** The W6 measurement
(`Hess S_W` min `≈ −1.12β`, 21.7% negative directions at β=2.5) shows this **fails at the head**; BE lives
only for `β ≲ Ric/1.1` (`W6`), the high-temperature window. **`[T]`** for the criterion; **`[NG]`** for
`U′`-uniform at β=2.5.

**(b) Brascamp–Lieb after the flat change of variables.** For a log-concave fiber (`K ⪰ 0`),
`Var_{μ_{U′}}(f) ≤ E[⟨K^{-1}∇_⊥ f, ∇_⊥ f⟩]`, i.e. `D ≤ 0` (`(1.7)`). The boundary enters `K` exactly as in
(a). **BL = the `D`-half of SCBL; it is `U′`-uniform iff `K ⪰ 0` is.** `[T]` (BL) / `[NG]` (uniform).

**(c) Prékopa for the Haar part.** `Φ` is convex (`Φ'' ⪰ (N/12)I`, `[T]`), and `U′`-independent. By
Prékopa, every linear marginal of `e^{−Φ}dX` has convex free energy, i.e. `W_geo = E[Φ''] − Var(Φ') ≥ 0`,
uniformly in `U′` and in volume. **This is the ONLY boundary-independent positive term**, and it is exactly
the `+W_j` of `(★)`. What Prékopa gives for the *marginal of the block* requires **joint** convexity of
`e^{−S_W−Φ}`, which fails (`R_{U′}` indefinite) — so Prékopa supplies `W_geo` but not `S ≥ δ`. **`[T]`.**

**(d) Tensorization over plaquettes / links.** Tensorization would need the action to be a sum of terms
each depending on *disjoint* variables. The Wilson action's terms share links (each link lies on `2(d−1)`
plaquettes), so tensorization is unavailable. The Dobrushin/one-site criterion (a high-temperature,
strong-coupling condition) gives a gap only at `β < O(1/deg)`, and in `d=4` this is the `β < 1/48`
regime of SZZ / Lemma SC — **off the tower's path** (`β_j ∈ [2.5,∞)`). At `β≈2.5` tensorization fails.
**`[T]`** (criterion) / **`[NG]`** (at 2.5).

**Summary of §2.** The three usable tools (BE, BL, Prékopa) reduce the conditional half of `(★)` to:
`K ⪰ 0` (from BE/BL) **and** `Φ'' ⪰ N/12` (Prékopa). The former is the `U′`-dependent one.

---

## 3 · Lemma S1 — the conditional log-concavity threshold, and its `U′`-uniform corollary

> **Lemma S1 (fiber modulus and its boundary decomposition) — `[T]`.**
> Let `ρ(U′) := inf_{X_B∈(𝔰𝔲(N))^{|B|}} λ_min( Hess_{X_B} V(X_B;U′) )`. Then:
>
> **(a) [Strong log-concavity ⇒ everything conditional.]** If `ρ(U′) > 0`, the conditional measure
> `μ_{U′}` in `(1.1)` is `ρ(U′)`-strongly log-concave; hence it satisfies `LSI(ρ(U′))` and `PI(ρ(U′))`,
> Brascamp–Lieb holds (`D ≤ 0` in `(1.7)`), and the SCBL bound `G″ ≥ S ≥ N/12`-part is valid.
>
> **(b) [Boundary decomposition.]** With `(1.4)`, `c_□ = c_□(B)` the box Dirichlet gap,
> `F_∂(U′) = sup_{p∈∂B} f_p(U′)` the boundary field strength (`f_p = 1 − (1/N)Re tr U_p`), and `C₁` an
> explicit plaquette-adjacency constant,
> ```
> ρ(U′)  ≥  min{ c_□ + N/12 − C₁ β F_∂(U′) ,  N/12 + R₀(U′) } ,                          (3.1)
> ```
> where the two arguments are the bounds on the non-constant modes (Laplacian-protected) and the
> constant/soft mode (`R₀` = boundary Savvidy curvature there). `C₁`, `c_□` are `U′`-independent apart
> from the (diamagnetically favourable) connection in `L_{U′}`.

*Proof.* (a) is the standard equivalence "strong convexity of `V` ⇔ `e^{−V}` strongly log-concave ⇔ LSI/PI
with constant `ρ`," plus BL (`K ⪰ 0 ⇒ D ≤ 0`, `(1.7)`, `INVENT` §3). (b): `Hess V = L + R + Φ''`; on
non-constant modes `L ⪰ c_□ Π_{>0}` (diamagnetic inequality) and `R ⪰ −C₁βF_∂ I` (Savvidy curvature
bounded by the field strength times the plaquette count sharing a link), `Φ'' ⪰ (N/12)I`; on the constant
mode `L = 0` and `R = R₀` (no Laplacian protection), giving the second argument. Take the min over modes. ∎

> **Corollary S1′ (`U′`-uniform small-field conditional LSI) — `[C]` (given `C₁`).**
> For every `t > 0` with `C₁ β t < c_□ + N/12`, there is `ρ₀(t) := c_□ + N/12 − C₁βt > 0` such that for all
> boundary fields in the small-field set `𝒮_t := {U′ : F_∂(U′) ≤ t}`, uniformly in `U′` and in the points
> of the fiber,
> ```
> μ_{U′}  satisfies LSI(ρ₀) and PI(ρ₀), and BL ⇒ D ≤ 0 .                                    (3.2)
> ```
> In particular the `D`-half and the Haar/W part of `(★)` are `U′`-uniform on `𝒮_t`, with the explicit
> constant `ρ₀(t)`. `[C]` because `C₁` (the plaquette-geometry constant) is not reduced to a decimal —
> exactly as the campaign's `C_g ≤ CN` (`ch05` §P3) is "existence proved, value not reduced."

**Reading.** Lemma S1/Corollary S1′ say: **the conditional measure is uniformly log-concave on the
small-field boundary set, with an explicit modulus, and not on the full boundary range.** The threshold is
`C₁βF_∂ < c_□ + N/12`. Whether this covers the boundary fields that survive Lemma B is the quantitative
question of §4.

---

## 4 · Boundary-dependence: no `U′`-uniform log-concavity, counterexample family, corrected forms

### 4.1 The impossibility (method no-go)

> **Proposition U′-NG — `[T-standard]` + `[V]`.**
> There are boundary configurations `U′` for which `Hess V` is **not** positive semidefinite, so
> `ρ(U′) ≤ 0`, `μ_{U′}` is not log-concave, `K` is indefinite, `D > 0` in `(1.7)`, and the
> BE/BL/PI/LSI-constant lower bounds are `< 0` (vacuous). Concretely: take a **constant chromomagnetic
> background** on the boundary (`f_∂ > 0`). The transverse charged mode (the nonabelian `σ₁,σ₂` sector of
> SU(2)) is then the lattice Savvidy operator `−D² + 2F`; on its constant-transverse mode the Laplacian
> eigenvalue is `0`, so the curvature `−2·c₂·F < 0` makes the Hessian indefinite for **any** nonzero `F`.
> As `F` grows the negative eigenvalue grows; the Haar floor `N/12` is far too small to restore positivity.

*Mechanism, and its `U′`-monotonicity.* `R_{U′}` is **linear in the boundary field strength**; `Φ''` and
`c_□` are bounded below independently of `U′`. Hence `ρ(U′)` is **not bounded below `0` uniformly** over
`U′ ∈ SU(N)^{∂B}`. **`[T-standard]`** (the Savvidy instability is the classical one of Savvidy 1977 /
Nielsen–Olesen 1978; W4/W6 re-measure it in this campaign: `Hess S_W` min `−1.12β`, 21.7% negative
directions at β=2.5). **Witness §5(B)** reproduces the crossing on a 4D SU(2) `L=2` block:

```
  theta   f=1-cos   #neg   lambda_min   lambda_min+N/12   log-concave?
   0.00    0.0000     0       0.0000         0.1667          YES
   0.50    0.1224    16      -1.4284        -1.2618       NO (indefinite)
   1.00    0.4597    16      -3.0703        -2.9037       NO
   2.50    1.8011    19      -8.0114        -7.8448       NO
```
`λ_min` leaves `0` and plunges negative as `f_∂` grows; the Haar floor `N/12 = 0.167` cannot compensate a
curvature that is already `−1.43` (`−0.57β`) at `f=0.12`. **`[V]`.**

> **`[NG]` (named class).** *No proof of `(★)` that proceeds through `U′`-uniform conditional
> log-concavity — Bakry–Émery, Brascamp–Lieb, Poincaré, log-Sobolev, Prékopa-convexity — can succeed at
> β≈2.5, because the conditional block measure is provably not log-concave for generic boundary fields.*

**Honest scope of the `[NG]`.** This is a no-go for a **method class**, not for `(★)` itself. A
non-log-concave measure can still satisfy a Poincaré inequality with a positive constant (e.g. a symmetric
double well `e^{−(|x|²−1)²}` does). So `K ⋡ 0` refutes *these tools*, not the *inequality*.

### 4.2 Does the true conditional PI constant degrade in `U′`? (model)

The soft direction of the negative mode is a 1-D effective potential `μ_{U′} ≈ e^{−W_{U′}(y)}` on the
Savvidy direction `y`, with `W_{U′}` a **double well** once the boundary field passes the critical value
(the tachyonic quadratic is stabilized by the quartic/quenching of SU(N)). Metastability gives a Poincaré
constant `≳ e^{h(U′)}`, where `h(U′)` is the barrier height, increasing with `|f_∂ − f_c|`. Thus the
**true** conditional PI constant is finite for each `U′` but **grows without bound** as `f_∂` ranges over
`SU(N)^{∂B}` — i.e. no `U′`-uniform positive constant from this picture. **`[C]`** (model normal form is
rigorous; embedding into the exact block measure is the standing approximation).

**Existence of a `U′`-uniform PI for the true measure is `[O]`** and is not decided here; only the
log-concavity route is closed.

### 4.3 The corrected forms — and which one L3 actually needs

The campaign applies `(★)` **after Lemma B removes large-field blocks** (`C2-BALABAN-UNIFORM` §2.2). So the
boundary that actually reaches `(★)` is already **small-field** in the deep-UV / LEMMA-B-surviving regime.
The three corrected statements:

- **`(★_sf)` small-field boundary — `[C]`.**
  `(★)` restricted to `𝒮_t = {U′ : F_∂(U′) ≤ t}` with `t = t_j ≍ κ g_j`, margin `δ_W(t) > O(g_j²) > 0`.
  By Corollary S1′ the fiber is uniformly log-concave there (`D ≤ 0`), and the `W_geo ≥ N/12` term is
  `U′`-uniform; only the coarse gap `S ≥ δ` remains (H1 at the head; certified deep UV by cone-emptiness).
  **This is the form to prove for L3.**
- **`(★_avg)` averaged over `U′` — `[C]`.**
  `E_{U′~ν_j}[\text{Var} − E[H] − W] ≤ −δ_W`. The large-field boundary set that makes `K` indefinite has
  **ν-measure `≤ e^{−c/g_j²}`** (Lemma B). So the average-margin is `≤ −δ_W + O(e^{−c/g_j²})` — a
  volume-free, constant-level statement of the same species as Lemma B.
- **`(★_{δ(U′)})` boundary-dependent margin — `[C]`.**
  `δ_W(U′) = δ₀ − C₁βF_∂(U′)`, positive on `𝒮_t` by Lemma S1. Summing `1/δ_W(U′)` over the tower
  converges by Lemma B's geometric suppression — this is the "tower min" bookkeeping of
  `G1-L3-TOWER-POINCARE` §4–5, now with the correct `U′`-dependence made explicit.

**Which is needed?** For L3 (fixed-volume existence) the honest target is **`(★_sf)`**; `(★_avg)` is the
natural companion when one wants the *density* form (as in Lemma B). **`(★)` with `δ_W` uniform in **all**
`U′` is false as a log-concavity statement and must be replaced by one of these three.** `[C]`.

---

## 5 · Acid tests

### (A) Free massless field — `(★)` MUST fail. It does. `[T]` + `[V]`

Conditional on `U′`, the free Gaussian measure on the block is Gaussian (`V = ½⟨X_B, A X_B⟩` + linear
coupling to `U′`, `A` the Dirichlet Laplacian). **Stein's identity** for any Gaussian gives
`Var_fiber(force) = E_fiber[Hess]` **exactly**. `Φ`-curvature `W = 0`. Thus `(★)` reads
```
E[H] ≤ E[H] + 0 − δ_W   ⟺   0 ≤ −δ_W   FALSE.
```
and the coarse curvature `G″ = E[H] − Var = 0` (no gap). Witness §5(A): `Var = E[H] = 0.25`, difference
`0.000e+00`; the lowest Dirichlet eigenvalue `→ 0` as the block grows, so `S = 1/Var(a) → 0`. **The acid
test is passed: the tools manufacture no false positivity on the free field.** `[T]` + `[V]`.

### (B) Interacting YM at β ≈ 2.5 — non-vacuous, but not proved. `[V]` + `[O]`

- **Non-vacuity (numerical).** `W7` measures the physical `G″(0) = +0.27 … +0.51 > 0` at β=2.5, with
  `flat E[H]−Var ≈ 0` (range `−0.03…+0.04`) and a large `W_geo ≈ +0.25`; `PROVE-Gpp` reports `D/S ≈ 2%`.
  So `(★)` **holds numerically** (SCBL: `G″ = S − D ≈ 0.98 S > 0`), and `W_geo ≥ 0` supplies the margin.
  **`[V]`** — but it is the *measured* object, not a theorem.
- **Why the conditional route gives no theorem at 2.5.** The fiber Hessian is **indefinite** there
  (`W6-BE`: 90–97% of negative BE directions lie in the fiber), so `K ⋡ 0`, `D > 0`, and **BL / LSI / PI
  constants are vacuous** (they require `K ⪰ 0`). The `D ≤ 0` half is **not available** at the head. So
  the conditional-Poincaré angle is silent at β=2.5; the measured positivity is exactly the open `D < S`
  of `INVENT-VARIANCE-INEQ` §5, i.e. H1. **`[O]`.**

---

## 6 · The structural caveat, made exact (why conditional Poincaré cannot close `(★)`)

Even granting a perfect `U′`-uniform conditional Poincaré/LSI (i.e. `K ⪰ ρ₀ > 0` everywhere), the SCBL
identity `(1.6)` gives `G″ = S − D` with `D ≤ 0` but says **nothing** about `S`. A conditional LSI controls
`Var_fiber(f)` by the **mixed** gradient `c = Π₀ Hess V v`:
```
Var_fiber(f) ≤ (1/ρ) E[ |∇_⊥ f|² ] = (1/ρ) E[ |c|² ] ,   c = Π₀ Hess V v  (mixed Hessian).   (6.1)
```
It bounds the fiber force-fluctuation by the **mixed** Hessian, **not** by the **coarse** Hessian; the
coarse (`S`) curvature is a separate Schur-complement object. Witness §5(D): the Gaussian block is
log-concave (`K>0`, `D=0`) for every size, yet the `k→0` mode has `S = 1/Var(a) = 0.198 → 0.058 → 0.016 →
0.004` (`n = 6→48`) — so `G″ = S → 0`. **Log-concavity/PI hold; `(★)` still fails.** Therefore:
```
(★)  ⟺  [fiber half: K ⪰ 0, i.e. D ≤ 0]   AND   [coarse half: S ≥ δ, the marginal spectral gap].
                     ↑ this file: provable on 𝒮_t                 ↑ H1: the gap itself.
```
This is the sharpest output of the angle: **it isolates the `U′`-uniformity of `(★)` into a provable
conditional log-concavity statement on the small-field boundary (Lemma S1, Corollary S1′), and it proves
that the *remaining* content of `(★)` is exactly the coarse spectral gap — the gap, not a Poincaré
estimate.** This agrees with `INVENT-VARIANCE-INEQ` (the residual is `D<S` / fiber `H_ff`) and with
`PROOF-COARSE-POINCARE` (the effective action is Savvidy-indefinite, so not convex) — and it sharpens both
by giving the boundary decomposition `(3.1)`, the threshold `C₁βF_∂ < c_□ + N/12`, and the corrected forms.

---

## 7 · Verdict and what it does to the campaign

**Verdict.** For the angle "conditional Poincaré / LSI for one block with arbitrary boundary":

1. **Setup and Hessian — done, explicit** (`§1`): the conditional block measure `(1.1)`, the flat-chart
   split `(1.3)`, the Hessian decomposition `(1.4)`, and the exact identities `(1.5)`–`(1.6)`. The boundary
   `U′` enters **only** through the Savvidy curvature `R_{U′}` of the boundary plaquettes. **`[T]`**
2. **Attempted lemma — Lemma S1 / Corollary S1′ — proved** (`§3`): the conditional modulus
   `ρ(U′) ≥ min{c_□ + N/12 − C₁βF_∂, N/12 + R₀}`, and uniform strong log-concavity / LSI / BL on the
   small-field boundary set with explicit `ρ₀(t)`. **`[T]`/`[C]`** (C₁ not decimalised).
3. **`U′`-uniformity is impossible for the method, with a counterexample family** (`§4.1`): chromomagnetic
   boundary fields make `K` indefinite for any `f_∂ > 0` on the soft mode; re-measured (`λ_min` to `−8`).
   **`[NG]` for the log-concavity/BE/BL/PI class.**
4. **Corrected forms stated** (`§4.3`): `(★_sf)` (small-field boundary — the L3-relevant one),
   `(★_avg)`, `(★_{δ(U′)})`. **`[C]`**
5. **Acid tests** (`§5`): free field fails `(★)` exactly (`Var=E[H]`, `G″=0`); β≈2.5 is non-vacuous
   (`G″>0` measured, `D/S≈2%`) but unproved. **`[T]`+`[V]` / `[V]`+`[O]`.
6. **The structural limit** (`§6`): conditional Poincaré controls the *mixed* Hessian (`D`), not the
   *coarse* Hessian (`S`); `(★)` needs **both**, and `S ≥ δ` is H1. **`[T]`+`[V]`.

**Net effect on the campaign.** The conditional-Poincaré angle **does not close `(★)` and produces no
no-go for `(★)` itself.** It (i) proves the `U′`-uniform conditional log-concavity/LSI *on the small-field
boundary set that survives Lemma B*, with an explicit modulus — so it **discharges the boundary-uniformity
half of L3 in the regime where Lemma B applies**; and (ii) proves that the remaining content of `(★)` is
**exactly the coarse spectral gap = H1**, not a conditional variance estimate. Concretely, L3 reduces from
`(★)∧(★dec)` to:
- **`(★_sf)`** — small-field conditional log-concavity + coarse coercivity, the first half now provided by
  Lemma S1 there; and
- the **head** `β≈2.5`, where the conditional measure is not log-concave and `(★)` is the measured-but-
  unproved `D<S` of `INVENT`, i.e. **H1**.
That is a genuine sharpening: the UV end of L3 is now computable with a named threshold
(`C₁β F_∂ < c_□ + N/12`), and the IR end is confirmed to be the gap. **No QED; Clay 0%.**

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| conditional block measure + flat-chart Haar split | eqs `(1.1)`–`(1.3)`; `Φ` is `U′`-independent | **`[T]`** |
| exact identities (marginal + SCBL) | `G″ = E[H]−Var`; `G″=S−D`; BL sign law `D≤0 ⇔ K⪰0` | **`[T]`** (carried from `INVENT`) |
| **Lemma S1** fiber modulus + boundary decomposition | `ρ(U′) ≥ min{c_□+N/12−C₁βF_∂, N/12+R₀}` | **`[T]`** |
| **Corollary S1′** uniform small-field conditional LSI/BL | `LSI(ρ₀(t))`, `D≤0` on `𝒮_t`, `ρ₀=c_□+N/12−C₁βt` | **`[C]`** (C₁ not decimalised) |
| `Φ'' ≥ N/12`, `W_geo ≥ 0` | Prékopa; `U′`-independent | **`[T]`** (re-verified §5C) |
| **Prop U′-NG** uniform log-concavity fails | chromomagnetic `U′` ⇒ `K` indefinite for `f_∂>0` (soft mode) | **`[T-standard]`+`[V]`** |
| **`[NG]` class** `(★)` via `U′`-uniform BE/BL/PI/LSI/Prékopa | closed at β≈2.5; log-concavity refuted | **`[NG]`** |
| log-concavity/PI **not sufficient** for `(★)` | controls `D` (mixed Hessian), not `S` (coarse) | **`[T]`+`[V]`** |
| corrected forms `(★_sf)`, `(★_avg)`, `(★_{δ(U′)})` | small-field / averaged / boundary-dependent margin | **`[C]`** |
| `(★)` true at β≈2.5 | `G″=+0.27…+0.51`, `D/S≈2%` (W7/PROVE-Gpp) | **`[V]`** |
| `(★)` proved | uniform in `j,U′,a` | **`[O]`** |
| free-field acid test | `Var=E[H]`, `G″=0`, `(★)` false | **`[T]`+`[V]`** |

**Sharpest line.** The conditional-Poincaré angle turns the `U′`-dependence of `(★)` into one matrix
inequality, `K = Π₀ Hess V Π₀ ⪰ 0`, governed by the competition
`(block Laplacian gap + Haar floor N/12) vs (boundary Savvidy curvature C₁βF_∂)`. It is **provable**
(uniform in `U′`) exactly on the small-field boundary set that Lemma B leaves behind — and **impossible**
(full range) because the boundary can realize the Savvidy instability. Crucially it also proves that even a
perfect conditional Poincaré controls only the **mixed**-Hessian defect `D`, never the **coarse** curvature
`S`; so `(★)` remains `[K⪰0 on small field]` **plus** `[S≥δ = H1]`. The conditional-variance half is
bankable in the UV; the gap is still the gap. **No QED.**

---

## 9 · Witness (`research/v3/w7/s1_conditional_poincare.py`)

Deterministic, numpy-only, runs under `py -3.12` (numpy 2.4.6), exit 0. Four parts:
- **(A)** free-field Dirichlet chain: `Var_fiber(force)=E_fiber[Hess]=0.25` (difference `0`), `c_n→0`;
  `(★)` fails, as required.
- **(B)** 4D SU(2) `L=2` fine Wilson Hessian at a U(1) constant field `f=1−cosθ`; `λ_min: 0 → −1.43 →
  −3.07 → … → −8.01` (`#neg` 0→16..19); Haar floor `N/12=0.167` cannot compensate. Savvidy threshold.
  (Caveat: `L=2` global soft mode, no block Dirichlet gap; a finite block's `c_□ ~ 1/a_j²` restores
  log-concavity deep UV — the UV/IR split.)
- **(C)** SU(2) Jacobian `g''(d) = ½csc²(d/2) − 2/d² ≥ 1/6` on `|d|<2π` (min `= 0.16667`); `W_geo≥0`
  by Prékopa, `U′`-independent.
- **(D)** Gaussian block: `D=0`, `S=1/Var(a)`; `k→0` gives `S = 0.198→0.058→0.016→0.004` — log-concavity
  holds for every size yet `G″→0`: log-concavity is not sufficient for `(★)`.

*Refs.* `C2-BALABAN-UNIFORM` §2.4; `V3-WAVE-SYNTHESIS` §2–3; `ch05-theorems.tex` §`sec:wgeo`, §`sec:feshbach`,
`eq:split`; `PROOF-COARSE-POINCARE` (fiber `H_ff` indefinite; Prékopa fails jointly); `PROOF-L3`; `G1-L3`;
`G1-L3-TOWER-POINCARE` §2–5; `G1-L3-UV-FORCE-BOUND` §1,§3; `G1-UNIFORM-ZA` §3; `INVENT-VARIANCE-INEQ`
(SCBL `G″=S−D`, `D≤0⇔K⪰0`); `CONDITIONAL-LATTICE-GAP` (H1/H2, `Ric≥ρ_G`); `W6-BE-at-2.5` (`Hess S_W` min
`−1.12β`, 21.7% neg); `W7-Gpp0` (`G″+0.27…+0.51`, flat `E[H]−Var≈0`); Savvidy 1977; Nielsen–Olesen 1978;
Brascamp–Lieb 1976; Prékopa 1971; Bakry–Émery; Holley–Stroock Ann. Probab. 15 (1987) 1541; Shen–Zhu–Zhu
CMP 400 (2023) 805, arXiv:2204.12737 (`β<1/48`, Lemma SC); Ledoux (LSI pushforward).

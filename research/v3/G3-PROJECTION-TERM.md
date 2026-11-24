# G3-PROJECTION-TERM — does `G'' = 2W_geo` survive the coarse projection `P`?

> 2026-11-24 · Attack surface G3. Generalizes the exact S2 one-plaquette identity
> `G''(v) = 2W_geo(v)` from a *single fixed fine direction* `v` to the *coarse mode* `v = P w`.
> Reads: `S2-BOUNDARY-COUNTEREXAMPLE.md` (the identity + its per-direction proof),
> `S3-UV-FESHBACH.md` (Schur/variance split `G'' = E[Σ] - D`, the defect `D`),
> `S5-N-UNIFORM.md` (N/a bookkeeping, boundary-plaquette counting),
> `docs/paper-v2/ch05-theorems.tex` eqs. (eq:wgeo, eq:split, eq:vardecomp).
> Witness: `research/v3/w7/g3_projection.py` (deterministic 3-D quadrature; `py -3.12`).
> *(G1-BLOCK-IDENTITY.md and G2-TWO-PLAQUETTE.md do not exist in the repo; this note builds
> directly on S2/S3/S5.)*
>
> **Labels (charter §2).** **[T]** true theorem (and machine-checked here) · **[C]** conditional /
> proof-sketch · **[V]** numerical (deterministic quadrature; identity residual quoted) ·
> **[NG]** no-go for a named class · **[O]** exact missing step.

---

## 0. Verdict (no QED)

**The identity survives a constant-coefficient (flat-chart-linear) projection exactly, and it
acquires a single, explicit, sign-indefinite extra term for a covariant (parallel-transport)
projection.** Precisely:

1. **Case (a) — linear `P` in the flat chart. [T]** If the coarse mode pulls back to a *fixed*
   fine direction `w = Pv` (constant coefficients: `(PX)_c = sum_l P_{cl} X_l`, `X_l = log U_l`),
   then `d_w = sum_l (Pv)_l d_{X_l}` is an ordinary directional derivative along a constant
   vector, the S2 IBP goes through verbatim, and
   ```
   G''(w) = E[H_P] - Var(G_P) + W_geo(w) = 2 W_geo(w)        EXACTLY,
   ```
   with `W_geo(w) = sum_l w_l^2 E[Phi''_l] - Var(sum_l w_l Phi'_l)`. **No obstruction.** The
   "per-link" structure is real but does *not* make `G''` a sum of per-link identities (Var is
   not additive): the curvature part `E[Phi''_w]` is a per-link sum, the variance part
   `Var(Phi'_w)` carries the cross-link correlations. Witness `s=0`: `G'' = 2W_geo` to `2e-7`,
   `R_proj = 0` to machine precision.

2. **Case (b) — covariant `P` (parallel transport). [T] exact identity.** The coarse direction
   pulls back to a *vector field* `w(X)_l = Ad_{Omega_l(X)^{-1}} v` (transporter `Omega_l` from
   link `l` to the block base). The S2 lemmas generalise to vector-field Stein identities
   ```
     (L2')  E_mu[ d_w V ]  =  E_mu[ div w ]          (= int div(Pw) e^{-S} J dX) ,
     (L1')  Var_mu( d_w V ) = E_mu[ w^T Hess V w ] + E_mu[ <(Dw)w, grad V> ] + Cov_mu( d_w V , div w ) ,
   ```
   with `V = S + Phi`. From these, **exactly**,
   ```
   G''_P = 2 W_geo(w) + R_proj ,
   R_proj = E[<(Dw)w, grad Phi>] - E[<(Dw)w, grad S>]
            + Cov(<w,grad Phi>, div w) - Cov(<w,grad S>, div w) .
   ```
   Equivalently `R_proj` is the "covariance form" of the task's `int div(Pw) e^{-S} J`; **its mean
   part is exactly the failure of S2 Lemma 2** (`E[d_w V] = E[div w] != 0`). Machine-verified to
   `<= 1.1e-7` at every scanned `(beta, theta, v, s)`.

3. **Sign — [V] sign-indefinite, [NG] for "controlled by `W_geo` at the head".** `R_proj` is
   negative at the worst head configuration (boundary at maximal holonomy `theta=pi`, the Savvidy
   regime) but positive in other directions/sizes:
   ```
    beta=2.5, theta=pi,  v=e3, s=0.3 :  2W_geo=+0.2618 ,  R_proj=-0.1143 , G''=+0.1475
    beta=5.0, theta=pi/2,v=e1, s=0.7 :  2W_geo=+0.3064 ,  R_proj=+0.2549 , G''=+0.5613
    beta=0.5, theta=pi,  v=e1, s=0.7 :  2W_geo=+0.00759,  R_proj=-0.00291, G''=+0.00469 (smallest)
   ```
   `|R_proj|/(2W_geo)` reaches **0.76** at the head (`beta=5, theta=pi, e1, s=0.7`) — the extra
   term is **the same order as the floor**, not a small correction. The scanned `G''_P` stays
   **positive** everywhere (no counterexample to (★)); the floor shrinks but does not cross `0`.

4. **Control — [C] in the UV tower, [O]/[NG] at the head.** In the UV `g_j -> 0`, both `div w` and
   `Dw` are `O(g_j)` (they are commutator/Maurer-Cartan objects that vanish when `w` becomes
   constant, i.e. `g_j=0`), while `W_geo = O(1)` is coupling-free. Hence
   ```
   R_proj = O(g_j)  ,   |R_proj| / W_geo = O(g_j) -> 0      (deep UV tower, conditional on the
                                                             moment bootstrap E|F_j| <= C g_j)
   ```
   At the **head** (`beta ~ 2.5`, `g ~ O(1)`) this suppression is gone: `R_proj = O(W_geo)`,
   sign-indefinite, **not bounded by `W_geo` alone** ([NG] for that class). A genuine head bound
   is [O] and joins `(★dec)`/H1.

5. **Acid test — [T].** Free flat field (`U(1)`, `Phi=0`, `J=1`): `W_geo=0`, `R_proj=0` (adjoint
   trivial => `w` constant => `div w = 0`, `Dw=0`), so `G''_P = 0`. **(★) fails, `delta_W=0`.** On
   a Gaussian free field the witness returns `G''_P = 1e-15`. The projection manufactures no floor
   where the free theory has none.

**Sharpest line.** A coarse mode `v = Pw` is answered in two regimes: **linear (flat-chart) block
averaging leaves S2's `G'' = 2W_geo` exact — the identity is a directional statement and a constant
`P` is just another fixed fine direction** — while **covariant (parallel-transport) block averaging
adds one explicit term `R_proj`, the divergence/connection of the transported direction, which is
provably `O(g_j)`-small in the UV tower but `O(W_geo)` and sign-indefinite at the IR head.** The
projection does **not** break the UV shell inequalities `(★)_j`; it moves the burden to the head,
exactly where `(★dec)` already lives. Clay 0%.

---

## 1. Setup, conventions, and what "projection" means

Work in the flat chart `X_l = log U_l in su(N)` with Haar `dU = J dX`, `Phi = -log J`, full
potential `V = S + Phi`, measure `dmu ∝ e^{-V} dX`. The S2/S5 convention (kept here, and the one
for which `G''=2W_geo` is exact): the force and Hessian in `G''` come from the **action** `S`,
the geometric curvature `W_geo` from the **Jacobian** `Phi`:
```
G_v = d_v S ,  H_vv = d_v^2 S ,  W_geo(v) = E[Phi''_vv] - Var(Phi'_v) ,
G''(v) = E[H_vv] - Var(G_v) + W_geo(v) .
```
(★) is `G''_j(v) >= delta_W > 0`, uniform in the environment `U'` and mode `v`.

The coarse mode is `v = P w`. Two readings of "the projected mode":

* **Reading 1 (direction / fibre — the campaign's, per S2).** `w = Pv` is the *fine* direction
  along which (★) is evaluated; the measure stays the fine fibre `mu(.|U')`. S2's own line: "★ uses
  `langle P^T grad S, v rangle = d_w S` with `w = Pv`, a single fixed direction". This is the
  reading the task's cases (a)/(b) probe (constant `P`, then covariant `P`).
* **Reading 2 (pushed-forward marginal — `nu = P_* mu`).** Then one wants the curvature of the
  *coarse effective action*. This is a **different** object (`nu`, not the fibre) and acquires a
  *law-of-total-variance* Fisher/entropy term, not the `R_proj` below. It is flagged as [O] in §8,
  because the campaign's `(★)` is Reading 1.

Throughout, `d_w f = <w, grad f> = sum_i w_i d_i f` is the pairing with the gradient (a scalar
directional derivative); `w^T Hess f w = sum_{ij} w_i (d_i d_j f) w_j` is the quadratic form (no
derivative of `w`); `div w = sum_i d_i w_i`; `(Dw)_{ij} = d_j w_i` is the Jacobian of the vector
field; `(Dw)w` is its action on `w`.

**Boundary/edge.** `e^{-V} = e^{-S} J` vanishes on the SU(N) cut locus (the flat-chart edge,
`rho = 2 pi` per link), and `V -> +inf` there (the Haar density `J -> 0`). Every IBP below drops a
flux `f w e^{-V}` that vanishes on that edge. This is the *only* global input; it is why S2's
lemmas are true and why they extend to vector fields.

---

## 2. Case (a): constant-coefficient `P` — exact, **[T]**

Let `(PX)_c = sum_l P_{cl} X_l` with **constant** `P_{cl}` (the flat-chart linear block average).
The coarse Hessian/projected force of (★) are
```
<P^T grad S, v> = sum_l (Pv)_l d_{X_l} S = d_w S ,   w := Pv  (a constant fine vector) ,
<P^T H P v, v>  = w^T H_S w = H_P .
```
`w = Pv` does not depend on `X`. So `d_w` is an ordinary partial derivative in the constant
direction `w`, and `H_P = d_w^2 S` with no `(Dw)` term. **S2 §2 applies verbatim** with the single
fixed direction `w`:

> **Proposition (a). [T]** For every constant-coefficient `P`, every `beta, U'` and every coarse
> direction `v`, `G''(Pv) = 2 W_geo(Pv)` exactly, with
> `W_geo(Pv) = E[(Pv)^T Hess Phi (Pv)] - Var(<Pv, grad Phi>)`.

*Proof.* `d_w V = d_w(S+Phi)`; `e^{-V}` vanishes at the chart edge; `int d_w(d_w V e^{-V}) = 0`
gives `Var(d_w V) = E[d_w^2 V]`; `int d_w e^{-V} = 0` gives `E[d_w V] = 0`; subtracting the
`Phi`-diagonal reproduces S2's corollary with `v -> w`. ∎

**The "per-link" refinement (what is and is not additive).** `Phi = sum_l Phi_l(X_l)` is a sum of
single-link potentials, so `d_w` acts link-by-link:
```
Phi'_w   = sum_l w_l Phi'_l(X_l)          (first derivative: additive) ,
Phi''_w  = sum_l w_l^2 Phi''_l(X_l)       (second derivative: DIAGONAL, no cross terms) .
```
Hence the curvature half of `W_geo` is a per-link sum, but the variance half is not:
```
W_geo(w) = sum_l w_l^2 E[Phi''_l]  -  Var( sum_l w_l Phi'_l ) .
```
`Var` of the sum carries all the cross-link correlations `Cov(Phi'_l, Phi'_m)`. So the IBP is
"per-link" in the sense that each `Phi_l` contributes its own term and each plaquette its own
`S_p`, but `G''(Pv)` is **not** `sum_l G''_l`; it is one identity for the aggregate direction `w`.
This is exactly why the identity survives any constant `P`: `P` only chooses *which* fixed direction
`w` to look along; it never forces you to integrate out a fine direction. Witness confirms `s=0`:
`G'' = 2W_geo` (e.g. `beta=2.5, theta=pi, e3`: `G''=0.2618240`, `2W_geo=0.2618240`), residual
`R_proj <= 2e-8`.

**Obstruction in case (a)?** None in Reading 1. (In Reading 2 there is one — see §8.) Moverover the
argument needs `w` supported on the *integrated* (fibre-interior) links; if `P` also moves the
*conditioned* boundary links `U'`, project onto the interior first — the boundary is pinned, so its
contribution to the fibre derivative is zero. This is automatic for a block-interior average.

---

## 3. The generalized (vector-field) Stein identities

For a smooth vector field `w(X)` and density `e^{-V}` vanishing on the chart edge:

**Lemma 2' (mean force). [T]**
```
  E_mu[ d_w V ] = E_mu[ div w ] .
```
*Proof.* `div(w e^{-V}) = div(w) e^{-V} - <w,grad V> e^{-V}`; integrate, the boundary flux of
`w e^{-V}` is zero. ∎  For constant `w` (`div w=0`): `E[d_w V]=0`, S2 Lemma 2. **The mean part of
`R_proj` is exactly `E[div w]` — the total-derivative density `int div(Pw) e^{-S} J` the task
names.** (Note `e^{-S}J = e^{-V}`, so `int div(Pw) e^{-S}J = E_mu[div w] = E_mu[d_w V]`.)

**Lemma 1' (second order). [T]**
```
  Var_mu( d_w V ) = E_mu[ w^T Hess V w ] + E_mu[ <(Dw)w, grad V> ] + Cov_mu( d_w V , div w ) .
```
*Proof.* `int d_w[d_w V e^{-V}] = -int d_w V e^{-V} div w` (total derivative along `w`; boundary
flux zero). LHS `= int[(d_w^2 V) - (d_w V)^2]e^{-V}` with
`d_w^2 V = <w,grad(d_w V)> = w^T Hess V w + <(Dw)w, grad V>`. Rearrange; subtract `E[d_w V]^2` and
use Lemma 2' to turn `E[d_w V . div w] - E[div w]^2` into `Cov(d_w V, div w)`. ∎  For constant `w`
(`Dw=0`, `div w=0`): `Var(d_w V)=E[w^T Hess V w]`, S2 Lemma 1.

These two are the whole content. Everything below is algebra.

---

## 4. Case (b): covariant `P` — the explicit extra term, **[T]**

Model the covariant coarse-graining (parallel transport to a block base point):
```
X_l = log U_l ,   w(X)_l = Ad_{Omega_l(X)^{-1}} v = Omega_l^{-1} v Omega_l ,   l in block ,
```
with `Omega_l` the ordered product of links from `l` to the base (a path of length `~` block
diameter). `w(X)` is a **vector field**: it vanishes as a rotation when `X -> 0` (so `w -> v` only
in the trivial-field limit), and `Dw`, `div w` are the Maurer-Cartan/structure-constant objects
that measure how the transported direction turns.

Project the definitions of §1 along `w`, and split `d_w V = d_w S + d_w Phi =: G + phi`, with
`H := w^T Hess_S w`, `Psi := w^T Hess_Phi w`. Define the projected geometric curvature
`W_geo(w) := E[Psi] - Var(phi)` — the exact same functional as S2's, evaluated on the vector field
`w`. Substituting the action/Phi split into Lemma 1':

**Theorem (projection identity). [T]**
```
  G''_P  :=  E[H] - Var(G) + W_geo(w)  =  2 W_geo(w) + R_proj ,
  R_proj = E[<(Dw)w, grad Phi>] - E[<(Dw)w, grad S>]
           + Cov(<w,grad Phi>, div w) - Cov(<w,grad S>, div w) .
```
*Proof (algebra from §3).* Lemma 1' on `V=S+Phi`:
`Var(G+phi) = E[H+Psi] + E[<(Dw)w,grad V>] + Cov(G+phi, div w)`. Insert `Var(G+phi) =
Var G + Var phi + 2Cov(G,phi)`, `E[H+Psi]`, and `W_geo=E[Psi]-Var phi`; use the cross-Stein
`Cov(G,phi) = W_geo + E[<(Dw)w,grad Phi>] + Cov(phi, div w)` (from `int div(phi w e^{-V})=0` plus
Lemma 2'), and `E[G]+E[phi]=E[div w]`. Collect. ∎

**Three readings of the same `R_proj`.**
- **Divergence form (task's form).** The mean piece is `E[div w] = int div(Pw) e^{-S}J` (Lemma 2'),
  the failure of S2 Lemma 2 for a covariant direction. The second-order piece is the covariance of
  the projected force with `div w`, `Cov(G+phi, div w)`.
- **Connection form.** The `E[<(Dw)w, grad V>]` terms survive *even when `div w = 0`*; they are the
  anholonomy of the transported direction, proportional to `Dw`. The witness shows both are active:
  for `v=e1` (where `psi`-symmetry forces `E[div w]=0`) `R_proj` is nevertheless nonzero, carried
  entirely by the `(Dw)` terms.
- **S2 limit.** Constant `w` (`Dw=0`, `div w=0`) gives `R_proj=0` and `G''_P = 2W_geo`, i.e. case
  (a). The witness `s -> 0` scan gives `R_proj -> 0` (like `O(s^2)` for the axial transport here —
  the leading `O(s)` is a total derivative and integrates out), recovering Proposition (a).

**Verification.** `g3_projection.py` computes `G''_P`, `W_geo(w)`, `R_num = G''_P - 2W_geo` and
`R_ana` (the formula above, with `Dw` by finite differences) by deterministic 3-D quadrature in
`(rho, u=cos theta, psi)`. Across all 72 scanned points `|R_num - R_ana| <= 1.1e-7`; representative:

| `beta` | `theta` | `v` | `s` | `2W_geo` | `R_num` | `R_ana` | `G''_P` |
|---|---|---|---|---|---|---|---|
| 2.50 | pi | e3 | 0.30 | 0.261824 | -0.114291 | -0.114291 | 0.147533 |
| 2.50 | pi | e1 | 0.70 | 0.131731 | -0.064903 | -0.064903 | 0.066828 |
| 2.50 | pi/2 | e3 | 0.70 | 0.274987 | +0.014777 | +0.014777 | 0.289764 |
| 5.00 | pi/2 | e1 | 0.70 | 0.306448 | +0.254863 | +0.254863 | 0.561311 |

`R_proj = 0` **exactly** whenever `theta = 0` (boundary `B = I`): full axial symmetry makes the
transport a symmetry of `mu`, so `div w` and the `(Dw)`-contractions are isotropic and integrate
out. The extra term is generated by **boundary holonomy breaking isotropy**, not by the projection
per se.

---

## 5. Sign and size — **[V]**, **[NG]** at the head

* **Sign-indefinite.** From the table: `R_proj = -0.114` at the targeted worst case
  (`theta=pi`, transverse-ish alignment) and `+0.255` at `theta=pi/2, e1, s=0.7`. It flips sign
  under changes of direction `v`, path size `s`, and boundary `theta`. **So `G''_P = 2W_geo +
  R_proj` is not automatically `>= 2W_geo`: at `theta=pi` the genuine floor is `2W_geo - |R_proj|`.**
* **Not small at the head.** `|R_proj|/(2W_geo)` reaches **0.76** (`beta=5, theta=pi, e1, s=0.7`),
  with typical values `0.4-0.5` at `O(1)` coupling. The extra term is the same order as the floor.
  (At `beta=0.5` the ratio is `0.2-0.4`, and the smallest observed `G''_P` is `0.00469` at
  `beta=0.5, theta=pi, e1, s=0.7` — still positive, but only `62%` of `2W_geo`.)
* **Simplest negative configuration.** Boundary at maximal holonomy `theta = pi` (`B -> -I`, the
  Savvidy/maximal-background regime) with the direction roughly aligned to the transport axis
  (`e3`), `s ~ 0.3`: `2W_geo = +0.2618`, `R_proj = -0.1143` (a `44%` bite), `G''_P = +0.1475 > 0`.
* **No crossed floor in the scan [V].** `G''_P > 0` at all 72 sampled points; the toy finds no
  counterexample to (★). But because `|R_proj| = O(W_geo)` at the head, a *generic* covariant
  transporter (the toy uses an axial one) is **not excluded** from crossing `0`. This is [O] and is
  the head burden.

**Size / block scaling.**
- `div w`, `Dw` vanish as `P -> const` (path length -> 0): `R_proj -> 0`, recovering case (a). The
  toy's `s` is a path-length proxy; `R_proj = O(s^2)` for the axial transport there (leading `O(s)`
  integrates out), generically `O(s)`.
- `E[div w] = int div(Pw) e^{-V}` is a **total-derivative / surface** quantity: for a translation-
  invariant density it would be `∮ w.n = 0` (w supported in the block interior, vanishing at its
  support edge), so it is a *boundary-correlation* with `grad e^{-V}`, i.e. a boundary effect. In
  lattice units it scales with the block **surface** `~(L_b/a)^{d-1}`, which is exactly S5's
  "boundary-plaquette" row normalized to `O(1)` per coarse mode. It is **not** an extensive
  `O(L_b^d)` bulk term. (The toy is 0-dimensional and cannot resolve `L_b`; the surface reading is
  the structural one.)

---

## 6. Control: UV tower safe, head open — **[C]** / **[O]**

The decisive point is the **coupling scaling** of `R_proj`'s ingredients versus `W_geo`:

| object | scaling | why |
|---|---|---|
| `w` vs constant | `w - v = O(g_j)` | `Ad_{Omega^{-1}} = 1 + O(g_j)`: rotation angle `~` transported field `~ g_j` |
| `Dw`, `div w` | `O(g_j)` | Maurer-Cartan / commutator of the transported direction; `= 0` at `g_j = 0` |
| `grad S` | `O(beta_j g_j)` | Wilson force; `beta_j = 4N beta_0 j` |
| `grad Phi` | `O(1)` | Haar force, coupling-free |
| `W_geo = E[Psi]-Var(phi)` | `O(1)` (indeed `>= c N`) | coupling-free Haar geometric curvature (S5: `W_geo/N ~ 0.125`) |
| `R_proj` | `O(g_j)` | every term carries one `Dw` or `div w` factor `= O(g_j)` |

Hence
```
  R_proj = O(g_j) ,    |R_proj| / W_geo = O(g_j) -> 0     (g_j = 1/sqrt(2 beta_0 j) -> 0 in the tower)
```
**So in the deep UV shell `(★)_j` is safe up to an explicit `O(g_j)` correction:** writing
`G''_{P,j} = 2W_geo(w)(1 + O(g_j))`, the certified floor `delta_W` survives with an `O(g_j)`-small
loss. This is **[C]**, conditional on the same gauge-invariant moment bootstrap `E|F_j| <= C g_j`
(=`(H-C)`) that already underlies Route B — the bound on the covariances `Cov(., div w)` needs
`Var(div w), Var(G)` controlled at `O(g_j^2)`, which is exactly that bootstrap. It is **not** a new
obstruction.

At the **head** (`beta ~ 2.5`, `g ~ O(1)`; and at the `beta=5` toy points where the saddle sits at
`rho=pi` so `w` is a large rotation), the suppression is absent, `R_proj = O(W_geo)`, and there is
**no bound of `R_proj` by `W_geo` alone** — it is sign-indefinite and the same size. **[NG]** for
"controlled by `W_geo`". A genuine head bound on `|R_proj|` (a clustering/`(★dec)`-type estimate of
the force-divergence correlation on the boundary-conditioned fibre) is **[O]** and joins H1.

> **Net.** The projection does **not** threaten the UV tower — it adds a provably `O(g_j)`-small
> correction there. It acts like a *boundary/holonomy* effect concentrated at the IR head, where it
> is `O(W_geo)` and indefinite, hence part of `(★dec)`/H1, not a new species of obstruction.

---

## 7. Acid test — **[T]**

Free flat field: `U(1)` (or `R`), `J = 1`, `Phi = 0`, `S` Gaussian. Then `W_geo = 0`. For the
abelian group the adjoint action is trivial, so the covariant direction degenerates to a constant
`w` (`Dw = 0`, `div w = 0`) and `R_proj = 0`. By Gaussian equipartition `Var(G) = E[H]`, so
`G''_P = 0`: **(★) fails with `delta_W = 0`.** Required. Witness: Gaussian free field, `Phi=0`,
`w = Ad_{exp(-sX)}v`, returns `G''_P = 1e-15` for all `a in {2,8}`, `v in {e3,e1}`, `s in {0,0.7}`.
The covariant correction adds no `Phi`-independent positivity — it cannot conjure a floor that the
free theory lacks.

---

## 8. What this means for (★), and the one thing left outside

* **Linear `P` (case a) is free.** If the campaign's coarse-graining is a flat-chart linear
  average (or if one works, as S2 does, per fixed fine direction `w=Pv`), the identity is exact for
  every `P`, `(★) ⟺ 2W_geo(w) >= delta_W`, and everything reduces to S2's positivity of the
  interacting Haar curvature — no projection term.
* **Covariant `P` (case b) costs one explicit term.** The genuine coarse-graining of lattice gauge
  theory uses parallel transport (covariant block spin), so case (b) is the physically relevant
  one. The correction `R_proj` is exact, explicit, and **small in the UV tower** (`O(g_j)`), which
  is where the shell inequalities live. `(★)_j` survives with an `O(g_j)` correction. At the head
  it is `O(W_geo)` and indefinite — the head burden grows by one named term, `(★)` itself is neither
  proved nor broken.
* **Reading 2 (marginal) is a different extra term — [O].** If instead one wants the curvature of
  the *pushed-forward coarse measure* `nu = P_* mu` (the coarse effective Hessian after integrating
  out the fibre), the law of total variance gives an additional `-E_nu[Var_mu(d_w V | Y)]`-type
  Fisher/entropy term, which is **not** `R_proj`. The campaign's `(★)` (per S2/S5) is Reading 1
  (fibre at fixed `U'`, not the marginal), so this term is out of scope here; whether Reading 2's
  object is what Clay ultimately needs is [O] and separate. (It is the same door S2 §6 and S3's
  defect `D` already flagged: an interior coarse constraint / decimation defect beyond the
  boundary-conditioned fibre.)

**Relation to S3.** S3's exact split is `G'' = E[Σ] - D`, `D = Var(f) - E[<K^{-1}c,c>]` the
Brascamp-Lieb fibre defect, with **constant** `P` implicit. `R_proj` is a *third*, independent
correction that appears only when `P` is covariant — it is disjoint from `D` (S3's `D` measures
fibre log-concavity; `R_proj` measures the divergence of the transported direction). For constant
`P`, `R_proj = 0` and S3's `D` is the whole story; for covariant `P`, both are present and must be
bounded at the head. Neither is controlled by the other.

---

## 9. Status ledger

| item | statement | label |
|---|---|---|
| Case (a): constant-coefficient `P` ⇒ `G''(Pv) = 2W_geo(Pv)` exactly | global IBP along fixed `w=Pv`; witness `s=0` residual `2e-7` | **[T]** (+[V]) |
| "per-link" structure: `Phi''_w` diagonal (per-link sum), `Var(Phi'_w)` carries cross-link correlations | `G''` is *not* `sum_l G''_l` | **[T]** |
| Case (b): exact `G''_P = 2W_geo(w) + R_proj`, explicit `R_proj` | vector-field Stein L1'+L2'; verified `<=1.1e-7` | **[T]** (+[V]) |
| `E[d_w V] = E[div w] = int div(Pw) e^{-S}J` (mean part = failure of S2 Lemma 2) | Lemma 2' | **[T]** |
| `R_proj` sign-indefinite; `<=0.76 * 2W_geo` at the head | explicit `+`/`-` points | **[V]** |
| `R_proj = O(g_j)`, `R_proj/W_geo = O(g_j) -> 0` in the UV tower | `Dw,div w = O(g_j)`, `W_geo=O(1)` | **[C]** (on `E|F_j|<=C g_j`) |
| `R_proj` bounded by `W_geo` at the head | **false**: `R_proj = O(W_geo)`, indefinite | **[NG]** |
| a genuine head bound on `|R_proj|` (force-divergence clustering on the fibre) | joins `(★dec)`/H1 | **[O]** |
| `G''_P > 0` at all 72 scanned `(beta,theta,v,s)`; no counterexample to (★) | witness | **[V]** |
| acid test: free Gaussian (`Phi=0`) gives `G''_P = 0`, (★) fails | witness `1e-15` | **[T]** |
| Reading 2 (marginal `P_*mu`) extra term (law of total variance / Fisher) | out of scope; separate door | **[O]** |
| **(★) itself** | UV tower safe up to `O(g_j)`; head unchanged (now carries `R_proj`) | **[C]/[O]** |

---

## 10. Verification ledger & witness

| item | check | status |
|---|---|---|
| case (a): `G''=2W_geo`, `R_proj=0` | `g3_projection.py` §0, residual `<=2e-7` | **[V]** |
| case (b): `R_num` vs analytic `R_proj` | §1, all 72 points `<=1.1e-7` | **[V]** |
| `R_proj` sign changes (`-0.114` at `theta=pi`; `+0.255` at `theta=pi/2`) | §1, §3 | **[V]** |
| `R_proj=0` at `theta=0` (isotropy) | §1 | **[V]** |
| `R_proj -> 0` as `s -> 0` (recovers case a) | §2b, `O(s^2)` for axial transport | **[V]** |
| free Gaussian `Phi=0`: `G''_P = 0` | §2, `1e-15` | **[T]**+[V] |
| `div w = 0` but `R_proj != 0` (v=e1) ⇒ `(Dw)` terms alone generate it | §2b (`E[div w]=0`, `R != 0`) | **[V]** |
| conventions match S2/S5 (`G` from `S`, `W` from `Phi`); `beta=2.5,theta=pi,e1` gives `W_geo=0.0659` | cross-checked against S2 table §3 | **[V]** |

Witness `research/v3/w7/g3_projection.py`: SU(2) one plaquette, one free link, boundary
`B = exp(i theta sigma_3/2)`; deterministic 3-D Gauss-Legendre quadrature in `(rho, u=cos theta,
psi)` (`90 x 56 x 56`), `Phi` radial with cancellation-free `Phi_r, Phi_rr` (S2 forms); `S`
gradient/Hessian by central differences (`h=1e-4`); `w(X) = Ad_{exp(-sX)}v` by Rodrigues, `Dw` by
finite differences. Run `py -3.12 g3_projection.py`. Output quoted in §2, §4, §5, §7.

---

## 11. Sharpest line (restated)

**A constant-coefficient projection leaves `G'' = 2W_geo` exact — `w=Pv` is just another fixed fine
direction and S2 is directional; a covariant (parallel-transport) projection adds exactly one term,
`R_proj = E[<(Dw)w,grad Phi> - <(Dw)w,grad S>] + Cov(<w,grad Phi>,div w) - Cov(<w,grad S>,div w)`,
which is the covariance form of `int div(Pw) e^{-S}J` and the failure of S2's Lemma 2. It is
sign-indefinite and `O(W_geo)` at the IR head (up to `76%` of the floor in the toy, never crossing
zero there), but `O(g_j)`-small and hence harmless in the UV tower.** The projection does not break
`(★)_j` in the UV; it adds a named head term to `(★dec)`/H1. Free field still gives `G''=0`, so the
acid test passes. Clay 0%.

---

*Refs.* `S2-BOUNDARY-COUNTEREXAMPLE.md` §2 (the directional identity, Lemma 1–Cor.),
§6 (interior-constraint caveat); `S3-UV-FESHBACH.md` §1–2 (Schur split `G''=E[Σ]-D`, defect `D`);
`S5-N-UNIFORM.md` §2 (N/a table; boundary-plaquette `O(1)`-per-mode counting), §3 (Haar reservoir);
`docs/paper-v2/ch05-theorems.tex` eqs. (eq:wgeo, eq:split, eq:vardecomp). Witness
`research/v3/w7/g3_projection.py`. No third-party sources downloaded; all external inputs are
standard (Stokes/total-derivative, law of total variance) and cited from the corpus.

# STAR2-CLASSIFICATION — (★) is classified: strong form false [NG], good-set theorem [T]/[C], head = H1 [O]; the 2W_geo route is dead

> **VERDICT (no QED, Clay 0%).** The `(★)` wave is now decided as a **classification**, not a
> proof of the gap:
>
> 1. **The strong form (uniform in the coarse field `a₀`) is FALSE `[NG]`** — a rigorous
>    counterexample in an exactly solvable SU(2) model: `G″=0.20511530−0.01105929 β²`, negative
>    for `β>4.30661`, with the longitudinal mode in `im P` and full-support boundary.
> 2. **On the small-field set `𝒮_t` the true `G″` has a genuine positive bound** `[T]/[C]`:
>    `G″_true ≥ c_□+N/12−C₁βF_∂ > 0`, linear in `N`, via marginal strong convexity
>    (`Hess V ⪰ ρI ⇒ G″≥ρ`), with `X` included.
> 3. **The head (typical field, beyond `𝒮_t`) is H1 `[O]`** — the bulk Savvidy modulus is `≤0`,
>    the true `G″` is measured `+0.13` at `β=2.5`, but proving it is the gap.
> 4. **The `2W_geo` reduction is dead `[NG]`**: G1/G3 computed the *unconditioned, X-less*
>    object; the true `(★)` is the **conditioned coarse curvature** `G″(v;a₀)=E_a[∂_v²V]−Var_a(∂_vV)`.
> 5. **The admissible projection class is linear/divergence-free `[T]`** (the acid test selects
>    it; covariant directions give false positivity or false negativity).
> 6. **The UV statement is renormalized** `(★_ren)`: `Z_j^{-1}G″ ≥ c_L`, with
>    `δ_W=Θ(β_j) ∝ log(1/a_j)`; a constant bare floor is false.
>
> This is a rigorous no-go plus a theorem plus a reduction; it does not prove the mass gap.

---

## §0 · The correct object

`(★)` is about the **coarse-field-conditioned marginal curvature**:
\[
G″_j(v;a_0)=E_a[\partial_v^2V]-Var_a(\partial_vV)=\Delta_{\rm dyn}+W_j=E[\Sigma]-D,
\qquad W_j=W_{\rm geo}+X,\quad X=-2\operatorname{Cov}(G,\Phi'),
\]
conditioned at `a=a₀` (ch05 `eq:split`/`eq:vardecomp`). It **includes `X`** and **must be
conditioned**. The object computed by S2/G1/G2/G3,
`G″_0=E_\mu[H]-Var_\mu(G)+W_{\rm geo}=2W_{\rm geo}`, is a true identity but **not `(★)`**.

---

## §1 · The strong form is false `[NG]`

Exactly solvable model: SU(2), one plaquette, boundary holonomy `B=e^{iθσ₃/2}`, one free link,
coarse field `a₀`. At `θ=π`, longitudinal `v=e₃`, `a₀=0` (the barrier):
\[
G″(v;a_0)=A-B\beta^2,\qquad A=0.20511530,\quad B=0.01105929,\quad \beta_*=4.30661 .
\]
Three independent methods agree to `≤10^{-6}` (this corrects G4's `0.2057/0.0123/4.09`, which
had a missing quadrature weight). Negative for `β>β_*`. Mechanically, conditioning sets
`E[H|a₀]=0`, exposing `−β²Var(G)`. Transverse modes and `θ=0` stay positive.

**Scope:** longitudinal modes lie in `im P`; `a₀=0` and `B→−I` are full-support. So the
counterexample is **inside the definition** of the true `(★)`. The strong "uniform in every `a₀`"
form is **`[NG]`**. The literal C2 variance inequality `2W_geo` is a different object and is not
refuted.

---

## §2 · The good-set theorem `[T]/[C]`

On the small-field set `𝒮_t` the true object is controlled exactly:
\[
Hess V\succeq \rho I\quad\Longrightarrow\quad G″_{\rm true}(v;a_0)\ \ge\ \rho(B,U′)
=\operatorname{ess\,inf}\lambda_{\min}(Hess V),
\]
with `X` included (marginal strong convexity). The modulus is
\[
\rho = N/12 + c_\Box - C_1\beta F_\partial - C_S\beta F_{\rm bulk}.
\]
On `𝒮_t` (`F_∂≤t`, small block), `G″_{\rm true} ≥ c_\Box+N/12-C_1\beta F_\partial > 0`, linear in
`N` `[C]+[V]`. One-plaquette at `β=2.5`: threshold `θ_c≈0.81`, `ρ_0(0)=+0.217`, true
`G″(e₃)=0.741≥0.217`. The Schur loss can exceed `N/12`; the corrected bound is the exact Schur
complement `Σ⪰ρ` (no loss).

---

## §3 · The head is H1 `[O]`

Beyond `𝒮_t` (general block, typical field), `ρ≤0` (bulk Savvidy). The true `G″` stays `+0.13`
at `β=2.5` (measure-carried defect `D`), but proving it is **H1 `[O]`**. This is the same wall as
before, now reached from the UV side with the correct object.

---

## §4 · Projection and normalization

- **Projection class `[T]`:** the acid test selects **linear/divergence-free** directions. For
  these, the X-less identity `G″_0=2W_geo` is exact for any block (S2/G1/G2). Covariant
  (right-invariant) directions give `R=G″-2W_geo≠0`, sign-indefinite, with false positivity or
  false negativity; they must be excluded.
- **Normalization `[T]/[C]`:** a constant bare `δ_W` is false at the UV end; the correct form is
  `(★_ren): Z_j^{-1}G″_j ≥ c_L=4π²/L²`, with `δ_W=Θ(β_j) ∝ log(1/a_j)`. The UV block must pass
  the acid test (`Q_c=4π²/L²>0`); the IR zero mode fails (`Q(0)=0`).

---

## §5 · The pincer, honestly

- **Strong jaw** (S6 Dobrushin `β<1/72`; SZZ `β<1/48` ⇒ `β_W≈0.19`): `[T]`, but **grips zero
  tower shells** (thresholds crossed at `j<1`; the tower sits at `β_j≥2.5`).
- **UV jaw** (BL + Lemma B / S1): covers `j≥j_UV` (`j_UV=1–3`), all but `O(1)` shells, but only
  `[C]` on a moment bootstrap and `(★dec)`, and only on `(★_sf)`.
- **Head**: `O(1)` shells, `β≈2.5`, `[O]=H1`.

No whole intermediate regime inside the tower; the head is the gap.

---

## §6 · What was won

1. **A rigorous no-go:** the uniform-in-`a₀` UV inequality is false; the `2W_geo` certificate
   mis-certifies exactly where it fails. This kills the "L3 via W_geo" route.
2. **A genuine theorem on a domain:** the true `G″` is bounded below by
   `c_□+N/12−C₁βF_∂>0`, linearly in `N`, on the small-field set, with `X` included.
3. **A precise reduction:** the head is H1, and the UV face is not a new route.
4. **Two selection principles:** the projection class is linear/divergence-free; the UV floor is
   `log(1/a)`-running after `Z_j`.

This is the victory available at this stage: **a classification with a theorem and a no-go, not a
proof of Clay.**

---

## §7 · Honest bottom line

`(★)` is **not** a single open inequality anymore. It is:
- **false** as a uniform statement in the coarse field `[NG]`;
- **true** on the small-field set, with an explicit positive modulus `[T]/[C]`;
- **H1** on the head `[O]`;
- correctly **renormalized** at the UV end `[T]/[C]`;
- restricted to the **linear/divergence-free** projection class `[T]`.

**Clay = 0%. No QED.** The wall did not move; the UV side of it is now fully mapped: one door is
locked (`[NG]`), one door opens onto a proved room (`𝒮_t`), and the last door opens onto the same
room as before (H1).

*Refs:* `RECONCILIATION-STAR.md`, `V1-UV-NORMALIZATION.md`, `V2-TRUE-BOUND.md`,
`V3-COUNTEREXAMPLE-SCOPE.md`, `V4-PINCER.md`, `G2-TWO-PLAQUETTE.md`; `S1`–`S6`; witnesses in `w7/`.

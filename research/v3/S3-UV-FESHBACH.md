# S3-UV-FESHBACH — the two-block/Feshbach decomposition of the UV shell inequality (★): what it derives, what W_j must contain, and why the IR TCI does not imply it

> 2026-11-24 · Attack surface S3: derive the live inequality (★) — the UV one-sided bound on the
> conditioned force variance, equivalently `G″_j ≥ δ_W` — from a Schur/Feshbach two-block split of
> the shell Hessian, and adjudicate C2's claim that the IR two-block lemma (TCI) is "same form,
> different object."
> Reads: `research/v3/C2-BALABAN-UNIFORM.md` §2.4 and §4; `docs/paper-v2/ch05-theorems.tex`
> (Feshbach–Bogoliubov lemma, Lemma NS, variance decomposition); `research/unconditional/FORGE-BLIND.md`
> §2–4; `D-BBS-INTEGRATE-NEG.md`; `INVENT-VARIANCE-INEQ.md`; `PROOF-JACOBIAN-SAVES.md`;
> `G1-L3-UV-FORCE-BOUND.md`; `PROOF-L3.md`.
> Witness: `research/v3/w7/s3_feshbach_witness.py` (deterministic; output quoted in §5).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step / open · **[R]** reported / unrefereed.

---

## 0. Verdict (no QED)

**(★) admits an exact two-block/Schur decomposition, and in the UV it is *equivalent* to a
Schur-complement curvature lower bound — given a spectral gap on the fibre.** That is the only
structural statement I can prove; it does **not** close (★). Precisely:

1. **Exact identity [T].** With `V = S + Φ` (`Φ = −log J` the Haar potential), the coarse/test mode
   `v`, and the conditioning on the coarse field `U′` modelled as `a = ⟨PᵀX,v⟩ = a₀`, the shell
   Hessian splits into a coarse block `A`, a **fibre** block `C = K` (the directions decimated away),
   and a coupling `c`. Its **Schur complement** `Σ = A − ⟨K⁻¹c,c⟩` is the Gaussian (one-loop) coarse
   curvature, and
   ```
   G″_j(v)  =  E_{μ(·|U′)}[Σ]  −  D ,      D := Var_{μ(·|U′)}(f) − E_{μ(·|U′)}[⟨K⁻¹c,c⟩] ,
   ```
   with `f = ⟨Pᵀ∇V,v⟩` the force. This is `INVENT-VARIANCE-INEQ` §1 in two-block language, and equals
   `G″ = Δ_dyn + W_geo + X` of `PROOF-JACOBIAN-SAVES` after substitution (§1).
2. **UV reduction [T corollary].** If the fibre has a gap `K ⪰ γ_j I > 0`, then `μ(·|U′)` is
   log-concave on the fibre (Prékopa), Brascamp–Lieb gives `D ≤ 0`, hence `G″_j ≥ E[Σ]`. Therefore
   ```
   (★)  ⟸  [ fibre gap  K ⪰ γ_j I ]  ∧  [ Schur bound  E_{μ(·|U′)}[Σ] ≥ δ_W ] .
   ```
   The variance has been traded for the **expected Schur complement** — a curvature (operator) object,
   not a fluctuation. This is the value of the decomposition: `(★)`'s hard content is *not* "bound a
   variance" but "the expected Gaussian coarse curvature stays positive" — where the positivity engines
   (Prékopa, BL, and at small `g_j` the LSI) can dock.
3. **`W_j` is exactly the Jacobian (Haar) contribution.** In the S-split, `W_j = W_geo + X`,
   `W_geo = E[Φ″] − Var(Φ′)`, `X = −2Cov(G,Φ′)`, `G = ⟨Pᵀ∇S,v⟩`. Its **provable content [T]** is the
   pointwise floor `Φ″ ⪰ (N/12)I` (`PROOF-JACOBIAN-SAVES` Lemma B), which enters `A` (and `K`)
   positively. Its job: **raise the coarse-curvature floor above the Schur loss `⟨K⁻¹c,c⟩` so that
   `E[Σ] ≥ δ_W` can hold.** It **cannot** by itself absorb the Schur loss or the cross term `X`
   (§2.3). So `W_j` "absorbs the geometric remainder" only in the narrow sense that `Φ″ ⪰ N/12` is a
   proven positive floor; it does **not** absorb the dynamical remainder.
4. **C2's "same form, different object": CONFIRMED and sharpened.** Both `(★)` and TCI are
   Schur-complement positivity of a *small* block against a *gapped* complement — the *same algebra*.
   But they are **different objects**: TCI's retained block is the *tachyonic* Savvidy block
   (`A ≺ 0`; TCI's truth *is* the gap), `(★)`'s is the *coarse UV mode* already `A ⪰ 0` in the UV
   (its truth is *certification*); different operator (`−D²+2adF` at `β≈2.5` vs the conditional shell
   marginal at `g_j→0`); TCI has no conditioning, no measure, no defect `D`, no Jacobian `W_j`.
   **TCI does not imply (★).** The sharpest form of the mismatch: TCI and `(★)` eliminate a
   **positive** complement (`C ⪰ γ`, so the Schur correction `−BC⁻¹Bᵀ ⪯ 0` enters with the sign that
   *hurts* the small block), whereas the corpus's other two-block tool (Lemma NS / `D-BBS`) eliminates
   a **negative** block (`C ≺ 0`, Schur correction `−BC⁻¹Bᵀ = +B|C|⁻¹Bᵀ ⪰ 0` *helps*). The two
   conventions are dual and must not be conflated; citing "the two-block lemma" as if there were one
   is the exact trap (§2.4, §3).
5. **Verdict [C] on `[fibre gap] ∧ [Schur-loss bound]`, [O] to close.** Not `[T]`, not `[NG]`.
   Clay 0%.

---

## 1. Setup and the exact Schur–variance identity  **[T]**

Fix shell `j` and coarse configuration `U′`. The conditional (fibre) measure is
`μ(·|U′) ∝ e^{−V_j(X;U′)} dX`, `V_j = S_j + Φ`, `Φ = −log J` the Haar potential (flat chart
`X_ℓ = log U_ℓ`). The coarse field `U′` enters the physical collective mode `v` through the gauge
projection `P` (the orbit-transverse projection `Q⊥`); along `v` we condition on
`a := ⟨PᵀX, v⟩ = a₀`, with `a₀` the value encoded by `U′` (this is the local model of "conditioning on
`U′`"; `INVENT-VARIANCE-INEQ` §1, `PROOF-JACOBIAN-SAVES` §1). On the slice `{a = a₀}`, with
`Π₀ = I − vvᵀ` the fibre projector, define the **coarse (A) block**, the fibre (**C**) block, and the
coupling:
```
A(X) = ⟨v, Pᵀ Hess V_j(X) P v⟩ ,        K(X) = Π₀ Pᵀ Hess V_j(X) P Π₀ ,     c(X) = Π₀ Pᵀ Hess V_j(X) P v ,
Σ(X) = A − ⟨K⁻¹ c, c⟩   (Schur complement) ,       f(X) = ⟨Pᵀ ∇V_j(X), v⟩   (force) .
```
This is exactly the `[[A, c],[cᵀ, K]]` block form of `FORGE-BLIND` §2 / ch05, with the *small* block
`A` the coarse direction and the *complement* `K` the fibre.

**Lemma S (Schur–variance identity).**
```
G″_j(v) = E_{μ(·|U′)}[A] − Var_{μ(·|U′)}(f) = E_{μ(·|U′)}[Σ] − D ,
D := Var_{μ(·|U′)}(f) − E_{μ(·|U′)}[⟨K⁻¹c,c⟩]  (Brascamp–Lieb defect of the conditioned fibre).
```
*Proof.* The marginal-curvature identity (`INVENT-VARIANCE-INEQ` (1); `PROOF-JACOBIAN-SAVES` (1)):
differentiating `−log ρ(a)` twice at `a=a₀`, with `f = ∂_n V`, gives `G″ = E[A] − Var(f)`. Add and
subtract `E[⟨K⁻¹c,c⟩]`: `E[A] − Var(f) = (E[A] − E[⟨K⁻¹c,c⟩]) − (Var(f) − E[⟨K⁻¹c,c⟩]) = E[Σ] − D`. ∎

**Reconciliation with the campaign's split.** Substituting `A = H_vv + Φ″_n`, `f = G + Φ′`,
`Σ = (H_vv + Φ″_n) − ⟨K⁻¹c,c⟩` reproduces `PROOF-JACOBIAN-SAVES`:
`G″ = Δ_dyn + W_geo + X`, `Δ_dyn = E[H] − Var(G)`, `W_geo = E[Φ″] − Var(Φ′)`, `X = −2Cov(G,Φ′)`; the
two identities are the same statement, `Σ`-form vs `W`-form. The witness (§5(1),(3)) verifies
`G″ = Σ` exactly for the Gaussian and `G″ = E[Σ] − D` to quadrature for the non-Gaussian model.

---

## 2. The UV reduction and the role of `W_j`

### 2.1 `(★)` from a fibre gap plus a Schur bound  **[T]**

**Corollary S.** Suppose `K(X) ⪰ γ_j I` on the typical slice (fibre gap — the spectral-gap assumption
on the complement). Then the conditioned fibre measure is log-concave (Prékopa: conditioning a
log-concave density on the linear constraint `a = a₀` stays log-concave), Brascamp–Lieb (1976) gives
`Var_{μ₀}(f) ≤ E_{μ₀}[⟨K⁻¹c,c⟩]`, i.e. `D ≤ 0`, so
```
G″_j(v) = E[Σ] − D  ≥  E[Σ] .
```
Hence `(★)` (i.e. `G″_j ≥ δ_W`) follows from `E_{μ(·|U′)}[Σ] ≥ δ_W`. ∎

The witness (§5(3)) checks this on a non-Gaussian 2D model: `K > 0` everywhere gives `D = −0.0208 < 0`
and `G″ = 0.260 ≥ E[Σ] = 0.240`; the concave-fibre contrast (`K` changes sign) gives `D = +7.29 > 0`,
so the sign of `D` is exactly the sign of the fibre's log-concavity. This is the precise sense in which
`(★)` **is** a statement about the fibre gap, not an unconstrained variance.

### 2.2 What the Schur bound says, and where `W_j` sits  **[T]/[C]**

Write `A = A_S + Φ″_n` and `K = K_S + K_Φ`, so `Σ = A_S + Φ″_n − ⟨(K_S+K_Φ)⁻¹ c, c⟩`. Using the
**proven pointwise floor** `Φ″ ⪰ c_Φ I`, `c_Φ = N/12` (SU2 `1/6`, SU3 `1/4`; `PROOF-JACOBIAN-SAVES`
§2, `ch05` Lemma B — under *any* measure, any `β`),
```
E[Σ] = E[A_S] + E[Φ″_n] − E[⟨K⁻¹ c,c⟩]  ≥  E[A_S] + c_Φ − E[⟨K⁻¹ c,c⟩] .
```
In the deep UV the physical coarse mode sits above the free box gap (`A_S ⪰ c_L = 4π²/L²`, Route B
cone-emptiness), the coupling graph angle is perturbatively small (`c` is interaction-induced; `v` is
a free eigenmode and the free Hessian is diagonal in the free basis), so the **Schur loss**
`⟨K⁻¹c,c⟩ = (graph angle)² · A` is `o(A)` (the squared angle is small; the exact power is unproven —
the campaign's `Σ_j g_j⁴ < ∞` counting, `G1-UNIFORM-ZA` §2, is the relevant series). Hence in the
deep UV
```
E[Σ]  ≳  (1 − o(1)) (c_L + c_Φ) ,        so  (★) holds with δ_W ≳ c_L + c_Φ − o(1) .
```
This is a **conditional** quantitative statement: it presupposes (i) the fibre gap `K ⪰ γ_j`
(§2.1), and (ii) that `c` is a perturbative correction, i.e. the graph angle is small. Both are the
open hard analysis of Route B (`G1-L3-UV-FORCE-BOUND` §3, §6):
```
(H-C)  fibre gap  K ⪰ γ_j I  (UV log-concavity; = cone-emptying + moment bootstrap E|F_j| ≤ C g_j)   [O]
(H-SL) Schur loss  E[⟨K⁻¹c,c⟩] ≤ c_L + c_Φ − δ_W  (graph-angle / equipartition-defect bound)            [O]
```
Given `(H-C) ∧ (H-SL)`, `(★)` follows. `A_S ⪰ c_L` is *not* an extra hypothesis in the deep UV — it is
the same cone-emptiness as `(H-C)` — but at the **head** (`j` small, `g_j = O(1)`, `β ≈ 2.5`) the
action block `A_S` can be Savvidy-indefinite, `K` can fail `⪰ 0`, and neither holds; that end needs the
cross-block clustering `(★dec)`.

### 2.3 What `W_j` must contain — and what it cannot absorb  **[T]/[O]**

In the S-split, `W_j := W_geo + X = E[Φ″] − Var(Φ′) − 2Cov(G,Φ′)`. Two facts:
- **Proven [T]:** `W_geo ≥ 0` (Prékopa on the pure-Haar marginal), with the pointwise floor
  `Φ″ ⪰ c_Φ = N/12`. This is what raises the tolerable Schur loss from `c_L − δ_W` (action-only) to
  `c_L + c_Φ − δ_W`: **`W_j` must contain at least `c_Φ = N/12` at the coarse block.**
- **Not sign-controlled [O]:** `X = −2Cov(G,Φ′)` is indefinite. `W_j` **cannot** absorb `X` without a
  bound on the force–Jacobian covariance, and the Schur loss `⟨K⁻¹c,c⟩` is not bounded by `W_geo`.
  This is the exact sense in which "the Jacobian carries the sign" (`PROOF-JACOBIAN-SAVES` §3) is a
  *measured* statement at `β=2.5` (`W_full ≈ +0.28` beats `−Δ_dyn ≈ +0.02`) and a *proven* statement
  only through the `Φ″ ≥ N/12` floor, not through `X`.

In the V-picture there is no separate `W_j`: `W_geo` and `X` are dissolved into `Σ` (via `Φ″`,
`Π₀ Hess Φ`) and `D` (via `Var(Φ′)`, `Cov`). The witness §5(5) shows the mechanism concretely: an
action-only model with strong coupling gives `E[Σ_S] = 0.077 < δ_W = 0.10`, while adding the Jacobian
floor `Φ″ = N/12 = 0.167` raises `E[Σ_V] = 0.399 > δ_W`. **`W_j`'s provable content is the `Φ″`
floor; its job is to pay the Schur loss `⟨K⁻¹c,c⟩`, and that payment is exactly sufficient only if
`(H-SL)` is certified.**

### 2.4 Two conventions — the sign of the Schur correction is the whole story

The corpus contains **two** inequivalent two-block manipulations, and the sign of the Schur term is
opposite:
```
(A) TCI / (★) convention  — eliminate the GAPPED complement  C ⪰ γ > 0 :
      Σ = A − B C⁻¹ Bᵀ = A − (positive)  ⪯  A      (coupling HURTS the small block).
    FORGE-BLIND §2 / ch05: "Σ ≥ μ" on the tachyonic small block A.

(B) Lemma NS / D-BBS convention — eliminate the NEGATIVE Savvidy block  C ≺ 0 :
      Σ = A − B C⁻¹ Bᵀ = A + B|C|⁻¹ Bᵀ  ⪰  A      (coupling HELPS the survivor).
    D-BBS-INTEGRATE-NEG (b); FORGE-FESHBACH-ANGLE §1; witness §5(2).
```
`(★)` is convention **(A)**: the eliminated object is the *fibre* (decimated directions), which in the
UV is positive (`K ⪰ γ_j`), and decimation therefore *softens* the coarse curvature (`Σ ⪯ A`). TCI is
also convention **(A)** (retained block = tachyonic, eliminated block = gapped Savvidy complement). The
`D-BBS`/`FORGE-FESHBACH` "eliminate the negatives" tool is convention **(B)**. **Conflating (A) and (B)
is exactly the category error** that makes "the two-block lemma" seem to close both the IR and the UV;
(B)'s "coupling only helps" does not apply to `(★)`, where the coupling hurts.

**The two-block lemma *per se* is not the tool.** The `FORGE-BLIND` bound
`A ≥ min{μ/(1+2‖B‖²/γ²), γ/2}` is a *full-Hessian spectral* statement; its leakage factor
`1/(1+2‖B‖²/γ²)` is the price of converting a Schur bound into a spectrum. `(★)` needs only the
*marginal* curvature, for which the Schur bound `E[Σ] ≥ δ_W` is direct and **leakage-free** (Jensen
is not even needed: `G″ = E[Σ] − D ≥ E[Σ]` when `D ≤ 0`). So the honest tool for `(★)` is Lemma S +
Brascamp–Lieb, and the two-block lemma enters only as the abstract reason `Σ ≥ μ ⇒ gap` — which the UV
already has from cone-emptiness and does not need.

---

## 3. Why the IR TCI does not imply the UV `(★)`

C2 §4 says "same form, different object." The decomposition sharpens *both* halves.

**Same form [T].** Both are convention-(A) Schur positivity: retain a small block, eliminate a gapped
complement `C ⪰ γ`, demand `Σ = A − BC⁻¹Bᵀ ≥ μ`.

**Different object [T] — five independent mismatches.**

| | UV `(★)` | IR TCI (`FORGE-BLIND` §2–3) |
|---|---|---|
| **operator** | conditional shell marginal of `V_j = S_j+Φ`, `g_j→0` | frozen quadratic `−D²+2adF` at `β≈2.5`, `g≈1.5` |
| **retained small block `A`** | coarse UV mode, **already `A_V ⪰ c_L + N/12 > 0`** (cone empty) | **tachyonic** Savvidy block, `A ≺ 0` |
| **what `Σ ≥ μ` means** | *certification* that decimation preserves an existing positivity | **the gap itself** (positivity manufactured by interaction + coupling) |
| **complement `C`** | UV fibre, gap `γ_j ∼ Λ_j²` (free Laplacian scale) | Savvidy gapped sector, `γ_0 ∼ F` (background-field scale) |
| **conditioning / defect** | conditioned on `U′`: introduces `a=a₀`, the fibre, and the Brascamp–Lieb defect `D` | *unconditional* quadratic form; no measure, no variance, no `D` |
| **Jacobian `W_j`** | present: `Φ = −log J`, floor `Φ″ ⪰ N/12` | **absent** (Hessian statement, no `Φ`, no Haar) |
| **boundary / decimation** | per-block, `U′`-boundary ⇒ boundary covariance `(★dec)`; Schur = **decimation** (coarse-graining map) | fixed torus, no boundary; Schur = **spectral** split (invariant subspace), not a decimation |
| **sign of eliminated block** | `C ⪰ γ > 0` ⇒ Schur term hurts (§2.4-A) | TCI also convention (A); but `D-BBS` (B) eliminates `C ≺ 0` and helps |

**Why the implication genuinely fails.** A proof of TCI would establish `Σ_IR ≥ μ > 0` for the
**Savvidy operator at `β≈2.5`**, a spectral statement about the *full* Hessian. It says nothing about:
(i) the **conditional marginal** `G″_j` of a *different* operator (the shell `V_j` at `g_j→0`), because
`G″` is `E[Σ] − D`, an expectation over the conditioned fibre, with a defect `D` TCI never sees;
(ii) the **decimation** Schur elimination (coarse-graining), as opposed to TCI's spectral split;
(iii) the **Jacobian** `W_j`, which has no TCI analogue and whose floor `Φ″ ⪰ N/12` is a component of
`(★)`'s positive curvature; (iv) the **`U′`-boundary covariance** `(★dec)`, absent from the
no-boundary TCI. Conversely, `(★)` at the UV end does not give TCI (C2 §4): the regimes are disjoint
(`g_j→0` vs `g≈1.5`) and the retained block's *sign* is opposite. The "one species" conjecture, if it
is to be cashed, must be a **single** object that degrades gracefully from certification (`A ⪰ 0`,
UV) to sign control (`A ≺ 0`, IR); the decomposition shows the two are the *same algebra on dual small
blocks*, which is the sharpest available statement of the conjecture — **and the reason it is only a
conjecture, not a proof: the duality flips the sign of the very term that does the work.**

---

## 4. Status: the conditional derivation and the failing step

**Succeeds conditionally.** `(★)` follows from:
```
(H-C)  fibre gap       K(X) ⪰ γ_j I  on the typical slice          [O]  (UV log-concavity / Route B (i))
(H-SL) Schur-loss bound E[⟨K⁻¹ c,c⟩] ≤ E[A_S] + c_Φ − δ_W   [O]  (graph angle small; Route B (ii))
```
Given `(H-C)` the defect is signed (`D ≤ 0`); given `(H-SL)` the expected Schur complement
`E[Σ] ≥ E[A_S] + c_Φ − E[loss] ≥ δ_W` clears
`δ_W`, with `W_j`'s `Φ″ ⪰ N/12` floor doing the positive work. Both are the campaign's existing Route
B items (`G1-L3-UV-FORCE-BOUND` §6), restated as *curvature* conditions rather than a variance bound.

**Failing step (exactly named).** Neither `(H-C)` nor `(H-SL)` is proven.
- `(H-C)` requires the gauge-invariant `a`-uniform moment/cone bootstrap `E|F_j| ≤ C g_j` — the same
  circular-with-the-flow object the campaign already flags (`G1-L3` §4(a), `G1-L3-UV-FORCE-BOUND` §6(1)).
  Without it the fibre can be Savvidy-indefinite and `D > 0` (witness §5(3) contrast), so `G″ ≥ E[Σ]`
  fails.
- `(H-SL)` requires bounding the graph angle `c` (equivalently the equipartition defect `Var − E[H]`),
  blocked by the `a⁻⁴` renormalon (`CONTINUUM-ASYMPTOTIC-BEYOND`); at the **head** (`β≈2.5`) `c_L` is
  replaced by a finite head constant and clustering `(★dec)` is needed.

**What changes vs the campaign's framing.** The decomposition shows the binding inequality is **not**
"`Var ≤ E[H]`" (one-loop equipartition, which gives `G″ ≈ 0`, not `≥ δ`) but
"`E[A] − E[⟨K⁻¹c,c⟩] ≥ δ_W`" — the **Jacobian floor plus free gap minus the Schur loss**. This is why
`W_j` is *necessary*: at `β=0` (`V = Φ` convex) `A_S = 0`, `K ⪰ c_Φ`, `D ≤ 0`, and Lemma S gives
`G″ ≥ E[Σ] ≥ c_Φ = N/12` **unconditionally [T]** (Schur-complement monotonicity: `Hess Φ ⪰ c_Φ I ⇒
Σ ⪰ c_Φ`) — the `β=0` theorem of `ch05`/`PROOF-JACOBIAN-SAVES`
is exactly the `Σ`-form at `A_S = 0`. At `β≈2.5` the action part turns indefinite and the fibre gap can
fail; that is the whole difference between a theorem (`β=0`) and the open inequality (`β≈2.5`, UV tower).

**Acid test.** For the *identity* the correct free-field check is that it manufactures **no** positivity:
`V = ½⟨X,AX⟩` massless gives `S = D = 0`, `G″ = Σ = 0` (witness §5(4): `G″ = E[Σ] = D = 0`). This
**passes** — as it must: for a *UV-certification* statement (`C2 §5`), the free field is UV-stable, so
`(★)`-by-Schur correctly returns `0`, the opposite of the gap statement H1 (which must fail on the free
field). Do not apply the H1 acid test (`f₁² = 0 ⇒ vacuous`) to `(★)`; the correct `(★)` test is
non-vacuity with explicit `δ_W`, which the decomposition supplies (`δ_W ≳ c_Φ = N/12` in the deep UV,
conditional on `(H-C)`).

---

## 5. Witness — `research/v3/w7/s3_feshbach_witness.py`

Deterministic, `py -3.12`, quoted output:

```
(1) GAUSSIAN, C > 0:  || G'' - Schur Sigma || = 2.2e-16 ; eig(Sigma - A) = [-0.413, -0.050]  => Sigma <= A
(2) SIGN FLIP, C < 0:  eig(Sigma - A) = [+0.019, +0.140]                                     => Sigma >= A   (Lemma NS)
(3) NON-GAUSSIAN 2D:  K_min = 1.05>0 ; G''=0.2605 ; E[Sigma]=0.2397 ; D=-0.0208 ; identity error 0.0 ; G''>=E[Sigma]
    concave fibre:      K_min=-120  ; D=+7.286  (D>0 when fibre not log-concave)
(4) ACID TEST (massless free):  G'' = E[Sigma] = D = 0   (no false positivity)
(5) W_j ROLE:  action only E[Sigma_S]=0.077 < delta=0.10 ; + Jacobian Phi''=N/12: E[Sigma_V]=0.399 > delta
```
Readings: (1) Gaussian coarse curvature **is** the Schur complement, and eliminating a *positive* fibre
*hurts*; (2) eliminating a *negative* block *helps* (Lemma NS) — the dual convention; (3) the exact
identity `G″ = E[Σ] − D` holds, `D ≤ 0` iff the fibre is log-concave, so `G″ ≥ E[Σ]`; (4) the free field
returns `0`; (5) the Jacobian `Φ″ = N/12` floor is the proven content of `W_j` that lifts `E[Σ]` over
`δ_W`. **Not a proof — a witness for the algebra and the sign structure.**

---

## 6. Status ledger

| item | statement | label |
|---|---|---|
| Schur–variance identity `G″ = E[Σ] − D` | exact, two-block form of `PROVE-Gpp` | **[T]** |
| fibre gap ⇒ `D ≤ 0` ⇒ `G″ ≥ E[Σ]` | Prékopa + Brascamp–Lieb | **[T]** |
| `(★) ⟸ (H-C) ∧ (H-SL)` | conditional reduction to a Schur (curvature) bound | **[C]** |
| `W_j` provable content = `Φ″ ⪰ N/12` floor | Lemma B; pays the Schur loss | **[T]** |
| `W_j` absorbs the cross term `X` / Schur loss | **false** without `(H-SL)` | **[O]** |
| `β=0`: `G″ ≥ E[Σ] ≥ N/12` | `A_S=0`, fibre convex — the `Σ`-form of the `β=0` theorem | **[T]** |
| `(H-C)` fibre gap `K ⪰ γ_j` (UV log-concavity/moment bootstrap) | Route B (i) | **[O]** |
| `(H-SL)` Schur-loss bound (graph angle / equipartition defect) | Route B (ii) | **[O]** |
| head `β≈2.5`: `A_S` indefinite, `K` may fail `⪰0`, needs `(★dec)` | clustering | **[O]** |
| TCI ⇒ `(★)` | **does not hold**: different operator/regime/sign/conditioning/no-`W` | **[T]** (structural) |
| TCI "same form, different object" (C2) | **confirmed, sharpened**: same (A)-algebra, dual small blocks | **[T]** |
| mixing convention (A) `C⪰γ` with (B) `C≺0` | category error | **[NG]** for the conflation |
| **(★) itself** | reduced to `(H-C) ∧ (H-SL)`; still open | **[C]/[O]** |

---

## 7. Sharpest line

**(★)`s two-block decomposition is exact, and in the UV it converts the one-sided *variance* bound into
a one-sided *Schur-complement* (curvature) bound — `G″ = E[Σ] − D`, `D ≤ 0` given a fibre gap — with
`W_j` identified as the Haar/Jacobian contribution whose only provable content is the pointwise floor
`Φ″ ⪰ N/12`.** That floor is *necessary*: it is what lifts the expected coarse curvature over the
decimation's Schur loss. The derivation **succeeds conditionally** on the fibre gap `(H-C)` and the
Schur-loss bound `(H-SL)` — the campaign's existing Route B items — and **fails exactly there**: the
moment bootstrap and the graph-angle/equipartition-defect bound are open, and at the head `β≈2.5` the
fibre can go Savvidy-indefinite so `D > 0`. C2's "same form, different object" is **right and sharper
than stated**: TCI and `(★)` are the same Schur algebra on *dual* small blocks, but the retained block's
sign flips (`A≺0` IR / `A≻0` UV), the operators, conditioning, defect, and Jacobian all differ, and the
coupling that *helps* in the `D-BBS`/Lemma-NS convention *hurts* in the TCI/`(★)` convention — so TCI
neither implies nor is implied by `(★)`. Not `[T]`, not `[NG]`: `[C]` on `(H-C) ∧ (H-SL)`, `[O]` to close.
Clay 0%.

---

## 8. Verification / references (internal, this session)

- `research/v3/C2-BALABAN-UNIFORM.md` §2.4 (the (★) statement), §4 (TCI "IR-native"), §5 (UV acid test).
- `docs/paper-v2/ch05-theorems.tex` §Feshbach–Bogoliubov (two-block lemma with `C ≥ γI`, `Σ ≥ μI`;
  no positive-Hessian hypothesis), Lemma NS (`Σ ⪰ A` for `C≺0`), variance decomposition `G″ = S − D`.
- `research/unconditional/FORGE-BLIND.md` §2 (two-block lemma, convention A), §3 (TCI), §4 (ledger).
- `research/unconditional/D-BBS-INTEGRATE-NEG.md` (b) (Lemma NS, convention B); `FORGE-FESHBACH-ANGLE.md` §1–3.
- `research/unconditional/INVENT-VARIANCE-INEQ.md` §1–3 (exact `G″ = S − D`, `K`, `c`, BL defect).
- `research/unconditional/PROOF-JACOBIAN-SAVES.md` §1–2 (`G″ = Δ_dyn + W_geo + X`; `Φ″ ⪰ N/12`
  pointwise; Prékopa `W_geo ≥ 0`).
- `research/unconditional/G1-L3-UV-FORCE-BOUND.md` §3 (Route B cone-emptying), §6 (open items (1)(2)(3));
  `PROOF-L3.md` §1–5; `G1-UNIFORM-ZA.md` §3; `PROOF-IR-H1.md` §1 (orthogonal `∇∥`/`∇⊥` split).
- Brascamp–Lieb (1976); Prékopa (1971); Schur complement / Haynsworth inertia additivity; Kato (1966).
- Witness `research/v3/w7/s3_feshbach_witness.py` (deterministic; run `py -3.12`), output quoted in §5.

*No third-party sources were downloaded this session; all external results above are standard and cited
from the corpus files, not re-fetched.*

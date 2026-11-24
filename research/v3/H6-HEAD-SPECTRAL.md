# H6-HEAD-SPECTRAL — the spectral / transfer-matrix attack on the head (β≈2.5): no volume-uniform spectral bound reaches it, and the obstruction is *typical-point positivity vs the essential infimum* (the Savvidy double well)

> 2026-11-24 · Sixth-wave angle **H6**: attack the head of the UV tower — `β≈2.5`, typical field,
> where the Hessian modulus fails but the measure may carry positivity — with **spectral /
> transfer-matrix** methods (gap, log-Sobolev, Cheeger, block-transfer-matrix).
> Reads: `STAR2-CLASSIFICATION.md`; `PROOF-GAUSS-INVARIANT-PF.md`; `PROOF-COARSE-POINCARE.md`;
> `PROOF-JACOBIAN-SAVES.md` (the file the charter calls `PROOF-JACOBSON-SAVES.md`); `S6-DOBRUSHIN.md`;
> `docs/paper-v2/ch05-theorems.tex`. Cross-checks: `G2G-TRANSFER-SPECTRAL.md`,
> `F2-KREIN-RUTMAN-QUANT.md`, `B3-FUNCTIONAL-BE.md`, `B3-ISOPERIMETRIC.md`, `CHEEGER-CORE.md`,
> `S1-CONDITIONAL-POINCARE.md`, `OTHER-H1-POINCARE.md`, `D3-SUSCEPTIBILITY-VARIATIONAL.md`,
> `A4-SPECTRAL-BOOTSTRAP.md`, `C1-SZZ-NSCALING.md`.
> Witness: `research/v3/w7/h6_head_spectral.py` (output quoted in §6).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[V]** numerical witness · **[R]** reported/unrefereed · **[std]** standard
> textbook tool used as such (not re-verified this session).

---

## 0 · Verdict (no QED, Clay 0%)

**No spectral bound in the bank — gap, log-Sobolev, Poincaré, Cheeger, or block-transfer-matrix —
reaches `β≈2.5` with a volume-uniform constant that lower-bounds `G″`.** The attack does not fail for
lack of technique; it fails at a single structural object, and the failure is now stated in the cleanest
form the campaign has:

1. **The exact transfer-matrix / coarse-mode relation is an identity, not a bound `[T]`.** The coarse
   marginal curvature is the curvature of the coarse transfer matrix's ground state:
   `G″(a) = V_eff''(a) = −2 ∂_a² log ψ₀(a)`, where `ν = |ψ₀|²` and `ψ₀` is the Perron ground state of the
   coarse-mode block transfer operator. This *expresses* `G″` through the transfer matrix; it does not
   estimate it.

2. **The one constructive one-sided spectral inequality reads the essential infimum, not the typical
   point `[T]`.** `λ₁(ν) ≥ ess inf_a V_eff''(a) = ess inf_a G″(a)` (Bakry–Émery, sharp in 1-D / Brascamp–Lieb;
   `ν` log-concave ⇒ equality iff Gaussian). The head supplies only a **typical-point** value
   `G″(a₀*) > 0`; the inequality needs the **essential infimum over the support of `ν`**, and at the head
   that infimum is negative (the Savvidy tail). **Typical-point positivity does not lift to the essential
   infimum.** This is the precise obstruction `[T]`+`[V]`.

3. **The block transfer matrix splits into two objects, and neither closes `[T]`.** (a) The **fast fiber**
   per-block gap is positive and volume-uniform (`~Λ_UV`) but **provably silent on `G″`** — two-coordinate
   independence / fast-≠-slow (`B3-ISOPERIMETRIC`, `PROOF-COARSE-POINCARE`(b); witness §6·D). (b) The
   **coarse / marginal** block-transfer gap *is* the coarse Poincaré constant `λ₁(ν) = 1/Var_ν(a)` — i.e.
   **`H1` by definition**, not a bound on it — and it degenerates as `k→0`. There is no third object.

4. **The named no-go `[NG]`.** *No gauge-invariant-PF / Krein–Rutman–Bushell, Bakry–Émery / log-Sobolev /
   Poincaré, Cheeger / isoperimetric, or per-block-coarse block-transfer-matrix spectral bound reaches
   `β≈2.5` with a volume-uniform constant controlling `G″`.* Each reduces either to the silent fast fiber
   gap or to the coarse marginal Poincaré constant `=H1`; the bridge between them ("typical point ⇒
   essential infimum") is exactly the Savvidy double-well barrier. `[NG]` for the **method class**, not for
   the inequality `G″>0` itself.

5. **What is missing `[O]`.** A *Savvidy-corrected uniform log-concavity* of the coarse marginal on the
   typical-field bulk (a lower bound on `ess inf G″` there, tail-controlled) — equivalently the campaign's
   `(◆)=H1`; **or** a second invariant cone / total-positivity (Karlin flag) separating bulk from tail —
   absent on the orderless quotient `𝒜/G`; **or** the two-scale Feshbach/Schur transfer quotient, whose
   hypothesis is `Σ≥μ` = TCI = `H1`. None exists.

The free-field acid test **passes cleanly and non-vacuously** (`λ₁ = G″ = 1/Var → 0` for the massless
field, witness §6·A), and the `β=0` Haar measure gives a genuine theorem (`G″ = W_geo ≥ N/12 ⇒ λ₁ ≥ N/12`,
`[T]`); the head is where both stop. **The wall did not move; it is now labelled in transfer-matrix
coordinates.**

---

## 1 · The spectral setup: four objects, three exact identities

**Context (head).** At `β≈2.5`, typical field, the marginal coarse object is the *conditioned coarse
curvature* `G″(v;a₀)=E_a[∂_v²V]−Var_a(∂_vV)` (STAR2 §0). Fix the softest coarse mode (`k→0`, block-uniform)
and let `a=⟨v,PX⟩` be its amplitude. The **coarse marginal** is `ν(da)`, and

```
V_eff(a) := −log(dν/da),          G″(a) = V_eff''(a)          [k0-CONSTRAINT-EFFPOT; T]
```

is the exact marginal curvature (ch05 `eq:split`; `OTHER-H1-POINCARE`). `(★)` is `G″≥δ_W` and, by Ward
`G″·Var_ν(a)=1` (Gaussian-exact, measured 2–5%), is the same as `Var_ν(a) ≤ 1/δ_W` = the **coarse spectral
gap** (`D3`).

**The four spectral objects and what each is.**

| # | object | definition | relation to `G″` |
|---|--------|-----------|------------------|
| (i) | full transfer matrix `T=e^{−a_0H}` | OS/PF kernel on `ker Gauss` | Bushell/KR gap contracts at rate `=G″` (`F2`, `G2G`) |
| (ii) | **coarse marginal transfer** | generator `L_ν=∂²−V_eff'∂` reversible w.r.t. `ν` | Poincaré constant `λ₁(ν)=inf ℰ_ν/Var_ν` **`=1/Var_ν(a)` = H1** |
| (iii) | **per-block transfer matrix** (all block modes) | `T_block=e^{−a_0 H_block}` | fast fiber gap `~Λ_UV`, **silent** on (ii) |
| (iv) | Cheeger / LSI / log-Sobolev | `λ₁ ≥ h(ν)²/4`, `Ent ≤ (2/λ)ℰ` | `h_binding∼√λ` (tautology); LSI `⇔` `G″⪰λ` |

**Three exact identities (all `[T]`, none an estimate).**

**(H6-A) Schrödinger / Riccati identity.** The coarse marginal density is the squared ground state,
`ν = |ψ₀|²`, `ψ₀>0`, where `ψ₀` is the Perron ground state of the self-adjoint form of the coarse-mode
block transfer operator, `−ψ₀''+Uψ₀=E₀ψ₀`. Then

```
G″(a) = V_eff''(a) = −2 ∂_a² log ψ₀(a) = 2( r(a)² − (U(a)−E₀) ),      r(a):=ψ₀'(a)/ψ₀(a),
r'(a) + r(a)² = U(a) − E₀        (Riccati).                                            (1.1)
```

This is the requested **transfer-matrix / coarse-mode relation**: `G″` *is* the log-curvature of the
transfer-matrix ground state. It is an identity — it repackages `G″` in the transfer variable, and by
itself gives no positivity (the RHS can have either sign).

**(H6-B) Bakry–Émery / Brascamp–Lieb one-sided inequality.** Let `λ₁(ν) := inf_φ ℰ_ν(φ)/Var_ν(φ)` be the
coarse Poincaré (spectral-gap) constant, `ℰ_ν(φ)=∫|∇φ|²dν`. Then

```
λ₁(ν) ≥ ess inf_a V_eff''(a) = ess inf_a G″(a),        equality iff Gaussian (constant V_eff'').   (1.2)
```

(`V_eff''≥λ` a.e. ⇒ `LSI(λ)`/`PI(λ)`; sharp in 1-D and for Gaussians.) In the log-concave regime one also
has the Brascamp–Lieb form `λ₁(ν) ≥ 1/E_ν[1/V_eff''] ≥ ess inf V_eff''`. **The constant is the essential
infimum** — this single word is the whole obstruction (§3).

**(H6-C) Cheeger tautology (banked).** `λ₁(ν) ≥ h(ν)²/4`; the binding competitor is the transverse `k→0`
saddle with `h_binding ∼ 1/σ = √λ`, `σ²=Var_ν(a)=1/λ`, so `λ ≥ (√λ)²/4 = λ/4` (`CHEEGER-CORE`, `B3-ISOPERIMETRIC`).
No independent perimeter floor; `h≥h₀ ⟺ H1`.

---

## 2 · The attempted transfer-matrix / coarse-mode inequality — and why each reading fails

**The inequality we would want (the target).**

```
∃ λ>0:  G″(a₀*) > 0  at the TYPICAL field a₀*   ⟹   λ₁(ν) ≥ λ  uniformly in (L,a).        (2.1)
```

This is the natural way "the measure-carried positivity at the head" would enter a proof: the measured
`G″(a₀*)=+0.13` (STAR2) / `+0.27…+0.51` (W7) at a typical field would seed a uniform spectral gap. **(2.1)
is false in the direction used, and (H6-B) says exactly why:** the left side is a *pointwise* value, the
right side is a *global* constant controlled by the **essential infimum** `(1.2)`. A positive value at one
typical point gives no lower bound on an infimum.

**Reading (i) — full transfer / Bushell / Krein–Rutman.** This is real and gives the *only* right-sign
lower bound on a gap in the corpus (`F2`, `G2G`): `λ₁/λ₀ ≤ tanh(Δ/4) < 1 ⇒ a·m_eff ≥ −ln c(Δ)`. But the
contraction rate is `−b′(0)=G″=S−D` (`D2`; `F2 §3`), so *uniform* `c<1 ⟺ G″(0)≥λ ⟺ H1`; and at finite `V`
the PF/KR constants scale like `e^{−cV}` (`PROOF-GAUSS-INVARIANT-PF (2)`, `P4`), so no **volume-uniform**
constant follows: `inf_V Δ_V>0` *is* the mass gap and is not obtained at `β=2.5`. `[T]` finite-`V`,
`[O]` uniform at the head.

**Reading (ii) — coarse marginal transfer (Poincaré `λ₁(ν)`).** The coarse Dirichlet form is exact,
`ℰ_ν(φ)=ℰ_μ(φ∘P)`, and `Hess V_eff = G″`, so on flat coarse space Bakry–Émery is *lossless*:
`λ₁(ν) ≥ λ ⟺ G″⪰λ ⟺ H1` (`B3-FUNCTIONAL-BE`, `OTHER-H1-POINCARE`). This is a **faithful rename** of `H1`,
not a bound on it. The only two constructions of a coarse form are (A) push the full LSI through `P` —
needs the `2.5` gap (circular); (B) Bakry–Émery on `ν` — needs `Hess V_eff⪰λ` (H1 itself) (`OTHER-H1-POINCARE`).
Exhausted.

**Reading (iii) — per-block (fiber) transfer matrix.** This gap is positive and volume-uniform
(`~Λ_UV ~ 1/a²`), the tempting object. It is **silent on `G″`**: the joint is `μ = ν(da)·μ(·|a)`, and the
marginal `ν` and the fiber conditionals are *independent coordinates*, so no fiber spectral property
constrains `ν` (`B3-ISOPERIMETRIC`, the two-coordinate obstruction). Quantitatively, the fiber gap is a
**fast relaxation rate**, not the slow coarse stiffness, and the only variance bound it supplies runs the
**wrong way** by `≈e^{5N_p}` (`PROOF-COARSE-POINCARE`(b), H2). Witness §6·D makes this exact: an
integrable Gaussian with slow precision `λ_s` and fast precision `Λ` has `Var(a)=1/λ_s` *independent of
`Λ`* — the fast gap can be sent to `∞` with **zero** effect on `G″_marg=λ_s`. So the per-block spectral gap,
though positive and uniform, certifies nothing about the coarse curvature.

**Reading (iv) — Cheeger / isoperimetry.** Curvature-free in *hypothesis* (a genuine advantage over
Bakry–Émery, which dies literally at `Γ₂<0`), but isoperimetry is a *bottleneck* functional, and the
Savvidy negativity is a **saddle**, not an expander: it *spreads* the measure along the very direction that
defines the cheap cut (`CHEEGER-CORE §2`). `h_binding∼√λ` ⇒ tautology; `h≥h₀ ⟺ Var(a)≤1/λ ⟺ H1`.

**Conclusion of §2.** Every spectral reading either (a) is the silent fast fiber gap, (b) is `H1` renamed
(coarse marginal Poincaré), or (c) is a tautology (Cheeger). None is a bound `[T]`. The failure is uniform
across the class — this is the named `[NG]` of §4.

---

## 3 · The obstruction, precisely: typical point vs essential infimum (Savvidy double well)

The common failure has one geometric name. At the head, the softest coarse mode is the lattice **Savvidy
operator** `−D²+2F` at `k→0`, whose quadratic part is tachyonic (`m²=k²−2g|F|<0` at the crossover), stabilised
at large amplitude by the SU(N) compactness. Its 1-D normal form is the **double well**

```
V_eff(y) = ½ m² y² + ¼ κ y⁴ ,   m² < 0  ⇒  wells y₀=√(−m²/κ),  barrier ΔV = m⁴/(4κ) > 0.
G″(well) = −2m² > 0  (TYPICAL positivity);   ess inf G″ = m² < 0  (barrier: density minimum).   (3.1)
```

This is the exact structure `S1-CONDITIONAL-POINCARE §4.2` uses for the conditioned soft mode ("the
tachyonic quadratic is stabilized by the quartic ⇒ a double well once the field passes the critical
value"). On it, the two sides of `(2.1)` decouple:

- **Typical point:** `G″(a₀*)=−2m²>0` at each well (the "measure-carried positivity": at typical field the
  conditioned mean `E[H|0]=+0.22` plus `W_full≈+0.54` make `G″` positive, `PROOF-JACOBIAN-SAVES`, W7).
- **Essential infimum:** `ess inf G″ = m²<0` at the barrier; the Haar floor `N/12` is `|F|`-free and does
  *not* cover the `∝|F|` well depth (`B3-FUNCTIONAL-BE §3`; `Φ″`'s convexity sits at the log-cut where
  `F≈0`, anti-localised from the Savvidy driver).
- **Spectral constant:** `λ₁(ν) ~ (prefactor)·e^{−ΔV} → 0` as the barrier grows (Kramers metastability).
  So `λ₁(ν)` is **exponentially smaller** than `G″(well)` — typical-point positivity and the spectral
  constant differ by the Kramers factor.

Witness §6·B quantifies it on `V_eff=(g/2)(x²−1)²` (the symmetric Savvidy normal form, `ΔV=g/2`):
`G″(well)=4g` grows, while `λ₁` falls `0.79→0.12` over `g=0.5→8` and the ratio `λ₁/G″(well)`
`0.40→0.0037 → 0` — a clean Kramers collapse at **fixed** positive typical curvature. In the physical
variable this is worse: as the mode approaches the exact crossover (`m²→0⁻`, `k→0`), `ΔV=m⁴/(4κ)→0` but so
does the typical curvature `−2m²→0`, and the *tower* descent (`k→0`) drives `λ₁(ν)→0` regardless
(`OTHER-H1-POINCARE`; `PROOF-COARSE-POINCARE`(a)). Either way: **no uniform `λ` survives the tail.**

**The precise obstruction, stated once.**
```
λ₁(ν)  ≥  ess inf_a  G″(a)                    [H6-B; the ONLY constructive one-sided bound]
λ₁(ν)  =  1/Var_ν(a)   (Gaussian)   =  H1      [coarse marginal transfer = Poincaré]
head:    G″(a₀*) = −2m² > 0  (typical point, a₀* in bulk of ν)
         ess inf_a G″(a) = m² < 0  (Savvidy barrier)   and   λ₁(ν) ~ e^{−ΔV} → 0
⇒ "typical point > 0"  does NOT imply  "ess inf > 0"  nor  "λ₁ ≥ λ>0".
```
The map **typical-point value → essential infimum** is exactly where the head lives, and it is blocked by
the Savvidy double well. This is the transfer-matrix face of `H1` — a fifth face, distinct from the banked
*fast-vs-slow* (`PROOF-COARSE-POINCARE`b), *fiber-vs-marginal* (`B3-ISOPERIMETRIC`), *harmonic-vs-arithmetic*
(`D3`), and *saddle-neck* (`CHEEGER-CORE`) faces; it specialises to each of them but its native statement is
`ess inf vs typical point`.

---

## 4 · The named no-go `[NG]`

> **`[NG]` (spectral method class).** *No bound in the class*
> *{gauge-invariant Perron–Frobenius / Krein–Rutman–Bushell; Bakry–Émery / log-Sobolev / Poincaré on the
> coarse marginal; Cheeger / isoperimetric / capacity; per-block or coarse block-transfer-matrix spectral
> gap}*
> *reaches `β≈2.5` with a constant uniform in `(L,a)` that lower-bounds `G″`.*

**Proof of the no-go (classification, not a new inequality).** Each member reduces to one of two objects
(§2): the **fast fiber** gap (positive, volume-uniform, but silent on `G″` by the two-coordinate
independence — `B3-ISOPERIMETRIC`), or the **coarse marginal** Poincaré constant `λ₁(ν)` (`=H1` by H6-B's
sharpness and `OTHER-H1-POINCARE`). The only bridge from a positive *fast* or *typical-point* constant to a
uniform lower bound on `λ₁(ν)` is the essential infimum `(1.2)`, which the Savvidy double well makes
negative while `k→0` makes `λ₁(ν)→0` (§3). Hence no member closes. ∎

**Honest scope.** This is `[NG]` for the **method class**, not for `(★)`. A non-log-concave measure can
still have a positive Poincaré constant (a symmetric double well with a *fixed* finite barrier does); what
the no-go refutes is that these *tool-derived* constants reach the head **uniformly**. Also, the `[NG]`
is *off the relevant scope* for the banked positive results: the spectral route **does** give
`[T]` volume-uniform constants at `β<1/48` (SZZ LSI, Lemma SC) and `[T]` at `β=0` (Haar), and `[C]` in the
deep UV (cone-emptiness, `S6 §5`) — all of which sit at `β_j<1/48` or `β_j→∞`, disjoint from the tower's
head (`STAR2 §5`, `S6 §4`).

**Why the class is exactly the one to name.** It is the union of every "spectral/energy" tool the campaign
has tried (PF, KR, Bushell, Bakry–Émery, Brascamp–Lieb, Cheeger, LSI, capacity, RCD): each was scored `2/5`
in its own note (`CHEEGER-CORE` scores all four isoperimetric siblings `2/5`; `G2G` 3/5; `F2` 3/5), and all
fail on the *same* two of five points (`#3 non-circular`, `#5 uniform`). The no-go unifies that pattern: the
missing step is always the same object, `ess inf G″ ≥ λ` on the head.

---

## 5 · Acid tests

**(A) Free field — the bound MUST return `0`. It does, cleanly and non-vacuously. `[T]`+`[V]`.**
Free massless Gaussian coarse marginal `V_eff=½ k² x²`: `G″=k²`, `Var=1/k²`, and the transfer gap
`λ₁(ν)=k²` **exactly** (Gaussian: `(1.2)` is equality). As `k→0`, `λ₁=G″→0` (witness §6·A: `λ=1,0.5,0.1,0.02`
give `λ₁=1.000,0.500,0.100,0.020`). The tool manufactures **no false positivity** on the gapless free field
— the charter acid test is passed. (Contrast: Cheeger also returns the correct `0`, in the *continuum*,
`CHEEGER-CORE §3`.)

**(B) `β=0` Haar — genuine theorem `[T]`.** `V=S=Φ`, `Φ` convex, `Φ''⪰N/12` pointwise (Prékopa /
`PROOF-JACOBIAN-SAVES` Lemma B), and the marginal is log-concave ⇒ `(1.2)` is `λ₁(ν) ≥ N/12`, i.e.
`Var_ν(a) ≤ 12/N`, unconditional, all `L`, all modes (`D3 §2`, the `χ≤12/N` theorem). This is where the
transfer route genuinely closes — `β=0`, not the head.

**(C) `β<1/48` — volume-uniform but off-tower `[T]`.** Shen–Zhu–Zhu give a volume-uniform LSI with
constant `K_S=N(½−24|β|)>0` (`C1-SZZ-NSCALING`; `[V]` source-verified), pushed through the 1-Lipschitz block
map to the coarse Poincaré (Lemma SC). This *is* a volume-uniform spectral bound at `β<1/48` — and it is
exactly where the full LSI already holds, disjoint from `β≈2.5` (`STAR2 §5`, `S6 §4`). It does not transport.

**(D) Savvidy double well — the head mechanism `[T]`+`[V]`.** Witness §6·B–C: `G″(well)=4g>0`,
`ess inf G″=−2g<0`, `λ₁∼e^{−g/2}`; the one-sided bound `(1.2)` is sharp for the Gaussian and **vacuous**
(RHS `<0`) for the double well. The exact mechanism by which typical-point positivity fails to lift.

---

## 6 · Witness (`research/v3/w7/h6_head_spectral.py`, deterministic, `py -3.12`, numpy+scipy)

Quoted output (numpy 2.4, scipy P1-FEM generalized eigenproblem `A c = λ M c` for the generator on `L²(ν)`):

**(A) free field** — `λ₁` vs `G″`:
```
   lam   lambda_1(fem)      G''=lam    Var=1/lam
 1.000        0.999995      1.000       1.000
 0.500        0.499999      0.500       2.000
 0.100        0.100000      0.100      10.000
 0.020        0.020000      0.020      50.000     => λ₁=G″=1/Var → 0 (ACID PASS)
```
**(B) double well** `V_eff=(g/2)(x²−1)²` — typical `G″` grows, `λ₁` collapses (Kramers):
```
   g    G''(well)   ess inf G''   lambda_1     Var(x)    lam1/G''well
 0.5      2.000      -1.000      7.921e-01    1.042      3.960e-01
 1.0      4.000      -2.000      8.459e-01    0.893      2.115e-01
 2.0      8.000      -4.000      7.487e-01    0.833      9.359e-02
 4.0     16.000      -8.000      4.545e-01    0.852      2.841e-02
 8.0     32.000     -16.000      1.191e-01    0.918      3.722e-03   => ratio → 0
```
**(C) the one-sided inequality** `λ₁ ≥ ess inf V_eff''`: Gaussian `λ₁=1.000=ess inf` (sharp);
double well `g=2`: `λ₁=0.749 ≥ −4` holds but **vacuous** (RHS `<0`).
**(D) two-coordinate block transfer** (exact Gaussian, `λ_s=0.4`): `G″_marg=1/Var=0.4` **independent of**
fast precision `Λ` (1→1000) and coupling `b`; the full 2D rate is `λ_s/(1+b²)≠G″` and also `Λ`-independent
as `Λ→∞`. The fast fiber gap is silent; the only block-`T` gap that *equals* `G″` is the marginal one `=H1`.
**(E) susceptibility direction**: test function `φ=a` gives only `λ₁ ≤ 1/Var` (an **upper** bound on the
gap); Gaussian equality `1/λ₁=Var`, double well `1/λ₁` grows while `Var=O(1)` ⇒ `χ` cannot lower-bound the
gap (confirms `D3`).

*Honest scope of the witness:* (A),(D),(E) are exact linear algebra; (B),(C) are 1-D FEM eigenvalue
computations of the two structural claims (Kramers collapse, essential-infimum vacuity). It is a
**witness, not a proof**; the theorem content is H6-A/H6-B (§1), and the open step is §3's head.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| exact marginal curvature | `G″(a)=V_eff''(a)`, `ν` coarse marginal | **[T]** (k0-CONSTRAINT; carried) |
| transfer/Schrödinger identity (H6-A) | `G″=−2∂²log ψ₀`, Riccati `r'+r²=U−E₀` | **[T]** |
| BE/BL one-sided (H6-B) | `λ₁(ν) ≥ ess inf_a G″(a)`; equality iff Gaussian | **[T]** (BE 1-D sharp / BL) |
| Cheeger tautology (H6-C) | `λ₁≥h²/4`, `h_binding∼√λ` ⇒ `λ≥λ/4` | **[T]** (carried `CHEEGER-CORE`) |
| coarse marginal transfer = Poincaré = H1 | lossless on flat coarse space | **[T]** (other-H1, B3-FE) |
| fast fiber block-`T` gap silent on `G″` | two-coordinate independence; wrong-direction H2 | **[T]** (`B3-ISO`, coarse-Poincaré b; §6·D) |
| `β=0`: `λ₁ ≥ N/12` theorem | `Φ` convex, `G″=W_geo≥N/12` | **[T]** (D3) |
| `β<1/48`: volume-uniform LSI (SZZ) → coarse Poincaré | Lemma SC | **[T]** (SZZ `[V]`) |
| typical point `G″(a₀*)>0` at head | `+0.13` (STAR2) / `+0.27…+0.51` (W7) | **[V]** |
| `ess inf_a G″<0` on Savvidy tail | double well `m²<0`; `Φ''` floor `|F|`-free | **[T-standard]** + **[V]** |
| **no spectral bound reaches `β≈2.5` volume-uniformly** | the class of §4 | **[NG]** |
| **uniform `ess inf G″ ≥ λ` on typical bulk** | the missing estimate (=H1) | **[O]** |
| free-field `λ₁=G″=0` | exact, non-vacuous | **[T]**+**[V]** |
| **`(★)` uniform in `(j,U′,a)`** | the campaign's L3 | **[O]** — not closed |

---

## 8 · What would be needed instead `[O]`

The transfer-matrix route is exhausted; what it asks for is one of three structurally distinct objects,
none of which exists:

1. **A Savvidy-corrected uniform log-concavity of the coarse marginal on the typical-field bulk** — a
   lower bound `ess inf_{a∈bulk} G″(a) ≥ λ`, *tail-controlled* so the quadratically-suppressed Savvidy cone
   (`ν`-measure `≤ e^{−c/g²}`, Lemma B) does not spoil the infimum. This is the campaign's `(◆)=H1`
   (`PROOF-JACOBIAN-SAVES §2`, `OTHER-H1-POINCARE`).
2. **A second invariant cone / total-positivity (Karlin flag)** on `𝒜/G` that places the typical bulk in
   the interior and the Savvidy tail on the boundary of a *second* cone below positivity — absent because
   `𝒜/G` carries **no order** (`PROOF-F1SQ iii`; `INVENT-SECOND-CONE`; `F2 §3`); the only invariant cone is
   the full positive one, on whose boundary the `0⁺⁺` excitation lives (`PROOF-GAUSS-INVARIANT-PF`).
3. **The two-scale Feshbach/Schur transfer quotient** — eliminating the fast (positive, gapped) fiber to
   raise the slow coarse operator. This *is* the Feshbach–Bogoliubov lemma (ch05 §`sec:feshbach`), and its
   hypothesis is `Σ≥μ` = the tachyonic-compensation inequality `(TCI)` = `H1`; the fast elimination
   certifies only the fast half (`Lemma NS`), never the slow variance (`S1 §6`).

So the transfer-matrix attack sharpens the wall into a single named object and supplies no residue of its
own.

---

## 9 · Verification ledger

| item | checked | status |
|---|---|---|
| Bakry–Émery `V''≥ρ ⇒ LSI(ρ)`, 1-D sharp | standard | **[std]** |
| Brascamp–Lieb `Var ≤ E[1/V_eff'']` (1-D log-concave) | standard (J. Funct. Anal. 22 (1976) 366) | **[std]** |
| Cheeger `λ₁≥h²/4`, Buser converse needs Ric | carried from `CHEEGER-CORE` | **[std]** |
| Krein–Rutman / Bushell `λ₁≤cλ₀`, `c=tanh(Δ/4)` | carried from `F2`/`G2G` | **[std]** |
| Kramers metastability `λ₁∼e^{−ΔV}` | standard (Bovier et al.); witness §6·B | **[std]**+**[V]** |
| Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737, `β<1/48`, `K_S=N(½−24|β|)` | carried, source-verified in `C1` | **[V]** (carried) |
| Prékopa; Lemma B `Φ''≥N/12`, `W_geo≥0` | carried, re-verified `S1 §5C` | **[T]** (carried) |
| Savvidy `−D²+2F`, `λ_min=−g|B|`, double-well normal form | carried (`S1 §4.2` exact model) | **[T-standard]** |
| witness `research/v3/w7/h6_head_spectral.py` | ran `py -3.12`, exit 0 | **[V]** |

**Notes.** (i) All external "standard" results (BE, BL, Cheeger, Buser, KR/Bushell, Kramers, Prékopa) are
used as textbook tools and marked `[std]`, not `[V]` — none downloaded this session (charter, and S6's
same convention). (ii) The normalisation throughout is the standard-Wilson one of the head `β≈2.5`
(`S6 §2`); the `1/48` threshold is the strong-coupling LSI radius, disjoint from the tower (`STAR2 §5`).
(iii) The double-well `(3.1)` is the exact normal form `S1-CONDITIONAL-POINCARE §4.2` already uses for the
conditioned soft mode; the witness's `λ₁` numbers are 1-D FEM, `G″(well)`/`ess inf` are exact.

---

## 10 · Bottom line (no QED)

The spectral / transfer-matrix attack **does not reach the head**, and it now says exactly why, in the
cleanest available coordinate:

> **`G″` *is* the curvature of the transfer-matrix ground state (identity H6-A). The only constructive
> spectral inequality reads the *essential infimum* of that curvature (H6-B). The head gives a positive
> curvature at a *typical point* but a negative essential infimum on the Savvidy *tail*; the fast per-block
> spectral gap is positive and uniform but provably silent on `G″`. Hence no gap/LSI/Cheeger/block-transfer
> bound is volume-uniform at `β≈2.5`. Typical-point positivity does not lift to a global spectral
> constant — that lift is `H1`.**

The free-field acid test is passed cleanly (`λ₁=G″=0`), `β=0` and `β<1/48` give genuine volume-uniform
theorems `[T]`, and the head is `[O]`. **Clay = 0%. The wall is unmoved, newly labelled.**

*Refs.* `STAR2-CLASSIFICATION.md`; `PROOF-GAUSS-INVARIANT-PF.md`; `PROOF-COARSE-POINCARE.md`;
`PROOF-JACOBIAN-SAVES.md`; `S6-DOBRUSHIN.md`; `docs/paper-v2/ch05-theorems.tex` (`eq:split`, `sec:feshbach`);
`G2G-TRANSFER-SPECTRAL.md`; `F2-KREIN-RUTMAN-QUANT.md`; `B3-FUNCTIONAL-BE.md`; `B3-ISOPERIMETRIC.md`;
`CHEEGER-CORE.md`; `S1-CONDITIONAL-POINCARE.md`; `OTHER-H1-POINCARE.md`; `D3-SUSCEPTIBILITY-VARIATIONAL.md`;
`A4-SPECTRAL-BOOTSTRAP.md`; `C1-SZZ-NSCALING.md`. External (standard): Bakry–Émery 1985; Brascamp–Lieb,
Adv. Math. **22** (1976) 366; Cheeger, *Problems in Analysis* (1970) 195; Buser, Ann. Sci. ENS **15**
(1982) 213; Krein–Rutman 1948; Birkhoff, Trans. AMS **85** (1957) 219; Bushell, Arch. Rational Mech. Anal.
**52** (1973) 330; Kramers 1940 / Bovier–Eckhoff–Gayrard–Klein (metastability); Savvidy, Phys. Lett. B **71**
(1977) 133; Nielsen–Olesen, Nucl. Phys. B **136** (1978) 45; Prékopa, Acta Sci. Math. (Szeged) **32** (1971)
301; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737. Witness: `research/v3/w7/h6_head_spectral.py`.

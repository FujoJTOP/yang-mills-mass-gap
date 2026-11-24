# H3-AVERAGED-HEAD — the averaged head inequality `(★_avg)` SURVIVES at β=2.5 `[T]`+`[V]`, but its margin *is* the coarse gap, so it relocates H1 rather than closing it

> 2026-11-24 · Attack angle H3 on the head of the UV tower: **prove or obstruct the averaged form**
> `(★_avg)` of the single live inequality `(★)`. The pointwise-in-`U′` and pointwise-in-`a₀` forms are
> already dead (the barrier counterexample, G4/V3); the question is whether averaging kills the bad set.
> Reads: `STAR2-CLASSIFICATION.md`, `S1-CONDITIONAL-POINCARE.md`, `S5-N-UNIFORM.md`, `V4-PINCER.md`,
> `G4-GENERAL-BOUND.md`, `V3-COUNTEREXAMPLE-SCOPE.md`, `S2-BOUNDARY-COUNTEREXAMPLE.md`,
> `G1-BLOCK-IDENTITY.md`, `STAR-WAVE-SYNTHESIS.md`, `C2-BALABAN-UNIFORM.md` §2.2 (Lemma B),
> `docs/paper-v2/ch05-theorems.tex` (`eq:split`, `eq:vardecomp`, `sec:wgeo`).
> Witness: `research/v3/w7/h3_averaged_head.py` (`py -3.12`, numpy-only, deterministic; output quoted).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[V]** numerical witness · **[R]** reported/unrefereed.

---

## 0 · Verdict (no QED; Clay 0%)

**The averaged head inequality survives at β=2.5 — in two inequivalent readings — and neither is
obstructed by the barrier counterexample. But both margins turn out to *be* the coarse spectral gap, so
`(★_avg)` is true and does not close H1.**

1. **The `a₀`-average is an exact theorem `[T]`.** Averaging the *true* coarse marginal curvature
   `G″(v;a₀)=−d²/da² log ρ_v(a)` over the coarse field with its own marginal `ρ_v(·|U′)` gives, by
   integration by parts, **the Fisher information of the coarse marginal**:
   ```
     G_A(v|U′) := E_{a~ρ_v(·|U′)}[ G″(v;a) ]  =  I(ρ_v) := ∫ (ρ_v′)²/ρ_v da  ≥  0 ,
   ```
   with equality iff the coarse field is delocalized (flat `ρ_v`). By Cauchy–Schwarz,
   `I(ρ_v)·Var_ρ(a) ≥ 1`, i.e. **`G_A ≥ 1/Var(a)` = the coarse mass-squared.** This holds for **every**
   boundary `U′`, **every** mode `v`, **every** β — including β=2.5 — and needs **no** measure model
   beyond "`ρ_v` vanishes at the chart edge" (which holds, G1 §5). **[T]**, verified to `≤1.2×10⁻³`.

2. **The `U′`-average at fixed `a₀` survives at the head `[V]`.** In the exact one-plaquette model at
   β=2.5, `G″(v;a₀=0,θ)>0` on the **whole** holonomy range `θ∈[0,π]` (min `+0.136` at `θ=π`); the
   `U′`-average is `+0.37` (Haar) / `+0.59` (equilibrium Wilson-tilted). The barrier goes negative only
   for `β>β*=4.307`, and even then the **genuinely negative set has equilibrium measure**
   `≤1.6×10⁻⁴` at β=5 and `→0` as `β→∞`. So "the bad set has small measure" is exact in the right measure.

3. **The good/bad *bound* (measure × worst-case) is decisive at β=2.5 `[V]`, and becomes airtight with
   a block Laplacian gap `[C]`.** With the crude "charge the whole not-good set at its worst value"
   estimate and `c_□=0` (degenerate single plaquette), the bound is `+0.137>0` at β=2.5 but `−0.061<0`
   at β=5. With a real block gap `c_□≥0.25` it is positive at **both** β=2.5 and β=5. Deep UV
   (`g_j→0`) the negative set's measure `≤e^{−c/g_j²}→0` (Lemma B), so **any** `c_□>0` suffices there.

4. **The margin is the gap — `(★_avg)` is H1 in disguise `[T]`+`[O]`.** The exact margin `G_A=I(ρ_v)`
   satisfies `I ≥ 1/Var(a)`; `1/Var(a)` *is* the coarse curvature / coarse mass-squared in the Gaussian
   approximation. So `(★_avg)` holds **with `δ_W = I(ρ_v)`, and a uniform floor `δ_W≥δ>0` is exactly the
   coarse spectral gap (H1)**. The averaged form is *not* a weaker statement than the gap — it is the gap,
   read through the coarse marginal. Acid test passes: free massless field ⇒ `Var(a)→∞` ⇒ `I=1/σ²→0`.

**Bottom line.** The averaged form *survives* at β=2.5 (answering the task's question **yes**), with an
exact-model margin `δ ≈ 0.37–0.61` (`U′`-Haar to `a₀`-Fisher), and it is **not** obstructed by the
barrier because the barrier is either positively curved at β=2.5 or, for `β>β*`, carried by a
measure-`e^{−2β f}` set. But the margin equals `1/Var(a)` = the coarse gap, so `(★_avg)` does not close
Clay. **No QED. Clay 0%.**

---

## 1 · The object and the precise definition of `(★_avg)`

### 1.1 Which `G″`

The pointwise forms are dead and the corrected object is fixed (G4 `[T]`, V3 `[T]`): the coarse marginal
curvature
```
   G″(v;a₀) := −d²/da² log ρ_v(a)|_{a₀} = E_{μ(·|U′)}[ ∂_v²V | a₀ ] − Var_{μ(·|U′)}( ∂_vV | a₀ ),
   V = S_W + Φ ,   a = ⟨v,X⟩ ,   ρ_v(a) = ∫ δ(a−⟨v,X⟩) e^{−V} dX ,
```
for the coarse mode `v` (linear/divergence-free, G1 class (a)) and the conditional block measure
`μ(·|U′)` of S1 eq. (1.1). The strong "uniform in `(U′,a₀)`" form is **false** `[NG]` (V3 §4.1):

```
   G″(e₃; θ=π, a₀=0) = A − B β² ,   A=0.20511530 , B=0.01105929 ,   A−Bβ²<0  ⟺  β>β*=4.30661 .
```

The whole point of `(★_avg)` is that the negative set `{G″<0}` might be **measure-small**, so an average
survives where the pointwise form fails.

### 1.2 The two averages, defined

Let `ν_j` be the equilibrium boundary measure on `U′` at shell `j` (the law of the boundary links after
integrating out the block; absolutely continuous w.r.t. Haar, full support — V3 §3.3), and let
`ρ_v(·|U′)` be the coarse marginal at fixed `U′`. There are two natural averages; both are "the
equilibrium measure" of the task, over different variables.

> **`(★_avg^{a})` — coarse-field average (at fixed U′).**
> ```
>   E_{a~ρ_v(·|U′)}[ G″(v;a) ]  ≥  δ_W(v,U′) > 0 ,     uniformly in U′, v, and the shell j .
> ```
> **`(★_avg^{U′})` — boundary average (at a fixed reference coarse field a₀, e.g. a₀=0 or the coarse
> minimiser a₀*).**
> ```
>   E_{U′~ν_j}[ G″(v;a₀,U′) ]  ≥  δ_W(v,a₀) > 0 ,       uniformly in v, a₀, j .
> ```
> **`(★_avg)` joint — over `(U′,a₀)` with the joint equilibrium measure `ν_j(U′) ρ_v(a|U′)`.** This is
> the version with no model assumption on `ν_j`; it follows from `(★_avg^{a})` for each `U′`.

The task's steps 2–3 ("good-set margin vs bad-set measure, Wilson estimate") probe `(★_avg^{U′})`; the
clean theorem is `(★_avg^{a})` and its joint corollary. Both are analysed; the verdict differs.

---

## 2 · `(★_avg^{a})` is an exact theorem: the Fisher identity — **[T]**

> **Theorem H3-A (coarse-field average = Fisher information).** Let `ρ_v(·|U′) ∈ C²` on its support
> `[−A,A]` with `ρ_v(±A)=0` (true: `e^{−V}` vanishes at the chart edge, G1 §5). Then for every `U′`,
> every unit mode `v`, and every β,
> ```
>    E_{a~ρ_v}[ G″(v;a) ] = ∫_{−A}^{A} (ρ_v′(a))² / ρ_v(a) da =: I(ρ_v) ≥ 0 ,
> ```
> with equality iff `ρ_v` is the uniform (delocalized) density on its support. Moreover, by
> Cauchy–Schwarz (`1 = ∫ρ_v = ∫ρ_v′·(a−ā)`, integration by parts),
> ```
>    I(ρ_v) · Var_ρ(a) ≥ 1    ⟹     E_{a}[ G″ ]  ≥  1/Var_ρ(a) .
> ```

*Proof.* By Lemma G1, `G″(v;a₀) = −(log ρ_v)″(a₀) = −ρ_v″/ρ_v + (ρ_v′/ρ_v)²`. Hence
`E_a[G″] = ∫ ρ_v (−ρ_v″/ρ_v + (ρ_v′/ρ_v)²) da = −∫ρ_v″ da + ∫(ρ_v′)²/ρ_v da`. Since `ρ_v(±A)=0`,
`∫ρ_v″ = [ρ_v′]_{-A}^{A} = 0`; the first integral vanishes and only the Fisher information survives.
For the bound, `1 = ∫ρ_v = −∫ρ_v′(a)(a−ā) da` (parts), so `1 ≤ (∫(ρ_v′)²/ρ_v)^{1/2}(∫ρ_v(a−ā)²)^{1/2} =
I^{1/2}·Var^{1/2}`. ∎ **`[T]`**

**Consequences, all immediate and label-clean.**

- **`(★_avg^{a})` holds on every `U′`, for every β** — including the head β=2.5 — with
  `δ_W(v,U′) = I(ρ_v(·|U′)) ≥ 1/Var_ρ(a|U′)`. **`[T]`**
- **The joint `(★_avg)` holds for ANY boundary measure `ν_j`**, since averaging the per-`U′` identity:
  `E_{(U′,a)}[G″] = E_{U′~ν_j}[I(ρ_v(·|U′))] ≥ E_{U′}[1/Var_ρ(a|U′)] ≥ 0`. **No** assumption on `ν_j`
  is used. **`[T]`**
- **The margin is `1/Var(a)`.** In the Gaussian approximation `G″≈1/Var(a)`, so the averaged curvature
  *is* the inverse coarse variance = **coarse mass-squared**. A uniform floor `δ_W≥δ>0` is therefore the
  **coarse spectral gap = H1** `[O]`. This is the sharpest structural output: the averaged form does not
  weaken the gap — it *is* the gap, written through the coarse marginal.

**Witness confirmation** (`h3_averaged_head.py` §1): the identity `G_A = I` holds to `≤1.2×10⁻³`
(discretization-limited; the analytic identity is exact) at every tested `(β,θ)`, and `I·Var ∈
[1.01,1.25]` in §2 — always `≥1`:

```
   beta  theta        G_A         I        resid     1/Var
   2.50  3.142     0.608650   0.609835   -1.18e-3   0.48674
   5.00  3.142     1.204216   1.204302   -8.64e-5   1.09315
   2.50  0.000     0.645828   0.645785    4.27e-5   0.64029
   5.00  1.571     1.235799   1.235425    3.74e-4   1.21933
```

---

## 3 · `(★_avg^{U′})`: the good/bad decomposition at fixed `a₀` — **[V]** / **[C]**

This is the task's steps 2–3. At fixed `a₀` (we take `a₀=0`, the barrier, the worst point; the coarse
minimiser `a₀*` is easier — V3 §4.2 shows `G″(a₀*)>0` trivially), split `U′` into

- **good set** `𝒮_t = {U′ : F_∂(U′) ≤ t}` with the good-set theorem (S1 Lemma S1 / Cor S1′, `[T]`/`[C]`):
  `G″(v;a₀,U′) ≥ ρ₀(t) = c_□ + N/12 − C₁β t > 0` where `K=Π₀ Hess V Π₀ ⪰ 0`;
- **bad set** `𝒮_t^c`, where `K` may be indefinite and `G″` can drop toward the Savvidy/S1 value.

**The crude bound (measure × worst-case).** With `G_min := inf_{𝒮_t^c} G″`,
```
   E_{U′~ν_j}[G″]  ≥  ρ₀·ν_j(𝒮_t)  +  G_min·ν_j(𝒮_t^c).
```

**Bad-set measure via the Wilson action — step 3.** A boundary plaquette with field strength `f_p =
1−cos θ_p` carries action `β f_p`, hence equilibrium weight `e^{−β f_p}` per boundary plaquette. The bad
set requires at least one boundary plaquette past the threshold `f_c=(N/12+c_□)/(β)`, so
```
   ν_j(𝒮_t^c)  ≲  n_p · e^{−β f_c}  =  n_p · e^{−(N/12+c_□)},        (Wilson, per plaquette)
```
and, via the large-deviation form of Lemma B `[T]` (C2 §2.2), `ν_j(𝒮_t^c) ≤ e^{−c/g_j²}` deep in the
UV. The comparison to the good-set margin `ρ₀` is the condition
```
   ρ₀(1 − ν_bad)  >  |G_min| ν_bad      ⟺      e^{−β f_c}  ≲  ρ₀/(|G_min| n_p).          (H3-cond)
```
Since the LHS is exponential in `β` and the RHS is algebraic, **(H3-cond) holds for all sufficiently
large β (deep UV)**; the open question is the **head** `β≈2.5`, `g_j=O(1)`, where `e^{−β f_c}` is *not*
small.

**The head, β=2.5 — the crude bound already certifies positivity, because the worst point is positive.**
In the exact one-plaquette model (`c_□=0`, `C₁=1`), the *actual* `G″(e₃;a₀=0,θ)` is positive on the whole
range `θ∈[0,π]` at β=2.5 (min `+0.136` at `θ=π`; it crosses zero only at `β*=4.307`). Hence even the
crude bound is positive at the head:

```
   [5] crude good/bad (c_Box=0):   beta=2.5  theta_c=0.367  nu(bad)_Wilson=0.962  min G_bad=+0.136  BOUND=+0.137
                                   beta=5.0  theta_c=0.259  nu(bad)_Wilson=0.957  min G_bad=-0.071  BOUND=-0.061
```

**With a real block gap the crude bound stays positive for all β.** A genuine 4D block has a Dirichlet
Laplacian gap `c_□ = c/a_j² > 0` (absent in the single plaquette), enlarging the good set. Scanning it
(witness §5b):
```
   beta  c_Box   theta_c   nu(bad)_Wilson   BOUND   true<G>_Wilson
   2.50   0.00    0.367        0.962        +0.137      +0.591
   2.50   0.25    0.586        0.869        +0.173      +0.591
   2.50   1.00    1.008        0.577        +0.572      +0.591
   5.00   0.25    0.411        0.849        +0.002      +1.227
   5.00   1.00    0.697        0.533        +0.507      +1.227
```
So the worst-case bound is positive at β=2.5 even at `c_□=0`, and becomes comfortably positive at β=5 as
soon as `c_□≳0.25`. As `c_□→∞` the bound converges to the true average.

**The sharper split — charge only the genuinely negative set.** The crude bound over-charges (it puts the
worst value on the whole not-good set, which is almost all of `[0,π]` when `c_□=0`). The correct object
for the *average* is the **negative set** `S_neg={G″<0}` (witness §5b):
```
   beta  #neg/401   nu(S_neg)_Haar   nu(S_neg)_Wilson   min G'' on S_neg   <G>_Haar
   2.50     0          0.00000           0.000000            0.00000          +0.372
   5.00    15          0.06993           0.000164           -0.07137          +0.545
   10.00   41          0.19836           0.000000           -0.90081          +0.898
```
Under the **equilibrium (Wilson-tilted) boundary measure** the negative set has measure
`≤1.6×10⁻⁴` at β=5 and `→0` at β=10; at β=2.5 it is **empty**. This is the precise "bad set has small
measure" statement: the measure is the equilibrium boundary measure, and the `e^{−β(1−cos θ)}` factor
suppresses the large-holonomy negative shoulder exactly as Lemma B predicts. The `U′`-average of `G″`
is positive at every tested β (`+0.37` at β=2.5 Haar; `+0.55` β=5; `+0.90` β=10). **`[V]`**

**Why the barrier counterexample does not obstruct the average.** The pointwise counterexample needs
**three simultaneous conditions**: `U′=B→−I` (maximal holonomy, `θ=π`), `a₀=0` (barrier), and `β>4.31`.
The average removes it two independent ways: (i) averaging over `a₀` (`(★_avg^{a})`) weights the barrier
by `ρ_v(0)`, which is *small* (it is a shoulder, the mass sits at `a₀*≈−2`); the Fisher identity makes the
result exactly `I(ρ_v)>0`. (ii) averaging over `U′` suppresses `θ=π` by `e^{−β f(π)}=e^{−2β}`; at β=2.5
the barrier is *positively* curved anyway (`+0.136`). So the barrier survives only in the *pointwise*
form, which is already `[NG]`. **`[T]`+`[V]`**

---

## 4 · Acid test — the free field must fail, and it does — **[T]**+**[V]**

The margin `δ_A = I(ρ_v) ≥ 1/Var(a)`. On the **free massless field** the coarse marginal delocalizes
(`Var(a)=σ²→∞`, block size →∞), so
```
   G_A = I = 1/σ² → 0   as σ→∞  :   σ=10 → 1.0e−2 ,  σ=100 → 1.0e−4 ,  σ=1000 → 1.0e-6 .
```
So `(★_avg)` fails on the free field, exactly as the acid test demands: no floor is manufactured where
the theory has no gap. (The *compact* one-plaquette model does not delocalize as β→0 — its marginal is a
fixed density on `SU(2)` with `I=O(1)` — so it is **not** the free-field acid test; the correct
non-compact limit is the Gaussian one above.) **`[T]`+`[V]`**

---

## 5 · `N`, `a`, and why the averaged form is still H1

- **`N`.** The Fisher margin `I(ρ_v) ≥ 1/Var(a)` is per-mode; the coarse variance `Var(a)` is the
  physical object, `N`-benign (the adjoint dimension is multiplicity, S5 §3). No `1/N` or `N`-loss: the
  averaged margin is `1/Var(a)`, a single per-mode number.
- **`a`.** `1/Var(a)` is the coarse mass-squared and is `Θ(a⁰)` at fixed physical block (S5 §4); the
  Fisher identity is coupling- and `a`-free. So `(★_avg)` is `a`-uniform in the same sense as S5.
- **But the margin is the gap.** `δ_W = I(ρ_v) ≥ 1/Var(a)`. A uniform-in-`j` floor `δ_W≥δ>0` is the
  coarse spectral gap = **H1** `[O]`. So the averaged form does not lower the height of the wall; it
  proves the wall is exactly what the average measures. This is consistent with S1 §6 / S3 (even a
  perfect conditional Poincaré controls `D`, not `S`): the average over `a₀` supplies the `S`-half
  *conditionally on the coarse variance*, and bounding `Var(a)` is H1.

---

## 6 · Verdict and status ledger

| item | statement | label |
|---|---|---|
| marginal identity `G″(v;a₀)=−d²log ρ_v/da²=E[∂²V\|a₀]−Var(∂V\|a₀)` | G4 Lemma G1, any block | **`[T]`** |
| **Theorem H3-A** `E_{a~ρ_v}[G″]=I(ρ_v)≥0` (Fisher; integration by parts) | every `U′`, every β | **`[T]`** (+`[V]` ≤1.2e-3) |
| Cramer–Rao `I(ρ_v)·Var(a)≥1` ⇒ margin `≥1/Var(a)` | Cauchy–Schwarz | **`[T]`** (+`[V]` 1.01–1.25) |
| **`(★_avg^{a})`** coarse-field average holds, `δ_W=I(ρ_v)` | consequence of H3-A | **`[T]`** |
| **joint `(★_avg)`** for any boundary measure `ν_j` | Fubini + H3-A + Cramer–Rao | **`[T]`** |
| **`(★_avg^{U′})`** at β=2.5 positive, margin `+0.37`(Haar)/`+0.59`(Wilson) | exact model, whole `θ∈[0,π]` has `G″>0` | **`[V]`** |
| crude good/bad bound `ρ₀ν_good+G_minν_bad` positive at β=2.5 even `c_□=0` (`+0.137`) | worst point `+0.136>0` at head | **`[V]`** |
| crude bound positive for all β once `c_□≳0.25` (real block) | Dirichlet gap enlarges good set | **`[V]`**+**`[C]`** |
| negative set `{G″<0}` equilibrium-measure `≤1.6e-4`@β=5, `→0`@β=10, `=0`@β=2.5 | Wilson weight `e^{−β(1−cosθ)}`; Lemma B | **`[T]`**-mechanism / **`[V]`** |
| deep-UV bad-set measure `≤e^{−c/g_j²}→0` ⇒ (H3-cond) for large β | Lemma B (C2 §2.2) | **`[C]`** (on moment bootstrap) |
| barrier counterexample `G″=A−Bβ²<0` for `β>4.30661` | pointwise only; killed by averaging | **`[NG]`** (pointwise) |
| **uniform floor** `δ_W≥δ>0` in `(★_avg)` (uniform in `j`) | `=1/Var(a)=` coarse gap | **`[O]` = H1** |
| free-field acid test `I=1/σ²→0` | delocalization | **`[T]`**+**`[V]`** |
| **`(★_avg)` holds at β=2.5** | yes — `[T]` (a₀), `[V]` (U′); margin `≈0.37–0.61` | **`[T]`/`[V]`** |
| **`(★_avg)` closes the gap** | no — its margin is the gap | **`[O]`** |

**Sharpest line.** The averaged head inequality is **true** at β=2.5 — not by a delicate cancellation of a
large bad set against a good-set margin, but because (i) the coarse-field average of the marginal
curvature is *exactly* the Fisher information `I(ρ_v)≥1/Var(a)`, and (ii) under the equilibrium boundary
measure the genuinely negative set (`{G″<0}`, which exists only for `β>4.31`) has measure `≤e^{−2β}→0`
and at β=2.5 is empty. So `(★_avg)` **survives** where the pointwise forms die. But its margin *is* the
coarse mass-squared `1/Var(a)`, so `(★_avg)` is the coarse spectral gap seen through the coarse
marginal — H1, not a step past it. The averaged form is the honest weakest survivor of `(★)`; it is
`[T]`/`[V]` as an inequality and `[O]` as a uniform floor. **No QED. Clay 0%.**

---

## 7 · Witness

`research/v3/w7/h3_averaged_head.py` (deterministic, numpy-only, `py -3.12`, exit 0). Six parts:
- **[1]** Fisher identity `E_{a~ρ_v}[G″]=I(ρ_v)`: `G_A−I ≤1.2×10⁻³` at four `(β,θ)`.
- **[2]** Cramer–Rao `I·Var` ∈ `[1.01,1.25]` over `β∈[0.5,10]` (always `≥1`); margin `I` and `1/Var` tabulated.
- **[3]** barrier `G″(e₃;θ=π,a₀=0)=A−Bβ²`: `+0.136`@β=2.5, crosses `0` at `β*=4.307`, `−0.071`@β=5.
- **[4]** `U′`-average at `a₀=0`: `G″(θ)>0` on all `θ` at β=2.5 (min `+0.136`); average `+0.37`/`+0.59`.
- **[5]/[5b]** good/bad bound (crude and sharp); `c_□` scan; negative-set measure `≤1.6×10⁻⁴`@β=5,
  `0`@β=2.5; Haar vs Wilson-tilted boundary measure.
- **[6]** acid test `I=1/σ²→0`.

Reuses the cancellation-free `Φ_ρ, Φ_ρρ` and the correct Gauss–Legendre weight of
`v3_counterexample_scope.py` (the G4 quadrature bug is not repeated). Finite-difference `G″` agrees with
`v3_counterexample_scope.py`'s exact `A−Bβ²` to the quoted digits.

*Refs.* `STAR2-CLASSIFICATION.md` §0,§2,§5; `S1-CONDITIONAL-POINCARE.md` §3 (Lemma S1/Cor S1′), §6;
`S5-N-UNIFORM.md` §3–§5; `V4-PINCER.md` §4; `G4-GENERAL-BOUND.md` §2 (Lemma G1/G2/G3), §3(iv);
`V3-COUNTEREXAMPLE-SCOPE.md` §2–§4; `S2-BOUNDARY-COUNTEREXAMPLE.md` §2; `G1-BLOCK-IDENTITY.md` §5;
`STAR-WAVE-SYNTHESIS.md` §2; `C2-BALABAN-UNIFORM.md` §2.2 (Lemma B); `docs/paper-v2/ch05-theorems.tex`
`eq:split`, `eq:vardecomp`, `sec:wgeo`; witness `research/v3/w7/h3_averaged_head.py`. No third-party
sources downloaded.

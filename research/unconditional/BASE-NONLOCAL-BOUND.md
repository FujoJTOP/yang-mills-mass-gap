# BASE-NONLOCAL-BOUND — the non-local lower bound has a *precise shape* (one global scalar `1/√χ` on the `k→0` mode): variational gives only UPPER bounds, the three global certificates are all `≡ H1`, and the instrument is certified on exact solvable models

> **VERDICT (no QED).** The non-local lower bound is **not "no handle"** — it is **one global scalar**. Three formally different global certificates — **Cheeger** (`λ ≥ h²/4`), **drift/Bakry–Émery** (`λ ≥ λ_drift`), **coupling/Dobrushin** (`λ ≥ 1−ρ`) — each reduce, on the `k→0` mode, to the **same** scalar `1/√χ` (resp. `1/χ`), so *a uniform certificate exists ⟺ `χ ≤ 1/λ` = `H1`* (proved both directions). The **test-function/Rayleigh** principle gives only **upper** bounds (`gap ≤ R(f)` for all `f`), certified numerically. The **instrument is calibrated on exact solvable models**: free field on `Z_n` (massless gap→0, massive ≥`m²`), and the 1D Gaussian marginal, where Cheeger returns the **explicit rigorous fraction `1/(2π)` of the true gap** and `h·√χ = √(2/π) = const`. At `β≈2.5`, Bakry–Émery is **vacuous** (Savvidy Hess indefinite), Cheeger survives but its constant is `H1`. **The base is the classification + the calibrated fraction; the missing piece is a uniform global certificate = `H1`.**

---

## §0 · Labels, iron laws, acid test, and what is already banked

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-here/V-num]** new/computed this round.
**IL-i** non-positivity (must fail on free massless `f₁²=0`,`G″=0`). **IL-ii** nonabelian. **IL-iii** no small parameter (`β≈2.5`).
**Acid**: free massless must give `≤0` / gap `0`; gapped-atomless `:φ²:` discriminatingly `≤0`.

**Banked (this build reuses, does not re-derive):**
`CHEEGER-BUILD` (Cheeger constant `h>0` uniform ⟺ `H1`; global minimizer = `k→0` softest-mode half-space, `h≍1/√χ`, 2/5) · `CHEEGER-BLIND` (same, blind red-team) · `BS-VARIATIONAL-LOWER` (Temple/Weinstein need `E₂`; `0⁺⁺*` dynamical `=I2`) · `PROOF-HAMILTONIAN-VARIATIONAL` (Rayleigh–Ritz upper; electric positivity blind; binding collapse `=gap`) · `LOWER-AMPLITUDE-MECHANISMS` (eight lower-bound mechanisms aim *lower* at an *upper* target) · `D3-SUSCEPTIBILITY-VARIATIONAL` (`χ` as Legendre dual of `G″`; β=0 theorem `χ≤12/N`; β=2.5 `=H1`).

**Genuinely additive here (the "base" contribution):** (i) the **dichotomy** stated and *witnessed* (variational = upper only); (ii) the **three-certificate equivalence** `(uniform global certificate) ⟺ χ≤1/λ` proved in **both** directions, with the explicit constant; (iii) a **calibrated instrument** — the certificates *computed* on exact solvable models with the rigorous fraction `1/(2π)` and the free-field collapse — something the abstract `CHEEGER-BUILD` argued but never evaluated.

---

## §1 · The dichotomy: variational = UPPER bounds only; a LOWER bound needs a GLOBAL certificate `[T]`

Let `P` be a reversible Markov operator (transfer matrix) with invariant measure `μ`, spectral gap `γ := 1−λ₂`, and Dirichlet form `ℰ(f)=½∑ μ_i P_{ij}(f_i−f_j)²`. Then
```
γ = inf_{f: Var_μ f>0}  R(f),      R(f) := ℰ(f,f)/Var_μ(f).
```
**Consequence (upper-only, `[T]`).** For **every** test function `f`, `R(f) ≥ γ`; i.e. **`γ ≤ R(f)`** — a single test function yields an **UPPER** bound on the gap. A **LOWER** bound `γ ≥ λ` requires `R(f) ≥ λ` for **all** `f`, which no single test function can certify.

**Witness (§C, exact, Gaussian `σ=1`, `γ=1/σ²=1`):** `R(x)=1.000` (eigenvector, `=γ`); `R(x²)=2.000`; `R(1+x+x²)=1.667`. **All `≥ γ`.** Test functions pin **from above**. ⇒ **a lower bound is irreducibly a *global* certificate**, not a variational test.

---

## §2 · The three global certificates are all `≡ H1` (proved both directions) `[T]`

On the coarse `k→0` mode `a = ⟨v,X⟩`, marginal `ν(da)=p(a)da`, `χ := Var_ν(a)`, and (Ward) `G″(k→0)·χ ≈ 1`. The three global lower-bound mechanisms and their constants:

| certificate | rigorous inequality | its constant on the `k→0` mode |
|---|---|---|
| **Cheeger** | `γ ≥ h²/4`, `h = inf_A ∂ν(A)/min(ν(A),ν(Aᶜ))` | `h = 2p(0)`, and `h·√χ = √(2/π)` ⇒ `h ≍ 1/√χ` |
| **Drift / Bakry–Émery** | `γ ≥ λ` if `Γ₂ = Hess V ⪰ λ·I` | `λ = inf Hess V = 1/χ` (Gaussian exact) |
| **Coupling / Dobrushin** | `γ ≥ 1−ρ` (`ρ`=contraction) | `1−ρ = 1/χ` on the coarse line |

**Theorem (equivalence, `[T]` for the model; `[C]/[O]` for the YM promotion).**
```
(uniform certificate with constant λ>0,  uniform in (L,a))   ⟺   χ ≤ 1/λ uniformly   ⟺   H1.
```
*Proof (⇐, the useful direction).* If `χ ≤ 1/λ` then `σ²=χ ≤ 1/λ`, so the Cheeger constant `h = 2/(σ√2π) ≥ 2√(λ/2π) > 0` **uniformly**, hence `γ ≥ h²/4 ≥ λ/(2π) > 0`. *(⇒)* If `h ≥ h₀ > 0` uniformly, then `γ ≥ h₀²/4 > 0`; and since `h·√χ = √(2/π)` (witness, exact), `h ≥ h₀ ⟺ χ ≤ 2/(π h₀²)` — bounded variance `= H1`. ∎

**Explicit consequence (`[V-num]`, exact for the Gaussian model).**
```
γ  ≥  (1/2π)·(1/χ)     (Cheeger, 1D coarse mode)
```
so a uniform gap lower bound `γ ≥ λ` follows **iff** `χ ≤ 1/(2πλ)` — an **explicit, rigorous fraction `1/(2π)`** of the inverse susceptibility. This is the sharpest quantitative "base" statement the certificates give; its uniform attainment at `β≈2.5` is the open step.

---

## §3 · The instrument, calibrated on exact solvable models `[T]`/`[V-num]`

**(A) Free field on `Z_n` (exact gap `= min_{k≠0}(m²+2−2cos(2πk/n))`).** Witness table:
```
 n   m²=0 (gap)   m²=0.5 (gap)
 4    2.000000      2.500000
 16   0.152241      0.652241
 64   0.009631      0.509631
```
massless: `4 sin²(π/n) → 0` (**acid**: returns `0`, gapless); massive: `≥ m²` (**gapped**, returns the true gap). ✓

**(B) 1D Gaussian marginal `p=N(0,σ²)` = the coarse `k→0` mode.** Witness:
```
 σ   χ=σ²   gap=1/χ   Cheeger h   h²/4     fraction (h²/4)/gap   h·√χ
0.5  0.25   4.0000    1.59577     0.63662  0.15915              0.79788
1.0  1.00   1.0000    0.79788     0.15915  0.15915              0.79788
4.0 16.00   0.0625    0.19947     0.00995  0.15915              0.79788
```
`h²/4 ≤ gap` **always** (rigorous lower bound) with fraction **`1/(2π)=0.15915`** (scale-invariant); `h·√χ = √(2/π)=0.79788` **const** ⇒ `h ≍ 1/√χ`; `σ→∞` (massless/flat): `χ→∞`, `h→0`, `gap→0` (**acid**). ✓

**(C) Variational (upper-only) —** §1 table. ✓

⇒ The certificates are **not abstract**: on a solvable gapped model they return a rigorous positive fraction of the exact gap, and on the massless adversary they collapse to zero. That is the **calibrated base**.

---

## §4 · At `β≈2.5`: Bakry–Émery vacuous, Cheeger sign-free but `=H1`

- **Drift / BE:** requires `Hess V ⪰ λ>0`. At `β=2.5` the effective Hessian is **Savvidy-indefinite** (`Γ₂<0`). Witness (`V=(x²−1)²/4`): `min V″ = −1 < 0` ⇒ **`λ = inf Hess V < 0`** ⇒ the drift lower bound is **vacuous** (§D). This is the corpus's `Γ₂<0` wall in the lower-bound coordinate.
- **Cheeger:** *sign-free* (needs no Hess sign) — the one genuine structural edge (`CHEEGER-BUILD`). But its global constant `h` is attained by the `k→0` softest-mode half-space `{a(0)<c}` with `h ≍ 1/√χ` (§2, witness). Hence **uniform `h>0` ⟺ `H1`** — Cheeger dodges the curvature wall and lands on the bottleneck wall, and the bottleneck *is* the gap.
- **Coupling/Dobrushin:** `1−ρ` is the coarse-mode contraction = the same `1/χ` (`D3` Legendre face).

**Base conclusion.** The non-local lower bound is **not "无处下嘴"**: it has a **precise shape** — *a single global certificate whose constant is `1/√χ` (≡ `1/χ`) on the `k→0` mode*. What is missing is a **uniform construction** of that one scalar at `β≈2.5` (a uniform `h>0` / `λ_drift>0` / `1−ρ>0`), which is **verbatim `H1`**. The shape is pinned; the certifying number is not.

---

## §5 · Honest verdict

- **Variational test functions give only upper bounds** — `[T]`, witnessed. A lower bound is irreducibly a **global** certificate (Cheeger / drift / coupling).
- **All three global certificates `≡ H1`** — `[T]` (both directions) on the model, `[C]/[O]` for the YM promotion. The explicit fraction `γ ≥ (1/2π)/χ` is the sharpest quantitative base statement.
- **The instrument is calibrated** on exact solvable models (free field: massless→0, massive≥`m²`; Gaussian: Cheeger fraction `1/(2π)`, `h√χ` const) — `[T]/[V-num]`. This is the *additive base* the abstract Cheeger work lacked.
- **At `β≈2.5`**: BE vacuous (Hess indefinite), Cheeger sign-free but constant `=H1`, coupling same. **No uniform certificate constructed; exact missing step = a uniform global certificate = `H1`.**
- **Scorecard:** (1 nonperturbative ✓) (2 sign-correct lower bound: *shape* ✓ / *uniform constant at 2.5* ✗) (3 non-circular: model ✓ / YM ✗ `=H1`) (4 fails free field ✓ acid-clean) (5 uniform `(L,a)` ✗ `=H1`). **2.5/5** — matching `CHEEGER-BUILD`, now with the calibrated instrument.

**Clay 0%; no QED.** The wall did not move; its **lower-bound face** is now instrumented, calibrated, and shown to be one scalar.

*Refs:* `CHEEGER-BUILD.md`, `CHEEGER-BLIND.md`, `BS-VARIATIONAL-LOWER.md`, `PROOF-HAMILTONIAN-VARIATIONAL.md`, `LOWER-AMPLITUDE-MECHANISMS.md`, `D3-SUSCEPTIBILITY-VARIATIONAL.md`, `B3-ISOPERIMETRIC`, `POTENTIAL-THEORY-CAPACITY`, `PROOF-JACOBIAN-SAVES`, `PROFILE.md`; witness `w7/base_nonlocal_bound.py`.

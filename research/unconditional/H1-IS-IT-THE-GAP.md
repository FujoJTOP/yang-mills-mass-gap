# H1-IS-IT-THE-GAP — the H1 curvature `G″(k=0)` is the **inverse zero-momentum susceptibility**, NOT the pole mass; the exact squeeze, and why `Var(a)≤1/λ` is *necessary but blind to a thin atom*

> 2026-11-23 · Conceptual settlement. Question under test: **is the H1 quantity *the* mass gap, or is it a susceptibility / condensate-fluctuation curvature — a different object?** Answers below are done as explicit Källén–Lehmann algebra, with the free-field normalisation fixed for the massive and massless cases *and* — the discriminating one the campaign never checked — the free **composite** `:φ²:` (gapped, `f₁²=0`). Iron rules: 真定理 / 条件定理 / no-go / 缺失步 are labelled at each step; nothing is called proved that is not.

**判词 (verdict, up front).**
1. **`G″(k=0)` is literally `1/χ = 1/∫s⁻¹dν = 1/C̃(0)`, the inverse zero-momentum susceptibility** (equivalently the curvature of the zero-momentum effective potential, `V_eff″(0)`, the Legendre `Γ″(0)`). It is **the reciprocal of** the spectral number `∫s⁻¹dν` — that integral *is* `χ`, the denomination, not `G″`. It equals the Källén–Lehmann pole mass `m₁²=inf supp ν` **if and only if ν is a single atom** (Gaussian marginal, `f₁²=M₀`). In every other case they are different functionals.
2. **The exact squeeze is `f₁²/χ ≤ m₁² ≤ M₀/χ`** (`M₀:=∫dν` = total spectral weight, `f₁²=ν({m₁²})` = edge atom). **H1 (`χ≤1/λ`, i.e. `G″=1/χ≥λ`) yields `m₁² ≥ f₁²·λ` — a gap only together with the residue bound `f₁² ≥ r₀>0`**; it yields *nothing* on the mass from above that helps, because the correct harmonic-mean relation is `1/χ ≥ m₁²` (harmonic mean ≥ min), so `G″` is an **upper** bound on `m₁²`, not a lower one.
3. **VERDICT: NO — the H1 quantity is a susceptibility, not the gap.** It is the inverse zero-momentum susceptibility; the gap is the **KL edge/pole datum** `m₁=inf supp ν` (lightest `0⁺⁺` glueball), equivalently the **edge atom `f₁²=ν({m₁²})` bounded below with isolation** — i.e. `H* = I1∧I2`, a two-point pole/decay-rate statement. H1 is a **curvature face** of that one object and is **necessary (up to `M₀`) but not sufficient**; the equivalence `H1 ⟺ gap` is the corpus's own **conjecture** (G2-ASSEMBLE-100 §7: "a structural identification … a conjecture, not a theorem"), not a theorem. The `NO` carries a precise **conditional** refinement in §3: `H1 + {f₁²≥r₀, isolation I2, the a²-ratio, a gauge-invariant carrier}` = `H*` = the gap.

---

## §1 · Definitions, fixed once and unambiguously (the overloaded `a` is the whole confusion)

A gauge-invariant scalar interpolator `O` (a composite — `tr F²`, a plaquette/Wilson-loop density, a glueball operator). All objects are transverse, `0⁺⁺`, colour-singlet. Euclidean.

| symbol | definition | meaning / status |
|---|---|---|
| `C(τ)` | `∑_x ⟨O(x,τ)O(0,0)⟩_c` | zero-**spatial**-momentum connected correlator |
| `ν` | KL measure: `C(τ)=∫_0^∞ e^{−mτ}dν(m)`, `dν≥0`, `supp ν⊂[m₁,∞)` | spectral measure (on **mass** `m`, equivalently on `s=m²`; I keep `m` and write `m₁²` for the edge of `ν` on `s`) — 真定理 (RP/OS) |
| `m₁` | `inf supp ν` | lightest mass (lightest `0⁺⁺`) — **the gap** |
| `f₁²` | `ν({m₁²})=|⟨0|O|G⟩|²` | edge **atom** (residue) |
| `M₀` | `∫_0^∞ dν` | total spectral weight (residues sum; UV-regularised on the lattice) |
| `χ` | `∫_0^∞ s⁻¹ dν(s) = C̃(0) = ∑_x∫dτ⟨O(x,τ)O(0)⟩_c` | **zero-momentum susceptibility** (the denomination) |
| `C̃(k)` | `∑_x e^{ikx}⟨O(x)O(0)⟩_c` | momentum-space connected propagator; `C̃(0)=χ` |
| `G″(k)` | `V_eff″(k)` = curvature of `−log p_k` | marginal / effective-potential curvature — 真定理 = `1/Var(a_k)` (Gaussian-exact; `G″·Var=1−R`, `R`=2–5% measured, H1-DIRECT-SPLIT (I3)) |

**The two normalisations of the soft mode `a_k` — the crux the campaign never disambiguated.** Write `Õ(k)=(1/V)∑_x e^{ikx}O(x)` (the *mean*) and `a^{(tot)}_k=(1/√V)∑_x e^{ikx}O(x)=√V Õ(k)` (the *total* amplitude). Then, exactly (CLT/translation invariance on the finite torus):
```
Var(a^{(tot)}_k) = C̃(k)            → C̃(0) = χ     at k→0   (NO V-suppression)
Var(Õ(k))        = C̃(k)/V          → χ/V → 0      at k→0   (self-averaging)
```
They differ by the factor `V` and **only** by it. *Which* `a` one writes decides whether H1 reads `χ<∞` (meaningful) or `χ/V→0` (vacuous). The lattice uses **both** without saying so: the gluon proxy `a=⟨v,log U⟩` is the *total*-normalised field mode (Var `≈0.27`, V-robust, W7); the gauge-invariant density zero-mode `q̄=(1/V)∑q` is the *mean* (Var `∝V⁻¹`, H1-NUMERIC-RIGHT §2). I call these the **propagator branch** and the **susceptibility branch**.

---

## §2 · The Källén–Lehmann computation and the squeeze (algebra done, not quoted)

### §2.1 `G″(k=0)` in terms of `ν` — all four candidate readings resolved

Start from `C(τ)=∫e^{−mτ}dν(m)`. The zero-momentum susceptibility (Laplace at `s=0`):
```
χ = ∫_0^∞ C(τ)dτ = ∫_0^∞ ∫_0^∞ e^{−mτ}dν(m) dτ = ∫_0^∞ (1/m)dν(m) = ∫ s⁻¹ dν(s).      (真定理)
```
So **`χ = ∫s⁻¹dν = C̃(0)`** — the susceptibility *is* the zero-momentum propagator *is* that integral. On the *total*-normalised mode, `Var(a^{(tot)}_0)=C̃(0)=χ`; on the *mean* mode, `Var(Õ(0))=χ/V`. The marginal curvature of `−log p_0` is, by the constrained-free-energy identity (H1-DIRECT-SPLIT (I2), 真定理),
```
G″(0) = V_eff″(0) = E[H_00|0] − Var(f|0) = 1/Var(a^{(tot)}_0) = 1/χ = 1/C̃(0) = 1/∫s⁻¹dν.   (Gaussian-exact)
```
`V_eff″` is the curvature of the **zero-momentum effective potential**; its Legendre dual `Γ″(0)` (curvature of the effective action w.r.t. a source) satisfies `Γ″(0)·C̃(0)=1`, so `Γ″(0)=1/χ` too. Therefore, term by term against the task's four candidates:

| candidate | is it `G″(k=0)`? | what it actually is |
|---|---|---|
| `1/χ` | **YES** | — |
| `∫ s⁻¹ dν` | **NO** | this is `χ` itself (= `C̃(0)`), the *denominator* |
| `1/C̃(0)` | **YES** | identical to `1/χ` |
| curvature of zero-momentum effective potential `V_eff″(0)` | **YES** | identical to `1/χ` |

**`G″(k=0) = 1/χ = 1/∫s⁻¹dν = 1/C̃(0) = V_eff″(0) = Γ″(0)`** — all one object; the lone impostor is "`∫s⁻¹dν`," which is `χ`, and "the inverse" of it is `G″`. **真定理**: `G″` *is* the inverse zero-momentum susceptibility. It is `m₁²` only if `ν` is a single atom.

### §2.2 The squeeze (the exact inequality, proved on `ν`)

`ν` on `s=m²≥m₁²`. Decompose `ν=f₁²δ_{m₁²}+ν_rest`, `ν_rest` total `M₀−f₁²`.

- **Lower bound on `m₁²` from the atom.** `χ=∫s⁻¹dν ≥ ∫_{atom}s⁻¹dν = f₁²/m₁²` (all terms `≥0`). Hence
  ```
  m₁²  ≥  f₁²/χ.                                                              (L)
  ```
- **Upper bound on `m₁²` from the total.** On `supp ν`, `s≥m₁² ⇒ s⁻¹≤1/m₁²`, so `χ=∫s⁻¹dν ≤ m₁⁻²M₀`. Hence
  ```
  m₁²  ≤  M₀/χ.                                                              (U)
  ```

**Squeeze (真定理):**  `f₁²/χ ≤ m₁² ≤ M₀/χ.`  The two sides are the edge-atom contribution (L) and the residue-weighted harmonic mean (U). They coincide (`f₁²=M₀`) iff `ν` is a single atom; otherwise they bracket `m₁²` by a factor `M₀/f₁²`.

### §2.3 Free-field normalisation — done for **three** cases (the third is the discriminating one)

**(i) Free MASSIVE, elementary `O=φ`.** `ν=δ_{m²}` (single atom, canonical residue 1).
`M₀=1`, `f₁²=1`, `m₁²=m²`, `χ=1/m²`, `G″=1/χ=m²`. Squeeze tight: `1/(1/m²)=m²=m²=m²`. **All readings agree — degenerate.**

**(ii) Free MASSLESS, elementary `O=φ`.** `ν=δ_0` (atom at `s=0`).
`m₁²=0`, `f₁²=1`, `χ=1/0=∞`, `G″=1/χ=0=k̂²|_{k=0}`. **`G″(0)=0`, `χ=∞`, `m₁²=0`** — again degenerate (every quantity vanishes at the gap edge), but the **acid test is clean**: a *finite* `G″` or a *finite* `χ` would be vacuous; here `χ=∞` correctly diverges, so **free massless violates H1 (`χ≤1/λ` fails, `∞≰1/λ`)** — H1 is non-vacuous and does not falsely certify the gapless theory. 真定理.

**(iii) Free MASSIVE, composite `O=:φ²:` — GAPPED, `f₁²=0`.** `ν` is a *cut* on `[4m²,∞)` (`:φ²:` couples only to two-particle states; `⟨0|:φ²:|1⟩=0` by `φ→−φ`). So
```
m₁² = 4m²  (threshold),   f₁² = 0,   χ_{:φ²:} = ∫_{4m²}^∞ s⁻¹ρ(s)ds  finite (O(1/16π²), m>0),
G″ = 1/χ_{:φ²:}  finite,  and  G″ ≠ 4m²  in general.
```
This is **the case that breaks the degeneracy**: a genuinely gapped theory (`m₁=2m>0`) with `f₁²=0`, whose H1 quantity `G″=1/χ_{:φ²:}` is **unrelated** to its gap `4m²`. Squeeze reads `0 ≤ 4m² ≤ ∞` — *no constraint*. **`H1` on `:φ²:` (`χ_{:φ²:}<∞`) is satisfied and says nothing about whether or for how large a gap exists.** It also exhibits the thin-atom mechanism (the edge atom is *absent*, weight 0). This is why the elementary free-field checks (i)–(ii) — both single-atom — are **insufficient** to test H1 and *why the campaign's "acid test" on the gluon never engaged the real question*.

---

## §3 · Sufficiency vs necessity — the crux, and why the naive "harmonic-mean" argument is *wrong*

The tension the task names: a naive "harmonic-mean ≤ min" would give `1/χ ≤ m₁²`, making H1 sufficient. **It is false on two independent counts.**

**(a) Direction of the harmonic mean is reversed.** `χ=∫s⁻¹dν`; the residue-weighted **harmonic mean** of the masses-squared is `H=(∫dν)/(∫s⁻¹dν)=M₀/χ`. The harmonic mean satisfies **`H ≥ min`** (not `H ≤ min`; check `{1,4}`, equal weights: `H=1.6>1=min`). Hence
```
M₀/χ = H ≥ m₁²    ⇒    1/χ ≥ m₁²/M₀     (so, for a normalised measure M₀=1:  1/χ = G″ ≥ m₁²).
```
So `G″` bounds `m₁²` from **above** (`m₁² ≤ G″`). The naive premise "`1/χ ≤ m₁²`" is exactly backwards: it would require `harmonic mean ≤ min`, which is the **root of the error**. Consequence: **H1 (`G″≥λ`) is logically compatible with an arbitrarily small `m₁`** (`m₁² ≤ G″` gives nothing below); it is *not* sufficient. REDUCTION-AUDIT is correct; the naive route is a sign error in a textbook inequality.

**(b) Even the correct lower bound carries the atom weight `f₁²`.** The genuine *lower* bound is `m₁² ≥ f₁²/χ` (L). Combine with H1 `1/χ≥λ`:
```
H1  (χ≤1/λ  ⟺  G″=1/χ ≥ λ)   ⇒   m₁² ≥ f₁²/χ ≥ f₁²·λ.        (真定理, one line)
```
- **Sufficient IF `f₁² ≥ r₀ > 0`:** `m₁² ≥ r₀·λ > 0` — a gap. (With the constant rescaled; the sharp statement is `m₁² ≥ r₀/χ`.)
- **Blind if `f₁² = 0` (thin/absent atom):** the bound collapses to `m₁² ≥ 0`. **A near-massless thin atom passes H1.** Concretely `ν = wδ_ε + (1−w)δ_S` with `w/ε ≤ 1/λ` gives `χ ≤ 1/λ + (1−w)/S ≤ 1/λ'` (H1 holds) while `m₁² = ε → 0`. This is the exact counterexample in REDUCTION-AUDIT (b); it is decisive.

**(c) Necessity *does* hold (the one real theorem in the chain).** If there is a gap `m₁² ≥ λ`, then `s ≥ λ` on `supp ν`, so `χ = ∫s⁻¹dν ≤ M₀/λ`. Thus
```
gap (m₁²≥λ)  ⇒  χ ≤ M₀/λ  ( = H1 with 1/λ' = M₀/λ ).        (真定理)
```
So H1 is a **necessary** condition (up to the `M₀` normalisation) — *exactly* what a susceptibility that is dominated by the light modes should be — but the converse fails per (a)+(b).

**Resolution.** The chain `gap ⟺ H1` is a one-directional theorem plus a conjecture. H1 is the **inverse susceptibility** face; sufficiency requires the **edge/residue** datum `f₁²≥r₀` (and isolation `I2`) — which *is* the gap's unprotected content, merely relocated. That is precisely the corpus's `H* = I1∧I2`, which RESIDUE-CONDITIONAL/G2-RESIDUE-ATOM prove is **strictly stronger than bare `m>0`**; and PROOF-F1SQ shows the gap alone does **not** manufacture the atom (massive free `:φ²:` is gapped, `f₁²=0`). **H1 ⟺ `H*` is the structural identification the corpus itself labels a conjecture** (G2-ASSEMBLE-100 §7).

---

## §4 · Physical identity of `G″(k=0)` — which one, and are they the same object?

| | **pole mass** `m₁²=inf supp ν` | **H1 curvature** `G″=1/χ=1/∫s⁻¹dν` |
|---|---|---|
| what it is | position of the lightest KL **pole** (lightest glueball) | **inverse zero-momentum susceptibility** (effective-potential/action curvature) |
| how it reads off `ν` | the edge of the support | a residue-weighted **inverse-harmonic** moment, dominated by light states |
| gauge-invariant? | yes (if `O` invariant) | yes (if `O` invariant) — but see §5: for composites `Var∝1/V` makes it V-divergent |
| free massive elementary `φ` | `m²` | `m²` (agree — single atom) |
| free massless `φ` | `0` | `0` (agree — both vanish) |
| free massive **composite** `:φ²:` | `4m²` | `1/χ_{:φ²:}` **≠** `4m²` (disagree; gapped, `f₁²=0`) |
| general | `m₁²` ∈ `[f₁²/χ, M₀/χ]` | `1/χ` = reciprocal harmonic mean |

**Determined conclusively:** `G″(k=0)` is **the inverse zero-momentum susceptibility**, not the pole mass. The two are the **same object only in the single-pole/Gaussian limit** (`f₁²=M₀`, no continuum/tail) — which is *exactly the elementary free field*, the case the campaign used. They are **different in every composite/multi-state case**, and every gauge-invariant operator in Yang–Mills is composite. The elementary-field free check is the degenerate case and therefore *cannot* vindicate "H1 = gap"; the free massive `:φ²:` check does the opposite — it falsifies it.

---

## §5 · Reconciliation with the numerics — it confirms the wrong-object verdict, *not* a patch

Two measured facts (H1-NUMERIC-RIGHT, W7, NUMERICS-AUDIT), understood through §1's normalisation split:

1. **Gauge-invariant composite `tr F²` density: `Var(q̄)=χ/V ∝ V^{−1.0}`, `G″_dens=1/Var∝V→+∞`.** This is the **susceptibility branch under the mean normalisation** (§1): the physically finite object is `χ=∫s⁻¹dν≈0.04`, V-uniform; the mode variance is `χ/V→0`; `G″_dens=V/χ` **diverges** — the *opposite* sign of scaling to a bounded inverse propagator. So `Var(a)≤1/λ` holds for **every** `λ` at large `V`: **H1 is vacuous on this object.** This is *not* a patchable normalisation slip — it is the CLT/self-averaging of any composite zero mode. **It realizes option (i): the k=0 curvature of a gauge-invariant composite is a susceptibility and carries no gap content.**

2. **Gluon `a=⟨v,log U⟩`: `G″≈0.26±0.11`, V-robust, but `G″/free<1` at every β (0.18–0.63), `LB_rig<0` everywhere.** This is the **propagator branch** (`a^{(tot)}`-normalised field mode): `Var=1/(k̂²+m²)`, `G″=k̂²+m²`, the **inverse gluon propagator** — the *right type* of object. But it (a) is the gauge-**dependent** transverse gluon (`1⁻` adjoint), not the `0⁺⁺`; (b) sits **below** free at every accessible `k`, so it fails the acid test (no IR mass **excess** — `G″` is a Haar/`W_geo` number the free field *beats*). **It realizes option (ii)-but-broken:** right type, but needs a *gauge-invariant* carrier with `Var↛0`, which H1-NUMERIC-RIGHT §3 proves **does not exist linearly** (no linear gauge-invariant function of the connection exists).

**Neither branch certifies the gap, and they are disjoint** (H1-NUMERIC-RIGHT §3 structural theorem: `a_{k→0}` is a composite ⇒ `Var→1/V` ⇒ vacuous, *or* the fundamental field ⇒ gauge-dependent; no third option).

**On the dimension correction `G″_a=a²m_phys²`** (PROOF-SKELETON §0, H1-DIRECT-SPLIT §2): it is real and correctly forces the **ratio** `G″_a/a²≥Δ₀²` instead of a constant `λ`. But it is strictly a fix *within* the **propagator branch** (it is the Wilson/field Hessian's engineering dimension). It is **necessary but not sufficient**: it removes the `a`-engineering error and does nothing about (a) gauge-invariance (the composite/field dichotomy) or (b) the composite `V/χ` divergence. Patching the ratio on the **gluon** still leaves a gauge-dependent number that the free field beats; patching it on the **density** still leaves `V/χ→∞`. **The dimension correction is not the missing ingredient; the object-type is.**

**Bottom line:** the numerics did not find H1 "slightly wrong-scaled" — they found it **wrong-type** on the gauge-invariant side and **gauge-dependent and acid-test-failing** on the field side. That is the signature of a susceptibility, not the gap.

---

## §6 · The correct object (since H1 is not it)

The Clay gap is a **gauge-invariant** statement. The only gauge-invariant carriers are composites (§5, structural theorem). For a composite, the k=0 curvature is a susceptibility (`1/χ` or `V/χ`), **not** the pole mass. Therefore the gap lives in the **two-point POLE / edge residue**, exactly:

> **X := the KL edge datum of the gauge-invariant `0⁺⁺` channel: `m₁ = inf supp ν > 0` (lightest glueball), equivalently the edge atom `f₁² = ν({m₁²})` bounded below with isolation `Λ−m₁²≥δ>0`.** This is `H* = I1∧I2` (RESIDUE-CONDITIONAL; PROOF-ISOLATION L1; GAP-G2-edge residue–level-gap lemma). It is a **decay-rate / pole-residue** datum of the gauge-invariant 2-point correlator — read numerically as the exponential decay `C(τ)∼e^{−m₁τ}` (the standard glueball mass), **proved** as `f₁²≥r₀` plus isolation, `(L,a)`-uniform.

- It is **not** a k=0 curvature sign; it is a **pole location + residue**.
- It is **one unprotected number** (RESIDUE-CONDITIONAL): trace anomaly/Polchinski `c`, marginal curvature `G″`, and LGP `γ` are its *vocabulary*, not its proof; and the identification with `G″` is the conjecture, not the theorem.
- The free massive `:φ²:` (§2.3-iii) is the acid model: gapped (`m₁=2m`), `f₁²=0` — showing why the *residue/isolation* is the irreducible content, strictly above bare `m>0`.

**Scorecard (5 points, honest).**

| test | H1 quantity (`1/χ`) | correct object `X` (edge pole/residue `f₁²`) |
|---|---|---|
| Nonperturbative | PASS (exact KL algebra) | PASS (OS/RP exact; residue nonperturbative) |
| Same as pole mass? | **FAIL** — inverse susceptibility; `=m₁²` iff single atom | **YES by definition** (`m₁=inf supp ν`) |
| Fails free massless (acid) | PASS (`χ=∞`) | PASS (`f₁²=0`, gap 0) |
| Fails free massive `:φ²:` (the real acid) | **FAIL** — `1/χ_{:φ²:}` finite, gap `4m²`, unrelated | PASS (`f₁²=0` there; `>0` only when bound) |
| Gauge-invariant, `Var↛0` carrier exists | **FAIL** — composite⇒`1/V` vacuous; field⇒gauge-dep | PASS (2-pt pole is gauge-invariant, no `Var` issue) |
| Sufficient for gap | **FAIL** alone (thin atom); needs `f₁²≥r₀` | **YES** (it *is* the gap, strictly stronger) |
| Uniform `(L,a)` | open (=H1) | open (=H*, residue unprotected) |

**Sharpest line.** `G″(k=0)=1/χ=1/∫s⁻¹dν=1/C̃(0)` is the **inverse zero-momentum susceptibility** — the reciprocal of a residue-weighted harmonic mean, which by "harmonic mean ≥ min" *upper*-bounds `m₁²`; it equals the pole mass only in the single-atom Gaussian case the campaign mistook for generic. The squeeze is `f₁²/χ ≤ m₁² ≤ M₀/χ`, so H1 buys exactly `m₁² ≥ f₁²λ` and **nothing more** — a gap the instant `f₁²≥r₀`, silence for a thin atom. **H1 is the curvature face of the gap, not the gap; the gap is the edge atom `f₁²=ν({m₁²})` at `m₁=inf supp ν` — one unprotected residue (H*), the corpus's own `H* = I1∧I2`, strictly stronger than `m>0`.**

---

## Exact answer / exact stuck step

```
QUESTION:  is the H1 quantity (Var(a_{k→0})≤1/λ  ⟺  G″=1/χ≥λ  ⟺  m_eff²(k→0)≥λ) the mass gap?
ANSWER:    NO — it is the inverse zero-momentum susceptibility   G″(k=0)=1/χ=1/∫s⁻¹dν=1/C̃(0)=V_eff″(0).
           It equals the pole mass m₁²=inf supp ν  IFF ν is a single atom (f₁²=M₀, Gaussian).
           Exact squeeze:   f₁²/χ  ≤  m₁²  ≤  M₀/χ     (M₀=∫dν, f₁²=ν({m₁²})).
           Necessity (真定理):  gap m₁²≥λ  ⇒  χ ≤ M₀/λ   (H1 up to M₀).
           Sufficiency (真定理, one line):  H1  ⇒  m₁² ≥ f₁²/χ ≥ f₁²·λ  — gap ONLY if f₁²≥r₀>0.
           Thin-atom witness: ν=wδ_ε+(1−w)δ_S, w/ε≤1/λ ⇒ χ≤1/λ but m₁²=ε→0.  (H1 blind.)
           Naive "1/χ≤m₁²" is FALSE: harmonic mean ≥ min ⇒ 1/χ ≥ m₁²/M₀ (G″ is an UPPER bound on m₁²).
GAP IS:    X = KL edge datum  m₁=inf supp ν  (lightest 0⁺⁺ pole), equivalently the edge atom
               f₁²=ν({m₁²}) bounded below with isolation Λ−m₁²≥δ>0   ( = H* = I1∧I2 ).
CONDITIONAL form:  H1 ⟺ gap  holds only as  H1 + {f₁²≥r₀, I2 isolation, a²-ratio, gauge-invariant carrier}
               = H* ; the identification H1≡H* is a CONJECTURE (G2-ASSEMBLE-100 §7), not a theorem.
NUMERICS:  gauge-invariant composite → Var(q̄)=χ/V→0, G″=V/χ→∞ : H1 VACUOUS (wrong object, option (i));
           gluon → G″=k̂²+m² right type but gauge-dependent and G″/free<1 (acid FAIL, option (ii) broken).
           The dimension fix G″_a/a² is a prop-branch patch — necessary, NOT sufficient; object-type is the wall.
STATUS:    H1 = the gap  →  necessary-not-sufficient (真定理 necessary + 条件定理 sufficient), a susceptibility.
           Clay gap = 0%.  The unprotected content is one edge residue f₁²≥r₀ (H*), strictly above m>0.
```

**What this settles (and does not).** It settles the *conceptual aim*: the campaign is aimed at the **curvature face** (`1/χ`) of the right object, and `Var(a)≤1/λ`/`G″≥λ` is a **susceptibility bound**, not the gap — because `G″=1/χ` is an inverse harmonic moment that (i) upper-bounds the mass and (ii) is blind to a thin edge atom. The correct target, in this vocabulary, is the **edge atom `f₁²=ν({m₁²})` bound below with isolation** (the KL pole mass `m₁=inf supp ν`), i.e. a **two-point pole/decay-rate** datum of a gauge-invariant `0⁺⁺` correlator — *not* a k=0 curvature. It does **not** prove the gap: the residue `f₁²≥r₀` is unprotected and uniform-in-`(L,a)` equally open (RESIDUE-CONDITIONAL; GAP-G2-edge: edge control is one rung above the gap). The re-aim is real, the wall is unmoved.

---

*Refs.* **Källén–Lehmann / OS:** Källén 1952; Lehmann 1954; Osterwalder–Schrader 1973–75; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Seiler, LNP **159** (1982). **Internal — definitions & the squeeze:** `PROOF-F1SQ.md` (`ν`, `M_k`, residue–level-gap lemma, massive `:φ²:` `f₁²=0`); `GAP-G2-edge.md` (moment-ratio *upper* bound; `f₁²≥(ΛM₀−M₁)/(Λ−m₁²)`; edge ≡ one rung above gap); `RESIDUE-CONDITIONAL.md` (`H*=I1∧I2` is the gap; one unprotected residue); `PROOF-ISOLATION.md` (L1: atom `f₁²≥r₀` + isolation `Λ−m₁²≥δ`); `G2-RESIDUE-ATOM.md` §9 / `REDUCTION-AUDIT.md` (atom strictly stronger than gap; `χ=∫s⁻¹dν` is a susceptibility; thin-atom counterexample `wδ_ε+(1−w)δ_S`); `D3-SUSCEPTIBILITY-VARIATIONAL.md` (`χ=Var(a)=∫dρ`, `G″·χ≈1`, free `χ=∞`); `H1-DIRECT-SPLIT.md` (I2 `G″=E[H|0]−Var(f|0)=V_eff″`; (I3) `G″·Var=1` iff Gaussian, 2–5%; dimension ratio `G″_a/a²`). **Internal — the equivalence is a conjecture:** `G2-ASSEMBLE-100.md` §7 ("structural identification … a conjecture, not a theorem"), §1 MASTER, §3 rows 1–24 (`H1=H*`, no third class). **Internal — numerics:** `H1-NUMERIC-RIGHT.md` (§1 gluon `G″/free<1`, `LB_rig<0`; §2 density `Var∝V⁻¹`, `G″=V/χ→∞`; §3 structural theorem no 3rd option; §4 glueball 2-pt unresolved at L=8); `NUMERICS-AUDIT.md` (finite-k free beats interacting; LB<0; not a continuum); `K0-OVERLAP-AUDIT.md` (measured mode is `1⁻` gluon, not `0⁺⁺`); `W7-Gpp0.md`, `W7-RERUN.md` (`G″≈0.27`, `1/Var≈0.27`, `W_geo`-carried, free at `k=π/4` is `0.586`). **Dimension:** `docs/PROOF-SKELETON.md` §0 (`G″_a=a²m_phys²`; ratio target `G″_a/a²≥Δ₀²`). **External:** Brascamp–Lieb, Adv. Math. **22** (1976) 366; Prékopa 1971/73.

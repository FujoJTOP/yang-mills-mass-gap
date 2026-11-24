# H1-DIRECT-SPLIT — the exact variance decomposition on the one inequality H1: `G″ = E[f²] − D`, the sharpest form, and the single missing estimate

> 2026-11-23 · Head-on attack on **H1** (no new "species"): untangle the banked identities, write `H1` in its sharpest exact-equivalent forms, and pin the one estimate that separates "D/S≈2% measured" from a theorem.
>
> **判词 (verdict, up front).** The identities are all **真定理** and they close the algebra exactly; they do **not** close the sign. The honest chain is: **(I1)** `E[H]=E[f²]` (IBP/SD, every measure, every β); **(I2)** `G″(0)=E[H|0]−Var(f|0)` (constrained free energy, exact); **(I5)** `G″=S−D` (SCBL, exact). Two corrections to the task's framing are load-bearing: (i) **`G″·Var(a)=1` is *not* an identity** — the exact Ward statement is `E[a f]=1`; `G″·Var(a)=1` holds **iff the marginal is Gaussian** (measured to 2–20%, `REDUCTION-AUDIT` (C)); (ii) the "`S=E[f²]`" of the task and the **SCBL `S`** (`=E[h−⟨K⁻¹c,c⟩|0]`, the Gaussian–Schur curvature, `=1/Var(a)` only for Gaussian) are **different objects** — the exact identity is `G″=E[f²]−D̃` with `D̃=(E[H]−E[H|0])+Var(f|0)`, whereas SCBL is `G″=S−D` with `D=Var(f|0)−E[⟨K⁻¹c,c⟩|0]`. **H1 in sharpest form** is one conditioned-variance inequality: `Var_μ(⟨v,∇V⟩|a=0) ≤ E_μ[H_vv|a=0] − λ` (⟺ `D ≤ S−λ` ⟺ `G″_a/a² ≥ Δ₀²`), `(L,a)`-uniform. **It is not `D ≤ (1−ε)S` alone** — that form silently needs `S ≥ λ` too, which is the *same* coarse-variance bound, and the `E[f²]` version is **dimensionally unsafe** (the measured `E[f²]≈2` is UV/equipartition-pinned and does not scale as `a²`). **Nothing here proves H1**; the exact missing step and its obstruction (`K` Savvidy-indefinite ⇒ BL reverses ⇒ `D` is a difference of two `O(2)` numbers with no structural sign) are named in §3. Clay: **0%**.

---

## §0 · Three "S"s and one "D_"s — the objects, exactly

Flat chart `X = log U ∈ su(N)` per link, Haar `dU = J dX = e^{−Φ}dX`, `Φ = −log J ⪰` (convex, `Φ″ ⪰ N/12`, Lemma B). `V = S_W + Φ`, `μ ∝ e^{−V}dX`. Fix a **coarse transverse plane-wave mode** `v` (unit in the algebra norm), amplitude `a = ⟨v,X⟩`. The three derived scalars:

| symbol | definition | meaning |
|---|---|---|
| `f` | `∂_n V := ⟨v,∇V⟩ = G + Φ′` | **force** (`G=∂_n S_W` dynamical, `Φ′=∂_nΦ` Jacobian) |
| `H_vv` | `∂_n f = ⟨v, Hess V · v⟩ = H^S_vv + Φ″` | **curvature** (Hessian along `v`) |
| `H` | `Hess V` | full Hessian; `H^S = −D²+2ad F_B` (Savvidy, §4) |

**Two distinct "defects", do not conflate:**
- **`D̃`** (the `E[f²]`-form): `D̃ := E[f²] − G″`. Identity `G″=E[f²]−D̃` is exact *given* (I1),(I2); it is the "unconditioning gap + force variance."
- **`D`** (the SCBL form): `D := Var(f|0) − E[⟨K⁻¹c,c⟩|0]`, `K=Π₀ Hess V Π₀`, `c=Π₀ Hess V v`. Identity `G″=S−D` is exact; this `D` is the **Brascamp–Lieb defect of the conditioned fiber** — the only one with a clean sign theory (`D≤0 ⇔ K⪰0`).

The task's slogan "`S=E[f²]`" upgrades (I1) to *define* the SCBL `S`; that is wrong (`S≠E[f²]` in general). The correct exact map is below.

---

## §1 · The identities, each with its exact status

**(I1) IBP / Schwinger–Dyson / Stein / "Stratonovich".** On the compact link group (boundary term vanishes by compactness/periodicity), for any `C¹` `g`:
```
E_μ[∂_n g] = E_μ[g f].                                  (IBP)
```
Two specialisations:
- `g=a` (`∂_n a=1`): **`E_μ[a f] = 1`** — the exact Ward identity. (真定理, every measure, every β.)
- `g=f` (`∂_n f=H_vv`): **`E_μ[H_vv] = E_μ[f²]`** — the "SD identity." (真定理.)

*This* is the identity the task calls "Stratonovich/self-dual." It is a genuine theorem for **every** measure (YM and free alike). **Acid test:** free massless `k→0` has `H_vv=0` so `E[H]=E[f²]=0` — finite — while `Var(a)=∞`; the identity is **gap-blind** (FORGE5-CONE-VARIANCE, SUB1).

**(I2) Constrained free energy.** With `p(a)` the marginal and `V_eff(a)=−log p(a)`:
```
V_eff′(a)=E_μ[f|a],   V_eff″(a)=E_μ[H_vv|a]−Var_μ(f|a),
G″(0):=V_eff″(0)=E_μ[H_vv|a=0]−Var_μ(f|a=0).            (真定理, exact, every β)
```
The numerics compute the right-hand side *directly* (`w7/gpp_core.py`: `G″=E[H_vv|a]−Var(G_v|a)`; `G_v`=`vᵀ∇S`), so the measured `G″>0` is this object.

**(I3) Ward ⇒ Gaussian relation — a correction.** From (I1) with `g=a`, `E[a m(a)]=1` where `m(a)=E[f|a]=V_eff′(a)`. At the reflection-symmetric point `m(0)=0`, `m′(0)=G″`:
```
E[a m(a)] = G″·Var(a) + R,   R = Σ_{k≥2} m^{(2k+1)}(0)/(2k+1)! · E[a^{2k+2}].
⇒  G″·Var(a) = 1 − R,   R = 0 ⟺ marginal Gaussian.
```
**So "`G″·Var(a)=1`" is NOT an exact identity** — it is the Gaussian statement (`V_eff` quadratic). The exact identity is `E[a m(a)]=1`; the deviation `R` (higher cumulants) is measured at 2–20% (`G″` vs `1/Var`: W7-RERUN `0.281` vs `0.272`; SU3 `0.425` vs `0.404`). [Matches `REDUCTION-AUDIT` (C): definition-or-false, not a theorem.]

**(I4) Jacobian / `W` split.** `f=G+Φ′`, `H_vv=H^S_vv+Φ″`; differentiating twice,
```
G″(0) = Δ_dyn + W_geo + X,
  Δ_dyn = E[H^S_vv|0] − Var(G|0),   W_geo = E[Φ″|0] − Var(Φ′|0),   X = −2Cov(G,Φ′|0).
```
(identity, all conditioned). **Proven pieces:** `E[Φ″|0] ≥ (N/12)‖v‖²` under *any* measure (Lemma B); `W_geo ≥ 0` (Prékopa, geometric marginal). **Uncontrolled:** signs of `Δ_dyn`, `X`.

**(I5) Savvidy-corrected BL (SCBL).** With `Π₀=I−vvᵀ`, add/subtract the Schur form `⟨K⁻¹c,c⟩` (it cancels):
```
G″ = S − D,   S = E_μ[h_vv − ⟨K⁻¹c,c⟩ | 0],   D = Var_μ(f|0) − E_μ[⟨K⁻¹c,c⟩|0].
D ≤ 0 ⇔ K ⪰ 0 (BL on μ₀);  D = 0 ⇔ μ₀ Gaussian;  D > 0 ⇔ K indefinite (Savvidy).
S = 1/Var(a)  exactly for Gaussian V (Schur complement).                        (真定理)
```
`S` is the **Gaussian–Schur curvature**; `D` is the **fiber BL defect**. This is the only decomposition with a *sign theory* for the defect (`D`'s sign = `K`'s definiteness).

**(I5′) Relation to (I1)** — the task's "`G″=E[f²]−D`":
```
G″ = E[f²] − D̃,   D̃ := E[f²] − G″ = (E[H] − E[H|0]) + Var(f|0).
```
Exact, but `D̃` is the *conditioning gap + force variance*, **not** the SCBL `D`. For a Gaussian marginal `E[H]=E[H|0]` and `D̃=Var(f|0)`, but in general `D̃ ≠ D`. **Use `D` (SCBL) when you want a sign; use `D̃` only to state the inequality against `E[f²]`.**

**`D` explicitly in terms of `Var(∇S|0)` and `W_geo`.** Using `Var(f|0)=Var(G+Φ'|0)` with `G=∂_n S_W` (so `Var(G|0)=Var(⟨v,∇S_W⟩|0)`), `Var(Φ'|0)=E[Φ''|0]−W_geo`, and `X=−2Cov(G,Φ'|0)`:
```
D = Var(G|0) + Var(Φ'|0) + 2Cov(G,Φ'|0) − E[⟨K⁻¹c,c⟩|0]
  = Var(G|0) + [E[Φ''|0] − W_geo] − X − E[⟨K⁻¹c,c⟩|0].       (exact, all conditioned)
```
So `D` is the conditioned dynamical-force variance `Var(G|0)`, **plus** the Jacobian force variance `Var(Φ'|0)=E[Φ'']−W_geo`, **plus** their covariance `−X`, **minus** the fiber Schur reference `E[⟨K⁻¹c,c⟩|0]`. The only *provably signed* ingredient is `W_geo≥0` (so `Var(Φ'|0)≤E[Φ'']`, an upper bound that pushes `D` **up**); `Var(G|0)`, `Cov`, and the Schur term carry the uncontrolled competition. This is the precise sense in which **the sign lives in `Var(∇S|0)` against `E[⟨K⁻¹c,c⟩|0]`**, not in the provable `W_geo` margin.

**Chain verdict.** (I1)(I2)(I4)(I5) are 真定理 (exact identities). (I3)'s Gaussian form is conditional. The `E[f²]=S` upgrade is **false** in general. The gap is in every case the same inequality `G″≥λ`.

---

## §2 · H1 in its sharpest equivalent forms; the two target inequalities

Given (I2), **H1 is exactly one inequality** and every "form" below is an exact algebraic relabel of `G″(0)≥λ`:

| form | statement | what it isolates | status |
|---|---|---|---|
| **primitive** | `Var_μ(f\|a=0) ≤ E_μ[H_vv\|a=0] − λ` | conditioned force variance | **the target** |
| **SCBL** | `D(v) ≤ S(v) − λ` | fiber BL defect `D` vs Schur `S` | **the target** |
| **ratio (continuum)** | `G″_a/a² ≥ Δ₀²`, any `Δ₀>0` | removes the `a²` engineering | **the target** |
| `E[f²]`-form | `D̃ ≤ E[f²] − λ` | unconditioned force moment | relabel (cumbrous) |
| susceptibility | `a² Var_μ(a) ≤ 1/Δ₀²` | coarse variance (`≈1/G″`) | relabel (via I3, 2–20%) |

**Dimension correction (mandatory).** `G″_a = a²m_phys²` (dimensionless Wilson Hessian), `m_lat=√G″_a=a·Δ₀`, `Δ=lim m_lat/a`. So the **dimensionless, RG-correct** target is the ratio
```
G″_a ≥ a²Δ₀²   ⟺   a²·Var_μ(a) ≤ 1/Δ₀²,   uniform in (L,a).
```
"`G″_a ≥ λ` with `λ` an `a`-independent constant" is **dimensionally wrong**: it forces `m_phys=lim√λ/a=+∞` — the free-box artifact `4π²/L²`, not a particle (G1-CONTINUUM-SCALING, P1-L0). The task's `Var(a)≤1/λ` must be read as the ratio `a²Var(a)≤1/Δ₀²`.

**The two "defect-small" forms the task names — and why neither closes:**

1. **`D ≤ E[f²] − λ`** (⟺ `G″≥λ`): exact relabel. To use it you must upper-bound `D̃` against the *unconditioned* `E[f²]`. **Fatal:** `E[f²]≈E[H]≈2` at β=2.5 is **equipartition/UV-pinned** (W7: `E[H]=2.46/3.18`; `Var≈E[H]`; B4: `2.01→2.05` L4→L6) and does **not** scale as `a²`. A constant `G″≥λ` from it is the free-box artifact. **The `E[f²]`-form is dimension-unsafe without the ratio.**

2. **`D̃ ≤ (1−ε)E[f²]`** (⟺ `G″ = E[f²] − D̃ ≥ εE[f²]`): the `E[f²]`-defect-small form. Same dimension warning as (1): if `E[f²]` is UV-pinned (`≈2`, not `∝a²`), `G″≥εE[f²]` is the free-box constant again. **Safe only in the ratio `D̃_phys ≤ (1−ε)E_phys[f²]` with `a²`-scaled moments.**

3. **`D ≤ (1−ε)S`** (SCBL defect fraction `< 1`): *necessary-looking, insufficient.* Since `G″=S−D`, `D≤(1−ε)S ⇒ G″≥εS`. This yields a **uniform gap only if `S ≥ λ_S` uniformly** — and `S=1/Var(a)` (Gaussian) is *itself* the coarse variance bound. So closing needs the **pair** `{S≥λ_S, D≤S−λ}`; the first conjunct is the *same* H1 at a possibly larger constant, not a weaker one. There is **no factorization of H1** into "denominator large" + "defect small": both are the single conditioned-variance inequality. (This is the sharpest structural statement in §2.)

**Cleanest sharpest line for H1:**
```
H1  ⟺  Var_μ(⟨v,∇V⟩ | a=0)  ≤  E_μ[H_vv | a=0] − λ,   λ ≥ a²Δ₀²,  (L,a)-uniform,  v⊥orbit coarse transverse.
    ⟺  D(v) ≤ S(v) − λ   (SCBL).    Not "D ≤ (1−ε)S" alone (needs S≥λ), not "D≤E[f²]−λ" (dimension-unsafe).
```

---

## §3 · Why it is hard despite `D/S≈2%`; the exact missing estimate

**Three independent reasons the 2% does not close it.**

**(a) The 2% is not uniform — one point, off the scaling window.** `D/S≈2%` is measured at `β=2.5`, `k=π/2,π/3`, `ξ/a≲0.3`, a single lattice spacing, no `a`-scan (CONTINUUM-VARIANCE-DEFECT). `D/S` is RG-**covariant**, not invariant: a dimensionless function of the running coupling, `→0` at the Gaussian point, `→O(1)` at the gap-relevant IR mode `k→0`. Passing `a→0` kills the UV numerator and proves nothing about the IR.

**(b) `D` is a difference of two `O(2)` numbers, not a small object.** `D = Var(f|0) − E[⟨K⁻¹c,c⟩|0]`. Measured `Var(f|0)≈2` (B4), and BL holds to 2% *because* `E[⟨K⁻¹c,c⟩|0]≈Var(f|0)` (almost-cancellation). There is **no structural reason** `D≪S`: on the Savvidy cone `K` is indefinite, `‖K⁻¹‖` blows up at the negative-mode zero-crossing, and `Var(f|0)` can be `≈E[f²]`. So **`D` can in principle reach `E[f²]`** — the 2% is the numerics of BL *barely* holding at β=2.5, not a bound.

**(c) D is UV, S is IR — the ratio is decided on the denominator.** `D` is fiber/UV-dominated (90–97% of Savvidy negatives in the fluctuation block, W6-BE) and AF-suppressed as `a→0`. But `S=1/Var(a)` is the **IR** coarse inverse propagator; `S≥a²Δ₀²` *is* `m_eff²(k→0)≥λ` *is* H1. A UV upper bound on `D` cannot bound `D/S` without an independent IR lower bound on `S`. Trace anomaly `θ=(β/2g)¼trF²` fixes `D`'s *scale* (`Z_O²`) but not its sign or size.

**The exact missing estimate (pinned, one line per face):**
```
∃Δ₀>0:  ∀ U′, ∀ coarse transverse v ⊥ orbit,  uniformly in (L,a):
   (◆)   Var_μ(⟨v,∇V⟩ | a=0)  ≤  E_μ[H_vv | a=0] − a²Δ₀²
   ⟺     D(v) ≤ S(v) − a²Δ₀²
   ⟺     G″_a(v)/a² ≥ Δ₀²      (the mass gap in curvature language).

Obstruction to proving (◆):
   the RHS needs  E_μ[⟨K⁻¹c,c⟩|0] ≤ Var(f|0) − a²Δ₀²,
   i.e. the conditioned fiber is BL-good to within the gap — but K is Savvidy-indefinite
   (D>0), ‖K⁻¹‖ is unbounded at the negative-mode edge (det K=0 ⟺ transverse massless mode),
   and no volume-uniform upper bound on Var(f|0) exists
   (only the β² Efron–Stein magnitude, largest on the cone — wrong sign).
```

**Why the acid test passes / the free field is excluded:** free massless `k=0`: `H_vv=0`, `Var(f|0)=0`, so (◆) reads `0 ≤ 0 − a²Δ₀²` — **unsatisfiable** for `Δ₀>0`; correctly rejected (`Var(a)=∞`, `G″=0`, `f₁²=0`). The inequality is non-vacuous.

---

## §4 · What can be proven now vs assumed; the `S`-split check; scorecard

### §4.1 Truth ledger

| statement | status | source |
|---|---|---|
| `E[H]=E[f²]` (IBP/SD), every measure/β | **真定理** | (I1); FORGE5-CONE-VARIANCE |
| `G″(0)=E[H\|0]−Var(f\|0)=V_eff″(0)` | **真定理** | (I2); k0-CONSTRAINT-EFFPOT |
| `G″=S−D` (SCBL); `D≤0⇔K⪰0`, `D>0⇔K` indef | **真定理** | (I5); INVENT-VARIANCE-INEQ |
| `W_geo≥0` (Prékopa), geometric marginal | **真定理** | (I4); PROOF-JACOBIAN-SAVES |
| `E[Φ″\|0]≥(N/12)‖v‖²`, any measure | **真定理** | Lemma B; A5 |
| `Tr H^S = Tr(−D²) ≥ 0` (trace-half-positivity) | **真定理** | PROOF-SAVVIDY-EXACT-SIGN |
| β=0: `G″≥N/12` (k→0, L-uniform) | **真定理** | k0-CONSTRAINT-EFFPOT |
| β<1/48: coarse Poincaré `Var_ν(φ)≤C_blk²λ_μ⁻¹∫\|∇φ\|²` | **真定理 (SZZ)** | PROOF-COARSE-POINCARE (c) |
| `G″·Var(a)=1` | **conditional** (Gaussian marginal; exact is `E[af]=1`) | (I3); REDUCTION-AUDIT (C) |
| `S=E[f²]` | **false in general** (only Gaussian/free) | §0, §4.2 |
| `D/S≈2%` at β=2.5 | **numerical, one point, not uniform** | W7-RERUN; CONTINUUM-VARIANCE-DEFECT |
| **H1 / (◆) at β=2.5, `(L,a)`-uniform** | **open (= the gap)** | this file; G2-ASSEMBLE-100 row 18 |

### §4.2 The `S = E[−D²]+E[2ad F_B]+E[Φ″]` check (task item 4)

The **operator** split is the Savvidy one: `H^S = Hess S_W = −D² + 2ad F_B` (PROOF-SAVVIDY-EXACT-SIGN §1, spectrum `λ_{n,s}=k_∥²+(2n+1)g|B|+2s gB`, `λ_min=−g|B|`). Hence the exact decomposition
```
E[H] = E[−D²] + E[2ad F_B] + E[Φ″],        and by (I1):   E[f²] = E[−D²] + E[2ad F_B] + E[Φ″].
```
**Check of the hypothesized `E[−D²] = E[f²]`:** this would need `E[2ad F_B] + E[Φ″] = 0`. It holds **only for the free field** (`F_B=0`, `Φ″=0`, so `E[f²]=E[−D²]=k̂²→0` at `k→0` — the acid test). In general it is **FALSE**: `E[f²]=E[H]` is the *whole* sum, not the kinetic piece.

**Which piece is uncontrolled — and why the split does not help.** The three pieces: `E[−D²]≥0` (kinetic, → free `k̂²`, dies at `k→0`); `E[Φ″]≥N/12` (Jacobian, provable positive); `E[2ad F_B]` (**Savvidy-indefinite**, `λ_min=−g|B|`, unbounded below, pointwise). One might hope `E[2adF_B]` is "the" hard piece — but the split **cannot localize H1**, for an exact reason: (I1) says the three pieces sum to `E[f²]≥0` *unconditionally* (a Schwarz fact, content-free). The gap is the **conditioned** object `G″=E[H|0]−Var(f|0)`; conditioning replaces the unconditioned `≥0` slack by the slice variance, and `K=Π₀ Hess V Π₀` is indefinite **whatever** operator decomposition you apply to `H`. So **no decomposition of `S` factors H1**; the uncontrolled object is invariant — it is `Var_μ(⟨v,∇S⟩|0)` (equivalently the BL defect `D>0` from the indefinite fiber `K`). The sharper decomposition exists as an *identity* but is **not a sharper sub-inequality**.

**The single sharpest sub-inequality remains (◆)** (`Var(f|0)≤E[H_vv|0]−λ`). If one insists on a two-piece split, the sharpest is `{S≥λ_S}` ∧ `{D≤S−λ}` — but both conjuncts are the same coarse-variance bound (§2).

### §4.3 Five-point scorecard (honest)

1. **Nonperturbative — PASS.** All identities are exact IBP/constrained-free-energy/Prékopa; no loop expansion; every β.
2. **Sign-controlling — FAIL.** The sign is `G″=E[H|0]−Var(f|0)≥λ`; the two terms are opposing and neither is bounded in the needed direction (`E[H|0]` Savvidy-indefinite, `Var(f|0)` un-upper-bounded).
3. **Non-circular — FAIL at closing.** (◆) = SCBL `D≤S−λ` = H1 verbatim; every decomposition is a relabel with zero slack.
4. **Non-SUSY + fails free field — PASS (cleanest).** Free massless `k→0` gives `0≤−a²Δ₀²` unsatisfiable; `G″=0`, `f₁²=0`. The Haar term `E[Φ″]≥N/12` is the discriminator.
5. **Uniform `m_eff²(k→0)≥a²Δ₀²` — FAIL at 2.5.** Only β=0 (`N/12`) and β<1/48 (SZZ) are uniform; at β=2.5 only `≥0` (uncertified, LB `−0.10`).

---

## Sharpest line

`H1` is exactly `G″(0)=E[f²]−D̃=E[H|0]−Var(f|0)=S−D≥a²Δ₀²`; every identity in the chain is a 真定理, and every one of them is a **relabel of the single conditioned-variance inequality** `Var_μ(⟨v,∇S⟩|a=0) ≤ E_μ[H_vv|a=0] − a²Δ₀²` — equivalently `D ≤ S − a²Δ₀²`. The `S`-split `E[f²]=E[−D²]+E[2ad F_B]+E[Φ″]` is exact but **does not localize** the gap, because the unconditioned sum is `≥0` for free reasons (Schwarz) and the sign lives in the **conditioned** fiber Hessian `K`, which is Savvidy-indefinite no matter how `H` is split. **The wall is one number:** the BL defect `D>0` of the conditioned fiber — a difference `Var(f|0)−E[⟨K⁻¹c,c⟩|0]` of two `O(2)` quantities with no structural sign, small (2%) only because BL *barely* holds at one non-scaling-window point.

## Exact stuck step

```
(◆)  Var_μ(⟨v,∇V⟩ | a=0)  ≤  E_μ[H_vv | a=0] − a²Δ₀²       [uniform (L,a,U′), v⊥orbit]
     ⟺ D ≤ S − a²Δ₀²  ⟺  G″_a/a² ≥ Δ₀²  ⟺  a²Var_μ(a) ≤ 1/Δ₀².
Needs:  E_μ[⟨K⁻¹c,c⟩ | 0] ≤ Var(f|0) − a²Δ₀²  with K⪰0 up to the gap.
Blocker: K = Π₀ Hess V Π₀  is Savvidy-indefinite (D>0), ‖K⁻¹‖→∞ at the negative-mode edge
         (det K=0 = a transverse massless mode = the wall); no volume-uniform upper bound on
         Var(f|0) exists.  This is H1 renamed, not weakened.  Clay = 0%.
```

**What would move it:** an `(L,a)`-uniform proof that the conditioned fiber measure `μ(·|a=0)` is log-concave **to within `a²Δ₀²`** — a per-scale, measure-weighted positivity of `Π₀ Hess V Π₀` (not pointwise: pointwise fails, `−1.12β`). Equivalently a defect estimate `D ≤ S − a²Δ₀²` whose two sides *both* carry the `a²` scaling (the ratio, not the constant). Neither is in the corpus; both are the single named species.

---

*Refs.* IBP/SD/Stein: FORGE5-CONE-VARIANCE (`E[H]=E[f²]`), SUB1-COV-BOUND (unconditioned, gap-blind). Ward/Gaussian: ABC-NEW-INEQ §1 (`E[af]=1⟺G″Var=1`), REDUCTION-AUDIT (C). Constrained free energy: k0-CONSTRAINT-EFFPOT, PROOF-COARSE-POINCARE (§"object"), `w7/gpp_core.py`. SCBL/`D`: INVENT-VARIANCE-INEQ (`G″=S−D`, `D≤0⇔K⪰0`). `W_geo`/Prékopa/Lemma B: PROOF-JACOBIAN-SAVES, A5-RIGOROUS-NUMERICS. Savvidy operator: PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, `Tr H=Tr(−D²)≥0`). SD/BPS: PROOF-SELFDUAL-BACKGROUND. Measurements: W7-Gpp0 (+0.51), W7-RERUN (identity, β=0 closure 0.248/0.246, `1/Var`), B4-VARIANCE-MAG (equipartition `2.01→2.05`), NUMERICS-AUDIT (LB `−0.10`). Scaling/ratio: G1-CONTINUUM-SCALING (`G″_a=a²m²`), CONTINUUM-VARIANCE-DEFECT (`D/S` RG-covariant), P1-L0 (free `4π²/L²`). Reduction: PROOF-IR-H1, OTHER-H1-POINCARE, FORGE5-SL-IR-SINGLE, G2-ASSEMBLE-100 (row 18 = H1). External: Brascamp–Lieb, JFA **22** (1976) 366; Prékopa 1971; Holley–Stroock, Ann. Probab. **15** (1987) 1541; Efron–Stein 1981; Savvidy, PLB **71** (1977) 133; Leutwyler, NPB **179** (1981) 129; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737.

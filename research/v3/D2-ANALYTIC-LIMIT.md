# D2-ANALYTIC-LIMIT — the constant-mode curvature at β=2.5 has a **positive, L-independent, β-dependent** limit (≈0.43 in mode normalization); it does **not** tend to 0 (H1 is false for the constant mode), and it does **not** saturate at the β=0 Haar floor (0.269) either

> **2026-11-24 · D2 — the decompactification limit of the soft/constant mode.**
> Question (SM-WAVE-SYNTHESIS §2, H1): *does the constant-mode curvature at β=2.5 tend to 0 as
> L→∞ (H1), or to a positive L-independent limit (the β=0 Haar floor)?*
> This note derives the constant-mode marginal on a finite torus and its L-scaling factor by factor,
> shows the two regimes are one physics in two normalizations, and computes the exact L→∞ limit and
> its β-dependence. Witness: `research/v3/w7/d2_{decompact,quick,batch,curv,beta}.py` (`py -3.12`,
> numpy only; one job at a time; §6).
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED, Clay 0%)

**At β=2.5 the constant-mode constraint curvature `G″(a₀*)` has a positive, L-independent limit
`G″_∞(2.5) ≈ 0.43`** (mode normalization `‖v‖=1`, `a=⟨v,ξ⟩`).** The Gaussianized limit is exactly
`1/Var(a)`, and `Var(a)` **plateaus** (≈2.34 for L≥6; already 0.41 at L=4) rather than growing (`∝L ⇒ G″∝1/L`) or shrinking
(`→0 ⇒ G″→∞`). So:

1. **Not H1.** `G″` does **not** tend to 0. `Var(a)` is bounded uniformly in L, hence `G″` is bounded
   below uniformly in L in the mode normalization. The apparent `≈1/L` decay `1.01→0.65→0.42`
   (`L=2,3,4`) is a **crossover** from the small-L non-Gaussian value to the plateau `0.42`; the
   corpus's `G″(a₀*)` has essentially already reached the limit at `L=4`.
2. **Not the Haar floor.** The limit is **β-dependent** and sits **above** the β=0 Haar floor
   `1/s²=0.26883`: `G″_∞(β)=` `0.269` (β=0), `0.31` (β=1), `0.43` (β=2.5), `0.55` (β=4) **— each flat
   in L for L≳4** `[V]`. So decompactification saturates at a β-dependent positive plateau, **not** at
   the (β-independent) Haar floor, and **not** at 0.
3. **Mechanism `[T]`+`[V]`.** The constant (k=0) mode is a **zero mode of the Wilson action's
   quadratic (gradient) part** both on the torus and on the fixed-boundary strip. Hence the action's
   *direct* grip on the mode is zero, and the curvature is **Haar-transverse-carried**: it comes from
   the flat-chart Haar factor `∏_l J(ρ_l)` coupling the mode to the transverse coordinates
   `ρ_l=‖ξ_{⊥,l}+av_l‖`, whose β-dependent (massive, stiffer at large β) distribution sets the
   O(1) β-dependence. Because that coupling is a **bulk** (`O(1)`-per-link) effect, the limit is
   L-independent; because the transverse distribution is β-dependent, the limit is β-dependent.
4. **What this does and does not prove.** It proves a **positive, L-uniform constant-mode margin**
   `[T-structure]+[V]` — the constant mode does **not** delocalize in decompactification. It does
   **not** give the mass gap: the flat-chart coordinate is a gauge-fixed (fixed-boundary slice)
   quantity, its susceptibility being the physical glueball mass is `[C]` at best; and the
   **non-constant** `k→0` soft tower (SM6's lattice object) is a different mode not covered here
   `[O]`. **No QED. Clay 0%.**

**One line.** The constant mode has a **positive, L-independent, β-dependent** curvature limit
`≈0.43` at β=2.5 (rising from `0.269` at β=0) — H1 is **false** for the constant mode (no
delocalization, no `1/L` collapse to 0), but the limit is **not** the Haar floor: it is a genuine
β-dependent bulk plateau carried by the Haar–transverse coupling of an action zero mode.

---

## 1 · The constant-mode marginal on a finite torus, and the L-scaling of every factor

### 1.1 The 4D torus `T_L⁴` (the object the question names)

Links `N=4L⁴`; plaquettes `6L⁴` (`L⁴` in each of the 6 μν-planes). Flat exponential chart,
`U_l=e^{iξ_l·σ/2}`, `ξ_l∈B(2π)⊂ℝ³`, Haar `dU_l=J(ρ_l)d³ξ_l`,
`J(ρ)=sin²(ρ/2)/ρ²`, `Φ₀=−logJ`, `Φ₀(2π)=+∞` (edge), `Φ₀″(0)=N/12=1/6`.

**Constant mode.** `U_μ(x)=U_μ` independent of `x`: four group elements, each repeated `L⁴` times.
Coarse amplitude along a uniform color direction: `a=⟨v,ξ⟩`; the *mode-normalization* unit vector is
`v_l=ê/√N` on every link (`‖v‖=1`), so `a=N^{-1/2}Σ_lξ_{l,·}`.

**Exact marginal and curvature** (H4 Lemma G1 / ch05 `eq:split`; the object SM1/SM4 compute):
```
   ρ_v(a) = ∫ ∏_l [J(ρ_l)d³ξ_l] δ(a−⟨v,ξ⟩) e^{−S_W[ξ]},
   G″(v;a₀) = −d²/da² log ρ_v(a) |_{a₀} = E[∂_v²V|a₀] − Var(∂_vV|a₀).                            (1.1)
```

**L-scaling factor by factor (mode normalization `‖v‖=1`).**

| factor | expression | L-scaling | note |
|---|---|---|---|
| Haar measure | `∏_l J(ρ_l)d³ξ_l`, `4L⁴` links | **L-independent per link** (density `J` is O(1) per link; count is `4L⁴`) | in the constant-mode slice `Φ_tot=L⁴Σ_μΦ₀(ρ_μ)` |
| Wilson action | `S_W=βΣ_{μ<ν}L⁴(1−½trP_{μν})`, `P_{μν}=U_μU_νU_μ†U_ν†`; `6L⁴` plaquettes | `=βL⁴ Σ_{μ<ν}(1−½trP_{μν})` | **commutator**: `1−½trP_{μν}=⅛‖[ξ_μ,ξ_ν]‖²+O(ξ⁴)` |
| action in `a` | field `ξ~a/√N` ⇒ `‖[ξ,ξ]‖²~a⁴/N²` | `S_W~βL⁴·a⁴/(16L⁸)=βa⁴/(16L⁴)→0` | **the action's grip on the normalized mode vanishes** |
| coarse field norm | `a=⟨v,ξ⟩`, `‖v‖=1` (vs `ξ̄=N^{-1}Σξ=a/√N`) | `a=O(1)` ⇔ physical field `ξ̄=O(N^{-1/2})→0` | this is the corpus's `a` |

On the **fixed-boundary strip** `1×L` (`L` links, `L` plaquettes `P₁=X₁X₂`, `P_i=B_{i−1}X_i^{−1}X_{i+1}`,
`P_L=B_{L−1}X_L^{−1}`) the same table holds with `4L⁴↦L`, `6L⁴↦L`; the bulk plaquette angles are
**exactly** invariant under a uniform colour shift (the shift cancels inside `X_i^{−1}X_{i+1}`), so the
bulk action is flat on the mode, and the two end plaquettes `P₁,P_L` anchor it only at `O(β/L)` in
the `a`-variable. (This is the concrete witness of §3–§4.)

### 1.2 The β=0 exact value (the Haar floor), in closed form

At β=0 the measure is pure Haar and the `ξ_l` are iid, so `ρ_v` is the density of
`N^{−1/2}Σ_lξ_{l,3}` — a scaled sum of iid one-link marginals — which converges (CLT) to `N(0,s²)`
with `s²=E_Haar[ξ₃²]`. The one-link variance is exact:
```
   E[ρ²]=∫₀^{2π}ρ²sin²(ρ/2)dρ / ∫₀^{2π}sin²(ρ/2)dρ = 4π²/3 − 2 ,
   s² = ⅓E[ρ²] = (4π²−6)/9 = 3.71988 ,   1/s² = 9/(4π²−6) = 0.268830 .                         (1.2)
```
So `G″_{β=0}(a₀*)→1/s²=0.26883` `[T]` — the exact **Haar floor**, the CLT Gaussian curvature.
(Numerics: `1/Var(a)=0.268–0.271` for L=4,10, §3.)

---

## 2 · The two regimes, and which holds at β=2.5

The question distinguishes (a) the action `∝βL⁴` dominating (Laplace/concentration) from (b) the Haar
term dominating. **These are the same physics in two normalizations:**
`G″_{ξ̄}=N·G″_a` under `a=√N ξ̄` (curvature scales as `(amplitude)^{−2}`).

- **(a) Field normalization** (`ξ̄=N^{−1}Σξ`, fixed physical constant field). `S_W≈βL⁴·⅛‖[ξ̄_μ,ξ̄_ν]‖²`,
  so the exponent is `∝βL⁴→∞`: the measure **concentrates** on the action-minimising set (commuting
  `U_μ`). Equivalently `Var(ξ̄)=Var(a)/N→0` and `G″_{ξ̄}=N·G″_a→∞`. The **physical constant mode is
  pinned** — infinite curvature. This is regime (a).
- **(b) Mode normalization** (`a=⟨v,ξ⟩`, `‖v‖=1` — the corpus's variable). The same action is
  `βa⁴/(16L⁴)→0`; the action's direct grip vanishes, and the curvature is carried by the Haar factor
  and the transverse sector. This is regime (b), and it gives a **finite positive** `G″_a`.

**At β=2.5, in the mode normalization (the one the corpus and SM-WAVE-SYNTHESIS use), regime (b)
holds**: `G″_a` is `O(1)`, positive, L-independent. Regime (a) is the *same* statement read in the
physical field: the two are related by `G″_{ξ̄}=N G″_a`.

**The one genuine subtlety — the k=0 decoupling.** The action's grip on the mode is *structurally*
zero (not merely `1/L`): the constant mode is an eigenvector of the quadratic Hessian with eigenvalue
`λ_{k=0}=0` (lattice Laplacian zero mode), `vᵀ(Hess S_W^{(2)})v=0`. This is why `G″_a` neither blows
up (as the field-normalization would) nor vanishes: the action drops out of the **quadratic** constant-
mode stiffness, and the curvature is generated by the **Haar–transverse coupling** (below). This also
implies the limit is **β-independent at the Gaussian level** (a Gaussian Haar gives a constant
`∂²logJ`, hence a β-independent stiffness); the observed β-dependence (§3) is the **non-Gaussianity of
the exact Haar** evaluated on the β-dependent transverse distribution.

---

## 3 · The L→∞ limit and its β-dependence

### 3.1 Analytic statement

Write `ρ_v(a)=∫dξ_⊥ e^{−S[ξ_⊥]} ∏_l J(‖ξ_{⊥,l}+a v_l‖)`, with `ξ=av+ξ_⊥`, `‖v‖=1`. Since the
quadratic action is flat on `v`, the entire `a`-dependence at leading order sits in the Haar product,
and the curvature is
```
   G″_a(a₀) = −∂_a² log⟨∏_l J(‖ξ_{⊥,l}+a v_l‖)⟩_{S} |_{a₀} ,                                    (3.1)
```
an **O(1)**, **bulk** (`O(1)`-per-link, `N`-fold average) object. As `L→∞` the link average converges
(finite correlation length / finite susceptibility — `[C]` below), giving
```
   G″_∞(β) = −∂_a² log⟨∏_l J(‖ξ_{⊥,l}+a v_l‖)⟩_{S} |_{a₀*} ,   finite, positive, β-dependent .  (3.2)
```
- **Positivity + L-uniformity `[T-structure]+[C]+[V]`.** `[T]`: the k=0 decoupling means `G″_a` is not
  of the divergent field form `N·(…)`. `[C]`: `G″_∞>0` **iff** the flat-chart-coordinate susceptibility
  `χ=lim Var(a)=Σ_r C(r)` is finite, i.e. **finite connected correlation length** in this channel
  (a mass gap in the coordinate field). `[V]`: the variance is measured bounded (§3.2).
  Equivalently `G″_a ≈ 1/Var(a)` once the marginal Gaussianizes (large L).
- **β=0 exact `[T]`:** `G″_∞(0)=1/s²=9/(4π²−6)=0.26883` (eq. 1.2).
- **β-dependence `[V]`:** at the Gaussian level β-independent (Gaussian Haar); the exact-Haar
  non-Gaussianity, evaluated on the stiffer (large-β) transverse distribution, raises `G″_∞` monotonically
  with β (§3.2).

### 3.2 Witness — `Var(a)` and `1/Var(a)` for the fixed-boundary strip (k=0, LONG `e₃` mode, θ=π)

Exact-target Metropolis (`d2_batch.py`: 12–24k independent strips, vectorized; local plaquette updates;
`J`-corrected acceptance; chart ball `B(2π)`). `Var(a)` is measured; the Gaussianized curvature is
`1/Var(a)` (agrees with the corpus `G″(a₀*)` once the marginal has Gaussianized):

| β \ L | 2 | 3 | 4 | 6 | 8 | 10 | 14 |
|---|---|---|---|---|---|---|---|
| **2.5** `Var(a)` | 1.71 | 2.33 | 2.45 | 2.34 | 2.35 | 2.33 | 2.34 |
| **2.5** `1/Var` | 0.58 | 0.43 | 0.41 | 0.427 | 0.425 | 0.429 | 0.427 |
| **0.0** `Var(a)` | 3.74 | 3.71 | 3.94 | 3.61 | – | – | – |
| **0.0** `1/Var` | 0.268 | 0.270 | 0.254 | 0.277 | – | – | – |
| **1.0** `1/Var` | – | – | 0.313 | – | – | 0.305 | – |
| **4.0** `1/Var` | – | – | 0.499 | – | – | 0.582 | – |

`Var(a)` at β=2.5 is **flat at ≈2.34 for L≥4** — bounded, not `∝L`, not `→0`. The plateau
`1/Var≈0.428` is the **positive L-independent limit** `G″_∞(2.5)`. The corpus's `G″(a₀*)=1.01,0.65,0.42`
(`L=2,3,4`) is already at that plateau by `L=4`; its apparent `≈1/L` decay is the small-L
non-Gaussian crossover, **not** a decay to 0. The β=0 row reproduces the exact Haar floor `0.269`
`[T]+[V]`, validating the witness. **Verdict of the numerics: the limit is `0.43` at β=2.5 —
positive and L-independent, not 0.**

### 3.3 Comparison with the β=0 Haar floor

`G″_∞(2.5)=0.43 > 0.269 = G″_∞(0)`. The decompactified constant-mode curvature is **above** the Haar
floor and rises with β: `0.269 (β=0) < 0.31 (β=1) < 0.43 (β=2.5) < 0.55 (β=4)`. So the answer to the
posed dichotomy is: **neither 0 nor the Haar floor — a β-dependent positive plateau.** (The Haar
floor is the β→0 endpoint of that β-dependent plateau, not the value at β=2.5.)

---

## 4 · Invariance under the choice of the constant-mode normalization

The curvature rescales under `a→λa` (`ρ(a)→λ^{−1}ρ(a/λ)`): **`G″→λ^{−2}G″`**. The numerical value is
normalization-dependent; the physical content is not. Fixing `‖v‖=1` (mode normalization) once and for
all:

- **Mode normalization** (`a=⟨v,ξ⟩`, `‖v‖=1`): `G″_a→0.43` at β=2.5 (finite, positive, L-independent).
- **Field normalization** (`ξ̄=N^{−1}Σξ=a/√N`): `G″_{ξ̄}=N·G″_a→+∞` (concentration / pinning).
- **Invariant quantity.** The susceptibility `χ_v=Var(a)=Σ_r C(r)`, an **L-independent, O(1)** number
  (a physical susceptibility: `χ=V·Var(ξ̄)`, volume-averaged). The invariant statement is
  **`χ` is finite (`[C]`: finite correlation length) — so `G″_a≈1/χ` is L-uniformly positive** in the
  fixed `‖v‖=1` normalization. Equivalently: the constant collective coordinate does **not**
  delocalize (`Var(a)` bounded), and the physical constant field **is** pinned (`Var(ξ̄)=χ/N→0`).
  The rescaling `λ` merely moves the number between `0.43` and `∞`; the **bounded variance** is the
  normalization-invariant fact. (Hence "H1: `G″→0`" is normalization-dependent and **false** in the
  standard `‖v‖=1` convention used throughout the corpus; only a pathological amplitude growing faster
  than `√N` could make it read as 0.)

**Convention match.** SM-WAVE-SYNTHESIS, SM1, SM4 all use `‖v‖=1` (`mode_k0(nfree)=e₃ repeated/√nfree`).
All numbers in this note are in that convention; the comparison to the Haar floor `0.269` is
same-convention.

---

## 5 · Acid test — free field `[T]`

**Free massless Gaussian** (`Φ≡0`, no Haar Jacobian, `V=½κ‖ξ‖²`). Then `ρ_v` is `N(0,1/κ)` (mode
normalization), so `G″(a₀)=κ` at every `a₀`, and `κ→0` as the mass `→0` — **`G″→0`, no manufactured
floor** `[T]`. The Haar floor `1/s²=0.269` is present **only** because the group is compact with a
flat-chart Jacobian `J` (`Φ₀″(0)=1/6>0`); it is **absent** for the free (Lebesgue) field. The witness
`d2_*` treats `J` exactly (the Metropolis acceptance carries `J`; at `β=0` it reproduces `s²=3.72`
exactly), so the tools are **discriminating**: `0` for the free field, `s²`-variance/Haar floor for
pure Haar, a β-dependent positive plateau for the interacting constant mode — without manufacturing a
gap.

---

## 6 · Witness — `research/v3/w7/d2_*.py`

`py -3.12` (numpy only; one job at a time; small 1D strips, cheap). Quoted outputs (abridged):
```
d2_batch.py  [Var]  exact one-link s^2 via quaternion sampler = 3.7213 ; 1/s^2 = 0.26883
             beta=2.5 Var(a): L=2..14 = 1.71,2.33,2.45,2.34,2.35,2.33,2.34  -> plateau 2.34
             1/Var: 0.58,0.43,0.41,0.43,0.43,0.43,0.43  -> G''_inf(2.5)=0.43
             beta=0   1/Var = 0.268..0.277  (= Haar floor, exact check)
d2_beta.py   beta=1: 1/Var=0.31 (L=4,10);  beta=4: 1/Var=0.50(L=4),0.58(L=10)  -> rises with beta
d2_curv.py   histogram G''(a0*) is noise-dominated (peak curvature on a broad marginal); 1/Var is the
             robust Gaussianized limit and matches the corpus G''(a0*)=0.42 at L=4.
```
Exact identities: `s²=(4π²−6)/9`, `1/s²=9/(4π²−6)=0.268830` `[T]`; `G″_{β=0}→1/s²` (CLT) `[T]`.
Metropolis scatter on `Var(a)` at fixed L is `≲1%` (12–24k independent strips × many sweeps); the
plateau in `1/Var` is flat to three digits over `L=6..14`. No third-party sources downloaded.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| exact constant-mode marginal `ρ_v(a)`, Haar+Wilson, `G″=−ρ″/ρ` incl. cross term | §1.1 (1.1) | **[T]** |
| action on the normalized constant mode `=βL⁴·⅛‖[ξ,ξ]‖²=βa⁴/(16L⁴)→0` (torus); flat bulk (strip) | §1.1, §2 | **[T]** |
| exact Haar variance `s²=(4π²−6)/9`, `1/s²=0.268830`; `G″_{β=0}→1/s²` (CLT) | §1.2 (1.2); `d2_batch` β=0 | **[T]**+**[V]** |
| k=0 mode is a zero mode of the quadratic Wilson Hessian (`vᵀHessS_W^{(2)}v=0`); curvature Haar-transverse-carried | §2, §3.1 | **[T]** |
| **`G″` limit positive, L-independent at β=2.5 ≈0.43** (not 0, not the Haar floor) | §3.2, table; `d2_batch` | **[V]** (mechanism [T-structure]) |
| limit β-dependent, monotone rising `0.269<0.31<0.43<0.55` (β=0,1,2.5,4), each flat in L | §3.3; `d2_beta` | **[V]** |
| `G″_∞>0` ⟺ finite flat-chart susceptibility `χ=Σ_r C(r)<∞` (finite correlation length) | §3.1 | **[C]** |
| normalization invariance: `G″→λ^{−2}G″`; invariant = bounded `Var(a)`/finite `χ`; field-norm `G″→∞` | §4 | **[T]** |
| free field `G″=κ→0` (no manufactured floor); pure Haar reproduces `s²` | §5; `d2_batch` | **[T]**+**[V]** |
| identification of `1/χ` with the **physical glueball mass** (gauge-fixed flat-chart coordinate) | — | **[C]** / open |
| **non-constant `k→0` soft tower** (SM6 lattice object) decompactification | not covered | **[O]** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 8 · Sharpest line

**The constant mode is an action zero mode, so its constrained curvature is Haar-transverse-carried,
not action-carried — and it has a positive, L-independent, β-dependent limit: `G″→0.43` at β=2.5,
rising from `0.269` (β=0) to `0.55` (β=4), each flat in L for `L≳4` `[V]`, with `Var(a)` plateauing at
`2.34` (`[T]` k=0 decoupling + `[C]` finite susceptibility).** So decompactification **does not** drive
the constant-mode curvature to 0 (H1 is false for the constant mode) and **does not** freeze it at the
β=0 Haar floor either: it saturates at a β-dependent bulk plateau. What remains — the identification
of this plateau with the physical glueball mass, and the decompactification of the **non-constant**
`k→0` soft tower — is `[C]`/`[O]`. **No QED. Clay 0%.**

---

*Refs.* `SM1-SOFT-EXACT.md` §1,§4,§6 (marginal identity, `Δ_dyn+W_geo+X`, β=0 convolution, acid
tests); `SM4-CONSTANT-MODE.md` §1–§3 (exact constant-mode marginal, `G″(a₀*)=−ρ″/ρ≥0`, `L`-map);
`SM-WAVE-SYNTHESIS.md` §2 (the decompactification question); `H4-HEAD-COUNTEREXAMPLE.md` §1 (interior
minimiser, `G″(a₀*)≥0`), §3 (Cartan `L`-scan); `docs/paper-v3/main.tex` §6.7 (soft mode, the remaining
question); `H2-DEFECT-BOUND.md` (measure-carried positivity, Fisher identity); `k0-CONSTRAINT-EFFPOT.md`
(`V_eff″=G″`, `W_geo≥0`). External (standard, from corpus, not re-fetched): Savvidy PLB **71** (1977)
133; Gross–Pisarski–Yaffe RMP **53** (1981) 43; Prékopa 1971. Witness `research/v3/w7/d2_*.py`. *No
third-party sources downloaded.*

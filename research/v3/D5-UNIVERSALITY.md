# D5-UNIVERSALITY — is the constant-mode curvature's L-dependence universal? **No common exponent and no universal value; at the twist boundary the whole group/action family *saturates* to a positive L-independent floor (`G″→1/χ`, the inverse amplitude susceptibility), the ≈1/L is a small-L transient, and the floor value is group/action/boundary-dependent.**

> **2026-11-24 · Wave D5** — universality audit of the soft/constant-mode decompactification claim.
> Question: is the `≈1/L` decay (or Haar-floor saturation) of the constant-mode curvature `G″(a₀*)`
> universal, or an artifact of SU(2) Wilson? Attack: recompute the L-dependence across **groups**
> (U(1), SU(2), SU(3)) and **actions** (Wilson, fundamental–adjoint), in a normalisation-independent
> way, with U(1) as the exactly solvable anchor.
> Reads: `research/v3/SM1-SOFT-EXACT.md`, `SM4-CONSTANT-MODE.md`, `SM-WAVE-SYNTHESIS.md`,
> `docs/paper-v3/main.tex §6.7`; machinery in `research/v3/w7/sm1_*.py`.
> Witness: `research/v3/w7/d5_universality.py` (`py -3.12 [quick|full]`, numpy only; exact
> Gauss–Legendre at L=1, Haar-importance MC for L≥2 with ESS quoted; U(1) exact product quadrature
> to L=4). **One job at a time.** No third-party sources.
>
> **Labels (charter §2).** **[T]** exact identity/analytic · **[V]** numerical witness · **[C]**
> conditional · **[NG]** no-go for a named class · **[O]** exact missing step. **No QED. Clay 0%.**

---

## 0 · Verdict (no QED; Clay 0%)

**The L-dependence of the constant-mode curvature is *not* universal in either exponent or value. What *is* universal (at the twist boundary the corpus uses) is the *structure*: the amplitude variance `Var(a)` saturates to a finite L-independent susceptibility `χ`, so `G″(a₀*) → 1/χ > 0` — a positive L-independent floor, not a `1/L→0` decay. The `≈1/L` seen in SM1 is the small-L transient (the ratio `R=G″·Var` runs `1.7→1.0`), not the asymptote. The floor *value* `1/χ` is group-, action-, β- and boundary-dependent, so no universal number exists; and at the trivial boundary (`b=0`) the same chains instead *grow* `Var(a)` (⇒ `G″` decays), so even the sign of the L-trend is boundary-dependent.**

Point by point:

1. **Exactly solvable anchor — U(1) L=1 [T].** On the minimal torus the U(1) marginal is the von Mises
   `ρ(a) ∝ e^{−β(1−cos(a+b))}`, and `G″(a₀*) = β` **exactly** (verified to 5 digits, all β). There is
   **no Haar floor**: U(1) Haar is flat (`J=1`, `Φ=0`), so at β=0, L=1 the density is flat and `G″=0`.
   The positive `β=0` "Haar floor" of SU(2) (`0.205 > N/12`) is thus a **non-abelian** (chart-Jacobian)
   phenomenon `[T]`. → the *floor at β=0* is **not** universal.

2. **The L-trend is a saturation, not a power law [V].** Measuring the knot-free, group-universal
   amplitude variance `Var(a) = E[a²]−E[a]² = (1/L)Σ_{ij}Cov(ξ_i,ξ_j)` (the amplitude susceptibility)
   at the twist boundary `θ=π`:

   | group·action | β | L=1 | 2 | 3 | 4 | 6 | 8 | `1/Var_∞` (floor) |
   |---|---|---|---|---|---|---|---|---|
   | U(1) Wilson | 1.0 | 5.19 | 3.49 | 2.68 | 2.66 | 2.47 | 2.40 | **≈0.42** |
   | U(1) Wilson | 2.5 | 6.81 | 3.92 | 2.50 | 2.76 | 2.48 | ~2.1 | **≈0.40** |
   | SU(2) Wilson | 1.0 | 3.31 | 2.82 | 3.11 | 3.18 | 3.22 | 3.25 | **≈0.31** |
   | SU(2) Wilson | 2.5 | 2.05 | 1.72 | 2.30 | 2.35 | 2.31 | 2.36 | **≈0.43** |
   | SU(3) Wilson | 1.0 | 2.05 | 1.96 | 1.98 | 1.99 | – | – | **≈0.50** |
   | SU(3) Wilson | 2.0 | 2.05 | 1.89 | 1.93 | 1.96 | – | – | **≈0.51** |

   (U(1) rows exact by product quadrature to L=4; validated against MC to 4 digits. L≥6 MC,
   ESS 40k–1×10⁶. SU(3) L≥2 importance-sampled, ESS ≥5×10⁵.) **All three groups saturate to a finite,
   positive, L-independent `χ` ⇒ `G″→1/χ>0`.** No `1/L→0`.

3. **Direct confirmation that `G″` plateaus (SU(2), β=2.5, twist) [V].** Re-computing the *direct*
   conditioned curvature (not the proxy) with the SM1 reweight estimator:

   | L | 2 | 3 | 4 | 5 | 6 |
   |---|---|---|---|---|---|
   | `G″(a₀*)` | 0.995 | 0.634 | 0.448 | 0.494 | **0.482** |
   | `Δ_dyn` | 0.769 | 0.381 | 0.193 | 0.235 | 0.202 |
   | `W_geo` | 0.264 | 0.294 | 0.276 | 0.275 | 0.270 |
   | `R=G″·Var` | 1.71 | 1.46 | 1.05 | 1.14 | 1.14 |

   `G″` **drops 0.99→0.45 then plateaus ≈0.48 for L≥4** (it does *not* continue as `1/L`: 1/L would
   give 0.34 at L=5, 0.28 at L=6). The decay is entirely in `Δ_dyn` (the dynamical part converging by
   L=4); **`W_geo` — the geometric/Haar part — is constant ≈0.27 for every L** `[T]+[V]`. The apparent
   local exponent `p = −dlnG″/dlnL` is `≈1.1–1.2` over L=2→4 but collapses to `≈0` over L=4→6:
   **the "≈1/L law" is a finite-L transient, not an asymptotic power law.**

4. **The `1/L` transient is the Gaussianisation transient [T]+[V].** `R=G″·Var`: L=2→1.71, L=4→1.05,
   L=5,6→1.14. As the `k=0` mode delocalises, by the central limit the marginal of `a=Σξ_i/√L`
   Gaussianises, so `G″→1/Var` (R→1, verified) and the remaining L-dependence is that of `1/Var=1/χ`.
   The corpus's own three points (1.01, 0.65, 0.42 — **reproduced here as 0.995, 0.634, 0.448**) trace
   this transient; fitting `1/L` through L=2,3,4 extrapolates the *wrong* asymptote.

5. **No universal floor value; it is action-dependent [V].** Exact L=1 `G″(a₀*)` (LONG, θ=π) moves
   monotonically with the action while staying positive:
   `Wilson(β=2.5)=0.693`, `FA(2,0.5)=0.775`, `FA(1.5,1)=0.888`, `FA(1,1.5)=1.029`, `FA(0.5,2)=1.190`,
   `Wilson(β=5)=1.325`. In the fundamental–adjoint family the variance saturation is equally present
   (β_F,β_A)=(2,0.5): Var = 1.83,2.42,2.48,2.41; (1.5,1): 2.06,2.66,2.68,2.59 — **the saturation
   *structure* is action-robust, the *value* is not.**

6. **Boundary-dependence kills a blanket universality claim [V].** At the **trivial** boundary `b=0`
   the same chains *grow*:
   U(1) exact β=1 Var = 1.60,1.43,2.03,**2.72** (L=1..4, rising); SU(2) β=1 Var = 2.40,3.08,3.63,**4.32**
   (L=2,3,4,6, rising). So at `b=0` `G″` decays over the accessible range (the constant mode is a
   near-zero mode there — the free-field/abelian flat-direction limit), while at the **twist `b=π`** it
   saturates. The corpus works at `θ=π`; the saturation statement is *that* boundary's.

7. **Normalisation-invariant statement.** Under a chart rescaling `ξ→cξ` (`a→ca`), `G″→G″/c²` and
   `Var→c²Var`, so `R:=G″·Var` is **dimensionless and chart-invariant**, and the floor is
   `G″_∞ = 1/χ` with `χ = lim_L Var_ρ(a)` a chart-covariant susceptibility. The invariant content is:
   **`R_L → 1` as `L→∞` (marginal Gaussianises), and `χ_L → χ_∞ ∈ (0,∞)` at the twist boundary**, so the
   floor `1/χ_∞ > 0` is a well-defined *inverse amplitude susceptibility*. Only the *structure*
   (positivity & L-independence of `1/χ_∞`) is universal; `χ_∞` (hence `1/χ_∞`) is not.

**One line.** The L-dependence is **not** universal: the exponent (`≈1.1` at L=2–4 → `0` at L≥4) and the
floor value (`1/χ_∞` ∈ {0.31,0.40,0.43,0.50,…}, action/boundary-dependent) both vary. What survives is
the **invariant structural statement** `G″(a₀*)·Var_ρ(a) → 1` and `Var_ρ(a) → χ_∞ > 0` at the twist
boundary — a positive L-independent floor carried by the amplitude susceptibility, not a `1/L` decay and
not a universal number. **Decompactification of the *value* is not universal; only the *saturation
structure* is.** Clay 0%. No QED.

---

## 1 · The object, the proxy, and why `Var(a)` is the right invariant

**Object (SM1/SM4).** Constant-mode curvature `G″(v;a₀*) = −d²/da² log ρ_v(a)|_{a₀*}`, `a=⟨v,ξ⟩`,
`v=(e₃,…,e₃)/√L`, `a = Σ_i ξ_{i,3}/√L`, on the `1×L` chain with fixed boundary holonomy `θ` (flat
exponential chart `U=e^{iξ·σ/2}`, Haar Jacobian `J(ρ)=sin²(ρ/2)/ρ²`). At an interior mode `a₀*`,
`G″=−ρ″/ρ ≥ 0` `[T]` (H4).

**Why `Var(a)`.** The marginal `ρ_v(a)` of a sum of `L` compact variables is a **spline with knots**
(where the support combinatorics change); at a knot `G″` is distributional. E.g. the exact U(1)
`L=2, b=0` marginal has a cusp at `a=0` and its 5-point `G″(0)` is an artifact (we measured 150 vs the
true `2β=1`). The **variance** `Var(a)=E[a²]−E[a]²` is an unrestricted moment of the measure (no
marginal, no constraint, no knot), smooth in `L`, and equals the *amplitude susceptibility*
`χ_L=(1/L)Σ_{ij}Cov(ξ_i,ξ_j)` — the object that controls whether `G″≈1/χ_L` decays (if `χ_L` grows) or
saturates (if `χ_L→const`). It is computed group-universally by Haar-importance sampling with **no chart
Jacobian needed** (sampling Haar *is* the Jacobian), which is why U(1)/SU(2)/SU(3) are on the same
footing. For a Gaussian marginal `G″=1/Var` exactly; the ratio `R=G″·Var` measures the departure. `[T]`

**Validation `[V]`.** U(1) exact product-quadrature `Var(a)` (L=1,2) matches the Haar-importance MC to
4 digits: β=1: exact 5.191/3.492 vs MC 5.187/3.494; β=2.5: 6.814/3.915 vs 6.814/3.911. The engine is
sound; no third-party input.

---

## 2 · Group universality — the table and the exponents

**U(1) — exactly solvable.** L=1: `G″=β` exactly `[T]`. The `1×L` chain marginal is a spline; its exact
`Var(a)` (product quadrature, `NX=140` in `full`) at `b=π`: β=1 → 5.19, 3.49, 2.68, 2.66 (L=1..4);
β=2.5 → 6.81, 3.92, 2.50, 2.76. **Saturates** (the decrease is a transient; L=4 already ≈ the L=8 MC
value 2.4/2.5). `G″≈1/Var→0.40`.

**SU(2) Wilson — the corpus model, extended.** Direct `G″` (§0.3) plateaus 0.48 for L≥4; `W_geo≈0.27`
constant. `Var(a)` at `b=π`: β=1 → 3.31,2.82,3.11,3.18,3.22,3.25 (L=1..8, saturating to 3.25); β=2.5 →
2.05,1.72,2.30,2.35,2.31,2.36 (saturating to 2.36). L=1 `G″` exact `0.694` (Cartan) / `0.362`
(transverse); β=0 Haar floor `0.205>N/12`.

**SU(3) Wilson.** `Var(a)` (diagonal Cartan direction, `λ_3`) at `b=π`: β=1 → 2.05,1.96,1.98,1.99
(L=1..4); β=2 → 2.05,1.89,1.93,1.96. **Saturates to ≈1.97** → floor `1/χ≈0.51`. (SU(3) is computed via
direct group-Haar sampling of `Var(a)`; the flat-chart `G″` needs the 8-dim chart Jacobian and is not
computed here — but `R→1` ties it to `1/Var`.)

**The exponent, honestly.** Fit `Var∝L^α` (equivalently `G″∝L^{−α}`):
- SU(2) β=2.5, over L=2→4: `α≈1.2` (Var 1.72→2.35); over L=4→8: `α≈0.0` (2.35→2.36).
- U(1) β=1, over L=1→4: `α≈−0.4` (Var *decreases*, 5.19→2.66); over L=4→8: `α≈0.0`.
So there is **no single power law**: the local exponent is group- and range-dependent and collapses to
`0` (saturation) once `L≳4`. The corpus's `≈1/L` is the `α≈1.2` piece of SU(2) at L=2–4 only.

---

## 3 · Action universality — same structure, different value

Exact L=1 `G″(a₀*)` (LONG `e₃`, θ=π), fixed plaquette weight `w(P)`:

| action | `(β_F,β_A)` | `a₀*` | `G″_LONG` |
|---|---|---|---|
| Wilson | (2.5, 0) | −2.14 | **0.693** |
| Fund–Adj soft | (2.0, 0.5) | −2.28 | 0.775 |
| Fund–Adj mid | (1.5, 1.0) | −2.38 | 0.888 |
| Fund–Adj adj | (1.0, 1.5) | −2.48 | 1.029 |
| Fund–Adj hard | (0.5, 2.0) | −2.56 | 1.190 |
| Wilson | (5.0, 0) | −2.57 | 1.325 |

Every action gives `G″>0` at the mode `[T]`; the **value** rises monotonically with the adjoint
admixture / β. The **saturation** persists in the FA family (Var at `b=π`, `(2,0.5)`: 1.83,2.42,2.48,2.41;
`(1.5,1)`: 2.06,2.66,2.68,2.59). So *action-independence of the structure, action-dependence of the value*.
(Scope: an explicit heat-kernel weight `w_t(P)=−log Σ_j(2j+1)e^{−t j(j+1)}χ_j(P)` is not separately run
— it is a truncated character sum and sits on the **same multi-character axis** the FA family already
spans; the FA scan is the tractable stand-in for the heat-kernel/alternative-action test.)

---

## 4 · Boundary dependence — the caveat that forbids a blanket statement

The twist `θ=π` (Z₂ centre) is the corpus's choice and the one where the floor appears. At the **trivial**
boundary `θ=0` the same chains *grow*:

| model | b | L=1 | 2 | 3 | 4 | 6 | trend |
|---|---|---|---|---|---|---|---|
| U(1), β=1 (exact) | 0 | 1.60 | 1.43 | 2.03 | 2.72 | – | **rising** |
| U(1), β=1 (exact) | π/2 | 2.41 | 2.10 | 2.58 | 2.56 | – | saturating |
| U(1), β=1 (exact) | π | 5.19 | 3.49 | 2.68 | 2.66 | – | saturating |
| SU(2), β=1 (MC) | 0 | – | 2.40 | 3.08 | 3.63 | 4.32 | **rising** |
| SU(2), β=1 (MC) | π | 3.31 | 2.82 | 3.11 | 3.18 | 3.22 | saturating |

At `b=0` the constant mode is a near-zero mode (the abelian/free flat direction survives), so `Var` grows
and `G″` decays over the accessible range — the **free-field acid-test regime** (§6). The twist lifts the
zero mode and produces the positive floor. **A claim of the form "the L-dependence is universal" is false
across boundary conditions; the correct statement must fix `θ=π`.** `[V]`

---

## 5 · Normalisation invariance — the invariant statement

The only chart-invariant combination of the two natural quantities is
```
   R(a₀*) := G″(a₀*) · Var_ρ(a)          (dimensionless; invariant under ξ→cξ),
```
and the floor is `G″_∞ = 1/χ_∞`, `χ_∞ = lim_L Var_ρ(a)`. Measured:
- SU(2) β=2.5: `R` = 1.71 (L2), 1.46 (L3), 1.05 (L4), 1.14 (L5,6) → **`R→1`** (Gaussianised marginal).
- U(1) β=2.5 L=1: `R = β·Var = 2.5·6.81 = 17.0` (strongly non-Gaussian von Mises) → decreases to 1 as L grows.

**Invariant content `[T]+[V]`.** *(i)* `R_L → 1` as `L→∞` at the twist boundary (CLT on `a=Σξ_i/√L`).
*(ii)* `χ_L → χ_∞ ∈ (0,∞)` (Var saturates). *(iii)* Hence `G″_∞ = 1/χ_∞ > 0` — a positive L-independent
floor, the inverse amplitude susceptibility. `χ_∞` is **not** universal (§2–§4); only *(i)+(ii)* are.

---

## 6 · Acid test (free field) — passes

- **Free massless Gaussian** (`Φ≡0`, `S=κ|ξ|²/2`): `ρ(a)` Gaussian `⇒ G″=κ→0`, `Var=1/κ`, `R=1`. No
  manufactured floor `[T]`.
- **U(1) Haar (β=0)** at L=1: flat density, `G″=0` `[T]`. The SU(2) β=0 "Haar floor" `0.205>N/12` is the
  **chart-Jacobian** curvature, absent for U(1) — the floor is non-abelian at β=0.
- **Trivial boundary `b=0`**: the near-zero mode ⇒ `Var` grows ⇒ `G″` decays (§4) — the free-field
  regime; the tool returns 0/decay there, non-vacuously.
So the tool gives 0 on the free/Gaussian field, the `N/12` Jacobian value on flat SU(2) Haar, and a
positive saturating `1/χ` at the interacting twist head — discriminatingly `[T]+[V]`.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| U(1) L=1 `G″(a₀*)=β`; β=0 flat `G″=0` | von Mises identity | **[T]** |
| interior-mode sign `G″(a₀*)=−ρ″/ρ≥0`, all groups/actions | H4 specialised | **[T]** |
| `G″≈1/Var` in Gaussian limit; `R=G″·Var` chart-invariant; `R→1` | §1,§5 | **[T]** |
| **`Var(a)` saturates (twist `θ=π`): U(1)≈2.4, SU(2)≈2.3–3.25, SU(3)≈1.97** | §2, tables | **[V]** |
| **direct SU(2) `G″` plateaus ≈0.48 for L≥4; `W_geo≈0.27` L-independent** | §0.3 | **[V]** |
| `≈1/L` is a small-L transient (`α≈1.2`→`0`), **not** asymptotic | §2 exponents | **[V]** |
| floor value `1/χ_∞` action-dependent (0.69…1.33 at L=1; Var-sat. FA) | §3 | **[V]** |
| floor **structure** action-robust; **value** not | §3 | **[V]** |
| at `b=0` `Var` *grows* (U(1), SU(2)) ⇒ `G″` decays | §4 | **[V]** |
| floor at β=0 is **non-abelian** (SU(2) `N/12`, U(1) `0`) | §6 | **[T]**+**[V]** |
| **uniform-in-L floor `χ_∞>0` for all β/boundary/group** (the invariant statement, proved) | — | **[O]** = H1 |
| transport of the 1D-chain floor to **4D `L→∞` decompactification** | SM6 hints softening | **[C]/[O]** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 8 · Sharpest line

**The constant-mode curvature's L-dependence is not universal — neither its exponent nor its value.** The
exactly solvable U(1) anchor gives `G″=β` at L=1 with **no** β=0 Haar floor (that floor is non-abelian);
across U(1), SU(2), SU(3) and Wilson vs fundamental–adjoint, the amplitude variance `Var(a)` *saturates*
at the twist boundary `θ=π` (to `χ_∞` ≈ 2.4, 2.3–3.3, 2.0 respectively), so `G″→1/χ_∞>0` — a positive
**L-independent** floor, and the corpus's `≈1/L` (`1.01,0.65,0.42`, reproduced as `0.995,0.634,0.448`) is
the small-L Gaussianisation transient (`R=G″·Var`: 1.7→1.0), *not* the asymptote (`G″` plateaus at 0.48
for L≥4, while `W_geo≈0.27` is exactly L-independent). The floor **value** `1/χ_∞` is group/action/β-
dependent, and at trivial boundary `b=0` the same chains instead grow `Var` (⇒ `G″` decays, the free-field
regime). The only **invariant** statement is `R=G″·Var→1` and `Var→χ_∞∈(0,∞)`: a positive L-independent
floor that is the *inverse amplitude susceptibility*, whose saturation is universal in **structure** but
whose value — and the 4D decompactification of it — remain **[O] = H1**. **Clay 0%. No QED.**

---

*Refs.* `SM1-SOFT-EXACT.md` §0,§4 (direct `G″` 1.01/0.65/0.42, `W_geo≈0.25`, L-scan); `SM4-CONSTANT-MODE.md`
§2–§3 (exact L=1 marginal, `G″(0)=A−Bβ²`); `SM-WAVE-SYNTHESIS.md` §2 (decompactification question);
`docs/paper-v3/main.tex` §6.7. Witness `research/v3/w7/d5_universality.py`; U(1) exact product
quadrature, Haar-importance MC with ESS, `sm1_blocks` reweight for the direct SU(2) `G″`. Standard inputs
(Brascamp–Lieb, CLT) cited from the corpus. **No third-party sources downloaded; no QED.**

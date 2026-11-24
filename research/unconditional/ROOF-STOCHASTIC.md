# ROOF-STOCHASTIC — stochastic quantization / regularity structures: the spectral gap of the 4D gauge Fokker–Planck operator

> **VERDICT (no QED).** Stochastic quantization (Parisi–Wu) turns the gap into the **spectral gap of the
> Fokker–Planck (FP) operator** of the Langevin SPDE — but that gap **is** the mass gap `[T]`, so a gap bound
> is *circular unless its origin is non-circular*. Modern singular-SPDE / regularity-structure theory
> (Hairer) supplies a **rigorous, non-circular** construction — and there it delivers the strut — **only up
> to the dimension where 4D Yang–Mills becomes critical.** In `d=4` the gauge coupling is **marginal**
> (`[g²]=4−d=0`, asymptotically free): regularity structures need subcriticality + finitely many
> counterterms, and the 4D counterterms **proliferate to all loop orders** = the **renormalon** = the same
> `H1` wall. **The SPDE framework does not give the strut in 4D; it relocates `H1` to "the 4D gauge SPDE is
> critical."** Clay 0%.

---

## §0 · Labels, acid, iron laws (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-here]** new.
**Acid**: the free massless field (gap `→0`) must be a solution; the gapped-atomless `:φ²:` must not be certified.
**IL-i** non-positivity · **IL-ii** nonabelian · **IL-iii** no small parameter at `α_s~O(1)`.

---

## §1 · What stochastic quantization actually says (the equivalence is a tautology)

Parisi–Wu: the Euclidean YM measure `e^{−S}dA` is the invariant measure of the Langevin equation
```
∂_t A_μ = −δS/δA_μ + ξ_μ = D_νF_{νμ} + ξ_μ (+ gauge fixing),     ⟨ξ_μ(x,t)ξ_ν(x′,t′)⟩ = 2δ_{μν}δ(x−x′)δ(t−t′).
```
The associated Fokker–Planck operator, conjugated to symmetric form, is the **Dirichlet-form generator**
```
L = −Δ_F + (∇S)·∇   on L²(e^{−S}dA),   equivalently  −Q†Q  in the BRST/SUSY grading.
```
Its spectral gap `= E₁ − E₀` **equals the mass gap** `[T]` (the FP semigroup is the transfer semigroup of the
Euclidean theory; `e^{−Lt}` = transfer). Consequences:
- **A bound on `gap(L)` is a bound on the mass gap** — and conversely. So "SPDE gives a rigorous gap lower
  bound" is, *by itself*, a restatement of the target.
- **Non-circularity must come from the SPDE analysis, not the identity.** The two standard non-circular
  handles are **hypercontractivity / log-Sobolev** (Gross) and **Bakry–Émery `Γ₂≥0`**. But `L` is exactly
  a **gradient-flow generator**, so `Γ₂`/BE on it is **Theorem BE** (`BASE-FIELDMEASURE` `[T]`:
  homogeneous in the Hessian / zero at equilibrium / sign-indefinite) — the BE handle is dead. The
  log-Sobolev/hypercontractivity handle for 4D YM is **open and equivalent to the gap** `[O]`.

---

## §2 · Where regularity structures works — and stops (the dimension count)

The rigorous tool (Hairer's regularity structures / paracontrolled distributions) builds the SPDE and its
**non-circular** spectral gap in the **subcritical** regime, where the nonlinearity is *super-renormalizable*
and finitely many counterterms suffice:

| `d` | `[g²] = 4 − d` | class | RS status |
|---|---|---|---|
| 2 | `2` | super-renorm. | `[T]` construction + gap |
| 3 | `1` | super-renorm. | `[T]` YM-Higgs: **Chandra–Chevyrev–Hairer–Shen 2022** (3D stochastic quantisation) |
| **4** | **`0`** | **marginal / critical, asymptotically free** | **`[O]`/`[NG]`: no global RS construction; counterterms to all orders** |

In `d=3` the strut is genuinely supplied (a rigorous gap for the 3D SPDE), matching the known 3D
strong-coupling gap. **In `d=4` the SPDE is critical**: `g²` dimensionless, `β(g)=−β₀g³` (asymptotic
freedom), so the required counterterms are the **full perturbative series** — one per loop order — with no
finite truncation. Regularity structures assume subcriticality precisely to avoid this; 4D YM sits *on*
the critical line.

---

## §3 · The 4D counterterms are the renormalon wall (precise)

Gauge fixing (Faddeev–Popov) + ghosts, then the renormalization constants
`Z_A, Z_c, Z_g (with Z_g Z_A^{1/2} Z_c = 1, STI), Z_α,` and the vacuum energy. In `d=4` these are needed to
**all orders**; their Borel transform has an **IR renormalon pole at `u = d/2 = 2`** (`'t Hooft 1977`),
ambiguity `~Λ⁴` — **the same `u=2` renormalon that `ROOF-CONDENSATE-FLOOR` `[NG]` proved blocks every
scheme-independent condensate floor.** So "make the 4D SPDE well-posed" and "floor the `0⁺⁺` anomaly" are the
**same obstruction in two coordinates** (`SLOT-SYNTHESIS` §5). The asymptotic series is non-summable
without a renormalon prescription; controlling it uniformly in `a→0` is the transport wall `= H1`.

---

## §4 · Verdict and the exact step

- **`[T]`** `gap(FP) = mass gap`; `d=2,3` super-renormalizable, `d=4` marginal; the FP generator is a
  gradient-flow generator (so BE is the dead handle).
- **`[T]`** regularity structures give a rigorous, non-circular gap through `d=3` (CCHS 2022 for YM-Higgs).
- **`[NG]`** SPDE/RS does **not** give the strut in 4D: the theory is **critical**, not subcritical — the
  hypothesis RS requires is exactly the one 4D YM violates.
- **`[O]`** exact missing step: a global (all-orders, renormalon-controlled) well-posedness of the 4D gauge
  SPDE — equivalently a non-gradient-flow (non-BE, non-order-parameter) hypercontractivity/log-Sobolev
  estimate for the 4D YM Fokker–Planck operator. This is `H1` in SPDE coordinates; **no such estimate exists.**

**Net:** stochastic quantization is the *cleanest existing-theory restatement* of the target — and the
dimension count shows precisely why it stops one dimension short: `d=3` super-renormalizable (strut
supplied), `d=4` critical (strut = the renormalon wall). It does **not** close Clay; it pins the wall at
"4D is the critical dimension."

Witness `w7/roof_stochastic.py` (numpy, `<1s`, `py -3.12 -I`): (A) free-field FP gap `= m → 0` as `m→0`
(acid); (B) coupling-dimension table `[g²]=4−d` (marginal at 4); (C) perturbative coefficient growth
`~β₀ⁿ n!` (renormalon, asymptotic).

*Refs:* Parisi–Wu 1981; Hairer, *Invent. Math.* 198 (2014); Chandra–Chevyrev–Hairer–Shen, "Stochastic
quantisation of Yang–Mills–Higgs in 3D" (2022); Gross 1975 (log-Sobolev ⟺ hypercontractivity); 't Hooft
1977 (renormalon); `BASE-FIELDMEASURE` (Thm BE), `ROOF-CONDENSATE-FLOOR` (u=2), `SLOT-SYNTHESIS`.

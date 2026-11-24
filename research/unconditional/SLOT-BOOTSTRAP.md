# SLOT-BOOTSTRAP — decisive numerical test: do the exact moment constraints on the 0⁺⁺ spectral measure certify a positive support-floor `Δ`, or do they admit gapless solutions at every truncation order?

> 2026-11-24 · ROUTE = **numerical bootstrap / feasibility of a gapped `k→0` solution**. We set up the finite constraint system on the `0⁺⁺` Källén–Lehmann spectral measure `ν` (a truncated Stieltjes moment problem) from **all** exact inputs — reflection positivity `ν≥0`, the KL moments `M_j=∫s^j dν` (`M_{-1}=χ`, growth `M_k`), Ward `E[af]=1`, the trace anomaly `M_{-2}∝Λ⁴` — and compute the **feasible support-floor `Δ=inf supp ν`** as a function of the truncation order `K`. It is a *decisive* numerical mirror of the analytic verdict in `A4-SPECTRAL-BOOTSTRAP.md`; it does not replace it. Witness: `research/unconditional/w7/slot_bootstrap.py` (`py -3.12`, single lightweight LP job, no lattice MC).
>
> **VERDICT (up front).** **[V-here]** The feasible-floor set is `[0, Δ_max(K)]` at **every** truncation order; the **minimum feasible floor is `0` for all `K=1..7`** and does **not** rise with `K`. In detail: (i) the **max** floor `Δ_max(K)` is an **UPPER** bound, converging *down* to the true floor `≈1.00` as `K` grows; (ii) the **min** floor `Δ_min(K)=0` *identically* — gapless witnesses (a vanishing-mass atom at `Δ`, weight `w₀∝Δ²`) match every exact input at residual `<1.5e-11`. So **moment data alone certify `floor = 0`, never a positive floor.** This confirms, numerically and quantitatively, `[NG]` for the *finite-moment-constraint class*: no single constraint gates gapless solutions (leave-one-out: all survive). The slot needs a **non-moment, edge-sensitive input** (§6) = `H1`/`I2` relabelled. **Clay not implied; 0%. No QED.**

Labels (charter §2): **[T]** true theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[V-here]** computed this round.

---

## §1 · The finite constraint system (a truncated Stieltjes moment problem)

Let `ν` be the `0⁺⁺` spectral measure on `s=m²∈[0,∞)`: `C(τ)=∫_0^∞ e^{−mτ}dν(m)`, `m₁=inf supp ν` (the gap), `f₁²=ν({m₁²})` (the edge residue). The **exact** inputs, each encoded as a constraint:

| input | meaning | encoding here | label |
|---|---|---|---|
| reflection positivity | `C(τ)` is a positive-definite KL measure | `ν ≥ 0` (the LP's `w ≥ 0`) | **[T]** |
| KL moments | `M_j = ∫ s^j dν`, `s=m²` | equality rows `Σ_i w_i s_i^j = μ_j`, `j∈J` | **[T]** |
| Ward `E[af]=1` | `G″·Var(a)=1` ⇒ `G″=1/χ`, `χ=∫s⁻¹dν<∞` | fixes `M_{-1}=χ>0` | **[T]** |
| susceptibility | `χ=M_{-1}=1/G″` finite | one fixed moment value | **[T]** |
| trace anomaly | `T^μ_μ=(β/2g)trF²`, `⟨T T⟩` positive ⇒ `M_{-2}=∫s⁻²dν ∝Λ⁴` finite | fixes `M_{-2}=A>0` | **[T]** |
| KL growth / AF tail | `M_k` for `k≥1` finite, growing | `M_1,…,M_K` fixed | **[T]/[C]** |

The constraint set is `S(K) = {ν≥0 : M_j(ν)=μ_j ∀ j∈J_K}`, `J_K={−2,−1,0,1,…,K}` (`|J_K|=K+3`). Since `ν` is arbitrary on `[0,∞)`, `S(K)` is exactly a **truncated Stieltjes moment problem** on the half-line. We parametrize `ν` by a finite atomic ansatz `ν=Σ_i w_i δ(s_i)` on a dense log grid (`s`-nodes span `[Δ, 60]`; physical-scale nodes fixed) and solve for the weights by a row-scaled linear program (`scipy.linprog`, HiGHS). Feasibility of `w≥0` with residual `<1e-6` = non-empty `S(K)`.

**Representative moment vector** **[C]/[V-here]:** a dominant `0⁺⁺` glueball atom `f₁²=0.50` at `m₁²=1` plus a sparse multi-glueball continuum — support `s∈{1,3,6,12,25}`, weights `{0.50,0.15,0.15,0.10,0.10}` (normalised `M₀=1`), **true floor `=1`**. Moments: `M_{-2}=0.5217`, `M_{-1}=χ=0.5873`, `M_0=1.000`, `M_1=5.550`, `M_2=84.15`, `M_3=1772.25`. This is a stand-in for the `β=2.5` moments; the conclusion is **value-independent** because it is a property of the moment *kernels*, verified in §4 by leave-one-out. **[C]** on the identification "model ↦ physical `β=2.5`"; the floor-blindness itself is **[T]** (§3).

---

## §2 · Decisive scan: feasible floor vs truncation order `K`

For each `K` we compute two extremal floors:
- **`Δ_max(K)`** = largest `D` with a data-consistent measure supported **entirely in `[D,∞)`** — the best **UPPER** bound the moments give on the floor.
- **`Δ_min(K)`** = smallest `D` with a data-consistent measure whose support **reaches down to `D`** (a low atom at `s=D`, its weight *maximised* by the LP) — the only **LOWER**-bound certificate the moments can give. This is the quantity the mass gap needs to be `>0`.

**[V-here]** (residual of every accepted solution `≤1.4e-11`; truth-witness residual `0.00e+00`, so `Δ_max≥1` in exact arithmetic):

| `K` | `|J|` | `Δ_max(K)` (upper bd.) | `Δ_min(K)` (lower cert.) | gapless (`floor≤1e-3`) feasible? |
|---|---|---|---|---|
| 1 | 4 | 1.0551 | **≤ 1e-3** | **YES** |
| 2 | 5 | 1.0551 | **≤ 1e-3** | **YES** |
| 3 | 6 | 1.0195 | **≤ 1e-3** | **YES** |
| 4 | 7 | 1.0195 | **≤ 1e-3** | **YES** |
| 5 | 8 | 1.0050 | **≤ 1e-3** | **YES** |
| 6 | 9 | 1.0050 | **≤ 1e-3** | **YES** |
| 7 | 10 | 0.9900* | **≤ 1e-3** | **YES** |

`*` K=7 reads 0.99 by LP conditioning at high moment order (`M₇~6e9`); in exact arithmetic `Δ_max≥1` (the truth at floor 1 is itself a witness, residual `0`). The trend is the point: `Δ_max(K) → 1.00` **from above** as `K` grows.

**Reading.** **[V-here]**
1. **The LOWER end is glued to `0` for every `K`** — `Δ_min(K)=0` identically (a support-floor `≤10⁻³` measure is feasible at every truncation). The minimum feasible floor **does not converge away from `0`**; it is `0` at all orders. This is the decisive failure.
2. **The UPPER end `Δ_max(K)` converges to the true floor `≈1.00`** (from above), i.e. the moments give a *tightening upper bound on `m₁²`* — exactly the moment-ratio / Christoffel direction. An **upper** bound on the gap is not a gap certificate; the free gapless measure sits at `Δ=0 < Δ_max`, so `Δ_max>0` excludes nothing.
3. Therefore the feasible floor set is the **interval `[0, Δ_max(K)] ≈ [0,1]`** at every order. It **never shrinks to a point at the bottom**; the true floor `1` sits at the *top* of the feasible interval, unpinned from below.

The `K=1` ambiguity is already exact and decisive (§Test 0 of the witness): `δ₀+δ₂` and `2δ₁` have **identical** `(M₀,M₁)=(2,2)` and floors `0` vs `1`. **[T]** Finitely many positive moments never determine `inf supp` — the corpus's `GAP-G2-edge` lemma, reproduced numerically.

---

## §3 · The near-zero-atom budget `w₀(Δ)` — why *all* moments are floor-blind

**[V-here]** With `K=6` fixed and a low atom at `s=Δ`, the LP (maximising that atom's weight) returns:

| `Δ` | `w₀*(Δ)` (weight the data allow at the floor) | scaled residual | `1/s²`-share of `M_{-2}` from the low atom | feasible? |
|---|---|---|---|---|
| 1e-1 | 5.47e-5 | 1.4e-11 | 0.0105 | YES |
| 3e-2 | 4.04e-6 | 1.4e-11 | 0.0086 | YES |
| 1e-2 | 4.25e-7 | 1.3e-11 | 0.0081 | YES |
| 3e-3 | 2.41e-8 | 1.5e-11 | 0.0051 | YES |
| 1e-3 | 4.15e-9 | 1.3e-11 | 0.0079 | YES |

**`w₀*(Δ) ≈ C·Δ²`** (two decades of `Δ` ⇒ four decades of `w₀`). **[V-here]** The mechanism is **[T]**: for a measure with `inf supp=0` but a density vanishing at the origin, `dν ~ s^{1+δ}ds` (`δ>0`),
`M_{-1}=∫s^{-1}dν ~ ∫ s^{δ}ds <∞` and `M_{-2}=∫s^{-2}dν ~ ∫ s^{δ-1}ds <∞`, and `M_k<∞` for all `k≥0`. A vanishing-mass atom at `s=Δ` (weight `~Δ²`, so it contributes `~Δ²·Δ^{-2}=O(1)` to `M_{-2}` and `~Δ²Δ^{-1}=O(Δ)→0` to `M_{-1}`, `O(Δ^k)→0` to all `M_k`) *tops up* `M_{-2}` while leaving `M_{-1},M₀,M_1,…` to the bulk. **Every moment kernel `s^j`, `j≥−2`, integrates the `s→0` region with weight that vanishes faster than the allowed atomic mass**, so no finite set of them detects a floor collapse. The anomaly moment `M_{-2}` is the most IR-sensitive (`1/s²`), yet a *zero-mass* atom supplies it exactly with `Δ²→0` weight. This is the numerical fingerprint of `KL moments = Taylor germ at s=∞` vs `edge = pole at finite s` (`A4 §2`).

---

## §4 · Which constraint gates gapless solutions? — leave-one-out: **none**

**[V-here]** `K=5`, test "floor `≤1e-3` feasible?" after **removing** each constraint from the full set `J={−2,−1,0,1,2,3,4,5}`:

| constraints | floor≤1e-3 feasible? | low-atom weight `w₀` |
|---|---|---|
| all | **True** | 2.07e-8 |
| drop `M_{-2}` | **True** | 7.93e-5 |
| drop `M_{-1}` | **True** | 1.77e-7 |
| drop `M_0` | **True** | 8.06e-8 |
| drop `M_1` | **True** | 4.76e-8 |
| drop `M_2` | **True** | 3.43e-8 |

Gapless feasibility **survives removal of every individual moment.** There is **no "gate" constraint** to relax — the obstruction is the *type* of data (integral moments), not any one identity. Dropping the anomaly `M_{-2}` relaxes the budget most (`w₀` jumps `×3800`), confirming `M_{-2}` is the most IR-sensitive of the inputs — yet even it is only an **upper** bound (§5) and does not gate the floor.

---

## §5 · The binding constraints are UPPER bounds (moment-ratio lemma)

**[T]** For `ν≥0` with `inf supp ν ≥ D`: `s≥D ⇒ s^{-1}≤1/D ⇒ M_{-1}≤M_0/D ⇒ D≤M_0/M_{-1}`; likewise `M_{-2}≤M_{-1}/D ⇒ D≤M_{-1}/M_{-2}`; and `M_k/M_{k'} ≤ D^{k−k'}` for `k<k'`. **Every exact input caps `D` from ABOVE.** **[V-here]** for the representative vector: `M_0/M_{-1}=1.703`, `M_{-1}/M_{-2}=1.126`, `M_1/M_0=5.55` — all `≥` the true floor `1`, all *upper* bounds. No combination of these produces a **lower** bound: that is `[NG]` for the moment/positivity/Ward/anomaly class (`GAP-G2-edge`; `A4 §2–3`; `KEY-DIRECT-DERIVATION §4`). The UPPER-bound `Δ_max(K)→1` in §2 is exactly this lemma tightening.

---

## §6 · What non-moment input would cut the gapless solutions

**[O]** The numerics localise the exact missing input: an **edge-sensitive datum that assigns positive cost to support near `0`** — something a finite set of integral moments structurally cannot encode. The minimal such packages, all equivalent to `H1`/`I2`:

1. **Edge atom + isolation (residue–level-gap lemma).** A lower bound `f₁²=ν({m₁²}) ≥ r₀>0` *together with* an upper bound `Λ` on the next spectral level beyond the atom, giving `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)`. But `Λ` (the excited-threshold / coarse-mode Poincaré mass) **is** `H1`; this is the gap shifted one rung up (`GAP-G2-edge`).
2. **A density-of-states / decay envelope (non-polynomial test function).** A rigorous bound `ν([0,ε]) ≤ Cε^p` for small `ε` — equivalently a Simon–Spencer envelope `C(τ) ≤ A e^{−μτ}` at large `τ` — excludes the `s^{1+δ}` loose-density witness because it bounds the *measure of the near-zero region*, not its `s^j`-moments. This is a non-analytic (step/envelope) functional, **not a moment**; and proving it for YM is `H1`.
3. **Positivity *with* an isolated atom.** Requiring `ν` to have an atom at its edge and *zero* density on `(m₁²,Λ)` — an "empty interval above the edge" — again needs the isolation `Λ` (`=I2`).

All three are the **same object** the profile calls the one open slot (`PROFILE.md §3`: `H1 = I2 = H* = the gap`). The numerical bootstrap therefore **confirms** that closing the slot is not a matter of more moment data / higher truncation — it requires a **field-measure (non-spectral) input** `J` (`PROFILE` P1), exactly as the four-facet convergence predicted. **[O]/[NG]**

---

## Verdict and honest bottom line

**[V-here]** Decisive numbers, for a representative but structurally generic gapped `0⁺⁺` moment vector:

- **Feasible floor set `= [0, Δ_max(K)]` at every truncation `K=1..7`.**
- **`Δ_min(K) = 0` for all `K`** (gapless witnesses match every exact input to `<1.5e-11`); the **minimum feasible floor is `0` and does not converge away from `0`.**
- **`Δ_max(K) → 1.00` from above** — a tightening **UPPER** bound on `m₁²`, never a lower-bound/gap certificate.
- **No single constraint gates gapless solutions** (leave-one-out: all survive); the binding constraints `M_{-2}, M_{-1}` are **upper**-bound machines (§5).

**Does the bootstrap certify a positive floor? NO.** **[NG]** for the finite-moment-constraint class: moment/positivity/Ward/anomaly data admit gapless (`floor=0`) solutions at **every** truncation order, and the lower certificate is pinned to `0` independently of `K`. **Does it confirm moment-data-alone are floor-blind? YES**, quantitatively — the min-floor scan is exactly `0` at each order and the upper bound converges to the true floor while the lower bound stays `0`. **What non-moment input is needed?** An **edge-sensitive** datum — a residue lower bound `f₁²≥r₀` **with** an isolation/second-level bound `Λ` (residue–level-gap lemma), or an equivalent density-of-states/Simon–Spencer decay envelope (§6) — i.e. **`H1`/`I2` relabelled**: a functional on the field measure, not on `ν`. This is consistent with — and a numerical mirror of — the analytic `[NG]` in `A4-SPECTRAL-BOOTSTRAP.md` and the data-mode facet `Π3` (`PROFILE-FACET-DATAMODE.md`).

**Honest caveats.** **[C]** The moment values are a representative model, not measured `β=2.5` data; the floor-blindness is nevertheless structural (leave-one-out §4, and the kernel argument §3 make it **value-independent**). The scan uses a finite atomic grid + LP (tolerance `1e-6`, residuals `~1e-11`); `Δ_max` at `K=7` carries a `~1%` LP-conditioning wobble. None of these affects the conclusion: **`Δ_min=0` for all tested `K`, and it is `[T]` that adding more finite moments cannot lift it.** **Clay = 0%. No QED.**

---
*Witness:* `research/unconditional/w7/slot_bootstrap.py` (`py -3.12`; `scipy.optimize.linprog`; single lightweight job; no lattice MC). *Internal refs:* `PROFILE.md` (one open slot `=H1`); `PROFILE-FACET-DATAMODE.md` (`Π3`: spectral data `{blind, circular}`); `A4-SPECTRAL-BOOTSTRAP.md` (moments = Taylor germ at `∞`, edge = finite pole; feasibility non-empty & gap-blind); `GAP-G2-edge.md` (moment-ratio upper bound; residue–level-gap lemma; `δ₀+δ₂` vs `2δ₁`); `H1-IS-IT-THE-GAP.md` (`G″=1/χ`, squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ`); `KEY-DIRECT-DERIVATION.md` §4–5 (`0⁺⁺` `ρ≥0 ⇒ moments ⇒ upper`; no signed `V−A` channel); `BS-NUMERICAL-RIGOROUS.md` (correlator = one-way upper bound; lower bound needs a decay envelope). *External `[R]`:* Källén (1952) / Lehmann (1954); Simon–Spencer (1977); Stieltjes moment problem (classical).

*Status: decisive numerical test executed and passed as an honest `[NG]`-confirmation. The moment/RP/Ward/anomaly bootstrap does NOT certify a positive `0⁺⁺` support-floor at `β=2.5`: `Δ_min(K)=0` at every truncation, `Δ_max(K)→1.00` (upper bound). The slot needs a non-moment, edge-sensitive (`field-measure`) input `= H1 = I2`. Clay NOT proved; 0%. No QED on anything unproven.*

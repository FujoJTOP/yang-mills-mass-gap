# ROOF-PEIERLS-0FORM — the 0-form Peierls escapes the *circularity* but not the *small parameter*: not the same wall, still a wall

> **VERDICT (no QED).** The `0⁺⁺` (`0`-form) contour is a **genuinely different object** from the
> `1`-form (Wilson-loop) surface — a **3-surface domain wall** (codim‑1 in 4D) vs a **2-surface magnetic
> string**, carrying a **domain-wall tension `σ_w`** vs the **string tension `σ_s`**, and targeting the
> **exponential decay of the `0⁺⁺` correlator (the gap itself)** vs the **area law**. So the 0-form
> Peierls **does escape the specific circularity** flagged in the base round (`σ_s>0 ⇔ confinement`):
> it never invokes `σ_s`. **But it is equally blocked at the crossover.** At `β=2.5` the defect activity
> is `a₁/a₀ = 0.565`, so `σ := −log(activity) = 0.570` — **2–3× below the `O(1)` entropic threshold
> `c ∈ [log 3, log 7] ≈ [1.10, 1.95]` for *either* channel.** The contour sum diverges in *both*
> channels, for the *same* root cause: **no small parameter at `α_s~O(1)`.** **[NG]** for the `0`-form
> Peierls/polymer class at the crossover; **[O]** the exact step = a rigorous `σ_w(2.5) > c₀`. The
> roofing move *removes one obstacle (circularity) and isolates the other (small parameter)* — pin
> progress, **not** closure. **Clay 0%.** *(The answer to the directive's either/or: neither "same
> circularity" nor "genuine escape" — a genuine escape from circularity that lands on the same wall.)*

---

## §0 · Labels, iron laws, acid

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**IL-i** non-positivity (must fail on free massless) · **IL-ii** nonabelian · **IL-iii** no small parameter.
**Acid**: the inequality must be *false* on the free massless field (must not spuriously floor a gapless field).

---

## §1 · The 0-form contour is a different object (⟹ not the same circularity)

| | 1-form (Wilson loop) | 0-form (`tr F²`, `0⁺⁺`) |
|---|---|---|
| defect | magnetic-string worldsheet | domain wall separating 0-form regions |
| geometry | **2-surface** | **3-surface** (codim‑1 in 4D) |
| order parameter | area law (perimeter/area) | **decay of `⟨tr F²(x) tr F²(0)⟩_c`** |
| tension | string tension `σ_s` | **domain-wall tension `σ_w`** |
| logical target | confinement | **the mass gap `m_{0⁺⁺}`** (transfer-matrix spectral theory) |

The 0-form Peierls **never invokes `σ_s`**; proving `σ_w>0` is *not* proving confinement. So it is **not
circular** in the base round's sense (`σ_s>0 ⇔ confinement`). **This is a real structural improvement:**
it removes the *area-law → gap* bridge that Cao–Nissim–Sheffield (and every 1-form argument) cannot cross,
because it targets the gap directly.

> **Honest caveat.** Elitzur forbids a *local* 0-form order parameter, so the 0-form "domain wall" must be
> built from a non-local / dual defect representation of the `tr F²` correlator — **less canonical and
> harder to define rigorously** than the Wilson-loop vortex. The object is different; it is also less
> well-defined. (This is the same Elitzur fact that the base round used to explain the eight SLOT no-gos.)

---

## §2 · The Peierls criterion and the numbers at `β=2.5`

Peierls/polymer convergence criterion:

```
σ(β) := −log(activity)  >  c := log(branching),
```

i.e. the *surface tension* must exceed the *entropic* constant (log of the number of ways a contour can
grow). Witness `w7/roof_peierls_0form.py` (SU(2) single-plaquette, weight `exp(β cosθ)`, coefficients
`a_r/a₀ = (2r+1) I_{2r+1}(β)/I₁(β)`; the numerical `β` is this witness's convention — the *structural*
conclusion is convention-independent):

| quantity | value at `β=2.5` |
|---|---|
| defect activity `a₁/a₀` | **0.5655** |
| `σ = −log(activity)` | **0.5701** |
| 1-form threshold `c₁ = log 5` (each link meets 6 plaquettes → 5 continuations) | 1.6094 |
| 0-form threshold `c₀ ∈ [log 3, log 7]` (3-surface entropy) | [1.0986, 1.9459] |
| **1-form: `σ − c₁`** | **−1.039 → FAILS** |
| **0-form: `σ` vs `c₀`** | **0.570 < 1.099 → FAILS** |

**Validity window** (where each Peierls argument proves the gap):

```
 beta    a1/a0    sigma    1form(>log5)   0form(>log7)
 0.020   0.0000   9.9035     ok             ok
 1.053   0.1297   2.0423     ok             ok
 1.260   0.1806   1.7112     ok             fail
 1.467   0.2373   1.4386     fail           fail
 2.500   0.5655   0.5701     fail           fail

 1-form Peierls validity ends at beta* = 1.333  (rigorous gap only for beta < beta*)
 crossover target beta = 2.5 : OUTSIDE -> Peierls cannot reach it
```

**Both channels stop ~1 order of magnitude in activity short of `β=2.5`.** Constant-optimization can lower
`c₀` to at best `log 3 ≈ 1.10` — still **~2× above `σ(2.5)=0.57`** — so **no tuning of constants closes it.**
(This is a *sufficient*-condition bound; even granting a generous constant, `σ(2.5)=0.57` is too small.)

---

## §3 · Is the 0-form defect free-energy cost controllable at `2.5`?

The 0-form "order parameter" would be a set of 3-surface domain walls; its free-energy cost per unit
3-volume is `σ_w(β)`. It is controllable at `β=2.5` iff `σ_w(2.5) > c₀`. **This is not the string tension
(so not circular), but it is still an `O(1)`-activity bound:** at `β=2.5` the activity is `0.565` and no
rigorous lower bound `σ_w(2.5) > 1.1` exists — indeed the *definition* of the crossover is that the
activity is `O(1)`. The only *genuine* small parameters available to the 0-form channel are:

- **`1/N`** (large‑N): asymptotic, and *not* a small parameter at fixed `N` on the crossover;
- **the UV** (asymptotic freedom): controls the UV, not the IR.

**Neither is a small parameter at `α_s~O(1)`.** Hence: **not controllable at 2.5.**

---

## §4 · What is genuinely improved, and what is not

- **Improved (real).** The 0-form Peierls removes the **area-law → gap bridge** (it targets the gap
  directly) and removes the **`σ_s>0 ⇔ confinement` circularity**. The roof's target is now *well-posed*:
  prove `σ_w(2.5) > c₀` for the 0-form domain wall. This is a strictly cleaner statement than the 1-form
  version, which had to detour through the string tension.
- **Not improved.** It does **not** supply the small parameter. The same `O(1)` activity that kills every
  *absolute* method (`BASE-PHASE-RIGOROUS`: cluster `13×` short, Dobrushin an order of magnitude,
  Pirogov–Sinai wrong object) kills this one too: `σ(2.5)=0.57` is below every `O(1)` threshold, and the
  validity window ends at `β≈1.3–1.5 ≪ 2.5`.

---

## §5 · Verdict + exact step

- **[T]** the 0-form target is well-posed (`0⁺⁺` correlator decay ⇔ gap) and the object is genuinely
  different from the 1-form (3-surface / domain-wall / `σ_w`, **not** 2-surface / string / `σ_s`).
- **[NG]** `0`-form Peierls/polymer at the crossover: at `β=2.5`, `σ=0.57 < c ∈ [1.10, 1.95]`; the
  contour sum diverges. **Not circular, but equally blocked.**
- **[O]** the exact step = a rigorous small-parameter control of the `0`-form contour/domain-wall activity
  at `β≈2.5` (equivalently `σ_w(2.5) > c₀`) — the **same "no small parameter" wall as
  `BASE-PHASE-RIGOROUS`**, relocated to the 0-form channel.
- **Net for the roof.** One obstacle removed (circularity), one isolated (small parameter). The roof's
  requirement is now *exactly* **"small-parameter control of the 0-form domain-wall activity at the
  crossover"** — pin progress, not closure.

**Iron-law / acid check.** IL-i: on the free massless field the activity is `O(1)` ⇒ the inequality
`σ>c` is *false* ⇒ Peierls gives no gap there (it does **not** spuriously floor a gapless field — acid
passes). IL-ii: the counting is nonabelian-**agnostic** (U(1) has the same vortices) ⇒ abelian-blind, as
in the base round. IL-iii: **fails** at the crossover — that *is* the exact step. **Clay 0%; no QED.**

---

*Refs:* `BASE-PHASE-RIGOROUS.md` (absolute methods stop short), `BASE-MINIMAL-MODEL.md` (magnetic defect
condensation in 3D Z₂; channel transfer to 4D blocked), `BASE-DUALITY.md` (character expansion; no small
parameter at `β=2.5`), `BASE-SYNTHESIS.md` §2/§6; Cao–Nissim–Sheffield (area law ≠ gap); witness
`w7/roof_peierls_0form.py`.

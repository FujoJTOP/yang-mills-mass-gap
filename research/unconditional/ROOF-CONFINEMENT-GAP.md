# ROOF-CONFINEMENT-GAP — does the 1-form area law imply the 0-form `0⁺⁺` gap?

> **VERDICT (no QED).** **Not kinematic. Not equivalent. True and provable at strong coupling; open at the crossover.** The `0⁺⁺` operators (`tr F²`) are local and center-neutral ⇒ they sit in the **trivial sector** of the center 1-form symmetry, while the Wilson/'t Hooft area law is a statement about **1-form-charged** extended operators; the two sectors are orthogonal, so **no kinematic implication exists in either direction** — any bridge is *dynamical* `[T]`. In the opposite direction **gap ⇏ confinement** (high-`T` deconfined pure YM is gapped but has perimeter law) `[T]`. In the direction we need — **confinement ⟹ gap** — the bridge is **TRUE and PROVEN at strong coupling**: area law (`σ>0`, Osterwalder–Seiler) and the `0⁺⁺` gap (SZZ, `β<1/48`) come from the **same polymer expansion**, yielding `m_{0⁺⁺} ≥ c√σ` with `c>0` `[C]`. It **fails to extend to the crossover** — same KP/polymer radius (`β* ≤ 8.9×10⁻³`, ≈2 orders below `2.5`), i.e. exactly `H1`. No model with (area law ∧ gapless `0⁺⁺`) is known. **Clay 0%.**

---

## §0 · Labels, the two properties, acid

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class.
- **(A) Confinement / area law**: the center-electric Wilson loop `W(C)` obeys `⟨W(C)⟩ ~ e^{−σ·Area}`, `σ>0`; equivalently the center 1-form symmetry is unbroken and the dual ('t Hooft) sector orders.
- **(B) `0⁺⁺` gap**: `Δ := inf spec(H) − 0 > 0` in the **vacuum / trivial** sector, with the lightest state having `f₁² = |⟨0|tr F²|G⟩|² > 0` and isolation (`H*`).
**Acid**: any claimed bridge must *fail* on the massless free field (no area law, no gap) and must be non-abelian-essential.

---

## §1 · `[T]` The kinematic no-go — why the bridge cannot be automatic

Let `R` be the center 1-form symmetry of pure `SU(N)` (electric) and `R*` its magnetic dual. The
color-singlet local operator
```
O_{0⁺⁺}(x) = tr F_{μν}F^{μν}(x)      (and tr F∧F for 0⁻⁺)
```
is **local** and **center-neutral**: it carries **trivial charge under `R`**. Hence the entire `0⁺⁺`
spectrum lives in the **`R`-trivial sector** `H_{q=0}`. In contrast:
```
W(C)  is R-charged (electric),       B(C) ('t Hooft)  is R*-charged (magnetic).
```
The area law is a statement about expectation values of these **charged, extended** operators — i.e. the
**1-form-charged sector**. The Hilbert space decomposes `H = ⊕_q H_q` by 1-form charge, and `H_{q=0}` is
**orthogonal** to every `H_{q≠0}`. **A statement about the charged sector therefore places no kinematic
constraint on the trivial sector.** ⇒

> **`[T]` "area law ⟹ `0⁺⁺` gap" is not a theorem of the definitions.** Any bridge must be *dynamical*
> (a transfer of order from the `q≠0` sector to `q=0`). **This is exactly the "channel mismatch" the base
> round (`BASE-ORDER-PARAMETER`) measured**, now stated as a structural fact.

**Consequence for the campaign:** the 't Hooft / center-vortex order parameter *cannot* be the `0⁺⁺` gap
detector (different sector). The `0⁺⁺` gap must be re-derived in the trivial sector — consistent with the
base round's finding that the useful object is a `0⁺⁺` anomaly (`⟨tr F²⟩`), not a 1-form order parameter.

---

## §2 · `[T]` The other direction is cleanly false: gap ⇏ confinement

At **high temperature** pure `SU(N)` YM is in the **deconfined** phase: `⟨W(C)⟩` obeys a **perimeter** law
(no string tension, `σ=0`), while the **spatial** transfer matrix is **gapped** — the screening masses
(electric `m_E ~ gT`, magnetic `m_M ~ g²T`) are positive. So

> **`[T]` (at the level of the high-`T` gapped deconfined phase) mass gap does NOT imply confinement.**

Hence the tempting "elegant equivalence" (confinement ⟺ mass gap) is **false**: the forward direction
(confinement ⟹ gap) is what is at issue, and the reverse is demonstrably not general.

---

## §3 · `[C]` The bridge that actually exists — strong coupling

At strong coupling both (A) and (B) are proven **from the same polymer / character expansion**, and one
can read off a *quantitative* bridge.

- **Area law `[T]`** (Osterwalder–Seiler; Wilson strong-coupling expansion): `σ a² = −ln u + O(u) > 0` for
  small `β` (`u = u(β)` the fundamental character weight).
- **`0⁺⁺` gap `[T]`/`[C]`** (P3, `PROOF-P3-UNCONDITIONAL.md`; SZZ `β<1/48`): the lightest `0⁺⁺` residue is
  a **four-walled tube** of length `t`, giving `m₁ = 4 ln(1/u)/a` and `f₁² ≥ r₀(β)>0`, `L`-uniform.

Combining the **same** `u`:
```
m₁ / √σ  =  [4 ln(1/u)/a] / [√(ln(1/u))/a]  =  4√(ln(1/u))  ≥  4   for  u ≤ e^{−1} (strong coupling).
```
⇒

> **`[C]` (strong-coupling bridge) `m_{0⁺⁺} ≥ 4√σ · √(ln(1/u)) ≥ 4√σ > 0`.** In the strong-coupling regime,
> **confinement (σ>0) implies the `0⁺⁺` gap, provably, and the constant is explicit.** The area law and the
> gap are two readings of the **same** polymer weight `u` — this is *why* they track each other physically.

This is genuine: at strong coupling the implication is **true and proven**, and the "bridge" is a concrete,
computable inequality (`m₁ ≳ 4√σ`), not a slogan. It also shows the *mechanism* of the transfer: the
`0⁺⁺` tube of length `t` is built from the **same** plaquette weight that gives the area law — the
`q=0` and `q≠0` sectors are correlated through `u`, not kinematically.

---

## §4 · `[O]` Where the bridge breaks — the same wall

The strong-coupling bridge is only as good as the polymer radius. Both `σ a²=−ln u+O(u)` and
`m₁=4 ln(1/u)/a` (with `f₁²≥r₀`) are controlled only while the expansion converges:
```
β*(N) ≤ 8.9×10⁻³   (optimised grading; 3.36×10⁻³ crude KP; benchmark 1/48 ≈ 2.08×10⁻²)
```
— **≈2 orders of magnitude below the crossover `β≈2.5`** (`BASE-INFINITE-VOL`). Past the KP boundary
`λ_aΦ(a)=½` the activity sum diverges and no `β`-uniform inequality replaces it `[NG]` for the class.
⇒

> **`[O]` The bridge "confinement ⟹ `0⁺⁺` gap" exists and is proven at strong coupling, and its missing
> span is exactly the IR crossover `α_s~O(1)` = `H1`.** The roof is a *two-ended* structure (strong coupling:
> bridge proven; `β→∞`: Balaban UV control) with the **same notch** in the middle.

---

## §5 · Counterexample search (for a no-go, in the confined phase)

A *no-go* for "confinement ⟹ gap" would need a model with **area law ∧ gapless `0⁺⁺`**. Status:

- **No such model is known.** Every solvable confining model is *also* gapped: 2D `SU(N)` YM (exact:
  area law + discrete gapped spectrum); 3D `Z₂` gauge ≡ dual 3D Ising (Peierls area law + exponential
  clustering `[T]`). Confinement and gap co-occur in all of them.
- **'t Hooft-anomaly gaplessness at `θ=π` is not a counterexample.** `CP¹`/`O(3)` at `θ=π` is forced *not*
  trivially gapped by the anomaly — but there the would-be **1-form symmetry is spontaneously broken**
  (deconfined), so the premise "confinement (area law)" **fails** at `θ=π`; the anomaly forces the theory
  out of the confined phase, not a gapless confined phase. This shows the bridge is **not universal**
  (anomalies can break it) but does **not** produce (area law ∧ gapless `0⁺⁺`) in 4D-like pure gauge.
- **Holographic "confinement without mass gap"** exists only in *non-QFT / tuned* soft-wall backgrounds and
  is not a counterexample for a local 4D gauge theory.

⇒ **`[NG]` for the "separate them in the confined phase" hope:** we found **no** rigorous model with area
law but a gapless `0⁺⁺`; the implication, where testable, holds. But (per §1) it is **not kinematic**, so
its 4D truth is dynamical = the wall, not a theorem.

---

## §6 · Verdict

| question | answer |
|---|---|
| Is "area law ⟹ `0⁺⁺` gap" kinematic? | **No `[T]`** — orthogonal 1-form charge sectors (`0⁺⁺` is trivial, area law is charged). |
| Is "area law ⟺ gap" an equivalence? | **No `[T]`** — gap ⇏ confinement (high-`T` deconfined gapped). |
| Does confinement ⟹ gap **at strong coupling**? | **Yes, provable `[C]`** — `m_{0⁺⁺} ≥ 4√σ`, same polymer weight. |
| Does it carry to the crossover? | **No `[O]`** — same KP radius, `β* ≤ 8.9×10⁻³` ≪ 2.5. |
| Is there a 4D counterexample (area law ∧ gapless `0⁺⁺`)? | **None known `[NG]`** — θ=π anomalies force *deconfinement*, not a gapless confined phase. |

**Net:** the implication is **dynamical, true and proven at strong coupling (with the explicit bridge
`m_{0⁺⁺} ≥ 4√σ`), false in the reverse direction, and open at the crossover = `H1`.** The roof has two
proven ends and the **same central notch** as the base. It does not "solve" the gap; it **locates the roof's
missing span precisely** and shows the transfer is real where controlled. **Clay 0%; no QED.**

*Refs:* `BASE-ORDER-PARAMETER.md` (channel mismatch `[NG]`), `BASE-PHASE-RIGOROUS.md` (two-ended base),
`BASE-INFINITE-VOL.md` (`β*`, `f₁²=W_tube·R_vac`), `PROOF-P3-UNCONDITIONAL.md` (`m₁=4ln(1/u)/a`, `f₁²≥r₀`),
Osterwalder–Seiler 1978 (area law + gap, strong coupling), Shen–Zhu–Zhu CMP 400 (2023) `β<1/48`.

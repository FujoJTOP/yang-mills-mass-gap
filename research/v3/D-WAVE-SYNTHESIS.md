# D-WAVE-SYNTHESIS — decompactification: the constant-mode curvature SATURATES in the strip (c>0, Haar floor 0.2688 exact); the 4D floor remains H1 because the d≥2 susceptibility can diverge

> **VERDICT (no QED, Clay 0%).** Six routes attacked the decompactification question: does the
> constant-mode curvature decay to zero as `L→∞` (then the uniform floor is `H1`), or saturate at
> the `L`-independent Haar floor? **In the strip model it saturates** — every fit gives `c>0`, the
> `β=0` Haar floor is exact and `L`-independent, and the transfer matrix is gapped. **But in `d≥2`
> the constant-mode susceptibility can diverge (Kosterlitz–Thouless/Goldstone), so the 4D uniform
> floor is not settled; it remains `H1`.** The result is a genuine positive theorem in 1D and a
> sharp reduction in 4D, not a proof of the gap.

---

## §0 · The six verdicts

| Route | Verdict | Result |
|---|---|---|
| **D1** exact strip `L=1..16` | `[T]`/`[V]` | Fits `c+A/L²`: `c=0.248` (`θ=0`), `c=0.410` (`θ=π`); `c+A q^L`: `0.272`/`0.416`; every form gives `c>0`. The `θ=0` value brackets the Haar floor `0.269`. Mechanism: the constant mode is a zero mode of the Wilson difference operator, so the one-loop Wilson part is exactly `0`; the floor is Haar pinning. |
| **D3** transfer matrix | `[T]`/`[V]` | Bulk gap `Δ=3.18` (`θ=π`), `0.64` (`θ=0`) at `β=2.5`, `L`-independent, `ξ<∞`; total `G″` falls `1.59→…→0.32` (`L=2..10`) and fits `c+A/L²` with `c=+0.22±0.03`. |
| **D4** Haar/Wilson split | `[T]`/`[V]` | Exact `G″=G″_Haar+G″_Wilson`; Haar part is a deterministic `L`-fold convolution with exact limit `1/Var(ξ₃)=0.26883>0`; Wilson part is positive `O(1)` on strips; total flattens at `≈0.46` (`L=4–6`). |
| **D2** analytic limit | `[T]`/`[C]`/`[V]` | The constant mode is a zero mode of the Wilson Hessian; `Var(a)` plateaus at `≈2.34`, so `G″→1/Var≈0.43` at `β=2.5`; `β`-monotone `0.269→0.31→0.43→0.55`. |
| **D5** universality | `[T]`/`[V]` | The L-dependence is **not universal**. The invariant statement is `R:=G″·Var_ρ(a)→1` (marginal Gaussianizes), `Var_ρ(a)→χ_∞∈(0,∞)`, hence `G″→1/χ_∞>0`. Values depend on group/action/boundary; U(1)'s `β=0` floor is zero, so `N/12` is non-abelian. Saturation is **twist-specific**: at `b=0` the variance grows and `G″` decays. |
| **D6** exact identity | `[T]`/`[V]`/`[C]` | `G″(a₀*)=1/χ`, `χ=Σ_d C(d)`; at `β=0`, `χ=(4π²−6)/9`, `G″→0.2688` exactly; the `1×L` chain has `χ<∞`. At `β=2.5` the fit gives `b≈0.266`. In 4D, `χ<∞` is exactly `H1`. |

---

## §1 · The invariant statement

The only normalization-independent statement the wave establishes is:

> In the strip (1D) model, and at twist boundary `θ=π`, the marginal Gaussianizes,
> `G″·Var_ρ(a)→1`, and `Var_ρ(a)→χ_∞∈(0,∞)`, so
> `G″(a₀*)→1/χ_∞>0`.
> The value of `χ_∞` is **not universal**: it depends on the group, the action, the twist, and
> `β` (SU(2) strip: `0.41`–`0.48`; `β=0` Haar: `0.2688`; `θ=0`: `0.248`).

The `≈1/L` seen in the first small-volume scans is a transient through the `L=2` dimer bump, not
the asymptote.

---

## §2 · What is proved, and what is not

**Proved `[T]`:**
- The constant mode is a zero mode of the Wilson difference operator; the one-loop Wilson part of
  `G″` is exactly zero.
- The Haar part is a deterministic `L`-fold convolution with exact `L`-independent limit
  `1/Var(ξ₃)=0.26883>0`.
- In the 1D strip, extensivity plus finite correlation length imply `χ<∞`, hence a positive
  `L`-independent floor.

**Not proved `[O]`:**
- In `d≥2`, the constant-mode susceptibility can diverge (Kosterlitz–Thouless/Goldstone), so
  `c=0` is possible; the 4D uniform floor is unchanged `=H1`.
- The periodic (untwisted) boundary case: at `b=0` the variance grows and `G″` decays; the physical
  periodic torus is not the twist model.
- The non-constant `k→0` soft tower.
- The identification `1/χ =` the pole/glueball mass is `[C]`; and a positive curvature floor is
  **not** an atom: at `β=0`, `G″=0.2688>0` while `f₁²=0`.

---

## §3 · What changes

- **The decay-to-zero route is dead in 1D.** The constant-mode curvature saturates; the small-volume
  `≈1/L` was a transient. This is the strongest positive decompactification result of the campaign.
- **The 4D question is sharpened:** not "does it decay?" but "does the `d≥2` susceptibility stay
  finite?" — which is `H1`, now with an explicit invariant (`χ<∞`) and a concrete value at `β=0`.
- **The positive floor is not the gap:** `G″>0` at `β=0` coexists with `f₁²=0`. The atom remains
  the missing object.

**Clay = 0%. No QED.**

*Refs:* `D1-EXACT-L.md`, `D2-ANALYTIC-LIMIT.md`, `D3-TRANSFER-CONSTANT.md`, `D4-HAAR-WILSON.md`,
`D5-UNIVERSALITY.md`, `D6-DECISION.md`; witnesses in `w7/`; `SM-WAVE-SYNTHESIS.md`.

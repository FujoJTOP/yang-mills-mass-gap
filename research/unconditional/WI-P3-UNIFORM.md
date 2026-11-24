# WI-P3-UNIFORM — the P3 approach-law's uniformity in `(L,a)` **is** the gap (`H1`); it is a finite-`L` diagnostic, not an independent route to `I2`

> **VERDICT (decisive).** `WI-P3-CONTINUUM` showed the approach-law is a robust, acid-clean isolator *given* `Δ=E_thr−e₁>0`. This file attacks the **uniformity** `Δ(L,a) ≥ δ > 0` across `(L,a)` at β≈2.5 (the step that would make it an `I2` certifier). **The uniformity is exactly `m₁ > 0`.** The isolation gap of the `0⁺⁺` channel is `Δ(L) = m₁ + O(L⁻²)` for an **atom** but `Δ(L) = O(L⁻²) → 0` for a **cut** (`:φ²:`), so `Δ(L,a) ≥ δ` uniformly **⟺** the floor is an isolated atom **⟺** the gap. Making P3 uniform is *proving the gap*, not a new handle. `[T]/[O]`, Clay 0%, no QED.

## §1 · The finite-size-scaling law `[T]`

On `T³_L`, the `0⁺⁺` spectrum is discrete, `E₁(L) < E₂(L) < …`.
- **Atom** (lightest `0⁺⁺` = glueball at `m₁`, below the two-glueball threshold `2m₁`): `E₂(L)=min(2m₁, m_{0⁺⁺*})` up to `O(L⁻²)`, so `Δ(L)=m₁+O(L⁻²) → m₁>0`; and as `a→0`, `m₁(a)→m₁>0`. **Uniform in `(L,a)`.**
- **Cut** (`:φ²:`: floor `2m>0` but a two-particle continuum, `f₁²=0`): `E₁=2m` (`p_rel=0`), `E₂=2√(m²+(2π/L)²) ≈ 2m+(2π/L)²/m`, so `Δ(L)=(2π/L)²/m → 0`. **No uniform positive bound.**

Witness (`w7/wi_p3_uniform.py`): atom `Δ=1.000000` at every `L`; cut `Δ=1.72, 0.543, 0.149, …` with `Δ·L² → (2π)²=39.48`. A **positive floor is strictly weaker than uniform isolation** — they agree only for an atom.

## §2 · Why P3 is only a finite-`L` diagnostic `[O]`

At fixed `(L,a)` the transfer matrix is compact, so `Δ>0` is **trivial** and the approach-law "works" **vacuously** (same finite-`L` triviality as `FALSIFY-FIELD`'s finite-`L` atom). The **content** is the `L→∞`, `a→0` survival of `Δ≥δ`, which by §1 is `m₁>0`. All three uniformity directions — in `L`, in `a`, and across `β` (`1/48 → 2.5`) — reduce to `H1`. **No monotonicity gives it for free**: `Δ(L)=m₁+c/(m₁L²)` only yields a floor if `m₁>0` is already known; the cut's `Δ→0`.

## §3 · Genuine yield `[C]`

A **concrete, checkable FSS law**: `Δ(L)=m₁+c/(m₁L²)` (two-particle threshold) or `m_{0⁺⁺*}−m₁` (excited glueball). A lattice measurement of `m_O(τ)` that shows **exponential + plateau** (atom) vs **`1/L²` decay** (cut) is a real finite-`L` isolation signature — but its `L→∞` verdict is `H1`.

**Decisive answer:** P3 is a **finite-`L` diagnostic**; uniformising it *is* proving the gap (`H1`). Not a route to `I2` independent of `H1`.

*Refs:* `WI-PLATEAU-RIGOR.md`, `WI-P3-CONTINUUM.md` (adjacent: continuum-robustness), `WI-EMBEDDED.md`, `FALSIFY-FIELD.md`.

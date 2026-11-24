# INVENT-FLUX — the center-flux twisted-gap functional: the 0⁺⁺ gap is **exactly blind** to the center flux by 1-form neutrality; the flux response is O(1/L) (Casimir / momentum shift), and the only flux-sensitive object is the 't Hooft loop = 1-form — **it inherits the channel wall**

> **VERDICT (no QED).** The invented species `Γ(𝔷)` — the 0⁺⁺ spectral gap as a function of a
> threaded discrete center flux — is **not a new handle.** **`[T]`** The 0⁺⁺ operator `tr F²` is
> center-neutral (zero N-ality), hence blind to the 1-form center-flux superselection sector; its
> two-point function — and therefore the 0⁺⁺ mass — is **sector-independent**. The residual flux
> dependence is **O(1/L)** (a finite-volume Casimir / momentum-shift effect that vanishes as `L→∞`),
> and the only genuinely flux-sensitive object is the **'t Hooft loop**, which is **1-form** and
> already channel-mismatched. So `INVENT-FLUX` **collapses onto the center-neutrality / channel-mismatch
> wall**: it is the 't Hooft loop in spectral clothing, not a 0⁺⁺ probe. Failing step pinned. **Clay 0%.**

---

## §0 · Labels, iron laws, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric. **IL-i** non-positivity (fail on free massless) · **IL-ii** nonabelian · **IL-iii** no small parameter. **Acid**: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The species

Thread a discrete center flux (twist) `𝔷 ∈ Z_N` through a spatial cycle of `T³_L`. **Define** the twisted
0⁺⁺ gap
```
Γ(𝔷) := m_{0⁺⁺}(𝔷) = −lim_{t→∞} (1/t) log ⟨0_𝔷| tr F²(x,t) tr F²(x,0) |0_𝔷⟩,   |0_𝔷⟩ = ground state in flux sector 𝔷.
```
The **flux response** is `δΓ(𝔷) = Γ(𝔷) − Γ(1)` (or `dΓ/d𝔷` in the continuous `α` relaxation). The hope: a
**global, non-order-parameter, flux-sensitive** functional that is inhomogeneous at `β=2.5` even where the
untwisted gap is not (since B2: no order parameter, and the untwisted gap is not floored).

---

## §2 · The theorem: the 0⁺⁺ is exactly flux-blind `[T]`

Pure `SU(N)` at `θ=0` has an **exact 1-form center symmetry** `G^{(1)} = Z_N` (non-anomalous at `θ=0`).
On `T³` the Hilbert space decomposes into **superselection sectors** labeled by the 't Hooft flux
`𝔷 ∈ H²(T³;Z_N) = (Z_N)³`, and the 1-form symmetry **maps sectors to sectors unitarily** while
**commuting with every neutral local operator** (i.e. every operator of zero 1-form charge). The 0⁺⁺
operator `tr F²` is:
- a **color singlet** (adjoint `Ad_{zI} = id ⇒` invariant, verified `|Ad_{zI}(X)−X| = 2.2×10⁻¹⁶` in
  `ROOF-0FORM-DUAL`), and
- **center-neutral** (zero N-ality, zero 1-form charge).

Hence `[tr F², U_𝔷] = 0` for the symmetry relating sectors, so
```
⟨0_𝔷| tr F²(x) tr F²(y) |0_𝔷⟩  is independent of 𝔷,   ⇒   Γ(𝔷) = Γ(1)   (in the mass). `[T]`
```
**The 0⁺⁺ gap cannot see the center flux.** Any `𝔷`-dependence of the *ground-state energy* `E_𝔷 − E_1` is a
**Casimir/finite-volume term** of order `1/L` (below), which **cancels in the mass** and **vanishes** as
`L→∞`. This is the same center-neutrality that killed the 0-form 't Hooft algebra (`ROOF-0FORM-DUAL`),
now shown to make the flux-response species vacuous.

---

## §3 · Acid test and witness (`w7/invent_flux_witness.py`, exact)

**(1) Local gap vs twist, massless `m=0` (acid: must → 0 as `L→∞`):** the twist lifts the zero mode only
at finite `L`, giving `gap(α) = π·dist(α,2πℤ)/L`:

| L | gap(α=0) | gap(α=π/2) |
|---|---|---|
| 4 | 0.00000 | 0.39270 |
| 8 | 0.00000 | 0.19635 |
| 16 | 0.00000 | 0.09817 |
| 32 | 0.00000 | 0.04909 |
| 64 | 0.00000 | 0.02454 |

→ `π/(2L) → 0`: the twist effect is **pure finite volume**; the physical (L→∞) gap is 0 for every α. **Acid passes** (returns 0).

**(2) Local gap vs twist, massive `m=1` (gapped):** → `m=1` for every α as `L→∞`
(`1.074, 1.019, 1.005, 1.001, 1.0003` at `L=4…64`). **The thermodynamic gap is twist-independent.**

**(3) Twisted partition function (Casimir) is α-dependent, but only at O(1/L):** `logZ` at `L=8` varies by
`~6×10⁻⁴` across α; at `L=16` it is flat to `6` digits. So even the *global* flux response is a
**finite-volume Casimir effect** that dies as `L→∞` — not an infinite-volume handle.

**(4) Neutral vs charged sector-independence (1-form toy, exact):** a neutral operator has phase `1` in
**every** flux sector `𝔷`; a charged (q=1) operator has phase `e^{2πi𝔷/N}` (`N=5`: sectors give
`1, e^{2πi/5}, …`). **The neutral 0⁺⁺ cannot distinguish sectors; only 1-form operators can.**

---

## §4 · What *is* flux-sensitive — and why it is not a 0⁺⁺ handle

The genuinely flux-sensitive objects are the **'t Hooft loop** `B(C)` (creates center flux through `C`) and
the **twisted partition function** `Z(𝔷)`:
- `B(C)` is the **1-form order parameter** (confinement), **Schur-orthogonal to the neutral 0⁺⁺** — exactly
  the channel mismatch proved in `ROOF-CONFINEMENT-GAP` / `ROOF-ORDER-PARAMETER`.
- `Z(𝔷)`'s `𝔷`-dependence is a **Casimir/ground-state term** `O(1/L)` (witness (3)), which cancels in the
  mass and vanishes at `L→∞`.

So "flux response" = the 't Hooft loop (1-form) + a finite-volume Casimir — **neither is a 0⁺⁺ gap floor.**

---

## §5 · Loopholes checked (all closed)

1. **Anomalous 1-form symmetry (θ=π):** would break center-neutrality and allow a sector dependence — but at
   `θ=π` the 1-form symmetry is spontaneously broken (deconfinement), the confined-phase premise fails, and
   `θ=0` is the physical theory. `[NG]`
2. **Finite-`L` `O(1/L)` mixing with the mass:** vanishes as `L→∞`; not an infinite-volume handle. `[T]`
3. **Mixed correlator `⟨B(C) tr F² tr F²⟩` (0⁺⁺ around a flux string):** a codim-2 defect; the leading 0⁺⁺
   mass is unchanged, corrections are power-law — no floored gap. `[C]`
4. **A non-neutral probe:** that is the 't Hooft loop = 1-form = channel-mismatched. `[NG]`
5. **Volume-reduced / twisted-EK large-N:** the twist stabilizes the center and the 0⁺⁺ mass becomes a
   matrix-model quantity — but this is the **large-N limit**, not a flux handle at `α_s~O(1)` fixed `N`
   (see `ROOF-LARGEN-MASTER`); it does not floor the finite-`N` crossover. `[O]`

---

## §6 · Verdict and exact failing step

**`INVENT-FLUX` fails — it inherits the wall.** The exact failing step is **[T] center-neutrality**:
`tr F²` has zero 1-form charge, so `Γ(𝔷) = Γ(1)` in the mass, and the residual flux response is `O(1/L)`
Casimir. The only flux-sensitive object is the **'t Hooft loop (1-form)**, already channel-mismatched.
The species is therefore **the 't Hooft loop in spectral clothing**, not a new 0⁺⁺ handle — it lands on the
same center-neutrality / channel-mismatch obstruction as `ROOF-0FORM-DUAL` and `ROOF-ORDER-PARAMETER`.

**Iron laws:** IL-i passes numerically **vacuously** (returns the gapless branch); IL-ii PASS (center
symmetry is nonabelian-specific); IL-iii **FAIL** (no small parameter: the surviving flux response is the
1-form/'t Hooft wall).

---

## §7 · Honest bottom line

- **The twisted-gap functional is well-defined and global, but it is center-blind and finite-volume: `[T]`.**
- **No new 0⁺⁺ handle.** The flux response is the 't Hooft loop (1-form) plus an `O(1/L)` Casimir.
- The strut (`χ_lattice ≲ a²` / `σ_w(2.5)>c₀` / `R_vac` monotone) is **untouched** by this species.
- **Clay = 0%; no QED.** This is a clean `[NG]+[T]` no-go with a named exact step, as the round requires.

*Refs:* `ROOF-0FORM-DUAL.md` (center-neutrality), `ROOF-CONFINEMENT-GAP.md` (1-form orthogonality),
`ROOF-ORDER-PARAMETER.md` (channel landscape), `ROOF-SYNTHESIS.md`; witness `w7/invent_flux_witness.py`.

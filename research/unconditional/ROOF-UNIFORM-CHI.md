# ROOF-UNIFORM-CHI — the roof needs ONE strut: a uniform upper bound on `χ_phys` at β≈2.5; every route to it is either channel-mismatched, small-β-absolute, or H1 verbatim

> **VERDICT (no QED).** The base gives `γ ≥ (1/2π)(1/χ)` on the `k→0` mode, calibrated on Gaussian/`Z_n`.
> The roof needs `1/√χ ≥ λ` **uniformly at β≈2.5**, i.e. an **upper bound on the physical coarse-mode
> variance `χ_phys`**. Five candidate routes (trivial boundedness, FKG/positive-association on the dual,
> Gaussian domination on χ, Haar/Jacobian positivity `Φ″≥N/12`, monotone/LSI) are evaluated. **None gives a
> uniform floor.** The **minimal missing input is exactly `χ_phys ≤ 1/λ` uniformly**, which **is `H1`
> verbatim**. There is **no route** to a uniform `1/χ` floor at `α_s~O(1)` that is not `H1`. **Clay 0%.**

---

## §0 · Labels, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric.
**Acid**: `γ=0` / `χ=∞` on the free massless field. Iron laws: non-positivity, nonabelian, no small parameter.

---

## §1 · The precise object, and why the scaling is the whole game

Cheeger (base): `γ ≥ (1/2π)(1/χ)`, `χ = Var_ν(a)` on the `k→0` 0⁺⁺ mode, `1/χ = G″(k→0)` (Ward `G″χ≈1`).
The roof needs `1/χ ≥ λ`, i.e.
```
χ_phys(β=2.5, L→∞) ≤ 1/λ ,   uniformly in (L,a).                          (★)
```
**Scaling is decisive.** `a = ⟨v, X⟩`, `X=logU` (dimensionless, `|X|≤π√(2(N²−1))`); the physical mode is
`a_phys = X/a_ℓ` (`a_ℓ` = lattice spacing), so `χ_phys = χ_lattice / a_ℓ²`. Hence `(★)` reads
```
χ_phys ≤ 1/λ   ⟺   χ_lattice ≤ a_ℓ²/λ   (an UPPER bound that SHRINKS like a_ℓ² as a_ℓ→0).
```
So the roof is **not** "χ_lattice bounded" (that is trivial and useless — see R1); it is
**"χ_lattice ≲ a_ℓ² non-perturbatively at β=2.5"** = the weak-coupling scaling = **dimensional transmutation**.

---

## §2 · The five routes, and the exact reason each dies

**R1 — Trivial boundedness.** `U∈SU(N)` compact ⇒ `|X|`, hence `|a|`, bounded ⇒ `χ_lattice ≤ A²` uniformly.
**Dies (acid):** this is a **kinematic** (compactness) bound. Uniform-on-`[−A,A]` witness: gap `γ=π²/A²`,
`χ=A²/3`; as `A→∞` (the `a_ℓ→0` limit) `γ→0`, `χ→∞`. Cheeger holds at every `A` but the floor `1/χ→0`.
**A bounded state space gives a lattice gap that vanishes in the continuum — never a physical floor.** `[NG]`

**R2 — FKG / positive association on the dual.** **Dies (channel mismatch + absent hypothesis):** FKG needs a
ferromagnetic (monotone) coupling; the full-representation Wilson measure is **not FKG** (BASE-INFINITE-VOL),
and the dual FKG that *does* exist acts on the center / 1-form sector, which is **Schur-orthogonal to the
0-form `a`** (BASE-ORDER-PARAMETER). It bounds the wrong observable. `[NG]`

**R3 — Gaussian domination (Fröhlich–Simon–Spencer) on χ.** FSS upper-bounds correlations, so in principle
`χ = Σ_x C(x) ≤ Σ_x (Gaussian bound)`. **Dies (small-β absolute):** the FSS rate is useful only at small `β`;
at `β=2.5` it degenerates and gives no useful decay (the same `β*` wall as BASE-INFINITE-VOL, `β*≤8.9×10⁻³`).
`[NG]`

**R4 — Haar/Jacobian positivity `Φ″≥N/12` (Prékopa).** Real and pointwise ⇒ the Haar part of `−log p` is
**strictly convex**. **Dies (acid-vacuous in the continuum):** `Φ` is a lattice-scale Jacobian; in the
continuum **linearization** — where the free massless field lives — `Φ→0` (`W_geo=0` on the acid field,
BASE-SD-FIXPOINT O3; Theorem BE). So it produces a lattice-scale curvature that **does not survive** into a
physical mass unit. It cannot give a uniform-in-`a` floor. (Also: the coarse 0⁺⁺ mode may be *convex alone*,
but "convex alone" ⇒ `κ≥0`, not `κ≥λ`; and `≥N/12` in lattice units `→0` physically.) `[NG]`

**R5 — Monotone quantum inequality / LSI.** An LSI with constant `c` ⇒ `χ ≤ 1/c`; but on the `k→0` mode the
LSI constant **is** the gap `≈κ` (circular), and there is **no proven monotone in β** (numerically the residue
*dissolves* with β, BASE-INFINITE-VOL). `[NG]`

---

## §3 · The minimal missing input (exact)

Beyond the base's `1/(2π)` fraction, the roof needs **exactly one strut**:
```
(★)  χ_phys(β=2.5, L→∞) ≤ 1/λ   uniformly,   ⟺   G″(k→0) ≥ λ uniformly.
```
This is **`H1` verbatim** (the coarse-mode Poincaré / inverse-susceptibility floor). No smaller or
"cleaner" object suffices: the Cheeger fraction already extracts everything the *shape* of the bound gives;
what remains is the *magnitude* of `1/χ`, and that magnitude is `H1`.

---

## §4 · Is there ANY route to a uniform `1/χ` floor that is not `H1`?

**No — and this is structural, not an accident of the five tried routes.** A uniform `1/χ` floor at
`α_s~O(1)` is the statement that a **mass scale is generated non-perturbatively** (dimensional
transmutation). Any hypothesis strong enough to yield it at the crossover already encodes
**clustering / the gap**; any hypothesis weaker than that reproduces one of:
- **channel-mismatch** (R2: controls 1-form/dual, not 0-form `a`);
- **small-β absoluteness** (R3, and all of BASE-PHASE-RIGOROUS: convergence radius stops short of 2.5);
- **kinematic/jacobian artifacts that vanish in the continuum** (R1, R4: acid-vacuous);
- **circularity** (R5: LSI = gap).

These three failure modes are exhaustive for lattice-gauge positivity tools, so **no sixth route escapes.**
(The only exceptions would be a genuinely new non-perturbative 4D-gauge RG — which does not exist and *is*
the Clay-scale object.)

---

## §5 · What CAN be built: the conditional roof

The roof is **complete modulo one certificate**. Given Cheeger `[T]` plus a **certified** `χ`:
```
γ ≥ (1/2π)·(1/χ_certified)    rigorously.                                  (ROOF|χ)
```
At `β=2.5` with the finite-volume measured `1/χ ≈ G″ ≈ 0.28`, this yields `γ ≥ (1/2π)(0.28) ≈ 0.045` — a
**rigorous 4.5% fraction of the measured inverse susceptibility**, *conditional on certifying* `χ ≤ 1/λ`
on an `L→∞` sequence. The **only missing strut is the χ certificate**; every other piece of the roof
(Cheeger, calibration, Ward `G″χ≈1`) is `[T]`.

---

## §6 · Honest bottom line

- **Roof not laid.** Five routes, five `[NG]`; no uniform `1/√χ` floor at β≈2.5.
- **Minimal missing input named exactly:** `χ_phys ≤ 1/λ` uniform `= H1`. The roof stands **iff** that one
  certificate exists.
- **The base's `1/(2π)` fraction is genuine but magnitude-blind** — it converts a certified `χ` into a
  rigorous gap fraction; it cannot manufacture the `χ`.
- **No `QED`. Clay 0%.** The roof is a **conditional** structure: complete given the `H1` certificate,
  empty without it.

*Refs:* `BASE-NONLOCAL-BOUND.md` (Cheeger, calibration), `BASE-INFINITE-VOL.md` (`β*`, no monotone),
`BASE-ORDER-PARAMETER.md` (channel mismatch), `BASE-FIELDMEASURE.md` (Theorem BE, Φ-vanishing),
`BASE-SD-FIXPOINT.md` (O3: `W_geo=0` on acid field); witness `w7/roof_uniform_chi.py`.

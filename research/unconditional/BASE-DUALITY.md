# BASE-DUALITY — electric–magnetic / character duality as a foundation: the right non-local order parameter, but **no dual small parameter at the crossover**

> **VERDICT (no QED).** Duality supplies a genuine **BASE** piece — it fixes the *type* of the global
> object (a non-local loop/surface order parameter, not a local observable) and the only convergent dual
> coordinate (the character / Fourier-on-the-group expansion). It does **not** supply a **KEY**: at the
> crossover `β≈2.5` the dual expansion has **no small parameter** (character coefficients `a_r/a_0 → 1`),
> and the 't Hooft-loop / center-vortex area-perimeter law *is* the confinement statement. **The difficulty
> does not transfer to a controlled dual side — it stays put = `H1`.** A true small-parameter dual exists
> only on a **deformed small circle** (`Unsal–Yaffe`), and decompactifying it is the same transport wall.
> **Labels:** `[T]` loop algebra & character expansion; `[V-num]` no-small-parameter at `β=2.5`;
> `[NG]` character/EM duality as a gap-forcing tool; `[O]` the transfer = decompactification = `H1`.

---

## §0 · Labels, iron laws, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric this round.
**IL-i** non-positivity (fail on free massless). **IL-ii** nonabelian. **IL-iii** no small parameter at `α_s~O(1)`.
**Acid**: `Φ(free massless)≤0` *and* `≤0` on the gapped-atomless `:φ²:`.

---

## §1 · The three duality coordinates and their rigorous status

1. **Character / Fourier-on-the-group expansion `[T]`** — the only exact dual coordinate available.
   `e^{βRe tr U} = Σ_r d_r a_r(β) χ_r(U)`; the dual (random-surface / strong-coupling) expansion weights a
   surface `S` by `∏_p a_{r_p}(β)`. Controlled ⟺ `a_{r≠0}/a_0 ≪ 1`. This is the standard strong-coupling
   expansion; the area law at `β<1/48` is **proven** (Osterwalder–Seiler 1978; Shen–Zhu–Zhu, *CMP* 400 (2023)).

2. **'t Hooft loop algebra `[T]`** — the rigorous content of "electric–magnetic duality" in pure `SU(N)`:
   the electric (Wilson) `W(C)` and magnetic ('t Hooft) `B(C′)` loop operators obey the **center algebra**
   `W(C)B(C′) = B(C′)W(C) e^{2πi\,link(C,C′)/N}` ('t Hooft 1978). This is a theorem and it is a *base*
   fact: it proves the order parameter **must be non-local** (a loop/surface) and that the relevant global
   structure is the center `Z_N`. But it is an **algebra, not a small parameter** — it fixes the object's
   type, it does not evaluate it.

3. **Deformed small-circle semiclassics `[C-partial]`** — on `ℝ³×S¹` with a center-symmetric deformation
   the theory is weakly coupled (a genuine small parameter; monopole-instantons/bions; calculable gap,
   Unsal–Yaffe). This is a real duality *with* a small parameter — **but on a deformed theory at small
   circle**. Removing the deformation and growing the circle (decompactification) is the transport wall.

> **What is NOT available.** A true exact electric–magnetic duality for pure 4D YM (the kind that
> maps strong coupling to weak) is proven only with supersymmetry (`N=4` Montonen–Olive; `N=2`
> Seiberg–Witten). Large-`N` (t Hooft) and the exact 2D solution are *different limits / different
> models*, not dual coordinates of the 4D finite-`N` continuum problem. So the honest duality inventory is
> exactly {character expansion, 't Hooft loop algebra, deformed small circle}.

---

## §2 · The measurement: no dual small parameter at the crossover

`w7/base_duality_witness.py` evaluates the character-coefficient ratios (`I_r(2β)/I_0(2β)`, the SU(2) /
abelian proxy; the robust convention-independent fact is `a_r/a_0 → 1` for every compact group as
`β→∞`):

| β | a₁/a₀ | a₂/a₀ | a₃/a₀ | regime |
|---|---|---|---|---|
| 0.020 (`<1/48`) | 0.020 | 0.0002 | 0.0000 | **strong coupling — gap PROVEN (OS/SZZ)** |
| 0.125 | 0.124 | 0.008 | 0.0003 | dual marginal |
| 0.500 | 0.446 | 0.107 | 0.018 | dual marginal |
| 1.000 | 0.698 | 0.302 | 0.093 | dual marginal |
| 2.000 | 0.864 | 0.568 | 0.295 | **crossover: NO small parameter** |
| **2.500** | **0.893** | **0.643** | **0.379** | **crossover: NO small parameter** |
| 5.000 | 0.949 | 0.810 | 0.625 | weak: dual expansion degenerate |

**Reading.** `β<1/48`: `a_r/a₀` is small, the dual surface expansion converges — but this is exactly the
regime where the gap is *already a theorem*, so the duality is a **backward certificate**, not a key.
`β=2.5`: `a₁/a₀≈0.89`, `a₂/a₀≈0.64`, `a₃/a₀≈0.38` — surface weights are `O(1)`, the dual series is
**uncontrolled exactly where the problem lives**.

---

## §3 · Vortex dilution (the order parameter's own control parameter)

If the order parameter is a center vortex / 't Hooft loop, a controlled *dilute* picture needs
core size `≪` lattice spacing. Core `~1/√β` (lattice units):

| β | core `~1/√β` | dilute? |
|---|---|---|
| 0.5 | 1.41 | no |
| 1.0 | 1.00 | no |
| 2.0 | 0.71 | no |
| **2.5** | **0.63** | **no — vortices overlap** |
| 3.0 | 0.58 | no |

At `β=2.5` the vortex core (`≈0.63 a`) is **not thin**: there is no dilute gas, no controlled
semiclassical/vortex dual picture. So even the correct *type* of order parameter has no small control
parameter at the crossover.

---

## §4 · Where the difficulty transfers (honest answer to the directive)

- **Character duality:** transfers a small-`β` expansion to *another small-`β` statement*. At `β=2.5` it
  degenerates (`a_r/a_0` `O(1)`) — the difficulty **does not move**; it is the same weak-coupling side.
- **'t Hooft loop:** the area/perimeter law at `β=2.5` *is* the confinement statement, i.e. the gap. The
  duality is a **dictionary**, not a reduction: Wilson area law ↔ 't Hooft perimeter law is a rewriting
  of the same unknown.
- **Deformed small circle:** genuinely trades strong coupling for weak coupling (`Unsal–Yaffe`), so the
  gap *is* calculable there. But the object computed is a **deformed** theory; the transfer to pure 4D
  YM is decompactification → continuum → the transport wall. **The difficulty transfers to exactly
  the `a→0`, large-circle limit `= I2 = H1`.**

---

## §5 · Verdict and what it contributes to the base

- **Base contribution `[T]` (real):** the order parameter is **non-local** ('t Hooft loop / center vortex
  / Wilson loop), and the relevant global structure is the **center `Z_N`** — a rigorous anchor (`'t Hooft
  algebra`) that any future construction must respect. The dual coordinate is the **character
  expansion**. This is a genuine plinth: it says *what kind of object* the creature's base carries.
- **No key `[NG]`:** character / electric–magnetic duality is **not a gap-forcing tool**. At `β=2.5` it
  has no small parameter, and its exact version requires SUSY (absent here). It passes the acid test only
  vacuously (it is an exact reorganization; it does not *force* a floor — it returns the same
  weak/gapless branch at the crossover). So it is acid-**neutral**, not acid-**clean**.
- **Exact missing step `[O]`:** wherever a real dual small parameter exists (deformed small circle), the
  transfer to pure `ℝ⁴` is decompactification = `I2 = H1`. **Difficulty stays at the crossover; not
  transferred.**
- **No QED. Clay 0%.**

*Refs:* 't Hooft 1978 (loop algebra); Osterwalder–Seiler 1978; Shen–Zhu–Zhu, *CMP* 400 (2023);
Unsal–Yaffe (deformed small-circle semiclassics). Related: `PROFILE.md` (P1: non-local field
measure), `SLOT-SYNTHESIS.md` (all local data classes are floor-blind). Witness:
`w7/base_duality_witness.py` (runs clean, `py -3.12 -I`).

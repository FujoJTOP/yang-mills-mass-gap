# STAR-NO-IRFREE — pure YM's IR is NOT IR-free (the "free massless IR" branch of (★) is excluded)

> **Verdict (link of the (★) chain).** The "IR-free gapless" branch — the IR being a *free massless* sector —
> is **excluded for pure SU(N) YM**. The load-bearing, non-circular reason is **asymptotic freedom alone**:
> `β(g)<0` for all `g>0` ⇒ `g` runs **away from 0** into the IR ⇒ the IR coupling is `O(1)` at the confinement
> scale `Λ`, never `0`. Supporting (with a confinement caveat): no massless gluon/photon, no massless
> composite/Goldstone, no moduli. **[T]** given AF + the definition of `Λ`; the "no massless sector" sub-links are
> `[C]`. No QED.

---

## §0 · The three realizations of "IR-free", and why each fails

**1. Trivial IR fixed point `g*=0` — EXCLUDED (this is the rigorous link).**
Pure YM: `β(g) = −β₀g³ − β₁g⁵ − …`, `β₀ = 11N/(3·16π²) > 0`, `β₁ = 34N²/(3·(16π²)²) > 0` (no matter ⇒ no
Banks–Zaks). So `β(g)<0` for **all** `g>0` (witness `star_no_irfree.py`: SU(2),SU(3), zero sign-changes, no
zero on `g∈(0,10]`). Integrating `dg/dlnμ = β(g)` drives `g` **up** as `μ` falls: the UV `g→0` is the *only*
free point, and the IR flows **away** from it. Hence `g_IR = O(1)` at `Λ` (the scale where `g` becomes strong),
**not `0`**. This uses only `β₀>0` (a theorem) + the definition of `Λ` — **it does not assume the gap.**

**2. Free massless gluon/photon — EXCLUDED `[C]`.** SU(N) is simple: no unbroken `U(1)` factor, hence no
Coulomb-phase photon; the gluon is not a free IR degree of freedom. (Relies on `θ=0` / no Coulomb phase.)

**3. Free massless composite / Goldstone — EXCLUDED `[C]`.** A massless composite needs either a spontaneously
broken **continuous** global symmetry (Goldstone theorem) or an exact protecting symmetry. Pure YM has **no
continuous global symmetry** (internal symmetries are gauge or the discrete center `Z_N`), and no moduli
(isolated vacuum). So no protected massless composite.

## §1 · The two-gluon-cut subtlety (the directive's warning)

The `0⁺⁺` channel is a **two-gluon cut**. But a *cut* is a continuum in the spectral function of **confined
constituents** — it is **not** a free massless gluon propagating in the IR. The existence of the cut does not
make the IR free; it is a feature of the (confined) spectral density. So "two-gluon cut" ≠ "IR-free".

## §2 · Acid test

A free massless field has `β(g)=0`: its IR **is** free. So the statement "IR not free" is **FALSE on the free
field** — exactly the discriminating behavior. The argument is **specific to the interacting AF gauge theory**
(uses `β<0`, `β₀>0`), not vacuous.

## §3 · Honest scope

- **[T]** (AF + `Λ`): the IR coupling is `O(1) ≠ 0` ⇒ **not a free** (weakly-coupled) theory. This link is clean.
- **[C]**: "no massless gluon/composite" — supporting, partly confinement-flavored; a strong-coupling IR fixed
  point (`β=0` at `g*~O(1)`, if it existed) would be an **interacting CFT**, still **not free** — so IR-free is
  excluded regardless of whether such a zero exists.
- **What this link does NOT do:** it excludes *free*, not *gapless-strongly-coupled* — the latter is the CFT
  branch left to the rest of (★). It is the correct, non-circular exclusion of the IR-free branch.

*Witness:* `w7/star_no_irfree.py` (runs clean). **Clay 0%; no QED.**

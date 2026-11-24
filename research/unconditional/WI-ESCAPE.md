# WI-ESCAPE — weight-invariance is **OUT of the OPE-RHO *integral* class** (the theorem does not forbid it) but **IN the A3 class** (a re-encoding of the Stieltjes pole), so it does NOT escape the wall: what it stabilizes at `L→∞` is **isolation `I2` = `H1`**

> **VERDICT (no QED).** WI (`𝒜=π/Im(τ_zero)` invariant under operator-reweighting) **is not
> forbidden by the OPE-RHO theorem**, because it is **not a finite measurable integral datum** — it is
> the **location of the nearest zero of an exponential sum**, a root of a transcendental equation, which
> is *not* weak-continuously sensitive to the measure (so it lies outside the lsc-integral class the
> theorem covers). In that precise sense WI is a **genuinely new type of handle**. **But it does not
> escape the wall.** WI reads the **analytic singularity type** of the Laplace transform — a *simple
> pole* (isolated atom) vs a *branch cut* (continuum) — which is **exactly A3** (`f₁² = −Res Π`), the
> **read-not-force** object. At `L→∞` the stable invariant is "the first level is an isolated pole,"
> i.e. **the isolation `I2`** — which is part of the gap (`H1`). **Out of one theorem's class; into the
> same wall, in the A3 coordinate. Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric (witness).

---

## §1 · What WI is, structurally

At finite `L`, `H` is discrete and the `0⁺⁺` two-point function is the exponential sum
`Z(τ) = w₀ + Σ_i f_i² e^{−E_i τ}` (`w₀=|⟨tr F²⟩|²` the VEV). Its zeros in complex `τ` encode the
**singularity type of the Laplace transform** `Π(z)=∫ρ(s)/(z+s)ds`:

- **isolated atom** `f₁²δ(s−m₁²)` ⇒ a **simple pole** at `s=−m₁²` ⇒ `Z(τ)=w₀+f₁²e^{−m₁²τ}` has zeros at
  `τ = [ln(f₁²/w₀) − iπ(2k+1)]/m₁²`. The **imaginary part `Im τ = −π/m₁²` is independent of `f₁²`**; the
  weight only shifts `Re τ = ln(f₁²/w₀)/m₁²`. Hence `𝒜=π/|Im τ| = m₁²` **exactly, weight-invariant**.
- **continuum** (branch cut) ⇒ `Z` is a sum/integral over a band; the nearest zero **moves** with the
  weight distribution and with the density.

Witness (`w7/wi_escape_witness.py`) confirms exactly this: atom `𝒜=1.000000` for `f₁²∈{0.05,1,20}`
(only `Re τ` moves); two close levels `(1,1.3)` give `𝒜=1.14, 1.03, 1.23, 1.06` as the weight split
varies; resolving a band `[1,1.3]` at `N=2,3,5,8` gives `𝒜=1.137,1.136,1.134,1.133` (moving with
density). **Weight-invariance ⟺ simple pole ⟺ atom; failure ⟺ cut/continuum.**

---

## §2 · The class question — WI is OUT of the OPE-RHO integral class `[T]`

OPE-RHO's obstruction is: *finite measurable **integral** data `∫ρφ` are weak-continuous (lsc) and cannot
pin a support property; a bump of weight `CΔ²→0` near the floor matches all of them with `Δ→0`.*

`𝒜` is **not** such an integral. It is defined through the **nearest zero of `Z`**, i.e. the root of a
transcendental equation in `τ`. The map `weights ↦ (zero location)` is **not weak-continuous**: under a
measure perturbation that changes the pole-to-cut structure, the nearest zero can appear, disappear, or
jump to a different `Im`. A functional built from a **zero location** is therefore **not** in the
lsc-integral class, and the OPE-RHO theorem's hypothesis does **not** apply to it. **Conclusion: OPE-RHO
does not forbid WI.** In this precise sense WI is a new type of datum — it exploits **analytic
singularity structure (pole vs cut)**, which no integral can see (integrals are blind to a thin atom
because it is `L¹`-small; the *pole location* is not `L¹`-controlled).

---

## §3 · But WI is IN the A3 class `[T]`

**Reading "pole vs cut" *is* reading the Stieltjes pole.** "There is a simple pole of `Π` at `s=−m₁²` with
residue `f₁²`" is verbatim A3 (`f₁² = −Res Π`), the corpus's **read-not-force** object: it is an *exact*
reading of the atom, but **evaluating it requires knowing the spectrum** (locating the pole = solving the
`0⁺⁺` eigenvalue problem). WI changes the *representation* of A3 (a zero-location instead of a residue) but
not its *content*. And the invariant WI stabilizes at `L→∞` is precisely **"the first level is an isolated
pole separated from the rest"** = the **isolation `I2`**, which with the atom weight is the gap (`H*`).
Hence at `L→∞`, **WI ≡ A3 ≡ `I2` ≡ `H1`** — the same wall in a new coordinate.

---

## §4 · Finite-`L` reality: a real diagnostic, not an atom certifier `[V-num]/[O]`

At finite `L`, **everything is a finite exponential sum**, so *every* spectral feature is a pole of the
finite-`L` `Z`. WI then distinguishes **"first level isolated from a packed ladder"** from **"first level
part of a dense ladder"** — a genuine, testable property, and **the raw `𝒜` alone is a scale detector, not
acid-clean** (it fires on `:φ²:` and on the free massless box, `CONSTRUCT-IR-ANALYTIC-NUM`). The
weight-invariance refines this to a real finite-`L` **pole-vs-ladder diagnostic**. But it does **not**
certify that the lowest level *survives as an isolated atom in the `L→∞` limit* — at each finite `L` the
"atom" could be a discretized continuum edge (it is, for free field / `:φ²:`). The witness shows the
continuum's `𝒜` keeps moving as `N` grows, so the finite-`L` "isolation" is only as sharp as the level
spacing — i.e. it is probing **existing isolation**, not manufacturing one.

**Exact failing step `[O]`:** proving that the weight-invariant `𝒜` is **stable uniformly in `(L,a)` at
`β≈2.5`** (with reweighting that preserves the coupling content) — equivalently that an isolated `0⁺⁺`
pole survives the continuum limit. That is `I2/H1`.

---

## §5 · Bottom line

- **Is WI in or out of the OPE-RHO class?** **OUT** — it is not a finite measurable integral datum, so the
  theorem does not apply; WI genuinely exploits analytic singularity structure that integrals cannot see.
- **Does it escape the wall?** **NO** — it is **A3 re-encoded** (Stieltjes pole = residue), the
  read-not-force object, and the invariant it stabilizes (`L→∞`) is the **isolation `I2` = `H1`**.
- **Net:** WI is a **genuine new finite-`L` diagnostic / representation** — a pole-vs-cut detector outside
  the OPE-RHO integral class — but it **pays the A3/H1 price**, so it re-locates the wall rather than
  crossing it. **A real new tool; not a gap certifier. Clay 0%; no QED.**

*Witness:* `w7/wi_escape_witness.py` (numpy). *Refs:* `PNT-OPE-RHO.md`, `PNT-IR-ANALYTIC.md`,
`CONSTRUCT-IR-ANALYTIC-NUM.md`, `PROOF-F1SQ.md` (A3), `SIGNED-FOURIER-*`, `PROFILE.md`.

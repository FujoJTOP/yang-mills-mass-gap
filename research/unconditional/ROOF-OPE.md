# ROOF-OPE — the OPE / SVZ sum-rule route: the anomalous Ward identity transports the ambiguity, dispersion + positivity cannot remove it, and the only renormalon-free object is the full spectral function (whose floor is `H1`)

> **VERDICT (no QED).** Tested the standard SVZ sum-rule setup as a way to pin the gluon condensate
> `renormalon-free` at the crossover, using OPE + anomalous (dilatation) Ward identity + dispersion
> relation + positivity. **It does not work.** The Wilson coefficient of the dimension-4 operator has the
> `u=2` IR renormalon; the upper/lower Borel prescriptions differ by `2πr·Λ⁴/Q⁴` — a **fixed multiple of
> the condensate's own term** at every scale. The anomalous Ward identity ties the condensate to the
> vacuum-energy density, which carries the **same** `Λ⁴` renormalon: it **transports** the ambiguity,
> it does not remove it. Dispersion + positivity constrain the **physical spectral function** ρ(s),
> which is renormalon-free — but the OPE *split* into "perturbative coefficient × condensate" is
> scheme-dependent, so every positivity **bound** on the condensate inherits the `±πrΛ⁴` ambiguity.
> The only renormalon-free object is ρ(s) itself; bounding the `0⁺⁺` atom/floor from ρ is exactly `H1`
> (finite-moment data are floor-blind). **`[NG]` for the OPE-pins-condensate class; `[O]` exact step =
> the `u=2` subtraction. Clay 0%.**

---

## §0 · Labels

**[T]** 真定理/standard textbook · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-here]** new this round.
**Acid test** (carried): a candidate must fail/return `0` on the free massless field *and discriminatingly* on the gapped-atomless `:φ²:`.
Convention note: "u = 2" is stated in the standard convention where the leading **dimension-4 IR renormalon** sits at Borel `u = 2`; the physical statement (the ambiguity is `O(Λ⁴/Q⁴)`) is convention-independent.

---

## §1 · The SVZ sum-rule setup

Two-point function of the `0⁺⁺` gluonic current at Euclidean `Q² = −k²`:

```
Π(Q²) = C₀(Q²) + C_G(Q²)·⟨G²⟩ + C₄(Q²)·⟨O₄⟩ + …        (OPE, SVZ 1979)
```

- **Current choice.** `J = tr G²` is *not* RG-invariant (mixes with the identity under renormalization);
  the RG-invariant `0⁺⁺` operator is the **trace anomaly** `θ = (β(g)/2g)·tr G²`. The `⟨G²⟩` condensate
  therefore means the vacuum matrix element of the trace-anomaly operator — the **vacuum-energy density**.
- **Dispersion.** `Π(Q²) = (1/π)∫ ds ρ(s)/(s+Q²)`, with `ρ(s) = Im Π(s)` the physical spectral density.
- **Positivity.** `ρ(s) ≥ 0` in the `0⁺⁺` channel (reflection positivity / Källén–Lehmann). `[T]`
- **Sum rule.** Borel-transforming in `Q²` (or a finite-energy sum rule over a duality window) trades the
  unknown ρ for the OPE series. The question is whether the **leading power correction** `⟨G²⟩` is thereby
  fixed *renormalon-free*.

This is the standard SVZ procedure, and the standard SVZ number is `⟨(α_s/π)G²⟩ ≈ 0.012 GeV⁴` (order of
magnitude; scheme- and truncation-dependent — see §3).

---

## §2 · The `u=2` IR renormalon: the ambiguity is the condensate's own size

The perturbative coefficient `C_G(Q²)` (and `C₀`) are computed as **asymptotic** series in `α_s`. Their
Borel transforms have an **infrared pole at `u = 2`**, the dimension-4 renormalon whose residue is exactly
the gluon-condensate ambiguity. `[T]` (standard renormalon calculus)

- The Borel sum along the real axis is ill-defined at the pole; upper vs lower contour prescriptions differ by
  `2πi·r·e^{−2t}`, giving a **real ambiguity `±πr·e^{−2t}`** (`r` = residue, `t` the Borel variable). `[T]`
- Standard one-loop correspondence `e^{−2t} ↔ (Λ²/Q²)²` ⇒ **ambiguity `= πr·Λ⁴/Q⁴`** — the **same
  `Λ⁴/Q⁴` power** as the condensate term `A·Λ⁴/Q⁴`. `[T]`

**Witness** (`w7/roof_ope_witness.py`, Part A): the ratio `ambiguity / (A·Λ⁴/Q⁴) = πr/A` is **exactly
constant in `Q²`** (`3.1416` at every scale for `r=A=1`). The ambiguity never shrinks relative to the
condensate — so there is **no renormalon-free *value*** of the condensate extractable from the OPE.

This is not a truncation artifact: **the condensate *is* the transseries coefficient of the `e^{−2s}`
term, and the renormalon ambiguity *is* the size of that same term.** Resurgence *identifies* the
condensate with the ambiguity — a diagnosis, not a construction. (Matches `ROOF-CONDENSATE-FLOOR`.)

---

## §3 · The anomalous Ward identity transports the ambiguity

The dilatation Ward identity `∂_μ J_D^μ = θ^μ_μ = (β/2g)·tr G²` is **exact** and makes the trace-anomaly
operator RG-invariant. `[T]` But its vacuum matrix element is the **vacuum-energy density**, which carries
the **same `Λ⁴` ambiguity** (the cosmological-constant renormalon `u = 2` for `d = 4`). `[T]`

⇒ The Ward identity **relocates** the ambiguity from the correlator to the vacuum energy; it supplies
**no independent handle** on the magnitude. The identity is a *constraint among ambiguous quantities*, not
a renormalon-free determination. (This is the OPE-sector restatement of `ROOF-CONDENSATE-FLOOR`'s
"the identity transports the ambiguity, does not remove it.")

---

## §4 · Dispersion + positivity: renormalon-free ρ, but an ambiguous *split*

- The **physical** spectral function `ρ(s) ≥ 0` is **unambiguous** — it is a property of the exact
  correlator, independent of any Borel prescription. `[T]`
- But the OPE **splits** `Π = C₀ + C_G⟨G²⟩ + …`; the split between "perturbative coefficient" and
  "condensate" is **prescription-dependent** at the renormalon level. Two prescriptions `s = ±1`
  reproduce the **same** physical `Π` (hence the **same** ρ) while giving condensate coefficients
  differing by `2πr·Λ⁴` (**witness Part B**: `A(+) − A(−) = 6.2832 = 2πrΛ⁴` at every `Q²`, both
  consistent with the identical positive ρ).
- **Positivity cannot choose between them.** ρ ≥ 0 is the *same object* for both prescriptions
  (**witness Part C**). Any positivity **bound** on the condensate — e.g. from `Π(Q²) > 0` —
  shifts by `s·πrΛ⁴` between prescriptions and therefore inherits the ambiguity. Positivity gives
  **bounds**, not a renormalon-free **value**, and the bounds are themselves ambiguous.

---

## §5 · The only renormalon-free object, and why it is `H1`

The one object free of the OPE split is `ρ(s)` itself. But to turn ρ into a statement about the `0⁺⁺`
**mass gap** one needs to locate the support floor `inf supp ν` / the atom weight — and **finite-moment
data from ρ are floor-blind**. `[T]` (base round)

- `SLOT-BOOTSTRAP`: `Δ_min(K) = 0` for **every** moment truncation `K` — the lower feasible support floor
  is pinned at `0` regardless of how many moments are supplied.
- `BASE-FIELDMEASURE`: every Markov-semigroup / gradient-flow functional of μ is homogeneous in the coarse
  mass (or zero, or sign-indefinite) — it is `κ` relabelled, not a source.
- ⇒ bounding the `0⁺⁺` atom/floor from the physical ρ (or from any finite moment data) **is exactly `H1`**.

---

## §6 · Verdict and the exact failing step

- **OPE does NOT pin the condensate renormalon-free.** `[NG]` for the class "SVZ sum rule + Ward +
  dispersion + positivity ⇒ renormalon-free `⟨G²⟩`."
- **Three independent reasons, any one sufficient:**
  1. the `u=2` IR renormalon makes the perturbative/condensate split ambiguous at `O(Λ⁴)` = the
     condensate's own scale (**§2**);
  2. the anomalous Ward identity routes RG-invariance to the vacuum energy, which carries the same
     `Λ⁴` renormalon (**§3**);
  3. positivity constrains the physical ρ, which is unambiguous, but every **bound** on the condensate
     inherits the prescription ambiguity (**§4**).
- **Exact missing step `[O]`:** a renormalon-free (scheme-independent) characterization of the `0⁺⁺`
  anomaly — equivalently a rigorous control of the `u=2` subtraction uniform in `(L, a)`. That is the
  continuum/transport wall in the OPE coordinate.
- **Renormalon-free survivor:** ρ(s) (and positivity bounds on it); but the floor/atom extracted from ρ
  is `H1`. **`H1` in the OPE coordinate.**

**Clay 0%; no QED.**

*Refs:* `ROOF-CONDENSATE-FLOOR.md`, `BASE-FIELDMEASURE.md`, `SLOT-BOOTSTRAP.md`, `BASE-SYNTHESIS.md`,
`ROOF-SYNTHESIS.md`. Witness: `w7/roof_ope_witness.py`. Standard inputs (SVZ OPE; `u=2` IR renormalon;
dispersion; RP positivity of ρ) are `[T]` and labelled as standard, not re-derived here.

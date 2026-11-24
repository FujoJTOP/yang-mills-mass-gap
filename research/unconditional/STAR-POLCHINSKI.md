# STAR-POLCHINSKI — the scale-invariance ⇒ conformal-invariance link for target (★)

> **VERDICT.** The Polchinski link is **[T] and solid**: in d=4, a unitary QFT with a conserved symmetric
> (improvable) stress tensor that is **scale-invariant** is **conformally invariant** (Polchinski, *Nucl. Phys.*
> **B303** (1988) 226). Pure YM satisfies the **structural** hypotheses (symmetric conserved improvable
> `T^{μν}`; unitarity). **But (★) does NOT close** — it needs two further, separately-open steps:
> **(a)** the hypothetical gapless IR of pure YM is **scale-invariant** (not automatic), and
> **(b)** pure YM has **no conformal IR**. The sharp handle for (b) is the **trace anomaly**:
> `T^μ_μ = (β/2g) tr F² + (Weyl², Euler)`; IR conformality requires this operator to vanish in the IR, i.e.
> `β(g_IR) = 0`. Pure YM's β-function is AF with **no Banks–Zaks zero**, so **no perturbative** IR fixed point;
> ruling out a **strongly-coupled** CFT IR needs more than β. Both (a) and (b) are open = **`H1`-adjacent**.
> So: `[T]` for the Polchinski link; `[O]` for (a) and (b). (★) reduces to `H1`, but in a **sharpened
> trace-anomaly coordinate** the campaign had not made explicit.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The theorem (Polchinski 1988), precisely

**Theorem (Polchinski, d=4).** Let a four-dimensional QFT satisfy:

- **(H1) Unitarity** — the Hilbert space has a positive-definite inner product.
- **(H2) Symmetric, conserved, improvable stress tensor** — `∃ T^{μν}=T^{νμ}`, `∂_μT^{μν}=0`, and `T^{μν}` is
  Belinfante/improvement-improvable so that the **trace** `T^μ_μ` is a well-defined local scalar and the
  **dilation current** `D^μ = x_ν T^{μν}` exists.
- **(H3) Scale invariance** — `∂_μD^μ = T^μ_μ = 0` as an operator equation. (In general `T^μ_μ = ∂_μ V^μ` for a
  local **virial current** `V^μ`; scale invariance means `D^μ` is conserved up to at most a trivial total
  derivative.)
- **(H4) Discrete spectrum** — no continuous spectrum of scaling dimensions / a mild gap, so the OPE and the
  unitarity bounds apply.

Then the theory is **conformally invariant**: the special-conformal currents
`K^{μν} = 2x^νx_λT^{μλ} − x²T^{μν}` (+ improvement) are conserved and the full algebra `so(4,2)` is realized.

**The Polchinski subtlety (virial obstruction).** In general `T^μ_μ = ∂_μ V^μ` with `V^μ` a vector operator of
dimension `d−1 = 3`. A non-trivial, non-conserved `V^μ` would give *scale without conformal*. Polchinski shows
that in d=4 with unitarity and symmetric `T^{μν}`, the 2-point function of `V^μ` together with the unitarity
bounds forces `V^μ` to be **trivial** (a total derivative / conserved) **unless** the theory is non-unitary
(negative-norm states) or has an exotic virial dimension. Hence, for a **unitary** theory: **scale ⇒ conformal.**

**Proof sketch.** Scale invariance ⇒ `D^μ` conserved. The conformal currents satisfy
`∂_μK^{μν} = 2x^ν T^μ_μ + (total derivatives)`; so if `T^μ_μ=0` the `K^{μν}` are conserved and the
conformal algebra closes. The only gap is whether `T^μ_μ = ∂_μV^μ` can be improved to zero; unitarity + the
OPE of `V` with itself (via the `so(4,2)`-compatible unitarity bounds) forces `V^μ` trivial.

---

## §2 · Check against pure YM's hypothetical gapless IR

Pure YM (4D, `SU(N)`, θ=0):

| hypothesis | pure YM |
|---|---|
| (H1) unitarity | **✓** (OS reconstruction; reflection positivity) |
| (H2) symmetric conserved improvable `T^{μν}` | **✓** (gauge-invariant Belinfante-improved stress tensor) |
| (H4) discrete spectrum | ~ (a gapless IR has a continuum, but the local OPE / unitarity bounds still apply) |
| (H3) **scale invariance** | **← THE CRUX** |

**The sharp handle is the trace anomaly.** For YM,
`T^μ_μ = (β(g)/2g)\, tr F² + (Weyl², Euler)`.
This is **not zero** as an operator in the interacting theory — `tr F²` is a nonzero local operator and
`β(g)/2g ≠ 0` for `g>0`. So the IR can be scale-invariant only if the **anomaly operator vanishes in the IR**:
`β(g_IR) = 0` — an IR fixed point.

- **Pure YM's β-function is asymptotically free**: `β(g) = −(β₀g³ + β₁g⁵)/(16π²) + …` with
  `β₀ = 11N/3 > 0`, `β₁ = 34N²/3 > 0` ⇒ **`β(g) < 0` for all `g>0` at 2 loops: no Banks–Zaks zero.**
  (Witness below: `β<0` on all `g>0` for `N=2,3,4`.)
- Hence **no perturbative IR fixed point**. A **strongly-coupled** CFT IR (not described by the perturbative β)
  is **not excluded by β alone**, but has no known realization and is not expected for pure YM.

⇒ **(★) decomposes into two open steps:**
1. **(a)** a hypothetical gapless AF IR is **scale-invariant** (so Polchinski applies). Not a theorem — a gapless
   IR could in principle be non-scale-invariant ("scale without conformal", or a dimensionful coupling).
2. **(b)** pure YM has **no conformal IR** — the trace-anomaly operator `(β/2g)tr F²` must vanish in the IR
   (`β(g_IR)=0`). The no-BZ β rules out **perturbative** fixed points; excluding a strongly-coupled CFT IR
   requires more (the `a`-theorem, the 't Hooft anomaly, or a nonperturbative β statement).

---

## §3 · Verdict

- **Polchinski link `[T]`:** scale-inv + unitarity + symmetric `T` ⇒ conformal, in d=4 — **solid**.
- **Pure YM structural hypotheses (H1,H2,H4):** satisfied.
- **(★) does NOT close:** it needs (a) scale invariance of the gapless IR and (b) absence of a conformal IR —
  both open. The trace anomaly `T^μ_μ = (β/2g)tr F²` is the sharpest statement of the obstruction: **IR
  conformality ⟺ the anomaly vanishes in the IR ⟺ an IR fixed point**, which pure YM lacks perturbatively.
- **Honest:** `[T]` Polchinski link; `[O]` steps (a) and (b). (★) is a genuine, well-posed reformulation, and
  the trace anomaly is a **new, sharper axis** on the same gap — but the β-function alone does **not** deliver
  the "no conformal IR" needed to close it. **Clay 0%; no QED.**

*References:* J. Polchinski, "Scale and conformal invariance in quantum field theory," *Nucl. Phys.* **B303**
(1988) 226–236. (Companion literature: A. Zamolodchikov, *JETP Lett.* **43** (1986) 730 — `c`-theorem; the
"scale without conformal" reviews of Nakayama; Jackiw–Pi trace-anomaly improvements.)

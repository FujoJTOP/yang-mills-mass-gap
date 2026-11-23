# SCHUR-CONVEX — can a protected total be bounded ABOVE by the edge atom in 4D YM?

> 2026-11-23 · Direct attempt to **construct** the one object EDGE-f1sq named: a coercivity `∫φ dν ≤ C(φ)·ν({m₁²})` — a *protected total bounded above by the edge atom*, the reverse of every RP/Cauchy–Schwarz inequality. **Verdict: not constructible from protected data.** A master continuity obstruction kills all four routes at once; per-route there is a further sign or regime wall. What the construction would *be*, if posited, is verbatim `H*(r₀)`. No fabricated proof.

**Normal forms.** `ν≥0` KL measure of a gauge-invariant interpolator, `supp ν⊂[m₁²,∞)`, `f₁²=ν({m₁²})`. A protected total is `T_φ=∫φ dν` (moments `M_k`, `χ_t=M_{−1}`, norm `M₀=1`). With `M₀` fixed, the coercivity is equivalent to

  **(RD)** `f₁² ≥ M₀/C`, i.e. edge fraction `ϑ:=f₁²/M₀ ≥ 1/C>0` (resonance/edge dominance).

The barrier `R:=T_φ − C·ν({m₁²})` *is* Schur-convex: `∫φ dν` is Schur-convex (φ convex, Jensen under mean-preserving spreads) and `−ν({m₁²})` is Schur-convex (the atom weight is Schur-concave). So `{R≤0}` is the majorization-**concentration** set `{ϑ≥…}`. The object is correctly *shaped*. The question is provability.

**Master obstruction (free-field no-go).** Take `O=tr F²`, the anomaly/`0⁺⁺` channel itself. Every protected datum — anomaly coefficient, `M₀` normalisation, KL positivity, Ward/BRST — is a consequence of structure that is *exact or continuous in `g` and holds at `g=0`*. The free theory (`g=0`) is admissible: `ν` is a two-photon **cut** on `[0,∞)`, purely absolutely continuous ⇒ **`f₁²=0`**, while `M₀=⟨O†O⟩=1>0`. Then (RD) reads `1 ≤ 0`. Hence **no inequality provable from protected, `g`-continuous data implies the coercivity**: it must *fail* exactly at `g=0` and be restored non-perturbatively as `g→g_phys`. That restoration *is* the gap. (Also: for `tr F²` the bare `M₂=∫s²dν` is UV-divergent; only the anomaly *log-coefficient* is scheme-independent — a UV number, not a finite IR total.)

**The four candidate routes — exact death.**
1. *Anomaly / `χ_t` total, proven Schur-convexity.* Schur-convexity is genuine (`s^k`, `1/s` convex on `[m₁²,∞)`; this is precisely why `χ_t≥f₁²m₁² ⇒ f₁²≤χ_t/m₁²` — atom ≤ total, the **wrong** direction). `χ_t=M_{−1}` lives in the `0⁻⁺` (`q=∂·K`) channel whose atom is **unprotected/scheme-dependent** (GAP-IN11); the `0⁺⁺` total is divergent. Schur-convex with a *usable* total: absent.
2. *Density-of-states / counting.* A δ atom has zero spectral *density*; every Weyl-type cumulative/counting bound is Schur-concave and admits `atom=ε`, continuum filling the rest. Counting is blind to atoms. Reinforced by GAP-G2-mass(c): the edge mode is a delocalized `k→0` collective mode, so every local-density functional carries `λ_min→0` — no density bound sees the edge atom.
3. *Lattice protected sum.* A genuine coercivity exists only for `β<β₀` (OS'78; SZZ `β<1/48`), where the gap is *already* a theorem; it cannot cross `β_c≈2.47` (15–120× beyond, T5) without the uniform-in-`a` RG control that is IN-9. Real sign, wrong regime.
4. *Nonlinear Ward/STI at `k=0`.* Ward/STI are exact all-coupling total/contact relations; the transverse residue has **no anomalous slot** (ghost#1, dim-5 empty — W3-marginal-BPHZ; GAP-G2-mass). Exactness forces them to hold at `g=0` where `f₁²=0` (master obstruction). Products add totals, not a residue slot.

**Identification.** Granted, the coercivity yields `f₁²≥r₀>0` — **verbatim `H*(r₀)`** of RESIDUE-CONDITIONAL, the campaign's single unprotected input, strictly stronger than a bare gap. The Schur-convex coercivity is not a new lever; it is `H*` in measure language.

**Sharpest line.** `f₁²` is a Schur-*concave* functional; the coercivity is a *lower* bound on it. All provable data give upper bounds on Schur-concave (lower on Schur-convex) functionals — the opposite cone — and the free theory sits legitimately on its boundary `ϑ=0`. A protected, `g`-continuous total that is Schur-convex-coercive in 4D YM does not exist; constructing one is proving edge dominance, i.e. the gap.

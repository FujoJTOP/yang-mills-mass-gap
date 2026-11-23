# INEQ-R6 — 4-point spectral positivity + 4-point Ward vs. the residue `f₁²`

> 2026-11-23 · The 4-point *was* the one place to hope the pole Ward stops being homogeneous: an STI on `⟨O O T T⟩` carries an anomalous contact that does not scale with `f₁`, so one might solve `f₁²·V=a` non-homogeneously. **Verdict: the homogeneity is not broken in a gauge-invariant connected correlator.** The only protected contact reachable in `⟨O O T T⟩` is the **T-norm** (central charge `c`), and it enters **as `f₁²·c`** — multiplying the atom exactly as the 3-point did. The **a-anomaly `E₄` requires ≥3 T-legs** and never meets an O-contraction without an `f₁²` factor. Re-collapses to aggregate. Constructed nothing; no fabrication.

**Setup.** Double spectral rep `G⁽⁴⁾=∫∫ ρ₂(s,t)/(…)`, `ρ₂≥0` by unitarity; `O=¼trF²`, `f₁²=ν({m₁²})`.

**Construction attempt (best faith).** In the 3-point every residue is `f₁·(vertex)·f₁`, and `g^{μν}Γ_μν=(β/2g)⟨OOO⟩` multiplies both sides by `f₁²` → cancels (INEQ-T4). Suppose a 4-point STI `∂·T[O O T]=…` adds an anomalous contact `a·C(p)` **not** scaling with `f₁`, giving `f₁²·V(p) + a·C(p) = (pole parts)`. A once-subtracted dispersion then gives `f₁²·V = a − ∫∫_{heavier} ρ₂`; a lower bound would follow if `V` and the heavier double-spectrum were bounded **above**.

**Why it dies — four exact points.**
1. **No f₁-free contact survives connectivity.** The `E₄`/`a` contact is a *gravitational* (TT) contact; inside `⟨O O T T⟩` it is `⟨OO⟩·(TT-contact)=f₁²·c` — homogeneous. The scalar a-anomaly proper lives in `⟨TTT⟩`/`⟨TTTT⟩` (needing a third T); the f₁-free pieces it can leave with one O are *disconnected* (`⟨O⟩×a`) or pure-T. None is the connected amplitude carrying `f₁²·V`, so `V` is never isolated — `f₁²` keeps its factor.
2. **The protected contact is a UV moment.** When it does appear it is the `x→0`/high-`s` OPE coefficient — an *integral* of `ρ₂` (aggregate), fixing a sum rule, never the atom (INEQ-T4 leg 2; INEQ-T5: UV vs edge).
3. **`V` is unbounded.** Granting `f₁²·V=a`: at `g=0`, `a=a_free≠0` (type-A/B anomalies are one-loop-exact) while `f₁²=0` ⇒ `V=∞`. The vertex diverges precisely where the atom vanishes; the free-field no-go is localized here.
4. **Positivity is the wrong cone.** `ρ₂≥0`; the Gram 2×2 on `|O⟩,|OO⟩` gives `|⟨O†OO⟩|²≤⟨O†O⟩⟨(OO)†OO⟩`, i.e. `f₁²g²≤M₀·(4pt)` ⇒ **upper**. Positivity never lower-bounds a diagonal entry beyond `≥0`.

**On `⟨TTTT⟩` vs `⟨O O T T⟩`.** `⟨TTTT⟩` supplies only the pure-gravity `a` (no O-legs); `⟨OOTT⟩` supplies `f₁²·c`; their dispersive ratio is `f₁²·(c/a)·(kinematics)`, and the kinematics is the unprotected vertex again — the same wall.

**Sharpest line.** The 4-point's extra spectral power is real but Schur-convex in the *aggregate* (double) measure, and every protected contact it exposes is either homogeneous (`f₁²·c`) or disconnected. A lower bound still needs the heavier double-spectrum **and** the off-shell `⟨G|…⟩` vertex bounded above — two-variable resonance/duality saturation, strictly harder than the 2-point assumption and equally an assumption. **Flip condition:** a proven upper bound on the off-shell vertex plus two-variable duality ⇒ `f₁²≥a/V_max>0`; none known. The edge theorem (`H*(r₀)`) is unchanged.

*Refs.* INEQ-T4 (homogeneous pole Ward); INEQ-T5 (UV vs edge; positivity = upper); SCHUR-CONVEX (wrong cone); B1 (free-field no-go); GLUEBALL-BOOTSTRAP (⟨TTTT⟩/a-anomaly ill-posed, protected sector only); GAP-IN11 (moment lemma); RESIDUE-CONDITIONAL (`H*`).

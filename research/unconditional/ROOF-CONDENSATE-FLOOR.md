# ROOF-CONDENSATE-FLOOR — "floor the gluon condensate uniformly in `a`" is **not a well-posed scheme-independent demand**: the `u=2` IR renormalon is the exact obstruction

> **VERDICT (no QED).** The base round converged on the gluon condensate `⟨tr F²⟩` as the unique `0⁺⁺` anomaly `= H1`. Attempting to *floor* it uniformly as `a→0` (the roof) runs into a sharp, structural obstruction: **`⟨tr F²⟩` is manifestly bounded below only at finite `a` (where the bare plaquette defect `≥0`), but the finite, renormalized "condensate" is defined only up to an IR-renormalon ambiguity at `u=2` whose size is `~Λ⁴` — the same order as the condensate itself.** The perturbative/nonperturbative split is therefore scheme-dependent and ambiguous at `O(1)` relative precision: there is no scheme-independent, `a`-uniform lower bound on "the condensate number." **The roof cannot be laid by flooring `⟨tr F²⟩` as a number.** A floor, if it exists, must come from a *different characterization of the same anomaly* (the trace anomaly as an exact identity, or the 3D-Z₂ magnetic-defect-condensation mechanism in the `0`-form channel) — not from the OPE-renormalized condensate. **Clay 0%; no QED.**

---

## §0 · Labels, task

Labels: **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** schematic numeric witness · **[R]** inherited.
Task: separate the perturbative (divergent) and nonperturbative (floored) parts of `⟨tr F²⟩`; find a scheme/definition where the nonperturbative part is provably floored uniformly in `a`; analyze the `u=2` IR renormalon and whether Borel/resurgence controls it.

---

## §1 · The separation problem (why "the condensate" is not a number)

The gluon condensate enters as the dimension-4 operator in the OPE, `⟨(α_s/π)G²⟩_μ`, with a known anomalous dimension. Three separate facts make the split non-canonical:

1. **The bare lattice object is pointwise nonnegative but divergent.** The plaquette-action density `1 − (1/N)Re tr U_□ ≥ 0` on every configuration `[T]`; but `⟨tr F²⟩ ∝ a^{-4}(1 − ⟨Re tr U_□⟩/N)` diverges as `a→0` (the perturbative series `= Σ c_n α_s^n` grows factorially).
2. **The subtraction is scheme-dependent.** Isolating a finite "condensate" requires subtracting the perturbative series to all orders — which is itself the renormalon-ambiguous object (§2).
3. **The trace anomaly route does not escape it.** `⟨θ^μ_μ⟩ = (β(g)/2g)⟨tr G²⟩` ties the condensate to the vacuum-energy density, but the vacuum-energy density carries the **same** `u=2` (Λ⁴) renormalon ambiguity (the cosmological-constant / bag-constant ambiguity), so the anomaly identity *transports* the ambiguity rather than removing it.

---

## §2 · The `u=2` IR renormalon is exactly the condensate's own scale

The Borel transform of the perturbative series for the condensate has an **IR renormalon at `u = d/2 = 2`** (dimension-4 operator; leading IR renormalon of the gluon condensate). Its consequences `[T]/[R]`:

- **Ambiguity `~ (Λ²/Q²)^2 = Λ⁴ / Q⁴`** in the dimensionless observable — i.e. an absolute ambiguity `~Λ⁴`, the condensate scale.
- Equivalently, **the position `u=2` is exactly where "perturbative" ends and "nonperturbative" begins**; the condensate *is* the coefficient of the transseries term `~e^{-2s}` (`s=ln(Q²/Λ²)`), and the renormalon ambiguity is the size of that same term.

**Schematic witness** (`w7/roof_condensate_floor_witness.py`, illustrative — *not* a YM calculation): for a Borel pole at `u = d/2`, the discontinuity across the cut scales as `e^{-(d/2)s}` with a **scale-independent** coefficient (`≈4π` full / `2π` half-width). Relative to a natural condensate term `c·e^{-2s}`, the ambiguity is `2π/c` **at every scale** — a fixed fraction, never shrinking. For `c < 2π` the ambiguity **exceeds** the condensate, so the interval `[PV − |Δ|, PV + |Δ|]` always straddles zero: *the "nonperturbative condensate" can be dialed to zero inside the ambiguity.* The `d=6` control (`u=3`) is parametrically smaller by `e^{s}`, confirming `d=4` is the marginal, leading case.

**This is the obstruction, not a technicality:** a quantity whose definition carries an ambiguity *of its own size* cannot be floored in a scheme-independent way.

---

## §3 · Why no rigorous uniform floor survives

- **Finite `a`:** a "floor" at finite `a` (`⟨plaquette defect⟩ ≥ 0`) is real but **trivial and divergent** — it is the bare `1/a⁴` quantity, not the physical condensate.
- **Continuum `a→0`:** the subtracted value approaches the **renormalon-ambiguous** continuum number; any finite-`a` "floor" is a lattice artifact, not a bound on the physical object.
- **Trace anomaly / vacuum energy:** transports to the same Λ⁴ renormalon (cosmological-constant ambiguity); the rigorous positivity `H ≥ 0` is *after* subtracting the vacuum energy to zero and does not floor the anomaly's magnitude.
- **Borel/resurgence:** *does not rescue* — resurgence **identifies** the condensate with the renormalon ambiguity (the condensate is the transseries coefficient). That is a *diagnosis* of why there is no floor, not a construction of one.

---

## §4 · The exact failing step, and what would actually lay the roof

**Exact failing step `[O]`:** the **`u=2` IR-renormalon subtraction** — the boundary between perturbative and nonperturbative parts of `⟨tr F²⟩`. Making it rigorous (a scheme-fixed Borel resummation with a *controlled* `u=2` ambiguity) is an open problem in renormalon theory, and controlling it uniformly in `a` is the same continuum/transport wall in disguise.

**Hence the roof must be laid via a different characterization of the same anomaly**, not the condensate number:

1. **Defect condensation in the `0`-form channel** (the base round's named route): prove the *3D-Z₂ magnetic-defect-condensation mechanism* transferred to 4D `0⁺⁺` at `α_s~O(1)` — a statement about the measure/defects, not the OPE-renormalized number, hence immune to the `u=2` ambiguity.
2. **Trace anomaly as an exact identity** with a *renormalization-group-invariant*, scheme-independent left-hand side (e.g. the vacuum energy density defined by a finite-volume gap, not the OPE condensate) — a floor on the invariant, if one exists, must be proven on that side.

Both remain **open**; both are the same `H1` in different coordinates.

---

## §5 · Acid test and iron laws

- **Acid** `[T]`: on the free massless field the perturbative series has no condensate (condensate `= 0`, only the divergent perturbative part); the witness's `a`-dependence (`∝ initial scale`) returns the free behavior. Acid passes — the obstruction is *specific to the interacting, gapped-measure* side.
- **IL-i** (non-positivity): the obstruction is the *absence* of a positivity-printable floor — consistent.
- **IL-ii** (nonabelian): the `u=2` renormalon with the correct `β` coefficient is a nonabelian (asymptotically-free) statement; U(1) has no condensate.
- **IL-iii** (no small parameter): the renormalon ambiguity is `O(1)` relative at `α_s~O(1)` — the very reason it cannot be controlled.

---

## §6 · Honest bottom line

- **The roof is not laid. The direct route is blocked at a precise step:** the `u=2` IR renormalon makes "the gluon condensate" scheme-ambiguous at `O(1)` relative precision, so **"floor `⟨tr F²⟩` uniformly in `a`" is not a well-posed scheme-independent demand** `[NG]` for the condensate-number class.
- **This is a structural finding, not a failure of effort:** it tells us the charter's identification `⟨¼trF²⟩ = H1` is a **physical** identification, but the *rigorous target* must be recast on a scheme-independent object — the defect-condensation statement or the exact trace-anomaly identity.
- **The physical nonzero condensate is real** (nonperturbative, measured); what is missing is a **scheme-independent, `a`-uniform floor**, and that is blocked by the renormalon at its own size.
- **Clay = 0%; no QED.** The next roof piece is the `0`-form defect-condensation route (§4.1), which sidesteps the `u=2` ambiguity entirely.

*Refs:* `BASE-SYNTHESIS.md` (§5, the condensate `=H1`), `BASE-FIELDMEASURE.md` (renormalon wall), `PROFILE.md`;
witness `w7/roof_condensate_floor_witness.py`. Renormalon background: 't Hooft (1977); Beneke, *Phys. Rept.*
317 (1999) 1 ("Renormalons"); SVZ, *Nucl. Phys.* B147 (1979) 385.

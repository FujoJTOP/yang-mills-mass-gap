# M1-DECISION-RULE — the three-branch rule for the L=8 result, fixed before the number arrives

> **Purpose.** The L=8 run (N=2000, β=2.5, SU(2)) is the first high-statistics extension of the
> L=4/6 plateau `χ_flat=3.807±0.074, 3.760±0.072`. This note fixes, *before* seeing L=8, what each
> possible outcome means and the exact wording the papers may use. It prevents post-hoc storytelling.
> **Labels:** [T]/[C]/[NG]/[O]/[V]. No QED. Clay 0%.

---

## §0 · The reference numbers

```
   L=4:  chi_flat = 3.807 ± 0.074   (N=2000, therm=800, gap=10, tau=1.03)
   L=6:  chi_flat = 3.760 ± 0.072   (same protocol)
   Haar: 3.7198                     (exact beta=0 value)
   free: chi(L) ∝ L^2               (the acid test; (6/4)^2 = 2.25)
```

The L=4/6 ratio is `0.987`; the gauge-invariant `chi_P` ratio is `0.935`; the plaquette is stable
(`0.65479 → 0.65285`). So the plateau is not an artifact of one estimator.

---

## §1 · The three branches

Let `chi8 ± se8` be the L=8 result, and define the deviation from the L=6 value in units of the
combined error,
```
   z = (chi8 − 3.760) / sqrt(se8^2 + 0.072^2) .
```

### Branch A — plateau (`|z| ≤ 1.5`, or `chi8 ∈ [3.55, 3.97]`)

**Meaning.** The variance face is flat from L=4 to L=8, consistent with the Haar value `3.72`.

**Allowed wording (v3 §6.9 / v2 ch15):**
> At β=2.5 the k=0 susceptibility is `χ=3.807±0.074` (L=4), `3.760±0.072` (L=6) and
> `χ8±se8` (L=8); the L=4→8 ratio is consistent with one, and a free massless field would have
> grown by `(8/4)^2=4`. Growth is excluded at more than `5σ`; the plateau is consistent with the
> Haar value `3.7198`. This is evidence for a finite k=0 susceptibility at β=2.5 on the lattice; it
> is **not** a proof of H1, because the carrier is gauge-invariant only up to the flat-chart
> normalization, the continuum limit `a→0` is not taken, and the atom/isolation pair is untouched.

**Do not write:** “H1 is proved”, “the mass gap follows”, “the continuum limit is under control”.

### Branch B — growth (`z ≥ 2`, or `chi8 ≥ 4.2`)

**Meaning.** The variance grows beyond L=6; the plateau is a finite-size transient.

**Allowed wording:**
> The L=8 value `χ8±se8` is `zσ` above the L=6 value; the L=4/6 plateau does **not** persist. The
> variance face of H1 is therefore not settled at β=2.5; the earlier positive reading is withdrawn.
> A growing `χ(L)` is consistent with a diverging k=0 susceptibility, which would be a no-go for the
> variance route. The decisive run is L=12 with the same protocol.

**Do not write:** “the plateau is confirmed”, “the discrepancy is a fluctuation” (unless `z<2`).

### Branch C — inconclusive (`1.5 < z < 2`, or the error bar spans both regions)

**Meaning.** L=8 does not decide.

**Allowed wording:**
> The L=8 value `χ8±se8` does not separate the plateau from growth; the L=4/6 trend is neither
> confirmed nor excluded. The decisive run is L=12 with the same protocol; no claim about the
> variance face is made until then.

**Do not write:** “consistent with a plateau” (that is Branch A language), “growing” (Branch B).

---

## §2 · What gets updated in each branch

| file | A (plateau) | B (growth) | C (inconclusive) |
|---|---|---|---|
| `docs/paper-v3/main.tex` §6.9 | add L=8 to the plateau table; keep caveats | replace plateau sentence with the growth statement | add “L=8 does not decide; L=12 needed” |
| `docs/paper-v2/ch15-uvface.tex` | same | same | same |
| `research/v3/Q-WAVE-SYNTHESIS.md` §3 | update the 4D pilot row | mark the pilot as withdrawn | mark as undecided |
| `docs/yang-mills-frontier.md` §18/§21 | add the L=8 number | add the withdrawal | add the undecided status |
| `docs/YANG-MILLS-MASTER.md` §11 | update the 4D pilot line | same | same |

---

## §3 · Invariants (true in every branch)

- `χ_flat` is a flat-chart coordinate; it is gauge-invariant only infinitesimally. The gauge-invariant
  cross-check is `χ_P` / the Polyakov correlator, which must be reported alongside.
- The atom/isolation pair (`A/W/I/G`) is untouched by M1; the L=8 result says nothing about `f₁²`,
  `w₁`, or `ρ²`.
- No continuum limit is taken; no OS reconstruction is performed.
- **Clay stays 0% and no step is labelled QED in any branch.**

*Refs:* `Q5-4D-CHI.md`, `Q-WAVE-SYNTHESIS.md` §3, `M1-4D-RUN.md` (when written),
`research/v3/w7/m1_4d_run_results.json`.

# IMPOSSIBLE-CONSISTENCY — the gap is NOT forced by global consistency: the free massless gauge theory is a consistent, gapless 4D YM; pure YM has no 0-form anomaly to match; the 1-form anomaly is satisfied by the confined vacuum

> **VERDICT (no QED).** The "consistency bootstrap" — unitarity + analyticity + crossing + a-theorem +
> 't Hooft anomaly matching + conformal-collider positivity — does **not** force the 0⁺⁺ massless sector of 4D
> pure YM to be empty. **Decisive acid test:** the **free massless gauge theory (`g=0`) is a valid, gapless
> 4D YM satisfying every consistency axiom** — so no consistency functional can exclude the gapless case
> without excluding a valid theory. A functional that fires only on the gapped theory must use the
> interacting-specific operator content, i.e. the gap itself (`H1`). Labels **[T]/[NG]/[O]**.

## §1 · Why it fails (four exact reasons)

1. **[T] The free massless field is consistent and gapless.** At `g=0` pure YM is a free (conformal) theory:
   `ρ(s)≥0` (unitarity), the Stieltjes transform is Herglotz (`Im Π/Im z≥0`, analyticity), crossing is
   automatic for a self-conjugate `0⁺⁺`, the a-theorem inequality `a_UV>a_IR` holds, and it satisfies the
   conformal-collider/positivity bounds exactly (it *is* a CFT). It is **gapless** (`:tr F²:` two-particle
   continuum from `s=0`). ⇒ consistency holds on a gapless theory ⇒ **consistency cannot force the gap**
   (acid test: the consistency functional also fires on the gapless free field ⇒ not discriminating).

2. **[T] Pure YM has no 0-form (chiral) 't Hooft anomaly to match.** No fermions, no chiral symmetry ⇒ the
   standard anomaly-matching argument (which requires the IR to carry the UV chiral anomaly) is **empty**.

3. **[T] The only anomaly is the 1-form center `Z_N` anomaly.** At `θ=0` it is matched by the **confined
   vacuum's `Z_N` topological order** — a gapped TQFT. Anomaly matching is *consistent with* the gap; it does
   **not** force a massless 0-form (`0⁺⁺`) state. (`θ=π` forces 1-form deconfinement, a different sector.)

4. **[T] a-theorem is an inequality, not a gap-forcer.** `a_UV>a_IR` holds for *both* a gapless CFT IR
   (`a_IR>0`) and a gapped IR (`a_IR=0`). Conformal-collider/positivity/crossing are satisfied by the free
   massless field *and* by interacting CFTs.

## §2 · The honest nuance

The modern bootstrap **proves gaps only *given* spectral data / *assuming* the gap it studies** (it bounds
couplings given a gap, or pins OPE data in a CFT). It does **not** derive a gap from consistency alone —
because consistency is satisfied by gapless theories (free fields, and many interacting CFTs). Forcing the
gap needs an **interacting-specific** input = `H1`.

## §3 · Witness (`w7/impossible_consistency.py`, runs clean)

`ρ_free=θ(s)` (gapless continuum from 0) vs `ρ_gap=δ(s−m²)+continuum`: both have `min ρ≥0` (unitarity),
both Herglotz (`min Im Π/Im z = 0.073` vs `0.010`, both `≥0`), both satisfy `a_UV>a_IR`, in the crossing
cone ⇒ **not discriminated**. The consistency conditions hold on the gapless theory.

**Verdict:** `[NG]` consistency bootstrap; `[T]` free-massless consistency + no 0-form anomaly; `[O]` the
interacting-specific input is `H1`. Clay `0%`; no QED.

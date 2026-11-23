# Method 3 — Stein / SLC: coarse-staple neighbour regression

> CPB-attack probe. Verdict: **INCONCLUSIVE — the instrument fails.**

## What was run

Scripts left in this directory: `neighbour_regression.py` (Coulomb-gauge-fixed, R²
gauge-orbit invariant) and `transported_test.py` (explicit Ad-transport of the 18
staple links). L=6, Python312.

## Results

- **Cold-perturbation control:** `f = 4βξ` exactly at an isolated link ⟹ gradient /
  algebra coordinates are correct; the Coulomb fixer works (functional 1.5→420; plaquette
  0.588 at β=6).
- b=2 (coarse), raw coords: R²(same + 18-neighbour) ≈ **0.007** (β=2.6, 3.0).
- b=1 (fine control), raw: R² ≈ **0.001**.
- b=1, explicit Ad-transport: R² ≈ **0.004** (β=6, 3.0).

## Verdict: instrument fails — inconclusive

The **fine-scale control fails** (R²≈0) although the force is *exactly* a function of
the staple. The force is the gauge-covariant **curl**, nonlinear in the links; the
linear regression (even transported) explains ≈0 of its variance. Therefore:

- "R²≈0 ⟹ genuinely nonlocal ⟹ operator target" is **NOT supported**;
- "local-but-off-site" is **not supported** either.

The entry is **neither confirmed nor closed**.

## Right instrument

A **nonlinear** decomposition (exact plaquette) or a gauge-covariant **transported link
blocking** for b=2 — a larger build.

## Honest boundary

A decision input, not a proof; do **not** conclude the operator route from these numbers.
**Correction to ENTRY.md §4.2:** the earlier "force is nonlocal (R²≈0)" reading was an
artifact of the linear regression, not a physical statement.

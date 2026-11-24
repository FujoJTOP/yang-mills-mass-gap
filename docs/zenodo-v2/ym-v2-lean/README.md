# Lean 4 supplement — Yang–Mills v2

Companion source to:

> **A Conditional Reduction of the 4D Pure Yang–Mills Mass Gap — Second Edition:
> The Wall in Two Coordinates** (2026-11-24, 81 pp.)

**This is not a proof of the Clay problem.** The Lean development machine-checks the *logic* of
the conditional reduction and the H1 sharpening; the analytic hypotheses remain explicit
hypotheses. `Clay = 0%`; no step is labelled QED.

## Layout

```
YMStatement/           the campaign's statement modules (6 files)
  Reduction.lean         conditional reduction: clay_gap_conditional; (A)/(C) proved, (B) sorry hypotheses
  H1Sharpening.lean      H1 does not imply an atom; m >= lambda f (five theorems, 0 sorry)
  Handoff.lean           handoff statements
  RGHandoff.lean         RG handoff statements
  Conditional.lean       conditional-theorem skeleton with #guard_msgs axiom audit
  Basic.lean             imports the external OSforGFF package (see below)
research/
  handoff-proof/         HandoffMargin.lean, axcheck*
  glue-a-limits/         GL3Limits, GL4L2, GL4Transfer, GL9Decay
  glue-b-blocking/lean/  GL2a, GL2bPrime
  glue-c-gauge-group/    GL6Fubini, GroupUniform
  glue-d-dual-handoff/   GL7Dual
  lean-infra/            YangMillsInfra library (14 modules)
  minimal-handoff/       MinimalHandoff.lean
  nontriviality/         NonTrivialCore.lean, axiom checks
  os2-rotation/lean/     SYMDim4.lean
scripts/
  build-core.sh          builds Reduction + H1Sharpening + Conditional
  build-all.sh           builds every library
  original/              the campaign's original check scripts, kept for provenance
lakefile.lean            pinned Mathlib requirement
lean-toolchain           leanprover/lean4:v4.33.0-rc1
MANIFEST.sha256          checksums of the files in this supplement
```

## Pins

| Item | Pin |
|---|---|
| Lean toolchain | `leanprover/lean4:v4.33.0-rc1` |
| Mathlib | `79d0395a1825a6264ad5d269e35e60537518955e` |
| OSforGFF | vendored in the campaign repository; **not redistributed here** |

The pin matches the toolchain recorded in the paper's Appendix A. The package
`OSforGFF` (and its transitive packages `GaussianField`, `BochnerMinlos`,
`kolmogorov_extension4`) is third-party and is **not** included in this
supplement; `YMStatement/Basic.lean` therefore does not build from this folder
alone. Every other module builds against Mathlib only.

## Quick start

```bash
# requires elan and network access
./scripts/build-core.sh
```

What it does:

```bash
lake update
lake exe cache get   # if this fails, Mathlib is compiled from source (slow)
lake build YMStatement.Reduction YMStatement.H1Sharpening YMStatement.Conditional
```

Expected: exit 0. `Reduction.lean` and `H1Sharpening.lean` contain **no `sorry`**.
`Conditional.lean` carries `#guard_msgs` blocks that freeze the axiom footprint; its
`(B)` inputs are deliberate hypotheses, not proved claims.

## Expected axiom audit

The theorems advertised in the paper have the standard footprint only:

```
[propext, Classical.choice, Quot.sound]
```

for `YMStatement.Reduction.clay_gap_conditional` and for the five theorems of
`YMStatement.H1Sharpening` (`lemmaG_twoSided`, `cert_necessity`, `cert_sufficiency`,
`H1_blind_to_thin_atom`, `H1_does_not_imply_atom`). To print the footprints:

```bash
lake env lean YMStatement/H1Sharpening.lean
```

## Mapping to the paper

| Paper | File |
|---|---|
| Appendix A, reduction theorem `clay_gap_conditional` | `YMStatement/Reduction.lean` |
| Appendix A, H1 sharpening (five theorems) | `YMStatement/H1Sharpening.lean` |
| Chapter 5, conditional theorem and its gluing lemmas | `YMStatement/Conditional.lean` + `research/glue-*` |
| Chapter 3, `m >= lambda f` and the countermodel | `YMStatement/H1Sharpening.lean` |

## Honesty notes

- The analytic hypotheses of `clay_gap_conditional` (H1, CP, W_geo ≥ 0, P3,
  no-dissolution) are `sorry` in `Reduction.lean`; their intended content is named
  in the paper and in `docs/LEAN-STATUS.md` of the repository.
- The `(B)` section of `Reduction.lean` is explicitly stated, not proved.
- `YMStatement/Basic.lean` imports the external `OSforGFF` package and is excluded
  from the default build roots.
- `scripts/original/*.sh` contain absolute Windows paths from the campaign machine.
  They are included for provenance only; use `build-core.sh` / `build-all.sh`.
- Numerical claims in the paper are measurements, not Lean theorems; no numerical
  result is formalized here.

## License / provenance

These files are part of the Zenodo record and of the public repository
`github.com/FujoJTOP/yang-mills-mass-gap`. The vendored `OSforGFF` package is
Apache-2.0 and is not redistributed in this supplement.

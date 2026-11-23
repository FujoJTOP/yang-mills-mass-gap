# Zenodo submission metadata

> Companion file: `CPB-preprint.md` (English, submission version).
> The Chinese draft `CPB-preprint-zh.md` is retained locally and is **not** part of the submission.

---

## Upload type

```
Preprint
```

## Title

```
Crossover Preservation and the Weak–Strong Handoff for Four-Dimensional
Pure Yang–Mills: A Reduction and a Map of Obstructions
```

## Authors

```
Jiang, Yuxuan
```

(Affiliation / ORCID optional; fill in the Zenodo form only, not in the manuscript body.)

## Publication date

```
2026-11-22
```

## Description (abstract)

```
We reduce the existence and mass-gap problem for four-dimensional pure Yang–Mills to a
weak–strong handoff inequality q0 ≥ q1, formulate the two-sided crossover preservation
statements CPB-W and CPB-S that would discharge it, and prove the handoff algebra in a
machine-checked (Lean 4) form. We then report, with proofs and falsification tests, that
every mechanism currently available for CPB-W/CPB-S fails at the crossover q ≈ 1:

  * Bianchi/sign cancellation is absent: dF = 0 is a topological (flux 2-cocycle)
    statement acting on cluster entropy, not on single-block weights; quadratic
    scaling is empty at q = 1; and κ_c = 1 is provable in the strong-coupling
    representation, with no additive defect charge in 4D non-abelian gauge theory.
  * Gaussian-weighted class preservation is blocked by a cluster constant
    (C_cum ≥ 2 against a required ≤ 1.6).
  * Numerically, κ_c is not directly measurable, the σ defect-fraction proxy is
    shown to be a Gaussian tail of a blocked field that is itself Gaussian, and
    ξ/a ≈ 7 is a feasibility boundary.

The outcome is not a proof of the mass gap. It is a precise statement of the
obstruction: at q ≈ 1 the single-block activity must be small, either because the
defect action is large or because signs cancel, and no known structure supplies
either. The paper also records the surrounding conditional theorem: ten named inputs,
nine glue lemmas, nine machine-checked blocks, and four explicit external declarations
(minimal assumptions are pushed into type signatures rather than hidden in prose).
```

## Keywords

```
Yang--Mills
mass gap
lattice gauge theory
renormalization group
crossover
Bianchi identity
Lean 4
formal verification
numerical falsification
```

## License

```
CC BY 4.0
```

## Language

```
eng
```

## Related identifiers

```
is supplemented by   <repository / code archive URL>   -- optional
```

## Files to upload

1. `CPB-preprint.md` — manuscript (English).
2. (Recommended) `CPB-preprint.pdf` — convert from Markdown.
   - `pandoc` / `xelatex` are **not installed** on this machine.
   - After installing:
     ```
     pandoc CPB-preprint.md -o CPB-preprint.pdf --pdf-engine=xelatex
     ```
3. (Optional) code and data archive: `lean/`, `research/`, `docs/conditional-theorem.md`.

## Pre-submission checklist

- [ ] Author filled as `Jiang, Yuxuan`
- [ ] Language set to `eng`
- [ ] Numbers in the abstract agree with the manuscript
      (`t_up^{2L}(6)=2.639`, `q₀=0.571`, `λ≈7.28`, `κ_c=1`, `C_cum≥2`)
- [ ] The three references marked *(to be verified)* are checked or removed
- [ ] Disclaimer retained ("does not claim a proof")
- [ ] AI-assistance statement retained (Deepseek v4.1 flash, deep participation)
- [ ] If code is uploaded: no personal data, no keys, no absolute paths
- [ ] Record the DOI after publication

## Notes

- The hexagonmath citation is written as `hexagon:2611.00492v1` (a site ID), **not** an arXiv number; `arXiv:2611.00492` is an unrelated paper.
- The 2026 OpenAI construction is cited as an unpublished technical report assessed by an independent committee, not as a peer-reviewed paper.
- The Chinese draft is not part of the submission ("English only").

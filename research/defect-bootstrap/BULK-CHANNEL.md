# Does the defect-4pt bootstrap reach the bulk gap?

> CPB-attack probe (sub-task 3 of round 3). Verdict: the link is **indirect and conditional** —
> a probe, not a proof.

## Bulk-channel decomposition

The defect 4-point (four operators on the line) has *two* OPE channels: the defect channel
(defect operators, hypergeometric blocks) and a **bulk channel** — gauge-invariant bulk
operators `O_b` entering the defect–defect OPE,
`O(τ)O(τ') ⊃ c̃ |τ−τ'|^{Δ_b−2Δ_O} O_b(x)`, with `c̃ ∝` the bulk one-point function in the
defect background. Crossing:

```
Σ_def c² 𝔽(χ) = Σ_def c² 𝔽(1−χ) + Σ_bulk c̃² 𝔾(χ) .
```

The bulk channel carries `tr F²` (`Δ=4`), `T_{μν}`, `tr F⁴` (`Δ=8`), …; the lightest is the
`0^{++}`.

## Logical link

The bootstrap constrains `Δ_b` only for bulk operators with a **nonzero bulk–defect OPE**. The
gap is the lightest bulk-scalar dimension, so the link holds **iff** (a) the lightest bulk state
couples to the defect — plausible (free Maxwell: `tr F²` couples to the displacement `D`),
unproven; and (b) the defect-4pt exclusion is strong enough to pin the bulk spectrum.

## Honest boundary

The link is **indirect and conditional**: excluding interacting *defect* solutions is *weaker*
than excluding an interacting *bulk* (the bulk could be interacting while the defect sector is
effectively free). The defect-4-point bootstrap is a **probe, not a proof**, of the bulk gap;
the 3-cross-ratio bulk–defect block sees the bulk operator directly and is the cleaner object.

Sources: arXiv:2501.06900 (defect-4pt blocks); arXiv:1601.02883 (defect CFT, bulk–defect OPE).

# Gap ↔ defect-spectrum dictionary

> CPB-attack probe (round 4). Verdict: the link is precise — theorem (a) + **dynamical assumption
> (b)**; no symmetry forces (b).

## Dictionary

For `φ = tr F²`, define the `φ`-channel gap `Δ_φ := inf{Δ_b : C^{φφ}_b ≠ 0}` (the lightest bulk
primary in the `φφ` OPE). The bulk–defect block `⟨φφ⟩_def`'s **bulk-channel** bootstrap bounds
exactly `Δ_φ`. The physical gap is `Δ_bulk = inf over all bulk primaries`. The defect spectrum
`{Δ_i}` (operators on the line) is **defect** data and does **not** bound the bulk gap directly;
the bound comes only from the bulk-operator channel of the bulk–defect correlator.

## Precise implication

`Δ_bulk ≥ δ` is established iff:

- **(a)** the bootstrap gives `Δ_φ ≥ δ` — rigorous up to interval arithmetic, and the OPE
  coefficient is automatic for an interpolating `φ`;
- **and (b)** `Δ_bulk = Δ_φ`, i.e. the lightest bulk state lies in the `φ = tr F²` channel (the
  `0^{++}` glueball). If the gap sits in a defect-invisible channel, the bound does not transfer.

## Status

(a) is a **theorem** (bootstrap rigor). (b) is a **dynamical assumption** — true on the lattice
(`0^{++}` lightest), **not** symmetry-protected. So the combined implication is **conditional**,
exactly as `BULK-CHANNEL.md` said, with the condition now explicit: *the gap channel must couple
to `φ`*.

## Honest boundary

No symmetry (`Z_N`, anomaly) forces (b); it is a physical input. The 3-cross-ratio block sees the
bulk operator directly and is cleaner than the defect-4-point's indirect bulk channel.

Sources: arXiv:2501.06900; arXiv:1601.02883; parent `BULK-CHANNEL.md`.

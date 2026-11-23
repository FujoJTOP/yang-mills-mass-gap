# Explicit SDP form (corrected)

> CPB-attack probe (round 5, reuse). Verdict: the **untwisted** 1-form constraint is **not** a
> new solver — it is **standard SDPB + a per-operator phase rotation + one linear equality**.

## Correction to the first report

The *untwisted* anomaly is **cheaper** than first said — not a shrunken positive cone, just a
**per-operator phase rotation**. Ray-positivity, not a new cone.

## Variables

Functional `Λ` (the SDPB object). The anomaly fixes a **phase `φ_i` per defect operator** (from
`ω(q_i,q_j)=e^{2πi k q_i q_j/N}`).

## SDP form

```
max  Λ(z)
s.t. Re(e^{−2iφ_i} Λ(F_i)) ≥ 0      for every allowed (Δ_i, j_i)   [ray positivity]
     Λ(G_O) ≥ 0                       for every bulk O                [real]
     Re Λ(crossing) = 0 ,  Im Λ(crossing) = 0   (per Z_N sector q)
     Λ(F_identity) = 1
```

## Solver

**Standard SDPB + a fixed per-operator phase rotation** of each block column
(`F_i ↦ e^{−2iφ_i}F_i`) plus **one extra linear equality** (`Im Λ(cross)=0`). Implementable as
preprocessing; **no new positivity cone, no new solver**.

## Cost

~`N ×` standard SDPB (one crossing per `Z_N` sector); phase rotation `O(#ops)`. Tractable for
small `N` (SU(3): ~`3×`).

## Genuinely new

The **twisted (line-changing) sector** — a two-line / domain-wall SDP, which *is* a larger solver
(see `TWISTED-2.md`: the line-changing-local reading is also standard SDPB; the point-defect-on-
line is the hard one).

## Honest boundary

`φ_i` depend on the level `k` / BF normalisation; untested.

Sources: GKSW arXiv:1412.5148; Simmons-Duffin SDPB arXiv:1502.02033; arXiv:2501.06900.

# SDP setup for the simplified defect bootstrap

> CPB-attack probe (sub-task 9). Verdict: the SDP is writable, but **the defect blocks do
> not exist**, so it cannot be instantiated; and with no fusion layer `Z_N` is only a
> grading — likely too weak.

## SDP form

Crossing: `Σ_i a_i² 𝔽_i(ξ) = Σ_O C²_O 𝒢_O(ξ)`, `ξ∈[0,1]` (defect cross-ratio). Positivity:
`a_i² ≥ 0`, `C²_O ≥ 0`; spectra on the unitarity bounds — defect `(Δ_i, s_i)` in the
`SO(2,1)` representation, bulk `Δ_O ≥ ℓ+1` (`ℓ≥1`), `≥1` (`ℓ=0`). `Z_N` grading: defect
operators carry `q∈Z_N`; the anomalous phase enters `𝔽_i` as a **fixed** twist, not a
variable.

## SDPB input

A vector of block functions sampled at points `ξ_j`; search a linear functional `Λ` with
`Λ(𝔽_identity)=1` and `Λ(𝔽_O) ≥ 0` for all allowed `(Δ_i, s_i)`. SDPB (Simmons-Duffin,
arXiv:1502.02033) solves the PSD program over the `(Δ,ℓ)` blocks.

## Island vs proof

An *island* = a bounded region of allowed CFT data, rigorous up to interval arithmetic.
"Only free" = exclusion of **all** interacting solutions — a **global** uniqueness statement
over a continuum, not a finite set of point-wise SDPs.

## Honest boundary

(i) With **no fusion layer**, `Z_N` is only a grading — likely too weak to exclude interacting
solutions. (ii) The **4D non-SUSY defect blocks do not exist**, so the SDP cannot even be
instantiated yet. (iii) A global "only free" needs an **analytic** bound, not a finite scan.

Sources: Simmons-Duffin, SDPB, arXiv:1502.02033; Poland–Rychkov–Simmons-Duffin, arXiv:1805.04405.

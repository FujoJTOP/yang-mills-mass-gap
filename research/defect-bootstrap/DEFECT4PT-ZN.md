# Defect 4-point bootstrap for pure-gauge `Z_N` (the shortcut)

> CPB-attack probe (sub-task 1 of round 3). Verdict: **runnable, but blind to the `Z_N`
> grading and the anomaly for neutral operators**; yields at most bulk-operator bounds.

## Crossing

Four defect operators on the line; invariant `χ = τ₁₂τ₃₄/(τ₁₃τ₂₄)`. With equal external
scalars,

```
G(χ) = Σ_k C₁₂^k C₃₄^k G_k(χ) = Σ_k C₁₄^k C₂₃^k G_k(1−χ)   (s↔t),
G_k(χ) = χ^{Δ_k} ₂F₁(Δ_k, Δ_k; 2Δ_k; χ)   — the 1D block already computed to 1e-15.
```

Blocks available; the only constraint is `χ ↔ 1−χ`.

## `Z_N` adaptation

Defect operators carry `q_i∈Z_N`; the OPE preserves `q_k = q_i + q_j mod N`. For **neutral**
external operators (`q_i=0` — e.g. the displacement `D`, adjoint / N-ality 0), both channels run
over the `q=0` sector: **the grading is trivial.** To see the 1-form structure one needs
**line-changing (twisted) defect operators** (`q→q'`), which are not the standard defect
4-point.

## Anomaly phase

`ω(q_i,q_j) = e^{2πi k q_i q_j/N}`; for `q_i=0`, `ω=1` — also **blind**.

## Bulk channel → bulk gap?

Via the bulk–defect OPE, the `D·D` bulk channel contains the bulk stress tensor and `tr F²` (the
`0^{++}` creator), so it *can* constrain `Δ(tr F²)` — but only as a **bound**, and only if
`D·D → tr F²` has nonzero coefficient.

## Verdict

Runnable, but **blind to both the `Z_N` grading and the anomaly** for neutral operators; it
reduces to the standard line-defect bootstrap (matching the SDPB warning "too weak"). The
`Z_N`-based exclusion needs the **twisted sector** or the 3-cross-ratio bulk–defect block. The
shortcut yields at most bulk-operator bounds.

Sources: arXiv:2501.06900; arXiv:2208.11715; GKSW arXiv:1412.5148.

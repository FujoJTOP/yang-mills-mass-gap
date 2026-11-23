# `Z_N` twisted (line-changing) defect sector

> CPB-attack probe (round 4). Verdict: the twisted sector is **where the 1-form symmetry is
> visible** — but it is a **selection rule**, not the dynamics, and not the bulk gap directly.

## Twisted operators

A line of charge `q∈Z_N` carries a defect Hilbert space `H_def = ⊕_{q_i∈Z_N} H^{(q_i)}`.
Ordinary (untwisted) defect operators change `q_i→q_i+q_j`. A **line-changing (twisted)
operator** `O_{q→q'}` sits at a point where the line's total charge jumps `q→q'`; its charge
label is `Δq = q'−q ∈ Z_N`, and it is what the ordinary defect OPE **cannot** see.

## Crossing

The twisted 4-point `⟨O_{q₁→q₂} O_{q₂→q₃} O_{q₃→q₄} O_{q₄→q₁}⟩` crosses with the **anomaly
2-cocycle** `ω(q_i,q_j) = e^{2πi k q_i q_j/N}` (`k` = level, θ-related, from the 5d BF symmetry
TQFT). Unlike the neutral 4-point, where `q_i=0` gives `ω=1`, the twisted sector has `q_i≠0`,
so **`ω≠1` and the crossing carries a non-trivial projective phase**. This is the only place the
`Z_N` 1-form structure enters the defect crossing: the phase is an **irremovable 2-cocycle**,
fixed by the anomaly, not a free parameter.

## `Z_N` visibility

**Yes** — the twisted sector is where the 1-form symmetry is **visible**. The neutral
(bulk-scalar) 4-point is blind (`ω=1`, trivial grading); the twisted 4-point sees the grading
and the anomalous phase. So the defect-4pt shortcut (which used neutral externals) was blind
**by construction**, not by accident.

## What it constrains

The twisted crossing constrains the `q_i`-graded defect spectrum and the relative OPE phases —
a **selection rule / anomaly-matching constraint**, not the dynamics (`Δ_i`, OPE coefficients).
It can **exclude** defect spectra inconsistent with the projective `Z_N` representation, but does
**not** by itself give the bulk gap. Whether a gapless phase matching the twisted sector exists
is the same open classification step (iii).

## Honest boundary

Requires the full 1-form-constraint layer / symmetry-TQFT machinery (the SDPB "needs-new-solver"
gap); the toy (3D codim-2) analogue exists, the 4D pure-gauge twisted bootstrap does not.

Sources: GKSW arXiv:1412.5148; Cordova–Tachikawa (symmetry TQFT); arXiv:2501.06900;
arXiv:2208.11715; Etingof–Nikshych–Ostrik (fusion categories).

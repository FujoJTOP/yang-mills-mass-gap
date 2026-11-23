# SDP machinery de-risked on a known 1D toy

> Toy: 1D CFT, one real scalar (Δφ=1/2), Z2. Known answer = generalized free field.
> Code: `sdpb_toy.py` (numpy 2.4.6, scipy 1.18.1; **no cvxpy** — scipy NNLS only).
> Run: `python sdpb_toy.py`. Real output below.

## Results (actual numbers)

| check | result |
|---|---|
| [1] blocks + crossing `max|G(u)−G(1−u)|` | `4.3e-14` |
| [2] positivity cone: GFF nonneg fit residual | `4.2e-16` |
| [2] recovered OPE `a²[0:4]` | `[2.0, 0.3333, 0.0286, 0.0022]`, all ≥0 |
| [3] gap scan: best-fit Δ_gap | `1.000` (=2Δφ), residual `3e-16` |
| [3] drop lightest op → residual | `0.569` (**excluded**) |
| [4] Z3 phase rotation `|ω|` | `1.000000000000000` |
| [4] Im(cross): trivial q=0 / one-form q=n | `0.0` / `4.31` (nonzero) |

The sharp zero of the residual at `Δ_gap=1.0` reproduces the known free-theory
gap `Δ_min=2Δφ=1.0`; fits are grid-stable (argmin 1.000 for u1∈[0.7,0.95], N∈[25,40]).

## What this validates (plumbing)

- 1D blocks `u^{Δ−2Δφ}₂F₁(Δ,Δ;2Δ;u)` and crossing vectors assemble correctly.
- Positivity cone (a²≥0, NNLS) contains the known GFF solution; OPE data recovered.
- Gap/exclusion arithmetic reproduces a known bound within truncation tolerance.
- Z_N per-operator phase rotation: unit modulus ⇒ ray-positivity cone unchanged;
  trivial sector stays real (standard SDPB), the one-form sector adds exactly one
  `Im Λ(cross)=0` linear equality (Im term is nonzero). Matches `SDPB-FORM.md`'s claim.

## What it does NOT validate

- **Only plumbing, not YM physics.** No 4D block, no defect kinematics, no pure-gauge Z_N.
- **No trustworthy dual bound.** A dual-functional LP (Λ(R₀)=1, Λ(R_Δ)≥0) was tried; its
  threshold is grid/dof-dependent (2.27 on one grid, none on a finer one) — a point-eval LP
  is not a faithful SOS/PSD test. The real PSD/SOS layer (SDPB) is absent. The `[3]` result is
  primal known-answer recovery, not a dual extremal bound.

**Verdict:** end-to-end SDP plumbing works on a known case; a real exclusion bound needs SDPB.

# Can the bootstrap close `(iii)`? — uniqueness vs islands

> CPB-attack probe (the crux). Verdict: **NO.** The bootstrap **cannot close `(iii)`** — no
> analytic uniqueness theorem exists for 4D CFTs; it gives at most a **partial exclusion**.

## (a) Upgrade path

An SDPB island is a rigorous **bound** (interval arithmetic) at each sampled point of spectrum
space — **not** a uniqueness statement. Excluding *all* interacting CFTs needs one of:

1. an **analytic extremal functional** — a closed-form inequality saturated by the free solution
   and violated by every interacting one (rare; exists only in analytic / large-N limits);
2. an **a priori compactness reduction** (prove the spectrum lies in a compact finite-dimensional
   set via a gap bound + cutoff), then a rigorous interval-arithmetic scan — but the reduction is
   itself an analytic input;
3. a **bootstrap+theorem hybrid**, showing the island's boundary collapses onto the free point.

A finite scan over sampled `(Δ_i, j_i)` covers only a discrete set; it gives point-wise
infeasibility, **never** global exclusion over the continuum. **There is no "more SDPB points"
route to uniqueness.**

## (b) Precedent

**None for 4D CFT *uniqueness* by bootstrap.** 4D bootstrap yields **islands** (`N=4`,
stress-tensor bootstrap) and analytic bounds (extremal functionals), but **no 4D CFT has been
proven unique** this way. A discrete `Z_N` symmetry is a *selection rule* (trims the spectrum),
**not** a uniqueness engine.

## (c) Strongest honest statement

A rigorous **partial exclusion**: "no interacting CFT with the pure-gauge `Z_N` structure whose
lightest operator in channel `X` lies below a computable `δ_bootstrap`." Strictly weaker than
`(iii)`: an interacting window **above** `δ_bootstrap` survives. Conditional on the bulk gap
(IN-b), the bootstrap bounds the defect spectrum — **the wrong direction** to prove the gap.

## Honest boundary

`(iii)` requires an analytic uniqueness theorem **absent for 4D CFTs**; the bootstrap can only
**bound**, not close `(iii)`. **This confirms the route's endpoint is a conditional theorem
(IN-b), not an unconditional one.**

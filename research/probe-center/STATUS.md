# STATUS — center-vortex probe (NP-7)

> 2026-11-23 · 4D pure SU(3) lattice, Z(3) maximal-centre gauge · full: `REPORT.md`

## Verdict

**`κ_c < 1` cannot come from center vortices for the object `CPB-W` needs.
NP-7 closed as a route to the mass gap.** (Mechanism-insufficient / no-known-
route, per `SURVEY.md` §1.2 — not "proved false".)

| test | result |
|---|---|
| projected **adjoint** loop | `= 1` to machine zero (`max|W_A−1|=0`), all β,L,R,T |
| projected adjoint Polyakov correlator | `~1e-16 (≡0)` |
| projected **fundamental** loop vs full (β=6) | reproduces area law (centre dominance, `χ_proj/χ_full ≈ 1.5–1.7` at (2,2)) |
| full adjoint Creutz (β=6) | `0.589 ≈ (9/4)·χ_F` — Casimir scaling, carried by non-centre configs |
| signed/absolute D1 (real phase) | `−0.5000` **constant**, every cluster/size |
| signed/absolute D2 (orientation) | `∈ [−0.004, +0.005] ≈ 0` (CP symmetry), size-independent |
| vortex ensemble | **percolating** (one spanning cluster) at all β = 1.7–6.0, L=6–10 |
| density (β=2.6–3.3, L6) | 0.56–0.59; L6 vs L8 agree <1% |

**Why `|D1|=0.5<0.8` is not a pass**: it is `|cos 120°|`, size-independent, and
applies to *fundamental* (charged) loops; a `CPB-W` block weight is neutral and
sees factor `1`.

## What was computed

* MCG (= adjoint-Landau functional `Σ|Tr_F U|²`) by Riemannian ascent, 2
  restarts; Z(3) projection; vortex plaquettes + link-adjacency clusters.
* Vortex density, cluster sizes, `k=1`/`k=2` counts, `z_signed/z_abs` under 3
  definitions, on 14 (β,L) points incl. L=8 β=2.6/3.0 generated here.
* Fundamental & adjoint Wilson loops (full + projected), fundamental & adjoint
  Polyakov correlators, Creutz ratios.

## Files

`mcg.py` `run_probe.py` `analyze.py` `gen_L8_warm.py` `probe_*.npz`
`runs_L8_b*.npz` `plots/` `analyze_out.txt`

## Not done / caveats

* Effective-action `κ_c` not directly measurable (STATUS §9 of `lattice-cpb`);
  probe shows the only sign channel (Z(3) phase) is neutral-blind.
* β=2.6–3.3 Creutz unstable (coarse); clean sector separation quoted at β=6.
* Snapshot errors statistical only (autocorrelation not corrected).

# SAVVIDY-LATTICE-CELL — empirical test of the softened χ-floor in the negative-mode region

> **VERDICT `[V-num]`.** Using EXISTING configs (no new MC), the S−D defect `D/S` in the
> negative-mode-active region is **BOUNDED (≤~1.1, no divergence)** — consistent with the *bounded
> softening* claim — **but it is O(1) (≈0.94–1.00), so `G″ = S−D` is a near-cancellation (≈0), NOT a
> robustly positive floor.** At the softest measured mode `k=2π/L` the floor is positive-but-thin
> (`G″=+0.088`, `D/S=0.937`); at `k=4π/L` it is `G″≈0` (`D/S=1.001`); **5/12 individual modes have `G″<0`**.
> **The actual `k→0` mode (`n=0`) is NOT in this dataset**, so the `k→0` positive floor is *not*
> established here. The negative mode is ubiquitous (100% of configs carry a negative plaquette at
> β=2.6). Clay 0%, no QED.

---

## §0 · Data and exact identity

Two existing files (no MC):

- `research/unconditional/w7/k0_su2_b2.4_L8_n2.npz` — SU(2) β=2.4, L=8, N=400, per-config per-mode
  `a=⟨v,x⟩`, `H=d²S/ds²`, `G=dS/ds`, `P1=d(−logJ)/ds`, `P2=d²(−logJ)/ds²`, at `k=2π/L (n=1)` and `4π/L (n=2)`.
- `research/lattice-cpb/runs/L6_b2.600_s51.npz` — SU(2) β=2.6, L=6, negative-plaquette "defect" series.

Exact identity (`k0_continuum_run.py` docstring):
```
G″(0) = (E[H] − Var(G)) + W ,   W = E[P2] − Var(P1) − 2·Cov(G,P1)
⇒ S = E[H] ,  D = Var(G) − W ,  G″ = S − D ,  D/S = defect fraction.
```

## §1 · Per-mode and aggregate results (β=2.4, L=8)

| k | S=E[H] | VarG | W | **G″=S−D** | D | **D/S** |
|---|---|---|---|---|---|---|
| `2π/L` (n=1) | 1.4042 | 1.3296 | +0.0133 | **+0.0879** | 1.3164 | **0.9374** |
| `4π/L` (n=2) | 1.4045 | 1.4367 | +0.0315 | **−0.0008** | 1.4053 | **1.0006** |
| ALL | 1.4043 | 1.3839 | +0.0227 | **+0.0432** | 1.3611 | **0.9692** |

Per-mode spread: `G″ ∈ [−0.144, +0.237]`, `D/S ∈ [0.831, 1.103]`; **5/12 modes have `G″<0` (`D/S>1`)**.
So `S` and `D` are both O(1.4); `G″` is their small difference — a near-cancellation, sign-fluctuating.

## §2 · Negative-mode ubiquity (β=2.6, L=6)

`frac(configs with a negative plaquette) = 1.000`, `plaq_min mean = −0.438`, `defect(thr) mean = [0.299, 0.421, 0.547, 0.704]`
for `thr = [0, 0.1, 0.2, 1/3]`. The Savvidy negative mode is **not rare** — it is present in every config.

## §3 · Interpretation (honest)

- **Bounded: YES.** `D/S ≤ ~1.1` everywhere — no divergence — consistent with the softened floor's
  "the negative mode does not destroy the floor by divergence".
- **"<1 and positive": NO (not uniformly).** `D/S ≈ 1` (0.94–1.00 aggregate), and exceeds 1 for 5/12
  modes. The floor is *softened to ≈0*, i.e. `G″` is a near-cancellation of two O(1.4) numbers, not a
  robust positive floor.
- **`k→0` not tested.** Softest mode is `k=2π/L`; the `n=0` (`k=0`) mode is absent. The positivity at
  the *softest finite k* (`G″=+0.088`) is real but thin and is the standard near-cancellation regime.
- **Volume lesson.** The corpus's earlier `G″≈+0.28` at β=2.5, L=4 is consistent with a **small-volume/
  kinematic** positive (`4sin²(π/L)` box gap); at L=8, β=2.4 it drops to `+0.04` — i.e. the "positive
  `G″`" weakens with volume, exactly the kinematic-artifact signature, not a genuine IR floor.

**Bottom line `[V-num]/[C]/[O]`.** The empirical defect is **bounded but O(1)**: `G″≈0` at finite k,
sign-fluctuating, with `D/S≈1`. This *supports* "no divergence" but *does not support* a robust
positive χ-floor at these modes; the `k→0` positive floor (the requirement) is **not measured here**
and would need the `n=0` mode, larger volume, and higher statistics. Far from a gap proof; consistent
with the `χ-floor ≈ near-cancellation` picture. **Clay 0%; no QED.**

*Witness:* `research/unconditional/w7/savvidy_lattice_cell.py` (runs clean, `py -3.12 -I`, numpy only).

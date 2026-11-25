# M1-FAST-CROSSCHECK — a Numba heat-bath version of the 4D run: L=12 in 233 s, `χ=4.16±0.20`, consistent with the slow L=12 `3.62±0.20` within 1.9σ

> **Purpose.** The M1 heat-bath run (`m1_4d_run.py`) is pure NumPy and takes ~58 min at L=8 and
> several hours at L=12. This note records a Numba-accelerated version (`m1_fast.py`) that runs
> **L=12, N=300 in 233 s** on the same measure, and its cross-check against the slow result.
> **Labels:** [V]/[O]. No QED. Clay 0%.

---

## §0 · Result

| run | algorithm | L | N | `χ_flat` | wall |
|---|---|---|---|---|---|
| M1 slow | NumPy checkerboard heat bath | 4 | 2000 | 3.807±0.074 | 221 s |
| M1 slow | same | 6 | 2000 | 3.760±0.072 | 958 s |
| M1 slow | same | 8 | 2000 | 3.823±0.069 | 3473 s |
| M1 slow | same | 12 | 300 | 3.615±0.198 | (feasibility) |
| **fast** | **Numba sequential heat bath** | **12** | **300** | **4.160±0.200** | **233 s** |

The fast and slow L=12 points differ by `0.545`, i.e. `1.9σ` on the combined error
`sqrt(0.200²+0.198²)=0.281`. Both are far below the free-massless expectation from L=8:
`3.823 × (12/8)² = 8.60`. **No growth is seen in either.**

Validation at small L (fast, N=200–300):
- L=4 (therm=200, gap=2): `4.08±0.26` (reference 3.807±0.074);
- L=6 (therm=200, gap=5): `3.62±0.18` (reference 3.760±0.072).

Both are within ~1σ of the slow values, so the fast kernel is sampling the same measure.

---

## §1 · What was optimized

The slow script spends its time in the NumPy `staple_sum2`/`su2_heatbath` path: large-array
`np.roll` and batched 2x2 matmuls, repeated for every checkerboard/`mu`. The fast script:

1. **Numba `@njit`** on the whole sweep, with explicit 2x2 complex arithmetic (no `np.roll`,
   no batched temporaries).
2. **Scalar exact SU(2) heat bath** per link (`sample_w0` rejection sampler), the same conditional
   distribution as the slow code.
3. **No gauge fixing, no glueball fits** in the timing run — only `χ_flat` and the Polyakov
   susceptibility. Those are separate diagnostics and are not needed for the plateau test.
4. **Sequential site/mu order** instead of checkerboard; the heat bath is exact, so the order does
   not change the stationary measure.

Timing: L=12 has `4·12^4 = 82944` links/sweep; 1700 sweeps at ~0.137 s/sweep = 233 s.
L=8 at the same settings would be ~3.16× cheaper (~74 s).

---

## §2 · What this does and does not establish

- **Does:** an independent, faster implementation reproduces the flat-chart plateau at L=12 within
  ~2σ of the slow run, and excludes free-field growth (`8.60`) by a wide margin. The Numba kernel
  is a reusable fast engine for future L=12/16 scans.
- **Does not:** change the carrier caveat. `χ_flat` is a flat-chart coordinate; the gauge-invariant
  Polyakov cross-check still drops (`χ_P=2.77` at L=8), and `χ_Landau` still grows with Gribov
  spread. The fast run does **not** measure the atom/glueball; that remains `[O]`.
- **Does not:** prove H1. It is a numerical cross-check at fixed β=2.5, not a uniform-in-`a` theorem.

*Refs:* `M1-4D-RUN.md`, `Q5-4D-CHI.md`, `m1_fast.py`, `m1_fast_results.json`,
`m1_4d_run_results.json`.

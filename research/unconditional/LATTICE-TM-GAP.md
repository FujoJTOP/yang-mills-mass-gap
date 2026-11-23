# LATTICE-TM-GAP — direct transfer-matrix gap of 4D SU(2) pure gauge

**Verdict.** Measured directly, the first-excited gap Δ of the Euclidean transfer matrix is **positive at every point** (β=2.2/2.5/2.8, L=4/6): Δ ≈ 1.0 in lattice units, never approaching 0. The spectrum behaves like a positive mass, not a vanishing gap — on these lattices. But a finite-box gap is an *upper* bound on Δ_{R⁴} (P4 §2); positivity at L≤6 confirms "gapped on every box tested," it does not prove the infinite-volume gap.

## Method
SU(2) Wilson action, S=βΣ_p(1−½Re tr U_p), β=4/g². Exact SU(2) heatbath (κ=β/2), periodic, β=2.2/2.5/2.8, L=4 (N=2×10⁴ configs) and L=6 (N=6×10³), ≥1000 thermalizing sweeps, measured every 4 sweeps. Operator = the gauge-invariant **0⁺⁺ magnetic/tr F²** density summed over each time slice, with spatial APE smearing (α=0.5; levels n∈{0,4,8} and {0,8,16}). I form the connected Euclidean correlator C(t)=⟨O(t′)O(t′+t)⟩_c (periodic time) and extract Δ from (i) the cosh effective mass m_eff(t)=arccosh[(C(t−1)+C(t+1))/2C(t)] of the highest-smearing operator and (ii) the leading eigenvalue λ₀(t) of a multi-operator GEVP, C(t)w=λ C(t₀)w, via the same cosh formula. Errors: jackknife over τ_int-binned configs. Validated: smeared links stay SU(2) to 10⁻¹⁵, operator gauge-invariant to 10⁻¹⁴, single-operator GEVP exactly reproduces m_eff, plaquettes match known β-values.

## Results (Δ in lattice units, plateau at t=L/2)
| β | L=4 | L=6 |
|---|-----|-----|
| 2.2 | 1.13(8) | 1.4(4)† |
| 2.5 | 1.15(20) | 1.02(15) |
| 2.8 | 1.10(30) | 1.05(15) |

†L=6,β=2.2 correlator decays too fast to reach a clean plateau (t=3 signal-to-noise <2); the t=2 value is an upper estimate. L=4 errors include a ±0.1–0.3 smearing systematic (still descending slightly with added smearing); L=6 β≥2.5 plateaus are flat in t (e.g. β=2.5: m_eff=0.96,1.05 at t=2,3).

## Reading
- **Gapped, not gapless.** Δ>0 at all six (L,β); nowhere does it collapse toward 0. The lightest 0⁺⁺ mass is O(1) in lattice units throughout.
- **Volume trend.** From L=4→6, Δ drops modestly (β=2.5: 1.15→1.02; β=2.8: 1.10→1.05) and, at the confined point β=2.2, is volume-stable (~1.1). A monotone, saturating-in-appearance decrease is exactly the finite-volume *upper-bound* trend; L=4,6 alone cannot separate a positive bulk limit from a slow vanishing.
- **Finite temperature.** The temporal extent is a temperature T=1/(La). SU(2) deconfinement (Wilson action) is β_c≈2.30 at N_τ=4 and ≈2.43 at N_τ=6, so β=2.5,2.8 sit *above* the transition at these small N_τ: those Δ are finite-T screening/thermal masses, and only β=2.2 is in the confined phase at both L. All remain gapped.

## vs the proxies
The G″≈−0.165 curvature, residue and induced-mass probes (W6) are effective-action **curvatures/signs**, not masses; W6 already notes curvature sign and the spectral gap are independent objects. The direct measurement supplies the actual number those proxies target: a real, positive m_{0⁺⁺}≈1/a. Nothing here contradicts G″<0; it shows the spectrum is gapped precisely where the marginal curvature probe was negative. Data/scripts: `tm_gap_run.py`, `tm_gap_gevp.py`, `tm_gap_results.csv`, `tm_gap_data/`.

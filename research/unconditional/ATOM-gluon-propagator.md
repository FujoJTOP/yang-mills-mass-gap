# ATOM — direct transverse gluon propagator at k→0 (SU(2), L≤8)

**Verdict.** D(k) is finite and positive at every momentum the small lattice resolves, and it falls *steeper than* 1/k² (D·k² decreasing) — unlike a free massless field. A simple massive pole D=Z/(k²+m²) is excluded (**m²<0** at every β,L); the shape is *consistent with* the shoulder of a Gribov–Stingl propagator whose IR plateau would sit *below* k_min (a finite-D(0) GS form fits L=8). The k→0 data do **not** resolve D(0): on L≤8 the transverse mode at k→0 is exactly the unprobed one — the edge/localization obstruction of `GAP-G2` in numbers.

**Setup.** SU(2) Wilson, heat-bath (cross-checked vs Metropolis); β=4/g², β_c≈2.5. Landau gauge by parity Gauss–Seidel (divergence/RMS ≈ 2–3×10⁻³; longitudinal power ≲10⁻⁵·D). Gauge field A_μ^a=2 Im tr(U_μτ^a); *exact* lattice-transverse projection matched to the forward-difference gauge condition (L_μ=1−e^{−ik_μ}): D(k)=⅓⟨Σ_{μ,a}|Ã_μ^a(k)|²⟩/V, lattice momentum k̂_μ=2 sin(πn_μ/L). **Control:** a free massless field reproduces D∝1/k² with D·k² flat (1.27–1.33 for k̂²≥2), validating the estimator; using the exact log field instead of the linear A shifts D by +17% but leaves the D·k² trend unchanged. Errors: jackknife + blocking. The k̂²=1 shell is unusable (axis anisotropy ×2, halves disagree) and excluded.

**Numbers (L=6, β=2.5, 100 cfgs; plaq β(1−⟨P⟩)=0.87 ⇒ ⟨P⟩=0.65):**

| k̂² | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| D(k) | 4.80 | 2.72 | 1.80 | 1.14 | 0.83 |
| D·k² | 9.61 | 8.15 | 7.21 | 6.85 | 6.61 |

Massive fit over k̂²≥2: **m²=−0.53±0.03** (intercept b<0 — D extrapolates *up*, not to a pole). Across β: m²=−1.37 (β=2.2), −0.53 (2.5), −0.34 (2.8); D(2)=7.63, 4.80, 4.02. The Gribov–Stingl form D=Z(k²+m₁²)/(k⁴+m₂²k²+m₃⁴) fits the *shape* but is volume-degenerate: D(0)→∞ (∼10¹⁸) at most points, finite (≈15–100) only where unconstrained — no robust D(0).

**Volume / β (β=2.5, fixed k̂²):** D(k̂²=2)=6.24, 4.80, 3.91 for L=4,6,8; D(k̂²=4)=2.02, 1.80, 1.71. Values shrink with L (finite-volume enhancement decaying) and with β; the L=6→8 change is small at high k̂² (approaching a continuum value). But the smallest resolved k̂² drops only to 0.59 (L=8), where D=24 and still *rising* steeply. D·k² is still increasing as k↓ over the entire reliable window, so the turnover where D(0) saturates lies below k_min=2π/L.

**Bottom line.** Over k̂²∈[0.6,16] the propagator is finite and monotone-decreasing; it is *not* massless-1/k² (D·k² falls 10→6.6, a genuine interaction shoulder), and the simple massive ansatz is excluded. But D(0) is neither measured (k=0 removed by gauge fixing) nor reliably extrapolated (massive→m²<0; GS diverges). Reaching the plateau — seen at 64⁴–128⁴ by Cucchieri–Mendes (PRL 100, 2008) and Bogolubsky et al. (PLB 676, 2009) — needs L≈24–32 at these β. On this evidence alone **m_eff(0)>0 is neither demonstrated nor excluded**: direct small-lattice measurement confirms instead that the k→0 transverse mode is precisely the edge a fixed small volume cannot see. *(Scripts: `research/lattice-cpb/gluon_prop.py`, `run_production.py`; raw output `gluon_prop_results.json`.)*

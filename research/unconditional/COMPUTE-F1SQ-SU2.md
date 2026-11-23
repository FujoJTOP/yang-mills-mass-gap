# COMPUTE-F1SQ-SU2 — direct lattice measurement of the 0⁺⁺ residue `f₁²=|⟨0|tr F²|G⟩|²` on SU(2)

> 2026-11-23 · The Clay-residue computation, run directly. Connected zero-momentum correlator `C(t)=⟨O(t)O(0)⟩_c`, `O_s(t)=Σ_x (1/3)Σ_{i<j∈spatial}(1/2)Re tr U_ij` from APE-smeared spatial links (levels s=0,2,4,8; α=0.5), SU(2) heat bath (CM), L=6, single process, memory-safe streaming (only N×L time-slice sums stored; peak RSS ≈ 60 MB). Two points: β=1.0 and β=2.2. Residue = amplitude `A` of the single-cosh fit `C(t)=A(e^{-mt}+e^{-m(T-t)})` over t=1,2 with `m` fixed to the measured plateau `m_G`, errors = 40-block jackknife. Code: `_f1sq_su2_mem.py` (run), `_f1sq_su2_residue.py` (table; built on `_f1sq_su2_analyze.py` correlator machinery). Data: `f1sq_b1.0_L6.npz`, `f1sq_b2.2_L6.npz` (each N=4000).

**Verdict.** `f₁²>0` is numerically **established at β=2.2** at high significance and **not resolvable at β=1.0** on this lattice. The task's expectation is inverted: the *finer* lattice (β=2.2, smaller `m·a`) resolves the state; the *coarse* strong-coupling lattice (β=1.0, `m·a` at the cutoff) does not. Both statements are physical, not statistical failures of the run.

**Results (L=6, N=4000 configs each; β=1.0: therm 2000, gap 6, seed 11, 665 s; β=2.2: therm 2500, gap 8, seed 22, 724 s).**

| β | `m_G` (plateau) | level s | C(1)/C(0) (S/N) | f₁² = A ± err (m fixed) | S/N |
|---|---|---|---|---|---|
| **2.2** | **1.51 ± 0.15** | 0 (bare tr F²) | +0.123 ± 0.007 (17σ) | **+7.14 ± 0.42** | 17 |
| | | 2 | +0.199 ± 0.010 (24σ) | +5.17 ± 0.26 | 20 |
| | | 4 | +0.197 ± 0.009 (25σ) | +2.72 ± 0.13 | 22 |
| | | 8 | +0.164 ± 0.007 (25σ) | +0.885 ± 0.041 | 22 |
| **1.0** | ≳3–6 (at cutoff) | 0 | +0.0013 ± 0.0068 (**0.2σ**) | not extractable | −0.1 |
| | | 2 | +0.0023 ± 0.0069 (0.3σ) | " | 0.2 |
| | | 4 | +0.0025 ± 0.0078 (0.3σ) | " | 0.3 |
| | | 8 | −0.0012 ± 0.0080 (0.1σ) | " | −0.1 |

`f₁²` is per-operator (the smearing levels are different operators with different normalizations); the bare `tr F²` value `7.14` also needs continuum renormalization, so its absolute value is scheme-dependent — the robust, scheme-independent content is **positivity**. A crude free-mass fit gives `f₁²_bare = 7.1 ± 0.4(stat) −1.3/+1.7(syst from m_G)`, still >4σ from zero.

**β=2.2 (confined point, below β_c≈2.30).** `C(t)/C(0)` is positive and monotone; effective mass forms a plateau `m_eff(t=2) = 1.53, 1.50, 1.49` (levels 2,4,8) → `m_G = 1.51 ± 0.15 a`, consistent with the known SU(2) 0⁺⁺ glueball (`m a ≈ 1.3–1.5` at β≈2.2–2.3). Every smearing level gives a strongly positive residue. **`f₁²>0` is numerically confirmed.**

**β=1.0 (strong coupling).** The connected correlator is a contact term: `C(0)=16` but `C(t>0)<0.35`, with `C(2),C(3)` slightly *negative* — a correlation length below one lattice spacing. `m_eff(t=1)≈6`, i.e. the lightest state sits at/beyond the cutoff (`m·a≳3`). The t=1 signal is `0.2σ`; `f₁²` is not extractable, consistent with zero. This is **not** evidence against P3: P3 is an analytic statement in the β→0 regime where `m·a` diverges, so the state is *provably* present but *numerically sub-lattice-spacing* — unresolvable at any statistics on T=6. More configs cannot help; one needs larger time extent (T≳12–16) at fixed L, not more N. Extrapolating the 0.2σ floor, a 5σ detection at β=1.0 on T=6 would need ~10⁶ configs — and would still sit at the cutoff.

**Is P3 numerically confirmed?** The positivity `f₁²>0` that P3 asserts is **confirmed where the state is resolvable** (β=2.2, >17σ). It is **not** confirmed at β=1.0 purely because strong coupling = coarse lattice = the state is at the cutoff on L=6. Numerical positivity and strong-coupling provability live at opposite ends of the coupling.

**Sharpest line.** The Clay residue is directly seen: at β=2.2, `f₁²_bare = 7.1(4)` with `m_G = 1.51(15)a` (17σ); at β=1.0 it dissolves into the cutoff (0.2σ) — strong coupling is analytically *provable* but numerically the *hardest* point, because there `m·a` exceeds the inverse lattice spacing. `f₁²>0` is a measured fact in the crossover region; its continuum value needs renormalization (`Z_a`, `a→0`), which remains the open wall, not the sign.

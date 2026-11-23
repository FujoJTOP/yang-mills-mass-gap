# E4-RMT-EDGE — the RMT / free-probability edge is the random-Schrödinger DOS edge; its theorems size the Savvidy tachyon and point **anti-gap**, and no moment/edge bound controls `Var(a)` non-circularly

> 2026-11-23 · NEW SPECIES E4. Treat `H_slow=−D²(Ā)+2ad F(Ā)` (Savvidy-indefinite, `λ_min=−g|B|`) as a random operator over the gauge ensemble; use RMT / free-probability edges (Wigner, Tracy–Widom, Marchenko–Pastur, the negative-spectrum density) to bound the negative spectrum's contribution to `Var(a)` / the backreaction `BC⁻¹Bᵀ`. **判词: does an RMT edge/moment bound control the variance non-circularly? — No.** `H_slow` is a *random Schrödinger operator* (deterministic `−D²`, sign-indefinite `V=2ad F_B`); the rigorous edge theory is the Anderson/Lifshitz DOS edge, whose theorems give the spectrum's *density/fluctuation*, never a positive floor. The one positive moment is exact and field-blind (`Tr H=Tr(−D²)≥0`); the one edge is the tachyon `|λ_min|=g|B|`, monotone ↑ in `|B|`; the only sign flip needs conditioning on a positive-Rayleigh-quotient slow mode — an exponentially rare large-deviation tilt whose rate is the effective action, H1. **Five-point: #1 PASS-in-principle(unbuilt); #2,#3,#5 FAIL; #4 PASS narrow, target-vacuous. Acid test passes trivially, but the YM signal is the wrong-sign tachyonic edge, not a mass.**

**§1 The three edges, applied.** *(Wigner)* the semicircle is symmetric about `0` — sign-blind to `G″>0` vs `<0`; and `H_slow` is a differential operator plus `ad F`, not an iid Wigner matrix. *(Marchenko–Pastur)* MP is the PSD sample-covariance edge `[(1−√γ)²σ²,(1+√γ)²σ²]`; the only PSD object is the stiff `C⁻¹` (SZZ-banked, B5 §2), whose floor *closes* as `γ→1` (soft-mode density `D~V`) ⇒ **no uniform floor (#5)**; the dangerous negatives sit in the slow block `A`. *(Tracy–Widom)* the `N^{−2/3}` fluctuation of the extreme eigenvalue, absolute location ensemble-set — a concentration statement, not a deterministic floor. *(Lieb–Thirring, banked)* `Tr(H_-)≤8C_N∫|F|³` is an exact **upper** bound on negative weight.

**§2 The edge is where the ESD premise fails — and it is a tachyon.** `H_slow` indefinite: extreme eigenvalue `λ_min=−g|B|`, count `∝(gB)²V`, depth `∝|B|` (PROOF-SAVVIDY-EXACT-SIGN §1). The one-loop ESD enters as `∫λ^{-1}d⟨a,E_λ a⟩`; at the negative edge this has an imaginary part (decay) — the ESD is the *one-loop* object and its edge is the *breakdown locus* of that expansion. RMT is controlled *away* from the edge; the gap lives *at* it. Reading `|λ_min|=g|B|` returns a tachyonic scale, F-dependent, increasing in the field — the **anti-gap**. The sign-carrier is the *conditioned* quotient `E[H|0]≈2>0` (B4-VARIANCE-MAG; COMPUTE-GPP-SU3), not the unconditioned ESD (90–97% negative in band, W6-BE).

**§3 Conditioning = large-deviation tilt = H1.** Positivity needs conditioning on `E[H|a_0]>0` — an exponentially rare event (rate `∼V`). The conditioned ESD is a large-deviation tilt (Voiculescu deformation / Ben-Arous–Guionnet); its free energy is the extensive effective-action difference, i.e. `W>Var−E[H]`. The transfer RMT cannot make — raw `90%`-negative ESD → `E[H|0]>0` — *is* H1. Circular.

**§4 `BC⁻¹Bᵀ` via MP — wrong sector, wrong structure.** `B` is the cubic slow↔fast coupling: banded (nearest momentum shell), **not** iid; `C⁻¹` is PSD only on the SZZ stiff block, while the negative band is in `A`. MP's hypotheses (iid, PSD) fail on both counts; even granted, MP gives `|BC⁻¹Bᵀ|≲‖B‖²/γ_fast`, and `‖B‖`-small / `γ_fast`-large is verbatim B5's missing estimate `|BC⁻¹Bᵀ|<A+W_geo` = H1.

**Five-point scorecard.**
1. **Nonperturbative — PASS in principle, unbuilt.** Exact `N=∞` (free CLT; Lifshitz/Wegner rigorous), but the 4D ensemble is the unbuilt master field (H1-T8-free-prob §1; free-prob proven only 0D/2D).
2. **Sign-controlling LOWER bound — FAIL.** Moments are sign-blind; `Tr H=Tr(−D²)≥0` is trace-half-positivity (`Tr(H_-)≤Tr(H_+)`, sat. `0.94`), positive weight field-proportional — no `λ>0`; the edge is the tachyon.
3. **Non-circular — FAIL.** Ensemble unbuilt; conditioned ESD = LD tilt, rate = effective action = H1; MP backreaction = B5's missing step.
4. **Non-SUSY + fails free field — PASS (narrow, target-vacuous).** Non-SUSY ✓; free `F=0⇒V=0⇒H=−Δ≥0`, edge `0`, gapless, correctly returns no gap ✓. But the YM signal is the negative tachyon edge — certifies no mass.
5. **Uniform `λ` — FAIL.** Tracy–Widom fluctuation at ensemble-set location; MP floor closes as `D~V`; no `λ` uniform over `[1/48,2.5]`.

**Exact failing estimate.**
```
H_slow = −D² + 2ad F_B   (random Schrödinger: −D² kinetic, V=2ad F_B sign-indefinite)
Uncond. edge: λ_min=−g|B|<0 , N_neg/V∝(gB)² , depth∝|B|        [Savvidy; PR§1]
Moment (n=1): Tr H=Tr(−D²)≥0 ⇒ Tr(H_-)≤Tr(H_+) , sat 0.94 , NO gap .
Edge:         |λ_min|=g|B| ↑ in |B| = tachyon scale = ANTI-gap, not m .
One-loop:     Var = ∫λ^{-1} dESD → Im part at negative edge = decay;
              the edge is the ESD's invalidity locus, not the physical floor.
Condition:    E[H|0]≈2>0  ≠  uncond. ESD (90–97% neg) ⇒ LD tilt, rate = eff. action = H1 .
MP route:     |BC⁻¹Bᵀ|≲‖B‖²/γ_fast < A+W_geo ; B banded≠iid, C⁻¹ PSD only stiff ⇒ = H1 .
Free field:   F=0 ⇒ V=0 ⇒ H=−Δ≥0 , edge 0 , gapless ✓ — YM signal is −g|B|, wrong sign.
```

**Sharpest line.** Correctly named, the RMT/free-probability edge is random-Schrödinger DOS-edge theory — and its theorems do the **opposite** of a gap: a sign-indefinite potential pulls the spectral bottom to `−g|B|` and fills a density of negative states there (`∝∫|F|³`); Lifshitz/Wegner describe the DOS *near that tachyon edge*, bound its regularity, never assert a positive floor. A moment bound returns only trace-half-positivity; an edge bound returns the fluctuation of the tachyon. The edge flips sign only by the exponentially rare conditioning whose rate function *is* the effective action, H1. **E4 does not merely fail to bound the variance non-circularly: it locates the gap at the very edge where its one-loop ESD premise is invalid — the edge is the wall** (k0-SPECTRAL-REDUCTION's "the tail's gap is the gap," now "the edge's tachyon is the anti-gap").

*Refs.* Wigner, *Ann. Math.* **62** (1955) 548; Marchenko–Pastur, *Mat. Sb.* **72** (1967) 507; Tracy–Widom, *CMP* **159** (1994) 151; Erdős–Yau–Yin (arXiv:1212.0833); Lifshitz, *Sov. Phys. JETP* **26** (1968) 462; Wegner, *Z. Phys. B* **44** (1981) 9; Ben-Arous–Guionnet, *CMP* **218** (2001) 177; Lieb–Thirring, *Stud. Appl. Math.* **55** (1976) 179; Voiculescu, *J. Funct. Anal.* **66** (1986) 323; VDN, *Free Random Variables* (AMS 1992); Biane, *Indiana Univ. Math. J.* **46** (1997); Savvidy, PLB **71** (1977) 133; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440. Internal: k0-FREEPROB-AVERAGE, k0-SPECTRAL-REDUCTION, PROOF-SAVVIDY-EXACT-SIGN, B5-CORRELATION-CONTROL, B4-VARIANCE-MAG, C7-SPECIES-CONSTRUCT, H1-T8-free-prob, GAP-G2-edge, PROOF-CONCENTRATION, W6-BE-at-2.5, COMPUTE-GPP-SU3, PROOF-COARSE-POINCARE (H1).

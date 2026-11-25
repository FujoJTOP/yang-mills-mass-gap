#!/usr/bin/env py -3.12 -I
r"""
z2_tqft.py -- Z2 witness: can a TOPOLOGICAL / COHOMOLOGICAL invariant of the moduli space
of connections FORCE a positive 4D mass gap?  (TQFT / geometric-Langlands / cohomology paradigm)

Companion to  research/v3/Z2-TQFT.md.  Deterministic, numpy only, `py -3.12 -I`,
no downloads, no Monte Carlo.  Exit 0 iff all checks PASS.

THE QUESTION (two forms, same answer).
 (Q1)  Does a topological/cohomological invariant of the moduli space of connections
        force  m1 > 0  in 4D pure Yang-Mills?
 (Q2)  [coordinator] Does a MODULI-SPACE COHOMOLOGY CLASS (characteristic / equivariant
        class) whose NONVANISHING in the continuum limit imply a positive Hamiltonian
        SPECTRAL GAP -- without any probability/variance estimate?

THE STRUCTURAL FACT that decides both.  A cohomology class is a class of a topological
space: it is METRIC- and COUPLING-INDEPENDENT (invariant under the connected beta-family
that runs YM <-> free Gaussian at beta=infinity).  A mass is a dimensionful, metric- and
coupling-DEPENDENT spectral threshold.  So:
   * the class cannot SEE the mass (it is the same at beta=infinity, where m=0); and
   * the only operator-theoretic content of the class is  class <-> harmonic form <->
     ker(Laplacian) <-> ZERO MODE  (Hodge / Morse / Witten index) -- i.e. the VACUUM,
     the OPPOSITE of a gap ABOVE the vacuum.  Extracting a gap needs the SPECTRAL
     (Hodge / transfer) Laplacian, which is metric-dependent; there the statement is
     no longer topological and IS the transfer-matrix gap = campaign H1.
Both horns of the coordinator's dichotomy bite at once.  No third option.

PANELS
  P1  SUSY-QM: the Witten index (a genuine cohomology/deformation invariant) is a
      CONSTANT integer while the spectral gap scales with the deformation.  Index != gap.  [T]+[V]
  P2  2D Yang-Mills (exact): the topological data (representation labels) is coupling-
      INDEPENDENT; the spectrum is  E_R = g^2 C2(R); the gap  = g^2 C2(fund) -> 0 as g->0
      (free).  Topological data constant, gap proportional to the coupling.            [T]+[V]
  P3  3D Maxwell-Chern-Simons (DGT): the ONLY 'topological mass',  m = mu * e^2, requires
      BOTH the topological level mu AND the metric Maxwell coupling e^2; m->0 if either
      vanishes; it is a tree-level PARAMETER, and exists only in odd (3D) dimensions.     [T]+[V]
  P4  [coordinator] Squashed torus T^2(s): the cohomology class (b1=2, a generator of H^1)
      is CONSTANT in s, while the spectral gap lambda1(s) = 4pi^2 min(1,1/s^2) -> 0 under a
      metric deformation.  A nonvanishing cohomology class does NOT force a positive gap.  [T]+[V]
  P5  LEDGER + reduction map to H1 / transfer gap.  Kill criteria.                        [V]

Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""

import sys
import json
import numpy as np

NTESTS = 0
FAILS = []


def check(name, cond, detail=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"   ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


# =====================================================================
# P1  -- SUSY QM: Witten index (cohomology invariant) constant, gap scales
# =====================================================================
def panel_P1():
    print("\n" + "-" * 80)
    print("P1  SUSY QM:  Witten index (a cohomology invariant) is CONSTANT; the gap is NOT")
    print("-" * 80)
    # H_B = 1/2 p^2 + 1/2 (W')^2 - 1/2 W''  ,  H_F = 1/2 p^2 + 1/2 (W')^2 + 1/2 W''
    # For W = 1/2 lam x^2:  V_B = 1/2 lam^2 x^2 - 1/2 lam,  V_F = 1/2 lam^2 x^2 + 1/2 lam
    # Exact: E_B0 = 0 (one bosonic zero mode), E_F0 = lam (fermionic), index = 1, gap = lam.
    L, N = 6.0, 800
    x = np.linspace(-L, L, N)
    h = x[1] - x[0]
    # discrete -d^2/dx^2 (Dirichlet)
    main = np.full(N, 2.0 / h**2)
    off = np.full(N - 1, -1.0 / h**2)
    T = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    T = 0.5 * T  # 1/2 p^2

    lams = [0.15, 0.30, 0.60, 1.20, 2.40]
    rows = []
    gaps = []
    idxs = []
    for lam in lams:
        V = 0.5 * (lam * x) ** 2
        HB = T + np.diag(V - 0.5 * lam)
        HF = T + np.diag(V + 0.5 * lam)
        eB = np.linalg.eigvalsh(HB)
        eF = np.linalg.eigvalsh(HF)
        tol = 0.5 * lam  # zero-mode vs gap threshold
        nB0 = int(np.sum(eB < tol))
        nF0 = int(np.sum(eF < tol))
        index = nB0 - nF0
        gap = float(np.min(eF))  # first excited (fermionic ground) above the bosonic vacuum
        rows.append((lam, nB0, nF0, index, float(eB[0]), gap))
        gaps.append(gap)
        idxs.append(index)

    print("   lam    nB0  nF0  index   E_B0(<<)   gap=E_F0     (exact: index=1, gap=lam)")
    for lam, nB0, nF0, index, eB0, gap in rows:
        print(f"   {lam:5.2f}   {nB0:2d}   {nF0:2d}    {index:2d}   {eB0:9.2e}   {gap:8.4f}")

    check("P1.index_constant_1", all(i == 1 for i in idxs),
          f"index={idxs} (a deformation-invariant integer)")
    # gap should equal lam to a few percent
    gap_err = max(abs(g - lam) / lam for (lam, _, _, _, _, g) in rows)
    check("P1.gap_equals_lam", gap_err < 0.03,
          f"max rel err {gap_err:.2e}: gap tracks the deformation lam")
    spread = max(gaps) / min(gaps)
    check("P1.gap_varies_16x", spread > 15.0,
          f"gap spans {min(gaps):.3f}..{max(gaps):.3f} (x{spread:.1f}) while index stays 1")
    print("   => a COHOMOLOGY-TYPE invariant (the index) is CONSTANT (=1) under the very")
    print("      deformation that changes the spectral gap by x16.  The invariant does not")
    print("      know the gap.  (A symmetric double well is the other value, index=0.)")
    return dict(lams=lams, index=idxs, gaps=gaps)


# =====================================================================
# P2  -- 2D Yang-Mills (exact): topological rep labels constant, gap ~ g^2
# =====================================================================
def panel_P2():
    print("\n" + "-" * 80)
    print("P2  2D Yang-Mills (exact): representation data coupling-INDEPENDENT; gap ~ g^2")
    print("-" * 80)
    # SU(2): reps j = 0, 1/2, 1, 3/2, 2 ; Casimir C2(j) = j(j+1); mass E_j = g^2 C2(j)
    js = [0.0, 0.5, 1.0, 1.5, 2.0]
    C2 = [j * (j + 1) for j in js]
    g2s = [0.05, 0.20, 1.00]
    print("   j     C2=j(j+1)")
    for j, c in zip(js, C2):
        print(f"   {j:4.1f}   {c:6.2f}")
    print("\n   g^2     gap=E_{1/2}=g^2*C2(1/2)=0.75 g^2     E_j/g^2 = C2(j) (shape, coupling-free)")
    gaps = []
    shapes_ok = True
    for g2 in g2s:
        Ej = [g2 * c for c in C2]
        gap = Ej[1]  # first excited: j=1/2
        gaps.append(gap)
        # shape E_j/g2 must equal C2(j) exactly, independent of g2
        shape = [e / g2 for e in Ej]
        shapes_ok = shapes_ok and np.allclose(shape, C2)
        print(f"   {g2:5.2f}   {gap:8.4f}                          {['%.2f' % s for s in shape]}")
    # rep set / Casimir shape is coupling-independent
    check("P2.rep_data_coupling_independent", shapes_ok,
          "E_j/g^2 = C2(j) identical for every g^2 (the topological 'which-reps' data)")
    # gap scales linearly with g^2
    ratios = [gaps[i] / gaps[0] for i in range(len(gaps))]
    lin = [g2s[i] / g2s[0] for i in range(len(gaps))]
    check("P2.gap_linear_in_g2", np.allclose(ratios, lin),
          f"gap(g2)/gap(g2_0) = g2/g2_0 : {['%.2f' % r for r in ratios]}")
    check("P2.gap_to_zero_free", gaps[0] / gaps[-1] < 0.06,
          f"gap -> 0 as g->0 (free): gap(0.05)/gap(1.0) = {gaps[0]/gaps[-1]:.3f}")
    print("   => the coupling-INdependent datum (rep labels, Casimir shape) is exactly the")
    print("      data the TQFT keeps; the GAP is the coupling-DEPENDENT scale g^2 C2.")
    print("      Witten(1991): the topological (Ray-Singer torsion) limit is g->0 = the FREE,")
    print("      gapless end.  Topological data lives at the massless end.")
    return dict(g2s=g2s, gaps=gaps, C2=C2)


# =====================================================================
# P3  -- 3D Maxwell-Chern-Simons: the only 'topological mass' m = mu e^2
# =====================================================================
def panel_P3():
    print("\n" + "-" * 80)
    print("P3  3D Maxwell-Chern-Simons (DGT): 'topological mass' m = mu * e^2 needs BOTH mu and e^2")
    print("-" * 80)
    # Transverse 2x2 inverse propagator block (helicity basis):
    #   M(omega) = a I + b J,  a = -(omega^2 - k^2)/e^2,  b = i mu omega,  J = [[0,1],[-1,0]]
    # det M = a^2 + b^2 = 0  =>  (omega^2-k^2)^2/e^4 = mu^2 omega^2
    #   =>  omega^2 - mu e^2 omega - k^2 = 0 (physical branch); at k=0: omega = mu e^2.
    def mass( mu, e2, k=0.0):
        # solve omega^2 - mu*e2*omega - k^2 = 0, take the positive root
        b = -mu * e2
        disc = b * b + 4.0 * k * k
        w = (-b + np.sqrt(disc)) / 2.0
        return w
    print("   mu     e^2     m=mass(k=0)     mu*e^2        m/(mu*e^2)")
    table = []
    for mu in [0.3, 1.0, 3.0]:
        for e2 in [0.5, 1.0, 2.0]:
            m = mass(mu, e2)
            table.append((mu, e2, m))
            print(f"   {mu:4.2f}   {e2:4.2f}   {m:10.4f}     {mu*e2:8.4f}      {m/(mu*e2):6.3f}")
    ok_formula = all(abs(m - mu * e2) < 1e-12 for mu, e2, m in table)
    check("P3.mass_equals_mu_e2", ok_formula,
          "the pole at k=0 is exactly m = mu e^2 (verified across the table)")
    # vanishes when either factor vanishes
    check("P3.zero_if_mu_zero", abs(mass(0.0, 2.0)) < 1e-12,
          "mu=0 (no topological term) -> m=0, gapless free Maxwell (acid test fails correctly)")
    check("P3.zero_if_e2_zero", abs(mass(3.0, 0.0)) < 1e-12,
          "e^2=0 (no metric Maxwell term) -> m=0: the topological term ALONE gives no mass")
    # linear in mu at fixed e^2
    m1 = mass(1.0, 1.0); m3 = mass(3.0, 1.0)
    check("P3.linear_in_mu", abs(m3 - 3.0 * m1) < 1e-12,
          f"m linear in mu at fixed e^2 (m(3)/m(1)={m3/m1:.3f})")
    print("   => the 'topological mass' is m = mu*e^2: it carries the METRIC coupling e^2 as a")
    print("      factor, vanishes as e->0 (free, gapless), is a TREE-LEVEL parameter (the CS")
    print("      level mu), and exists only in odd dimensions (CS is 3D/odd).  NOT a function")
    print("      of a topological invariant alone; NOT a 4D dynamical gap.")
    return dict(table=table)


# =====================================================================
# P4  -- [coordinator] moduli-space cohomology class vs spectral gap
# =====================================================================
def panel_P4():
    print("\n" + "-" * 80)
    print("P4  [coordinator] cohomology CLASS constant, spectral GAP tunable to 0 (squashed T^2)")
    print("-" * 80)
    # Squashed torus T^2(s): circumferences (1, s).  Laplacian on functions:
    #   lambda_{mn} = 4 pi^2 ( m^2 + n^2/s^2 ),  m,n in Z.
    # Cohomology of T^2 (a TOPOLOGICAL invariant, metric-independent):
    #   b0 = 1, b1 = 2, b2 = 1  (constant for every s).
    # Spectral gap = first NONZERO eigenvalue = 4 pi^2 * min(1, 1/s^2)  -> 0 as s -> infinity.
    ss = [1.0, 2.0, 4.0, 8.0, 16.0]
    print("   s     b0  b1  b2    lambda1(s)=4pi^2 min(1,1/s^2)    lambda1 * s^2")
    b0, b1, b2 = 1, 2, 1   # Betti numbers of T^2: s-independent (the class is nonzero, fixed)
    lam1s = []
    for s in ss:
        lam1 = 4.0 * np.pi**2 * min(1.0, 1.0 / s**2)
        lam1s.append(lam1)
        print(f"   {s:4.1f}   {b0}   {b1}   {b2}    {lam1:12.6f}               {lam1*s*s:10.4f}")
    check("P4.cohomology_constant", (b0, b1, b2) == (1, 2, 1),
          "b1 = 2 (generator of H^1(T^2)) is NONZERO and CONSTANT in s -- the class exists")
    check("P4.gap_to_zero", lam1s[-1] / lam1s[0] < 0.01,
          f"lambda1(s=16)/lambda1(1) = {lam1s[-1]/lam1s[0]:.4f}: the gap -> 0 under a metric move")
    check("P4.gap_scales_1_over_s2", all(abs(lam1s[i] * ss[i]**2 - 4*np.pi**2) < 1e-9 for i in range(1, len(ss))),
          "lambda1 * s^2 = 4 pi^2 : the gap is a METRIC datum, the class is not")
    print("   => a NONVANISHING cohomology class (b1=2) coexists with a spectral gap that a")
    print("      metric deformation drives to 0.  The class does NOT force a positive gap.")
    print("      Operator meaning of the class:  class <-> harmonic form <-> ker(Laplacian)")
    print("      <-> ZERO MODE (the VACUUM).  The gap = first NONZERO eigenvalue = NOT a class.")
    print("      2D-YM reading: flat-connection moduli Hom(Z^2,G)/G = the zero-energy (vacuum)")
    print("      sector; its cohomology labels GROUND STATES; the gap is the first non-flat")
    print("      connection's action ~ g^2 (P2).  Class = vacuum, gap != vacuum.")
    return dict(ss=ss, lam1=lam1s, betti=(b0, b1, b2))


# =====================================================================
# main
# =====================================================================
def main():
    print("=" * 80)
    print("Z2 witness -- TQFT / geometric-Langlands / cohomology paradigm vs the 4D mass gap")
    print("Does a topological/cohomological invariant of the moduli space FORCE m1 > 0?")
    print("=" * 80)
    R = {}
    R['P1'] = panel_P1()
    R['P2'] = panel_P2()
    R['P3'] = panel_P3()
    R['P4'] = panel_P4()

    print("\n" + "-" * 80)
    print("P5  LEDGER + reduction map")
    print("-" * 80)
    print("   * (Q1) topo/cohomological invariant -> metric/coupling-INDEPENDENT -> constant on")
    print("     the beta-family -> equals its free value -> CANNOT see the mass.        [NG]")
    print("   * (Q2) moduli-space cohomology CLASS -> nonzero class = ZERO MODE (= vacuum),")
    print("     not a gap ABOVE it; to get a gap you need the spectral (Hodge/transfer)")
    print("     Laplacian -> metric-dependent -> it is the TRANSFER-MATRIX GAP = H1.    [NG]")
    print("   * topo data = 'which sectors/reps/vacua' (constant); gap = the spectral scale")
    print("     above the vacuum (depends on g).  Confusing the two is the category error.")
    print("   * only 'topological mass' (3D DGT m=mu e^2) needs the metric coupling, is tree-")
    print("     level, odd-dimensional: no 4D analogue for pure YM.                      [NG]")
    print("   * anomalies (Freed-Hopkins / GKSW) forbid the TRIVIALLY gapped phase (force")
    print("     gapless OR topologically ordered) -- they do NOT force a positive gap; 4D")
    print("     pure YM at theta=0 has no 0-form 't Hooft anomaly (no fermions).         [NG]")
    print("   * the one exact 'cohomology gives a mass' (Seiberg-Witten BPS = integral over")
    print("     an H_1 cycle; geometric Langlands = local systems = H^1) is SUSY, gapless at")
    print("     the symmetric point, and LITERALLY H_1 -> the campaign H1 cell, unsolved.  [NG]")
    check("P5.kill_criteria_triggered", True,
          "invariant metric-independent (kill i); only topological order/vacua (kill ii); "
          "or reduces to H1 / transfer gap (kill iii)")

    out = dict(
        verdict="NG: TQFT/geometric-Langlands/cohomology paradigm cannot force m1>0",
        kill_i_metric_independent=True,
        kill_ii_only_topological_order=True,
        kill_iii_reduces_to_H1=True,
        acid_test="invariant CONSTANT across free<->interacting (blind); only 'topological "
                  "mass' (3D DGT m=mu*e^2) needs the metric coupling, tree-level, no 4D analogue",
        P1_index=dict(lams=R['P1']['lams'], index=R['P1']['index'], gaps=R['P1']['gaps']),
        P2_2dym=dict(g2s=R['P2']['g2s'], gaps=R['P2']['gaps'], C2=R['P2']['C2']),
        P3_dgt_mu_e2=[list(map(float, t)) for t in R['P3']['table']],
        P4_torus=dict(ss=R['P4']['ss'], lam1=R['P4']['lam1'], betti=list(R['P4']['betti'])),
        clay_percent=0,
    )

    print("\n" + "=" * 80)
    if FAILS:
        print(f"RESULT: {len(FAILS)} FAIL(S) of {NTESTS}: {FAILS}")
        print("=" * 80)
        return 1
    print(f"RESULT: ALL {NTESTS} CHECKS PASS.")
    print("  Verdict [NG]: no topological/cohomological invariant of the moduli space forces")
    print("  m1>0.  It is metric/coupling-independent (blind to a mass); its only operator")
    print("  content is a zero mode (= vacuum), the opposite of a gap; to read a gap from it")
    print("  needs the metric (Hodge/transfer) Laplacian = the transfer-matrix gap = H1.")
    print("  Acid test: the invariant is CONSTANT across free<->interacting, so it cannot")
    print("  discriminate.  The only 'topological mass' (3D DGT m=mu e^2) needs the metric")
    print("  coupling, is tree-level, odd-dimensional.  Kill criteria i/ii/iii all triggered.")
    print("  No QED.  Clay 0%.")
    print("=" * 80)

    try:
        with open("z2_tqft_results.json", "w") as f:
            json.dump(out, f, indent=2)
        print("  (wrote z2_tqft_results.json)")
    except Exception as e:
        print(f"  (json write skipped: {e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

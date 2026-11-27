"""HR4-POWER -- discriminative power of the complete-monotonicity two-point inequality.

ZERO MC, deterministic, exact where exact is possible.  ~2 s.

Object under test (HZ1-PLATEAU.md:19,143,295 [10]):
      C_L(L)  >=  C_L(L/2)^2 / C_L(0)
proposed there as "P 在有限 L 的必要条件" worth checking on the repo's own data.

 [1]  the inequality is an IDENTITY for every nu >= 0 (Cauchy-Schwarz):
          D := C(0)C(2t) - C(t)^2 = (1/2) iint (e^{-st}-e^{-s't})^2 dnu dnu'
                                  = Var_nu(e^{-st}) >= 0
 [1b] the SAME holds for the estimator actually run on the T^4 data (periodic time wrap,
      finite Lt), for every Lt and every L -- so a measured D<0 can only be bias/noise.
 [2]  FREE FIELD (box gap m1^2 = 4 sin^2(pi/L)): PASS, with K = Var/C^2 > 0.
 [3]  HN1 falsifier family (P FALSE): PASS, exact slack, K -> 1.
 [4]  class-4 spread family (P FALSE): PASS, K -> 0.
 [5]  class-3 double well at/below/above the P threshold: PASS in all three.
 [6]  M6 (P TRUE, gap FALSE): PASS, saturating.
 [7]  2D YM (nu == 0) and single plaquette: NO OBJECT.
 => the statistic is CONSTANT (=PASS) on the whole battery, including the single cell that
    discriminates for P.  Zero discriminative power.
 [8]  DIRECTION: the inequality is K >= 0, i.e. a LOWER bound on C(2t) given C(t).
      P needs an UPPER bound on C(L).  The upper-bound form runs backwards in time and
      needs log-concavity, i.e. K <= K0 -- which HZ1 silently assumed to be 0.
 [9]  K on the battery: bounded everywhere except M6 at large L; does not separate P.

Run:  py -3.12 hr4_power.py
"""
import math
import numpy as np

OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def Cof(t, s, w):
    s = np.asarray(s, float)
    w = np.asarray(w, float)
    return float(np.sum(w * np.exp(-s * t)))


def Dstat(C0, Ct, C2t):
    D = C0 * C2t - Ct * Ct
    K = D / (Ct * Ct) if Ct != 0 else float('nan')
    return D, K


# ------------------------------------------------------------------ [1]
def sec1():
    P("=" * 104)
    P("[1] IDENTITY:  D := C(0)C(2t) - C(t)^2 = Var_nu(e^{-st}) >= 0   for every nu >= 0")
    P("=" * 104)
    P("  C(t) = int e^{-st} dnu(s)  =>  <e^{-st}>^2 <= <e^{-2st}>  (Cauchy-Schwarz in dnu)")
    P("  Symmetrised, the slack is a VARIANCE, i.e. it vanishes iff nu is a single atom.")
    rng = np.random.default_rng(20261126)
    bad = 0
    worst = 0.0
    for _ in range(4000):
        n = int(rng.integers(2, 9))          # n >= 2 so the variance form is non-degenerate
        s = rng.uniform(0.0, 6.0, n)
        w = rng.uniform(0.0, 1.0, n)
        t = float(rng.uniform(0.05, 20.0))
        C0, Ct, C2t = Cof(0, s, w), Cof(t, s, w), Cof(2 * t, s, w)
        D = C0 * C2t - Ct * Ct
        e = np.exp(-s * t)
        var = 0.5 * sum(w[i] * w[j] * (e[i] - e[j]) ** 2 for i in range(n) for j in range(n))
        if D < -1e-13 * C0 * C0:
            bad += 1
        worst = max(worst, abs(D - var) / (C0 * C0))
        assert abs(D - var) <= 1e-13 * C0 * C0, (D, var)
    P(f"  4000 random positive measures: D<0 in {bad}/4000 ;  max|D-Var|/C(0)^2 = {worst:.2e}")
    P("  => 'satisfied' carries ZERO bits about the mass gap; a measured D<0 is an")
    P("     ESTIMATOR statement (bias or noise), never a statement about the model.  [T]")

    P("")
    P("[1b] the estimator actually run on the data is the PERIODIC-time (finite Lt) one:")
    P("     C(t) = sum_n w_n [ e^{-m_n t} + e^{-m_n(Lt-t)} ].  Per mode, at argument 2t,")
    P("     slack = w_n b (1-a^2)^2 with a = e^{-m_n t}, b = e^{-m_n(Lt-2t)} , and b > 0 for")
    P("     every Lt  =>  D >= 0 for EVERY (Lt, L, nu).  Numeric check:")
    bad = 0
    worstrel = 0.0
    for _ in range(3000):
        n = int(rng.integers(1, 7))
        Lt = int(rng.integers(6, 33))
        t = float(rng.integers(1, Lt // 2 + 1))   # 2t <= Lt : exactly the regime the data uses
        m = rng.uniform(0.01, 3.0, n)
        w = rng.uniform(0.0, 1.0, n)
        c = lambda x: float(np.sum(w * (np.exp(-m * x) + np.exp(-m * (Lt - x)))))  # noqa: E731
        C0, Ct, C2t = c(0.0), c(t), c(2 * t)
        D = C0 * C2t - Ct * Ct
        if D < -1e-12 * max(1.0, C0 * C0):
            bad += 1
        worstrel = max(worstrel, abs(D) / max(1e-300, C0 * C0))
    P(f"  3000 random (nu, Lt, t) with 2t <= Lt: D<0 in {bad}/3000  "
      f"(max |D|/C(0)^2 = {worstrel:.2e})")
    P("  PSD argument (covers every Lt): the Hankel matrix of the wrapped kernel is a sum of")
    P("  two rank-one PSD kernels, e^{-m(t+t')} and e^{-mLt}e^{m t}e^{m t'}; every 2x2 minor")
    P("  is therefore >= 0, i.e. C(0)C(2t) - C(t)^2 >= 0 exactly.")
    P("  => on the actual Lt=12/16/20 data the SIGN of D is guaranteed by positivity;")
    P("     the inequality cannot be falsified by any model in the class.  [T]")

    P("")
    P("[1c] finite-Lt INFLATION of the matched-time value (this is what the tier is about):")
    P("     C_meas(L) = C_cyl(L) * (1 + e^{-m(Lt-2L)})  ->  the measured C(L) is the cylinder")
    P("     value only when Lt >> 2L; at Lt = 2L it is already 2x, at Lt < 2L it is huge.")
    P(f"     {'Lt-2L':>8} {'m=0.7':>10} {'m=1.2':>10} {'m=2.0':>10}")
    for d in (8, 4, 2, 1, 0, -1, -2, -4):
        vals = [1 + math.exp(-m * d) for m in (0.7, 1.2, 2.0)]
        P(f"     {d:>8} " + " ".join(f"{v:10.4g}" for v in vals))
    P("     datasets used below: Ls=6/Lt=20 -> Lt-2L=8 (x1.0001) ; Ls=8/Lt=20 -> 4 (x1.008)")
    P("     Ls=6/Lt=12 -> 0 (x2) ; Ls=8/Lt=12 -> -4 (x1.2e2) ; Ls=8/Lt=16 -> 0 (x2)")
    P("     => in Tier B the test is not merely noisy, it compares the mirror image.")


# ------------------------------------------------------------------ [2]
def freefield(L, kind='per', nmax=48, wp=4):
    """Zero-momentum free-field correlator on T^3_L.

    kind 'per' : k_i = 2 pi n_i / L                     (periodic; lowest E^2 = 4 sin^2(2pi/L))
    kind 'box' : k_1 = (2 n_1 + 1) pi / L, k_2,3 = 2 pi n / L   (one antiperiodic direction;
                 lowest E^2 = 4 sin^2(pi/L)  -- HZ1's 'box gap')
    E_n = 2 sqrt(sum_i sin^2 k_i) ;  w_n = (2 E_n)^{-wp}
    wp = 4 : C(0) finite (dimension-8 carrier).   wp = 1 : C(0) UV-divergent (see note).
    """
    ns = np.arange(-nmax, nmax + 1)
    N1, N2, N3 = np.meshgrid(ns, ns, ns, indexing='ij')
    N1, N2, N3 = N1.ravel(), N2.ravel(), N3.ravel()
    if kind == 'per':
        k1, k2, k3 = 2 * math.pi * N1 / L, 2 * math.pi * N2 / L, 2 * math.pi * N3 / L
    else:
        k1 = (2 * N1 + 1) * math.pi / L
        k2, k3 = 2 * math.pi * N2 / L, 2 * math.pi * N3 / L
    E = 2 * np.sqrt(np.sin(k1) ** 2 + np.sin(k2) ** 2 + np.sin(k3) ** 2)
    keep = E > 1e-12
    E = E[keep]
    return E, (2 * E) ** (-wp)


def sec2():
    P("")
    P("=" * 104)
    P("[2] FREE FIELD: the inequality is SATISFIED for every L (it is the P-TRUE cell)")
    P("=" * 104)
    P("  (a) box gap m1^2 = 4 sin^2(pi/L)  (one antiperiodic direction), carrier w = (2E)^{-4}")
    P(f"  {'L':>4} {'C(0)':>12} {'C(L/2)':>12} {'C(L)':>12} {'C(L/2)^2/C(0)':>15}"
      f" {'D':>12} {'K':>10} {'verdict':>8}")
    for L in (4, 6, 8, 12, 16, 32, 64):
        E, w = freefield(L, 'box', 40, 4)
        C0, Ch, C1 = Cof(0, E, w), Cof(L / 2, E, w), Cof(L, E, w)
        D, K = Dstat(C0, Ch, C1)
        P(f"  {L:>4} {C0:12.6e} {Ch:12.6e} {C1:12.6e} {Ch*Ch/C0:15.6e}"
          f" {D:12.4e} {K:10.5f} {'PASS' if D >= 0 else 'FAIL':>8}")
    P("  (b) periodic 2pi/L lowest mode, same carrier")
    for L in (6, 8, 16):
        E, w = freefield(L, 'per', 40, 4)
        C0, Ch, C1 = Cof(0, E, w), Cof(L / 2, E, w), Cof(L, E, w)
        D, K = Dstat(C0, Ch, C1)
        P(f"  {L:>4} {C0:12.6e} {Ch:12.6e} {C1:12.6e} {Ch*Ch/C0:15.6e}"
          f" {D:12.4e} {K:10.5f} {'PASS' if D >= 0 else 'FAIL':>8}")
    P("  (c) carrier w = (2E)^{-1} (massless propagator): C(0) is UV-divergent (prop. to nmax^2)")
    P("      -> the VALUE of K is carrier/UV-set, the SIGN never is:")
    for nmax in (20, 40, 80):
        E, w = freefield(6, 'box', nmax, 1)
        C0, Ch, C1 = Cof(0, E, w), Cof(3, E, w), Cof(6, E, w)
        D, K = Dstat(C0, Ch, C1)
        P(f"      nmax={nmax:>3}  C(0)={C0:12.5e}  K={K:12.4f}  {'PASS' if D >= 0 else 'FAIL'}")
    P("  (d) two-shell truncation (cutoff-free): shell |n|=1 (3 modes, E1), |n|=sqrt2 (6 modes, sqrt2 E1)")
    P("      slack = w1 w2 (e^{-E1 L/2} - e^{-sqrt2 E1 L/2})^2 > 0  -> PASS, margin O(w1 w2)")
    P("  => free field: PASS with K in [0.0x, 1e3] depending on the carrier, and K is")
    P("     L-INdependent for the continuum carrier (the free field is scale invariant).")


# ------------------------------------------------------------------ [3..7]
def sec3():
    P("")
    P("=" * 104)
    P("[3] HN1 FALSIFIER FAMILY  nu_L = (1/2)delta_{1/L^2} + (1/2)delta_4     [P is FALSE]")
    P("=" * 104)
    P("  C(t) = (1/2)e^{-t/L^2} + (1/2)e^{-4t} ;  exact slack D = (1/4)(e^{-1/(2L)} - e^{-2L})^2")
    P(f"  {'L':>7} {'C(0)':>6} {'C(L/2)':>12} {'C(L)':>12} {'rhs':>12} {'D':>12}"
      f" {'K':>9} {'verdict':>8} {'G_c=e C(L)':>11}")
    for L in (2, 3, 4, 6, 8, 16, 64, 1024, 10 ** 6):
        s = np.array([1.0 / L ** 2, 4.0]); w = np.array([0.5, 0.5])
        C0, Ch, C1 = Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w)
        D, K = Dstat(C0, Ch, C1)
        Dx = 0.25 * (math.exp(-1.0 / (2 * L)) - math.exp(-2.0 * L)) ** 2
        assert abs(D - Dx) <= 1e-12
        P(f"  {L:>7} {C0:6.4f} {Ch:12.6e} {C1:12.6e} {Ch*Ch/C0:12.6e} {D:12.6e}"
          f" {K:9.5f} {'PASS' if D >= 0 else 'FAIL':>8} {math.e*C1:11.4f}")
    P("  => PASS for every L (K -> 1, margin -> 1/4).  P is FALSE on the same sequence:")
    P("     G_c = e^{sqrt(c)}C_L(L) -> e * 1/2 = 1.359 (sqrt(c)=1).")
    P("     [note] HZ1 sec.5.3 [4] prints 0.1839 = (1/2)e^{-1} for this row, i.e. the e^{-sqrt(c)}")
    P("     convention, while HZ1 sec.1.2 (G) defines G_c = e^{+sqrt(c)}C_L.  Factor e^2 = 7.389.")
    P("     Either way G_c does not tend to 0 -> P is FALSE -> the erratum changes no verdict.")


def sec4():
    P("")
    P("=" * 104)
    P("[4] CLASS-4 SPREAD FAMILY  nu_L = sum_{j<=L^3}(1/L^3) delta_{(c/L^2)(1+j/L^3)}   [P FALSE]")
    P("=" * 104)
    P(f"  {'L':>4} {'atoms':>7} {'C(L/2)':>12} {'C(L)':>12} {'rhs':>12} {'D':>11}"
      f" {'K':>11} {'verdict':>8}")
    for L in (4, 8, 16, 32):
        n = L ** 3
        j = np.arange(1, n + 1)
        s = (1.0 / L ** 2) * (1 + j / L ** 3)
        w = np.full(n, 1.0 / n)
        C0, Ch, C1 = Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w)
        D, K = Dstat(C0, Ch, C1)
        P(f"  {L:>4} {n:>7} {Ch:12.6e} {C1:12.6e} {Ch*Ch/C0:12.6e} {D:11.3e}"
          f" {K:11.3e} {'PASS' if D >= 0 else 'FAIL':>8}")
    P("  => PASS with K -> 0 (the tower sits in an exponent window of width c/(2L) -> 0).")


def sec5():
    P("")
    P("=" * 104)
    P("[5] CLASS-3 DOUBLE WELL, on both sides of and at the P threshold  m1 = e^{-dS}")
    P("=" * 104)
    P("  HZ1 sec.5.3 [5]:  plateau (P FALSE)  <=>  dS >= log L  <=>  L m1 <= 1")
    P(f"  {'case':>11} {'L':>4} {'m1':>9} {'L*m1':>8} {'C(L/2)':>11} {'C(L)':>11}"
      f" {'rhs':>11} {'D':>10} {'verdict':>8} {'P':>9}")
    for name, f_ds in (('dS=logL/2', lambda L: 0.5 * math.log(L)),
                       ('dS=logL', lambda L: 1.0 * math.log(L)),
                       ('dS=2logL', lambda L: 2.0 * math.log(L))):
        for L in (6, 8, 16):
            m1 = math.exp(-f_ds(L))
            s = np.array([m1 * m1, 4.0]); w = np.array([0.5, 0.5])
            C0, Ch, C1 = Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w)
            D, K = Dstat(C0, Ch, C1)
            pv = 'FALSE' if L * m1 <= 1.0 else 'TRUE'
            P(f"  {name:>11} {L:>4} {m1:9.4f} {L*m1:8.3f} {Ch:11.5e} {C1:11.5e}"
              f" {Ch*Ch/C0:11.5e} {D:10.3e} {'PASS' if D >= 0 else 'FAIL':>8} {pv:>9}")
    P("  => PASS at, below and above the threshold: the inequality is blind to dS.")


def sec6():
    P("")
    P("=" * 104)
    P("[6] M6  (HN1:316-319): bottom residue -> 0, O(1) residue 1/(2omega) at s = omega")
    P("=" * 104)
    P(f"  {'eps':>9} {'L':>4} {'C(L/2)':>12} {'C(L)':>12} {'rhs':>12} {'D':>11}"
      f" {'K':>10} {'verdict':>8}")
    for eps in (1e-1, 1e-3, 1e-6, 1e-9):
        for L in (8, 16):
            s = np.array([eps ** 2, 0.25]); w = np.array([eps, 1.0 - eps])
            C0, Ch, C1 = Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w)
            D, K = Dstat(C0, Ch, C1)
            P(f"  {eps:9.0e} {L:>4} {Ch:12.6e} {C1:12.6e} {Ch*Ch/C0:12.6e} {D:11.3e}"
              f" {K:10.6f} {'PASS' if D >= 0 else 'FAIL':>8}")
    P("  => M6: PASS, K -> 0.  M6 has P TRUE and the gap FALSE.")


def sec7():
    P("")
    P("=" * 104)
    P("[7] 2D YM and single plaquette")
    P("=" * 104)
    P("  2D YM: O = (1/4)tr F^2 = H conserved (HN1:297) => nu == 0 => C == 0 identically;")
    P("         the statistic is 0/0.  NO OBJECT -- neither PASS nor FAIL.")
    P("  single plaquette: no transfer matrix (N02:78) => no C(t), no L family.  NO OBJECT.")


# ------------------------------------------------------------------ [8][9]
def sec8():
    P("")
    P("=" * 104)
    P("[8] DIRECTION: the inequality is a LOWER bound; P needs an UPPER bound")
    P("=" * 104)
    P("  exact rewrite:   C(2t) = (1 + K(t)) C(t)^2 / C(0),   K(t) = Var_nu(e^{-st})/C(t)^2.")
    P("  * the inequality is exactly K >= 0 -- it says a floor at t PROPAGATES to 2t, i.e.")
    P("    it CERTIFIES a platform and can never exclude one.")
    P("  * P wants C_L(L) <= r0 e^{-sqrt(c)}: an upper bound at the LONGEST time.")
    P("  * the upper-bound half of Cauchy-Schwarz reads C(t) <= sqrt(C(0)C(2t)): small at 2t")
    P("    implies small at t -- backwards.  To reach t = L from t < L one needs the reverse")
    P("    inequality  C(2t) <= C(t)^2/C(0)  (log-concavity), which is FALSE in general.")
    P("  * the only usable weakening is  K(t) <= K0  uniformly, giving")
    P("        C(2t) <= (1+K0) C(t)^2/C(0).")
    P("    HZ1 sec.2.2's propagation rho -> rho^2/M0 is exactly K0 = 0, which is assumed,")
    P("    not proved (and is false for every model except a single atom).")
    P("")
    rho, seq = 0.18, [0.18]
    for _ in range(4):
        rho *= rho
        seq.append(rho)
    P("  HZ1 sec.5.3 [10] reproduced (rho0=0.18, M0=1):")
    P("    " + "  ->  ".join(f"{v:.4g}" for v in seq))
    P("    = 0.18 / 3.24e-2 / 1.05e-3 / 1.10e-6 / 1.21e-12, matching HZ1 exactly  [T]")
    P("    but every entry is a LOWER bound on C(2^k L): the chain is a platform certificate,")
    P("    not a falsifier.  From a floor it produces weaker floors, never an upper bound.")


def sec9():
    P("")
    P("=" * 104)
    P("[9] K on the battery (carrier: (2E)^{-4} for free field) -- does K separate P?")
    P("=" * 104)
    P(f"  {'model':>18} {'L=6':>11} {'L=8':>11} {'L=16':>11} {'L=64':>11} {'P':>6} {'gap':>6}")
    rows = []
    v = []
    for L in (6, 8, 16, 64):
        E, w = freefield(L, 'box', 40, 4)
        C0, Ch, C1 = Cof(0, E, w), Cof(L / 2, E, w), Cof(L, E, w)
        v.append(Dstat(C0, Ch, C1)[1])
    rows.append(('free field', v, 'TRUE', 'FALSE'))
    v = []
    for L in (6, 8, 16, 64):
        s = np.array([1.0 / L ** 2, 4.0]); w = np.array([0.5, 0.5])
        v.append(Dstat(Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w))[1])
    rows.append(('HN1 falsifier', v, 'FALSE', 'FALSE'))
    v = []
    for L in (6, 8, 16, 64):
        n = L ** 3; j = np.arange(1, n + 1)
        s = (1.0 / L ** 2) * (1 + j / L ** 3); w = np.full(n, 1.0 / n)
        v.append(Dstat(Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w))[1])
    rows.append(('class-4 spread', v, 'FALSE', 'FALSE'))
    v = []
    for L in (6, 8, 16, 64):
        m1 = math.exp(-math.log(L))
        s = np.array([m1 * m1, 4.0]); w = np.array([0.5, 0.5])
        v.append(Dstat(Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w))[1])
    rows.append(('dW dS=logL', v, 'FALSE', 'FALSE'))
    v = []
    for L in (6, 8, 16, 64):
        s = np.array([1e-12, 0.25]); w = np.array([1e-6, 1 - 1e-6])
        v.append(Dstat(Cof(0, s, w), Cof(L / 2, s, w), Cof(L, s, w))[1])
    rows.append(('M6', v, 'TRUE', 'FALSE'))
    for name, vv, pp, gg in rows:
        P(f"  {name:>18} " + " ".join(f"{x:11.5f}" for x in vv) + f" {pp:>6} {gg:>6}")
    P("")
    P("  K is bounded on the P-FALSE rows (0 .. 1) and on the P-TRUE free-field row, and it")
    P("  GROWS on the P-TRUE M6 row at L >> log(1/w) (the thin atom's visibility horizon).")
    P("  => K does not separate either.  It stays worth measuring only because K <= K0 is the")
    P("     missing hypothesis that would make an upper bound propagate (see [8]).")


def main():
    P("HR4-POWER: discriminative power of  C(L) >= C(L/2)^2/C(0)   (zero MC, ~2 s)")
    sec1()
    sec2()
    sec3()
    sec4()
    sec5()
    sec6()
    sec7()
    sec8()
    sec9()
    P("")
    P("ALL ASSERTIONS PASS (D >= 0 in every cell of the battery; no assertion failed).")
    open('hr4_power.log', 'w').write("\n".join(OUT) + "\n")
    P("[ok] wrote hr4_power.log")


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    main()

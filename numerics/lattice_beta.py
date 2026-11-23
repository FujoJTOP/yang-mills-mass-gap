"""
lattice_beta.py  --  Task 4.

The lattice beta axis  beta = 2N/g^2  (= 6/g^2 for SU(3) Wilson),
2-loop asymptotic scaling, and the two candidate handoff crossings.

--------------------------------------------------------------------------
LITERATURE INPUTS (NOT derived here; flagged explicitly)
--------------------------------------------------------------------------
  * Lambda_MSbar / Lambda_lat = 28.8     (SU(3) Wilson action; standard
    2-loop / nonperturbative estimate.  Quoted range in the literature is
    roughly 25-30; we use 28.8 as instructed.)
  * Delta / Lambda_MSbar  in  [4, 8]     (Delta = the physical handoff /
    mass-gap scale.  This number is itself the NONPERTURBATIVE output one
    would like to prove -- we use it only as an illustrative INPUT range to
    see whether the two candidate windows are close.  It is CIRCULAR as a
    proof, but legitimate as a consistency / plausibility check.)

--------------------------------------------------------------------------
2-loop asymptotic scaling
--------------------------------------------------------------------------
Convention (same as rg_handoff.py):  beta = -b0 g^3/(16pi^2) - b1 g^5/(16pi^2)^2,
b0 = 11N/3 = 11, b1 = 34N^2/3 = 102 for N=3.  Solving the 2-loop flow at
mu = 1/a gives the standard asymptotic-scaling relation

    a^{-1} / Lambda_lat = exp(8 pi^2 / (b0 g^2)) * (g^2 / (16 pi^2))^{b1/(2 b0^2)}
                        =: F_asym(g^2).

(The scheme constant is absorbed into Lambda_lat; using the exact 2-loop
antiderivative exp(A_2loop) instead changes F_asym by ~2% at g^2=1, which
shifts beta by <0.02 -- far below the O(1) effects we are tracking.)

--------------------------------------------------------------------------
Quantities
--------------------------------------------------------------------------
  a*Delta       = Delta / a^{-1}
                = (Delta/Lambda_MS)*(Lambda_MS/Lambda_lat)/F_asym
                = 28.8 * (Delta/Lambda_MS) / F_asym
  xi/a          = 1/(a*Delta).
  q             = a/xi_phys = a*Delta.   So "q ~ 1"  <==>  "a*Delta ~ 1"
                <==>  "xi/a ~ 1"  (the strong-coupling boundary).

Two loci are compared on the SAME beta axis:
  (W)  g^2(a^{-1}) = 1     ->  beta = 6  (trivial; the formal weak threshold);
  (S)  a*Delta      = 1     ->  beta determined by F_asym = 28.8*(Delta/Lambda).
Also reported: the MS-bar running coupling g^2(mu=Delta), to show that the
*physical* handoff scale (mu ~ Delta) has g^2 ~ 1.8-2.2, NOT 1 -- i.e. the
'coupling=1' threshold is not an invariant.

--------------------------------------------------------------------------
Honest caveat
--------------------------------------------------------------------------
At beta ~ 6 the Wilson lattice spacing is a ~ 0.1-0.15 fm (UV), where
2-loop asymptotic scaling has O(30-50%) violations.  So the coincidence of
the two beta windows is an INDICATION of an overlap, not a controlled one.
The MS-bar g^2(Delta) number is cleaner (continuum perturbation theory) but
still uses the literature Delta/Lambda input.
"""
import math

import rg_handoff as R

MPI2 = R.MPI2


# ---------------- literature inputs (flagged) ----------------
LAM_MS_OVER_LAT = 28.8          # SU(3) Wilson, literature input
DELTA_LO = 4.0                  # Delta/Lambda_MS, literature input range
DELTA_HI = 8.0
LAM_MS_GEV = 0.22               # illustrative anchor for a GeV scale only


def coeffs(N=3):
    return R.coeffs(N)


def F_asym(g2, N=3):
    """a^{-1}/Lambda_lat  = 2-loop asymptotic scaling factor."""
    b0, b1, _ = coeffs(N)
    k = b1 / (2.0 * b0 * b0)
    return math.exp(8.0 * math.pi ** 2 / (b0 * g2)) * (g2 / MPI2) ** k


def invA(t, N=3, loops=3, lo=1e-12, hi=0.5):
    """a = g^2/16pi^2 at ln(mu/Lambda)=t, on the weak branch."""
    f = (lambda a: R.A(a, N, loops))
    if f(lo) < t or f(hi) > t:
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > t:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def g2_at_t(t, N=3, loops=3):
    a = invA(t, N, loops)
    return MPI2 * a if a == a else float("nan")


def beta_for_aDelta_one(D, N=3, lo=0.2, hi=3.0):
    """Bare g^2 where a*Delta=1, i.e. F_asym(g^2)=28.8*D.  Return beta=2N/g^2."""
    target = LAM_MS_OVER_LAT * D
    # F_asym is decreasing in g^2
    glo, ghi = lo, hi
    if not (F_asym(glo, N) > target > F_asym(ghi, N)):
        return float("nan"), float("nan")
    for _ in range(200):
        mid = 0.5 * (glo + ghi)
        if F_asym(mid, N) > target:
            glo = mid
        else:
            ghi = mid
    g2 = 0.5 * (glo + ghi)
    return 2.0 * N / g2, g2


def run(outdir="."):
    import os
    import _plot as P

    N = 3
    b0, b1, b2 = coeffs(N)
    res = {}

    # ---------------- beta scan table ----------------
    rows = []
    betas = [4.0 + 0.25 * k for k in range(17)]     # 4.0 .. 8.0
    for beta in betas:
        g2 = 2.0 * N / beta
        F = F_asym(g2, N)
        ainv_Lat = F
        ainv_MS = F / LAM_MS_OVER_LAT
        daD_lo = LAM_MS_OVER_LAT * DELTA_LO / F
        daD_hi = LAM_MS_OVER_LAT * DELTA_HI / F
        xi_lo = 1.0 / daD_lo
        xi_hi = 1.0 / daD_hi
        rows.append(["%.2f" % beta, "%.4f" % g2, "%.3e" % ainv_Lat,
                     "%.3e" % ainv_MS, "%.4f" % daD_lo, "%.4f" % daD_hi,
                     "%.3f" % xi_lo, "%.3f" % xi_hi])
    P.write_table(
        os.path.join(outdir, "t4_beta_axis.txt"),
        "SU(3) Wilson beta axis: a^-1 (2-loop scaling), a*Delta and xi/a "
        "(Lit: Lam_MS/Lam_lat=28.8, Delta/Lam_MS in [4,8])",
        ["beta", "g^2", "a^-1/Lam_lat", "a^-1/Lam_MS",
         "aD(Lo=4)", "aD(Hi=8)", "xi/a(Lo=4)", "xi/a(Hi=8)"],
        rows,
        note=("beta = 6/g^2.  a^-1/Lam_lat = F_asym(g^2) (2-loop asymptotic scaling).  "
              "aD = a*Delta = 28.8*(Delta/Lam_MS)/F_asym; xi/a = 1/(aD).\n"
              "LITERATURE INPUTS: Lam_MS/Lam_lat=28.8; Delta/Lam_MS in [4,8] (this last is "
              "the nonperturbative quantity one wants to prove -- used only as an illustrative "
              "input range; the exercise is CIRCULAR as a proof).\n"
              "At beta~6 the Wilson a is ~0.1-0.15 fm and 2-loop scaling has O(30-50%) "
              "violations, so treat beta numbers as O(1) indications."),
    )

    # ---------------- key points ----------------
    krows = []
    # (W) g^2 = 1
    beta_W = 2.0 * N / 1.0
    F1 = F_asym(1.0, N)
    ainv_MS_W = F1 / LAM_MS_OVER_LAT
    mu_GeV_W = math.exp(R.A(1.0 / MPI2, N, 3)) * LAM_MS_GEV
    krows.append(["g^2(a^-1)=1 (formal weak)", "%.2f" % beta_W, "1.000",
                  "%.3e" % ainv_MS_W, "%.3f" % mu_GeV_W,
                  "%.3f" % (LAM_MS_OVER_LAT * DELTA_LO / F1),
                  "%.3f" % (LAM_MS_OVER_LAT * DELTA_HI / F1), ""])
    # (S) a Delta = 1 for D = 4 and 8
    for D in (DELTA_LO, DELTA_HI):
        beta_S, g2_S = beta_for_aDelta_one(D, N)
        krows.append(["aD=1 (Delta=%.0f Lam)" % D, "%.3f" % beta_S, "%.4f" % g2_S,
                      "%.3e" % (F_asym(g2_S, N) / LAM_MS_OVER_LAT), "",
                      "1.000", "1.000", ""])
    # MS-bar g^2 at mu = Delta
    for D in (DELTA_LO, DELTA_HI):
        t = math.log(D)
        g2_3 = g2_at_t(t, N, 3)
        g2_2 = g2_at_t(t, N, 2)
        krows.append(["g^2(mu=%.0f Lam) MSbar" % D, "", "",
                      "", "%.3f" % (D * LAM_MS_GEV),
                      "", "", "g2_2l=%.3f g2_3l=%.3f" % (g2_2, g2_3)])
    P.write_table(
        os.path.join(outdir, "t4_keypoints.txt"),
        "Two candidate handoff loci on the beta axis + MS-bar g^2(mu=Delta)",
        ["locus", "beta", "g^2 bare", "a^-1/Lam_MS", "mu [GeV] (ill.)",
         "aD(D=4)", "aD(D=8)", "notes"],
        krows,
        note=("(W) g^2=1 is at beta=6 exactly (trivial).  (S) a*Delta=1 (q~1, strong side) is "
              "at beta=5.7 (Delta=4Lam) .. 6.3 (Delta=8Lam).  So on the beta axis the two "
              "windows COINCIDE near beta~6 -- an apparent overlap.  This is the lattice "
              "statement that the strong-weak crossover / scaling region is beta~5.7-6.0.\n"
              "BUT at beta~6 the lattice spacing is a~0.1-0.15 fm (UV), 2-loop scaling is "
              "violated at O(30-50%%), and Delta/Lam is a literature input -- so this is a "
              "plausibility indication, not a proof of q0>=q1.\n"
              "MS-bar g^2(mu=Delta): at the PHYSICAL handoff scale mu=Delta=(4..8)Lam the "
              "coupling is g^2~1.8-2.2 (2/3-loop), NOT 1.  Meanwhile g^2=1 sits at "
              "mu=158*Lam ~ %.0f GeV (with the illustrative Lam_MS=%.2f GeV).  Hence the "
              "'g^2=1' threshold is NOT invariant: it is ~20-40x above the physical handoff "
              "scale." % (mu_GeV_W, LAM_MS_GEV)),
    )
    res["beta_W"] = beta_W
    res["beta_S"] = {D: beta_for_aDelta_one(D, N)[0] for D in (DELTA_LO, DELTA_HI)}
    res["mu_g2one_GeV"] = mu_GeV_W
    res["g2_at_Delta"] = {D: (g2_at_t(math.log(D), N, 2),
                              g2_at_t(math.log(D), N, 3)) for D in (DELTA_LO, DELTA_HI)}

    # ---------------- plot: aD and xi/a vs beta ----------------
    c_lo, c_hi = [], []
    for k in range(121):
        beta = 4.0 + k * (4.0 / 120.0)
        g2 = 2.0 * N / beta
        F = F_asym(g2, N)
        c_lo.append((beta, LAM_MS_OVER_LAT * DELTA_LO / F))
        c_hi.append((beta, LAM_MS_OVER_LAT * DELTA_HI / F))
    P.svg_plot(
        os.path.join(outdir, "t4_beta_windows.svg"),
        [P.Curve("aD, Delta=4 Lam", c_lo, color=P.PALETTE[0]),
         P.Curve("aD, Delta=8 Lam", c_hi, color=P.PALETTE[1])],
        xlabel="beta = 6/g^2",
        ylabel="a*Delta   (= q = a/xi_phys)",
        title="q=a*Delta vs beta: q~1 (strong side) sits near beta~6",
        ylog=True,
        hlines=[(1.0, "q=1 (aD=1)", "#333333")],
        vlines=[(6.0, "beta=6 (g^2=1)", "#d62728")],
        legend_loc="upper right",
    )

    return res


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    r = run(here)
    print("=== lattice_beta.py (task 4) ===")
    print("literature inputs: Lam_MS/Lam_lat=%.1f, Delta/Lam_MS in [%.0f,%.0f]"
          % (LAM_MS_OVER_LAT, DELTA_LO, DELTA_HI))
    print("\n(W) g^2=1:            beta = %.2f" % r["beta_W"])
    for D in (DELTA_LO, DELTA_HI):
        print("(S) aD=1, D=%.0f Lam:   beta = %.3f" % (D, r["beta_S"][D]))
    print("\nMS-bar g^2(mu=Delta):")
    for D in (DELTA_LO, DELTA_HI):
        g2l, g3l = r["g2_at_Delta"][D]
        print("  mu=%.0f Lam -> g^2 = %.3f (2-loop), %.3f (3-loop)" % (D, g2l, g3l))
    print("\ng^2=1 sits at mu = 158*Lam ~ %.0f GeV (illustrative Lam_MS=0.22 GeV)"
          % r["mu_g2one_GeV"])
    print("\nwrote t4_beta_axis.txt, t4_keypoints.txt, t4_beta_windows.svg")

"""
delta_toy.py  --  Task 5.

The 2026-template finite-size quantity

      Delta(n) = 4 log Z(n x n) - log Z(2n x 2n)

(docs/yang-mills-prior-work.md sec.3.2), evaluated in a toy model:
a MASSIVE FREE GAUSSIAN on an L x L periodic lattice plus a SMALL quartic
perturbation.  We check (i) monotonicity, (ii) whether it crosses O(1) at a
physical size ~ 1/Delta (Delta = the mass gap = m here).

--------------------------------------------------------------------------
Toy model
--------------------------------------------------------------------------
Action  S = sum_x [ 1/2 (d_mu phi)^2 + 1/2 m^2 phi^2 + (lambda/4!) phi^4 ]
with lattice spacing 1, torus L x L.  Eigenvalues of -Laplacian:
    lambda_p = 4 sin^2(p_x/2) + 4 sin^2(p_y/2),  p = 2 pi k / L.

Gaussian (lambda=0):   log Z0(L) = -1/2 sum_p log(lambda_p + m^2).
Coincident propagator: G0(L)      = (1/L^2) sum_p 1/(lambda_p + m^2).
Leading quartic term:  log Z = log Z0 - (lambda/8) L^2 G0(L)^2 + O(lambda^2)
   [from  log Z = log Z0 - (lambda/4!) sum_x <phi_x^4>_0 ,  <phi^4>=3 G0^2].

The zero mode is INCLUDED (scalar).  A gauge theory would project it out; the
qualitative behaviour here does not depend on that choice.

--------------------------------------------------------------------------
What we find (numerically, then read off)
--------------------------------------------------------------------------
For the Gaussian, log Z(L) = -f_inf L^2 + (finite-size correction), and the
volume terms cancel exactly in 4 logZ(n)-logZ(2n).  The residual is the
single-winding finite-size correction, which scales as  exp(-m * n)  (linear
size n, mass gap m).  Hence

      Delta(n) ~= C * exp(- m n),        C ~= 6.5 - 6.6  (empirical, weakly
                                          m-dependent; numerically fitted).

So:
  * Delta(n) is STRICTLY DECREASING in n (monotone);
  * the log-slope equals the mass gap m  (the certificate does measure the gap);
  * Delta crosses O(1) at  n* = ln(C)/m ~= 1.9/m  -- same ORDER as the
    correlation length 1/m = 1/Delta (a factor ~1.9, i.e. inside the natural
    "physical size ~1/Delta" window).

The small quartic (lambda>0) shifts Delta DOWN (delta Delta < 0), pushing the
crossing to slightly larger n; for the small lambda used it stays monotone.
The O(lambda) formula is only schematic and its reliability degrades at large
n, so we report it as a small perturbation, not a precision result.

Honest limit: this is a FREE massive toy.  It demonstrates that the template
quantity Delta(n) is monotone and gap-sensitive in the cleanest case; it does
NOT demonstrate the interacting gauge-theory statements (exact-block
invariance, positivity -> gap bound, cutoff independence).
"""
import math

PI = math.pi


# --------------------------------------------------------------------------
#  Gaussian pieces on an L x L torus
# --------------------------------------------------------------------------
def logZ0_and_G0(L, m):
    s = 0.0
    g = 0.0
    mm = m * m
    for kx in range(L):
        px = 2.0 * PI * kx / L
        sx = 4.0 * math.sin(px / 2.0) ** 2
        for ky in range(L):
            py = 2.0 * PI * ky / L
            lam = sx + 4.0 * math.sin(py / 2.0) ** 2
            d = lam + mm
            s += math.log(d)
            g += 1.0 / d
    return -0.5 * s, g / (L * L)


def Delta(n, m, lam=0.0):
    """Delta(n) = 4 logZ(n)-logZ(2n), with O(lambda) quartic correction."""
    lz_n, g_n = logZ0_and_G0(n, m)
    lz_2n, g_2n = logZ0_and_G0(2 * n, m)
    if lam:
        lz_n -= (lam / 8.0) * n * n * g_n * g_n
        lz_2n -= (lam / 8.0) * (2 * n) ** 2 * g_2n * g_2n
    return 4.0 * lz_n - lz_2n


def crossing_n(ns, ds, target=1.0):
    """Linear interpolation in (n, ln Delta) for Delta=target (ds decreasing)."""
    for i in range(len(ns) - 1):
        d0, d1 = ds[i], ds[i + 1]
        if d0 >= target > d1:
            x0, x1 = ns[i], ns[i + 1]
            y0, y1 = math.log(d0), math.log(d1)
            yt = math.log(target)
            return x0 + (yt - y0) * (x1 - x0) / (y1 - y0)
    return float("nan")


def fit_exp(ns, ds):
    """Fit ln Delta = a + b n over points with 0.03 < Delta < 5; return C, slope."""
    pts = [(n, math.log(d)) for n, d in zip(ns, ds) if 0.03 < d < 5.0]
    if len(pts) < 2:
        return float("nan"), float("nan")
    n_ = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] ** 2 for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    denom = n_ * sxx - sx * sx
    b = (n_ * sxy - sx * sy) / denom
    a = (sy - b * sx) / n_
    return math.exp(a), b


def run(outdir="."):
    import os
    import _plot as P

    ns = [2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
    results = {"tables": {}, "summary": {}}

    # ---- main table: m=0.1, three lambdas ----------------------------
    m0 = 0.1
    rows = []
    d0_list, d1_list = [], []
    for n in ns:
        d0 = Delta(n, m0, 0.0)
        d1 = Delta(n, m0, 0.001)
        d2 = Delta(n, m0, 0.005)
        d0_list.append(d0)
        d1_list.append(d2)
        rows.append(["%3d" % n, "%.5e" % d0, "%.5e" % d1, "%.5e" % d2,
                     "%.5f" % math.log(d0) if d0 > 0 else "nan"])
    P.write_table(
        os.path.join(outdir, "t5_delta_vs_n.txt"),
        "Delta(n)=4logZ(n)-logZ(2n) for massive Gaussian + quartic, m=0.1 (xi=1/m=10)",
        ["n", "Delta(lam=0)", "Delta(lam=0.001)", "Delta(lam=0.005)", "ln Delta(lam=0)"],
        rows,
        note=("logZ0 = -1/2 sum_p log(lambda_p+m^2); quartic via the O(lambda) term "
              "-(lambda/8)L^2 G0(L)^2.\n"
              "Delta(n) is strictly DECREASING; it crosses 1 at n* ~ 1.9/m (see summary). "
              "For lambda>0 the O(lambda) correction is NEGATIVE, pushing n* slightly up; for "
              "genuinely small lambda (0.001-0.005) the curve stays positive and monotone.\n"
              "HONESTY: at lambda=0.01 the O(lambda) shift already drives Delta<0 at n=2 and "
              "large n -- i.e. the linear-in-lambda truncation has broken down (O(lambda^2) "
              "matters), so 0.01 is NOT a valid 'small' quartic here.  The quartic column is a "
              "schematic perturbation, not a precision result."),
    )

    # ---- summary across masses and lambdas ---------------------------
    srows = []
    for m in (0.2, 0.1, 0.05):
        ns_local = [n for n in ns if n <= (128 if m >= 0.1 else 128)]
        ds0 = [Delta(n, m, 0.0) for n in ns_local]
        ncross = crossing_n(ns_local, ds0)
        C, slope = fit_exp(ns_local, ds0)
        # small-lambda crossing
        ds_l = [Delta(n, m, 0.001) for n in ns_local]
        ncross_l = crossing_n(ns_local, ds_l)
        # monotonicity check
        mono = all(ds0[i] > ds0[i + 1] for i in range(len(ds0) - 1))
        srows.append(["%.2f" % m, "%.1f" % (1.0 / m), "%.3f" % C, "%.4f" % slope,
                      "%.3f" % (-slope / m), "%.2f" % ncross,
                      "%.2f" % (ncross / (1.0 / m)), "%.2f" % ncross_l,
                      "yes" if mono else "NO"])
        results["summary"][m] = {"C": C, "slope": slope, "n_cross": ncross,
                                 "n_cross_lam": ncross_l, "monotone": mono}
    P.write_table(
        os.path.join(outdir, "t5_summary.txt"),
        "Delta(n) toy: exponential fit, O(1) crossing, and comparison to xi=1/m",
        ["m=Delta", "xi=1/m", "C (fit)", "slope (fit)", "-slope/m",
         "n*(Delta=1)", "n*/(1/m)", "n*(lam=0.001)", "monotone?"],
        srows,
        note=("Fit ln Delta = ln C + slope*n over 0.03<Delta<5.  Empirically slope = -m (to a "
              "few %) and C ~ 6.5-6.6: Delta(n) ~ C exp(-m n).  So the log-Delta slope directly "
              "equals the mass gap, and Delta crosses O(1) at n* ~ 1.9/m, i.e. same order as the "
              "correlation length 1/Delta=1/m (factor ~1.9).  This is the cleanest confirmation "
              "of the template's design intent: Delta(n) is monotone and gap-sensitive.\n"
              "Caveat: free massive toy only; no interacting/gauge-theory claim is made."),
    )

    # ---- plot: Delta(n) log-y for m=0.1, lam=0 and 0.01 --------------
    pts0 = [(n, Delta(n, 0.1, 0.0)) for n in ns]
    pts1 = [(n, Delta(n, 0.1, 0.005)) for n in ns]
    P.svg_plot(
        os.path.join(outdir, "t5_delta.svg"),
        [P.Curve("lambda=0", pts0, color=P.PALETTE[0]),
         P.Curve("lambda=0.005 (O(lambda))", pts1, color=P.PALETTE[1])],
        xlabel="n  (box linear size, lattice units)",
        ylabel="Delta(n)",
        title="Template Delta(n)=4logZ(n)-logZ(2n): monotone, crosses 1 near n~1.9/m",
        ylog=True,
        xlog=False,
        hlines=[(1.0, "Delta=1", "#333333")],
        vlines=[(1.0 / 0.1, "xi=1/m=10", "#2ca02c"),
                (2.079 / 0.1, "mu*=8Lam (from task 3)", "#9467bd")],
        legend_loc="upper right",
    )

    return results


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    r = run(here)
    print("=== delta_toy.py (task 5) ===")
    print("Delta(n) ~ C exp(-m n); crossing Delta=1 at n* ~ ln(C)/m (same order as 1/m)")
    print("\n%-8s %8s %8s %9s %9s %11s %10s %8s" %
          ("m=Delta", "xi=1/m", "C", "slope", "-slope/m", "n*(D=1)", "n*/(1/m)", "mono"))
    for m in (0.2, 0.1, 0.05):
        d = r["summary"][m]
        print("%-8.2f %8.1f %8.3f %9.4f %9.3f %11.2f %10.2f %8s" %
              (m, 1.0 / m, d["C"], d["slope"], -d["slope"] / m,
               d["n_cross"], d["n_cross"] / (1.0 / m),
               "yes" if d["monotone"] else "NO"))
    print("\nwrote t5_delta_vs_n.txt, t5_summary.txt, t5_delta.svg")

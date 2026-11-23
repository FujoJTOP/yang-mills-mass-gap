"""
rg_handoff.py  --  Pure Yang-Mills SU(N) beta function, RG flow, and the
RG-invariant scale Lambda (task 1).

--------------------------------------------------------------------------
CONVENTION (stated explicitly; checked against standard QCD nf=0 values)
--------------------------------------------------------------------------
We use

    beta(g) = dg/d ln(mu)
            = - b0 g^3/(16 pi^2)
              - b1 g^5/(16 pi^2)^2
              - b2 g^7/(16 pi^2)^3 - ...

with the MS-bar coefficients for a pure SU(N) gauge theory

    b0 = 11 N / 3
    b1 = 34 N^2 / 3
    b2 = 2857 N^3 / 54

For N=3 this gives b0=11, b1=102, b2=2857/2 = 1428.5, which are exactly the
textbook QCD values for nf=0 in this normalisation.  (b0 and b1 are
scheme-independent for a single coupling in a mass-independent scheme; b2
is scheme-dependent in general, and this is used in scheme_dependence.py.)

--------------------------------------------------------------------------
Working variable
--------------------------------------------------------------------------
It is numerically cleanest to use   a = g^2/(16 pi^2)  (so g^2 = 16 pi^2 a).
Then

    da/dt = - 2 a^2 ( b0 + b1 a + b2 a^2 ),      t = ln(mu).

--------------------------------------------------------------------------
RG-invariant scale (derived, not looked up)
--------------------------------------------------------------------------
The running is  dt = dg/beta = - da / (2 a^2 P(a)),  P(a)=b0+b1 a+b2 a^2.
Integrating gives the RG invariant  Lambda = mu * exp(-int^g dg'/beta(g')):

    ln(mu/Lambda) = A(a)
      A_2loop(a) = 1/(2 b0 a)
                   + (b1/(2 b0^2)) ln[ a / (1 + (b1/b0) a) ]
      A_3loop(a) = A_2loop(a)
                   + (b2/2) * int_0^a da' / [ P(a') (b0 + b1 a') ]

(The 3-loop correction follows from dt/da - dt_2loop/da = b2/(2 P P2).)

The identity d A/dt = 1 is self-checked here by integrating the ODE
da/dt = -2 a^2 P(a) and confirming that A(a(t)) - t is constant.
--------------------------------------------------------------------------
"""
import math

PI = math.pi
MPI2 = 16.0 * PI * PI          # 16 pi^2


# --------------------------------------------------------------------------
#  coefficients / beta function
# --------------------------------------------------------------------------
def coeffs(N):
    b0 = 11.0 * N / 3.0
    b1 = 34.0 * N * N / 3.0
    b2 = 2857.0 * N ** 3 / 54.0
    return b0, b1, b2


def beta_g(g, N, loops=3):
    """beta(g) = dg/d ln mu  in the convention above."""
    b0, b1, b2 = coeffs(N)
    g2 = g * g
    val = -b0 * g * g2 / MPI2
    if loops >= 2:
        val -= b1 * g * g2 * g2 / MPI2 ** 2
    if loops >= 3:
        val -= b2 * g * g2 * g2 * g2 / MPI2 ** 3
    return val


def beta_a(a, N, loops=3):
    """da/dt for a = g^2/16pi^2."""
    b0, b1, b2 = coeffs(N)
    P = b0
    if loops >= 2:
        P += b1 * a
    if loops >= 3:
        P += b2 * a * a
    return -2.0 * a * a * P


# --------------------------------------------------------------------------
#  quadrature (Simpson, stdlib only)
# --------------------------------------------------------------------------
def simpson(f, x0, x1, n=400):
    if x1 == x0:
        return 0.0
    if n % 2:
        n += 1
    h = (x1 - x0) / n
    s = f(x0) + f(x1)
    for i in range(1, n):
        s += (4.0 if i % 2 else 2.0) * f(x0 + i * h)
    return s * h / 3.0


# --------------------------------------------------------------------------
#  RG-invariant scale:  ln(mu/Lambda) = A(a)
# --------------------------------------------------------------------------
def A_2loop(a, N):
    b0, b1, _ = coeffs(N)
    ca = (b1 / b0) * a
    return 1.0 / (2.0 * b0 * a) + (b1 / (2.0 * b0 * b0)) * math.log(a / (1.0 + ca))


def A(a, N, loops=3):
    """ln(mu/Lambda) as a function of a = g^2/16pi^2 (weak-coupling branch)."""
    if loops == 1:
        b0, _, _ = coeffs(N)
        return 1.0 / (2.0 * b0 * a)
    b0, b1, b2 = coeffs(N)
    val = A_2loop(a, N)
    if loops >= 3:
        def integrand(x):
            P = b0 + b1 * x + b2 * x * x
            P2 = b0 + b1 * x
            return b2 / (P * P2)
        val += 0.5 * simpson(integrand, 0.0, a, 800)
    return val


def lambda_at_threshold(target_g2, N, loops=3):
    """Return t = ln(mu/Lambda) at which g^2 = target_g2."""
    a = target_g2 / MPI2
    return A(a, N, loops)


# --------------------------------------------------------------------------
#  self-check:  integrate ODE and confirm A(a(t)) - t = const
# --------------------------------------------------------------------------
def selfcheck(N=3, loops=3, a_lo=1e-3, a_hi=0.09, dt=-0.05):
    """Integrate da/dt = -2 a^2 P(a) and confirm A(a(t)) - t is constant.

    The analytic A(a) is the exact antiderivative of dt/da, so along the
    numerical flow the quantity A(a(t)) - t must be the flow constant
    (-ln Lambda).  Any drift measures numerical error only.
    """
    # start at a_lo, integrate downward in t until a exceeds a_hi
    a = a_lo
    t = 0.0
    vals = []
    for _ in range(20000):
        def f(aa):
            return beta_a(aa, N, loops)
        k1 = f(a)
        k2 = f(a + dt / 2 * k1)
        k3 = f(a + dt / 2 * k2)
        k4 = f(a + dt * k3)
        a += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        if a <= 0:
            break
        if a > a_hi:
            break
        vals.append(A(a, N, loops) - t)
    if not vals:
        return 0.0
    return max(vals) - min(vals)


# --------------------------------------------------------------------------
#  driver
# --------------------------------------------------------------------------
def run(outdir="."):
    import os
    import _plot as P

    results = {"thresholds": {}, "selfcheck": {}}

    # thresholds of interest, expressed as g^2
    targets = [
        ("g2=1",        1.0),
        ("g2=4pi",      4.0 * PI),
        ("lambda=g2*N=1", None),    # filled per N -> g^2 = 1/N
    ]

    for N in (2, 3):
        rows = []
        results["thresholds"][N] = {}
        for name, g2 in targets:
            tg2 = (1.0 / N) if g2 is None else g2
            row = [name, "%.5f" % tg2]
            tvals = {}
            for loops in (1, 2, 3):
                t = lambda_at_threshold(tg2, N, loops)
                tvals[loops] = t
                row.append("%.4f" % t)
            row.append("%.3f" % tvals[3] if 3 in tvals else "")
            rows.append(row)
            results["thresholds"][N][name] = {
                "g2": tg2,
                "t_1loop": tvals[1],
                "t_2loop": tvals[2],
                "t_3loop": tvals[3],
            }
        P.write_table(
            os.path.join(outdir, "t1_thresholds_SU%d.txt" % N),
            "SU(%d) pure YM: t = ln(mu/Lambda) at coupling thresholds" % N,
            ["criterion", "g^2", "t(1-loop)", "t(2-loop)", "t(3-loop)"],
            rows,
            note=("Convention: beta = -b0 g^3/(16pi^2) - b1 g^5/(16pi^2)^2 - b2 g^7/(16pi^2)^3; "
                  "b0=11N/3, b1=34N^2/3, b2=2857N^3/54.\n"
                  "t = ln(mu/Lambda) is the RG-invariant logarithm; it is the SAME for every "
                  "initial condition (Lambda is fixed by the theory).\n"
                  "The spread across loop orders (esp. at g^2=4pi) is a truncation diagnostic."),
        )

    # sanity: print the b0,b1,b2 for N=3 and compare with QCD nf=0
    b0, b1, b2 = coeffs(3)
    results["b_coeffs_SU3"] = (b0, b1, b2)

    # self-check identity
    for N in (2, 3):
        dev = selfcheck(N, 3)
        results["selfcheck"][N] = dev

    # plot g^2 vs t for N=2,3 (3-loop)
    curves = []
    for i, N in enumerate((2, 3)):
        pts = []
        # sample a logarithmically from tiny to 0.12
        amax = 0.12
        npts = 500
        for k in range(npts + 1):
            frac = k / npts
            a = 10.0 ** (math.log10(1e-6) + frac * (math.log10(amax) - math.log10(1e-6)))
            t = A(a, N, 3)
            g2 = MPI2 * a
            pts.append((t, g2))
        curves.append(P.Curve("SU(%d), 3-loop" % N, pts, color=P.PALETTE[i]))

    vlines = [
        (0.0, "mu=Lambda", "#555555"),
    ]
    hlines = [
        (1.0, "g^2=1", "#d62728"),
        (4.0 * PI, "g^2=4pi (alpha=1)", "#2ca02c"),
        (1.0 / 3.0, "lambda=g^2 N=1 (N=3)", "#9467bd"),
    ]
    P.svg_plot(
        os.path.join(outdir, "t1_running_g2.svg"),
        curves,
        xlabel="t = ln(mu / Lambda)",
        ylabel="g^2",
        title="Pure YM running coupling g^2(mu) vs t=ln(mu/Lambda)  (MS-bar, 3-loop)",
        ylog=True,
        vlines=vlines,
        hlines=hlines,
        legend_loc="upper right",
    )

    # also write a raw data table for the N=3 curve
    rows = []
    for k in range(41):
        frac = k / 40.0
        a = 10.0 ** (math.log10(1e-5) + frac * (math.log10(0.12) - math.log10(1e-5)))
        rows.append(["%.6e" % (MPI2 * a), "%.5f" % A(a, 3, 3),
                     "%.5f" % A(a, 3, 2), "%.5f" % A(a, 3, 1)])
    P.write_table(
        os.path.join(outdir, "t1_flow_SU3.txt"),
        "SU(3) g^2 -> t=ln(mu/Lambda) (3/2/1-loop)",
        ["g^2", "t(3-loop)", "t(2-loop)", "t(1-loop)"],
        rows,
    )

    return results


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    r = run(here)
    print("=== rg_handoff.py (task 1) ===")
    print("b-coefficients SU(3): b0=%.4f b1=%.4f b2=%.4f" % r["b_coeffs_SU3"])
    for N in (2, 3):
        print("\nSU(%d): t=ln(mu/Lambda) at thresholds (3-loop; [1-loop])" % N)
        for name, d in r["thresholds"][N].items():
            print("  %-14s g^2=%.5f  t=%.4f   [1-loop %.4f]"
                  % (name, d["g2"], d["t_3loop"], d["t_1loop"]))
    print("\nself-check max|A(a(t))-t| deviation: N=2 %.2e  N=3 %.2e"
          % (r["selfcheck"][2], r["selfcheck"][3]))
    print("\nwrote t1_thresholds_SU2.txt, t1_thresholds_SU3.txt, t1_flow_SU3.txt, t1_running_g2.svg")

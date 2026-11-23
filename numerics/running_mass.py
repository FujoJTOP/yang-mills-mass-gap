"""
running_mass.py  --  Task 3.

The control quantity q(mu) := m(mu)/mu  and its unique crossing of 1.

--------------------------------------------------------------------------
RG equation and the one stated assumption
--------------------------------------------------------------------------
Pure YM has no bare mass.  Following research/physics_handoff.md we DEFINE
m(mu) := 1/xi(mu), the scale-mu effective two-point decay rate (inverse
correlation length).  Its RG equation is

      mu dm/dmu = - gamma_m m ,        0 < gamma_m < 1   (physical range).

This gives  m(mu) ~ mu^{-gamma_m}  and hence

      q(mu) = m(mu)/mu ~ mu^{-(1+gamma_m)} ,
      dq/dmu = -(1+gamma_m) m/mu^2 < 0   (strictly decreasing in mu,
                                          strictly increasing as mu -> 0).

ASSUMPTION (illustrative, flagged everywhere it matters): we take
gamma_m = 0.4 (constant) as a representative value in (0,1).  Real QCD-like
anomalous dimensions run with g; a constant gamma is the minimal toy that
still has the qualitative behaviour.  The monotonicity and the uniqueness of
the q=1 crossing do NOT depend on the value of gamma_m, only on gamma_m>0.

--------------------------------------------------------------------------
What is computed
--------------------------------------------------------------------------
* We integrate the coupled system (a(t), m(t)), t = ln mu, by RK4:
        da/dt = beta_a(a)  (3-loop, SU(3), MS-bar)
        dm/dt = -gamma_m m
  from a UV point t0 = 8 (mu = Lambda e^8) downward.
* m is anchored so that the crossing sits at mu* = c Lambda (c = 1, 2, 4, 8);
  the anchor is the single free O(1) number m_gap/Lambda of dimensional
  transmutation.
* We locate the crossing q = 1 numerically, verify mu* = m(mu*) there, verify
  q is strictly monotone on a log axis, and report q at the places that
  matter: the formal weak boundary g^2 = 1 (t=5.06) and the PHYSICAL handoff
  scale mu ~ Delta = (4..8) Lambda (where g^2 ~ 1.8-2.2).

--------------------------------------------------------------------------
The honest framing
--------------------------------------------------------------------------
q is a COORDINATE: it is built from m, and m/ Lambda is itself the unknown
nonperturbative number.  So q does not by itself prove anything.  What is
actually to be proved is the quantitative overlap

        q0 := q(weak effective boundary)  >=  q1 := q(strong boundary) ~ 1.

This script shows where q=1 sits and how far q is from 1 at the nominal
boundaries; it does NOT establish q0 >= q1.
"""
import math

import rg_handoff as R

MPI2 = R.MPI2
GAMMA = 0.4                      # <-- the stated illustrative assumption


# --------------------------------------------------------------------------
#  local 2-loop inverse of A  (A is strictly monotone for a>0 at 2 loops)
# --------------------------------------------------------------------------
def A2(a, N=3):
    return R.A_2loop(a, N)


def invA2(t, N=3, lo=1e-12, hi=0.5):
    if A2(lo, N) < t or A2(hi, N) > t:
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if A2(mid, N) > t:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def g2_of_t_2loop(t, N=3):
    a = invA2(t, N)
    return MPI2 * a if a == a else float("nan")


# --------------------------------------------------------------------------
#  RK4 integration of (a, m) from t0 downward
# --------------------------------------------------------------------------
def integrate(gamma, t0=8.0, t1=-1.0, dt=-0.002, m_at_t0=1.0, N=3):
    def f(t, a, m):
        return R.beta_a(a, N, 3), -gamma * m

    a0 = invA2(t0, N)
    a, m = a0, m_at_t0
    ts, aas, ms = [t0], [a0], [m_at_t0]
    t = t0
    nsteps = int(round((t0 - t1) / (-dt)))
    for _ in range(nsteps):
        k1a, k1m = f(t, a, m)
        k2a, k2m = f(t + dt / 2, a + dt / 2 * k1a, m + dt / 2 * k1m)
        k3a, k3m = f(t + dt / 2, a + dt / 2 * k2a, m + dt / 2 * k2m)
        k4a, k4m = f(t + dt, a + dt * k3a, m + dt * k3m)
        a += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        m += dt / 6.0 * (k1m + 2 * k2m + 2 * k3m + k4m)
        t += dt
        ts.append(t)
        aas.append(a)
        ms.append(m)
    return ts, aas, ms


def q_of(ts, ms):
    return [m / math.exp(t) for t, m in zip(ts, ms)]


def find_crossing(ts, qs):
    """Locate q=1 by bisection on the (t-decreasing) grid; q increasing as t falls."""
    # ts is decreasing from t0 to t1; q increases as t decreases
    for i in range(len(ts) - 1):
        q0, q1 = qs[i], qs[i + 1]     # q1 > q0 (t_{i+1} < t_i)
        if q0 < 1.0 <= q1:
            t0, t1 = ts[i], ts[i + 1]
            for _ in range(80):
                tm = 0.5 * (t0 + t1)
                # interpolate q at tm linearly between endpoints (grid is fine)
                f = (tm - t0) / (t1 - t0)
                qm = q0 + f * (q1 - q0)
                if qm < 1.0:
                    t0 = tm
                    q0 = qm
                else:
                    t1 = tm
                    q1 = qm
            return 0.5 * (t0 + t1)
    return float("nan")


def run(outdir="."):
    import os
    import _plot as P

    N = 3
    results = {"gamma": GAMMA, "crossings": {}, "tables": {}}

    # anchor scales c = mu*/Lambda  (crossing at mu*=c Lambda)
    anchors = [1.0, 2.0, 4.0, 8.0]
    # t(g^2=1) in 3-loop MS-bar (from rg_handoff)
    t_g2one = R.A(1.0 / MPI2, N, 3)

    for c in anchors:
        tstar = math.log(c)
        # anchor m(t0) so that m(tstar)=e^{tstar}
        m0 = math.exp(tstar) * math.exp(GAMMA * (tstar - 8.0))
        ts, aas, ms = integrate(GAMMA, t0=8.0, t1=-0.6, dt=-0.002,
                                m_at_t0=m0, N=N)
        qs = q_of(ts, ms)
        tc = find_crossing(ts, qs)
        # value of m / exp(t) at the crossing
        mc = m0 * math.exp(-GAMMA * (tc - 8.0))
        muc = math.exp(tc)
        # monotonicity of q: grid runs t=8 -> -1 (decreasing t), so q must INCREASE
        dq = [qs[i + 1] - qs[i] for i in range(len(qs) - 1)]  # all >0
        # q at the formal weak boundary g^2=1
        q_at_g2one = m0 * math.exp(-GAMMA * (t_g2one - 8.0)) / math.exp(t_g2one)
        # q at the physical handoff mu = 4 Lambda and 8 Lambda
        q_at_4 = m0 * math.exp(-GAMMA * (math.log(4) - 8.0)) / 4.0
        q_at_8 = m0 * math.exp(-GAMMA * (math.log(8) - 8.0)) / 8.0
        # q at g^2 = 1.8 (the physical handoff coupling found in task 4)
        t_g2_18 = R.A(1.8 / MPI2, N, 3)
        q_at_g2_18 = m0 * math.exp(-GAMMA * (t_g2_18 - 8.0)) / math.exp(t_g2_18)
        results["crossings"][c] = {
            "tstar": tstar, "t_cross": tc, "mu_star_over_Lam": math.exp(tc),
            "m_at_cross": mc, "mu_at_cross": muc,
            "q_at_g2one": q_at_g2one,
            "q_at_4Lam": q_at_4, "q_at_8Lam": q_at_8,
            "q_at_g2_18": q_at_g2_18,
            "monotone_ok": all(d > 0 for d in dq),
            "m_minus_mu_at_cross": mc - muc,
        }
        if c == 8.0:
            # keep the full curve for the table / plot
            results["tables"]["ts"] = ts
            results["tables"]["g2s"] = [MPI2 * a for a in aas]
            results["tables"]["ms"] = ms
            results["tables"]["qs"] = qs

    # ---- main table: q(t) curve for the c=8 anchor --------------------
    ts = results["tables"]["ts"]
    g2s = results["tables"]["g2s"]
    ms = results["tables"]["ms"]
    qs = results["tables"]["qs"]
    # sample every 500 steps (dt=-0.002 -> 1 unit of t per 500 steps)
    rows = []
    last_t = None
    for i in range(0, len(ts), 500):
        t = ts[i]
        region = "UV (q<<1)" if qs[i] < 0.5 else ("handoff (q~1)" if qs[i] < 2.0 else "IR (q>1)")
        rows.append(["%+.2f" % t, "%.4e" % math.exp(t), "%.4e" % g2s[i],
                     "%.4e" % ms[i], "%.4e" % qs[i], region])
        last_t = t
    P.write_table(
        os.path.join(outdir, "t3_running_mass.txt"),
        "SU(3), gamma_m=%.1f (ASSUMED const): q=m/mu along the flow (anchor c=mu*/Lam=8)" % GAMMA,
        ["t=ln(mu/Lam)", "mu/Lam", "g^2(mu) 2l", "m/Lam", "q=m/mu", "region"],
        rows,
        note=("m(mu) from mu dm/dmu = -gamma_m m with gamma_m=%.1f (illustrative constant, "
              "0<gamma_m<1).  g^2 from 2-loop MS-bar running (reference only).\n"
              "q=m/mu is strictly increasing as mu decreases (strictly decreasing in t); it "
              "crosses q=1 exactly once, at mu*=m(mu*).  For this anchor the crossing is at "
              "mu*=8*Lambda (t*=ln8=2.079).\n"
              "q is a COORDINATE (it contains the unknown m/Lambda); the real open problem is "
              "the quantitative overlap q0>=q1, not the existence of this crossing." % GAMMA),
    )

    # ---- crossing / boundary table ------------------------------------
    rows = []
    for c in anchors:
        d = results["crossings"][c]
        rows.append(["%.1f" % c, "%.4f" % d["tstar"], "%.4f" % d["t_cross"],
                     "%.4f" % d["mu_star_over_Lam"], "%.2e" % abs(d["m_minus_mu_at_cross"]),
                     "%.3e" % d["q_at_g2one"], "%.3f" % d["q_at_g2_18"],
                     "%.3f" % d["q_at_4Lam"], "%.3f" % d["q_at_8Lam"],
                     "yes" if d["monotone_ok"] else "NO"])
    P.write_table(
        os.path.join(outdir, "t3_crossing.txt"),
        "q=1 crossing and q at the candidate boundaries (SU(3), gamma_m=%.1f assumed)" % GAMMA,
        ["mu*/Lam", "t*=ln(mu*/Lam)", "t_cross(num)", "mu*/Lam(num)",
         "|m-mu|@cross", "q@(g^2=1)", "q@(g^2=1.8)", "q@(mu=4Lam)", "q@(mu=8Lam)", "q monotone?"],
        rows,
        note=("Rows: anchor mu*=c*Lambda sets the crossing (m_gap/Lambda = c is the one free "
              "O(1) number).  t_cross(num) is found by RK4 integration + bisection and matches "
              "ln(c); |m-mu|@cross confirms mu*=m(mu*) at the crossing.\n"
              "The column q@(g^2=1) is the value of q at the FORMAL weak boundary t=5.06 "
              "(g^2=1).  It is tiny (1e-4 .. 1e-2) for every anchor: the g^2=1 line lies deep "
              "in the q<<1 region and CANNOT be the handoff boundary.\n"
              "The column q@(g^2=1.8) is q at the PHYSICAL handoff scale g^2~1.8 (task 4).  "
              "For the literature anchor mu*=(4..8)Lam there q is O(0.4..1): the physical "
              "handoff lives where q~1, i.e. at mu~Delta~m_gap, "
              "consistent with the dimensional-transmutation picture.  This is why q, not g^2, "
              "is the right control quantity.\n"
              "Caveat: with gamma_m constant the power law m~mu^{-gamma} is only the weak-side "
              "behaviour; in the deep IR m freezes at m_gap, which makes q grow even faster, so "
              "the q>1 strong side is unaffected.  These numbers are schematic, not a proof."),
    )

    # ---- plot: q(t) on a log axis, several gamma ----------------------
    curves = []
    for idx, gam in enumerate((0.2, 0.4, 0.7)):
        tstar = math.log(8.0)
        m0 = math.exp(tstar) * math.exp(gam * (tstar - 8.0))
        ts2, _, ms2 = integrate(gam, t0=8.0, t1=-1.0, dt=-0.005,
                                m_at_t0=m0, N=N)
        qs2 = q_of(ts2, ms2)
        pts = [(t, q) for t, q in zip(ts2, qs2) if 1e-6 < q < 1e6]
        curves.append(P.Curve("gamma_m=%.1f" % gam, pts,
                              color=P.PALETTE[idx]))
    P.svg_plot(
        os.path.join(outdir, "t3_q_log.svg"),
        curves,
        xlabel="t = ln(mu / Lambda)",
        ylabel="q = m(mu)/mu",
        title="Running-mass ratio q=m/mu is monotone and crosses 1 once (anchor mu*=8 Lambda)",
        ylog=True,
        hlines=[(1.0, "q=1", "#333333")],
        vlines=[(0.0, "mu=Lambda", "#555555"),
                (math.log(4), "mu=4Lam", "#2ca02c"),
                (math.log(8), "mu=8Lam (q=1)", "#9467bd"),
                (t_g2one, "g^2=1", "#d62728")],
        legend_loc="upper right",
    )

    return results


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    r = run(here)
    print("=== running_mass.py (task 3) ===")
    print("gamma_m = %.1f  (assumed constant, 0<gamma<1)" % r["gamma"])
    print("\n%-8s %10s %10s %12s %12s %10s" %
          ("mu*/Lam", "t*", "t_cross", "mu*/Lam(num)", "|m-mu|@cross", "monotone"))
    for c in (1.0, 2.0, 4.0, 8.0):
        d = r["crossings"][c]
        print("%-8.1f %10.4f %10.4f %12.4f %12.2e %10s" %
              (c, d["tstar"], d["t_cross"], d["mu_star_over_Lam"],
               abs(d["m_minus_mu_at_cross"]), "yes" if d["monotone_ok"] else "NO"))
    d = r["crossings"][8.0]
    print("\nFor anchor mu*=8 Lambda: q at g^2=1      = %.3e  (formal weak boundary, q<<1)" % d["q_at_g2one"])
    print("                          q at g^2=1.8    = %.3f   (PHYSICAL handoff, q~1)" % d["q_at_g2_18"])
    print("                          q at mu=4 Lam   = %.3f" % d["q_at_4Lam"])
    print("                          q at mu=8 Lam   = %.3f" % d["q_at_8Lam"])
    print("\nwrote t3_running_mass.txt, t3_crossing.txt, t3_q_log.svg")

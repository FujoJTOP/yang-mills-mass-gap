"""kt_e_identity.py -- KT-RHO2: the Wronskian identity that localises the DD-class bound.

  Identity [T-standard]:  for  -u'' + q u = lam u  on (0,1), Dirichlet, with u1 > 0 a ground
  state and u2 an eigenfunction with lam2 > lam1,

        lam2 - lam1 = ( INT (W/u1)^2 dx ) / ( INT u2^2 dx ),      W = u1 u2' - u1' u2 = const.

  Hence  rho^2 <= 4  <=>  INT (W/u1)^2 <= 3 lam1 INT u2^2  -- i.e. the whole DD-class bound is
  ONE inequality about the ground-state log-derivative (u1'/u1).
Verified numerically here on random q >= 0 (DD and ND), and on a two-step string.
Writes kt_e_identity.json + .log.
"""
import json, numpy as np
from scipy.linalg import eigh_tridiagonal

JSON, LOG = "kt_e_identity.json", []
OUT = {}


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


def flush():
    with open(JSON, "w", encoding="utf-8") as f:
        json.dump(OUT, f, indent=1)
    with open("kt_e_identity.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def solve(N=800, q=None, w=None, bc="DD"):
    h = 1.0 / N
    x = np.linspace(0, 1, N + 1)
    qn = np.zeros(N + 1) if q is None else np.asarray(q(x), float)
    idx = {"DD": np.arange(1, N), "ND": np.arange(0, N), "NN": np.arange(0, N + 1)}[bc]
    m = len(idx)
    diag = np.zeros(m); off = np.zeros(m - 1)
    for a, j in enumerate(idx):
        if j == 0 and bc in ("NN", "ND"):
            diag[a] = 1 / h ** 2 + qn[0]; off[a] = -1 / h ** 2; continue
        d = qn[j] + (1 / h ** 2 if j >= 1 else 0.0)
        if j < N:
            d += 1 / h ** 2
        elif bc == "NN":
            d = qn[N] + 1 / h ** 2
        diag[a] = d
        if a < m - 1:
            off[a] = -1 / h ** 2
    if w is not None:
        ww = np.asarray(w(x), float)[idx]
        diag = diag / ww; off = off / np.sqrt(ww[:-1] * ww[1:])
    lam, V = eigh_tridiagonal(diag, off)
    return x[idx], lam, V


P("=" * 100)
P("E1  verify  lam2-lam1 = INT (W/u1)^2 / INT u2^2      (W = const Wronskian)")
P("=" * 100)
rng = np.random.default_rng(11)
res = {}
for trial in range(6):
    kind = ["free", "parab", "smooth", "bump", "two-step"][trial % 5]
    if kind == "free":
        q = None
    elif kind == "parab":
        q = lambda t: 5.0 * (t - 0.5) ** 2
    elif kind == "smooth":
        a = rng.exponential(1.0, 5); ph = rng.uniform(0, 2 * np.pi, 5)
        q = (lambda t, a=a, ph=ph: sum(a[i] * (1 + np.cos(2 * np.pi * (i + 1) * t + ph[i]))
                                      for i in range(5)))
    elif kind == "bump":
        q = lambda t: 40.0 * np.exp(-((t - 0.3) / 0.05) ** 2)
    else:
        q = lambda t: 30.0 * ((t < 0.2) | (t > 0.8))
    x, lam, V = solve(N=800, q=q)
    u1, u2 = V[:, 0], V[:, 1]
    W = u1 * np.gradient(u2, x) - np.gradient(u1, x) * u2
    lhs = lam[1] - lam[0]
    rhs = np.trapezoid((W / u1) ** 2, x) / np.trapezoid(u2 ** 2, x)
    d = dict(kind=kind, lam1=float(lam[0]), lam2=float(lam[1]),
             lhs=float(lhs), rhs=float(rhs), rel_err=float(abs(rhs - lhs) / lhs),
             W_spread=float(np.ptp(W) / (abs(np.mean(W)) + 1e-30)))
    res[f"trial{trial}_{kind}"] = d
    P(f"  {kind:<9s} lam1={lam[0]:.6f} lam2={lam[1]:.6f}  lhs={lhs:.8f} rhs={rhs:.8f}"
      f"  rel.err={d['rel_err']:.2e}   (W is constant to {d['W_spread']:.1e})")
OUT["E1_identity"] = res
flush()

P("=" * 100)
P("E2  what the identity says about the constant:   rho^2 <= 4  <=>  INT(W/u1)^2 <= 3 lam1 INT u2^2")
P("=" * 100)
P("    equivalently  3 lam1 >= INT (W/u1)^2 / INT u2^2 = lam2-lam1 ,   i.e. lam2 <= 4 lam1")
for bc in ["DD", "ND"]:
    x, lam, V = solve(N=1200, bc=bc)
    u1, u2 = V[:, 0], V[:, 1]
    W = u1 * np.gradient(u2, x) - np.gradient(u1, x) * u2
    r = np.trapezoid((W / u1) ** 2, x) / np.trapezoid(u2 ** 2, x) / lam[0]
    P(f"    bc={bc}: (lam2-lam1)/lam1 = {r:.6f}   3 - this = {3 - r:+.6f}"
      f"   -> rho^2 = {1 + r:.6f}")
    OUT.setdefault("E2", {})[bc] = dict(rho2=float(1 + r), deficit=float(3 - r))
flush()

P("=" * 100)
P("E3  the identity on a two-step string (w=1 on [0,1/2], 10 on [1/2,1]) -- generalised form")
P("=" * 100)
w = lambda t: np.where(t < 0.5, 1.0, 10.0)
x, lam, V = solve(N=1600, w=w)
u1, u2 = V[:, 0], V[:, 1]
W = u1 * np.gradient(u2, x) - np.gradient(u1, x) * u2
lhs = lam[1] - lam[0]
rhs = np.trapezoid((W / u1) ** 2, x) / np.trapezoid(w(x[x >= 0] if False else x) * u2 ** 2, x) if False else None
num = np.trapezoid((W / u1) ** 2, x)
den = np.trapezoid(w(x) * u2 ** 2, x)
P(f"    lam1={lam[0]:.6f} lam2={lam[1]:.6f}  rho^2={lam[1]/lam[0]:.6f}")
P(f"    generalised identity: lam2-lam1 = {lhs:.8f} ;  INT(W/u1)^2/INT(w u2^2) = {num/den:.8f}")
OUT["E3_string"] = dict(rho2=float(lam[1] / lam[0]), lhs=float(lhs), rhs=float(num / den),
                        rel_err=float(abs(num / den - lhs) / lhs))
flush()
P("=" * 100)
P("E4  EXECUTABLE CRITERION (T3) on three models:")
P("    Q := INT (W/u1)^2 / ( 3 * e1 * INT u2^2 )      [rho^2 <= 4  <=>  Q <= 1]")
P("    W = u1 u2' - u1' u2 ;  e1 = lowest tower eigenvalue ;  u1 nodeless is the HYPOTHESIS")
P("=" * 100)


def crit(x, lam, V, ref=0, probe=1):
    """Q for the pair (u_ref, u_probe); e1 = lam[probe]-lam[ref] of the tower is supplied by
    the caller through lam_tower (see below). Here: returns the identity's ratio and the
    two integrals, with u_ref assumed nodeless."""
    u1, u2 = V[:, ref], V[:, probe]
    W = u1 * np.gradient(u2, x) - np.gradient(u1, x) * u2
    A = np.trapezoid((W / u1) ** 2, x)
    B = np.trapezoid(u2 ** 2, x)
    return A, B, A / B


E4 = {}
# ---- (a) DD free string: e1 = pi^2/L^2, e2 = 4 pi^2/L^2  (no vacuum in the tower)
x, lam, V = solve(N=1200, bc="DD")
A, B, r = crit(x, lam, V, 0, 1)
e1 = lam[0]
Q = A / (3 * e1 * B)
P(f"  (a) DD free string    e1={e1:.6f} e2={lam[1]:.6f} rho^2={lam[1]/e1:.8f}")
P(f"      INT(W/u1)^2={A:.8f}  INT u2^2={B:.6f}  ratio={r:.8f} (= e2-e1 = {lam[1]-e1:.8f})")
P(f"      Q = A/(3 e1 B) = {Q:.8f}      (analytic: 1 exactly -> saturated)")
E4["DD_free"] = dict(e1=float(e1), e2=float(lam[1]), A=float(A), B=float(B), Q=float(Q))
# scale invariance check: u2 -> 3 u2
u1c, u2c = V[:, 0], 3.0 * V[:, 1]
Wc = u1c * np.gradient(u2c, x) - np.gradient(u1c, x) * u2c
Qc = np.trapezoid((Wc / u1c) ** 2, x) / (3 * e1 * np.trapezoid(u2c ** 2, x))
P(f"      scale check (u2 -> 3 u2): Q = {Qc:.8f}  (must be unchanged)")
E4["scale_check_Q"] = float(Qc)

# ---- (b) ND half-string
x, lam, V = solve(N=1200, bc="ND")
A, B, r = crit(x, lam, V, 0, 1)
e1 = lam[0]
Q = A / (3 * e1 * B)
P(f"  (b) ND half-string    e1={e1:.6f} e2={lam[1]:.6f} rho^2={lam[1]/e1:.8f}")
P(f"      INT(W/u1)^2={A:.8f}  ratio={r:.8f} (= e2-e1 = {lam[1]-e1:.8f})   Q = {Q:.8f}  (>1 !)")
E4["ND"] = dict(e1=float(e1), e2=float(lam[1]), A=float(A), B=float(B), Q=float(Q))

# ---- (c) double well Lam=20 (vacuum subtracted): u0 nodeless, tower = (E1-E0, E2-E0)
N = 8000
xw = np.linspace(-3, 3, N)
hw = xw[1] - xw[0]
dw = 2.0 / hw ** 2 + 20.0 * (xw ** 2 - 1) ** 2
ew = np.full(N - 1, -1.0 / hw ** 2)
ev, VV = eigh_tridiagonal(dw, ew, select="i", select_range=(0, 3))
A1, B1, r1 = crit(xw, ev, VV, 0, 1)
A2, B2, r2 = crit(xw, ev, VV, 0, 2)
rho2 = r2 / r1
P(f"  (c) double well Lam=20 E0={ev[0]:.6f} E1={ev[1]:.6f} E2={ev[2]:.6f}")
P(f"      two-integral form (nodeless reference u0):")
P(f"        E1-E0={ev[1]-ev[0]:.8e}  vs  INT(W01/u0)^2/INT u1^2 = {r1:.8e}")
P(f"        E2-E0={ev[2]-ev[0]:.8f}  vs  INT(W02/u0)^2/INT u2^2 = {r2:.8f}")
P(f"        rho^2 = r2/r1 = {rho2:.6f}   (direct: {(ev[2]-ev[0])/(ev[1]-ev[0]):.6f})"
  f"   Q = (rho^2-1)/3 = {(rho2-1)/3:.4f}  (>1 !)")
E4["DW20"] = dict(E0=float(ev[0]), E1=float(ev[1]), E2=float(ev[2]),
                  r1=float(r1), r2=float(r2), rho2=float(rho2), Q=float((rho2 - 1) / 3))
# why u1 must be nodeless: single-integral form with the NODAL u1 diverges like 1/eps
u1n, u2n = VV[:, 1], VV[:, 2]
Wn = u1n * np.gradient(u2n, xw) - np.gradient(u1n, xw) * u2n
mask = np.abs(u1n) > 1e-12
P("      single-integral form with the NODAL u1 (hypothesis violated) -> 1/eps divergence:")
E4["node_div"] = {}
for eps in [3e-2, 1e-2, 3e-3, 1e-3, 3e-4]:
    ok = np.abs(u1n) > eps
    Aq = np.trapezoid(np.where(ok, (Wn / np.where(ok, u1n, 1.0)) ** 2, 0.0), xw)
    Qq = Aq / (3 * (ev[1] - ev[0]) * np.trapezoid(u2n ** 2, xw))
    E4["node_div"][f"eps{eps:g}"] = float(Qq)
    P(f"        eps={eps:8.1e}  Q_eps={Qq:12.4e}")
OUT["E4_criterion"] = E4
flush()

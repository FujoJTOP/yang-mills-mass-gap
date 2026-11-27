#!/usr/bin/env python3
"""gm2edge_check.py -- deterministic checks for GM2-EDGE (no Monte Carlo, no sampling).

[1] two-level exact expansion  E(t) = D + t*kappa - t^2*beta^2/(D2-D) + O(t^3)
    => c1 = kappa/D (gap-ratio independent), c2 = -beta^2/(D*(D2-D)) <= 0.
[2] degenerate edge: first-order shift = lambda_min(P1 K2 P1), NOT <psi1|K2|psi1>.
[3] free field: omega_a := (2/a) sinh(a*Delta_a/2) == m exactly (Delta_a from cosh(aD)=1+a^2 m^2/2).
[4] off-diagonal independence of the a^2 coefficient.
[5] concavity in t = a^2 of the exact edge eigenvalue (=> c2 <= 0 at eigenvalue level).
"""
from mpmath import mp, mpf, matrix, eigsy, sinh, cosh, acosh, sqrt, nstr

mp.dps = 40
OUT = []


def p(s=""):
    OUT.append(s)
    print(s)


def evals(H):
    E = eigsy(H, eigvals_only=True)
    return [E[i, 0] for i in range(E.rows)]


def emin(H):
    return min(evals(H))


# ---------------------------------------------------------------- [1]
p("[1] two-level: H(t) = [[D,0],[0,D2]] + t*[[kappa,beta],[beta,kappa2]]")
D, D2, ka, kb, be = mpf("1.7"), mpf("3.1"), mpf("-0.41"), mpf("0.93"), mpf("0.37")


def e1(t):
    return emin(matrix([[D + t * ka, t * be], [t * be, D2 + t * kb]]))


e0 = D  # t=0 exactly (H(0) diagonal)
p("    t          (E(t)-D)/t              kappa - t*beta^2/D         err")
for t in [mpf("1e-2"), mpf("1e-3"), mpf("1e-4"), mpf("1e-6")]:
    est = (e1(t) - e0) / t
    pred = ka - t * be**2 / (D2 - D)
    p("    %-10s %-24s %-24s %s" % (nstr(t, 4), nstr(est, 18), nstr(pred, 18),
                                    nstr(est - pred, 4)))
p("    => c1 = kappa/D = %s   (kappa = %s)" % (nstr(ka / D, 18), nstr(ka, 6)))
p("    => delta1 = c1*D + D^3/24 = %s ; (SM-2) iff delta1 = 0 iff kappa = -D^3/24 = %s"
  % (nstr(ka + D**3 / 24, 18), nstr(-D**3 / 24, 18)))
p("    a^4 coefficient: (E(t)-D-t*kappa)/t^2 at t=1e-6 :")
t = mpf("1e-6")
c2_est = (e1(t) - e0 - t * ka) / t**2
p("      numeric %s   vs  -beta^2/(D2-D) = %s   (c2 = that /D = %s)"
  % (nstr(c2_est, 18), nstr(-be**2 / (D2 - D), 18), nstr(-be**2 / (D2 - D) / D, 18)))
p("    gap-ratio scan (kappa fixed): D2, c1, c2")
for D2v in [mpf("1.9"), mpf("2.5"), mpf("5.0")]:
    e = lambda tt: emin(matrix([[D + tt * ka, tt * be], [tt * be, D2v + tt * kb]]))
    c1v = (e(mpf("1e-8")) - D) / mpf("1e-8")
    c2v = (e(mpf("1e-4")) - D - mpf("1e-4") * ka) / mpf("1e-4") ** 2
    p("      D2=%-5s c1=%s  c2=%s" % (nstr(D2v, 3), nstr(c1v, 12), nstr(c2v, 12)))
p("")

# ---------------------------------------------------------------- [2]
p("[2] degenerate edge (Delta1 = Delta2 = 1.0): first order = lambda_min(P1 K2 P1)")
H0 = matrix([[mpf(1), 0, 0], [0, mpf(1), 0], [0, 0, mpf("2.4")]])
K2 = matrix([[mpf("0.3"), mpf("0.5"), mpf("0.2")],
             [mpf("0.5"), mpf("-0.7"), mpf("0.1")],
             [mpf("0.2"), mpf("0.1"), mpf("0.0")]])
blk = matrix([[mpf("0.3"), mpf("0.5")], [mpf("0.5"), mpf("-0.7")]])
lam_min = min(evals(blk))
lam_max = max(evals(blk))
p("    lambda_min(P1 K2 P1) = %s ; lambda_max = %s ; <e1|K2|e1> = %s"
  % (nstr(lam_min, 18), nstr(lam_max, 18), nstr(mpf("0.3"), 6)))
for t in [mpf("1e-6"), mpf("1e-9")]:
    Ht = H0 + t * K2
    est = (emin(Ht) - 1) / t
    p("      t=%-8s (Emin(t)-1)/t = %s   err vs lambda_min = %s"
      % (nstr(t, 3), nstr(est, 18), nstr(est - lam_min, 4)))
p("    => <psi1|K2|psi1> = +0.3 would give c1 = +0.3/1 (wrong sign side);")
p("       the correct edge-only object is the MIN over the edge eigenspace.")
p("")

# ---------------------------------------------------------------- [3]
p("[3] free field: Delta_a from cosh(a*D)=1+a^2 m^2/2 ; omega_a=(2/a)sinh(aD/2)")
for m in [mpf(1), mpf("0.5"), mpf(3)]:
    for a in [mpf("1e-1"), mpf("1e-2"), mpf("1e-3")]:
        Dv = acosh(1 + a**2 * m**2 / 2) / a
        om = (2 / a) * sinh(a * Dv / 2)
        p("    m=%-4s a=%-7s omega_a - m = %s" % (nstr(m, 3), nstr(a, 3), nstr(om - m, 4)))
p("    => omega_a == m identically (delta1 = 0 by operator identity K2 = -H^3/24).")
p("")

# ---------------------------------------------------------------- [4]
p("[4] off-diagonal independence of the a^2 coefficient (4x4, simple edge)")
import itertools
Dg = [mpf(1), mpf("1.3"), mpf(2.0), mpf("2.7")]
h = matrix(4, 4)
for i in range(4):
    h[i, i] = Dg[i]
k_a = matrix(4, 4)
for i in range(4):
    k_a[i, i] = mpf("0.1") * (i + 1)
for (i, j), val in {}.items():
    pass
k_a[0, 1] = k_a[1, 0] = mpf("0.0")
k_b = matrix(k_a)
k_b[0, 1] = k_b[1, 0] = mpf("0.9")
k_b[1, 2] = k_b[2, 1] = mpf("-0.6")
for tag, K in (("no off-diag", k_a), ("with off-diag", k_b)):
    e = lambda tt: emin(h + tt * K)
    c1v = (e(mpf("1e-8")) - Dg[0]) / mpf("1e-8")
    c2v = (e(mpf("1e-5")) - Dg[0] - mpf("1e-5") * K[0, 0]) / mpf("1e-5") ** 2
    p("    %-14s c1 = %s (K[0,0]=%s)   c2 = %s"
      % (tag, nstr(c1v, 12), nstr(K[0, 0], 4), nstr(c2v, 12)))
p("    => same c1, different c2: c1 sees only P1 K2 P1.")
p("")

# ---------------------------------------------------------------- [5]
p("[5] concavity in t=a^2 of the exact edge eigenvalue (2-level, beta>0)")
t1, t2, t3 = mpf("0.01"), mpf("0.02"), mpf("0.03")
f1, f2, f3 = e1(t1), e1(t2), e1(t3)
p("    second difference (t2 midpoint) = %s  (<=0 expected)"
  % nstr(f3 - 2 * f2 + f1, 6))
p("    => c2 <= 0 at eigenvalue level; c1 unconstrained (sign of kappa free).")

with open("gm2edge_check.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n[log written: gm2edge_check.log]")

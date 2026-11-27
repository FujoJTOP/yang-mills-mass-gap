# -*- coding: utf-8 -*-
"""MMMM1-UBETA — u_P3(beta) = 2*I2(beta)/I1(beta): table, root beta*, verdict.

Deterministic numerics only (NO Monte Carlo). Witness for MMMM1-UBETA.md.
  u_P3(beta) = 2 I2(beta)/I1(beta)        (matches corpus 2I2/I1 = 1.014390200094 at beta=2.5)
  gate u < 1/(e*Delta), Delta = 20        -> 0.0183939720586...
Run: py -3 research/v3/species/MMMM1-UBETA.py
"""
import math

GATE = 1.0 / (math.e * 20.0)                 # 1/(e*Delta)
BETAS = [0.02, 0.036, 0.05, 0.5, 1.0, 2.5]
SZZ = 1.0 / 48.0
TGT = 2.5


def I_series(n, x):
    """I_n(x) = sum_k (x/2)^(2k+n) / (k! (k+n)!), integer n >= 0, x >= 0."""
    h = 0.5 * x
    t = h ** n / math.factorial(n)
    s = 0.0
    k = 0
    while True:
        s += t
        k += 1
        t *= (h * h) / (k * (k + n))
        if t < 1e-19 * s or k > 20000:
            return s


try:
    from scipy.special import iv as _iv

    def u(beta):
        return 0.0 if beta <= 0.0 else 2.0 * _iv(2, beta) / _iv(1, beta)

    SRC = "scipy.special.iv"
except Exception as _e:                       # pragma: no cover
    def u(beta):
        return 0.0 if beta <= 0.0 else 2.0 * I_series(2, beta) / I_series(1, beta)

    SRC = "power series (scipy absent: %r)" % (_e,)


def bisect(f, a, b, rel=1e-16, itmax=400):
    fa = f(a)
    for _ in range(itmax):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0.0:
            b = m
        else:
            a, fa = m, fm
        if b - a <= rel * max(1.0, abs(b)):
            break
    return 0.5 * (a + b)


print("== S0  gate = 1/(e*Delta), Delta=20 ==")
print("GATE = %.16f   (1/GATE = e*Delta = %.10f)" % (GATE, 1.0 / GATE))
print("src  = %s" % SRC)
print("u(2.5) = %.15f  [corpus 1.014390200094]" % u(TGT))
print("series cross-check u(2.5) = %.15f" % (2.0 * I_series(2, TGT) / I_series(1, TGT)))
print("excess  u(2.5)/GATE = %.6f x   [corpus 55.148x]" % (u(TGT) / GATE))

print("\n== S1  table: u_P3(beta) = 2 I2/I1 ==")
print("%8s %18s %10s %10s %10s %8s" % ("beta", "u_P3", "u/GATE", "beta/beta*", "u/(beta/2)", "gate"))
bst = bisect(lambda b: u(b) - GATE, 1e-12, 1.0)
for b in BETAS:
    print("%8.4f %18.15f %10.5f %10.5f %10.9f %8s"
          % (b, u(b), u(b) / GATE, b / bst, u(b) / (0.5 * b), "OK" if u(b) < GATE else "NO"))

print("\n== S2  root beta* of u_P3(beta) = GATE ==")
print("beta* (bisect, dp) = %.16f" % bst)
print("2/(e*Delta)        = %.16f   (leading small-beta law)" % (2.0 * GATE))
print("beta*/(2/eDelta)   = %.16f   (+%.3e rel; predicted +beta^2/24 = +%.3e)"
      % (bst / (2.0 * GATE), bst / (2.0 * GATE) - 1.0, bst * bst / 24.0))
try:
    import mpmath as mp

    mp.mp.dps = 40
    um = lambda b: 2 * mp.besseli(2, b) / mp.besseli(1, b)
    r = mp.findroot(lambda b: um(b) - mp.mpf(1) / (20 * mp.e), mp.mpf(0.0368))
    print("beta* (mpmath dps=40) = %s" % mp.nstr(r, 30))
    print("residual              = %s" % mp.nstr(um(r) - mp.mpf(1) / (20 * mp.e), 5))
except Exception as _e:
    print("mpmath unavailable (%r); dp bisection stands" % (_e,))

print("\n== S3  monotonicity + window structure ==")
# u/2 = I2/I1 = (d/dbeta ln I1) - 1/beta ; log I_n concave => ratio I_{n}/I_{n-1} increasing
prev, mono = None, True
for i in range(1, 4001):
    b = 2.5 * i / 4000.0
    cur = u(b)
    if prev is not None and cur <= prev:
        mono = False
    prev = cur
print("strictly increasing on (0,2.5], 4000-node grid: %s" % mono)
b = 1.0
print("identity check (u/2) vs (I1'/I1 - 1/beta) at beta=1: %.15f vs %.15f"
      % (0.5 * u(b), 0.5 * (_iv(0, b) + _iv(2, b)) / _iv(1, b) - 1.0 / b))

print("\n== S4  ratio verdicts ==")
print("1/48           = %.15f  -> 1/48 < beta* ? %s   (%.6f x of window)"
      % (SZZ, SZZ < bst, SZZ / bst))
print("u(1/48)        = %.15f  (%.6f x of GATE)" % (u(SZZ), u(SZZ) / GATE))
print("beta=2.5 vs beta*: %.6f x" % (TGT / bst))
print("beta=2.5 vs 1/48 : %.6f x" % (TGT / SZZ))
print("window width in beta: (0, %.16f)   [1/48 sits at %.4f%% of it]"
      % (bst, 100.0 * SZZ / bst))
# K-dependence of the endpoint: window = u^{-1}(1/K)
for K, tag in ((math.e * 20.0, "K = e*Delta (task gate)"),
               (8.0 * math.e * 20.0, "K = 8*e*Delta (corpus Cbar route, linear u=beta/8)")):
    K = float(K) if False else K
    print("  %-52s 1/K = %.15f" % (tag, 1.0 / K))
print("ADJ2 beta_1 = 1/(e*Delta) = %.15f ; 1/48 exceeds it by %.6f x"
      % (GATE, SZZ / GATE))

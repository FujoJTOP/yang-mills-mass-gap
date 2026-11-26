"""lc3_ratios.py -- recompute, digit by digit, the two ratio tables in dispute.

  Gamma_n = #{connected n-plaquette sets through p0}      (animal count, P3:89 object)
  m_n     = sum_{Gamma, p0 in Gamma, |Gamma|=n} |phi^U|   (mass/multiplicity weighted)

Both from the SAME enumeration: w7/ur4_enum_n6.json (pool R=10, n<=6).
Exact integer ratios (Fraction), then decimals. Nothing is quoted from TN2/UR4.
"""
import json
import math
from fractions import Fraction as F

d = json.load(open("w7/ur4_enum_n6.json"))
G = [d["res"][str(n)]["sets"] for n in range(1, 7)]
M = [d["res"][str(n)]["m_n"] for n in range(1, 7)]

print("source: w7/ur4_enum_n6.json   pool=%d  R=%d  elapsed=%.1fs"
      % (d["pool"], d["R"], d["elapsed"]))
print("Gamma (=-#connected sets) :", G)
print("m_n   (=sum |phi^U|)      :", M)
print("check sum(dist) == Gamma  :", [sum(d["res"][str(n)]["dist"].values()) == G[n - 1]
                                     for n in range(1, 7)])
print("check sum(value*count)    :", [sum(int(k) * v for k, v in d["res"][str(n)]["dist"].items()) == M[n - 1]
                                     for n in range(1, 7)])


def table(name, xs):
    print("\n== %s : ratios %s_{n}/%s_{n-1} ==" % (name, name, name))
    print(" n   value      EXACT            decimal      ^(1/n)      increment")
    prev = None
    prev_root = None
    for n, x in enumerate(xs, 1):
        r = "" if prev is None else str(F(x, prev))
        dec = float("nan") if prev is None else x / prev
        root = x ** (1.0 / n)
        inc = "" if prev_root is None else "%+.4f" % (root - prev_root)
        print("%2d %9d  %-16s %12.6f %10.6f  %s" % (n, x, r, dec, root, inc))
        prev, prev_root = x, root


table("Gamma", G)
table("m", M)

print("\n== rounding test: which column reproduces the two quoted strings? ==")
q1 = [20, 25.70, 28.69, 30.53, 31.77]
q2 = [20, 22.9, 24.3, 25.17]
rG = [G[i] / G[i - 1] for i in range(1, 6)]
rM = [M[i] / M[i - 1] for i in range(1, 6)]
for lab, q, r in (("20->25.70->28.69->30.53->31.77", q1, rM),
                  ("20->22.9->24.3->25.17", q2, rG)):
    hit = all(abs(a - b) < 0.006 for a, b in zip(q, r))
    print("%-34s vs %s ratios -> %s" % (lab, "m_n" if r is rM else "Gamma", "MATCH" if hit else "MISMATCH"))
    print("     quoted :", q)
    print("     m_n    :", ["%.4f" % v for v in rM])
    print("     Gamma  :", ["%.4f" % v for v in rG])

print("\n== is it the same sequence? ==")
print("TN2 a_n (1,20,458,11132,280175) == Gamma_1..5 :",
      [1, 20, 458, 11132, 280175] == G[:5])
print("first step identical for both (n=2 all |phi^U|=1): m_2=G_2=%d" % G[1])
print("m_n/G_n (mean |phi^U|) =", ["%.6f" % (M[i] / G[i]) for i in range(6)])

print("\n== reciprocal cross-check against 1/48 (the SZZ wall) ==")
for lab, lam in (("Gamma ratio n=5", G[4] / G[3]), ("Gamma ratio n=6", G[5] / G[4]),
                 ("m_n ratio n=6", M[5] / M[4]), ("m_n^(1/n) at n=6", M[5] ** (1 / 6))):
    print("%-18s lam=%9.4f   1/lam=%.6f   vs 1/48=%.6f   ratio=%.3f"
          % (lab, lam, 1 / lam, 1 / 48, (1 / lam) / (1 / 48)))

print("\n== mass-rate extrapolations (for the lambda_Ur^mass slot only) ==")
print("m_n^(1/n)            :", ["%.6f" % (M[n - 1] ** (1 / n)) for n in range(1, 7)])
print("increments           :", ["%+.4f" % (M[n - 1] ** (1 / n) - M[n - 2] ** (1 / (n - 1))) for n in range(2, 7)])
xs = [1 / n for n in (4, 5, 6)]
ys = [M[n - 1] ** (1 / n) for n in (4, 5, 6)]
nn = 3
sx, sy = sum(xs), sum(ys)
sxx = sum(x * x for x in xs)
sxy = sum(x * y for x, y in zip(xs, ys))
b = (nn * sxy - sx * sy) / (nn * sxx - sx * sx)
a = (sy - b * sx) / nn
print("affine fit m_n^(1/n)=lam+b/n on n=4,5,6 (verify UR4's 24.54): lam=%.3f b=%.2f" % (a, b))
print("\nno float ambiguity: all ratios above are exact rationals of integers.")

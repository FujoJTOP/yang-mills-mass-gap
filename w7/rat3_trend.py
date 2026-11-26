# rat3_trend.py -- arithmetic on a_n (n<=7) : ratios, roots, trend fits, 48-test, lambda_* bounds
# exact integers; roots/limits via mpmath (50 dps).  No fitting is called a theorem.
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 50

A = [None, 1, 20, 458, 11132, 280175, 7218768, 189215964]   # a_n, n = 1..7  (rat3_a7.log)
NMAX = 7

P = print
P("=" * 78)
P("a_n (rooted at p0, 4D plaquette graph, deg 20)   n = 1..%d" % NMAX)
P("=" * 78)
for n in range(1, NMAX + 1):
    P("  a_%d = %d" % (n, A[n]))

P("\n-- adjacent ratios (exact fractions) --")
rat = [None, None]
for n in range(2, NMAX + 1):
    r = Fraction(A[n], A[n - 1])
    rat.append(r)
    P("  r_%d = a_%d/a_%d = %d/%d = %.6f" % (n, n, n - 1, r.numerator, r.denominator, float(r)))

P("\n-- increments of the ratio sequence --")
d = []
for n in range(3, NMAX + 1):
    inc = float(rat[n] - rat[n - 1])
    d.append((n, inc))
    P("  r_%d - r_%d = %+.6f" % (n, n - 1, inc))
rr = [(d[i][0], d[i][1] / d[i - 1][1]) for i in range(1, len(d))]
P("  increment ratios: " + ", ".join("%d:%.4f" % (k, v) for k, v in rr))

P("\n-- roots a_n^(1/n) --")
roots = [None]
for n in range(1, NMAX + 1):
    v = mp.mpf(A[n]) ** (mp.mpf(1) / n)
    roots.append(v)
    P("  a_%d^(1/%d) = %s" % (n, n, mp.nstr(v, 15)))

# ---- lambda_* bounds -------------------------------------------------------
P("\n-- lambda_* bounds (lambda_* := lim a_n^(1/n); if the limit exists it equals limsup) --")
P("  [T] lambda_* = limsup_n a_n^(1/n) >= a_7^(1/7) = %s   (limsup >= any term: unconditional)"
  % mp.nstr(roots[7], 15))
P("  [T] was: >= a_6^(1/6) = %s" % mp.nstr(roots[6], 15))
P("  [T] upper (Penrose/tree envelope, e*Delta, Delta=20): lambda_* <= e*20 = %s"
  % mp.nstr(mp.e * 20, 15))
P("  [T] local-ratio form: IF a_{n+1} <= 48 a_n for all n THEN lambda_* <= 48")
P("  observed sup_{n<=7} r_n = %s  -> 48/sup = %.4f (slack)" % (mp.nstr(rat[7], 12), 48 / float(rat[7])))

P("\n-- 48-test: C_48(n) = 48^n / a_n  (margin of the local-ratio envelope) --")
for n in range(2, NMAX + 1):
    c = mp.mpf(48) ** n / mp.mpf(A[n])
    P("  n=%d  48^n/a_n = %.4e" % (n, float(c)))

# ---- extrapolation of the ratio -------------------------------------------------
P("\n-- trend of r_n --")
P("  r_2..r_7 = " + ", ".join("%.4f" % float(rat[n]) for n in range(2, NMAX + 1)))
P("  monotone increasing: %s" % all(rat[n] > rat[n - 1] for n in range(3, NMAX + 1)))
# power-law fit  Delta_n = C n^{-p}  on the increments Delta_n = r_{n+1}-r_n, n = 2..6
import math
xs = [math.log(n) for n, _ in d]
ys = [math.log(v) for _, v in d]
m = len(xs)
sx, sy = sum(xs), sum(ys)
sxx = sum(x * x for x in xs)
sxy = sum(x * y for x, y in zip(xs, ys))
p = -(m * sxy - sx * sy) / (m * sxx - sx * sx)
C = math.exp((sy + p * sx) / m)
P("  fit  Dr_n = C n^(-p):  p = %.4f   C = %.4f  (least squares on log-log, n=2..6)" % (p, C))
tail = sum(C * (NMAX + k) ** (-p) for k in range(1, 4000))
P("  => sum_{n>=7} Dr_n ~ %.4f  =>  r_inf ~ %.3f" % (tail, float(rat[NMAX]) + tail))
# geometric fit on increments (conservative: ratio -> last observed)
g = d[-1][1] / d[-2][1]
tail_g = d[-1][1] * g / (1 - g)
P("  geometric fit (ratio fixed at %.4f):  r_inf ~ %.3f" % (g, float(rat[NMAX]) + tail_g))
P("  naive step-1 extrapolation: a_8 ~ a_7 * r_inf ~ %.3e" % (float(A[NMAX]) * (float(rat[NMAX]) + tail)))

P("\n-- level-6 -> level-7 generation stats (from rat3_a7.log, 4 shard passes) --")
print_raw = 2240297184            # printed as `cand = ...` (sum over the 4 identical passes)
per_pass = print_raw // 4
P("  raw candidate pairs per pass = %d ; avg boundary |dS| at |S|=6 = %.4f"
  % (per_pass, per_pass / A[6]))
P("  distinct 7-sets = %d ; dedup factor = %.4f" % (A[7], per_pass / A[7]))
P("  => a_8 scale: raw pairs ~ a_7 * |dS|(7) ~ %.2e  (|dS| slope ~11.1/step => |dS|~89)"
  % (A[7] * 89.0))
P("  => a_8 needs >= %.0f shard-passes at 16 B/key and 3 GB/pass" % (A[7] * 89.0 * 16 / 3e9))

P("\n-- verdict slots --")
P("  (a) exact: a_7 = %d   a_7/a_6 = %.6f   a_7^(1/7) = %s" % (A[7], float(rat[7]), mp.nstr(roots[7], 12)))
P("  ratio <= 48 at n=7 ?  %s   (48 - r_7 = %.4f)" % (float(rat[7]) <= 48, 48 - float(rat[7])))
P("  ratio still increasing ? %s   (r_7 - r_6 = %+.4f)" % (rat[7] > rat[6], float(rat[7] - rat[6])))
P("  any sign of approaching 48 ? r_7/48 = %.4f ; extrapolated r_inf/48 ~ %.4f"
  % (float(rat[7]) / 48, (float(rat[NMAX]) + tail) / 48))

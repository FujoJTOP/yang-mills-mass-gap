import numpy as np
from scipy.special import iv

# FIX: c_d must be evaluated at real d (j-units), which can be half-integer.
# Previous script indexed a precomputed integer array -> wrong whenever d is half-integer.
#   c(d) = exp(-b) * ( I_{2d}(b) - I_{2d+2}(b) )
#   M[x,y] = sum_{t=|x-y|}^{x+y} c(t),  t stepping by 1 in j-units (half-integers allowed)
#   minor  = M[j1,k1] M[j2,k2] - M[j1,k2] M[j2,k1]

def cfun(d, b):
    return np.exp(-b) * (iv(2 * d, b) - iv(2 * d + 2, b))

def Mxy(x, y, b):
    lo, hi = abs(x - y), x + y
    n = int(round(hi - lo))
    return sum(cfun(lo + t, b) for t in range(n + 1))

def minor(j1, j2, k1, k2, b):
    return Mxy(j1, k1, b) * Mxy(j2, k2, b) - Mxy(j1, k2, b) * Mxy(j2, k1, b)

print("=== sanity: Mxy vs Haar quadrature (must match) ===")
om = np.linspace(0, np.pi, 600001); meas = (2 / np.pi) * np.sin(om) ** 2
b = 1.5
f = np.exp(-b * (1 - np.cos(om)))
def ch(j, o):
    s = np.sin(o); num = np.sin((2 * j + 1) * o)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))
for (jp, j) in [(0.0, 0.0), (0.5, 0.0), (1.5, 0.0), (1.5, 0.5), (2.0, 1.5)]:
    q = np.trapezoid(ch(jp, om) * f * ch(j, om) * meas, om)
    print("  (j'=%3.1f,j=%3.1f): quad=%.10f  Mxy=%.10f  diff=%.2e" % (jp, j, q, Mxy(jp, j, b), abs(q - Mxy(jp, j, b))))

print()
print("=== all 2x2 minors, j in {0,1/2,...,5} ===")
js = [k / 2.0 for k in range(0, 11)]
mn = np.inf; arg = None; neg = 0; tot = 0
for i1 in range(len(js)):
    for i2 in range(i1 + 1, len(js)):
        for k1i in range(len(js)):
            for k2 in range(k1i + 1, len(js)):
                v = minor(js[i1], js[i2], js[k1i], js[k2], b)
                tot += 1
                if v < -1e-12:
                    neg += 1
                if v < mn:
                    mn = v; arg = (js[i1], js[i2], js[k1i], js[k2])
print("  total = %d ; negative = %d ; min = %+.6e at %s" % (tot, neg, mn, arg))
print()
print("=== same for several b ===")
for bb in (0.5, 1.0, 1.5, 3.0, 6.0):
    mn = np.inf; neg = 0
    for i1 in range(len(js)):
        for i2 in range(i1 + 1, len(js)):
            for k1i in range(len(js)):
                for k2 in range(k1i + 1, len(js)):
                    v = minor(js[i1], js[i2], js[k1i], js[k2], bb)
                    if v < -1e-12: neg += 1
                    mn = min(mn, v)
    print("  b=%.1f : negative minors = %d ; min = %+.4e" % (bb, neg, mn))
